#!/usr/bin/env nextflow

/*
 * =============================================================================
 * test_haddock.nf — Isolated test for the HADDOCK3 module
 * =============================================================================
 * Per Session 7 restructure (notes/design_audit.md A149), one combined
 * test on real biological inputs exercises both restraint modes
 * simultaneously:
 *   - contact-pair mode (haddock_contact_pairs)
 *   - active-residues mode (haddock_receptor_active_residues
 *     + haddock_effector_active_residues)
 *
 * Pipeline exercised:
 *   HADDOCK3_PREPARE → HADDOCK3_DOCK → HADDOCK3_PLOTS →
 *   HADDOCK_CLUSTER_METRICS → SELECT_HADDOCK_CLUSTER → BUILD_CONTIGS
 *
 * (EXTRACT_HOTSPOTS deleted in Session 7 / commit 3 of the restructure.)
 *
 * Usage:
 *   sbatch tests/haddock/run_test_haddock.slurm.sh
 *
 * =============================================================================
 */

nextflow.enable.dsl = 2

// ---------------------------------------------------------------------------
// Parameter defaults — override via params.yml or --param on the command line
// ---------------------------------------------------------------------------

params.receptor_input    = "${projectDir}/data/sr50_3bi_lrr.pdb"
params.effector_input    = "${projectDir}/data/pwl2.pdb"
params.receptor_chain    = "A"
params.effector_chain    = "B"
params.contigs           = "B A1-390/20-40/A421-438"

// HADDOCK sampling reduced for test runs; production default 10000.
params.haddock_sampling  = 100
params.haddock_seletop   = 20
// Reduced from production default 4: with sampling=100 the clusters
// produced by clustfcc are naturally smaller (1 cluster of ~3 models is
// typical at this scale).  The test exists to exercise the pipeline
// plumbing end-to-end; 2 is the floor that still demands real clustering
// signal.  Production keeps haddock_min_cluster_size = 4 via nextflow.config.
params.haddock_min_cluster_size = 2

// Session 7 restraint params — exercise BOTH modes in the same test.
// Contact-pair mode: 2 hard CA-CA pins on plausible interface residues
// (chosen to exercise the plumbing, not to recover a known interface —
// the sr50/pwl2 interface biology isn't required for this test).
params.haddock_contact_pairs            = "A395-B45 A415-B70"
// Active-residues mode: receptor design region + effector face.  The
// receptor active list covers residues 391-420 (the de novo gap region
// from the contig string) plus a few flanking anchors.  The effector
// list covers a contiguous run on chain B.
params.haddock_receptor_active_residues = "391-420"
params.haddock_effector_active_residues = "40-80"
params.haddock_pair_distance            = "2,2,4"
params.haddock_chosen_cluster           = null   // auto-pick

params.rfdiff_contact_cutoff = 8.0
params.project_name          = "test_haddock"
params.outdir                = "${projectDir}/results"

// Infrastructure — params.rfdiff_container is inherited from nextflow.config.

// ---------------------------------------------------------------------------
// Includes
// ---------------------------------------------------------------------------

include { HADDOCK3_PREPARE         } from '../../modules/haddock'
include { HADDOCK3_DOCK            } from '../../modules/haddock'
include { HADDOCK3_PLOTS           } from '../../modules/haddock'
include { HADDOCK_CLUSTER_METRICS  } from '../../modules/haddock'
include { SELECT_HADDOCK_CLUSTER   } from '../../modules/haddock'
include { BUILD_CONTIGS            } from '../../modules/haddock'
include { WRITE_DUMMY_MAPPING as WRITE_DUMMY_MAPPING_REC } from '../../modules/preprocessing'
include { WRITE_DUMMY_MAPPING as WRITE_DUMMY_MAPPING_EFF } from '../../modules/preprocessing'

// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------

workflow {

    receptor_ch = Channel.fromPath(params.receptor_input, checkIfExists: true)
    effector_ch = Channel.fromPath(params.effector_input, checkIfExists: true)

    // Generate empty trim mappings (placeholder for the BUILD_CONTIGS API).
    WRITE_DUMMY_MAPPING_REC(Channel.value("receptor"))
    WRITE_DUMMY_MAPPING_EFF(Channel.value("effector"))

    HADDOCK3_PREPARE(
        receptor_ch,
        effector_ch,
        params.receptor_chain,
        params.effector_chain,
        params.haddock_contact_pairs,
        params.haddock_receptor_active_residues,
        params.haddock_effector_active_residues,
        params.haddock_pair_distance,
        Channel.value(file("${projectDir}/bin/haddock3_prepare.py"))
    )

    HADDOCK3_DOCK(
        HADDOCK3_PREPARE.out.receptor_pdb_out,
        HADDOCK3_PREPARE.out.effector_pdb_out,
        HADDOCK3_PREPARE.out.ambig_restraints,
        HADDOCK3_PREPARE.out.unambig_restraints,
        params.haddock_sampling,
        params.haddock_seletop,
        Channel.value(file("${projectDir}/bin/collect_haddock3_dock.py"))
    )

    HADDOCK3_PLOTS(
        HADDOCK3_DOCK.out.capri_scores,
        HADDOCK3_DOCK.out.cluster_summary,
        HADDOCK3_DOCK.out.run_dir,
        params.haddock_receptor_active_residues,
        params.haddock_effector_active_residues,
        Channel.value(file("${projectDir}/bin/haddock3_plots.py"))
    )

    HADDOCK_CLUSTER_METRICS(
        HADDOCK3_DOCK.out.haddock_report,
        HADDOCK3_DOCK.out.cluster_models,
        HADDOCK3_PREPARE.out.restraints_summary,
        HADDOCK3_PREPARE.out.ambig_restraints,
        HADDOCK3_PREPARE.out.unambig_restraints,
        Channel.value(file("${projectDir}/bin/haddock_cluster_metrics.py")),
        Channel.value(file("${projectDir}/bin"))
    )

    SELECT_HADDOCK_CLUSTER(
        HADDOCK3_DOCK.out.haddock_report,
        HADDOCK_CLUSTER_METRICS.out.cluster_metrics,
        HADDOCK3_PREPARE.out.restraints_summary,
        HADDOCK3_DOCK.out.cluster_models,
        receptor_ch,
        effector_ch,
        params.haddock_contact_pairs,
        params.haddock_receptor_active_residues,
        params.haddock_effector_active_residues,
        params.haddock_pair_distance,
        params.haddock_chosen_cluster,
        Channel.value(file("${projectDir}/bin/select_haddock_cluster.py"))
    )

    BUILD_CONTIGS(
        SELECT_HADDOCK_CLUSTER.out.selected_pdb,
        params.receptor_chain,
        params.effector_chain,
        params.contigs,
        WRITE_DUMMY_MAPPING_REC.out.mapping,
        WRITE_DUMMY_MAPPING_EFF.out.mapping,
        Channel.value(file("${projectDir}/bin/build_contigs.py"))
    )
}

// ---------------------------------------------------------------------------
// On completion
// ---------------------------------------------------------------------------

workflow.onComplete {
    log.info """
    =============================================================
    HADDOCK3 Test Complete
    =============================================================
    Output:   ${params.outdir}
    Duration: ${workflow.duration}
    Success:  ${workflow.success}
    =============================================================
    """.stripIndent()
}

workflow.onError {
    log.error "HADDOCK3 test failed: ${workflow.errorMessage}"
}
