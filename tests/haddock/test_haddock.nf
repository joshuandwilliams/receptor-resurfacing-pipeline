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

// Real biological inputs: Pikp-1_HMA (78-residue HMA domain, chain A)
// and avr-Pia (68-residue MAX-fold effector, chain B) as monomer PDBs.
// Both renumbered to start at residue 1 in ChimeraX before saving.
params.receptor_input    = "${projectDir}/data/Pikp-1_HMA.pdb"
params.effector_input    = "${projectDir}/data/avr-pia.pdb"
params.receptor_chain    = "A"
params.effector_chain    = "B"
// Contig: A1-32 + A50-68 are fixed (51 native residues kept), with two
// de novo regions (33-49 → length 10-30, 69-78 → length 10).  The contig's
// fixed segments are what HADDOCK_PREPARE uses to derive the receptor
// design region (residues 33-49 + 69-78, 27 total) for clash bookkeeping.
params.contigs           = "A1-32/10-30/A50-68/10-10 B"

// Production-length sampling — the short (sampling=100) test only ever
// produced one qualifying cluster, hiding multi-cluster behaviour.
// At 10000 / 400 / min 4 we typically see 3-6 qualifying clusters
// with distinct poses, which exercises the auto-pick logic and the
// per-cluster plots properly.  Estimated runtime: ~30–60 min on
// 64 CPUs at the jic-medium queue's typical wait time.
params.haddock_sampling  = 10000
params.haddock_seletop   = 400
params.haddock_min_cluster_size = 4

// Session 7 restraint params (post-commit-3 amendment: effector-only AIRs
// are now valid, and the receptor design region comes from the contig).
//
// Contact-pair mode: 3 hard CA-CA pins on antiparallel beta strand
// contacts the user wants locked in place (core of the intended
// receptor:effector interface; the receptor halves fall in the second
// de novo region, which is the expected "HADDOCK pins to native coords;
// RFDiffusion redesigns the residue identities" workflow).
params.haddock_contact_pairs            = "A73-B31 A72-B32 A71-B33"
// Receptor design region 1 (33-49) on the active list — combined with
// the effector hotspot list below this becomes a two-sided AIR that
// pulls DR1 residues toward B20-26.  DR2 (residues 69-78) is intentionally
// EXCLUDED here because the contact pairs already pin DR2 to the existing
// AvrPikF beta-strand interface; adding DR2 to the active list would
// double-restrain.  See notes/design_audit.md "DR1 rotation experiment".
params.haddock_receptor_active_residues = "33-49"
params.haddock_effector_active_residues = "20-26"
params.haddock_pair_distance            = "2,2,4"
params.haddock_chosen_cluster           = null   // auto-pick
// Strip the receptor's design-region sidechains (33-49 + 69-78) before
// HADDOCK runs — leaves backbone-only GLY in those positions so the
// effector can pack tight against the receptor backbone without
// being deflected by sidechains that RFDiffusion will redesign anyway.
// Per A147 future-consideration in notes/design_audit.md.
params.haddock_strip_design_sidechains  = true

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
include { HADDOCK_CLUSTER_SC       } from '../../modules/haddock'
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
        params.contigs,
        params.haddock_strip_design_sidechains,
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

    HADDOCK_CLUSTER_METRICS(
        HADDOCK3_DOCK.out.haddock_report,
        HADDOCK3_DOCK.out.cluster_models,
        HADDOCK3_PREPARE.out.restraints_summary,
        HADDOCK3_PREPARE.out.ambig_restraints,
        HADDOCK3_PREPARE.out.unambig_restraints,
        Channel.value(file("${projectDir}/bin/haddock_cluster_metrics.py")),
        Channel.value(file("${projectDir}/bin"))
    )

    HADDOCK_CLUSTER_SC(
        HADDOCK3_DOCK.out.haddock_report,
        HADDOCK3_DOCK.out.cluster_models,
        Channel.value(file("${projectDir}/bin/haddock_cluster_sc.py")),
        Channel.value(file("${projectDir}/bin"))
    )

    // PLOTS depends on cluster_metrics + cluster_sc, so it runs after them.
    HADDOCK3_PLOTS(
        HADDOCK3_DOCK.out.capri_scores,
        HADDOCK3_DOCK.out.cluster_summary,
        HADDOCK3_DOCK.out.run_dir,
        HADDOCK3_PREPARE.out.restraints_summary,
        HADDOCK_CLUSTER_METRICS.out.cluster_metrics,
        HADDOCK_CLUSTER_SC.out.cluster_sc,
        params.haddock_contact_pairs,
        params.haddock_effector_active_residues,
        Channel.value(file("${projectDir}/bin/haddock3_plots.py"))
    )

    SELECT_HADDOCK_CLUSTER(
        HADDOCK3_DOCK.out.haddock_report,
        HADDOCK_CLUSTER_METRICS.out.cluster_metrics,
        HADDOCK_CLUSTER_SC.out.cluster_sc,
        HADDOCK3_PREPARE.out.restraints_summary,
        HADDOCK3_DOCK.out.cluster_models,
        receptor_ch,
        effector_ch,
        params.haddock_contact_pairs,
        params.haddock_receptor_active_residues,
        params.haddock_effector_active_residues,
        params.haddock_pair_distance,
        // Nextflow refuses to bind null to a `val` input — string-sentinel.
        (params.haddock_chosen_cluster != null ? "${params.haddock_chosen_cluster}" : "null"),
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
