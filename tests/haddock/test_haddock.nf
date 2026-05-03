#!/usr/bin/env nextflow

/*
 * =============================================================================
 * test_haddock.nf — Isolated test for the HADDOCK3 module
 * =============================================================================
 * Runs HADDOCK3_PREPARE → HADDOCK3_DOCK → HADDOCK3_PLOTS → EXTRACT_HOTSPOTS
 * → BUILD_CONTIGS in isolation from the rest of the pipeline.
 *
 * Usage:
 *   sbatch tests/haddock/run_test_haddock.sh
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
params.haddock_sampling  = 100    // reduced from 10000 for faster test runs
params.haddock_seletop   = 20
params.rfdiff_contact_cutoff = 8.0
params.effector_active_residues = ""   // Comma-separated effector residues for HADDOCK AIRs
params.receptor_seq      = null   // Optional: reference sequence for chain disambiguation
params.effector_seq      = null   // Optional: reference sequence for chain disambiguation
params.project_name      = "test_haddock"
params.outdir            = "${projectDir}/results"

// Infrastructure — params.rfdiff_container is inherited from nextflow.config
// (single source of truth across the main pipeline and per-module tests).

// ---------------------------------------------------------------------------
// Includes
// ---------------------------------------------------------------------------

include { HADDOCK3_PREPARE     } from '../../modules/haddock'
include { HADDOCK3_DOCK        } from '../../modules/haddock'
include { HADDOCK3_PLOTS       } from '../../modules/haddock'
include { EXTRACT_HOTSPOTS     } from '../../modules/haddock'
include { BUILD_CONTIGS        } from '../../modules/haddock'
include { WRITE_DUMMY_MAPPING as WRITE_DUMMY_MAPPING_REC } from '../../modules/preprocessing'
include { WRITE_DUMMY_MAPPING as WRITE_DUMMY_MAPPING_EFF } from '../../modules/preprocessing'

// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------

workflow {

    receptor_ch    = Channel.fromPath(params.receptor_input, checkIfExists: true)
    effector_ch    = Channel.fromPath(params.effector_input, checkIfExists: true)

    // Generate empty trim mappings (placeholder for the BUILD_CONTIGS API).
    WRITE_DUMMY_MAPPING_REC(Channel.value("receptor"))
    WRITE_DUMMY_MAPPING_EFF(Channel.value("effector"))

    HADDOCK3_PREPARE(
        receptor_ch,
        effector_ch,
        params.receptor_chain,
        params.effector_chain,
        params.contigs,
        params.effector_active_residues,
        Channel.value(file("${projectDir}/bin/haddock3_prepare.py"))
    )

    HADDOCK3_DOCK(
        HADDOCK3_PREPARE.out.receptor_pdb_out,
        HADDOCK3_PREPARE.out.effector_pdb_out,
        HADDOCK3_PREPARE.out.restraints,
        params.haddock_sampling,
        params.haddock_seletop,
        Channel.value(file("${projectDir}/bin/collect_haddock3_dock.py"))
    )

    HADDOCK3_PLOTS(
        HADDOCK3_DOCK.out.capri_scores,
        HADDOCK3_DOCK.out.cluster_summary,
        HADDOCK3_DOCK.out.run_dir,
        params.contigs,
        params.receptor_chain,
        params.effector_active_residues,
        Channel.value(file("${projectDir}/bin/haddock3_plots.py"))
    )

    EXTRACT_HOTSPOTS(
        HADDOCK3_DOCK.out.best_model,
        params.receptor_chain,
        params.effector_chain,
        params.rfdiff_contact_cutoff,
        params.receptor_seq ?: "",
        params.effector_seq ?: "",
        Channel.value(file("${projectDir}/bin/extract_hotspots.py"))
    )

    BUILD_CONTIGS(
        HADDOCK3_DOCK.out.best_model,
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
