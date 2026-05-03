#!/usr/bin/env nextflow

/*
 * =============================================================================
 * test_rfdiffusion.nf — Isolated test for the RFDiffusion module
 * =============================================================================
 * Runs RFDIFFUSION → RFDIFFUSION_FILTER → RFDIFFUSION_PLOTS in isolation,
 * using a pre-docked complex PDB (best_model.pdb from a HADDOCK3 run) as
 * input.
 *
 * Contigs and hotspot are provided directly — no upstream HADDOCK or
 * BUILD_CONTIGS steps are needed.
 *
 * Usage:
 *   sbatch tests/rfdiffusion/run_test_rfdiffusion.sh
 *
 * =============================================================================
 */

nextflow.enable.dsl = 2

// ---------------------------------------------------------------------------
// Parameter defaults — override via params.yml or --param on the command line
// ---------------------------------------------------------------------------

params.pdb_file          = "${projectDir}/data/af3_pikp1_native_avrpikf_complex.pdb"
params.receptor_chain    = "A"
params.effector_chain    = "C"
params.contigs           = "A1-32/10-20/A46-72/6-6 C"
params.hotspot           = ""              // leave blank to test without hotspots
params.num_designs       = 8               // small number for fast test
params.rfdiff_iterations = 50
params.rfdiff_contact_cutoff = 8.0
params.min_hotspot_frac  = 0.0             // 0.0 = no filtering (test all designs pass)
params.project_name      = "test_rfdiffusion"
params.outdir            = "${projectDir}/results"

// Infrastructure — params.rfdiff_container is inherited from nextflow.config
// (single source of truth across the main pipeline and per-module tests).

// ---------------------------------------------------------------------------
// Includes
// ---------------------------------------------------------------------------

include { RFDIFFUSION        } from '../../modules/rfdiffusion'
include { RFDIFFUSION_FILTER } from '../../modules/rfdiffusion'
include { RFDIFFUSION_PLOTS  } from '../../modules/rfdiffusion'

// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------

workflow {

    pdb_ch     = Channel.fromPath(params.pdb_file, checkIfExists: true)
    contigs_ch = Channel.value(params.contigs)
    hotspot_ch = Channel.value(params.hotspot)

    RFDIFFUSION(
        pdb_ch,
        contigs_ch,
        hotspot_ch,
        params.num_designs,
        params.rfdiff_iterations,
        Channel.value(file("${projectDir}/bin/rfdiffusion_contigs.py"))
    )

    RFDIFFUSION_FILTER(
        RFDIFFUSION.out.design_pdbs,
        pdb_ch,
        contigs_ch,
        hotspot_ch,
        params.receptor_chain,
        params.effector_chain,
        params.rfdiff_contact_cutoff,
        params.min_hotspot_frac,
        Channel.value(file("${projectDir}/bin/rfdiffusion_filter.py"))
    )

    RFDIFFUSION_PLOTS(
        RFDIFFUSION_FILTER.out.metrics,
        Channel.value(file("${projectDir}/bin/rfdiffusion_plots.py"))
    )
}

// ---------------------------------------------------------------------------
// On completion
// ---------------------------------------------------------------------------

workflow.onComplete {
    log.info """
    =============================================================
    RFDiffusion Test Complete
    =============================================================
    Output:   ${params.outdir}
    Duration: ${workflow.duration}
    Success:  ${workflow.success}
    =============================================================
    """.stripIndent()
}

workflow.onError {
    log.error "RFDiffusion test failed: ${workflow.errorMessage}"
}
