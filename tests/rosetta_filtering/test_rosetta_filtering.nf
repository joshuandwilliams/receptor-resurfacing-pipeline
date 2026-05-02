#!/usr/bin/env nextflow

/*
 * =============================================================================
 * test_rosetta_filtering.nf — Isolated test for the Rosetta filtering module
 * =============================================================================
 * Runs ROSETTA_SC → ROSETTA_FILTER → ROSETTA_FILTER_PLOTS in isolation,
 * using split two-chain PDBs from ``tests/rosetta_filtering/data/`` as
 * input.
 *
 * The split PDBs must have chain A (receptor) and chain B (effector),
 * renumbered from 1 — these are the polyvaline backbone structures
 * produced by RFDiffusion and split by rfdiffusion_filter.py.
 *
 * Inputs (canonical, no upstream chaining)
 * ----------------------------------------
 *   data/rfdiffusion_split/design_*.pdb
 *       Curated split-PDB fixtures captured from a discovery run's
 *       RFDIFFUSION_FILTER ``split/`` output.  See data/README.md.
 *
 * Usage:
 *   sbatch tests/rosetta_filtering/run_test_rosetta_filtering_slurm.sh
 * =============================================================================
 */

nextflow.enable.dsl = 2

// ---------------------------------------------------------------------------
// Parameter defaults — override via params.yml or --param on the command line
// ---------------------------------------------------------------------------

// Canonical input path: this test's own data/ directory.  No fallback to
// upstream test outputs — per-module tests run from committed fixtures.
params.design_pdbs = "${projectDir}/data/rfdiffusion_split/design_*.pdb"

params.sc_threshold      = 0.5
params.project_name      = "test_rosetta_filtering"
params.outdir            = "${projectDir}/results"

// Infrastructure — params.rosetta_container is inherited from nextflow.config
// (single source of truth across the main pipeline and per-module tests).

// ---------------------------------------------------------------------------
// Includes
// ---------------------------------------------------------------------------

include { ROSETTA_SC           } from '../../modules/rosetta_filtering'
include { ROSETTA_FILTER       } from '../../modules/rosetta_filtering'
include { ROSETTA_FILTER_PLOTS } from '../../modules/rosetta_filtering'

// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------

workflow {

    log.info "Rosetta filtering test — design_pdbs glob: ${params.design_pdbs}"

    design_pdbs_ch = Channel.fromPath(params.design_pdbs, checkIfExists: true)

    ROSETTA_SC(
        design_pdbs_ch
    )

    ROSETTA_FILTER(
        ROSETTA_SC.out.pdb_and_scores.collect(),
        params.sc_threshold,
        Channel.value(file("${projectDir}/bin/rosetta_filter_collect.py"))
    )

    ROSETTA_FILTER_PLOTS(
        ROSETTA_FILTER.out.metrics,
        Channel.value(file("${projectDir}/bin/rosetta_filter_plots.py"))
    )
}

// ---------------------------------------------------------------------------
// On completion
// ---------------------------------------------------------------------------

workflow.onComplete {
    log.info """
    =============================================================
    Rosetta Filtering Test Complete
    =============================================================
    Output:   ${params.outdir}
    Duration: ${workflow.duration}
    Success:  ${workflow.success}
    =============================================================
    """.stripIndent()
}

workflow.onError {
    log.error "Rosetta filtering test failed: ${workflow.errorMessage}"
}
