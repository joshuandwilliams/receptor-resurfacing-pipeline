#!/usr/bin/env nextflow

/*
 * =============================================================================
 * test_rosetta_filtering.nf — Isolated test for the Rosetta filtering module
 * =============================================================================
 * Runs ROSETTA_SC → ROSETTA_FILTER → ROSETTA_FILTER_PLOTS in isolation,
 * using split two-chain PDBs from a previous RFDIFFUSION_FILTER run as
 * input.
 *
 * The split PDBs must have chain A (receptor) and chain B (effector),
 * renumbered from 1 — these are the polyvaline backbone structures
 * produced by RFDiffusion and split by rfdiffusion_filter.py.
 *
 * Test chaining
 * -------------
 * By default this consumes the output of ``test_rfdiffusion`` so the
 * tests form a cascade when run in order:
 *
 *   test_rfdiffusion/results/rfdiffusion/split/design_*.pdb   (produced)
 *                            │
 *                            ▼
 *   test_rosetta_filtering    (consumed here)
 *
 * Point ``--upstream_rfdiff_outdir`` elsewhere to consume a different
 * upstream run's output, or override ``--design_pdbs`` directly with a
 * glob to the ``tests/rosetta_filtering/data/design_*.pdb`` test data
 * for standalone runs.
 *
 * Usage:
 *   sbatch tests/rosetta_filtering/run_test_rosetta_filtering_slurm.sh
 * =============================================================================
 */

nextflow.enable.dsl = 2

// ---------------------------------------------------------------------------
// Parameter defaults — override via params.yml or --param on the command line
// ---------------------------------------------------------------------------

// ── Test-chaining: prefer upstream test output, fall back to cached data/ ──
// If test_rfdiffusion has been run, rfdiffusion_filter.py splits every
// passing design into separate receptor+effector chains under
// <outdir>/rfdiffusion/split/.  We look there first; if the glob is
// empty (test_rfdiffusion hasn't run, or published elsewhere), we fall
// back to the stand-alone cached test data in this test's data/ dir.
//
// Override ``--upstream_rfdiff_outdir`` to point at a different upstream
// run, or ``--design_pdbs`` directly to skip the chaining logic entirely.
params.upstream_rfdiff_outdir = "${projectDir}/../rfdiffusion/receptor_resurfacing_results"

def _upstream_design_glob = "${params.upstream_rfdiff_outdir}/rfdiffusion/split/design_*.pdb"
def _fallback_design_glob = "${projectDir}/data/design_*.pdb"

// files(glob) returns a (possibly empty) list of matching paths.  Use
// cached data/ only if the upstream test hasn't produced output yet.
params.design_pdbs = files(_upstream_design_glob).size() > 0 \
    ? _upstream_design_glob                                  \
    : _fallback_design_glob

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

    // Log which input source was chosen — upstream chain or cached data/.
    def _selected = params.design_pdbs.startsWith(params.upstream_rfdiff_outdir) \
        ? "upstream (${params.upstream_rfdiff_outdir})"                          \
        : "cached test data (${projectDir}/data/)"
    log.info "Rosetta filtering test — input source: ${_selected}"
    log.info "design_pdbs glob: ${params.design_pdbs}"

    design_pdbs_ch = Channel.fromPath(params.design_pdbs, checkIfExists: true)

    ROSETTA_SC(
        design_pdbs_ch
    )

    ROSETTA_FILTER(
        ROSETTA_SC.out.pdb_and_scores.collect(),
        params.sc_threshold
    )

    ROSETTA_FILTER_PLOTS(
        ROSETTA_FILTER.out.metrics
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
