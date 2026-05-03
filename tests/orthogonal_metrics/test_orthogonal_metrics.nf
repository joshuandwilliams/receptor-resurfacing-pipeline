#!/usr/bin/env nextflow

/*
 * =============================================================================
 * test_orthogonal_metrics.nf
 * =============================================================================
 * P0-31 test harness.  Runs the full post-negative-steering orthogonal-
 * metrics cascade on a completed test_negative_steering run:
 *
 *   NEGSTEER_CROSS_SEQUENCE output (cached)
 *      │
 *      ▼
 *   NEGSTEER_INTERFACE_METRICS               (P0-29 — iRMSD/fnat/DockQ/iPSAE_15/intact_core)
 *      │
 *      ▼
 *   EXTRACT_SURVIVOR_MANIFEST                (manifest of per-survivor tuples)
 *      │
 *      ├──→  AF3_SETUP_DB  ──→  AF3_NOMSA_ON_SURVIVORS  ──→  AF3_PARSE_OUTPUT
 *      ├──→  NEGSTEER_BIOPHYSICAL_METRICS                                       (fan-out)
 *      └──→  NEGSTEER_ROSETTA_METRICS                                           (fan-out)
 *      │
 *      ▼
 *   NEGSTEER_ORTHOGONAL_METRICS              (merge + filter cascade)
 *
 * Inputs (canonical, no upstream chaining)
 * ----------------------------------------
 *   data/negsteer_run/cross_sequence_summary.csv
 *       Cohort table from a curated negative-steering fixture run.
 *   data/negsteer_run/runs/<seq_name>/
 *       Per-sequence workdirs containing plan.json,
 *       effector_template.cif, cycle_0/..., etc.  See data/README.md
 *       for the expected per-workdir layout.
 *
 * Usage:
 *   sbatch tests/orthogonal_metrics/run_test_orthogonal_metrics_slurm.sh
 * =============================================================================
 */

nextflow.enable.dsl = 2

// ---------------------------------------------------------------------------
// Parameter defaults
// ---------------------------------------------------------------------------

params.project_name   = "test_orthogonal_metrics"
params.outdir         = "${projectDir}/results"

// Canonical input path: this test's own data/ directory.  No fallback to
// upstream test outputs — per-module tests run from committed fixtures.
params.input_dir = "${projectDir}/data/negsteer_run"

// ── Chain conventions ────────────────────────────────────────────────
// Per pipeline_notes10 §"chain-param wiring bug":
//   receptor_chain / effector_chain      — INPUT PDB convention.
//   rfdiff_output_*_chain                — prediction-PDB convention
//                                          (always A/B; hardcoded by
//                                          rfdiffusion_filter.py).
// Every orthogonal-metrics module operates on prediction PDBs (Boltz
// canonical_pdb + RFDiffusion split ground_truth — both A/B), so all
// five — NEGSTEER_INTERFACE_METRICS, EXTRACT_SURVIVOR_MANIFEST,
// AF3_NOMSA_ON_SURVIVORS, NEGSTEER_BIOPHYSICAL_METRICS,
// NEGSTEER_ROSETTA_METRICS — read params.rfdiff_output_*_chain.
// receptor_chain / effector_chain are also set to A/B because this
// test consumes pre-baked prediction-PDB workdirs only — there is no
// input PDB present in the fixture.
params.receptor_chain               = "A"
params.effector_chain               = "B"
params.rfdiff_output_receptor_chain = "A"
params.rfdiff_output_effector_chain = "B"

// ── P0-29 defaults (forwarded to NEGSTEER_INTERFACE_METRICS) ──────────
params.interface_plddt_trim_threshold = 50.0
params.interface_intact_threshold     = 5.0

// ── P0-38 default (also forwarded to NEGSTEER_INTERFACE_METRICS) ──────
// See main.nf for the rationale on hard-coded μ/σ vs CLI cutoff.
params.weighted_jaccard_pair_cutoff   = 8.0

// ── P0-31 defaults ────────────────────────────────────────────────────
params.af3_nomsa_seeds              = [42, 123, 456]
params.orthogonal_filter_sc_min     = 0.55
params.orthogonal_filter_bsa_min    = 600
params.orthogonal_filter_plddt_min  = 0.75
params.orthogonal_filter_af3_ra_max = 5.0

// Forwarded to biophysical metrics (contact cutoff for interface).
params.negsteer_postprocess_contact_cutoff = 5.0

// ── Concurrency defaults ──────────────────────────────────────────────
params.max_af3_parallel = 30

// ---------------------------------------------------------------------------
// Includes
// ---------------------------------------------------------------------------

include { NEGSTEER_INTERFACE_METRICS   } from '../../modules/negsteer_interface_metrics'
include { EXTRACT_SURVIVOR_MANIFEST    } from '../../modules/negsteer_manifest'
include { AF3_SETUP_DB                 } from '../../modules/negsteer_af3_nomsa'
include { AF3_NOMSA_ON_SURVIVORS       } from '../../modules/negsteer_af3_nomsa'
include { AF3_PARSE_OUTPUT             } from '../../modules/negsteer_af3_nomsa'
include { NEGSTEER_BIOPHYSICAL_METRICS } from '../../modules/negsteer_biophysical_metrics'
include { NEGSTEER_ROSETTA_METRICS     } from '../../modules/negsteer_rosetta_metrics'
include { NEGSTEER_ORTHOGONAL_METRICS  } from '../../modules/negsteer_orthogonal_metrics'


// ---------------------------------------------------------------------------
// EXTRACT_SURVIVOR_MANIFEST is imported from modules/negsteer_manifest.nf
// (promoted out of this file by P0-34 so production main.nf and this test
// share one definition).  See the module for documentation.
// ---------------------------------------------------------------------------


// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------

workflow {

    log.info "Orthogonal-metrics test — input_dir: ${params.input_dir}"

    // ── Inputs ─────────────────────────────────────────────────────
    cross_csv_ch = Channel
        .fromPath("${params.input_dir}/cross_sequence_summary.csv",
                  checkIfExists: true)

    workdirs_ch = Channel
        .fromPath("${params.input_dir}/runs/*", type: 'dir', checkIfExists: true)
        .collect()

    // ── P0-29: interface-restricted metrics ────────────────────────
    NEGSTEER_INTERFACE_METRICS(
        cross_csv_ch,
        workdirs_ch,
        Channel.value(file("${projectDir}/bin/compute_interface_metrics.py"))
    )
    extended_csv_ch = NEGSTEER_INTERFACE_METRICS.out.extended_csv

    // ── Build per-survivor manifest ────────────────────────────────
    EXTRACT_SURVIVOR_MANIFEST(
        extended_csv_ch,
        workdirs_ch,
        Channel.value(file("${projectDir}/bin/extract_survivor_manifest.py"))
    )

    // ── Fan-out: one record per survivor ───────────────────────────
    // The manifest CSV is consumed via splitCsv.  Each record becomes
    // three downstream tuples — one per orthogonal-metrics stream.
    manifest_records_ch = EXTRACT_SURVIVOR_MANIFEST.out.manifest
        .splitCsv(header: true)
        .map { r ->
            [
                r.seq_name,
                r.canonical_pdb_abs,
                r.ground_truth_abs,
                r.effector_template_cif_abs,
                r.receptor_seq,
                r.effector_seq,
            ]
        }

    // ── AF3-no-MSA stream ──────────────────────────────────────────
    AF3_SETUP_DB()

    af3_input_ch = manifest_records_ch.map { rec ->
        tuple(
            rec[0],              // seq_name
            rec[4],              // receptor_seq
            rec[5],              // effector_seq
            file(rec[3]),        // effector_template_cif path
            file(rec[2]),        // ground_truth_pdb path
        )
    }

    AF3_NOMSA_ON_SURVIVORS(
        af3_input_ch,
        AF3_SETUP_DB.out.db_dir,
        AF3_SETUP_DB.out.ready_flag,
    )
    AF3_PARSE_OUTPUT(
        AF3_NOMSA_ON_SURVIVORS.out.prediction,
        Channel.value(file("${projectDir}/bin/parse_af3_output.py")),
    )

    // ── Biophysical stream ─────────────────────────────────────────
    biophysical_input_ch = manifest_records_ch.map { rec ->
        tuple(
            rec[0],              // seq_name
            file(rec[1]),        // canonical_pdb
            file(rec[2]),        // ground_truth_pdb (unused but kept symmetric)
        )
    }
    NEGSTEER_BIOPHYSICAL_METRICS(
        biophysical_input_ch,
        Channel.value(file("${projectDir}/bin/run_biophysical_metrics.py"))
    )

    // ── Rosetta stream ─────────────────────────────────────────────
    rosetta_input_ch = manifest_records_ch.map { rec ->
        tuple(
            rec[0],              // seq_name
            file(rec[1]),        // canonical_pdb
        )
    }
    // FastRelax XML — this test was previously calling
    // NEGSTEER_ROSETTA_METRICS with only the per-survivor tuple,
    // which had been silently inconsistent with the production
    // module's two-input signature (tuple + fastrelax_xml).  Adding
    // fastrelax_xml here as a side effect of the cache-busting
    // refactor; the test now matches main.nf's invocation shape.
    fastrelax_xml_ch = Channel.value(
        file("${projectDir}/bin/fastrelax_for_ia.xml")
    )
    NEGSTEER_ROSETTA_METRICS(
        rosetta_input_ch,
        fastrelax_xml_ch,
        Channel.value(file("${projectDir}/bin/run_rosetta_metrics.py"))
    )

    // ── Merge the three streams into the final CSV ─────────────────
    NEGSTEER_ORTHOGONAL_METRICS(
        extended_csv_ch,
        AF3_PARSE_OUTPUT.out.summary_csv.collect(),
        NEGSTEER_BIOPHYSICAL_METRICS.out.summary_csv.collect(),
        NEGSTEER_ROSETTA_METRICS.out.summary_csv.collect(),
        Channel.value(file("${projectDir}/bin/merge_orthogonal_metrics.py"))
    )
}

// ---------------------------------------------------------------------------
// On completion
// ---------------------------------------------------------------------------

workflow.onComplete {
    log.info """
    =============================================================
    Orthogonal-Metrics Test Complete
    =============================================================
    Output:   ${params.outdir}
    Duration: ${workflow.duration}
    Success:  ${workflow.success}
    =============================================================
    """.stripIndent()
}

workflow.onError {
    log.error "Orthogonal-metrics test failed: ${workflow.errorMessage}"
}
