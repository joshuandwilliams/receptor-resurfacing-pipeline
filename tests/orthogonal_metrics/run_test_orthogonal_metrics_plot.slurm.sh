#!/bin/bash
#SBATCH --job-name="orthog_plots_test"
#SBATCH -p jic-medium
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=test_plots_%j.out
#SBATCH --error=test_plots_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jowillia@nbi.ac.uk
#SBATCH --chdir=/hpc-home/jowillia/receptor_design/receptor-resurfacing-pipeline/tests/orthogonal_metrics

# =============================================================================
# run_test_orthogonal_metrics_plots_slurm.sh
# -----------------------------------------------------------------------------
# Re-renders the orthogonal-metrics plots from a previous test_orthogonal_metrics
# run's survivors_with_orthogonal_metrics.csv joined to the upstream negsteer
# cross_sequence_summary.csv, without invoking Nextflow.  Designed for fast
# iteration of the plot code itself against cached test outputs.
#
# As of 2026-04-28 BOTH csvs are required — the metrics-vs-composite scatter
# and the combined cohort+orthogonal summary both join on mpnn_sequence to
# source upstream composite scores and cohort-wide context.
#
# Outputs (under <outdir>):
#   orthogonal_af3_vs_boltz.png             cross-model ra_eff agreement
#                                           (the headline orthogonal plot).
#   orthogonal_filter_cascade.png           waterfall through the GATING
#                                           cascade (sc → bsa → interface_plddt),
#                                           pre-filtered to rows that actually
#                                           entered the orthogonal stage.
#                                           AF3 disagreement appears as a
#                                           trailing yellow bar (informational,
#                                           not gating).
#   orthogonal_metrics_vs_composite.png     2x3 scatter: composite (Y) vs
#                                           each orthogonal metric (X) — Sc,
#                                           BSA, interface_plddt, ΔΔG,
#                                           AF3 ra_eff (one slot empty).
#                                           Replaces the old strip-plot
#                                           distributions and standalone
#                                           ΔΔG plot.
#   orthogonal_combined_cohort_summary.png  the big combined view: one row
#                                           per cross_summary sequence (steered
#                                           only), sorted by composite desc.
#                                           Tier as left-edge stripe; columns
#                                           span both negsteer + orthogonal
#                                           cascades; missing orthogonal
#                                           cells grey-hatched.
#
# Usage:
#   sbatch tests/orthogonal_metrics/run_test_orthogonal_metrics_plots_slurm.sh
#
# Optional CLI flags forwarded to the python script (override defaults):
#   --survivors-csv     PATH   override survivors_with_orthogonal_metrics.csv
#   --cross-summary-csv PATH   override upstream cross_sequence_summary.csv
#   --outdir            PATH   override output dir
# =============================================================================

set -euo pipefail

# ── Paths ─────────────────────────────────────────────────────────────────
PIPELINE_DIR="/hpc-home/jowillia/receptor_design/receptor-resurfacing-pipeline"
TEST_DIR="${PIPELINE_DIR}/tests/orthogonal_metrics"

# nextflow.config overrides params.outdir, so the real results dir is
# ${projectDir}/${params.project_name}_results = receptor_resurfacing_results.
RESULTS_DIR="${TEST_DIR}/receptor_resurfacing_results"
SURVIVORS_CSV="${RESULTS_DIR}/orthogonal_metrics/survivors_with_orthogonal_metrics.csv"

# Upstream negsteer context — REQUIRED.  Carries cross_tier, composite
# score, mutation strings, and the full upstream cohort (including
# sequences that didn't reach the orthogonal stage, which the combined
# summary plot needs to render with grey-hatched cells).
CROSS_SUMMARY_CSV="${PIPELINE_DIR}/tests/negative_steering/receptor_resurfacing_results/negative_steering/cross_sequence_summary.csv"

OUTDIR="${RESULTS_DIR}/plots_iter"

TEST_SCRIPT="${TEST_DIR}/test_orthogonal_metrics_plots.py"

# ── Container ─────────────────────────────────────────────────────────────
# Same image as production ORTHOGONAL_METRICS_PLOTS process.  Plotting
# only depends on matplotlib + numpy + stdlib csv — no GPU needed.
RFDIFF_CONTAINER="${PIPELINE_DIR}/containers/HADDOCK_RFDiffusion_ProteinMPNN_MMseqs2.img"

# ── Sanity checks ─────────────────────────────────────────────────────────
if [ ! -f "${TEST_SCRIPT}" ]; then
    echo "ERROR: test script not found: ${TEST_SCRIPT}" >&2
    exit 1
fi

if [ ! -f "${SURVIVORS_CSV}" ]; then
    echo "ERROR: survivors_with_orthogonal_metrics.csv not found: ${SURVIVORS_CSV}" >&2
    echo "       Run test_orthogonal_metrics first to generate it." >&2
    exit 1
fi

if [ ! -f "${CROSS_SUMMARY_CSV}" ]; then
    echo "ERROR: cross_sequence_summary.csv not found: ${CROSS_SUMMARY_CSV}" >&2
    echo "       The metrics-vs-composite and combined-summary plots both" >&2
    echo "       require it to source composite scores and the upstream" >&2
    echo "       cohort.  Run test_negative_steering first." >&2
    exit 1
fi

mkdir -p "${OUTDIR}"

# ── Launch ────────────────────────────────────────────────────────────────
echo "============================================================"
echo "Orthogonal-Metrics Plots Iteration — variant test"
echo "============================================================"
echo "Pipeline dir:      ${PIPELINE_DIR}"
echo "Test dir:          ${TEST_DIR}"
echo "Survivors CSV:     ${SURVIVORS_CSV}"
echo "Cross-summary CSV: ${CROSS_SUMMARY_CSV}"
echo "Output dir:        ${OUTDIR}"
echo "Date:              $(date)"
echo "Node:              $(hostname)"
echo "============================================================"

singularity exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --env MPLCONFIGDIR=/tmp \
    "${RFDIFF_CONTAINER}" \
    python "${TEST_SCRIPT}" \
        --survivors-csv     "${SURVIVORS_CSV}" \
        --cross-summary-csv "${CROSS_SUMMARY_CSV}" \
        --outdir            "${OUTDIR}" \
        "$@"

echo ""
echo "============================================================"
echo "Plots written to: ${OUTDIR}"
echo "Finished: $(date)"
echo "============================================================"
