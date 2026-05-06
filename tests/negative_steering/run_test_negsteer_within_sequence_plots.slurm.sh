#!/bin/bash
#SBATCH --job-name="negsteer_within_seq_test"
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
#SBATCH --chdir=/hpc-home/jowillia/receptor_design/receptor-resurfacing-pipeline/tests/negative_steering

# =============================================================================
# run_test_negsteer_within_sequence_plots_slurm.sh
# -----------------------------------------------------------------------------
# Re-renders the within-sequence (per-seed) negsteer plots from a previous
# test_negative_steering run's per-sequence raw_per_seed_results.csv files.
# Companion to run_test_negsteer_plots_slurm.sh — that wrapper handles the
# cohort-level (one-row-per-sequence) plots; this wrapper handles the
# within-sequence (one-row-per-seed) plots.
#
# Outputs (under <outdir>):
#   negsteer_per_seed_dispersion_overview.png
#       single scatter, one point per (MPNN sequence, Boltz seed) of the
#       representative sg.  Shape = stage that produced the final
#       prediction (cold-start / steering / reversion); colour = tier.
#   negsteer_per_seed_dispersion_grid.png
#       small-multiples grid: one panel per MPNN sequence, sorted by
#       composite descending.  Shows ALL seeds across all sgs in grey,
#       representative-sg seeds highlighted green (pass) or red (fail).
#       Vertical threshold line at ra_eff = 5 Å.
#   negsteer_per_design_stage_trajectories.png
#       one panel per RFDiffusion design.  Line per (sequence, sg, seed)
#       showing ra_eff at each stage (cold-start → steered → reverted).
#       Marker on the FINAL stage.  Lines coloured by MPNN sequence
#       number (so seq 3 is the same colour wherever it appears).
#   negsteer_weighted_vs_true_jaccard.png
#       single scatter, weighted vs true jaccard per representative seed.
#       y = x diagonal — points below indicate weighted understates,
#       points above indicate weighted overstates.
#   negsteer_composite_vs_confidence.png
#       multi-panel scatter: composite (y) vs each confidence metric (x).
#       One point per MPNN sequence (median across rep seeds).
#       Bounded [0, 1] metrics shown on a fixed 0-1 axis; pLDDT shown
#       on its native 0-100 scale; ipae and pae_mean autoscale.
#
# Usage:
#   sbatch tests/negative_steering/run_test_negsteer_within_sequence_plots_slurm.sh
#
# Optional CLI flags forwarded to the python script (override defaults):
#   --runs-dir          PATH   per-sequence runs directory
#   --cross-summary-csv PATH   cross_sequence_summary.csv path
#   --outdir            PATH   output dir
# =============================================================================

set -euo pipefail

# ── Paths ─────────────────────────────────────────────────────────────────
PIPELINE_DIR="/hpc-home/jowillia/receptor_design/receptor-resurfacing-pipeline"
TEST_DIR="${PIPELINE_DIR}/tests/negative_steering"

RESULTS_DIR="${TEST_DIR}/receptor_resurfacing_results"
RUNS_DIR="${RESULTS_DIR}/negative_steering/runs"
CROSS_CSV="${RESULTS_DIR}/negative_steering/cross_sequence_summary.csv"

OUTDIR="${RESULTS_DIR}/plots_iter"

TEST_SCRIPT="${TEST_DIR}/test_negsteer_within_sequence_plots.py"

# ── Container ─────────────────────────────────────────────────────────────
# Same image as the cohort-level plots — matplotlib + numpy only.
RFDIFF_CONTAINER="/hpc-home/jowillia/singularity/HADDOCK_RFDiffusion_ProteinMPNN_MMseqs2/HADDOCK_RFDiffusion_ProteinMPNN_MMseqs2.img"

# ── Sanity checks ─────────────────────────────────────────────────────────
if [ ! -f "${TEST_SCRIPT}" ]; then
    echo "ERROR: test script not found: ${TEST_SCRIPT}" >&2
    exit 1
fi

if [ ! -d "${RUNS_DIR}" ]; then
    echo "ERROR: runs/ directory not found: ${RUNS_DIR}" >&2
    echo "       Run test_negative_steering first to generate it." >&2
    exit 1
fi

mkdir -p "${OUTDIR}"

# Cross-summary CSV is now REQUIRED (the new plots use it for the
# representative sequence_group lookup and for composite-score sorting).
if [ ! -f "${CROSS_CSV}" ]; then
    echo "ERROR: cross_sequence_summary.csv not found: ${CROSS_CSV}" >&2
    echo "       Run test_negative_steering first to generate it." >&2
    exit 1
fi

# ── Launch ────────────────────────────────────────────────────────────────
echo "============================================================"
echo "Negative-Steering Within-Sequence Plots — variant test"
echo "============================================================"
echo "Pipeline dir:  ${PIPELINE_DIR}"
echo "Test dir:      ${TEST_DIR}"
echo "Runs dir:      ${RUNS_DIR}"
echo "Cross CSV:     ${CROSS_CSV}"
echo "Output dir:    ${OUTDIR}"
echo "Date:          $(date)"
echo "Node:          $(hostname)"
echo "============================================================"

singularity exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --env MPLCONFIGDIR=/tmp \
    "${RFDIFF_CONTAINER}" \
    python "${TEST_SCRIPT}" \
        --runs-dir          "${RUNS_DIR}" \
        --cross-summary-csv "${CROSS_CSV}" \
        --outdir            "${OUTDIR}" \
        "$@"

echo ""
echo "============================================================"
echo "Plots written to: ${OUTDIR}"
echo "Finished: $(date)"
echo "============================================================"
