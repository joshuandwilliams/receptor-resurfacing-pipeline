#!/bin/bash
#SBATCH --job-name="mpnn_plots_test"
# No -p: the partition is site-specific.  tests/run_tests.sh passes
# RECEPTOR_PARTITION_CPU from site.conf as `sbatch -p`, where it beats both
# the environment and a header line.  Submitting this script directly takes
# the cluster's default partition.
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=test_plots_%j.out
#SBATCH --error=test_plots_%j.err
#SBATCH --mail-type=END,FAIL
# No --mail-user: sbatch(1) mails the submitting user when only --mail-type
# is given, so an address here would be site-specific and redundant.
# No --chdir: the path is site-specific.  tests/run_tests.sh passes it as
# `sbatch --chdir`, which is what sets the $PWD this script resolves its
# inputs from.  Submit from this directory if you run it by hand.

# =============================================================================
# run_test_proteinmpnn_plots_slurm.sh
# -----------------------------------------------------------------------------
# Re-renders the ProteinMPNN plots from a previous test_proteinmpnn run's
# scored_metadata.csv without invoking Nextflow.  Designed for fast
# iteration of the plot code itself against cached test outputs.
#
# Outputs (under <outdir>):
#   mpnn_score_distribution.png   (panels flipped, ranked by design-region)
#   mpnn_aa_composition.png       (with native reference overlay)
#   mpnn_sequence_diversity.png   (per-region pairwise distance heatmaps)
#   mpnn_physicochem.png          (NEW: hydrophobicity + net charge)
#
# Usage:
#   sbatch tests/proteinmpnn/run_test_proteinmpnn_plots_slurm.sh
#
# Optional CLI flags forwarded to the python script:
#   --metadata    PATH   override the default metadata CSV path
#   --outdir      PATH   override the default output dir
#   --cluster-csv PATH   (kept for CLI parity; unused by new diversity plot)
# =============================================================================

set -euo pipefail

# ── Paths ─────────────────────────────────────────────────────────────────
PIPELINE_DIR="${RECEPTOR_PIPELINE_DIR:-$(cd "../.." && pwd)}"
TEST_DIR="${PIPELINE_DIR}/tests/proteinmpnn"

# nextflow.config overrides params.outdir, so the real results dir is
# ${projectDir}/${params.project_name}_results = receptor_resurfacing_results.
RESULTS_DIR="${TEST_DIR}/receptor_resurfacing_results"
METADATA="${RESULTS_DIR}/sequences/scored_metadata.csv"
CLUSTER_CSV="${RESULTS_DIR}/sequences/mpnn_cluster_counts.csv"
OUTDIR="${RESULTS_DIR}/plots_iter"

TEST_SCRIPT="${TEST_DIR}/iterate_mpnn_plots.py"

# ── Container ─────────────────────────────────────────────────────────────
# Same image as production MPNN_PLOTS process — kept in sync via
# nextflow.config (params.rfdiff_container; ProteinMPNN ships in the
# same container as RFDiffusion for this pipeline).
RFDIFF_CONTAINER="${PIPELINE_DIR}/containers/HADDOCK_RFDiffusion_ProteinMPNN_MMseqs2.img"

# ── Sanity checks ─────────────────────────────────────────────────────────
if [ ! -f "${TEST_SCRIPT}" ]; then
    echo "ERROR: test script not found: ${TEST_SCRIPT}" >&2
    exit 1
fi

mkdir -p "${OUTDIR}"

# ── Launch ────────────────────────────────────────────────────────────────
echo "============================================================"
echo "ProteinMPNN Plots Iteration — variant test"
echo "============================================================"
echo "Pipeline dir:  ${PIPELINE_DIR}"
echo "Test dir:      ${TEST_DIR}"
echo "Metadata:      ${METADATA}"
echo "Cluster CSV:   ${CLUSTER_CSV}"
echo "Output dir:    ${OUTDIR}"
echo "Date:          $(date)"
echo "Node:          $(hostname)"
echo "============================================================"

singularity exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --env MPLCONFIGDIR=/tmp \
    "${RFDIFF_CONTAINER}" \
    python "${TEST_SCRIPT}" \
        --metadata    "${METADATA}" \
        --cluster-csv "${CLUSTER_CSV}" \
        --outdir      "${OUTDIR}" \
        "$@"

echo ""
echo "============================================================"
echo "Plots written to: ${OUTDIR}"
echo "Finished: $(date)"
echo "============================================================"
