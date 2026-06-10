#!/bin/bash
#SBATCH --job-name="rosetta_plots_test"
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
#SBATCH --chdir=/hpc-home/jowillia/receptor_design/receptor-resurfacing-pipeline/tests/rosetta_filtering

# =============================================================================
# run_test_rosetta_filtering_plots_slurm.sh
# -----------------------------------------------------------------------------
# Re-renders the Rosetta filtering plot from a previous
# test_rosetta_filtering run's rosetta_filter_metrics.json without
# invoking Nextflow.  Designed for fast iteration of the plot code
# itself against cached test outputs.
#
# Output: <outdir>/rosetta_sc_histogram.png
#
# Usage:
#   sbatch tests/rosetta_filtering/run_test_rosetta_filtering_plots_slurm.sh
#
# Optional CLI flags forwarded to the python script:
#   --metrics PATH   override the default metrics JSON path
#   --outdir PATH    override the default output dir
# =============================================================================

set -euo pipefail

# ── Paths ─────────────────────────────────────────────────────────────────
PIPELINE_DIR="/hpc-home/jowillia/receptor_design/receptor-resurfacing-pipeline"
TEST_DIR="${PIPELINE_DIR}/tests/rosetta_filtering"

# nextflow.config overrides params.outdir, so the real results dir is
# ${projectDir}/${params.project_name}_results = receptor_resurfacing_results.
RESULTS_DIR="${TEST_DIR}/receptor_resurfacing_results"
METRICS="${RESULTS_DIR}/rosetta_filtering/rosetta_filter_metrics.json"
OUTDIR="${RESULTS_DIR}/plots_iter"

TEST_SCRIPT="${TEST_DIR}/test_rosetta_filtering_plots.py"

# ── Container ─────────────────────────────────────────────────────────────
# Same image as production ROSETTA_FILTER_PLOTS process — kept in sync
# via nextflow.config (params.rosetta_container).
ROSETTA_CONTAINER="${PIPELINE_DIR}/containers/Rosetta.img"

# ── Sanity checks ─────────────────────────────────────────────────────────
if [ ! -f "${TEST_SCRIPT}" ]; then
    echo "ERROR: test script not found: ${TEST_SCRIPT}" >&2
    exit 1
fi

mkdir -p "${OUTDIR}"

# ── Launch ────────────────────────────────────────────────────────────────
echo "============================================================"
echo "Rosetta Filtering Plots Iteration — variant test"
echo "============================================================"
echo "Pipeline dir:  ${PIPELINE_DIR}"
echo "Test dir:      ${TEST_DIR}"
echo "Metrics:       ${METRICS}"
echo "Output dir:    ${OUTDIR}"
echo "Date:          $(date)"
echo "Node:          $(hostname)"
echo "============================================================"

singularity exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --env MPLCONFIGDIR=/tmp \
    "${ROSETTA_CONTAINER}" \
    python "${TEST_SCRIPT}" \
        --metrics "${METRICS}" \
        --outdir  "${OUTDIR}" \
        "$@"

echo ""
echo "============================================================"
echo "Plot written to: ${OUTDIR}/rosetta_sc_histogram.png"
echo "Finished: $(date)"
echo "============================================================"
