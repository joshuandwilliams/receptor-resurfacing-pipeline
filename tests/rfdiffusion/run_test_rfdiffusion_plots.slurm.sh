#!/bin/bash
#SBATCH --job-name="rfdiff_plots_test"
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
#SBATCH --chdir=/hpc-home/jowillia/receptor_design/receptor-resurfacing-pipeline/tests/rfdiffusion

# =============================================================================
# run_test_rfdiffusion_plots_slurm.sh
# -----------------------------------------------------------------------------
# Re-renders the RFDiffusion plots from a previous test_rfdiffusion run's
# rfdiffusion_metrics.json without invoking Nextflow or the GPU.
#
# The five candidate dendrogram-placement variants for the two problem
# plots (contact_map and design_clustering) are dropped under
#   <outdir>/contactmap_variants/
#   <outdir>/clustering_variants/
# so all five can be eyeballed side-by-side.
#
# The three already-good plots (specificity_coverage, design_lengths,
# com_displacement) are reproduced from the production rfdiffusion_plots.py
# unchanged for reference.
#
# Usage:
#   sbatch tests/rfdiffusion/run_test_rfdiffusion_plots_slurm.sh
#
# Optional CLI flags forwarded to the python script:
#   --metrics PATH   override the default metrics JSON path
#   --outdir PATH    override the default output dir
#   --variants ...   subset of B,C,F,G,H to render
# =============================================================================

set -euo pipefail

# ── Paths ─────────────────────────────────────────────────────────────────
PIPELINE_DIR="/hpc-home/jowillia/receptor_design/receptor-resurfacing-pipeline"
TEST_DIR="${PIPELINE_DIR}/tests/rfdiffusion"

# Default metrics location: the canonical test_rfdiffusion output.
METRICS="${TEST_DIR}/receptor_resurfacing_results/rfdiffusion/rfdiffusion_metrics.json"

# Default outdir: a sibling of the existing plots/ dir so we don't
# clobber the production-test plot output.
OUTDIR="${TEST_DIR}/receptor_resurfacing_results/plots_iter"

# Production plotting script (used for the three "good" plots).
PROD_SCRIPT="${PIPELINE_DIR}/bin/rfdiffusion_plots.py"

# Test-side python script (this run's iteration target).
TEST_SCRIPT="${TEST_DIR}/test_rfdiffusion_plots.py"

# ── Container ─────────────────────────────────────────────────────────────
# Reuse the rfdiff container — has matplotlib + scipy + numpy.  Uses
# the same image as the production RFDIFFUSION_PLOTS process (kept in
# sync via nextflow.config, see params.rfdiff_container).
RFDIFF_CONTAINER="/hpc-home/jowillia/singularity/HADDOCK_RFDiffusion_ProteinMPNN_MMseqs2/LRR_Pipeline.img"

# ── Sanity checks ─────────────────────────────────────────────────────────
if [ ! -f "${TEST_SCRIPT}" ]; then
    echo "ERROR: test script not found: ${TEST_SCRIPT}" >&2
    exit 1
fi
if [ ! -f "${PROD_SCRIPT}" ]; then
    echo "WARNING: production plot script not found: ${PROD_SCRIPT}" >&2
    echo "  → variants will still render, but the three good plots will be skipped." >&2
    PROD_FLAG=""
else
    PROD_FLAG="--prod-script ${PROD_SCRIPT}"
fi

mkdir -p "${OUTDIR}"

# ── Launch ────────────────────────────────────────────────────────────────
echo "============================================================"
echo "RFDiffusion Plots Iteration — variant test"
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
    "${RFDIFF_CONTAINER}" \
    python "${TEST_SCRIPT}" \
        --metrics "${METRICS}" \
        --outdir  "${OUTDIR}" \
        ${PROD_FLAG} \
        "$@"

echo ""
echo "============================================================"
echo "Plots written to: ${OUTDIR}"
echo "  Contact-map variants:   ${OUTDIR}/contactmap_variants/"
echo "  Clustering variants:    ${OUTDIR}/clustering_variants/"
echo "  Other plots (good):     ${OUTDIR}/"
echo "Finished: $(date)"
echo "============================================================"
