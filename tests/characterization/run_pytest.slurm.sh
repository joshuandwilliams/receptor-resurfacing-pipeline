#!/bin/bash
#SBATCH --job-name="pytest_characterisation"
#SBATCH -p jic-short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jowillia@nbi.ac.uk
#SBATCH --chdir=/hpc-home/jowillia/receptor_design/receptor-resurfacing-pipeline/tests/characterization

# =============================================================================
# run_pytest.slurm.sh
# -----------------------------------------------------------------------------
# Runs the hpc-tier characterisation test suite against the per-module test
# outputs produced by the five run_test_<module>.slurm.sh scripts.
#
# Prerequisites:
#   - All five per-module tests have completed successfully, producing
#     tests/<module>/receptor_resurfacing_results/ on disk.
#   - tests/<module>/example_output_files/ (the committed reference sets)
#     are present in the repo.
#
# The test suite resolves all paths automatically via the stage_reference_root
# and stage_output_root fixtures in conftest.py — no environment variables
# needed.
#
# Usage:
#   sbatch tests/characterization/run_pytest.slurm.sh
#
# Optional: pass pytest flags through $@, e.g.:
#   sbatch tests/characterization/run_pytest.slurm.sh -- -v
# =============================================================================

set -euo pipefail

# ── Paths ─────────────────────────────────────────────────────────────────
PIPELINE_DIR="/hpc-home/jowillia/receptor_design/receptor-resurfacing-pipeline"
PYTEST_IMG="/hpc-home/jowillia/singularity/pytest/pytest_runner.img"

TEST_DIR="${PIPELINE_DIR}/tests/characterization"

# ── Launch ────────────────────────────────────────────────────────────────
echo "============================================================"
echo "Characterisation Test Suite — pytest"
echo "============================================================"
echo "Pipeline dir: ${PIPELINE_DIR}"
echo "Test dir:     ${TEST_DIR}"
echo "Image:        ${PYTEST_IMG}"
echo "Date:         $(date)"
echo "Node:         $(hostname)"
echo "============================================================"

singularity exec \
    --bind "${PIPELINE_DIR}:${PIPELINE_DIR}" \
    --env MPLCONFIGDIR=/tmp \
    "${PYTEST_IMG}" \
    python -m pytest \
        -m hpc \
        --tb=short \
        -q \
        --rootdir="${PIPELINE_DIR}" \
        "${TEST_DIR}" \
        "$@"

echo ""
echo "============================================================"
echo "pytest finished: $(date)"
echo "============================================================"
