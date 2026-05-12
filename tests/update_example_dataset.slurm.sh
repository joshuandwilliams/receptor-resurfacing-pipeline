#!/bin/bash
#SBATCH --job-name="update_examples"
#SBATCH -p jic-short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=update_examples_%j.out
#SBATCH --error=update_examples_%j.err
#SBATCH --chdir=/hpc-home/jowillia/receptor_design/receptor-resurfacing-pipeline/tests
# -----------------------------------------------------------------------------
# update_example_dataset.slurm.sh
# -----------------------------------------------------------------------------
# Replace tests/<module>/example_output_files/ with the subset of a freshly-
# run output folder that matches tests/<module>/reference_manifest.txt.
#
# Usage:
#   sbatch tests/update_example_dataset.slurm.sh \
#         --module <module> \
#         --updated-output-folder <path-to-receptor_resurfacing_results>
#
# Example:
#   sbatch tests/update_example_dataset.slurm.sh \
#         --module rfdiffusion \
#         --updated-output-folder tests/rfdiffusion/receptor_resurfacing_results
#
# Notes:
# - The script runs non-interactively (no confirmation prompt).  Diff is
#   printed to stdout so you can review in slurm_<job-id>.out / .err.
# - A manifest pattern that matches ZERO files in the updated output
#   folder is a hard error: either the test run didn't produce that
#   output, or the manifest is stale.  Fix one before re-running.
# - To preview without writing, run locally on Mac:
#       python3 tests/_update_example_dataset_impl.py \
#           --module rfdiffusion --updated-output-folder <path> --dry-run
# - After this job completes, sync example_output_files/ back to Mac via
#   scripts/sync_from_hpc.sh, then commit + push.
# -----------------------------------------------------------------------------

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IMPL="${SCRIPT_DIR}/_update_example_dataset_impl.py"

if [ ! -f "${IMPL}" ]; then
    echo "ERROR: implementation not found: ${IMPL}" >&2
    exit 2
fi

echo "============================================================"
echo "update_example_dataset.slurm.sh"
echo "  Date: $(date)"
echo "  Node: $(hostname)"
echo "  CWD:  $(pwd)"
echo "  Args: $*"
echo "============================================================"

python3 "${IMPL}" "$@"

echo ""
echo "============================================================"
echo "Update complete: $(date)"
echo "============================================================"
echo "Next step (on Mac):"
echo "  ./scripts/sync_from_hpc.sh"
echo "  git diff --stat tests/<module>/example_output_files/"
echo "  git add tests/<module>/example_output_files/"
echo "  git commit -m \"tests/<module>: regenerate example_output_files\""
echo "  git push"
echo "  ./scripts/sync_to_hpc.sh"
