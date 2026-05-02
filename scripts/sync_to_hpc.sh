#!/bin/bash
#
# Sync the local repo to HPC via the mounted /Volumes/HPC-Home.
#
# Usage:
#   ./scripts/sync_to_hpc.sh           # real sync
#   ./scripts/sync_to_hpc.sh --dry     # dry run, no changes
#
# Excludes:
#   - .git/, Python/IDE caches, macOS metadata, Word lockfiles
#   - run artefacts (work/, tmp/, results/, .nextflow/)
#   - slurm log files in tests/*/
#
# Includes (despite being run artefacts in some sense):
#   - tests/full_test_run/example_output_files/  -- curated reference set
#
# The sync uses --delete, so files removed from the local repo are
# also removed from HPC. Always run with --dry first when in doubt.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
HPC_DEST="/Volumes/HPC-Home/receptor_design/receptor-resurfacing-pipeline/"

if [ ! -d "/Volumes/HPC-Home" ]; then
    echo "ERROR: /Volumes/HPC-Home is not mounted." >&2
    echo "Mount the HPC home directory first, then re-run." >&2
    exit 1
fi

DRY_RUN=""
if [ "${1:-}" = "--dry" ]; then
    DRY_RUN="--dry-run"
    echo "=== DRY RUN -- no files will be changed ==="
fi

cd "$REPO_ROOT"

rsync -av --delete $DRY_RUN \
    --exclude='.git/' \
    --exclude='__pycache__/' \
    --exclude='.pytest_cache/' \
    --exclude='.DS_Store' \
    --exclude='*.pyc' \
    --exclude='.ruff_cache/' \
    --exclude='.mypy_cache/' \
    --exclude='nxf_home/' \
    --exclude='~$*' \
    --exclude='tests/*/work/' \
    --exclude='tests/*/tmp/' \
    --exclude='tests/*/results/' \
    --exclude='tests/*/receptor_resurfacing_results/' \
    --exclude='tests/*/.nextflow*' \
    --exclude='tests/*/*.out' \
    --exclude='tests/*/*.err' \
    --exclude='tests/full_test_run/results/' \
    --exclude='tests/full_test_run/work/' \
    --exclude='tests/full_test_run/tmp/' \
    --exclude='tests/full_test_run/.nextflow*' \
    --exclude='rsync_dryrun*.txt' \
    --exclude='rsync_deletions*.txt' \
    "$REPO_ROOT/" "$HPC_DEST"

echo
echo "Sync complete."
