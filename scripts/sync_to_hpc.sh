#!/bin/bash
#
# Sync the local repo to HPC via SSH (uses ~/.ssh/config alias 'slurm').
#
# Usage:
#   ./scripts/sync_to_hpc.sh           # real sync
#   ./scripts/sync_to_hpc.sh --dry     # dry run, no changes
#
# Transport:
#   rsync over SSH, using the 'slurm' host alias defined in
#   ~/.ssh/config.  Requires passwordless key auth -- the
#   id_ed25519_nbi key is mapped to slurm.nbi.ac.uk in ~/.ssh/config.
#
# Excludes:
#   - .git/, Python/IDE caches, macOS metadata, Word lockfiles
#   - run artefacts in tests/*/ and experiments/campaigns/*/runs/*/
#     (work/, tmp/, results/, .nextflow*, receptor_resurfacing_results/)
#   - slurm log files in tests/*/ and experiments/campaigns/*/runs/*/
#   - tests/curation_staging/ -- temporary tarball staging area used
#     during fixture curation; safe to recreate locally without syncing
#   - tests/full_test_run/reference_data_helpers/ -- ad-hoc helpers
#     pulled from HPC for path-coverage analysis
#
# Includes (despite being run artefacts in some sense):
#   - tests/full_test_run/example_output_files/  -- curated reference set
#     (will retire at phase 2.9)
#
# The sync uses --delete, so files removed from the local repo are
# also removed from HPC. Always run with --dry first when in doubt.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
HPC_USER_HOST="slurm"
HPC_DEST="receptor_design/receptor-resurfacing-pipeline/"

DRY_RUN=""
if [ "${1:-}" = "--dry" ]; then
    DRY_RUN="--dry-run"
    echo "=== DRY RUN -- no files will be changed ==="
fi

cd "$REPO_ROOT"

rsync -av --delete $DRY_RUN -e ssh \
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
    --exclude='experiments/campaigns/*/runs/*/results/' \
    --exclude='experiments/campaigns/*/runs/*/work/' \
    --exclude='experiments/campaigns/*/runs/*/tmp/' \
    --exclude='experiments/campaigns/*/runs/*/.nextflow*' \
    --exclude='experiments/campaigns/*/runs/*/*.out' \
    --exclude='experiments/campaigns/*/runs/*/*.err' \
    --exclude='tests/full_test_run/reference_data_helpers/' \
    --exclude='tests/curation_staging/' \
    --exclude='rsync_dryrun*.txt' \
    --exclude='rsync_deletions*.txt' \
    "$REPO_ROOT/" "${HPC_USER_HOST}:${HPC_DEST}"

echo
echo "Sync complete."
