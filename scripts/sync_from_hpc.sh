#!/bin/bash
#
# Sync test fixtures back from HPC to Mac (the reverse of sync_to_hpc.sh).
#
# Use this after running tests/update_example_dataset.slurm.sh on HPC.  The
# update script writes the new tests/<module>/example_output_files/ on HPC;
# this script pulls those updated fixtures back to Mac so they can be
# committed and pushed.
#
# Usage:
#   ./scripts/sync_from_hpc.sh                      # all modules, real sync
#   ./scripts/sync_from_hpc.sh --dry                # all modules, dry-run
#   ./scripts/sync_from_hpc.sh --module rfdiffusion # single module
#   ./scripts/sync_from_hpc.sh --module negative_steering --dry
#
# Transport:
#   rsync over SSH using the 'slurm' host alias defined in ~/.ssh/config.
#
# Scope:
#   Pulls tests/<module>/example_output_files/ only.  Source-tree code,
#   data/, run artefacts (work/, results/, *.out, *.err) are NOT pulled —
#   those belong on Mac as authoritative.
#
# Uses --delete on the per-module example_output_files/ subtree so files
# removed by the manifest update on HPC are also removed on Mac.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# Site settings: the ssh host alias and the path on the cluster are specific
# to one person's setup, so they come from site.conf rather than from here.
_SELF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=scripts/site_conf.sh
source "${_SELF_DIR}/site_conf.sh"

HPC_USER_HOST="${RECEPTOR_HPC_HOST:-slurm}"
HPC_BASE="${RECEPTOR_HPC_DEST:-receptor_design/receptor-resurfacing-pipeline}"
HPC_BASE="${HPC_BASE%/}"   # trailing slash is optional in site.conf

VALID_MODULES=(
    rfdiffusion
    proteinmpnn
    rosetta_filtering
    negative_steering
    orthogonal_metrics
    pose_solver
)

DRY_RUN=""
SELECTED=""

usage() {
    cat <<EOF
Usage: $(basename "$0") [--module MOD] [--dry]

  --module MOD   Sync only the named module (one of: ${VALID_MODULES[*]}).
                 Omit to sync all modules.
  --dry          rsync dry-run; show what would change without copying.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry|--dry-run) DRY_RUN="--dry-run"; shift ;;
        --module)        SELECTED="$2"; shift 2 ;;
        -h|--help)       usage; exit 0 ;;
        *)               echo "ERROR: unknown argument: $1" >&2; usage; exit 2 ;;
    esac
done

if [ -n "${SELECTED}" ]; then
    if [[ ! " ${VALID_MODULES[*]} " == *" ${SELECTED} "* ]]; then
        echo "ERROR: unknown module: ${SELECTED}" >&2
        echo "       valid: ${VALID_MODULES[*]}" >&2
        exit 2
    fi
    MODULES=("${SELECTED}")
else
    MODULES=("${VALID_MODULES[@]}")
fi

cd "${REPO_ROOT}"

if [ -n "${DRY_RUN}" ]; then
    echo "=== DRY RUN — no files will be changed ==="
fi

for m in "${MODULES[@]}"; do
    REL="tests/${m}/example_output_files"
    SRC="${HPC_USER_HOST}:${HPC_BASE}/${REL}/"
    DST="${REPO_ROOT}/${REL}/"

    # Skip silently if the module has never had a fixture.
    # Probe by listing the remote dir; rsync would also handle missing
    # source but we want a friendly skip message.
    if ! ssh -q -o BatchMode=yes "${HPC_USER_HOST}" \
            "[ -d ${HPC_BASE}/${REL} ]" 2>/dev/null; then
        echo "[${m}] SKIP: no remote fixture at ${HPC_BASE}/${REL}"
        continue
    fi

    mkdir -p "${DST}"
    echo "[${m}] rsync ${SRC} → ${DST}"
    rsync -av --delete ${DRY_RUN} -e ssh \
        --exclude='.DS_Store' \
        "${SRC}" "${DST}"
    echo ""
done

echo "Sync from HPC complete."
if [ -n "${DRY_RUN}" ]; then
    echo "(dry-run — nothing was actually pulled)"
else
    echo "Next: review with 'git diff --stat tests/', then commit + push."
fi
