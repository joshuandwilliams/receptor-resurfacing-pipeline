#!/bin/bash
# -----------------------------------------------------------------------------
# scripts/refresh_orthog_fixture.sh
# -----------------------------------------------------------------------------
# Refresh the orthogonal_metrics test's input fixture from the latest
# negative_steering test outputs.  Mirrors what would happen in production:
# the negsteer cohort emerges from the negsteer pipeline, then feeds the
# orthogonal_metrics processing.
#
# After this script runs, tests/orthogonal_metrics/data/negsteer_run/ will
# contain every steered design's workdir (and the matching
# cross_sequence_summary.csv) from the most recent negsteer test run.  The
# orthogonal_metrics test will then exercise the full cohort instead of the
# previously-staged 2-survivor mini-fixture.
#
# Stale absolute paths embedded in plan.json (e.g. /Users/... vs /hpc-home/...)
# are handled at runtime by extract_survivor_manifest.py's remap logic — no
# pre-processing step is needed.
#
# Usage:
#   ./scripts/refresh_orthog_fixture.sh           # real refresh
#   ./scripts/refresh_orthog_fixture.sh --dry     # dry-run, no changes
#
# Pre-requisites: a completed negative_steering test run.  The script reads
# from tests/negative_steering/receptor_resurfacing_results/.
# -----------------------------------------------------------------------------

set -euo pipefail

trap 'rc=$?; printf >&2 "\nERROR: %s exited %d at line %d:\n  %s\n" \
    "$(basename "$0")" "$rc" "${LINENO}" "${BASH_COMMAND}"' ERR

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

SRC_RESULTS="${REPO_ROOT}/tests/negative_steering/receptor_resurfacing_results/negative_steering"
DST_FIXTURE="${REPO_ROOT}/tests/orthogonal_metrics/data/negsteer_run"

DRY=""
if [ "${1:-}" = "--dry" ]; then
    DRY="--dry-run"
    echo "=== DRY RUN — no files will be changed ==="
fi

if [ ! -d "${SRC_RESULTS}/runs" ] || [ ! -f "${SRC_RESULTS}/cross_sequence_summary.csv" ]; then
    echo "ERROR: negsteer test outputs not found at ${SRC_RESULTS}" >&2
    echo "       run ./tests/run_tests.sh --modules negative_steering --with-plots first" >&2
    exit 2
fi

echo "Source:      ${SRC_RESULTS}"
echo "Destination: ${DST_FIXTURE}"
echo

# Copy the CSV.  Overwrites; the previous fixture's 2-row CSV had the same
# structure as the 8-row negsteer output.
echo "→ cross_sequence_summary.csv"
rsync -a $DRY \
    "${SRC_RESULTS}/cross_sequence_summary.csv" \
    "${DST_FIXTURE}/cross_sequence_summary.csv"

# Mirror the runs/ tree.  --delete so stale workdirs (e.g. controls that
# aren't relevant to orthogonal_metrics) don't accumulate.
echo "→ runs/ (mirror with --delete)"
rsync -a --delete $DRY \
    "${SRC_RESULTS}/runs/" \
    "${DST_FIXTURE}/runs/"

if [ -z "$DRY" ]; then
    echo
    echo "Refreshed.  New fixture cohort:"
    awk -F, 'NR==1{for(i=1;i<=NF;i++){h[$i]=i}; next} {
        printf "  %-25s cross_tier=%-5s row_type=%s\n", \
            $h["mpnn_sequence"], $h["cross_tier"], $h["row_type"]
    }' "${DST_FIXTURE}/cross_sequence_summary.csv"
    echo
    echo "Run the orthogonal_metrics test against the refreshed fixture:"
    echo "  ./tests/run_tests.sh --modules orthogonal_metrics --with-plots"
fi
