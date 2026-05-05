#!/bin/bash
#SBATCH --job-name="cleanup_finished_runs"
#SBATCH -p jic-medium
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=1G
#SBATCH --time=00:30:00
#SBATCH --output=cleanup_%A_%a.out
#SBATCH --error=cleanup_%A_%a.err
#SBATCH --array=0-7
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jowillia@nbi.ac.uk
#
# clean_all_finished_runs.slurm.sh — array job that runs
# clean_finished_run.sh on each eligible campaign run identified by the
# storage bloat audit (notes/inventory/17_storage_bloat_audit.md §8).
#
# Submit from the pipeline root on the HPC login node:
#   sbatch scripts/clean_all_finished_runs.slurm.sh
#
# Each array task is independent; failure of one does not affect the
# others. To run a subset, override --array on the command line, e.g.:
#   sbatch --array=0,7 scripts/clean_all_finished_runs.slurm.sh
# To throttle concurrency (max N tasks at once), append %N:
#   sbatch --array=0-7%4 scripts/clean_all_finished_runs.slurm.sh
#
# Logs land in cwd as cleanup_<arrayJobId>_<taskId>.{out,err}.

set -euo pipefail

# SLURM copies the submitted script into /var/spool/slurmd/jobN/ before
# executing, so BASH_SOURCE[0] does NOT point to the original location in
# scripts/. Use SLURM_SUBMIT_DIR (set by SLURM to the directory sbatch was
# called from) as the pipeline root. Submit from the pipeline root, e.g.:
#   cd /hpc-home/jowillia/receptor_design/receptor-resurfacing-pipeline
#   sbatch scripts/clean_all_finished_runs.slurm.sh
PIPELINE_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
CLEAN="$PIPELINE_ROOT/scripts/clean_finished_run.sh"

if [[ ! -x "$CLEAN" ]]; then
    echo "ERROR: $CLEAN is not executable" >&2
    echo "  PIPELINE_ROOT=$PIPELINE_ROOT" >&2
    echo "  Submit from the pipeline root, or run: chmod +x $CLEAN" >&2
    exit 1
fi

# Eligible runs as of 2026-05-05 (see notes/inventory/17_storage_bloat_audit.md §8).
# Index in this array corresponds to the SLURM array task id.
# Status (OK / ERR) shown as a comment — the per-run cleaner re-checks
# .nextflow/history at runtime so an out-of-date comment is harmless.
RUNS=(
    "experiments/campaigns/pikp1_avrpia/runs/v1_4a"     # 0  OK,  62 negsteer designs
    "experiments/campaigns/pikp1_avrpia/runs/v1_4b"     # 1  OK, 102 negsteer designs (representative)
    "experiments/campaigns/pikp1_avrpikf/runs/v1_4a"    # 2  OK, 118 negsteer designs
    "experiments/campaigns/pikp1_pby2/runs/v1_4a"       # 3  OK, 114 negsteer designs
    "experiments/campaigns/pikp1_pby2/runs/v1_4e"       # 4  OK, 116 negsteer designs
    "experiments/campaigns/pikp1_pwt3/runs/v1_4e"       # 5  OK,  26 negsteer designs
    "experiments/campaigns/pikp1_pwt7/runs/v1_4a"       # 6  OK,  30 negsteer designs
    "experiments/campaigns/pikp1_pwt7/runs/v1_4e"       # 7  ERR, 0 designs passed; clean for work/
)

# Excluded — never executed (no .nextflow/history); preserved for future runs:
#   experiments/campaigns/pikp1_avrpikf/runs/v2_4b
#   experiments/campaigns/pikp1_pby2/runs/v1_4b
#   experiments/campaigns/pikp1_pwt3/runs/v1_4b
#   experiments/campaigns/pikp1_pwt7/runs/v1_4b
# When new runs complete, append them to RUNS and bump --array=0-N.

cd "$PIPELINE_ROOT"

IDX="${SLURM_ARRAY_TASK_ID:-}"
if [[ -z "$IDX" ]]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID not set — this script must be submitted as an array job (sbatch …)" >&2
    exit 1
fi
if (( IDX < 0 || IDX >= ${#RUNS[@]} )); then
    echo "ERROR: SLURM_ARRAY_TASK_ID=$IDX out of range [0, ${#RUNS[@]})" >&2
    exit 1
fi

RUN_DIR="${RUNS[$IDX]}"

echo "============================================================"
echo "cleanup task $IDX / $((${#RUNS[@]} - 1))"
echo "run dir:      $RUN_DIR"
echo "pipeline:     $PIPELINE_ROOT"
echo "node:         $(hostname)"
echo "date:         $(date)"
echo "slurm job:    ${SLURM_ARRAY_JOB_ID:-?}_${SLURM_ARRAY_TASK_ID:-?}"
echo "============================================================"

if [[ ! -d "$RUN_DIR" ]]; then
    echo "ERROR: $RUN_DIR not found relative to $PIPELINE_ROOT" >&2
    exit 1
fi

"$CLEAN" "$RUN_DIR"

echo "============================================================"
echo "done: $(date)"
echo "============================================================"
