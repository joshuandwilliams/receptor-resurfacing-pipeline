#!/usr/bin/env bash
# clean_finished_run.sh — reclaim disk + inodes from one finished pipeline run.
#
# Usage:
#   scripts/clean_finished_run.sh <run_dir>
#
# Where <run_dir> is e.g. experiments/campaigns/pikp1_avrpia/runs/v1_4b/
#
# Performs, in order:
#   1. Refuse if .nextflow/history is missing (run was never executed).
#   2. Read status from .nextflow/history; clean both OK and ERR runs
#      (warn loudly on ERR), skip anything else.
#   3. Delete <run>/work, <run>/.nextflow/cache, <run>/.nextflow/plr
#      (cache and plr are orphaned once work/ is gone — see
#      notes/inventory/17_storage_bloat_audit.md §7).
#   4. Inside results/negative_steering/runs/, prune unused Boltz diffusion
#      samples (model_1..9 and rank_1..9 + their PAE/PDE/plddt/confidence
#      sidecars) — only model_0/rank_0 is consumed downstream.
#   5. Inside results/negative_steering/runs/, prune Boltz scratch dirs
#      (processed/, lightning_logs/, msa/) that are never read after Boltz
#      returns.
#
# All file operations tolerate missing paths — the script is safe to re-run
# on a partially-cleaned tree.
#
# Run on the HPC login node, NOT over SSHFS — find/unlink over SSHFS is
# orders of magnitude slower than local FS.

set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $(basename "$0") <run_dir>" >&2
    exit 2
fi

RUN="${1%/}"

if [[ ! -d "$RUN" ]]; then
    echo "[$RUN] not a directory; aborting" >&2
    exit 1
fi

HIST="$RUN/.nextflow/history"
if [[ ! -f "$HIST" ]]; then
    echo "[$RUN] SKIP — no .nextflow/history (run was never executed)"
    exit 0
fi

# .nextflow/history is tab-separated, one line per launch:
#   timestamp  duration  run_name  status  session_hash  uuid  command
# Take the status from the last line; if --resume was used, the most recent
# launch's status is what matters.
STATUS=$(awk -F'\t' 'END{print $4}' "$HIST")

case "$STATUS" in
    OK)
        echo "[$RUN] status=OK; cleaning"
        ;;
    ERR)
        echo "[$RUN] status=ERR; cleaning anyway (work/ is the largest waste regardless)"
        ;;
    "")
        echo "[$RUN] SKIP — could not parse .nextflow/history" >&2
        exit 0
        ;;
    *)
        echo "[$RUN] SKIP — unrecognised status '$STATUS'" >&2
        exit 0
        ;;
esac

# ── 1. work/ and orphan Nextflow state ─────────────────────────────────
# These are the largest single contributors. Cache/plr are useful only
# while work/ exists, so they go together.
for p in work .nextflow/cache .nextflow/plr; do
    if [[ -e "$RUN/$p" ]]; then
        echo "[$RUN]   removing $p"
        rm -rf -- "$RUN/$p"
    fi
done

NEGSTEER="$RUN/results/negative_steering/runs"

if [[ -d "$NEGSTEER" ]]; then
    # ── 2. Prune unused Boltz diffusion samples + their sidecars ───────
    # Only model_0 / rank_0 (the top-confidence sample chosen by Boltz's
    # internal ranking) is read by compute_metrics.py and the downstream
    # negsteer post-processing. Samples 1-9 are dead weight.
    #
    # The glob *_model_[1-9]* matches:
    #   input_model_1.pdb           (and 2..4 in current params)
    #   pae_input_model_1.npz       (and 2..4)
    #   pde_input_model_1.npz       (and 2..4)
    #   plddt_input_model_1.npz     (and 2..4)
    #   confidence_input_model_1.json (and 2..4)
    # *_rank_[1-9]* covers older Boltz output naming defensively.
    echo "[$RUN]   pruning unused Boltz samples (*_model_[1-9]*, *_rank_[1-9]*)"
    find "$NEGSTEER" \
        -path '*/boltz_results_input/predictions/*' \
        \( -name '*_model_[1-9]*' -o -name '*_rank_[1-9]*' \) \
        -type f -delete

    # ── 3. Prune Boltz scratch dirs ────────────────────────────────────
    # processed/, lightning_logs/, msa/ inside boltz_results_input/ are
    # Boltz's internal staging — never read by negsteer post-processing.
    # -prune stops find descending into a dir we're about to rm.
    echo "[$RUN]   pruning Boltz scratch dirs (processed, lightning_logs, msa)"
    find "$NEGSTEER" \
        -type d \
        \( -name processed -o -name lightning_logs -o -name msa \) \
        -path '*/boltz_results_input/*' \
        -prune -exec rm -rf {} +
else
    echo "[$RUN]   no results/negative_steering/runs/ — skipping in-results pruning"
fi

# ── 4. Prune RFDiffusion trajectory files ──────────────────────────────────
# RFDiffusion writes a traj/ directory alongside its design PDBs.
# This directory is published to results/rfdiffusion/traj/ but the
# trajectories channel is never consumed downstream (grep main.nf confirms).
# The design PDBs in results/rfdiffusion/ are kept; only traj/ is removed.
RFDIFF_TRAJ="$RUN/results/rfdiffusion/traj"
if [[ -d "$RFDIFF_TRAJ" ]]; then
    echo "[$RUN]   removing results/rfdiffusion/traj/"
    rm -rf -- "$RFDIFF_TRAJ"
else
    echo "[$RUN]   no results/rfdiffusion/traj/ — skipping"
fi

echo "[$RUN] done"
