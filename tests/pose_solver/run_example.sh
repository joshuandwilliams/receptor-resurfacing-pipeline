#!/bin/bash
# Constraint-driven pose solver example — historical / hand-debugging.
#
# This script invokes bin/pose_solver.py directly outside the Nextflow
# pipeline so you can iterate on pair choices in ChimeraX without
# spinning up a SLURM job.  The Nextflow path lives in:
#   modules/pose_solver.nf::POSE_SOLVE
#   tests/pose_solver/test_pose_solver.nf
#
# Pair constraints (from ChimeraX visual inspection):
#   A73-B31  }  orientate the two proteins and line up the anti-parallel
#   A71-B33  }  beta strands correctly
#   A8-B22   — pulls the receptor rotation down to build shape complementarity
#
# Defaults match the Nextflow process defaults in main.nf (max distance
# 6.0 Å per pipeline_notes16 §4.1 — natural CA-CA range).
#
# Usage:
#   bash run_example.sh                   # 1000 restarts (default)
#   bash run_example.sh --n-restarts 3000 # more restarts for tighter convergence

set -e
cd "$(dirname "$0")"

REPO_ROOT="$(cd ../.. && pwd)"

python3 "${REPO_ROOT}/bin/pose_solver.py" \
    --binder               ./data/Pikp-1_HMA.pdb \
    --target               ./data/avr-pia.pdb    \
    --binder-chain         A \
    --target-chain         B \
    --pairs                "A73-B31 A71-B33 A8-B22" \
    --exclusions           "A8-B33@4.5" \
    --min-pair-distance    3.5 \
    --max-pair-distance    6.0 \
    --contig-design-region "33-49,69-78" \
    --clash-cutoff         2.0 \
    --contact-cutoff       8.0 \
    --n-restarts           1000 \
    --global-interp \
    --output-prefix        solved_pose \
    "$@"

echo ""
echo "Outputs:"
echo "  solved_pose_posed.pdb    — open in ChimeraX: open solved_pose_posed.pdb"
echo "  solved_pose_heatmap.png  — contact / clash / design region plot"
echo "  solved_pose_results.json — pair distances, clashes, three contig strings"
