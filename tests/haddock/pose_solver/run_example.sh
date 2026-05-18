#!/bin/bash
# Constraint-driven pose solver example
#
# Pair constraints (from ChimeraX visual inspection):
#   A73-B48  }  orientate the two proteins and line up the anti-parallel
#   A71-B50  }  beta strands correctly
#   A8-B39   — pulls the receptor rotation down to build shape complementarity
#
# Usage:
#   bash run_example.sh                  # 300 restarts (default)
#   bash run_example.sh --n-restarts 500 # more restarts for tighter convergence

set -e
cd "$(dirname "$0")"

python3 pose_solver.py \
    --binder               ../data/Pikp-1_HMA.pdb \
    --target               ../data/avr-pia.pdb    \
    --binder-chain         A \
    --target-chain         B \
    --pairs                "A73-B31 A71-B33 A8-B22" \
    --exclusions           "A8-B33@4.5" \
    --max-pair-distance    4.0 \
    --contig-design-region "33-49,69-78" \
    --clash-cutoff         2.0 \
    --contact-cutoff       8.0 \
    --output-prefix        solved_pose \
    "$@"

echo ""
echo "Outputs:"
echo "  solved_pose_posed.pdb    — open in ChimeraX: open solved_pose_posed.pdb"
echo "  solved_pose_heatmap.png  — contact / clash / design region plot"
echo "  solved_pose_results.json — pair distances, clashes, three contig strings"
