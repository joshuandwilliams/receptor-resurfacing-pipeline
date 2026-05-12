#!/usr/bin/env python3
"""
parse_af3_output.py
-------------------
P0-31 · Post-process a single AF3-no-MSA prediction directory.

CLI wrapper around bin/af3_confidence.AF3ConfidenceAggregate (the
Phase 4 deep-module type at Tier 0.4).  The wrapper preserves the
existing CLI contract so the AF3_PARSE_OUTPUT Nextflow process keeps
working unchanged.

Output columns (matches AF3ConfidenceAggregate.to_summary_row plus the
threshold-dependent n_correct_interface column):
  seq_name
  af3_nomsa_best_ra_eff         min across all seed×sample predictions
  af3_nomsa_mean_ra_eff         mean across all predictions
  af3_nomsa_best_iptm           max ipTM across predictions
  af3_nomsa_mean_iptm           mean ipTM
  af3_nomsa_n_correct_interface count of predictions where ra_eff < threshold
  af3_nomsa_total_predictions   total predictions parsed
  af3_nomsa_failures            comma-separated failure tags

AF3 output layout (NBI build, Dec 2024):
  output/<name>/
      <name>_model.cif                  top-ranked structure
      seed-<N>_sample-<M>/model.cif     all seed × sample predictions
      summary_confidences.json          per-prediction ipTM/pTM/ranking_score
"""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List

_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

from af3_confidence import AF3ConfidenceAggregate  # noqa: E402


@dataclass
class _Thresholds:
    """Tiny shim so AF3ConfidenceAggregate.n_correct_interface() can be
    called with an object exposing the threshold attribute it expects,
    without pulling in the full PipelineInternalThresholds (which would
    add an unneeded dependency for the standalone CLI)."""
    orthogonal_af3_ra_max: float


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seq-name", required=True)
    ap.add_argument("--af3-output-dir", required=True, type=Path)
    ap.add_argument("--ground-truth", required=True, type=Path)
    ap.add_argument("--receptor-chain", default="A")
    ap.add_argument("--effector-chain", default="B")
    ap.add_argument("--ra-eff-threshold", type=float, default=5.0)
    ap.add_argument("--output-csv", required=True, type=Path)
    args = ap.parse_args()

    aggregate = AF3ConfidenceAggregate.from_output_dir(
        af3_output_dir=args.af3_output_dir,
        ground_truth=args.ground_truth,
        receptor_chain=args.receptor_chain,
        effector_chain=args.effector_chain,
    )

    row = aggregate.to_summary_row(args.seq_name)
    thresholds = _Thresholds(orthogonal_af3_ra_max=args.ra_eff_threshold)
    row["af3_nomsa_n_correct_interface"] = str(
        aggregate.n_correct_interface(thresholds)
    )

    fieldnames: List[str] = [
        "seq_name",
        "af3_nomsa_best_ra_eff",
        "af3_nomsa_mean_ra_eff",
        "af3_nomsa_best_iptm",
        "af3_nomsa_mean_iptm",
        "af3_nomsa_n_correct_interface",
        "af3_nomsa_total_predictions",
        "af3_nomsa_failures",
    ]
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        w.writerow(row)

    print(
        f"[{args.seq_name}] best ra_eff = {row['af3_nomsa_best_ra_eff']} Å, "
        f"n_correct = {row['af3_nomsa_n_correct_interface']}/"
        f"{row['af3_nomsa_total_predictions']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
