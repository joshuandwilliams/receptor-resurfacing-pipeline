#!/usr/bin/env python3
"""
merge_orthogonal_metrics.py
---------------------------
P0-31 · Merge the three per-survivor orthogonal-metrics streams
(AF3-no-MSA, biophysical, Rosetta) back onto the P0-29 extended CSV,
and apply the filter cascade.

Inputs (all CSV):
  --input-csv                       cross_sequence_summary_with_interface_metrics.csv
  --af3-summaries-glob              glob for per-survivor AF3 summary CSVs
  --biophysical-summaries-glob      glob for per-survivor biophysical CSVs
  --rosetta-summaries-glob          glob for per-survivor Rosetta CSVs

Filter thresholds:
  --filter-af3-ra-max     AF3-no-MSA best ra_eff (< threshold passes).
                          INFORMATIONAL — flag emitted but does NOT gate
                          passes_orthogonal_filters (demoted 2026-04-28).
  --filter-sc-min         Rosetta Sc (>= passes).  Flag only, not a hard drop.
  --filter-bsa-min        FreeSASA BSA (Å², >= passes).  Flag only.
  --filter-plddt-min      Interface pLDDT (>= passes).  Flag only.

Output: survivors_with_orthogonal_metrics.csv with columns:
  <all original columns from P0-29>,
  af3_nomsa_best_ra_eff, af3_nomsa_mean_ra_eff,
  af3_nomsa_best_iptm, af3_nomsa_mean_iptm,
  af3_nomsa_n_correct_interface, af3_nomsa_total_predictions,
  af3_nomsa_failures,
  bsa, interface_plddt, interface_hbonds, biophysical_failures,
  sc, rosetta_ddg, rosetta_failures,
  orthogonal_flags, passes_orthogonal_filters

passes_orthogonal_filters = 1 iff every NON-AF3 flag is absent.
AF3 flags (af3_nomsa_missing, af3_nomsa_ra_eff_too_high) are
informational only.
"""

from __future__ import annotations

import argparse
import csv
import glob
import sys
from pathlib import Path
from typing import Dict, List, Optional


AF3_COLS = [
    "af3_nomsa_best_ra_eff", "af3_nomsa_mean_ra_eff",
    "af3_nomsa_best_iptm", "af3_nomsa_mean_iptm",
    "af3_nomsa_n_correct_interface", "af3_nomsa_total_predictions",
    "af3_nomsa_failures",
]
BIO_COLS = ["bsa", "interface_hbonds", "biophysical_failures"]
ROSETTA_COLS = ["sc", "rosetta_ddg", "rosetta_failures"]
SUMMARY_COLS = ["orthogonal_flags", "passes_orthogonal_filters"]
ALL_NEW = AF3_COLS + BIO_COLS + ROSETTA_COLS + SUMMARY_COLS


def _load_summaries(pattern: str, key_col: str = "seq_name") -> Dict[str, Dict[str, str]]:
    """Load every CSV matching pattern and index by key_col."""
    index: Dict[str, Dict[str, str]] = {}
    for path in sorted(glob.glob(pattern)):
        with open(path, newline="") as f:
            for row in csv.DictReader(f):
                key = row.get(key_col, "")
                if key:
                    index[key] = row
    return index


def _as_float(value: str) -> Optional[float]:
    if value is None or value == "":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def _apply_filters(
    row: Dict[str, str],
    af3_ra_max: float,
    sc_min: float,
    bsa_min: float,
    plddt_min: float,
) -> None:
    """Populate orthogonal_flags and passes_orthogonal_filters in-place.

    AF3 is INFORMATIONAL — its flag is emitted (so plots can report
    cross-model disagreement) but it does NOT gate
    passes_orthogonal_filters.  Negative steering optimises against
    Boltz; AF3 is a useful sanity check but disagreement is not by
    itself a kill criterion.
    """
    flags: List[str] = []

    # AF3 — INFORMATIONAL flag only (was previously gating; demoted
    # 2026-04-28 so AF3 disagreement is a warning rather than a hard
    # drop).  The flag is still emitted so the plot's filter cascade
    # can show how many survivors AF3 disputes.
    af3_ra = _as_float(row.get("af3_nomsa_best_ra_eff", ""))
    if af3_ra is None:
        flags.append("af3_nomsa_missing")
    elif af3_ra >= af3_ra_max:
        flags.append(f"af3_nomsa_ra_eff_too_high:{af3_ra:.2f}")

    # Sc — flag only.
    sc_val = _as_float(row.get("sc", ""))
    if sc_val is None:
        flags.append("sc_missing")
    elif sc_val < sc_min:
        flags.append(f"sc_too_low:{sc_val:.3f}")

    # BSA — flag only.
    bsa_val = _as_float(row.get("bsa", ""))
    if bsa_val is None:
        flags.append("bsa_missing")
    elif bsa_val < bsa_min:
        flags.append(f"bsa_too_low:{bsa_val:.0f}")

    # Interface pLDDT — flag only.  Sourced from
    # rep_interface_plddt_median, which comes from
    # cross_sequence_summary.csv (computed per-prediction by
    # compute_metrics.py and propagated through the aggregator).
    # Falls back to the legacy interface_plddt column if the row
    # came from an older biophysical_summary.csv (defensive — should
    # not happen post-migration).
    plddt_val = _as_float(
        row.get("rep_interface_plddt_median", "")
        or row.get("interface_plddt", "")
    )
    if plddt_val is None:
        flags.append("interface_plddt_missing")
    elif plddt_val < plddt_min:
        flags.append(f"interface_plddt_too_low:{plddt_val:.3f}")

    row["orthogonal_flags"] = ",".join(flags)
    # Pass requires: every NON-AF3 flag absent.  AF3 flags
    # (af3_nomsa_missing, af3_nomsa_ra_eff_too_high) are informational
    # only and do not gate passing.  Missing non-AF3 metrics still
    # count as flags — strict "every column populated" interpretation.
    gating_flags = [f for f in flags if not f.startswith("af3_nomsa")]
    passes = 1 if not gating_flags else 0
    row["passes_orthogonal_filters"] = str(passes)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-csv", required=True, type=Path)
    ap.add_argument("--af3-summaries-glob", required=True, type=str)
    ap.add_argument("--biophysical-summaries-glob", required=True, type=str)
    ap.add_argument("--rosetta-summaries-glob", required=True, type=str)
    ap.add_argument("--filter-af3-ra-max", type=float, required=True)
    ap.add_argument("--filter-sc-min", type=float, required=True)
    ap.add_argument("--filter-bsa-min", type=float, required=True)
    ap.add_argument("--filter-plddt-min", type=float, required=True)
    ap.add_argument("--output-csv", required=True, type=Path)
    args = ap.parse_args()

    af3 = _load_summaries(args.af3_summaries_glob)
    bio = _load_summaries(args.biophysical_summaries_glob)
    ros = _load_summaries(args.rosetta_summaries_glob)
    print(f"Loaded AF3={len(af3)}  biophysical={len(bio)}  rosetta={len(ros)} summaries")

    with open(args.input_csv, newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            print("ERROR: input CSV has no header.", file=sys.stderr)
            return 1
        in_fieldnames = list(reader.fieldnames)
        rows = list(reader)

    out_fieldnames = list(in_fieldnames)
    for col in ALL_NEW:
        if col not in out_fieldnames:
            out_fieldnames.append(col)

    # Survivors are keyed by mpnn_sequence (the same workdir name used
    # everywhere else in the pipeline).
    n_passing = 0
    for row in rows:
        key = row.get("mpnn_sequence", "")
        af3_row = af3.get(key, {})
        bio_row = bio.get(key, {})
        ros_row = ros.get(key, {})
        for c in AF3_COLS:
            row[c] = af3_row.get(c, "")
        for c in BIO_COLS:
            row[c] = bio_row.get(c, "")
        for c in ROSETTA_COLS:
            row[c] = ros_row.get(c, "")
        _apply_filters(
            row,
            af3_ra_max=args.filter_af3_ra_max,
            sc_min=args.filter_sc_min,
            bsa_min=args.filter_bsa_min,
            plddt_min=args.filter_plddt_min,
        )
        if row.get("passes_orthogonal_filters") == "1":
            n_passing += 1

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=out_fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    print(f"Wrote {args.output_csv}")
    print(f"  {len(rows)} total rows, {n_passing} pass all orthogonal filters.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
