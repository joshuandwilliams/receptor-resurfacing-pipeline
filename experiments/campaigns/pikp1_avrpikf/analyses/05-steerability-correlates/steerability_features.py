#!/usr/bin/env python3
"""Do any pre-Boltz-2 features separate the designs the pipeline kept?

This is the analysis behind @tbl-steerability-features in Chapter 5. An earlier
version of this folder answered the same question with a LightGBM model read
through SHAP, which is 05-steerability-correlates.qmd. That was dropped: a
gradient boosted model needs a page of Methods to explain, and with 8 features
and 128 designs there is nothing it can find that a per-feature test cannot.

WHAT THE TWO GROUPS ARE. Only the designs negative steering was actually run
on, split by whether it corrected them: 18 corrected against 88 not. The 22
designs the pipeline never steered are excluded, because "steering corrected
it" is undefined for them.

An earlier version split all 128 designs on helpers.load_run's `hit` column,
`cross_tier` in {A, B, C}, giving 38 against 90. That is whether a design left
the whole pipeline in a passing tier, not whether steering worked, and it
included the 22 unsteered designs. It answered a different question.

Per feature: a two-sided Mann-Whitney U test, with the AUC as the effect size.
Mann-Whitney rather than a t-test because these features are not normal and
several are bounded, and its U statistic is a monotone function of the AUC, so
the test and the effect size are the same quantity reported two ways.

Six features were considered and dropped before testing:

  dG_separated, dG_dSASA_density, packstat, delta_unsatHbonds
      Rosetta terms computed on the RFdiffusion backbone before any sequence
      exists, so they describe the backbone, not the designed sequence.
  nres_int
      Spearman 0.887 with design region 1 length, so redundant with it.
  mpnn_dr_ratio
      A ratio of the two ProteinMPNN scores already tested separately.

loop_has_proline was also dropped: it is binary, so a rank test on it is a
two-by-two contingency and the AUC is not interpretable as separability.

Benjamini-Hochberg q values are reported but no correction is needed for the
conclusion: no raw p reaches 0.05, and a correction can only raise it.

Usage:
    steerability_features.py [--out steerability_features.csv]
"""
import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_auc_score

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import helpers as H  # noqa: E402

# helpers reads the campaign from the HPC mount by default. Point it at the
# run tree copied into this repo so the table can be rebuilt offline.
H.HPC_BASE = Path(HERE).parents[1] / "runs"

# Column, and the label as it appears in the thesis table.
FEATURES = [
    ("sc_value", "Shape complementarity"),
    ("loop_net_charge", "DR1 net charge"),
    ("loop_hydrophob_frac", "DR1 hydrophobic fraction"),
    ("loop_aromatic_frac", "DR1 aromatic fraction"),
    ("mpnn_score", "ProteinMPNN score, whole receptor"),
    ("design_region_score", "ProteinMPNN score, both design regions"),
    ("loop1_len", "DR1 length"),
    ("dSASA_int", "Buried interface area (Å²)"),
    # The pose before steering. Asks whether how badly a design starts out
    # predicts whether steering can rescue it.
    ("cold_ra_eff", "Receptor-aligned effector RMSD before steering (Å)"),
]


def benjamini_hochberg(p):
    """BH q values, step-up from the largest p."""
    p = np.asarray(p, dtype=float)
    order = np.argsort(p)
    n = len(p)
    q = np.empty(n)
    running = 1.0
    for rank in range(n - 1, -1, -1):
        i = order[rank]
        running = min(running, p[i] * n / (rank + 1))
        q[i] = running
    return q


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out",
                    default=os.path.join(HERE, "steerability_features.csv"))
    args = ap.parse_args()

    df = H.load_run("crystal", H.RUNS["crystal"])
    # Steering outcome and pre-steering pose, from the pipeline's own record.
    chk = pd.read_csv(os.path.join(HERE, "..", "04-negative-steering-outcomes",
                                   "campaign_interface_check.csv"))
    df = df.merge(chk[["sequence", "was_steered", "steered_correct",
                       "cold_ra_eff"]],
                  left_on="mpnn_sequence", right_on="sequence", how="left")
    df = df[df.was_steered.astype(bool)].copy()
    y = df["steered_correct"].values.astype(bool)
    print(f"steered = {len(df)}   corrected = {y.sum()}   "
          f"not corrected = {(~y).sum()}   features = {len(FEATURES)}\n")

    rows = []
    for col, label in FEATURES:
        v = pd.to_numeric(df[col], errors="coerce")
        ok = v.notna().values
        rows.append({
            "Feature": label,
            "Median success": np.median(v[ok][y[ok]]),
            "Median failure": np.median(v[ok][~y[ok]]),
            "AUC": roc_auc_score(y[ok], v[ok]),
            "p": mannwhitneyu(v[ok][y[ok]], v[ok][~y[ok]]).pvalue,
        })
    out = pd.DataFrame(rows)
    out["q"] = benjamini_hochberg(out.p.values)
    out = out.sort_values("p").reset_index(drop=True)

    print(out.to_string(index=False, float_format=lambda x: f"{x:,.4g}"))
    out.to_csv(args.out, index=False)
    print(f"\nwrote {args.out}")
    print(f"lowest raw p = {out.p.min():.3f}, lowest q = {out.q.min():.3f}")


if __name__ == "__main__":
    main()
