#!/usr/bin/env python3
"""The constraint level-and-rank figure, campaign sequences only.

Migrated from structure-negative-steering analysis/06-confidence-synthesis,
which built it alongside the benchmark cohorts. The benchmark half is not
reproduced here: that document itself notes that only 3 of the benchmark's 12
pose-matched targets have the two arms on the same sequence, so the thesis
figure was already campaign-only. Everything the benchmark needs stays in the
benchmarking repository.

Reads the two campaign cohorts, unconstrained and constrained, keeps the
pose-matched pairs, and splits them by whether both arms' representative
designs are the same sequence.

Usage:
    rebuild_constraint_figure.py
"""

import os

import numpy as np
import pandas as pd

import confidence_common as cc

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
OUT = os.path.join(HERE, "thesis-figures", "constraint_level_and_rank.png")


def main():
    cc.apply_style()
    unconstrained = cc.build_representative_table(
        pd.read_csv(os.path.join(DATA, "campaign_per_seed.csv")))
    constrained = cc.build_representative_table(
        pd.read_csv(os.path.join(DATA, "campaign_per_seed_constrained.csv")))

    u, c = cc.pose_matched_pairs(unconstrained, constrained)
    match = (pd.read_csv(os.path.join(HERE, "sequence_match.csv"))
             .set_index("unit")["sequence_matched"].astype(bool).reindex(u.index))
    groups = pd.Series(np.where(match, *cc.SEQ_GROUPS), index=u.index)
    print(groups.value_counts().reindex(cc.SEQ_GROUPS).rename("units")
          .to_string())

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    cc.fig_constraint_level_and_rank(u, c, groups, OUT)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
