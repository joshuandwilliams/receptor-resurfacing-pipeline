#!/usr/bin/env python3
"""Replacements for the survivor-metrics and negative-steering-controls figures.

Both originals had the same kind of problem: they showed a quantity that the
plotting choice had already determined.

survivor_metrics plotted the representative design's confidence against how
many of its three seeds were correctly posed. The representative is the median
of those seeds, so the pass count mechanically decided whether the median seed
passed and confidence followed from that. The rise was the selection rule, not
a finding.

The fix is not to abandon the representative, only the seed-count axis. Each
design already has a single representative prediction, and once the split is a
plain correct against incorrect there is no circularity left: the
representative is chosen on pose and the quantity plotted is confidence, which
took no part in choosing it. So both columns here are one point per design,
split the same way, on the same threshold, which also lets the confidence and
the orthogonal metrics be compared like for like.

negsteer_controls plotted interface Jaccard against pose error for the 83
actually-steered sequences. Jaccard is a second pose measure sitting beside a
pose measure and the paragraph never used it, and excluding the designs whose
initial prediction already passed made steering look like the only route to a
correct pose when it accounts for 20 of the 41. Here all 128 designs are shown,
split by how they got there.

    seed_confidence     confidence per seed by that seed's own pose outcome
    outcome_histogram   all 128 designs by pose error, stacked by route

Usage:
    make_figures.py [--outdir figure-options]
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")

INK, MUTED, GRID = "#222222", "#666666", "#DDDDDD"
CORRECT, WRONG = "#0072B2", "#BBBBBB"   # incorrect is grey, not a second hue
THRESHOLD_LINE = "#D55E00"
CTRL_POLYA, CTRL_SCRAMBLED = "#D62728", "#8E44AD"
INITIAL, STEERED, FAILED = "#7FB3D5", "#0072B2", "#CCCCCC"
RA_EFF_THRESHOLD = 5.0
RA_TICKS = [0.5, 1, 2, 5, 10, 20, 40]

CONF = [("complex_plddt", "Complex pLDDT"), ("iptm", "ipTM"),
        ("ipsae_min", "ipSAE (min)")]   # prefixed rep_ / suffixed _median
# ddG is the quantity, REU the unit it is reported in: Rosetta Energy Units,
# the arbitrary scale of the Rosetta score function, not kcal/mol.
ORTHO = [("sc", "Shape complementarity"), ("bsa", "Buried surface area (Å²)"),
         ("rosetta_ddg", "Rosetta binding ΔΔG (REU)")]


def apply_style():
    plt.rcParams.update({
        "figure.dpi": 130, "savefig.dpi": 300, "savefig.bbox": "tight",
        "font.size": 9, "axes.titlesize": 9.5, "axes.labelsize": 9,
        "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
        "xtick.color": MUTED, "ytick.color": MUTED, "axes.grid": False,
        "axes.spines.top": False, "axes.spines.right": False,
        "legend.frameon": False})


def readable_log(axis, ticks):
    axis.set_ticks(ticks)
    axis.set_ticklabels([f"{v:g}" for v in ticks])
    axis.set_minor_locator(plt.NullLocator())


def save(fig, outdir, name):
    path = os.path.join(outdir, f"{name}.png")
    fig.savefig(path)
    plt.close(fig)
    print("wrote", path)


def _box(ax, groups, colours, labels):
    bp = ax.boxplot(groups, patch_artist=True, widths=0.55,
                    medianprops=dict(color=INK, linewidth=1.4),
                    flierprops=dict(marker="o", markersize=2.5,
                                    markerfacecolor=MUTED,
                                    markeredgecolor="none", alpha=0.45))
    for patch, c in zip(bp["boxes"], colours):
        patch.set_facecolor(c)
        patch.set_alpha(0.65)
        patch.set_edgecolor(INK)
        patch.set_linewidth(0.6)
    for w in bp["whiskers"] + bp["caps"]:
        w.set_color(MUTED)
        w.set_linewidth(0.8)
    ax.set_xticks(range(1, len(labels) + 1))
    ax.set_xticklabels(labels, fontsize=8.5)


def survivor_metric_histograms(outdir):
    surv = pd.read_csv(os.path.join(
        DATA, "crystal_full_test_contig.survivors.csv"))
    surv = surv[surv.row_type == "steered"].copy()
    # One representative per design, split on the same 5 A pose threshold for
    # both columns. "Validated" was a different and stricter thing, tier A/B/C,
    # which also requires no misfolding and no residual contamination, and
    # using it on one column while the other used pose made the two panels
    # answer slightly different questions.
    surv["correct"] = surv.rep_ra_eff_vs_truth_median <= RA_EFF_THRESHOLD

    # Three rows, two columns: confidence down the left, the orthogonal
    # metrics down the right. Taller than it is wide, which suits a page.
    fig, axes = plt.subplots(3, 2, figsize=(6.6, 8.2))
    for row, ((c_col, c_title), (o_col, o_title)) in enumerate(zip(CONF, ORTHO)):
        for ax, col, title in [(axes[row][0], f"rep_{c_col}_median", c_title),
                               (axes[row][1], o_col, o_title)]:
            g = [surv[surv.correct][col].dropna(),
                 surv[~surv.correct][col].dropna()]
            _box(ax, g, [CORRECT, WRONG],
                 [f"Correct\nn={len(g[0])}", f"Incorrect\nn={len(g[1])}"])
            ax.set_title(title)
    axes[0][0].annotate("Boltz-2 confidence", (0.5, 1.22),
                        xycoords="axes fraction", ha="center", fontsize=10)
    axes[0][1].annotate("Orthogonal metrics", (0.5, 1.22),
                        xycoords="axes fraction", ha="center", fontsize=10)
    fig.tight_layout()
    save(fig, outdir, "survivor_metric_histograms")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=os.path.join(HERE, "figure-options"))
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    apply_style()
    survivor_metric_histograms(args.outdir)


if __name__ == "__main__":
    main()
