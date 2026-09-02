#!/usr/bin/env python3
"""AlphaFold3 pose accuracy on the designed sequences, two views.

af3_ra_eff_distribution   the same stacked-histogram form as the Boltz-2
                          figure in the results chapter, but for AlphaFold3,
                          splitting designs by whether the pipeline steered
                          them.
af3_vs_boltz2_scatter     each design's Boltz-2 representative RMSD against
                          its AlphaFold3 best, on matched log axes.

AlphaFold3 values are the best of a sequence's 15 predictions, so they flatter
AlphaFold3 relative to Boltz-2's median over three seeds. That asymmetry is
deliberate: it is the most generous reading of AlphaFold3 available, and it
still places nothing at the intended interface.

Usage:
    make_af3_figures.py
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pipeline_truth import never_steered

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "figure-options")
INK, MUTED = "#222222", "#666666"
UNSTEERED, STEERED = "#56B4E9", "#0072B2"
THRESHOLD_LINE = "#D55E00"
RA_EFF_THRESHOLD = 5.0
TICKS = [1, 2, 5, 10, 20, 40]


def style():
    plt.rcParams.update({
        "figure.dpi": 130, "savefig.dpi": 300, "savefig.bbox": "tight",
        "font.size": 9, "axes.labelsize": 9, "axes.edgecolor": MUTED,
        "axes.labelcolor": INK, "text.color": INK, "xtick.color": MUTED,
        "ytick.color": MUTED, "axes.grid": False, "axes.spines.top": False,
        "axes.spines.right": False, "legend.frameon": False})


def load():
    af3 = pd.read_csv(os.path.join(HERE, "data",
                                   "survivors_with_orthogonal_metrics.csv"))
    boltz = pd.read_csv(os.path.join(
        HERE, "data", "crystal_full_test_contig.survivors.csv"))
    d = (af3.set_index("mpnn_sequence")[["af3_nomsa_best_ra_eff"]]
         .join(boltz.set_index("mpnn_sequence")[["rep_ra_eff_vs_truth_median"]]))
    d = d[~d.index.str.contains("control")].dropna()
    d["steered"] = ~d.index.isin(never_steered())
    return d


def distribution(d):
    lo = min(d.af3_nomsa_best_ra_eff.min(), 1.0)
    hi = d.af3_nomsa_best_ra_eff.max() * 1.05
    bins = np.geomspace(lo, hi, 22)
    fig, ax = plt.subplots(figsize=(6.6, 3.9))
    groups = [(~d.steered, "Never steered", UNSTEERED),
              (d.steered, "Steered", STEERED)]
    ax.hist([d.af3_nomsa_best_ra_eff[m] for m, _, _ in groups], bins=bins,
            stacked=True, color=[c for _, _, c in groups],
            label=[f"{lab} ({int(m.sum())})" for m, lab, _ in groups],
            edgecolor="black", linewidth=0.5)
    ax.axvline(RA_EFF_THRESHOLD, color=THRESHOLD_LINE, ls="--", lw=1.4)
    ax.set_xscale("log")
    ax.set_xticks(TICKS)
    ax.set_xticklabels([f"{t:g}" for t in TICKS])
    ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set_xlabel("AlphaFold3 receptor-aligned effector RMSD (Å), best of 15")
    ax.set_ylabel("Designs")
    ax.legend(fontsize=8.5, loc="upper left")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "af3_ra_eff_distribution.png"))
    plt.close(fig)


def scatter(d):
    fig, ax = plt.subplots(figsize=(5.2, 5.0))
    for m, lab, c in [(~d.steered, "Never steered", UNSTEERED),
                      (d.steered, "Steered", STEERED)]:
        ax.scatter(d.rep_ra_eff_vs_truth_median[m], d.af3_nomsa_best_ra_eff[m],
                   s=34, color=c, alpha=0.75, edgecolor="white", linewidth=0.6,
                   label=f"{lab} ({int(m.sum())})", zorder=3)
    lim = (1.0, max(d.max().max(), 40) * 1.1)
    ax.plot(lim, lim, ls="--", lw=1.0, color=MUTED, zorder=1)
    for f in (ax.axvline, ax.axhline):
        f(RA_EFF_THRESHOLD, color=THRESHOLD_LINE, ls="--", lw=1.2, zorder=2)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(lim); ax.set_ylim(lim)
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_locator(matplotlib.ticker.FixedLocator(TICKS))
        axis.set_major_formatter(matplotlib.ticker.FixedFormatter(
            [f"{t:g}" for t in TICKS]))
        axis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set_xlabel("Boltz-2 RMSD (Å), representative median")
    ax.set_ylabel("AlphaFold3 RMSD (Å), best of 15")
    ax.legend(fontsize=8.5, loc="lower right")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "af3_vs_boltz2_scatter.png"))
    plt.close(fig)


def main():
    os.makedirs(OUT, exist_ok=True)
    style()
    d = load()
    print(f"{len(d)} designs   steered {int(d.steered.sum())}   "
          f"never steered {int((~d.steered).sum())}")
    print(f"AlphaFold3 below {RA_EFF_THRESHOLD} A: "
          f"{int((d.af3_nomsa_best_ra_eff <= RA_EFF_THRESHOLD).sum())}")
    print(f"Boltz-2 below {RA_EFF_THRESHOLD} A: "
          f"{int((d.rep_ra_eff_vs_truth_median <= RA_EFF_THRESHOLD).sum())}")
    print("Spearman Boltz-2 vs AlphaFold3: "
          f"{d.rep_ra_eff_vs_truth_median.corr(d.af3_nomsa_best_ra_eff, method='spearman'):.3f}")
    distribution(d)
    scatter(d)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
