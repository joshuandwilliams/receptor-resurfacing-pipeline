#!/usr/bin/env python3
"""Distribution of post-steering receptor-aligned effector RMSD.

Every one of the 128 designed sequences, stacked by the route it took to its
outcome, with the polyA and scrambled sequence controls stacked alongside.

Split out of the campaign metric figures so that it sits with the rest of the
negative-steering results rather than with the confidence comparison.

Usage:
    make_ra_eff_distribution.py
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pipeline_truth import never_steered

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



def ra_eff_distribution_by_outcome(outdir):
    seeds = pd.read_csv(os.path.join(DATA, "campaign_per_seed.csv"))
    cold = (seeds[seeds.kind == "cold_start"]
            .groupby("unit").ra_eff.median().rename("cold"))
    surv = pd.read_csv(os.path.join(
        DATA, "crystal_full_test_contig.survivors.csv"))
    surv = surv[surv.row_type == "steered"].copy()
    d = surv.set_index("mpnn_sequence").join(cold)
    d["final"] = d.rep_ra_eff_vs_truth_median

    never = never_steered()

    def route(r):
        if r.final > RA_EFF_THRESHOLD:
            return "failed"
        return ("correct without steering" if r.name in never
                else "corrected by steering")
    d["route"] = d.apply(route, axis=1)

    # Bins anchored ON the 5 A threshold so no bar straddles it, and spaced by
    # a constant ratio so every bar is the same width on the log axis. Building
    # the two sides separately, as an earlier version did, gave them different
    # widths.
    ratio = 10 ** 0.0875
    bins = np.array([RA_EFF_THRESHOLD * ratio ** k for k in range(-9, 12)])

    ctrl = pd.read_csv(os.path.join(DATA, "campaign_per_seed_controls.csv"))
    ctrl = ctrl[ctrl.kind == "steered"].groupby("unit").ra_eff.median()
    ctrl_colour = {"polyA": CTRL_POLYA, "scrambled": CTRL_SCRAMBLED}

    order = ["correct without steering", "corrected by steering", "failed"]
    colours = {"correct without steering": INITIAL,
               "corrected by steering": STEERED, "failed": FAILED}
    fig, ax = plt.subplots(figsize=(7.4, 4.4))
    series = [d[d.route == k].final for k in order]
    labels = [f"{k[0].upper()}{k[1:]} ({(d.route == k).sum()})" for k in order]
    # The controls ride on the same bars as the designs they sit among, so
    # they read as members of the distribution rather than as annotations.
    for name, v in ctrl.items():
        short = name.replace("input_control_", "")
        series.append(pd.Series([v]))
        labels.append(f"{short} control")
    ax.hist(series, bins=bins, stacked=True,
            color=[colours[k] for k in order]
                  + [ctrl_colour[n.replace("input_control_", "")]
                     for n in ctrl.index],
            edgecolor="black", linewidth=0.5, label=labels)

    ax.axvline(RA_EFF_THRESHOLD, ls="--", lw=1.2, color=THRESHOLD_LINE)
    ax.set_xscale("log")
    readable_log(ax.xaxis, RA_TICKS)
    ax.set_xlabel("Receptor-aligned effector RMSD (Å)")
    ax.set_ylabel("Designs")
    ax.legend(fontsize=8.5, loc="upper left")
    fig.tight_layout()
    save(fig, outdir, "ra_eff_distribution_by_outcome")




def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=os.path.join(HERE, "figure-options"))
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    apply_style()
    ra_eff_distribution_by_outcome(args.outdir)


if __name__ == "__main__":
    main()
