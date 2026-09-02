#!/usr/bin/env python3
"""Diversity of the two design regions, both curves on one axis.

Both curves count clusters against a sequence identity threshold. Design
region 1 varies from 10 to 20 residues, so pairs are aligned before identity
is taken. Design region 2 is a fixed 6 residues, so identity is exactly the
number of matching positions over 6 and no alignment is involved. That leaves
only seven attainable thresholds, so its points sit at sixths of the axis
rather than on region 1's grid.

Usage:
    make_figure.py [--outdir figure-options]
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
INK, MUTED = "#222222", "#666666"
R1, R2 = "#0072B2", "#D55E00"


def style():
    plt.rcParams.update({
        "figure.dpi": 130, "savefig.dpi": 300, "savefig.bbox": "tight",
        "font.size": 9, "axes.labelsize": 9, "axes.edgecolor": MUTED,
        "axes.labelcolor": INK, "text.color": INK, "xtick.color": MUTED,
        "ytick.color": MUTED, "axes.grid": False, "axes.spines.top": False,
        "axes.spines.right": False, "legend.frameon": False})


def label_points(ax, x, y, colour, dy):
    for xi, yi in zip(x, y):
        ax.annotate(str(int(yi)), (xi, yi), xytext=(0, dy),
                    textcoords="offset points", ha="center", va="center",
                    fontsize=8.5, color=colour)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=os.path.join(HERE, "figure-options"))
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    style()

    r1 = pd.read_csv(os.path.join(HERE, "region1_cluster_counts.csv"))
    r2 = pd.read_csv(os.path.join(HERE, "region2_cluster_counts.csv"))
    n = int(r1.n_sequences.iloc[0])

    fig, ax = plt.subplots(figsize=(5.8, 3.9))
    x1, y1 = r1.threshold * 100, r1.n_clusters
    x2, y2 = r2.threshold * 100, r2.n_clusters
    ax.plot(x1, y1, marker="o", ms=5, lw=1.8, color=R1, zorder=3,
            label="Design region 1, 10 to 20 residues")
    ax.plot(x2, y2, marker="o", ms=5, lw=1.8, color=R2, zorder=2,
            label="Design region 2, 6 residues")
    label_points(ax, x1, y1, R1, 11)
    label_points(ax, x2, y2, R2, -11)
    ax.set_xlabel("Sequence identity threshold (%)")
    ax.set_ylabel("Clusters")
    ax.set_xlim(-4, 104)
    ax.set_xticks(range(0, 101, 10))
    ax.set_ylim(-16, n * 1.16)
    ax.legend(fontsize=8.5, loc="upper left")
    fig.tight_layout()
    path = os.path.join(args.outdir, "diversity_curves.png")
    fig.savefig(path)
    plt.close(fig)
    print("wrote", path)


if __name__ == "__main__":
    main()
