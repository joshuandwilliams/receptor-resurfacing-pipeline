#!/usr/bin/env python3
"""Two alignment-free views of what ProteinMPNN produced in design region 1.

Both figures exist because the region cannot be aligned across designs. It is a
hairpin whose loop length RFdiffusion varies, so a residue at a given index
sits at a different place in the fold depending on how long that design's loop
is. Anything positional is therefore unsafe across the campaign, which rules
out per-position conservation profiles, sequence logos and alignment views.

Composition and clustering both sidestep that. Neither needs residues to
correspond between designs.

    diversity_curve         clusters against identity threshold
    composition_grouped     per-residue frequency, amino acids in class order
    composition_panels      the same split into one panel per class

Frequencies are per SEQUENCE, so n = 128 rather than the 64 backbones, because
composition is a property of the sequence ProteinMPNN wrote rather than of the
backbone it wrote it on.

Usage:
    make_figures.py [--outdir figure-options]
"""

import argparse
import collections
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")

INK, MUTED, GRID = "#222222", "#666666", "#DDDDDD"
NATIVE = "#D55E00"

# Standard side-chain classes. Cysteine is put with the polar group rather than
# given its own, since ProteinMPNN never uses it here.
GROUPS = [
    ("Aliphatic", ["A", "V", "L", "I", "M"], "#0072B2"),
    ("Aromatic", ["F", "W", "Y"], "#56B4E9"),
    ("Polar", ["S", "T", "N", "Q", "C"], "#009E73"),
    ("Positive", ["K", "R", "H"], "#CC79A7"),
    ("Negative", ["D", "E"], "#E69F00"),
    ("Special", ["G", "P"], "#999999"),
]


def apply_style():
    plt.rcParams.update({
        "figure.dpi": 130, "savefig.dpi": 300, "savefig.bbox": "tight",
        "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
        "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
        "xtick.color": MUTED, "ytick.color": MUTED,
        "axes.grid": False, "axes.spines.top": False,
        "axes.spines.right": False, "legend.frameon": False,
    })


def load():
    df = pd.read_csv(os.path.join(
        DATA, "crystal_full_test_contig.scored_metadata.csv"))
    designed = [str(v).split("|")[0] for v in df["designed_residues"]]
    native = str(df["native_residues"].iloc[0]).split("|")[0]

    order = [aa for _, aas, _ in GROUPS for aa in aas]
    freqs = {aa: [] for aa in order}
    for seq in designed:
        c = collections.Counter(seq)
        for aa in order:
            freqs[aa].append(100.0 * c[aa] / len(seq))
    nat = collections.Counter(native)
    native_freq = {aa: 100.0 * nat[aa] / len(native) for aa in order}
    return freqs, native_freq, len(designed)


def save(fig, outdir, name):
    path = os.path.join(outdir, f"{name}.png")
    fig.savefig(path)
    plt.close(fig)
    print("wrote", path)


def diversity_curve(outdir):
    d = pd.read_csv(os.path.join(
        DATA, "crystal_full_test_contig.mpnn_cluster_counts.csv"))
    fig, ax = plt.subplots(figsize=(5.6, 3.8))
    ax.plot(d.threshold * 100, d.n_clusters, marker="o", ms=5, lw=1.8,
            color="#0072B2")
    n = int(d.n_sequences.iloc[0])
    ax.axhline(n, ls="--", lw=1.0, color=MUTED)
    ax.annotate(f"{n} sequences", (100, n), xytext=(-4, 4),
                textcoords="offset points", ha="right", va="bottom",
                fontsize=8.5, color=MUTED)
    for t in (90, 50):
        row = d[np.isclose(d.threshold * 100, t)]
        if len(row):
            v = int(row.n_clusters.iloc[0])
            ax.annotate(str(v), (t, v), xytext=(0, 7),
                        textcoords="offset points", ha="center",
                        fontsize=8.5, color=INK)
    ax.set_xlabel("sequence identity threshold (%)")
    ax.set_ylabel("clusters")
    ax.set_ylim(0, n * 1.12)
    fig.tight_layout()
    save(fig, outdir, "diversity_curve")


def _box(ax, aas, freqs, native_freq, colour):
    bp = ax.boxplot([freqs[aa] for aa in aas], patch_artist=True, widths=0.6,
                    medianprops=dict(color=INK, linewidth=1.4),
                    flierprops=dict(marker="o", markersize=2.5,
                                    markerfacecolor=MUTED,
                                    markeredgecolor="none", alpha=0.5))
    for patch in bp["boxes"]:
        patch.set_facecolor(colour)
        patch.set_alpha(0.65)
        patch.set_edgecolor(INK)
        patch.set_linewidth(0.6)
    for w in bp["whiskers"] + bp["caps"]:
        w.set_color(MUTED)
        w.set_linewidth(0.8)
    ax.scatter(range(1, len(aas) + 1), [native_freq[aa] for aa in aas],
               marker="D", s=26, color=NATIVE, zorder=4,
               edgecolor="white", linewidth=0.5)
    ax.set_xticks(range(1, len(aas) + 1))
    ax.set_xticklabels(aas)


def composition_grouped(freqs, native_freq, n, outdir):
    """One axis, amino acids ordered by class with the classes marked."""
    order = [aa for _, aas, _ in GROUPS for aa in aas]
    fig, ax = plt.subplots(figsize=(9.0, 4.2))
    _box(ax, order, freqs, native_freq, "#BBBBBB")
    # Recolour each box by its class and rule between the classes.
    boxes = [p for p in ax.patches]
    i = 0
    edge = 0.5
    for name, aas, colour in GROUPS:
        for _ in aas:
            if i < len(boxes):
                boxes[i].set_facecolor(colour)
            i += 1
        edge += len(aas)
        if edge < len(order) + 0.5:
            ax.axvline(edge, color=GRID, lw=1.0)
        ax.annotate(name, (edge - len(aas) / 2, 1.005),
                    xycoords=("data", "axes fraction"), ha="center",
                    va="bottom", fontsize=8.5, color=MUTED)
    ax.set_ylabel("frequency in design region 1 (%)")
    ax.set_xlabel("amino acid")
    ax.legend(handles=[Line2D([], [], marker="D", linestyle="none",
                              color=NATIVE, markersize=6,
                              label="wild-type frequency")],
              fontsize=8.5, loc="upper right")
    fig.tight_layout()
    save(fig, outdir, "composition_grouped")


def composition_panels(freqs, native_freq, n, outdir):
    """One panel per class, which keeps the x axis short in each."""
    widths = [len(aas) for _, aas, _ in GROUPS]
    fig, axes = plt.subplots(1, len(GROUPS), figsize=(9.6, 3.6), sharey=True,
                             gridspec_kw={"width_ratios": widths})
    for ax, (name, aas, colour) in zip(axes, GROUPS):
        _box(ax, aas, freqs, native_freq, colour)
        ax.set_title(name, fontsize=9.5)
    axes[0].set_ylabel("frequency in design region 1 (%)")
    fig.legend(handles=[Line2D([], [], marker="D", linestyle="none",
                               color=NATIVE, markersize=6,
                               label="wild-type frequency")],
               fontsize=8.5, loc="lower center", bbox_to_anchor=(0.5, 0.0),
               ncol=1)
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    save(fig, outdir, "composition_panels")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=os.path.join(HERE, "figure-options"))
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    apply_style()
    freqs, native_freq, n = load()
    print(f"{n} sequences")
    diversity_curve(args.outdir)
    composition_grouped(freqs, native_freq, n, args.outdir)
    composition_panels(freqs, native_freq, n, args.outdir)


if __name__ == "__main__":
    main()
