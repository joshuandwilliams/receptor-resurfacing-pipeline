#!/usr/bin/env python3
"""What ProteinMPNN wrote, plotted against depth into the hairpin.

Three views of the same table, differing in how much they commit to a binning.

    composition_stack       the whole composition, in equal-occupancy depth bins
    composition_stack_fixed the same, in fixed-width bins instead
    composition_stack_wt    fixed bins, with the wild-type residues marked
    composition_stack_wt_bars   the same, wild type as a paired narrow bar
    class_trends            one line per side-chain class, sliding window
    depth_by_amino_acid     the depth distribution of each amino acid, no bins

The bins are equal-occupancy rather than equal-width because depth is not
uniformly sampled. Every design contributes residues near the anchors and only
the long ones reach the far end, so fixed-width bins hold anywhere from 8 to 266
residues and the sparse ones swing wildly. Quantile bins put the same number of
residues in each, which makes the error bars comparable across the axis.

Five of them, because the bin has to be wider than one strand period. A beta
strand advances about 3.3 A per residue and alternate residues face the buried
core, so composition oscillates with a period near 6.6 A along the depth axis.
Eight bins are 2 to 3 A wide, which samples that oscillation rather than
averaging over it, and the aliphatic frequency comes out as 92, 49, 83, 54, 13
per cent: a sawtooth that is the strand register showing through, not a trend in
depth. Five bins are 4 to 10 A wide and give 74, 69, 35, 24, 16 per cent. The
choice is set by the periodicity, not by which picture is tidier.

The fixed-width version is the same picture drawn the other way, with bins of
equal width and unequal counts, so the depth axis is linear and the bars are
directly comparable in position. It costs the balanced counts: at 5 A the bins
hold 256, 380, 388, 490, 252, 126 and 8 residues, so the last bar rests on eight
residues from the few longest loops and should not be read as a level. The n on
each bar is there to make that visible rather than to be looked up.

class_trends drops bins for a window that slides along the depth axis in small
steps, which gives a curve rather than five joined points and removes the
dependence on where a particular set of edges happened to fall. The window is
6 A wide, chosen the same way as before: a beta strand advances about 3.3 A per
residue and alternate residues face the buried core, so a window narrower than
the resulting 6.6 A period tracks the strand register instead of averaging over
it. Adjacent points share most of their residues, so the curve is smooth by
construction and the intervals on neighbouring points are not independent. They
say how well each point is pinned down, not how much the curve could move as a
whole.

Two versions carry the wild type. It has 13 residues in the region against 1900,
and they fall 3, 5, 4, 1 and 0 across the five bins, so a wild-type composition
per bin is built on single figures: the 20 to 26 A bin would read 100 per cent
aliphatic off one leucine. composition_stack_wt therefore marks the 13 residues
individually at the depths they actually occupy, coloured by class and labelled,
which is a statement about where they are rather than a frequency. The
paired-bar version draws the direct comparison instead, a narrow faded bar
beside each design bar, and carries its n above each one because those
percentages move in steps of a third or a fifth. The 20 to 26 A bar is one
residue and the 26 to 32 A bin has none, so those two are labels rather than
measurements.

The last figure commits to no binning at all, which is the point of including
it: if the trend only exists for one choice of bin edges it is not a trend.

Strand and coil come from DSSP on the same backbones. The band marks where the
two assignments cross over, so the compositional change can be read against the
secondary structure rather than asserted alongside it.

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
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

HERE = os.path.dirname(os.path.abspath(__file__))

INK, MUTED, GRID = "#222222", "#666666", "#DDDDDD"
NATIVE = "#D55E00"
BAND = "#EFEFEF"
N_BINS = 5

# The same classes and colours as the composition figure, so the two can be read
# together. Cysteine sits with the polar group; ProteinMPNN never uses it here.
GROUPS = [
    ("Aliphatic", ["A", "V", "L", "I", "M"], "#0072B2"),
    ("Aromatic", ["F", "W", "Y"], "#56B4E9"),
    ("Polar", ["S", "T", "N", "Q", "C"], "#009E73"),
    ("Positive", ["K", "R", "H"], "#CC79A7"),
    ("Negative", ["D", "E"], "#E69F00"),
    ("Special", ["G", "P"], "#999999"),
]
CLASS_OF = {aa: name for name, aas, _ in GROUPS for aa in aas}
COLOUR_OF = {name: c for name, _, c in GROUPS}


def apply_style():
    plt.rcParams.update({
        "figure.dpi": 130, "savefig.dpi": 300, "savefig.bbox": "tight",
        "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
        "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
        "xtick.color": MUTED, "ytick.color": MUTED,
        "axes.grid": False, "axes.spines.top": False,
        "axes.spines.right": False, "legend.frameon": False,
    })


def save(fig, outdir, name):
    path = os.path.join(outdir, f"{name}.png")
    fig.savefig(path)
    plt.close(fig)
    print("wrote", path)


def load(path):
    df = pd.read_csv(path)
    df["group"] = df.aa.map(CLASS_OF)
    return df[df.design != "wild-type"].copy(), df[df.design == "wild-type"].copy()


def strand_coil_band(designs):
    """Depth range over which DSSP flips from mostly strand to mostly coil."""
    edges = np.arange(0, designs.depth.max() + 1, 1.0)
    mid, frac = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        w = designs[(designs.depth >= lo) & (designs.depth < hi)]
        if len(w) >= 20:
            mid.append((lo + hi) / 2)
            frac.append((w.ss == "coil").mean())
    mid, frac = np.array(mid), np.array(frac)
    inside = mid[(frac > 0.25) & (frac < 0.75)]
    return (inside.min(), inside.max()) if len(inside) else (np.nan, np.nan)


def draw_band(ax, band, where="top"):
    """Shade the depths over which DSSP flips from mostly strand to mostly coil."""
    lo, hi = band
    if np.isnan(lo):
        return
    ax.axvspan(lo, hi, color=BAND, zorder=0)
    y, dy, va = (1.0, -10, "top") if where == "top" else (0.0, 6, "bottom")
    ax.annotate("strand", (lo, y), xytext=(-4, dy), ha="right", va=va,
                xycoords=("data", "axes fraction"),
                textcoords="offset points", fontsize=8, color=MUTED)
    ax.annotate("loop", (hi, y), xytext=(4, dy), ha="left", va=va,
                xycoords=("data", "axes fraction"),
                textcoords="offset points", fontsize=8, color=MUTED)


def wilson(k, n, z=1.96):
    """Wilson interval, which behaves at the 0% and 100% ends where the normal
    approximation does not, and several classes sit at 0% in the first bin."""
    if n == 0:
        return 0.0, 0.0
    p = k / n
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, centre - half), min(1.0, centre + half)


BIN_WIDTH = 6.0
MIN_BIN = 30
WINDOW, STEP, MIN_WINDOW = 6.0, 0.25, 50
# The grid starts at 2 A rather than 0. Depth cannot approach zero here: the
# two anchor CA atoms are 7.6 to 7.8 A apart, so the midpoint between them is
# already 3 to 4 A from any residue next to an anchor, and the shallowest of
# the 1900 is 3.10 A. A grid from 0 spends its first bin mostly on depths that
# cannot occur and leaves a sliver past 30; from 2 it is five bins of 6 A
# covering 2 to 32, which is the occupied range and nothing else.
BIN_ORIGIN = 2.0


def binned(designs, width=None, min_bin=MIN_BIN, origin=BIN_ORIGIN):
    """Equal-occupancy bins by default, fixed-width ones if a width is given."""
    d = designs.copy()
    if width is None:
        d["bin"] = pd.qcut(d.depth, N_BINS, labels=False, duplicates="drop")
        edges = None
    else:
        d["bin"] = ((d.depth - origin) // width).astype(int)
        # Only the longest few loops reach the far end, so the top bin can hold a
        # handful of residues from a handful of designs. Folding it into the one
        # below keeps the bars comparable; the label says where it now ends.
        while d.bin.nunique() > 1 and (d.bin == d.bin.max()).sum() < min_bin:
            d.loc[d.bin == d.bin.max(), "bin"] = d.bin.max() - 1
        top = d.bin.max()
        edges = {b: (origin + b * width, origin + (b + 1) * width)
                 for b in d.bin.unique()}
        edges[top] = (origin + top * width,
                      min(origin + (top + 1) * width,
                          float(np.ceil(d.depth.max()))))
    stats = d.groupby("bin").agg(lo=("depth", "min"), hi=("depth", "max"),
                                 mid=("depth", "median"), n=("depth", "size"))
    if edges is not None:
        # Label by the bin's own edges, not by the residues that landed in it,
        # or a sparse bin would be labelled narrower than it is.
        stats["lo"] = [edges[b][0] for b in stats.index]
        stats["hi"] = [edges[b][1] for b in stats.index]
    return d, stats


def _bar_x(depth, stats, half=0.41):
    """Map a depth onto its bar, so a marker sits where the residue really is."""
    for i, (b, row) in enumerate(stats.iterrows()):
        if row.lo <= depth <= row.hi:
            return i - half + 2 * half * (depth - row.lo) / (row.hi - row.lo)
    return None


def wt_residue_strip(ax, wt, stats):
    """The 13 wild-type residues as points at their own depths, not as a rate.

    Several sit within a fraction of an angstrom of each other, so a marker that
    would land on top of its neighbour is stepped to a second row instead of
    being nudged along the axis, which would misplace it in depth.
    """
    ax.set_ylim(0, 1)
    # Three residues can fall within a tenth of an angstrom, so a marker takes
    # the first row that is clear rather than simply alternating.
    ROWS, GAP = 3, 0.075
    occupied = [None] * ROWS
    placed = []
    for r in wt.sort_values("depth").itertuples():
        x = _bar_x(r.depth, stats)
        if x is None:
            continue
        row = next((i for i in range(ROWS)
                    if occupied[i] is None or x - occupied[i] >= GAP), 0)
        occupied[row] = x
        placed.append((x, row, r.aa))
    ys = [0.80, 0.50, 0.20]
    for x, row, aa in placed:
        y = ys[row]
        ax.scatter([x], [y], s=52, color=COLOUR_OF[CLASS_OF[aa]],
                   edgecolor=INK, linewidth=0.6, zorder=4, clip_on=False)
        ax.annotate(aa, (x, y), ha="center", va="center", fontsize=6.5,
                    color="white", zorder=5, fontweight="bold")
    ax.set_yticks([])
    for side in ("left", "bottom", "right", "top"):
        ax.spines[side].set_visible(False)
    ax.set_ylabel("wild type", fontsize=8, rotation=0, ha="right", va="center")


def wt_paired_bars(ax, wt, stats, xs, offset, width):
    """A narrow faded bar beside each design bar, from the 13 wild-type residues.

    Fading alone is not enough to tell the two bars apart, because a 50 per cent
    aliphatic blue is close to a full-strength aromatic blue. The wild-type bar
    therefore also carries a dark outline, and its count goes under the axis
    beside the design count rather than above the bar, where it would fall in the
    gap under the coil strip and be painted over.
    """
    wt = wt.copy()
    wt["bin"] = [next((i for i, (_, r) in enumerate(stats.iterrows())
                       if r.lo <= v <= r.hi), None) for v in wt.depth]
    counts = []
    for i in xs:
        w = wt[wt.bin == i]
        counts.append(len(w))
        bottom = 0.0
        for name, _, colour in GROUPS:
            if not len(w):
                continue
            v = 100.0 * (w.group == name).sum() / len(w)
            ax.bar([i + offset], [v], bottom=[bottom], width=width,
                   color=colour, alpha=0.5, edgecolor=INK, linewidth=0.7)
            bottom += v
    return counts


def composition_stack(designs, band, outdir, width=None, origin=BIN_ORIGIN,
                      wt=None, wt_style="dots", outname="composition_stack"):
    """Everything at once: the full class composition of each depth bin."""
    d, stats = binned(designs, width, origin=origin)
    counts = (d.groupby(["bin", "group"]).size().unstack(fill_value=0)
              .reindex(columns=[g for g, _, _ in GROUPS], fill_value=0))
    pct = 100.0 * counts.div(counts.sum(axis=1), axis=0)

    if wt is not None and wt_style == "dots":
        fig, (top, ax, foot) = plt.subplots(
            3, 1, figsize=(8.0, 6.0), sharex=True,
            gridspec_kw={"height_ratios": [1, 7, 1.1], "hspace": 0.08})
    else:
        foot = None
        fig, (top, ax) = plt.subplots(
            2, 1, figsize=(8.0, 5.0), sharex=True,
            gridspec_kw={"height_ratios": [1, 7], "hspace": 0.08})
    x = np.arange(len(pct))
    # Paired bars need room beside each design bar; otherwise it takes the lot.
    paired = wt is not None and wt_style == "bars"
    main_x = x - 0.11 if paired else x
    main_w = 0.58 if paired else 0.82
    bottom = np.zeros(len(pct))
    for name, _, colour in GROUPS:
        ax.bar(main_x, pct[name], bottom=bottom, width=main_w, color=colour,
               edgecolor="white", linewidth=0.7, label=name)
        bottom += pct[name].values

    # DSSP for the same residues, so the compositional change can be read
    # against the secondary structure instead of being asserted beside it.
    coil = 100.0 * d.groupby("bin").ss.apply(lambda s: (s == "coil").mean())
    top.bar(main_x, coil, width=main_w, color="#BBBBBB", edgecolor="white",
            linewidth=0.7)
    top.set_ylim(0, 100)
    top.set_yticks([0, 100])
    top.set_yticklabels(["0", "100"], fontsize=7.5)
    top.set_ylabel("% coil", fontsize=8, rotation=0, ha="right", va="center")
    top.spines["left"].set_visible(True)
    for i, v in enumerate(coil):
        top.annotate(f"{v:.0f}", (main_x[i], v), xytext=(0, 2), ha="center",
                     va="bottom", textcoords="offset points", fontsize=7.5,
                     color=MUTED)
    ax.set_xticks(x)
    wt_counts = None
    if wt is not None:
        if wt_style == "dots":
            wt_residue_strip(foot, wt, stats)
        else:
            wt_counts = wt_paired_bars(ax, wt, stats, x, 0.29, 0.22)
    if width is None:
        labels = [f"{lo:.0f}–{hi:.0f}" for lo, hi in zip(stats.lo, stats.hi)]
    else:
        # Counts go under each bar rather than over it, where they would fall in
        # the gap beneath the coil strip and be painted over.
        labels = [f"{lo:.0f}–{hi:.0f}\nn={int(n)}"
                  for lo, hi, n in zip(stats.lo, stats.hi, stats.n)]
        if wt_counts is not None:
            labels = [f"{lab}\nwt n={c}" for lab, c in zip(labels, wt_counts)]
    ax.set_xticklabels(labels, fontsize=8.5)
    ax.set_xlim(-0.6, len(pct) - 0.4)
    ax.set_ylim(0, 100)
    if width is None:
        # Equal occupancy means every bin holds the same count by construction,
        # so the count belongs in the axis label rather than over each bar.
        ax.set_xlabel("Depth from the hairpin anchors (Å), equal-occupancy bins "
                      f"of {int(stats.n.iloc[0])} residues")
    else:
        ax.set_xlabel("Depth from the hairpin anchors (Å), "
                      f"{width:g} Å bins")
    ax.set_ylabel("Composition of design region 1 (%)")
    if foot is None:
        ax.legend(ncol=6, fontsize=8.5, loc="lower center",
                  bbox_to_anchor=(0.5, -0.30))
        fig.tight_layout()
    else:
        # The strip sits between the bars and the axis label, so the legend goes
        # on the figure rather than hanging off the bar axes into it.
        handles, labels = ax.get_legend_handles_labels()
        fig.legend(handles, labels, ncol=6, fontsize=8.5, loc="lower center",
                   bbox_to_anchor=(0.5, 0.0))
        fig.tight_layout(rect=(0, 0.08, 1, 1))
    save(fig, outdir, outname)


def class_trends(designs, outdir, window=WINDOW, step=STEP, min_n=MIN_WINDOW):
    """One curve per class, from a window sliding along the depth axis."""
    depth = designs.depth.values
    centres = np.arange(depth.min() + window / 2.0,
                        depth.max() - window / 2.0 + step, step)
    keep, members = [], []
    for c in centres:
        sel = (depth >= c - window / 2.0) & (depth < c + window / 2.0)
        if sel.sum() >= min_n:
            keep.append(c)
            members.append(sel)
    print(f"class_trends: {len(keep)} windows of {window:g} A, "
          f"{min(s.sum() for s in members)}-{max(s.sum() for s in members)} "
          f"residues each")

    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    for name, _, colour in GROUPS:
        is_class = (designs.group == name).values
        ys, los, his = [], [], []
        for sel in members:
            k, n = int((is_class & sel).sum()), int(sel.sum())
            lo, hi = wilson(k, n)
            ys.append(100.0 * k / n)
            los.append(100.0 * lo)
            his.append(100.0 * hi)
        ax.fill_between(keep, los, his, color=colour, alpha=0.16, linewidth=0)
        ax.plot(keep, ys, lw=2.0, color=colour, label=name)
    ax.set_xlabel("Depth from the hairpin anchors (Å)")
    ax.set_ylabel("Frequency in design region 1 (%)")
    ax.set_ylim(bottom=0)
    ax.legend(ncol=3, fontsize=8.5, loc="upper center",
              bbox_to_anchor=(0.5, -0.16))
    fig.tight_layout()
    save(fig, outdir, "class_trends")


def depth_by_amino_acid(designs, wt, band, outdir):
    """No bins at all: the depths at which each amino acid was actually used."""
    order = [aa for _, aas, _ in GROUPS for aa in aas]
    used = [aa for aa in order if (designs.aa == aa).sum() > 0]
    used.sort(key=lambda aa: designs[designs.aa == aa].depth.median())

    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    draw_band(ax, band, where="bottom")
    data = [designs[designs.aa == aa].depth.values for aa in used]
    bp = ax.boxplot(data, orientation="horizontal", patch_artist=True,
                    widths=0.66,
                    medianprops=dict(color=INK, linewidth=1.3),
                    flierprops=dict(marker="o", markersize=2.2,
                                    markerfacecolor=MUTED,
                                    markeredgecolor="none", alpha=0.4))
    for patch, aa in zip(bp["boxes"], used):
        patch.set_facecolor(COLOUR_OF[CLASS_OF[aa]])
        patch.set_alpha(0.7)
        patch.set_edgecolor(INK)
        patch.set_linewidth(0.6)
    for w in bp["whiskers"] + bp["caps"]:
        w.set_color(MUTED)
        w.set_linewidth(0.8)
    # The wild-type residues sit on the same axis, one marker per residue, so
    # the native hairpin can be compared without it being a distribution.
    for aa in used:
        v = wt[wt.aa == aa].depth.values
        if len(v):
            ax.scatter(v, np.full(len(v), used.index(aa) + 1), marker="D",
                       s=22, color=NATIVE, zorder=5, edgecolor="white",
                       linewidth=0.5)
    for i, aa in enumerate(used):
        ax.annotate(f"n={(designs.aa == aa).sum()}", (designs.depth.max(), i + 1),
                    xytext=(6, 0), textcoords="offset points", va="center",
                    fontsize=7.5, color=MUTED)
    ax.set_yticks(range(1, len(used) + 1))
    ax.set_yticklabels(used)
    ax.set_ylim(0.2, len(used) + 0.8)
    ax.set_xlabel("Depth from the hairpin anchors (Å)")
    ax.set_ylabel("Amino acid, ordered by median depth")
    ax.set_xlim(0, designs.depth.max() * 1.13)
    handles = [Patch(facecolor=c, alpha=0.7, edgecolor=INK, linewidth=0.6,
                     label=n) for n, _, c in GROUPS]
    handles.append(Line2D([], [], marker="D", linestyle="none", color=NATIVE,
                          markersize=6, label="wild-type residue"))
    ax.legend(handles=handles, ncol=4, fontsize=8.5, loc="upper center",
              bbox_to_anchor=(0.5, -0.13))
    fig.tight_layout()
    save(fig, outdir, "depth_by_amino_acid")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default=os.path.join(HERE, "loop_depth.csv"))
    ap.add_argument("--outdir", default=os.path.join(HERE, "figure-options"))
    ap.add_argument("--width", type=float, default=BIN_WIDTH,
                    help="fixed bin width in angstroms")
    ap.add_argument("--origin", type=float, default=BIN_ORIGIN,
                    help="depth the fixed-width grid starts at")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    apply_style()

    designs, wt = load(args.csv)
    band = strand_coil_band(designs)
    print(f"{len(designs)} residues, strand/coil crossover {band[0]:.0f}-{band[1]:.0f} A")
    _, stats = binned(designs)
    print(stats.round(1).to_string())

    composition_stack(designs, band, args.outdir)
    composition_stack(designs, band, args.outdir, width=args.width,
                      origin=args.origin, outname="composition_stack_fixed")
    composition_stack(designs, band, args.outdir, width=args.width,
                      origin=args.origin, wt=wt, wt_style="dots",
                      outname="composition_stack_wt")
    composition_stack(designs, band, args.outdir, width=args.width,
                      origin=args.origin, wt=wt, wt_style="bars",
                      outname="composition_stack_wt_bars")
    class_trends(designs, args.outdir)
    depth_by_amino_acid(designs, wt, band, args.outdir)


if __name__ == "__main__":
    main()
