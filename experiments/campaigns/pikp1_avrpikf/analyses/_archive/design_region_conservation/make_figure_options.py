#!/usr/bin/env python3
"""Four candidate thesis figures for the design-region conservation result.

The result is that ProteinMPNN keeps the strand anchors at each end of design
region 1 and rebuilds the middle, and that the extra length is inserted at the
same two positions that diverge most. Each option below foregrounds a different
part of that.

    option1_bars        mean BLOSUM62 per position, bars shaded by how often
                        the position stayed in the alignment
    option2_lollipop    the same score as a stem plot with the aligned fraction
                        on its own track beneath, so the two are not conflated
    option3_heatstrip   both quantities as a two-row heat strip, compact enough
                        to sit beside another panel
    option4_annotated   design region 1 only, diverging bars with the anchors
                        and the insertion point called out

Usage:
    make_figure_options.py [--outdir figure-options]
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize

HERE = os.path.dirname(os.path.abspath(__file__))

INK, MUTED, GRID = "#222222", "#666666", "#DDDDDD"
LENGTH_COLOUR = {"shorter": "#D55E00", "equal": "#666666", "longer": "#0072B2"}
KEEP, LOSE = "#0072B2", "#D55E00"      # conserved / diverged
FILL = LinearSegmentedColormap.from_list(
    "keep_lose", [LOSE, "#F0F0F0", KEEP])


def apply_style():
    plt.rcParams.update({
        "figure.dpi": 130, "savefig.dpi": 300, "savefig.bbox": "tight",
        "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
        "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
        "xtick.color": MUTED, "ytick.color": MUTED,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
        "axes.axisbelow": True, "axes.spines.top": False,
        "axes.spines.right": False, "legend.frameon": False,
    })


def labels(b):
    return [f"{r}{p}" for r, p in zip(b.wt_residue, b.wt_position)]


def save(fig, outdir, name):
    path = os.path.join(outdir, f"{name}.png")
    fig.savefig(path)
    plt.close(fig)
    print("wrote", path)


def option1_bars(cons, outdir):
    """Score as bar height, aligned fraction as fill intensity."""
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.8),
                             gridspec_kw={"width_ratios": [13, 6]})
    norm = Normalize(0.6, 1.0)
    for ax, region in zip(axes, [1, 2]):
        b = cons[cons.design_region == region]
        shade = plt.get_cmap("Blues")(0.25 + 0.6 * norm(b.aligned_fraction))
        ax.bar(b.wt_position, b.mean_blosum62, color=shade,
               edgecolor="black", linewidth=0.6)
        ax.axhline(0, color=INK, lw=0.9)
        ax.set_xticks(b.wt_position)
        ax.set_xticklabels(labels(b), fontsize=8)
        ax.set_title(f"Design region {region}", fontsize=10)
        ax.set_xlabel("wild-type position")
    axes[0].set_ylabel("mean BLOSUM62 score")
    sm = plt.cm.ScalarMappable(cmap="Blues", norm=norm)
    cb = fig.colorbar(sm, ax=axes, shrink=0.75, pad=0.02)
    cb.set_label("designs with a residue at this position", fontsize=8)
    save(fig, outdir, "option1_bars")


def option2_lollipop(cons, outdir):
    """Score and aligned fraction on separate tracks, nothing conflated."""
    fig, axes = plt.subplots(
        2, 2, figsize=(9.6, 4.6), sharex="col",
        gridspec_kw={"width_ratios": [13, 6], "height_ratios": [2.4, 1]})
    for col, region in enumerate([1, 2]):
        b = cons[cons.design_region == region]
        top, bot = axes[0][col], axes[1][col]
        colour = [KEEP if v >= 0 else LOSE for v in b.mean_blosum62]
        top.vlines(b.wt_position, 0, b.mean_blosum62, color=colour, lw=2.2)
        top.scatter(b.wt_position, b.mean_blosum62, s=38, color=colour,
                    zorder=3, edgecolor="white", linewidth=0.8)
        top.axhline(0, color=INK, lw=0.9)
        top.set_title(f"Design region {region}", fontsize=10)
        bot.bar(b.wt_position, b.aligned_fraction, color=MUTED, width=0.65)
        bot.set_ylim(0, 1.05)
        bot.set_xticks(b.wt_position)
        bot.set_xticklabels(labels(b), fontsize=8)
        bot.set_xlabel("wild-type position")
    axes[0][0].set_ylabel("mean BLOSUM62")
    axes[1][0].set_ylabel("designs with\na residue here", fontsize=8)
    fig.tight_layout()
    save(fig, outdir, "option2_lollipop")


def option3_heatstrip(cons, outdir):
    """Compact two-row strip, for pairing with a structure panel."""
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 2.3),
                             gridspec_kw={"width_ratios": [13, 6]})
    for ax, region in zip(axes, [1, 2]):
        b = cons[cons.design_region == region]
        M = np.vstack([b.mean_blosum62.values / 4.0, b.aligned_fraction.values])
        ax.imshow(M, cmap=FILL, vmin=-1, vmax=1, aspect="auto")
        for j in range(M.shape[1]):
            ax.text(j, 0, f"{b.mean_blosum62.values[j]:+.1f}", ha="center",
                    va="center", fontsize=7)
            ax.text(j, 1, f"{b.aligned_fraction.values[j]:.2f}", ha="center",
                    va="center", fontsize=7)
        ax.set_yticks([0, 1])
        ax.set_yticklabels(["BLOSUM62", "has residue"], fontsize=8)
        ax.set_xticks(range(len(b)))
        ax.set_xticklabels(labels(b), fontsize=8)
        ax.set_title(f"Design region {region}", fontsize=10)
        ax.grid(False)
    fig.tight_layout()
    save(fig, outdir, "option3_heatstrip")


def option4_annotated(cons, outdir):
    """Design region 1 alone, with the reading written onto the figure."""
    b = cons[cons.design_region == 1]
    fig, ax = plt.subplots(figsize=(7.4, 4.2))
    colour = [KEEP if v >= 0 else LOSE for v in b.mean_blosum62]
    ax.bar(b.wt_position, b.mean_blosum62, color=colour, edgecolor="black",
           linewidth=0.6, width=0.72)
    ax.axhline(0, color=INK, lw=0.9)
    gapped = b[b.aligned_fraction < 0.9]
    ax.scatter(gapped.wt_position, [-2.3] * len(gapped), marker="v", s=44,
               color=MUTED, clip_on=False)
    for _, r in gapped.iterrows():
        ax.annotate(f"{r.aligned_fraction:.2f}", (r.wt_position, -2.55),
                    ha="center", fontsize=7.5, color=MUTED)
    top = b.nlargest(2, "mean_blosum62")
    for _, r in top.iterrows():
        ax.annotate(f"{r.mean_blosum62:+.2f}",
                    (r.wt_position, r.mean_blosum62 + 0.15),
                    ha="center", fontsize=8, color=KEEP, weight="bold")
    ax.set_xticks(b.wt_position)
    ax.set_xticklabels(labels(b), fontsize=9)
    ax.set_ylim(-2.7, 4.6)
    ax.set_ylabel("mean BLOSUM62 score")
    ax.set_xlabel("wild-type position in design region 1")
    # Placed low-left, where the negative bars leave the panel empty, so the
    # value annotations above the two anchors are not written over.
    ax.text(0.015, 0.20, "strand anchors kept", transform=ax.transAxes,
            fontsize=8.5, color=KEEP, va="top")
    ax.text(0.015, 0.145, "loop rebuilt", transform=ax.transAxes,
            fontsize=8.5, color=LOSE, va="top")
    ax.text(0.015, 0.09, "\u25bc  position the alignment gaps out",
            transform=ax.transAxes, fontsize=8, color=MUTED, va="top")
    fig.tight_layout()
    save(fig, outdir, "option4_annotated")


def option5_by_length(full, outdir):
    """The aligned fraction split by design length, which is what drives it.

    A design shorter than the wild-type 13 must lose a position somewhere, so
    the dip at Leu8 and Arg9 is mostly those 34 designs. Pooling them with the
    84 longer designs makes it look like a single insertion point, which it is
    not.
    """
    b = full[(full.design_region == 1) & (full.length_class != "all")]
    fig, axes = plt.subplots(1, 2, figsize=(9.8, 3.8), sharex=True)
    for cls in ["shorter", "equal", "longer"]:
        sub = b[b.length_class == cls].sort_values("wt_position")
        n = int(sub.n_designs.iloc[0])
        axes[0].plot(sub.wt_position, sub.aligned_fraction, marker="o", ms=4,
                     color=LENGTH_COLOUR[cls], lw=1.6,
                     label=f"{cls} than wild type (n={n})" if cls != "equal"
                     else f"equal to wild type (n={n})")
        axes[1].plot(sub.wt_position, sub.mean_blosum62, marker="o", ms=4,
                     color=LENGTH_COLOUR[cls], lw=1.6)
    ref = b[b.length_class == "longer"].sort_values("wt_position")
    for ax, ylab in [(axes[0], "designs with a residue here"),
                     (axes[1], "mean BLOSUM62 score")]:
        ax.set_xticks(ref.wt_position)
        ax.set_xticklabels(labels(ref), fontsize=8)
        ax.set_ylabel(ylab)
        ax.set_xlabel("wild-type position in design region 1")
    axes[1].axhline(0, color=INK, lw=0.9)
    axes[0].set_ylim(0, 1.05)
    axes[0].legend(fontsize=8, loc="lower left")
    fig.tight_layout()
    save(fig, outdir, "option5_by_length")


def option6_design_coordinates(prof, outdir):
    """Design region 1 in its own coordinates, which is where the extension is.

    The wild-type-indexed options cannot show the extra residues at all, since
    a 20-mer has 7 positions with no wild-type counterpart. Here every design
    position gets a bar for how many designs reach it, split by whether the
    residue there matches a wild-type position or is an insertion.
    """
    fig, axes = plt.subplots(2, 1, figsize=(8.0, 5.0), sharex=True,
                             gridspec_kw={"height_ratios": [2.0, 1.0]})
    x = prof.design_position
    axes[0].bar(x, prof.n_matched, color=KEEP, edgecolor="black",
                linewidth=0.6, label="aligns to a wild-type residue")
    axes[0].bar(x, prof.n_inserted, bottom=prof.n_matched, color=LOSE,
                edgecolor="black", linewidth=0.6, label="inserted residue")
    for xi, n in zip(x, prof.n_designs):
        axes[0].annotate(str(int(n)), (xi, n + 2.5), ha="center", fontsize=7.5,
                         color=MUTED)
    axes[0].set_ylabel("designs with a residue here")
    axes[0].set_ylim(0, 140)
    axes[0].legend(fontsize=8.5, loc="upper right")
    axes[1].bar(x, prof.entropy_bits, color=MUTED, width=0.65)
    axes[1].axhline(np.log2(20), ls="--", color=INK, lw=0.9)
    axes[1].annotate("uniform over 20 residues", (0.99, np.log2(20)),
                     xycoords=("axes fraction", "data"), ha="right",
                     va="bottom", fontsize=7.5, color=MUTED)
    axes[1].set_ylabel("entropy (bits)")
    axes[1].set_xlabel("position in design region 1, as built")
    axes[1].set_xticks(x)
    fig.tight_layout()
    save(fig, outdir, "option6_design_coordinates")


def option7_structural_matrix(long_df, outdir):
    """Every design as a row, so the gap in the middle is literally visible.

    Sequence alignment could not place the insertion consistently, because
    ProteinMPNN recovers little of the wild-type sequence in this region. This
    uses the structural correspondence instead: each design is superimposed on
    the wild type by its fixed scaffold, and a residue counts as matched when a
    wild-type Ca sits within 2.5 A of it. Rows are sorted by design region
    length, so the widening band of unmatched residues in the middle is the
    extension.
    """
    piv = (long_df.pivot_table(index=["region_length", "design"],
                               columns="design_position", values="matched")
                  .sort_index(level=["region_length", "design"]))
    M = np.full(piv.shape, np.nan)
    M[piv.notna().values] = piv.values[piv.notna().values].astype(float)

    fig, ax = plt.subplots(figsize=(7.6, 6.2))
    cmap = matplotlib.colors.ListedColormap([LOSE, KEEP])
    cmap.set_bad("#F2F2F2")
    ax.imshow(np.ma.masked_invalid(M), cmap=cmap, vmin=0, vmax=1,
              aspect="auto", interpolation="nearest")
    ax.set_xticks(range(piv.shape[1]))
    ax.set_xticklabels(piv.columns, fontsize=8)
    ax.set_xlabel("position in design region 1, as built")
    lengths = [i[0] for i in piv.index]
    ticks = [lengths.index(v) for v in sorted(set(lengths))]
    ax.set_yticks(ticks)
    ax.set_yticklabels(sorted(set(lengths)), fontsize=8)
    ax.set_ylabel("design region 1 length (one row per backbone)")
    ax.grid(False)
    handles = [Patch(facecolor=KEEP, edgecolor="black", linewidth=0.4,
                     label="matches a wild-type residue"),
               Patch(facecolor=LOSE, edgecolor="black", linewidth=0.4,
                     label="no wild-type counterpart"),
               Patch(facecolor="#F2F2F2", edgecolor="black", linewidth=0.4,
                     label="design is shorter than this")]
    ax.legend(handles=handles, fontsize=8.5, ncol=3, loc="lower center",
              bbox_to_anchor=(0.5, 1.01))
    fig.tight_layout()
    save(fig, outdir, "option7_structural_matrix")


def option8_structural_alignment(long_df, outdir):
    """The structural correspondence laid out as a real alignment.

    option7 answers only "does this residue have a wild-type counterpart", so
    every row starts flush left and the unmatched residues trail off as a
    staircase. Here each residue sits in the column of the wild-type residue it
    is superimposed on, so structurally equivalent residues line up.

    Two things decide the layout, and both were wrong in the first version.

    Wild-type positions 7, 8 and 9 are matched by no design at all, because
    they are the wild-type loop and every design replaces it. Drawn as ordinary
    empty columns they read as holes splitting the middle of the figure, so
    they get their own colour: present in the wild type, absent from the
    designs.

    Insertions are CENTRED in their block rather than left-aligned. A design
    with 4 insertions and one with 10 both belong in the middle of the rebuilt
    region, and left-aligning them put the short ones hard against the
    left-hand block and left a ragged hole on the right.
    """
    rows, max_ins = {}, 0
    for design, sub in long_df.groupby("design"):
        sub = sub.sort_values("design_position")
        matched, inserted = [], 0
        for r in sub.itertuples():
            if not np.isnan(r.wt_position):
                matched.append(int(r.wt_position))
            else:
                inserted += 1
        rows[design] = (matched, inserted, int(sub.region_length.iloc[0]))
        max_ins = max(max_ins, inserted)

    # The two strands are the wild-type positions nearly every design lands
    # on. 1-5 and 11-13 are matched by all 64, 6 by 46 and 10 by 60, while
    # 7-9 are matched by only 3 to 8. Splitting on a simple majority puts the
    # rebuilt loop in the middle without hard-coding which positions it is.
    n_wt = int(long_df.wt_position.max())
    n_designs = long_df.design.nunique()
    hits = (long_df.wt_position.dropna().astype(int)
            .value_counts().reindex(range(1, n_wt + 1), fill_value=0))
    core = [w for w in range(1, n_wt + 1) if hits[w] >= 0.5 * n_designs]
    absent = [w for w in range(1, n_wt + 1) if w not in core]
    split = min(absent) if absent else n_wt
    left = [w for w in core if w < split]
    right = [w for w in core if w > split]

    # columns: matched wild-type block, the rebuilt middle, matched block
    cols = ([("wt", w) for w in left]
            + [("wtonly", w) for w in absent]
            + [("ins", k) for k in range(max_ins)]
            + [("wt", w) for w in right])
    index = {c: i for i, c in enumerate(cols)}

    order = sorted(rows, key=lambda d: (rows[d][2], d))
    M = np.full((len(order), len(cols)), np.nan)
    for i, design in enumerate(order):
        matched, inserted, _ = rows[design]
        for w in matched:
            key = ("wtonly", w) if w in absent else ("wt", w)
            M[i, index[key]] = 1.0
        start = (max_ins - inserted) // 2          # centred, not left-aligned
        for k in range(inserted):
            M[i, index[("ins", start + k)]] = 0.0
    # The wild type has these residues whatever the design does, so the column
    # is filled for every row unless that design actually lands on it.
    for w in absent:
        col = index[("wtonly", w)]
        M[np.isnan(M[:, col]), col] = 0.5

    fig, ax = plt.subplots(figsize=(10.0, 6.4))
    cmap = matplotlib.colors.ListedColormap([LOSE, "#E8D7C3", KEEP])
    cmap.set_bad("#FFFFFF")
    ax.imshow(np.ma.masked_invalid(M), cmap=cmap, vmin=0, vmax=1,
              aspect="auto", interpolation="nearest")
    ticks = [index[("wt", w)] for w in left + right]
    ticks += [index[("wtonly", w)] for w in absent]
    labels = [str(w) for w in left + right] + [str(w) for w in absent]
    o = np.argsort(ticks)
    ax.set_xticks([ticks[k] for k in o])
    ax.set_xticklabels([labels[k] for k in o], fontsize=8)
    ax.set_xlabel("wild-type position in design region 1, with the rebuilt "
                  "middle between")
    lengths = [rows[d][2] for d in order]
    yt = [lengths.index(v) for v in sorted(set(lengths))]
    ax.set_yticks(yt)
    ax.set_yticklabels(sorted(set(lengths)), fontsize=8)
    ax.set_ylabel("design region 1 length (one row per backbone)")
    ax.grid(False)
    handles = [Patch(facecolor=KEEP, edgecolor="black", linewidth=0.4,
                     label="design residue on a wild-type residue"),
               Patch(facecolor="#E8D7C3", edgecolor="black", linewidth=0.4,
                     label="wild-type loop residue, unmatched"),
               Patch(facecolor=LOSE, edgecolor="black", linewidth=0.4,
                     label="inserted design residue"),
               Patch(facecolor="#FFFFFF", edgecolor="black", linewidth=0.4,
                     label="design is shorter than this")]
    ax.legend(handles=handles, fontsize=8.5, ncol=2, loc="lower center",
              bbox_to_anchor=(0.5, 1.01))
    fig.tight_layout()
    save(fig, outdir, "option8_structural_alignment")


def option9_hairpin_layout(lay, outdir):
    """Design region 1 laid out by depth down each arm of the hairpin.

    option8 anchored the rows to the wild-type sequence and had to decide where
    the insertions went, which is a layout choice rather than a measurement.
    Here there is no choice left. The region is a U with both ends fixed, so
    residues index inward from the left anchor down the N arm and inward from
    the right anchor down the C arm. Every row is as wide as the longest design
    and the hole in its middle is exactly how many residues shorter it is,
    which is the only thing RFdiffusion varied.

    The turn comes from the design's own geometry, the residue whose Ca is
    furthest from the line joining the two flanking anchors. Nothing about the
    layout depends on the wild type. Colour does: it marks which residues have
    a wild-type counterpart, so the rebuilt part of each hairpin is visible
    against a layout that did not use that information.
    """
    # Length first, then where the rebuilt block starts and ends. Ordering
    # only by length left the rows inside each band in file order, which made
    # the edges of the rebuilt region look ragged for no reason.
    key = []
    for design, sub in lay.groupby("design"):
        loop = sub[sub.ss != "strand"].column
        key.append((design, int(sub.region_length.iloc[0]),
                    int(loop.min()) if len(loop) else 99,
                    -int(loop.max()) if len(loop) else 0))
    order = [k[0] for k in sorted(key, key=lambda k: (k[1], k[2], k[3]))]
    ncol = int(lay.column.max())
    M = np.full((len(order), ncol), np.nan)
    turns = []
    for i, design in enumerate(order):
        sub = lay[lay.design == design]
        for r in sub.itertuples():
            M[i, int(r.column) - 1] = 1.0 if bool(r.matched) else 0.0
            if r.is_turn:
                turns.append((int(r.column) - 1, i))

    fig, ax = plt.subplots(figsize=(9.0, 6.0))
    cmap = matplotlib.colors.ListedColormap([LOSE, KEEP])
    cmap.set_bad("#FFFFFF")
    ax.imshow(np.ma.masked_invalid(M), cmap=cmap, vmin=0, vmax=1,
              aspect="auto", interpolation="nearest")
    tx, ty = zip(*turns)
    ax.scatter(tx, ty, s=6, color=INK, marker="|", linewidths=0.8)
    ax.set_xticks(range(ncol))
    ax.set_xticklabels(range(1, ncol + 1), fontsize=8)
    ax.set_xlabel("position, counted in from each anchor of the hairpin")
    lengths = [lay[lay.design == d].region_length.iloc[0] for d in order]
    # Rule between each design-length band, so the rows a tick label refers to
    # are bounded rather than left to be counted.
    for i in range(1, len(lengths)):
        if lengths[i] != lengths[i - 1]:
            ax.axhline(i - 0.5, color=INK, lw=0.6, ls=(0, (3, 3)))
    ticks = [lengths.index(v) for v in sorted(set(lengths))]
    ax.set_yticks(ticks)
    ax.set_yticklabels(sorted(set(lengths)), fontsize=8)
    ax.set_ylabel("design region 1 length (one row per backbone)")
    wt_row = order.index("wild-type") if "wild-type" in order else None
    if wt_row is not None:
        ax.annotate("wild type", xy=(1, wt_row), xycoords=("axes fraction",
                    "data"), xytext=(6, 0), textcoords="offset points",
                    ha="left", va="center", fontsize=8, color="#5A5A5A")
    ax.grid(False)
    handles = [Patch(facecolor=KEEP, edgecolor="black", linewidth=0.4,
                     label="has a wild-type counterpart"),
               Patch(facecolor=LOSE, edgecolor="black", linewidth=0.4,
                     label="rebuilt, no wild-type counterpart"),
               Patch(facecolor="#FFFFFF", edgecolor="black", linewidth=0.4,
                     label="shorter than the longest design"),
               Line2D([], [], color=INK, marker="|", linestyle="none",
                      markersize=6, label="turn of the hairpin")]
    ax.legend(handles=handles, fontsize=8.5, ncol=4, loc="lower center",
              bbox_to_anchor=(0.5, 1.01))
    fig.tight_layout()
    save(fig, outdir, "option9_hairpin_layout")


def option9_hairpin_layout(lay, outdir):
    """Design region 1 laid out by depth down each arm of the hairpin.

    option8 anchored the rows to the wild-type sequence and had to decide where
    the insertions went, which is a layout choice rather than a measurement.
    Here there is no choice left. The region is a U with both ends fixed, so
    residues index inward from the left anchor down the N arm and inward from
    the right anchor down the C arm. Every row is as wide as the longest design
    and the hole in its middle is exactly how many residues shorter it is,
    which is the only thing RFdiffusion varied.

    The turn comes from the design's own geometry, the residue whose Ca is
    furthest from the line joining the two flanking anchors. Nothing about the
    layout depends on the wild type. Colour does: it marks which residues have
    a wild-type counterpart, so the rebuilt part of each hairpin is visible
    against a layout that did not use that information.
    """
    # Length first, then where the rebuilt block starts and ends. Ordering
    # only by length left the rows inside each band in file order, which made
    # the edges of the rebuilt region look ragged for no reason.
    key = []
    for design, sub in lay.groupby("design"):
        loop = sub[sub.ss != "strand"].column
        key.append((design, int(sub.region_length.iloc[0]),
                    int(loop.min()) if len(loop) else 99,
                    -int(loop.max()) if len(loop) else 0))
    order = [k[0] for k in sorted(key, key=lambda k: (k[1], k[2], k[3]))]
    ncol = int(lay.column.max())
    M = np.full((len(order), ncol), np.nan)
    turns = []
    for i, design in enumerate(order):
        sub = lay[lay.design == design]
        for r in sub.itertuples():
            M[i, int(r.column) - 1] = 1.0 if bool(r.matched) else 0.0
            if r.is_turn:
                turns.append((int(r.column) - 1, i))

    fig, ax = plt.subplots(figsize=(9.0, 6.0))
    cmap = matplotlib.colors.ListedColormap([LOSE, KEEP])
    cmap.set_bad("#FFFFFF")
    ax.imshow(np.ma.masked_invalid(M), cmap=cmap, vmin=0, vmax=1,
              aspect="auto", interpolation="nearest")
    tx, ty = zip(*turns)
    ax.scatter(tx, ty, s=6, color=INK, marker="|", linewidths=0.8)
    ax.set_xticks(range(ncol))
    ax.set_xticklabels(range(1, ncol + 1), fontsize=8)
    ax.set_xlabel("position, counted in from each anchor of the hairpin")
    lengths = [lay[lay.design == d].region_length.iloc[0] for d in order]
    # Rule between each design-length band, so the rows a tick label refers to
    # are bounded rather than left to be counted.
    for i in range(1, len(lengths)):
        if lengths[i] != lengths[i - 1]:
            ax.axhline(i - 0.5, color=INK, lw=0.6, ls=(0, (3, 3)))
    ticks = [lengths.index(v) for v in sorted(set(lengths))]
    ax.set_yticks(ticks)
    ax.set_yticklabels(sorted(set(lengths)), fontsize=8)
    ax.set_ylabel("design region 1 length (one row per backbone)")
    wt_row = order.index("wild-type") if "wild-type" in order else None
    if wt_row is not None:
        ax.annotate("wild type", xy=(1, wt_row), xycoords=("axes fraction",
                    "data"), xytext=(6, 0), textcoords="offset points",
                    ha="left", va="center", fontsize=8, color="#5A5A5A")
    ax.grid(False)
    handles = [Patch(facecolor=KEEP, edgecolor="black", linewidth=0.4,
                     label="has a wild-type counterpart"),
               Patch(facecolor=LOSE, edgecolor="black", linewidth=0.4,
                     label="rebuilt, no wild-type counterpart"),
               Patch(facecolor="#FFFFFF", edgecolor="black", linewidth=0.4,
                     label="shorter than the longest design"),
               Line2D([], [], color=INK, marker="|", linestyle="none",
                      markersize=6, label="turn of the hairpin")]
    ax.legend(handles=handles, fontsize=8.5, ncol=4, loc="lower center",
              bbox_to_anchor=(0.5, 1.01))
    fig.tight_layout()
    save(fig, outdir, "option9_hairpin_layout")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=os.path.join(HERE, "figure-options"))
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    apply_style()
    full = pd.read_csv(os.path.join(HERE, "design_region_conservation.csv"))
    cons = full[full.length_class == "all"]
    option1_bars(cons, args.outdir)
    option2_lollipop(cons, args.outdir)
    option3_heatstrip(cons, args.outdir)
    option4_annotated(cons, args.outdir)
    option5_by_length(full, args.outdir)
    prof = pd.read_csv(os.path.join(HERE, "design_coordinate_profile.csv"))
    option6_design_coordinates(prof, args.outdir)
    long_df = pd.read_csv(os.path.join(
        HERE, "structural_correspondence_per_design.csv"))
    option7_structural_matrix(long_df, args.outdir)
    option8_structural_alignment(long_df, args.outdir)
    lay = pd.read_csv(os.path.join(HERE, "hairpin_layout.csv"))
    option9_hairpin_layout(lay, args.outdir)
    lay = pd.read_csv(os.path.join(HERE, "hairpin_layout.csv"))
    option9_hairpin_layout(lay, args.outdir)








def option9_hairpin_layout(lay, outdir):
    """Design region 1 laid out by depth down each arm of the hairpin.

    option8 anchored the rows to the wild-type sequence and had to decide where
    the insertions went, which is a layout choice rather than a measurement.
    Here there is no choice left. The region is a U with both ends fixed, so
    residues index inward from the left anchor down the N arm and inward from
    the right anchor down the C arm. Every row is as wide as the longest design
    and the hole in its middle is exactly how many residues shorter it is,
    which is the only thing RFdiffusion varied.

    The turn comes from the design's own pairing register, so nothing about the
    layout depends on the wild type.

    Colour is DSSP secondary structure rather than wild-type correspondence.
    Correspondence was misleading here: it is a property of the comparison, not
    of the design, and its 2.5 A nearest-neighbour rule scattered isolated
    matches through the rebuilt loop that were coincidences of proximity rather
    than real equivalences. Strand against loop is a property of the design
    itself, computed from its own backbone hydrogen bonds.
    """
    # Length first, then where the rebuilt block starts and ends. Ordering
    # only by length left the rows inside each band in file order, which made
    # the edges of the rebuilt region look ragged for no reason.
    key = []
    for design, sub in lay.groupby("design"):
        loop = sub[sub.ss != "strand"].column
        key.append((design, int(sub.region_length.iloc[0]),
                    int(loop.min()) if len(loop) else 99,
                    -int(loop.max()) if len(loop) else 0))
    order = [k[0] for k in sorted(key, key=lambda k: (k[1], k[2], k[3]))]
    ncol = int(lay.column.max())
    M = np.full((len(order), ncol), np.nan)
    # The wild type is laid out by the same rules and sorted in with the other
    # 13-residue rows rather than pinned to the top, so it reads as one of the
    # set. Grey only marks which row it is.
    for i, design in enumerate(order):
        sub = lay[lay.design == design]
        wt = design == "wild-type"
        for r in sub.itertuples():
            strand = r.ss == "strand"
            M[i, int(r.column) - 1] = (2 + strand) if wt else strand

    fig, ax = plt.subplots(figsize=(9.0, 6.0))
    WT_LOOP, WT_STRAND = "#C8C8C8", "#5A5A5A"
    cmap = matplotlib.colors.ListedColormap([LOSE, KEEP, WT_LOOP, WT_STRAND])
    cmap.set_bad("#FFFFFF")
    norm = matplotlib.colors.BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], cmap.N)
    ax.imshow(np.ma.masked_invalid(M), cmap=cmap, norm=norm,
              aspect="auto", interpolation="nearest")
    ax.set_xticks(range(ncol))
    ax.set_xticklabels(range(1, ncol + 1), fontsize=8)
    ax.set_xlabel("position, counted in from each anchor of the hairpin")
    lengths = [lay[lay.design == d].region_length.iloc[0] for d in order]
    # Rule between each design-length band, so the rows a tick label refers to
    # are bounded rather than left to be counted.
    for i in range(1, len(lengths)):
        if lengths[i] != lengths[i - 1]:
            ax.axhline(i - 0.5, color=INK, lw=0.6, ls=(0, (3, 3)))
    ticks = [lengths.index(v) for v in sorted(set(lengths))]
    ax.set_yticks(ticks)
    ax.set_yticklabels(sorted(set(lengths)), fontsize=8)
    ax.set_ylabel("design region 1 length (one row per backbone)")
    wt_row = order.index("wild-type") if "wild-type" in order else None
    if wt_row is not None:
        ax.annotate("wild type", xy=(1, wt_row), xycoords=("axes fraction",
                    "data"), xytext=(6, 0), textcoords="offset points",
                    ha="left", va="center", fontsize=8, color="#5A5A5A")
    ax.grid(False)
    handles = [Patch(facecolor=KEEP, edgecolor="black", linewidth=0.4,
                     label=r"$\beta$-strand"),
               Patch(facecolor=LOSE, edgecolor="black", linewidth=0.4,
                     label="loop"),
               Patch(facecolor=WT_STRAND, edgecolor="black", linewidth=0.4,
                     label=r"wild-type $\beta$-strand"),
               Patch(facecolor=WT_LOOP, edgecolor="black", linewidth=0.4,
                     label="wild-type loop"),
               Patch(facecolor="#FFFFFF", edgecolor="black", linewidth=0.4,
                     label="shorter than the longest design")]
    ax.legend(handles=handles, fontsize=8.5, ncol=5, loc="lower center",
              bbox_to_anchor=(0.5, 1.01))
    fig.tight_layout()
    save(fig, outdir, "option9_hairpin_layout")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=os.path.join(HERE, "figure-options"))
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    apply_style()
    full = pd.read_csv(os.path.join(HERE, "design_region_conservation.csv"))
    cons = full[full.length_class == "all"]
    option1_bars(cons, args.outdir)
    option2_lollipop(cons, args.outdir)
    option3_heatstrip(cons, args.outdir)
    option4_annotated(cons, args.outdir)
    option5_by_length(full, args.outdir)
    prof = pd.read_csv(os.path.join(HERE, "design_coordinate_profile.csv"))
    option6_design_coordinates(prof, args.outdir)
    long_df = pd.read_csv(os.path.join(
        HERE, "structural_correspondence_per_design.csv"))
    option7_structural_matrix(long_df, args.outdir)
    option8_structural_alignment(long_df, args.outdir)
    lay = pd.read_csv(os.path.join(HERE, "hairpin_layout.csv"))
    option9_hairpin_layout(lay, args.outdir)
    lay = pd.read_csv(os.path.join(HERE, "hairpin_layout.csv"))
    option9_hairpin_layout(lay, args.outdir)







if __name__ == "__main__":
    main()
