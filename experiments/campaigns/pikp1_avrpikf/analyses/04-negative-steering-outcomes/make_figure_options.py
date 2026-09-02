#!/usr/bin/env python3
"""Candidate thesis figures for the AVR-Pia contact result.

The paragraph's claim is that negative steering REMOVES AVR-Pia contact. An
earlier version of this script plotted the proportion of incorrect poses still
touching AVR-Pia, which foregrounds the residue that survives rather than the
reduction, and its largest bar was a post-steering one. These options put the
reduction first.

The reduction, over all 128 designs:

    total AVR-Pia-specific contacts   207 -> 111
    designs contacting the surface     90 ->  48
    cleared entirely                   52
    gained contact                     10
    median contacts per design          2 ->   0

    option4_distribution    how many designs contact how many residues, before
                            and after steering
    option5_scatter         pose error against contacts, one point per design

Usage:
    make_figure_options.py [--outdir figure-options]
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))

INK, MUTED, GRID = "#222222", "#666666", "#DDDDDD"
COLD, STEERED = "#999999", "#0072B2"
REMOVED, KEPT, GAINED = "#0072B2", "#D55E00", "#8C564B"
NEVER = "#CCCCCC"
RA_EFF_THRESHOLD = 5.0
# Same tick set as the benchmarking chapter's log axes.
RA_TICKS = [0.5, 1, 2, 5, 10, 20, 40]


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


def readable_log(axis, ticks):
    """Plain-number ticks on a log axis; the decade defaults are unreadable."""
    axis.set_ticks(ticks)
    axis.set_ticklabels([f"{v:g}" for v in ticks])
    axis.set_minor_locator(plt.NullLocator())


def outcome(row):
    c, s = row.cold_n_pia_specific, row.steered_n_pia_specific
    if c > 0 and s == 0:
        return "cleared"
    if c > 0 and s > 0:
        return "reduced" if s < c else "still contacting"
    if c == 0 and s > 0:
        return "gained"
    return "never contacted"


def save(fig, outdir, name):
    path = os.path.join(outdir, f"{name}.png")
    fig.savefig(path)
    plt.close(fig)
    print("wrote", path)


def option4_distribution(df, outdir):
    """How many designs contact how many residues, before and after."""
    fig, ax = plt.subplots(figsize=(6.8, 4.2))
    hi = int(df[["cold_n_pia_specific", "steered_n_pia_specific"]].max().max())
    bins = np.arange(-0.5, hi + 1.5)
    counts, _, patches = ax.hist(
        [df.cold_n_pia_specific, df.steered_n_pia_specific], bins=bins,
        color=[COLD, STEERED], edgecolor="black", linewidth=0.6,
        label=["cold start", "after steering"])
    # Counts on the bars, so the shift is readable without measuring heights.
    for series in patches:
        for bar in series:
            h = int(bar.get_height())
            if h:
                ax.annotate(str(h),
                            (bar.get_x() + bar.get_width() / 2, h),
                            xytext=(0, 2), textcoords="offset points",
                            ha="center", va="bottom", fontsize=7.5,
                            color=INK)
    ax.set_xticks(np.arange(0, hi + 1))
    ax.set_xlabel("AVR-Pia interface-specific residues contacted")
    ax.set_ylabel("Frequency")
    ax.set_ylim(0, max(counts.max() for counts in [counts[0], counts[1]]) * 1.12)
    ax.legend(fontsize=9, ncol=2, loc="lower center",
              bbox_to_anchor=(0.5, 1.01))
    save(fig, outdir, "option4_distribution")


def option5_scatter(df, outdir):
    """Per design, with the log axis given plain-number ticks."""
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.3), sharey=True)
    for ax, (stage, label, colour) in zip(
            axes, [("cold", "cold start", COLD),
                   ("steered", "after steering", STEERED)]):
        # Points overlap heavily at integer contact counts, so identical
        # (RMSD band, count) cells are drawn once at a size set by how many
        # designs fall in them. Plain jitter hid the mass.
        d = df.copy()
        d["band"] = pd.cut(d[f"{stage}_ra_eff"], bins=np.geomspace(0.4, 45, 16))
        g = (d.groupby(["band", f"{stage}_n_pia_specific"], observed=True)
               .size().reset_index(name="n"))
        g["x"] = [b.mid for b in g.band]
        ax.scatter(g.x, g[f"{stage}_n_pia_specific"], s=14 + 26 * g.n,
                   alpha=0.65, color=colour, edgecolor="black", linewidth=0.4)
        ax.axvline(RA_EFF_THRESHOLD, ls="--", color=KEPT, lw=1.2)
        ax.set_xscale("log")
        readable_log(ax.xaxis, RA_TICKS)
        ax.set_xlabel("receptor-aligned effector RMSD (Å)")
        ax.set_title(label, fontsize=10)
    axes[0].set_ylabel("AVR-Pia-specific residues contacted")
    save(fig, outdir, "option5_scatter")


STAGES2 = [("cold", "Initial", COLD), ("steered", "Steered", STEERED)]


def _counts(df, stage, hi):
    """(correctly posed, incorrectly posed) counts per contact value."""
    ok = df[df[f"{stage}_correct"].astype(bool)][f"{stage}_n_pia_specific"]
    bad = df[~df[f"{stage}_correct"].astype(bool)][f"{stage}_n_pia_specific"]
    return (np.array([int((ok == v).sum()) for v in range(hi + 1)]),
            np.array([int((bad == v).sum()) for v in range(hi + 1)]))


def _hi(df):
    return int(df[["cold_n_pia_specific", "steered_n_pia_specific"]].max().max())


def option6_stacked_by_pose(df, outdir):
    """Grouped bars, each split at the 5 A pose threshold.

    Every segment carries its own count rather than only the column total, so
    the pose composition can be read without measuring segment heights. Zero
    segments are left unlabelled.
    """
    hi = _hi(df)
    x = np.arange(hi + 1)
    width = 0.4
    fig, ax = plt.subplots(figsize=(7.4, 4.4))
    for k, (stage, label, base) in enumerate(STAGES2):
        h_ok, h_bad = _counts(df, stage, hi)
        xs = x + (k - 0.5) * width
        ax.bar(xs, h_ok, width, color=base, edgecolor="black", linewidth=0.6,
               label=f"{label}, correctly posed")
        ax.bar(xs, h_bad, width, bottom=h_ok, color=base, alpha=0.42,
               edgecolor="black", linewidth=0.6,
               label=f"{label}, incorrectly posed")
        # Each count sits just above the top edge of its own segment rather
        # than inside it, so a thin segment still gets a readable label.
        for xi, a, b in zip(xs, h_ok, h_bad):
            if a:
                ax.annotate(str(a), (xi, a), xytext=(0, 2),
                            textcoords="offset points", ha="center",
                            va="bottom", fontsize=7.5, color=INK)
            if b:
                ax.annotate(str(b), (xi, a + b), xytext=(0, 2),
                            textcoords="offset points", ha="center",
                            va="bottom", fontsize=7.5, color=INK)
    ax.set_xticks(x)
    # No gridlines at all. Every segment carries its own count, so a rule
    # only competes with the numbers it would otherwise help estimate.
    ax.grid(False)
    ax.set_xlabel("AVR-Pia interface-specific residues contacted")
    ax.set_ylabel("Frequency")
    # Head room so the topmost label is not clipped by the axes.
    ax.set_ylim(0, ax.get_ylim()[1] * 1.08)
    ax.legend(fontsize=8.5, ncol=2, loc="lower center",
              bbox_to_anchor=(0.5, 1.01))
    save(fig, outdir, "option6_stacked_by_pose")


def option9_rmsd_by_contacts(df, outdir):
    """Pose error as a distribution, split by how much AVR-Pia is contacted.

    Reads the relationship the other way round from option 5. If contacting
    the competing surface goes with being badly posed, these should climb.
    """
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.2), sharey=True)
    for ax, (stage, title, colour) in zip(
            axes, [("cold", "cold start", COLD),
                   ("steered", "after steering", STEERED)]):
        groups, labels = [], []
        for v in sorted(df[f"{stage}_n_pia_specific"].unique()):
            sel = df[df[f"{stage}_n_pia_specific"] == v][f"{stage}_ra_eff"]
            if len(sel) >= 3:
                groups.append(sel.values)
                labels.append(f"{int(v)}\nn={len(sel)}")
        bp = ax.boxplot(groups, patch_artist=True, widths=0.6,
                        medianprops=dict(color=INK, linewidth=1.5))
        for patch in bp["boxes"]:
            patch.set_facecolor(colour)
            patch.set_alpha(0.65)
        ax.set_xticklabels(labels, fontsize=8)
        ax.axhline(RA_EFF_THRESHOLD, ls="--", color=KEPT, lw=1.2)
        ax.set_yscale("log")
        readable_log(ax.yaxis, RA_TICKS)
        ax.set_xlabel("AVR-Pia interface-specific residues contacted")
        ax.set_title(title, fontsize=10)
    axes[0].set_ylabel("receptor-aligned effector RMSD (Å)")
    save(fig, outdir, "option9_rmsd_by_contacts")


def option15_change(df, outdir):
    """The change between Initial and Steered, rather than both levels.

    option14 shows two distributions and leaves the reader to subtract them.
    Here the bar IS the subtraction: how many more or fewer designs sit at each
    contact count after steering. Green above the axis is a gain, red below is
    a loss, so the result reads as designs draining out of the contacting rows
    and piling up at zero.
    """
    hi = _hi(df)
    x = np.arange(hi + 1)
    GAIN, LOSS = "#1B7F4B", "#C0392B"
    fs_val, fs_lab, fs_title, lw = 6.5, 8.0, 9.0, 0.45
    fig, axes = plt.subplots(2, 1, figsize=(6.8, 4.4), sharex=True)
    for ax, (correct, title) in zip(
            axes, [(True, "Correctly posed"), (False, "Incorrectly posed")]):
        cold_ok, cold_bad = _counts(df, "cold", hi)
        st_ok, st_bad = _counts(df, "steered", hi)
        delta = (st_ok - cold_ok) if correct else (st_bad - cold_bad)
        ax.bar(x, delta, 0.62, color=[GAIN if v >= 0 else LOSS for v in delta],
               edgecolor="black", linewidth=lw)
        for xi, v in zip(x, delta):
            if v:
                ax.annotate(f"{v:+d}", (xi, v),
                            xytext=(0, 2 if v > 0 else -2),
                            textcoords="offset points", ha="center",
                            va="bottom" if v > 0 else "top",
                            fontsize=fs_val, color=INK)
        ax.axhline(0, color=INK, lw=0.8)
        ax.set_title(title, fontsize=fs_title)
        ax.set_xticks(x)
        ax.tick_params(labelsize=fs_lab - 0.5, width=0.6, length=3)
        ax.grid(axis="x", visible=False)
        for sp in ax.spines.values():
            sp.set_linewidth(0.6)
        pad = max(abs(delta).max(), 1) * 0.28
        ax.set_ylim(delta.min() - pad, delta.max() + pad)
    axes[-1].set_xlabel("AVR-Pia interface-specific residues contacted",
                        fontsize=fs_lab)
    fig.supylabel("Change in number of designs after steering", fontsize=fs_lab)
    handles = [Patch(facecolor=GAIN, edgecolor="black", linewidth=lw,
                     label="more designs after steering"),
               Patch(facecolor=LOSS, edgecolor="black", linewidth=lw,
                     label="fewer designs after steering")]
    fig.legend(handles=handles, fontsize=fs_lab, ncol=2, loc="lower center",
               bbox_to_anchor=(0.5, 0.0), frameon=False)
    fig.tight_layout(rect=(0.02, 0.06, 1, 1))
    save(fig, outdir, "option15_change")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=os.path.join(HERE, "figure-options"))
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    apply_style()
    df = pd.read_csv(os.path.join(HERE, "campaign_interface_check.csv"))
    option4_distribution(df, args.outdir)
    option5_scatter(df, args.outdir)
    option6_stacked_by_pose(df, args.outdir)
    option15_change(df, args.outdir)
    option9_rmsd_by_contacts(df, args.outdir)



if __name__ == "__main__":
    main()
