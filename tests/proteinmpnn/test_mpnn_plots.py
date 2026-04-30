#!/usr/bin/env python3
"""
test_mpnn_plots.py
------------------
Standalone iteration harness for ``mpnn_plots.py``.

Reads the QC metadata CSV produced by a previous ``test_proteinmpnn``
run and re-emits the diagnostic plots in an output directory of
choice.  No GPU, no Nextflow — designed for fast iteration of the
plot code itself against cached test outputs.

Plots produced
--------------

mpnn_score_distribution.png
    Two stacked per-design boxplots, both ranked by **median
    design-region score** (this iteration: was global before).
    Top:    design-region score (the newly designed bit — primary).
    Bottom: global score (sanity check that global also looks fine).

mpnn_sequence_diversity.png
    MMseqs2 cluster count vs identity threshold (line chart, kept
    from production).  At scale (hundreds-to-thousands of sequences)
    a pairwise distance heatmap would be unreadable and expensive to
    compute, so this plot stays as a curve.

mpnn_aa_composition.png
    Per-design AA frequency heatmap + ranked boxplot, with a
    **native AA composition reference line overlaid on the boxplot**
    (this iteration: was design-only before).  Shows whether MPNN
    is producing AA frequencies similar to what the original gap
    residues had.

mpnn_physicochem.png  (NEW)
    Per-design mean hydrophobicity (Kyte-Doolittle) and net charge
    of the design region, with native value as horizontal reference.
    Shows whether designs systematically differ from native in
    interface physicochemistry.

Usage
-----
    python test_mpnn_plots.py \\
        --metadata    tests/proteinmpnn/receptor_resurfacing_results/sequences/scored_metadata.csv \\
        --cluster-csv tests/proteinmpnn/receptor_resurfacing_results/sequences/mpnn_cluster_counts.csv \\
        --outdir      tests/proteinmpnn/receptor_resurfacing_results/plots_iter

If ``--metadata`` is omitted the script falls back to the canonical
test output path.  ``--cluster-csv`` is used by the diversity plot
(it reads from the precomputed ``mpnn_cluster_counts.csv``).
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MaxNLocator


AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"

# Kyte-Doolittle hydrophobicity index.
KD_HYDROPHOBICITY = {
    "A":  1.8, "C":  2.5, "D": -3.5, "E": -3.5, "F":  2.8,
    "G": -0.4, "H": -3.2, "I":  4.5, "K": -3.9, "L":  3.8,
    "M":  1.9, "N": -3.5, "P": -1.6, "Q": -3.5, "R": -4.5,
    "S": -0.8, "T": -0.7, "V":  4.2, "W": -0.9, "Y": -1.3,
}

# Net charge at neutral pH.  K, R = +1; D, E = -1; H ~ +0.1 (skip).
CHARGE = {"K": +1, "R": +1, "D": -1, "E": -1}


# ── Styling ────────────────────────────────────────────────────────────────
COLOUR_DESIGN = "#4C72B0"   # blue (designed sequences)
COLOUR_NATIVE = "#D62728"   # red  (native reference)
COLOUR_REGION = "#DD8452"   # orange (design-region score boxplots)


def _design_label(design_idx):
    return f"d{design_idx}"


# ═══════════════════════════════════════════════════════════════════════════
# Plot 1: Score distribution (flipped, ranked by design-region)
# ═══════════════════════════════════════════════════════════════════════════

def plot_score_distribution(rows, out_path):
    """Two stacked per-design boxplots.

    Iteration vs production:
    - Top panel is now design-region score (was global).
    - Bottom panel is global score (was design-region).
    - Ranking is by median design-region score (was global).

    Rationale: design-region score is what we're actually evaluating
    (the newly designed bit) — global is dominated by long fixed
    regions and barely moves between designs, so it's a sanity check.

    Falls back to global-ranked / global-only display if no
    design-region scores are present.
    """
    global_scores = defaultdict(list)
    region_scores = defaultdict(list)

    for row in rows:
        design = row.get("design", "?")
        try:
            gs = float(row.get("mpnn_score", "nan"))
            if np.isfinite(gs):
                global_scores[design].append(gs)
        except (ValueError, TypeError):
            pass
        try:
            rs = float(row.get("design_region_score", "nan"))
            if np.isfinite(rs):
                region_scores[design].append(rs)
        except (ValueError, TypeError):
            pass

    if not global_scores:
        _make_empty_plot("No MPNN scores available", out_path)
        return False

    has_region = bool(region_scores)

    # Rank by design-region median if available, else fall back to global.
    if has_region:
        sorted_designs = sorted(
            region_scores.keys(),
            key=lambda d: np.median(region_scores[d]) if region_scores[d] else float("inf"),
        )
        # Designs that have global but not region scores get appended at end.
        for d in sorted(global_scores.keys()):
            if d not in region_scores:
                sorted_designs.append(d)
        rank_label = "ranked by median design-region score"
    else:
        sorted_designs = sorted(
            global_scores.keys(),
            key=lambda d: np.median(global_scores[d]),
        )
        rank_label = "ranked by median global score"

    n = len(sorted_designs)
    labels = [_design_label(d) for d in sorted_designs]

    fig_w = max(6, min(20, n * 0.35 + 2))
    fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(fig_w, 8), sharex=True)

    # ── Top: design-region score (primary) ─────────────────────────────
    if has_region:
        data_region = [region_scores.get(d, []) for d in sorted_designs]
        nonempty = [(i, d) for i, d in enumerate(data_region) if d]
        if nonempty:
            positions = [x[0] + 1 for x in nonempty]
            vals = [x[1] for x in nonempty]
            bp_top = ax_top.boxplot(vals, positions=positions,
                                     patch_artist=True, showfliers=True,
                                     flierprops={"markersize": 2, "alpha": 0.4},
                                     widths=0.6)
            for patch in bp_top["boxes"]:
                patch.set_facecolor(COLOUR_REGION)
                patch.set_alpha(0.7)
        ax_top.set_ylabel("Design Region MPNN Score", fontsize=11)
    else:
        ax_top.text(0.5, 0.5,
                    "Design-region scores not available\n"
                    "(showing global only below)",
                    ha="center", va="center", transform=ax_top.transAxes,
                    fontsize=11, color="grey")
        ax_top.set_ylabel("Design Region MPNN Score", fontsize=11)
    ax_top.yaxis.grid(True, linestyle="-", linewidth=0.3, alpha=0.5,
                      color="grey")
    ax_top.set_axisbelow(True)

    # ── Bottom: global score (sanity check) ────────────────────────────
    data_global = [global_scores[d] for d in sorted_designs]
    bp_bot = ax_bot.boxplot(data_global, patch_artist=True, showfliers=True,
                             flierprops={"markersize": 2, "alpha": 0.4},
                             widths=0.6)
    for patch in bp_bot["boxes"]:
        patch.set_facecolor(COLOUR_DESIGN)
        patch.set_alpha(0.7)
    ax_bot.set_ylabel("Global MPNN Score", fontsize=11)
    ax_bot.yaxis.grid(True, linestyle="-", linewidth=0.3, alpha=0.5,
                      color="grey")
    ax_bot.set_axisbelow(True)

    ax_bot.set_xticks(range(1, n + 1))
    ax_bot.set_xticklabels(labels, fontsize=max(5, 8 - n // 20), rotation=0)
    ax_bot.set_xlabel(f"Design ({rank_label})", fontsize=10)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    return True


# ═══════════════════════════════════════════════════════════════════════════
# Plot 2: AA composition with native reference
# ═══════════════════════════════════════════════════════════════════════════

def _aa_composition_pct(seq_string):
    """Return a 20-element AA frequency array (%) for a single string,
    ignoring '|' region separators and any non-AA characters."""
    counts = np.zeros(len(AA_ORDER))
    total = 0
    for ch in seq_string:
        if ch in AA_ORDER:
            counts[AA_ORDER.index(ch)] += 1
            total += 1
    if total == 0:
        return counts
    return (counts / total) * 100.0


def plot_aa_composition(rows, out_path):
    """Per-design AA composition heatmap + boxplot, with the **native
    AA composition** (pooled across all rows) drawn as a red overlay
    on the boxplot.

    Iteration vs production:
    - Adds a native-composition reference (red dots and thin vertical
      lines) so you can tell at a glance whether MPNN is producing
      sequences with similar AA frequencies to the original gap
      residues — informative for receptor resurfacing where the
      native interface was under selection.
    """
    design_compositions = defaultdict(lambda: defaultdict(int))
    design_totals = defaultdict(int)
    # Native should be constant across rows (all rows share the same
    # input PDB and the same contig spec, so the gap residues that
    # the design replaces are identical across designs).  Defensively
    # collect the unique values — if we somehow get more than one
    # (e.g. concatenated runs), pick the longest and warn.
    native_unique = set()
    for row in rows:
        designed = row.get("designed_residues", "")
        design = row.get("design", "?")
        for aa in designed:
            if aa in AA_ORDER:
                design_compositions[design][aa] += 1
                design_totals[design] += 1
        nat = row.get("native_residues", "")
        if nat:
            native_unique.add(nat)

    if len(native_unique) > 1:
        print(f"WARNING: AA composition: found {len(native_unique)} distinct "
              f"native_residues strings across rows; using the longest.")
    native_seq = max(native_unique, key=len) if native_unique else ""

    if not design_totals:
        _make_empty_plot("No designed residues", out_path)
        return False

    designs = sorted(design_compositions.keys())
    n_designs = len(designs)
    n_aa = len(AA_ORDER)

    pct_matrix = np.zeros((n_designs, n_aa))
    for i, design in enumerate(designs):
        total = design_totals[design]
        if total > 0:
            for j, aa in enumerate(AA_ORDER):
                pct_matrix[i, j] = (design_compositions[design].get(aa, 0)
                                    / total * 100)

    # Native composition from the single input-structure gap residues.
    native_pct = _aa_composition_pct(native_seq)
    has_native = native_seq != ""

    # Sort AAs by median DESIGN frequency descending — same as before.
    median_freq = np.median(pct_matrix, axis=0)
    aa_sort_idx = np.argsort(-median_freq)
    sorted_aa = [AA_ORDER[i] for i in aa_sort_idx]
    sorted_pct = pct_matrix[:, aa_sort_idx]
    sorted_native = native_pct[aa_sort_idx]

    # ── Layout with dedicated colourbar column ────────────────────────
    heatmap_height = max(1.5, n_designs * 0.25)
    boxplot_height = 3.5
    total_height = heatmap_height + boxplot_height + 0.5
    fig_width = max(8, n_aa * 0.45 + 2)

    fig = plt.figure(figsize=(fig_width, total_height))
    gs = GridSpec(2, 2, width_ratios=[30, 1],
                  height_ratios=[heatmap_height, boxplot_height],
                  hspace=0.08, wspace=0.05, figure=fig)

    ax_heat = fig.add_subplot(gs[0, 0])
    ax_cbar = fig.add_subplot(gs[0, 1])
    ax_box = fig.add_subplot(gs[1, 0])
    ax_empty = fig.add_subplot(gs[1, 1])
    ax_empty.axis("off")

    # ── Heatmap ───────────────────────────────────────────────────────
    im = ax_heat.imshow(sorted_pct, aspect="auto", cmap="YlOrRd",
                        extent=[-0.5, n_aa - 0.5, n_designs - 0.5, -0.5])
    ax_heat.set_yticks(range(n_designs))
    ax_heat.set_yticklabels([_design_label(d) for d in designs], fontsize=7)
    ax_heat.set_ylabel("Design\n(design region)", fontsize=10)

    ax_heat.set_xticks(range(n_aa))
    ax_heat.set_xticklabels([])
    ax_heat.tick_params(axis="x", length=0)

    plt.colorbar(im, cax=ax_cbar, label="Frequency (%)")

    # ── Boxplot (design distribution) ─────────────────────────────────
    data_for_box = [sorted_pct[:, j] for j in range(n_aa)]
    bp = ax_box.boxplot(data_for_box, positions=range(n_aa),
                         patch_artist=True, showfliers=True,
                         flierprops={"markersize": 2, "alpha": 0.4},
                         widths=0.6)
    for patch in bp["boxes"]:
        patch.set_facecolor(COLOUR_DESIGN)
        patch.set_alpha(0.7)

    # ── Native overlay: red diamond markers on top of each AA's box ───
    if has_native:
        ax_box.scatter(range(n_aa), sorted_native, marker="D",
                       color=COLOUR_NATIVE, s=40, zorder=5,
                       edgecolors="black", linewidths=0.5,
                       label="Native")
        ax_box.legend(loc="upper right", fontsize=9, framealpha=0.95)

    ax_box.set_xlim(-0.5, n_aa - 0.5)
    ax_box.set_xticks(range(n_aa))
    ax_box.set_xticklabels(sorted_aa, fontsize=9, fontfamily="monospace")
    ax_box.set_xlabel("Amino Acid (ranked by median design frequency)",
                       fontsize=10)
    ax_box.set_ylabel("Frequency (%)\n(design region)", fontsize=10)

    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    return True



# ═══════════════════════════════════════════════════════════════════════════
# Plot 3: Sequence diversity — MMseqs2 cluster count vs identity threshold
# ═══════════════════════════════════════════════════════════════════════════
# Kept as the production line chart: at scale (hundreds-to-thousands of
# sequences) a pairwise distance heatmap becomes both unreadable and
# expensive to compute.  The line chart trades per-pair detail for
# scalability and is computed by the upstream ``mpnn_cluster_sequences.py``
# step which writes ``mpnn_cluster_counts.csv``.

def plot_sequence_diversity(cluster_csv, out_path):
    """Sequence diversity: cluster count vs identity threshold.

    Reads the precomputed ``mpnn_cluster_counts.csv`` (one row per
    threshold/cluster-count pair) emitted by ``mpnn_cluster_sequences.py``
    and plots the curve.  The curve's shape tells you how aggressively
    sequences merge as the identity threshold drops — a sharp drop
    means MPNN is producing distinct sequence families; a slow drop
    means many near-duplicates.
    """
    if not cluster_csv or not Path(cluster_csv).exists():
        _make_empty_plot("No clustering data available", out_path)
        return False

    thresholds = []
    n_clusters = []
    with open(cluster_csv) as f:
        for row in csv.DictReader(f):
            thresholds.append(float(row["threshold"]) * 100)
            n_clusters.append(int(row["n_clusters"]))

    if not thresholds:
        _make_empty_plot("No clustering data", out_path)
        return False

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(thresholds, n_clusters, "o-", color=COLOUR_DESIGN,
            linewidth=2, markersize=5)

    ax.set_xlabel("Sequence Identity Threshold (Design Region) (%)", fontsize=11)
    ax.set_ylabel("Number of Clusters", fontsize=11)
    ax.set_xlim(min(thresholds) - 2, max(thresholds) + 2)
    ax.set_ylim(bottom=0)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    return True


# ═══════════════════════════════════════════════════════════════════════════
# Plot 4 (NEW): Per-design hydrophobicity and net charge
# ═══════════════════════════════════════════════════════════════════════════

def _mean_hydrophobicity(seq):
    if not seq:
        return float("nan")
    vals = [KD_HYDROPHOBICITY[c] for c in seq if c in KD_HYDROPHOBICITY]
    return float(np.mean(vals)) if vals else float("nan")


def _net_charge(seq):
    if not seq:
        return float("nan")
    return float(sum(CHARGE.get(c, 0) for c in seq))


def plot_physicochem(rows, out_path):
    """Per-design hydrophobicity and net charge of the design region,
    with native value as a horizontal reference.

    Bypasses the variable-length-per-design problem by reducing each
    sequence to one number per axis.  Tells you whether designs are
    systematically more/less hydrophobic or more/less charged than
    native — important for surface-exposed interfaces.
    """
    per_design_hydro = defaultdict(list)
    per_design_charge = defaultdict(list)
    # See note in plot_aa_composition — native should be constant.
    native_unique = set()

    for row in rows:
        designed = row.get("designed_residues", "").replace("|", "")
        design = row.get("design", "?")
        if designed:
            per_design_hydro[design].append(_mean_hydrophobicity(designed))
            per_design_charge[design].append(_net_charge(designed))
        nat = row.get("native_residues", "").replace("|", "")
        if nat:
            native_unique.add(nat)

    if len(native_unique) > 1:
        print(f"WARNING: physicochem: found {len(native_unique)} distinct "
              f"native_residues strings across rows; using the longest.")
    native_seq = max(native_unique, key=len) if native_unique else ""

    if not per_design_hydro:
        _make_empty_plot("No designed residues for physicochem plot", out_path)
        return False

    designs = sorted(per_design_hydro.keys())
    n = len(designs)
    labels = [_design_label(d) for d in designs]

    # Native reference from the single input-structure gap residues.
    native_hydro = _mean_hydrophobicity(native_seq) if native_seq else None
    native_charge_per_len = (_net_charge(native_seq) / len(native_seq)
                              if native_seq else None)

    fig_w = max(6, min(20, n * 0.35 + 2))
    fig, (ax_h, ax_c) = plt.subplots(2, 1, figsize=(fig_w, 7), sharex=True)

    # ── Hydrophobicity ────────────────────────────────────────────────
    data_h = [per_design_hydro[d] for d in designs]
    bp_h = ax_h.boxplot(data_h, patch_artist=True, showfliers=True,
                         flierprops={"markersize": 2, "alpha": 0.4},
                         widths=0.6)
    for patch in bp_h["boxes"]:
        patch.set_facecolor(COLOUR_DESIGN)
        patch.set_alpha(0.7)

    if native_hydro is not None:
        ax_h.axhline(native_hydro, color=COLOUR_NATIVE, linestyle="--",
                      linewidth=1.2, zorder=4,
                      label=f"Native = {native_hydro:.2f}")
        ax_h.legend(loc="upper right", fontsize=9, framealpha=0.95)
    ax_h.set_ylabel("Mean Kyte–Doolittle\nhydrophobicity", fontsize=10)
    ax_h.yaxis.grid(True, linestyle="-", linewidth=0.3, alpha=0.5,
                     color="grey")
    ax_h.set_axisbelow(True)

    # ── Net charge per residue (so design-region length differences
    # don't dominate) ─────────────────────────────────────────────────
    data_c = []
    for d in designs:
        # Each entry in per_design_charge[d] is a raw net charge for
        # one designed sequence.  Normalise by the sequence's length
        # so designs with longer regions aren't penalised.
        per_seq_norm = []
        for raw_charge, seq in zip(per_design_charge[d],
                                    [r.get("designed_residues", "").replace("|", "")
                                     for r in rows if r.get("design") == d]):
            if seq:
                per_seq_norm.append(raw_charge / len(seq))
        data_c.append(per_seq_norm if per_seq_norm else [0])

    bp_c = ax_c.boxplot(data_c, patch_artist=True, showfliers=True,
                         flierprops={"markersize": 2, "alpha": 0.4},
                         widths=0.6)
    for patch in bp_c["boxes"]:
        patch.set_facecolor(COLOUR_DESIGN)
        patch.set_alpha(0.7)

    if native_charge_per_len is not None:
        ax_c.axhline(native_charge_per_len, color=COLOUR_NATIVE, linestyle="--",
                      linewidth=1.2, zorder=4,
                      label=f"Native = {native_charge_per_len:+.3f}")
        ax_c.legend(loc="upper right", fontsize=9, framealpha=0.95)
    ax_c.axhline(0, color="black", linewidth=0.5, zorder=2)
    ax_c.set_ylabel("Net charge per residue\n(at neutral pH)", fontsize=10)
    ax_c.yaxis.grid(True, linestyle="-", linewidth=0.3, alpha=0.5,
                     color="grey")
    ax_c.set_axisbelow(True)

    ax_c.set_xticks(range(1, n + 1))
    ax_c.set_xticklabels(labels, fontsize=max(5, 8 - n // 20), rotation=0)
    ax_c.set_xlabel("Design", fontsize=10)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    return True


# ═══════════════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════════════

def _make_empty_plot(message, out_path):
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.text(0.5, 0.5, message, ha="center", va="center",
            transform=ax.transAxes, fontsize=14, color="grey")
    ax.set_axis_off()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════
# Top-level orchestration
# ═══════════════════════════════════════════════════════════════════════════

def _resolve_metadata_path(supplied):
    if supplied:
        p = Path(supplied)
        if not p.exists():
            raise SystemExit(f"--metadata path does not exist: {p}")
        return p
    fallback = Path(
        "tests/proteinmpnn/receptor_resurfacing_results/"
        "sequences/scored_metadata.csv"
    )
    if fallback.exists():
        return fallback
    raise SystemExit(
        f"No --metadata passed and fallback {fallback} does not exist. "
        f"Run test_proteinmpnn first or pass --metadata explicitly."
    )


def _parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--metadata", default=None,
        help="QC metadata CSV (scored_metadata.csv) from a test_proteinmpnn run.")
    parser.add_argument("--cluster-csv", default="",
        help="Cluster counts CSV (mpnn_cluster_counts.csv) emitted by "
             "mpnn_cluster_sequences.py.  Used by the diversity plot.")
    parser.add_argument("--outdir", required=True,
        help="Directory to write plots into; created if missing.")
    return parser.parse_args()


def main():
    args = _parse_args()
    metadata_path = _resolve_metadata_path(args.metadata)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Reading metadata: {metadata_path}")
    rows = []
    with open(metadata_path) as f:
        for row in csv.DictReader(f):
            rows.append(row)
    print(f"  rows: {len(rows)}")

    if not rows:
        print("No rows in metadata — nothing to plot.")
        return 0

    plots = [
        ("Score distribution (flipped, ranked by design-region)",
         "mpnn_score_distribution.png",
         lambda p: plot_score_distribution(rows, p)),
        ("AA composition (with native reference)",
         "mpnn_aa_composition.png",
         lambda p: plot_aa_composition(rows, p)),
        ("Sequence diversity (cluster count vs identity threshold)",
         "mpnn_sequence_diversity.png",
         lambda p: plot_sequence_diversity(args.cluster_csv, p)),
        ("Physicochem (hydrophobicity + net charge)",
         "mpnn_physicochem.png",
         lambda p: plot_physicochem(rows, p)),
    ]

    for desc, fname, fn in plots:
        out = outdir / fname
        print(f"Building {desc}...")
        try:
            ok = fn(str(out))
            if ok:
                print(f"  → {out}")
            else:
                print(f"  (skipped — no data)")
        except Exception as e:
            print(f"  FAILED: {type(e).__name__}: {e}")

    print(f"\nDone.  Inspect: {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
