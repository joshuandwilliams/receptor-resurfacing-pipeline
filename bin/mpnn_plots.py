#!/usr/bin/env python3
"""
mpnn_plots.py
-------------
Diagnostic plots from ProteinMPNN sequence design.

    mpnn_score_distribution.png  - Two stacked subplots: design-region
                                   MPNN score (top, primary) and global
                                   MPNN score (bottom, sanity check).
                                   Ranked by median design-region score.
    mpnn_sequence_diversity.png  - MMseqs2 clusters vs identity threshold
    mpnn_aa_composition.png      - Per-design AA heatmap + ranked boxplot
                                   with the input PDB's native gap-residue
                                   AA composition overlaid as red diamonds.
    mpnn_physicochem.png         - Per-design mean hydrophobicity (Kyte-
                                   Doolittle) and net charge per residue,
                                   with native value as horizontal reference.
"""

import argparse
import csv
import os
import sys
from collections import defaultdict

import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    HAS_MPL = True
except ImportError:
    HAS_MPL = False


AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"

# Kyte-Doolittle hydrophobicity index.  Used by the physicochem plot.
KD_HYDROPHOBICITY = {
    "A":  1.8, "C":  2.5, "D": -3.5, "E": -3.5, "F":  2.8,
    "G": -0.4, "H": -3.2, "I":  4.5, "K": -3.9, "L":  3.8,
    "M":  1.9, "N": -3.5, "P": -1.6, "Q": -3.5, "R": -4.5,
    "S": -0.8, "T": -0.7, "V":  4.2, "W": -0.9, "Y": -1.3,
}

# Net charge at neutral pH.  H is ~+0.1 — skipped as commonly done.
CHARGE = {"K": +1, "R": +1, "D": -1, "E": -1}

# ── Styling ──────────────────────────────────────────────────────────────
COLOUR_DESIGN = "#4C72B0"   # blue (designed sequences)
COLOUR_NATIVE = "#D62728"   # red  (native reference)
COLOUR_REGION = "#DD8452"   # orange (design-region score boxplots)

ALL_PLOT_FILES = [
    "mpnn_score_distribution.png",
    "mpnn_sequence_diversity.png",
    "mpnn_aa_composition.png",
    "mpnn_physicochem.png",
]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True,
                        help="QC metadata CSV (scored_metadata.csv)")
    parser.add_argument("--cluster-csv", default="",
                        help="Cluster counts CSV from mpnn_cluster_sequences.py")
    parser.add_argument("--receptor-seq", default="",
                        help="(unused, kept for CLI compatibility)")
    parser.add_argument("--contigs", default="",
                        help="(unused, kept for CLI compatibility)")
    return parser.parse_args()


def make_empty_plot(message, path):
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.text(0.5, 0.5, message, ha="center", va="center",
            transform=ax.transAxes, fontsize=14, color="grey")
    ax.set_axis_off()
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()


def save_fallback_plots(message):
    for name in ALL_PLOT_FILES:
        make_empty_plot(message, name)


def _design_label(design_idx):
    return f"d{design_idx}"


# ==============================================================================
# Plot 1: MPNN score — two stacked subplots (global + design-region)
# ==============================================================================

def plot_score_distribution(rows):
    """Two stacked per-design boxplots sharing the same x-axis order.

    Top:    design-region MPNN score (primary — the newly designed bit
            that we're actually evaluating).
    Bottom: global MPNN score (sanity check that global also looks
            reasonable; dominated by the long fixed regions, so barely
            moves between designs).

    Both are ranked left-to-right by median design-region score (best
    first).  Falls back to global-only / global-ranked display if no
    design-region scores are present in the metadata.
    """
    out_path = "mpnn_score_distribution.png"

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
        make_empty_plot("No MPNN scores available", out_path)
        return

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

    # ── Top: design-region score (primary) ────────────────────────────
    if has_region:
        data_region = [region_scores.get(d, []) for d in sorted_designs]
        nonempty = [(i, d) for i, d in enumerate(data_region) if d]
        if nonempty:
            positions = [x[0] + 1 for x in nonempty]   # boxplot is 1-indexed
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

    # ── Bottom: global score (sanity check) ───────────────────────────
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
    print(f"Saved {out_path}")


# ==============================================================================
# Plot 2: Sequence diversity — cluster count vs identity threshold
# ==============================================================================

def plot_sequence_diversity(cluster_csv):
    out_path = "mpnn_sequence_diversity.png"

    if not cluster_csv or not os.path.exists(cluster_csv):
        make_empty_plot("No clustering data available", out_path)
        return

    thresholds = []
    n_clusters = []

    with open(cluster_csv) as f:
        for row in csv.DictReader(f):
            thresholds.append(float(row["threshold"]) * 100)
            n_clusters.append(int(row["n_clusters"]))

    if not thresholds:
        make_empty_plot("No clustering data", out_path)
        return

    fig, ax = plt.subplots(figsize=(8, 5))

    ax.plot(thresholds, n_clusters, "o-", color="#4C72B0",
            linewidth=2, markersize=5)

    ax.set_xlabel("Sequence Identity Threshold (Design Region) (%)", fontsize=11)
    ax.set_ylabel("Number of Clusters", fontsize=11)

    ax.set_xlim(min(thresholds) - 2, max(thresholds) + 2)
    ax.set_ylim(bottom=0)
    ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved {out_path}")


# ==============================================================================
# Plot 3: AA composition — heatmap + boxplot, fully aligned
# ==============================================================================

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


def plot_aa_composition(rows):
    """Per-design AA composition heatmap + boxplot, with the native
    AA composition (gap residues from the input PDB) drawn as red
    diamonds on top of the boxplot.

    The native string spans all design regions, joined with '|'
    separators by ``pipeline_correct_sequences.py``.  Every row in
    one run carries the same native string (gap positions are
    determined by the contig spec, which doesn't change between
    designs), so we collect the unique values defensively and pick
    the longest if more than one is found (a sign of concatenated
    runs — print a warning in that case).
    """
    out_path = "mpnn_aa_composition.png"

    design_compositions = defaultdict(lambda: defaultdict(int))
    design_totals = defaultdict(int)
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
              f"native_residues strings across rows; using the longest.",
              file=sys.stderr)
    native_seq = max(native_unique, key=len) if native_unique else ""

    if not design_totals:
        make_empty_plot("No designed residues", out_path)
        return

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

    # Native composition.
    native_pct = _aa_composition_pct(native_seq)
    has_native = native_seq != ""

    # Sort AAs by median DESIGN frequency descending — same ordering
    # as before, so heatmap and boxplot stay aligned.
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

    # 2 rows, 2 columns: [heatmap, cbar] / [boxplot, empty]
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

    # ── Native overlay: red diamond markers per AA ────────────────────
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
    print(f"Saved {out_path}")


# ==============================================================================
# Plot 4: Hydrophobicity & net charge per design
# ==============================================================================

def _mean_hydrophobicity(seq):
    if not seq:
        return float("nan")
    vals = [KD_HYDROPHOBICITY[c] for c in seq if c in KD_HYDROPHOBICITY]
    return float(np.mean(vals)) if vals else float("nan")


def _net_charge(seq):
    if not seq:
        return 0.0
    return float(sum(CHARGE.get(c, 0) for c in seq))


def plot_physicochem(rows):
    """Per-design mean hydrophobicity (Kyte–Doolittle) and net charge
    per residue across the design region(s), with the input PDB's
    native gap-residue value drawn as a horizontal red dashed
    reference line on each panel.

    Sidesteps the variable-length-per-design problem of positional
    analysis by reducing each sequence to one number per axis.  Tells
    you whether designs systematically differ from native in interface
    physicochemistry — important for surface-exposed interfaces.

    Native values are computed from the unique (defensively collected)
    ``native_residues`` string in the metadata; the longest is picked
    if more than one is found, with a warning.
    """
    out_path = "mpnn_physicochem.png"

    per_design_hydro = defaultdict(list)
    per_design_charge_per_res = defaultdict(list)
    native_unique = set()

    for row in rows:
        designed = row.get("designed_residues", "").replace("|", "")
        design = row.get("design", "?")
        if designed:
            per_design_hydro[design].append(_mean_hydrophobicity(designed))
            # Normalise charge by sequence length so designs with
            # longer regions aren't artificially penalised.
            per_design_charge_per_res[design].append(
                _net_charge(designed) / len(designed)
            )
        nat = row.get("native_residues", "").replace("|", "")
        if nat:
            native_unique.add(nat)

    if not per_design_hydro:
        make_empty_plot("No designed residues for physicochem plot", out_path)
        return

    if len(native_unique) > 1:
        print(f"WARNING: physicochem: found {len(native_unique)} distinct "
              f"native_residues strings across rows; using the longest.",
              file=sys.stderr)
    native_seq = max(native_unique, key=len) if native_unique else ""
    native_hydro = _mean_hydrophobicity(native_seq) if native_seq else None
    native_charge_per_len = (_net_charge(native_seq) / len(native_seq)
                             if native_seq else None)

    designs = sorted(per_design_hydro.keys())
    n = len(designs)
    labels = [_design_label(d) for d in designs]

    fig_w = max(6, min(20, n * 0.35 + 2))
    fig, (ax_h, ax_c) = plt.subplots(2, 1, figsize=(fig_w, 7), sharex=True)

    # ── Hydrophobicity ───────────────────────────────────────────────
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

    # ── Net charge per residue ───────────────────────────────────────
    data_c = [per_design_charge_per_res[d] for d in designs]
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
    print(f"Saved {out_path}")


# ==============================================================================
# Main
# ==============================================================================

def main():
    args = parse_args()

    if not HAS_MPL:
        for name in ALL_PLOT_FILES:
            with open(name, "w") as f:
                f.write("matplotlib not available")
        sys.exit(0)

    if not os.path.exists(args.metadata) or os.path.getsize(args.metadata) == 0:
        save_fallback_plots("No metadata")
        sys.exit(0)

    rows = []
    with open(args.metadata) as f:
        for row in csv.DictReader(f):
            rows.append(row)

    if not rows:
        save_fallback_plots("No sequences in metadata")
        sys.exit(0)

    plot_score_distribution(rows)
    plot_sequence_diversity(args.cluster_csv)
    plot_aa_composition(rows)
    plot_physicochem(rows)


if __name__ == "__main__":
    main()