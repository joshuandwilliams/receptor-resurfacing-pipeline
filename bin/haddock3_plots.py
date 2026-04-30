#!/usr/bin/env python3
"""
haddock3_plots.py
-----------------
Diagnostic plots from a HADDOCK3 docking run.

Produces three panels:
    haddock_score_vs_bsa.png       - Score vs BSA scatter
    haddock_cluster_sizes.png      - Cluster size bar chart by mean score
    haddock_interface_heatmap.png  - Per-cluster receptor & effector contact heatmaps
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
import numpy as np

from haddock_utils import (
    parse_capri_tsv, parse_clustfcc_tsv, get_numeric_col,
    model_stem, extract_heavy_atoms, contacted_residues,
    collect_pdb_index, cluster_mean_score,
)

# ── Styling ──────────────────────────────────────────────────────────────────
COLOUR_ALL      = "#4C72B0"
COLOUR_BEST     = "#DD4444"
COLOUR_WARN     = "#FFCCCC"
COLOUR_DENOVO   = "#FF8C00"
BSA_WARN_CUTOFF = 1000.0
CONTACT_CUTOFF  = 8.0
MIN_CLUSTER_SIZE = 4  # Default; overridden by --min-cluster-size CLI flag.
                      # Single source of truth lives in nextflow.config as
                      # params.haddock_min_cluster_size, plumbed via haddock.nf.


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--capri-scores", required=True, help="Path to capri_scores.tsv")
    parser.add_argument("--cluster-summary", default=None, help="Path to clustfcc.tsv")
    parser.add_argument("--run-dir", default=None,
                        help="HADDOCK3 run directory (auto-discovers PDBs)")
    parser.add_argument("--complex-dir", default=None,
                        help="Directory with docked PDBs (fallback if no --run-dir)")
    parser.add_argument("--receptor-chain", default="A")
    parser.add_argument("--effector-chain", default="B")
    parser.add_argument("--contact-cutoff", type=float, default=CONTACT_CUTOFF,
                        help=f"Heavy-atom distance cutoff in Å (default: {CONTACT_CUTOFF})")
    parser.add_argument("--contigs", default=None,
                        help="RFDiffusion contig string for design-region shading")
    parser.add_argument("--min-cluster-size", type=int, default=MIN_CLUSTER_SIZE,
                        help="Minimum cluster size to qualify (default: %(default)s). "
                             "Should match `min_population` in haddock.nf clustfcc block.")
    parser.add_argument("--effector-active-residues", default="",
                        help="Comma-separated effector active residues for AIR annotation bar")
    return parser.parse_args()


# ═══════════════════════════════════════════════════════════════════════════════
# Utility
# ═══════════════════════════════════════════════════════════════════════════════

def make_empty_plot(message, path):
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.text(0.5, 0.5, message, ha="center", va="center",
            transform=ax.transAxes, fontsize=12, color="grey")
    ax.set_axis_off()
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()


def save_fallback_plots(message):
    for name in ["haddock_score_vs_bsa.png", "haddock_cluster_sizes.png",
                  "haddock_interface_heatmap.png"]:
        make_empty_plot(message, name)


# ═══════════════════════════════════════════════════════════════════════════════
# Plot 1: Score vs BSA
# ═══════════════════════════════════════════════════════════════════════════════

def plot_score_vs_bsa(scores, bsas, best_score, best_bsa):
    """HADDOCK score vs BSA scatter. Ideal region is bottom-right."""
    fig, ax = plt.subplots(figsize=(7, 5))

    ax.scatter(bsas, scores, color=COLOUR_ALL, alpha=0.6, s=30,
               edgecolors="white", linewidths=0.5, label="All models", zorder=3)
    ax.scatter(best_bsa, best_score, color=COLOUR_BEST, s=120, zorder=5,
               edgecolors="black", linewidths=0.8,
               label=f"Best model (BSA={best_bsa:.0f} Å²)", marker="*")

    saved_xlim = ax.get_xlim()
    ax.axvspan(saved_xlim[0], BSA_WARN_CUTOFF, color=COLOUR_WARN, alpha=0.5,
               label=f"BSA < {BSA_WARN_CUTOFF:.0f} Å² (weak interface)", zorder=1)
    ax.set_xlim(saved_xlim)

    ax.set_xlabel("Buried Surface Area (Å²)", fontsize=12)
    ax.set_ylabel("HADDOCK Score", fontsize=12)

    y_min, y_max = ax.get_ylim()
    ax.set_ylim(top=y_max + (y_max - y_min) * 0.25)
    ax.legend(fontsize=10, loc="upper right")

    plt.tight_layout()
    plt.savefig("haddock_score_vs_bsa.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("Saved haddock_score_vs_bsa.png")


# ═══════════════════════════════════════════════════════════════════════════════
# Plot 2: Cluster sizes
# ═══════════════════════════════════════════════════════════════════════════════

def plot_cluster_sizes(clusters, best_cluster, best_model_cluster):
    """Cluster size bar chart coloured by mean HADDOCK score."""
    if not clusters:
        make_empty_plot("No cluster data available", "haddock_cluster_sizes.png")
        print("Saved haddock_cluster_sizes.png (no cluster data)")
        return

    means = {cid: (np.mean([s for (_, s) in m if s is not None])
                    if any(s is not None for (_, s) in m) else 0.0)
             for cid, m in clusters.items()}

    sorted_clusters = sorted(clusters.items(), key=lambda x: (x[0] == "-", -len(x[1])))
    cluster_ids   = [c[0] for c in sorted_clusters]
    cluster_sizes = [len(c[1]) for c in sorted_clusters]
    cluster_means = [means[c[0]] for c in sorted_clusters]

    from matplotlib.colors import TwoSlopeNorm
    vmin = min(cluster_means)
    vmax = max(max(cluster_means), abs(vmin))
    if vmin >= 0: vmin = -1.0
    if vmax <= 0: vmax = 1.0
    norm = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
    cmap = plt.cm.RdBu
    colours = [cmap(norm(m)) for m in cluster_means]

    fig, ax = plt.subplots(figsize=(max(7, len(cluster_ids) * 0.8 + 2), 5))
    bars = ax.bar(range(len(cluster_ids)), cluster_sizes, color=colours,
                  edgecolor="white", linewidth=0.8)

    for i, (bar, cid) in enumerate(zip(bars, cluster_ids)):
        if cid == best_model_cluster:
            ax.text(i, bar.get_height() + 0.3, "★", ha="center", va="bottom",
                    color="black", fontsize=14)

    tick_labels = ["Unclustered" if cid == "-" else f"Cluster {cid}" for cid in cluster_ids]
    ax.set_xticks(range(len(cluster_ids)))
    ax.set_xticklabels(tick_labels, rotation=45, ha="right", fontsize=9)
    ax.set_xlabel("FCC Cluster", fontsize=12)
    ax.set_ylabel("Number of Models", fontsize=12)
    ax.set_ylim(top=ax.get_ylim()[1] * 1.25)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_label("Mean HADDOCK Score", fontsize=10)

    ax.text(0.97, 0.95, "★ Pipeline-selected cluster",
            transform=ax.transAxes, ha="right", va="top", fontsize=10)

    plt.tight_layout()
    plt.savefig("haddock_cluster_sizes.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("Saved haddock_cluster_sizes.png")


# ═══════════════════════════════════════════════════════════════════════════════
# Plot 3: Interface contact heatmap
# ═══════════════════════════════════════════════════════════════════════════════

def _compute_contact_frequencies(clusters, pdb_index, rec_chain, eff_chain, cutoff):
    """
    For each named cluster, compute per-residue contact frequency for BOTH
    chains in a single pass per PDB (reads each file once for both chains).

    Returns for each chain:
        cluster_freqs : {cid: {resnum: frequency}}
        all_resnums   : sorted list of all residue numbers seen
        models_found  : {cid: int}
        models_total  : {cid: int}
    """
    rec_freqs, eff_freqs = {}, {}
    rec_all, eff_all = set(), set()
    rec_found, eff_found = {}, {}
    rec_total, eff_total = {}, {}

    for cid, members in clusters.items():
        if cid == "-":
            continue

        rec_total[cid] = eff_total[cid] = len(members)
        rec_counts, eff_counts = {}, {}
        found = 0

        for model_name, _score in members:
            stem = model_stem(model_name)
            fpath = pdb_index.get(stem)
            if fpath is None:
                continue

            # Single-pass: extract both chains at once
            both = extract_heavy_atoms(fpath, (rec_chain, eff_chain))
            r_atoms, e_atoms = both[rec_chain], both[eff_chain]
            if not r_atoms or not e_atoms:
                continue
            found += 1

            for r in contacted_residues(r_atoms, e_atoms, cutoff):
                rec_counts[r] = rec_counts.get(r, 0) + 1
            for r in contacted_residues(e_atoms, r_atoms, cutoff):
                eff_counts[r] = eff_counts.get(r, 0) + 1

            rec_all.update(r_atoms.keys())
            eff_all.update(e_atoms.keys())

        rec_found[cid] = eff_found[cid] = found
        rec_freqs[cid] = {r: c / found for r, c in rec_counts.items()} if found else {}
        eff_freqs[cid] = {r: c / found for r, c in eff_counts.items()} if found else {}

    return (rec_freqs, sorted(rec_all), rec_found, rec_total,
            eff_freqs, sorted(eff_all), eff_found, eff_total)


def _parse_fixed_residues(contigs, rec_chain):
    """
    Parse the contig string to extract the set of **fixed** receptor PDB
    residue numbers.

    Fixed segments are explicit ranges like A1-32, A46-72.  Everything
    else in the receptor (de novo segments) will be replaced by RFDiffusion,
    but those residues still exist in the pre-RFDiffusion PDB.

    Returns a set of PDB residue numbers that are fixed, or None if no
    contig is provided.
    """
    if not contigs:
        return None

    blocks = contigs.replace(",", " ").split()
    rec_block = None
    for block in blocks:
        if any(seg.strip() and seg.strip()[0].upper() == rec_chain.upper()
               for seg in block.split("/") if seg.strip() and seg.strip()[0].isalpha()):
            rec_block = block
            break
    if rec_block is None:
        return None

    fixed = set()
    for seg in rec_block.split("/"):
        seg = seg.strip()
        if not seg:
            continue
        if seg[0].isalpha() and seg[0].upper() == rec_chain.upper():
            rest = seg[1:]
            if "-" in rest:
                parts = rest.split("-")
                fixed.update(range(int(parts[0]), int(parts[1]) + 1))
            else:
                fixed.add(int(rest))
    return fixed


def _design_ranges_from_fixed(fixed_residues, all_pdb_residues):
    """
    Given the set of fixed residue numbers and the full set of receptor
    residue numbers present in the PDB, return contiguous (start, end)
    ranges for all non-fixed (design) residues.
    """
    if not fixed_residues or not all_pdb_residues:
        return []

    design_residues = sorted(set(all_pdb_residues) - fixed_residues)
    if not design_residues:
        return []

    # Group into contiguous runs
    ranges = []
    run_start = design_residues[0]
    prev = run_start
    for r in design_residues[1:]:
        if r != prev + 1:
            ranges.append((run_start, prev))
            run_start = r
        prev = r
    ranges.append((run_start, prev))
    return ranges


def _build_heatmap_matrix(sorted_cids, cluster_freqs, all_resnums, res_range,
                          best_model_contacts=None):
    """Build a 2-D matrix (rows × residue columns) for an interface heatmap."""
    all_set = set(all_resnums)
    n_res = len(res_range)
    has_best = best_model_contacts is not None
    n_rows = len(sorted_cids) + (1 if has_best else 0)

    matrix = np.full((n_rows, n_res), np.nan)
    row_offset = 0

    if has_best:
        for col_idx, resnum in enumerate(res_range):
            if resnum in best_model_contacts:
                matrix[0, col_idx] = 1.0
            elif resnum in all_set:
                matrix[0, col_idx] = 0.0
        row_offset = 1

    for row_idx, cid in enumerate(sorted_cids):
        freqs = cluster_freqs.get(cid, {})
        for col_idx, resnum in enumerate(res_range):
            if resnum in freqs:
                matrix[row_idx + row_offset, col_idx] = freqs[resnum]
            elif resnum in all_set:
                matrix[row_idx + row_offset, col_idx] = 0.0

    return matrix, row_offset


def _render_heatmap_panel(ax, matrix, res_range, n_rows, cmap, has_best_row):
    """Render a contact-frequency heatmap on the given axes."""
    n_res = len(res_range)
    im = ax.imshow(
        matrix, aspect="auto", cmap=cmap, vmin=0.0, vmax=1.0,
        interpolation="nearest", origin="upper",
        extent=[-0.5, n_res - 0.5, n_rows - 0.5, -0.5],
    )
    if has_best_row:
        ax.axhline(0.5, color="#888888", linewidth=1.0, linestyle="-", zorder=3)
    return im


def _set_heatmap_yticks(ax, sorted_cids, mean_scores, models_found, models_total,
                        clusters, row_offset, has_best_row, best_model_name=None):
    """Set Y-axis tick labels for a heatmap panel."""
    labels, positions = [], []

    if has_best_row:
        labels.append("★ Best model")
        positions.append(0)

    for row_idx, cid in enumerate(sorted_cids):
        size = models_total.get(cid, len(clusters.get(cid, [])))
        mean = mean_scores.get(cid)
        found = models_found.get(cid, 0)
        label = f"Cluster {cid} (N={size}, H={mean:.1f})" if mean is not None else f"Cluster {cid} (N={size})"
        if found < size:
            label += f" [{found} PDBs]"
        labels.append(label)
        positions.append(row_idx + row_offset)

    ax.set_yticks(positions)
    ax.set_yticklabels(labels, fontsize=9)


def _set_heatmap_xticks(ax, res_range, xlabel):
    """Set X-axis ticks at regular intervals, always showing the last residue."""
    n_res = len(res_range)
    tick_step = max(10, round(n_res / 20 / 10) * 10)
    positions = list(range(0, n_res, tick_step))

    last = n_res - 1
    if last not in positions:
        min_gap = max(tick_step // 2, 5)
        if positions and (last - positions[-1]) < min_gap:
            positions[-1] = last
        else:
            positions.append(last)

    ax.set_xticks(positions)
    ax.set_xticklabels([str(res_range[i]) for i in positions], fontsize=8, rotation=45, ha="right")
    ax.set_xlabel(xlabel, fontsize=11)


def plot_interface_heatmap(clusters, pdb_index, rec_chain, eff_chain, cutoff,
                           fixed_residues=None, best_model_name=None,
                           effector_active_residues=None):
    """Combined receptor + effector interface contact frequency heatmap."""
    out_path = "haddock_interface_heatmap.png"

    if not clusters:
        make_empty_plot("No cluster data available", out_path)
        print("Saved haddock_interface_heatmap.png (no cluster data)")
        return

    if not pdb_index:
        make_empty_plot("No complex PDBs found.\nProvide --complex-dir to enable\nthis plot.",
                        out_path)
        print("Saved haddock_interface_heatmap.png (no PDBs found)")
        return

    # ── Compute per-cluster contact frequencies (both chains, single pass) ──
    (rec_freqs, rec_all, rec_found, rec_total,
     eff_freqs, eff_all, eff_found, eff_total) = \
        _compute_contact_frequencies(clusters, pdb_index, rec_chain, eff_chain, cutoff)

    named = {cid: v for cid, v in clusters.items() if cid != "-"}
    if not named:
        make_empty_plot("No named clusters found", out_path)
        print("Saved haddock_interface_heatmap.png (no named clusters)")
        return
    if not rec_all:
        make_empty_plot("No receptor residues found in PDB files.\n"
                        "Check --receptor-chain and --complex-dir.", out_path)
        print("Saved haddock_interface_heatmap.png (no receptor residues)")
        return

    # ── Sort clusters by mean score ──────────────────────────────────────
    mean_scores = {}
    for cid, members in clusters.items():
        scored = [s for (_, s) in members if s is not None]
        mean_scores[cid] = np.mean(scored) if scored else None

    sorted_cids = sorted(named.keys(),
                         key=lambda c: mean_scores[c] if mean_scores[c] is not None else float("inf"))

    # ── Best-model contacts (single-pass PDB read) ──────────────────────
    best_rec_contacts = best_eff_contacts = None
    if best_model_name:
        stem = model_stem(best_model_name)
        fpath = pdb_index.get(stem)
        if fpath:
            both = extract_heavy_atoms(fpath, (rec_chain, eff_chain))
            r_atoms, e_atoms = both[rec_chain], both[eff_chain]
            if r_atoms and e_atoms:
                best_rec_contacts = contacted_residues(r_atoms, e_atoms, cutoff)
                best_eff_contacts = contacted_residues(e_atoms, r_atoms, cutoff)

    has_best_row = best_rec_contacts is not None

    # ── Receptor matrix ──────────────────────────────────────────────────
    rec_res_range = list(range(rec_all[0], rec_all[-1] + 1))
    rec_matrix, rec_row_offset = _build_heatmap_matrix(
        sorted_cids, rec_freqs, rec_all, rec_res_range, best_rec_contacts)
    n_rec_res = len(rec_res_range)
    n_rec_rows = rec_matrix.shape[0]

    # ── Effector matrix ──────────────────────────────────────────────────
    has_eff = bool(eff_all)
    if has_eff:
        eff_res_range = list(range(eff_all[0], eff_all[-1] + 1))
        eff_matrix, eff_row_offset = _build_heatmap_matrix(
            sorted_cids, eff_freqs, eff_all, eff_res_range, best_eff_contacts)
        n_eff_res = len(eff_res_range)
        n_eff_rows = eff_matrix.shape[0]
    else:
        n_eff_res = n_eff_rows = 0

    # ── Figure layout ────────────────────────────────────────────────────
    fig_w = min(30, max(10, max(n_rec_res, n_eff_res if has_eff else 0) * 0.04 + 3))
    rec_res_to_col = {r: i for i, r in enumerate(rec_res_range)}

    # Design regions = PDB residues that exist but are NOT fixed
    denovo_ranges = _design_ranges_from_fixed(fixed_residues, set(rec_res_range)) \
                    if fixed_residues is not None else []
    has_denovo = bool(denovo_ranges)
    if has_denovo:
        print(f"Design regions (PDB numbering): {denovo_ranges}")

    has_eff_active = bool(effector_active_residues) and has_eff

    BAR_H_IN = 0.30
    GAP_H_IN = 1.20
    rec_heat_h = max(2.5, n_rec_rows * 0.55 + 0.6)
    eff_heat_h = max(2.5, n_eff_rows * 0.55 + 0.6) if has_eff else 0.0

    # Build gridspec rows: [rec_bar, rec_heat, gap, (eff_bar), eff_heat]
    if has_eff:
        if has_eff_active:
            h_ratios = [BAR_H_IN, rec_heat_h, GAP_H_IN, BAR_H_IN, eff_heat_h]
            n_gs_rows = 5
        else:
            h_ratios = [BAR_H_IN, rec_heat_h, GAP_H_IN, eff_heat_h]
            n_gs_rows = 4
        total_h = sum(h_ratios) + 1.6
        fig = plt.figure(figsize=(fig_w, total_h))
        gs = fig.add_gridspec(n_gs_rows, 1, height_ratios=h_ratios, hspace=0.0)
        ax_bar = fig.add_subplot(gs[0])
        ax_rec = fig.add_subplot(gs[1])
        if has_eff_active:
            ax_eff_bar = fig.add_subplot(gs[3])
            ax_eff = fig.add_subplot(gs[4])
            all_axes = [ax_bar, ax_rec, ax_eff_bar, ax_eff]
        else:
            ax_eff_bar = None
            ax_eff = fig.add_subplot(gs[3])
            all_axes = [ax_bar, ax_rec, ax_eff]
    else:
        fig = plt.figure(figsize=(fig_w, BAR_H_IN + rec_heat_h + 1.6))
        gs = fig.add_gridspec(2, 1, height_ratios=[BAR_H_IN, rec_heat_h], hspace=0.06)
        ax_bar = fig.add_subplot(gs[0])
        ax_rec = fig.add_subplot(gs[1])
        ax_eff = None
        ax_eff_bar = None
        all_axes = [ax_bar, ax_rec]

    # ── Design-region annotation bar ─────────────────────────────────────
    bar_data = np.zeros((1, n_rec_res))
    if has_denovo:
        for dn_start, dn_end in denovo_ranges:
            for r in range(dn_start, dn_end + 1):
                if r in rec_res_to_col:
                    bar_data[0, rec_res_to_col[r]] = 1.0

    from matplotlib.colors import ListedColormap
    bar_cmap = ListedColormap(["#DDDDDD", COLOUR_DENOVO])
    ax_bar.imshow(bar_data, aspect="auto", cmap=bar_cmap, vmin=0, vmax=1,
                  interpolation="nearest", extent=[-0.5, n_rec_res - 0.5, 0, 1])
    ax_bar.set_yticks([])
    ax_bar.set_ylabel("", visible=False)
    for spine in ["top", "left", "right"]:
        ax_bar.spines[spine].set_visible(False)
    ax_bar.spines["bottom"].set_linewidth(0.5)
    ax_bar.spines["bottom"].set_color("#AAAAAA")
    ax_bar.tick_params(bottom=False, labelbottom=False)
    ax_bar.set_xlim(-0.5, n_rec_res - 0.5)

    if has_denovo:
        n_regions = len(denovo_ranges)
        ax_bar.text(-0.01, 0.5,
                    "Design region" if n_regions == 1 else "Design regions",
                    ha="right", va="center", transform=ax_bar.transAxes, fontsize=9)

    # ── Shared colourmap ─────────────────────────────────────────────────
    cmap = plt.cm.Blues.copy()
    cmap.set_bad(color="#EEEEEE")

    # ── Receptor heatmap ─────────────────────────────────────────────────
    im_rec = _render_heatmap_panel(ax_rec, rec_matrix, rec_res_range, n_rec_rows, cmap, has_best_row)

    if has_denovo:
        for dn_start, dn_end in denovo_ranges:
            cols = [rec_res_to_col[r] for r in range(dn_start, dn_end + 1) if r in rec_res_to_col]
            if cols:
                for xv in (min(cols) - 0.5, max(cols) + 0.5):
                    ax_rec.axvline(xv, color=COLOUR_DENOVO, linewidth=1.2,
                                   linestyle="--", alpha=0.95, zorder=2)

    _set_heatmap_yticks(ax_rec, sorted_cids, mean_scores, rec_found, rec_total,
                        clusters, rec_row_offset, has_best_row, best_model_name)
    ax_rec.set_ylabel("FCC Cluster", fontsize=11)
    _set_heatmap_xticks(ax_rec, rec_res_range, "Receptor residue")
    ax_rec.set_xlim(-0.5, n_rec_res - 0.5)

    # ── Effector active residues annotation bar ──────────────────────────
    if ax_eff_bar is not None and has_eff:
        eff_res_to_col = {r: i for i, r in enumerate(eff_res_range)}
        eff_bar_data = np.zeros((1, n_eff_res))
        for r in effector_active_residues:
            if r in eff_res_to_col:
                eff_bar_data[0, eff_res_to_col[r]] = 1.0

        COLOUR_EFF_ACTIVE = "#8B5CF6"  # purple to distinguish from receptor design bar
        eff_bar_cmap = ListedColormap(["#DDDDDD", COLOUR_EFF_ACTIVE])
        ax_eff_bar.imshow(eff_bar_data, aspect="auto", cmap=eff_bar_cmap, vmin=0, vmax=1,
                          interpolation="nearest",
                          extent=[-0.5, n_eff_res - 0.5, 0, 1])
        ax_eff_bar.set_yticks([])
        ax_eff_bar.set_ylabel("", visible=False)
        for spine in ["top", "left", "right"]:
            ax_eff_bar.spines[spine].set_visible(False)
        ax_eff_bar.spines["bottom"].set_linewidth(0.5)
        ax_eff_bar.spines["bottom"].set_color("#AAAAAA")
        ax_eff_bar.tick_params(bottom=False, labelbottom=False)
        ax_eff_bar.set_xlim(-0.5, n_eff_res - 0.5)
        ax_eff_bar.text(-0.01, 0.5, "Active residues",
                        ha="right", va="center", transform=ax_eff_bar.transAxes, fontsize=9)

    # ── Effector heatmap ─────────────────────────────────────────────────
    if has_eff and ax_eff is not None:
        _render_heatmap_panel(ax_eff, eff_matrix, eff_res_range, n_eff_rows, cmap, has_best_row)
        _set_heatmap_yticks(ax_eff, sorted_cids, mean_scores, eff_found, eff_total,
                            clusters, eff_row_offset, has_best_row, best_model_name)
        ax_eff.set_ylabel("FCC Cluster", fontsize=11)
        _set_heatmap_xticks(ax_eff, eff_res_range, "Effector residue")
        ax_eff.set_xlim(-0.5, n_eff_res - 0.5)

    # ── Shared colourbar ─────────────────────────────────────────────────
    cbar = fig.colorbar(im_rec, ax=all_axes, pad=0.01, fraction=0.02, aspect=30)
    cbar.set_label("Contact frequency", fontsize=10)
    cbar.set_ticks([0, 0.25, 0.5, 0.75, 1.0])
    cbar.set_ticklabels(["0", "0.25", "0.5", "0.75", "1"])

    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved {out_path}")


# ═══════════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    args = parse_args()
    # Re-bind the module-level constant from the CLI flag so all existing
    # references pick up the override without threading through signatures.
    global MIN_CLUSTER_SIZE
    MIN_CLUSTER_SIZE = args.min_cluster_size

    if not os.path.exists(args.capri_scores) or os.path.getsize(args.capri_scores) == 0:
        save_fallback_plots("No HADDOCK results")
        sys.exit(0)

    data, header = parse_capri_tsv(args.capri_scores)
    if not data:
        save_fallback_plots("No data in CAPRI scores")
        sys.exit(0)

    scores = get_numeric_col(data, ["score", "haddock-score", "total"])
    bsas = get_numeric_col(data, ["bsa", "BSA", "buried_surface_area"])

    if scores is None:
        save_fallback_plots("No HADDOCK score column found")
        sys.exit(0)

    best_score = scores[0]

    # ── Plot 1: Score vs BSA ─────────────────────────────────────────────
    if bsas is not None:
        plot_score_vs_bsa(scores, bsas, best_score, bsas[0])
    else:
        make_empty_plot("BSA column not found\nin capri_ss.tsv", "haddock_score_vs_bsa.png")
        print("WARNING: BSA column not found in capri_ss.tsv")

    # ── Plot 2: Cluster sizes ────────────────────────────────────────────
    clusters = {}
    best_qual_cluster = None
    if args.cluster_summary and os.path.exists(args.cluster_summary):
        clusters, _, best_qual_cluster = parse_clustfcc_tsv(args.cluster_summary, MIN_CLUSTER_SIZE)
    else:
        print("WARNING: No cluster summary file provided or found")

    # Largest named cluster (for backwards compat — not used by star placement)
    named = {c: m for c, m in clusters.items() if c != "-"}
    best_cluster = max(named, key=lambda c: len(named[c])) if named else None

    plot_cluster_sizes(clusters, best_cluster, best_qual_cluster)

    # ── Plot 3: Interface heatmap ────────────────────────────────────────
    pdb_index = collect_pdb_index(complex_dir=args.complex_dir, run_dir=args.run_dir)
    if not pdb_index:
        src = args.run_dir or args.complex_dir
        if src:
            print(f"WARNING: No docked PDB files found under '{src}'")
        else:
            print("WARNING: --run-dir or --complex-dir not provided; heatmap will be empty.")

    fixed_residues = _parse_fixed_residues(args.contigs, args.receptor_chain)
    if fixed_residues is not None:
        print(f"Fixed receptor residues from contig: {len(fixed_residues)}")
    else:
        print("No contig string provided; design region shading disabled.")

    # Best model: lowest-scoring member of pipeline-selected cluster
    best_model_name = None
    if best_qual_cluster and best_qual_cluster in clusters:
        scored = [(m, s) for m, s in clusters[best_qual_cluster] if s is not None]
        if scored:
            best_model_name = min(scored, key=lambda x: x[1])[0]
    if best_model_name is None:
        best_model_name = data[0].get("model", data[0].get("structure", data[0].get("pdb")))

    # Parse effector active residues for annotation bar
    eff_active = []
    if args.effector_active_residues:
        for token in args.effector_active_residues.split(","):
            token = token.strip()
            if not token:
                continue
            if "-" in token:
                parts = token.split("-", 1)
                try:
                    eff_active.extend(range(int(parts[0]), int(parts[1]) + 1))
                except ValueError:
                    pass
            else:
                try:
                    eff_active.append(int(token))
                except ValueError:
                    pass
    if eff_active:
        print(f"Effector active residues for annotation: {len(eff_active)}")

    plot_interface_heatmap(clusters, pdb_index, args.receptor_chain, args.effector_chain,
                           args.contact_cutoff, fixed_residues=fixed_residues,
                           best_model_name=best_model_name,
                           effector_active_residues=sorted(set(eff_active)) if eff_active else None)


if __name__ == "__main__":
    main()