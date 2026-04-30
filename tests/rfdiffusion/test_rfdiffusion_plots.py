#!/usr/bin/env python3
"""
test_rfdiffusion_plots.py
-------------------------
Standalone iteration harness for ``rfdiffusion_plots.py``.

Reads the metrics JSON produced by a previous ``test_rfdiffusion`` run
and re-emits the five RFDiffusion plots in an output directory of
choice.  No GPU, no Nextflow — designed for fast iteration of the
plot code itself against cached test outputs.

The three plots that already look correct are produced by importing
the production functions unchanged:

    rfdiff_specificity_coverage.png
    rfdiff_design_lengths.png
    rfdiff_com_displacement.png

The two plots with the dendrogram-placement issue
(``rfdiff_contact_map`` and ``rfdiff_design_clustering``) are rendered
using a single calibrated layout, picked after a multi-variant
exploration phase.  Strategy:

    - Heatmap drawn at an explicit fig-coord position with normal
      left-side y-tick labels.
    - Dendrogram overlaid as an inset axis to the left of those
      labels.  Inset position is calibrated to the MEASURED width of
      the longest label, so the visual gap between labels and
      dendrogram is consistent regardless of label-text length and
      regardless of how many panels share the figure.
    - The clustering heatmap uses ``aspect='equal'`` so cells are
      square, with a dedicated colorbar axis to the right (rather
      than the auto-colorbar that steals space from the heatmap).
    - Figure dimensions are absolute (inches) and ``bbox_inches`` is
      not set on save, so layout is fully under our control and
      nothing crops.

Usage
-----
    python3 test_rfdiffusion_plots.py \\
        --metrics tests/rfdiffusion/receptor_resurfacing_results/rfdiffusion/rfdiffusion_metrics.json \\
        --outdir  tests/rfdiffusion/receptor_resurfacing_results/plots_iter

If ``--metrics`` is omitted the script falls back to the canonical
test output path.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import sys
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import MaxNLocator

from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform


# ─── Colours (must match production rfdiffusion_plots.py) ───────────────────
COLOUR_PASS     = "#1F77B4"
COLOUR_FAIL     = "#D62728"
COLOUR_DESIGN   = "#FFD54F"
COLOUR_SUMMARY  = "#4C72B0"
COLOUR_COVERAGE = "#FF7F0E"


# ═══════════════════════════════════════════════════════════════════════════
# Production-script import — used to render the three "good" plots
# (specificity_coverage, design_lengths, com_displacement) unchanged.
# ═══════════════════════════════════════════════════════════════════════════

def _import_production_plots(prod_path: Path):
    """Import rfdiffusion_plots.py as a module by file path.

    Imported lazily so the script can still produce the variant
    plots even if the production file isn't on the same machine.
    """
    spec = importlib.util.spec_from_file_location("rfdiffusion_plots_prod", prod_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"could not load {prod_path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ═══════════════════════════════════════════════════════════════════════════
# Helpers (lifted from rfdiffusion_plots.py — keep in sync)
# ═══════════════════════════════════════════════════════════════════════════

def _design_label(d):
    return d["design"].replace("design_", "").replace(".pdb", "")


def _get_per_design_residues(d, global_design_residues):
    per = d.get("per_design_design_residues")
    if per and len(per) > 0:
        return set(per)
    return global_design_residues


def _get_per_region_lengths(d):
    prl = d.get("per_region_lengths")
    return prl if prl and len(prl) > 0 else []


def _label_with_region_len(d, region_idx, global_design_residues):
    lbl = _design_label(d)
    prl = _get_per_region_lengths(d)
    if prl and region_idx < len(prl):
        region_len = prl[region_idx]
    else:
        region_len = len(_get_per_design_residues(d, global_design_residues))
    return f"d{lbl} ({region_len})"


def _get_dendrogram_order(designs):
    """Average-linkage clustering on flat design_region_coords."""
    valid = [(i, d) for i, d in enumerate(designs)
             if d.get("design_region_coords") and len(d["design_region_coords"]) > 0]
    if len(valid) < 2:
        return list(range(len(designs))), None, None

    coord_arrays = [np.array(d["design_region_coords"]) for _, d in valid]
    min_len = min(len(c) for c in coord_arrays)
    if min_len == 0:
        return list(range(len(designs))), None, None
    trimmed = [c[:min_len] for c in coord_arrays]
    n = len(trimmed)

    rmsd_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            diff = trimmed[i] - trimmed[j]
            rmsd = np.sqrt(np.mean(np.sum(diff ** 2, axis=1)))
            rmsd_matrix[i, j] = rmsd
            rmsd_matrix[j, i] = rmsd

    condensed = squareform(rmsd_matrix)
    Z = linkage(condensed, method="average")

    temp_fig, temp_ax = plt.subplots()
    dendro_ref = dendrogram(Z, ax=temp_ax, no_labels=True)
    order = dendro_ref["leaves"]
    plt.close(temp_fig)

    orig_order = [valid[k][0] for k in order]
    return orig_order, Z, rmsd_matrix


def _per_region_clustering_data(designs, region_idx):
    """Return {valid, Z, rmsd_matrix, order, n} or None for one region."""
    def _get_region_coords(d, r_idx):
        prc = d.get("per_region_coords")
        if prc and r_idx < len(prc):
            return prc[r_idx]
        if r_idx == 0:
            return d.get("design_region_coords", [])
        return []

    valid = [(i, d) for i, d in enumerate(designs)
             if len(_get_region_coords(d, region_idx)) > 0]
    if len(valid) < 2:
        return None

    coord_arrays = [np.array(_get_region_coords(d, region_idx)) for _, d in valid]
    min_len = min(len(c) for c in coord_arrays)
    if min_len == 0:
        return None
    trimmed = [c[:min_len] for c in coord_arrays]
    n = len(trimmed)

    rmsd_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            diff = trimmed[i] - trimmed[j]
            rmsd = np.sqrt(np.mean(np.sum(diff ** 2, axis=1)))
            rmsd_matrix[i, j] = rmsd
            rmsd_matrix[j, i] = rmsd

    condensed = squareform(rmsd_matrix)
    Z = linkage(condensed, method="average")

    temp_fig, temp_ax = plt.subplots()
    dendro_ref = dendrogram(Z, ax=temp_ax, no_labels=True)
    order = dendro_ref["leaves"]
    plt.close(temp_fig)

    return {"valid": valid, "Z": Z, "rmsd_matrix": rmsd_matrix,
            "order": order, "n": n}


def _n_regions(designs):
    """How many design regions are present in the metrics."""
    for d in designs:
        prc = d.get("per_region_coords")
        if prc and len(prc) > 0:
            return len(prc)
    return 1


def _build_contact_map_data(designs, global_design_residues, fixed_residues):
    """Compute the per-design data needed for plot_contact_map.

    Returned in dendrogram order.  Mirrors the production
    rfdiffusion_plots.plot_contact_map preparation block exactly so
    the variant renderers can share one data path.
    """
    dendro_order, Z, _ = _get_dendrogram_order(designs)
    sorted_designs = [designs[i] for i in dendro_order]
    n_designs = len(sorted_designs)

    design_contact_positions = []
    design_region_positions = []
    design_n_res = []

    for d in sorted_designs:
        per_design_dr = _get_per_design_residues(d, global_design_residues)

        pos_order = d.get("receptor_position_order")
        if pos_order and len(pos_order) > 0:
            all_resnums = pos_order
        else:
            all_resnums = sorted(fixed_residues | per_design_dr)

        resnum_to_pos = {r: i for i, r in enumerate(all_resnums)}
        n_res = len(all_resnums)
        design_n_res.append(n_res)

        dr_positions = set(resnum_to_pos[r] for r in per_design_dr if r in resnum_to_pos)
        design_region_positions.append(dr_positions)

        contacts = {}
        for r in d.get("receptor_contact_residues", []):
            if r in resnum_to_pos:
                pos = resnum_to_pos[r]
                contacts[pos] = "in" if r in per_design_dr else "out"
        design_contact_positions.append(contacts)

    return {
        "sorted_designs":    sorted_designs,
        "n_designs":         n_designs,
        "design_n_res":      design_n_res,
        "design_region_positions": design_region_positions,
        "design_contact_positions": design_contact_positions,
        "Z":                 Z,
    }


def _build_effector_freq(designs, n_designs):
    all_eff_res = set()
    for d in designs:
        all_eff_res.update(d.get("effector_contact_residues", []))
    if not all_eff_res:
        return None

    eff_min_r = min(all_eff_res)
    eff_max_r = max(all_eff_res)
    eff_range = list(range(eff_min_r, eff_max_r + 1))
    n_eff = len(eff_range)
    eff_to_idx = {e: i for i, e in enumerate(eff_range)}
    eff_freq = np.zeros(n_eff)
    for d in designs:
        for e in d.get("effector_contact_residues", []):
            if e in eff_to_idx:
                eff_freq[eff_to_idx[e]] += 1
    if n_designs > 0:
        eff_freq /= n_designs

    return {"range": eff_range, "freq": eff_freq, "n_eff": n_eff}


def _build_design_labels(sorted_designs, global_design_residues):
    """Y-axis labels for the contact map: 'd<i> (l1+l2+...)'."""
    labels = []
    for d in sorted_designs:
        lbl = _design_label(d)
        prl = _get_per_region_lengths(d)
        if prl:
            len_str = "+".join(str(l) for l in prl)
        else:
            len_str = str(len(_get_per_design_residues(d, global_design_residues)))
        labels.append(f"d{lbl} ({len_str})")
    return labels


# ═══════════════════════════════════════════════════════════════════════════
# Shared drawing routines — all variants call these so the heatmap content
# itself is identical across variants; only the dendrogram + label
# placement differs.
# ═══════════════════════════════════════════════════════════════════════════

def _draw_contact_heatmap_content(ax_heat, data, max_n_res):
    """Draw the contact-map heatmap content (hatching, design region
    backgrounds, contact blocks) into the given axis.  Sets x/y limits
    and inverts y-axis but does NOT touch tick labels — variants set
    those themselves."""
    n_designs = data["n_designs"]
    design_n_res = data["design_n_res"]
    design_region_positions = data["design_region_positions"]
    design_contact_positions = data["design_contact_positions"]

    # Beyond-protein-end hatching per design
    for i in range(n_designs):
        n_res = design_n_res[i]
        if n_res < max_n_res:
            ax_heat.add_patch(plt.Rectangle(
                (n_res - 0.5, i - 0.5), (max_n_res - n_res), 1.0,
                facecolor="#E0E0E0", alpha=0.5, edgecolor="#888888",
                linewidth=0, zorder=0, hatch="///",
            ))

    # Design-region backgrounds with black boundary lines
    for i in range(n_designs):
        dr_pos = sorted(design_region_positions[i])
        if not dr_pos:
            continue
        runs = []
        run_start = dr_pos[0]
        prev = run_start
        for p in dr_pos[1:]:
            if p != prev + 1:
                runs.append((run_start, prev))
                run_start = p
            prev = p
        runs.append((run_start, prev))
        for rs, re in runs:
            ax_heat.add_patch(plt.Rectangle(
                (rs - 0.5, i - 0.5), (re - rs + 1), 1.0,
                facecolor=COLOUR_DESIGN, alpha=0.18, edgecolor="black",
                linewidth=0.8, zorder=0,
            ))

    # Contact blocks
    for i in range(n_designs):
        for pos, cat in design_contact_positions[i].items():
            colour = COLOUR_PASS if cat == "in" else COLOUR_FAIL
            ax_heat.add_patch(plt.Rectangle(
                (pos - 0.5, i - 0.5), 1.0, 1.0,
                facecolor=colour, edgecolor="none", zorder=1,
            ))

    ax_heat.set_xlim(-0.5, max_n_res - 0.5)
    ax_heat.set_ylim(-0.5, n_designs - 0.5)
    ax_heat.invert_yaxis()


def _set_contact_xticks(ax_heat, data, max_n_res, global_design_residues, fixed_residues):
    """Receptor-position x-axis ticks for the contact-map heatmap."""
    sorted_designs = data["sorted_designs"]
    design_n_res = data["design_n_res"]
    longest_idx = design_n_res.index(max_n_res)
    ref_design = sorted_designs[longest_idx]
    ref_pos_order = ref_design.get("receptor_position_order")
    ref_dr = _get_per_design_residues(ref_design, global_design_residues)
    if ref_pos_order and len(ref_pos_order) > 0:
        ref_all = ref_pos_order
    else:
        ref_all = sorted(fixed_residues | ref_dr)

    tick_step = max(1, len(ref_all) // 25)
    xtick_positions = list(range(0, len(ref_all), tick_step))
    # Only append the final tick if it's at least one full tick_step
    # away from the previous one — otherwise the last two labels
    # collide (this was the "82/83" / "131/132" overlap).  When too
    # close, replace the previous-to-last so we still show the
    # endpoint without doubling up.
    last_idx = len(ref_all) - 1
    if last_idx not in xtick_positions:
        if not xtick_positions or (last_idx - xtick_positions[-1]) >= tick_step:
            xtick_positions.append(last_idx)
        else:
            xtick_positions[-1] = last_idx
    xtick_labels = [str(i + 1) for i in xtick_positions]

    ax_heat.set_xticks(xtick_positions)
    ax_heat.set_xticklabels(xtick_labels, fontsize=7, rotation=0, ha="center")
    ax_heat.set_xlabel("Receptor position", fontsize=9)


def _draw_effector_freq_bar(ax_eff, eff_data):
    if eff_data is None:
        return
    eff_range = eff_data["range"]
    eff_freq = eff_data["freq"]
    n_eff = eff_data["n_eff"]
    ax_eff.bar(range(n_eff), eff_freq, width=1.0,
               color=COLOUR_SUMMARY, edgecolor="none", alpha=0.8)
    ax_eff.set_xlim(-0.5, n_eff - 0.5)
    max_freq = max(eff_freq) if len(eff_freq) > 0 else 1.0
    ax_eff.set_ylim(0, max(1.05, max_freq * 1.1))
    tick_step = max(1, max(5, round(n_eff / 15 / 5) * 5))
    xticks = list(range(0, n_eff, tick_step))
    last_idx = n_eff - 1
    if last_idx not in xticks:
        if not xticks or (last_idx - xticks[-1]) >= tick_step:
            xticks.append(last_idx)
        else:
            xticks[-1] = last_idx
    ax_eff.set_xticks(xticks)
    ax_eff.set_xticklabels([str(eff_range[i]) for i in xticks], fontsize=7)
    ax_eff.set_xlabel("Effector residue", fontsize=9)
    ax_eff.set_ylabel("Freq.", fontsize=8)
    ax_eff.tick_params(axis="y", labelsize=7)


def _contact_legend_handles():
    return [
        mpatches.Patch(color=COLOUR_PASS, label="Contact in design region"),
        mpatches.Patch(color=COLOUR_FAIL, label="Contact outside design region"),
        mpatches.Patch(facecolor=COLOUR_DESIGN, alpha=0.3, edgecolor="black",
                       linewidth=0.8, label="Design region"),
        mpatches.Patch(facecolor="#E0E0E0", alpha=0.5, hatch="///",
                       edgecolor="#999999", label="Beyond protein end"),
    ]


def _draw_clustering_heatmap_content(ax_heat, ax_dendro, rd, n_regions, r_idx,
                                      global_design_residues, draw_labels=True,
                                      label_side="left", aspect="auto",
                                      draw_colorbar=True):
    """Draw clustering-heatmap content into ax_heat.  If draw_labels is
    True, also sets y/x ticklabels with pass/fail colour-coding on the
    given side ('left' or 'right').

    aspect:
        'auto'  — heatmap fills the axis (default, used by old variants)
        'equal' — square cells, used by the calibrated variant

    draw_colorbar:
        True  — create the colorbar inside the helper using
                ax.figure.colorbar (steals space from ax_heat)
        False — skip; the caller is expected to add a dedicated cbar
                axis via add_axes()  (returns im for caller use)
    """
    n = rd["n"]
    order = rd["order"]
    Z_panel = rd["Z"]
    rmsd_matrix = rd["rmsd_matrix"]
    valid = rd["valid"]

    reordered = rmsd_matrix[np.ix_(order, order)]
    reordered_flipped = reordered[:, ::-1]
    flipped_col_labels = list(reversed(order))

    heat_size = max(3, n * 0.5)
    label_fontsize = max(7, min(12, heat_size * 1.2))
    cell_fontsize = max(5, min(10, heat_size * 0.9 - n // 5))
    legend_fontsize = max(8, min(11, heat_size * 1.0))
    cbar_fontsize = max(9, min(12, heat_size * 1.0))

    vmax = max(1.0, np.max(reordered_flipped))
    im = ax_heat.imshow(reordered_flipped, cmap="viridis_r", vmin=0, vmax=vmax,
                        aspect=aspect, origin="lower")

    if ax_dendro is not None:
        dendrogram(Z_panel, ax=ax_dendro, orientation="left",
                   no_labels=True, color_threshold=0,
                   above_threshold_color="#555555")
        ax_dendro.set_axis_off()

    if n <= 15:
        for i in range(n):
            for j in range(n):
                val = reordered_flipped[i, j]
                colour = "white" if val > vmax * 0.6 else "black"
                ax_heat.text(j, i, f"{val:.1f}", ha="center", va="center",
                             fontsize=cell_fontsize, color=colour)

    ordered_designs_rows = [valid[k][1] for k in order]
    row_labels = [_label_with_region_len(d, r_idx, global_design_residues)
                  for d in ordered_designs_rows]
    row_colours = [COLOUR_PASS if d["passes_filter"] else COLOUR_FAIL
                   for d in ordered_designs_rows]

    ordered_designs_cols = [valid[k][1] for k in flipped_col_labels]
    col_labels = [_label_with_region_len(d, r_idx, global_design_residues)
                  for d in ordered_designs_cols]
    col_colours = [COLOUR_PASS if d["passes_filter"] else COLOUR_FAIL
                   for d in ordered_designs_cols]

    if draw_labels:
        ax_heat.set_yticks(range(n))
        ax_heat.set_yticklabels(row_labels, fontsize=label_fontsize)
        for i, ytl in enumerate(ax_heat.get_yticklabels()):
            ytl.set_color(row_colours[i])
        if label_side == "right":
            ax_heat.yaxis.tick_right()
            ax_heat.yaxis.set_label_position("right")
    else:
        ax_heat.set_yticks([])

    ax_heat.set_xticks(range(n))
    rotate = n > 8 or max((len(l) for l in col_labels), default=0) > 8
    ax_heat.set_xticklabels(col_labels,
                            rotation=45 if rotate else 0,
                            ha="right" if rotate else "center",
                            fontsize=label_fontsize)
    for i, xtl in enumerate(ax_heat.get_xticklabels()):
        xtl.set_color(col_colours[i])

    if draw_colorbar:
        cbar = ax_heat.figure.colorbar(im, ax=ax_heat, pad=0.04, fraction=0.046,
                                        aspect=20, shrink=0.85)
        cbar.set_label("Design-region Cα RMSD (Å)", fontsize=cbar_fontsize)
        cbar.ax.tick_params(labelsize=max(7, cbar_fontsize - 2))

    region_label = f"Design region {r_idx + 1}" if n_regions > 1 else "Design region"
    ax_heat.set_title(region_label, fontsize=10, loc="left", pad=8)

    legend_handles = [
        mpatches.Patch(color=COLOUR_PASS, label="Pass"),
        mpatches.Patch(color=COLOUR_FAIL, label="Filtered"),
    ]
    # Sit above heatmap, right-aligned (matches contact-map legend style).
    ax_heat.legend(handles=legend_handles, fontsize=legend_fontsize,
                   loc="lower right", bbox_to_anchor=(1.0, 1.01),
                   ncol=2, frameon=True, fancybox=True,
                   edgecolor="#cccccc")

    return {
        "row_labels":  row_labels,
        "row_colours": row_colours,
        "col_labels":  col_labels,
        "col_colours": col_colours,
        "n":           n,
        "order":       order,
        "im":          im,
        "vmax":        vmax,
        "label_fontsize": label_fontsize,
        "cbar_fontsize": cbar_fontsize,
    }




# ═══════════════════════════════════════════════════════════════════════════
# Label-width measurement
# ═══════════════════════════════════════════════════════════════════════════

def _measure_label_width_inches(labels, fontsize, dpi=100):
    """Render the longest label to a throwaway figure and measure its
    width in inches.  Used to calibrate the dendrogram-inset offset
    and the heatmap left margin so labels neither overlap the
    dendrogram nor leave a wasteful gap.

    The throwaway figure uses dpi=100 so pixel-to-inch conversion is
    independent of the final figure dpi (we want a true label width).
    """
    if not labels:
        return 0.0
    longest = max(labels, key=len)
    fig = plt.figure(figsize=(2, 1), dpi=dpi)
    ax = fig.add_subplot(111)
    txt = ax.text(0.0, 0.0, longest, fontsize=fontsize)
    fig.canvas.draw()
    bbox_px = txt.get_window_extent()
    width_in = bbox_px.width / dpi
    plt.close(fig)
    return width_in


# ═══════════════════════════════════════════════════════════════════════════
# Calibrated rendering — single variant for both plots, replacing
# B/C/F/G/H.  Strategy: G-style (heatmap with normal left-side labels,
# dendrogram overlaid as an inset axis) but with the dendrogram-inset
# horizontal position calibrated to the measured longest-label width
# so there's a consistent visual gap between labels and dendrogram
# regardless of label length.
# ═══════════════════════════════════════════════════════════════════════════

# Tunables — single source of truth for both plots.
DENDRO_W_INCHES        = 1.6     # width of the dendrogram inset
LABEL_DENDRO_GAP_IN    = 0.20    # visual gap between labels and dendrogram
RIGHT_MARGIN_IN        = 1.6     # space for colour bar + label, contact map plot's legend, etc.
TOP_MARGIN_IN          = 0.7
BOTTOM_MARGIN_IN       = 0.7


def plot_contact_map(data, eff_data, max_n_res, labels,
                      global_design_residues, fixed_residues, out_path):
    """Calibrated contact-map render — heatmap with normal left-side
    labels and dendrogram overlaid as an inset axis to the left,
    spaced consistently from the labels."""
    n_designs = data["n_designs"]
    Z = data["Z"]
    has_dendro = Z is not None
    has_eff = eff_data is not None

    label_fontsize = max(5, 8 - n_designs // 15)
    label_w_in = _measure_label_width_inches(labels, label_fontsize)

    # Heatmap dimensions in inches.
    heat_w_in = max(8, max_n_res * 0.08 + 2)
    heat_h_in = max(3, 0.8 + n_designs * 0.35)
    eff_h_in = 1.5 if has_eff else 0
    eff_gap_in = 0.4 if has_eff else 0

    # Left margin = dendrogram width + dendro-label gap + label width.
    left_margin_in = (DENDRO_W_INCHES + LABEL_DENDRO_GAP_IN + label_w_in
                      if has_dendro else label_w_in + 0.2)

    fig_w_in = left_margin_in + heat_w_in + RIGHT_MARGIN_IN
    fig_h_in = TOP_MARGIN_IN + heat_h_in + eff_gap_in + eff_h_in + BOTTOM_MARGIN_IN

    fig = plt.figure(figsize=(fig_w_in, fig_h_in))

    # Heatmap axis — explicit fig-coord position.
    heat_y0 = (BOTTOM_MARGIN_IN + eff_h_in + eff_gap_in) / fig_h_in
    heat_h = heat_h_in / fig_h_in
    heat_x0 = left_margin_in / fig_w_in
    heat_w = heat_w_in / fig_w_in
    ax_heat = fig.add_axes([heat_x0, heat_y0, heat_w, heat_h])

    if has_eff:
        eff_y0 = BOTTOM_MARGIN_IN / fig_h_in
        eff_h_frac = eff_h_in / fig_h_in
        ax_eff = fig.add_axes([heat_x0, eff_y0, heat_w, eff_h_frac])
    else:
        ax_eff = None

    _draw_contact_heatmap_content(ax_heat, data, max_n_res)
    _set_contact_xticks(ax_heat, data, max_n_res, global_design_residues, fixed_residues)
    ax_heat.set_yticks(range(n_designs))
    ax_heat.set_yticklabels(labels, fontsize=label_fontsize)

    # Dendrogram inset, positioned so its right edge sits LABEL_DENDRO_GAP_IN
    # to the left of the labels' left edge.
    if has_dendro:
        dendro_w_frac = DENDRO_W_INCHES / fig_w_in
        # Dendrogram's right edge sits at:
        #   heat_x0 - (label_w + LABEL_DENDRO_GAP_IN) / fig_w_in
        dendro_x1 = heat_x0 - (label_w_in + LABEL_DENDRO_GAP_IN) / fig_w_in
        dendro_x0 = dendro_x1 - dendro_w_frac
        ax_dendro = fig.add_axes([dendro_x0, heat_y0, dendro_w_frac, heat_h])
        dendrogram(Z, ax=ax_dendro, orientation="left",
                   no_labels=True, color_threshold=0,
                   above_threshold_color="#555555")
        ax_dendro.set_axis_off()

    ax_heat.legend(handles=_contact_legend_handles(), fontsize=7,
                   loc="lower center", bbox_to_anchor=(0.5, 1.01),
                   ncol=4, frameon=True, fancybox=True)
    if ax_eff is not None:
        _draw_effector_freq_bar(ax_eff, eff_data)

    plt.savefig(out_path, dpi=150)  # bbox_inches=None — fixed layout
    plt.close(fig)


def plot_design_clustering(designs, global_design_residues, out_path):
    """Calibrated design-clustering render — multi-panel, square
    heatmaps with dedicated colorbar axes and dendrogram insets
    sized by the actual label width.

    Each panel is laid out in absolute inches:
      [dendro | gap | labels | square heatmap | gap | colorbar]
    so that:
      - heatmap is square (cell aspect = 1)
      - labels never collide with the dendrogram (gap is consistent)
      - colorbar label has space (right margin reserves room)
    """
    n_regions = _n_regions(designs)
    region_data = [_per_region_clustering_data(designs, r) for r in range(n_regions)]
    plottable = [(i, rd) for i, rd in enumerate(region_data) if rd is not None]
    if not plottable:
        # Caller (orchestration) handles the empty case.
        return False

    # First pass: compute label widths per panel and panel sizes.
    panel_specs = []
    for r_idx, rd in plottable:
        n = rd["n"]
        # Match _draw_clustering_heatmap_content's font-size logic.
        heat_size = max(3, n * 0.5)
        label_fontsize = max(7, min(12, heat_size * 1.2))
        # Build labels here just to measure.  Cheap (just string ops).
        ordered = [rd["valid"][k][1] for k in rd["order"]]
        labels = [_label_with_region_len(d, r_idx, global_design_residues)
                  for d in ordered]
        label_w_in = _measure_label_width_inches(labels, label_fontsize)

        # Heatmap is square: side length proportional to n.
        cell_in = 0.5
        heat_side_in = max(3.0, n * cell_in)

        panel_specs.append({
            "r_idx":          r_idx,
            "rd":             rd,
            "label_w_in":     label_w_in,
            "label_fontsize": label_fontsize,
            "heat_side_in":   heat_side_in,
        })

    # Use the WIDEST label width across panels to fix the heatmap left
    # edge — heatmaps line up vertically across panels for easy
    # cross-region comparison.  Each panel's dendrogram, however, is
    # positioned relative to ITS OWN labels (right edge of dendro =
    # GAP inches from that panel's longest label), so a panel with
    # shorter labels gets its dendrogram drifted rightwards, keeping
    # the visual label-to-dendro gap consistent everywhere.
    common_label_w_in = max(p["label_w_in"] for p in panel_specs)

    # Use the LARGEST heat side so all panels share the same width.
    common_heat_side_in = max(p["heat_side_in"] for p in panel_specs)

    # Per-panel x-layout in inches.  Tightened from initial pass to
    # match the user's "less whitespace" feedback.
    panel_left_margin_in  = 0.2
    cbar_gap_in           = 0.2    # heatmap → cbar (was 0.4)
    cbar_w_in             = 0.25
    cbar_label_w_in       = 0.5    # rotated label "Design-region Cα RMSD (Å)" (was 0.9)
    panel_right_margin_in = 0.2

    panel_w_in = (panel_left_margin_in + DENDRO_W_INCHES + LABEL_DENDRO_GAP_IN
                  + common_label_w_in + common_heat_side_in + cbar_gap_in
                  + cbar_w_in + cbar_label_w_in + panel_right_margin_in)

    # Vertical stacking — also tightened.
    xtick_room_in   = 0.4    # was 0.5; xticks aren't rotated, plenty of room
    title_room_in   = 0.25   # was 0.3
    legend_room_in  = 0.25   # was 0.35; pass/filtered legend is small
    panel_h_in      = (common_heat_side_in + xtick_room_in + title_room_in
                       + legend_room_in)
    panel_gap_in    = 0.2    # was 0.3

    # Local margins (override the module-level constants which are
    # tuned for the contact-map plot).
    top_margin_in    = 0.3
    bottom_margin_in = 0.3

    fig_w_in = panel_w_in
    fig_h_in = (top_margin_in + panel_h_in * len(panel_specs)
                + panel_gap_in * (len(panel_specs) - 1) + bottom_margin_in)

    fig = plt.figure(figsize=(fig_w_in, fig_h_in))

    # Heatmap x-position is fixed across panels (uses common_label_w_in)
    heat_x_in = (panel_left_margin_in + DENDRO_W_INCHES + LABEL_DENDRO_GAP_IN
                 + common_label_w_in)
    cbar_x_in = heat_x_in + common_heat_side_in + cbar_gap_in

    from matplotlib.ticker import FormatStrFormatter

    for k, spec in enumerate(panel_specs):
        rd = spec["rd"]
        r_idx = spec["r_idx"]
        n = rd["n"]
        panel_label_w_in = spec["label_w_in"]

        # Vertical position: top-most panel first.
        panel_top_in = fig_h_in - top_margin_in - k * (panel_h_in + panel_gap_in)
        panel_y_in = panel_top_in - panel_h_in   # bottom of panel
        # Heatmap sits above the x-tick room.
        heat_y_in = panel_y_in + xtick_room_in
        heat_y_frac = heat_y_in / fig_h_in
        heat_h_frac = common_heat_side_in / fig_h_in

        heat_x_frac = heat_x_in / fig_w_in
        heat_w_frac = common_heat_side_in / fig_w_in
        ax_heat = fig.add_axes([heat_x_frac, heat_y_frac, heat_w_frac, heat_h_frac])

        drawn = _draw_clustering_heatmap_content(
            ax_heat, None, rd, n_regions, r_idx,
            global_design_residues, draw_labels=True, label_side="left",
            aspect="equal", draw_colorbar=False)

        # Dendrogram inset — right edge sits LABEL_DENDRO_GAP_IN to the
        # left of THIS PANEL's labels (per-panel calibration).
        dendro_w_frac = DENDRO_W_INCHES / fig_w_in
        dendro_x_in_panel = (heat_x_in - panel_label_w_in
                             - LABEL_DENDRO_GAP_IN - DENDRO_W_INCHES)
        ax_dendro = fig.add_axes([
            dendro_x_in_panel / fig_w_in, heat_y_frac,
            dendro_w_frac, heat_h_frac,
        ])
        dendrogram(rd["Z"], ax=ax_dendro, orientation="left",
                   no_labels=True, color_threshold=0,
                   above_threshold_color="#555555")
        ax_dendro.set_axis_off()

        # Colorbar in its own axis to the right of the heatmap.
        cbar_w_frac = cbar_w_in / fig_w_in
        cax = fig.add_axes([cbar_x_in / fig_w_in, heat_y_frac,
                            cbar_w_frac, heat_h_frac])
        cbar = fig.colorbar(drawn["im"], cax=cax)
        cbar.set_label("Design-region Cα RMSD (Å)",
                       fontsize=drawn["cbar_fontsize"])
        cbar.ax.tick_params(labelsize=max(7, drawn["cbar_fontsize"] - 2))
        # Format ticks as 0.0, 1.0, 2.0 instead of 0, 1, 2.
        cbar.ax.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))

    plt.savefig(out_path, dpi=150)
    plt.close(fig)
    return True


# ═══════════════════════════════════════════════════════════════════════════
# Design lengths — copied from production for iteration in this script
# ═══════════════════════════════════════════════════════════════════════════

def plot_design_lengths(designs, global_design_residues, fixed_residues, out_path):
    """Multi-panel bar chart of design-region lengths, one panel per
    contiguous design region.

    Differences from production:
    - Bars get thin black borders.
    - A universal bar width (in inches) is used; each panel's width
      scales to the number of unique length values, so a region with
      one unique value gets a narrow panel rather than one massive
      stretched bar.
    - Per-panel x-ticks are set explicitly without MaxNLocator (which
      previously over-rode the explicit set_xticks and produced too
      many ticks).
    """
    if not designs:
        return False

    all_region_lengths = []
    for d in designs:
        prl = d.get("per_region_lengths")
        if prl and len(prl) > 0:
            for k, length in enumerate(prl):
                while len(all_region_lengths) <= k:
                    all_region_lengths.append([])
                all_region_lengths[k].append(length)

    if not all_region_lengths:
        return False

    n_regions = len(all_region_lengths)

    # Per-region pre-computation so we can size each panel by its
    # unique-value count.
    panel_specs = []
    for k in range(n_regions):
        lengths = np.array(all_region_lengths[k])
        unique_lengths, counts = np.unique(lengths, return_counts=True)
        # Show every integer in [min, max]; one-bar regions get one
        # tick total (no padding ints).
        lo, hi = int(unique_lengths.min()), int(unique_lengths.max())
        all_ints = list(range(lo, hi + 1))
        panel_specs.append({
            "lengths": lengths,
            "unique":  unique_lengths,
            "counts":  counts,
            "all_ints": all_ints,
            "n_bars":  len(all_ints),
        })

    # Layout constants (same target physical bar width is used in
    # the per-panel sizing block below).
    target_bar_w_in_setup = 0.6
    bar_width_data_setup  = 0.8
    edge_pad_data_setup   = 0.4
    min_panel_w_in_setup  = 3.0
    panel_padding_in = 1.2   # horizontal margin per panel for ylabel/cell pad
    panel_h_in = 3.0
    panel_gap_in = 0.4

    # Per-panel natural widths: target_data_range * bar_w / data_w.
    panel_widths_in = []
    for s in panel_specs:
        span_data = s["unique"][-1] - s["unique"][0]
        target_data_range = span_data + bar_width_data_setup + 2 * edge_pad_data_setup
        natural = target_data_range * target_bar_w_in_setup / bar_width_data_setup
        panel_widths_in.append(max(min_panel_w_in_setup, natural) + panel_padding_in)
    fig_w_in = max(panel_widths_in)
    fig_h_in = panel_h_in * n_regions + panel_gap_in * (n_regions - 1) + 1.0

    fig = plt.figure(figsize=(fig_w_in, fig_h_in))

    # Each panel gets its own gridspec row, but the panel doesn't fill
    # the row — it's left-anchored to a width matching n_bars.
    # Each panel is centered horizontally within the figure (not
    # left-anchored), so a one-bar panel sits in the middle rather
    # than at the left with a runaway xlabel.
    bottom_margin_in = 0.5
    top_margin_in = 0.3
    # Reserve a minimum panel width so the xlabel "Design-region length
    # (residues)" — about 2.6 inches wide — has room to fit.
    min_panel_w_in = 3.0
    # Target physical bar width.  We derive ax_w_in from the desired
    # data range, not the other way around — otherwise the first/last
    # bars get half-clipped by the xlim because their edges fall
    # outside it.
    target_bar_w_in = 0.6
    bar_width_data  = 0.8     # data-units; matches plt.bar(..., width=0.8)
    edge_pad_data   = 0.4     # half-bar of breathing room on each side

    for k in range(n_regions):
        spec = panel_specs[k]
        panel_top_in = fig_h_in - top_margin_in - k * (panel_h_in + panel_gap_in)
        panel_y_in = panel_top_in - panel_h_in
        ax_y_frac = (panel_y_in + 0.5) / fig_h_in
        ax_h_frac = (panel_h_in - 0.8) / fig_h_in

        # Data-coord range that exactly fits all bars + edge padding.
        span_data = spec["unique"][-1] - spec["unique"][0]   # 0 for one-bar
        target_data_range = span_data + bar_width_data + 2 * edge_pad_data
        # Inches needed to render each bar at target_bar_w_in physical.
        natural_ax_w_in = target_data_range * target_bar_w_in / bar_width_data
        ax_w_in = max(min_panel_w_in, natural_ax_w_in)
        # If the min-width floor inflated the axis above natural, the
        # bar would render wider than target.  Scale target_data_range
        # up to match — this gives more empty space around the bar
        # but keeps the bar's physical width identical across panels.
        if ax_w_in > natural_ax_w_in:
            target_data_range = bar_width_data * ax_w_in / target_bar_w_in
        ax_x_in = (fig_w_in - ax_w_in) / 2.0   # centered
        ax_x_frac = ax_x_in / fig_w_in
        ax_w_frac = ax_w_in / fig_w_in

        ax = fig.add_axes([ax_x_frac, ax_y_frac, ax_w_frac, ax_h_frac])

        ax.bar(spec["unique"], spec["counts"], width=bar_width_data,
               color=COLOUR_SUMMARY, edgecolor="black", linewidth=0.5,
               alpha=0.85)

        # X-ticks: every integer between min and max unique value
        # (inclusive), so gap-values get a tick + label too — but
        # never extend past the data range.
        lo, hi = int(spec["unique"][0]), int(spec["unique"][-1])
        all_ints_between = list(range(lo, hi + 1))
        ax.set_xticks(all_ints_between)
        ax.set_xticklabels([str(v) for v in all_ints_between])
        # xlim sized to fit all bars exactly with edge_pad on each side.
        center = (lo + hi) / 2.0
        ax.set_xlim(center - target_data_range / 2.0,
                    center + target_data_range / 2.0)

        ax.set_xlabel("Design-region length (residues)", fontsize=10)
        ax.set_ylabel("Count", fontsize=10)
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        ax.set_title(f"Design region {k + 1}", fontsize=11)

    plt.savefig(out_path, dpi=150)
    plt.close(fig)
    return True


# ═══════════════════════════════════════════════════════════════════════════
# COM displacement — copied from production for iteration in this script
# ═══════════════════════════════════════════════════════════════════════════

def plot_com_displacement(designs, out_path):
    """Per-design centre-of-mass displacement of each design region
    from the input gap residues.

    Differences from production:
    - Title removed (plot speaks for itself; the y-axis label is
      now expanded to carry the meaning).
    - Y-axis label is more verbose to stand alone without a title.
    """
    if not designs:
        return False

    plotted = []
    skipped = 0
    for d in designs:
        disps = d.get("design_region_com_displacement")
        if not disps:
            skipped += 1
            continue
        clean = [float(v) if (v is not None and not (isinstance(v, float) and np.isnan(v)))
                 else None for v in disps]
        if all(v is None for v in clean):
            skipped += 1
            continue
        plotted.append((d.get("design", "?"), clean))

    if not plotted:
        return False

    def _max_disp(item):
        vals = [v for v in item[1] if v is not None]
        return max(vals) if vals else -1.0
    plotted.sort(key=_max_disp)

    n_designs = len(plotted)
    n_regions = max(len(p[1]) for p in plotted)

    fig, ax = plt.subplots(figsize=(max(8, n_designs * 0.4 + 2), 5))

    bar_width = 0.8 / n_regions
    x_base = np.arange(n_designs)
    region_colours = plt.cm.viridis(np.linspace(0.15, 0.85, max(n_regions, 1)))

    for region_idx in range(n_regions):
        heights = []
        for _name, disps in plotted:
            v = disps[region_idx] if region_idx < len(disps) else None
            heights.append(v if v is not None else 0.0)
        offsets = x_base + (region_idx - (n_regions - 1) / 2) * bar_width
        label = f"region {region_idx + 1}" if n_regions > 1 else "design region"
        ax.bar(offsets, heights, width=bar_width,
               color=region_colours[region_idx], edgecolor="black",
               linewidth=0.5, label=label)

    ax.set_xticks(x_base)
    ax.set_xticklabels([p[0].replace(".pdb", "") for p in plotted],
                       rotation=45, ha="right", fontsize=8)
    # Verbose y-axis label (no title above).
    ax.set_ylabel("Design region COM displacement\nfrom input gap, scaffold-aligned (Å)",
                  fontsize=10)
    # Title removed deliberately.
    ax.grid(axis="y", linestyle=":", alpha=0.5)
    ax.set_axisbelow(True)
    if n_regions > 1:
        ax.legend(fontsize=8, frameon=False)

    if skipped:
        ax.text(0.99, 0.98,
                f"{skipped} design{'s' if skipped != 1 else ''} skipped (no data)",
                transform=ax.transAxes, ha="right", va="top",
                fontsize=8, color="#888")

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return True


# ═══════════════════════════════════════════════════════════════════════════
# Top-level orchestration
# ═══════════════════════════════════════════════════════════════════════════

def _resolve_metrics_path(supplied):
    if supplied:
        p = Path(supplied)
        if not p.exists():
            raise SystemExit(f"--metrics path does not exist: {p}")
        return p
    fallback = Path(
        "tests/rfdiffusion/receptor_resurfacing_results/rfdiffusion/rfdiffusion_metrics.json"
    )
    if fallback.exists():
        return fallback
    raise SystemExit(
        f"No --metrics passed and fallback {fallback} does not exist. "
        f"Run test_rfdiffusion first or pass --metrics explicitly."
    )


def _parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--metrics", default=None,
        help="Path to rfdiffusion_metrics.json from a test_rfdiffusion run.")
    parser.add_argument("--outdir", required=True,
        help="Directory to write plots into; created if missing.")
    parser.add_argument("--prod-script", default=None,
        help="Path to the production rfdiffusion_plots.py.  Used for the "
             "three already-good plots (specificity_coverage, design_lengths, "
             "com_displacement).  If omitted, those three plots are skipped.")
    return parser.parse_args()


def main():
    args = _parse_args()
    metrics_path = _resolve_metrics_path(args.metrics)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Reading metrics: {metrics_path}")
    data = json.loads(metrics_path.read_text())
    designs = data.get("designs", [])
    global_design_residues = set(data.get("design_residues", []))
    fixed_residues = set(data.get("fixed_residues", []))
    threshold = data.get("min_hotspot_frac", 0.0)
    print(f"  designs: {len(designs)}")

    if not designs:
        print("No designs in metrics — nothing to plot.")
        return 0

    # ── Render the three "good" plots via the production script ─────────
    if args.prod_script:
        prod_path = Path(args.prod_script)
        if not prod_path.exists():
            print(f"WARNING: --prod-script {prod_path} does not exist; "
                  f"skipping specificity_coverage plot.")
        else:
            print(f"Importing production plots from: {prod_path}")
            prod = _import_production_plots(prod_path)
            cwd_save = Path.cwd()
            os.chdir(outdir)
            try:
                # specificity_coverage looks good in production — reuse as-is.
                prod.plot_specificity_coverage(designs, global_design_residues,
                                                threshold, fixed_residues)
            finally:
                os.chdir(cwd_save)

    # ── Local design-lengths render (modified from production) ──────────
    print("Building design-lengths plot...")
    out = outdir / "rfdiff_design_lengths.png"
    try:
        ok = plot_design_lengths(designs, global_design_residues,
                                 fixed_residues, str(out))
        if ok:
            print(f"  → {out}")
        else:
            print(f"  (no plottable data)")
    except Exception as e:
        print(f"  FAILED: {type(e).__name__}: {e}")

    # ── Local com-displacement render (modified from production) ────────
    print("Building com-displacement plot...")
    out = outdir / "rfdiff_com_displacement.png"
    try:
        ok = plot_com_displacement(designs, str(out))
        if ok:
            print(f"  → {out}")
        else:
            print(f"  (no plottable data)")
    except Exception as e:
        print(f"  FAILED: {type(e).__name__}: {e}")

    # ── Calibrated contact-map render ───────────────────────────────────
    print("Building contact-map plot (calibrated G-style)...")
    cm_data = _build_contact_map_data(designs, global_design_residues, fixed_residues)
    eff_data = _build_effector_freq(designs, cm_data["n_designs"])
    max_n_res = max(cm_data["design_n_res"]) if cm_data["design_n_res"] else 0
    labels = _build_design_labels(cm_data["sorted_designs"], global_design_residues)
    if max_n_res == 0:
        print("  No residue data — skipping contact map.")
    else:
        out = outdir / "rfdiff_contact_map.png"
        try:
            plot_contact_map(cm_data, eff_data, max_n_res, labels,
                             global_design_residues, fixed_residues, str(out))
            print(f"  → {out}")
        except Exception as e:
            print(f"  FAILED: {type(e).__name__}: {e}")

    # ── Calibrated design-clustering render ─────────────────────────────
    print("Building design-clustering plot (calibrated G-style, square)...")
    out = outdir / "rfdiff_design_clustering.png"
    try:
        ok = plot_design_clustering(designs, global_design_residues, str(out))
        if ok:
            print(f"  → {out}")
        else:
            print(f"  (no plottable regions)")
    except Exception as e:
        print(f"  FAILED: {type(e).__name__}: {e}")

    print(f"\nDone.  Inspect: {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
