#!/usr/bin/env python3
"""
rfdiffusion_plots.py
--------------------
Diagnostic plots from an RFDiffusion backbone generation run.

    rfdiff_contact_map.png           - Per-design receptor contacts with
                                       per-design variable-length design regions,
                                       dendrogram-ordered, plus effector contact
                                       frequency bar below
    rfdiff_design_clustering.png     - Pairwise RMSD heatmap with dendrograms
    rfdiff_specificity_coverage.png  - Two-panel: specificity + coverage bars
    rfdiff_design_lengths.png        - Histogram of design-region lengths
    rfdiff_com_displacement.png      - Per-design centre-of-mass displacement
                                       of each design region from the input
                                       gap it replaces (Å), in the
                                       scaffold-aligned frame.
"""

import argparse
import json
import os
import sys

import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, BoundaryNorm
    from matplotlib.ticker import MaxNLocator, FormatStrFormatter
    import matplotlib.patches as mpatches
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

COLOUR_PASS     = "#1F77B4"   # blue  (colourblind-safe)
COLOUR_FAIL     = "#D62728"   # red
COLOUR_DESIGN   = "#FFD54F"
COLOUR_SUMMARY  = "#4C72B0"
COLOUR_COVERAGE = "#FF7F0E"

ALL_PLOT_FILES = [
    "rfdiff_contact_map.png",
    "rfdiff_design_clustering.png",
    "rfdiff_specificity_coverage.png",
    "rfdiff_design_lengths.png",
    "rfdiff_com_displacement.png",
]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metrics", required=True)
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


def _design_label(d):
    return d["design"].replace("design_", "").replace(".pdb", "")


def _get_per_design_residues(d, global_design_residues):
    """Get the design residues for a specific design, falling back to global."""
    per = d.get("per_design_design_residues")
    if per and len(per) > 0:
        return set(per)
    return global_design_residues


def _split_design_into_regions(per_design_dr, fixed_residues, per_region_lengths=None):
    """Split design resnums into contiguous groups (one per design region).

    If per_region_lengths is provided (from the metrics JSON), use it directly
    to partition the design residues — this avoids the sorting bug where
    synthetic resnums all land after fixed resnums.

    Falls back to the old sort-based heuristic only when per_region_lengths
    is not available (legacy metrics).

    Returns list of sets, one per contiguous design region.
    """
    if per_region_lengths and len(per_region_lengths) > 0:
        sorted_dr = sorted(per_design_dr)
        regions = []
        offset = 0
        for length in per_region_lengths:
            region = set(sorted_dr[offset:offset + length])
            if region:
                regions.append(region)
            offset += length
        if offset < len(sorted_dr):
            leftover = set(sorted_dr[offset:])
            if regions:
                regions[-1] |= leftover
            else:
                regions.append(leftover)
        return regions

    # Legacy fallback
    all_resnums = sorted(fixed_residues | per_design_dr)
    regions = []
    current_region = set()
    for r in all_resnums:
        if r in per_design_dr:
            current_region.add(r)
        else:
            if current_region:
                regions.append(current_region)
                current_region = set()
    if current_region:
        regions.append(current_region)
    return regions


def _get_dendrogram_order(designs):
    """Compute average-linkage clustering order from design-region coords."""
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

    from scipy.cluster.hierarchy import linkage, dendrogram
    from scipy.spatial.distance import squareform

    condensed = squareform(rmsd_matrix)
    Z = linkage(condensed, method="average")

    temp_fig, temp_ax = plt.subplots()
    dendro_ref = dendrogram(Z, ax=temp_ax, no_labels=True)
    order = dendro_ref["leaves"]
    plt.close(temp_fig)

    orig_order = [valid[k][0] for k in order]
    return orig_order, Z, rmsd_matrix


def _get_per_region_lengths(d):
    """Return per_region_lengths from metrics, or empty list."""
    prl = d.get("per_region_lengths")
    return prl if prl and len(prl) > 0 else []


def _label_with_region_len(d, region_idx, global_design_residues):
    """Build a 'd<n> (<len>)' label for a single region of one design."""
    lbl = _design_label(d)
    prl = _get_per_region_lengths(d)
    if prl and region_idx < len(prl):
        region_len = prl[region_idx]
    else:
        region_len = len(_get_per_design_residues(d, global_design_residues))
    return f"d{lbl} ({region_len})"


def _per_region_clustering_data(designs, region_idx):
    """Compute pairwise design-region Cα RMSD + linkage + dendrogram order
    for a single region.  Returns dict {valid, Z, rmsd_matrix, order, n}
    or None if the region has fewer than 2 designs with usable coords."""
    from scipy.cluster.hierarchy import linkage, dendrogram
    from scipy.spatial.distance import squareform

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
    """Compute per-design data needed for the contact map, in
    dendrogram order.  Returns a dict containing the sorted designs,
    per-design receptor-position arrays, contact assignments, and the
    linkage matrix Z (None if dendrogram couldn't be built)."""
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
        "sorted_designs":           sorted_designs,
        "n_designs":                n_designs,
        "design_n_res":             design_n_res,
        "design_region_positions":  design_region_positions,
        "design_contact_positions": design_contact_positions,
        "Z":                        Z,
    }


def _build_effector_freq(designs, n_designs):
    """Compute per-effector-residue contact frequencies across designs.
    Returns dict {range, freq, n_eff} or None if no effector contacts."""
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


def _draw_contact_heatmap_content(ax_heat, data, max_n_res):
    """Render the contact-map heatmap interior into ax_heat: hatched
    'beyond protein end' zones, design-region backgrounds, and contact
    blocks.  Sets x/y limits and inverts y-axis but does NOT touch
    tick labels."""
    n_designs = data["n_designs"]
    design_n_res = data["design_n_res"]
    design_region_positions = data["design_region_positions"]
    design_contact_positions = data["design_contact_positions"]

    # Beyond-protein-end hatching per design.  edgecolor must be set
    # (not 'none') for matplotlib to actually render the hatch lines;
    # linewidth=0 keeps the rectangle border itself invisible.
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
    # collide.  When too close, replace the previous-to-last so we
    # still show the endpoint without doubling up.
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
    """Render the effector-residue contact-frequency bar chart below
    the contact map."""
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
    """Return the four-patch legend for the contact-map plot."""
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
    """Draw clustering-heatmap content into ax_heat with optional
    dendrogram into ax_dendro.  When draw_colorbar=False the caller is
    responsible for adding a dedicated cbar axis (and gets the
    QuadMesh handle back via the return dict's 'im' key).  aspect can
    be 'auto' or 'equal' (square cells)."""
    from scipy.cluster.hierarchy import dendrogram

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
        "row_labels":     row_labels,
        "row_colours":    row_colours,
        "col_labels":     col_labels,
        "col_colours":    col_colours,
        "n":              n,
        "order":          order,
        "im":             im,
        "vmax":           vmax,
        "label_fontsize": label_fontsize,
        "cbar_fontsize":  cbar_fontsize,
    }


def _measure_label_width_inches(labels, fontsize, dpi=100):
    """Render the longest label to a throwaway figure and measure its
    width in inches.  Used to calibrate dendrogram-inset placement so
    the gap between labels and dendrogram is consistent regardless of
    label-text length."""
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


# ──────────────────────────────────────────────────────────────────────────
# Layout tunables — used by both contact-map and clustering plots
# ──────────────────────────────────────────────────────────────────────────
DENDRO_W_INCHES     = 1.6   # width of the dendrogram inset
LABEL_DENDRO_GAP_IN = 0.20  # visual gap between labels and dendrogram
RIGHT_MARGIN_IN     = 1.6   # space for cbar/legend on the right
TOP_MARGIN_IN       = 0.7
BOTTOM_MARGIN_IN    = 0.7


# ==============================================================================
# Plot 1: Per-design contact map with effector frequency bar below
# ==============================================================================

def plot_contact_map(designs, global_design_residues, fixed_residues):
    """
    Per-design receptor contact bar chart with dendrogram ordering.

    Layout strategy ("calibrated G"): heatmap drawn with normal
    left-side y-tick labels, dendrogram overlaid as an inset axis to
    the LEFT of those labels.  Inset position is calibrated to the
    measured width of the longest label, so the visual gap between
    labels and dendrogram is consistent regardless of label length.
    Figure dimensions are absolute (inches) and bbox_inches is not
    set on save, so layout is fully under our control.

    X-axis uses sequential position numbers (1..N) so the protein
    reads as a continuous chain.  Design region boundaries are marked
    with black lines.  Positions beyond a design's actual length are
    hatched (visible diagonal lines on a grey background).
    """
    out_path = "rfdiff_contact_map.png"

    if not designs:
        make_empty_plot("No designs to plot", out_path)
        return

    from scipy.cluster.hierarchy import dendrogram

    data = _build_contact_map_data(designs, global_design_residues, fixed_residues)
    if not data["design_n_res"]:
        make_empty_plot("No residue data in metrics", out_path)
        return

    eff_data = _build_effector_freq(designs, data["n_designs"])
    max_n_res = max(data["design_n_res"])
    labels = _build_design_labels(data["sorted_designs"], global_design_residues)

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
    print(f"Saved {out_path}")


# ==============================================================================
# Plot 2: Design Diversity Clustering
# ==============================================================================

def plot_design_clustering(designs, global_design_residues):
    """
    Multi-panel pairwise design-region Cα RMSD heatmaps with
    dendrograms (one panel per design region).

    Layout strategy: same as plot_contact_map but with per-panel
    dendrogram calibration — each panel's dendrogram is positioned
    LABEL_DENDRO_GAP_IN from THAT PANEL's labels' right edge, so a
    panel with shorter labels gets its dendrogram drifted rightwards,
    keeping the visual label-to-dendro gap consistent everywhere.
    Heatmaps stay aligned across panels (anchored on widest-panel
    label width) for easy cross-region comparison.

    Heatmaps are square (aspect='equal'), with a dedicated colorbar
    axis to the right and the Pass/Filtered legend above the heatmap
    right-aligned.
    """
    out_path = "rfdiff_design_clustering.png"

    if not designs:
        make_empty_plot("No designs to plot", out_path)
        return

    from scipy.cluster.hierarchy import dendrogram

    n_regions = _n_regions(designs)
    region_data = [_per_region_clustering_data(designs, r) for r in range(n_regions)]
    plottable = [(i, rd) for i, rd in enumerate(region_data) if rd is not None]
    if not plottable:
        make_empty_plot(
            "Not enough designs with usable region coords for clustering.",
            out_path,
        )
        return

    # First pass: per-panel label width and heatmap side length.
    panel_specs = []
    for r_idx, rd in plottable:
        n = rd["n"]
        heat_size = max(3, n * 0.5)
        label_fontsize = max(7, min(12, heat_size * 1.2))
        ordered = [rd["valid"][k][1] for k in rd["order"]]
        labels = [_label_with_region_len(d, r_idx, global_design_residues)
                  for d in ordered]
        label_w_in = _measure_label_width_inches(labels, label_fontsize)

        cell_in = 0.5
        heat_side_in = max(3.0, n * cell_in)

        panel_specs.append({
            "r_idx":          r_idx,
            "rd":             rd,
            "label_w_in":     label_w_in,
            "label_fontsize": label_fontsize,
            "heat_side_in":   heat_side_in,
        })

    # Heatmap left edge is fixed across panels (widest label width);
    # each panel's dendrogram floats horizontally to maintain its
    # own label-to-dendro gap.
    common_label_w_in = max(p["label_w_in"] for p in panel_specs)
    common_heat_side_in = max(p["heat_side_in"] for p in panel_specs)

    # Per-panel x-layout in inches.
    panel_left_margin_in  = 0.2
    cbar_gap_in           = 0.2
    cbar_w_in             = 0.25
    cbar_label_w_in       = 0.5    # rotated cbar label
    panel_right_margin_in = 0.2

    panel_w_in = (panel_left_margin_in + DENDRO_W_INCHES + LABEL_DENDRO_GAP_IN
                  + common_label_w_in + common_heat_side_in + cbar_gap_in
                  + cbar_w_in + cbar_label_w_in + panel_right_margin_in)

    # Vertical stacking — tight.
    xtick_room_in   = 0.4
    title_room_in   = 0.25
    legend_room_in  = 0.25
    panel_h_in      = (common_heat_side_in + xtick_room_in + title_room_in
                       + legend_room_in)
    panel_gap_in    = 0.2

    # Local margins (clustering plot is more compact than contact map).
    top_margin_in    = 0.3
    bottom_margin_in = 0.3

    fig_w_in = panel_w_in
    fig_h_in = (top_margin_in + panel_h_in * len(panel_specs)
                + panel_gap_in * (len(panel_specs) - 1) + bottom_margin_in)

    fig = plt.figure(figsize=(fig_w_in, fig_h_in))

    # Heatmap x-position is fixed across panels.
    heat_x_in = (panel_left_margin_in + DENDRO_W_INCHES + LABEL_DENDRO_GAP_IN
                 + common_label_w_in)
    cbar_x_in = heat_x_in + common_heat_side_in + cbar_gap_in

    for k, spec in enumerate(panel_specs):
        rd = spec["rd"]
        r_idx = spec["r_idx"]
        n = rd["n"]
        panel_label_w_in = spec["label_w_in"]

        panel_top_in = fig_h_in - top_margin_in - k * (panel_h_in + panel_gap_in)
        panel_y_in = panel_top_in - panel_h_in
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
    print(f"Saved {out_path}")


# ==============================================================================
# Plot 3: Specificity & Coverage (two side-by-side subplots)
# ==============================================================================

def plot_specificity_coverage(designs, global_design_residues, threshold, fixed_residues):
    out_path = "rfdiff_specificity_coverage.png"

    sorted_designs = sorted(designs, key=lambda d: d["frac_contacts_in_design"],
                            reverse=True)

    names = [_design_label(d) for d in sorted_designs]
    specificity = [d["frac_contacts_in_design"] for d in sorted_designs]
    passes = [d["passes_filter"] for d in sorted_designs]
    n = len(names)

    x = np.arange(n)
    bar_w = 0.7

    fig, (ax_spec, ax_cov) = plt.subplots(2, 1,
        figsize=(max(10, n * 0.8 + 3), 10), sharex=True)

    # ── Top: Specificity ─────────────────────────────────────────────────
    spec_colours = [COLOUR_PASS if p else COLOUR_FAIL for p in passes]
    ax_spec.bar(x, specificity, bar_w, color=spec_colours,
                edgecolor="white", linewidth=0.5)

    fontsize = max(5, 8 - n // 10)
    for i, d in enumerate(sorted_designs):
        ax_spec.text(x[i], specificity[i] + 0.02,
                     f"{d['n_contacts_in_design']}/{d['n_contact_pairs']}",
                     ha="center", va="bottom", fontsize=fontsize, color="black")

    if threshold > 0:
        ax_spec.axhline(threshold, color="black", linewidth=1.5,
                        linestyle="--", alpha=0.8, zorder=4)

    ax_spec.set_xticks(x)
    ax_spec.set_xticklabels(names, rotation=0, ha="center",
                            fontsize=max(6, 9 - n // 10))
    ax_spec.set_ylabel("Contact Specificity\n(Fraction)", fontsize=11)
    ax_spec.set_ylim(0, 1.15)

    spec_handles = [
        mpatches.Patch(color=COLOUR_PASS, label="Pass"),
        mpatches.Patch(color=COLOUR_FAIL, label="Filtered"),
    ]
    if threshold > 0:
        spec_handles.append(plt.Line2D([0], [0], color="black", linewidth=1.5,
                                        linestyle="--",
                                        label=f"Threshold ({threshold:.0%})"))
    ax_spec.legend(handles=spec_handles, fontsize=8, loc="upper right")

    # ── Bottom: Design-region coverage, split by region ──────────────────
    n_regions_per_design = []
    region_coverages = []
    for d in sorted_designs:
        d_design_res = _get_per_design_residues(d, global_design_residues)
        regions = _split_design_into_regions(d_design_res, fixed_residues,
                                             d.get("per_region_lengths"))
        n_regions_per_design.append(len(regions))
        contacted = set(d.get("receptor_contact_residues", []))
        coverages = []
        for region_set in regions:
            n_in_region = len(region_set)
            n_contacted = len(contacted & region_set)
            coverages.append(n_contacted / n_in_region if n_in_region > 0 else 0.0)
        region_coverages.append(coverages)

    max_n_regions = max(n_regions_per_design) if n_regions_per_design else 1
    region_colours = [COLOUR_COVERAGE, "#E377C2", "#17BECF", "#BCBD22"][:max_n_regions]

    group_width = 0.8
    sub_bar_w = group_width / max_n_regions
    for k in range(max_n_regions):
        offsets = x - group_width / 2 + sub_bar_w * (k + 0.5)
        vals = []
        for i in range(n):
            if k < len(region_coverages[i]):
                vals.append(region_coverages[i][k])
            else:
                vals.append(0.0)
        ax_cov.bar(offsets, vals, sub_bar_w * 0.9,
                   color=region_colours[k % len(region_colours)],
                   edgecolor="white", linewidth=0.3, alpha=0.85,
                   label=f"Region {k + 1}" if k < max_n_regions else None)

        for i in range(n):
            if k < len(region_coverages[i]):
                d_design_res = _get_per_design_residues(sorted_designs[i], global_design_residues)
                regions = _split_design_into_regions(d_design_res, fixed_residues,
                                                     sorted_designs[i].get("per_region_lengths"))
                if k < len(regions):
                    n_reg = len(regions[k])
                    n_cont = round(region_coverages[i][k] * n_reg)
                    ax_cov.text(offsets[i], vals[i] + 0.02,
                                f"{n_cont}/{n_reg}", ha="center", va="bottom",
                                fontsize=max(4, fontsize - 1), color="black")

    ax_cov.set_xticks(x)
    ax_cov.set_xticklabels(names, rotation=0, ha="center",
                           fontsize=max(6, 9 - n // 10))
    ax_cov.set_xlabel("Design", fontsize=11)
    ax_cov.set_ylabel("Design-Region Coverage\n(Fraction)", fontsize=11)
    ax_cov.set_ylim(0, 1.15)
    if max_n_regions > 1:
        ax_cov.legend(fontsize=8, loc="upper right")

    plt.subplots_adjust(hspace=0.15)
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved {out_path}")


# ==============================================================================
# Plot 4: Design-region length distribution
# ==============================================================================

def plot_design_lengths(designs, global_design_residues, fixed_residues):
    """
    Multi-panel bar chart of design-region lengths, one panel per
    contiguous design region.  Discrete integer bars, not histograms.

    Each panel's width scales with the number of unique length values
    so that bars render at a consistent physical width across regions
    (a region with one unique length gets a narrow panel rather than
    one massive stretched bar).  X-ticks include every integer between
    the min and max unique value (so gap-values get a tick + label too)
    but never extend past the data range.
    """
    out_path = "rfdiff_design_lengths.png"

    if not designs:
        make_empty_plot("No designs to plot", out_path)
        return

    all_region_lengths = []
    designs_missing_prl = 0
    for d in designs:
        prl = d.get("per_region_lengths")
        if prl and len(prl) > 0:
            for k, length in enumerate(prl):
                while len(all_region_lengths) <= k:
                    all_region_lengths.append([])
                all_region_lengths[k].append(length)
        else:
            # Refuse to fall back to the legacy sort-based heuristic for
            # size reporting — synthetic de novo resnums (filter line ~306)
            # all land after fixed resnums, so the heuristic would collapse
            # all regions into one and produce actively misleading sizes.
            designs_missing_prl += 1

    if designs_missing_prl > 0:
        print(f"WARNING: {designs_missing_prl}/{len(designs)} designs have no "
              f"per_region_lengths in metrics; their region sizes are not "
              f"plotted.  Re-run rfdiffusion_filter.py to regenerate metrics.",
              file=sys.stderr)
        if not all_region_lengths:
            make_empty_plot(
                "No per_region_lengths in metrics.\n"
                "Re-run rfdiffusion_filter.py to regenerate.",
                out_path,
            )
            return

    n_regions = len(all_region_lengths)
    if n_regions == 0:
        make_empty_plot("No design regions found", out_path)
        return

    # Per-region pre-computation so we can size each panel by its
    # unique-value count.
    panel_specs = []
    for k in range(n_regions):
        lengths = np.array(all_region_lengths[k])
        unique_lengths, counts = np.unique(lengths, return_counts=True)
        lo, hi = int(unique_lengths.min()), int(unique_lengths.max())
        all_ints = list(range(lo, hi + 1))
        panel_specs.append({
            "lengths":  lengths,
            "unique":   unique_lengths,
            "counts":   counts,
            "all_ints": all_ints,
            "n_bars":   len(all_ints),
        })

    # Layout constants — bars render at target_bar_w_in physical width
    # regardless of how many bars share the panel.  The min_panel_w_in
    # floor stops the panel collapsing to a sliver when there are few
    # unique values; when the floor kicks in, target_data_range scales
    # up too so the bar stays the target width.
    target_bar_w_in       = 0.6
    bar_width_data        = 0.8
    edge_pad_data         = 0.4
    min_panel_w_in        = 3.0
    panel_padding_in      = 1.2
    panel_h_in            = 3.0
    panel_gap_in          = 0.4
    bottom_margin_in      = 0.5
    top_margin_in         = 0.3

    panel_widths_in = []
    for s in panel_specs:
        span_data = s["unique"][-1] - s["unique"][0]
        target_data_range = span_data + bar_width_data + 2 * edge_pad_data
        natural = target_data_range * target_bar_w_in / bar_width_data
        panel_widths_in.append(max(min_panel_w_in, natural) + panel_padding_in)
    fig_w_in = max(panel_widths_in)
    fig_h_in = panel_h_in * n_regions + panel_gap_in * (n_regions - 1) + 1.0

    fig = plt.figure(figsize=(fig_w_in, fig_h_in))

    for k in range(n_regions):
        spec = panel_specs[k]
        panel_top_in = fig_h_in - top_margin_in - k * (panel_h_in + panel_gap_in)
        panel_y_in = panel_top_in - panel_h_in
        ax_y_frac = (panel_y_in + 0.5) / fig_h_in
        ax_h_frac = (panel_h_in - 0.8) / fig_h_in

        # Data-coord range that exactly fits all bars + edge padding.
        span_data = spec["unique"][-1] - spec["unique"][0]
        target_data_range = span_data + bar_width_data + 2 * edge_pad_data
        natural_ax_w_in = target_data_range * target_bar_w_in / bar_width_data
        ax_w_in = max(min_panel_w_in, natural_ax_w_in)
        # If the min-width floor inflated the axis above natural, scale
        # target_data_range up to match — extra empty space around the
        # bar but bar physical width stays constant across panels.
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
        center = (lo + hi) / 2.0
        ax.set_xlim(center - target_data_range / 2.0,
                    center + target_data_range / 2.0)

        ax.set_xlabel("Design-region length (residues)", fontsize=10)
        ax.set_ylabel("Count", fontsize=10)
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        ax.set_title(f"Design region {k + 1}", fontsize=11)

    plt.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Saved {out_path}")


# ==============================================================================
# Plot 5: Per-design COM displacement of each design region from the input gap
# ==============================================================================

def plot_com_displacement(designs):
    """
    Per-design centre-of-mass displacement of each design region from the
    input gap residues it replaced, in the scaffold-aligned frame (Å).

    For multi-region contigs, each design region gets its own grouped bar
    within the design's bar group.  Designs are ordered left-to-right by
    the maximum displacement across their regions, so the most-changed
    designs sit at the right of the plot.

    A value of 0 means RFDiffusion built the new region centred at exactly
    the same place in space as the original gap residues sat; larger
    values mean the new region's centre of mass is somewhere different.
    Length-independent — works the same way for fixed-length and
    variable-length contigs.

    Designs missing the per-region field, or with all-NaN region values,
    are skipped from the plot (and reported via the in-plot annotation).
    """
    out_path = "rfdiff_com_displacement.png"

    if not designs:
        make_empty_plot("No designs to plot", out_path)
        return

    # Collect (design_name, [region_disp, ...]) pairs for designs that
    # have any non-null displacement values.  None and NaN both count as
    # "missing"; treat them uniformly here.
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
        make_empty_plot(
            "No design_region_com_displacement values in metrics.\n"
            "Re-run rfdiffusion_filter.py to regenerate.",
            out_path,
        )
        return

    # Sort designs by maximum (non-None) region displacement, ascending
    # so the rightmost bars are the most-changed designs.
    def _max_disp(item):
        vals = [v for v in item[1] if v is not None]
        return max(vals) if vals else -1.0
    plotted.sort(key=_max_disp)

    n_designs = len(plotted)
    n_regions = max(len(p[1]) for p in plotted)

    fig, ax = plt.subplots(figsize=(max(8, n_designs * 0.4 + 2), 5))

    # Group bars: one slot per design, n_regions sub-bars per slot.
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
    # Verbose y-axis label carries the meaning so the title can be omitted.
    ax.set_ylabel("Design region COM displacement\nfrom input gap, scaffold-aligned (Å)",
                  fontsize=10)
    # No title — the y-axis label says it all.
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

    if not os.path.exists(args.metrics) or os.path.getsize(args.metrics) == 0:
        save_fallback_plots("No metrics data")
        sys.exit(0)

    with open(args.metrics) as f:
        data = json.load(f)

    designs = data.get("designs", [])
    global_design_residues = set(data.get("design_residues", []))
    fixed_residues = set(data.get("fixed_residues", []))
    threshold = data.get("min_hotspot_frac", 0.0)

    if not designs:
        save_fallback_plots("No designs in metrics")
        sys.exit(0)

    plot_contact_map(designs, global_design_residues, fixed_residues)
    plot_design_clustering(designs, global_design_residues)
    plot_specificity_coverage(designs, global_design_residues, threshold, fixed_residues)
    plot_design_lengths(designs, global_design_residues, fixed_residues)
    plot_com_displacement(designs)


if __name__ == "__main__":
    main()