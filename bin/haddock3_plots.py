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
import numpy as np

from haddock_utils import (
    parse_capri_tsv, parse_clustfcc_tsv, get_numeric_col,
    model_stem, extract_heavy_atoms, contacted_residues,
    collect_pdb_index,
)

# ── Styling ──────────────────────────────────────────────────────────────────
COLOUR_ALL      = "#4C72B0"
COLOUR_BEST     = "#DD4444"
COLOUR_WARN     = "#FFCCCC"
COLOUR_DENOVO       = "#FF8C00"
COLOUR_EFF_ACTIVE   = "#8B5CF6"
COLOUR_CLASH        = "#DC2626"
BSA_WARN_CUTOFF = 700.0
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
    parser.add_argument("--min-cluster-size", type=int, default=MIN_CLUSTER_SIZE,
                        help="Minimum cluster size to qualify (default: %(default)s). "
                             "Should match `min_population` in haddock.nf clustfcc block.")
    parser.add_argument("--effector-active-residues", default="",
                        help="Comma-separated effector active residues for AIR annotation bar")
    parser.add_argument("--contact-pairs", default="",
                        help="Space-separated CA-CA pin pairs, e.g. 'A73-B31 A72-B32'. "
                             "Each pair's receptor and effector residues are marked with "
                             "a triangle on the heatmap annotation bars.")
    parser.add_argument("--restraints-summary", default=None,
                        help="Path to restraints_summary.json from haddock3_prepare.py. "
                             "If provided, the 'contig_design_region' field is used to "
                             "shade the receptor design region on the heatmap.")
    parser.add_argument("--cluster-metrics", default=None,
                        help="Path to cluster_metrics.json from haddock_cluster_metrics.py. "
                             "Required for the cluster-overview / ranking / pair-satisfaction "
                             "/ clash-breakdown plots.  When absent those plots are skipped.")
    parser.add_argument("--cluster-sc", default=None,
                        help="Path to cluster_sc.json from haddock_cluster_sc.py. "
                             "Optional — used as a fallback when cluster_metrics.json has "
                             "sc=null for some clusters.")
    return parser.parse_args()


def _parse_residue_spec(spec):
    """Parse "25,40-44" into a set of residue numbers.  Empty -> empty set."""
    if not spec or not spec.strip():
        return set()
    out = set()
    for token in spec.split(","):
        token = token.strip()
        if not token:
            continue
        if "-" in token:
            lo, hi = token.split("-", 1)
            try:
                out.update(range(int(lo), int(hi) + 1))
            except ValueError:
                pass
        else:
            try:
                out.add(int(token))
            except ValueError:
                pass
    return out


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
                  "haddock_interface_heatmap.png",
                  "haddock_cluster_overview.png",
                  "haddock_cluster_ranking.png",
                  "haddock_pair_satisfaction.png",
                  "haddock_clash_breakdown.png"]:
        make_empty_plot(message, name)


# ═══════════════════════════════════════════════════════════════════════════════
# Plot 1: Score vs BSA
# ═══════════════════════════════════════════════════════════════════════════════

def plot_score_vs_bsa(scores, bsas, best_score, best_bsa,
                      cluster_ids=None, best_cluster_id=None):
    """HADDOCK score vs BSA scatter, coloured by cluster.

    With ``cluster_ids`` supplied (one per data point), points are
    coloured by their cluster_id using a discrete tab10 palette; the
    star marker still highlights the lowest-score model.  Without
    cluster_ids, falls back to the old single-colour behaviour.
    """
    fig, ax = plt.subplots(figsize=(7.5, 5))

    if cluster_ids is not None and len(cluster_ids) == len(scores):
        # Discrete colour per cluster.  Sort the distinct IDs so cluster 1
        # always gets the first colour (consistent across plots).  '-' is
        # treated as 'unclustered' and drawn in grey.
        distinct = sorted({c for c in cluster_ids if c != "-"},
                          key=lambda x: (str(x).isdigit() and int(x), str(x)))
        cmap = plt.get_cmap("tab10")
        colour_by_cid = {cid: cmap(i % 10) for i, cid in enumerate(distinct)}
        colour_by_cid["-"] = "#999999"

        # Draw cluster by cluster so each gets a separate legend entry.
        for cid in distinct + (["-"] if any(c == "-" for c in cluster_ids) else []):
            xs = [b for b, c in zip(bsas, cluster_ids) if c == cid]
            ys = [s for s, c in zip(scores, cluster_ids) if c == cid]
            if not xs:
                continue
            label = f"Cluster {cid} (n={len(xs)})" if cid != "-" else f"Unclustered (n={len(xs)})"
            ax.scatter(xs, ys, color=colour_by_cid[cid], alpha=0.8, s=45,
                       edgecolors="white", linewidths=0.6, label=label, zorder=3)
    else:
        ax.scatter(bsas, scores, color=COLOUR_ALL, alpha=0.6, s=30,
                   edgecolors="white", linewidths=0.5,
                   label="Top-N per cluster (seletopclusts output)", zorder=3)

    # Star: lowest-score model overall.  When best_cluster_id is given,
    # outline the star in the cluster's colour so the link is visible.
    star_edge = "black"
    if cluster_ids is not None and best_cluster_id is not None and best_cluster_id != "-":
        # Recompute palette to find the star's edge colour.
        distinct = sorted({c for c in cluster_ids if c != "-"},
                          key=lambda x: (str(x).isdigit() and int(x), str(x)))
        cmap = plt.get_cmap("tab10")
        for i, cid in enumerate(distinct):
            if cid == best_cluster_id:
                star_edge = cmap(i % 10)
                break
    ax.scatter(best_bsa, best_score, color="#FFD700", s=180, zorder=5,
               edgecolors=star_edge, linewidths=1.5,
               label=f"Best model (BSA={best_bsa:.0f} Å²)", marker="*")

    saved_xlim = ax.get_xlim()
    ax.axvspan(saved_xlim[0], BSA_WARN_CUTOFF, color=COLOUR_WARN, alpha=0.5,
               label=f"BSA < {BSA_WARN_CUTOFF:.0f} Å² (weak interface)", zorder=1)
    ax.set_xlim(saved_xlim)

    ax.set_xlabel("Buried Surface Area (Å²)", fontsize=12)
    ax.set_ylabel("HADDOCK Score", fontsize=12)

    y_min, y_max = ax.get_ylim()
    ax.set_ylim(top=y_max + (y_max - y_min) * 0.25)
    ax.legend(fontsize=9, loc="upper right", framealpha=0.9)

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
    # Y axis is a count → integer ticks only.
    from matplotlib.ticker import MaxNLocator
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))

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
# Plot 3: Cluster metrics overview (heatmap-table)
# ═══════════════════════════════════════════════════════════════════════════════


# Metric columns shown in the cluster overview + their colour-mapping
# rules.  Each tuple: (column_header, json_key, colourmap_name,
# higher_is_better).  When higher_is_better is False the colourmap is
# inverted so red = bad in both directions.
_OVERVIEW_COLUMNS = [
    ("BSA (Å²)",         "bsa",                            "Blues",   True),
    ("Sc",               "sc",                             "Blues",   True),
    ("Pair contacts",    "pair_contact_fraction",          "Greens",  True),
    ("AIR sat.",         "air_sat_frac",                   "Greens",  True),
    ("COM (Å)",          "com_distance",                   "Purples", False),
    ("H-bonds",          "interface_hbonds",               "Blues",   True),
    ("Clashes in DR",    "clashes_in_design_region",       "Oranges", False),
    ("Clashes out DR",   "clashes_outside_design_region",  "Reds",    False),
    ("Members",          "size",                           "Greys",   True),
]


def _norm01(values):
    """Min-max normalise a list of numeric values to [0,1].  Returns
    a list of normalised floats; None values stay as None.  All-equal
    inputs return 0.5 for each (avoids divide-by-zero)."""
    nums = [v for v in values if v is not None]
    if not nums:
        return [None] * len(values)
    lo, hi = min(nums), max(nums)
    if hi - lo < 1e-9:
        return [None if v is None else 0.5 for v in values]
    return [None if v is None else (v - lo) / (hi - lo) for v in values]


def plot_cluster_overview(metrics_by_cid, sc_by_cid, sorted_cids,
                          selected_cid, cluster_sizes_by_cid,
                          out_path="haddock_cluster_overview.png"):
    """Heatmap-style table: rows = qualifying clusters (sorted by auto-pick
    key), cols = metrics.  Cell shading = relative value within the column;
    cell text = actual value.  Selected cluster row is marked with a star.
    """
    if not sorted_cids:
        make_empty_plot("No qualifying clusters", out_path)
        print(f"Saved {out_path} (no qualifying clusters)")
        return

    # Pull values column by column (so colour scaling is per-column).
    rows = []
    for cid in sorted_cids:
        m = dict(metrics_by_cid.get(cid, {}))
        # Compose derived columns.
        sat = m.get("air_satisfaction_count")
        total = m.get("air_total_count")
        m["air_sat_frac"] = (sat / total) if (sat is not None and total) else None
        # Pull sc from sidecar if missing here.
        if m.get("sc") is None and cid in sc_by_cid:
            m["sc"] = sc_by_cid[cid].get("sc")
        m["size"] = cluster_sizes_by_cid.get(cid)
        rows.append(m)

    # Per-column value array.
    col_values = []
    for _, key, _, _ in _OVERVIEW_COLUMNS:
        col_values.append([r.get(key) for r in rows])

    n_rows = len(rows)
    n_cols = len(_OVERVIEW_COLUMNS)
    fig_w = max(8.5, n_cols * 1.05 + 2.0)
    fig_h = max(2.2, n_rows * 0.55 + 1.4)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    from matplotlib.patches import Rectangle
    ax.set_xlim(-0.5, n_cols - 0.5)
    ax.set_ylim(n_rows - 0.5, -0.5)
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels([h for (h, _, _, _) in _OVERVIEW_COLUMNS],
                       rotation=30, ha="right", fontsize=9)
    row_labels = []
    for cid in sorted_cids:
        star = "★ " if cid == selected_cid else "  "
        row_labels.append(f"{star}Cluster {cid}")
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(row_labels, fontsize=10)

    for j, (header, key, cmap_name, _higher_better) in enumerate(_OVERVIEW_COLUMNS):
        # Shade by magnitude only — deeper colour = more of the metric.
        # The column header / colourmap choice tells the reader whether
        # "more of this" is good (Blues/Greens) or bad (Oranges/Reds).
        # This avoids the brain-bender of "low BSA is bad so its colour
        # should be deep red" — instead, low BSA = light blue, and the
        # reader sees the absolute amount of BSA in the cell.
        values = col_values[j]
        norm = _norm01(values)
        cmap = plt.get_cmap(cmap_name)
        for i, (val, nv) in enumerate(zip(values, norm)):
            if val is None:
                facecolor = "#EEEEEE"
                text = "—"
            else:
                shade = 0.15 + 0.65 * (nv if nv is not None else 0.5)
                facecolor = cmap(shade)
                if isinstance(val, float):
                    text = f"{val:.2f}" if abs(val) < 1000 else f"{val:.0f}"
                else:
                    text = f"{val}"
            ax.add_patch(Rectangle(
                (j - 0.5, i - 0.5), 1.0, 1.0,
                facecolor=facecolor, edgecolor="white", linewidth=1.0,
            ))
            ax.text(j, i, text, ha="center", va="center",
                    fontsize=9, color="#222222")

    # Highlight selected row with a thick outer border.
    if selected_cid in sorted_cids:
        sel_idx = sorted_cids.index(selected_cid)
        ax.add_patch(Rectangle(
            (-0.5, sel_idx - 0.5), n_cols, 1.0,
            facecolor="none", edgecolor="#222222", linewidth=2.2,
            zorder=5,
        ))

    ax.set_title(
        "Per-cluster metrics overview "
        "(★ = auto-pick / chosen cluster; deeper colour = larger value "
        "in that column — read the column header for whether high or low is good)",
        fontsize=10, pad=12,
    )
    for spine in ("top", "right", "bottom", "left"):
        ax.spines[spine].set_visible(False)
    ax.tick_params(left=False, bottom=False)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved {out_path}")


# ═══════════════════════════════════════════════════════════════════════════════
# Plot 4: Cluster ranking scatter
# ═══════════════════════════════════════════════════════════════════════════════


def plot_cluster_ranking(metrics_by_cid, sc_by_cid, sorted_cids, selected_cid,
                         cluster_sizes_by_cid,
                         out_path="haddock_cluster_ranking.png"):
    """2D ranking scatter: x = BSA, y = pair contact fraction.

    Each cluster is a dot.  Dot size = cluster member count.  Dot colour
    = Sc (diverging Blues — higher = more complementary).  Selected
    cluster is a star.  Top-right = best on the auto-pick axes.
    """
    if not sorted_cids:
        make_empty_plot("No qualifying clusters", out_path)
        print(f"Saved {out_path} (no qualifying clusters)")
        return

    xs, ys, sizes, scs, labels = [], [], [], [], []
    for cid in sorted_cids:
        m = metrics_by_cid.get(cid, {}) or {}
        bsa = m.get("bsa")
        pf = m.get("pair_contact_fraction")
        sc = m.get("sc")
        if sc is None and cid in sc_by_cid:
            sc = sc_by_cid[cid].get("sc")
        if bsa is None or pf is None:
            continue
        xs.append(bsa)
        ys.append(pf)
        sizes.append(120 + 30 * cluster_sizes_by_cid.get(cid, 1))
        scs.append(sc if sc is not None else float("nan"))
        labels.append(str(cid))

    fig, ax = plt.subplots(figsize=(7.5, 5.5))

    # Sc colourbar — drop NaNs so vmin/vmax are well-defined.
    valid_sc = [s for s in scs if s == s]  # NaN != NaN trick
    if valid_sc:
        vmin = min(valid_sc)
        vmax = max(valid_sc)
        if vmax - vmin < 1e-9:
            vmin, vmax = max(0.0, vmin - 0.05), vmax + 0.05
        sc_cmap = plt.get_cmap("viridis")
    else:
        vmin, vmax, sc_cmap = 0.0, 1.0, plt.get_cmap("viridis")

    for cid, x, y, sz, sc in zip(labels, xs, ys, sizes, scs):
        is_selected = (cid == selected_cid)
        color = sc_cmap((sc - vmin) / (vmax - vmin)) if sc == sc else "#BBBBBB"
        marker = "*" if is_selected else "o"
        edge = "#222222" if is_selected else "white"
        lw = 1.5 if is_selected else 0.6
        adjusted_size = sz * 1.8 if is_selected else sz
        ax.scatter(x, y, s=adjusted_size, c=[color], marker=marker,
                   edgecolors=edge, linewidths=lw, zorder=5,
                   label=f"Cluster {cid}" + (" ★" if is_selected else ""))
        ax.text(x, y - 0.04, f"  {cid}", fontsize=9, ha="left", va="top",
                color="#333333")

    ax.set_xlabel("Buried Surface Area (Å²)", fontsize=12)
    ax.set_ylabel("Pair contact fraction (user pins satisfied)", fontsize=12)
    ax.set_ylim(-0.05, 1.10)
    ax.axhline(1.0, color="#888888", linestyle=":", linewidth=1.0, alpha=0.6,
               zorder=1)
    ax.axvline(BSA_WARN_CUTOFF, color=COLOUR_WARN, linestyle=":",
               linewidth=1.0, alpha=0.7, zorder=1)
    ax.set_title(
        "Cluster ranking — top-right = best on the auto-pick "
        "(pair_contact_fraction first, BSA tie-breaker).  "
        "Dot size = cluster members; colour = Sc.",
        fontsize=9, pad=10,
    )
    # Sc colourbar.
    if valid_sc:
        sm = plt.cm.ScalarMappable(
            cmap=sc_cmap,
            norm=plt.Normalize(vmin=vmin, vmax=vmax),
        )
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, pad=0.02, fraction=0.04)
        cbar.set_label("Sc (Lawrence–Colman shape complementarity)", fontsize=9)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved {out_path}")


# ═══════════════════════════════════════════════════════════════════════════════
# Plot 5: Pair contact satisfaction per cluster
# ═══════════════════════════════════════════════════════════════════════════════


def _pair_ca_distance(pdb_path, rec_chain, rec_resnum, eff_chain, eff_resnum):
    """CA-CA distance between two residues on the same PDB.  Returns
    None if either CA is missing.  Transparently handles gzipped PDBs.
    """
    import gzip as _gzip
    import numpy as _np
    rec_ca = eff_ca = None
    opener = _gzip.open if str(pdb_path).endswith(".gz") else open
    with opener(pdb_path, "rt") as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            if line[12:16].strip() != "CA":
                continue
            ch = line[21]
            try:
                rn = int(line[22:26].strip())
            except ValueError:
                continue
            if ch == rec_chain and rn == rec_resnum:
                rec_ca = _np.array([float(line[30:38]), float(line[38:46]),
                                    float(line[46:54])])
            elif ch == eff_chain and rn == eff_resnum:
                eff_ca = _np.array([float(line[30:38]), float(line[38:46]),
                                    float(line[46:54])])
            if rec_ca is not None and eff_ca is not None:
                break
    if rec_ca is None or eff_ca is None:
        return None
    return float(_np.linalg.norm(rec_ca - eff_ca))


def plot_pair_satisfaction(contact_pairs, pair_distance, pdb_index,
                           clusters, sorted_cids, selected_cid,
                           rec_chain, eff_chain,
                           out_path="haddock_pair_satisfaction.png"):
    """Per cluster, per user contact pair: scatter the actual CA-CA
    distance on each cluster's best model PDB.  Vertical band marks the
    restraint window [target-lo_dev, target+hi_dev] in the pair_distance
    triple.  Cluster row is starred if it's the auto-pick.
    """
    if not contact_pairs:
        make_empty_plot(
            "No contact pairs specified.\n"
            "(pass --contact-pairs to enable this plot)",
            out_path,
        )
        print(f"Saved {out_path} (no contact pairs)")
        return
    if not sorted_cids:
        make_empty_plot("No qualifying clusters", out_path)
        print(f"Saved {out_path} (no clusters)")
        return

    target, lo_dev, hi_dev = pair_distance
    lower = max(0.0, target - lo_dev)
    upper = target + hi_dev

    n_pairs = len(contact_pairs)
    n_clusters = len(sorted_cids)
    fig, ax = plt.subplots(figsize=(max(7, n_pairs * 1.0 + 4),
                                    max(2.5, n_clusters * 0.55 + 1.8)))

    # Compute distance per (cluster, pair).
    distances = {}  # (cid, pair_idx) -> float or None
    for cid in sorted_cids:
        members = clusters.get(cid, [])
        if not members:
            continue
        scored = [(m, s) for (m, s) in members if s is not None]
        best_member = (min(scored, key=lambda x: x[1])[0]
                       if scored else members[0][0])
        stem = model_stem(best_member)
        pdb = pdb_index.get(stem)
        if pdb is None:
            continue
        for pi, (rc, rn, ec, en) in enumerate(contact_pairs):
            d = _pair_ca_distance(pdb, rc, rn, ec, en)
            distances[(cid, pi)] = d

    # Layout: one ROW per cluster, one X position per pair.  Markers
    # coloured by satisfaction status.
    ax.axvspan(lower, upper, color="#A7F3D0", alpha=0.5, zorder=1,
               label=f"Restraint window ({lower:.1f}–{upper:.1f} Å)")
    ax.axvline(target, color="#10B981", linestyle="--", linewidth=1.0,
               alpha=0.7, zorder=2, label=f"Target ({target:.1f} Å)")

    pair_marker_offset = {}
    for pi, (rc, rn, ec, en) in enumerate(contact_pairs):
        pair_marker_offset[pi] = pi  # used in legend, just for label
    # Stagger label vertical offsets by pair index so they don't pile up
    # when multiple pairs land at similar distances.
    n_pairs_local = max(1, len(contact_pairs))
    label_offsets = [-0.18 - 0.10 * pi for pi in range(n_pairs_local)]
    for i, cid in enumerate(sorted_cids):
        for pi, (rc, rn, ec, en) in enumerate(contact_pairs):
            d = distances.get((cid, pi))
            if d is None:
                continue
            satisfied = d <= upper
            color = "#22C55E" if satisfied else "#DC2626"
            ax.scatter(d, i, s=80, color=color, edgecolors="black",
                       linewidths=0.6, zorder=4)
            # Small connector line from marker to label so the pairing
            # is obvious even when the label is staggered upward.
            label_y = i + label_offsets[pi]
            ax.plot([d, d], [i, label_y + 0.04], color="#888888",
                    linewidth=0.5, alpha=0.6, zorder=3)
            ax.text(d, label_y, f"{rc}{rn}-{ec}{en}",
                    ha="center", va="bottom", fontsize=7,
                    color="#444444")

    ax.set_xlabel("CA–CA distance on cluster best model (Å)", fontsize=11)
    row_labels = []
    for cid in sorted_cids:
        star = "★ " if cid == selected_cid else "  "
        row_labels.append(f"{star}Cluster {cid}")
    ax.set_yticks(range(n_clusters))
    ax.set_yticklabels(row_labels, fontsize=10)
    ax.set_ylim(n_clusters - 0.5, -0.7)

    # Custom legend.
    from matplotlib.lines import Line2D
    legend_elems = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#22C55E",
               markeredgecolor="black", markersize=9, label="Satisfied"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#DC2626",
               markeredgecolor="black", markersize=9, label="Outside window"),
    ]
    ax.legend(handles=legend_elems, loc="lower right", fontsize=9, framealpha=0.9)

    ax.set_title(
        f"Pair contact satisfaction — green band = restraint window "
        f"[{lower:.1f}, {upper:.1f}] Å; markers = CA-CA distance per pair "
        f"on each cluster's best model",
        fontsize=9, pad=10,
    )
    ax.grid(True, axis="x", linestyle=":", alpha=0.3, zorder=0)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved {out_path}")


# ═══════════════════════════════════════════════════════════════════════════════
# Plot 6: Clash bookkeeping breakdown
# ═══════════════════════════════════════════════════════════════════════════════


def plot_clash_breakdown(metrics_by_cid, sorted_cids, selected_cid,
                         out_path="haddock_clash_breakdown.png"):
    """Per cluster, stacked horizontal bar:
       green segment = clashes inside design region (tolerated — these
           residues will be redesigned by RFDiffusion);
       red segment   = clashes outside design region (real geometric
           problems that survive into the RFDiffusion input).

    Selected cluster row is starred.  Bars sorted by auto-pick order.
    """
    if not sorted_cids:
        make_empty_plot("No qualifying clusters", out_path)
        print(f"Saved {out_path} (no clusters)")
        return

    in_dr  = [metrics_by_cid.get(c, {}).get("clashes_in_design_region", 0)
              for c in sorted_cids]
    out_dr = [metrics_by_cid.get(c, {}).get("clashes_outside_design_region", 0)
              for c in sorted_cids]

    n_clusters = len(sorted_cids)
    fig, ax = plt.subplots(figsize=(8, max(2.5, n_clusters * 0.55 + 1.5)))

    y_positions = list(range(n_clusters))
    bar_in  = ax.barh(y_positions, in_dr, color="#86EFAC",
                      edgecolor="#15803D", linewidth=0.8,
                      label="Clashes inside design region (tolerated)",
                      zorder=3)
    bar_out = ax.barh(y_positions, out_dr, left=in_dr, color="#FCA5A5",
                      edgecolor="#B91C1C", linewidth=0.8,
                      label="Clashes outside design region (real problem)",
                      zorder=3)

    # Numeric annotations on each segment.
    for i, (a, b) in enumerate(zip(in_dr, out_dr)):
        if a > 0:
            ax.text(a / 2, i, f"{a}", ha="center", va="center",
                    fontsize=9, color="#14532D", fontweight="bold", zorder=4)
        if b > 0:
            ax.text(a + b / 2, i, f"{b}", ha="center", va="center",
                    fontsize=9, color="#7F1D1D", fontweight="bold", zorder=4)
        if a + b == 0:
            ax.text(0.5, i, "no clashes", ha="left", va="center",
                    fontsize=8, color="#777777", style="italic", zorder=4)

    row_labels = []
    for cid in sorted_cids:
        star = "★ " if cid == selected_cid else "  "
        row_labels.append(f"{star}Cluster {cid}")
    ax.set_yticks(y_positions)
    ax.set_yticklabels(row_labels, fontsize=10)
    ax.set_ylim(n_clusters - 0.5, -0.7)
    ax.set_xlabel("Heavy-atom clash count (< 2.0 Å cross-chain)", fontsize=11)
    from matplotlib.ticker import MaxNLocator
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    if max(a + b for a, b in zip(in_dr, out_dr)) == 0:
        ax.set_xlim(0, 1)
    ax.legend(loc="lower right", fontsize=9, framealpha=0.9)
    ax.set_title(
        "Clash bookkeeping — clashes inside the design region get "
        "redesigned away by RFDiffusion; clashes outside survive into the "
        "RFDiffusion input.",
        fontsize=9, pad=10,
    )
    ax.grid(True, axis="x", linestyle=":", alpha=0.3, zorder=0)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved {out_path}")


# ═══════════════════════════════════════════════════════════════════════════════
# Plot 7: Interface contact heatmap
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


def _per_residue_clashes(pdb_path, rec_chain, eff_chain, cutoff=2.0):
    """Return (set_of_clashing_receptor_resnums, set_of_clashing_effector_resnums)
    for heavy-atom pairs under ``cutoff`` Å across the two chains.

    Used by the heatmap to draw red borders around residues with at
    least one cross-chain heavy-atom clash.  Mirrors the bookkeeping
    that ``structure_metrics.clash_count`` does but exposes the per-
    residue identity rather than just the count.
    """
    import numpy as _np
    import gzip as _gzip
    rec_atoms = []  # list of (resnum, [x,y,z])
    eff_atoms = []
    opener = _gzip.open if str(pdb_path).endswith(".gz") else open
    with opener(pdb_path, "rt") as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            ch = line[21]
            atom_name = line[12:16].strip()
            element = line[76:78].strip()
            if element == "H" or (not element and atom_name.startswith("H")):
                continue
            try:
                resnum = int(line[22:26].strip())
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            if ch == rec_chain:
                rec_atoms.append((resnum, x, y, z))
            elif ch == eff_chain:
                eff_atoms.append((resnum, x, y, z))
    if not rec_atoms or not eff_atoms:
        return set(), set()
    eff_coords = _np.array([(x, y, z) for _, x, y, z in eff_atoms])
    eff_res    = _np.array([r for r, _, _, _ in eff_atoms])
    cutoff_sq = cutoff * cutoff
    rec_clash = set()
    eff_clash = set()
    for rn, x, y, z in rec_atoms:
        diffs = eff_coords - _np.array([x, y, z])
        sq = _np.einsum("ij,ij->i", diffs, diffs)
        hits = sq < cutoff_sq
        if hits.any():
            rec_clash.add(rn)
            eff_clash.update(int(r) for r in eff_res[hits])
    return rec_clash, eff_clash


def _contiguous_ranges(residues):
    """Group a set/list of residue numbers into contiguous (start, end) ranges."""
    if not residues:
        return []
    sorted_r = sorted(set(residues))
    ranges = []
    run_start = sorted_r[0]
    prev = run_start
    for r in sorted_r[1:]:
        if r != prev + 1:
            ranges.append((run_start, prev))
            run_start = r
        prev = r
    ranges.append((run_start, prev))
    return ranges


def _design_ranges(design_residues, all_pdb_residues):
    """Restrict ``design_residues`` to those actually present in the PDB and
    return contiguous (start, end) ranges for heatmap shading.

    Post-Session 7: ``design_residues`` is the HADDOCK design region
    (= receptor_active_residues + receptor halves of contact_pairs).
    Previously this function complemented a fixed-residues set; both
    framings shade the same residues, just inverted in their input
    semantics.
    """
    if not design_residues or not all_pdb_residues:
        return []
    return _contiguous_ranges(set(design_residues) & set(all_pdb_residues))


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


def _overlay_clash_borders(ax, res_to_col, row_idx, clash_residues):
    """Draw a red border around each (clash residue, row) cell.

    The cell fill colour is NOT touched — only the rectangle outline is
    drawn.  zorder=5 sits above the heatmap imshow (zorder=0/1) but
    below the design-region / active-residue dashed lines (zorder=10),
    matching the user's requirement that dashed lines remain on top.
    """
    from matplotlib.patches import Rectangle
    for r in clash_residues:
        col = res_to_col.get(r)
        if col is None:
            continue
        rect = Rectangle(
            (col - 0.5, row_idx - 0.5), 1.0, 1.0,
            facecolor="none", edgecolor=COLOUR_CLASH,
            linewidth=1.5, zorder=5, clip_on=True,
        )
        ax.add_patch(rect)


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
                           effector_active_residues=None,
                           contact_pairs=None):
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
    best_rec_clashes  = best_eff_clashes  = set()
    best_model_pdb_path = None
    if best_model_name:
        stem = model_stem(best_model_name)
        fpath = pdb_index.get(stem)
        if fpath:
            both = extract_heavy_atoms(fpath, (rec_chain, eff_chain))
            r_atoms, e_atoms = both[rec_chain], both[eff_chain]
            if r_atoms and e_atoms:
                best_rec_contacts = contacted_residues(r_atoms, e_atoms, cutoff)
                best_eff_contacts = contacted_residues(e_atoms, r_atoms, cutoff)
            best_model_pdb_path = fpath
            best_rec_clashes, best_eff_clashes = _per_residue_clashes(
                fpath, rec_chain, eff_chain, cutoff=2.0,
            )

    has_best_row = best_rec_contacts is not None

    # ── Per-cluster clash sets (one set per cluster's best model) ───────
    # Used for red borders on the heatmap.  Cluster best model = lowest-
    # score member that we have a PDB for.
    rec_clashes_by_cid = {}
    eff_clashes_by_cid = {}
    for cid in sorted_cids:
        members = clusters.get(cid, [])
        if not members:
            continue
        scored = [(m, s) for (m, s) in members if s is not None]
        if scored:
            best_member_name, _ = min(scored, key=lambda x: x[1])
        else:
            best_member_name = members[0][0]
        cstem = model_stem(best_member_name)
        cpath = pdb_index.get(cstem)
        if cpath is None:
            continue
        r_c, e_c = _per_residue_clashes(cpath, rec_chain, eff_chain, cutoff=2.0)
        rec_clashes_by_cid[cid] = r_c
        eff_clashes_by_cid[cid] = e_c

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

    # Design regions = user-supplied design residues intersected with PDB residues
    denovo_ranges = _design_ranges(fixed_residues, set(rec_res_range)) \
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

    # Receptor contact-pair markers: small black triangles above the bar
    # pointing down to the receptor halves of each user-defined pair.
    rec_pair_residues = []
    if contact_pairs:
        rec_pair_residues = [rn for (_rc, rn, _ec, _en) in contact_pairs
                             if rn in rec_res_to_col]
        for rn in rec_pair_residues:
            ax_bar.scatter(rec_res_to_col[rn], 1.15, marker="v",
                           color="#000000", s=40, clip_on=False, zorder=5)
    if rec_pair_residues:
        ax_bar.text(0.5, 1.45,
                    f"▼ user contact-pair anchors ({len(rec_pair_residues)})",
                    ha="center", va="bottom", transform=ax_bar.transAxes,
                    fontsize=8, color="#333333")

    # ── Shared colourmap ─────────────────────────────────────────────────
    cmap = plt.cm.Blues.copy()
    cmap.set_bad(color="#EEEEEE")

    # ── Receptor heatmap ─────────────────────────────────────────────────
    im_rec = _render_heatmap_panel(ax_rec, rec_matrix, rec_res_range, n_rec_rows, cmap, has_best_row)

    # Clash borders on receptor heatmap — drawn BEFORE the dashed
    # design-region lines so the dashes overlay any red box edges.
    if has_best_row and best_rec_clashes:
        _overlay_clash_borders(ax_rec, rec_res_to_col, 0, best_rec_clashes)
    for row_idx, cid in enumerate(sorted_cids):
        if cid in rec_clashes_by_cid:
            _overlay_clash_borders(
                ax_rec, rec_res_to_col, row_idx + rec_row_offset,
                rec_clashes_by_cid[cid],
            )

    if has_denovo:
        for dn_start, dn_end in denovo_ranges:
            cols = [rec_res_to_col[r] for r in range(dn_start, dn_end + 1) if r in rec_res_to_col]
            if cols:
                for xv in (min(cols) - 0.5, max(cols) + 0.5):
                    ax_rec.axvline(xv, color=COLOUR_DENOVO, linewidth=1.2,
                                   linestyle="--", alpha=0.95, zorder=10)

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

        # Effector contact-pair markers (mirror of the receptor markers).
        eff_pair_residues = []
        if contact_pairs:
            eff_pair_residues = [en for (_rc, _rn, _ec, en) in contact_pairs
                                 if en in eff_res_to_col]
            for en in eff_pair_residues:
                ax_eff_bar.scatter(eff_res_to_col[en], 1.15, marker="v",
                                   color="#000000", s=40, clip_on=False, zorder=5)
        if eff_pair_residues:
            ax_eff_bar.text(0.5, 1.45,
                            f"▼ user contact-pair anchors "
                            f"({len(eff_pair_residues)})",
                            ha="center", va="bottom",
                            transform=ax_eff_bar.transAxes,
                            fontsize=8, color="#333333")

    # ── Effector heatmap ─────────────────────────────────────────────────
    if has_eff and ax_eff is not None:
        _render_heatmap_panel(ax_eff, eff_matrix, eff_res_range, n_eff_rows, cmap, has_best_row)
        _set_heatmap_yticks(ax_eff, sorted_cids, mean_scores, eff_found, eff_total,
                            clusters, eff_row_offset, has_best_row, best_model_name)
        ax_eff.set_ylabel("FCC Cluster", fontsize=11)
        _set_heatmap_xticks(ax_eff, eff_res_range, "Effector residue")
        ax_eff.set_xlim(-0.5, n_eff_res - 0.5)

        # Clash borders on effector heatmap (best-model row + per-cluster).
        eff_res_to_col_map = {r: i for i, r in enumerate(eff_res_range)}
        if has_best_row and best_eff_clashes:
            _overlay_clash_borders(ax_eff, eff_res_to_col_map, 0, best_eff_clashes)
        for row_idx, cid in enumerate(sorted_cids):
            if cid in eff_clashes_by_cid:
                _overlay_clash_borders(
                    ax_eff, eff_res_to_col_map, row_idx + eff_row_offset,
                    eff_clashes_by_cid[cid],
                )

        # Dashed lines at the boundaries of effector active-residue runs —
        # mirrors the receptor design-region dashed lines above.
        if effector_active_residues:
            eff_active_set = set(effector_active_residues)
            in_pdb = sorted(r for r in eff_active_set if r in eff_res_to_col)
            if in_pdb:
                runs = []
                run_start = prev = in_pdb[0]
                for r in in_pdb[1:]:
                    if r != prev + 1:
                        runs.append((run_start, prev))
                        run_start = r
                    prev = r
                runs.append((run_start, prev))
                for r_start, r_end in runs:
                    cols = [eff_res_to_col[r] for r in range(r_start, r_end + 1)
                            if r in eff_res_to_col]
                    if cols:
                        for xv in (min(cols) - 0.5, max(cols) + 0.5):
                            ax_eff.axvline(xv, color=COLOUR_EFF_ACTIVE,
                                           linewidth=1.2, linestyle="--",
                                           alpha=0.95, zorder=10)

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
    # cluster_id is a string column ('1', '2', '-', ...) — keep raw.
    cluster_ids_per_model = None
    for name in ("cluster_id", "cluster", "cluster-id"):
        if data and name in data[0]:
            cluster_ids_per_model = [row[name].strip() or "-" for row in data]
            break

    if scores is None:
        save_fallback_plots("No HADDOCK score column found")
        sys.exit(0)

    best_score = scores[0]
    best_cluster_id = (cluster_ids_per_model[0]
                       if cluster_ids_per_model else None)

    # Parse contact pairs (used by both score-vs-BSA's downstream and by
    # the per-cluster plots / heatmap markers).  Single source of truth.
    contact_pairs = []
    if args.contact_pairs:
        import re as _re
        for token in args.contact_pairs.split():
            m = _re.match(r"^([A-Za-z])(\d+)-([A-Za-z])(\d+)$", token.strip())
            if m:
                contact_pairs.append((
                    m.group(1).upper(), int(m.group(2)),
                    m.group(3).upper(), int(m.group(4)),
                ))
    if contact_pairs:
        print(f"Contact pairs parsed: {len(contact_pairs)}")

    # ── Plot 1: Score vs BSA ─────────────────────────────────────────────
    if bsas is not None:
        plot_score_vs_bsa(scores, bsas, best_score, bsas[0],
                          cluster_ids=cluster_ids_per_model,
                          best_cluster_id=best_cluster_id)
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

    # ── Plots 3-6 require cluster_metrics.json from haddock_cluster_metrics.py.
    metrics_by_cid: dict = {}
    sc_by_cid: dict = {}
    if args.cluster_metrics and os.path.exists(args.cluster_metrics):
        import json as _json
        try:
            with open(args.cluster_metrics) as f:
                raw = _json.load(f)
            metrics_by_cid = {str(k): v for k, v in raw.items()}
        except (ValueError, OSError) as e:
            print(f"WARNING: could not read --cluster-metrics: {e}",
                  file=sys.stderr)
    if args.cluster_sc and os.path.exists(args.cluster_sc):
        import json as _json
        try:
            with open(args.cluster_sc) as f:
                raw = _json.load(f)
            sc_by_cid = {str(k): v for k, v in raw.items()}
        except (ValueError, OSError) as e:
            print(f"WARNING: could not read --cluster-sc: {e}",
                  file=sys.stderr)

    # Cluster sort order — match SELECT_HADDOCK_CLUSTER's auto-pick key.
    overview_sorted_cids: list = []
    selected_cid: str = None
    if metrics_by_cid:
        named_ids = [c for c in clusters.keys() if c != "-"]
        qualifying_for_overview = [c for c in named_ids
                                   if str(c) in metrics_by_cid]
        def _auto_pick_key(cid):
            m = metrics_by_cid[str(cid)]
            pf = m.get("pair_contact_fraction", 1.0)
            bsa = m.get("bsa", 0.0)
            return (-pf, -bsa)
        overview_sorted_cids = sorted(qualifying_for_overview, key=_auto_pick_key)
        if overview_sorted_cids:
            selected_cid = overview_sorted_cids[0]

    sizes_by_cid = {cid: len(members) for cid, members in clusters.items()}

    # ── Plot 3: Cluster overview heatmap-table ───────────────────────────
    if metrics_by_cid:
        plot_cluster_overview(metrics_by_cid, sc_by_cid, overview_sorted_cids,
                              selected_cid, sizes_by_cid)
    else:
        make_empty_plot(
            "cluster_metrics.json not provided\n"
            "(pass --cluster-metrics to enable this plot)",
            "haddock_cluster_overview.png",
        )

    # ── Plot 4: Cluster ranking scatter ──────────────────────────────────
    if metrics_by_cid:
        plot_cluster_ranking(metrics_by_cid, sc_by_cid, overview_sorted_cids,
                             selected_cid, sizes_by_cid)
    else:
        make_empty_plot(
            "cluster_metrics.json not provided",
            "haddock_cluster_ranking.png",
        )

    # ── Plot 5: Pair contact satisfaction ───────────────────────────────
    # Pull pair_distance from restraints_summary.json if present;
    # otherwise default to (2, 2, 4).
    pair_distance = (2.0, 2.0, 4.0)
    if args.restraints_summary and os.path.exists(args.restraints_summary):
        try:
            import json as _json
            with open(args.restraints_summary) as f:
                _rs = _json.load(f)
            pd_field = _rs.get("pair_distance")
            if isinstance(pd_field, list) and len(pd_field) == 3:
                pair_distance = tuple(float(x) for x in pd_field)
        except (ValueError, OSError):
            pass

    # The contact-pairs structure is the parsed tuple list from earlier.
    if contact_pairs and overview_sorted_cids:
        # Build a pdb_index that may not have been built yet if we're
        # ahead of the heatmap path — defer to the existing pdb_index
        # construction.  We need it here too.
        if 'pdb_index' not in locals() or not pdb_index:
            pdb_index = collect_pdb_index(complex_dir=args.complex_dir,
                                          run_dir=args.run_dir)
        plot_pair_satisfaction(contact_pairs, pair_distance, pdb_index,
                               clusters, overview_sorted_cids, selected_cid,
                               args.receptor_chain, args.effector_chain)
    elif not contact_pairs:
        make_empty_plot(
            "No contact pairs specified", "haddock_pair_satisfaction.png",
        )
    else:
        make_empty_plot(
            "No qualifying clusters", "haddock_pair_satisfaction.png",
        )

    # ── Plot 6: Clash breakdown ─────────────────────────────────────────
    if metrics_by_cid:
        plot_clash_breakdown(metrics_by_cid, overview_sorted_cids, selected_cid)
    else:
        make_empty_plot(
            "cluster_metrics.json not provided",
            "haddock_clash_breakdown.png",
        )

    # ── Plot 7: Interface heatmap ────────────────────────────────────────
    pdb_index = collect_pdb_index(complex_dir=args.complex_dir, run_dir=args.run_dir)
    if not pdb_index:
        src = args.run_dir or args.complex_dir
        if src:
            print(f"WARNING: No docked PDB files found under '{src}'")
        else:
            print("WARNING: --run-dir or --complex-dir not provided; heatmap will be empty.")

    # Design-region shading: read contig_design_region from
    # restraints_summary.json (post-commit-4 amendment).  The contig is
    # the authoritative source for the receptor design region; HADDOCK
    # PREPARE pre-computes it and stores the residue list there.
    design_residues: set = set()
    if args.restraints_summary and os.path.exists(args.restraints_summary):
        try:
            import json as _json
            with open(args.restraints_summary) as f:
                _rs = _json.load(f)
            design_residues = set(int(r) for r in
                                  _rs.get("contig_design_region", []))
        except (ValueError, OSError) as e:
            print(f"WARNING: could not read --restraints-summary "
                  f"{args.restraints_summary}: {e}", file=sys.stderr)
    fixed_residues = design_residues if design_residues else None
    if fixed_residues:
        print(f"Design-region residues for shading "
              f"(from contig_design_region): {len(fixed_residues)}")
    else:
        print("No design region available; heatmap design-region shading disabled.")

    # (contact_pairs parsed earlier in main() — see above.)

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
                           effector_active_residues=sorted(set(eff_active)) if eff_active else None,
                           contact_pairs=contact_pairs)


if __name__ == "__main__":
    main()