#!/usr/bin/env python3
"""
pose_solver_plots.py
--------------------
Render the five diagnostic plots from a completed POSE_SOLVE run:

  1. pose_loss_breakdown.png       — per-term loss contributions at the
                                     final pose (validity, dist_upper,
                                     dist_lower, excl, interp,
                                     pair_sc_clash) with weight
                                     constants annotated.
  2. pose_pair_distances.png       — one whisker row per pair constraint
                                     showing [min, max] window + achieved
                                     distance.
  3. pose_interface_dashboard.png  — 4-panel interface metrics summary
                                     (BSA, contact-distance histogram,
                                     gap-index, residue composition).
  4. pose_restart_loss_curve.png   — per-restart loss + best-so-far
                                     trajectory; convergence diagnostic.
  5. pose_contig_comparison.png    — three contig options (gap-0, gap-1,
                                     gap-2) shown as residue bars so the
                                     user can pick the right join_gap.

Sized for sparse data: 2-4 pairs, single solved pose.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Dict, List

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np


# ── Common style ─────────────────────────────────────────────────────────────

PALETTE = {
    "validity":      "#7F0000",  # very serious penalty
    "dist_upper":    "#DC2626",
    "dist_lower":    "#F97316",
    "excl":          "#A855F7",
    "interp":        "#0EA5E9",
    "pair_sc_clash": "#22C55E",
    "bsa_a":         "#2563EB",
    "bsa_b":         "#0891B2",
    "bsa_total":     "#1E40AF",
    "satisfied":     "#16A34A",
    "unsatisfied":   "#DC2626",
    "fixed":         "#E5E7EB",
    "design":        "#FB923C",
    "window":        "#BFDBFE",
    "reference":     "#9CA3AF",
    "hydrophobic":   "#FBBF24",
    "polar":         "#34D399",
    "charged":       "#F87171",
    "aromatic":      "#A78BFA",
}


# ── Plot 1: loss breakdown ──────────────────────────────────────────────────

def plot_loss_breakdown(results: Dict, out_path: Path) -> None:
    """Six horizontal bars (one per penalty term) with the weight
    constant annotated at the end of each bar.  Drops the 'total' key.
    """
    breakdown = results.get("loss_breakdown") or {}
    weights = results.get("loss_weights") or {}
    weight_key = {
        "validity": "W_VALIDITY", "dist_upper": "W_DIST", "dist_lower": "W_LOWER",
        "excl": "W_EXCL", "interp": "W_INTERP", "pair_sc_clash": "W_PAIR_SC_CLASH",
    }
    terms = [k for k in
             ("validity", "dist_upper", "dist_lower",
              "excl", "interp", "pair_sc_clash")
             if k in breakdown]
    if not terms:
        _empty(out_path, "No loss breakdown available.")
        return

    values = [breakdown[k] for k in terms]
    colours = [PALETTE.get(k, "#888888") for k in terms]

    fig, ax = plt.subplots(figsize=(6.5, 3.4))
    y = np.arange(len(terms))
    ax.barh(y, values, color=colours, edgecolor="black", linewidth=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels(terms, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("Weighted loss contribution", fontsize=9)
    ax.set_title(f"Loss breakdown at final pose  "
                 f"(total = {breakdown.get('total', sum(values)):.2f})",
                 fontsize=10)
    ax.grid(axis="x", alpha=0.3, linestyle=":")

    # Force a non-zero x-axis so 0-bars are still visible as ticks.
    if max(values) <= 0:
        ax.set_xlim(0, 1)
    else:
        ax.set_xlim(0, max(values) * 1.25)

    # Annotate value + weight on each bar.
    for yi, term, v in zip(y, terms, values):
        wk = weight_key.get(term, "?")
        wv = weights.get(wk, "?")
        wv_str = f"{wv:.0e}" if isinstance(wv, (int, float)) else "?"
        txt = f"{v:.2f}  ({wk}={wv_str})"
        ax.text(v + max(values, default=1) * 0.01, yi, txt,
                fontsize=7.5, va="center", ha="left", color="black")

    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ── Plot 2: pair-constraint distance whiskers ───────────────────────────────

def plot_pair_distances(results: Dict, binder_chain: str, target_chain: str,
                        out_path: Path) -> None:
    pairs = results.get("pair_results") or []
    if not pairs:
        _empty(out_path, "No pair results to plot.")
        return

    # Tight rows but readable.
    fig_h = max(2.0, 0.6 * len(pairs) + 1.2)
    fig, ax = plt.subplots(figsize=(6.5, fig_h))

    # Fixed x-axis 0-10 Å so the window isn't auto-zoomed away.
    ax.set_xlim(0, 10)
    ax.set_xlabel("CA–CA distance (Å)", fontsize=9)

    labels = []
    for i, pr in enumerate(pairs):
        lo, hi = pr["min_d"], pr["max_d"]
        achieved = pr["achieved_d"]
        satisfied = pr.get("satisfied", lo <= achieved <= hi)

        # Shaded allowed window
        ax.barh(i, hi - lo, left=lo, height=0.55,
                color=PALETTE["window"], edgecolor="#3B82F6", linewidth=0.7)
        # Achieved marker
        ax.scatter(achieved, i, marker="D", s=80,
                   color=PALETTE["satisfied"] if satisfied else PALETTE["unsatisfied"],
                   edgecolors="black", linewidths=0.6, zorder=5)
        labels.append(f"{binder_chain}{pr['binder_res']} – {target_chain}{pr['target_res']}")

    ax.set_yticks(np.arange(len(pairs)))
    ax.set_yticklabels(labels, fontsize=9)
    ax.invert_yaxis()
    ax.grid(axis="x", alpha=0.3, linestyle=":")

    # Reference line at typical β-strand CA–CA (~5 Å) — annotate via
    # the x-axis transform so the label sits at the bottom of the axes
    # regardless of how many rows are stacked above it.
    ax.axvline(5.0, color=PALETTE["reference"], linestyle="--",
               linewidth=0.8, alpha=0.7)
    ax.text(5.05, 0.02, "β-strand ~5 Å",
            transform=ax.get_xaxis_transform(),
            fontsize=7, color=PALETTE["reference"], va="bottom")

    # Legend
    legend = [
        mpatches.Patch(facecolor=PALETTE["window"], edgecolor="#3B82F6",
                       label="Allowed window"),
        plt.Line2D([0], [0], marker="D", color="w",
                   markerfacecolor=PALETTE["satisfied"], markeredgecolor="black",
                   markersize=8, label="Satisfied"),
        plt.Line2D([0], [0], marker="D", color="w",
                   markerfacecolor=PALETTE["unsatisfied"], markeredgecolor="black",
                   markersize=8, label="Unsatisfied"),
    ]
    ax.legend(handles=legend, fontsize=7, framealpha=0.9, loc="lower right")
    ax.set_title("Pair-constraint distances at final pose", fontsize=10)

    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ── Plot 3: interface metrics dashboard ─────────────────────────────────────

def plot_interface_dashboard(metrics_list: List[Dict],
                             binder_chain: str, target_chain: str,
                             out_path: Path) -> None:
    """Single-pose 4-panel summary.  metrics_list is biopython-shaped:
    a list with one entry per PDB (pose_interface_metrics.py emits a
    list).  We render the first entry.
    """
    if not metrics_list:
        _empty(out_path, "No interface metrics to plot.")
        return
    m = metrics_list[0]

    fig, axes = plt.subplots(2, 2, figsize=(10, 7.5))
    fig.suptitle("Interface metrics — solved pose", fontsize=11, y=1.00)

    # ── (a) BSA bars ────────────────────────────────────────────────
    ax = axes[0, 0]
    bsa = m.get("bsa_A2", {})
    bars = [
        (f"Chain {binder_chain}", bsa.get("chain_a", 0.0), PALETTE["bsa_a"]),
        (f"Chain {target_chain}", bsa.get("chain_b", 0.0), PALETTE["bsa_b"]),
        ("Total",                  bsa.get("total",   0.0), PALETTE["bsa_total"]),
    ]
    x = np.arange(len(bars))
    vals = [b[1] for b in bars]
    ax.bar(x, vals, color=[b[2] for b in bars], edgecolor="black", linewidth=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels([b[0] for b in bars], fontsize=9)
    ax.set_ylabel("BSA (Å²)", fontsize=9)
    ax.set_title("(a) Buried surface area", fontsize=9.5)
    ax.grid(axis="y", alpha=0.3, linestyle=":")
    # Annotate values above each bar.
    for xi, v in zip(x, vals):
        ax.text(xi, v + max(vals, default=1) * 0.02, f"{v:.0f}",
                ha="center", va="bottom", fontsize=8.5)
    # Reference dashed line — typical native interface BSA (~2700 Å² for the
    # validated AvrPikF pose, per pipeline_notes16 §3.1).
    ax.axhline(2700, color=PALETTE["reference"], linestyle="--",
               linewidth=0.8, alpha=0.7)
    ax.text(len(bars) - 0.4, 2700, "  native ≈2700 Å²",
            fontsize=7, color=PALETTE["reference"], va="bottom")

    # ── (b) Contact-distance histogram ──────────────────────────────
    ax = axes[0, 1]
    hist = (m.get("heavy_contacts_5A", {}).get("distance_histogram")
            or {})
    if hist:
        bins = list(hist.keys())
        counts = list(hist.values())
        # Colour bars by clash severity.
        cols = []
        for b in bins:
            lo = float(b.split("-")[0])
            if lo < 2.0:
                cols.append("#7F0000")     # clash zone
            elif lo < 2.5:
                cols.append("#DC2626")
            elif lo < 3.0:
                cols.append("#F97316")
            elif lo < 3.5:
                cols.append("#FACC15")
            else:
                cols.append("#16A34A")     # comfortable
        ax.bar(range(len(bins)), counts, color=cols,
               edgecolor="black", linewidth=0.4)
        ax.set_xticks(range(len(bins)))
        ax.set_xticklabels(bins, rotation=45, ha="right", fontsize=7.5)
        ax.set_xlabel("Heavy-atom distance bin (Å)", fontsize=9)
        ax.set_ylabel("Atom-pair count", fontsize=9)
    else:
        ax.text(0.5, 0.5, "No histogram", ha="center", va="center",
                transform=ax.transAxes)
    ax.set_title("(b) Inter-chain heavy-atom distances", fontsize=9.5)
    ax.grid(axis="y", alpha=0.3, linestyle=":")

    # ── (c) Gap-index ───────────────────────────────────────────────
    ax = axes[1, 0]
    gi = m.get("shape_complementarity_proxy", {}).get("gap_index_mean_A")
    gi = float(gi) if gi is not None and not np.isnan(gi) else None
    # Reference: AF3-native PikF gap-index = 3.41 Å (pipeline_notes16 §3.1).
    native = 3.41
    if gi is not None:
        ax.barh([0], [gi], height=0.35,
                color="#0EA5E9", edgecolor="black", linewidth=0.5)
        ax.set_xlim(0, max(5.0, gi * 1.2, native * 1.2))
        ax.set_ylim(-0.6, 0.6)
        ax.text(gi + 0.05, 0, f"{gi:.2f} Å",
                fontsize=9, va="center", ha="left")
    else:
        ax.text(0.5, 0.5, "No gap-index", ha="center", va="center",
                transform=ax.transAxes)
        ax.set_xlim(0, 5.0)
        ax.set_ylim(-0.6, 0.6)
    ax.axvline(native, color=PALETTE["reference"], linestyle="--",
               linewidth=0.8, alpha=0.7)
    ax.text(native + 0.05, -0.4, f"native ≈{native} Å",
            fontsize=7, color=PALETTE["reference"], va="top")
    ax.set_yticks([0])
    ax.set_yticklabels(["this pose"], fontsize=9)
    ax.set_xlabel("Mean nearest-other-chain heavy-atom distance (Å)",
                  fontsize=8.5)
    ax.set_title("(c) Gap-index (lower = tighter packing)", fontsize=9.5)
    ax.grid(axis="x", alpha=0.3, linestyle=":")

    # ── (d) Interface composition ──────────────────────────────────
    ax = axes[1, 1]
    comp_a = m.get("interface_residues", {}).get("composition_a") or {}
    comp_b = m.get("interface_residues", {}).get("composition_b") or {}
    categories = ["hydrophobic", "polar", "charged", "aromatic"]
    width = 0.35
    xpos = np.arange(len(categories))
    a_vals = [comp_a.get(c, 0.0) for c in categories]
    b_vals = [comp_b.get(c, 0.0) for c in categories]
    ax.bar(xpos - width/2, a_vals, width, label=f"Chain {binder_chain}",
           color=PALETTE["bsa_a"], edgecolor="black", linewidth=0.4)
    ax.bar(xpos + width/2, b_vals, width, label=f"Chain {target_chain}",
           color=PALETTE["bsa_b"], edgecolor="black", linewidth=0.4)
    ax.set_xticks(xpos)
    ax.set_xticklabels(categories, fontsize=9)
    ax.set_ylim(0, max(1.0, max(a_vals + b_vals, default=0.0) * 1.2))
    ax.set_ylabel("Fraction of interface residues", fontsize=9)
    ax.set_title("(d) Interface residue composition", fontsize=9.5)
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(axis="y", alpha=0.3, linestyle=":")

    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ── Plot 4: restart loss convergence curve ──────────────────────────────────

def plot_restart_loss_curve(csv_path: Path, out_path: Path) -> None:
    """Scatter of per-restart loss + step plot of best-so-far.  Highlights
    the restart at which the global best was achieved.
    """
    restarts: List[int] = []
    losses: List[float] = []
    bests: List[float] = []
    with open(csv_path) as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            restarts.append(int(row["restart_idx"]))
            losses.append(float(row["restart_loss"]))
            bests.append(float(row["best_loss_so_far"]))

    if not restarts:
        _empty(out_path, "No restart history.")
        return

    fig, ax = plt.subplots(figsize=(8, 3.8))
    if len(restarts) == 1:
        # Single point (DE path) — just show the value.
        ax.scatter(restarts, losses, s=80, color="#2563EB",
                   edgecolors="black", linewidths=0.5, label="DE final loss")
        ax.set_xlim(-0.5, 0.5)
    else:
        ax.scatter(restarts, losses, s=12, color="#94A3B8", alpha=0.5,
                   label="Per-restart loss")
        ax.step(restarts, bests, where="post", color="#DC2626", linewidth=1.5,
                label="Best-so-far")
        # Mark the first restart that hit the global minimum.
        best_global = min(bests)
        first_hit = next(i for i, b in enumerate(bests) if b == best_global)
        ax.axvline(first_hit, color="#16A34A", linestyle="--", linewidth=0.8,
                   alpha=0.8)
        ax.text(first_hit + max(1, len(restarts) * 0.01),
                best_global * 1.05 if best_global > 0 else 1,
                f"converged at restart {first_hit}\n(loss = {best_global:.3f})",
                fontsize=8, color="#16A34A", va="bottom")

    # Log y-axis: validity penalties can dwarf everything else.
    if min(losses) > 0:
        ax.set_yscale("log")
    ax.set_xlabel("Restart index", fontsize=9)
    ax.set_ylabel("Loss", fontsize=9)
    ax.set_title("Restart-loss convergence", fontsize=10)
    ax.grid(alpha=0.3, linestyle=":")
    ax.legend(fontsize=8, loc="upper right")
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ── Plot 5: contig comparison ───────────────────────────────────────────────

def plot_contig_comparison(results: Dict, posed_pdb: Path,
                           binder_chain: str, out_path: Path) -> None:
    """Three stacked rows showing the three contig options (no joining,
    gap ≤ 1, gap ≤ 2) as residue bars over the binder sequence.  Fixed
    residues grey, design-region residues orange.
    """
    contigs = results.get("contigs") or {}
    if not contigs:
        _empty(out_path, "No contig options to compare.")
        return

    # Get the full binder residue list from the posed PDB.
    binder_resnums = sorted({
        int(line[22:26])
        for line in posed_pdb.read_text().splitlines()
        if line.startswith("ATOM") and line[21] == binder_chain
    })
    if not binder_resnums:
        _empty(out_path, f"No chain {binder_chain} residues in {posed_pdb.name}.")
        return

    # Parse the design region for each label.  The contig key stores a
    # `design_region` field like "33-49, 69-78".
    rows = []
    for label, info in contigs.items():
        dr_str = info.get("design_region", "")
        dr_set: set[int] = set()
        for tok in (dr_str or "").split(","):
            tok = tok.strip()
            if not tok or tok == "(none)":
                continue
            if "-" in tok:
                lo, hi = tok.split("-", 1)
                dr_set.update(range(int(lo), int(hi) + 1))
            else:
                dr_set.add(int(tok))
        rows.append((label, dr_set, info.get("n_residues", 0)))

    fig, ax = plt.subplots(figsize=(8, max(2.5, 0.55 * len(rows) + 1.2)))
    n_res = len(binder_resnums)
    # Draw each row as a series of unit-width coloured cells.
    for i, (label, dr, n_design) in enumerate(rows):
        for j, rn in enumerate(binder_resnums):
            colour = PALETTE["design"] if rn in dr else PALETTE["fixed"]
            ax.add_patch(plt.Rectangle((j, i - 0.4), 1, 0.8,
                                       facecolor=colour, edgecolor="white",
                                       linewidth=0.2))
        ax.text(n_res + 0.5, i, f"{n_design} de novo",
                fontsize=8, va="center", ha="left", color="#4B5563")

    ax.set_xlim(0, n_res + 8)
    ax.set_ylim(-0.6, len(rows) - 0.4)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r[0] for r in rows], fontsize=9)
    ax.invert_yaxis()

    # X-ticks: every 10th residue.
    step = max(1, n_res // 12)
    ticks = list(range(0, n_res, step))
    if ticks[-1] != n_res - 1:
        ticks.append(n_res - 1)
    ax.set_xticks([t + 0.5 for t in ticks])
    ax.set_xticklabels([str(binder_resnums[t]) for t in ticks], fontsize=7.5)
    ax.set_xlabel(f"Binder residue ({binder_chain})", fontsize=9)
    ax.set_title("Contig options — join-gap sensitivity", fontsize=10)

    legend = [
        mpatches.Patch(facecolor=PALETTE["fixed"],  label="Fixed"),
        mpatches.Patch(facecolor=PALETTE["design"], label="De novo / design"),
    ]
    ax.legend(handles=legend, fontsize=8, loc="upper right",
              bbox_to_anchor=(1.0, -0.12), ncol=2, frameon=False)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ── Helpers ─────────────────────────────────────────────────────────────────

def _empty(out_path: Path, message: str) -> None:
    """Write a 'no data' placeholder so the Nextflow output pattern matches."""
    fig, ax = plt.subplots(figsize=(5, 2))
    ax.text(0.5, 0.5, message, ha="center", va="center", fontsize=10)
    ax.axis("off")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ── CLI ─────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results-json",   required=True, type=Path,
                    help="solved_pose_results.json from POSE_SOLVE")
    ap.add_argument("--metrics-json",   required=True, type=Path,
                    help="interface_metrics.json from POSE_INTERFACE_METRICS")
    ap.add_argument("--restart-losses", required=True, type=Path,
                    help="solved_pose_restart_losses.csv from POSE_SOLVE")
    ap.add_argument("--posed-pdb",      required=True, type=Path,
                    help="solved_pose_posed.pdb (needed for the contig "
                         "comparison plot to enumerate binder residues)")
    ap.add_argument("--receptor-chain", default="A")
    ap.add_argument("--effector-chain", default="B")
    return ap.parse_args()


def main() -> int:
    args = parse_args()
    bc = args.receptor_chain.upper()
    tc = args.effector_chain.upper()

    results = json.loads(args.results_json.read_text())
    metrics = json.loads(args.metrics_json.read_text())
    # pose_interface_metrics.py emits a list (one entry per input PDB).
    if isinstance(metrics, dict):
        metrics = [metrics]

    plot_loss_breakdown(results,                          Path("pose_loss_breakdown.png"))
    plot_pair_distances(results, bc, tc,                  Path("pose_pair_distances.png"))
    plot_interface_dashboard(metrics, bc, tc,             Path("pose_interface_dashboard.png"))
    plot_restart_loss_curve(args.restart_losses,          Path("pose_restart_loss_curve.png"))
    plot_contig_comparison(results, args.posed_pdb, bc,   Path("pose_contig_comparison.png"))

    print("Wrote:")
    for fn in ("pose_loss_breakdown.png", "pose_pair_distances.png",
               "pose_interface_dashboard.png", "pose_restart_loss_curve.png",
               "pose_contig_comparison.png"):
        print(f"  {fn}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
