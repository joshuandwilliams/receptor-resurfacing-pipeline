"""Loading and plotting for the design-region novelty comparison.

Three runs of the same campaign differ only in how the design regions were
chosen: one large hand-picked contig, and two auto-derived from effector contact
at 5 A and 3 A.  The question is whether that choice changes what RFDiffusion is
able to build.

Everything here reads `data/<run>.json`, the per-run `rfdiffusion_metrics.json`
the pipeline already writes, so nothing is recomputed from structures.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

DATA = Path(__file__).parent / "data"
PLOTS = Path(__file__).parent / "plots"

# Ordered from fewest, largest design regions to most, smallest.
RUNS = {
    "full (2 regions)": "crystal_full_test_contig",
    "auto 5 A (4 regions)": "crystal_auto5A_contig",
    "auto 3 A (6 regions)": "crystal_auto3A_contig",
}
RUN_COLOR = {"full (2 regions)": "#0072B2",
             "auto 5 A (4 regions)": "#E69F00",
             "auto 3 A (6 regions)": "#D55E00"}

INK, MUTED, GRID = "#222222", "#666666", "#DDDDDD"


def apply_style():
    plt.rcParams.update({
        "figure.dpi": 130, "savefig.dpi": 200, "savefig.bbox": "tight",
        "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
        "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
        "xtick.color": MUTED, "ytick.color": MUTED,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
        "axes.axisbelow": True, "axes.spines.top": False,
        "axes.spines.right": False, "legend.frameon": False,
    })


def load_regions() -> pd.DataFrame:
    """One row per (run, design, design region).

    `com_displacement` is the distance between the built region's Ca centroid
    and the centroid of the input residues it replaced, in the scaffold-aligned
    frame. `span_design` and `span_input` are the corresponding first-to-last
    Ca distances, which is what tells you whether a region was stretched or
    pulled tight rather than merely moved.
    """
    rows = []
    for label, run in RUNS.items():
        d = json.loads((DATA / f"{run}.json").read_text())
        for des in d["designs"]:
            com = des.get("design_region_com_displacement") or []
            lens = des.get("per_region_lengths") or []
            sp_d = des.get("endpoint_distance_design") or []
            sp_i = des.get("endpoint_distance_input") or []
            for i, c in enumerate(com):
                rows.append({
                    "run": label, "design": des["design"], "region": i + 1,
                    "com_displacement": np.nan if c is None else float(c),
                    "length": lens[i] if i < len(lens) else np.nan,
                    "span_design": _num(sp_d[i] if i < len(sp_d) else None),
                    "span_input": _num(sp_i[i] if i < len(sp_i) else None),
                    "passes_filter": bool(des.get("passes_filter")),
                })
    df = pd.DataFrame(rows)
    df["com_per_residue"] = df.com_displacement / df.length
    df["span_ratio"] = df.span_design / df.span_input
    return df


def load_designs() -> pd.DataFrame:
    """One row per (run, design), for the run-level summaries."""
    rows = []
    for label, run in RUNS.items():
        d = json.loads((DATA / f"{run}.json").read_text())
        for des in d["designs"]:
            lens = des.get("per_region_lengths") or []
            rows.append({
                "run": label, "design": des["design"],
                "n_regions": len(lens), "total_designed": int(np.sum(lens)) if lens else 0,
                "motif_rmsd": _num(des.get("motif_rmsd")),
                "design_region_coverage": _num(des.get("design_region_coverage")),
                "frac_contacts_in_design": _num(des.get("frac_contacts_in_design")),
                "passes_filter": bool(des.get("passes_filter")),
            })
    return pd.DataFrame(rows)


def region_com_spread() -> pd.DataFrame:
    """Diversity: how far apart are different designs' versions of one region?

    Mean pairwise distance between the region centroids of every pair of designs
    within a run. Pairwise Ca RMSD would be the obvious diversity measure but is
    undefined here, because variable-length contigs mean two designs' versions of
    the same region often have different residue counts.
    """
    rows = []
    for label, run in RUNS.items():
        d = json.loads((DATA / f"{run}.json").read_text())
        per_region: dict[int, list] = {}
        for des in d["designs"]:
            for i, coords in enumerate(des.get("per_region_coords") or []):
                if coords:
                    per_region.setdefault(i, []).append(
                        np.asarray(coords, float).mean(axis=0))
        for i, cents in per_region.items():
            c = np.asarray(cents)
            if len(c) < 2:
                continue
            dist = np.linalg.norm(c[:, None, :] - c[None, :, :], axis=-1)
            iu = np.triu_indices(len(c), k=1)
            rows.append({"run": label, "region": i + 1, "n_designs": len(c),
                         "mean_pairwise_com_dist": float(dist[iu].mean()),
                         "max_pairwise_com_dist": float(dist[iu].max())})
    return pd.DataFrame(rows)


def _num(v):
    try:
        f = float(v)
        return np.nan if f != f else f
    except (TypeError, ValueError):
        return np.nan


# ── Plot helpers ─────────────────────────────────────────────────────────────

def _box_by_run(ax, df, value, by_region=True):
    """Boxes grouped by run, optionally split by region index within each run."""
    labels, data, colors = [], [], []
    for run in RUNS:
        sub = df[df.run == run]
        if by_region:
            for reg in sorted(sub.region.unique()):
                v = sub[sub.region == reg][value].dropna().values
                if v.size:
                    labels.append(f"{reg}")
                    data.append(v)
                    colors.append(RUN_COLOR[run])
        else:
            v = sub[value].dropna().values
            if v.size:
                labels.append(run)
                data.append(v)
                colors.append(RUN_COLOR[run])
    bp = ax.boxplot(data, patch_artist=True, widths=0.62, showfliers=False,
                    medianprops=dict(color=INK, linewidth=1.6))
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c)
        patch.set_alpha(0.45)
        patch.set_edgecolor(c)
    for i, (v, c) in enumerate(zip(data, colors), start=1):
        jit = np.random.default_rng(i).uniform(-0.13, 0.13, len(v))
        ax.scatter(i + jit, v, s=12, color=c, alpha=0.55, zorder=3,
                   edgecolor="white", linewidth=0.3)
    ax.set_xticks(range(1, len(labels) + 1))
    ax.set_xticklabels(labels)
    ax.grid(axis="x", visible=False)
    return labels, colors


# ── Contig map ───────────────────────────────────────────────────────────────
# The three strategies differ in where the design regions sit as well as how big
# they are, and neither summary table shows placement. This draws the contig
# strings themselves, so the length and the position arguments can be read off
# the same figure.

CONTIGS = {
    "full (2 regions)": "A1-32/10-20/A46-72/6-6",
    "auto 5 A (4 regions)": "A1-2/2-5/A6-31/6-12/A40-42/5-11/A50-67/8-17",
    "auto 3 A (6 regions)": "A1-2/1-1/A4-31/1-3/A34-38/1-1/A40-42/1-1/A44-67/2-5/A71-73/4-8",
}
FIXED_COLOR = "#BBBBBB"


def parse_contig(contig: str):
    """-> [(kind, start, end)] over WILD-TYPE receptor numbering.

    A fixed block `A<i>-<j>` keeps wild-type residues i..j. A design block
    `<lo>-<hi>` replaces whatever wild-type residues lie between the flanking
    fixed blocks, which is what has to be drawn: the design region's footprint
    on the input, not its length in the output.
    """
    segs = contig.split()[0].split("/")
    out, cursor = [], None
    for k, seg in enumerate(segs):
        if seg.startswith("A"):
            lo, hi = (int(v) for v in seg[1:].split("-"))
            out.append(("fixed", lo, hi, seg))
            cursor = hi
        else:
            lo, hi = (int(v) for v in seg.split("-"))
            nxt = next((s for s in segs[k + 1:] if s.startswith("A")), None)
            start = (cursor or 0) + 1
            end = int(nxt[1:].split("-")[0]) - 1 if nxt else None
            out.append(("design", start, end, f"{lo}-{hi}"))
    return out


def fig_contig_map(wt_len: int = 78, figsize=(8.2, 2.6)):
    """One row per strategy, drawn along wild-type receptor numbering."""
    fig, ax = plt.subplots(figsize=figsize)
    for row, (label, contig) in enumerate(CONTIGS.items()):
        y = len(CONTIGS) - 1 - row
        for kind, start, end, txt in parse_contig(contig):
            end = end if end is not None else wt_len
            width = max(end - start + 1, 0.8)
            ax.barh(y, width, left=start - 0.5, height=0.55,
                    color=FIXED_COLOR if kind == "fixed" else RUN_COLOR[label],
                    edgecolor="white", linewidth=0.6, zorder=3)
            if kind == "design":
                ax.text(start - 0.5 + width / 2, y + 0.38, txt, ha="center",
                        va="bottom", fontsize=7.5, color=RUN_COLOR[label])
    ax.set_yticks(range(len(CONTIGS)))
    ax.set_yticklabels(list(CONTIGS)[::-1])
    ax.set_xlabel("wild-type Pikp-1 HMA residue")
    ax.set_xlim(0, wt_len + 1)
    ax.set_ylim(-0.6, len(CONTIGS) - 0.2)
    ax.grid(axis="y", visible=False)
    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=FIXED_COLOR),
               plt.Rectangle((0, 0), 1, 1, facecolor=MUTED)]
    ax.legend(handles, ["retained", "rebuilt (label = allowed length)"],
              loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, frameon=False)
    fig.tight_layout()
    return fig
