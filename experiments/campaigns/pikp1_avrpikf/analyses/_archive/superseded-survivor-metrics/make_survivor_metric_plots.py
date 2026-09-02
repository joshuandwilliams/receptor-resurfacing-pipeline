#!/usr/bin/env python3
"""Metric distributions for the designs that survived negative steering.

Run: crystal_full_test_contig.  Rows whose cross_tier is "none" are dropped, as
are the two sequence controls, leaving the designs that steering actually
validated.  Bars are stacked by how many of the design's three seeds passed,
which is the quantity the tier letters encode.

Deep-learning confidence (pLDDT, ipTM, ipSAE) and biophysical interface metrics
(shape complementarity, Rosetta ddG, buried surface area) are shown together
because the point of the orthogonal-metrics stage is that they need not agree.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = Path(__file__).parent
DATA = HERE / "data" / "crystal_full_test_contig_survivors.csv"
PLOTS = HERE / "plots"

INK, MUTED, GRID = "#222222", "#666666", "#DDDDDD"
# Sequential, so "more seeds passing" reads as "darker" without a colour lookup.
PASS_COLOR = {1: "#C6DBEF", 2: "#6BAED6", 3: "#08519C"}

# (column, panel title, x label, scale factor)
METRICS = [
    ("rep_complex_plddt_median", "Complex pLDDT", "complex pLDDT (0–100)", 100.0),
    ("rep_iptm_median", "ipTM", "ipTM", 1.0),
    ("rep_ipsae_min_median", "ipSAE", "ipSAE (min)", 1.0),
    ("sc", "Shape complementarity", "sc", 1.0),
    ("rosetta_ddg", "Rosetta ΔΔG", "ddG (REU)", 1.0),
    ("bsa", "Buried surface area", "BSA (Å²)", 1.0),
]


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


def main():
    apply_style()
    PLOTS.mkdir(exist_ok=True)

    d = pd.read_csv(DATA)
    d = d[(d.row_type == "steered") & (d.cross_tier.astype(str).str.lower() != "none")]
    d = d.copy()
    d["n_pass"] = d.rep_n_pass.astype(int)
    groups = sorted(d.n_pass.unique())

    fig, axes = plt.subplots(2, 3, figsize=(10.0, 6.0))
    for ax, (col, title, xlab, scale) in zip(axes.ravel(), METRICS):
        vals = pd.to_numeric(d[col], errors="coerce") * scale
        ok = vals.notna()
        # One shared bin edge set per panel, so the stack segments line up.
        edges = np.histogram_bin_edges(vals[ok], bins=12)
        stack = [vals[ok & (d.n_pass == g)].values for g in groups]
        ax.hist(stack, bins=edges, stacked=True,
                color=[PASS_COLOR[g] for g in groups],
                edgecolor="white", linewidth=0.6, zorder=3)
        ax.set_title(title)
        ax.set_xlabel(xlab)
        ax.set_ylabel("designs")
        ax.grid(axis="x", visible=False)
        ax.yaxis.get_major_locator().set_params(integer=True)

    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=PASS_COLOR[g],
                             edgecolor="white", label=f"{g}/3 passing seeds")
               for g in groups]
    fig.legend(handles=handles, loc="upper center", ncol=len(groups),
               frameon=False, bbox_to_anchor=(0.5, 1.0))
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    out = PLOTS / "survivor_metric_histograms.png"
    fig.savefig(out)
    print(f"{out}  n = {len(d)} designs, groups {dict(d.n_pass.value_counts().sort_index())}")

    summary = (d.assign(**{"complex pLDDT": d.rep_complex_plddt_median * 100})
               .groupby("n_pass")[["complex pLDDT", "rep_iptm_median",
                                   "rep_ipsae_min_median", "sc", "rosetta_ddg", "bsa"]]
               .median().round(3))
    summary.index = [f"{i}/3" for i in summary.index]
    summary.to_csv(HERE / "survivor_metric_medians.csv")
    print()
    print(summary.to_string())


if __name__ == "__main__":
    main()
