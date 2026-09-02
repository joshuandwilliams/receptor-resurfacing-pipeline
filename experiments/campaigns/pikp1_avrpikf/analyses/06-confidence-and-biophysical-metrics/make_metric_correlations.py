#!/usr/bin/env python3
"""Rank correlation between the two metric families, within the correct set.

Selection only ever chooses among designs Boltz-2 posed correctly, so the
correlations are computed on those 39 and not on all 128. Pooling the two pose
groups inflates every confidence-against-biophysical correlation, because
confidence differs by pose and the biophysical metrics do not, so mixing the
groups induces association that is absent within either.

ra_eff is included so that residual pose quality inside the correct set can be
read against both families. Rosetta binding energy is favourable when
negative, and is left unflipped here so the matrix shows raw Spearman.

Usage:
    make_metric_correlations.py
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "figure-options")
RA_EFF_THRESHOLD = 5.0
INK, MUTED = "#222222", "#666666"

COLUMNS = {
    "Pose": "rep_ra_eff_vs_truth_median",
    "pLDDT": "rep_complex_plddt_median",
    "ipTM": "rep_iptm_median",
    "ipSAE": "rep_ipsae_min_median",
    "Sc": "sc",
    "BSA": "bsa",
    "Rosetta ddG": "rosetta_ddg",
}


def load(correct_only=True):
    d = pd.read_csv(os.path.join(
        HERE, "data", "crystal_full_test_contig_survivors.csv"))
    d = d[~d.mpnn_sequence.str.contains("control")]
    if correct_only:
        return d[d.rep_ra_eff_vs_truth_median <= RA_EFF_THRESHOLD]
    return d


def matrices(d):
    names = list(COLUMNS)
    n = len(names)
    rho = np.zeros((n, n))
    pval = np.ones((n, n))
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            r, p = spearmanr(d[COLUMNS[a]], d[COLUMNS[b]])
            rho[i, j], pval[i, j] = r, p
    return names, rho, pval


def heatmap(names, rho, pval, n_designs, title, stem):
    plt.rcParams.update({
        "figure.dpi": 130, "savefig.dpi": 300, "savefig.bbox": "tight",
        "font.size": 9, "axes.labelsize": 9, "text.color": INK,
        "xtick.color": MUTED, "ytick.color": MUTED})
    masked = np.array(rho, dtype=float)
    np.fill_diagonal(masked, np.nan)
    fig, ax = plt.subplots(figsize=(5.8, 5.0))
    im = ax.imshow(masked, cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(len(names)))
    ax.set_yticks(range(len(names)))
    ax.set_xticklabels(names, rotation=45, ha="right")
    ax.set_yticklabels(names)
    for i in range(len(names)):
        for j in range(len(names)):
            if i == j:
                continue
            # Bold only where the correlation would survive at 0.05, so the
            # eye is not drawn to colour that n = 39 cannot support.
            weight = "bold" if pval[i, j] < 0.05 else "normal"
            ax.text(j, i, f"{rho[i, j]:.2f}", ha="center", va="center",
                    fontsize=8, fontweight=weight,
                    color="white" if abs(rho[i, j]) > 0.55 else INK)
    cb = fig.colorbar(im, ax=ax, shrink=0.8)
    cb.set_label(r"Spearman $\rho$", fontsize=9)
    fig.tight_layout()
    path = os.path.join(OUT, f"{stem}.png")
    fig.savefig(path)
    plt.close(fig)
    print("wrote", path)


# Rows are everything read off the Boltz-2 prediction, including ra_eff RMSD,
# which is measured on that predicted structure against the RFdiffusion
# reference. Columns are computed from the structure independently of Boltz-2.
CONF = ["Pose", "pLDDT", "ipTM", "ipSAE"]
ORTHO = ["Sc", "BSA", "Rosetta ddG"]


def cross_heatmap(d, title, stem):
    """Confidence against the orthogonal metrics only.

    The square matrix is dominated by the strong within-family blocks, which
    are not the question. This keeps only the cross-family cells.
    """
    rho = np.zeros((len(CONF), len(ORTHO)))
    pval = np.ones_like(rho)
    for i, a in enumerate(CONF):
        for j, b in enumerate(ORTHO):
            rho[i, j], pval[i, j] = spearmanr(d[COLUMNS[a]], d[COLUMNS[b]])
    fig, ax = plt.subplots(figsize=(4.4, 3.4))
    im = ax.imshow(rho, cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(len(ORTHO))); ax.set_yticks(range(len(CONF)))
    ax.set_xticklabels(ORTHO, rotation=45, ha="right")
    ax.set_yticklabels(CONF)
    for i in range(len(CONF)):
        for j in range(len(ORTHO)):
            ax.text(j, i, f"{rho[i, j]:.2f}", ha="center", va="center",
                    fontsize=8.5,
                    fontweight="bold" if pval[i, j] < 0.05 else "normal",
                    color="white" if abs(rho[i, j]) > 0.55 else INK)
    cb = fig.colorbar(im, ax=ax, shrink=0.9)
    cb.set_label(r"Spearman $\rho$", fontsize=9)
    ax.set_ylabel("Boltz-2", fontsize=9)
    ax.set_xlabel("Biophysical", fontsize=9)
    fig.tight_layout()
    path = os.path.join(OUT, f"{stem}.png")
    fig.savefig(path); plt.close(fig)
    print("wrote", path)


def report(d, title, stem):
    names, rho, pval = matrices(d)
    tab = pd.DataFrame(rho, index=names, columns=names).round(2)
    tab.to_csv(os.path.join(HERE, f"{stem}.csv"))
    print(f"=== {title}, n = {len(d)} ===")
    print(tab.to_string())
    print("significant at 0.05:")
    hits = [(names[i], names[j], rho[i, j], pval[i, j])
            for i in range(len(names)) for j in range(i + 1, len(names))
            if pval[i, j] < 0.05]
    for a, b, r, p in sorted(hits, key=lambda x: -abs(x[2])):
        print(f"  {a} vs {b}: rho = {r:.2f}, p = {p:.4f}")
    print()
    heatmap(names, rho, pval, len(d), title, stem)
    cross_heatmap(d, title, stem + "_cross")


def main():
    os.makedirs(OUT, exist_ok=True)
    report(load(correct_only=True), "Correctly posed designs only",
           "metric_correlations")
    report(load(correct_only=False), "All designs",
           "metric_correlations_all")


if __name__ == "__main__":
    main()
