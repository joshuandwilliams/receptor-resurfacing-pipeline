#!/usr/bin/env python3
"""ClustalW alignment + visualisation of the design-region sequences.

For each of the two chain-A design regions (region 1 = variable-length,
region 2 = the trailing 6 residues), across the 10 representative negsteer
poses + the crystal input:

  1. Extract the sequence straight from each PDB (same ranges as the
     ChimeraX colouring script).
  2. Align with ClustalW 2.1.
  3. Render: coloured MSA grid, per-column sequence logo, the ClustalW
     guide-tree dendrogram, and a pairwise %-identity heatmap.

Run with the clustalw conda env's python (or any python3 with biopython,
logomaker, matplotlib); the clustalw2 binary path is set below.
"""
from pathlib import Path
import re
import subprocess
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import logomaker
import pandas as pd
from Bio import Phylo, AlignIO

CLUSTALW = "/Users/jowillia/miniforge3/envs/clustalw/bin/clustalw2"
HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1]  # .../crystal_full_test_contig/results
NS = RESULTS / "negative_steering" / "runs"
INPUT = Path("/Users/jowillia/Documents/GitHub/receptor-resurfacing-pipeline/"
             "experiments/campaigns/pikp1_avrpikf/inputs/"
             "pikp1_avrpikf_pik_interface_complex.pdb")

AA3TO1 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E',
          'GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F',
          'PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}

# (id, pdb, region1 (lo,hi), region2 (lo,hi))
MODELS = [
    ("d39_s0", NS/"design_39_seq_0/cycle_0/initial_prediction_s2.pdb", (33,48),(76,81)),
    ("d17_s0", NS/"design_17_seq_0/cycle_0/initial_prediction_s1.pdb", (33,46),(74,79)),
    ("d47_s0", NS/"design_47_seq_0/cycle_0/initial_prediction_s1.pdb", (33,48),(76,81)),
    ("d39_s1", NS/"design_39_seq_1/cycle_0/steered/design_00_s1/prediction.pdb", (33,48),(76,81)),
    ("d26_s1", NS/"design_26_seq_1/cycle_0/steered/design_01_s0/prediction.pdb", (33,51),(79,84)),
    ("d18_s0", NS/"design_18_seq_0/cycle_0/steered/design_00_s0/prediction.pdb", (33,42),(70,75)),
    ("d49_s0", NS/"design_49_seq_0/cycle_0/initial_prediction_s2.pdb", (33,46),(74,79)),
    ("d40_s0", NS/"design_40_seq_0/cycle_0/initial_prediction_s1.pdb", (33,47),(75,80)),
    ("d6_s1",  NS/"design_6_seq_1/cycle_0/initial_prediction_s2.pdb", (33,47),(75,80)),
    ("d14_s0", NS/"design_14_seq_0/cycle_0/steered/design_02_s0/prediction.pdb", (33,47),(75,80)),
    ("crystal_input", INPUT, (33,45),(73,78)),
]

# Clustal-X-ish residue colours (by physicochemical group).
RESIDUE_COLOUR = {
    **{a: "#80a0f0" for a in "AILMFWV"},   # hydrophobic - blue
    **{a: "#f01505" for a in "KR"},          # positive - red
    **{a: "#c048c0" for a in "ED"},          # negative - magenta
    **{a: "#15c015" for a in "STNQ"},        # polar - green
    "G": "#f09048", "P": "#c0c000",          # glycine orange, proline yellow
    "C": "#f08080", "H": "#15a4a4", "Y": "#15a4a4",
    "-": "#ffffff",
}

# Pastel version of above (50% blend with white) — for slides.
RESIDUE_COLOUR_PALE = {
    **{a: "#bfcff7" for a in "AILMFWV"},
    **{a: "#f78a82" for a in "KR"},
    **{a: "#dfa3df" for a in "ED"},
    **{a: "#8adf8a" for a in "STNQ"},
    "G": "#f7c7a3", "P": "#dfdf7f",
    "C": "#f7bfbf", "H": "#8ad1d1", "Y": "#8ad1d1",
    "-": "#ffffff",
}


def _rename_label(mid):
    """crystal_input -> 'Input Design'; d{N}_s{M} -> 'Design N Sequence M'."""
    if mid == "crystal_input":
        return "Input Design"
    m = re.match(r"d(\d+)_s(\d+)$", mid)
    if m:
        return f"Design {m.group(1)} Sequence {m.group(2)}"
    return mid


def extract(path, lo, hi):
    res = {}
    for line in open(path):
        if line.startswith("ATOM") and line[12:16].strip() == "CA" and line[21] == "A":
            res[int(line[22:26])] = AA3TO1.get(line[17:20].strip(), "X")
    return "".join(res.get(i, "-") for i in range(lo, hi + 1))


def run_clustalw(fasta):
    aln = fasta.with_suffix(".aln")
    subprocess.run([CLUSTALW, f"-INFILE={fasta}", "-ALIGN", "-TYPE=PROTEIN",
                    f"-OUTFILE={aln}", "-OUTPUT=CLUSTAL"],
                   check=True, capture_output=True, text=True)
    return aln, fasta.with_suffix(".dnd")


def plot_msa_grid(records, out, title, colour_map=None, show_title=True):
    if colour_map is None:
        colour_map = RESIDUE_COLOUR
    names = [r[0] for r in records]
    seqs = [r[1] for r in records]
    ncol = len(seqs[0]); nrow = len(seqs)
    fig, ax = plt.subplots(figsize=(max(6, ncol * 0.42), nrow * 0.42 + 0.6))
    for i, seq in enumerate(seqs):
        y = nrow - 1 - i
        for j, aa in enumerate(seq):
            ax.add_patch(Rectangle((j, y), 1, 1,
                         facecolor=colour_map.get(aa, "#dddddd"),
                         edgecolor="white", linewidth=0.5))
            if aa != "-":
                ax.text(j + 0.5, y + 0.5, aa, ha="center", va="center",
                        fontsize=8, family="monospace")
    ax.set_xlim(0, ncol); ax.set_ylim(0, nrow)
    ax.set_yticks([nrow - 0.5 - i for i in range(nrow)])
    ax.set_yticklabels(names, fontsize=8, family="monospace")
    ax.set_xticks([j + 0.5 for j in range(ncol)])
    ax.set_xticklabels(range(1, ncol + 1), fontsize=7)
    ax.set_xlabel("Alignment column")
    if show_title:
        ax.set_title(title, fontsize=11)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    plt.tight_layout(); plt.savefig(out, dpi=150, bbox_inches="tight"); plt.close()


def plot_logo(records, out, title):
    seqs = [r[1] for r in records]
    ncol = len(seqs[0])
    aas = "ACDEFGHIKLMNPQRSTVWY"
    counts = np.zeros((ncol, len(aas)))
    for seq in seqs:
        for j, aa in enumerate(seq):
            if aa in aas:
                counts[j, aas.index(aa)] += 1
    df = pd.DataFrame(counts, columns=list(aas))
    # information-content logo (bits)
    info = logomaker.transform_matrix(df, from_type="counts",
                                       to_type="information")
    fig, ax = plt.subplots(figsize=(max(6, ncol * 0.5), 2.8))
    logomaker.Logo(info, ax=ax, color_scheme="chemistry")
    ax.set_ylabel("bits"); ax.set_xlabel("Alignment column")
    ax.set_title(title, fontsize=11)
    plt.tight_layout(); plt.savefig(out, dpi=150, bbox_inches="tight"); plt.close()


def plot_tree(dnd, out, title):
    tree = Phylo.read(dnd, "newick")
    tree.ladderize()
    fig, ax = plt.subplots(figsize=(7, 5))
    Phylo.draw(tree, axes=ax, do_show=False,
               label_func=lambda c: c.name if c.is_terminal() else "")
    ax.set_title(title, fontsize=11)
    plt.tight_layout(); plt.savefig(out, dpi=150, bbox_inches="tight"); plt.close()


def plot_identity_heatmap(records, out, title):
    names = [r[0] for r in records]
    seqs = [r[1] for r in records]
    n = len(seqs)
    M = np.zeros((n, n))
    for i in range(n):
        for k in range(n):
            cols = [(a, b) for a, b in zip(seqs[i], seqs[k]) if a != "-" or b != "-"]
            same = sum(1 for a, b in cols if a == b and a != "-")
            M[i, k] = 100.0 * same / len(cols) if cols else 0.0
    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.imshow(M, cmap="viridis", vmin=0, vmax=100)
    ax.set_xticks(range(n)); ax.set_xticklabels(names, rotation=90, fontsize=7)
    ax.set_yticks(range(n)); ax.set_yticklabels(names, fontsize=7)
    for i in range(n):
        for k in range(n):
            ax.text(k, i, f"{M[i,k]:.0f}", ha="center", va="center",
                    fontsize=6, color="white" if M[i, k] < 60 else "black")
    fig.colorbar(im, ax=ax, label="% identity")
    ax.set_title(title, fontsize=11)
    plt.tight_layout(); plt.savefig(out, dpi=150, bbox_inches="tight"); plt.close()


def process(region_idx, tag, label):
    fasta = HERE / f"{tag}.fasta"
    with open(fasta, "w") as fh:
        for mid, pdb, *regs in MODELS:
            seq = extract(pdb, *regs[region_idx])
            fh.write(f">{mid}\n{seq}\n")
    aln, dnd = run_clustalw(fasta)
    msa = AlignIO.read(aln, "clustal")
    records = [(rec.id, str(rec.seq)) for rec in msa]
    if region_idx == 0:
        msa_records = sorted(records, key=lambda r: (0 if r[0] == "crystal_input" else 1))
        msa_records = [(_rename_label(r[0]), r[1]) for r in msa_records]
        plot_msa_grid(msa_records, HERE / f"{tag}_msa.png", "",
                      colour_map=RESIDUE_COLOUR_PALE, show_title=False)
    else:
        plot_msa_grid(records, HERE / f"{tag}_msa.png", f"{label} - ClustalW alignment")
    plot_logo(records, HERE / f"{tag}_logo.png", f"{label} - sequence logo")
    plot_identity_heatmap(records, HERE / f"{tag}_identity.png",
                          f"{label} - pairwise % identity")
    if dnd.exists():
        plot_tree(dnd, HERE / f"{tag}_tree.png", f"{label} - ClustalW guide tree")
    print(f"[{tag}] {len(records)} seqs, {msa.get_alignment_length()} cols -> "
          f"{tag}_msa.png, _logo.png, _identity.png, _tree.png")


def main():
    process(0, "region1", "Design region 1 (red)")
    process(1, "region2", "Design region 2 (purple)")
    print("Done. Outputs in", HERE)


if __name__ == "__main__":
    main()
