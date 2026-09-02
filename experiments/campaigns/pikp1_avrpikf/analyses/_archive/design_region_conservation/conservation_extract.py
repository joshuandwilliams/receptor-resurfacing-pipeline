#!/usr/bin/env python3
"""Per-position conservation of the two design regions against the wild type.

The designs differ in length, so a positional comparison cannot be made by
index. Each design region is aligned to its wild-type counterpart with
BLOSUM62 under the gap settings the thesis uses elsewhere (open -10, extend
-0.5), and every wild-type position is scored by the mean substitution score
of whatever aligned to it across the campaign.

Two numbers come out per wild-type position. The mean BLOSUM62 score says
whether ProteinMPNN kept the residue or something chemically like it. The
aligned fraction says in how many of the designs the alignment placed a residue
opposite that position rather than a gap.

Read the aligned fraction carefully. It is a property of the alignment, not of
the designs. Designs SHORTER than the wild-type 13 must have a gap somewhere,
and the aligner puts it at whichever position costs least, so a low aligned
fraction marks the position the alignment sacrifices rather than a biological
insertion point. Splitting by length makes that plain, which is why
`length_class` is reported alongside.

The profile also only describes the 13 wild-type positions. A design sampling
20 residues has 7 with no wild-type counterpart, and they contribute nothing
here.

Usage:
    conservation_extract.py [--csv data/<run>.scored_metadata.csv]
                            [--output design_region_conservation.csv]
"""

import argparse
import os

import numpy as np
import pandas as pd
from Bio import Align
from Bio.Align import substitution_matrices

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_CSV = os.path.join(
    HERE, "data", "crystal_full_test_contig.scored_metadata.csv")

# Matches the pose-accuracy scoring in the benchmarking repo, which took them
# from the defaults of ChimeraX's matchmaker.
GAP_OPEN = -10.0
GAP_EXTEND = -0.5


def build_aligner():
    al = Align.PairwiseAligner()
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score = GAP_OPEN
    al.extend_gap_score = GAP_EXTEND
    al.mode = "global"
    return al


def split_regions(value):
    """scored_metadata packs the two design regions into one pipe-joined cell."""
    return [part.strip() for part in str(value).split("|")]


def profile(native, designs, aligner, matrix):
    """Mean BLOSUM62 score and aligned fraction per native position.

    `designs` is the list of design-region strings to profile. Callers pass
    either all of them or one length class at a time.
    """
    scores = [[] for _ in native]
    for query in designs:
        alignment = aligner.align(native, query)[0]
        for (ns, ne), (qs, qe) in zip(*alignment.aligned):
            for k in range(ne - ns):
                scores[ns + k].append(float(matrix[native[ns + k], query[qs + k]]))
    n = len(designs)
    return [
        {"wt_residue": aa,
         "wt_position": i + 1,
         "n_aligned": len(v),
         "aligned_fraction": len(v) / n,
         "mean_blosum62": float(np.mean(v)) if v else np.nan,
         "sd_blosum62": float(np.std(v, ddof=1)) if len(v) > 1 else np.nan}
        for i, (aa, v) in enumerate(zip(native, scores))
    ]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default=DEFAULT_CSV)
    ap.add_argument("--output", default=os.path.join(
        HERE, "design_region_conservation.csv"))
    args = ap.parse_args()

    df = pd.read_csv(args.csv)
    aligner = build_aligner()
    matrix = aligner.substitution_matrix

    native = split_regions(df["native_residues"].iloc[0])
    designed = [split_regions(v) for v in df["designed_residues"]]

    rows = []
    for region_idx, wt in enumerate(native, start=1):
        all_q = [d[region_idx - 1] for d in designed]
        classes = {"all": all_q,
                   "shorter": [q for q in all_q if len(q) < len(wt)],
                   "equal": [q for q in all_q if len(q) == len(wt)],
                   "longer": [q for q in all_q if len(q) > len(wt)]}
        for name, queries in classes.items():
            if not queries:
                continue
            for row in profile(wt, queries, aligner, matrix):
                row["design_region"] = region_idx
                row["length_class"] = name
                row["n_designs"] = len(queries)
                rows.append(row)

    out = pd.DataFrame(rows)[
        ["design_region", "length_class", "wt_position", "wt_residue",
         "n_designs", "n_aligned", "aligned_fraction", "mean_blosum62",
         "sd_blosum62"]]
    out.to_csv(args.output, index=False)
    print(f"wrote {args.output} ({len(out)} positions, "
          f"{len(designed)} designs)")


if __name__ == "__main__":
    main()
