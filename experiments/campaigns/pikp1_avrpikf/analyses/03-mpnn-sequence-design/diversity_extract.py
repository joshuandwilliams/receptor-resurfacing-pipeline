#!/usr/bin/env python3
"""Sequence diversity of the designed regions, from full pairwise identity.

The pipeline's own clustering (bin/mpnn_cluster_sequences.py) runs MMseqs2
easy-cluster with --min-seq-id T and -c 0.8 --cov-mode 1. Those two settings
measure different things: the identity applies only to the part that was
aligned, and the coverage rule lets a fifth of the shorter sequence go
unaligned. On design regions of 16 to 26 residues that means 3 to 5 residues
are simply not examined, so a pair reported at "100% identity" can differ at
several positions. The pipeline output is left alone; this is a second opinion
computed a different way.

The two design regions are handled separately rather than concatenated. All
the length variation lives in design region 1, so concatenating lets a global
aligner absorb that variation into the design region 2 block and align
residues of one region against the other. It does so in 873 of the 8128 pairs.

Design region 1 varies from 10 to 20 residues, so pairs are aligned globally,
end to end, and identity is the number of identical aligned positions divided
by the alignment length. A length difference counts against the pair rather
than being excused. There is no coverage threshold and no k-mer prefilter.

Design region 2 is a fixed 6 residues in every design, so a gap can never be
justified and no alignment is needed. Pairs are compared position by position
and the result is a count out of 6. The aligner opened a gap in 4 of the 8128
pairs when it was allowed to, returning identities of 3/7, 4/7 and 5/7.

Clustering uses greedy set cover, the same rule as MMseqs2's --cluster-mode 0:
take the sequence with the most neighbours as a representative, remove it and
its neighbours, repeat. Ties are broken on the sequence name so the count is
reproducible. Connected components are reported alongside, since single-linkage
chains and greedy does not, and the gap between them says how much the answer
depends on that choice.

Usage:
    diversity_extract.py [--csv data/<run>.scored_metadata.csv]
"""

import argparse
import itertools
import os

import numpy as np
import pandas as pd
from Bio import Align

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_CSV = os.path.join(
    HERE, "..", "..", "runs", "crystal_full_test_contig", "results",
    "sequences", "scored_metadata.csv")
THRESHOLDS = [i / 10 for i in range(11)]


def build_aligner():
    al = Align.PairwiseAligner()
    al.mode = "global"
    al.match_score = 1
    al.mismatch_score = 0
    al.open_gap_score = -1
    al.extend_gap_score = -0.5
    return al


def identity(a, b, al):
    """Identical aligned positions over alignment length, gaps included."""
    if a == b:
        return 1.0
    aln = al.align(a, b)[0]
    top, bottom = str(aln[0]), str(aln[1])
    same = sum(1 for x, y in zip(top, bottom) if x == y and x != "-")
    return same / len(top)


def matches(a, b):
    """Identical positions between two sequences of the same length."""
    return sum(1 for x, y in zip(a, b) if x == y)


def greedy_cover(names, neighbours):
    """MMseqs2 --cluster-mode 0: repeatedly take the best-connected sequence."""
    remaining, clusters = set(names), 0
    while remaining:
        rep = max(sorted(remaining),
                  key=lambda n: len(neighbours[n] & remaining))
        clusters += 1
        remaining -= ({rep} | neighbours[rep])
    return clusters


def components(names, neighbours):
    seen, count = set(), 0
    for n in names:
        if n in seen:
            continue
        count += 1
        stack = [n]
        while stack:
            x = stack.pop()
            if x in seen:
                continue
            seen.add(x)
            stack.extend(neighbours[x] - seen)
    return count


def within_between(pair_score, scale=1):
    """Median pair score for sequences on one backbone against all others.

    Names are design_<d>_seq_<s>, so the design field says whether a pair
    shares an RFdiffusion backbone. The two medians separate what the
    sequence designer contributes from what the backbone set contributes.
    """
    within, between = [], []
    for (a, b), v in pair_score.items():
        design = a.rsplit("_seq_", 1)[0]
        (within if b.startswith(design + "_seq_") else between).append(
            v / scale)
    return float(np.median(within)), float(np.median(between))


def cluster_curve(names, pair_score, thresholds, scale=1):
    """Cluster and component counts at each threshold on a pair score.

    scale converts a threshold expressed as a fraction into the units the
    pair score is held in, so a match count out of 6 can be thresholded on
    identity without being rescaled first.
    """
    rows = []
    for t in thresholds:
        nb = {n: set() for n in names}
        for (a, b), v in pair_score.items():
            if v >= t * scale - 1e-9:
                nb[a].add(b)
                nb[b].add(a)
        rows.append({"threshold": t,
                     "n_clusters": greedy_cover(names, nb),
                     "n_components": components(names, nb),
                     "n_sequences": len(names)})
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default=DEFAULT_CSV)
    ap.add_argument("--outdir", default=HERE)
    args = ap.parse_args()

    df = pd.read_csv(args.csv)
    regions = {}
    for i in (0, 1):
        regions[i + 1] = {
            f"design_{r.design}_seq_{r.seq}":
            str(r.designed_residues).split("|")[i] for r in df.itertuples()}

    r1 = regions[1]
    names = sorted(r1)
    al = build_aligner()
    pair_id = {(a, b): identity(r1[a], r1[b], al)
               for a, b in itertools.combinations(names, 2)}
    vals = np.array(list(pair_id.values()))
    print(f"design region 1: {len(names)} sequences, "
          f"{len(set(r1.values()))} distinct, "
          f"lengths {min(map(len, r1.values()))}-{max(map(len, r1.values()))}")
    print(f"  pairwise identity: median {np.median(vals):.3f}, "
          f"max {vals.max():.3f}, pairs at 1.000: {(vals == 1).sum()}")
    w, b = within_between(pair_id)
    print(f"  median identity: within backbone {w:.3f}, "
          f"between backbones {b:.3f}")
    curve = cluster_curve(names, pair_id, THRESHOLDS)
    curve_path = os.path.join(args.outdir, "region1_cluster_counts.csv")
    curve.to_csv(curve_path, index=False)
    print(curve.to_string(index=False))

    r2 = regions[2]
    lengths = {len(v) for v in r2.values()}
    if len(lengths) != 1:
        raise SystemExit(f"design region 2 is not a fixed length: {lengths}")
    width = lengths.pop()
    pair_match = {(a, b): matches(r2[a], r2[b])
                  for a, b in itertools.combinations(names, 2)}
    counts = pd.Series(list(pair_match.values())).value_counts()
    hist = pd.DataFrame({"n_matching": range(width + 1)})
    hist["n_pairs"] = [int(counts.get(k, 0)) for k in hist.n_matching]
    hist["pct_pairs"] = 100 * hist.n_pairs / hist.n_pairs.sum()
    hist["region_length"] = width
    hist_path = os.path.join(args.outdir, "region2_match_counts.csv")
    hist.to_csv(hist_path, index=False)
    print(f"\ndesign region 2: {len(names)} sequences, "
          f"{len(set(r2.values()))} distinct, fixed length {width}")
    print(hist.to_string(index=False))

    # Identity is exactly matches / width here, so the curve is the same
    # quantity as design region 1 and belongs on the same axis. Only the
    # width + 1 attainable values exist, which is what makes it a step.
    w2, b2 = within_between(pair_match, scale=width)
    print(f"  median identity: within backbone {w2:.3f}, "
          f"between backbones {b2:.3f}")
    r2_curve = cluster_curve(names, pair_match,
                             [k / width for k in range(width + 1)],
                             scale=width)
    r2_curve.insert(1, "n_matching", range(width + 1))
    r2_curve_path = os.path.join(args.outdir, "region2_cluster_counts.csv")
    r2_curve.to_csv(r2_curve_path, index=False)
    print(r2_curve.to_string(index=False))

    print(f"\nwrote {curve_path}")
    print(f"wrote {hist_path}")
    print(f"wrote {r2_curve_path}")


if __name__ == "__main__":
    main()
