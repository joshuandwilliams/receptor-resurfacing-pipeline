#!/usr/bin/env python3
"""Profile design region 1 in DESIGN coordinates, positions 1 to 20.

The companion extract profiles the 13 wild-type positions. That answers "what
happened to each wild-type residue" but it cannot describe the extensions,
because a design sampling 20 residues has 7 with no wild-type counterpart and
they contribute nothing to a wild-type-indexed profile. The extensions are the
interesting part, so this extract indexes the other way.

For each design position 1..20:

    n_designs       how many of the 128 designs are long enough to have one
    n_matched       of those, how many put a wild-type residue there
    n_inserted      of those, how many are aligned to a gap in the wild type,
                    which is what "an extra residue" means concretely
    frac_inserted   n_inserted / n_designs
    entropy_bits    Shannon entropy of the residue identity at that position
    top_residue     the most common residue, and its share

The alignment is the same one the wild-type profile uses, BLOSUM62 at gap open
-10 and extend -0.5, so the two extracts agree about which residues are extra.

Usage:
    design_coordinate_extract.py [--region 1] [--max-length 20]
"""

import argparse
import collections
import math
import os

import numpy as np
import pandas as pd
from Bio import Align
from Bio.Align import substitution_matrices

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_CSV = os.path.join(
    HERE, "data", "crystal_full_test_contig.scored_metadata.csv")

GAP_OPEN = -10.0
GAP_EXTEND = -0.5


def build_aligner():
    al = Align.PairwiseAligner()
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score = GAP_OPEN
    al.extend_gap_score = GAP_EXTEND
    al.mode = "global"
    return al


def entropy(letters):
    n = len(letters)
    counts = collections.Counter(letters)
    return -sum(v / n * math.log2(v / n) for v in counts.values())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default=DEFAULT_CSV)
    ap.add_argument("--region", type=int, default=1)
    ap.add_argument("--max-length", type=int, default=20)
    ap.add_argument("--output", default=os.path.join(
        HERE, "design_coordinate_profile.csv"))
    args = ap.parse_args()

    df = pd.read_csv(args.csv)
    al = build_aligner()
    native = str(df["native_residues"].iloc[0]).split("|")[args.region - 1]
    designs = [str(v).split("|")[args.region - 1] for v in df["designed_residues"]]

    L = args.max_length
    present = np.zeros(L, dtype=int)
    matched = np.zeros(L, dtype=int)
    letters = [[] for _ in range(L)]

    for q in designs:
        a = al.align(native, q)[0]
        # design index -> wild-type index, for the pairs the alignment matched
        to_wt = {}
        for (ns, ne), (qs, qe) in zip(*a.aligned):
            for k in range(ne - ns):
                to_wt[qs + k] = ns + k
        for i, ch in enumerate(q[:L]):
            present[i] += 1
            letters[i].append(ch)
            if i in to_wt:
                matched[i] += 1

    rows = []
    for i in range(L):
        if not present[i]:
            continue
        ins = present[i] - matched[i]
        top, count = collections.Counter(letters[i]).most_common(1)[0]
        rows.append({
            "design_region": args.region,
            "design_position": i + 1,
            "n_designs": int(present[i]),
            "n_matched": int(matched[i]),
            "n_inserted": int(ins),
            "frac_inserted": ins / present[i],
            "entropy_bits": entropy(letters[i]),
            "top_residue": top,
            "top_residue_frac": count / present[i],
        })

    out = pd.DataFrame(rows)
    out.to_csv(args.output, index=False)
    print(f"wrote {args.output}")
    print(f"wild-type region length {len(native)}, "
          f"{len(designs)} designs, {int(out.n_inserted.sum())} inserted "
          f"residues in total")


if __name__ == "__main__":
    main()
