#!/usr/bin/env python3
"""Secondary structure of the two design regions, per design.

Reads the ChimeraX DSSP assignment written by dssp_extract.py and the run's
scored metadata, which carries the observed length of each design region as
"<dr1>|<dr2>". Design region 1 starts at residue 33 and runs for its observed
length. Design region 2 is the final 6 residues of the receptor chain.

An earlier pass used pydssp 0.9.1 rather than ChimeraX and reported 57 designs
with a two-residue strand at design region 2. pydssp implements only the
backbone hydrogen-bond term and emits three states, and it disagrees with
ChimeraX on whether a third residue at the strand edge counts. The numbers
here are the ChimeraX ones.

Usage:
    design_region_summary.py
"""

import argparse
import collections
import os

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DR1_START = 33
DR2_LENGTH = 6


def code(ss):
    return {"strand": "S", "helix": "H"}.get(ss, "C")


def segments(pattern, char="S"):
    """Number of runs of char in pattern."""
    return sum(1 for i, c in enumerate(pattern)
               if c == char and (i == 0 or pattern[i - 1] != char))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dssp", default=os.path.join(
        HERE, "design_secondary_structure.csv"))
    ap.add_argument("--metadata", default=os.path.join(
        HERE, "..", "..", "runs", "crystal_full_test_contig",
        "results", "sequences", "scored_metadata.csv"))
    ap.add_argument("--output", default=os.path.join(
        HERE, "design_region_summary.csv"))
    args = ap.parse_args()

    ss = pd.read_csv(args.dssp).sort_values(["design", "resnum"])
    meta = pd.read_csv(args.metadata).drop_duplicates("design")
    lengths = {f"design_{r.design}":
               int(str(r.design_region_length_observed).split("|")[0])
               for r in meta.itertuples()}

    rows = []
    for design, grp in ss[ss.design != "wild-type"].groupby("design"):
        pattern = "".join(code(v) for v in grp.ss)
        length = lengths[design]
        dr1 = pattern[DR1_START - 1:DR1_START - 1 + length]
        dr2 = pattern[-DR2_LENGTH:]
        rows.append({"design": design,
                     "dr1_length": length,
                     "dr1_pattern": dr1,
                     "dr1_strand_residues": dr1.count("S"),
                     "dr1_strand_segments": segments(dr1),
                     "dr2_pattern": dr2,
                     "dr2_strand_residues": dr2.count("S")})
    out = pd.DataFrame(rows).sort_values("dr1_length").reset_index(drop=True)
    out.to_csv(args.output, index=False)

    wt = ss[ss.design == "wild-type"]
    wt_pattern = "".join(code(v) for v in wt.ss)
    print(f"wild type, {len(wt_pattern)} residues")
    print(f"  design region 1 equivalent: {wt_pattern[DR1_START - 1:DR1_START - 1 + 13]}")
    print(f"  design region 2 equivalent: {wt_pattern[-DR2_LENGTH:]}")

    print(f"\n{len(out)} designs")
    two_seg = (out.dr1_strand_segments == 2).sum()
    print(f"  design region 1 with two strand segments: {two_seg}")
    print("  median design region 1 strand residues by length band:")
    band = pd.cut(out.dr1_length, [9, 12, 15, 17, 20],
                  labels=["10-12", "13-15", "16-17", "18-20"])
    print(out.groupby(band, observed=True)
          .dr1_strand_residues.agg(["size", "median", "min", "max"])
          .to_string())

    print("\n  design region 2 patterns:")
    for pat, n in collections.Counter(out.dr2_pattern).most_common():
        print(f"    {pat}  x{n}")
    print(f"\nwrote {args.output}")


if __name__ == "__main__":
    main()
