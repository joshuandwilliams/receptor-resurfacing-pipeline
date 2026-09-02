#!/usr/bin/env python3
"""Index design region 1 by depth down each arm of the hairpin.

Design region 1 is a U. Two beta strands run down from fixed anchors either
side and a loop turns through 180 degrees at the bottom, and the only thing
RFdiffusion varies is how deep the U goes. Position counted from the start of
the region is therefore the wrong coordinate: it makes a residue's index depend
on how much length was added below it.

The right coordinate is which arm a residue is on and how many steps down from
that arm's anchor it sits. Residues at the same depth on the same arm are
equivalent between designs whatever the loop length, and the extra residues in
a longer design are the ones nearest the turn.

Finding the turn is the only judgement, and it is made on the design's own
geometry so nothing external can shift it. It is found from the pairing
register of the hairpin: in a U the two arms run antiparallel, so a residue i
sits opposite some residue j and i + j is the same constant for every pair,
namely twice the turn. Taking the median of (i + j) / 2 over all residues with
a close cross-arm partner averages that estimate over the whole hairpin.

The obvious alternative, taking the residue whose Ca is furthest from the line
joining the two anchors, does not work. The depth profile is not a single clean
peak but a plateau that wanders by several Angstrom, so argmax lands anywhere
along it: across the six designs of length 17 it put the turn at positions 6,
8, 9, 9, 9 and 12, implying one arm six residues longer than the other. The
pairing register gives 10.0 for all six. Within every length band its spread is
under half a residue, and the estimate rises by almost exactly 0.5 per residue
added, which is what a hairpin extending both arms alternately should do. That
line-distance version is kept below as a fallback for any design with no
cross-arm contacts.

Writes one row per residue with its arm, its depth, and the column it should
occupy in a layout anchored at both ends and as wide as the longest design.

Usage:
    hairpin_layout_extract.py [--max-length 20]
"""

import argparse
import glob
import os
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")

SCAFFOLD_END = 32          # last fixed residue before design region 1
DR1_START = 33
TAIL_LEN = 78 - 45         # fixed residues after the region, from the wild type


def ca_map(path, chain="A"):
    model = PDBParser(QUIET=True).get_structure("s", path)[0]
    return {r.id[1]: np.asarray(r["CA"].coord, dtype=float)
            for r in model[chain] if r.id[0] == " " and "CA" in r}


def distance_to_line(points, a, b):
    """Perpendicular distance from each point to the line through a and b."""
    d = b - a
    d = d / np.linalg.norm(d)
    v = points - a
    proj = np.outer(v @ d, d)
    return np.linalg.norm(v - proj, axis=1)


def turn_from_pairing(points, min_sep=3, max_dist=8.0):
    """Turn index from the antiparallel register, or None if it cannot be set.

    For a residue i and its closest partner j at least min_sep apart in
    sequence, i + j is twice the turn index in an ideal hairpin. The median
    over all such pairs is the estimate. max_dist keeps it to genuine
    cross-arm contacts rather than anything that happens to be nearest.
    """
    n = len(points)
    D = np.linalg.norm(points[:, None, :] - points[None, :, :], axis=-1)
    mids = []
    for i in range(n):
        far = [j for j in range(n) if abs(i - j) >= min_sep]
        if not far:
            continue
        j = min(far, key=lambda j: D[i, j])
        if D[i, j] <= max_dist:
            mids.append((i + j) / 2.0)
    return float(np.median(mids)) if mids else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-length", type=int, default=20)
    ap.add_argument("--output", default=os.path.join(HERE, "hairpin_layout.csv"))
    args = ap.parse_args()

    # The wild type is laid out by exactly the same rules, so its row is a
    # reference produced by the method rather than drawn by hand.
    paths = sorted(glob.glob(os.path.join(DATA, "rfdiffusion", "design_*.pdb")))
    paths.append(os.path.join(DATA, "wt_reference.pdb"))

    rows = []
    for path in paths:
        ca = ca_map(path)
        n_res = max(ca)
        length = n_res - SCAFFOLD_END - TAIL_LEN
        nums = [DR1_START + k for k in range(length)]
        if not all(n in ca for n in nums):
            continue

        anchor_n, anchor_c = ca[SCAFFOLD_END], ca[DR1_START + length]
        pts = np.array([ca[n] for n in nums])
        depth = distance_to_line(pts, anchor_n, anchor_c)
        mid = turn_from_pairing(pts)
        # The N arm takes everything up to and including the turn, so a
        # half-integer register splits the two arms evenly.
        turn = int(np.floor(mid)) if mid is not None else int(np.argmax(depth))

        design = os.path.basename(path).replace(".pdb", "")
        if design == "wt_reference":
            design = "wild-type"
        for k, n in enumerate(nums):
            if k <= turn:
                arm, arm_depth = "N", k + 1
                column = k + 1                       # in from the left anchor
            else:
                arm, arm_depth = "C", length - k
                column = args.max_length - (length - 1 - k)   # in from the right
            rows.append({
                "design": design,
                "region_length": length,
                "design_position": k + 1,
                "arm": arm,
                "arm_depth": arm_depth,
                "is_turn": k == turn,
                "turn_register": mid,
                "dist_to_anchor_line": float(depth[k]),
                "column": column,
            })

    out = pd.DataFrame(rows)

    # Carry the wild-type correspondence across so the layout can be coloured
    # by it. The layout itself does not depend on it.
    corr_path = os.path.join(HERE, "structural_correspondence_per_design.csv")
    if os.path.isfile(corr_path):
        corr = pd.read_csv(corr_path)[
            ["design", "design_position", "matched", "wt_position"]]
        out = out.merge(corr, on=["design", "design_position"], how="left")

    # And the DSSP assignment, which is the more honest thing to colour by.
    # Whether a residue has a wild-type counterpart is a property of the
    # comparison; whether it is strand or loop is a property of the design.
    ss_path = os.path.join(HERE, "design_secondary_structure.csv")
    if os.path.isfile(ss_path):
        ss = pd.read_csv(ss_path)
        ss["design_position"] = ss.resnum - DR1_START + 1
        out = out.merge(ss[["design", "design_position", "ss"]],
                        on=["design", "design_position"], how="left")

    out.to_csv(args.output, index=False)
    print(f"wrote {args.output}")
    n_arm = out.groupby("design").arm.apply(lambda s: (s == "N").sum())
    print(f"{out.design.nunique()} designs, region lengths "
          f"{out.region_length.min()}-{out.region_length.max()}")
    print(f"N-arm length: median {n_arm.median():.0f}, "
          f"range {n_arm.min()}-{n_arm.max()}")
    print(out.groupby("region_length").apply(
        lambda g: pd.Series({
            "designs": g.design.nunique(),
            "N-arm": (g.arm == "N").sum() / g.design.nunique(),
            "C-arm": (g.arm == "C").sum() / g.design.nunique()}),
        include_groups=False).round(1).to_string())


if __name__ == "__main__":
    main()
