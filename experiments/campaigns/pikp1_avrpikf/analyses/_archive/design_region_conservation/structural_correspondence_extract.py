#!/usr/bin/env python3
"""Match design-region residues to wild-type residues by STRUCTURE, not sequence.

Sequence alignment is the wrong instrument here. ProteinMPNN recovers about 12%
of the wild-type sequence in the designed region, so BLOSUM62 has almost no
signal to work with and the gap it opens lands wherever the local scoring
happens to prefer. Per design that is often sensible, but the gap position
wanders from design to design, so pooling 128 of them smears insertions across
every position and hides the pattern the fold actually has.

The fold gives a better correspondence for free. Every design keeps the same
fixed scaffold either side of design region 1, and the region is two beta
strands with a loop between them. Superimposing a design on the wild type
using the fixed N-terminal block puts the two strands on top of each other, and
then a design residue either has a wild-type Ca close enough to be its
structural equivalent or it does not. Residues with no equivalent are the
extension, and they fall where the loop is rather than being spread out.

For each design position this writes:

    n_designs     designs long enough to have a residue there
    n_matched     of those, how many sit within MATCH_CUTOFF of a wild-type Ca
    n_inserted    the rest, which have no wild-type counterpart
    median_dist   median distance to the nearest wild-type Ca

Usage:
    structural_correspondence_extract.py [--cutoff 2.5]
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

# Chain A residue numbering is sequential from 1 in both the reference and the
# designs. The scaffold before design region 1 is identical in every design.
SCAFFOLD = list(range(1, 33))
WT_DR1 = list(range(33, 46))       # the wild-type 13-mer
DR1_START = 33
MATCH_CUTOFF = 2.5                 # Angstrom between superimposed Ca atoms
MAX_LENGTH = 20


def ca_map(path, chain="A"):
    """residue number -> Ca coordinate."""
    model = PDBParser(QUIET=True).get_structure("s", path)[0]
    return {r.id[1]: r["CA"].coord for r in model[chain]
            if r.id[0] == " " and "CA" in r}


def kabsch(P, Q):
    """Rotation and translation putting P onto Q."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, qc - R @ pc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cutoff", type=float, default=MATCH_CUTOFF)
    ap.add_argument("--output", default=os.path.join(
        HERE, "structural_correspondence.csv"))
    args = ap.parse_args()

    wt = ca_map(os.path.join(DATA, "wt_reference.pdb"))
    wt_scaffold = np.array([wt[i] for i in SCAFFOLD])
    wt_region = np.array([wt[i] for i in WT_DR1])

    present = np.zeros(MAX_LENGTH, dtype=int)
    matched = np.zeros(MAX_LENGTH, dtype=int)
    dists = [[] for _ in range(MAX_LENGTH)]
    per_design = []   # long format, one row per design per position

    paths = sorted(glob.glob(os.path.join(DATA, "rfdiffusion", "design_*.pdb")))
    lengths = []
    for path in paths:
        d = ca_map(path)
        if not all(i in d for i in SCAFFOLD):
            continue
        R, t = kabsch(np.array([d[i] for i in SCAFFOLD]), wt_scaffold)

        # Design region 1 runs from DR1_START to wherever the fixed C-terminal
        # block resumes. Its length is what RFdiffusion sampled, so it is read
        # off the model rather than assumed.
        n_res = max(d)
        length = n_res - len(SCAFFOLD) - (78 - 45)
        lengths.append(length)
        for k in range(length):
            num = DR1_START + k
            if num not in d or k >= MAX_LENGTH:
                continue
            xyz = R @ d[num] + t
            d_all = np.linalg.norm(wt_region - xyz, axis=1)
            j = int(np.argmin(d_all))
            dist = float(d_all[j])
            present[k] += 1
            dists[k].append(dist)
            if dist <= args.cutoff:
                matched[k] += 1
            per_design.append({
                "design": os.path.basename(path).replace(".pdb", ""),
                "region_length": length,
                "design_position": k + 1,
                "dist": dist,
                "matched": bool(dist <= args.cutoff),
                # Which wild-type residue this one sits on, so the rows can be
                # laid out as a real alignment with the insertions in their own
                # columns rather than every row starting flush left.
                "wt_position": (WT_DR1[j] - DR1_START + 1
                                if dist <= args.cutoff else None)})

    rows = []
    for k in range(MAX_LENGTH):
        if not present[k]:
            continue
        rows.append({
            "design_position": k + 1,
            "n_designs": int(present[k]),
            "n_matched": int(matched[k]),
            "n_inserted": int(present[k] - matched[k]),
            "frac_inserted": (present[k] - matched[k]) / present[k],
            "median_dist": float(np.median(dists[k])),
        })

    out = pd.DataFrame(rows)
    out.to_csv(args.output, index=False)
    long_path = args.output.replace(".csv", "_per_design.csv")
    pd.DataFrame(per_design).to_csv(long_path, index=False)
    print(f"wrote {long_path}")
    print(f"wrote {args.output}")
    print(f"{len(paths)} design backbones, region lengths "
          f"{min(lengths)}-{max(lengths)}, cutoff {args.cutoff} A")
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
