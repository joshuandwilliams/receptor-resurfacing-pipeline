#!/usr/bin/env python3
"""Where each designed residue sits in the hairpin, as a distance rather than an index.

Design region 1 cannot be aligned across the campaign. RFdiffusion varied its
length from 10 to 20 residues, so position 7 of a short design and position 7 of
a long one are at different places in the fold, which rules out per-position
profiles, logos and any alignment view.

Distance sidesteps that. Design region 1 is a beta hairpin anchored to fixed
scaffold at residue 32 and again at the residue immediately after the region, and
those two anchors are in the same place in every design because the scaffold was
never redesigned. The midpoint between their CA atoms is therefore a fixed
reference at the base of the hairpin, and the distance from it to a residue's CA
says how far out towards the tip that residue sits. It is a continuous coordinate,
so it needs no correspondence between designs and no superposition: each design is
measured in its own frame against its own anchors.

Depth answers a different question from position. Position asks which residue of
another design a residue lines up with, which has no answer here. Depth asks how
far from the scaffold it is, which every design can answer on the same scale.

Frequencies downstream are per SEQUENCE, n = 128, since the amino acid is what
ProteinMPNN chose; the depth comes from the backbone, so the 64 backbones each
contribute their two sequences at identical depths.

The wild type goes through the same calculation from the same reference structure,
so its residues land on the same axis.

Writes design, sequence, region length, residue number, index within the region,
amino acid, depth in angstroms, and the DSSP assignment for that residue.

Usage:
    depth_extract.py [--output loop_depth.csv]
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
RUN = os.path.join(
    HERE, "..", "..", "runs", "crystal_full_test_contig", "results")
INPUTS = os.path.join(HERE, "..", "..", "inputs")
DSSP = os.path.join(HERE, "..", "01-design-region-backbones",
                    "design_secondary_structure.csv")

# The scaffold either side of design region 1 is fixed in every design: 32
# residues before it and 33 after, so a backbone of N residues carries a region
# of N - 65 and the region runs from 33 to 32 + (N - 65).
N_BEFORE = 32
N_AFTER = 33


def alpha_carbons(path):
    model = PDBParser(QUIET=True).get_structure("s", path)[0]
    return {r.id[1]: np.asarray(r["CA"].coord, dtype=float)
            for r in model["A"] if r.id[0] == " " and "CA" in r}


def depths(path):
    """Distance from the anchor midpoint to each design region 1 CA."""
    ca = alpha_carbons(path)
    length = max(ca) - N_BEFORE - N_AFTER
    # One anchor each side. Their midpoint sits at the base of the hairpin, on
    # the scaffold, so it does not move when the loop is lengthened.
    anchor = (ca[N_BEFORE] + ca[N_BEFORE + 1 + length]) / 2.0
    return length, {n: float(np.linalg.norm(ca[n] - anchor))
                    for n in range(N_BEFORE + 1, N_BEFORE + 1 + length)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", default=os.path.join(HERE, "loop_depth.csv"))
    args = ap.parse_args()

    dssp = pd.read_csv(DSSP)
    ss = {(r.design, r.resnum): r.ss for r in dssp.itertuples()}

    meta = pd.read_csv(os.path.join(
        RUN, "sequences", "scored_metadata.csv"))
    meta["region1"] = meta.designed_residues.str.split("|").str[0]
    sequences = {}
    for r in meta.itertuples():
        sequences.setdefault(f"design_{r.design}", []).append((r.seq, r.region1))

    # The wild type is one structure with one sequence, read off the same file
    # the designs were built from.
    native = str(meta.native_residues.iloc[0]).split("|")[0]
    jobs = [(p, os.path.basename(p)[:-4], sequences.get(os.path.basename(p)[:-4], []))
            for p in sorted(glob.glob(os.path.join(RUN, "rfdiffusion", "design_*.pdb")))]
    jobs.append((os.path.join(INPUTS, "pikp1_avrpikf_complex_ab.pdb"), "wild-type",
                 [("wt", native)]))

    rows = []
    for path, design, seqs in jobs:
        length, depth = depths(path)
        for seq_id, seq in seqs:
            if len(seq) != length:
                raise ValueError(
                    f"{design} seq {seq_id}: sequence is {len(seq)} residues but "
                    f"the backbone region is {length}")
            for k, aa in enumerate(seq):
                resnum = N_BEFORE + 1 + k
                rows.append({"design": design, "seq": seq_id,
                             "region_length": length, "resnum": resnum,
                             "index_in_region": k + 1, "aa": aa,
                             "depth": depth[resnum],
                             "ss": ss.get((design, resnum), "unknown")})

    out = pd.DataFrame(rows)
    out.to_csv(args.output, index=False)
    designs = out[out.design != "wild-type"]
    print(f"wrote {args.output}")
    print(f"{designs.design.nunique()} backbones, {len(designs.groupby(['design','seq']))} "
          f"sequences, {len(designs)} residues")
    print(f"region length {designs.region_length.min()}-{designs.region_length.max()}, "
          f"depth {designs.depth.min():.1f}-{designs.depth.max():.1f} A")
    print("\ndepth by DSSP assignment:")
    print(designs.groupby("ss").depth.describe()[["count", "min", "50%", "max"]]
          .round(1).to_string())


if __name__ == "__main__":
    main()
