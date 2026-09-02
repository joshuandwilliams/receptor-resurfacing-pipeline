#!/usr/bin/env python3
"""Re-align design region 1 so the sequences are flush at BOTH ends.

ClustalW puts every gap at the right-hand end, so the sequences are left
aligned and the longest design's last three residues hang off unaligned. That
is wrong here: design region 1 is a beta hairpin anchored to fixed scaffold at
both ends, so its last residues are the second strand and belong in a column
with everybody else's. Every one of these sequences has isoleucine two from
the end, which is the conserved strand position, and ClustalW leaves it
scattered across three different columns.

No gap penalty fixes it. Sweeping GAPOPEN over 10 to 30 and GAPEXT over 0.1 to
0.5, with and without -ENDGAPS, leaves exactly one of eleven sequences flush at
the right in every one of the 24 combinations. That is a property of the
algorithm rather than the parameters: terminal gaps are free in ClustalW's
model, so packing them all at one end is always optimal or tied, and
-ENDGAPS turns the separation penalty off rather than on. -TERMINALGAP applies
only to structure-guided alignment and needs masks this has no way to supply.

So the alignment is done here instead, with an aligner whose end gaps can be
charged. Each sequence is aligned to the longest one, which is the reference
that spans every column, with end gaps penalised at -8 so gaps are forced
inside the sequence rather than onto its ends. Open -20 and extend -0.5 were
chosen by a grid search as the setting that leaves all eleven sequences flush
at both ends AND gives each exactly one internal gap, rather than breaking the
loop into two or three fragments.

Region 2 is a fixed six residues in every design, so it needs no realignment
and its ClustalW output is passed through unchanged.

Usage:
    anchor_alignment.py [--out region1_anchored.aln]
"""

import argparse
import os
import warnings

warnings.filterwarnings("ignore")

from Bio import Align
from Bio.Align import substitution_matrices

HERE = os.path.dirname(os.path.abspath(__file__))

END_GAP = -8.0
OPEN_GAP = -20.0
EXTEND_GAP = -0.5


def read_aln(path):
    seqs = {}
    for line in open(path):
        if not line.strip() or line.startswith(("CLUSTAL", " ")):
            continue
        parts = line.split()
        if len(parts) == 2:
            seqs[parts[0]] = seqs.get(parts[0], "") + parts[1]
    return seqs


def build_aligner():
    al = Align.PairwiseAligner()
    al.mode = "global"
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score = OPEN_GAP
    al.extend_gap_score = EXTEND_GAP
    # Charging the end gaps is the whole point. Left at 0 they are free and the
    # aligner behaves like ClustalW.
    al.end_insertion_score = END_GAP
    al.end_deletion_score = END_GAP
    return al


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--aln", default=os.path.join(HERE, "region1.aln"))
    ap.add_argument("--out", default=os.path.join(HERE, "region1_anchored.aln"))
    args = ap.parse_args()

    ungapped = {k: v.replace("-", "") for k, v in read_aln(args.aln).items()}
    ref = max(ungapped, key=lambda k: len(ungapped[k]))
    al = build_aligner()

    out = {}
    for name, seq in ungapped.items():
        out[name] = seq if name == ref else str(al.align(ungapped[ref], seq)[0][1])

    width = len(out[ref])
    flush = sum(1 for v in out.values()
                if not v.startswith("-") and not v.endswith("-"))
    runs = sum(v.strip("-").count("-") > 0 for v in out.values())
    print(f"reference {ref}, {len(out)} sequences, width {width}")
    print(f"flush at both ends: {flush}/{len(out)}")

    with open(args.out, "w") as fh:
        fh.write("CLUSTAL 2.1 multiple sequence alignment\n\n\n")
        for name, seq in out.items():
            fh.write(f"{name:<18s} {seq}\n")
    print(f"wrote {args.out}")
    for name, seq in out.items():
        print(f"  {name:<15s} {seq}")


if __name__ == "__main__":
    main()
