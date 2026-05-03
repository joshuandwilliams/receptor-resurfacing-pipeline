#!/usr/bin/env python3
"""
extract_hotspots.py
-------------------
Identify interface residues from a docked complex PDB.

Finds effector residues within a distance cutoff of any receptor heavy atom
and writes the hotspot string in RFDiffusion format (e.g. "B24,B25,B26").

Outputs:
    hotspot_string.txt      - RFDiffusion-format hotspot string
    interface_analysis.json  - full interface residue lists and metadata
"""

import argparse
import json
import sys

from contig_utils import THREE_TO_ONE
from haddock_utils import extract_heavy_atoms, find_interface_residues


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--complex", required=True, help="Docked complex PDB")
    parser.add_argument("--receptor-chain", default="A", help="Receptor chain ID")
    parser.add_argument("--effector-chain", default="B", help="Effector chain ID")
    parser.add_argument("--cutoff", type=float, default=8.0,
                        help="Distance cutoff in Angstroms (default: 8.0)")
    parser.add_argument("--receptor-seq", default=None,
                        help="Reference receptor sequence (one-letter). Used "
                             "to disambiguate chains if HADDOCK has relabelled "
                             "them.  Optional but strongly recommended.")
    parser.add_argument("--effector-seq", default=None,
                        help="Reference effector sequence (one-letter).  Same "
                             "purpose as --receptor-seq.")
    return parser.parse_args()


def _read_chain_sequences(pdb_path):
    """Return {chain_id: one_letter_seq} for every protein chain in the PDB.

    First atom record per (chain, resnum) wins (matches HADDOCK's first-altloc
    convention).  Non-standard residues become 'X'.
    """
    chains_seen = {}  # chain_id -> {resnum: resname}
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            ch = line[21]
            try:
                resnum = int(line[22:26].strip())
            except ValueError:
                continue
            resname = line[17:20].strip()
            chains_seen.setdefault(ch, {})
            if resnum not in chains_seen[ch]:
                chains_seen[ch][resnum] = resname
    return {
        ch: ''.join(THREE_TO_ONE.get(rn, 'X') for _, rn in sorted(d.items()))
        for ch, d in chains_seen.items()
    }


def _best_identity(query, reference):
    """Return the best identity fraction (0-1) between query and reference,
    sliding query over reference at every offset (positive and negative) and
    taking the maximum match-fraction over the overlapping region.

    Adequate for the chain-disambiguation case, where the candidate sequence
    is either identical to a reference (intact chain), terminally truncated
    (HADDOCK trimmed termini), or unrelated (it's the other chain).
    """
    if not query or not reference:
        return 0.0
    q_len, r_len = len(query), len(reference)
    best = 0.0
    # Slide query over reference: offsets where q starts at -q_len+1 .. r_len-1
    for offset in range(-q_len + 1, r_len):
        q_start = max(0, -offset)
        r_start = max(0, offset)
        overlap = min(q_len - q_start, r_len - r_start)
        if overlap < 10:  # ignore tiny overlaps to avoid spurious 100%
            continue
        matches = sum(
            1 for i in range(overlap)
            if query[q_start + i] == reference[r_start + i]
        )
        frac = matches / overlap
        if frac > best:
            best = frac
    return best


def resolve_chains(pdb_path, rec_chain, eff_chain, rec_ref=None, eff_ref=None):
    """
    Identify which chain in the docked complex is the receptor and which is
    the effector.  Strategy:

      1. If both expected chain IDs (rec_chain, eff_chain) are present in
         the PDB, use them directly.  This is the 99% case.
      2. Otherwise, fall back to sequence-identity matching: read every
         chain's sequence from the PDB and assign each chain to whichever
         reference (rec_ref or eff_ref) it matches best.  Requires
         --receptor-seq and --effector-seq to be supplied; fails loudly if
         not.
      3. The match must be >=80% identity, otherwise fail loudly — chain
         identity should be unambiguous.
    """
    atoms = extract_heavy_atoms(pdb_path, (rec_chain, eff_chain))
    if atoms[rec_chain] and atoms[eff_chain]:
        return atoms[rec_chain], atoms[eff_chain], rec_chain, eff_chain

    # Expected chain IDs missing — fall back to sequence matching.
    all_chains = sorted(_read_chain_sequences(pdb_path).keys())
    print(f"WARNING: Expected chains {rec_chain},{eff_chain} not both present "
          f"in {pdb_path}; found {all_chains}.  Falling back to "
          f"sequence-identity matching.")

    if rec_ref is None or eff_ref is None:
        sys.exit("ERROR: cannot disambiguate chains — expected chain IDs "
                 "missing and reference sequences not provided.  Pass "
                 "--receptor-seq and --effector-seq to enable fallback.")

    chain_seqs = _read_chain_sequences(pdb_path)
    if len(chain_seqs) < 2:
        sys.exit(f"ERROR: need at least two chains in {pdb_path}; found {len(chain_seqs)}.")

    # Score every chain against both references.
    scores = {}
    for ch, seq in chain_seqs.items():
        scores[ch] = (_best_identity(seq, rec_ref), _best_identity(seq, eff_ref))
        print(f"  chain {ch} ({len(seq)} aa): receptor-id={scores[ch][0]:.2f} "
              f"effector-id={scores[ch][1]:.2f}")

    # Pick the chain best matching receptor reference, then the chain best
    # matching effector reference (excluding the one already chosen).
    rec_pick = max(scores, key=lambda c: scores[c][0])
    rec_id = scores[rec_pick][0]
    remaining = {c: s for c, s in scores.items() if c != rec_pick}
    eff_pick = max(remaining, key=lambda c: remaining[c][1])
    eff_id = remaining[eff_pick][1]

    if rec_id < 0.8 or eff_id < 0.8:
        sys.exit(f"ERROR: chain disambiguation failed.  Best receptor match "
                 f"{rec_pick} at {rec_id:.0%} identity; best effector match "
                 f"{eff_pick} at {eff_id:.0%} identity.  Both must be >=80%.  "
                 f"Check --receptor-seq / --effector-seq are correct and the "
                 f"docked complex contains the expected chains.")

    print(f"Sequence-based assignment: receptor={rec_pick} ({rec_id:.0%}), "
          f"effector={eff_pick} ({eff_id:.0%})")
    fallback = extract_heavy_atoms(pdb_path, (rec_pick, eff_pick))
    return fallback[rec_pick], fallback[eff_pick], rec_pick, eff_pick


def main():
    args = parse_args()

    rec_atoms, eff_atoms, _, eff_chain = resolve_chains(
        args.complex, args.receptor_chain, args.effector_chain,
        rec_ref=args.receptor_seq, eff_ref=args.effector_seq,
    )

    interface_eff, interface_rec = find_interface_residues(
        rec_atoms, eff_atoms, args.cutoff
    )

    if interface_eff:
        print(f"Interface effector residues ({len(interface_eff)}): "
              f"{interface_eff[0]}-{interface_eff[-1]}")
    else:
        print("WARNING: No effector interface residues found")

    if interface_rec:
        print(f"Interface receptor residues ({len(interface_rec)}): "
              f"{interface_rec[0]}-{interface_rec[-1]}")
    else:
        print("WARNING: No receptor interface residues found")

    hotspot_str = ",".join(f"{eff_chain}{r}" for r in interface_eff)
    print(f"Hotspot string: {hotspot_str}")

    with open("hotspot_string.txt", "w") as f:
        f.write(hotspot_str)

    with open("interface_analysis.json", "w") as f:
        json.dump({
            "effector_interface_residues": interface_eff,
            "receptor_interface_residues": interface_rec,
            "n_effector_contacts": len(interface_eff),
            "n_receptor_contacts": len(interface_rec),
            "contact_cutoff": args.cutoff,
            "hotspot_string": hotspot_str,
        }, f, indent=2)


if __name__ == "__main__":
    main()