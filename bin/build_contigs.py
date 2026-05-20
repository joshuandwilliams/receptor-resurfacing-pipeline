#!/usr/bin/env python3
"""
build_contigs.py
----------------
Update the RFDiffusion contig string after the docking step.

Given the user's original contig string and the docked complex PDB:
  - Detects renumbering + pLDDT trimming offsets
  - Replaces bare effector chain references with actual PDB length
  - Adjusts receptor residue numbers for the (possibly renumbered) complex
  - Extracts receptor/effector sequences from the complex

Outputs:
    updated_contigs.txt   - corrected contig string for RFDiffusion
    updated_params.json   - receptor/effector sequences, lengths, contigs
    rfdiffusion_input.pdb - copy of the complex PDB for RFDiffusion input
"""

import argparse
import json
import os
import shutil

from contig_utils import THREE_TO_ONE


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--complex", required=True, help="Docked complex PDB from the pose solver")
    parser.add_argument("--receptor-chain", default="A", help="Receptor chain ID")
    parser.add_argument("--effector-chain", default="B", help="Effector chain ID")
    parser.add_argument("--contigs", required=True, help="Original user contig string")
    parser.add_argument("--rec-trim-mapping", default=None,
                        help="JSON trim mapping for receptor (optional)")
    parser.add_argument("--eff-trim-mapping", default=None,
                        help="JSON trim mapping for effector (optional)")
    return parser.parse_args()


def load_mapping(path):
    """Load a trim mapping JSON. Returns None if absent, empty, or placeholder ({})."""
    if path and os.path.exists(path) and os.path.getsize(path) > 2:
        with open(path) as f:
            return json.load(f)
    return None


def extract_chain(pdb_path, chain_id):
    """Return (sequence_str, [resnum, ...]) for a given chain."""
    residues = {}
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM") or line[21] != chain_id:
                continue
            resnum = int(line[22:26].strip())
            resname = line[17:20].strip()
            if resnum not in residues:
                residues[resnum] = resname
    resnums = sorted(residues.keys())
    seq = "".join(THREE_TO_ONE.get(residues[r], "X") for r in resnums)
    return seq, resnums


def resolve_chains(pdb_path, rec_chain, eff_chain):
    """
    Extract both chains, falling back to the first two found alphabetically
    if expected chains are missing.
    """
    rec_seq, rec_resnums = extract_chain(pdb_path, rec_chain)
    eff_seq, eff_resnums = extract_chain(pdb_path, eff_chain)

    if not rec_seq or not eff_seq:
        all_chains = set()
        with open(pdb_path) as f:
            for line in f:
                if line.startswith("ATOM"):
                    all_chains.add(line[21])
        all_chains = sorted(all_chains)
        if len(all_chains) >= 2:
            if not rec_seq:
                rec_seq, rec_resnums = extract_chain(pdb_path, all_chains[0])
                print(f"Using chain {all_chains[0]} as receptor")
            if not eff_seq:
                eff_seq, eff_resnums = extract_chain(pdb_path, all_chains[1])
                print(f"Using chain {all_chains[1]} as effector")

    return rec_seq, rec_resnums, eff_seq, eff_resnums


def parse_contig_receptor_refs(contigs, rec_chain):
    """
    Extract all receptor residue numbers referenced in fixed segments.
    Returns sorted boundary residue numbers, e.g. [1, 400, 426, 435].
    """
    blocks = contigs.replace(",", " ").split()
    refs = []
    for block in blocks:
        if not any(seg.strip() and seg.strip()[0].isalpha()
                   and seg.strip()[0].upper() == rec_chain.upper()
                   for seg in block.split("/")):
            continue
        for seg in block.split("/"):
            seg = seg.strip()
            if not seg or not (seg[0].isalpha() and seg[0].upper() == rec_chain.upper()):
                continue
            rest = seg[1:]
            if "-" in rest:
                parts = rest.split("-")
                refs.extend([int(parts[0]), int(parts[1])])
            elif rest:
                refs.append(int(rest))
    return sorted(refs)


def compute_offset(contigs, rec_chain, rec_resnums, rec_map):
    """
    Compute the total offset to translate contig residue numbers to complex
    PDB numbering, combining renumbering and pLDDT trimming.

    The offset is (contig_last_ref - pdb_last_residue), clamped >= 0, plus
    any N-terminal residues removed by pLDDT trimming.
    """
    contig_refs = parse_contig_receptor_refs(contigs, rec_chain)
    if not contig_refs or not rec_resnums:
        return 0

    contig_last = contig_refs[-1]
    pdb_last = rec_resnums[-1]
    renumber_offset = max(0, contig_last - pdb_last)

    trim_offset = 0
    if rec_map:
        trim_offset = rec_map.get("n_trimmed_nterm", 0)

    total_offset = renumber_offset + trim_offset

    if renumber_offset > 0:
        print(f"Detected renumbering offset: {renumber_offset}")
        print(f"  Contig references residues {contig_refs[0]}-{contig_last}")
        print(f"  PDB receptor spans residues {rec_resnums[0]}-{pdb_last}")
    if trim_offset > 0:
        print(f"pLDDT trim offset: {trim_offset} residues removed from N-term")
    if total_offset > 0:
        print(f"Total offset: {total_offset}")

    return total_offset


def update_contigs(user_contigs, rec_chain, eff_chain, eff_len, rec_len, rec_offset):
    """
    Update the contig string: replace effector blocks with actual length,
    adjust receptor residue numbers by offset, clamp to receptor length.
    """
    blocks = user_contigs.replace(",", " ").split()
    new_blocks = []

    for block in blocks:
        # Bare effector chain reference
        if block.strip().upper() == eff_chain.upper() or (
            len(block) > 1 and block[0].upper() == eff_chain.upper()
        ):
            new_blocks.append(f"{eff_chain}1-{eff_len}")
            continue

        # Receptor block — adjust residue numbers
        segments = block.split("/")
        new_segs = []
        for seg in segments:
            seg = seg.strip()
            if not seg or seg == "0":
                new_segs.append(seg)
                continue
            if seg[0].isalpha() and seg[0].upper() == rec_chain.upper():
                chain = seg[0]
                rest = seg[1:]
                if "-" in rest:
                    parts = rest.split("-")
                    old_start, old_end = int(parts[0]), int(parts[1])
                else:
                    old_start = old_end = int(rest)
                new_start = max(1, min(old_start - rec_offset, rec_len))
                new_end = max(1, min(old_end - rec_offset, rec_len))
                new_segs.append(f"{chain}{new_start}-{new_end}")
            else:
                new_segs.append(seg)
        new_blocks.append("/".join(new_segs))

    return " ".join(new_blocks)


def main():
    args = parse_args()
    rec_map = load_mapping(args.rec_trim_mapping)

    rec_seq, rec_resnums, eff_seq, eff_resnums = resolve_chains(
        args.complex, args.receptor_chain, args.effector_chain
    )

    rec_len, eff_len = len(rec_seq), len(eff_seq)
    print(f"Receptor: {rec_len} residues (PDB range {rec_resnums[0]}-{rec_resnums[-1]})")
    print(f"Effector: {eff_len} residues (PDB range {eff_resnums[0]}-{eff_resnums[-1]})")

    rec_offset = compute_offset(args.contigs, args.receptor_chain, rec_resnums, rec_map)

    updated_contigs = update_contigs(
        args.contigs, args.receptor_chain, args.effector_chain,
        eff_len, rec_len, rec_offset,
    )
    print(f"Original contigs : {args.contigs}")
    print(f"Updated contigs  : {updated_contigs}")

    with open("updated_contigs.txt", "w") as f:
        f.write(updated_contigs)

    shutil.copy(args.complex, "rfdiffusion_input.pdb")

    with open("updated_params.json", "w") as f:
        json.dump({
            "receptor_seq": rec_seq,
            "effector_seq": eff_seq,
            "receptor_len": rec_len,
            "effector_len": eff_len,
            "contigs": updated_contigs,
            "receptor_start_pdb": rec_resnums[0] if rec_resnums else 1,
            "offset_applied": rec_offset,
        }, f, indent=2)


if __name__ == "__main__":
    main()