#!/usr/bin/env python3
"""
haddock3_prepare.py
-------------------
Prepare inputs for HADDOCK3 docking.

Parses the RFDiffusion contig string to identify de novo regions on the
receptor as HADDOCK active residues. Generates ambiguous interaction
restraints (AIRs) and copies PDBs with HADDOCK-friendly names.

Outputs:
    receptor_haddock.pdb  - receptor PDB (renamed copy)
    effector_haddock.pdb  - effector PDB (renamed copy)
    ambig_restraints.tbl  - HADDOCK AIR restraints file
    active_residues.json  - summary of active residues
"""

import argparse
import json
import shutil
import sys


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--receptor", required=True, help="Receptor PDB file")
    parser.add_argument("--effector", required=True, help="Effector PDB file")
    parser.add_argument("--contigs", required=True, help="RFDiffusion contig string")
    parser.add_argument("--receptor-chain", default="A", help="Receptor chain ID")
    parser.add_argument("--effector-chain", default="B", help="Effector chain ID")
    parser.add_argument("--effector-active-residues", default="",
                        help="Comma-separated effector residues for two-sided AIRs "
                             "(e.g. '24,25,26,40,41'). If empty, AIRs use entire effector.")
    return parser.parse_args()


def parse_contig_segments(contigs, rec_chain):
    """
    Parse the receptor block from a contig string into an ordered list of
    (type, start, end) tuples where type is 'fixed' or 'denovo'.

    Thin adapter around :class:`contig_spec.ContigSpec` (Phase 4 Tier 0).
    Fixed segments start with the receptor chain letter (e.g. A1-400).
    De novo segments are bare numbers (e.g. 20-30) representing lengths.
    Commas are tolerated for back-compat.
    """
    import sys as _sys
    from pathlib import Path as _Path
    _sys.path.insert(0, str(_Path(__file__).resolve().parent))
    from contig_spec import ContigSpec, FixedSegment, DeNovoSegment  # noqa: E402

    normalised = contigs.replace(",", " ")
    spec = ContigSpec.from_string(normalised)

    if rec_chain.upper() not in (c.upper() for c in spec.chain_ids):
        sys.exit(f"ERROR: receptor chain '{rec_chain}' not found in contig "
                 f"string '{contigs}'")

    parsed = []
    chain = next(c for c in spec.chains if c.chain_id.upper() == rec_chain.upper())
    for seg in chain.segments:
        if isinstance(seg, FixedSegment):
            parsed.append(("fixed", seg.start, seg.end))
        elif isinstance(seg, DeNovoSegment):
            parsed.append(("denovo", seg.min_len, seg.max_len))
    return parsed


def find_denovo_residues(parsed_segments, receptor_pdb=None, rec_chain="A"):
    """
    Walk parsed segments sequentially to compute output PDB residue numbers
    for de novo regions. De novo tokens specify LENGTHS, not residue numbers.
    The output PDB is numbered contiguously from 1.

    Variable-length segments (min != max) are resolved from the receptor PDB
    total residue count for the receptor chain.
    """
    # Count fixed residues and identify variable-length de novo segments
    total_fixed = 0
    known_denovo = 0
    variable_indices = []

    for i, (seg_type, seg_min, seg_max) in enumerate(parsed_segments):
        if seg_type == "fixed":
            total_fixed += seg_max - seg_min + 1
        elif seg_min == seg_max:
            known_denovo += seg_min
        else:
            variable_indices.append(i)

    # Resolve variable-length de novo segments via PDB residue count
    resolved_lengths = {}
    if variable_indices:
        if receptor_pdb is None:
            sys.exit("ERROR: variable-length de novo segment(s) detected but "
                     "no receptor PDB provided to determine actual length.")
        pdb_residues = set()
        with open(receptor_pdb) as fh:
            for line in fh:
                # Filter by receptor chain to avoid overcounting in PDBs
                # that contain additional chains (e.g. cofactors).
                if line.startswith("ATOM") and line[21] == rec_chain.upper():
                    pdb_residues.add(int(line[22:26].strip()))
        remaining = len(pdb_residues) - total_fixed - known_denovo

        if len(variable_indices) == 1:
            idx = variable_indices[0]
            seg_min, seg_max = parsed_segments[idx][1], parsed_segments[idx][2]
            if not (seg_min <= remaining <= seg_max):
                print(f"WARNING: inferred de novo length {remaining} outside "
                      f"contig range {seg_min}-{seg_max}", file=sys.stderr)
            resolved_lengths[idx] = remaining
        else:
            print(f"WARNING: {len(variable_indices)} variable-length de novo segments; "
                  f"using minimum lengths.", file=sys.stderr)
            for idx in variable_indices:
                resolved_lengths[idx] = parsed_segments[idx][1]

    # Walk segments and assign output residue numbers
    pos = 1
    active_residues = []
    for i, (seg_type, seg_min, seg_max) in enumerate(parsed_segments):
        if seg_type == "fixed":
            pos += seg_max - seg_min + 1
        else:
            length = seg_min if seg_min == seg_max else resolved_lengths[i]
            active_residues.extend(range(pos, pos + length))
            pos += length

    return sorted(active_residues)


def format_ranges(residue_list):
    """Convert a sorted list of residue numbers to compact range strings."""
    if not residue_list:
        return []
    ranges = []
    start = prev = residue_list[0]
    for r in residue_list[1:]:
        if r != prev + 1:
            ranges.append(f"{start}-{prev}" if start != prev else str(start))
            start = r
        prev = r
    ranges.append(f"{start}-{prev}" if start != prev else str(start))
    return ranges


def write_air_restraints(path, active_residues, rec_chain, eff_chain, eff_active_residues=None):
    """
    Write HADDOCK3 ambiguous interaction restraint (AIR) file.

    If eff_active_residues is provided, each receptor active residue is
    restrained to only those effector residues. Otherwise, restrained
    to the entire effector chain.
    """
    with open(path, "w") as f:
        f.write("! Ambiguous Interaction Restraints\n")
        if eff_active_residues:
            f.write(f"! Effector active residues ({len(eff_active_residues)}):\n")
            for i in range(0, len(eff_active_residues), 20):
                chunk = eff_active_residues[i:i+20]
                f.write(f"!   {','.join(str(r) for r in chunk)}\n")
            for resnum in active_residues:
                f.write(f"assign (resid {resnum} and segid {rec_chain})\n")
                f.write("       (\n")
                for i, eff_res in enumerate(eff_active_residues):
                    connector = "or" if i < len(eff_active_residues) - 1 else "  "
                    f.write(f"        (resid {eff_res} and segid {eff_chain}) {connector}\n")
                f.write("       ) 2.0 2.0 2.0\n")
        else:
            for resnum in active_residues:
                f.write(f"assign (resid {resnum} and segid {rec_chain})\n")
                f.write(f"       ((segid {eff_chain})) 2.0 2.0 0.0\n")


def parse_effector_residue_spec(spec):
    """Parse '1-80,90,95-100' into a sorted list of residue numbers."""
    if not spec or not spec.strip():
        return []
    residues = set()
    for token in spec.split(","):
        token = token.strip()
        if not token:
            continue
        if "-" in token:
            start, end = token.split("-", 1)
            residues.update(range(int(start), int(end) + 1))
        else:
            residues.add(int(token))
    return sorted(residues)


def main():
    args = parse_args()

    # ── Parse contig string ──────────────────────────────────────────────
    parsed_segments = parse_contig_segments(args.contigs, args.receptor_chain)
    active_residues = find_denovo_residues(parsed_segments,
                                           receptor_pdb=args.receptor,
                                           rec_chain=args.receptor_chain)

    # Diagnostics
    all_fixed = set()
    n_denovo = 0
    for seg_type, seg_min, seg_max in parsed_segments:
        if seg_type == "fixed":
            all_fixed.update(range(seg_min, seg_max + 1))
        else:
            n_denovo += 1

    print("Contig parse summary:")
    print(f"  Total fixed receptor residues : {len(all_fixed)}")
    print(f"  De novo segments              : {n_denovo}")
    print(f"  Active residues for HADDOCK   : {len(active_residues)}")
    if active_residues:
        print(f"  Active segments               : {', '.join(format_ranges(active_residues))}")

    # ── Parse effector active residues ───────────────────────────────────
    eff_active_residues = parse_effector_residue_spec(args.effector_active_residues)
    if eff_active_residues:
        print(f"  Effector active residues    : {len(eff_active_residues)}")
        print(f"  Effector active segments    : {', '.join(format_ranges(eff_active_residues))}")

    if not active_residues:
        print("ERROR: No active residues identified. Check contig string.", file=sys.stderr)
        sys.exit(1)

    # ── Write restraints ─────────────────────────────────────────────────
    write_air_restraints(
        "ambig_restraints.tbl", active_residues,
        args.receptor_chain, args.effector_chain,
        eff_active_residues if eff_active_residues else None,
    )
    if eff_active_residues:
        print(f"Generated {len(active_residues)} two-sided AIR restraints "
              f"(receptor ↔ {len(eff_active_residues)} effector residues) -> ambig_restraints.tbl")
    else:
        print(f"Generated {len(active_residues)} AIR restraints -> ambig_restraints.tbl")

    # ── Save active residues JSON ────────────────────────────────────────
    with open("active_residues.json", "w") as f:
        json.dump({
            "receptor_chain": args.receptor_chain,
            "active_residues": active_residues,
            "n_active": len(active_residues),
            "all_fixed_count": len(all_fixed),
            "effector_active_residues": eff_active_residues,
            "n_effector_active": len(eff_active_residues),
        }, f, indent=2)

    # ── Copy PDBs ────────────────────────────────────────────────────────
    shutil.copy(args.receptor, "receptor_haddock.pdb")
    shutil.copy(args.effector, "effector_haddock.pdb")
    print(f"Copied {args.receptor} -> receptor_haddock.pdb")
    print(f"Copied {args.effector} -> effector_haddock.pdb")


if __name__ == "__main__":
    main()