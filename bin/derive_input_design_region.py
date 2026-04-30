#!/usr/bin/env python3
"""
derive_input_design_region.py
-----------------------------
Task 6 revision (v7) · Derive design-region and true-interface indices
for the INPUT STRUCTURE — the receptor+effector complex going into
RFDiffusion — rather than for an RFDiffusion output design.

Used only by the Task 6 negative-control path: the scrambled and polyA
controls run against the INPUT structure (answering "if RFDiffusion did
a terrible job, what metrics would those sequences get?") rather than
against individual RFDiffusion designs.

Why a new script
================
``derive_design_region.py`` reads rfdiffusion_metrics.json which is a
per-design output.  The input structure has no such file — its "design
region" is whatever the USER's contigs string marks for de novo
rebuilding.  Likewise ``derive_true_interface.py`` reads
rfdiffusion_metrics.json; for the input structure we compute contacts
directly off the input PDB.

Outputs (matching the formats the downstream pipeline expects)
==============================================================
--design-region-output
    1-based comma-separated positional indices on the receptor chain,
    same format as derive_design_region.py produces.

--true-interface-output
    0-based, one per line, same format as derive_true_interface.py
    produces.

Design-region semantics
========================
The RFDiffusion contigs string ``A1-50 5-15 A60-100 B`` means:
  - Anchor receptor positions 1..50 in the input PDB (kept verbatim)
  - Build 5..15 de novo receptor residues next
  - Anchor receptor positions 60..100 (kept verbatim)
  - Effector chain B (any length)

On the INPUT PDB's receptor chain, positions 1..50 and 60..100 are
anchors.  Positions 51..59 (the gap between anchor 1 end and anchor 2
start) are the part of the NATIVE receptor that RFDiffusion is being
told to REPLACE with de novo chemistry.  Those positions are what
we call the "design region on the input structure" — scrambling or
polyA-ing THEM answers "what if we wreck the native target sequence
where RFDiffusion was going to rebuild?"

True-interface semantics
========================
Receptor residues with at least one heavy atom within ``contact_cutoff``
(default 5.0 Å) of any effector heavy atom in the input complex PDB.
Uses the same find_contact_residues_heavy helper from compute_metrics.py
that the per-design path uses.

Usage
=====
  python derive_input_design_region.py \\
      --input-pdb rfdiffusion_input.pdb \\
      --contigs "A1-50 5-15 A60-100 B" \\
      --receptor-chain A \\
      --effector-chain B \\
      --design-region-output input_design_region.txt \\
      --true-interface-output input_true_interface.txt
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import List, Tuple


# ── Contig parsing ──────────────────────────────────────────────────

def _parse_contigs(
    contigs_str: str,
    receptor_chain: str,
    effector_chain: str,
) -> List[Tuple[str, int, int, bool]]:
    """
    Parse an RFDiffusion contigs string into a flat sequence of segments.

    Returns a list of (chain_id, start, end, is_denovo) tuples in the
    order they appear in the contigs.  is_denovo is True for numeric-only
    segments (design region — positions TO BE built, not positions IN
    the input), False for chain-letter-prefixed anchor segments.

    The ``start`` and ``end`` are:
      - for anchor segments: the receptor/effector position RANGE that
        must be present in the input PDB (inclusive);
      - for de novo segments: the min/max of the length sampling range
        (these are NOT positions — they're lengths that RFDiffusion
        samples from).
    """
    # Strip slashes — contigs blocks can be grouped with "/" but for our
    # purposes every segment is independent.
    normalised = contigs_str.replace(",", " ").replace("/", " ")
    tokens = [t.strip() for t in normalised.split() if t.strip()]

    out: List[Tuple[str, int, int, bool]] = []
    for tok in tokens:
        # Pure effector-chain tokens: "B" alone.
        if tok.upper() == effector_chain.upper():
            out.append((effector_chain, 0, 0, False))
            continue
        # Chain-letter-prefixed anchor segment, e.g. "A1-50" or "A5".
        if tok[0].isalpha():
            chain = tok[0]
            rest = tok[1:]
            if not rest:
                raise ValueError(f"Malformed contig token: {tok!r}")
            if "-" in rest:
                try:
                    s, e = [int(x) for x in rest.split("-")]
                except ValueError:
                    raise ValueError(f"Malformed contig token: {tok!r}")
            else:
                s = e = int(rest)
            out.append((chain, s, e, False))
            continue
        # Numeric-only de novo / flexible segment, e.g. "5-15" or "10".
        if "-" in tok:
            try:
                s, e = [int(x) for x in tok.split("-")]
            except ValueError:
                raise ValueError(f"Malformed contig token: {tok!r}")
        else:
            try:
                s = e = int(tok)
            except ValueError:
                raise ValueError(f"Malformed contig token: {tok!r}")
        out.append(("_denovo", s, e, True))
    return out


# ── Input-structure receptor position walker ──────────────────────────

def _receptor_positions_by_contigs(
    segments: List[Tuple[str, int, int, bool]],
    receptor_chain: str,
) -> List[Tuple[int, int]]:
    """
    Walk the contigs and return (anchor_start, anchor_end) pairs for
    consecutive receptor anchor segments — the RECEPTOR's position
    ranges that are preserved.  The GAPS between consecutive anchors
    are the positions we treat as the "design region on the input
    receptor."

    Segments for the effector chain and for de novo blocks are skipped
    in this walk — we only care about the receptor anchors' positions.
    """
    anchors: List[Tuple[int, int]] = []
    for chain, s, e, is_denovo in segments:
        if is_denovo:
            continue
        if chain.upper() != receptor_chain.upper():
            continue
        anchors.append((s, e))
    # Ensure anchors are sorted by start — malformed contigs that list
    # anchors out of order would otherwise produce nonsense gaps.
    anchors.sort()
    return anchors


def _design_region_positions_1b(
    anchors: List[Tuple[int, int]],
    receptor_len: int,
) -> List[int]:
    """
    Given the anchor ranges on the input receptor (1-based inclusive,
    as parsed from the contigs) and the total receptor length, return
    the 1-based positions that fall OUTSIDE any anchor — these are the
    positions the design region would replace.

    For the classic ``A1-50 5-15 A60-100`` contig on a 100-residue
    receptor, anchors are [(1,50), (60,100)] and the gap positions are
    [51, 52, 53, ..., 59].
    """
    if not anchors:
        # Edge case: no receptor anchors at all (unusual but possible —
        # entire receptor is being rebuilt).  Return every position as
        # design-region.
        return list(range(1, receptor_len + 1))

    anchor_set = set()
    for s, e in anchors:
        for p in range(s, e + 1):
            anchor_set.add(p)

    return [p for p in range(1, receptor_len + 1) if p not in anchor_set]


# ── True interface via heavy-atom contacts ────────────────────────────

def _compute_true_interface_0b(
    input_pdb: Path,
    receptor_chain: str,
    effector_chain: str,
    contact_cutoff: float,
    expected_rec_seq: str,
) -> List[int]:
    """
    Compute receptor residues at the true interface via heavy-atom
    contact detection on the input complex PDB.  Returns 0-based
    positional indices on the receptor chain.

    Re-uses find_contact_residues_heavy from compute_metrics.py to
    match the per-design path's contact-detection semantics exactly.
    """
    # Use the same import pattern as boltz2_negative_steering.py line 2249
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from compute_metrics import find_contact_residues_heavy  # noqa: E402

    contacts = find_contact_residues_heavy(
        str(input_pdb),
        receptor_chain,
        effector_chain,
        contact_cutoff,
        expected_rec_seq=expected_rec_seq,
    )
    # contacts is a list of (positional_index, closest_distance) tuples.
    return sorted(i for i, _ in contacts)


# ── Receptor sequence extraction ──────────────────────────────────────

def _extract_receptor_sequence(
    input_pdb: Path,
    receptor_chain: str,
) -> str:
    """Extract the receptor chain sequence using the existing helper."""
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from boltz2_negative_steering import get_chain_sequence  # noqa: E402
    return get_chain_sequence(input_pdb, receptor_chain)


# ── Output writers ────────────────────────────────────────────────────

def _write_design_region(
    output_path: Path,
    positions_1b: List[int],
    contigs_str: str,
    input_pdb: Path,
    receptor_chain: str,
) -> None:
    """
    Emit the design-region file in the exact format that
    derive_design_region.py writes, so downstream tooling reads it
    without caring whether it came from an RFDiffusion output or the
    input structure.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as f:
        f.write(f"# derived from input structure {input_pdb} + contigs\n")
        f.write(f"# contigs: {contigs_str}\n")
        f.write(f"# receptor chain: {receptor_chain}\n")
        f.write(f"# design region size: {len(positions_1b)} residues\n")
        f.write(f"# coordinate system: 1-based positional indices on the receptor\n")
        f.write(f"#                    chain, matching input PDB sequence order\n")
        if positions_1b:
            f.write(",".join(str(p) for p in positions_1b) + "\n")
        else:
            f.write("\n")


def _write_true_interface(
    output_path: Path,
    positions_0b: List[int],
    input_pdb: Path,
    receptor_chain: str,
    effector_chain: str,
    contact_cutoff: float,
) -> None:
    """
    Emit the true-interface file in the exact format that
    derive_true_interface.py writes.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as f:
        f.write(f"# derived from input structure {input_pdb}\n")
        f.write(f"# heavy-atom contacts <= {contact_cutoff} Angstroms\n")
        f.write(f"# receptor chain: {receptor_chain}\n")
        f.write(f"# effector chain: {effector_chain}\n")
        f.write(f"# n contacts: {len(positions_0b)}\n")
        f.write(f"# coordinate system: 0-based positional indices on the "
                f"receptor chain,\n")
        f.write(f"#                    matching input PDB sequence order\n")
        for p in positions_0b:
            f.write(f"{p}\n")


# ── Main ──────────────────────────────────────────────────────────────

def main() -> int:
    ap = argparse.ArgumentParser(
        description="Derive design-region and true-interface indices "
                    "for the INPUT STRUCTURE (used by Task 6 negative "
                    "controls).",
    )
    ap.add_argument("--input-pdb", required=True, type=Path,
                    help="Input complex PDB going into RFDiffusion.")
    ap.add_argument("--contigs", required=True, type=str,
                    help="RFDiffusion contigs string (e.g. 'A1-50 5-15 "
                         "A60-100 B').")
    ap.add_argument("--receptor-chain", default="A",
                    help="Receptor chain ID (default A).")
    ap.add_argument("--effector-chain", default="B",
                    help="Effector chain ID (default B).")
    ap.add_argument("--contact-cutoff", type=float, default=5.0,
                    help="Heavy-atom contact cutoff in Angstroms "
                         "(default 5.0 — matches the per-design path).")
    ap.add_argument("--design-region-output", required=True, type=Path,
                    help="Output file for 1-based design-region positions.")
    ap.add_argument("--true-interface-output", required=True, type=Path,
                    help="Output file for 0-based true-interface positions.")
    args = ap.parse_args()

    # ── Extract receptor sequence from input PDB ───────────────────
    try:
        rec_seq = _extract_receptor_sequence(args.input_pdb, args.receptor_chain)
    except Exception as e:  # noqa: BLE001
        print(f"ERROR: receptor sequence extraction failed: {e}",
              file=sys.stderr)
        return 2
    if not rec_seq:
        print(f"ERROR: receptor chain '{args.receptor_chain}' empty in "
              f"{args.input_pdb}", file=sys.stderr)
        return 2
    receptor_len = len(rec_seq)
    print(f"Receptor: {receptor_len} residues, chain {args.receptor_chain}")

    # ── Parse contigs → extract anchor ranges → compute design region ─
    try:
        segments = _parse_contigs(
            args.contigs, args.receptor_chain, args.effector_chain,
        )
    except ValueError as e:
        print(f"ERROR: contigs parse failed: {e}", file=sys.stderr)
        return 2
    anchors = _receptor_positions_by_contigs(segments, args.receptor_chain)
    design_region_1b = _design_region_positions_1b(anchors, receptor_len)
    print(f"Contigs: {args.contigs}")
    print(f"Receptor anchors: {anchors}")
    print(f"Design region (input receptor): {len(design_region_1b)} positions")
    if design_region_1b:
        print(f"  first: {design_region_1b[0]}  last: {design_region_1b[-1]}  "
              f"(1-based)")

    # ── Compute true interface via heavy-atom contacts ──────────────
    try:
        true_iface_0b = _compute_true_interface_0b(
            args.input_pdb,
            args.receptor_chain,
            args.effector_chain,
            args.contact_cutoff,
            expected_rec_seq=rec_seq,
        )
    except Exception as e:  # noqa: BLE001
        print(f"ERROR: true-interface computation failed: {e}",
              file=sys.stderr)
        return 2
    print(f"True interface: {len(true_iface_0b)} receptor residues "
          f"within {args.contact_cutoff} Å of effector")

    # ── Write outputs ──────────────────────────────────────────────
    _write_design_region(
        args.design_region_output, design_region_1b,
        args.contigs, args.input_pdb, args.receptor_chain,
    )
    _write_true_interface(
        args.true_interface_output, true_iface_0b,
        args.input_pdb, args.receptor_chain, args.effector_chain,
        args.contact_cutoff,
    )
    print(f"Wrote {args.design_region_output}")
    print(f"Wrote {args.true_interface_output}")

    # Sanity: the design region and true interface should NOT overlap
    # completely — if they do, the contigs are telling RFDiffusion to
    # rebuild the exact binding interface, which would defeat the
    # purpose of the whole design philosophy.  Warn advisorily.
    design_set = set(p - 1 for p in design_region_1b)   # 0-based
    true_set = set(true_iface_0b)
    overlap = design_set & true_set
    if design_region_1b and true_iface_0b:
        overlap_frac = len(overlap) / len(true_set) if true_set else 0
        if overlap_frac > 0.8:
            print(f"WARNING: true interface is {overlap_frac:.0%} inside "
                  f"the design region — controls will fully destroy "
                  f"the binding interface.  Make sure this is intended.",
                  file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())