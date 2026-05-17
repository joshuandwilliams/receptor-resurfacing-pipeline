#!/usr/bin/env python3
"""
haddock3_prepare.py
-------------------
Prepare inputs for HADDOCK3 docking.

Per Session 7 grill-me (notes/design_audit.md Q131, Q133, Q137, Q138,
Q142, Q147) this script no longer parses the RFDiffusion contig string
to derive HADDOCK active residues.  At HADDOCK time the contig is
abstract (the de novo regions don't exist on the input PDB yet), so
restraints come from explicit user params instead.

Two restraint files may be written:

- ``ambig_restraints.tbl`` — ambiguous interaction restraints (AIRs).
  Generated from ``--receptor-active-residues`` (and optionally
  ``--effector-active-residues``).  AND-over-receptor / OR-over-effector
  structure with 50% nrest tolerance (HADDOCK's [airs] block treats
  this as "at least 50% must be satisfied").  Soft constraint.

- ``unambig_restraints.tbl`` — unambiguous distance restraints.
  Generated from ``--contact-pairs`` (format ``"A25-C42 A13-C94"``).
  Each pair becomes one CA-CA distance restraint with the global
  ``--pair-distance`` target/tolerance triple.  Hard pin.

Either or both restraint sets may be empty; HADDOCK reads whichever
files exist.  The pipeline validator (validate_params.py) enforces that
at least one of contact pairs or receptor active residues is supplied
in Branch A.

Chain-ID contract (per A142):
- Caller passes ``--receptor-chain`` / ``--effector-chain`` referring to
  the chain IDs in the user's INPUT PDBs.
- This script writes ``receptor_haddock.pdb`` / ``effector_haddock.pdb``
  with the chain ID FORCED to A (receptor) and B (effector) regardless
  of input chain letters.  Downstream HADDOCK + post-processing assume
  A/B.

Outputs:
    receptor_haddock.pdb     - receptor PDB, chain forced to A
    effector_haddock.pdb     - effector PDB, chain forced to B
    ambig_restraints.tbl     - AIRs from active residues (may be empty)
    unambig_restraints.tbl   - distance restraints from pairs (may be empty)
    restraints_summary.json  - what was written and from which inputs
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from pathlib import Path

# Local imports — share parsers with haddock_run.HaddockRun so the user-
# facing string formats round-trip identically.
_BIN_DIR = os.path.dirname(os.path.abspath(__file__))
if _BIN_DIR not in sys.path:
    sys.path.insert(0, _BIN_DIR)
from contig_spec import ContigSpec, FixedSegment  # noqa: E402
from haddock_run import (  # noqa: E402
    parse_active_residues,
    parse_contact_pairs,
    parse_pair_distance,
)


# Receptor/effector chains inside the HADDOCK workspace (post-relabel).
HADDOCK_RECEPTOR_CHAIN = "A"
HADDOCK_EFFECTOR_CHAIN = "B"

# AIR satisfaction fraction — 50% per A146, hardcoded constant.
AIR_SATISFACTION_FRACTION = 0.5


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--receptor", required=True,
                        help="Receptor monomer PDB file")
    parser.add_argument("--effector", required=True,
                        help="Effector monomer PDB file")
    parser.add_argument("--receptor-chain", default="A",
                        help="Receptor chain ID in input PDB (default: A)")
    parser.add_argument("--effector-chain", default="B",
                        help="Effector chain ID in input PDB (default: B)")
    # Contact-pair mode (case a) — primary workhorse per A136
    parser.add_argument("--contact-pairs", default="",
                        help='Pairwise CA-CA distance restraints, format '
                             '"A25-C42 A13-C94".  Hard pin.  Empty by default.')
    # Active-residues mode (case b)
    parser.add_argument("--receptor-active-residues", default="",
                        help='Comma-separated receptor residue numbers / ranges, '
                             'e.g. "25,35,40-44".  Empty by default.')
    parser.add_argument("--effector-active-residues", default="",
                        help='Comma-separated effector residue numbers / ranges. '
                             'Empty = restrain to entire effector chain.')
    # Global pair distance
    parser.add_argument("--pair-distance", default="2,2,4",
                        help='Global distance triple "target,lo_dev,hi_dev" '
                             'for all contact pairs (default 2,2,4).')
    # Contig string — for receptor design region derivation (clash bookkeeping
    # only; AIRs still come from --receptor-active-residues).  Per Session 7
    # post-commit-3 amendment in notes/design_audit.md.
    parser.add_argument("--contigs", default="",
                        help='RFDiffusion contig string (e.g. '
                             '"A1-32/10-30/A50-68/10-10 B").  Used at HADDOCK '
                             'time to derive the receptor design region for '
                             'clash bookkeeping — receptor residues NOT in '
                             'fixed contig segments are considered design '
                             'region.  Does NOT drive AIR generation.')
    return parser.parse_args()


# ── PDB relabelling ─────────────────────────────────────────────────


def write_pdb_with_chain(in_path: Path, out_path: Path,
                         src_chain: str, dst_chain: str) -> int:
    """Copy ``in_path`` to ``out_path``, restricting to ATOM records on
    ``src_chain`` and rewriting the chain column to ``dst_chain``.
    Returns the number of ATOM lines written.

    Other record types (HEADER, TITLE, TER, END) are preserved but
    re-written with the dst chain letter where applicable (TER lines).
    """
    src_chain = src_chain.upper()
    dst_chain = dst_chain.upper()
    n_written = 0
    with open(in_path) as fin, open(out_path, "w") as fout:
        for line in fin:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                if line[21].upper() != src_chain:
                    continue
                # Rewrite chain column (index 21).
                line = line[:21] + dst_chain + line[22:]
                fout.write(line)
                n_written += 1
            elif line.startswith("TER"):
                if len(line) > 21 and line[21].upper() == src_chain:
                    line = line[:21] + dst_chain + line[22:]
                fout.write(line)
            elif line.startswith("END"):
                fout.write(line)
                break
            else:
                # HEADER / TITLE / SEQRES etc. — pass through.
                fout.write(line)
    return n_written


# ── Contig-derived design region (for clash bookkeeping) ───────────


def read_pdb_residues(pdb_path: Path, chain: str) -> set:
    """Return the set of receptor residue numbers present in the chain."""
    out = set()
    chain = chain.upper()
    with open(pdb_path) as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            if line[21].upper() != chain:
                continue
            try:
                out.add(int(line[22:26].strip()))
            except ValueError:
                continue
    return out


def receptor_design_region_from_contig(
    contig_string: str,
    receptor_pdb: Path,
    receptor_chain_in_pdb: str,
) -> list:
    """Compute the set of receptor PDB residues that fall in de novo
    (non-fixed) regions of the contig.  Returns a sorted list.

    Algorithm: parse contig → collect FixedSegment residue ranges on
    the receptor chain → subtract from the set of residues actually
    present in the receptor PDB.

    Notes:
    - The chain letter inside the contig MUST match the receptor's
      chain in the user's input PDB (e.g. "A" if the receptor PDB
      uses chain A).  HADDOCK_PREPARE relabels to A/B for HADDOCK
      input, but the contig still references the USER's chain letters.
    - Returns [] when the contig is empty (caller decides whether
      that's an error — the validator should have rejected empty in
      Branch A already).
    - Returns the entire receptor PDB residue list if no FixedSegment
      references the receptor chain (degenerate "everything is de novo"
      case).
    """
    contig_string = (contig_string or "").strip()
    if not contig_string:
        return []
    try:
        spec = ContigSpec.from_string(contig_string)
    except ValueError as e:
        print(f"WARNING: cannot parse contig {contig_string!r}: {e}; "
              f"design region empty.", file=sys.stderr)
        return []

    chain_letter = receptor_chain_in_pdb.upper()
    fixed_residues: set = set()
    chain_present = False
    for ch in spec.chains:
        if ch.chain_id.upper() == chain_letter:
            chain_present = True
            for seg in ch.segments:
                if isinstance(seg, FixedSegment):
                    fixed_residues.update(range(seg.start, seg.end + 1))
    if not chain_present:
        print(f"WARNING: receptor chain {chain_letter!r} not in contig "
              f"{contig_string!r}; design region empty.", file=sys.stderr)
        return []

    pdb_residues = read_pdb_residues(receptor_pdb, chain_letter)
    return sorted(pdb_residues - fixed_residues)


# ── AIR generation ──────────────────────────────────────────────────


def write_air_restraints(
    path: Path,
    receptor_active: list,
    receptor_chain: str,
    effector_chain: str,
    effector_active: list,
) -> None:
    """Write the ambiguous interaction restraints file.

    Four cases (per Session 7 A133, plus the post-commit-3 amendment
    that adds the effector-only case):

    - **Both sides specified**: two-sided AIR, distance 3.0 +- 3.0 +- 5.0 A.
      Each receptor active residue must contact (OR-of-every) effector
      active residue.
    - **Receptor only**: one-sided AIR, distance 5.0 +- 5.0 +- 5.0 A.
      Each receptor active residue must contact any residue on the
      effector chain.
    - **Effector only**: one-sided AIR, distance 5.0 +- 5.0 +- 5.0 A.
      Each effector active residue must contact any residue on the
      receptor chain.  The "I want this effector face involved but
      I don't know which receptor residues" mode.
    - **Neither**: empty file.
    """
    if not receptor_active and not effector_active:
        path.write_text("")
        return

    lines = ["! Ambiguous Interaction Restraints (HADDOCK AIRs)\n"]

    if receptor_active and effector_active:
        # Two-sided: each receptor active → OR of every effector active.
        lines.append(
            f"! Receptor active residues ({len(receptor_active)}): "
            f"{','.join(str(r) for r in receptor_active)}\n"
        )
        lines.append(
            f"! Effector active residues ({len(effector_active)}): "
            f"{','.join(str(r) for r in effector_active)}\n"
        )
        for resnum in receptor_active:
            lines.append(
                f"assign (resid {resnum} and segid {receptor_chain})\n"
            )
            lines.append("       (\n")
            for i, eff_res in enumerate(effector_active):
                connector = "or" if i < len(effector_active) - 1 else "  "
                lines.append(
                    f"        (resid {eff_res} and segid {effector_chain}) "
                    f"{connector}\n"
                )
            lines.append("       ) 3.0 3.0 5.0\n")
    elif receptor_active:
        # Receptor-only: each receptor active → entire effector chain.
        lines.append(
            f"! Receptor active residues ({len(receptor_active)}): "
            f"{','.join(str(r) for r in receptor_active)}\n"
        )
        lines.append(f"! Effector side: entire chain {effector_chain}\n")
        for resnum in receptor_active:
            lines.append(
                f"assign (resid {resnum} and segid {receptor_chain})\n"
                f"       ((segid {effector_chain})) 5.0 5.0 5.0\n"
            )
    else:
        # Effector-only: each effector active → entire receptor chain.
        lines.append(f"! Receptor side: entire chain {receptor_chain}\n")
        lines.append(
            f"! Effector active residues ({len(effector_active)}): "
            f"{','.join(str(r) for r in effector_active)}\n"
        )
        for resnum in effector_active:
            lines.append(
                f"assign (resid {resnum} and segid {effector_chain})\n"
                f"       ((segid {receptor_chain})) 5.0 5.0 5.0\n"
            )
    path.write_text("".join(lines))


# ── Unambig (contact-pair) restraints ───────────────────────────────


def write_unambig_restraints(
    path: Path,
    contact_pairs: list,
    pair_distance: tuple,
    receptor_chain: str,
    effector_chain: str,
) -> None:
    """Write unambiguous CA-CA distance restraints, one per pair.

    Format per HADDOCK CNS convention:
        assign (name CA and resid R and segid A)
               (name CA and resid E and segid B) target lo_dev hi_dev

    No-op when contact_pairs is empty: nothing to restrain.
    """
    if not contact_pairs:
        path.write_text("")
        return
    target, lo_dev, hi_dev = pair_distance
    lines = [
        "! Unambiguous CA-CA distance restraints (contact-pair mode)\n",
        f"! Target {target} A, lower dev {lo_dev} A, upper dev {hi_dev} A "
        f"=> distance window [{max(0.0, target - lo_dev):.2f}, "
        f"{target + hi_dev:.2f}] A\n",
    ]
    for _src_rec_chain, rec_resnum, _src_eff_chain, eff_resnum in contact_pairs:
        # NOTE: the chain letters in the pair string refer to the USER'S
        # input PDB chains.  The .tbl file uses HADDOCK_RECEPTOR_CHAIN /
        # HADDOCK_EFFECTOR_CHAIN (A/B) because we relabel to A/B above.
        lines.append(
            f"assign (name CA and resid {rec_resnum} and segid {receptor_chain})\n"
            f"       (name CA and resid {eff_resnum} and segid {effector_chain}) "
            f"{target:.2f} {lo_dev:.2f} {hi_dev:.2f}\n"
        )
    path.write_text("".join(lines))


# ── nrest tolerance hint ────────────────────────────────────────────


def air_nrest(receptor_active: list, effector_active: list) -> int:
    """How many of the N AIRs HADDOCK should require satisfied.

    50% per A146; HADDOCK's [airs] block reads this as "at least this
    many must satisfy."  The total AIR count is whichever side has
    active residues (one assign-block per receptor active in two-sided
    or receptor-only mode; one per effector active in effector-only).
    """
    total = len(receptor_active) if receptor_active else len(effector_active)
    if total <= 0:
        return 0
    return max(1, math.ceil(AIR_SATISFACTION_FRACTION * total))


# ── Main ────────────────────────────────────────────────────────────


def main():
    args = parse_args()

    # ── Parse restraint params ───────────────────────────────────────
    contact_pairs = parse_contact_pairs(args.contact_pairs)
    receptor_active = list(parse_active_residues(args.receptor_active_residues))
    effector_active = list(parse_active_residues(args.effector_active_residues))
    pair_distance = parse_pair_distance(args.pair_distance)

    # Validation — at least one restraint mode must be active.
    # (The pipeline validator should have caught this already; defensive.)
    if not contact_pairs and not receptor_active and not effector_active:
        sys.exit(
            "ERROR: HADDOCK requires at least one of --contact-pairs, "
            "--receptor-active-residues, or --effector-active-residues to be "
            "non-empty.  Blind docking is not supported by this pipeline; see "
            "notes/design_audit.md A135 for the workflow when you don't have "
            "a target interface in mind."
        )

    # ── Derive receptor design region from contig (clash bookkeeping) ────
    # Used only by haddock_cluster_metrics.py to classify clashes as
    # "in design region" (tolerated; will be redesigned by RFDiffusion)
    # vs "outside design region" (real geometric problems).  Does NOT
    # affect AIR generation.  Per Session 7 post-commit-3 amendment.
    contig_design_region = receptor_design_region_from_contig(
        args.contigs, Path(args.receptor), args.receptor_chain,
    )

    # ── Relabel input PDBs to chains A / B ───────────────────────────
    n_rec_atoms = write_pdb_with_chain(
        Path(args.receptor),
        Path("receptor_haddock.pdb"),
        src_chain=args.receptor_chain,
        dst_chain=HADDOCK_RECEPTOR_CHAIN,
    )
    n_eff_atoms = write_pdb_with_chain(
        Path(args.effector),
        Path("effector_haddock.pdb"),
        src_chain=args.effector_chain,
        dst_chain=HADDOCK_EFFECTOR_CHAIN,
    )
    if n_rec_atoms == 0:
        sys.exit(
            f"ERROR: no ATOM records found on receptor chain "
            f"'{args.receptor_chain}' in {args.receptor}"
        )
    if n_eff_atoms == 0:
        sys.exit(
            f"ERROR: no ATOM records found on effector chain "
            f"'{args.effector_chain}' in {args.effector}"
        )
    print(f"Relabelled receptor chain {args.receptor_chain} -> "
          f"{HADDOCK_RECEPTOR_CHAIN} ({n_rec_atoms} atoms) "
          f"-> receptor_haddock.pdb")
    print(f"Relabelled effector chain {args.effector_chain} -> "
          f"{HADDOCK_EFFECTOR_CHAIN} ({n_eff_atoms} atoms) "
          f"-> effector_haddock.pdb")

    # ── Write restraint files ────────────────────────────────────────
    write_air_restraints(
        Path("ambig_restraints.tbl"),
        receptor_active,
        HADDOCK_RECEPTOR_CHAIN,
        HADDOCK_EFFECTOR_CHAIN,
        effector_active,
    )
    write_unambig_restraints(
        Path("unambig_restraints.tbl"),
        contact_pairs,
        pair_distance,
        HADDOCK_RECEPTOR_CHAIN,
        HADDOCK_EFFECTOR_CHAIN,
    )

    # ── Summary JSON for haddock.nf to consume ───────────────────────
    nrest = air_nrest(receptor_active, effector_active)
    summary = {
        "receptor_chain": HADDOCK_RECEPTOR_CHAIN,
        "effector_chain": HADDOCK_EFFECTOR_CHAIN,
        "src_receptor_chain": args.receptor_chain,
        "src_effector_chain": args.effector_chain,
        "n_receptor_atoms": n_rec_atoms,
        "n_effector_atoms": n_eff_atoms,
        "contact_pairs": [
            {"rec_chain": rc, "rec_resnum": rn,
             "eff_chain": ec, "eff_resnum": en}
            for rc, rn, ec, en in contact_pairs
        ],
        "n_contact_pairs": len(contact_pairs),
        "receptor_active_residues": receptor_active,
        "effector_active_residues": effector_active,
        "pair_distance": list(pair_distance),
        "air_nrest": nrest,
        "air_satisfaction_fraction": AIR_SATISFACTION_FRACTION,
        "ambig_restraints_present": bool(receptor_active) or bool(effector_active),
        "unambig_restraints_present": bool(contact_pairs),
        # Contig-derived receptor design region for clash bookkeeping.
        # The receptor PDB chain letter under USER's numbering — not the
        # HADDOCK A/B relabel — because contig + native PDB share a frame.
        "contig": args.contigs,
        "receptor_chain_in_pdb": args.receptor_chain,
        "contig_design_region": contig_design_region,
    }
    Path("restraints_summary.json").write_text(json.dumps(summary, indent=2))

    # ── Diagnostic output ────────────────────────────────────────────
    if contact_pairs:
        print(f"Wrote {len(contact_pairs)} unambig pair restraint(s) "
              f"-> unambig_restraints.tbl  "
              f"(target {pair_distance[0]}, "
              f"window [{max(0.0, pair_distance[0] - pair_distance[1]):.2f}, "
              f"{pair_distance[0] + pair_distance[2]:.2f}] A)")
    else:
        print("No contact pairs -> unambig_restraints.tbl is empty")
    if receptor_active and effector_active:
        print(f"Wrote AIRs: two-sided ({len(receptor_active)} receptor "
              f"active x {len(effector_active)} effector active) -> "
              f"ambig_restraints.tbl  "
              f"(nrest = {nrest}, "
              f"{AIR_SATISFACTION_FRACTION*100:.0f}% satisfaction required)")
    elif receptor_active:
        print(f"Wrote AIRs: receptor-only ({len(receptor_active)} receptor "
              f"active -> entire effector chain) -> ambig_restraints.tbl  "
              f"(nrest = {nrest})")
    elif effector_active:
        print(f"Wrote AIRs: effector-only ({len(effector_active)} effector "
              f"active -> entire receptor chain) -> ambig_restraints.tbl  "
              f"(nrest = {nrest})")
    else:
        print("No active residues -> ambig_restraints.tbl is empty")
    if contig_design_region:
        ranges = []
        run_start = prev = contig_design_region[0]
        for r in contig_design_region[1:]:
            if r != prev + 1:
                ranges.append(f"{run_start}-{prev}" if run_start != prev else str(run_start))
                run_start = r
            prev = r
        ranges.append(f"{run_start}-{prev}" if run_start != prev else str(run_start))
        print(f"Receptor design region from contig ({len(contig_design_region)} "
              f"residues): {','.join(ranges)}")
    elif args.contigs:
        print("Contig provided but receptor design region is empty "
              "(check that the contig's receptor chain matches "
              f"--receptor-chain '{args.receptor_chain}')")
    else:
        print("No contig provided -> contig_design_region is empty "
              "-> clash bookkeeping falls back to pair + receptor-active set")


if __name__ == "__main__":
    main()
