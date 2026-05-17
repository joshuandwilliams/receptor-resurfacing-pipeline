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


# ── AIR generation ──────────────────────────────────────────────────


def write_air_restraints(
    path: Path,
    receptor_active: list,
    receptor_chain: str,
    effector_chain: str,
    effector_active: list,
) -> None:
    """Write the ambiguous interaction restraints file.

    AND-over-receptor / OR-over-effector with a 2-sided 3.0 +- 3.0 +- 5.0 A
    distance restraint when both sides are specified, or a one-sided
    "receptor residue to entire effector chain" 5.0 +- 5.0 +- 5.0 A
    restraint when only the receptor side is.

    No-op when receptor_active is empty: nothing to restrain.
    """
    if not receptor_active:
        path.write_text("")
        return

    lines = ["! Ambiguous Interaction Restraints (HADDOCK AIRs)\n",
             f"! Receptor active residues ({len(receptor_active)}): "
             f"{','.join(str(r) for r in receptor_active)}\n"]
    if effector_active:
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
    else:
        lines.append(f"! Effector active residues: entire chain {effector_chain}\n")
        for resnum in receptor_active:
            lines.append(
                f"assign (resid {resnum} and segid {receptor_chain})\n"
                f"       ((segid {effector_chain})) 5.0 5.0 5.0\n"
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


def air_nrest(n_active: int) -> int:
    """How many of the N AIRs HADDOCK should require satisfied.

    50% per A146; HADDOCK's [airs] block reads this as "at least this
    many must satisfy."  Returned for callers to plumb into docking.cfg
    (the HADDOCK3 noecv / nrest knob).
    """
    if n_active <= 0:
        return 0
    return max(1, math.ceil(AIR_SATISFACTION_FRACTION * n_active))


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
    if not contact_pairs and not receptor_active:
        sys.exit(
            "ERROR: HADDOCK requires at least one of --contact-pairs or "
            "--receptor-active-residues to be non-empty.  Blind docking is "
            "not supported by this pipeline; see notes/design_audit.md A135 "
            "for the workflow when you don't have a target interface in mind."
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
    nrest = air_nrest(len(receptor_active))
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
        "ambig_restraints_present": bool(receptor_active),
        "unambig_restraints_present": bool(contact_pairs),
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
    if receptor_active:
        print(f"Wrote AIRs for {len(receptor_active)} receptor active "
              f"residue(s) -> ambig_restraints.tbl "
              f"(nrest = {nrest}, "
              f"{AIR_SATISFACTION_FRACTION*100:.0f}% satisfaction required)")
        if effector_active:
            print(f"  Effector active residues: {len(effector_active)}")
        else:
            print(f"  Effector active residues: entire chain "
                  f"{HADDOCK_EFFECTOR_CHAIN}")
    else:
        print("No receptor active residues -> ambig_restraints.tbl is empty")


if __name__ == "__main__":
    main()
