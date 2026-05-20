#!/usr/bin/env python3
"""
pose_solver_prepare.py
----------------------
Normalise the two input PDBs for the pose solver: rewrite the
receptor chain to ``--receptor-chain`` and the effector chain to
``--effector-chain`` so the rest of the pipeline can rely on the
user-declared chain IDs being present in every downstream PDB.

``build_contigs.py`` falls back to the first two alphabetical chains
when the requested chain IDs are absent, but the contig string it
emits still references the user's claimed chain IDs.  Without this
normalisation step, RFDiffusion would receive a PDB with one set of
chain IDs and a contig string naming a different set.

Outputs:
    receptor_pose.pdb   - input receptor with chain ID forced to
                          ``--receptor-chain``
    effector_pose.pdb   - input effector with chain ID forced to
                          ``--effector-chain``
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path


def relabel_chain(in_pdb: Path, out_pdb: Path, new_chain: str) -> int:
    """Copy in_pdb to out_pdb with column 22 (chain ID) overwritten to
    new_chain on every ATOM / HETATM / TER record.  Other record types
    pass through unchanged.  Returns the number of atom records written.
    """
    n_atoms = 0
    with open(in_pdb) as f_in, open(out_pdb, "w") as f_out:
        for line in f_in:
            rec = line[:6]
            if rec in ("ATOM  ", "HETATM", "TER   ", "ANISOU"):
                if len(line) < 22:
                    f_out.write(line)
                    continue
                f_out.write(line[:21] + new_chain + line[22:])
                if rec in ("ATOM  ", "HETATM"):
                    n_atoms += 1
            else:
                f_out.write(line)
    return n_atoms


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--receptor", required=True, type=Path)
    ap.add_argument("--effector", required=True, type=Path)
    ap.add_argument("--receptor-chain", default="A",
                    help="Chain ID to assign to receptor (default A)")
    ap.add_argument("--effector-chain", default="B",
                    help="Chain ID to assign to effector (default B)")
    ap.add_argument("--receptor-out", default="receptor_pose.pdb", type=Path)
    ap.add_argument("--effector-out", default="effector_pose.pdb", type=Path)
    return ap.parse_args()


def main() -> int:
    args = parse_args()
    rc = args.receptor_chain.upper()
    ec = args.effector_chain.upper()
    if rc == ec:
        sys.exit(f"ERROR: receptor-chain and effector-chain are both "
                 f"'{rc}'.  Pass distinct chain IDs.")
    if len(rc) != 1 or len(ec) != 1:
        sys.exit(f"ERROR: chain IDs must be single characters, got "
                 f"receptor={rc!r}, effector={ec!r}.")

    n_rec = relabel_chain(args.receptor, args.receptor_out, rc)
    n_eff = relabel_chain(args.effector, args.effector_out, ec)
    print(f"Wrote {args.receptor_out} ({n_rec} atoms, chain {rc})")
    print(f"Wrote {args.effector_out} ({n_eff} atoms, chain {ec})")
    if n_rec == 0:
        sys.exit(f"ERROR: receptor PDB {args.receptor} had no atom records.")
    if n_eff == 0:
        sys.exit(f"ERROR: effector PDB {args.effector} had no atom records.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
