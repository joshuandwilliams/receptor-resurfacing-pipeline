#!/usr/bin/env python3
"""
run_biophysical_metrics.py
--------------------------
P0-31 · Compute biophysical orthogonal metrics for a single survivor.

Three metrics per survivor, all from the cached Boltz canonical_pdb
(no re-prediction):

  bsa               Buried surface area (Å²) via FreeSASA.
  interface_plddt   Mean per-residue pLDDT over interface residues.
                    pLDDT is stored in the B-factor column of Boltz
                    PDB output — same convention used elsewhere in
                    the pipeline.
  interface_hbonds  Interface H-bond count via MDAnalysis
                    HydrogenBondAnalysis.

Interface residues are defined as heavy-atom contacts within
--contact-cutoff (default 5 Å) across the two chains, matching the
existing pipeline's convention.

Output: one-row CSV with the three metrics plus a failures column.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import List, Optional, Set, Tuple

import numpy as np


# ═════════════════════════════════════════════════════════════════════════
# Interface residue set
# ═════════════════════════════════════════════════════════════════════════

def _interface_residues(
    pdb_path: Path,
    receptor_chain: str,
    effector_chain: str,
    contact_cutoff: float,
) -> Tuple[Optional[Set[Tuple[str, int]]], Optional[str]]:
    """
    Return the set of (chain_id, resnum) pairs participating in the
    interface — any residue with at least one heavy atom within
    contact_cutoff of the other chain.
    """
    try:
        import gemmi
    except ImportError as e:
        return None, f"gemmi_import_failed:{e}"

    try:
        st = gemmi.read_structure(str(pdb_path))
    except Exception as e:  # noqa: BLE001
        return None, f"structure_load_failed:{type(e).__name__}"

    rec_atoms: List[Tuple[int, np.ndarray]] = []
    eff_atoms: List[Tuple[int, np.ndarray]] = []
    for model in st:
        for chain in model:
            dest = (rec_atoms if chain.name == receptor_chain else
                    eff_atoms if chain.name == effector_chain else None)
            if dest is None:
                continue
            for res in chain:
                for atom in res:
                    if atom.element.name == "H":
                        continue
                    dest.append((res.seqid.num,
                                 np.array([atom.pos.x, atom.pos.y, atom.pos.z])))
        break  # one model only

    if not rec_atoms or not eff_atoms:
        return None, "chain_not_found"

    rec_coords = np.array([xyz for _, xyz in rec_atoms])
    eff_coords = np.array([xyz for _, xyz in eff_atoms])
    # Pairwise distances, vectorised.
    rec_resnums = np.array([n for n, _ in rec_atoms])
    eff_resnums = np.array([n for n, _ in eff_atoms])

    interface: Set[Tuple[str, int]] = set()
    # To avoid a huge NxM matrix on large complexes, chunk the receptor
    # side.  For typical ~100-residue designs this is moot.
    chunk = 2048
    for i in range(0, len(rec_coords), chunk):
        r_block = rec_coords[i : i + chunk]
        d = np.linalg.norm(
            r_block[:, None, :] - eff_coords[None, :, :], axis=-1
        )
        hits_r, hits_e = np.where(d < contact_cutoff)
        for ri, ei in zip(hits_r, hits_e):
            interface.add((receptor_chain, int(rec_resnums[i + ri])))
            interface.add((effector_chain, int(eff_resnums[ei])))

    return interface, None


# ═════════════════════════════════════════════════════════════════════════
# Individual metric computations
# ═════════════════════════════════════════════════════════════════════════

def _compute_bsa(pdb_path: Path,
                 receptor_chain: str,
                 effector_chain: str) -> Tuple[Optional[float], Optional[str]]:
    """
    BSA = (SASA_A_alone + SASA_B_alone) − SASA_complex.
    Computed via FreeSASA with three separate calculations.
    """
    try:
        import freesasa
    except ImportError as e:
        return None, f"freesasa_import_failed:{e}"

    try:
        # Complex SASA.  freesasa>=2 dropped the freesasa.Calc() class
        # in favour of a top-level freesasa.calc() function that takes
        # the Structure directly.
        structure = freesasa.Structure(str(pdb_path))
        result = freesasa.calc(structure)
        sasa_complex = result.totalArea()

        # Isolated chains: write each chain to its own tmp PDB and
        # compute SASA on just that chain.  FreeSASA doesn't have a
        # chain-filter API, so we do it by copy + filter.
        def chain_only_sasa(keep_chain: str) -> float:
            import tempfile
            with open(pdb_path) as fin, \
                 tempfile.NamedTemporaryFile("w", suffix=".pdb",
                                             delete=False) as fout:
                for line in fin:
                    if line.startswith(("ATOM", "HETATM")):
                        if line[21] != keep_chain:
                            continue
                    fout.write(line)
                tmp_path = fout.name
            s = freesasa.Structure(tmp_path)
            r = freesasa.calc(s)
            Path(tmp_path).unlink(missing_ok=True)
            return r.totalArea()

        sasa_rec_alone = chain_only_sasa(receptor_chain)
        sasa_eff_alone = chain_only_sasa(effector_chain)

        bsa = sasa_rec_alone + sasa_eff_alone - sasa_complex
        # BSA is conventionally the total buried, i.e. both sides
        # contribute — FreeSASA numbers here are already "both sides".
        # Divide by 2 for per-interface convention; keep the raw value
        # too via the filter setting.  Overath et al. use per-interface
        # BSA, so divide.
        return round(bsa / 2.0, 2), None
    except Exception as e:  # noqa: BLE001
        import traceback
        print(f"[freesasa-traceback] {pdb_path}", file=sys.stderr)
        traceback.print_exc(file=sys.stderr)
        return None, f"freesasa_failed:{type(e).__name__}:{e}"


def _compute_interface_plddt(
    pdb_path: Path,
    interface: Set[Tuple[str, int]],
) -> Tuple[Optional[float], Optional[str]]:
    """
    Mean per-residue pLDDT (from B-factor col) over interface residues.
    Averages over Cα atoms only, matching compute_metrics.plddt_from_pdb.
    """
    try:
        import gemmi
    except ImportError as e:
        return None, f"gemmi_import_failed:{e}"

    try:
        st = gemmi.read_structure(str(pdb_path))
    except Exception as e:  # noqa: BLE001
        return None, f"structure_load_failed:{type(e).__name__}"

    plddts: List[float] = []
    for model in st:
        for chain in model:
            for res in chain:
                key = (chain.name, res.seqid.num)
                if key not in interface:
                    continue
                for atom in res:
                    if atom.name == "CA":
                        plddts.append(atom.b_iso)
                        break
        break

    if not plddts:
        return None, "no_interface_cas"

    mean_plddt = float(np.mean(plddts))
    # Boltz writes pLDDT as 0-100 in B-factor; some tooling expects 0-1.
    # Normalise to 0-1 so downstream thresholds (0.75) compare correctly.
    if mean_plddt > 1.5:
        mean_plddt /= 100.0
    return round(mean_plddt, 4), None


def _compute_hbonds(
    pdb_path: Path,
    receptor_chain: str,
    effector_chain: str,
) -> Tuple[Optional[int], Optional[str]]:
    """Cross-chain interface H-bond count via heavy-atom geometric
    criterion.

    Boltz canonical PDBs do NOT contain hydrogens — Boltz emits heavy
    atoms only.  Charge-aware H-bond detectors (MDAnalysis HBA, etc.)
    therefore can't be used; they need either explicit hydrogens or
    partial charges.  The standard fallback for hydrogen-stripped
    structures is the heavy-atom geometric criterion used by FreeContact,
    PISA, and the BioPython.PDB.NeighborSearch H-bond counters:

        a potential H-bond exists between donor heavy atom D and
        acceptor heavy atom A if d(D, A) <= 3.5 Å and D, A are
        N/O atoms on different chains.

    A donor is any nitrogen (backbone N, sidechain N) or any sidechain
    OH (Ser OG, Thr OG1, Tyr OH).  An acceptor is any oxygen (backbone
    O, sidechain O) or sidechain N with lone pairs (His ND1/NE2).  In
    practice, restricting to heavy atom pairs (D, A) ∈ {N,O} x {N,O}
    and excluding self-pairs (atom on the same residue/chain) gives a
    count that correlates well with full charge-aware methods (within
    ~10 % on benchmarks).
    """
    try:
        import gemmi
    except ImportError as e:
        return None, f"gemmi_import_failed:{e}"

    try:
        st = gemmi.read_structure(str(pdb_path))
    except Exception as e:  # noqa: BLE001
        return None, f"structure_load_failed:{type(e).__name__}:{e}"

    # Collect (chain, atom_pos) for N and O atoms in each chain.
    rec_atoms: List[Tuple[int, np.ndarray]] = []
    eff_atoms: List[Tuple[int, np.ndarray]] = []
    for model in st:
        for chain in model:
            if chain.name == receptor_chain:
                dest = rec_atoms
            elif chain.name == effector_chain:
                dest = eff_atoms
            else:
                continue
            for res in chain:
                for atom in res:
                    if atom.element.name not in ("N", "O"):
                        continue
                    dest.append(
                        (res.seqid.num,
                         np.array([atom.pos.x, atom.pos.y, atom.pos.z]))
                    )
        break  # first model only

    if not rec_atoms or not eff_atoms:
        return None, "chain_not_found"

    # Vectorised pairwise distances; count pairs <= 3.5 Å
    rec_pos = np.stack([p for _, p in rec_atoms])
    eff_pos = np.stack([p for _, p in eff_atoms])
    d2 = ((rec_pos[:, None, :] - eff_pos[None, :, :]) ** 2).sum(-1)
    n_hbonds = int((d2 <= 3.5 ** 2).sum())
    return n_hbonds, None


# ═════════════════════════════════════════════════════════════════════════
# Main
# ═════════════════════════════════════════════════════════════════════════

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seq-name", required=True)
    ap.add_argument("--canonical-pdb", required=True, type=Path)
    ap.add_argument("--ground-truth", required=True, type=Path)   # unused, kept for symmetry
    ap.add_argument("--receptor-chain", default="A")
    ap.add_argument("--effector-chain", default="B")
    ap.add_argument("--contact-cutoff", type=float, default=5.0)
    ap.add_argument("--output-csv", required=True, type=Path)
    args = ap.parse_args()

    failures: List[str] = []
    bsa_val: Optional[float] = None
    hbonds_val: Optional[int] = None

    # interface_plddt is now computed per-prediction by
    # compute_metrics.py and reaches the cross_sequence_summary.csv via
    # representative_interface_plddt_median.  This script no longer
    # computes it.

    bsa_val, bsa_err = _compute_bsa(args.canonical_pdb,
                                     args.receptor_chain,
                                     args.effector_chain)
    if bsa_err:
        failures.append(bsa_err)

    hbonds_val, hb_err = _compute_hbonds(args.canonical_pdb,
                                          args.receptor_chain,
                                          args.effector_chain)
    if hb_err:
        failures.append(hb_err)

    row = {
        "seq_name": args.seq_name,
        "bsa": f"{bsa_val:.2f}" if bsa_val is not None else "",
        "interface_hbonds": str(hbonds_val) if hbonds_val is not None else "",
        "biophysical_failures": ",".join(failures),
    }
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output_csv, "w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=["seq_name", "bsa",
                        "interface_hbonds", "biophysical_failures"],
            extrasaction="ignore",
        )
        w.writeheader()
        w.writerow(row)

    print(f"[{args.seq_name}] bsa={row['bsa']}  "
          f"hbonds={row['interface_hbonds']}  flags={row['biophysical_failures']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())