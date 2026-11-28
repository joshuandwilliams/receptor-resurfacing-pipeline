"""
structure_metrics.py
--------------------
Geometric / structural primitives shared across the pipeline.

Holds the centre-of-mass and clash-counting code that used to live
inline in ``rfdiffusion_filter.py``, which now imports it from here.

Contract
========
- All functions are pure.  No globals, and no I/O beyond reading the
  path passed in.
- Functions that take a PDB read the first model only and consider ATOM
  records exclusively, so HETATMs are ignored.
- The coordinate readers key residues on resseq alone, so two residues
  sharing a resseq with different insertion codes collapse into one
  pick.  This matches rfdiffusion_filter's pre-existing behaviour.
  ``get_chain_sequence`` is the exception and keys on (resnum, icode).

Public API
==========
- ``centroid(coords)`` — mean along axis 0.  Raises ValueError on empty
  input or on any non-finite coordinate, rather than returning NaN.
- ``centroid_distance(c1, c2)`` — Euclidean distance between two centroids.
- ``read_chain_ca_coords(pdb_path, chain_id)`` — Nx3 numpy array of CA
  coordinates for one chain, in residue order.
- ``chain_centre_of_mass(pdb_path, chain_id)`` — composes the above two.
  Returns ``None`` if the chain has no CAs.
- ``chain_com_distance(pdb_path, chain_a, chain_b)`` — composes the above.
- ``read_chain_heavy_atoms(pdb_path, chain_id)`` — list of
  ``(resnum, atom_name, x, y, z)`` for all non-hydrogen ATOM records.
- ``clash_count(pdb_path, chain_a, chain_b, design_region_residues_a,
  cutoff=2.0)`` — ``(in_design, outside_design)`` heavy-atom pair counts
  under the cutoff.  A pair counts as "in design" when its chain_a
  residue number is in ``design_region_residues_a``.
- ``get_chain_sequence(pdb_path, chain_id)`` — one-letter sequence.
- ``find_contact_residues_heavy(pdb_path, receptor_chain,
  effector_chain, cutoff=5.0)`` — 0-based positional indices of
  receptor residues contacting the effector.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, List, Optional, Set, Tuple

import numpy as np

# ── Numpy-level primitives ──────────────────────────────────────────


def centroid(coords: np.ndarray) -> np.ndarray:
    """Return ``coords.mean(axis=0)`` as a 1D length-3 array.

    Empty coords, or coords containing NaN or Inf, raise ValueError.
    Otherwise callers would silently propagate a NaN centroid.
    """
    if coords.size == 0:
        raise ValueError("centroid: empty coordinate array")
    if not np.isfinite(coords).all():
        raise ValueError("centroid: coordinates contain NaN/Inf")
    return coords.mean(axis=0)


def centroid_distance(c1: np.ndarray, c2: np.ndarray) -> float:
    """Euclidean distance between two centroids."""
    return float(np.linalg.norm(c1 - c2))


# ── PDB readers ─────────────────────────────────────────────────────


def _iter_atom_lines(pdb_path: Path) -> Iterable[str]:
    """Yield ATOM lines from the first model in the PDB.

    Stops at ENDMDL, so multi-model files such as NMR ensembles or
    docking cluster outputs emit only the first model.
    """
    with open(pdb_path) as fh:
        for line in fh:
            if line.startswith("ENDMDL"):
                return
            if line.startswith("ATOM"):
                yield line


def read_chain_ca_coords(pdb_path: Path, chain_id: str) -> np.ndarray:
    """Return an Nx3 numpy array of CA coordinates for ``chain_id``,
    in PDB residue-record order.  Returns shape (0, 3) when the chain
    has no CAs or doesn't exist in the PDB.
    """
    coords: List[Tuple[float, float, float]] = []
    seen_residues: Set[int] = set()
    for line in _iter_atom_lines(Path(pdb_path)):
        if line[21] != chain_id:
            continue
        if line[12:16].strip() != "CA":
            continue
        try:
            resnum = int(line[22:26].strip())
        except ValueError:
            continue
        if resnum in seen_residues:
            # Defensive: ignore duplicate CAs from alternate locations.
            continue
        seen_residues.add(resnum)
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        coords.append((x, y, z))
    if not coords:
        return np.zeros((0, 3), dtype=float)
    return np.array(coords, dtype=float)


def chain_centre_of_mass(
    pdb_path: Path, chain_id: str
) -> Optional[np.ndarray]:
    """Centre of mass (CA-only) for one chain.  Returns ``None`` if the
    chain has no CAs in the PDB.

    This is the geometric centroid of the CAs, not a mass-weighted COM.
    For comparing two domains or measuring the distance between chains
    the difference does not matter, and CA-only is robust to side-chain
    differences across redesigns.
    """
    coords = read_chain_ca_coords(Path(pdb_path), chain_id)
    if coords.shape[0] == 0:
        return None
    return centroid(coords)


def chain_com_distance(
    pdb_path: Path, chain_a: str, chain_b: str
) -> Optional[float]:
    """Distance between the CA-centroids of two chains in one PDB.
    Returns ``None`` if either chain is empty.
    """
    a = chain_centre_of_mass(Path(pdb_path), chain_a)
    b = chain_centre_of_mass(Path(pdb_path), chain_b)
    if a is None or b is None:
        return None
    return centroid_distance(a, b)


# ── Heavy-atom + clash counting ─────────────────────────────────────


def read_chain_heavy_atoms(
    pdb_path: Path, chain_id: str
) -> List[Tuple[int, str, float, float, float]]:
    """List of ``(resnum, atom_name, x, y, z)`` for every non-hydrogen
    ATOM record on the chain.  Used by clash counting.
    """
    out: List[Tuple[int, str, float, float, float]] = []
    for line in _iter_atom_lines(Path(pdb_path)):
        if line[21] != chain_id:
            continue
        atom_name = line[12:16].strip()
        # Element column is 76:78; fall back to name-prefix when missing.
        element = line[76:78].strip()
        if element == "H" or (not element and atom_name.startswith("H")):
            continue
        try:
            resnum = int(line[22:26].strip())
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except ValueError:
            continue
        out.append((resnum, atom_name, x, y, z))
    return out


def clash_count(
    pdb_path: Path,
    chain_a: str,
    chain_b: str,
    design_region_residues_a: Set[int],
    cutoff: float = 2.0,
) -> Tuple[int, int]:
    """Count heavy-atom pairs within ``cutoff`` Å, one atom from chain_a
    and one from chain_b.  Splits the count by whether the chain_a
    residue is in the design region.

    Returns ``(clashes_in_design_region, clashes_outside_design_region)``.

    A pair counts as "in design region" iff the chain_a (receptor)
    residue number is in ``design_region_residues_a``.  Chain B atoms
    have no notion of design region at this stage.

    Uses a naive O(N*M) loop.  Typical 2-chain interface PDBs are one
    ~100-residue receptor against one ~100-residue effector, so ~10^4
    atom pairs, which is fast enough.  Swap in a KDTree if a future
    consumer has much larger inputs.
    """
    atoms_a = read_chain_heavy_atoms(Path(pdb_path), chain_a)
    atoms_b = read_chain_heavy_atoms(Path(pdb_path), chain_b)
    if not atoms_a or not atoms_b:
        return (0, 0)
    coords_b = np.array([(x, y, z) for _, _, x, y, z in atoms_b], dtype=float)

    in_design = 0
    outside = 0
    cutoff_sq = cutoff * cutoff
    for resnum_a, _atom_a, xa, ya, za in atoms_a:
        diffs = coords_b - np.array([xa, ya, za], dtype=float)
        sq_dists = np.einsum("ij,ij->i", diffs, diffs)
        n_clashes = int((sq_dists < cutoff_sq).sum())
        if n_clashes == 0:
            continue
        if resnum_a in design_region_residues_a:
            in_design += n_clashes
        else:
            outside += n_clashes
    return (in_design, outside)


# ── Sequence and interface contacts ──────────────────────────────────
# Both functions below were this repository's only reason to reach into the
# negative-steering engine, which is a separate repository installed
# separately.  They are reimplemented here on this module's own PDB reader
# rather than copied across, so nothing in this repository depends on engine
# internals.
#
# tests/characterization/test_structure_metrics_ported.py asserts the
# properties the pipeline depends on, on the tracked reference structures.  The
# engine's own copies were compared directly before they were deleted.  80 of
# 80 chain sequences and 120 of 120 contact sets agreed, the latter over three
# cutoffs on 40 structures.

def get_chain_sequence(pdb_path: Path, chain_id: str) -> str:
    """One-letter sequence of `chain_id`, in the order the file lists it.

    One residue per (resnum, icode), taken from the first ATOM record of each,
    so an altloc does not duplicate a position.  A non-standard three-letter
    name becomes 'X', keeping the sequence the same length as the residue list
    and so keeping every downstream positional index valid.

    Raises ValueError when the chain is absent, naming the chains that are
    present.  The alternative is a silent empty string that reads as an empty
    chain.
    """
    from contig_utils import THREE_TO_ONE

    seen: set = set()
    seq: List[str] = []
    found = False
    present: List[str] = []
    for line in _iter_atom_lines(Path(pdb_path)):
        ch = line[21]
        if ch not in present:
            present.append(ch)
        if ch != chain_id:
            continue
        found = True
        key = (line[22:26].strip(), line[26])      # resnum plus insertion code
        if key in seen:
            continue
        seen.add(key)
        seq.append(THREE_TO_ONE.get(line[17:20].strip().upper(), "X"))
    if not found:
        raise ValueError(
            f"Chain {chain_id} not found in {pdb_path}. "
            f"Available chains: {present}")
    return "".join(seq)


def find_contact_residues_heavy(
    pdb_path: Path,
    receptor_chain: str,
    effector_chain: str,
    cutoff: float = 5.0,
) -> List[int]:
    """0-based positional indices of receptor residues contacting the effector.

    A contact is any heavy atom of the residue within `cutoff` angstroms of any
    heavy atom of the effector chain.  read_chain_heavy_atoms drops hydrogens,
    so 5 Å is a heavy-atom-to-heavy-atom distance.

    The index is positional within the receptor chain, counted over distinct
    resnums in file order.  That matches the design-region and true-interface
    files, which are also positional.  Returning residue numbers instead would
    break on any structure that does not start at 1.  Note that this walk
    ignores insertion codes, whereas get_chain_sequence above keys on
    (resnum, icode), so the two agree only on structures without icodes.
    """
    import numpy as _np

    rec = read_chain_heavy_atoms(Path(pdb_path), receptor_chain)
    eff = read_chain_heavy_atoms(Path(pdb_path), effector_chain)
    if not rec or not eff:
        return []

    order: List[int] = []
    for resnum, _atom, _x, _y, _z in rec:
        if resnum not in order:
            order.append(resnum)
    index_of = {rn: i for i, rn in enumerate(order)}

    eff_xyz = _np.asarray([(x, y, z) for _r, _a, x, y, z in eff], dtype=float)
    hits: set = set()
    for resnum, _atom, x, y, z in rec:
        if resnum in hits:
            continue
        d = _np.linalg.norm(eff_xyz - _np.asarray((x, y, z), dtype=float), axis=1)
        if float(d.min()) <= cutoff:
            hits.add(resnum)
    return sorted(index_of[rn] for rn in hits)
