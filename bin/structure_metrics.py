"""
structure_metrics.py
--------------------
Geometric / structural primitives shared across the pipeline.

Consolidates inline centre-of-mass + clash-counting code that previously
lived only in ``rfdiffusion_filter.py`` so other consumers can reuse
the same primitives without duplication.

Contract
========
- All functions are pure (no globals, no I/O side effects beyond reading
  the path passed in).
- Functions that take a PDB read the first model only and consider ATOM
  records exclusively (HETATMs ignored).
- Insertion codes are NOT distinguished — two residues with the same
  resseq but different icodes collapse into one CA pick.  This matches
  rfdiffusion_filter's pre-existing behaviour; refine here if a future
  consumer needs icode-aware logic.

Public API
==========
- ``centroid(coords)`` — numpy.mean along axis 0 with NaN handling.
- ``centroid_distance(c1, c2)`` — Euclidean distance between two centroids.
- ``read_chain_ca_coords(pdb_path, chain_id)`` — Nx3 numpy array of CA
  coordinates for one chain, in residue order.
- ``chain_centre_of_mass(pdb_path, chain_id)`` — convenience composing the
  above two.  Returns ``None`` if the chain has no CAs.
- ``chain_com_distance(pdb_path, chain_a, chain_b)`` — convenience composing
  the above.
- ``read_chain_heavy_atoms(pdb_path, chain_id)`` — list of
  ``(resnum, atom_name, x, y, z)`` for all non-hydrogen ATOM records.
- ``clash_count(pdb_path, chain_a, chain_b, design_region, cutoff=2.0)`` —
  returns ``(in_design, outside_design)`` heavy-atom pair counts under the
  cutoff.  ``design_region`` is a set of receptor residue numbers; a pair
  counts as "in design" iff EITHER atom belongs to a residue in that set.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, List, Optional, Set, Tuple

import numpy as np

# ── Numpy-level primitives ──────────────────────────────────────────


def centroid(coords: np.ndarray) -> np.ndarray:
    """Return ``coords.mean(axis=0)`` as a 1D length-3 array.

    Empty coords or coords containing NaN raise ValueError so callers
    don't silently get NaN centroids and propagate them.
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

    Stops at ENDMDL so multi-model PDBs (e.g. NMR ensembles, docking
    cluster outputs) only emit the first model.
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

    Note: this is geometric centroid of CAs, not mass-weighted COM.  For
    the use cases here (compare two domains, distance between chains)
    the distinction doesn't matter and CA-only is robust to side-chain
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
    """Count heavy-atom pairs (one atom from chain_a, one from chain_b)
    within ``cutoff`` Å of each other.  Splits the count by whether the
    chain_a residue is in the design region.

    Returns ``(clashes_in_design_region, clashes_outside_design_region)``.

    A pair counts as "in design region" iff the chain_a (receptor)
    residue number is in ``design_region_residues_a``.  Chain B atoms
    don't have a notion of design region at this stage.

    Note: this uses a naive O(N*M) loop.  For typical 2-chain interface
    PDBs (one receptor ~100 residues × one effector ~100 residues)
    that's ~10^4 atom pairs per chain pair — fast enough.  If a future
    consumer has much larger inputs, swap in a KDTree.
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
# Both functions below were this repository's only reason to import from
# bin/boltz2_negative_steering.py and bin/compute_metrics.py, which belong to
# the negative-steering engine and are installed separately.  They are
# reimplemented here on this module's own PDB reader rather than copied, so
# nothing in this repository depends on engine internals.
#
# tests/characterization/test_structure_metrics_ported.py checks both against
# the engine's originals on every reference PDB in the test fixtures.

def get_chain_sequence(pdb_path: Path, chain_id: str) -> str:
    """One-letter sequence of `chain_id`, in the order the file lists it.

    One residue per (resnum, icode), taken from the first ATOM record of each,
    so an altloc does not duplicate a position.  A residue whose three-letter
    name is not standard becomes 'X', which keeps the sequence the same length
    as the residue list and so keeps every positional index downstream valid.

    Raises ValueError when the chain is absent, naming the chains that are
    present, because the alternative is a silent empty string that reads as an
    empty chain.
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
    heavy atom of the effector chain.  Hydrogens are excluded by
    read_chain_heavy_atoms, which is why the cutoff can be as tight as 5 A.

    The index is positional within the receptor chain, counted over
    (resnum, icode) in file order, so it lines up with get_chain_sequence above
    and with the design-region and true-interface files, which are also
    positional.  Returning residue numbers instead would break on any structure
    that does not start at 1.
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
