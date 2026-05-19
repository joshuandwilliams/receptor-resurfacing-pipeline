"""
interface_metrics.py

Compute interface metrics for a two-chain protein complex.

Metrics produced:
  - CA-CA contact pairs at <=8 A (count + distance distribution stats)
  - Minimum heavy-atom contact pairs at <=5 A (count + distance distribution)
  - Heavy-atom clash counts at thresholds 2.0, 2.5, 3.0 A
  - Buried surface area (BSA) via Shrake-Rupley SASA difference (complex vs separated)
  - Interface residue count per chain (any heavy atom within 5 A of other chain)
  - Interface composition: hydrophobic / polar / charged / aromatic fractions per chain
  - Gap-index proxy for shape complementarity:
      mean nearest-other-chain distance over interface heavy atoms (lower = tighter packing)
  - Hydrogen-bond candidate count (donor->acceptor heavy-atom distance <=3.5 A)

CLI:
  python interface_metrics.py <pdb> [--chain-a A] [--chain-b B] [--json-out OUT.json]
  python interface_metrics.py --batch <pdb1> <pdb2> ...  [--chains A:B A:C ...]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np

try:
    from Bio.PDB import PDBParser
    from Bio.PDB.SASA import ShrakeRupley
except ImportError as exc:  # pragma: no cover
    sys.exit(f"biopython is required: {exc}")


HYDROPHOBIC = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "GLY", "CYS"}
POLAR       = {"SER", "THR", "ASN", "GLN", "TYR", "HIS"}
CHARGED     = {"ASP", "GLU", "LYS", "ARG"}
AROMATIC    = {"PHE", "TYR", "TRP", "HIS"}

HBOND_DONOR_ATOMS    = {"N", "ND1", "ND2", "NE", "NE1", "NE2", "NH1", "NH2", "NZ", "OG", "OG1", "OH", "SG"}
HBOND_ACCEPTOR_ATOMS = {"O", "OD1", "OD2", "OE1", "OE2", "OG", "OG1", "OH", "ND1", "NE2", "SD"}


def _load_chains(pdb_path: Path, chain_a: str, chain_b: str):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_path.stem, str(pdb_path))
    model = next(structure.get_models())
    chains = {c.id: c for c in model.get_chains()}
    if chain_a not in chains or chain_b not in chains:
        raise ValueError(
            f"{pdb_path.name}: requested chains {chain_a}/{chain_b}, found {sorted(chains)}"
        )
    return structure, model, chains[chain_a], chains[chain_b]


def _heavy_atoms(chain) -> List:
    atoms = []
    for res in chain.get_residues():
        if res.id[0] != " ":
            continue
        for atom in res.get_atoms():
            if atom.element != "H":
                atoms.append(atom)
    return atoms


def _ca_atoms(chain) -> List:
    out = []
    for res in chain.get_residues():
        if res.id[0] != " ":
            continue
        if "CA" in res:
            out.append(res["CA"])
    return out


def _pairwise(coords_a: np.ndarray, coords_b: np.ndarray) -> np.ndarray:
    diff = coords_a[:, None, :] - coords_b[None, :, :]
    return np.sqrt((diff ** 2).sum(-1))


def _dist_stats(values: np.ndarray) -> Dict[str, float]:
    if values.size == 0:
        return {"n": 0}
    return {
        "n": int(values.size),
        "min": float(values.min()),
        "max": float(values.max()),
        "mean": float(values.mean()),
        "median": float(np.median(values)),
        "p25": float(np.percentile(values, 25)),
        "p75": float(np.percentile(values, 75)),
    }


def _bin_distribution(values: np.ndarray, edges: Sequence[float]) -> Dict[str, int]:
    out = {}
    for lo, hi in zip(edges[:-1], edges[1:]):
        out[f"{lo:.1f}-{hi:.1f}"] = int(((values >= lo) & (values < hi)).sum())
    return out


def _residue_composition(residue_ids: List[Tuple[str, int]]) -> Dict[str, float]:
    if not residue_ids:
        return {"hydrophobic": 0.0, "polar": 0.0, "charged": 0.0, "aromatic": 0.0}
    n = len(residue_ids)
    counts = {"hydrophobic": 0, "polar": 0, "charged": 0, "aromatic": 0}
    for resname, _ in residue_ids:
        if resname in HYDROPHOBIC:
            counts["hydrophobic"] += 1
        if resname in POLAR:
            counts["polar"] += 1
        if resname in CHARGED:
            counts["charged"] += 1
        if resname in AROMATIC:
            counts["aromatic"] += 1
    return {k: v / n for k, v in counts.items()}


def compute_metrics(pdb_path: Path, chain_a: str, chain_b: str) -> Dict:
    structure, model, ca, cb = _load_chains(pdb_path, chain_a, chain_b)

    ca_atoms_a = _ca_atoms(ca)
    ca_atoms_b = _ca_atoms(cb)
    heavy_a = _heavy_atoms(ca)
    heavy_b = _heavy_atoms(cb)

    ca_coords_a = np.array([a.coord for a in ca_atoms_a])
    ca_coords_b = np.array([a.coord for a in ca_atoms_b])
    heavy_coords_a = np.array([a.coord for a in heavy_a])
    heavy_coords_b = np.array([a.coord for a in heavy_b])

    # --- CA-CA contacts ----------------------------------------------------
    d_ca = _pairwise(ca_coords_a, ca_coords_b)
    ca_contact_mask = d_ca <= 8.0
    ca_contact_distances = d_ca[ca_contact_mask]

    # --- heavy-atom contacts / clashes ------------------------------------
    d_heavy = _pairwise(heavy_coords_a, heavy_coords_b)
    contact_mask_5a = d_heavy <= 5.0
    contact_distances_5a = d_heavy[contact_mask_5a]

    clash_counts = {
        f"<{thr:.1f}A": int((d_heavy < thr).sum())
        for thr in (2.0, 2.5, 3.0, 3.5)
    }

    # --- per-residue interface identification -----------------------------
    interface_res_a: Dict[Tuple[str, int], str] = {}
    interface_res_b: Dict[Tuple[str, int], str] = {}
    contact_pairs = np.argwhere(d_heavy <= 5.0)
    for i, j in contact_pairs:
        ra = heavy_a[i].get_parent()
        rb = heavy_b[j].get_parent()
        interface_res_a[(ra.get_resname(), ra.id[1])] = ca.id
        interface_res_b[(rb.get_resname(), rb.id[1])] = cb.id

    comp_a = _residue_composition(list(interface_res_a.keys()))
    comp_b = _residue_composition(list(interface_res_b.keys()))

    # --- BSA via Shrake-Rupley -------------------------------------------
    sr = ShrakeRupley(probe_radius=1.4, n_points=960)
    sr.compute(structure, level="A")
    sasa_complex_a = sum(a.sasa for a in heavy_a)
    sasa_complex_b = sum(a.sasa for a in heavy_b)

    # Recompute SASA after detaching each chain in turn to get free SASA
    # We use temporary structures to avoid mutating the model in place.
    def _chain_only_sasa(target_chain) -> float:
        from Bio.PDB.Structure import Structure
        from Bio.PDB.Model import Model
        from Bio.PDB.Chain import Chain
        s = Structure(target_chain.id + "_iso")
        m = Model(0)
        s.add(m)
        new_chain = Chain(target_chain.id)
        for res in target_chain.get_residues():
            new_chain.add(res.copy())
        m.add(new_chain)
        sr_iso = ShrakeRupley(probe_radius=1.4, n_points=960)
        sr_iso.compute(s, level="A")
        return float(sum(a.sasa for a in new_chain.get_atoms() if a.element != "H"))

    sasa_free_a = _chain_only_sasa(ca)
    sasa_free_b = _chain_only_sasa(cb)
    bsa_a = sasa_free_a - sasa_complex_a
    bsa_b = sasa_free_b - sasa_complex_b
    bsa_total = bsa_a + bsa_b

    # --- shape-complementarity proxy --------------------------------------
    # For each interface heavy atom in A, distance to nearest heavy atom in B.
    # Lower mean = tighter packing; very low minima paired with low mean ->
    # well-mated interface.
    interface_atom_mask_a = (d_heavy <= 5.0).any(axis=1)
    interface_atom_mask_b = (d_heavy <= 5.0).any(axis=0)
    nearest_a_to_b = d_heavy[interface_atom_mask_a].min(axis=1)
    nearest_b_to_a = d_heavy[:, interface_atom_mask_b].min(axis=0)
    gap_index = float(
        np.concatenate([nearest_a_to_b, nearest_b_to_a]).mean()
        if nearest_a_to_b.size and nearest_b_to_a.size
        else float("nan")
    )

    # --- hydrogen-bond candidate count ------------------------------------
    def _atom_pairs(donors, acceptors) -> int:
        if not donors or not acceptors:
            return 0
        d_coords = np.array([a.coord for a in donors])
        a_coords = np.array([a.coord for a in acceptors])
        dd = _pairwise(d_coords, a_coords)
        return int((dd <= 3.5).sum())

    donors_a    = [a for a in heavy_a if a.get_name() in HBOND_DONOR_ATOMS]
    accept_a    = [a for a in heavy_a if a.get_name() in HBOND_ACCEPTOR_ATOMS]
    donors_b    = [a for a in heavy_b if a.get_name() in HBOND_DONOR_ATOMS]
    accept_b    = [a for a in heavy_b if a.get_name() in HBOND_ACCEPTOR_ATOMS]
    hbond_a_to_b = _atom_pairs(donors_a, accept_b)
    hbond_b_to_a = _atom_pairs(donors_b, accept_a)

    return {
        "pdb": str(pdb_path),
        "chains": {"a": chain_a, "b": chain_b},
        "size": {
            "residues_a": len(ca_atoms_a),
            "residues_b": len(ca_atoms_b),
            "heavy_atoms_a": len(heavy_a),
            "heavy_atoms_b": len(heavy_b),
        },
        "ca_contacts_8A": {
            "count": int(ca_contact_mask.sum()),
            "distance_stats": _dist_stats(ca_contact_distances),
            "distance_histogram": _bin_distribution(
                ca_contact_distances, [0, 4, 5, 6, 7, 8]
            ),
        },
        "heavy_contacts_5A": {
            "count": int(contact_mask_5a.sum()),
            "distance_stats": _dist_stats(contact_distances_5a),
            "distance_histogram": _bin_distribution(
                contact_distances_5a, [0, 2, 2.5, 3, 3.5, 4, 4.5, 5]
            ),
        },
        "clash_counts": clash_counts,
        "interface_residues": {
            "count_a": len(interface_res_a),
            "count_b": len(interface_res_b),
            "composition_a": comp_a,
            "composition_b": comp_b,
        },
        "bsa_A2": {
            "chain_a": float(bsa_a),
            "chain_b": float(bsa_b),
            "total": float(bsa_total),
        },
        "shape_complementarity_proxy": {
            "gap_index_mean_A": gap_index,
            "interface_heavy_atoms_a": int(interface_atom_mask_a.sum()),
            "interface_heavy_atoms_b": int(interface_atom_mask_b.sum()),
        },
        "hbond_candidates_3p5A": {
            "donor_a_to_acceptor_b": hbond_a_to_b,
            "donor_b_to_acceptor_a": hbond_b_to_a,
            "total": hbond_a_to_b + hbond_b_to_a,
        },
    }


def _pretty(metrics: Dict) -> str:
    lines = []
    p = metrics
    lines.append(f"=== {Path(p['pdb']).name}  chains {p['chains']['a']}/{p['chains']['b']} ===")
    lines.append(
        f"  size:  A {p['size']['residues_a']} res / {p['size']['heavy_atoms_a']} heavy   "
        f"B {p['size']['residues_b']} res / {p['size']['heavy_atoms_b']} heavy"
    )
    ca = p["ca_contacts_8A"]
    lines.append(
        f"  CA contacts <=8A: {ca['count']}  "
        f"(mean {ca['distance_stats'].get('mean', float('nan')):.2f} A, "
        f"median {ca['distance_stats'].get('median', float('nan')):.2f} A)"
    )
    lines.append(f"    histogram: {ca['distance_histogram']}")
    hv = p["heavy_contacts_5A"]
    lines.append(
        f"  heavy contacts <=5A: {hv['count']}  "
        f"(mean {hv['distance_stats'].get('mean', float('nan')):.2f} A)"
    )
    lines.append(f"    histogram: {hv['distance_histogram']}")
    lines.append(f"  clashes: {p['clash_counts']}")
    ir = p["interface_residues"]
    lines.append(
        f"  interface residues:  A {ir['count_a']}  B {ir['count_b']}"
    )
    lines.append(
        f"    composition A: hydroph {ir['composition_a']['hydrophobic']:.2f}  "
        f"polar {ir['composition_a']['polar']:.2f}  "
        f"charged {ir['composition_a']['charged']:.2f}  "
        f"aromatic {ir['composition_a']['aromatic']:.2f}"
    )
    lines.append(
        f"    composition B: hydroph {ir['composition_b']['hydrophobic']:.2f}  "
        f"polar {ir['composition_b']['polar']:.2f}  "
        f"charged {ir['composition_b']['charged']:.2f}  "
        f"aromatic {ir['composition_b']['aromatic']:.2f}"
    )
    bsa = p["bsa_A2"]
    lines.append(
        f"  BSA: total {bsa['total']:.0f} A^2  (A {bsa['chain_a']:.0f}  B {bsa['chain_b']:.0f})"
    )
    sc = p["shape_complementarity_proxy"]
    lines.append(
        f"  gap-index (mean nearest-other-chain d, interface atoms): {sc['gap_index_mean_A']:.3f} A   "
        f"(A iface atoms {sc['interface_heavy_atoms_a']}, B iface atoms {sc['interface_heavy_atoms_b']})"
    )
    hb = p["hbond_candidates_3p5A"]
    lines.append(
        f"  H-bond candidates (<=3.5A, geom. only): A->B {hb['donor_a_to_acceptor_b']}  "
        f"B->A {hb['donor_b_to_acceptor_a']}  total {hb['total']}"
    )
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pdb", nargs="+", help="PDB file(s)")
    ap.add_argument("--chains", nargs="+", default=None,
                    help="Per-PDB chain pair like A:B  (one per PDB). Defaults to A:B for all.")
    ap.add_argument("--json-out", default=None,
                    help="If set, write full metrics to this JSON file.")
    args = ap.parse_args()

    if args.chains and len(args.chains) != len(args.pdb):
        sys.exit("--chains must have same length as pdb list")
    chain_pairs = args.chains or ["A:B"] * len(args.pdb)

    all_metrics = []
    for pdb_str, chain_pair in zip(args.pdb, chain_pairs):
        a, b = chain_pair.split(":")
        m = compute_metrics(Path(pdb_str), a, b)
        all_metrics.append(m)
        print(_pretty(m))
        print()

    if args.json_out:
        Path(args.json_out).write_text(json.dumps(all_metrics, indent=2))
        print(f"Wrote {args.json_out}")


if __name__ == "__main__":
    main()
