#!/usr/bin/env python3
"""
Analysis 1 — Sidechain clash check for RMSD-passing sequences.

For sequences where Boltz-2 correctly reproduces the RFDiffusion binder backbone
(receptor RMSD ≤ 6 Å) but Ra_eff is still high (effector not found correctly):
check whether MPNN-assigned binder sidechains sterically clash with the effector
when both are placed in their RFDiffusion-intended positions.

Method:
  1. Align Boltz-2 chain A Cα → RFD chain A Cα  (binder in correct position)
  2. Align Boltz-2 chain B Cα → RFD chain B Cα  (effector in correct position)
  3. Count heavy-atom clashes (< 2.0 Å) between chain A and chain B atoms

Both runs are processed: pose_solved_3a (12 passing) and pose_solved_5a (47 passing).

Usage:
  python3 analysis1_sidechain_clashes.py
"""

import csv
import math
import os
from collections import defaultdict

import numpy as np

# ── Path configuration ──────────────────────────────────────────────────────
HPC_PREFIX = "/hpc-home/jowillia/"
MOUNT_PREFIX = "/Volumes/HPC-Home/"
REPO_HPC = MOUNT_PREFIX + "receptor_design/receptor-resurfacing-pipeline"

RMSD_THRESHOLD = 6.0
CLASH_DIST = 2.0  # Å — heavy-atom clash threshold

# Contig-derived fixed-segment lengths for each run
# Structure: alternating fixed / gap / fixed / gap / ... / fixed
# fixed_lengths[i] = length of i-th fixed block
# gap lengths come from designed_residues in the summary CSV
RUN_CONFIG = {
    "pose_solved_3a": {
        "fixed_lengths": [6, 27, 2, 24, 6],   # A1-6, A10-36, A40-41, A44-67, A73-78
    },
    "pose_solved_5a": {
        "fixed_lengths": [4, 25, 8, 12, 2, 2], # A1-4, A11-35, A46-53, A55-66, A74-75, A77-78
    },
}


def remap(path: str) -> str:
    if path and path.startswith(HPC_PREFIX):
        return MOUNT_PREFIX + path[len(HPC_PREFIX):]
    return path


def parse_atoms(pdb_path: str, chain: str) -> dict:
    """Return {resnum: {atom_name: np.array([x,y,z])}} for the given chain."""
    atoms = defaultdict(dict)
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            if line[21] != chain:
                continue
            name = line[12:16].strip()
            resnum = int(line[22:26])
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            atoms[resnum][name] = np.array([x, y, z])
    return dict(atoms)


def kabsch(P: np.ndarray, Q: np.ndarray):
    """
    Kabsch algorithm: find rotation R and translation t such that
    R @ P[i] + t ≈ Q[i].  Returns (R, t, rmsd).
    P and Q are (N, 3) arrays of corresponding points.
    """
    pc = P - P.mean(0)
    qc = Q - Q.mean(0)
    H = pc.T @ qc
    U, S, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    t = Q.mean(0) - R @ P.mean(0)
    aligned = (R @ P.T).T + t
    rmsd = math.sqrt(((aligned - Q) ** 2).sum(1).mean())
    return R, t, rmsd


def apply_transform(atoms: dict, R: np.ndarray, t: np.ndarray) -> dict:
    """Apply rotation R and translation t to all atoms in the dict."""
    result = {}
    for resnum, res_atoms in atoms.items():
        result[resnum] = {
            name: R @ coord + t for name, coord in res_atoms.items()
        }
    return result


def get_ca_array(atoms: dict, residues: list) -> np.ndarray:
    """Extract CA coordinates for the given residue list (ordered)."""
    return np.array([atoms[r]["CA"] for r in residues if r in atoms and "CA" in atoms[r]])


def get_ca_residues(atoms: dict, residues: list) -> list:
    """Return residues that have CA in atoms dict."""
    return [r for r in residues if r in atoms and "CA" in atoms[r]]


def compute_residue_positions(gap_lengths: list, fixed_lengths: list):
    """
    Compute fixed and designed residue numbers in the output PDB.
    Structure: fixed_0, gap_0, fixed_1, gap_1, ..., fixed_n
    Returns: (fixed_residues, designed_residues) as sorted lists.
    """
    fixed, designed = [], []
    pos = 1
    for i, fl in enumerate(fixed_lengths):
        fixed.extend(range(pos, pos + fl))
        pos += fl
        if i < len(gap_lengths):
            gl = gap_lengths[i]
            designed.extend(range(pos, pos + gl))
            pos += gl
    return fixed, designed


def get_initial_pdb(run_dir: str, seq_name: str) -> tuple[str, float]:
    """Read raw_per_seed_results.csv and return (initial_pdb_path, rmsd)."""
    csv_path = os.path.join(run_dir, "results", "negative_steering", "runs",
                            seq_name, "raw_per_seed_results.csv")
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            if row.get("design") == "initial":
                pdb = remap(row["pdb"])
                rmsd = float(row["steered_independent_receptor_rmsd"])
                return pdb, rmsd
    raise ValueError(f"No initial row in {csv_path}")


def clash_check(boltz_a: dict, boltz_b: dict) -> list[tuple]:
    """
    Find heavy-atom clashes (< CLASH_DIST Å) between chain A and chain B atoms.
    Returns list of (resnum_a, atom_a, resnum_b, atom_b, distance).
    Uses numpy broadcasting for speed (~300K pairs per sequence).
    """
    a_coords, a_labels = [], []
    for ra, atoms in boltz_a.items():
        for an, ac in atoms.items():
            a_coords.append(ac)
            a_labels.append((ra, an))

    b_coords, b_labels = [], []
    for rb, atoms in boltz_b.items():
        for bn, bc in atoms.items():
            b_coords.append(bc)
            b_labels.append((rb, bn))

    A = np.array(a_coords)  # (N_A, 3)
    B = np.array(b_coords)  # (N_B, 3)
    dists = np.sqrt(((A[:, None, :] - B[None, :, :]) ** 2).sum(-1))  # (N_A, N_B)
    pairs = np.argwhere(dists < CLASH_DIST)
    return [(a_labels[i][0], a_labels[i][1], b_labels[j][0], b_labels[j][1],
             round(float(dists[i, j]), 3)) for i, j in pairs]


def analyse_run(run_name: str):
    run_dir = os.path.join(REPO_HPC, "experiments", "campaigns", "pikp1_avrpia",
                           "runs", run_name)
    rfd_dir = os.path.join(run_dir, "results", "rfdiffusion")
    summary_csv = os.path.join(run_dir, "results", "negative_steering",
                               "cross_sequence_summary_with_interface_metrics.csv")
    fixed_lengths = RUN_CONFIG[run_name]["fixed_lengths"]

    passing_seqs = []
    with open(summary_csv) as f:
        for row in csv.DictReader(f):
            if row["row_type"] != "steered":
                continue
            rmsd = float(row["rep_independent_receptor_rmsd_median"])
            if rmsd <= RMSD_THRESHOLD:
                passing_seqs.append({
                    "seq": row["mpnn_sequence"],
                    "rmsd": rmsd,
                    "ra_eff": float(row["rep_ra_eff_vs_truth_median"]),
                    "designed_residues": row["designed_residues"],
                })

    print(f"\n{'='*60}")
    print(f"Run: {run_name}  ({len(passing_seqs)} RMSD-passing sequences)")
    print(f"{'='*60}")

    results = []
    for item in passing_seqs:
        seq_name = item["seq"]
        design_id = "_".join(seq_name.split("_")[:2])  # e.g. design_47
        rfd_pdb = os.path.join(rfd_dir, f"{design_id}.pdb")

        try:
            init_pdb, init_rmsd = get_initial_pdb(run_dir, seq_name)
        except Exception as e:
            print(f"  SKIP {seq_name}: {e}")
            continue

        if not os.path.exists(rfd_pdb):
            print(f"  SKIP {seq_name}: RFD PDB missing at {rfd_pdb}")
            continue
        if not os.path.exists(init_pdb):
            print(f"  SKIP {seq_name}: Boltz-2 PDB missing at {init_pdb}")
            continue

        # Parse gap lengths from designed_residues (e.g. "GVG|SPD|RQ|ANYGY")
        gap_lengths = [len(s) for s in item["designed_residues"].split("|")]
        if len(gap_lengths) != len(fixed_lengths) - 1:
            print(f"  WARN {seq_name}: expected {len(fixed_lengths)-1} gaps, "
                  f"got {len(gap_lengths)} from designed_residues={item['designed_residues']!r}")
            continue

        # Compute fixed/designed residue positions
        fixed_res, designed_res = compute_residue_positions(gap_lengths, fixed_lengths)

        # Load PDBs
        rfd_a = parse_atoms(rfd_pdb, "A")
        rfd_b = parse_atoms(rfd_pdb, "B")
        boltz_a = parse_atoms(init_pdb, "A")
        boltz_b = parse_atoms(init_pdb, "B")

        # Step 1: Align Boltz-2 chain A → RFD chain A on all common CA
        common_a = sorted(set(rfd_a.keys()) & set(boltz_a.keys()))
        if len(common_a) < 6:
            print(f"  SKIP {seq_name}: only {len(common_a)} common chain-A residues")
            continue

        P_a = np.array([boltz_a[r]["CA"] for r in common_a])
        Q_a = np.array([rfd_a[r]["CA"] for r in common_a])
        R1, t1, rmsd1 = kabsch(P_a, Q_a)

        # Apply to all Boltz-2 atoms (chain A and B together)
        boltz_a_t = apply_transform(boltz_a, R1, t1)
        boltz_b_t = apply_transform(boltz_b, R1, t1)

        # Step 2: Align Boltz-2 chain B (now rotated) → RFD chain B on CA
        rfd_b_res = sorted(rfd_b.keys())
        boltz_b_res = sorted(boltz_b_t.keys())
        n_effector = min(len(rfd_b_res), len(boltz_b_res))
        if n_effector < 6:
            print(f"  SKIP {seq_name}: only {n_effector} effector residues")
            continue

        rfd_b_aligned = rfd_b_res[:n_effector]
        boltz_b_aligned = boltz_b_res[:n_effector]
        P_b = np.array([boltz_b_t[r]["CA"] for r in boltz_b_aligned])
        Q_b = np.array([rfd_b[r]["CA"] for r in rfd_b_aligned])
        R2, t2, rmsd2 = kabsch(P_b, Q_b)

        # Apply R2/t2 only to Boltz-2 chain B
        boltz_b_final = apply_transform(boltz_b_t, R2, t2)

        # Step 3: Count heavy-atom clashes between Boltz-2 chain A and chain B
        # (both now in their RFD-intended positions)
        clashes = clash_check(boltz_a_t, boltz_b_final)

        # Summarise clashing binder residues
        clashing_binder_res = sorted(set(c[0] for c in clashes))
        is_designed = set(designed_res)
        designed_clashes = [r for r in clashing_binder_res if r in is_designed]
        fixed_clashes = [r for r in clashing_binder_res if r not in is_designed]

        results.append({
            "seq": seq_name,
            "rmsd": item["rmsd"],
            "ra_eff": item["ra_eff"],
            "rmsd_A_align": round(rmsd1, 3),
            "rmsd_B_align": round(rmsd2, 3),
            "n_clashes": len(clashes),
            "n_clashing_binder_res": len(clashing_binder_res),
            "designed_clashing_res": designed_clashes,
            "fixed_clashing_res": fixed_clashes,
        })

    # ── Print results ──────────────────────────────────────────────────────
    if not results:
        print("  No results computed.")
        return

    n_with_clashes = sum(1 for r in results if r["n_clashes"] > 0)
    print(f"\nClash summary ({run_name}):")
    print(f"  Sequences analysed:      {len(results)}")
    print(f"  With ≥1 clash:           {n_with_clashes} "
          f"({100*n_with_clashes//len(results)}%)")
    print(f"  Mean clash count (all):  {sum(r['n_clashes'] for r in results)/len(results):.1f}")
    print(f"  Mean clash count (if any): "
          f"{sum(r['n_clashes'] for r in results if r['n_clashes']>0) / max(n_with_clashes,1):.1f}")

    # Per-residue clash frequency
    binder_clash_freq = defaultdict(int)
    for r in results:
        for res in r["designed_clashing_res"] + r["fixed_clashing_res"]:
            binder_clash_freq[res] += 1
    if binder_clash_freq:
        top = sorted(binder_clash_freq.items(), key=lambda x: -x[1])[:10]
        print(f"\n  Most frequently clashing binder residues (position, count):")
        for res, cnt in top:
            in_des = res in set(sum((r["designed_clashing_res"] for r in results), []))
            tag = "designed" if res in set(sum((r["designed_clashing_res"] for r in results if res in r["designed_clashing_res"]), [])) else "fixed"
            print(f"    residue {res:3d} ({tag:8s}): {cnt}/{len(results)} sequences")

    print(f"\n  Per-sequence details (RMSD ≤ 6 Å sequences):")
    print(f"  {'Sequence':<25} {'RecRMSD':>8} {'Ra_eff':>7} {'A-align':>8} {'B-align':>8} "
          f"{'Clashes':>8} {'Clashing_res'}")
    for r in sorted(results, key=lambda x: -x["n_clashes"]):
        all_clashing = sorted(set(r["designed_clashing_res"] + r["fixed_clashing_res"]))
        print(f"  {r['seq']:<25} {r['rmsd']:>8.2f} {r['ra_eff']:>7.2f} "
              f"{r['rmsd_A_align']:>8.3f} {r['rmsd_B_align']:>8.3f} "
              f"{r['n_clashes']:>8d}  {all_clashing}")

    return results


def main():
    print("Analysis 1: Sidechain clash check for RMSD-passing sequences")
    print(f"Clash distance threshold: {CLASH_DIST} Å")
    print(f"RMSD pass threshold: {RMSD_THRESHOLD} Å")

    all_results = {}
    for run_name in ["pose_solved_3a", "pose_solved_5a"]:
        all_results[run_name] = analyse_run(run_name)

    # Cross-run summary
    print("\n" + "=" * 60)
    print("Cross-run interpretation")
    print("=" * 60)
    for run_name, results in all_results.items():
        if not results:
            continue
        n_with = sum(1 for r in results if r["n_clashes"] > 0)
        n_sidechain = sum(1 for r in results if r["n_clashing_binder_res"] > 0)
        print(f"\n{run_name}: {n_with}/{len(results)} sequences have sidechain clashes "
              f"at the intended interface")
        if n_with / max(len(results), 1) > 0.5:
            print("  → MPNN sidechain clashes are a COMMON feature of Level-2 failures")
            print("    in this run. Redesigning MPNN at the interface residues may improve")
            print("    effector placement.")
        else:
            print("  → Sidechain clashes are NOT the dominant cause of Level-2 failures.")
            print("    The effector displacement is likely from backbone geometry mismatch,")
            print("    not steric exclusion by binder sidechains.")


if __name__ == "__main__":
    main()
