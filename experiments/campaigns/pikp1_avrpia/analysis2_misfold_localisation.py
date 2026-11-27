#!/usr/bin/env python3
"""
Analysis 2 — Receptor misfold localisation for RMSD-failing sequences (3a only).

For sequences where Boltz-2 cannot fold the MPNN sequence into the RFDiffusion
backbone (receptor RMSD > 6 Å): align on fixed regions only and measure per-residue
Cα displacement across the whole binder to localise the misfold.

Hypothesis A: backbone pathology — even fixed regions show large displacement after
              fixed-region alignment (the whole scaffold deforms).
Hypothesis B: sequence-specific misfold — only designed residues diverge; fixed
              regions stay close to the RFDiffusion backbone.

Method:
  1. Align Boltz-2 chain A Cα → RFD chain A Cα using ONLY fixed-region residues
  2. Compute per-residue Cα displacement for ALL chain A residues
  3. Aggregate across all RMSD-failing sequences for 3a

Usage:
  python3 analysis2_misfold_localisation.py
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

# pose_solved_3a contig: A1-6/2-5/A10-36/2-5/A40-41/1-3/A44-67/4-8/A73-78 B
FIXED_LENGTHS_3A = [6, 27, 2, 24, 6]  # A1-6, A10-36, A40-41, A44-67, A73-78


def remap(path: str) -> str:
    if path and path.startswith(HPC_PREFIX):
        return MOUNT_PREFIX + path[len(HPC_PREFIX):]
    return path


def parse_ca(pdb_path: str, chain: str) -> dict[int, np.ndarray]:
    """Return {resnum: Cα_coords} for the given chain."""
    ca = {}
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA" and line[21] == chain:
                resnum = int(line[22:26])
                ca[resnum] = np.array([float(line[30:38]), float(line[38:46]),
                                       float(line[46:54])])
    return ca


def kabsch(P: np.ndarray, Q: np.ndarray):
    """Kabsch alignment. Returns (R, t, rmsd)."""
    pc = P - P.mean(0)
    qc = Q - Q.mean(0)
    H = pc.T @ qc
    U, S, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    t = Q.mean(0) - R @ P.mean(0)
    aligned = (R @ P.T).T + t
    rmsd = math.sqrt(((aligned - Q) ** 2).sum(1).mean())
    return R, t, rmsd


def compute_residue_positions(gap_lengths: list, fixed_lengths: list):
    """Return (fixed_residues, designed_residues) in the output PDB."""
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
    """Return (initial_pdb_path, rmsd) from raw_per_seed_results.csv."""
    csv_path = os.path.join(run_dir, "results", "negative_steering", "runs",
                            seq_name, "raw_per_seed_results.csv")
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            if row.get("design") == "initial":
                return remap(row["pdb"]), float(row["steered_independent_receptor_rmsd"])
    raise ValueError(f"No initial row in {csv_path}")


def main():
    run_name = "pose_solved_3a"
    run_dir = os.path.join(REPO_HPC, "experiments", "campaigns", "pikp1_avrpia",
                           "runs", run_name)
    rfd_dir = os.path.join(run_dir, "results", "rfdiffusion")
    summary_csv = os.path.join(run_dir, "results", "negative_steering",
                               "cross_sequence_summary_with_interface_metrics.csv")

    print(f"Analysis 2: Receptor misfold localisation ({run_name})")
    print(f"RMSD fail threshold: > {RMSD_THRESHOLD} Å\n")

    # Load failing sequences
    failing_seqs = []
    with open(summary_csv) as f:
        for row in csv.DictReader(f):
            if row["row_type"] != "steered":
                continue
            rmsd = float(row["rep_independent_receptor_rmsd_median"])
            if rmsd > RMSD_THRESHOLD:
                failing_seqs.append({
                    "seq": row["mpnn_sequence"],
                    "rmsd": rmsd,
                    "designed_residues": row["designed_residues"],
                })

    print(f"Found {len(failing_seqs)} RMSD-failing sequences")

    # Accumulate per-residue displacements across sequences
    # Key: residue position label (residue number in output PDB)
    # We use relative position labels: F1..Fn for fixed, D1..Dm for designed
    # For aggregation, store raw per-residue displacements indexed by output residue number
    per_res_displacements = defaultdict(list)   # resnum → [displacement, ...]
    res_category = {}  # resnum → "fixed" or "designed"
    sequence_level = []  # per-sequence summary

    skipped = 0
    for item in failing_seqs:
        seq_name = item["seq"]
        design_id = "_".join(seq_name.split("_")[:2])
        rfd_pdb = os.path.join(rfd_dir, f"{design_id}.pdb")

        try:
            init_pdb, init_rmsd = get_initial_pdb(run_dir, seq_name)
        except Exception as e:
            skipped += 1
            continue

        if not os.path.exists(rfd_pdb) or not os.path.exists(init_pdb):
            skipped += 1
            continue

        gap_lengths = [len(s) for s in item["designed_residues"].split("|")]
        if len(gap_lengths) != len(FIXED_LENGTHS_3A) - 1:
            skipped += 1
            continue

        fixed_res, designed_res = compute_residue_positions(gap_lengths, FIXED_LENGTHS_3A)
        fixed_set = set(fixed_res)
        designed_set = set(designed_res)

        # Update category map (may overwrite but should be consistent per design)
        for r in fixed_res:
            res_category[r] = "fixed"
        for r in designed_res:
            res_category[r] = "designed"

        rfd_ca = parse_ca(rfd_pdb, "A")
        boltz_ca = parse_ca(init_pdb, "A")

        # Align on fixed-region CA only
        shared_fixed = sorted(fixed_set & set(rfd_ca.keys()) & set(boltz_ca.keys()))
        if len(shared_fixed) < 6:
            skipped += 1
            continue

        P = np.array([boltz_ca[r] for r in shared_fixed])
        Q = np.array([rfd_ca[r] for r in shared_fixed])
        R, t, fixed_rmsd = kabsch(P, Q)

        # Apply rotation/translation to all Boltz-2 CA coords
        all_res = sorted(set(rfd_ca.keys()) & set(boltz_ca.keys()))
        P_all = np.array([boltz_ca[r] for r in all_res])
        P_aligned = (R @ P_all.T).T + t

        # Per-residue Cα displacement
        Q_all = np.array([rfd_ca[r] for r in all_res])
        disp = np.sqrt(((P_aligned - Q_all) ** 2).sum(1))

        # Store
        for r, d in zip(all_res, disp):
            per_res_displacements[r].append(d)

        # Compute region-level stats for this sequence
        fixed_disp = [d for r, d in zip(all_res, disp) if r in fixed_set]
        designed_disp = [d for r, d in zip(all_res, disp) if r in designed_set]

        sequence_level.append({
            "seq": seq_name,
            "total_rmsd": item["rmsd"],
            "fixed_region_rmsd_post_align": round(fixed_rmsd, 3),
            "fixed_mean_disp": round(sum(fixed_disp) / len(fixed_disp), 3) if fixed_disp else None,
            "designed_mean_disp": round(sum(designed_disp) / len(designed_disp), 3) if designed_disp else None,
            "fixed_max_disp": round(max(fixed_disp), 3) if fixed_disp else None,
            "designed_max_disp": round(max(designed_disp), 3) if designed_disp else None,
        })

    print(f"Processed: {len(sequence_level)} sequences  (skipped: {skipped})\n")

    if not sequence_level:
        print("No results computed.")
        return

    # ── Sequence-level summary ─────────────────────────────────────────────
    all_fixed_mean = [s["fixed_mean_disp"] for s in sequence_level if s["fixed_mean_disp"] is not None]
    all_des_mean = [s["designed_mean_disp"] for s in sequence_level if s["designed_mean_disp"] is not None]
    all_fixed_max = [s["fixed_max_disp"] for s in sequence_level if s["fixed_max_disp"] is not None]
    all_des_max = [s["designed_max_disp"] for s in sequence_level if s["designed_max_disp"] is not None]

    print("Region-level Cα displacement (after fixed-region alignment):")
    print(f"  Fixed regions   — mean disp: {sum(all_fixed_mean)/len(all_fixed_mean):.2f} Å  "
          f"median: {sorted(all_fixed_mean)[len(all_fixed_mean)//2]:.2f} Å  "
          f"max: {max(all_fixed_max):.2f} Å")
    print(f"  Designed regions — mean disp: {sum(all_des_mean)/len(all_des_mean):.2f} Å  "
          f"median: {sorted(all_des_mean)[len(all_des_mean)//2]:.2f} Å  "
          f"max: {max(all_des_max):.2f} Å")

    # ── Hypothesis adjudication ────────────────────────────────────────────
    overall_fixed_mean = sum(all_fixed_mean) / len(all_fixed_mean)
    overall_des_mean = sum(all_des_mean) / len(all_des_mean)
    ratio = overall_des_mean / max(overall_fixed_mean, 0.01)

    print(f"\n  Designed/Fixed displacement ratio: {ratio:.1f}×")
    if overall_fixed_mean > 5.0:
        print("\n  HYPOTHESIS A SUPPORTED: Fixed regions also diverge significantly after")
        print("  fixed-region alignment. The RFDiffusion backbone is globally incompatible")
        print("  with the MPNN sequence — Boltz-2 predicts an entirely different fold.")
        print("  This is not fixable by resampling MPNN sequences for the same backbone.")
    elif ratio >= 3.0 and overall_fixed_mean < 3.0:
        print("\n  HYPOTHESIS B SUPPORTED: Displacement is concentrated in the designed")
        print("  regions while fixed regions stay close (< 3 Å). The MPNN sequences")
        print("  misfold at the specific de novo segments. Resampling MPNN sequences")
        print("  (different temperatures, different seeds) may find sequences that fold.")
    else:
        print("\n  MIXED SIGNAL: Both fixed and designed regions show elevated displacement.")
        print("  The backbone is partially incompatible but designed segments contribute more.")

    # ── Per-residue aggregated displacement table ──────────────────────────
    print(f"\nPer-residue mean Cα displacement (across {len(sequence_level)} failing sequences):")
    print(f"  {'Res':>4}  {'Category':>10}  {'Mean Å':>8}  {'Median Å':>10}  {'Max Å':>8}  {'N':>4}")

    all_residues = sorted(per_res_displacements.keys())
    for r in all_residues:
        disps = per_res_displacements[r]
        mean_d = sum(disps) / len(disps)
        sorted_d = sorted(disps)
        med_d = sorted_d[len(sorted_d) // 2]
        max_d = max(disps)
        cat = res_category.get(r, "?")
        flag = " *** HIGH ***" if mean_d > 10.0 else (" * elevated" if mean_d > 5.0 else "")
        print(f"  {r:>4}  {cat:>10}  {mean_d:>8.2f}  {med_d:>10.2f}  {max_d:>8.2f}  "
              f"{len(disps):>4}{flag}")

    # ── Per-sequence detail ────────────────────────────────────────────────
    print(f"\nPer-sequence summary (sorted by designed region mean displacement):")
    print(f"  {'Sequence':<25} {'TotRMSD':>8} {'FixedMean':>10} {'DesMean':>9} {'DesMean/FixMean':>15}")
    for s in sorted(sequence_level, key=lambda x: -(x["designed_mean_disp"] or 0)):
        ratio_s = (s["designed_mean_disp"] or 0) / max(s["fixed_mean_disp"] or 0.01, 0.01)
        print(f"  {s['seq']:<25} {s['total_rmsd']:>8.2f} "
              f"{s['fixed_mean_disp']:>10.3f} {s['designed_mean_disp']:>9.3f} "
              f"{ratio_s:>15.2f}×")


if __name__ == "__main__":
    main()
