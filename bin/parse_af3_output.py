#!/usr/bin/env python3
"""
parse_af3_output.py
-------------------
P0-31 · Post-process a single AF3-no-MSA prediction directory.

For each (seed × diffusion-sample) mmCIF written by AF3, compute the
receptor-aligned effector RMSD against the RFDiffusion reference PDB.
Emit a one-row CSV per survivor keyed by seq_name.

Output columns:
  seq_name
  af3_nomsa_best_ra_eff         Minimum across all seed×sample predictions.
  af3_nomsa_mean_ra_eff         Mean across all predictions.
  af3_nomsa_best_iptm           Max ipTM across predictions (from AF3
                                summary_confidences.json).
  af3_nomsa_mean_iptm           Mean ipTM.
  af3_nomsa_n_correct_interface Count of predictions where ra_eff
                                < ra_eff_threshold (3 seeds × 5 samples
                                = 15 by default; 1 pass fails the filter).
  af3_nomsa_total_predictions   Total number of predictions parsed.
  af3_nomsa_failures            Comma-separated failure reasons, or empty.

AF3 output layout (NBI build, Dec 2024):
  output/<name>/
      <name>_model.cif                  top-ranked structure
      seed-<N>_sample-<M>/<name>_model.cif   all seed × sample predictions
      summary_confidences.json          per-prediction ipTM/pTM/ranking_score
      ...
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np


# ═════════════════════════════════════════════════════════════════════════
# Structure loading and RMSD
# ═════════════════════════════════════════════════════════════════════════

def _chain_ca_coords(structure, chain_id: str):
    """Return an ordered list of CA (x, y, z) coords for the given chain,
    in the order they appear in the structure (NOT keyed by residue
    number).

    Why ordered list instead of {resnum: xyz} dict:
      Different sources use different numbering conventions for chain B.
      RFDiffusion / Boltz / our cleaned ground-truth PDBs typically use
      CONTINUOUS numbering across chains (chain A is 1-82, chain B is
      83-164).  AF3 always emits per-chain numbering starting from 1
      (chain A is 1-82, chain B is 1-82).  Intersecting by residue
      number gives an empty (or near-empty) effector overlap and the
      RMSD computation fails with `too_few_effector_residues:0`.

      Position-pairing is correct because the effector sequence is
      identical between pred and ref (the manifest extracts the same
      effector that was passed to AF3 as input), and the receptor
      sequence in both pred and ref is the same length (steering
      mutates identities, not lengths).  So pred[i] should pair with
      ref[i] for both chains.

      This was the cause of the 2026-04-29 'half the cohort missing
      AF3 metrics' incident — design_19 happened to have 5 residues
      of overlap by coincidence (78-82) so it slipped through with a
      barely-meaningful RMSD; design_20's chain B numbering shifted
      one position further and the intersection was empty.  The bare
      length check below now catches this consistently.
    """
    coords = []
    for model in structure:
        for chain in model:
            if chain.name != chain_id:
                continue
            for res in chain:
                for atom in res:
                    if atom.name == "CA":
                        coords.append(np.array(
                            [atom.pos.x, atom.pos.y, atom.pos.z]
                        ))
                        break
    return coords


def _receptor_aligned_effector_rmsd(
    pred_cif_path: Path,
    ref_pdb_path: Path,
    receptor_chain: str,
    effector_chain: str,
) -> Tuple[Optional[float], Optional[str]]:
    """Return receptor-aligned effector Cα RMSD or (None, error)."""
    import gemmi

    try:
        pred = gemmi.read_structure(str(pred_cif_path))
        ref = gemmi.read_structure(str(ref_pdb_path))
    except Exception as e:  # noqa: BLE001
        return None, f"structure_load_failed:{type(e).__name__}"

    try:
        pred_rec = _chain_ca_coords(pred, receptor_chain)
        pred_eff = _chain_ca_coords(pred, effector_chain)
        ref_rec  = _chain_ca_coords(ref,  receptor_chain)
        ref_eff  = _chain_ca_coords(ref,  effector_chain)
    except Exception as e:  # noqa: BLE001
        return None, f"chain_parse_failed:{type(e).__name__}"

    # Pair by position up to the shorter of the two chains.  The
    # sequences should match in length for the effector (identical
    # input) and the receptor (steering doesn't alter length); a
    # mismatch indicates an upstream issue worth flagging.
    rec_n = min(len(pred_rec), len(ref_rec))
    eff_n = min(len(pred_eff), len(ref_eff))
    if rec_n < 3:
        return None, (f"too_few_receptor_residues:{rec_n} "
                      f"(pred={len(pred_rec)}, ref={len(ref_rec)})")
    if len(pred_rec) != len(ref_rec):
        # Don't fail — pair on the shorter — but flag in stderr so
        # the user knows there's a length mismatch on the receptor.
        print(f"[{pred_cif_path.parent.name}] receptor length mismatch: "
              f"pred={len(pred_rec)} ref={len(ref_rec)}; "
              f"pairing on first {rec_n}", file=sys.stderr)

    # Kabsch superposition on receptor.
    P = np.array(pred_rec[:rec_n])
    Q = np.array(ref_rec[:rec_n])
    P_c, P_m = P - P.mean(axis=0), P.mean(axis=0)
    Q_c, Q_m = Q - Q.mean(axis=0), Q.mean(axis=0)
    H = P_c.T @ Q_c
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T

    if eff_n < 3:
        return None, (f"too_few_effector_residues:{eff_n} "
                      f"(pred={len(pred_eff)}, ref={len(ref_eff)})")
    if len(pred_eff) != len(ref_eff):
        print(f"[{pred_cif_path.parent.name}] effector length mismatch: "
              f"pred={len(pred_eff)} ref={len(ref_eff)}; "
              f"pairing on first {eff_n}", file=sys.stderr)

    P_eff = np.array(pred_eff[:eff_n])
    Q_eff = np.array(ref_eff[:eff_n])
    P_eff_aligned = (P_eff - P_m) @ R.T + Q_m
    diff = P_eff_aligned - Q_eff
    return float(np.sqrt((diff * diff).sum(axis=1).mean())), None


def _read_confidence(cif_path: Path) -> Tuple[Optional[float], Optional[str]]:
    """Read ipTM from the AF3 summary_confidences.json that sits next to a CIF."""
    for candidate in [
        cif_path.parent / "summary_confidences.json",
        cif_path.parent / f"{cif_path.stem.replace('_model', '')}_summary_confidences.json",
    ]:
        if candidate.is_file():
            try:
                with open(candidate) as f:
                    data = json.load(f)
                iptm = data.get("iptm")
                if iptm is not None:
                    return float(iptm), None
            except (json.JSONDecodeError, OSError) as e:
                return None, f"confidence_load_failed:{type(e).__name__}"
    return None, "confidence_json_not_found"


# ═════════════════════════════════════════════════════════════════════════
# Main
# ═════════════════════════════════════════════════════════════════════════

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seq-name", required=True)
    ap.add_argument("--af3-output-dir", required=True, type=Path)
    ap.add_argument("--ground-truth", required=True, type=Path)
    ap.add_argument("--receptor-chain", default="A")
    ap.add_argument("--effector-chain", default="B")
    ap.add_argument("--ra-eff-threshold", type=float, default=5.0)
    ap.add_argument("--output-csv", required=True, type=Path)
    args = ap.parse_args()

    failures: List[str] = []
    fieldnames = [
        "seq_name",
        "af3_nomsa_best_ra_eff",
        "af3_nomsa_mean_ra_eff",
        "af3_nomsa_best_iptm",
        "af3_nomsa_mean_iptm",
        "af3_nomsa_n_correct_interface",
        "af3_nomsa_total_predictions",
        "af3_nomsa_failures",
    ]

    # Find every per-sample CIF under the AF3 output dir.
    #
    # AF3's actual output convention (verified empirically on the
    # 2026-04-29 run):
    #   output/<name>/seed-<N>_sample-<M>/model.cif    ← 15 per-sample
    #                                                    predictions
    #                                                    (3 seeds × 5
    #                                                    diffusion
    #                                                    samples by
    #                                                    default)
    #   output/<name>/<name>_model.cif                  ← duplicate of
    #                                                    AF3's highest-
    #                                                    ranked sample
    #                                                    (excluded from
    #                                                    aggregation to
    #                                                    avoid double-
    #                                                    counting and
    #                                                    biasing the
    #                                                    median toward
    #                                                    AF3's own pick)
    #
    # The earlier `*_model.cif` glob matched ONLY the top-level
    # <name>_model.cif (because the 15 per-sample files are bare
    # "model.cif", no name prefix), so the aggregation was running on
    # 1 CIF instead of 15.  The current pattern explicitly targets the
    # per-sample directory layout.
    if not args.af3_output_dir.is_dir():
        failures.append(f"af3_output_dir_missing:{args.af3_output_dir}")
        _write_single(args.output_csv, fieldnames,
                      {"seq_name": args.seq_name,
                       "af3_nomsa_failures": ",".join(failures)})
        return 0

    # Match seed-*/sample-*/model.cif rather than the broader **/model.cif
    # so we don't accidentally pull in any future top-level files AF3
    # may add.  The seed-*_sample-* directory naming has been stable
    # since the AF3 v2 release.
    cifs = sorted(args.af3_output_dir.rglob("seed-*_sample-*/model.cif"))
    if not cifs:
        # Fallback for older / non-standard AF3 output layouts: any
        # model.cif under the tree, excluding the top-level
        # <name>_model.cif (which is a renamed duplicate of one of the
        # per-sample CIFs).
        cifs = [p for p in sorted(args.af3_output_dir.rglob("*.cif"))
                if p.name == "model.cif" or
                   (p.name.endswith(".cif") and not p.name.endswith("_model.cif"))]
    if not cifs:
        failures.append("no_cif_files_found")
        _write_single(args.output_csv, fieldnames,
                      {"seq_name": args.seq_name,
                       "af3_nomsa_failures": ",".join(failures)})
        return 0

    print(f"[{args.seq_name}] found {len(cifs)} AF3 CIF predictions")
    ra_effs: List[float] = []
    iptms: List[float] = []
    per_prediction_errors = 0

    for cif in cifs:
        ra_eff, rmsd_err = _receptor_aligned_effector_rmsd(
            cif, args.ground_truth,
            args.receptor_chain, args.effector_chain,
        )
        iptm, iptm_err = _read_confidence(cif)
        if ra_eff is None:
            per_prediction_errors += 1
            # Echo the actual error so .command.err carries diagnostic
            # info.  Without this, the per_prediction_rmsd_errors:N
            # tag in the CSV is the only signal that anything went
            # wrong, with no clue as to WHY.  rmsd_err comes from
            # _receptor_aligned_effector_rmsd; iptm_err is independent
            # (a prediction can fail RMSD parsing but still yield an
            # iptm, so log them separately).
            err_msg = rmsd_err or "(no error message returned)"
            print(f"[{args.seq_name}] RMSD failed for {cif.parent.name}: "
                  f"{err_msg}", file=sys.stderr)
            continue
        ra_effs.append(ra_eff)
        if iptm is not None:
            iptms.append(iptm)
        elif iptm_err:
            print(f"[{args.seq_name}] iPTM read failed for "
                  f"{cif.parent.name}: {iptm_err}", file=sys.stderr)

    if per_prediction_errors:
        failures.append(f"per_prediction_rmsd_errors:{per_prediction_errors}")

    if not ra_effs:
        failures.append("no_successful_rmsd")
        _write_single(args.output_csv, fieldnames,
                      {"seq_name": args.seq_name,
                       "af3_nomsa_failures": ",".join(failures)})
        return 0

    row = {
        "seq_name": args.seq_name,
        "af3_nomsa_best_ra_eff":  f"{min(ra_effs):.3f}",
        "af3_nomsa_mean_ra_eff":  f"{sum(ra_effs) / len(ra_effs):.3f}",
        "af3_nomsa_best_iptm":    f"{max(iptms):.4f}" if iptms else "",
        "af3_nomsa_mean_iptm":    f"{sum(iptms) / len(iptms):.4f}" if iptms else "",
        "af3_nomsa_n_correct_interface":
            str(sum(1 for r in ra_effs if r < args.ra_eff_threshold)),
        "af3_nomsa_total_predictions": str(len(ra_effs)),
        "af3_nomsa_failures": ",".join(failures),
    }
    _write_single(args.output_csv, fieldnames, row)
    print(f"[{args.seq_name}] best ra_eff = {row['af3_nomsa_best_ra_eff']} Å, "
          f"n_correct = {row['af3_nomsa_n_correct_interface']}/{row['af3_nomsa_total_predictions']}")
    return 0


def _write_single(path: Path, fieldnames: List[str], row: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        w.writerow(row)


if __name__ == "__main__":
    raise SystemExit(main())