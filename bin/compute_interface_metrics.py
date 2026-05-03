#!/usr/bin/env python3
"""
compute_interface_metrics.py
----------------------------
P0-29 · Add interface-restricted metrics to cross_sequence_summary.csv.
P0-38 · Add a smooth distance-weighted Jaccard contact overlap.

For each row in the input CSV, adds eight columns:

  - irmsd          Interface-restricted effector Cα RMSD (via DockQ).
  - fnat           Fraction of native contacts preserved (via DockQ).
  - dockq          CAPRI-standard 0-1 score (via DockQ API).
  - ipsae_ab_15    ipSAE chain A→B at PAE cutoff 15 Å (Dunbrack 2025
                   recommends 10-15 Å; 10 Å already in source data).
  - ipsae_ba_15    ipSAE chain B→A at PAE cutoff 15 Å.
  - ipsae_min_15   min(ipsae_ab_15, ipsae_ba_15).
  - intact_core    1/0 flag — whole-chain receptor-aligned RMSD after
                   trimming residues below a pLDDT threshold (core-only
                   replacement for the current whole-chain intact
                   filter).
  - weighted_jaccard
                   Smooth distance-weighted contact overlap between the
                   model and the native (P0-38).  Generalises the
                   existing residue-set true_jaccard by:
                     (a) operating on (receptor i, effector j) PAIRS
                         rather than receptor residues alone — strictly
                         more discriminating;
                     (b) weighting each pair by its Cβ–Cβ distance
                         under exp(−(d − 4)² / 2.25) instead of a hard
                         5 Å cutoff — smooths the boundary so that a
                         single-Å predictor displacement no longer
                         flips contacts in/out and amplifies noise.
                   See _compute_weighted_jaccard() for the full
                   formula.  Co-exists with true_jaccard which
                   continues to flow through unchanged from the
                   per-prediction inner loop, so Task 8 AUROC can
                   compare the two head-to-head before any composite-
                   score migration.

The model is the `representative_canonical_pdb` cached Boltz prediction.
The native is the parent RFDiffusion design PDB (receptor+effector),
the same ground-truth structure used by the existing ra_eff metric.
Its path is read from the per-sequence workdir's `plan.json` under the
`ground_truth` key.

Inputs
------
--input-csv        cross_sequence_summary.csv from NEGSTEER_CROSS_SEQUENCE.
                   Must have a `representative_canonical_pdb` column.
--output-csv       Where to write the extended CSV.
--workdirs-glob    Glob pattern matching per-sequence workdirs, e.g.
                   '<outdir>/runs/design_*_seq_*'.  Used to locate
                   plan.json for each row.
--plddt-threshold  Per-residue pLDDT cutoff for the intact_core filter.
                   Residues with pLDDT below this value are dropped
                   before the receptor-aligned RMSD is computed.
                   Default 50.
--intact-threshold Å cutoff for intact_core.  Default 5.0 (matches the
                   pre-existing whole-chain `intact` filter).
--receptor-chain   Usually 'A'.  Must match the Boltz prediction chain.
--effector-chain   Usually 'B'.
--weighted-jaccard-pair-cutoff
                   Hard distance ceiling (Å) for residue pairs included
                   in the weighted_jaccard sum.  Pairs above this are
                   dropped to keep the all-vs-all loop O(N_iface²)
                   rather than O(N_total²).  Default 8.0 — at d=8 the
                   Gaussian weight is ≈ 8e-4, so the cutoff is
                   numerically inert (it discards weight that would
                   round to zero anyway).

Usage
-----
  compute_interface_metrics.py \\
      --input-csv  cross_sequence_summary.csv \\
      --output-csv cross_sequence_summary_with_interface_metrics.csv \\
      --workdirs-glob 'runs/design_*_seq_*' \\
      --receptor-chain A --effector-chain B

Notes
-----
- The existing 10 Å iPSAE values (ipsae_ab, ipsae_ba, ipsae_min) are
  left as-is and NOT recomputed; they come through as
  `representative_ipsae_*` columns from the per-row passing_summary.
  This keeps the existing ranker behaviour unchanged.
- Rows where DockQ, ipSAE, or pLDDT fail have the corresponding
  columns set to NaN and a comma-separated reason written to a new
  `interface_metrics_failures` column.  Failures are flagged, not
  silently dropped (todo_list3 P0-29 acceptance criterion).
- DockQ's Python API finds the best chain mapping between model and
  native automatically.  For the pipeline's (A, B) two-chain
  complexes, this is trivial.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

# Reuse the existing ipSAE machinery rather than reimplementing it.
# compute_metrics.py is shipped alongside this script in bin/.
sys.path.insert(0, str(Path(__file__).resolve().parent))


# DockQ-only column set.  All other interface metrics (15Å iPSAE,
# intact_core, weighted_jaccard, interface_plddt) are now computed
# per-prediction inside compute_metrics.py and propagated through
# passing_summary.csv → cross_sequence_summary.csv as
# `representative_*_median` — see boltz2_iterate_steering.py's
# METRIC_KEYS_RAW for the authoritative list.  This script only
# adds the DockQ family, which genuinely needs a post-hoc external
# tool invocation.
NEW_COLUMNS = [
    "irmsd",
    "fnat",
    "dockq",
    "dockq_failures",
]


# ═════════════════════════════════════════════════════════════════════════
# Path resolution
# ═════════════════════════════════════════════════════════════════════════

def _find_sidecar(pdb_path: Path, prefix: str, ext: str) -> Optional[Path]:
    """
    Locate Boltz sidecar files (confidence JSON, PAE npz) for a given
    canonical_pdb.  Mirrors the glob pattern used by compute_metrics.py's
    existing logic so path resolution is consistent.

    Boltz emits these files with a base name matching the model name.
    They typically live under one of:
        <pdb_dir>/<prefix>_<name>.<ext>
        <pdb_dir>/../<prefix>/<prefix>_<name>.<ext>
    """
    pdb_dir = pdb_path.parent
    name = pdb_path.stem
    candidates = [
        pdb_dir / f"{prefix}_{name}.{ext}",
        pdb_dir.parent / prefix / f"{prefix}_{name}.{ext}",
    ]
    for c in candidates:
        if c.is_file():
            return c
    # Fallback: recursive search from pred_dir two levels up.
    pred_dir = pdb_dir.parent.parent
    if pred_dir.is_dir():
        matches = list(pred_dir.rglob(f"{prefix}_{name}.{ext}"))
        if matches:
            return matches[0]
    return None


def _find_ground_truth(seq_name: str, workdirs_glob: str) -> Optional[Path]:
    """
    Locate the ground-truth PDB for a row by reading its per-sequence
    workdir's plan.json.  `seq_name` is the row's
    representative_canonical_pdb workdir suffix (e.g. "design_3_seq_1").
    """
    from glob import glob
    for wd in glob(workdirs_glob):
        wd_path = Path(wd)
        if wd_path.name != seq_name:
            continue
        # plan.json lives at the cycle_0 level typically; try a couple
        # of common layouts.
        for plan_candidate in [
            wd_path / "plan.json",
            wd_path / "cycle_0" / "plan.json",
        ]:
            if plan_candidate.is_file():
                try:
                    with open(plan_candidate) as f:
                        plan = json.load(f)
                    gt = plan.get("ground_truth")
                    if gt:
                        return Path(gt)
                except (json.JSONDecodeError, OSError):
                    pass
    return None


# ═════════════════════════════════════════════════════════════════════════
# Per-row metric computation
# ═════════════════════════════════════════════════════════════════════════

def _run_dockq(
    model_pdb: Path,
    native_pdb: Path,
    receptor_chain: str,
    effector_chain: str,
) -> Tuple[Optional[float], Optional[float], Optional[float], Optional[str]]:
    """
    Compute iRMSD, fnat, DockQ via the DockQ Python API.

    Returns (irmsd, fnat, dockq, error) — the first three may be None
    if the call failed; error is a short string on failure or None
    on success.
    """
    try:
        from DockQ.DockQ import load_PDB, run_on_all_native_interfaces
    except ImportError as e:
        return None, None, None, f"dockq_import_failed:{e}"

    try:
        model = load_PDB(str(model_pdb))
        native = load_PDB(str(native_pdb))
    except Exception as e:  # noqa: BLE001 — DockQ raises various errors
        return None, None, None, f"dockq_load_failed:{type(e).__name__}"

    try:
        # run_on_all_native_interfaces auto-picks the best chain mapping.
        # Result is a dict keyed by chain-pair tuple.
        results, _total = run_on_all_native_interfaces(model, native)
    except Exception as e:  # noqa: BLE001
        return None, None, None, f"dockq_run_failed:{type(e).__name__}"

    # For two-chain A/B complexes there's exactly one interface.
    # Prefer the (receptor, effector) pair if present; otherwise take
    # whatever DockQ found (it may have flipped the mapping).
    interface_key = None
    for key in results:
        if set(key) == {receptor_chain, effector_chain}:
            interface_key = key
            break
    if interface_key is None and results:
        interface_key = next(iter(results))

    if interface_key is None:
        return None, None, None, "dockq_no_interface"

    d = results[interface_key]
    return (
        float(d.get("iRMSD", float("nan"))),
        float(d.get("fnat", float("nan"))),
        float(d.get("DockQ", float("nan"))),
        None,
    )


def _load_pae(
    canonical_pdb: Path,
) -> Tuple[Optional[np.ndarray], Optional[Tuple[int, ...]], Optional[str]]:
    """
    Load the PAE matrix and per-chain lengths for a cached Boltz
    prediction.

    Returns (pae_matrix, chain_lengths, error).
    """
    pae_npz = _find_sidecar(canonical_pdb, prefix="pae", ext="npz")
    if pae_npz is None:
        return None, None, "pae_npz_not_found"
    conf_json = _find_sidecar(canonical_pdb, prefix="confidence", ext="json")
    if conf_json is None:
        return None, None, "confidence_json_not_found"

    try:
        data = np.load(pae_npz)
        # Boltz stores under 'pae' key (same convention as AlphaFold).
        pae_matrix = data["pae"] if "pae" in data else data[data.files[0]]
    except (OSError, KeyError) as e:
        return None, None, f"pae_load_failed:{type(e).__name__}"

    try:
        with open(conf_json) as f:
            conf = json.load(f)
    except (json.JSONDecodeError, OSError) as e:
        return None, None, f"confidence_load_failed:{type(e).__name__}"

    # Boltz confidence json carries chain lengths either under
    # 'chain_lengths' directly, or can be derived from 'token_chain_ids'
    # counts.  Prefer the explicit field if present.
    chain_lengths = conf.get("chain_lengths")
    if chain_lengths is None:
        token_ids = conf.get("token_chain_ids")
        if token_ids is not None:
            # Count consecutive occurrences per chain.
            from collections import Counter
            counts = Counter(token_ids)
            # Preserve order of first appearance.
            seen: List[str] = []
            for cid in token_ids:
                if cid not in seen:
                    seen.append(cid)
            chain_lengths = [counts[c] for c in seen]
    if chain_lengths is None:
        return pae_matrix, None, "chain_lengths_unresolved"

    return pae_matrix, tuple(chain_lengths), None


def _compute_intact_core(
    model_pdb: Path,
    native_pdb: Path,
    plddt_threshold: float,
    intact_threshold: float,
    receptor_chain: str,
    effector_chain: str,
) -> Tuple[Optional[int], Optional[float], Optional[str]]:
    """
    Compute the core-only intact filter: align receptor cores (trimmed
    by per-residue pLDDT), compute whole-complex Cα RMSD on the trimmed
    set, emit 1 if RMSD < intact_threshold else 0.

    Returns (intact_flag, rmsd_angstroms, error).
    """
    try:
        import gemmi
    except ImportError as e:
        return None, None, f"gemmi_import_failed:{e}"

    # Per-residue pLDDT is stored in the B-factor field of Boltz PDBs.
    # This matches the convention in compute_metrics.plddt_from_pdb.
    try:
        model = gemmi.read_structure(str(model_pdb))
        native = gemmi.read_structure(str(native_pdb))
    except Exception as e:  # noqa: BLE001
        return None, None, f"structure_load_failed:{type(e).__name__}"

    def chain_ca_records(st, chain_id):
        """Yield (resnum, coord, bfactor) for every CA in the chain."""
        for m in st:
            for ch in m:
                if ch.name != chain_id:
                    continue
                for res in ch:
                    for atom in res:
                        if atom.name == "CA":
                            yield (res.seqid.num, np.array([
                                atom.pos.x, atom.pos.y, atom.pos.z,
                            ]), atom.b_iso)
                            break

    try:
        model_recs = {num: (xyz, bf) for (num, xyz, bf)
                      in chain_ca_records(model, receptor_chain)}
        native_recs = {num: (xyz, bf) for (num, xyz, bf)
                       in chain_ca_records(native, receptor_chain)}
    except Exception as e:  # noqa: BLE001
        return None, None, f"chain_parse_failed:{type(e).__name__}"

    # Align on residues present in both, with pLDDT >= threshold in the
    # model (the native is the RFDiffusion reference — no pLDDT).
    common_resnums = sorted(set(model_recs) & set(native_recs))
    kept = [r for r in common_resnums
            if model_recs[r][1] is not None
            and model_recs[r][1] >= plddt_threshold]
    if len(kept) < 3:
        return None, None, f"too_few_core_residues:{len(kept)}"

    # Kabsch alignment on receptor core.
    M = np.array([model_recs[r][0] for r in kept])
    N = np.array([native_recs[r][0] for r in kept])
    M_c = M - M.mean(axis=0)
    N_c = N - N.mean(axis=0)
    H = M_c.T @ N_c
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T

    # Apply transform to the whole model (receptor + effector).
    model_translate = N.mean(axis=0) - M.mean(axis=0) @ R.T
    # Gather all CAs from both chains in the model for RMSD.
    try:
        all_model = list(chain_ca_records(model, receptor_chain)) + \
                    list(chain_ca_records(model, effector_chain))
        all_native = list(chain_ca_records(native, receptor_chain)) + \
                     list(chain_ca_records(native, effector_chain))
    except Exception as e:  # noqa: BLE001
        return None, None, f"whole_chain_parse_failed:{type(e).__name__}"

    # Match residue numbers across both chains for RMSD.
    model_full = {(chain_id, num): xyz
                  for chain_id, recs in [(receptor_chain, chain_ca_records(model, receptor_chain)),
                                          (effector_chain, chain_ca_records(model, effector_chain))]
                  for (num, xyz, _bf) in recs}
    native_full = {(chain_id, num): xyz
                   for chain_id, recs in [(receptor_chain, chain_ca_records(native, receptor_chain)),
                                           (effector_chain, chain_ca_records(native, effector_chain))]
                   for (num, xyz, _bf) in recs}

    shared = sorted(set(model_full) & set(native_full))
    if len(shared) < 3:
        return None, None, f"too_few_shared_residues:{len(shared)}"

    M_all = np.array([model_full[k] for k in shared])
    N_all = np.array([native_full[k] for k in shared])
    M_aligned = (M_all - M.mean(axis=0)) @ R.T + N.mean(axis=0)
    diff = M_aligned - N_all
    rmsd = float(np.sqrt((diff * diff).sum(axis=1).mean()))
    return (1 if rmsd < intact_threshold else 0), round(rmsd, 3), None


# ═════════════════════════════════════════════════════════════════════════
# P0-38 · Weighted Jaccard contact overlap
# ═════════════════════════════════════════════════════════════════════════

# Gaussian weight parameters from the v5/v6 todo_list spec:
#
#     w(d) = exp( -(d - mu)^2 / (2 * sigma^2) )
#          = exp( -(d - 4)^2 / 2.25 )
#
# So mu = 4.0 Å (the empirical mode of close-contact Cβ–Cβ distances) and
# 2 * sigma^2 = 2.25  =>  sigma^2 = 1.125  =>  sigma ≈ 1.061 Å.
#
# These are deliberately HARD-CODED rather than CLI-tunable so that
# weighted_jaccard values are directly comparable across cohorts.  Re-
# tuning would invalidate cross-target AUROC analysis (Task 8) without
# anyone noticing — the column would still exist with the same name but
# different semantics.  If a future cohort needs different parameters,
# emit a NEW column (e.g. weighted_jaccard_v2) rather than re-tuning
# this one.
_WJ_MU = 4.0
_WJ_TWO_SIGMA_SQ = 2.25


def _gaussian_contact_weight(distance_angstroms: float) -> float:
    """Per-pair Gaussian weight: exp(-(d - 4)^2 / 2.25)."""
    delta = distance_angstroms - _WJ_MU
    return float(np.exp(-(delta * delta) / _WJ_TWO_SIGMA_SQ))


def _collect_chain_cb_positions(
    structure,
    chain_id: str,
) -> List[Tuple[int, np.ndarray, str]]:
    """
    Return a list of (positional_index_0b, cb_xyz, residue_name) for
    every residue in the named chain, walking residues in PDB file
    order.  Cβ is taken from the side chain when present, falling back
    to Cα for glycine (and any other residue lacking a Cβ atom — e.g.
    the rare missing-side-chain case in a Cα-only RFDiffusion model).

    Positional indexing — NOT seqid.num — is what makes the (i, j) keys
    comparable between the prediction PDB and the ground-truth PDB even
    when their numbering schemes differ.  This is the same convention
    used by find_contact_residues_heavy() in compute_metrics.py and by
    the residue-set true_jaccard upstream.
    """
    out: List[Tuple[int, np.ndarray, str]] = []
    pos_idx = 0
    for model in structure:
        for chain in model:
            if chain.name != chain_id:
                continue
            for res in chain:
                cb_atom = None
                ca_atom = None
                for atom in res:
                    if atom.name == "CB":
                        cb_atom = atom
                    elif atom.name == "CA":
                        ca_atom = atom
                # Cα fallback handles glycine and any missing-Cβ case
                # (e.g. Cα-only ground-truth PDBs from RFDiffusion).
                anchor = cb_atom if cb_atom is not None else ca_atom
                if anchor is None:
                    pos_idx += 1
                    continue
                xyz = np.array(
                    [anchor.pos.x, anchor.pos.y, anchor.pos.z],
                    dtype=np.float64,
                )
                out.append((pos_idx, xyz, res.name))
                pos_idx += 1
            # Only the first matching chain — gemmi may return
            # duplicates if a structure has multiple models.
            break
        # First model only.
        break
    return out


def _build_weighted_pair_map(
    pdb_path: Path,
    receptor_chain: str,
    effector_chain: str,
    pair_cutoff: float,
) -> Dict[Tuple[int, int], float]:
    """
    For one structure, return {(receptor_pos_0b, effector_pos_0b): weight}
    over all pairs within `pair_cutoff` Å Cβ–Cβ.  Pairs with d > cutoff
    are dropped — at the spec's σ they contribute weight ≈ 0 anyway and
    keeping them only inflates the union denominator with negligible
    summands.
    """
    try:
        import gemmi
    except ImportError as e:
        raise RuntimeError(f"gemmi_import_failed:{e}") from e

    try:
        st = gemmi.read_structure(str(pdb_path))
    except Exception as e:  # noqa: BLE001
        raise RuntimeError(
            f"structure_load_failed:{type(e).__name__}"
        ) from e

    rec = _collect_chain_cb_positions(st, receptor_chain)
    eff = _collect_chain_cb_positions(st, effector_chain)
    if not rec or not eff:
        # Empty chain → empty pair map.  Caller treats this as a
        # legitimate "no contacts" rather than an error.
        return {}

    # Vectorised distance matrix: |rec| × |eff|.  Both chains are
    # ~75-150 residues in this pipeline; the matrix is ~10k entries
    # max, so the all-vs-all Cβ comparison is comfortably fast.
    rec_xyz = np.array([r[1] for r in rec])  # (R, 3)
    eff_xyz = np.array([e[1] for e in eff])  # (E, 3)
    diffs = rec_xyz[:, None, :] - eff_xyz[None, :, :]  # (R, E, 3)
    dists = np.sqrt((diffs * diffs).sum(axis=2))       # (R, E)

    cutoff_sq = float(pair_cutoff)
    pair_map: Dict[Tuple[int, int], float] = {}
    # Vectorise the threshold mask, then iterate only over the small
    # number of pairs that pass.  np.where keeps memory bounded.
    keep_i, keep_j = np.where(dists <= cutoff_sq)
    for ii, jj in zip(keep_i.tolist(), keep_j.tolist()):
        d = float(dists[ii, jj])
        rec_pos = rec[ii][0]
        eff_pos = eff[jj][0]
        pair_map[(rec_pos, eff_pos)] = _gaussian_contact_weight(d)
    return pair_map


def _compute_weighted_jaccard(
    model_pdb: Path,
    native_pdb: Path,
    receptor_chain: str,
    effector_chain: str,
    pair_cutoff: float = 8.0,
) -> Tuple[Optional[float], Optional[str]]:
    """
    Smooth distance-weighted contact overlap between model and native.

    Per the v5/v6 todo_list P0-38 spec, this is the generalised
    (weighted) Jaccard index over the union of (receptor, effector)
    Cβ–Cβ residue PAIRS:

        weighted_jaccard
            = sum_pairs min(w_model(d), w_native(d))
            / sum_pairs max(w_model(d), w_native(d))

    where w(d) = exp(-(d - 4)^2 / 2.25), and pairs absent from one
    structure contribute 0 there.

    This co-exists with the existing residue-set true_jaccard (which
    flows through cross_sequence_summary unchanged from the per-
    prediction inner loop).  Both columns are emitted side-by-side so
    Task 8 AUROC can determine which one better discriminates correct-
    vs-wrong interface classification before any composite-score
    migration.

    Returns
    -------
    (weighted_jaccard, error)
        weighted_jaccard is in [0, 1] when both pair maps are non-
        empty; NaN when both are empty (no interface to compare —
        matches the convention of jaccard() in boltz2_negative_
        steering.py); 0.0 when exactly one is empty (legitimate "no
        overlap").  error is a short string on failure or None.
    """
    try:
        model_pairs = _build_weighted_pair_map(
            model_pdb, receptor_chain, effector_chain, pair_cutoff,
        )
    except RuntimeError as e:
        return None, f"wj_model:{e}"
    try:
        native_pairs = _build_weighted_pair_map(
            native_pdb, receptor_chain, effector_chain, pair_cutoff,
        )
    except RuntimeError as e:
        return None, f"wj_native:{e}"

    if not model_pairs and not native_pairs:
        return float("nan"), None
    if not model_pairs or not native_pairs:
        return 0.0, None

    union_keys = set(model_pairs) | set(native_pairs)
    num = 0.0
    den = 0.0
    for k in union_keys:
        wm = model_pairs.get(k, 0.0)
        wn = native_pairs.get(k, 0.0)
        num += min(wm, wn)
        den += max(wm, wn)
    if den == 0.0:
        # All pairs in the union have zero weight in both structures —
        # only possible if every kept pair sits at d ≫ 4 in both.
        # Treat as no-meaningful-overlap rather than a divide-by-zero.
        return 0.0, None
    return num / den, None


# ═════════════════════════════════════════════════════════════════════════
# Main
# ═════════════════════════════════════════════════════════════════════════

def process_row(
    row: Dict[str, str],
    workdirs_glob: str,
    receptor_chain: str,
    effector_chain: str,
) -> Dict[str, str]:
    """Populate the DockQ family for one row.

    Other interface metrics (15Å iPSAE, intact_core, weighted_jaccard,
    interface_plddt) are now computed per-prediction by
    compute_metrics.py and reach this CSV as `representative_*_median`
    columns through cross_sequence_summary.py — they are NOT recomputed
    here.
    """
    out = dict(row)
    for col in NEW_COLUMNS:
        out[col] = ""
    failures: List[str] = []

    canonical_pdb_str = row.get("representative_canonical_pdb", "")
    if not canonical_pdb_str:
        out["dockq_failures"] = "no_canonical_pdb"
        return out
    model_pdb = Path(canonical_pdb_str)
    if not model_pdb.is_file():
        out["dockq_failures"] = f"canonical_pdb_missing:{model_pdb.name}"
        return out

    # Ground truth: read from plan.json under the per-sequence workdir.
    seq_name = row.get("mpnn_sequence", "")
    native_pdb = _find_ground_truth(seq_name, workdirs_glob)
    if native_pdb is None or not native_pdb.is_file():
        out["dockq_failures"] = f"ground_truth_not_found_for:{seq_name}"
        return out

    # --- DockQ: irmsd, fnat, dockq --------------------------------------
    irmsd, fnat, dockq, dq_err = _run_dockq(
        model_pdb, native_pdb, receptor_chain, effector_chain,
    )
    if dq_err:
        failures.append(dq_err)
    else:
        out["irmsd"] = f"{irmsd:.3f}" if not math.isnan(irmsd) else ""
        out["fnat"] = f"{fnat:.4f}" if not math.isnan(fnat) else ""
        out["dockq"] = f"{dockq:.4f}" if not math.isnan(dockq) else ""

    out["dockq_failures"] = ",".join(failures) if failures else ""
    return out


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Add P0-29 interface-restricted metrics to a "
                    "cross_sequence_summary.csv.",
    )
    ap.add_argument("--input-csv", required=True, type=Path)
    ap.add_argument("--output-csv", required=True, type=Path)
    ap.add_argument("--workdirs-glob", required=True, type=str,
                    help="Glob matching per-sequence workdirs (contain plan.json).")
    # The following flags are accepted for backward compatibility only.
    # The metrics they used to govern (intact_core, weighted_jaccard) are
    # now computed per-prediction in compute_metrics.py — this script
    # only emits the DockQ family.  Keeping the args declared so callers
    # passing them on the CLI don't break.
    ap.add_argument("--plddt-threshold", type=float, default=50.0,
                    help="(unused, accepted for compat) intact_core is "
                         "now computed per-prediction by compute_metrics.py.")
    ap.add_argument("--intact-threshold", type=float, default=5.0,
                    help="(unused, accepted for compat) intact_core is "
                         "now computed per-prediction by compute_metrics.py.")
    ap.add_argument("--receptor-chain", type=str, default="A")
    ap.add_argument("--effector-chain", type=str, default="B")
    ap.add_argument(
        "--weighted-jaccard-pair-cutoff", type=float, default=8.0,
        help="(unused, accepted for compat) weighted_jaccard is now "
             "computed per-prediction by compute_metrics.py.",
    )
    args = ap.parse_args()

    with open(args.input_csv, newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            print("ERROR: input CSV has no header.", file=sys.stderr)
            return 1
        in_fieldnames = list(reader.fieldnames)
        rows = list(reader)

    # Append new columns, preserving any that already exist (idempotent
    # re-runs: if re-run on an already-extended CSV, we just overwrite
    # the same columns rather than duplicating them).
    out_fieldnames = list(in_fieldnames)
    for col in NEW_COLUMNS:
        if col not in out_fieldnames:
            out_fieldnames.append(col)

    print(f"Processing {len(rows)} rows from {args.input_csv.name}")
    out_rows: List[Dict[str, str]] = []
    n_ok = 0
    n_flagged = 0
    for i, row in enumerate(rows, 1):
        enriched = process_row(
            row, args.workdirs_glob,
            args.receptor_chain, args.effector_chain,
        )
        if enriched.get("dockq_failures"):
            n_flagged += 1
        else:
            n_ok += 1
        out_rows.append(enriched)
        seq = row.get("mpnn_sequence", f"row_{i}")
        dockq = enriched.get("dockq", "—")
        irmsd = enriched.get("irmsd", "—")
        flags = enriched.get("dockq_failures", "")
        status = f"flags={flags}" if flags else "ok"
        print(f"  [{i:3d}/{len(rows)}] {seq:35s}  dockq={dockq:>6s}  "
              f"irmsd={irmsd:>7s}  {status}")

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=out_fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in out_rows:
            writer.writerow(row)

    print(f"\nWrote {args.output_csv}")
    print(f"  {n_ok} rows fully populated, {n_flagged} flagged with failures.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())