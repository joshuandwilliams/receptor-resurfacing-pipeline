#!/usr/bin/env python3
"""
compute_interface_metrics.py
----------------------------
Add the DockQ family of interface metrics to cross_sequence_summary.csv.

For each row in the input CSV, appends four columns:

  - irmsd           Interface-restricted effector Cα RMSD (via DockQ).
  - fnat            Fraction of native contacts preserved (via DockQ).
  - dockq           CAPRI-standard 0-1 score (via DockQ API).
  - dockq_failures  Comma-separated failure reasons, blank on success.

The model is the `rep_canonical_pdb` cached Boltz prediction.
The native is the parent RFDiffusion design PDB (receptor+effector),
the same ground-truth structure used by the existing ra_eff metric.
Its path is read from the per-sequence workdir's `plan.json` under the
`ground_truth` key.

The other interface metrics this script used to emit (15 Å iPSAE,
intact_core, weighted_jaccard, interface_plddt) are now computed per
prediction by the negative-steering engine.  They reach this CSV as
`rep_*_median` columns through passing_summary.csv and
cross_sequence_summary.py.  The helpers that computed them are still
defined below, but nothing calls them.

Inputs
------
--input-csv        cross_sequence_summary.csv from NEGSTEER_CROSS_SEQUENCE.
                   Must have a `rep_canonical_pdb` column.
--output-csv       Where to write the extended CSV.
--workdirs-glob    Glob pattern matching per-sequence workdirs, e.g.
                   '<outdir>/runs/design_*_seq_*'.  Used to locate
                   plan.json for each row.
--receptor-chain   Usually 'A'.  Must match the Boltz prediction chain.
--effector-chain   Usually 'B'.

--plddt-threshold (default 50), --intact-threshold (default 5.0) and
--weighted-jaccard-pair-cutoff (default 8.0) are parsed but unused.
They stay declared because the Nextflow process passes all three.

Usage
-----
  compute_interface_metrics.py \\
      --input-csv  cross_sequence_summary.csv \\
      --output-csv cross_sequence_summary_with_interface_metrics.csv \\
      --workdirs-glob 'runs/design_*_seq_*' \\
      --receptor-chain A --effector-chain B

Notes
-----
- The 10 Å iPSAE values are not recomputed here.  ipsae_min arrives as
  the `rep_ipsae_min_median` column from passing_summary.  That keeps
  the existing ranker behaviour unchanged.
- The engine's 15 Å iPSAE cutoff follows Dunbrack 2025, which
  recommends 10-15 Å.  The 10 Å value was already in the source data.
- A row whose DockQ call fails gets blank irmsd/fnat/dockq and a reason
  in `dockq_failures`.  Failures are flagged, never silently dropped.
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

# Put bin/ on sys.path so sibling modules are importable.  Vestigial
# now — the ipSAE machinery this once reused moved to the
# negative-steering engine, so nothing here imports a sibling.
sys.path.insert(0, str(Path(__file__).resolve().parent))


# DockQ-only column set.  The other interface metrics (15 Å iPSAE,
# intact_core, weighted_jaccard, interface_plddt) are computed per
# prediction by the negative-steering engine and propagated through
# passing_summary.csv → cross_sequence_summary.csv as `rep_*_median`.
# This script only adds the DockQ family, which genuinely needs a
# post-hoc external tool invocation.
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
    canonical_pdb.  Mirrors the layout the negative-steering engine
    writes, so path resolution stays consistent.

    Boltz names these files after the model.  They usually live under
    one of:
        <pdb_dir>/<prefix>_<name>.<ext>
        <pdb_dir>/../<prefix>/<prefix>_<name>.<ext>
    Failing both, the search recurses from two levels above <pdb_dir>.
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
    workdir's plan.json.  `seq_name` is the row's `mpnn_sequence` value
    (e.g. "design_3_seq_1"), matched against the workdir basename.
    """
    from glob import glob
    for wd in glob(workdirs_glob):
        wd_path = Path(wd)
        if wd_path.name != seq_name:
            continue
        # plan.json sits either at the workdir root or under cycle_0.
        # Try both.
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
    Compute iRMSD, fnat and DockQ via the DockQ Python API.

    Returns (irmsd, fnat, dockq, error).  The first three are None when
    the call failed.  `error` is a short string on failure, None on
    success.
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


# ═════════════════════════════════════════════════════════════════════════
# Weighted Jaccard contact overlap
# ═════════════════════════════════════════════════════════════════════════

# Gaussian weight parameters:
#
#     w(d) = exp( -(d - mu)^2 / (2 * sigma^2) )
#          = exp( -(d - 4)^2 / 2.25 )
#
# So mu = 4.0 Å (the empirical mode of close-contact Cβ–Cβ distances) and
# 2 * sigma^2 = 2.25, giving sigma^2 = 1.125 and sigma ≈ 1.061 Å.
#
# These are deliberately HARD-CODED rather than CLI-tunable so that
# weighted_jaccard values stay comparable across cohorts.  Re-tuning
# would silently invalidate the planned cross-target AUROC comparison of
# the overlap metrics.  The column would keep its name while changing
# meaning.  A cohort needing different parameters should get a NEW
# column (e.g. weighted_jaccard_v2) rather than a re-tune of this one.
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
    to Cα for glycine and any other residue lacking a Cβ atom (e.g. the
    rare missing-side-chain case in a Cα-only RFDiffusion model).

    Positional indexing — NOT seqid.num — is what makes the (i, j) keys
    comparable between the prediction PDB and the ground-truth PDB even
    when their numbering schemes differ.  find_contact_residues_heavy()
    in structure_metrics.py and the residue-set true_jaccard upstream
    use the same convention.
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
    over all pairs within `pair_cutoff` Å Cβ–Cβ.  Pairs beyond the cutoff
    are dropped.  At σ ≈ 1.06 Å they contribute weight ≈ 0 anyway, and
    keeping them only inflates the union denominator.
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
        # Empty chain → empty pair map.  The caller treats this as a
        # legitimate "no contacts" rather than an error.
        return {}

    # Vectorised distance matrix: |rec| × |eff|.  Both chains run
    # ~75-150 residues in this pipeline, so the matrix is ~10k entries
    # at most and the all-vs-all Cβ comparison is fast.
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

    The other interface metrics (15 Å iPSAE, intact_core,
    weighted_jaccard, interface_plddt) are computed per prediction by
    the negative-steering engine.  They reach this CSV as `rep_*_median`
    columns through cross_sequence_summary.py and are NOT recomputed
    here.
    """
    out = dict(row)
    for col in NEW_COLUMNS:
        out[col] = ""
    failures: List[str] = []

    canonical_pdb_str = row.get("rep_canonical_pdb", "")
    if not canonical_pdb_str:
        out["dockq_failures"] = "no_canonical_pdb"
        return out
    model_pdb = Path(canonical_pdb_str)
    if not model_pdb.is_file():
        out["dockq_failures"] = f"canonical_pdb_missing:{model_pdb.name}"
        return out

    # Ground truth comes from plan.json under the per-sequence workdir.
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
    # The next three flags are accepted for backward compatibility only.
    # The metrics they used to govern (intact_core, weighted_jaccard) are
    # computed per prediction by the negative-steering engine, and this
    # script only emits the DockQ family.  The args stay declared so
    # callers that still pass them on the CLI don't break.
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

    # Append the new columns, preserving any that already exist.  A
    # re-run on an already-extended CSV overwrites the same columns
    # rather than duplicating them.
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
