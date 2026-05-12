#!/usr/bin/env python3
"""
extract_survivor_manifest.py
----------------------------
P0-31 · Given the P0-29 extended CSV and the set of per-sequence
workdirs, emit a one-row-per-survivor manifest that NextFlow can
splitCsv over to fan out the three orthogonal streams.

"Survivor" here is any row with:
  - a non-empty rep_canonical_pdb
  - a workdir under the runs/ glob containing plan.json

For each survivor, the manifest provides:
  seq_name, canonical_pdb_abs, ground_truth_abs,
  effector_template_cif_abs, receptor_seq, effector_seq

Sequences are extracted from chain A and chain B of the canonical_pdb
itself — the final designed construct, same one Boltz predicted on.
This guarantees AF3 sees exactly the sequence that will go to wet-lab.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Dict, List, Optional

# _AA3TO1 dict removed in Phase 4 free-function consolidation —
# extract_chain_seq below now delegates to the canonical
# boltz2_negative_steering.get_chain_sequence (which uses the
# THREE_TO_ONE map in contig_utils.py).


def _extract_chain_seq(pdb_path: Path, chain_id: str) -> Optional[str]:
    """Read Cα sequence for a chain from a PDB file.

    Thin wrapper around the canonical
    ``boltz2_negative_steering.get_chain_sequence`` (Phase 4
    free-function consolidation).  Preserves the previous behaviour of
    returning ``None`` on read failure or on a missing chain (callers
    depend on the Optional return; the canonical helper raises on
    missing chain).
    """
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from boltz2_negative_steering import get_chain_sequence  # noqa: E402
    try:
        seq = get_chain_sequence(pdb_path, chain_id)
    except (OSError, ValueError):
        return None
    return seq or None


def _find_workdir(seq_name: str, workdirs_glob: str) -> Optional[Path]:
    from glob import glob
    for wd in glob(workdirs_glob):
        if Path(wd).name == seq_name:
            return Path(wd)
    return None


def _resolve_plan_json(workdir: Path) -> Optional[Path]:
    for candidate in [
        workdir / "plan.json",
        workdir / "cycle_0" / "plan.json",
    ]:
        if candidate.is_file():
            return candidate
    return None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-csv", required=True, type=Path,
                    help="cross_sequence_summary_with_interface_metrics.csv")
    ap.add_argument("--workdirs-glob", required=True, type=str,
                    help="Glob for per-sequence workdirs containing plan.json")
    ap.add_argument("--receptor-chain", default="A")
    ap.add_argument("--effector-chain", default="B")
    ap.add_argument("--output-manifest", required=True, type=Path)
    args = ap.parse_args()

    with open(args.input_csv, newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    out_rows: List[Dict[str, str]] = []
    skipped: Dict[str, int] = {}

    for row in rows:
        seq_name = row.get("mpnn_sequence", "")
        if not seq_name:
            skipped["no_mpnn_sequence"] = skipped.get("no_mpnn_sequence", 0) + 1
            continue

        canonical_pdb = row.get("rep_canonical_pdb", "")
        if not canonical_pdb or not Path(canonical_pdb).is_file():
            skipped["canonical_pdb_missing"] = skipped.get("canonical_pdb_missing", 0) + 1
            continue

        workdir = _find_workdir(seq_name, args.workdirs_glob)
        if workdir is None:
            skipped["workdir_missing"] = skipped.get("workdir_missing", 0) + 1
            continue

        plan_json = _resolve_plan_json(workdir)
        if plan_json is None:
            skipped["plan_json_missing"] = skipped.get("plan_json_missing", 0) + 1
            continue

        try:
            with open(plan_json) as f:
                plan = json.load(f)
        except (json.JSONDecodeError, OSError):
            skipped["plan_json_unreadable"] = skipped.get("plan_json_unreadable", 0) + 1
            continue

        ground_truth = plan.get("ground_truth")
        effector_template = plan.get("effector_template_cif")
        if not ground_truth or not Path(ground_truth).is_file():
            skipped["ground_truth_missing"] = skipped.get("ground_truth_missing", 0) + 1
            continue
        if not effector_template or not Path(effector_template).is_file():
            skipped["effector_template_missing"] = skipped.get("effector_template_missing", 0) + 1
            continue

        receptor_seq = _extract_chain_seq(Path(canonical_pdb), args.receptor_chain)
        effector_seq = _extract_chain_seq(Path(canonical_pdb), args.effector_chain)
        if not receptor_seq or not effector_seq:
            skipped["seq_extraction_failed"] = skipped.get("seq_extraction_failed", 0) + 1
            continue

        out_rows.append({
            "seq_name": seq_name,
            "canonical_pdb_abs":     str(Path(canonical_pdb).resolve()),
            "ground_truth_abs":      str(Path(ground_truth).resolve()),
            "effector_template_cif_abs": str(Path(effector_template).resolve()),
            "receptor_seq": receptor_seq,
            "effector_seq": effector_seq,
        })

    args.output_manifest.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output_manifest, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["seq_name", "canonical_pdb_abs", "ground_truth_abs",
                        "effector_template_cif_abs", "receptor_seq", "effector_seq"],
        )
        writer.writeheader()
        for r in out_rows:
            writer.writerow(r)

    print(f"Manifest: {len(out_rows)} survivors written to {args.output_manifest}")
    if skipped:
        print("Skipped:")
        for reason, n in sorted(skipped.items()):
            print(f"  {reason}: {n}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())