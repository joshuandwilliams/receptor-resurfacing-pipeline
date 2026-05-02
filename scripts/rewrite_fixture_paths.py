#!/usr/bin/env python3
"""
rewrite_fixture_paths.py
------------------------
Rewrite absolute HPC paths inside a staged orthogonal_metrics fixture
tree to point at the local (Mac/repo) staging location.

Background
==========
The discovery run on HPC writes per-sequence workdirs under
`tests/full_test_run/results/negative_steering/runs/<seq>/`.  Files
inside those workdirs (notably `cycle_0/plan.json` and the cohort
`cross_sequence_summary.csv`) carry absolute paths to:

  - the per-sequence workdir itself (`/hpc-home/.../runs/<seq>/...`)
  - the parent RFDiffusion design PDB used as `ground_truth`
  - the per-sequence `effector_template_cif`
  - per-cold-start-seed PDB paths
  - the singularity container image
  - Nextflow scratch `work/<hash>/` paths (read-only refs)

When the fixture is staged into the repo at
`tests/orthogonal_metrics/data/negsteer_run/runs/<seq>/`, those
absolute paths point at HPC locations that don't exist on Mac.  This
script rewrites them.

Rewrite rules
=============
For each per-sequence workdir under `--staged-runs-dir`:

  1. plan.json fields rewritten:
       - ground_truth: pointed at the staged parent design PDB
         under <design_pdbs_dir>/design_<N>.pdb (read from the
         `_RFDIFF_DESIGN_STEM` filename token if present, otherwise
         inferred from the current ground_truth basename).
       - effector_template_cif: pointed at
         <staged-runs-dir>/<seq>/cycle_0/effector_template.cif
       - cold_start_seeds[*].pdb_path: pointed at
         <staged-runs-dir>/<seq>/cycle_0/initial_prediction[_s<N>].pdb
       - boltz_container, receptor_fasta_override, effector_fasta_override,
         designs[*].dir, designs[*].yaml: blanked (set to "").  These
         refer to scratch work dirs or HPC-only resources that the
         orthogonal_metrics test never reads.

  2. Inside cross_sequence_summary.csv (passed via --cross-csv):
       - representative_canonical_pdb: rewrite the workdir prefix
         /hpc-home/.../runs/<seq>/  →  <staged-runs-dir>/<seq>/

The rewriter is idempotent: paths already pointing at the local
staging are left alone.

Usage
=====
    python3 scripts/rewrite_fixture_paths.py \\
        --staged-runs-dir tests/orthogonal_metrics/data/negsteer_run/runs \\
        --cross-csv       tests/orthogonal_metrics/data/negsteer_run/cross_sequence_summary.csv \\
        --design-pdbs-dir tests/negative_steering/data/design_pdbs \\
        --hpc-prefix      /hpc-home/jowillia/receptor_design/receptor-resurfacing-pipeline/tests/full_test_run/results

Stdlib only — no third-party imports.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional


# ── Helpers ──────────────────────────────────────────────────────────


def _rewrite_plan_json(
    plan_path: Path,
    seq_name: str,
    staged_runs_dir: Path,
    design_pdbs_dir: Path,
    hpc_prefix: str,
) -> Dict[str, str]:
    """Rewrite paths inside one cycle_0/plan.json.  Returns a summary
    dict mapping field-name → before/after for logging."""
    with plan_path.open() as f:
        plan = json.load(f)

    summary: Dict[str, str] = {}
    seq_runs_local = staged_runs_dir / seq_name
    cycle_local = seq_runs_local / "cycle_0"

    def _rec(field: str, before: str, after: str) -> None:
        if before != after:
            summary[field] = f"{before!s:.80} → {after!s:.80}"

    # 1. ground_truth — the parent RFDiffusion design PDB.  Infer
    # design stem from the current ground_truth basename (e.g.
    # design_15.pdb → design_15) and point at the staged copy.
    gt_old = plan.get("ground_truth", "") or ""
    if gt_old:
        gt_basename = Path(gt_old).name  # e.g. "design_15.pdb"
        gt_new = str((design_pdbs_dir / gt_basename).resolve())
        plan["ground_truth"] = gt_new
        _rec("ground_truth", gt_old, gt_new)

    # 2. effector_template_cif — lives inside cycle_0/.
    et_old = plan.get("effector_template_cif", "") or ""
    if et_old:
        et_new = str((cycle_local / "effector_template.cif").resolve())
        plan["effector_template_cif"] = et_new
        _rec("effector_template_cif", et_old, et_new)

    # 3. cold_start_seeds[*].pdb_path
    css = plan.get("cold_start_seeds") or []
    for i, seed in enumerate(css):
        if not isinstance(seed, dict):
            continue
        old = seed.get("pdb_path", "") or ""
        if not old:
            continue
        # Determine local target.  Seed 0 → initial_prediction.pdb;
        # seed N>0 → initial_prediction_s<N>.pdb (per
        # boltz2_negative_steering.py:1992).
        offset = seed.get("seed_offset", i)
        if offset == 0:
            local = cycle_local / "initial_prediction.pdb"
        else:
            local = cycle_local / f"initial_prediction_s{offset}.pdb"
        seed["pdb_path"] = str(local.resolve())
        _rec(f"cold_start_seeds[{i}].pdb_path", old, str(local))

    # 4. Blank out fields that point at HPC scratch / containers.
    # These are not read by the orthogonal_metrics test cascade.
    for key in (
        "boltz_container",
        "receptor_fasta_override",
        "effector_fasta_override",
    ):
        old = plan.get(key, "")
        if old and isinstance(old, str) and old.startswith("/"):
            plan[key] = ""
            _rec(key, old, "")

    designs = plan.get("designs") or []
    for i, d in enumerate(designs):
        if not isinstance(d, dict):
            continue
        for key in ("dir", "yaml"):
            old = d.get(key, "")
            if old and isinstance(old, str) and old.startswith("/"):
                d[key] = ""
                _rec(f"designs[{i}].{key}", old, "")

    plan_path.write_text(json.dumps(plan, indent=2))
    return summary


def _rewrite_cross_csv(
    csv_path: Path,
    staged_runs_dir: Path,
    hpc_prefix: str,
) -> int:
    """Rewrite representative_canonical_pdb in cross_sequence_summary.csv.
    Returns the number of rows whose path was rewritten."""
    with csv_path.open(newline="") as f:
        rows = list(csv.DictReader(f))
        fieldnames = list(rows[0].keys()) if rows else []

    rewritten = 0
    runs_hpc_root = (
        f"{hpc_prefix.rstrip('/')}/negative_steering/runs/"
    )
    for r in rows:
        old = r.get("representative_canonical_pdb", "") or ""
        if not old:
            continue
        if old.startswith(runs_hpc_root):
            tail = old[len(runs_hpc_root):]
            new = str((staged_runs_dir / tail).resolve())
            r["representative_canonical_pdb"] = new
            rewritten += 1

    if not fieldnames:
        return 0

    with csv_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)
    return rewritten


# ── CLI ──────────────────────────────────────────────────────────────


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--staged-runs-dir", required=True, type=Path,
        help="Local dir holding the staged per-sequence workdirs "
             "(e.g. tests/orthogonal_metrics/data/negsteer_run/runs)",
    )
    ap.add_argument(
        "--cross-csv", required=True, type=Path,
        help="Path to the staged cross_sequence_summary.csv",
    )
    ap.add_argument(
        "--design-pdbs-dir", required=True, type=Path,
        help="Local dir holding parent RFDiffusion design PDBs "
             "(e.g. tests/negative_steering/data/design_pdbs).  "
             "Used as the rewrite target for plan.ground_truth.",
    )
    ap.add_argument(
        "--hpc-prefix", required=True, type=str,
        help="HPC absolute path prefix to strip, e.g. "
             "'/hpc-home/jowillia/receptor_design/receptor-resurfacing-pipeline/tests/full_test_run/results'",
    )
    args = ap.parse_args()

    staged = args.staged_runs_dir.resolve()
    design_pdbs = args.design_pdbs_dir.resolve()
    cross_csv = args.cross_csv.resolve()

    if not staged.is_dir():
        print(f"ERROR: --staged-runs-dir not found: {staged}", file=sys.stderr)
        return 2
    if not cross_csv.is_file():
        print(f"ERROR: --cross-csv not found: {cross_csv}", file=sys.stderr)
        return 2
    if not design_pdbs.is_dir():
        print(f"ERROR: --design-pdbs-dir not found: {design_pdbs}", file=sys.stderr)
        return 2

    # Per-sequence plan.json rewrite
    n_seqs = 0
    for seq_dir in sorted(staged.iterdir()):
        if not seq_dir.is_dir():
            continue
        plan_path = seq_dir / "cycle_0" / "plan.json"
        if not plan_path.is_file():
            print(f"  {seq_dir.name}: SKIP (no cycle_0/plan.json)")
            continue
        summary = _rewrite_plan_json(
            plan_path,
            seq_dir.name,
            staged,
            design_pdbs,
            args.hpc_prefix,
        )
        n_seqs += 1
        print(f"  {seq_dir.name}: rewrote {len(summary)} plan.json fields")
        for field, change in summary.items():
            print(f"    {field}: {change}")

    # Cross-summary rewrite
    n_rewritten = _rewrite_cross_csv(cross_csv, staged, args.hpc_prefix)
    print(f"  cross_sequence_summary.csv: rewrote "
          f"representative_canonical_pdb on {n_rewritten} row(s)")

    print(f"Done.  {n_seqs} per-sequence plan.json file(s) updated.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
