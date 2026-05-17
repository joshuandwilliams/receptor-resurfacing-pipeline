#!/usr/bin/env python3
"""
haddock_cluster_sc.py
---------------------
Compute Rosetta shape complementarity (Sc) for every qualifying HADDOCK
cluster's best model and emit ``cluster_sc.json`` (keyed by cluster_id).

Sibling to ``haddock_cluster_metrics.py``: that script runs in
``boltz2_container`` (needs numpy + freesasa + MDAnalysis); this one
runs in ``rosetta_container`` (needs the Rosetta binaries) because
neither container has both stacks.  The two output JSONs are merged
downstream by ``HaddockRun.from_workdir``.

Reads:
- ``<workdir>/haddock_report.json`` — cluster list + best_cluster*.pdb names.
- ``<workdir>/best_cluster<id>.pdb`` — per-cluster best model PDBs.

For each cluster:
- Calls ``run_rosetta_metrics.py`` as a subprocess for Sc (the ΔΔG
  output is discarded; HADDOCK's geometric-placement objective doesn't
  use it).

Outputs:
    cluster_sc.json  - {"<cluster_id>": {"sc": float|null, "failures": [...]}}
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

_BIN_DIR = os.path.dirname(os.path.abspath(__file__))


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workdir", required=True, type=Path,
                    help="HADDOCK3_DOCK output dir containing haddock_report.json")
    ap.add_argument("--receptor-chain", default="A",
                    help="Receptor chain inside HADDOCK output (default A)")
    ap.add_argument("--effector-chain", default="B",
                    help="Effector chain inside HADDOCK output (default B)")
    ap.add_argument("--rosetta-script",
                    default=str(Path(_BIN_DIR) / "run_rosetta_metrics.py"),
                    help="Path to run_rosetta_metrics.py")
    ap.add_argument("--fastrelax-xml",
                    default=str(Path(_BIN_DIR) / "fastrelax_for_ia.xml"),
                    help="Path to fastrelax_for_ia.xml for Rosetta InterfaceAnalyzer")
    return ap.parse_args()


def _read_csv_single_row(path: Path) -> Dict[str, str]:
    if not path.is_file():
        return {}
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            return dict(row)
    return {}


def _run_rosetta_sc(
    rosetta_script: Path, fastrelax_xml: Path, pdb: Path,
    receptor_chain: str, effector_chain: str, cluster_id: int,
) -> Tuple[Optional[float], List[str]]:
    with tempfile.TemporaryDirectory(prefix=f"rosetta_c{cluster_id}_") as tmp:
        out_csv = Path(tmp) / "out.csv"
        # sys.executable (not 'python3') so the subprocess uses the same
        # interpreter the parent is running under — guarantees consistent
        # site-packages across container environments.
        cmd = [
            sys.executable, str(rosetta_script),
            "--seq-name", f"cluster_{cluster_id}",
            "--canonical-pdb", str(pdb),
            "--receptor-chain", receptor_chain,
            "--effector-chain", effector_chain,
            "--fast-relax-xml", str(fastrelax_xml),
            "--output-csv", str(out_csv),
        ]
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            print(f"[cluster {cluster_id}] rosetta subprocess failed: "
                  f"{e.stderr}", file=sys.stderr)
            return None, ["rosetta_subprocess_failed"]
        row = _read_csv_single_row(out_csv)
    raw_sc = row.get("sc", "")
    sc = None
    if raw_sc not in ("", None):
        try:
            sc = float(raw_sc)
        except (TypeError, ValueError):
            sc = None
    failures: List[str] = []
    raw_fail = row.get("rosetta_failures", "")
    if raw_fail:
        failures.extend(raw_fail.split(","))
    return sc, failures


def main() -> int:
    args = parse_args()
    workdir = args.workdir.resolve()
    report_path = workdir / "haddock_report.json"
    if not report_path.is_file():
        print(f"ERROR: {report_path} not found", file=sys.stderr)
        return 1
    report = json.loads(report_path.read_text())
    cluster_models = report.get("cluster_models", {})
    if not cluster_models:
        print("WARNING: no qualifying clusters in haddock_report.json",
              file=sys.stderr)
        (workdir / "cluster_sc.json").write_text("{}\n")
        return 0

    out: Dict[str, dict] = {}
    for cid_str, info in cluster_models.items():
        cid = int(cid_str)
        pdb = workdir / info["filename"]
        if not pdb.is_file():
            print(f"[cluster {cid}] PDB {pdb} not found; skipping",
                  file=sys.stderr)
            continue
        sc, failures = _run_rosetta_sc(
            Path(args.rosetta_script), Path(args.fastrelax_xml),
            pdb, args.receptor_chain, args.effector_chain, cid,
        )
        out[str(cid)] = {"sc": sc, "failures": failures}
        print(f"[cluster {cid}] Sc={'' if sc is None else f'{sc:.3f}'}")

    (workdir / "cluster_sc.json").write_text(json.dumps(out, indent=2))
    print(f"Wrote {len(out)} cluster Sc record(s) -> "
          f"{workdir / 'cluster_sc.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
