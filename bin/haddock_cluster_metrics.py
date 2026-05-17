#!/usr/bin/env python3
"""
haddock_cluster_metrics.py
--------------------------
Compute per-cluster metrics for a completed HADDOCK3 docking run and
emit ``cluster_metrics.json`` (keyed by cluster_id) for HaddockRun.

Reads:
- ``<workdir>/haddock_report.json`` — cluster list + best_cluster*.pdb names.
- ``<workdir>/restraints_summary.json`` — contact pairs + active residues
  used to generate the restraints (so the metric set knows the
  denominators).
- ``<workdir>/best_cluster<id>.pdb`` — per-cluster best model PDBs.

For each cluster:
- BSA + interface H-bonds via ``run_biophysical_metrics.py`` subprocess.
- Sc + ΔΔG via ``run_rosetta_metrics.py`` subprocess (ΔΔG discarded; the
  HADDOCK stage doesn't use it).  Per A150 — adding ~3–12 min of
  FastRelax overhead is acceptable.
- COM distance between receptor and effector via structure_metrics.
- Pair contact fraction: count of contact pairs with CA-CA distance
  within ``pair_distance[0] + pair_distance[2]`` (the restraint upper bound).
- AIR satisfaction count: how many of the receptor active residues have
  ANY CA within the AIR upper bound of an effector active residue (or
  any effector CA when effector_active is empty).
- Clash counts in/outside the design region via structure_metrics.

Bridge distances are NOT computed here — they require the contig string,
which is intentionally not a HaddockRun field (per A139 / Q152(i)).  The
Nextflow wrapper that has access to both this output and the contig
fills bridge distances onto HaddockCluster post-construction.

Outputs:
    cluster_metrics.json  - {"<cluster_id>": {metric: value, ...}, ...}
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np

_BIN_DIR = os.path.dirname(os.path.abspath(__file__))
if _BIN_DIR not in sys.path:
    sys.path.insert(0, _BIN_DIR)
import structure_metrics  # noqa: E402


# Heavy-atom clash cutoff and AIR upper-bound default mirror
# haddock3_prepare.py constants; one shared source could be a follow-up
# refactor but small enough to repeat here for now.
CLASH_CUTOFF_A = 2.0
# When AIRs are two-sided (eff active residues present) we use the
# 3 +- 3 +- 5 distance triple => up to 8 A.  When one-sided (entire
# effector chain), we use 5 +- 5 +- 5 => up to 10 A.
AIR_TWO_SIDED_UPPER = 8.0
AIR_ONE_SIDED_UPPER = 10.0


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workdir", required=True, type=Path,
                    help="HADDOCK3_DOCK output dir containing haddock_report.json")
    ap.add_argument("--receptor-chain", default="A",
                    help="Receptor chain inside HADDOCK output (default A)")
    ap.add_argument("--effector-chain", default="B",
                    help="Effector chain inside HADDOCK output (default B)")
    ap.add_argument("--biophysical-script",
                    default=str(Path(_BIN_DIR) / "run_biophysical_metrics.py"),
                    help="Path to run_biophysical_metrics.py")
    ap.add_argument("--rosetta-script",
                    default=str(Path(_BIN_DIR) / "run_rosetta_metrics.py"),
                    help="Path to run_rosetta_metrics.py")
    ap.add_argument("--fastrelax-xml",
                    default=str(Path(_BIN_DIR) / "fastrelax_for_ia.xml"),
                    help="Path to fastrelax_for_ia.xml for Rosetta InterfaceAnalyzer")
    ap.add_argument("--skip-sc", action="store_true",
                    help="Skip Rosetta Sc computation (debugging only)")
    return ap.parse_args()


# ── Subprocess invocations ──────────────────────────────────────────


def _read_csv_single_row(path: Path) -> Dict[str, str]:
    """Read a single-row CSV emitted by run_*_metrics.py.  Returns the
    row as a dict; empty dict if file missing or empty.
    """
    if not path.is_file():
        return {}
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            return dict(row)
    return {}


def _run_biophys(
    biophys_script: Path, pdb: Path, receptor_chain: str, effector_chain: str,
    cluster_id: int,
) -> Tuple[Optional[float], Optional[int], List[str]]:
    """Call run_biophysical_metrics.py; return (bsa, hbonds, failures)."""
    with tempfile.TemporaryDirectory(prefix=f"biophys_c{cluster_id}_") as tmp:
        out_csv = Path(tmp) / "out.csv"
        cmd = [
            "python3", str(biophys_script),
            "--seq-name", f"cluster_{cluster_id}",
            "--canonical-pdb", str(pdb),
            "--ground-truth", str(pdb),  # unused per the script's comment
            "--receptor-chain", receptor_chain,
            "--effector-chain", effector_chain,
            "--output-csv", str(out_csv),
        ]
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            print(f"[cluster {cluster_id}] biophysical subprocess failed: "
                  f"{e.stderr}", file=sys.stderr)
            return None, None, ["biophysical_subprocess_failed"]
        row = _read_csv_single_row(out_csv)
    bsa = _safe_float(row.get("bsa", ""))
    hbonds = _safe_int(row.get("interface_hbonds", ""))
    failures = []
    raw_fail = row.get("biophysical_failures", "")
    if raw_fail:
        failures.extend(raw_fail.split(","))
    return bsa, hbonds, failures


def _run_rosetta_sc(
    rosetta_script: Path, fastrelax_xml: Path, pdb: Path,
    receptor_chain: str, effector_chain: str, cluster_id: int,
) -> Tuple[Optional[float], List[str]]:
    """Call run_rosetta_metrics.py; return (sc, failures).  ΔΔG is discarded."""
    with tempfile.TemporaryDirectory(prefix=f"rosetta_c{cluster_id}_") as tmp:
        out_csv = Path(tmp) / "out.csv"
        cmd = [
            "python3", str(rosetta_script),
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
    sc = _safe_float(row.get("sc", ""))
    failures = []
    raw_fail = row.get("rosetta_failures", "")
    if raw_fail:
        failures.extend(raw_fail.split(","))
    return sc, failures


def _safe_float(s) -> Optional[float]:
    if s in ("", None):
        return None
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def _safe_int(s) -> Optional[int]:
    if s in ("", None):
        return None
    try:
        return int(s)
    except (TypeError, ValueError):
        return None


# ── Distance-based metrics ──────────────────────────────────────────


def _ca_by_resnum(pdb: Path, chain: str) -> Dict[int, np.ndarray]:
    """Return {resnum: CA_xyz_array} for the chain."""
    out: Dict[int, np.ndarray] = {}
    for line in open(pdb):
        if not line.startswith("ATOM"):
            continue
        if line[21] != chain:
            continue
        if line[12:16].strip() != "CA":
            continue
        try:
            rn = int(line[22:26].strip())
        except ValueError:
            continue
        if rn in out:
            continue
        x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
        out[rn] = np.array([x, y, z], dtype=float)
    return out


def _pair_contact_fraction(
    pdb: Path,
    contact_pairs: List[dict],
    pair_upper_bound: float,
    receptor_chain: str, effector_chain: str,
) -> Tuple[float, int]:
    """Fraction of user contact pairs whose CA-CA distance on this PDB is
    within ``pair_upper_bound``.  Returns ``(fraction, n_total)``.

    Per A154's "vacuous truth" convention: when there are no pairs the
    fraction reports 1.0 so the auto-pick lexicographic sort still works.
    """
    n_total = len(contact_pairs)
    if n_total == 0:
        return 1.0, 0
    rec_cas = _ca_by_resnum(pdb, receptor_chain)
    eff_cas = _ca_by_resnum(pdb, effector_chain)
    n_satisfied = 0
    for pair in contact_pairs:
        rn = int(pair["rec_resnum"])
        en = int(pair["eff_resnum"])
        a = rec_cas.get(rn)
        b = eff_cas.get(en)
        if a is None or b is None:
            continue
        if float(np.linalg.norm(a - b)) <= pair_upper_bound:
            n_satisfied += 1
    return n_satisfied / n_total, n_total


def _air_satisfaction(
    pdb: Path,
    receptor_active: List[int],
    effector_active: List[int],
    receptor_chain: str, effector_chain: str,
) -> Tuple[int, int]:
    """Per-AIR satisfaction count.  Each receptor active residue
    contributes one AIR; it's "satisfied" iff any CA within the AIR upper
    bound is found on the effector side (active list if non-empty,
    otherwise any effector residue).

    Returns ``(satisfied, total)``.
    """
    total = len(receptor_active)
    if total == 0:
        return 0, 0
    rec_cas = _ca_by_resnum(pdb, receptor_chain)
    eff_cas = _ca_by_resnum(pdb, effector_chain)
    if effector_active:
        eff_subset = [eff_cas[rn] for rn in effector_active if rn in eff_cas]
        upper = AIR_TWO_SIDED_UPPER
    else:
        eff_subset = list(eff_cas.values())
        upper = AIR_ONE_SIDED_UPPER
    if not eff_subset:
        return 0, total
    eff_arr = np.array(eff_subset)
    satisfied = 0
    for rn in receptor_active:
        ca = rec_cas.get(rn)
        if ca is None:
            continue
        diffs = eff_arr - ca
        min_d = float(np.min(np.sqrt(np.einsum("ij,ij->i", diffs, diffs))))
        if min_d <= upper:
            satisfied += 1
    return satisfied, total


def _design_region(restraints: dict) -> Set[int]:
    """Receptor residues counted as "in design region" for clash bookkeeping.

    Per Session 7 post-commit-3 amendment: prefer the contig-derived
    design region (computed by haddock3_prepare.py from the contig
    string + receptor PDB).  Falls back to pair receptor halves +
    receptor_active_residues when the contig wasn't provided (e.g. a
    test invocation that hand-rolled a restraints_summary.json).
    """
    contig_dr = restraints.get("contig_design_region")
    if contig_dr:
        return set(int(r) for r in contig_dr)
    rec_from_pairs = {int(p["rec_resnum"]) for p in restraints.get("contact_pairs", [])}
    rec_active = set(int(r) for r in restraints.get("receptor_active_residues", []))
    return rec_from_pairs | rec_active


# ── Main ────────────────────────────────────────────────────────────


def main() -> int:
    args = parse_args()
    workdir = args.workdir.resolve()

    report_path = workdir / "haddock_report.json"
    restraints_path = workdir / "restraints_summary.json"
    if not report_path.is_file():
        print(f"ERROR: {report_path} not found", file=sys.stderr)
        return 1
    if not restraints_path.is_file():
        print(f"ERROR: {restraints_path} not found", file=sys.stderr)
        return 1

    report = json.loads(report_path.read_text())
    restraints = json.loads(restraints_path.read_text())

    cluster_models = report.get("cluster_models", {})
    if not cluster_models:
        print("WARNING: no qualifying clusters in haddock_report.json", file=sys.stderr)
        (workdir / "cluster_metrics.json").write_text("{}\n")
        return 0

    pair_upper = restraints["pair_distance"][0] + restraints["pair_distance"][2]
    design_region = _design_region(restraints)

    out: Dict[str, dict] = {}
    for cid_str, info in cluster_models.items():
        cid = int(cid_str)
        pdb = workdir / info["filename"]
        if not pdb.is_file():
            print(f"[cluster {cid}] PDB {pdb} not found; skipping", file=sys.stderr)
            continue

        # BSA + hbonds via biophysical subprocess.
        bsa, hbonds, biophys_fail = _run_biophys(
            Path(args.biophysical_script), pdb,
            args.receptor_chain, args.effector_chain, cid,
        )

        # Sc via Rosetta subprocess.
        if args.skip_sc:
            sc, ros_fail = None, []
        else:
            sc, ros_fail = _run_rosetta_sc(
                Path(args.rosetta_script), Path(args.fastrelax_xml),
                pdb, args.receptor_chain, args.effector_chain, cid,
            )

        # COM distance.
        com_dist = structure_metrics.chain_com_distance(
            pdb, args.receptor_chain, args.effector_chain,
        )

        # Pair satisfaction.
        pair_frac, pair_total = _pair_contact_fraction(
            pdb, restraints.get("contact_pairs", []),
            pair_upper, args.receptor_chain, args.effector_chain,
        )

        # AIR satisfaction.
        air_sat, air_total = _air_satisfaction(
            pdb,
            restraints.get("receptor_active_residues", []),
            restraints.get("effector_active_residues", []),
            args.receptor_chain, args.effector_chain,
        )

        # Clash counts.
        clashes_in, clashes_out = structure_metrics.clash_count(
            pdb, args.receptor_chain, args.effector_chain,
            design_region, cutoff=CLASH_CUTOFF_A,
        )

        out[str(cid)] = {
            "bsa": bsa if bsa is not None else 0.0,
            "interface_hbonds": hbonds,
            "sc": sc,
            "com_distance": com_dist if com_dist is not None else 0.0,
            "air_satisfaction_count": air_sat,
            "air_total_count": air_total,
            "pair_contact_fraction": pair_frac,
            "pair_total_count": pair_total,
            "clashes_in_design_region": clashes_in,
            "clashes_outside_design_region": clashes_out,
            "failures": biophys_fail + ros_fail,
        }
        print(f"[cluster {cid}] BSA={out[str(cid)]['bsa']:.2f}  "
              f"COM={out[str(cid)]['com_distance']:.2f}  "
              f"AIR={air_sat}/{air_total}  pairs={pair_frac*100:.0f}%  "
              f"clashes={clashes_in}/{clashes_out}  "
              f"Sc={'' if sc is None else f'{sc:.3f}'}")

    (workdir / "cluster_metrics.json").write_text(json.dumps(out, indent=2))
    print(f"Wrote {len(out)} cluster metric record(s) -> "
          f"{workdir / 'cluster_metrics.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
