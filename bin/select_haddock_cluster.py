#!/usr/bin/env python3
"""
select_haddock_cluster.py
-------------------------
Pick the HADDOCK cluster that downstream stages consume.

Per Session 7 grill-me (notes/design_audit.md Q134, Q141, Q154):
- When ``--chosen-cluster-id N`` is set, select cluster N (validator
  upstream guarantees it's in the qualifying set).
- Otherwise auto-pick using ``HaddockCluster.auto_pick_key``:
  lexicographic ``(-pair_contact_fraction, -bsa)``.  Pair contact fraction
  is the primary key because user-specified pair contacts encode
  biological intent that RFDiffusion cannot recover; BSA breaks ties
  because RFDiffusion can improve BSA but can't redirect contacts.

Inputs (workdir):
- ``haddock_report.json``      - emitted by collect_haddock3_dock.py
- ``cluster_metrics.json``     - emitted by haddock_cluster_metrics.py
- ``restraints_summary.json``  - emitted by haddock3_prepare.py
- ``best_cluster{N}.pdb``      - per-cluster best models

Outputs (current dir):
- ``selected_complex.pdb``     - the chosen cluster's best model PDB
- ``selected_cluster_id.txt``  - the chosen cluster_id as a line
- ``cluster_metrics_table.csv``- all clusters with metrics, sorted by
                                 auto-pick key for easy inspection
"""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import sys
from pathlib import Path

_BIN_DIR = os.path.dirname(os.path.abspath(__file__))
if _BIN_DIR not in sys.path:
    sys.path.insert(0, _BIN_DIR)
from haddock_run import HaddockRun  # noqa: E402


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workdir", required=True, type=Path,
                    help="HADDOCK3 workdir (contains haddock_report.json + "
                         "cluster_metrics.json + restraints_summary.json + "
                         "best_cluster*.pdb)")
    ap.add_argument("--receptor-pdb", required=True, type=Path,
                    help="Input receptor monomer PDB (used by HaddockRun "
                         "construction; not modified)")
    ap.add_argument("--effector-pdb", required=True, type=Path,
                    help="Input effector monomer PDB")
    ap.add_argument("--contact-pairs", default="",
                    help="Pair string (must match what was passed to "
                         "haddock3_prepare.py).")
    ap.add_argument("--receptor-active-residues", default="",
                    help="Active residue string (receptor side).")
    ap.add_argument("--effector-active-residues", default="",
                    help="Active residue string (effector side).")
    ap.add_argument("--pair-distance", default="2,2,4",
                    help='Pair distance triple "target,lo_dev,hi_dev".')
    ap.add_argument("--chosen-cluster-id", default=None,
                    help="Cluster ID to use; auto-pick when blank/null.")
    return ap.parse_args()


def _parse_chosen(value):
    """Nextflow passes 'null' as a string when params.haddock_chosen_cluster
    is unset.  Normalise to None.
    """
    if value is None:
        return None
    s = str(value).strip()
    if s in ("", "null", "None"):
        return None
    return int(s)


def _write_cluster_metrics_table(run, out_csv: Path) -> None:
    """Write cluster_metrics_table.csv sorted by auto-pick key.

    Top row = cluster the auto-pick would select.  The chosen cluster
    is highlighted by the ``is_selected`` column.
    """
    clusters_sorted = sorted(run.qualifying_clusters, key=lambda c: c.auto_pick_key)
    selected_id = run.selected.cluster_id
    fieldnames = [
        "rank", "cluster_id", "is_selected", "size", "mean_haddock_score",
        "bsa", "sc", "com_distance",
        "air_satisfaction_count", "air_total_count",
        "pair_contact_fraction", "pair_total_count",
        "bridge_distance_min", "bridge_distance_max",
        "clashes_in_design_region", "clashes_outside_design_region",
        "best_model_pdb",
    ]
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for rank, c in enumerate(clusters_sorted, start=1):
            row = c.to_summary_row()
            row["rank"] = rank
            row["is_selected"] = int(c.cluster_id == selected_id)
            w.writerow(row)


def main() -> int:
    args = parse_args()
    chosen = _parse_chosen(args.chosen_cluster_id)

    run = HaddockRun.from_workdir(
        workdir=args.workdir,
        receptor_pdb=args.receptor_pdb,
        effector_pdb=args.effector_pdb,
        contact_pairs=args.contact_pairs,
        receptor_active_residues=args.receptor_active_residues,
        effector_active_residues=args.effector_active_residues,
        pair_distance=args.pair_distance,
        chosen_cluster_id=chosen,
    )
    if not run.qualifying_clusters:
        print("ERROR: HaddockRun has no qualifying clusters", file=sys.stderr)
        return 1

    selected = run.selected
    src_pdb = selected.best_model_pdb
    shutil.copy(src_pdb, "selected_complex.pdb")
    Path("selected_cluster_id.txt").write_text(f"{selected.cluster_id}\n")
    _write_cluster_metrics_table(run, Path("cluster_metrics_table.csv"))

    mode = (
        f"USER-SPECIFIED (--chosen-cluster-id {chosen})"
        if chosen is not None else
        "AUTO-PICK by (pair_contact_fraction desc, BSA desc)"
    )
    print(f"Selected cluster {selected.cluster_id}: {mode}")
    print(f"  best_model_pdb = {src_pdb.name}")
    print(f"  BSA = {selected.bsa:.2f} A^2")
    print(f"  pair_contact_fraction = {selected.pair_contact_fraction:.2f} "
          f"({int(selected.pair_contact_fraction * selected.pair_total_count)}"
          f"/{selected.pair_total_count})")
    print(f"  AIR satisfaction = {selected.air_satisfaction_count}/"
          f"{selected.air_total_count}")
    print(f"  COM distance = {selected.com_distance:.2f} A")
    print(f"  clashes (in/out design region) = "
          f"{selected.clashes_in_design_region}/{selected.clashes_outside_design_region}")
    if selected.sc is not None:
        print(f"  Sc = {selected.sc:.3f}")
    print(f"Wrote: selected_complex.pdb, selected_cluster_id.txt, "
          f"cluster_metrics_table.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
