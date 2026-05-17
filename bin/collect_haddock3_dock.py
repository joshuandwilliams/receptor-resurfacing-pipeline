#!/usr/bin/env python3
"""
collect_haddock3_dock.py
------------------------
Post-process a completed HADDOCK3 docking run.

Locates the best-scoring model PDB for every qualifying cluster
(>= ``--min-cluster-size`` members), copies CAPRI and cluster summary
files, and writes ``haddock_report.json``.

Per Session 7 restructure (notes/design_audit.md A129/A130/A134):
- No "best model" is selected here.  Cluster ranking by HADDOCK score has
  been replaced by ranking on BSA + pair contact fraction, which is
  computed by the downstream ``haddock_cluster_metrics.py`` step and
  used by ``select_haddock_cluster.py`` to pick the chosen cluster.
- The fail-if-best-score-positive guard has also been removed: with
  multiple qualifying clusters reported and BSA driving the choice, a
  single positive HADDOCK score isn't grounds to fail the whole run.
  (Pipeline error if NO clusters qualify is preserved.)

Outputs:
    best_cluster{N}.pdb   - best model PDB from each qualifying cluster
    capri_scores.tsv      - CAPRI evaluation scores (copy)
    cluster_summary.txt   - cluster membership table (copy, optional)
    haddock_report.json   - run summary with per-cluster metadata
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import shutil
import sys
import tarfile

from haddock_utils import (
    parse_capri_tsv, parse_clustfcc_tsv,
    get_score, get_model_name, model_stem, copy_pdb,
    cluster_mean_score,
)


MIN_CLUSTER_SIZE = 4  # Default; overridden by --min-cluster-size CLI flag.
                      # Single source of truth lives in nextflow.config as
                      # params.haddock_min_cluster_size, plumbed via haddock.nf.


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", default="run/run-haddock",
                        help="Path to the HADDOCK3 run directory")
    parser.add_argument("--min-cluster-size", type=int, default=MIN_CLUSTER_SIZE,
                        help="Minimum cluster size to qualify (default: %(default)s). "
                             "Should match `min_population` in haddock.nf clustfcc block.")
    return parser.parse_args()


def find_file(run_dir, analysis_pattern, step_pattern, filename):
    """Search for a file in analysis/, then step dir, then recursively."""
    for pattern in [
        os.path.join(run_dir, "analysis", analysis_pattern, filename),
        os.path.join(run_dir, step_pattern, filename),
    ]:
        paths = glob.glob(pattern)
        if paths:
            return paths[0]
    paths = glob.glob(os.path.join(run_dir, "**", filename), recursive=True)
    return paths[0] if paths else None


def select_all_qualifying_clusters(clusters):
    """
    Return a list of dicts for every qualifying cluster (>= MIN_CLUSTER_SIZE,
    excluding '-'), sorted by mean HADDOCK score (best first).

    Each dict: {cluster_id, size, mean_score, best_model_name, best_model_score}

    The mean-HADDOCK-score sort here is purely cosmetic for the report
    output order — the actual cluster ranking happens downstream on BSA +
    pair contact fraction.
    """
    results = []
    for cid, members in clusters.items():
        if cid == "-" or len(members) < MIN_CLUSTER_SIZE:
            continue
        mean = cluster_mean_score(members)
        scored = [(m, s) for (m, s) in members if s is not None]
        if scored:
            best_name, best_score = min(scored, key=lambda x: x[1])
        else:
            best_name, best_score = members[0][0], None

        results.append({
            "cluster_id": cid,
            "size": len(members),
            "mean_score": mean,
            "best_model_name": best_name,
            "best_model_score": best_score,
        })

    results.sort(key=lambda r: r["mean_score"])
    return results


def extract_model_by_name(run_dir, target_model_name, output_name):
    """
    Locate a model PDB by stem match in seletopclusts/emref directories.
    Falls back to first available model, then summary.tgz.
    Returns True on success.
    """
    target = model_stem(target_model_name)

    search_dirs = [
        sorted(glob.glob(os.path.join(run_dir, "*seletopclusts*"))),
        sorted(glob.glob(os.path.join(run_dir, "*emref*"))),
    ]

    # Pass 1: exact stem match
    for dir_list in search_dirs:
        for d in dir_list:
            if not os.path.isdir(d):
                continue
            for pattern in ("*.pdb.gz", "*.pdb"):
                for fpath in sorted(glob.glob(os.path.join(d, pattern))):
                    if model_stem(fpath) == target:
                        copy_pdb(fpath, output_name)
                        print(f"Extracted {output_name} ({target}) from: {fpath}")
                        return True

    # Pass 2: first available model
    print(f"WARNING: No exact match for '{target}'; using first available model.",
          file=sys.stderr)
    for dir_list in search_dirs:
        for d in dir_list:
            if not os.path.isdir(d):
                continue
            for pattern in ("*.pdb.gz", "*.pdb"):
                files = sorted(glob.glob(os.path.join(d, pattern)))
                if files:
                    copy_pdb(files[0], output_name)
                    print(f"Extracted fallback model to {output_name} from: {files[0]}")
                    return True

    # Pass 3: summary.tgz
    tgz = os.path.join(run_dir, "analysis", "summary.tgz")
    if os.path.exists(tgz):
        with tarfile.open(tgz) as tar:
            members = [m for m in tar.getmembers() if m.name.endswith(".pdb")]
            if members:
                tar.extract(members[0])
                shutil.move(members[0].name, output_name)
                print(f"Extracted fallback model to {output_name} from summary.tgz")
                return True

    return False


def _fail(report, message):
    """Write failure report, then exit 1."""
    print(f"ERROR: {message}", file=sys.stderr)
    with open("haddock_report.json", "w") as f:
        json.dump(report, f, indent=2)
    sys.exit(1)


def main():
    args = parse_args()
    # Re-bind the module-level constant from the CLI flag so all existing
    # references (including from helper functions) pick up the override
    # without needing to thread the value through every function signature.
    global MIN_CLUSTER_SIZE
    MIN_CLUSTER_SIZE = args.min_cluster_size
    run_dir = args.run_dir

    report = {
        "success": False, "n_clusters": 0,
        "n_qualifying_clusters": 0, "cluster_models": {},
    }

    # ── CAPRI scores ─────────────────────────────────────────────────────
    capri_file = find_file(run_dir, "*caprieval*", "*caprieval*", "capri_ss.tsv")
    if not capri_file:
        _fail(report, "No capri_ss.tsv found")

    print(f"Using CAPRI file: {capri_file}")
    shutil.copy(capri_file, "capri_scores.tsv")

    capri_rows, _ = parse_capri_tsv(capri_file)
    if not capri_rows:
        _fail(report, "capri_ss.tsv is empty")

    # ── Cluster summary ──────────────────────────────────────────────────
    clust_file = find_file(run_dir, "*clustfcc*", "*clustfcc*", "clustfcc.tsv")
    clusters = {}
    if clust_file:
        shutil.copy(clust_file, "cluster_summary.txt")
        clusters, n_clusters, _ = parse_clustfcc_tsv(clust_file, MIN_CLUSTER_SIZE)
        report["n_clusters"] = n_clusters
        print(f"Found {n_clusters} cluster(s) (excluding unclustered models)")
    else:
        print("WARNING: No clustfcc.tsv found; cannot perform cluster-aware "
              "selection — pipeline will halt", file=sys.stderr)
        _fail(report, "No clustfcc.tsv found; HADDOCK3 must have crashed during "
                      "clustering")

    qualifying = select_all_qualifying_clusters(clusters)
    if not qualifying:
        sizes = {cid: len(m) for cid, m in clusters.items() if cid != "-"}
        size_str = ", ".join(f"cluster {c}: {s}" for c, s in sorted(sizes.items()))
        _fail(report, f"No cluster has >= {MIN_CLUSTER_SIZE} models. "
                      f"Cluster sizes: {size_str or 'none'}")

    report["n_qualifying_clusters"] = len(qualifying)
    print(f"{len(qualifying)} qualifying cluster(s) "
          f"(min size {MIN_CLUSTER_SIZE}); extracting best model from each")

    # ── Extract best model from every qualifying cluster ─────────────────
    cluster_models = {}
    for info in qualifying:
        cid = info["cluster_id"]
        filename = f"best_cluster{cid}.pdb"
        if extract_model_by_name(run_dir, info["best_model_name"], output_name=filename):
            cluster_models[cid] = {
                "filename": filename,
                "model_name": info["best_model_name"],
                "score": info["best_model_score"],
                "mean_score": round(info["mean_score"], 3),
                "size": info["size"],
            }
            print(f"  Cluster {cid}: {filename} "
                  f"(score={info['best_model_score']}, "
                  f"mean={info['mean_score']:.2f}, N={info['size']})")
        else:
            print(f"  WARNING: Could not extract best model for cluster {cid}",
                  file=sys.stderr)

    if not cluster_models:
        _fail(report, "Failed to extract any cluster's best-model PDB")

    report["success"] = True
    report["cluster_models"] = cluster_models

    with open("haddock_report.json", "w") as f:
        json.dump(report, f, indent=2)

    print(f"SUCCESS: {len(cluster_models)} per-cluster best model(s) "
          f"saved; downstream haddock_cluster_metrics.py will rank by BSA + "
          f"pair contact fraction")


if __name__ == "__main__":
    main()
