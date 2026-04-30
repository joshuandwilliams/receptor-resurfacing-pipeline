#!/usr/bin/env python3
"""
collect_haddock3_dock.py
------------------------
Post-process a completed HADDOCK3 docking run.

Locates the best-scoring model PDB from the best qualifying cluster, copies
CAPRI and cluster summary files, and writes a JSON report.

Cluster selection:
    1. Only clusters with >= 4 models qualify.
    2. Best cluster = lowest mean HADDOCK score among qualifying clusters.
    3. Best model = lowest-scoring individual model within that cluster.
    4. Fails if no cluster qualifies or best model score >= 0.

Outputs:
    best_model.pdb        - best-scoring docked complex from best cluster
    best_cluster{N}.pdb   - best model from each qualifying cluster
    capri_scores.tsv      - CAPRI evaluation scores (copy)
    cluster_summary.txt   - cluster membership table (copy, optional)
    haddock_report.json   - run summary with per-cluster metadata
"""

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


def select_best_cluster(clusters):
    """
    Select the single best qualifying cluster.
    Returns (cluster_id, model_name, model_score, mean_score, size)
    or exits if no cluster qualifies.
    """
    all_qual = select_all_qualifying_clusters(clusters)
    if not all_qual:
        sizes = {cid: len(m) for cid, m in clusters.items() if cid != "-"}
        size_str = ", ".join(f"cluster {c}: {s}" for c, s in sorted(sizes.items()))
        print(f"ERROR: No cluster has >= {MIN_CLUSTER_SIZE} models. "
              f"Cluster sizes: {size_str or 'none'}", file=sys.stderr)
        sys.exit(1)

    best = all_qual[0]
    return (best["cluster_id"], best["best_model_name"], best["best_model_score"],
            best["mean_score"], best["size"])


def extract_model_by_name(run_dir, target_model_name, output_name="best_model.pdb"):
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
    """Write failure report + placeholder PDB, then exit 1."""
    print(f"ERROR: {message}", file=sys.stderr)
    with open("haddock_report.json", "w") as f:
        json.dump(report, f, indent=2)
    with open("best_model.pdb", "w") as f:
        f.write(f"REMARK HADDOCK3 failed - {message}\n")
    sys.exit(1)


def _validate_score(score, report):
    """Fail if score is missing or non-negative (repulsive dock)."""
    if score is None:
        _fail(report, "Best model has no HADDOCK score — cannot validate dock quality")
    if score >= 0:
        _fail(report, f"Best model HADDOCK score is {score:.2f} (>= 0). "
              f"Positive scores indicate repulsive/non-specific docking poses.")


def main():
    args = parse_args()
    # Re-bind the module-level constant from the CLI flag so all existing
    # references (including from helper functions) pick up the override
    # without needing to thread the value through every function signature.
    global MIN_CLUSTER_SIZE
    MIN_CLUSTER_SIZE = args.min_cluster_size
    run_dir = args.run_dir

    report = {
        "success": False, "best_model": None, "best_score": None,
        "best_cluster_id": None, "best_cluster_size": 0,
        "best_cluster_mean_score": None, "n_clusters": 0,
        "n_qualifying_clusters": 0,
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
        print("WARNING: No clustfcc.tsv found; cannot perform cluster-aware selection",
              file=sys.stderr)

    # ── Select best model ────────────────────────────────────────────────
    if clusters:
        (best_cid, best_model_name, best_model_score,
         best_mean, best_size) = select_best_cluster(clusters)

        all_qual = select_all_qualifying_clusters(clusters)
        report["n_qualifying_clusters"] = len(all_qual)
        report["best_cluster_id"] = best_cid
        report["best_cluster_size"] = best_size
        report["best_cluster_mean_score"] = round(best_mean, 3)

        print(f"Best cluster: {best_cid} (size={best_size}, mean score={best_mean:.2f})")
        print(f"Best model in cluster: {best_model_name} (score={best_model_score})")

        _validate_score(best_model_score, report)
        report["best_model"] = best_model_name
        report["best_score"] = best_model_score
    else:
        # No cluster file — fall back to best row in capri_ss.tsv
        print("WARNING: No cluster data; selecting globally best model from capri_ss.tsv")
        best_model_name = get_model_name(capri_rows[0])
        best_model_score = get_score(capri_rows[0])
        _validate_score(best_model_score, report)
        report["best_model"] = best_model_name
        report["best_score"] = best_model_score

    # ── Extract best model PDB ───────────────────────────────────────────
    if not extract_model_by_name(run_dir, best_model_name):
        _fail(report, "Could not find best model PDB on disk")

    # ── Extract best model from every qualifying cluster ─────────────────
    cluster_models = {}
    if clusters:
        for info in select_all_qualifying_clusters(clusters):
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

    report["success"] = True
    report["cluster_models"] = cluster_models

    with open("haddock_report.json", "w") as f:
        json.dump(report, f, indent=2)

    print(f"SUCCESS: best_model.pdb written "
          f"(cluster={report.get('best_cluster_id', 'N/A')}, "
          f"score={report['best_score']:.2f}); "
          f"{len(cluster_models)} per-cluster model(s) saved")


if __name__ == "__main__":
    main()