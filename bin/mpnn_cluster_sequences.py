#!/usr/bin/env python3
"""
mpnn_cluster_sequences.py
-------------------------
Cluster ProteinMPNN-designed sequences at multiple identity thresholds
using MMseqs2 easy-cluster, and report the number of clusters at each.

Clusters the design-region sequences (variable length) directly.
MMseqs2 handles length differences via local alignment.

Outputs:
    mpnn_cluster_counts.csv  - threshold, n_clusters, n_sequences
"""

import argparse
import csv
import os
import shutil
import subprocess
import sys
import tempfile


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata-csv", required=True,
                        help="QC metadata CSV with designed_residues column")
    parser.add_argument("--output-csv", default="mpnn_cluster_counts.csv",
                        help="Output CSV (default: mpnn_cluster_counts.csv)")
    parser.add_argument("--thresholds", default="0.30,0.40,0.50,0.60,0.70,0.80,0.90,1.00",
                        help="Comma-separated identity thresholds (default: 0.30-1.00)")
    parser.add_argument("--mmseqs-bin", default="mmseqs",
                        help="Path to mmseqs binary (default: mmseqs)")
    return parser.parse_args()


def write_design_fasta(rows, fasta_path):
    """Write design-region sequences to a FASTA file.

    Each sequence is named by design_seq for traceability.
    """
    n_written = 0
    with open(fasta_path, "w") as f:
        for row in rows:
            designed = row.get("designed_residues", "")
            if not designed:
                continue
            # Strip region separators for sequence clustering
            designed_flat = designed.replace("|", "")
            if not designed_flat:
                continue
            name = f"design_{row['design']}_seq_{row['seq']}"
            f.write(f">{name}\n{designed_flat}\n")
            n_written += 1
    return n_written


def run_mmseqs_cluster(fasta_path, threshold, mmseqs_bin, work_dir):
    """Run mmseqs easy-cluster at a given identity threshold.

    Returns the number of clusters.
    """
    prefix = os.path.join(work_dir, "clust")
    tmp_dir = os.path.join(work_dir, "tmp")
    os.makedirs(tmp_dir, exist_ok=True)

    # Clean previous results
    for ext in ["_cluster.tsv", "_rep_seq.fasta", "_all_seqs.fasta"]:
        p = prefix + ext
        if os.path.exists(p):
            os.remove(p)

    cmd = [
        mmseqs_bin, "easy-cluster",
        fasta_path,
        prefix,
        tmp_dir,
        "--min-seq-id", f"{threshold:.2f}",
        "-c", "0.8",           # coverage threshold
        "--cov-mode", "1",     # coverage of shorter sequence
        "--cluster-mode", "0", # greedy set cover
        "-v", "0",             # quiet
    ]

    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"  WARNING: mmseqs failed at threshold {threshold}: {e.stderr}",
              file=sys.stderr)
        return None

    # Count clusters from the cluster TSV
    cluster_tsv = prefix + "_cluster.tsv"
    if not os.path.exists(cluster_tsv):
        return None

    representatives = set()
    with open(cluster_tsv) as f:
        for line in f:
            parts = line.strip().split("\t")
            if parts:
                representatives.add(parts[0])

    return len(representatives)


def main():
    args = parse_args()

    # Read metadata
    rows = []
    with open(args.metadata_csv) as f:
        for row in csv.DictReader(f):
            rows.append(row)

    if not rows:
        print("WARNING: No sequences in metadata")
        with open(args.output_csv, "w") as f:
            f.write("threshold,n_clusters,n_sequences\n")
        return

    thresholds = [float(t) for t in args.thresholds.split(",")]
    thresholds.sort()

    # Write FASTA
    work_dir = tempfile.mkdtemp(prefix="mpnn_cluster_")
    fasta_path = os.path.join(work_dir, "design_regions.fasta")
    n_seqs = write_design_fasta(rows, fasta_path)

    if n_seqs < 2:
        print(f"Only {n_seqs} sequences — skipping clustering")
        with open(args.output_csv, "w") as f:
            f.write("threshold,n_clusters,n_sequences\n")
            for t in thresholds:
                f.write(f"{t:.2f},{n_seqs},{n_seqs}\n")
        shutil.rmtree(work_dir, ignore_errors=True)
        return

    print(f"Clustering {n_seqs} design-region sequences at "
          f"{len(thresholds)} thresholds")

    results = []
    for threshold in thresholds:
        # MMseqs2 needs a clean tmp dir each run
        cluster_work = os.path.join(work_dir, f"t{int(threshold*100)}")
        os.makedirs(cluster_work, exist_ok=True)

        n_clusters = run_mmseqs_cluster(
            fasta_path, threshold, args.mmseqs_bin, cluster_work)

        if n_clusters is None:
            print(f"  {threshold:.0%}: FAILED")
            n_clusters = n_seqs  # fallback: assume all unique
        else:
            print(f"  {threshold:.0%}: {n_clusters} clusters")

        results.append({
            "threshold": f"{threshold:.2f}",
            "n_clusters": n_clusters,
            "n_sequences": n_seqs,
        })

    # Write output
    with open(args.output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["threshold", "n_clusters", "n_sequences"])
        writer.writeheader()
        writer.writerows(results)

    # Cleanup
    shutil.rmtree(work_dir, ignore_errors=True)

    print(f"\nCluster counts written to {args.output_csv}")


if __name__ == "__main__":
    main()