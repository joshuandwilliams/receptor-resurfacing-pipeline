#!/usr/bin/env python3
"""
mpnn_select_top.py
------------------
Select the top N ProteinMPNN sequences by score (lower = better).

Outputs:
    top_fastas/           - FASTA files for the top N sequences
    top_metadata.csv      - Metadata for selected sequences
    selection_report.txt  - Human-readable selection summary
"""

import argparse
import csv
import os
import shutil


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--qc-fastas-dir", required=True)
    parser.add_argument("--qc-metadata", required=True)
    parser.add_argument("--top-n", type=int, required=True)
    parser.add_argument("--output-dir", default=".")
    return parser.parse_args()


def main():
    args = parse_args()
    top_fastas_dir = os.path.join(args.output_dir, "top_fastas")
    os.makedirs(top_fastas_dir, exist_ok=True)

    rows = []
    with open(args.qc_metadata) as f:
        for row in csv.DictReader(f):
            rows.append(row)

    rows.sort(key=lambda r: float(r.get("mpnn_score", "999")))
    selected = rows[:args.top_n]

    for row in selected:
        name = f"design_{row['design']}_seq_{row['seq']}"
        src = os.path.join(args.qc_fastas_dir, f"{name}.fasta")
        if os.path.exists(src):
            shutil.copy(src, os.path.join(top_fastas_dir, f"{name}.fasta"))
        else:
            print(f"WARNING: {src} not found")

    meta_path = os.path.join(args.output_dir, "top_metadata.csv")
    if selected:
        with open(meta_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=selected[0].keys())
            writer.writeheader(); writer.writerows(selected)
    else:
        with open(meta_path, "w") as f:
            f.write("design,seq,mpnn_score\n")

    with open(os.path.join(args.output_dir, "selection_report.txt"), "w") as f:
        f.write(f"MPNN Top-N Selection Report\n{'=' * 60}\n")
        f.write(f"Total QC-passed: {len(rows)}\nSelected top {args.top_n}: {len(selected)}\n")
        if selected:
            scores = [float(r.get("mpnn_score", "0")) for r in selected]
            f.write(f"Score range: {min(scores):.3f} - {max(scores):.3f}\n")
        f.write("\nSelected:\n")
        for i, row in enumerate(selected):
            f.write(f"  {i+1}. design_{row['design']}_seq_{row['seq']} "
                    f"(score: {row.get('mpnn_score', 'N/A')})\n")

    print(f"Selected {len(selected)} / {len(rows)} (top {args.top_n} by MPNN score)")


if __name__ == "__main__":
    main()