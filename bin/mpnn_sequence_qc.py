#!/usr/bin/env python3
"""
mpnn_sequence_qc.py
-------------------
Quality-control filter for ProteinMPNN-designed sequences.

Checks:
  - Poly-X homopolymer runs exceeding a configurable length
  - Unusual amino acid compositions (>30% single identity, or non-standard AAs)
  - Percentage identity to native within configurable bounds
  - Duplicate receptor design sequences across all designs (keeps best MPNN score)

Outputs:
    qc_fastas/              - FASTA files that passed all checks
    qc_metadata.csv         - Metadata for passing sequences
    qc_report.txt           - Human-readable pass/fail summary
"""

import argparse
import csv
import os
import shutil


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--af2-fastas-dir", required=True)
    parser.add_argument("--metadata-csv", required=True)
    parser.add_argument("--receptor-seq", required=True)
    parser.add_argument("--max-poly-x", type=int, default=5)
    parser.add_argument("--min-pct-identity", type=float, default=0.0)
    parser.add_argument("--max-pct-identity", type=float, default=100.0)
    parser.add_argument("--max-single-aa-frac", type=float, default=0.30)
    parser.add_argument("--no-dedup", action="store_true",
                        help="Skip duplicate receptor sequence removal")
    parser.add_argument("--output-dir", default=".")
    return parser.parse_args()


def has_poly_x(seq, max_run):
    for aa in set(seq):
        if aa * (max_run + 1) in seq:
            return True, aa
    return False, None


def unusual_aa_composition(seq, max_single_aa_frac=0.30):
    common = set("AVILMFYWKRHDENQSTGPC")
    for aa in set(seq):
        if aa not in common:
            return True, aa
    total = len(seq)
    if total == 0:
        return False, None
    for aa in set(seq):
        frac = seq.count(aa) / total
        if frac > max_single_aa_frac and aa not in "ALG":
            return True, aa
    return False, None


def main():
    args = parse_args()
    qc_fastas_dir = os.path.join(args.output_dir, "qc_fastas")
    os.makedirs(qc_fastas_dir, exist_ok=True)

    rows = []
    with open(args.metadata_csv) as f:
        for row in csv.DictReader(f):
            rows.append(row)

    passed, failed = [], []
    report_lines = ["Sequence QC Report", "=" * 60]

    for row in rows:
        design, seq_n = row["design"], row["seq"]
        name = f"design_{design}_seq_{seq_n}"
        fasta_path = os.path.join(args.af2_fastas_dir, f"{name}.fasta")

        if not os.path.exists(fasta_path):
            failed.append((name, "FASTA file missing")); continue

        corrected_receptor = None
        with open(fasta_path) as f:
            header = None
            for line in f:
                line = line.strip()
                if line.startswith(">"): header = line
                elif header and "receptor" in header.lower():
                    corrected_receptor = line; break
        if not corrected_receptor:
            corrected_receptor = row.get("corrected_receptor", "")

        has_poly, poly_aa = has_poly_x(corrected_receptor, args.max_poly_x)
        if has_poly:
            failed.append((name, f"Poly-{poly_aa} run > {args.max_poly_x}")); continue

        unusual, unusual_aa = unusual_aa_composition(corrected_receptor, args.max_single_aa_frac)
        if unusual:
            report_lines.append(f"  WARNING: {name} has unusual AA composition ({unusual_aa})")

        designed = row.get("designed_residues", "")
        native = row.get("native_residues", "")
        # Strip region separators for identity comparison
        designed_flat = designed.replace("|", "")
        native_flat = native.replace("|", "")
        if designed_flat and native_flat and len(designed_flat) == len(native_flat):
            matches = sum(1 for d, n in zip(designed_flat, native_flat) if d == n)
            pct_id = (matches / len(native_flat)) * 100 if native_flat else 100
            if pct_id < args.min_pct_identity:
                failed.append((name, f"% identity {pct_id:.1f}% < {args.min_pct_identity}%")); continue
            if pct_id >= args.max_pct_identity:
                failed.append((name, f"% identity {pct_id:.1f}% >= {args.max_pct_identity}%")); continue

        passed.append(row)
        shutil.copy(fasta_path, os.path.join(qc_fastas_dir, f"{name}.fasta"))

    # ── Deduplicate identical receptor design sequences ────────────────
    # Two designs can independently produce the same corrected receptor
    # sequence.  Keep only the one with the best (lowest) MPNN score.
    n_before_dedup = len(passed)
    duplicates = []
    if not args.no_dedup and passed:
        seen_seqs = {}  # corrected_receptor -> (index_in_passed, mpnn_score)
        dedup_indices_to_remove = set()
        for idx, row in enumerate(passed):
            rec_seq = row.get("corrected_receptor", "")
            score = float(row.get("mpnn_score", "999"))
            if rec_seq in seen_seqs:
                prev_idx, prev_score = seen_seqs[rec_seq]
                # Keep the one with the lower (better) MPNN score
                if score < prev_score:
                    # Current is better — remove previous
                    dedup_indices_to_remove.add(prev_idx)
                    seen_seqs[rec_seq] = (idx, score)
                    dup_name = f"design_{passed[prev_idx]['design']}_seq_{passed[prev_idx]['seq']}"
                    kept_name = f"design_{row['design']}_seq_{row['seq']}"
                    duplicates.append((dup_name, f"Duplicate receptor seq (kept {kept_name}, score {score:.3f})"))
                else:
                    # Previous is better — remove current
                    dedup_indices_to_remove.add(idx)
                    kept_name = f"design_{passed[prev_idx]['design']}_seq_{passed[prev_idx]['seq']}"
                    dup_name = f"design_{row['design']}_seq_{row['seq']}"
                    duplicates.append((dup_name, f"Duplicate receptor seq (kept {kept_name}, score {prev_score:.3f})"))
            else:
                seen_seqs[rec_seq] = (idx, score)

        # Remove duplicate FASTAs from qc_fastas/
        for idx in dedup_indices_to_remove:
            dup_row = passed[idx]
            dup_fasta = os.path.join(
                qc_fastas_dir,
                f"design_{dup_row['design']}_seq_{dup_row['seq']}.fasta")
            if os.path.exists(dup_fasta):
                os.remove(dup_fasta)

        # Filter the passed list
        passed = [row for idx, row in enumerate(passed)
                  if idx not in dedup_indices_to_remove]

    n_dedup_removed = n_before_dedup - len(passed)

    qc_meta_path = os.path.join(args.output_dir, "qc_metadata.csv")
    if passed:
        with open(qc_meta_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=passed[0].keys())
            writer.writeheader(); writer.writerows(passed)
    else:
        with open(qc_meta_path, "w") as f:
            f.write("design,seq,mpnn_score,native_residues,designed_residues,num_changes,corrected_receptor\n")

    report_lines.extend([f"\nTotal input sequences: {len(rows)}",
                         f"Passed QC: {n_before_dedup}",
                         f"Duplicates removed: {n_dedup_removed}",
                         f"Final unique sequences: {len(passed)}",
                         f"Failed QC: {len(failed)}",
                         "\nFailed sequences:"])
    for name, reason in failed:
        report_lines.append(f"  {name}: {reason}")
    if duplicates:
        report_lines.append("\nDuplicate sequences removed:")
        for name, reason in duplicates:
            report_lines.append(f"  {name}: {reason}")
    with open(os.path.join(args.output_dir, "qc_report.txt"), "w") as f:
        f.write("\n".join(report_lines) + "\n")

    if not passed:
        print("ERROR: No sequences passed QC!")
        with open(os.path.join(qc_fastas_dir, "NONE_PASSED_QC.txt"), "w") as f:
            f.write("No sequences passed quality control.\n")
    print(f"QC complete: {len(passed)} passed ({n_dedup_removed} duplicates removed), {len(failed)} failed")


if __name__ == "__main__":
    main()