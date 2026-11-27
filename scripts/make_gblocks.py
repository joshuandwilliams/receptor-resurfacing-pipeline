#!/usr/bin/env python3
"""
Reverse translate amino acid sequences to E. coli codon-optimised DNA
and append Golden Gate primer tails for gBlock ordering.

Uses the E. coli K12 codon usage table from:
  https://www.bioinformatics.org/sms2/rev_trans.html
  (Kazusa codon usage database)

Picks the most likely codon (highest fraction) for each amino acid.

Usage:
    python make_gblocks.py <input_csv> [output_csv]

Input CSV must contain a 'corrected_receptor' column with amino acid sequences.
Adds a 'gblock_sequence' column with the full gBlock sequence ready to order.
"""

import sys

import pandas as pd

# E. coli K12 most likely codons (highest fraction per amino acid)
# Source: Kazusa codon usage database via SMS2 rev_trans tool
ECOLI_BEST_CODONS = {
    'A': 'gcg',  # Ala: 0.34
    'C': 'tgc',  # Cys: 0.55
    'D': 'gat',  # Asp: 0.63
    'E': 'gaa',  # Glu: 0.68
    'F': 'ttt',  # Phe: 0.57
    'G': 'ggc',  # Gly: 0.39
    'H': 'cat',  # His: 0.57
    'I': 'att',  # Ile: 0.50
    'K': 'aaa',  # Lys: 0.75
    'L': 'ctg',  # Leu: 0.49
    'M': 'atg',  # Met: 1.00
    'N': 'aac',  # Asn: 0.54
    'P': 'ccg',  # Pro: 0.51
    'Q': 'cag',  # Gln: 0.66
    'R': 'cgc',  # Arg: 0.38
    'S': 'agc',  # Ser: 0.26
    'T': 'acc',  # Thr: 0.42
    'V': 'gtg',  # Val: 0.36
    'W': 'tgg',  # Trp: 1.00
    'Y': 'tat',  # Tyr: 0.58
    '*': 'taa',  # Stop: 0.62
}

# Golden Gate primer tails for pICSL system (AATG / GCTT overhangs)
FORWARD_TAIL = 'gacggtctcgaatg'
REVERSE_TAIL = 'gcttcgagaccctt'


def reverse_translate(aa_sequence: str) -> str:
    """Convert amino acid sequence to E. coli codon-optimised DNA (no stop codon)."""
    dna = ''
    for i, aa in enumerate(aa_sequence.upper()):
        if aa not in ECOLI_BEST_CODONS:
            raise ValueError(f"Unknown amino acid '{aa}' at position {i+1}")
        dna += ECOLI_BEST_CODONS[aa]
    return dna


def make_gblock(aa_sequence: str) -> str:
    """Full gBlock sequence: forward tail + CDS + reverse tail. No stop codon."""
    cds = reverse_translate(aa_sequence)
    return FORWARD_TAIL + cds + REVERSE_TAIL


def main():
    if len(sys.argv) < 2:
        print("Usage: python make_gblocks.py <input_csv> [output_csv]")
        sys.exit(1)

    input_csv = sys.argv[1]
    output_csv = sys.argv[2] if len(sys.argv) > 2 else input_csv.replace('.csv', '_gblocks.csv')

    df = pd.read_csv(input_csv)

    if 'corrected_receptor' not in df.columns:
        print("Error: 'corrected_receptor' column not found.")
        print(f"Available columns: {list(df.columns)}")
        sys.exit(1)

    df['gblock_sequence'] = df['corrected_receptor'].apply(
        lambda x: make_gblock(x) if pd.notna(x) else ''
    )

    df.to_csv(output_csv, index=False)
    print(f"Done. {len(df)} sequences processed.")
    print(f"Output: {output_csv}")

    if 'mpnn_sequence' not in df.columns:
        print("Warning: 'mpnn_sequence' column not found — skipping sequences-only output.")
        return

    sequences_csv = output_csv.replace('_gblocks.csv', '_sequences.csv')
    df[['mpnn_sequence', 'gblock_sequence']].to_csv(sequences_csv, index=False)
    print(f"Sequences output: {sequences_csv}")


if __name__ == '__main__':
    main()
