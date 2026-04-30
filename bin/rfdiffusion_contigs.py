#!/usr/bin/env python3
"""
rfdiffusion_contigs.py
----------------------
Preprocess user-friendly contig notation into RFDiffusion native format.

Transformations:
    - "B"      → "B1-123" (resolved from PDB)
    - "A403"   → "A403-403"
    - "5"      → "5-5"
    - Fixed-segment residue numbers remapped to match input PDB

Outputs:
    processed_contigs.txt  - corrected contig string for RFDiffusion
"""

import argparse

from contig_utils import resolve_contigs


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contigs", required=True, help="Raw contig string")
    parser.add_argument("--pdb", required=True, help="Input PDB file")
    parser.add_argument("--output", default="processed_contigs.txt",
                        help="Output file (default: processed_contigs.txt)")
    return parser.parse_args()


def main():
    args = parse_args()
    processed = resolve_contigs(args.contigs, args.pdb)

    print(f"Raw contigs:       {args.contigs}")
    print(f"Processed contigs: {processed}")

    with open(args.output, "w") as f:
        f.write(processed)


if __name__ == "__main__":
    main()