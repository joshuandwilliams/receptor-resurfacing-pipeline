#!/usr/bin/env python3
"""
Fetch Pikp1 target effector structures from RCSB PDB.
Extracts effector chains and saves as individual PDB files.
"""

import sys
from pathlib import Path

import requests

# Target definitions: (allele_name, pdb_id, chain_id, description)
TARGETS = [
    ("pia", "6Q76", "B", "AVR-Pia (Varden 2019)"),
    ("pikd", "6G10", "C", "AVR-PikD (De la Concepcion 2018) - highest res"),
    ("pike", "6G11", "C", "AVR-PikE (De la Concepcion 2018)"),
    ("pikc", "7QPX", "C", "AVR-PikC (Maidment 2023)"),
    ("pikf", "7QZD", "C", "AVR-PikF (Maidment 2023)"),
    ("co39", "5ZNG", "C", "AVR1-CO39 (Guo 2018)"),
]

def download_pdb(pdb_id: str) -> str:
    """Download PDB file from RCSB."""
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    print(f"  Downloading {pdb_id} from {url}...", end=" ", flush=True)
    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        print("✓")
        return response.text
    except Exception as e:
        print(f"✗ Error: {e}")
        return None

def extract_chain(pdb_content: str, chain_id: str) -> str:
    """Extract specific chain from PDB content."""
    lines = pdb_content.split("\n")
    output_lines = []

    # Copy header/remarks
    for line in lines:
        if line.startswith("REMARK") or line.startswith("HEADER"):
            output_lines.append(line)

    # Extract chain atoms
    for line in lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            if len(line) >= 22 and line[21] == chain_id:
                output_lines.append(line)

    # Add END
    output_lines.append("END")
    return "\n".join(output_lines)

def main():
    if len(sys.argv) < 2:
        print("Usage: python fetch_pikp1_targets.py <output_directory>")
        print("Example: python fetch_pikp1_targets.py '/Users/jowillia/Documents/GitHub/receptor-resurfacing-pipeline/experiments/inputs/structures/Pikp1 targets'")
        sys.exit(1)

    output_dir = Path(sys.argv[1])
    output_dir.mkdir(parents=True, exist_ok=True)

    print("\nFetching Pikp1 target effector structures...")
    print(f"Output directory: {output_dir}\n")

    success_count = 0
    for allele, pdb_id, chain_id, description in TARGETS:
        print(f"{allele.upper():6s} ({description})")

        # Download
        pdb_content = download_pdb(pdb_id)
        if not pdb_content:
            continue

        # Extract chain
        print(f"  Extracting chain {chain_id}...", end=" ", flush=True)
        try:
            chain_content = extract_chain(pdb_content, chain_id)

            # Write output
            output_file = output_dir / f"avr_{allele}.pdb"
            output_file.write_text(chain_content)
            print(f"✓ → {output_file.name}")
            success_count += 1
        except Exception as e:
            print(f"✗ Error: {e}")

    print(f"\n✓ Fetched {success_count}/{len(TARGETS)} structures")

if __name__ == "__main__":
    main()
