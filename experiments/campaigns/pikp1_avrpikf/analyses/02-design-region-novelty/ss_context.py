"""What secondary structure does each design region replace?

Region size and anchoring are only two of the three things that plausibly set
how much freedom RFdiffusion has. The third is what the region is cut out of.
A region replacing a strand of the central sheet is held by backbone hydrogen
bonds to its neighbouring strands, whereas a region replacing a loop is not.

The receptor's secondary structure is taken from the HELIX and SHEET records of
the input complex rather than recomputed, and the residues each region replaces
follow from the contig alone, so this annotation is the same for every design.

Pikp-1 HMA is a ferredoxin-like babbab fold. The records give
    b1 3-10, a1 14-26, b2 31-37, b3 43-49, a2 53-63, b4 67-74
so everything outside those spans is loop or terminus.
"""
from __future__ import annotations

import os
from pathlib import Path

import pandas as pd

import motif_rigidity as mr
import novelty_common as nc

WT_LEN = 78


def parse_ss(pdb: Path, chain: str = "A") -> dict[int, str]:
    """resnum -> 'H' | 'E' | 'C', from the PDB's own HELIX/SHEET records."""
    ss = {}
    for line in pdb.read_text().splitlines():
        if line.startswith("HELIX") and line[19] == chain:
            lo, hi, code = int(line[21:25]), int(line[33:37]), "H"
        elif line.startswith("SHEET") and line[21] == chain:
            lo, hi, code = int(line[22:26]), int(line[33:37]), "E"
        else:
            continue
        for r in range(lo, hi + 1):
            ss[r] = code
    return {r: ss.get(r, "C") for r in range(1, WT_LEN + 1)}


def replaced_residues(contig: str) -> list[list[int]]:
    """Input residues each design region replaces, one list per region.

    A region between two retained blocks replaces the gap between them. A
    trailing region replaces everything from the last retained residue to the
    end of the wild-type domain.
    """
    blocks = mr.fixed_blocks(contig)
    segs = contig.split()[0].split("/")
    out, b = [], 0
    for seg in segs:
        if seg.startswith("A"):
            b += 1
            continue
        left = blocks[b - 1][1] if b > 0 else 0
        right = blocks[b][0] if b < len(blocks) else None
        out.append(list(range(left + 1, right if right is not None else WT_LEN + 1)))
    return out


def annotate() -> pd.DataFrame:
    """One row per (run, region) with the composition of what it replaces."""
    ss = parse_ss(mr.PDB_ROOT / "inputs" / mr.INPUT_PDB)
    rows = []
    for label, contig in nc.CONTIGS.items():
        n_regions = sum(1 for s in contig.split()[0].split("/")
                        if not s.startswith("A"))
        for i, res in enumerate(replaced_residues(contig), start=1):
            codes = [ss[r] for r in res if r in ss]
            n = len(codes) or 1
            frac = {c: codes.count(c) / n for c in "HEC"}
            rows.append({
                "run": label, "region": i,
                "region_id": f"{label}|{i}",
                "terminal": int(i == n_regions),
                "replaced_first": res[0] if res else None,
                "replaced_last": res[-1] if res else None,
                "n_replaced": len(res),
                "frac_strand": frac["E"], "frac_helix": frac["H"],
                "frac_loop": frac["C"],
                "ss_class": _classify(frac),
            })
    return pd.DataFrame(rows)


def _classify(frac: dict[str, float]) -> str:
    """Majority call, with anything under 60% of one type called mixed."""
    code, f = max(frac.items(), key=lambda kv: kv[1])
    if f < 0.6:
        return "mixed"
    return {"E": "strand", "H": "helix", "C": "loop"}[code]


if __name__ == "__main__":
    pd.set_option("display.width", 170)
    df = annotate()
    print(df.to_string(index=False))
    df.to_csv(Path(__file__).parent / "region_ss_context.csv", index=False)
