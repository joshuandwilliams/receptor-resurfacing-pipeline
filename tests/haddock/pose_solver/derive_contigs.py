"""
Derive RFDiffusion contigs from a 2-chain complex at a given heavy-atom
contact cutoff, using the same conventions as pose_solved_5a:
  - design region = binder residues with any heavy-atom contact <= cutoff to other chain
  - join-gap <= 1 (segments separated by a single fixed residue are merged)
  - strategy-4a length-range scaling: per segment of length N -> [round(N*0.7), round(N*1.5)],
    with N=1 staying 1-1.

Usage:
  python derive_contigs.py PDB CHAIN_A CHAIN_B CUTOFF [CUTOFF ...]
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np


def read_heavy(pdb: Path, chain: str) -> dict[int, np.ndarray]:
    data: dict[int, list] = {}
    for ln in pdb.read_text().splitlines():
        if not ln.startswith("ATOM"): continue
        if ln[21] != chain: continue
        atom = ln[12:16].strip()
        elem = ln[76:78].strip()
        if elem == "H" or (not elem and atom.startswith("H")): continue
        rn = int(ln[22:26])
        data.setdefault(rn, []).append([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
    return {r: np.array(v) for r, v in data.items()}


def all_residue_numbers(pdb: Path, chain: str) -> list[int]:
    seen = set()
    out = []
    for ln in pdb.read_text().splitlines():
        if not ln.startswith("ATOM"): continue
        if ln[21] != chain: continue
        rn = int(ln[22:26])
        if rn not in seen:
            seen.add(rn); out.append(rn)
    return sorted(out)


def contacts_at_cutoff(binder_heavy, target_heavy, cutoff):
    """Binder residue numbers with any heavy-atom distance <= cutoff to any target heavy atom."""
    t_all = np.vstack(list(target_heavy.values()))
    hit = set()
    c2 = cutoff * cutoff
    for rn, atoms in binder_heavy.items():
        d2 = ((atoms[:, None, :] - t_all[None, :, :]) ** 2).sum(-1)
        if (d2 <= c2).any():
            hit.add(rn)
    return hit


def join_gap1(residues: set[int], all_rn: list[int]) -> set[int]:
    """Add any residue that sits between two design residues separated by exactly 1 fixed slot."""
    out = set(residues)
    idx = {r: i for i, r in enumerate(all_rn)}
    for r in residues:
        i = idx[r]
        # check if r+1 (the next residue in sequence) plus r+2 is also design
        if i + 2 < len(all_rn):
            r1, r2 = all_rn[i + 1], all_rn[i + 2]
            if r2 in residues and r1 not in residues:
                out.add(r1)
    return out


def segments(residues: set[int], all_rn: list[int]) -> list[tuple[int, int]]:
    """Return list of (start_resnum, end_resnum) for each contiguous design segment."""
    idx = {r: i for i, r in enumerate(all_rn)}
    sorted_res = sorted(residues)
    segs = []
    cur_start = cur_end = None
    for r in sorted_res:
        if cur_start is None:
            cur_start = cur_end = r
        elif idx[r] == idx[cur_end] + 1:
            cur_end = r
        else:
            segs.append((cur_start, cur_end))
            cur_start = cur_end = r
    if cur_start is not None:
        segs.append((cur_start, cur_end))
    return segs


def strat4a_range(n: int) -> tuple[int, int]:
    """Strategy 4a: low = round(N*0.7), high = round(N*1.5), N=1 stays 1-1.
    Uses round-half-up to match the convention in pose_solved_5a/params.yml
    (N=7 -> 5-11, where 7*1.5=10.5 rounds up to 11)."""
    if n == 1:
        return 1, 1
    lo = max(1, int(n * 0.7 + 0.5))
    hi = max(lo, int(n * 1.5 + 0.5))
    return lo, hi


def make_contig(design_residues: set[int], all_rn: list[int],
                binder_chain: str, target_chain: str) -> tuple[str, int]:
    """Build RFDiffusion contig (with strategy-4a scaling on de novo lengths).
    Returns (contig_string, total_design_residues_count).
    """
    design_residues = join_gap1(design_residues, all_rn)
    design_segs = segments(design_residues, all_rn)

    # all-fixed segments are runs of non-design residues sandwiched between design ones
    # (or before the first / after the last design segment).
    parts: list[str] = []
    if not design_segs:
        return f"{binder_chain}{all_rn[0]}-{all_rn[-1]} {target_chain}", 0

    first_design_start = design_segs[0][0]
    last_design_end = design_segs[-1][1]

    # leading fixed region
    if first_design_start > all_rn[0]:
        parts.append(f"{binder_chain}{all_rn[0]}-{first_design_start - 1}")

    for i, (ds, de) in enumerate(design_segs):
        n = de - ds + 1
        lo, hi = strat4a_range(n)
        parts.append(f"{lo}-{hi}")
        # trailing fixed (between this design segment and the next, or after last)
        if i + 1 < len(design_segs):
            next_start = design_segs[i + 1][0]
            parts.append(f"{binder_chain}{de + 1}-{next_start - 1}")
        else:
            if last_design_end < all_rn[-1]:
                parts.append(f"{binder_chain}{de + 1}-{all_rn[-1]}")

    contig = "/".join(parts) + f" {target_chain}"
    total = sum(de - ds + 1 for ds, de in design_segs)
    return contig, total


def main():
    if len(sys.argv) < 5:
        sys.exit("usage: derive_contigs.py PDB CHAIN_A CHAIN_B CUTOFF [CUTOFF ...]")
    pdb = Path(sys.argv[1])
    ca, cb = sys.argv[2], sys.argv[3]
    cutoffs = [float(x) for x in sys.argv[4:]]

    binder_heavy = read_heavy(pdb, ca)
    target_heavy = read_heavy(pdb, cb)
    binder_rn = all_residue_numbers(pdb, ca)

    print(f"=== {pdb.name}  binder chain {ca} ({len(binder_rn)} res)  target chain {cb} ===")
    for cut in cutoffs:
        hits = contacts_at_cutoff(binder_heavy, target_heavy, cut)
        contig, n = make_contig(hits, binder_rn, ca, cb)
        joined = sorted(join_gap1(hits, binder_rn))
        raw    = sorted(hits)
        segs = segments(join_gap1(hits, binder_rn), binder_rn)
        seg_lengths = [de - ds + 1 for ds, de in segs]
        print()
        print(f"--- cutoff {cut} A ---")
        print(f"  raw design residues ({len(raw)}): {raw}")
        print(f"  after join-gap=1 ({n} res, {len(segs)} segments, lengths {seg_lengths}): {joined}")
        print(f"  contig: {contig}")


if __name__ == "__main__":
    main()
