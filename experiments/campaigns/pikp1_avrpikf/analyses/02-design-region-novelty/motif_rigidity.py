"""Is the retained scaffold internally distorted, or just repositioned?

`motif_rmsd` is the residual after Kabsch on every fixed residue at once, so a
large value has two possible causes that the single number cannot separate.
Either each fixed block is internally deformed, or each block is still rigid but
the blocks have moved relative to one another.

The decomposition is direct. Align one fixed block onto its input counterpart on
its own and the residual is that block's internal deformation. Compare that with
the all-blocks-at-once residual, and the gap is inter-block rigid-body motion.

Needs the RFdiffusion output PDBs, which live on the HPC rather than in the repo.
Point PDB_ROOT at a local copy, laid out as
    <root>/inputs/<input>.pdb
    <root>/runs/<run>/results/rfdiffusion/design_*.pdb
"""
from __future__ import annotations

import json
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd

import novelty_common as nc

PDB_ROOT = Path(os.environ.get("PDB_ROOT", "")).expanduser()
INPUT_PDB = "pikp1_avrpikf_pik_interface_complex.pdb"


def ca_coords(pdb: Path, chain: str = "A") -> dict[int, np.ndarray]:
    """resnum -> Ca coordinate for one chain."""
    out = {}
    for line in pdb.read_text().splitlines():
        if not line.startswith("ATOM") or line[12:16] != " CA ":
            continue
        if line[21] != chain:
            continue
        out[int(line[22:26])] = np.array(
            [float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return out


def fixed_blocks(contig: str) -> list[tuple[int, int]]:
    """The A<i>-<j> segments of a contig, in order."""
    return [(int(a), int(b)) for a, b in
            re.findall(r"A(\d+)-(\d+)", contig.split()[0])]


def design_to_input(contig: str, region_lengths: list[int]) -> dict[int, int]:
    """Design residue number -> input residue number, for fixed residues only.

    RFdiffusion renumbers its output 1..N contiguously, so the mapping is
    recovered by walking the contig and consuming the design regions at the
    lengths this particular design actually used.
    """
    mapping, d_idx, region = {}, 1, 0
    for seg in contig.split()[0].split("/"):
        if seg.startswith("A"):
            lo, hi = (int(v) for v in seg[1:].split("-"))
            for rn in range(lo, hi + 1):
                mapping[d_idx] = rn
                d_idx += 1
        else:
            d_idx += region_lengths[region]
            region += 1
    return mapping


def _kabsch_rmsd(P: np.ndarray, Q: np.ndarray) -> float:
    """RMSD of P onto Q after optimal superposition, as the pipeline does it."""
    Pc, Qc = P - P.mean(axis=0), Q - Q.mean(axis=0)
    U, _, Vt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return float(np.sqrt(np.mean(np.sum((Pc @ R.T - Qc) ** 2, axis=1))))


def build() -> pd.DataFrame:
    """One row per (run, design) with the global and per-block residuals."""
    inp = ca_coords(PDB_ROOT / "inputs" / INPUT_PDB)
    rows = []

    for label, run in nc.RUNS.items():
        contig = nc.CONTIGS[label]
        blocks = fixed_blocks(contig)
        meta = json.loads((nc.DATA / f"{run}.json").read_text())
        pdb_dir = PDB_ROOT / "runs" / run / "results" / "rfdiffusion"

        for des in meta["designs"]:
            pdb = pdb_dir / des["design"]
            if not pdb.exists():
                continue
            dc = ca_coords(pdb)
            m = design_to_input(contig, des["per_region_lengths"])

            pairs = [(dc[d], inp[i]) for d, i in m.items()
                     if d in dc and i in inp]
            if len(pairs) < 3:
                continue
            P = np.array([p for p, _ in pairs])
            Q = np.array([q for _, q in pairs])
            global_rmsd = _kabsch_rmsd(P, Q)

            # Same residues, but each block aligned on its own.
            per_block = []
            for lo, hi in blocks:
                bp = [(dc[d], inp[i]) for d, i in m.items()
                      if lo <= i <= hi and d in dc and i in inp]
                if len(bp) < 3:
                    per_block.append(np.nan)
                    continue
                per_block.append(_kabsch_rmsd(
                    np.array([p for p, _ in bp]), np.array([q for _, q in bp])))

            rows.append({
                "run": label, "design": des["design"],
                "n_blocks": len(blocks),
                "reported_motif_rmsd": nc._num(des.get("motif_rmsd")),
                "global_motif_rmsd": global_rmsd,
                "mean_within_block_rmsd": float(np.nanmean(per_block)),
                "max_within_block_rmsd": float(np.nanmax(per_block)),
                **{f"block{k + 1}_rmsd": v for k, v in enumerate(per_block)},
            })

    return pd.DataFrame(rows)
