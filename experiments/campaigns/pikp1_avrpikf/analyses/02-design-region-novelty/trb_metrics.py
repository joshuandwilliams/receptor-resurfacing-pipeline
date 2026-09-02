"""Rebuild this analysis's geometry from RFdiffusion's own .trb mapping.

WHY THIS EXISTS

`rfdiffusion_metrics.json` is built by bin/rfdiffusion_filter.py, which does not
know how long each design region actually came out. It knows only the total de
novo budget, so it splits that total across the regions in proportion to their
max_len. The guess is exact only when the budget is forced, which is the case for
a contig with a single variable-length region and not otherwise.

`crystal_full_test_contig` has one variable region, so its mapping is right. The
two auto contigs have several, so their mapping is wrong for most designs, and
every residue downstream of the first mis-sized region is shifted by one or two.
That shift is what the reported motif RMSD of 2.6 and 4.0 A was measuring. With
the correct mapping the retained scaffold superposes to about 0.1 A, so nothing
had actually drifted.

RFdiffusion records the truth in the .trb it writes next to each design.
`sampled_mask` holds the realised contig and `con_ref_pdb_idx`/`con_hal_pdb_idx`
hold the exact input-to-output residue correspondence. The pipeline discards the
.trb, so these are recovered from the Nextflow work directory.

This module is a local fix for this analysis only. bin/rfdiffusion_filter.py
still has the bug and carries a note pointing here.

Set TRB_ROOT and PDB_ROOT to local copies laid out as
    $PDB_ROOT/inputs/<input>.pdb
    $PDB_ROOT/runs/<run>/results/rfdiffusion/design_*.pdb
    $TRB_ROOT/<run>/**/design_*.trb
"""
from __future__ import annotations

import json
import os
import re
from pathlib import Path

import numpy as np

import motif_rigidity as mr
import novelty_common as nc

TRB_ROOT = Path(os.environ.get("TRB_ROOT", "")).expanduser()
OUT_SUFFIX = ".trb_corrected.json"


def trb_paths(run: str) -> dict[str, Path]:
    """design name -> .trb path.

    A retried Nextflow task can leave a superseded execution behind with the
    same design names. Those are resolved upstream by checksumming the work
    directory's PDB against the published one, so anything left here is
    expected to be unique.
    """
    out: dict[str, Path] = {}
    for p in sorted((TRB_ROOT / run).rglob("*.trb")):
        name = p.stem + ".pdb"
        if name in out:
            raise RuntimeError(
                f"{run}: two .trb files for {name}. Resolve which execution was "
                f"published before rebuilding.")
        out[name] = p
    return out


def realised_lengths(sampled_mask: str) -> list[int]:
    """Design region lengths RFdiffusion actually sampled, from sampled_mask.

    The realised mask writes every segment with equal bounds, so `4-4` is a
    four-residue design region and `A44-67` is retained scaffold.
    """
    return [int(seg.split("-")[0])
            for seg in sampled_mask.split("/") if not seg.startswith("A")
            and re.fullmatch(r"\d+-\d+", seg)]


def rebuild(run: str, label: str) -> dict:
    """Recompute the geometry this analysis reads, using the .trb mapping."""
    inp = mr.ca_coords(mr.PDB_ROOT / "inputs" / mr.INPUT_PDB)
    meta = json.loads((nc.DATA / f"{run}.json").read_text())
    paths = trb_paths(run)
    pdb_dir = mr.PDB_ROOT / "runs" / run / "results" / "rfdiffusion"
    blocks = mr.fixed_blocks(nc.CONTIGS[label])

    designs = []
    for des in meta["designs"]:
        name = des["design"]
        if name not in paths or not (pdb_dir / name).exists():
            continue
        trb = np.load(paths[name], allow_pickle=True)
        dc = mr.ca_coords(pdb_dir / name)

        # Authoritative fixed-residue correspondence, receptor chain only.
        pairs = [(int(h[1]), int(r[1]))
                 for r, h in zip(trb["con_ref_pdb_idx"], trb["con_hal_pdb_idx"])
                 if r[0] == "A" and int(h[1]) in dc and int(r[1]) in inp]
        if len(pairs) < 3:
            continue
        P = np.array([dc[h] for h, _ in pairs])
        Q = np.array([inp[i] for _, i in pairs])
        motif_rmsd = mr._kabsch_rmsd(P, Q)

        # Put the design into the input frame, as the pipeline does, so region
        # displacements are comparable with the originals.
        Pc, Qc = P - P.mean(axis=0), Q - Q.mean(axis=0)
        U, _, Vt = np.linalg.svd(Pc.T @ Qc)
        d = np.sign(np.linalg.det(Vt.T @ U.T))
        R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
        aligned = {rn: (c - P.mean(axis=0)) @ R.T + Q.mean(axis=0)
                   for rn, c in dc.items()}

        lens = realised_lengths(trb["sampled_mask"][0])
        fixed_hal = sorted(h for h, _ in pairs)

        # Walk the realised mask to find each design region's output residues.
        region_coords, region_lengths = [], []
        cursor, region = 1, 0
        for seg in trb["sampled_mask"][0].split("/"):
            if seg.startswith("A"):
                lo, hi = (int(v) for v in seg[1:].split("-"))
                cursor += hi - lo + 1
            elif re.fullmatch(r"\d+-\d+", seg):
                n = lens[region]
                coords = [aligned[cursor + k] for k in range(n)
                          if (cursor + k) in aligned]
                region_coords.append([c.tolist() for c in coords])
                region_lengths.append(n)
                cursor += n
                region += 1

        # Input residues each region replaced, taken between the flanking
        # retained blocks in input numbering.
        com_disp, ep_design, ep_input = [], [], []
        gaps = []
        prev_hi = None
        for k, (lo, hi) in enumerate(blocks):
            if prev_hi is not None:
                gaps.append((prev_hi, lo))
            prev_hi = hi
        # A trailing design region has a left anchor only.
        anchors = list(gaps) + [(prev_hi, None)] if len(region_coords) > len(gaps) \
            else list(gaps)

        for k, coords in enumerate(region_coords):
            arr = np.array(coords) if coords else np.zeros((0, 3))
            ep_design.append(float(np.linalg.norm(arr[-1] - arr[0]))
                             if len(arr) >= 2 else float("nan"))
            if k >= len(anchors) or not len(arr):
                com_disp.append(float("nan"))
                ep_input.append(float("nan"))
                continue
            left, right = anchors[k]
            if right is not None:
                gap = [inp[r] for r in range(left + 1, right) if r in inp]
            else:
                gap = [inp[r] for r in range(left + 1, left + 1 + len(arr))
                       if r in inp]
            if gap:
                g = np.array(gap)
                com_disp.append(float(np.linalg.norm(arr.mean(axis=0) - g.mean(axis=0))))
                ep_input.append(float(np.linalg.norm(g[-1] - g[0]))
                                if right is not None and len(g) >= 2 else float("nan"))
            else:
                com_disp.append(float("nan"))
                ep_input.append(float("nan"))

        designs.append({
            "design": name,
            "motif_rmsd": round(motif_rmsd, 4),
            "reported_motif_rmsd": des.get("motif_rmsd"),
            "per_region_lengths": region_lengths,
            "reported_per_region_lengths": des.get("per_region_lengths"),
            "lengths_agree": region_lengths == des.get("per_region_lengths"),
            "per_region_coords": region_coords,
            "design_region_com_displacement": com_disp,
            "endpoint_distance_design": ep_design,
            "endpoint_distance_input": ep_input,
            "passes_filter": bool(des.get("passes_filter")),
            "sampled_mask": trb["sampled_mask"][0],
            # Not recomputed. These depend on contact classification, which is
            # also downstream of the bad mapping, so the originals are wrong for
            # the auto runs and are carried only to keep the schema.
            "design_region_coverage": None,
            "frac_contacts_in_design": None,
        })

    return {"run": run, "label": label, "n_designs": len(designs),
            "source": "rebuilt from RFdiffusion .trb", "designs": designs}


def main():
    for label, run in nc.RUNS.items():
        out = rebuild(run, label)
        path = nc.DATA / f"{run}{OUT_SUFFIX}"
        path.write_text(json.dumps(out))
        bad = sum(1 for d in out["designs"] if not d["lengths_agree"])
        rms = np.median([d["motif_rmsd"] for d in out["designs"]])
        old = np.median([d["reported_motif_rmsd"] for d in out["designs"]
                         if d["reported_motif_rmsd"] is not None])
        print(f"{label}: {out['n_designs']} designs, "
              f"{bad} with wrong per_region_lengths, "
              f"motif RMSD {old:.3f} -> {rms:.3f} A")


if __name__ == "__main__":
    main()
