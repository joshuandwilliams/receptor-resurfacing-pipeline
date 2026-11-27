#!/usr/bin/env python3
"""Does the AVR-Pia-specific contact test survive Chapter 4's 8 A Ca cutoff?

interface_extract.py scores contacts as closest heavy atom within 5 A, because
that is how the pipeline reports design contacts and the pose contact sets are
read straight out of its tables. The benchmarking chapter instead used Ca-Ca
within 8 A. This rebuilds the whole test under both criteria from the
prediction PDBs, so the only thing that changes between the two runs is the
contact definition.

Both the reference surfaces and the pose contacts are recomputed here. Reading
the pose contacts from the pipeline CSV would fix them at 5 A heavy atom and
make the comparison meaningless.

Usage:
    interface_cutoff_comparison.py [--output interface_cutoff_comparison.csv]
"""
import argparse
import os
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import protein_letters_3to1

from pipeline_truth import never_steered

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
REPO = os.path.abspath(os.path.join(HERE, "..", "..", "..", "..", ".."))
HPC_ROOT = "/hpc-home/jowillia/receptor_design/receptor-resurfacing-pipeline"
RUNS = os.path.join(HERE, "..", "..", "runs", "crystal_full_test_contig",
                    "results", "negative_steering", "runs")

NEVER_STEERED = never_steered()
RA_EFF_THRESHOLD = 5.0
GAP_OPEN, GAP_EXTEND = -10.0, -0.5
DR1 = set(range(33, 46))
DR2 = set(range(73, 79))

MODES = {"heavy5": ("heavy", 5.0), "ca8": ("ca", 8.0)}


def aligner():
    al = Align.PairwiseAligner()
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score = GAP_OPEN
    al.extend_gap_score = GAP_EXTEND
    al.mode = "global"
    return al


def read_chain(path, chain_id):
    """(sequence, heavy-atom coords per residue, Ca coord per residue)."""
    model = PDBParser(QUIET=True).get_structure("s", path)[0]
    seq, heavy, ca = [], [], []
    for res in model[chain_id]:
        if res.id[0] != " ":
            continue
        try:
            aa = protein_letters_3to1[res.get_resname()]
        except KeyError:
            continue
        h = [a.coord for a in res if a.element != "H"]
        if not h or "CA" not in res:
            continue
        seq.append(aa)
        heavy.append(np.asarray(h, dtype=float))
        ca.append(np.asarray(res["CA"].coord, dtype=float))
    return "".join(seq), heavy, np.asarray(ca, dtype=float)


def contact_positions(rec, eff, mode, cutoff):
    """1-based receptor positions in contact with the effector."""
    if mode == "ca":
        d = np.linalg.norm(rec[2][:, None, :] - eff[2][None, :, :], axis=-1)
        return {int(i) + 1 for i in np.where(d.min(axis=1) <= cutoff)[0]}
    eff_atoms = np.vstack(eff[1])
    out = set()
    for i, a in enumerate(rec[1]):
        d = np.linalg.norm(a[:, None, :] - eff_atoms[None, :, :], axis=-1)
        if d.min() <= cutoff:
            out.add(i + 1)
    return out


def map_to_frame(query_seq, frame_seq, al):
    a = al.align(frame_seq, query_seq)[0]
    out = {}
    for (fs, fe), (qs, qe) in zip(*a.aligned):
        for k in range(fe - fs):
            out[qs + k + 1] = fs + k + 1
    return out


def local(path):
    return path.replace(HPC_ROOT, REPO) if isinstance(path, str) else path


def cold_pdb(unit, per_seed):
    """Median-ra_eff cold-start seed for a unit, mirroring the pipeline."""
    g = per_seed[(per_seed.unit == unit) & (per_seed.design == "initial")]
    g = g.dropna(subset=["ra_eff"]).sort_values("ra_eff")
    if not len(g):
        return None
    row = g.iloc[(len(g) - 1) // 2]
    idx = int(row.seed_index) if not pd.isna(row.seed_index) else 0
    name = "initial_prediction.pdb" if idx == 0 else f"initial_prediction_s{idx}.pdb"
    p = os.path.join(RUNS, unit, "cycle_0", name)
    return p if os.path.isfile(p) else None


def build(mode_name, al, surv, per_seed):
    mode, cutoff = MODES[mode_name]
    frame = read_chain(os.path.join(DATA, "pikp1_avrpikf_reference.pdb"), "A")
    pikf = read_chain(os.path.join(DATA, "pikp1_avrpikf_reference.pdb"), "B")
    pik_site = contact_positions(frame, pikf, mode, cutoff)

    pia_rec = read_chain(os.path.join(DATA, "6Q76.pdb"), "A")
    pia_eff = read_chain(os.path.join(DATA, "6Q76.pdb"), "B")
    pia_map = map_to_frame(pia_rec[0], frame[0], al)
    pia_site = {pia_map[p] for p in contact_positions(pia_rec, pia_eff, mode, cutoff)
                if p in pia_map}

    shared = pik_site & pia_site
    pik_only, pia_only = pik_site - shared, pia_site - shared
    pia_spec = pia_only - (DR1 | DR2)

    rows = []
    for r in surv.itertuples():
        if not isinstance(r.corrected_receptor, str):
            continue
        rec_pdb = local(r.rep_canonical_pdb)
        cpdb = cold_pdb(r.mpnn_sequence, per_seed)
        row = {"sequence": r.mpnn_sequence,
               "steered_correct": r.rep_ra_eff_vs_truth_median <= RA_EFF_THRESHOLD,
               "cold_correct": (r.mpnn_sequence in NEVER_STEERED
                                and r.rep_ra_eff_vs_truth_median <= RA_EFF_THRESHOLD)}
        for tag, path in (("cold", cpdb), ("steered", rec_pdb)):
            if not path or not os.path.isfile(path):
                row[f"{tag}_n_pia_specific"] = np.nan
                row[f"{tag}_touches"] = np.nan
                continue
            rec = read_chain(path, "A")
            eff = read_chain(path, "B")
            m = map_to_frame(rec[0], frame[0], al)
            hit = {m[p] for p in contact_positions(rec, eff, mode, cutoff) if p in m}
            row[f"{tag}_n_pia_specific"] = len(hit & pia_spec)
            row[f"{tag}_touches"] = len(hit & pia_spec) > 0
        rows.append(row)
    return pd.DataFrame(rows), dict(
        pik_only=len(pik_only), shared=len(shared), pia_only=len(pia_only),
        pia_specific=sorted(pia_spec))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", default=os.path.join(
        HERE, "interface_cutoff_comparison.csv"))
    args = ap.parse_args()
    al = aligner()
    surv = pd.read_csv(os.path.join(DATA, "crystal_full_test_contig.survivors.csv"))
    per_seed = pd.read_csv(os.path.join(DATA, "campaign_per_seed.csv"))

    frames = {}
    for name in MODES:
        df, sets = build(name, al, surv, per_seed)
        frames[name] = df
        print(f"\n=== {name} ===")
        print(f"  AVR-Pik only {sets['pik_only']} | shared {sets['shared']} "
              f"| AVR-Pia only {sets['pia_only']}")
        print(f"  AVR-Pia-specific set (outside design regions): {sets['pia_specific']}")
        for stage in ("cold", "steered"):
            sub = df.dropna(subset=[f"{stage}_touches"])
            ok = sub[sub[f"{stage}_correct"].astype(bool)]
            bad = sub[~sub[f"{stage}_correct"].astype(bool)]
            print(f"  {stage}: correct {int(ok[f'{stage}_touches'].sum())}/{len(ok)} "
                  f"touch | incorrect {int(bad[f'{stage}_touches'].sum())}/{len(bad)} touch")

    merged = frames["heavy5"].merge(frames["ca8"], on="sequence",
                                    suffixes=("_h5", "_ca8"))
    merged.to_csv(args.output, index=False)
    print(f"\nwrote {args.output} ({len(merged)} designs)")


if __name__ == "__main__":
    main()
