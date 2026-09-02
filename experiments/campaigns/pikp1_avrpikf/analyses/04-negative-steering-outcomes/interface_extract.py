#!/usr/bin/env python3
"""Do de novo designs mis-pose onto the AVR-Pia surface, before and after steering?

On natural Pik HMAs the benchmarking chapter found that Boltz-2 mis-poses AVR-Pik
effectors onto the AVR-Pia surface. Its test was not which surface a pose
overlaps more, since the AVR-Pik surface is the larger of the two and the shared
positions are contacted either way. The test was whether the pose contacts any
AVR-Pia-SPECIFIC residue at all. Poses within 8 A of the correct position
contacted none, and every pose beyond 21 A contacted some.

This repeats that test on the 128 resurfaced receptors of the AVR-PikF campaign,
for the cold start and for the representative steered design, so the effect of
steering on AVR-Pia-specific contact can be read directly.

Three disjoint sets are built on a common frame, as in the benchmarking chapter:

* frame          chain A of the campaign reference complex, 78 residues, 1-78
* AVR-Pik site   frame residues contacting AVR-PikF, closest heavy atom < 5 A
* AVR-Pia site   6Q76 chain A residues contacting AVR-Pia on the same criterion,
                 mapped on by global BLOSUM62 alignment
* shared         the intersection, subtracted from both to give the -only sets

Three of the AVR-Pia-only positions fall inside design region 1, which sits
between the two surfaces. A correctly posed effector contacts them anyway,
because alternating side chains on the same strand face opposite ways, so the
AVR-Pia-only set is further restricted to positions outside both design
regions. Without that restriction the test cannot separate correct from
incorrect poses on this system at all.

Contact sets come from the pipeline's own tables. `steered_contact_residues` on
the `initial` row of each run gives the cold start, and
`rep_contact_residues_majority` in the survivor table gives the representative
steered design. Both are per-chain 0-based indices, shifted to 1-based here.

Usage:
    interface_extract.py [--output campaign_interface_check.csv]
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

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")

# The pipeline reports design contacts as closest heavy atom within 5 A, so the
# two reference surfaces are defined the same way. Using the benchmarking
# chapter's 8 A Ca criterion here instead makes correct poses appear to touch
# AVR-Pia-specific residues, because the two definitions disagree rather than
# because the effector moved.
CONTACT_CUTOFF = 5.0
from pipeline_truth import never_steered

NEVER_STEERED = never_steered()
RA_EFF_THRESHOLD = 5.0

GAP_OPEN = -10.0
GAP_EXTEND = -0.5

# Where the two design regions sit on the frame, used only to report how much
# of each reference surface the designer was free to overwrite.
DR1 = set(range(33, 46))
DR2 = set(range(73, 79))


def aligner():
    al = Align.PairwiseAligner()
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score = GAP_OPEN
    al.extend_gap_score = GAP_EXTEND
    al.mode = "global"
    return al


def read_chain(path, chain_id):
    """Ordered (sequence, per-residue heavy-atom coords) for one chain."""
    model = PDBParser(QUIET=True).get_structure("s", path)[0]
    seq, atoms = [], []
    for res in model[chain_id]:
        if res.id[0] != " ":
            continue
        try:
            aa = protein_letters_3to1[res.get_resname()]
        except KeyError:
            continue
        coords = [a.coord for a in res if a.element != "H"]
        if not coords:
            continue
        seq.append(aa)
        atoms.append(np.asarray(coords, dtype=float))
    return "".join(seq), atoms


def contact_positions(rec_atoms, eff_atoms, cutoff=CONTACT_CUTOFF):
    """1-based receptor positions whose closest heavy atom is within cutoff."""
    eff = np.vstack(eff_atoms)
    out = set()
    for i, a in enumerate(rec_atoms):
        d = np.linalg.norm(a[:, None, :] - eff[None, :, :], axis=-1)
        if d.min() <= cutoff:
            out.add(i + 1)
    return out


def map_to_frame(query_seq, frame_seq, al):
    """query 1-based position -> frame 1-based position, for aligned pairs."""
    a = al.align(frame_seq, query_seq)[0]
    out = {}
    for (fs, fe), (qs, qe) in zip(*a.aligned):
        for k in range(fe - fs):
            out[qs + k + 1] = fs + k + 1
    return out


def parse_contacts(cell):
    """Pipeline contact lists are 0-based; return 1-based positions."""
    if not isinstance(cell, str) or not cell.strip():
        return set()
    return {int(tok) + 1 for tok in cell.split(",") if tok.strip()}


def describe(contacts, pik_only, shared, pia_only, pia_undesigned):
    def frac(s):
        return len(contacts & s) / len(s) if s else np.nan
    return {
        "n_contacts": len(contacts),
        "n_pia_specific": len(contacts & pia_undesigned),
        "touches_pia_specific": len(contacts & pia_undesigned) > 0,
        "frac_pik_only": frac(pik_only),
        "frac_shared": frac(shared),
        "frac_pia_only": frac(pia_only),
        "frac_pia_specific": frac(pia_undesigned),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", default=os.path.join(
        HERE, "campaign_interface_check.csv"))
    args = ap.parse_args()

    al = aligner()

    frame_seq, frame_atoms = read_chain(
        os.path.join(DATA, "pikp1_avrpikf_reference.pdb"), "A")
    _, pikf_atoms = read_chain(
        os.path.join(DATA, "pikp1_avrpikf_reference.pdb"), "B")
    pik_site = contact_positions(frame_atoms, pikf_atoms)

    pia_seq, pia_atoms = read_chain(os.path.join(DATA, "6Q76.pdb"), "A")
    _, avrpia_atoms = read_chain(os.path.join(DATA, "6Q76.pdb"), "B")
    pia_map = map_to_frame(pia_seq, frame_seq, al)
    pia_site = {pia_map[p] for p in contact_positions(pia_atoms, avrpia_atoms)
                if p in pia_map}

    shared = pik_site & pia_site
    pik_only = pik_site - shared
    pia_only = pia_site - shared
    print(f"frame {len(frame_seq)} residues | AVR-Pik only {len(pik_only)} "
          f"| shared {len(shared)} | AVR-Pia only {len(pia_only)}")
    pia_undesigned = pia_only - (DR1 | DR2)
    print(f"AVR-Pia-only positions: {sorted(pia_only)}")
    print(f"  {sorted(pia_only & (DR1 | DR2))} lie inside a design region "
          f"and are excluded")
    print(f"  AVR-Pia-specific set used for the test: {sorted(pia_undesigned)}")

    surv = pd.read_csv(
        os.path.join(DATA, "crystal_full_test_contig.survivors.csv"))
    cold = pd.read_csv(
        os.path.join(DATA, "crystal_full_test_contig.cold_start_contacts.csv")
    ).set_index("sequence")

    rows = []
    for r in surv.itertuples():
        # The two sequence controls carry no corrected receptor.
        if not isinstance(r.corrected_receptor, str):
            continue
        m = map_to_frame(r.corrected_receptor, frame_seq, al)

        def to_frame(cell):
            return {m[p] for p in parse_contacts(cell) if p in m}

        row = {"sequence": r.mpnn_sequence,
               "cross_tier": r.cross_tier,
               "steered_ra_eff": r.rep_ra_eff_vs_truth_median,
               "steered_correct": (r.rep_ra_eff_vs_truth_median
                                   <= RA_EFF_THRESHOLD),
               "was_steered": r.rep_design != "initial"}

        c = cold.loc[r.mpnn_sequence] if r.mpnn_sequence in cold.index else None
        if c is not None:
            row["cold_ra_eff"] = c.cold_ra_eff
            # Already correctly posed before steering: the pipeline left it
            # unsteered (cold_start_all_clean) AND its representative clears
            # the 5 A pose threshold. The two differ for one design, whose
            # cold seeds all clear the 6 A skip bar but not the 5 A pass bar.
            row["cold_correct"] = (r.mpnn_sequence in NEVER_STEERED
                                   and row["steered_correct"])
            for k, v in describe(to_frame(c.cold_contact_residues), pik_only,
                                 shared, pia_only, pia_undesigned).items():
                row[f"cold_{k}"] = v
        for k, v in describe(to_frame(r.rep_contact_residues_majority),
                             pik_only, shared, pia_only,
                             pia_undesigned).items():
            row[f"steered_{k}"] = v
        rows.append(row)

    out = pd.DataFrame(rows)
    out.to_csv(args.output, index=False)
    print(f"\nwrote {args.output} ({len(out)} designs)")

    for stage in ("cold", "steered"):
        sub = out.dropna(subset=[f"{stage}_touches_pia_specific"])
        ok = sub[sub[f"{stage}_correct"]]
        bad = sub[~sub[f"{stage}_correct"].astype(bool)]
        print(f"\n{stage}: {len(ok)} correct, {len(bad)} incorrect")
        print(f"  correct   touching an AVR-Pia-specific residue: "
              f"{int(ok[f'{stage}_touches_pia_specific'].sum())} of {len(ok)}")
        print(f"  incorrect touching an AVR-Pia-specific residue: "
              f"{int(bad[f'{stage}_touches_pia_specific'].sum())} of {len(bad)}")


if __name__ == "__main__":
    main()
