#!/usr/bin/env python3
"""Where the six design region 1 differences between d39_s0 and d39_s1 point.

The two sequences sit on the same RFdiffusion backbone, so every difference is
a side chain swap in a fixed loop. The question is which of the six face the
effector and which face back onto the HMA body, because only the second kind
can plausibly perturb the domain itself.

For each position, in each structure, this reports:

  eff       closest heavy-atom approach from the side chain (CB outwards) to
            AVR-PikF, chain B
  body      closest approach to the HMA body, taken as chain A outside the
            design region, so contacts with neighbours inside the rebuilt loop
            do not count as packing against the domain
  rel SASA  side-chain solvent accessibility relative to a Gly-X-Gly reference,
            computed on the receptor alone so the effector does not bury it
  partner   nearest charged or polar group on the HMA body, for the buried
            charge check

A buried charge with no counter-charge is the thing worth finding: it is the
classic way a loop substitution destabilises a small domain.
"""
from pathlib import Path

from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

HERE = Path(__file__).parent
DATA = HERE / "data"

S0 = DATA / "d39_s0_coldstart_s2.pdb"
S1 = DATA / "d39_s1_coldstart_s1.pdb"
WT = DATA / "wt_pikp1_avrpikf.pdb"

DR1 = range(33, 49)          # design region 1 in the two designs
DR1_WT = range(33, 46)       # shorter in the wild type
POSITIONS = [36, 39, 40, 44, 45, 46]

# Tien et al. 2013, theoretical maxima from Gly-X-Gly.
MAXA = {"ALA": 129, "ARG": 274, "ASN": 195, "ASP": 193, "CYS": 167, "GLN": 225,
        "GLU": 223, "GLY": 104, "HIS": 224, "ILE": 197, "LEU": 201, "LYS": 236,
        "MET": 224, "PHE": 240, "PRO": 159, "SER": 155, "THR": 172, "TRP": 285,
        "TYR": 263, "VAL": 174}
CHARGED = {"ARG": "+", "LYS": "+", "ASP": "-", "GLU": "-", "HIS": "+/0"}
POLAR_AT = {"N", "O", "S"}
BACKBONE = {"N", "CA", "C", "O", "OXT"}


def load(path):
    return PDBParser(QUIET=True).get_structure(path.stem, str(path))[0]


def sidechain(res):
    """Heavy atoms past the backbone. Gly falls back to CA."""
    at = [a for a in res if a.element != "H" and a.get_name() not in BACKBONE]
    return at or [a for a in res if a.get_name() == "CA"]


def min_dist(atoms, others):
    if not atoms or not others:
        return float("inf")
    return min((a - b) for a in atoms for b in others)


def rel_sasa(model, chain_id, resid):
    """Side-chain SASA of one residue, receptor only, relative to Gly-X-Gly."""
    import copy
    m = copy.deepcopy(model)
    for ch in list(m):
        if ch.id != chain_id:
            m.detach_child(ch.id)
    ShrakeRupley().compute(m, level="A")
    res = m[chain_id][resid]
    tot = sum(a.sasa for a in res if a.get_name() not in BACKBONE and a.element != "H")
    return 100.0 * tot / MAXA.get(res.get_resname(), 200)


def report(path, dr1, label):
    model = load(path)
    A = model["A"]
    B = model["B"] if "B" in model else None
    body = [a for r in A for a in r
            if r.id[1] not in dr1 and a.element != "H"]
    eff = [a for r in B for a in r if a.element != "H"] if B else []

    print(f"\n=== {label} ===")
    print(f"{'res':>5} {'aa':>4} {'eff A':>7} {'body A':>7} {'relSASA':>8}  note")
    for pos in POSITIONS:
        if pos not in [r.id[1] for r in A]:
            continue
        res = A[pos]
        sc = sidechain(res)
        d_eff = min_dist(sc, eff)
        d_body = min_dist(sc, body)
        sas = rel_sasa(model, "A", pos)
        name = res.get_resname()

        notes = []
        if d_eff <= 4.5:
            notes.append("contacts effector")
        if d_body <= 4.0:
            notes.append("packs on HMA body")
        if sas < 20:
            notes.append("buried")
        elif sas < 40:
            notes.append("partly buried")
        if name in CHARGED and sas < 30:
            polar = [a for a in body if a.element in POLAR_AT]
            near = min_dist([a for a in sc if a.element in POLAR_AT] or sc, polar)
            notes.append(f"BURIED CHARGE, nearest body polar {near:.2f} A")
        print(f"{pos:>5} {name:>4} {d_eff:>7.2f} {d_body:>7.2f} {sas:>7.1f}%  "
              f"{'; '.join(notes) if notes else 'solvent facing'}")


if __name__ == "__main__":
    report(S0, DR1, "d39_s0  (cold start, not autoactive)")
    report(S1, DR1, "d39_s1  (cold start, autoactive)")
