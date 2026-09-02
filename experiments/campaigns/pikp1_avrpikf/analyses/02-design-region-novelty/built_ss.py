"""Secondary structure of what RFdiffusion actually built, by length.

The length-diversity curve has a minimum near the wild-type gap. One reading is
that the residue budget near 13-14 is what lets the native strand-turn-strand
element be rebuilt, and that both shorter and longer budgets force something
non-native. That predicts strand content in the built region to peak in the
middle of the allowance and fall off on both sides.

RFdiffusion writes N, CA, C and O, so DSSP is computable directly. This is the
Kabsch and Sander electrostatic hydrogen-bond definition, three-state only.
Validated against the HELIX and SHEET records of the input complex.
"""
from __future__ import annotations

import re

import numpy as np

BACKBONE = (" N  ", " CA ", " C  ", " O  ")
Q1Q2F = 0.42 * 0.20 * 332.0
HB_CUTOFF = -0.5


def backbone_coords(pdb, chain: str = "A"):
    """(resnums, dict of N/CA/C/O arrays) for one chain, ordered by residue."""
    per: dict[int, dict[str, np.ndarray]] = {}
    for line in pdb.read_text().splitlines():
        if not line.startswith("ATOM") or line[21] != chain:
            continue
        name = line[12:16]
        if name not in BACKBONE:
            continue
        rn = int(line[22:26])
        per.setdefault(rn, {})[name.strip()] = np.array(
            [float(line[30:38]), float(line[38:46]), float(line[46:54])])
    resnums = sorted(r for r in per if len(per[r]) == 4)
    return resnums, {a: np.array([per[r][a] for r in resnums])
                     for a in ("N", "CA", "C", "O")}


def hbond_map(bb: dict[str, np.ndarray]) -> np.ndarray:
    """Boolean [donor, acceptor]. Donor i is the NH, acceptor j the C=O."""
    N, C, O = bb["N"], bb["C"], bb["O"]
    n = len(N)

    # Amide H sits one angstrom from N, along the previous residue's C-to-O
    # direction. Residue 0 has no preceding carbonyl and cannot donate.
    H = np.full_like(N, np.nan)
    v = C[:-1] - O[:-1]
    H[1:] = N[1:] + v / np.linalg.norm(v, axis=1, keepdims=True)

    def d(a, b):
        return np.linalg.norm(a[:, None, :] - b[None, :, :], axis=-1) + 1e-9

    E = Q1Q2F * (1.0 / d(N, O) + 1.0 / d(H, C) - 1.0 / d(H, O) - 1.0 / d(N, C))
    hb = E < HB_CUTOFF
    hb[0, :] = False                      # no H on the first residue
    idx = np.arange(n)
    hb[np.abs(idx[:, None] - idx[None, :]) < 2] = False
    return hb


def assign(bb: dict[str, np.ndarray]) -> str:
    """Three-state string of 'H', 'E' and 'C', one character per residue."""
    hb = hbond_map(bb)
    n = len(bb["N"])
    ss = np.array(["C"] * n)

    # Strands first, so helices overwrite them where both are called.
    for i in range(1, n - 1):
        for j in range(i + 3, n - 1):
            anti = ((hb[i, j] and hb[j, i])
                    or (hb[i - 1, j + 1] and hb[j - 1, i + 1]))
            para = ((hb[i, j - 1] and hb[j + 1, i])
                    or (hb[j, i - 1] and hb[i + 1, j]))
            if anti or para:
                ss[i] = ss[j] = "E"

    # A 4-turn at i and at i-1 makes residues i..i+3 helical.
    turn = np.zeros(n, bool)
    for i in range(n - 4):
        turn[i] = hb[i + 4, i]
    for i in range(1, n - 4):
        if turn[i] and turn[i - 1]:
            ss[i:i + 4] = "H"

    return "".join(ss)


def spans(ss: str, code: str) -> list[tuple[int, int]]:
    """Contiguous 1-based runs of one code, for printing next to PDB records."""
    return [(m.start() + 1, m.end()) for m in re.finditer(f"{code}+", ss)]


def demote_isolated_bridges(ss: str) -> str:
    """Recode single-residue 'E' runs as 'B', following DSSP.

    A residue can satisfy the bridge test on its own, but DSSP only calls a
    strand when two consecutive residues do. A lone bridge is 'B'. Counting
    those as strand segments splits one loop into two and inflates the
    apparent number of elements.
    """
    return re.sub(r"(?<!E)E(?!E)", "B", ss)


def strand_segments(ss: str) -> int:
    """Number of true strand segments, isolated bridges excluded."""
    return len(re.findall(r"E{2,}", demote_isolated_bridges(ss)))
