#!/usr/bin/env python3
"""
pose_solver.py
--------------
Rigid-body placement of target protein relative to binder,
driven by pair CA-CA contact constraints.

Objective (strict priority):
  a) Pair validity  : pair residues must not be inside the opposite protein.
                      Violation → astronomically large penalty (pose is invalid).
  b) Pair distances : all pair CA-CA distances must be < max_pair_distance.
                      Violation → very large quadratic penalty.
  c) Interpenetration: minimise depth of ALL CA atoms inside the opposite hull.
                      Small coefficient; only meaningful once a+b are satisfied.

Surface detection uses the convex hull of each protein's CA positions.

Usage:
    python3 pose_solver.py \\
        --binder ../data/receptor.pdb --target ../data/effector.pdb \\
        --pairs "A73-B48 A71-B50 A8-B39" --max-pair-distance 4.0
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
from scipy.optimize import minimize, differential_evolution
from scipy.spatial import ConvexHull
from scipy.spatial.transform import Rotation


# ── Objective weights ─────────────────────────────────────────────────────────
W_VALIDITY = 1e6   # pair residue buried deep inside opposite hull → invalid pose
VALID_TOL  = -1.0  # Å — pair residue must be this far OUTSIDE the opposite hull (negative = outside)
W_DIST     = 1e4   # penalty for pair CA-CA > max_pair_distance
W_LOWER    = 1e4   # penalty for pair CA-CA < min_pair_distance (too close to be realistic)
W_EXCL     = 1e4   # penalty for exclusion pair CA-CA < min_excl_distance (too close)
W_INTERP   = 10.0  # overall CA interpenetration depth — overridden by --interp-weight
W_PAIR_SC_CLASH = 1e4  # heavy-atom clash between sidechains of paired residues
PAIR_SC_CLASH_TOL = 2.0  # Å — paired sidechain heavy atoms within this distance count as a clash

# Backbone atoms excluded when extracting sidechain-only heavy atoms.
_BACKBONE_ATOMS = frozenset({'N', 'CA', 'C', 'O', 'OXT'})


# ── CLI ───────────────────────────────────────────────────────────────────────
def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--binder',  required=True, type=Path)
    ap.add_argument('--target',  required=True, type=Path)
    ap.add_argument('--binder-chain', default='A')
    ap.add_argument('--target-chain',  default='B')
    ap.add_argument('--pairs', required=True,
                    help='Space-separated ARES-BRES or ARES-BRES@DIST, e.g. "A73-B48 A71-B50"')
    ap.add_argument('--min-pair-distance', type=float, default=3.5,
                    help='Minimum allowed CA-CA distance for pair constraints (Å); below this is unrealistic')
    ap.add_argument('--max-pair-distance', type=float, default=4.0,
                    help='Maximum allowed CA-CA distance for pair constraints (Å)')
    ap.add_argument('--exclusions', default='',
                    help='Space-separated ARES-BRES@DIST pairs that must stay FURTHER than DIST Å, '
                         'e.g. "A8-B33@4.5" penalises if A8 and B33 come within 4.5 Å')
    ap.add_argument('--clash-cutoff', type=float, default=2.0,
                    help='Heavy-atom distance counted as a clash (Å)')
    ap.add_argument('--pair-sc-clash-cutoff', type=float, default=PAIR_SC_CLASH_TOL,
                    help='Min heavy-atom distance between sidechains of paired residues (Å); '
                         'enforced during optimisation. Set to 0 to disable.')
    ap.add_argument('--contact-cutoff', type=float, default=8.0,
                    help='CA-CA distance shown on heatmap (Å)')
    ap.add_argument('--contig-design-region', default='',
                    help='Binder residue ranges already planned for redesign, e.g. "33-49,69-78"')
    ap.add_argument('--n-restarts', type=int, default=300,
                    help='Number of random restarts (L-BFGS-B) or population×generations (DE)')
    ap.add_argument('--global-interp', action='store_true',
                    help='Penalise ALL CA atoms inside the opposite hull (not just pair residues)')
    ap.add_argument('--interp-weight', type=float, default=10.0,
                    help='Weight for interpenetration penalty (default 10.0)')
    ap.add_argument('--use-de', action='store_true',
                    help='Use differential evolution instead of random restarts')
    ap.add_argument('--output-prefix', default='solved_pose')
    return ap.parse_args()


# ── PDB I/O ──────────────────────────────────────────────────────────────────
def read_ca(pdb: Path, chain: str) -> dict[int, np.ndarray]:
    cas: dict[int, np.ndarray] = {}
    with open(pdb) as f:
        for line in f:
            if line[:4] != 'ATOM': continue
            if line[21] != chain: continue
            if line[12:16].strip() != 'CA': continue
            try:
                cas[int(line[22:26])] = np.array([float(line[30:38]),
                                                   float(line[38:46]),
                                                   float(line[46:54])])
            except ValueError:
                pass
    return cas


def read_heavy(pdb: Path, chain: str) -> dict[int, np.ndarray]:
    data: dict[int, list] = {}
    with open(pdb) as f:
        for line in f:
            if line[:4] != 'ATOM': continue
            if line[21] != chain: continue
            atom = line[12:16].strip()
            elem = line[76:78].strip()
            if elem == 'H' or (not elem and atom.startswith('H')): continue
            try:
                rn = int(line[22:26])
                data.setdefault(rn, []).append([float(line[30:38]),
                                                float(line[38:46]),
                                                float(line[46:54])])
            except ValueError:
                pass
    return {rn: np.array(v) for rn, v in data.items()}


def read_sidechain_heavy(pdb: Path, chain: str) -> dict[int, np.ndarray]:
    """Read sidechain heavy atoms only (excludes backbone N, CA, C, O, OXT)."""
    data: dict[int, list] = {}
    with open(pdb) as f:
        for line in f:
            if line[:4] != 'ATOM': continue
            if line[21] != chain: continue
            atom = line[12:16].strip()
            elem = line[76:78].strip()
            if elem == 'H' or (not elem and atom.startswith('H')): continue
            if atom in _BACKBONE_ATOMS: continue
            try:
                rn = int(line[22:26])
                data.setdefault(rn, []).append([float(line[30:38]),
                                                float(line[38:46]),
                                                float(line[46:54])])
            except ValueError:
                pass
    return {rn: np.array(v) for rn, v in data.items()}


def write_posed_pdb(out: Path, binder_pdb: Path, target_pdb: Path,
                    binder_chain: str, target_chain: str,
                    R: np.ndarray, t: np.ndarray) -> None:
    lines = []
    with open(binder_pdb) as f:
        for line in f:
            if line[:4] == 'ATOM' and line[21] == binder_chain:
                lines.append(line)
    lines.append('TER\n')
    with open(target_pdb) as f:
        for line in f:
            if line[:4] == 'ATOM' and line[21] == target_chain:
                try:
                    xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                    nxyz = R @ xyz + t
                    # Preserve the target chain ID end-to-end.  Callers
                    # must ensure binder_chain != target_chain to avoid
                    # producing a complex with two same-named chains.
                    lines.append(line[:21] + target_chain
                                 + line[22:30]
                                 + f'{nxyz[0]:8.3f}{nxyz[1]:8.3f}{nxyz[2]:8.3f}'
                                 + line[54:])
                except ValueError:
                    lines.append(line)
    lines.append('TER\nEND\n')
    # Renumber atom serials sequentially to avoid duplicate serial warnings
    serial = 1
    renumbered = []
    for line in lines:
        if line[:4] in ('ATOM', 'HETA'):
            line = f'{line[:6]}{serial:5d}{line[11:]}'
            serial += 1
        renumbered.append(line)
    with open(out, 'w') as f:
        f.writelines(renumbered)


# ── Pair parsing ──────────────────────────────────────────────────────────────
def parse_pairs(s: str, default_dist: float) -> list[tuple[int, int, float]]:
    out = []
    for tok in s.split():
        m = re.match(r'[A-Za-z]?(\d+)-[A-Za-z]?(\d+)(?:@([\d.]+))?$', tok.strip())
        if not m:
            sys.exit(f"Cannot parse pair {tok!r}. Expected e.g. A73-B48 or A73-B48@4.0")
        a, b, d = m.groups()
        out.append((int(a), int(b), float(d) if d else default_dist))
    return out


def parse_residue_ranges(spec: str) -> list[int]:
    res: list[int] = []
    for tok in spec.split(','):
        tok = tok.strip()
        if not tok: continue
        if '-' in tok:
            lo, hi = tok.split('-', 1)
            res.extend(range(int(lo), int(hi) + 1))
        else:
            res.append(int(tok))
    return sorted(set(res))


# ── Convex hull depth ─────────────────────────────────────────────────────────
def hull_depths(points: np.ndarray, hull_equations: np.ndarray) -> np.ndarray:
    """Unsigned depth: 0 if outside, positive if inside. Used for W_INTERP."""
    scores = points @ hull_equations[:, :3].T + hull_equations[:, 3]
    return np.maximum(0.0, -scores.max(axis=1))


def hull_signed_depths(points: np.ndarray, hull_equations: np.ndarray) -> np.ndarray:
    """Signed depth: negative = outside (how far), positive = inside (how deep).
    Used for validity check — negative values are good, positive are bad."""
    scores = points @ hull_equations[:, :3].T + hull_equations[:, 3]
    return -scores.max(axis=1)


# ── Objective function ────────────────────────────────────────────────────────
def _obj(params: np.ndarray,
         b_arr: np.ndarray,          # (N_b, 3) binder CA — fixed
         t_arr: np.ndarray,          # (N_t, 3) target CA — original coords
         pair_b_idx: np.ndarray,     # indices into b_arr
         pair_t_idx: np.ndarray,     # indices into t_arr
         min_pair_dist: float,
         max_pair_dist: float,
         b_hull_eq: np.ndarray,      # binder hull equations
         t_hull_eq: np.ndarray,      # target hull equations (original frame)
         excl_b_idx: np.ndarray,     # exclusion binder residue indices
         excl_t_idx: np.ndarray,     # exclusion target residue indices
         excl_min_dists: np.ndarray, # minimum allowed distances for exclusions
         global_interp: bool,        # if True, penalise all CAs; if False, pair residues only
         interp_weight: float,       # weight for interpenetration penalty
         b_pair_sc_arr: np.ndarray,  # (N_b_sc, 3) flat sidechain heavy atoms of binder pair residues
         t_pair_sc_arr: np.ndarray,  # (N_t_sc, 3) flat sidechain heavy atoms of target pair residues
         pair_sc_offsets: list,      # list of (b_start, b_end, t_start, t_end) per pair
         pair_sc_clash_cutoff: float, # min allowed heavy-atom distance for paired sidechains
         return_components: bool = False,  # dict of per-term losses (end-of-run only)
         ):

    R = Rotation.from_rotvec(params[:3]).as_matrix()
    t = params[3:]

    # Transformed target CAs in world frame
    target_t = (R @ t_arr.T).T + t

    # Binder CAs in target original frame (inverse transform).
    # Checking binder[i] inside transformed target hull is equivalent to
    # checking R^T(binder[i]-t) inside the original target hull.
    b_in_t_frame = (R.T @ (b_arr - t).T).T

    # Pair validity: penalise unless pair residue is at least |VALID_TOL| Å OUTSIDE opposite hull.
    # Signed depth: negative = outside (good), positive = inside (bad).
    # Penalty triggers when signed_depth > VALID_TOL (e.g. VALID_TOL=-2 means must be 2 Å outside).
    sdepth_b = hull_signed_depths(b_in_t_frame[pair_b_idx], t_hull_eq)
    sdepth_t = hull_signed_depths(target_t[pair_t_idx], b_hull_eq)
    validity_loss = W_VALIDITY * float(
        np.maximum(0.0, sdepth_b - VALID_TOL).sum() +
        np.maximum(0.0, sdepth_t - VALID_TOL).sum()
    )

    # Pair distances: penalise outside the window [min_pair_dist, max_pair_dist].
    # Too far (> max): unreachable contact.  Too close (< min): physically unrealistic.
    dists = np.linalg.norm(b_arr[pair_b_idx] - target_t[pair_t_idx], axis=1)
    upper_violations = np.maximum(0.0, dists - max_pair_dist)
    lower_violations = np.maximum(0.0, min_pair_dist - dists)
    dist_loss = (W_DIST  * float((upper_violations ** 2).sum()) +
                 W_LOWER * float((lower_violations ** 2).sum()))

    # Exclusion: penalise if excluded pair CA-CA < min distance (must stay apart)
    if len(excl_b_idx):
        excl_dists = np.linalg.norm(b_arr[excl_b_idx] - target_t[excl_t_idx], axis=1)
        excl_violations = np.maximum(0.0, excl_min_dists - excl_dists)
        excl_loss = W_EXCL * float((excl_violations ** 2).sum())
    else:
        excl_loss = 0.0

    # Interpenetration — global (all CAs) or pair residues only
    if global_interp:
        interp_loss = interp_weight * float(
            hull_depths(b_in_t_frame, t_hull_eq).sum() +
            hull_depths(target_t, b_hull_eq).sum()
        )
    else:
        interp_loss = interp_weight * float(
            hull_depths(b_in_t_frame[pair_b_idx], t_hull_eq).sum() +
            hull_depths(target_t[pair_t_idx], b_hull_eq).sum()
        )

    # Pair sidechain clash: minimum heavy-atom distance between the sidechain
    # heavy atoms of each binder pair-residue and its partner target pair-residue.
    # Catches the over-packing case where CA-CA satisfies the [min, max] window
    # but bulky sidechains overlap.  Pairs involving GLY (no sidechain) skipped.
    pair_sc_loss = 0.0
    if pair_sc_clash_cutoff > 0.0 and pair_sc_offsets and t_pair_sc_arr.shape[0] > 0:
        t_pair_sc_t = (R @ t_pair_sc_arr.T).T + t
        for bs, be, ts, te in pair_sc_offsets:
            if be <= bs or te <= ts:
                continue
            b_at = b_pair_sc_arr[bs:be]
            t_at = t_pair_sc_t[ts:te]
            diff = b_at[:, None, :] - t_at[None, :, :]
            min_d = float(np.sqrt((diff * diff).sum(-1)).min())
            if min_d < pair_sc_clash_cutoff:
                pair_sc_loss += (pair_sc_clash_cutoff - min_d) ** 2
        pair_sc_loss *= W_PAIR_SC_CLASH

    # Split dist_loss into upper/lower for the breakdown view.  Recomputing
    # only when explicitly asked keeps the hot optimisation path scalar.
    if return_components:
        dist_upper = W_DIST  * float((upper_violations ** 2).sum())
        dist_lower = W_LOWER * float((lower_violations ** 2).sum())
        return {
            'validity':      validity_loss,
            'dist_upper':    dist_upper,
            'dist_lower':    dist_lower,
            'excl':          excl_loss,
            'interp':        interp_loss,
            'pair_sc_clash': pair_sc_loss,
            'total':         (validity_loss + dist_loss + excl_loss
                              + interp_loss + pair_sc_loss),
        }

    return validity_loss + dist_loss + excl_loss + interp_loss + pair_sc_loss


# ── Solver ────────────────────────────────────────────────────────────────────
def solve(binder_ca: dict, target_ca: dict,
          pairs: list[tuple[int, int, float]],
          exclusions: list[tuple[int, int, float]],
          n_restarts: int, min_pair_dist: float, max_pair_dist: float,
          use_de: bool = False, global_interp: bool = False,
          interp_weight: float = 10.0,
          binder_sc: dict | None = None,
          target_sc: dict | None = None,
          pair_sc_clash_cutoff: float = 0.0,
          ) -> tuple[np.ndarray, np.ndarray, float, list[dict], list[tuple], tuple]:

    b_rn = sorted(binder_ca); t_rn = sorted(target_ca)
    b_arr = np.array([binder_ca[r] for r in b_rn])
    t_arr = np.array([target_ca[r] for r in t_rn])
    b_idx = {r: i for i, r in enumerate(b_rn)}
    t_idx = {r: i for i, r in enumerate(t_rn)}

    valid = [(rb, rt, d) for rb, rt, d in pairs if rb in b_idx and rt in t_idx]
    missing = [(rb, rt) for rb, rt, _ in pairs if rb not in b_idx or rt not in t_idx]
    if missing:
        print(f"  WARNING: residues not found, skipping: {missing}")
    if not valid:
        sys.exit("ERROR: no valid pairs found in PDB files.")

    pb = np.array([b_idx[rb] for rb, _, _ in valid])
    pt = np.array([t_idx[rt] for _, rt, _ in valid])

    b_hull = ConvexHull(b_arr)
    t_hull = ConvexHull(t_arr)

    # Build exclusion arrays
    excl_valid = [(rb, rt, d) for rb, rt, d in exclusions if rb in b_idx and rt in t_idx]
    if excl_valid:
        excl_b_idx   = np.array([b_idx[rb] for rb, _, _ in excl_valid])
        excl_t_idx   = np.array([t_idx[rt] for _, rt, _ in excl_valid])
        excl_min_d   = np.array([d          for _, _, d  in excl_valid])
    else:
        excl_b_idx = excl_t_idx = excl_min_d = np.array([], dtype=int)

    # Build flat sidechain heavy-atom arrays for pair residues, with per-pair offsets.
    b_pair_sc_list: list = []
    t_pair_sc_list: list = []
    pair_sc_offsets: list[tuple[int, int, int, int]] = []
    pair_sc_skipped: list[tuple[int, int, str]] = []
    if binder_sc is not None and target_sc is not None and pair_sc_clash_cutoff > 0.0:
        for rb, rt, _ in valid:
            b_at = binder_sc.get(rb)
            t_at = target_sc.get(rt)
            bs = len(b_pair_sc_list)
            ts = len(t_pair_sc_list)
            if b_at is None or t_at is None or len(b_at) == 0 or len(t_at) == 0:
                pair_sc_skipped.append((rb, rt, 'no sidechain (likely GLY)'))
                pair_sc_offsets.append((bs, bs, ts, ts))
                continue
            b_pair_sc_list.extend(b_at.tolist())
            t_pair_sc_list.extend(t_at.tolist())
            pair_sc_offsets.append((bs, bs + len(b_at), ts, ts + len(t_at)))
    b_pair_sc_arr = (np.array(b_pair_sc_list, dtype=float)
                     if b_pair_sc_list else np.empty((0, 3), dtype=float))
    t_pair_sc_arr = (np.array(t_pair_sc_list, dtype=float)
                     if t_pair_sc_list else np.empty((0, 3), dtype=float))

    obj_args = (b_arr, t_arr, pb, pt, min_pair_dist, max_pair_dist,
                b_hull.equations, t_hull.equations,
                excl_b_idx, excl_t_idx, excl_min_d, global_interp, interp_weight,
                b_pair_sc_arr, t_pair_sc_arr, pair_sc_offsets, pair_sc_clash_cutoff)

    # Smart initialisation: for every restart, start with the target's pair
    # residue centroid placed near the binder's pair residue centroid.
    # This seeds the optimiser in the right neighbourhood regardless of rotation.
    pair_b_centroid = b_arr[pb].mean(0)
    pair_t_centroid = t_arr[pt].mean(0)

    # Outward normal at the binder pair region (approximate surface normal)
    outward = pair_b_centroid - b_arr.mean(0)
    outward /= np.linalg.norm(outward) + 1e-12

    # Kabsch alignment: finds the rotation/translation that minimises the sum
    # of squared distances between all pair residues simultaneously.
    # Used as restart 0 — places all pairs at distance ~0, giving the
    # optimiser the best possible starting point.
    def kabsch(A: np.ndarray, B: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        cA, cB = A.mean(0), B.mean(0)
        H = (B - cB).T @ (A - cA)
        U, _, Vt = np.linalg.svd(H)
        d = np.sign(np.linalg.det(Vt.T @ U.T))
        R_k = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
        return R_k, cA - R_k @ cB

    R_k, t_k = kabsch(b_arr[pb], t_arr[pt])
    rotvec_k = Rotation.from_matrix(R_k).as_rotvec()

    # Interface normal: perpendicular to the plane spanned by the three binder pair positions.
    # In the Kabsch pose the pairs coincide, so moving the target by n*d gives all pairs at distance d.
    if len(pb) >= 3:
        p1, p2, p3 = b_arr[pb[0]], b_arr[pb[1]], b_arr[pb[2]]
        n = np.cross(p2 - p1, p3 - p1)
        n = n / (np.linalg.norm(n) + 1e-12)
    else:
        n = outward  # fallback for < 3 pairs

    best_loss = np.inf
    best_x: np.ndarray | None = None

    np.random.seed(42)
    interp_scope = 'global' if global_interp else 'pair-only'
    print(f"  Weights: W_validity={W_VALIDITY:.0e} tol={VALID_TOL}Å  "
          f"W_dist={W_DIST:.0e}  W_lower={W_LOWER:.0e}  W_interp={interp_weight} ({interp_scope})")
    if pair_sc_clash_cutoff > 0.0:
        print(f"  Pair sidechain-clash: W={W_PAIR_SC_CLASH:.0e}  cutoff={pair_sc_clash_cutoff:.2f}Å  "
              f"(active pairs: {sum(1 for bs, be, ts, te in pair_sc_offsets if be > bs and te > ts)}"
              f"/{len(valid)})")
        for rb, rt, reason in pair_sc_skipped:
            print(f"    skipped {rb}-{rt}: {reason}")

    # Per-restart loss history.  Recorded so POSE_SOLVER_PLOTS can render
    # the restart-loss convergence curve.  In the DE path we get one
    # final value (no restart-level history); a single-row CSV is fine.
    restart_history: list[tuple[int, float, float]] = []

    if use_de:
        # Differential evolution: global optimiser, no gradients required.
        # Search bounds: rotation in [-π, π]³, translation ±50 Å around binder pair centroid.
        lo = pair_b_centroid - 50.0
        hi = pair_b_centroid + 50.0
        bounds = [(-np.pi, np.pi)] * 3 + list(zip(lo, hi))
        popsize = 15  # population = popsize × 6
        maxiter = max(1, n_restarts // popsize)
        print(f"  Differential evolution: popsize={popsize}  maxiter={maxiter}  "
              f"(effective evaluations ≈ {popsize*6*maxiter})")

        def cb(xk, convergence):
            pass  # silent

        res = differential_evolution(
            _obj, bounds, args=obj_args,
            maxiter=maxiter, popsize=popsize,
            seed=42, tol=1e-9, mutation=(0.5, 1.0), recombination=0.9,
            workers=1, polish=True, init='sobol',
            callback=cb,
        )
        best_x = res.x
        best_loss = res.fun
        restart_history.append((0, float(res.fun), float(res.fun)))
        print(f"  DE converged: loss = {best_loss:.4f}  (success={res.success})")
    else:
        for i in range(n_restarts):
            if i == 0:
                x0 = np.r_[rotvec_k, t_k + n * 30.0]
            elif i == 1:
                x0 = np.r_[rotvec_k, t_k - n * 30.0]
            elif i < 20:
                rand_dir = np.random.randn(3); rand_dir /= np.linalg.norm(rand_dir) + 1e-12
                x0 = np.r_[rotvec_k, t_k + rand_dir * 30.0]
            else:
                angle = np.random.uniform(0, 2 * np.pi)
                axis = np.random.randn(3); axis /= np.linalg.norm(axis) + 1e-12
                rotvec = axis * angle
                R_init = Rotation.from_rotvec(rotvec).as_matrix()
                rand_dir = np.random.randn(3); rand_dir /= np.linalg.norm(rand_dir) + 1e-12
                t_init = (pair_b_centroid + rand_dir * 30.0
                          - R_init @ pair_t_centroid)
                x0 = np.r_[rotvec, t_init]

            res = minimize(_obj, x0, args=obj_args, method='L-BFGS-B',
                           options={'maxiter': 2000, 'ftol': 1e-12, 'gtol': 1e-8})
            if res.fun < best_loss:
                best_loss = res.fun
                best_x = res.x.copy()
            restart_history.append((i, float(res.fun), float(best_loss)))

            if (i + 1) % 50 == 0 or i == n_restarts - 1:
                print(f"    {i+1:3d}/{n_restarts}  best loss = {best_loss:.4f}")

    R = Rotation.from_rotvec(best_x[:3]).as_matrix()
    t = best_x[3:]

    target_t = (R @ t_arr.T).T + t
    pair_results = []
    for (rb, rt, _), bi, ti in zip(valid, pb, pt):
        achieved = float(np.linalg.norm(b_arr[bi] - target_t[ti]))
        pair_results.append({
            'binder_res': rb, 'target_res': rt,
            'achieved_d': round(achieved, 3),
            'min_d': min_pair_dist, 'max_d': max_pair_dist,
            'satisfied': min_pair_dist - 0.01 <= achieved <= max_pair_dist + 0.01,
        })

    # Per-term loss breakdown at the final pose (one extra _obj() call;
    # negligible overhead compared to optimisation).
    loss_breakdown = _obj(best_x, *obj_args, return_components=True)

    return R, t, best_loss, pair_results, restart_history, loss_breakdown


# ── Clash / contact analysis ──────────────────────────────────────────────────
def analyse_contacts(binder_heavy: dict, target_heavy: dict,
                     R: np.ndarray, t: np.ndarray,
                     contact_cutoff: float) -> dict[tuple[int, int], float]:
    contacts: dict[tuple[int, int], float] = {}
    for br, bc in binder_heavy.items():
        for tr, tc_orig in target_heavy.items():
            tc = (R @ tc_orig.T).T + t
            diffs = bc[:, None, :] - tc[None, :, :]
            d = float(np.sqrt((diffs ** 2).sum(axis=2)).min())
            if d <= contact_cutoff:
                contacts[(br, tr)] = d
    return contacts


# ── Contig generation ─────────────────────────────────────────────────────────
def expand_with_gap_join(design: set[int], all_rn: list[int], max_gap: int) -> set[int]:
    expanded = set(design)
    if max_gap == 0:
        return expanded
    n = len(all_rn)
    changed = True
    while changed:
        changed = False
        i = 0
        while i < n:
            if all_rn[i] not in expanded:
                j = i
                while j < n and all_rn[j] not in expanded:
                    j += 1
                before = i > 0 and all_rn[i - 1] in expanded
                after  = j < n and all_rn[j] in expanded
                if (j - i) <= max_gap and before and after:
                    for k in range(i, j):
                        expanded.add(all_rn[k])
                    changed = True
                i = j
            else:
                i += 1
    return expanded


def make_contig(all_binder_rn: list[int], design_region: set[int],
                binder_chain: str, target_chain: str, join_gap: int) -> str:
    dr = expand_with_gap_join(design_region, all_binder_rn, join_gap)
    parts: list[str] = []
    in_denovo: bool | None = None
    current: list[int] = []
    for rn in sorted(all_binder_rn):
        is_denovo = rn in dr
        if is_denovo != in_denovo:
            if current:
                parts.append(f'{len(current)}-{len(current)}' if in_denovo
                              else f'{binder_chain}{current[0]}-{current[-1]}')
            current = [rn]; in_denovo = is_denovo
        else:
            current.append(rn)
    if current:
        parts.append(f'{len(current)}-{len(current)}' if in_denovo
                     else f'{binder_chain}{current[0]}-{current[-1]}')
    return '/'.join(parts) + f' {target_chain}'


def describe_design_region(all_binder_rn: list[int], design_region: set[int],
                            join_gap: int) -> tuple[str, int]:
    dr = expand_with_gap_join(design_region, all_binder_rn, join_gap)
    denovo = sorted(r for r in all_binder_rn if r in dr)
    if not denovo:
        return '(none)', 0
    runs, s, p = [], denovo[0], denovo[0]
    for r in denovo[1:]:
        if r != p + 1:
            runs.append(f'{s}-{p}' if s != p else str(s))
            s = r
        p = r
    runs.append(f'{s}-{p}' if s != p else str(p))
    return ', '.join(runs), len(denovo)


# ── Heatmap ───────────────────────────────────────────────────────────────────
_CMAP = LinearSegmentedColormap.from_list('contact', [
    (0.00, '#7F0000'), (0.25, '#DC2626'), (0.50, '#FDE68A'),
    (0.75, '#93C5FD'), (1.00, '#F0F9FF'),
])


def make_heatmap(contacts: dict, pairs: list[tuple[int, int, float]],
                 contig_dr: list[int], clash_binder: set[int], clash_target: set[int],
                 binder_chain: str, target_chain: str,
                 clash_cutoff: float, contact_cutoff: float,
                 pair_results: list[dict], out_path: Path) -> None:
    from matplotlib.colors import ListedColormap as LCM

    b_rns = sorted({b for b, _ in contacts})
    t_rns = sorted({t for _, t in contacts})
    if not b_rns:
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, 'No contacts', ha='center', va='center')
        plt.savefig(out_path, dpi=150); plt.close()
        return

    bi = {r: i for i, r in enumerate(b_rns)}
    ti = {r: i for i, r in enumerate(t_rns)}
    grid = np.full((len(t_rns), len(b_rns)), np.nan)
    for (br, tr), d in contacts.items():
        if br in bi and tr in ti:
            grid[ti[tr], bi[br]] = d

    pair_b = {rb for rb, _, _ in pairs}
    pair_t = {rt for _, rt, _ in pairs}
    design_set = set(contig_dr) | clash_binder

    fw = max(11, len(b_rns) * 0.075 + 4)
    fh = max(7, len(t_rns) * 0.13 + 3.5)
    fig, (ax_bar, ax_t_bar, ax_heat) = plt.subplots(
        3, 1, figsize=(fw, fh),
        gridspec_kw={'height_ratios': [0.22, 0.22, 1], 'hspace': 0.02})

    # Binder annotation bar
    bbar = np.zeros((1, len(b_rns)))
    for i, r in enumerate(b_rns):
        bbar[0, i] = 2 if r in clash_binder else (1 if r in set(contig_dr) else 0)
    ax_bar.imshow(bbar, aspect='auto',
                  cmap=LCM(['#DDDDDD', '#FF8C00', '#DC2626']),
                  vmin=0, vmax=2, interpolation='nearest',
                  extent=[-0.5, len(b_rns) - 0.5, 0, 1])
    ax_bar.set_xlim(-0.5, len(b_rns) - 0.5)
    ax_bar.set_yticks([]); ax_bar.tick_params(bottom=False, labelbottom=False)
    for sp in ['top', 'left', 'right']: ax_bar.spines[sp].set_visible(False)
    ax_bar.set_ylabel('Binder\nregions', fontsize=8, rotation=0,
                       ha='right', va='center', labelpad=30)
    for r in pair_b:
        if r in bi:
            ax_bar.scatter(bi[r], 1.25, marker='v', color='black', s=35, zorder=6, clip_on=False)
    ax_bar.legend(handles=[
        mpatches.Patch(color='#DDDDDD', label='Fixed'),
        mpatches.Patch(color='#FF8C00', label='Contig design region'),
        mpatches.Patch(color='#DC2626', label='Clash-derived'),
    ], fontsize=7, framealpha=0.9, loc='upper left',
        bbox_to_anchor=(1.01, 1.0), borderaxespad=0)

    # Target annotation bar
    tbar = np.zeros((1, len(t_rns)))
    for i, r in enumerate(t_rns):
        tbar[0, i] = 2 if r in clash_target else (1 if r in pair_t else 0)
    ax_t_bar.imshow(tbar, aspect='auto',
                    cmap=LCM(['#DDDDDD', '#8B5CF6', '#DC2626']),
                    vmin=0, vmax=2, interpolation='nearest',
                    extent=[-0.5, len(t_rns) - 0.5, 0, 1])
    ax_t_bar.set_xlim(-0.5, len(t_rns) - 0.5)
    ax_t_bar.set_yticks([]); ax_t_bar.tick_params(bottom=False, labelbottom=False)
    for sp in ['top', 'left', 'right']: ax_t_bar.spines[sp].set_visible(False)
    ax_t_bar.set_ylabel('Target\nregions', fontsize=8, rotation=0,
                         ha='right', va='center', labelpad=30)
    for r in pair_t:
        if r in ti:
            ax_t_bar.scatter(ti[r], 1.25, marker='v', color='black', s=35, zorder=6, clip_on=False)
    ax_t_bar.legend(handles=[
        mpatches.Patch(color='#DDDDDD', label='Free'),
        mpatches.Patch(color='#8B5CF6', label='Pair anchors'),
        mpatches.Patch(color='#DC2626', label='Clashing'),
    ], fontsize=7, framealpha=0.9, loc='upper left',
        bbox_to_anchor=(1.01, 1.0), borderaxespad=0)

    # Heatmap
    _cmap = _CMAP.copy(); _cmap.set_bad('#F8F8F8')
    im = ax_heat.imshow(grid, aspect='auto', cmap=_cmap,
                        vmin=0, vmax=contact_cutoff, origin='lower',
                        interpolation='nearest',
                        extent=[-0.5, len(b_rns) - 0.5, -0.5, len(t_rns) - 0.5])
    for rb, rt, _ in pairs:
        if rb in bi and rt in ti:
            ax_heat.scatter(bi[rb], ti[rt], marker='*', color='gold',
                            s=300, zorder=10, edgecolors='black', linewidths=0.8)
    for (br, tr), d in contacts.items():
        if d < clash_cutoff and br in bi and tr in ti:
            ax_heat.add_patch(Rectangle((bi[br] - 0.5, ti[tr] - 0.5), 1, 1,
                lw=1.8, edgecolor='#DC2626', facecolor='none', zorder=8))

    def _ticks(items, mx=20):
        step = max(1, len(items) // mx)
        pos = list(range(0, len(items), step))
        if pos[-1] != len(items) - 1: pos.append(len(items) - 1)
        return pos

    tp = _ticks(b_rns)
    ax_heat.set_xticks(tp)
    ax_heat.set_xticklabels([b_rns[i] for i in tp], rotation=45, ha='right', fontsize=7)
    ax_heat.set_xlabel(f'Binder residue (chain {binder_chain})', fontsize=10)
    tp_t = _ticks(t_rns)
    ax_heat.set_yticks(tp_t)
    ax_heat.set_yticklabels([t_rns[i] for i in tp_t], fontsize=7)
    ax_heat.set_ylabel(f'Target residue (chain {target_chain})', fontsize=10)
    ax_heat.set_xlim(-0.5, len(b_rns) - 0.5)
    ax_heat.set_ylim(-0.5, len(t_rns) - 0.5)

    cbar = fig.colorbar(im, ax=[ax_bar, ax_t_bar, ax_heat],
                        fraction=0.018, pad=0.01, aspect=30)
    cbar.set_label('Min heavy-atom distance (Å)', fontsize=9)
    cbar.ax.axhline(clash_cutoff, color='#DC2626', lw=1.5, ls='--')

    ax_heat.legend(handles=[
        plt.Line2D([0], [0], marker='*', color='w', markerfacecolor='gold',
                   markeredgecolor='black', markersize=13, label='Pair pin'),
        mpatches.Patch(edgecolor='#DC2626', facecolor='none', lw=2,
                       label=f'Clash (<{clash_cutoff} Å)'),
    ], fontsize=8, framealpha=0.9, loc='lower right')

    if pair_results:
        txt = ['Pair results:']
        for pr in pair_results:
            flag = '✓' if pr['satisfied'] else '!'
            txt.append(f"  B{pr['binder_res']}-T{pr['target_res']}: "
                       f"{pr['achieved_d']:.2f} Å (max {pr['max_d']} Å) {flag}")
        ax_heat.text(0.01, 0.98, '\n'.join(txt), transform=ax_heat.transAxes,
                     fontsize=7.5, va='top', ha='left',
                     bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                               alpha=0.88, edgecolor='#CCCCCC'))

    fig.suptitle('Constraint-solved pose — contact/clash heatmap', fontsize=10, y=1.01)
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved heatmap: {out_path}")


# ── Main ─────────────────────────────────────────────────────────────────────
def main() -> None:
    args = parse_args()
    bc = args.binder_chain.upper()
    tc = args.target_chain.upper()
    if bc == tc:
        sys.exit(f"ERROR: binder-chain and target-chain are both '{bc}'.  "
                 f"Output would contain two same-named chains.  Pass distinct "
                 f"chain IDs (rename inputs first if necessary).")
    pairs      = parse_pairs(args.pairs, args.max_pair_distance)
    exclusions = parse_pairs(args.exclusions, 4.5) if args.exclusions else []
    contig_dr  = parse_residue_ranges(args.contig_design_region)

    print('=' * 64)
    print('CONSTRAINT-DRIVEN POSE SOLVER')
    print('=' * 64)
    print(f"Binder : {args.binder}  (chain {bc})")
    print(f"Target : {args.target}  (chain {tc})")
    for rb, rt, d in pairs:
        print(f"  Pair      {bc}{rb} — {tc}{rt}  max {d:.1f} Å")
    for rb, rt, d in exclusions:
        print(f"  Exclusion {bc}{rb} — {tc}{rt}  min {d:.1f} Å")
    print()

    binder_ca    = read_ca(args.binder, bc)
    target_ca    = read_ca(args.target, tc)
    binder_heavy = read_heavy(args.binder, bc)
    target_heavy = read_heavy(args.target, tc)
    binder_sc    = read_sidechain_heavy(args.binder, bc)
    target_sc    = read_sidechain_heavy(args.target, tc)

    print(f"Binder: {len(binder_ca)} residues   Target: {len(target_ca)} residues")
    print()

    print("Solving pose...")
    R, t, final_loss, pair_results, restart_history, loss_breakdown = solve(
        binder_ca, target_ca, pairs, exclusions, args.n_restarts,
        args.min_pair_distance, args.max_pair_distance,
        use_de=args.use_de, global_interp=args.global_interp,
        interp_weight=args.interp_weight,
        binder_sc=binder_sc, target_sc=target_sc,
        pair_sc_clash_cutoff=args.pair_sc_clash_cutoff)
    print()

    all_satisfied = all(pr['satisfied'] for pr in pair_results)
    print(f"Final loss: {final_loss:.4f}  —  all pairs satisfied: {all_satisfied}")
    for pr in pair_results:
        flag = '✓' if pr['satisfied'] else '!'
        print(f"  {bc}{pr['binder_res']:3d} — {tc}{pr['target_res']:3d}: "
              f"achieved {pr['achieved_d']:.3f} Å  (max {pr['max_d']:.1f} Å)  {flag}")
    print()

    # Per-term loss breakdown (per project memory feedback_loss_breakdown:
    # report unprompted after every run so the user can judge which term
    # dominates and whether weights need rebalancing).
    print("Loss breakdown:")
    for term, val in loss_breakdown.items():
        if term == 'total':
            continue
        print(f"  {term:<14s} {val:>12.4f}")
    print(f"  {'total':<14s} {loss_breakdown['total']:>12.4f}")
    print()

    print("Analysing contacts and clashes...")
    contacts = analyse_contacts(binder_heavy, target_heavy, R, t, args.contact_cutoff)
    clash_pairs = {k for k, d in contacts.items() if d < args.clash_cutoff}
    clash_b = sorted({b for b, _ in clash_pairs})
    clash_t = sorted({t_ for _, t_ in clash_pairs})
    design_set = set(clash_b) | set(contig_dr)
    all_binder_rn = sorted(binder_ca)

    print(f"  Contacts ≤ {args.contact_cutoff} Å : {len(contacts)} residue pairs")
    print(f"  Clashes  <  {args.clash_cutoff} Å  : "
          f"{len(clash_pairs)} pairs   {len(clash_b)} binder res   {len(clash_t)} target res")
    print(f"  Clashing binder residues: {clash_b}")
    print()

    print('-' * 64)
    print('RFDiffusion CONTIG STRINGS')
    print('-' * 64)
    contigs_out: dict[str, dict] = {}
    for label, gap in [
        ('A) No joining', 0),
        ('B) Join gap ≤ 1', 1),
        ('C) Join gap ≤ 2', 2),
    ]:
        contig = make_contig(all_binder_rn, design_set, bc, tc, gap)
        dr_desc, dr_total = describe_design_region(all_binder_rn, design_set, gap)
        contigs_out[label] = {'contig': contig, 'design_region': dr_desc,
                              'n_residues': dr_total, 'join_gap': gap}
        print(f"\n  {label}")
        print(f"    Design region : {dr_desc}  [{dr_total} residues]")
        print(f"    Contig string : {contig}")

    posed_pdb = Path(f'{args.output_prefix}_posed.pdb')
    write_posed_pdb(posed_pdb, args.binder, args.target, bc, tc, R, t)
    print(f"\nSaved posed complex : {posed_pdb}")

    results = {
        'pair_results': pair_results,
        'all_pairs_satisfied': all_satisfied,
        'final_loss': round(final_loss, 4),
        'loss_breakdown': {k: round(v, 4) for k, v in loss_breakdown.items()},
        'loss_weights': {
            'W_VALIDITY':       W_VALIDITY,
            'W_DIST':           W_DIST,
            'W_LOWER':          W_LOWER,
            'W_EXCL':           W_EXCL,
            'W_INTERP':         args.interp_weight,
            'W_PAIR_SC_CLASH':  W_PAIR_SC_CLASH,
        },
        'n_contact_pairs': len(contacts),
        'n_clash_pairs': len(clash_pairs),
        'clash_binder_residues': clash_b,
        'clash_target_residues': clash_t,
        'contig_design_region_input': contig_dr,
        'clash_cutoff_A': args.clash_cutoff,
        'contact_cutoff_A': args.contact_cutoff,
        'contigs': contigs_out,
    }
    json_out = Path(f'{args.output_prefix}_results.json')
    json_out.write_text(json.dumps(results, indent=2))
    print(f"Saved results JSON  : {json_out}")

    # Per-restart loss history CSV — consumed by POSE_SOLVER_PLOTS for
    # the restart-loss convergence curve.  Always emitted (DE path
    # writes one row, L-BFGS-B path writes n_restarts rows).
    csv_out = Path(f'{args.output_prefix}_restart_losses.csv')
    with open(csv_out, 'w') as fh:
        fh.write('restart_idx,restart_loss,best_loss_so_far\n')
        for i, restart_loss, best_so_far in restart_history:
            fh.write(f'{i},{restart_loss:.6f},{best_so_far:.6f}\n')
    print(f"Saved restart loss  : {csv_out}")

    make_heatmap(contacts, pairs, contig_dr,
                 set(clash_b), set(clash_t), bc, tc,
                 args.clash_cutoff, args.contact_cutoff,
                 pair_results, Path(f'{args.output_prefix}_heatmap.png'))

    print()
    print('=' * 64)
    print(f"Open in ChimeraX:  open {posed_pdb}")
    print('=' * 64)


if __name__ == '__main__':
    main()
