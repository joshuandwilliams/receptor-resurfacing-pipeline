"""
tests/characterization/test_pose_solver.py
------------------------------------------
Logic-tier tests for bin/pose_solver.py, the rigid-body constraint solver
that places a target monomer against a fixed binder.

Every geometry here is stated before the code runs, so the expected value is
arithmetic rather than a blessed previous output.

Two fixture shapes are used, for two different reasons:

  * `fixtures_runtree.complex_pdb` writes two chains on parallel lines a
    stated distance apart.  That makes a contact or clash count at a given
    cutoff a matter of counting lattice offsets, so it drives
    `analyse_contacts` and the PDB readers.
  * `fixtures_runtree.lattice_pdb` writes a cuboid point cloud.  Everything
    in the solver that takes a convex hull needs one, because Qhull rejects
    the collinear chains `complex_pdb` produces.

Loss-term note: the breakdown `_obj` returns has keys `validity`,
`dist_upper`, `dist_lower`, `excl`, `interp`, `pair_sc_clash` and `total`.
There is no centre-of-mass term.  COM-COM separation is a reported
diagnostic, not a penalty, and `test_breakdown_has_no_com_term` pins that so
a future COM penalty cannot be added silently.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.spatial import ConvexHull

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT / "bin") not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT / "bin"))

import pose_solver as ps  # noqa: E402

from tests.characterization.fixtures_runtree import (  # noqa: E402
    complex_pdb,
    lattice_pdb,
    lattice_resnum,
)

# ── Shared geometry ──────────────────────────────────────────────────
#
# A 4 Å cube: the eight corners of [0, 4]^3 plus its centre (2, 2, 2).
# The hull is the cube itself, so the signed depth of any point is
# readable off the face planes by hand.  The centre sits at index 8 and is
# 2 Å inside every face.
_CUBE_CORNERS = np.array([[x, y, z]
                          for x in (0.0, 4.0)
                          for y in (0.0, 4.0)
                          for z in (0.0, 4.0)])
_CUBE = np.vstack([_CUBE_CORNERS, [[2.0, 2.0, 2.0]]])
_CUBE_CENTRE_IDX = 8
_CUBE_EQ = ConvexHull(_CUBE).equations

# Lattice used for the solver and CLI tests: 3x3x3 points, 4 Å apart.
_SHAPE = (3, 3, 3)
_STEP = 4.0


def _rn(i: int, j: int, k: int) -> int:
    return lattice_resnum(i, j, k, shape=_SHAPE, start=1)


def _two_lattices(tmp_path: Path) -> tuple[Path, Path]:
    """A binder and a target, each a 3x3x3 lattice at the origin.

    Both occupy [0, 8]^3 in their own frame.  Residue numbering follows
    `lattice_resnum`, so residue `_rn(0, 0, 2)` is the corner (0, 0, 8) and
    residue `_rn(0, 0, 0)` is the corner (0, 0, 0).
    """
    binder = lattice_pdb(tmp_path / "binder.pdb", chain="A",
                         shape=_SHAPE, step=_STEP)
    target = lattice_pdb(tmp_path / "target.pdb", chain="B",
                         shape=_SHAPE, step=_STEP)
    return binder, target


def _obj_breakdown(*, pair_b=(_CUBE_CENTRE_IDX,), pair_t=(_CUBE_CENTRE_IDX,),
                   dz=10.0, min_d=3.5, max_d=4.0, b_arr=None, t_arr=None,
                   b_eq=None, t_eq=None, excl=(), global_interp=False,
                   interp_weight=10.0, b_sc=None, t_sc=None,
                   sc_offsets=(), sc_cutoff=0.0):
    """Evaluate `_obj` with the target translated by (0, 0, dz), no rotation."""
    b_arr = _CUBE if b_arr is None else b_arr
    t_arr = _CUBE if t_arr is None else t_arr
    b_eq = _CUBE_EQ if b_eq is None else b_eq
    t_eq = _CUBE_EQ if t_eq is None else t_eq
    eb = np.array([e[0] for e in excl], dtype=int)
    et = np.array([e[1] for e in excl], dtype=int)
    ed = np.array([e[2] for e in excl], dtype=float)
    return ps._obj(
        np.r_[0.0, 0.0, 0.0, 0.0, 0.0, dz],
        b_arr, t_arr, np.array(pair_b), np.array(pair_t),
        min_d, max_d, b_eq, t_eq, eb, et, ed,
        global_interp, interp_weight,
        np.empty((0, 3)) if b_sc is None else np.asarray(b_sc, dtype=float),
        np.empty((0, 3)) if t_sc is None else np.asarray(t_sc, dtype=float),
        list(sc_offsets), sc_cutoff,
        return_components=True,
    )


# ── Pair-token parsing ───────────────────────────────────────────────
#
# The repo's working notes record this token format as a recurring source of
# error, so it is pinned token by token rather than by example.


@pytest.mark.logic
def test_parse_pairs_single_token():
    """`A73-B31` is one pair that inherits the default distance."""
    assert ps.parse_pairs("A73-B31", 4.0) == [(73, 31, 4.0)]


@pytest.mark.logic
def test_parse_pairs_per_token_distance_overrides_default():
    """`@DIST` sets that pair's distance and leaves the default unused."""
    assert ps.parse_pairs("A73-B48@4.5", 4.0) == [(73, 48, 4.5)]
    assert ps.parse_pairs("A73-B48@12", 4.0) == [(73, 48, 12.0)]


@pytest.mark.logic
def test_parse_pairs_space_separated_preserves_order():
    """Several tokens are split on whitespace, in the order written.

    Mixed per-token and default distances coexist in one string.
    """
    out = ps.parse_pairs("A73-B31 A71-B33@4.5 A8-B22", 4.0)
    assert out == [(73, 31, 4.0), (71, 33, 4.5), (8, 22, 4.0)]


@pytest.mark.logic
def test_parse_pairs_tolerates_surrounding_and_repeated_whitespace():
    assert ps.parse_pairs("   A73-B31    A71-B33  ", 4.0) == [
        (73, 31, 4.0), (71, 33, 4.0)]


@pytest.mark.logic
def test_parse_pairs_chain_letters_are_optional_and_discarded():
    """Only the residue numbers survive parsing.

    The chain letter is accepted for readability and thrown away, so
    `A73-B31`, `73-31` and even `Z73-Q31` parse identically.  The chains the
    numbers are looked up in come from --binder-chain / --target-chain, not
    from the token.
    """
    expected = [(73, 31, 4.0)]
    assert ps.parse_pairs("A73-B31", 4.0) == expected
    assert ps.parse_pairs("73-31", 4.0) == expected
    assert ps.parse_pairs("Z73-Q31", 4.0) == expected


@pytest.mark.logic
@pytest.mark.parametrize("bad", [
    # Comma separation.  The notes record this as the recurring mistake:
    # the pair list is space-separated, and a comma makes the whole string
    # one unparseable token.
    "A73-B31,A71-B33",
    "A73-B31, A71-B33",
    "A73B31",          # no hyphen
    "A73-",            # no second residue
    "-B31",            # no first residue
    "A73-B31@",        # '@' with no distance
    "A73-B31@x",       # non-numeric distance
    "AA73-B31",        # two-letter chain id
    "A73-B31-C5",      # three members
    "A73--B31",        # doubled separator
    "A73:B31",         # wrong separator
])
def test_parse_pairs_rejects_malformed_tokens(bad):
    """Anything that is not `[CHAIN]NUM-[CHAIN]NUM[@DIST]` exits non-zero.

    `parse_pairs` calls `sys.exit` with a message rather than raising, so a
    typo in a pair string stops the run instead of silently dropping a
    restraint.
    """
    with pytest.raises(SystemExit) as exc:
        ps.parse_pairs(bad, 4.0)
    assert "Cannot parse pair" in str(exc.value)


@pytest.mark.logic
def test_parse_pairs_empty_string_is_no_pairs():
    """An empty string yields no pairs rather than an error."""
    assert ps.parse_pairs("", 4.0) == []
    assert ps.parse_pairs("   ", 4.0) == []


@pytest.mark.logic
def test_parse_residue_ranges_expands_and_dedupes():
    """Design-region ranges are comma-separated and inclusive at both ends.

    This is the one place in the solver where commas are the separator, in
    contrast to the space-separated pair list.
    """
    assert ps.parse_residue_ranges("33-35,69-70") == [33, 34, 35, 69, 70]
    assert ps.parse_residue_ranges("5, 7-9 ,5") == [5, 7, 8, 9]
    assert ps.parse_residue_ranges("") == []
    assert ps.parse_residue_ranges(" , ") == []


# ── Hull depth ───────────────────────────────────────────────────────


@pytest.mark.logic
def test_hull_depths_are_zero_outside_and_positive_inside():
    """Unsigned depth: 0 outside the hull, the distance to the nearest face
    inside it."""
    pts = np.array([[2.0, 2.0, 2.0],    # centre, 2 Å from every face
                    [2.0, 2.0, 3.5],    # 0.5 Å below the z = 4 face
                    [2.0, 2.0, 4.0],    # on the face
                    [2.0, 2.0, 9.0]])   # 5 Å outside
    got = ps.hull_depths(pts, _CUBE_EQ)
    assert got == pytest.approx([2.0, 0.5, 0.0, 0.0])


@pytest.mark.logic
def test_hull_signed_depths_are_negative_outside():
    """Signed depth: negative outside (how far out), positive inside."""
    pts = np.array([[2.0, 2.0, 2.0],
                    [2.0, 2.0, 4.0],
                    [2.0, 2.0, 9.0]])
    got = ps.hull_signed_depths(pts, _CUBE_EQ)
    assert got == pytest.approx([2.0, 0.0, -5.0])


# ── Loss breakdown, term by term ─────────────────────────────────────


@pytest.mark.logic
def test_breakdown_keys_are_the_documented_terms():
    got = _obj_breakdown(dz=10.0)
    assert set(got) == {"validity", "dist_upper", "dist_lower", "excl",
                        "interp", "pair_sc_clash", "total"}


@pytest.mark.logic
def test_breakdown_has_no_com_term():
    """Centre-of-mass separation is a diagnostic, not a penalty.

    Nothing in `_obj` references the centres of mass, so translating the
    target along the pair-restraint sphere changes the COM-COM distance
    without changing any loss term.  Both poses below put the pair residues
    exactly 4 Å apart, on opposite sides, and both score zero.
    """
    assert "com" not in _obj_breakdown(dz=10.0)
    at_plus = _obj_breakdown(dz=4.0)
    at_minus = _obj_breakdown(dz=-4.0)
    com_plus = np.linalg.norm(_CUBE.mean(0) - (_CUBE.mean(0) + [0, 0, 4.0]))
    com_minus = np.linalg.norm(_CUBE.mean(0) - (_CUBE.mean(0) + [0, 0, -4.0]))
    assert com_plus == pytest.approx(com_minus)       # symmetric by design
    assert at_plus["total"] == pytest.approx(0.0)
    assert at_minus["total"] == pytest.approx(0.0)


@pytest.mark.logic
def test_distance_term_is_quadratic_in_the_upper_overshoot():
    """Over `max_pair_distance`, the penalty is W_DIST x overshoot^2.

    The two cube centres are the only pair, so translating the target by
    (0, 0, dz) sets the pair distance to exactly dz.
    """
    assert ps.W_DIST == 1e4
    for dz, overshoot in [(10.0, 6.0), (12.0, 8.0), (4.5, 0.5)]:
        got = _obj_breakdown(dz=dz)
        assert got["dist_upper"] == pytest.approx(ps.W_DIST * overshoot ** 2)
        assert got["dist_lower"] == pytest.approx(0.0)


@pytest.mark.logic
def test_distance_term_is_zero_inside_the_window():
    """Anywhere in [min_pair_distance, max_pair_distance] costs nothing.

    Both bounds are inclusive: the penalty is max(0, d - max) and
    max(0, min - d), so landing exactly on a bound is free.
    """
    for dz in (3.5, 3.75, 4.0):
        got = _obj_breakdown(dz=dz)
        assert got["dist_upper"] == pytest.approx(0.0)
        assert got["dist_lower"] == pytest.approx(0.0)
        assert got["total"] == pytest.approx(0.0)


@pytest.mark.logic
def test_distance_term_is_quadratic_in_the_lower_shortfall():
    """Under `min_pair_distance`, the penalty is W_LOWER x shortfall^2."""
    assert ps.W_LOWER == 1e4
    got = _obj_breakdown(dz=3.4)
    assert got["dist_lower"] == pytest.approx(ps.W_LOWER * 0.1 ** 2)
    closer = _obj_breakdown(dz=3.0)
    assert closer["dist_lower"] == pytest.approx(ps.W_LOWER * 0.5 ** 2)
    assert closer["dist_lower"] > got["dist_lower"]


@pytest.mark.logic
def test_validity_term_fires_once_a_pair_residue_reaches_the_opposite_hull():
    """A pair residue must sit |VALID_TOL| Å clear of the opposite hull.

    VALID_TOL is -1.0, so the penalty is W_VALIDITY x (signed_depth + 1)
    summed over both sides, and it is still zero while a residue is more
    than 1 Å outside.

    At dz = 2.0 each centre lands exactly on the opposite hull's face
    (signed depth 0), giving 1e6 per side.  At dz = 1.0 each is 1 Å inside
    (signed depth 1), giving 2e6 per side.
    """
    assert (ps.W_VALIDITY, ps.VALID_TOL) == (1e6, -1.0)
    assert _obj_breakdown(dz=10.0)["validity"] == pytest.approx(0.0)
    assert _obj_breakdown(dz=4.0)["validity"] == pytest.approx(0.0)
    assert _obj_breakdown(dz=2.0)["validity"] == pytest.approx(2e6)
    assert _obj_breakdown(dz=1.0)["validity"] == pytest.approx(4e6)
    assert _obj_breakdown(dz=0.0)["validity"] == pytest.approx(6e6)


@pytest.mark.logic
def test_validity_term_grows_monotonically_as_the_pair_buries():
    seq = [_obj_breakdown(dz=dz)["validity"] for dz in (3.0, 2.0, 1.0, 0.0)]
    assert seq == sorted(seq)
    assert seq[0] == pytest.approx(0.0)
    assert seq[-1] > 0.0


@pytest.mark.logic
def test_exclusion_term_is_zero_while_the_pair_stays_apart():
    """An exclusion costs nothing at or beyond its stated distance."""
    far = _obj_breakdown(dz=10.0, excl=[(_CUBE_CENTRE_IDX,
                                         _CUBE_CENTRE_IDX, 4.5)])
    assert far["excl"] == pytest.approx(0.0)
    on_the_bound = _obj_breakdown(
        dz=10.0, excl=[(_CUBE_CENTRE_IDX, _CUBE_CENTRE_IDX, 10.0)])
    assert on_the_bound["excl"] == pytest.approx(0.0)


@pytest.mark.logic
def test_exclusion_term_is_quadratic_in_the_shortfall():
    """Inside the exclusion distance, the penalty is W_EXCL x shortfall^2.

    The excluded pair here is the same two cube centres, so their separation
    is dz.  Requiring 12 Å while they sit 10 Å apart is a 2 Å shortfall.
    """
    assert ps.W_EXCL == 1e4
    got = _obj_breakdown(dz=10.0,
                         excl=[(_CUBE_CENTRE_IDX, _CUBE_CENTRE_IDX, 12.0)])
    assert got["excl"] == pytest.approx(ps.W_EXCL * 2.0 ** 2)
    worse = _obj_breakdown(dz=10.0,
                           excl=[(_CUBE_CENTRE_IDX, _CUBE_CENTRE_IDX, 14.0)])
    assert worse["excl"] == pytest.approx(ps.W_EXCL * 4.0 ** 2)
    assert worse["excl"] > got["excl"]


@pytest.mark.logic
def test_interpenetration_term_is_the_summed_hull_depth_of_pair_residues():
    """Pair-only interpenetration is interp_weight x summed unsigned depth.

    At dz = 1.0 each centre is 1 Å inside the opposite hull, so the sum is
    2 Å and the term is 2 x interp_weight.
    """
    assert _obj_breakdown(dz=10.0)["interp"] == pytest.approx(0.0)
    assert _obj_breakdown(dz=1.0)["interp"] == pytest.approx(20.0)
    assert _obj_breakdown(dz=1.0, interp_weight=2.5)["interp"] == (
        pytest.approx(5.0))
    assert _obj_breakdown(dz=0.0)["interp"] == pytest.approx(40.0)


@pytest.mark.logic
def test_global_interpenetration_counts_non_pair_atoms_too():
    """--global-interp penalises every CA, not only the pair residues.

    The target here is a 2 Å cube fully swallowed by an 8 Å binder, so its
    non-pair corners are buried.  Pair-only scoring sees just the two
    centres; global scoring sees all of it, so global is strictly larger.
    """
    big = np.vstack([[[x, y, z] for x in (0.0, 8.0)
                      for y in (0.0, 8.0) for z in (0.0, 8.0)],
                     [[4.0, 4.0, 4.0]]])
    small = np.vstack([[[x, y, z] for x in (0.0, 2.0)
                        for y in (0.0, 2.0) for z in (0.0, 2.0)],
                       [[1.0, 1.0, 1.0]]])
    big_eq = ConvexHull(big).equations
    small_eq = ConvexHull(small).equations
    kwargs = dict(b_arr=big, t_arr=small, b_eq=big_eq, t_eq=small_eq, dz=3.0)
    pair_only = _obj_breakdown(global_interp=False, **kwargs)["interp"]
    everything = _obj_breakdown(global_interp=True, **kwargs)["interp"]
    assert everything > pair_only > 0.0


@pytest.mark.logic
def test_pair_sidechain_clash_term_is_quadratic_in_the_shortfall():
    """Paired sidechain heavy atoms closer than the cutoff are penalised.

    One heavy atom per side, 1.5 Å apart after the translation, against a
    2.0 Å cutoff: a 0.5 Å shortfall at weight W_PAIR_SC_CLASH.
    """
    assert ps.W_PAIR_SC_CLASH == 1e4
    got = _obj_breakdown(dz=1.5, b_sc=[[2.0, 2.0, 2.0]],
                         t_sc=[[2.0, 2.0, 2.0]],
                         sc_offsets=[(0, 1, 0, 1)], sc_cutoff=2.0)
    assert got["pair_sc_clash"] == pytest.approx(
        ps.W_PAIR_SC_CLASH * 0.5 ** 2)
    clear = _obj_breakdown(dz=10.0, b_sc=[[2.0, 2.0, 2.0]],
                           t_sc=[[2.0, 2.0, 2.0]],
                           sc_offsets=[(0, 1, 0, 1)], sc_cutoff=2.0)
    assert clear["pair_sc_clash"] == pytest.approx(0.0)


@pytest.mark.logic
def test_pair_sidechain_clash_is_disabled_at_cutoff_zero():
    """Cutoff 0 turns the term off, as --pair-sc-clash-cutoff 0 documents."""
    got = _obj_breakdown(dz=1.5, b_sc=[[2.0, 2.0, 2.0]],
                         t_sc=[[2.0, 2.0, 2.0]],
                         sc_offsets=[(0, 1, 0, 1)], sc_cutoff=0.0)
    assert got["pair_sc_clash"] == pytest.approx(0.0)


@pytest.mark.logic
def test_total_is_the_sum_of_the_terms():
    got = _obj_breakdown(dz=1.0,
                         excl=[(_CUBE_CENTRE_IDX, _CUBE_CENTRE_IDX, 4.0)])
    assert got["total"] == pytest.approx(
        got["validity"] + got["dist_upper"] + got["dist_lower"]
        + got["excl"] + got["interp"] + got["pair_sc_clash"])


@pytest.mark.logic
def test_obj_returns_a_scalar_when_components_are_not_requested():
    """The optimiser's hot path gets a float, not the breakdown dict."""
    scalar = ps._obj(
        np.r_[0.0, 0.0, 0.0, 0.0, 0.0, 10.0],
        _CUBE, _CUBE, np.array([_CUBE_CENTRE_IDX]),
        np.array([_CUBE_CENTRE_IDX]), 3.5, 4.0, _CUBE_EQ, _CUBE_EQ,
        np.array([], dtype=int), np.array([], dtype=int),
        np.array([], dtype=float), False, 10.0,
        np.empty((0, 3)), np.empty((0, 3)), [], 0.0,
    )
    assert isinstance(scalar, float)
    assert scalar == pytest.approx(_obj_breakdown(dz=10.0)["total"])


# ── PDB readers ──────────────────────────────────────────────────────


@pytest.mark.logic
def test_read_ca_keys_on_residue_number_and_filters_by_chain(tmp_path):
    """Only CA atoms of the requested chain are read, keyed by residue
    number.

    `complex_pdb` puts chain A along x at z = 0 and chain B along x at
    z = separation, so every coordinate is known.
    """
    pdb = complex_pdb(tmp_path / "c.pdb", n_receptor=4, n_effector=3,
                      spacing=3.8, separation=5.0)
    a = ps.read_ca(pdb, "A")
    b = ps.read_ca(pdb, "B")
    assert sorted(a) == [1, 2, 3, 4]
    assert sorted(b) == [1, 2, 3]
    assert a[1] == pytest.approx([0.0, 0.0, 0.0])
    assert a[4] == pytest.approx([3 * 3.8, 0.0, 0.0])
    assert b[1] == pytest.approx([0.0, 0.0, 5.0])
    assert ps.read_ca(pdb, "Z") == {}


@pytest.mark.logic
def test_read_heavy_groups_every_non_hydrogen_atom_by_residue(tmp_path):
    """read_heavy keeps backbone atoms; read_sidechain_heavy drops them.

    The residue below carries N, CA, C, O and CB plus one hydrogen.  Heavy
    reading keeps five atoms, sidechain-only reading keeps the CB alone, and
    the hydrogen is dropped by both.
    """
    pdb = tmp_path / "res.pdb"
    pdb.write_text("\n".join([
        "ATOM      1  N   ALA A   7       0.000   0.000   0.000  1.00  0.00           N",
        "ATOM      2  CA  ALA A   7       1.000   0.000   0.000  1.00  0.00           C",
        "ATOM      3  C   ALA A   7       2.000   0.000   0.000  1.00  0.00           C",
        "ATOM      4  O   ALA A   7       3.000   0.000   0.000  1.00  0.00           O",
        "ATOM      5  CB  ALA A   7       1.000   1.500   0.000  1.00  0.00           C",
        "ATOM      6  HA  ALA A   7       1.000  -1.000   0.000  1.00  0.00           H",
        "ATOM      7  CB  GLY B   7       9.000   9.000   9.000  1.00  0.00           C",
        "END",
    ]) + "\n")
    heavy = ps.read_heavy(pdb, "A")
    assert list(heavy) == [7]
    assert heavy[7].shape == (5, 3)
    side = ps.read_sidechain_heavy(pdb, "A")
    assert list(side) == [7]
    assert side[7] == pytest.approx(np.array([[1.0, 1.5, 0.0]]))
    assert ps.read_sidechain_heavy(pdb, "B")[7] == pytest.approx(
        np.array([[9.0, 9.0, 9.0]]))


@pytest.mark.logic
def test_read_sidechain_heavy_is_empty_for_a_ca_only_structure(tmp_path):
    """A CA-only fixture has no sidechain atoms, which is what makes the
    solver report every pair as `no sidechain (likely GLY)`."""
    pdb = lattice_pdb(tmp_path / "ca.pdb", chain="A", shape=(2, 2, 2))
    assert ps.read_sidechain_heavy(pdb, "A") == {}
    assert len(ps.read_heavy(pdb, "A")) == 8


# ── Posed-PDB writing ────────────────────────────────────────────────


@pytest.mark.logic
def test_write_posed_pdb_applies_the_transform_to_the_target_only(tmp_path):
    """The binder is copied verbatim and the target gets R x + t.

    Serials are renumbered from 1 across the whole file, and the target
    keeps its own chain id.
    """
    binder = lattice_pdb(tmp_path / "b.pdb", chain="A", shape=(2, 2, 2),
                         step=4.0)
    target = lattice_pdb(tmp_path / "t.pdb", chain="B", shape=(2, 2, 2),
                         step=4.0)
    out = tmp_path / "posed.pdb"
    shift = np.array([0.0, 0.0, 20.0])
    ps.write_posed_pdb(out, binder, target, "A", "B", np.eye(3), shift)

    posed_a = ps.read_ca(out, "A")
    posed_b = ps.read_ca(out, "B")
    orig_b = ps.read_ca(target, "B")
    assert posed_a[1] == pytest.approx(ps.read_ca(binder, "A")[1])
    for rn, xyz in orig_b.items():
        assert posed_b[rn] == pytest.approx(xyz + shift)

    serials = [int(ln[6:11]) for ln in out.read_text().splitlines()
               if ln.startswith("ATOM")]
    assert serials == list(range(1, 17))


@pytest.mark.logic
def test_write_posed_pdb_applies_rotation(tmp_path):
    """A 90 degree rotation about z maps (x, y, z) to (-y, x, z)."""
    binder = lattice_pdb(tmp_path / "b.pdb", chain="A", shape=(2, 2, 2))
    target = lattice_pdb(tmp_path / "t.pdb", chain="B", shape=(2, 2, 2))
    R = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    out = tmp_path / "posed.pdb"
    ps.write_posed_pdb(out, binder, target, "A", "B", R, np.zeros(3))
    orig = ps.read_ca(target, "B")
    posed = ps.read_ca(out, "B")
    for rn, (x, y, z) in orig.items():
        assert posed[rn] == pytest.approx([-y, x, z])


# ── Clash and contact counting ───────────────────────────────────────


@pytest.mark.logic
def test_analyse_contacts_counts_pairs_inside_the_cutoff(tmp_path):
    """Two parallel 5-residue chains 2 Å apart, CA only, identity transform.

    Residue i of the receptor sits at (3.8 i, 0, 0) and residue j of the
    effector at (3.8 j, 0, 2), so their separation is
    sqrt(2^2 + (3.8 (i - j))^2):

        |i - j| = 0 -> 2.000 Å
        |i - j| = 1 -> 4.294 Å
        |i - j| = 2 -> 7.861 Å
        |i - j| = 3 -> 11.571 Å

    At an 8 Å contact cutoff the pairs with |i - j| <= 2 count, which is
    5 + 2x4 + 2x3 = 19 pairs.
    """
    pdb = complex_pdb(tmp_path / "c.pdb", n_receptor=5, n_effector=5,
                      spacing=3.8, separation=2.0)
    heavy_a = ps.read_heavy(pdb, "A")
    heavy_b = ps.read_heavy(pdb, "B")
    contacts = ps.analyse_contacts(heavy_a, heavy_b, np.eye(3), np.zeros(3),
                                   contact_cutoff=8.0)
    assert len(contacts) == 19
    assert contacts[(1, 1)] == pytest.approx(2.0)
    assert contacts[(1, 2)] == pytest.approx(np.hypot(2.0, 3.8))
    assert (1, 4) not in contacts


@pytest.mark.logic
def test_clash_count_is_the_number_of_pairs_under_the_clash_cutoff(tmp_path):
    """Only the five directly-stacked pairs are within 2.5 Å.

    Same geometry as above.  The clash set is the subset of the contact
    dictionary below `clash_cutoff`, which main() computes as
    `{k for k, d in contacts.items() if d < clash_cutoff}`.
    """
    pdb = complex_pdb(tmp_path / "c.pdb", n_receptor=5, n_effector=5,
                      spacing=3.8, separation=2.0)
    contacts = ps.analyse_contacts(ps.read_heavy(pdb, "A"),
                                   ps.read_heavy(pdb, "B"),
                                   np.eye(3), np.zeros(3), contact_cutoff=8.0)
    clashes = {k for k, d in contacts.items() if d < 2.5}
    assert clashes == {(i, i) for i in range(1, 6)}
    assert len(clashes) == 5
    # Raising the bar past the stacked distance catches the next shell too.
    assert len({k for k, d in contacts.items() if d < 5.0}) == 5 + 8


@pytest.mark.logic
def test_analyse_contacts_sees_no_clashes_once_the_chains_separate(tmp_path):
    pdb = complex_pdb(tmp_path / "c.pdb", n_receptor=5, n_effector=5,
                      spacing=3.8, separation=30.0)
    contacts = ps.analyse_contacts(ps.read_heavy(pdb, "A"),
                                   ps.read_heavy(pdb, "B"),
                                   np.eye(3), np.zeros(3), contact_cutoff=8.0)
    assert contacts == {}


@pytest.mark.logic
def test_analyse_contacts_applies_the_transform_before_measuring(tmp_path):
    """The same structures 30 Å apart become contacts once translated back."""
    pdb = complex_pdb(tmp_path / "c.pdb", n_receptor=3, n_effector=3,
                      spacing=3.8, separation=30.0)
    pulled_in = ps.analyse_contacts(
        ps.read_heavy(pdb, "A"), ps.read_heavy(pdb, "B"),
        np.eye(3), np.array([0.0, 0.0, -28.0]), contact_cutoff=8.0)
    assert pulled_in[(1, 1)] == pytest.approx(2.0)


# ── Contig construction ──────────────────────────────────────────────


@pytest.mark.logic
def test_expand_with_gap_join_closes_only_interior_gaps_within_the_limit():
    """A gap is joined when it is flanked on both sides and short enough.

    Leading and trailing stretches have no flank on one side, so they are
    never absorbed however small.
    """
    all_rn = list(range(1, 11))
    assert ps.expand_with_gap_join({3, 5}, all_rn, 0) == {3, 5}
    assert ps.expand_with_gap_join({3, 5}, all_rn, 1) == {3, 4, 5}
    assert ps.expand_with_gap_join({3, 6}, all_rn, 1) == {3, 6}
    assert ps.expand_with_gap_join({3, 6}, all_rn, 2) == {3, 4, 5, 6}
    # Residue 1 is not flanked below, residue 10 not flanked above.
    assert ps.expand_with_gap_join({3}, all_rn, 5) == {3}


@pytest.mark.logic
def test_make_contig_is_slash_separated_with_chain_prefixed_fixed_runs():
    """The contig format: `/`-separated parts, then a space and the target
    chain.

    A fixed run is written chain-prefixed as `A<first>-<last>`.  A de novo
    run is written as a `<len>-<len>` length range with no chain prefix.
    Commas never appear.  Residues 5-8 of a 1-12 binder are designed, so the
    string is `A1-4/4-4/A9-12 B`.
    """
    all_rn = list(range(1, 13))
    got = ps.make_contig(all_rn, {5, 6, 7, 8}, "A", "B", join_gap=0)
    assert got == "A1-4/4-4/A9-12 B"
    assert "," not in got


@pytest.mark.logic
def test_make_contig_handles_a_design_region_at_each_end():
    all_rn = list(range(1, 9))
    assert ps.make_contig(all_rn, {1, 2}, "A", "B", 0) == "2-2/A3-8 B"
    assert ps.make_contig(all_rn, {7, 8}, "A", "B", 0) == "A1-6/2-2 B"
    assert ps.make_contig(all_rn, set(), "A", "B", 0) == "A1-8 B"
    assert ps.make_contig(all_rn, set(all_rn), "A", "B", 0) == "8-8 B"


@pytest.mark.logic
def test_make_contig_gap_joining_merges_two_segments_into_one():
    """Joining a 1-residue gap turns two de novo runs into a single longer
    one, which is the whole point of the A/B/C variants main() prints."""
    all_rn = list(range(1, 13))
    assert ps.make_contig(all_rn, {5, 7}, "A", "B", 0) == (
        "A1-4/1-1/A6-6/1-1/A8-12 B")
    assert ps.make_contig(all_rn, {5, 7}, "A", "B", 1) == (
        "A1-4/3-3/A8-12 B")


@pytest.mark.logic
def test_make_contig_respects_the_target_chain_id():
    assert ps.make_contig([1, 2, 3], {2}, "X", "Y", 0) == "X1-1/1-1/X3-3 Y"


@pytest.mark.logic
def test_describe_design_region_lists_runs_and_counts_residues():
    """Runs are comma-separated ranges, singletons printed bare."""
    all_rn = list(range(1, 13))
    assert ps.describe_design_region(all_rn, {5, 6, 7, 8}, 0) == ("5-8", 4)
    assert ps.describe_design_region(all_rn, {3, 7, 8}, 0) == ("3, 7-8", 3)
    assert ps.describe_design_region(all_rn, set(), 0) == ("(none)", 0)
    assert ps.describe_design_region(all_rn, {5, 7}, 1) == ("5-7", 3)
    # A design position outside the binder's own numbering is ignored.
    assert ps.describe_design_region(all_rn, {99}, 0) == ("(none)", 0)


# ── The solver itself ────────────────────────────────────────────────


@pytest.mark.logic
def test_solve_brings_each_restrained_pair_inside_the_max_pair_distance(
        tmp_path):
    """A pair restraint must end with those two CA atoms close enough.

    Both monomers are 3x3x3 lattices filling [0, 8]^3 in their own frame.
    Three restraints pin the binder's +z face to the target's -z face.  The
    restraint is satisfiable by a pure translation of (0, 0, 8 + d), so a
    zero loss is reachable, and the distance is re-measured here from the
    returned rotation and translation rather than trusted from
    `pair_results`.
    """
    binder, target = _two_lattices(tmp_path)
    binder_ca = ps.read_ca(binder, "A")
    target_ca = ps.read_ca(target, "B")
    pairs = [(_rn(0, 0, 2), _rn(0, 0, 0), 4.0),
             (_rn(2, 0, 2), _rn(2, 0, 0), 4.0),
             (_rn(0, 2, 2), _rn(0, 2, 0), 4.0)]

    R, t, loss, pair_results, history, breakdown = ps.solve(
        binder_ca, target_ca, pairs, [], n_restarts=5,
        min_pair_dist=3.5, max_pair_dist=4.0)

    assert loss == pytest.approx(0.0, abs=1e-6)
    assert breakdown["total"] == pytest.approx(0.0, abs=1e-6)
    assert len(history) == 5
    assert all(pr["satisfied"] for pr in pair_results)

    b_rn = sorted(binder_ca)
    t_rn = sorted(target_ca)
    moved = {rn: R @ target_ca[rn] + t for rn in t_rn}
    for rb, rt, _ in pairs:
        measured = float(np.linalg.norm(binder_ca[rb] - moved[rt]))
        assert 3.5 - 0.01 <= measured <= 4.0 + 0.01
        reported = next(pr["achieved_d"] for pr in pair_results
                        if pr["binder_res"] == rb and pr["target_res"] == rt)
        assert measured == pytest.approx(reported, abs=1e-3)
    assert len(b_rn) == len(t_rn) == 27


@pytest.mark.logic
def test_solve_keeps_an_exclusion_pair_at_least_its_distance_apart(tmp_path):
    """An exclusion must end up no closer than its stated distance.

    One restraint holds target residue `_rn(0, 0, 0)` 3.5-4 Å off the
    binder's (0, 0, 8) corner.  The exclusion additionally requires the
    binder's top-face centre (4, 4, 8) to stay at least 10 Å from target
    residue `_rn(0, 0, 1)`, which the restraint alone would allow much
    closer.  Both constraints are then checked on the returned pose.
    """
    binder, target = _two_lattices(tmp_path)
    binder_ca = ps.read_ca(binder, "A")
    target_ca = ps.read_ca(target, "B")
    pairs = [(_rn(0, 0, 2), _rn(0, 0, 0), 4.0)]
    excl_b, excl_t, excl_d = _rn(1, 1, 2), _rn(0, 0, 1), 10.0

    R, t, loss, pair_results, _, breakdown = ps.solve(
        binder_ca, target_ca, pairs, [(excl_b, excl_t, excl_d)],
        n_restarts=8, min_pair_dist=3.5, max_pair_dist=4.0)

    assert breakdown["excl"] == pytest.approx(0.0, abs=1e-6)
    assert pair_results[0]["satisfied"]
    moved_excl = R @ target_ca[excl_t] + t
    separation = float(np.linalg.norm(binder_ca[excl_b] - moved_excl))
    assert separation >= excl_d - 0.01


@pytest.mark.logic
def test_solve_skips_pairs_whose_residues_are_absent(tmp_path, capsys):
    """A restraint naming a residue neither chain has is dropped with a
    warning, and the remaining restraints still solve."""
    binder, target = _two_lattices(tmp_path)
    binder_ca = ps.read_ca(binder, "A")
    target_ca = ps.read_ca(target, "B")
    pairs = [(_rn(0, 0, 2), _rn(0, 0, 0), 4.0), (9999, 9999, 4.0)]
    _, _, _, pair_results, _, _ = ps.solve(
        binder_ca, target_ca, pairs, [], n_restarts=2,
        min_pair_dist=3.5, max_pair_dist=4.0)
    assert [pr["binder_res"] for pr in pair_results] == [_rn(0, 0, 2)]
    assert "residues not found" in capsys.readouterr().out


@pytest.mark.logic
def test_solve_exits_when_no_restraint_names_a_real_residue(tmp_path):
    binder, target = _two_lattices(tmp_path)
    with pytest.raises(SystemExit) as exc:
        ps.solve(ps.read_ca(binder, "A"), ps.read_ca(target, "B"),
                 [(9999, 9999, 4.0)], [], n_restarts=1,
                 min_pair_dist=3.5, max_pair_dist=4.0)
    assert "no valid pairs" in str(exc.value)


@pytest.mark.logic
def test_solve_with_differential_evolution_returns_one_history_row(tmp_path):
    """The DE path reports a single final value, not a per-restart curve."""
    binder, target = _two_lattices(tmp_path)
    _, _, loss, pair_results, history, _ = ps.solve(
        ps.read_ca(binder, "A"), ps.read_ca(target, "B"),
        [(_rn(0, 0, 2), _rn(0, 0, 0), 4.0)], [], n_restarts=30,
        min_pair_dist=3.5, max_pair_dist=4.0, use_de=True)
    assert len(history) == 1
    assert history[0][0] == 0
    assert loss == pytest.approx(0.0, abs=1e-3)
    assert pair_results[0]["satisfied"]


@pytest.mark.logic
def test_solve_reports_a_pair_with_no_sidechain_as_skipped(tmp_path, capsys):
    """CA-only inputs leave the sidechain-clash term with nothing to do.

    The solver says so rather than silently scoring zero, because a real
    GLY pair and a stripped structure look the same to it.
    """
    binder, target = _two_lattices(tmp_path)
    binder_ca = ps.read_ca(binder, "A")
    target_ca = ps.read_ca(target, "B")
    _, _, _, _, _, breakdown = ps.solve(
        binder_ca, target_ca, [(_rn(0, 0, 2), _rn(0, 0, 0), 4.0)], [],
        n_restarts=2, min_pair_dist=3.5, max_pair_dist=4.0,
        binder_sc=ps.read_sidechain_heavy(binder, "A"),
        target_sc=ps.read_sidechain_heavy(target, "B"),
        pair_sc_clash_cutoff=2.0)
    out = capsys.readouterr().out
    assert "no sidechain (likely GLY)" in out
    assert "active pairs: 0/1" in out
    assert breakdown["pair_sc_clash"] == pytest.approx(0.0)


@pytest.mark.logic
def test_solve_with_global_interp_still_satisfies_the_restraint(tmp_path):
    binder, target = _two_lattices(tmp_path)
    _, _, _, pair_results, _, _ = ps.solve(
        ps.read_ca(binder, "A"), ps.read_ca(target, "B"),
        [(_rn(0, 0, 2), _rn(0, 0, 0), 4.0)], [], n_restarts=4,
        min_pair_dist=3.5, max_pair_dist=4.0,
        global_interp=True, interp_weight=1.0)
    assert pair_results[0]["satisfied"]


# ── Heatmap ──────────────────────────────────────────────────────────


@pytest.mark.logic
def test_make_heatmap_writes_a_png(tmp_path):
    pdb = complex_pdb(tmp_path / "c.pdb", n_receptor=6, n_effector=6,
                      spacing=3.8, separation=2.0)
    contacts = ps.analyse_contacts(ps.read_heavy(pdb, "A"),
                                   ps.read_heavy(pdb, "B"),
                                   np.eye(3), np.zeros(3), contact_cutoff=8.0)
    pair_results = [{"binder_res": 1, "target_res": 1, "achieved_d": 2.0,
                     "min_d": 3.5, "max_d": 4.0, "satisfied": False}]
    out = tmp_path / "heat.png"
    ps.make_heatmap(contacts, [(1, 1, 4.0)], [2, 3], {1, 2}, {1},
                    "A", "B", 2.5, 8.0, pair_results, out)
    assert out.exists() and out.stat().st_size > 0


@pytest.mark.logic
def test_make_heatmap_renders_a_placeholder_when_there_are_no_contacts(
        tmp_path):
    out = tmp_path / "empty.png"
    ps.make_heatmap({}, [], [], set(), set(), "A", "B", 2.0, 8.0, [], out)
    assert out.exists() and out.stat().st_size > 0


# ── CLI ──────────────────────────────────────────────────────────────


@pytest.mark.logic
def test_parse_args_defaults_match_the_documented_window(monkeypatch):
    monkeypatch.setattr(sys, "argv", [
        "pose_solver.py", "--binder", "b.pdb", "--target", "t.pdb",
        "--pairs", "A1-B2"])
    args = ps.parse_args()
    assert (args.min_pair_distance, args.max_pair_distance) == (3.5, 4.0)
    assert (args.binder_chain, args.target_chain) == ("A", "B")
    assert args.clash_cutoff == 2.0
    assert args.contact_cutoff == 8.0
    assert args.n_restarts == 300
    assert args.interp_weight == 10.0
    assert args.pair_sc_clash_cutoff == ps.PAIR_SC_CLASH_TOL
    assert args.global_interp is False
    assert args.use_de is False
    assert args.output_prefix == "solved_pose"


@pytest.mark.logic
def test_main_writes_every_output_and_the_restraints_are_satisfied(
        tmp_path, monkeypatch):
    """End-to-end: four artefacts, a satisfied restraint set, a contig.

    main() writes relative to the working directory, so the test runs in
    tmp_path.  The clash count is zero because the solved pose holds the
    pair CA atoms just under 4 Å apart, which is above the 2 Å clash
    cutoff.
    """
    binder, target = _two_lattices(tmp_path)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", [
        "pose_solver.py",
        "--binder", str(binder), "--target", str(target),
        "--binder-chain", "a", "--target-chain", "b",
        "--pairs", (f"A{_rn(0, 0, 2)}-B{_rn(0, 0, 0)} "
                    f"A{_rn(2, 0, 2)}-B{_rn(2, 0, 0)} "
                    f"A{_rn(0, 2, 2)}-B{_rn(0, 2, 0)}"),
        "--n-restarts", "3",
        "--contig-design-region", "5-8",
        "--output-prefix", "out",
    ])
    ps.main()

    for name in ("out_posed.pdb", "out_results.json",
                 "out_restart_losses.csv", "out_heatmap.png"):
        assert (tmp_path / name).exists(), name

    results = json.loads((tmp_path / "out_results.json").read_text())
    assert results["all_pairs_satisfied"] is True
    assert len(results["pair_results"]) == 3
    for pr in results["pair_results"]:
        assert 3.5 - 0.01 <= pr["achieved_d"] <= 4.0 + 0.01
    assert results["final_loss"] == pytest.approx(0.0, abs=1e-3)
    assert results["n_clash_pairs"] == 0
    assert results["clash_cutoff_A"] == 2.0
    assert results["contig_design_region_input"] == [5, 6, 7, 8]
    assert set(results["loss_breakdown"]) == {
        "validity", "dist_upper", "dist_lower", "excl", "interp",
        "pair_sc_clash", "total"}
    assert results["loss_weights"]["W_VALIDITY"] == 1e6
    # The de novo span is the 4 designed residues; everything else is fixed
    # and chain-prefixed.  '/' separates, never a comma.
    assert results["contigs"]["A) No joining"]["contig"] == (
        "A1-4/4-4/A9-27 B")
    assert results["contigs"]["A) No joining"]["n_residues"] == 4

    # Chain ids are upper-cased from the lowercase flags, and the posed
    # complex carries both chains.
    posed = tmp_path / "out_posed.pdb"
    assert len(ps.read_ca(posed, "A")) == 27
    assert len(ps.read_ca(posed, "B")) == 27

    csv_rows = (tmp_path / "out_restart_losses.csv").read_text().splitlines()
    assert csv_rows[0] == "restart_idx,restart_loss,best_loss_so_far"
    assert len(csv_rows) == 1 + 3


@pytest.mark.logic
def test_main_rejects_identical_binder_and_target_chains(
        tmp_path, monkeypatch):
    """Two same-named chains in one output file would be unusable, so the
    run stops before any work."""
    binder, target = _two_lattices(tmp_path)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", [
        "pose_solver.py", "--binder", str(binder), "--target", str(target),
        "--binder-chain", "A", "--target-chain", "a",
        "--pairs", "A1-B1", "--n-restarts", "1"])
    with pytest.raises(SystemExit) as exc:
        ps.main()
    assert "both 'A'" in str(exc.value)


@pytest.mark.logic
def test_main_accepts_an_exclusion_and_records_the_weights(
        tmp_path, monkeypatch):
    """--exclusions defaults each token to 4.5 Å when no `@DIST` is given."""
    binder, target = _two_lattices(tmp_path)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", [
        "pose_solver.py", "--binder", str(binder), "--target", str(target),
        "--pairs", f"A{_rn(0, 0, 2)}-B{_rn(0, 0, 0)}",
        "--exclusions", f"A{_rn(1, 1, 2)}-B{_rn(0, 0, 1)}@10.0",
        "--n-restarts", "3", "--interp-weight", "1.0",
        "--pair-sc-clash-cutoff", "0",
        "--output-prefix", "excl"])
    ps.main()
    results = json.loads((tmp_path / "excl_results.json").read_text())
    assert results["loss_breakdown"]["excl"] == pytest.approx(0.0, abs=1e-6)
    assert results["loss_weights"]["W_INTERP"] == 1.0
    assert ps.parse_pairs("A8-B33", 4.5) == [(8, 33, 4.5)]
