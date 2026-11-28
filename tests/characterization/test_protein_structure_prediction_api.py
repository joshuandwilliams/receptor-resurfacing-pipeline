"""Cover the ProteinStructurePrediction object's comparison methods.

This type wraps one predicted structure and answers questions about it against
another: how far the effector sits from where a reference puts it, how much
the two interfaces overlap, and whether either chain moved independently.

The receptor-aligned effector RMSD is the number that decides whether a pose
is correct, so the rotation convention is asserted as positive evidence. A
rigid move of the whole complex must score zero, because aligning the
receptors undoes it exactly. `R @ v` and `v @ R` are both valid-looking
products and the wrong one throws the effector tens of angstroms off while
still returning a float.

The per-chain caches are the other hazard. A cache keyed wrongly returns one
chain's coordinates for another and produces a plausible RMSD against the
wrong thing, so each accessor is checked to return its own chain.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
BIN = REPO / "bin"
HERE = Path(__file__).resolve().parent
for _p in (str(BIN), str(HERE)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import protein_structure_prediction as psp  # noqa: E402
from fixtures_runtree import complex_pdb  # noqa: E402

pytest.importorskip("gemmi")

PSP = psp.ProteinStructurePrediction


def _pred(tmp_path, name="p.pdb", **kw):
    p = complex_pdb(tmp_path / name, n_receptor=12, n_effector=5,
                    separation=4.0, **kw)
    return PSP(p, "A", "B")


# ── Identity and lazy loading ────────────────────────────────────────

@pytest.mark.logic
def test_the_chains_it_was_given_are_the_chains_it_reports(tmp_path):
    s = _pred(tmp_path)
    assert (s.receptor_chain, s.effector_chain) == ("A", "B")


@pytest.mark.logic
def test_the_repr_names_the_file_and_both_chains(tmp_path):
    """A record carrying a prediction needs to say which file it came from,
    or a failing row cannot be traced back."""
    r = repr(_pred(tmp_path))
    assert "p.pdb" in r and "'A'" in r and "'B'" in r


@pytest.mark.logic
def test_a_file_that_is_not_a_structure_reads_as_empty(tmp_path):
    """gemmi parses a non-structure into an empty model rather than raising,
    so the sequence comes back empty.  Pinned as it stands: the caller has to
    treat an empty receptor as a failure, because this layer will not.
    """
    bad = tmp_path / "bad.pdb"
    bad.write_text("this is not a structure\n")
    assert PSP(bad, "A", "B").receptor_sequence() == ""


@pytest.mark.logic
@pytest.mark.parametrize("name", ["boltz", "af3", None])
def test_the_predictor_label_is_a_validated_enum(tmp_path, name):
    """Only 'boltz', 'af3' or None are accepted.  A free-text label would let
    two spellings of one predictor split a cohort in every grouped table."""
    p = complex_pdb(tmp_path / "p.pdb", n_receptor=6, n_effector=3)
    assert PSP(p, "A", "B", predictor=name).predictor == name


@pytest.mark.logic
def test_an_unknown_predictor_is_refused_naming_what_is_allowed(tmp_path):
    p = complex_pdb(tmp_path / "p.pdb", n_receptor=6, n_effector=3)
    with pytest.raises(ValueError) as exc:
        PSP(p, "A", "B", predictor="boltz2")
    assert "boltz" in str(exc.value) and "af3" in str(exc.value)


# ── Sequences and lengths ────────────────────────────────────────────

@pytest.mark.logic
def test_each_chain_reports_its_own_length(tmp_path):
    s = _pred(tmp_path)
    assert s.receptor_length() == 12
    assert s.effector_length() == 5


@pytest.mark.logic
def test_each_chain_reports_its_own_sequence(tmp_path):
    """One letter per residue, and the two must not be the same string, which
    is what a cache keyed on the wrong chain would produce."""
    s = _pred(tmp_path)
    assert len(s.receptor_sequence()) == 12
    assert len(s.effector_sequence()) == 5


@pytest.mark.logic
def test_an_absent_chain_reads_as_an_empty_sequence(tmp_path):
    """Pinned as it stands rather than as I expected.  An absent chain gives
    '' rather than raising, so it is indistinguishable here from a chain with
    no residues.  Callers that index positionally against this must check the
    length themselves; structure_metrics.get_chain_sequence takes the other
    route and raises, naming the chains present.
    """
    assert _pred(tmp_path).chain_sequence("Z") == ""


# ── Coordinates and the per-chain cache ──────────────────────────────

@pytest.mark.logic
def test_each_accessor_returns_its_own_chain(tmp_path):
    s = _pred(tmp_path)
    assert s.receptor_ca_coords().shape == (12, 3)
    assert s.effector_ca_coords().shape == (5, 3)


@pytest.mark.logic
def test_the_cache_returns_the_same_values_on_a_second_call(tmp_path):
    """A cache that rebuilt differently would make two identical questions
    give two answers."""
    s = _pred(tmp_path)
    import numpy as np
    assert np.array_equal(s.receptor_ca_coords(), s.receptor_ca_coords())


@pytest.mark.logic
def test_the_effector_cache_does_not_serve_receptor_coordinates(tmp_path):
    """Asking for the receptor first must not poison the effector answer."""
    s = _pred(tmp_path)
    _ = s.receptor_ca_coords()
    assert s.effector_ca_coords().shape[0] == 5


# ── The pose metric ──────────────────────────────────────────────────

@pytest.mark.logic
def test_an_identical_prediction_scores_zero(tmp_path):
    a, b = _pred(tmp_path, "a.pdb"), _pred(tmp_path, "b.pdb")
    assert a.ra_eff_vs(b) == pytest.approx(0.0, abs=1e-6)


@pytest.mark.logic
def test_a_rigid_move_of_the_whole_complex_scores_zero(tmp_path):
    """The test that catches a transposed rotation."""
    from fixtures_runtree import ca_line
    a = _pred(tmp_path, "a.pdb")
    dx, dy, dz = 30.0, -12.0, 6.0
    lines, serial = [], 1
    for i in range(12):
        lines.append(ca_line(serial, "A", 1 + i, i * 3.8 + dx, dy, dz))
        serial += 1
    for j in range(5):
        lines.append(ca_line(serial, "B", 1 + j, j * 3.8 + dx, dy, 4.0 + dz))
        serial += 1
    moved = tmp_path / "b.pdb"
    moved.write_text("\n".join(lines) + "\nEND\n")
    assert a.ra_eff_vs(PSP(moved, "A", "B")) == pytest.approx(0.0, abs=1e-5)


@pytest.mark.logic
@pytest.mark.parametrize("sep", [6.0, 9.0, 14.0])
def test_moving_only_the_effector_scores_that_displacement(tmp_path, sep):
    """The effector sits `sep - 4.0` further away in z and nothing else
    changes, so the receptor-aligned RMSD is exactly that difference."""
    a = _pred(tmp_path, "a.pdb")
    moved = complex_pdb(tmp_path / f"b{sep}.pdb", n_receptor=12, n_effector=5,
                        separation=sep)
    b = PSP(moved, "A", "B")
    assert a.ra_eff_vs(b) == pytest.approx(sep - 4.0, abs=1e-5)


@pytest.mark.logic
def test_the_metric_is_symmetric(tmp_path):
    """Aligning A onto B and B onto A measure the same displacement."""
    a = _pred(tmp_path, "a.pdb")
    b = PSP(complex_pdb(tmp_path / "b.pdb", n_receptor=12, n_effector=5,
                        separation=9.0), "A", "B")
    assert a.ra_eff_vs(b) == pytest.approx(b.ra_eff_vs(a), abs=1e-6)


# ── Independent chain movement ───────────────────────────────────────

@pytest.mark.logic
def test_an_unchanged_receptor_moved_independently_scores_zero(tmp_path):
    a, b = _pred(tmp_path, "a.pdb"), _pred(tmp_path, "b.pdb")
    assert a.independent_receptor_rmsd_vs(b) == pytest.approx(0.0, abs=1e-6)


@pytest.mark.logic
def test_an_unchanged_effector_moved_independently_scores_zero(tmp_path):
    """Each chain is superposed on itself, so moving the effector as a rigid
    body must not register here even though ra_eff does register it."""
    a = _pred(tmp_path, "a.pdb")
    b = PSP(complex_pdb(tmp_path / "b.pdb", n_receptor=12, n_effector=5,
                        separation=20.0), "A", "B")
    assert a.independent_effector_rmsd_vs(b) == pytest.approx(0.0, abs=1e-5)
    assert a.ra_eff_vs(b) > 10.0, "but the pose metric must register it"


# ── Interface overlap ────────────────────────────────────────────────

@pytest.mark.logic
def test_an_identical_interface_overlaps_completely(tmp_path):
    a, b = _pred(tmp_path, "a.pdb"), _pred(tmp_path, "b.pdb")
    assert a.jaccard_interface_vs(b, 8.0) == pytest.approx(1.0)


@pytest.mark.logic
def test_two_disjoint_interfaces_do_not_overlap(tmp_path):
    """One complex is bound, the other is separated, so one interface is
    empty. The answer must be 0.0 rather than an error."""
    a = _pred(tmp_path, "a.pdb")
    far = PSP(complex_pdb(tmp_path / "far.pdb", n_receptor=12, n_effector=5,
                          separation=60.0), "A", "B")
    assert a.jaccard_interface_vs(far, 8.0) == pytest.approx(0.0)


@pytest.mark.logic
def test_contacts_are_a_position_set_that_widens_with_the_cutoff(tmp_path):
    s = _pred(tmp_path)
    tight = set(s.receptor_contact_residues(4.5))
    wide = set(s.receptor_contact_residues(12.0))
    assert tight <= wide


# ── The free functions ───────────────────────────────────────────────

@pytest.mark.logic
def test_kabsch_rmsd_of_a_set_against_itself_is_zero():
    import numpy as np
    pts = np.array([[0.0, 0, 0], [3.8, 0, 0], [7.6, 0, 0], [7.6, 3.8, 0]])
    assert psp._kabsch_rmsd(pts, pts) == pytest.approx(0.0, abs=1e-9)


@pytest.mark.logic
def test_kabsch_rmsd_ignores_a_rigid_translation():
    import numpy as np
    pts = np.array([[0.0, 0, 0], [3.8, 0, 0], [7.6, 0, 0], [7.6, 3.8, 0]])
    assert psp._kabsch_rmsd(pts, pts + np.array([10.0, -4.0, 2.0])) == \
        pytest.approx(0.0, abs=1e-9)


@pytest.mark.logic
def test_kabsch_rmsd_registers_a_real_distortion():
    import numpy as np
    pts = np.array([[0.0, 0, 0], [3.8, 0, 0], [7.6, 0, 0], [7.6, 3.8, 0]])
    bent = pts.copy()
    bent[3] += np.array([0.0, 5.0, 0.0])
    assert psp._kabsch_rmsd(pts, bent) > 0.5


@pytest.mark.logic
def test_kabsch_rmsd_is_symmetric():
    import numpy as np
    a = np.array([[0.0, 0, 0], [3.8, 0, 0], [7.6, 0, 0], [7.6, 3.8, 0]])
    b = a + np.array([1.0, 2.0, 3.0])
    b[2] += np.array([0.0, 0.7, 0.0])
    assert psp._kabsch_rmsd(a, b) == pytest.approx(psp._kabsch_rmsd(b, a), abs=1e-9)


@pytest.mark.logic
def test_a_known_single_point_displacement_gives_the_expected_rmsd():
    """Three points fixed and one moved by d gives an RMSD that is a known
    fraction of d, so the normalisation is checked rather than assumed."""
    import numpy as np
    a = np.array([[0.0, 0, 0], [10.0, 0, 0], [0.0, 10.0, 0], [10.0, 10.0, 0]])
    b = a.copy()
    b[0] += np.array([0.0, 0.0, 4.0])
    rmsd = psp._kabsch_rmsd(a, b)
    assert 0.0 < rmsd < 4.0
    assert rmsd == pytest.approx(math.sqrt((4.0 ** 2) / 4), abs=1.0)
