"""Unit tests for bin/protein_structure_prediction.py (Phase 4 Tier 2).

Construction-time + geometry-helper coverage runs without gemmi.

The tests against a committed Boltz output skip when that fixture is absent,
which it is in this repository, so the second half of this file drives the
same accessors on structures built by tests/characterization/fixtures_runtree
at a geometry stated at the call site. A built structure states the answer
first and then constructs the input that must produce it, which a blessed
real output cannot do.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

_REPO_ROOT = Path(__file__).resolve().parents[2]
_HERE = Path(__file__).resolve().parent
for _p in (str(_REPO_ROOT / "bin"), str(_HERE)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import protein_structure_prediction as psp  # noqa: E402
from fixtures_runtree import ca_line, complex_pdb  # noqa: E402


@pytest.mark.logic
class TestConstruction:
    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            psp.ProteinStructurePrediction(
                path=tmp_path / "nope.pdb",
                receptor_chain="A",
                effector_chain="B",
            )

    def test_empty_receptor_chain_raises(self, tmp_path):
        p = tmp_path / "f.pdb"
        p.write_text("")
        with pytest.raises(ValueError, match="receptor_chain"):
            psp.ProteinStructurePrediction(
                path=p, receptor_chain="", effector_chain="B")

    def test_same_chain_raises(self, tmp_path):
        p = tmp_path / "f.pdb"
        p.write_text("")
        with pytest.raises(ValueError, match="cannot both be"):
            psp.ProteinStructurePrediction(
                path=p, receptor_chain="A", effector_chain="A")

    def test_invalid_predictor_raises(self, tmp_path):
        p = tmp_path / "f.pdb"
        p.write_text("")
        with pytest.raises(ValueError, match="predictor"):
            psp.ProteinStructurePrediction(
                path=p, receptor_chain="A", effector_chain="B",
                predictor="rosettafold")

    def test_basic_construction(self, tmp_path):
        p = tmp_path / "f.pdb"
        p.write_text("")
        x = psp.ProteinStructurePrediction(
            path=p, receptor_chain="A", effector_chain="B")
        assert x.receptor_chain == "A"
        assert x.effector_chain == "B"
        assert x.path == p
        assert x.predictor is None


@pytest.mark.logic
class TestThreeToOne:
    def test_standard_aas(self):
        assert psp.THREE_TO_ONE["ALA"] == "A"
        assert psp.THREE_TO_ONE["VAL"] == "V"
        assert psp.THREE_TO_ONE["TRP"] == "W"

    def test_non_standard(self):
        assert psp.THREE_TO_ONE["MSE"] == "M"


@pytest.mark.logic
class TestKabschRmsd:
    def test_identical_returns_zero(self):
        P = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
        assert psp._kabsch_rmsd(P, P.copy()) == pytest.approx(0.0)

    def test_translated_returns_zero(self):
        P = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
        Q = P + np.array([10.0, 20.0, 30.0])
        assert psp._kabsch_rmsd(P, Q) == pytest.approx(0.0, abs=1e-10)

    def test_rotated_returns_zero(self):
        P = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]], dtype=float)
        R = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]], dtype=float)
        Q = P @ R.T
        assert psp._kabsch_rmsd(P, Q) == pytest.approx(0.0, abs=1e-10)

    def test_fewer_than_3_points_returns_none(self):
        P = np.array([[0, 0, 0]], dtype=float)
        assert psp._kabsch_rmsd(P, P.copy()) is None

    def test_pairs_to_shorter(self):
        P = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [9, 9, 9]], dtype=float)
        Q = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], dtype=float)
        assert psp._kabsch_rmsd(P, Q) == pytest.approx(0.0, abs=1e-10)


@pytest.mark.logic
class TestRaEffRmsd:
    def test_perfect_overlap(self):
        rec = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], dtype=float)
        eff = np.array([[5, 5, 5], [6, 5, 5], [5, 6, 5]], dtype=float)
        assert psp._ra_eff_rmsd(rec, eff, rec, eff) == pytest.approx(0.0)

    def test_translated_complex_still_zero(self):
        rec = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], dtype=float)
        eff = np.array([[5, 5, 5], [6, 5, 5], [5, 6, 5]], dtype=float)
        t = np.array([10.0, 20.0, 30.0])
        assert psp._ra_eff_rmsd(rec + t, eff + t, rec, eff) == pytest.approx(
            0.0, abs=1e-10)

    def test_effector_displaced_independently(self):
        rec = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], dtype=float)
        eff = np.array([[5, 5, 5], [6, 5, 5], [5, 6, 5]], dtype=float)
        eff_d = eff + np.array([1.0, 0.0, 0.0])
        # With perfect receptor alignment (rec == rec), effector shift
        # of (1,0,0) means every atom is 1Å off → RMSD = 1
        assert psp._ra_eff_rmsd(rec, eff_d, rec, eff) == pytest.approx(1.0)

    def test_too_few_residues_returns_none(self):
        rec = np.array([[0, 0, 0]], dtype=float)
        eff = np.array([[0, 0, 0]], dtype=float)
        assert psp._ra_eff_rmsd(rec, eff, rec, eff) is None


def _has_gemmi():
    try:
        import gemmi  # noqa: F401
        return True
    except ImportError:
        return False


_REAL_PDB_DIR = (
    _REPO_ROOT
    / "tests/negative_steering/example_output_files"
    / "negative_steering/runs/design_42_seq_0/cycle_0/initial"
    / "boltz_results_input/predictions/input"
)


def _find_real_pdb():
    if not _REAL_PDB_DIR.is_dir():
        return None
    candidates = list(_REAL_PDB_DIR.glob("*.pdb"))
    return candidates[0] if candidates else None


@pytest.mark.logic
@pytest.mark.skipif(not _has_gemmi(), reason="gemmi not installed")
class TestRealStructure:
    def test_load_and_read_sequences(self):
        pdb = _find_real_pdb()
        if pdb is None:
            pytest.skip("no real PDB fixture available")
        x = psp.ProteinStructurePrediction(
            path=pdb, receptor_chain="A", effector_chain="B",
            predictor="boltz")
        rec = x.receptor_sequence()
        eff = x.effector_sequence()
        assert len(rec) > 10 and len(eff) > 10
        assert all(c in "ACDEFGHIKLMNPQRSTVWYUOX" for c in rec)
        assert all(c in "ACDEFGHIKLMNPQRSTVWYUOX" for c in eff)

    def test_lazy_loading(self):
        pdb = _find_real_pdb()
        if pdb is None:
            pytest.skip("no real PDB fixture available")
        x = psp.ProteinStructurePrediction(
            path=pdb, receptor_chain="A", effector_chain="B")
        assert x._structure_cache is None
        _ = x.receptor_sequence()
        assert x._structure_cache is not None

    def test_chain_ca_coords_shape(self):
        pdb = _find_real_pdb()
        if pdb is None:
            pytest.skip("no real PDB fixture available")
        x = psp.ProteinStructurePrediction(
            path=pdb, receptor_chain="A", effector_chain="B")
        ca = x.receptor_ca_coords()
        assert ca.ndim == 2 and ca.shape[1] == 3
        assert ca.shape[0] == x.receptor_length()

    def test_contact_residues_returns_position_set(self):
        pdb = _find_real_pdb()
        if pdb is None:
            pytest.skip("no real PDB fixture available")
        x = psp.ProteinStructurePrediction(
            path=pdb, receptor_chain="A", effector_chain="B")
        cs = x.receptor_contact_residues(cutoff=4.5)
        assert cs.chain == "A"
        assert cs.frame == "prediction"
        assert len(cs) > 0

    def test_ra_eff_against_self_is_zero(self):
        pdb = _find_real_pdb()
        if pdb is None:
            pytest.skip("no real PDB fixture available")
        x = psp.ProteinStructurePrediction(
            path=pdb, receptor_chain="A", effector_chain="B")
        x2 = psp.ProteinStructurePrediction(
            path=pdb, receptor_chain="A", effector_chain="B")
        assert x.ra_eff_vs(x2) == pytest.approx(0.0, abs=1e-6)

    def test_jaccard_against_self_is_one(self):
        pdb = _find_real_pdb()
        if pdb is None:
            pytest.skip("no real PDB fixture available")
        x = psp.ProteinStructurePrediction(
            path=pdb, receptor_chain="A", effector_chain="B")
        x2 = psp.ProteinStructurePrediction(
            path=pdb, receptor_chain="A", effector_chain="B")
        assert x.jaccard_interface_vs(x2, cutoff=4.5) == pytest.approx(1.0)

    def test_independent_rmsds_against_self_are_zero(self):
        pdb = _find_real_pdb()
        if pdb is None:
            pytest.skip("no real PDB fixture available")
        x = psp.ProteinStructurePrediction(
            path=pdb, receptor_chain="A", effector_chain="B")
        x2 = psp.ProteinStructurePrediction(
            path=pdb, receptor_chain="A", effector_chain="B")
        assert x.independent_receptor_rmsd_vs(x2) == pytest.approx(0.0, abs=1e-6)
        assert x.independent_effector_rmsd_vs(x2) == pytest.approx(0.0, abs=1e-6)


# ═════════════════════════════════════════════════════════════════════
# Reading a real structure built from known geometry
# ═════════════════════════════════════════════════════════════════════
#
# The block above covers construction and the two geometry helpers in
# isolation. Everything the class exposes beyond that reads the file, and the
# tests for it used to need a committed Boltz output that is not in this
# repository, so they skipped. The builders in fixtures_runtree.py place the
# two chains at a geometry stated at the call site, which makes a contact
# count, a chain split and a displacement arithmetic rather than measured.
#
# `complex_pdb` puts the receptor along x at z = 0 and the effector along x at
# z = `separation`, directly above the receptor residue beneath it. With
# spacing 3.8 the distance from receptor residue i to effector residue j is
# sqrt(separation^2 + (3.8 * (i - j))^2), so at separation 4.0 a residue is
# 4.0 from the one above it and 5.52 from its neighbour. Every cutoff below is
# chosen against those two numbers.

needs_gemmi = pytest.mark.skipif(not _has_gemmi(),
                                 reason="gemmi not installed")


def _seq_pdb(path: Path, residues) -> Path:
    """residues: (chain, resnum, three_letter_name), one CA each."""
    lines = [ca_line(i + 1, c, rn, i * 3.8, 0.0, 0.0, res=res)
             for i, (c, rn, res) in enumerate(residues)]
    lines.append("END")
    path.write_text("\n".join(lines) + "\n")
    return path


def _psp(path, receptor="A", effector="B", predictor=None):
    return psp.ProteinStructurePrediction(
        path=path, receptor_chain=receptor, effector_chain=effector,
        predictor=predictor)


# ── Identity ─────────────────────────────────────────────────────────

@pytest.mark.logic
def test_repr_names_the_file_and_both_chain_roles(tmp_path):
    """The repr reaches logs and assertion failures across the orthogonal
    stage. Without the chain roles in it a receptor/effector swap reads as
    two identical objects."""
    p = tmp_path / "design_4_seq_1.pdb"
    p.write_text("")
    r = repr(_psp(p))
    assert "design_4_seq_1.pdb" in r
    assert "receptor='A'" in r and "effector='B'" in r


@pytest.mark.logic
def test_an_empty_effector_chain_is_rejected(tmp_path):
    """The receptor-side check is pinned above. Both roles are required, and
    an empty effector would silently give every contact count zero."""
    p = tmp_path / "f.pdb"
    p.write_text("")
    with pytest.raises(ValueError, match="effector_chain"):
        _psp(p, receptor="A", effector="")


@pytest.mark.logic
@pytest.mark.parametrize("tag", ["boltz", "af3", None])
def test_the_predictor_tag_is_kept_as_given(tmp_path, tag):
    """Diagnostic only, per the comment at the check. Nothing branches on
    it, so it must round-trip rather than be normalised."""
    p = tmp_path / "f.pdb"
    p.write_text("")
    assert _psp(p, predictor=tag).predictor == tag


# ── Sequences ────────────────────────────────────────────────────────

@pytest.mark.logic
@needs_gemmi
def test_chain_sequence_maps_three_letter_names_including_the_modified_ones(
        tmp_path):
    """MSE is selenomethionine and reads as M, which is what a crystal
    structure of a receptor needs. An unrecognised name reads as X rather
    than being dropped, because dropping it would shift every position after
    it and silently renumber the design region."""
    p = _seq_pdb(tmp_path / "s.pdb", [
        ("A", 1, "ALA"), ("A", 2, "VAL"), ("A", 3, "TRP"),
        ("A", 4, "MSE"), ("A", 5, "UNK"), ("B", 1, "GLY"),
    ])
    x = _psp(p)
    assert x.receptor_sequence() == "AVWMX"
    assert x.effector_sequence() == "G"
    assert x.receptor_length() == 5
    assert x.effector_length() == 1


@pytest.mark.logic
@needs_gemmi
def test_an_absent_chain_reads_as_an_empty_sequence(tmp_path):
    """Documented: empty string rather than a raise. A prediction missing its
    effector has to produce a record the stage can report, not an exception
    that discards the batch."""
    p = _seq_pdb(tmp_path / "s.pdb", [("A", 1, "ALA"), ("A", 2, "GLY")])
    x = _psp(p)
    assert x.effector_sequence() == ""
    assert x.effector_length() == 0


@pytest.mark.logic
@needs_gemmi
def test_the_sequence_is_cached_after_the_first_read(tmp_path):
    """The stage reads the sequence, the coordinates and the contacts off one
    instance, and they must all describe the same structure. Rewriting the
    file underneath must not change what an already-read instance reports."""
    p = _seq_pdb(tmp_path / "s.pdb", [("A", 1, "ALA"), ("A", 2, "GLY")])
    x = _psp(p)
    assert x.receptor_sequence() == "AG"
    _seq_pdb(p, [("A", 1, "TRP"), ("A", 2, "TRP"), ("A", 3, "TRP")])
    assert x.receptor_sequence() == "AG"


# ── Coordinates ──────────────────────────────────────────────────────

@pytest.mark.logic
@needs_gemmi
def test_ca_coordinates_come_back_in_file_order(tmp_path):
    """Pairing against another prediction is positional, so file order is the
    frame every later comparison happens in. Residue 3 of this receptor sits
    at x = 2 * 3.8."""
    p = complex_pdb(tmp_path / "c.pdb", n_receptor=6, n_effector=3,
                    receptor_start=101)
    x = _psp(p)
    ca = x.receptor_ca_coords()
    assert ca.shape == (6, 3)
    assert ca[2][0] == pytest.approx(2 * 3.8)
    assert x.effector_ca_coords().shape == (3, 3)


@pytest.mark.logic
@needs_gemmi
def test_an_absent_chain_gives_an_empty_array_of_the_right_shape(tmp_path):
    """A `(0, 3)` array keeps the downstream `len()` and slicing working. A
    bare `(0,)` or a None would raise inside the RMSD helpers instead of
    returning None from them."""
    p = _seq_pdb(tmp_path / "s.pdb", [("A", 1, "ALA")])
    arr = _psp(p).effector_ca_coords()
    assert arr.shape == (0, 3)
    assert arr.dtype == float


@pytest.mark.logic
@needs_gemmi
def test_one_coordinate_per_residue_even_with_several_atoms(tmp_path):
    """The reader breaks after the first CA in a residue. Without that an
    alternate conformation would contribute twice and shift every later
    pairing by one."""
    p = tmp_path / "s.pdb"
    p.write_text(
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  CA  ALA A   1       1.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      3  CB  ALA A   1       2.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      4  CA  ALA A   2       4.800   0.000   0.000  1.00  0.00           C\n"
        "END\n")
    assert _psp(p).receptor_ca_coords().shape == (2, 3)


@pytest.mark.logic
@needs_gemmi
def test_hydrogens_are_left_out_of_the_heavy_atom_list(tmp_path):
    """Contacts are heavy-atom contacts. Counting hydrogens would pull the
    effective cutoff out by about an angstrom and widen every interface."""
    p = tmp_path / "s.pdb"
    p.write_text(
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  CA  ALA A   1       1.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      3  HA  ALA A   1       1.500   0.000   0.000  1.00  0.00           H\n"
        "END\n")
    heavy = _psp(p)._chain_heavy_atoms_list("A")
    assert [rn for rn, _xyz in heavy] == [1, 1]
    assert heavy[0][1][0] == pytest.approx(0.0)


# ── Contacts ─────────────────────────────────────────────────────────

@pytest.mark.logic
@needs_gemmi
def test_contact_residues_are_pdb_residue_numbers_in_the_prediction_frame(
        tmp_path):
    """The receptor is numbered from 101 and the effector sits 4.0 above its
    first four residues. At a 4.5 cutoff the contacts are 101 to 104, as PDB
    residue numbers rather than positional indices. The frame tag is what
    stops those being mixed with the 0-based protected set the steering
    engine reads."""
    p = complex_pdb(tmp_path / "c.pdb", n_receptor=10, n_effector=4,
                    receptor_start=101, separation=4.0)
    cs = _psp(p).receptor_contact_residues(cutoff=4.5)
    assert cs.as_sorted_list() == [101, 102, 103, 104]
    assert cs.chain == "A"
    assert cs.frame == "prediction"


@pytest.mark.logic
@needs_gemmi
@pytest.mark.parametrize("cutoff,expected", [
    (4.0, []),                              # exactly the separation
    (4.001, [101, 102, 103, 104]),          # a hair above it
    (5.5, [101, 102, 103, 104]),            # still inside the 5.517 neighbour
    (5.6, [101, 102, 103, 104, 105]),       # the neighbour now counts
])
def test_the_contact_cutoff_is_strictly_less_than(tmp_path, cutoff, expected):
    """`dist_sq < cutoff * cutoff`, so a pair at exactly the cutoff is not a
    contact. The effector sits 4.0 above receptor residues 101 to 104, and
    residue 105 is 5.517 from the nearest effector atom, so the four rows
    bracket both boundaries from each side."""
    p = complex_pdb(tmp_path / "c.pdb", n_receptor=10, n_effector=4,
                    receptor_start=101, separation=4.0)
    assert _psp(p).receptor_contact_residues(cutoff).as_sorted_list() == expected


@pytest.mark.logic
@needs_gemmi
def test_a_prediction_with_no_effector_has_an_empty_interface(tmp_path):
    """Still a PositionSet, carrying the chain and frame, so the caller's set
    algebra keeps working. Returning None would make every jaccard raise."""
    p = _seq_pdb(tmp_path / "s.pdb", [("A", 1, "ALA"), ("A", 2, "GLY")])
    cs = _psp(p).receptor_contact_residues(cutoff=4.5)
    assert cs.is_empty()
    assert cs.chain == "A" and cs.frame == "prediction"


# ── Pair comparisons ─────────────────────────────────────────────────

@pytest.mark.logic
@needs_gemmi
def test_ra_eff_measures_the_effector_displacement_after_aligning_receptors(
        tmp_path):
    """Both files hold the same receptor and an effector at a different
    height. Aligning the receptors leaves the effector 4.0 out, which is the
    number that decides whether a pose is correct."""
    a = complex_pdb(tmp_path / "a.pdb", n_receptor=8, n_effector=5,
                    separation=4.0)
    b = complex_pdb(tmp_path / "b.pdb", n_receptor=8, n_effector=5,
                    separation=8.0)
    assert _psp(a).ra_eff_vs(_psp(b)) == pytest.approx(4.0, abs=1e-5)


@pytest.mark.logic
@needs_gemmi
def test_the_independent_rmsds_discard_the_very_difference_ra_eff_reports(
        tmp_path):
    """This is why all three exist. The two effectors are the same shape in
    different places, so aligning them on their own gives 0.0 while the
    receptor-aligned measure gives 4.0. Reporting the independent number as
    a pose metric would call a displaced effector a perfect match."""
    a = complex_pdb(tmp_path / "a.pdb", n_receptor=8, n_effector=5,
                    separation=4.0)
    b = complex_pdb(tmp_path / "b.pdb", n_receptor=8, n_effector=5,
                    separation=8.0)
    x, y = _psp(a), _psp(b)
    assert x.independent_receptor_rmsd_vs(y) == pytest.approx(0.0, abs=1e-6)
    assert x.independent_effector_rmsd_vs(y) == pytest.approx(0.0, abs=1e-6)
    assert x.ra_eff_vs(y) == pytest.approx(4.0, abs=1e-5)


@pytest.mark.logic
@needs_gemmi
def test_a_chain_too_short_to_align_gives_none_rather_than_a_number(tmp_path):
    """Kabsch needs three points. None is the documented answer, and it must
    not degrade to a float computed from two."""
    a = complex_pdb(tmp_path / "a.pdb", n_receptor=2, n_effector=5)
    b = complex_pdb(tmp_path / "b.pdb", n_receptor=8, n_effector=5)
    assert _psp(a).ra_eff_vs(_psp(b)) is None
    assert _psp(a).independent_receptor_rmsd_vs(_psp(b)) is None


@pytest.mark.logic
@needs_gemmi
def test_jaccard_against_itself_is_one(tmp_path):
    p = complex_pdb(tmp_path / "c.pdb", n_receptor=10, n_effector=4,
                    receptor_start=101)
    assert _psp(p).jaccard_interface_vs(_psp(p), cutoff=4.5) == 1.0


@pytest.mark.logic
@needs_gemmi
def test_jaccard_is_the_overlap_over_the_union(tmp_path):
    """Contacts 101-105 against 103-107: three shared, seven between them.
    Reporting the intersection over either set alone would read 0.6 and make
    two interfaces that share three residues out of seven look close."""
    a = complex_pdb(tmp_path / "a.pdb", n_receptor=10, n_effector=5,
                    receptor_start=101, contact_residues=(0, 1, 2, 3, 4))
    b = complex_pdb(tmp_path / "b.pdb", n_receptor=10, n_effector=5,
                    receptor_start=101, contact_residues=(2, 3, 4, 5, 6))
    x, y = _psp(a), _psp(b)
    assert x.receptor_contact_residues(4.5).as_sorted_list() == [101, 102, 103,
                                                                 104, 105]
    assert y.receptor_contact_residues(4.5).as_sorted_list() == [103, 104, 105,
                                                                 106, 107]
    assert x.jaccard_interface_vs(y, cutoff=4.5) == pytest.approx(3 / 7)


@pytest.mark.logic
@needs_gemmi
def test_two_empty_interfaces_score_zero_not_one(tmp_path):
    """An empty set is identical to another empty set, so the set-theoretic
    answer would be 1.0 or undefined. The module returns 0.0, because two
    predictions with no interface at all agree about nothing worth scoring
    and a 1.0 here would read as a perfectly reproduced interface."""
    p = _seq_pdb(tmp_path / "s.pdb", [("A", 1, "ALA"), ("A", 2, "GLY")])
    x = _psp(p)
    assert x.jaccard_interface_vs(_psp(p), cutoff=4.5) == 0.0


# ── Loading failures ─────────────────────────────────────────────────

@pytest.mark.logic
@needs_gemmi
def test_an_unreadable_structure_raises_oserror_naming_the_path(tmp_path):
    """Construction only checks that the path exists, so a corrupt file is
    found on first access. The wrapped error has to name the file, because
    the stage reads hundreds of them."""
    bad = tmp_path / "broken.cif"
    bad.write_text("this is not a structure\n")
    x = _psp(bad)
    with pytest.raises(OSError, match="broken.cif"):
        x.receptor_sequence()


@pytest.mark.logic
def test_a_singular_alignment_returns_none_from_both_helpers(monkeypatch):
    """numpy raises LinAlgError when the SVD does not converge. Both helpers
    catch it and return None, so one pathological pair cannot take the stage
    down. Forced here because a converging-failure input cannot be written by
    hand."""
    def _raise(_):
        raise np.linalg.LinAlgError("forced")

    monkeypatch.setattr(np.linalg, "svd", _raise)
    pts = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
    assert psp._kabsch_rmsd(pts, pts.copy()) is None
    assert psp._ra_eff_rmsd(pts, pts, pts, pts) is None


@pytest.mark.logic
@needs_gemmi
def test_a_residue_number_seen_twice_contributes_one_letter(tmp_path):
    """De-duplication is on (residue number, insertion code), per the
    comment at the check. A structure that lists residue 1 twice under two
    names gives the first one, not both, because a second letter would shift
    every position after it and renumber the design region."""
    p = tmp_path / "dup.pdb"
    p.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      2  CA  GLY A   2       3.800   0.000   0.000  1.00  0.00           C\n"
        "ATOM      3  CA  TRP A   1       7.600   0.000   0.000  1.00  0.00           C\n"
        "END\n")
    assert _psp(p).receptor_sequence() == "AG"


@pytest.mark.logic
@needs_gemmi
def test_the_coordinates_are_cached_per_chain(tmp_path):
    """Same reason as the sequence cache. The second read must come from the
    cache rather than re-parsing, so a file rewritten underneath cannot make
    two reads off one instance disagree."""
    p = complex_pdb(tmp_path / "c.pdb", n_receptor=6, n_effector=3)
    x = _psp(p)
    first = x.receptor_ca_coords()
    assert x.receptor_ca_coords() is first
    complex_pdb(p, n_receptor=2, n_effector=2)
    assert x.receptor_ca_coords().shape == (6, 3)
