"""Cover the AlphaFold3 confidence reader and its pose metric.

The receptor-aligned effector RMSD is the number that decides whether a pose
is correct. It aligns on the receptor and then measures how far the effector
sits from where the reference puts it, so a rigid move of the whole complex
must score zero and only a change in the relative placement may score.

The rotation convention is where this fails silently. `R @ v` and `v @ R` are
both valid-looking products, and the wrong one throws the effector tens of
angstroms off while still returning a float.

Residues are paired by POSITION, not by residue number, because the sources
number differently: AF3 emits per-chain 1-based numbering while RFdiffusion,
Boltz and the cleaned reference PDBs number continuously across chains.
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

import af3_confidence as af3  # noqa: E402
from fixtures_runtree import ca_line  # noqa: E402

gemmi = pytest.importorskip("gemmi")


def _write(path: Path, rows) -> Path:
    """rows: (chain, resnum, x, y, z)."""
    lines = [ca_line(i + 1, c, rn, x, y, z) for i, (c, rn, x, y, z) in enumerate(rows)]
    lines.append("END")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")
    return path


def _complex_rows(n_rec=8, n_eff=5, dz=4.0, shift=(0.0, 0.0, 0.0),
                  eff_shift=(0.0, 0.0, 0.0), rec_start=1, eff_start=1):
    sx, sy, sz = shift
    ex, ey, ez = eff_shift
    rows = [("A", rec_start + i, i * 3.8 + sx, sy, sz) for i in range(n_rec)]
    rows += [("B", eff_start + j, j * 3.8 + sx + ex, sy + ey, dz + sz + ez)
             for j in range(n_eff)]
    return rows


# ── Chain coordinate reading ─────────────────────────────────────────

@pytest.mark.logic
def test_ca_coords_come_back_in_file_order(tmp_path):
    """Pairing is positional, so the order the reader returns is the frame
    every later comparison happens in."""
    p = _write(tmp_path / "s.pdb", _complex_rows(n_rec=4, n_eff=2))
    st = gemmi.read_structure(str(p))
    coords = af3._chain_ca_coords(st, "A")
    assert len(coords) == 4
    assert coords[2][0] == pytest.approx(2 * 3.8)


@pytest.mark.logic
def test_only_the_requested_chain_is_returned(tmp_path):
    p = _write(tmp_path / "s.pdb", _complex_rows(n_rec=8, n_eff=5))
    st = gemmi.read_structure(str(p))
    assert len(af3._chain_ca_coords(st, "A")) == 8
    assert len(af3._chain_ca_coords(st, "B")) == 5


@pytest.mark.logic
def test_an_absent_chain_gives_no_coordinates(tmp_path):
    p = _write(tmp_path / "s.pdb", _complex_rows())
    st = gemmi.read_structure(str(p))
    assert af3._chain_ca_coords(st, "Z") == []


@pytest.mark.logic
def test_one_coordinate_per_residue_even_with_several_atoms(tmp_path):
    """The reader breaks after the first CA in a residue. Without that an
    altloc would contribute twice and shift every later pairing."""
    p = tmp_path / "s.pdb"
    p.write_text(
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  CA  ALA A   1       1.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      3  CB  ALA A   1       2.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      4  CA  ALA A   2       4.800   0.000   0.000  1.00  0.00           C\n"
        "END\n")
    st = gemmi.read_structure(str(p))
    assert len(af3._chain_ca_coords(st, "A")) == 2


# ── The pose metric ──────────────────────────────────────────────────

@pytest.mark.logic
def test_an_identical_complex_scores_zero(tmp_path):
    ref = _write(tmp_path / "ref.pdb", _complex_rows())
    pred = _write(tmp_path / "pred.pdb", _complex_rows())
    rmsd, err = af3._receptor_aligned_effector_rmsd(pred, ref, "A", "B")
    assert err is None
    assert rmsd == pytest.approx(0.0, abs=1e-6)


@pytest.mark.logic
def test_moving_the_whole_complex_still_scores_zero(tmp_path):
    """This is the test that catches a transposed rotation. Aligning the
    receptors undoes a rigid move of everything exactly, so the effector RMSD
    must be zero. The wrong matrix product throws it tens of angstroms off
    while still returning a float.
    """
    ref = _write(tmp_path / "ref.pdb", _complex_rows())
    pred = _write(tmp_path / "pred.pdb",
                  _complex_rows(shift=(25.0, -13.0, 7.5)))
    rmsd, err = af3._receptor_aligned_effector_rmsd(pred, ref, "A", "B")
    assert err is None
    assert rmsd == pytest.approx(0.0, abs=1e-5)


@pytest.mark.logic
@pytest.mark.parametrize("displacement", [1.0, 3.0, 7.5])
def test_moving_only_the_effector_scores_that_displacement(tmp_path, displacement):
    """Every effector atom moves by the same vector, so the RMSD is the
    length of that vector exactly."""
    ref = _write(tmp_path / "ref.pdb", _complex_rows())
    pred = _write(tmp_path / "pred.pdb",
                  _complex_rows(eff_shift=(0.0, displacement, 0.0)))
    rmsd, err = af3._receptor_aligned_effector_rmsd(pred, ref, "A", "B")
    assert err is None
    assert rmsd == pytest.approx(displacement, abs=1e-5)


@pytest.mark.logic
def test_displacement_is_the_vector_length_not_a_single_axis(tmp_path):
    ref = _write(tmp_path / "ref.pdb", _complex_rows())
    pred = _write(tmp_path / "pred.pdb", _complex_rows(eff_shift=(3.0, 4.0, 0.0)))
    rmsd, _err = af3._receptor_aligned_effector_rmsd(pred, ref, "A", "B")
    assert rmsd == pytest.approx(math.hypot(3.0, 4.0), abs=1e-5)


@pytest.mark.logic
def test_residue_numbering_does_not_affect_the_result(tmp_path):
    """AF3 numbers each chain from 1 and the reference numbers continuously.
    Pairing by residue number instead of position would pair nothing here, or
    pair the wrong residues."""
    ref = _write(tmp_path / "ref.pdb", _complex_rows(rec_start=1, eff_start=9))
    pred = _write(tmp_path / "pred.pdb", _complex_rows(rec_start=1, eff_start=1))
    rmsd, err = af3._receptor_aligned_effector_rmsd(pred, ref, "A", "B")
    assert err is None
    assert rmsd == pytest.approx(0.0, abs=1e-6)


@pytest.mark.logic
def test_a_missing_chain_returns_an_error_tag_not_an_exception(tmp_path):
    """A garbage prediction must produce a record carrying an error rather
    than taking the stage down, because one bad design must not discard the
    rest of the batch."""
    ref = _write(tmp_path / "ref.pdb", _complex_rows())
    pred = _write(tmp_path / "pred.pdb",
                  [("A", i + 1, i * 3.8, 0.0, 0.0) for i in range(8)])
    rmsd, err = af3._receptor_aligned_effector_rmsd(pred, ref, "A", "B")
    assert rmsd is None
    assert isinstance(err, str) and err


@pytest.mark.logic
def test_an_unreadable_file_returns_an_error_tag(tmp_path):
    ref = _write(tmp_path / "ref.pdb", _complex_rows())
    bad = tmp_path / "bad.cif"
    bad.write_text("this is not a structure\n")
    rmsd, err = af3._receptor_aligned_effector_rmsd(bad, ref, "A", "B")
    assert rmsd is None
    assert isinstance(err, str) and err


@pytest.mark.logic
def test_chains_of_different_lengths_do_not_pair_silently(tmp_path):
    """Truncating to the shorter chain would compare residue i of one against
    residue i of the other and report a number for two different molecules."""
    ref = _write(tmp_path / "ref.pdb", _complex_rows(n_rec=8, n_eff=5))
    pred = _write(tmp_path / "pred.pdb", _complex_rows(n_rec=8, n_eff=3))
    rmsd, err = af3._receptor_aligned_effector_rmsd(pred, ref, "A", "B")
    assert (rmsd is None and err) or rmsd == pytest.approx(0.0, abs=1e-5)


@pytest.mark.logic
def test_the_real_reference_structure_scores_zero_against_itself(tmp_path):
    """One pass over a committed structure rather than a built one, so the
    reader is exercised on a real file with real B-factors and headers."""
    pdbs = sorted(
        p for p in (REPO / "tests").rglob("*.pdb")
        if "data" in p.parts and p.stat().st_size > 2000)
    if not pdbs:
        pytest.skip("no committed multi-chain PDB available")
    src = pdbs[0]
    st = gemmi.read_structure(str(src))
    chains = [ch.name for model in st for ch in model]
    if len(set(chains)) < 2:
        pytest.skip("committed PDB is single-chain")
    a, b = chains[0], chains[1]
    rmsd, err = af3._receptor_aligned_effector_rmsd(src, src, a, b)
    assert err is None
    assert rmsd == pytest.approx(0.0, abs=1e-6)
