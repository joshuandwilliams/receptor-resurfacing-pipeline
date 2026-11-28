"""Cover the Rosetta and biophysical metric runners, minus the binaries.

Both stages shell out to a tool that is not installed here, Rosetta's
InterfaceAnalyzer and FastRelax. What is testable without them is everything
either side of that call: deciding where the binary is, parsing what it wrote,
and turning a failure into a recorded row rather than an exception.

That parsing is where these fail quietly. A column read by the wrong index
returns a plausible number from a different metric, and every downstream table
still looks healthy. `notes/things_to_fix.md` already records that shape
complementarity never fires as a gate because it scores a bare backbone, so a
wrong Sc would not be noticed by eye either.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
BIN = REPO / "bin"
HERE = Path(__file__).resolve().parent
for _p in (str(BIN), str(HERE)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import run_biophysical_metrics as bio  # noqa: E402
import run_rosetta_metrics as rose  # noqa: E402
from fixtures_runtree import complex_pdb  # noqa: E402

# ── InterfaceAnalyzer score file ─────────────────────────────────────

def _scorefile(tmp_path: Path, header: str, row: str, seq: bool = True) -> Path:
    """A score.sc in the shape InterfaceAnalyzer writes it."""
    lines = []
    if seq:
        lines.append("SEQUENCE:")
    lines.append(f"SCORE: {header}")
    lines.append(f"SCORE: {row}")
    p = tmp_path / "score.sc"
    p.write_text("\n".join(lines) + "\n")
    return p


@pytest.mark.logic
def test_sc_and_ddg_are_read_from_their_named_columns(tmp_path):
    """Read by column NAME, not by position. InterfaceAnalyzer's column order
    changes with its options, so a positional read would return a different
    metric under a plausible name."""
    p = _scorefile(
        tmp_path,
        "total_score dG_separated sc_value description",
        "-123.4 -59.739 0.665 design_1")
    sc, ddg, err = rose._parse_ia_scorefile(p)
    assert err is None
    assert sc == pytest.approx(0.665)
    assert ddg == pytest.approx(-59.739)


@pytest.mark.logic
def test_column_order_does_not_change_the_result(tmp_path):
    """Same values, columns swapped. A positional read returns them the wrong
    way round and nothing downstream objects, because both are plausible."""
    p = _scorefile(
        tmp_path,
        "total_score sc_value dG_separated description",
        "-123.4 0.665 -59.739 design_1")
    sc, ddg, err = rose._parse_ia_scorefile(p)
    assert err is None
    assert sc == pytest.approx(0.665)
    assert ddg == pytest.approx(-59.739)


@pytest.mark.logic
def test_a_missing_column_is_reported_by_name(tmp_path):
    """The token says which column was absent, so a changed Rosetta option
    set is diagnosable from the row rather than from the log."""
    p = _scorefile(tmp_path, "total_score dG_separated description",
                   "-123.4 -59.739 design_1")
    sc, ddg, err = rose._parse_ia_scorefile(p)
    assert sc is None
    assert ddg == pytest.approx(-59.739)
    assert err is not None and "sc" in err


@pytest.mark.logic
def test_a_non_numeric_value_is_reported_rather_than_coerced(tmp_path):
    """Coercing to 0.0 would read as a real measurement of a terrible
    interface instead of as a parse failure."""
    p = _scorefile(tmp_path, "total_score dG_separated sc_value description",
                   "-123.4 nan_or_junk 0.665 design_1")
    sc, ddg, err = rose._parse_ia_scorefile(p)
    assert sc == pytest.approx(0.665)
    assert ddg is None
    assert err is not None


@pytest.mark.logic
def test_a_file_with_no_data_row_is_reported(tmp_path):
    p = tmp_path / "score.sc"
    p.write_text("SEQUENCE:\nSCORE: total_score sc_value dG_separated description\n")
    sc, ddg, err = rose._parse_ia_scorefile(p)
    assert (sc, ddg) == (None, None)
    assert err == "scorefile_no_data_row"


@pytest.mark.logic
def test_an_absent_file_is_reported_not_raised(tmp_path):
    """FastRelax can die without writing a score file. One design failing
    must not take the stage down."""
    sc, ddg, err = rose._parse_ia_scorefile(tmp_path / "nope.sc")
    assert (sc, ddg) == (None, None)
    assert err is not None and "read_failed" in err


@pytest.mark.logic
def test_the_first_data_row_is_the_one_taken(tmp_path):
    """One structure in, one row expected. Taking the last would silently
    report a different pose if Rosetta ever appended."""
    p = tmp_path / "score.sc"
    p.write_text(
        "SEQUENCE:\n"
        "SCORE: total_score dG_separated sc_value description\n"
        "SCORE: -1.0 -10.0 0.100 first\n"
        "SCORE: -2.0 -20.0 0.200 second\n")
    sc, ddg, err = rose._parse_ia_scorefile(p)
    assert err is None
    assert (sc, ddg) == (pytest.approx(0.100), pytest.approx(-10.0))


# ── Finding the Rosetta binary ───────────────────────────────────────

@pytest.mark.logic
def test_a_binary_that_is_not_installed_resolves_to_none():
    """Returning a path that does not exist would turn a missing install into
    a confusing exec failure three stages later."""
    assert rose._find_rosetta_binary("definitely_not_a_rosetta_binary_xyz") is None


@pytest.mark.logic
def test_a_binary_on_path_is_found(tmp_path, monkeypatch):
    fake = tmp_path / "InterfaceAnalyzer.linuxgccrelease"
    fake.write_text("#!/bin/sh\nexit 0\n")
    fake.chmod(0o755)
    monkeypatch.setenv("PATH", str(tmp_path))
    assert rose._find_rosetta_binary("InterfaceAnalyzer") is not None


# ── Interface residues ───────────────────────────────────────────────

@pytest.mark.logic
def test_interface_residues_are_pairs_of_chain_and_residue_number(tmp_path):
    """The set keys on (chain, resnum). Keying on resnum alone would merge a
    receptor and an effector residue that share a number."""
    pytest.importorskip("gemmi")
    p = complex_pdb(tmp_path / "c.pdb", n_receptor=10, n_effector=4,
                    separation=4.0, receptor_start=1, effector_start=1)
    residues, err = bio._interface_residues(p, "A", "B", 8.0)
    assert err is None
    assert residues
    assert all(isinstance(k, tuple) and len(k) == 2 for k in residues)
    assert {c for c, _n in residues} == {"A", "B"}


@pytest.mark.logic
def test_a_separated_complex_has_no_interface(tmp_path):
    pytest.importorskip("gemmi")
    p = complex_pdb(tmp_path / "far.pdb", separation=60.0)
    residues, err = bio._interface_residues(p, "A", "B", 8.0)
    assert err is None
    assert residues == set()


@pytest.mark.logic
def test_the_interface_grows_monotonically_with_the_cutoff(tmp_path):
    pytest.importorskip("gemmi")
    p = complex_pdb(tmp_path / "c.pdb", separation=4.0)
    sets = []
    for cutoff in (4.5, 6.0, 10.0):
        residues, err = bio._interface_residues(p, "A", "B", cutoff)
        assert err is None
        sets.append(residues)
    assert sets[0] <= sets[1] <= sets[2]


@pytest.mark.logic
def test_an_unreadable_structure_is_reported_not_raised(tmp_path):
    pytest.importorskip("gemmi")
    bad = tmp_path / "bad.pdb"
    bad.write_text("not a structure at all\n")
    residues, err = bio._interface_residues(bad, "A", "B", 8.0)
    assert residues is None or residues == set()
    if residues is None:
        assert isinstance(err, str) and err


@pytest.mark.logic
def test_an_absent_chain_is_named_rather_than_scored_as_no_interface(tmp_path):
    """The stronger of the two possible answers, and the one the code gives.

    An empty set would be indistinguishable from a real complex whose chains
    do not touch, so a typo in the chain id would read as a design that failed
    to bind. 'chain_not_found' cannot be mistaken for a measurement.
    """
    pytest.importorskip("gemmi")
    p = complex_pdb(tmp_path / "c.pdb", separation=4.0)
    residues, err = bio._interface_residues(p, "A", "Z", 8.0)
    assert residues is None
    assert err == "chain_not_found"
