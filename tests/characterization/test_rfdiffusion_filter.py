"""Drive rfdiffusion_filter against built structures with known answers.

This stage decides which RFDiffusion backbones go forward. It reads a PDB that
may be one chain or two, finds the chain break, works out which segment is the
receptor, counts the interface, and applies the filters.

Every fixture here is constructed, so the expected contact count, chain split
and displacement are arithmetic rather than something read off a previous run.
The geometry is two parallel lines of CA atoms a stated distance apart, which
makes a contact at cutoff c exactly the pairs whose separation is under c.
"""
from __future__ import annotations

import importlib.util
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

import rfdiffusion_filter as rf  # noqa: E402
from fixtures_runtree import complex_pdb  # noqa: E402

# ── Reading ──────────────────────────────────────────────────────────

@pytest.mark.logic
def test_read_ca_atoms_reads_one_atom_per_residue(tmp_path):
    p = complex_pdb(tmp_path / "c.pdb", n_receptor=20, n_effector=10)
    assert len(rf.read_ca_atoms(p)) == 30


@pytest.mark.logic
def test_read_ca_atoms_filters_by_chain(tmp_path):
    p = complex_pdb(tmp_path / "c.pdb", n_receptor=20, n_effector=10)
    assert len(rf.read_ca_atoms(p, chain="A")) == 20
    assert len(rf.read_ca_atoms(p, chain="B")) == 10


@pytest.mark.logic
def test_read_ca_atoms_preserves_file_order_and_coordinates(tmp_path):
    """Order is the frame every positional index is counted in, so a reader
    that sorted would silently renumber the design region."""
    p = complex_pdb(tmp_path / "c.pdb", n_receptor=5, n_effector=1,
                    spacing=3.8, receptor_start=101)
    rec = rf.read_ca_atoms(p, chain="A")
    assert [a["resnum"] for a in rec] == [101, 102, 103, 104, 105]
    assert rec[2]["x"] == pytest.approx(2 * 3.8)
    assert rec[0]["z"] == 0.0


@pytest.mark.logic
def test_count_chains(tmp_path):
    two = complex_pdb(tmp_path / "two.pdb")
    one = complex_pdb(tmp_path / "one.pdb", single_chain=True)
    assert rf.count_chains(two) == ["A", "B"]
    assert rf.count_chains(one) == ["A"]


# ── Chain break ──────────────────────────────────────────────────────

@pytest.mark.logic
def test_break_is_found_where_the_geometry_jumps(tmp_path):
    """RFDiffusion writes both chains under one id with continuous numbering.
    The only thing distinguishing them is the distance at the join, so the
    detector works on geometry rather than on numbering.
    """
    p = complex_pdb(tmp_path / "one.pdb", n_receptor=20, n_effector=10,
                    single_chain=True, separation=30.0)
    atoms = rf.read_ca_atoms(p)
    seg1, seg2, break_dist = rf.split_at_chain_break(atoms)
    assert len(seg1) == 20
    assert len(seg2) == 10
    # The join is the only large gap: 30 A in z plus 19*3.8 in x.
    assert break_dist == pytest.approx(math.hypot(19 * 3.8, 30.0), abs=1e-6)


@pytest.mark.logic
def test_a_single_atom_has_no_break_to_find():
    """A one-residue chain has no pair to measure, so the detector reports
    index 0 and distance 0 rather than indexing past the end."""
    assert rf._find_break_index([{"x": 1.0, "y": 2.0, "z": 3.0}]) == (0, 0.0)
    assert rf._find_break_index([]) == (0, 0.0)


@pytest.mark.logic
def test_the_split_keeps_every_atom(tmp_path):
    p = complex_pdb(tmp_path / "one.pdb", n_receptor=14, n_effector=7,
                    single_chain=True, separation=30.0)
    atoms = rf.read_ca_atoms(p)
    seg1, seg2, break_dist = rf.split_at_chain_break(atoms)
    assert len(seg1) + len(seg2) == len(atoms)


# ── Which segment is the receptor ────────────────────────────────────

@pytest.mark.logic
def test_segments_are_identified_by_the_effector_length(tmp_path):
    """The effector block of a contig is entirely fixed, so its length is
    exact and is the one to match against. Matching on the receptor would
    fail whenever a de novo segment resolved to a different length."""
    p = complex_pdb(tmp_path / "one.pdb", n_receptor=20, n_effector=10,
                    single_chain=True, separation=30.0)
    seg1, seg2, _ = rf.split_at_chain_break(rf.read_ca_atoms(p))
    rec, eff, rec_is_seg1 = rf.identify_segments(seg1, seg2, 20, 10)
    assert len(rec) == 20 and len(eff) == 10 and rec_is_seg1 is True


@pytest.mark.logic
def test_segments_are_identified_when_they_arrive_reversed(tmp_path):
    """The effector can come first in the file. Assuming the receptor leads
    would swap the two and every later metric would describe the wrong chain.
    """
    p = complex_pdb(tmp_path / "one.pdb", n_receptor=10, n_effector=20,
                    single_chain=True, separation=30.0)
    seg1, seg2, _ = rf.split_at_chain_break(rf.read_ca_atoms(p))
    rec, eff, rec_is_seg1 = rf.identify_segments(seg1, seg2, 20, 10)
    assert len(rec) == 20 and len(eff) == 10 and rec_is_seg1 is False


@pytest.mark.logic
def test_a_de_novo_length_that_resolved_short_still_identifies(tmp_path):
    """The receptor is 18 where 20 was expected, because a variable de novo
    segment resolved short. The effector is still exact, so the fuzzy match
    must pick on that."""
    p = complex_pdb(tmp_path / "one.pdb", n_receptor=18, n_effector=10,
                    single_chain=True, separation=30.0)
    seg1, seg2, _ = rf.split_at_chain_break(rf.read_ca_atoms(p))
    rec, eff, _ = rf.identify_segments(seg1, seg2, 20, 10)
    assert len(rec) == 18 and len(eff) == 10


@pytest.mark.logic
def test_the_fuzzy_match_picks_the_leading_segment_as_the_effector(tmp_path, capsys):
    """Neither segment matches exactly: 11 and 19 against an expected 20 and
    10. The effector length is the exact one, so the segment nearer to 10 is
    the effector, which here is the one that came first. Falling back on
    'the receptor leads' would describe the wrong chain in every metric.
    """
    p = complex_pdb(tmp_path / "one.pdb", n_receptor=11, n_effector=19,
                    single_chain=True, separation=30.0)
    seg1, seg2, _ = rf.split_at_chain_break(rf.read_ca_atoms(p))
    assert (len(seg1), len(seg2)) == (11, 19)

    rec, eff, rec_is_seg1 = rf.identify_segments(seg1, seg2, 20, 10)
    assert len(rec) == 19 and len(eff) == 11 and rec_is_seg1 is False
    assert "Fuzzy match" in capsys.readouterr().out


# ── Contacts ─────────────────────────────────────────────────────────

@pytest.mark.logic
def test_contact_count_is_arithmetic_at_a_given_cutoff(tmp_path):
    """Effector residue j sits 4.0 A above receptor residue j. A receptor
    residue k places away is sqrt(4^2 + (k*3.8)^2) distant. At cutoff 8.0
    that admits k = 0 and k = 1 on each side, because sqrt(16 + 14.44) = 5.52
    is under 8 and sqrt(16 + 57.76) = 8.59 is not.
    """
    p = complex_pdb(tmp_path / "c.pdb", n_receptor=20, n_effector=10,
                    spacing=3.8, separation=4.0)
    rec = rf.read_ca_atoms(p, chain="A")
    eff = rf.read_ca_atoms(p, chain="B")

    assert math.isclose(math.hypot(4.0, 3.8), 5.517, abs_tol=1e-3)
    assert math.isclose(math.hypot(4.0, 2 * 3.8), 8.588, abs_tol=1e-3)

    contacts = rf.find_contacts(rec, eff, 8.0)
    expected = sum(
        1
        for j in range(10)
        for k in range(20)
        if math.hypot(4.0, abs(k - j) * 3.8) < 8.0
    )
    assert len(contacts) == expected


@pytest.mark.logic
def test_a_separated_complex_has_no_contacts(tmp_path):
    p = complex_pdb(tmp_path / "far.pdb", separation=40.0)
    rec = rf.read_ca_atoms(p, chain="A")
    eff = rf.read_ca_atoms(p, chain="B")
    assert rf.find_contacts(rec, eff, 8.0) == []


@pytest.mark.logic
def test_contacts_grow_monotonically_with_the_cutoff(tmp_path):
    """A wider cutoff can only add pairs. A frame-of-reference bug usually
    does not preserve that ordering."""
    p = complex_pdb(tmp_path / "c.pdb", separation=4.0)
    rec = rf.read_ca_atoms(p, chain="A")
    eff = rf.read_ca_atoms(p, chain="B")
    counts = [len(rf.find_contacts(rec, eff, c)) for c in (5.0, 8.0, 12.0, 20.0)]
    assert counts == sorted(counts)


@pytest.mark.logic
def test_contact_distances_are_all_under_the_cutoff(tmp_path):
    p = complex_pdb(tmp_path / "c.pdb", separation=4.0)
    rec = rf.read_ca_atoms(p, chain="A")
    eff = rf.read_ca_atoms(p, chain="B")
    for _r, _e, d in rf.find_contacts(rec, eff, 8.0):
        assert d < 8.0


@pytest.mark.logic
def test_contact_indices_address_the_lists_they_came_from(tmp_path):
    """The pair is (receptor index, effector index). Returning residue
    numbers, or swapping the two, would make every downstream lookup read the
    wrong atom."""
    p = complex_pdb(tmp_path / "c.pdb", n_receptor=20, n_effector=10,
                    separation=4.0, receptor_start=101, effector_start=501)
    rec = rf.read_ca_atoms(p, chain="A")
    eff = rf.read_ca_atoms(p, chain="B")
    for r, e, d in rf.find_contacts(rec, eff, 8.0):
        assert 0 <= r < len(rec)
        assert 0 <= e < len(eff)
        dx = rec[r]["x"] - eff[e]["x"]
        dz = rec[r]["z"] - eff[e]["z"]
        assert math.isclose(math.hypot(dx, dz), d, abs_tol=1e-6)


@pytest.mark.logic
def test_an_empty_chain_contacts_nothing(tmp_path):
    p = complex_pdb(tmp_path / "c.pdb")
    rec = rf.read_ca_atoms(p, chain="A")
    assert rf.find_contacts(rec, [], 8.0) == []
    assert rf.find_contacts([], rec, 8.0) == []


# ── Coordinate array ─────────────────────────────────────────────────

@pytest.mark.logic
def test_coords_array_shape_and_order(tmp_path):
    p = complex_pdb(tmp_path / "c.pdb", n_receptor=7, n_effector=1)
    rec = rf.read_ca_atoms(p, chain="A")
    arr = rf.ca_coords_array(rec)
    assert arr.shape == (7, 3)
    assert arr[3][0] == pytest.approx(3 * 3.8)


@pytest.mark.logic
def test_coords_array_of_nothing_is_empty_not_an_error():
    assert rf.ca_coords_array([]).shape == (0, 3)


# ── Writing the split ────────────────────────────────────────────────

@pytest.mark.logic
def test_write_split_pdb_produces_two_chains_from_one(tmp_path):
    """The stage's job: take RFDiffusion's single-chain output and emit a
    two-chain complex the rest of the pipeline can read."""
    src = complex_pdb(tmp_path / "one.pdb", n_receptor=20, n_effector=10,
                      single_chain=True, separation=30.0)
    out = tmp_path / "split.pdb"
    rf.write_split_pdb(src, out, 20, 10)
    assert out.exists()
    assert rf.count_chains(out) == ["A", "B"]
    assert len(rf.read_ca_atoms(out, chain="A")) == 20
    assert len(rf.read_ca_atoms(out, chain="B")) == 10


@pytest.mark.logic
def test_the_written_split_assigns_chain_a_to_the_receptor_whichever_came_first(tmp_path):
    """The file holds 11 residues then 19, and the effector is expected to be
    10, so the leading 11 are the effector. Chain A must still be the
    receptor, because that is what every later stage reads as the binder.
    """
    src = complex_pdb(tmp_path / "one.pdb", n_receptor=11, n_effector=19,
                      single_chain=True, separation=30.0)
    out = tmp_path / "split.pdb"
    rf.write_split_pdb(src, out, 19, 10)
    assert len(rf.read_ca_atoms(out, chain="A")) == 19
    assert len(rf.read_ca_atoms(out, chain="B")) == 11


@pytest.mark.logic
def test_the_split_renumbers_each_chain_from_one(tmp_path):
    """split/design_*.pdb is 1-based output space, whatever the input
    numbering was. The module docstring says the metrics JSON is in
    contig space and these files are not, and downstream joins depend on
    exactly that."""
    src = complex_pdb(tmp_path / "one.pdb", n_receptor=12, n_effector=6,
                      single_chain=True, separation=30.0, receptor_start=201)
    out = tmp_path / "split.pdb"
    rf.write_split_pdb(src, out, 12, 6)
    assert [a["resnum"] for a in rf.read_ca_atoms(out, chain="A")] == list(range(1, 13))
    assert [a["resnum"] for a in rf.read_ca_atoms(out, chain="B")] == list(range(1, 7))


@pytest.mark.logic
def test_the_split_moves_no_atom(tmp_path):
    """Only chain ids and numbering change. A coordinate that moved would
    make every metric computed afterwards describe a different structure."""
    src = complex_pdb(tmp_path / "one.pdb", n_receptor=12, n_effector=6,
                      single_chain=True, separation=30.0)
    out = tmp_path / "split.pdb"
    rf.write_split_pdb(src, out, 12, 6)
    before = [(a["x"], a["y"], a["z"]) for a in rf.read_ca_atoms(src)]
    after = [(a["x"], a["y"], a["z"]) for a in rf.read_ca_atoms(out)]
    assert before == after


# ── Importing the file on its own ────────────────────────────────────

@pytest.mark.logic
def test_the_module_adds_its_own_directory_to_the_path_when_loaded_by_file():
    """A Nextflow stage runs this script by absolute path with no PYTHONPATH,
    so `import structure_metrics` and `from contig_utils import ...` resolve
    only because the module puts bin/ on sys.path itself. Every test file
    here adds bin/ before importing, which hides that line, and losing it
    would break the stage while the suite stayed green.
    """
    saved = list(sys.path)
    try:
        sys.path[:] = [p for p in sys.path
                       if Path(p or ".").resolve() != BIN.resolve()]
        assert not any(Path(p or ".").resolve() == BIN.resolve() for p in sys.path)

        spec = importlib.util.spec_from_file_location(
            "rfdiffusion_filter_loaded_by_path", BIN / "rfdiffusion_filter.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        assert any(Path(p or ".").resolve() == BIN.resolve() for p in sys.path)
        assert module.read_ca_atoms is not rf.read_ca_atoms, "a second load"
    finally:
        sys.path[:] = saved
