"""Cover derivation of the design region and true interface from the input.

These two files are the contract between the pipeline and the steering engine.
The design region says which receptor positions were rebuilt. The true
interface says which must never be mutated. An index one out in either tells
the engine to mutate a residue it was told to keep, and the run completes.

The two use different frames on purpose, which is the thing most likely to be
collapsed by a refactor:

    design region    1-based positions on the receptor
    true interface   0-based positional indices

Both are positional within the receptor chain rather than PDB residue numbers,
so a structure that does not start at residue 1 still lines up.
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

import derive_input_design_region as didr  # noqa: E402
from fixtures_runtree import complex_pdb  # noqa: E402

# ── Contig parsing to the flat form ──────────────────────────────────

@pytest.mark.logic
def test_anchors_and_denovo_segments_come_back_distinguishable():
    """A de novo segment is marked with the sentinel chain id "_denovo" and
    carries its length bounds where an anchor carries residue numbers.
    Confusing the two would treat a length as a residue range."""
    segs = didr._parse_contigs("A1-10/5-15/A20-30 B", "A", "B")
    assert ("A", 1, 10, False) in segs
    assert ("A", 20, 30, False) in segs
    denovo = [s for s in segs if s[0] == "_denovo"]
    assert denovo == [("_denovo", 5, 15, True)]


@pytest.mark.logic
def test_a_bare_effector_chain_parses_to_a_zero_range():
    """`B` alone means carry the chain through. It has no residue range, so
    it is recorded as 0-0 rather than being dropped."""
    segs = didr._parse_contigs("A1-10/5-5/A16-25 B", "A", "B")
    assert ("B", 0, 0, False) in segs


@pytest.mark.logic
def test_commas_are_tolerated_but_are_not_the_canonical_form():
    """The canonical grammar is '/' within a chain and a space between
    chains. Commas are accepted for back-compatibility only, and this pins
    that they parse to the same thing rather than silently to something
    else."""
    canonical = didr._parse_contigs("A1-10/5-5/A16-25 B", "A", "B")
    with_commas = didr._parse_contigs("A1-10/5-5/A16-25,B", "A", "B")
    assert canonical == with_commas


# ── Receptor anchor walk ─────────────────────────────────────────────

@pytest.mark.logic
def test_only_receptor_anchors_are_walked():
    """The effector's own block must not contribute anchors, or the design
    region would be computed against the wrong chain."""
    segs = didr._parse_contigs("A1-10/5-5/A16-25 B", "A", "B")
    anchors = didr._receptor_positions_by_contigs(segs, "A")
    assert anchors == [(1, 10), (16, 25)]


@pytest.mark.logic
def test_an_all_fixed_receptor_gives_one_anchor():
    segs = didr._parse_contigs("A1-30 B", "A", "B")
    assert didr._receptor_positions_by_contigs(segs, "A") == [(1, 30)]


# ── Design region positions ──────────────────────────────────────────

@pytest.mark.logic
def test_the_documented_example_reproduces():
    """The docstring states: anchors [(1,50), (60,100)] on a 100-residue
    receptor give positions 51..59."""
    got = didr._design_region_positions_1b([(1, 50), (60, 100)], 100)
    assert got == list(range(51, 60))


@pytest.mark.logic
def test_positions_are_one_based():
    """Position 1 is the first residue. Returning 0 would shift the whole
    region by one against the file the engine reads."""
    got = didr._design_region_positions_1b([(4, 10)], 10)
    assert got == [1, 2, 3]
    assert min(got) == 1


@pytest.mark.logic
def test_a_gap_at_each_end_is_found():
    got = didr._design_region_positions_1b([(4, 7)], 10)
    assert got == [1, 2, 3, 8, 9, 10]


@pytest.mark.logic
def test_two_gaps_both_appear():
    got = didr._design_region_positions_1b([(1, 3), (7, 8)], 10)
    assert got == [4, 5, 6, 9, 10]


@pytest.mark.logic
def test_an_all_anchor_receptor_has_no_design_region():
    assert didr._design_region_positions_1b([(1, 10)], 10) == []


@pytest.mark.logic
def test_no_anchors_means_the_whole_receptor_is_rebuilt():
    """Documented edge case. An empty list here would mean nothing is
    rebuilt, which is the opposite of what no anchors means."""
    assert didr._design_region_positions_1b([], 5) == [1, 2, 3, 4, 5]


@pytest.mark.logic
def test_overlapping_anchors_do_not_double_count():
    got = didr._design_region_positions_1b([(1, 5), (3, 7)], 10)
    assert got == [8, 9, 10]


# ── True interface ───────────────────────────────────────────────────

@pytest.mark.logic
def test_the_interface_is_zero_based_positional(tmp_path):
    """Effector residue j sits 4 A above receptor residue j, so a 5 A cutoff
    gives the directly-opposite residues. The receptor starts at 101 here,
    and the answer must still be 0, 1, 2 rather than 101, 102, 103.
    """
    pdb = complex_pdb(tmp_path / "c.pdb", n_receptor=10, n_effector=3,
                      spacing=3.8, separation=4.0, receptor_start=101)
    seq = didr._extract_receptor_sequence(pdb, "A")
    idx = didr._compute_true_interface_0b(pdb, "A", "B", 5.0, seq)
    assert idx == [0, 1, 2]


@pytest.mark.logic
def test_the_interface_widens_with_the_cutoff(tmp_path):
    pdb = complex_pdb(tmp_path / "c.pdb", n_receptor=12, n_effector=4,
                      separation=4.0)
    seq = didr._extract_receptor_sequence(pdb, "A")
    tight = set(didr._compute_true_interface_0b(pdb, "A", "B", 4.5, seq))
    wide = set(didr._compute_true_interface_0b(pdb, "A", "B", 9.0, seq))
    assert tight <= wide


@pytest.mark.logic
def test_a_separated_complex_has_an_empty_interface(tmp_path):
    pdb = complex_pdb(tmp_path / "far.pdb", separation=50.0)
    seq = didr._extract_receptor_sequence(pdb, "A")
    assert didr._compute_true_interface_0b(pdb, "A", "B", 5.0, seq) == []


@pytest.mark.logic
def test_indices_are_sorted_and_unique(tmp_path):
    pdb = complex_pdb(tmp_path / "c.pdb", n_receptor=12, n_effector=5,
                      separation=4.0)
    seq = didr._extract_receptor_sequence(pdb, "A")
    idx = didr._compute_true_interface_0b(pdb, "A", "B", 8.0, seq)
    assert idx == sorted(set(idx))


# ── The two files, and the frame they disagree on ────────────────────

@pytest.mark.logic
def test_the_design_region_file_carries_its_provenance(tmp_path):
    out = tmp_path / "design_region.txt"
    pdb = complex_pdb(tmp_path / "c.pdb", n_receptor=10, n_effector=2)
    didr._write_design_region(out, [5, 6, 7], "A1-4/3-3/A8-10 B", pdb, "A")
    text = out.read_text()
    body = [ln for ln in text.splitlines() if ln and not ln.startswith("#")]
    assert body, "the positions themselves must be written"
    nums = [int(x) for ln in body for x in ln.replace(",", " ").split()]
    assert nums == [5, 6, 7]


@pytest.mark.logic
def test_the_true_interface_file_carries_its_provenance(tmp_path):
    out = tmp_path / "true_interface.txt"
    pdb = complex_pdb(tmp_path / "c.pdb", n_receptor=10, n_effector=2)
    didr._write_true_interface(out, [0, 2, 4], pdb, "A", "B", 5.0)
    body = [ln for ln in out.read_text().splitlines()
            if ln and not ln.startswith("#")]
    nums = [int(x) for ln in body for x in ln.replace(",", " ").split()]
    assert nums == [0, 2, 4]


@pytest.mark.logic
def test_the_two_frames_differ_by_one_on_the_same_residue(tmp_path):
    """The seam, stated as a test. For a receptor where residue 1 is both an
    interface residue and outside every anchor, the true interface records 0
    and the design region records 1.
    """
    pdb = complex_pdb(tmp_path / "c.pdb", n_receptor=10, n_effector=1,
                      spacing=3.8, separation=4.0, receptor_start=1)
    seq = didr._extract_receptor_sequence(pdb, "A")
    interface = didr._compute_true_interface_0b(pdb, "A", "B", 5.0, seq)
    assert interface == [0]
    design = didr._design_region_positions_1b([(2, 10)], 10)
    assert design == [1]
