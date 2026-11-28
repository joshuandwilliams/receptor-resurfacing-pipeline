"""Cover contig derivation from a posed complex.

This builds the RFdiffusion contig from the residues that contact the target.
Its output decides which stretches of the binder are rebuilt, so an error here
is an error in every design that follows, and the run completes normally.

The format is pinned deliberately, because it is a recurring source of
mistakes: segments are `/`-separated, a fixed segment carries its chain letter
as in `A1-30`, a de novo segment is a bare length range as in `5-11`, and the
target chain follows after a space. Commas are not separators.
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

import derive_pose_contigs as dpc  # noqa: E402
from fixtures_runtree import complex_pdb  # noqa: E402

# ── Reading ──────────────────────────────────────────────────────────

@pytest.mark.logic
def test_heavy_atoms_are_grouped_by_residue(tmp_path):
    p = complex_pdb(tmp_path / "c.pdb", n_receptor=6, n_effector=3,
                    receptor_start=1)
    heavy = dpc.read_heavy(p, "A")
    assert sorted(heavy) == [1, 2, 3, 4, 5, 6]
    assert heavy[1].shape[1] == 3


@pytest.mark.logic
def test_hydrogens_are_excluded(tmp_path):
    """The cutoffs below are heavy-atom distances. Including hydrogens would
    shrink every contact distance and widen the design region."""
    p = tmp_path / "h.pdb"
    p.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      2  HA  ALA A   1       0.500   0.000   0.000  1.00  0.00           H\n"
        "ATOM      3  CA  ALA A   2       3.800   0.000   0.000  1.00  0.00           C\n"
        "END\n")
    heavy = dpc.read_heavy(p, "A")
    assert len(heavy[1]) == 1, "the hydrogen must not be counted"


@pytest.mark.logic
def test_residue_numbers_come_back_in_file_order_without_duplicates(tmp_path):
    p = complex_pdb(tmp_path / "c.pdb", n_receptor=5, n_effector=2,
                    receptor_start=101)
    assert dpc.all_residue_numbers(p, "A") == [101, 102, 103, 104, 105]


# ── Contacts ─────────────────────────────────────────────────────────

@pytest.mark.logic
def test_contacts_are_the_residues_within_the_cutoff(tmp_path):
    """Effector residue j sits 4 A above receptor residue j, so at a 5 A
    cutoff exactly those directly-opposite residues contact."""
    p = complex_pdb(tmp_path / "c.pdb", n_receptor=10, n_effector=3,
                    spacing=3.8, separation=4.0, receptor_start=1)
    binder = dpc.read_heavy(p, "A")
    target = dpc.read_heavy(p, "B")
    assert dpc.contacts_at_cutoff(binder, target, 5.0) == {1, 2, 3}


@pytest.mark.logic
def test_a_separated_complex_has_no_contacts(tmp_path):
    p = complex_pdb(tmp_path / "far.pdb", separation=40.0)
    binder = dpc.read_heavy(p, "A")
    target = dpc.read_heavy(p, "B")
    assert dpc.contacts_at_cutoff(binder, target, 5.0) == set()


@pytest.mark.logic
def test_contacts_grow_monotonically_with_the_cutoff(tmp_path):
    p = complex_pdb(tmp_path / "c.pdb", n_receptor=12, n_effector=4,
                    separation=4.0)
    binder = dpc.read_heavy(p, "A")
    target = dpc.read_heavy(p, "B")
    sets = [dpc.contacts_at_cutoff(binder, target, c) for c in (4.5, 6.0, 10.0)]
    assert sets[0] <= sets[1] <= sets[2]


# ── Gap joining ──────────────────────────────────────────────────────

@pytest.mark.logic
def test_a_single_fixed_residue_between_two_design_residues_is_absorbed():
    """Leaving one fixed residue stranded inside a rebuilt stretch gives
    RFdiffusion a one-residue anchor it cannot usefully hold."""
    all_rn = list(range(1, 11))
    assert dpc.join_gap1({3, 5}, all_rn) == {3, 4, 5}


@pytest.mark.logic
def test_a_two_residue_gap_is_left_alone():
    """Only a gap of exactly one is closed. Closing two would start merging
    genuinely separate design regions."""
    all_rn = list(range(1, 11))
    assert dpc.join_gap1({3, 6}, all_rn) == {3, 6}


@pytest.mark.logic
def test_joining_is_idempotent():
    all_rn = list(range(1, 11))
    once = dpc.join_gap1({3, 5}, all_rn)
    assert dpc.join_gap1(once, all_rn) == once


# ── Segments ─────────────────────────────────────────────────────────

@pytest.mark.logic
def test_contiguous_residues_form_one_segment():
    assert dpc.segments({4, 5, 6}, list(range(1, 11))) == [(4, 6)]


@pytest.mark.logic
def test_separated_runs_form_separate_segments():
    assert dpc.segments({2, 3, 7, 8}, list(range(1, 11))) == [(2, 3), (7, 8)]


@pytest.mark.logic
def test_adjacency_follows_the_residue_list_not_arithmetic():
    """The binder may be numbered with gaps. Two residues whose numbers
    differ by more than one are still adjacent if nothing sits between them
    in the structure."""
    all_rn = [1, 2, 10, 11, 12]
    assert dpc.segments({2, 10}, all_rn) == [(2, 10)]


@pytest.mark.logic
def test_no_design_residues_gives_no_segments():
    assert dpc.segments(set(), list(range(1, 6))) == []


# ── The 4a length scaling ────────────────────────────────────────────

@pytest.mark.logic
def test_the_documented_example_reproduces():
    """The docstring states N=7 gives 5-11, because 7*1.5 = 10.5 rounds half
    up. Python's round() would give 10, so the rule is written as
    int(x + 0.5)."""
    assert dpc.strat4a_range(7) == (5, 11)


@pytest.mark.logic
def test_a_single_residue_stays_one_to_one():
    """Scaling 1 would ask RFdiffusion to rebuild a stretch that is not
    there."""
    assert dpc.strat4a_range(1) == (1, 1)


@pytest.mark.logic
@pytest.mark.parametrize("n", [2, 3, 5, 8, 13, 20])
def test_the_range_brackets_the_observed_length(n):
    lo, hi = dpc.strat4a_range(n)
    assert 1 <= lo <= n <= hi


@pytest.mark.logic
def test_the_range_never_inverts():
    for n in range(1, 60):
        lo, hi = dpc.strat4a_range(n)
        assert lo <= hi


# ── The contig itself ────────────────────────────────────────────────

@pytest.mark.logic
def test_a_contig_with_one_design_segment_in_the_middle():
    """Fixed head, de novo middle, fixed tail, then the target chain. The
    de novo part is a bare range with no chain letter, which is what makes it
    de novo rather than fixed."""
    all_rn = list(range(1, 31))
    contig, total = dpc.make_contig({10, 11, 12}, all_rn, "A", "B")
    assert contig == "A1-9/2-5/A13-30 B"
    assert total == 3


@pytest.mark.logic
def test_segments_are_slash_separated_and_the_target_follows_a_space():
    all_rn = list(range(1, 31))
    contig, _total = dpc.make_contig({10, 11, 12}, all_rn, "A", "B")
    assert "," not in contig, "commas are not separators in this grammar"
    head, target = contig.rsplit(" ", 1)
    assert target == "B"
    assert head.count("/") == 2


@pytest.mark.logic
def test_two_design_segments_keep_the_fixed_stretch_between_them():
    all_rn = list(range(1, 31))
    contig, total = dpc.make_contig({5, 6, 20, 21}, all_rn, "A", "B")
    assert contig.startswith("A1-4/")
    assert "/A7-19/" in contig, "the fixed stretch between the two must survive"
    assert contig.endswith("/A22-30 B")
    assert total == 4


@pytest.mark.logic
def test_a_design_segment_at_the_start_emits_no_leading_fixed_part():
    all_rn = list(range(1, 31))
    contig, _total = dpc.make_contig({1, 2, 3}, all_rn, "A", "B")
    assert not contig.startswith("A1-0")
    assert contig.split("/")[0] == "2-5"


@pytest.mark.logic
def test_a_design_segment_at_the_end_emits_no_trailing_fixed_part():
    all_rn = list(range(1, 31))
    contig, _total = dpc.make_contig({28, 29, 30}, all_rn, "A", "B")
    assert contig.endswith("2-5 B")


@pytest.mark.logic
def test_no_design_residues_gives_an_all_fixed_contig():
    """Nothing to rebuild means the binder passes through whole. An empty
    contig would make RFdiffusion rebuild everything."""
    all_rn = list(range(1, 31))
    contig, total = dpc.make_contig(set(), all_rn, "A", "B")
    assert contig == "A1-30 B"
    assert total == 0


@pytest.mark.logic
def test_the_total_counts_design_residues_after_gap_joining():
    """A residue absorbed by join_gap1 is rebuilt, so it counts."""
    all_rn = list(range(1, 31))
    _contig, total = dpc.make_contig({10, 12}, all_rn, "A", "B")
    assert total == 3


@pytest.mark.logic
def test_the_contig_parses_back_to_the_segments_it_describes():
    """The round trip through the repository's own parser, so this module and
    contig_utils cannot drift apart."""
    import contig_utils as cu
    all_rn = list(range(1, 31))
    contig, _total = dpc.make_contig({10, 11, 12}, all_rn, "A", "B")
    descs, chain = cu.parse_block_segments(contig.split(" ")[0])
    assert chain == "A"
    assert [d[0] for d in descs] == ["fixed", "denovo", "fixed"]
    assert descs[0] == ("fixed", "A", 1, 9)
    assert descs[2] == ("fixed", "A", 13, 30)
