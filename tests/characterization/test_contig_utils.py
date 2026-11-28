"""Cover the contig parser, which decides which residues survive redesign.

A contig string says which stretches of the receptor RFDiffusion keeps and
which it rebuilds. Every failure in it is quiet. A contig resolved against the
wrong residue numbers fixes the wrong stretch and the run completes normally,
producing designs that kept the segment they were meant to replace.

So the assertions below are on the residues themselves, not on the shape of
the output.

Format, which these tests also pin because it is easy to get wrong: segments
are separated by `/`, a fixed segment carries its chain letter (`A1-390`), and
a de novo segment is a bare length range (`20-40`). Commas are not separators.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
BIN = REPO / "bin"
if str(BIN) not in sys.path:
    sys.path.insert(0, str(BIN))

import contig_utils as cu  # noqa: E402


def _pdb(tmp_path: Path, chains: dict[str, tuple[int, int]]) -> Path:
    """A minimal PDB with one CA per residue, for the ranges given."""
    lines = []
    serial = 1
    for chain, (lo, hi) in chains.items():
        for rn in range(lo, hi + 1):
            lines.append(
                f"ATOM  {serial:5d}  CA  ALA {chain}{rn:4d}    "
                f"{0.0:8.3f}{0.0:8.3f}{float(rn):8.3f}  1.00  0.00           C")
            serial += 1
    lines.append("END")
    p = tmp_path / "x.pdb"
    p.write_text("\n".join(lines) + "\n")
    return p


# ── THREE_TO_ONE ─────────────────────────────────────────────────────

@pytest.mark.logic
def test_three_to_one_covers_the_twenty_standard_residues():
    """A gap here turns a real residue into X and shifts nothing, but a wrong
    letter silently changes a sequence that is compared downstream."""
    assert len(set(cu.THREE_TO_ONE.values()) & set("ACDEFGHIKLMNPQRSTVWY")) == 20
    assert cu.THREE_TO_ONE["GLY"] == "G"
    assert cu.THREE_TO_ONE["TRP"] == "W"


# ── Residue ranges off a PDB ─────────────────────────────────────────

@pytest.mark.logic
def test_chain_residue_range_reads_the_real_numbers(tmp_path):
    """Not 1..N. A structure numbered from 101 must report 101."""
    pdb = _pdb(tmp_path, {"A": (101, 180)})
    assert cu.get_chain_residue_range(pdb, "A") == (101, 180)


@pytest.mark.logic
def test_chain_residue_range_is_none_for_an_absent_chain(tmp_path):
    pdb = _pdb(tmp_path, {"A": (1, 10)})
    assert cu.get_chain_residue_range(pdb, "Z") == (None, None)


@pytest.mark.logic
def test_chain_residues_sorted_is_ordered_and_unique(tmp_path):
    pdb = _pdb(tmp_path, {"A": (5, 12)})
    got = cu.get_chain_residues_sorted(pdb, "A")
    assert got == sorted(set(got)) == list(range(5, 13))


@pytest.mark.logic
def test_an_unparseable_residue_number_is_skipped_rather_than_crashing(tmp_path):
    """Columns 23-26 of an ATOM record are not always an integer. A structure
    carrying `ABCD` there must not take the whole stage down, and the residues
    around it must still be read."""
    pdb = tmp_path / "bad.pdb"
    pdb.write_text(
        "ATOM      1  CA  ALA A   7       0.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      2  CA  ALA AABCD       3.800   0.000   0.000  1.00  0.00           C\n"
        "ATOM      3  CA  ALA A   9       7.600   0.000   0.000  1.00  0.00           C\n"
        "END\n")
    assert cu.get_chain_residue_range(pdb, "A") == (7, 9)
    assert cu.get_chain_residues_sorted(pdb, "A") == [7, 9]


# ── Block parsing ────────────────────────────────────────────────────

@pytest.mark.logic
def test_fixed_and_denovo_segments_parse_in_order():
    descs, chain = cu.parse_block_segments("A1-390/20-40/A421-438")
    assert chain == "A"
    assert descs == [
        ("fixed", "A", 1, 390),
        ("denovo", "20-40"),
        ("fixed", "A", 421, 438),
    ]


@pytest.mark.logic
def test_a_bare_chain_letter_is_a_passthrough():
    """`B` alone means 'keep this chain as it is', which is how the effector
    is carried through without being named residue by residue."""
    descs, chain = cu.parse_block_segments("B")
    assert chain == "B"
    assert descs == [("passthrough", "B")]


@pytest.mark.logic
def test_a_bare_chain_resolves_to_a_fixed_range_when_given_a_pdb(tmp_path):
    pdb = _pdb(tmp_path, {"B": (7, 96)})
    descs, _ = cu.parse_block_segments("B", pdb_path=pdb)
    assert descs == [("fixed", "B", 7, 96)]


@pytest.mark.logic
def test_a_chain_break_marker_becomes_its_own_descriptor():
    """`0` tells RFDiffusion the flanking residues are not contiguous. It
    carries no length, so it has to survive as a descriptor of its own rather
    than being dropped or counted as a residue."""
    descs, chain = cu.parse_block_segments("A1-10/0/A15-20")
    assert chain == "A"
    assert descs == [
        ("fixed", "A", 1, 10),
        ("break",),
        ("fixed", "A", 15, 20),
    ]


@pytest.mark.logic
def test_a_bare_chain_letter_inside_a_block_is_a_passthrough():
    """`A1-10/B/A15-20` keeps chain B whole in the middle of the receptor
    block. Without a PDB its range is unknown, so it stays a letter."""
    descs, chain = cu.parse_block_segments("A1-10/B/A15-20")
    assert chain == "A"
    assert descs[1] == ("passthrough", "B")


@pytest.mark.logic
def test_an_in_block_passthrough_resolves_against_the_pdb(tmp_path):
    pdb = _pdb(tmp_path, {"A": (1, 30), "B": (7, 96)})
    descs, _chain = cu.parse_block_segments("A1-10/B/A15-20", pdb_path=pdb)
    assert descs[1] == ("fixed", "B", 7, 96)


@pytest.mark.logic
def test_an_in_block_passthrough_for_an_absent_chain_stays_a_letter(tmp_path):
    """The PDB cannot supply a range, so the letter is carried through for
    RFDiffusion to expand. Inventing a range here would fix residues that are
    not in the structure."""
    pdb = _pdb(tmp_path, {"A": (1, 30)})
    descs, _chain = cu.parse_block_segments("A1-10/Z/A15-20", pdb_path=pdb)
    assert descs[1] == ("passthrough", "Z")


@pytest.mark.logic
def test_a_malformed_block_raises_naming_the_block():
    """RFDiffusion's own failure for a bad contig does not mention the contig,
    so the message here has to."""
    with pytest.raises(ValueError) as exc:
        cu.parse_block_segments("A1-2-3-4/??")
    assert "A1-2-3-4" in str(exc.value)


# ── Round trip ───────────────────────────────────────────────────────

@pytest.mark.logic
@pytest.mark.parametrize("block", [
    "A1-390/20-40/A421-438",
    "A1-100/15-15/A150-200",
    "A5-9",
    "A1-10/0/A15-20",       # the chain-break marker survives the round trip
    "A1-10/B/A15-20",       # so does an unresolved passthrough letter
])
def test_segments_to_string_round_trips(block):
    descs, _ = cu.parse_block_segments(block)
    assert cu.segments_to_string(descs) == block


# ── Expected chain lengths ───────────────────────────────────────────

@pytest.mark.logic
def test_expected_chain_lengths_returns_receptor_and_effector():
    """Returns a (rec_len, eff_len) pair, and a de novo segment contributes
    the MIDPOINT of its range, not its maximum.

    For 5-5 that is 5 either way. For 20-40 it is 30, so a caller that assumed
    the maximum would size the receptor 10 residues long and every positional
    check against it would be wrong at the tail.
    """
    assert cu.get_expected_chain_lengths("B A1-10/5-5/A21-30", "A", "B") == (25, 1)


@pytest.mark.logic
def test_denovo_length_uses_the_midpoint_of_the_range():
    rec, _eff = cu.get_expected_chain_lengths("B A1-10/20-40/A21-30", "A", "B")
    assert rec == 10 + (20 + 40) // 2 + 10


@pytest.mark.logic
def test_an_already_resolved_denovo_length_counts_once():
    """A bare `5` with no range is a length RFDiffusion has already sampled.
    Averaging it with itself, or reading it as residue 5, would both be
    wrong: the receptor is 10 + 5 + 10 residues long."""
    assert cu.get_expected_chain_lengths("B A1-10/5/A21-30", "A", "B") == (25, 1)


@pytest.mark.logic
def test_chain_break_marker_contributes_nothing():
    """`0` is RFDiffusion's chain-break marker, not a one-residue segment."""
    with_break = cu.get_expected_chain_lengths("B A1-10/0/A21-30", "A", "B")
    without = cu.get_expected_chain_lengths("B A1-10/A21-30", "A", "B")
    assert with_break == without


# ── Design region ────────────────────────────────────────────────────

@pytest.mark.logic
def test_design_region_is_the_native_gap_between_fixed_ranges():
    """Returns (design, fixed) as sets of NATIVE residue numbers.

    The de novo segment's own length does not decide the design region. The
    gap between the fixed ranges does. Here the contig asks for 3 new residues
    between A1-10 and A21-30, but the gap in native numbering is 11..20, so
    the design region is ten residues. Reading the segment length instead
    would protect ten residues while believing three were being rebuilt.
    """
    design, fixed = cu.parse_design_region("B A1-10/3-3/A21-30", "A")
    assert design == set(range(11, 21))
    assert fixed == set(range(1, 11)) | set(range(21, 31))
    assert not (design & fixed), "a residue cannot be both rebuilt and kept"


@pytest.mark.logic
def test_design_region_is_empty_with_a_single_fixed_range():
    """One fixed segment has no gap after it, so nothing is rebuilt."""
    design, fixed = cu.parse_design_region("B A1-30", "A")
    assert design == set()
    assert fixed == set(range(1, 31))


@pytest.mark.logic
def test_two_gaps_both_appear_in_the_design_region():
    design, fixed = cu.parse_design_region("B A1-5/2-2/A11-15/3-3/A21-25", "A")
    assert design == set(range(6, 11)) | set(range(16, 21))
    assert fixed == set(range(1, 6)) | set(range(11, 16)) | set(range(21, 26))


@pytest.mark.logic
def test_receptor_chain_lookup_is_case_insensitive():
    assert cu.parse_design_region("B a1-10/3-3/a21-30", "A")[0] == set(range(11, 21))


@pytest.mark.logic
def test_absent_receptor_chain_gives_empty_sets_rather_than_raising():
    """A caller that got an exception here would lose the whole run. Empty
    sets mean nothing is protected and nothing is rebuilt, which the design
    stage rejects on its own terms."""
    assert cu.parse_design_region("B A1-10", "Z") == (set(), set())


@pytest.mark.logic
def test_unparseable_contig_gives_empty_sets():
    assert cu.parse_design_region("not a contig at all ???", "A") == (set(), set())


@pytest.mark.logic
def test_a_chain_carried_through_whole_has_no_design_region():
    """The effector block is a bare letter, so it has no fixed segments to
    find gaps between. Nothing of it is rebuilt and nothing of it is pinned
    by residue number, which is not the same as it being absent."""
    assert cu.parse_design_region("A1-10/3-3/A21-30 B", "B") == (set(), set())


# ── De novo length resolution ────────────────────────────────────────

@pytest.mark.logic
def test_one_variable_segment_absorbs_whatever_is_left():
    """The PDB's residue count is the budget. A fixed total of 20 against a
    30-residue chain leaves 10 for the single variable segment, whatever its
    declared range said."""
    descs = [("fixed", "A", 1, 10), ("denovo", "5-20"), ("fixed", "A", 21, 30)]
    assert cu._resolve_denovo_lengths(descs, 30) == {1: 10}


@pytest.mark.logic
def test_a_fixed_length_denovo_segment_is_not_resolved():
    """5-5 is already decided, so it is subtracted from the budget rather
    than being handed the remainder."""
    descs = [("fixed", "A", 1, 10), ("denovo", "5-5"), ("fixed", "A", 16, 25)]
    assert cu._resolve_denovo_lengths(descs, 25) == {}


@pytest.mark.logic
def test_two_variable_segments_each_take_their_minimum():
    """With more than one unknown the split is underdetermined, so each takes
    its minimum rather than the remainder being guessed between them."""
    descs = [("fixed", "A", 1, 5), ("denovo", "3-9"),
             ("fixed", "A", 10, 14), ("denovo", "4-8")]
    assert cu._resolve_denovo_lengths(descs, 40) == {1: 3, 3: 4}


@pytest.mark.logic
def test_a_fixed_segment_is_subtracted_before_the_remainder():
    descs = [("fixed", "A", 1, 10), ("denovo", "2-2"), ("denovo", "1-50")]
    assert cu._resolve_denovo_lengths(descs, 30) == {2: 18}


# ── Remapping onto real PDB numbering ────────────────────────────────

@pytest.mark.logic
def test_segments_already_matching_the_pdb_are_left_alone(tmp_path):
    """No remap means no change. Rewriting numbers that already resolve would
    be a silent no-op at best and a shift at worst."""
    pdb = _pdb(tmp_path, {"A": (1, 30)})
    descs = [("fixed", "A", 1, 10), ("denovo", "5-5"), ("fixed", "A", 16, 30)]
    assert cu.remap_segments_to_pdb(descs, "A", pdb) == descs


@pytest.mark.logic
def test_segments_are_walked_onto_the_pdb_when_the_numbers_do_not_exist(tmp_path):
    """A contig written against 1..30 applied to a structure numbered
    101..130. The fixed stretches must land on real residues, keeping their
    lengths, or RFDiffusion fails with '(A, 1) is not in pdb file!' which says
    nothing about the contig that caused it.
    """
    pdb = _pdb(tmp_path, {"A": (101, 130)})
    descs = [("fixed", "A", 1, 10), ("denovo", "5-5"), ("fixed", "A", 16, 30)]
    out = cu.remap_segments_to_pdb(descs, "A", pdb)
    assert out[0] == ("fixed", "A", 101, 110)
    assert out[1] == ("denovo", "5-5")
    assert out[2][0] == "fixed"
    assert out[2][2] == 116, "the second anchor starts after the 5-residue gap"


@pytest.mark.logic
def test_remapping_preserves_every_fixed_length(tmp_path):
    pdb = _pdb(tmp_path, {"A": (501, 540)})
    descs = [("fixed", "A", 1, 12), ("denovo", "4-4"), ("fixed", "A", 17, 40)]
    out = cu.remap_segments_to_pdb(descs, "A", pdb)
    for before, after in zip(descs, out):
        if before[0] == "fixed":
            assert after[3] - after[2] == before[3] - before[2]


@pytest.mark.logic
def test_a_chain_break_marker_survives_remapping(tmp_path):
    """Remapping walks the PDB residue by residue. A break carries no length,
    so it must pass through untouched without consuming a residue: A1-10
    lands on 101-110 and the 15-residue anchor after the break lands on
    111-125, not on 112-126.
    """
    pdb = _pdb(tmp_path, {"A": (101, 130)})
    descs = [("fixed", "A", 1, 10), ("break",), ("fixed", "A", 16, 30)]
    assert cu.remap_segments_to_pdb(descs, "A", pdb) == [
        ("fixed", "A", 101, 110),
        ("break",),
        ("fixed", "A", 111, 125),
    ]


@pytest.mark.logic
def test_a_contig_naming_an_absent_chain_exits_with_the_reason(tmp_path):
    """RFDiffusion's own failure for this does not mention the contig, so the
    message here has to name the chain and the segments that referenced it."""
    pdb = _pdb(tmp_path, {"A": (1, 20)})
    descs = [("fixed", "Z", 1, 10)]
    with pytest.raises(SystemExit) as exc:
        cu.remap_segments_to_pdb(descs, "Z", pdb)
    assert "Z" in str(exc.value)


@pytest.mark.logic
def test_an_absent_chain_with_no_fixed_segments_passes_through(tmp_path):
    """Nothing references a residue number, so there is nothing to validate."""
    pdb = _pdb(tmp_path, {"A": (1, 20)})
    descs = [("denovo", "5-5")]
    assert cu.remap_segments_to_pdb(descs, "Z", pdb) == descs


@pytest.mark.logic
def test_no_block_chain_means_nothing_to_remap(tmp_path):
    pdb = _pdb(tmp_path, {"A": (1, 20)})
    descs = [("denovo", "5-5")]
    assert cu.remap_segments_to_pdb(descs, None, pdb) == descs


# ── resolve_contigs, end to end ──────────────────────────────────────

@pytest.mark.logic
def test_a_bare_chain_resolves_to_its_real_range(tmp_path):
    pdb = _pdb(tmp_path, {"A": (1, 30), "B": (7, 96)})
    assert cu.resolve_contigs("B", pdb) == "B7-96"


@pytest.mark.logic
def test_a_bare_chain_that_is_absent_is_passed_through_with_a_warning(tmp_path, capsys):
    """Passing it through keeps the failure where the user can see it, rather
    than inventing a range."""
    pdb = _pdb(tmp_path, {"A": (1, 30)})
    assert cu.resolve_contigs("Z", pdb) == "Z"
    assert "WARNING" in capsys.readouterr().err


@pytest.mark.logic
def test_a_full_contig_resolves_both_blocks(tmp_path):
    pdb = _pdb(tmp_path, {"A": (1, 30), "B": (1, 12)})
    out = cu.resolve_contigs("A1-10/5-5/A16-30 B", pdb)
    assert out.endswith("B1-12")
    assert out.split()[0].count("/") == 2
    assert "," not in out


@pytest.mark.logic
def test_commas_and_colons_are_accepted_as_block_separators(tmp_path):
    """Tolerated for back-compatibility. The canonical form is a space, and
    the output is canonical either way."""
    pdb = _pdb(tmp_path, {"A": (1, 30), "B": (1, 12)})
    space = cu.resolve_contigs("A1-30 B", pdb)
    assert cu.resolve_contigs("A1-30,B", pdb) == space
    assert cu.resolve_contigs("A1-30:B", pdb) == space
