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
