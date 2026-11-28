"""Cover the MPNN sequence correction stage.

This stage takes ProteinMPNN's output and rebuilds the receptor sequence so the
fixed segments carry the native residues and only the de novo positions carry
designed ones. Everything downstream, including the sequence that reaches
Boltz-2 and the one that is eventually ordered as a gBlock, is what this stage
writes.

Two of its functions have a recorded failure in their own docstrings, and both
are pinned below. `split_sequence` once assumed ProteinMPNN wrote the effector
first, and silently produced a "corrected receptor" built from effector
residues. `align_to_native_by_anchors` distributes de novo length across
segments, and getting that wrong shifts every residue after the first gap.
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

import pipeline_correct_sequences as pcs  # noqa: E402
from fixtures_runtree import complex_pdb  # noqa: E402

# ── Contig segments ──────────────────────────────────────────────────

@pytest.mark.logic
def test_segments_come_back_in_contig_order():
    """Order is what the alignment walks, so a reordering would put the
    designed residues in the wrong stretch."""
    segs = pcs.parse_contig_segments("A1-10/5-5/A21-30 B", "A")
    assert [s["type"] for s in segs] == ["fixed", "denovo", "fixed"]
    assert segs[0] == {"type": "fixed", "start": 1, "end": 10}
    assert segs[2] == {"type": "fixed", "start": 21, "end": 30}
    assert segs[1]["min_len"] == 5 and segs[1]["max_len"] == 5


@pytest.mark.logic
def test_a_variable_denovo_segment_keeps_both_bounds():
    segs = pcs.parse_contig_segments("A1-10/20-40/A21-30 B", "A")
    assert segs[1]["min_len"] == 20 and segs[1]["max_len"] == 40


@pytest.mark.logic
def test_receptor_chain_lookup_is_case_insensitive():
    assert pcs.parse_contig_segments("a1-10/5-5/a21-30 B", "A")


@pytest.mark.logic
@pytest.mark.parametrize("contigs,chain", [
    ("not a contig ???", "A"),
    ("A1-10/5-5/A21-30 B", "Z"),
])
def test_an_unusable_contig_gives_no_segments_rather_than_raising(contigs, chain):
    """An exception here would lose the whole run. An empty list makes the
    caller fall back to the length-based split, which it reports."""
    assert pcs.parse_contig_segments(contigs, chain) == []


# ── PDB residues ─────────────────────────────────────────────────────

@pytest.mark.logic
def test_chain_residues_are_per_chain_in_file_order(tmp_path):
    p = complex_pdb(tmp_path / "c.pdb", n_receptor=5, n_effector=3,
                    receptor_start=101, effector_start=501)
    chains = pcs.get_pdb_chain_residues(p)
    assert chains["A"] == [101, 102, 103, 104, 105]
    assert chains["B"] == [501, 502, 503]


@pytest.mark.logic
def test_a_residue_with_several_atoms_is_counted_once(tmp_path):
    """complex_pdb writes one CA per residue, so add a second atom by hand and
    check the reader de-duplicates on (chain, resnum). Counting atoms instead
    of residues would make every expected length wrong."""
    p = tmp_path / "c.pdb"
    p.write_text(
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  CA  ALA A   1       1.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      3  CA  ALA A   2       4.800   0.000   0.000  1.00  0.00           C\n"
        "END\n")
    assert pcs.get_pdb_chain_residues(p) == {"A": [1, 2]}


# ── FASTA ────────────────────────────────────────────────────────────

@pytest.mark.logic
def test_fasta_entries_keep_their_headers_and_join_wrapped_lines(tmp_path):
    f = tmp_path / "m.fa"
    f.write_text(">seq1, score=1.23\nAAAA\nCCCC\n>seq2, score=4.56\nDDDD\n")
    entries = pcs.parse_mpnn_fasta(f)
    assert len(entries) == 2
    assert entries[0][1] == "AAAACCCC"
    assert entries[1][0].startswith(">seq2")


@pytest.mark.logic
def test_an_empty_fasta_yields_no_entries(tmp_path):
    f = tmp_path / "empty.fa"
    f.write_text("")
    assert pcs.parse_mpnn_fasta(f) == []


@pytest.mark.logic
@pytest.mark.parametrize("header,expected", [
    (">d, score=1.2345, global_score=2.0", 1.2345),
    (">d, global_score=2.0, score=0.5", 0.5),
    (">d, no score here", 0.0),
    (">d, score=notanumber", 0.0),
])
def test_mpnn_score_is_read_or_defaults_to_zero(header, expected):
    """A header without a parsable score returns 0.0 rather than raising,
    because one malformed header must not discard a whole design set. 0.0
    sorts best for MPNN score, which is why selection must not rank on it
    alone."""
    assert pcs.extract_mpnn_score(header) == expected


# ── The chain-order bug ──────────────────────────────────────────────

@pytest.mark.logic
def test_split_assigns_chains_by_length_not_by_position():
    """ProteinMPNN's chain order depends on how the PDB was written. The
    previous implementation hardcoded effector-first and produced a
    "corrected receptor" built from effector residues, which then flowed
    downstream as if it were MPNN's receptor design.

    Receptor 20, effector 5. Whichever way round MPNN writes them, the
    5-residue part must come back as the effector.
    """
    rec, eff = "R" * 20, "E" * 5
    assert pcs.split_sequence(f"{rec}/{eff}", 5) == (eff, rec)
    assert pcs.split_sequence(f"{eff}/{rec}", 5) == (eff, rec)


@pytest.mark.logic
def test_both_separators_are_accepted():
    rec, eff = "R" * 20, "E" * 5
    assert pcs.split_sequence(f"{rec}:{eff}", 5) == (eff, rec)


@pytest.mark.logic
def test_without_a_separator_the_split_is_positional():
    """The documented fallback. It is not reached in this pipeline, because
    MPNN always emits a separator, and it is pinned so a change to it is
    deliberate."""
    assert pcs.split_sequence("EEEEERRRRR", 5) == ("EEEEE", "RRRRR")


# ── Anchor alignment ─────────────────────────────────────────────────

@pytest.mark.logic
def test_fixed_positions_take_native_and_denovo_takes_mpnn():
    """The point of the stage. Fixed stretches must come back as the native
    residues even though MPNN wrote something else there, and only the gap
    may carry designed residues.
    """
    segments = [
        {"type": "fixed", "start": 1, "end": 5},
        {"type": "denovo", "min_len": 3, "max_len": 3},
        {"type": "fixed", "start": 9, "end": 13},
    ]
    native = "NNNNN" + "xxx" + "CCCCC"          # 5 + 3 + 5
    mpnn = "mmmmm" + "DDD" + "mmmmm"            # MPNN rewrote the anchors too
    corrected, designed, native_at_design, observed = \
        pcs.align_to_native_by_anchors(mpnn, native, segments, 1)

    assert corrected.startswith("NNNNN"), "fixed head must be native"
    assert corrected.endswith("CCCCC"), "fixed tail must be native"
    assert designed == "DDD"
    assert "m" not in corrected, "no MPNN residue may survive at a fixed position"
    assert observed == "3"


@pytest.mark.logic
def test_the_designed_stretch_is_reported_at_its_observed_length():
    """A variable de novo segment resolves to whatever MPNN produced. The
    observed length is what the metrics table records, so it must be the real
    one rather than the contig's maximum."""
    # The residue numbering has to agree with the native string: fixed 1-4,
    # a five-residue native gap at 5-9, then fixed 10-13. Numbering the second
    # anchor 11-14 against a 13-residue native makes it read past the end and
    # come back a residue short, which is a broken fixture rather than a
    # finding.
    segments = [
        {"type": "fixed", "start": 1, "end": 4},
        {"type": "denovo", "min_len": 2, "max_len": 8},
        {"type": "fixed", "start": 10, "end": 13},
    ]
    native = "NNNN" + "xxxxx" + "CCCC"          # 4 + 5 + 4 = 13
    mpnn = "mmmm" + "DDDDDD" + "mmmm"           # six designed, gap was five
    corrected, designed, _nat, observed = \
        pcs.align_to_native_by_anchors(mpnn, native, segments, 1)
    assert designed == "DDDDDD"
    assert observed == "6"
    # The insertion is longer than the gap it replaces, and both anchors
    # survive intact. The corrected sequence grows rather than the tail being
    # truncated to hold the native length.
    assert corrected == "NNNN" + "DDDDDD" + "CCCC"


@pytest.mark.logic
def test_two_denovo_segments_are_reported_separately():
    """The lengths are '|'-joined per segment. Collapsing them to one number
    would hide which gap grew."""
    # Numbering agrees with the native string throughout: 1-3, gap 4-7,
    # 8-10, gap 11-14, 15-17.
    segments = [
        {"type": "fixed", "start": 1, "end": 3},
        {"type": "denovo", "min_len": 2, "max_len": 4},
        {"type": "fixed", "start": 8, "end": 10},
        {"type": "denovo", "min_len": 2, "max_len": 4},
        {"type": "fixed", "start": 15, "end": 17},
    ]
    native = "NNN" + "xxxx" + "MMM" + "xxxx" + "CCC"   # 3+4+3+4+3 = 17
    mpnn = "mmm" + "DDD" + "mmm" + "EEE" + "mmm"
    corrected, designed, _nat, observed = \
        pcs.align_to_native_by_anchors(mpnn, native, segments, 1)
    assert observed == "3|3"
    assert corrected == "NNN" + "DDD" + "MMM" + "EEE" + "CCC"
    # `designed` is per-segment and "|"-joined too, not one concatenated run.
    # Reading it as a single sequence would splice two unrelated stretches.
    assert designed == "DDD|EEE"


@pytest.mark.logic
def test_an_all_fixed_contig_designs_nothing():
    segments = [{"type": "fixed", "start": 1, "end": 6}]
    corrected, designed, _nat, _obs = \
        pcs.align_to_native_by_anchors("mmmmmm", "NNNNNN", segments, 1)
    assert corrected == "NNNNNN"
    assert designed == ""


@pytest.mark.logic
def test_an_mpnn_sequence_shorter_than_the_fixed_total_does_not_raise():
    """The module prints a warning and falls back rather than raising, because
    one short design must not take down the batch. Pinned so the fallback
    stays reachable."""
    segments = [
        {"type": "fixed", "start": 1, "end": 10},
        {"type": "denovo", "min_len": 3, "max_len": 3},
    ]
    corrected, designed, _nat, _obs = \
        pcs.align_to_native_by_anchors("short", "N" * 13, segments, 1)
    assert isinstance(corrected, str)
    assert isinstance(designed, str)
