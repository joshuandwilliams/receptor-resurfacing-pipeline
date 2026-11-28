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

import csv
import json
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
from fixtures_runtree import ca_line, complex_pdb  # noqa: E402

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


# ── Anchor alignment at its boundaries ───────────────────────────────
#
# Three cases the walk handles explicitly and the block above does not reach:
# a de novo segment with no fixed segment before it, one with no fixed segment
# either side, and one whose flanking anchors are adjacent in the native
# numbering so there is no native stretch for it to have replaced.

@pytest.mark.logic
def test_a_leading_denovo_segment_takes_the_native_head_as_its_counterpart():
    """No fixed segment precedes it, so the native stretch it replaced is
    read from position 0. Without that case the first region's native
    composition comes back empty and the plots under-report it."""
    segments = [
        {"type": "denovo", "min_len": 3, "max_len": 3},
        {"type": "fixed", "start": 4, "end": 8},
    ]
    native = "xxx" + "CCCCC"
    corrected, designed, native_at_design, observed = \
        pcs.align_to_native_by_anchors("DDD" + "mmmmm", native, segments, 1)
    assert corrected == "DDDCCCCC"
    assert designed == "DDD"
    assert native_at_design == "xxx"
    assert observed == "3"


@pytest.mark.logic
def test_a_receptor_that_is_entirely_denovo_takes_the_whole_native():
    """Pathological but representable, per the comment on the boundary
    cases. With no anchor on either side the native counterpart is the whole
    receptor."""
    segments = [{"type": "denovo", "min_len": 4, "max_len": 4}]
    corrected, designed, native_at_design, _obs = \
        pcs.align_to_native_by_anchors("DDDD", "NNNN", segments, 1)
    assert corrected == "DDDD"
    assert designed == "DDDD"
    assert native_at_design == "NNNN"


@pytest.mark.logic
def test_an_insertion_between_adjacent_anchors_has_no_native_counterpart():
    """Anchors 1-5 and 6-10 are adjacent in the native numbering, so the
    three designed residues are a pure insertion and replaced nothing. The
    native counterpart is the empty string rather than a stretch borrowed
    from one of the anchors."""
    segments = [
        {"type": "fixed", "start": 1, "end": 5},
        {"type": "denovo", "min_len": 3, "max_len": 3},
        {"type": "fixed", "start": 6, "end": 10},
    ]
    corrected, designed, native_at_design, _obs = \
        pcs.align_to_native_by_anchors("NNNNNDDDCCCCC", "NNNNNCCCCC",
                                       segments, 1)
    assert corrected == "NNNNNDDDCCCCC"
    assert designed == "DDD"
    assert native_at_design == ""


@pytest.mark.logic
def test_a_variable_segment_clamped_to_its_bounds_says_so(capsys):
    """The contig asks for 20-40 de novo residues and the MPNN receptor only
    has room for 3. The allocation is clamped to the minimum and the walk
    then reads past the end, so the warning is the only signal that the
    corrected sequence is suspect."""
    segments = [
        {"type": "fixed", "start": 1, "end": 5},
        {"type": "denovo", "min_len": 20, "max_len": 40},
        {"type": "fixed", "start": 9, "end": 13},
    ]
    pcs.align_to_native_by_anchors("m" * 13, "N" * 13, segments, 1)
    out = capsys.readouterr().out
    assert "clamping to [20, 40]" in out
    assert "suggests 3" in out
    assert "does not match contig accounting" in out


# ═════════════════════════════════════════════════════════════════════
# gen_fixed_positions, driven through main()
# ═════════════════════════════════════════════════════════════════════
#
# The JSONL this mode writes is what ProteinMPNN is told to hold constant.
# Everything in it is a PDB residue number, not a positional index, and a
# residue wrongly listed there is a residue MPNN will never redesign. The
# tests below assert the file's exact bytes on layouts where every number is
# worked out at the call site.


def _merged_pdb(path, *, effector_len, receptor_len, start=1):
    """The single-chain layout RFDiffusion writes: the effector residues
    first, then the receptor, all under one chain with continuous numbering.

    `generate_fixed_positions` reads the first `effector_len` residue numbers
    as the effector and everything after as the receptor, so this order is
    the assumption under test rather than an arbitrary fixture choice.
    """
    n = effector_len + receptor_len
    lines = [ca_line(i + 1, "A", start + i, i * 3.8, 0.0, 0.0)
             for i in range(n)]
    lines.append("END")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")
    return path


def _split_pdb(path, chains):
    """chains: ordered (chain_id, first_resnum, n_residues) tuples, written
    in that order so the reader's file-order output is known."""
    lines, serial = [], 1
    for chain_id, first, n in chains:
        for i in range(n):
            lines.append(ca_line(serial, chain_id, first + i,
                                 i * 3.8, 0.0, 0.0))
            serial += 1
    lines.append("END")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")
    return path


def _run(monkeypatch, **kw):
    argv = ["pipeline_correct_sequences.py"]
    for key, value in kw.items():
        argv += [f"--{key}", str(value)]
    monkeypatch.setattr(sys, "argv", argv)
    return pcs.main()


@pytest.mark.logic
def test_fixed_positions_for_a_merged_chain_are_the_exact_resnums(
        tmp_path, monkeypatch, capsys):
    """Chain A holds four effector residues numbered 1-4 and then thirteen
    receptor residues numbered 5-17. The contig fixes receptor positions 1-5
    and 9-13 with three de novo residues between them, which consumes the
    receptor as 5-9 fixed, 10-12 de novo, 13-17 fixed.

    So the whole effector plus those ten receptor residues are held, and
    exactly the three de novo residues are free. The file is one JSON line
    keyed on the PDB stem, and it is asserted whole because every number in
    it is a residue MPNN will or will not be allowed to touch.
    """
    pdb = _merged_pdb(tmp_path / "design_0.pdb", effector_len=4,
                      receptor_len=13)
    out = tmp_path / "out"
    out.mkdir()
    _run(monkeypatch, mode="gen_fixed_positions", pdb_file=pdb,
         output_dir=out, receptor_seq="N" * 13, effector_seq="EEEE",
         contigs="A1-5/3-3/A9-13 B", num_designs=1)

    written = out / "fixed_positions_design_0.jsonl"
    assert written.read_text() == (
        '{"design_0": {"A": [1, 2, 3, 4, 5, 6, 7, 8, 9, '
        '13, 14, 15, 16, 17]}}\n'
    )
    assert "design_0: 3 free / 14 fixed (single chain A, 17 total)" in \
        capsys.readouterr().out


@pytest.mark.logic
def test_two_denovo_segments_split_the_slack_by_their_maximum(
        tmp_path, monkeypatch):
    """Fifteen fixed receptor residues and twenty-two in the PDB leave seven
    de novo to place across two segments with maxima 4 and 5. Seven split by
    maximum rounds to 3 and 4, so the receptor consumes as 5-9 fixed, three
    free, 13-17 fixed, four free, 22-26 fixed.

    Splitting them evenly instead would hand the second segment three
    residues and shift the last anchor by one, which is exactly the class of
    error the fixed/variable split in the sibling function exists to avoid.
    """
    pdb = _merged_pdb(tmp_path / "design_0.pdb", effector_len=4,
                      receptor_len=22)
    out = tmp_path / "out"
    out.mkdir()
    _run(monkeypatch, mode="gen_fixed_positions", pdb_file=pdb,
         output_dir=out, receptor_seq="N" * 22, effector_seq="EEEE",
         contigs="A1-5/2-4/A9-13/3-5/A18-22 B", num_designs=1)

    assert json.loads((out / "fixed_positions_design_0.jsonl").read_text()) == {
        "design_0": {"A": [1, 2, 3, 4, 5, 6, 7, 8, 9,
                           13, 14, 15, 16, 17,
                           22, 23, 24, 25, 26]}
    }


@pytest.mark.logic
def test_fixed_positions_for_a_split_pdb_name_both_chains(
        tmp_path, monkeypatch, capsys):
    """Receptor chain A numbered 101-113 and effector chain B numbered 1-4.
    The effector is held whole and the receptor is walked the same way as the
    merged case, so 101-105 and 109-113 are held and 106-108 are free.

    The effector chain is listed first because it is written first, and the
    dict order is what reaches the JSONL.
    """
    pdb = _split_pdb(tmp_path / "design_7.pdb",
                     [("A", 101, 13), ("B", 1, 4)])
    out = tmp_path / "out"
    out.mkdir()
    _run(monkeypatch, mode="gen_fixed_positions", pdb_file=pdb,
         output_dir=out, receptor_seq="N" * 13, effector_seq="EEEE",
         contigs="A1-5/3-3/A9-13 B", num_designs=1)

    assert (out / "fixed_positions_design_7.jsonl").read_text() == (
        '{"design_7": {"B": [1, 2, 3, 4], '
        '"A": [101, 102, 103, 104, 105, 109, 110, 111, 112, 113]}}\n'
    )
    assert "design_7: 3 free / 10 fixed receptor / 4 fixed effector" in \
        capsys.readouterr().out


@pytest.mark.logic
def test_the_effector_chain_is_chosen_by_length_not_by_its_chain_id(
        tmp_path, monkeypatch):
    """Here chain A is the four-residue effector and chain B is the receptor,
    the opposite of the pipeline's usual split, and --receptor_chain still
    says A because that names the chain in the contig rather than in the PDB.

    The module picks the chain whose length is closest to the effector
    sequence, so it holds all of A and walks B. That is deliberate: the
    receptor length varies with the insertion and cannot be matched, so only
    the effector can be identified by length. It also means the chain ids in
    the output come from the PDB's geometry and not from --receptor_chain.
    """
    pdb = _split_pdb(tmp_path / "design_0.pdb",
                     [("A", 1, 4), ("B", 101, 13)])
    out = tmp_path / "out"
    out.mkdir()
    _run(monkeypatch, mode="gen_fixed_positions", pdb_file=pdb,
         output_dir=out, receptor_seq="N" * 13, effector_seq="EEEE",
         contigs="A1-5/3-3/A9-13 B", num_designs=1)

    assert json.loads((out / "fixed_positions_design_0.jsonl").read_text()) == {
        "design_0": {"A": [1, 2, 3, 4],
                     "B": [101, 102, 103, 104, 105,
                           109, 110, 111, 112, 113]}
    }


@pytest.mark.logic
def test_a_third_chain_is_left_out_of_the_file_entirely(tmp_path, monkeypatch):
    """Per the comment at the branch, more than two chains means best
    effector match plus largest remaining chain as receptor. Chain C is
    neither, so no residue of it is listed, which leaves every one of them
    free for MPNN to redesign. Pinned because the pipeline's own PDBs have
    two chains and a third would pass through silently.
    """
    pdb = _split_pdb(tmp_path / "design_0.pdb",
                     [("A", 101, 13), ("B", 1, 4), ("C", 301, 6)])
    out = tmp_path / "out"
    out.mkdir()
    _run(monkeypatch, mode="gen_fixed_positions", pdb_file=pdb,
         output_dir=out, receptor_seq="N" * 13, effector_seq="EEEE",
         contigs="A1-5/3-3/A9-13 B", num_designs=1)

    written = json.loads(
        (out / "fixed_positions_design_0.jsonl").read_text())
    assert sorted(written["design_0"]) == ["A", "B"]
    assert 301 not in written["design_0"]["A"]


@pytest.mark.logic
def test_the_design_loop_skips_the_designs_that_are_not_there(
        tmp_path, monkeypatch):
    """RFDiffusion's filter leaves gaps in the design numbering, so the loop
    tests for each file rather than assuming a dense range. A missing
    design_1 must not stop design_2 from being written."""
    rfdiff = tmp_path / "rfdiffusion"
    _merged_pdb(rfdiff / "design_0.pdb", effector_len=4, receptor_len=13)
    _merged_pdb(rfdiff / "design_2.pdb", effector_len=4, receptor_len=13)
    out = tmp_path / "out"
    out.mkdir()
    _run(monkeypatch, mode="gen_fixed_positions", rfdiff_dir=rfdiff,
         output_dir=out, receptor_seq="N" * 13, effector_seq="EEEE",
         contigs="A1-5/3-3/A9-13 B", num_designs=3)

    assert sorted(p.name for p in out.glob("*.jsonl")) == [
        "fixed_positions_design_0.jsonl",
        "fixed_positions_design_2.jsonl",
    ]


@pytest.mark.logic
def test_gen_fixed_positions_needs_somewhere_to_read_pdbs_from(
        tmp_path, monkeypatch):
    """Neither --rfdiff_dir nor --pdb_file would otherwise give an empty file
    list and a silent no-op, which in a Nextflow stage reads as success."""
    out = tmp_path / "out"
    out.mkdir()
    with pytest.raises(AssertionError, match="rfdiff_dir"):
        _run(monkeypatch, mode="gen_fixed_positions", output_dir=out,
             receptor_seq="N" * 13, effector_seq="EEEE",
             contigs="A1-5/3-3/A9-13 B", num_designs=1)


@pytest.mark.logic
def test_a_contig_naming_an_absent_chain_frees_the_whole_receptor(
        tmp_path, monkeypatch, capsys):
    """--receptor_chain Z matches no chain in the contig, so the segment list
    is empty and the walk has nothing to hold. The effector is still fixed
    and every receptor residue is left free.

    This is the current behaviour and it is pinned rather than fixed. It
    matters because the failure is silent: the stage exits 0, writes a
    well-formed JSONL, and ProteinMPNN then redesigns the receptor's fixed
    anchors as well as its de novo block. The sibling `correct` mode at
    least prints a warning on the same input; this mode prints the empty
    segment list and nothing else.
    """
    pdb = _merged_pdb(tmp_path / "design_0.pdb", effector_len=4,
                      receptor_len=13)
    out = tmp_path / "out"
    out.mkdir()
    _run(monkeypatch, mode="gen_fixed_positions", pdb_file=pdb,
         output_dir=out, receptor_seq="N" * 13, effector_seq="EEEE",
         contigs="A1-5/3-3/A9-13 B", receptor_chain="Z", num_designs=1)

    written = json.loads(
        (out / "fixed_positions_design_0.jsonl").read_text())
    assert written == {"design_0": {"A": [1, 2, 3, 4]}}
    assert "Contig segments: []" in capsys.readouterr().out


@pytest.mark.logic
@pytest.mark.xfail(strict=False, reason=(
    "Reported rather than fixed: an empty segment list is indistinguishable "
    "from a contig that fixes nothing, so gen_fixed_positions writes a JSONL "
    "freeing the entire receptor and exits 0.  A contig that resolves to no "
    "segments for the named receptor chain is always a specification error "
    "and should abort or at least warn.  Adding the check would change the "
    "stage's exit behaviour on inputs that currently 'succeed', so it is "
    "pinned here instead."))
def test_an_empty_segment_list_should_be_reported_as_an_error(
        tmp_path, monkeypatch, capsys):
    pdb = _merged_pdb(tmp_path / "design_0.pdb", effector_len=4,
                      receptor_len=13)
    out = tmp_path / "out"
    out.mkdir()
    _run(monkeypatch, mode="gen_fixed_positions", pdb_file=pdb,
         output_dir=out, receptor_seq="N" * 13, effector_seq="EEEE",
         contigs="A1-5/3-3/A9-13 B", receptor_chain="Z", num_designs=1)
    assert "WARNING" in capsys.readouterr().out


# ═════════════════════════════════════════════════════════════════════
# correct mode, driven through main()
# ═════════════════════════════════════════════════════════════════════
#
# This mode writes the sequences that reach Boltz-2 and, for a hit, the
# gBlock that gets ordered. Three artefacts come out of one run and they have
# to agree: a per-sequence AF2 FASTA, a combined FASTA, and the metadata CSV
# the plots read. All three are asserted whole below.
#
# The run used throughout: native receptor NNNNN xxx CCCCC numbered from 1,
# native effector EEEE, and the contig A1-5/3-3/A9-13 B. So anchors 1-5 and
# 9-13 come back native, the three residues between them are MPNN's, and the
# native stretch they replaced is "xxx".

_NATIVE_RECEPTOR = "NNNNN" + "xxx" + "CCCCC"
_CONTIG = "A1-5/3-3/A9-13 B"


def _mpnn_fasta(path, designed):
    """The FASTA ProteinMPNN writes: a native entry first, which the stage
    skips by slicing it off, then one entry per sampled sequence.

    `designed` is a list of (score, sequence) pairs, each sequence already
    carrying MPNN's '/' chain separator.
    """
    lines = [">design_0, score=0.0000, global_score=0.0000, seq_recovery=1.0",
             _NATIVE_RECEPTOR + "/EEEE"]
    for i, (score, seq) in enumerate(designed):
        lines.append(f">T=0.1, sample={i + 1}, score={score}, "
                     f"global_score=2.0000, seq_recovery=0.40")
        lines.append(seq)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")
    return path


def _mpnn_tree(tmp_path, designs):
    """designs: {design index: [(score, sequence), ...]}."""
    root = tmp_path / "mpnn"
    for m, entries in designs.items():
        _mpnn_fasta(root / f"design_{m}" / "seqs" / f"design_{m}.fa", entries)
    return root


def _correct(monkeypatch, tmp_path, mpnn_dir, **extra):
    out = tmp_path / "out"
    kw = dict(mpnn_dir=mpnn_dir, output_dir=out,
              receptor_seq=_NATIVE_RECEPTOR, effector_seq="EEEE",
              contigs=_CONTIG, num_designs=1)
    kw.update(extra)
    _run(monkeypatch, **kw)
    return out


@pytest.mark.logic
def test_the_af2_fasta_holds_the_corrected_receptor_and_the_native_effector(
        tmp_path, monkeypatch):
    """MPNN rewrote the anchors and the effector. The anchors come back
    native and the effector is replaced wholesale by the native one, because
    the effector is the target and is never a design target. Keeping MPNN's
    effector here would hand Boltz a mutated target and make every pose
    metric meaningless.
    """
    mpnn = _mpnn_tree(tmp_path, {0: [("1.234", "mmmmmDDDmmmmm/eeee")]})
    out = _correct(monkeypatch, tmp_path, mpnn)
    assert (out / "af2_fastas" / "design_0_seq_0.fasta").read_text() == (
        ">receptor\nNNNNNDDDCCCCC\n>effector\nEEEE\n"
    )


@pytest.mark.logic
def test_the_combined_fasta_carries_the_per_sequence_provenance(
        tmp_path, monkeypatch):
    """One header line per sequence, carrying the MPNN score, the designed
    residues, the length this design actually produced, the length the contig
    asked for, and the change count. The spec and the observed length are
    separate fields on purpose, so a design that came out at a different
    length than requested is visible rather than inferred."""
    mpnn = _mpnn_tree(tmp_path, {0: [("1.234", "mmmmmDDDmmmmm/eeee")]})
    out = _correct(monkeypatch, tmp_path, mpnn)
    assert (out / "mpnn_corrected.fasta").read_text() == (
        ">design_0_seq_0|mpnn=1.234|designed=DDD"
        "|design_region_length_observed=3|design_region_spec=3-3|changes=3\n"
        "NNNNNDDDCCCCC:EEEE\n"
    )


@pytest.mark.logic
def test_the_metadata_csv_is_one_row_per_sequence_with_its_columns_in_order(
        tmp_path, monkeypatch):
    """The plots and the cross-sequence summary read this file by column
    name, so both the header and the values are pinned. `num_changes` is 3
    because all three designed residues differ from the native 'xxx' they
    replaced."""
    mpnn = _mpnn_tree(tmp_path, {0: [("1.234", "mmmmmDDDmmmmm/eeee")]})
    out = _correct(monkeypatch, tmp_path, mpnn)
    rows = list(csv.DictReader((out / "sequence_metadata.csv").open()))
    assert list(rows[0]) == [
        "design", "seq", "mpnn_score", "native_residues", "designed_residues",
        "num_changes", "design_region_length_observed", "design_region_spec",
        "corrected_receptor", "receptor_total_length",
    ]
    assert rows[0] == {
        "design": "0", "seq": "0", "mpnn_score": "1.234",
        "native_residues": "xxx", "designed_residues": "DDD",
        "num_changes": "3", "design_region_length_observed": "3",
        "design_region_spec": "3-3",
        "corrected_receptor": "NNNNNDDDCCCCC",
        "receptor_total_length": "13",
    }


@pytest.mark.logic
def test_changes_counts_only_the_designed_residues_that_differ(
        tmp_path, monkeypatch):
    """The native stretch is 'xxx' and the design is 'xDx', so two of the
    three positions match and the count is 1. Counting the anchors too would
    report a change for every position MPNN touched and then had reverted,
    which is the whole receptor."""
    mpnn = _mpnn_tree(tmp_path, {0: [("1.000", "mmmmmxDxmmmmm/eeee")]})
    out = _correct(monkeypatch, tmp_path, mpnn)
    rows = list(csv.DictReader((out / "sequence_metadata.csv").open()))
    assert rows[0]["designed_residues"] == "xDx"
    assert rows[0]["num_changes"] == "1"


@pytest.mark.logic
def test_an_insertion_with_no_native_counterpart_counts_every_residue(
        tmp_path, monkeypatch):
    """Anchors 1-5 and 6-10 are adjacent, so the three designed residues
    replaced nothing and `native_residues` is empty. With no native to
    compare against, every designed residue counts as a change."""
    native = "NNNNN" + "CCCCC"
    mpnn = tmp_path / "mpnn"
    _mpnn_fasta(mpnn / "design_0" / "seqs" / "design_0.fa",
                [("1.000", "mmmmmDDDmmmmm/eeee")])
    out = tmp_path / "out"
    _run(monkeypatch, mpnn_dir=mpnn, output_dir=out, receptor_seq=native,
         effector_seq="EEEE", contigs="A1-5/3-3/A6-10 B", num_designs=1)
    rows = list(csv.DictReader((out / "sequence_metadata.csv").open()))
    assert rows[0]["native_residues"] == ""
    assert rows[0]["designed_residues"] == "DDD"
    assert rows[0]["num_changes"] == "3"


@pytest.mark.logic
def test_only_the_first_num_seqs_sequences_of_a_design_are_written(
        tmp_path, monkeypatch):
    """ProteinMPNN is run at a higher sample count than the pipeline carries
    forward, so the cap is what keeps the two in step. Taking them all would
    put sequences into Boltz that no downstream table expects."""
    mpnn = _mpnn_tree(tmp_path, {0: [
        ("1.000", "mmmmmAAAmmmmm/eeee"),
        ("2.000", "mmmmmBBBmmmmm/eeee"),
        ("3.000", "mmmmmCCCmmmmm/eeee"),
    ]})
    out = _correct(monkeypatch, tmp_path, mpnn, num_seqs=2)
    names = sorted(p.name for p in (out / "af2_fastas").iterdir())
    assert names == ["design_0_seq_0.fasta", "design_0_seq_1.fasta"]
    rows = list(csv.DictReader((out / "sequence_metadata.csv").open()))
    assert [r["designed_residues"] for r in rows] == ["AAA", "BBB"]


@pytest.mark.logic
def test_a_missing_design_is_reported_and_the_rest_still_run(
        tmp_path, monkeypatch, capsys):
    """MPNN can fail on one design out of sixty-four. The run has to carry
    the rest, and it has to say which one it dropped, because a short
    metadata CSV is otherwise indistinguishable from a short design set."""
    mpnn = _mpnn_tree(tmp_path, {1: [("1.000", "mmmmmDDDmmmmm/eeee")]})
    out = _correct(monkeypatch, tmp_path, mpnn, num_designs=2)
    rows = list(csv.DictReader((out / "sequence_metadata.csv").open()))
    assert [r["design"] for r in rows] == ["1"]
    err = capsys.readouterr().out
    assert "design_0.fa not found, skipping" in err


@pytest.mark.logic
def test_a_run_that_corrects_nothing_writes_an_empty_combined_fasta(
        tmp_path, monkeypatch, capsys):
    """The documented empty result. With no MPNN output at all the combined
    FASTA is written and holds a single newline, the metadata CSV is not
    written, and the af2_fastas directory is created and left empty.

    A reader downstream therefore sees an existing but empty FASTA and an
    absent CSV, which is the shape to test for rather than a non-zero exit.
    """
    mpnn = tmp_path / "mpnn"
    mpnn.mkdir()
    out = _correct(monkeypatch, tmp_path, mpnn)
    assert (out / "mpnn_corrected.fasta").read_text() == "\n"
    assert not (out / "sequence_metadata.csv").exists()
    assert (out / "af2_fastas").is_dir()
    assert list((out / "af2_fastas").iterdir()) == []
    assert "Done: 0 corrected sequences" in capsys.readouterr().out


@pytest.mark.logic
def test_correct_mode_needs_an_mpnn_directory(tmp_path, monkeypatch):
    with pytest.raises(AssertionError, match="mpnn_dir"):
        _run(monkeypatch, output_dir=tmp_path / "out",
             receptor_seq=_NATIVE_RECEPTOR, effector_seq="EEEE",
             contigs=_CONTIG, num_designs=1)


@pytest.mark.logic
def test_correct_is_the_default_mode(tmp_path, monkeypatch):
    """No --mode on the command line runs the correction, because that is
    what the stage invokes in the common case. A default of
    gen_fixed_positions would silently produce no sequences."""
    mpnn = _mpnn_tree(tmp_path, {0: [("1.234", "mmmmmDDDmmmmm/eeee")]})
    out = _correct(monkeypatch, tmp_path, mpnn)
    assert (out / "mpnn_corrected.fasta").is_file()


@pytest.mark.logic
def test_a_contig_naming_an_absent_chain_corrects_to_an_empty_receptor(
        tmp_path, monkeypatch, capsys):
    """--receptor_chain Z resolves to no segments, so the walk has nothing to
    consume and the corrected receptor is the empty string. The AF2 FASTA is
    still written, with an empty receptor record.

    Unlike gen_fixed_positions this mode does warn, because the accounting
    check fires: zero accounted residues against a thirteen-residue MPNN
    receptor. That warning is the only signal, and it goes to stdout rather
    than changing the exit code, so a Nextflow stage reads the run as a
    success. Pinned so the warning cannot quietly disappear.
    """
    mpnn = _mpnn_tree(tmp_path, {0: [("1.234", "mmmmmDDDmmmmm/eeee")]})
    out = _correct(monkeypatch, tmp_path, mpnn, receptor_chain="Z")
    assert (out / "af2_fastas" / "design_0_seq_0.fasta").read_text() == (
        ">receptor\n\n>effector\nEEEE\n"
    )
    captured = capsys.readouterr().out
    assert "does not match contig accounting" in captured
    assert "diff=+13" in captured


@pytest.mark.logic
def test_the_receptor_start_offset_shifts_which_native_residues_are_anchors(
        tmp_path, monkeypatch):
    """The contig numbers residues as the PDB does and the native sequence is
    a bare string, so --receptor_start_pdb is what converts between them. The
    receptor here starts at PDB residue 101, and the contig fixes 101-105 and
    109-113, which are offsets 0-4 and 8-12 of the native string.

    Dropping the offset would read the native from position 101 of a
    thirteen-residue string and hand Boltz an empty receptor.
    """
    mpnn = _mpnn_tree(tmp_path, {0: [("1.000", "mmmmmDDDmmmmm/eeee")]})
    out = tmp_path / "out"
    _run(monkeypatch, mpnn_dir=mpnn, output_dir=out,
         receptor_seq=_NATIVE_RECEPTOR, effector_seq="EEEE",
         contigs="A101-105/3-3/A109-113 B", receptor_start_pdb=101,
         num_designs=1)
    rows = list(csv.DictReader((out / "sequence_metadata.csv").open()))
    assert rows[0]["corrected_receptor"] == "NNNNNDDDCCCCC"
    assert rows[0]["native_residues"] == "xxx"


@pytest.mark.logic
def test_two_denovo_segments_split_the_slack_in_a_split_pdb_too(
        tmp_path, monkeypatch):
    """The merged-chain and split-chain branches each carry their own copy of
    the proportional allocation, so the same contig has to give the same
    answer through both. Receptor chain A numbered 101-122 with the effector
    on B: seven de novo residues split 3 and 4 by maximum, leaving 101-105,
    109-113 and 118-122 held.
    """
    pdb = _split_pdb(tmp_path / "design_0.pdb",
                     [("A", 101, 22), ("B", 1, 4)])
    out = tmp_path / "out"
    out.mkdir()
    _run(monkeypatch, mode="gen_fixed_positions", pdb_file=pdb,
         output_dir=out, receptor_seq="N" * 22, effector_seq="EEEE",
         contigs="A1-5/2-4/A9-13/3-5/A18-22 B", num_designs=1)

    assert json.loads((out / "fixed_positions_design_0.jsonl").read_text()) == {
        "design_0": {"B": [1, 2, 3, 4],
                     "A": [101, 102, 103, 104, 105,
                           109, 110, 111, 112, 113,
                           118, 119, 120, 121, 122]}
    }
