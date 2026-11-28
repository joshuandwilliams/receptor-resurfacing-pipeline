"""Drive the MPNN top-N selection stage end to end.

`main()` is the MPNN_SELECT_TOP process. It reads the QC-passed metadata CSV,
keeps the N sequences with the lowest MPNN score, and writes the three things
the structure-prediction stage downstream reads:

    top_fastas/           one FASTA per selected sequence
    top_metadata.csv      the selected rows, in selection order
    selection_report.txt  the human-readable summary

The failure mode worth catching is a selection that is well-formed and wrong.
MPNN score is better when lower, so a sort in the wrong direction would hand
the campaign its N worst sequences and every downstream file would look
exactly as it does now. The same goes for the FASTA copies: a sequence present
in top_metadata.csv but missing from top_fastas/ is never predicted, and
nothing in the stage's own output says so.

The tests state the intended selection at the call site and build the cohort
that must produce it.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
BIN = REPO / "bin"
HERE = Path(__file__).resolve().parent
for _p in (str(BIN), str(HERE)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import mpnn_select_top as mst  # noqa: E402
from fixtures_mpnn_stage import af2_fasta, qc_row, sequence_metadata_csv  # noqa: E402


def _cohort(tmp_path: Path, scores, *, with_fastas=True, columns=None):
    """A QC-passed cohort of len(scores) sequences, one per score.

    `scores` maps a (design, seq) pair to that sequence's MPNN score, so the
    ranking a test expects is written out in the call rather than derived.
    Every sequence gets a FASTA in qc_fastas/ unless `with_fastas` says
    otherwise, which is the only state in which a selection can be complete.
    """
    fastas = tmp_path / "qc_fastas"
    rows = []
    for (design, seq), score in scores.items():
        rows.append(qc_row(design=design, seq=seq, mpnn_score=score))
        if with_fastas:
            af2_fasta(fastas / f"design_{design}_seq_{seq}.fasta",
                      receptor=f"MKT{design}{seq}")
    fastas.mkdir(parents=True, exist_ok=True)
    meta = sequence_metadata_csv(tmp_path / "qc_metadata.csv", rows,
                                 columns=columns)
    return fastas, meta


def _run(monkeypatch, fastas: Path, meta: Path, out: Path, top_n: int) -> int:
    monkeypatch.setattr(sys, "argv", [
        "mpnn_select_top.py",
        "--qc-fastas-dir", str(fastas),
        "--qc-metadata", str(meta),
        "--top-n", str(top_n),
        "--output-dir", str(out),
    ])
    try:
        mst.main()
    except SystemExit as exc:
        return int(exc.code or 0)
    return 0


def _selected(out: Path):
    """The (design, seq) pairs in top_metadata.csv, in file order."""
    with open(out / "top_metadata.csv", newline="") as f:
        return [(r["design"], r["seq"]) for r in csv.DictReader(f)]


def _metadata(out: Path):
    with open(out / "top_metadata.csv", newline="") as f:
        reader = csv.DictReader(f)
        return list(reader.fieldnames or []), list(reader)


def _copied(out: Path):
    return sorted(p.name for p in (out / "top_fastas").iterdir())


# ── Selection by score ───────────────────────────────────────────────

@pytest.mark.logic
def test_the_lowest_scoring_sequences_are_the_ones_selected(monkeypatch, tmp_path):
    """MPNN score is a negative log likelihood, so lower is better. Sorting
    the other way would hand the campaign its worst two sequences and every
    output file would still be well-formed."""
    fastas, meta = _cohort(tmp_path, {
        ("0", "0"): 1.50,
        ("0", "1"): 0.50,   # best
        ("1", "0"): 2.50,
        ("1", "1"): 0.90,   # second best
    })
    out = tmp_path / "out"
    assert _run(monkeypatch, fastas, meta, out, 2) == 0
    assert _selected(out) == [("0", "1"), ("1", "1")]


@pytest.mark.logic
def test_the_selection_is_written_in_ascending_score_order(monkeypatch, tmp_path):
    """top_metadata.csv is the ranking, not just the membership. The report
    numbers the rows 1..N off this order, so an unsorted file would print a
    rank that does not match the score beside it."""
    fastas, meta = _cohort(tmp_path, {("0", "0"): 3.0, ("0", "1"): 1.0,
                                      ("0", "2"): 2.0})
    out = tmp_path / "out"
    _run(monkeypatch, fastas, meta, out, 3)
    _fields, rows = _metadata(out)
    assert [float(r["mpnn_score"]) for r in rows] == [1.0, 2.0, 3.0]


@pytest.mark.logic
def test_a_tie_keeps_the_order_the_metadata_had(monkeypatch, tmp_path):
    """Three sequences on the same score, so only a stable sort can give a
    reproducible selection. An unstable one would make a rerun of the same
    cohort pick a different sequence."""
    fastas, meta = _cohort(tmp_path, {("0", "0"): 1.0, ("0", "1"): 1.0,
                                      ("0", "2"): 1.0})
    out = tmp_path / "out"
    _run(monkeypatch, fastas, meta, out, 2)
    assert _selected(out) == [("0", "0"), ("0", "1")]


@pytest.mark.logic
def test_asking_for_more_than_the_cohort_holds_selects_everything(monkeypatch, tmp_path):
    """--top-n is a cap, not a quota. A cohort of 2 asked for 10 must give 2
    rather than padding or failing."""
    fastas, meta = _cohort(tmp_path, {("0", "0"): 1.0, ("0", "1"): 2.0})
    out = tmp_path / "out"
    assert _run(monkeypatch, fastas, meta, out, 10) == 0
    assert len(_selected(out)) == 2


@pytest.mark.logic
def test_asking_for_none_selects_none_from_a_full_cohort(monkeypatch, tmp_path):
    fastas, meta = _cohort(tmp_path, {("0", "0"): 1.0, ("0", "1"): 2.0})
    out = tmp_path / "out"
    assert _run(monkeypatch, fastas, meta, out, 0) == 0
    assert _selected(out) == []
    assert _copied(out) == []


@pytest.mark.logic
def test_a_cohort_with_no_score_column_at_all_still_selects(monkeypatch, tmp_path):
    """The sort key falls back to 999 for an absent mpnn_score column, so a
    metadata CSV written before the scoring stage ran still produces a
    selection instead of taking the stage down."""
    fastas, meta = _cohort(tmp_path, {("0", "0"): 1.0, ("0", "1"): 2.0},
                           columns=["design", "seq"])
    out = tmp_path / "out"
    assert _run(monkeypatch, fastas, meta, out, 1) == 0
    assert _selected(out) == [("0", "0")]


# ── The FASTA copies ─────────────────────────────────────────────────

@pytest.mark.logic
def test_one_fasta_is_copied_per_selected_sequence_and_no_others(monkeypatch, tmp_path):
    """A FASTA copied for an unselected sequence would be predicted and
    reported as a top-N design that was never chosen."""
    fastas, meta = _cohort(tmp_path, {("0", "0"): 1.0, ("0", "1"): 2.0,
                                      ("0", "2"): 3.0})
    out = tmp_path / "out"
    _run(monkeypatch, fastas, meta, out, 2)
    assert _copied(out) == ["design_0_seq_0.fasta", "design_0_seq_1.fasta"]


@pytest.mark.logic
def test_the_copy_is_byte_identical_to_the_qc_fasta(monkeypatch, tmp_path):
    """The stage copies rather than regenerates, so the sequence predicted is
    the sequence QC passed."""
    fastas, meta = _cohort(tmp_path, {("0", "0"): 1.0})
    out = tmp_path / "out"
    _run(monkeypatch, fastas, meta, out, 1)
    src = (fastas / "design_0_seq_0.fasta").read_text()
    assert (out / "top_fastas" / "design_0_seq_0.fasta").read_text() == src


@pytest.mark.logic
def test_a_missing_fasta_warns_and_the_stage_still_finishes(monkeypatch, tmp_path,
                                                            capsys):
    """The row stays in top_metadata.csv with no FASTA behind it, so the
    warning on stdout is the only record that the sequence will never be
    predicted. Exiting instead would lose the rest of the selection."""
    fastas, meta = _cohort(tmp_path, {("0", "0"): 1.0, ("0", "1"): 2.0},
                           with_fastas=False)
    out = tmp_path / "out"
    assert _run(monkeypatch, fastas, meta, out, 2) == 0
    printed = capsys.readouterr().out
    assert printed.count("WARNING") == 2
    assert "design_0_seq_0.fasta not found" in printed
    assert _copied(out) == []
    assert len(_selected(out)) == 2


@pytest.mark.logic
def test_the_output_tree_is_created_rather_than_assumed(monkeypatch, tmp_path):
    """Nextflow runs the stage in a bare work directory."""
    fastas, meta = _cohort(tmp_path, {("0", "0"): 1.0})
    out = tmp_path / "a" / "b" / "c"
    assert _run(monkeypatch, fastas, meta, out, 1) == 0
    assert (out / "top_fastas").is_dir()


# ── top_metadata.csv ─────────────────────────────────────────────────

@pytest.mark.logic
def test_every_input_column_survives_into_the_selection(monkeypatch, tmp_path):
    """The downstream stages read corrected_receptor and designed_residues off
    this file. Dropping a column here would make the selected sequences
    unusable even though the selection itself was right."""
    fastas, meta = _cohort(tmp_path, {("0", "0"): 1.0})
    out = tmp_path / "out"
    _run(monkeypatch, fastas, meta, out, 1)
    with open(meta, newline="") as f:
        expected = list(csv.DictReader(f).fieldnames or [])
    fields, rows = _metadata(out)
    assert fields == expected
    assert rows[0]["corrected_receptor"]


# ── selection_report.txt ─────────────────────────────────────────────

@pytest.mark.logic
def test_the_report_counts_both_the_cohort_and_the_selection(monkeypatch, tmp_path):
    """4 QC-passed, top 2 requested, 2 selected. Reporting the cohort size as
    the selection size is how an over-wide --top-n would hide itself."""
    fastas, meta = _cohort(tmp_path, {("0", "0"): 1.0, ("0", "1"): 2.0,
                                      ("0", "2"): 3.0, ("0", "3"): 4.0})
    out = tmp_path / "out"
    _run(monkeypatch, fastas, meta, out, 2)
    report = (out / "selection_report.txt").read_text()
    assert "Total QC-passed: 4" in report
    assert "Selected top 2: 2" in report


@pytest.mark.logic
def test_the_report_states_the_score_range_of_the_selection_not_the_cohort(
        monkeypatch, tmp_path):
    """The cohort spans 0.5 to 9.0 and the selected two span 0.5 to 1.5, so a
    report quoting 0.500 - 9.000 would be describing sequences it did not
    select. Three decimal places, as the format string asks for."""
    fastas, meta = _cohort(tmp_path, {("0", "0"): 0.5, ("0", "1"): 1.5,
                                      ("0", "2"): 9.0})
    out = tmp_path / "out"
    _run(monkeypatch, fastas, meta, out, 2)
    assert "Score range: 0.500 - 1.500" in (out / "selection_report.txt").read_text()


@pytest.mark.logic
def test_the_report_lists_the_selection_numbered_from_one_in_rank_order(
        monkeypatch, tmp_path):
    fastas, meta = _cohort(tmp_path, {("0", "0"): 2.0, ("1", "3"): 1.0})
    out = tmp_path / "out"
    _run(monkeypatch, fastas, meta, out, 2)
    report = (out / "selection_report.txt").read_text()
    assert "  1. design_1_seq_3 (score: 1.0)" in report
    assert "  2. design_0_seq_0 (score: 2.0)" in report


# ── The empty cohort ─────────────────────────────────────────────────

@pytest.mark.logic
def test_an_empty_cohort_still_writes_every_published_file(monkeypatch, tmp_path):
    """Nextflow publishes a fixed file list, so a stage that exits without
    writing fails with nothing in the log explaining why."""
    fastas, meta = _cohort(tmp_path, {})
    out = tmp_path / "out"
    assert _run(monkeypatch, fastas, meta, out, 3) == 0
    assert (out / "top_fastas").is_dir()
    assert (out / "top_metadata.csv").is_file()
    assert (out / "selection_report.txt").is_file()


@pytest.mark.logic
def test_the_empty_report_omits_the_score_range_rather_than_inventing_one(
        monkeypatch, tmp_path):
    """min() of nothing raises, so the range line is skipped. A zero range
    would read as a cohort of identical sequences rather than no cohort."""
    fastas, meta = _cohort(tmp_path, {})
    out = tmp_path / "out"
    _run(monkeypatch, fastas, meta, out, 3)
    report = (out / "selection_report.txt").read_text()
    assert "Total QC-passed: 0" in report
    assert "Selected top 3: 0" in report
    assert "Score range" not in report


@pytest.mark.logic
def test_the_empty_cohort_writes_a_three_column_header_not_the_input_header(
        monkeypatch, tmp_path):
    """Today's behaviour, pinned. The populated path copies the input's ten
    columns through; the empty path writes a hardcoded three. A reader that
    indexes top_metadata.csv by column therefore sees a different schema
    depending on whether anything was selected."""
    fastas, meta = _cohort(tmp_path, {})
    out = tmp_path / "out"
    _run(monkeypatch, fastas, meta, out, 3)
    assert (out / "top_metadata.csv").read_text() == "design,seq,mpnn_score\n"


@pytest.mark.logic
@pytest.mark.xfail(strict=False, reason=(
    "Defect: the empty-selection branch writes a hardcoded "
    "'design,seq,mpnn_score' header instead of the input CSV's columns, so "
    "top_metadata.csv changes schema when nothing is selected. Not fixed "
    "here because these are characterization tests and bin/ is unchanged; "
    "the fix is to carry the input fieldnames through both branches."))
def test_the_empty_cohort_should_keep_the_input_header(monkeypatch, tmp_path):
    fastas, meta = _cohort(tmp_path, {})
    out = tmp_path / "out"
    _run(monkeypatch, fastas, meta, out, 3)
    with open(meta, newline="") as f:
        expected = list(csv.DictReader(f).fieldnames or [])
    assert _metadata(out)[0] == expected


# ── Defects pinned ───────────────────────────────────────────────────

@pytest.mark.logic
def test_a_blank_score_takes_the_stage_down(monkeypatch, tmp_path):
    """Today's behaviour, pinned. The 999 fallback only applies when the
    mpnn_score KEY is absent. A present-but-empty value reaches float("")
    and raises, which is the shape a sequence whose scoring step failed
    actually has, because the QC stage writes the column either way."""
    fastas, meta = _cohort(tmp_path, {("0", "0"): 1.0, ("0", "1"): ""})
    out = tmp_path / "out"
    with pytest.raises(ValueError):
        _run(monkeypatch, fastas, meta, out, 1)


@pytest.mark.logic
@pytest.mark.xfail(strict=False, reason=(
    "Defect: a blank mpnn_score raises ValueError and loses the whole "
    "cohort, where an absent column is tolerated via the '999' default. "
    "Not fixed here because bin/ is unchanged under a characterization "
    "pass; the fix is to coerce a blank to the same 999 sentinel so one "
    "unscored sequence sorts last instead of failing the stage."))
def test_a_blank_score_should_sort_last_like_an_absent_column(monkeypatch, tmp_path):
    fastas, meta = _cohort(tmp_path, {("0", "0"): 1.0, ("0", "1"): ""})
    out = tmp_path / "out"
    assert _run(monkeypatch, fastas, meta, out, 2) == 0
    assert _selected(out) == [("0", "0"), ("0", "1")]


@pytest.mark.logic
def test_a_negative_top_n_selects_all_but_the_worst(monkeypatch, tmp_path):
    """Today's behaviour, pinned. --top-n is sliced as rows[:n] with no
    validation, so -1 is a Python negative slice: it drops the worst
    sequence and keeps the other two. The report then prints
    'Selected top -1: 2'."""
    fastas, meta = _cohort(tmp_path, {("0", "0"): 1.0, ("0", "1"): 2.0,
                                      ("0", "2"): 3.0})
    out = tmp_path / "out"
    assert _run(monkeypatch, fastas, meta, out, -1) == 0
    assert _selected(out) == [("0", "0"), ("0", "1")]
    assert "Selected top -1: 2" in (out / "selection_report.txt").read_text()


@pytest.mark.logic
@pytest.mark.xfail(strict=False, reason=(
    "Defect: a negative --top-n is interpreted as a negative slice and "
    "quietly selects all-but-N sequences rather than being rejected. Not "
    "fixed here because bin/ is unchanged under a characterization pass; "
    "the fix is to reject top_n < 0 in parse_args."))
def test_a_negative_top_n_should_be_rejected(monkeypatch, tmp_path):
    fastas, meta = _cohort(tmp_path, {("0", "0"): 1.0, ("0", "1"): 2.0,
                                      ("0", "2"): 3.0})
    out = tmp_path / "out"
    assert _run(monkeypatch, fastas, meta, out, -1) == 2
