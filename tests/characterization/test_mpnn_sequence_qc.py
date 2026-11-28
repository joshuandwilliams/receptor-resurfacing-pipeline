"""Cover the MPNN sequence QC predicates and the filter they drive.

The two predicates decide which designed sequences are thrown away before
anything is paid for in GPU time. Both fail in the expensive direction rather
than the loud one: a predicate that under-reports lets a degenerate sequence
through to structure prediction, and one that over-reports silently discards a
design that was fine.

`main()` wires them to the SEQUENCE_QC stage's three outputs, and adds two
decisions of its own that nothing else in the pipeline revisits: the percent
identity window, and the deduplication that keeps one of two designs whose
corrected receptor sequences came out identical. Which one it keeps is a
choice, so the tests below say which and why.

The thresholds themselves are the contract, so the tests sit on the boundary
rather than well inside it. Row shapes come from
tests/proteinmpnn/example_output_files/sequences/sequence_metadata.csv and the
FASTA shape from that run's qc_fastas/.
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

import mpnn_sequence_qc as qc  # noqa: E402
from fixtures_mpnn_stage import (  # noqa: E402
    af2_fasta,
    qc_row,
    sequence_metadata_csv,
)

# A receptor sequence that clears every check: no run longer than two, no
# residue above 18% of the chain, nothing outside the twenty canonical AAs.
CLEAN_RECEPTOR = "MKTAYIAKQRQISFVKSHFSRQ"

# ── Poly-X runs ──────────────────────────────────────────────────────

@pytest.mark.logic
def test_a_run_longer_than_the_limit_is_caught_and_the_residue_named():
    """Naming the residue is what makes the rejection actionable. 'failed QC'
    on its own does not say whether the sequence is poly-alanine or has one
    stray tract."""
    hit, aa = qc.has_poly_x("MKT" + "A" * 6 + "GG", max_run=5)
    assert hit is True
    assert aa == "A"


@pytest.mark.logic
def test_a_run_exactly_at_the_limit_is_allowed():
    """max_run is the longest run that passes, not the first that fails. The
    predicate tests for aa * (max_run + 1), so 5 A's pass at max_run=5."""
    hit, aa = qc.has_poly_x("MKT" + "A" * 5 + "GG", max_run=5)
    assert hit is False
    assert aa is None


@pytest.mark.logic
def test_one_over_the_limit_fails():
    assert qc.has_poly_x("A" * 6, max_run=5)[0] is True


@pytest.mark.logic
def test_runs_are_detected_for_any_residue_not_just_alanine():
    hit, aa = qc.has_poly_x("MK" + "Q" * 8 + "T", max_run=5)
    assert hit is True and aa == "Q"


@pytest.mark.logic
def test_a_run_split_by_another_residue_is_not_a_run():
    """AAAxAAA is two runs of three, not one of six. Counting the total
    occurrences instead would reject a perfectly ordinary sequence."""
    assert qc.has_poly_x("AAA" + "K" + "AAA", max_run=5)[0] is False


@pytest.mark.logic
def test_an_empty_sequence_has_no_run():
    assert qc.has_poly_x("", max_run=5) == (False, None)


# ── Composition ──────────────────────────────────────────────────────

@pytest.mark.logic
def test_a_non_standard_residue_is_rejected_and_named():
    """X, B, Z and U reach here from a parse that could not resolve a residue.
    Passing one through puts an unresolvable letter into a FASTA that a
    structure predictor will later reject or silently reinterpret."""
    hit, aa = qc.unusual_aa_composition("MKTAX" + "G" * 10)
    assert hit is True and aa == "X"


@pytest.mark.logic
def test_an_ordinary_sequence_passes():
    seq = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"
    assert qc.unusual_aa_composition(seq) == (False, None)


@pytest.mark.logic
def test_a_dominant_residue_above_the_fraction_is_rejected():
    """Half the sequence is tryptophan, which no designed surface should be."""
    seq = "W" * 10 + "MKTAYIAKQR"
    hit, aa = qc.unusual_aa_composition(seq, max_single_aa_frac=0.30)
    assert hit is True and aa == "W"


@pytest.mark.logic
@pytest.mark.parametrize("aa", ["A", "L", "G"])
def test_alanine_leucine_and_glycine_are_exempt_from_the_fraction(aa):
    """The exemption is deliberate. RFdiffusion writes glycine at every de
    novo position and the scaffold is alanine-rich, so applying the fraction
    to these three would reject essentially every design before MPNN has
    assigned a sequence.
    """
    seq = aa * 20 + "MKTYIQR"
    assert qc.unusual_aa_composition(seq, max_single_aa_frac=0.30)[0] is False


@pytest.mark.logic
def test_the_fraction_is_strict_so_exactly_at_the_threshold_passes():
    """The comparison is `frac > max_single_aa_frac`, so a residue at exactly
    the threshold is kept. 3 of 10 is 0.30 and passes."""
    seq = "KKK" + "MTYIQRV"
    assert len(seq) == 10
    assert qc.unusual_aa_composition(seq, max_single_aa_frac=0.30)[0] is False


@pytest.mark.logic
def test_just_over_the_fraction_fails():
    seq = "KKKK" + "MTYIQR"
    assert len(seq) == 10
    hit, aa = qc.unusual_aa_composition(seq, max_single_aa_frac=0.30)
    assert hit is True and aa == "K"


@pytest.mark.logic
def test_an_empty_sequence_does_not_divide_by_zero():
    assert qc.unusual_aa_composition("") == (False, None)


@pytest.mark.logic
def test_a_stricter_fraction_rejects_what_the_default_allows():
    """The flag has to actually reach the comparison."""
    seq = "KKK" + "MTYIQRV"
    assert qc.unusual_aa_composition(seq, max_single_aa_frac=0.30)[0] is False
    assert qc.unusual_aa_composition(seq, max_single_aa_frac=0.10)[0] is True


# ── main(): the filter the predicates drive ──────────────────────────

MISSING = object()


def _run_qc(tmp_path, rows, monkeypatch, *, fastas=None, args=()):
    """Write the stage's two inputs, run main(), return the output directory.

    One af2 FASTA per metadata row, each carrying that row's
    corrected_receptor, which is the arrangement AF2_FASTAS produces.
    `fastas` overrides a named one: a string puts a different receptor in the
    file than in the CSV, None writes a FASTA with no receptor record at all,
    and MISSING leaves the file out.
    """
    meta = sequence_metadata_csv(tmp_path / "sequence_metadata.csv", rows)
    fasta_dir = tmp_path / "af2_fastas"
    fasta_dir.mkdir(exist_ok=True)
    for row in rows:
        name = f"design_{row['design']}_seq_{row['seq']}"
        receptor = (fastas or {}).get(name, row["corrected_receptor"])
        if receptor is MISSING:
            continue
        af2_fasta(fasta_dir / f"{name}.fasta", receptor=receptor)
    out = tmp_path / "out"
    monkeypatch.setattr(sys, "argv", [
        "mpnn_sequence_qc.py",
        "--af2-fastas-dir", str(fasta_dir),
        "--metadata-csv", str(meta),
        "--receptor-seq", CLEAN_RECEPTOR,
        "--output-dir", str(out),
        *args,
    ])
    qc.main()
    return out


def _passed(out: Path) -> list[dict]:
    with open(out / "qc_metadata.csv") as handle:
        return list(csv.DictReader(handle))


def _names(out: Path) -> list[str]:
    return sorted(p.name for p in (out / "qc_fastas").iterdir())


def _report(out: Path) -> str:
    return (out / "qc_report.txt").read_text()


@pytest.mark.logic
def test_a_clean_sequence_passes_and_its_fasta_is_copied(tmp_path, monkeypatch):
    """The three outputs the stage publishes, on the simplest input. The FASTA
    copy is the one that matters downstream: AF2_PREDICT globs qc_fastas/, so a
    row in qc_metadata.csv with no file beside it is a sequence that is
    reported as passing and then never predicted.
    """
    out = _run_qc(tmp_path, [qc_row(design="0", seq="0",
                                    corrected_receptor=CLEAN_RECEPTOR)],
                  monkeypatch)
    assert [(r["design"], r["seq"]) for r in _passed(out)] == [("0", "0")]
    assert _names(out) == ["design_0_seq_0.fasta"]
    assert "Passed QC: 1" in _report(out)
    assert "Failed QC: 0" in _report(out)


@pytest.mark.logic
def test_the_input_columns_are_carried_through_to_qc_metadata(
        tmp_path, monkeypatch):
    """MPNN_DESIGN_REGION_SCORE reads qc_metadata.csv and appends one column to
    it, so every column SEQUENCE_QC was given has to survive."""
    out = _run_qc(tmp_path, [qc_row(corrected_receptor=CLEAN_RECEPTOR)],
                  monkeypatch)
    with open(out / "qc_metadata.csv") as handle:
        header = next(csv.reader(handle))
    assert header == [
        "design", "seq", "mpnn_score", "native_residues", "designed_residues",
        "num_changes", "design_region_length_observed", "design_region_spec",
        "corrected_receptor", "receptor_total_length",
    ]


@pytest.mark.logic
def test_a_missing_fasta_fails_the_sequence_and_names_the_reason(
        tmp_path, monkeypatch):
    """A design whose FASTA never arrived is a plumbing failure, not a bad
    sequence, so the report has to distinguish it from a QC rejection."""
    out = _run_qc(tmp_path,
                  [qc_row(design="0", seq="0"), qc_row(design="0", seq="1")],
                  monkeypatch, fastas={"design_0_seq_0": MISSING})
    assert [(r["design"], r["seq"]) for r in _passed(out)] == [("0", "1")]
    assert "design_0_seq_0: FASTA file missing" in _report(out)


@pytest.mark.logic
def test_the_receptor_checked_is_the_one_in_the_fasta_not_the_one_in_the_csv(
        tmp_path, monkeypatch):
    """The FASTA is what AF2 will actually fold, so it is the sequence QC has
    to judge. Here the CSV column is clean and the FASTA carries a run of six
    alanines, and the sequence must fail. Reading the CSV instead would pass
    a sequence that is about to be predicted as poly-alanine.
    """
    out = _run_qc(
        tmp_path,
        [qc_row(design="0", seq="0", corrected_receptor=CLEAN_RECEPTOR)],
        monkeypatch,
        fastas={"design_0_seq_0": "MKTYIQR" + "A" * 6 + "VSFVKSH"})
    assert _passed(out) == []
    assert "Poly-A run > 5" in _report(out)


@pytest.mark.logic
def test_a_fasta_with_no_receptor_record_falls_back_to_the_csv_column(
        tmp_path, monkeypatch):
    """`af2_fasta(receptor=None)` writes the effector record alone, so the
    header scan never matches and corrected_receptor stays None. The fallback
    is the CSV column, and it has to be the one judged: a dirty CSV value with
    no FASTA receptor must still fail rather than skip the check.
    """
    dirty = "MKTYIQR" + "A" * 6 + "VSFVKSH"
    out = _run_qc(tmp_path,
                  [qc_row(design="0", seq="0", corrected_receptor=dirty)],
                  monkeypatch, fastas={"design_0_seq_0": None})
    assert _passed(out) == []
    assert "Poly-A run > 5" in _report(out)


@pytest.mark.logic
def test_a_run_exactly_at_max_poly_x_passes_through_main(tmp_path, monkeypatch):
    """The boundary, through the CLI rather than the predicate: --max-poly-x 5
    is the longest run that passes, because the predicate searches for
    `aa * (max_run + 1)`."""
    out = _run_qc(
        tmp_path,
        [qc_row(design="0", seq="0",
                corrected_receptor="MKTYIQR" + "A" * 5 + "VSFVKSH")],
        monkeypatch, args=("--max-poly-x", "5"))
    assert len(_passed(out)) == 1


@pytest.mark.logic
def test_one_residue_over_max_poly_x_fails_through_main(tmp_path, monkeypatch):
    out = _run_qc(
        tmp_path,
        [qc_row(design="0", seq="0",
                corrected_receptor="MKTYIQR" + "A" * 6 + "VSFVKSH")],
        monkeypatch, args=("--max-poly-x", "5"))
    assert _passed(out) == []
    assert "Poly-A run > 5" in _report(out)


@pytest.mark.logic
def test_an_unusual_composition_warns_but_does_not_reject(tmp_path, monkeypatch):
    """Composition is advisory and identity is not. A chain that is 31%
    lysine is flagged in the report and still goes to AF2, which is the
    deliberate asymmetry: the fraction has no principled cutoff, so it informs
    rather than filters.
    """
    out = _run_qc(
        tmp_path,
        [qc_row(design="0", seq="0", corrected_receptor="KQKQKQKQKQMTYIVF")],
        monkeypatch)
    assert len(_passed(out)) == 1
    assert "has unusual AA composition" in _report(out)


@pytest.mark.logic
def test_identity_exactly_at_the_lower_bound_passes(tmp_path, monkeypatch):
    """The comparison is `pct_id < min_pct_identity`, so it is strict and a
    sequence exactly at the bound is kept. AAKK against AADD is 2 of 4, so
    50% passes --min-pct-identity 50.
    """
    out = _run_qc(
        tmp_path,
        [qc_row(design="0", seq="0", corrected_receptor=CLEAN_RECEPTOR,
                designed_residues="AAKK", native_residues="AADD")],
        monkeypatch, args=("--min-pct-identity", "50"))
    assert len(_passed(out)) == 1


@pytest.mark.logic
def test_identity_below_the_lower_bound_fails_and_the_report_says_by_how_much(
        tmp_path, monkeypatch):
    """AKKK against ADDD is 1 of 4. Too little identity means the design has
    drifted off the native fold, which is the failure RFDiffusion's contigs
    were supposed to prevent."""
    out = _run_qc(
        tmp_path,
        [qc_row(design="0", seq="0", corrected_receptor=CLEAN_RECEPTOR,
                designed_residues="AKKK", native_residues="ADDD")],
        monkeypatch, args=("--min-pct-identity", "50"))
    assert _passed(out) == []
    assert "% identity 25.0% < 50.0%" in _report(out)


@pytest.mark.logic
def test_identity_exactly_at_the_upper_bound_fails(tmp_path, monkeypatch):
    """The upper comparison is `pct_id >= max_pct_identity`, so it is
    inclusive and the boundary is a rejection, the opposite of the lower one.
    50% identity fails --max-pct-identity 50. The asymmetry is deliberate:
    --max-pct-identity 100 has to reject a sequence identical to the native,
    which is a design that changed nothing.
    """
    out = _run_qc(
        tmp_path,
        [qc_row(design="0", seq="0", corrected_receptor=CLEAN_RECEPTOR,
                designed_residues="AAKK", native_residues="AADD")],
        monkeypatch, args=("--max-pct-identity", "50"))
    assert _passed(out) == []
    assert "% identity 50.0% >= 50.0%" in _report(out)


@pytest.mark.logic
def test_identity_below_the_upper_bound_passes(tmp_path, monkeypatch):
    out = _run_qc(
        tmp_path,
        [qc_row(design="0", seq="0", corrected_receptor=CLEAN_RECEPTOR,
                designed_residues="AKKK", native_residues="ADDD")],
        monkeypatch, args=("--max-pct-identity", "50"))
    assert len(_passed(out)) == 1


@pytest.mark.logic
def test_the_default_window_rejects_only_an_unchanged_design(
        tmp_path, monkeypatch):
    """0.0 to 100.0 is the default, and 100% identity is `>= 100.0`, so a
    design identical to its native is rejected even with the window left
    open. That is the only thing the default window excludes.
    """
    unchanged = qc_row(design="0", seq="0", corrected_receptor=CLEAN_RECEPTOR,
                       designed_residues="AADD", native_residues="AADD")
    changed = qc_row(design="1", seq="0", corrected_receptor=CLEAN_RECEPTOR,
                     designed_residues="AAKK", native_residues="AADD")
    out = _run_qc(tmp_path, [unchanged, changed], monkeypatch)
    assert [r["design"] for r in _passed(out)] == ["1"]
    assert "% identity 100.0% >= 100.0%" in _report(out)


@pytest.mark.logic
def test_the_separator_is_stripped_before_the_identity_comparison(
        tmp_path, monkeypatch):
    """Both strings are pipe-joined per region, and the separators always
    match, so leaving them in would count each one as an identity and inflate
    the percentage. "AA|KK" against "AA|DD" is 2 of 4 and not 3 of 5.
    """
    out = _run_qc(
        tmp_path,
        [qc_row(design="0", seq="0", corrected_receptor=CLEAN_RECEPTOR,
                designed_residues="AA|KK", native_residues="AA|DD")],
        monkeypatch, args=("--min-pct-identity", "50"))
    assert len(_passed(out)) == 1
    assert "% identity" not in _report(out)


@pytest.mark.logic
def test_identity_is_not_checked_when_the_two_strings_differ_in_length(
        tmp_path, monkeypatch):
    """A length-range contig gives a design region shorter or longer than the
    native stretch it replaced, so there is no position-wise comparison to
    make. The check is skipped rather than guessed at, and the sequence
    passes on the other two checks alone.
    """
    out = _run_qc(
        tmp_path,
        [qc_row(design="0", seq="0", corrected_receptor=CLEAN_RECEPTOR,
                designed_residues="AAKKKK", native_residues="AADD")],
        monkeypatch, args=("--min-pct-identity", "99"))
    assert len(_passed(out)) == 1


@pytest.mark.logic
def test_identity_is_not_checked_when_the_native_string_is_empty(
        tmp_path, monkeypatch):
    """native_residues is blank when the input PDB had no gap residues to read
    from, which is the de-novo-only case."""
    out = _run_qc(
        tmp_path,
        [qc_row(design="0", seq="0", corrected_receptor=CLEAN_RECEPTOR,
                designed_residues="AAKK", native_residues="")],
        monkeypatch, args=("--min-pct-identity", "99"))
    assert len(_passed(out)) == 1


# ── Deduplication ────────────────────────────────────────────────────

@pytest.mark.logic
def test_dedup_keeps_the_better_scoring_duplicate_when_it_comes_second(
        tmp_path, monkeypatch):
    """Two designs can independently converge on the same receptor sequence,
    and MPNN score is lower-is-better, so the keeper is the 1.0 and not the
    2.0. The file on disk is the assertion that matters: a kept row whose
    FASTA was deleted, or a dropped row whose FASTA survived, both send the
    wrong sequence to AF2.
    """
    rows = [qc_row(design="0", seq="0", mpnn_score=2.0,
                   corrected_receptor=CLEAN_RECEPTOR),
            qc_row(design="1", seq="0", mpnn_score=1.0,
                   corrected_receptor=CLEAN_RECEPTOR)]
    out = _run_qc(tmp_path, rows, monkeypatch)
    assert [(r["design"], r["mpnn_score"]) for r in _passed(out)] == [("1", "1.0")]
    assert _names(out) == ["design_1_seq_0.fasta"]
    assert "design_0_seq_0: Duplicate receptor seq (kept design_1_seq_0" in _report(out)


@pytest.mark.logic
def test_dedup_keeps_the_better_scoring_duplicate_when_it_comes_first(
        tmp_path, monkeypatch):
    """The mirror case, which takes the other branch of the comparison. The
    answer must not depend on input order."""
    rows = [qc_row(design="0", seq="0", mpnn_score=1.0,
                   corrected_receptor=CLEAN_RECEPTOR),
            qc_row(design="1", seq="0", mpnn_score=2.0,
                   corrected_receptor=CLEAN_RECEPTOR)]
    out = _run_qc(tmp_path, rows, monkeypatch)
    assert [r["design"] for r in _passed(out)] == ["0"]
    assert _names(out) == ["design_0_seq_0.fasta"]
    assert "design_1_seq_0: Duplicate receptor seq (kept design_0_seq_0" in _report(out)


@pytest.mark.logic
def test_dedup_breaks_a_score_tie_in_favour_of_the_earlier_row(
        tmp_path, monkeypatch):
    """`score < prev_score` is strict, so an exact tie keeps the row already
    seen. Equal scores are reachable: two designs that converged on the same
    sequence often converged on the same score too.
    """
    rows = [qc_row(design="0", seq="0", mpnn_score=1.5,
                   corrected_receptor=CLEAN_RECEPTOR),
            qc_row(design="1", seq="0", mpnn_score=1.5,
                   corrected_receptor=CLEAN_RECEPTOR)]
    out = _run_qc(tmp_path, rows, monkeypatch)
    assert [r["design"] for r in _passed(out)] == ["0"]


@pytest.mark.logic
def test_dedup_counts_three_copies_as_two_removals(tmp_path, monkeypatch):
    """The best of the three is kept whatever position it is in. 0.5 is in the
    middle, so both branches of the comparison run in one pass."""
    rows = [qc_row(design="0", seq="0", mpnn_score=2.0,
                   corrected_receptor=CLEAN_RECEPTOR),
            qc_row(design="1", seq="0", mpnn_score=0.5,
                   corrected_receptor=CLEAN_RECEPTOR),
            qc_row(design="2", seq="0", mpnn_score=1.0,
                   corrected_receptor=CLEAN_RECEPTOR)]
    out = _run_qc(tmp_path, rows, monkeypatch)
    assert [r["design"] for r in _passed(out)] == ["1"]
    assert "Duplicates removed: 2" in _report(out)
    assert "Final unique sequences: 1" in _report(out)


@pytest.mark.logic
def test_different_receptor_sequences_are_not_duplicates(tmp_path, monkeypatch):
    rows = [qc_row(design="0", seq="0", mpnn_score=2.0,
                   corrected_receptor=CLEAN_RECEPTOR),
            qc_row(design="1", seq="0", mpnn_score=1.0,
                   corrected_receptor=CLEAN_RECEPTOR[::-1])]
    out = _run_qc(tmp_path, rows, monkeypatch)
    assert sorted(r["design"] for r in _passed(out)) == ["0", "1"]
    assert "Duplicates removed: 0" in _report(out)
    assert "Duplicate sequences removed:" not in _report(out)


@pytest.mark.logic
def test_no_dedup_keeps_both_copies(tmp_path, monkeypatch):
    """--no-dedup exists because the duplicate rows are still distinct
    backbones, and a campaign that wants per-backbone statistics needs them.
    Both FASTAs stay on disk."""
    rows = [qc_row(design="0", seq="0", mpnn_score=2.0,
                   corrected_receptor=CLEAN_RECEPTOR),
            qc_row(design="1", seq="0", mpnn_score=1.0,
                   corrected_receptor=CLEAN_RECEPTOR)]
    out = _run_qc(tmp_path, rows, monkeypatch, args=("--no-dedup",))
    assert sorted(r["design"] for r in _passed(out)) == ["0", "1"]
    assert _names(out) == ["design_0_seq_0.fasta", "design_1_seq_0.fasta"]
    assert "Duplicates removed: 0" in _report(out)


@pytest.mark.logic
@pytest.mark.xfail(strict=False, reason=(
    "Defect: the dedup score read is float(row.get('mpnn_score', '999')), so "
    "the 999 default only applies when the column is absent entirely. A "
    "present-but-empty cell reaches float('') and raises ValueError out of "
    "main(), failing the whole SEQUENCE_QC stage over one unparsed score. "
    "Treating an unreadable score as worst-possible would match the intent of "
    "the 999 default. Not fixed here because it changes which of two "
    "duplicates the stage keeps."))
def test_an_unparseable_mpnn_score_should_not_fail_the_stage(
        tmp_path, monkeypatch):
    rows = [qc_row(design="0", seq="0", mpnn_score=1.0,
                   corrected_receptor=CLEAN_RECEPTOR),
            qc_row(design="1", seq="0", mpnn_score="",
                   corrected_receptor=CLEAN_RECEPTOR)]
    out = _run_qc(tmp_path, rows, monkeypatch)
    assert [r["design"] for r in _passed(out)] == ["0"]


@pytest.mark.logic
def test_an_unparseable_mpnn_score_raises_today(tmp_path, monkeypatch):
    """Pins the behaviour the xfail above describes, so the defect cannot be
    half-fixed without one of the two tests changing."""
    rows = [qc_row(design="0", seq="0", mpnn_score=1.0,
                   corrected_receptor=CLEAN_RECEPTOR),
            qc_row(design="1", seq="0", mpnn_score="",
                   corrected_receptor=CLEAN_RECEPTOR)]
    with pytest.raises(ValueError):
        _run_qc(tmp_path, rows, monkeypatch)


# ── Degenerate input ─────────────────────────────────────────────────

@pytest.mark.logic
def test_nothing_passing_writes_a_marker_file_and_a_header_only_csv(
        tmp_path, monkeypatch):
    """qc_fastas/ is a Nextflow output path, so an empty directory would not
    be staged and the stage would fail on a missing output. The marker file
    keeps the directory non-empty and says why, and AF2_PREDICT's glob then
    finds nothing to fold rather than crashing.
    """
    out = _run_qc(tmp_path, [qc_row(design="0", seq="0")], monkeypatch,
                  fastas={"design_0_seq_0": MISSING})
    assert _names(out) == ["NONE_PASSED_QC.txt"]
    assert _passed(out) == []
    assert "Failed QC: 1" in _report(out)


@pytest.mark.logic
def test_the_empty_qc_metadata_header_is_shorter_than_the_input_header(
        tmp_path, monkeypatch):
    """Characterising a wart rather than asserting an intention. With nothing
    passing, the header is a hardcoded seven columns instead of the ten the
    input carried, so a downstream reader that indexes by name finds
    design_region_length_observed, design_region_spec and
    receptor_total_length missing. Nothing downstream reads those from an
    empty file today, which is why this is pinned rather than changed.
    """
    out = _run_qc(tmp_path, [qc_row(design="0", seq="0")], monkeypatch,
                  fastas={"design_0_seq_0": MISSING})
    assert (out / "qc_metadata.csv").read_text() == (
        "design,seq,mpnn_score,native_residues,designed_residues,"
        "num_changes,corrected_receptor\n")


@pytest.mark.logic
def test_zero_input_rows(tmp_path, monkeypatch):
    """A header-only metadata CSV: ProteinMPNN produced no sequences at all.
    Every output still has to appear."""
    out = _run_qc(tmp_path, [], monkeypatch)
    assert _names(out) == ["NONE_PASSED_QC.txt"]
    assert "Total input sequences: 0" in _report(out)
    assert "Passed QC: 0" in _report(out)


@pytest.mark.logic
def test_main_with_an_absent_metadata_file_raises(tmp_path, monkeypatch):
    """No guard on the input path, which is right for a missing upstream
    file."""
    fasta_dir = tmp_path / "af2_fastas"
    fasta_dir.mkdir()
    monkeypatch.setattr(sys, "argv", [
        "mpnn_sequence_qc.py",
        "--af2-fastas-dir", str(fasta_dir),
        "--metadata-csv", str(tmp_path / "absent.csv"),
        "--receptor-seq", CLEAN_RECEPTOR,
        "--output-dir", str(tmp_path / "out"),
    ])
    with pytest.raises(FileNotFoundError):
        qc.main()


# ── CLI ──────────────────────────────────────────────────────────────

@pytest.mark.logic
@pytest.mark.parametrize("omit", ["--af2-fastas-dir", "--metadata-csv",
                                  "--receptor-seq"])
def test_parse_args_requires_its_three_paths(monkeypatch, omit):
    argv = ["mpnn_sequence_qc.py", "--af2-fastas-dir", "d",
            "--metadata-csv", "m.csv", "--receptor-seq", "MKT"]
    index = argv.index(omit)
    del argv[index:index + 2]
    monkeypatch.setattr(sys, "argv", argv)
    with pytest.raises(SystemExit):
        qc.parse_args()


@pytest.mark.logic
def test_parse_args_defaults_match_the_predicate_defaults(monkeypatch):
    """SEQUENCE_QC passes max_poly_x and the identity window from params but
    never --max-single-aa-frac, so that one comes from here. It has to agree
    with `unusual_aa_composition`'s own default or the two paths into the
    predicate would disagree.
    """
    monkeypatch.setattr(sys, "argv", [
        "mpnn_sequence_qc.py", "--af2-fastas-dir", "d",
        "--metadata-csv", "m.csv", "--receptor-seq", "MKT"])
    args = qc.parse_args()
    assert args.max_poly_x == 5
    assert args.min_pct_identity == 0.0
    assert args.max_pct_identity == 100.0
    assert args.max_single_aa_frac == 0.30
    assert args.no_dedup is False
    assert args.output_dir == "."


@pytest.mark.logic
def test_the_required_receptor_seq_argument_is_never_read():
    """--receptor-seq is required by parse_args and referenced nowhere in
    main(). SEQUENCE_QC interpolates the whole native receptor sequence into
    the command line to supply it, so the wart costs a long argument on every
    invocation. Asserted on the source rather than on behaviour, because
    there is no behaviour to assert.
    """
    source = (BIN / "mpnn_sequence_qc.py").read_text()
    assert '"--receptor-seq", required=True' in source
    assert "args.receptor_seq" not in source
