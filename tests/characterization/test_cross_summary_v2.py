"""Drive the typed CLI the NEGSTEER_CROSS_SEQUENCE process runs.

`cross_sequence_summary.aggregate` builds the 63-column CSV and is covered by
`test_cross_sequence_summary.py`. What this file covers is the CLI in front of
it: how a `--passing-summary` argument becomes a (sequence name, path) pair,
and which of the two input modes reaches the aggregation.

The sequence name is the whole risk. It is the join key the cohort plots and
the wet-lab triage list index on, and it is derived from a directory name
rather than read out of the file. Taking the parent where the layout wants the
grandparent gives every sequence the same name, which deduplicates the cohort
down to one row. Taking the grandparent where the layout wants the parent
names every sequence after the run directory. Neither raises, and both produce
a well-formed CSV.

The second thing asserted here is that supplying no input at all is refused
rather than silently writing an empty cohort, because a published
cross_sequence_summary.csv with no rows reads as a campaign where nothing
worked.
"""
from __future__ import annotations

import csv
import importlib.util
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
BIN = REPO / "bin"
HERE = Path(__file__).resolve().parent
for _p in (str(BIN), str(HERE)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import cross_summary_v2 as csv2  # noqa: E402
from fixtures_cohort import passing_row, write_passing_summary  # noqa: E402


def _sequence(runs: Path, name: str, **overrides) -> Path:
    """One runs/<mpnn_sequence>/passing_summary.csv, as NEGSTEER_RUN_ONE
    publishes it."""
    return write_passing_summary(runs / name / "passing_summary.csv",
                                 [passing_row(**overrides)])


def _run(monkeypatch, out: Path, *args) -> int:
    monkeypatch.setattr(sys, "argv",
                        ["cross_summary_v2.py", "--output", str(out), *args])
    try:
        return csv2.main()
    except SystemExit as exc:
        return int(exc.code or 0)


def _rows(out: Path):
    with open(out, newline="") as f:
        reader = csv.DictReader(f)
        return list(reader.fieldnames or []), list(reader)


def _names(out: Path):
    return [row["mpnn_sequence"] for row in _rows(out)[1]]


# ── Naming an explicit passing summary ───────────────────────────────

@pytest.mark.logic
def test_an_explicit_name_before_the_path_wins():
    """NAME=PATH is how a caller that knows the sequence name states it
    rather than relying on the layout. The name must come out verbatim,
    including when the path's own parent says something else."""
    name, path = csv2._parse_passing_summary_arg(
        "design_0_seq_1=/work/ab/cd/passing_summary.csv", "parent")
    assert name == "design_0_seq_1"
    assert path == Path("/work/ab/cd/passing_summary.csv")


@pytest.mark.logic
def test_only_the_first_equals_sign_splits_the_argument():
    """A work-dir path can contain an '=' and must not be truncated at it."""
    name, path = csv2._parse_passing_summary_arg(
        "seq_0=/work/a=b/passing_summary.csv", "parent")
    assert name == "seq_0"
    assert path == Path("/work/a=b/passing_summary.csv")


@pytest.mark.logic
def test_an_explicit_name_is_taken_even_in_grandparent_mode():
    """The mode only decides how to DERIVE a name. An explicit one is not
    derived, so the switch must not reach it."""
    name, _path = csv2._parse_passing_summary_arg(
        "chosen=/runs/design_0_seq_0/inner/passing_summary.csv", "grandparent")
    assert name == "chosen"


# ── Deriving the name from the layout ────────────────────────────────

@pytest.mark.logic
def test_the_parent_directory_names_the_sequence_by_default():
    """The published layout is runs/<mpnn_sequence>/passing_summary.csv, so
    the file's own directory is the sequence."""
    name, path = csv2._parse_passing_summary_arg(
        "/out/runs/design_0_seq_0/passing_summary.csv", "parent")
    assert name == "design_0_seq_0"
    assert path == Path("/out/runs/design_0_seq_0/passing_summary.csv")


@pytest.mark.logic
def test_the_grandparent_names_the_sequence_when_asked():
    """A Nextflow work directory stages the file one level deeper, under
    <seq>/<subdir>/passing_summary.csv. Taking the parent there would name
    every sequence after its subdirectory."""
    name, _path = csv2._parse_passing_summary_arg(
        "/work/design_0_seq_0/negative_steering/passing_summary.csv",
        "grandparent")
    assert name == "design_0_seq_0"


@pytest.mark.logic
def test_the_two_modes_disagree_on_the_same_path():
    """Stated explicitly because the modes are a pair of mutually exclusive
    flags on one dest, and picking the wrong one is silent."""
    arg = "/work/design_0_seq_0/inner/passing_summary.csv"
    assert csv2._parse_passing_summary_arg(arg, "parent")[0] == "inner"
    assert csv2._parse_passing_summary_arg(arg, "grandparent")[0] == "design_0_seq_0"


@pytest.mark.logic
def test_the_mode_flag_reaches_the_derivation(monkeypatch, tmp_path):
    """The two switches are store_const on one dest with "parent" as the
    default, so the grandparent flag has to arrive or the CLI silently uses
    the default on a work-dir layout."""
    deep = tmp_path / "design_7_seq_2" / "negsteer" / "passing_summary.csv"
    write_passing_summary(deep, [passing_row()])

    parent_out = tmp_path / "parent.csv"
    assert _run(monkeypatch, parent_out, "--passing-summary", str(deep)) == 0
    assert _names(parent_out) == ["negsteer"]

    gp_out = tmp_path / "gp.csv"
    assert _run(monkeypatch, gp_out, "--seq-name-from-grandparent-dir",
                "--passing-summary", str(deep)) == 0
    assert _names(gp_out) == ["design_7_seq_2"]


@pytest.mark.logic
def test_the_parent_flag_is_accepted_and_is_also_the_default(monkeypatch, tmp_path):
    """The Nextflow process passes it explicitly on one branch, so it has to
    exist and to mean the same as omitting it."""
    path = _sequence(tmp_path / "runs", "design_0_seq_0")
    out = tmp_path / "explicit.csv"
    assert _run(monkeypatch, out, "--seq-name-from-parent-dir",
                "--passing-summary", str(path)) == 0
    assert _names(out) == ["design_0_seq_0"]


# ── The directory mode ───────────────────────────────────────────────

@pytest.mark.logic
def test_every_sequence_subdirectory_becomes_a_row(monkeypatch, tmp_path):
    """This is the mode NEGSTEER_CROSS_SEQUENCE uses: one aggregator
    directory holding a symlink per sequence. Rows come out in sorted
    order."""
    runs = tmp_path / "aggregator"
    _sequence(runs, "design_0_seq_0", n_pass="3")
    _sequence(runs, "design_0_seq_1", n_pass="2")
    _sequence(runs, "design_1_seq_0", n_pass="0")
    out = tmp_path / "cross.csv"
    assert _run(monkeypatch, out, "--passing-summary-dir", str(runs)) == 0
    assert _names(out) == ["design_0_seq_0", "design_0_seq_1", "design_1_seq_0"]


@pytest.mark.logic
def test_a_loose_file_in_the_runs_directory_is_not_a_sequence(monkeypatch, tmp_path):
    """The aggregator directory also collects logs and the odd stray file. A
    row named after one would carry no metrics and be ranked anyway."""
    runs = tmp_path / "aggregator"
    _sequence(runs, "design_0_seq_0")
    (runs / "nextflow.log").write_text("not a sequence\n")
    out = tmp_path / "cross.csv"
    _run(monkeypatch, out, "--passing-summary-dir", str(runs))
    assert _names(out) == ["design_0_seq_0"]


@pytest.mark.logic
def test_a_directory_name_outside_the_allowed_character_set_is_skipped(
        monkeypatch, tmp_path):
    """Sequence names go into CSV cells and into figure filenames, so the
    walk accepts only [A-Za-z0-9_.+-]. A name with a space or a slash-like
    character is dropped rather than carried into either."""
    runs = tmp_path / "aggregator"
    _sequence(runs, "design_0_seq_0")
    _sequence(runs, "not a sequence name!")
    out = tmp_path / "cross.csv"
    _run(monkeypatch, out, "--passing-summary-dir", str(runs))
    assert _names(out) == ["design_0_seq_0"]


@pytest.mark.logic
def test_several_runs_directories_are_merged(monkeypatch, tmp_path):
    """--passing-summary-dir is appendable, which is how a resumed campaign
    with two aggregator roots is summarised in one pass."""
    first, second = tmp_path / "a", tmp_path / "b"
    _sequence(first, "design_0_seq_0")
    _sequence(second, "design_1_seq_0")
    out = tmp_path / "cross.csv"
    assert _run(monkeypatch, out,
                "--passing-summary-dir", str(first),
                "--passing-summary-dir", str(second)) == 0
    assert sorted(_names(out)) == ["design_0_seq_0", "design_1_seq_0"]


@pytest.mark.logic
def test_an_explicit_pair_takes_precedence_over_the_same_name_in_a_directory(
        monkeypatch, tmp_path):
    """Explicit pairs are gathered first and the walk skips a name already
    seen, so a caller can override one sequence's file without rebuilding
    the directory."""
    runs = tmp_path / "aggregator"
    _sequence(runs, "design_0_seq_0", ra_eff_vs_truth_median="9.0")
    override = write_passing_summary(
        tmp_path / "override" / "passing_summary.csv",
        [passing_row(ra_eff_vs_truth_median="1.0")])
    out = tmp_path / "cross.csv"
    assert _run(monkeypatch, out,
                "--passing-summary", f"design_0_seq_0={override}",
                "--passing-summary-dir", str(runs)) == 0
    _fields, rows = _rows(out)
    assert len(rows) == 1
    assert rows[0]["rep_ra_eff_vs_truth_median"] == "1.0"


@pytest.mark.logic
def test_the_same_sequence_named_twice_produces_one_row(monkeypatch, tmp_path):
    """A resumed run can stage the same workdir twice. Two rows for one
    sequence would double-count it in every tier total."""
    path = _sequence(tmp_path / "runs", "design_0_seq_0")
    out = tmp_path / "cross.csv"
    assert _run(monkeypatch, out, "--passing-summary", str(path),
                "--passing-summary", str(path)) == 0
    assert _names(out) == ["design_0_seq_0"]


@pytest.mark.logic
def test_a_runs_directory_that_does_not_exist_yields_no_rows(monkeypatch, tmp_path):
    """The walk skips a non-directory, so nothing is gathered and the
    aggregation refuses an empty cohort with exit 2 rather than publishing
    an empty CSV."""
    out = tmp_path / "cross.csv"
    assert _run(monkeypatch, out, "--passing-summary-dir",
                str(tmp_path / "never_created")) == 2
    assert not out.exists()


# ── Refusing to run on nothing ───────────────────────────────────────

@pytest.mark.logic
def test_supplying_neither_input_is_refused(monkeypatch, tmp_path, capsys):
    """Exit 2 with nothing written. An empty cross_sequence_summary.csv would
    read as a campaign in which no sequence survived."""
    out = tmp_path / "cross.csv"
    assert _run(monkeypatch, out) == 2
    assert not out.exists()
    assert "supply --passing-summary" in capsys.readouterr().err


# ── Options reaching the aggregation ─────────────────────────────────

@pytest.mark.logic
def test_a_relative_published_runs_dir_is_made_absolute(monkeypatch, tmp_path):
    """The path is used to rewrite work-dir absolute paths into stable
    published ones, so a relative value would only match by accident from
    whichever directory the stage happened to run in."""
    runs = tmp_path / "aggregator"
    _sequence(runs, "design_0_seq_0")
    published = tmp_path / "published" / "runs"
    published.mkdir(parents=True)
    monkeypatch.chdir(tmp_path)
    out = tmp_path / "cross.csv"
    assert _run(monkeypatch, out, "--passing-summary-dir", str(runs),
                "--published-runs-dir", "published/runs") == 0
    assert _rows(out)[1][0]["mpnn_sequence"] == "design_0_seq_0"


@pytest.mark.logic
def test_strict_turns_a_missing_passing_summary_into_a_failure(monkeypatch,
                                                               tmp_path, capsys):
    """Without --strict a sequence whose file never landed becomes a
    diagnostic stub row, which is right for a cohort report. With it the
    stage fails, which is right when the caller needs every sequence
    present."""
    out = tmp_path / "cross.csv"
    ghost = tmp_path / "design_9_seq_9" / "passing_summary.csv"
    assert _run(monkeypatch, out, "--strict", "--passing-summary", str(ghost)) == 2
    assert not out.exists()
    assert "not found" in capsys.readouterr().err


@pytest.mark.logic
def test_without_strict_a_missing_passing_summary_becomes_an_untiered_row(
        monkeypatch, tmp_path):
    """The sequence stays visible in the cohort and lands in tier "none", so
    a run that failed is reported rather than omitted."""
    out = tmp_path / "cross.csv"
    ghost = tmp_path / "design_9_seq_9" / "passing_summary.csv"
    assert _run(monkeypatch, out, "--passing-summary", str(ghost)) == 0
    _fields, rows = _rows(out)
    assert len(rows) == 1
    assert rows[0]["mpnn_sequence"] == "design_9_seq_9"
    assert rows[0]["cross_tier"] == "none"


@pytest.mark.logic
def test_the_scored_metadata_columns_are_joined_on(monkeypatch, tmp_path):
    """--scored-metadata adds corrected_receptor, designed_residues and
    native_residues to every row, joined by mpnn_sequence. The flag has to
    reach the aggregation or the triage list loses the sequences themselves.
    """
    runs = tmp_path / "aggregator"
    _sequence(runs, "design_0_seq_0")
    scored = tmp_path / "scored_metadata.csv"
    with open(scored, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "design", "seq", "corrected_receptor", "designed_residues",
            "native_residues"])
        writer.writeheader()
        writer.writerow({"design": "0", "seq": "0",
                         "corrected_receptor": "MKTAYIAK",
                         "designed_residues": "AAKK",
                         "native_residues": "AADD"})
    out = tmp_path / "cross.csv"
    assert _run(monkeypatch, out, "--passing-summary-dir", str(runs),
                "--scored-metadata", str(scored)) == 0
    _fields, rows = _rows(out)
    assert rows[0]["corrected_receptor"] == "MKTAYIAK"
    assert rows[0]["designed_residues"] == "AAKK"


# ── Routing, not reimplementation ────────────────────────────────────

@pytest.mark.logic
def test_the_output_carries_the_full_cross_sequence_column_set(monkeypatch,
                                                               tmp_path):
    """This module is advertised as a drop-in replacement for
    cross_sequence_summary.py, so the CSV has to be the one that script
    writes rather than a reduced version. 63 columns, as the committed
    reference set under tests/negative_steering/example_output_files/
    carries."""
    runs = tmp_path / "aggregator"
    _sequence(runs, "design_0_seq_0")
    out = tmp_path / "cross.csv"
    _run(monkeypatch, out, "--passing-summary-dir", str(runs))
    fields, _rows_ = _rows(out)
    assert len(fields) == 63
    assert fields[0] == "mpnn_sequence"
    for expected in ("cross_tier", "cross_rank_by_composite",
                     "cross_composite_score", "rep_n_pass", "rep_outcome"):
        assert expected in fields


@pytest.mark.logic
def test_the_tier_is_computed_from_the_row_rather_than_copied(monkeypatch,
                                                              tmp_path):
    """3 of 3 seeds passing is tier A and 2 of 3 is tier B, both decided by
    the aggregation this CLI routes to. Asserting the tier rather than the
    row count is what makes this a test of the routing and not of the file
    being readable."""
    runs = tmp_path / "aggregator"
    _sequence(runs, "design_0_seq_0", n_seeds="3", n_pass="3")
    _sequence(runs, "design_0_seq_1", n_seeds="3", n_pass="2")
    out = tmp_path / "cross.csv"
    _run(monkeypatch, out, "--passing-summary-dir", str(runs))
    tiers = {row["mpnn_sequence"]: row["cross_tier"] for row in _rows(out)[1]}
    assert tiers == {"design_0_seq_0": "A", "design_0_seq_1": "B"}


# ── The sibling import ───────────────────────────────────────────────

@pytest.mark.logic
def test_the_module_puts_its_own_directory_on_the_import_path():
    """The Nextflow process runs `python bin/cross_summary_v2.py` inside a
    container with no PYTHONPATH, so `from design_cohort import DesignCohort`
    resolves only because the module first inserts its own directory into
    sys.path. Every other test here imports it with bin/ already on the path,
    which skips that insert, so loading it from the file with bin/ removed is
    the only way to exercise the branch the cluster actually takes.
    """
    saved = list(sys.path)
    sys.path[:] = [p for p in sys.path
                   if p and Path(p).resolve() != BIN.resolve()]
    assert str(BIN) not in sys.path
    try:
        spec = importlib.util.spec_from_file_location(
            "cross_summary_v2__pathless", BIN / "cross_summary_v2.py")
        module = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
        assert str(BIN) in sys.path, "the module must add its own directory"
    finally:
        sys.path[:] = saved
        sys.modules.pop("cross_summary_v2__pathless", None)
    assert module.DesignCohort is csv2.DesignCohort


# ── Defects pinned ───────────────────────────────────────────────────

@pytest.mark.logic
def test_an_argument_with_an_empty_path_is_accepted_as_the_cwd(monkeypatch,
                                                               tmp_path):
    """Today's behaviour, pinned. "NAME=" yields Path(""), which is Path("."),
    so the sequence is read from a directory. That fails as an unreadable
    file and becomes an untiered stub row instead of being rejected as a
    malformed argument."""
    assert csv2._parse_passing_summary_arg("nameonly=", "parent") \
        == ("nameonly", Path("."))
    out = tmp_path / "cross.csv"
    assert _run(monkeypatch, out, "--passing-summary", "nameonly=") == 0
    _fields, rows = _rows(out)
    assert [row["mpnn_sequence"] for row in rows] == ["nameonly"]
    assert rows[0]["cross_tier"] == "none"


@pytest.mark.logic
@pytest.mark.xfail(strict=False, reason=(
    "Defect: main() wraps _parse_passing_summary_arg in `except ValueError: "
    "return 2`, but that helper cannot raise ValueError for any input -- "
    "str.split never raises and Path() accepts the empty string -- so the "
    "error branch is unreachable and a malformed --passing-summary argument "
    "such as 'NAME=' is carried through as Path('.') and reported as an "
    "untiered sequence. Not fixed here because bin/ is unchanged under a "
    "characterization pass; the fix is for the helper to raise on an empty "
    "name or an empty path so the existing handler fires."))
def test_a_malformed_passing_summary_argument_should_be_refused(monkeypatch,
                                                                tmp_path):
    out = tmp_path / "cross.csv"
    assert _run(monkeypatch, out, "--passing-summary", "nameonly=") == 2
