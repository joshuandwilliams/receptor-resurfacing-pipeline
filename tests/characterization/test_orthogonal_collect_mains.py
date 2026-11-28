"""Drive the Rosetta score collector and the orthogonal merge end to end.

`tests/characterization/test_orthogonal_merge_and_collect.py` covers the filter
cascade, the float coercion and the score-file parser. What it leaves is
`main()` for both: the four files the collector writes into its working
directory, and the sixteen columns the merge appends.

The two stages have opposite failure styles, and both are pinned here.

The collector gates on shape complementarity and copies the winners into
`passing/`. Its quiet failure is a design with no usable Sc being carried
forward, so the nan path is asserted as a fail rather than as a pass.

The merge drops nothing. Every input row is written out, and a survivor with no
summary CSV behind it comes out flagged `sc_missing,bsa_missing` and gated,
because a missing Rosetta run is not evidence that the design was good.
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

import merge_orthogonal_metrics as mom  # noqa: E402
import rosetta_filter_collect as rfc  # noqa: E402


def _write_csv(path: Path, fieldnames, rows) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    return path


def _read_csv(path: Path):
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        return list(reader.fieldnames or []), list(reader)


# ═════════════════════════════════════════════════════════════════════
# rosetta_filter_collect
# ═════════════════════════════════════════════════════════════════════

_FULL_COLUMNS = ("sc_value dG_separated dSASA_int dG_separated/dSASAx100 "
                 "packstat delta_unsatHbonds nres_int description")


def _score_file(score_dir: Path, stem: str, values: str,
                columns: str = _FULL_COLUMNS) -> Path:
    score_dir.mkdir(parents=True, exist_ok=True)
    p = score_dir / f"interface_scores_{stem}.sc"
    p.write_text(f"SEQUENCE: AAAA\nSCORE: {columns}\nSCORE: {values}\n")
    return p


def _run_rfc(monkeypatch, cwd: Path, score_dir: Path, pdb_dir: Path,
             *, sc_threshold: float = 0.62) -> int:
    """Run the collector in `cwd`, which is where it writes its four outputs."""
    monkeypatch.chdir(cwd)
    monkeypatch.setattr(sys, "argv", [
        "rosetta_filter_collect.py",
        "--score-dir", str(score_dir),
        "--pdb-dir", str(pdb_dir),
        "--sc-threshold", str(sc_threshold),
    ])
    try:
        rfc.main()
    except SystemExit as e:
        return int(e.code or 0)
    return 0


def _collector_outputs(cwd: Path):
    metrics = json.loads((cwd / "rosetta_filter_metrics.json").read_text())
    summary = json.loads((cwd / "rosetta_filter_summary.json").read_text())
    passing = [ln for ln in
               (cwd / "rosetta_passing_designs.txt").read_text().splitlines() if ln]
    return metrics, summary, passing


@pytest.mark.logic
def test_one_design_collects_into_one_fully_known_metrics_entry(monkeypatch, tmp_path):
    """Every field is a rounding of a value written into the score file, at
    the per-column precision the committed rosetta_filter_metrics.json
    carries. nres_int is an int rather than a float because the downstream
    plots label it as a residue count."""
    scores, pdbs, cwd = tmp_path / "s", tmp_path / "p", tmp_path / "work"
    cwd.mkdir()
    _score_file(scores, "design_0",
                "0.70001 -45.6789 1234.56 -3.69876 0.65432 3.44 27 design_0")
    pdbs.mkdir()
    (pdbs / "design_0.pdb").write_text("ATOM  design 0\n")

    assert _run_rfc(monkeypatch, cwd, scores, pdbs) == 0
    metrics, summary, passing = _collector_outputs(cwd)

    assert metrics["sc_threshold"] == 0.62
    assert metrics["designs"] == [{
        "design": "design_0.pdb",
        "design_stem": "design_0",
        "sc_value": 0.7,
        "dG_separated": -45.679,
        "dSASA_int": 1234.6,
        "dG_dSASA_density": -3.6988,
        "packstat": 0.6543,
        "delta_unsatHbonds": 3.4,
        "nres_int": 27,
        "passes_filter": True,
    }]
    assert passing == ["design_0.pdb"]
    assert summary == {"n_total": 1, "n_passing": 1, "n_filtered": 0,
                       "sc_threshold": 0.62, "filtered_designs": []}
    assert (cwd / "passing" / "design_0.pdb").read_text() == "ATOM  design 0\n"


@pytest.mark.logic
@pytest.mark.parametrize("sc,passes", [
    ("0.619", False), ("0.62", True), ("0.621", True),
])
def test_the_sc_threshold_is_not_strict_so_the_boundary_passes(
        monkeypatch, tmp_path, sc, passes):
    """The comparison is `sc_value >= threshold`, so a design sitting exactly
    on the published 0.62 is kept. A strict comparison would make the
    documented threshold one the design can never attain."""
    scores, pdbs, cwd = tmp_path / "s", tmp_path / "p", tmp_path / "work"
    cwd.mkdir()
    pdbs.mkdir()
    _score_file(scores, "design_0", f"{sc} -40.0 900.0 -4.0 0.6 2.0 25 design_0")
    (pdbs / "design_0.pdb").write_text("x\n")

    assert _run_rfc(monkeypatch, cwd, scores, pdbs, sc_threshold=0.62) == 0
    metrics, summary, passing = _collector_outputs(cwd)
    assert metrics["designs"][0]["passes_filter"] is passes
    assert bool(passing) is passes
    assert summary["n_passing"] == (1 if passes else 0)
    assert (cwd / "passing" / "design_0.pdb").exists() is passes


@pytest.mark.logic
def test_a_design_with_no_sc_value_fails_rather_than_passes(monkeypatch, tmp_path):
    """Sc is the only gate at this stage, so a nan Sc is no evidence at all.
    `nan >= threshold` is False either way, but the code says so explicitly
    rather than relying on that, and the entry records None so the plots do
    not draw a zero."""
    scores, pdbs, cwd = tmp_path / "s", tmp_path / "p", tmp_path / "work"
    cwd.mkdir()
    pdbs.mkdir()
    _score_file(scores, "design_0", "-40.0 900.0 25 design_0",
                columns="dG_separated dSASA_int nres_int description")
    (pdbs / "design_0.pdb").write_text("x\n")

    assert _run_rfc(monkeypatch, cwd, scores, pdbs) == 0
    metrics, summary, passing = _collector_outputs(cwd)
    entry = metrics["designs"][0]
    assert entry["sc_value"] is None
    assert entry["passes_filter"] is False
    assert entry["packstat"] is None, "an absent column is None, never 0.0"
    assert passing == []
    assert summary["filtered_designs"] == ["design_0.pdb"]


@pytest.mark.logic
def test_a_mixed_cohort_splits_into_passing_and_filtered(monkeypatch, tmp_path):
    """Only the passing PDBs are copied. The filtered ones stay named in the
    summary so the next stage can report what was dropped."""
    scores, pdbs, cwd = tmp_path / "s", tmp_path / "p", tmp_path / "work"
    cwd.mkdir()
    pdbs.mkdir()
    for stem, sc in (("design_0", "0.70"), ("design_1", "0.30"),
                     ("design_2", "0.80")):
        _score_file(scores, stem, f"{sc} -40.0 900.0 -4.0 0.6 2.0 25 {stem}")
        (pdbs / f"{stem}.pdb").write_text(stem)

    assert _run_rfc(monkeypatch, cwd, scores, pdbs) == 0
    metrics, summary, passing = _collector_outputs(cwd)
    assert [d["design"] for d in metrics["designs"]] == [
        "design_0.pdb", "design_1.pdb", "design_2.pdb"]
    assert passing == ["design_0.pdb", "design_2.pdb"]
    assert summary["filtered_designs"] == ["design_1.pdb"]
    assert sorted(p.name for p in (cwd / "passing").iterdir()) == [
        "design_0.pdb", "design_2.pdb"]


@pytest.mark.logic
def test_no_score_files_writes_empty_outputs_and_exits_clean(monkeypatch, tmp_path):
    """Nextflow would fail on a missing output file, so the empty case still
    writes all four. Exiting 0 is what lets a cohort where every Rosetta job
    died be diagnosed from the summary rather than from the scheduler."""
    scores, pdbs, cwd = tmp_path / "s", tmp_path / "p", tmp_path / "work"
    scores.mkdir()
    pdbs.mkdir()
    cwd.mkdir()

    assert _run_rfc(monkeypatch, cwd, scores, pdbs) == 0
    metrics, summary, passing = _collector_outputs(cwd)
    assert metrics == {"designs": [], "sc_threshold": 0.62}
    assert summary == {"n_total": 0, "n_passing": 0, "n_filtered": 0,
                       "sc_threshold": 0.62}
    assert passing == []
    assert (cwd / "passing").is_dir()


@pytest.mark.logic
def test_a_score_file_with_no_data_row_is_skipped_not_counted(
        monkeypatch, tmp_path, capsys):
    """A header with nothing under it is what a killed InterfaceAnalyzer
    leaves. Counting it as a design would put a row of nulls in the metrics
    JSON and a zero in every plot."""
    scores, pdbs, cwd = tmp_path / "s", tmp_path / "p", tmp_path / "work"
    cwd.mkdir()
    pdbs.mkdir()
    scores.mkdir()
    (scores / "interface_scores_design_0.sc").write_text(
        f"SEQUENCE: AAAA\nSCORE: {_FULL_COLUMNS}\n")
    _score_file(scores, "design_1", "0.70 -40.0 900.0 -4.0 0.6 2.0 25 design_1")
    (pdbs / "design_1.pdb").write_text("x\n")

    assert _run_rfc(monkeypatch, cwd, scores, pdbs) == 0
    metrics, summary, _passing = _collector_outputs(cwd)
    assert [d["design_stem"] for d in metrics["designs"]] == ["design_1"]
    assert summary["n_total"] == 1
    assert "No scores in interface_scores_design_0.sc" in capsys.readouterr().out


@pytest.mark.logic
def test_a_filename_with_no_design_stem_is_skipped(monkeypatch, tmp_path, capsys):
    """`interface_scores_.sc` matches the glob but yields no stem, so there is
    no PDB name to attribute the scores to. Guessing one would copy the wrong
    structure into passing/."""
    scores, pdbs, cwd = tmp_path / "s", tmp_path / "p", tmp_path / "work"
    cwd.mkdir()
    pdbs.mkdir()
    scores.mkdir()
    (scores / "interface_scores_.sc").write_text(
        f"SCORE: {_FULL_COLUMNS}\nSCORE: 0.9 -40.0 900.0 -4.0 0.6 2.0 25 x\n")

    assert _run_rfc(monkeypatch, cwd, scores, pdbs) == 0
    metrics, _summary, _passing = _collector_outputs(cwd)
    assert metrics["designs"] == []
    assert "Could not parse design name" in capsys.readouterr().out


@pytest.mark.logic
def test_a_passing_design_whose_pdb_is_absent_is_warned_about(
        monkeypatch, tmp_path, capsys):
    """The score file and the PDB come from two different Nextflow channels,
    so they can disagree. The design stays in the metrics and the passing
    list, and the missing copy is reported rather than crashing the stage."""
    scores, pdbs, cwd = tmp_path / "s", tmp_path / "p", tmp_path / "work"
    cwd.mkdir()
    pdbs.mkdir()
    _score_file(scores, "design_0", "0.70 -40.0 900.0 -4.0 0.6 2.0 25 design_0")

    assert _run_rfc(monkeypatch, cwd, scores, pdbs) == 0
    _metrics, summary, passing = _collector_outputs(cwd)
    assert passing == ["design_0.pdb"]
    assert summary["n_passing"] == 1
    assert list((cwd / "passing").iterdir()) == []
    assert "not found for passing design" in capsys.readouterr().out


@pytest.mark.logic
def test_a_cohort_where_nothing_passes_says_so_on_stderr(monkeypatch, tmp_path, capsys):
    """The whole point of the stage is to hand on survivors. Zero of them is
    the operator's problem to see, so it goes to stderr where .command.err
    picks it up."""
    scores, pdbs, cwd = tmp_path / "s", tmp_path / "p", tmp_path / "work"
    cwd.mkdir()
    pdbs.mkdir()
    _score_file(scores, "design_0", "0.10 -40.0 900.0 -4.0 0.6 2.0 25 design_0")

    assert _run_rfc(monkeypatch, cwd, scores, pdbs) == 0
    assert "No designs passed the Sc filter" in capsys.readouterr().err


# ═════════════════════════════════════════════════════════════════════
# merge_orthogonal_metrics
# ═════════════════════════════════════════════════════════════════════

AF3_SUMMARY = {
    "seq_name": "design_1_seq_1",
    "af3_nomsa_best_ra_eff": "2.100",
    "af3_nomsa_mean_ra_eff": "3.400",
    "af3_nomsa_best_iptm": "0.780",
    "af3_nomsa_mean_iptm": "0.710",
    "af3_nomsa_n_correct_interface": "4",
    "af3_nomsa_total_predictions": "5",
    "af3_nomsa_failures": "",
}
BIO_SUMMARY = {
    "seq_name": "design_1_seq_1",
    "bsa": "912.50",
    "interface_hbonds": "11",
    "biophysical_failures": "",
}
ROSETTA_SUMMARY = {
    "seq_name": "design_1_seq_1",
    "sc": "0.6651",
    "rosetta_ddg": "-59.739",
    "rosetta_failures": "",
}


def _summaries(tmp_path: Path, *, af3=AF3_SUMMARY, bio=BIO_SUMMARY,
               ros=ROSETTA_SUMMARY) -> tuple[str, str, str]:
    for name, row in (("af3", af3), ("bio", bio), ("ros", ros)):
        _write_csv(tmp_path / name / "s.csv", list(row), [row])
    return (str(tmp_path / "af3" / "*.csv"),
            str(tmp_path / "bio" / "*.csv"),
            str(tmp_path / "ros" / "*.csv"))


def _run_mom(monkeypatch, in_csv: Path, out_csv: Path, globs,
             *, af3_ra_max: float = 5.0, sc_min: float = 0.300,
             bsa_min: float = 400.0, plddt_min: float = 0.75) -> int:
    monkeypatch.setattr(sys, "argv", [
        "merge_orthogonal_metrics.py",
        "--input-csv", str(in_csv),
        "--af3-summaries-glob", globs[0],
        "--biophysical-summaries-glob", globs[1],
        "--rosetta-summaries-glob", globs[2],
        "--filter-af3-ra-max", str(af3_ra_max),
        "--filter-sc-min", str(sc_min),
        "--filter-bsa-min", str(bsa_min),
        "--filter-plddt-min", str(plddt_min),
        "--output-csv", str(out_csv),
    ])
    return mom.main()


@pytest.mark.logic
def test_the_three_streams_merge_into_one_fully_known_row(monkeypatch, tmp_path):
    """Every value traces to a named input: the first two columns from the
    input CSV, thirteen from the three summary CSVs, and the last two from
    the filter cascade. The row passes with no flags, so the merge is
    exercised without the filter logic masking a mismapped column."""
    globs = _summaries(tmp_path)
    in_csv = _write_csv(tmp_path / "in.csv", ["mpnn_sequence", "cross_tier"],
                        [{"mpnn_sequence": "design_1_seq_1", "cross_tier": "A"}])
    out_csv = tmp_path / "out" / "survivors.csv"
    assert _run_mom(monkeypatch, in_csv, out_csv, globs) == 0

    header, rows = _read_csv(out_csv)
    assert header == ["mpnn_sequence", "cross_tier"] + mom.ALL_NEW
    assert rows == [{
        "mpnn_sequence": "design_1_seq_1",
        "cross_tier": "A",
        "af3_nomsa_best_ra_eff": "2.100",
        "af3_nomsa_mean_ra_eff": "3.400",
        "af3_nomsa_best_iptm": "0.780",
        "af3_nomsa_mean_iptm": "0.710",
        "af3_nomsa_n_correct_interface": "4",
        "af3_nomsa_total_predictions": "5",
        "af3_nomsa_failures": "",
        "bsa": "912.50",
        "interface_hbonds": "11",
        "biophysical_failures": "",
        "sc": "0.6651",
        "rosetta_ddg": "-59.739",
        "rosetta_failures": "",
        "orthogonal_flags": "",
        "passes_orthogonal_filters": "1",
    }]


@pytest.mark.logic
def test_a_survivor_with_no_summaries_behind_it_is_gated_not_passed(
        monkeypatch, tmp_path):
    """A Rosetta or FreeSASA job that died leaves no summary CSV. Absent is
    not the same as passing, so the row comes out flagged on both gating
    metrics and gated. Carrying it forward as a pass would put an unmeasured
    design into the gBlock selection."""
    globs = _summaries(tmp_path)
    in_csv = _write_csv(tmp_path / "in.csv", ["mpnn_sequence"],
                        [{"mpnn_sequence": "design_99_seq_1"}])
    out_csv = tmp_path / "out.csv"
    assert _run_mom(monkeypatch, in_csv, out_csv, globs) == 0

    _header, rows = _read_csv(out_csv)
    assert rows[0]["sc"] == "" and rows[0]["bsa"] == ""
    assert rows[0]["passes_orthogonal_filters"] == "0"
    flags = rows[0]["orthogonal_flags"].split(",")
    assert set(flags) == {"af3_nomsa_missing", "sc_missing", "bsa_missing"}


@pytest.mark.logic
def test_no_row_is_ever_dropped_whatever_its_flags(monkeypatch, tmp_path):
    """The merge is the last table before the cohort views, which report
    failures as well as survivors. Dropping the failing rows would make the
    hitrate denominator wrong."""
    globs = _summaries(tmp_path)
    in_csv = _write_csv(
        tmp_path / "in.csv", ["mpnn_sequence"],
        [{"mpnn_sequence": "design_1_seq_1"},
         {"mpnn_sequence": "design_2_seq_1"},
         {"mpnn_sequence": "design_3_seq_1"}])
    out_csv = tmp_path / "out.csv"
    assert _run_mom(monkeypatch, in_csv, out_csv, globs) == 0

    _header, rows = _read_csv(out_csv)
    assert [r["mpnn_sequence"] for r in rows] == [
        "design_1_seq_1", "design_2_seq_1", "design_3_seq_1"]
    assert [r["passes_orthogonal_filters"] for r in rows] == ["1", "0", "0"]


@pytest.mark.logic
def test_a_summary_for_a_sequence_not_in_the_input_is_ignored(monkeypatch, tmp_path):
    """The summaries are keyed on seq_name and the input on mpnn_sequence.
    They name the same workdir, so a summary with no matching input row must
    not be attached to whichever row happened to be first."""
    globs = _summaries(tmp_path, ros={**ROSETTA_SUMMARY,
                                      "seq_name": "design_77_seq_1"})
    in_csv = _write_csv(tmp_path / "in.csv", ["mpnn_sequence"],
                        [{"mpnn_sequence": "design_1_seq_1"}])
    out_csv = tmp_path / "out.csv"
    assert _run_mom(monkeypatch, in_csv, out_csv, globs) == 0

    _header, rows = _read_csv(out_csv)
    assert rows[0]["sc"] == ""
    assert "sc_missing" in rows[0]["orthogonal_flags"]


@pytest.mark.logic
def test_the_boundary_values_of_both_gates_are_kept(monkeypatch, tmp_path):
    """Sc and BSA are compared with `value < minimum`, so a design sitting
    exactly on either published threshold passes. The AF3 comparison is the
    other way round, `value >= maximum`, so a ra_eff exactly on its limit is
    flagged, and the flag still does not gate."""
    globs = _summaries(
        tmp_path,
        af3={**AF3_SUMMARY, "af3_nomsa_best_ra_eff": "5.0"},
        bio={**BIO_SUMMARY, "bsa": "400"},
        ros={**ROSETTA_SUMMARY, "sc": "0.300"})
    in_csv = _write_csv(tmp_path / "in.csv", ["mpnn_sequence"],
                        [{"mpnn_sequence": "design_1_seq_1"}])
    out_csv = tmp_path / "out.csv"
    assert _run_mom(monkeypatch, in_csv, out_csv, globs,
                    af3_ra_max=5.0, sc_min=0.300, bsa_min=400.0) == 0

    _header, rows = _read_csv(out_csv)
    assert rows[0]["passes_orthogonal_filters"] == "1"
    assert rows[0]["orthogonal_flags"] == "af3_nomsa_ra_eff_too_high:5.00"


@pytest.mark.logic
def test_a_rerun_does_not_duplicate_the_appended_columns(monkeypatch, tmp_path):
    """The stage is idempotent over its own output. A second `sc` column
    would make the CSV unreadable by every downstream DictReader."""
    globs = _summaries(tmp_path)
    fieldnames = ["mpnn_sequence"] + mom.ALL_NEW
    stale = {c: "stale" for c in mom.ALL_NEW}
    stale["mpnn_sequence"] = "design_1_seq_1"
    in_csv = _write_csv(tmp_path / "in.csv", fieldnames, [stale])
    out_csv = tmp_path / "out.csv"
    assert _run_mom(monkeypatch, in_csv, out_csv, globs) == 0

    header, rows = _read_csv(out_csv)
    assert header == fieldnames
    assert rows[0]["sc"] == "0.6651", "the stale value must be overwritten"
    assert rows[0]["passes_orthogonal_filters"] == "1"


@pytest.mark.logic
def test_a_header_only_input_gives_a_header_only_output(monkeypatch, tmp_path):
    globs = _summaries(tmp_path)
    in_csv = _write_csv(tmp_path / "in.csv", ["mpnn_sequence"], [])
    out_csv = tmp_path / "out.csv"
    assert _run_mom(monkeypatch, in_csv, out_csv, globs) == 0

    header, rows = _read_csv(out_csv)
    assert rows == []
    assert header == ["mpnn_sequence"] + mom.ALL_NEW


@pytest.mark.logic
def test_an_input_with_no_header_at_all_is_an_error_exit(monkeypatch, tmp_path):
    """An empty file means the upstream stage wrote nothing. Exiting 1 stops
    the cascade rather than writing a survivors table with no survivors in
    it, which reads as a cohort that simply failed."""
    globs = _summaries(tmp_path)
    in_csv = tmp_path / "empty.csv"
    in_csv.write_text("")
    out_csv = tmp_path / "out.csv"
    assert _run_mom(monkeypatch, in_csv, out_csv, globs) == 1
    assert not out_csv.exists()


@pytest.mark.logic
def test_globs_matching_nothing_leave_every_metric_blank(monkeypatch, tmp_path):
    """The three streams are optional in the sense that the stage still runs
    without them. What it must not do is report the result as a pass."""
    globs = (str(tmp_path / "none_a" / "*.csv"),
             str(tmp_path / "none_b" / "*.csv"),
             str(tmp_path / "none_c" / "*.csv"))
    in_csv = _write_csv(tmp_path / "in.csv", ["mpnn_sequence"],
                        [{"mpnn_sequence": "design_1_seq_1"}])
    out_csv = tmp_path / "out.csv"
    assert _run_mom(monkeypatch, in_csv, out_csv, globs) == 0

    _header, rows = _read_csv(out_csv)
    assert all(rows[0][c] == "" for c in mom.AF3_COLS + mom.BIO_COLS
               + mom.ROSETTA_COLS)
    assert rows[0]["passes_orthogonal_filters"] == "0"
