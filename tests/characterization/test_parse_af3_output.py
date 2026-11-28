"""Drive the AF3-no-MSA post-processing CLI the AF3_PARSE_OUTPUT process runs.

`af3_confidence.AF3ConfidenceAggregate` walks the prediction tree and is
covered by `test_af3_confidence.py`. What this file covers is the CLI around
it: the eight-column CSV that is the only thing the orthogonal-metrics merge
downstream ever sees.

Two failure modes, both invisible once the CSV is written.

`af3_nomsa_n_correct_interface` is the one column the aggregate cannot fill on
its own, because the count needs a threshold. The CLI builds a `_Thresholds`
shim for it, and a shim carrying the wrong attribute would make the count
silently zero while every other column stayed right.

And the column set is fixed with `extrasaction="ignore"`, so a renamed column
upstream would be dropped rather than raising, leaving a blank where a metric
was. Both are asserted by value here.

Every expected number is arithmetic over `af3_output_tree`, which places each
prediction's effector at a stated displacement from the ground truth, so a
prediction's receptor-aligned effector RMSD is known before the parser runs.
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

pytest.importorskip("gemmi")

import parse_af3_output as pao  # noqa: E402
from fixtures_runtree import af3_output_tree  # noqa: E402

# The column order the AF3_PARSE_OUTPUT process publishes, read by
# merge_orthogonal_metrics.
EXPECTED_COLUMNS = [
    "seq_name",
    "af3_nomsa_best_ra_eff",
    "af3_nomsa_mean_ra_eff",
    "af3_nomsa_best_iptm",
    "af3_nomsa_mean_iptm",
    "af3_nomsa_n_correct_interface",
    "af3_nomsa_total_predictions",
    "af3_nomsa_failures",
]


def _run(monkeypatch, truth: Path, out_dir: Path, csv_path: Path, **extra) -> int:
    argv = ["parse_af3_output.py",
            "--seq-name", extra.pop("seq_name", "design_0_seq_0"),
            "--af3-output-dir", str(out_dir),
            "--ground-truth", str(truth),
            "--output-csv", str(csv_path)]
    for key, value in extra.items():
        argv += [f"--{key.replace('_', '-')}", str(value)]
    monkeypatch.setattr(sys, "argv", argv)
    try:
        return pao.main()
    except SystemExit as exc:
        return int(exc.code or 0)


def _row(csv_path: Path):
    """The single data row, plus the header, from the summary CSV."""
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        fields = list(reader.fieldnames or [])
        rows = list(reader)
    assert len(rows) == 1, "one survivor means exactly one row"
    return fields, rows[0]


# ── The aggregate columns ────────────────────────────────────────────

@pytest.mark.logic
def test_every_column_of_a_two_prediction_run_is_a_known_number(monkeypatch,
                                                                tmp_path):
    """Ground truth puts the effector 4 A above the receptor. Sample 0 is the
    same complex, so its ra_eff is 0.0; sample 1 moves the effector to 7 A,
    so its ra_eff is exactly 7 - 4 = 3.0. Best is 0.000 and the mean is
    (0.0 + 3.0) / 2 = 1.500. iPTMs are 0.80 and 0.60, so best is 0.8000 and
    the mean is 0.7000. The threshold is 5 A and both are below it, so
    n_correct is 2 of 2 and no failure tag is emitted.
    """
    truth, out_dir = af3_output_tree(
        tmp_path / "af3", predictions=[(4.0, 0.80), (7.0, 0.60)],
        truth_separation=4.0)
    summary = tmp_path / "af3_nomsa_summary.csv"
    assert _run(monkeypatch, truth, out_dir, summary,
                ra_eff_threshold=5.0) == 0

    fields, row = _row(summary)
    assert fields == EXPECTED_COLUMNS
    assert row == {
        "seq_name": "design_0_seq_0",
        "af3_nomsa_best_ra_eff": "0.000",
        "af3_nomsa_mean_ra_eff": "1.500",
        "af3_nomsa_best_iptm": "0.8000",
        "af3_nomsa_mean_iptm": "0.7000",
        "af3_nomsa_n_correct_interface": "2",
        "af3_nomsa_total_predictions": "2",
        "af3_nomsa_failures": "",
    }


@pytest.mark.logic
def test_the_best_ra_eff_is_the_minimum_and_the_best_iptm_the_maximum(monkeypatch,
                                                                     tmp_path):
    """ra_eff is an error, so best is lowest; iPTM is a confidence, so best is
    highest. Reducing both the same way is the mistake this catches. Three
    predictions at 4, 10 and 6 A against a 4 A truth give ra_eff 0, 6 and 2,
    so best is 0.000 and mean is (0 + 6 + 2) / 3 = 2.667.
    """
    truth, out_dir = af3_output_tree(
        tmp_path / "af3",
        predictions=[(4.0, 0.30), (10.0, 0.90), (6.0, 0.60)],
        truth_separation=4.0)
    summary = tmp_path / "s.csv"
    _run(monkeypatch, truth, out_dir, summary)
    _fields, row = _row(summary)
    assert row["af3_nomsa_best_ra_eff"] == "0.000"
    assert row["af3_nomsa_mean_ra_eff"] == "2.667"
    assert row["af3_nomsa_best_iptm"] == "0.9000"
    assert row["af3_nomsa_mean_iptm"] == "0.6000"   # (0.3 + 0.9 + 0.6) / 3


@pytest.mark.logic
def test_the_seq_name_is_carried_through_as_given(monkeypatch, tmp_path):
    """The merge joins on this column, so a transformed name joins to
    nothing and the design loses its AF3 metrics silently."""
    truth, out_dir = af3_output_tree(tmp_path / "af3",
                                     predictions=[(4.0, 0.50)])
    summary = tmp_path / "s.csv"
    _run(monkeypatch, truth, out_dir, summary, seq_name="input_control_polyA")
    assert _row(summary)[1]["seq_name"] == "input_control_polyA"


# ── The threshold-dependent count ────────────────────────────────────

@pytest.mark.logic
def test_the_threshold_decides_the_correct_interface_count(monkeypatch, tmp_path):
    """Three predictions with ra_eff 0, 3 and 6 A. The count is strictly
    below the threshold, so at 5 A two qualify and at 2 A only one does. A
    shim whose attribute name did not match would give zero in both."""
    truth, out_dir = af3_output_tree(
        tmp_path / "af3",
        predictions=[(4.0, 0.50), (7.0, 0.50), (10.0, 0.50)],
        truth_separation=4.0)

    wide = tmp_path / "wide.csv"
    _run(monkeypatch, truth, out_dir, wide, ra_eff_threshold=5.0)
    assert _row(wide)[1]["af3_nomsa_n_correct_interface"] == "2"

    tight = tmp_path / "tight.csv"
    _run(monkeypatch, truth, out_dir, tight, ra_eff_threshold=2.0)
    assert _row(tight)[1]["af3_nomsa_n_correct_interface"] == "1"


@pytest.mark.logic
def test_the_default_threshold_is_five_angstroms(monkeypatch, tmp_path):
    """The default has to match the pipeline's own
    orthogonal_filter_af3_ra_max, because a run that omits the flag must not
    count differently from one that passes the campaign's value. A prediction
    at 4.5 A counts and one at 5.5 A does not."""
    truth, out_dir = af3_output_tree(
        tmp_path / "af3", predictions=[(8.5, 0.50), (9.5, 0.50)],
        truth_separation=4.0)   # ra_eff 4.5 and 5.5
    summary = tmp_path / "s.csv"
    _run(monkeypatch, truth, out_dir, summary)
    _fields, row = _row(summary)
    assert row["af3_nomsa_n_correct_interface"] == "1"
    assert row["af3_nomsa_total_predictions"] == "2"


@pytest.mark.logic
def test_the_count_is_strictly_below_the_threshold(monkeypatch, tmp_path):
    """A prediction exactly at the threshold does not count. The orthogonal
    merge treats this column as evidence the pose was reproduced, so the
    boundary belongs on the conservative side."""
    truth, out_dir = af3_output_tree(tmp_path / "af3",
                                     predictions=[(9.0, 0.50)],
                                     truth_separation=4.0)   # ra_eff 5.0
    summary = tmp_path / "s.csv"
    _run(monkeypatch, truth, out_dir, summary, ra_eff_threshold=5.0)
    _fields, row = _row(summary)
    assert row["af3_nomsa_best_ra_eff"] == "5.000"
    assert row["af3_nomsa_n_correct_interface"] == "0"


# ── Absent and broken input ──────────────────────────────────────────

@pytest.mark.logic
def test_a_missing_output_directory_writes_a_flagged_row_and_exits_zero(
        monkeypatch, tmp_path):
    """An AF3 job that never produced a directory must still publish a CSV,
    because Nextflow expects the file and the merge expects the row. The
    metric columns are blank rather than zero, so a failed prediction cannot
    read as a perfect one, and the failure tag names the directory."""
    truth, _out_dir = af3_output_tree(tmp_path / "af3",
                                      predictions=[(4.0, 0.50)])
    ghost = tmp_path / "never_ran"
    summary = tmp_path / "s.csv"
    assert _run(monkeypatch, truth, ghost, summary) == 0

    _fields, row = _row(summary)
    assert row["af3_nomsa_best_ra_eff"] == ""
    assert row["af3_nomsa_mean_ra_eff"] == ""
    assert row["af3_nomsa_best_iptm"] == ""
    assert row["af3_nomsa_total_predictions"] == "0"
    assert row["af3_nomsa_n_correct_interface"] == "0"
    assert row["af3_nomsa_failures"] == f"af3_output_dir_missing:{ghost}"


@pytest.mark.logic
def test_an_empty_output_directory_is_flagged_as_finding_no_structures(
        monkeypatch, tmp_path):
    """A directory that exists but holds no CIF is a different failure from
    one that is absent, and the two tags distinguish them in the merged CSV."""
    truth, _out_dir = af3_output_tree(tmp_path / "af3",
                                      predictions=[(4.0, 0.50)])
    empty = tmp_path / "empty"
    empty.mkdir()
    summary = tmp_path / "s.csv"
    assert _run(monkeypatch, truth, empty, summary) == 0
    assert _row(summary)[1]["af3_nomsa_failures"] == "no_cif_files_found"


@pytest.mark.logic
def test_a_prediction_with_no_confidence_json_still_contributes_its_ra_eff(
        monkeypatch, tmp_path):
    """ra_eff comes from the structure and iPTM from a sibling JSON, so the
    two can fail independently. Dropping the prediction entirely because its
    confidence write did not land would throw away the pose evidence, which
    is the number the campaign actually gates on."""
    truth, out_dir = af3_output_tree(
        tmp_path / "af3", predictions=[(4.0, 0.80), (7.0, None)],
        truth_separation=4.0)
    summary = tmp_path / "s.csv"
    _run(monkeypatch, truth, out_dir, summary)
    _fields, row = _row(summary)
    assert row["af3_nomsa_total_predictions"] == "2"
    assert row["af3_nomsa_mean_ra_eff"] == "1.500"   # (0.0 + 3.0) / 2
    # Only one prediction supplied an iPTM, so both reductions are that one.
    assert row["af3_nomsa_best_iptm"] == "0.8000"
    assert row["af3_nomsa_mean_iptm"] == "0.8000"


# ── The chain flags ──────────────────────────────────────────────────

@pytest.mark.logic
def test_the_chain_flags_reach_the_pose_computation(monkeypatch, tmp_path):
    """Naming a chain the structures do not hold must make every prediction
    fail rather than quietly scoring against whichever chain came first. The
    tag counts the failures, and the metric columns stay blank."""
    truth, out_dir = af3_output_tree(
        tmp_path / "af3", predictions=[(4.0, 0.80), (7.0, 0.60)])
    summary = tmp_path / "s.csv"
    assert _run(monkeypatch, truth, out_dir, summary, effector_chain="Z") == 0
    _fields, row = _row(summary)
    assert row["af3_nomsa_total_predictions"] == "0"
    assert row["af3_nomsa_best_ra_eff"] == ""
    assert "per_prediction_rmsd_errors:2" in row["af3_nomsa_failures"]
    assert "no_successful_rmsd" in row["af3_nomsa_failures"]


@pytest.mark.logic
def test_the_default_chains_are_the_pipeline_s_a_and_b(monkeypatch, tmp_path):
    """`af3_output_tree` writes the receptor as A and the effector as B, which
    is what the pipeline's rfdiff_output_*_chain params hold, so a run with
    neither flag must resolve the pose rather than failing."""
    truth, out_dir = af3_output_tree(tmp_path / "af3",
                                     predictions=[(4.0, 0.80)])
    summary = tmp_path / "s.csv"
    _run(monkeypatch, truth, out_dir, summary)
    assert _row(summary)[1]["af3_nomsa_total_predictions"] == "1"


# ── The sibling import ───────────────────────────────────────────────

@pytest.mark.logic
def test_the_module_puts_its_own_directory_on_the_import_path():
    """A Nextflow stage invokes this script by absolute path with no
    PYTHONPATH set, so `from af3_confidence import AF3ConfidenceAggregate`
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
            "parse_af3_output__pathless", BIN / "parse_af3_output.py")
        module = importlib.util.module_from_spec(spec)
        # @dataclass resolves its own module out of sys.modules, so the
        # loading module has to be registered before it is executed.
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
        assert str(BIN) in sys.path, "the module must add its own directory"
    finally:
        sys.path[:] = saved
        sys.modules.pop("parse_af3_output__pathless", None)
    assert module.AF3ConfidenceAggregate is pao.AF3ConfidenceAggregate


# ── The output file ──────────────────────────────────────────────────

@pytest.mark.logic
def test_the_output_directory_is_created_rather_than_assumed(monkeypatch, tmp_path):
    """The flag is passed an absolute path in some callers, so the parent may
    not exist when the stage starts."""
    truth, out_dir = af3_output_tree(tmp_path / "af3",
                                     predictions=[(4.0, 0.80)])
    summary = tmp_path / "a" / "b" / "af3_nomsa_summary.csv"
    assert _run(monkeypatch, truth, out_dir, summary) == 0
    assert summary.is_file()


@pytest.mark.logic
def test_the_run_reports_its_headline_numbers_on_stdout(monkeypatch, tmp_path,
                                                        capsys):
    """The Nextflow log line is how a failed AF3 job is spotted without
    opening the CSV."""
    truth, out_dir = af3_output_tree(
        tmp_path / "af3", predictions=[(4.0, 0.80), (7.0, 0.60)],
        truth_separation=4.0)
    summary = tmp_path / "s.csv"
    _run(monkeypatch, truth, out_dir, summary, seq_name="design_3_seq_2",
         ra_eff_threshold=2.0)
    printed = capsys.readouterr().out
    assert "[design_3_seq_2] best ra_eff = 0.000" in printed
    assert "n_correct = 1/2" in printed
