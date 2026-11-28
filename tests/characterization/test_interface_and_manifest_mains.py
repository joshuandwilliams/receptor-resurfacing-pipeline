"""Drive the DockQ enrichment stage and the survivor manifest stage end to end.

`tests/characterization/test_interface_and_manifest.py` covers the helpers in
these two modules. What it cannot reach is the part a Nextflow process actually
invokes, which is `main()`: argument handling, the per-row loop, the output
header, and the row a failure produces.

Both stages are built so one design failing does not take the stage down.
`compute_interface_metrics` writes a blank metric plus a token in
`dockq_failures`, and `extract_survivor_manifest` drops the row and counts the
reason in a histogram. Those two behaviours are the ones worth pinning. A stage
that raised instead would lose a whole cohort to one bad PDB, and a stage that
dropped silently would hand the next stage a shorter channel than it expects.

DockQ is not installed here. Where a test needs the call to succeed it installs
a stub `DockQ.DockQ` into `sys.modules`, so the chain-mapping choice and the
formatting of the four output columns are exercised against values known before
the code runs.
"""
from __future__ import annotations

import csv
import json
import sys
import types
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
BIN = REPO / "bin"
HERE = Path(__file__).resolve().parent
for _p in (str(BIN), str(HERE)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import compute_interface_metrics as cim  # noqa: E402
import extract_survivor_manifest as esm  # noqa: E402
from fixtures_runtree import complex_pdb  # noqa: E402

gemmi = pytest.importorskip("gemmi")


# ── Shared CSV helpers ───────────────────────────────────────────────

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


# ── Stubbing DockQ ───────────────────────────────────────────────────

# Deliberately more digits than the output columns carry, so the per-column
# precision below is a truncation the test can see rather than a copy of the
# input.  Each rounds unambiguously: 1.23449 -> 1.234, 0.567849 -> 0.5678,
# 0.765451 -> 0.7655.
_GOOD_RESULT = {"iRMSD": 1.23449, "fnat": 0.567849, "DockQ": 0.765451}


def _install_fake_dockq(monkeypatch, *, results=None, load_raises=None,
                        run_raises=None):
    """Put a stub `DockQ.DockQ` on sys.modules.

    `_run_dockq` imports it inside the function, so the stub only has to be
    in place for the duration of the call. Returning the results dict
    directly is what makes the expected irmsd/fnat/dockq known in advance.
    """
    pkg = types.ModuleType("DockQ")
    mod = types.ModuleType("DockQ.DockQ")

    def load_pdb(path):
        if load_raises is not None:
            raise load_raises
        return f"structure:{path}"

    def run_on_all_native_interfaces(model, native):
        if run_raises is not None:
            raise run_raises
        return (results if results is not None else {}), 0.0

    mod.load_PDB = load_pdb
    mod.run_on_all_native_interfaces = run_on_all_native_interfaces
    pkg.DockQ = mod
    monkeypatch.setitem(sys.modules, "DockQ", pkg)
    monkeypatch.setitem(sys.modules, "DockQ.DockQ", mod)


# ── Sidecar discovery ────────────────────────────────────────────────

@pytest.mark.logic
def test_a_sidecar_beside_the_pdb_is_found(tmp_path):
    pdb = tmp_path / "pred" / "model_0.pdb"
    pdb.parent.mkdir(parents=True)
    pdb.write_text("ATOM\n")
    want = pdb.parent / "confidence_model_0.json"
    want.write_text("{}")
    assert cim._find_sidecar(pdb, "confidence", "json") == want


@pytest.mark.logic
def test_a_sidecar_in_its_own_prefix_directory_is_found(tmp_path):
    """Boltz writes the PAE arrays into a sibling directory named after the
    prefix. Looking only beside the PDB would miss every one of them."""
    pdb = tmp_path / "pred" / "model_0.pdb"
    pdb.parent.mkdir(parents=True)
    pdb.write_text("ATOM\n")
    want = tmp_path / "pae" / "pae_model_0.npz"
    want.parent.mkdir()
    want.write_text("x")
    assert cim._find_sidecar(pdb, "pae", "npz") == want


@pytest.mark.logic
def test_a_sidecar_further_up_the_tree_is_found_by_recursion(tmp_path):
    """The last resort walks down from two levels above the PDB's directory,
    which is where a Boltz prediction root sits."""
    pdb = tmp_path / "pred" / "records" / "model_0.pdb"
    pdb.parent.mkdir(parents=True)
    pdb.write_text("ATOM\n")
    want = tmp_path / "elsewhere" / "deeper" / "confidence_model_0.json"
    want.parent.mkdir(parents=True)
    want.write_text("{}")
    assert cim._find_sidecar(pdb, "confidence", "json") == want


@pytest.mark.logic
def test_a_sidecar_that_is_not_anywhere_is_none(tmp_path):
    pdb = tmp_path / "a" / "b" / "model_0.pdb"
    pdb.parent.mkdir(parents=True)
    pdb.write_text("ATOM\n")
    assert cim._find_sidecar(pdb, "confidence", "json") is None


# ── Ground truth from plan.json ──────────────────────────────────────

def _workdir_with_plan(root: Path, seq_name: str, plan: dict,
                       *, under_cycle_0: bool = False) -> Path:
    wd = root / "runs" / seq_name
    target = wd / "cycle_0" / "plan.json" if under_cycle_0 else wd / "plan.json"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(plan))
    return wd


@pytest.mark.logic
@pytest.mark.parametrize("under_cycle_0", [False, True])
def test_the_ground_truth_is_read_from_either_plan_location(tmp_path, under_cycle_0):
    gt = tmp_path / "truth.pdb"
    gt.write_text("ATOM\n")
    _workdir_with_plan(tmp_path, "design_1_seq_1", {"ground_truth": str(gt)},
                       under_cycle_0=under_cycle_0)
    got = cim._find_ground_truth("design_1_seq_1", str(tmp_path / "runs" / "*"))
    assert got == gt


@pytest.mark.logic
def test_only_the_workdir_whose_name_matches_is_read(tmp_path):
    """The glob matches every sequence's workdir. Taking the first would hand
    one design another design's ground truth, and every metric after that
    would be a comparison against the wrong structure."""
    mine = tmp_path / "mine.pdb"
    mine.write_text("ATOM\n")
    _workdir_with_plan(tmp_path, "design_1_seq_1", {"ground_truth": "/other.pdb"})
    _workdir_with_plan(tmp_path, "design_2_seq_1", {"ground_truth": str(mine)})
    got = cim._find_ground_truth("design_2_seq_1", str(tmp_path / "runs" / "*"))
    assert got == mine


@pytest.mark.logic
def test_a_plan_without_a_ground_truth_key_is_none(tmp_path):
    _workdir_with_plan(tmp_path, "design_1_seq_1", {"something_else": 1})
    assert cim._find_ground_truth(
        "design_1_seq_1", str(tmp_path / "runs" / "*")) is None


@pytest.mark.logic
def test_an_unparsable_plan_is_none_rather_than_an_exception(tmp_path):
    """A truncated plan.json is what a killed engine run leaves behind. The
    row has to come back flagged, not take the stage down."""
    wd = tmp_path / "runs" / "design_1_seq_1"
    wd.mkdir(parents=True)
    (wd / "plan.json").write_text("{not json")
    assert cim._find_ground_truth(
        "design_1_seq_1", str(tmp_path / "runs" / "*")) is None


@pytest.mark.logic
def test_no_matching_workdir_gives_no_ground_truth(tmp_path):
    (tmp_path / "runs").mkdir()
    assert cim._find_ground_truth(
        "design_1_seq_1", str(tmp_path / "runs" / "*")) is None


# ── The DockQ call ───────────────────────────────────────────────────

@pytest.mark.logic
def test_dockq_not_being_installed_is_reported_not_raised(tmp_path):
    """DockQ is an optional dependency. Without it every row must come back
    flagged rather than the stage dying on import."""
    irmsd, fnat, dockq, err = cim._run_dockq(
        tmp_path / "m.pdb", tmp_path / "n.pdb", "A", "B")
    assert (irmsd, fnat, dockq) == (None, None, None)
    assert err is not None and err.startswith("dockq_import_failed")


@pytest.mark.logic
def test_dockq_values_come_from_the_receptor_effector_interface(monkeypatch, tmp_path):
    """A multi-interface result is searched for the (receptor, effector) pair
    by chain set, so the mapping DockQ happens to return first does not
    decide which interface is reported."""
    _install_fake_dockq(monkeypatch, results={
        ("C", "D"): {"iRMSD": 99.0, "fnat": 0.0, "DockQ": 0.0},
        ("B", "A"): _GOOD_RESULT,
    })
    irmsd, fnat, dockq, err = cim._run_dockq(
        tmp_path / "m.pdb", tmp_path / "n.pdb", "A", "B")
    assert err is None
    assert irmsd == pytest.approx(1.23449)
    assert fnat == pytest.approx(0.567849)
    assert dockq == pytest.approx(0.765451)


@pytest.mark.logic
def test_a_renamed_mapping_is_still_used_when_no_ab_pair_exists(monkeypatch, tmp_path):
    """DockQ may rename the chains while searching for the best mapping. The
    only interface it found is the one to report."""
    _install_fake_dockq(monkeypatch, results={("C", "D"): _GOOD_RESULT})
    irmsd, _fnat, _dockq, err = cim._run_dockq(
        tmp_path / "m.pdb", tmp_path / "n.pdb", "A", "B")
    assert err is None
    assert irmsd == pytest.approx(1.23449)


@pytest.mark.logic
def test_no_interface_at_all_is_reported(monkeypatch, tmp_path):
    _install_fake_dockq(monkeypatch, results={})
    irmsd, fnat, dockq, err = cim._run_dockq(
        tmp_path / "m.pdb", tmp_path / "n.pdb", "A", "B")
    assert (irmsd, fnat, dockq) == (None, None, None)
    assert err == "dockq_no_interface"


@pytest.mark.logic
def test_a_structure_dockq_cannot_load_is_reported_by_exception_type(monkeypatch, tmp_path):
    _install_fake_dockq(monkeypatch, load_raises=ValueError("bad atom"))
    _irmsd, _fnat, _dockq, err = cim._run_dockq(
        tmp_path / "m.pdb", tmp_path / "n.pdb", "A", "B")
    assert err == "dockq_load_failed:ValueError"


@pytest.mark.logic
def test_a_dockq_run_that_blows_up_is_reported_by_exception_type(monkeypatch, tmp_path):
    _install_fake_dockq(monkeypatch, run_raises=RuntimeError("no common chains"))
    _irmsd, _fnat, _dockq, err = cim._run_dockq(
        tmp_path / "m.pdb", tmp_path / "n.pdb", "A", "B")
    assert err == "dockq_run_failed:RuntimeError"


# ── One row ──────────────────────────────────────────────────────────

def _row(tmp_path: Path, **overrides) -> dict:
    model = tmp_path / "model.pdb"
    complex_pdb(model)
    row = {"mpnn_sequence": "design_1_seq_1", "rep_canonical_pdb": str(model)}
    row.update(overrides)
    return row


@pytest.mark.logic
def test_a_row_with_no_canonical_pdb_is_flagged_and_keeps_its_other_columns(tmp_path):
    out = cim.process_row(
        {"mpnn_sequence": "design_1_seq_1", "rep_canonical_pdb": "",
         "cross_tier": "A"},
        str(tmp_path / "runs" / "*"), "A", "B")
    assert out["dockq_failures"] == "no_canonical_pdb"
    assert out["irmsd"] == out["fnat"] == out["dockq"] == ""
    assert out["cross_tier"] == "A", "the input columns must survive the failure"


@pytest.mark.logic
def test_a_canonical_pdb_that_is_not_on_disk_is_named_in_the_flag(tmp_path):
    """The token carries the filename, so a stale absolute path is
    diagnosable from the output CSV without opening a log."""
    out = cim.process_row(
        {"mpnn_sequence": "d", "rep_canonical_pdb": "/gone/model_7.pdb"},
        str(tmp_path / "runs" / "*"), "A", "B")
    assert out["dockq_failures"] == "canonical_pdb_missing:model_7.pdb"


@pytest.mark.logic
def test_a_row_with_no_ground_truth_names_the_sequence(tmp_path):
    out = cim.process_row(_row(tmp_path), str(tmp_path / "runs" / "*"), "A", "B")
    assert out["dockq_failures"] == "ground_truth_not_found_for:design_1_seq_1"


@pytest.mark.logic
def test_a_successful_row_formats_each_metric_to_its_own_precision(monkeypatch, tmp_path):
    """irmsd is an angstrom distance written to 3 dp, fnat and dockq are
    0-1 scores written to 4 dp. The precisions are part of the committed
    output, so a change to either shows up here."""
    gt = tmp_path / "truth.pdb"
    complex_pdb(gt)
    _workdir_with_plan(tmp_path, "design_1_seq_1", {"ground_truth": str(gt)})
    _install_fake_dockq(monkeypatch, results={("A", "B"): _GOOD_RESULT})

    out = cim.process_row(_row(tmp_path), str(tmp_path / "runs" / "*"), "A", "B")
    assert out["dockq_failures"] == ""
    assert out["irmsd"] == "1.234"
    assert out["fnat"] == "0.5678"
    assert out["dockq"] == "0.7655"


@pytest.mark.logic
def test_a_dockq_failure_blanks_the_metrics_and_records_the_reason(monkeypatch, tmp_path):
    gt = tmp_path / "truth.pdb"
    complex_pdb(gt)
    _workdir_with_plan(tmp_path, "design_1_seq_1", {"ground_truth": str(gt)})
    _install_fake_dockq(monkeypatch, run_raises=RuntimeError("boom"))

    out = cim.process_row(_row(tmp_path), str(tmp_path / "runs" / "*"), "A", "B")
    assert out["dockq_failures"] == "dockq_run_failed:RuntimeError"
    assert out["irmsd"] == out["fnat"] == out["dockq"] == ""


@pytest.mark.logic
@pytest.mark.xfail(strict=False, reason=(
    "Known defect, reported rather than fixed: DockQ returning nan for a "
    "metric leaves the column blank and dockq_failures empty, so main() "
    "counts the row as fully populated.  Emitting a token would change the "
    "dockq_failures column of committed "
    "cross_sequence_summary_with_interface_metrics.csv output."))
def test_a_nan_metric_must_not_be_reported_as_a_success(monkeypatch, tmp_path):
    """A blank irmsd with no failure token is indistinguishable from a
    metric that was never requested.

    `_run_dockq` substitutes nan for any key DockQ omits, and `process_row`
    turns nan into an empty string without adding a reason. The row then
    reads as a clean success in the stage's "N rows fully populated" count
    while carrying no measurement at all.
    """
    gt = tmp_path / "truth.pdb"
    complex_pdb(gt)
    _workdir_with_plan(tmp_path, "design_1_seq_1", {"ground_truth": str(gt)})
    _install_fake_dockq(monkeypatch, results={("A", "B"): {}})

    out = cim.process_row(_row(tmp_path), str(tmp_path / "runs" / "*"), "A", "B")
    assert out["irmsd"] == ""
    assert out["dockq_failures"] != ""


# ── The weighted pair map ────────────────────────────────────────────

@pytest.mark.logic
def test_the_pair_map_keys_on_positional_index_and_carries_the_gaussian(tmp_path):
    """Two residues per chain, 4 A apart in z, so the facing pairs sit on the
    peak and weigh exactly 1.0 while the diagonal pairs are further off."""
    p = complex_pdb(tmp_path / "c.pdb", n_receptor=2, n_effector=2,
                    spacing=3.8, separation=4.0)
    got = cim._build_weighted_pair_map(p, "A", "B", 8.0)
    assert set(got) == {(0, 0), (0, 1), (1, 0), (1, 1)}
    assert got[(0, 0)] == pytest.approx(1.0)
    assert got[(1, 1)] == pytest.approx(1.0)
    assert got[(0, 1)] < 1.0 and got[(1, 0)] < 1.0


@pytest.mark.logic
def test_a_pair_exactly_on_the_cutoff_is_kept(tmp_path):
    """The comparison is `distance <= cutoff`, not strict, so the documented
    cutoff is an attainable value rather than one that is always excluded."""
    p = complex_pdb(tmp_path / "c.pdb", n_receptor=1, n_effector=1,
                    separation=6.0)
    on_it = cim._build_weighted_pair_map(p, "A", "B", 6.0)
    assert set(on_it) == {(0, 0)}
    assert on_it[(0, 0)] == pytest.approx(cim._gaussian_contact_weight(6.0))
    assert cim._build_weighted_pair_map(p, "A", "B", 5.999) == {}


@pytest.mark.logic
def test_a_residue_with_no_anchor_atom_still_advances_the_index(tmp_path):
    """A residue carrying neither CB nor CA is dropped from the list but not
    from the numbering, so the positional keys stay aligned with the other
    structure's. Skipping the increment instead would shift every key after
    it and compare two different residues."""
    lines = [
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C",
        "ATOM      2  SG  CYS A   2       3.800   0.000   0.000  1.00  0.00           S",
        "ATOM      3  CA  ALA A   3       7.600   0.000   0.000  1.00  0.00           C",
    ]
    p = tmp_path / "gap.pdb"
    p.write_text("\n".join(lines) + "\nEND\n")
    got = cim._collect_chain_cb_positions(gemmi.read_structure(str(p)), "A")
    assert [g[0] for g in got] == [0, 2]


@pytest.mark.logic
def test_an_absent_chain_gives_an_empty_pair_map_rather_than_an_error(tmp_path):
    """An empty chain is treated as no contacts, which is what the caller's
    union denominator needs."""
    p = complex_pdb(tmp_path / "c.pdb")
    assert cim._build_weighted_pair_map(p, "A", "Z", 8.0) == {}


@pytest.mark.logic
def test_a_structure_that_will_not_load_raises_a_tagged_error(tmp_path):
    """The token is what the caller turns into a row flag, so it has to say
    which step failed rather than surfacing gemmi's own message."""
    with pytest.raises(RuntimeError) as exc:
        cim._build_weighted_pair_map(tmp_path / "absent.pdb", "A", "B", 8.0)
    assert str(exc.value).startswith("structure_load_failed:")


# ── compute_interface_metrics main() ─────────────────────────────────

def _run_cim(monkeypatch, tmp_path, in_csv: Path, out_csv: Path) -> int:
    monkeypatch.setattr(sys, "argv", [
        "compute_interface_metrics.py",
        "--input-csv", str(in_csv),
        "--output-csv", str(out_csv),
        "--workdirs-glob", str(tmp_path / "runs" / "*"),
        "--receptor-chain", "A", "--effector-chain", "B",
    ])
    return cim.main()


@pytest.mark.logic
def test_the_dockq_columns_are_appended_and_one_full_row_is_known(monkeypatch, tmp_path):
    """Every value in this row is fixed before the stage runs. The three
    input columns come from the input CSV and the four new ones from the
    stub's result dict."""
    model = tmp_path / "model.pdb"
    complex_pdb(model)
    gt = tmp_path / "truth.pdb"
    complex_pdb(gt)
    _workdir_with_plan(tmp_path, "design_1_seq_1", {"ground_truth": str(gt)})
    _install_fake_dockq(monkeypatch, results={("A", "B"): _GOOD_RESULT})

    in_csv = _write_csv(
        tmp_path / "in.csv", ["mpnn_sequence", "cross_tier", "rep_canonical_pdb"],
        [{"mpnn_sequence": "design_1_seq_1", "cross_tier": "A",
          "rep_canonical_pdb": str(model)}])
    out_csv = tmp_path / "out" / "extended.csv"
    assert _run_cim(monkeypatch, tmp_path, in_csv, out_csv) == 0

    header, rows = _read_csv(out_csv)
    assert header == ["mpnn_sequence", "cross_tier", "rep_canonical_pdb",
                      "irmsd", "fnat", "dockq", "dockq_failures"]
    assert rows == [{
        "mpnn_sequence": "design_1_seq_1",
        "cross_tier": "A",
        "rep_canonical_pdb": str(model),
        "irmsd": "1.234",
        "fnat": "0.5678",
        "dockq": "0.7655",
        "dockq_failures": "",
    }]


@pytest.mark.logic
def test_one_failing_row_does_not_cost_the_others(monkeypatch, tmp_path):
    """The failing row carries a token and the good rows keep their metrics.
    Raising instead would lose a whole cohort to one missing PDB."""
    model = tmp_path / "model.pdb"
    complex_pdb(model)
    gt = tmp_path / "truth.pdb"
    complex_pdb(gt)
    for name in ("design_1_seq_1", "design_3_seq_1"):
        _workdir_with_plan(tmp_path, name, {"ground_truth": str(gt)})
    _install_fake_dockq(monkeypatch, results={("A", "B"): _GOOD_RESULT})

    in_csv = _write_csv(
        tmp_path / "in.csv", ["mpnn_sequence", "rep_canonical_pdb"],
        [{"mpnn_sequence": "design_1_seq_1", "rep_canonical_pdb": str(model)},
         {"mpnn_sequence": "design_2_seq_1", "rep_canonical_pdb": "/gone/x.pdb"},
         {"mpnn_sequence": "design_3_seq_1", "rep_canonical_pdb": str(model)}])
    out_csv = tmp_path / "out.csv"
    assert _run_cim(monkeypatch, tmp_path, in_csv, out_csv) == 0

    _header, rows = _read_csv(out_csv)
    assert [r["mpnn_sequence"] for r in rows] == [
        "design_1_seq_1", "design_2_seq_1", "design_3_seq_1"]
    assert rows[1]["dockq_failures"] == "canonical_pdb_missing:x.pdb"
    assert rows[1]["dockq"] == ""
    assert rows[0]["dockq"] == rows[2]["dockq"] == "0.7655"


@pytest.mark.logic
def test_a_rerun_overwrites_the_columns_rather_than_duplicating_them(monkeypatch, tmp_path):
    """The stage is idempotent over its own output. Appending a second
    `dockq` column would make the CSV unreadable by DictReader."""
    in_csv = _write_csv(
        tmp_path / "in.csv",
        ["mpnn_sequence", "rep_canonical_pdb", "irmsd", "fnat", "dockq",
         "dockq_failures"],
        [{"mpnn_sequence": "design_1_seq_1", "rep_canonical_pdb": "",
          "irmsd": "9.9", "fnat": "0.1", "dockq": "0.2",
          "dockq_failures": "stale"}])
    out_csv = tmp_path / "out.csv"
    assert _run_cim(monkeypatch, tmp_path, in_csv, out_csv) == 0

    header, rows = _read_csv(out_csv)
    assert len(header) == 6 and header.count("dockq") == 1
    assert rows[0]["irmsd"] == ""
    assert rows[0]["dockq_failures"] == "no_canonical_pdb"


@pytest.mark.logic
def test_a_header_only_csv_gives_a_header_only_output(monkeypatch, tmp_path):
    """Zero rows is a legitimate cohort, not an error. The header still has
    to carry the four new columns so the next stage's reader works."""
    in_csv = _write_csv(tmp_path / "in.csv",
                        ["mpnn_sequence", "rep_canonical_pdb"], [])
    out_csv = tmp_path / "out.csv"
    assert _run_cim(monkeypatch, tmp_path, in_csv, out_csv) == 0

    header, rows = _read_csv(out_csv)
    assert rows == []
    assert header[-4:] == cim.NEW_COLUMNS


@pytest.mark.logic
def test_an_input_csv_with_no_header_at_all_is_an_error_exit(monkeypatch, tmp_path):
    """An empty file means the upstream stage wrote nothing. Exiting 1 stops
    the cascade rather than producing an empty extended CSV that reads as a
    cohort with no survivors."""
    in_csv = tmp_path / "empty.csv"
    in_csv.write_text("")
    out_csv = tmp_path / "out.csv"
    assert _run_cim(monkeypatch, tmp_path, in_csv, out_csv) == 1
    assert not out_csv.exists()


# ── Chain sequence extraction ────────────────────────────────────────

@pytest.mark.logic
def test_a_chain_sequence_is_read_at_its_own_length(tmp_path):
    """The manifest's receptor_seq is what AF3 is handed, so it has to be the
    designed construct's own chain rather than the whole file."""
    p = complex_pdb(tmp_path / "c.pdb", n_receptor=7, n_effector=4)
    assert esm._extract_chain_seq(p, "A") == "A" * 7
    assert esm._extract_chain_seq(p, "B") == "A" * 4


@pytest.mark.logic
def test_an_absent_chain_is_none_rather_than_an_exception(tmp_path):
    """The canonical helper raises on a missing chain. This wrapper converts
    that to None because every caller branches on the Optional."""
    p = complex_pdb(tmp_path / "c.pdb")
    assert esm._extract_chain_seq(p, "Z") is None


@pytest.mark.logic
def test_an_unreadable_pdb_is_none(tmp_path):
    assert esm._extract_chain_seq(tmp_path / "absent.pdb", "A") is None


# ── extract_survivor_manifest main() ─────────────────────────────────

MANIFEST_COLUMNS = ["seq_name", "canonical_pdb_abs", "ground_truth_abs",
                    "effector_template_cif_abs", "receptor_seq", "effector_seq"]


def _run_esm(monkeypatch, tmp_path, in_csv: Path, manifest: Path,
             *, tier_filter: str = "all") -> int:
    monkeypatch.setattr(sys, "argv", [
        "extract_survivor_manifest.py",
        "--input-csv", str(in_csv),
        "--workdirs-glob", str(tmp_path / "runs" / "*"),
        "--output-manifest", str(manifest),
        "--tier-filter", tier_filter,
    ])
    return esm.main()


def _survivor(tmp_path: Path, seq_name: str = "design_1_seq_1",
              *, plan_extra: dict | None = None) -> dict:
    """One row whose every manifest value is known before the stage runs."""
    model = tmp_path / "runs" / seq_name / "cycle_0" / "model.pdb"
    complex_pdb(model, n_receptor=7, n_effector=4)
    gt = tmp_path / "truth.pdb"
    complex_pdb(gt)
    cif = tmp_path / "effector.cif"
    cif.write_text("data_x\n")
    plan = {"ground_truth": str(gt), "effector_template_cif": str(cif)}
    plan.update(plan_extra or {})
    _workdir_with_plan(tmp_path, seq_name, plan)
    return {"mpnn_sequence": seq_name, "cross_tier": "A",
            "rep_canonical_pdb": str(model)}


@pytest.mark.logic
def test_one_survivor_produces_one_fully_known_manifest_row(monkeypatch, tmp_path):
    """The two sequences are the receptor and effector chains of the
    canonical PDB itself, 7 and 4 residues of alanine, so the row is
    arithmetic rather than something to read off a real run."""
    row = _survivor(tmp_path)
    in_csv = _write_csv(tmp_path / "in.csv", list(row), [row])
    manifest = tmp_path / "out" / "manifest.csv"
    assert _run_esm(monkeypatch, tmp_path, in_csv, manifest) == 0

    header, rows = _read_csv(manifest)
    assert header == MANIFEST_COLUMNS
    assert rows == [{
        "seq_name": "design_1_seq_1",
        "canonical_pdb_abs": str(Path(row["rep_canonical_pdb"]).resolve()),
        "ground_truth_abs": str((tmp_path / "truth.pdb").resolve()),
        "effector_template_cif_abs": str((tmp_path / "effector.cif").resolve()),
        "receptor_seq": "A" * 7,
        "effector_seq": "A" * 4,
    }]


@pytest.mark.logic
def test_zero_input_rows_gives_a_header_only_manifest_and_a_clean_exit(monkeypatch, tmp_path):
    """No rows at all is not the empty-manifest failure. There was nothing to
    skip, so there is nothing to diagnose."""
    in_csv = _write_csv(tmp_path / "in.csv",
                        ["mpnn_sequence", "rep_canonical_pdb"], [])
    manifest = tmp_path / "manifest.csv"
    assert _run_esm(monkeypatch, tmp_path, in_csv, manifest) == 0
    header, rows = _read_csv(manifest)
    assert header == MANIFEST_COLUMNS
    assert rows == []


@pytest.mark.logic
def test_every_row_being_skipped_is_an_error_exit_with_the_reason_named(
        monkeypatch, tmp_path, capsys):
    """Previously this exited 0 and the downstream fan-out ran over an empty
    channel, which Nextflow treats as a successful no-op. Exit 2 plus the
    histogram is what makes a whole-cohort path failure visible."""
    in_csv = _write_csv(
        tmp_path / "in.csv", ["mpnn_sequence", "rep_canonical_pdb"],
        [{"mpnn_sequence": "design_1_seq_1", "rep_canonical_pdb": "/gone.pdb"}])
    manifest = tmp_path / "manifest.csv"
    assert _run_esm(monkeypatch, tmp_path, in_csv, manifest) == 2
    out = capsys.readouterr()
    assert "workdir_missing: 1" in out.out
    assert "0 survivors" in out.err


@pytest.mark.logic
@pytest.mark.parametrize("mutate,reason", [
    (lambda r, t: r.update(mpnn_sequence=""), "no_mpnn_sequence"),
    (lambda r, t: r.update(rep_canonical_pdb=""), "canonical_pdb_missing"),
    (lambda r, t: r.update(rep_canonical_pdb="/elsewhere/model.pdb"),
     "canonical_pdb_missing"),
])
def test_a_row_the_stage_cannot_use_is_dropped_with_its_reason_counted(
        monkeypatch, tmp_path, capsys, mutate, reason):
    """Each reason is its own histogram key. One bucket called 'skipped'
    would make a path problem indistinguishable from a tier problem."""
    good = _survivor(tmp_path, "design_1_seq_1")
    bad = _survivor(tmp_path, "design_2_seq_1")
    mutate(bad, tmp_path)
    in_csv = _write_csv(tmp_path / "in.csv", list(good), [good, bad])
    manifest = tmp_path / "manifest.csv"
    assert _run_esm(monkeypatch, tmp_path, in_csv, manifest) == 0

    _header, rows = _read_csv(manifest)
    assert [r["seq_name"] for r in rows] == ["design_1_seq_1"]
    assert f"{reason}: 1" in capsys.readouterr().out


@pytest.mark.logic
def test_a_stale_canonical_pdb_is_remapped_onto_the_local_workdir(monkeypatch, tmp_path):
    """The CSV stamps absolute paths at aggregation time. A fixture built on
    one host and consumed on another has to find the file by its
    /runs/<seq_name>/ tail rather than be skipped."""
    row = _survivor(tmp_path)
    row["rep_canonical_pdb"] = (
        "/built/on/another/host/runs/design_1_seq_1/cycle_0/model.pdb")
    in_csv = _write_csv(tmp_path / "in.csv", list(row), [row])
    manifest = tmp_path / "manifest.csv"
    assert _run_esm(monkeypatch, tmp_path, in_csv, manifest) == 0

    _header, rows = _read_csv(manifest)
    assert rows[0]["canonical_pdb_abs"] == str(
        (tmp_path / "runs" / "design_1_seq_1" / "cycle_0" / "model.pdb").resolve())


@pytest.mark.logic
def test_a_stale_ground_truth_is_remapped_against_this_checkout(monkeypatch, tmp_path):
    """plan.json carries the repo path the engine saw. Anchoring on the first
    canonical top-level directory is what lets a committed fixture resolve on
    a different machine."""
    here = Path(__file__).resolve()
    row = _survivor(tmp_path, plan_extra={
        "ground_truth": f"/some/other/checkout/tests/{here.relative_to(REPO / 'tests')}"})
    in_csv = _write_csv(tmp_path / "in.csv", list(row), [row])
    manifest = tmp_path / "manifest.csv"
    assert _run_esm(monkeypatch, tmp_path, in_csv, manifest) == 0

    _header, rows = _read_csv(manifest)
    assert rows[0]["ground_truth_abs"] == str(here)


@pytest.mark.logic
def test_a_stale_effector_template_is_remapped_through_the_workdir(monkeypatch, tmp_path):
    """The template is written inside the workdir under cycle_0, so the
    /runs/<seq>/ remap is tried before the repo-wide one."""
    template = tmp_path / "runs" / "design_1_seq_1" / "cycle_0" / "effector.cif"
    row = _survivor(tmp_path, plan_extra={
        "effector_template_cif": "/old/host/runs/design_1_seq_1/cycle_0/effector.cif"})
    template.write_text("data_x\n")
    in_csv = _write_csv(tmp_path / "in.csv", list(row), [row])
    manifest = tmp_path / "manifest.csv"
    assert _run_esm(monkeypatch, tmp_path, in_csv, manifest) == 0

    _header, rows = _read_csv(manifest)
    assert rows[0]["effector_template_cif_abs"] == str(template.resolve())


@pytest.mark.logic
@pytest.mark.parametrize("plan_extra,reason", [
    ({"ground_truth": None}, "ground_truth_missing"),
    ({"ground_truth": "/nowhere/truth.pdb"}, "ground_truth_missing"),
    ({"effector_template_cif": None}, "effector_template_missing"),
    ({"effector_template_cif": "/nowhere/eff.cif"}, "effector_template_missing"),
])
def test_an_unresolvable_plan_path_drops_the_row_with_its_reason(
        monkeypatch, tmp_path, capsys, plan_extra, reason):
    row = _survivor(tmp_path, plan_extra=plan_extra)
    in_csv = _write_csv(tmp_path / "in.csv", list(row), [row])
    manifest = tmp_path / "manifest.csv"
    assert _run_esm(monkeypatch, tmp_path, in_csv, manifest) == 2
    assert f"{reason}: 1" in capsys.readouterr().out


@pytest.mark.logic
def test_a_workdir_with_no_plan_is_skipped_for_that_reason(monkeypatch, tmp_path, capsys):
    model = tmp_path / "runs" / "design_1_seq_1" / "model.pdb"
    complex_pdb(model)
    in_csv = _write_csv(
        tmp_path / "in.csv", ["mpnn_sequence", "rep_canonical_pdb"],
        [{"mpnn_sequence": "design_1_seq_1", "rep_canonical_pdb": str(model)}])
    manifest = tmp_path / "manifest.csv"
    assert _run_esm(monkeypatch, tmp_path, in_csv, manifest) == 2
    assert "plan_json_missing: 1" in capsys.readouterr().out


@pytest.mark.logic
def test_an_unparsable_plan_is_skipped_for_its_own_reason(monkeypatch, tmp_path, capsys):
    """Distinct from plan_json_missing. A truncated plan means the engine run
    died mid-write, which is a different thing to fix."""
    row = _survivor(tmp_path)
    (tmp_path / "runs" / "design_1_seq_1" / "plan.json").write_text("{truncated")
    in_csv = _write_csv(tmp_path / "in.csv", list(row), [row])
    manifest = tmp_path / "manifest.csv"
    assert _run_esm(monkeypatch, tmp_path, in_csv, manifest) == 2
    assert "plan_json_unreadable: 1" in capsys.readouterr().out


@pytest.mark.logic
def test_a_canonical_pdb_with_only_one_chain_is_skipped(monkeypatch, tmp_path, capsys):
    """Both sequences are required, because AF3 is handed a two-chain job.
    Writing a row with a blank effector_seq would fail three stages later."""
    row = _survivor(tmp_path)
    complex_pdb(Path(row["rep_canonical_pdb"]), n_receptor=7, n_effector=0)
    in_csv = _write_csv(tmp_path / "in.csv", list(row), [row])
    manifest = tmp_path / "manifest.csv"
    assert _run_esm(monkeypatch, tmp_path, in_csv, manifest) == 2
    assert "seq_extraction_failed: 1" in capsys.readouterr().out


@pytest.mark.logic
def test_the_default_tier_filter_keeps_a_tier_none_row(monkeypatch, tmp_path):
    """'all' is the default so a failed design still gets AF3, biophysical
    and Rosetta numbers. Those rows are the diagnostic for why it failed."""
    row = _survivor(tmp_path)
    row["cross_tier"] = ""
    in_csv = _write_csv(tmp_path / "in.csv", list(row), [row])
    manifest = tmp_path / "manifest.csv"
    assert _run_esm(monkeypatch, tmp_path, in_csv, manifest) == 0
    _header, rows = _read_csv(manifest)
    assert [r["seq_name"] for r in rows] == ["design_1_seq_1"]


@pytest.mark.logic
@pytest.mark.parametrize("tier,kept", [
    ("A", True), ("B", True), ("C", True), ("", False), ("None", False),
    ("D", False), (" A ", True),
])
def test_the_abc_tier_filter_admits_exactly_a_b_and_c(monkeypatch, tmp_path, tier, kept):
    """The membership test is on the stripped value against {A, B, C}, so a
    whitespace-padded tier is kept and anything else is excluded. Letting
    'None' through would spend GPU time on rows the operator asked to skip."""
    row = _survivor(tmp_path)
    row["cross_tier"] = tier
    in_csv = _write_csv(tmp_path / "in.csv", list(row), [row])
    manifest = tmp_path / "manifest.csv"
    rc = _run_esm(monkeypatch, tmp_path, in_csv, manifest, tier_filter="abc")
    _header, rows = _read_csv(manifest)
    assert bool(rows) is kept
    assert rc == (0 if kept else 2)
