"""Drive bin/mpnn_cluster_sequences.py against constructed metadata and a stub mmseqs.

The MPNN_CLUSTER stage produces one number per identity threshold, and that
number is the pipeline's only measure of whether a design round produced 128
distinct sequences or 128 copies of three. MMseqs2 does the alignment; what
this module owns is three things, and each is asserted on input where the
answer is obvious:

  - which rows reach the FASTA at all, since a row with no design region has
    nothing to cluster and must not be counted as a sequence;
  - the cluster count, which is the number of DISTINCT representatives in
    MMseqs2's two-column TSV and not the number of lines in it;
  - what happens when MMseqs2 does not answer, because the fallback writes a
    number that the diversity plot then presents as a measurement.

MMseqs2 is not installed here and is not meant to be. The module reaches it
through `--mmseqs-bin`, so `fixtures_mpnn_stage.fake_mmseqs` puts an
executable on disk and the module's own `subprocess.run` finds it. Nothing in
the module is patched.

Column names come from bin/pipeline_correct_sequences.py and from
tests/proteinmpnn/example_output_files/sequences/sequence_metadata.csv; the
output header from tests/proteinmpnn/example_output_files/sequences/
mpnn_cluster_counts.csv.
"""
from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
BIN = REPO / "bin"
HERE = Path(__file__).resolve().parent
for _p in (str(BIN), str(HERE)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import mpnn_cluster_sequences as mc  # noqa: E402
from fixtures_mpnn_stage import (  # noqa: E402
    fake_mmseqs,
    mmseqs_calls,
    qc_row,
    sequence_metadata_csv,
)


def _read_counts(path: Path) -> list[dict]:
    with open(path) as handle:
        return list(csv.DictReader(handle))


def _run_main(monkeypatch, *args):
    monkeypatch.setattr(sys, "argv", ["mpnn_cluster_sequences.py", *args])
    mc.main()


# ── Which rows reach the FASTA ───────────────────────────────────────

@pytest.mark.logic
def test_the_fasta_carries_one_record_per_row_named_design_and_seq(tmp_path):
    """The record name is the only link back from a cluster representative to
    the sequence it came from, so it has to be design_D_seq_S exactly."""
    rows = [qc_row(design="0", seq="0", designed_residues="AVSVTEKGGWSDIT"),
            qc_row(design="28", seq="3", designed_residues="KKDDEEWWFFGGHH")]
    fasta = tmp_path / "design_regions.fasta"
    assert mc.write_design_fasta(rows, fasta) == 2
    assert fasta.read_text() == (
        ">design_0_seq_0\nAVSVTEKGGWSDIT\n"
        ">design_28_seq_3\nKKDDEEWWFFGGHH\n")


@pytest.mark.logic
def test_the_region_separator_is_stripped_before_clustering(tmp_path):
    """designed_residues is pipe-joined across design regions. MMseqs2 would
    treat the pipe as an unknown residue and let it count towards both the
    identity and the coverage, so a two-region design would cluster slightly
    differently from a one-region design of the same length."""
    rows = [qc_row(design="0", seq="0", designed_residues="AVSV|GDIV")]
    fasta = tmp_path / "f.fasta"
    assert mc.write_design_fasta(rows, fasta) == 1
    assert fasta.read_text() == ">design_0_seq_0\nAVSVGDIV\n"


@pytest.mark.logic
def test_a_row_with_no_design_region_is_skipped_and_not_counted(tmp_path):
    """A blank designed_residues is what a design with no de novo region
    leaves behind. Writing it as an empty FASTA record would make MMseqs2
    fail on the whole file, and counting it would inflate n_sequences."""
    rows = [qc_row(design="0", seq="0", designed_residues=""),
            qc_row(design="1", seq="0", designed_residues="AVSV")]
    fasta = tmp_path / "f.fasta"
    assert mc.write_design_fasta(rows, fasta) == 1
    assert fasta.read_text() == ">design_1_seq_0\nAVSV\n"


@pytest.mark.logic
def test_a_row_that_is_only_separators_is_also_skipped(tmp_path):
    """Two empty regions give "|", which is truthy before the strip and empty
    after it. The second emptiness check is what catches it."""
    rows = [qc_row(design="0", seq="0", designed_residues="|")]
    fasta = tmp_path / "f.fasta"
    assert mc.write_design_fasta(rows, fasta) == 0
    assert fasta.read_text() == ""


@pytest.mark.logic
def test_zero_rows_write_an_empty_fasta(tmp_path):
    fasta = tmp_path / "f.fasta"
    assert mc.write_design_fasta([], fasta) == 0
    assert fasta.read_text() == ""


# ── Counting clusters ────────────────────────────────────────────────

@pytest.mark.logic
def test_the_count_is_distinct_representatives_not_tsv_lines(tmp_path):
    """MMseqs2's _cluster.tsv is one line per sequence, representative in
    column 1. Four sequences in two clusters give four lines, and the answer
    is 2. Counting lines would report every design round as maximally diverse.
    """
    stub = fake_mmseqs(tmp_path / "stub", {"0.50": 2})
    fasta = tmp_path / "f.fasta"
    mc.write_design_fasta(
        [qc_row(design=str(i), designed_residues="AVSV") for i in range(4)],
        fasta)
    assert mc.run_mmseqs_cluster(fasta, 0.50, str(stub),
                                 str(tmp_path / "work")) == 2


@pytest.mark.logic
def test_every_sequence_its_own_cluster_is_the_maximum(tmp_path):
    stub = fake_mmseqs(tmp_path / "stub", {"0.90": 3})
    fasta = tmp_path / "f.fasta"
    mc.write_design_fasta(
        [qc_row(design=str(i), designed_residues="AVSV") for i in range(3)],
        fasta)
    assert mc.run_mmseqs_cluster(fasta, 0.90, str(stub),
                                 str(tmp_path / "work")) == 3


@pytest.mark.logic
def test_the_documented_mmseqs_flags_reach_the_binary(tmp_path):
    """The coverage mode and cluster mode are the difference between
    clustering the design regions and clustering whichever pairs happened to
    align. `--min-seq-id` is formatted to two decimals, which is also the
    format the plan file and the output CSV are keyed on.
    """
    root = tmp_path / "stub"
    stub = fake_mmseqs(root, {"0.30": 1})
    fasta = tmp_path / "f.fasta"
    mc.write_design_fasta(
        [qc_row(design=str(i), designed_residues="AVSV") for i in range(2)],
        fasta)
    mc.run_mmseqs_cluster(fasta, 0.30, str(stub), str(tmp_path / "work"))
    argv = mmseqs_calls(root)[0]["argv"]
    assert argv[0] == "easy-cluster"
    assert argv[argv.index("--min-seq-id") + 1] == "0.30"
    assert argv[argv.index("-c") + 1] == "0.8"
    assert argv[argv.index("--cov-mode") + 1] == "1"
    assert argv[argv.index("--cluster-mode") + 1] == "0"
    assert mmseqs_calls(root)[0]["names"] == ["design_0_seq_0", "design_1_seq_0"]


# ── When MMseqs2 does not answer ─────────────────────────────────────

@pytest.mark.logic
def test_a_non_zero_exit_returns_none_and_names_the_threshold(
        tmp_path, capsys):
    """None is the module's "no answer", distinct from a count of 0."""
    stub = fake_mmseqs(tmp_path / "stub", {"1.00": "fail"})
    fasta = tmp_path / "f.fasta"
    mc.write_design_fasta([qc_row(designed_residues="AVSV")], fasta)
    assert mc.run_mmseqs_cluster(fasta, 1.00, str(stub),
                                 str(tmp_path / "work")) is None
    assert "mmseqs failed at threshold 1.0" in capsys.readouterr().err


@pytest.mark.logic
def test_a_clean_exit_with_no_tsv_also_returns_none(tmp_path):
    """MMseqs2 can exit 0 having written nothing, e.g. when every sequence
    fails the coverage filter. The TSV existence check is the second guard."""
    stub = fake_mmseqs(tmp_path / "stub", {"0.70": "no_tsv"})
    fasta = tmp_path / "f.fasta"
    mc.write_design_fasta([qc_row(designed_residues="AVSV")], fasta)
    assert mc.run_mmseqs_cluster(fasta, 0.70, str(stub),
                                 str(tmp_path / "work")) is None


@pytest.mark.logic
def test_a_stale_tsv_from_an_earlier_run_is_removed_first(tmp_path):
    """The three result files are deleted before mmseqs is invoked. Without
    that, a failed run would be read as whatever the previous threshold left
    behind — here a count of 9 instead of no answer at all, which is the
    worst kind of wrong because it looks like a measurement.
    """
    work = tmp_path / "work"
    work.mkdir()
    stale = work / "clust_cluster.tsv"
    stale.write_text("".join(f"rep_{i}\tmember_{i}\n" for i in range(9)))
    stub = fake_mmseqs(tmp_path / "stub", {"0.60": "fail"})
    fasta = tmp_path / "f.fasta"
    mc.write_design_fasta([qc_row(designed_residues="AVSV")], fasta)
    assert mc.run_mmseqs_cluster(fasta, 0.60, str(stub), str(work)) is None
    assert not stale.exists()


@pytest.mark.logic
@pytest.mark.xfail(strict=False, reason=(
    "Defect: only CalledProcessError is caught, so an mmseqs binary that is "
    "not on PATH raises FileNotFoundError out of main() and fails the "
    "MPNN_CLUSTER stage. The module's own design is to warn and fall back to "
    "n_seqs, and a missing binary is the likelier of the two failures in a "
    "container that was built without MMseqs2. Not fixed here because the "
    "fix changes what the stage writes."))
def test_an_absent_mmseqs_binary_should_fall_back_rather_than_raise(tmp_path):
    fasta = tmp_path / "f.fasta"
    mc.write_design_fasta([qc_row(designed_residues="AVSV")], fasta)
    assert mc.run_mmseqs_cluster(fasta, 0.50, str(tmp_path / "not_installed"),
                                 str(tmp_path / "work")) is None


@pytest.mark.logic
def test_an_absent_mmseqs_binary_raises_today(tmp_path):
    """Pins the behaviour the xfail above describes, so the defect cannot be
    half-fixed without one of the two tests changing."""
    fasta = tmp_path / "f.fasta"
    mc.write_design_fasta([qc_row(designed_residues="AVSV")], fasta)
    with pytest.raises(FileNotFoundError):
        mc.run_mmseqs_cluster(fasta, 0.50, str(tmp_path / "not_installed"),
                              str(tmp_path / "work"))


# ── main() ───────────────────────────────────────────────────────────

@pytest.mark.logic
def test_main_writes_one_row_per_threshold_in_ascending_order(
        tmp_path, monkeypatch):
    """The thresholds are sorted before use, so the diversity plot's x-axis is
    monotonic whatever order the flag listed them in. 1.00,0.30,0.60 in must
    come out 0.30,0.60,1.00 — and the counts must travel with their own
    threshold rather than with their input position.
    """
    stub = fake_mmseqs(tmp_path / "stub", {"0.30": 1, "0.60": 2, "1.00": 4})
    metadata = sequence_metadata_csv(tmp_path / "meta.csv", [
        qc_row(design=str(i), designed_residues="AVSV" * (i + 1))
        for i in range(4)])
    out = tmp_path / "counts.csv"
    _run_main(monkeypatch, "--metadata-csv", str(metadata),
              "--output-csv", str(out), "--thresholds", "1.00,0.30,0.60",
              "--mmseqs-bin", str(stub))
    assert _read_counts(out) == [
        {"threshold": "0.30", "n_clusters": "1", "n_sequences": "4"},
        {"threshold": "0.60", "n_clusters": "2", "n_sequences": "4"},
        {"threshold": "1.00", "n_clusters": "4", "n_sequences": "4"},
    ]


@pytest.mark.logic
def test_main_falls_back_to_all_unique_when_a_threshold_fails(
        tmp_path, monkeypatch):
    """A failed threshold is written as n_clusters == n_sequences, which the
    diversity plot draws as "maximally diverse". That is the optimistic
    direction, so it is worth having pinned: 0.90 failed and reports 3, while
    0.50 really did cluster to 1.
    """
    stub = fake_mmseqs(tmp_path / "stub", {"0.50": 1, "0.90": "fail"})
    metadata = sequence_metadata_csv(tmp_path / "meta.csv", [
        qc_row(design=str(i), designed_residues="AVSV") for i in range(3)])
    out = tmp_path / "counts.csv"
    _run_main(monkeypatch, "--metadata-csv", str(metadata),
              "--output-csv", str(out), "--thresholds", "0.50,0.90",
              "--mmseqs-bin", str(stub))
    assert _read_counts(out) == [
        {"threshold": "0.50", "n_clusters": "1", "n_sequences": "3"},
        {"threshold": "0.90", "n_clusters": "3", "n_sequences": "3"},
    ]


@pytest.mark.logic
def test_main_clusters_at_two_sequences_which_is_the_boundary(
        tmp_path, monkeypatch):
    """`n_seqs < 2` is the skip, so two sequences is the smallest input that
    actually invokes mmseqs. One invocation per threshold, no more."""
    root = tmp_path / "stub"
    stub = fake_mmseqs(root, {"0.50": 1})
    metadata = sequence_metadata_csv(tmp_path / "meta.csv", [
        qc_row(design="0", designed_residues="AVSV"),
        qc_row(design="1", designed_residues="AVSV")])
    out = tmp_path / "counts.csv"
    _run_main(monkeypatch, "--metadata-csv", str(metadata),
              "--output-csv", str(out), "--thresholds", "0.50",
              "--mmseqs-bin", str(stub))
    assert len(mmseqs_calls(root)) == 1
    assert _read_counts(out) == [
        {"threshold": "0.50", "n_clusters": "1", "n_sequences": "2"}]


@pytest.mark.logic
def test_main_with_one_sequence_skips_clustering_and_reports_one_cluster(
        tmp_path, monkeypatch, capsys):
    """One sequence is one cluster at every threshold by definition, and
    mmseqs is never invoked. The row still has to be written, because the
    diversity plot reads this file and the stage publishes it."""
    root = tmp_path / "stub"
    stub = fake_mmseqs(root, {})
    metadata = sequence_metadata_csv(tmp_path / "meta.csv",
                                     [qc_row(designed_residues="AVSV")])
    out = tmp_path / "counts.csv"
    _run_main(monkeypatch, "--metadata-csv", str(metadata),
              "--output-csv", str(out), "--thresholds", "0.30,1.00",
              "--mmseqs-bin", str(stub))
    assert mmseqs_calls(root) == []
    assert _read_counts(out) == [
        {"threshold": "0.30", "n_clusters": "1", "n_sequences": "1"},
        {"threshold": "1.00", "n_clusters": "1", "n_sequences": "1"},
    ]
    assert "Only 1 sequences" in capsys.readouterr().out


@pytest.mark.logic
def test_main_with_rows_but_no_design_regions_reports_zero(
        tmp_path, monkeypatch):
    """Rows exist, so the early return does not fire, but nothing reached the
    FASTA. n_sequences is 0 rather than the row count, which is the honest
    answer: there was nothing to cluster."""
    metadata = sequence_metadata_csv(tmp_path / "meta.csv", [
        qc_row(design="0", designed_residues=""),
        qc_row(design="1", designed_residues="|")])
    out = tmp_path / "counts.csv"
    _run_main(monkeypatch, "--metadata-csv", str(metadata),
              "--output-csv", str(out), "--thresholds", "0.50")
    assert _read_counts(out) == [
        {"threshold": "0.50", "n_clusters": "0", "n_sequences": "0"}]


@pytest.mark.logic
def test_main_with_a_header_only_metadata_file_writes_a_header_only_csv(
        tmp_path, monkeypatch, capsys):
    """Zero rows: the design round produced no sequences at all. The output is
    the header alone, which is what the diversity plot's no-data path reads."""
    metadata = sequence_metadata_csv(tmp_path / "meta.csv", [])
    out = tmp_path / "counts.csv"
    _run_main(monkeypatch, "--metadata-csv", str(metadata),
              "--output-csv", str(out))
    assert out.read_text() == "threshold,n_clusters,n_sequences\n"
    assert "No sequences in metadata" in capsys.readouterr().out


@pytest.mark.logic
def test_main_with_an_absent_metadata_file_raises(tmp_path, monkeypatch):
    """No guard on the input path. The stage fails loudly, which is right for
    a missing upstream file, so this is pinned rather than wished away."""
    with pytest.raises(FileNotFoundError):
        _run_main(monkeypatch, "--metadata-csv", str(tmp_path / "absent.csv"),
                  "--output-csv", str(tmp_path / "counts.csv"))


@pytest.mark.logic
def test_the_default_thresholds_are_the_eight_in_the_committed_counts_file(
        tmp_path, monkeypatch):
    """mpnn_cluster_counts.csv in the reference set has rows 0.30 through 1.00
    in steps of 0.10, and the diversity plot's x-axis is built from them, so
    the default is part of the contract."""
    root = tmp_path / "stub"
    stub = fake_mmseqs(root, {})
    metadata = sequence_metadata_csv(tmp_path / "meta.csv", [
        qc_row(design="0", designed_residues="AVSV"),
        qc_row(design="1", designed_residues="KKDD")])
    out = tmp_path / "counts.csv"
    _run_main(monkeypatch, "--metadata-csv", str(metadata),
              "--output-csv", str(out), "--mmseqs-bin", str(stub))
    assert [r["threshold"] for r in _read_counts(out)] == [
        "0.30", "0.40", "0.50", "0.60", "0.70", "0.80", "0.90", "1.00"]


@pytest.mark.logic
def test_main_cleans_up_its_temporary_work_directory(tmp_path, monkeypatch):
    """The FASTA and the mmseqs tmp trees go to a mkdtemp directory, which the
    Nextflow work dir never sees. Leaving them behind fills /tmp over a
    campaign of 64 designs times 128 sequences."""
    import tempfile

    made: list[str] = []
    real_mkdtemp = tempfile.mkdtemp
    monkeypatch.setattr(tempfile, "mkdtemp",
                        lambda *a, **k: made.append(real_mkdtemp(*a, **k)) or made[-1])
    stub = fake_mmseqs(tmp_path / "stub", {"0.50": 1})
    metadata = sequence_metadata_csv(tmp_path / "meta.csv", [
        qc_row(design="0", designed_residues="AVSV"),
        qc_row(design="1", designed_residues="KKDD")])
    _run_main(monkeypatch, "--metadata-csv", str(metadata),
              "--output-csv", str(tmp_path / "counts.csv"),
              "--thresholds", "0.50", "--mmseqs-bin", str(stub))
    assert made and not Path(made[0]).exists()


@pytest.mark.logic
def test_parse_args_requires_the_metadata_csv(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["mpnn_cluster_sequences.py"])
    with pytest.raises(SystemExit):
        mc.parse_args()


@pytest.mark.logic
def test_parse_args_defaults_match_the_nextflow_invocation(monkeypatch):
    """MPNN_CLUSTER passes only --metadata-csv and --output-csv, so the
    threshold list and the binary name come from these defaults."""
    monkeypatch.setattr(sys, "argv", ["mpnn_cluster_sequences.py",
                                      "--metadata-csv", "m.csv"])
    args = mc.parse_args()
    assert args.output_csv == "mpnn_cluster_counts.csv"
    assert args.mmseqs_bin == "mmseqs"
    assert [float(t) for t in args.thresholds.split(",")] == [
        0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00]


@pytest.mark.logic
def test_the_stub_mmseqs_really_is_invoked_as_a_subprocess(tmp_path):
    """Guards the fixture itself. If the wrapper were not executable the
    module's subprocess.run would raise PermissionError and every test above
    would be asserting against a stub that never ran."""
    stub = fake_mmseqs(tmp_path / "stub", {})
    result = subprocess.run([str(stub), "easy-cluster", str(tmp_path / "x.fasta"),
                             str(tmp_path / "p"), str(tmp_path / "t"),
                             "--min-seq-id", "0.50"],
                            capture_output=True, text=True)
    assert result.returncode != 0  # x.fasta does not exist
    assert "No such file" in result.stderr
