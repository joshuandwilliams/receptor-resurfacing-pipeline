"""
tests/characterization/test_extract_passing.py
----------------------------------------------
Logic-tier tests for bin/extract_passing.py, which turns one sequence's
``aggregated_results.csv`` into its ``passing_summary.csv``.

Column names come from the producing code and from the committed reference
set under ``tests/negative_steering/example_output_files/``: the input
carries ``steered_*`` and ``reverted_*`` aggregates plus ``outcome``,
``n_pass``, ``n_seeds`` and ``rank_by_composite_score``, and the output is
exactly ``extract_passing.OUTPUT_FIELDS`` in that order.

Threshold strictness, stated once here and asserted below:

  * ``pae_pass_frac_median``, ``iptm_median`` and ``complex_plddt_median``
    trigger on ``value < minimum``.  A value exactly on the minimum passes.
  * ``ipae_median`` triggers on ``value > maximum``.  A value exactly on
    the maximum passes.
  * The reverted-mutation majority threshold is ``(n // 2) + 1``, i.e. a
    strict majority of seeds.  For n = 3 that is 2 seeds; for n = 2 it is
    both seeds.

Tier assignment itself lives in ``cross_sequence_summary._tier_for_row`` and
is tested in ``test_cross_sequence_summary.py``.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT / "bin") not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT / "bin"))

import extract_passing as ep  # noqa: E402

# A confidence block that trips nothing, so each test below can move one
# metric at a time.  Defaults come from
# PipelineInternalThresholds.default(): pae_pass_frac_min 0.10,
# iptm_min 0.30, complex_plddt_min 0.70, ipae_max 15.0.
_CLEAN_CONFIDENCE = {
    "pae_pass_frac_median": "0.50",
    "iptm_median": "0.90",
    "complex_plddt_median": "0.90",
    "ipae_median": "5.0",
}


def _aggregated_row(**overrides) -> dict:
    """One ``aggregated_results.csv`` row with both metric stages present.

    `steered_*` and `reverted_*` carry deliberately different values so a
    test can tell which stage `extract_row` sourced a column from.
    """
    row = {
        "cycle": "0",
        "pathway": "cycle_0",
        "sequence_group": "0",
        "design": "initial",
        "canonical_pdb": "/work/ab/cd/design_1_seq_0/cycle_0/p.pdb",
        "canonical_seed_index": "0",
        "canonical_seed_choice_warning": "",
        "n_seeds": "3",
        "n_pass": "3",
        "outcome": "no_reversion",
        "rank_by_composite_score": "1",
        "rank_by_ra_eff": "1",
        "steered_total_mutations_median": "3.0",
        "reverted_total_mutations_median": "1.0",
        "steered_ra_eff_vs_truth_median": "9.9",
        "reverted_ra_eff_vs_truth_median": "1.1",
        "steered_independent_receptor_rmsd_median": "1.5",
        "reverted_independent_receptor_rmsd_median": "0.5",
        "steered_independent_effector_rmsd_median": "0.6",
        "reverted_independent_effector_rmsd_median": "0.2",
        "steered_true_jaccard_median": "0.80",
        "reverted_true_jaccard_median": "0.20",
        "steered_wrong_jaccard_median": "0.10",
        "steered_n_shared_true_median": "12",
        "steered_n_shared_wrong_median": "1",
        "steered_n_contact_residues_median": "20",
        "steered_contact_residues_majority": "2,5,33",
        "steered_mutated_contact_positions_majority": "5",
        "steered_intact_core_majority": "1",
    }
    for key, value in _CLEAN_CONFIDENCE.items():
        row[f"steered_{key}"] = value
        row[f"reverted_{key}"] = value
    row.update(overrides)
    return row


# ── _get ─────────────────────────────────────────────────────────────


@pytest.mark.logic
def test_get_prefers_the_primary_key_and_falls_back_when_blank():
    """An empty string counts as absent, so the fallback key is used."""
    assert ep._get({"a": "3", "b": "7"}, "a", "b") == "3"
    assert ep._get({"a": "", "b": "7"}, "a", "b") == "7"
    assert ep._get({"a": None, "b": "9"}, "a", "b") == "9"
    assert ep._get({}, "a", "b") == ""
    assert ep._get({}, "a") == ""
    # Zero is a real value, not a blank, so it is not replaced.
    assert ep._get({"a": "0", "b": "7"}, "a", "b") == "0"


# ── Confidence flag thresholds ───────────────────────────────────────


@pytest.mark.logic
def test_confidence_flag_is_ok_when_every_median_is_comfortable():
    assert ep._compute_confidence_flag(dict(_CLEAN_CONFIDENCE)) == "ok"


@pytest.mark.logic
def test_confidence_flag_minimum_thresholds_are_strict():
    """A median exactly on its minimum passes; a hair below trips.

    The comparison is `value < minimum`, so equality is a pass.
    """
    on_the_bound = {"pae_pass_frac_median": "0.10", "iptm_median": "0.30",
                    "complex_plddt_median": "0.70", "ipae_median": "5.0"}
    assert ep._compute_confidence_flag(on_the_bound) == "ok"
    assert ep._compute_confidence_flag(
        {**_CLEAN_CONFIDENCE, "iptm_median": "0.2999"}) == "low_iptm"
    assert ep._compute_confidence_flag(
        {**_CLEAN_CONFIDENCE, "pae_pass_frac_median": "0.0999"}
    ) == "low_pass_frac"
    assert ep._compute_confidence_flag(
        {**_CLEAN_CONFIDENCE, "complex_plddt_median": "0.6999"}
    ) == "low_plddt"


@pytest.mark.logic
def test_confidence_flag_ipae_maximum_is_strict():
    """ipae trips on `value > maximum`, so 15.0 exactly is a pass."""
    assert ep._compute_confidence_flag(
        {**_CLEAN_CONFIDENCE, "ipae_median": "15.0"}) == "ok"
    assert ep._compute_confidence_flag(
        {**_CLEAN_CONFIDENCE, "ipae_median": "15.0001"}) == "high_ipae"


@pytest.mark.logic
def test_confidence_flag_reports_multiple_when_more_than_one_trips():
    """One trigger is named; two or more collapse to `multiple`."""
    assert ep._compute_confidence_flag(
        {**_CLEAN_CONFIDENCE, "iptm_median": "0.1",
         "ipae_median": "99.0"}) == "multiple"


@pytest.mark.logic
def test_confidence_flag_treats_a_missing_metric_as_passing():
    """Absent columns default to a passing value rather than tripping.

    An empty dict reads every metric through its `.get(key, default)`
    fallback, so a sequence whose metrics were never computed is reported
    `ok` rather than flagged.
    """
    assert ep._compute_confidence_flag({}) == "ok"


@pytest.mark.logic
def test_confidence_flag_ignores_an_unparseable_metric():
    """A non-numeric value is skipped rather than crashing the run."""
    assert ep._compute_confidence_flag(
        {**_CLEAN_CONFIDENCE, "iptm_median": "abc"}) == "ok"


@pytest.mark.logic
def test_confidence_flag_thresholds_can_be_overridden():
    """Explicit thresholds bypass PipelineInternalThresholds entirely."""
    assert ep._compute_confidence_flag(
        dict(_CLEAN_CONFIDENCE), iptm_min=0.95, pass_frac_min=0.0,
        complex_plddt_min=0.0, ipae_max=1e9) == "low_iptm"


@pytest.mark.logic
@pytest.mark.xfail(strict=False, reason=(
    "Known defect, reported rather than fixed: nan compares False against "
    "every threshold, so a nan median is flagged 'ok'.  Changing it would "
    "change the confidence_flag column of committed passing_summary.csv "
    "output."))
def test_a_nan_metric_must_not_be_reported_as_confident():
    """A nan median is not evidence of confidence, so it should flag.

    `float('nan') < 0.30` is False, as is `float('nan') > 15.0`, so every
    one of the four comparisons in `_compute_confidence_flag` skips a nan
    and the row comes out `ok`.  That reads as "the typical prediction
    passed the confidence filter" when in fact no usable number exists.

    Contrast `cross_sequence_summary._try_float`, which rejects nan
    explicitly so a nan ra_eff cannot be ranked.  This test asserts the
    same standard here.
    """
    assert ep._compute_confidence_flag(
        {**_CLEAN_CONFIDENCE, "iptm_median": "nan"}) != "ok"
    assert ep._compute_confidence_flag(
        {**_CLEAN_CONFIDENCE, "ipae_median": "nan"}) != "ok"


# ── Reverted-mutation summarisation ──────────────────────────────────


@pytest.mark.logic
def test_mutation_majority_needs_a_strict_majority_of_seeds():
    """Threshold is (n // 2) + 1, so n = 3 needs 2 seeds.

    Position 5 is in all three seeds, 7 in two, and 22 and 30 in one each.
    Majority is {5, 7}; the union is {5, 7, 22, 30}.  The seeds differ, so
    the identical flag is "0".
    """
    rows = [{"reverted_mutations_chimerax": "/A:5,7,22"},
            {"reverted_mutations_chimerax": "/A:5,7"},
            {"reverted_mutations_chimerax": "/A:5,30"}]
    assert ep._summarize_mutations_across_seeds(
        rows, "reverted_mutations_chimerax") == (
            "/A:5,7", "/A:5,7,22,30", "0")


@pytest.mark.logic
def test_mutation_majority_for_two_seeds_requires_both():
    """n = 2 gives threshold (2 // 2) + 1 = 2, so a 1-of-2 position is out.

    This is the case where "majority" is strictly more than half rather
    than at least half.
    """
    rows = [{"reverted_mutations_chimerax": "/B:5,7"},
            {"reverted_mutations_chimerax": "/B:7,9"}]
    assert ep._summarize_mutations_across_seeds(
        rows, "reverted_mutations_chimerax") == ("/B:7", "/B:5,7,9", "0")


@pytest.mark.logic
def test_mutation_summary_chain_prefix_comes_from_the_first_non_empty_seed():
    rows = [{"reverted_mutations_chimerax": ""},
            {"reverted_mutations_chimerax": "/C:4"}]
    maj, every, _ = ep._summarize_mutations_across_seeds(
        rows, "reverted_mutations_chimerax")
    assert maj.startswith("/C:")
    assert every.startswith("/C:")


@pytest.mark.logic
def test_mutation_summary_identical_flag():
    """"1" when every seed reverted identically, "0" when they differ, and
    "" when there is only one seed with data (nothing to compare)."""
    same = [{"reverted_mutations_chimerax": "/A:5,7"}] * 3
    assert ep._summarize_mutations_across_seeds(
        same, "reverted_mutations_chimerax") == ("/A:5,7", "/A:5,7", "1")
    single = [{"reverted_mutations_chimerax": "/A:5"}]
    assert ep._summarize_mutations_across_seeds(
        single, "reverted_mutations_chimerax") == ("/A:5", "/A:5", "")


@pytest.mark.logic
def test_mutation_summary_of_no_usable_seed_is_three_blanks():
    """Empty values and strings without a ':' are skipped entirely."""
    assert ep._summarize_mutations_across_seeds(
        [], "reverted_mutations_chimerax") == ("", "", "")
    assert ep._summarize_mutations_across_seeds(
        [{"reverted_mutations_chimerax": ""}],
        "reverted_mutations_chimerax") == ("", "", "")
    assert ep._summarize_mutations_across_seeds(
        [{"reverted_mutations_chimerax": "no-colon-here"}],
        "reverted_mutations_chimerax") == ("", "", "")


@pytest.mark.logic
def test_muts_aa_merges_tokens_across_seeds_and_sorts_by_position():
    """No single seed lists every retained position, so all are merged.

    Seed 1 kept 5 and 7, seed 2 kept 22.  Asking for [5, 7, 22] returns
    all three tokens, ordered by residue number whatever order they were
    requested in.
    """
    rows = [{"reverted_mutations_aa": "K5E V7D",
             "steered_mutations_aa": "K5E V7D M22R"},
            {"reverted_mutations_aa": "M22R",
             "steered_mutations_aa": "K5E V7D M22R"}]
    assert ep._muts_aa_for_positions(
        rows, "reverted_mutations_aa", [5, 7, 22]) == "K5E V7D M22R"
    assert ep._muts_aa_for_positions(
        rows, "reverted_mutations_aa", [22, 5]) == "K5E M22R"
    assert ep._muts_aa_for_positions(rows, "reverted_mutations_aa", []) == ""


@pytest.mark.logic
def test_muts_aa_falls_back_to_the_steered_set_for_unknown_positions():
    """steered_mutations_aa always carries every mutated position."""
    rows = [{"steered_mutations_aa": "K5E V7D M22R"}]
    assert ep._muts_aa_for_positions(
        rows, "reverted_mutations_aa", [7]) == "V7D"


@pytest.mark.logic
def test_muts_aa_skips_tokens_it_cannot_parse():
    rows = [{"reverted_mutations_aa": "XX K5E 12"}]
    assert ep._muts_aa_for_positions(
        rows, "reverted_mutations_aa", [5]) == "K5E"


# ── extract_row stage selection ──────────────────────────────────────


@pytest.mark.logic
def test_extract_row_emits_exactly_the_output_schema():
    """Every OUTPUT_FIELDS column except the runtime is produced here.

    run_one_runtime_sec is stamped by main() from a sidecar file, not by
    extract_row, so it is the one documented omission.
    """
    out = ep.extract_row(_aggregated_row())
    assert set(out) == set(ep.OUTPUT_FIELDS) - {"run_one_runtime_sec"}


@pytest.mark.logic
@pytest.mark.parametrize("outcome", ["pose_holds", "pose_collapses",
                                     "new_contamination"])
def test_reversion_attempted_outcomes_source_metrics_from_reverted(outcome):
    """When reversion ran, the reverted prediction is the final design.

    Structural and confidence columns come from `reverted_*`.  Interface
    similarity still comes from `steered_*`, because the designed interface
    is a property of the sequence rather than of the reverted pose.
    """
    out = ep.extract_row(_aggregated_row(outcome=outcome))
    assert out["ra_eff_vs_truth_median"] == "1.1"
    assert out["independent_receptor_rmsd_median"] == "0.5"
    assert out["total_mutations_median"] == "1.0"
    assert out["true_jaccard_median"] == "0.80"
    assert out["wrong_jaccard_median"] == "0.10"


@pytest.mark.logic
@pytest.mark.parametrize("outcome", ["no_reversion", "", "singleton"])
def test_outcomes_without_reversion_source_metrics_from_steered(outcome):
    """No reversion means the steered prediction IS the final design.

    `no_reversion` is not a cold-start marker, it only means reversion was
    never triggered, so the steered metrics are the real ones.
    """
    out = ep.extract_row(_aggregated_row(outcome=outcome))
    assert out["ra_eff_vs_truth_median"] == "9.9"
    assert out["independent_receptor_rmsd_median"] == "1.5"
    assert out["total_mutations_median"] == "3.0"
    assert out["true_jaccard_median"] == "0.80"


@pytest.mark.logic
def test_extract_row_falls_back_to_steered_when_reverted_is_blank():
    """A pose_holds row missing its reverted aggregate still gets a value.

    Only the columns with a documented fallback behave this way;
    ra_eff_vs_truth_median has no fallback and stays blank.
    """
    row = _aggregated_row(outcome="pose_holds",
                          reverted_total_mutations_median="",
                          reverted_ra_eff_vs_truth_median="",
                          reverted_iptm_median="")
    out = ep.extract_row(row)
    assert out["total_mutations_median"] == "3.0"
    assert out["iptm_median"] == "0.90"
    assert out["ra_eff_vs_truth_median"] == ""


@pytest.mark.logic
def test_extract_row_copies_n_pass_and_outcome_through_untouched():
    """n_pass is path-agnostic and authoritative, so it is not recomputed.

    A `no_reversion` row with n_pass below n_seeds keeps both numbers, and
    the tier decision downstream is left to read them.
    """
    out = ep.extract_row(_aggregated_row(outcome="no_reversion", n_pass="2"))
    assert (out["n_pass"], out["n_seeds"]) == ("2", "3")
    assert out["outcome"] == "no_reversion"


@pytest.mark.logic
def test_extract_row_leaves_mutation_columns_blank_without_per_seed_rows():
    out = ep.extract_row(_aggregated_row(), per_seed_rows=None)
    for col in ("steered_mutations_chimerax", "steered_mutations_aa",
                "reverted_mutations_majority_chimerax",
                "reverted_mutations_majority_aa",
                "reverted_mutations_all_chimerax",
                "reverted_mutations_all_aa",
                "reverted_mutations_seeds_identical"):
        assert out[col] == "", col


@pytest.mark.logic
def test_extract_row_populates_mutation_columns_from_per_seed_rows():
    """Steered mutations are identical across seeds, so any seed supplies
    them.  Reverted mutations are summarised across seeds."""
    per_seed = [
        {"steered_mutations_chimerax": "/A:5,7,22",
         "steered_mutations_aa": "K5E V7D M22R",
         "reverted_mutations_chimerax": "/A:5,7",
         "reverted_mutations_aa": "K5E V7D"},
        {"steered_mutations_chimerax": "/A:5,7,22",
         "steered_mutations_aa": "K5E V7D M22R",
         "reverted_mutations_chimerax": "/A:5,22",
         "reverted_mutations_aa": "K5E M22R"},
    ]
    out = ep.extract_row(_aggregated_row(), per_seed_rows=per_seed)
    assert out["steered_mutations_chimerax"] == "/A:5,7,22"
    assert out["steered_mutations_aa"] == "K5E V7D M22R"
    # n = 2, so the majority threshold is 2: only position 5 qualifies.
    assert out["reverted_mutations_majority_chimerax"] == "/A:5"
    assert out["reverted_mutations_majority_aa"] == "K5E"
    assert out["reverted_mutations_all_chimerax"] == "/A:5,7,22"
    assert out["reverted_mutations_all_aa"] == "K5E V7D M22R"
    assert out["reverted_mutations_seeds_identical"] == "0"


@pytest.mark.logic
def test_extract_row_blanks_the_majority_columns_when_no_position_qualifies():
    """Two seeds reverting disjoint sets leave no strict majority.

    The majority chimerax string is then empty, and its AA companion has no
    positions to look up, so both come out blank while the union keeps
    everything.
    """
    per_seed = [
        {"steered_mutations_chimerax": "/A:5,22",
         "steered_mutations_aa": "K5E M22R",
         "reverted_mutations_chimerax": "/A:5",
         "reverted_mutations_aa": "K5E"},
        {"steered_mutations_chimerax": "/A:5,22",
         "steered_mutations_aa": "K5E M22R",
         "reverted_mutations_chimerax": "/A:22",
         "reverted_mutations_aa": "M22R"},
    ]
    out = ep.extract_row(_aggregated_row(), per_seed_rows=per_seed)
    assert out["reverted_mutations_majority_chimerax"] == ""
    assert out["reverted_mutations_majority_aa"] == ""
    assert out["reverted_mutations_all_chimerax"] == "/A:5,22"


@pytest.mark.logic
def test_extract_row_computes_the_confidence_flag_from_its_own_output():
    """The flag reads the medians already selected for this row, so the
    stage choice feeds the flag."""
    out = ep.extract_row(
        _aggregated_row(outcome="pose_holds", reverted_iptm_median="0.05"))
    assert out["iptm_median"] == "0.05"
    assert out["confidence_flag"] == "low_iptm"


# ── main() ───────────────────────────────────────────────────────────


def _write_csv(path: Path, rows: list[dict], fieldnames=None) -> Path:
    fieldnames = fieldnames or sorted({k for r in rows for k in r})
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return path


def _run_main(monkeypatch, *argv_tail) -> int:
    monkeypatch.setattr(sys, "argv", ["extract_passing.py", *argv_tail])
    return ep.main()


@pytest.mark.logic
def test_main_filters_on_outcome_and_on_a_populated_rank(
        tmp_path, monkeypatch):
    """Two gates, both required.

      * `outcome` must be "", "no_reversion" or "pose_holds".
      * `rank_by_composite_score` must be non-empty.

    So the `pose_collapses` row and the unranked `no_reversion` row are
    both dropped, and the survivors come out sorted by rank ascending.
    """
    rows = [
        _aggregated_row(sequence_group="0", outcome="no_reversion",
                        rank_by_composite_score="2"),
        _aggregated_row(sequence_group="1", outcome="pose_holds",
                        rank_by_composite_score="1"),
        _aggregated_row(sequence_group="2", outcome="pose_collapses",
                        rank_by_composite_score="3"),
        _aggregated_row(sequence_group="3", outcome="new_contamination",
                        rank_by_composite_score="4"),
        _aggregated_row(sequence_group="4", outcome="no_reversion",
                        rank_by_composite_score=""),
        _aggregated_row(sequence_group="5", outcome="",
                        rank_by_composite_score="5"),
    ]
    agg = _write_csv(tmp_path / "aggregated_results.csv", rows)
    assert _run_main(monkeypatch, "--input", str(agg)) == 0

    with open(tmp_path / "passing_summary.csv") as f:
        got = list(csv.DictReader(f))
    assert [r["sequence_group"] for r in got] == ["1", "0", "5"]
    assert [r["rank_by_composite_score"] for r in got] == ["1", "2", "5"]


@pytest.mark.logic
def test_main_writes_the_output_schema_in_order(tmp_path, monkeypatch):
    agg = _write_csv(tmp_path / "aggregated_results.csv", [_aggregated_row()])
    _run_main(monkeypatch, "--input", str(agg))
    with open(tmp_path / "passing_summary.csv") as f:
        assert next(csv.reader(f)) == ep.OUTPUT_FIELDS


@pytest.mark.logic
def test_main_writes_a_header_only_file_when_nothing_passes(
        tmp_path, monkeypatch, capsys):
    """An empty passing_summary.csv is a valid result, not an error.

    cross_sequence_summary's tier-none fallback exists precisely to put
    something in the cohort table for these sequences.
    """
    agg = _write_csv(tmp_path / "aggregated_results.csv",
                     [_aggregated_row(outcome="pose_collapses")])
    assert _run_main(monkeypatch, "--input", str(agg)) == 0
    with open(tmp_path / "passing_summary.csv") as f:
        rows = list(csv.DictReader(f))
    assert rows == []
    assert "no ranked rows found" in capsys.readouterr().out


@pytest.mark.logic
def test_main_joins_per_seed_rows_on_cycle_pathway_and_sequence_group(
        tmp_path, monkeypatch, capsys):
    """The join key is the (cycle, pathway, sequence_group) triple.

    The second per-seed row below belongs to sequence_group 9, so it must
    not leak into sequence_group 0's mutation columns.
    """
    agg = _write_csv(tmp_path / "aggregated_results.csv",
                     [_aggregated_row(sequence_group="0")])
    _write_csv(tmp_path / "raw_per_seed_results.csv", [
        {"cycle": "0", "pathway": "cycle_0", "sequence_group": "0",
         "steered_mutations_chimerax": "/A:5,7",
         "steered_mutations_aa": "K5E V7D",
         "reverted_mutations_chimerax": "/A:5",
         "reverted_mutations_aa": "K5E"},
        {"cycle": "0", "pathway": "cycle_0", "sequence_group": "9",
         "steered_mutations_chimerax": "/A:99",
         "steered_mutations_aa": "W99A",
         "reverted_mutations_chimerax": "/A:99",
         "reverted_mutations_aa": "W99A"},
    ])
    _run_main(monkeypatch, "--input", str(agg))
    with open(tmp_path / "passing_summary.csv") as f:
        got = list(csv.DictReader(f))
    assert got[0]["steered_mutations_chimerax"] == "/A:5,7"
    assert got[0]["reverted_mutations_all_chimerax"] == "/A:5"
    assert "99" not in got[0]["steered_mutations_aa"]
    assert "Loaded 2 per-seed rows" in capsys.readouterr().out


@pytest.mark.logic
def test_main_warns_and_blanks_mutations_without_a_per_seed_file(
        tmp_path, monkeypatch, capsys):
    agg = _write_csv(tmp_path / "aggregated_results.csv", [_aggregated_row()])
    _run_main(monkeypatch, "--input", str(agg))
    with open(tmp_path / "passing_summary.csv") as f:
        got = list(csv.DictReader(f))
    assert got[0]["steered_mutations_chimerax"] == ""
    assert "mutation columns" in capsys.readouterr().out


@pytest.mark.logic
def test_main_accepts_an_explicit_per_seed_path(tmp_path, monkeypatch):
    agg = _write_csv(tmp_path / "aggregated_results.csv", [_aggregated_row()])
    per_seed = _write_csv(tmp_path / "elsewhere" / "seeds.csv", [
        {"cycle": "0", "pathway": "cycle_0", "sequence_group": "0",
         "steered_mutations_chimerax": "/A:11",
         "steered_mutations_aa": "A11G",
         "reverted_mutations_chimerax": "/A:11",
         "reverted_mutations_aa": "A11G"},
    ])
    _run_main(monkeypatch, "--input", str(agg), "--per-seed", str(per_seed))
    with open(tmp_path / "passing_summary.csv") as f:
        got = list(csv.DictReader(f))
    assert got[0]["steered_mutations_chimerax"] == "/A:11"


@pytest.mark.logic
def test_main_stamps_the_runtime_sidecar_onto_every_row(tmp_path, monkeypatch):
    """One `negsteer run` invocation covers one sequence, so every row of
    that sequence's file carries the same wall-clock value."""
    agg = _write_csv(tmp_path / "aggregated_results.csv", [
        _aggregated_row(sequence_group="0", rank_by_composite_score="1"),
        _aggregated_row(sequence_group="1", rank_by_composite_score="2"),
    ])
    (tmp_path / "run_one_runtime_sec.txt").write_text("231.5\n")
    _run_main(monkeypatch, "--input", str(agg))
    with open(tmp_path / "passing_summary.csv") as f:
        got = list(csv.DictReader(f))
    assert [r["run_one_runtime_sec"] for r in got] == ["231.5", "231.5"]


@pytest.mark.logic
def test_main_blanks_a_corrupt_runtime_sidecar(tmp_path, monkeypatch):
    """The sidecar is validated as a number so garbage is not propagated."""
    agg = _write_csv(tmp_path / "aggregated_results.csv", [_aggregated_row()])
    (tmp_path / "run_one_runtime_sec.txt").write_text("not-a-number\n")
    _run_main(monkeypatch, "--input", str(agg))
    with open(tmp_path / "passing_summary.csv") as f:
        got = list(csv.DictReader(f))
    assert got[0]["run_one_runtime_sec"] == ""


@pytest.mark.logic
def test_main_sorts_unrankable_rows_last(tmp_path, monkeypatch):
    """A non-integer rank sorts behind every real rank via the 10**9
    sentinel, rather than raising."""
    agg = _write_csv(tmp_path / "aggregated_results.csv", [
        _aggregated_row(sequence_group="0", rank_by_composite_score="junk"),
        _aggregated_row(sequence_group="1", rank_by_composite_score="7"),
    ])
    _run_main(monkeypatch, "--input", str(agg))
    with open(tmp_path / "passing_summary.csv") as f:
        got = list(csv.DictReader(f))
    assert [r["sequence_group"] for r in got] == ["1", "0"]


@pytest.mark.logic
def test_main_honours_an_explicit_output_path(tmp_path, monkeypatch):
    agg = _write_csv(tmp_path / "aggregated_results.csv", [_aggregated_row()])
    out = tmp_path / "nested" / "custom.csv"
    assert _run_main(monkeypatch, "--input", str(agg),
                     "--output", str(out)) == 0
    assert out.exists()
    assert not (tmp_path / "passing_summary.csv").exists()


@pytest.mark.logic
def test_main_returns_two_when_the_input_is_missing(
        tmp_path, monkeypatch, capsys):
    rc = _run_main(monkeypatch, "--input", str(tmp_path / "absent.csv"))
    assert rc == 2
    assert "not found" in capsys.readouterr().err
