"""
tests/characterization/test_cross_sequence_summary.py
-----------------------------------------------------
Logic-tier tests for bin/cross_sequence_summary.py, which stacks one
``passing_summary.csv`` per MPNN sequence into one cross-sequence table.

Column names are taken from ``extract_passing.OUTPUT_FIELDS`` (the producing
code) and from the committed reference set under
``tests/negative_steering/example_output_files/``, not invented here.

Threshold strictness, stated once and asserted below:

  * Tier A is ``n_pass == n_seeds``.  Exact equality, so 3/3 is A and 4/3
    is "none" rather than A.
  * Tier B is ``1 < n_pass < n_seeds``, strict on both sides.  2/3 is B;
    1/3 is not B, and 3/3 is not B.
  * Tier C is ``n_pass == 1``, exact.  1/1 is A, not C, because the
    ``n_pass == n_seeds`` test is checked first.
  * Anything else, including ``n_pass == 0``, an unparseable count or
    ``n_seeds <= 0``, is "none".

n_pass is path-agnostic, so ``outcome == "no_reversion"`` is not a tier
shortcut in either direction.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT / "bin") not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT / "bin"))

import cross_sequence_summary as css  # noqa: E402
import extract_passing as ep  # noqa: E402


def _passing_row(**overrides) -> dict:
    """One ``passing_summary.csv`` row, every OUTPUT_FIELDS column present."""
    row = {col: "" for col in ep.OUTPUT_FIELDS}
    row.update({
        "rank_by_composite_score": "1",
        "rank_by_ra_eff": "1",
        "cycle": "0",
        "pathway": "cycle_0",
        "sequence_group": "0",
        "design": "initial",
        "n_seeds": "3",
        "n_pass": "3",
        "outcome": "no_reversion",
        "ra_eff_vs_truth_median": "2.0",
        "true_jaccard_median": "0.80",
        "confidence_flag": "ok",
    })
    row.update(overrides)
    return row


def _write_passing_summary(directory: Path, rows: list[dict]) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / "passing_summary.csv"
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=ep.OUTPUT_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    return path


def _read_output(path: Path) -> list[dict]:
    with open(path) as f:
        return list(csv.DictReader(f))


# ── Numeric parsing ──────────────────────────────────────────────────


@pytest.mark.logic
def test_try_float_parses_numbers_and_rejects_non_finite():
    """nan and both infinities are not valid ranking values.

    This is the gate that keeps a nan metric out of every rank and score,
    so it is asserted explicitly rather than inferred.
    """
    assert css._try_float("2.5") == 2.5
    assert css._try_float(2) == 2.0
    assert css._try_float(" 3.0 ") == 3.0
    assert css._try_float(None) is None
    assert css._try_float("") is None
    assert css._try_float("   ") is None
    assert css._try_float("abc") is None
    assert css._try_float("nan") is None
    assert css._try_float(float("nan")) is None
    assert css._try_float("inf") is None
    assert css._try_float("-inf") is None


@pytest.mark.logic
def test_try_int_tolerates_a_float_string_by_truncating():
    assert css._try_int("3") == 3
    assert css._try_int(5) == 5
    assert css._try_int("3.0") == 3
    assert css._try_int("3.7") == 3
    assert css._try_int("") is None
    assert css._try_int(None) is None
    assert css._try_int("abc") is None


# ── Tier assignment ─────────────────────────────────────────────────


@pytest.mark.logic
@pytest.mark.parametrize("n_pass,n_seeds,expected", [
    # Tier A: exact equality, at any n_seeds.
    (3, 3, "A"),
    (2, 2, "A"),
    (1, 1, "A"),
    (5, 5, "A"),
    # Tier B: strictly between 1 and n_seeds.  2/3 is the first case that
    # can exist, since B needs n_seeds >= 3.
    (2, 3, "B"),
    (3, 4, "B"),
    (2, 4, "B"),
    # Tier C: exactly one passing seed, and not all of them.
    (1, 3, "C"),
    (1, 2, "C"),
    (1, 5, "C"),
    # None: no passing seed, or more passes than seeds (inconsistent).
    (0, 3, "none"),
    (0, 1, "none"),
    (4, 3, "none"),
])
def test_tier_for_row_on_each_side_of_each_threshold(n_pass, n_seeds,
                                                     expected):
    """Every boundary of the n_pass / n_seeds tiering, stated as a case.

    The checks run in order A, B, C, so 1/1 reaches the A test first and
    comes out A rather than C.  B is strict on both sides: `1 < n_pass`
    excludes 1/3 and `n_pass < n_seeds` excludes 3/3.
    """
    assert css._tier_for_row(
        {"n_pass": str(n_pass), "n_seeds": str(n_seeds)}) == expected


@pytest.mark.logic
def test_tier_is_none_when_the_counts_are_unusable():
    """Missing, blank, non-numeric or non-positive counts give up on a tier.

    "Insufficient data to place in any tier" is a deliberate "none" rather
    than a guess, so cross-ranking skips the row.
    """
    assert css._tier_for_row({}) == "none"
    assert css._tier_for_row({"n_pass": "", "n_seeds": "3"}) == "none"
    assert css._tier_for_row({"n_pass": "3", "n_seeds": ""}) == "none"
    assert css._tier_for_row({"n_pass": "abc", "n_seeds": "3"}) == "none"
    assert css._tier_for_row({"n_pass": "0", "n_seeds": "0"}) == "none"
    assert css._tier_for_row({"n_pass": "3", "n_seeds": "0"}) == "none"
    assert css._tier_for_row({"n_pass": "3", "n_seeds": "-1"}) == "none"


@pytest.mark.logic
def test_a_nan_pass_count_must_not_be_treated_as_passing():
    """A nan n_pass is not a passing seed count.

    `_try_int` cannot parse "nan" through either int() or int(float()), so
    the row lands in tier "none" and is kept out of the cross-rank.
    """
    assert css._tier_for_row({"n_pass": "nan", "n_seeds": "3"}) == "none"
    assert css._tier_for_row({"n_pass": "nan", "n_seeds": "nan"}) == "none"


@pytest.mark.logic
def test_no_reversion_is_not_a_tier_shortcut():
    """The outcome label is ignored; only n_pass / n_seeds decide.

    A wrong-placement group can carry outcome "no_reversion" with every
    seed in no_data or pose_collapses, so n_pass is 0 and the tier is
    "none".  A cold-start group where every seed passed carries the same
    label with n_pass == n_seeds and is tier A.
    """
    wrong_placement = {"outcome": "no_reversion", "n_pass": "0",
                       "n_seeds": "3"}
    cold_start_all_pass = {"outcome": "no_reversion", "n_pass": "3",
                           "n_seeds": "3"}
    assert css._tier_for_row(wrong_placement) == "none"
    assert css._tier_for_row(cold_start_all_pass) == "A"
    # A failing outcome label with passing counts still tiers on the counts.
    assert css._tier_for_row({"outcome": "pose_collapses", "n_pass": "2",
                              "n_seeds": "3"}) == "B"


@pytest.mark.logic
def test_tier_order_puts_a_before_b_before_c_before_none():
    assert css._TIER_ORDER == {"A": 0, "B": 1, "C": 2, "none": 3}


# ── Score and rank helpers ──────────────────────────────────────────


@pytest.mark.logic
def test_composite_score_is_true_jaccard_minus_a_twentieth_of_ra_eff():
    """composite = true_jaccard_median - 0.05 x ra_eff_vs_truth_median."""
    assert css._composite_score({"ra_eff_vs_truth_median": "2.0",
                                 "true_jaccard_median": "0.80"}) == (
        pytest.approx(0.7))
    assert css._composite_score({"ra_eff_vs_truth_median": "10.0",
                                 "true_jaccard_median": "0.50"}) == (
        pytest.approx(0.0))


@pytest.mark.logic
def test_composite_score_is_none_when_either_component_is_unusable():
    """A nan or blank component yields no score, so the row cannot rank."""
    assert css._composite_score({"true_jaccard_median": "0.8"}) is None
    assert css._composite_score({"ra_eff_vs_truth_median": "2.0"}) is None
    assert css._composite_score({"ra_eff_vs_truth_median": "nan",
                                 "true_jaccard_median": "0.8"}) is None
    assert css._composite_score({"ra_eff_vs_truth_median": "2.0",
                                 "true_jaccard_median": "nan"}) is None


@pytest.mark.logic
def test_median_ra_eff_sorts_an_unusable_value_last():
    """Lower is better, so a missing or nan median becomes +inf."""
    assert css._median_ra_eff({"ra_eff_vs_truth_median": "2.0"}) == 2.0
    assert css._median_ra_eff({"ra_eff_vs_truth_median": ""}) == float("inf")
    assert css._median_ra_eff({"ra_eff_vs_truth_median": "nan"}) == (
        float("inf"))
    assert css._median_ra_eff({}) == float("inf")


@pytest.mark.logic
def test_within_sequence_rank_sentinel_sorts_unrankable_rows_last():
    assert css._within_sequence_rank({"rank_by_composite_score": "3"}) == 3
    assert css._within_sequence_rank({}) == 10 ** 9
    assert css._within_sequence_rank(
        {"rank_by_composite_score": "x"}) == 10 ** 9


# ── Representative pick ─────────────────────────────────────────────


@pytest.mark.logic
def test_pick_representative_prefers_tier_over_a_better_median():
    """Reliability first, accuracy second.

    The tier-B row has a much better (lower) ra_eff than the tier-A row and
    still loses, because more seeds independently landing the pose beats a
    single lucky median.
    """
    tier_a = _passing_row(sequence_group="A", n_pass="3", n_seeds="3",
                          ra_eff_vs_truth_median="4.0")
    tier_b = _passing_row(sequence_group="B", n_pass="2", n_seeds="3",
                          ra_eff_vs_truth_median="0.5")
    rep = css._pick_representative([tier_b, tier_a])
    assert rep["sequence_group"] == "A"
    assert rep["_picked_tier"] == "A"


@pytest.mark.logic
def test_pick_representative_breaks_a_tier_tie_on_median_not_composite():
    """Within one tier, the lower ra_eff wins even with a worse jaccard.

    The losing row has the higher composite score (0.95 - 0.05 x 3.0 =
    0.80 against 0.10 - 0.05 x 1.0 = 0.05), which is exactly the case the
    switch from composite to median was made for.
    """
    better_pose = _passing_row(sequence_group="pose",
                               ra_eff_vs_truth_median="1.0",
                               true_jaccard_median="0.10")
    better_overlap = _passing_row(sequence_group="overlap",
                                  ra_eff_vs_truth_median="3.0",
                                  true_jaccard_median="0.95")
    rep = css._pick_representative([better_overlap, better_pose])
    assert rep["sequence_group"] == "pose"


@pytest.mark.logic
def test_pick_representative_returns_none_for_no_rows_or_only_none_tier():
    assert css._pick_representative([]) is None
    assert css._pick_representative(
        [_passing_row(n_pass="0", n_seeds="3")]) is None


@pytest.mark.logic
def test_pick_representative_does_not_mutate_its_input():
    """The pick is returned as an annotated copy."""
    row = _passing_row()
    rep = css._pick_representative([row])
    assert "_picked_tier" in rep
    assert "_picked_tier" not in row


# ── Tier-none fallback from aggregated_results.csv ──────────────────


def _aggregated_csv(directory: Path, rows: list[dict]) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / "aggregated_results.csv"
    fieldnames = sorted({k for r in rows for k in r})
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return path


@pytest.mark.logic
def test_fallback_picks_the_best_raw_composite_ignoring_the_eligibility_gate(
        tmp_path):
    """With no passing rows, the best aggregated row stands in.

    rank_by_composite_score is blank on all three rows below, which is what
    the engine's eligibility gate does, and the fallback ignores that.  The
    winner is the highest `true_jaccard - 0.05 x ra_eff`:

        row 0:  0.80 - 0.05 x 20.0 = -0.20
        row 1:  0.60 - 0.05 x  6.0 =  0.30   <- best
        row 2:  0.10 - 0.05 x  6.0 = -0.20
    """
    rows = [
        {"outcome": "no_reversion", "sequence_group": "0",
         "steered_ra_eff_vs_truth_median": "20.0",
         "steered_true_jaccard_median": "0.80"},
        {"outcome": "no_reversion", "sequence_group": "1",
         "steered_ra_eff_vs_truth_median": "6.0",
         "steered_true_jaccard_median": "0.60"},
        {"outcome": "no_reversion", "sequence_group": "2",
         "steered_ra_eff_vs_truth_median": "6.0",
         "steered_true_jaccard_median": "0.10"},
    ]
    _aggregated_csv(tmp_path, rows)
    rep = css._pick_rep_from_aggregated(tmp_path / "passing_summary.csv")
    assert rep["sequence_group"] == "1"
    assert rep["_picked_tier"] == "none"
    assert rep["_fallback_composite"] == pytest.approx(0.30)


@pytest.mark.logic
def test_fallback_is_stage_aware_on_ra_eff_but_not_on_true_jaccard(tmp_path):
    """pose_holds reads reverted ra_eff; true_jaccard is always steered.

    The row's reverted ra_eff is 1.0 and its steered ra_eff 50.0, so the
    composite is 0.80 - 0.05 x 1.0 = 0.75.  Had the steered value been used
    the composite would be -1.70.
    """
    _aggregated_csv(tmp_path, [{
        "outcome": "pose_holds",
        "steered_ra_eff_vs_truth_median": "50.0",
        "reverted_ra_eff_vs_truth_median": "1.0",
        "steered_true_jaccard_median": "0.80",
        "reverted_true_jaccard_median": "0.01",
    }])
    rep = css._pick_rep_from_aggregated(tmp_path / "passing_summary.csv")
    assert rep["_fallback_composite"] == pytest.approx(0.75)
    assert rep["true_jaccard_median"] == "0.80"


@pytest.mark.logic
def test_fallback_skips_rows_whose_components_do_not_parse(tmp_path):
    """A nan or blank component disqualifies a row from the fallback too.

    Only the third row has both components, so it wins despite the first
    two carrying a far better nominal jaccard.
    """
    _aggregated_csv(tmp_path, [
        {"outcome": "no_reversion",
         "steered_ra_eff_vs_truth_median": "nan",
         "steered_true_jaccard_median": "0.99", "sequence_group": "0"},
        {"outcome": "no_reversion",
         "steered_ra_eff_vs_truth_median": "1.0",
         "steered_true_jaccard_median": "", "sequence_group": "1"},
        {"outcome": "no_reversion",
         "steered_ra_eff_vs_truth_median": "30.0",
         "steered_true_jaccard_median": "0.10", "sequence_group": "2"},
    ])
    rep = css._pick_rep_from_aggregated(tmp_path / "passing_summary.csv")
    assert rep["sequence_group"] == "2"


@pytest.mark.logic
def test_fallback_returns_none_without_a_usable_aggregated_file(tmp_path):
    assert css._pick_rep_from_aggregated(
        tmp_path / "passing_summary.csv") is None
    _aggregated_csv(tmp_path, [{"outcome": "no_reversion"}])
    assert css._pick_rep_from_aggregated(
        tmp_path / "passing_summary.csv") is None


@pytest.mark.logic
def test_fallback_returns_none_for_a_header_only_aggregated_file(tmp_path):
    (tmp_path).mkdir(parents=True, exist_ok=True)
    (tmp_path / "aggregated_results.csv").write_text("outcome\n")
    assert css._pick_rep_from_aggregated(
        tmp_path / "passing_summary.csv") is None


# ── Sidecars and path rewriting ─────────────────────────────────────


@pytest.mark.logic
def test_rewrite_workdir_path_swaps_everything_up_to_the_sequence_segment():
    """The segment match is on a whole path component, not a substring."""
    published = Path("/out/negative_steering/runs")
    raw = "/w/b6/afc/design_7_seq_3/cycle_0/steered/p.pdb"
    assert css._rewrite_workdir_path_to_published(
        raw, "design_7_seq_3", published) == (
        "/out/negative_steering/runs/design_7_seq_3/cycle_0/steered/p.pdb")


@pytest.mark.logic
def test_rewrite_workdir_path_is_a_no_op_in_every_unsafe_case():
    """Returned unchanged with no published dir, no path, no matching
    segment, or nothing after the segment."""
    published = Path("/out/runs")
    raw = "/w/b6/design_7_seq_3/cycle_0/p.pdb"
    assert css._rewrite_workdir_path_to_published(
        raw, "design_7_seq_3", None) == raw
    assert css._rewrite_workdir_path_to_published(
        "", "design_7_seq_3", published) == ""
    assert css._rewrite_workdir_path_to_published(
        raw, "", published) == raw
    assert css._rewrite_workdir_path_to_published(
        raw, "design_9_seq_9", published) == raw
    # A name that only appears as a substring of a component is not a match.
    substring_only = "/w/b6/design_7_seq_30/cycle_0/p.pdb"
    assert css._rewrite_workdir_path_to_published(
        substring_only, "design_7_seq_3", published) == substring_only
    ends_at_segment = "/w/b6/design_7_seq_3"
    assert css._rewrite_workdir_path_to_published(
        ends_at_segment, "design_7_seq_3", published) == ends_at_segment


@pytest.mark.logic
def test_row_type_defaults_to_steered_without_a_sidecar(tmp_path):
    assert css._read_row_type_sidecar(tmp_path / "passing_summary.csv") == (
        "steered")
    (tmp_path / "row_type.txt").write_text("control_polyA\n")
    assert css._read_row_type_sidecar(tmp_path / "passing_summary.csv") == (
        "control_polyA")
    (tmp_path / "row_type.txt").write_text("   \n")
    assert css._read_row_type_sidecar(tmp_path / "passing_summary.csv") == (
        "steered")


@pytest.mark.logic
def test_runtime_prefers_the_representative_column_over_the_sidecar(tmp_path):
    """The rep column is canonical; the sidecar covers tier-none
    sequences whose passing_summary.csv has no row to read from."""
    source = tmp_path / "passing_summary.csv"
    (tmp_path / "run_one_runtime_sec.txt").write_text("999\n")
    assert css._resolve_run_one_runtime_sec(
        {"run_one_runtime_sec": "231.5"}, source) == "231.5"
    assert css._resolve_run_one_runtime_sec(
        {"run_one_runtime_sec": ""}, source) == "999"
    assert css._resolve_run_one_runtime_sec(None, source) == "999"


@pytest.mark.logic
def test_runtime_is_blank_when_the_sidecar_is_corrupt_or_absent(tmp_path):
    source = tmp_path / "passing_summary.csv"
    assert css._resolve_run_one_runtime_sec(None, source) == ""
    (tmp_path / "run_one_runtime_sec.txt").write_text("not-a-number\n")
    assert css._resolve_run_one_runtime_sec(None, source) == ""


# ── Sequence-name inference ─────────────────────────────────────────


@pytest.mark.logic
def test_infer_seq_name_from_the_parent_or_grandparent_directory(tmp_path):
    nested = tmp_path / "design_12_seq_3" / "cycle_0"
    nested.mkdir(parents=True)
    path = nested / "passing_summary.csv"
    path.touch()
    assert css._infer_seq_name(path, "parent") == "cycle_0"
    assert css._infer_seq_name(path, "grandparent") == "design_12_seq_3"
    with pytest.raises(ValueError, match="unknown seq-name inference mode"):
        css._infer_seq_name(path, "greatgrandparent")


@pytest.mark.logic
def test_infer_seq_name_rejects_a_directory_name_with_odd_characters(
        tmp_path):
    odd = tmp_path / "design 1 seq 2"
    odd.mkdir()
    path = odd / "passing_summary.csv"
    path.touch()
    with pytest.raises(ValueError, match="unexpected characters"):
        css._infer_seq_name(path, "parent")


@pytest.mark.logic
def test_parse_passing_summary_arg_accepts_name_equals_path(tmp_path):
    """An explicit NAME=PATH wins; a bare path falls back to inference.

    The split is on the first '=' only, because HPC paths can contain one.
    """
    seq_dir = tmp_path / "design_4_seq_1"
    path = _write_passing_summary(seq_dir, [])
    assert css._parse_passing_summary_arg(
        f"explicit_name={path}", "parent") == ("explicit_name", path)
    name, got = css._parse_passing_summary_arg(str(path), "parent")
    assert name == "design_4_seq_1"
    assert got == path


@pytest.mark.logic
def test_parse_passing_summary_arg_falls_back_when_the_head_is_not_a_name(
        tmp_path):
    """`NAME=` is only honoured when the path after it exists, so a real
    path containing '=' is still read as a bare path."""
    seq_dir = tmp_path / "design_4_seq_1"
    weird = seq_dir / "a=b"
    weird.mkdir(parents=True)
    path = weird / "passing_summary.csv"
    path.touch()
    name, got = css._parse_passing_summary_arg(str(path), "grandparent")
    assert name == "design_4_seq_1"
    assert got == path


# ── Input reading ───────────────────────────────────────────────────


@pytest.mark.logic
def test_read_inputs_adopts_the_first_fieldnames_and_unions_a_mismatch(
        tmp_path, capsys):
    """Version drift warns rather than aborting, and emits the superset."""
    first = _write_passing_summary(tmp_path / "s1", [_passing_row()])
    second_dir = tmp_path / "s2"
    second_dir.mkdir()
    second = second_dir / "passing_summary.csv"
    with open(second, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=[*ep.OUTPUT_FIELDS, "brand_new_column"])
        writer.writeheader()
    rows, paths, fieldnames = css._read_passing_summary_inputs(
        [("s1", first), ("s2", second)], strict=False)
    assert list(rows) == ["s1", "s2"]
    assert paths["s1"] == first
    assert fieldnames[:len(ep.OUTPUT_FIELDS)] == ep.OUTPUT_FIELDS
    assert fieldnames[-1] == "brand_new_column"
    assert "differ from the first" in capsys.readouterr().err


@pytest.mark.logic
def test_read_inputs_rejects_a_duplicate_sequence_name(tmp_path):
    path = _write_passing_summary(tmp_path / "s1", [])
    with pytest.raises(css._AggregateInputError, match="duplicate"):
        css._read_passing_summary_inputs(
            [("s1", path), ("s1", path)], strict=False)


@pytest.mark.logic
def test_read_inputs_stubs_a_missing_file_unless_strict(tmp_path, capsys):
    missing = tmp_path / "gone" / "passing_summary.csv"
    rows, paths, _ = css._read_passing_summary_inputs(
        [("s1", missing)], strict=False)
    assert rows == {"s1": []}
    assert "not found" in capsys.readouterr().err
    with pytest.raises(css._AggregateInputError, match="not found"):
        css._read_passing_summary_inputs([("s1", missing)], strict=True)


@pytest.mark.logic
def test_read_inputs_stubs_an_unreadable_file_unless_strict(tmp_path, capsys):
    """A directory where a CSV is expected is unreadable, not missing."""
    as_dir = tmp_path / "s1" / "passing_summary.csv"
    as_dir.mkdir(parents=True)
    rows, _, _ = css._read_passing_summary_inputs(
        [("s1", as_dir)], strict=False)
    assert rows == {"s1": []}
    assert "failed to read" in capsys.readouterr().err
    with pytest.raises(css._AggregateInputError, match="failed to read"):
        css._read_passing_summary_inputs([("s1", as_dir)], strict=True)


# ── Cross-sequence ranking ──────────────────────────────────────────


def _record(name: str, tier: str, ra: str, tj: str,
            row_type: str = "steered") -> dict:
    return {
        "mpnn_sequence": name,
        "row_type": row_type,
        "cross_tier": tier,
        "cross_composite_score": "",
        "cross_rank_by_composite": "",
        "cross_rank_by_ra_eff": "",
        "rep_ra_eff_vs_truth_median": ra,
        "rep_true_jaccard_median": tj,
    }


@pytest.mark.logic
def test_cross_ranks_put_tier_first_then_score():
    """Tier A beats tier B even with a worse composite inside the tier.

    Composites: a1 0.80 - 0.05 x 4.0 = 0.60, a2 0.70 - 0.05 x 2.0 = 0.60...
    so use clearly distinct values instead: a1 0.50, a2 0.70, b1 0.95.
    b1 has the best raw composite and still ranks third.
    """
    a1 = _record("a1", "A", "2.0", "0.60")   # 0.60 - 0.10 = 0.50
    a2 = _record("a2", "A", "2.0", "0.80")   # 0.80 - 0.10 = 0.70
    b1 = _record("b1", "B", "1.0", "1.00")   # 1.00 - 0.05 = 0.95
    css._assign_cross_ranks([a1, a2, b1])
    assert a2["cross_rank_by_composite"] == 1
    assert a1["cross_rank_by_composite"] == 2
    assert b1["cross_rank_by_composite"] == 3
    assert a2["cross_composite_score"] == "0.700000"


@pytest.mark.logic
def test_cross_rank_by_ra_eff_is_tier_first_then_ascending_ra_eff():
    a1 = _record("a1", "A", "4.0", "0.10")
    a2 = _record("a2", "A", "1.0", "0.10")
    c1 = _record("c1", "C", "0.1", "0.99")
    css._assign_cross_ranks([a1, a2, c1])
    assert a2["cross_rank_by_ra_eff"] == 1
    assert a1["cross_rank_by_ra_eff"] == 2
    assert c1["cross_rank_by_ra_eff"] == 3


@pytest.mark.logic
def test_tier_none_and_control_rows_are_not_ranked():
    """Controls are diagnostic, and tier-none rows are not promotable."""
    steered = _record("s", "A", "2.0", "0.80")
    none_tier = _record("n", "none", "1.0", "0.99")
    control = _record("c", "A", "1.0", "0.99", row_type="control_scrambled")
    css._assign_cross_ranks([steered, none_tier, control])
    assert steered["cross_rank_by_composite"] == 1
    assert none_tier["cross_rank_by_composite"] == ""
    assert none_tier["cross_composite_score"] == ""
    assert control["cross_rank_by_composite"] == ""


@pytest.mark.logic
def test_a_nan_metric_must_not_be_ranked_as_passing():
    """A nan ra_eff or true_jaccard keeps a tier-A row out of the ranking.

    The row still has tier A from its seed counts, but neither rank is
    assigned because `_try_float` rejects nan, so it cannot outrank a row
    with real numbers.
    """
    nan_ra = _record("nan_ra", "A", "nan", "0.80")
    nan_tj = _record("nan_tj", "A", "2.0", "nan")
    real = _record("real", "A", "9.0", "0.10")
    css._assign_cross_ranks([nan_ra, nan_tj, real])
    assert nan_ra["cross_rank_by_composite"] == ""
    assert nan_ra["cross_rank_by_ra_eff"] == ""
    assert nan_ra["cross_composite_score"] == ""
    # A nan jaccard blocks the composite rank; the ra_eff rank is still
    # available because that component parses.
    assert nan_tj["cross_rank_by_composite"] == ""
    assert nan_tj["cross_composite_score"] == ""
    assert nan_tj["cross_rank_by_ra_eff"] in (1, 2)
    assert real["cross_rank_by_composite"] == 1


@pytest.mark.logic
@pytest.mark.xfail(strict=False, reason=(
    "Known defect, reported rather than fixed: the rank keys in "
    "_assign_cross_ranks end in id(rec), so tied records are ordered by "
    "heap address instead of input order.  Fixing it would change "
    "committed cross_sequence_summary.csv output."))
def test_tied_cross_ranks_must_not_depend_on_object_identity():
    """Three exactly-tied records should rank by input order, every time.

    `_composite_rank_key` returns `(tier_order, -comp, id(rec))` and
    `_ra_rank_key` returns `(tier_order, ra, id(rec))`.  `sorted` is
    already stable, so the `id(rec)` component adds nothing but
    non-determinism: for an exact tie it orders the records by heap
    address, which depends on allocation history rather than on the input.

    Each cohort below is built with the dicts ALLOCATED in a different
    order and then assembled into the SAME input order, so a correct
    implementation must give the same rank to the same sequence all three
    times.  It does not: the ranks follow the allocation order.
    """
    def cohort(allocation_order):
        built = {name: _record(name, "A", "2.0", "0.80")
                 for name in allocation_order}
        return [built[name] for name in ("s1", "s2", "s3")]

    observed = []
    for allocation_order in [("s1", "s2", "s3"),
                             ("s3", "s2", "s1"),
                             ("s2", "s1", "s3")]:
        records = cohort(allocation_order)
        css._assign_cross_ranks(records)
        observed.append({r["mpnn_sequence"]: r["cross_rank_by_composite"]
                         for r in records})

    assert observed[0] == observed[1] == observed[2], (
        f"tied ranks varied within one process: {observed}")
    assert observed[0] == {"s1": 1, "s2": 2, "s3": 3}


# ── scored_metadata join ────────────────────────────────────────────


@pytest.mark.logic
def test_scored_metadata_join_key_is_design_and_seq(tmp_path):
    path = tmp_path / "scored_metadata.csv"
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "design", "seq", "corrected_receptor", "designed_residues",
            "native_residues"])
        writer.writeheader()
        writer.writerows([
            {"design": "12", "seq": "3", "corrected_receptor": "MKV",
             "designed_residues": "AB|CD", "native_residues": "EF|GH"},
            {"design": "", "seq": "3", "corrected_receptor": "skipped"},
            {"design": "9", "seq": "", "corrected_receptor": "skipped"},
        ])
    lookup = css._load_scored_metadata(path)
    assert list(lookup) == ["design_12_seq_3"]
    assert lookup["design_12_seq_3"]["designed_residues"] == "AB|CD"


@pytest.mark.logic
def test_scored_metadata_is_an_empty_lookup_when_not_supplied(
        tmp_path, capsys):
    assert css._load_scored_metadata(None) == {}
    assert css._load_scored_metadata(tmp_path / "absent.csv") == {}
    assert "not found" in capsys.readouterr().err


# ── aggregate() ─────────────────────────────────────────────────────


@pytest.mark.logic
def test_aggregate_writes_one_row_per_sequence_with_the_meta_columns(
        tmp_path):
    """Three sequences, three tiers, ranked A before B before C.

    Tier counts are per-sequence counts over that file's rows, not over the
    representative alone, so sequence s1 reports two tier-A rows.
    """
    s1 = _write_passing_summary(tmp_path / "s1", [
        _passing_row(n_pass="3", n_seeds="3", ra_eff_vs_truth_median="2.0",
                     rank_by_composite_score="1"),
        _passing_row(n_pass="3", n_seeds="3", ra_eff_vs_truth_median="3.0",
                     rank_by_composite_score="2"),
    ])
    s2 = _write_passing_summary(tmp_path / "s2", [
        _passing_row(n_pass="2", n_seeds="3", ra_eff_vs_truth_median="0.5"),
    ])
    s3 = _write_passing_summary(tmp_path / "s3", [
        _passing_row(n_pass="1", n_seeds="3", ra_eff_vs_truth_median="0.1"),
    ])
    out = tmp_path / "cross.csv"
    rc = css.aggregate([("s1", s1), ("s2", s2), ("s3", s3)], out,
                       strict=False)
    assert rc == 0

    got = _read_output(out)
    assert [r["mpnn_sequence"] for r in got] == ["s1", "s2", "s3"]
    assert [r["cross_tier"] for r in got] == ["A", "B", "C"]
    assert [r["cross_rank_by_composite"] for r in got] == ["1", "2", "3"]
    assert got[0]["passing_rows_in_sequence"] == "2"
    assert got[0]["n_tier_A_in_sequence"] == "2"
    assert got[1]["n_tier_B_in_sequence"] == "1"
    assert got[2]["n_tier_C_in_sequence"] == "1"
    # The representative is the better of s1's two tier-A rows.
    assert got[0]["rep_ra_eff_vs_truth_median"] == "2.0"
    assert got[0]["within_sequence_rank"] == "1"
    assert got[0]["row_type"] == "steered"


@pytest.mark.logic
def test_aggregate_column_order_is_meta_then_prefixed_rep_columns(tmp_path):
    s1 = _write_passing_summary(tmp_path / "s1", [_passing_row()])
    out = tmp_path / "cross.csv"
    css.aggregate([("s1", s1)], out, strict=False)
    with open(out) as f:
        header = next(csv.reader(f))
    assert header[:13] == [
        "mpnn_sequence", "row_type", "cross_rank_by_composite",
        "cross_rank_by_ra_eff", "cross_tier", "cross_composite_score",
        "within_sequence_rank", "passing_rows_in_sequence",
        "n_tier_A_in_sequence", "n_tier_B_in_sequence",
        "n_tier_C_in_sequence", "run_one_runtime_sec",
        "source_passing_summary"]
    assert header[13:16] == ["corrected_receptor", "designed_residues",
                             "native_residues"]
    assert header[16:] == [f"rep_{c}" for c in ep.OUTPUT_FIELDS]


@pytest.mark.logic
def test_aggregate_gives_an_empty_sequence_a_tier_none_stub(tmp_path):
    """An empty passing_summary.csv with no aggregated sibling still gets a
    row, with blank rep_* columns and both cross-rank columns blank."""
    empty = _write_passing_summary(tmp_path / "s_empty", [])
    out = tmp_path / "cross.csv"
    css.aggregate([("s_empty", empty)], out, strict=False)
    got = _read_output(out)
    assert len(got) == 1
    assert got[0]["cross_tier"] == "none"
    assert got[0]["cross_rank_by_composite"] == ""
    assert got[0]["cross_rank_by_ra_eff"] == ""
    assert got[0]["cross_composite_score"] == ""
    assert got[0]["within_sequence_rank"] == ""
    assert got[0]["rep_ra_eff_vs_truth_median"] == ""


@pytest.mark.logic
def test_aggregate_fills_a_tier_none_row_from_aggregated_results(tmp_path):
    """The fallback stamps real metrics and a composite onto a tier-none
    row, so the cohort plot has something to show, while both rank columns
    stay blank."""
    seq_dir = tmp_path / "s_none"
    empty = _write_passing_summary(seq_dir, [])
    _aggregated_csv(seq_dir, [{
        "outcome": "no_reversion",
        "n_seeds": "3", "n_pass": "0",
        "steered_ra_eff_vs_truth_median": "12.0",
        "steered_true_jaccard_median": "0.40",
        "steered_iptm_median": "0.85",
    }])
    out = tmp_path / "cross.csv"
    css.aggregate([("s_none", empty)], out, strict=False)
    got = _read_output(out)[0]
    assert got["cross_tier"] == "none"
    # 0.40 - 0.05 x 12.0 = -0.20
    assert got["cross_composite_score"] == "-0.200000"
    assert got["cross_rank_by_composite"] == ""
    assert got["rep_ra_eff_vs_truth_median"] == "12.0"
    assert got["rep_iptm_median"] == "0.85"


@pytest.mark.logic
def test_aggregate_rewrites_the_canonical_pdb_path_when_asked(tmp_path):
    workdir_path = "/w/ab/cd/s1/cycle_0/steered/pred.pdb"
    s1 = _write_passing_summary(tmp_path / "s1", [
        _passing_row(canonical_pdb=workdir_path)])
    out = tmp_path / "cross.csv"
    published = tmp_path / "published" / "runs"
    css.aggregate([("s1", s1)], out, strict=False,
                  published_runs_dir=published)
    got = _read_output(out)[0]
    assert got["rep_canonical_pdb"] == str(
        published / "s1" / "cycle_0" / "steered" / "pred.pdb")

    css.aggregate([("s1", s1)], out, strict=False)
    assert _read_output(out)[0]["rep_canonical_pdb"] == workdir_path


@pytest.mark.logic
def test_aggregate_excludes_controls_from_the_ranking_but_keeps_the_row(
        tmp_path):
    """A control row is emitted with every cross_rank column blank.

    main.nf's control-miscalibration check reads them back from the final
    CSV, so they must be present and identifiable by row_type.
    """
    steered = _write_passing_summary(tmp_path / "s1", [
        _passing_row(ra_eff_vs_truth_median="5.0")])
    control_dir = tmp_path / "input_control_polyA"
    control = _write_passing_summary(control_dir, [
        _passing_row(ra_eff_vs_truth_median="0.1")])
    (control_dir / "row_type.txt").write_text("control_polyA\n")

    out = tmp_path / "cross.csv"
    css.aggregate([("s1", steered), ("input_control_polyA", control)],
                  out, strict=False)
    by_name = {r["mpnn_sequence"]: r for r in _read_output(out)}
    assert by_name["s1"]["cross_rank_by_composite"] == "1"
    assert by_name["input_control_polyA"]["row_type"] == "control_polyA"
    assert by_name["input_control_polyA"]["cross_rank_by_composite"] == ""
    assert by_name["input_control_polyA"]["cross_rank_by_ra_eff"] == ""
    # The control still carries its representative's metrics.
    assert by_name["input_control_polyA"]["rep_ra_eff_vs_truth_median"] == (
        "0.1")


@pytest.mark.logic
def test_aggregate_joins_scored_metadata_by_mpnn_sequence(tmp_path):
    s1 = _write_passing_summary(tmp_path / "design_12_seq_3",
                                [_passing_row()])
    s2 = _write_passing_summary(tmp_path / "design_99_seq_9",
                                [_passing_row()])
    meta = tmp_path / "scored_metadata.csv"
    with open(meta, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "design", "seq", "corrected_receptor", "designed_residues",
            "native_residues"])
        writer.writeheader()
        writer.writerow({"design": "12", "seq": "3",
                         "corrected_receptor": "MKV",
                         "designed_residues": "AB|CD",
                         "native_residues": "EF|GH"})
    out = tmp_path / "cross.csv"
    css.aggregate([("design_12_seq_3", s1), ("design_99_seq_9", s2)],
                  out, strict=False, scored_metadata_path=meta)
    by_name = {r["mpnn_sequence"]: r for r in _read_output(out)}
    assert by_name["design_12_seq_3"]["corrected_receptor"] == "MKV"
    assert by_name["design_12_seq_3"]["native_residues"] == "EF|GH"
    # An unmatched sequence gets blanks, not a missing column.
    assert by_name["design_99_seq_9"]["corrected_receptor"] == ""
    assert by_name["design_99_seq_9"]["designed_residues"] == ""


@pytest.mark.logic
def test_aggregate_reads_the_runtime_from_the_rep_row_or_the_sidecar(
        tmp_path):
    with_column = _write_passing_summary(tmp_path / "s1", [
        _passing_row(run_one_runtime_sec="100.0")])
    sidecar_dir = tmp_path / "s2"
    sidecar_only = _write_passing_summary(sidecar_dir, [_passing_row()])
    (sidecar_dir / "run_one_runtime_sec.txt").write_text("200.0\n")
    out = tmp_path / "cross.csv"
    css.aggregate([("s1", with_column), ("s2", sidecar_only)], out,
                  strict=False)
    by_name = {r["mpnn_sequence"]: r for r in _read_output(out)}
    assert by_name["s1"]["run_one_runtime_sec"] == "100.0"
    assert by_name["s2"]["run_one_runtime_sec"] == "200.0"


@pytest.mark.logic
def test_aggregate_reports_two_with_no_inputs(capsys):
    assert css.aggregate([], Path("unused.csv"), strict=False) == 2
    assert "no inputs supplied" in capsys.readouterr().err


@pytest.mark.logic
def test_aggregate_reports_two_on_a_duplicate_sequence_name(tmp_path, capsys):
    path = _write_passing_summary(tmp_path / "s1", [_passing_row()])
    rc = css.aggregate([("s1", path), ("s1", path)], tmp_path / "out.csv",
                       strict=False)
    assert rc == 2
    assert "duplicate sequence name" in capsys.readouterr().err


@pytest.mark.logic
def test_aggregate_prints_the_tier_and_row_type_breakdown(tmp_path, capsys):
    steered = _write_passing_summary(tmp_path / "s1", [_passing_row()])
    control_dir = tmp_path / "ctrl"
    control = _write_passing_summary(control_dir, [_passing_row()])
    (control_dir / "row_type.txt").write_text("control_scrambled\n")
    out = tmp_path / "cross.csv"
    css.aggregate([("s1", steered), ("ctrl", control)], out, strict=False)
    err = capsys.readouterr().err
    assert "2 sequences, 1 placed in tiers A/B/C" in err
    assert "A=2 B=0 C=0 none=0" in err
    assert "control_scrambled=1" in err
    assert "steered=1" in err


# ── CLI ─────────────────────────────────────────────────────────────


def _run_main(monkeypatch, *argv_tail) -> int:
    monkeypatch.setattr(sys, "argv",
                        ["cross_sequence_summary.py", *argv_tail])
    return css.main()


@pytest.mark.logic
def test_main_infers_names_from_passing_summary_dir(tmp_path, monkeypatch):
    """--passing-summary-dir enumerates one subdirectory per sequence."""
    runs = tmp_path / "runs"
    _write_passing_summary(runs / "design_1_seq_0", [_passing_row()])
    _write_passing_summary(runs / "design_2_seq_1", [
        _passing_row(n_pass="1", n_seeds="3")])
    out = tmp_path / "cross.csv"
    assert _run_main(monkeypatch, "--passing-summary-dir", str(runs),
                     "--output", str(out)) == 0
    got = _read_output(out)
    assert sorted(r["mpnn_sequence"] for r in got) == [
        "design_1_seq_0", "design_2_seq_1"]


@pytest.mark.logic
def test_main_skips_a_subdirectory_whose_name_is_not_a_sequence_id(
        tmp_path, monkeypatch, capsys):
    runs = tmp_path / "runs"
    _write_passing_summary(runs / "design_1_seq_0", [_passing_row()])
    _write_passing_summary(runs / "not a name", [_passing_row()])
    (runs / "loose_file.txt").write_text("ignored\n")
    out = tmp_path / "cross.csv"
    assert _run_main(monkeypatch, "--passing-summary-dir", str(runs),
                     "--output", str(out)) == 0
    assert [r["mpnn_sequence"] for r in _read_output(out)] == [
        "design_1_seq_0"]
    assert "not a valid sequence identifier" in capsys.readouterr().err


@pytest.mark.logic
def test_main_accepts_explicit_name_equals_path_pairs(tmp_path, monkeypatch):
    path = _write_passing_summary(tmp_path / "whatever", [_passing_row()])
    out = tmp_path / "cross.csv"
    assert _run_main(monkeypatch, "--passing-summary",
                     f"design_7_seq_2={path}", "--output", str(out)) == 0
    assert [r["mpnn_sequence"] for r in _read_output(out)] == [
        "design_7_seq_2"]


@pytest.mark.logic
def test_main_uses_the_grandparent_name_when_asked(tmp_path, monkeypatch):
    path = _write_passing_summary(
        tmp_path / "design_5_seq_5" / "cycle_0", [_passing_row()])
    out = tmp_path / "cross.csv"
    assert _run_main(monkeypatch, "--passing-summary", str(path),
                     "--seq-name-from-grandparent-dir",
                     "--output", str(out)) == 0
    assert [r["mpnn_sequence"] for r in _read_output(out)] == [
        "design_5_seq_5"]


@pytest.mark.logic
def test_main_does_not_re_add_a_sequence_already_given_explicitly(
        tmp_path, monkeypatch):
    runs = tmp_path / "runs"
    path = _write_passing_summary(runs / "design_1_seq_0", [_passing_row()])
    out = tmp_path / "cross.csv"
    assert _run_main(monkeypatch, "--passing-summary", str(path),
                     "--passing-summary-dir", str(runs),
                     "--output", str(out)) == 0
    assert len(_read_output(out)) == 1


@pytest.mark.logic
def test_main_resolves_published_runs_dir_to_an_absolute_path(
        tmp_path, monkeypatch):
    """A relative --published-runs-dir would resolve against the wrong CWD
    in a downstream container task, so it is made absolute first."""
    runs = tmp_path / "runs"
    _write_passing_summary(runs / "design_1_seq_0", [
        _passing_row(canonical_pdb="/w/ab/design_1_seq_0/cycle_0/p.pdb")])
    monkeypatch.chdir(tmp_path)
    out = tmp_path / "cross.csv"
    assert _run_main(monkeypatch, "--passing-summary-dir", str(runs),
                     "--published-runs-dir", "published/runs",
                     "--output", str(out)) == 0
    got = _read_output(out)[0]["rep_canonical_pdb"]
    assert Path(got).is_absolute()
    assert got.endswith("published/runs/design_1_seq_0/cycle_0/p.pdb")


@pytest.mark.logic
def test_main_reports_two_on_a_duplicate_explicit_name(tmp_path, monkeypatch,
                                                       capsys):
    path = _write_passing_summary(tmp_path / "s1", [_passing_row()])
    rc = _run_main(monkeypatch, "--passing-summary", f"dup={path}",
                   "--passing-summary", f"dup={path}",
                   "--output", str(tmp_path / "out.csv"))
    assert rc == 2
    assert "duplicate sequence name" in capsys.readouterr().err


@pytest.mark.logic
def test_main_reports_two_on_an_unusable_name_or_directory(
        tmp_path, monkeypatch, capsys):
    odd = tmp_path / "bad name"
    path = _write_passing_summary(odd, [_passing_row()])
    rc = _run_main(monkeypatch, "--passing-summary", str(path),
                   "--output", str(tmp_path / "out.csv"))
    assert rc == 2
    assert "unexpected characters" in capsys.readouterr().err

    rc = _run_main(monkeypatch, "--passing-summary-dir",
                   str(tmp_path / "absent"),
                   "--output", str(tmp_path / "out.csv"))
    assert rc == 2
    assert "not a directory" in capsys.readouterr().err


@pytest.mark.logic
def test_main_reports_two_when_a_directory_holds_no_sequences(
        tmp_path, monkeypatch, capsys):
    empty = tmp_path / "runs"
    empty.mkdir()
    rc = _run_main(monkeypatch, "--passing-summary-dir", str(empty),
                   "--output", str(tmp_path / "out.csv"))
    assert rc == 2
    assert "no inputs found" in capsys.readouterr().err


@pytest.mark.logic
def test_main_strict_aborts_on_a_missing_input(tmp_path, monkeypatch, capsys):
    runs = tmp_path / "runs"
    (runs / "design_1_seq_0").mkdir(parents=True)
    rc = _run_main(monkeypatch, "--passing-summary-dir", str(runs),
                   "--strict", "--output", str(tmp_path / "out.csv"))
    assert rc == 2
    assert "not found" in capsys.readouterr().err
