"""Stage 5 (cohort figures) — bin/negsteer_plots.py.

Producer: ``NEGSTEER_PLOTS`` → ``bin/negsteer_plots.py``, which reads
cross_sequence_summary.csv plus the per-sequence runs/ tree and writes the
seven cohort PNGs that Nextflow publishes.

What these tests hold in place, beyond "it ran":

* Every ``plot_*`` function returns True when it drew and False when it had no
  data, and on the False path it still writes a placeholder PNG. The process
  publishes a fixed figure list, so a missing path fails the stage with nothing
  in the log to explain why.
* The PNG is not blank. A caught exception leaves a file behind too.
* ``classify_seed`` sorts a per-seed row into its stage and outcome, and its
  copy in bin/negsteer_within_sequence_plots.py agrees with it for every
  outcome. The two modules carry separate copies; a drift between them would
  make the cohort figures and the within-sequence figures disagree about the
  same seed.
* Which prediction is final. A sequence that went through reversion has both a
  steered and a reverted prediction, and the reverted one is the last made.
  Reading the steered value unconditionally is what made a collapsed pose look
  like a tier-A hit.
* A nan metric is absent data, not a failed check. A prediction that never
  landed writes nan, and nan fails every comparison it is put through.

Fixture cohort: tests/characterization/fixtures_cohort.py, which states each
figure's content at the call site and writes the CSV rows that produce it.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

from tests.characterization import fixtures_cohort as fc
from tests.characterization.helpers.png_compare import count_distinct_colours

_BIN = Path(__file__).resolve().parents[2] / "bin"
if str(_BIN) not in sys.path:
    sys.path.insert(0, str(_BIN))

# Both modules select the Agg backend themselves, before importing pyplot, so
# importing them here is what makes this file headless. Nothing above may
# import pyplot.
import negsteer_plots as cohort  # noqa: E402
import negsteer_within_sequence_plots as within  # noqa: E402

THRESHOLDS = {
    "complex_plddt_min": 0.70,
    "ipae_max": 15.0,
    "pae_pass_frac_min": 0.10,
    "iptm_min": 0.30,
    "control_ra_eff_min": 5.0,
}


def assert_drawn(out: Path) -> None:
    """The figure landed and has more than one colour in it."""
    assert out.exists(), f"{out.name} was not written"
    assert out.stat().st_size > 0, f"{out.name} is empty"
    n = count_distinct_colours(out)
    assert n > 1, f"{out.name} is a single flat colour ({n} distinct)"


def steered_rows(n: int = 3) -> list[dict]:
    """`n` steered cohort rows, tier A first, every metric present."""
    rows = []
    for i in range(n):
        rows.append({
            "mpnn_sequence": f"design_{i}_seq_0",
            "row_type": "steered",
            "cross_tier": ["A", "B", "C", "none"][min(i, 3)],
            "cross_rank_by_composite": str(i + 1),
            "cross_composite_score": f"{0.6 - 0.1 * i:.3f}",
            "rep_ra_eff_vs_truth_median": "2.0",
            "rep_true_jaccard_median": f"{0.7 - 0.1 * i:.2f}",
            "rep_complex_plddt_median": "0.92",
            "rep_ipae_median": "4.0",
            "rep_pae_pass_frac_median": "0.80",
            "rep_iptm_median": "0.72",
            "rep_outcome": "no_reversion",
            "rep_steered_mutations_chimerax": f"/A:{5 + i}",
        })
    return rows


# ═══════════════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════════════

@pytest.mark.logic
@pytest.mark.parametrize("raw, expected", [
    ("2.5", 2.5),
    (" 2.5 ", 2.5),
    (3, 3.0),
    ("", None),
    ("   ", None),
    (None, None),
    ("not_a_number", None),
    ("nan", None),
    ("NaN", None),
    ("inf", None),
    ("-inf", None),
])
def test_try_float_treats_nan_and_inf_as_absent(raw, expected):
    """nan is what a prediction that never landed writes. Parsed as a float it
    would survive into every plot and then fail every comparison silently, so
    both modules must hand it back as absent."""
    assert cohort._try_float(raw) == expected
    assert within._try_float(raw) == expected


@pytest.mark.logic
@pytest.mark.parametrize("raw, expected", [
    ("3", 3),
    ("3.0", 3),
    ("3.7", 3),
    ("", None),
    (None, None),
    ("x", None),
])
def test_try_int_accepts_the_float_strings_the_engine_writes(raw, expected):
    assert cohort._try_int(raw) == expected
    assert within._try_int(raw) == expected


@pytest.mark.logic
@pytest.mark.parametrize("raw, expected", [
    ("/A:5,7,22", [5, 7, 22]),
    ("/A:5", [5]),
    ("/A: 5 , 7 ", [5, 7]),
    ("/A:5,x,7", [5, 7]),
    ("", []),
    ("A5,A7", []),
])
def test_parse_chimerax_positions(raw, expected):
    assert cohort._parse_chimerax_positions(raw) == expected


@pytest.mark.logic
@pytest.mark.parametrize("raw, expected", [
    ("design_0_seq_1", "d0_s1"),
    ("design_27_seq_0", "d27_s0"),
    ("input_control_polyA", "ctrl_poly"),
    ("input_control_scrambled", "ctrl_scra"),
    ("something_else", "something_else"),
])
def test_short_name(raw, expected):
    assert cohort._short_name(raw) == expected
    assert within._short_name(raw) == expected


@pytest.mark.logic
def test_row_type_defaults_to_steered():
    assert cohort._row_type({}) == "steered"
    assert cohort._row_type({"row_type": ""}) == "steered"
    assert cohort._row_type({"row_type": " control_polyA "}) == "control_polyA"
    assert cohort._is_steered({"row_type": "steered"}) is True
    assert cohort._is_steered({"row_type": "control_polyA"}) is False


@pytest.mark.logic
def test_composite_is_jaccard_minus_a_twentieth_of_the_rmsd():
    row = {"rep_ra_eff_vs_truth_median": "4.0", "rep_true_jaccard_median": "0.60"}
    assert cohort._composite_from_row(row) == pytest.approx(0.60 - 0.05 * 4.0)


@pytest.mark.logic
@pytest.mark.parametrize("row", [
    {"rep_true_jaccard_median": "0.60"},
    {"rep_ra_eff_vs_truth_median": "4.0"},
    {"rep_ra_eff_vs_truth_median": "nan", "rep_true_jaccard_median": "0.60"},
])
def test_composite_is_none_when_a_component_is_missing(row):
    assert cohort._composite_from_row(row) is None


# ═══════════════════════════════════════════════════════════════════════════
# classify_seed — and the agreement between the two modules' copies
# ═══════════════════════════════════════════════════════════════════════════

_CLASSIFY_CASES = [
    ("cold_start", "pass"),
    ("cold_start", "off_target"),
    ("cold_start", "poor_prediction"),
    ("cold_start", "multiple_failures"),
    ("cold_start", "no_data"),
    ("steering", "pass"),
    ("steering", "off_target"),
    ("steering", "poor_prediction"),
    ("steering", "multiple_failures"),
    ("steering", "no_data"),
    ("reversion", "pass"),
    ("reversion", "off_target"),
    ("reversion", "poor_prediction"),
    ("reversion", "multiple_failures"),
    ("reversion", "no_data"),
    ("reversion", "new_contamination"),
    ("no_data", "no_data"),
]


@pytest.mark.logic
@pytest.mark.parametrize("stage, outcome", _CLASSIFY_CASES,
                          ids=[f"{s}-{o}" for s, o in _CLASSIFY_CASES])
def test_classify_seed_agrees_across_modules(stage, outcome):
    """Both modules carry their own copy of classify_seed. Each copy must put
    the seed in the stated stage and outcome, and both must say the same thing
    about the same row — a drift would make the cohort figures and the
    within-sequence figures disagree about one seed."""
    row = fc.seed_row(stage=stage, outcome=outcome)
    assert cohort.classify_seed(row) == (stage, outcome)
    assert within.classify_seed(row) == cohort.classify_seed(row)


@pytest.mark.logic
def test_classify_seed_reads_the_reverted_prediction_not_the_steered_one():
    """A reverted seed has both predictions. The reverted one is the last made,
    so a seed whose steering looked healthy and whose reversion collapsed is an
    off-target seed, not a passing one."""
    row = fc.seed_row(stage="reversion", outcome="pass",
                      reverted_ra_eff_vs_truth=30.5)
    assert row["steered_ra_eff_vs_truth"] == "2.5"
    assert cohort.classify_seed(row) == ("reversion", "off_target")
    assert within.classify_seed(row) == ("reversion", "off_target")


@pytest.mark.logic
@pytest.mark.parametrize("missing", [
    "steered_ra_eff_vs_truth",
    "steered_receptor_intact",
    "steered_complex_plddt",
    "steered_ipae",
    "steered_pae_pass_frac",
    "steered_iptm",
])
def test_classify_seed_calls_a_nan_metric_no_data_not_a_failure(missing):
    """nan fails every comparison, so a prediction that never landed would
    otherwise be reported as a structural or confidence failure that the run
    never actually measured."""
    row = fc.seed_row(stage="steering", outcome="pass", **{missing: "nan"})
    assert cohort.classify_seed(row) == ("steering", "no_data")
    assert within.classify_seed(row) == ("steering", "no_data")


@pytest.mark.logic
def test_classify_seed_keeps_the_stage_when_the_outcome_is_unknown():
    """An unreadable metric loses the outcome but not the stage: the seed still
    reached reversion, which is what the heatmap's marker shape shows."""
    row = fc.seed_row(stage="reversion", outcome="no_data")
    stage, outcome = cohort.classify_seed(row)
    assert (stage, outcome) == ("reversion", "no_data")
    assert stage in cohort.STAGE_MARKER


# ═══════════════════════════════════════════════════════════════════════════
# Which prediction is final
# ═══════════════════════════════════════════════════════════════════════════

@pytest.mark.logic
@pytest.mark.parametrize("outcome, expected_ra", [
    ("pose_holds", "30.5"),
    ("pose_collapses", "30.5"),
    ("new_contamination", "30.5"),
    ("no_reversion", "2.0"),
    ("singleton", "2.0"),
    ("", "2.0"),
])
def test_final_prediction_follows_the_outcome(outcome, expected_ra):
    """The three outcomes that mean reversion ran take the reverted metric. The
    rest never produced a reverted prediction, so the steered column holds the
    only prediction there is."""
    row = fc.agg_row_for(outcome, steered_ra=2.0, steered_tj=0.70,
                         reverted_ra=30.5, reverted_tj=0.05)
    got = cohort._final_prediction_value(row, "rep_ra_eff_vs_truth_median")
    assert got == expected_ra


@pytest.mark.logic
def test_final_prediction_falls_back_when_reversion_made_no_prediction():
    """Reversion can run and fail to produce a prediction at all. The steered
    value is then the only number available."""
    row = fc.agg_row_for("pose_collapses", steered_ra=2.0, steered_tj=0.70)
    assert row["reverted_ra_eff_vs_truth_median"] == ""
    got = cohort._final_prediction_value(row, "rep_ra_eff_vs_truth_median")
    assert got == "2.0"


@pytest.mark.logic
def test_final_prediction_of_a_field_with_no_reverted_twin_is_the_steered_one():
    row = fc.agg_row(n_pass="3", outcome="pose_collapses")
    assert cohort._final_prediction_value(row, "rep_n_pass") == "3"
    assert cohort._final_prediction_value(row, "not_a_field") == ""


@pytest.mark.logic
def test_project_agg_row_gives_a_tierless_cross_shaped_row():
    row = cohort._project_agg_row(
        "design_9_seq_0",
        fc.agg_row_for("no_reversion", steered_ra=3.0, steered_tj=0.40),
        row_type="control_polyA",
    )
    assert row["mpnn_sequence"] == "design_9_seq_0"
    assert row["row_type"] == "control_polyA"
    assert row["cross_tier"] == "none"
    assert row["cross_rank_by_composite"] == ""
    assert row["rep_ra_eff_vs_truth_median"] == "3.0"


# ═══════════════════════════════════════════════════════════════════════════
# _pick_best_agg_row
# ═══════════════════════════════════════════════════════════════════════════

@pytest.mark.logic
def test_pick_best_agg_row_returns_none_for_no_rows():
    assert cohort._pick_best_agg_row([]) is None


@pytest.mark.logic
def test_pick_best_agg_row_honours_the_cross_summary_pointer():
    """A tiered sequence names its representative row, and that is the row the
    cross-rank was computed on — not the best-looking one."""
    rows = [
        fc.agg_row_for("no_reversion", steered_ra=1.0, steered_tj=0.90,
                        sequence_group="0"),
        fc.agg_row_for("no_reversion", steered_ra=9.0, steered_tj=0.10,
                        sequence_group="2"),
    ]
    got = cohort._pick_best_agg_row(rows, prefer_cycle="0",
                                    prefer_pathway="cycle_0",
                                    prefer_sequence_group="2")
    assert got["sequence_group"] == "2"


@pytest.mark.logic
def test_pick_best_agg_row_falls_back_to_the_best_composite():
    """Tier-none sequences leave the pointer blank. The best attempt is the
    most informative single summary of a sequence nothing passed on."""
    rows = [
        fc.agg_row_for("no_reversion", steered_ra=9.0, steered_tj=0.10,
                        sequence_group="0"),
        fc.agg_row_for("no_reversion", steered_ra=1.0, steered_tj=0.90,
                        sequence_group="1"),
    ]
    assert cohort._pick_best_agg_row(rows)["sequence_group"] == "1"


@pytest.mark.logic
def test_pick_best_agg_row_skips_singletons_unless_they_are_all_there_is():
    """Singleton rows carry no real seed counts."""
    singleton = fc.agg_row_for("singleton", steered_ra=0.1, steered_tj=0.99,
                                sequence_group="9")
    real = fc.agg_row_for("no_reversion", steered_ra=9.0, steered_tj=0.10,
                           sequence_group="0")
    assert cohort._pick_best_agg_row([singleton, real])["sequence_group"] == "0"
    assert cohort._pick_best_agg_row([singleton])["sequence_group"] == "9"


@pytest.mark.logic
def test_pick_best_agg_row_prefers_a_row_with_numbers_over_one_without():
    blank = fc.agg_row(sequence_group="0")
    real = fc.agg_row_for("no_reversion", steered_ra=40.0, steered_tj=0.0,
                           sequence_group="1")
    assert cohort._pick_best_agg_row([blank, real])["sequence_group"] == "1"


# ═══════════════════════════════════════════════════════════════════════════
# load_unified_cohort
# ═══════════════════════════════════════════════════════════════════════════

@pytest.mark.logic
def test_load_unified_cohort_resources_metrics_from_the_final_prediction(tmp_path):
    """design_0_seq_1 collapsed on reversion. cross_summary's own rep_* columns
    describe its steered intermediate (2.0 Å, Jaccard 0.70); the merge must
    replace them with the reverted prediction (30.5 Å, Jaccard 0.05)."""
    paths = fc.cohort_tree(tmp_path)
    rows = cohort.load_unified_cohort(paths["cross_csv"], paths["runs_dir"],
                                      verbose=True)
    by_name = {r["mpnn_sequence"]: r for r in rows}

    collapsed = by_name["design_0_seq_1"]
    assert collapsed["rep_outcome"] == "pose_collapses"
    assert float(collapsed["rep_ra_eff_vs_truth_median"]) == 30.5
    assert float(collapsed["rep_true_jaccard_median"]) == 0.05
    assert cohort._composite_from_row(collapsed) == pytest.approx(0.05 - 0.05 * 30.5)

    held = by_name["design_0_seq_0"]
    assert float(held["rep_ra_eff_vs_truth_median"]) == 2.0
    assert float(held["rep_true_jaccard_median"]) == 0.70


@pytest.mark.logic
def test_load_unified_cohort_fills_blank_seed_counts_but_not_metrics(tmp_path):
    """Seed counts and mutation strings have no reverted twin, so they are
    filled only where cross_summary left them blank."""
    paths = fc.cohort_tree(tmp_path)
    rows = cohort.load_unified_cohort(paths["cross_csv"], paths["runs_dir"],
                                      verbose=False)
    by_name = {r["mpnn_sequence"]: r for r in rows}
    assert by_name["design_0_seq_0"]["rep_n_pass"] == "3"
    assert by_name["design_1_seq_0"]["rep_n_seeds"] == "3"


@pytest.mark.logic
def test_load_unified_cohort_synthesises_sequences_absent_from_the_summary(tmp_path):
    """An older workdir or a partial run leaves sequence dirs the summary never
    saw. They still belong in the cohort, with the row type read off the
    directory name so controls stay labelled."""
    paths = fc.cohort_tree(tmp_path)
    rows = cohort.load_unified_cohort(paths["cross_csv"], paths["runs_dir"],
                                      verbose=False)
    by_name = {r["mpnn_sequence"]: r for r in rows}
    assert "design_2_seq_0" in by_name
    assert by_name["design_2_seq_0"]["row_type"] == "steered"
    assert by_name["design_2_seq_0"]["cross_tier"] == "none"


@pytest.mark.logic
@pytest.mark.parametrize("name, expected", [
    ("input_control_polyA_x", "control_polyA"),
    ("input_control_scrambled_x", "control_scrambled"),
    ("design_9_seq_9", "steered"),
])
def test_load_unified_cohort_labels_a_synthesised_control(tmp_path, name, expected):
    cross = tmp_path / "cross_sequence_summary.csv"
    fc.write_cross_summary(cross, [fc.cross_row("design_0_seq_0")])
    runs = tmp_path / "runs"
    fc.write_sequence_dir(runs, name, agg_rows=[
        fc.agg_row_for("no_reversion", steered_ra=5.0, steered_tj=0.2)])
    rows = cohort.load_unified_cohort(cross, runs, verbose=False)
    by_name = {r["mpnn_sequence"]: r for r in rows}
    assert by_name[name]["row_type"] == expected


@pytest.mark.logic
def test_load_unified_cohort_without_a_runs_dir_returns_the_summary_verbatim(tmp_path):
    paths = fc.cohort_tree(tmp_path)
    rows = cohort.load_unified_cohort(paths["cross_csv"], None, verbose=True)
    assert len(rows) == 5
    collapsed = {r["mpnn_sequence"]: r for r in rows}["design_0_seq_1"]
    # No runs/ means no backfill, so the steered intermediate is all there is.
    assert collapsed["rep_ra_eff_vs_truth_median"] == "2.0"

    missing = cohort.load_unified_cohort(paths["cross_csv"],
                                         tmp_path / "absent", verbose=False)
    assert len(missing) == 5


@pytest.mark.logic
def test_load_unified_cohort_passes_through_a_sequence_with_no_runs_dir(tmp_path):
    cross = tmp_path / "cross_sequence_summary.csv"
    fc.write_cross_summary(cross, [
        fc.cross_row("design_0_seq_0", rep_ra_eff_vs_truth_median=2.0,
                      rep_true_jaccard_median=0.7),
        fc.cross_row("design_7_seq_7", rep_ra_eff_vs_truth_median=3.0,
                      rep_true_jaccard_median=0.3),
    ])
    runs = tmp_path / "runs"
    fc.write_sequence_dir(runs, "design_0_seq_0", agg_rows=[
        fc.agg_row_for("no_reversion", steered_ra=2.0, steered_tj=0.70)])
    # A loose file beside the sequence dirs, and a dir with no
    # aggregated_results.csv, are both skipped rather than fatal.
    (runs / "stray.txt").write_text("not a sequence dir\n")
    (runs / "design_8_seq_8").mkdir()

    rows = cohort.load_unified_cohort(cross, runs, verbose=True)
    names = {r["mpnn_sequence"] for r in rows}
    assert names == {"design_0_seq_0", "design_7_seq_7"}


# ═══════════════════════════════════════════════════════════════════════════
# Seed outcomes
# ═══════════════════════════════════════════════════════════════════════════

@pytest.mark.logic
def test_load_seed_outcomes_keeps_the_representative_group_by_default(tmp_path):
    """design_0_seq_0 has three passing seeds in sequence_group 0 and three
    off-target seeds in group 1, and cross_summary names group 0. The heatmap
    sees three cells, the bars see all six."""
    paths = fc.cohort_tree(tmp_path)
    rep, comps = cohort.load_seed_outcomes(paths["runs_dir"],
                                           paths["cross_csv"], all_sgs=False)
    all_sg, _ = cohort.load_seed_outcomes(paths["runs_dir"],
                                          paths["cross_csv"], all_sgs=True)

    assert [o for _s, o, _i in rep["design_0_seq_0"]] == ["pass"] * 3
    assert len(all_sg["design_0_seq_0"]) == 6
    assert sorted(o for _s, o, _i in all_sg["design_0_seq_0"]) == (
        ["off_target"] * 3 + ["pass"] * 3)
    # The composite ranks the heatmap's rows, best at top.
    assert comps["design_0_seq_0"] == pytest.approx(0.70 - 0.05 * 2.0)


@pytest.mark.logic
def test_load_seed_outcomes_skips_controls_and_the_cold_start_baseline(tmp_path):
    """The cycle-0 initial row has no sequence_group, so it belongs to no
    group; controls are not what the per-seed figures are about."""
    paths = fc.cohort_tree(tmp_path)
    rep, _ = cohort.load_seed_outcomes(paths["runs_dir"], paths["cross_csv"])
    assert not [n for n in rep if n.startswith("input_control")]
    # Four seed rows on disk for design_0_seq_1, one of them the baseline.
    assert len(rep["design_0_seq_1"]) == 3


@pytest.mark.logic
def test_load_seed_outcomes_classifies_a_collapsed_reversion_as_off_target(tmp_path):
    paths = fc.cohort_tree(tmp_path)
    rep, _ = cohort.load_seed_outcomes(paths["runs_dir"], paths["cross_csv"])
    assert {(s, o) for s, o, _i in rep["design_0_seq_1"]} == {
        ("reversion", "off_target")}


@pytest.mark.logic
def test_rep_sg_uses_the_cross_summary_pointer(tmp_path):
    runs = tmp_path / "runs"
    seq = fc.write_sequence_dir(runs, "design_0_seq_0", agg_rows=[
        fc.agg_row_for("no_reversion", steered_ra=1.0, steered_tj=0.9,
                        sequence_group="0"),
        fc.agg_row_for("no_reversion", steered_ra=9.0, steered_tj=0.1,
                        sequence_group="3"),
    ])
    lookup = {"design_0_seq_0": fc.cross_row("design_0_seq_0",
                                              rep_sequence_group="3")}
    assert cohort._rep_sg_for_outcomes(seq, lookup) == "3"
    assert within._rep_sg(seq, lookup) == "3"


@pytest.mark.logic
def test_rep_sg_falls_back_to_the_best_composite_for_a_tier_none_sequence(tmp_path):
    """A tier-none sequence has no pointer. The reverted side decides the
    composite where reversion ran, so the group whose steering looked best can
    lose to one whose reversion held."""
    runs = tmp_path / "runs"
    seq = fc.write_sequence_dir(runs, "design_0_seq_0", agg_rows=[
        fc.agg_row_for("pose_collapses", steered_ra=1.0, steered_tj=0.95,
                        reverted_ra=40.0, reverted_tj=0.01,
                        sequence_group="0"),
        fc.agg_row_for("no_reversion", steered_ra=6.0, steered_tj=0.40,
                        sequence_group="1"),
    ])
    assert cohort._rep_sg_for_outcomes(seq, {}) == "1"
    assert within._rep_sg(seq, {}) == "1"


@pytest.mark.logic
def test_rep_sg_is_none_without_an_aggregated_results_file(tmp_path):
    seq = tmp_path / "runs" / "design_0_seq_0"
    seq.mkdir(parents=True)
    assert cohort._rep_sg_for_outcomes(seq, {}) is None
    assert within._rep_sg(seq, {}) is None


# ═══════════════════════════════════════════════════════════════════════════
# Plot 1 — tier landscape
# ═══════════════════════════════════════════════════════════════════════════

@pytest.mark.logic
def test_plot_tier_landscape_draws(tmp_path):
    out = tmp_path / "negsteer_tier_landscape.png"
    assert cohort.plot_tier_landscape(steered_rows(4), str(out)) is True
    assert_drawn(out)


@pytest.mark.logic
def test_plot_tier_landscape_marks_a_sequence_with_no_composite(tmp_path):
    """A blank Jaccard leaves no composite to plot. The bar would be
    indistinguishable from a genuine zero, so it gets an × instead."""
    rows = steered_rows(2)
    rows[1]["rep_true_jaccard_median"] = ""
    rows[1]["cross_composite_score"] = ""
    rows[1]["cross_rank_by_composite"] = ""
    out = tmp_path / "landscape_with_gap.png"
    assert cohort.plot_tier_landscape(rows, str(out)) is True
    assert_drawn(out)


@pytest.mark.logic
@pytest.mark.parametrize("rows", [
    [],
    [{"mpnn_sequence": "input_control_polyA", "row_type": "control_polyA"}],
], ids=["no-rows", "controls-only"])
def test_plot_tier_landscape_writes_a_placeholder_with_no_steered_rows(tmp_path, rows):
    """Controls sit outside the cross-rank, so a cohort of nothing but controls
    has no bar to draw. Nextflow still publishes the path."""
    out = tmp_path / "negsteer_tier_landscape.png"
    assert cohort.plot_tier_landscape(rows, str(out)) is False
    assert_drawn(out)


# ═══════════════════════════════════════════════════════════════════════════
# Plot 2 — seed outcomes heatmap and bars
# ═══════════════════════════════════════════════════════════════════════════

def _seed_outcomes():
    return (
        {
            "design_0_seq_0": [("steering", "pass", 0), ("steering", "pass", 1),
                                ("reversion", "new_contamination", 2)],
            "design_1_seq_0": [("cold_start", "off_target", 0),
                                ("no_data", "no_data", 1),
                                ("reversion", "multiple_failures", 2)],
        },
        {"design_0_seq_0": 0.60, "design_1_seq_0": -1.2},
    )


@pytest.mark.logic
def test_plot_seed_outcomes_heatmap_draws(tmp_path):
    outcomes, comps = _seed_outcomes()
    out = tmp_path / "negsteer_seed_outcomes_heatmap.png"
    assert cohort.plot_seed_outcomes_heatmap(outcomes, comps, str(out)) is True
    assert_drawn(out)


@pytest.mark.logic
def test_plot_seed_outcomes_heatmap_writes_a_placeholder_with_no_seeds(tmp_path):
    out = tmp_path / "negsteer_seed_outcomes_heatmap.png"
    assert cohort.plot_seed_outcomes_heatmap({}, {}, str(out)) is False
    assert_drawn(out)


@pytest.mark.logic
def test_plot_seed_outcomes_bars_draws(tmp_path):
    outcomes, comps = _seed_outcomes()
    out = tmp_path / "negsteer_seed_outcomes_bars.png"
    assert cohort.plot_seed_outcomes_bars(outcomes, comps, str(out)) is True
    assert_drawn(out)


@pytest.mark.logic
def test_plot_seed_outcomes_bars_writes_a_placeholder_with_no_seeds(tmp_path):
    out = tmp_path / "negsteer_seed_outcomes_bars.png"
    assert cohort.plot_seed_outcomes_bars({}, {}, str(out)) is False
    assert_drawn(out)


# ═══════════════════════════════════════════════════════════════════════════
# Plot 3 — ra_eff vs jaccard
# ═══════════════════════════════════════════════════════════════════════════

@pytest.mark.logic
def test_plot_ra_eff_vs_jaccard_draws_with_controls_and_an_off_scale_failure(tmp_path):
    rows = steered_rows(3)
    rows.append({
        "mpnn_sequence": "design_9_seq_0", "row_type": "steered",
        "cross_tier": "none",
        "rep_ra_eff_vs_truth_median": "40.0",   # far past the scatter clip
        "rep_true_jaccard_median": "0.02",
    })
    rows.append({
        "mpnn_sequence": "input_control_polyA", "row_type": "control_polyA",
        "rep_ra_eff_vs_truth_median": "25.0", "rep_true_jaccard_median": "0.01",
    })
    rows.append({
        "mpnn_sequence": "input_control_scrambled", "row_type": "control_scrambled",
        "rep_ra_eff_vs_truth_median": "22.0", "rep_true_jaccard_median": "0.02",
    })
    out = tmp_path / "negsteer_ra_eff_vs_jaccard.png"
    assert cohort.plot_ra_eff_vs_jaccard(rows, str(out)) is True
    assert_drawn(out)


@pytest.mark.logic
def test_plot_ra_eff_vs_jaccard_draws_with_only_tier_none(tmp_path):
    """No tier A/B/C point means nothing to anchor the scatter's x range to."""
    rows = [{
        "mpnn_sequence": "design_9_seq_0", "row_type": "steered",
        "cross_tier": "none",
        "rep_ra_eff_vs_truth_median": "18.0", "rep_true_jaccard_median": "1.4",
    }]
    out = tmp_path / "ra_none_only.png"
    assert cohort.plot_ra_eff_vs_jaccard(rows, str(out)) is True
    assert_drawn(out)


@pytest.mark.logic
def test_plot_ra_eff_vs_jaccard_draws_with_only_controls(tmp_path):
    rows = [{
        "mpnn_sequence": "input_control_polyA", "row_type": "control_polyA",
        "rep_ra_eff_vs_truth_median": "25.0", "rep_true_jaccard_median": "0.01",
    }]
    out = tmp_path / "ra_controls_only.png"
    assert cohort.plot_ra_eff_vs_jaccard(rows, str(out)) is True
    assert_drawn(out)


@pytest.mark.logic
@pytest.mark.parametrize("rows", [
    [],
    [{"mpnn_sequence": "design_0_seq_0", "rep_ra_eff_vs_truth_median": "nan",
      "rep_true_jaccard_median": "nan"}],
], ids=["no-rows", "nan-metrics"])
def test_plot_ra_eff_vs_jaccard_writes_a_placeholder_with_no_numbers(tmp_path, rows):
    """A nan metric is no data. It must not become a point at the origin."""
    out = tmp_path / "negsteer_ra_eff_vs_jaccard.png"
    assert cohort.plot_ra_eff_vs_jaccard(rows, str(out)) is False
    assert_drawn(out)


# ═══════════════════════════════════════════════════════════════════════════
# Plot 4 — confidence-filter cascade
# ═══════════════════════════════════════════════════════════════════════════

def _cascade_rows():
    """Six steered sequences, each failing exactly one step of the cascade, in
    cascade order. So the waterfall must read 6 → 5 → 4 → 3 → 2 → 1."""
    base = {
        "row_type": "steered",
        "rep_complex_plddt_median": "0.92", "rep_ipae_median": "4.0",
        "rep_pae_pass_frac_median": "0.80", "rep_iptm_median": "0.72",
        "rep_ra_eff_vs_truth_median": "2.0", "rep_true_jaccard_median": "0.70",
    }
    failures = [
        {},                                        # survives everything
        {"rep_complex_plddt_median": "0.50"},      # fails step 1
        {"rep_ipae_median": "30.0"},               # fails step 2
        {"rep_pae_pass_frac_median": "0.01"},      # fails step 3
        {"rep_iptm_median": "0.10"},               # fails step 4
        {"rep_ra_eff_vs_truth_median": "30.0"},    # fails step 5
    ]
    return [{**base, "mpnn_sequence": f"design_{i}_seq_0", **f}
            for i, f in enumerate(failures)]


@pytest.mark.logic
def test_plot_filter_cascade_counts_one_drop_per_step(tmp_path, monkeypatch):
    out = tmp_path / "negsteer_filter_cascade.png"
    heights = _capture_bar_heights(monkeypatch)
    assert cohort.plot_filter_cascade(_cascade_rows(), THRESHOLDS, str(out)) is True
    assert_drawn(out)
    # First six bars are the survivors, the next six the drops stacked on top.
    assert heights[:6] == [6.0, 5.0, 4.0, 3.0, 2.0, 1.0]
    assert heights[6:12] == [0.0, 1.0, 1.0, 1.0, 1.0, 1.0]


@pytest.mark.logic
def test_plot_filter_cascade_treats_a_missing_metric_as_a_drop(tmp_path, monkeypatch):
    """A sequence whose prediction never landed cannot be said to have passed a
    filter it was never measured against."""
    rows = [{"mpnn_sequence": "design_0_seq_0", "row_type": "steered"}]
    out = tmp_path / "cascade_blank.png"
    heights = _capture_bar_heights(monkeypatch)
    assert cohort.plot_filter_cascade(rows, THRESHOLDS, str(out)) is True
    assert heights[:6] == [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]


@pytest.mark.logic
def test_plot_filter_cascade_writes_a_placeholder_with_no_steered_rows(tmp_path):
    out = tmp_path / "negsteer_filter_cascade.png"
    rows = [{"mpnn_sequence": "input_control_polyA", "row_type": "control_polyA"}]
    assert cohort.plot_filter_cascade(rows, THRESHOLDS, str(out)) is False
    assert_drawn(out)


# ═══════════════════════════════════════════════════════════════════════════
# Plot 5 — controls diagnostic
# ═══════════════════════════════════════════════════════════════════════════

@pytest.mark.logic
def test_plot_controls_diagnostic_draws(tmp_path):
    rows = steered_rows(3) + [
        {"mpnn_sequence": "input_control_polyA", "row_type": "control_polyA",
         "rep_ra_eff_vs_truth_median": "25.0", "rep_true_jaccard_median": "0.01"},
        {"mpnn_sequence": "input_control_scrambled", "row_type": "control_scrambled",
         "rep_ra_eff_vs_truth_median": "22.0", "rep_true_jaccard_median": "0.02"},
    ]
    out = tmp_path / "negsteer_controls_diagnostic.png"
    assert cohort.plot_controls_diagnostic(rows, THRESHOLDS, str(out)) is True
    assert_drawn(out)


@pytest.mark.logic
def test_plot_controls_diagnostic_draws_with_an_empty_group(tmp_path):
    """Only one of the four groups has data. The other three keep their x
    positions so the panel still reads as four groups."""
    rows = [{"mpnn_sequence": "design_0_seq_0", "row_type": "steered",
             "cross_tier": "none", "rep_ra_eff_vs_truth_median": "9.0",
             "rep_true_jaccard_median": "0.10"}]
    out = tmp_path / "controls_one_group.png"
    assert cohort.plot_controls_diagnostic(rows, THRESHOLDS, str(out)) is True
    assert_drawn(out)


@pytest.mark.logic
def test_plot_controls_diagnostic_writes_a_placeholder_with_no_numbers(tmp_path):
    rows = [{"mpnn_sequence": "design_0_seq_0", "row_type": "steered",
             "rep_ra_eff_vs_truth_median": "", "rep_true_jaccard_median": "nan"}]
    out = tmp_path / "negsteer_controls_diagnostic.png"
    assert cohort.plot_controls_diagnostic(rows, THRESHOLDS, str(out)) is False
    assert_drawn(out)


# ═══════════════════════════════════════════════════════════════════════════
# Plot 6 — mutation impact
# ═══════════════════════════════════════════════════════════════════════════

def _capture_bar_heights(monkeypatch) -> list:
    """Record the bar heights of the figure being saved.

    The plot functions close their figure straight after savefig, so the only
    place to read what was drawn is between the two.
    """
    heights: list = []
    real_savefig = cohort.plt.savefig

    def recording_savefig(*args, **kwargs):
        fig = cohort.plt.gcf()
        for ax in fig.axes:
            heights.extend(p.get_height() for p in ax.patches)
        return real_savefig(*args, **kwargs)

    monkeypatch.setattr(cohort.plt, "savefig", recording_savefig)
    return heights


def _capture_bar_positions(monkeypatch) -> list:
    """Record which x positions carried a bar on the saved figure."""
    positions: list = []
    real_savefig = cohort.plt.savefig

    def recording_savefig(*args, **kwargs):
        fig = cohort.plt.gcf()
        for ax in fig.axes:
            for p in ax.patches:
                if abs(p.get_width() - 0.8) < 1e-9:
                    positions.append(round(p.get_x() + p.get_width() / 2))
        return real_savefig(*args, **kwargs)

    monkeypatch.setattr(cohort.plt, "savefig", recording_savefig)
    return positions


@pytest.mark.logic
def test_plot_mutation_impact_uses_the_reverted_mutation_set(tmp_path, monkeypatch):
    """Reversion put position 5 back to wild type, so the final structure
    carries mutations at 7 and 26 only. Counting the steered set would credit
    position 5 with a composite score no structure ever had."""
    rows = [{
        "mpnn_sequence": "design_0_seq_0", "row_type": "steered",
        "rep_ra_eff_vs_truth_median": "2.0", "rep_true_jaccard_median": "0.70",
        "rep_outcome": "pose_holds",
        "rep_steered_mutations_chimerax": "/A:5,7,26",
        "rep_reverted_mutations_majority_chimerax": "/A:7,26",
    }]
    out = tmp_path / "negsteer_mutation_impact.png"
    positions = _capture_bar_positions(monkeypatch)
    assert cohort.plot_mutation_impact(rows, None, None, str(out)) is True
    assert sorted(positions) == [7, 26]
    assert_drawn(out)


@pytest.mark.logic
def test_plot_mutation_impact_uses_the_steered_set_when_reversion_never_ran(
        tmp_path, monkeypatch):
    rows = [{
        "mpnn_sequence": "design_0_seq_0", "row_type": "steered",
        "rep_ra_eff_vs_truth_median": "2.0", "rep_true_jaccard_median": "0.70",
        "rep_outcome": "no_reversion",
        "rep_steered_mutations_chimerax": "/A:5,7,26",
        "rep_reverted_mutations_majority_chimerax": "",
    }]
    out = tmp_path / "mutation_impact_steered.png"
    positions = _capture_bar_positions(monkeypatch)
    assert cohort.plot_mutation_impact(rows, None, None, str(out)) is True
    assert sorted(positions) == [5, 7, 26]


@pytest.mark.logic
def test_plot_mutation_impact_falls_back_when_the_reverted_set_is_blank(
        tmp_path, monkeypatch):
    """Reversion ran but recorded no majority mutation set. The steered set is
    then the only record of which positions were touched."""
    rows = [{
        "mpnn_sequence": "design_0_seq_0", "row_type": "steered",
        "rep_ra_eff_vs_truth_median": "2.0", "rep_true_jaccard_median": "0.70",
        "rep_outcome": "pose_collapses",
        "rep_steered_mutations_chimerax": "/A:31",
        "rep_reverted_mutations_majority_chimerax": "",
    }]
    out = tmp_path / "mutation_impact_fallback.png"
    positions = _capture_bar_positions(monkeypatch)
    assert cohort.plot_mutation_impact(rows, None, None, str(out)) is True
    assert positions == [31]


@pytest.mark.logic
def test_plot_mutation_impact_draws_the_design_region_and_the_full_axis(tmp_path):
    """The design-region footprint is context, not data: steering mutations
    cannot fall there by construction. Two runs of positions exercise the
    gap-splitting in the shading."""
    rows = steered_rows(4)
    rows[3]["rep_steered_mutations_chimerax"] = "/A:8"
    out = tmp_path / "mutation_impact_shaded.png"
    assert cohort.plot_mutation_impact(
        rows, {33, 34, 35, 40, 41}, 72, str(out)) is True
    assert_drawn(out)


@pytest.mark.logic
def test_plot_mutation_impact_writes_a_placeholder_with_no_mutations(tmp_path):
    rows = steered_rows(2)
    for r in rows:
        r["rep_steered_mutations_chimerax"] = ""
    out = tmp_path / "negsteer_mutation_impact.png"
    assert cohort.plot_mutation_impact(rows, None, None, str(out)) is False
    assert_drawn(out)


# ═══════════════════════════════════════════════════════════════════════════
# CLI plumbing
# ═══════════════════════════════════════════════════════════════════════════

@pytest.mark.logic
def test_resolve_csv_path_requires_an_existing_file(tmp_path):
    with pytest.raises(SystemExit):
        cohort._resolve_csv_path(None)
    with pytest.raises(SystemExit):
        cohort._resolve_csv_path(str(tmp_path / "absent.csv"))
    present = tmp_path / "cross_sequence_summary.csv"
    present.write_text("mpnn_sequence\n")
    assert cohort._resolve_csv_path(str(present)) == present


@pytest.mark.logic
def test_load_input_design_region_parses_comma_separated_positions(tmp_path):
    p = tmp_path / "input_design_region.txt"
    p.write_text("# contigs: A1-32/10-20/A46-72/6-6 C\n"
                 "\n"
                 "33,34,35\n"
                 "36, not_a_number ,37\n")
    positions, pdb_length = cohort._load_input_design_region(str(p))
    assert positions == {33, 34, 35, 36, 37}
    # The file lists only the design region, so its maximum is all the length
    # this plot can infer.
    assert pdb_length == 37


@pytest.mark.logic
@pytest.mark.parametrize("supplied, body", [
    (False, None),
    (True, None),
    (True, "# nothing but a comment\n"),
], ids=["not-supplied", "missing-file", "no-positions"])
def test_load_input_design_region_degrades_to_no_shading(tmp_path, supplied, body):
    """The mutation-impact plot still renders without the design region."""
    if not supplied:
        arg = None
    else:
        p = tmp_path / "input_design_region.txt"
        if body is not None:
            p.write_text(body)
        arg = str(p)
    assert cohort._load_input_design_region(arg) == (None, None)


@pytest.mark.logic
def test_main_writes_every_published_figure(tmp_path, monkeypatch, capsys):
    """The seven paths in the NEGSTEER_PLOTS publish list, all present and all
    with a figure in them."""
    paths = fc.cohort_tree(tmp_path)
    outdir = tmp_path / "plots"
    monkeypatch.setattr(sys, "argv", [
        "negsteer_plots.py",
        "--csv", str(paths["cross_csv"]),
        "--runs-dir", str(paths["runs_dir"]),
        "--input-design-region", str(paths["design_region"]),
        "--outdir", str(outdir),
    ])
    assert cohort.main() == 0
    expected = [
        "negsteer_tier_landscape.png",
        "negsteer_seed_outcomes_heatmap.png",
        "negsteer_seed_outcomes_bars.png",
        "negsteer_ra_eff_vs_jaccard.png",
        "negsteer_filter_cascade.png",
        "negsteer_controls_diagnostic.png",
        "negsteer_mutation_impact.png",
    ]
    for name in expected:
        assert_drawn(outdir / name)
    assert "FAILED" not in capsys.readouterr().out


@pytest.mark.logic
def test_main_finds_runs_beside_the_csv_and_honours_pdb_length(tmp_path, monkeypatch):
    """--runs-dir is optional; the process can also point at a csv with runs/
    beside it."""
    paths = fc.cohort_tree(tmp_path)
    outdir = tmp_path / "plots"
    monkeypatch.setattr(sys, "argv", [
        "negsteer_plots.py",
        "--csv", str(paths["cross_csv"]),
        "--pdb-length", "72",
        "--complex-plddt-min", "0.60",
        "--outdir", str(outdir),
    ])
    assert cohort.main() == 0
    assert_drawn(outdir / "negsteer_tier_landscape.png")


@pytest.mark.logic
def test_main_warns_and_carries_on_when_runs_dir_is_absent(tmp_path, monkeypatch, capsys):
    cross = tmp_path / "cross_sequence_summary.csv"
    fc.write_cross_summary(cross, [
        fc.cross_row("design_0_seq_0", cross_tier="A", cross_rank_by_composite=1,
                      cross_composite_score=0.6, rep_ra_eff_vs_truth_median=2.0,
                      rep_true_jaccard_median=0.70,
                      rep_steered_mutations_chimerax="/A:5")])
    outdir = tmp_path / "plots"
    monkeypatch.setattr(sys, "argv", [
        "negsteer_plots.py", "--csv", str(cross),
        "--runs-dir", str(tmp_path / "absent"),
        "--input-design-region", str(tmp_path / "absent.txt"),
        "--outdir", str(outdir),
    ])
    assert cohort.main() == 0
    out = capsys.readouterr().out
    assert "WARNING: --runs-dir not found" in out
    # The two per-seed figures have no raw_per_seed_results.csv to read, so
    # they are the placeholder, and every path is still on disk.
    for name in ("negsteer_seed_outcomes_heatmap.png",
                  "negsteer_seed_outcomes_bars.png",
                  "negsteer_tier_landscape.png"):
        assert_drawn(outdir / name)


@pytest.mark.logic
def test_main_returns_zero_on_an_empty_cohort(tmp_path, monkeypatch, capsys):
    cross = tmp_path / "cross_sequence_summary.csv"
    fc.write_cross_summary(cross, [])
    monkeypatch.setattr(sys, "argv", [
        "negsteer_plots.py", "--csv", str(cross),
        "--outdir", str(tmp_path / "plots"),
    ])
    assert cohort.main() == 0
    assert "nothing to plot" in capsys.readouterr().out.lower()
