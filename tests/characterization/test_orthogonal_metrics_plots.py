"""Drive bin/orthogonal_metrics_plots.py against constructed cohort CSVs.

Four figures from two joined CSVs. Two of them make decisions a reader cannot
check from the image: the filter cascade computes a cumulative-AND survival
count, and the combined cohort summary sorts every steered sequence by a
three-key ranking. Both are asserted here on inputs where the answer is
countable by hand, and every threshold that colours a cell is tested on its
boundary.

Column names come from bin/merge_orthogonal_metrics.py (which columns it
appends to the upstream CSV) and bin/cross_sequence_summary.py, with the
committed headers under tests/orthogonal_metrics/ as the check. Thresholds come
from bin/pipeline_thresholds.py through the module's own constants, so a
threshold change moves the test with the code instead of breaking it.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
BIN = REPO / "bin"
HERE = Path(__file__).resolve().parent
for _p in (str(BIN), str(HERE)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import orthogonal_metrics_plots as op  # noqa: E402
from fixtures_plots import (  # noqa: E402
    BOLTZ_CONF_ALL_PASS,
    assert_png_not_blank,
    figures_kept,
    orthogonal_cross_row,
    orthogonal_cross_summary_csv,
    orthogonal_survivor_row,
    orthogonal_survivors_csv,
)
from matplotlib.colors import to_rgba  # noqa: E402

pytestmark = pytest.mark.logic


def _ticklabels(ax, axis="x"):
    getter = ax.get_xticklabels if axis == "x" else ax.get_yticklabels
    return [t.get_text() for t in getter()]


def _rgb(artist):
    return tuple(artist.get_facecolor()[:3])


def _cell(ax, row_idx: int, column: str):
    """The Rectangle for one cell of the combined summary.

    Each row contributes its tier stripe first, then one patch per entry of
    COMBINED_SUMMARY_COLUMNS in order, so the cell index is 1 + the column's
    position within that list.
    """
    n_cols = len(op.COMBINED_SUMMARY_COLUMNS)
    col_idx = [c[0] for c in op.COMBINED_SUMMARY_COLUMNS].index(column)
    return ax.patches[row_idx * (n_cols + 1) + 1 + col_idx]


# A survivor row that passes everything, used as the baseline a test perturbs.
def _clean_survivor(name: str, **overrides) -> dict:
    kwargs = dict(
        mpnn_sequence=name, cross_tier="A", boltz_ra_eff=2.5,
        af3_best_ra_eff=3.0, sc=0.70, bsa=900.0, rosetta_ddg=-40.0,
        interface_plddt=0.80, af3_total_predictions=5,
    )
    kwargs.update(overrides)
    return orthogonal_survivor_row(**kwargs)


# ── Scalar helpers ───────────────────────────────────────────────────

@pytest.mark.parametrize(
    ("raw", "expected"),
    [(None, None), ("", None), ("   ", None), ("abc", None),
     ("1.5", 1.5), (" 2 ", 2.0), (3, 3.0), (4.25, 4.25),
     ("nan", None), ("inf", None), ("-inf", None), (float("nan"), None)],
)
def test_try_float_rejects_blanks_and_non_finite(raw, expected):
    """A CSV cell is a string, and a missing metric is the empty string. NaN and
    the infinities also have to come back as None, because a plot that accepts
    them draws a point at the axis edge instead of omitting it."""
    assert op._try_float(raw) == expected


@pytest.mark.parametrize(
    ("raw", "expected"),
    [(None, None), ("", None), ("3", 3), ("3.0", 3), ("3.7", 3), ("x", None)],
)
def test_try_int_truncates_a_float_string(raw, expected):
    """Mutation counts are written as "3.0" by the median aggregation upstream,
    so the int parse has to fall through to float and truncate."""
    assert op._try_int(raw) == expected


def test_parse_flags_strips_the_value_suffix():
    """orthogonal_flags carries the failing value inline, e.g. sc_too_low:0.512.
    The lookups in FILTER_ORDER match on the bare name, so the suffix goes."""
    assert op._parse_flags("sc_too_low:0.512,bsa_missing") == [
        "sc_too_low", "bsa_missing"]
    assert op._parse_flags("") == []
    assert op._parse_flags(" , sc_missing , ") == ["sc_missing"]


@pytest.mark.parametrize(("raw", "expected"),
                         [("1", True), ("0", False), ("", False), (None, False)])
def test_passes_reads_the_literal_one(raw, expected):
    assert op._passes({"passes_orthogonal_filters": raw}) is expected


def test_short_name_compacts_a_sequence_id():
    assert op._short_name("design_0_seq_1") == "d0_s1"
    assert op._short_name("design_28_seq_0") == "d28_s0"
    assert op._short_name("native_control") == "native_control"


@pytest.mark.parametrize(
    ("value", "fmt", "expected"),
    [(3.0, "int", "3"), (3.6, "int", "4"), (0.5, "0p2", "0.50"),
     (0.5, "0p3", "0.500"), (912.4, "0f", "912"), (-40.26, "1f", "-40.3"),
     (1234.0, "auto", "1234"), (12.34, "auto", "12.3"),
     (0.1234, "auto", "0.123")],
)
def test_format_cell_per_format(value, fmt, expected):
    assert op._format_cell(value, fmt) == expected


def test_boltz_conf_thresholded_columns_comes_from_the_group_span():
    """The tie-breaker's column set is derived from the "Boltz-2 confidence"
    group, so adding a column to that block automatically joins the ranking.
    The five with thresholds are the ones listed; ipsae_min_15 and
    weighted_jaccard have none and are excluded."""
    assert op._boltz_conf_thresholded_columns() == [
        ("rep_complex_plddt_median", 0.70, ">="),
        ("rep_iptm_median", 0.30, ">="),
        ("rep_ipae_median", 15.0, "<="),
        ("rep_pae_pass_frac_median", 0.10, ">="),
        ("rep_interface_plddt_median", 0.75, ">="),
    ]
    assert set(BOLTZ_CONF_ALL_PASS) == {
        c for c, _t, _o in op._boltz_conf_thresholded_columns()}


def test_entered_orthogonal_stage_needs_one_non_empty_presence_column():
    """merge_orthogonal_metrics.py carries forward every upstream row, including
    rows the orthogonal stage never ran on. Those have all five presence columns
    blank and must not pad the cascade's starting count."""
    assert op._entered_orthogonal_stage({"sc": "0.6"}) is True
    assert op._entered_orthogonal_stage({"bsa": "0"}) is True
    assert op._entered_orthogonal_stage(
        {c: "" for c in op.ORTHOG_PRESENCE_COLS}) is False
    assert op._entered_orthogonal_stage({}) is False


# ── _rep_n_mutations ─────────────────────────────────────────────────

def test_rep_n_mutations_cold_start_is_zero():
    """Cold start means no reversion and zero steered mutations, so the final
    prediction carries none. Zero is a real answer here, not a missing one."""
    assert op._rep_n_mutations(orthogonal_cross_row(
        mpnn_sequence="design_0_seq_0", total_mutations=0)) == 0.0


def test_rep_n_mutations_steering_reads_the_steered_count():
    assert op._rep_n_mutations(orthogonal_cross_row(
        mpnn_sequence="design_0_seq_0", total_mutations="3.0")) == 3.0


def test_rep_n_mutations_reversion_reads_the_post_reversion_count():
    """Documented branch, unreachable from the real header:
    cross_sequence_summary.py emits no rep_reversion_verdict or
    rep_reverted_total_mutations column, because extract_passing.py already
    switches total_mutations_median to the post-reversion count when reversion
    ran. The branch is tested on a hand-built row so the documented behaviour is
    pinned, and the committed cross-summary header is the evidence that the
    fallback is what runs in practice.
    """
    row = orthogonal_cross_row(mpnn_sequence="design_0_seq_0",
                               total_mutations="9.0",
                               reversion_verdict="pose_holds",
                               reverted_total_mutations="2")
    assert op._rep_n_mutations(row) == 2.0


def test_rep_n_mutations_reversion_with_no_count_is_missing():
    """A verdict with no count means reversion ran but the number was lost. That
    is missing, not zero, so the cell must hatch rather than claim 0 mutations."""
    row = orthogonal_cross_row(mpnn_sequence="design_0_seq_0",
                               total_mutations="9.0",
                               reversion_verdict="pose_holds")
    assert op._rep_n_mutations(row) is None


def test_rep_n_mutations_is_none_when_nothing_is_populated():
    assert op._rep_n_mutations(
        orthogonal_cross_row(mpnn_sequence="design_0_seq_0")) is None


# ── Plot 1: AF3 vs Boltz ─────────────────────────────────────────────

def test_af3_vs_boltz_groups_points_by_tier(tmp_path):
    """One scatter per non-empty tier, and the legend reports each tier's n. A
    tier string the colour map does not know falls back to "no tier" rather than
    raising, because cross_tier is free text in the CSV."""
    out = tmp_path / "af3.png"
    rows = [
        _clean_survivor("design_0_seq_0", cross_tier="A", af3_best_ra_eff=2.0,
                        boltz_ra_eff=2.0),
        _clean_survivor("design_1_seq_0", cross_tier="A", af3_best_ra_eff=3.0,
                        boltz_ra_eff=3.0),
        _clean_survivor("design_2_seq_0", cross_tier="B", af3_best_ra_eff=4.0,
                        boltz_ra_eff=4.0),
        _clean_survivor("design_3_seq_0", cross_tier="Z", af3_best_ra_eff=9.0,
                        boltz_ra_eff=9.0),
    ]
    with figures_kept() as figs:
        assert op.plot_af3_vs_boltz(rows, str(out)) is True
        texts = [t.get_text() for t in figs()[-1].axes[0].get_legend().get_texts()]
    assert "Tier A (n=2)" in texts
    assert "Tier B (n=1)" in texts
    assert "No tier (n=1)" in texts
    assert_png_not_blank(out)


def test_af3_vs_boltz_axis_reaches_the_data(tmp_path):
    """AF3 routinely returns 30-40 Å on a bad candidate, so the axis cannot stop
    at the threshold. The floor is ORTHOG_AF3_RA_MAX * 1.5 and the data wins when
    it is larger: 40 Å of data gives 44."""
    out = tmp_path / "af3.png"
    rows = [_clean_survivor("design_0_seq_0", af3_best_ra_eff=40.0,
                            boltz_ra_eff=3.0)]
    with figures_kept() as figs:
        op.plot_af3_vs_boltz(rows, str(out))
        ax = figs()[-1].axes[0]
    assert ax.get_xlim() == pytest.approx((0.0, 44.0))
    assert ax.get_ylim() == ax.get_xlim()


def test_af3_vs_boltz_axis_floor_applies_when_the_data_is_small(tmp_path):
    out = tmp_path / "af3.png"
    rows = [_clean_survivor("design_0_seq_0", af3_best_ra_eff=1.0,
                            boltz_ra_eff=1.0)]
    with figures_kept() as figs:
        op.plot_af3_vs_boltz(rows, str(out))
        lim = figs()[-1].axes[0].get_xlim()
    assert lim == pytest.approx((0.0, op.ORTHOG_AF3_RA_MAX * 1.5))


def test_af3_vs_boltz_draws_rows_without_af3_at_the_right_edge(tmp_path):
    """AF3 can be missing while Boltz is present. Those rows become right-arrow
    markers at the axis edge with a count, instead of vanishing."""
    out = tmp_path / "af3.png"
    rows = [
        _clean_survivor("design_0_seq_0", af3_best_ra_eff=2.0, boltz_ra_eff=2.0),
        _clean_survivor("design_1_seq_0", af3_best_ra_eff="", boltz_ra_eff=4.0),
    ]
    with figures_kept() as figs:
        assert op.plot_af3_vs_boltz(rows, str(out)) is True
        texts = [t.get_text() for t in figs()[-1].axes[0].texts]
    assert any("AF3 missing\n(n=1)" == t for t in texts)
    assert_png_not_blank(out)


def test_af3_vs_boltz_with_only_missing_af3_rows_still_plots(tmp_path):
    """No AF3 anywhere, which is what a run with the AF3 stream switched off
    produces. With no x data to reach for, the x range falls back to twice the
    threshold, and the axes are square so both ends come out at 10."""
    out = tmp_path / "af3.png"
    rows = [_clean_survivor("design_0_seq_0", af3_best_ra_eff="",
                            boltz_ra_eff=4.0)]
    with figures_kept() as figs:
        assert op.plot_af3_vs_boltz(rows, str(out)) is True
        lim = figs()[-1].axes[0].get_xlim()
    assert lim == pytest.approx((0.0, op.ORTHOG_AF3_RA_MAX * 2))
    assert_png_not_blank(out)


def test_af3_vs_boltz_skips_rows_with_no_boltz_value(tmp_path):
    """Boltz ra_eff is the y-axis, so a row without it cannot be placed. A row
    with neither value is skipped before the tier lookup."""
    out = tmp_path / "af3.png"
    rows = [
        _clean_survivor("design_0_seq_0", af3_best_ra_eff=2.0, boltz_ra_eff=""),
        _clean_survivor("design_1_seq_0", af3_best_ra_eff="", boltz_ra_eff=""),
        _clean_survivor("design_2_seq_0", af3_best_ra_eff=2.0, boltz_ra_eff=2.0),
    ]
    with figures_kept() as figs:
        op.plot_af3_vs_boltz(rows, str(out))
        texts = [t.get_text() for t in figs()[-1].axes[0].get_legend().get_texts()]
    assert "Tier A (n=1)" in texts


def test_af3_vs_boltz_with_no_ra_eff_data_returns_false_and_writes_a_file(tmp_path):
    """No-data path. The return value is the documented False, and the
    placeholder still has to appear because ORTHOG_PLOTS publishes a fixed
    figure list."""
    out = tmp_path / "af3.png"
    rows = [_clean_survivor("design_0_seq_0", af3_best_ra_eff="",
                            boltz_ra_eff="")]
    assert op.plot_af3_vs_boltz(rows, str(out)) is False
    assert_png_not_blank(out)


def test_af3_vs_boltz_with_no_rows_returns_false_and_writes_a_file(tmp_path):
    out = tmp_path / "af3.png"
    assert op.plot_af3_vs_boltz([], str(out)) is False
    assert_png_not_blank(out)


# ── Plot 2: filter cascade ───────────────────────────────────────────

def test_filter_cascade_counts_survivors_at_each_step(tmp_path):
    """Four rows that entered the stage. The sc step drops the sc_too_low row,
    the bsa step drops the bsa_missing row, the interface_plddt step drops nobody
    because nothing in the pipeline emits that flag, and the trailing
    informational AF3 bar drops the row AF3 disagreed with. So the green bars
    read 4, 3, 2, 2, 1 and the stacked drops read 0, 1, 1, 0, 1.
    """
    out = tmp_path / "cascade.png"
    rows = [
        _clean_survivor("design_0_seq_0"),
        _clean_survivor("design_1_seq_0", sc=0.40,
                        flags=["sc_too_low:0.400"]),
        _clean_survivor("design_2_seq_0", bsa="",
                        flags=["bsa_missing"]),
        _clean_survivor("design_3_seq_0", af3_best_ra_eff=30.0,
                        flags=["af3_nomsa_ra_eff_too_high:30.00"]),
    ]
    with figures_kept() as figs:
        assert op.plot_filter_cascade(rows, str(out)) is True
        ax = figs()[-1].axes[0]
        passing = [p.get_height() for p in ax.containers[0].patches]
        dropped = [p.get_height() for p in ax.containers[1].patches]
    assert passing == [4, 3, 2, 2, 1]
    assert dropped == [0, 1, 1, 0, 1]
    assert _ticklabels(ax)[0] == "All survivors"
    assert "informational" in _ticklabels(ax)[-1]
    assert_png_not_blank(out)


def test_filter_cascade_charges_a_doubly_failing_row_to_the_first_filter(tmp_path):
    """The cascade is a cumulative AND, so a row failing both sc and bsa is
    dropped once, at sc. Counting it twice would make the drops sum above the
    cohort size."""
    out = tmp_path / "cascade.png"
    rows = [
        _clean_survivor("design_0_seq_0"),
        _clean_survivor("design_1_seq_0", sc=0.1, bsa=10.0,
                        flags=["sc_too_low:0.100", "bsa_too_low:10"]),
    ]
    with figures_kept() as figs:
        op.plot_filter_cascade(rows, str(out))
        ax = figs()[-1].axes[0]
        passing = [p.get_height() for p in ax.containers[0].patches]
        dropped = [p.get_height() for p in ax.containers[1].patches]
    assert passing == [2, 1, 1, 1, 1]
    assert dropped == [0, 1, 0, 0, 0]


def test_filter_cascade_excludes_rows_that_never_entered_the_stage(tmp_path):
    """Two of the three rows have every presence column blank, so the cascade
    starts from one, not three. Starting from three would make the stage look
    like it filtered rows it never saw."""
    out = tmp_path / "cascade.png"
    rows = [
        _clean_survivor("design_0_seq_0"),
        orthogonal_survivor_row(mpnn_sequence="design_1_seq_0",
                                flags=["sc_missing", "bsa_missing",
                                       "af3_nomsa_missing"]),
        orthogonal_survivor_row(mpnn_sequence="design_2_seq_0",
                                flags=["sc_missing", "bsa_missing",
                                       "af3_nomsa_missing"]),
    ]
    with figures_kept() as figs:
        op.plot_filter_cascade(rows, str(out))
        passing = [p.get_height()
                   for p in figs()[-1].axes[0].containers[0].patches]
    assert passing == [1, 1, 1, 1, 1]


def test_filter_cascade_drops_a_row_whose_metric_is_missing(tmp_path):
    """A missing metric is a failure, not a pass: sc_missing drops the row at the
    sc step exactly as sc_too_low would."""
    out = tmp_path / "cascade.png"
    rows = [_clean_survivor("design_0_seq_0", sc="", flags=["sc_missing"])]
    with figures_kept() as figs:
        op.plot_filter_cascade(rows, str(out))
        passing = [p.get_height()
                   for p in figs()[-1].axes[0].containers[0].patches]
    assert passing == [1, 0, 0, 0, 0]


def test_filter_cascade_with_no_rows_returns_false_and_writes_a_file(tmp_path):
    out = tmp_path / "cascade.png"
    assert op.plot_filter_cascade([], str(out)) is False
    assert_png_not_blank(out)


def test_filter_cascade_with_nothing_entering_the_stage_returns_false(tmp_path):
    """No-data path that is not an empty CSV: every row was carried forward from
    upstream with blank orthogonal columns."""
    out = tmp_path / "cascade.png"
    rows = [orthogonal_survivor_row(mpnn_sequence="design_0_seq_0")]
    assert op.plot_filter_cascade(rows, str(out)) is False
    assert_png_not_blank(out)


# ── Plot 3: metrics vs composite ─────────────────────────────────────

def test_metrics_vs_composite_draws_one_panel_per_metric(tmp_path):
    """Five metrics in a 3x2 grid, so the sixth cell is switched off rather than
    left as empty axes with ticks."""
    out = tmp_path / "mvc.png"
    rows = [_clean_survivor("design_0_seq_0")]
    cross = {"design_0_seq_0": orthogonal_cross_row(
        mpnn_sequence="design_0_seq_0", composite=0.5)}
    with figures_kept() as figs:
        assert op.plot_metrics_vs_composite(rows, cross, str(out)) is True
        axes = figs()[-1].axes
    assert len(axes) == 6
    assert len(op.METRIC_VS_COMPOSITE_PANELS) == 5
    assert axes[-1].axison is False
    assert_png_not_blank(out)


def test_metrics_vs_composite_bakes_the_threshold_into_the_axis_label(tmp_path):
    """Per-panel legends were replaced by one shared figure legend, so each
    panel's threshold is carried in its x-label with a unicode comparison glyph
    that states the direction."""
    out = tmp_path / "mvc.png"
    rows = [_clean_survivor("design_0_seq_0")]
    cross = {"design_0_seq_0": orthogonal_cross_row(
        mpnn_sequence="design_0_seq_0", composite=0.5)}
    with figures_kept() as figs:
        op.plot_metrics_vs_composite(rows, cross, str(out))
        labels = [ax.get_xlabel() for ax in figs()[-1].axes]
    assert labels[0] == (
        f"Sc (shape complementarity)  [thr ≥ {op.ORTHOG_SC_MIN:g}]")
    assert labels[4] == (
        f"AF3-no-MSA best ra_eff  [thr ≤ {op.ORTHOG_AF3_RA_MAX:g} Å]")


def test_metrics_vs_composite_bounds_the_fraction_metrics_to_zero_one(tmp_path):
    """Sc and interface pLDDT are [0,1] by construction, so their panels are
    pinned to that range however narrow the cohort's spread is. BSA is unbounded
    and scales to its data."""
    out = tmp_path / "mvc.png"
    rows = [_clean_survivor("design_0_seq_0", sc=0.70, bsa=900.0)]
    cross = {"design_0_seq_0": orthogonal_cross_row(
        mpnn_sequence="design_0_seq_0", composite=0.5)}
    with figures_kept() as figs:
        op.plot_metrics_vs_composite(rows, cross, str(out))
        axes = figs()[-1].axes
    assert axes[0].get_xlim() == (0.0, 1.0)
    assert axes[1].get_xlim()[1] > 900.0


def test_metrics_vs_composite_joins_the_composite_from_cross_summary(tmp_path):
    """survivors_with_orthogonal_metrics.csv does not carry the composite score,
    so the y value comes from cross_summary joined on mpnn_sequence. A survivor
    with no cross row has no y and is dropped rather than plotted at zero."""
    out = tmp_path / "mvc.png"
    rows = [_clean_survivor("design_0_seq_0", sc=0.70),
            _clean_survivor("design_9_seq_9", sc=0.90)]
    cross = {"design_0_seq_0": orthogonal_cross_row(
        mpnn_sequence="design_0_seq_0", composite=0.42)}
    with figures_kept() as figs:
        op.plot_metrics_vs_composite(rows, cross, str(out))
        points = figs()[-1].axes[0].collections[0].get_offsets()
    assert len(points) == 1
    assert tuple(points[0]) == pytest.approx((0.70, 0.42))


def test_metrics_vs_composite_annotates_names_only_at_small_n(tmp_path):
    """Names are readable at 12 points and noise above it. The comparison is
    non-strict, so exactly 12 is annotated and 13 is not."""
    out = tmp_path / "mvc.png"
    for n, expect_labels in ((12, True), (13, False)):
        rows = [_clean_survivor(f"design_{i}_seq_0", sc=0.60 + i * 0.01)
                for i in range(n)]
        cross = {r["mpnn_sequence"]: orthogonal_cross_row(
            mpnn_sequence=r["mpnn_sequence"], composite=0.5) for r in rows}
        with figures_kept() as figs:
            op.plot_metrics_vs_composite(rows, cross, str(out))
            annotated = bool(len(figs()[-1].axes[0].texts))
        assert annotated is expect_labels


def test_metrics_vs_composite_marks_an_empty_panel(tmp_path):
    """A metric absent across the whole cohort gets a "No <col> data" note and
    its axes switched off, so an empty panel is distinguishable from a panel
    whose points happen to sit off-screen."""
    out = tmp_path / "mvc.png"
    rows = [_clean_survivor("design_0_seq_0", rosetta_ddg="")]
    cross = {"design_0_seq_0": orthogonal_cross_row(
        mpnn_sequence="design_0_seq_0", composite=0.5)}
    with figures_kept() as figs:
        assert op.plot_metrics_vs_composite(rows, cross, str(out)) is True
        ddg_panel = figs()[-1].axes[3]
        assert [t.get_text() for t in ddg_panel.texts] == ["No rosetta_ddg\ndata"]
        assert ddg_panel.axison is False
    assert_png_not_blank(out)


def test_metrics_vs_composite_with_no_rows_returns_false_and_writes_a_file(tmp_path):
    out = tmp_path / "mvc.png"
    assert op.plot_metrics_vs_composite([], {}, str(out)) is False
    assert_png_not_blank(out)


def test_metrics_vs_composite_with_no_joinable_data_returns_false(tmp_path):
    """No-data path. Rows exist but none joins to a composite score, so every
    panel is empty and the documented False comes back with a placeholder."""
    out = tmp_path / "mvc.png"
    rows = [_clean_survivor("design_0_seq_0")]
    assert op.plot_metrics_vs_composite(rows, {}, str(out)) is False
    assert_png_not_blank(out)


# ── Plot 4: combined cohort summary ──────────────────────────────────

def test_combined_summary_sorts_by_tier_then_confidence_then_composite(tmp_path):
    """Three keys, in order. Tier wins outright, so the tier-B row sits below
    both tier-A rows despite holding the highest composite. Within tier A, the
    row that passes every thresholded Boltz-2 confidence metric is promoted above
    the one that fails one, despite the failing row's higher composite. So the
    rows read A/clean (0.10), A/failing (0.90), B (0.99).
    """
    out = tmp_path / "summary.png"
    failing_conf = dict(BOLTZ_CONF_ALL_PASS, rep_iptm_median="0.10")
    cross = [
        orthogonal_cross_row(mpnn_sequence="design_1_seq_0", cross_tier="A",
                             composite=0.90, boltz_conf=failing_conf),
        orthogonal_cross_row(mpnn_sequence="design_2_seq_0", cross_tier="B",
                             composite=0.99),
        orthogonal_cross_row(mpnn_sequence="design_3_seq_0", cross_tier="A",
                             composite=0.10),
    ]
    with figures_kept() as figs:
        assert op.plot_combined_cohort_orthogonal_summary(
            cross, {}, str(out)) is True
        labels = _ticklabels(figs()[-1].axes[0], "y")
    assert labels == ["A  d3_s0", "A  d1_s0", "B  d2_s0"]
    assert_png_not_blank(out)


def test_combined_summary_sorts_a_missing_composite_to_the_bottom(tmp_path):
    """Within a tier-and-confidence group, a row with no composite has not earned
    a position, so it goes last rather than ranking as zero."""
    out = tmp_path / "summary.png"
    cross = [
        orthogonal_cross_row(mpnn_sequence="design_1_seq_0", cross_tier="A",
                             composite=""),
        orthogonal_cross_row(mpnn_sequence="design_2_seq_0", cross_tier="A",
                             composite=-5.0),
    ]
    with figures_kept() as figs:
        op.plot_combined_cohort_orthogonal_summary(cross, {}, str(out))
        labels = _ticklabels(figs()[-1].axes[0], "y")
    assert labels == ["A  d2_s0", "A  d1_s0"]


def test_combined_summary_ranks_an_unknown_tier_with_the_untiered(tmp_path):
    """cross_tier is free text. An unrecognised value sorts as "none", so it goes
    last despite its higher composite, and its stripe takes the no-tier grey.
    The row label still shows the raw string rather than the "·" marker, because
    only the literal "none" is rewritten.
    """
    out = tmp_path / "summary.png"
    cross = [
        orthogonal_cross_row(mpnn_sequence="design_1_seq_0", cross_tier="Z",
                             composite=0.99),
        orthogonal_cross_row(mpnn_sequence="design_2_seq_0", cross_tier="C",
                             composite=0.10),
    ]
    with figures_kept() as figs:
        op.plot_combined_cohort_orthogonal_summary(cross, {}, str(out))
        ax = figs()[-1].axes[0]
        labels = _ticklabels(ax, "y")
        stripes = [_rgb(ax.patches[0]),
                   _rgb(ax.patches[len(op.COMBINED_SUMMARY_COLUMNS) + 1])]
    assert labels == ["C  d2_s0", "Z  d1_s0"]
    assert stripes[1] == pytest.approx(to_rgba(op.COLOUR_TIER["none"])[:3])


def test_combined_summary_excludes_control_rows(tmp_path):
    """cross_sequence_summary.py keeps controls in the CSV with row_type set to
    something other than 'steered', and they are not ranked. A blank row_type
    defaults to steered, which is what an older CSV without the column has."""
    out = tmp_path / "summary.png"
    cross = [
        orthogonal_cross_row(mpnn_sequence="design_1_seq_0", cross_tier="A",
                             composite=0.5),
        orthogonal_cross_row(mpnn_sequence="native_control", cross_tier="A",
                             composite=0.9, row_type="control"),
        orthogonal_cross_row(mpnn_sequence="design_3_seq_0", cross_tier="A",
                             composite=0.4, row_type=""),
    ]
    with figures_kept() as figs:
        op.plot_combined_cohort_orthogonal_summary(cross, {}, str(out))
        labels = _ticklabels(figs()[-1].axes[0], "y")
    assert labels == ["A  d1_s0", "A  d3_s0"]


@pytest.mark.parametrize(
    ("column", "value", "expect_pass"),
    [("sc", op.ORTHOG_SC_MIN, True),
     ("sc", op.ORTHOG_SC_MIN - 0.001, False),
     ("bsa", op.ORTHOG_BSA_MIN, True),
     ("bsa", op.ORTHOG_BSA_MIN - 1.0, False),
     ("af3_nomsa_best_ra_eff", op.ORTHOG_AF3_RA_MAX, True),
     ("af3_nomsa_best_ra_eff", op.ORTHOG_AF3_RA_MAX + 0.001, False),
     ("rosetta_ddg", op.DDG_BENNETT_REFERENCE, True),
     ("rosetta_ddg", op.DDG_BENNETT_REFERENCE + 0.001, False)],
)
def test_combined_summary_cell_thresholds_are_non_strict(tmp_path, column,
                                                         value, expect_pass):
    """Both comparisons are non-strict: a ">=" cell passes at exactly the
    threshold and a "<=" cell passes at exactly the threshold. So a value sitting
    on the line is green, and only a value past it is red.
    """
    out = tmp_path / "summary.png"
    kwargs = {"sc": 0.70, "bsa": 900.0, "af3_best_ra_eff": 3.0,
              "rosetta_ddg": -40.0}
    rename = {"sc": "sc", "bsa": "bsa",
              "af3_nomsa_best_ra_eff": "af3_best_ra_eff",
              "rosetta_ddg": "rosetta_ddg"}
    kwargs[rename[column]] = value
    survivor = _clean_survivor("design_1_seq_0", **kwargs)
    cross = [orthogonal_cross_row(mpnn_sequence="design_1_seq_0",
                                  cross_tier="A", composite=0.5)]
    with figures_kept() as figs:
        op.plot_combined_cohort_orthogonal_summary(
            cross, {"design_1_seq_0": survivor}, str(out))
        cell = _cell(figs()[-1].axes[0], 0, column)
    want = op.COLOUR_PASS if expect_pass else op.COLOUR_FAIL
    assert _rgb(cell) == pytest.approx(to_rgba(want)[:3])


def test_combined_summary_hatches_cells_for_a_sequence_that_never_reached_orthogonal(
        tmp_path):
    """Most of the cohort never reaches the orthogonal stage. Those rows get
    grey-hatched cells with an em dash in the four orthogonal columns, which is
    distinguishable from a red failing cell. Hatching is per cell, not per row:
    the populated columns of the same row keep their fill.

    This fixture also leaves six cross columns blank, so the em dashes number
    ten, not four: the four orthogonal columns plus n_mut, Boltz ra_eff, rec
    RMSD, true_jaccard, ipsae_min_15 and weighted_jacc.
    """
    out = tmp_path / "summary.png"
    cross = [orthogonal_cross_row(mpnn_sequence="design_1_seq_0",
                                  cross_tier="A", composite=0.5)]
    with figures_kept() as figs:
        op.plot_combined_cohort_orthogonal_summary(cross, {}, str(out))
        ax = figs()[-1].axes[0]
        for column in ("af3_nomsa_best_ra_eff", "sc", "bsa", "rosetta_ddg"):
            assert _cell(ax, 0, column).get_hatch() == "///"
        assert _cell(ax, 0, "cross_composite_score").get_hatch() is None
        assert [t.get_text() for t in ax.texts].count("—") == 10
    assert_png_not_blank(out)


def test_combined_summary_uses_a_neutral_fill_for_unthresholded_metrics(tmp_path):
    """true_jaccard has no threshold, so its cell is light blue and still carries
    the value. Colouring it red or green would imply a gate that does not exist."""
    out = tmp_path / "summary.png"
    cross = [orthogonal_cross_row(mpnn_sequence="design_1_seq_0",
                                  cross_tier="A", composite=0.5,
                                  rep_true_jaccard_median=0.62)]
    with figures_kept() as figs:
        op.plot_combined_cohort_orthogonal_summary(cross, {}, str(out))
        ax = figs()[-1].axes[0]
        cell = _cell(ax, 0, "rep_true_jaccard_median")
        assert _rgb(cell) == pytest.approx(to_rgba(op.COLOUR_NEUTRAL)[:3])
        assert "0.62" in [t.get_text() for t in ax.texts]


def test_combined_summary_draws_the_three_axis_group_headers(tmp_path):
    """The column blocks carry a header band each, except the descriptive
    columns, whose group label is None."""
    out = tmp_path / "summary.png"
    cross = [orthogonal_cross_row(mpnn_sequence="design_1_seq_0",
                                  cross_tier="A", composite=0.5)]
    with figures_kept() as figs:
        op.plot_combined_cohort_orthogonal_summary(cross, {}, str(out))
        texts = [t.get_text() for t in figs()[-1].axes[0].texts]
    assert "Boltz-2 structure" in texts
    assert "Boltz-2 confidence" in texts
    assert "Orthogonal" in texts
    assert "tier" in texts


def test_combined_summary_tier_stripe_takes_the_tier_colour(tmp_path):
    out = tmp_path / "summary.png"
    cross = [orthogonal_cross_row(mpnn_sequence="design_1_seq_0",
                                  cross_tier="B", composite=0.5)]
    with figures_kept() as figs:
        op.plot_combined_cohort_orthogonal_summary(cross, {}, str(out))
        stripe = figs()[-1].axes[0].patches[0]
    assert _rgb(stripe) == pytest.approx(to_rgba(op.COLOUR_TIER["B"])[:3])


def test_combined_summary_with_no_steered_rows_returns_false_and_writes_a_file(
        tmp_path):
    """No-data path: a cohort of controls only. The documented False comes back
    and the placeholder is written."""
    out = tmp_path / "summary.png"
    cross = [orthogonal_cross_row(mpnn_sequence="native_control",
                                  row_type="control")]
    assert op.plot_combined_cohort_orthogonal_summary(cross, {}, str(out)) is False
    assert_png_not_blank(out)


def test_combined_summary_with_no_rows_returns_false_and_writes_a_file(tmp_path):
    out = tmp_path / "summary.png"
    assert op.plot_combined_cohort_orthogonal_summary([], {}, str(out)) is False
    assert_png_not_blank(out)


def test_combined_summary_single_row(tmp_path):
    """One steered sequence is the realistic small-cohort case and the figure
    height is computed from the row count, so a single row is where a degenerate
    height would show up."""
    out = tmp_path / "summary.png"
    cross = [orthogonal_cross_row(mpnn_sequence="design_1_seq_0",
                                  cross_tier="A", composite=0.5)]
    with figures_kept() as figs:
        op.plot_combined_cohort_orthogonal_summary(cross, {}, str(out))
        assert len(_ticklabels(figs()[-1].axes[0], "y")) == 1
    assert_png_not_blank(out)


# ── CLI ──────────────────────────────────────────────────────────────

PUBLISHED = ("orthogonal_af3_vs_boltz.png", "orthogonal_filter_cascade.png",
             "orthogonal_metrics_vs_composite.png",
             "orthogonal_combined_cohort_summary.png")


def _write_pair(tmp_path, survivors, cross):
    return (orthogonal_survivors_csv(tmp_path / "survivors.csv", survivors),
            orthogonal_cross_summary_csv(tmp_path / "cross.csv", cross))


def test_main_writes_all_four_published_figures(tmp_path, monkeypatch):
    """modules/negsteer_orthogonal_metrics.nf publishes these four names, so one
    run has to produce every one of them."""
    survivors = [
        _clean_survivor("design_0_seq_0"),
        _clean_survivor("design_1_seq_0", cross_tier="B", sc=0.40,
                        flags=["sc_too_low:0.400"]),
    ]
    cross = [
        orthogonal_cross_row(mpnn_sequence="design_0_seq_0", cross_tier="A",
                             composite=0.58, total_mutations="0.0",
                             rep_true_jaccard_median=0.60,
                             rep_ra_eff_vs_truth_median=2.72),
        orthogonal_cross_row(mpnn_sequence="design_1_seq_0", cross_tier="B",
                             composite=0.39, total_mutations="3.0",
                             rep_ra_eff_vs_truth_median=4.07),
        orthogonal_cross_row(mpnn_sequence="native_control", row_type="control",
                             composite=0.99),
    ]
    surv_csv, cross_csv = _write_pair(tmp_path, survivors, cross)
    outdir = tmp_path / "plots"
    monkeypatch.setattr(sys, "argv", [
        "orthogonal_metrics_plots.py",
        "--survivors-csv", str(surv_csv),
        "--cross-summary-csv", str(cross_csv),
        "--outdir", str(outdir),
    ])
    assert op.main() == 0
    for name in PUBLISHED:
        assert_png_not_blank(outdir / name)


def test_main_creates_the_output_directory(tmp_path, monkeypatch):
    surv_csv, cross_csv = _write_pair(
        tmp_path, [_clean_survivor("design_0_seq_0")],
        [orthogonal_cross_row(mpnn_sequence="design_0_seq_0", cross_tier="A",
                              composite=0.5)])
    outdir = tmp_path / "a" / "b" / "plots"
    monkeypatch.setattr(sys, "argv", [
        "orthogonal_metrics_plots.py", "--survivors-csv", str(surv_csv),
        "--cross-summary-csv", str(cross_csv), "--outdir", str(outdir)])
    assert op.main() == 0
    assert outdir.is_dir()


def test_main_on_an_empty_survivors_csv_returns_early(tmp_path, monkeypatch):
    """Header but no rows. main() returns 0 without writing any figure, so the
    Nextflow stage's output pattern has to tolerate the absence. This is the one
    path where the published list is not produced, which is worth knowing.
    """
    surv_csv, cross_csv = _write_pair(tmp_path, [], [])
    outdir = tmp_path / "plots"
    monkeypatch.setattr(sys, "argv", [
        "orthogonal_metrics_plots.py", "--survivors-csv", str(surv_csv),
        "--cross-summary-csv", str(cross_csv), "--outdir", str(outdir)])
    assert op.main() == 0
    assert list(outdir.glob("*.png")) == []


def test_main_rejects_a_missing_survivors_csv(tmp_path, monkeypatch):
    _surv, cross_csv = _write_pair(tmp_path, [], [])
    monkeypatch.setattr(sys, "argv", [
        "orthogonal_metrics_plots.py",
        "--survivors-csv", str(tmp_path / "absent.csv"),
        "--cross-summary-csv", str(cross_csv),
        "--outdir", str(tmp_path / "plots")])
    with pytest.raises(SystemExit) as exc:
        op.main()
    assert "--survivors-csv path does not exist" in str(exc.value)


def test_main_rejects_a_missing_cross_summary_csv(tmp_path, monkeypatch):
    surv_csv, _cross = _write_pair(
        tmp_path, [_clean_survivor("design_0_seq_0")], [])
    monkeypatch.setattr(sys, "argv", [
        "orthogonal_metrics_plots.py", "--survivors-csv", str(surv_csv),
        "--cross-summary-csv", str(tmp_path / "absent.csv"),
        "--outdir", str(tmp_path / "plots")])
    with pytest.raises(SystemExit) as exc:
        op.main()
    assert "--cross-summary-csv path does not exist" in str(exc.value)


def test_main_reports_a_failing_plot_without_aborting_the_rest(
        tmp_path, monkeypatch, capsys):
    """Each plot is wrapped, so one raising does not cost the other three. That
    also means a broken plot leaves its published path missing while the stage
    still exits 0, which is the trade-off being pinned here.
    """
    surv_csv, cross_csv = _write_pair(
        tmp_path, [_clean_survivor("design_0_seq_0")],
        [orthogonal_cross_row(mpnn_sequence="design_0_seq_0", cross_tier="A",
                              composite=0.5)])
    outdir = tmp_path / "plots"

    def _boom(*_a, **_k):
        raise RuntimeError("synthetic failure")

    monkeypatch.setattr(op, "plot_af3_vs_boltz", _boom)
    monkeypatch.setattr(sys, "argv", [
        "orthogonal_metrics_plots.py", "--survivors-csv", str(surv_csv),
        "--cross-summary-csv", str(cross_csv), "--outdir", str(outdir)])
    assert op.main() == 0
    assert "FAILED: RuntimeError: synthetic failure" in capsys.readouterr().out
    assert not (outdir / "orthogonal_af3_vs_boltz.png").exists()
    assert_png_not_blank(outdir / "orthogonal_filter_cascade.png")
