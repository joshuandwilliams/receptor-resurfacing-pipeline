"""Drive bin/rosetta_filter_plots.py against constructed rosetta_filter_metrics.json.

One figure, and everything in it is a claim about the Sc filter that gates the
whole RFDiffusion-to-MPNN handover. Three of those claims are only ever visible
inside the PNG, so each is asserted on an input whose right answer is
arithmetic:

  - the filtered count in the legend, which is the count of designs the
    pipeline is about to throw away, and whose boundary case is a design
    sitting exactly on the threshold;
  - the shaded region's left edge, which the module goes out of its way to pin
    to the axis spine and which silently drifts if `set_xlim` ever moves after
    the `axvspan`;
  - the bin count, which switches between Freedman-Diaconis and a fixed count
    at twelve designs and is clamped at both ends.

`plot_sc_histogram` and `make_empty_plot` write to bare relative filenames, so
every test chdirs into tmp_path first. Key names come from
tests/rosetta_filtering/example_output_files/rosetta_filtering/
rosetta_filter_metrics.json.
"""
from __future__ import annotations

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

import matplotlib  # noqa: E402

matplotlib.use("Agg")

import rosetta_filter_plots as rp  # noqa: E402
from fixtures_plots import (  # noqa: E402
    ROSETTA_DEFAULT_SC_THRESHOLD,
    assert_png_not_blank,
    figures_kept,
    rosetta_design,
    rosetta_filter_metrics_json,
)
from matplotlib.colors import to_hex  # noqa: E402

HISTOGRAM = "rosetta_sc_histogram.png"


def _legend_texts(ax):
    legend = ax.get_legend()
    return [] if legend is None else [t.get_text() for t in legend.get_texts()]


def _patches_coloured(ax, colour):
    """Patches whose face is `colour`, which separates the histogram bars from
    the shaded filtered region without depending on which patch class
    `axvspan` happens to return in this matplotlib."""
    want = to_hex(colour).lower()
    return [p for p in ax.patches if to_hex(p.get_facecolor()).lower() == want]


def _x_extent(patch):
    """(left, right) of a patch in data coordinates, for either a Polygon
    (identity transform, vertices already in data space) or a Rectangle
    (unit square mapped by its patch transform)."""
    verts = patch.get_patch_transform().transform(patch.get_path().vertices)
    return float(verts[:, 0].min()), float(verts[:, 0].max())


def _bars(ax):
    return _patches_coloured(ax, rp.COLOUR_ALL)


# ── The filtered count ───────────────────────────────────────────────

@pytest.mark.logic
def test_the_threshold_is_inclusive_so_a_design_exactly_on_it_is_kept(
        tmp_path, monkeypatch):
    """The count is `np.sum(sc_vals >= threshold)`, so the comparison is
    non-strict and a design at exactly the threshold passes. Three designs at
    0.49, 0.50 and 0.51 against a threshold of 0.50 therefore give one
    filtered, not two. Getting this wrong discards a design the pipeline
    kept, and the legend is the only place the number is ever written down.
    """
    monkeypatch.chdir(tmp_path)
    designs = [rosetta_design(i, sc_value=v)
               for i, v in enumerate((0.49, 0.50, 0.51))]
    with figures_kept() as figs:
        rp.plot_sc_histogram(designs, 0.50)
        texts = _legend_texts(figs()[-1].axes[0])
    assert texts == ["Sc < 0.50 (filtered, n=1 of 3)", "Threshold = 0.50"]
    assert_png_not_blank(tmp_path / HISTOGRAM)


@pytest.mark.logic
def test_designs_with_no_sc_value_are_left_out_of_both_counts(
        tmp_path, monkeypatch):
    """A design Rosetta failed on carries sc_value null. It is neither passed
    nor filtered, so the legend's denominator is 2 and not 3 — reporting it as
    filtered would overstate how many designs the Sc gate rejected."""
    monkeypatch.chdir(tmp_path)
    designs = [rosetta_design(0, sc_value=0.40),
               rosetta_design(1, sc_value=None),
               rosetta_design(2, sc_value=0.70)]
    with figures_kept() as figs:
        rp.plot_sc_histogram(designs, 0.50)
        assert _legend_texts(figs()[-1].axes[0])[0] == (
            "Sc < 0.50 (filtered, n=1 of 2)")


@pytest.mark.logic
def test_every_design_filtered_reports_n_equal_to_the_total(
        tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    designs = [rosetta_design(i, sc_value=0.10) for i in range(4)]
    with figures_kept() as figs:
        rp.plot_sc_histogram(designs, 0.62)
        assert _legend_texts(figs()[-1].axes[0])[0] == (
            "Sc < 0.62 (filtered, n=4 of 4)")


# ── The shaded region and the axis limits ────────────────────────────

@pytest.mark.logic
def test_the_shaded_region_runs_flush_from_the_left_spine_to_the_threshold(
        tmp_path, monkeypatch):
    """`set_xlim` is called before `axvspan` so the span's left edge is the
    explicit left limit rather than whatever autoscaling settles on. The
    failure this guards against is cosmetic but misleading: a white strip at
    the left of the plot reads as "no designs down here" when it is really
    "the shading started late".
    """
    monkeypatch.chdir(tmp_path)
    designs = [rosetta_design(i, sc_value=v) for i, v in enumerate((0.40, 0.60))]
    with figures_kept() as figs:
        rp.plot_sc_histogram(designs, 0.50)
        ax = figs()[-1].axes[0]
        left, right = ax.get_xlim()
        span_left, span_right = _x_extent(
            _patches_coloured(ax, rp.COLOUR_WARN)[0])
    assert span_left == pytest.approx(left)
    assert span_right == pytest.approx(0.50)
    assert right > 0.50


@pytest.mark.logic
def test_the_x_axis_stretches_to_include_a_threshold_no_design_reaches(
        tmp_path, monkeypatch):
    """Every design at 0.90 against a threshold of 0.50 gives a zero-width
    data range, so the margin falls back to its 0.02 floor and the left limit
    comes from the threshold rather than from the data: min(0.90-0.02,
    0.50-0.02) = 0.48. Without the threshold term the dashed line would sit
    off the axes and the figure would not show the gate at all.
    """
    monkeypatch.chdir(tmp_path)
    designs = [rosetta_design(i, sc_value=0.90) for i in range(3)]
    with figures_kept() as figs:
        rp.plot_sc_histogram(designs, 0.50)
        left, right = figs()[-1].axes[0].get_xlim()
    assert left == pytest.approx(0.48)
    assert right == pytest.approx(0.92)


@pytest.mark.logic
def test_the_x_axis_stretches_right_for_a_threshold_above_every_design(
        tmp_path, monkeypatch):
    """The mirror case, and the one where the margin comes from the data
    rather than from its 0.02 floor. Designs at 0.10 and 0.50 span 0.40, so
    the margin is 0.04: the limits are 0.10 - 0.04 on the left and
    0.90 + 0.04 on the right, the right one driven by the threshold.
    """
    monkeypatch.chdir(tmp_path)
    designs = [rosetta_design(i, sc_value=v) for i, v in enumerate((0.10, 0.50))]
    with figures_kept() as figs:
        rp.plot_sc_histogram(designs, 0.90)
        left, right = figs()[-1].axes[0].get_xlim()
    assert left == pytest.approx(0.06)
    assert right == pytest.approx(0.94)


@pytest.mark.logic
def test_the_y_axis_is_padded_above_the_tallest_bar(tmp_path, monkeypatch):
    """The legend sits in the upper right, so the top limit is lifted 20%
    above the autoscaled top. Four identical designs give a single bar of
    height 4, and the top must therefore clear 4."""
    monkeypatch.chdir(tmp_path)
    designs = [rosetta_design(i, sc_value=0.55) for i in range(4)]
    with figures_kept() as figs:
        rp.plot_sc_histogram(designs, 0.50)
        ax = figs()[-1].axes[0]
        tallest = max(b.get_height() for b in _bars(ax))
        top = ax.get_ylim()[1]
    assert tallest == 4
    assert top > 4


# ── Bin selection ────────────────────────────────────────────────────

@pytest.mark.logic
def test_under_twelve_designs_uses_a_fixed_bin_count(tmp_path, monkeypatch):
    """Below twelve there is not enough data for Freedman-Diaconis, so the
    count is max(5, min(n, 10)): three designs give 5 bins, not 3."""
    monkeypatch.chdir(tmp_path)
    designs = [rosetta_design(i, sc_value=v)
               for i, v in enumerate((0.40, 0.50, 0.60))]
    with figures_kept() as figs:
        rp.plot_sc_histogram(designs, 0.50)
        assert len(_bars(figs()[-1].axes[0])) == 5


@pytest.mark.logic
def test_eleven_designs_is_still_the_fixed_count_and_caps_at_ten(
        tmp_path, monkeypatch):
    """Eleven is the last count below the Freedman-Diaconis switch, and
    min(n, 10) caps it, so eleven designs give 10 bins."""
    monkeypatch.chdir(tmp_path)
    designs = [rosetta_design(i, sc_value=0.40 + 0.01 * i) for i in range(11)]
    with figures_kept() as figs:
        rp.plot_sc_histogram(designs, 0.50)
        assert len(_bars(figs()[-1].axes[0])) == 10


@pytest.mark.logic
def test_twelve_designs_switches_to_freedman_diaconis(tmp_path, monkeypatch):
    """Twelve is the boundary: `n_total >= 12` is non-strict, so twelve
    designs take the Freedman-Diaconis branch. On 0.40-0.51 in steps of 0.01
    the IQR is 0.055, the bin width 2*0.055/12^(1/3) = 0.0480 and the raw
    count ceil(0.11/0.0480) = 3 — below the floor, so the answer is 8. That
    the floor rather than the formula decides it here is the point: a narrow
    distribution would otherwise be drawn as three fat bars.
    """
    monkeypatch.chdir(tmp_path)
    designs = [rosetta_design(i, sc_value=0.40 + 0.01 * i) for i in range(12)]
    with figures_kept() as figs:
        rp.plot_sc_histogram(designs, 0.50)
        assert len(_bars(figs()[-1].axes[0])) == 8


@pytest.mark.logic
def test_freedman_diaconis_bin_count_is_used_when_it_clears_the_floor(
        tmp_path, monkeypatch):
    """Ten designs clustered at 0.48-0.57 with outliers at 0.0 and 1.0: the
    IQR stays at 0.055 while the range grows to 1.0, so the raw count is
    ceil(1.0/0.0480) = 21 and the formula, not the clamp, decides."""
    monkeypatch.chdir(tmp_path)
    values = [0.0] + [0.48 + 0.01 * i for i in range(10)] + [1.0]
    designs = [rosetta_design(i, sc_value=v) for i, v in enumerate(values)]
    with figures_kept() as figs:
        rp.plot_sc_histogram(designs, 0.50)
        assert len(_bars(figs()[-1].axes[0])) == 21


@pytest.mark.logic
def test_freedman_diaconis_is_capped_at_forty_bins(tmp_path, monkeypatch):
    """The same shape with the cluster ten times tighter asks for 209 bins
    over twelve designs, which would be one bar per design and a lot of empty
    air. The cap holds it at 40."""
    monkeypatch.chdir(tmp_path)
    values = [0.0] + [0.500 + 0.001 * i for i in range(10)] + [1.0]
    designs = [rosetta_design(i, sc_value=v) for i, v in enumerate(values)]
    with figures_kept() as figs:
        rp.plot_sc_histogram(designs, 0.50)
        assert len(_bars(figs()[-1].axes[0])) == 40


@pytest.mark.logic
def test_twelve_identical_designs_fall_back_to_twelve_bins(
        tmp_path, monkeypatch):
    """A zero IQR would divide by zero in the Freedman-Diaconis width, so the
    `iqr > 0` guard sends it to a flat 12. Twelve designs all at 0.55 is the
    input that gets there."""
    monkeypatch.chdir(tmp_path)
    designs = [rosetta_design(i, sc_value=0.55) for i in range(12)]
    with figures_kept() as figs:
        rp.plot_sc_histogram(designs, 0.50)
        assert len(_bars(figs()[-1].axes[0])) == 12


# ── Degenerate input ─────────────────────────────────────────────────

@pytest.mark.logic
def test_a_single_design_still_renders(tmp_path, monkeypatch):
    """One design is a degenerate histogram — zero-width data range, a single
    bar — but ROSETTA_FILTER_PLOTS publishes a fixed figure list, so it has to
    produce a figure rather than a stack trace."""
    monkeypatch.chdir(tmp_path)
    with figures_kept() as figs:
        rp.plot_sc_histogram([rosetta_design(0, sc_value=0.62)], 0.62)
        ax = figs()[-1].axes[0]
        assert _legend_texts(ax)[0] == "Sc < 0.62 (filtered, n=0 of 1)"
        assert len(_bars(ax)) == 5
    assert_png_not_blank(tmp_path / HISTOGRAM)


@pytest.mark.logic
def test_no_sc_values_at_all_writes_the_placeholder(tmp_path, monkeypatch):
    """Rosetta failed on every design. The placeholder carries the reason, and
    it still has to be a real PNG because the stage's output glob is
    rosetta_*.png."""
    monkeypatch.chdir(tmp_path)
    designs = [rosetta_design(i, sc_value=None) for i in range(3)]
    with figures_kept() as figs:
        rp.plot_sc_histogram(designs, 0.50)
        texts = [t.get_text() for t in figs()[-1].axes[0].texts]
    assert texts == ["No Sc data available"]
    assert_png_not_blank(tmp_path / HISTOGRAM)


@pytest.mark.logic
def test_an_empty_design_list_writes_the_placeholder(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    rp.plot_sc_histogram([], 0.50)
    assert_png_not_blank(tmp_path / HISTOGRAM)


# ── Fallbacks ────────────────────────────────────────────────────────

@pytest.mark.logic
def test_make_empty_plot_writes_the_message_with_no_axes(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with figures_kept() as figs:
        rp.make_empty_plot("nothing here", "custom.png")
        ax = figs()[-1].axes[0]
        assert [t.get_text() for t in ax.texts] == ["nothing here"]
        assert ax.axison is False
    assert_png_not_blank(tmp_path / "custom.png")


@pytest.mark.logic
def test_save_fallback_plots_covers_every_published_figure(
        tmp_path, monkeypatch):
    """modules/rosetta_filtering.nf publishes `rosetta_*.png`, and
    ALL_PLOT_FILES is the single source of truth for what that is. The
    docstring still mentions a second dG_separated plot that the module
    deliberately does not draw, so the list is the thing to assert.
    """
    monkeypatch.chdir(tmp_path)
    rp.save_fallback_plots("No metrics data")
    assert rp.ALL_PLOT_FILES == [HISTOGRAM]
    for name in rp.ALL_PLOT_FILES:
        assert_png_not_blank(tmp_path / name)


# ── CLI ──────────────────────────────────────────────────────────────

@pytest.mark.logic
def test_parse_args_requires_the_metrics_path(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["rosetta_filter_plots.py"])
    with pytest.raises(SystemExit):
        rp.parse_args()


@pytest.mark.logic
def test_main_draws_the_histogram_from_the_file(tmp_path, monkeypatch):
    """End to end on the committed shape: four designs, one of them below the
    0.5 threshold, which is the reference set's own arrangement."""
    monkeypatch.chdir(tmp_path)
    metrics = rosetta_filter_metrics_json(
        tmp_path / "rosetta_filter_metrics.json",
        [rosetta_design(0, sc_value=0.576), rosetta_design(1, sc_value=0.504),
         rosetta_design(14, sc_value=0.542), rosetta_design(59, sc_value=0.433)],
        sc_threshold=0.5,
    )
    monkeypatch.setattr(sys, "argv",
                        ["rosetta_filter_plots.py", "--metrics", str(metrics)])
    with figures_kept() as figs:
        rp.main()
        assert _legend_texts(figs()[-1].axes[0])[0] == (
            "Sc < 0.50 (filtered, n=1 of 4)")
    assert_png_not_blank(tmp_path / HISTOGRAM)


@pytest.mark.logic
def test_main_falls_back_to_a_threshold_of_0_62_when_the_key_is_absent(
        tmp_path, monkeypatch):
    """A metrics file written before sc_threshold existed. The default has to
    be visible in the legend, because a plot drawn against the wrong gate is
    worse than no plot."""
    monkeypatch.chdir(tmp_path)
    metrics = rosetta_filter_metrics_json(
        tmp_path / "metrics.json",
        [rosetta_design(0, sc_value=0.60), rosetta_design(1, sc_value=0.70)],
        sc_threshold=None,
    )
    monkeypatch.setattr(sys, "argv",
                        ["rosetta_filter_plots.py", "--metrics", str(metrics)])
    with figures_kept() as figs:
        rp.main()
        assert _legend_texts(figs()[-1].axes[0])[1] == (
            f"Threshold = {ROSETTA_DEFAULT_SC_THRESHOLD:.2f}")


@pytest.mark.logic
def test_main_with_an_absent_metrics_file_writes_the_fallback(
        tmp_path, monkeypatch):
    """main() exits 0 so the Nextflow stage succeeds, which is only safe
    because the published figure was written first."""
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", ["rosetta_filter_plots.py", "--metrics",
                                      str(tmp_path / "absent.json")])
    with pytest.raises(SystemExit) as exc:
        rp.main()
    assert exc.value.code == 0
    assert_png_not_blank(tmp_path / HISTOGRAM)


@pytest.mark.logic
def test_main_with_a_zero_byte_metrics_file_writes_the_fallback(
        tmp_path, monkeypatch):
    """What an upstream crash leaves behind. The size check comes before
    json.load, so this is a placeholder rather than a JSONDecodeError."""
    monkeypatch.chdir(tmp_path)
    empty = tmp_path / "metrics.json"
    empty.write_text("")
    monkeypatch.setattr(sys, "argv",
                        ["rosetta_filter_plots.py", "--metrics", str(empty)])
    with pytest.raises(SystemExit) as exc:
        rp.main()
    assert exc.value.code == 0
    assert_png_not_blank(tmp_path / HISTOGRAM)


@pytest.mark.logic
def test_main_with_no_designs_writes_the_fallback(tmp_path, monkeypatch):
    """Valid JSON, zero designs: every RFDiffusion backbone failed before
    Rosetta ran."""
    monkeypatch.chdir(tmp_path)
    metrics = rosetta_filter_metrics_json(tmp_path / "metrics.json", [])
    monkeypatch.setattr(sys, "argv",
                        ["rosetta_filter_plots.py", "--metrics", str(metrics)])
    with pytest.raises(SystemExit) as exc:
        rp.main()
    assert exc.value.code == 0
    assert_png_not_blank(tmp_path / HISTOGRAM)


@pytest.mark.logic
def test_main_with_a_designs_key_missing_altogether_writes_the_fallback(
        tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    metrics = tmp_path / "metrics.json"
    metrics.write_text(json.dumps({"sc_threshold": 0.5}))
    monkeypatch.setattr(sys, "argv",
                        ["rosetta_filter_plots.py", "--metrics", str(metrics)])
    with pytest.raises(SystemExit) as exc:
        rp.main()
    assert exc.value.code == 0
    assert_png_not_blank(tmp_path / HISTOGRAM)


@pytest.mark.logic
def test_main_without_matplotlib_writes_placeholder_files(
        tmp_path, monkeypatch):
    """When the matplotlib import failed, main() still creates each published
    path so the stage does not fail on a missing output. The files carry a
    text marker rather than PNG bytes, which Nextflow accepts because it
    matches on the filename, but anything that opens them as images will not.
    """
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(rp, "HAS_MPL", False)
    monkeypatch.setattr(sys, "argv", ["rosetta_filter_plots.py", "--metrics",
                                      str(tmp_path / "metrics.json")])
    with pytest.raises(SystemExit) as exc:
        rp.main()
    assert exc.value.code == 0
    for name in rp.ALL_PLOT_FILES:
        assert (tmp_path / name).read_text() == "matplotlib not available"
