"""Drive bin/rfdiffusion_plots.py against constructed rfdiffusion_metrics.json.

Five figures, and four of them order their rows before drawing: the contact map
and the clustering heatmaps use a dendrogram leaf order, the specificity panel
sorts by contact specificity, and the COM-displacement bars sort by the largest
region displacement. None of those orders is recorded anywhere but in the
figure, so each is asserted on an input whose answer is arithmetic.

The fixtures build the real metrics shape, including the part that trips people
up: rfdiffusion_filter numbers de novo residues from a synthetic base above
every real receptor resnum, so sorting the design residues cannot recover the
region boundaries. That is what per_region_lengths is for, and the legacy
sort-based fallback is tested here to show what it does instead.

Design-region coordinates are pure translations of one another, so the pairwise
design-region RMSD between two designs is exactly the difference of their
offsets. The clustering distance matrix is then checkable by hand.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[2]
BIN = REPO / "bin"
HERE = Path(__file__).resolve().parent
for _p in (str(BIN), str(HERE)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import rfdiffusion_plots as rp  # noqa: E402
from fixtures_plots import (  # noqa: E402
    DE_NOVO_BASE,
    assert_png_not_blank,
    figures_kept,
    receptor_layout,
    rfdiff_design,
    rfdiffusion_metrics_json,
)
from matplotlib import pyplot as plt  # noqa: E402
from matplotlib.colors import to_rgba  # noqa: E402

pytestmark = pytest.mark.logic


def _ticklabels(ax, axis="x"):
    getter = ax.get_xticklabels if axis == "x" else ax.get_yticklabels
    return [t.get_text() for t in getter()]


def _legend_texts(ax):
    legend = ax.get_legend()
    return [] if legend is None else [t.get_text() for t in legend.get_texts()]


def _clustering_panels(fig):
    """The clustering heatmaps, in panel order.

    Each panel contributes three axes in creation order — heatmap, dendrogram
    inset, colourbar — and only the heatmap carries a title. The title is set
    with loc="left", so get_title() with no argument returns the empty centre
    title and would find nothing.
    """
    return [ax for ax in fig.axes if ax.get_title(loc="left")]


# A two-region layout: 5 fixed, 4 de novo, 3 fixed, 2 de novo, 5 fixed.
# receptor_position_order is 1-5, 1000-1003, 6-8, 1004-1005, 9-13.
TWO_REGION = {"fixed_runs": (5, 3, 5), "region_lengths": (4, 2)}


# ── Pure helpers ─────────────────────────────────────────────────────

def test_design_label_strips_the_prefix_and_suffix():
    assert rp._design_label({"design": "design_12.pdb"}) == "12"


def test_get_per_design_residues_prefers_the_per_design_list():
    """Design regions are per-design because contigs can be variable-length, so
    the global list is only a fallback for legacy metrics."""
    d = {"per_design_design_residues": [7, 8]}
    assert rp._get_per_design_residues(d, {99}) == {7, 8}


@pytest.mark.parametrize("value", [None, []])
def test_get_per_design_residues_falls_back_to_global(value):
    d = {} if value is None else {"per_design_design_residues": value}
    assert rp._get_per_design_residues(d, {99}) == {99}


def test_split_design_into_regions_uses_per_region_lengths():
    """Six de novo residues numbered 1000-1005 split 4 then 2, which is what the
    contig asked for. The sort-based route cannot see that boundary."""
    regions = rp._split_design_into_regions(
        set(range(1000, 1006)), {1, 2, 3}, [4, 2])
    assert regions == [set(range(1000, 1004)), {1004, 1005}]


def test_split_design_into_regions_absorbs_a_length_shortfall():
    """per_region_lengths summing below the residue count means the metrics and
    the PDB disagree. The leftover is attached to the last region rather than
    dropped, so no residue goes missing from the coverage denominator."""
    regions = rp._split_design_into_regions(set(range(1000, 1006)), set(), [2, 2])
    assert regions == [{1000, 1001}, {1002, 1003, 1004, 1005}]


def test_split_design_into_regions_legacy_route_collapses_synthetic_resnums():
    """The documented failure of the sort-based fallback. De novo resnums all
    sort after every fixed resnum, so the two regions read as one contiguous run
    of six. This is exactly why plot_design_lengths refuses the fallback and
    warns instead of reporting a wrong region size.
    """
    regions = rp._split_design_into_regions(
        set(range(1000, 1006)), {1, 2, 3, 4, 5}, None)
    assert regions == [set(range(1000, 1006))]


def test_split_design_into_regions_legacy_route_splits_interleaved_resnums():
    """With real interleaved numbering, which is what legacy metrics carried, the
    fallback does separate the regions: 6-7 and 11-12 are divided by 8-10."""
    regions = rp._split_design_into_regions({6, 7, 11, 12},
                                            {1, 2, 3, 4, 5, 8, 9, 10}, None)
    assert regions == [{6, 7}, {11, 12}]


def test_get_per_region_lengths_normalises_absent_and_empty():
    assert rp._get_per_region_lengths({"per_region_lengths": [3, 4]}) == [3, 4]
    assert rp._get_per_region_lengths({"per_region_lengths": []}) == []
    assert rp._get_per_region_lengths({}) == []


def test_label_with_region_len_names_the_single_region():
    d = rfdiff_design(3, **TWO_REGION)
    assert rp._label_with_region_len(d, 0, set()) == "d3 (4)"
    assert rp._label_with_region_len(d, 1, set()) == "d3 (2)"


def test_label_with_region_len_falls_back_to_the_total_design_length():
    """Without per_region_lengths the label can only report the whole design
    region, so a two-region design reads as its total of six."""
    d = rfdiff_design(3, **TWO_REGION, include_per_region_lengths=False)
    assert rp._label_with_region_len(d, 0, set()) == "d3 (6)"


def test_n_regions_reads_the_coordinate_blocks():
    assert rp._n_regions([rfdiff_design(0, **TWO_REGION)]) == 2
    assert rp._n_regions([rfdiff_design(0, region_lengths=(4,))]) == 1
    assert rp._n_regions([rfdiff_design(0, include_coords=False)]) == 1


def test_build_design_labels_joins_the_per_region_lengths():
    labels = rp._build_design_labels([rfdiff_design(7, **TWO_REGION)], set())
    assert labels == ["d7 (4+2)"]


def test_build_design_labels_falls_back_to_a_single_total():
    labels = rp._build_design_labels(
        [rfdiff_design(7, **TWO_REGION, include_per_region_lengths=False)], set())
    assert labels == ["d7 (6)"]


def test_measure_label_width_grows_with_the_longest_label():
    """Used to calibrate the dendrogram inset, so the only property that matters
    is that a longer label measures wider and an empty list measures zero."""
    assert rp._measure_label_width_inches([], 8) == 0.0
    short = rp._measure_label_width_inches(["d1 (4)"], 8)
    long = rp._measure_label_width_inches(["d1 (4)", "d123 (44+22+11)"], 8)
    assert 0 < short < long


# ── Clustering and ordering ──────────────────────────────────────────

def test_per_region_clustering_rmsd_is_the_coordinate_offset_difference():
    """The fixtures translate each design's region along x, so the design-region
    Cα RMSD between two designs is the difference of their offsets. Offsets 0, 2
    and 5 therefore give a matrix of 2, 5 and 3.
    """
    designs = [rfdiff_design(i, **TWO_REGION, coord_offset=off)
               for i, off in enumerate((0.0, 2.0, 5.0))]
    rd = rp._per_region_clustering_data(designs, 0)
    m = rd["rmsd_matrix"]
    assert m[0, 1] == pytest.approx(2.0)
    assert m[0, 2] == pytest.approx(5.0)
    assert m[1, 2] == pytest.approx(3.0)
    assert np.allclose(m, m.T)
    assert np.allclose(np.diag(m), 0.0)


def test_per_region_clustering_needs_two_designs_with_coords():
    designs = [rfdiff_design(0, **TWO_REGION)]
    assert rp._per_region_clustering_data(designs, 0) is None


def test_per_region_clustering_returns_none_for_an_absent_region():
    """A single-region metrics set has no region 2, so panel 2 is not plottable
    rather than being drawn from whatever coords happen to be first."""
    designs = [rfdiff_design(i, region_lengths=(4,)) for i in range(3)]
    assert rp._per_region_clustering_data(designs, 1) is None


def test_per_region_clustering_region_zero_falls_back_to_the_flat_coords():
    """Legacy metrics carry design_region_coords but no per_region_coords, and
    region 0 reads from it so the single-region case still clusters."""
    designs = []
    for i, off in enumerate((0.0, 3.0)):
        d = rfdiff_design(i, region_lengths=(4,), coord_offset=off)
        d.pop("per_region_coords")
        designs.append(d)
    rd = rp._per_region_clustering_data(designs, 0)
    assert rd["rmsd_matrix"][0, 1] == pytest.approx(3.0)


def test_dendrogram_order_puts_the_outlier_at_an_end():
    """Designs 0 and 1 sit 0.5 Å apart and design 2 sits 50 Å away, so the two
    near-identical designs must end up adjacent and the outlier at one edge.
    Asserting the exact permutation would pin scipy's tie-breaking rather than
    the clustering, so the two structural facts are asserted instead.
    """
    designs = [rfdiff_design(i, **TWO_REGION, coord_offset=off)
               for i, off in enumerate((0.0, 0.5, 50.0))]
    order, Z, m = rp._get_dendrogram_order(designs)
    assert sorted(order) == [0, 1, 2]
    assert abs(order.index(0) - order.index(1)) == 1
    assert order.index(2) in (0, 2)
    assert Z is not None
    assert m[0, 1] == pytest.approx(0.5)


def test_dendrogram_order_is_the_identity_below_two_designs():
    """One design cannot be clustered, so the order is the input order and the
    linkage is absent, which is what suppresses the dendrogram inset."""
    order, Z, m = rp._get_dendrogram_order([rfdiff_design(0, **TWO_REGION)])
    assert (order, Z, m) == ([0], None, None)


def test_dendrogram_order_is_the_identity_without_coords():
    designs = [rfdiff_design(i, include_coords=False) for i in range(3)]
    assert rp._get_dendrogram_order(designs) == ([0, 1, 2], None, None)


# ── Contact-map data ─────────────────────────────────────────────────

def test_build_contact_map_data_splits_contacts_by_design_region():
    """A contact's colour depends only on whether its resnum is a de novo one.
    Here resnum 3 is fixed and 1000 is de novo, so one contact is 'out' and the
    other 'in', and their x positions are their indices in
    receptor_position_order: position 2 and position 5.
    """
    order, design_res, _ = receptor_layout(**TWO_REGION)
    assert order[2] == 3 and order[5] == DE_NOVO_BASE
    d = rfdiff_design(0, **TWO_REGION, contact_resnums=(3, DE_NOVO_BASE))
    data = rp._build_contact_map_data([d], set(), set())
    assert data["design_contact_positions"][0] == {2: "out", 5: "in"}
    assert data["design_region_positions"][0] == {5, 6, 7, 8, 12, 13}
    assert data["design_n_res"] == [len(order)]


def test_build_contact_map_data_ignores_contacts_outside_the_position_order():
    """A resnum absent from receptor_position_order cannot be placed on the
    x-axis, so it is dropped rather than plotted at a guessed position."""
    d = rfdiff_design(0, **TWO_REGION, contact_resnums=(3, 9999))
    data = rp._build_contact_map_data([d], set(), set())
    assert data["design_contact_positions"][0] == {2: "out"}


def test_build_contact_map_data_falls_back_to_sorted_resnums():
    """Without receptor_position_order the axis is the sorted union of fixed and
    design residues, which for synthetic de novo numbering puts the whole design
    region at the right-hand end."""
    d = rfdiff_design(0, **TWO_REGION, contact_resnums=(3,))
    d.pop("receptor_position_order")
    data = rp._build_contact_map_data([d], set(), {1, 2, 3})
    assert data["design_n_res"] == [9]            # 3 fixed + 6 de novo
    assert data["design_contact_positions"][0] == {2: "out"}
    assert data["design_region_positions"][0] == {3, 4, 5, 6, 7, 8}


def test_build_contact_map_data_returns_designs_in_dendrogram_order():
    designs = [rfdiff_design(i, **TWO_REGION, coord_offset=off)
               for i, off in enumerate((0.0, 0.5, 50.0))]
    data = rp._build_contact_map_data(designs, set(), set())
    names = [d["design"] for d in data["sorted_designs"]]
    assert sorted(names) == ["design_0.pdb", "design_1.pdb", "design_2.pdb"]
    assert abs(names.index("design_0.pdb") - names.index("design_1.pdb")) == 1


def test_build_effector_freq_is_the_fraction_of_designs_contacting_each_residue():
    """Residue 95 is contacted by both designs, 97 by one of two, and 96 by
    neither — but 96 still gets a bar at 0 because the range is contiguous from
    the lowest to the highest contacted residue.
    """
    designs = [
        rfdiff_design(0, effector_contact_residues=(95, 97)),
        rfdiff_design(1, effector_contact_residues=(95,)),
    ]
    eff = rp._build_effector_freq(designs, 2)
    assert eff["range"] == [95, 96, 97]
    assert eff["n_eff"] == 3
    assert list(eff["freq"]) == [1.0, 0.0, 0.5]


def test_build_effector_freq_is_none_without_effector_contacts():
    """No effector contacts means no frequency bar below the contact map, which
    is what shrinks the figure rather than drawing an empty strip."""
    assert rp._build_effector_freq([rfdiff_design(0)], 1) is None


def test_set_contact_xticks_uses_the_longest_design_as_the_reference():
    """Designs have different lengths with a variable-length contig, so the axis
    has to come from the longest one or the shorter designs' hatched tails would
    sit past the last tick. The longest here is design 1 at 19 positions, and the
    final tick label is its position count.
    """
    short = rfdiff_design(0, fixed_runs=(5, 5), region_lengths=(2,))
    long = rfdiff_design(1, fixed_runs=(5, 5), region_lengths=(9,))
    data = rp._build_contact_map_data([short, long], set(), set())
    fig, ax = plt.subplots()
    try:
        rp._set_contact_xticks(ax, data, max(data["design_n_res"]), set(), set())
        assert _ticklabels(ax)[-1] == "19"
        assert ax.get_xlabel() == "Receptor position"
    finally:
        plt.close(fig)


def test_contact_legend_has_four_entries():
    labels = [h.get_label() for h in rp._contact_legend_handles()]
    assert labels == ["Contact in design region", "Contact outside design region",
                      "Design region", "Beyond protein end"]


# ── Plot 1: contact map ──────────────────────────────────────────────

def test_contact_map_labels_rows_in_dendrogram_order(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    designs = [rfdiff_design(i, **TWO_REGION, coord_offset=off,
                             contact_resnums=(3, DE_NOVO_BASE),
                             effector_contact_residues=(95, 96))
               for i, off in enumerate((0.0, 0.5, 50.0))]
    with figures_kept() as figs:
        rp.plot_contact_map(designs, set(), set())
        ax_heat = figs()[-1].axes[0]
        labels = _ticklabels(ax_heat, "y")
    assert sorted(labels) == ["d0 (4+2)", "d1 (4+2)", "d2 (4+2)"]
    assert abs(labels.index("d0 (4+2)") - labels.index("d1 (4+2)")) == 1
    assert_png_not_blank(tmp_path / "rfdiff_contact_map.png")


def test_contact_map_hatches_the_tail_of_a_shorter_design(tmp_path, monkeypatch):
    """Variable-length contigs give designs different receptor lengths. The short
    design's missing positions are drawn as a hatched 'beyond protein end' block
    so an absent contact is not read as a failed one."""
    monkeypatch.chdir(tmp_path)
    designs = [
        rfdiff_design(0, fixed_runs=(5, 5), region_lengths=(2,), coord_offset=0.0),
        rfdiff_design(1, fixed_runs=(5, 5), region_lengths=(9,), coord_offset=4.0),
    ]
    with figures_kept() as figs:
        rp.plot_contact_map(designs, set(), set())
        hatched = [p for p in figs()[-1].axes[0].patches if p.get_hatch() == "///"]
    assert len(hatched) == 1
    assert hatched[0].get_width() == pytest.approx(7.0)   # 19 - 12 positions
    assert_png_not_blank(tmp_path / "rfdiff_contact_map.png")


def test_contact_map_without_an_effector_bar_drops_that_axis(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    designs = [rfdiff_design(i, **TWO_REGION, coord_offset=float(i))
               for i in range(2)]
    with figures_kept() as figs:
        rp.plot_contact_map(designs, set(), set())
        fig = figs()[-1]
    # heatmap + dendrogram inset only; no frequency bar.
    assert len(fig.axes) == 2
    assert_png_not_blank(tmp_path / "rfdiff_contact_map.png")


def test_contact_map_single_design_has_no_dendrogram(tmp_path, monkeypatch):
    """One design cannot be clustered, so the inset is suppressed and the figure
    is the heatmap alone."""
    monkeypatch.chdir(tmp_path)
    with figures_kept() as figs:
        rp.plot_contact_map([rfdiff_design(0, **TWO_REGION)], set(), set())
        assert len(figs()[-1].axes) == 1
    assert_png_not_blank(tmp_path / "rfdiff_contact_map.png")


def test_contact_map_with_no_designs_still_writes_a_file(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    rp.plot_contact_map([], set(), set())
    assert_png_not_blank(tmp_path / "rfdiff_contact_map.png")


# ── Plot 2: clustering heatmaps ──────────────────────────────────────

def test_design_clustering_draws_one_panel_per_region(tmp_path, monkeypatch):
    """Two regions in the metrics means two heatmaps, each with its own
    dendrogram and colourbar, so a region that diversified and a region that did
    not can be told apart."""
    monkeypatch.chdir(tmp_path)
    designs = [rfdiff_design(i, **TWO_REGION, coord_offset=off)
               for i, off in enumerate((0.0, 1.0, 4.0))]
    with figures_kept() as figs:
        rp.plot_design_clustering(designs, set())
        fig = figs()[-1]
        heatmaps = _clustering_panels(fig)
    assert [ax.get_title(loc="left") for ax in heatmaps] == ["Design region 1",
                                                              "Design region 2"]
    assert len(fig.axes) == 6          # 2 x (heatmap + dendrogram + colourbar)
    assert_png_not_blank(tmp_path / "rfdiff_design_clustering.png")


def test_design_clustering_titles_a_single_region_without_a_number(
        tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    designs = [rfdiff_design(i, region_lengths=(4,), coord_offset=float(i))
               for i in range(2)]
    with figures_kept() as figs:
        rp.plot_design_clustering(designs, set())
        titles = [ax.get_title(loc="left")
                  for ax in _clustering_panels(figs()[-1])]
    assert titles == ["Design region"]


def test_design_clustering_colours_labels_by_filter_outcome(tmp_path, monkeypatch):
    """Row and column labels are recoloured to carry pass/fail, which is the only
    place the clustering plot shows the filter verdict."""
    monkeypatch.chdir(tmp_path)
    designs = [
        rfdiff_design(0, region_lengths=(4,), coord_offset=0.0, passes_filter=True),
        rfdiff_design(1, region_lengths=(4,), coord_offset=3.0, passes_filter=False),
    ]
    with figures_kept() as figs:
        rp.plot_design_clustering(designs, set())
        heat = _clustering_panels(figs()[-1])[0]
        by_label = {t.get_text(): t.get_color() for t in heat.get_yticklabels()}
    assert by_label["d0 (4)"] == rp.COLOUR_PASS
    assert by_label["d1 (4)"] == rp.COLOUR_FAIL


def test_design_clustering_annotates_cells_only_up_to_fifteen_designs(
        tmp_path, monkeypatch):
    """A 16x16 grid of numbers is unreadable, so the in-cell RMSD text is drawn
    only at 15 designs or fewer. The comparison is non-strict: exactly 15 is
    still annotated."""
    monkeypatch.chdir(tmp_path)
    for n, expect_text in ((15, True), (16, False)):
        designs = [rfdiff_design(i, region_lengths=(4,), coord_offset=float(i))
                   for i in range(n)]
        with figures_kept() as figs:
            rp.plot_design_clustering(designs, set())
            heat = _clustering_panels(figs()[-1])[0]
            assert bool(len(heat.texts)) is expect_text


def test_design_clustering_with_one_design_still_writes_a_file(
        tmp_path, monkeypatch):
    """A single design has nothing to cluster against. The placeholder names the
    reason rather than leaving a blank heatmap."""
    monkeypatch.chdir(tmp_path)
    rp.plot_design_clustering([rfdiff_design(0, **TWO_REGION)], set())
    assert_png_not_blank(tmp_path / "rfdiff_design_clustering.png")


def test_design_clustering_with_no_designs_still_writes_a_file(
        tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    rp.plot_design_clustering([], set())
    assert_png_not_blank(tmp_path / "rfdiff_design_clustering.png")


# ── Plot 3: specificity and coverage ─────────────────────────────────

def test_specificity_ranks_designs_by_descending_contact_specificity(
        tmp_path, monkeypatch):
    """Design 1 is the most specific at 0.9 and design 0 the least at 0.2, so the
    axis must read 1, 2, 0. Input order 0, 1, 2 is the wrong answer.

    The two panels share an x-axis, so the visible tick labels belong to the
    lower one; the upper panel's labels are set and then hidden.
    """
    monkeypatch.chdir(tmp_path)
    designs = [
        rfdiff_design(i, **TWO_REGION, frac_contacts_in_design=f,
                      n_contacts_in_design=int(f * 10), n_contact_pairs=10)
        for i, f in enumerate((0.2, 0.9, 0.5))
    ]
    with figures_kept() as figs:
        rp.plot_specificity_coverage(designs, set(), 0.5, set())
        assert _ticklabels(figs()[-1].axes[1]) == ["1", "2", "0"]
    assert_png_not_blank(tmp_path / "rfdiff_specificity_coverage.png")


def test_specificity_bars_are_coloured_by_the_filter_verdict(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    designs = [
        rfdiff_design(0, **TWO_REGION, frac_contacts_in_design=0.9,
                      n_contact_pairs=10, n_contacts_in_design=9,
                      passes_filter=True),
        rfdiff_design(1, **TWO_REGION, frac_contacts_in_design=0.1,
                      n_contact_pairs=10, n_contacts_in_design=1,
                      passes_filter=False),
    ]
    with figures_kept() as figs:
        rp.plot_specificity_coverage(designs, set(), 0.5, set())
        bars = figs()[-1].axes[0].patches[:2]
        colours = [tuple(b.get_facecolor()) for b in bars]
    assert colours[0] == pytest.approx(to_rgba(rp.COLOUR_PASS))
    assert colours[1] == pytest.approx(to_rgba(rp.COLOUR_FAIL))


@pytest.mark.parametrize(("threshold", "expect_line"),
                         [(0.0, False), (0.0001, True), (0.5, True)])
def test_specificity_threshold_line_appears_only_above_zero(
        tmp_path, monkeypatch, threshold, expect_line):
    """The comparison is strict: `threshold > 0`. A min_hotspot_frac of exactly 0
    means no hotspot requirement was set, and the committed reference metrics
    carry 0, so the no-line case is the common one. Any positive value, however
    small, draws the dashed line and adds its legend entry.
    """
    monkeypatch.chdir(tmp_path)
    designs = [rfdiff_design(i, **TWO_REGION, frac_contacts_in_design=0.5,
                             n_contact_pairs=2, n_contacts_in_design=1)
               for i in range(2)]
    with figures_kept() as figs:
        rp.plot_specificity_coverage(designs, set(), threshold, set())
        texts = _legend_texts(figs()[-1].axes[0])
    assert (len(texts) == 3) is expect_line
    assert texts[:2] == ["Pass", "Filtered"]


def test_specificity_coverage_is_contacted_residues_over_region_size(
        tmp_path, monkeypatch):
    """Region 1 holds four de novo residues and two are contacted, so its bar is
    0.5. Region 2 holds two and both are contacted, so its bar is 1.0. The bars
    come from separate containers, one per region.
    """
    monkeypatch.chdir(tmp_path)
    contacts = (DE_NOVO_BASE, DE_NOVO_BASE + 1,          # 2 of region 1's 4
                DE_NOVO_BASE + 4, DE_NOVO_BASE + 5)      # 2 of region 2's 2
    d = rfdiff_design(0, **TWO_REGION, contact_resnums=contacts)
    with figures_kept() as figs:
        rp.plot_specificity_coverage([d], set(), 0.0, set())
        ax_cov = figs()[-1].axes[1]
        heights = [c.patches[0].get_height() for c in ax_cov.containers]
    assert heights == [pytest.approx(0.5), pytest.approx(1.0)]
    assert "Region 1" in _legend_texts(ax_cov)


def test_specificity_coverage_of_an_uncontacted_region_is_zero(
        tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    d = rfdiff_design(0, **TWO_REGION, contact_resnums=(3,))
    with figures_kept() as figs:
        rp.plot_specificity_coverage([d], set(), 0.0, set())
        heights = [c.patches[0].get_height()
                   for c in figs()[-1].axes[1].containers]
    assert heights == [0.0, 0.0]


def test_specificity_single_region_has_no_coverage_legend(tmp_path, monkeypatch):
    """With one region there is nothing to distinguish, so the Region legend is
    suppressed."""
    monkeypatch.chdir(tmp_path)
    d = rfdiff_design(0, region_lengths=(4,), contact_resnums=(DE_NOVO_BASE,))
    with figures_kept() as figs:
        rp.plot_specificity_coverage([d], set(), 0.0, set())
        assert _legend_texts(figs()[-1].axes[1]) == []
    assert_png_not_blank(tmp_path / "rfdiff_specificity_coverage.png")


# ── Plot 4: design-region lengths ────────────────────────────────────

def test_design_lengths_draws_one_panel_per_region(tmp_path, monkeypatch):
    """Region 1 lengths are 4 and 5, one design each. Region 2 is 2 for both
    designs, so its single bar has a count of 2."""
    monkeypatch.chdir(tmp_path)
    designs = [
        rfdiff_design(0, fixed_runs=(5, 3, 5), region_lengths=(4, 2)),
        rfdiff_design(1, fixed_runs=(5, 3, 5), region_lengths=(5, 2)),
    ]
    with figures_kept() as figs:
        rp.plot_design_lengths(designs, set(), set())
        panels = figs()[-1].axes
    assert [ax.get_title() for ax in panels] == ["Design region 1",
                                                  "Design region 2"]
    assert [(p.get_x() + p.get_width() / 2, p.get_height())
            for p in panels[0].patches] == [(pytest.approx(4.0), 1),
                                            (pytest.approx(5.0), 1)]
    assert [p.get_height() for p in panels[1].patches] == [2]
    assert_png_not_blank(tmp_path / "rfdiff_design_lengths.png")


def test_design_lengths_ticks_every_integer_in_the_observed_span(
        tmp_path, monkeypatch):
    """Lengths of 4 and 7 with nothing between them still get ticks at 5 and 6,
    so the gap in the distribution is readable. The ticks never run past the
    data, so there is no tick at 3 or 8."""
    monkeypatch.chdir(tmp_path)
    designs = [rfdiff_design(0, region_lengths=(4,)),
               rfdiff_design(1, region_lengths=(7,))]
    with figures_kept() as figs:
        rp.plot_design_lengths(designs, set(), set())
        assert _ticklabels(figs()[-1].axes[0]) == ["4", "5", "6", "7"]


def test_design_lengths_warns_and_skips_designs_without_per_region_lengths(
        tmp_path, monkeypatch, capsys):
    """Refusing the sort-based fallback is deliberate: it would collapse the
    regions and report a wrong size. The design is dropped from the counts and
    the omission is reported, so a count of 1 from two designs is explained."""
    monkeypatch.chdir(tmp_path)
    designs = [
        rfdiff_design(0, region_lengths=(4,)),
        rfdiff_design(1, region_lengths=(4,), include_per_region_lengths=False),
    ]
    with figures_kept() as figs:
        rp.plot_design_lengths(designs, set(), set())
        assert [p.get_height() for p in figs()[-1].axes[0].patches] == [1]
    assert "1/2 designs have no per_region_lengths" in capsys.readouterr().err


def test_design_lengths_with_no_per_region_lengths_at_all_writes_a_placeholder(
        tmp_path, monkeypatch):
    """No-data path. Every design missing the field leaves nothing plottable, and
    the placeholder tells the user to regenerate the metrics."""
    monkeypatch.chdir(tmp_path)
    designs = [rfdiff_design(i, include_per_region_lengths=False)
               for i in range(2)]
    rp.plot_design_lengths(designs, set(), set())
    assert_png_not_blank(tmp_path / "rfdiff_design_lengths.png")


def test_design_lengths_with_no_designs_still_writes_a_file(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    rp.plot_design_lengths([], set(), set())
    assert_png_not_blank(tmp_path / "rfdiff_design_lengths.png")


# ── Plot 5: COM displacement ─────────────────────────────────────────

def test_com_displacement_sorts_designs_by_their_largest_region_shift(
        tmp_path, monkeypatch):
    """Ordering is ascending by the maximum displacement across a design's
    regions, so the most-moved design sits at the right-hand edge. Design 1 peaks
    at 1.0, design 0 at 5.0 and design 2 at 9.0, giving 1, 0, 2 — and design 0's
    second region at 8.0 must not promote it above design 2.
    """
    monkeypatch.chdir(tmp_path)
    designs = [
        rfdiff_design(0, com_displacement=[5.0, 8.0]),
        rfdiff_design(1, com_displacement=[1.0, 0.5]),
        rfdiff_design(2, com_displacement=[9.0, 2.0]),
    ]
    with figures_kept() as figs:
        rp.plot_com_displacement(designs)
        assert _ticklabels(figs()[-1].axes[0]) == ["design_1", "design_0",
                                                    "design_2"]
    assert_png_not_blank(tmp_path / "rfdiff_com_displacement.png")


def test_com_displacement_treats_null_and_nan_as_missing(tmp_path, monkeypatch):
    """A region with no measurable displacement comes through as null or NaN. The
    design is still plotted from its other region, and the missing region draws
    at zero height rather than dropping the whole bar group."""
    monkeypatch.chdir(tmp_path)
    designs = [
        rfdiff_design(0, com_displacement=[float("nan"), 3.0]),
        rfdiff_design(1, com_displacement=[None, 4.0]),
    ]
    with figures_kept() as figs:
        rp.plot_com_displacement(designs)
        ax = figs()[-1].axes[0]
        first_region = [p.get_height() for p in ax.containers[0].patches]
        second_region = [p.get_height() for p in ax.containers[1].patches]
    assert first_region == [0.0, 0.0]
    assert second_region == [3.0, 4.0]


def test_com_displacement_reports_how_many_designs_were_skipped(
        tmp_path, monkeypatch):
    """A design whose every region value is missing cannot be plotted, so it is
    skipped and the count is annotated on the axes. Silently dropping it would
    make the figure understate the cohort."""
    monkeypatch.chdir(tmp_path)
    designs = [
        rfdiff_design(0, com_displacement=[2.0]),
        rfdiff_design(1, com_displacement=[None]),
        rfdiff_design(2),                                   # field absent
    ]
    with figures_kept() as figs:
        rp.plot_com_displacement(designs)
        ax = figs()[-1].axes[0]
        assert _ticklabels(ax) == ["design_0"]
        assert "2 designs skipped (no data)" in [t.get_text() for t in ax.texts]
    assert_png_not_blank(tmp_path / "rfdiff_com_displacement.png")


def test_com_displacement_singular_wording_for_one_skipped_design(
        tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    designs = [rfdiff_design(0, com_displacement=[2.0]),
               rfdiff_design(1, com_displacement=[])]
    with figures_kept() as figs:
        rp.plot_com_displacement(designs)
        texts = [t.get_text() for t in figs()[-1].axes[0].texts]
    assert "1 design skipped (no data)" in texts


def test_com_displacement_single_region_has_no_legend(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with figures_kept() as figs:
        rp.plot_com_displacement([rfdiff_design(0, com_displacement=[2.0])])
        assert _legend_texts(figs()[-1].axes[0]) == []


def test_com_displacement_labels_each_region_when_there_are_several(
        tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with figures_kept() as figs:
        rp.plot_com_displacement([rfdiff_design(0, com_displacement=[2.0, 3.0])])
        assert _legend_texts(figs()[-1].axes[0]) == ["region 1", "region 2"]


def test_com_displacement_with_no_values_still_writes_a_file(tmp_path, monkeypatch):
    """No-data path: every design missing the field, which is what a metrics file
    from before the displacement metric existed looks like."""
    monkeypatch.chdir(tmp_path)
    rp.plot_com_displacement([rfdiff_design(0), rfdiff_design(1)])
    assert_png_not_blank(tmp_path / "rfdiff_com_displacement.png")


def test_com_displacement_with_no_designs_still_writes_a_file(
        tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    rp.plot_com_displacement([])
    assert_png_not_blank(tmp_path / "rfdiff_com_displacement.png")


# ── Fallbacks and CLI ────────────────────────────────────────────────

def test_save_fallback_plots_covers_every_published_figure(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    rp.save_fallback_plots("No metrics data")
    assert rp.ALL_PLOT_FILES == [
        "rfdiff_contact_map.png", "rfdiff_design_clustering.png",
        "rfdiff_specificity_coverage.png", "rfdiff_design_lengths.png",
        "rfdiff_com_displacement.png",
    ]
    for name in rp.ALL_PLOT_FILES:
        assert_png_not_blank(tmp_path / name)


def test_main_writes_all_five_figures(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    order, design_res, fixed_res = receptor_layout(**TWO_REGION)
    designs = [
        rfdiff_design(i, **TWO_REGION, coord_offset=off,
                      contact_resnums=(3, DE_NOVO_BASE, DE_NOVO_BASE + 4),
                      effector_contact_residues=(95, 96, 97),
                      com_displacement=[5.5 + i, 8.9],
                      passes_filter=(i != 2))
        for i, off in enumerate((0.0, 0.5, 6.0))
    ]
    metrics = rfdiffusion_metrics_json(
        tmp_path / "rfdiffusion_metrics.json", designs,
        fixed_residues=fixed_res, design_residues=design_res,
        min_hotspot_frac=0.5)
    monkeypatch.setattr(sys, "argv",
                        ["rfdiffusion_plots.py", "--metrics", str(metrics)])
    rp.main()
    for name in rp.ALL_PLOT_FILES:
        assert_png_not_blank(tmp_path / name)


def test_main_with_a_missing_metrics_file_writes_fallbacks(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv",
                        ["rfdiffusion_plots.py", "--metrics", "absent.json"])
    with pytest.raises(SystemExit) as exc:
        rp.main()
    assert exc.value.code == 0
    for name in rp.ALL_PLOT_FILES:
        assert_png_not_blank(tmp_path / name)


def test_main_with_an_empty_metrics_file_writes_fallbacks(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    empty = tmp_path / "rfdiffusion_metrics.json"
    empty.write_text("")
    monkeypatch.setattr(sys, "argv",
                        ["rfdiffusion_plots.py", "--metrics", str(empty)])
    with pytest.raises(SystemExit) as exc:
        rp.main()
    assert exc.value.code == 0
    for name in rp.ALL_PLOT_FILES:
        assert_png_not_blank(tmp_path / name)


def test_main_with_zero_designs_writes_fallbacks(tmp_path, monkeypatch):
    """A valid metrics file whose design list is empty: every backbone failed
    before metrics were computed. All five published paths still have to exist."""
    monkeypatch.chdir(tmp_path)
    metrics = rfdiffusion_metrics_json(tmp_path / "m.json", [])
    monkeypatch.setattr(sys, "argv",
                        ["rfdiffusion_plots.py", "--metrics", str(metrics)])
    with pytest.raises(SystemExit) as exc:
        rp.main()
    assert exc.value.code == 0
    for name in rp.ALL_PLOT_FILES:
        assert_png_not_blank(tmp_path / name)


def test_main_without_matplotlib_writes_placeholder_files(tmp_path, monkeypatch):
    """When the matplotlib import failed, main() still creates each published
    path, carrying a text marker rather than PNG bytes. Nextflow matches on the
    filename so the stage passes, but the files are not readable as images."""
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(rp, "HAS_MPL", False)
    monkeypatch.setattr(sys, "argv",
                        ["rfdiffusion_plots.py", "--metrics", "absent.json"])
    with pytest.raises(SystemExit) as exc:
        rp.main()
    assert exc.value.code == 0
    for name in rp.ALL_PLOT_FILES:
        assert (tmp_path / name).read_text() == "matplotlib not available"
