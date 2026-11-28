"""Stage 5 (within-sequence figures) — bin/negsteer_within_sequence_plots.py.

Producer: ``NEGSTEER_WITHIN_SEQUENCE_PLOTS`` → ``bin/negsteer_within_sequence_plots.py``,
which reads the per-sequence raw_per_seed_results.csv files and writes the five
per-seed PNGs that Nextflow publishes.

Where the cohort script asks how the sequences compare, this one asks whether a
sequence's cohort median is trustworthy or the work of one lucky seed. So the
assertions here are about the per-seed record: which column family a seed's
final metrics come from, which seeds belong to the representative
sequence_group, and where the cold-start baseline comes from.

``classify_seed`` lives in both this module and bin/negsteer_plots.py as two
separate copies. The parametrised agreement test over every outcome is in
test_negsteer_plots.py, since it has to import both modules either way.

Fixture cohort: tests/characterization/fixtures_cohort.py.
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

# The module selects the Agg backend itself, before importing pyplot, so
# importing it here is what makes this file headless. Nothing above may import
# pyplot.
import negsteer_within_sequence_plots as wsp  # noqa: E402


def assert_drawn(out: Path) -> None:
    """The figure landed and has more than one colour in it."""
    assert out.exists(), f"{out.name} was not written"
    assert out.stat().st_size > 0, f"{out.name} is empty"
    n = count_distinct_colours(out)
    assert n > 1, f"{out.name} is a single flat colour ({n} distinct)"


def seed_rec(*, seed=0, sg="0", stage="steering", outcome="pass",
             ra=2.0, tj=0.70, wj=0.66, steered_ra=None, reverted_ra=None,
             **metrics) -> dict:
    """One record in the shape ``load_cohort_seeds`` hands the plot functions.

    Built directly rather than through the CSV so a plot test can state the one
    value it is about. The loader's own output is pinned separately.
    """
    rec = {
        "seed": seed, "sg": sg, "stage": stage, "outcome": outcome,
        "ra_final": ra, "tj_final": tj, "wj_final": wj,
        "complex_plddt": 0.92, "avg_plddt": 92.0, "interface_plddt": 0.80,
        "iptm": 0.72, "ptm": 0.80, "actifptm": 0.65,
        "ipae": 4.0, "pae_pass_frac": 0.80, "pae_mean": 8.5,
        "ipsae_min": 0.44, "ipsae_ab": 0.46, "ipsae_ba": 0.42,
        "ipsae_min_15": 0.48, "ipsae_ab_15": 0.50, "ipsae_ba_15": 0.46,
        "steered_ra_raw": ra if steered_ra is None else steered_ra,
        "reverted_ra_raw": reverted_ra,
        "n_mut": 0 if stage == "cold_start" else 3,
    }
    rec.update(metrics)
    return rec


def capture_figure(monkeypatch) -> dict:
    """Record what the figure carried at the moment it was saved.

    The plot functions close their figure straight after savefig, so this is
    the only point at which the drawn artists can be read back.
    """
    seen: dict = {}
    real_savefig = wsp.plt.savefig

    def recording_savefig(*args, **kwargs):
        fig = wsp.plt.gcf()
        seen["titles"] = [ax.get_title() for ax in fig.axes]
        seen["xlabels"] = [ax.get_xlabel() for ax in fig.axes]
        seen["texts"] = [t.get_text() for ax in fig.axes for t in ax.texts]
        seen["n_scatters"] = sum(len(ax.collections) for ax in fig.axes)
        seen["n_lines"] = sum(len(ax.lines) for ax in fig.axes)
        return real_savefig(*args, **kwargs)

    monkeypatch.setattr(wsp.plt, "savefig", recording_savefig)
    return seen


# ═══════════════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════════════

@pytest.mark.logic
@pytest.mark.parametrize("name, expected", [
    ("design_0_seq_2", "d0"),
    ("design_27_seq_1", "d27"),
    ("input_control_polyA", "input_control_polyA"),
])
def test_design_id_groups_the_sequences_of_one_backbone(name, expected):
    """The trajectory plot puts one panel per RFDiffusion design, so every MPNN
    sequence of that backbone has to map onto the same panel."""
    assert wsp._design_id(name) == expected


@pytest.mark.logic
@pytest.mark.parametrize("name, expected", [
    ("design_0_seq_3", 3),
    ("design_12_seq_0", 0),
    ("not_a_design", 0),
])
def test_seq_num_is_stable_across_designs(name, expected):
    """Sequence 3 keeps one colour wherever it appears, so the colour is read
    off the sequence number rather than assigned per panel."""
    assert wsp._seq_num(name) == expected


@pytest.mark.logic
@pytest.mark.parametrize("n, expected", [
    (0, (1, 1)),
    (1, (1, 1)),
    (3, (1, 3)),
    (4, (2, 2)),
    (5, (2, 3)),
    (9, (3, 3)),
    (16, (4, 4)),
    (17, (5, 4)),
])
def test_grid_shape_stays_within_four_columns(n, expected):
    assert wsp._grid_shape(n) == expected
    rows, cols = expected
    assert rows * cols >= n


@pytest.mark.logic
def test_make_empty_plot_writes_a_readable_placeholder(tmp_path):
    out = tmp_path / "placeholder.png"
    wsp._make_empty_plot("No per-seed RMSD data", str(out))
    assert_drawn(out)


# ═══════════════════════════════════════════════════════════════════════════
# load_cohort_seeds
# ═══════════════════════════════════════════════════════════════════════════

@pytest.mark.logic
def test_load_cohort_seeds_takes_the_reverted_metrics_as_final(tmp_path):
    """design_0_seq_1's seeds all went through reversion. Their steered
    prediction still reads 2.5 Å; the reverted one, the last made, reads 30 Å.
    The record's final metrics must be the reverted ones, with the steered
    value kept separately for the trajectory plot's middle point."""
    paths = fc.cohort_tree(tmp_path)
    rep, all_seeds, comps, tiers, cold = wsp.load_cohort_seeds(
        paths["runs_dir"], paths["cross_csv"])

    reverted = rep["design_0_seq_1"]
    assert len(reverted) == 3
    assert {s["stage"] for s in reverted} == {"reversion"}
    assert {s["ra_final"] for s in reverted} == {30.0}
    assert {s["steered_ra_raw"] for s in reverted} == {2.5}
    assert {s["reverted_ra_raw"] for s in reverted} == {30.0}


@pytest.mark.logic
def test_load_cohort_seeds_keeps_the_steered_metrics_when_reversion_never_ran(tmp_path):
    paths = fc.cohort_tree(tmp_path)
    rep, _all, _c, _t, _cb = wsp.load_cohort_seeds(paths["runs_dir"],
                                                   paths["cross_csv"])
    held = rep["design_0_seq_0"]
    assert {s["stage"] for s in held} == {"steering"}
    assert {s["ra_final"] for s in held} == {2.0}
    assert {s["tj_final"] for s in held} == {0.70}
    assert {s["reverted_ra_raw"] for s in held} == {None}


@pytest.mark.logic
def test_load_cohort_seeds_separates_the_representative_group(tmp_path):
    """design_0_seq_0 ran two sequence_groups and cross_summary names group 0.
    The dispersion overview sees that group's three seeds; the grid sees all
    six, with the other group's three in grey."""
    paths = fc.cohort_tree(tmp_path)
    rep, all_seeds, _c, _t, _cb = wsp.load_cohort_seeds(paths["runs_dir"],
                                                        paths["cross_csv"])
    assert {s["sg"] for s in rep["design_0_seq_0"]} == {"0"}
    assert len(all_seeds["design_0_seq_0"]) == 6
    assert {s["sg"] for s in all_seeds["design_0_seq_0"]} == {"0", "1"}


@pytest.mark.logic
def test_load_cohort_seeds_reads_the_cold_start_baseline_from_the_initial_row(tmp_path):
    """The cycle-0 initial row has no sequence_group and no seed index: it is
    the single pre-steering prediction, shared by every seed, and it is the
    trajectory plot's leftmost point."""
    paths = fc.cohort_tree(tmp_path)
    _rep, all_seeds, _c, _t, cold = wsp.load_cohort_seeds(paths["runs_dir"],
                                                           paths["cross_csv"])
    assert cold["design_0_seq_0"] == 30.1
    assert cold["design_0_seq_1"] == 29.4
    # The baseline row is not itself a seed.
    assert all(s["sg"] != "" for s in all_seeds["design_0_seq_0"])
    # design_2_seq_0's fixture has no initial row at all.
    assert "design_2_seq_0" not in cold


@pytest.mark.logic
def test_load_cohort_seeds_skips_the_controls(tmp_path):
    """The within-sequence figures are about per-seed dispersion of steered
    sequences. Control behaviour is the cohort script's controls diagnostic."""
    paths = fc.cohort_tree(tmp_path)
    rep, all_seeds, _c, _t, _cb = wsp.load_cohort_seeds(paths["runs_dir"],
                                                        paths["cross_csv"])
    assert not [n for n in rep if n.startswith("input_control")]
    assert not [n for n in all_seeds if n.startswith("input_control")]


@pytest.mark.logic
def test_load_cohort_seeds_takes_tier_and_composite_from_the_summary(tmp_path):
    paths = fc.cohort_tree(tmp_path)
    _rep, _all, comps, tiers, _cb = wsp.load_cohort_seeds(paths["runs_dir"],
                                                           paths["cross_csv"])
    assert tiers["design_0_seq_0"] == "A"
    assert tiers["design_0_seq_1"] == "B"
    assert comps["design_0_seq_0"] == pytest.approx(0.70 - 0.05 * 2.0)


@pytest.mark.logic
def test_load_cohort_seeds_falls_back_to_the_best_aggregate_composite(tmp_path):
    """design_2_seq_0 is in runs/ but not in the summary, so it has no tier and
    no rep_* medians. Its composite comes from its best aggregated row
    (Jaccard 0.50, RMSD 4.0 Å) so it can still be ranked in the grid."""
    paths = fc.cohort_tree(tmp_path)
    _rep, _all, comps, tiers, _cb = wsp.load_cohort_seeds(paths["runs_dir"],
                                                           paths["cross_csv"])
    assert tiers["design_2_seq_0"] == "none"
    assert comps["design_2_seq_0"] == pytest.approx(0.50 - 0.05 * 4.0)


@pytest.mark.logic
def test_load_cohort_seeds_leaves_a_nan_metric_absent(tmp_path):
    """A prediction that never landed writes nan. It must arrive as None so the
    scatter drops the point, rather than as a number the axes scale to."""
    cross = tmp_path / "cross_sequence_summary.csv"
    fc.write_cross_summary(cross, [fc.cross_row("design_0_seq_0",
                                                 rep_sequence_group="0")])
    runs = tmp_path / "runs"
    fc.write_sequence_dir(
        runs, "design_0_seq_0",
        agg_rows=[fc.agg_row_for("no_reversion", steered_ra=2.0, steered_tj=0.7)],
        seed_rows=[fc.seed_row(stage="steering", outcome="no_data")],
    )
    rep, _all, _c, _t, _cb = wsp.load_cohort_seeds(runs, cross)
    assert rep["design_0_seq_0"][0]["ra_final"] is None
    assert rep["design_0_seq_0"][0]["outcome"] == "no_data"


@pytest.mark.logic
def test_load_cohort_seeds_skips_a_sequence_with_no_per_seed_file(tmp_path):
    cross = tmp_path / "cross_sequence_summary.csv"
    fc.write_cross_summary(cross, [fc.cross_row("design_0_seq_0")])
    runs = tmp_path / "runs"
    fc.write_sequence_dir(runs, "design_0_seq_0", agg_rows=[
        fc.agg_row_for("no_reversion", steered_ra=2.0, steered_tj=0.7)])
    rep, all_seeds, _c, _t, _cb = wsp.load_cohort_seeds(runs, cross)
    assert rep == {}
    assert all_seeds == {}


@pytest.mark.logic
def test_load_cohort_seeds_defaults_an_unparseable_seed_index_to_zero(tmp_path):
    cross = tmp_path / "cross_sequence_summary.csv"
    fc.write_cross_summary(cross, [fc.cross_row("design_0_seq_0",
                                                 rep_sequence_group="0")])
    runs = tmp_path / "runs"
    fc.write_sequence_dir(
        runs, "design_0_seq_0",
        agg_rows=[fc.agg_row_for("no_reversion", steered_ra=2.0, steered_tj=0.7)],
        seed_rows=[fc.seed_row(stage="steering", outcome="pass",
                                seed_index="not_an_index")],
    )
    rep, _all, _c, _t, _cb = wsp.load_cohort_seeds(runs, cross)
    assert rep["design_0_seq_0"][0]["seed"] == 0


# ═══════════════════════════════════════════════════════════════════════════
# Plot 1a — per-seed dispersion overview
# ═══════════════════════════════════════════════════════════════════════════

@pytest.mark.logic
def test_plot_dispersion_overview_draws_one_point_per_usable_seed(tmp_path, monkeypatch):
    """Three seeds, one of which never produced an RMSD. Two points."""
    rep = {"design_0_seq_0": [
        seed_rec(seed=0, stage="steering", ra=2.0),
        seed_rec(seed=1, stage="reversion", ra=3.0, reverted_ra=3.0),
        seed_rec(seed=2, stage="steering", ra=None),
    ]}
    out = tmp_path / "negsteer_per_seed_dispersion_overview.png"
    seen = capture_figure(monkeypatch)
    assert wsp.plot_dispersion_overview(rep, {"design_0_seq_0": "A"}, str(out)) is True
    assert seen["n_scatters"] == 2
    assert_drawn(out)


@pytest.mark.logic
def test_plot_dispersion_overview_still_writes_a_figure_with_no_seeds(tmp_path):
    """Nothing to scatter, but the axes and labels are drawn and the path
    Nextflow publishes is on disk."""
    out = tmp_path / "negsteer_per_seed_dispersion_overview.png"
    assert wsp.plot_dispersion_overview({}, {}, str(out)) is True
    assert_drawn(out)


@pytest.mark.logic
def test_plot_dispersion_overview_marks_an_unclassifiable_seed(tmp_path):
    """A seed whose stage could not be determined gets the fallback marker
    rather than being dropped."""
    rep = {"design_0_seq_0": [seed_rec(stage="no_data", outcome="no_data")]}
    out = tmp_path / "overview_no_data.png"
    assert wsp.plot_dispersion_overview(rep, {}, str(out)) is True
    assert_drawn(out)


# ═══════════════════════════════════════════════════════════════════════════
# Plot 1b — per-seed dispersion grid
# ═══════════════════════════════════════════════════════════════════════════

@pytest.mark.logic
def test_plot_dispersion_grid_draws_a_panel_per_sequence(tmp_path, monkeypatch):
    """Four sequences, so a 2 × 2 grid, each panel titled with the short name
    and ranked by composite."""
    rep, all_seeds, tiers, comps = {}, {}, {}, {}
    for i in range(4):
        name = f"design_{i}_seq_0"
        rep[name] = [seed_rec(seed=0, sg="0", ra=2.0 + i)]
        all_seeds[name] = rep[name] + [seed_rec(seed=0, sg="1", ra=9.0)]
        tiers[name] = ["A", "B", "C", "none"][i]
        comps[name] = 0.6 - 0.1 * i
    out = tmp_path / "negsteer_per_seed_dispersion_grid.png"
    seen = capture_figure(monkeypatch)
    assert wsp.plot_dispersion_grid(rep, all_seeds, tiers, comps, str(out)) is True
    assert [t for t in seen["titles"] if t] == ["d0_s0", "d1_s0", "d2_s0", "d3_s0"]
    assert_drawn(out)


@pytest.mark.logic
def test_plot_dispersion_grid_draws_a_single_panel(tmp_path):
    """One sequence means matplotlib hands back a bare Axes rather than an
    array of them."""
    rep = {"design_0_seq_0": [seed_rec(ra=2.0)]}
    out = tmp_path / "grid_single.png"
    assert wsp.plot_dispersion_grid(rep, rep, {}, {}, str(out)) is True
    assert_drawn(out)


@pytest.mark.logic
def test_plot_dispersion_grid_blanks_the_unused_panels(tmp_path):
    """Five sequences in a 2 × 3 grid leaves one panel with no sequence."""
    rep = {f"design_{i}_seq_0": [seed_rec(ra=2.0 + i)] for i in range(5)}
    out = tmp_path / "grid_five.png"
    assert wsp.plot_dispersion_grid(rep, rep, {}, {}, str(out)) is True
    assert_drawn(out)


@pytest.mark.logic
@pytest.mark.parametrize("rep, all_seeds", [
    ({}, {}),
    ({"design_0_seq_0": [seed_rec(ra=None, tj=None)]},
     {"design_0_seq_0": [seed_rec(ra=None, tj=None)]}),
], ids=["no-sequences", "no-rmsd"])
def test_plot_dispersion_grid_writes_a_placeholder_with_no_rmsd(tmp_path, rep, all_seeds):
    """With no RMSD anywhere there is nothing to scale the axes to. The
    placeholder keeps the published path populated and says why."""
    out = tmp_path / "negsteer_per_seed_dispersion_grid.png"
    assert wsp.plot_dispersion_grid(rep, all_seeds, {}, {}, str(out)) is False
    assert_drawn(out)


# ═══════════════════════════════════════════════════════════════════════════
# Plot 2 — per-design stage trajectories
# ═══════════════════════════════════════════════════════════════════════════

@pytest.mark.logic
def test_plot_stage_trajectories_draws_a_panel_per_design(tmp_path, monkeypatch):
    all_seeds = {
        "design_0_seq_0": [seed_rec(seed=i, stage="steering", ra=2.0)
                            for i in range(3)],
        "design_0_seq_1": [seed_rec(seed=i, stage="reversion", ra=30.0,
                                     steered_ra=2.5, reverted_ra=30.0)
                            for i in range(3)],
        "design_1_seq_0": [seed_rec(seed=0, stage="cold_start", ra=28.0)],
    }
    cold = {"design_0_seq_0": 30.1, "design_0_seq_1": 29.4,
            "design_1_seq_0": 28.0}
    out = tmp_path / "negsteer_per_design_stage_trajectories.png"
    seen = capture_figure(monkeypatch)
    assert wsp.plot_stage_trajectories(all_seeds, {}, cold, str(out)) is True
    assert [t for t in seen["titles"] if t] == ["Design 0", "Design 1"]
    assert_drawn(out)


@pytest.mark.logic
def test_plot_stage_trajectories_drops_the_per_sequence_legend_when_crowded(tmp_path):
    """A design with more than six sequences would spend the whole panel on its
    legend, so the panel keeps the lines and loses the key."""
    all_seeds = {f"design_0_seq_{i}": [seed_rec(seed=0, ra=2.0 + i)]
                 for i in range(7)}
    out = tmp_path / "trajectories_crowded.png"
    assert wsp.plot_stage_trajectories(all_seeds, {}, {}, str(out)) is True
    assert_drawn(out)


@pytest.mark.logic
def test_plot_stage_trajectories_plots_a_cold_start_only_sequence(tmp_path):
    """A cold-start seed has no steering point, so its trajectory is the single
    baseline point."""
    all_seeds = {"design_0_seq_0": [seed_rec(seed=0, stage="cold_start",
                                              ra=28.0, steered_ra=None)]}
    out = tmp_path / "trajectories_cold.png"
    assert wsp.plot_stage_trajectories(all_seeds, {},
                                       {"design_0_seq_0": 28.0}, str(out)) is True
    assert_drawn(out)


@pytest.mark.logic
def test_plot_stage_trajectories_writes_a_placeholder_with_no_rmsd(tmp_path):
    all_seeds = {"design_0_seq_0": [seed_rec(ra=None, steered_ra=None)]}
    out = tmp_path / "negsteer_per_design_stage_trajectories.png"
    assert wsp.plot_stage_trajectories(all_seeds, {}, {}, str(out)) is False
    assert_drawn(out)


# ═══════════════════════════════════════════════════════════════════════════
# Plot 3 — weighted vs true jaccard
# ═══════════════════════════════════════════════════════════════════════════

@pytest.mark.logic
def test_plot_weighted_vs_true_jaccard_draws(tmp_path, monkeypatch):
    rep = {"design_0_seq_0": [seed_rec(seed=i, tj=0.7, wj=0.6)
                               for i in range(3)]}
    out = tmp_path / "negsteer_weighted_vs_true_jaccard.png"
    seen = capture_figure(monkeypatch)
    assert wsp.plot_weighted_vs_true_jaccard(
        rep, {"design_0_seq_0": "A"}, str(out)) is True
    assert seen["n_scatters"] == 3
    assert "3 points" in seen["titles"][0]
    assert_drawn(out)


@pytest.mark.logic
def test_plot_weighted_vs_true_jaccard_says_so_when_there_is_no_weighted_set(
        tmp_path, monkeypatch):
    """weighted_jaccard arrived with this stage's metrics and an older run has
    none. The panel is drawn and labelled rather than left out."""
    rep = {"design_0_seq_0": [seed_rec(wj=None)]}
    out = tmp_path / "wj_absent.png"
    seen = capture_figure(monkeypatch)
    assert wsp.plot_weighted_vs_true_jaccard(rep, {}, str(out)) is True
    assert any("No weighted_jaccard data" in t for t in seen["texts"])
    assert_drawn(out)


# ═══════════════════════════════════════════════════════════════════════════
# Plot 4 — composite vs confidence metrics
# ═══════════════════════════════════════════════════════════════════════════

@pytest.mark.logic
def test_plot_composite_vs_confidence_draws_a_panel_per_metric(tmp_path, monkeypatch):
    rep = {f"design_{i}_seq_0": [seed_rec(seed=j) for j in range(3)]
           for i in range(3)}
    comps = {n: 0.6 for n in rep}
    tiers = {n: "A" for n in rep}
    out = tmp_path / "negsteer_composite_vs_confidence.png"
    seen = capture_figure(monkeypatch)
    assert wsp.plot_composite_vs_confidence(rep, comps, tiers, str(out)) is True
    titles = [t for t in seen["titles"] if t]
    assert "ipTM" in titles
    assert "ipae" in titles
    assert "average pLDDT (0-100)" in titles
    assert_drawn(out)


@pytest.mark.logic
def test_plot_composite_vs_confidence_drops_a_metric_no_sequence_carries(
        tmp_path, monkeypatch):
    """An older run has no ipSAE columns. The panel is dropped rather than
    drawn empty, so the grid stays readable."""
    rep = {"design_0_seq_0": [seed_rec(ipsae_min=None, ipsae_min_15=None,
                                        pae_mean=None)]}
    out = tmp_path / "confidence_fewer_panels.png"
    seen = capture_figure(monkeypatch)
    assert wsp.plot_composite_vs_confidence(
        rep, {"design_0_seq_0": 0.6}, {}, str(out)) is True
    titles = [t for t in seen["titles"] if t]
    assert "ipSAE_min (10 Å)" not in titles
    assert "pae_mean" not in titles
    assert "ipTM" in titles


@pytest.mark.logic
def test_plot_composite_vs_confidence_skips_a_sequence_with_no_composite(tmp_path):
    """A sequence with no computable composite has no y value, so it cannot be
    placed against a confidence metric."""
    rep = {"design_0_seq_0": [seed_rec()]}
    out = tmp_path / "confidence_no_composite.png"
    assert wsp.plot_composite_vs_confidence(
        rep, {"design_0_seq_0": float("-inf")}, {}, str(out)) is True
    assert_drawn(out)


# ═══════════════════════════════════════════════════════════════════════════
# CLI plumbing
# ═══════════════════════════════════════════════════════════════════════════

@pytest.mark.logic
def test_main_writes_every_published_figure(tmp_path, monkeypatch, capsys):
    """The five paths in the NEGSTEER_WITHIN_SEQUENCE_PLOTS publish list, all
    present and all with a figure in them."""
    paths = fc.cohort_tree(tmp_path)
    outdir = tmp_path / "plots"
    monkeypatch.setattr(sys, "argv", [
        "negsteer_within_sequence_plots.py",
        "--runs-dir", str(paths["runs_dir"]),
        "--cross-summary-csv", str(paths["cross_csv"]),
        "--outdir", str(outdir),
    ])
    assert wsp.main() == 0
    for name in (
        "negsteer_per_seed_dispersion_overview.png",
        "negsteer_per_seed_dispersion_grid.png",
        "negsteer_per_design_stage_trajectories.png",
        "negsteer_weighted_vs_true_jaccard.png",
        "negsteer_composite_vs_confidence.png",
    ):
        assert_drawn(outdir / name)
    assert "FAILED" not in capsys.readouterr().out


@pytest.mark.logic
def test_main_refuses_a_missing_runs_dir_or_summary(tmp_path, monkeypatch):
    """Both inputs are produced by an upstream process. A missing one is a
    broken run, not an empty cohort, so it fails loudly."""
    paths = fc.cohort_tree(tmp_path)
    monkeypatch.setattr(sys, "argv", [
        "negsteer_within_sequence_plots.py",
        "--runs-dir", str(tmp_path / "absent"),
        "--cross-summary-csv", str(paths["cross_csv"]),
        "--outdir", str(tmp_path / "plots"),
    ])
    with pytest.raises(SystemExit):
        wsp.main()

    monkeypatch.setattr(sys, "argv", [
        "negsteer_within_sequence_plots.py",
        "--runs-dir", str(paths["runs_dir"]),
        "--cross-summary-csv", str(tmp_path / "absent.csv"),
        "--outdir", str(tmp_path / "plots"),
    ])
    with pytest.raises(SystemExit):
        wsp.main()


@pytest.mark.logic
def test_main_returns_zero_when_no_sequence_has_seed_data(tmp_path, monkeypatch, capsys):
    cross = tmp_path / "cross_sequence_summary.csv"
    fc.write_cross_summary(cross, [fc.cross_row("design_0_seq_0")])
    runs = tmp_path / "runs"
    fc.write_sequence_dir(runs, "design_0_seq_0", agg_rows=[
        fc.agg_row_for("no_reversion", steered_ra=2.0, steered_tj=0.7)])
    monkeypatch.setattr(sys, "argv", [
        "negsteer_within_sequence_plots.py",
        "--runs-dir", str(runs), "--cross-summary-csv", str(cross),
        "--outdir", str(tmp_path / "plots"),
    ])
    assert wsp.main() == 0
    assert "Nothing to plot" in capsys.readouterr().out
