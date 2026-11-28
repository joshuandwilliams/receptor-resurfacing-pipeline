"""Drive bin/pose_solver_plots.py against constructed POSE_SOLVE output.

The five figures are the only readout of a pose-solver run, so a figure that
silently renders nothing is a silent loss of the diagnosis. Each test here
asserts a figure was written AND that it has more than one pixel colour, then
pins the decisions that are not visible from the file alone: which restart the
convergence curve calls the converged one, which side of a distance window
counts as satisfied, and which colour a histogram bin on a boundary gets.

Every input is built by tests/characterization/fixtures_plots.py from the keys
bin/pose_solver.py writes and bin/pose_interface_metrics.py emits.
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
import pose_solver_plots as psp  # noqa: E402
from fixtures_plots import (  # noqa: E402
    POSE_LOSS_WEIGHTS,
    assert_png_not_blank,
    figures_kept,
    pose_contig,
    pose_interface_metrics,
    pose_metrics_json,
    pose_pair,
    pose_restart_losses_csv,
    pose_results,
    pose_results_json,
)
from fixtures_runtree import complex_pdb  # noqa: E402
from matplotlib.colors import to_rgba  # noqa: E402

FULL_BREAKDOWN = {
    "validity": 0.0, "dist_upper": 2.5, "dist_lower": 0.0,
    "excl": 1.0, "interp": 4.0, "pair_sc_clash": 0.5, "total": 8.0,
}


def _facecolours(artists):
    return [tuple(a.get_facecolor()) for a in artists]


# ── Plot 1: loss breakdown ───────────────────────────────────────────

@pytest.mark.logic
def test_loss_breakdown_draws_one_bar_per_term_present(tmp_path):
    """Only the six known penalty terms become bars; `total` goes in the title."""
    out = tmp_path / "loss.png"
    with figures_kept() as figs:
        psp.plot_loss_breakdown(
            pose_results(loss_breakdown=FULL_BREAKDOWN,
                         loss_weights=POSE_LOSS_WEIGHTS), out)
        ax = figs()[-1].axes[0]
        assert [t.get_text() for t in ax.get_yticklabels()] == [
            "validity", "dist_upper", "dist_lower", "excl", "interp",
            "pair_sc_clash",
        ]
        assert "total = 8.00" in ax.get_title()
    assert_png_not_blank(out)


@pytest.mark.logic
def test_loss_breakdown_skips_absent_terms(tmp_path):
    """A breakdown carrying two terms draws two bars, in FILTER-order not dict
    order, so a reordered JSON cannot reorder the figure."""
    out = tmp_path / "loss.png"
    with figures_kept() as figs:
        psp.plot_loss_breakdown(
            pose_results(loss_breakdown={"interp": 3.0, "validity": 1.0}), out)
        labels = [t.get_text() for t in figs()[-1].axes[0].get_yticklabels()]
    assert labels == ["validity", "interp"]
    assert_png_not_blank(out)


@pytest.mark.logic
def test_loss_breakdown_with_no_breakdown_still_writes_a_file(tmp_path):
    """No-data path. The three committed solved_pose_results.json fixtures carry
    no loss_breakdown at all, so this is the common case, not an edge case.
    Nextflow publishes a fixed figure list, so the placeholder has to exist."""
    out = tmp_path / "loss.png"
    psp.plot_loss_breakdown(pose_results(), out)
    assert_png_not_blank(out)


@pytest.mark.logic
def test_loss_breakdown_all_zero_forces_a_unit_x_axis(tmp_path):
    """A perfect pose has every term at 0. Scaling to max(values) would give a
    degenerate (0, 0) axis, so the code pins x to (0, 1) instead."""
    out = tmp_path / "loss.png"
    zeros = {k: 0.0 for k in ("validity", "dist_upper", "excl")}
    with figures_kept() as figs:
        psp.plot_loss_breakdown(pose_results(loss_breakdown=zeros), out)
        assert figs()[-1].axes[0].get_xlim() == (0.0, 1.0)
    assert_png_not_blank(out)


@pytest.mark.logic
def test_loss_breakdown_annotates_unknown_weight_as_question_mark(tmp_path):
    """loss_weights is absent from every committed results JSON, so the bar
    annotation has to degrade to '?' rather than raise."""
    out = tmp_path / "loss.png"
    with figures_kept() as figs:
        psp.plot_loss_breakdown(pose_results(loss_breakdown={"interp": 3.0}), out)
        texts = [t.get_text() for t in figs()[-1].axes[0].texts]
    assert texts == ["3.00  (W_INTERP=?)"]


# ── Plot 2: pair-constraint distances ────────────────────────────────

@pytest.mark.logic
def test_pair_distances_marks_rows_in_input_order(tmp_path):
    """Row labels are chain-prefixed and follow pair_results order, not a sort,
    so a label can be traced back to the pair that produced it."""
    out = tmp_path / "pairs.png"
    results = pose_results(pair_results=[
        pose_pair(binder_res=73, target_res=76, achieved_d=4.008),
        pose_pair(binder_res=71, target_res=78, achieved_d=5.730),
    ])
    with figures_kept() as figs:
        psp.plot_pair_distances(results, "A", "B", out)
        labels = [t.get_text() for t in figs()[-1].axes[0].get_yticklabels()]
    assert labels == ["A73 – B76", "A71 – B78"]
    assert_png_not_blank(out)


@pytest.mark.logic
@pytest.mark.parametrize(
    ("achieved", "expect_satisfied"),
    [(4.0, True), (6.0, True), (3.999, False), (6.001, False), (5.0, True)],
)
def test_pair_distances_window_boundary_is_inclusive(tmp_path, achieved,
                                                     expect_satisfied):
    """The fallback satisfaction test is `min_d <= achieved <= max_d`, so BOTH
    comparisons are non-strict: a pair sitting exactly on min_d or on max_d is
    satisfied. Only used when `satisfied` is absent from the pair dict.
    """
    out = tmp_path / "pairs.png"
    results = pose_results(pair_results=[
        pose_pair(binder_res=10, target_res=20, achieved_d=achieved,
                  min_d=4.0, max_d=6.0),
    ])
    with figures_kept() as figs:
        psp.plot_pair_distances(results, "A", "B", out)
        marker = figs()[-1].axes[0].collections[0]
        colour = tuple(marker.get_facecolor()[0])
    want = to_rgba(psp.PALETTE["satisfied" if expect_satisfied else "unsatisfied"])
    assert colour == pytest.approx(want)


@pytest.mark.logic
def test_pair_distances_trusts_an_explicit_satisfied_flag(tmp_path):
    """An explicit `satisfied` from the solver wins over the window arithmetic,
    because the solver applies its own tolerance."""
    out = tmp_path / "pairs.png"
    results = pose_results(pair_results=[
        pose_pair(binder_res=10, target_res=20, achieved_d=5.0,
                  min_d=4.0, max_d=6.0, satisfied=False),
    ])
    with figures_kept() as figs:
        psp.plot_pair_distances(results, "A", "B", out)
        colour = tuple(figs()[-1].axes[0].collections[0].get_facecolor()[0])
    assert colour == pytest.approx(to_rgba(psp.PALETTE["unsatisfied"]))


@pytest.mark.logic
def test_pair_distances_with_no_pairs_still_writes_a_file(tmp_path):
    out = tmp_path / "pairs.png"
    psp.plot_pair_distances(pose_results(pair_results=[]), "A", "B", out)
    assert_png_not_blank(out)


# ── Plot 3: interface dashboard ──────────────────────────────────────

@pytest.mark.logic
def test_interface_dashboard_renders_four_panels(tmp_path):
    out = tmp_path / "dash.png"
    with figures_kept() as figs:
        psp.plot_interface_dashboard([pose_interface_metrics()], "A", "B", out)
        assert len(figs()[-1].axes) == 4
    assert_png_not_blank(out)


@pytest.mark.logic
def test_interface_dashboard_reads_only_the_first_entry(tmp_path):
    """pose_interface_metrics.py emits one entry per input PDB and the dashboard
    is documented as single-pose. The BSA bars must come from entry 0."""
    out = tmp_path / "dash.png"
    entries = [
        pose_interface_metrics(bsa_chain_a=100.0, bsa_chain_b=200.0),
        pose_interface_metrics(bsa_chain_a=9000.0, bsa_chain_b=9000.0),
    ]
    with figures_kept() as figs:
        psp.plot_interface_dashboard(entries, "A", "B", out)
        bsa_panel = figs()[-1].axes[0]
        heights = [p.get_height() for p in bsa_panel.patches]
    assert heights == [100.0, 200.0, 300.0]


@pytest.mark.logic
@pytest.mark.parametrize(
    ("bin_label", "expected"),
    [("1.9-2.0", "#7F0000"), ("2.0-2.5", "#DC2626"), ("2.4-2.5", "#DC2626"),
     ("2.5-3.0", "#F97316"), ("3.0-3.5", "#FACC15"), ("3.5-4.0", "#16A34A"),
     ("4.5-5.0", "#16A34A")],
)
def test_interface_dashboard_histogram_colour_boundaries_are_strict_less_than(
        tmp_path, bin_label, expected):
    """A bin's colour comes from its LOWER edge, compared with strict `<`:
    lo < 2.0 is the clash colour, lo < 2.5 the next, and so on. So the bin
    starting exactly at 2.0 is NOT the clash colour, and the bin starting
    exactly at 3.5 falls through to the comfortable green.
    """
    out = tmp_path / "dash.png"
    metrics = pose_interface_metrics(distance_histogram={bin_label: 7})
    with figures_kept() as figs:
        psp.plot_interface_dashboard([metrics], "A", "B", out)
        bar = figs()[-1].axes[1].patches[0]
    assert tuple(bar.get_facecolor()) == pytest.approx(to_rgba(expected))


@pytest.mark.logic
def test_interface_dashboard_without_a_histogram_says_so(tmp_path):
    """A degenerate metrics entry whose histogram is empty must annotate the
    panel rather than draw zero bars with no explanation."""
    out = tmp_path / "dash.png"
    with figures_kept() as figs:
        psp.plot_interface_dashboard(
            [pose_interface_metrics(distance_histogram={})], "A", "B", out)
        hist_panel = figs()[-1].axes[1]
        assert "No histogram" in [t.get_text() for t in hist_panel.texts]
        assert len(hist_panel.patches) == 0
    assert_png_not_blank(out)


@pytest.mark.logic
@pytest.mark.parametrize("gap_index", [None, float("nan")])
def test_interface_dashboard_missing_gap_index_says_so(tmp_path, gap_index):
    """null and NaN are both treated as missing. NaN reaches the plot because
    json.dumps writes bare NaN and json.loads accepts it."""
    out = tmp_path / "dash.png"
    with figures_kept() as figs:
        psp.plot_interface_dashboard(
            [pose_interface_metrics(gap_index=gap_index)], "A", "B", out)
        gap_panel = figs()[-1].axes[2]
        assert "No gap-index" in [t.get_text() for t in gap_panel.texts]
        assert len(gap_panel.patches) == 0
    assert_png_not_blank(out)


@pytest.mark.logic
def test_interface_dashboard_gap_index_bar_carries_the_value(tmp_path):
    out = tmp_path / "dash.png"
    with figures_kept() as figs:
        psp.plot_interface_dashboard(
            [pose_interface_metrics(gap_index=4.25)], "A", "B", out)
        gap_panel = figs()[-1].axes[2]
        assert gap_panel.patches[0].get_width() == pytest.approx(4.25)
        assert "4.25 Å" in [t.get_text() for t in gap_panel.texts]


@pytest.mark.logic
def test_interface_dashboard_with_no_entries_still_writes_a_file(tmp_path):
    out = tmp_path / "dash.png"
    psp.plot_interface_dashboard([], "A", "B", out)
    assert_png_not_blank(out)


# ── Plot 4: restart-loss convergence ─────────────────────────────────

@pytest.mark.logic
def test_restart_loss_curve_marks_the_first_restart_reaching_the_minimum(tmp_path):
    """The convergence annotation must name the EARLIEST restart that reached the
    global best, not the last one holding it. best_loss_so_far goes
    10, 10, 5, 5, 5, so restart 2 is the answer and restart 4 is not.
    """
    csv_path = pose_restart_losses_csv(tmp_path / "r.csv", [
        (0, 10.0, 10.0), (1, 12.0, 10.0), (2, 5.0, 5.0),
        (3, 9.0, 5.0), (4, 5.0, 5.0),
    ])
    out = tmp_path / "curve.png"
    with figures_kept() as figs:
        psp.plot_restart_loss_curve(csv_path, out)
        texts = " ".join(t.get_text() for t in figs()[-1].axes[0].texts)
    assert "converged at restart 2" in texts
    assert "loss = 5.000" in texts
    assert_png_not_blank(out)


@pytest.mark.logic
def test_restart_loss_curve_marks_the_row_position_not_the_restart_id(tmp_path):
    """Documented behaviour, not an endorsement: the marked restart is the ROW
    INDEX of the first best-so-far minimum, not the restart_idx value in that
    row. bin/pose_solver.py always writes restart_idx 0..n-1 contiguously, so
    the two coincide on real input and this only bites a hand-edited CSV.
    """
    csv_path = pose_restart_losses_csv(tmp_path / "r.csv", [
        (100, 10.0, 10.0), (200, 5.0, 5.0), (300, 6.0, 5.0),
    ])
    out = tmp_path / "curve.png"
    with figures_kept() as figs:
        psp.plot_restart_loss_curve(csv_path, out)
        texts = " ".join(t.get_text() for t in figs()[-1].axes[0].texts)
    assert "converged at restart 1" in texts


@pytest.mark.logic
def test_restart_loss_curve_single_row_uses_the_de_branch(tmp_path):
    """One row means the differential-evolution path ran, so there is no curve to
    draw and no convergence restart to name."""
    csv_path = pose_restart_losses_csv(tmp_path / "r.csv", [(0, 3.5, 3.5)])
    out = tmp_path / "curve.png"
    with figures_kept() as figs:
        psp.plot_restart_loss_curve(csv_path, out)
        ax = figs()[-1].axes[0]
        assert ax.get_xlim() == (-0.5, 0.5)
        assert len(ax.texts) == 0
    assert_png_not_blank(out)


@pytest.mark.logic
def test_restart_loss_curve_zero_loss_stays_linear(tmp_path):
    """A solved pose has loss 0, which has no logarithm. The y-scale switches to
    log only when every loss is strictly positive."""
    zero = pose_restart_losses_csv(tmp_path / "zero.csv",
                                   [(0, 1.0, 1.0), (1, 0.0, 0.0)])
    positive = pose_restart_losses_csv(tmp_path / "pos.csv",
                                       [(0, 1.0, 1.0), (1, 0.5, 0.5)])
    with figures_kept() as figs:
        psp.plot_restart_loss_curve(zero, tmp_path / "a.png")
        assert figs()[-1].axes[0].get_yscale() == "linear"
    with figures_kept() as figs:
        psp.plot_restart_loss_curve(positive, tmp_path / "b.png")
        assert figs()[-1].axes[0].get_yscale() == "log"


@pytest.mark.logic
def test_restart_loss_curve_header_only_csv_still_writes_a_file(tmp_path):
    """Zero rows. The CSV always has its header, so the reader succeeds and the
    emptiness only shows up after the loop."""
    csv_path = pose_restart_losses_csv(tmp_path / "r.csv", [])
    out = tmp_path / "curve.png"
    psp.plot_restart_loss_curve(csv_path, out)
    assert_png_not_blank(out)


# ── Plot 5: contig comparison ────────────────────────────────────────

def _three_contigs(design_region: str = "3-5") -> dict:
    return {
        "A) No joining": pose_contig(contig="A1-2/3-3/A6-10 B",
                                     design_region=design_region,
                                     n_residues=3, join_gap=0),
        "B) Join gap ≤ 1": pose_contig(contig="A1-2/3-3/A6-10 B",
                                       design_region=design_region,
                                       n_residues=3, join_gap=1),
        "C) Join gap ≤ 2": pose_contig(contig="A1-2/3-3/A6-10 B",
                                       design_region=design_region,
                                       n_residues=3, join_gap=2),
    }


@pytest.mark.logic
def test_contig_comparison_colours_exactly_the_design_region(tmp_path):
    """Ten binder residues, three contig options, design region 3-5 in each.
    So 3 of the 10 cells per row are orange and 7 are grey, 9 and 21 overall.
    An off-by-one in the range parse shows up here as a wrong count.
    """
    pdb = complex_pdb(tmp_path / "posed.pdb", n_receptor=10, n_effector=4)
    out = tmp_path / "contigs.png"
    with figures_kept() as figs:
        psp.plot_contig_comparison(
            pose_results(contigs=_three_contigs("3-5")), pdb, "A", out)
        cells = _facecolours(figs()[-1].axes[0].patches)
    design = to_rgba(psp.PALETTE["design"])
    fixed = to_rgba(psp.PALETTE["fixed"])
    assert sum(1 for c in cells if c == pytest.approx(design)) == 9
    assert sum(1 for c in cells if c == pytest.approx(fixed)) == 21
    assert_png_not_blank(out)


@pytest.mark.logic
def test_contig_comparison_parses_multi_range_and_bare_resnums(tmp_path):
    """design_region is comma-separated and each token is either a range or a
    single resnum, so "2-3, 7" means residues 2, 3 and 7: three cells."""
    pdb = complex_pdb(tmp_path / "posed.pdb", n_receptor=10, n_effector=4)
    out = tmp_path / "contigs.png"
    contigs = {"A) No joining": pose_contig(contig="A1 B",
                                            design_region="2-3, 7",
                                            n_residues=3, join_gap=0)}
    with figures_kept() as figs:
        psp.plot_contig_comparison(pose_results(contigs=contigs), pdb, "A", out)
        cells = _facecolours(figs()[-1].axes[0].patches)
    assert sum(1 for c in cells
               if c == pytest.approx(to_rgba(psp.PALETTE["design"]))) == 3


@pytest.mark.logic
def test_contig_comparison_treats_none_as_an_empty_region(tmp_path):
    """pose_solver.py writes the literal "(none)" for an empty design region, so
    every cell stays fixed-grey and nothing is parsed as a resnum."""
    pdb = complex_pdb(tmp_path / "posed.pdb", n_receptor=6, n_effector=3)
    out = tmp_path / "contigs.png"
    contigs = {"A) No joining": pose_contig(contig="A1-6 B",
                                            design_region="(none)",
                                            n_residues=0, join_gap=0)}
    with figures_kept() as figs:
        psp.plot_contig_comparison(pose_results(contigs=contigs), pdb, "A", out)
        cells = _facecolours(figs()[-1].axes[0].patches)
    assert sum(1 for c in cells
               if c == pytest.approx(to_rgba(psp.PALETTE["design"]))) == 0
    assert len(cells) == 6
    assert_png_not_blank(out)


@pytest.mark.logic
def test_contig_comparison_with_no_contigs_still_writes_a_file(tmp_path):
    pdb = complex_pdb(tmp_path / "posed.pdb", n_receptor=6, n_effector=3)
    out = tmp_path / "contigs.png"
    psp.plot_contig_comparison(pose_results(contigs={}), pdb, "A", out)
    assert_png_not_blank(out)


@pytest.mark.logic
def test_contig_comparison_with_an_absent_binder_chain_still_writes_a_file(tmp_path):
    """Wrong --receptor-chain is the realistic failure. The PDB parses, the chain
    filter matches nothing, and a placeholder naming the file has to appear."""
    pdb = complex_pdb(tmp_path / "posed.pdb", n_receptor=6, n_effector=3)
    out = tmp_path / "contigs.png"
    psp.plot_contig_comparison(pose_results(contigs=_three_contigs()),
                               pdb, "Z", out)
    assert_png_not_blank(out)


@pytest.mark.logic
def test_contig_comparison_single_residue_binder(tmp_path):
    """One row of one cell. The x-tick logic divides by the binder length, so a
    length of 1 is where an integer-division step of 0 would raise."""
    pdb = complex_pdb(tmp_path / "posed.pdb", n_receptor=1, n_effector=1)
    out = tmp_path / "contigs.png"
    contigs = {"A) No joining": pose_contig(contig="A1 B", design_region="1",
                                            n_residues=1, join_gap=0)}
    with figures_kept() as figs:
        psp.plot_contig_comparison(pose_results(contigs=contigs), pdb, "A", out)
        assert len(figs()[-1].axes[0].patches) == 1
    assert_png_not_blank(out)


@pytest.mark.logic
def test_contig_comparison_always_labels_the_last_binder_residue(tmp_path):
    """The x-axis aims for about twelve ticks by stepping len//12, which skips the
    final residue whenever the length is not a multiple of the step. A 24-residue
    binder steps by 2 and so stops at index 22, and the endpoint is appended.
    The last tick label must be the last binder resnum, 24, not 23.
    """
    pdb = complex_pdb(tmp_path / "posed.pdb", n_receptor=24, n_effector=4)
    out = tmp_path / "contigs.png"
    contigs = {"A) No joining": pose_contig(contig="A1-24 B",
                                            design_region="5-9",
                                            n_residues=5, join_gap=0)}
    with figures_kept() as figs:
        psp.plot_contig_comparison(pose_results(contigs=contigs), pdb, "A", out)
        labels = [t.get_text() for t in figs()[-1].axes[0].get_xticklabels()]
    assert labels[0] == "1"
    assert labels[-1] == "24"
    assert_png_not_blank(out)


# ── CLI ──────────────────────────────────────────────────────────────

@pytest.mark.logic
def test_main_writes_all_five_published_figures(tmp_path, monkeypatch):
    """modules/pose_solver.nf publishes a fixed list of five PNGs. A missing path
    fails the stage, so main() has to write every one of them from one run."""
    results = pose_results_json(
        tmp_path / "solved_pose_results.json",
        loss_breakdown=FULL_BREAKDOWN, loss_weights=POSE_LOSS_WEIGHTS,
        pair_results=[pose_pair(binder_res=73, target_res=76, achieved_d=4.0),
                      pose_pair(binder_res=71, target_res=78, achieved_d=7.2)],
        contigs=_three_contigs("3-5"),
    )
    metrics = pose_metrics_json(tmp_path / "interface_metrics.json",
                                [pose_interface_metrics()])
    restarts = pose_restart_losses_csv(tmp_path / "restarts.csv",
                                       [(0, 9.0, 9.0), (1, 2.0, 2.0)])
    pdb = complex_pdb(tmp_path / "solved_pose_posed.pdb",
                      n_receptor=10, n_effector=4)

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", [
        "pose_solver_plots.py",
        "--results-json", str(results),
        "--metrics-json", str(metrics),
        "--restart-losses", str(restarts),
        "--posed-pdb", str(pdb),
        "--receptor-chain", "a", "--effector-chain", "b",
    ])
    assert psp.main() == 0
    for name in ("pose_loss_breakdown.png", "pose_pair_distances.png",
                 "pose_interface_dashboard.png", "pose_restart_loss_curve.png",
                 "pose_contig_comparison.png"):
        assert_png_not_blank(tmp_path / name)


@pytest.mark.logic
def test_main_accepts_a_bare_dict_from_the_metrics_json(tmp_path, monkeypatch):
    """pose_interface_metrics.py emits a list, but main() normalises a single
    dict too, so a hand-run single-pose metrics file still plots."""
    results = pose_results_json(tmp_path / "res.json",
                                contigs=_three_contigs("2-3"))
    metrics = pose_metrics_json(tmp_path / "m.json",
                                pose_interface_metrics(bsa_chain_a=500.0))
    restarts = pose_restart_losses_csv(tmp_path / "r.csv", [(0, 1.0, 1.0)])
    pdb = complex_pdb(tmp_path / "p.pdb", n_receptor=8, n_effector=3)
    assert json.loads(metrics.read_text()).__class__ is dict

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", [
        "pose_solver_plots.py", "--results-json", str(results),
        "--metrics-json", str(metrics), "--restart-losses", str(restarts),
        "--posed-pdb", str(pdb),
    ])
    assert psp.main() == 0
    assert_png_not_blank(tmp_path / "pose_interface_dashboard.png")


@pytest.mark.logic
def test_main_on_an_empty_run_still_writes_all_five(tmp_path, monkeypatch):
    """Nothing solved: no pairs, no contigs, no loss breakdown, no restarts, and
    an empty metrics list. All five published paths must still exist."""
    results = pose_results_json(tmp_path / "res.json")
    metrics = pose_metrics_json(tmp_path / "m.json", [])
    restarts = pose_restart_losses_csv(tmp_path / "r.csv", [])
    pdb = complex_pdb(tmp_path / "p.pdb", n_receptor=4, n_effector=2)

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", [
        "pose_solver_plots.py", "--results-json", str(results),
        "--metrics-json", str(metrics), "--restart-losses", str(restarts),
        "--posed-pdb", str(pdb),
    ])
    assert psp.main() == 0
    for name in ("pose_loss_breakdown.png", "pose_pair_distances.png",
                 "pose_interface_dashboard.png", "pose_restart_loss_curve.png",
                 "pose_contig_comparison.png"):
        assert_png_not_blank(tmp_path / name)
