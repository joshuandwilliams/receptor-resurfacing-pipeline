"""Drive bin/mpnn_plots.py against constructed scored_metadata.csv rows.

The four figures are the only readout of a ProteinMPNN design round, and three
of the four rank something before they draw it: the score plot ranks designs by
median design-region score, the composition plot ranks amino acids by median
design frequency, and both composition and physicochem pick one native
reference string out of however many the rows carry. Those rankings exist
nowhere but inside the figure, so each is asserted on an input whose right
answer is arithmetic.

Every plot function writes to a bare relative filename, so each test chdirs into
tmp_path first. Row shapes come from
tests/proteinmpnn/example_output_files/sequences/scored_metadata.csv.
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

import mpnn_plots as mp  # noqa: E402
from fixtures_plots import (  # noqa: E402
    assert_png_not_blank,
    figures_kept,
    mpnn_cluster_csv,
    mpnn_metadata_csv,
    mpnn_row,
)


def _xticklabels(ax):
    return [t.get_text() for t in ax.get_xticklabels()]


def _legend_texts(ax):
    legend = ax.get_legend()
    return [] if legend is None else [t.get_text() for t in legend.get_texts()]


# ── Plot 1: score distribution ───────────────────────────────────────

@pytest.mark.logic
def test_score_distribution_ranks_designs_by_median_design_region_score(
        tmp_path, monkeypatch):
    """Lower MPNN score is better, so "best first" means ascending median. Design
    1 has median 1.0, design 2 has 3.0 and design 0 has 5.0, so the x-axis must
    read d1, d2, d0 — the input order 0, 1, 2 is the wrong answer.
    """
    monkeypatch.chdir(tmp_path)
    rows = [
        mpnn_row(design="0", mpnn_score=9.0, design_region_score=5.0),
        mpnn_row(design="1", mpnn_score=9.0, design_region_score=1.0),
        mpnn_row(design="2", mpnn_score=9.0, design_region_score=3.0),
    ]
    with figures_kept() as figs:
        mp.plot_score_distribution(rows)
        ax_bot = figs()[-1].axes[1]
        assert _xticklabels(ax_bot) == ["d1", "d2", "d0"]
        assert "design-region score" in ax_bot.get_xlabel()
    assert_png_not_blank(tmp_path / "mpnn_score_distribution.png")


@pytest.mark.logic
def test_score_distribution_takes_the_median_not_the_best_replicate(
        tmp_path, monkeypatch):
    """Ranking is by median across a design's sequences, so one outstanding
    sequence does not promote its design. Design 0's scores are 1, 9, 9 (median
    9) and design 1's are 4, 5, 6 (median 5), so design 1 ranks first even though
    design 0 holds the single best sequence.
    """
    monkeypatch.chdir(tmp_path)
    rows = (
        [mpnn_row(design="0", seq=str(i), design_region_score=v)
         for i, v in enumerate((1.0, 9.0, 9.0))]
        + [mpnn_row(design="1", seq=str(i), design_region_score=v)
           for i, v in enumerate((4.0, 5.0, 6.0))]
    )
    with figures_kept() as figs:
        mp.plot_score_distribution(rows)
        assert _xticklabels(figs()[-1].axes[1]) == ["d1", "d0"]


@pytest.mark.logic
def test_score_distribution_appends_designs_with_no_region_score_last(
        tmp_path, monkeypatch):
    """A design carrying a global score but no design-region score cannot be
    ranked against the others, so it is appended after the ranked block rather
    than dropped. d1 has the better global score and still comes last."""
    monkeypatch.chdir(tmp_path)
    rows = [
        mpnn_row(design="0", mpnn_score=9.0, design_region_score=4.0),
        mpnn_row(design="1", mpnn_score=0.1, design_region_score="nan"),
    ]
    with figures_kept() as figs:
        mp.plot_score_distribution(rows)
        assert _xticklabels(figs()[-1].axes[1]) == ["d0", "d1"]


@pytest.mark.logic
def test_score_distribution_falls_back_to_global_ranking(tmp_path, monkeypatch):
    """A whole column of nan — every design_region_score missing — switches the
    ranking to the global score and replaces the top panel with a note."""
    monkeypatch.chdir(tmp_path)
    rows = [
        mpnn_row(design="0", mpnn_score=5.0, design_region_score="nan"),
        mpnn_row(design="1", mpnn_score=2.0, design_region_score=""),
    ]
    with figures_kept() as figs:
        mp.plot_score_distribution(rows)
        fig = figs()[-1]
        assert _xticklabels(fig.axes[1]) == ["d1", "d0"]
        assert "ranked by median global score" in fig.axes[1].get_xlabel()
        assert any("not available" in t.get_text() for t in fig.axes[0].texts)
    assert_png_not_blank(tmp_path / "mpnn_score_distribution.png")


@pytest.mark.logic
def test_score_distribution_with_no_scores_still_writes_a_file(
        tmp_path, monkeypatch):
    """No-data path. Non-numeric scores leave both dictionaries empty, and the
    MPNN_PLOTS stage publishes a fixed figure list, so a placeholder PNG has to
    appear rather than nothing."""
    monkeypatch.chdir(tmp_path)
    mp.plot_score_distribution([
        mpnn_row(design="0", mpnn_score="", design_region_score=""),
    ])
    assert_png_not_blank(tmp_path / "mpnn_score_distribution.png")


@pytest.mark.logic
def test_score_distribution_single_row(tmp_path, monkeypatch):
    """One design, one sequence. A boxplot of one value is degenerate but must
    still render."""
    monkeypatch.chdir(tmp_path)
    with figures_kept() as figs:
        mp.plot_score_distribution([
            mpnn_row(design="7", mpnn_score=1.2, design_region_score=1.5)])
        assert _xticklabels(figs()[-1].axes[1]) == ["d7"]
    assert_png_not_blank(tmp_path / "mpnn_score_distribution.png")


# ── Plot 2: sequence diversity ───────────────────────────────────────

@pytest.mark.logic
def test_sequence_diversity_converts_the_threshold_to_a_percentage(
        tmp_path, monkeypatch):
    """mpnn_cluster_counts.csv stores identity as a fraction and the x-axis is a
    percentage, so 0.30-1.00 in the file is 30-100 on the axis, padded by 2."""
    monkeypatch.chdir(tmp_path)
    csv_path = mpnn_cluster_csv(tmp_path / "clusters.csv",
                               [(0.30, 1), (0.60, 2), (1.00, 4)])
    with figures_kept() as figs:
        mp.plot_sequence_diversity(str(csv_path))
        ax = figs()[-1].axes[0]
        assert ax.get_xlim() == (28.0, 102.0)
        xs, ys = ax.lines[0].get_data()
        assert list(xs) == [30.0, 60.0, 100.0]
        assert list(ys) == [1, 2, 4]
    assert_png_not_blank(tmp_path / "mpnn_sequence_diversity.png")


@pytest.mark.logic
def test_sequence_diversity_single_threshold(tmp_path, monkeypatch):
    """One row means min == max, which would give a zero-width axis without the
    padding."""
    monkeypatch.chdir(tmp_path)
    csv_path = mpnn_cluster_csv(tmp_path / "clusters.csv", [(0.50, 3)])
    with figures_kept() as figs:
        mp.plot_sequence_diversity(str(csv_path))
        assert figs()[-1].axes[0].get_xlim() == (48.0, 52.0)
    assert_png_not_blank(tmp_path / "mpnn_sequence_diversity.png")


@pytest.mark.logic
@pytest.mark.parametrize("arg", ["", "does_not_exist.csv"])
def test_sequence_diversity_without_clustering_still_writes_a_file(
        tmp_path, monkeypatch, arg):
    """--cluster-csv defaults to the empty string and the file is optional, so
    both an unset flag and a missing path must produce the placeholder."""
    monkeypatch.chdir(tmp_path)
    mp.plot_sequence_diversity(arg)
    assert_png_not_blank(tmp_path / "mpnn_sequence_diversity.png")


@pytest.mark.logic
def test_sequence_diversity_header_only_csv_still_writes_a_file(
        tmp_path, monkeypatch):
    """Zero rows. The file exists and parses, so the emptiness is only visible
    after the read loop."""
    monkeypatch.chdir(tmp_path)
    csv_path = mpnn_cluster_csv(tmp_path / "clusters.csv", [])
    mp.plot_sequence_diversity(str(csv_path))
    assert_png_not_blank(tmp_path / "mpnn_sequence_diversity.png")


# ── Plot 3: AA composition ───────────────────────────────────────────

@pytest.mark.logic
def test_aa_composition_pct_counts_only_amino_acids():
    """The design-region string is pipe-joined across regions, so the separator
    and anything else non-canonical must not reach the denominator. "AA|CC" is
    four residues, 50% A and 50% C, not five."""
    pct = mp._aa_composition_pct("AA|CC")
    assert pct[mp.AA_ORDER.index("A")] == pytest.approx(50.0)
    assert pct[mp.AA_ORDER.index("C")] == pytest.approx(50.0)
    assert pct.sum() == pytest.approx(100.0)


@pytest.mark.logic
def test_aa_composition_pct_of_an_empty_string_is_all_zero():
    assert mp._aa_composition_pct("").sum() == 0.0


@pytest.mark.logic
def test_aa_composition_ranks_amino_acids_by_median_design_frequency(
        tmp_path, monkeypatch):
    """Ten designed residues, AAAACCCGGW: A 40%, C 30%, G 20%, W 10%. The x-axis
    is ranked by descending median design frequency, so the first four labels
    must be A, C, G, W. The remaining sixteen are all tied at 0% and their
    relative order is not defined by the code, so they are not asserted.
    """
    monkeypatch.chdir(tmp_path)
    rows = [mpnn_row(design="0", designed_residues="AAAACCCGGW",
                     native_residues="AAAACCCGGW")]
    with figures_kept() as figs:
        mp.plot_aa_composition(rows)
        ax_box = figs()[-1].axes[2]
        assert _xticklabels(ax_box)[:4] == ["A", "C", "G", "W"]
    assert_png_not_blank(tmp_path / "mpnn_aa_composition.png")


@pytest.mark.logic
def test_aa_composition_native_diamonds_sit_at_the_native_percentage(
        tmp_path, monkeypatch):
    """The native overlay is computed from native_residues, not from the designs.
    An all-D native against an all-A design puts the red diamond at 100% above
    D's column and at 0% above A's."""
    monkeypatch.chdir(tmp_path)
    rows = [mpnn_row(design="0", designed_residues="AAAA",
                     native_residues="DDDD")]
    with figures_kept() as figs:
        mp.plot_aa_composition(rows)
        ax_box = figs()[-1].axes[2]
        labels = _xticklabels(ax_box)
        offsets = {labels[int(x)]: y
                   for x, y in ax_box.collections[-1].get_offsets()}
        assert offsets["D"] == pytest.approx(100.0)
        assert offsets["A"] == pytest.approx(0.0)
        assert "Native" in _legend_texts(ax_box)


@pytest.mark.logic
def test_aa_composition_picks_the_longest_native_string(tmp_path, monkeypatch):
    """Every row of one run carries the same native string, because the gap
    positions come from the contig spec. Two distinct strings mean concatenated
    runs, and the documented resolution is to warn and use the longest. Here the
    longest is the all-E one, so the diamond sits over E and not over D.
    """
    monkeypatch.chdir(tmp_path)
    rows = [
        mpnn_row(design="0", designed_residues="AAAA", native_residues="DD"),
        mpnn_row(design="1", designed_residues="AAAA", native_residues="EEEE"),
    ]
    with figures_kept() as figs:
        mp.plot_aa_composition(rows)
        ax_box = figs()[-1].axes[2]
        labels = _xticklabels(ax_box)
        offsets = {labels[int(x)]: y
                   for x, y in ax_box.collections[-1].get_offsets()}
    assert offsets["E"] == pytest.approx(100.0)
    assert offsets["D"] == pytest.approx(0.0)


@pytest.mark.logic
def test_aa_composition_warns_when_native_strings_disagree(
        tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    mp.plot_aa_composition([
        mpnn_row(design="0", designed_residues="AAAA", native_residues="DD"),
        mpnn_row(design="1", designed_residues="AAAA", native_residues="EEEE"),
    ])
    assert "found 2 distinct native_residues" in capsys.readouterr().err


@pytest.mark.logic
def test_aa_composition_without_a_native_string_omits_the_overlay(
        tmp_path, monkeypatch):
    """native_residues is blank when the input PDB had no gap residues to read.
    The boxplot still draws; only the diamonds and their legend go away."""
    monkeypatch.chdir(tmp_path)
    with figures_kept() as figs:
        mp.plot_aa_composition([
            mpnn_row(design="0", designed_residues="AAAACCCC",
                     native_residues="")])
        ax_box = figs()[-1].axes[2]
        assert _legend_texts(ax_box) == []
    assert_png_not_blank(tmp_path / "mpnn_aa_composition.png")


@pytest.mark.logic
def test_aa_composition_with_no_designed_residues_still_writes_a_file(
        tmp_path, monkeypatch):
    """No-data path: every designed_residues blank, so there is nothing to count
    and the placeholder must be written."""
    monkeypatch.chdir(tmp_path)
    mp.plot_aa_composition([mpnn_row(design="0", designed_residues="",
                                     native_residues="")])
    assert_png_not_blank(tmp_path / "mpnn_aa_composition.png")


# ── Plot 4: physicochemistry ─────────────────────────────────────────

@pytest.mark.logic
def test_mean_hydrophobicity_averages_the_kyte_doolittle_index():
    """A is +1.8 and I is +4.5, so AI averages to 3.15. Non-index characters are
    skipped, which is how the pipe separator survives."""
    assert mp._mean_hydrophobicity("AI") == pytest.approx(3.15)
    assert mp._mean_hydrophobicity("A|I") == pytest.approx(3.15)
    assert mp._mean_hydrophobicity("") != mp._mean_hydrophobicity("")  # nan
    assert mp._mean_hydrophobicity("XX") != mp._mean_hydrophobicity("XX")


@pytest.mark.logic
def test_net_charge_counts_kr_as_plus_one_and_de_as_minus_one():
    """Histidine is about +0.1 at neutral pH and is deliberately excluded, so
    KRDEH is net zero rather than slightly positive."""
    assert mp._net_charge("KR") == pytest.approx(2.0)
    assert mp._net_charge("DE") == pytest.approx(-2.0)
    assert mp._net_charge("KRDEH") == pytest.approx(0.0)
    assert mp._net_charge("") == 0.0


@pytest.mark.logic
def test_physicochem_native_reference_lines_carry_the_native_values(
        tmp_path, monkeypatch):
    """Both panels draw the native value as a dashed reference whose legend text
    is the number itself. Native IIII is KD 4.50; native KKKK is +1.000 charge
    per residue. Charge is normalised by length so a longer region is not
    penalised, which is why the label is per-residue rather than a total.
    """
    monkeypatch.chdir(tmp_path)
    rows = [mpnn_row(design="0", designed_residues="AAAA",
                     native_residues="IIII")]
    with figures_kept() as figs:
        mp.plot_physicochem(rows)
        fig = figs()[-1]
        assert _legend_texts(fig.axes[0]) == ["Native = 4.50"]
    rows = [mpnn_row(design="0", designed_residues="KKDD",
                     native_residues="KKKK")]
    with figures_kept() as figs:
        mp.plot_physicochem(rows)
        fig = figs()[-1]
        assert _legend_texts(fig.axes[1]) == ["Native = +1.000"]
    assert_png_not_blank(tmp_path / "mpnn_physicochem.png")


@pytest.mark.logic
def test_physicochem_charge_is_per_residue_across_both_regions(
        tmp_path, monkeypatch):
    """The separator is stripped before the division, so "KK|DDDD" is six
    residues with net -2, giving -1/3 per residue. Dividing by the string length
    instead of the residue count would give -0.286.
    """
    monkeypatch.chdir(tmp_path)
    with figures_kept() as figs:
        mp.plot_physicochem([mpnn_row(design="0", designed_residues="KK|DDDD",
                                      native_residues="AAAAAA")])
        charge_panel = figs()[-1].axes[1]
        medians = [ln.get_ydata()[0] for ln in charge_panel.lines
                   if len(ln.get_ydata()) and ln.get_linestyle() == "-"]
    assert any(v == pytest.approx(-1 / 3) for v in medians)


@pytest.mark.logic
def test_physicochem_orders_designs_by_name(tmp_path, monkeypatch):
    """Unlike the score plot, physicochem does not rank: designs appear in sorted
    key order so the panel can be read against the AA-composition heatmap."""
    monkeypatch.chdir(tmp_path)
    rows = [mpnn_row(design=d, designed_residues="AAAA", native_residues="AAAA")
            for d in ("2", "0", "1")]
    with figures_kept() as figs:
        mp.plot_physicochem(rows)
        assert _xticklabels(figs()[-1].axes[1]) == ["d0", "d1", "d2"]


@pytest.mark.logic
def test_physicochem_with_no_designed_residues_still_writes_a_file(
        tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    mp.plot_physicochem([mpnn_row(design="0", designed_residues="",
                                  native_residues="")])
    assert_png_not_blank(tmp_path / "mpnn_physicochem.png")


@pytest.mark.logic
def test_physicochem_warns_when_native_strings_disagree(
        tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    mp.plot_physicochem([
        mpnn_row(design="0", designed_residues="AAAA", native_residues="DD"),
        mpnn_row(design="1", designed_residues="AAAA", native_residues="EEEE"),
    ])
    assert "found 2 distinct native_residues" in capsys.readouterr().err


# ── Fallbacks and CLI ────────────────────────────────────────────────

@pytest.mark.logic
def test_save_fallback_plots_covers_every_published_figure(tmp_path, monkeypatch):
    """modules/proteinmpnn.nf publishes these four names, so the fallback has to
    write all four or the stage fails with nothing in the log."""
    monkeypatch.chdir(tmp_path)
    mp.save_fallback_plots("No metadata")
    assert mp.ALL_PLOT_FILES == [
        "mpnn_score_distribution.png", "mpnn_sequence_diversity.png",
        "mpnn_aa_composition.png", "mpnn_physicochem.png",
    ]
    for name in mp.ALL_PLOT_FILES:
        assert_png_not_blank(tmp_path / name)


@pytest.mark.logic
def test_main_writes_all_four_figures(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    metadata = mpnn_metadata_csv(tmp_path / "scored_metadata.csv", [
        mpnn_row(design="0", seq="0", mpnn_score=1.35, design_region_score=1.51,
                 designed_residues="AVSVTEKGGWSDIT|GDIVVP",
                 native_residues="SVALVGDLRDKIE|SQANKD"),
        mpnn_row(design="0", seq="1", mpnn_score=1.32, design_region_score=1.48,
                 designed_residues="AVYAEKEDGWWRIT|GEITVP",
                 native_residues="SVALVGDLRDKIE|SQANKD"),
        mpnn_row(design="28", seq="0", mpnn_score=1.40, design_region_score=1.20,
                 designed_residues="KKDDEEWWFFGGHH|AAAAAA",
                 native_residues="SVALVGDLRDKIE|SQANKD"),
    ])
    clusters = mpnn_cluster_csv(tmp_path / "mpnn_cluster_counts.csv",
                               [(t / 10, n) for t, n in
                                ((3, 1), (5, 2), (8, 3), (10, 3))])
    monkeypatch.setattr(sys, "argv", [
        "mpnn_plots.py", "--metadata", str(metadata),
        "--cluster-csv", str(clusters),
    ])
    mp.main()
    for name in mp.ALL_PLOT_FILES:
        assert_png_not_blank(tmp_path / name)


@pytest.mark.logic
def test_main_with_a_missing_metadata_file_writes_fallbacks(tmp_path, monkeypatch):
    """No-data path through the CLI. main() exits 0 so the Nextflow stage
    succeeds, which is only safe because all four files were written."""
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", [
        "mpnn_plots.py", "--metadata", str(tmp_path / "absent.csv")])
    with pytest.raises(SystemExit) as exc:
        mp.main()
    assert exc.value.code == 0
    for name in mp.ALL_PLOT_FILES:
        assert_png_not_blank(tmp_path / name)


@pytest.mark.logic
def test_main_with_an_empty_metadata_file_writes_fallbacks(tmp_path, monkeypatch):
    """A zero-byte CSV is what an upstream crash leaves behind."""
    monkeypatch.chdir(tmp_path)
    empty = tmp_path / "scored_metadata.csv"
    empty.write_text("")
    monkeypatch.setattr(sys, "argv", ["mpnn_plots.py", "--metadata", str(empty)])
    with pytest.raises(SystemExit) as exc:
        mp.main()
    assert exc.value.code == 0
    for name in mp.ALL_PLOT_FILES:
        assert_png_not_blank(tmp_path / name)


@pytest.mark.logic
def test_main_with_a_header_only_metadata_file_writes_fallbacks(
        tmp_path, monkeypatch):
    """Header but zero rows: the design round produced no sequences."""
    monkeypatch.chdir(tmp_path)
    metadata = mpnn_metadata_csv(tmp_path / "scored_metadata.csv", [])
    monkeypatch.setattr(sys, "argv",
                        ["mpnn_plots.py", "--metadata", str(metadata)])
    with pytest.raises(SystemExit) as exc:
        mp.main()
    assert exc.value.code == 0
    for name in mp.ALL_PLOT_FILES:
        assert_png_not_blank(tmp_path / name)


@pytest.mark.logic
def test_main_without_matplotlib_writes_placeholder_files(tmp_path, monkeypatch):
    """When the import of matplotlib failed, main() still creates each published
    path so the stage does not fail on a missing output. The files carry a text
    marker rather than PNG bytes, which Nextflow accepts because it matches on
    the filename, but anything that opens them as images will not.
    """
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(mp, "HAS_MPL", False)
    monkeypatch.setattr(sys, "argv",
                        ["mpnn_plots.py", "--metadata", str(tmp_path / "x.csv")])
    with pytest.raises(SystemExit) as exc:
        mp.main()
    assert exc.value.code == 0
    for name in mp.ALL_PLOT_FILES:
        assert (tmp_path / name).read_text() == "matplotlib not available"
