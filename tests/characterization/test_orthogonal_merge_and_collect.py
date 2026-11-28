"""Cover the orthogonal-metric merge and the Rosetta score collector.

`_apply_filters` decides `passes_orthogonal_filters`, which is the column the
cohort tables and the gBlock selection read. Which metrics gate and which only
annotate is a research decision written into the code, so the tests state it
rather than restating the implementation:

    Sc and BSA       gate
    Rosetta ddG      carried through, no threshold applied
    AF3 ra_eff       informational only, gating until 2026-04-28
    interface_pLDDT  never checked, and --filter-plddt-min has no effect

That last one is what `notes/things_to_fix.md` records under "Interface pLDDT
never reaches the output table, and its threshold is dead". The test here
makes it a property the code must keep rather than a note someone has to find.
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

import merge_orthogonal_metrics as mom  # noqa: E402
import rosetta_filter_collect as rfc  # noqa: E402

PASSING = {
    "af3_nomsa_best_ra_eff": "2.0",
    "sc": "0.700",
    "bsa": "900",
    "rosetta_ddg": "-50.0",
    "interface_plddt": "0.10",
}
LIMITS = dict(af3_ra_max=5.0, sc_min=0.300, bsa_min=400.0, plddt_min=0.75)


def _filtered(**overrides):
    row = dict(PASSING)
    row.update(overrides)
    mom._apply_filters(row, **LIMITS)
    return row


# ── Float coercion ───────────────────────────────────────────────────

@pytest.mark.logic
@pytest.mark.parametrize("raw", ["", "   ", "NA", "none", None])
def test_an_unusable_value_becomes_none_not_zero(raw):
    """0.0 would read as a real measurement of a terrible interface and gate
    the design out for the wrong reason."""
    assert mom._as_float(raw) is None


@pytest.mark.logic
def test_a_real_number_is_read():
    assert mom._as_float("0.665") == pytest.approx(0.665)
    assert mom._as_float("-59.739") == pytest.approx(-59.739)


# ── What gates ───────────────────────────────────────────────────────

@pytest.mark.logic
def test_a_fully_passing_row_passes_with_no_flags():
    row = _filtered()
    assert row["passes_orthogonal_filters"] == "1"
    assert row["orthogonal_flags"] == ""


@pytest.mark.logic
def test_shape_complementarity_below_the_minimum_gates():
    row = _filtered(sc="0.100")
    assert row["passes_orthogonal_filters"] == "0"
    assert "sc_too_low" in row["orthogonal_flags"]


@pytest.mark.logic
def test_buried_surface_area_below_the_minimum_gates():
    row = _filtered(bsa="100")
    assert row["passes_orthogonal_filters"] == "0"
    assert "bsa_too_low" in row["orthogonal_flags"]


@pytest.mark.logic
@pytest.mark.parametrize("col,flag", [("sc", "sc_missing"), ("bsa", "bsa_missing")])
def test_a_missing_gating_metric_gates(col, flag):
    """Absent is not the same as passing. A design whose Rosetta run died
    must not be carried forward as though it had passed."""
    row = _filtered(**{col: ""})
    assert row["passes_orthogonal_filters"] == "0"
    assert flag in row["orthogonal_flags"]


@pytest.mark.logic
def test_the_thresholds_are_not_strict_so_exact_equality_passes():
    """The comparison is `value < minimum`, so a value sitting exactly on the
    threshold is kept. That makes the published threshold an attainable
    value."""
    assert _filtered(sc="0.300")["passes_orthogonal_filters"] == "1"
    assert _filtered(bsa="400")["passes_orthogonal_filters"] == "1"


# ── What does not gate ───────────────────────────────────────────────

@pytest.mark.logic
def test_af3_disagreement_is_flagged_but_does_not_gate():
    """It gated until 2026-04-28. Negative steering optimises against Boltz,
    so AF3 disagreement is a sanity-check warning rather than a kill
    criterion."""
    row = _filtered(af3_nomsa_best_ra_eff="25.0")
    assert "af3_nomsa_ra_eff_too_high" in row["orthogonal_flags"]
    assert row["passes_orthogonal_filters"] == "1"


@pytest.mark.logic
def test_a_missing_af3_value_is_flagged_but_does_not_gate():
    row = _filtered(af3_nomsa_best_ra_eff="")
    assert "af3_nomsa_missing" in row["orthogonal_flags"]
    assert row["passes_orthogonal_filters"] == "1"


@pytest.mark.logic
def test_rosetta_ddg_is_carried_but_never_thresholded():
    """A wildly favourable and a wildly unfavourable ddG both pass, because
    no threshold is applied to it."""
    for ddg in ("-500.0", "500.0", ""):
        row = _filtered(rosetta_ddg=ddg)
        assert row["passes_orthogonal_filters"] == "1"
        assert "ddg" not in row["orthogonal_flags"].lower()


@pytest.mark.logic
def test_interface_plddt_is_never_checked_and_its_threshold_is_dead():
    """notes/things_to_fix.md records this. The flag is required by the CLI
    so the Nextflow caller need not change, and it has no effect.

    Both a pLDDT far below the minimum and a changed minimum leave the row
    passing with no pLDDT flag.
    """
    row = _filtered(interface_plddt="0.01")
    assert row["passes_orthogonal_filters"] == "1"
    assert "plddt" not in row["orthogonal_flags"].lower()

    strict = dict(PASSING)
    mom._apply_filters(strict, **{**LIMITS, "plddt_min": 0.99})
    assert strict["passes_orthogonal_filters"] == "1"
    assert "plddt" not in strict["orthogonal_flags"].lower()


@pytest.mark.logic
def test_several_failures_are_all_reported():
    """The flags are the diagnosis. Stopping at the first would hide the
    rest."""
    row = _filtered(sc="0.1", bsa="10")
    assert "sc_too_low" in row["orthogonal_flags"]
    assert "bsa_too_low" in row["orthogonal_flags"]


# ── Summary loading ──────────────────────────────────────────────────

@pytest.mark.logic
def test_summaries_are_keyed_by_sequence_name(tmp_path):
    a = tmp_path / "a.csv"
    a.write_text("seq_name,sc\ndesign_1_seq_1,0.7\ndesign_2_seq_1,0.4\n")
    got = mom._load_summaries(str(tmp_path / "*.csv"))
    assert set(got) == {"design_1_seq_1", "design_2_seq_1"}
    assert got["design_1_seq_1"]["sc"] == "0.7"


@pytest.mark.logic
def test_no_matching_files_gives_an_empty_mapping(tmp_path):
    assert mom._load_summaries(str(tmp_path / "nothing_*.csv")) == {}


# ── The Rosetta score file ───────────────────────────────────────────

@pytest.mark.logic
def test_a_score_file_parses_to_header_and_rows(tmp_path):
    p = tmp_path / "score.sc"
    p.write_text(
        "SEQUENCE: AAAA\n"
        "SCORE: total_score sc_value dG_separated description\n"
        "SCORE: -100.0 0.600 -40.0 design_0\n"
        "SCORE: -110.0 0.650 -45.0 design_1\n")
    header, rows = rfc.parse_rosetta_scorefile(p)
    assert header[0] == "total_score"
    assert len(rows) == 2
    assert rows[1]["description"] == "design_1"
    assert rows[0]["sc_value"] == "0.600"


@pytest.mark.logic
def test_a_short_row_pairs_only_what_is_present(tmp_path):
    """Rosetta sometimes writes fewer fields than headers. zip stops at the
    shorter of the two, which is deliberate: the alternative is pairing a
    value with the wrong column name."""
    p = tmp_path / "score.sc"
    p.write_text(
        "SCORE: total_score sc_value dG_separated description\n"
        "SCORE: -100.0 0.600\n")
    _header, rows = rfc.parse_rosetta_scorefile(p)
    assert rows[0] == {"total_score": "-100.0", "sc_value": "0.600"}
    assert "description" not in rows[0]


@pytest.mark.logic
def test_sequence_lines_are_skipped(tmp_path):
    p = tmp_path / "score.sc"
    p.write_text("SEQUENCE: AAAA\nSEQUENCE: BBBB\n"
                 "SCORE: total_score description\nSCORE: -1.0 d0\n")
    header, rows = rfc.parse_rosetta_scorefile(p)
    assert header == ["total_score", "description"]
    assert len(rows) == 1


@pytest.mark.logic
def test_a_file_with_only_a_header_has_no_rows(tmp_path):
    p = tmp_path / "score.sc"
    p.write_text("SCORE: total_score description\n")
    header, rows = rfc.parse_rosetta_scorefile(p)
    assert header == ["total_score", "description"]
    assert rows == []


@pytest.mark.logic
@pytest.mark.parametrize("raw", ["", "NA", None, "not_a_number"])
def test_an_unparsable_score_becomes_nan_not_zero(raw):
    """nan fails every `< threshold` comparison downstream, so it cannot be
    mistaken for a good score. 0.0 would read as a real measurement."""
    import math
    assert math.isnan(rfc.safe_float(raw))


@pytest.mark.logic
def test_a_real_score_is_read():
    assert rfc.safe_float("-59.739") == pytest.approx(-59.739)


@pytest.mark.logic
def test_an_explicit_default_is_honoured():
    assert rfc.safe_float("junk", default=0.0) == 0.0
