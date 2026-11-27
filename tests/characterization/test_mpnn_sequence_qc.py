"""Cover the MPNN sequence QC predicates.

These two decide which designed sequences are thrown away before anything is
paid for in GPU time. Both fail in the expensive direction rather than the
loud one: a predicate that under-reports lets a degenerate sequence through to
structure prediction, and one that over-reports silently discards a design
that was fine.

The thresholds themselves are the contract, so the tests sit on the boundary
rather than well inside it.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
BIN = REPO / "bin"
if str(BIN) not in sys.path:
    sys.path.insert(0, str(BIN))

import mpnn_sequence_qc as qc  # noqa: E402

# ── Poly-X runs ──────────────────────────────────────────────────────

@pytest.mark.logic
def test_a_run_longer_than_the_limit_is_caught_and_the_residue_named():
    """Naming the residue is what makes the rejection actionable. 'failed QC'
    on its own does not say whether the sequence is poly-alanine or has one
    stray tract."""
    hit, aa = qc.has_poly_x("MKT" + "A" * 6 + "GG", max_run=5)
    assert hit is True
    assert aa == "A"


@pytest.mark.logic
def test_a_run_exactly_at_the_limit_is_allowed():
    """max_run is the longest run that passes, not the first that fails. The
    predicate tests for aa * (max_run + 1), so 5 A's pass at max_run=5."""
    hit, aa = qc.has_poly_x("MKT" + "A" * 5 + "GG", max_run=5)
    assert hit is False
    assert aa is None


@pytest.mark.logic
def test_one_over_the_limit_fails():
    assert qc.has_poly_x("A" * 6, max_run=5)[0] is True


@pytest.mark.logic
def test_runs_are_detected_for_any_residue_not_just_alanine():
    hit, aa = qc.has_poly_x("MK" + "Q" * 8 + "T", max_run=5)
    assert hit is True and aa == "Q"


@pytest.mark.logic
def test_a_run_split_by_another_residue_is_not_a_run():
    """AAAxAAA is two runs of three, not one of six. Counting the total
    occurrences instead would reject a perfectly ordinary sequence."""
    assert qc.has_poly_x("AAA" + "K" + "AAA", max_run=5)[0] is False


@pytest.mark.logic
def test_an_empty_sequence_has_no_run():
    assert qc.has_poly_x("", max_run=5) == (False, None)


# ── Composition ──────────────────────────────────────────────────────

@pytest.mark.logic
def test_a_non_standard_residue_is_rejected_and_named():
    """X, B, Z and U reach here from a parse that could not resolve a residue.
    Passing one through puts an unresolvable letter into a FASTA that a
    structure predictor will later reject or silently reinterpret."""
    hit, aa = qc.unusual_aa_composition("MKTAX" + "G" * 10)
    assert hit is True and aa == "X"


@pytest.mark.logic
def test_an_ordinary_sequence_passes():
    seq = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"
    assert qc.unusual_aa_composition(seq) == (False, None)


@pytest.mark.logic
def test_a_dominant_residue_above_the_fraction_is_rejected():
    """Half the sequence is tryptophan, which no designed surface should be."""
    seq = "W" * 10 + "MKTAYIAKQR"
    hit, aa = qc.unusual_aa_composition(seq, max_single_aa_frac=0.30)
    assert hit is True and aa == "W"


@pytest.mark.logic
@pytest.mark.parametrize("aa", ["A", "L", "G"])
def test_alanine_leucine_and_glycine_are_exempt_from_the_fraction(aa):
    """The exemption is deliberate. RFdiffusion writes glycine at every de
    novo position and the scaffold is alanine-rich, so applying the fraction
    to these three would reject essentially every design before MPNN has
    assigned a sequence.
    """
    seq = aa * 20 + "MKTYIQR"
    assert qc.unusual_aa_composition(seq, max_single_aa_frac=0.30)[0] is False


@pytest.mark.logic
def test_the_fraction_is_strict_so_exactly_at_the_threshold_passes():
    """The comparison is `frac > max_single_aa_frac`, so a residue at exactly
    the threshold is kept. 3 of 10 is 0.30 and passes."""
    seq = "KKK" + "MTYIQRV"
    assert len(seq) == 10
    assert qc.unusual_aa_composition(seq, max_single_aa_frac=0.30)[0] is False


@pytest.mark.logic
def test_just_over_the_fraction_fails():
    seq = "KKKK" + "MTYIQR"
    assert len(seq) == 10
    hit, aa = qc.unusual_aa_composition(seq, max_single_aa_frac=0.30)
    assert hit is True and aa == "K"


@pytest.mark.logic
def test_an_empty_sequence_does_not_divide_by_zero():
    assert qc.unusual_aa_composition("") == (False, None)


@pytest.mark.logic
def test_a_stricter_fraction_rejects_what_the_default_allows():
    """The flag has to actually reach the comparison."""
    seq = "KKK" + "MTYIQRV"
    assert qc.unusual_aa_composition(seq, max_single_aa_frac=0.30)[0] is False
    assert qc.unusual_aa_composition(seq, max_single_aa_frac=0.10)[0] is True
