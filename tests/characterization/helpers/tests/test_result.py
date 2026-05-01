from __future__ import annotations

from pathlib import Path

import pytest

from tests.characterization.helpers.result import ComparisonResult


@pytest.mark.local_unit
def test_comparison_result_truthy_when_passed():
    """bool(result) is True when passed=True."""
    r = ComparisonResult(passed=True, reference=Path("/r"), actual=Path("/a"), strategy="X")
    assert bool(r) is True
    assert r


@pytest.mark.local_unit
def test_comparison_result_falsy_when_failed():
    """bool(result) is False when passed=False."""
    r = ComparisonResult(passed=False, reference=Path("/r"), actual=Path("/a"), strategy="X")
    assert bool(r) is False
    assert not r


@pytest.mark.local_unit
def test_comparison_result_assert_passed_noop_on_pass():
    """assert_passed() on a passed result returns without raising."""
    r = ComparisonResult(passed=True, reference=Path("/r"), actual=Path("/a"), strategy="X")
    r.assert_passed()


@pytest.mark.local_unit
def test_comparison_result_assert_passed_raises_with_full_message_on_failure():
    """assert_passed() on a failed result raises with strategy, both paths, message, and differences."""
    r = ComparisonResult(
        passed=False,
        reference=Path("/refpath"),
        actual=Path("/actpath"),
        strategy="MY-STRATEGY",
        message="things went wrong",
        differences=["row 0 col x: 1 vs 2", "row 3 col y: nan vs 4"],
    )
    with pytest.raises(AssertionError) as exc_info:
        r.assert_passed()
    msg = str(exc_info.value)
    assert "MY-STRATEGY" in msg
    assert "/refpath" in msg
    assert "/actpath" in msg
    assert "things went wrong" in msg
    assert "row 0 col x: 1 vs 2" in msg
    assert "row 3 col y: nan vs 4" in msg
