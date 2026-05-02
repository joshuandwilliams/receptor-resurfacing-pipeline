"""Stage 3 — Rosetta filtering.

Pins the outputs of ``modules/rosetta_filtering.nf``::ROSETTA_FILTER, which
runs ``bin/rosetta_filter.py``.
"""
from __future__ import annotations

import pytest

from tests.characterization.helpers.json_compare import compare_json_deep
from tests.characterization.helpers.text_compare import compare_text_exact


@pytest.mark.hpc
@pytest.mark.wave1
def test_rosetta_filter_metrics_from_rosetta_filter_py(reference_root, output_root):
    """bin/rosetta_filter.py — JSON-DEEP."""
    rel = "rosetta_filtering/rosetta_filter_metrics.json"
    compare_json_deep(reference_root / rel, output_root / rel).assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
def test_rosetta_filter_summary_from_rosetta_filter_py(reference_root, output_root):
    """bin/rosetta_filter.py — JSON-DEEP."""
    rel = "rosetta_filtering/rosetta_filter_summary.json"
    compare_json_deep(reference_root / rel, output_root / rel).assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
def test_rosetta_passing_designs_from_rosetta_filter_py(reference_root, output_root):
    """bin/rosetta_filter.py — TEXT-EXACT."""
    rel = "rosetta_filtering/rosetta_passing_designs.txt"
    compare_text_exact(reference_root / rel, output_root / rel).assert_passed()
