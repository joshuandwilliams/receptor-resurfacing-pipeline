"""Stage 3 — Rosetta filtering.

Pins the outputs of ``modules/rosetta_filtering.nf``::ROSETTA_FILTER, which
runs ``bin/rosetta_filter.py``.

The per-module fixture curates four designs spanning the sc_threshold gate
(design_0 strong pass, design_1 marginal pass, design_14 boundary fail,
design_59 clear fail — see notes/inventory/15_discovery_run_path_coverage.md
§Rosetta filtering). The cohort-level metrics JSON contains all four
designs' rows, so the existing cohort-shaped tests cover the four-design
parametrisation implicitly.
"""
from __future__ import annotations

import pytest

from tests.characterization.helpers.json_compare import compare_json_deep
from tests.characterization.helpers.text_compare import compare_text_exact


@pytest.fixture
def ref(stage_reference_root):
    return stage_reference_root("rosetta_filtering")


@pytest.fixture
def out(stage_output_root):
    return stage_output_root("rosetta_filtering")


@pytest.mark.hpc
@pytest.mark.wave1
def test_rosetta_filter_metrics_from_rosetta_filter_py(ref, out):
    """bin/rosetta_filter.py — JSON-DEEP."""
    rel = "rosetta_filtering/rosetta_filter_metrics.json"
    compare_json_deep(ref / rel, out / rel).assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
def test_rosetta_filter_summary_from_rosetta_filter_py(ref, out):
    """bin/rosetta_filter.py — JSON-DEEP."""
    rel = "rosetta_filtering/rosetta_filter_summary.json"
    compare_json_deep(ref / rel, out / rel).assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
def test_rosetta_passing_designs_from_rosetta_filter_py(ref, out):
    """bin/rosetta_filter.py — TEXT-EXACT."""
    rel = "rosetta_filtering/rosetta_passing_designs.txt"
    compare_text_exact(ref / rel, out / rel).assert_passed()
