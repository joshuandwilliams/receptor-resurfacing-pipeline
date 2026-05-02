"""Stage 2 — RFDiffusion + filter.

Pins the outputs of ``modules/rfdiffusion.nf``::RFDIFFUSION_FILTER, which
runs ``bin/rfdiffusion_filter.py``.

NOTE: ``rfdiffusion_metrics.json`` carries ``design_region_coords[][3]`` — Cα
coordinates of the de novo region — which depend on RFDiffusion's RNG seed.
If the first HPC round-trip fails on this file because the seed is not
honoured end-to-end, demote to JSON-STRUCT (keys/types/ranges only). For
this prompt: pin JSON-DEEP and let the round-trip tell us. (Plan §2.2.)
"""
from __future__ import annotations

import pytest

from tests.characterization.helpers.json_compare import compare_json_deep
from tests.characterization.helpers.text_compare import compare_text_exact


@pytest.mark.hpc
@pytest.mark.wave1
def test_rfdiffusion_metrics_from_rfdiffusion_filter_py(reference_root, output_root):
    """bin/rfdiffusion_filter.py — JSON-DEEP."""
    rel = "rfdiffusion/rfdiffusion_metrics.json"
    compare_json_deep(reference_root / rel, output_root / rel).assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
def test_filter_summary_from_rfdiffusion_filter_py(reference_root, output_root):
    """bin/rfdiffusion_filter.py — JSON-DEEP."""
    rel = "rfdiffusion/filter_summary.json"
    compare_json_deep(reference_root / rel, output_root / rel).assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
def test_passing_designs_from_rfdiffusion_filter_py(reference_root, output_root):
    """bin/rfdiffusion_filter.py — TEXT-EXACT."""
    rel = "rfdiffusion/passing_designs.txt"
    compare_text_exact(reference_root / rel, output_root / rel).assert_passed()
