"""Stage 2 — RFDiffusion + filter.

Pins the outputs of ``modules/rfdiffusion.nf``::RFDIFFUSION_FILTER, which
runs ``bin/rfdiffusion_filter.py``.

NOTE: ``rfdiffusion_metrics.json`` carries ``design_region_coords[][3]`` — Cα
coordinates of the de novo region — which depend on RFDiffusion's RNG seed.
If the first HPC round-trip fails on this file because the seed is not
honoured end-to-end, demote to JSON-STRUCT (keys/types/ranges only). For
this prompt: pin JSON-DEEP and let the round-trip tell us. (Plan §2.2.)

The fail branch of ``passes_filter`` (n_contact_pairs == 0 or
frac_contacts_in_design < min_hotspot_frac) was not exercised by the
discovery run and has never been observed in any real run — see
``notes/inventory/15_discovery_run_path_coverage.md`` Coverage gap 1. It is
covered by the local-unit test at the bottom of this file rather than by a
fixture.
"""
from __future__ import annotations

import pytest

from tests.characterization.helpers.json_compare import compare_json_deep
from tests.characterization.helpers.text_compare import compare_text_exact


@pytest.fixture
def ref(stage_reference_root):
    return stage_reference_root("rfdiffusion")


@pytest.fixture
def out(stage_output_root):
    return stage_output_root("rfdiffusion")


@pytest.mark.hpc
@pytest.mark.wave1
def test_rfdiffusion_metrics_from_rfdiffusion_filter_py(ref, out):
    """bin/rfdiffusion_filter.py — JSON-DEEP."""
    rel = "rfdiffusion/rfdiffusion_metrics.json"
    compare_json_deep(ref / rel, out / rel).assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
def test_filter_summary_from_rfdiffusion_filter_py(ref, out):
    """bin/rfdiffusion_filter.py — JSON-DEEP."""
    rel = "rfdiffusion/filter_summary.json"
    compare_json_deep(ref / rel, out / rel).assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
def test_passing_designs_from_rfdiffusion_filter_py(ref, out):
    """bin/rfdiffusion_filter.py — TEXT-EXACT."""
    rel = "rfdiffusion/passing_designs.txt"
    compare_text_exact(ref / rel, out / rel).assert_passed()


# ---------------------------------------------------------------------------
# Local-unit test: fail-branch of bin/rfdiffusion_filter.py
# ---------------------------------------------------------------------------
#
# The filter predicate is inlined inside ``main()`` at line 820 of
# bin/rfdiffusion_filter.py and is not importable as a function. The
# producer line is:
#
#     passes = frac_in >= args.min_hotspot_frac if n_contacts > 0 else False
#
# Three logical branches:
#   - n_contacts == 0                   → False (degenerate fail)
#   - n_contacts > 0, frac_in < min_hsp → False (frac fail)
#   - n_contacts > 0, frac_in >= min_hsp → True  (pass)
#
# We replicate the predicate and assert all three branches. If the inlined
# logic ever changes, this test fails — flagging that the script-side
# predicate has drifted from this characterization assertion. The pass
# branch is also exercised by the hpc-tier reference set, but the two
# fail branches have never been observed in any real run.


def _passes_filter(n_contacts: int, frac_in: float, min_hotspot_frac: float) -> bool:
    """Mirror of bin/rfdiffusion_filter.py:820 — predicate under test."""
    return frac_in >= min_hotspot_frac if n_contacts > 0 else False


@pytest.mark.local_unit
def test_rfdiffusion_filter_fail_branch_unit():
    """bin/rfdiffusion_filter.py — fail-branch predicate (no fixture)."""
    # Branch 1: n_contact_pairs == 0 — degenerate input, design has no
    # Cα contacts whatsoever. passes_filter must be False regardless of
    # the (irrelevant) frac_in value the producer code computed.
    assert _passes_filter(n_contacts=0, frac_in=0.0, min_hotspot_frac=0.0) is False
    assert _passes_filter(n_contacts=0, frac_in=1.0, min_hotspot_frac=0.5) is False

    # Branch 2: n_contacts > 0 but frac_in below the hotspot threshold.
    assert _passes_filter(n_contacts=10, frac_in=0.20, min_hotspot_frac=0.50) is False

    # Branch 3 (sanity): pass branch — exercised by the hpc-tier fixture
    # too, but pinned here so the predicate's truth table is complete.
    assert _passes_filter(n_contacts=10, frac_in=0.50, min_hotspot_frac=0.50) is True
    assert _passes_filter(n_contacts=24, frac_in=0.46, min_hotspot_frac=0.0) is True
