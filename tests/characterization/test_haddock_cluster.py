"""Unit tests for bin/haddock_cluster.py."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO_ROOT / "bin"))

import haddock_cluster as hc  # noqa: E402


def _make_cluster(**overrides):
    defaults = dict(
        cluster_id=1,
        size=5,
        mean_haddock_score=-50.0,
        best_model_pdb=Path("best_cluster1.pdb"),
        bsa=800.0,
        com_distance=15.5,
        air_satisfaction_count=10,
        air_total_count=20,
        pair_contact_fraction=1.0,
        pair_total_count=3,
        clashes_in_design_region=0,
        clashes_outside_design_region=0,
    )
    defaults.update(overrides)
    return hc.HaddockCluster(**defaults)


@pytest.mark.local_unit
class TestConstruction:
    def test_basic(self):
        c = _make_cluster()
        assert c.cluster_id == 1
        assert c.bsa == 800.0
        assert c.sc is None
        assert c.bridge_distance_min is None

    def test_bad_cluster_id_raises(self):
        with pytest.raises(ValueError, match="cluster_id"):
            _make_cluster(cluster_id=0)

    def test_bad_size_raises(self):
        with pytest.raises(ValueError, match="size"):
            _make_cluster(size=0)

    def test_pair_fraction_out_of_range(self):
        with pytest.raises(ValueError, match="pair_contact_fraction"):
            _make_cluster(pair_contact_fraction=1.5)

    def test_air_sat_greater_than_total_raises(self):
        with pytest.raises(ValueError, match="air_satisfaction_count"):
            _make_cluster(air_satisfaction_count=30, air_total_count=20)

    def test_partial_bridge_distance_raises(self):
        with pytest.raises(ValueError, match="bridge_distance"):
            _make_cluster(
                bridge_distance_min=3.0,
                # _max and _per_anchor still None
            )


@pytest.mark.local_unit
class TestAutoPickKey:
    def test_lexicographic_pair_first(self):
        # Higher pair fraction should sort lower (earlier) by min().
        c1 = _make_cluster(cluster_id=1, pair_contact_fraction=2/3, bsa=900.0)
        c2 = _make_cluster(cluster_id=2, pair_contact_fraction=3/3, bsa=750.0)
        # auto_pick_key smaller = "better"; c2 has higher pair fraction
        # so its key is smaller.
        assert c2.auto_pick_key < c1.auto_pick_key

    def test_bsa_breaks_ties(self):
        # Same pair fraction → BSA decides.
        c1 = _make_cluster(cluster_id=1, pair_contact_fraction=1.0, bsa=800.0)
        c2 = _make_cluster(cluster_id=2, pair_contact_fraction=1.0, bsa=850.0)
        assert c2.auto_pick_key < c1.auto_pick_key

    def test_no_pairs_collapses_to_bsa(self):
        # pair_total_count == 0 → caller should set pair_contact_fraction=1.0
        # (the vacuous case) so the lex sort still works on BSA only.
        c1 = _make_cluster(pair_total_count=0, pair_contact_fraction=1.0, bsa=700.0)
        c2 = _make_cluster(pair_total_count=0, pair_contact_fraction=1.0, bsa=800.0)
        assert c2.auto_pick_key < c1.auto_pick_key


@pytest.mark.local_unit
class TestImmutableUpdates:
    def test_with_sc(self):
        c = _make_cluster()
        c2 = c.with_sc(0.75)
        assert c.sc is None
        assert c2.sc == 0.75
        assert c2.cluster_id == c.cluster_id

    def test_with_bridge_distances(self):
        c = _make_cluster()
        c2 = c.with_bridge_distances((3.5, 8.1, 12.7))
        assert c.bridge_distance_per_anchor is None
        assert c2.bridge_distance_per_anchor == (3.5, 8.1, 12.7)
        assert c2.bridge_distance_min == 3.5
        assert c2.bridge_distance_max == 12.7

    def test_with_bridge_distances_empty_raises(self):
        c = _make_cluster()
        with pytest.raises(ValueError, match="non-empty"):
            c.with_bridge_distances(())


@pytest.mark.local_unit
class TestSummaryRow:
    def test_basic_fields(self):
        c = _make_cluster()
        row = c.to_summary_row()
        assert row["cluster_id"] == 1
        assert row["bsa"] == 800.0
        assert row["sc"] is None
        assert row["bridge_distance_min"] is None
        assert row["pair_contact_fraction"] == 1.0

    def test_with_optionals_filled(self):
        c = _make_cluster().with_sc(0.62).with_bridge_distances((2.5, 7.5))
        row = c.to_summary_row()
        assert row["sc"] == 0.62
        assert row["bridge_distance_min"] == 2.5
        assert row["bridge_distance_max"] == 7.5
        assert row["bridge_distance_per_anchor"] == [2.5, 7.5]
