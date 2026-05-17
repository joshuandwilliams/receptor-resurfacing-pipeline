"""Unit tests for bin/haddock_run.py."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO_ROOT / "bin"))

import haddock_cluster as hc  # noqa: E402
import haddock_run as hr  # noqa: E402


# ── Parsers ──────────────────────────────────────────────────────────


@pytest.mark.local_unit
class TestParseContactPairs:
    def test_basic(self):
        assert hr.parse_contact_pairs("A25-C42 A13-C94") == (
            ("A", 25, "C", 42),
            ("A", 13, "C", 94),
        )

    def test_empty(self):
        assert hr.parse_contact_pairs("") == ()
        assert hr.parse_contact_pairs("   ") == ()

    def test_lowercase_chain_uppercased(self):
        assert hr.parse_contact_pairs("a25-c42") == (("A", 25, "C", 42),)

    def test_bad_format_raises(self):
        with pytest.raises(ValueError, match="bad pair"):
            hr.parse_contact_pairs("A25-42")
        with pytest.raises(ValueError, match="bad pair"):
            hr.parse_contact_pairs("25-C42")


@pytest.mark.local_unit
class TestParseActiveResidues:
    def test_singles(self):
        assert hr.parse_active_residues("25,35,40") == (25, 35, 40)

    def test_range(self):
        assert hr.parse_active_residues("40-44") == (40, 41, 42, 43, 44)

    def test_mixed(self):
        assert hr.parse_active_residues("25,35,40-44") == (25, 35, 40, 41, 42, 43, 44)

    def test_empty(self):
        assert hr.parse_active_residues("") == ()

    def test_dedup_and_sort(self):
        assert hr.parse_active_residues("5,3,3,1-2") == (1, 2, 3, 5)

    def test_bad_range_raises(self):
        with pytest.raises(ValueError, match="lo > hi"):
            hr.parse_active_residues("10-5")


@pytest.mark.local_unit
class TestParsePairDistance:
    def test_basic(self):
        assert hr.parse_pair_distance("2,2,4") == (2.0, 2.0, 4.0)

    def test_floats(self):
        assert hr.parse_pair_distance("1.5,0.5,3.5") == (1.5, 0.5, 3.5)

    def test_wrong_count(self):
        with pytest.raises(ValueError, match="3 comma-separated"):
            hr.parse_pair_distance("2,2")

    def test_negative_rejected(self):
        with pytest.raises(ValueError, match=">= 0"):
            hr.parse_pair_distance("-1,2,4")


# ── Helpers ──────────────────────────────────────────────────────────


def _make_cluster(cluster_id=1, **overrides):
    defaults = dict(
        cluster_id=cluster_id,
        size=5,
        mean_haddock_score=-50.0,
        best_model_pdb=Path(f"best_cluster{cluster_id}.pdb"),
        bsa=800.0,
        com_distance=15.5,
        air_satisfaction_count=10,
        air_total_count=20,
        pair_contact_fraction=1.0,
        pair_total_count=2,
        clashes_in_design_region=0,
        clashes_outside_design_region=0,
    )
    defaults.update(overrides)
    return hc.HaddockCluster(**defaults)


def _make_run(tmp_path, **overrides):
    receptor = tmp_path / "receptor.pdb"
    effector = tmp_path / "effector.pdb"
    receptor.write_text("")
    effector.write_text("")
    defaults = dict(
        receptor_pdb=receptor,
        effector_pdb=effector,
        contact_pairs=(("A", 25, "C", 42), ("A", 13, "C", 94)),
        receptor_active_residues=(35, 40, 41, 42),
        effector_active_residues=(),
        pair_distance=(2.0, 2.0, 4.0),
        qualifying_clusters=(_make_cluster(1, bsa=900.0),),
        chosen_cluster_id=None,
    )
    defaults.update(overrides)
    return hr.HaddockRun(**defaults)


# ── HaddockRun ───────────────────────────────────────────────────────


@pytest.mark.local_unit
class TestConstruction:
    def test_basic(self, tmp_path):
        run = _make_run(tmp_path)
        assert run.chosen_cluster_id is None
        assert run.pair_distance == (2.0, 2.0, 4.0)
        assert len(run.qualifying_clusters) == 1

    def test_invalid_chosen_cluster_raises(self, tmp_path):
        with pytest.raises(ValueError, match="chosen_cluster_id"):
            _make_run(tmp_path, chosen_cluster_id=99)

    def test_wrong_cluster_type_raises(self, tmp_path):
        with pytest.raises(TypeError, match="HaddockCluster"):
            _make_run(tmp_path, qualifying_clusters=("not a cluster",))

    def test_pair_distance_arity(self, tmp_path):
        with pytest.raises(ValueError, match="3-tuple"):
            _make_run(tmp_path, pair_distance=(2.0, 2.0))


@pytest.mark.local_unit
class TestDesignRegion:
    def test_union_of_pairs_and_active(self, tmp_path):
        run = _make_run(
            tmp_path,
            contact_pairs=(("A", 25, "C", 42),),
            receptor_active_residues=(35, 36, 37),
        )
        assert run.design_region == {25, 35, 36, 37}

    def test_no_restraints(self, tmp_path):
        run = _make_run(
            tmp_path,
            contact_pairs=(),
            receptor_active_residues=(),
        )
        assert run.design_region == set()
        assert run.has_restraints is False

    def test_has_restraints_true_with_pairs_only(self, tmp_path):
        run = _make_run(
            tmp_path,
            contact_pairs=(("A", 25, "C", 42),),
            receptor_active_residues=(),
        )
        assert run.has_restraints is True


@pytest.mark.local_unit
class TestSelected:
    def test_explicit_chosen_wins(self, tmp_path):
        clusters = (
            _make_cluster(1, bsa=900.0, pair_contact_fraction=1.0),
            _make_cluster(2, bsa=850.0, pair_contact_fraction=1.0),
        )
        run = _make_run(tmp_path, qualifying_clusters=clusters, chosen_cluster_id=2)
        assert run.selected.cluster_id == 2

    def test_auto_pick_prefers_pair_fraction(self, tmp_path):
        # cluster 1: BSA 900, pairs 2/3.  cluster 2: BSA 800, pairs 3/3.
        # Auto-pick should choose 2 (higher pair fraction).
        clusters = (
            _make_cluster(1, bsa=900.0, pair_contact_fraction=2/3, pair_total_count=3),
            _make_cluster(2, bsa=800.0, pair_contact_fraction=3/3, pair_total_count=3),
        )
        run = _make_run(tmp_path, qualifying_clusters=clusters)
        assert run.selected.cluster_id == 2

    def test_auto_pick_bsa_breaks_ties(self, tmp_path):
        clusters = (
            _make_cluster(1, bsa=750.0, pair_contact_fraction=1.0, pair_total_count=2),
            _make_cluster(2, bsa=900.0, pair_contact_fraction=1.0, pair_total_count=2),
        )
        run = _make_run(tmp_path, qualifying_clusters=clusters)
        assert run.selected.cluster_id == 2

    def test_no_clusters_raises(self, tmp_path):
        run = _make_run(tmp_path, qualifying_clusters=())
        with pytest.raises(ValueError, match="no qualifying"):
            _ = run.selected


@pytest.mark.local_unit
class TestToRfdiffusionInput:
    def test_returns_selected_pdb(self, tmp_path):
        clusters = (_make_cluster(1, bsa=800.0),)
        run = _make_run(tmp_path, qualifying_clusters=clusters)
        assert run.to_rfdiffusion_input().name == "best_cluster1.pdb"


# ── from_workdir factory ─────────────────────────────────────────────


@pytest.mark.local_unit
class TestFromWorkdir:
    def test_round_trip(self, tmp_path):
        receptor = tmp_path / "receptor.pdb"
        effector = tmp_path / "effector.pdb"
        receptor.write_text("")
        effector.write_text("")
        # Fake HADDOCK output fixture.
        report = {
            "success": True,
            "cluster_models": {
                "1": {
                    "filename": "best_cluster1.pdb",
                    "size": 6,
                    "mean_score": -45.5,
                    "model_name": "rigidbody_42",
                    "score": -52.1,
                },
                "2": {
                    "filename": "best_cluster2.pdb",
                    "size": 4,
                    "mean_score": -30.2,
                    "model_name": "rigidbody_88",
                    "score": -33.7,
                },
            },
        }
        metrics = {
            "1": {
                "bsa": 950.0,
                "com_distance": 14.5,
                "air_satisfaction_count": 18,
                "air_total_count": 20,
                "pair_contact_fraction": 1.0,
                "pair_total_count": 2,
                "clashes_in_design_region": 1,
                "clashes_outside_design_region": 0,
                "sc": 0.62,
            },
            "2": {
                "bsa": 700.0,
                "com_distance": 18.5,
                "air_satisfaction_count": 12,
                "air_total_count": 20,
                "pair_contact_fraction": 0.5,
                "pair_total_count": 2,
                "clashes_in_design_region": 0,
                "clashes_outside_design_region": 2,
                "sc": 0.55,
            },
        }
        (tmp_path / "haddock_report.json").write_text(json.dumps(report))
        (tmp_path / "cluster_metrics.json").write_text(json.dumps(metrics))
        # The best_cluster*.pdb files don't need real content for the factory.
        (tmp_path / "best_cluster1.pdb").write_text("")
        (tmp_path / "best_cluster2.pdb").write_text("")

        run = hr.HaddockRun.from_workdir(
            workdir=tmp_path,
            receptor_pdb=receptor,
            effector_pdb=effector,
            contact_pairs="A25-C42 A13-C94",
            receptor_active_residues="35,40-44",
            pair_distance="2,2,4",
        )
        assert len(run.qualifying_clusters) == 2
        assert run.contact_pairs == (("A", 25, "C", 42), ("A", 13, "C", 94))
        assert run.receptor_active_residues == (35, 40, 41, 42, 43, 44)
        # Auto-pick: cluster 1 has higher pair_contact_fraction.
        assert run.selected.cluster_id == 1
        assert run.selected.bsa == 950.0
        assert run.selected.sc == 0.62
        assert run.design_region == {25, 13, 35, 40, 41, 42, 43, 44}

    def test_missing_report_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="haddock_report"):
            hr.HaddockRun.from_workdir(
                workdir=tmp_path,
                receptor_pdb=tmp_path / "r.pdb",
                effector_pdb=tmp_path / "e.pdb",
            )
