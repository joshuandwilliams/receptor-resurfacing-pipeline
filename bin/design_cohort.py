"""
design_cohort.py
----------------
DesignCohort — the set of all MPNN sequences in one pipeline run,
each carrying its NegativeSteeringRun.

Per Phase 4 spec §2.10.  Tier 6 type — depends on NegativeSteeringRun,
OrthogonalMetrics (optional, attached separately), and
PipelineInternalThresholds.  Top of the hierarchy.

Knows tier breakdowns, the composite-score ranking across sequences,
and which sequences are survivors eligible for the orthogonal stage.
"""

from __future__ import annotations

import sys
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

from negative_steering_run import NegativeSteeringRun  # noqa: E402


@dataclass(frozen=True)
class DesignCohort:
    """All MPNN sequences from one pipeline run.

    Holds NegativeSteeringRun instances keyed by mpnn_sequence_id.
    Tier 6 constructor takes a list of runs directly; the workdir-
    walking factory (from_runs_directory) is implemented at migration
    time when the existing cross_sequence_summary.py absorbs.
    """
    runs: Tuple[NegativeSteeringRun, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.runs, tuple):
            object.__setattr__(self, "runs", tuple(self.runs))
        seen_ids = set()
        for r in self.runs:
            if not isinstance(r, NegativeSteeringRun):
                raise TypeError(
                    f"every run must be a NegativeSteeringRun, "
                    f"got {type(r).__name__}"
                )
            if r.mpnn_sequence_id in seen_ids:
                raise ValueError(
                    f"duplicate mpnn_sequence_id {r.mpnn_sequence_id!r}"
                )
            seen_ids.add(r.mpnn_sequence_id)

    # ── Lookup / filter ─────────────────────────────────────────────

    def __len__(self) -> int:
        return len(self.runs)

    def __iter__(self):
        return iter(self.runs)

    def get(self, mpnn_sequence_id: str) -> Optional[NegativeSteeringRun]:
        for r in self.runs:
            if r.mpnn_sequence_id == mpnn_sequence_id:
                return r
        return None

    def steered_runs(self) -> List[NegativeSteeringRun]:
        """Only steered rows (excludes controls)."""
        return [r for r in self.runs if r.row_type == "steered"]

    def control_runs(self) -> List[NegativeSteeringRun]:
        """Only the negative controls."""
        return [r for r in self.runs if r.row_type != "steered"]

    def by_tier(self, tier: str, thresholds) -> List[NegativeSteeringRun]:
        """Steered runs at the given tier (A/B/C/none)."""
        return [r for r in self.steered_runs() if r.tier(thresholds) == tier]

    def survivors(self, thresholds) -> List[NegativeSteeringRun]:
        """Tier A/B/C steered rows — the set fed to the orthogonal stage."""
        return [
            r for r in self.steered_runs()
            if r.tier(thresholds) in ("A", "B", "C")
        ]

    def tier_breakdown(self, thresholds) -> Dict[str, int]:
        """Per-tier count of STEERED runs (controls excluded)."""
        counts = Counter(
            r.tier(thresholds) for r in self.steered_runs()
        )
        return {
            "A": counts.get("A", 0),
            "B": counts.get("B", 0),
            "C": counts.get("C", 0),
            "none": counts.get("none", 0),
        }

    # ── Ranking ─────────────────────────────────────────────────────

    _TIER_ORDER = {"A": 0, "B": 1, "C": 2, "none": 3}

    def ranked_by_composite(
        self, thresholds
    ) -> List[NegativeSteeringRun]:
        """Sort steered runs by (tier_order, -composite_score).

        composite_score = `n_pass / n_seeds_total` (a simple ratio for
        Tier 6; the existing composite that weights jaccard − 0.05 ×
        ra_eff requires per-prediction metrics not yet accessible here
        and is computed in `to_cross_summary_csv` at migration time).
        """
        def key(r):
            tier = r.tier(thresholds)
            tier_idx = self._TIER_ORDER.get(tier, 99)
            n_seeds = r.n_seeds_total()
            ratio = r.n_pass(thresholds) / n_seeds if n_seeds > 0 else 0.0
            return (tier_idx, -ratio)
        return sorted(self.steered_runs(), key=key)

    # ── Stats ─────────────────────────────────────────────────────

    def n_steered(self) -> int:
        return sum(1 for r in self.runs if r.row_type == "steered")

    def n_controls(self) -> int:
        return sum(1 for r in self.runs if r.row_type != "steered")

    def __repr__(self) -> str:
        return (
            f"DesignCohort(n_steered={self.n_steered()}, "
            f"n_controls={self.n_controls()})"
        )
