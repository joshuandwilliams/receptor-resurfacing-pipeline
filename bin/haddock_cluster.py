"""
haddock_cluster.py
------------------
HaddockCluster — per-cluster representative for one HADDOCK3 docking run.

Per Session 7 grill-me (notes/design_audit.md Q153–Q155).  Composed
inside HaddockRun (one HaddockCluster per qualifying cluster).  All
metric fields are flat scalars; no nested Metrics sub-object (matches
OrthogonalMetrics style).

Some metric fields are Optional because they depend on inputs unavailable
at HaddockRun.from_workdir construction time:
- ``sc`` is None when Rosetta InterfaceAnalyzer has not been run yet.
- ``bridge_distance_*`` are None when the receptor contig is unavailable
  (per Session 7 Q152(i): bridge distance is filled by the downstream
  Nextflow wrapper that has access to both HaddockRun and the contig).
- ``pair_contact_fraction`` is reported as 1.0 (vacuous) when
  ``pair_total_count == 0`` so the auto-pick lexicographic sort
  ``(-pair_contact_fraction, -bsa)`` still produces a well-defined order
  with no contact-pair restraints.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Tuple


@dataclass(frozen=True)
class HaddockCluster:
    """One qualifying HADDOCK3 cluster, summarised by its best model.

    "Best model" within a cluster means the lowest HADDOCK-score member
    (carried as ``mean_haddock_score`` at the cluster level for
    diagnostics; the per-model score is not stored here).
    """
    cluster_id: int
    size: int                                 # members in the cluster
    mean_haddock_score: float                 # informational only
    best_model_pdb: Path                      # best_cluster{id}.pdb

    # ── Ranking + reporting metrics ─────────────────────────────────
    bsa: float                                # buried surface area, Å²
    com_distance: float                       # CA-centroid distance, Å
    air_satisfaction_count: int               # AIRs satisfied on best model
    air_total_count: int                      # AIRs in restraints file
    pair_contact_fraction: float              # 0.0–1.0; 1.0 when no pairs
    pair_total_count: int                     # user-defined pairs (denominator)
    clashes_in_design_region: int             # heavy-atom pairs <2 Å, in set
    clashes_outside_design_region: int        # everything else

    # ── Optional metrics filled post-construction ───────────────────
    sc: Optional[float] = None                # shape complementarity
    bridge_distance_per_anchor: Optional[Tuple[float, ...]] = None
    bridge_distance_min: Optional[float] = None
    bridge_distance_max: Optional[float] = None

    def __post_init__(self) -> None:
        if self.cluster_id < 1:
            raise ValueError(
                f"cluster_id must be >= 1, got {self.cluster_id}"
            )
        if self.size < 1:
            raise ValueError(f"size must be >= 1, got {self.size}")
        if not 0.0 <= self.pair_contact_fraction <= 1.0:
            raise ValueError(
                f"pair_contact_fraction must be in [0, 1], "
                f"got {self.pair_contact_fraction}"
            )
        if self.pair_total_count < 0:
            raise ValueError("pair_total_count must be >= 0")
        if self.air_total_count < 0:
            raise ValueError("air_total_count must be >= 0")
        if self.air_satisfaction_count > self.air_total_count:
            raise ValueError(
                f"air_satisfaction_count {self.air_satisfaction_count} > "
                f"air_total_count {self.air_total_count}"
            )
        if self.clashes_in_design_region < 0 or self.clashes_outside_design_region < 0:
            raise ValueError("clash counts must be >= 0")
        # Bridge distance fields move together: either all three filled
        # or all three None.  Catches the half-populated case.
        bd_fields = (
            self.bridge_distance_per_anchor,
            self.bridge_distance_min,
            self.bridge_distance_max,
        )
        n_set = sum(1 for f in bd_fields if f is not None)
        if n_set not in (0, 3):
            raise ValueError(
                "bridge_distance_per_anchor / _min / _max must all be set "
                "or all be None"
            )

    @property
    def auto_pick_key(self) -> Tuple[float, float]:
        """Sort key for the auto-pick when ``stop_after_haddock=false``
        and ``haddock_chosen_cluster`` is unset.

        Per Session 7 A154 — lexicographic ``(-pair_contact_fraction, -bsa)``.
        ``min(clusters, key=auto_pick_key)`` picks the cluster with the
        highest pair contact fraction; ties on pair fraction are broken
        by highest BSA.
        """
        return (-self.pair_contact_fraction, -self.bsa)

    def with_sc(self, sc: float) -> "HaddockCluster":
        """Return a copy with ``sc`` filled in (immutable update)."""
        return self.__class__(
            cluster_id=self.cluster_id,
            size=self.size,
            mean_haddock_score=self.mean_haddock_score,
            best_model_pdb=self.best_model_pdb,
            bsa=self.bsa,
            com_distance=self.com_distance,
            air_satisfaction_count=self.air_satisfaction_count,
            air_total_count=self.air_total_count,
            pair_contact_fraction=self.pair_contact_fraction,
            pair_total_count=self.pair_total_count,
            clashes_in_design_region=self.clashes_in_design_region,
            clashes_outside_design_region=self.clashes_outside_design_region,
            sc=sc,
            bridge_distance_per_anchor=self.bridge_distance_per_anchor,
            bridge_distance_min=self.bridge_distance_min,
            bridge_distance_max=self.bridge_distance_max,
        )

    def with_bridge_distances(
        self, per_anchor: Tuple[float, ...]
    ) -> "HaddockCluster":
        """Return a copy with bridge distance fields filled from a
        per-anchor tuple.  Min and max are derived automatically.
        """
        if not per_anchor:
            raise ValueError(
                "with_bridge_distances: per_anchor must be non-empty; "
                "for empty case keep the field None"
            )
        return self.__class__(
            cluster_id=self.cluster_id,
            size=self.size,
            mean_haddock_score=self.mean_haddock_score,
            best_model_pdb=self.best_model_pdb,
            bsa=self.bsa,
            com_distance=self.com_distance,
            air_satisfaction_count=self.air_satisfaction_count,
            air_total_count=self.air_total_count,
            pair_contact_fraction=self.pair_contact_fraction,
            pair_total_count=self.pair_total_count,
            clashes_in_design_region=self.clashes_in_design_region,
            clashes_outside_design_region=self.clashes_outside_design_region,
            sc=self.sc,
            bridge_distance_per_anchor=tuple(per_anchor),
            bridge_distance_min=min(per_anchor),
            bridge_distance_max=max(per_anchor),
        )

    def to_summary_row(self) -> dict:
        """Flat dict for cluster_metrics.csv / haddock_report.json output."""
        return {
            "cluster_id": self.cluster_id,
            "size": self.size,
            "mean_haddock_score": round(self.mean_haddock_score, 3),
            "best_model_pdb": self.best_model_pdb.name,
            "bsa": round(self.bsa, 2),
            "com_distance": round(self.com_distance, 3),
            "air_satisfaction_count": self.air_satisfaction_count,
            "air_total_count": self.air_total_count,
            "pair_contact_fraction": round(self.pair_contact_fraction, 3),
            "pair_total_count": self.pair_total_count,
            "clashes_in_design_region": self.clashes_in_design_region,
            "clashes_outside_design_region": self.clashes_outside_design_region,
            "sc": None if self.sc is None else round(self.sc, 4),
            "bridge_distance_min": (
                None if self.bridge_distance_min is None
                else round(self.bridge_distance_min, 3)
            ),
            "bridge_distance_max": (
                None if self.bridge_distance_max is None
                else round(self.bridge_distance_max, 3)
            ),
            "bridge_distance_per_anchor": (
                None if self.bridge_distance_per_anchor is None
                else [round(v, 3) for v in self.bridge_distance_per_anchor]
            ),
        }
