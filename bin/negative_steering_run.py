"""
negative_steering_run.py
------------------------
NegativeSteeringRun — one MPNN sequence's full negsteer experiment.

Per Phase 4 spec §2.9.  Tier 5 type — depends on StageResult,
DesignedSequence, PipelineInternalThresholds.

Read-only view: constructed from already-completed stages (the
orchestration stays in `negative_steering_run_one.sh` + the Nextflow
process body per Q10).  Owns the cross-stage per-seed verdict
aggregation and the outcome / n_pass / tier decisions.
"""

from __future__ import annotations

import math
import sys
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

from designed_sequence import DesignedSequence  # noqa: E402
from position_set import PositionSet  # noqa: E402
from protein_structure_prediction import ProteinStructurePrediction  # noqa: E402
from stage_result import StageResult  # noqa: E402


VALID_ROW_TYPES = frozenset(
    ("steered", "control_scrambled", "control_polyA")
)


@dataclass(frozen=True)
class NegativeSteeringRun:
    """One MPNN sequence's complete negsteer chain.

    Tier 5 constructor takes pre-built StageResults; the workdir-
    walking factory (``from_workdir``) is implemented at migration time
    when the existing parsers in cross_sequence_summary.py are absorbed.
    """
    mpnn_sequence_id: str
    workdir: Path
    row_type: str
    cold_start: StageResult
    steered: Dict[str, StageResult] = field(default_factory=dict)
    reversion: Dict[str, StageResult] = field(default_factory=dict)
    truth: Optional[ProteinStructurePrediction] = None
    contamination_positions: Optional[PositionSet] = None
    designed_sequence: Optional[DesignedSequence] = None
    run_one_runtime_sec: Optional[float] = None

    def __post_init__(self) -> None:
        if not isinstance(self.mpnn_sequence_id, str) or not self.mpnn_sequence_id:
            raise ValueError("mpnn_sequence_id must be a non-empty string")
        if self.row_type not in VALID_ROW_TYPES:
            raise ValueError(
                f"row_type {self.row_type!r} not in {sorted(VALID_ROW_TYPES)}"
            )
        if not isinstance(self.cold_start, StageResult):
            raise TypeError("cold_start must be a StageResult")
        if self.cold_start.stage_type != "cold_start":
            raise ValueError(
                f"cold_start stage must have stage_type='cold_start', "
                f"got {self.cold_start.stage_type!r}"
            )
        for design_id, st in self.steered.items():
            if not isinstance(st, StageResult):
                raise TypeError(
                    f"steered[{design_id!r}] must be a StageResult"
                )
            if st.stage_type != "steered":
                raise ValueError(
                    f"steered[{design_id!r}] must have stage_type='steered'"
                )
        for rev_id, st in self.reversion.items():
            if not isinstance(st, StageResult):
                raise TypeError(
                    f"reversion[{rev_id!r}] must be a StageResult"
                )
            if st.stage_type != "reversion":
                raise ValueError(
                    f"reversion[{rev_id!r}] must have stage_type='reversion'"
                )

    @property
    def num_seeds(self) -> int:
        return self.cold_start.num_seeds

    # ── Reversion mapping (per spec §2.9 + Step 4.3 addition) ──────

    def _reversion_for_design(self, design_id: str) -> Optional[StageResult]:
        """Walk reversion stages, return the one whose
        ``applies_to_designs`` contains ``design_id``, or None."""
        for st in self.reversion.values():
            if st.applies_to_designs and design_id in st.applies_to_designs:
                return st
        return None

    # ── Per-seed cross-stage aggregation ────────────────────────────

    def per_seed_final_verdicts(self, thresholds) -> List[Dict]:
        """Final verdict per (design, seed) across all stages.

        Rule (per spec §2.9):
        - If steering didn't run (no `steered` entries), use the
          cold-start stage's per-seed verdicts.  This is the cold-start
          path; outcome will be ``no_reversion``.
        - Otherwise, for each (design, seed) in `steered`:
            - If reversion ran for that design (matched by
              applies_to_designs), use the matching reversion verdict
              (paired by seed_index — current convention).
            - Otherwise, use the steered verdict.

        Returns list of dicts: ``{design_id, seed_index, stage_run,
        verdict}``.
        """
        if self._truth_required():
            self._ensure_truth()

        # Cold-start path
        if not self.steered:
            verdicts = self.cold_start.per_seed_verdicts(thresholds, self.truth)
            return [
                {
                    "design_id": "initial",
                    "seed_index": si,
                    "stage_run": "cold_start",
                    "verdict": v,
                }
                for si, v in zip(self.cold_start.seed_indices, verdicts)
            ]

        # Normal path: steered (and possibly reversion)
        out: List[Dict] = []
        for design_id, steered_stage in self.steered.items():
            rev_stage = self._reversion_for_design(design_id)
            if rev_stage is not None:
                # Pair by seed_index — current convention
                rev_verdicts = rev_stage.per_seed_verdicts(
                    thresholds, self.truth, self.contamination_positions,
                )
                rev_by_seed = dict(zip(rev_stage.seed_indices, rev_verdicts))
                for si in steered_stage.seed_indices:
                    if si in rev_by_seed:
                        out.append({
                            "design_id": design_id,
                            "seed_index": si,
                            "stage_run": "reverted",
                            "verdict": rev_by_seed[si],
                        })
                    else:
                        # Reversion missing this seed — keep the steered verdict
                        steered_verdicts = steered_stage.per_seed_verdicts(
                            thresholds, self.truth, self.contamination_positions,
                        )
                        st_by_seed = dict(
                            zip(steered_stage.seed_indices, steered_verdicts)
                        )
                        out.append({
                            "design_id": design_id,
                            "seed_index": si,
                            "stage_run": "steered",
                            "verdict": st_by_seed.get(si, "no_data"),
                        })
            else:
                # No reversion → steered verdict stands
                steered_verdicts = steered_stage.per_seed_verdicts(
                    thresholds, self.truth, self.contamination_positions,
                )
                for si, v in zip(steered_stage.seed_indices, steered_verdicts):
                    out.append({
                        "design_id": design_id,
                        "seed_index": si,
                        "stage_run": "steered",
                        "verdict": v,
                    })
        return out

    def _truth_required(self) -> bool:
        # Truth always needed for verdicts.
        return True

    def _ensure_truth(self) -> None:
        if self.truth is None:
            raise ValueError(
                "NegativeSteeringRun verdict methods require a 'truth' "
                "ProteinStructurePrediction at construction"
            )

    # ── Outcome + n_pass + tier ─────────────────────────────────────

    PASS_EQUIVALENT_VERDICTS = frozenset(
        ("clean", "clean_steered", "pose_holds")
    )

    def n_pass(self, thresholds) -> int:
        """Count seeds whose final verdict is pass-equivalent."""
        verdicts = self.per_seed_final_verdicts(thresholds)
        return sum(
            1 for e in verdicts if e["verdict"] in self.PASS_EQUIVALENT_VERDICTS
        )

    def n_seeds_total(self) -> int:
        """Total seed count for this MPNN sequence.

        For the cold-start path, this is the cold-start stage's
        num_seeds.  For the steered path, it's `n_designs × num_seeds`
        summed across all steered designs (matches the existing per-
        design × per-seed accounting in the codebase).
        """
        if not self.steered:
            return self.cold_start.num_seeds
        return sum(st.num_seeds for st in self.steered.values())

    def outcome(self, thresholds) -> str:
        """Aggregate outcome label per spec §2.9.

        Returns one of: ``no_reversion`` | ``pose_holds`` |
        ``pose_collapses`` | ``new_contamination`` | ``singleton``.

        Logic mirrors the existing `_classify_outcome` rules with the
        post-rename labels.
        """
        if not self.steered:
            # No steering ran — cold-start path
            return "no_reversion"
        # Did any reversion run?
        if not self.reversion:
            return "no_reversion"  # steering ran, no contamination triggered reversion
        verdicts = self.per_seed_final_verdicts(thresholds)
        counts = Counter(e["verdict"] for e in verdicts)
        nph = counts.get("pose_holds", 0)
        npc = counts.get("pose_collapses", 0)
        nnc = counts.get("new_contamination", 0)
        ncs = counts.get("clean_steered", 0)
        # Per the existing logic: if any pose_holds AND failure modes
        # don't dominate, return pose_holds.  Otherwise pick by
        # plurality between failure modes.
        nph_eff = nph + ncs
        if nph_eff > 0 and nph_eff >= max(npc, nnc):
            return "pose_holds"
        if nnc > npc:
            return "new_contamination"
        return "pose_collapses"

    def outcome_reason(self, thresholds) -> str:
        """Short diagnostic explanation of the outcome decision."""
        if not self.steered:
            return "cold-start path: no steering ran"
        if not self.reversion:
            return ("steering ran but no design triggered reversion under "
                    "the CL-3 majority rule")
        verdicts = self.per_seed_final_verdicts(thresholds)
        counts = Counter(e["verdict"] for e in verdicts)
        parts = [f"{k}={v}" for k, v in sorted(counts.items())]
        return "; ".join(parts)

    def tier(self, thresholds) -> str:
        """Tier classification (A / B / C / none) from n_pass / n_seeds.

        Per the post-rename rule (no `outcome=='no_reversion'` shortcut):
            A: n_pass == n_seeds (e.g. 3/3)
            B: 1 < n_pass < n_seeds (e.g. 2/3)
            C: n_pass == 1
            none: 0
        """
        n_pass = self.n_pass(thresholds)
        n_seeds = self.n_seeds_total()
        if n_seeds <= 0:
            return "none"
        if n_pass == n_seeds:
            return "A"
        if 1 < n_pass < n_seeds:
            return "B"
        if n_pass == 1:
            return "C"
        return "none"

    # ── Identification / passthrough ────────────────────────────────

    def __repr__(self) -> str:
        return (
            f"NegativeSteeringRun(mpnn_sequence_id="
            f"{self.mpnn_sequence_id!r}, row_type={self.row_type!r}, "
            f"n_steered_designs={len(self.steered)}, "
            f"n_reversions={len(self.reversion)})"
        )
