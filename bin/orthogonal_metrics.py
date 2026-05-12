"""
orthogonal_metrics.py
---------------------
OrthogonalMetrics — the AF3 + biophysical + Rosetta measurements for
one cohort survivor.

Per Phase 4 spec §2.11.  Tier 4 type — depends on
ProteinStructurePrediction (Tier 2), AF3ConfidenceAggregate (Tier 0),
and PipelineInternalThresholds (Tier 0).

Independent of Boltz (the original prediction tool) so disagreement
between Boltz and these metrics flags potential prediction artefacts.
Strictly Sc + BSA + ΔΔG in the orthogonal_filters gate — `interface_plddt`
is Boltz-derived and explicitly NOT in the gate (per Session 5 Q118 and
commit `47cb9f2`).
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

from af3_confidence import AF3ConfidenceAggregate  # noqa: E402
from protein_structure_prediction import ProteinStructurePrediction  # noqa: E402


@dataclass(frozen=True)
class OrthogonalMetrics:
    """Orthogonal validation summary for one cohort survivor.

    Holds the AF3-no-MSA aggregate, the biophysical metrics (BSA,
    H-bonds), and the Rosetta metrics (Sc, ΔΔG).  Knows the
    passes_orthogonal_filters gate.
    """
    mpnn_sequence_id: str
    canonical_prediction: ProteinStructurePrediction
    af3: AF3ConfidenceAggregate
    bsa: Optional[float] = None
    hbonds: Optional[int] = None
    sc: Optional[float] = None
    ddg: Optional[float] = None

    def __post_init__(self) -> None:
        if not isinstance(self.mpnn_sequence_id, str) or not self.mpnn_sequence_id:
            raise ValueError("mpnn_sequence_id must be non-empty string")
        if not isinstance(self.canonical_prediction, ProteinStructurePrediction):
            raise TypeError(
                "canonical_prediction must be a ProteinStructurePrediction"
            )
        if not isinstance(self.af3, AF3ConfidenceAggregate):
            raise TypeError("af3 must be an AF3ConfidenceAggregate")

    # ── The orthogonal gate ──────────────────────────────────────────

    def passes_orthogonal_filters(self, thresholds) -> bool:
        """Sc + BSA + ΔΔG only — the audit-confirmed gate.  Returns
        False if any required metric is missing."""
        if self.sc is None or self.bsa is None or self.ddg is None:
            return False
        return (
            self.sc >= thresholds.orthogonal_sc_min
            and self.bsa >= thresholds.orthogonal_bsa_min
            and self.ddg <= thresholds.orthogonal_ddg_max
        )

    def failed_filter_names(self, thresholds) -> List[str]:
        """Names of the metrics that failed the gate.  Empty list iff
        ``passes_orthogonal_filters`` is True.  ``missing_X`` for any
        metric whose value is None."""
        failed: List[str] = []
        if self.sc is None:
            failed.append("missing_sc")
        elif self.sc < thresholds.orthogonal_sc_min:
            failed.append("sc")
        if self.bsa is None:
            failed.append("missing_bsa")
        elif self.bsa < thresholds.orthogonal_bsa_min:
            failed.append("bsa")
        if self.ddg is None:
            failed.append("missing_ddg")
        elif self.ddg > thresholds.orthogonal_ddg_max:
            failed.append("ddg")
        return failed

    # ── Informational (NOT in the gate) ──────────────────────────────

    def af3_disagrees(self, thresholds) -> bool:
        """True iff AF3 best ra_eff is at or above the disagreement
        threshold.  Diagnostic only; never part of the gate."""
        return self.af3.af3_disagrees(thresholds)

    # ── Output ─────────────────────────────────────────────────────

    def to_summary_row(self) -> dict:
        """Flat dict matching the existing
        ``survivors_with_orthogonal_metrics.csv`` schema for drop-in
        compatibility with downstream consumers."""
        row = {
            "mpnn_sequence": self.mpnn_sequence_id,
            "sc": "" if self.sc is None else f"{self.sc:.4f}",
            "bsa": "" if self.bsa is None else f"{self.bsa:.2f}",
            "hbonds": "" if self.hbonds is None else str(self.hbonds),
            "rosetta_ddg": "" if self.ddg is None else f"{self.ddg:.2f}",
        }
        row.update(self.af3.to_summary_row(self.mpnn_sequence_id))
        return row
