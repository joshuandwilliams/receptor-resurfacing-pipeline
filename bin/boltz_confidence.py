"""
boltz_confidence.py
-------------------
BoltzConfidenceMetrics — wrapper around one Boltz prediction sample's
confidence sidecar JSON (and optional pLDDT npz).

Pairs 1:1 with one Boltz `ProteinStructurePrediction`.  Each sample
produces one PDB and one confidence JSON.  This type wraps the
non-structural side.

Scope is the JSON-direct metrics plus mean pLDDT.  PAE-derived metrics
(`ipae`, `pae_pass_frac`, `ipsae*`) are absent, because they need the
chain lengths from the matching PDB.  They are computed instead by the
negative-steering engine and by `compute_interface_metrics.py`, both of
which read the `pae_*.npz` sidecar directly.

Imports nothing from the other typed modules in `bin/`, only json,
pathlib and typing.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Optional


class BoltzConfidenceMetrics:
    """One sample's confidence metrics from a Boltz prediction.

    Construction takes the path to a `confidence_<name>.json` file and
    optionally a matching `plddt_<name>.npz`.  A missing confidence JSON
    raises FileNotFoundError.  The JSON is read lazily on first access,
    then cached for the instance's lifetime.

    Scalar properties return `float | None`.  `None` means the JSON was
    unparseable or the field was absent.  Otherwise the value is rounded
    to 4 dp.  `avg_plddt` is the exception and rounds to 2.
    """

    def __init__(
        self,
        confidence_json: Path,
        plddt_npz: Optional[Path] = None,
    ):
        self._confidence_json = Path(confidence_json)
        self._plddt_npz = Path(plddt_npz) if plddt_npz is not None else None
        if not self._confidence_json.is_file():
            raise FileNotFoundError(
                f"confidence JSON not found: {self._confidence_json}"
            )
        self._cached_json: Optional[Dict[str, Any]] = None
        self._cached_avg_plddt: Optional[float] = None
        self._avg_plddt_computed = False

    # ── Properties (lazy from JSON) ──────────────────────────────────

    def _load_json(self) -> Dict[str, Any]:
        if self._cached_json is None:
            try:
                with open(self._confidence_json) as f:
                    self._cached_json = json.load(f)
            except (json.JSONDecodeError, OSError):
                self._cached_json = {}
        return self._cached_json

    @property
    def confidence_score(self) -> Optional[float]:
        return self._field_float("confidence_score")

    @property
    def ptm(self) -> Optional[float]:
        return self._field_float("ptm")

    @property
    def iptm(self) -> Optional[float]:
        """Prefers the `iptm` field Boltz writes directly.  Falls back to
        `protein_iptm` when only the protein-protein variant is
        present."""
        val = self._field_float("iptm")
        if val is not None:
            return val
        return self._field_float("protein_iptm")

    @property
    def protein_iptm(self) -> Optional[float]:
        return self._field_float("protein_iptm")

    @property
    def ligand_iptm(self) -> Optional[float]:
        return self._field_float("ligand_iptm")

    @property
    def complex_plddt(self) -> Optional[float]:
        """Whole-complex pLDDT on a 0-1 scale (Boltz convention)."""
        return self._field_float("complex_plddt")

    @property
    def complex_iplddt(self) -> Optional[float]:
        return self._field_float("complex_iplddt")

    @property
    def complex_pde(self) -> Optional[float]:
        return self._field_float("complex_pde")

    @property
    def complex_ipde(self) -> Optional[float]:
        return self._field_float("complex_ipde")

    @property
    def chains_ptm(self) -> Dict[str, float]:
        """Per-chain pTM dict keyed by chain index (string)."""
        raw = self._load_json().get("chains_ptm", {})
        out: Dict[str, float] = {}
        if isinstance(raw, dict):
            for k, v in raw.items():
                try:
                    out[str(k)] = float(v)
                except (TypeError, ValueError):
                    continue
        return out

    @property
    def pair_chains_iptm(self) -> Dict[str, Dict[str, float]]:
        """Pairwise inter-chain ipTM dict, keyed by chain index strings."""
        raw = self._load_json().get("pair_chains_iptm", {})
        out: Dict[str, Dict[str, float]] = {}
        if isinstance(raw, dict):
            for k, inner in raw.items():
                if not isinstance(inner, dict):
                    continue
                inner_out: Dict[str, float] = {}
                for k2, v in inner.items():
                    try:
                        inner_out[str(k2)] = float(v)
                    except (TypeError, ValueError):
                        continue
                out[str(k)] = inner_out
        return out

    # ── Properties (from PLDDT npz, if provided) ─────────────────────

    @property
    def avg_plddt(self) -> Optional[float]:
        """Mean pLDDT across all residues, rescaled to 0-100.

        Boltz writes the per-residue .npz on a 0-1 scale.  Values are
        multiplied by 100 when the array max is <= 1.5, matching the
        0-100 PDB B-factor convention used elsewhere in the pipeline.
        `complex_plddt` by contrast stays on Boltz's 0-1 scale.
        Returns None when no .npz was supplied."""
        if self._avg_plddt_computed:
            return self._cached_avg_plddt
        self._avg_plddt_computed = True
        if self._plddt_npz is None or not self._plddt_npz.is_file():
            self._cached_avg_plddt = None
            return None
        try:
            import numpy as np  # local — keep top-level import light
            data = np.load(self._plddt_npz)
            arr = data[list(data.keys())[0]]
            if arr.ndim > 1:
                arr = arr.flatten()
            mean = float(arr.mean())
            if float(arr.max()) <= 1.5:
                mean *= 100.0
            self._cached_avg_plddt = round(mean, 2)
        except (OSError, ValueError, KeyError, ImportError):
            self._cached_avg_plddt = None
        return self._cached_avg_plddt

    # ── Derived ──────────────────────────────────────────────────────

    def confidence_flag(self, thresholds) -> str:
        """Returns one of: ok | low_iptm | low_plddt | multiple.

        Checks iptm and complex_plddt only, against
        `thresholds.iptm_min` and `thresholds.complex_plddt_min`.  The
        low_pass_frac and high_ipae triggers that
        `extract_passing._compute_confidence_flag` also applies need the
        PAE matrix, which this type does not read.
        """
        triggers: List[str] = []
        if self.iptm is not None and self.iptm < thresholds.iptm_min:
            triggers.append("low_iptm")
        if self.complex_plddt is not None and self.complex_plddt < thresholds.complex_plddt_min:
            triggers.append("low_plddt")
        if not triggers:
            return "ok"
        if len(triggers) == 1:
            return triggers[0]
        return "multiple"

    # ── Helpers ──────────────────────────────────────────────────────

    def _field_float(self, key: str) -> Optional[float]:
        raw = self._load_json().get(key)
        if raw is None:
            return None
        try:
            return round(float(raw), 4)
        except (TypeError, ValueError):
            return None

    @classmethod
    def from_prediction_dir(
        cls, prediction_dir: Path, base_name: str = "input"
    ) -> List["BoltzConfidenceMetrics"]:
        """Return one BoltzConfidenceMetrics per per-sample confidence
        JSON in `prediction_dir`, sorted by filename.

        Boltz writes one `confidence_<base_name>_model_<sample>.json`
        per sample, where sample runs 0..N-1 for N diffusion samples.
        The sort is lexicographic, so model_10 precedes model_2 once
        there are ten or more samples.  `base_name="input"` matches the
        Boltz output layout this pipeline uses.

        Each instance is linked to the matching `plddt_*.npz` when one
        is present in the same directory.
        """
        prediction_dir = Path(prediction_dir)
        if not prediction_dir.is_dir():
            raise FileNotFoundError(
                f"prediction dir not found: {prediction_dir}"
            )
        confs = sorted(
            prediction_dir.glob(f"confidence_{base_name}_model_*.json")
        )
        out: List["BoltzConfidenceMetrics"] = []
        for cj in confs:
            # Match plddt npz by replacing the prefix.
            stem = cj.name[len("confidence_"):-len(".json")]
            plddt = prediction_dir / f"plddt_{stem}.npz"
            out.append(cls(cj, plddt if plddt.is_file() else None))
        return out
