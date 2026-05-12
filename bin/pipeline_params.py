"""
pipeline_params.py
------------------
PipelineParams — runtime wrapper around the user-facing parameter set
validated by bin/validate_params.py PARAM_SPECS.

The validator file (validate_params.py) is the catalogue of WHAT to
check.  This module is the runtime carrier of THE VALUES.  Construction
runs full validation, collects every error, raises ParamValidationError
if any field is bad.  Successfully constructed → all values are valid.

Type design (Phase 4 spec §2.12):
- Immutable after construction (frozen dataclass-like behaviour).
- Construction via from_dict / from_nextflow_json; both routes funnel
  through validate_params.validate_params for spec enforcement.
- Carries USER-facing params only.  Research thresholds (ra_eff
  cutoffs, intact thresholds, composite weights) live in
  bin/pipeline_thresholds.py.

Tier 0 type — no upstream dependencies on other Phase 4 types.
"""

from __future__ import annotations

import json
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Mapping

_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

import validate_params  # noqa: E402  (bin/validate_params.py)


class ParamValidationError(ValueError):
    """Raised when one or more parameters fail validation at construction.

    Carries the full list of error messages (not just the first) so
    callers can surface every problem in one go.
    """

    def __init__(self, errors: List[str]):
        self.errors = errors
        super().__init__(self._format(errors))

    @staticmethod
    def _format(errors: List[str]) -> str:
        joined = "\n".join(errors)
        return (
            f"Parameter validation failed ({len(errors)} error(s)):\n{joined}"
        )


@dataclass(frozen=True)
class PipelineParams:
    """The validated, user-facing parameter set for one pipeline run.

    Constructed via classmethods (``from_dict`` / ``from_nextflow_json``)
    rather than directly so that validation always runs.  Values are
    stored in a single dict (``values``) keyed by param name; access
    via attribute-like ``params["key"]`` or via ``params.values["key"]``.

    A dict rather than one field per param keeps the type small and
    avoids drift between two catalogues — `validate_params.PARAM_SPECS`
    is the single source of truth for which params exist.  Accessor
    methods on this type return the same dict view.
    """

    values: Mapping[str, Any] = field(default_factory=dict)

    def __getitem__(self, key: str) -> Any:
        return self.values[key]

    def __contains__(self, key: str) -> bool:
        return key in self.values

    def get(self, key: str, default: Any = None) -> Any:
        return self.values.get(key, default)

    def as_dict(self) -> Dict[str, Any]:
        """Flat snapshot of all stored values."""
        return dict(self.values)

    def validate(self) -> List[str]:
        """Re-run every spec against the current values.

        Useful as a defence-in-depth check (e.g. before passing the
        params to a downstream stage that might silently tolerate a bad
        value).  Returns empty list on success.
        """
        return validate_params.validate_params(dict(self.values))

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "PipelineParams":
        """Construct from a flat dict (e.g. the JSON Nextflow dumps).

        Runs validation; raises ParamValidationError if any field is
        invalid (with all errors collected, not just the first).
        """
        errors = validate_params.validate_params(dict(d))
        if errors:
            raise ParamValidationError(errors)
        return cls(values=dict(d))

    @classmethod
    def from_nextflow_json(cls, path: Path) -> "PipelineParams":
        """Construct from a JSON file dumped by Nextflow's
        groovy.json.JsonOutput.toJson(params).  Matches the wiring
        already in main.nf workflow head."""
        path = Path(path)
        if not path.is_file():
            raise FileNotFoundError(f"params JSON not found: {path}")
        try:
            raw = json.loads(path.read_text())
        except json.JSONDecodeError as e:
            raise ValueError(f"malformed params JSON ({path}): {e}") from e
        if not isinstance(raw, dict):
            raise TypeError(
                f"params JSON must be an object, got {type(raw).__name__}"
            )
        return cls.from_dict(raw)
