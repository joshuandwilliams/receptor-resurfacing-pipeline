#!/usr/bin/env python3
"""
test_orthogonal_metrics_plots.py
--------------------------------
Iteration harness for orthogonal-metrics cohort plots.  Thin wrapper
around ``bin/orthogonal_metrics_plots.py``.

Previously a near-verbatim duplicate of the production script
(~1140 LOC).  Phase 4 architecture migration eliminated the
duplication; all plot logic lives in `bin/orthogonal_metrics_plots.py`
and this file is a thin shim that adds the canonical-test-path
fallbacks for ``--survivors-csv`` and ``--cross-summary-csv``.
"""

from __future__ import annotations

import sys
from pathlib import Path


def _has_arg(name: str) -> bool:
    return any(a == name or a.startswith(name + "=") for a in sys.argv)


if not _has_arg("--survivors-csv"):
    canonical = Path(
        "tests/orthogonal_metrics/receptor_resurfacing_results/"
        "orthogonal_metrics/survivors_with_orthogonal_metrics.csv"
    )
    if not canonical.exists():
        raise SystemExit(
            f"--survivors-csv not supplied and canonical test path "
            f"does not exist:\n  {canonical}\n"
            f"Run test_orthogonal_metrics first or pass --survivors-csv "
            f"explicitly."
        )
    sys.argv.extend(["--survivors-csv", str(canonical)])

# cross-summary-csv has its own canonical fallback (from the upstream
# negsteer test); only inject if user didn't pass and it exists.
if not _has_arg("--cross-summary-csv"):
    canonical = Path(
        "tests/negative_steering/receptor_resurfacing_results/"
        "negative_steering/cross_sequence_summary.csv"
    )
    if canonical.exists():
        sys.argv.extend(["--cross-summary-csv", str(canonical)])

_BIN_DIR = Path(__file__).resolve().parents[2] / "bin"
if str(_BIN_DIR) not in sys.path:
    sys.path.insert(0, str(_BIN_DIR))

import orthogonal_metrics_plots  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(orthogonal_metrics_plots.main())
