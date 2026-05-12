#!/usr/bin/env python3
"""
test_negsteer_within_sequence_plots.py
--------------------------------------
Iteration harness for within-sequence (per-seed) negsteer plots.  Thin
wrapper around ``bin/negsteer_within_sequence_plots.py``.

Previously a verbatim duplicate of the production script (~980 LOC).
Phase 4 architecture migration eliminated the duplication; all plot
logic lives in `bin/negsteer_within_sequence_plots.py` and this file
is a thin shim.
"""

from __future__ import annotations

import sys
from pathlib import Path


def _has_arg(name: str) -> bool:
    return any(a == name or a.startswith(name + "=") for a in sys.argv)


# Canonical test-output paths for the two main inputs.
if not _has_arg("--csv"):
    canonical = Path(
        "tests/negative_steering/receptor_resurfacing_results/"
        "negative_steering/cross_sequence_summary.csv"
    )
    if not canonical.exists():
        raise SystemExit(
            f"--csv not supplied and canonical test path does not exist:\n"
            f"  {canonical}\n"
            f"Run test_negative_steering first or pass --csv explicitly."
        )
    sys.argv.extend(["--csv", str(canonical)])

_BIN_DIR = Path(__file__).resolve().parents[2] / "bin"
if str(_BIN_DIR) not in sys.path:
    sys.path.insert(0, str(_BIN_DIR))

import negsteer_within_sequence_plots  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(negsteer_within_sequence_plots.main())
