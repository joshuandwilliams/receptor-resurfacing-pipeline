#!/usr/bin/env python3
"""
test_negsteer_plots.py
----------------------
Iteration harness for negsteer cohort plots.  Thin wrapper around
``bin/negsteer_plots.py``.

Previously this file was a near-verbatim duplicate of the production
script (~1900 LOC).  Phase 4 architecture migration eliminated the
duplication; all plot logic now lives in `bin/negsteer_plots.py` and
this file is a 20-line shim that adds the iteration-time conveniences:

- ``--csv`` is OPTIONAL; falls back to the canonical per-module test
  output path when omitted.
- Same flag set, same outputs, same plotting code as production.

Same CLI as before — anything that worked before still works.
"""

from __future__ import annotations

import sys
from pathlib import Path


def _has_arg(name: str) -> bool:
    return any(a == name or a.startswith(name + "=") for a in sys.argv)


# Inject the canonical test-output path for --csv when not supplied.
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

# Pull in the production plot script (it lives in bin/ alongside the
# rest of the pipeline scripts).
_BIN_DIR = Path(__file__).resolve().parents[2] / "bin"
if str(_BIN_DIR) not in sys.path:
    sys.path.insert(0, str(_BIN_DIR))

import negsteer_plots  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(negsteer_plots.main())
