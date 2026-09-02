#!/usr/bin/env python3
"""Whether the pipeline steered each unit, read from its own plan.

`boltz2_negative_steering.py` skips steering only when every cold-start seed
passes, at ra_eff <= rmsd_threshold with receptor and effector intact, and
records that decision as `cold_start_all_clean` in cycle_0/plan.json. Any
label derived here from a median or from a single cold ra_eff disagrees with
what the run actually did, so the plan is the only source used.
"""

import glob
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, "..", "..", "runs", "crystal_full_test_contig",
                    "results", "negative_steering", "runs")


def never_steered():
    """Units the pipeline left unsteered because the cold start was all clean."""
    out = set()
    for path in glob.glob(os.path.join(RUNS, "*", "cycle_0", "plan.json")):
        unit = os.path.basename(os.path.dirname(os.path.dirname(path)))
        with open(path) as fh:
            if json.load(fh).get("cold_start_all_clean"):
                out.add(unit)
    if not out:
        raise SystemExit(f"no plan.json found under {RUNS}; sync results/ first")
    return out
