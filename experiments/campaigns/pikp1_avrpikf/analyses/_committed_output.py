"""Refuse to replace a committed CSV with a smaller one.

Three extractors in this campaign read their inputs from
``runs/crystal_full_test_contig/results/``, which is gitignored and therefore
absent from any fresh checkout, and write a CSV that *is* committed. The glob
that finds the inputs returns an empty list when that directory is missing, and
nothing downstream objects, so the script exits 0 having overwritten a complete
table with a near-empty one. `design_secondary_structure.csv` is 5188 rows and
would come back as 1, the lone hardcoded wild-type entry.

The committed tables are the provenance of figures and numbers in the thesis,
so losing one is not recoverable by re-running. The input it was built from is
not in version control either.

So a write goes through `guard`, which compares the new row count against the
committed file and stops before touching it when the count has fallen.
`--allow-partial` is the deliberate override, for the case where a shorter
table is the intended result.

This is the pattern structure-negative-steering adopted in its PR #11, after
the same shape of bug silently rewrote a 452-row table with 324 rows.
"""

from __future__ import annotations

import os
import sys


class WouldShrink(SystemExit):
    """Raised instead of writing. SystemExit so an unguarded caller still stops."""


def committed_row_count(path: str) -> int | None:
    """Data rows in an existing CSV, or None when there is nothing to protect.

    Counts lines and subtracts the header rather than parsing, because the
    comparison only needs a magnitude and a parse would need the dialect.
    """
    if not os.path.exists(path):
        return None
    with open(path, newline="") as fh:
        lines = sum(1 for _ in fh)
    return max(0, lines - 1)


def guard(path: str, new_rows: int, *, allow_partial: bool = False,
          label: str | None = None) -> None:
    """Stop before replacing `path` with fewer rows than it already holds.

    Returns normally when the file is new, when the count held or grew, or when
    allow_partial is set. Raises WouldShrink otherwise, leaving the file as it
    was.
    """
    what = label or os.path.basename(path)
    have = committed_row_count(path)
    if have is None:
        print(f"  {what}: new file, {new_rows} rows")
        return
    if new_rows >= have:
        print(f"  {what}: {have} -> {new_rows} rows")
        return
    if allow_partial:
        print(f"  {what}: {have} -> {new_rows} rows "
              f"(SHRINKING, allowed by --allow-partial)", file=sys.stderr)
        return
    raise WouldShrink(
        f"\nREFUSING TO WRITE {path}\n"
        f"  committed: {have} data rows\n"
        f"  this run:  {new_rows} data rows\n\n"
        f"The committed table is longer than what this run produced, which is\n"
        f"what a missing input directory looks like. This campaign's inputs\n"
        f"live under runs/crystal_full_test_contig/results/, which is\n"
        f"gitignored, so a fresh checkout does not have them.\n\n"
        f"The file on disk is untouched. Pass --allow-partial if a shorter\n"
        f"table really is what you meant.\n")
