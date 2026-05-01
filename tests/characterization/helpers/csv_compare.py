from __future__ import annotations

import math
from collections.abc import Iterable, Sequence
from pathlib import Path

import pandas as pd

from tests.characterization.helpers.result import ComparisonResult

_MAX_DIFFS = 10
_STRATEGY_EXACT = "CSV-EXACT"
_STRATEGY_STRUCT = "CSV-STRUCT"


def compare_csv_exact(
    reference: Path,
    actual: Path,
    *,
    abs_tol: float = 1e-6,
    rel_tol: float = 1e-9,
    string_columns: Iterable[str] | None = None,
) -> ComparisonResult:
    """Strategy CSV-EXACT.

    Header (column names + order) must match exactly. Row count and row order
    must match. String columns are compared byte-equal; numeric columns within
    `abs(a - b) <= abs_tol + rel_tol * abs(b)`. NaN equals NaN; +/-Inf equal
    only to themselves. ``string_columns`` overrides the default dtype-based
    string/numeric classification (an explicit name list forces those columns
    to byte-equal comparison regardless of inferred dtype).

    On mismatch, ``differences`` carries up to ~10 row/column locations so
    failure output stays digestible.
    """
    early = _file_check(reference, actual, _STRATEGY_EXACT)
    if early is not None:
        return early
    ref_df = pd.read_csv(reference)
    act_df = pd.read_csv(actual)
    return _compare_dataframes(
        ref_df, act_df,
        reference=reference, actual=actual, strategy=_STRATEGY_EXACT,
        abs_tol=abs_tol, rel_tol=rel_tol, string_columns=string_columns,
    )


def compare_csv_struct(
    reference: Path,
    actual: Path,
    *,
    sort_by: Sequence[str],
    abs_tol: float = 1e-6,
    rel_tol: float = 1e-9,
    string_columns: Iterable[str] | None = None,
) -> ComparisonResult:
    """Strategy CSV-STRUCT.

    Same contract as :func:`compare_csv_exact` except both frames are sorted
    by ``sort_by`` before comparison. Use this when the producer does not
    guarantee row order. ``sort_by`` must reference existing columns in both
    frames; otherwise a clear failure is returned without attempting cell
    comparison.
    """
    early = _file_check(reference, actual, _STRATEGY_STRUCT)
    if early is not None:
        return early
    ref_df = pd.read_csv(reference)
    act_df = pd.read_csv(actual)

    sort_keys = list(sort_by)
    missing_ref = [c for c in sort_keys if c not in ref_df.columns]
    missing_act = [c for c in sort_keys if c not in act_df.columns]
    if missing_ref or missing_act:
        return ComparisonResult(
            passed=False,
            reference=reference,
            actual=actual,
            strategy=_STRATEGY_STRUCT,
            message="invalid sort_by: column(s) not present in both frames",
            differences=[
                f"missing in reference: {missing_ref}",
                f"missing in actual:    {missing_act}",
            ],
        )

    ref_sorted = ref_df.sort_values(by=sort_keys, kind="mergesort").reset_index(drop=True)
    act_sorted = act_df.sort_values(by=sort_keys, kind="mergesort").reset_index(drop=True)
    return _compare_dataframes(
        ref_sorted, act_sorted,
        reference=reference, actual=actual, strategy=_STRATEGY_STRUCT,
        abs_tol=abs_tol, rel_tol=rel_tol, string_columns=string_columns,
    )


def _file_check(reference: Path, actual: Path, strategy: str) -> ComparisonResult | None:
    if not Path(reference).exists():
        return ComparisonResult(
            passed=False, reference=reference, actual=actual, strategy=strategy,
            message=f"reference file not found: {reference}",
        )
    if not Path(actual).exists():
        return ComparisonResult(
            passed=False, reference=reference, actual=actual, strategy=strategy,
            message=f"actual file not found: {actual}",
        )
    return None


def _string_column_set(df: pd.DataFrame, override: Iterable[str] | None) -> set[str]:
    if override is not None:
        return set(override)
    return {c for c in df.columns if not pd.api.types.is_numeric_dtype(df[c])}


def _values_equal(a, b, *, abs_tol: float, rel_tol: float, treat_as_string: bool) -> bool:
    a_nan = isinstance(a, float) and math.isnan(a)
    b_nan = isinstance(b, float) and math.isnan(b)
    if a_nan and b_nan:
        return True
    if a_nan != b_nan:
        return False
    if treat_as_string:
        return str(a) == str(b)
    try:
        af = float(a)
        bf = float(b)
    except (TypeError, ValueError):
        return a == b
    if math.isinf(af) or math.isinf(bf):
        return af == bf
    return abs(af - bf) <= abs_tol + rel_tol * abs(bf)


def _compare_dataframes(
    ref_df: pd.DataFrame,
    act_df: pd.DataFrame,
    *,
    reference: Path,
    actual: Path,
    strategy: str,
    abs_tol: float,
    rel_tol: float,
    string_columns: Iterable[str] | None,
) -> ComparisonResult:
    ref_cols = list(ref_df.columns)
    act_cols = list(act_df.columns)
    if ref_cols != act_cols:
        return ComparisonResult(
            passed=False, reference=reference, actual=actual, strategy=strategy,
            message="header mismatch",
            differences=[
                f"reference columns: {ref_cols}",
                f"actual columns:    {act_cols}",
            ],
        )
    if len(ref_df) != len(act_df):
        return ComparisonResult(
            passed=False, reference=reference, actual=actual, strategy=strategy,
            message=f"row count mismatch: reference={len(ref_df)} actual={len(act_df)}",
        )

    string_cols = _string_column_set(ref_df, string_columns)
    diffs: list[str] = []
    n_rows = len(ref_df)
    for ri in range(n_rows):
        for col in ref_cols:
            rv = ref_df.iat[ri, ref_cols.index(col)]
            av = act_df.iat[ri, ref_cols.index(col)]
            if not _values_equal(
                rv, av, abs_tol=abs_tol, rel_tol=rel_tol,
                treat_as_string=col in string_cols,
            ):
                diffs.append(f"row {ri}, column {col!r}: reference={rv!r} actual={av!r}")
                if len(diffs) >= _MAX_DIFFS:
                    break
        if len(diffs) >= _MAX_DIFFS:
            break

    if diffs:
        return ComparisonResult(
            passed=False, reference=reference, actual=actual, strategy=strategy,
            message=f"{len(diffs)}{'+' if len(diffs) >= _MAX_DIFFS else ''} cell mismatch(es)",
            differences=diffs,
        )
    return ComparisonResult(passed=True, reference=reference, actual=actual, strategy=strategy)
