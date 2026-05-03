from __future__ import annotations

import json
import math
from collections.abc import Callable
from pathlib import Path

from tests.characterization.helpers.result import ComparisonResult

_MAX_DIFFS = 10
_STRATEGY_DEEP = "JSON-DEEP"
_STRATEGY_MODULO = "JSON-MODULO-PATHS"


def compare_json_deep(
    reference: Path,
    actual: Path,
    *,
    abs_tol: float = 1e-6,
    rel_tol: float = 1e-9,
    list_orders_significant: bool = True,
) -> ComparisonResult:
    """Strategy JSON-DEEP.

    Recursive structural comparison. Numeric scalars are compared with
    ``abs(a - b) <= abs_tol + rel_tol * abs(b)``; NaN equals NaN; +/-Inf
    equal only to themselves. Lists are compared positionally unless
    ``list_orders_significant=False``, in which case both lists are sorted
    by their JSON repr before pairwise comparison.

    Booleans are treated as a distinct type from integers — JSON ``true``
    does NOT match JSON ``1`` even though ``isinstance(True, int)`` is
    True in Python.

    Mismatches report nested locations like ``designs[3].motif_rmsd``;
    the ``differences`` list is capped at ~10 entries.
    """
    early = _file_check(reference, actual, _STRATEGY_DEEP)
    if early is not None:
        return early
    ref_data = json.loads(Path(reference).read_text())
    act_data = json.loads(Path(actual).read_text())
    diffs = _diff(ref_data, act_data, "", abs_tol, rel_tol, list_orders_significant)
    return _build_result(diffs, reference, actual, _STRATEGY_DEEP)


def compare_json_modulo_paths(
    reference: Path,
    actual: Path,
    path_normalizer: Callable[[str], str],
    *,
    abs_tol: float = 1e-6,
    rel_tol: float = 1e-9,
    list_orders_significant: bool = True,
) -> ComparisonResult:
    """Strategy JSON-MODULO-PATHS.

    Same as :func:`compare_json_deep` but every string *value* is passed
    through ``path_normalizer`` before comparison. Object keys are not
    normalised. Use this for files that contain absolute Nextflow workdir
    paths or user-home prefixes that vary between runs.
    """
    early = _file_check(reference, actual, _STRATEGY_MODULO)
    if early is not None:
        return early
    ref_data = json.loads(Path(reference).read_text())
    act_data = json.loads(Path(actual).read_text())
    ref_n = _normalize_strings(ref_data, path_normalizer)
    act_n = _normalize_strings(act_data, path_normalizer)
    diffs = _diff(ref_n, act_n, "", abs_tol, rel_tol, list_orders_significant)
    return _build_result(diffs, reference, actual, _STRATEGY_MODULO)


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


def _normalize_strings(data, fn: Callable[[str], str]):
    if isinstance(data, dict):
        return {k: _normalize_strings(v, fn) for k, v in data.items()}
    if isinstance(data, list):
        return [_normalize_strings(v, fn) for v in data]
    if isinstance(data, str):
        return fn(data)
    return data


def _path(parent: str, key) -> str:
    if parent:
        return f"{parent}.{key}"
    return str(key)


def _here(parent: str) -> str:
    return parent if parent else "<root>"


def _diff(ref, act, path: str, abs_tol: float, rel_tol: float, list_orders_significant: bool) -> list[str]:
    ref_is_bool = isinstance(ref, bool)
    act_is_bool = isinstance(act, bool)
    if ref_is_bool != act_is_bool:
        return [f"{_here(path)}: type mismatch (reference={type(ref).__name__}, actual={type(act).__name__})"]
    if ref_is_bool and act_is_bool:
        if ref != act:
            return [f"{_here(path)}: reference={ref!r} actual={act!r}"]
        return []

    if ref is None or act is None:
        if ref is act:
            return []
        return [f"{_here(path)}: reference={ref!r} actual={act!r}"]

    if isinstance(ref, dict) and isinstance(act, dict):
        return _diff_dicts(ref, act, path, abs_tol, rel_tol, list_orders_significant)

    if isinstance(ref, list) and isinstance(act, list):
        return _diff_lists(ref, act, path, abs_tol, rel_tol, list_orders_significant)

    if isinstance(ref, (int, float)) and isinstance(act, (int, float)):
        return _diff_numbers(ref, act, path, abs_tol, rel_tol)

    if isinstance(ref, str) and isinstance(act, str):
        if ref != act:
            return [f"{_here(path)}: reference={ref!r} actual={act!r}"]
        return []

    return [f"{_here(path)}: type mismatch (reference={type(ref).__name__}, actual={type(act).__name__})"]


def _diff_dicts(ref: dict, act: dict, path: str, abs_tol: float, rel_tol: float, list_orders_significant: bool) -> list[str]:
    out: list[str] = []
    ref_keys = set(ref.keys())
    act_keys = set(act.keys())
    for k in sorted(ref_keys - act_keys):
        out.append(f"{_path(path, k)}: missing in actual")
        if len(out) >= _MAX_DIFFS:
            return out
    for k in sorted(act_keys - ref_keys):
        out.append(f"{_path(path, k)}: unexpected in actual")
        if len(out) >= _MAX_DIFFS:
            return out
    for k in sorted(ref_keys & act_keys):
        out.extend(_diff(ref[k], act[k], _path(path, k), abs_tol, rel_tol, list_orders_significant))
        if len(out) >= _MAX_DIFFS:
            return out[:_MAX_DIFFS]
    return out


def _diff_lists(ref: list, act: list, path: str, abs_tol: float, rel_tol: float, list_orders_significant: bool) -> list[str]:
    if len(ref) != len(act):
        return [f"{_here(path)}: list length mismatch (reference={len(ref)} actual={len(act)})"]
    if list_orders_significant:
        ref_iter, act_iter = ref, act
    else:
        key = lambda x: json.dumps(x, sort_keys=True, default=str)
        ref_iter = sorted(ref, key=key)
        act_iter = sorted(act, key=key)
    out: list[str] = []
    for i, (r, a) in enumerate(zip(ref_iter, act_iter)):
        out.extend(_diff(r, a, f"{path}[{i}]", abs_tol, rel_tol, list_orders_significant))
        if len(out) >= _MAX_DIFFS:
            return out[:_MAX_DIFFS]
    return out


def _diff_numbers(ref, act, path: str, abs_tol: float, rel_tol: float) -> list[str]:
    rf = float(ref)
    af = float(act)
    if math.isnan(rf) and math.isnan(af):
        return []
    if math.isnan(rf) or math.isnan(af):
        return [f"{_here(path)}: reference={ref!r} actual={act!r}"]
    if math.isinf(rf) or math.isinf(af):
        if rf == af:
            return []
        return [f"{_here(path)}: reference={ref!r} actual={act!r}"]
    if abs(rf - af) <= abs_tol + rel_tol * abs(af):
        return []
    return [f"{_here(path)}: reference={ref!r} actual={act!r}"]


def _build_result(diffs: list[str], reference: Path, actual: Path, strategy: str) -> ComparisonResult:
    if not diffs:
        return ComparisonResult(passed=True, reference=reference, actual=actual, strategy=strategy)
    capped = diffs[:_MAX_DIFFS]
    suffix = "+" if len(diffs) > _MAX_DIFFS else ""
    return ComparisonResult(
        passed=False, reference=reference, actual=actual, strategy=strategy,
        message=f"{len(capped)}{suffix} difference(s)",
        differences=capped,
    )
