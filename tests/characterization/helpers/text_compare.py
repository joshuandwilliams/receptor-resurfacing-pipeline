from __future__ import annotations

from pathlib import Path

from tests.characterization.helpers.result import ComparisonResult

_STRATEGY_TEXT = "TEXT-EXACT"


def compare_text_exact(reference: Path, actual: Path) -> ComparisonResult:
    """Strategy TEXT-EXACT.

    Byte-equal comparison of two text files. Used for outputs whose contents
    are line-deterministic but not strictly CSV/JSON (e.g. ``processed_contigs.txt``,
    ``passing_designs.txt``, FASTA files, ``summary.txt``).
    """
    early = _file_check(reference, actual)
    if early is not None:
        return early
    ref_bytes = Path(reference).read_bytes()
    act_bytes = Path(actual).read_bytes()
    if ref_bytes == act_bytes:
        return ComparisonResult(
            passed=True, reference=reference, actual=actual, strategy=_STRATEGY_TEXT,
        )
    return ComparisonResult(
        passed=False, reference=reference, actual=actual, strategy=_STRATEGY_TEXT,
        message=f"byte mismatch (reference={len(ref_bytes)} bytes, actual={len(act_bytes)} bytes)",
        differences=_first_line_diffs(ref_bytes, act_bytes),
    )


def _file_check(reference: Path, actual: Path) -> ComparisonResult | None:
    if not Path(reference).exists():
        return ComparisonResult(
            passed=False, reference=reference, actual=actual, strategy=_STRATEGY_TEXT,
            message=f"reference file not found: {reference}",
        )
    if not Path(actual).exists():
        return ComparisonResult(
            passed=False, reference=reference, actual=actual, strategy=_STRATEGY_TEXT,
            message=f"actual file not found: {actual}",
        )
    return None


def _first_line_diffs(ref_bytes: bytes, act_bytes: bytes, max_diffs: int = 10) -> list[str]:
    ref_lines = ref_bytes.splitlines()
    act_lines = act_bytes.splitlines()
    diffs: list[str] = []
    for i in range(max(len(ref_lines), len(act_lines))):
        rl = ref_lines[i] if i < len(ref_lines) else None
        al = act_lines[i] if i < len(act_lines) else None
        if rl != al:
            diffs.append(f"line {i + 1}: reference={rl!r} actual={al!r}")
            if len(diffs) >= max_diffs:
                break
    return diffs
