from __future__ import annotations

from pathlib import Path

import numpy as np
from PIL import Image, UnidentifiedImageError
from skimage.metrics import structural_similarity

from tests.characterization.helpers.result import ComparisonResult

_STRATEGY_PERCEPTUAL = "PNG-PERCEPTUAL"
_STRATEGY_EXISTS = "PNG-EXISTS"


def compare_png_perceptual(
    reference: Path,
    actual: Path,
    *,
    ssim_threshold: float = 0.95,
    size_tolerance_pct: float = 5.0,
) -> ComparisonResult:
    """Strategy PNG-PERCEPTUAL.

    Three-stage check:

    1. Both files exist, are non-empty, and decode as valid PNGs.
    2. File sizes are within ``size_tolerance_pct`` of each other, computed
       as ``abs(ref - act) / max(ref, act) * 100``.
    3. Structural similarity (SSIM) on the RGB pixel arrays is at least
       ``ssim_threshold``. If actual's dimensions differ from reference's,
       actual is resized via PIL bilinear interpolation before SSIM.

    The 5% default for ``size_tolerance_pct`` is calibrated for real
    matplotlib PNGs (typically 50-500 KB). On very small procedurally
    generated test images, percentage size differences inflate because
    PNG compression is sensitive to small perturbations in compact images;
    this is a property of small images, not a defect. Tests that exercise
    the SSIM path on small synthetic images should pass an explicit
    relaxed ``size_tolerance_pct``.
    """
    early = _file_check(reference, actual, _STRATEGY_PERCEPTUAL)
    if early is not None:
        return early
    ref_size = Path(reference).stat().st_size
    act_size = Path(actual).stat().st_size
    if ref_size == 0:
        return ComparisonResult(
            passed=False, reference=reference, actual=actual, strategy=_STRATEGY_PERCEPTUAL,
            message="reference file is empty",
        )
    if act_size == 0:
        return ComparisonResult(
            passed=False, reference=reference, actual=actual, strategy=_STRATEGY_PERCEPTUAL,
            message="actual file is empty",
        )

    ref_img = _open_png(reference)
    if isinstance(ref_img, ComparisonResult):
        return ref_img
    act_img = _open_png(actual)
    if isinstance(act_img, ComparisonResult):
        return act_img

    larger = max(ref_size, act_size)
    pct = (abs(ref_size - act_size) / larger) * 100.0
    if pct > size_tolerance_pct:
        return ComparisonResult(
            passed=False, reference=reference, actual=actual, strategy=_STRATEGY_PERCEPTUAL,
            message=f"file size differs by {pct:.2f}% (tolerance {size_tolerance_pct}%)",
            differences=[f"reference size: {ref_size} bytes", f"actual size:    {act_size} bytes"],
        )

    if act_img.size != ref_img.size:
        act_img = act_img.resize(ref_img.size, Image.Resampling.BILINEAR)

    ref_arr = np.asarray(ref_img)
    act_arr = np.asarray(act_img)
    ssim_val = float(
        structural_similarity(
            ref_arr, act_arr,
            channel_axis=-1, data_range=255,
        )
    )
    if ssim_val < ssim_threshold:
        return ComparisonResult(
            passed=False, reference=reference, actual=actual, strategy=_STRATEGY_PERCEPTUAL,
            message=f"SSIM={ssim_val:.4f} below threshold {ssim_threshold}",
        )
    return ComparisonResult(
        passed=True, reference=reference, actual=actual, strategy=_STRATEGY_PERCEPTUAL,
        message=f"SSIM={ssim_val:.4f}",
    )


def compare_png_exists(reference: Path, actual: Path) -> ComparisonResult:
    """Strategy PNG-EXISTS — escape hatch for plots too stochastic for SSIM.

    Both files must exist and be non-empty. No pixel comparison.
    """
    early = _file_check(reference, actual, _STRATEGY_EXISTS)
    if early is not None:
        return early
    if Path(reference).stat().st_size == 0:
        return ComparisonResult(
            passed=False, reference=reference, actual=actual, strategy=_STRATEGY_EXISTS,
            message="reference file is empty",
        )
    if Path(actual).stat().st_size == 0:
        return ComparisonResult(
            passed=False, reference=reference, actual=actual, strategy=_STRATEGY_EXISTS,
            message="actual file is empty",
        )
    return ComparisonResult(
        passed=True, reference=reference, actual=actual, strategy=_STRATEGY_EXISTS,
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


def _open_png(path: Path) -> Image.Image | ComparisonResult:
    try:
        with Image.open(path) as probe:
            probe.verify()
        return Image.open(path).convert("RGB")
    except (UnidentifiedImageError, OSError) as exc:
        return ComparisonResult(
            passed=False, reference=path, actual=path, strategy=_STRATEGY_PERCEPTUAL,
            message=f"not a valid PNG ({path.name}): {exc}",
        )
