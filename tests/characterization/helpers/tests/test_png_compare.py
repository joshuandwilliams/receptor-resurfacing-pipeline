from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from PIL import Image, ImageDraw

from tests.characterization.helpers.png_compare import (
    compare_png_exists,
    compare_png_perceptual,
)


def _plot_like_image(size=(128, 96), seed: int = 0) -> Image.Image:
    rng = np.random.default_rng(seed)
    img = Image.new("RGB", size, "white")
    draw = ImageDraw.Draw(img)
    w, h = size
    draw.line((10, h - 10, w - 10, h - 10), fill="black", width=1)
    draw.line((10, 10, 10, h - 10), fill="black", width=1)
    for _ in range(12):
        x = int(rng.integers(15, w - 12))
        y = int(rng.integers(12, h - 15))
        draw.ellipse((x - 2, y - 2, x + 2, y + 2), fill="black")
    return img


def _save(img: Image.Image, path: Path) -> Path:
    img.save(path, format="PNG")
    return path


@pytest.mark.local_unit
def test_compare_png_perceptual_identity_passes(tmp_path):
    """Two byte-identical PNGs pass with high SSIM."""
    img = _plot_like_image()
    ref = _save(img, tmp_path / "r.png")
    act = _save(img, tmp_path / "a.png")
    assert compare_png_perceptual(ref, act).passed


@pytest.mark.local_unit
def test_compare_png_perceptual_minor_noise_passes_with_relaxed_size_tolerance(tmp_path):
    """Slight pixel perturbation on a plot-like image passes SSIM with relaxed size tolerance.

    Procedurally generated PNGs are small enough that PNG-compression-driven
    size differences inflate beyond the 5% default — see the docstring on
    compare_png_perceptual for the calibration note. The SSIM check on the
    decoded RGB arrays is what we are exercising here.
    """
    base = _plot_like_image(seed=0)
    ref = _save(base, tmp_path / "r.png")
    arr = np.asarray(base).copy()
    rng = np.random.default_rng(1)
    for _ in range(5):
        x = int(rng.integers(0, arr.shape[1]))
        y = int(rng.integers(0, arr.shape[0]))
        arr[y, x] = (240, 240, 240)
    perturbed = Image.fromarray(arr)
    act = _save(perturbed, tmp_path / "a.png")
    result = compare_png_perceptual(ref, act, size_tolerance_pct=50.0)
    assert result.passed


@pytest.mark.local_unit
def test_compare_png_perceptual_completely_different_images_fail(tmp_path):
    """A solid-black image vs a plot-like image fails SSIM (size tolerance disabled)."""
    ref = _save(_plot_like_image(seed=0), tmp_path / "r.png")
    act = _save(Image.new("RGB", (128, 96), "black"), tmp_path / "a.png")
    result = compare_png_perceptual(ref, act, size_tolerance_pct=99.0)
    assert not result.passed


@pytest.mark.local_unit
def test_compare_png_perceptual_size_outside_tolerance_fails(tmp_path):
    """File size difference exceeding size_tolerance_pct fails before SSIM is computed."""
    small = _plot_like_image(size=(50, 50), seed=0)
    big = _plot_like_image(size=(500, 500), seed=0)
    ref = _save(small, tmp_path / "r.png")
    act = _save(big, tmp_path / "a.png")
    result = compare_png_perceptual(ref, act, size_tolerance_pct=5.0)
    assert not result.passed
    assert "size" in result.message.lower()


@pytest.mark.local_unit
def test_compare_png_perceptual_ssim_threshold_override(tmp_path):
    """A high ssim_threshold causes a perceptibly-different image to fail."""
    base = _plot_like_image(seed=0)
    ref = _save(base, tmp_path / "r.png")
    arr = np.asarray(base).copy()
    arr[20:60, 20:80] = (200, 0, 0)
    act = _save(Image.fromarray(arr), tmp_path / "a.png")
    result = compare_png_perceptual(ref, act, ssim_threshold=0.999, size_tolerance_pct=99.0)
    assert not result.passed


@pytest.mark.local_unit
def test_compare_png_perceptual_missing_reference_fails(tmp_path):
    """Nonexistent reference file fails cleanly."""
    act = _save(_plot_like_image(), tmp_path / "a.png")
    result = compare_png_perceptual(tmp_path / "missing.png", act)
    assert not result.passed
    assert "reference" in result.message.lower()


@pytest.mark.local_unit
def test_compare_png_perceptual_empty_actual_fails(tmp_path):
    """Empty (zero-byte) actual file fails."""
    ref = _save(_plot_like_image(), tmp_path / "r.png")
    act = tmp_path / "a.png"
    act.write_bytes(b"")
    result = compare_png_perceptual(ref, act)
    assert not result.passed


@pytest.mark.local_unit
def test_compare_png_exists_passes_when_present(tmp_path):
    """compare_png_exists returns passed=True when both files exist and are non-empty."""
    ref = _save(_plot_like_image(), tmp_path / "r.png")
    act = _save(_plot_like_image(), tmp_path / "a.png")
    assert compare_png_exists(ref, act).passed


@pytest.mark.local_unit
def test_compare_png_exists_fails_when_empty(tmp_path):
    """compare_png_exists fails when actual is zero bytes."""
    ref = _save(_plot_like_image(), tmp_path / "r.png")
    act = tmp_path / "a.png"
    act.write_bytes(b"")
    assert not compare_png_exists(ref, act).passed


@pytest.mark.local_unit
def test_compare_png_exists_fails_when_missing(tmp_path):
    """compare_png_exists fails when actual file does not exist."""
    ref = _save(_plot_like_image(), tmp_path / "r.png")
    assert not compare_png_exists(ref, tmp_path / "missing.png").passed
