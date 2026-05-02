"""Stage 7 — Plots.

Pins each PNG in ``plots/`` with the default PNG-PERCEPTUAL strategy
(``ssim_threshold=0.95``, ``size_tolerance_pct=5.0``). Per-plot threshold
overrides are deferred — the first HPC round-trip will surface which (if
any) plots need them.

Producers (by filename prefix; full mapping in plan §2.7):
- ``rfdiff_*.png``     → ``RFDIFFUSION_PLOTS`` / ``bin/rfdiffusion_plots.py``
- ``rosetta_*.png``    → ``ROSETTA_FILTER_PLOTS`` / ``bin/rosetta_filter_plots.py``
- ``mpnn_*.png``       → ``MPNN_PLOTS`` / ``bin/mpnn_plots.py``
- ``negsteer_*.png``   → ``NEGSTEER_PLOTS`` / ``NEGSTEER_WITHIN_SEQUENCE_PLOTS``
                         (``bin/negsteer_plots.py`` / within-sequence variants)
- ``orthogonal_*.png`` → ``ORTHOG_PLOTS`` / ``bin/orthogonal_plots.py``
"""
from __future__ import annotations

import os
from pathlib import Path

import pytest

from tests.characterization.helpers.png_compare import compare_png_perceptual


def _safe_ref_root() -> Path | None:
    """Resolve reference root at import time so parametrize ids are stable."""
    override = os.environ.get("RECEPTOR_REF_ROOT")
    if override:
        candidate = Path(override)
    else:
        here = Path(__file__).resolve()
        candidate = None
        for parent in (here, *here.parents):
            if (parent / "nextflow.config").exists() or (parent / "main.nf").exists():
                candidate = parent / "tests" / "full_test_run" / "example_output_files"
                break
        if candidate is None:
            return None
    return candidate if candidate.exists() else None


def _list_plot_pngs(ref_root: Path | None) -> list[str]:
    if ref_root is None:
        return []
    plots = ref_root / "plots"
    if not plots.exists():
        return []
    return sorted(p.name for p in plots.iterdir() if p.is_file() and p.suffix == ".png")


_REF = _safe_ref_root()
_PNGS = _list_plot_pngs(_REF)


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("filename", _PNGS)
def test_plot_from_plot_module(reference_root, output_root, filename):
    """matplotlib plot module (see filename prefix in module docstring) — PNG-PERCEPTUAL."""
    rel = f"plots/{filename}"
    compare_png_perceptual(reference_root / rel, output_root / rel).assert_passed()
