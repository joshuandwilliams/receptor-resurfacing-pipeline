"""Stage 7 — Plots.

Pins each PNG in the per-stage ``plots/`` directories with the default
PNG-PERCEPTUAL strategy (``ssim_threshold=0.95``,
``size_tolerance_pct=5.0``). Per-plot threshold overrides are deferred —
the first HPC round-trip will surface which (if any) plots need them.

Each parametrised case carries an explicit ``stage`` field so the
reference and output paths resolve under the per-module roots:

    rfdiffusion        → tests/rfdiffusion/example_output_files/plots/
    rosetta_filtering  → tests/rosetta_filtering/example_output_files/plots/
    proteinmpnn        → tests/proteinmpnn/example_output_files/plots/
    negative_steering  → tests/negative_steering/example_output_files/plots/
    orthogonal_metrics → tests/orthogonal_metrics/example_output_files/plots/

Producers (by filename prefix; full mapping in plan §2.7):
- ``rfdiff_*.png``     → ``RFDIFFUSION_PLOTS`` / ``bin/rfdiffusion_plots.py``
- ``rosetta_*.png``    → ``ROSETTA_FILTER_PLOTS`` / ``bin/rosetta_filter_plots.py``
- ``mpnn_*.png``       → ``MPNN_PLOTS`` / ``bin/mpnn_plots.py``
- ``negsteer_*.png``   → ``NEGSTEER_PLOTS`` / ``NEGSTEER_WITHIN_SEQUENCE_PLOTS``
                         (``bin/negsteer_plots.py`` / within-sequence variants)
- ``orthogonal_*.png`` → ``ORTHOG_PLOTS`` / ``bin/orthogonal_plots.py``

Plots emitted by the orthogonal_metrics module are absent from the
per-module reference set (the AF3/biophys/rosetta streams that feed those
plots did not run on the per-module fixture — see test_orthogonal_metrics.py).
A skip-only placeholder row keeps each ORTHOG_PLOTS plot in the
parametrize list so the reference can be extended without adding new
tests.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from tests.characterization.conftest import STAGE_ROOTS
from tests.characterization.helpers.png_compare import compare_png_perceptual


def _safe_repo_root() -> Path | None:
    here = Path(__file__).resolve()
    for parent in (here, *here.parents):
        if (parent / "nextflow.config").exists() or (parent / "main.nf").exists():
            return parent
    return None


def _enumerate_pngs() -> list[tuple[str, str]]:
    """Walk each per-stage example_output_files/plots/ and return
    (stage, filename) tuples for every PNG present. Resolved at import
    time so parametrize ids are stable across collection runs.
    """
    repo_root = _safe_repo_root()
    if repo_root is None:
        return []
    pairs: list[tuple[str, str]] = []
    for stage, stage_dir in STAGE_ROOTS.items():
        plots = repo_root / stage_dir / "example_output_files" / "plots"
        if not plots.exists():
            continue
        for p in sorted(plots.iterdir()):
            if p.is_file() and p.suffix == ".png":
                pairs.append((stage, p.name))
    return pairs


# Plots known to be produced by ORTHOG_PLOTS but absent from the per-module
# orthogonal_metrics reference set (the AF3/biophys/rosetta streams that
# feed them did not run on the per-module fixture). Kept as parametrize
# placeholders so the reference can be extended without new tests.
_ORTHOGONAL_PLOTS_PLACEHOLDERS: list[tuple[str, str]] = [
    ("orthogonal_metrics", "orthogonal_af3_vs_boltz.png"),
    ("orthogonal_metrics", "orthogonal_combined_cohort_summary.png"),
    ("orthogonal_metrics", "orthogonal_filter_cascade.png"),
    ("orthogonal_metrics", "orthogonal_metrics_vs_composite.png"),
]


_PNGS_PRESENT = _enumerate_pngs()
_PNGS = sorted(set(_PNGS_PRESENT) | set(_ORTHOGONAL_PLOTS_PLACEHOLDERS))


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize(
    ("stage", "filename"),
    _PNGS,
    ids=[f"{s}-{f}" for s, f in _PNGS],
)
def test_plot_from_plot_module(stage_reference_root, stage_output_root, stage, filename):
    """matplotlib plot module (see filename prefix in module docstring) — PNG-PERCEPTUAL."""
    ref_root = stage_reference_root(stage)
    rel = f"plots/{filename}"
    target = ref_root / rel
    if not target.exists():
        pytest.skip(
            f"absent from per-module reference set ({target}); "
            "the producer chain feeding this plot did not run on the per-module fixture."
        )
    out_root = stage_output_root(stage)
    compare_png_perceptual(target, out_root / rel).assert_passed()
