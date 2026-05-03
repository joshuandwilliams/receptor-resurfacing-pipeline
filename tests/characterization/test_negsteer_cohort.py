"""Stage 5 (cohort) — Negative steering cross-sequence summary + input derivations.

Per-sequence outputs are pinned in ``test_negsteer_per_sequence.py``. This
module pins the cohort-level summaries plus the input-derivation artefacts
under ``negative_steering/controls_inputs/`` and ``negative_steering/indices/``.

The per-module fixture's ``indices/`` directory is curated to the eight
parent designs whose sequences appear in the negsteer fixture cohort
(design_0, design_3, design_27, design_28, design_42, design_44, design_55,
design_62) — each contributing two files (``_design_region.txt`` and
``_true_interface.txt``). The parametrize list resolves at import time from
whatever files are on disk.

``cross_sequence_summary_with_interface_metrics.csv`` is produced by
NEGSTEER_INTERFACE_METRICS, which runs as part of the per-module
*orthogonal_metrics* test (not the per-module negative_steering test). The
file therefore appears in the orthogonal_metrics reference set, not the
negative_steering one — see test_orthogonal_metrics.py for the actual pin.
The test below stays as a placeholder and skips.

Producers:
- ``NEGSTEER_CROSS_SEQUENCE`` → ``bin/cross_sequence_summary.py::aggregate``
  (cross_sequence_summary.csv)
- ``NEGSTEER_INTERFACE_METRICS`` → ``bin/compute_interface_metrics.py``
  (cross_sequence_summary_with_interface_metrics.csv — emitted by the
  orthogonal_metrics per-module test, not this one)
- ``DERIVE_INPUT_INDICES`` → ``bin/derive_input_indices.py``
  (controls_inputs/*.txt, indices/*.txt)
"""
from __future__ import annotations

from pathlib import Path

import pytest

from tests.characterization.conftest import STAGE_ROOTS
from tests.characterization.helpers.csv_compare import compare_csv_exact
from tests.characterization.helpers.text_compare import compare_text_exact


@pytest.fixture
def ref(stage_reference_root):
    return stage_reference_root("negative_steering")


@pytest.fixture
def out(stage_output_root):
    return stage_output_root("negative_steering")


@pytest.mark.hpc
@pytest.mark.wave1
def test_cross_sequence_summary_from_aggregate(ref, out):
    """bin/cross_sequence_summary.py::aggregate — CSV-EXACT."""
    rel = "negative_steering/cross_sequence_summary.csv"
    compare_csv_exact(ref / rel, out / rel).assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
def test_cross_sequence_summary_with_interface_metrics_from_compute_interface_metrics_py(
    ref, out,
):
    """bin/compute_interface_metrics.py — CSV-EXACT.

    Skipped — absent from the per-module negative_steering reference set.
    NEGSTEER_INTERFACE_METRICS runs in the per-module orthogonal_metrics
    test; the produced file lives at
    ``tests/orthogonal_metrics/example_output_files/negative_steering/cross_sequence_summary_with_interface_metrics.csv``.
    test_orthogonal_metrics.py owns that pin.
    """
    rel = "negative_steering/cross_sequence_summary_with_interface_metrics.csv"
    if not (ref / rel).exists():
        pytest.skip(
            f"absent from per-module reference set ({ref / rel}); "
            "produced by the per-module orthogonal_metrics test instead."
        )
    compare_csv_exact(ref / rel, out / rel).assert_passed()


def _files_in_subdir(ref_root: Path, sub: str) -> list[str]:
    base = ref_root / sub
    if not base.exists():
        return []
    return sorted(p.name for p in base.iterdir() if p.is_file() and not p.name.startswith("."))


def _safe_negsteer_ref_root() -> Path | None:
    """Resolve the negative_steering reference root at import time so
    parametrize ids are stable. Returns None if the directory is absent on
    a fresh checkout (data files are gitignored)."""
    here = Path(__file__).resolve()
    for parent in (here, *here.parents):
        if (parent / "nextflow.config").exists() or (parent / "main.nf").exists():
            candidate = parent / STAGE_ROOTS["negative_steering"] / "example_output_files"
            return candidate if candidate.exists() else None
    return None


_REF = _safe_negsteer_ref_root()
_CONTROLS_INPUTS = _files_in_subdir(_REF, "negative_steering/controls_inputs") if _REF else []
_INDICES = _files_in_subdir(_REF, "negative_steering/indices") if _REF else []


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("filename", _CONTROLS_INPUTS)
def test_controls_inputs_from_derive_input_indices_py(ref, out, filename):
    """bin/derive_input_indices.py — TEXT-EXACT (controls_inputs/)."""
    rel = f"negative_steering/controls_inputs/{filename}"
    compare_text_exact(ref / rel, out / rel).assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("filename", _INDICES)
def test_indices_from_derive_input_indices_py(ref, out, filename):
    """bin/derive_input_indices.py — TEXT-EXACT (indices/)."""
    rel = f"negative_steering/indices/{filename}"
    compare_text_exact(ref / rel, out / rel).assert_passed()
