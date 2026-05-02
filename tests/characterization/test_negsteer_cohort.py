"""Stage 5 (cohort) — Negative steering cross-sequence summary + input derivations.

Per-sequence outputs are pinned in ``test_negsteer_per_sequence.py``. This
module pins the cohort-level summaries plus the input-derivation artefacts
under ``negative_steering/controls_inputs/`` and ``negative_steering/indices/``.

Producers:
- ``NEGSTEER_CROSS_SEQUENCE`` → ``bin/cross_sequence_summary.py::aggregate``
  (cross_sequence_summary.csv)
- ``NEGSTEER_INTERFACE_METRICS`` → ``bin/compute_interface_metrics.py``
  (cross_sequence_summary_with_interface_metrics.csv)
- ``DERIVE_INPUT_INDICES`` → ``bin/derive_input_indices.py``
  (controls_inputs/*.txt, indices/*.txt)
"""
from __future__ import annotations

from pathlib import Path

import pytest

from tests.characterization.helpers.csv_compare import compare_csv_exact
from tests.characterization.helpers.text_compare import compare_text_exact


@pytest.mark.hpc
@pytest.mark.wave1
def test_cross_sequence_summary_from_aggregate(reference_root, output_root):
    """bin/cross_sequence_summary.py::aggregate — CSV-EXACT."""
    rel = "negative_steering/cross_sequence_summary.csv"
    compare_csv_exact(reference_root / rel, output_root / rel).assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
def test_cross_sequence_summary_with_interface_metrics_from_compute_interface_metrics_py(
    reference_root, output_root,
):
    """bin/compute_interface_metrics.py — CSV-EXACT."""
    rel = "negative_steering/cross_sequence_summary_with_interface_metrics.csv"
    compare_csv_exact(reference_root / rel, output_root / rel).assert_passed()


def _controls_inputs_files(ref_root: Path) -> list[str]:
    base = ref_root / "negative_steering" / "controls_inputs"
    if not base.exists():
        return []
    return sorted(p.name for p in base.iterdir() if p.is_file() and not p.name.startswith("."))


def _indices_files(ref_root: Path) -> list[str]:
    base = ref_root / "negative_steering" / "indices"
    if not base.exists():
        return []
    return sorted(p.name for p in base.iterdir() if p.is_file() and not p.name.startswith("."))


# Resolve at import time so parametrize ids are stable. The reference root may
# not exist on a fresh checkout (the data files are gitignored); in that case
# we leave the parametrize lists empty and pytest skips with "no tests
# collected", consistent with how the reference_root fixture behaves.
def _safe_ref_root() -> Path | None:
    import os

    override = os.environ.get("RECEPTOR_REF_ROOT")
    if override:
        candidate = Path(override)
    else:
        here = Path(__file__).resolve()
        for parent in (here, *here.parents):
            if (parent / "nextflow.config").exists() or (parent / "main.nf").exists():
                candidate = parent / "tests" / "full_test_run" / "example_output_files"
                break
        else:
            return None
    return candidate if candidate.exists() else None


_REF = _safe_ref_root()
_CONTROLS_INPUTS = _controls_inputs_files(_REF) if _REF else []
_INDICES = _indices_files(_REF) if _REF else []


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("filename", _CONTROLS_INPUTS)
def test_controls_inputs_from_derive_input_indices_py(reference_root, output_root, filename):
    """bin/derive_input_indices.py — TEXT-EXACT (controls_inputs/)."""
    rel = f"negative_steering/controls_inputs/{filename}"
    compare_text_exact(reference_root / rel, output_root / rel).assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("filename", _INDICES)
def test_indices_from_derive_input_indices_py(reference_root, output_root, filename):
    """bin/derive_input_indices.py — TEXT-EXACT (indices/)."""
    rel = f"negative_steering/indices/{filename}"
    compare_text_exact(reference_root / rel, output_root / rel).assert_passed()
