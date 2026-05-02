"""Stage 1 — Preprocessing.

Pins the outputs of ``modules/preprocessing.nf``:
- ``RESOLVE_CONTIGS`` → ``bin/rfdiffusion_contigs.py`` produces ``processed_contigs.txt``.
- ``EXTRACT_SEQUENCES`` (inline Python in ``modules/preprocessing.nf``) produces ``sequences.json``.
"""
from __future__ import annotations

import pytest

from tests.characterization.helpers.json_compare import compare_json_deep
from tests.characterization.helpers.text_compare import compare_text_exact


@pytest.mark.hpc
@pytest.mark.wave1
def test_processed_contigs_from_rfdiffusion_contigs_py(reference_root, output_root):
    """bin/rfdiffusion_contigs.py — TEXT-EXACT."""
    rel = "preprocessing/processed_contigs.txt"
    compare_text_exact(reference_root / rel, output_root / rel).assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
def test_sequences_json_from_extract_sequences(reference_root, output_root):
    """EXTRACT_SEQUENCES inline Python in modules/preprocessing.nf — JSON-DEEP."""
    rel = "preprocessing/sequences.json"
    compare_json_deep(reference_root / rel, output_root / rel).assert_passed()
