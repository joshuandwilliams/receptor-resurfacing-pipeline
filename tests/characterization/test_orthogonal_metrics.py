"""Stage 6 — Orthogonal metrics.

Pins the cohort-level orthogonal-metrics outputs for the survivors that
passed Stage 5. Per-survivor summary CSVs (``af3_nomsa/<seq>/af3_nomsa_summary.csv``,
``biophysical/<seq>/biophysical_summary.csv``, ``rosetta/<seq>/rosetta_summary.csv``)
are *not* pinned individually in Wave 1 — they are captured indirectly by
``survivors_with_orthogonal_metrics.csv`` (which merges them via
``bin/merge_orthogonal_metrics.py``).

Producers:
- ``EXTRACT_SURVIVOR_MANIFEST`` → ``bin/extract_survivor_manifest.py``
  (survivor_manifest.csv)
- ``NEGSTEER_ORTHOGONAL_METRICS`` → ``bin/merge_orthogonal_metrics.py``
  (survivors_with_orthogonal_metrics.csv)

Both files contain absolute Nextflow workdir/results paths in cells. The
shared ``compare_csv_exact_modulo_paths`` helper applies
``canonicalize_workdir_paths`` to every string cell before comparison —
analogous to JSON-MODULO-PATHS for CSVs. Plan §2.6 flagged this exact need.
"""
from __future__ import annotations

import pytest

from tests.characterization.helpers.csv_compare import compare_csv_exact_modulo_paths


@pytest.mark.hpc
@pytest.mark.wave1
def test_survivor_manifest_from_extract_survivor_manifest_py(
    reference_root, output_root, path_normalizer,
):
    """bin/extract_survivor_manifest.py — CSV-EXACT-MODULO-PATHS.

    Columns ``canonical_pdb_abs``, ``ground_truth_abs``,
    ``effector_template_cif_abs`` carry absolute paths with workdir hashes.
    """
    rel = "orthogonal_metrics/survivor_manifest.csv"
    compare_csv_exact_modulo_paths(
        reference_root / rel, output_root / rel, path_normalizer,
    ).assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
def test_survivors_with_orthogonal_metrics_from_merge_orthogonal_metrics_py(
    reference_root, output_root, path_normalizer,
):
    """bin/merge_orthogonal_metrics.py — CSV-EXACT-MODULO-PATHS.

    Carries ``representative_canonical_pdb`` (workdir-hashed path) plus
    several ``source_*`` columns. Plan §2.6 flagged this file's path content
    explicitly.
    """
    rel = "orthogonal_metrics/survivors_with_orthogonal_metrics.csv"
    compare_csv_exact_modulo_paths(
        reference_root / rel, output_root / rel, path_normalizer,
    ).assert_passed()
