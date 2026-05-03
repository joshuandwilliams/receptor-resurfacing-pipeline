"""Stage 6 — Orthogonal metrics.

Pins the cohort-level orthogonal-metrics outputs for the survivors that
passed Stage 5. Per-survivor summary CSVs (``af3_nomsa/<seq>/af3_nomsa_summary.csv``,
``biophysical/<seq>/biophysical_summary.csv``, ``rosetta/<seq>/rosetta_summary.csv``)
are *not* pinned individually in Wave 1 — they are captured indirectly by
``survivors_with_orthogonal_metrics.csv`` (which merges them via
``bin/merge_orthogonal_metrics.py``).

The per-module orthogonal_metrics test on the curated fixture only ran
NEGSTEER_INTERFACE_METRICS, EXTRACT_SURVIVOR_MANIFEST, and AF3_SETUP_DB —
the AF3/biophysical/rosetta streams need GPU resources that the per-module
test node does not provide. As a result,
``orthogonal_metrics/survivors_with_orthogonal_metrics.csv`` is absent from
the per-module reference set; the test for it skips with a clear message
and remains as a placeholder for when the reference set is extended.

NEGSTEER_INTERFACE_METRICS publishes its output back into
``negative_steering/cross_sequence_summary_with_interface_metrics.csv``
(the per-module orthogonal_metrics test runs the merge with negsteer
inputs). That file therefore appears under this stage's reference root
even though its publish prefix says ``negative_steering/``.

Producers:
- ``EXTRACT_SURVIVOR_MANIFEST`` → ``bin/extract_survivor_manifest.py``
  (survivor_manifest.csv)
- ``NEGSTEER_INTERFACE_METRICS`` → ``bin/compute_interface_metrics.py``
  (cross_sequence_summary_with_interface_metrics.csv — published under
  negative_steering/ prefix)
- ``NEGSTEER_ORTHOGONAL_METRICS`` → ``bin/merge_orthogonal_metrics.py``
  (survivors_with_orthogonal_metrics.csv — absent from the per-module
  fixture; pin will activate once AF3/biophys/rosetta streams run)

Both files contain absolute Nextflow workdir/results paths in cells. The
shared ``compare_csv_exact_modulo_paths`` helper applies
``canonicalize_workdir_paths`` to every string cell before comparison —
analogous to JSON-MODULO-PATHS for CSVs. Plan §2.6 flagged this exact need.
"""
from __future__ import annotations

import pytest

from tests.characterization.helpers.csv_compare import (
    compare_csv_exact,
    compare_csv_exact_modulo_paths,
)


@pytest.fixture
def ref(stage_reference_root):
    return stage_reference_root("orthogonal_metrics")


@pytest.fixture
def out(stage_output_root):
    return stage_output_root("orthogonal_metrics")


@pytest.mark.hpc
@pytest.mark.wave1
def test_survivor_manifest_from_extract_survivor_manifest_py(ref, out, path_normalizer):
    """bin/extract_survivor_manifest.py — CSV-EXACT-MODULO-PATHS.

    Columns ``canonical_pdb_abs``, ``ground_truth_abs``,
    ``effector_template_cif_abs`` carry absolute paths with workdir hashes.
    """
    rel = "orthogonal_metrics/survivor_manifest.csv"
    compare_csv_exact_modulo_paths(
        ref / rel, out / rel, path_normalizer,
    ).assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
def test_cross_sequence_summary_with_interface_metrics_from_compute_interface_metrics_py(
    ref, out,
):
    """bin/compute_interface_metrics.py — CSV-EXACT.

    NEGSTEER_INTERFACE_METRICS publishes under the ``negative_steering/``
    prefix even though the producer runs as part of the orthogonal_metrics
    per-module test. The reference set lives at this stage's root.
    """
    rel = "negative_steering/cross_sequence_summary_with_interface_metrics.csv"
    target = ref / rel
    if not target.exists():
        pytest.skip(
            f"absent from per-module reference set ({target}); "
            "NEGSTEER_INTERFACE_METRICS did not publish under this prefix."
        )
    compare_csv_exact(target, out / rel).assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
def test_survivors_with_orthogonal_metrics_from_merge_orthogonal_metrics_py(
    ref, out, path_normalizer,
):
    """bin/merge_orthogonal_metrics.py — CSV-EXACT-MODULO-PATHS.

    Skipped — absent from per-module reference set. AF3/biophys/rosetta
    streams require GPU resources that the per-module test node does not
    provide; ``survivors_with_orthogonal_metrics.csv`` is not produced. Pin
    activates once the reference set is extended (see
    notes/inventory/15_discovery_run_path_coverage.md §Orthogonal metrics
    Coverage gap 3).
    """
    rel = "orthogonal_metrics/survivors_with_orthogonal_metrics.csv"
    target = ref / rel
    if not target.exists():
        pytest.skip(
            f"absent from per-module reference set ({target}); "
            "AF3/biophys/rosetta streams did not run on the per-module fixture."
        )
    compare_csv_exact_modulo_paths(
        target, out / rel, path_normalizer,
    ).assert_passed()
