from __future__ import annotations

import os
from pathlib import Path
from typing import Callable

import pytest

from tests.characterization.helpers.path_normalize import canonicalize_workdir_paths


def _discover_repo_root() -> Path:
    markers = ("nextflow.config", "main.nf")
    here = Path(__file__).resolve()
    for candidate in (here, *here.parents):
        if candidate.is_dir() and any((candidate / m).exists() for m in markers):
            return candidate
    raise RuntimeError(
        "Could not locate repo root: walked upward from "
        f"{here} without finding any of {markers}."
    )


_REPO_ROOT = _discover_repo_root()


# Per-stage directories, relative to repo_root. Each stage owns
# `<stage_dir>/example_output_files/` (the reference set) and
# `<stage_dir>/receptor_resurfacing_results/` (the per-module test's output).
STAGE_ROOTS: dict[str, str] = {
    "rfdiffusion":        "tests/rfdiffusion",
    "rosetta_filtering":  "tests/rosetta_filtering",
    "proteinmpnn":        "tests/proteinmpnn",
    "negative_steering":  "tests/negative_steering",
    "orthogonal_metrics": "tests/orthogonal_metrics",
}


@pytest.fixture(scope="session")
def repo_root() -> Path:
    return _REPO_ROOT


@pytest.fixture(scope="session")
def stage_reference_root() -> Callable[[str], Path]:
    """Factory fixture: stage name → per-stage reference_root.

    The returned callable resolves ``tests/<stage>/example_output_files/``
    and pytest.skip's with a clear message if the directory is absent.
    Each test file binds its own stage via a thin local fixture, e.g.::

        @pytest.fixture
        def ref(stage_reference_root):
            return stage_reference_root("negative_steering")
    """

    def resolve(stage: str) -> Path:
        if stage not in STAGE_ROOTS:
            raise KeyError(
                f"Unknown stage {stage!r}; valid stages: {sorted(STAGE_ROOTS)}"
            )
        path = _REPO_ROOT / STAGE_ROOTS[stage] / "example_output_files"
        if not path.exists():
            pytest.skip(
                f"Per-stage reference root for {stage!r} does not exist: {path}. "
                f"Populate {STAGE_ROOTS[stage]}/example_output_files/ from a per-module "
                "test run before exercising hpc-tier tests."
            )
        return path

    return resolve


@pytest.fixture(scope="session")
def stage_output_root() -> Callable[[str], Path]:
    """Factory fixture: stage name → per-stage output_root.

    The returned callable resolves ``tests/<stage>/receptor_resurfacing_results/``
    and pytest.skip's with a clear message if the directory is absent
    (output_root is only present after the per-module test has run on HPC).
    """

    def resolve(stage: str) -> Path:
        if stage not in STAGE_ROOTS:
            raise KeyError(
                f"Unknown stage {stage!r}; valid stages: {sorted(STAGE_ROOTS)}"
            )
        path = _REPO_ROOT / STAGE_ROOTS[stage] / "receptor_resurfacing_results"
        if not path.exists():
            pytest.skip(
                f"Per-stage output root for {stage!r} does not exist: {path}. "
                f"Run tests/{stage}/run_test_{stage}.slurm.sh on HPC and sync the "
                f"resulting receptor_resurfacing_results/ directory before exercising "
                "hpc-tier tests."
            )
        return path

    return resolve


# ---------------------------------------------------------------------------
# Deprecated aliases — kept only so files not yet migrated to the per-stage
# fixtures don't ImportError. To be removed once
# tests/full_test_run/example_output_files/ is retired.
# ---------------------------------------------------------------------------


@pytest.fixture(scope="session")
def reference_root() -> Path:
    override = os.environ.get("RECEPTOR_REF_ROOT")
    path = (
        Path(override)
        if override
        else _REPO_ROOT / "tests" / "full_test_run" / "example_output_files"
    )
    if not path.exists():
        pytest.skip(
            f"Deprecated reference_root does not exist: {path}. "
            "Migrate the test to stage_reference_root or set RECEPTOR_REF_ROOT."
        )
    return path


@pytest.fixture(scope="session")
def output_root() -> Path:
    override = os.environ.get("RECEPTOR_OUTPUT_ROOT")
    if not override:
        pytest.skip(
            "RECEPTOR_OUTPUT_ROOT is not set; deprecated output_root fixture "
            "needs an explicit path. Migrate the test to stage_output_root."
        )
    path = Path(override)
    if not path.exists():
        pytest.skip(f"RECEPTOR_OUTPUT_ROOT points to a missing path: {path}.")
    return path


@pytest.fixture(scope="session")
def numeric_tolerance() -> tuple[float, float]:
    return (1e-6, 1e-9)


@pytest.fixture(scope="session")
def path_normalizer():
    return canonicalize_workdir_paths


@pytest.fixture(scope="session")
def ssim_threshold() -> float:
    return 0.95
