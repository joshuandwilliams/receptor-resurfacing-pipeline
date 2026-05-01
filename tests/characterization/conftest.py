from __future__ import annotations

import os
from pathlib import Path

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


@pytest.fixture(scope="session")
def repo_root() -> Path:
    return _REPO_ROOT


@pytest.fixture(scope="session")
def reference_root() -> Path:
    override = os.environ.get("RECEPTOR_REF_ROOT")
    path = Path(override) if override else _REPO_ROOT / "tests" / "full_test_run" / "example_output_files"
    if not path.exists():
        pytest.skip(
            f"Reference root does not exist: {path}. "
            "Set RECEPTOR_REF_ROOT or populate the default location."
        )
    return path


@pytest.fixture(scope="session")
def output_root() -> Path:
    override = os.environ.get("RECEPTOR_OUTPUT_ROOT")
    if not override:
        pytest.skip(
            "RECEPTOR_OUTPUT_ROOT is not set; HPC-tier tests need a fresh pipeline run directory."
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
