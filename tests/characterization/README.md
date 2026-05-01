# Characterization Test Suite

This directory holds the Phase 2 characterization tests for the receptor-resurfacing pipeline. The suite pins the pipeline's *current* observable behavior so that any divergence introduced by Phase 3+ refactoring becomes visible immediately.

## What this suite is — and is not

These are **characterization tests** in the sense of Feathers, *Working Effectively with Legacy Code*. They assert **stability**, not **correctness**.

The reference outputs live under `tests/full_test_run/example_output_files/` and were captured from a real pipeline run (the supervisor-demo: 32 designs × 4 MPNN sequences + 2 controls). **They are not certified correct.** Some entries are known or suspected to encode latent bugs (e.g. inverted `scaffold` semantics from glossary §F3, the production/test divergence in `merge_orthogonal_metrics.py` flagged in finding A4, the suspected `--n-cycles 1` reversion-skip from notes6). These tests will not catch those bugs. They will catch any *change* in behavior — which during refactoring is the property that matters.

A test failing during refactoring does **not** by itself mean the refactor is wrong. It means behavior changed. Evaluate the change on its merits; a refactor that *fixes* a latent bug will produce a failure that should be welcomed, investigated, and (after confirmation) used to update the reference output via the manual procedure below.

The full framing is in `notes/inventory/11_phase_2_plan.md` §1.

## Quick-start invocations

Three tiers, mutually exclusive — every test carries exactly one:

| Tier | Where it runs | Speed | Typical use |
|---|---|---|---|
| `local_unit` | Mac / anywhere with Python | <1 s each | Comparator/helper sanity check during framework development |
| `local_integration` | Mac (no HPC needed) | seconds–minutes | Pure-Python characterization on saved inputs (Wave 2) |
| `hpc` | HPC only | hours (full pipeline run + comparison) | After every batch of changes synced to HPC |

```bash
# Mac — fast comparator sanity check
pytest -m local_unit

# Mac — everything that can run locally (before HPC transfer)
pytest -m "local_unit or local_integration"

# HPC — after a fresh pipeline run
pytest -m hpc

# HPC — full check
pytest

# Filter further by wave
pytest -m "local_unit and wave1"
pytest -m "hpc and wave1 and not slow"
```

For HPC-tier tests, point the suite at a fresh run directory:

```bash
export RECEPTOR_OUTPUT_ROOT=/path/to/fresh/pipeline/run
pytest -m hpc
```

## Layout

```
tests/
├── __init__.py
└── characterization/
    ├── __init__.py
    ├── conftest.py             # session-scope fixtures
    ├── README.md               # this file
    └── helpers/
        ├── __init__.py
        ├── result.py           # ComparisonResult dataclass
        ├── path_normalize.py   # canonicalize_workdir_paths
        └── tests/
            └── __init__.py     # unit tests for the helpers (Prompt 2)
```

Comparator modules (`csv_compare.py`, `json_compare.py`, `png_compare.py`) and the stage-by-stage characterization tests (`test_preprocessing.py`, …) are added in subsequent prompts.

## Fixtures (session-scope, defined in `conftest.py`)

| Fixture | Returns | Notes |
|---|---|---|
| `repo_root` | `Path` | Discovered by walking upward until a directory containing `nextflow.config` or `main.nf` is found. |
| `reference_root` | `Path` | Defaults to `<repo_root>/tests/full_test_run/example_output_files/`. Override with `RECEPTOR_REF_ROOT`. Skips the test if the path does not exist. |
| `output_root` | `Path` | Set via `RECEPTOR_OUTPUT_ROOT` (no default). Skips the test if unset or missing. |
| `numeric_tolerance` | `(1e-6, 1e-9)` | Default `(abs, rel)`. Per-test override: `@pytest.mark.tolerance(abs=…, rel=…)`. |
| `path_normalizer` | callable | The `canonicalize_workdir_paths` function. |
| `ssim_threshold` | `0.95` | Default; per-plot override via marker. |

## Comparison strategies (Wave 1)

| Strategy | Used for | What it checks |
|---|---|---|
| `EXACT` | Small deterministic text files | Byte-equal. |
| `CSV-EXACT` | Most CSVs | Header names + order; row count; row order; string columns byte-equal; numeric columns within `(abs=1e-6, rel=1e-9)`. NaN == NaN. |
| `CSV-STRUCT` | CSVs whose row order is not guaranteed | Same as `CSV-EXACT` but rows sorted by a designated key column first. |
| `JSON-DEEP` | Most JSONs | Recursive structural compare; numerics within tolerance; lists positional unless flagged. |
| `JSON-MODULO-PATHS` | JSONs containing absolute workdir paths | `JSON-DEEP` after applying `canonicalize_workdir_paths` to string values. |
| `PNG-PERCEPTUAL` | Most plots | Both files exist and parse; sizes within ±5%; SSIM ≥ `ssim_threshold` (default 0.95). |
| `PNG-EXISTS` | Highly-stochastic plots | Both files exist and are non-empty (escape hatch from `PNG-PERCEPTUAL`). |

Per-file mappings live in `notes/inventory/11_phase_2_plan.md` §2.

## Manual reference-update procedure

When a Phase 3+ change legitimately alters behavior (e.g. fixing the Stage 6 AF3 divergence), the reference set must be updated. The procedure is **fully manual** so each update is a deliberate decision, not an automated overwrite.

1. On Mac, run `pytest -m "local_unit or local_integration"`. Confirm the pre-update state is green for everything that does not depend on the affected stage.
2. Sync the change to HPC and run the pipeline.
3. On HPC, run `pytest -m hpc` (with `RECEPTOR_OUTPUT_ROOT` pointing at the fresh run). Note which tests fail; confirm the failures match the *expected* behavior change.
4. For each expected failure, manually copy the new output file from the HPC run into `tests/full_test_run/example_output_files/`, replacing the old reference.
5. Update `tests/full_test_run/example_output_files/README.md`:
   - Add an entry to a "Reference set change log" section noting the date, the Phase 3+ commit hash, the affected files, and a one-line summary of why the reference changed.
   - Update any file-description text that has changed (e.g. column added, field renamed).
6. Re-run `pytest -m hpc`. The previously-failing tests should now pass. Anything still failing is an unexpected regression and needs investigation before commit.
7. Commit the reference update and the README update together. Commit message format: `chore(reference): update for <commit hash> — <summary>`.

## Troubleshooting

- **"unknown marker" warnings.** All markers are registered in the root `pyproject.toml` under `[tool.pytest.ini_options].markers`. If you add a marker, register it there or pytest will fail under `--strict-markers`.
- **Tests skipped with "Reference root does not exist".** Either populate `tests/full_test_run/example_output_files/` (the data files are gitignored) or point `RECEPTOR_REF_ROOT` at a populated copy.
- **Tests skipped with "RECEPTOR_OUTPUT_ROOT is not set".** Expected on Mac: HPC-tier tests need a fresh pipeline run directory and have no useful default.
- **Comparator passes locally but fails on HPC round-trip.** This is *real* signal — either non-determinism not yet accounted for (e.g. an unset random seed), a Nextflow cache poisoning issue (cache hashes the command string rather than the script content; refactored scripts produce cached old outputs), or a comparator bug. Investigate before treating as a tolerance problem. Per plan §7.1, the response to repeated benign failures on a single column is a per-column tolerance override + Verification Queue entry, **not** a blanket loosening.
- **Plot SSIM is too noisy.** Demote the offending plot to `PNG-EXISTS`, log the demotion in the Verification Queue, and rely on the hand-eye plot review at the end of each Phase 3 commit (plan §6.3).
