# 11 — Phase 2 Plan: Characterization Test Suite

A planning document, written before any test code exists. The goal is to produce a behavioral safety net that makes Phase 3 mechanical cleanup possible: a suite that captures the pipeline's *current* output behavior so that any divergence introduced by refactoring is visible immediately.

This document fixes the scope, structure, comparison strategy, and delivery sequence of that suite. It does *not* yet contain test code; it is the artefact to review and amend before code is written.

---

## 1. Framing

### 1.1 What these tests are

Characterization tests in the sense of [Feathers, *Working Effectively with Legacy Code*]: tests that pin current behavior so changes become visible. They assert **stability**, not **correctness**.

The reference outputs we compare against are the artefacts of a real pipeline run (the supervisor-demo run, 32 designs × 4 MPNN sequences + 2 controls, files already on disk in `tests/full_test_run/example_output_files/`). Some of those outputs almost certainly contain undetected logical errors — the inverted `scaffold` semantics flagged in glossary §F3, the `merge_orthogonal_metrics.py` semantic divergence flagged in finding A4, and the suspected `--n-cycles 1` reversion-skip bug from notes6 are all candidates. These tests will not detect those errors. They will detect any *change* in behavior, which during refactoring is the property that matters.

A test failing during refactoring does not by itself indicate the refactor is wrong. It indicates behavior changed. The change must be evaluated on its merits — and a refactor that *fixes* a latent bug will produce a failure that should be welcomed, investigated, and (after confirmation) used to update the reference output.

### 1.2 What these tests are not

- They are not unit tests. Phase 5 owns those.
- They are not correctness tests. The reference outputs are not certified.
- They are not a substitute for the cache-busting audit. If the Nextflow process cache hashes the interpolated command string rather than the script content, refactored scripts will silently produce cached old outputs and the tests will pass falsely. This is recorded as an entry condition (§6.1), not a test step.
- They do not cover Branch A (HADDOCK). Branch A is out of scope for the entire remediation pass and will be redesigned separately.

### 1.3 Two-wave structure

Phase 2 ships in two waves:

- **Wave 1 — Cohort + intermediates.** All comparison-friendly file outputs the pipeline emits, plus the per-stage intermediates that downstream stages consume. Round-trip tested against a fresh re-run on the HPC. This is the working safety net before any Phase 3 code change.
- **Wave 2 — Python-internal state.** Targeted instrumentation around the small set of complex functions whose behavior is too subtle to capture from output files alone (`align_to_native_by_anchors` is the canonical example). Layered on after Wave 1 is green.

The plan below specifies Wave 1 in full detail. Wave 2 is sketched (§5) — the precise list of functions and instrumentation strategies will be decided after Wave 1 surfaces what's actually load-bearing.

---

## 2. Scope of comparison — Wave 1

The following table enumerates every reference output that Wave 1's tests will pin. Each entry specifies the producing process/script (so the test's failure points to a known location), the comparison strategy, and whether the file exists in the current example set or needs to be sourced from HPC.

Comparison strategies are abbreviated:

- **EXACT** — byte-equal.
- **CSV-EXACT** — column names and order identical; row order identical; string columns byte-equal; numeric columns equal within `1e-6` absolute and `1e-9` relative tolerance.
- **CSV-STRUCT** — column names and order identical; row count identical; types per column identical; row-order tolerated (sort by a designated key column before comparison). Used where the producing code does not guarantee row order.
- **JSON-DEEP** — recursive comparison with numeric tolerance (same `1e-6` / `1e-9`); list element order significant unless flagged otherwise.
- **JSON-MODULO-PATHS** — JSON-DEEP after rewriting any absolute path strings (workdir hashes, HPC-specific prefixes) to a canonical form.
- **PNG-PERCEPTUAL** — file exists; size within ±5%; structural similarity (SSIM) ≥ 0.95 against reference.
- **PNG-EXISTS** — file exists and is non-empty. For the cases where pinning even SSIM is too brittle (very stochastic plot content).

### 2.1 Stage 1 — Preprocessing

| File | Producer | Strategy | Status |
|---|---|---|---|
| `preprocessing/sequences.json` | `EXTRACT_SEQUENCES` (inline Python in `modules/preprocessing.nf`) | JSON-DEEP | ✅ on disk |
| `preprocessing/processed_contigs.txt` | `RESOLVE_CONTIGS` → `bin/rfdiffusion_contigs.py` | EXACT | ✅ on disk |

**Gap to fill:** the `WRITE_DUMMY_MAPPING` placeholder JSON. Trivially reproducible (2-byte `{}`); decision: **skip pinning** — its content is degenerate and any divergence would mean the script started writing real data, which a downstream consumer would catch.

### 2.2 Stage 2 — RFDiffusion + filter

| File | Producer | Strategy | Status |
|---|---|---|---|
| `rfdiffusion/rfdiffusion_metrics.json` | `RFDIFFUSION_FILTER` → `bin/rfdiffusion_filter.py` | JSON-DEEP | ✅ on disk |
| `rfdiffusion/filter_summary.json` | same | JSON-DEEP | ✅ on disk |
| `rfdiffusion/passing_designs.txt` | same | EXACT | ✅ on disk |

**Note:** `rfdiffusion_metrics.json` contains `design_region_coords[][3]` — Cα coordinates of the de novo region. These are deterministic given fixed RFDiffusion seed *and* the upstream contigs, but RFDiffusion itself is non-deterministic across runs. Reference comparison for Wave 1 assumes the test driver pins the random seed (see §6.3). If seed-pinning proves impractical, this becomes a Wave 2 problem and falls back to JSON-STRUCT-only checks (keys exist, types correct, value ranges sane).

### 2.3 Stage 3 — Rosetta filtering

| File | Producer | Strategy | Status |
|---|---|---|---|
| `rosetta_filtering/rosetta_filter_summary.csv` | `ROSETTA_FILTER` → `bin/rosetta_filter.py` (per finding inventory) | CSV-EXACT | ⚠️ **verify on disk** |
| Per-design Rosetta `.sc` files | same | (skip — captured indirectly via summary) | n/a |

**Gap to fill:** confirm the rosetta-filtering CSV is in the example set. If not, copy from HPC.

### 2.4 Stage 4 — ProteinMPNN

| File | Producer | Strategy | Status |
|---|---|---|---|
| Per-design `mpnn_results.json` | `MPNN_RUN` → `bin/run_proteinmpnn.py` | JSON-DEEP | ⚠️ **verify on disk** |
| `top_fastas/*.fasta` | `EXTRACT_TOP` → `bin/extract_top_sequences.py` | EXACT (per file) | ⚠️ partial on disk |

**Gap to fill:** the README notes only one sample FASTA from `af2_fastas/` is present. Need the full `top_fastas/` directory and the full `af2_fastas/` directory copied from HPC.

### 2.5 Stage 5 — Negative steering core

This is the most reference-data-rich stage and the one where the gap-filling matters most.

| File | Producer | Strategy | Status |
|---|---|---|---|
| Per-sequence `aggregated_results.csv` (×34: 32 designs + 2 controls) | `bin/boltz2_iterate_steering.py cmd_collect` | CSV-EXACT | ⚠️ **verify on disk** |
| Per-sequence `passing_summary.csv` (×34) | `bin/extract_passing.py` | CSV-EXACT (with empty-body case allowed) | ⚠️ partial — empty-body example present |
| Per-sequence `raw_per_seed_results.csv` (×34) | reversion harvest | CSV-EXACT | ⚠️ **verify on disk** |
| Per-sequence `pathways.json` (×34) | iterate-steering orchestrator | JSON-MODULO-PATHS | ⚠️ partial — one example present |
| Per-sequence `contaminated.json` (×34) | reversion stage | JSON-DEEP | ⚠️ partial — only empty case present |
| Per-sequence `reversion_plan.json` (×34) | reversion stage | JSON-DEEP | ⚠️ partial — only empty case present |
| Per-sequence `reversion_results_per_seed.json` (×34) | reversion stage | JSON-DEEP | ⚠️ partial — only empty case present |
| `cross_sequence_summary.csv` | `bin/cross_sequence_summary.py::aggregate` | CSV-EXACT | ✅ on disk (presumed — confirm) |
| `cross_sequence_summary_with_interface_metrics.csv` | `bin/compute_interface_metrics.py` | CSV-EXACT | ⚠️ **verify on disk** |
| `survivor_manifest.csv` | `EXTRACT_SURVIVOR_MANIFEST` → `bin/extract_survivor_manifest.py` | CSV-EXACT modulo absolute paths | ✅ on disk |

**Gaps to fill (priority for HPC hunt):**

1. **A populated-reversion per-sequence example.** All four reversion-related JSONs in the current set are empty (`n_contaminated=0`). At least one MPNN sequence whose cold-start prediction *did* produce contaminating mutations needs to be added — this exercises the entire reversion harvest path that is currently uncovered.
2. **The two control sequences' rows.** The current `cross_sequence_summary.csv` does not include `control_scrambled` / `control_polyA`. These exercise the cold-start-only path with no contact residues → silent skip → 3-row stub (per glossary §5).
3. **All 34 per-sequence subdirectories' aggregated/passing/raw CSVs.** Confirm whether the example set has all 34 or only a subset.

### 2.6 Stage 6 — Orthogonal metrics

| File | Producer | Strategy | Status |
|---|---|---|---|
| `orthogonal_metrics/af3_nomsa_summary.csv` (cohort) | `AF3_PARSE_OUTPUT` → `bin/parse_af3_output.py` | CSV-EXACT | ✅ on disk |
| `orthogonal_metrics/biophysical_summary.csv` (cohort) | `NEGSTEER_BIOPHYSICAL_METRICS` → `bin/run_biophysical_metrics.py` | CSV-EXACT | ✅ on disk |
| `orthogonal_metrics/rosetta_summary.csv` (cohort) | `NEGSTEER_ROSETTA_METRICS` → `bin/run_rosetta_metrics.py` | CSV-EXACT | ✅ on disk |
| `orthogonal_metrics/survivors_with_orthogonal_metrics.csv` | `bin/merge_orthogonal_metrics.py` | CSV-EXACT | ⚠️ **verify on disk** |
| `merged_orthogonal_metrics.csv` | same | CSV-EXACT | ⚠️ **verify on disk** |
| Per-survivor AF3/biophys/rosetta subdirectories | per-stream modules | (skip in Wave 1 — captured indirectly via cohort merge) | gap |

**Gap to fill:** confirm `merged_orthogonal_metrics.csv` (or its actual filename — finding A4 records semantic divergence between production and test versions).

**Note flagged for later:** finding A4 — the production vs test `merge_orthogonal_metrics.py` divergence. The reference output reflects whichever version actually ran in the supervisor-demo run. The test pin will lock that behavior, *including* the bug if the production version was used. The Verification Queue gets an entry for this (§7).

### 2.7 Stage 7 — Plots

| File pattern | Count | Producer | Strategy |
|---|---:|---|---|
| `plots/rfdiff_*.png` | 5 | `RFDIFFUSION_PLOTS` | PNG-PERCEPTUAL |
| `plots/rosetta_sc_histogram.png` | 1 | `ROSETTA_FILTER_PLOTS` | PNG-PERCEPTUAL |
| `plots/mpnn_*.png` | 4 | `MPNN_PLOTS` | PNG-PERCEPTUAL |
| `plots/negsteer_*.png` (cohort) | 7 | `NEGSTEER_PLOTS` | PNG-PERCEPTUAL |
| `plots/negsteer_*.png` (within-sequence) | 5 | `NEGSTEER_WITHIN_SEQUENCE_PLOTS` | PNG-PERCEPTUAL |
| `plots/orthogonal_*.png` | 4 | `ORTHOG_PLOTS` | PNG-PERCEPTUAL |

26 PNGs total. SSIM threshold of 0.95 is a starting calibration; the test framework will allow per-plot override (some plots are inherently more stochastic than others — e.g., plots that depend on random tie-breaking in clustering).

If SSIM proves too noisy in practice, individual plots demote to PNG-EXISTS and a Verification Queue entry records the demotion. The principle: plots are weak evidence; a plot regression should prompt investigation but rarely block a refactor.

---

## 3. Test framework architecture

### 3.1 Layout

```
tests/
└── characterization/
    ├── __init__.py
    ├── conftest.py                    # pytest fixtures: paths, tolerances, helpers
    ├── pyproject.toml                 # or pytest.ini at repo root
    ├── README.md                      # how to run, how to update references
    ├── helpers/
    │   ├── __init__.py
    │   ├── csv_compare.py             # CSV-EXACT, CSV-STRUCT comparators
    │   ├── json_compare.py            # JSON-DEEP, JSON-MODULO-PATHS comparators
    │   ├── png_compare.py             # PNG-PERCEPTUAL, PNG-EXISTS comparators
    │   └── path_normalize.py          # workdir-hash → canonical-form rewriter
    ├── test_preprocessing.py
    ├── test_rfdiffusion.py
    ├── test_rosetta_filtering.py
    ├── test_mpnn.py
    ├── test_negsteer_per_sequence.py  # parameterized over the 34 sequences
    ├── test_negsteer_cohort.py
    ├── test_orthogonal_metrics.py
    └── test_plots.py
```

The reference data continues to live at `tests/full_test_run/example_output_files/`; tests resolve it via a fixture, not a hardcoded path. This separation keeps the framework portable: when the reference set is updated, only the fixture's resolution logic changes.

### 3.2 pytest conventions

- **Fixtures in `conftest.py`:**
  - `reference_root` — resolves to `tests/full_test_run/example_output_files/`. Configurable via env var `RECEPTOR_REF_ROOT` for running against alternative reference sets.
  - `output_root` — resolves to the directory of a fresh pipeline run. Configurable via env var `RECEPTOR_OUTPUT_ROOT`. Defaults to a sensible HPC location.
  - `numeric_tolerance` — returns `(abs=1e-6, rel=1e-9)`. Per-test override via marker: `@pytest.mark.tolerance(abs=1e-4, rel=1e-6)`.
  - `path_normalizer` — returns the canonicalization function for absolute paths.
  - `ssim_threshold` — returns `0.95` by default; per-plot override via marker.

- **Parameterization:** `test_negsteer_per_sequence.py` parameterizes over the 34 sequences. Each sequence's failure is a separate test, so a single-sequence regression doesn't mask others.

- **Markers:**
  - **Tier markers** (mutually exclusive — every test has exactly one):
    - `@pytest.mark.local_unit` — tests of individual functions/comparators in isolation. Fast (<1 sec each). Runnable anywhere with Python. Wave 1 comparator unit tests live here.
    - `@pytest.mark.local_integration` — tests of pure-Python characterization on saved inputs (e.g., Wave 2 snapshot tests on `align_to_native_by_anchors`). Runnable on Mac without HPC, Singularity, or Nextflow. Wave 1 has no tests here; Wave 2 populates it.
    - `@pytest.mark.hpc` — tests requiring a fresh pipeline run. HPC-only. Wave 1's stage-by-stage characterization tests live here.
  - **Selection markers** (orthogonal, multiple allowed):
    - `@pytest.mark.wave1` / `@pytest.mark.wave2` — selectable test subsets by remediation phase.
    - `@pytest.mark.requires_populated_reversion` — tests that need the populated-reversion reference (skipped until that gap is filled).
    - `@pytest.mark.requires_controls` — tests that need control rows (skipped until those are filled).
    - `@pytest.mark.slow` — for tests that compute SSIM over large images.

- **Reporting:** failures emit a unified-diff-style report for CSV/JSON, and side-by-side diff plus SSIM score for PNG. The comparator helpers own that formatting; tests just assert.

### 3.3 Invocation

Tier-aware invocation patterns. Note the three tiers correspond to where each test can run:

| Tier | Runs where | Speed | Used when |
|---|---|---|---|
| `local_unit` | Mac, anywhere with Python | <1 sec each | Every few minutes during framework development |
| `local_integration` | Mac (no HPC) | Seconds–minutes | Before transferring code to HPC; whenever Python-internal logic changes |
| `hpc` | HPC only (needs Singularity, Nextflow, GPUs) | Hours (full pipeline run + comparison) | After every batch of changes synced to HPC |

```bash
# Mac, fast comparator sanity check:
pytest -m local_unit

# Mac, before HPC transfer (everything that can run locally):
pytest -m "local_unit or local_integration"

# HPC, after a fresh pipeline run:
pytest -m hpc

# HPC, full check (everything):
pytest

# Filtering further by wave:
pytest -m "local_unit and wave1"
pytest -m "hpc and wave1 and not slow"
```

Exit code 0 means all selected tests passed. The marker scheme means a Mac invocation never accidentally triggers an HPC-only test (it would be skipped, not failed), and an HPC invocation never accidentally skips a regression-relevant local test.

### 3.4 Dependencies

Added to whatever dependency manifest the project uses (assuming `requirements.txt` or `pyproject.toml`):

- `pytest>=7.0`
- `pandas` (likely already a dependency; used for CSV comparison)
- `numpy` (likely already a dependency)
- `Pillow` (PNG loading)
- `scikit-image` (SSIM) — only one needed beyond what's likely already present

---

## 4. Comparison-helper specifications

Each comparator lives in `tests/characterization/helpers/`. These are the contracts they must satisfy.

### 4.1 `csv_compare.py`

```python
def compare_csv_exact(
    reference: Path,
    actual: Path,
    *,
    abs_tol: float = 1e-6,
    rel_tol: float = 1e-9,
    string_columns: Iterable[str] | None = None,
) -> ComparisonResult:
    """
    Strategy: CSV-EXACT.
    - Header must match exactly (names + order).
    - Row count must match.
    - Row order must match.
    - String columns: byte-equal.
    - Numeric columns: within tolerances.
    - NaN equality: NaN == NaN treated as match.
    On mismatch, ComparisonResult contains a bounded diff (first 10
    differing rows, columns highlighted).
    """

def compare_csv_struct(
    reference: Path,
    actual: Path,
    *,
    sort_by: list[str],
    abs_tol: float = 1e-6,
    rel_tol: float = 1e-9,
) -> ComparisonResult:
    """Strategy: CSV-STRUCT. Same as exact but rows sorted by sort_by first."""
```

### 4.2 `json_compare.py`

```python
def compare_json_deep(
    reference: Path,
    actual: Path,
    *,
    abs_tol: float = 1e-6,
    rel_tol: float = 1e-9,
    list_orders_significant: bool = True,
) -> ComparisonResult:
    """
    Strategy: JSON-DEEP.
    Recursive structural comparison.
    Numeric scalars compared with tolerances.
    Lists compared positionally unless list_orders_significant=False.
    Nested paths reported on mismatch (e.g., "designs[3].scaffold_rmsd").
    """

def compare_json_modulo_paths(
    reference: Path,
    actual: Path,
    path_normalizer: Callable[[str], str],
    **kwargs,
) -> ComparisonResult:
    """
    Strategy: JSON-MODULO-PATHS.
    Wraps compare_json_deep, applying path_normalizer to all string values
    that look like absolute paths before comparison.
    """
```

### 4.3 `png_compare.py`

```python
def compare_png_perceptual(
    reference: Path,
    actual: Path,
    *,
    ssim_threshold: float = 0.95,
    size_tolerance_pct: float = 5.0,
) -> ComparisonResult:
    """
    Strategy: PNG-PERCEPTUAL.
    1. Both files exist and are valid PNG.
    2. File sizes within size_tolerance_pct.
    3. SSIM (computed on decoded image arrays, RGB channel) >= ssim_threshold.
    """

def compare_png_exists(reference: Path, actual: Path) -> ComparisonResult:
    """Strategy: PNG-EXISTS. Both files exist and are non-empty."""
```

### 4.4 `path_normalize.py`

```python
def canonicalize_workdir_paths(s: str) -> str:
    """
    Replace patterns like '/work/ab/cdef1234.../' with '<WORKDIR>/'.
    Replace user-home prefixes ('/hpc-home/<user>/') with '<HOME>/'.
    Idempotent.
    """
```

The exact regex set is small but matters for false-positive avoidance. The first iteration will likely need refinement once the actual paths in the reference data are inspected.

---

## 5. Wave 2 — Python-internal state (sketched)

Wave 2 is layered on after Wave 1 is green and the cache-busting audit is complete. It targets functions whose behavior is too subtle to capture from output files alone.

**Tier classification:** Wave 2 tests are explicitly `local_integration` — they pin Python-internal state from saved inputs and run on the Mac without needing the HPC. This is one of the main payoffs of the tiered scheme: as soon as Wave 2 ships, you can iterate on the highest-complexity functions (e.g., `align_to_native_by_anchors`) with Mac-speed feedback rather than HPC-round-trip-speed feedback.

The synthesis identified one canonical case explicitly (`align_to_native_by_anchors` in `pipeline_correct_sequences.py`, complexity E=36). Wave 2 will survey the remaining complexity hotspots and decide which deserve internal-state pinning. Provisional candidates from §2 of the synthesis:

| Candidate | File | Rationale |
|---|---|---|
| `align_to_native_by_anchors` | `pipeline_correct_sequences.py` | Most subtle upstream logic; correctness hinges on it. |
| `harvest_reversion_results` | `reversion.py` | F=74 complexity; full state cycle through reversion. |
| `cmd_plan` (Boltz scheduler) | `boltz2_negative_steering.py` | F=109 — highest single-function complexity. |
| `aggregate` (cohort) | `cross_sequence_summary.py` | F=58; tier classification + composite score live here. |
| `extract_effector_template_cif` | `boltz2_negative_steering.py` | E=33; per-design template generation. |

Wave 2's *strategy* will be one of:

- **Snapshot testing** of the function's return value or side-effect file, captured on a representative input from the supervisor-demo cohort.
- **Property assertions** for functions where the exact value is non-deterministic but invariants are not (e.g., "anchor alignment preserves contact residues at fixed positions").
- **Instrumentation** — temporary stdout / file dumps added to capture intermediate state, removed after Wave 2 reference outputs are collected.

The exact list and strategies are deferred to a post-Wave-1 follow-up planning pass.

---

## 6. Entry conditions and execution sequence

### 6.1 Entry conditions (must hold before Phase 2 execution begins)

1. **Cache-busting audit complete.** Every `bin/X.py` invocation in `modules/*.nf` must use the `path script_input` pattern (or equivalent) so refactors invalidate the Nextflow cache. This is a standalone pre-Phase-2 task per the most recent direction; tracked separately. **Without this, characterization tests will silently pass on cached old outputs.**
2. **Reference set gaps inventoried.** §2 above contains a checklist of files marked ⚠️ — each must be confirmed present on the example set or sourced from HPC before its corresponding test ships.
3. **Empty `main` file removed and `bin/sequence_registry.py` flagged for deletion.** Both are housekeeping; the test suite will be run before and after the deletions in §6.3.

### 6.2 Wave 1 build sequence

1. **Scaffolding.** Create `tests/characterization/` with `conftest.py`, `pyproject.toml` (or update existing), and the four helper modules. Tests are no-ops; goal is that `pytest tests/characterization/` runs and discovers the test files.
2. **Helper unit tests.** Each comparator gets a small test verifying it correctly identifies identity (reference == reference returns pass), single-cell numeric drift (returns fail with correct location), and edge cases (empty CSV, NaN, missing file). These are real unit tests on the comparators themselves, not characterization tests.
3. **Stage-by-stage characterization tests.** Order: preprocessing → RFDiffusion → Rosetta filtering → MPNN → negsteer cohort → orthogonal metrics → plots → negsteer per-sequence. Per-sequence is last because it's the largest parameterization and shaking out the comparators on simpler tests first reduces churn.
4. **Local self-check.** Run the suite on the Mac against the reference set as both reference *and* actual. All tests must pass. (This catches comparator bugs before round-trip cost is incurred.)
5. **HPC round-trip.** Copy the suite to HPC. Run against a fresh re-run of the same `params_test.yml`. All tests should pass — failures here are *real* signal: either a non-determinism not yet accounted for, or a bug in cache-busting, or a flawed comparator. Iterate until clean.

### 6.3 First Phase 3 commit gated on the safety net

Once the suite is green on a fresh HPC run, Phase 3 begins with the smallest possible code change to validate the safety net:

1. Delete `bin/sequence_registry.py`.
2. Delete the empty `main` file at repo root.
3. Re-run pipeline on HPC.
4. Re-run characterization tests.
5. **Hand-eye plot review.** Open the 26 reference plots and the 26 fresh plots side-by-side. Visual inspection is a sanity check on the SSIM thresholds and a chance to catch regressions that pass numeric comparison but look wrong. Plots are the most human-readable view into the pipeline's outputs; suspicious changes will surface here even when CSVs pass.

If green and the plots look right, the safety net works. If not green, the failure must be diagnosed before any further Phase 3 commits. This is the canonical "first refactor commit" — chosen because both deletions are demonstrably unreachable from the current pipeline and should produce zero output diffs.

**On the pass criterion broadly:** a clean `pytest` exit on a fresh HPC run is the formal pass criterion. The hand-eye plot review is the informal sanity check that complements it. Both are required to declare a wave's tests green.

---

## 7. Risks and Verification Queue entries

The following entries should be added to `remediation_state.md`'s Verification Queue at the conclusion of this planning pass.

### 7.1 Risks specific to Phase 2

| Risk | Mitigation |
|---|---|
| Reference outputs bake in a latent bug; the fix later registers as a "regression" | Verification Queue entries (below); commit message convention noting when a Phase 3+ failure is a confirmed-bug-fix vs. a real regression. |
| Default tolerances (`abs=1e-6, rel=1e-9`) trip on legitimately stable but noisier columns (e.g., Boltz iPTM, where benign drift is closer to 1e-3) | Per-column tolerance override via marker / fixture. **Principle:** if a column repeatedly fails on benign re-runs, the response is per-column override + Verification Queue entry recording the override, **not** blanket tolerance loosening. |
| Plot SSIM threshold is wrong for some plots and produces false positives or negatives | Per-plot override + demote-to-EXISTS escape hatch + Verification Queue entry on demotion. Hand-eye plot review (per §6.3) provides an independent sanity check. |
| RFDiffusion non-determinism leaks into `rfdiffusion_metrics.json` despite seed pinning | If first HPC round-trip fails on this file, demote to JSON-STRUCT (keys + types + ranges, not values) and record the demotion. |
| The cache-busting audit regresses or is incomplete | Tests pass falsely on cached outputs. Mitigation: any post-Phase-3 commit that does not produce expected diff in tests is a red flag for cache poisoning, not a green light. |

### 7.2 Verification Queue entries to add

> - **What:** `merge_orthogonal_metrics.py` semantic divergence (production gates on AF3; test demotes to flag-only). Reference output reflects whichever version ran in the supervisor-demo.
> - **Why suspicious:** Finding A4 explicitly flags this as a latent bug. The Phase 2 test pin will lock production behavior — including the bug — until Phase 3.
> - **How to verify:** Compare the production `merge_orthogonal_metrics.py` against the test-tree copy; check the cohort's AF3 column distribution (synthesis §5 open question); if production was used, expect this test to fail when Phase 3 fixes the divergence, and update the reference at that point.

> - **What:** Suspected `--n-cycles 1` silent-skip-reversion bug from notes6.
> - **Why suspicious:** Synthesis §5 records the bug as not verified during inventory. The supervisor-demo ran with `--n-cycles 1`, so reference outputs may bake in this skipped behavior.
> - **How to verify:** Inspect a per-sequence run's `pathways.json` and reversion JSONs; cross-check whether reversion was actually attempted on contaminated sequences. If skipped where it should not have been, the reversion-related JSONs need to be regenerated from a `--n-cycles >= 2` run before Phase 3.

> - **What:** `scaffold_rmsd` field in `rfdiffusion_metrics.json` is actually motif RMSD (glossary §F3 inversion).
> - **Why suspicious:** Field name is inverted from Baker-lab convention; semantics not in question, naming is.
> - **How to verify:** Naming-only fix in Phase 3.3. Test will need to be updated at that point — reference file's field name changes, not its values.

---

## 8. What this plan is committing to

The deliverable of Phase 2, in order:

1. **This plan document** (`notes/inventory/11_phase_2_plan.md`), reviewed and amended.
2. **The reference-set gap-fill checklist from §2** completed by an HPC round-trip; README in `example_output_files/` updated with the new files.
3. **Wave 1 test suite under `tests/characterization/`,** green on the Mac (`local_unit` tier) and on the HPC against a fresh re-run (`hpc` tier). Pass criterion: clean `pytest` exit *and* hand-eye plot review.
4. **`tests/characterization/README.md`** documenting how to run the suite at each tier, the marker scheme, and the **manual reference-update procedure** (see §8.1 below).
5. **First Phase 3 commit (delete dead files), green tests post-commit,** as the safety-net validation.
6. **Wave 2 planning addendum**, written after Wave 1 ships, listing concrete Wave 2 instrumentation targets (`local_integration` tier).
7. **`remediation_state.md` updated** with completion of Phase 2 wave 1, Verification Queue entries from §7.2 added, tag `phase-2.1-wave1-complete` applied.

Wave 2 has its own deliverable cycle, which this document does not commit to in detail.

### 8.1 Manual reference-update procedure

When a Phase 3+ change legitimately alters behavior (e.g., the Stage 6 AF3 divergence fix), the reference set must be updated. The procedure is fully manual to ensure each update is a deliberate decision, not an automated overwrite:

1. Run `pytest -m "local_unit or local_integration"` on Mac. Confirm pre-update state is green for everything that doesn't depend on the affected stage.
2. Sync the change to HPC and run the pipeline.
3. Run `pytest -m hpc` on HPC. Note which tests fail and confirm the failures match the *expected* behavior change.
4. For each expected failure: manually copy the new output file from the HPC pipeline run into `tests/full_test_run/example_output_files/`, replacing the old reference.
5. Update `tests/full_test_run/example_output_files/README.md`:
   - Add an entry to a "Reference set change log" section noting the date, the Phase 3+ commit hash, the affected files, and a one-line summary of why the reference changed.
   - Update any file-description text in the README that has changed (e.g., column added, field renamed).
6. Re-run `pytest -m hpc`. The previously-failing tests should now pass against the updated reference. Any tests that *still* fail are unexpected regressions and need investigation before commit.
7. Commit the reference update and the README update together; commit message format `chore(reference): update for <commit hash> — <summary>`.

This procedure is documented in `tests/characterization/README.md` and is the canonical way to acknowledge a behavior change.

---

## 9. Decisions made during planning

The four open questions raised during planning were resolved as follows. Recorded here so this document stands alone for future readers.

1. **Tolerance calibration.** Default `abs=1e-6, rel=1e-9`. If a column repeatedly fails on benign re-runs, response is per-column override + Verification Queue entry — *not* blanket loosening. Boltz columns are expected to need overrides (drift typically <1e-3); this will be confirmed empirically on the first HPC round-trip.

2. **Pass criterion.** Clean `pytest` exit on a fresh HPC run, plus a hand-eye review of the 26 reference vs fresh plots side-by-side. Both are required to declare a wave's tests green. The plot review is the easiest human-readable check on whether numbers that pass numeric comparison still *look* right.

3. **CI scope.** No GitHub Actions. The marginal value over disciplined local invocation is low for a solo-developer remediation pass. The framework is CI-ready — `pytest -m local_unit` is single-command and would add to a workflow file in 10 minutes — so this can be revisited without refactoring if the project's collaboration profile changes.

4. **Reference-update workflow.** Fully manual, per the procedure in §8.1. The README in `example_output_files/` doubles as the change log for reference-set updates.

The three-tier marker scheme (`local_unit` / `local_integration` / `hpc`) was adopted in preference to a flat local/hpc split, so that pure-Python characterization tests on saved inputs (Wave 2 territory) have a clear home distinct from comparator unit tests.

---

## 10. Document version history

- **v1** (initial draft) — Sections 1–8 specifying scope, framework architecture, comparators, waves, entry conditions, risks, and deliverables. Section 9 listed four open questions for review.
- **v2** (current) — Open questions resolved into §9 "Decisions made." §3.2 marker scheme expanded to three tiers (`local_unit` / `local_integration` / `hpc`). §3.3 invocation examples updated. §5 Wave 2 explicitly classified as `local_integration`. §6.3 pass criterion adds hand-eye plot review. §7.1 risk table adds tolerance-handling row. §8 deliverables list adds manual reference-update procedure (§8.1) and `tests/characterization/README.md` deliverable. GitHub Actions / CI integration not adopted.
