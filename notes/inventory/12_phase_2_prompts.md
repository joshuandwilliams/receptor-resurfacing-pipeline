# Phase 2 Wave 1 Scaffolding — Prompts for Claude Code

Two prompts, executed in order on the `remediation` branch. Prompt 1 produces a working but empty framework; Prompt 2 fills in the comparators and their unit tests and establishes `TRACEABILITY.md`.

Two further prompts (3 and 4) are placeholders for later sessions:
- **Prompt 3** will write the hpc-tier characterization tests, one per producer file, populating `TRACEABILITY.md` as it goes.
- **Prompt 4** will be a standalone suspicion-finding pass that reads each producer in `TRACEABILITY.md` and writes `SUSPICIONS.md`. Splitting suspicion-finding from test-writing keeps each task focused — reading code for "what does this do" is a different cognitive job from reading code for "does this match what it claims to do."

Both prompts below treat `notes/inventory/11_phase_2_plan.md` as the authoritative reference. Claude Code should read it before writing anything.

---

## Prompt 1 — Framework foundation

Goal: a working but empty pytest framework. No comparators, no characterization tests yet — just the layout, configuration, fixtures, and small foundational modules. Acceptance: `pytest --collect-only` runs cleanly with zero tests collected.

```
We are starting Phase 2 Wave 1 of the codebase remediation. The plan
document at notes/inventory/11_phase_2_plan.md is the authoritative
reference. Please read §1 (framing), §3 (test framework architecture),
and §4 (comparison-helper specifications) before writing anything. The
state document at remediation_state.md gives broader context.

Your job in this prompt is to scaffold the empty framework. No
comparators or characterization tests yet — those come in Prompt 2.

Create the following, all derived from the plan document:

1. pyproject.toml at the repo root:
   - Declares test dependencies under an optional-dependencies "test"
     extra: pytest>=7.0, pandas>=2.0, numpy>=1.24, Pillow>=10.0,
     scikit-image>=0.21.
   - Sets requires-python = ">=3.10". The codebase targets the lowest
     container Python: 3.10 in haddock/rfdiffusion/proteinmpnn
     (3.11 in boltz2_negsteer).
   - Configures pytest with testpaths = ["tests/characterization"],
     --strict-markers and -ra in addopts, and the full marker
     registry from plan §3.2: local_unit, local_integration, hpc,
     wave1, wave2, requires_populated_reversion, requires_controls,
     slow. Each marker registered with a one-line description per
     pytest convention.

2. tests/characterization/ directory layout per plan §3.1, with empty
   __init__.py files at:
     tests/__init__.py
     tests/characterization/__init__.py
     tests/characterization/helpers/__init__.py
     tests/characterization/helpers/tests/__init__.py

3. tests/characterization/conftest.py with the five fixtures from plan
   §3.2:
     - reference_root (session scope): defaults to
       <repo_root>/tests/full_test_run/example_output_files/,
       overridable via RECEPTOR_REF_ROOT env var. pytest.skip with a
       clear message if the resolved path doesn't exist.
     - output_root (session scope): env-var-only via
       RECEPTOR_OUTPUT_ROOT, pytest.skip if unset or missing. No
       default — fresh runs live in run-specific HPC paths.
     - numeric_tolerance (session scope): returns (1e-6, 1e-9).
     - path_normalizer (session scope): returns the
       canonicalize_workdir_paths function from
       helpers.path_normalize.
     - ssim_threshold (session scope): returns 0.95.
   Add a sixth fixture, repo_root (session scope), exposing the
   discovered repo root path for tests that need it.
   Repo-root discovery: walk upward from __file__ until you find a
   directory containing either nextflow.config or main.nf. Cache the
   result at module import time. Raise RuntimeError with a clear
   message if no marker is found before the filesystem root.

4. tests/characterization/helpers/result.py: a frozen=False
   ComparisonResult dataclass with fields
     passed: bool
     reference: Path
     actual: Path
     strategy: str  # short identifier like "CSV-EXACT", "JSON-DEEP"
     message: str
     differences: list[str]  # default_factory=list, capped at ~10
                             # entries by callers
   Implement __bool__ delegating to .passed and an assert_passed()
   method that raises AssertionError with a multi-line formatted
   message on failure (strategy header, both paths, message,
   differences indented). Comparators in Prompt 2 will all return
   this type; their unit tests will both check .passed and call
   .assert_passed().

5. tests/characterization/helpers/path_normalize.py exporting one
   function:
     canonicalize_workdir_paths(s: str) -> str
   The function replaces:
     - Nextflow workdir prefixes …/work/AB/HASH (where AB is two
       hex chars and HASH is 32 hex chars) with <WORKDIR>.
     - /hpc-home/<user>, /Users/<user>, /home/<user> with <HOME>.
   The function must be idempotent.

   IMPORTANT corrections to the plan document for this module:
   - Plan §4.4 wording mentions a "30-char hash" — that is wrong.
     Nextflow workdir hashes are 32 hex characters. Use {32} in the
     regex.
   - Pattern ordering matters: the workdir-hash regex MUST run
     BEFORE the /hpc-home, /Users, /home patterns. The workdir
     pattern is the more specific match and can include a /hpc-home
     prefix as part of its match; if the home-stripping patterns
     fire first, they leave behind a string that no longer contains
     the /hpc-home/.../work/AB/HASH shape, and the workdir pattern
     can no longer fire.

6. tests/characterization/README.md covering:
   - What the suite is (characterization tests, not correctness
     tests — see plan §1.1 for the framing). Include the caveat
     about reference outputs not being certified-correct.
   - Quick-start invocation for each tier (local_unit,
     local_integration, hpc).
   - Layout map.
   - Fixture summary.
   - Comparison strategy table from plan §2.
   - The manual reference-update procedure from plan §8.1.
   - Brief troubleshooting section.

Acceptance criteria:
   - `pytest --collect-only -q` from the repo root runs without
     errors and reports zero tests collected (correct — no tests
     written yet).
   - `pytest -m local_unit` and `pytest -m hpc` both run cleanly,
     reporting that all collected tests were deselected (zero
     selected).
   - No "unknown marker" warnings.
   - These imports succeed:
       python -c "from tests.characterization.helpers.result import ComparisonResult"
       python -c "from tests.characterization.helpers.path_normalize import canonicalize_workdir_paths"
   - canonicalize_workdir_paths is verified manually to be
     idempotent on at least one example string with both a workdir
     hash and an /hpc-home prefix.

Stop after acceptance is met. Commit with message
"Phase 2 Wave 1: scaffold characterization framework". Do not
proceed to comparators in this prompt.
```

---

## Prompt 2 — Comparators and unit tests

Run after Prompt 1 is committed and acceptance verified.

Goal: the four comparator modules with full unit test coverage, plus an empty `TRACEABILITY.md` ready for Prompt 3 to populate. Acceptance: `pytest -m local_unit` passes; `TRACEABILITY.md` exists with the documented schema and zero rows.

```
Continuing Phase 2 Wave 1. Prompt 1 created the framework foundation;
this prompt fills in the four comparator modules, writes their unit
tests, and establishes TRACEABILITY.md.

Read plan §4 (comparison-helper specifications) for the contracts each
comparator must satisfy. Tier all the unit tests in this prompt as
@pytest.mark.local_unit. Place them in
tests/characterization/helpers/tests/ — one test file per comparator
module.

A note on test organisation. These unit tests test the *comparators
themselves*, not the pipeline. They exist so that when Prompt 3 uses
the comparators to pin pipeline outputs, we already know the
comparators detect the things they claim to detect. Coverage targets:
identity passes, drift inside tolerance passes, drift outside fails,
NaN equality, missing files (reference and actual), header mismatches,
type mismatches where applicable, edge cases (empty CSVs, list-order
sensitivity, etc).

Test naming convention for this prompt:
  test_<comparator_function>_<scenario>
e.g. test_compare_csv_exact_drift_outside_tolerance_fails. Keep the
test name unambiguous — pytest's failure output uses it directly.
Each test gets a one-line docstring explaining what scenario it
exercises.

Create the following:

1. tests/characterization/helpers/csv_compare.py exporting:
     compare_csv_exact(reference, actual, *, abs_tol=1e-6,
                       rel_tol=1e-9, string_columns=None)
                       -> ComparisonResult
     compare_csv_struct(reference, actual, *, sort_by, abs_tol=1e-6,
                       rel_tol=1e-9, string_columns=None)
                       -> ComparisonResult
   Numeric tolerance: abs(a - b) <= abs_tol + rel_tol * abs(b).
   NaN equals NaN. Inf equals Inf of same sign. Use pandas to load.
   On failure, populate differences with row/column locations of the
   first ~10 mismatches; cap the list to keep failure output
   digestible. Header check (names + order) is its own failure mode,
   distinct from cell-mismatch.

2. tests/characterization/helpers/json_compare.py exporting:
     compare_json_deep(reference, actual, *, abs_tol=1e-6,
                       rel_tol=1e-9, list_orders_significant=True)
                       -> ComparisonResult
     compare_json_modulo_paths(reference, actual, path_normalizer, *,
                              abs_tol=1e-6, rel_tol=1e-9,
                              list_orders_significant=True)
                              -> ComparisonResult
   Recursive comparison. Report nested locations on mismatch like
   "designs[3].scaffold_rmsd". When list_orders_significant is False,
   compare lists as unordered (sort by JSON repr of each element
   before pairwise compare).

   Critical detail to handle correctly: Python booleans are a subclass
   of int (isinstance(True, int) is True), but JSON true must NOT
   match JSON 1. Type-check bool separately and treat bool-vs-non-bool
   as a type mismatch.

3. tests/characterization/helpers/png_compare.py exporting:
     compare_png_perceptual(reference, actual, *, ssim_threshold=0.95,
                           size_tolerance_pct=5.0) -> ComparisonResult
     compare_png_exists(reference, actual) -> ComparisonResult
   Three-stage check for PERCEPTUAL: (a) both files exist and are
   non-empty valid PNG; (b) file sizes within size_tolerance_pct
   (computed as percentage of the larger size); (c) SSIM on RGB pixel
   arrays >= threshold. Use scikit-image's structural_similarity with
   channel_axis=-1 and data_range=255. If actual's dimensions differ
   from reference's, resize via Pillow's BILINEAR before SSIM.

   Calibration note for the docstring: the 5% size-tolerance default
   is calibrated for real matplotlib PNGs (typically 50-500 KB). On
   very small procedurally-generated test images, percentage size
   differences inflate because PNG compression is sensitive to
   perturbation in compact images. This is a property of small
   images, not a defect. Unit tests in this prompt that need to
   verify the SSIM path on small images should pass an explicit
   relaxed size_tolerance_pct.

4. tests/characterization/helpers/tests/test_csv_compare.py — unit
   tests for compare_csv_exact and compare_csv_struct. Generate
   tiny CSVs in tmp_path; do not depend on any reference data.
   Cover at minimum: identity-passes, drift-within-tolerance-passes,
   drift-outside-tolerance-fails, missing-column-header-fails,
   reordered-header-fails, row-count-mismatch-fails,
   row-order-significant-for-exact, NaN-equals-NaN,
   missing-reference-file-fails, missing-actual-file-fails,
   empty-body (header-only CSV), assert_passed-raises-on-failure,
   csv_struct-tolerates-row-reordering,
   csv_struct-still-detects-value-drift,
   csv_struct-invalid-sort-key-fails.

5. tests/characterization/helpers/tests/test_json_compare.py — unit
   tests for compare_json_deep and compare_json_modulo_paths. Cover
   at minimum: identity-passes, numeric-drift-within/outside-
   tolerance, nested-path-reported-on-mismatch, missing-key-fails,
   unexpected-key-fails, list-length-mismatch-fails,
   list-order-significant-by-default,
   list-order-insensitive-when-disabled, type-mismatch-fails,
   bool-distinct-from-int, null-handling, null-vs-value-fails,
   modulo-paths-workdir-hash-normalised,
   modulo-paths-real-data-diff-still-caught, missing-files.

6. tests/characterization/helpers/tests/test_png_compare.py — unit
   tests for compare_png_perceptual and compare_png_exists. Generate
   tiny PNGs procedurally with numpy + Pillow. Cover at minimum:
   identity-passes, minor-noise-passes-with-relaxed-size-tolerance
   (use a plot-like base image: white background, a few axes-like
   lines, scattered marker dots — closer to matplotlib output than
   a smooth gradient or random noise),
   completely-different-images-fail, size-outside-tolerance-fails,
   ssim-threshold-override, missing-reference-fails,
   empty-actual-fails, png_exists-passes-when-present,
   png_exists-fails-when-empty, png_exists-fails-when-missing.

7. tests/characterization/helpers/tests/test_path_normalize.py — unit
   tests for canonicalize_workdir_paths. Cover: workdir-hash-
   replaced (expect <WORKDIR>/...), workdir-with-hpc-home-prefix
   (the hpc-home prefix gets absorbed into the workdir match,
   leaving <WORKDIR>/...), hpc-home-replaced, macos-home-replaced,
   linux-home-replaced, idempotent (apply twice == apply once),
   non-path-strings-pass-through (test on contig strings, sequence
   names, plain numbers — they must come back unchanged),
   short-hash-not-matched (a 30-char hex string in a workdir-shaped
   path must NOT match — guards against the regex matching short
   hex-like substrings).

8. tests/characterization/helpers/tests/test_result.py — small
   tests for ComparisonResult: truthy-when-passed, falsy-when-
   failed, assert_passed-noop-on-pass, assert_passed-raises-with-
   full-message-on-failure (verify message contains strategy,
   both paths, the message text, and the differences).

9. tests/characterization/TRACEABILITY.md — a markdown file with
   the schema documented but zero data rows. Format:

     # Characterization Test Traceability

     This document maps each characterization test to the pipeline
     producer it pins. It is the canonical answer to "if I touch
     <file>, which tests should I expect to flicker?"

     This file is updated as hpc-tier and local_integration-tier
     tests are added. Comparator unit tests (local_unit) are NOT
     listed here — they test the comparators themselves, not the
     pipeline.

     ## Tests

     | Test (file::function) | Producer (file::function) | Output pinned | Strategy | Notes |
     |---|---|---|---|---|

     (rows added as tests are written)

     ## Conventions

     - Test naming: `test_<output>_from_<producer>` where producer
       is the .py file or function name.
     - Producer column: file path relative to repo root, optionally
       with ::function suffix.
     - Output pinned: relative path from reference_root.
     - Strategy: one of CSV-EXACT, CSV-STRUCT, JSON-DEEP,
       JSON-MODULO-PATHS, PNG-PERCEPTUAL, PNG-EXISTS.

Acceptance criteria:
   - `pytest -m local_unit` from the repo root passes with zero
     failures and zero warnings.
   - `pytest -m local_unit -v` shows test names that follow the
     naming convention (test_<comparator>_<scenario>).
   - Each comparator function has a clear docstring with its
     contract.
   - tests/characterization/TRACEABILITY.md exists with the schema
     above; the table has zero data rows (correct — no
     pipeline-pinning tests written yet).
   - Before commit, run the full suite and confirm `pytest` (no
     marker filter) gives the same passing count as
     `pytest -m local_unit` (no other tiers have tests yet, so
     they should match).

Stop after acceptance is met. Commit with message
"Phase 2 Wave 1: comparators, unit tests, traceability skeleton".

Do not proceed to writing hpc-tier characterization tests in this
prompt — that is Prompt 3.
```

---

## Notes for the human (you)

**On running these prompts.** Prompt 1 is small and should produce a clean acceptance check on the first run; if it doesn't, the framework foundation has a structural issue worth diagnosing before Prompt 2. Prompt 2 is larger and Claude Code may iterate against `pytest -m local_unit` before all tests pass — that's expected.

**On the two corrections to the plan document.** Prompt 1 calls these out explicitly so Claude Code uses 32-char workdir hashes and the right pattern ordering. After both prompts succeed, the plan document at §4.4 should be amended to fix the "30-char hash" wording. That's a one-line edit best done as part of the same commit that updates `remediation_state.md`.

**On verification.** When Claude Code finishes each prompt, you should:
- Read the diff yourself before accepting it. The acceptance criteria are necessary but not sufficient — they verify the code runs, not that it does what it should.
- Run `pytest -m local_unit -v` and skim the test names. If a test name doesn't match its docstring, or the failure message format on a deliberately-broken test is unhelpful, push back.
- Spot-check one comparator implementation against plan §4 to make sure the contract is faithfully implemented (not just plausible-looking).

**On Prompts 3 and 4.** Both are deferred to future sessions:
- **Prompt 3** writes hpc-tier characterization tests. Should be drafted only after the cache-busting audit and the gap-fill HPC hunt (per plan §6.1) are complete, because both inform what tests can actually run.
- **Prompt 4** is the suspicion-finding pass. Reads each producer in `TRACEABILITY.md` (populated by Prompt 3) and writes `SUSPICIONS.md`. Standalone task; no test code touched. Producing this in its own prompt — not bundled into Prompt 3 — keeps the cognitive job clean: writing tests is "what does this code do," finding suspicions is "does it match what it claims."

**On the broader plan amendment.** Two updates to the plan document worth doing alongside Prompt 1:
1. §4.4 — fix "30-char hash" to "32-char hash" (corrected in the prompt).
2. §6 entry conditions — add "Prompts 1 and 2 complete (framework + comparators) before Prompt 3 begins."

Both are small enough that you might prefer to do them in the same commit that updates `remediation_state.md` after Prompt 2 passes acceptance.
