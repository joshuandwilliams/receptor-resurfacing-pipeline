# Phase 2.9 — Prompts for Claude Code

Two prompts, executed in order on the `remediation` branch.

**Prompt A** is a small housekeeping fix: add `#SBATCH --chdir` to all five
`run_test_*.slurm.sh` scripts and delete the stray `slurm_*.out/.err` files that
landed at the repo root during earlier test runs. Fast and low-risk; commit before
Prompt B.

**Prompt B** is the main 2.9 work: restructure the 184 hpc-tier characterisation
tests to use the per-module reference sets, expand the negsteer sequence parametrize
list from 3 to 10, add the RFDiffusion fail-branch unit test, and update
TRACEABILITY.md throughout.

The active plan is `notes/inventory/14_phase_2_revision_per_module_tests.md`. The
original plan at `notes/inventory/11_phase_2_plan.md` is a historical record only.
The state document is `notes/remediation_state.md`.

---

## Prompt A — SLURM housekeeping

Goal: SLURM log files always land inside the test directory regardless of where
`sbatch` is invoked. Stray `slurm_*.out/.err` files at the repo root are deleted.

```
Small housekeeping task before the main Phase 2.9 work.

Context: when a `run_test_<module>.slurm.sh` script is submitted with
`sbatch` from the repo root rather than from inside the test directory,
SLURM resolves the relative `--output=slurm_%j.out` path against the
submission directory and writes the log to the repo root. This has happened;
there are stray slurm log files at the repo root that should not be there.

The fix is to add `#SBATCH --chdir=<test_dir>` to each script so SLURM
always changes into the test directory before resolving output paths,
regardless of where `sbatch` is called from.

STEP 1 — Add `#SBATCH --chdir` to each of the five SLURM scripts.

The five scripts and their test directories:

  tests/rfdiffusion/run_test_rfdiffusion.slurm.sh
    → /hpc-home/jowillia/receptor_design/receptor-resurfacing-pipeline/tests/rfdiffusion

  tests/rosetta_filtering/run_test_rosetta_filtering.slurm.sh
    → /hpc-home/jowillia/receptor_design/receptor-resurfacing-pipeline/tests/rosetta_filtering

  tests/proteinmpnn/run_test_proteinmpnn.slurm.sh
    → /hpc-home/jowillia/receptor_design/receptor-resurfacing-pipeline/tests/proteinmpnn

  tests/negative_steering/run_test_negative_steering.slurm.sh
    → /hpc-home/jowillia/receptor_design/receptor-resurfacing-pipeline/tests/negative_steering

  tests/orthogonal_metrics/run_test_orthogonal_metrics.slurm.sh
    → /hpc-home/jowillia/receptor_design/receptor-resurfacing-pipeline/tests/orthogonal_metrics

In each script, read the file first, then insert the `#SBATCH --chdir=<path>`
directive immediately after the existing `#SBATCH --mail-user` line. Do not
change any other part of any script.

STEP 2 — Delete stray SLURM log files from the repo root.

Find and delete all `slurm_*.out` and `slurm_*.err` files at the repo root
only (not inside subdirectories):

  find . -maxdepth 1 \( -name 'slurm_*.out' -o -name 'slurm_*.err' \)

Delete only files matching this pattern at maxdepth 1. Do not touch any
slurm log files inside tests/*/ directories.

ACCEPTANCE CRITERIA

  - Each of the five SLURM scripts contains exactly one `#SBATCH --chdir`
    line pointing at its test directory.
  - No `slurm_*.out` or `slurm_*.err` files remain at the repo root.
  - `grep -r "chdir" tests/*/run_test_*.slurm.sh` lists exactly five hits,
    one per script.
  - No other content in any script has changed.

DO NOT run git add, git commit, or git push. Staging and committing is the
human's responsibility.
```

---

## Prompt B — Restructure characterisation tests to per-module reference sets

Run after Prompt A is committed. Goal: all 184 hpc-tier characterisation tests
point at per-module reference sets, SEQUENCE_TRIO is replaced by SEQUENCE_COHORT
(10 sequences), the RFDiffusion fail-branch is covered by a new local_unit test,
and TRACEABILITY.md is consistent throughout. Acceptance: `pytest -m local_unit`
green, `pytest --collect-only -m hpc` clean, all hpc tests skip gracefully when
RECEPTOR_OUTPUT_ROOT is unset.

```
Continuing Phase 2. Prompts 1–3 produced the framework, comparators, and 184
hpc-tier characterisation tests. Those tests are currently tied to:

  - A single reference root: tests/full_test_run/example_output_files/
  - A fixed sequence trio: input_control_polyA, design_0_seq_0, design_13_seq_2

Both must change. Per-module reference sets now exist at
tests/<module>/example_output_files/, each containing the outputs of that
module's per-module test run on curated fixtures. This prompt restructures the
existing tests to use those per-module roots and the sequence names the new
fixtures actually produced.

Read these in order before writing anything:

  1. notes/inventory/14_phase_2_revision_per_module_tests.md — the active spec.
     Focus on §2.3 (reference set structure), §2.5 (fixture curation procedure),
     and §4.2 step 7 (what this restructure entails).
  2. notes/inventory/15_discovery_run_path_coverage.md — the fixture rationale.
     This tells you which sequence names exist in each per-module fixture and
     which pathway classes each exercises.
  3. notes/remediation_state.md — confirms current state and lists the five
     per-module reference sets now on disk.
  4. tests/characterization/conftest.py — the current fixture definitions you
     will be changing.
  5. Each tests/<module>/example_output_files/ directory — use `find` to confirm
     which files exist before updating any test. Do NOT pin files that are absent.

The five per-module reference roots are:

  tests/rfdiffusion/example_output_files/
  tests/rosetta_filtering/example_output_files/
  tests/proteinmpnn/example_output_files/
  tests/negative_steering/example_output_files/
  tests/orthogonal_metrics/example_output_files/

The corresponding per-module output roots (where per-module tests write their
results) are:

  tests/rfdiffusion/receptor_resurfacing_results/
  tests/rosetta_filtering/receptor_resurfacing_results/
  tests/proteinmpnn/receptor_resurfacing_results/
  tests/negative_steering/receptor_resurfacing_results/
  tests/orthogonal_metrics/receptor_resurfacing_results/

─────────────────────────────────────────────────────────────────────────────
PART 1 — conftest.py restructure
─────────────────────────────────────────────────────────────────────────────

Replace the single session-scoped `reference_root` and `output_root` fixtures
with per-stage variants. The new design:

  - A module-level mapping `STAGE_ROOTS` in conftest.py maps each stage name
    to its directory, relative to repo_root:

      STAGE_ROOTS = {
          "rfdiffusion":        "tests/rfdiffusion",
          "rosetta_filtering":  "tests/rosetta_filtering",
          "proteinmpnn":        "tests/proteinmpnn",
          "negative_steering":  "tests/negative_steering",
          "orthogonal_metrics": "tests/orthogonal_metrics",
      }

    For each stage:
      reference_root = <stage_dir>/example_output_files/
      output_root    = <stage_dir>/receptor_resurfacing_results/

  - A session-scoped fixture `stage_reference_root` that returns a callable
    accepting a stage name and returning the resolved Path. pytest.skip with
    a clear message if the directory doesn't exist.

  - A session-scoped fixture `stage_output_root` that returns a callable
    accepting a stage name and returning the resolved Path. pytest.skip with
    a clear message if the directory doesn't exist (output_root is only
    present after an HPC run).

  Both are factory fixtures (they return a callable, not a Path directly).
  Each test file binds its own stage:

      @pytest.fixture
      def ref(stage_reference_root):
          return stage_reference_root("negative_steering")

  Keep `reference_root` and `output_root` as deprecated aliases pointing at
  tests/full_test_run/example_output_files/ — this prevents import errors in
  any file not yet migrated. The aliases are deleted in a later commit when
  tests/full_test_run/example_output_files/ is retired.

  Keep all other fixtures unchanged: numeric_tolerance, path_normalizer,
  ssim_threshold, repo_root.

─────────────────────────────────────────────────────────────────────────────
PART 2 — Update each test file
─────────────────────────────────────────────────────────────────────────────

Eight test files need updating. For each, the changes are:

  (a) Replace reference_root / output_root fixture usage with the per-stage
      factory fixtures.
  (b) Update any SEQUENCE_TRIO or parametrize lists to match the sequence
      names that exist in the per-module fixture. Walk the per-module
      example_output_files/ to confirm exact names on disk.
  (c) Update any hardcoded subdirectory paths that assumed the old
      tests/full_test_run/example_output_files/ layout, if the per-module
      layout differs.
  (d) For any file pinned in the old reference set but absent from the
      per-module set, add pytest.skip() with a message explaining the absence.
      Do NOT delete the test — it remains as a placeholder for when the
      reference set is extended.

Work through the files in this order:

1. tests/characterization/test_rfdiffusion.py
   Stage name: "rfdiffusion"
   Walk tests/rfdiffusion/example_output_files/ to confirm which files are
   present. RFDiffusion is not parametrised by sequence; no SEQUENCE_TRIO
   update needed. The per-module test runs a single fixed input PDB and
   produces one set of outputs.

2. tests/characterization/test_rosetta_filtering.py
   Stage name: "rosetta_filtering"
   Walk tests/rosetta_filtering/example_output_files/.
   Rosetta filtering is parametrised by design PDB (design_0, design_1,
   design_14, design_59 from the curated fixture — 2 pass, 2 fail the sc
   threshold). Update the parametrize list to reflect these four design names.

3. tests/characterization/test_mpnn.py
   Stage name: "proteinmpnn"
   Walk tests/proteinmpnn/example_output_files/.
   The curated proteinmpnn fixture has 2 parent designs (design_0, design_28),
   each producing sequences. Update any design- or sequence-level parametrize
   lists to match what the walk reveals on disk.

4. tests/characterization/test_negsteer_cohort.py
   Stage name: "negative_steering"
   Walk tests/negative_steering/example_output_files/negative_steering/ for
   cohort-level files (cross_sequence_summary.csv, indices/, controls_inputs/).
   These are not parametrised by sequence; no SEQUENCE_TRIO update needed
   beyond confirming the files exist.

5. tests/characterization/test_negsteer_per_sequence.py
   Stage name: "negative_steering"
   This is the most changed file. The old SEQUENCE_TRIO was:
     ["input_control_polyA", "design_0_seq_0", "design_13_seq_2"]
   The new sequence set is all 10 sequences in the per-module fixture:
     ["design_28_seq_1", "design_62_seq_0", "design_0_seq_0", "design_42_seq_0",
      "design_3_seq_1",  "design_55_seq_1", "design_44_seq_1", "design_27_seq_0",
      "input_control_polyA", "input_control_scrambled"]
   Call this constant SEQUENCE_COHORT (replacing SEQUENCE_TRIO) so the name
   reflects the full curated cohort, not a trio.

   For each per-sequence file pinned in the existing test, walk the per-module
   reference set to confirm the file exists for each sequence before keeping
   the test. Cold-start sequences (design_28_seq_1) will lack cycle_0/steered/
   and some reversion outputs — use pytest.skip() where files are absent for
   specific sequences, consistent with the existing pattern.

6. tests/characterization/test_orthogonal_metrics.py
   Stage name: "orthogonal_metrics"
   Walk tests/orthogonal_metrics/example_output_files/.
   The per-module orthogonal_metrics test ran only NEGSTEER_INTERFACE_METRICS,
   EXTRACT_SURVIVOR_MANIFEST, and AF3_SETUP_DB (AF3/biophys/rosetta processes
   did not run — no GPU on the test node). Confirm which output files are
   actually present and pin only those. Tests for absent files get
   pytest.skip().

7. tests/characterization/test_preprocessing.py
   Stage name: None — preprocessing has no per-module test.
   Add a module-level skip on all tests in the file with message:
   "preprocessing has no per-module reference set; tests deferred until the
   full pipeline run reference is retired or a preprocessing per-module test
   is added."
   Do not restructure the test logic — just skip the whole file cleanly.

8. tests/characterization/test_plots.py
   Stage name: multiple — plots are emitted by several modules.
   Walk each per-module example_output_files/ for .png files. Update the
   parametrize list to cover the PNGs that exist across all five per-module
   reference sets. Each PNG test should record which stage it came from (add
   a `stage` parameter to the parametrize tuple alongside the filename).
   Plots absent from all per-module sets get pytest.skip().

─────────────────────────────────────────────────────────────────────────────
PART 3 — Small RFDiffusion fail-branch unit test
─────────────────────────────────────────────────────────────────────────────

The discovery run confirmed that the RFDiffusion fail-branch (passes_filter
== false) has never been observed in any real run and cannot be exercised by
the per-module fixture. This path is covered by a unit test instead.

Add one test to tests/characterization/test_rfdiffusion.py:

  test_rfdiffusion_filter_fail_branch_unit

  Mark: @pytest.mark.local_unit (NOT hpc — this test does not need HPC
  outputs).
  What it does: construct a minimal in-memory design dict with
  passes_filter=False (n_contact_pairs=0, frac_contacts_in_design=0.0,
  passes_filter=False), pass it through bin/rfdiffusion_filter.py's filter
  predicate logic directly (import the relevant function), and assert the
  design is excluded from the passing set. Do not invoke Nextflow or any
  external tool. If the filter logic is not easily importable (e.g. buried
  in a script with side effects), write the test against the logical
  condition itself with a clear comment explaining the indirection.

─────────────────────────────────────────────────────────────────────────────
PART 4 — TRACEABILITY.md update
─────────────────────────────────────────────────────────────────────────────

Update tests/characterization/TRACEABILITY.md:

  - For each restructured test (reference path changed, parametrize list
    changed), update the "Output pinned" column to reflect the new per-module
    reference path.
  - For each skipped test (file absent from per-module set), update the Notes
    column to "skipped — absent from per-module reference set".
  - For the new RFDiffusion fail-branch unit test, add a row with marker
    local_unit and strategy "import + assert".
  - For tests/characterization/test_preprocessing.py entries, update Notes
    to "file-level skip — no per-module reference set".

─────────────────────────────────────────────────────────────────────────────
ACCEPTANCE CRITERIA
─────────────────────────────────────────────────────────────────────────────

  - `pytest -m local_unit` passes with zero failures. The new
    test_rfdiffusion_filter_fail_branch_unit must appear and pass.
  - `pytest --collect-only -m hpc` succeeds without errors. Count reported
    should be >= the previous count (SEQUENCE_COHORT expands the parametrize
    matrix).
  - `pytest -m hpc` with RECEPTOR_OUTPUT_ROOT unset reports all hpc tests as
    skipped, not errored. Skip messages must be informative (missing env var
    or missing reference directory — not a bare "skipped").
  - No test file imports a hardcoded path to
    tests/full_test_run/example_output_files/ except via the deprecated alias
    fixtures in conftest.py.
  - tests/characterization/test_preprocessing.py is present but all its tests
    skip cleanly (no errors, no failures).
  - TRACEABILITY.md rows are consistent with the test files (no phantom rows,
    no missing rows).

DO NOT run git add, git commit, or git push. Staging and committing is the
human's responsibility.

In the final summary, report:
  - Total test function count and parametrised instance count before and after.
  - Which files were found absent from one or more per-module reference sets
    and were therefore skipped.
  - Any judgment calls on comparison strategy or path handling that the human
    should review.
  - Whether the RFDiffusion fail-branch function was directly importable or
    required an indirect approach.
```
