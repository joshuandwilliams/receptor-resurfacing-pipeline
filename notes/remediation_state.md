# Remediation State

A living document tracking where the codebase remediation effort currently stands. Read this at the start of every session; update it at the end of every session.

**Last updated:** 2026-05-02

---

## Current Phase

**Phase 2 (Establish Behavioral Tests) — in progress, with a strategic revision.**

Plan written, framework scaffolded, comparators tested, Nextflow caching fixed, reference set rebuilt, 184 hpc-tier characterization tests written, JAVA_HOME fixed, repo synced to HPC.

**A pivot has been adopted** before the first HPC round-trip: per-module tests will become the primary characterization safety net, with full pipeline runs reserved for milestone verification only. See `notes/inventory/14_phase_2_revision_per_module_tests.md` for the revision plan and rationale.

---

## Just Completed

**Phase 2 work to date:**

- `notes/inventory/11_phase_2_plan.md` — original Phase 2 plan (still in repo as historical record; superseded by the revision below).
- `notes/inventory/14_phase_2_revision_per_module_tests.md` — pivot to per-module-tests-as-safety-net. Active spec.
- `tests/characterization/` — pytest framework with conftest, fixtures, README. Includes ComparisonResult dataclass and four+ comparators (CSV-EXACT, CSV-STRUCT, CSV-EXACT-MODULO-PATHS, JSON-DEEP, JSON-MODULO-PATHS, PNG-PERCEPTUAL, PNG-EXISTS, TEXT-EXACT, plus path normalisation).
- `tests/characterization/helpers/tests/` — 67 unit tests, all passing under `pytest -m local_unit`.
- `tests/characterization/TRACEABILITY.md` — 184 rows mapping characterization tests to producers.
- `pyproject.toml` — test dependencies, pytest config, marker registry with `--strict-markers`.
- `containers/` — Singularity definition files moved into the repo with a README.
- Nextflow cache-busting — orthogonal-metrics pattern (`path X_script` input, content-hashed) propagated across modules; corresponding call-site updates in `main.nf` and `tests/*/test_*.nf`.
- `tests/full_test_run/example_output_files/` — supervisor-demo reference set (currently the canonical reference; will be retired once per-module references are in place per the revision plan).
- 184 hpc-tier characterization tests across 8 stage files, parametrized over the trio `input_control_polyA, design_0_seq_0, design_13_seq_2`.
- Pre-HPC-roundtrip audit (`audit_pre_hpc_roundtrip.md`); JAVA_HOME path corrected in 7 slurm launchers.
- `scripts/sync_to_hpc.sh` — sync script with documented excludes; `.gitignore` tightened for Word lockfiles and slurm log files.
- Repo synced to HPC at `/Volumes/HPC-Home/receptor_design/receptor-resurfacing-pipeline/`.

**Phase 1 deliverables** remain in `notes/inventory/` (`01_module_map.md` through `10_phase_1_synthesis.md`).

---

## Next Concrete Step

**Begin execution of the revised Phase 2 plan.**

Per `14_phase_2_revision_per_module_tests.md` §4.2:

1. **Decouple per-module tests from upstream chaining.** Each `tests/<module>/test_<module>.nf` should resolve its inputs from `tests/<module>/data/` as the canonical path, not from a previous module's results. The existing fallback logic in `test_negative_steering.nf` already supports this; other tests may need similar adjustment.

2. **Design and run the discovery run on HPC.** Params per revision §2.4: `rfdiffusion_n_designs=64`, `mpnn_seqs_per_design=2`, all 128 forward to negsteer, `negsteer_n_designs=4`, Branch B only. Estimated runtime 2–4 hours. Capture the full `results/` tree.

3. **Inventory the discovery run** to identify per-stage path coverage; pick fixture candidates.

(Steps 4–9 follow per the revision plan.)

The 184 existing characterization tests will need their `reference_root` resolution and parametrize lists updated once per-module reference sets are in place.

---

## Important Context Not Captured Elsewhere

These are decisions, framings, or observations that matter for future sessions but aren't formalized in the inventory or synthesis documents.

### Why we pivoted

The original Phase 2 plan assumed the supervisor-demo run (~6 hours) would be both the source of reference outputs and the verification target during Phase 3+ iteration. Three problems:

1. 6-hour iterations are unworkable across many Phase 3 commits.
2. Path coverage was retroactively-discovered, not designed-in.
3. Per-module tests were silently chained to each other, defeating their isolation.

The pivot uses per-module tests as primary safety net with curated input fixtures, and reserves full pipeline runs for milestone verification. See `14_phase_2_revision_per_module_tests.md` §1 for full rationale.

### Negative-steering path semantics (recorded for clarity)

- **Cold-start**: predicting where the MPNN sequence's effector lands. Can pass (correct placement) or fail (incorrect placement).
- **Steering**: introducing mutations. May produce contamination if mutations interact with the target.
- **Reversion**: rolling back contaminating mutations to test whether correct placement persists.

The path-coverage outcomes the discovery run aims to capture: cold-start fail; cold-start pass + clean steering; cold-start pass + contaminated + reversion succeeds; cold-start pass + contaminated + reversion fails. Plus the two controls (polyA, scrambled) which the test workflow generates from `params.input_pdb`.

### "Golden master" is the wrong framing

Phase 2's reference outputs are **current pipeline behavior**, not verified-correct outputs. Some outputs almost certainly contain undetected logical errors. Numbers can look realistic without being correct.

This means:

- Characterization tests assert *stability*, not *correctness*. Their job is to make changes visible, not to certify behavior as right.
- A test failing during refactoring does not necessarily indicate the refactor is wrong. It indicates behavior changed.
- Some refactoring may *fix* latent bugs, which will appear as test failures. Investigate, don't roll back automatically.

### File-by-file simplification belongs to Phase 4

The original plan included consolidation of duplicated functionality in Phase 3. A deeper "rewrite for clarity" pass is intentionally deferred to Phase 4 (deep modules), for two reasons:

- Phase 2 characterization tests need to exist first to make rewrite-for-clarity safe.
- A file-by-file deep simplification pass is more rigorous than a bird's-eye redundancy sweep across the whole codebase. It's also more naturally combined with the architectural restructuring of Phase 4.

Phase 3 will still consolidate obvious duplication that the inventory has already surfaced (e.g., the four identical blocks in `boltz2_iterate_steering.py`), but ambitious simplification waits.

### Dead-code deletion deferred for safety-net validation

`bin/sequence_registry.py` (vulture-confirmed unreachable) is intentionally **not** being deleted yet. Per the revision plan §4.2 step 9, it is reserved as the canonical safety-net validation commit — the smallest possible change to validate that the per-module tests + characterization framework correctly detect (or correctly don't detect) behavior changes. The empty `main` file at the repo root has already been removed.

### Specific known issues to track

- **`bin/boltz2_iterate_steering.py`** is the largest file (5,500+ lines), contains four identical code blocks that build the same `cmd` and call `subprocess.run` on `submit_script`, and references a `submit_boltz2_iterate_steering.sh` that does not exist on disk. The script has fallback handling for this missing file (`if not submit_script.exists(): print(WARN...)`), so the dead path doesn't break execution — but it represents an entire orphaned self-resubmission feature from a pre-Nextflow workflow.
- **`bin/reversion.py::harvest_reversion_results`** has cyclomatic complexity 74 — pathological. Top refactoring target.
- **`bin/cross_sequence_summary.py::aggregate`** has cyclomatic complexity 58. Second-worst offender.
- **`bin/rfdiffusion_plots.py`** contains several D-rated functions (CC 21–24). Plotting code, but worth attention.
- **Vulture's 8 high-confidence findings** (`08_vulture_high_confidence.txt`) are the safest dead-code candidates to remove first.

### HPC workflow constraint

The pipeline runs on HPC. The Mac is for development with Claude Code only — code cannot be tested there. Workflow:

1. Develop changes locally with Claude Code.
2. Sync via `./scripts/sync_to_hpc.sh` (the HPC home is mounted at `/Volumes/HPC-Home/`).
3. Run pipeline / tests on HPC.
4. Iterate.

Round-trips are slow. Plan thoroughly before each transfer; batch related changes; lean on `pytest -m local_unit` and static analysis (which run locally) for fast feedback during development.

### Tooling baseline

Local Mac has:

- Homebrew (installed during Phase 1 setup)
- Node.js + npm (via Homebrew)
- Claude Code (`@anthropic-ai/claude-code`, npm-global, runs as `claude`)
- Conda (miniforge base) with `ruff`, `vulture`, `radon` installed
- Conda env `receptor-tests` (Python 3.10) with the test framework's dependencies — created via `pyproject.toml`'s `[test]` extra
- Git, with project pushed to https://github.com/joshuandwilliams/receptor-resurfacing-pipeline
- HPC home mounted at `/Volumes/HPC-Home/` (read+write; treat as read-only when prompting Claude Code)

The HPC does not have git installed. Work happens on Mac; HPC receives synced code only.

### Branch and tag state

- Default branch: `main` (pristine baseline; do not modify).
- Active branch: `remediation` (all work happens here).
- Tags so far:
  - `baseline-pre-remediation`
  - `phase-1.1-inventory-complete`
  - `phase-1.2-glossary-complete`
  - `phase-1.3-static-analysis-complete`
  - `phase-1-complete`
  - `phase-2.1-plan-complete`
  - `phase-2.2-scaffolding-complete`
  - `phase-2.3-comparators-complete`
  - `phase-2.4-cache-busting-complete`
  - `phase-2.5-reference-rebuild-complete`

Tag at the completion of each phase or significant sub-step.

### Anticipated remaining Phase 2 tag sequence (under the revision)

- `phase-2.6-revision-plan-complete` — when the revision plan and updated state document are committed (this commit).
- `phase-2.7-discovery-run-complete` — after the discovery run executes and is curated.
- `phase-2.8-fixtures-and-references-complete` — after per-module fixtures and reference sets are in place.
- `phase-2.9-tests-restructured-complete` — after the 184 characterization tests are updated to the per-module structure.
- `phase-2-complete` — after the safety-net validation deletion commit produces zero diffs.

The earlier-anticipated `phase-2.6-hpc-tier-complete` is superseded — the equivalent verification now happens at `phase-2.9`.

---

## Verification Queue

Outputs or behaviors to scrutinize for correctness, but not investigated immediately. Revisited during or after refactoring.

> - **What:** `merge_orthogonal_metrics.py` semantic divergence between production and test versions.
> - **Why suspicious:** Finding A4 in `05_findings.md` flags this as a latent bug — production gates on AF3 presence; test demotes it to a flag-only column. Reference output reflects whichever version actually ran in the supervisor-demo cohort, and the Phase 2 test pin will lock that behavior — including the bug — until Phase 3.
> - **How to verify:** Compare the production `merge_orthogonal_metrics.py` against the test-tree copy. Inspect the cohort's AF3 column distribution. If the production version was used, expect this test to fail when Phase 3 fixes the divergence; update the reference at that point.

> - **What:** Suspected `--n-cycles 1` silent-skip-reversion bug from notes6.
> - **Why suspicious:** Synthesis §5 records the bug as not verified during Phase 1. The supervisor-demo ran with `--n-cycles 1`, so reference outputs may bake in skipped reversion behavior.
> - **How to verify:** Inspect a per-sequence run's `pathways.json` and reversion JSONs; cross-check whether reversion was actually attempted on contaminated sequences. If skipped where it should not have been, the reversion-related JSONs need regenerating from a `--n-cycles >= 2` run before Phase 3.

> - **What:** `scaffold_rmsd` field in `rfdiffusion_metrics.json` is actually motif RMSD.
> - **Why suspicious:** Glossary §F3 — field name is inverted from Baker-lab convention. Semantics likely fine; naming is misleading.
> - **How to verify:** Naming-only fix in Phase 3.3. Test will need to be updated at that point — reference file's field name changes, not its values.

> - **What:** Cache-busting residual gaps — sub-scripts loaded indirectly are not content-hashed by Nextflow.
> - **Why suspicious:** `bin/negative_steering_run_one.sh` loads sub-scripts via `--bin-dir` at runtime; `cross_sequence_summary.py` and the `NEGSTEER_CONTROLS` heredoc import helpers via `sys.path.insert`. Edits to these indirectly-loaded scripts will not invalidate the corresponding process cache.
> - **How to verify:** When testing whether a Phase 3 edit invalidates cache correctly, edit the orchestrator or a directly-tracked script to be safe. Closing this gap requires a new pattern (declaring all sub-scripts as `path` inputs, or staging the entire `bin/` directory) and is deferred.

> - **What:** Pre-existing `fastrelax_xml` mismatch in `tests/orthogonal_metrics/test_orthogonal_metrics.nf`.
> - **Why suspicious:** The test harness was calling `NEGSTEER_ROSETTA_METRICS(rosetta_input_ch)` without the `fastrelax_xml` input the production module has required since the initial baseline. Corrected in passing during the cache-busting commit.
> - **How to verify:** Confirm nothing else in the test tree relied on the old shape.

> - **What:** 184 existing hpc-tier characterization tests are tied to the supervisor-demo reference paths and the trio `input_control_polyA, design_0_seq_0, design_13_seq_2`.
> - **Why suspicious:** The revision plan replaces both the reference set location and the trio. The tests still encode comparator strategies and per-stage logic correctly, but their resolution paths and parametrize lists are stale.
> - **How to verify:** Per revision §4.2 step 7, restructure tests once per-module references are built. The framework code (conftest.py, fixtures) needs updating so `reference_root` and `output_root` are per-stage. Each test file then resolves via per-stage fixtures.

---

## Open Questions

Genuine uncertainties that may need resolution at some point.

1. **HADDOCK branch in the per-module strategy.** Branch A is out of Wave 1 scope. Should `tests/haddock/` get a per-module reference set under the new strategy at all, or stay deferred until Branch A is redesigned? Suggested: defer until the redesign.
2. **Disk-space cost of per-module reference sets.** Each `tests/<module>/example_output_files/` mirrors that module's run outputs. Worth measuring after the first per-module rebuild to see whether further trimming is needed.
3. **Stochasticity in per-module outputs.** Even with fixed inputs, some Boltz-2 prediction outputs may vary slightly run-to-run if seeds aren't fully pinned. The first per-module round-trip will reveal which outputs need the comparator strategy demoting.

---

## How to Use This Document

**At the start of a session:**

1. Read this document.
2. Read `notes/codebase_remediation_plan.md` if it's been a while.
3. Read `notes/inventory/10_phase_1_synthesis.md` for the substantive Phase 1 findings.
4. Read `notes/inventory/14_phase_2_revision_per_module_tests.md` for the active Phase 2 spec.
5. Begin work on the "Next Concrete Step" listed above.

**During a session:**

- If a decision is made, note it.
- If a doubt or suspicion arises, add it to the Verification Queue.
- If a question can't be resolved now, add it to Open Questions.

**At the end of a session:**

- Update "Just Completed" with what was done.
- Update "Next Concrete Step" with what should happen next.
- Update "Last Updated" date.
- Commit and push.

This document is the canonical source of "where am I in this work." If it disagrees with another document, this one wins (and the other should be updated).
