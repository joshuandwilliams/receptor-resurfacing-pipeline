# Remediation State

A living document tracking where the codebase remediation effort currently stands. Read this at the start of every session; update it at the end of every session.

**Last updated:** 2026-05-01

---

## Current Phase

**Phase 2 (Establish Behavioral Tests) — in progress.**

Plan written, framework scaffolded, comparators tested, Nextflow caching fixed. Remaining: reference set gap-fill, hpc-tier characterization tests, suspicions pass, safety-net validation.

---

## Just Completed

**Phase 2 work to date:**

- `notes/inventory/11_phase_2_plan.md` — characterization test suite plan (two-wave structure, comparator strategies, manual reference-update procedure)
- `tests/characterization/` — pytest framework with conftest, fixtures, README
- `tests/characterization/helpers/` — four comparators (CSV, JSON, PNG, path normalisation) + `ComparisonResult` dataclass
- `tests/characterization/helpers/tests/` — 50+ unit tests, all passing under `pytest -m local_unit`
- `tests/characterization/TRACEABILITY.md` — schema for test-to-producer mapping (zero rows; populated when hpc-tier tests are written)
- `pyproject.toml` — test dependencies, pytest config, marker registry with `--strict-markers`
- `containers/` — Singularity definition files moved into the repo with a README
- Nextflow cache-busting — orthogonal-metrics pattern (`path X_script` input, content-hashed) propagated to all modules that invoke `bin/*.py` scripts; corresponding call-site updates in `main.nf` and `tests/*/test_*.nf`

**Phase 1 deliverables** remain in `notes/inventory/` (`01_module_map.md` through `10_phase_1_synthesis.md`).

---

## Next Concrete Step

**Reference set gap-fill on HPC.**

Locate and copy missing files into `tests/full_test_run/example_output_files/` per the ⚠️ checklist in plan §2:

- `rosetta_filter_summary.csv`
- Per-design `mpnn_results.json`; full `top_fastas/` and `af2_fastas/` directories
- All 34 per-sequence subdirectories under `negative_steering/` (32 designs + 2 controls), each containing the seven expected files
- A populated-reversion per-sequence example (currently all reversion JSONs on disk are empty)
- The two control sequence rows in `cross_sequence_summary.csv`
- `cross_sequence_summary_with_interface_metrics.csv`
- `merged_orthogonal_metrics.csv` / `survivors_with_orthogonal_metrics.csv` (confirm exact filenames)

For each file copied in, update `tests/full_test_run/example_output_files/README.md` with a one-line note. Genuine gaps (files that don't exist on HPC) get noted in the Verification Queue.

This is a manual HPC task — no Claude Code prompt. Required before Prompt 3 (hpc-tier characterization tests) can be drafted.

---

## Important Context Not Captured Elsewhere

These are decisions, framings, or observations that matter for future sessions but aren't formalized in the inventory or synthesis documents.

### "Golden master" is the wrong framing

Phase 2's reference outputs are **current pipeline behavior**, not verified-correct outputs. Some outputs — particularly large summary CSVs — almost certainly contain undetected logical errors. Numbers can look realistic without being correct.

This means:

- Characterization tests assert *stability*, not *correctness*. Their job is to make changes visible, not to certify behavior as right.
- A test failing during refactoring does not necessarily indicate the refactor is wrong. It indicates behavior changed. The change must then be evaluated on its merits.
- Some refactoring may *fix* latent bugs, which will appear as test failures. These are good outcomes that need investigation rather than rollback.

The README in `example_output_files/` carries this caveat. It should be re-read at the start of Phase 2 work.

### File-by-file simplification belongs to Phase 4

The original plan included consolidation of duplicated functionality in Phase 3. A deeper "could this be rewritten more cleanly while preserving behavior" pass is intentionally deferred to Phase 4 (deep modules), for two reasons:

- Phase 2 characterization tests need to exist first to make rewrite-for-clarity safe.
- A file-by-file deep simplification pass is more rigorous than a bird's-eye redundancy sweep across the whole codebase. It's also more naturally combined with the architectural restructuring of Phase 4.

Phase 3 will still consolidate obvious duplication that the inventory has already surfaced (e.g., the four identical blocks in `boltz2_iterate_steering.py`), but ambitious simplification waits.

### Dead-code deletion deferred for safety-net validation

`bin/sequence_registry.py` (vulture-confirmed unreachable) and the empty `main` file at the repo root are intentionally **not** being deleted yet. Per plan §6.3, they are reserved as the canonical first Phase 3 commit — the smallest possible change to validate the safety net. Deleting them now would use up the cleanest validation case for nothing.

### Specific known issues to track

- **`bin/boltz2_iterate_steering.py`** is the largest file (5,500+ lines), contains four identical code blocks that build the same `cmd` and call `subprocess.run` on `submit_script`, and references a `submit_boltz2_iterate_steering.sh` that does not exist on disk. The script has fallback handling for this missing file (`if not submit_script.exists(): print(WARN...)`), so the dead path doesn't break execution — but it represents an entire orphaned self-resubmission feature from a pre-Nextflow workflow.
- **`bin/reversion.py::harvest_reversion_results`** has cyclomatic complexity 74 — pathological. Top refactoring target.
- **`bin/cross_sequence_summary.py::aggregate`** has cyclomatic complexity 58. Second-worst offender.
- **`bin/rfdiffusion_plots.py`** contains several D-rated functions (CC 21–24). Plotting code, but worth attention.
- **Vulture's 8 high-confidence findings** (`08_vulture_high_confidence.txt`) are the safest dead-code candidates to remove first.

### HPC workflow constraint

The pipeline runs on an airgapped HPC. The Mac is for development with Claude Code only — code cannot be tested there. Workflow:

1. Develop changes locally with Claude Code.
2. Sync modified files to HPC (manual transfer; airgap means no rsync over network).
3. Run pipeline / tests on HPC.
4. Sync results back to Mac.
5. Iterate.

Round-trips are slow. Plan thoroughly before each transfer; batch related changes; lean on static analysis (which runs locally) for early-phase work.

### Tooling baseline

Local Mac has:

- Homebrew (installed during Phase 1 setup)
- Node.js + npm (via Homebrew)
- Claude Code (`@anthropic-ai/claude-code`, npm-global, runs as `claude`)
- Conda (miniforge base) with `ruff`, `vulture`, `radon` installed
- Conda env `receptor-tests` (Python 3.10) with the test framework's dependencies — created via `pyproject.toml`'s `[test]` extra
- Git, with project pushed to https://github.com/joshuandwilliams/receptor-resurfacing-pipeline

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

Tag at the completion of each phase or significant sub-step.

### Anticipated Phase 2 tag sequence

- `phase-2.5-gap-fill-complete` — after the HPC reference-set gap-fill
- `phase-2.6-hpc-tier-complete` — after Prompt 3 (hpc-tier characterization tests, green on a fresh HPC run)
- `phase-2.7-suspicions-complete` — after Prompt 4 (suspicion-finding pass producing `SUSPICIONS.md`)
- `phase-2-complete` — after the dead-code deletion (`bin/sequence_registry.py`, empty `main` file) produces zero diffs in the test suite, validating the safety net

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

---

## Open Questions

Genuine uncertainties that may need resolution at some point.

*(None currently.)*

---

## How to Use This Document

**At the start of a session:**

1. Read this document.
2. Read `notes/codebase_remediation_plan.md` if it's been a while.
3. Read `notes/inventory/10_phase_1_synthesis.md` for the substantive findings.
4. Read `notes/inventory/11_phase_2_plan.md` for the Phase 2 design.
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
