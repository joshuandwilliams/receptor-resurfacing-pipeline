# Remediation State

A living document tracking where the codebase remediation effort currently stands. Read this at the start of every session; update it at the end of every session.

**Last updated:** 2026-04-30

---

## Current Phase

**Phase 1 — complete.**
**Phase 2 (Establish Behavioral Tests) — about to begin.**

The boundary between Phase 1 and Phase 2 marks the transition from diagnosis (read-only analysis) to active intervention. From this point forward, code changes are possible, though Phase 2's changes are limited to *adding* characterization tests rather than modifying the existing pipeline.

---

## Just Completed

**Phase 1: full diagnostic pass.** All deliverables are in `notes/inventory/`:

- `01_module_map.md` — every code file with one-line summaries
- `02_function_inventory.md` — Python functions and Nextflow processes catalogued
- `03_dependency_graph.md` — import and invocation relationships
- `04_functional_categorization.md` — modules grouped by purpose (the codebase's "design concept" made explicit)
- `05_findings.md` — anomalies and issues spotted during the inventory
- `06_ubiquitous_language.md` — canonical glossary of domain terms (interactively reviewed and curated)
- `07_ruff_report.txt`, `07_ruff_summary.txt` — linting (~194 findings)
- `08_vulture_report.txt`, `08_vulture_high_confidence.txt` — dead code candidates (66 medium-confidence, 8 high-confidence)
- `09_radon_complexity.txt`, `09_radon_maintainability.txt`, `09_radon_loc.txt` — complexity and size metrics
- `10_phase_1_synthesis.md` — synthesis briefing summarizing findings, priorities, and refinements to the plan

**Reference output sampling.** Hand-picked example output files from a full pipeline test run have been copied to `tests/full_test_run/example_output_files/`. A README in that directory describes each file's structure and pipeline-stage origin (committed; the data files themselves are gitignored).

---

## Next Concrete Step

**Begin Phase 2: characterization test setup.**

Specifically, the next session should:

1. Identify a tractable initial set of pipeline outputs (likely 5–15) to use as references for characterization tests. Use the synthesis document and `example_output_files/README.md` to choose representative outputs across the major pipeline stages.
2. For each chosen output, decide on the appropriate comparison strategy:
   - Exact equality (for deterministic, low-volume outputs like config dumps or small summary files).
   - Structural assertion (for tabular outputs where row counts, column names, types matter but exact values may vary).
   - Tolerance-based numerical comparison (for floating-point metrics where bitwise equality is too brittle).
   - Partial / property-based assertions (for large outputs where full validation is impossible — e.g., asserting summary statistics rather than exact contents).
3. Sketch the characterization test framework: where the test scripts live, how they are invoked, what the comparison harness looks like. Tests will need to run on the HPC (since that's where the pipeline runs), but the comparison code itself can be lightweight Python that runs anywhere.
4. Produce a *plan document* before writing tests, so the approach can be reviewed before significant work happens.

The first concrete deliverable of Phase 2 is therefore a planning document, not yet test code.

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

Tag at the completion of each phase or significant sub-step.

---

## Verification Queue

Outputs or behaviors I want to scrutinize for correctness, but don't have time/clarity to investigate immediately. To be revisited during or after refactoring.

*(Empty for now — populate as suspicious outputs are noticed during Phase 2 setup or later phases.)*

Suggested template for entries:

> - **What:** Brief description of the output or behavior.
> - **Why suspicious:** What gave me pause.
> - **How to verify:** Reference data, paper, manual calculation, or alternative implementation that could be used to check.

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
4. Begin work on the "Next Concrete Step" listed above.

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
