# Remediation State

A living document tracking where the codebase remediation effort currently stands. Read this at the start of every session; update it at the end of every session.

**Last updated:** 2026-05-03 (Phase 2 complete)

---

## Current Phase

**Phase 2 (Establish Behavioral Tests) — COMPLETE.**

Tag: `phase-2-complete`.

The full per-module characterisation safety net is in place: five per-module Nextflow tests with curated fixtures, five per-module reference sets, 272 passing hpc-tier pytest tests (298 parametrized instances, 26 expected skips, 0 failures), and a confirmed-dead `bin/sequence_registry.py` deleted as the safety-net validation commit.

---

## Just Completed (this session)

- **Phase 2.7 closed.** All five per-module tests ran successfully on HPC against curated fixtures (`Success: true`; negsteer 57m, others <10 min each).
- **Spot-check passed.** Negsteer cohort summary confirmed: 10 rows (8 steered + 2 controls), all 8 curated pathway classes present with correct `aggregated_verdict` and `cross_tier`. Cold-start sentinel for `design_28_seq_1` confirmed (single `initial` row, `no_reversion`).
- **Subtractive rebuild.** Five `tests/<module>/example_output_files/` tarballs built on HPC (rfdiffusion 1.6M, rosetta_filtering 1.4M, proteinmpnn 1.5M, negative_steering 3.5M, orthogonal_metrics 1.4M), pulled to Mac, unpacked, committed per-module. Tags: `phase-2.7-discovery-run-complete`, `phase-2.8-fixtures-and-references-complete`.
- **Phase 2.9 complete.** Two Claude Code prompts:
  - Prompt A: added `#SBATCH --chdir` to five per-module SLURM scripts; deleted stray `slurm_*.out/.err` from repo root.
  - Prompt B: restructured 184 → 298 parametrized hpc-tier tests to per-module reference roots; expanded negsteer `SEQUENCE_TRIO` (3) → `SEQUENCE_COHORT` (10); added `test_rfdiffusion_filter_fail_branch_unit` (local_unit, indirect approach with clear producer comment); updated `TRACEABILITY.md` throughout; `test_preprocessing.py` file-level skipped (no per-module fixture). Tag: `phase-2.9-tests-restructured-complete`.
- **`pytest_runner.def` container added.** Lightweight AlmaLinux 9 / Python 3.11 container for running pytest (pytest, pandas, numpy, Pillow, scikit-image only). Built at `/hpc-home/jowillia/singularity/pytest/pytest_runner.img`. Added to `containers/README.md`.
- **`run_pytest.slurm.sh` added** at `tests/characterization/`. Runs `pytest -m hpc` inside the pytest_runner container with correct bind mount.
- **Characterisation suite green on HPC.** `272 passed, 26 skipped, 0 failed in 20.75s`. All skips expected and documented (effector_template.cif excluded from subtractive rebuild; cold-start per-sequence files absent for design_28_seq_1; orthogonal plots/survivors absent — AF3/biophys/rosetta didn't run on per-module fixture; preprocessing deferred).
- **`bin/sequence_registry.py` deleted.** Safety-net validation: vulture confirmed all five functions dead (no other file imported them). Pytest green after deletion confirms zero behavioural change.
- **`tests/full_test_run/example_output_files/` retired.** Deleted from disk; tracked files (`README.md`) removed via `git rm`. Per doc 14 §5 order of operations completed.
- **`#SBATCH --chdir` added to all 14 SLURM scripts.** Comprehensive fix across all `*.slurm.sh` files in the repo (8 were missing it, including `run_pipeline.slurm.sh` and all plotting helpers).
- **Prompt document `notes/inventory/16_prompt_phase_2_9.md` committed** (two-prompt structure: Prompt A = SLURM housekeeping, Prompt B = test restructure).

---

## Phase 2 Deliverables (cumulative)

- `notes/inventory/11_phase_2_plan.md` — original plan (historical record).
- `notes/inventory/14_phase_2_revision_per_module_tests.md` — active spec (supersedes §2, §6, §8 of original).
- `notes/inventory/15_discovery_run_path_coverage.md` — fixture rationale, pathway taxonomy, 8 observable Classes.
- `notes/inventory/16_prompt_phase_2_9.md` — Claude Code prompts for Phase 2.9.
- `tests/characterization/` — pytest framework, conftest (per-stage factory fixtures), README, TRACEABILITY.md.
- `tests/characterization/helpers/` — six comparators (CSV-EXACT, CSV-STRUCT, CSV-EXACT-MODULO-PATHS, JSON-DEEP, JSON-MODULO-PATHS, PNG-PERCEPTUAL/EXISTS, TEXT-EXACT) + `ComparisonResult` dataclass.
- 68 `local_unit` tests passing (67 comparator unit tests + 1 RFDiffusion fail-branch unit test).
- 272 `hpc` tests passing (298 parametrized instances, 26 expected skips) across 8 stage files, parametrized over `SEQUENCE_COHORT` (10 sequences covering Classes 1–7 + 2 controls).
- Five per-module fixture sets in `tests/<module>/data/`.
- Five per-module reference sets in `tests/<module>/example_output_files/`.
- `containers/pytest_runner.def` + built image at `/hpc-home/jowillia/singularity/pytest/pytest_runner.img`.
- `tests/characterization/run_pytest.slurm.sh`.
- `scripts/sync_to_hpc.sh`, `scripts/rewrite_fixture_paths.py`.
- Nextflow cache-busting (`path X_script` content-hashed).
- `#SBATCH --chdir` on all 14 SLURM scripts.

**Phase 1 deliverables** remain in `notes/inventory/` (`01_module_map.md` through `10_phase_1_synthesis.md`).

---

## Next Concrete Steps — Phase 3

Phase 3 is **targeted simplification of the highest-complexity producers**, using the Phase 2 safety net to verify each change produces zero characterisation diffs.

The priority order follows the complexity findings from Phase 1:

1. **`bin/boltz2_iterate_steering.py`** — 5,500+ lines, four identical code blocks, references a non-existent submit script. Top refactoring target. Start here.
2. **`bin/reversion.py::harvest_reversion_results`** — cyclomatic complexity 74. Extract sub-functions.
3. **`bin/cross_sequence_summary.py::aggregate`** — CC 58. Second-worst.
4. **`bin/rfdiffusion_plots.py`** — several D-rated functions (CC 21–24).
5. **Docstring fix:** `bin/cross_sequence_summary.py` tier-rule docstring (~lines 17–23) claims `n_seeds_pose_holds == n_seeds`; implementation uses `n_pass = n_seeds_pose_holds + n_seeds_clean_steered`. One-line fix.
6. **Vulture high-confidence findings** (`08_vulture_high_confidence.txt`) — 8 items, safest dead-code candidates now that `sequence_registry.py` (the reserved one) is gone.

**Per-commit workflow during Phase 3:**

- Run the four cheap per-module tests (rfdiffusion, rosetta_filtering, proteinmpnn, orthogonal_metrics) after each substantive change — total <10 min.
- Run negsteer per-module test only for changes that touch negsteer producers — ~57 min.
- Run `pytest -m hpc` after each module sweep to confirm zero diffs.
- Reserve full pipeline runs for end-of-phase milestones.

---

## Important Context Not Captured Elsewhere

### "Golden master" framing

Phase 2 reference outputs are *current pipeline behavior*, not verified-correct outputs. Tests assert *stability*, not *correctness*. A test failing during Phase 3 refactoring may indicate a legitimate behavior change or a fix to a latent bug — investigate, don't roll back automatically.

### Negative-steering path semantics

The full pathway taxonomy lives in `notes/inventory/15_discovery_run_path_coverage.md` §Negative steering — four orthogonal axes (cold-start outcome / per-seed verdict / aggregated verdict / cohort tier) and the eight observable per-MPNN-sequence outcome Classes. That document is the authority on negsteer semantics for fixture purposes.

### Class 8 (`new_contamination`) coverage gap

The discovery run produced no `new_contamination` examples at either per-seed or aggregated resolution. This path is **not currently covered** by any fixture or unit test. Deferred from Phase 2 (the safety-net validation was green without it). Options for Phase 3: (a) manufacture an input that triggers `new_contamination` aggregation via a small targeted negsteer run; (b) write a Python unit test against `classify_reversion_verdict` + `_classify_aggregated_verdict` with a hand-constructed reversion-result blob. Option (b) is cheaper; option (a) is more representative.

### Known expected skips in pytest

These 26 skips are permanent fixtures of the current reference set, not regressions:
- `cycle_0/effector_template.cif` — `.cif` excluded from subtractive rebuild across all 10 sequences.
- Cold-start per-sequence files — `design_28_seq_1` lacks `cycle_statistics.csv`, `cycle_0/passing.json`, `cycle_0/summary.txt`, `steered_results_aggregate.csv`, `steered_results.csv`, `{true,wrong}_interface_residues.txt` (correct — cold-start path never runs steering).
- `inputs/receptor.fasta` — both controls use `control_*_receptor.fasta` variants instead.
- `cross_sequence_summary_with_interface_metrics.csv` — pinned in orthogonal_metrics, not negsteer_cohort.
- `survivors_with_orthogonal_metrics.csv` — AF3/biophys/rosetta did not run on per-module fixture.
- 4 `orthogonal_*.png` plots — same reason.
- 2 preprocessing tests — no per-module preprocessing fixture.

### HPC workflow constraint

Mac is authoritative. HPC has no git. Workflow: edit on Mac → `./scripts/sync_to_hpc.sh` → run on HPC → iterate. Round-trips are slow.

### Tooling

- Mac: Homebrew, Node.js, Claude Code, miniforge conda, `receptor-tests` env (Python 3.10), Nextflow v26.04.0.0.
- HPC: Nextflow via Singularity (legacy parser default). JAVA_HOME: `/hpc-home/jowillia/singularity/jdk-17.0.2`.
- HPC home mounted at `/Volumes/HPC-Home/` (treat as read-only when prompting Claude Code).
- Containers at `/hpc-home/jowillia/singularity/`.

### Branch and tag state

- Default branch: `main` (pristine baseline).
- Active branch: `remediation`.
- Tags: `baseline-pre-remediation`, `phase-1.1` through `phase-1-complete`, `phase-2.1-plan-complete` through `phase-2-complete`.

### Phase 3+ work currently deferred

- File-by-file simplification (the main Phase 3 job).
- Strict-syntax migration of `main.nf` and test workflows (deferred until HPC's Nextflow forces it).
- Closing cache-busting residual gaps for indirectly-loaded sub-scripts (`bin/negative_steering_run_one.sh`, `cross_sequence_summary.py` via `--bin-dir` / `sys.path.insert`).
- Renaming `scaffold_rmsd` field (semantically motif RMSD — Baker-lab convention inverted).
- Investigating `merge_orthogonal_metrics.py` test-vs-production divergence.
- Investigating suspected `--n-cycles 1` silent-skip-reversion bug.
- Class 8 (`new_contamination`) coverage — see above.
- HADDOCK branch (Branch A) per-module test — deferred until Branch A is redesigned.

---

## Verification Queue

Outputs or behaviors to scrutinize for correctness, but not investigated immediately. Revisited during or after refactoring.

> - **What:** `merge_orthogonal_metrics.py` semantic divergence between production and test versions.
> - **Why suspicious:** Finding A4 in `05_findings.md`. Reference output may bake in this divergence.
> - **How to verify:** Compare scripts and inspect AF3 column distribution in the discovery run's outputs.

> - **What:** Suspected `--n-cycles 1` silent-skip-reversion bug from notes6.
> - **Why suspicious:** Synthesis §5 records this as not verified during Phase 1.
> - **How to verify:** Inspect `pathways.json` and reversion JSONs for cohort sequences; check whether reversion was actually attempted on contaminated cases.

> - **What:** `scaffold_rmsd` field in `rfdiffusion_metrics.json` is actually motif RMSD.
> - **Why suspicious:** Glossary §F3 — field name inverted from Baker-lab convention.
> - **How to verify:** Naming-only fix in Phase 3.

> - **What:** Cache-busting residual gaps — `bin/negative_steering_run_one.sh` and `cross_sequence_summary.py` use `--bin-dir` / `sys.path.insert` for sub-scripts; not content-hashed by Nextflow.
> - **Why suspicious:** Edits to indirectly-loaded scripts won't invalidate Nextflow cache.
> - **How to verify:** During Phase 3, prefer editing directly-tracked scripts; closing the gap requires declaring all sub-scripts as `path` inputs.

> - **What:** Legacy Nextflow syntax (`workflow.onComplete {}` handlers, top-level `if` in `main.nf:43`). `NXF_SYNTAX_PARSER=v1` pinned in SLURM wrappers as stopgap.
> - **Why suspicious:** When HPC's Nextflow eventually upgrades, v1 parser removal will break compilation.
> - **How to verify:** Migrate event handlers to `nextflow.config`; lift top-level `if` into workflow body. Test against strict-default Nextflow on Mac.

> - **What:** Other potential `set -euo pipefail` + piped command sites with the SIGPIPE-on-141 bug.
> - **Why suspicious:** Pattern is easy to miss; failures are intermittent.
> - **How to verify:** Audit all `.nf` process bodies combining `pipefail` with piped commands.

> - **What:** `contact_cutoff` parameter documented under HADDOCK section but consumed by `rfdiffusion_filter.py` regardless of mode.
> - **Why suspicious:** Misleading section caused discovery run's first attempt to use `contact_cutoff: 1.0`, breaking Cα contact detection.
> - **How to verify:** Phase 3 — rename to `rfdiff_contact_cutoff`, move to its own section in `params_example.yml`.

> - **What:** Tier rule docstring at `bin/cross_sequence_summary.py` lines ~17–23 is stale.
> - **Why suspicious:** Claims `n_seeds_pose_holds == n_seeds`; implementation uses `n_pass = n_seeds_pose_holds + n_seeds_clean_steered`.
> - **How to verify:** One-line docstring fix in Phase 3.

> - **What:** Production survivor manifest does not gate on `cross_tier`.
> - **Why suspicious:** AF3 + biophys + rosetta computed for tier-none representatives. Meaningful GPU time in production.
> - **How to verify:** Product decision — whether diagnostic completeness is wanted or a tier gate would save time without losing signal.

> - **What:** `interface_plddt_median` column blank in `survivors_with_orthogonal_metrics.csv` but `interface_plddt_too_low` gate fired on 70 sequences.
> - **Why suspicious:** Gate may read a differently-named column than what gets emitted in the cross_summary join.
> - **How to verify:** Trace `merge_orthogonal_metrics.py` — which input column the gate consumes vs. which output column it writes.

> - **What:** `cross_sequence_summary_with_interface_metrics.csv` columns `irmsd`, `fnat`, `dockq` may be silently wrong in production runs predating the `negsteer_interface_metrics.nf` chain-param fix.
> - **Why suspicious:** Same shape as the bug pipeline_notes10 fixed for the other four orthogonal-metrics modules — `--effector-chain ${params.effector_chain}` (input-PDB convention, "C") instead of `--effector-chain ${params.rfdiff_output_effector_chain}` (prediction-PDB convention, "B").
> - **How to verify:** After next per-module orthogonal_metrics sweep, spot-check survivors against any prior production CSV. If materially different, alert users of historical interface-metrics data.

> - **What:** Pipeline output verbosity — work/ and results/ trees contain many thousands of tiny files per run.
> - **Why suspicious:** Deleting a full run takes minutes-to-hours on HPC filesystem. The fixture tarball exclude list (*.npz, *.a3m, msa/, predictions/, contamination_scratch/) is itself evidence of the problem.
> - **How to verify:** Phase 3 audit — identify files emitted but never consumed downstream. Likely fixes: scratch-only outputs into work/; coalesced per-seed outputs; explicit deletion at end of negsteer cycles.

> - **What:** Per-module test cost asymmetry — negsteer ~57 min vs. <10 min for other four.
> - **Why suspicious:** Design constraint, not a bug. Current params (negsteer_n_designs=4, num_seeds=3) are at the floor for fixture path coverage.
> - **How to verify:** Treat four cheap tests as the default per-commit safety net during Phase 3; reserve negsteer for substantive changes. Class 5's 57-minute outlier (design_3_seq_1, 12 contaminated mutations) is the cost-driving case.

---

## Open Questions

1. **Class 8 (`new_contamination`) coverage.** Manufactured fixture on HPC vs. Python unit test? See above for options.
2. **HADDOCK branch per-module test.** Deferred until Branch A is redesigned.
3. **Stochasticity in per-module outputs.** First round-trip showed no issues, but future Boltz-2 version changes may require comparator strategy demotions (e.g. JSON-DEEP → JSON-STRUCT for prediction outputs).

---

## How to Use This Document

**At the start of a session:**

1. Read this document.
2. Read `notes/codebase_remediation_plan.md` if it's been a while.
3. Read `notes/inventory/10_phase_1_synthesis.md` for Phase 1 findings.
4. Begin work on the first item under "Next Concrete Steps — Phase 3".

**During a session:**

- If a decision is made, note it.
- If a doubt or suspicion arises, add it to the Verification Queue.
- If a question can't be resolved now, add it to Open Questions.

**At the end of a session:**

- Update "Just Completed" with what was done.
- Update "Next Concrete Steps".
- Update "Last Updated" date.
- Commit and push.

This document is the canonical "where am I" source. If it disagrees with another document, this one wins (and the other should be updated).