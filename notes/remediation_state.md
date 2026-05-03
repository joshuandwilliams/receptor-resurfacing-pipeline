# Remediation State

A living document tracking where the codebase remediation effort currently stands. Read this at the start of every session; update it at the end of every session.

**Last updated:** 2026-05-02 (late evening, post-curation)

---

## Current Phase

**Phase 2 (Establish Behavioral Tests) — in progress, with a strategic revision in execution.**

Plan written, framework scaffolded, comparators tested, Nextflow caching fixed, supervisor-demo reference set built (now being retired), 184 hpc-tier characterization tests written, JAVA_HOME fixed, repo synced to HPC, per-module test decoupling complete, discovery run executed and curated, per-module fixtures staged into each `tests/<module>/data/`.

The pivot to per-module tests as the primary safety net is described in `notes/inventory/14_phase_2_revision_per_module_tests.md` — that's the active spec.

---

## Just Completed (since last update)

- **Discovery run completed.** `params_full_test.yml` ran to completion on HPC after the SIGPIPE fix, producing 122 cross-sequence rows (120 steered + 2 controls) covering 7 of 8 observable per-MPNN-sequence outcome classes.
- **Path-coverage analysis.** `notes/inventory/15_discovery_run_path_coverage.md` derives the negsteer pathway taxonomy from the producer code's vocabulary (4 orthogonal axes: cold-start outcome / per-seed verdict / aggregated verdict / cohort tier), tabulates the 8 observable classes, identifies discriminating signals per class, and proposes minimum-coverage fixture candidates per stage.
- **Verification round on HPC.** All verification commands from report 15 run on HPC; results applied back into the report. Resolved Open Q3 (design_44_seq_1 tier-A anomaly: representative `n_pass = n_seeds_pose_holds + n_seeds_clean_steered = 1 + 2 = 3`). Surfaced and corrected: MPNN path was wrong (`top_metadata.csv` lives in `sequences/`, not `mpnn/`); `singleton` verdict is systematic (one per sequence, not "n/a here"); Class 3 splits into 3a (63 sequences, all-zeros) and 3b (6 sequences, mixed pose_collapses).
- **Class 8 (`new_contamination`) confirmed unreachable from this discovery run** at both per-seed and aggregated resolution. Will be addressed as a standalone task before Phase 3.
- **RFDiffusion fail branch confirmed not exercisable** by any real run, ever. Defensive code path; will be covered by a small unit test added during Phase 2.9, not by a fixture.
- **Orthogonal-metrics 0/122 pass rate downgraded** from "coverage gap" to "expected behavior" — the af3_nomsa_ra_eff gate is a warning-style filter; pass branch exercisable via param override in the per-module test.
- **Per-module fixtures curated.** Five `tests/<module>/data/` directories populated:
  - `tests/rfdiffusion/data/` — no-op (already correctly set up from decoupling).
  - `tests/rosetta_filtering/data/rfdiffusion_split/` — 4 split PDBs (design_0, 1, 14, 59) covering sc_threshold both sides.
  - `tests/proteinmpnn/data/rosetta_passing/` — 2 parent PDBs (design_0, design_28); `input_complex.pdb` populated from rfdiffusion data (was 0-byte placeholder).
  - `tests/negative_steering/data/` — 8 FASTAs spanning Classes 1, 2, 3a, 4, 5, 6, 7 + 8 deduped parent design PDBs + `rfdiffusion_metrics.json` + `input_complex.pdb`.
  - `tests/orthogonal_metrics/data/negsteer_run/` — 2 per-sequence workdirs (design_28_seq_1, design_62_seq_0) + trimmed `cross_sequence_summary.csv`. HPC absolute paths rewritten to local Mac/repo paths via `scripts/rewrite_fixture_paths.py`.
- **Path-rewriter script added** at `scripts/rewrite_fixture_paths.py` — stdlib-only, idempotent, re-runnable when fixtures are re-pulled.
- **`sync_to_hpc.sh` continuation-line fix** committed (the missing `\` after `nxf_home/` exclude).

**Phase 2 work to date (cumulative):**

- `notes/inventory/11_phase_2_plan.md` — original plan (historical record).
- `notes/inventory/14_phase_2_revision_per_module_tests.md` — active spec.
- `notes/inventory/15_discovery_run_path_coverage.md` — discovery-run path coverage and fixture candidates.
- `tests/characterization/` — pytest framework, conftest, fixtures, README.
- `tests/characterization/helpers/` — six comparators (CSV-EXACT, CSV-STRUCT, CSV-EXACT-MODULO-PATHS, JSON-DEEP, JSON-MODULO-PATHS, PNG-PERCEPTUAL/EXISTS, TEXT-EXACT) + `ComparisonResult` dataclass.
- 67 `local_unit` tests passing.
- 184 `hpc-tier` characterization tests across 8 stage files, parametrized over `input_control_polyA, design_0_seq_0, design_13_seq_2`. **Restructuring required once per-module references are built (next concrete step #4).**
- Nextflow cache-busting (`path X_script` content-hashed).
- `containers/` Singularity definitions in repo with README.
- `scripts/sync_to_hpc.sh` + `.gitignore` hygiene + `scripts/rewrite_fixture_paths.py`.
- `audit_pre_hpc_roundtrip.md` and JAVA_HOME fix.
- Per-module fixtures in each `tests/<module>/data/`.

**Phase 1 deliverables** remain in `notes/inventory/` (`01_module_map.md` through `10_phase_1_synthesis.md`).

---

## Right Now (state mid-session)

Curation work complete on Mac. Per-module fixtures are staged but not yet exercised end-to-end on HPC. Eleven commits pending push (SIGPIPE fix, sync script fix, gitignore update, per-module test decoupling + module READMEs, path-coverage report, path-rewriter script, four per-module fixture commits, this state-doc update).

Anticipated tag for these commits: `phase-2.7-discovery-run-complete`. Tagging may be deferred until per-module tests have been exercised on HPC against the curated fixtures (next concrete step #1) — that's the first place the curation could turn out to be subtly wrong.

---

## Next Concrete Steps (in order)

1. **Sync curated fixtures to HPC and run each per-module test.** Use `./scripts/sync_to_hpc.sh`; then run `tests/<module>/run_test_<module>.slurm.sh` for each of the four populated modules (rosetta_filtering, proteinmpnn, negative_steering, orthogonal_metrics). Confirm each test produces a clean output tree without missing-file errors. The orthogonal_metrics test is the most likely to surface issues because of the path-rewriting and the cross-fixture reference into `negative_steering/data/design_pdbs/`.

2. **Subtractive rebuild** of each `tests/<module>/example_output_files/` from per-module test output. Same approach as the supervisor-demo rebuild (mirror structure, exclude PDBs/NPZs/MSA/predictions/ subtrees), but per-module rather than cohort-wide.

3. **Tag `phase-2.8-fixtures-and-references-complete`** once all per-module reference sets are committed.

4. **Restructure 184 characterization tests** to point at per-module reference paths and update parametrize lists. Fold in the small RFDiffusion fail-branch unit test as part of this work (the fail branch is a defensive code path the discovery run can't exercise; covering it requires a hand-constructed input + threshold tweak rather than a real fixture). Tag `phase-2.9-tests-restructured-complete`.

5. **Class 8 (`new_contamination`) standalone task.** Manufactured fixture or unit test against `classify_reversion_verdict` + `_classify_aggregated_verdict` to cover the path the discovery run didn't produce. Sequencing: ideally before phase-2-complete; can be slotted in around the test-restructuring work.

6. **Safety-net validation** — delete `bin/sequence_registry.py`, expect zero diffs. This is `phase-2-complete`.

7. **`tests/old_full_test/` deletion** on HPC and `tests/full_test_run/example_output_files/` retirement on Mac. Both stay until step 4 is done.

---

## Important Context Not Captured Elsewhere

### Why we pivoted

The original plan assumed the supervisor-demo run (~6 hours) would be the iteration target. Three problems forced a revision: (1) 6-hour iterations unworkable across many Phase 3 commits; (2) path coverage was retroactively-discovered, not designed-in; (3) per-module tests were silently chained, defeating their isolation. The pivot uses per-module tests as primary safety net with curated fixtures, full pipeline runs reserved for milestones. See `14_phase_2_revision_per_module_tests.md` §1.

### Negative-steering path semantics

The full pathway taxonomy lives in `notes/inventory/15_discovery_run_path_coverage.md` §Negative steering — four orthogonal axes (cold-start outcome / per-seed verdict / aggregated verdict / cohort tier) and the eight observable per-MPNN-sequence outcome classes the discovery run produced (or failed to produce, in the case of Class 8). That document is the authority on negsteer semantics for fixture purposes; the producer code at `bin/boltz2_iterate_steering.py` and `bin/cross_sequence_summary.py` is the authority for code-level questions.

### Reference-set sequence trio (currently)

Trio used by the existing 184 characterization tests:
- `input_control_polyA` — control, cold-start-only path
- `design_0_seq_0` — empty reversion (n_contaminated=0)
- `design_13_seq_2` — populated reversion (n_contaminated=33)

These names are baked into `tests/characterization/test_negsteer_per_sequence.py` and `TRACEABILITY.md`. **The new fixtures from the discovery run produce different names** — see `notes/inventory/15_discovery_run_path_coverage.md` for the locked candidate list (8 sequences across Classes 1, 2, 3a, 4, 5, 6, 7). Test restructuring is next concrete step #4.

### "Golden master" framing

Phase 2 reference outputs are *current pipeline behavior*, not verified-correct outputs. Tests assert *stability*, not *correctness*. A test failing during refactoring may indicate a legitimate behavior change or a fix to a latent bug — investigate, don't roll back automatically.

### Phase 3+ work currently deferred

- File-by-file simplification (Phase 4, after characterization tests are fully in place).
- Strict-syntax migration of `main.nf` and test workflows (deferred until HPC's Nextflow forces it).
- Closing cache-busting residual gaps for indirectly-loaded sub-scripts.
- Renaming `scaffold_rmsd` field (semantically motif RMSD).
- Investigating `merge_orthogonal_metrics.py` test-vs-production divergence.
- Investigating suspected `--n-cycles 1` silent-skip-reversion bug.

### Specific known issues

- **`bin/boltz2_iterate_steering.py`** is the largest file (5,500+ lines), four identical code blocks, references a non-existent submit script. Top refactoring target.
- **`bin/reversion.py::harvest_reversion_results`** has cyclomatic complexity 74. Pathological.
- **`bin/cross_sequence_summary.py::aggregate`** has CC 58. Second-worst.
- **`bin/rfdiffusion_plots.py`** several D-rated functions (CC 21–24).
- **Vulture's 8 high-confidence findings** (`08_vulture_high_confidence.txt`) are safest dead-code candidates; `bin/sequence_registry.py` is reserved for safety-net validation commit.

### HPC workflow constraint

Mac is authoritative. HPC has no git. Workflow: edit on Mac, sync via `./scripts/sync_to_hpc.sh`, run on HPC, iterate. Round-trips are slow.

### Tooling

- Mac: Homebrew, Node.js, Claude Code, miniforge conda, `receptor-tests` env (Python 3.10), Nextflow v26.04.0.0 (Mac-only syntax checks; `-stub-run` works).
- HPC home mounted at `/Volumes/HPC-Home/` (treat as read-only when prompting Claude Code; explicit guardrail).
- HPC: older Nextflow via Singularity (legacy parser default — that's why the codebase still works there despite v26 strictness on Mac).

### Branch and tag state

- Default branch: `main` (pristine baseline).
- Active branch: `remediation`.
- Tags: `baseline-pre-remediation`, `phase-1.1` through `phase-1-complete`, `phase-2.1-plan-complete` through `phase-2.5-reference-rebuild-complete`.
- **Next anticipated tag**: `phase-2.7-discovery-run-complete` — applied once the curated fixtures have been confirmed working on HPC (next concrete step #1). The previously-anticipated `phase-2.6-revision-plan-complete` was never tagged; the revision plan + decoupling work was done but the planned tag boundary was overtaken by the curation work and is being skipped.

### Anticipated remaining Phase 2 tag sequence

- `phase-2.7-discovery-run-complete` — once the curated fixtures are confirmed working on HPC (next concrete step #1). The eleven commits in this session sit at the boundary of this tag.
- `phase-2.8-fixtures-and-references-complete` — after per-module reference sets are subtractive-rebuilt and committed.
- `phase-2.9-tests-restructured-complete` — after the 184 characterization tests are updated and the small RFDiffusion fail-branch unit test is added.
- `phase-2-complete` — after Class 8 standalone task plus safety-net validation deletion produces zero diffs.

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
> - **How to verify:** Naming-only fix in Phase 3.3.

> - **What:** Cache-busting residual gaps — `bin/negative_steering_run_one.sh` and `cross_sequence_summary.py` use `--bin-dir` / `sys.path.insert` for sub-scripts; not content-hashed by Nextflow.
> - **Why suspicious:** Edits to indirectly-loaded scripts won't invalidate Nextflow cache.
> - **How to verify:** During Phase 3, prefer editing directly-tracked scripts; closing the gap requires declaring all sub-scripts as `path` inputs.

> - **What:** Pre-existing `fastrelax_xml` mismatch in `tests/orthogonal_metrics/test_orthogonal_metrics.nf`, corrected during cache-busting.
> - **How to verify:** Confirm nothing else relied on the old shape.

> - **What:** 184 existing hpc-tier characterization tests are tied to the supervisor-demo reference paths and the trio `input_control_polyA, design_0_seq_0, design_13_seq_2`.
> - **Why suspicious:** The revision plan replaces both the reference set location and the trio. Resolution paths and parametrize lists are stale.
> - **How to verify:** Restructure tests once per-module references are built.

> - **What:** Codebase uses legacy Nextflow syntax (top-level `workflow.onComplete { }` handlers, top-level `if` blocks in `main.nf:43`). `NXF_SYNTAX_PARSER=v1` pinned in slurm wrappers as stopgap.
> - **Why suspicious:** When HPC's Nextflow eventually upgrades and v1 parser is no longer available, the codebase will fail to compile.
> - **How to verify:** Migrate event handlers to `nextflow.config` or appropriate workflow scopes; lift top-level `if` into a workflow body. Test against strict-default Nextflow on Mac before declaring done.

> - **What:** Other potential `set -euo pipefail` + piped command sites that might have the same SIGPIPE-on-141 bug as the two we just fixed.
> - **Why suspicious:** The pattern is easy to miss because failures are intermittent (only triggers when output is large enough).
> - **How to verify:** Audit all `.nf` process bodies that combine `pipefail` with piped commands. Apply `|| true` to the diagnostic ones; use explicit error handling for any that should propagate failure.

> - **What:** `contact_cutoff` parameter is documented under HADDOCK section but consumed by `rfdiffusion_filter.py` regardless of mode.
> - **Why suspicious:** Misleading section header caused the discovery run's first attempt to use `contact_cutoff: 1.0`, breaking all Cα contact detection. Naming + grouping is bug-prone.
> - **How to verify:** During Phase 3, consider renaming to `rfdiff_contact_cutoff` and moving to its own section in `params_example.yml`. (Note: this entry is on the lighter side of "verification queue" and could just be a Phase 3 todo.)

> - **What:** Tier rule docstring at the top of `bin/cross_sequence_summary.py` (lines ~17–23) is stale relative to the implementation.
> - **Why suspicious:** Docstring claims `n_seeds_pose_holds == n_seeds`; implementation at line 200 uses `n_pass = n_seeds_pose_holds + n_seeds_clean_steered`. `pipeline_notes/pipeline_notes3.md` agrees with the implementation. Anyone reading the docstring will form the wrong mental model.
> - **How to verify:** Phase 3 docstring fix; one-line edit. Source: report 15 Open Q1.

> - **What:** Production survivor manifest does not gate on `cross_tier`.
> - **Why suspicious:** AF3 + biophys + rosetta are computed for tier-none representatives in the discovery cohort (95 of 122). In production with `negsteer_n_designs=20` this is meaningful GPU time. Whether the diagnostic-completeness behaviour is wanted, or whether a tier gate would save time without losing signal, is a product question.
> - **How to verify:** Decision rather than verification. Source: report 15 Open Q5.

> - **What:** `interface_plddt_median` column was blank in the `survivors_with_orthogonal_metrics.csv` row Claude Code inspected, but the `interface_plddt_too_low` gate fired on 70 sequences cohort-wide.
> - **Why suspicious:** The gate must be reading a differently-named column than what gets emitted in the cross_summary join. Possibly a column-rename inconsistency between the gate input and the published output.
> - **How to verify:** Trace `merge_orthogonal_metrics.py` to find which input column the gate consumes vs which output column it writes. Source: report 15 final VC row.

> - **What:** `cross_sequence_summary_with_interface_metrics.csv` columns `irmsd`, `fnat`, `dockq` may be silently wrong in production runs predating the `negsteer_interface_metrics.nf` chain-param fix.
> - **Why suspicious:** Same shape as the bug pipeline_notes10 fixed for the four other orthogonal-metrics modules. The module called `compute_interface_metrics.py` with `--effector-chain ${params.effector_chain}` (input-PDB convention, "C" in production) instead of `--effector-chain ${params.rfdiff_output_effector_chain}` (prediction-PDB convention, "B"). Script defaults to A/B and the script docstring states `--receptor-chain Must match the Boltz prediction chain`. With chain "C" passed, `_iface_pair == {receptor_chain, effector_chain}` filtering would fail to match the prediction PDB's actual {A,B} interface, so DockQ family columns would either fail with an obscure reason or silently produce wrong values.
> - **How to verify:** After the next per-module orthogonal_metrics test sweep produces correct columns, spot-check a survivor or two against any prior production CSV containing those columns. If they differ materially, alert anyone using historical interface-metrics data and consider whether downstream analyses need re-running.

> - **What:** Pipeline output verbosity — work/ and results/ trees contain many thousands of tiny files per run, dominated by per-seed Boltz internals, intermediate JSONs, and per-design subtrees that are not consumed downstream. Deleting a full run (e.g. tests/old_full_test/) takes minutes-to-hours of filesystem time on HPC; tarballing is slow; backups are inflated.
> - **Why suspicious:** Likely a mix of (a) Boltz emitting full prediction trees per seed (legitimately many files, possibly excessive), (b) producer scripts writing intermediate JSONs that no downstream stage reads, (c) per-design subtrees duplicated across cycles, and (d) workdirs not pruned because Nextflow defaults retain everything. The path-coverage curation already required excluding *.npz, *.a3m, msa/, predictions/, contamination_scratch/ from fixture tarballs to keep them tractable — that exclude list is itself evidence of the problem.
> - **How to verify:** During Phase 3, audit each producer script and Nextflow process for files emitted but never consumed downstream. Cross-reference against the orthogonal_metrics fixture exclude list (those subtrees are confirmed unread by the orthogonal_metrics stack). Likely fixes: scratch-only outputs into work/ rather than published; coalesced per-seed outputs; explicit deletion at the end of negsteer cycles. Phase 3 work; not blocking.

> - **What:** Per-module test cost asymmetry — negsteer takes ~55–60 minutes wall-clock; the other four take <10 minutes each.
> - **Why suspicious:** Not a bug, but a design constraint worth recording. Per-sequence wall-clock from the phase-2.7 sweep ranged from ~4 min (Class 1, cold_start_all_clean, skips steering) to ~57 min (Class 5, 12 contaminated mutations to revert). The median heavy sequence is ~30 min. Cost is dominated by the steering+reversion Boltz prediction loop. Current params (negsteer_n_designs=4, num_seeds=3) are at the floor for fixture path coverage — reducing either compromises which Classes 1–7 are exercised.
> - **How to verify:** Treat the four cheap tests as the default per-commit safety net during Phase 3 refactoring; reserve negsteer for substantive changes (new modules, Phase 3 deliveries, pre-tag verification). If a future Phase 3 change requires faster negsteer iteration, audit which Boltz predictions are necessary to exercise each path class — there may be opportunities to reduce work without losing coverage. Class 5's 57-minute outlier is the cost-driving case; if iteration is needed and Class 5 coverage can be deferred, swapping design_3_seq_1 for a different Class 5 candidate with fewer contaminations would help.
---

## Open Questions

1. **HADDOCK branch in the per-module strategy.** Branch A is out of Wave 1 scope. Should `tests/haddock/` get a per-module reference set under the new strategy, or stay deferred? *Suggested: defer until Branch A is redesigned.*
2. **Disk-space cost of per-module reference sets.** Worth measuring after the first per-module rebuild.
3. **Stochasticity in per-module outputs.** Even with fixed inputs, some Boltz-2 outputs may vary slightly. The first per-module round-trip will reveal which need strategy demoting.
4. **Class 8 standalone task — manufactured fixture vs unit test.** The discovery run produced no `new_contamination` examples at either per-seed or aggregated resolution. The path needs coverage before phase-2-complete. Two options: (a) manufacture an input that will trigger `new_contamination` aggregation through a small targeted negsteer run on HPC; (b) write a Python unit test against `classify_reversion_verdict` + `_classify_aggregated_verdict` with a hand-constructed reversion-result blob. (b) is cheaper but covers less; (a) is more representative but slower. Decision needed when this task is picked up.

---

## How to Use This Document

**At the start of a session:**

1. Read this document.
2. Read `notes/codebase_remediation_plan.md` if it's been a while.
3. Read `notes/inventory/10_phase_1_synthesis.md` for Phase 1 findings.
4. Read `notes/inventory/14_phase_2_revision_per_module_tests.md` for the active Phase 2 spec.
5. Begin work on the "Next Concrete Step" listed above.

**During a session:**

- If a decision is made, note it.
- If a doubt or suspicion arises, add it to the Verification Queue.
- If a question can't be resolved now, add it to Open Questions.

**At the end of a session:**

- Update "Just Completed" with what was done.
- Update "Next Concrete Step".
- Update "Last Updated" date.
- Commit and push.

This document is the canonical "where am I" source. If it disagrees with another document, this one wins (and the other should be updated).
