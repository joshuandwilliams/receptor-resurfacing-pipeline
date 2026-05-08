# Remediation State

A living document tracking where the codebase remediation effort currently stands. Read this at the start of every session; update it at the end of every session.

**Last updated:** 2026-05-08 (experiments branch merged; RFDiffusion improvements; Phase 4 next)

---

## Current Phase

**Phase 3 (Mechanical Cleanup) — COMPLETE.**

Tag: `phase-3-complete`.

Dead code removed, key naming inconsistencies fixed, critical duplication consolidated, two scientific correctness bugs fixed, and the three worst complexity hotspots reduced from F-rated to C/D/E. Safety net confirmed green (272 passed, 26 skipped, 0 failed) after all changes.

---

## Just Completed (this session)

### experiments branch merged into remediation (2026-05-08)

All commits from the `experiments` worktree have been fast-forwarded into `remediation`. The following pipeline and infrastructure changes are now in both branches:

**RFDiffusion improvements (`b68e66c`)**
- `params.rfdiff_checkpoint` exposed as a first-class parameter (default `Complex_beta_ckpt.pt` — the PPI-optimised model). Old code used the container default; this makes the choice explicit and reproducible.
- Guiding potentials wired: `add_potential`, `rfdiff_guide_scale/decay`, `rfdiff_interface_weight`, `rfdiff_rog_weight`, `rfdiff_rog_min_dist` added to `main.nf`, `modules/rfdiffusion.nf`, `params_example.yml`.  Default `add_potential = false` so existing campaigns are unaffected unless they set it.
- `stop_after_rosetta` param added: halts pipeline cleanly after Sc filtering; `-resume` picks up from MPNN.
- `notes/inventory/rfdiffusion_parameter_audit.md` added — covers checkpoint landscape, guiding potentials guide, and changelog.

**Container rebuild (`b4780ee`, `756dcb8`, `cda2bbc`, `daf5c8a`)**
- `containers/LRR_Pipeline.def` renamed to `containers/HADDOCK_RFDiffusion_ProteinMPNN_MMseqs2.def`.
- `Complex_beta_ckpt.pt` download added; container rebuilt as `HADDOCK_RFDiffusion_ProteinMPNN_MMseqs2.img`.
- ColabDesign install switched from `git+https` to tarball (URL pointed to a branch, not a tag).

**Test data replacement (`06d5746`)**
- `tests/rfdiffusion/data/af3_pikp1_native_avrpikf_complex.pdb` replaced by `tests/rfdiffusion/data/pikp1_avrpikf_complex.pdb`.
- `tests/rfdiffusion/test_rfdiffusion.nf` and `tests/full_test_run/params_full_test.yml` updated to reference the new filename (stale reference was a broken-path bug — test would have failed on HPC).

**Negsteer bug fix (`7ce321c`)**
- `bin/negative_steering_run_one.sh`: `skip_steering` plan exit is now a soft failure, not a hard abort. Plan stage exits 1 with `skip_steering=true` in `plan.json`; old code was calling `fail()` before reaching the skip-steering check, aborting the pipeline and scancelling in-flight tasks.

**Infrastructure and cleanup**
- `scripts/clean_finished_run.sh` / `clean_all_finished_runs.slurm.sh` added (storage bloat audit tooling).
- `notes/inventory/17_storage_bloat_audit.md` added.
- All test SLURM plot scripts updated to reference new container image name.
- `scripts/sync_to_hpc.sh` updated to include experimental outputs.

**Documentation fix**
- `notes/inventory/rfdiffusion_parameter_audit.md` "Changes made" section had three stale `Complex_base_ckpt.pt` references (draft artifact); corrected to `Complex_beta_ckpt.pt` to match the actual code.

---

### Phase 2 closing (earlier today)

- Phase 2.7/2.8 spot-check passed. Five per-module tests green on HPC.
- Subtractive rebuild complete. Five per-module reference sets committed.
- Phase 2.9 complete: 184 → 298 parametrized hpc-tier tests restructured to per-module roots; SEQUENCE_COHORT (10 sequences); RFDiffusion fail-branch unit test added; TRACEABILITY.md updated.
- `pytest_runner.def` container built at `/hpc-home/jowillia/singularity/pytest/pytest_runner.img`.
- `run_pytest.slurm.sh` added at `tests/characterization/`.
- Characterisation suite green on HPC: 272 passed, 26 skipped, 0 failed.
- `bin/sequence_registry.py` deleted (safety-net validation — confirmed dead by vulture).
- `tests/full_test_run/example_output_files/` retired.
- `#SBATCH --chdir` added to all 14 SLURM scripts.

### Phase 3 — full account

**3.1a — Dead code removal**
- `ruff check --fix`: 125 fixes across 28 files (F401 unused imports, F541 f-string placeholders).
- 8 vulture high-confidence items resolved: 3 unused imports, 2 unused parameters (with call-site cleanup), 2 unreachable post-return blocks, 1 unused variable.
- Ruff errors: 196 → 70 remaining (E702/E741/F821/E402 — addressed in later steps or deferred to Phase 4).

**3.1b** — `sequence_registry.py` already deleted in Phase 2.

**3.2 — Comment cleanup** — deferred. Too high risk of removing load-bearing rationale in a domain-heavy codebase. Revisit per-file in Phase 4.

**3.3 — Naming fixes**
- `scaffold_rmsd` → `motif_rmsd`, `scaffold_str` → `motif_str`, `calc_scaffold_and_region_metrics` → `calc_motif_and_region_metrics`, `scaffold_diffs` → `motif_diffs` (Baker-lab convention, previously inverted). Reference JSONs updated.
- `THREE_TO_ONE` consolidated from 5 definitions across 5 files into single canonical definition in `bin/contig_utils.py` (27-entry dict). All consumers import from there. `boltz2_negative_steering.py` header comment updated (no longer standalone-droppable).
- `contact_cutoff` → `rfdiff_contact_cutoff` in params/main.nf/modules/rfdiffusion chain. Unrelated `contact_cutoff` parameters in negsteer/compute_metrics/reversion scripts left as-is (independent semantics).
- Tier-rule docstring in `bin/cross_sequence_summary.py` updated to accurately describe `n_pass = n_seeds_pose_holds + n_seeds_clean_steered`.
- `experiments/` added to `.gitignore`.

**3.4a — merge_orthogonal_metrics.py AF3 gating bug (scientific correctness fix)**
- Production `bin/merge_orthogonal_metrics.py` was gating designs on AF3 disagreement. Intended behaviour (pipeline_notes13, test copy) demotes AF3 to informational flag only. Fixed: `gating_flags = [f for f in flags if not f.startswith("af3_nomsa")]`. 8 differences between production and test copy — all class-a, all applied.

**3.4b — resubmit_pathway_chains helper**
- Extracted shared helper from four identical subprocess blocks in `bin/boltz2_iterate_steering.py`. Complexity reductions: cmd_iterate_collect D(21)→C(18), cmd_kickoff C(19)→C(16), both finalize F(50)→F(47).

**3.4c — _REVERTED_CONFIDENCE_FIELDS single source of truth**
- 24-entry constant moved from `bin/boltz2_iterate_steering.py` to `bin/reversion.py`. boltz2_iterate_steering now imports it. No mismatches found. Docstring explains deliberately-excluded fields and six-bug history.

**3.5a+b — Sub-function extraction**
- `reversion.py::harvest_reversion_results`: CC 74 → 20. Nine helpers extracted, all grade A–B.
- `cross_sequence_summary.py::aggregate`: CC 58 → 21. Six helpers + sentinel exception extracted, all grade A–C.

**3.5c — Sub-function extraction (cmd_compute_final_metrics)**
- `boltz2_iterate_steering.py::cmd_compute_final_metrics`: CC 72 → 31 (E). Six helpers extracted. Part of dormant multi-cycle chain — not exercised by per-module tests.

**Verification queue investigations**
- **VQ: interface_plddt_median blank** — resolved as non-bug. Gate correctly reads `representative_interface_plddt_median`. Stale docstring in `merge_orthogonal_metrics.py` line 29 fixed.
- **VQ: irmsd/fnat/dockq chain-param** — resolved as already fixed (commit `a4ef1f7`, 2026-05-03). **Historical runs before 2026-05-03 may have wrong irmsd/fnat/dockq values** — populated but wrong (fallback used wrong chain). Alert downstream consumers.

**HPC pytest sweep — green**
- 272 passed, 26 skipped, 0 failed. One intermediate failure on stale rfdiffusion cache after 3.3 rename; fixed by resubmitting rfdiffusion per-module test.

---

## Phase 3 Deliverables (cumulative, on top of Phase 2)

- `bin/merge_orthogonal_metrics.py` — AF3 gating bug fixed; stale docstring fixed.
- `bin/boltz2_iterate_steering.py` — resubmit_pathway_chains helper; _REVERTED_CONFIDENCE_FIELDS import; cmd_compute_final_metrics CC 72→31; six remaining F-rated cmd_* functions (CC 45–50, deferred to Phase 4).
- `bin/reversion.py` — _REVERTED_CONFIDENCE_FIELDS canonical definition; harvest_reversion_results CC 74→20.
- `bin/cross_sequence_summary.py` — aggregate CC 58→21; tier docstring fixed.
- `bin/contig_utils.py` — canonical THREE_TO_ONE definition (27 entries).
- `bin/rfdiffusion_filter.py` — scaffold→motif rename; rfdiff_contact_cutoff rename.
- `bin/rfdiffusion_plots.py` — scaffold→motif rename; BoundaryNorm unused import removed.
- `bin/negsteer_plots.py` — unused parameter + unreachable block removed.
- All 14 SLURM scripts — `#SBATCH --chdir` added.
- `experiments/` in `.gitignore`.
- Reference JSONs updated for motif_rmsd and rfdiff_contact_cutoff key renames.
- Tags: `phase-3.1a-complete`, `phase-3.3-complete`, `phase-3.4a-complete`, `phase-3.4b-complete`, `phase-3.4c-complete`, `phase-3.5-complete`, `phase-3-complete`.

---

## Next Concrete Steps — Phase 4

Phase 4 is **structural refactoring toward deep modules**. Interface design decisions are required before touching code.

**Before writing any code:**
1. Run the grill-me skill to design the target module architecture. Key inputs: `notes/inventory/04_functional_categorization.md`, `notes/inventory/06_ubiquitous_language.md`, `notes/inventory/10_phase_1_synthesis.md`. Output: a short architecture document committed to `notes/`.

**Targets in priority order:**
1. **`bin/boltz2_negative_steering.py`** — already partly structured as library + CLI. Extract the public interface explicitly. First real deep module.
2. **`bin/boltz2_iterate_steering.py`** — 6 remaining F-rated functions (CC 45–50). Decompose by `cmd_*` stage boundaries; each stage communicates via JSON files whose interfaces are already implicit.
3. **`reversion.py`** — CC now 20, but mixes reading/classifying/building/error-handling concerns. Phase 4 extracts those into separate modules.
4. **`bin/merge_orthogonal_metrics.py`** — eliminate test/production divergence permanently by making the test import from `bin/`.
5. **Remaining ruff errors** — E702/E701 (semicolons/colons), E741 (ambiguous names), F821 (undefined names — review for latent bugs first), E402 (import order).

**Per-commit workflow — same as Phase 3:**
- `pytest -m local_unit` locally after every change.
- Four cheap per-module tests on HPC for non-negsteer changes.
- Full pytest sweep after each sub-phase.
- negsteer per-module test only for negsteer-producer changes.

---

## Important Context Not Captured Elsewhere

### "Golden master" framing

Phase 2 reference outputs are *current pipeline behavior*, not verified-correct outputs. Tests assert *stability*, not *correctness*. A test failing during Phase 4 refactoring may indicate a legitimate behavior change or a fix to a latent bug — investigate, don't roll back automatically.

### Negative-steering path semantics

The full pathway taxonomy lives in `notes/inventory/15_discovery_run_path_coverage.md` §Negative steering — four orthogonal axes (cold-start outcome / per-seed verdict / aggregated verdict / cohort tier) and eight observable per-MPNN-sequence outcome Classes. That document is the authority on negsteer semantics for fixture purposes.

### Class 8 (`new_contamination`) coverage gap

The discovery run produced no `new_contamination` examples. Options: (a) manufacture an input via a small targeted negsteer run; (b) write a Python unit test against `classify_reversion_verdict` + `_classify_aggregated_verdict` with a hand-constructed reversion-result blob. Option (b) is cheaper; option (a) is more representative.

### Known expected skips in pytest (26 permanent)

- `cycle_0/effector_template.cif` — `.cif` excluded from subtractive rebuild across all 10 sequences.
- Cold-start per-sequence files — `design_28_seq_1` lacks `cycle_statistics.csv`, `cycle_0/passing.json`, `cycle_0/summary.txt`, `steered_results_aggregate.csv`, `steered_results.csv`, `{true,wrong}_interface_residues.txt`.
- `inputs/receptor.fasta` — both controls use `control_*_receptor.fasta` variants.
- `cross_sequence_summary_with_interface_metrics.csv` — pinned in orthogonal_metrics, not negsteer_cohort.
- `survivors_with_orthogonal_metrics.csv` — AF3/biophys/rosetta did not run on per-module fixture.
- 4 `orthogonal_*.png` plots — same reason.
- 2 preprocessing tests — no per-module preprocessing fixture.

### Historical data warning

Any `cross_sequence_summary_with_interface_metrics.csv` produced before 2026-05-03 (commit `a4ef1f7`) may have wrong `irmsd`/`fnat`/`dockq` values. Values look populated (fallback path used wrong chain) rather than blank. Alert downstream consumers of historical interface-metrics data.

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
- Tags: `baseline-pre-remediation`, `phase-1.1` through `phase-1-complete`, `phase-2.1-plan-complete` through `phase-2-complete`, `phase-3.1a-complete` through `phase-3-complete`.

### Phase 4+ work currently deferred

- Deep module extraction (`boltz2_negative_steering.py` library/CLI split, `boltz2_iterate_steering.py` stage decomposition).
- Comment cleanup (3.2 deferred — too high risk without per-file domain review).
- Strict-syntax migration of `main.nf` and test workflows.
- Closing cache-busting residual gaps (`bin/negative_steering_run_one.sh`, `cross_sequence_summary.py` via `--bin-dir` / `sys.path.insert`).
- Class 8 (`new_contamination`) coverage.
- HADDOCK branch (Branch A) per-module test.
- `rfdiff_contact_cutoff` / `haddock_hotspot_cutoff` split.
- F821 undefined-name ruff errors — review for latent runtime bugs before fixing.
- Six remaining F-rated `cmd_*` functions in `boltz2_iterate_steering.py` (CC 45–50).

---

## Verification Queue

> - **What:** RFDiffusion per-module reference set needs regeneration on HPC.
> - **Why:** Three things changed simultaneously: (1) test PDB replaced (`af3_pikp1_native_avrpikf_complex.pdb` → `pikp1_avrpikf_complex.pdb`); (2) explicit checkpoint now passed (`Complex_beta_ckpt.pt` — different from whatever the old container default was); (3) seven new `val` inputs added to the RFDIFFUSION process, invalidating the Nextflow cache. The `rfdiffusion_metrics.json` reference uses JSON-DEEP with exact Cα coordinates — any of these three changes will cause the comparison to fail. Also consider whether to demote `rfdiffusion_metrics.json` from JSON-DEEP to JSON-STRUCT (keys/types/ranges only) given that RFDiffusion is stochastic and seed honoring is not guaranteed end-to-end.
> - **How to verify:** On HPC: delete `tests/rfdiffusion/work/`, run `sbatch tests/rfdiffusion/run_test_rfdiffusion.slurm.sh`, copy the `receptor_resurfacing_results/` output to `tests/rfdiffusion/example_output_files/`, commit, re-run `pytest -m hpc` against the new reference.

> - **What:** Suspected `--n-cycles 1` silent-skip-reversion bug from notes6.
> - **Why suspicious:** Synthesis §5 records this as not verified during Phase 1.
> - **How to verify:** Inspect `pathways.json` and reversion JSONs for cohort sequences; check whether reversion was actually attempted on contaminated cases.

> - **What:** Cache-busting residual gaps — `bin/negative_steering_run_one.sh` and `cross_sequence_summary.py` use `--bin-dir` / `sys.path.insert`; not content-hashed by Nextflow.
> - **Why suspicious:** Edits to indirectly-loaded scripts won't invalidate Nextflow cache.
> - **How to verify:** During Phase 4, prefer editing directly-tracked scripts; closing the gap requires declaring all sub-scripts as `path` inputs.

> - **What:** Legacy Nextflow syntax (`workflow.onComplete {}` handlers, top-level `if` in `main.nf:43`). `NXF_SYNTAX_PARSER=v1` pinned in SLURM wrappers as stopgap.
> - **Why suspicious:** When HPC's Nextflow eventually upgrades, v1 parser removal will break compilation.
> - **How to verify:** Migrate event handlers to `nextflow.config`; lift top-level `if` into workflow body. Test against strict-default Nextflow on Mac.

> - **What:** Other potential `set -euo pipefail` + piped command sites with the SIGPIPE-on-141 bug.
> - **Why suspicious:** Pattern is easy to miss; failures are intermittent.
> - **How to verify:** Audit all `.nf` process bodies combining `pipefail` with piped commands.

> - **What:** `params.rfdiff_contact_cutoff` consumed by both RFDIFFUSION_FILTER and EXTRACT_HOTSPOTS (HADDOCK). The `rfdiff_` prefix is mildly misleading for the HADDOCK consumer.
> - **Why suspicious:** Should arguably be two separate params sharing the same default.
> - **How to verify:** Phase 4 parameter interface review — split into `rfdiff_contact_cutoff` and `haddock_hotspot_cutoff`.

> - **What:** Production survivor manifest does not gate on `cross_tier`.
> - **Why suspicious:** AF3 + biophys + rosetta computed for tier-none representatives. Meaningful GPU time in production.
> - **How to verify:** Product decision — whether diagnostic completeness is wanted or a tier gate would save time without losing signal.

> - **What:** Pipeline output verbosity — work/ and results/ trees contain many thousands of tiny files per run.
> - **Why suspicious:** Deleting a full run takes minutes-to-hours on HPC filesystem.
> - **How to verify:** Phase 4 audit — identify files emitted but never consumed downstream.

> - **What:** Per-module test cost asymmetry — negsteer ~57 min vs. <10 min for other four.
> - **Why suspicious:** Design constraint, not a bug. Class 5's 57-minute outlier (design_3_seq_1, 12 contaminated mutations) is the cost-driving case.
> - **How to verify:** Treat four cheap tests as the default per-commit safety net; reserve negsteer for substantive changes touching negsteer producers.

> - **What:** Historical `irmsd`/`fnat`/`dockq` values wrong in runs before 2026-05-03.
> - **Why suspicious:** Chain-param bug fixed in commit a4ef1f7; pre-fix runs used input-PDB effector chain "C" instead of prediction-PDB "B". Values look populated rather than blank.
> - **How to verify:** After next per-module orthogonal_metrics sweep, diff irmsd/fnat/dockq against any prior production CSV. Alert downstream consumers if materially different.

---

## Open Questions

1. **Class 8 (`new_contamination`) coverage.** Manufactured fixture on HPC vs. Python unit test?
2. **HADDOCK branch per-module test.** Deferred until Branch A is redesigned.
3. **Stochasticity in per-module outputs.** First round-trip showed no issues, but future Boltz-2 version changes may require comparator strategy demotions (e.g. JSON-DEEP → JSON-STRUCT).

---

## How to Use This Document

**At the start of a session:**

1. Read this document.
2. Read `notes/codebase_remediation_plan.md` if it's been a while.
3. Read `notes/inventory/10_phase_1_synthesis.md` for Phase 1 findings.
4. Begin work on the first item under "Next Concrete Steps".

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