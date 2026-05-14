# Session Handoff — 2026-05-14 → next session

A focused starting brief for the next chat.  Read this first.  For the full picture, follow the cross-references to `notes/remediation_state.md`, `notes/phase4_architecture_spec.md`, and `notes/design_audit.md`.

## Read these first (in order)

1. **This file** — orientation and immediate next step.
2. **`notes/remediation_state.md`** — full session log + branch state.  Authoritative.
3. **`notes/phase4_architecture_spec.md`** — the architectural contract.  Especially §"User clarifications" (CL-1 to CL-5) and the §2.2 note about `BreakSegment` + `PassthroughSegment`.
4. **`notes/design_audit.md`** — six-session grill-me record.  **Session 7 (HADDOCK restructure) is the next planned grill-me** — research already begun this session, see §"HADDOCK code surface" below for current state.
5. **`~/.claude/.../memory/MEMORY.md`** — two recurring-failure memory files indexed (`project_no_reversion_semantics.md`, `project_contig_string_format.md`).  Read both before touching negsteer outcomes or contig strings.

## Where the code is

- `main` at commit `18718db`: safe-fallback baseline.  Phase 4 spec + pre-architecture work.  No Phase 4 implementation code.
- `phase4-impl` at commit `be593db`: **active branch.**  38 commits ahead of `main`.  All Phase 4 implementation, deep-form migrations, polish, fixture regeneration, full-pipeline test integration done.
- Working directory should be on `phase4-impl`.  Verify with `git branch -vv`.

## What happened in the last session (2026-05-13 → 2026-05-14)

17 commits on top of the prior `ea6dd71` (notes update).  Three logical groupings:

### A. Phase 4 polish (post-implementation, pre-verification)

- **`1641c76`** — `extract_survivor_manifest.py`: remap stale absolute paths (`/Users/...` → `/hpc-home/...` survives via fallback to discovered workdir + repo-root); exit non-zero on 0-survivor manifest (was silently exiting 0, letting Nextflow fan out the orthogonal cascade over an empty channel and report "Success" with no orthogonal data).
- **`e3c2752`** — polish pass: ERR trap on `tests/run_tests.sh` (prints failing line on `set -e` exit); F821 typing fixes (`Optional` import in `boltz2_negative_steering.py`); `NEGSTEER_CROSS_SEQUENCE` cache-busting via passing the entire `bin/` directory as a path input (closes the indirect-import cache-busting gap).

### B. Orthogonal cascade enhancements

- **`f5f9167`** → **`373e164`** → **`0bf94f6`**: orthog combined-cohort plot ranking churn.  Final state: `(cross_tier, boltz_conf_passes, -composite_score)` — tier first, rows that pass every thresholded Boltz-2 confidence metric next, then descending composite within each group.  Legends side-by-side instead of stacked.  ORTHOG_PLOTS wired into `tests/orthogonal_metrics/test_orthogonal_metrics.nf` so plots land in `plots/` (matching other modules' patterns; `--with-plots` iterator still publishes to `plots_iter/`).
- **`20f9ed0`** — `params.orthogonal_tier_filter` introduced.  Default `'all'` runs AF3 + biophysical + Rosetta on every steered design including tier-none failed designs (useful diagnostic).  `'abc'` restricts to survivors.  Plumbed through `extract_survivor_manifest.py` (new `--tier-filter` CLI flag) and `modules/negsteer_manifest.nf`.
- **`scripts/refresh_orthog_fixture.sh`** (also `20f9ed0`) — one-shot refresh of `tests/orthogonal_metrics/data/negsteer_run/` from the latest `tests/negative_steering/receptor_resurfacing_results/`.  Brings the orthogonal_metrics test fixture from a 2-sequence mini-cohort up to all 8 steered designs (+ 2 controls).

### C. Test infrastructure + fixture regeneration

- **`06c0157`** — `tests/haddock/test_haddock.nf`: override `params.haddock_min_cluster_size = 2` for small test cohorts (production default of 4 is too strict at `haddock_sampling=100`).
- **`ba95ede`**, **`fa14d3a`** — `tests/update_example_dataset.slurm.sh` and its impl: two sbatch-specific bugs fixed.  (1) `${BASH_SOURCE[0]}` resolves to `/var/spool/slurmd/job<id>/` under sbatch, not the real `tests/` dir; switched to `${PWD}` (set by `#SBATCH --chdir`).  (2) Path doubling: when CWD is already `tests/` and the user passes `tests/<module>/...` relative to project root, the path resolved as `tests/tests/<module>/...`.  Impl now auto-strips a leading `tests/` segment with a stderr NOTE.
- **`e5f3be0`** — **all 6 module fixtures regenerated** against `phase4-impl`.  657 files updated.  By module: negative_steering 627 (CL-3 + schema rename touches per-design × per-seed × per-stage), rfdiffusion 9 (RFD container `complex_beta` update changed designs), proteinmpnn 9 (small `design_region_score` shifts), rosetta_filtering 5 (`design[3].dG_separated` 33→57 — unexplained, blessed), orthogonal_metrics 6 (10-sequence cohort + ranking + ORTHOG_PLOTS), **haddock NEW** (first time HADDOCK has a fixture — Nextflow run-artefacts only; scientific outputs are placeholders pending the restructure).
- **`279f078`** + **`8e21c7b`** + **`be593db`** — `full_test_run` is now a first-class module in `tests/run_tests.sh`.  Dispatcher special-cases it to invoke the project-root `run_pipeline.slurm.sh` with `tests/full_test_run/params_full_test.yml`.  Two fixes to the params file: `haddock_sampling: 1 → 100` (validator's min is 100; value is unused in mode 2 but the validator runs unconditionally); `pdb_file` corrected to `af3_pikp1_native_avrpikf_complex.pdb` (the actual filename on disk).  `orthogonal_tier_filter` added to `validate_params.py:PARAM_SPECS`.

### Verification state at end of session

- **272 / 272 hpc-marked characterization tests passing.**  386 local_unit tests deselected by `-m hpc`.  26 expected skips.
- **All 6 per-module test runs reported `Success: true`** prior to fixture regeneration.
- **Full pipeline test (`./tests/run_tests.sh --modules full_test_run`) running on HPC at end-of-session.**  ETA ~24h+ depending on GPU queue.  Submission required two iterations to clear validator + a stale `pdb_file` path; the running attempt is `nextflow_pipeline_<latest jobid>` and should be checked first thing.

## What needs to happen FIRST in the next session

1. **Check the full pipeline test status:**
   ```
   ssh slurm
   squeue -u $USER
   tail -100 nextflow_pipeline_<latest-jobid>.out
   ```
   Expect a long-running process tree.  If it finished successfully: skim the output for any tier-none orthogonal data (the `orthogonal_tier_filter: all` setting should have produced AF3 + biophys + Rosetta for every steered design).  If it failed mid-run: bring the failing process work dir into the next session.

2. **If full test passed, start the HADDOCK Session 7 grill-me.**  Background research has been done; see §"HADDOCK code surface" below.  No notes added to `design_audit.md` yet — the next chat should open with Batch 1 of Q129–Q136.

3. **If full test failed and isn't a quick fix, fix and resubmit.**  Then revisit HADDOCK after the long run lands.

## What's left to do (in priority order)

### High priority

1. **HADDOCK module restructure** (Session 7 grill-me + iterative redesign).  User-flagged motivations:
   - Most energetically favourable HADDOCK pose never matched the intended pose.
   - Suspected AIR table issues.
   - Suspected cluster-selection issues (HADDOCK ranks by score, but the "right" cluster may not be the top one).
   - Wants more targeted contact restraints to hold the input pose.

2. **`rfdiff_contact_cutoff` / `haddock_hotspot_cutoff` parameter split** (deferred since pre-Phase-4) — fits naturally with the HADDOCK restructure since `EXTRACT_HOTSPOTS` is the HADDOCK consumer with the odd name.

### Medium priority

3. **Investigate the Rosetta `dG_separated` 33→57 drift** — fixture was blessed in `e5f3be0` but the cause wasn't isolated.  Inputs are frozen PDBs, container should be unrelated to the RFD container update.  May be Rosetta non-determinism, may be a real change.

4. **Validator coverage gaps** — 9 params still flagged as "no spec entry": `af2_data_dir`, `af3_db_v3`, `af3_model_dir`, `af3_package_id`, `boltz2_container`, `colabfold_container`, `rfdiff_container`, `max_af3_parallel`, `haddock_min_cluster_size`.  Infrastructure-type params from `nextflow.config`.  Add specs as part of a threshold audit.

### Low priority / deferred

5. **`contig_utils.parse_design_region` / `resolve_contigs` migration** to thin ContigSpec adapters.  External contract preserved; cosmetic refactor.

6. **`NegativeSteeringRun.from_workdir` against a real production workdir** — synthetic-fixture tests pass; full-data verification only possible once a real negsteer run is observed on `phase4-impl`.  The full pipeline test run should provide this.

7. **`read_ca_atoms` consolidation** — intentionally NOT consolidated; cross-reference docstrings added.  No functional issue.

8. **Class 8 (`new_contamination`) coverage gap** — user has a separate plan to surface a Class 8 example from the full test run.

9. **Stale Nextflow process selectors in `nextflow.config`** — `BOLTZ2_PREPARE`, `COLABFOLD_SEARCH_PER_DESIGN`, `BOLTZ2_PREDICT`, `BOLTZ2_VERIFY_BINDING`, `BOLTZ2_FILTER_AND_RANK`, `BOLTZ2_PLOTS`, `AGGREGATE_RESULTS` — config has labels for these but no workflow imports/calls them.  Warns at every pipeline start.  Cosmetic cleanup.

## HADDOCK code surface (research for Session 7 grill-me)

Already-read files:

- `bin/haddock3_prepare.py` (258 LOC) — parses contig to identify de novo regions as HADDOCK active residues; generates ambig_restraints.tbl.  **AIR generation is one big AND-of-OR: every receptor active residue must contact OR(every effector active residue).**  No active/passive distinction.  One-sided default (receptor → entire effector chain) when `effector_active_residues=""`; two-sided when given.
- `bin/extract_hotspots.py` (199 LOC) — post-docking interface-residue extraction from the docked complex; outputs RFDiffusion-format hotspot string.  Uses `find_interface_residues` from `haddock_utils.py`.
- `bin/collect_haddock3_dock.py` (291 LOC) — post-processes HADDOCK run dir; rejects clusters below `--min-cluster-size`.
- `bin/haddock3_plots.py` (695 LOC) — plotting code (genuine HADDOCK plots, untested by characterization).
- `bin/build_contigs.py` (229 LOC) — builds RFDiffusion contig string from selected docked output.
- `bin/haddock_utils.py` (283 LOC) — shared helpers.

The Nextflow processes (`modules/haddock.nf`):

- `HADDOCK3_PREPARE` — runs `haddock3_prepare.py`.
- `HADDOCK3_DOCK` — runs HADDOCK3.  docking.cfg has 8 modules: topoaa → rigidbody → seletop → flexref → emref → clustfcc → seletopclusts → caprieval.  Sampling = `params.haddock_sampling` (production: 10000), `min_population = params.haddock_min_cluster_size` (production: 4), seletop = 20.
- `HADDOCK3_PLOTS` — runs `haddock3_plots.py`.
- `EXTRACT_HOTSPOTS` — runs `extract_hotspots.py` on the chosen docked complex; uses `params.rfdiff_contact_cutoff` (currently shared with RFDiffusion's filter).
- `BUILD_CONTIGS` — runs `build_contigs.py`.

Initial grill-me question candidates (NOT YET CAPTURED IN design_audit.md — open them in Batch 1):

- **VERIFY: AIR generation has no passive residues; every receptor active residue is restrained to OR-of-every-effector-active-residue with distance 2.0 ± 2.0 ± 0.0 (one-sided) or 2.0 ± 2.0 ± 2.0 (two-sided).  Is this what you want, or is the lack of passive residues the underlying cause of poor cluster discrimination?**
- **VERIFY: `effector_active_residues` is an empty string by default (params_example.yml line where this lives).  When empty, every receptor active residue is restrained to "any effector residue" — effectively a very loose "must contact effector somewhere" constraint.  Is the empty default why intended-pose recovery is poor?**
- **VERIFY: Cluster selection in `collect_haddock3_dock.py` is by HADDOCK score (default cluster `1` is the best-scoring).  No structural similarity check against the intended pose.  Is this the "wrong cluster picked" mechanism you described?**
- **UNKNOWN: What does "intended pose" mean operationally — do you have a reference complex (e.g. AF3 prediction) that's the target, or is it more abstract (a region of receptor surface you want contacted)?**
- **UNKNOWN: Have you tried HADDOCK's pairwise distance restraints (not AIRs) — `unambig_restraints.tbl` with explicit residue-pair distances?**
- **UNKNOWN: Is the goal to RECOVER an existing input pose, or to GENERATE a pose that holds specific contacts?**

## Critical context the next session must know

### CL-3 reversion-gating rule

**Old rule** (pre-Phase-4): any single contaminated (design, seed) triggers reversion for that design.

**New rule** (post-Phase-4, live in `cmd_build_contaminated`): reversion runs for a design IFF
```
n_correctly_placed > 0  AND  n_contaminated >= ceil(n_correctly_placed / 2)
```

Documented in `notes/phase4_architecture_spec.md` §"User clarifications" CL-3.  Encoded both in `bin/stage_result.py:triggers_next_stage` (the type-level encoding) and in `bin/boltz2_iterate_steering.py:cmd_build_contaminated` (the live pipeline call site).  Six dedicated tests in `tests/characterization/test_stage_result.py::TestCL3SteeredGating`.

### `no_reversion` is NOT cold-start

Recurring confusion documented in memory file `project_no_reversion_semantics.md`.  The outcome label `no_reversion` fires in TWO disjoint paths:
- (a) Cold-start `skip_steering`: all seeds passed at cycle 0, never went through steering.
- (b) Steering ran but no seed had contamination on mutated positions.

`n_pass` is path-agnostic: ALWAYS `pose_holds_count + clean_steered_count`.

### Contig string format

Canonical: `A1-10/5/A15-20 B`.  Slash-separated within a chain, space-separated between chains, never commas.  Fixed segments chain-prefixed; denovo segments are bare lengths (resolved) or length ranges (constraint).  As of Phase 4, `DeNovoSegment` always carries `(min_len, max_len)` — `min_len == max_len` for the resolved case.  RFDiffusion chain-break marker `0` parses as `BreakSegment` (no length, ignored by position-math methods).  Bare chain letters within a block (e.g. `A1-10/B/A15-20`) parse as `PassthroughSegment(chain)`.

### Orthogonal cascade scope

`params.orthogonal_tier_filter` controls which `cross_tier` values receive the AF3 + biophysical + Rosetta cascade:
- `'all'` (default): every steered design including tier-none failed designs.  Useful diagnostic.  More expensive (per-design AF3).
- `'abc'`: restrict to cross_tier in (A, B, C).  Lower GPU cost.

Controls (row_type != steered) are not currently filtered out by `extract_survivor_manifest.py` — so with `'all'` they also get the cascade.  Documented as a follow-up choice in the previous session's notes; carried forward.

### The 13 deep-module types and their tiers

```
T0 (no Phase-4 deps):
  ContigSpec, BoltzConfidenceMetrics, AF3ConfidenceAggregate,
  PipelineParams, PipelineInternalThresholds
T1: PositionSet                       (deps: ContigSpec)
T2: ProteinStructurePrediction        (deps: PositionSet)
T3: DesignedBackbone                  (deps: PSP, ContigSpec, PositionSet)
T4: DesignedSequence, StageResult, OrthogonalMetrics
T5: NegativeSteeringRun
T6: DesignCohort                      (top of hierarchy)
```

Each type's full spec is in `notes/phase4_architecture_spec.md` §2.1 through §2.13.

### Two memory files to be aware of

`~/.claude/projects/-Users-jowillia-Documents-GitHub-receptor-resurfacing-pipeline/memory/`:
- `project_no_reversion_semantics.md` — read before touching negsteer outcomes
- `project_contig_string_format.md` — read before writing any contig example

Both indexed in `MEMORY.md`.

## How to verify the session's work locally before the next one

```bash
# On Mac, on phase4-impl:
git status                                                   # should be clean
git log --oneline 18718db..HEAD | wc -l                      # expect 38
python3 -m pytest tests/characterization/ -q -m local_unit \
    --ignore=tests/characterization/test_plots.py \
    --ignore=tests/characterization/helpers/tests/test_png_compare.py
# Expect: 386 passed, 7 skipped (gemmi-dependent), 272 deselected
```

If any of those checks fail, something has changed since end-of-session — investigate before continuing.

## When to update this file

Update at the END of every session, not as you work.  Replace the "What needs to happen FIRST" section with the new immediate next step.  Update the branch state if branches move.  Preserve "What happened in the last session" so the next reader can trace progress.
