# Session Handoff — 2026-05-12 → next session

A focused starting brief for the next chat.  Read this first.  For the full picture, follow the cross-references to `notes/remediation_state.md`, `notes/phase4_architecture_spec.md`, and `notes/design_audit.md`.

## Read these first (in order)

1. **This file** — orientation and immediate next step.
2. **`notes/remediation_state.md`** — full session log + branch state.  Authoritative.
3. **`notes/phase4_architecture_spec.md`** — the architectural contract.  Especially §"User clarifications" (CL-1 to CL-5).
4. **`~/.claude/.../memory/MEMORY.md`** — two recurring-failure memory files indexed (`project_no_reversion_semantics.md`, `project_contig_string_format.md`).  Read both before touching negsteer outcomes or contig strings.

## Where the code is

- `main` at commit `18718db`: safe-fallback baseline.  Has the full Phase 4 spec + all pre-architecture work.  No Phase 4 implementation code.
- `phase4-impl` at commit `54dec6f` (or further if more landed): **active branch.**  All 13 deep-module types implemented + six caller migrations applied.
- Working directory should be on `phase4-impl`.  Verify with `git branch -vv`.

## What happened in the last session

Two major chunks:

**Phase 1 — pre-architecture work** (committed on `remediation`, merged to `main`).
Six items from the audit's "Audit close → Immediate actions" list:
- Schema rename: `aggregated_verdict → outcome`, `representative_ → rep_`, five `n_seeds_*` → `n_pass`.  60 files changed including 37 fixture CSVs.  Rides with a latent-bug fix: tier is now derived purely from `n_pass / n_seeds` (no more `no_reversion → tier A` shortcut).
- Cohort summary visual groups, MPNN-sequence join, `validate_params.py`, glossary updates, test runner / update-fixture scripts.

**Phase 2 — Phase 4 implementation** (on `phase4-impl`).
- All 13 deep-module types built (Tiers 0–6) — ~4,000 LOC new code + 350 local_unit tests.
- **CL-3 fix live in `cmd_build_contaminated`**: the headline behaviour change.  Reversion now uses the majority-of-correctly-placed rule per design instead of the per-(design, seed) gate.
- Six caller migrations: test/prod plot-script dedup (~4,200 LOC removed), `get_chain_sequence` dedup (2 of 4), `THREE_TO_ONE`/`_AA3TO1` dedup, `parse_af3_output.py` rewritten as 104-line wrapper, `extract_passing._compute_confidence_flag` reads from `PipelineInternalThresholds`.

Full detail in `notes/remediation_state.md` § "What landed in this session".

## What needs to happen FIRST in the next session

**Verify CL-3 against real cohort data before doing any more work.**  The HPC tests for rfdiffusion and negsteer were queued at end-of-session and had not run yet.

Concrete steps:

1. **Check HPC queue status:**
   ```
   ssh slurm
   squeue -u $USER
   ```
   Look for the `nf_test_rfdiff` and `nf_test_negsteer` jobs.

2. **If they've completed: read the logs.**
   ```
   ls -lt tests/negative_steering/slurm_*.out | head -3
   tail -50 tests/negative_steering/slurm_<latest-jobid>.out
   ```
   Look for `Success: true` and check the new `[build-contaminated] CL-3 gating: X / Y designs trigger reversion (Z contaminated entries queued)` line.  X should be lower than Y (some designs not triggering reversion under the new rule).

3. **If the test passed but with diffs, regenerate the fixture:**
   ```
   sbatch tests/update_example_dataset.slurm.sh \
       --module negative_steering \
       --updated-output-folder tests/negative_steering/receptor_resurfacing_results
   ```
   Then pull back to Mac:
   ```
   ./scripts/sync_from_hpc.sh --module negative_steering
   git diff --stat tests/negative_steering/example_output_files/
   ```
   The expected diff shape: some files removed (designs that no longer trigger reversion lose their `reversions/` subtrees), some designs that were tier C/none under the old rule may now be tier A/B.

4. **Run characterization tests against `phase4-impl`:**
   ```
   sbatch tests/characterization/run_pytest.slurm.sh
   ```
   Expect some failures driven by CL-3 changes.  Each failure is either (a) legitimate new behaviour we should bless, or (b) an unintended regression.  Investigate before piling on more migrations.

5. **If the test FAILED**: the most likely culprits are the CL-3 fix (`b67b4b3`) and the parse_af3_output rewrite (`8230ee5`).  Bisect against the pre-implementation baseline `1f735c0` if needed.

## What's left to do (in priority order)

These are all on `phase4-impl`; commit each separately, ideally with HPC validation between.

### High priority (architectural meat)

1. **`cross_sequence_summary.py` → `DesignCohort.to_cross_summary_csv`** — the biggest remaining migration.  ~900 LOC of production aggregator to replace.  Requires implementing `DesignCohort.from_runs_directory` and `NegativeSteeringRun.from_workdir` (currently stubs).  HIGH RISK.  Dedicated session.

2. **ContigSpec parser consolidation** — 4 existing parsers (`derive_input_design_region._parse_contigs`, `haddock3_prepare.parse_contig_segments`, `pipeline_correct_sequences.parse_contig_segments`, plus `contig_utils.parse_block_segments`) → `ContigSpec`.  Requires resolving the constraint-form-vs-resolved-form question:
   - Option A: revise `ContigSpec` (§2.2 of the architecture spec) to accept constraint-form denovo (range like `5-7`) alongside resolved-form (`5`).
   - Option B: introduce a separate `ContigConstraint` type for the pre-resolution form.
   Discuss with user before implementing.

### Medium priority

3. **Full threshold migration** — `orthogonal_metrics_plots.py` (`ORTHOG_SC_MIN`, etc.), `compute_metrics.py` defaults, `main.nf` inline values.  Each migration is small but they're scattered across ~10–15 files.  Pattern established in commit `54dec6f` (lazy import + fallback to literals).

4. **`find_contact_residues_heavy` / `read_ca_atoms` consolidation** — DIFFERENT semantic contracts between implementations.  Safe consolidation requires renaming for disambiguation, not merging.  Coordinate with user on naming.

### Low priority

5. **`pipeline_correct_sequences.get_pdb_sequence`** — different signature `(seq, resnums)` tuple.  Migrating requires changes to every caller.

6. **THREE_TO_ONE dedup in other files** — `boltz2_negative_steering.py`, `build_contigs.py`, `extract_hotspots.py`, `pipeline_correct_sequences.py` may still have local references; spot-check.

## Critical context the next session must know

### CL-3 reversion-gating rule (the headline behaviour change)

**Old rule** (pre-Phase-4): any single contaminated (design, seed) triggers reversion for that design.

**New rule** (post-Phase-4, live in `cmd_build_contaminated`): reversion runs for a design IFF
```
n_correctly_placed > 0  AND  n_contaminated >= ceil(n_correctly_placed / 2)
```

Examples with `num_seeds=3`:
- 3/3 correctly placed, 1 contaminated → **no reversion** (Tier-B-shape result holds; the contaminated seed counts as a failure in `n_pass` aggregation but isn't reverted).
- 3/3 correctly placed, 2 contaminated → **reversion runs**.
- 2/3 correctly placed, 1 contaminated → **reversion runs** (50%+ gate).
- 0 correctly placed → no reversion (nothing to revert).

Documented in `notes/phase4_architecture_spec.md` §"User clarifications" CL-3.  Encoded both in `bin/stage_result.py:triggers_next_stage` (the type-level encoding) and in `bin/boltz2_iterate_steering.py:cmd_build_contaminated` (the live pipeline call site).  Six dedicated tests in `tests/characterization/test_stage_result.py::TestCL3SteeredGating`.

### `no_reversion` is NOT cold-start

Recurring confusion documented in memory file `project_no_reversion_semantics.md`.  The outcome label `no_reversion` fires in TWO disjoint paths:
- (a) Cold-start `skip_steering`: all seeds passed at cycle 0, never went through steering.
- (b) Steering ran but no seed had contamination on mutated positions.

`n_pass` is path-agnostic: ALWAYS `pose_holds_count + clean_steered_count`.  A `no_reversion`-labeled group with `n_pass < n_seeds` correctly lands in a lower tier (this was the latent bug the renames fixed).

### Contig string format

Canonical: `A1-10/5/A15-20 B`.  Slash-separated within a chain, space-separated between chains, never commas.  Fixed segments chain-prefixed; denovo segments are bare lengths (resolved form) or length ranges (constraint form, `5-7`).  Documented in memory file `project_contig_string_format.md`.

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

### Two new memory files to be aware of

`~/.claude/projects/-Users-jowillia-Documents-GitHub-receptor-resurfacing-pipeline/memory/`:
- `project_no_reversion_semantics.md` — read before touching negsteer outcomes
- `project_contig_string_format.md` — read before writing any contig example

Both indexed in `MEMORY.md`.

## How to verify the session's work locally before the next one

```bash
# On Mac, on phase4-impl:
git status                                                   # should be clean
git log --oneline 18718db..HEAD | wc -l                      # expect 19
python3 -m pytest -m local_unit -q \
    --ignore=tests/characterization/test_plots.py \
    --ignore=tests/characterization/helpers/tests/test_png_compare.py
# Expect: 350 passed, 7 skipped (gemmi-dependent), 272 deselected
```

If any of those checks fail, something has changed since end-of-session — investigate before continuing.

## When to update this file

Update at the END of every session, not as you work.  Replace the "What needs to happen FIRST" section with the new immediate next step.  Update the branch state if branches move.  Preserve "What happened in the last session" so the next reader can trace progress.
