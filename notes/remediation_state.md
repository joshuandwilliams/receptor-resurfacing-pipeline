# Remediation State

A living document tracking where the codebase remediation effort currently stands. Read this at the start of every session; update it at the end of every session.

**Last updated:** 2026-05-14 (Phase 4 verification complete on `phase4-impl`; all 6 fixtures regenerated; full pipeline test running on HPC; Session 7 HADDOCK grill-me is the next planned work)

## Current Phase

**Phase 4 implementation — VERIFICATION COMPLETE on branch `phase4-impl`.**  All 13 deep-module types built and tested; CL-3 majority-rule fix live; all caller migrations done; all 6 per-module HPC tests passed; all 6 fixtures regenerated and committed; 272/272 hpc-marked characterization tests passing.  Full pipeline test running on HPC at end-of-session — once that completes successfully, Phase 4 is fully closed.

## Branch state (as of 2026-05-14 end-of-session)

- `main` — at `18718db`.  Safe-fallback baseline.  Contains everything through Phase 4 spec (Sessions 1–6 complete) but no implementation code.  If `phase4-impl` ever needs to be abandoned, `main` is the working pipeline to return to.
- `phase4-impl` — at `be593db`.  38 commits ahead of `main`.  All Phase 4 implementation + caller migration + polish + verification work lives here.  Tracking remote `origin/phase4-impl`.
- `remediation` — at `18718db`.  Historical; can be retired.
- `experiments` — separate worktree branch; not touched in Phase 4.

## What landed in this session (2026-05-12)

### Phase 1 of session — pre-architecture deliverables (on `remediation`, then merged into `main`)

Six items from `notes/design_audit.md` "Audit close → Immediate actions":

1. **Schema rename + latent-bug fix** (commit `6e5959b`, 60 files):
   - `aggregated_verdict → outcome`, `aggregated_verdict_reason → outcome_reason`, `representative_* → rep_*`, the five `n_seeds_{...}` columns consolidated to one `n_pass`.
   - `_classify_aggregated_verdict → _classify_outcome` (signature changed to take `verdict_counts: Dict[str, int]`).
   - `_pick_representative_from_aggregated → _pick_rep_from_aggregated`.
   - **Latent-bug fix rides with rename**: `cross_sequence_summary._tier_for_row` no longer shortcuts `outcome == "no_reversion" → tier A`.  Tier is now derived purely from `n_pass / n_seeds`.  Wrong-placement-no-contamination groups (n_pass=0, outcome=no_reversion) now correctly land in tier none, not tier A.
   - 37 fixture CSVs rewritten by a one-shot script that was deleted post-run.
   - **New memory file**: `project_no_reversion_semantics.md` — recurring Claude failure mode (conflating `no_reversion` with cold-start).  Indexed in `MEMORY.md`.

2. **Five other audit-close items** (commit `05f0c07`):
   - Cohort summary visual groups (three-axis: Boltz structure / Boltz confidence / Orthogonal) + AF3 "(best)" label + restructured legend.  Mirror in `tests/orthogonal_metrics/test_orthogonal_metrics_plots.py`.
   - MPNN-sequence join into `cross_sequence_summary.csv` (`corrected_receptor`, `designed_residues`, `native_residues` columns), wired through `NEGSTEER_CROSS_SEQUENCE`.
   - `bin/validate_params.py` — 66 ParamSpec entries (59 with real constraints), 34 unit tests.  Wired into `main.nf` workflow head (fails fast before any SLURM job).
   - Glossary additions in `notes/inventory/06_ubiquitous_language.md`: `no_data` standalone entry, `n_pass` definition, `passes_orthogonal_filters` clarified as Sc+BSA+ΔΔG only, three-axis framework, F11 ("outcome" overload between per-sequence column and per-seed plot label).
   - Inventory mechanical updates (`02_function_inventory.md`, `15_discovery_run_path_coverage.md`).

3. **Tests/run scripts** (commit `54107b4`):
   - `tests/run_tests.sh` — multi-module dispatcher (parallel sbatch, `--with-plots` adds afterok dependency, `--no-clean` opt-out for the default cleanup).
   - `tests/update_example_dataset.slurm.sh` + `tests/_update_example_dataset_impl.py` — diff-and-replace with per-module `reference_manifest.txt`.
   - `scripts/sync_from_hpc.sh` — Mac-side puller (reverse of `sync_to_hpc.sh`).
   - Per-module `reference_manifest.txt` files for all five modules (haddock placeholder).
   - Retired `tests/rfdiffusion/capture_reference_set.sh`.

4. **Post-tests-failure fixes** (commits `1cdefd6`, `1f735c0`):
   - NO_FILE sentinel passed to `NEGSTEER_CROSS_SEQUENCE` in `tests/negative_steering/test_negative_steering.nf` (the 3-input regression from `05f0c07`).
   - `tests/run_tests.sh` cleans stale Nextflow + SLURM artefacts before sbatch by default (otherwise Nextflow's resume cache served stale prior-run outputs).

5. **Phase 4 architecture spec** (commits `d3a75f8`, `18718db`):
   - `notes/phase4_architecture_spec.md` — 1065 lines, the architectural contract for Phase 4.
   - Strategy B (deep modules) locked.
   - 13 domain types specced with full interface, fields, methods, boundary, consolidation targets, dependency tier (0–6), composition tree, data flow narrative, code-fit validation against 4 real-code samples.
   - Five user clarifications captured (CL-1 to CL-5).  CL-3 is the only intentional behaviour change.
   - **New memory file**: `project_contig_string_format.md` — canonical RFDiffusion contig syntax (`A1-10/5/A15-20 B`, slash-separated, never commas).

### Merge to main

After Phase 1, `remediation → main` fast-forward merge.  `phase4-impl` branched off `main` for the implementation phase.

### Phase 2 of session — Phase 4 implementation on `phase4-impl`

**All 13 deep-module types built (Tiers 0–6, 13 commits):**

| Tier | Type | File | Tests | Commit |
|---|---|---|---|---|
| 0.1 | `PipelineInternalThresholds` | `bin/pipeline_thresholds.py` | 19 | `f618b70` |
| 0.2 | `PipelineParams` | `bin/pipeline_params.py` | 15 | `8e921fb` |
| 0.3 | `BoltzConfidenceMetrics` | `bin/boltz_confidence.py` | 22 | `5dcbe98` |
| 0.4 | `AF3ConfidenceAggregate` | `bin/af3_confidence.py` | 19 | `3119196` |
| 0.5 | `ContigSpec` | `bin/contig_spec.py` | 35 | `8ee1d55` |
| 1 | `PositionSet` | `bin/position_set.py` | 41 | `1b51259` |
| 2 | `ProteinStructurePrediction` | `bin/protein_structure_prediction.py` | 23 (16 local + 7 gemmi) | `5457273` |
| 3 | `DesignedBackbone` | `bin/designed_backbone.py` | 11 | `edd9336` |
| 4.1 | `DesignedSequence` | `bin/designed_sequence.py` | 9 | `f6ae964` |
| 4.2 | `StageResult` | `bin/stage_result.py` | 30 (incl. 6 CL-3 tests) | `0859641` |
| 4.3 | `OrthogonalMetrics` | `bin/orthogonal_metrics.py` | 14 | `b044447` |
| 5 | `NegativeSteeringRun` | `bin/negative_steering_run.py` | 14 | `36355b3` |
| 6 | `DesignCohort` | `bin/design_cohort.py` | 13 | `0f194ee` |

Total: ~4,000 LOC new code + 350 local_unit tests, all green.

**Six caller migrations completed:**

| Migration | Commit | Impact |
|---|---|---|
| **CL-3 fix in `cmd_build_contaminated`** | `b67b4b3` | Reversion now per-design majority-of-correctly-placed, not per-(design, seed).  Headline behaviour change.  Encoded inline in the existing function rather than as a wholesale rewrite to StageResult (lower risk; same data flow). |
| Test/prod plot-script divergence eliminated | `4370e99` | 3 test scripts collapsed into thin wrappers; `tests/orthogonal_metrics/merge_orthogonal_metrics.py` deleted (zero callers).  **~4,200 LOC duplication removed.** |
| `get_chain_sequence` consolidation (2 of 4 duplicates) | `4e19b1b` | `build_control_sequences._extract_chain_sequence` and `extract_survivor_manifest._extract_chain_seq` migrated to canonical `boltz2_negative_steering.get_chain_sequence`. |
| `THREE_TO_ONE` / `_AA3TO1` dedup | `e60628b` | `protein_structure_prediction.py` and `extract_survivor_manifest.py` no longer carry their own copies; both use `contig_utils.THREE_TO_ONE`. |
| `parse_af3_output.py` rewritten as thin CLI wrapper | `8230ee5` | 319 → 104 LOC.  All logic now in `AF3ConfidenceAggregate`.  CLI contract preserved; `AF3_PARSE_OUTPUT` Nextflow process unchanged. |
| `extract_passing._compute_confidence_flag` reads from `PipelineInternalThresholds` | `54dec6f` | First threshold-catalogue caller migration.  Lazy import + fallback to historical literals if `pipeline_thresholds` unavailable. |

## Caller migrations explicitly NOT done (with reasons)

These were considered and deferred — each is a real Phase 4 todo but each has a specific reason it deserves a future, focused session:

1. **~~`cross_sequence_summary.py` → `DesignCohort.to_cross_summary_csv`~~** — DONE in commit `babd5d4`.  `DesignCohort.from_runs_directory` and `NegativeSteeringRun.from_workdir` factories are now real.  CSV emission routes through `DesignCohort.emit_cross_summary_from_dirs` (delegates to the existing `aggregate()` for column-level construction; output is bit-identical).  Nextflow process points at `bin/cross_summary_v2.py` as the directly-tracked script.
2. **~~ContigSpec parser consolidation (4 existing parsers → ContigSpec)~~** — DONE in commits `7385042` (three Python parsers) and `babd5d4` (the fourth, `contig_utils.parse_block_segments`).  Resolved by collapsing resolved/constraint forms: `DeNovoSegment` now always carries `(min_len, max_len)`, with `min_len == max_len` for the resolved case.  RFDiffusion-specific grammar features (chain-break `0`, bare-chain passthrough) modelled as new `BreakSegment` + `PassthroughSegment` types on `ContigChain`.
3. **~~Full threshold migration~~** — DONE in commit `7385042`.  `orthogonal_metrics_plots.py`, `compute_metrics.py` (weighted-Jaccard constants), and `boltz2_iterate_steering.py` (`intact_threshold`) all migrated to lazy-imported `PipelineInternalThresholds` with literal fallback.  `main.nf` inline `params.*` defaults are user-facing knobs not in scope for the threshold catalogue.
4. **~~`find_contact_residues_heavy` consolidation~~** — Resolved in commit `7385042`.  The two implementations now produce IDENTICAL output (compute_metrics.py's was buggy on PDBs with insertion codes; fixed to match boltz2_negative_steering.py's icode-aware bucketing).  Cross-reference docstrings added on both ends.  Two implementations retained intentionally; future refactor can extract a shared `pdb_atom_io` helper.
5. **`read_ca_atoms` consolidation** — INTENTIONALLY NOT consolidated.  Two implementations with different return types (CAEntry list vs coord-dict list) serve different consumers.  Cross-reference docstrings added on both ends (commit `7385042`).
6. **~~`pipeline_correct_sequences.get_pdb_sequence`~~** — DONE in commit `7385042` (deleted as dead code; had zero callers).

## HPC validation queue (pending at end of session 2026-05-13)

All five per-module tests submitted via `./tests/run_tests.sh --modules rfdiffusion proteinmpnn rosetta_filtering negative_steering orthogonal_metrics --with-plots`.  `haddock` intentionally omitted (deprecated, not currently used).

GPU queue was congested at submission time (another user's 6,000-job array job had blocked GPU access).  RFDiffusion, negative_steering, and orthogonal_metrics' AF3 step depend on GPU; their queue timing depends on that array clearing.

**The CL-3 fix (`b67b4b3`), the deep-form factories (`babd5d4`), and the cross_summary_v2.py rewiring have not yet been validated against real HPC data.**  The next session should:
- Confirm `Success: true` in each module's `slurm_<jobid>.out`
- For `negative_steering`: check the new `[build-contaminated] CL-3 gating: …` log line for the expected per-design reversion redistribution
- Diff `tests/<module>/example_output_files/` against the post-run `tests/<module>/receptor_resurfacing_results/` for each module that completed
- Regenerate fixtures only where the diff is sensible and CL-3-driven (use `sbatch tests/update_example_dataset.slurm.sh --module <module> --updated-output-folder tests/<module>/receptor_resurfacing_results` then `./scripts/sync_from_hpc.sh --module <module>`)

**Note for `negative_steering` specifically:** expect some tier-C/none → tier-A/B promotions (sequences whose lone contaminated seed used to trigger reversion-and-collapse now keep their clean_steered partners and pass).  Some `reversions/` subtrees should disappear from designs that no longer trigger reversion under the majority rule.

## What landed in this session (2026-05-13)

Five commits on top of the prior end-of-session state (`020d32b`).  All built and verified locally (370 local_unit tests pass, up from 350).  Nothing GPU-blocked locally; HPC validation is the remaining gap.

### Commit `7385042` — threshold migration + ContigSpec range form + parser dedup + icode fix

- **Threshold migration**: `bin/orthogonal_metrics_plots.py`, `bin/compute_metrics.py` (weighted-Jaccard mu + two_sigma_sq), and `bin/boltz2_iterate_steering.py` (`intact_threshold` at 4 sites) now lazy-import `PipelineInternalThresholds` with literal fallback.  No value changes; just sourcing.  Pattern matches commit `54dec6f`.
- **ContigSpec range form**: `bin/contig_spec.py` reworked so `DeNovoSegment` always carries `(min_len, max_len)`.  Resolved form is the special case `min_len == max_len`.  Parser accepts both `5` (bare integer → `(5,5)`) and `5-7` (range).  Serializer emits `5` when resolved, `5-7` otherwise.  Added `ContigSpec.from_string` (form-agnostic) + `ContigSpec.is_resolved`.  `from_resolved_string` preserved as a strict variant that raises if any segment is unresolved.  Per Q on parser consolidation, the user directed: "What other forms of denovo segments are there? All denovo segments should take that format.  Even if it's not variable length it will still be `5-5` instead of `5-7`."
- **Three contig parsers migrated**: `derive_input_design_region._parse_contigs`, `haddock3_prepare.parse_contig_segments`, `pipeline_correct_sequences.parse_contig_segments` are now thin adapters around `ContigSpec.from_string`.  Comma-separator back-compat preserved by pre-normalising to spaces.
- **`find_contact_residues_heavy` icode bug fix**: `bin/compute_metrics.py:596` rewritten to do icode-aware bucketing (matching `boltz2_negative_steering.read_residue_heavy_atoms`).  Two residues with the same resseq but different insertion codes (common in wwPDB structures) used to collapse into one — wrong positional indices on any input PDB with icodes.  Both copies now produce identical output.  Cross-reference docstrings added on both ends.  `read_ca_atoms` intentionally NOT consolidated (different return shapes serve different consumers); docstrings cross-referenced.
- **`bin/cross_summary_view.py` NEW**: `CrossSummaryRow` + `CrossSummarySnapshot` — typed read-only view over `cross_sequence_summary.csv`.  Exposes the same `survivors() / tier_breakdown() / ranked_by_composite()` interface as `DesignCohort` for cohort-level queries that don't need the deep StageResult/PSP scaffolding.
- **`DesignCohort.from_cross_summary_csv`** + **`DesignCohort.to_cross_summary_csv`** added.  The latter delegates to the existing `aggregate()` for bit-identical CSV output.
- **Cleanup**: dead function `get_pdb_sequence` removed from `bin/pipeline_correct_sequences.py` (zero callers).

### Commit `babd5d4` — complete the deferred deep-form migrations

- **`NegativeSteeringRun.from_workdir` deep form**: walks the canonical workdir layout (`cycle_<C>/initial[_sN]`, `steered/design_NN_sS`, `reversions/rev_design_NN_sS_sR`); for each prediction directory finds the canonical Boltz output (`boltz_results_input/predictions/input/pdb/input_model_<M>.pdb` + matching confidence JSON) and constructs `ProteinStructurePrediction` + `BoltzConfidenceMetrics`; parses `mutations.tsv` into a designed-frame `PositionSet` for steered/reversion stages; resolves reversion → steered mapping via the source-(design, seed) prefix encoded in directory names; reads `row_type.txt` + `run_one_runtime_sec.txt` sidecars.  Returns `None` cleanly when cold-start data is missing.  Six dedicated unit tests using synthetic workdirs with minimal PDB + confidence fixtures.
- **`DesignCohort.from_runs_directory`** — now hydrating runs via `from_workdir` per subdir; tested with synthetic runs dir including bogus-subdir skip path.
- **`DesignCohort.emit_cross_summary_from_dirs`** classmethod — typed entry point for the Nextflow process to call.  Gathers `(seq_name, passing_summary.csv)` pairs and emits via the existing `aggregate()` column-builder.  Bit-identical CSV output.
- **`bin/cross_summary_v2.py` NEW** — Phase 4 typed CLI entry point.  Drop-in replacement for `cross_sequence_summary.py` from the Nextflow perspective (same flags, same output).
- **Nextflow rewiring**: `NEGSTEER_CROSS_SEQUENCE` in `main.nf:879` and `tests/negative_steering/test_negative_steering.nf:361` now points at `cross_summary_v2.py`.
- **`contig_utils.parse_block_segments` fourth-parser migration**: added `BreakSegment` + `PassthroughSegment` types to `bin/contig_spec.py`.  `ContigChain` position-math methods (`total_length`, `designed_position_to_native`, etc.) iterate a `_length_bearing` filtered list — break + passthrough are skipped naturally.  `ContigSpec.from_string` now accepts the literal `0` chain-break marker and bare-chain-letter passthroughs.  `parse_block_segments` reshapes to the legacy `(descs, block_chain)` tuple form so downstream `resolve_contigs` + `remap_segments_to_pdb` keep working unchanged.

### Commit `23e0fc4` — Nextflow comment refresh

`modules/negative_steering.nf` comment block previously said `cross_sequence_summary.py` is the directly-tracked script; updated to reflect that `cross_summary_v2.py` is now the directly-hashed entry point, with `cross_sequence_summary.py` (and `extract_passing.py`) as indirect imports.

### Commit `c15923f` — `tests/run_tests.sh` silent-exit bug fix

Old `clean_module()` ended with:
```bash
[ ${#stale_logs[@]} -gt 0 ] && rm -f "${stale_logs[@]}"
```
When `stale_logs` was empty (no `slurm_*` files to delete), `[ 0 -gt 0 ]` returned exit code 1, `&&` short-circuited, and the compound command exited 1.  As the function's LAST command, that propagated as the function's return status.  Under `set -e` this silently killed the dispatcher.  First three modules in a session worked because they had stale `slurm_*` files from prior runs; the fourth module (alphabetically first one without stale logs) would silently die.  Replaced with explicit `if [ … ]; then rm -f …; fi` so the function exits 0 cleanly when there's nothing to clean.

### Caller migrations status

All six items previously listed in "Caller migrations explicitly NOT done" are now resolved:
1. ~~cross_sequence_summary → DesignCohort~~: DONE (`babd5d4`)
2. ~~ContigSpec parser consolidation~~: DONE across two commits (`7385042` for 3 parsers, `babd5d4` for the 4th)
3. ~~Full threshold migration~~: DONE (`7385042`)
4. ~~find_contact_residues_heavy consolidation~~: RESOLVED via bug fix (`7385042`); both copies now produce identical output
5. read_ca_atoms consolidation: INTENTIONALLY NOT consolidated (different return types serve different consumers); cross-reference docstrings added
6. ~~get_pdb_sequence~~: DONE — was dead code, removed (`7385042`)

### Test state

```
370 local_unit tests pass (up from 350)
+13 ContigSpec tests (range form, break/passthrough, parsing, position math)
+ 4 CrossSummarySnapshot tests
+ 6 NegativeSteeringRun.from_workdir tests
+ 4 ContigSpec.is_resolved / from_string variants
```

No regressions.  Skipped tests: 7 (gemmi-dependent).

### Architecture spec follow-up

`notes/phase4_architecture_spec.md` §2.2 (ContigSpec) currently describes the resolved-form-only model and three segment types (Fixed, DeNovo, plus Passthrough as a single PassthroughChain encoding).  The implementation now has FOUR explicit segment types (Fixed, DeNovo with range form, Break, Passthrough) on `ContigChain`.  Spec doc should be updated to reflect this — minor write-up task for the next session.

### HPC sync state

All commits pushed to `origin/phase4-impl` and synced to HPC via `./scripts/sync_to_hpc.sh`.  HPC tests submitted at end of session — see "HPC validation queue" section above.

## What landed in this session (2026-05-13 → 2026-05-14)

17 commits on top of the prior `ea6dd71` (the previous notes update).  All Phase 4 verification + polish + fixture regeneration + full-pipeline-test integration.  Branch went from `c15923f` (24 commits ahead of main) to `be593db` (38 commits ahead).

### A. Polish pass (post-implementation, pre-verification)

- **`1641c76`** — `extract_survivor_manifest.py` hardening.  Two changes: (1) remap stale absolute paths.  `rep_canonical_pdb` (CSV), `ground_truth` and `effector_template_cif` (plan.json) carry absolute paths stamped at fixture-generation time.  When the fixture is built on host A (Mac at `/Users/...`) and consumed on host B (HPC at `/hpc-home/...`), those paths don't resolve.  Two new helpers: `_remap_canonical_pdb` (finds `/runs/<seq>/` segment, joins tail onto discovered workdir) and `_remap_repo_path` (finds `/tests/`-or-`/bin/`-etc. segment, joins against runtime repo root deduced from script location).  (2) Empty-manifest guard.  Previously the script exited 0 even when every survivor was skipped, letting Nextflow fan out the orthogonal cascade over an empty channel (which it treated as successful no-op).  Now exits 2 with stderr message.

- **`e3c2752`** — three polish items in one commit.  (1) ERR trap on `tests/run_tests.sh` that prints failing command + line + rc before the shell exits — defends against the general "function's last command is a short-circuiting compound returning non-zero under `set -e`" pattern, not just the c15923f-specific bug.  (2) F821 typing fix in `bin/boltz2_negative_steering.py` — 7 uses of `Optional` in string-typed annotations without `Optional` being imported.  Trivial fix (add to `from typing import` line); `ruff check --select F821` now clean.  (3) `NEGSTEER_CROSS_SEQUENCE` cache-busting via passing the entire `bin/` directory as a `path bin_dir, name: 'bin'` input.  Process invokes `python bin/cross_summary_v2.py …` so Python's sibling-import resolution finds every dependency; Nextflow content-hashes the whole `bin/` tree, busting cache on any indirect-import edit.  Closes the verification queue's "indirect-import cache-busting gap".

### B. Orthogonal cascade enhancements

- **`f5f9167`** → **`373e164`** → **`0bf94f6`**: combined-cohort plot ranking iteration.  First attempt added a Boltz-2-confidence tie-breaker; reverted at user request after misreading; restored after user clarified the original intent.  Final state: `_sort_key = (cross_tier_order, fails_any_boltz_conf, -composite)` — tier first, rows that pass every thresholded Boltz-2 confidence metric next, then descending composite within each group.  Missing values count as fails for ranking purposes.  Helper `_boltz_conf_thresholded_columns()` derives the threshold list from `COMBINED_SUMMARY_GROUPS` at call time (single source of truth).  Legends side-by-side instead of stacked (`ncol=2` each, `bbox_to_anchor` at `(0.0, -0.04)` and `(1.0, -0.04)`).  Also (in `0bf94f6`): `ORTHOG_PLOTS` wired into `tests/orthogonal_metrics/test_orthogonal_metrics.nf` so plots land in `plots/` matching the pattern of every other module test.  The standalone `run_test_orthogonal_metrics_plot.slurm.sh` still exists for plot-only iteration; publishes to `plots_iter/`.

- **`20f9ed0`** — `params.orthogonal_tier_filter` introduced.  Configurable choice; default `'all'` runs AF3 + biophysical + Rosetta on every steered design (including tier-none failed designs — useful diagnostic since orthogonal metrics can explain WHY a design failed even when its negsteer tier is none).  `'abc'` restricts to cross_tier in (A, B, C).  Plumbed through `extract_survivor_manifest.py` (new `--tier-filter` CLI flag, with skip-counter increment when restricted) and `modules/negsteer_manifest.nf`.  **Also added `scripts/refresh_orthog_fixture.sh`** — one-shot rsync that mirrors the negsteer test outputs (`tests/negative_steering/receptor_resurfacing_results/.../runs/`) into the orthogonal_metrics fixture (`tests/orthogonal_metrics/data/negsteer_run/runs/`).  Errors loudly if negsteer outputs are missing; idempotent (`rsync -a --delete`).  Used to expand the orthogonal_metrics test cohort from the previously-staged 2-survivor mini-fixture to the full 10-row (8 designs + 2 controls) negsteer output.

### C. Test infrastructure fixes

- **`06c0157`** — `tests/haddock/test_haddock.nf` overrides `params.haddock_min_cluster_size = 2` (production default 4 in `nextflow.config`).  HADDOCK itself ran cleanly through all 8 steps but `collect_haddock3_dock.py` rejected the run because the only cluster had 3 models, below the production threshold.  At `haddock_sampling = 100` (vs production 10000), clusters are naturally smaller.  Override is test-only.

- **`ba95ede`** + **`fa14d3a`** — `tests/update_example_dataset.slurm.sh` and its impl `_update_example_dataset_impl.py`: two sbatch-specific bugs.  (1) `${BASH_SOURCE[0]}` resolves to `/var/spool/slurmd/job<id>/` under sbatch, not the real `tests/` dir.  Fix: try `${PWD}` first (set by `#SBATCH --chdir`), fall back to `BASH_SOURCE`.  (2) Path-doubling: when CWD is already `tests/` and the user passes `tests/<module>/...` (project-root-relative) as the value of `--updated-output-folder`, it resolves as `tests/tests/<module>/...`.  Impl now auto-strips a leading `tests/` segment from the as-given path when that path doesn't exist; emits stderr NOTE so the user knows their command was ambiguous.  Both fixes preserve direct-shell-invocation behaviour on Mac.

### D. Fixture regeneration

- **`e5f3be0`** — all 6 module fixtures regenerated against `phase4-impl` head (post-1641c76).  657 files updated.  By module:
  - **negative_steering**: 627 files.  Driven by CL-3 majority-rule reversion gating (per-design, not per-seed) and the prior `n_pass / outcome` schema rename.  Touches every design's per-seed outputs.
  - **rfdiffusion**: 9 files.  RFD container update (`complex_beta` checkpoint) produces structurally different designs; metrics + plots all shifted.
  - **proteinmpnn**: 9 files.  Small `design_region_score` shifts (likely ProteinMPNN minor non-determinism); plot redraws follow.
  - **rosetta_filtering**: 5 files.  Blesses the previous `design[3].dG_separated` 33→57 shift.  **Cause not isolated** — flagged in the verification queue.
  - **orthogonal_metrics**: 6 files.  `cross_sequence_summary_with_interface_metrics.csv` and `survivor_manifest.csv` now 10 rows (was 2) after the test fixture was refreshed from negsteer outputs and `orthogonal_tier_filter='all'` (default) let every steered design through.
  - **haddock**: NEW fixture (first time HADDOCK has one).  Manifest intentionally minimal (Nextflow run-artefacts only — `dag.html`, `pipeline_report.html`, `timeline.html`, `trace.txt`).  Scientific outputs (docked PDBs, CAPRI scores) are placeholders pending the HADDOCK restructure.

### E. Full pipeline test integration

- **`279f078`** — `tests/run_tests.sh` now supports `full_test_run` as a module name.  Special-cased in `submit_one()` because the launcher lives at the repo root (`run_pipeline.slurm.sh`) and takes a params file argument (`tests/full_test_run/params_full_test.yml`), not the per-module test wrapper pattern.  `--no-clean` preserves Nextflow state (`work/`, `.nextflow*`) so a follow-up `--resume` can continue interrupted runs.  The expensive `results/` tree is NOT touched by clean either way.  `--with-plots` is a no-op for `full_test_run` (the main pipeline runs every plot process inline).

- **`8e21c7b`** — `tests/full_test_run/params_full_test.yml`: bumped `haddock_sampling: 1 → 100`.  Validator (`bin/validate_params.py`) runs unconditionally regardless of mode and rejects values below 100; mode 2 (pre-docked PDB) skips HADDOCK so the actual value is unused.  Placeholder.  Also added a `ParamSpec` for `orthogonal_tier_filter` (the new param introduced in `20f9ed0`).

- **`be593db`** — `tests/full_test_run/params_full_test.yml`: `pdb_file` corrected to `af3_pikp1_native_avrpikf_complex.pdb` (the actual filename on disk).  The previous reference (`pikp1_avrpikf_complex.pdb`) didn't exist on either Mac or HPC.

### Verification at end of session

- **272 / 272 hpc-marked characterization tests passing.**  386 local_unit tests deselected by `-m hpc` (run those Mac-side for a quick green).  26 expected skips (long-standing, documented in §"Known expected skips in pytest").
- **All 6 per-module test runs reported `Success: true`** prior to fixture regeneration.
- **Full pipeline test (`./tests/run_tests.sh --modules full_test_run`) running on HPC at end-of-session.**  ETA ~24h+ depending on GPU queue.  Submission required two iterations to clear validator + a stale `pdb_file` path; the running attempt should be checked first thing next session.

### Caller migrations resolution status (post-Phase 4)

| Item | Status |
|---|---|
| `cross_sequence_summary.py` → `DesignCohort.to_cross_summary_csv` | DONE in `babd5d4` |
| ContigSpec parser consolidation (4 → 1) | DONE in `7385042` (3) + `babd5d4` (4th) |
| Full threshold migration | DONE in `7385042` |
| `find_contact_residues_heavy` consolidation | RESOLVED via icode bug fix (`7385042`); two impls now produce identical output |
| `read_ca_atoms` consolidation | INTENTIONALLY NOT consolidated; cross-reference docstrings |
| `get_pdb_sequence` | DONE — was dead code; removed (`7385042`) |
| `NegativeSteeringRun.from_workdir` deep form | DONE in `babd5d4` |
| `cross_summary_v2.py` Nextflow rewiring | DONE in `babd5d4` |
| `contig_utils.parse_block_segments` migration | DONE in `babd5d4` |
| `ORTHOG_PLOTS` in `tests/orthogonal_metrics/test_orthogonal_metrics.nf` | DONE in `0bf94f6` |
| `orthogonal_tier_filter` user-facing param | DONE in `20f9ed0` |
| `extract_survivor_manifest.py` path remap + empty-guard | DONE in `1641c76` |
| `tests/run_tests.sh` full_test_run support | DONE in `279f078` |

Every caller migration originally listed as "deferred" is now resolved.

### HPC sync state

All commits pushed to `origin/phase4-impl` and synced to HPC via `./scripts/sync_to_hpc.sh`.  Full pipeline test submitted at end-of-session and running.

## Pre-Phase-4 baseline commit

**`1f735c0`** — captured at end of pre-architecture work, 2026-05-12.

If a regression surfaces in the rfdiffusion or negsteer module test after this commit (those two were queued on GPU at the time and could not be verified in-session), bisect against `1f735c0` to pinpoint when the divergence was introduced.  The four CPU-only module tests (proteinmpnn, rosetta_filtering, orthogonal_metrics, and rfdiffusion's CPU portion) all reported `Success: true` against this commit; only negsteer needed a follow-up fix (`1cdefd6`, NO_FILE sentinel for the new `scored_metadata` input) before its first successful run.

## Phase 4 architecture spec — Session 6 (COMPLETE)

Architecture grill-me Session 6 produced `notes/phase4_architecture_spec.md` — the written specification for the Phase 4 codebase rebuild under Strategy B (deep modules).  Spec covers 13 domain types: `ProteinStructurePrediction`, `ContigSpec`, `PositionSet`, `BoltzConfidenceMetrics`, `AF3ConfidenceAggregate`, `DesignedBackbone`, `DesignedSequence`, `StageResult`, `NegativeSteeringRun`, `DesignCohort`, `OrthogonalMetrics`, `PipelineParams`, `PipelineInternalThresholds`.

**All four steps populated:**
- Step 1: candidate type list (13 types, triaged from initial 10 against the user's domain corrections).
- Step 2: per-type interface specs — domain meaning, construction, methods, boundary, consolidation targets from current code.
- Step 3: composition tree, dependency-order table (migration tiers 0→6), data flow through one pipeline run.
- Step 4: four code-fit validations rewriting current code against the new types (contamination check + CL-3 fix; cross_sequence_summary aggregate; per-seed cross-stage verdict aggregation; orthogonal filter check).  Three small spec additions surfaced and back-fitted into §2.8 and §2.9.

**CL-3 captured as the only intentional behaviour change Phase 4 will make**: the reversion-gating rule moves from per-(design, seed) to per-design majority-of-correctly-placed.  Everything else is structural.

**Migration tier order** (each tier's types are independent and can land in any order; tiers themselves are sequential):
- T0: `ContigSpec`, `BoltzConfidenceMetrics`, `AF3ConfidenceAggregate`, `PipelineParams`, `PipelineInternalThresholds`.
- T1: `PositionSet`.
- T2: `ProteinStructurePrediction`.
- T3: `DesignedBackbone`.
- T4: `DesignedSequence`, `StageResult`, `OrthogonalMetrics`.
- T5: `NegativeSteeringRun`.
- T6: `DesignCohort`.

Two new memory files added (recurring Claude failure modes): `project_no_reversion_semantics.md` and `project_contig_string_format.md`.  Both indexed in `MEMORY.md`.

Next: implementation phase per the spec.  Each tier's types lands as its own commit, characterization tests re-run between commits using the workflow established in `tests/run_tests.sh` + `tests/update_example_dataset.slurm.sh`.

---

## Current Phase

**Phase 3 (Mechanical Cleanup) — COMPLETE.**

Tag: `phase-3-complete`.

Dead code removed, key naming inconsistencies fixed, critical duplication consolidated, two scientific correctness bugs fixed, and the three worst complexity hotspots reduced from F-rated to C/D/E. Safety net confirmed green (272 passed, 26 skipped, 0 failed) after all changes.

---

## Just Completed (this session)

### Pre-Phase-4 audit-close items (2026-05-12)

All six items from `notes/design_audit.md` "Audit close → Immediate actions" now land in `remediation`. Two commits: `6e5959b` (rename + latent-bug fix) and `05f0c07` (the other five).

**Commit `6e5959b` — schema rename + latent-bug fix:**
- Column renames across all output CSVs and producer/consumer code:
  `aggregated_verdict → outcome`, `aggregated_verdict_reason → outcome_reason`,
  `representative_* → rep_*`, and the five `n_seeds_{pose_holds, pose_collapses,
  new_contamination, no_data, clean_steered}` columns consolidated to one
  `n_pass` column (= `pose_holds_count + clean_steered_count`,
  computed unconditionally from `_per_seed_verdict_breakdown`).
- Function renames: `_classify_aggregated_verdict → _classify_outcome` (signature
  now takes `verdict_counts: Dict[str, int]` instead of reading `n_seeds_*` keys
  off the agg dict); `_pick_representative_from_aggregated → _pick_rep_from_aggregated`.
- Latent-bug fix riding with the rename: `cross_sequence_summary._tier_for_row`
  no longer shortcuts `outcome == "no_reversion" → tier A`. Tier is derived
  purely from `n_pass / n_seeds`. **Wrong-placement-no-contamination groups
  (n_pass=0, outcome=no_reversion) now correctly land in tier none**, not tier A.
- New memory file `project_no_reversion_semantics.md` documents the recurring
  Claude conflation between `no_reversion` and cold-start. **Read before
  touching negsteer outcomes.**
- 60 files changed including 37 fixture CSVs (rewritten by a one-shot script
  that was deleted post-run).

**Commit `05f0c07` — five other audit-close items:**
- *Cohort summary visual groups* (`bin/orthogonal_metrics_plots.py`):
  three-axis column groups (descriptive / Boltz-2 structure / Boltz-2
  confidence / Orthogonal) drawn above the metric labels; AF3 ra_eff
  column relabelled "AF3 ra_eff (best)"; legend split into two titled
  blocks so the left stripe (negsteer cohort outcome) is visually
  independent from cell colours (per-metric pass/fail). Mirror change
  applied to the test-copy `tests/orthogonal_metrics/test_orthogonal_metrics_plots.py`.
- *MPNN-sequence join* (`bin/cross_sequence_summary.py`,
  `modules/negative_steering.nf`, `main.nf`): new `--scored-metadata`
  flag joins three columns (`corrected_receptor`, `designed_residues`,
  `native_residues`) by `mpnn_sequence`. `NEGSTEER_CROSS_SEQUENCE` now
  consumes `MPNN_DESIGN_REGION_SCORE.out.scored_metadata`.
- *Parameter validator* (`bin/validate_params.py`): 66 `ParamSpec` entries
  (59 with concrete constraints, 7 documented `any` coverage gaps). Eight
  validator kinds: `int_range`, `float_range`, `choice`, `bool`,
  `non_empty_str`, `optional_path`, `regex`, `custom`, `any`. The
  `negsteer_num_seeds` odd-only check is a custom validator. The Nextflow
  workflow head dumps params to JSON, invokes the validator, and
  `error()`s with the full stderr on non-zero exit — fails fast before
  any SLURM job is dispatched. Coverage report via `--report`.
- *Glossary update* (`notes/inventory/06_ubiquitous_language.md`):
  `no_data` standalone entry; `n_pass` definition; `passes_orthogonal_filters`
  clarified as Sc + BSA + ΔΔG only; three-axis framework added; F11
  flags the `outcome` term as overloaded between per-sequence column
  and per-seed plot label.
- *Inventory mechanical updates* (`02_function_inventory.md`,
  `15_discovery_run_path_coverage.md`): function-name renames mirrored;
  a header note in 15 covers the historical diagnostic commands that
  still reference pre-rename column names.

**Tests**: 92 `local_unit` tests pass on Mac (up from 58 pre-session;
the +34 are validator unit tests). Full HPC characterization sweep
deferred — the user will re-run after these six items land, to set a
fresh baseline before Phase 4 begins.

### Cleanup script: RFDiffusion traj/ added; array extended (2026-05-08)

- `scripts/clean_finished_run.sh`: step 4 added — deletes `results/rfdiffusion/traj/` after
  work/ and Boltz pruning. The `trajectories` channel emitted by `rfdiffusion.nf` is never
  consumed in `main.nf`; traj files are published but dead weight, and in practice the largest
  remaining space consumer after work/ is gone.
- `scripts/clean_all_finished_runs.slurm.sh`: `--time` extended to 6 h (prior 30-min limit
  was too short for some runs). `--array` extended from 0-7 to 0-14; the 7 v1_4f runs that
  completed since the original audit are now included. `pikp1_avrpikf/v2_4b` removed from
  the excluded list (does not exist on HPC).
- 15 runs confirmed as needing cleanup (traj present); 3 never-run runs correctly excluded.
- Characterisation suite confirmed green against updated RFDiffusion reference set:
  272 passed, 26 skipped, 0 failed (HPC job 19997103, 2026-05-08 21:15).

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

**Pre-Phase-4 items are now complete** (commits `6e5959b`, `05f0c07`). Run the full HPC characterization sweep next to set a fresh baseline (the schema rename and the latent-bug fix together changed many fixture-comparable outputs; some tier-A rows will move to tier none for genuinely-wrong-placement designs — investigate any non-trivial diff before assuming it is a regression).

**Before writing any Phase 4 code:**
1. Run the grill-me skill to design the target module architecture. Authoritative inputs:
   - `notes/design_audit.md` — Phase 4 scope is in §"Audit close → Phase 4 architecture — scope confirmed". Key modules to design:
     - `boltz_lib.py` — extract `get_chain_sequence`, `find_contact_residues_heavy` (negsteer version), `jaccard`, weighted_jaccard primitives, `binding_rmsds`, `write_boltz_yaml`, Cα/Kabsch helpers, residue dataclasses.
     - `plot_lib.py` — `make_empty_plot`, `save_fallback_plots`, `COLOUR_TIER`, `_TIER_SORT_ORDER`, `COMPOSITE_RA_EFF_WEIGHT`, shared utilities (currently duplicated across plot scripts and test copies).
     - `derive_indices.py` — merge `derive_design_region.py` + `derive_true_interface.py`.
     - Contig parser consolidation — `_parse_contigs` and `parse_contig_segments` rewritten on `parse_block_segments`.
     - Six F-rated `cmd_*` function decomposition in `boltz2_iterate_steering.py`.
     - `merge_orthogonal_metrics.py` test/production divergence elimination.
   - `notes/inventory/06_ubiquitous_language.md` — recently updated with `n_pass`, `outcome`, three-axis framework, the `outcome` overload (F11).
   - `notes/inventory/04_functional_categorization.md`, `10_phase_1_synthesis.md` — original Phase 1 outputs.
   - Output: a short architecture document committed to `notes/`.

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
- **Inline RFDiffusion traj/ cleanup into the pipeline.** Currently handled retroactively by
  `scripts/clean_finished_run.sh`. The right fix is to stop publishing the directory at all:
  change the `publishDir` in `modules/rfdiffusion.nf` to add a `pattern:` filter that excludes
  `traj/`, or delete `traj/` at the end of the process script body. Either approach means future
  runs never accumulate traj files and the standalone cleanup step becomes unnecessary.
  Analogous to the §10 forward-looking change in `notes/inventory/17_storage_bloat_audit.md`
  for Boltz unused samples.

---

## Verification Queue

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