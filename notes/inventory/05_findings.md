# 05 — Findings

A running list of discoveries that didn't fit cleanly into the other documents. Roughly grouped by category. Severity is a rough triage label, not a formal scale.

---

## A. Notes-vs-code discrepancies

### A1. The "two codebases" framing in `notes/project_map.md` is now historical

**Severity**: low (informational, but important to internalise).

`notes/project_map.md` (dated 15 Apr) describes "two related but distinct codebases": the receptor-resurfacing Nextflow pipeline and the standalone negative-steering SLURM/Python pipeline, with the latter "being prepared for integration into the Nextflow pipeline." That integration has happened — the negsteer chain is now `modules/negative_steering.nf` + `modules/negsteer_*.nf` (8 modules total) running inside one workflow. A useful next step would be to update or replace `project_map.md` with an updated overview, but that's a documentation task beyond the scope of this inventory.

### A2. `project_map.md` lists 14 negsteer files that no longer exist

The project_map's "negative-steering pipeline (standalone)" section enumerates files that are absent from the current tree:

- `submit_boltz2_negative_steering.sh`, `submit_boltz2_iterate_steering.sh`, `submit_steering_experiments.sh` — SLURM submitters replaced by Nextflow-managed `bin/negative_steering_run_one.sh`.
- `compare_versions.py` — cross-version validator, never ported.
- `aggregate_results.sh` — multi-experiment aggregator; the equivalent functionality is now in `boltz2_iterate_steering.py:cmd_aggregate`.
- `generate_waypoints.py`, `generate_waypoints_7QPX_chainB.py`, `combine_waypoints.py` — diagnostic-experiment scripts from the 7QPX session; not in current tree.
- `all_results_multicycle_with_metrics.csv` — example output file; not in current tree.

The user's prompt explicitly mentioned `submit_boltz2_iterate_steering.sh` as already-known-missing. The above list is the complete set of files referenced in `project_map.md` that are absent. None appear to be missed integrations — they're all artefacts of the standalone-pipeline era that didn't make the cut into the integrated pipeline.

### A3. `project_map.md` lists modules that no longer exist

- `aggregate.nf`, `boltz2.nf`, `msa.nf` are referenced as current modules in `project_map.md §1.2`. None exist on disk. Their functionality has been absorbed into the negsteer chain modules (`negative_steering.nf`, `negsteer_*.nf`, `negsteer_orthogonal_metrics.nf` for the aggregate role). This is consistent with the project's intended end-state.

### A4. The `modules/negsteer_orthogonal_metrics.nf` notes claim AF3 demoted, but production code still gates on AF3

**Severity**: medium (semantic; affects whether designs are filtered on AF3 disagreement).

Pipeline_notes13 records: "AF3 demoted from gating to flag-only across `merge_orthogonal_metrics.py`". The test-tree copy at `tests/orthogonal_metrics/merge_orthogonal_metrics.py` indeed has the demotion (`gating_flags = [f for f in flags if not f.startswith("af3_nomsa")]`, line ~143). The production `bin/merge_orthogonal_metrics.py` does NOT have the demotion — it still computes `passes = 1 if (af3_pass and not flags) else 0` (line 133) where `af3_pass` is `True` only when the AF3 metric is present and below threshold. The two files diverge on this single semantic point. This is the canonical "patch applied to test, never mirrored to prod" failure of the test-script-first iteration workflow.

### A5. `project_map.md` claims `pipeline_correct_sequences.py` is "the only file in the project whose ownership crosses the two codebases"

That claim is consistent with the file's history (lives in receptor-resurfacing tree, patched in negsteer sessions). With the codebases merged, the framing is moot, but the file remains the most cross-cutting Python module — patched again in pipeline_notes11 (the trailing-denovo native-residues bug fix, lines 382-411). Worth knowing: any change here ripples to both the upstream MPNN flow AND any negsteer integration test.

---

## B. Files referenced but missing / files present but unused

### B1. `bin/sequence_registry.py` is dead code

**Severity**: low (no impact, but ~175 LOC of cognitive overhead).

Zero importers in `bin/`, zero invocations from `.nf` or `.sh`. Already flagged as "confirmed dead code" in `notes/notes/notes9.md`. Safe to delete; would be a Phase 3.1 candidate per the remediation plan.

### B2. Empty file `main` at the repo root

**Severity**: low (trivial).

`/Users/jowillia/Documents/GitHub/receptor-resurfacing-pipeline/main` is a 0-byte file with mtime 2026-04-30 13:00 (today). Almost certainly an accidental shell redirect (`> main` or similar) — the `nf` was probably the intended target. Safe to delete; pre-baseline scrub. Worth checking with `git log` if it's tracked.

### B3. `submit_boltz2_iterate_steering.sh` referenced in `bin/boltz2_iterate_steering.py:1998` (`cmd_kickoff`) but absent from disk

**Severity**: high (the pipeline's `cmd_kickoff` subcommand may try to invoke a submission script that no longer exists).

The user pre-flagged this. The `cmd_kickoff` function in `boltz2_iterate_steering.py` was responsible for chaining cycle 1+ submitters in the old standalone pipeline. With the Nextflow integration, `negative_steering_run_one.sh` handles single-cycle, but multi-cycle stages 1+ are still represented in the orchestrator code via `cmd_iterate_plan` / `cmd_kickoff` / `cmd_kickoff_distances` / `cmd_kickoff_prefilter` / `cmd_kickoff_finalize`. Whether any of these are reachable in the current Nextflow flow is worth a separate audit — main.nf currently only invokes the single-cycle path (`NEGSTEER_RUN_ONE`). If multi-cycle is intentionally deferred (as `project_map.md §3` and `pipeline_notes13` suggest with the "Step 5 (multi-cycle) becomes the next foreground item" remark), then a substantial portion of `boltz2_iterate_steering.py` (the `cmd_kickoff*` family, ~600 LOC) is currently dormant code paths.

### B4. Several test-only Python files duplicate production files

Listed in §C below.

---

## C. Surprising duplication

### C1. `bin/merge_orthogonal_metrics.py` ↔ `tests/orthogonal_metrics/merge_orthogonal_metrics.py` (divergent)

Already covered in A4. Two copies of the same script that are now semantically different.

### C2. The six test/production plot-script pairs (verbatim except for header + fallback)

Documented as the deliberate "test iterates, prod mirrors" workflow. But it does mean every plot fix touches two files:

| Production | Test (source of truth) | Diff size |
|---|---|---|
| `bin/rfdiffusion_plots.py` | `tests/rfdiffusion/test_rfdiffusion_plots.py` | ~1073 lines |
| `bin/mpnn_plots.py` | `tests/proteinmpnn/test_mpnn_plots.py` | ~657 lines |
| `bin/rosetta_filter_plots.py` | `tests/rosetta_filtering/test_rosetta_filtering_plots.py` | ~285 lines |
| `bin/negsteer_plots.py` | `tests/negative_steering/test_negsteer_plots.py` | ~30 lines |
| `bin/negsteer_within_sequence_plots.py` | `tests/negative_steering/test_negsteer_within_sequence_plots.py` | ~15 lines |
| `bin/orthogonal_metrics_plots.py` | `tests/orthogonal_metrics/test_orthogonal_metrics_plots.py` | ~25 lines |

The first three pairs have substantial diffs because the test scripts predate the verbatim-mirror discipline (they contain alternative implementations of the iterating-target plots, plus an `_import_production_plots` runtime importer for the "good" plots). The last three pairs (negsteer + orthogonal) are the cleaner application of the discipline — small diffs, clear handoff. A future refactor that promotes the plot bodies into a shared library and leaves both copies as thin entry points is a clean target.

### C3. `derive_design_region.py` and `derive_true_interface.py` share ~80% structure

Both files have the same four functions: `normalise_design_id`, `find_design_entry`, `derive_positional_indices`, `cross_check_*`, `format_output`. Only `derive_positional_indices` differs in body (one emits 1-based design-region indices, the other emits 0-based contact residues). Consolidatable into a single `derive_indices.py --mode design-region|true-interface`.

### C4. Single-chain sequence extraction implemented four times

- `bin/boltz2_negative_steering.py:get_chain_sequence` (canonical, 4 importers)
- `bin/build_control_sequences.py:_extract_chain_sequence` (defined locally; the file ALSO imports get_chain_sequence — has both)
- `bin/extract_survivor_manifest.py:_extract_chain_seq`
- `bin/pipeline_correct_sequences.py:get_pdb_sequence`

### C5. Cα reading implemented several times

- `bin/boltz2_negative_steering.py:read_ca_atoms` (CAEntry dataclass list)
- `bin/rfdiffusion_filter.py:read_ca_atoms`
- `bin/parse_af3_output.py:_chain_ca_coords` (positional-pair list)
- `bin/reversion.py:_read_ca_chain` (length-only guard helper)

### C6. Contig parsing implemented in four places

- `bin/contig_utils.py:parse_block_segments`
- `bin/derive_input_design_region.py:_parse_contigs`
- `bin/haddock3_prepare.py:parse_contig_segments`
- `bin/pipeline_correct_sequences.py:parse_contig_segments`

### C7. THREE_TO_ONE amino-acid dict defined in 3 files (plus a variant)

- `bin/build_contigs.py:THREE_TO_ONE`
- `bin/boltz2_negative_steering.py:THREE_TO_ONE` AND `_THREE_TO_ONE_RMSD`
- `bin/extract_hotspots.py:THREE_TO_ONE`

### C8. Weighted-Jaccard primitives defined twice

`_WJ_MU = 4.0`, `_WJ_TWO_SIGMA_SQ = 2.25`, `_gaussian_contact_weight`, `_collect_chain_cb_positions`, `_build_weighted_pair_map` exist in BOTH `bin/compute_metrics.py` AND `bin/compute_interface_metrics.py` with identical bodies. Comments in `compute_interface_metrics.py:421-422` indicate awareness; the consolidation never happened.

### C9. `make_empty_plot`/`save_fallback_plots` reimplemented per plot script

Found in `bin/haddock3_plots.py`, `bin/rfdiffusion_plots.py`, `bin/mpnn_plots.py`, `bin/rosetta_filter_plots.py` (and as `_make_empty_plot` in the negsteer/orthogonal plots). Each is ~10 lines, near-identical bodies. Easy consolidation candidate.

### C10. `classify_seed(row)` is reimplemented in two negsteer plot scripts

`bin/negsteer_plots.py:737` and `bin/negsteer_within_sequence_plots.py:223` both define `classify_seed(row) -> Tuple[str, str]`. Different bodies (one uses `_is_steered`, the other doesn't). Conceptually the same operation; semantically distinct. Each was independently iterated.

### C11. Tier constants defined identically in 3+ places

`COLOUR_TIER`, `_TIER_SORT_ORDER`, `COMPOSITE_RA_EFF_WEIGHT` defined in `bin/negsteer_plots.py`, `bin/negsteer_within_sequence_plots.py`, `bin/orthogonal_metrics_plots.py` (plus all three test copies). One source of truth would be cleaner.

---

## D. Architectural drift / structural smells

### D1. `bin/boltz2_negative_steering.py` is a hub library disguised as a CLI

Two responsibilities in one 3416-LOC file: (a) the `cmd_plan` / `cmd_predict_one` / `cmd_collect` CLI for single-cycle steering, and (b) a dependency-graph hub from which 4 other scripts import named symbols. There's no `__all__` declaration, so the public API is implicit. A natural deep-modules refactor: extract the shared primitives (`get_chain_sequence`, `run_boltz`, `find_contact_residues_heavy`, `jaccard`, `binding_rmsds`, `write_boltz_yaml`, `extract_sequences`, the Boltz I/O helpers) into a pure library module (e.g. `bin/_boltz_lib.py` or a `bin/boltz/` package), leaving the CLI thin.

### D2. `bin/compute_metrics.py` is a hub library + subprocess-target

Imported by 3 files (`compute_interface_metrics.py`, `derive_input_design_region.py`, `boltz2_negative_steering.py`) AND invoked via subprocess by 2 (`boltz2_iterate_steering.py`, `reversion.py`). Move/rename requires touching both import sites and CLI default paths in two callers' argparse blocks. Same recommendation as D1: clean library/CLI split.

### D3. The single biggest file is ~6000 LOC

`bin/boltz2_iterate_steering.py` at 5854 LOC has 14 subcommands. The remediation plan's Phase 4 ("Address giant files thoughtfully") explicitly notes some 6000-line files are legitimately one module if cohesive. This file is split into clear regions (pathway helpers / mutation handling / maximin selection / 14 cmd_* dispatch / aggregation helpers / final-metrics utilities) that suggest natural module boundaries. The 14 subcommands also represent a particular shape — the SLURM-array-stage pattern — that may be re-evaluable once the Nextflow integration is done (some `cmd_kickoff*` family members may be dormant; see B3).

### D4. Inline narrative comments in main.nf load-bearing

Roughly half of `main.nf`'s 1061 lines is narrative comments documenting invariants (the post-RFDiffusion chain hardcoding, the params.outdir path leak, the parse_script cache pattern). Any structural restructure of main.nf will need to preserve these. They're effectively the workflow's design documentation; moving them to a separate doc would lose the at-the-call-site context.

### D5. `bin/negative_steering_run_one.sh` is ~450 LOC of orchestration

It chains ~13 subcommand invocations of `boltz2_negative_steering.py`/`boltz2_iterate_steering.py`/`extract_passing.py` through SLURM-array primitives. This is the kind of glue code that's expensive to test end-to-end (requires SLURM + GPU + Boltz container) and easy to break when any of the 14 subcommands changes its argparse. A future refactor could either (a) absorb the chain into Nextflow processes (one process per stage, getting Nextflow's caching + parallelism) or (b) move the orchestration into a Python wrapper (testable in isolation). The current pattern is justified by the 64+ MPNN-sequence scaling concern, but the tradeoff cost is real.

### D6. Three "derive" scripts in negsteer (`derive_design_region`, `derive_true_interface`, `derive_input_design_region`)

The first two share ~80% structure (see C3). The third is conceptually doing the same job for the controls flow but takes a fundamentally different input (input complex PDB instead of `rfdiffusion_metrics.json`) — so it's not a verbatim duplicate, but the three together represent one concept ("compute the design region + true interface for this design") with two inputs (post-RFDiffusion or pre-RFDiffusion). A unified API with mode flags is a natural Phase 4 deep-modules opportunity.

---

## E. Observations worth flagging

### E1. The Nextflow caching gotcha is partially mitigated, not fully fixed

Pipeline_notes14 documents the cache-key gotcha (Nextflow hashes the interpolated command string, not the contents of `${projectDir}/bin/X.py`). The mitigation pattern — `path script_input` in the process input list — has been applied to `AF3_PARSE_OUTPUT`, `ORTHOG_PLOTS`, and `NEGSTEER_ROSETTA_METRICS` (with the FastRelax XML). Notes 14 explicitly tracks "Task 61 — audit every `${projectDir}/bin/<script>.py` invocation in the .nf modules and convert to `path script_input` pattern" as outstanding. Until done, every script edit needs a corresponding .nf body touch to invalidate the cache. Worth flagging as a remediation prerequisite — refactoring without this fix means edits silently no-op when re-running with `-resume`.

### E2. `params.outdir` workaround in `tests/full_test_run/params.yml` is fragile

Pipeline_notes14 Bug 2 documents that relative `outdir: "./results"` paths leak through to container CWDs and break path-resolution. The workaround was making `outdir` absolute in the test params file; the real fix (`.toAbsolutePath()` in main.nf or argparse validation) was deferred to "Task 60." Worth knowing during refactoring — anyone editing the example params.yml could re-trigger this trap.

### E3. The bare `except: pass` audit is incomplete

Pipeline_notes13 specifically: "the broader lesson is to audit the codebase for similar patterns. **Added to Task 46 audit list.**" One specific instance was instrumented in `reversion.py:866`. A grep across `bin/` for bare `except: pass` patterns is a free-lunch Phase 1 task.

### E4. Several `cmd_*` functions in `boltz2_iterate_steering.py` are 200-500 LOC each

Examples: `cmd_iterate_plan` (~322 LOC), `cmd_iterate_collect_finalize` (~356 LOC), `cmd_aggregate` (~405 LOC), `cmd_aggregate_per_sequence` (~457 LOC), `cmd_compute_final_metrics` (~609 LOC). These are the longest functions in the codebase by far. Each represents one phase of the negsteer chain; many have grown organically as bugs were fixed. They're prime candidates for a complexity-report pass with `radon` (per remediation plan §1.4).

### E5. `_REVERTED_CONFIDENCE_FIELDS` (in `boltz2_iterate_steering.py`) is a 30-element string tuple

The propagation list for reverted-prediction columns. Pipeline_notes12 describes 6 separate bugs that were caused by this list being out of sync with `reversion.py`'s per_seed_records dict. The two have to be edited in lockstep — changing one without the other silently drops columns. Worth a single source-of-truth refactor: either the list lives in `reversion.py` and is imported, or both sides derive from a shared schema constant.

### E6. The `--n-cycles 1` silent-skip-reversion bug from notes6.md may not have been fixed

Notes 6.md (`notes/notes/notes6.md`) lists "post-implementation bugs" including "`--n-cycles 1` silently skipping reversion." I did not verify whether this is still present. Worth confirming during the Phase 2 characterization-test phase. The prompt was to NOT modify code; this is just a flag for the next pass.

### E7. The `tests/full_test_run/` directory is ambiguously scoped

It contains a `params_test.yml` (165 LOC) and a single-row `rosetta_summary.csv`. There's no `.nf` workflow here — the directory exists to hold the params file for full end-to-end runs of `main.nf`. That's fine, but the directory name suggests test infrastructure that lives elsewhere; future readers may go looking for a workflow file that isn't there.

### E8. The "test iterates, production mirrors" pattern for plots assumes discipline

Each pair of plot scripts (e.g. `bin/negsteer_plots.py` and `tests/negative_steering/test_negsteer_plots.py`) carries a header comment instructing future editors to edit the test version first. The discipline has held in some places (negsteer + orthogonal pairs are nearly identical) but not others (the rfdiff/mpnn/rosetta pairs have substantial drift). Without an automated check (a pre-commit hook diffing the two files modulo header + fallback), the pattern will continue to drift.

### E9. The pipeline has zero unit tests and that's by design

Per the remediation plan §2: "Don't try to write comprehensive unit tests yet — that's a Phase 5 activity." The existing `tests/` directory is module-level integration tests (each runs a sub-workflow on real test inputs), not pytest-style unit tests. This is consistent with the plan; recording explicitly so it's not mis-read as a defect.

### E10. The `tests/full_test_run/` overnight cohort run completed today

Pipeline_notes14 (today's session) describes a successful full-pipeline run in `tests/full_test_run/` producing 32 designs × 4 mpnn = 32 sequences + 2 controls = 34 total, with 4 tier-A / 3 tier-B / 27 tier-none. The four tier-A candidates (d19_s1, d0_s2, d13_s1, d13_s2) are described as the supervisor demo highlights. **This means the codebase as-is produces the demo result the user is showing, so the Phase 2 characterization-test step has a recent golden-master to anchor against** (the produced cohort summary plot, the cross-sequence CSV, the per-sequence aggregated CSVs). Worth preserving these outputs before any Phase 3 changes touch the code.

### E11. The XML at `bin/fastrelax_for_ia.xml` is a Rosetta protocol, not a Python/Nextflow file

It's a 59-line RosettaScripts protocol pinned for the Bennett-2023 ΔΔG calibration. It's correctly placed in `bin/` (where `run_rosetta_metrics.py` reads it), and it's wired through main.nf via `Channel.value(file("${projectDir}/bin/fastrelax_for_ia.xml"))` — the same path-input pattern used for cache-busting Python scripts. Worth noting as the only non-Python/non-shell asset in `bin/`; if a future structural reorg moves things into subdirectories, this XML's path-input flow needs to follow.

### E12. The `compare_versions.py` cross-version validator is missing

Mentioned in `project_map.md §2.3` as central to the production criterion (`intact AND ra_eff < 3.0 AND true_jaccard >= 0.7`) used to score designs into production / borderline / failed tiers. Not in current tree. The functionality may have been absorbed into `cross_sequence_summary.py:_tier_for_row` and `extract_passing.py`'s tier classification — it'd be worth confirming during the characterization-test phase that the current production tier rules match what `compare_versions.py` would have produced.

---

## F. Things I expected but didn't find

### F1. No `requirements.txt` / `pyproject.toml` / `environment.yml`

Python dependencies are implicit (everything runs in Singularity containers per `nextflow.config`). For local-dev or pre-commit static-analysis (`vulture`, `ruff`, `radon` per remediation plan §1.3), the dev environment will have to be set up by hand. Worth adding a `dev-requirements.txt` or similar for Phase 1 tooling.

### F2. No `.pre-commit-config.yaml` or CI configuration

No git hooks, no GitHub Actions / GitLab CI / Jenkins file. All linting / testing is manual. A `.pre-commit-config.yaml` running `ruff` + a simple consistency check (e.g. diff `bin/X_plots.py` against `tests/<m>/test_X_plots.py` modulo expected diffs) would be a high-value Phase 6 (going-forward discipline) addition.

### F3. No `.gitattributes` for the large CSVs / data files

Several test data files live in the repo (`tests/rfdiffusion/data/af3_pikp1_native_avrpikf_complex.pdb`, etc.). If any grow large (overnight runs sometimes produce multi-MB summary CSVs), an `.gitattributes` with LFS pointers would prevent repo bloat. Not urgent but worth a Phase 6 consideration.

### F4. No top-level README.md

Repo root has `params_example.yml` and `run_pipeline.slurm.sh` but no README. The codebase's design concept is currently spread across 14 `pipeline_notes*.md` + a `project_map.md` + this inventory. A short README pointing future-self (and any collaborator) at the right starting documents would close the loop.
