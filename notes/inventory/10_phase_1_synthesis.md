# 10 — Phase 1 Synthesis

A focused briefing distilled from the inventory documents (`01`–`09`) and the remediation plan. This is the document to re-read at the start of every later phase.

---

## 1. Summary of the codebase

This is a single Nextflow DSL2 pipeline (`main.nf`, 1061 LOC) for de novo binder design. Given a receptor + effector PDB pair (or a pre-docked complex), it generates novel receptor backbones via RFDiffusion, redesigns their sequences with ProteinMPNN, validates each candidate with Boltz negative-steering, and cross-checks the cohort with three orthogonal-method streams (AF3-no-MSA, biophysical, Rosetta InterfaceAnalyzer). The codebase that the `project_map.md` calls "two related but distinct codebases" is now one — the standalone negsteer pipeline has been folded into the Nextflow workflow as the `NEGSTEER_*` modules. ~42k LOC total (~35k Python in `bin/` + `tests/`, ~2.7k Nextflow, the rest YAML/shell/XML).

Organisation is conventional Nextflow: `main.nf` orchestrates 13 process modules under `modules/`, each module shells out to one or more Python CLIs under `bin/`. Seven natural functional categories emerge from the topology: pipeline scaffolding, preprocessing, HADDOCK Branch A, RFDiffusion + Rosetta filtering + ProteinMPNN, negsteer core, orthogonal-metrics validation, and plotting. Cross-cutting helpers live as ad-hoc imports between `bin/` files — the most load-bearing being `boltz2_negative_steering.py` (3416 LOC) and `compute_metrics.py` (1658 LOC), each both a top-level CLI *and* a hub library that 3–4 other scripts pull named symbols from.

State of health: the pipeline produced a successful supervisor-demo run on 2026-04-30 (4 tier-A / 3 tier-B / 27 tier-none across 32 designs × 4 MPNN sequences + 2 controls), so functional behaviour is currently correct. But the codebase has all the predictable scars of long AI-assisted iteration: two files exceed 3000 LOC and seven functions register cyclomatic complexity F (>40), four independent contig-parsers, four independent chain-sequence extractors, six pairs of test/production plot scripts (some near-verbatim, some drifted), one fully dead module (`sequence_registry.py`), ~600 LOC of dormant multi-cycle code, and one production/test semantic divergence (`merge_orthogonal_metrics.py`) that means the production code disagrees with what the notes say it does. Average maintainability index is A-grade overall; a small set of files (the negsteer cores, `negsteer_plots.py`, `compute_metrics.py`) drag to C.

---

## 2. Priority targets for remediation

Ranked by combined signal — radon complexity, file size, vulture findings, and the cross-cutting roles each plays. The top 7 are unambiguous; the next 3 are second-order.

| Rank | File | LOC | Why it's priority | Work needed |
|---:|---|---:|---|---|
| 1 | `bin/boltz2_iterate_steering.py` | 5854 | Largest file in the project; 7 cmd_* functions at complexity F (max 72); MI rank C; 14 subcommands cover both single-cycle pipeline (live) and multi-cycle orchestration (dormant, ~600 LOC of `cmd_kickoff*`). `_REVERTED_CONFIDENCE_FIELDS` 30-tuple silently couples to `reversion.py`. | Phase 4 decomposition: split by phase boundary (live cmd_* family vs dormant kickoff family vs aggregation helpers vs final-metrics utilities). Confirm dormancy of `cmd_kickoff*` first, then either delete or quarantine. |
| 2 | `bin/boltz2_negative_steering.py` | 3416 | Second-largest file; MI rank C; `cmd_plan` complexity F=109 (the highest single function in the project); also a hub library that 4 other scripts import named symbols from with no `__all__`. Two roles in one file. | Phase 4 deep-modules refactor: extract the shared primitives (`get_chain_sequence`, `run_boltz`, `find_contact_residues_heavy`, `jaccard`, `binding_rmsds`, `write_boltz_yaml`, the Boltz I/O helpers, the Cα/Kabsch geometry) into a pure library module; leave a thin CLI behind. |
| 3 | `bin/compute_metrics.py` | 1658 | Hub library + subprocess target (3 importers, 2 subprocess callers). MI rank C. Five functions at D complexity. Owns the second copy of weighted-Jaccard primitives (also defined in `compute_interface_metrics.py`). Owns one of two `find_contact_residues_heavy` implementations. | Phase 3 consolidation: pick canonical home for weighted-Jaccard + `find_contact_residues_heavy`, then library/CLI split (same pattern as #2). |
| 4 | `bin/negsteer_plots.py` (+ test twin) | 1912 + 1916 | Two files, near-identical bodies. MI rank C. 5 functions at E complexity. Constants (`COLOUR_TIER`, `COMPOSITE_RA_EFF_WEIGHT`, `_VERDICTS_WITH_REVERSION`) duplicated identically. Vulture flags `unreachable code after 'return'` at line 1260 and `unused 'cross_tier'` at 261 (both copies). | Phase 4 (or end-of-Phase-3): promote plot bodies into a shared library; leave both files as thin entry points differing only in `_resolve_csv_path`. Closes the synchronisation burden the test/prod-mirror pattern carries. |
| 5 | `bin/rfdiffusion_plots.py` (+ test twin) | 1248 + 1163 | Largest module-stage plot script. Test twin has substantial drift (~1073-line diff per finding C2) — the discipline of "test iterates, prod mirrors" broke down here. Five functions at C/D complexity. | Same library extraction as #4, but with a pre-step of reconciling the two copies (audit which differences are intentional). |
| 6 | `bin/cross_sequence_summary.py` | 904 | `aggregate()` at complexity F=58 — the cohort-aggregation function in one giant block. Holds the single source of truth for tier classification and composite score (`_tier_for_row`, `_composite_score`) yet those constants are re-defined in 3+ plot scripts. | Phase 3 consolidation: extract tier + composite into a shared module that the plots import; Phase 4 decomposition of `aggregate()` into smaller passes. |
| 7 | `bin/rfdiffusion_filter.py` | 926 | `main()` at F=47 and `calc_scaffold_and_region_metrics` at D=30. Owns the chain-A=receptor / chain-B=effector hardcoded invariant that every downstream stage depends on; uses `scaffold` for the *motif* (inverted Baker-lab terminology, glossary §F3). Has its own `read_ca_atoms` and `find_contacts` implementations. | Phase 3 mechanical: rename `scaffold_*` → `motif_*`; promote the chain-A/B invariant into a module docstring. Phase 4: split `main` into orchestration + per-design metric computation. |
| 8 | `bin/pipeline_correct_sequences.py` | 770 | `generate_fixed_positions` E=39, `align_to_native_by_anchors` E=36. Cross-pipeline file (lives in receptor-resurfacing tree, has been patched in negsteer sessions). Contains its own `parse_contig_segments` (4th independent copy) and `get_pdb_sequence` (4th chain-extractor copy). The `align_to_native_by_anchors` is the single most subtle piece of logic in upstream design. | Phase 3: replace local contig parser with `contig_utils`. Phase 5 (deferred): unit-test the anchor alignment in isolation before any restructuring. |
| 9 | `bin/extract_passing.py` | 678 | `_summarize_mutations_across_seeds` D=25, `_muts_aa_for_positions` C=14. Stage-selection rule (Bug-4 fix from notes14) is implicit and easy to break. | Phase 3 cleanup: extract stage-selection predicate into a named helper; document the per-seed-vs-aggregated semantics inline. |
| 10 | `bin/haddock3_plots.py` | 697 | `plot_interface_heatmap` F=47 — the highest individual plot-function complexity. HADDOCK work was explicitly de-scoped in pipeline_notes1 ("HADDOCK work is going to come later"), so this is lower priority than its complexity score implies. | Defer to a post-cleanup HADDOCK pass; don't touch in Phase 3/4 unless it blocks something. |

Special case (not ranked because it's a one-line action, not a refactor):

- `bin/sequence_registry.py` — 175 LOC, zero importers, zero `.nf`/`.sh` invocations, already documented as dead in pipeline_notes9 and confirmed in `03_dependency_graph.md`. Phase 3.1 deletion target.
- `main` (empty 0-byte file at repo root) — accidental shell redirect from 2026-04-30, already flagged. Delete pre-baseline.

---

## 3. Patterns observed across the codebase

The findings cluster into ~7 issue categories. Some warrant a single global pass; others are inherently file-by-file.

**Cross-cutting duplication (best fixed globally in one Phase-3 pass).** The same primitive operation is implemented in 3–4 places:
- single-chain sequence extraction (4 implementations: `get_chain_sequence`, `_extract_chain_sequence`, `_extract_chain_seq`, `get_pdb_sequence`)
- contig parsing (4 implementations: `contig_utils.parse_block_segments` + 3 local parsers in `derive_input_design_region`, `haddock3_prepare`, `pipeline_correct_sequences`)
- `THREE_TO_ONE` amino-acid dict (3 files plus a `_THREE_TO_ONE_RMSD` variant)
- weighted-Jaccard primitives (`_WJ_MU`, `_WJ_TWO_SIGMA_SQ`, `_gaussian_contact_weight`, `_collect_chain_cb_positions`, `_build_weighted_pair_map` — identical in `compute_metrics.py` and `compute_interface_metrics.py`; the awareness comment is at `compute_interface_metrics.py:421-422`)
- `find_contact_residues_heavy` (2 distinct bodies under the same name in `compute_metrics.py` and `boltz2_negative_steering.py`)
- `make_empty_plot` / `save_fallback_plots` plot-fallback idiom (4 copies)
- tier constants `COLOUR_TIER`, `_TIER_SORT_ORDER`, `COMPOSITE_RA_EFF_WEIGHT` (3+ copies)
- `classify_seed` (2 distinct bodies in `negsteer_plots.py` and `negsteer_within_sequence_plots.py`)

**Naming inconsistencies (Phase 3, glossary-driven).** Five terminology issues are flagged in glossary §10: `survivor` no longer matches its connotation post-Bug-3-fix; `scaffold` is used for the *motif* (inverse of Baker-lab usage); `backbone` is used in prose for Cα-only operations; `design_region` is 1-based while `true_interface` is 0-based (a code accident, not a domain distinction); `construct` is overloaded between the plasmid noun and `construct_reliance_flag` predicate. Each is a coordinated rename; doing them all in one pass is cleaner than file-by-file.

**Complexity hotspots (Phase 4, per-file).** Concentrated in five files: `boltz2_iterate_steering.py` (7 F-rated functions), `boltz2_negative_steering.py` (`cmd_plan` at F=109), `reversion.py` (`harvest_reversion_results` at F=74), `cross_sequence_summary.py` (`aggregate` at F=58), `compute_metrics.py` (5 D-rated functions). These are not amenable to global treatment — each requires understanding the local domain logic.

**Test/production drift (Phase 3 audit + Phase 4 library extraction).** Six plot-script pairs follow the documented "test iterates, prod mirrors" workflow with varying discipline. The negsteer/orthogonal pairs (~15–30 line diffs) are clean; the rfdiff/mpnn/rosetta-filter pairs have substantial drift. One non-plot pair (`merge_orthogonal_metrics.py`) is *semantically* divergent — the test version has an AF3-flag-only demotion that pipeline_notes13 says was applied, but the production version still gates on AF3. That single divergence affects whether designs are filtered on AF3 disagreement, so it counts as both a duplication issue and a latent semantic bug.

**Dead and dormant code.** Three buckets: confirmed-dead (`sequence_registry.py`, the empty `main` file), high-confidence vulture findings (8 items: 3 unused imports, 2 unused variables, 2 unreachable post-return, 1 in plot test twin), and dormant multi-cycle code (`cmd_kickoff*` family in `boltz2_iterate_steering.py`, ~600 LOC, never reached by current `main.nf` per finding B3 — but possibly intentionally dormant for the multi-cycle roadmap). Bucket one is safe to delete in Phase 3.1; bucket two as soon as ruff `--fix` runs; bucket three needs a decision before deletion.

**Configuration sprawl.** Less severe than expected. Constants (cluster sizes, weighted-Jaccard parameters, tier weightings, intact cutoffs) are typically defined adjacent to use, not centralised, but they're not scattered across config files — they're scattered across Python modules. The mitigation is consolidation into a small constants module rather than a YAML; the project already has `nextflow.config` and `params_example.yml` for runtime parameters.

**Inline narrative comments.** Ruff finds 4614 comment lines (13% of total LOC). Concentrated in `main.nf` (~50%), `boltz2_iterate_steering.py` (1111 lines), `boltz2_negative_steering.py` (640 lines), and `orthogonal_metrics_plots.py` (212 lines). Many are load-bearing — the chain A/B invariant, the params.outdir leak workaround, the parse_script cache pattern, the post-Bug-12 propagation rules — and need preservation, not deletion. Phase 3.2 ("strip over-verbose comments") needs to be careful here: delete change-history commentary, keep invariant-documenting commentary.

---

## 4. Refinements to the remediation plan

The plan was written before the inventory; with the inventory in hand, several phases sharpen.

**Phase 2 — capture golden masters from the 2026-04-30 demo run.** Finding E10 records that the supervisor-demo run completed today and produced 32 sequences + 2 controls with 4 tier-A / 3 tier-B / 27 tier-none. **The exact outputs from that run should become the Phase 2 golden masters** rather than re-running the pipeline. Specifically: the cross-sequence CSV, every per-sequence `aggregated_results.csv` and `passing_summary.csv`, the 16 cohort PNGs (7 negsteer cohort + 5 within-sequence + 4 orthogonal), the survivor manifest, the merged orthogonal metrics CSV, and the `rfdiffusion_metrics.json` files for each design. Comparison tolerance: numerical CSV columns to ~1e-6, plots by structural identity (file size + a sampled-pixel hash) since matplotlib output is non-deterministic across runs.

**Phase 3 — sequence the cross-cutting consolidations first.** The plan lists 3.5 ("consolidate duplicated functionality") last. The inventory suggests the order should be: 3.1 (delete dead — quick win), 3.5 partial (consolidate the most-imported primitives: chain-sequence extraction, contig parsing, THREE_TO_ONE, weighted-Jaccard), 3.4 (centralise constants — covers the tier constants), 3.3 (naming — the glossary §10 renames), 3.2 (strip comments — last, because the consolidations may obviate some comment blocks). The reason is that the duplicated primitives sit at the top of the import graph; consolidating them once means the rename pass touches each canonical name in one place, not four.

**Phase 3 — fix the `merge_orthogonal_metrics.py` divergence as a Phase 3 task, not Phase 4.** Finding A4: production gates on AF3, test demotes to flag-only. This is one ~10-line semantic edit that brings production in line with the test (which is the documented intended behaviour per pipeline_notes13). It should happen before any plotting refactor that consumes `passes_orthogonal_filters`.

**Phase 3 — confirm `cmd_kickoff*` dormancy before any boltz2_iterate_steering decomposition.** Finding B3 strongly suggests these ~600 LOC are unreachable from the current Nextflow flow. A static-trace from `main.nf` and from `negative_steering_run_one.sh` will confirm. If dormant, they go into a Phase 3.1-style deletion (or quarantine into a `bin/_dormant/` directory if multi-cycle is on the roadmap). This single decision shrinks the largest file by ~10% and simplifies the Phase 4 split.

**Phase 4 — natural starting point is the `boltz2_negative_steering.py` library/CLI split.** The plan's "design the target architecture" step calls for the grill-me skill; before grilling, the obvious deep-module candidate is already visible. Extract `get_chain_sequence`, `run_boltz`, `find_contact_residues_heavy`, `jaccard`, `binding_rmsds`, `write_boltz_yaml`, `extract_sequences*`, the Cα/Kabsch geometry helpers, the residue dataclasses, and `THREE_TO_ONE` into a `bin/boltz_lib.py` (or a `bin/boltz/` package). Expose them via a real `__all__`. Update the four importers (`boltz2_iterate_steering`, `reversion`, `build_control_sequences`, `derive_input_design_region`). The `cmd_plan` / `cmd_predict_one` / `cmd_collect` CLI stays in a thin `boltz2_negative_steering.py`. Same pattern then applies to `compute_metrics.py`.

**Phase 4 — second target is the plot-script library.** All six plot pairs share the same shape: `_try_float` / `_try_int` / `_make_empty_plot` / `_short_name` plus per-script plot functions that differ by what data they consume. A `bin/plot_lib.py` module with the shared primitives and the `make_empty_plot` / `save_fallback_plots` / `ALL_PLOT_FILES` idiom would shrink each plot script substantially and let test/prod scripts collapse to thin entry points.

**Phase 4 — `derive_*.py` consolidation (small but clean).** `derive_design_region.py` and `derive_true_interface.py` share 4 functions verbatim. Combine into one `derive_indices.py --mode design-region|true-interface`. `derive_input_design_region.py` is conceptually the same operation for the controls flow but with a different input; whether it folds into the same module is a Phase 4 design decision.

**Phase 6 — add the test/prod-mirror consistency check as a pre-commit hook.** Finding E8: the test/prod mirror discipline holds when there's social pressure but drifts otherwise. A pre-commit hook diffing each pair (modulo expected diffs in `_resolve_csv_path`) is a one-time setup that prevents regression. This is a Phase 6 task per the plan; recording it now so it's not forgotten.

**Cross-phase prerequisite — the Nextflow cache-busting audit (notes14 Task 61).** Finding E1: the `path script_input` pattern has been applied to `AF3_PARSE_OUTPUT`, `ORTHOG_PLOTS`, and `NEGSTEER_ROSETTA_METRICS`, but every other process still hashes the interpolated command string rather than the `bin/X.py` content. Until this is fixed across the board, refactoring a script and re-running with `-resume` will silently use the cached old result. **This needs to happen before Phase 3 begins** or the characterisation tests will compare unchanged outputs and pass falsely.

---

## 5. Open questions and risks

**Genuine uncertainties from the inventory:**

- **Is `cmd_kickoff*` truly dead, or scheduled for revival?** Glossary §2 says multi-cycle is "currently dormant" and "on the roadmap." If revived, deletion is wrong. Need a yes/no decision before Phase 3.1 deletion of dormant code.
- **Does the production `merge_orthogonal_metrics.py` AF3 gating actually filter any rows in practice?** The semantic divergence with the test copy is real, but if AF3 always passes (or always present) on the current cohort, no rows are visible at the cohort plot level. Worth checking the 2026-04-30 cohort's AF3 column distribution before Phase 3 fixes it — the fix is correct either way, but the urgency depends.
- **Is the `--n-cycles 1` silent-skip-reversion bug from notes6.md still latent?** Finding E6: not verified during inventory. Worth confirming during Phase 2 characterisation tests. If still latent, the demo-run golden masters may *bake in* a bug.
- **The `pipeline_correct_sequences.py:align_to_native_by_anchors` function** at E=36 complexity is the single most subtle piece of upstream logic. Behavioural correctness hinges on it. Phase 2 characterisation needs to capture intermediate state from this function (the corrected sequence per design), not just the eventual MPNN output.
- **The `notes/full_test_run/` cohort outputs** referenced in finding E10 — are they checked into git, or only on the HPC? If only on the HPC, Phase 2 needs an explicit transfer step before any code touching can begin.

**Behaviours the characterisation tests will need to capture but that we haven't fully scoped:**

- Per-seed reversion verdict classification (is `pose_holds` deterministic across seeds for a given input?).
- The "clean steered" path where reversion is correctly skipped (per glossary §5: same underlying state as `no_reversion`, different granularity).
- The negative-controls path: cold-start prediction → no contact residues → silent skip → 3 cold-start-only rows per control. Easy to break with the steering-set / steering-mode logic.
- The path-rewriting in `_rewrite_workdir_path_to_published` — work-dir → published-tree path conversion is depended on by every downstream container task and is invisible to the import graph.
- Effector-template CIF generation in `extract_effector_template_cif` (E=33 complexity, 33-LOC function in `boltz2_negative_steering.py`).

**Risks specific to this remediation:**

- The HPC airgap means each round-trip is expensive. Sticking to Phase 1 + Phase 3.1 + glossary renames before any architectural restructure conserves round-trips by batching mechanical work.
- The narrative comments in `main.nf` and the negsteer cores are load-bearing. A Phase 3.2 comment-strip pass that doesn't preserve them risks losing the only documentation of the chain-A/B invariant and the cache-busting pattern.
- The `boltz2_iterate_steering.py` ↔ `reversion.py` `_REVERTED_CONFIDENCE_FIELDS` coupling (finding E5) means the propagation list is silently editable on either side. Until it's a single source of truth, *any* refactor that touches reverted-prediction columns risks reintroducing the notes12 bug family.

---

## 6. Recommended Phase 2 starting point

The smallest viable characterisation test suite — concretely, what to capture and what to compare.

**Inputs to capture (from the 2026-04-30 demo run already on disk):**

1. `tests/full_test_run/params_test.yml` — the parameter file that produced the demo cohort. Pin the exact contents.
2. The HADDOCK-output PDB or pre-docked-complex PDB used as input.
3. The `rfdiffusion_metrics.json` file for each design.

**Outputs to capture as golden masters (from the published tree of the 2026-04-30 run):**

| Tier | Artefact | Comparison |
|---|---|---|
| 1 (strongest) | `cross_sequence_summary.csv` | Exact CSV diff (column-by-column, row-ordered). Numerical columns at 1e-6 tolerance. |
| 1 | Each per-sequence `aggregated_results.csv` (32 + 2 = 34 files) | Same. |
| 1 | Each per-sequence `passing_summary.csv` (34 files) | Same. |
| 1 | `merged_orthogonal_metrics.csv` (post-orthogonal merge) | Same. |
| 2 | `survivor_manifest.csv` | Exact diff modulo path rewriting (workdir hashes will change between runs). |
| 2 | `rfdiffusion_metrics.json` per design | Exact diff. |
| 2 | The 16 cohort PNGs (7 negsteer cohort + 5 within-sequence + 4 orthogonal) | Image-similarity comparison: file size within 5%, structural-similarity index ≥ 0.95. Plots are matplotlib non-deterministic; exact-bytes comparison will fail. |
| 3 (weakest) | Per-design Boltz prediction CIF/PDB sidecars | Skip in initial pass — prediction non-determinism makes these unstable. Capture only if a metric regression is detected and we need to bisect. |

**Smallest viable test suite (the comparison-script set):**

1. `tests/characterization/compare_cross_sequence_csv.py` — diff `cross_sequence_summary.csv` against the golden master. ~50 LOC.
2. `tests/characterization/compare_per_sequence_csvs.py` — same for the 34 `aggregated_results.csv` + 34 `passing_summary.csv` files. ~80 LOC.
3. `tests/characterization/compare_orthogonal_merge.py` — same for `merged_orthogonal_metrics.csv` and `survivor_manifest.csv`. ~40 LOC.
4. `tests/characterization/compare_plots.py` — image-similarity over the 16 PNGs. ~60 LOC. Uses Pillow + a structural-similarity helper.
5. `tests/characterization/run_all.sh` — driver that walks an output directory and runs all four comparators against `golden/`. Returns 0 if all pass.

**Concrete first action items (in order):**

1. **Audit the cache-busting Task 61** (finding E1). One pass through the `.nf` modules to convert every `${projectDir}/bin/X.py` invocation to the `path script_input` pattern. Without this, characterisation tests are unreliable. Local + Mac, no HPC needed.
2. **Copy the 2026-04-30 demo outputs from HPC into `tests/characterization/golden/`.** One round-trip.
3. **Write the four comparator scripts above on the Mac.** Test them locally against a copy of the golden output (so they all pass on identity).
4. **Sync the comparators to HPC; run them against a fresh re-run of the same params file.** This validates that the comparators produce 0 false positives on a known-good rerun (modulo plot non-determinism, which the SSIM tolerance should absorb).
5. **Once green: delete `bin/sequence_registry.py` and the empty `main` file as the first Phase 3.1 commit.** Re-run characterisation. If green, the safety net works and the rest of Phase 3 can begin.

The key insight is that we *already have* the golden master — the demo run completed successfully today. We don't need to design synthetic inputs or re-run from scratch; we just need to capture and compare. That makes Phase 2 dramatically cheaper than the plan implies.
