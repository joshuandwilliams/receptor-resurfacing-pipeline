# 04 — Functional Categorization

The codebase's design concept, made explicit. The categories below emerged bottom-up from grouping `.nf` modules with the `bin/` scripts each invokes. That is the categorization the code itself implies, cross-checked against the domain narrative in `notes/project_map.md` and `notes/pipeline_notes*.md`.

The pipeline is one Nextflow DSL2 workflow (`main.nf`). It takes either a pre-docked complex PDB or a receptor+effector PDB pair. It generates candidate receptor mutants that bind the effector, then evaluates each candidate with structural prediction + multiple orthogonal metrics. There is no longer a separate "negative-steering pipeline". The project_map's intended end-state, negsteer integrated as the Boltz2 validation stage, has been realised. What was once two codebases is now one.

The natural seven-category grouping below corresponds roughly to the stages of the workflow.

---

## Category 1 — Pipeline scaffolding

**Files**: `main.nf`, `nextflow.config`, `params_example.yml`, `run_pipeline.slurm.sh`.

The top-level orchestrator. `main.nf` is a single workflow block. It includes 31 processes from 13 module files and wires them in topological order. Entry is via optional Branch-A (HADDOCK docking from two PDBs) or Branch-B (pre-docked complex). `nextflow.config` is the canonical home for container paths (`rfdiff_container`, `rosetta_container`, `boltz2_container`, `colabfold_container`, etc.) and per-process SLURM resources. Pipeline_notes1 carried out that consolidation. `params_example.yml` is the user-copyable parameter starter. `run_pipeline.slurm.sh` is the SLURM submission wrapper.

**Internal coherence**: clean. The only smell is in `main.nf` itself. At 1061 LOC it is the longest top-level Nextflow file in the project, and it carries substantial inline narrative comments. The file is half code, half explanation. That is not a defect on its own. The comments document the post-RFDiffusion chain-convention invariant (A=receptor / B=effector hardcoded by `write_split_pdb`) and the params.outdir relative-path leak workaround. But any restructure of main.nf must preserve that documentation, which is non-trivial.

---

## Category 2 — Preprocessing

**Files**: `modules/preprocessing.nf`, `bin/rfdiffusion_contigs.py`, `bin/contig_utils.py`.

This category covers three steps. Sequence extraction from input PDBs (`EXTRACT_SEQUENCES`, inline Python in the .nf). Contigs resolution (`RESOLVE_CONTIGS` → `bin/rfdiffusion_contigs.py`, a 42-LOC CLI wrapper around `contig_utils.resolve_contigs`). A 2-byte dummy-mapping JSON writer (`WRITE_DUMMY_MAPPING`, inline `printf`).

`bin/contig_utils.py` holds the shared contig-parsing primitives: `parse_block_segments`, `remap_segments_to_pdb`, `resolve_contigs`, `get_expected_chain_lengths`, `parse_design_region`. It is imported by `rfdiffusion_contigs.py` and `rfdiffusion_filter.py`.

**Internal coherence**: mostly clean. But contig parsing is one of the canonical examples of duplication in the codebase. There are at least four independent contig parsers across `bin/`: `contig_utils.parse_block_segments`, `derive_input_design_region._parse_contigs`, `haddock3_prepare.parse_contig_segments`, `pipeline_correct_sequences.parse_contig_segments`. Their responsibilities are similar but not identical. Some resolve denovo lengths, others extract specific subsets. Consolidation candidate (flagged in `05_findings.md`).

---

## Category 3 — Branch A: HADDOCK docking

**Files**: `modules/haddock.nf`, `bin/haddock3_prepare.py`, `bin/collect_haddock3_dock.py`, `bin/haddock3_plots.py`, `bin/extract_hotspots.py`, `bin/build_contigs.py`, `bin/haddock_utils.py`.

Branch-A entry point. When the user supplies receptor + effector as separate PDBs, HADDOCK3 docks them. The best cluster's best model becomes the input PDB to RFDiffusion. The five processes are `HADDOCK3_PREPARE`, `HADDOCK3_DOCK`, `HADDOCK3_PLOTS`, `EXTRACT_HOTSPOTS` and `BUILD_CONTIGS`. They form a linear chain: AIR restraints from contigs → docking → cluster picking → diagnostic plots → hotspot extraction → contigs reconstruction. `haddock_utils.py` is the module's library (CAPRI/clustfcc parsers, heavy-atom extraction, cluster scoring), imported by 3 scripts.

**Internal coherence**: this category was the deepest cleanup target in pipeline_notes1. Twelve fixes were applied, including a sequence-identity-fallback chain disambiguation (the H10 fix). Notes flagged HADDOCK plot work as out of scope for the current cleanup pass ("HADDOCK work is going to come later on"). Contigs handling here (`build_contigs.py`'s `parse_contig_receptor_refs`, `compute_offset`, `update_contigs`) is yet another contig-parsing implementation that bypasses `contig_utils`.

---

## Category 4 — RFDiffusion + Rosetta filtering + ProteinMPNN

**Files**:
- `modules/rfdiffusion.nf`, `bin/rfdiffusion_filter.py`, `bin/rfdiffusion_plots.py`
- `modules/rosetta_filtering.nf`, `bin/rosetta_filter_collect.py`, `bin/rosetta_filter_plots.py`
- `modules/proteinmpnn.nf`, `bin/pipeline_correct_sequences.py`, `bin/mpnn_sequence_qc.py`, `bin/mpnn_design_region_score.py`, `bin/mpnn_cluster_sequences.py`, `bin/mpnn_select_top.py`, `bin/mpnn_plots.py`

The "design-generation" middle of the pipeline. Take a complex PDB → RFDiffusion generates N backbones → Rosetta InterfaceAnalyzer filters by Sc → ProteinMPNN designs sequences for the surviving backbones → QC, region-scoring, clustering, top-N selection. After this stage the pipeline holds a cohort of designed receptor sequences. These are what the negsteer stage validates.

One internal seam matters. `rfdiffusion_filter.py:write_split_pdb` HARDCODES output chain labels A=receptor, B=effector, regardless of input chain naming. Every downstream consumer treats those labels as fixed (documented in main.nf). This is the canonical "untracked invariant", the kind of thing the remediation plan calls out.

`pipeline_correct_sequences.py` (770 LOC) is the cross-pipeline file. It lives in the receptor-resurfacing tree but was patched in negsteer sessions. Its `split_sequence` chain-order bug was producing systematically-wrong receptor sequences in negsteer's RFDiffusion-design test runs. Its anchor-based `align_to_native_by_anchors` is the single most subtle piece of logic in this category.

**Internal coherence**: each sub-module has a clear linear chain. The plot scripts are `rfdiffusion_plots.py` (1248 LOC), `mpnn_plots.py` (550 LOC) and `rosetta_filter_plots.py` (206 LOC). They are the most-iterated files in the project, and pipeline_notes11 reworked all three. All three carry the same pattern: `parse_args` / `main` / module-level `ALL_PLOT_FILES` / `make_empty_plot` / `save_fallback_plots` / one function per plot. Each has an iteration-script counterpart in `tests/<module>/`. The `make_empty_plot` / `save_fallback_plots` idiom is reimplemented separately in each plot script (4 copies). Mild observation, not a critical defect.

---

## Category 5 — Negative-steering core

**Files**:
- `modules/negative_steering.nf`, `modules/negsteer_controls.nf`
- `bin/negative_steering_run_one.sh` (449 LOC orchestrator), `bin/fastrelax_for_ia.xml`
- `bin/boltz2_negative_steering.py` (3416 LOC), `bin/boltz2_iterate_steering.py` (5854 LOC), `bin/reversion.py` (1152 LOC)
- `bin/extract_passing.py`, `bin/cross_sequence_summary.py`
- `bin/derive_design_region.py`, `bin/derive_true_interface.py`, `bin/derive_input_design_region.py`
- `bin/build_control_sequences.py`
- `bin/sequence_registry.py` (dead code)

The heart of the validation stage and by far the largest category. For each MPNN-designed sequence, `NEGSTEER_RUN_ONE` runs the full single-cycle steering experiment. The steering phase: cold-start Boltz prediction → identify wrong-interface residues → emit N steered design YAMLs → predict each → rank → kickoff distances/prefilter/finalize. The reversion phase: identify mutations contacting effector inside the protected set (`design_region ∪ true_interface`) → plan reversions of those contaminated mutations → re-predict reverted designs. Then harvest reversion outputs → classify each design's reversion verdict → aggregate per-sequence → extract passing summary. The whole chain is collapsed into one SLURM GPU job per MPNN sequence. That avoids the submission explosion the project_map's older SLURM-array chain would create at 64+ MPNN sequences. The two negative controls (`scrambled`, `polyA`) go through the same chain.

`bin/cross_sequence_summary.py` aggregates every MPNN sequence's `passing_summary.csv` into one cross-sequence CSV with tier ranking and composite-score sorting. Tier-A/B/C/none classification + the `true_jaccard − 0.05·ra_eff` composite drive the cohort plots and the hand-off to the orthogonal-metrics stage.

The three "derive" scripts are `derive_design_region.py`, `derive_true_interface.py` and `derive_input_design_region.py`. Each is a CLI computing one of three closely-related index sets. The source is either RFDiffusion metrics or the input complex, depending on whether the input is a real MPNN design or a control.

**Internal coherence**: this is the category most affected by accumulated complexity. Several specific observations:

1. `boltz2_iterate_steering.py` (5854 LOC, 14 subcommands) and `boltz2_negative_steering.py` (3416 LOC, 3 subcommands) together contain ~9300 LOC of negsteer logic. The split is by phase, single-cycle versus multi-cycle orchestration. But `boltz2_negative_steering.py` is also a hub library that 4 other scripts import named symbols from, which makes its public API implicit. Promoting those exported helpers into a dedicated library module is a natural deep-modules candidate. The helpers are `get_chain_sequence`, `run_boltz`, `find_contact_residues_heavy`, `jaccard`, `binding_rmsds` and `write_boltz_yaml`.

2. `derive_design_region.py` and `derive_true_interface.py` share an identical four-function quartet: `normalise_design_id`, `find_design_entry`, `derive_positional_indices`, `format_output`. Only the `derive_positional_indices` body differs. They are a candidate for a single `derive_indices.py` taking a `--mode design-region|true-interface` flag.

3. `bin/sequence_registry.py` is dead code (zero importers, zero invocations from `.nf` or `.sh`). Already flagged as such in pipeline_notes9.

4. `bin/negative_steering_run_one.sh` (449 LOC) is the negsteer chain in shell form. Each phase of the chain runs one of `boltz2_negative_steering.py`/`boltz2_iterate_steering.py`/`extract_passing.py`. The boundary between this shell script and the Nextflow process body in `modules/negative_steering.nf` is mostly arbitrary. A future refactor could move some or all of the shell logic into inline Nextflow (parallelisable per-stage) or a Python orchestrator (more testable). The current "one big GPU job" was a deliberate choice to avoid SLURM submission explosion at scale (notes commentary in `negative_steering.nf`).

5. The reversion-pass added in pipeline_notes6 (notes/notes/notes6.md) introduced `reversion.py` and propagation of `_REVERTED_CONFIDENCE_FIELDS`. Pipeline_notes12 surfaced 6 latent bugs in this propagation chain (jaccard gating, per-seed jaccard, truth detection, reverted interface metrics, aggregator blanking, reverted extended metrics). Most were fixed. One was in flight, the `expected_rec_seq` hypothesis, and was proven + fixed in pipeline_notes13.

---

## Category 6 — Orthogonal-metrics validation

**Files**:
- `modules/negsteer_manifest.nf`, `bin/extract_survivor_manifest.py`
- `modules/negsteer_interface_metrics.nf`, `bin/compute_interface_metrics.py`
- `modules/negsteer_af3_nomsa.nf`, `bin/parse_af3_output.py`
- `modules/negsteer_biophysical_metrics.nf`, `bin/run_biophysical_metrics.py`
- `modules/negsteer_rosetta_metrics.nf`, `bin/run_rosetta_metrics.py`, `bin/fastrelax_for_ia.xml` (shared with Category 5)
- `modules/negsteer_orthogonal_metrics.nf`, `bin/merge_orthogonal_metrics.py`
- `bin/compute_metrics.py` (1657 LOC, library role here)

After negsteer finishes, the survivors (designs that passed tier-A/B/C) are validated using "orthogonal" metrics. These are methods independent of Boltz, so disagreement between Boltz and them flags potential prediction artefacts. The `EXTRACT_SURVIVOR_MANIFEST` step fans the survivors out to three parallel streams. Those are AF3-no-MSA prediction, biophysical metrics (BSA / interface_plddt / interface_hbonds), and Rosetta metrics (Sc / ΔΔG via FastRelax → InterfaceAnalyzer). `NEGSTEER_INTERFACE_METRICS` separately appends iRMSD/fnat/DockQ/15Å iPSAE/intact_core/weighted_jaccard to the cross-sequence CSV (single CPU process). Finally `NEGSTEER_ORTHOGONAL_METRICS` merges the three stream summaries and computes `passes_orthogonal_filters`.

`compute_metrics.py` plays a dual role. It runs as a subprocess under `boltz2_iterate_steering.py` (for confidence-metric extraction during finalize) and `reversion.py` (for reverted-prediction confidence). It is also imported as a library by `compute_interface_metrics.py` (`compute_ipsae`), `derive_input_design_region.py` (`find_contact_residues_heavy`), and `boltz2_negative_steering.py` (`compute_effector_interface_residues`). That dual mode makes it hard to classify as Category 5 or Category 6. It is load-bearing for both.

**Internal coherence**: the four streams (interface metrics, AF3, biophysical, Rosetta) have clean per-stream module boundaries. Each is one short `.nf` file plus one Python script, all keyed on `seq_name`. The merge logic in `merge_orthogonal_metrics.py` is also small (204 LOC) and well-scoped. Two issues stand out:

1. **`merge_orthogonal_metrics.py` has a divergent test-tree copy** (`tests/orthogonal_metrics/merge_orthogonal_metrics.py`). The test version has the AF3-flag-only demotion. The production version still gates on AF3. Pipeline_notes13 says the demotion was applied, but only the test copy was updated.

2. The weighted-Jaccard primitive is implemented twice, in `compute_metrics.py` and `compute_interface_metrics.py`. Both copies use the identical constants `_WJ_MU=4.0` and `_WJ_TWO_SIGMA_SQ=2.25`, and the identical helpers `_gaussian_contact_weight`, `_collect_chain_cb_positions`, `_build_weighted_pair_map`. Comments in `compute_interface_metrics.py:421-422` indicate awareness. Consolidation never happened.

---

## Category 7 — Plot scripts (cohort + per-seed + orthogonal)

**Files**:
- `bin/negsteer_plots.py` (1911 LOC), `bin/negsteer_within_sequence_plots.py` (987 LOC), `bin/orthogonal_metrics_plots.py` (1146 LOC)
- All paired `tests/<module>/test_<module>_plots.py` iteration scripts
- (Module-level plot scripts `rfdiffusion_plots.py`, `mpnn_plots.py`, `rosetta_filter_plots.py`, `haddock3_plots.py` are listed in their parent categories above.)

The three negsteer-stage plot scripts are the largest individual files in the project after the negsteer cores. They render the cohort-level summary plots that are the supervisor-demo-facing outputs, for example `negsteer_tier_landscape.png` and `negsteer_seed_outcomes_heatmap.png`. The totals are 7 cohort + 5 per-seed + 4 orthogonal = 16 PNGs. Each was iterated in `tests/<module>/` first, under the pipeline_notes11 "test-script-first iteration" workflow. It was then mirrored to `bin/` verbatim minus a header and a test-path fallback in `_resolve_csv_path`.

**Internal coherence**: the test/production-mirror pattern is documented in each script's header. It is disciplined enough that the diff between the two copies of each script is small (header + ~10 lines of fallback path). It is still literally duplicated source code, and every plot fix has to be made twice. A future refactor could promote the plot bodies into a shared library. The test/prod scripts would become thin entry points, with the only real diff being which `_resolve_csv_path` is used. That would eliminate the synchronisation burden. The cost is one more layer of indirection.

The constants are also duplicated identically across the test and production copies, among them `COLOUR_TIER`, `COMPOSITE_RA_EFF_WEIGHT`, `STAGE_MARKER` and `_VERDICTS_WITH_REVERSION`. They are the canonical example of "the same constant defined in three places" that the remediation plan calls out.

---

## Cross-cutting concerns

A few responsibilities don't fit cleanly into any single category. They are cross-cutting, and worth recording so the remediation knows where to look.

### Single-chain sequence extraction from a PDB

There are at least **four independent implementations** of "given a PDB and a chain ID, return the 1-letter amino-acid sequence":

- `bin/boltz2_negative_steering.py:get_chain_sequence` (the one most other scripts import)
- `bin/build_control_sequences.py:_extract_chain_sequence` (also imports the above — has both)
- `bin/extract_survivor_manifest.py:_extract_chain_seq`
- `bin/pipeline_correct_sequences.py:get_pdb_sequence`

### Cα reading

Multiple Cα atom readers exist. `bin/boltz2_negative_steering.py:read_ca_atoms` returns a CAEntry dataclass list. `bin/rfdiffusion_filter.py:read_ca_atoms` is its own implementation. `bin/parse_af3_output.py:_chain_ca_coords` returns a positional-pair list, post Bug-6 fix. `bin/reversion.py:_read_ca_chain` is a length-only guard helper added in notes13. Each is "almost" the same but reflects subtly different domain needs.

### Three-letter→one-letter AA dict

`THREE_TO_ONE` is defined in `bin/build_contigs.py`, `bin/boltz2_negative_steering.py` and `bin/extract_hotspots.py`. `boltz2_negative_steering.py` also holds `_THREE_TO_ONE_RMSD`, a separate variant. Trivially consolidatable.

### Plot fallback idiom

`make_empty_plot(message, path)` + `save_fallback_plots(message)` exists with near-identical bodies in `bin/haddock3_plots.py`, `bin/rfdiffusion_plots.py`, `bin/mpnn_plots.py`, `bin/rosetta_filter_plots.py`. Each plot script also defines `ALL_PLOT_FILES` as a module-level list to drive the fallback.

### Tier classification + composite score

`cross_sequence_summary.py` defines `_tier_for_row(row)`. The `COLOUR_TIER` / `_TIER_SORT_ORDER` constants are defined identically in `bin/negsteer_plots.py`, `bin/negsteer_within_sequence_plots.py` and `bin/orthogonal_metrics_plots.py`, plus their test copies. `_composite_score(row)` is defined in cross_sequence_summary and reimplemented as `_composite_from_row` in the plot scripts. One concept, several implementations.

### Contig parsing

Four implementations enumerated above: `contig_utils.parse_block_segments`, `derive_input_design_region._parse_contigs`, `haddock3_prepare.parse_contig_segments`, `pipeline_correct_sequences.parse_contig_segments`. Each is scoped slightly differently. But a unified parser with mode flags is feasible.

These cross-cutting duplications are the cleanest targets for the Phase 3 mechanical-cleanup step in the remediation plan. They don't require architectural change to consolidate.
