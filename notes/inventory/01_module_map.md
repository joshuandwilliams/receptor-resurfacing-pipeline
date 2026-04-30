# 01 — Module Map

A complete inventory of every code file in the project, grouped by directory. Line counts are wc -l. One-line summaries describe what the file actually does (not what notes claim it should do — see `05_findings.md` for code/notes discrepancies).

Skipped: `notes/` (used for context, not inventoried), data files (`.pdb`, `.csv`, data `.yml`), and the `.git/` tree.

## Top level

| Path | Type | LOC | Summary |
|---|---|---:|---|
| `main.nf` | Nextflow DSL2 | 1061 | Top-level workflow. Wires preprocessing → HADDOCK (optional Branch A) → RFDiffusion → Rosetta → ProteinMPNN → negative-steering → orthogonal metrics → plots. |
| `nextflow.config` | Nextflow config | 387 | Single source of truth for container paths, per-process SLURM resources, default params (overrides values in params.yml for tested defaults). |
| `params_example.yml` | YAML | 223 | Annotated example parameter file users copy and edit. |
| `run_pipeline.slurm.sh` | shell | 104 | SLURM wrapper that submits the full Nextflow pipeline. |
| `main` | empty | 0 | **Empty file at repo root, mtime 2026-04-30. Almost certainly an accidental `>` redirect or `touch` mistake — flagged in `05_findings.md`.** |

## bin/ — backing scripts called by Nextflow processes (and by each other)

### Negative-steering core (largest files in the codebase)

| Path | Type | LOC | Summary |
|---|---|---:|---|
| `bin/boltz2_iterate_steering.py` | Python | 5854 | Multi-cycle negative-steering orchestrator. 14 subcommands: `iterate-plan`, `compute-distances`, `iterate-collect`, `iterate-collect-prefilter`, `build-contaminated`, `plan-reversions`, `predict-reversion-one`, `harvest-reversions`, `kickoff-distances`, `kickoff-prefilter`, `kickoff-finalize`, `aggregate`, `aggregate-per-sequence`, `compute-final-metrics`. |
| `bin/boltz2_negative_steering.py` | Python | 3416 | Single-cycle negative-steering core. Subcommands: `plan` (initial Boltz prediction → identify wrong-interface → emit N steered design YAMLs), `predict-one` (one design per array task), `collect` (rank). Holds shared helpers (Boltz YAML/run wrappers, contact-residue + Jaccard logic, Kabsch, sequence I/O) used by ALL other negsteer Python modules via direct `from boltz2_negative_steering import …` calls. |
| `bin/reversion.py` | Python | 1152 | Pure-logic helpers used by the multi-cycle orchestrator: `build_reverted_sequence`, `write_reversion_plan`, `harvest_reversion_results`, `classify_reversion_verdict`. Also computes per-seed structural-jaccard metrics on the reverted prediction. |
| `bin/extract_passing.py` | Python | 678 | Reads the per-sequence aggregated CSV, applies tier classification (A/B/C/none) and the eligibility predicate, and emits `passing_summary.csv` (one canonical row per MPNN sequence). Imported by `cross_sequence_summary.py`; invoked once by `negative_steering_run_one.sh` per MPNN sequence. |
| `bin/cross_sequence_summary.py` | Python | 904 | Cohort aggregator. Combines every per-sequence `passing_summary.csv` + `aggregated_results.csv` into one cross-sequence CSV with tier ranking, composite-score sorting, fallback row picker for tier-none sequences, and work-dir → published-tree path rewriting. |
| `bin/sequence_registry.py` | Python | 175 | Hash-based sequence deduplication helper (`load_registry`, `save_registry`, `lookup`, `register`, `count_existing_seeds`, `deduplicate_sequences`). **Confirmed dead code — no Python file imports it; no .nf or .sh references it; pipeline_notes9.md already flagged it as "confirmed dead code".** |

### Negative-steering shell orchestrator + Rosetta XML

| Path | Type | LOC | Summary |
|---|---|---:|---|
| `bin/negative_steering_run_one.sh` | shell | 449 | The single-job orchestrator that NEGSTEER_RUN_ONE shells out to per MPNN sequence. Chains: `plan` → array of `predict-one` → `collect` → `kickoff-distances` → `kickoff-prefilter` → `build-contaminated` → `plan-reversions` → array of `predict-reversion-one` → `harvest-reversions` → `kickoff-finalize` → `aggregate` → `compute-final-metrics` → `aggregate-per-sequence` → `extract_passing.py`. Replaces the old SLURM-array submitter pattern; collapses the whole chain into one Nextflow task. |
| `bin/fastrelax_for_ia.xml` | Rosetta XML | 59 | RosettaScripts protocol for FastRelax (1 repeat, Cartesian, ref2015_cart) prior to InterfaceAnalyzer. Pinned for Bennett-2023 ΔΔG calibration. |

### Negative-steering metrics & post-processing

| Path | Type | LOC | Summary |
|---|---|---:|---|
| `bin/compute_metrics.py` | Python | 1657 | Boltz-output metrics calculator: ipSAE / iPTM / actifPTM / pLDDT / PAE-derived metrics, contact analysis, weighted-Jaccard, intact-core, interface-pLDDT. **Invoked via `subprocess.run` from boltz2_iterate_steering.py and reversion.py** (not via Python import); also imported directly by `compute_interface_metrics.py` and `derive_input_design_region.py`. |
| `bin/compute_interface_metrics.py` | Python | 739 | Per-row interface metrics on cohort-aggregated rows: iRMSD, fnat, DockQ (via DockQ python API), 15Å iPSAE, intact_core, weighted_jaccard. Run by `NEGSTEER_INTERFACE_METRICS`. |
| `bin/run_biophysical_metrics.py` | Python | 333 | One of three orthogonal-metrics streams. Computes BSA (FreeSASA), interface_plddt (cached Boltz B-factors), interface_hbonds (MDAnalysis). Run by `NEGSTEER_BIOPHYSICAL_METRICS`. |
| `bin/run_rosetta_metrics.py` | Python | 411 | One of three orthogonal-metrics streams. Two-stage Rosetta call: FastRelax → InterfaceAnalyzer; emits Sc and ΔΔG. Run by `NEGSTEER_ROSETTA_METRICS`. |
| `bin/parse_af3_output.py` | Python | 320 | One of three orthogonal-metrics streams. Parses AlphaFold3 output (per-seed/per-sample CIFs), pairs CA atoms by position with reference, computes receptor-aligned-effector RMSD and ipTM medians. Run by `AF3_PARSE_OUTPUT`. |
| `bin/merge_orthogonal_metrics.py` | Python | 204 | Joins the three orthogonal-metrics summaries onto the per-row CSV; computes `passes_orthogonal_filters` flag. Run by `NEGSTEER_ORTHOGONAL_METRICS`. |
| `bin/extract_survivor_manifest.py` | Python | 174 | Reads cross-sequence CSV + per-sequence workdirs, emits a fan-out manifest CSV (one row per survivor) for the orthogonal-metrics streams. Run by `EXTRACT_SURVIVOR_MANIFEST`. |
| `bin/build_control_sequences.py` | Python | 340 | Generates two negative controls (`scrambled` permutes the design region under deterministic seed; `polyA` substitutes alanines). Each goes through the same negsteer chain as a real MPNN sequence. Run by `NEGSTEER_CONTROLS`. |
| `bin/derive_design_region.py` | Python | 379 | Reads `rfdiffusion_metrics.json`, emits 1-based positional indices of the ProteinMPNN-designed region for one design. Used by `NEGSTEER_DERIVE_INDICES`. |
| `bin/derive_true_interface.py` | Python | 326 | Reads `rfdiffusion_metrics.json`, emits 0-based positional indices of the design's intended-binding-site contact residues for one design. Used by `NEGSTEER_DERIVE_INDICES`. |
| `bin/derive_input_design_region.py` | Python | 398 | Computes the design region + true-interface indices for the **input complex** (used by negative controls — controls have no RFDiffusion ancestor to read from). Run by `DERIVE_INPUT_INDICES`. |

### Negative-steering plot scripts (production copies of test scripts)

| Path | Type | LOC | Summary |
|---|---|---:|---|
| `bin/negsteer_plots.py` | Python | 1911 | Cohort-level negsteer plots (7 PNGs): tier landscape, seed-outcomes heatmap + bars, ra_eff vs jaccard, filter cascade, controls diagnostic, mutation impact. Run by `NEGSTEER_PLOTS`. **Verbatim lift of `tests/negative_steering/test_negsteer_plots.py`** minus a test-path fallback in `_resolve_csv_path`. |
| `bin/negsteer_within_sequence_plots.py` | Python | 987 | Within-sequence (per-seed) negsteer plots (5 PNGs): per-seed dispersion overview + grid, per-design stage trajectories, weighted-vs-true jaccard, composite-vs-confidence. Run by `NEGSTEER_WITHIN_SEQUENCE_PLOTS`. **Verbatim lift of test counterpart.** |
| `bin/orthogonal_metrics_plots.py` | Python | 1146 | Orthogonal-metrics cohort plots (4 PNGs): AF3 vs Boltz scatter, filter cascade, metrics-vs-composite multi-panel, combined cohort+orthogonal summary. Run by `ORTHOG_PLOTS`. **Verbatim lift of test counterpart.** |

### Receptor-resurfacing upstream (HADDOCK, RFDiffusion, Rosetta filtering, MPNN)

| Path | Type | LOC | Summary |
|---|---|---:|---|
| `bin/haddock3_prepare.py` | Python | 273 | Generates HADDOCK3 AIR restraints from de novo regions in the contig string. Run by `HADDOCK3_PREPARE`. |
| `bin/collect_haddock3_dock.py` | Python | 291 | Post-HADDOCK: parses CAPRI / clustfcc TSVs, picks best cluster, extracts best model PDB. Run by `HADDOCK3_DOCK`. |
| `bin/haddock3_plots.py` | Python | 697 | HADDOCK diagnostic plots: score-vs-BSA, cluster sizes, interface contact heatmap. Run by `HADDOCK3_PLOTS`. |
| `bin/haddock_utils.py` | Python | 283 | Shared helpers for HADDOCK scripts: parse_capri_tsv, parse_clustfcc_tsv, model_stem, extract_heavy_atoms, find_interface_residues, copy_pdb, etc. Imported by `collect_haddock3_dock.py`, `haddock3_plots.py`, `extract_hotspots.py`. |
| `bin/extract_hotspots.py` | Python | 207 | After HADDOCK best model selected, identifies receptor hotspot residues (with sequence-identity chain disambiguation). Run by `EXTRACT_HOTSPOTS`. |
| `bin/build_contigs.py` | Python | 237 | Builds RFDiffusion contigs string from receptor hotspots / chain mapping. Run by `BUILD_CONTIGS`. |
| `bin/contig_utils.py` | Python | 326 | Shared contig-parsing helpers: `parse_block_segments`, `remap_segments_to_pdb`, `resolve_contigs`, `get_expected_chain_lengths`, `parse_design_region`. Imported by `rfdiffusion_contigs.py`, `rfdiffusion_filter.py`. |
| `bin/rfdiffusion_contigs.py` | Python | 42 | Thin CLI wrapper around `contig_utils.resolve_contigs` — used by `RESOLVE_CONTIGS` (preprocessing). |
| `bin/rfdiffusion_filter.py` | Python | 926 | Post-RFDiffusion filter. Splits chain-break design PDBs back into receptor+effector chains (hardcoded A=receptor, B=effector); computes per-design metrics (contacts, COM displacement, per-region clustering); emits `rfdiffusion_metrics.json`. Run by `RFDIFFUSION_FILTER`. |
| `bin/rfdiffusion_plots.py` | Python | 1248 | Five RFDiffusion plots: contact map, design clustering, specificity coverage, design lengths, COM displacement. Run by `RFDIFFUSION_PLOTS`. |
| `bin/rosetta_filter_collect.py` | Python | 227 | Per-design Rosetta scorefile parser; aggregates Sc/dG_separated/dSASA_int across all RFDiffusion designs and filters by Sc threshold. Run by `ROSETTA_FILTER`. |
| `bin/rosetta_filter_plots.py` | Python | 206 | Sc histogram + filter-region plot. Run by `ROSETTA_FILTER_PLOTS`. |
| `bin/pipeline_correct_sequences.py` | Python | 770 | Reconstructs full chimeric receptor sequences from MPNN's per-region output (via anchor-based alignment to the native sequence); also generates fixed-positions JSONL for ProteinMPNN. Run by `SEQUENCE_CORRECTION`. **Cross-pipeline file** — patched in negsteer sessions because its `split_sequence` bug broke RFDiffusion-design test runs. |
| `bin/mpnn_sequence_qc.py` | Python | 195 | Filters MPNN sequences on poly-X runs and unusual amino-acid composition. Run by `SEQUENCE_QC`. |
| `bin/mpnn_design_region_score.py` | Python | 319 | Re-scores each MPNN sequence using the ProteinMPNN model on the design region only (vs the global score MPNN emits). Run by `MPNN_DESIGN_REGION_SCORE`. |
| `bin/mpnn_cluster_sequences.py` | Python | 175 | MMseqs2 sequence clustering on MPNN designs. Run by `MPNN_CLUSTER`. |
| `bin/mpnn_select_top.py` | Python | 72 | Selects top-N designs after QC + clustering. Run by `MPNN_SELECT_TOP`. |
| `bin/mpnn_plots.py` | Python | 550 | Four MPNN plots: score distribution (design-region top, global bottom), AA composition with native overlay, sequence diversity, physicochem (KD hydrophobicity + net charge). Run by `MPNN_PLOTS`. |

## modules/ — per-stage Nextflow modules

| Path | Type | LOC | Summary |
|---|---|---:|---|
| `modules/preprocessing.nf` | Nextflow | 161 | Processes: `EXTRACT_SEQUENCES`, `RESOLVE_CONTIGS`, `WRITE_DUMMY_MAPPING`. PDB-only entry-point preprocessing. |
| `modules/haddock.nf` | Nextflow | 287 | Processes: `HADDOCK3_PREPARE`, `HADDOCK3_DOCK`, `HADDOCK3_PLOTS`, `EXTRACT_HOTSPOTS`, `BUILD_CONTIGS`. Branch-A docking flow. |
| `modules/rfdiffusion.nf` | Nextflow | 186 | Processes: `RFDIFFUSION`, `RFDIFFUSION_FILTER`, `RFDIFFUSION_PLOTS`. Backbone generation + per-design metrics. |
| `modules/rosetta_filtering.nf` | Nextflow | 165 | Processes: `ROSETTA_SC`, `ROSETTA_FILTER`, `ROSETTA_FILTER_PLOTS`. Pre-MPNN physics filter on Sc. |
| `modules/proteinmpnn.nf` | Nextflow | 330 | Processes: `MPNN_FIXED_POSITIONS`, `PROTEINMPNN`, `SEQUENCE_CORRECTION`, `SEQUENCE_QC`, `MPNN_DESIGN_REGION_SCORE`, `MPNN_CLUSTER`, `MPNN_SELECT_TOP`, `MPNN_PLOTS`. Sequence design + QC + selection. |
| `modules/negative_steering.nf` | Nextflow | 424 | Processes: `NEGSTEER_DERIVE_INDICES`, `NEGSTEER_RUN_ONE`, `NEGSTEER_CROSS_SEQUENCE`, `NEGSTEER_PLOTS`, `NEGSTEER_WITHIN_SEQUENCE_PLOTS`. The single-job-per-MPNN-sequence negsteer flow. |
| `modules/negsteer_controls.nf` | Nextflow | 222 | Processes: `DERIVE_INPUT_INDICES`, `NEGSTEER_CONTROLS`. Generates and runs negsteer for the two control variants (`scrambled`, `polyA`). |
| `modules/negsteer_manifest.nf` | Nextflow | 66 | Process: `EXTRACT_SURVIVOR_MANIFEST`. Single CPU step that fans the survivors out to the orthogonal-metrics streams. |
| `modules/negsteer_interface_metrics.nf` | Nextflow | 69 | Process: `NEGSTEER_INTERFACE_METRICS`. iRMSD/fnat/DockQ/15Å iPSAE/intact_core/weighted_jaccard. |
| `modules/negsteer_af3_nomsa.nf` | Nextflow | 280 | Processes: `AF3_SETUP_DB`, `AF3_NOMSA_ON_SURVIVORS`, `AF3_PARSE_OUTPUT`. AlphaFold3-no-MSA orthogonal validation. |
| `modules/negsteer_biophysical_metrics.nf` | Nextflow | 53 | Process: `NEGSTEER_BIOPHYSICAL_METRICS`. BSA/interface_plddt/interface_hbonds. |
| `modules/negsteer_rosetta_metrics.nf` | Nextflow | 67 | Process: `NEGSTEER_ROSETTA_METRICS`. FastRelax + InterfaceAnalyzer → Sc + ΔΔG. |
| `modules/negsteer_orthogonal_metrics.nf` | Nextflow | 143 | Processes: `NEGSTEER_ORTHOGONAL_METRICS`, `ORTHOG_PLOTS`. Joins the three streams + emits cohort plots. |

## tests/ — per-module test harnesses

### tests/full_test_run/

| Path | Type | LOC | Summary |
|---|---|---:|---|
| `tests/full_test_run/params_test.yml` | YAML | 165 | Parameter file for a full end-to-end pipeline test run. (Data file `rosetta_summary.csv` skipped.) |

### tests/haddock/

| Path | Type | LOC | Summary |
|---|---|---:|---|
| `tests/haddock/test_haddock.nf` | Nextflow | 127 | Module-test workflow that runs HADDOCK3_PREPARE → HADDOCK3_DOCK → HADDOCK3_PLOTS → EXTRACT_HOTSPOTS → BUILD_CONTIGS on small fixed inputs. |
| `tests/haddock/run_test_haddock.slurm.sh` | shell | 71 | SLURM wrapper to submit `test_haddock.nf`. |

### tests/rfdiffusion/

| Path | Type | LOC | Summary |
|---|---|---:|---|
| `tests/rfdiffusion/test_rfdiffusion.nf` | Nextflow | 101 | Module-test workflow that runs RFDIFFUSION → RFDIFFUSION_FILTER → RFDIFFUSION_PLOTS on a small fixed input. |
| `tests/rfdiffusion/run_test_rfdiffusion.slurm.sh` | shell | 71 | SLURM wrapper to submit `test_rfdiffusion.nf`. |
| `tests/rfdiffusion/test_rfdiffusion_plots.py` | Python | 1163 | Standalone iteration harness for `bin/rfdiffusion_plots.py`. Imports the prod module for three "good" plots; redefines the two iteration-target plots locally. |
| `tests/rfdiffusion/run_test_rfdiffusion_plots.slurm.sh` | shell | 108 | SLURM wrapper for the standalone-plot iteration script. |

### tests/rosetta_filtering/

| Path | Type | LOC | Summary |
|---|---|---:|---|
| `tests/rosetta_filtering/test_rosetta_filtering.nf` | Nextflow | 123 | Module-test workflow for the ROSETTA_SC → ROSETTA_FILTER → ROSETTA_FILTER_PLOTS chain. |
| `tests/rosetta_filtering/run_test_rosetta_filtering.slurm.sh` | shell | 71 | SLURM wrapper. |
| `tests/rosetta_filtering/test_rosetta_filtering_plots.py` | Python | 224 | Standalone iteration harness for `bin/rosetta_filter_plots.py`. |
| `tests/rosetta_filtering/run_test_rosetta_filtering_plots.slurm.sh` | shell | 84 | SLURM wrapper for the iteration script. |

### tests/proteinmpnn/

| Path | Type | LOC | Summary |
|---|---|---:|---|
| `tests/proteinmpnn/test_proteinmpnn.nf` | Nextflow | 232 | Module-test workflow for the full ProteinMPNN sub-chain (FIXED_POSITIONS → PROTEINMPNN → SEQUENCE_CORRECTION → SEQUENCE_QC → DESIGN_REGION_SCORE → CLUSTER → SELECT_TOP → PLOTS). |
| `tests/proteinmpnn/run_test_proteinmpnn.slurm.sh` | shell | 65 | SLURM wrapper. |
| `tests/proteinmpnn/test_mpnn_plots.py` | Python | 625 | Standalone iteration harness for `bin/mpnn_plots.py`. |
| `tests/proteinmpnn/run_test_mpnn_plots.slurm.sh` | shell | 92 | SLURM wrapper for the iteration script. |

### tests/negative_steering/

| Path | Type | LOC | Summary |
|---|---|---:|---|
| `tests/negative_steering/test_negative_steering.nf` | Nextflow | 433 | Module-test workflow including controls and both plot processes. Defines a local `EMPTY_DESIGN_REGION_PLACEHOLDER` process for when controls are off. |
| `tests/negative_steering/run_test_negative_steering.slurm.sh` | shell | 84 | SLURM wrapper. |
| `tests/negative_steering/test_negsteer_plots.py` | Python | 1916 | Standalone iteration harness for cohort negsteer plots. **Source of truth for `bin/negsteer_plots.py` — production is a verbatim copy minus header + test-path fallback.** |
| `tests/negative_steering/run_test_negsteer_plots.slurm.sh` | shell | 170 | SLURM wrapper for the cohort plot iteration script. |
| `tests/negative_steering/test_negsteer_within_sequence_plots.py` | Python | 981 | Standalone iteration harness for within-sequence negsteer plots. **Source of truth for `bin/negsteer_within_sequence_plots.py`.** |
| `tests/negative_steering/run_test_negsteer_within_sequence_plots.slurm.sh` | shell | 124 | SLURM wrapper. |

### tests/orthogonal_metrics/

| Path | Type | LOC | Summary |
|---|---|---:|---|
| `tests/orthogonal_metrics/test_orthogonal_metrics.nf` | Nextflow | 223 | Module-test workflow for the orthogonal-metrics fan-out (manifest → AF3/biophysical/Rosetta in parallel → merge → plots). |
| `tests/orthogonal_metrics/run_test_orthogonal_metrics.slurm.sh` | shell | 97 | SLURM wrapper. |
| `tests/orthogonal_metrics/merge_orthogonal_metrics.py` | Python | 216 | **Second copy of `bin/merge_orthogonal_metrics.py`. The two have diverged (see `05_findings.md`)** — the test version implements the AF3-flag-only demotion that pipeline_notes13 says was applied to merge_orthogonal_metrics.py; the production version still gates on AF3. |
| `tests/orthogonal_metrics/test_orthogonal_metrics_plots.py` | Python | 1170 | Standalone iteration harness for `bin/orthogonal_metrics_plots.py`. **Source of truth for the production copy.** |
| `tests/orthogonal_metrics/run_test_orthogonal_metrics_plot.slurm.sh` | shell | 135 | SLURM wrapper. |

## Totals

- 38 Python files in `bin/` (~28,500 LOC) + 1 shell + 1 XML
- 13 Nextflow modules in `modules/` (~2,700 LOC)
- 7 Python test scripts in `tests/*/` (~7,100 LOC) + 6 test `.nf` workflows + 11 SLURM wrappers
- 1 top-level `main.nf` (1,061 LOC), `nextflow.config`, slurm wrapper, params example
- Total: ~42,000 LOC across the source tree
