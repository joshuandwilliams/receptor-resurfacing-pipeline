# Example pipeline outputs — reference set for characterization tests

This directory holds a hand-picked sample of files produced by an end-to-end run of the pipeline. They are kept here as **reference outputs for the Phase 2 characterization tests** of the codebase remediation: each test will compare a fresh run's output against the corresponding file here and flag any divergence.

The cohort underlying these files is the `tests/full_test_run/params_test.yml` test run that was completed on 2026-04-28/29 (paths embedded in the JSONs trace back to `/hpc-home/jowillia/receptor_design/receptor-resurfacing-pipeline/tests/full_test_run/`). Branch B (pre-docked complex) was the entry path; Branch A (HADDOCK docking) was not exercised.

## Important caveat — these are *current behaviour*, not *verified-correct behaviour*

Characterization tests built from this reference set will detect **changes** in pipeline behaviour. They will not detect logical errors that are already present in these outputs. Several files in this directory may be subtly or grossly wrong — undetected bugs in the producing scripts, miscalibrated thresholds, or the inverted `scaffold` semantics flagged in the glossary (§F3) will all serialize directly into these CSVs and JSONs without any signal that something is off.

Use this set as a **pin** during refactoring, not a **proof of correctness**. Anything in the inventory's `05_findings.md` or the glossary's §10 cleanup list should be treated as suspect until the post-cleanup outputs are validated independently.

The negative-steering subdirectory is a particularly load-bearing example of this caveat: many of the files there are intermediates inside one MPNN sequence's per-cycle workdir (`negative_steering/runs/<seq_name>/cycle_0/`). They reflect what the orchestrator *currently* writes at each phase, not a reviewed schema.

---

## Directory layout

```
example_output_files/
├── README.md                       (this file)
├── dag.html, pipeline_report.html, timeline.html, trace.txt   (Nextflow execution reports)
├── rosetta_summary.csv             (anomalous top-level — see notes)
├── preprocessing/
├── rfdiffusion/
├── rosetta_filtering/
├── sequences/
├── negative_steering/
├── orthogonal_metrics/
└── plots/
```

Stages map onto the categories in `notes/inventory/04_functional_categorization.md`. Branch A (HADDOCK) outputs are not present — see the "What's missing" section.

---

## Top-level Nextflow execution reports

These four files are emitted by Nextflow itself, not by any pipeline script. They are useful as a record of *which* tasks ran and how long they took, but they are not characterization-test targets — they contain timestamps and workdir hashes that change every run.

| File | Producer | Notes |
|---|---|---|
| `dag.html` | Nextflow `-with-dag` | Workflow DAG visualisation. |
| `pipeline_report.html` | Nextflow `-with-report` | Per-task resource usage, status, exit codes. |
| `timeline.html` | Nextflow `-with-timeline` | Gantt-style timeline of task execution. |
| `trace.txt` | Nextflow `-with-trace` | TSV: one row per task with hash, status, duration, RSS, etc. Header on row 1; columns include `task_id`, `hash`, `name`, `status`, `duration`, `peak_rss`, `cpus`, `memory`. |

## Top-level `rosetta_summary.csv` (anomalous placement)

Schema: `seq_name,sc,rosetta_ddg,rosetta_failures` — single data row `design_0_seq_2,,,fastrelax_xml_not_found`.

This is a per-survivor Rosetta-orthogonal-metrics summary emitted by `NEGSTEER_ROSETTA_METRICS` (`bin/run_rosetta_metrics.py` via `modules/negsteer_rosetta_metrics.nf`). In a normal run it would be published to `orthogonal_metrics/rosetta/<seq_name>/rosetta_summary.csv` (one per survivor) and then merged into the cohort `orthogonal_metrics/rosetta_summary.csv`. The presence of this **failure-case row** at the top level appears intentional — likely a hand-picked example showing what the file looks like when the FastRelax XML can't be located. **Flagged for confirmation.**

---

## Stage 1 — Preprocessing  (`preprocessing/`)

Producer: `modules/preprocessing.nf` (Branch B entry path).

| File | Producing process | Producing script | Description |
|---|---|---|---|
| `sequences.json` | `EXTRACT_SEQUENCES` | inline Python in `modules/preprocessing.nf` | Receptor + effector amino-acid sequences and PDB residue numbers extracted from the input complex. Schema: `{receptor_seq, effector_seq, receptor_len, effector_len, receptor_resnums[], effector_resnums[], receptor_start_pdb}`. |
| `processed_contigs.txt` | `RESOLVE_CONTIGS` | `bin/rfdiffusion_contigs.py` (CLI wrapper around `bin/contig_utils.py:resolve_contigs`) | Single-line resolved RFDiffusion contig string. Example here: `A1-32/10-20/A46-72/6-6 C32-113`. Anchor segments use `A`/`C` chain prefixes; bare ranges are de novo regions. |

Not present: the dummy mapping JSON written by `WRITE_DUMMY_MAPPING` (a 2-byte `{}` placeholder) — trivially reproducible.

---

## Stage 2 — RFDiffusion backbone generation  (`rfdiffusion/`)

Producer: `modules/rfdiffusion.nf` → `RFDIFFUSION_FILTER` (`bin/rfdiffusion_filter.py`).

| File | Description |
|---|---|
| `rfdiffusion_metrics.json` | Per-design metrics for all 32 RFDiffusion outputs. Top-level key `designs` is an array; each entry has `design` (PDB filename), `n_receptor_residues`, `n_contact_pairs`, `n_unique_receptor_residues`, `n_unique_effector_residues`, `n_contacts_in_design`, `n_contacts_outside_design`, `frac_contacts_in_design`, `design_region_coverage`, `scaffold_rmsd` (note: per glossary §F3 this value is actually *motif* RMSD — the field name is inverted), `endpoint_distance_input`, `endpoint_distance_design`, `design_region_com_displacement`, `receptor_contact_residues[]`, `effector_contact_residues[]`, `contact_pairs[][]`, and `design_region_coords[][3]` (Cα coordinates of the de novo region). Consumed downstream by all three `bin/derive_*.py` scripts and by `NEGSTEER_DERIVE_INDICES`. |
| `filter_summary.json` | One-blob summary of the filter pass: `n_total`, `n_passing`, `n_filtered`, `min_hotspot_frac`, `filtered_designs[]`. Here all 32 designs pass. |
| `passing_designs.txt` | Newline-separated list of passing design PDB filenames (`design_0.pdb` … `design_31.pdb`). |

Not present: the actual `design_*.pdb` backbone files, the `passing/` and `split/` PDB subdirectories, and any RFDiffusion trajectory output (`traj/`).

---

## Stage 3 — Rosetta filtering (pre-MPNN)  (`rosetta_filtering/`)

Producer: `modules/rosetta_filtering.nf` → `ROSETTA_FILTER` (`bin/rosetta_filter_collect.py`). Pre-MPNN physics filter on shape complementarity (Sc).

| File | Description |
|---|---|
| `rosetta_filter_metrics.json` | Per-design Rosetta InterfaceAnalyzer scores. `designs[]` with fields `design`, `design_stem`, `sc_value`, `dG_separated`, `dSASA_int`, `dG_dSASA_density`, `packstat`, `delta_unsatHbonds`, `nres_int`, `passes_filter`. |
| `rosetta_filter_summary.json` | Pass/fail counts and the threshold used: `{n_total, n_passing, n_filtered, sc_threshold, filtered_designs[]}`. Here `sc_threshold=0.5`, all 32 pass. |
| `rosetta_passing_designs.txt` | Newline-separated list of passing PDBs. |

Not present: per-design Rosetta scorefiles (`interface_scores_*.sc` from `ROSETTA_SC`), or the `passing/*.pdb` symlinked subdirectory.

---

## Stage 4 — ProteinMPNN sequence design  (`sequences/`)

Producer: `modules/proteinmpnn.nf`. The files here cover the full SEQUENCE_CORRECTION → SEQUENCE_QC → MPNN_DESIGN_REGION_SCORE → MPNN_CLUSTER → MPNN_SELECT_TOP chain.

| File | Producing process | Producing script | Description |
|---|---|---|---|
| `mpnn_corrected.fasta` | `SEQUENCE_CORRECTION` | `bin/pipeline_correct_sequences.py` | All MPNN sequences, anchor-aligned and reconstructed into full chimeric receptor sequences. Each entry's header carries `mpnn=<score>`, `designed=<de-novo-residues>`, `design_region_length_observed=<obs>`, `design_region_spec=<from-contig>`, `changes=<num-changes>`. The body line is `<receptor>:<effector>` colon-joined. |
| `sequence_metadata.csv` | `SEQUENCE_CORRECTION` | `bin/pipeline_correct_sequences.py` | Tabular form of the same. Columns: `design`, `seq`, `mpnn_score`, `native_residues`, `designed_residues`, `num_changes`, `design_region_length_observed`, `design_region_spec`, `corrected_receptor`, `receptor_total_length`. |
| `qc_metadata.csv` | `SEQUENCE_QC` | `bin/mpnn_sequence_qc.py` | Same schema as `sequence_metadata.csv`, filtered to QC-passing sequences. (Currently this run's QC removed 4 duplicates and zero failures, so the row schema is identical.) |
| `qc_report.txt` | `SEQUENCE_QC` | `bin/mpnn_sequence_qc.py` | Human-readable QC summary: total in/out, duplicates removed (with the kept sequence and its score), failures. |
| `scored_metadata.csv` | `MPNN_DESIGN_REGION_SCORE` | `bin/mpnn_design_region_score.py` | `qc_metadata.csv` plus a `design_region_score` column (ProteinMPNN re-score on the design region only — distinct from the global `mpnn_score`). |
| `mpnn_cluster_counts.csv` | `MPNN_CLUSTER` | `bin/mpnn_cluster_sequences.py` | MMseqs2 cluster counts at sweep of identity thresholds. Columns: `threshold,n_clusters,n_sequences`. Diagnostic only — not used downstream. |
| `top_metadata.csv` | `MPNN_SELECT_TOP` | `bin/mpnn_select_top.py` | The top-N sequences after QC + clustering + scoring. Same schema as `scored_metadata.csv`. These are the sequences that go into the negative-steering stage. |
| `design_0_seq_0.fasta` | `SEQUENCE_CORRECTION` (publishes `af2_fastas/`) | `bin/pipeline_correct_sequences.py` | One per-sequence FASTA from the `af2_fastas/` directory: two records `>receptor` and `>effector`, full chains. This is the input file `NEGSTEER_RUN_ONE` consumes per MPNN sequence. Only one example included; in a real run there is one per `(design, sequence)` pair. |

Not present: the `af2_fastas/`, `qc_fastas/`, `top_fastas/` directories themselves; the raw `mpnn_results/` subdirectory; the fixed-positions JSONL (`MPNN_FIXED_POSITIONS`); the `selection_report.txt` from `MPNN_SELECT_TOP`.

---

## Stage 5 — Negative-steering validation  (`negative_steering/`)

This is the most complex stage and the example set here reflects that. **Most of the files in this subdirectory are intermediates that normally live inside one MPNN sequence's per-cycle workdir** — i.e. they would normally be at `negative_steering/runs/<seq_name>/cycle_0/<file>` — but here they have been hand-picked from a single representative sequence (likely `design_0_seq_0`, judging by the embedded paths) and placed flat. The two cohort-level files (`cross_sequence_summary*.csv`) are at the cohort level in a real run.

The producing scripts are largely subcommands of `bin/boltz2_negative_steering.py` and `bin/boltz2_iterate_steering.py`, sequenced by the shell orchestrator `bin/negative_steering_run_one.sh` and run by the `NEGSTEER_RUN_ONE` Nextflow process. The cohort step is `NEGSTEER_CROSS_SEQUENCE` (`bin/cross_sequence_summary.py`); the interface-metrics extension is `NEGSTEER_INTERFACE_METRICS` (`bin/compute_interface_metrics.py`).

### 5a. Plan phase outputs

| File | Producing subcommand | Description |
|---|---|---|
| `plan.json` | `boltz2_negative_steering.py plan` | The full plan blob for one MPNN sequence's cold-start: ground-truth PDB path, chain conventions (`pred_receptor_chain="A"`, `pred_effector_chain="B"`), wild-type receptor sequence, initial RMSD metrics (`initial_receptor_aligned_effector_rmsd`, `initial_independent_receptor_rmsd`, `initial_independent_effector_rmsd`), RMSD/contact thresholds, steering `mode` (here `mild`), candidate-pool size, `n_designs`, `num_seeds`, container path, recycling/diffusion-sample counts, and the `designs[]` array — one entry per steered design × seed combination, each with `index`, `name` (e.g. `design_00_s0`), `dir`, `yaml`, `total_mutations`, `mutated_positions[]`, `sequence_group`, `seed_index`, `boltz_seed`. |
| `wrong_interface_residues.txt` | `boltz2_negative_steering.py plan` | TSV of residues at the cold-start (wrong) predicted interface. Header comment `# 1-based receptor residue, wt one-letter, min_heavy_dist_to_eff (Å), in_true_site, in_protected_set, surface_exposed, in_pool` followed by one row per residue; the candidate pool for steering is the subset where `in_pool=1`. |
| `true_interface_residues.txt` | `boltz2_negative_steering.py plan` (or `NEGSTEER_DERIVE_INDICES` upstream) | TSV of residues in the true (intended) binding site. Header `# 1-based receptor residue, wt one-letter, min_heavy_dist_to_eff (Å) — protected from mutation`. (Note: 0-based vs 1-based convention conflict flagged in glossary §F4 — the protected-set logic re-bases internally.) |

### 5b. Steered-prediction collect phase

| File | Producing subcommand | Description |
|---|---|---|
| `steered_results.csv` | `boltz2_negative_steering.py collect` | Per-design ranked steered-prediction summary (one row per `(design, sequence_group, seed_index)`). Columns: `rank`, `design`, `sequence_group`, `seed_index`, `total_mutations`, `receptor_aligned_effector_rmsd`, `delta_receptor_aligned_effector_rmsd`, `independent_receptor_rmsd`, `independent_effector_rmsd`, `receptor_intact`, `wrong_jaccard`, `n_shared_wrong`, `true_jaccard`, `n_shared_true`, `n_design_interface_residues`, `status`, `pdb`. |
| `steered_results_aggregate.csv` | `boltz2_negative_steering.py collect` | Per-`sequence_group` aggregation across seeds: `sequence_group`, `n_seeds`, `total_mutations`, `representative_design`, plus mean/std/min/max for the three RMSD families and `n_intact_seeds`, `any_seed_intact`. |
| `summary.txt` | `boltz2_negative_steering.py collect` | Human-readable run summary: ground truth, chain conventions, contact cutoff, steering mode, max mutations, candidate pool size, seeds per sequence, n unique sequences, n predicted, n intact, the residue-set sizes (true site, wrong site, protected, pool), initial RMSDs, and best-overall + best-with-intact-proteins designs. |

### 5c. Multi-cycle/kickoff intermediates  (single-cycle in the current pipeline)

| File | Producing subcommand | Description |
|---|---|---|
| `kickoff_distances.json` | `boltz2_iterate_steering.py kickoff-distances` | Pose-novelty distances from each candidate to every upstream pose. Schema: `{n_generated, candidates[]}` where each candidate has `design`, `design_idx`, `sequence_group`, `seed_index`, `pdb`, the three RMSD metrics, `dists_to_upstream[]`, `min_dist_to_upstream`, and a nested `result{}` block with the original collect-phase row. |
| `prefilter.json` | `boltz2_iterate_steering.py kickoff-prefilter` | Same shape as `kickoff_distances.json` filtered to `intact_candidates[]`; adds `cycle`, `parent_pathway`, `n_generated`, `n_intact`. |
| `passing.json` | `boltz2_iterate_steering.py kickoff-finalize` | Selection of poses for cycle-N+1: `cycle`, `parent_pathway`, `novelty_cutoff`, `max_passing`, `n_generated`, `n_intact`, `n_novel`, `n_selected`, the three verdict counters, `reverted_labels[]`, `dropped_labels[]`, `all_candidates[]`. In single-cycle runs (current pipeline) this drives the cohort summary; in multi-cycle it drives the next cycle's planning. |
| `pathways.json` | `boltz2_iterate_steering.py aggregate` | Pathway lineage list: `[{label, cycle, workdir, n_designs}, ...]`. Currently only `cycle_0` since multi-cycle is dormant. |
| `cycle_statistics.csv` | `boltz2_iterate_steering.py aggregate` | One-row-per-cycle cohort stats: `cycle, pathway, n_generated, n_intact, n_novel, n_selected, min_dist_to_upstream_values, receptor_aligned_effector_rmsd_vs_truth_values`. |

### 5d. Reversion phase

The reversion pass identifies steering mutations inside the protected set and tries to revert them back to the MPNN sequence's residues. In this example the cold-start had no contaminating mutations, so the reversion machinery emitted "empty" outputs.

| File | Producing subcommand | Description |
|---|---|---|
| `contaminated.json` | `boltz2_iterate_steering.py build-contaminated` | List of mutations to revert. Schema: `{cycle, parent_pathway, n_checked, n_contaminated, contaminated[]}`. Empty here (`n_contaminated=0`). |
| `reversion_plan.json` | `boltz2_iterate_steering.py plan-reversions` | Per-design reversion YAMLs to predict. Schema: `{workdir, n_contaminated, n_staged, pred_receptor_chain, pred_effector_chain, entries[]}`. Empty here. |
| `reversion_results.json` | `boltz2_iterate_steering.py harvest-reversions` (per-design) | Per-design reverted-prediction harvest. Empty `{}` here. |
| `reversion_results_per_seed.json` | `boltz2_iterate_steering.py harvest-reversions` (per-seed) | Per-seed reverted-prediction harvest. Empty `{}` here. |

(Schemas for the populated case are not represented in this example set — flagged under "What's missing" below.)

### 5e. Aggregation and per-sequence finalisation

| File | Producing subcommand | Description |
|---|---|---|
| `raw_per_seed_results.csv` | `boltz2_iterate_steering.py aggregate` | The flat per-prediction CSV for one MPNN sequence. ~70 columns split into a `steered_*` family and a `reverted_*` family covering: mutation strings (chimerax + AA forms), the RMSD triple, `wrong_jaccard`, `true_jaccard`, `n_shared_wrong`, `n_shared_true`, contact-residue counts and lists, `mutated_contact_positions`, all the Boltz confidence metrics (`avg_plddt`, `complex_plddt`, `ptm`, `iptm`, `pae_mean`, `ipae`, `pae_pass_frac`, `ipsae_*`, `actifptm`, `interface_plddt`, `weighted_jaccard`), `intact_core`, `confidence_flag`, plus `min_dist_to_upstream`, `selected_for_next_cycle`, ranks, and `reversion_verdict`. The `steered_*` columns are populated for every row; the `reverted_*` columns are populated only if reversion ran for that seed. |
| `all_results_multicycle_with_metrics.csv` | `boltz2_iterate_steering.py compute-final-metrics` | `raw_per_seed_results.csv` schema, recomputed with the full extended metrics set. The basis the per-sequence aggregator reduces. |
| `aggregated_results.csv` | `boltz2_iterate_steering.py aggregate-per-sequence` | Per-`sequence_group` aggregate. ~140 columns: every metric from `raw_per_seed_results.csv` reappears as `<metric>_median`, `<metric>_mad`, `<metric>_n_used` triples; categorical columns become `<col>_majority` + `<col>_per_seed_counts`; counters `n_seeds_pose_holds`/`pose_collapses`/`new_contamination`/`no_data`/`clean_steered`; `aggregated_verdict` (per glossary §5: one of `pose_holds`/`pose_collapses`/`new_contamination`/`no_reversion`/`no_data`); `aggregated_verdict_reason`. |
| `passing_summary.csv` | `bin/extract_passing.py` | One canonical row per MPNN sequence after tier-A/B/C/none classification and the eligibility predicate. Filtered to `PASSING_VERDICTS = {"", "no_reversion", "pose_holds"}`. **Empty body in this example — the header is present but no rows survived for the picked sequence under `PASSING_VERDICTS`. This is unusual given the cohort summary has 4 tier-A rows; the file may have been picked from a sequence that did not pass.** Flagged for confirmation. |
| `run_one_runtime_sec.txt` | `bin/negative_steering_run_one.sh` | Single integer: per-MPNN-sequence wallclock seconds. Threaded into the cross-sequence CSV's `run_one_runtime_sec` column. |

### 5f. Cohort-level outputs

| File | Producing process | Producing script | Description |
|---|---|---|---|
| `cross_sequence_summary.csv` | `NEGSTEER_CROSS_SEQUENCE` | `bin/cross_sequence_summary.py` | The cohort table. One row per MPNN sequence (plus the two negative controls in a full run; not present here). 50 columns: cohort identification (`mpnn_sequence`, `row_type`, `cross_rank_by_composite`, `cross_rank_by_ra_eff`, `cross_tier`, `cross_composite_score`, `within_sequence_rank`, tier counts), source pointer (`source_passing_summary`), and a wide `representative_*` block carrying the chosen representative seed's metrics, mutation strings, canonical PDB path, all Boltz confidence metrics, the intact-core/confidence flags, and the per-sequence runtime. |
| `cross_sequence_summary_with_interface_metrics.csv` | `NEGSTEER_INTERFACE_METRICS` | `bin/compute_interface_metrics.py` | `cross_sequence_summary.csv` plus four trailing columns: `irmsd`, `fnat`, `dockq`, `dockq_failures`. The base for the orthogonal-metrics fan-out. |

---

## Stage 6 — Orthogonal-metrics validation  (`orthogonal_metrics/`)

Producer: the four per-stream Nextflow modules (`negsteer_manifest.nf`, `negsteer_af3_nomsa.nf`, `negsteer_biophysical_metrics.nf`, `negsteer_rosetta_metrics.nf`) merged by `negsteer_orthogonal_metrics.nf`.

| File | Producing process | Producing script | Description |
|---|---|---|---|
| `survivor_manifest.csv` | `EXTRACT_SURVIVOR_MANIFEST` | `bin/extract_survivor_manifest.py` | Fan-out manifest: one row per cohort row that has a `representative_canonical_pdb` on disk. Columns: `seq_name`, `canonical_pdb_abs`, `ground_truth_abs`, `effector_template_cif_abs`, `receptor_seq`, `effector_seq`. **Note**: per glossary §F1 the term "survivor" overstates what this manifest filters on — it is currently *any* row with a representative PDB, not a tier-passing row. |
| `af3_nomsa_summary.csv` | `AF3_PARSE_OUTPUT` | `bin/parse_af3_output.py` | Per-sequence AF3-no-MSA summary. Columns: `seq_name`, `af3_nomsa_best_ra_eff`, `af3_nomsa_mean_ra_eff`, `af3_nomsa_best_iptm`, `af3_nomsa_mean_iptm`, `af3_nomsa_n_correct_interface`, `af3_nomsa_total_predictions`, `af3_nomsa_failures`. (The per-stream summary is a single row in the cohort-merged file at this path; in the per-survivor working tree it would live under `orthogonal_metrics/af3/<seq_name>/`.) |
| `biophysical_summary.csv` | `NEGSTEER_BIOPHYSICAL_METRICS` | `bin/run_biophysical_metrics.py` | Per-sequence BSA + interface H-bonds. Columns: `seq_name`, `bsa`, `interface_hbonds`, `biophysical_failures`. Note `interface_plddt` (named in the inventory) does not appear in this row — flagged for confirmation; may be merged into the cohort row elsewhere. |
| `rosetta_summary.csv` | `NEGSTEER_ROSETTA_METRICS` | `bin/run_rosetta_metrics.py` | Per-sequence Rosetta FastRelax → InterfaceAnalyzer. Columns: `seq_name`, `sc`, `rosetta_ddg`, `rosetta_failures`. |
| `survivors_with_orthogonal_metrics.csv` | `NEGSTEER_ORTHOGONAL_METRICS` | `bin/merge_orthogonal_metrics.py` | The final cohort table. `cross_sequence_summary_with_interface_metrics.csv` plus all per-stream summary columns plus two trailing columns: `orthogonal_flags` (semicolon-joined predicate strings, e.g. `af3_nomsa_ra_eff_too_high:17.23`) and `passes_orthogonal_filters` (0/1). **Note**: per `notes/inventory/05_findings.md` A4, the production `bin/merge_orthogonal_metrics.py` and the test-tree copy have diverged on whether AF3 contributes to the final filter; this CSV reflects whichever version actually ran. Treat the `passes_orthogonal_filters` column as the canonical record of *what was computed*, not *what is correct*. |

---

## Stage 7 — Plots  (`plots/`)

All PNGs from the four plot scripts plus the two pre-MPNN ones, copied here verbatim. These are the supervisor-demo-facing outputs.

| Filename pattern | Producing process | Producing script | Description |
|---|---|---|---|
| `rfdiff_*.png` (5 files) | `RFDIFFUSION_PLOTS` | `bin/rfdiffusion_plots.py` | `contact_map`, `design_clustering`, `specificity_coverage`, `design_lengths`, `com_displacement`. |
| `rosetta_sc_histogram.png` | `ROSETTA_FILTER_PLOTS` | `bin/rosetta_filter_plots.py` | Sc score distribution + filter cutoff. (The script normally also emits a filter-region plot; only the histogram is included here.) |
| `mpnn_*.png` (4 files) | `MPNN_PLOTS` | `bin/mpnn_plots.py` | `score_distribution` (design-region top, global bottom), `aa_composition`, `sequence_diversity`, `physicochem`. |
| `negsteer_tier_landscape.png`, `negsteer_seed_outcomes_heatmap.png`, `negsteer_seed_outcomes_bars.png`, `negsteer_ra_eff_vs_jaccard.png`, `negsteer_filter_cascade.png`, `negsteer_controls_diagnostic.png`, `negsteer_mutation_impact.png` (7 files) | `NEGSTEER_PLOTS` | `bin/negsteer_plots.py` | Cohort-level negsteer plots. |
| `negsteer_per_seed_dispersion_overview.png`, `negsteer_per_seed_dispersion_grid.png`, `negsteer_per_design_stage_trajectories.png`, `negsteer_weighted_vs_true_jaccard.png`, `negsteer_composite_vs_confidence.png` (5 files) | `NEGSTEER_WITHIN_SEQUENCE_PLOTS` | `bin/negsteer_within_sequence_plots.py` | Per-seed (within-sequence) negsteer plots. |
| `orthogonal_af3_vs_boltz.png`, `orthogonal_filter_cascade.png`, `orthogonal_metrics_vs_composite.png`, `orthogonal_combined_cohort_summary.png` (4 files) | `ORTHOG_PLOTS` | `bin/orthogonal_metrics_plots.py` | Orthogonal-metrics cohort plots. |

PNGs are unsuitable for line-by-line characterization assertions. Pin them by file existence + size band, or render a deterministic test cohort and pin to a perceptual hash, but not pixel-equal.

---

## What's missing

Compared with the full pipeline as documented in `notes/inventory/01_module_map.md` and `04_functional_categorization.md`, the following kinds of outputs are not represented in this set. Decide before writing characterization tests whether each gap is acceptable, or whether a representative file should be added.

### Whole stages absent

- **Branch A — HADDOCK docking** (`modules/haddock.nf`, `bin/haddock3_*.py`, `bin/extract_hotspots.py`, `bin/build_contigs.py`, `bin/haddock_utils.py`). Nothing under a `haddock/` directory. AIR restraints, CAPRI/clustfcc tables, best-cluster/best-model PDB, hotspot lists, contigs reconstruction, and the four HADDOCK diagnostic plots from `bin/haddock3_plots.py` are all absent. Acceptable if the characterization tests intentionally only target Branch B; otherwise this is a gap.
- **Negative-steering controls** (`modules/negsteer_controls.nf`, `bin/build_control_sequences.py`, `bin/derive_input_design_region.py`, `DERIVE_INPUT_INDICES`, `NEGSTEER_CONTROLS`). The two synthetic non-design sequences (`scrambled`, `polyA`) and their per-sequence workdirs and output rows are not represented. The cohort summary in this set has 4 rows but the controls' rows are not visible; in a full run they would appear in `cross_sequence_summary.csv` as `mpnn_sequence=control_scrambled` / `control_polyA`.
- **Multi-cycle negative steering**. By design — the current Nextflow integration only runs `cycle_0`. The `pathways.json` here reflects this. Add multi-cycle examples only if/when the dormant `boltz2_iterate_steering.py` orchestration is wired in.

### Domain artefacts absent

- **PDBs**: no actual structure files anywhere in this set. RFDiffusion `design_*.pdb`, the split `passing/*.pdb`, the per-prediction `prediction.pdb` files inside `negative_steering/runs/<seq_name>/cycle_0/{initial_prediction,steered/,reverted/}/...`, the canonical PDBs referenced from `representative_canonical_pdb` in the cohort table, the AF3 per-sample CIFs, and the FastRelax-relaxed PDBs. Pinning these in characterization tests is expensive (file size + nondeterminism in random seeds); the inventory's recommendation is to pin metric outputs derived from them, not the structures themselves.
- **Boltz I/O**: no `input.yaml` files, no per-design `effector_template.cif`, no a3m/fasta inputs, no Boltz output CIFs.
- **ProteinMPNN raw outputs**: no `<pdb>_mpnn/` raw output directories, no fixed-positions JSONLs, no `selection_report.txt`.

### Negative-steering populated-reversion case

The reversion-related files in this example set are all empty (`contaminated.json` has `n_contaminated=0`, `reversion_plan.json` has zero entries, both `reversion_results*.json` are `{}`). A representative MPNN sequence whose cold-start *did* produce contaminating mutations would exercise:

- `contaminated.json` with populated `contaminated[]` entries.
- `reversion_plan.json` with `entries[]` carrying the per-design reversion YAMLs.
- `reversion_results_per_seed.json` with full per-seed reverted records.
- `aggregated_results.csv` rows where the `reverted_*` half-schema is populated.
- A non-trivial `reversion_verdict` distribution in `raw_per_seed_results.csv`.
- A `passing_summary.csv` with at least one populated row.

The current `passing_summary.csv` is empty-body, which is itself notable. A second per-sequence example with a populated reversion is recommended for the test set — without one, the reversion code path is untested by characterization.

### Per-survivor orthogonal-metrics structure

`orthogonal_metrics/` here holds the cohort-merged summaries. In a real run, each per-stream module also publishes per-survivor subdirectories (`orthogonal_metrics/af3/<seq_name>/`, `orthogonal_metrics/biophysical/<seq_name>/`, `orthogonal_metrics/rosetta/<seq_name>/`) with that survivor's raw outputs. None of those subtrees are present here.

### Other intermediates not represented

- `NEGSTEER_DERIVE_INDICES` per-design output (the `design_region.json` / `true_interface.json` files emitted by `bin/derive_design_region.py` and `bin/derive_true_interface.py`).
- `WRITE_DUMMY_MAPPING` JSON (trivial; reproducible).
- The `top_fastas/` and `af2_fastas/` directories in full (only one sample FASTA from `af2_fastas/` is included).

---

## Items flagged for confirmation

Three items above I could not resolve unambiguously from inventory + glossary alone:

1. **Top-level `rosetta_summary.csv`** — placement is anomalous; appears to be a hand-picked failure-case example normally found at `orthogonal_metrics/rosetta/<seq_name>/rosetta_summary.csv`. Confirm that's the intent (vs. it being misplaced).
2. **Empty-body `passing_summary.csv`** in `negative_steering/` — header-only file. Confirm whether this is intentional (an example of the empty case) or whether the intended file from a passing sequence got missed.
3. **Missing `interface_plddt` in `biophysical_summary.csv`** — `notes/inventory/01_module_map.md` lists `interface_plddt` as one of the three biophysical metrics, but the example file only has `bsa` and `interface_hbonds`. Confirm whether `interface_plddt` is computed and merged elsewhere, or whether the inventory's description is out of date.
