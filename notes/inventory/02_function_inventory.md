# 02 — Function Inventory

Per-file enumeration of top-level functions, classes, processes and workflows. Function summaries are derived from signatures, surrounding code, and docstrings. Module-level constants of structural significance are noted parenthetically. Duplication is called out inline; deeper findings live in `05_findings.md`.

The two largest files (`boltz2_iterate_steering.py` 5854 LOC, `boltz2_negative_steering.py` 3416 LOC) were sampled at function granularity. Bodies of long `cmd_*` subcommands were not read line-by-line; their summaries reflect docstrings and surrounding comments.

---

## Python — bin/

### bin/boltz2_iterate_steering.py (5854 LOC)

Multi-cycle negsteer orchestrator. Single argparse dispatch over 14 subcommands. Imports `boltz2_negative_steering` lazily (`_lazy_bns()`) and per-function (`from boltz2_negative_steering import get_chain_sequence`, `from boltz2_negative_steering import run_boltz`).

**Pathway helpers**
- `_lazy_bns()` — lazy-load `boltz2_negative_steering` module to avoid import-time cycle (the latter pulls numpy + Bio).
- `make_pathway_label(parent_label, cycle, design_idx)` → str — encode "cycle-N design-M" as a dot-separated parent.cNdMM label.
- `parse_pathway_label(label)` → list[(int, int)] — inverse.
- `pathway_workdir(experiment_root, label)` → Path — resolve a pathway label to its on-disk workdir.
- `parent_pathway_label(label)` → Optional[str] — strip leaf segment.
- `read_upstream_chain(experiment_root, leaf_label)` → list[Dict] — walk leaf → cycle 0 collecting per-cycle metadata.

**Mutation handling**
- `parse_mutations_tsv(path)` → list[int] — parse positions only.
- `read_mutations_tsv_full(path)` → list[(pos, wt_aa, mut_aa)] — full mutation tuples.
- `format_mutation_strings(...)` — build human-readable mutation labels.
- `read_cumulative_mutations(...)` — collect mutations across all ancestors.
- `accumulated_mutated_positions(chain)` → set — flatten upstream chain's mutated positions.

**Pose-distance maximin selection** (used by iterate-plan to choose which cycle-N+1 candidates to pursue)
- `compute_pose_distance(candidate_pdb, reference_pdb, ...)` — receptor-aligned RMSD.
- `maximin_select(candidates, k)` — greedy farthest-from-anything-seen selection.
- `append_cycle_statistics(experiment_root, ...)` — log per-cycle summary stats.

**Subcommand entry-points** (each ~50–500 LOC)
- `cmd_iterate_plan(args)` — pick cycle-N+1 candidates from cycle-N results.
- `cmd_compute_distances(args)` — populate min-dist columns on candidate pool.
- `cmd_iterate_collect(args)` — gather per-cycle prediction results.
- `cmd_iterate_collect_prefilter(args)` — apply eligibility predicate before reversion.
- `cmd_build_contaminated(args)` — identify mutations contacting effector inside protected set.
- `cmd_iterate_collect_finalize(args)` — produce finalized per-cycle CSV after reversion.
- `cmd_kickoff_distances/prefilter/finalize(args)` — submit downstream jobs.
- `cmd_kickoff(args)` — top-level cycle-kickoff for cycle 1+.
- `cmd_plan_reversions(args)` — emit reversion YAMLs per design.
- `cmd_predict_reversion_one(args)` — run Boltz on one reverted design.
- `cmd_harvest_reversions(args)` — collect reversion outputs.
- `cmd_aggregate(args)` — build flat per-prediction CSV.
- `cmd_aggregate_per_sequence(args)` — group flat rows by sequence.
- `cmd_compute_final_metrics(args)` — derived metrics + tier classification.

**Reversion-confidence forwarding**
- `_REVERTED_CONFIDENCE_FIELDS` (tuple of ~30 column names) — propagation list reversion.py → aggregator.
- `_copy_reverted_confidence(dst, rev)` — apply.

**Aggregation helpers** (block-B/C metric reduction)
- `_AGG_RENAME`, `_AGG_DROP`, `_AGG_CONTINUOUS_METRICS`, `_AGG_BINARY_METRICS`, `_AGG_INTEGER_METRICS`, `_AGG_POSITION_SET_METRICS` — schema dicts/tuples.
- `_translate_aggregate_row(r)` — column rename + drop.
- `_populate_reverted_mutations(rows, ...)` — annotate rows with reverted-mutation strings.
- `_row_is_clean_steered(r)`, `_row_is_pose_holds(r)` — eligibility predicates.
- `_ranking_metric(r, name)`, `_ranking_composite(r)` — rank-key computation.
- `_compute_unified_ranks(rows)` — final rank assignment.
- `_agg_continuous`, `_agg_majority_binary`, `_agg_majority_positions` — aggregators.
- `_parse_position_csv(s)` — parse "/A:5,7,18" position-set strings.
- `_aggregate_per_sequence(...)` — group + reduce per sequence.
- `_classify_aggregated_verdict(...)` — derive verdict from aggregated reverted-* columns.
- `_per_seed_verdict_breakdown(...)` — per-seed verdict counts.
- `_aggregate_csv_fieldnames(rows)` — final fieldname order.

**Final-metrics utilities**
- `_count_ca_per_chain(pdb_path)` — Cα count per chain.
- `_float_or_none(x)`, `_is_nan(x)` — coercion helpers.
- `_locate_boltz_sidecar_pdb(prediction_pdb)` — find the canonical-pdb sidecar to a prediction.
- `build_parser()` → argparse.ArgumentParser — all 14 subcommands.
- `main()` — dispatch.

### bin/boltz2_negative_steering.py (3416 LOC)

Single-cycle negsteer core; ALSO the central library that every other negsteer Python file imports from (`get_chain_sequence`, `run_boltz`, `find_contact_residues_heavy`, `jaccard`, `binding_rmsds`, `write_boltz_yaml`, `compute_effector_interface_residues`).

**Module-level constants**: `RECEPTOR_INTACT_CUTOFF=5.0`, `EFFECTOR_INTACT_CUTOFF=5.0`, `THREE_TO_ONE`, `_THREE_TO_ONE_RMSD`, `STEERING_SETS`, `CONSERVATIVE_SUBSTITUTIONS`, `PRED_REC_CHAIN="A"`, `PRED_EFF_CHAIN="B"`.

**Cα geometry / Kabsch alignment**
- `_read_ca_seq_by_chain(pdb_path)` — parse Cα sequences from a PDB.
- `_get_pairwise_aligner()` (singleton via `_PAIRWISE_ALIGNER`) — Bio.Align aligner.
- `_seqalign_pair_indices(seq_pred, seq_ref)` — sequence-alignment paired indices.
- `_kabsch_align(P, Q)` — classic Kabsch rotation matrix + translation.
- `compute_binding_rmsds(pred_pdb, design_pdb, ...)` — receptor-aligned effector RMSD vs design + truth.
- `extract_sequences_gemmi(path)` / `extract_sequences_biopython(path)` / `extract_sequences(path)` — three sequence-extraction routes (gemmi preferred, BioPython fallback, dispatcher).
- `get_chain_sequence(pdb_path, chain_id)` — single-chain sequence extraction. **Imported by 4 other bin/ scripts.**

**Boltz I/O**
- `write_single_seq_a3m(out_path, header, seq)` — write a single-sequence A3M.
- `extract_effector_template_cif(...)` — write an effector-only template CIF.
- `write_boltz_yaml(...)` — Boltz prediction YAML composer. **Imported by reversion.py.**
- `run_boltz(...)` — `boltz predict` subprocess wrapper. **Imported by boltz2_iterate_steering.py.**

**Residue / atom data classes**
- `class CAEntry` (dataclass) — Cα atom record.
- `read_ca_atoms(pdb_path)` → list[CAEntry].
- `_read_ca_seq_from_chain(pdb_path, chain_id)` — Cα 1-letter sequence.
- `class ResidueAtoms` (dataclass) — heavy-atom record.
- `read_residue_heavy_atoms(pdb_path)` → list[ResidueAtoms].

**Contact / interface analysis** — the heart of the steering logic
- `find_contact_residues_heavy(pdb_path, ...)` → set[int]. **Imported by 3 other bin/ files.**
- `_load_true_interface_indices(...)` — load 0-based interface indices from file or derive.
- `_load_design_region_indices(...)` — load 1-based design-region indices.
- `_second_shell_fallback(...)` — expand contact set by shell-2.
- `is_surface_exposed(...)` — SASA-based surface filter.
- `jaccard(a, b)` → float. **Imported by reversion.py.**
- `binding_rmsds(...)` — wrap compute_binding_rmsds with multi-chain logic. **Imported by reversion.py.**

**Steered design generation**
- `pick_mutant_residue(...)` — choose a substitution from a steering set.
- `make_steered_sequence(...)` — apply mutations to a sequence.

**CSV writers**
- `_write_initial_only_csv(...)` — when steering skipped (clean cold-start), emit single-row CSV.
- `_write_initial_multiseed_csv(...)` — when wrong-interface present but steering skipped, emit multi-seed cold-start CSV.

**Subcommand entry-points**
- `cmd_plan(args)` — full plan-phase workflow: cold-start prediction → wrong-interface detection → emit N steered design YAMLs.
- `cmd_predict_one(args)` — run Boltz on one steered design (one SLURM array task).
- `cmd_collect(args)` — rank predictions, emit `steered_results.csv`.
- `build_parser()`, `main()` — argparse dispatch.

### bin/build_contigs.py (237 LOC)

Used by `BUILD_CONTIGS` in HADDOCK module to convert hotspot residues + chain lengths into RFDiffusion contigs.

- `THREE_TO_ONE` — AA dict (also defined in `boltz2_negative_steering.py`, `extract_hotspots.py`, and elsewhere — see `05_findings.md`).
- `parse_args()`, `main()`.
- `load_mapping(path)` — read sequence-to-PDB mapping JSON.
- `extract_chain(pdb_path, chain_id)` — chain residue list.
- `resolve_chains(pdb_path, rec_chain, eff_chain)` — disambiguate chains.
- `parse_contig_receptor_refs(contigs, rec_chain)` — parse user-provided contigs.
- `compute_offset(contigs, rec_chain, rec_resnums, rec_map)` — figure out residue-numbering offset.
- `update_contigs(...)` — rewrite contigs with denovo lengths derived from hotspot count.

### bin/build_control_sequences.py (340 LOC)

Generates scrambled / polyA negative controls. Run by `NEGSTEER_CONTROLS`.

- `_extract_chain_sequence(pdb_path, chain_id)` → str. **Functionally equivalent to `boltz2_negative_steering.get_chain_sequence` and `extract_survivor_manifest._extract_chain_seq`.**
- `_load_design_region_indices_1b(design_region_file)` → list[int].
- `_load_effector_sequence(effector_fasta)` → str.
- `_build_scrambled_receptor(...)` — permute design-region residues under a deterministic seed.
- `_build_polyA_receptor(...)` — substitute alanines (preserving glycine in some cases).
- `_write_single_record_fasta(...)` — single-record FASTA writer.
- `_default_seed_from_pdb(source_pdb)` → int — deterministic seed from PDB content hash.
- `main()`.

Imports `get_chain_sequence` from `boltz2_negative_steering` for receptor-sequence extraction (so it has TWO ways to do the same thing, see `05_findings.md`).

### bin/collect_haddock3_dock.py (291 LOC)

Used by `HADDOCK3_DOCK`. Imports many helpers from `haddock_utils`.

- `MIN_CLUSTER_SIZE = 4` (also defined in haddock.nf, haddock3_plots.py, haddock_utils.py — flagged in pipeline_notes1).
- `parse_args()`, `main()`.
- `find_file(run_dir, analysis_pattern, step_pattern, filename)` — locate HADDOCK output.
- `select_all_qualifying_clusters(clusters)` — clusters above min size.
- `select_best_cluster(clusters)` — best by mean score.
- `extract_model_by_name(run_dir, target_model_name, output_name)` — extract the chosen PDB.
- `_fail(report, message)`, `_validate_score(score, report)` — error helpers.

### bin/compute_interface_metrics.py (739 LOC)

Per-row interface metrics for `NEGSTEER_INTERFACE_METRICS`. Imports `compute_ipsae` from `compute_metrics` and DockQ Python API.

- `NEW_COLUMNS` (list of column-name strings).
- `_find_sidecar(pdb_path, prefix, ext)` — locate companion files.
- `_find_ground_truth(seq_name, workdirs_glob)` — resolve plan.json → ground-truth PDB.
- `_run_dockq(...)` — DockQ subprocess wrapper.
- `_load_pae(...)` — load PAE matrix from Boltz output.
- `_compute_intact_core(...)` — pLDDT-trimmed RMSD-filter intact-core.
- `_WJ_MU=4.0`, `_WJ_TWO_SIGMA_SQ=2.25` — weighted-Jaccard parameters. **Duplicated in `compute_metrics.py` with the same constants.**
- `_gaussian_contact_weight(distance_angstroms)` — Gaussian weight at distance.
- `_collect_chain_cb_positions(...)` — Cβ collection.
- `_build_weighted_pair_map(...)` — pair-map for weighted-Jaccard.
- `_compute_weighted_jaccard(...)` — distance-weighted Jaccard.
- `process_row(...)` — per-row processing (one survivor's metrics).
- `main()` → int.

### bin/compute_metrics.py (1657 LOC)

Boltz output metrics. Run via subprocess by `boltz2_iterate_steering.py` (`compute-final-metrics`, `iterate-collect-finalize`) and `reversion.py` (`harvest_reversion_results`); imported directly by `compute_interface_metrics.py`, `derive_input_design_region.py`, and `boltz2_negative_steering.py`.

**iPSAE / actifPTM / iPAE / PAE-pass-frac**
- `_ipsae_d0(n)`, `compute_ipsae_one_direction(pae_ij, cutoff)`.
- `_pae_to_ptm_score(pae_ij, n_interface)`.
- `compute_actifptm(pae_matrix, chain_lengths, cutoff)`.
- `compute_ipsae(pae_matrix, chain_lengths, cutoff)` — exported.
- `compute_ipae(pae_matrix, chain_lengths)`.
- `compute_pae_pass_frac(pae_matrix, chain_lengths, cutoff)`.
- `_pae_derived_metrics(pae_matrix, chain_lengths, cutoff)` — bundle.

**Heavy-atom / contact analysis**
- `_parse_heavy_atoms_by_chain(pdb_path)` — load heavy atoms.
- `compute_effector_interface_residues(...)` — exported, called by `boltz2_negative_steering.py`.
- `compute_effector_interface_atom_mask(...)`.
- `compute_interface_contacts(...)`.
- `find_contact_residues_heavy(...)` — exported. **Distinct from `boltz2_negative_steering.find_contact_residues_heavy` despite the same name; both exist (see `05_findings.md`).**

**Mutation reliance / position parsing**
- `parse_mutated_positions(spec)`, `parse_position_list(spec)`, `load_positions_file(path)`.
- `classify_mutation_reliance(...)` — derives `mutation_reliance_flag` and `construct_reliance_flag`.
- `plddt_from_pdb(pdb_path)`.
- `_finalize_entry(entry, pae_matrix, chain_lengths, pae_cutoff)`.

**Weighted-Jaccard / intact-core**
- `_WJ_MU=4.0`, `_WJ_TWO_SIGMA_SQ=2.25` (same constants as compute_interface_metrics.py).
- `_gaussian_contact_weight(distance)`, `_collect_chain_cb_positions(structure, chain_id)`, `_build_weighted_pair_map(...)`.
- `compute_weighted_jaccard(model_pdb, native_pdb, ...)`.
- `compute_intact_core(model_pdb, native_pdb, plddt_threshold, ...)`.
- `compute_interface_plddt(model_pdb, ...)`.

**Boltz parser**
- `_dedup_paths(paths)`, `_is_reference(path)`.
- `parse_boltz2(...)` — read all relevant Boltz outputs into a single record.
- `PARSERS = {"boltz2": parse_boltz2}` — dispatch table.
- `CSV_FIELDS` — output CSV column order.
- `main()` — CLI entry-point.

### bin/contig_utils.py (326 LOC)

Shared contig parsing. Imported by `rfdiffusion_contigs.py` and `rfdiffusion_filter.py`.

- `get_chain_residue_range(pdb_path, chain_id)` — first/last resnum.
- `get_chain_residues_sorted(pdb_path, chain_id)` — sorted resnums.
- `parse_block_segments(block, rec_chain, pdb_path)` — parse one contigs block.
- `_resolve_denovo_lengths(seg_descs, pdb_total)` — resolve `N-N` ranges.
- `remap_segments_to_pdb(seg_descs, block_chain, pdb_path)` — fail loudly if a referenced chain has no atoms.
- `segments_to_string(seg_descs)` — round-trip.
- `resolve_contigs(raw_contigs, pdb_path)` — top-level resolution.
- `get_expected_chain_lengths(contigs, rec_chain, eff_chain)` — receptor / effector lengths.
- `parse_design_region(contigs, rec_chain)` — extract de novo region indices.

### bin/cross_sequence_summary.py (904 LOC)

Cohort aggregator run by `NEGSTEER_CROSS_SEQUENCE`. Imports `extract_passing.extract_row`.

- `_rewrite_workdir_path_to_published(raw_path, seq_name, published_runs_dir)` — rewrites Nextflow `work/<hash>/...` paths to `${params.outdir}/.../runs/<seq>/...` so downstream container tasks can bind the file.
- `_try_float(v)`, `_try_int(v)` — safe coercion.
- `_tier_for_row(row)` — tier-A/B/C/none classification.
- `_TIER_ORDER = {"A":0, "B":1, "C":2, "none":3}`.
- `_composite_score(row)` — `true_jaccard − 0.05·ra_eff`.
- `_within_sequence_rank(row)` — within-MPNN-sequence rank.
- `_SEQ_NAME_RE` — sequence-name validation.
- `_infer_seq_name(path, mode)`.
- `_parse_passing_summary_arg(arg, seq_name_mode)`.
- `_pick_representative(rows)` — pick best tier-A/B/C row.
- `_pick_representative_from_aggregated(...)` — Bug-3 fallback for tier-none sequences with empty passing_summary; reads aggregated_results.csv directly and uses `extract_passing.extract_row` to shape the output.
- `aggregate(...)` — main aggregation.
- `main()`.

### bin/derive_design_region.py (379 LOC)

Used by `NEGSTEER_DERIVE_INDICES`.

- `normalise_design_id(raw)` → str.
- `find_design_entry(metrics, design_id)` — locate the design block in `rfdiffusion_metrics.json`.
- `derive_positional_indices(entry)` — convert resnum ranges to 1-based positional indices.
- `cross_check_with_plan(plan_path, design_indices_0based, ...)` — sanity-check vs plan.
- `format_output(indices_1based, metrics_path, design_id, ranges)` — CLI output.
- `main()`.

### bin/derive_input_design_region.py (398 LOC)

Used by `DERIVE_INPUT_INDICES` (controls flow). Computes design-region + true-interface for the input complex (pre-RFDiffusion). Imports `find_contact_residues_heavy` from `compute_metrics` and `get_chain_sequence` from `boltz2_negative_steering`.

- `_parse_contigs(...)` — local contigs parser. **Conceptually duplicates `contig_utils.parse_block_segments`.**
- `_receptor_positions_by_contigs(...)` — derive 1-based positions.
- `_design_region_positions_1b(...)` — emit design-region indices.
- `_compute_true_interface_0b(...)` — derive true-interface indices.
- `_extract_receptor_sequence(...)`.
- `_write_design_region(...)`, `_write_true_interface(...)`.
- `main()`.

### bin/derive_true_interface.py (326 LOC)

Used by `NEGSTEER_DERIVE_INDICES`. **Sibling to `derive_design_region.py` with very similar shape — same `normalise_design_id`, `find_design_entry`, `derive_positional_indices`, `format_output` quartet.** Differs only in which subset of indices is emitted.

- `normalise_design_id(raw)` — same as in derive_design_region.py.
- `find_design_entry(metrics, design_id)` — same.
- `derive_positional_indices(entry)` — same name, different body (emits 0-based contact residues).
- `cross_check_design_region(design_region_path, interface_indices_0based)`.
- `format_output(indices_0based, metrics_path, design_id, n_contact_resnums)`.
- `main()`.

### bin/extract_hotspots.py (207 LOC)

Used by `EXTRACT_HOTSPOTS`. Imports `extract_heavy_atoms`, `find_interface_residues` from `haddock_utils`.

- `THREE_TO_ONE` — yet another copy.
- `parse_args()`, `main()`.
- `_read_chain_sequences(pdb_path)`.
- `_best_identity(query, reference)` — sliding-window identity for chain disambiguation.
- `resolve_chains(pdb_path, rec_chain, eff_chain, rec_ref, eff_ref)` — sequence-identity chain disambiguation (the H10 fix from pipeline_notes1).

### bin/extract_passing.py (678 LOC)

Used by `negative_steering_run_one.sh`; also imported by `cross_sequence_summary.py`.

- `OUTPUT_FIELDS` — explicit ordered field list.
- `_get(row, key, fallback_key)` — column lookup with fallback.
- `_compute_confidence_flag(...)` — combined confidence pass/fail.
- `_summarize_mutations_across_seeds(...)`.
- `_muts_aa_for_positions(...)` — convert position list to mutation strings.
- `extract_row(r, per_seed_rows)` — translate one aggregated row → passing-summary row. **Stage selection rule (Bug 4 fix in pipeline_notes14): use reverted_* whenever reversion was attempted, steered_* only for no_reversion.**
- `main()` — CLI; filters input to PASSING_VERDICTS.

### bin/extract_survivor_manifest.py (174 LOC)

Used by `EXTRACT_SURVIVOR_MANIFEST`.

- `_extract_chain_seq(pdb_path, chain_id)` → Optional[str]. **Third independent implementation of single-chain sequence extraction** (cf. `boltz2_negative_steering.get_chain_sequence`, `build_control_sequences._extract_chain_sequence`).
- `_find_workdir(seq_name, workdirs_glob)`.
- `_resolve_plan_json(workdir)`.
- `main()`.

### bin/haddock_utils.py (283 LOC)

Shared helpers; imported by `collect_haddock3_dock.py`, `haddock3_plots.py`, `extract_hotspots.py`.

- `parse_capri_tsv(capri_file)` — CAPRI scoring table.
- `parse_clustfcc_tsv(clust_file, min_cluster_size=4)` — cluster table.
- `_first_col(row, candidates)`, `_first_float(row, candidates)`.
- `get_score(row)`, `get_model_name(row)`, `get_numeric_col(data, name_options)`, `cluster_mean_score(members)`, `model_stem(model_path)`.
- `iter_pdb_lines(pdb_path)` — yield PDB lines (uncompresses .pdb.gz transparently).
- `extract_heavy_atoms(pdb_path, chain_ids)` — yield heavy-atom records.
- `copy_pdb(src, dst)`.
- `contacted_residues(query_atoms, target_atoms, cutoff)` — pair-distance contact set.
- `find_interface_residues(rec_atoms, eff_atoms, cutoff)`.
- `collect_pdb_index(complex_dir, run_dir)` — locate per-cluster PDBs.

### bin/haddock3_plots.py (697 LOC)

Used by `HADDOCK3_PLOTS`. Imports several helpers from `haddock_utils`.

- `BSA_WARN_CUTOFF=1000.0`, `MIN_CLUSTER_SIZE=4`.
- `parse_args()`, `main()`.
- `make_empty_plot(message, path)`, `save_fallback_plots(message)` — fallback plot helpers used in every plot script in the repo (each has its own copy).
- `plot_score_vs_bsa(...)`, `plot_cluster_sizes(...)`.
- `_compute_contact_frequencies(...)`, `_parse_fixed_residues(...)`, `_design_ranges_from_fixed(...)`, `_build_heatmap_matrix(...)`, `_render_heatmap_panel(...)`, `_set_heatmap_yticks(...)`, `_set_heatmap_xticks(...)`.
- `plot_interface_heatmap(...)` — main interface contact heatmap.

### bin/haddock3_prepare.py (273 LOC)

Used by `HADDOCK3_PREPARE`. Generates AIR restraints.

- `parse_args()`, `main()`.
- `parse_contig_segments(contigs, rec_chain)` — yet another contig parser. **Conceptually duplicates `contig_utils.parse_block_segments` and `derive_input_design_region._parse_contigs`.**
- `find_denovo_residues(parsed_segments, receptor_pdb, rec_chain)`.
- `format_ranges(residue_list)`.
- `write_air_restraints(path, active_residues, rec_chain, eff_chain, eff_active_residues)`.
- `parse_effector_residue_spec(spec)`.

### bin/merge_orthogonal_metrics.py (204 LOC)

Used by `NEGSTEER_ORTHOGONAL_METRICS`. **Has a divergent copy at `tests/orthogonal_metrics/merge_orthogonal_metrics.py` — see `05_findings.md`.**

- `BIO_COLS`, `ROSETTA_COLS`, `SUMMARY_COLS`, `ALL_NEW` — column lists (AF3_COLS defined in same block).
- `_load_summaries(pattern, key_col)` → Dict[str, Dict].
- `_as_float(value)` → Optional[float].
- `_apply_filters(...)` — populate orthogonal_flags + passes_orthogonal_filters. **Production version: AF3 is gating. Test-tree version: AF3 is informational only.**
- `main()`.

### bin/mpnn_cluster_sequences.py (175 LOC)

Used by `MPNN_CLUSTER`.

- `parse_args()`, `main()`.
- `write_design_fasta(rows, fasta_path)`.
- `run_mmseqs_cluster(fasta_path, threshold, mmseqs_bin, work_dir)`.

### bin/mpnn_design_region_score.py (319 LOC)

Used by `MPNN_DESIGN_REGION_SCORE`. Imports `torch` + `protein_mpnn_utils` (only file in the codebase that does).

- `MPNN_DIR="/opt/ProteinMPNN"`.
- `ALPHABET="ACDEFGHIKLMNPQRSTVWYX"`, `AA_TO_IDX`.
- `parse_args()`, `main()`.
- `find_checkpoint()` → Path.
- `load_model(checkpoint_path, device)`.
- `parse_mpnn_fasta(fasta_path)`, `load_fixed_positions(jsonl_path)`.
- `score_sequences_for_design(pdb_path, fasta_path, fixed_dict, ...)`.

### bin/mpnn_plots.py (550 LOC)

Used by `MPNN_PLOTS`. Plot script with `make_empty_plot` / `save_fallback_plots` / `ALL_PLOT_FILES` idiom.

- `AA_ORDER`, `KD_HYDROPHOBICITY`, `CHARGE`, `COLOUR_*` constants. **`KD_HYDROPHOBICITY`, `CHARGE`, `COLOUR_DESIGN`, `COLOUR_NATIVE`, `COLOUR_REGION` also defined identically in `tests/proteinmpnn/test_mpnn_plots.py`.**
- `parse_args()`, `main()`.
- `make_empty_plot(message, path)`, `save_fallback_plots(message)`.
- `_design_label(design_idx)`.
- `plot_score_distribution(rows)`.
- `plot_sequence_diversity(cluster_csv)`.
- `_aa_composition_pct(seq_string)`.
- `plot_aa_composition(rows)`.
- `_mean_hydrophobicity(seq)`, `_net_charge(seq)`.
- `plot_physicochem(rows)`.

### bin/mpnn_select_top.py (72 LOC)

Used by `MPNN_SELECT_TOP`. Two functions:

- `parse_args()`, `main()`.

### bin/mpnn_sequence_qc.py (195 LOC)

Used by `SEQUENCE_QC`.

- `parse_args()`, `main()`.
- `has_poly_x(seq, max_run)` — flag designs with long poly-X runs.
- `unusual_aa_composition(seq, max_single_aa_frac=0.30)` — flag mono-AA-dominated designs.

### bin/negsteer_plots.py (1911 LOC) — production lift of test_negsteer_plots.py

7-plot cohort suite for negsteer. **Verbatim copy of the test script** modulo header and `_resolve_csv_path` fallback.

- Module-level constants: `COMPOSITE_RA_EFF_WEIGHT=0.05`, `MUT_IMPACT_MIN_N=3`, `COLOUR_TIER`, `COLOUR_DESIGN_REGION`, `_AGG_FINAL_FIELD_MAP`, `_VERDICTS_WITH_REVERSION`, `OUTCOME_COLOUR`, `OUTCOME_LABEL`, `OUTCOME_ORDER`, `STAGE_MARKER`, `STAGE_COLOUR`.
- `_final_prediction_value(agg_row, dest_field)`, `_project_agg_row(...)`, `_pick_best_agg_row(...)`.
- `load_unified_cohort(...)` — cohort CSV loader.
- `_try_float`, `_try_int`, `_parse_chimerax_positions`, `_composite_from_row`, `_row_type`, `_is_steered`, `_make_empty_plot`, `_short_name`.
- `plot_tier_landscape(rows, out_path)`.
- `classify_seed(row)` → (stage, outcome).
- `_representative_sg_for_outcomes(...)`.
- `load_seed_outcomes(...)`.
- `plot_seed_outcomes_heatmap(...)`, `plot_seed_outcomes_bars(...)`.
- `plot_ra_eff_vs_jaccard(rows, out_path)`.
- `plot_filter_cascade(rows, thresholds, out_path)`.
- `plot_controls_diagnostic(rows, thresholds, out_path)`.
- `plot_mutation_impact(...)`.
- `_resolve_csv_path(supplied)`, `_load_input_design_region(...)`, `_parse_args()`, `main()`.

### bin/negsteer_within_sequence_plots.py (987 LOC)

5-plot per-seed suite. **Verbatim copy of the test script.**

- Module-level: `COLOUR_TIER`, `COLOUR_VERDICT_HOLD`, `COMPLEX_PLDDT_MIN`, `PAE_PASS_FRAC_MIN`, `STAGE_MARKER`, `_VERDICTS_WITH_REVERSION`.
- `_try_float`, `_try_int`, `_short_name`, `_design_id`, `_seq_num`, `_make_empty_plot`, `_is_control`, `_grid_shape`.
- `classify_seed(row)` → (stage, outcome). **Independent reimplementation of the same-named function in `negsteer_plots.py` — distinct code path on `_is_steered`.**
- `_representative_sg(seq_dir, cs_lookup)`.
- `load_cohort_seeds(...)`.
- `plot_dispersion_overview(rep_seeds, tiers, out_path)`.
- `plot_dispersion_grid(rep_seeds, all_seeds, tiers, composites, out_path)`.
- `plot_stage_trajectories(all_seeds, tiers, cold_baseline, out_path)`.
- `plot_weighted_vs_true_jaccard(rep_seeds, tiers, out_path)`.
- `plot_composite_vs_confidence(rep_seeds, composites, tiers, out_path)`.
- `_parse_args()`, `main()`.

### bin/orthogonal_metrics_plots.py (1146 LOC)

4-plot cohort suite for orthogonal metrics. **Verbatim copy of the test script.**

- `DDG_BENNETT_REFERENCE=-30.0`, `COMPOSITE_RA_EFF_WEIGHT=0.05`, `COLOUR_TIER`, `_TIER_SORT_ORDER`.
- `_try_float`, `_try_int`, `_parse_flags`, `_passes`, `_make_empty_plot`, `_short_name`.
- `plot_af3_vs_boltz(rows, out_path)`.
- `FILTER_ORDER` (list), `ORTHOG_PRESENCE_COLS` (tuple).
- `_entered_orthogonal_stage(row)`.
- `plot_filter_cascade(rows, out_path)`.
- `METRIC_VS_COMPOSITE_PANELS` — panel definitions list.
- `plot_metrics_vs_composite(...)`.
- `COMBINED_SUMMARY_COLUMNS` (list), `COLOUR_NEUTRAL`.
- `_format_cell(v, fmt)`, `_representative_n_mutations(cross_row)`.
- `plot_combined_cohort_orthogonal_summary(...)`.
- `_load_rows(path, label)`, `_parse_args()`, `main()`.

### bin/parse_af3_output.py (320 LOC)

Used by `AF3_PARSE_OUTPUT`.

- `_chain_ca_coords(structure, chain_id)` — pair-by-position list of xyz (after Bug-6 fix in pipeline_notes14).
- `_receptor_aligned_effector_rmsd(...)` — Kabsch-aligned RMSD.
- `_read_confidence(cif_path)` → (Optional[float], Optional[str]).
- `main()`.
- `_write_single(path, fieldnames, row)`.

### bin/pipeline_correct_sequences.py (770 LOC)

Used by `SEQUENCE_CORRECTION`. Cross-pipeline file (in receptor-resurfacing tree, patched in negsteer sessions).

- `parse_contig_segments(contigs, receptor_chain="A")` — yet another contig parser. **Distinct from `contig_utils.parse_block_segments`, `derive_input_design_region._parse_contigs`, `haddock3_prepare.parse_contig_segments`.**
- `get_fixed_residue_set(segments, receptor_start_pdb=1)`.
- `get_native_anchor_regions(segments)`.
- `get_pdb_chain_residues(pdb_path)` — one of three implementations.
- `get_pdb_sequence(pdb_path, chain_id)` — fourth implementation of single-chain sequence extraction.
- `parse_mpnn_fasta(fasta_path)`, `extract_mpnn_score(header)`.
- `split_sequence(seq, effector_len)` — the historic chain-order bug location (notes 4).
- `align_to_native_by_anchors(mpnn_receptor, native_receptor, segments, receptor_start_pdb)` — anchor-based alignment back to native.
- `generate_fixed_positions(args)` — emit fixed-positions JSONL.
- `correct_sequences(args)` — main correction logic.
- `main()`.

### bin/reversion.py (1152 LOC)

Used by the orchestrator's harvest-reversions step. Imports `write_boltz_yaml`, `binding_rmsds`, `find_contact_residues_heavy`, `jaccard` from `boltz2_negative_steering`.

- `build_reverted_sequence(...)` — apply reversions to a steered sequence.
- `_read_receptor_fasta(fasta_path)` → str.
- `_read_ca_chain(pdb_path, chain)` — Bug-fix helper from pipeline_notes13 (length-only guard).
- `write_reversion_plan(...)` — emit reversion-plan JSON.
- `harvest_reversion_results(workdir, plan_meta, compute_metrics_script, ...)` — collect predictions, run compute_metrics.py via subprocess, populate per-seed records (incl. structural-jaccard fields).
- `classify_reversion_verdict(...)` — derive `pose_holds` / `pose_collapses` / `new_contamination` / `no_data`.

### bin/rfdiffusion_contigs.py (42 LOC)

Used by `RESOLVE_CONTIGS`. Thin CLI wrapper around `contig_utils.resolve_contigs`.

- `parse_args()`, `main()`.

### bin/rfdiffusion_filter.py (926 LOC)

Used by `RFDIFFUSION_FILTER`. Imports `resolve_contigs`, `get_expected_chain_lengths`, `parse_design_region` from `contig_utils`.

- `parse_args()`, `main()`.
- `read_ca_atoms(pdb_path, chain=None)` — local Cα reader. **Distinct implementation from `boltz2_negative_steering.read_ca_atoms`.**
- `ca_coords_array(atoms)`, `count_chains(pdb_path)`.
- `_find_break_index(coords_or_atoms)`, `split_at_chain_break(ca_atoms)`.
- `identify_segments(seg1_atoms, seg2_atoms, expected_rec_len, expected_eff_len)`.
- `write_split_pdb(input_pdb_path, output_path, expected_rec_len, expected_eff_len)` — hardcodes A=receptor, B=effector.
- `_parse_receptor_segments(contigs, rec_chain)`.
- `build_receptor_resnum_map(rec_atoms, contigs, rec_chain)`.
- `find_contacts(rec_atoms, eff_atoms, cutoff)` — local contact finder.
- `calc_scaffold_and_region_metrics(...)` — per-design metric bundle.

### bin/rfdiffusion_plots.py (1248 LOC)

Used by `RFDIFFUSION_PLOTS`.

- `COLOUR_COVERAGE`, `ALL_PLOT_FILES`, `LABEL_DENDRO_GAP_IN`, `DENDRO_W_INCHES`, `RIGHT_MARGIN_IN`, `TOP_MARGIN_IN`, `BOTTOM_MARGIN_IN` (layout tunables).
- `parse_args()`, `main()`.
- `make_empty_plot(message, path)`, `save_fallback_plots(message)`.
- `_design_label(d)`.
- `_get_per_design_residues(d, global_design_residues)`, `_split_design_into_regions(...)`, `_get_dendrogram_order(designs)`, `_get_per_region_lengths(d)`, `_label_with_region_len(...)`, `_per_region_clustering_data(...)`, `_n_regions(designs)`.
- `_build_contact_map_data(...)`, `_build_effector_freq(...)`, `_build_design_labels(...)`.
- `_draw_contact_heatmap_content(...)`, `_set_contact_xticks(...)`, `_draw_effector_freq_bar(...)`, `_contact_legend_handles()`.
- `_draw_clustering_heatmap_content(...)`, `_measure_label_width_inches(labels, fontsize, dpi)`.
- `plot_contact_map(...)`, `plot_design_clustering(...)`, `plot_specificity_coverage(...)`, `plot_design_lengths(...)`, `plot_com_displacement(...)`.

### bin/rosetta_filter_collect.py (227 LOC)

Used by `ROSETTA_FILTER`.

- `parse_args()`, `main()`.
- `parse_rosetta_scorefile(path)`.
- `safe_float(val, default=float("nan"))`.

### bin/rosetta_filter_plots.py (206 LOC)

Used by `ROSETTA_FILTER_PLOTS`.

- `COLOUR_WARN`, `COLOUR_FAIL`, `ALL_PLOT_FILES`.
- `parse_args()`, `main()`.
- `make_empty_plot(message, path)`, `save_fallback_plots(message)` — yet another copy.
- `plot_sc_histogram(designs, threshold)`.

### bin/run_biophysical_metrics.py (333 LOC)

Used by `NEGSTEER_BIOPHYSICAL_METRICS`.

- `_interface_residues(...)` — interface residue selection.
- `_compute_bsa(pdb_path, ...)` — BSA via FreeSASA.
- `_compute_interface_plddt(...)` — interface pLDDT from cached B-factors.
- `_compute_hbonds(...)` — H-bonds via MDAnalysis.
- `main()`.

### bin/run_rosetta_metrics.py (411 LOC)

Used by `NEGSTEER_ROSETTA_METRICS`.

- `ROSETTA_FLAVOURS = ("mpi", "default", "static", "linuxgccrelease")`.
- `_find_rosetta_binary(stem)` → Optional[str].
- `_parse_ia_scorefile(...)` — InterfaceAnalyzer scorefile parser.
- `_run_fast_relax(...)` — FastRelax stage with stdout/stderr capture (Bug 1 fix).
- `_run_interface_analyzer(...)` — IA stage.
- `_emit_row(...)`, `_emit_failure_row(...)`.
- `main()`.

### bin/sequence_registry.py (175 LOC) — DEAD CODE

Never imported, never invoked from .nf or .sh. Already flagged as dead in pipeline notes.

- `_seq_hash(sequence)` → str.
- `load_registry(experiment_root)` → Dict.
- `save_registry(experiment_root, registry)` → Path.
- `lookup(registry, sequence)` → Optional[Dict].
- `register(...)`.
- `count_existing_seeds(registry, sequence)` → int.
- `get_existing_dirs(registry, sequence)` → List[str].
- `deduplicate_sequences(...)`.

---

## Python — tests/

### tests/orthogonal_metrics/merge_orthogonal_metrics.py (216 LOC)

**Divergent copy of `bin/merge_orthogonal_metrics.py`.** Same function signatures, divergent body for `_apply_filters` (AF3 demoted to flag-only).

### tests/orthogonal_metrics/test_orthogonal_metrics_plots.py (1170 LOC)

Source of truth for `bin/orthogonal_metrics_plots.py`. Same function set; only diff is the `_resolve_csv_path` fallback (test-path defaults). Includes a dual-csv arg parser.

### tests/proteinmpnn/test_mpnn_plots.py (625 LOC)

Iteration script for `bin/mpnn_plots.py`. Same plot functions (`plot_score_distribution`, `plot_aa_composition`, `plot_sequence_diversity`, `plot_physicochem`) and same constants (`AA_ORDER`, `KD_HYDROPHOBICITY`, `CHARGE`, `COLOUR_*`); test version takes `out_path` per-plot rather than relying on cwd.

- `_design_label`, `plot_score_distribution`, `_aa_composition_pct`, `plot_aa_composition`, `plot_sequence_diversity`, `_mean_hydrophobicity`, `_net_charge`, `plot_physicochem`, `_make_empty_plot`, `_resolve_metadata_path`, `_parse_args`, `main`.

### tests/rfdiffusion/test_rfdiffusion_plots.py (1163 LOC)

Iteration script. Imports the production module to render the three "good" plots (specificity_coverage, design_lengths, com_displacement) and reimplements the two iteration-target plots locally.

- `_import_production_plots(prod_path)` — runtime import of `rfdiffusion_plots.py`.
- All same `_design_label`, `_get_per_design_residues`, `_get_per_region_lengths`, `_label_with_region_len`, `_get_dendrogram_order`, `_per_region_clustering_data`, `_n_regions`, `_build_contact_map_data`, `_build_effector_freq`, `_build_design_labels`, `_draw_contact_heatmap_content`, `_set_contact_xticks`, `_draw_effector_freq_bar`, `_contact_legend_handles`, `_draw_clustering_heatmap_content`, `_measure_label_width_inches`, `plot_contact_map`, `plot_design_clustering`, `plot_design_lengths`, `plot_com_displacement` functions as production. **Roughly 80% identical code that was lifted.**
- `_resolve_metrics_path(supplied)`, `_parse_args`, `main`.

### tests/rosetta_filtering/test_rosetta_filtering_plots.py (224 LOC)

- `COLOUR_WARN`, `COLOUR_FAIL`.
- `plot_sc_histogram(designs, threshold, out_path)`, `_make_empty_plot`, `_resolve_metrics_path`, `_parse_args`, `main`.

### tests/negative_steering/test_negsteer_plots.py (1916 LOC)

**Source of truth for `bin/negsteer_plots.py`.** Same constants, same functions, same logic; production differs only in header and `_resolve_csv_path` fallback.

### tests/negative_steering/test_negsteer_within_sequence_plots.py (981 LOC)

**Source of truth for `bin/negsteer_within_sequence_plots.py`.** Identical except for header.

---

## Nextflow — modules/ and tests/

### modules/preprocessing.nf

- `process EXTRACT_SEQUENCES` — extract receptor/effector sequences + lengths from input PDB; emit `sequences.json`.
- `process RESOLVE_CONTIGS` — invoke `bin/rfdiffusion_contigs.py` to resolve `N-N` contig ranges against the PDB.
- `process WRITE_DUMMY_MAPPING` — emit a 2-byte `{}` JSON for downstream consumers that need a mapping file.

### modules/haddock.nf

- `process HADDOCK3_PREPARE` — run `bin/haddock3_prepare.py` to write AIR restraints + cfg.
- `process HADDOCK3_DOCK` — run HADDOCK3 + `bin/collect_haddock3_dock.py` to pick the best cluster's best model.
- `process HADDOCK3_PLOTS` — `bin/haddock3_plots.py`.
- `process EXTRACT_HOTSPOTS` — `bin/extract_hotspots.py`.
- `process BUILD_CONTIGS` — `bin/build_contigs.py`.

### modules/rfdiffusion.nf

- `process RFDIFFUSION` — invoke RFDiffusion (Singularity) to generate N backbones.
- `process RFDIFFUSION_FILTER` — `bin/rfdiffusion_filter.py`.
- `process RFDIFFUSION_PLOTS` — `bin/rfdiffusion_plots.py`.

### modules/rosetta_filtering.nf

- `process ROSETTA_SC` — per-design Rosetta InterfaceAnalyzer (parallelised).
- `process ROSETTA_FILTER` — `bin/rosetta_filter_collect.py`.
- `process ROSETTA_FILTER_PLOTS` — `bin/rosetta_filter_plots.py`.

### modules/proteinmpnn.nf

- `process MPNN_FIXED_POSITIONS` — `bin/pipeline_correct_sequences.py` (subcommand path) to emit fixed-positions JSONL.
- `process PROTEINMPNN` — invoke ProteinMPNN.
- `process SEQUENCE_CORRECTION` — `bin/pipeline_correct_sequences.py`.
- `process SEQUENCE_QC` — `bin/mpnn_sequence_qc.py`.
- `process MPNN_DESIGN_REGION_SCORE` — `bin/mpnn_design_region_score.py`.
- `process MPNN_CLUSTER` — `bin/mpnn_cluster_sequences.py`.
- `process MPNN_SELECT_TOP` — `bin/mpnn_select_top.py`.
- `process MPNN_PLOTS` — `bin/mpnn_plots.py`.

### modules/negative_steering.nf

- `process NEGSTEER_DERIVE_INDICES` — `bin/derive_design_region.py` + `bin/derive_true_interface.py` per design.
- `process NEGSTEER_RUN_ONE` — `bash bin/negative_steering_run_one.sh` per MPNN sequence (collapsed-chain GPU job).
- `process NEGSTEER_CROSS_SEQUENCE` — `bin/cross_sequence_summary.py`.
- `process NEGSTEER_PLOTS` — `bin/negsteer_plots.py` (passed as `path script_input` after Bug-7 cache fix).
- `process NEGSTEER_WITHIN_SEQUENCE_PLOTS` — `bin/negsteer_within_sequence_plots.py`.

### modules/negsteer_controls.nf

- `process DERIVE_INPUT_INDICES` — `bin/derive_input_design_region.py`.
- `process NEGSTEER_CONTROLS` — `bin/build_control_sequences.py` to generate scrambled+polyA, then `bash bin/negative_steering_run_one.sh` for each.

### modules/negsteer_manifest.nf

- `process EXTRACT_SURVIVOR_MANIFEST` — `bin/extract_survivor_manifest.py`.

### modules/negsteer_interface_metrics.nf

- `process NEGSTEER_INTERFACE_METRICS` — `bin/compute_interface_metrics.py`.

### modules/negsteer_af3_nomsa.nf

- `process AF3_SETUP_DB` — idempotent build of AF3 DB symlink farm.
- `process AF3_NOMSA_ON_SURVIVORS` — per-survivor AF3-no-MSA prediction (Singularity).
- `process AF3_PARSE_OUTPUT` — `bin/parse_af3_output.py` (passed as `path parse_script` after Bug-7 cache fix).

### modules/negsteer_biophysical_metrics.nf

- `process NEGSTEER_BIOPHYSICAL_METRICS` — `bin/run_biophysical_metrics.py`.

### modules/negsteer_rosetta_metrics.nf

- `process NEGSTEER_ROSETTA_METRICS` — `bin/run_rosetta_metrics.py` with `bin/fastrelax_for_ia.xml` passed as `path fastrelax_xml`.

### modules/negsteer_orthogonal_metrics.nf

- `process NEGSTEER_ORTHOGONAL_METRICS` — `bin/merge_orthogonal_metrics.py`.
- `process ORTHOG_PLOTS` — `bin/orthogonal_metrics_plots.py` (passed as `path plot_script`).

### main.nf (top-level workflow)

Single `workflow {}` block. Imports 31 processes (counted from include statements) and invokes them in topological order:

1. Preprocessing branch: `EXTRACT_SEQUENCES` / `RESOLVE_CONTIGS` / `WRITE_DUMMY_MAPPING_REC` / `WRITE_DUMMY_MAPPING_EFF`.
2. Branch A (HADDOCK): `HADDOCK3_PREPARE` → `HADDOCK3_DOCK` → `HADDOCK3_PLOTS` → `EXTRACT_HOTSPOTS` → `BUILD_CONTIGS`.
3. RFDiffusion: `RFDIFFUSION` → `RFDIFFUSION_FILTER` → `RFDIFFUSION_PLOTS`.
4. Rosetta: `ROSETTA_SC` → `ROSETTA_FILTER` → `ROSETTA_FILTER_PLOTS`.
5. ProteinMPNN: `MPNN_FIXED_POSITIONS` → `PROTEINMPNN` → `SEQUENCE_CORRECTION` → `SEQUENCE_QC` → `MPNN_DESIGN_REGION_SCORE` → `MPNN_CLUSTER` → `MPNN_PLOTS` → `MPNN_SELECT_TOP`.
6. Negsteer: `NEGSTEER_DERIVE_INDICES` → `NEGSTEER_RUN_ONE` (fan-out per MPNN sequence) → controls (`DERIVE_INPUT_INDICES` + `NEGSTEER_CONTROLS`) → `NEGSTEER_CROSS_SEQUENCE` → plots.
7. Orthogonal: `NEGSTEER_INTERFACE_METRICS` → `EXTRACT_SURVIVOR_MANIFEST` → fan-out (`AF3_*` / `NEGSTEER_BIOPHYSICAL_METRICS` / `NEGSTEER_ROSETTA_METRICS`) → `NEGSTEER_ORTHOGONAL_METRICS` → `ORTHOG_PLOTS`.

### Test workflows

Each test `.nf` re-uses the production `process` blocks (via include), wires them into a small workflow, and runs against fixed test inputs.

- `tests/haddock/test_haddock.nf` — preprocessing dummy mapping + HADDOCK chain.
- `tests/rfdiffusion/test_rfdiffusion.nf` — RFDIFFUSION + FILTER + PLOTS.
- `tests/rosetta_filtering/test_rosetta_filtering.nf` — ROSETTA_SC + FILTER + PLOTS.
- `tests/proteinmpnn/test_proteinmpnn.nf` — full MPNN sub-chain incl. preprocessing.
- `tests/negative_steering/test_negative_steering.nf` — negsteer chain incl. controls + plots; defines a local `EMPTY_DESIGN_REGION_PLACEHOLDER` process for when controls are off.
- `tests/orthogonal_metrics/test_orthogonal_metrics.nf` — manifest + AF3 + biophysical + Rosetta + merge.
