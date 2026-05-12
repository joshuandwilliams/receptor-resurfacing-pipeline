# 06 — Ubiquitous Language

The single canonical reference for domain terminology used in this project. Each entry defines a concept at the domain level (what the thing *is*), lists the code variants used to refer to it, and notes synonyms, cross-references, and known ambiguities.

This is the language the pipeline *should* use consistently. Where the current code is inconsistent — using one name for two concepts, or two names for one concept — those are flagged in §10. Phase 3 cleanup should resolve them.

**Scope**: terms come from the domain (computational structural biology of de novo binder design + Boltz validation), not from code conventions. Code variants are listed only as recognition aids; they should not be read as endorsing the names. Where a term has multiple casings or punctuations (e.g. `design_region`, `designRegion`, `design-region`), only the underscored Python form is listed unless the variant matters semantically.

**Terminology authority**: where Baker-lab / RFDiffusion-paper terminology (Watson et al. 2023, "motif scaffolding") differs from informal project usage, this glossary adopts the published terminology. Deviations from field standard are flagged.

---

## Contents

1. Tools and prediction engines
2. Pipeline stages
3. Biological and structural concepts
4. Identifiers and granularity
5. Verdicts, tiers, and outcomes
6. Metrics and scores
7. Mechanics and orchestration
8. File artefacts
9. Conventions and invariants
10. Flagged ambiguities for cleanup

---

## 1. Tools and prediction engines

### RFDiffusion
Generative diffusion model for protein backbones. Used here to generate novel receptor backbones around a fixed effector pose, producing N candidate designs per run.

- **Variants in code**: `rfdiffusion`, `RFDIFFUSION`, `RFDIFFUSION_FILTER`, `rfdiff_*`

### ProteinMPNN (often abbreviated MPNN)
Sequence-design model that proposes amino-acid sequences for a fixed backbone. Run on each surviving RFDiffusion design to produce K candidate sequences.

- **Variants in code**: `proteinmpnn`, `PROTEINMPNN`, `mpnn_*`, `MPNN`

### HADDOCK / HADDOCK3
Restraint-driven docking program. Used in Branch A to dock receptor + effector PDBs into an input complex when no pre-docked complex is available.

- **Variants in code**: `haddock`, `HADDOCK`, `haddock3_*`, `HADDOCK3_*`

### Boltz / Boltz2
Structural prediction model used as the core predictor for the negative-steering validation stage. Every Boltz prediction yields a complex PDB plus confidence metrics (pLDDT, PAE, ipTM, etc.).

- **Variants in code**: `boltz`, `boltz2`, `Boltz`, `BOLTZ2_*`

### AlphaFold3 (AF3)
Independent structure predictor used as an orthogonal cross-check against Boltz on cohort survivors. **In this pipeline AF3 is always run without an MSA** — a deliberate choice to ensure orthogonality with Boltz. The "AF3-no-MSA" qualifier is therefore implicit; bare "AF3" in this codebase always means the no-MSA configuration.

- **Variants in code**: `af3`, `AF3`, `alphafold3`, `AF3_NOMSA_ON_SURVIVORS`, `AF3_PARSE_OUTPUT`

### Rosetta
Umbrella term for the Rosetta protein modelling suite. Three distinct Rosetta sub-tools are used (FastRelax, InterfaceAnalyzer, RosettaScripts).

- **Variants in code**: `rosetta`, `ROSETTA_*`

### FastRelax (Rosetta protocol)
Rosetta protocol that minimises a structure under a force field prior to scoring. Pinned at 1 repeat Cartesian `ref2015_cart` to match the calibration target of Bennett et al. 2023.

- **Variants in code**: `fastrelax`, `_run_fast_relax`, `fastrelax_for_ia.xml`

### InterfaceAnalyzer (Rosetta module)
Rosetta module that computes interface scores: shape complementarity (Sc), `dG_separated`, `dSASA_int`, and ΔΔG of binding. Run after FastRelax on cohort survivors as one of the orthogonal-metrics streams.

- **Variants in code**: `InterfaceAnalyzer`, `_run_interface_analyzer`, `_parse_ia_scorefile`

### DockQ
Tool/metric that scores predicted vs reference complex structure quality. Computed per cohort row in the interface-metrics step.

- **Variants in code**: `dockq`, `DockQ`, `_run_dockq`

### MMseqs2
Fast sequence clustering tool. Used to cluster MPNN-designed sequences before top-N selection.

- **Variants in code**: `mmseqs`, `mmseqs_bin`, `run_mmseqs_cluster`

### CAPRI / clustfcc
Output-table formats from HADDOCK. CAPRI is the per-model scoring table; clustfcc is the cluster table. Both are parsed by `haddock_utils`.

- **Variants in code**: `parse_capri_tsv`, `parse_clustfcc_tsv`

### ColabFold (historical)
MSA-search tool referenced in `nextflow.config` and `project_map.md` as part of an earlier MSA-based Boltz path. **Not wired into the current workflow**; mention in code is historical.

- **Variants in code**: `colabfold`, `ColabFold`, `COLABFOLD_*`

---

## 2. Pipeline stages

Stages listed in workflow order. The pipeline is one Nextflow DSL2 workflow (`main.nf`) that runs preprocessing → optional HADDOCK docking → backbone generation → Rosetta filtering → sequence design → negative-steering validation → orthogonal-metrics validation → plotting.

### Preprocessing
Sequence extraction from input PDBs and resolution of contigs before backbone generation.

- **Processes**: `EXTRACT_SEQUENCES`, `RESOLVE_CONTIGS`, `WRITE_DUMMY_MAPPING`

### Branch A (HADDOCK docking)
Optional entry path: dock two separate PDBs (receptor + effector) using HADDOCK3 to produce the input complex. Contigs, AIR restraints, hotspot extraction, and re-construction of contigs all happen in this branch.

- **Processes**: `HADDOCK3_PREPARE`, `HADDOCK3_DOCK`, `HADDOCK3_PLOTS`, `EXTRACT_HOTSPOTS`, `BUILD_CONTIGS`

### Branch B (pre-docked complex)
Entry path that takes a pre-docked complex PDB directly, bypassing HADDOCK.

### Backbone generation (RFDiffusion stage)
RFDiffusion generates N candidate receptor backbones around the fixed effector pose. Each design is then split (chain A = receptor, chain B = effector) and scored.

- **Processes**: `RFDIFFUSION`, `RFDIFFUSION_FILTER`, `RFDIFFUSION_PLOTS`

### Rosetta filtering (pre-MPNN)
Pre-MPNN physics filter: drops backbones whose interface shape complementarity (Sc) is below threshold.

- **Processes**: `ROSETTA_SC`, `ROSETTA_FILTER`, `ROSETTA_FILTER_PLOTS`

### Sequence design (ProteinMPNN stage)
ProteinMPNN proposes K sequences per surviving backbone. Includes correction (reconstructing the full chimeric receptor sequence from MPNN's per-region output), QC, design-region rescoring, clustering, and top-N selection.

- **Processes**: `MPNN_FIXED_POSITIONS`, `PROTEINMPNN`, `SEQUENCE_CORRECTION`, `SEQUENCE_QC`, `MPNN_DESIGN_REGION_SCORE`, `MPNN_CLUSTER`, `MPNN_SELECT_TOP`, `MPNN_PLOTS`

### Negative-steering validation
Per-MPNN-sequence validation. Predict the design with Boltz, identify residues at the "wrong" predicted interface, mutate them to push Boltz off the wrong site, re-predict, then run a reversion pass to test whether the pose holds when the steering mutations are reverted. The full chain (plan → predict-one array → collect → prefilter → build-contaminated → plan-reversions → predict-reversion array → harvest-reversions → finalize → aggregate → compute-final-metrics → aggregate-per-sequence → extract-passing) runs as one SLURM GPU job per MPNN sequence.

- **Processes**: `NEGSTEER_DERIVE_INDICES`, `NEGSTEER_RUN_ONE`, `NEGSTEER_CROSS_SEQUENCE`, `NEGSTEER_PLOTS`, `NEGSTEER_WITHIN_SEQUENCE_PLOTS`
- **Variants in code**: `negative_steering`, `negsteer`

### Multi-cycle steering (currently dormant)
Repeated steering cycles where each cycle starts from the validated reverted sequence of the previous cycle. Implemented in `boltz2_iterate_steering.py` as a 14-subcommand orchestrator but **not wired into the current Nextflow integration** — the current pipeline only runs single-cycle (`--n-cycles 1`). Multi-cycle is on the roadmap.

- **Variants in code**: `iterate_steering`, `iterate-plan`, `kickoff_*`, `--n-cycles`

### Reversion pass (sub-stage of steering)
Sub-stage that takes each steered design, identifies steering mutations *inside the protected set* (design region ∪ true interface), reverts them back to the MPNN design's original residues, and re-predicts with Boltz. The reverted prediction is the *final* design — even if reversion fails, its metrics are what the cohort summary reports. Tests whether the pose holds without the load-bearing steering mutations.

- **Variants in code**: `reversion`, `plan-reversions`, `harvest-reversions`, `reverted_*`, `build_reverted_sequence`

### Negative controls
Two synthetic non-design sequences (`scrambled`, `polyA`) generated and run through the same negsteer chain as a real MPNN sequence, to test that the framework rejects nonsense receptors. Note: cold-start predictions on controls typically have no contact residues to mutate, so the steering phase silently skips — controls effectively become "3 cold-start predictions, no steered designs" and never reach tier A. This is semantically correct for testing the prediction-and-ranking framework, but the name "negative steering control" is mildly misleading because steering rarely fires.

- **Processes**: `DERIVE_INPUT_INDICES`, `NEGSTEER_CONTROLS`
- **Variants in code**: `controls`, `control_scrambled`, `control_polyA`

### Cross-sequence aggregation
Cohort-level join of every per-sequence `passing_summary.csv` into one tier-ranked CSV with composite-score sorting.

- **Process**: `NEGSTEER_CROSS_SEQUENCE`
- **Variants in code**: `cross_sequence_summary`

### Orthogonal-metrics validation
Three independent-method checks (AF3-no-MSA, biophysical, Rosetta) on cohort rows that have a representative PDB on disk, plus a per-row interface-metrics step (iRMSD, fnat, DockQ, weighted Jaccard).

- **Processes**: `EXTRACT_SURVIVOR_MANIFEST`, `NEGSTEER_INTERFACE_METRICS`, `AF3_*`, `NEGSTEER_BIOPHYSICAL_METRICS`, `NEGSTEER_ROSETTA_METRICS`, `NEGSTEER_ORTHOGONAL_METRICS`, `ORTHOG_PLOTS`
- **See also**: §10 — Survivor (term flagged for rename)

---

## 3. Biological and structural concepts

### Target
A protein one wishes to bind. General term in the binder-design field.

### Effector
The specific target protein this pipeline binds against. **Effector is an instance of target.** Used overwhelmingly in the codebase. After RFDiffusion, hardcoded to chain B in the output PDB convention.

- **Variants in code**: `effector`, `eff`, `eff_chain`, `PRED_EFF_CHAIN`
- **See also**: Target (parent concept)
- **Note**: "ligand" is a synonym for effector in some literature but is **not used** in this project.

### Binder
A protein designed (or selected) to bind a target. General term in the binder-design field.

### Receptor
The specific binder protein this pipeline redesigns. **Receptor is an instance of binder.** Used overwhelmingly in the codebase. After RFDiffusion, hardcoded to chain A in the output PDB convention.

- **Variants in code**: `receptor`, `rec`, `rec_chain`, `PRED_REC_CHAIN`
- **See also**: Binder (parent concept), §9 chain-A/B convention

### Native (≡ wild-type)
The protein as it arose through evolution — i.e. *not* designed. In this project, the native receptor is the starting structure that RFDiffusion generates designs around. "Native" and "wild-type" are synonyms.

- **Variants in code**: `native_*`, `wild_type_*`, `wt_*`

### Complex
A PDB containing receptor + effector together. Either pre-docked (Branch B input) or HADDOCK-docked (Branch A output).

- **Variants in code**: `complex`, `complex_pdb`, `input_pdb`, `input_complex`

### Backbone
The polypeptide main chain: **N, Cα, C, and O atoms** along the polymer. Side chains are everything else hanging off Cα. RFDiffusion generates backbone atoms in this sense (N/Cα/C/O).

- **Variants in code**: `backbone` (in narrative comments and prose)
- **See also**: Cα trace (the Cα-only abstraction)
- **Audit note**: prose in this codebase that says "backbone" while operating on `read_ca_atoms` output is computing Cα-trace operations, not full-backbone operations. Flagged in §10.

### Cα trace
The sequence of Cα atoms only — a coarse representation that ignores the rest of the backbone (N/C/O) and all side-chain atoms. **Most geometry helpers in this codebase operate on Cα traces** (`read_ca_atoms`, `_chain_ca_coords`, `_kabsch_align`), not full backbones.

- **Variants in code**: `read_ca_atoms`, `read_ca_seq_*`, `ca_coords`, `_chain_ca_coords`, `CAEntry`

### Design
One RFDiffusion-generated backbone (with its derived metadata: contigs, motif boundaries, design-region indices, true-interface indices, per-design metrics). Identified by `design_id` / `design_idx` (e.g. `d19`).

- **Variants in code**: `design`, `design_id`, `design_idx`, `RFDiffusion design`, `d{N}`

### Motif (Baker-lab framing)
Receptor residues **preserved** from the native sequence/structure during RFDiffusion. The fixed anchor regions of the contig — the parts that are NOT regenerated. Aligned with the "motif scaffolding" terminology of Watson et al. 2023.

- **Audit note**: this codebase has historically used `scaffold` for this concept (e.g. `calc_scaffold_and_region_metrics`, `scaffold_rmsd`). Under the Baker-lab framing adopted here, that usage is **inverted** — the code's `scaffold_*` actually measures the motif. **Strongly flagged for rename in §10.**

### Scaffold (≡ Design region)
Receptor residues **generated de novo** by RFDiffusion — the new structure built around the motif. Synonym of design region; "design region" is the preferred term in this glossary because it is unambiguous in code, while "scaffold" is the upstream Baker-lab term and may appear in upstream documentation.

- **See also**: Design region (preferred name), Motif (the complement)
- **Note**: existing code uses `scaffold` for the *motif* (the inverse meaning). When updating code, use `design_region` for the de novo region; reserve `motif` for the preserved region. **Strongly flagged in §10.**

### Design region
The set of receptor residue positions that ProteinMPNN was allowed to redesign. These are the de novo-generated positions, expressed as 1-based residue indices into the receptor sequence.

- **Variants in code**: `design_region`, `design_region_indices`, `design_region_positions`
- **Synonym**: Scaffold (Baker-lab framing)
- **Complement within receptor**: Motif

### Contig
RFDiffusion's specification string for a design: a comma-separated list of segments where each segment is either an *anchor* (residue range copied from native, contributing to the motif) or *de novo* (range to be generated, contributing to the design region). Resolved at preprocessing time.

- **Variants in code**: `contig`, `contigs`, `contigmap`, `parse_block_segments`, `resolve_contigs`

### De novo region
A segment in a contig string that RFDiffusion generates from scratch. The union of de novo segments in a design is the design region.

- **Variants in code**: `denovo`, `de novo`, `find_denovo_residues`

### Anchor
A native-sequence residue range that is pinned in place during design — i.e. part of the motif. The MPNN-corrected sequence is aligned back to native via these anchors.

- **Variants in code**: `anchor`, `align_to_native_by_anchors`, `native_anchor_regions`

### Hotspot
Receptor residues identified by HADDOCK as participating in the docked interface. Consumed downstream to build RFDiffusion contigs (the de novo region is sized around the hotspots).

- **Variants in code**: `hotspot`, `hotspot_residues`, `extract_hotspots`

### True interface
The set of receptor residue positions that contact the effector in the design's *intended* binding mode (i.e. as Boltz would predict if it predicted the design correctly). Stored as 0-based indices. Used as the gold-standard target the steering machinery is trying to hit.

- **Variants in code**: `true_interface`, `true_interface_idx`, `true_interface_indices`
- **Note**: stored 0-based by historical convention; design region is stored 1-based. The index-base difference is a code accident, not a domain distinction. See §9 and §10.

### Wrong interface
The set of receptor residue positions Boltz predicts as contacting the effector in the cold-start (pre-steering) prediction — i.e. the residues at the *unintended* binding mode the steering pass aims to push Boltz away from.

- **Variants in code**: `wrong_interface`, `wrong_interface_residues`, `initial_wrong_interface_idx`

### Protected set
The union (design region ∪ true interface): receptor positions where mutation contamination matters. The reversion pass reverts steering mutations only when they sit inside the protected set; mutations outside it are accepted as legitimate steering contacts.

- **Variants in code**: `protected_set`, `contamination_gating_positions`

### Chimeric receptor
A receptor sequence assembled by splicing native (anchored) and MPNN-designed (de novo) fragments — i.e. the corrected output of `pipeline_correct_sequences.py`.

- **Variants in code**: `chimera`, `chimeric`, `corrected_receptor`

### Pose
A specific predicted 3D placement of the effector relative to the receptor.

- **Variants in code**: `pose`, `pose_distance`, `pose_holds`, `pose_collapses`

### Cold-start prediction
A Boltz prediction made on the un-mutated MPNN sequence (no steering mutations applied yet). The starting reference point of every negsteer cycle.

- **Variants in code**: `cold_start`, `cold_baseline`, "initial Boltz prediction"

### Steered design
A receptor sequence with mutations applied at wrong-interface positions to push Boltz off the wrong site.

- **Variants in code**: `steered`, `steered_*`, `make_steered_sequence`

### Reverted design
A steered sequence with the protected-set steering mutations reverted back to the MPNN design's original residues. The reverted prediction is the *final* design — Boltz's prediction on this sequence is what's reported in cohort summaries when reversion was attempted.

- **Variants in code**: `reverted`, `reverted_*`, `build_reverted_sequence`

### Steering set / steering mode
A named pool of allowed substitutions used by the mutator to pick steering residues. **Four modes are valid and distinct**:

- **strong** — large, strongly charged, or bulky side chains. Pool: `WYFRKEDP`.
- **mild** — charge changes without large or awkward side-chain reshapes that might break the receptor. Pool: `DEKR`.
- **conservative** — fixed-dictionary substitutions that preserve approximate size/charge while flipping chemistry (e.g. K→E, A→S). Defined in `CONSERVATIVE_SUBSTITUTIONS`.
- **alanine** — substitute the residue with alanine.

- **Variants in code**: `STEERING_SETS`, `CONSERVATIVE_SUBSTITUTIONS`, `steering_mode`

### Construct (the plasmid)
The plasmid containing the designed receptor sequence, used for agroinfiltration in the wet-lab follow-up. A domain noun referring to a physical reagent.

- **Note**: distinct from `construct_reliance_flag` (see Reliance flags). The names overlap but refer to different things.

### Reliance flags (mutation predicates)
Computed predicates on cohort rows. Two granularities, both produced by `classify_mutation_reliance`:

- **`mutation_reliance_flag`** — the predicted pose depends on at least one mutated position contacting the effector.
- **`construct_reliance_flag`** — the predicted pose depends on a mutated position *inside the protected set* — i.e. on a mutation that would not exist in the wet-lab construct (because the reversion pass would have removed it). This is the wet-lab-relevant flag.

- **Variants in code**: `mutation_reliance_flag`, `construct_reliance_flag`, `classify_mutation_reliance`
- **Note**: the `construct_` prefix here is unrelated to Construct (the plasmid). Two different uses of the same root.

---

## 4. Identifiers and granularity

The pipeline operates at four distinct granularities, finest first.

### Diffusion sample
One AF3 per-sample output (5 per seed, organised as `seed-N_sample-M/` directories). Specific to AF3 output structure; not a unit elsewhere in the pipeline.

- **Variants in code**: `diffusion_samples`, `negsteer_diffusion_samples`

### Seed
One Boltz (or AF3) prediction trial — a single random-seeded run of the predictor. Typically 3 seeds per `(sequence, sequence_group)` pair.

- **Variants in code**: `seed`, `seed_idx`, `seed_index`, `num_seeds`

### Sequence group (sg)
A computational grouping of seeds that share the same input receptor sequence under one design. **No biological meaning** — the term comes from the standalone negsteer pipeline as a bookkeeping device. Aggregators group rows by the `sequence_group` integer.

- **Variants in code**: `sg`, `sg_idx`, `sequence_group`, `sgs`

### Sequence (cohort sense)
One MPNN-designed sequence for one RFDiffusion design — the unit one row of `cross_sequence_summary.csv` represents. Identified in plots and CSVs by names like `d19_s1` (design 19, sequence 1).

- **Variants in code**: `seq_name`, `mpnn_sequence`, `MPNN sequence`, `s{N}`

### Representative seed
The single seed chosen to represent a (sequence, sg) pair in cohort plots and the cross-sequence CSV — typically the seed with the best composite score.

- **Variants in code**: `rep_sg`, `rep_*` column prefix, `_rep_sg_for_outcomes`

### Cohort
The full set of MPNN sequences (plus the two negative controls) processed together in one pipeline run — the rows of `cross_sequence_summary.csv`.

- **Variants in code**: `cohort`, `load_unified_cohort`, `cross_sequence_*`

### Cycle (multi-cycle only — currently dormant)
One iteration of the multi-cycle steering loop. The current Nextflow integration runs only `cycle_0` (single-cycle); cycles 1+ exist in `boltz2_iterate_steering.py` as dormant orchestration code.

- **Variants in code**: `cycle`, `n_cycles`, `cycle_0`

### Pathway (multi-cycle only — currently dormant)
The lineage of steering mutations across multiple cycles for one design. Encoded as a dot-separated label `c0d05.c1d03` where each `cN dM` segment names one cycle's design choice.

- **Variants in code**: `pathway`, `pathway_label`, `make_pathway_label`

### Leaf label vs parent label (multi-cycle only)
A pathway label is itself a leaf label — the full chain identifying the tip of the lineage (e.g. `c0d05.c1d03`). The parent label is the prefix obtained by dropping the last segment (`c0d05`); it names the directory the leaf's design lives inside (`pathway_<parent>`).

- **Variants in code**: `leaf_label`, `parent_pathway_label`

### Cumulative mutations (multi-cycle only)
The union of mutated positions accumulated across all ancestors in a pathway — the full set of receptor positions that have been mutated to reach a given leaf prediction.

- **Variants in code**: `accumulated_mutated_positions`, `read_cumulative_mutations`

---

## 5. Verdicts, tiers, and outcomes

The pipeline applies three layered classifiers: a **per-seed reversion verdict** (after the reversion pass), a **per-sequence outcome**, and a **cohort tier**. All three feed into the final ranking.

### Reversion verdict (per-seed)
The classification produced by `classify_reversion_verdict` for one seed's reverted prediction.

- **`pose_holds`** — reverted prediction is intact and on-target with no new contamination. The successful outcome.
- **`pose_collapses`** — reverted prediction failed the structural filter (intact + ra_eff threshold). Reversion broke the pose.
- **`new_contamination`** — reverted prediction has a mutated residue contacting the effector inside the protected set. Reversion didn't fix the contamination.
- **`no_data`** — reverted prediction is missing or unparseable.

- **Variants in code**: `verdict`, `reversion_verdict`, `classify_reversion_verdict`

### `clean_steered` (per-seed)
A seed whose **steered** prediction had zero mutated-position contamination AND was structurally OK (intact + ra_eff < 5Å). Reversion was correctly skipped for this seed because there was nothing to revert. Counts as pass-equivalent in tier classification.

- **Variants in code**: `clean_steered`, `_row_is_clean_steered`; contributes to `n_pass`
- **See also**: `no_reversion` (the per-sequence outcome label that fires when every seed in the group was clean_steered or cold-start-skipped — same underlying state at different granularities; outcome=`no_reversion` is NOT a tier proxy)

### Outcome (per-sequence, renamed from `aggregated_verdict`)
The sequence-level label computed from per-seed reversion verdicts. Possible values include the four reversion-verdict classes above plus:

- **`no_reversion`** — no seed in this group needed reversion. Fires either when (a) all seeds passed cold-start `skip_steering`, or (b) steering ran but every seed's steered prediction had zero contamination. Note: **outcome `no_reversion` does NOT imply n_pass == n_seeds.** A wrong-placement-no-contamination group is `no_reversion` with n_pass == 0 → tier none. Use `n_pass / n_seeds` for tiering.
- **`singleton`** — passthrough row (cycle-0 initial baseline); not aggregated.

- **Variants in code**: `outcome`, `outcome_reason`, `_classify_outcome` (formerly `aggregated_verdict` / `_classify_aggregated_verdict`)

### `n_pass` (per-sequence)
The number of seeds in a sequence-group that passed the structural+contamination filters: `n_pass = (# pose_holds seeds) + (# clean_steered seeds)`. **Path-agnostic** — computed the same way regardless of whether reversion ran. Drives cohort tier directly. Replaces the five `n_seeds_*` per-verdict count columns from the pre-rename schema.

### `PASSING_OUTCOMES`
The set of outcome labels allowed to populate `passing_summary.csv`: `{pose_holds, no_reversion, ""}`. The empty string covers rows that haven't been outcome-classified.

### Tier (cohort A/B/C/none)
The cohort-level promotion class for one MPNN sequence, derived **only** from `n_pass / n_seeds`:

- **Tier A** — every seed passes: `n_pass == n_seeds`.
- **Tier B** — strict majority but not all: `1 < n_pass < n_seeds`.
- **Tier C** — exactly one seed passes: `n_pass == 1`.
- **Tier none** — `n_pass == 0`.

The outcome label is informational only for tier assignment.

- **Variants in code**: `tier`, `cross_tier`, `_tier_for_row`, `_TIER_ORDER`

### Stage (in seed-outcome plots)
The negsteer-chain phase that supplied a given seed's metrics: `cold_start`, `steering`, or `reversion`. Determined by which column family (`steered_*` vs `reverted_*`) and which row state populated the plot row.

- **Variants in code**: `stage`, `STAGE_MARKER`, `STAGE_COLOUR`

### Outcome (in seed-outcome plots)
The per-seed pass/fail decomposition along two axes: **structural pass** (intact + ra_eff < threshold) × **confidence pass** (pLDDT, iPAE, PAE-pass-frac, iPTM all above thresholds), plus the contamination/no-data buckets. Outcomes: `pass`, `off_target`, `poor_prediction`, `multiple_failures`, `new_contamination`, `no_data`.

- **Variants in code**: `outcome`, `OUTCOME_LABEL`, `OUTCOME_COLOUR`, `OUTCOME_ORDER`

### Composite score
`true_jaccard_median − 0.05 × ra_eff_vs_truth_median`. The headline ranking score (higher is better). Used to rank within tiers and to break ties in cohort plots.

- **Variants in code**: `composite`, `composite_score`, `_composite_score`, `COMPOSITE_RA_EFF_WEIGHT=0.05`

### `passes_orthogonal_filters`
The final gate after the orthogonal-metrics merge. True iff biophysical and Rosetta orthogonal metrics are above their thresholds. AF3 is computed but **demoted to flag-only** in the merge step (production code in `bin/merge_orthogonal_metrics.py` had not been updated to match `tests/orthogonal_metrics/merge_orthogonal_metrics.py` as of this writing — see `05_findings.md` A4).

- **Variants in code**: `passes_orthogonal_filters`, `_apply_filters`, `orthogonal_flags`

---

## 6. Metrics and scores

### 6.1 Confidence metrics (predictor-internal)

These come from Boltz / AF3 output per prediction.

- **pLDDT** — per-residue confidence score; aggregated as `avg_plddt`, `complex_plddt`, `interface_plddt`.
- **PAE** — predicted aligned error matrix between residue pairs.
- **`pae_pass_frac` (paepf)** — fraction of inter-chain residue pairs whose PAE is below cutoff.
- **iPAE** — interface-restricted PAE summary.
- **iPSAE** — interface pSAE-style score. A 15 Å variant is also computed in the interface-metrics step.
- **iPTM** — interface predicted-TM score.
- **pTM** — whole-structure predicted-TM score.
- **actifPTM** — active-interface-restricted pTM variant.

### 6.2 Geometry metrics

- **Receptor-aligned effector RMSD (`ra_eff`)** — after Kabsch-aligning predicted receptor onto reference receptor (Cα trace), RMSD of the predicted effector vs reference effector. The headline geometry metric of the steering pass; threshold typically 5 Å. Variants: `ra_eff`, `receptor_aligned_effector_rmsd`, `ra_eff_vs_truth`.
- **Independent receptor RMSD** — receptor-only RMSD vs reference (Cα trace). Used to detect whether the receptor's own fold survived steering.
- **Independent effector RMSD** — effector-only RMSD vs reference.
- **Receptor intact / effector intact** — boolean: independent RMSD below 5 Å (`RECEPTOR_INTACT_CUTOFF=5.0`, `EFFECTOR_INTACT_CUTOFF=5.0`). Used as a structural-pass gate.
- **iRMSD** — CAPRI-style interface RMSD.
- **fnat** — fraction of native interface contacts recovered.
- **DockQ** — CAPRI-derived combined metric of interface quality.
- **Pose distance** — receptor-aligned RMSD between two predicted poses; used by maximin selection in multi-cycle planning.
- **COM displacement** — centre-of-mass shift of the design region between native and design.
- **Kabsch alignment** — optimal rigid-body superposition used for all RMSD computations (`_kabsch_align`).

### 6.3 Interface metrics

- **Contact residues** — heavy-atom-distance-defined residues at the predicted interface. Computed by `find_contact_residues_heavy`. Two definitions of this name exist (see §10).
- **True jaccard** — Jaccard overlap between predicted contact residues and the true-interface set. Stored medianised over seeds.
- **Wrong jaccard** — Jaccard overlap between predicted contacts and the wrong-interface set (the cold-start interface).
- **Weighted jaccard** — distance-weighted Jaccard variant using Gaussian contact weights (`_WJ_MU=4.0`, `_WJ_TWO_SIGMA_SQ=2.25`).
- **Intact core** — pLDDT-trimmed RMSD-filter check on the receptor core.
- **BSA** — buried surface area (FreeSASA).
- **Sc (shape complementarity)** — Rosetta interface shape-complementarity score; used as the pre-MPNN filter.
- **dG_separated, dSASA_int** — Rosetta InterfaceAnalyzer interface-energy and ΔSASA outputs.
- **ΔΔG (`rosetta_ddg`)** — Rosetta InterfaceAnalyzer ΔΔG of binding (post-FastRelax). Calibrated against the −30 REU threshold of Bennett et al. 2023.
- **Interface H-bonds** — MDAnalysis-derived hydrogen-bond count at the predicted interface.

---

## 7. Mechanics and orchestration

### Plan / plan phase
First phase of a steering cycle: predict cold-start → identify wrong-interface residues → emit N steered design YAMLs.

- **Variants in code**: `cmd_plan`, `iterate-plan`, `plan_meta`, `reversion_plan`, `plan.json`

### Predict-one
Per-array-task Boltz invocation for one steered (or one reverted) design.

- **Variants in code**: `cmd_predict_one`, `predict-reversion-one`

### Collect
Phase that gathers per-design predictions and ranks them, emitting `steered_results.csv`.

- **Variants in code**: `cmd_collect`, `iterate-collect`

### Prefilter
Eligibility gate run before reversion: drops candidates that fail the upstream structural / confidence filters.

- **Variants in code**: `prefilter`, `iterate-collect-prefilter`, `kickoff-prefilter`

### Build-contaminated
Identifies which steering mutations contact the effector inside the protected set — these are the mutations the reversion pass will attempt to revert.

- **Variants in code**: `build-contaminated`, `cmd_build_contaminated`, `contaminated.json`

### Plan-reversions
Emits per-design reversion YAMLs from the contaminated mutations list.

- **Variants in code**: `plan-reversions`, `cmd_plan_reversions`

### Harvest reversions
Collects reverted predictions, runs `compute_metrics.py` on each, populates per-seed reverted records, and runs `classify_reversion_verdict`.

- **Variants in code**: `harvest-reversions`, `harvest_reversion_results`

### Finalize
Phase that produces the finalised per-cycle CSV with verdicts merged in.

- **Variants in code**: `finalize`, `iterate-collect-finalize`, `kickoff-finalize`

### Aggregate (per-experiment)
Reduces per-prediction rows to per-sequence rows. Two granularities exist within the same chain: `cmd_aggregate` builds the per-prediction flat CSV; `cmd_aggregate_per_sequence` reduces that to per-sequence rows.

- **Variants in code**: `cmd_aggregate`, `cmd_aggregate_per_sequence`, `_aggregate_per_sequence`
- **Note**: an older standalone `aggregate_results.sh` (cross-experiment combiner) existed in the pre-Nextflow era and has been retired (see `project_map.md` and `05_findings.md` A2).

### Compute-final-metrics
Derives the final metric columns + tier classification on the per-sequence CSV.

- **Variants in code**: `cmd_compute_final_metrics`

### Kickoff (multi-cycle only — currently dormant)
Submitter for downstream stages of multi-cycle runs. Five `cmd_kickoff*` family members exist in `boltz2_iterate_steering.py` (~600 LOC); none are reachable via the current Nextflow integration.

- **Variants in code**: `kickoff`, `cmd_kickoff*`

### Maximin selection (multi-cycle only)
Greedy farthest-from-anything-seen pose selection used to pick cycle-N+1 candidates.

- **Variants in code**: `maximin_select`, `compute_pose_distance`

---

## 8. File artefacts

### `rfdiffusion_metrics.json`
Per-design metric blob emitted by RFDiffusion filter; consumed by all three "derive" scripts (`derive_design_region.py`, `derive_true_interface.py`, `derive_input_design_region.py`) to extract design-region and true-interface index sets.

### `plan.json`
Per-sequence plan blob written by negsteer's plan phase. Holds the cold-start prediction metadata, true-interface indices, ground-truth path, effector-template path, and the list of steered designs to predict.

### `contaminated.json`
Per-cycle list of mutations to revert (mutations contacting the effector inside the protected set).

### `steered_results.csv`
Per-sequence ranked steered-prediction summary from `cmd_collect`.

### `aggregated_results.csv`
Per-sequence aggregate of all seed predictions (steered + reverted), with median/majority columns.

### `passing_summary.csv`
One canonical row per MPNN sequence after tier classification — the unit `cross_sequence_summary` aggregates from. Filtered to `PASSING_VERDICTS` when written.

### `cross_sequence_summary.csv`
Cohort-level join of all per-sequence passing summaries. The canonical "cohort table" — driver of plots and the survivor manifest.

### `survivor_manifest.csv`
Fan-out CSV with one row per cohort row that has a representative PDB on disk. Drives the orthogonal-metrics streams. **Term flagged in §10 — current usage of "survivor" does not match its connotation.**

### `all_results_multicycle.csv`
Per-experiment flat CSV of every prediction across all cycles (multi-cycle output). Currently only `cycle_0` rows are produced.

### AIR restraints
HADDOCK ambiguous interaction restraints derived from de novo regions in the contig string.

### Effector template CIF
Effector-only mmCIF passed to Boltz to anchor the effector at a known pose.

### Fixed-positions JSONL
ProteinMPNN input file specifying which residues are NOT to be redesigned (i.e. the motif positions).

### A3M / FASTA
Sequence-format files consumed by Boltz. A single-sequence A3M is used as the "MSA" when no real MSA is available.

### Canonical PDB
The single PDB chosen as the published representative for one cohort row. Path is stored in `rep_canonical_pdb` and is rewritten from work-dir form to published-tree form by `_rewrite_workdir_path_to_published`.

### Sidecar
A companion file accompanying a prediction (e.g. the canonical PDB next to a Boltz CIF output).

---

## 9. Conventions and invariants

### A=receptor / B=effector chain convention
After RFDiffusion, every output PDB is **always** labelled chain A = receptor, chain B = effector, regardless of the chain naming in the input. Hardcoded by `write_split_pdb`. Every downstream consumer (Boltz YAML composer, sequence extractors, contact finders) treats these labels as fixed.

- **Constants**: `PRED_REC_CHAIN="A"`, `PRED_EFF_CHAIN="B"`

### Index-base convention (1-based vs 0-based)
- Design region indices are stored **1-based** (e.g. `derive_design_region.py`, `design_region_positions`).
- True interface indices are stored **0-based** (e.g. `derive_true_interface.py`, `true_interface_idx`).

This split is a code accident, not a domain distinction. Any consumer that uses both must re-base one. **Flagged in §10.**

### Bennett 2023 calibration
The −30 REU threshold for `rosetta_ddg` is calibrated against the FastRelax → InterfaceAnalyzer flow described in Bennett et al. 2023. Changes to the FastRelax XML must preserve the calibration (1 repeat, Cartesian, `ref2015_cart`).

### `passing_summary.csv` is filtered
Only verdicts in `PASSING_VERDICTS = {"", "no_reversion", "pose_holds"}` reach `passing_summary.csv` via the normal flow. Code that reads `aggregated_results.csv` directly will see the full verdict set including `pose_collapses` and `new_contamination`.

---

## 10. Flagged ambiguities for cleanup

These are terms whose current code usage is ambiguous, contradictory, or non-standard. Each is a Phase 3 mechanical-cleanup target.

### F1. Survivor — term does not match its connotation (strong flag)

**Domain expectation**: "a cohort row that passed a quality threshold."

**Code reality**: any cohort row whose `rep_canonical_pdb` is a file on disk, with no metric/tier check (`extract_survivor_manifest.py:9-11, :110`).

A failed-metric row with an output PDB on disk is a "survivor" under the current definition. Pre-2026-04-29 the term effectively did imply tier-passage (because tier-none rows had no representative populated), but Bug 3's tier-none-fallback fix in pipeline_notes14 means tier-none rows now also acquire representatives for display, and so also become survivors. A bad prediction with a parseable PDB is now a survivor.

**Recommendation**: rename the file/symbol family (`survivor_manifest`, `EXTRACT_SURVIVOR_MANIFEST`, `_find_workdir`'s "survivor" vocabulary) to a tier-agnostic term such as `cohort_manifest` / `EXTRACT_COHORT_MANIFEST`, OR add an explicit tier-pass filter to the manifest extractor and keep the name. Pending decision.

### F2. Backbone usage in narrative comments (audit)

**Domain definition**: N/Cα/C/O main-chain atoms.

**Code reality**: most "backbone" prose in this codebase refers to operations on `read_ca_atoms` output, which is Cα-only — a Cα trace, not a backbone.

**Recommendation**: audit narrative comments and rename Cα-only operations as `ca_*` consistently. Reserve `backbone_*` for operations that genuinely use N/Cα/C/O.

### F3. Scaffold — inverted usage (strong flag)

**Baker-lab framing (adopted here)**: scaffold = the de novo-generated region (synonym of design region).

**Code reality**: `calc_scaffold_and_region_metrics`, `scaffold_rmsd`, `scaffold_str` use `scaffold` to mean the *non-designed* (motif) region — the inverse meaning.

**Recommendation**: rename code uses of `scaffold` (in their current sense) to `motif`. Avoid `scaffold` as an identifier altogether; prefer `design_region` for the de novo region. Visiting these renames during Phase 3 will require a coordinated touch across `rfdiffusion_filter.py`, `rfdiffusion_plots.py`, and any test scripts.

### F4. Index-base convention (1-based vs 0-based)

**Issue**: design region is 1-based, true interface is 0-based. Consumers that combine the two (e.g. for the protected set) must re-base.

**Recommendation**: pick one convention (suggest 1-based to match the design-region storage and to align with PDB residue numbering) and migrate. Or store both as 0-based and convert at presentation.

### F5. Contig-parsing duplication

Four independent contig parsers exist (`contig_utils.parse_block_segments`, `derive_input_design_region._parse_contigs`, `haddock3_prepare.parse_contig_segments`, `pipeline_correct_sequences.parse_contig_segments`). Each is scoped slightly differently but they share the core grammar. Already flagged in `05_findings.md` C6.

### F6. `find_contact_residues_heavy` defined twice

One in `bin/boltz2_negative_steering.py`, one in `bin/compute_metrics.py`. Different bodies, same name, both exported as named imports. Already flagged in `05_findings.md`.

### F7. Single-chain sequence extraction defined four times

`get_chain_sequence` (`boltz2_negative_steering.py`), `_extract_chain_sequence` (`build_control_sequences.py`), `_extract_chain_seq` (`extract_survivor_manifest.py`), `get_pdb_sequence` (`pipeline_correct_sequences.py`). Same operation, four implementations. See `05_findings.md` C4.

### F8. THREE_TO_ONE amino-acid dict defined four times

In `build_contigs.py`, `boltz2_negative_steering.py` (twice — `THREE_TO_ONE` plus `_THREE_TO_ONE_RMSD`), `extract_hotspots.py`. See `05_findings.md` C7.

### F9. Tier constants defined identically in 3+ places

`COLOUR_TIER`, `_TIER_SORT_ORDER`, `COMPOSITE_RA_EFF_WEIGHT` are duplicated across plot scripts and test copies. See `05_findings.md` C11.

### F10. Construct (plasmid) vs `construct_reliance_flag`

The root `construct` is overloaded: a domain noun (the wet-lab plasmid) vs the prefix on a flag predicate (`construct_reliance_flag`). The flag's name comes from "the wet-lab construct's predicted pose"; readability suffers because the plasmid sense is rarely surfaced in code while the flag is. Recommend keeping the flag name but adding a glossary cross-reference (already present here) and ensuring any new "construct"-named symbol is the flag, never the plasmid.
