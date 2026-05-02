# Project map: Receptor-Resurfacing pipeline + Negative-Steering validator

This project context contains files from **two related but distinct
codebases**, both belonging to the same overall research effort: a
computational pipeline for designing protein binders by **receptor
resurfacing** (RFDiffusion + ProteinMPNN), with a **Boltz2-based
validation stage** at the end. The strategic goal of the project is to
**replace the existing Boltz2 validation stage in the resurfacing
pipeline with the new negative-steering workflow**, because the
negative-steering workflow can detect when a Boltz2 binding-site
prediction is being held up by mutations that won't exist in the
final wet-lab construct.

The two codebases are at different points in their lifecycle:

- The **receptor-resurfacing pipeline** is a mature Nextflow DSL2
  workflow that is currently in a cleanup/audit phase (see
  `pipeline_notes1.md`, `pipeline_notes2.md`).
- The **negative-steering pipeline** is a standalone SLURM/Python
  workflow that has been validated on native benchmark complexes and on
  RFDiffusion designs and is being prepared for integration into the
  Nextflow pipeline (see `notes1.md` through `notes6.md`).

---

## 1. The receptor-resurfacing pipeline (Nextflow)

A Nextflow DSL2 pipeline for designing receptor mutants that bind a
given effector. Two entry branches:

- **Branch A:** two separate PDBs (receptor + effector) → HADDOCK3
  docking → complex PDB.
- **Branch B:** a pre-docked complex PDB.

Both branches feed into RFDiffusion → Rosetta → ProteinMPNN → (optional
ColabFold MSA) → Boltz2 validation → aggregate.

### 1.1 Top-level orchestration

| File | Role |
|---|---|
| `main.nf` | Top-level workflow. Wires every module together, routes branches A/B, includes all process definitions from `modules/`. |
| `nextflow.config` | Single source of truth for container paths (`rfdiff_container`, `rosetta_container`, `boltz2_container`, `colabfold_container`) and per-process SLURM resources. |
| `params_example.yml` | Example parameter file users copy and edit. |
| `run_pipeline_slurm.sh` | SLURM submission wrapper for the full pipeline. |

### 1.2 Per-block Nextflow modules

Each `.nf` module defines one or more Nextflow processes that shell out
to scripts in `bin/` (the Python files listed below).

| Module | Processes | Backing Python scripts |
|---|---|---|
| `preprocessing.nf` | `EXTRACT_SEQUENCES`, `RESOLVE_CONTIGS`, `WRITE_DUMMY_MAPPING` | (inline + helpers) |
| `haddock.nf` | `HADDOCK3_PREPARE`, `HADDOCK3_DOCK`, `HADDOCK3_PLOTS`, `EXTRACT_HOTSPOTS`, `BUILD_CONTIGS` | `haddock3_prepare.py`, `collect_haddock3_dock.py`, `haddock3_plots.py`, `haddock_utils.py`, `extract_hotspots.py`, `build_contigs.py`, `contig_utils.py` |
| `rfdiffusion.nf` | `RFDIFFUSION`, `RFDIFFUSION_FILTER`, `RFDIFFUSION_PLOTS` | `rfdiffusion_contigs.py`, `rfdiffusion_filter.py`, `rfdiffusion_plots.py`, `contig_utils.py` |
| `rosetta_filtering.nf` | `ROSETTA_SC`, `ROSETTA_FILTER`, `ROSETTA_FILTER_PLOTS` | `rosetta_filter_collect.py`, `rosetta_filter_plots.py` |
| `proteinmpnn.nf` | `MPNN_FIXED_POSITIONS`, `PROTEINMPNN`, `SEQUENCE_CORRECTION`, `SEQUENCE_QC`, `MPNN_DESIGN_REGION_SCORE`, `MPNN_CLUSTER`, `MPNN_SELECT_TOP`, `MPNN_PLOTS` | `pipeline_correct_sequences.py`, `mpnn_sequence_qc.py`, `mpnn_design_region_score.py`, `mpnn_cluster_sequences.py`, `mpnn_select_top.py`, `mpnn_plots.py` |
| `msa.nf` | `COLABFOLD_SEARCH_PER_DESIGN` | (inline; runs colabfold_search per MPNN design when `boltz2_use_real_msa=true`) |
| `boltz2.nf` | `BOLTZ2_PREPARE`, `BOLTZ2_PREDICT`, `BOLTZ2_VERIFY_BINDING`, `BOLTZ2_FILTER_AND_RANK`, `BOLTZ2_PLOTS` | `boltz2_prepare.py`, `boltz2_verify_binding.py`, `boltz2_filter.py`, `boltz2_plots.py` |
| `aggregate.nf` | `AGGREGATE_RESULTS` | (inline Python — merges MPNN metadata with Boltz2 ranked results into `results_summary.csv`) |

**This is the Boltz2 stage that the negative-steering pipeline is
intended to replace.** The current stage runs Boltz2 once per MPNN
sequence, computes ipSAE/iPTM/binding-RMSD vs the RFDiffusion design
PDB, filters on ipSAE, and ranks. It cannot detect contamination from
mutations that contact the effector in the predicted pose.

### 1.3 Per-module test infrastructure

Each module has an isolated test harness — a `test_<module>.nf`
workflow plus a `run_test_<module>_slurm.sh` SLURM wrapper that runs
just that module on small fixed inputs. Test scripts include production
modules directly (no inline duplicates) so the test exercises the same
code paths as production.

Pairs:
- `test_boltz2.nf` / `run_test_boltz2_slurm.sh`
- `test_haddock.nf` / `run_test_haddock_slurm.sh`
- `test_proteinmpnn.nf` / `run_test_proteinmpnn_slurm.sh`
- `test_rfdiffusion.nf` / `run_test_rfdiffusion_slurm.sh`
- `test_rosetta_filtering.nf` / `run_test_rosetta_filtering_slurm.sh`

### 1.4 Pipeline session notes

- `pipeline_notes1.md` — Cleanup-pass session 1: AF2-monomer branch
  removal, container-path consolidation, preprocessing/HADDOCK/
  RFDiffusion/Rosetta audits and bug fixes.
- `pipeline_notes2.md` — Cleanup-pass session 2: ProteinMPNN block.
  Includes the `pipeline_correct_sequences.py` `split_sequence` bug fix
  (chain-order swap) and the design-region-score replacement for the
  broken RMSD function.

Both notes flag remaining work on the MSA, Boltz2, and aggregate blocks,
plus a final whole-pipeline pass on `main.nf`.

---

## 2. The negative-steering pipeline (standalone SLURM + Python)

A multi-cycle pipeline that, given a complex where Boltz2 predicts the
wrong binding site, mutates the receptor residues at the wrong
interface, re-predicts, and iterates. Each cycle pushes the effector
away from where it landed before, with the goal of having Boltz2
converge on the true binding mode (`receptor_aligned_effector_rmsd <
3 Å` vs ground truth).

The pipeline has a dual purpose:

1. As a **diagnostic tool** for native benchmark complexes — does
   negative steering rescue Boltz2's mistakes?
2. As a **validation tool** for RFDiffusion+MPNN-designed binders — for
   each candidate, can Boltz2 be steered onto the design's intended
   binding site, and does the rescued pose hold up when the steering
   mutations are reverted?

Use case (2) is what motivates replacing the current `boltz2.nf`
validation stage.

### 2.1 Single-cycle core

| File | Role |
|---|---|
| `boltz2_negative_steering.py` | Single-cycle core (~2360 lines). Subcommands: `plan` (initial Boltz2 prediction → identify wrong-interface residues → emit N steered design YAMLs), `predict-one` (one design per SLURM array task), `collect` (rank and summarise). Stable across all multi-cycle work. |
| `submit_boltz2_negative_steering.sh` | Top-level submitter. Runs single-cycle if `--n-cycles 1`, otherwise writes `params.yml`/`submission.log` and chains a kickoff job after cycle 0. Now a 6-stage chain (was 3) since the reversion pass was added. |

### 2.2 Multi-cycle orchestration

| File | Role |
|---|---|
| `boltz2_iterate_steering.py` | Multi-cycle orchestrator (~2280 lines). Subcommands: `iterate-plan`, `compute-distances`, `iterate-collect-prefilter`, `build-contaminated`, `plan-reversions`, `predict-reversion-one`, `harvest-reversions`, `iterate-collect-finalize`, `kickoff-distances`, `kickoff-prefilter`, `kickoff-finalize`, `aggregate`, `compute-final-metrics`. |
| `submit_boltz2_iterate_steering.sh` | Per-cycle submitter for cycles 1+. 6-stage chain (was 3) as of the reversion-pass addition. |
| `reversion.py` | Pure-logic helpers used by the orchestrator: `build_reverted_sequence`, `write_reversion_plan`, `harvest_reversion_results`, `classify_reversion_verdict`. |

### 2.3 Metrics, filtering, and post-processing

| File | Role |
|---|---|
| `compute_metrics.py` | Computes Boltz confidence metrics (pTM, ipTM, ipSAE, actifPTM, pLDDT, PAE) from Boltz output JSONs. Also computes contact-residue analysis and the two reliance flags: `mutation_reliance_flag` (any contact on a mutated position) and `construct_reliance_flag` (mutated contact inside the protected set = design region ∪ true interface). The latter is the wet-lab-relevant flag. |
| `derive_design_region.py` | Standalone helper. Reads `rfdiffusion_metrics.json` and emits the 1-based positional indices of the ProteinMPNN-designed region for a given design. Output is consumed by `compute-final-metrics --design-region-positions-file`. |
| `compare_versions.py` | Cross-version validator. Scores designs into production / borderline / failed tiers using the locked-in criterion (`intact AND ra_eff < 3.0 AND true_jaccard >= 0.7`). Adds a `tier_reason` column explaining each rejection. |

### 2.4 Multi-experiment aggregation

| File | Role |
|---|---|
| `aggregate_results.sh` | Combines `all_results_multicycle.csv` files from many experiment dirs (default glob `runs/7QPX*/`) into a single CSV. With `--with-metrics`, runs `compute-final-metrics` per experiment first. **Note:** the orchestrator's `aggregate` subcommand and this shell script have confusingly overlapping names but do different things — orchestrator `aggregate` produces the per-experiment CSV; this shell script combines those across experiments. |

### 2.5 Diagnostic experiments (sequence-walk waypoints)

These were built during the notes-2 diagnostic session investigating why
6G10 converged but 7QPX did not, by walking the receptor/effector
sequence between the two complexes:

| File | Role |
|---|---|
| `generate_waypoints.py` | Emits 8 7QPX-anchored + 2 6G10-anchored chimera waypoints, plus an auto-generated submission script per anchor. |
| `generate_waypoints_7QPX_chainB.py` | Regeneration of the 7QPX-anchored waypoints at correct 73-aa chain-B length, after the chain-assignment bug was discovered. |
| `combine_waypoints.py` | Pure-stdlib post-processor. Combines per-waypoint `all_results_multicycle.csv` files into one combined CSV plus a one-row-per-waypoint summary. |
| `submit_steering_experiments.sh` | Multi-experiment wrapper that runs the same target through many steering-mode/parameter combinations (`strong_default`, `strong_few`, `mild_default`, `mild_pool`, `conservative`, `alanine_scan`, `alanine_pairs`) in one go. |

### 2.6 Negative-steering session notes

- `notes1.md` — Project overview, pipeline state at v5, file
  inventory.
- `notes2.md` — 7QPX diagnostic session: sequence-override mechanism
  (`--receptor-fasta` / `--effector-fasta` / `--skip-steering`),
  effector-override propagation, Singularity bind-path fix, waypoint
  generators.
- `notes3.md` — Production criterion lock-in (`compare_versions.py`
  refactor, dropped `tJ - wJ` constraint), skip-steering plumbing,
  aggregate retrofit, sidecar PDB locator fix. Pipeline declared
  production-ready for native benchmark complexes.
- `notes4.md` — RFDiffusion + ProteinMPNN extension session:
  `--true-interface-indices` / `--true-interface-indices-file` flags,
  the `split_sequence` chain-order bug fix in `pipeline_correct_sequences.py`,
  end-to-end Pikp1 design_1 run. Conclusion: workflow validated on
  RFDiffusion inputs; ipSAE confirmed broken on RFDiffusion data;
  HADDOCK identified as the upstream blocker for the receptor-resurfacing
  test case.
- `notes5.md` — Pikp1 design_3 single-design run with native ground
  truth (no HADDOCK confound). Discovered that high-confidence steered
  designs were contaminated by mutations contacting the effector
  inside the design region, motivating the reversion pass.
- `notes6.md` — Reversion-validation pass. Built the full 6-stage
  reversion chain (`iterate-collect-prefilter` → `build-contaminated` →
  `plan-reversions` → `predict-reversion-one` → `harvest-reversions`
  → `iterate-collect-finalize`), with the protected-set definition
  `design_region ∪ true_interface`. For `pose_holds` designs, the
  steered receptor.fasta is destructively overwritten with the
  reverted sequence so cycle N+1 starts from the validated sequence.
  Documents bugs found post-implementation (cycle-0 label parsing,
  reverted confidence not forwarded to aggregate CSV, `--n-cycles 1`
  silently skipping reversion).

### 2.7 Example output

| File | Role |
|---|---|
| `all_results_multicycle_with_metrics.csv` | Example of the final per-experiment output. One row per Boltz prediction across the multi-cycle tree, with structural RMSDs, confidence metrics, contact analysis, and (for designs that went through the reversion pass) reverted-prediction confidence and verdict columns. |

---

## 3. How the two pipelines connect (the integration goal)

The current `boltz2.nf` stage in the resurfacing pipeline runs each
MPNN design through Boltz2 once and ranks by ipSAE. This is
inadequate for two reasons surfaced in `notes5.md` and `notes6.md`:

1. **ipSAE is uniformly near-zero on RFDiffusion designs**, so it has
   no discriminative power.
2. **Mutations introduced by ProteinMPNN can directly contact the
   effector in Boltz's predicted pose**, meaning the predicted binding
   is conditioned on residues that may behave differently in the
   wet-lab construct.

The negative-steering workflow addresses both:

1. It uses **structural agreement with the design's intended binding
   site** (`receptor_aligned_effector_rmsd` vs the RFDiffusion design
   PDB) as the primary ranker, not confidence.
2. It runs the **reversion pass**, which reverts any mutation
   inside the protected set that contacts the effector and re-predicts.
   Designs whose pose holds without the load-bearing mutations
   (`pose_holds` verdict) are the ones safe to promote to wet-lab
   triage.

The intended end state is for `boltz2.nf` to be replaced with a
Nextflow port of the negative-steering chain (`plan` →
`predict array` → `prefilter` → `build-contaminated` →
`plan-reversions` → `predict-reversion array` → `harvest-reversions` →
`finalize`), iterated over `--n-cycles` cycles, with
`compute-final-metrics` and `compare_versions.py` producing the final
ranked output. Per `notes6.md` open work item 4, this Nextflow
integration is deferred until the negative-steering codebase is
restructured.

---

## 4. File-to-pipeline assignment quick reference

### Receptor-resurfacing pipeline (Nextflow)
**Modules:** `main.nf`, `aggregate.nf`, `boltz2.nf`, `haddock.nf`,
`msa.nf`, `preprocessing.nf`, `proteinmpnn.nf`, `rfdiffusion.nf`,
`rosetta_filtering.nf`
**Config / submission:** `nextflow.config`, `params_example.yml`,
`run_pipeline_slurm.sh`
**Test harnesses:** `test_boltz2.nf`, `test_haddock.nf`,
`test_proteinmpnn.nf`, `test_rfdiffusion.nf`,
`test_rosetta_filtering.nf` and their `run_test_*_slurm.sh` wrappers
**Backing scripts (bin/):** `boltz2_filter.py`, `boltz2_plots.py`,
`boltz2_prepare.py`, `boltz2_verify_binding.py`, `build_contigs.py`,
`collect_haddock3_dock.py`, `contig_utils.py`, `extract_hotspots.py`,
`haddock_utils.py`, `haddock3_plots.py`, `haddock3_prepare.py`,
`mpnn_cluster_sequences.py`, `mpnn_design_region_score.py`,
`mpnn_plots.py`, `mpnn_select_top.py`, `mpnn_sequence_qc.py`,
`pipeline_correct_sequences.py`, `rfdiffusion_contigs.py`,
`rfdiffusion_filter.py`, `rfdiffusion_plots.py`,
`rosetta_filter_collect.py`, `rosetta_filter_plots.py`
**Notes:** `pipeline_notes1.md`, `pipeline_notes2.md`

### Negative-steering pipeline (standalone)
**Core:** `boltz2_negative_steering.py`, `boltz2_iterate_steering.py`,
`reversion.py` (referenced in notes6 — may not be in this snapshot)
**Submission:** `submit_boltz2_negative_steering.sh`,
`submit_boltz2_iterate_steering.sh`, `submit_steering_experiments.sh`
**Metrics / post-processing:** `compute_metrics.py`,
`derive_design_region.py`, `compare_versions.py` (referenced in
notes — may not be in this snapshot), `aggregate_results.sh`
**Diagnostic experiments:** `generate_waypoints.py`,
`generate_waypoints_7QPX_chainB.py`, `combine_waypoints.py`
**Example output:** `all_results_multicycle_with_metrics.csv`
**Notes:** `notes1.md`, `notes2.md`, `notes3.md`, `notes4.md`,
`notes5.md`, `notes6.md`

### Cross-pipeline dependency
`pipeline_correct_sequences.py` lives in the **resurfacing pipeline**
(`bin/` for `proteinmpnn.nf`'s `SEQUENCE_CORRECTION` process) but was
**patched in a negative-steering session** (notes 4) because its
`split_sequence` bug was producing systematically wrong receptor
sequences that broke negative-steering's RFDiffusion-design test runs.
This is the only file in the project whose ownership crosses the two
codebases.
