# Reference outputs — `example_output_files/`

A trimmed copy of one real pipeline run's `results/` directory, used as
the reference data for the Phase 2 characterization-test suite. The
directory layout mirrors production exactly so tests can compare
fresh-run files to references by path symmetry: a fresh run's
`results/<area>/<file>` is compared directly against
`tests/full_test_run/example_output_files/<area>/<file>`.

**These are not certified-correct outputs.** They are the artefacts of
the supervisor-demo run and almost certainly contain undetected logical
errors (the inverted `scaffold` semantics, the `merge_orthogonal_metrics.py`
divergence, the suspected `--n-cycles 1` reversion-skip bug — see plan
§1.1 for the full caveat). The characterization tests assert
**stability**, not **correctness**: any fresh-run divergence from these
references is *visible*, but a green run does not certify the references
themselves. A refactor that fixes a latent bug should be expected to
*fail* these tests; that failure is signal, and the reference must then
be updated via the procedure in §6 below.

**Source:** 2026-04-29 supervisor-demo run (`params_test.yml`,
32 RFDiffusion designs × 4 MPNN sequences = 128 cohort sequences plus
2 negative controls). The full run was substantially trimmed; the
remaining files are documented per area below.

---

## 1. Per-area contents

### `preprocessing/`
Outputs of the preprocessing stage (`modules/preprocessing.nf`).
`sequences.json` is the receptor + effector sequence/length/resnum blob
emitted by the `EXTRACT_SEQUENCES` inline-Python process.
`processed_contigs.txt` is the resolved contig string written by
`RESOLVE_CONTIGS` → `bin/rfdiffusion_contigs.py`. The
`WRITE_DUMMY_MAPPING` placeholder JSON is intentionally not pinned (plan
§2.1).

### `rfdiffusion/`
Outputs of the backbone-generation stage
(`modules/rfdiffusion.nf`, `bin/rfdiffusion_filter.py`).
`rfdiffusion_metrics.json` holds per-design metrics (including
`design_region_coords`); `filter_summary.json` holds the cohort-level
filter counts; `passing_designs.txt` lists the design IDs that survived.
The `passing/`, `split/`, and `traj/` subdirectories are present but
empty — their PDB and trajectory contents were excluded (§5).

### `rosetta_filtering/`
Outputs of the pre-MPNN Rosetta-Sc filter stage
(`modules/rosetta_filtering.nf`, `bin/rosetta_filter_collect.py`).
`rosetta_filter_metrics.json` is the per-design Sc-and-friends blob;
`rosetta_filter_summary.json` is the cohort summary
(`n_total`/`n_passing`/`n_filtered`/`sc_threshold`);
`rosetta_passing_designs.txt` lists the surviving design IDs. The
`passing/` subdirectory is present but empty (PDBs excluded).

### `mpnn/`
ProteinMPNN per-design outputs (`modules/proteinmpnn.nf`). One
`design_<N>/` directory per RFDiffusion-passing design (32 in this run).
Each contains `fixed_positions_design_<N>.jsonl` (the motif/anchor
specification emitted by `MPNN_FIXED_POSITIONS`) and
`design_<N>_mpnn/seqs/design_<N>.fa` (the raw MPNN per-region sequences
emitted by the `PROTEINMPNN` invocation). The companion
`design_<N>_mpnn/scores/` directories are present but empty — the
per-design `.npz` score arrays were excluded.

### `sequences/`
Outputs of the post-MPNN sequence chain
(`modules/proteinmpnn.nf`, downstream of the MPNN run).
`mpnn_corrected.fasta` is the chimeric-receptor output of
`SEQUENCE_CORRECTION` → `bin/pipeline_correct_sequences.py`.
`qc_metadata.csv` / `qc_report.txt` / `qc_fastas/` (124 per-sequence
FASTAs) are produced by `SEQUENCE_QC` → `bin/mpnn_sequence_qc.py`.
`scored_metadata.csv` is from `MPNN_DESIGN_REGION_SCORE` →
`bin/mpnn_design_region_score.py`. `mpnn_cluster_counts.csv` is from
`MPNN_CLUSTER` → `bin/mpnn_cluster_sequences.py`. `top_metadata.csv` and
`top_fastas/` (32 selected sequences) are from `MPNN_SELECT_TOP` →
`bin/mpnn_select_top.py`. `sequence_metadata.csv` is the joined per-row
metadata table that flows through this chain.

### `negative_steering/`
Outputs of the negative-steering validation stage
(`modules/negative_steering.nf`, `modules/negsteer_controls.nf`).
At the top level: `cross_sequence_summary.csv` (cohort join produced by
`NEGSTEER_CROSS_SEQUENCE` → `bin/cross_sequence_summary.py`) and
`cross_sequence_summary_with_interface_metrics.csv` (the same table with
iRMSD / fnat / DockQ / 15 Å iPSAE / weighted-Jaccard appended by
`NEGSTEER_INTERFACE_METRICS` → `bin/compute_interface_metrics.py`).
`indices/` holds the 64 per-design index files (32 ×
`design_<N>_design_region.txt` + 32 × `design_<N>_true_interface.txt`)
emitted by `NEGSTEER_DERIVE_INDICES` → `bin/derive_design_region.py` /
`bin/derive_true_interface.py`. `controls_inputs/` holds the two
control-input index files (`input_design_region.txt`,
`input_true_interface.txt`) emitted by `DERIVE_INPUT_INDICES` →
`bin/derive_input_design_region.py`. `runs/` holds the per-sequence
`NEGSTEER_RUN_ONE` outputs (orchestrated by
`bin/negative_steering_run_one.sh` and the
`bin/boltz2_negative_steering.py` / `bin/boltz2_iterate_steering.py` /
`bin/reversion.py` chain) — only the three sequences listed in §2 are
pinned.

### `orthogonal_metrics/`
Outputs of the orthogonal-metrics validation stage
(`modules/negsteer_orthogonal_metrics.nf` and the four per-stream
modules). At the top level: `survivor_manifest.csv` (the per-survivor
fan-out from `EXTRACT_SURVIVOR_MANIFEST` →
`bin/extract_survivor_manifest.py`) and
`survivors_with_orthogonal_metrics.csv` (the merged AF3 / biophysical /
Rosetta table from `NEGSTEER_ORTHOGONAL_METRICS` →
`bin/merge_orthogonal_metrics.py`). Three per-stream subdirectories,
each with one entry per survivor (34 each in this run): `af3_nomsa/`
holds `<seq>/af3_nomsa_summary.csv` + `<seq>/input.json` from
`AF3_PARSE_OUTPUT` → `bin/parse_af3_output.py`; `biophysical/` holds
`<seq>/biophysical_summary.csv` from `NEGSTEER_BIOPHYSICAL_METRICS` →
`bin/run_biophysical_metrics.py`; `rosetta/` holds
`<seq>/rosetta_summary.csv` from `NEGSTEER_ROSETTA_METRICS` →
`bin/run_rosetta_metrics.py`.

### `plots/`
The 26 PNGs of the cohort run, produced by the six plot processes:
`rfdiff_*.png` (×5, `RFDIFFUSION_PLOTS` → `bin/rfdiffusion_plots.py`),
`rosetta_sc_histogram.png` (×1, `ROSETTA_FILTER_PLOTS` →
`bin/rosetta_filter_plots.py`), `mpnn_*.png` (×4, `MPNN_PLOTS` →
`bin/mpnn_plots.py`), `negsteer_*.png` cohort + within-sequence (×12
total, `NEGSTEER_PLOTS` / `NEGSTEER_WITHIN_SEQUENCE_PLOTS` →
`bin/negsteer_plots.py` / `bin/negsteer_within_sequence_plots.py`), and
`orthogonal_*.png` (×4, `ORTHOG_PLOTS` →
`bin/orthogonal_metrics_plots.py`). The 45-byte
`pipeline_plots.tar.gz` is the packaged-plots stub published by the
final tar step.

---

## 2. Sequence trio in `negative_steering/runs/`

Only three of the run's 34 per-sequence directories are pinned, chosen
to exercise the three meaningfully-distinct paths through the
negative-steering core. Every entry shares the same per-cycle skeleton:
`cycle_0/{plan.json, kickoff_distances.json, prefilter.json,
passing.json, contaminated.json, reversion_plan.json, reversion_results.json,
reversion_results_per_seed.json, steered_results.csv,
steered_results_aggregate.csv, summary.txt, true_interface_residues.txt,
wrong_interface_residues.txt, effector_template.cif, initial/, initial_s1/,
initial_s2/, steered/}` plus a sibling `inputs/` (FASTAs) and `logs/`
(empty stub). What differs is which JSONs are populated and whether
`reversions/` and `contamination_scratch/` exist.

### `input_control_polyA/` — control, cold-start-only path
The poly-alanine negative control. Cold-start Boltz runs but the
steering phase produces no effective contamination
(`contaminated.json` records `n_checked=0`, `n_contaminated=0`), so no
reversion is triggered. `cycle_0/steered/` holds only three entries
(`design_00_s0`, `design_00_s1`, `design_00_s2` — one design × three
seeds) rather than the ~60 of a real sequence; `reversion_plan.json` is
the 287-byte degenerate stub; `reversion_results.json` and
`reversion_results_per_seed.json` are the 2-byte `{}` empty case; there
is no `reversions/` subdirectory and no `contamination_scratch/`. The
sibling `inputs/` directory holds the receptor and effector FASTAs for
both controls (`control_polyA_receptor.fasta`,
`control_scrambled_receptor.fasta`, `effector.fasta`,
`source_effector.fasta`).

### `design_0_seq_0/` — empty reversion (`n_contaminated=0`)
A real MPNN sequence whose steering pass ran fully (60 steered designs
in `cycle_0/steered/`, populated `plan.json` / `kickoff_distances.json`
/ `passing.json` / `prefilter.json`) but for which the
`build-contaminated` step found no steering mutations inside the
protected set (`contaminated.json` records `n_contaminated=0`). The
reversion-related JSONs (`reversion_plan.json`,
`reversion_results.json`, `reversion_results_per_seed.json`) are
therefore the same empty stubs as the control, and there is no
`reversions/` subdirectory. `contamination_scratch/` exists but is
empty. This pins the "steering ran, reversion correctly skipped" path.

### `design_13_seq_2/` — populated reversion (`n_contaminated=33`)
The only sequence in the pinned set whose reversion pass actually fired.
`contaminated.json` is 37 KB and records `n_contaminated=33`;
`reversion_plan.json` (90 KB), `reversion_results.json` (54 KB), and
`reversion_results_per_seed.json` (105 KB) are all populated.
`cycle_0/reversions/` holds 66 reverted-design subdirectories named
`rev_design_<NN>_s<seed>_s<seed>/`, each containing `input.yaml`,
`receptor.fasta`, `metrics.csv`, `reversion_metadata.json`, and the
`boltz_results_input/` Boltz workdir stub.
`contamination_scratch/` is populated with per-design contamination
CSVs. This is the only entry in the reference set that exercises the
full reversion harvest path.

---

## 3. Exclusions

Deliberately removed during the subtractive rebuild:

- **All PDB files (`*.pdb`).** The cohort run produced thousands of
  per-design and per-prediction PDBs; none are pinned. Tests compare
  CSV/JSON metric outputs only. Empty `passing/`, `split/`, and `traj/`
  directories are left in place to preserve the production path layout.
- **All NPZ arrays (`*.npz`).** Boltz model-internal logits and
  embeddings are excluded (e.g. the empty `mpnn/.../scores/`
  directories, and inside Boltz run trees).
- **All MSA files (`*.a3m`).** Single-sequence A3M inputs and any MSA
  search outputs are excluded.
- **`predictions/` and `msa/` subdirectories** within Boltz run trees
  (the `boltz_results_input/processed/` and `lightning_logs/`
  internals). These are predictor scratch space and have no place in a
  reference set.
- **31 of 34 per-sequence directories in `negative_steering/runs/`.**
  Only the trio above is pinned.

---

## 4. Reference set change log

### 2026-05-01 — Subtractive rebuild
Replaced the prior flat-ish reference set with a structure mirroring
production. Source as named in the header (2026-04-29 supervisor-demo
run). Excluded as listed in §3. Per-sequence trio chosen to exercise
the three meaningfully-distinct paths through the negative-steering
core (control / empty-reversion / populated-reversion).

---

## 5. Updating the reference set

When a Phase 3+ change legitimately alters pipeline behavior (e.g. the
`merge_orthogonal_metrics.py` AF3-demotion fix, the `scaffold`/`motif`
rename, or any other deliberate behavior change), the reference set
must be updated by hand. The full procedure — pre-update green check,
fresh HPC run, expected-failure verification, file-by-file copy,
README change-log entry, post-update green check, joint commit — is
specified in **plan §8.1** (`notes/inventory/11_phase_2_plan.md`).
The procedure is intentionally manual so that each reference update is
a deliberate, reviewable decision rather than an automated overwrite.
