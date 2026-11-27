---
status: AUTHORITATIVE
scope: pipeline-overview
last_verified: 2026-05-19
branch_at_verification: phase4-impl
---

# Pipeline Overview

This document is the single canonical reference for what the receptor
resurfacing pipeline does, in order, with maximum technical detail.  It
supersedes the high-level summaries in `README.md`,
`notes/inventory/04_functional_categorization.md`, and the scattered
narrative comments in `main.nf`.

When this document and code disagree, code wins.  Open an issue and
update this document.  All file:line citations are against branch
`phase4-impl` at the verification date above.

The single most error-prone area is **negative steering** (§4).  Read
§4 in full before touching anything in `bin/boltz2_iterate_steering.py`,
`bin/boltz2_negative_steering.py`, `bin/reversion.py`,
`bin/stage_result.py`, or `bin/cross_sequence_summary.py`.

---

## 0. How to read this document

- §1 — what the pipeline produces.
- §2 — input modes.
- §3 — stage-by-stage flow.  One subsection per stage, with the
  Nextflow process name, the bin script(s) it calls, the inputs/outputs,
  and the params it reads.
- §4 — negative steering, in maximal detail.  This is the longest
  section because it is the most complicated and most-confused part of
  the pipeline.
- §5 — orthogonal validators (AF3 / biophysical / Rosetta).
- §6 — cohort table format.
- §7 — conventions (chains, frames, indices, contigs, thresholds).
- §8 — file:line map for every behavior described in §4.
- §9 — known divergences between this doc, `notes/design_audit.md`,
  and `notes/phase4_architecture_spec.md`.
- §10 — glossary cross-reference.

---

## 1. What the pipeline produces

A ranked cohort table of candidate binder sequences (`receptor`
sequences) for a fixed target (`effector`).  Each row is one
MPNN-designed receptor sequence on one RFDiffusion-generated backbone.
Columns include:

- the designed receptor sequence (full + design-region slice),
- a Boltz-2 tier (A / B / C / none) and per-seed verdicts,
- Boltz-2 confidence metrics (pLDDT, ipTM, ipSAE, PAE),
- orthogonal validator results (AF3 no-MSA `ra_eff`, Rosetta `Sc`,
  BSA, ΔΔG, interface H-bonds),
- per-row interface metrics (iRMSD, fnat, DockQ, weighted Jaccard),
- a composite score and within-cohort ranks,
- the path to the representative predicted complex PDB
  (`rep_canonical_pdb`).

The wet-lab triage list is `passing_summary.csv` (a per-MPNN-sequence
subset) joined into `cross_sequence_summary.csv` (the one row per MPNN
sequence cohort table).

---

## 2. Input modes

The pipeline accepts the receptor–effector pair in two mutually
exclusive ways:

- **Branch A (HADDOCK)** — supply two **monomer** PDBs
  (`params.receptor_input` + `params.effector_input`).  HADDOCK3 docks
  them.  Internal stages rename the resulting chains to `A` = receptor,
  `B` = effector regardless of input chain IDs.
  `params.receptor_chain` / `params.effector_chain` describe the
  *input* chain letters so the HADDOCK prepare step can extract the
  right atoms.

- **Branch B (pre-docked)** — supply one **complex** PDB
  (`params.pdb_file`) already containing receptor and effector in their
  intended pose.  This is the canonical entry point in current
  production runs.  HADDOCK is acknowledged dead-but-preserved code
  pending the Session 7 restructure (see
  `notes/haddock_protocol_experiments.md`).

In either case, **post-RFDiffusion all stages assume chain A =
receptor, chain B = effector**, hardcoded in
`bin/rfdiffusion_filter.py:write_split_pdb`.

---

## 3. Stage-by-stage flow

Order of execution (top-level orchestration in `main.nf:351-1127`):

### 3.0 Parameter validation

`main.nf:358-371` launches `bin/validate_params.py` against a JSON
dump of `params`.  `PARAM_SPECS` enumerates ~66 expected parameters
with type, range, and dependency constraints.  Non-zero exit aborts
the workflow before any process schedules.  Add new params here.
Known gaps (post-Session 6, deferred): infrastructure params like
`af3_db_v3`, `boltz2_container`, etc., do not yet have specs.

### 3.1 Preprocessing

| Process | Module | Script | Purpose |
| --- | --- | --- | --- |
| `EXTRACT_SEQUENCES` | `preprocessing.nf` | inline + `extract_sequences` (in `rfdiffusion_contigs.py`) | Extract receptor + effector chain sequences from the input PDB; record `receptor_start_pdb` (first residue number). |
| `RESOLVE_CONTIGS` | `preprocessing.nf` | `rfdiffusion_contigs.py` | Resolve user contig notation to RFDiffusion's native format (range expansion etc.). |
| `WRITE_DUMMY_MAPPING_REC/_EFF` | `preprocessing.nf` | inline | Branch A only: emit 2-byte `{}` placeholder mapping files so the BUILD_CONTIGS step downstream has a concrete file input. |

### 3.2 Branch A — HADDOCK3 docking (skipped on Branch B)

| Process | Module | Script | Purpose |
| --- | --- | --- | --- |
| `HADDOCK3_PREPARE` | `haddock.nf` | `haddock3_prepare.py` | Write receptor/effector PDBs renamed to A/B; emit `ambig_restraints.tbl` (active residues mode) and/or `unambig_restraints.tbl` (contact-pair mode); `restraints_summary.json`. |
| `HADDOCK3_DOCK` | `haddock.nf` | inline (HADDOCK3) + `collect_haddock3_dock.py` | Run HADDOCK3 (topoaa → rigidbody → seletop → flexref → emref → clustfcc → seletopclusts → caprieval); post-process cluster outputs.  Sampling = `params.haddock_sampling` (prod 10000, validator min 100), `seletop = params.haddock_seletop = 20`, `min_population = params.haddock_min_cluster_size` (prod 4). |
| `HADDOCK_CLUSTER_METRICS` | `haddock.nf` | `haddock_cluster_metrics.py` | Per-cluster: BSA, COM distance, AIR satisfaction count/total, pair_contact_fraction, clashes inside vs outside design region (2.0 Å heavy-atom cutoff). |
| `HADDOCK_CLUSTER_SC` | `haddock.nf` | `haddock_cluster_sc.py` | Per-cluster Rosetta shape complementarity. |
| `HADDOCK3_PLOTS` | `haddock.nf` | `haddock3_plots.py` | Diagnostic plots. |
| `SELECT_HADDOCK_CLUSTER` | `haddock.nf` | `select_haddock_cluster.py` | Auto-pick by lexicographic key `(-pair_contact_fraction, -BSA)` OR honor `params.haddock_chosen_cluster: N` user override.  Stop gate: `params.stop_after_haddock: true` halts the run for manual inspection. |
| `BUILD_CONTIGS` | `haddock.nf` | `build_contigs.py` | Re-derive the contig string in the chosen cluster's coordinate frame; emit `updated_params.json`. |

HADDOCK's role here is **geometric placement, not energetics**.
Cluster ranking is by `(pair_contact_fraction, BSA)` not by HADDOCK
score.  See `notes/haddock_protocol_experiments.md` for the rationale.
The Session 7 restructure removed the old auto-derived hotspot
extraction (the `EXTRACT_HOTSPOTS` process).  `params.hotspot` is now
user-optional.

### 3.3 RFDiffusion backbone generation

| Process | Module | Script | Purpose |
| --- | --- | --- | --- |
| `RFDIFFUSION` | `rfdiffusion.nf` | RFDiffusion inline | Generate `params.num_designs` candidate receptor backbones around the fixed effector pose.  Uses `Complex_beta_ckpt.pt`; `add_potential: true`; `rfdiff_iterations` (prod 100, test 50); `diffusion_samples = 5` per seed. |
| `RFDIFFUSION_FILTER` | `rfdiffusion.nf` | `rfdiffusion_filter.py` | Filter by interface contact (Cα 8 Å cutoff vs effector hotspots).  **Hardcodes chain A=receptor / B=effector** in `write_split_pdb`.  Currently maximally permissive (`hotspot=""`, `min_hotspot_frac=0.0`) — exploratory phase. |
| `RFDIFFUSION_PLOTS` | `rfdiffusion.nf` | `rfdiffusion_plots.py` | Diagnostic plots. |

RFDiffusion outputs Cα-only PDBs.  Geometry metrics
(`motif_rmsd`, `endpoint_distance`, COM displacement) are
**diagnostic only — not gating**.

Contig format is the canonical RFDiffusion string.  Segments are
`/`-separated within a chain and space-separated between chains, with
no commas.  Fixed segments are chain-prefixed (`A1-10`), and de novo
segments are bare lengths or length ranges (`5-7`).  Details in §7.4.

### 3.4 Rosetta pre-MPNN filtering

| Process | Module | Script | Purpose |
| --- | --- | --- | --- |
| `ROSETTA_SC` | `rosetta_filtering.nf` | Rosetta InterfaceAnalyzer inline | Per-design Sc, dG_separated, dSASA_int on the polyvaline-backbone outputs of RFDiffusion. |
| `ROSETTA_FILTER` | `rosetta_filtering.nf` | `rosetta_filter_collect.py` | Apply `params.sc_threshold` (default 0.5).  dG_separated is computed but **not used for filtering** at this stage — values are meaningless on polyvaline. |
| `ROSETTA_FILTER_PLOTS` | `rosetta_filtering.nf` | `rosetta_filter_plots.py` | Plots. |

Stop gate: `params.stop_after_rosetta: true` halts after this stage
for manual ChimeraX inspection.

### 3.5 ProteinMPNN sequence design

| Process | Module | Script | Purpose |
| --- | --- | --- | --- |
| `MPNN_FIXED_POSITIONS` | `proteinmpnn.nf` | `pipeline_correct_sequences.py` | Emit JSONL of fixed positions per design (where to keep native, where to design). |
| `PROTEINMPNN` | `proteinmpnn.nf` | ProteinMPNN inline | `num_seqs` sequences per surviving design.  Production: `num_seqs=2`, `mpnn_sampling_temp=0.25`, `rm_aa="C"`. |
| `SEQUENCE_CORRECTION` | `proteinmpnn.nf` | `pipeline_correct_sequences.py` | Overlay each MPNN-designed slice onto the native sequence at fixed positions to produce the "corrected receptor" — the full sequence Boltz will predict. |
| `SEQUENCE_QC` | `proteinmpnn.nf` | `mpnn_sequence_qc.py` | Filter: `max_poly_x=5` (default), single-AA fraction <30%, optional pct-identity bounds. |
| `MPNN_DESIGN_REGION_SCORE` | `proteinmpnn.nf` | `mpnn_design_region_score.py` | Re-score designed positions only. |
| `MPNN_CLUSTER` | `proteinmpnn.nf` | `mpnn_cluster_sequences.py` | MMseqs2 clusters at multiple identity thresholds. |
| `MPNN_PLOTS` | `proteinmpnn.nf` | `mpnn_plots.py` | Plots. |
| `MPNN_SELECT_TOP` | `proteinmpnn.nf` | `mpnn_select_top.py` | (Optional, gated on `mpnn_top_n>0`) Pick top-N MPNN sequences **globally** by score (lower = better).  Multiple sequences per design are allowed if their MPNN scores rank top globally. |

Production defaults (Session 3 of the grill-me):
`num_designs=64`, `num_seqs=2`, `mpnn_top_n=128`,
`negsteer_n_designs=3`.  Rationale: more backbone diversity + fewer
steering attempts per sequence produces more tier-A hits than forcing
many steering attempts on fewer backbones.

### 3.6 Negative steering

See §4 in full.  Summary: every selected MPNN sequence (plus two
controls) runs through one GPU SLURM job (`NEGSTEER_RUN_ONE`).  That
job performs cold-start prediction → steered design generation →
steered prediction → contamination check → optional reversion →
per-seed verdict assignment → per-sequence aggregation.  Then
`extract_passing.py` emits `passing_summary.csv`, and a single CPU
`NEGSTEER_CROSS_SEQUENCE` task joins all per-sequence summaries into
`cross_sequence_summary.csv`.

### 3.7 Cross-sequence aggregation

| Process | Module | Script | Purpose |
| --- | --- | --- | --- |
| `NEGSTEER_DERIVE_INDICES` | `negative_steering.nf` | `derive_design_region.py` + `derive_true_interface.py` | Per-design `_design_region.txt` (1-based) and `_true_interface.txt` (0-based).  `--cross-check-design-region` validates consistency. |
| `NEGSTEER_RUN_ONE` | `negative_steering.nf` | `negative_steering_run_one.sh` | One GPU job per MPNN sequence; runs the 7-stage negsteer chain (§4.2). |
| `DERIVE_INPUT_INDICES` | `negsteer_controls.nf` | `derive_input_design_region.py` | Same for the pre-RFDiffusion input PDB (used by controls). |
| `NEGSTEER_CONTROLS` | `negsteer_controls.nf` | `build_control_sequences.py` + `negative_steering_run_one.sh` | Scrambled + polyA controls through the same negsteer loop. |
| `NEGSTEER_CROSS_SEQUENCE` | `negative_steering.nf` | `cross_summary_v2.py` → `cross_sequence_summary.aggregate` | Symlink all per-sequence workdirs under `aggregator/`; emit `cross_sequence_summary.csv`. |
| `NEGSTEER_PLOTS` | `negative_steering.nf` | `negsteer_plots.py` | 7 cohort plots. |
| `NEGSTEER_WITHIN_SEQUENCE_PLOTS` | `negative_steering.nf` | `negsteer_within_sequence_plots.py` | 5 per-seed plots. |
| `NEGSTEER_INTERFACE_METRICS` | `negsteer_interface_metrics.nf` | `compute_interface_metrics.py` | Append iRMSD / fnat / DockQ / ipSAE_15 / intact_core / weighted_jaccard columns to `cross_sequence_summary.csv`. |

### 3.8 Orthogonal validators

See §5.  Three parallel streams run per survivor, plus a merge.
`params.orthogonal_tier_filter` selects the survivors.  `'all'` runs
the cascade on every steered row including tier-none, and `'abc'`
restricts to A/B/C.

---

## 4. Negative steering — full reference

### 4.1 Purpose

Given an MPNN sequence (one chosen receptor sequence for one
RFDiffusion design), test whether Boltz-2 predicts that the receptor
binds the effector in the intended pose.

Boltz can fail in three structurally distinct ways:

1. **Wrong receptor fold** — receptor structure collapses, ra_eff
   unmeaningful.
2. **Wrong pose** — receptor folds but binds at the wrong surface
   (`ra_eff_vs_truth` is large).  *Negative steering* mutates the
   residues mediating the wrong pose to push Boltz toward the
   intended pose.
3. **Apparent right pose but with off-interface contacts driving it**
   — Boltz placed the effector on the intended interface, but only
   because steering mutations happen to contact the effector.  When
   we mutate those mutated positions back to wild type, the pose
   collapses.  *Reversion* is the second test.  Revert the
   contaminating mutations and predict again.  If the pose still
   holds, the sequence is a true positive.

### 4.2 The unit of work

One **`NEGSTEER_RUN_ONE` job** = one MPNN sequence (or one control).
Inside that single GPU SLURM job, the orchestrator
(`bin/negative_steering_run_one.sh`, ~449 LOC) runs seven stages in
order.  Cross-MPNN-sequence parallelism comes from Nextflow scheduling
many `NEGSTEER_RUN_ONE` jobs concurrently.  There is **no
shared-state cooperation between sequences inside negsteer**.

Within one job the seven stages are:

1. **Plan (GPU)** — cold-start prediction(s) + steering plan
   (`boltz2_negative_steering plan`, `negative_steering_run_one.sh:222-247`).
2. **Predict-one × N (GPU)** — one Boltz call per steered design slot
   (`boltz2_negative_steering predict-one`, `negative_steering_run_one.sh:282-303`).
3. **Collect (CPU)** — rank per-design Boltz outputs into
   `steered_results.csv` (`negative_steering_run_one.sh:308-311`).
4. **Reversion prep (CPU)** — `kickoff-distances` →
   `kickoff-prefilter` → `build-contaminated` (CL-3 gate, §4.6) →
   `plan-reversions` (writes `reversion_plan.json` of *unique*
   reverted sequences) (`negative_steering_run_one.sh:316-343`).
5. **Predict-reversion-one × M (GPU)** — one Boltz call per staged
   reversion (`negative_steering_run_one.sh:346-377`).
6. **Harvest + finalize (CPU)** — `harvest-reversions` computes
   reverted metrics and assigns the per-design reversion verdict.
   `kickoff-finalize --max-cycles 0` writes `cycle_0/passing.json`
   (`negative_steering_run_one.sh:380-396`).
7. **Aggregate (CPU)** — `aggregate` (per-seed CSV) →
   `compute-final-metrics --populate-all` (steered confidence on
   every non-dropped row) → `aggregate-per-sequence` →
   `extract_passing.py` (`negative_steering_run_one.sh:408-443`).
   `aggregate-per-sequence` collapses per-seed into
   per-unique-sequence rows and derives `outcome` + `n_pass`, and
   `extract_passing.py` writes `passing_summary.csv`.

### 4.3 Cycle numbering

Cycles are **0-indexed**.  In current production only **cycle 0**
runs.  The orchestrator structure supports multi-cycle propagation
(cycle-1+ host-side paths in
`bin/boltz2_iterate_steering.py:cmd_iterate_*`) but those cmds are
unused by Nextflow.  Multi-cycle is live roadmap, not dead code.

The phrasing "cycle 0" everywhere refers to:

- the *first and only* steering pass on a sequence,
- the cold-start prediction(s) of the wild-type MPNN sequence,
- the steered design pool generated from the cold-start observation,
- the reversion staged off those steered designs.

When you see `cycle_0/plan.json` or `cycle_0/passing.json` in
workdirs, all of that content belongs to the one production cycle.

`n_cycles=1` in production.  The historical `cmd_kickoff*` family
(~600 LOC in `boltz2_iterate_steering.py:2055+`) is live multi-cycle
roadmap, not dead.  The only path Nextflow exercises is
`kickoff-distances` → `kickoff-prefilter` → `kickoff-finalize
--max-cycles 0`.

### 4.4 Cold-start

The **cold-start prediction** is Boltz-2's prediction of the
**unmodified MPNN sequence** (no steering mutations applied) against
the effector.  `num_seeds` cold-start predictions are run per MPNN
sequence (production `num_seeds=3`, which must be odd, enforced by
`bin/validate_params.py:_check_num_seeds_odd`).

Cold-start is performed by `cmd_plan` in
`bin/boltz2_negative_steering.py:1818`, specifically the seed loop at
`:1936-2015`.

A cold-start prediction is a "pass" iff:

- the receptor is *intact* (independent receptor RMSD < `intact_threshold`, default 5.0 Å), AND
- the effector is *intact* (independent effector RMSD < same threshold), AND
- `ra_eff_vs_truth < structural_ra_eff_cutoff` (default 5.0 Å, from `--rmsd-threshold` / `PipelineInternalThresholds.steered_ra_eff`).

**Skip-steering** (`cold_start_all_clean`, `boltz2_negative_steering.py:2017-2035`):
if **ALL** `num_seeds` cold-start predictions pass, the entire
steering apparatus is skipped.  The plan writes
`skip_steering=true`, produces no steered design pool, and the
sequence proceeds directly to per-seed verdict assignment based on
cold-start.

If even **one** cold-start fails, **all** seeds — including the ones
that already passed cold-start — proceed through steering.  The
rationale is recorded in Session 1 Q32, design_audit.md:804-822.
Median confidence metrics require the full seed set, and partial seed
sets corrupt median aggregation.

### 4.5 Steering

When cold-start does not all-pass, the plan stage selects a pool of
positions to mutate:

```
candidate_pool = wrong_interface − protected_set
protected_set  = design_region ∪ true_interface
```

- `wrong_interface` = receptor residues contacting the effector in
  the cold-start prediction (5 Å heavy-atom contact, computed by
  `find_contact_residues_heavy` in
  `bin/boltz2_negative_steering.py`).  These are the "wrong" because
  they were not supposed to be touching the effector.
- `true_interface` = receptor residues contacting the effector in
  the *ground-truth* pose (the user-supplied complex), pre-computed
  by `derive_true_interface.py`.  Stored 0-based.
- `design_region` = the de novo-generated positions (RFDiffusion's
  new residues).  Stored 1-based.

`true_interface` and `design_region` are both **excluded** from the
mutable pool.  Mutating `true_interface` would defeat the purpose,
and `design_region` carries the MPNN sequence we want to evaluate.
The `protected_set` is the architectural guarantee, not a runtime
check.

The plan generates `params.negsteer_n_designs` (production 3)
distinct steered receptor sequences.  Each steered sequence has up
to `--max-mutations` substitutions from the candidate pool, drawn
from a per-residue substitution mode (`--mode` one of
`strong=WYFRKEDP | mild=DEKR | conservative | alanine`).  Each
steered sequence is written to
`<workdir>/cycle_0/design_<NN>/receptor.fasta` along with a
`mutations.tsv` log and a Boltz YAML input file.

Each steered design then receives `num_seeds` Boltz predictions via
`cmd_predict_one` (`bin/boltz2_negative_steering.py:2788`).  Each
prediction yields a PDB and a `result.json` containing receptor and
effector RMSDs, `ra_eff_vs_truth`, and the Jaccard interface
overlaps against the true and wrong interfaces.

### 4.6 The CL-3 reversion gate

After all steered predictions complete and pass the structural
filter (intact + `ra_eff < cutoff`), each steered design is checked
for **contamination**.  The question is whether any steered mutation
position is in direct contact (5 Å heavy-atom) with the effector in
the steered prediction.

The **CL-3 gating rule** — current canonical behavior, encoded in
two places:

```
reversion runs for a steered design iff
    n_correctly_placed > 0
    AND n_contaminated >= ceil(n_correctly_placed / 2)
```

where, for *this one steered design*:

- `n_correctly_placed` = number of seeds whose steered prediction was
  intact AND ra_eff under threshold,
- `n_contaminated` = number of those correctly-placed seeds whose
  contacts include any mutated position (intersected with
  `contamination_positions = design_region ∪ true_interface`).

The predicate is enforced at two sites:

- **Type-level** (`bin/stage_result.py:278-282`):

  ```python
  n_correct = self.n_correctly_placed(thresholds, truth)
  if n_correct == 0:
      return False  # nothing to revert
  n_contam = self.n_contaminated(thresholds, truth, contamination_positions)
  return n_contam >= math.ceil(n_correct / 2)
  ```

- **Live pipeline** (`bin/boltz2_iterate_steering.py:cmd_build_contaminated`,
  lines 1532-1568): groups `per_seed_evaluations` by `design_slot`,
  then per design:

  ```python
  n_correctly_placed = len(evals)
  n_contaminated = sum(1 for e in evals if e["is_contaminated"])
  majority_threshold = math.ceil(n_correctly_placed / 2)
  if n_contaminated < majority_threshold:
      continue  # gate did not pass — no reversion for this design
  ```

**This replaced the pre-Phase-4 per-(design, seed) rule** which
queued reversion any time any single seed showed contamination.  The
new rule reduces unnecessary reversion when the contamination is a
minority phenomenon.

Worked examples for `num_seeds = 3` per design:

| n_correctly_placed | n_contaminated | Reversion runs? | Final outcome (assuming reverted-OK) |
| --- | --- | --- | --- |
| 3 | 0 | No (clean_steered × 3) | tier A, outcome `no_reversion` |
| 3 | 1 | No — `1 < ceil(3/2)=2` | tier B, outcome `no_reversion` (one seed = `contaminated`, no pass-equivalent) |
| 3 | 2 | Yes | tier depends on reverted result |
| 3 | 3 | Yes | tier depends on reverted result |
| 2 | 1 | Yes (`1 >= ceil(2/2)=1`) | tier depends on reverted result |
| 0 | 0 | No | tier none (all seeds wrong-placement, `no_data`) |

The "n_correctly_placed=3, n_contaminated=1" case is the subtle one.
The rule says reversion does NOT run, but the contaminated seed
still does not count as `pose_holds`/`clean_steered` in `n_pass`.
It buckets as `contaminated_skipped` (CL-3 §92 tentative name) or,
in the live code, as a non-pass-equivalent verdict that reduces
`n_pass` below `n_seeds`.

### 4.7 Reversion

For each design whose CL-3 gate fired, a **reverted sequence** is
built by `bin/reversion.py:build_reverted_sequence` (line 147).  The
reverted sequence is the steered sequence with the contaminating
mutations substituted back to wild type.  Mutations that did NOT
contaminate are *kept*, and only the contacts-with-effector
positions are reverted.

`bin/reversion.py:write_reversion_plan` (line 304) **dedups by
reverted-sequence hash**.  Two steered designs whose contaminations
happen to yield the same reverted sequence are predicted **once**.
An `applies_to` list maps the reverted prediction back to every
steered design it serves.  This is encoded in `StageResult.
applies_to_designs` for the reversion stage type (Phase 4 §2.8).

The reverted Boltz YAMLs are written under
`<workdir>/reversions/<reverted_id>/` and predicted by
`cmd_predict_reversion_one`
(`bin/boltz2_iterate_steering.py:2649`).

`bin/reversion.py:harvest_reversion_results` (line 832) walks the
predicted reverted directories and computes metrics.  It then routes
each **design–seed pair** that triggered reversion to a reversion
verdict via `bin/reversion.py:classify_reversion_verdict` (line
1204):

- `pose_collapses` — reverted is not intact OR reverted_ra_eff >= cutoff (`:1245`).
- `new_contamination` — reverted is structurally OK AND the mutated
  positions in the reverted prediction contact the effector inside
  the gating set `design_region ∪ true_interface` (`:1267`).
- `pose_holds` — reverted is structurally OK AND no new contamination
  (`:1297`).
- `no_data` — input metrics missing (`:1304`).

The reverted prediction is the **final** prediction for that
(design, seed) pair.  Its PDB and its confidence metrics override
the steered ones in downstream aggregation.

### 4.8 Per-seed verdicts (full table)

`bin/stage_result.py:per_seed_verdicts` (line 117) defines the
vocabulary per stage.  Below is the *production-relevant* version
(stage = `reversion` when reversion ran, else `steered`):

| Final per-seed verdict | Path | Defining condition | Counts in n_pass? |
| --- | --- | --- | --- |
| `clean_steered` | steered (reversion correctly skipped) | intact + ra_eff < cutoff + no contamination on mutated positions | **YES** |
| `pose_holds` | reversion ran, reverted ok | reverted intact + reverted_ra_eff < cutoff + no `new_contamination` on mutated positions ∩ gating set | **YES** |
| `pose_collapses` | reversion ran, reverted bad | reverted not intact OR reverted_ra_eff ≥ cutoff | NO |
| `new_contamination` | reversion ran, reverted re-contaminated | reverted intact + reverted_ra_eff ok BUT mutated positions still contact effector inside gating set | NO |
| `contaminated_skipped` (tentative CL-3 §92) | steered with contamination but CL-3 gate did not fire | correctly placed + contaminated, but `n_contaminated < ceil(n_correctly_placed/2)` for this design | NO |
| `no_data` | steered structural fail OR missing inputs | intact false OR ra_eff ≥ cutoff OR inputs unparseable | NO |

In the live aggregator, the per-seed verdict is reconstructed by
`_per_seed_verdict_breakdown`
(`bin/boltz2_iterate_steering.py:3728`).  It reads each row's
`reversion_verdict` column when reversion ran.  Otherwise it applies
a heuristic for the cold-start-all-clean / clean-steered path at
`:3813-3835`.

**Important.**  The "wrong placement, no contamination" case lands as
`no_data` (not `clean_steered`).  Earlier audit-era language used
"clean_steered" loosely to mean "no contamination regardless of
placement".  The live code's definition is stricter and is the
canonical one.

### 4.9 Aggregate outcomes (full table)

`_classify_outcome` in `bin/boltz2_iterate_steering.py:3555-3725`
emits one of these labels into the `outcome` column on
`aggregated_results.csv`:

| Outcome | Defining condition | Assignment site |
| --- | --- | --- |
| `no_reversion` | no reversion ran for any seed in the group: either skip_steering (cold-start all-clean) OR steering ran but no seed needed reversion | `:4742` (else branch) |
| `pose_holds` | reversion ran; aggregated intact-majority=1; reverted_ra_eff median < cutoff; no majority contamination at gating positions; `n_pose_holds + n_clean_steered >= 1`; `non_failing >= max(n_pose_collapses, n_new_contamination)` | `:3678-3725` |
| `pose_collapses` | reverted not intact by majority, OR reverted_ra_eff median missing/unparseable, OR median ≥ cutoff; also a downgrade target when no seed individually passes | `:3622-3697` |
| `new_contamination` | reverted_mutated_contact_positions_majority contains a position inside the gating set; also downgrade when `n_new_contamination > n_pose_collapses` | `:3654-3717` |
| `singleton` | sequence_group missing (pre-multi-seed rows); pass-through with `n_pass = ""` | `:4764` |

**Critical:** `no_reversion` is **NOT** synonymous with cold-start
all-clean.  It is path-agnostic and fires in TWO disjoint cases:

(a) **Cold-start all-clean.** All seeds passed cold-start, so
    `skip_steering=true` and steering never ran.  Every per-seed
    verdict is `clean_steered`, and `n_pass == n_seeds` → tier A.

(b) **Steering ran but no reversion needed.** Cold-start failed for
    at least one seed, so steering ran on every seed, but for every
    steered seed the result was either:
    - structurally OK with zero contamination on mutated positions →
      `clean_steered` (counts toward `n_pass`), OR
    - wrong placement / unintact → `no_data` (does NOT count).

    `n_pass` reflects only the `clean_steered` count and may be
    anywhere from 0 to `n_seeds`.  A "wrong placement, no
    contamination" group has every seed = `no_data`, `n_pass = 0`,
    outcome = `no_reversion`, **tier none** — NOT tier A.

This is the case that the memory file
`project_no_reversion_semantics.md` repeatedly flags.

### 4.10 n_pass — path-agnostic

`bin/boltz2_iterate_steering.py:4748-4749`:

```python
agg["n_pass"] = (verdict_counts["pose_holds"]
                + verdict_counts["clean_steered"])
```

`n_pass` is *always* this sum, regardless of which path got each
seed there.  A seed contributes 1 iff its final per-seed verdict is
either `pose_holds` or `clean_steered`.  `pose_holds` means reversion
ran and the reverted pose held, and `clean_steered` means reversion
correctly did not run because the steered prediction had zero
contamination.  Anything else — including `no_data`,
`pose_collapses`, `new_contamination`, `contaminated_skipped` —
contributes 0.

Denominator is `n_seeds` (also written into the aggregated row).
`n_pass` can be strictly less than `n_seeds` even when
`outcome == "no_reversion"` (the case-(b) above).

### 4.11 Tiers (A/B/C/none)

`bin/cross_sequence_summary.py:_tier_for_row` (line 205), called from
the cross-sequence join after `passing_summary.csv` is produced:

```
Tier A: n_pass == n_seeds        (all seeds pass-equivalent)
Tier B: 1 < n_pass < n_seeds     (some but not all)
Tier C: n_pass == 1              (exactly one)
none:   n_pass == 0  OR insufficient data
```

For `n_seeds = 3`: A=3, B=2, C=1, none=0.

For `n_seeds = 5`: A=5, B=2-4, C=1, none=0.

For `n_seeds = 7`: A=7, B=2-6, C=1, none=0.

The function explicitly rejects the outdated
`outcome=='no_reversion' → tier A` shortcut at its docstring
(`bin/cross_sequence_summary.py:220-225`).  Tier is derived from
`n_pass` only.

`num_seeds` MUST be odd (`PipelineInternalThresholds` invariant +
`validate_params._check_num_seeds_odd`) for unambiguous median
aggregation of confidence metrics.

Cross-rank ordering is in `_assign_cross_ranks`
(`bin/cross_sequence_summary.py:553`).  The primary sort is by tier
(A before B before C).  Tier-none rows are excluded from
`cross_rank_by_composite` and `cross_rank_by_ra_eff` but still
included in the CSV.

### 4.12 Representative prediction

Two layers of selection produce the "representative prediction" the
cohort plot displays:

**Layer 1 — per unique sequence (canonical_pdb).**
`_canonical_seed` in
`bin/boltz2_iterate_steering.py:4654-4683` picks one seed per
aggregated group:

- Sort the group's seeds by `steered_ra_eff_vs_truth` ascending
  (NaN/missing sorted to end).
- Odd N: pick index `(N-1)//2` — the median.
- Even N: pick index `N//2` — the *upper-middle*, i.e. the worse of
  the two middle seeds.  Emits a `even_n=N_chose_worse_of_two_middle_seeds`
  warning.
- N=1: that seed.

The chosen seed's prediction PDB becomes `canonical_pdb`.

Rationale: the median represents the typical seed honestly.  If the
median is failing, the conclusion is that more seeds are failing than
passing.  Best-seed PDBs are still available in
`raw_per_seed_results.csv` for manual inspection but are not the
default.

**Layer 2 — per MPNN sequence in cross-summary
(rep_canonical_pdb).**
`_pick_representative` in
`bin/cross_sequence_summary.py:324-351` picks one row from
`passing_summary.csv` minimizing `(tier_order, within_sequence_rank)`
— i.e. preferring tier A → B → C, then within-tier preferring lower
composite rank.  Its `canonical_pdb` is rewritten by
`_rewrite_workdir_path_to_published` to
`<published_runs_dir>/<seq_name>/...` and stored as
`rep_canonical_pdb`.

When `passing_summary.csv` is empty (every prediction failed the
ra_eff gate), the **tier-none fallback** `_pick_rep_from_aggregated`
(`bin/cross_sequence_summary.py:354-443`) takes over.  It reads
`aggregated_results.csv` directly and picks the best composite for
cohort-plot visibility, and the row's `_picked_tier` is `"none"`.

**Bottom line:** the representative prediction for a (design,
sequence) cohort row can come from cold-start, steered, or reverted
predictions depending on which path each seed took.  The provenance
is preserved in the underlying columns (`steered_*`, `reverted_*`)
and the `rep_*` prefix marks summary-of-representative-seed.

### 4.13 Skip-steering (cold-start all-clean)

If all `num_seeds` cold-start predictions pass, the plan sets
`skip_steering=true` and the entire steering+reversion apparatus is
bypassed for that MPNN sequence.  Every seed's per-seed verdict
becomes `clean_steered`, `n_pass == n_seeds`, outcome
`no_reversion`, tier A.

This is the best computational outcome.  No off-interface mutations
exist to inflate any confidence metric, and Boltz chose the
designed interface as its primary (not fallback) prediction.

It is currently ranked equivalently to a steered tier A.  Visual
indicators (a star or snowflake on the cohort plot's tier cell)
were discussed but not added (Session 1 Q71, design_audit.md:L1751-1764).

### 4.14 Controls

Two controls run alongside every cohort:

- `control_scrambled` — the input PDB sequence with the design-region
  residues shuffled (preserves composition, destroys order).
- `control_polyA` — the input PDB sequence with the design region
  replaced by polyalanine.

Both are run through the *same* negsteer pipeline.  They should
*never* reach tier A or B.  If they do, that flags miscalibration of
thresholds.  `workflow.onComplete` in `main.nf:1133-1229` reads the
final cross-summary and warns on any `control_*` row that breaches
either guard.  The guards are
`rep_ipsae_min_median > params.controls_warning_ipsae_max` (default
0.5) and
`rep_ra_eff_vs_truth_median < params.controls_warning_ra_eff_min`
(default 5.0).

### 4.15 Important non-obvious subtleties

- **One seed = one Boltz prediction**, not a thread.  Seeds are RNG
  controls that produce independent predictions.  Stages do not
  carry "this seed's trajectory".  Each stage runs all `num_seeds`
  predictions afresh on its current sequence.
- **The reverted prediction is the final prediction.**  When
  reversion ran for a (design, seed) pair, the reverted PDB and its
  confidence metrics override the steered ones.  The `steered_*`
  columns still carry the steered observation for diagnostics.
- **Reversion always runs all `num_seeds` predictions** of a reverted
  sequence (CL-4).  Never partial.
- **AF3 is informational only** for tier assignment.  AF3
  `best_ra_eff` is reported in the cross-summary and may turn the
  AF3 cell red, but it does NOT change tier.  Tier is purely
  Boltz-driven.
- **Three-axis interpretation** of the cohort plot.  Tier (Boltz
  placement) is the row colour, the orthogonal filter result is the
  cell red/green, and Boltz confidence flags are a third dimension.
  They are intentionally independent.
- **MPNN sequence ID** is the join key everywhere downstream.  The
  cross-summary, the orthogonal cascade, the manifest — all join on
  `mpnn_sequence` or its hash.

---

## 5. Orthogonal validators

Orthogonal = different model/method from Boltz.  Goal: independent
evidence that a tier-A/B/C sequence is real, not a Boltz-specific
artifact.

### 5.1 Survivor manifest

`bin/extract_survivor_manifest.py` reads `cross_sequence_summary.csv`
and emits one row per cohort entry with paths to:
- canonical PDB (`rep_canonical_pdb`),
- ground-truth complex PDB,
- effector template CIF,
- the FASTA for the receptor.

`params.orthogonal_tier_filter` (default `'all'`) controls which
rows the manifest includes:
- `'all'`: every steered design including tier-none failed designs.
  Useful diagnostic, but more expensive (per-design AF3).
- `'abc'`: restrict to `cross_tier ∈ {A, B, C}`.  Lower GPU cost.

Control rows (`row_type != 'steered'`) are currently **not**
filtered by tier, and they always flow through.

### 5.2 AF3 no-MSA

`AF3_SETUP_DB` builds the AF3 database symlink farm once per run.
`AF3_NOMSA_ON_SURVIVORS` predicts each survivor with AlphaFold3
**without an MSA**.  `params.af3_nomsa_seeds` (default 3) controls
the number of seed × diffusion-sample pairs (15 predictions per
survivor).  `bin/parse_af3_output.py` post-processes one output into
a per-survivor CSV row.

AF3 aggregation in the cross-summary is by **best** (minimum)
`ra_eff` across the 15 predictions.  The column is labelled
`(best)` to make the inconsistency with the median-aggregated Boltz
columns visible.  This is intentional leniency because AF3 currently
fails on every production sequence.  Aggregation will be revisited
when the AF3 hit-rate is non-zero.

### 5.3 Biophysical metrics

`NEGSTEER_BIOPHYSICAL_METRICS` → `bin/run_biophysical_metrics.py`.
Per survivor on the **Boltz** representative PDB:
- BSA (FreeSASA),
- interface H-bonds (geometric N/O within 3.5 Å cross-chain).

Informational unless thresholded.

### 5.4 Rosetta InterfaceAnalyzer

`NEGSTEER_ROSETTA_METRICS` → `bin/run_rosetta_metrics.py` +
`bin/fastrelax_for_ia.xml`.  Per survivor:
1. FastRelax with 1 repeat, Cartesian, `ref2015_cart` (Bennett 2023
   calibration).
2. InterfaceAnalyzer: Sc (shape complementarity), `dG_separated`
   (≈ ΔΔG of binding in REU), `dSASA_int`, `delta_unsatHbonds`.

ΔΔG ≤ −30 REU floor is the calibrated threshold from Bennett et al.
2023 (`PipelineInternalThresholds.orthogonal_ddg_max`).

### 5.5 Interface metrics

`NEGSTEER_INTERFACE_METRICS` → `bin/compute_interface_metrics.py`
appends per-row metrics to the cross-summary:
- iRMSD,
- fnat,
- DockQ (combined CAPRI: >0.23 acceptable, >0.49 medium, >0.80
  high),
- ipSAE_15 (interface-restricted at 15 Å),
- intact_core (≤ `interface_intact_threshold` = 5 Å),
- weighted_jaccard (Gaussian-decayed interface overlap, `mu=4.0 Å`,
  `2σ²=2.25`).

These are diagnostic.  Only weighted_jaccard is currently considered
for promotion into the composite.

### 5.6 Filter cascade

`NEGSTEER_ORTHOGONAL_METRICS` → `bin/merge_orthogonal_metrics.py`
joins all three streams on `mpnn_sequence` and computes:

```
passes_orthogonal_filters = 1 iff:
    Sc  >= params.orthogonal_filter_sc_min  (default 0.55)
AND BSA >= params.orthogonal_filter_bsa_min (default 600.0 Å²)
AND ΔΔG <= orthogonal_ddg_max               (calibrated -30 REU)
```

**AF3 is NOT in this gate** (Session 4 Q97, commit `47cb9f2`).  AF3
disagreement is reported separately as a flag, not a gate.

`interface_plddt` was historically in this cascade and was
**removed** because it is a Boltz metric, not an independent
orthogonal one.  It still appears under "Boltz confidence" in the
cohort summary.

Production pass rate (cohort run notes12): ~17% of cohort rows
(146/845) pass the non-AF3 filters.  BSA alone passes 85% of rows,
so that threshold is a floor, not a discriminator.

---

## 6. The cohort table

`cross_sequence_summary.csv` — one row per MPNN sequence + one row
per control (`control_scrambled`, `control_polyA`).

Column groups (three-axis framework, Session 5 Q121):

1. **Boltz-2 structure (drives tier)** — `cross_composite_score`,
   `ra_eff_vs_truth_median`, `independent_receptor_rmsd_median`,
   `true_jaccard_median`, `n_mut`, `n_pass`, `n_seeds`,
   `cross_tier`.

2. **Boltz-2 confidence (informational)** — `complex_plddt_median`,
   `iptm_median`, `ipsae_min_15_median`, `pae_pass_frac_median`,
   `interface_plddt_median`.

3. **Orthogonal (red/green display, independent of tier)** —
   `af3_best_ra_eff_(best)`, `orthog_sc`, `orthog_bsa`,
   `orthog_ddg`, `weighted_jaccard_median`,
   `passes_orthogonal_filters`.

The remaining columns carry provenance.  They are `mpnn_sequence`,
`corrected_receptor` (full MPNN sequence), `designed_residues`
(design-region slice) and `native_residues` (parallel native slice
for comparison).  Then come `design_id`, `seq_id`,
`rep_canonical_pdb` (full path), `rep_seed_index`,
`rep_aggregated_verdict` / `outcome` and `row_type`
(`steered | control_scrambled | control_polyA`).

Composite score: `cross_composite_score = true_jaccard_median − 0.05
× ra_eff_vs_truth_median` (`_composite_score`,
`bin/cross_sequence_summary.py:245`).  Coefficients are arbitrary
placeholders.  The threshold audit will revisit them.

`passing_summary.csv` — the per-MPNN-sequence wet-lab triage list.
Filtered to `PASSING_VERDICTS = {"", "no_reversion", "pose_holds"}`.
The empty string and `no_reversion` are kept because the test
fixture's pre-multi-seed rows had blank verdicts.  This is
path-agnostic by design.  Tier-A rows from skip_steering are
included because they are `no_reversion`.

---

## 7. Reference data and conventions

### 7.1 Chain naming

- **Input PDB** chains can be anything.  `params.receptor_chain` /
  `params.effector_chain` tell the preprocessing where to look.
- **Post-RFDiffusion / through negsteer / orthogonal**: chain
  A = receptor, chain B = effector.  Hardcoded in
  `bin/rfdiffusion_filter.py:write_split_pdb` and in the
  `PRED_REC_CHAIN`/`PRED_EFF_CHAIN` constants.
- **HADDOCK internal**: relabels to A/B before docking regardless
  of input chain IDs.

### 7.2 Coordinate frames

The `PositionSet` typed wrapper (Phase 4 §2.3) carries an explicit
frame label:
- `native` — input PDB numbering, may have gaps,
- `designed` — RFDiffusion 1-based contiguous,
- `prediction` — Boltz/AF3 1-based per-chain.

Set operations require matched frame + matched chain.  Cross-frame
conversion via `ContigSpec.designed_to_native` /
`native_to_designed` / `prediction_to_designed`.  Drops unmappable
positions by default, while `strict=True` raises.

### 7.3 Index bases

Recurring footgun.  Two scripts emit text-file lists used as input to
the negsteer plan:

- `bin/derive_design_region.py` emits **1-based** indices
  (`_design_region.txt`).
- `bin/derive_true_interface.py` emits **0-based** indices
  (`_true_interface.txt`).

Consumers combining both must re-base.

### 7.4 Contig string format

Canonical: `A1-10/5-7/A15-20 B`.

- **Forward-slash separated within a chain** (no commas, ever).
- **Space-separated between chains.**
- **Fixed segment**: chain-prefixed range, e.g. `A1-10`.  Inclusive.
- **De novo segment**: a length (resolved) or a length range
  (constraint).  As of Phase 4 always carries `(min_len, max_len)`
  with `min_len == max_len` for the resolved case.  Parser accepts
  bare integers (`5` → `(5,5)`) or ranges (`5-7` → `(5,7)`).
- **`0`** = RFDiffusion chain-break marker.  It parses as
  `BreakSegment` (no length, position-math methods skip).
- **Bare chain letter inside a block** (e.g. `A1-10/B/A15-20`) =
  `PassthroughSegment(chain)` (no length at spec level, resolved
  later by callers with PDB access).

Each `DesignedBackbone` has exactly one `ContigSpec`.

Memory file `project_contig_string_format.md` flags this as a
recurring Claude-mistake source.  Never use commas.  Always
chain-prefix fixed segments, and use bare lengths for de novo
segments.

### 7.5 Internal threshold catalogue

`PipelineInternalThresholds` (`bin/pipeline_thresholds.py`,
Phase 4 §2.13) — single source of truth for every hard-coded
threshold.  Selected production defaults:

| Threshold | Default | Used by |
| --- | --- | --- |
| `steered_ra_eff` | 5.0 Å | per-seed structural gate |
| `reverted_ra_eff` | 5.0 Å | reversion structural gate |
| `cold_start_rmsd_threshold` | 6.0 Å | cold-start receptor RMSD |
| `intact_threshold` | 5.0 Å | independent receptor/effector RMSD |
| `contact_cutoff` | 4.5 Å | mutated-position contamination check |
| `weighted_jaccard_pair_cutoff` | 8.0 Å | interface overlap pair window |
| `weighted_jaccard_sigma` | 1.5 Å | Gaussian decay σ |
| `weighted_jaccard_mu` | 4.0 Å | Gaussian decay μ |
| `composite_ra_eff_weight` | 0.05 | composite formula |
| `iptm_min` | 0.30 | Boltz confidence flag |
| `complex_plddt_min` | 0.70 | Boltz confidence flag |
| `ipae_max` | 15.0 | Boltz confidence flag |
| `pae_pass_frac_min` | 0.10 | Boltz confidence flag |
| `orthogonal_sc_min` | 0.55 | orthogonal filter |
| `orthogonal_bsa_min` | 600.0 Å² | orthogonal filter |
| `orthogonal_ddg_max` | -30.0 REU | orthogonal filter (Bennett 2023) |
| `orthogonal_af3_ra_max` | 5.0 Å | AF3 informational flag |
| `interface_plddt_trim_threshold` | 50.0 | interface metric trimming |
| `controls_warning_ipsae_max` | 0.5 | `onComplete` control warning |
| `controls_warning_ra_eff_min` | 5.0 Å | `onComplete` control warning |

Immutable.  `with_overrides(**kwargs)` returns a new instance.  The
Phase 4 §2.13 threshold-audit work item plans to drive these from a
YAML file with CLI overrides for per-campaign tuning.

---

## 8. File:line map (where each behavior is implemented)

| Behavior | File | Lines |
| --- | --- | --- |
| Pipeline orchestration | `main.nf` | 351-1127 |
| Branch A entry | `main.nf` | 376-523 |
| Branch B entry | `main.nf` | 527-554 |
| Common-path join | `main.nf` | 556 |
| RFD chain hardcoding | `bin/rfdiffusion_filter.py:write_split_pdb` | — |
| Cold-start prediction loop | `bin/boltz2_negative_steering.py:cmd_plan` | 1818, 1936-2015 |
| Skip-steering decision | `bin/boltz2_negative_steering.py:cmd_plan` | 2017-2035 |
| Steered design pool generation | `bin/boltz2_negative_steering.py:cmd_plan` | 1818-2159 |
| Steered prediction (one design) | `bin/boltz2_negative_steering.py:cmd_predict_one` | 2788 |
| Steered ranking | `bin/boltz2_negative_steering.py:cmd_collect` | 2912 |
| CL-3 gate (type-level) | `bin/stage_result.py:triggers_next_stage` | 278-282 |
| CL-3 gate (live pipeline) | `bin/boltz2_iterate_steering.py:cmd_build_contaminated` | 1532-1568 |
| Reverted-sequence builder | `bin/reversion.py:build_reverted_sequence` | 147 |
| Reverted-plan writer (with dedup) | `bin/reversion.py:write_reversion_plan` | 304 |
| Reverted prediction (one) | `bin/boltz2_iterate_steering.py:cmd_predict_reversion_one` | 2649 |
| Reverted-result harvester | `bin/reversion.py:harvest_reversion_results` | 832 |
| Per-design reversion verdict | `bin/reversion.py:classify_reversion_verdict` | 1204-1304 |
| Per-stage per-seed verdict vocabulary | `bin/stage_result.py:per_seed_verdicts` | 117-186 |
| Cross-stage per-seed verdict breakdown | `bin/boltz2_iterate_steering.py:_per_seed_verdict_breakdown` | 3728-3866 |
| `clean_steered` heuristic (cold-start-clean / no-reversion-needed) | `bin/boltz2_iterate_steering.py:_per_seed_verdict_breakdown` | 3813-3835 |
| Aggregate outcome derivation | `bin/boltz2_iterate_steering.py:_classify_outcome` | 3555-3725 |
| `no_reversion` outcome assignment | `bin/boltz2_iterate_steering.py:cmd_aggregate_per_sequence` | 4736-4743 |
| `n_pass` formula | `bin/boltz2_iterate_steering.py:cmd_aggregate_per_sequence` | 4748-4749 |
| `singleton` outcome assignment | `bin/boltz2_iterate_steering.py:cmd_aggregate_per_sequence` | 4764 |
| Per-stage `n_pass` (type) | `bin/stage_result.py:n_pass` | 222-246 |
| Canonical seed picker | `bin/boltz2_iterate_steering.py:_canonical_seed` | 4654-4683 |
| Tier derivation | `bin/cross_sequence_summary.py:_tier_for_row` | 205-239 |
| Composite score | `bin/cross_sequence_summary.py:_composite_score` | 245 |
| Cross-rank assignment | `bin/cross_sequence_summary.py:_assign_cross_ranks` | 553-601 |
| Representative pick | `bin/cross_sequence_summary.py:_pick_representative` | 324-351 |
| Tier-none fallback | `bin/cross_sequence_summary.py:_pick_rep_from_aggregated` | 354-443 |
| Orthogonal filter cascade | `bin/merge_orthogonal_metrics.py` | — |
| `passes_orthogonal_filters` predicate | `bin/orthogonal_metrics.py:passes_orthogonal_filters` | — |
| Controls warning | `main.nf:workflow.onComplete` | 1133-1229 |

---

## 9. Known divergences (audit vs spec vs code)

The repo carries three layers of design intent that occasionally
disagree.  When they do, **code wins**.  This document follows code.

**Divergence 1 — Reversion gating rule.**

- `notes/design_audit.md` Session 1 (Q31, Q33, ~L789-895):
  describes the pre-Phase-4 rule — "all steered residues that
  contact effector are reverted; reversion runs whenever any seed
  has any contact."
- `notes/phase4_architecture_spec.md` §"CL-3" (L72-92): introduces
  the new majority gate, marks it "TO BE IMPLEMENTED."
- **Code (current, phase4-impl)**: CL-3 is **implemented** at
  `bin/stage_result.py:278-282` and `bin/boltz2_iterate_steering.py:1532-1568`.
- This document, §4.6: CL-3 is the canonical behavior.

**Divergence 2 — `clean_steered` semantics.**

- `notes/design_audit.md` Session 1 (Q41) uses "clean_steered"
  loosely to mean "no contamination regardless of placement."
- `notes/design_audit.md` Session 2 (Q60, L1542-1555) corrects this:
  a wrong-placement seed is `no_data`, not `clean_steered`.
- **Code (current)**: a seed is `clean_steered` iff intact AND
  ra_eff < cutoff AND no contamination on mutated positions
  (`bin/stage_result.py:163-180`).  Wrong-placement seeds are
  `pose_collapses` (in `StageResult.per_seed_verdicts`) or
  `no_data` (in `_per_seed_verdict_breakdown`).
- This document, §4.8: follows the code definition.

**Divergence 3 — Tier-A shortcut on `no_reversion`.**

- Historical code shortcut `outcome == "no_reversion" → tier A`,
  which would mis-label the wrong-placement-no-contamination case
  as tier A.
- **Code (current)**: tier is derived from `n_pass / n_seeds` only
  (`bin/cross_sequence_summary.py:205-239`).  The docstring at
  `:220-225` explicitly addresses this case.
- This document, §4.11: tier derivation is `n_pass`-only.

**Divergence 4 — HADDOCK auto-pick key.**

- `notes/design_audit.md` Session 7 (Q134, L3117-3126): initially
  recorded "largest BSA cluster representative."
- Session 7 (Q154, L3935-3949) amends: lexicographic
  `(-pair_contact_fraction, -bsa)`.
- The Q134 row at L3117-3126 was not updated to match the
  amendment.
- **Code**: see `bin/select_haddock_cluster.py`.
- This document, §3.2: follows the amendment.

**Divergence 5 — `n_cycles` and multi-cycle behavior.**

- `notes/design_audit.md` Session 1 (Q36): multi-cycle is "live
  roadmap," not dead code.
- **Code**: `cmd_iterate_*` and `cmd_kickoff*` (~600 LOC) exist but
  are dormant.  Nextflow only invokes the cycle-0 path
  (`kickoff-distances` → `kickoff-prefilter` → `kickoff-finalize
  --max-cycles 0`).
- Production runs `n_cycles=1`.
- This document, §4.3: treats cycle 0 as the only production cycle
  and notes the multi-cycle scaffolding.

---

## 10. Glossary cross-reference

The canonical glossary is
`notes/inventory/06_ubiquitous_language.md`.  This document uses
its terms consistently.  Quick pointers to that file for the most
recurring terms (line ranges below are in `06_ubiquitous_language.md`):

| Term | Definition source |
| --- | --- |
| Backbone vs Cα trace | §3, L191-201 |
| Motif vs scaffold (Baker-lab framing) | §3, L208-217 |
| Design region | §3, L219-224 |
| Anchor / hotspot / true interface / wrong interface / protected set | §3, L236-260 |
| Cold-start prediction | §3, L272-275 |
| Steered design / reverted design | §3, L277-285 |
| Steering mode | §3, L287-295 |
| Reliance flags | §3, L302-309 |
| Seed / sequence group / sequence | §4, L322-335 |
| Representative seed | §4, L337-340 |
| Per-seed reversion verdicts | §5, L373-381 |
| `n_pass` | §5, L402-403 |
| `PASSING_OUTCOMES` | §5, L405-406 |
| Tier | §5, L408-418 |
| Composite score | §5, L430-433 |
| `passes_orthogonal_filters` | §5, L435-444 |
| `ra_eff` | §6, L476 |
| Intact predicate | §6, L477-485 |
| iRMSD / fnat / DockQ | §6, L477-485 |
| Weighted Jaccard | §6, L489-498 |
| BSA / Sc / ΔΔG | §6, L495-498 |

§10 of `06_ubiquitous_language.md` lists 11 flagged ambiguities
(F1-F11).  These are work items for future rename / consolidation
passes, and they do not affect runtime behavior.

---

## Maintenance

When code changes any behavior described here, update this document
in the same commit.  At minimum keep §4 and §8 synchronized with the
negsteer code.  Bump the `last_verified` date in the frontmatter on
every substantive update so future Claude sessions can detect
staleness.

When this document and `notes/phase4_architecture_spec.md` agree,
they are both canonical.  When they disagree, prefer this document
for *operational* questions ("what does the code do today?") and the
architecture spec for *design-intent* questions ("what should the
code do?").
