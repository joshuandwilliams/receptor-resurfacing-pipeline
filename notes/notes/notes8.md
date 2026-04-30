# Notes 8 — Deduplication, multi-seeding, and metric reliability

## Context

Notes 7 identified two immediate next features: (1) multi-seed validation of
final candidates to measure confidence reproducibility, and (2) reversion
deduplication to avoid predicting the same sequence multiple times. This session
built both, discovered and fixed three reporting bugs in the aggregate code,
uncovered a fundamental finding about Boltz confidence metric reliability through
a 200-prediction variance test, and deployed a configurable `--num-seeds` parameter
that unifies the exploration-vs-confidence tradeoff into a single dial.

## What was built

### 1. Sequence deduplication in design generation

**Problem:** `make_steered_sequence` samples randomly from the candidate pool.
With small pools (e.g. 8 positions, 3 mutations, mild mode = 4 residue choices),
two designs can land on the same sequence by chance. Each duplicate wastes a GPU
prediction on an identical input.

**Solution:** The design loop in `cmd_plan` (cycle 0) and `cmd_iterate_plan`
(cycle 1+) now generates unique sequences with retry-on-collision. A set of seen
sequences is maintained; if a collision occurs, the RNG draws again. If the
combinatorial space is exhausted before reaching `--n-designs`, the code warns
and caps at the maximum number of unique sequences. `plan.json` records
`n_unique_sequences` and `total_predictions`.

### 2. Configurable `--num-seeds` parameter

**Problem:** A single Boltz prediction per sequence gives one sample from a
stochastic process. Confidence metrics (ipSAE, iPTM) and even structural
placement (ra_eff) can vary across seeds. Without multi-seed data, a single
lucky or unlucky draw determines whether a design passes or fails.

**Solution:** New `--num-seeds N` CLI parameter (default 1). Each unique sequence
gets N independent Boltz predictions with different seeds. The total SLURM array
size is `n_unique_sequences × num_seeds`.

- `--n-designs 100 --num-seeds 1` = maximum exploration (100 unique sequences,
  each predicted once)
- `--n-designs 20 --num-seeds 5` = same GPU budget, 20 sequences with 5-fold
  confidence intervals
- `--n-designs 5 --num-seeds 20` = deep validation of a tiny set

Each design entry in `plan.json` carries `sequence_group`, `seed_index`, and
`boltz_seed` fields. `steered_results.csv` gains `sequence_group` and
`seed_index` columns. When `num_seeds > 1`, a separate
`steered_results_aggregate.csv` is written with per-group mean/std/min/max
of structural metrics.

Design naming: `design_00` when `num_seeds=1`; `design_00_s0`, `design_00_s1`,
... when `num_seeds>1`. Backward compatible — `predict-one` falls back to the
old seed formula when `boltz_seed` is absent from plan.json.

### 3. Reversion deduplication

**Problem:** When `--max-mutations 3` and all 3 mutations contact the effector,
the reversion pass reverts them all back to the MPNN wild-type. In the v3 run
(notes 7), ~40 designs each triggered a separate Boltz prediction of the same
sequence — pure waste.

**Solution:** `write_reversion_plan` in `reversion.py` now groups contaminated
designs by their actual reverted sequence (hash-based). If 40 designs all revert
to the same sequence, only one set of predictions is staged (× `num_seeds`).
The `reversion_plan.json` manifest carries an `all_labels` field listing which
designs share each prediction. `harvest_reversion_results` copies results to all
shared labels after computing metrics for the canonical prediction.

**Verified on cluster:**
- v1 test (10 designs, num_seeds=1): 5 contaminated → 4 unique → 4 predictions
  (saved 1 GPU job). design_01 and design_02 shared results correctly.
- v2 test (10 designs, num_seeds=3): 15 contaminated → 5 unique → 15 predictions
  (5 × 3 seeds). Without dedup would have been 45 predictions — saved 30 GPU jobs.

### 4. Multi-seed variance test

**Standalone script** (`multiseed_variance_test.py` + `submit_multiseed_variance_test.sh`)
that predicts a single sequence N times with different Boltz seeds and reports
per-seed metrics plus aggregate statistics. Used to characterise the baseline
variance of the MPNN wild-type sequence before committing to a `num_seeds` strategy.

### 5. Final validation for `num_seeds=1` runs

**Standalone script** (`submit_final_validation.sh`) that reads `passing_summary.csv`,
extracts unique survivor sequences, and runs multi-seed validation on each. This is
the "option E" step: when `num_seeds=1` was used for broad exploration, the
survivors get deep validation post-hoc. Seeds use a 200000+ offset to avoid
collision with the main pipeline.

### 6. Sequence registry infrastructure

`sequence_registry.py` provides hash-based tracking of which sequences have been
predicted, with how many seeds, and where results live. Designed for cross-phase
dedup (e.g. not re-predicting the cold-start sequence during reversion), but not
yet wired into the SLURM chain. The within-phase dedup (steered dedup, reversion
dedup) handles the most impactful cases.

## Bugs discovered and fixed

### Bug 1: `_populate_reverted_mutations` — wrong lookup key for cycle-0 rows

**What:** The aggregate code used `r["pathway"]` as the lookup key for
contaminated.json and reversion_results.json. For cycle-0 rows, `pathway` is
the constant `"cycle_0"`, not a per-design identifier. The per-design identity
is in `r["design"]` (e.g. `"design_04"`).

**Impact:** Three downstream failures:
- `steered_mutated_contact_positions` was blank for all cycle-0 reversion rows
  (notes 7 §3 bug — now explained and fixed)
- `positions_to_revert` lookup returned empty set → `reverted_total_mutations`
  was always 0 (accidentally correct for fully-reverted designs, wrong for
  partially-reverted ones)
- `read_cumulative_mutations("design_04")` crashed on `parse_pathway_label`
  (expects `cNdMM` format), silently caught → empty mutations list

**The critical consequence:** `passing_summary.csv` reported `total_mutations=0`
for pose_holds designs that actually carried 1 remaining mutation. This led to
the false belief that the MPNN wild-type sequence could find the correct interface
(4 out of 45 predictions), when in fact those predictions were of single-point
mutants (e.g. D32K).

**Fix:** For cycle-0, derive `leaf_label` from `r["design"]` instead of
`r["pathway"]`. Read `positions_to_revert` from `contaminated.json` (the
authoritative source) instead of `reversion_results.json` (which doesn't carry
it). Read mutations.tsv directly from the design's PDB directory for cycle-0
instead of calling `read_cumulative_mutations`.

### Bug 2: `positions_to_revert` read from wrong file

**What:** The aggregate code at line 2837 read `positions_to_revert` from
`reversion_results.json`, but `harvest_reversion_results` never writes that key
into the results dict. The authoritative source is `contaminated.json`.

**Fix:** Look up from `contaminated.json` first, fall back to
`reversion_plan.json` entries.

### Bug 3: reversion.py intact cutoff still 3.0 Å

**What:** When the intact cutoff was relaxed from 3.0 to 5.0 Å across the
codebase, `harvest_reversion_results` in `reversion.py` line 593 was missed.

**Fix:** Changed to 5.0 Å.

## Key experimental findings

### 200-prediction variance test of the MPNN wild-type

Ran the MPNN receptor sequence (design_3) 100 times with `diffusion_samples=1`
and 100 times with `diffusion_samples=5`, each with different seeds.

| | ds=1 | ds=5 |
|---|---|---|
| Intact (rec+eff RMSD ≤ 3 Å) | 63/100 | 97/100 |
| Correct interface (ra_eff < 5 Å) | 0/100 | 0/100 |
| ra_eff mean ± std | 30.0 ± 0.8 Å | 29.8 ± 0.7 Å |
| ipSAE mean ± std | 0.689 ± 0.099 | 0.787 ± 0.055 |
| iPTM mean ± std | 0.845 ± 0.050 | 0.893 ± 0.020 |

**Key findings:**

1. **The MPNN wild-type NEVER finds the correct interface** across 200
   independent predictions. The wrong binding site is the overwhelming Boltz
   solution for this sequence. ra_eff is rock-solid at ~30 Å with <1 Å std.

2. **`diffusion_samples=5` dramatically stabilises predictions** — effector
   deformation drops from 37% to 3%, and confidence metric variance roughly
   halves. Boltz picks better models when it has 5 internal candidates. This
   does NOT help with interface placement (still 0/100), only with fold quality.

3. **This resolves the notes 7 mystery.** Notes 7 observed 4/45 reversion
   predictions of "the MPNN sequence" finding the correct interface. We showed
   those predictions were NOT the MPNN wild-type — they carried a remaining
   D32K mutation (the `reverted_total_mutations=0` was a reporting bug). The
   single point mutation D32K is sufficient to shift the effector toward the
   correct interface in ~10% of predictions.

### ipSAE does not distinguish correct from incorrect interface

| Sequence | ra_eff | ipSAE | iPTM | Interface |
|---|---|---|---|---|
| MPNN wild-type (100 seeds, ds=5) | 29.8 | 0.787 | 0.893 | WRONG |
| design_24 (10 seeds, validated) | 2.5 | 0.755 | 0.862 | CORRECT |
| design_12_s0 (reverted, D32K) | 2.6 | 0.522 | 0.754 | CORRECT |
| design_10_s2 (reverted, D32K) | 3.6 | 0.245 | 0.649 | CORRECT |

Boltz is MORE confident about the wrong interface (ipSAE 0.787) than the correct
one (ipSAE 0.755). The pose_holds survivors with D32K have even lower confidence
(0.25–0.52) despite being structurally correct.

ipSAE and iPTM measure Boltz's self-consistency (how well-packed it thinks the
interface is), NOT whether the prediction is biologically correct. The wrong
interface is a perfectly reasonable protein surface that Boltz has seen in
training data; the correct interface (designed de novo) is novel, so Boltz is
less confident even when it finds it.

**This finding needs deeper investigation.** Notes 5 reported design_3 as the
first case where ipSAE showed meaningful signal (0.86), and concluded the
notes-3 decision to keep ipSAE as the production ranker was correct. The current
session's data complicates that conclusion — the 0.86 was a single seed that
happened to score high, the 10-seed mean is 0.755, and the wild-type wrong-site
prediction scores comparably. A thorough literature review and analysis of the
available data is needed to determine whether ipSAE, iPTM, or another metric
is appropriate for ranking in this pipeline. See "What to do next" §1.

### design_24 validated: 10/10 seeds find correct interface

Multi-seed validation of design_24 (K5E M22E L36D — the only clean steered
survivor from v3):

| Metric | Mean | Std | Min | Max |
|---|---|---|---|---|
| ra_eff | 2.51 | 0.53 | 1.88 | 3.34 |
| ipSAE_min | 0.755 | 0.055 | 0.679 | 0.846 |
| iPTM | 0.862 | 0.019 | 0.846 | 0.898 |
| rec RMSD | 1.89 | 0.60 | 1.03 | 3.02 |
| eff RMSD | 0.51 | 0.07 | 0.41 | 0.61 |

10/10 seeds find the correct interface (all ra_eff < 5 Å). The three mild
mutations (K→E, M→E, L→D — all charge-introducing) on the wrong-interface
surface reliably redirect Boltz to the correct binding site without any of
the mutations directly contacting the effector. This is a genuine,
reproducible steering hit.

### v3 vs v4 comparison: exploration vs confidence

| | v3 (100×1) | v4 (30×3) |
|---|---|---|
| Total predictions | 100 | 90 |
| Unique sequences | 100 | 30 |
| ra_eff < 5 Å | 57 | 48 |
| Clean steered | 1 (design_24) | 0 |
| pose_holds | 0 | 2 (D32K designs) |

v3's broad exploration found the only clean steered design (design_24) — a hit
v4 could have missed with only 30 sequences. v4's multi-seeding caught D32K
pose_holds that v3 missed because the single reversion seed collapsed. Both
approaches found valid hits, but different ones.

For production, `num_seeds=1` with post-hoc final validation (option E) appears
to be the better default — it finds the cleanest hits and validates them
cheaply afterward.

## Intact cutoff change: 3.0 → 5.0 Å

The `RECEPTOR_INTACT_CUTOFF` and `EFFECTOR_INTACT_CUTOFF` were hardcoded at
3.0 Å throughout the pipeline. This was too conservative for RFDiffusion
designs with floppy termini — design_16 in v4 had ra_eff 1.0–1.6 Å (best
ever) but was rejected because receptor RMSD was 3.2–3.4 Å.

Changed to 5.0 Å in:
- `boltz2_negative_steering.py` (module-level constants + cmd_collect)
- `boltz2_iterate_steering.py` (12+ locations across all filter stages)
- `reversion.py` (harvest_reversion_results)
- `compute-final-metrics` default `--rmsd-threshold`

Note: this threshold is intended as a temporary measure. The proper fix
is to replace whole-chain RMSD with a metric that ignores disordered
regions (core-only RMSD, iRMSD, or DockQ). See "What to do next" §3.

## Files changed

| File | Changes |
|---|---|
| `boltz2_negative_steering.py` | `--num-seeds` CLI; dedup design loop; `sequence_group`/`seed_index`/`boltz_seed` in plan.json; aggregate CSV; intact cutoff 3→5 |
| `boltz2_iterate_steering.py` | Cycle-0 aggregate bug fixes (leaf_label, positions_to_revert source, read_cumulative_mutations); cycle-1+ dedup design loop + num_seeds; num_seeds forwarding in all 4 resubmission paths; intact cutoff 3→5 across 12+ locations |
| `reversion.py` | Reversion dedup (group by sequence hash, predict once, share results); num_seeds for reversions; intact cutoff 3→5 in harvest |
| `sequence_registry.py` | New file — cross-phase sequence tracking (not yet wired into SLURM chain) |
| `submit_boltz2_negative_steering.sh` | `--num-seeds` arg; array size = n_designs × num_seeds |
| `submit_boltz2_iterate_steering.sh` | `--num-seeds` arg; array sizing; forwarding through collect/finalize |
| `multiseed_variance_test.py` | New file — plan/predict/harvest for variance testing |
| `submit_multiseed_variance_test.sh` | New file — SLURM submission for variance test |
| `submit_final_validation.sh` | New file — post-pipeline multi-seed validation for num_seeds=1 |
| `extract_passing.py` | No changes this session (ranking question deferred) |
| `compute_metrics.py` | No changes this session |

## Run commands (current recipe)

```bash
# Submit (num_seeds=1 for broad exploration)
./submit_boltz2_negative_steering.sh \
    --workdir         runs/design_N_vX \
    --ground-truth    resurface_pipeline_test/design_N.pdb \
    --receptor        A \
    --effector        B \
    --receptor-fasta  resurface_pipeline_test/design_N_seq_M_receptor.fasta \
    --mode            mild \
    --max-mutations   3 \
    --candidate-pool-size 12 \
    --n-designs       100 \
    --num-seeds       1 \
    --n-cycles        1 \
    --diffusion-samples 5 \
    --true-interface-indices-file resurface_pipeline_test/design_N_true_interface.txt \
    --design-region-indices-file  resurface_pipeline_test/design_N_design_region.txt

# After queue drains:
python3 boltz2_iterate_steering.py aggregate \
    --experiment-root runs/design_N_vX

singularity exec $CONTAINER python3 boltz2_iterate_steering.py compute-final-metrics \
    --experiment-root runs/design_N_vX \
    --rmsd-threshold 5.0 \
    --metric-column steered_ra_eff_vs_truth \
    --contact-cutoff 5.0

python3 extract_passing.py \
    --input runs/design_N_vX/all_results_multicycle_with_metrics.csv

# Final validation of survivors (option E):
./submit_final_validation.sh \
    --experiment-root runs/design_N_vX \
    --n-seeds 10 --diffusion-samples 5
```

## Known limitation: reversion dedup audit trail

When multiple contaminated designs revert to the same sequence, the harvest
computes metrics once using the canonical design's `reversion_metadata.json`.
The `reverted_mutated_contact_positions` column in shared results reflects the
canonical design's remaining mutations, not the non-canonical label's. This
doesn't affect verdicts or structural metrics (those depend only on the
sequence, which is shared), but the audit column may be misleading for
non-canonical labels. Low priority — the verdict is authoritative.

## What to do next

### 1. Literature review and metric analysis (HIGH PRIORITY)

Conduct a thorough analysis of which metrics are appropriate for ranking
negative-steering designs. The current session showed ipSAE doesn't
distinguish correct from incorrect interface placement on this system, but
notes 5 showed it had signal on an earlier run. Questions to resolve:

- Is ipSAE reliable for ranking designs that ALL find the correct interface
  (i.e. within the correct-interface regime, does higher ipSAE correlate
  with better binding)?
- What does the Boltz/AF2 literature say about ipSAE / iPTM reliability
  on de novo designed interfaces vs native interfaces?
- Should `extract_passing.py` rank by ra_eff instead of ipSAE, or by a
  composite metric?
- Are there alternative metrics (DockQ, lDDT-PPI, fnat, interface-only
  RMSD) that would be more appropriate?

### 2. Reduce contamination rate (HIGH PRIORITY)

In v3, 49/57 designs that found the correct interface had contaminating
mutations (mutated residues contacting the effector). This means the
steering mutations are doing the work of redirecting the effector but then
landing in the new interface, invalidating the prediction. Possible
approaches:

- **Expand the protected set.** Currently, only the true interface is
  protected from mutation. If we also protected residues within a
  second-shell radius of the true interface, fewer mutations would land
  near the correct binding site.
- **Better candidate pool selection.** Currently, the candidate pool is
  the top-N closest wrong-interface contacts that are surface-exposed and
  not in the true interface. Selecting positions that are maximally
  DISTANT from the true interface (not just absent from it) would reduce
  contamination.
- **Post-hoc analysis of which positions contaminate.** Position 32 (D→K/R)
  appeared consistently as a contaminating contact. Understanding which
  candidate-pool positions tend to contaminate would allow pre-filtering.

### 3. Replace ra_eff with a disordered-region-aware metric (MEDIUM PRIORITY)

Carried from notes 3, 4, 6. `ra_eff` is a whole-effector Cα RMSD that is
inflated by floppy termini and disordered loops. Options:

- **iRMSD** (interface-only RMSD): restrict to effector residues within
  a cutoff of the receptor in ground truth. Easiest drop-in.
- **DockQ**: the community standard for docking quality assessment.
  Combines iRMSD, LRMSD, and fnat into a single score. Biggest refactor
  but most standard and publishable.
- **Core-only RMSD**: trim low-confidence regions (low pLDDT or
  high PAE) before computing RMSD.
- **lDDT-PPI**: local distance difference test restricted to interface
  residues.

### 4. Threshold audit (MEDIUM PRIORITY)

The pipeline has accumulated multiple thresholds across different files and
stages, some inherited from early development, some changed during this
session. A thorough review should:

- Document every threshold, where it's used, and what its current value is
- Check consistency (e.g. the intact cutoff was 3.0 in some places and 5.0
  in others until this session)
- Determine correct values for production use
- Consider making key thresholds CLI-configurable where they aren't already

Key thresholds to audit:
- `RECEPTOR_INTACT_CUTOFF` / `EFFECTOR_INTACT_CUTOFF` (now 5.0 Å)
- `--rmsd-threshold` for compute-final-metrics (now 5.0 Å)
- `--contact-cutoff` for interface detection (4.5 Å heavy-atom)
- `--novelty-cutoff` for cycle-to-cycle novelty filter (10.0 Å)
- `build-contaminated` ra_eff threshold for sending to reversion (5.0 Å)
- `classify_reversion_verdict` structural cutoff (5.0 Å)
- `pae_cutoff` for ipSAE computation (10.0 default)
- Confidence flag thresholds in extract_passing.py

### 5. Wire sequence registry into SLURM chain (LOW PRIORITY)

The `sequence_registry.py` module is built and tested but not connected
to the pipeline. Integration points:

- `cmd_plan` should register the cold-start prediction and all steered
  sequences after plan.json is written
- `plan-reversions` should check the registry before staging reversion
  predictions, skipping sequences that were already predicted with
  sufficient seeds
- Each cycle's plan should check the registry for sequences predicted
  in earlier cycles

Currently, within-phase dedup (steered collisions, reversion grouping)
handles the most impactful cases. Cross-phase dedup would save additional
GPU time on the cold-start-as-reversion case.

### 6. Remaining carried-forward items

From notes 6:
- **Production tier update.** The wet-lab triage criterion needs
  formalising with the reversion pass: `reversion_verdict == 'pose_holds'
  OR (clean_steered AND intact AND ra_eff < threshold)`.
- **Nextflow migration.** `boltz2_iterate_steering.py` is now ~4500 lines.
  Natural split points exist. Deferred until Nextflow DSL2 migration.
- **Silent pass-through of unvalidated contaminated designs.** If
  `plan-reversions` crashes, designs should be dropped, not silently kept.
  Guard partially implemented (Bug 1 guard in finalize), but the general
  case (n_contaminated vs n_verdicts mismatch) isn't checked.

From notes 3/4:
- **RFDiffusion production validator wrapper** — iterate over a directory
  of designs, run the full pipeline per design, produce one summary CSV.

## Things to remember when working on this code

All items from notes 6/7 still apply, plus:

- **`--num-seeds` is inherited from cycle-0 plan.json.** Cycle-1+ plans
  read `num_seeds` from `cycle0_plan.get("num_seeds", 1)` and use the
  same value throughout. You cannot change num_seeds mid-experiment.
- **Design naming changes with num_seeds.** `design_00` when num_seeds=1;
  `design_00_s0`/`design_00_s1`/... when num_seeds>1. Code that parses
  design names must handle both formats.
- **`steered_results.csv` now has `sequence_group` and `seed_index` columns.**
  Downstream consumers (aggregate, extract_passing) should be aware of
  these. The aggregate code currently passes them through without
  multi-seed-aware grouping — this is acceptable because each seed is
  treated as an independent design for the reversion pipeline.
- **The intact cutoff is now 5.0 Å everywhere.** If you re-run
  aggregate + compute-final-metrics on an OLD experiment, you must also
  re-run `collect` (which writes `steered_results.csv` with the intact
  flag) to pick up the new threshold. Alternatively, delete the stale
  `all_results_multicycle_with_metrics.csv` to prevent resume from
  preserving old intact values.
- **`--n-cycles 1` works correctly for reversion.** Fixed in notes 7.
  No need for `--n-cycles 2` workaround.
