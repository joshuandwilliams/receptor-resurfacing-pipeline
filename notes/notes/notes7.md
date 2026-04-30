# Notes 7 — Contamination fix, schema cleanup, and validation runs

## Context

Notes 6 built the full reversion validation pass: for any steered design whose
rescue depends on mutations at "protected-set positions" (design region ∪ true
interface), revert those mutations and re-run Boltz to check if the pose holds.
This session discovered that the protected-set concept itself was fundamentally
wrong, fixed it, ran three validation experiments (v1 smoke, v2 full, v3 with
diffusion samples), and identified the next feature needed: multi-seed validation
of final candidates.

## The critical bug: protected-set contamination filter

### What the code was doing

The contamination check gated on `mutated_positions ∩ contact_residues ∩
protected_set`. A mutation at a position *outside* the protected set (design
region ∪ true interface) that was contacting the effector in the steered prediction
was considered acceptable — the code treated these as "proposed point-mutants the
construct can carry."

### Why that was wrong

The wet-lab construct is synthesised from the original ProteinMPNN sequence. It
does not carry *any* steering mutations — not at protected positions, not at
unprotected positions, not anywhere. A steering mutation at position 9 (native
scaffold, outside the design region, outside the true interface) that contacts the
effector is exactly as invalidating as a mutation at position 44 (inside the design
region): neither residue will exist in the expressed protein. The binding metrics
measured in the presence of that contact are describing a structure that will never
exist.

### Evidence from the v1 run data

The original 50-design run (pre-fix) had 33 designs marked "clean" that all had
mutations at positions 5, 9, 32, 36, 44, 46, 62, or 69 directly contacting the
effector. Only 2 designs were sent to reversion (those with mutations at position
44, which happened to be inside the design region). The other 31 were passed through
with invalid metrics. `n_mutated_contact_on_protected_positions = 0` for all 33
because their contaminating mutations were at positions outside the protected set.

### The fix

```
positions_to_revert = mutated_positions ∩ contact_residues
```

No protected-set filter. Any mutated residue contacting the effector is
contamination. Applied consistently in three places:

1. `cmd_build_contaminated` — flags designs for reversion
2. `harvest_reversion_results` — computes `reverted_mutated_contact_positions`
3. `classify_reversion_verdict` — fires `new_contamination` when remaining
   mutations contact the effector in the reverted prediction

### What was removed

The entire protected-set machinery:

- `classify_construct_reliance` function (compute_metrics.py)
- `construct_reliance_flag`, `n_contact_on_protected_positions`,
  `n_mutated_contact_on_protected_positions`, `mutated_contact_protected_positions`
  columns from CSV schema and all emitters
- `--protected-positions` and `--protected-positions-file` CLI flags from
  compute_metrics.py
- `--protected-positions-file` CLI flag from `cmd_harvest_reversions`
- `--design-region-positions-file` CLI flag from `cmd_compute_final_metrics`
- ~100-line protected-set resolution block in `cmd_compute_final_metrics`
- ~70-line protected-set resolution block in `cmd_harvest_reversions`
- `protected_positions` parameter from `harvest_reversion_results` signature
- `protected_positions` field from `contaminated.json` output
- `reverted_construct_reliance_flag` from all finalize functions and aggregate

### What was added

- `steered_mutated_contact_positions` — comma-separated list of positions where a
  steering mutation contacts the effector (the contamination receipt; the exact
  positions that get reverted)
- `reverted_mutated_contact_positions` — same for the reverted prediction (should
  be empty for `pose_holds` by definition; populated for `new_contamination`)
- `classify_mutation_reliance` now returns `(count, list)` — the list is the
  authoritative contamination identifier
- `steered_mutation_reliance_flag` was also removed as redundant —
  `n_contacts_on_mutated_positions > 0` is the contamination signal

## Additional bug: build-contaminated filtering on intact only, not ra_eff

### The bug

`build-contaminated` checked all `intact_candidates` from `prefilter.json`
regardless of their `steered_ra_eff_vs_truth`. A design with ra_eff = 28 Å but
intact = 1 (both chains fold fine, just docked completely wrong) would be sent
through the full reversion pipeline. If the reversion prediction happened to
find the correct interface by Boltz stochasticity, it would be marked
`pose_holds` — but this is not validating a steered rescue, it's just a lucky
re-roll of the dice on the same (or similar) sequence.

The v2 run had one such case: cycle 0 design_31 with `steered_ra_eff = 28.354`
ended up as `pose_holds` with `reverted_ra_eff = 3.599`.

### The fix

`build-contaminated` now filters candidates on `ra_eff < rmsd_threshold`
(default 5.0 Å) in addition to intact. Only designs whose steered prediction
actually found the correct interface are checked for contamination and sent to
reversion.

## Additional bug: n-cycles 1 skipped reversion entirely

### The bug

The reversion chain (stages 4a/b/c in `submit_boltz2_negative_steering.sh`) was
gated on `[[ "$N_CYCLES" -gt 1 ]]`. With `--n-cycles 1`, no reversion ever ran.
The `--allow-no-reversion` guard added in this session (before the root cause was
understood) forced users to use `--n-cycles 2`, but that was a workaround, not a
fix.

Reversion is not a "next cycle" feature — it validates cycle 0's own results.
It should always run.

### The fix

- Reversion chain always submits, regardless of `--n-cycles`
- `kickoff-finalize` receives `--max-cycles $((N_CYCLES - 1))` — when
  `--n-cycles 1`, max-cycles = 0, so no further cycles are spawned but
  reversion still runs
- Workdir always uses nested layout (`<workdir>/cycle_0/`) so reversion paths
  are consistent
- `--allow-no-reversion` flag removed entirely

## Additional bug: junk rows getting ranked

### The bug

`_row_is_clean_steered` only checked `steered_receptor_intact == 1` and
`reversion_verdict == ""`. Designs with intact chains but terrible ra_eff
(25+ Å) were being ranked by `rank_by_ra_eff`, then `extract_passing.py`
included any row with either `rank_by_ra_eff` or `rank_by_ipsae_min`.
This produced 13 junk rows in `passing_summary.csv` with ra_eff values of
5–34 Å and no confidence metrics.

### The fix

- `_row_is_clean_steered` now requires `steered_ipsae_min` to be populated
  (i.e. `compute-final-metrics` actually ran on the row). Designs that failed
  the ra_eff threshold won't have ipsae populated, so they won't be ranked.
- `extract_passing.py` now only includes rows with `rank_by_ipsae_min`
  (not `rank_by_ra_eff` alone)

## CSV schema changes

### Removed columns

| Column | Reason |
|---|---|
| `steered_construct_reliance_flag` | Protected-set concept removed |
| `steered_n_contact_on_protected_positions` | Protected-set concept removed |
| `steered_n_mutated_contact_on_protected_positions` | Protected-set concept removed |
| `steered_mutated_contact_protected_positions` | Protected-set concept removed |
| `reverted_construct_reliance_flag` | Protected-set concept removed |
| `steered_mutation_reliance_flag` | Redundant with `n_contacts_on_mutated_positions > 0` |

### Added columns

| Column | Meaning |
|---|---|
| `steered_mutated_contact_positions` | Comma-separated positions where a mutation contacts the effector (triggers reversion) |
| `reverted_mutated_contact_positions` | Same for reverted prediction (empty for pose_holds, populated for new_contamination) |
| `reverted_n_contact_residues` | Contact count on the reverted prediction |
| `reverted_contact_residues` | Contact residue list on the reverted prediction |
| `reverted_contact_cutoff_used` | Contact cutoff used (same as steered; from plan.json) |

### New file: `extract_passing.py`

Standalone script that reads `all_results_multicycle_with_metrics.csv` and produces
`passing_summary.csv` — a clean, flat table of only the ranked survivors with
unified columns (no `steered_`/`reverted_` prefix split). The `reverted` column
(0/1) indicates the source. For `pose_holds` rows, metrics come from the reverted
prediction; for clean steered rows, from the steered prediction. Computes
`confidence_flag` on the fly when the steered flag is blank (e.g. rows where the
steered ra_eff was too high for compute-final-metrics but the reverted ra_eff
passed).

Usage:
```bash
python3 extract_passing.py \
    --input  <experiment>/all_results_multicycle_with_metrics.csv \
    --output <experiment>/passing_summary.csv  # default
```

## Validation runs

### v1 smoke test (4 designs, 2 mutations, 2 cycles)

All 3 intact designs had mutation at position 44 contacting the effector →
all 3 sent to reversion → all 3 `pose_collapses` (reverted ra_eff 27–29 Å) →
0 survivors → cycle 1 never launched. Pipeline correctly killed all candidates.

### v2 full run (50 designs, 6 mutations, 2 cycles)

| Metric | Old run (broken) | v2 (fixed) |
|---|---|---|
| Sent to reversion | 2 | 67 |
| pose_holds | 0 | 13 |
| pose_collapses | 2 | 48 |
| new_contamination | 0 | 6 |
| Clean steered | 33 (most falsely clean) | 1 |
| Ranked survivors | 33 (invalid metrics) | 14 (genuine) |

Top survivors: reverted ipsae_min 0.17–0.33, ra_eff 1.5–2.4 Å. Moderate
but credible rescue signals.

### v3 run (100 designs, 3 mutations, 1 cycle, 5 diffusion samples)

Tested the `--n-cycles 1` fix (reversion now runs), wider candidate pool
(--candidate-pool-size 12), and `--diffusion-samples 5`.

| | Count |
|---|---|
| Total designs | 100 |
| Sent to reversion | 48 |
| pose_holds | 4 |
| pose_collapses | 40 |
| new_contamination | 4 |
| Clean steered | 1 |
| Ranked survivors | 5 |

Results:
- **Rank 1 (design_24):** Clean steered (no contamination). K5E M22E L36D.
  ipsae = 0.895, iptm = 0.91, pae_pass_frac = 0.97, ipae = 3.4 Å. Strongest
  rescue signal seen across all runs. These 3 mutations don't contact the
  effector — the interface is held entirely by native MPNN residues.
- **Ranks 2–4 (designs 73, 64, 04):** pose_holds with 0 remaining mutations.
  ipsae 0.38–0.49, iptm 0.77–0.79, pae_pass_frac > 0.8. Good confidence,
  but the scored sequence is identical to the cold start / MPNN sequence.
- **Rank 5 (design_90):** pose_holds, 0 mutations, ipsae = 0.0. No confidence.

### Key observation: redundant MPNN-sequence predictions

Many designs in v3 had all 3 mutations contacting the effector → all 3 reverted →
reverted sequence = original MPNN sequence. Combined with the cold start, the
MPNN sequence was predicted ~45 times (each with 5 diffusion samples = ~225 Boltz
trajectories). Of those, ~4 found the correct interface and ~41 didn't. This
gives a ~10% hit rate for the MPNN sequence finding the true interface —
it's accessible but not the dominant Boltz solution.

This redundancy is wasteful and should be deduplicated in future runs.

### Key observation: diffusion-samples vs seeds

`--diffusion-samples 5` generates 5 structure samples within a single Boltz
prediction; Boltz picks its top-ranked model (model_0) and the pipeline uses
that. This is NOT the same as running 5 independent predictions with different
seeds — it's Boltz's internal diversity within one run. The user expected
multiple independent cold starts; this was clarified.

## What to build next

### 1. Multi-seed validation of final candidates

**Goal:** For the survivors in `passing_summary.csv`, re-predict each unique
sequence multiple times with different seeds to measure confidence reproducibility.

**Why:** A single Boltz prediction (even with diffusion-samples > 1) gives one
confidence measurement. Two designs with ipsae 0.9 and 0.4 might both reproduce
at those levels across 10 seeds, or the 0.9 might scatter between 0.1–0.9 and
the 0.4 might be stable at 0.4. Without multi-seed data, we can't distinguish a
reliable prediction from a lucky trajectory.

**Design:**
- Input: `passing_summary.csv` (or a list of sequences to validate)
- Deduplicate sequences — e.g. all 0-mutation pose_holds rows are the same
  MPNN sequence; only predict it once
- Run each unique sequence N times (e.g. 10) with different seeds
- Collect ipSAE, iPTM, ipAE, ra_eff from each run
- Report: mean, std, and distribution of confidence metrics per sequence
- This is cheap: ~2 unique sequences × 10 seeds × 5 diffusion samples
  = 100 Boltz trajectories (vs 500+ in the main pipeline)

**Implementation options:**
- Standalone script that takes sequences + ground truth → submits Boltz array →
  harvests metrics → writes a summary CSV with per-seed rows and aggregate stats
- Alternatively, integrate into `extract_passing.py` as a `--validate` mode
  that auto-submits the validation jobs

### 2. Reversion deduplication

**Goal:** Avoid predicting the same sequence multiple times during the reversion
pass.

**Why:** When `max-mutations = 3` and all 3 mutations contact the effector (common
on design_3), every such design gets fully reverted to the MPNN sequence. 40+
designs each triggering a fresh Boltz prediction of the same sequence is pure waste.

**Design:**
- Before submitting reversion predictions, hash each reverted sequence
- If two designs produce the same reverted sequence, predict it once and share
  the result
- The reversion metadata (which positions were reverted, which design it came
  from) stays per-design; only the Boltz prediction is deduplicated
- Care needed around the `predict-reversion-one` array job structure: currently
  each array task predicts one design's reversion. With dedup, the array size
  shrinks to the number of unique sequences, and a mapping file tells
  `harvest-reversions` which designs share which prediction.

### 3. `steered_mutated_contact_positions` surfacing bug

The column is blank for all reversion rows in the aggregate CSV. The
`_populate_reverted_mutations` function reads `positions_to_revert` from
`reversion_results.json` but `steered_mutated_contact_positions` should come from
`contaminated.json`. The code to surface it was added but appears not to be
working — needs debugging. Low priority since `reversion_verdict` and
`reverted_total_mutations` carry the essential information, but the column should
be populated for audit completeness.

## Files changed

| File | Changes |
|---|---|
| `boltz2_iterate_steering.py` | Contamination fix (build-contaminated, finalize functions, harvest-reversions, aggregate fieldnames, METRIC_KEYS_RAW); ra_eff filter in build-contaminated; always-nested workdir layout; _row_is_clean_steered ipsae check; reverted contact columns wired through; all protected-set code removed |
| `reversion.py` | Contamination fix (classify_reversion_verdict uses reverted_mutated_contact_positions); harvest emits reverted_mutated_contact_positions; protected_positions parameter removed; module docstring updated |
| `compute_metrics.py` | classify_mutation_reliance returns (count, list); classify_construct_reliance removed; --protected-positions CLI removed; protected_positions parameter removed from parse_boltz2; CSV_FIELDS updated |
| `submit_boltz2_negative_steering.sh` | Reversion always runs (not gated on n-cycles); --allow-no-reversion removed; always-nested workdir layout |
| `extract_passing.py` | New file; produces passing_summary.csv |
| `columns.md` | Full rewrite for corrected contamination rule |

## Run commands (current recipe)

```bash
# Submit
./submit_boltz2_negative_steering.sh \
    --workdir         resurface_pipeline_test/runs/design_N_vX \
    --ground-truth    resurface_pipeline_test/design_N.pdb \
    --receptor        A \
    --effector        B \
    --receptor-fasta  resurface_pipeline_test/design_N_seq_M_receptor.fasta \
    --mode            mild \
    --max-mutations   3 \
    --candidate-pool-size 12 \
    --n-designs       100 \
    --n-cycles        1 \
    --diffusion-samples 5 \
    --true-interface-indices-file resurface_pipeline_test/design_N_true_interface.txt \
    --design-region-indices-file  resurface_pipeline_test/design_N_design_region.txt

# After queue drains:
singularity exec $CONTAINER python3 boltz2_iterate_steering.py aggregate \
    --experiment-root ./resurface_pipeline_test/runs/design_N_vX

singularity exec $CONTAINER python3 boltz2_iterate_steering.py compute-final-metrics \
    --experiment-root ./resurface_pipeline_test/runs/design_N_vX \
    --rmsd-threshold 5.0 \
    --metric-column steered_ra_eff_vs_truth \
    --contact-cutoff 5.0

singularity exec $CONTAINER python3 extract_passing.py \
    --input ./resurface_pipeline_test/runs/design_N_vX/all_results_multicycle_with_metrics.csv
```

Note: `--n-cycles 1` now works correctly — reversion runs as part of cycle 0's
finalization. `--design-region-positions-file` is no longer needed on
compute-final-metrics (the protected-set machinery is gone).
