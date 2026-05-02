# Notes 9 — P0 audit: per-sequence aggregation, the multi-seed reversion bug, and one-shot post-processing

## Context

Notes 8 introduced multi-seed prediction via `--num-seeds N` and reversion
deduplication. Both worked at the predict-Boltz level, but neither was wired
through to the metric-aggregation or filtering layer. The pipeline still
reasoned in terms of per-seed rows: each unique sequence's N seeds appeared
as N independent rows in `all_results_multicycle.csv`, each with its own
contamination decision, ranking position, and (for contaminated rows) its
own reversion attempt.

This session was a P0 audit pass before the v5 NextFlow integration. The
audit refactored the pipeline to operate on **unique sequences** (median
+ spread across seeds) rather than per-seed rows, uncovered a real
correctness bug in the multi-seed reversion harvest that was silently
discarding 2/3 of the GPU work for `num_seeds=3` runs, and packaged the
manual four-step post-processing chain into a single submit script.

## Mental model shift

The pipeline now thinks in three stages, each with the same shape:

1. **Steering generation** — K unique steered sequences × `num_seeds`
   Boltz predictions = K · N raw per-seed rows.
2. **Per-sequence aggregation (steered)** — collapse N seeds per sequence
   into one summary record (median + MAD on continuous, majority vote on
   categorical). K aggregated rows. Contamination decision now made on the
   aggregate.
3. **Reversion + per-sequence aggregation (reverted)** — for each
   contaminated steered sequence, build the reverted sequence. Dedup
   across designs → R unique reverted sequences (R typically ≪ K because
   many revert to the same MPNN wild-type). Run R · N reverted predictions.
   Aggregate per unique reverted sequence. Map each aggregated reverted
   record back to every contaminated steered sequence that produced it.

Three output CSVs now correspond to these stages:

- `raw_per_seed_results.csv` — every seed's raw prediction (unchanged
  schema from `all_results_multicycle.csv`, just renamed for clarity).
- `aggregated_results.csv` — one row per unique steered sequence with
  median + MAD aggregates, majority votes, and verdict re-derived from
  the aggregate.
- `passing_summary.csv` — the filtered + ranked subset of (b) for wet-lab
  triage.

## What was built

### 1. Per-sequence aggregation helpers (`boltz2_iterate_steering.py`)

A new block of helpers between `_compute_unified_ranks` and
`_aggregate_csv_fieldnames` implements the aggregation primitives.

**Continuous metrics → median + MAD.** `_agg_continuous` returns
`(median, mad, n_used)` for a list of per-seed values. None / NaN / Inf
are dropped. MAD (median absolute deviation from the median) was chosen
over SD: with N=3 seeds, SD is dominated by whichever seed is the
outlier and gives wildly variable estimates run-to-run; MAD is the
robust counterpart that pairs naturally with the median. For N=1 the
MAD is reported as 0.0 (no spread visible). For N=2 it's the median of
the two absolute deviations (which equals half the range).

Continuous columns aggregated by median + MAD:

```
ra_eff_vs_truth, independent_receptor_rmsd, independent_effector_rmsd,
wrong_jaccard, true_jaccard, avg_plddt, complex_plddt, ptm, iptm,
pae_mean, ipae, pae_pass_frac, ipsae_ab, ipsae_ba, ipsae_min, actifptm
```

Integer-counted metrics (e.g. `n_contact_residues`) get the same
median + MAD treatment — they're discrete but ordered, and median across
N=3 still makes sense.

**Binary categoricals → majority vote.** `_agg_majority_binary` returns
`(majority, n_positive, n_used)` for 0/1 values. Threshold is
`≥ ceil(N/2)`. Currently used only for `receptor_intact`.

**Position-set categoricals → per-position majority.** This is the rule
for `mutated_contact_positions` and `contact_residues`. Each per-seed row
holds a comma-separated list of 1-based positions. The aggregator reads
all N lists, counts how many seeds reported each position, and the
aggregate set contains exactly those positions reported by `≥ ceil(N/2)`
seeds. Per-position counts are also exposed as a JSON map for
diagnostics. This means:

- A position contaminating all 3 of 3 seeds appears in the aggregate.
- A position contaminating only 1 of 3 seeds drops out.
- A position contaminating 2 of 3 seeds appears (matches the majority
  rule used elsewhere).

### 2. Canonical-seed picker

Each aggregated row needs to surface ONE PDB for diagnostics (ChimeraX
viewing, manual triage). The rule:

- **Odd N**: pick the seed at the median ra_eff (sort ascending, take
  index `(N-1)//2`).
- **Even N (degraded — only happens when some seeds failed)**: pick the
  WORSE of the two middle seeds (sort ascending, take index `N//2`).
  Set `canonical_seed_choice_warning =
  "even_n=N_chose_worse_of_two_middle_seeds"` so the choice is
  self-documenting.
- **N=1**: that seed, no warning.

The conservative even-N rule was a deliberate decision: if you're
visually evaluating a structure that's worse than the median you ranked
on, you'll trust the design at most as much as the structure shows —
which is the safer error direction. Picking the better seed and ranking
on the median would let visual triage overestimate quality.

The canonical seed's PDB is exposed as `canonical_pdb` in the aggregate
CSV.

### 3. Aggregated verdict + per-seed verdict counts

Two complementary pieces of verdict information now ride on each
aggregated row.

**`aggregated_verdict`** is re-derived from the aggregated reverted
metrics, not from the per-seed verdicts. The classification rule mirrors
`classify_reversion_verdict` but on the aggregate:

- `no_reversion` if the row had no reverted predictions (sequence wasn't
  flagged as contaminated by the steered aggregate).
- `pose_collapses` if `receptor_intact_majority != 1` or median ra_eff
  ≥ 5 Å.
- `new_contamination` if the per-position majority of
  `reverted_mutated_contact_positions` is non-empty (any position
  contaminating ≥ ceil(N/2) reverted seeds).
- `pose_holds` otherwise.

**Per-seed verdict counts** (`n_seeds_pose_holds`, `n_seeds_pose_collapses`,
`n_seeds_new_contamination`, `n_seeds_no_data`) are computed by running
the SINGLE-SEED verdict logic on each per-seed row independently and
counting outcomes. This was added at the user's request: when triaging
between two otherwise-similar designs, a 3/3 pose_holds is more
trustworthy than a 2/3 pose_holds even though both have the same
aggregated verdict. The verdict acts as the initial filter, the counts
as the tie-breaker.

### 4. New subcommand `aggregate-per-sequence` (`boltz2_iterate_steering.py`)

Run after `aggregate` and (ideally) `compute-final-metrics`:

```
python3 boltz2_iterate_steering.py aggregate-per-sequence \
    --experiment-root runs/design_3_v5
```

Reads `all_results_multicycle_with_metrics.csv` (or the unmetriced
fallback), groups rows by `(cycle, pathway, sequence_group)`, runs
aggregation on each group, picks the canonical seed, classifies the
verdict, recomputes ranks on the aggregated rows, and writes:

- `raw_per_seed_results.csv` (verbatim copy of input)
- `aggregated_results.csv`

Initial-baseline rows (cycle-0 wild-type predictions with empty
`sequence_group`) and any pre-multi-seed rows pass through as singletons
with `aggregated_verdict = "singleton"`.

Ranking is done on the aggregated rows: `rank_by_ra_eff` (ascending) and
`rank_by_composite_score` (descending). The composite is unchanged from
P0.3 — `true_jaccard − 0.05 · ra_eff_vs_truth` — but the values are now
medians.

### 5. The multi-seed reversion harvest bug

**The bug.** `harvest_reversion_results` in `reversion.py` had a latent
correctness bug for `num_seeds > 1`. The reversion plan correctly stages
`num_unique_reverted_sequences × num_seeds` Boltz predictions; the
manifest contains one entry per `(canonical_label, seed_index)` pair.
But the harvest loop wrote into `results[label]` keyed by `label` only,
so for each canonical_label, the loop overwrote the result on each of
its N iterations and only the **last seed**'s metrics survived. The
other seeds' GPU work was discarded silently.

This was masked in earlier testing because:

1. The notes 8 v2 test ("15 predictions, 5 × 3 seeds") confirmed 15
   GPU jobs were dispatched, but didn't audit whether all 15 results
   propagated through the harvest.
2. The downstream readers (`cmd_iterate_collect_finalize` etc.) accept
   a scalar reverted record per label and don't notice that they're
   getting one seed's data instead of an aggregate.

The fix preserves all per-seed reversion results and aligns each
steered per-seed row with its matching reverted seed.

**The fix.** Three-part change:

(a) **Contaminated entries carry `seed_index`.** `cmd_build_contaminated`
in `boltz2_iterate_steering.py` now reads `seed_index` from the candidate
dict (which carries it from the plan via the new backfill — see §6) and
includes it in each contaminated entry.

(b) **Manifest entries carry `all_label_seeds`.**
`reversion.write_reversion_plan` now records, for each unique reverted
sequence, a list `[{label, seed_index}, …]` mapping every steered design
that shares this reverted sequence to the seed_index it had on the
steered side. The legacy `all_labels` list-of-strings is kept for
back-compat.

(c) **Harvest stores per-seed records, then pairs by seed_index.**
`reversion.harvest_reversion_results` now maintains an internal
`per_seed_records[(canonical_label, seed_index)] = result_dict` map
during the inner loop instead of writing into a flat `results[label]`
dict. After the loop, each steered per-seed label is paired with the
matching reverted seed_index — `design_05_s0 ↔ rev_seed_0`,
`design_05_s1 ↔ rev_seed_1`, etc. The legacy `reversion_results.json`
shape (`{label: result_dict}`) is preserved so the three downstream
readers don't need changes — but each label now gets its own
seed-matched reverted record, not the canonical's last-seed result.

A new sidecar `reversion_results_per_seed.json` exposes the per-seed
structure (`{label: {str(seed_index): result_dict}}`) for any post-hoc
analysis.

**Pairing rule.** Steered seed S of unique sequence X gets reverted seed
S of X's reverted sequence. If steered seed 1 of X wasn't contaminated
(no reversion needed), reverted seed 1 of X's reverted sequence is still
computed but isn't paired with any X-row — it may be paired with a Y-row
if Y also reverted to the same sequence and Y's seed 1 was contaminated.
Truly wasted reverted predictions only occur when no contaminated steered
row at any seed_index uses a particular reverted seed_index — vanishingly
rare in practice.

**Implication for the per-seed CSV.** With the fix, each contaminated
per-seed steered row has its own matching reverted record, so the per-seed
CSV's `reverted_*` columns now contain N distinct values per sequence
(one per seed_index). The new aggregator can then median across them
meaningfully.

### 6. `sequence_group` and `seed_index` plumbing

Aggregation needs to know which per-seed rows belong to the same unique
sequence. The plan files (`plan.json`) wrote `sequence_group` and
`seed_index` into design entries already, but they didn't survive into
`all_results_multicycle.csv` for the cycle-N+1 path (cycle-0 read them
via DictReader from `steered_results.csv` so they were already fine).

Backfill applied at four sites in `boltz2_iterate_steering.py`:

- `cmd_compute_distances` (cycle-N+1 candidate dict construction):
  reads `design.get("sequence_group")` + `design.get("seed_index")`
  from the plan and adds them to the candidate dict.
- `cmd_kickoff_distances` (cycle-0 kickoff candidate construction): same
  treatment.
- Cycle-N+1 row construction in `cmd_aggregate`: reads from the
  candidate dict `c` and writes into the row.
- Cycle-0 reconstructed-initial-row block (when `steered_results.csv`
  is missing): writes `sequence_group: "initial"`, `seed_index: 0`.
  The initial baseline is treated as a 1-of-1 pseudo-group.

`_aggregate_csv_fieldnames` now lists `sequence_group` and `seed_index`
in the preferred fieldnames so they appear in the CSV in a sensible
position (right after `design`).

### 7. Rewritten `extract_passing.py`

`extract_passing.py` now reads `aggregated_results.csv` instead of
`all_results_multicycle_with_metrics.csv`. The output schema reflects
that per-seed values are gone and only the aggregates remain:

- All numeric columns are medians (`<metric>_median`).
- Two ranking-driver metrics (`ra_eff_vs_truth`, `true_jaccard`) get
  median + MAD + `n_used`. Spread on these directly affects how
  trustworthy the rank is.
- Other confidence metrics get only the median; full spread is in
  `aggregated_results.csv` if needed.
- New columns expose the canonical PDB and the seed-choice warning:
  `canonical_pdb`, `canonical_seed_index`,
  `canonical_seed_choice_warning`.
- Per-seed verdict counts visible: `n_seeds_pose_holds`,
  `n_seeds_pose_collapses`, `n_seeds_new_contamination`,
  `n_seeds_no_data`.
- `n_seeds` exposes how many seeds the aggregate is built from
  (matches `n_used` in the simple case; can be smaller if any seeds
  failed mid-pipeline).

The filter is now `aggregated_verdict in {"" / "no_reversion" /
"pose_holds"} AND rank_by_composite_score is populated`, mirroring the
old logic but on the aggregated record.

The confidence flag is recomputed from the median values:
`pae_pass_frac_median < 0.10 → low_pass_frac`, etc. The intent is
that a sequence is flagged when its TYPICAL prediction fails the
confidence filter, not when a single rogue seed does.

### 8. New `submit_postprocess.sh` — one-shot post-processing

Previously, after the main pipeline's SLURM chain drained, the user had
to manually run four steps on the submission node:

```
python3 boltz2_iterate_steering.py aggregate ...
singularity exec <container> python3 boltz2_iterate_steering.py compute-final-metrics ...
python3 boltz2_iterate_steering.py aggregate-per-sequence ...
python3 extract_passing.py ...
```

`submit_postprocess.sh` packages these four into a single SLURM CPU job
that switches between container (compute-final-metrics needs numpy +
Bio.Align) and host (everything else).

Required flag: `--experiment-root`. Defaults match the v5 command:
`--rmsd-threshold 5.0`, `--metric-column steered_ra_eff_vs_truth`,
`--contact-cutoff 5.0`, `--populate-all` on. Overrides:
`--rmsd-threshold`, `--metric-column`, `--contact-cutoff`,
`--no-populate-all`.

The script is **standalone, not auto-appended** to any other chain.
Submit it manually after the main queue drains. For multi-cycle
runs (`--n-cycles ≥ 2`), wait until the full cycle tree has finished
(no outstanding `b2is_*` / `b2ns_*` jobs) before invoking it.

### 9. Updated closing messages

`submit_boltz2_negative_steering.sh` and `submit_boltz2_iterate_steering.sh`
now print a clearer post-pipeline reminder pointing at
`submit_postprocess.sh` and the actual final outputs
(`aggregated_results.csv`, `passing_summary.csv`) rather than the
intermediate `all_results_multicycle.csv` that the old messages claimed
was the final result.

## Decisions explicitly made (and the rejected alternatives)

A few choices in this audit had non-obvious tradeoffs and were
discussed before being implemented; recording them so future-me knows
why:

**MAD over SD for spread.** With N=3 seeds, SD is dominated by single
outliers and gives wildly variable estimates run-to-run. MAD is the
robust counterpart to median.

**Per-position majority on `mutated_contact_positions`, not
union-on-contamination.** The conservative alternative was: if any
seed's reverted prediction shows a mutated contact, treat as
new_contamination. The chosen majority rule treats a single dissenting
seed as noise rather than evidence. Per-seed verdict counts surface
the dissent for tie-breaking.

**Aggregated verdict drives filtering, per-seed counts drive
tie-breaking.** Verdict is a hard yes/no; counts are the soft signal
for distinguishing 5/5 from 3/5 outcomes within the same verdict
category.

**Even-N canonical seed picks the WORSE of the two middle seeds.**
Conservative — if the structure shown is worse than the rank metric
suggests, the user underestimates rather than overestimates quality.
Plus a warning column so the choice is self-documenting.

**Steered-side aggregation collapses per-seed rows; per-seed CSV
preserves them.** The two coexist: `aggregated_results.csv` is what
ranking + filtering operates on; `raw_per_seed_results.csv` is the
audit trail with all the underlying values.

**Reversion seed pairing by `seed_index`, not by Boltz random seed.**
The Boltz random seeds on the steered (`base + design_idx + 1`) and
reverted (`base + 100000 + prediction_idx + 1`) sides are deliberately
distinct so the reverted prediction isn't coupled to its steered parent.
Pairing by `seed_index` is just a bookkeeping convenience — "the Nth
of N independent seeds requested on each side."

**`submit_postprocess.sh` is standalone, not auto-appended to the main
chain.** Auto-appending would require either restricting to
`--n-cycles == 1` runs or making `kickoff-finalize` /
`iterate-collect-finalize` themselves submit the post-processing job.
Both options trade complexity now for what NextFlow will replace anyway.
Standalone keeps the cycle-tree-completion gate as the user's
responsibility, which is correct for the manual-submission era.

## Bugs found and fixed

### Bug 1: multi-seed reversion harvest discarded 2/3 of GPU work

See §5 above. This was the central correctness issue uncovered in the
audit. Affected every run with `num_seeds > 1` that triggered any
reversion.

### Bug 2: `sequence_group` / `seed_index` not propagated through cycle-N+1

These columns existed in `plan.json` but were dropped at the
`cmd_compute_distances` → `cmd_aggregate` boundary because the row
construction in `cmd_aggregate` didn't pull them from the candidate
dict. Same issue at the kickoff path. Backfilled at four sites.

This wasn't a correctness bug in the legacy single-row-per-seed model
(nothing depended on these columns for filtering or ranking) but it
made per-sequence aggregation impossible without parsing design names.

### Pre-existing partial aggregator left in place

`boltz2_negative_steering.py:cmd_collect` already had a multi-seed
aggregate block that wrote `steered_results_aggregate.csv`. That code
only aggregated 3 RMSDs (using mean ± SD, not median + MAD), didn't
handle confidence metrics or contamination, and wasn't read by anything
downstream. Left in place as a harmless diagnostic — can be ripped out
in the NextFlow chat.

## Files audited but unchanged

Several files were reviewed and intentionally not modified:

- `compute_metrics.py` — clean for our purposes; pre-existing fragility
  around 1-based PDB residue numbering noted but pre-existing.
- `compare_versions.py` — no multi-seed implications, no changes needed.
- `boltz2_negative_steering.py` — earlier-session edits only; no audit
  changes this session.
- `derive_design_region.py` — pure preprocessing utility, 1-based /
  0-based handling is correct, no multi-seed implications.
- `multiseed_variance_test.py` — standalone diagnostic tool, two minor
  issues found (`INTACT_CUTOFF=3.0` stale per notes 8; non-standard
  quartile indexing for small N) but explicitly left alone since the
  tool isn't in the main pipeline path.
- `sequence_registry.py` — confirmed dead code (nothing imports it).
  Two minor schema quirks noted for if it ever gets wired up.
- `aggregate_results.sh` — cross-experiment batch summary tool. Will
  need a complete rewrite for NextFlow (read `passing_summary.csv` per
  experiment, take top N per experiment, re-rank globally, two-level
  identification). Deferred to NextFlow integration since the directory
  layout isn't settled.
- `submit_steering_experiments.sh` — wrapper that loops over a hardcoded
  EXPERIMENTS array submitting one steering mode per experiment. Will
  be replaced by NextFlow's parameter sweep.
- `submit_final_validation.sh` — post-pipeline tool that re-validates
  survivors with multiple seeds. The "option E" use case (run main
  pipeline with `--num-seeds 1`, validate after) becomes irrelevant
  when `num_seeds ≥ 3` is the v5 default.
- `submit_multiseed_variance_test_slurm.sh` — wrapper for the
  diagnostic variance-test tool, used ad-hoc not in the main pipeline.

## Synthetic end-to-end verification

A synthetic per-seed CSV (`/tmp/p0test/build_synth.py`) was built to
exercise every code path in the new aggregator. Six scenarios:

- Clean winner (3 tightly-clustered seeds): ranks #1, MAD ≈ 0.1.
- All-pose_holds (3/3 seeds passing): ranks correctly, n_seeds_pose_holds=3.
- Outlier seed at 12 Å (2 of 3 intact): median is robust to the outlier
  (2.5 Å, MAD 0.1), majority intact = 1.
- Mixed verdicts (2 pose_holds, 1 new_contamination, dissenting position
  in only 1 seed): aggregated verdict = pose_holds (position fails
  per-position majority), per-seed counts (2/0/1) preserved for
  tie-breaking.
- Degraded N=2 (one seed failed): even-N warning fires correctly,
  canonical seed = worse of the two.
- All-new_contamination (3/3 seeds report position 44): aggregated
  verdict = new_contamination, position 44 survives majority, no rank
  assigned, correctly excluded from passing_summary.

Plus the cycle-0 initial baseline as a singleton passthrough.

`extract_passing.py` produced the expected 5 passing rows in correct
ranking order with the correct columns populated. The MAD column
revealed the wide-spread degraded case (`design_04`: MAD 0.3 with
n_used=2). Per-seed verdict counts let the user distinguish 3/3 wins
from 2/3 wins within the pose_holds bucket.

## Files changed

| File | What |
|------|------|
| `boltz2_iterate_steering.py` | New aggregation helpers (median+MAD, majority vote, position-set majority). New `aggregate-per-sequence` subcommand. Backfill of `sequence_group` + `seed_index` into 4 row-construction sites. Contaminated entries carry `seed_index`. |
| `reversion.py` | Per-seed harvest fix: `per_seed_records[(canonical_label, seed_index)]` instead of `results[label]`. New `all_label_seeds` field in manifest entries. Pairing by matching `seed_index`. New `reversion_results_per_seed.json` sidecar. |
| `extract_passing.py` | Rewritten for new schema: reads `aggregated_results.csv`, surfaces medians + MAD + per-seed verdict counts + canonical PDB + warnings. |
| `submit_postprocess.sh` | NEW. Single SLURM CPU job running aggregate → compute-final-metrics → aggregate-per-sequence → extract_passing. Standalone, manual submission after the main queue drains. |
| `submit_boltz2_negative_steering.sh` | Closing message updated to point at `submit_postprocess.sh` and the new outputs (`aggregated_results.csv`, `passing_summary.csv`). |
| `submit_boltz2_iterate_steering.sh` | Closing message updated with the same `submit_postprocess.sh` reminder. |

## v5 command (unchanged) + post-processing

```bash
./submit_boltz2_negative_steering.sh \
    --ground-truth resurface_pipeline_test/design_3.pdb \
    --receptor A --effector B \
    --receptor-fasta resurface_pipeline_test/design_3_seq_3_receptor.fasta \
    --workdir runs/design_3_dedup_test_v5 \
    --mode mild --max-mutations 3 --candidate-pool-size 12 \
    --n-designs 30 --num-seeds 3 --n-cycles 1 \
    --diffusion-samples 5 \
    --protected-set-source design_region_union \
    --true-interface-indices-file resurface_pipeline_test/design_3_true_interface.txt \
    --design-region-indices-file resurface_pipeline_test/design_3_design_region.txt
```

After the queue fully drains:

```bash
./submit_postprocess.sh --experiment-root runs/design_3_dedup_test_v5
```

That single submission runs all four post-processing phases as one
SLURM CPU job and produces:

- `all_results_multicycle.csv`
- `all_results_multicycle_with_metrics.csv`
- `raw_per_seed_results.csv`
- `aggregated_results.csv`
- `passing_summary.csv`

## Outstanding work (next chat)

1. **Diagnose whether this scaffold is a dead end.** Three consecutive
   live runs (v5, v6, v7) produced zero genuine passing designs.  The
   bugs uncovered along the way (sections below) were real and worth
   fixing, but fixing them only made the result LESS optimistic, not
   more: v5 and v7 had false-positive pose_holds that the fixes
   correctly eliminated.  This suggests the problem is not a pipeline
   bug at this point — it's the biology.  See "Outlook and next
   direction" at the end of this file.

2. **NextFlow integration.** The original next-chat plan, deferred from
   this session. Will subsume:
   - `submit_steering_experiments.sh` (parameter sweep → NextFlow channels)
   - `submit_final_validation.sh` (becomes a NextFlow process or removed
     if `num_seeds ≥ 3` is always the default)
   - `aggregate_results.sh` (rewrite as a Python tool that reads
     `passing_summary.csv` per experiment, takes top N per experiment,
     re-ranks globally; two-level identification by
     `rfdiffusion_design` + `mpnn_seq`; pure Python, no SBATCH
     scaffolding).

3. **Optional cleanup of the legacy partial aggregator** in
   `boltz2_negative_steering.py:cmd_collect`. The
   `steered_results_aggregate.csv` writer is dead code now that
   `aggregate-per-sequence` does the real aggregation. Safe to delete
   in the NextFlow chat.

4. **Minor cleanup**: when `extract_passing.py` finds zero passing
   rows, it should still write the CSV header (currently produces an
   empty file with 0 columns, which trips downstream tools).

---

# Addendum — v5 / v6 / v7 live shakedown and the three follow-up bugs

The audit section above covered design and synthetic verification.
This addendum covers what happened when that code was exercised
against real GPU-produced data, and the three additional bugs that
only surfaced on the full postprocess chain.

## v5 — shakedown: 30 designs × 3 seeds, 3 muts / pool 12

First live run.  Produced the expected output shapes (91 raw rows →
31 aggregated rows → 1 passing: design_01).  The multi-seed reversion
fix (§5 of main notes) was visibly working: design_01's three seeds
had distinct reverted_ra_eff values (4.276, 4.512, 4.424), confirming
the per-seed harvest pairing; under the old bug all three would have
been identical.

However, on closer inspection, design_01 was a false positive.  Seeds
0 and 2 had empty `reverted_mutated_contact_positions` in the
per-seed CSV but `reversion_verdict='new_contamination'` (set by the
harvest's authoritative classification).  The aggregator's
`_per_seed_verdict_breakdown` was re-deriving verdicts from the
blanked column, which made it count 3/3 as pose_holds when the
authoritative count was 1/3.  This cascaded into the aggregated
verdict being pose_holds instead of new_contamination.

### Bug A — `_populate_reverted_mutations` blanked contact columns for non-pose_holds rows

The legacy display choice: for dropped rows (pose_collapses,
new_contamination) the aggregate CSV showed blank `reverted_*` contact
columns on the principle that dropped rows "don't need this data."
Pre-audit that was harmless because single-seed rows had a direct 1:1
mapping to their verdict.  Post-audit the aggregator needs the
unblanked per-seed contamination data to compute the
per-position majority.

**Fix**: pulled the reverted contact population out of the
pose_holds-only branch of `_populate_reverted_mutations`.  Now runs
for ANY row with a non-empty `reversion_verdict`.
`reverted_mutated_contact_positions`, `reverted_contact_residues`,
`reverted_n_contact_residues`, `reverted_contact_cutoff_used` all
populate for pose_collapses and new_contamination rows too.  The
pose_holds-only branch retains its exclusive handling of
`reverted_mutations_aa` / `reverted_mutations_chimerax` /
`reverted_total_mutations` (those genuinely only make sense for
rows that continue to the next cycle).

### Bug B — `_per_seed_verdict_breakdown` re-derived instead of trusting the authoritative column

Same root cause: the helper ran the single-seed verdict logic on the
raw CSV columns, which saw the blanked contact fields and classified
as pose_holds.

**Fix**: `_per_seed_verdict_breakdown` now prefers the authoritative
`reversion_verdict` column (set by `classify_reversion_verdict` at
harvest time, before any display blanking) and only falls back to
re-derivation when the column is empty.  Maps
`pose_holds_not_intact` / `missing_verdict` / `unknown` all to
`pose_collapses` for the count.

### Bug C — `compute-final-metrics --populate-all` still excluded dropped rows

Separate from the display-blanking issue: the `_passes` filter in
`cmd_compute_final_metrics` skipped rows with `reversion_dropped == 1`
even when `--populate-all` was set.  Consequence: 33 rows in v5 had
`steered_ra_eff < 5`, `steered_receptor_intact == 1`, but blank
`steered_ipsae_min` / `_iptm` / `_complex_plddt` etc.  With per-seed
aggregation this biases the aggregate — the median across seeds is
wrong if some seeds' confidence metrics are missing entirely.

**Fix**: `--populate-all` now truly means "populate every row" —
the `reversion_dropped` gate is bypassed when the flag is set.
Without the flag, the default behaviour (skip dropped rows) is
preserved.

## v6 — 50 designs × 3 seeds, 4 muts / pool 8

With bugs A, B, C fixed, v6 reprocessed cleanly.  150 per-seed rows
→ 51 aggregated rows (50 unique steered + 1 initial singleton).
Aggregated verdict tally:

- 21 no_reversion (clean steered; mostly effector at wrong site,
  ra_eff ≈ 28 Å)
- 25 pose_collapses (steered ra_eff < 5 but reverted collapses to
  ra_eff ≈ 28 Å — rescue depended entirely on mutations)
- 4 new_contamination (design_00, design_17, design_25, design_36 —
  both steered AND reverted ra_eff < 5, but persistent mutation
  contacts at positions 22, 26, 32 or 5, 32)
- 0 passing

Of 50 unique steered sequences, 29 triggered reversion → 26 unique
reverted sequences after dedup.  Size distribution of reverted
sequences: 1 MPNN wild-type, 7 with 1 remaining mut, 11 with 2, 7
with 3.  **Most reverted sequences retain 1-3 mutations, not
collapsing to wild-type** — refuting a prior guess that low
`max_mutations` would collapse the reverted space.

### Bug D — `compute-final-metrics` read the pre-rename column

Noticed during v6 inspection that `steered_mutated_contact_positions`
was empty and `steered_n_contacts_on_mutated_positions = 0` for all
73 reverted rows, despite obvious intersections between the
steered mutation positions and steered contact residues (e.g.
design_02_s1 mutations `{5, 22, 26, 32}` and contacts including
`{5, 32}` — a clear intersection of `{5, 32}`).

Root cause: `cmd_compute_final_metrics` at line 4941 read
`row.get("mutations_chimerax", "")` — the **pre-rename** column name.
After `_translate_aggregate_row` has run in `cmd_aggregate`, the
column is named `steered_mutations_chimerax`.  Consequence:
`compute_metrics.py` was always called without `--mutated-positions`,
could never compute the intersection, and blanked
`mutated_contact_positions` / `n_contacts_on_mutated_positions` on
every row.  Worse, these blanks then **overwrote** the correct
values that `_populate_reverted_mutations` had written during
`cmd_aggregate` (because compute-final-metrics runs after aggregate).

**Fix**: `mutations_cx = row.get("steered_mutations_chimerax", "") or
row.get("mutations_chimerax", "") or ""`.  Reads the post-rename key
with fallback to the legacy name.

Why synthetic tests didn't catch this: they fed pre-built per-seed
CSVs directly into `aggregate-per-sequence`, bypassing both
`cmd_aggregate` and `cmd_compute_final_metrics`.  The full four-step
chain had never been exercised end-to-end against real data.

After fix: 146/151 rows correctly populated
`steered_mutated_contact_positions` (the 5 blanks are the initial
baseline and genuine no-contact rows).  `steered_n_contacts_on_mutated_positions`
values now meaningfully distributed from 1 to 4.

### New `--force-recompute` flag on `submit_postprocess.sh`

Fixing bug D required re-running `compute-final-metrics` on existing
v6 output.  But its default is `--skip-existing`, which would reuse
the stale (buggy) cached rows.  Added `--force-recompute` flag that
passes `--no-skip-existing` downstream.  Use when the postprocess
code has changed and cached rows are stale.

Usage: `./submit_postprocess.sh --experiment-root <root> --force-recompute`

## v7 — 50 designs × 3 seeds, 6 muts / pool 10

Rationale: after v6 showed reverted sequences keeping 1-3 mutations,
tried the OPPOSITE direction from v5 — more steering mutations, so
that contaminating ones can be reverted while leaving several
non-contaminating ones to hold the pose.

Initial v7 postprocess (before bug E fix) surfaced 2 passing
designs: design_08 and design_13.  Per-seed inspection showed both
were false positives:

- **design_08**: 1 pose_holds, 1 pose_collapses, 1 new_contamination
  (seeds fail in three different ways)
- **design_13**: 0 pose_holds, 1 pose_collapses, 2 new_contamination
  (no single seed says pose_holds; reverted_ra_eff MAD = 1.844
  indicates massive per-seed disagreement, not convergence)

### Bug E — aggregated verdict didn't require per-seed agreement

Root cause: `_classify_aggregated_verdict` checked three things:

1. `reverted_receptor_intact_majority == 1`
2. `reverted_ra_eff_vs_truth_median < 5`
3. `reverted_mutated_contact_positions_majority == ""` (no position
   reached majority contamination)

and if all three passed, returned pose_holds.  But condition 3 has a
subtle gap: if each seed contaminates at a DIFFERENT position
(seed 1 at 32, seed 2 at 46, seed 3 at 5), no position reaches
majority, so the majority set is empty — the check passes — but EVERY
seed individually said new_contamination.  Similarly, if one seed
pose_collapses and two contaminate at different positions, the
aggregate said pose_holds even though 0/3 seeds individually did.

**Fix**: added a per-seed-majority gate.  Aggregated pose_holds now
requires `n_seeds_pose_holds >= ceil(N/2)` in addition to the three
aggregated checks.  When the aggregated rules pass but per-seed
disagreement is too high, the verdict is downgraded to whichever
failure mode has plurality — with ties preferring the more
conservative pose_collapses.

Under the new rule:
- design_08 (1/1/1 tie) → pose_collapses (conservative tiebreak)
- design_13 (0/1/2) → new_contamination (plurality)
- 3/3 pose_holds rows (real wins) → pose_holds (unchanged)

After fix: v7 re-postprocessed with 0 passing rows, matching v5 and
v6.  Verified post-redeployment: design_08 and design_13 both
correctly downgraded (design_08 → pose_collapses, design_13 →
new_contamination), and `passing_summary.csv` contains 0 rows.
`aggregated_verdict_reason` now reads "only X/3 seeds individually
agree" on the downgraded rows — confirming the new logic ran.

**A note on deployment gotcha**: after editing
`boltz2_iterate_steering.py` locally and re-running
`submit_postprocess.sh`, the passing_summary kept showing the old
two false positives.  Root cause: only `submit_postprocess.sh` had
been redeployed to the cluster, not `boltz2_iterate_steering.py`.
Both files need to ship together when the fix is in the Python.
Diagnostic: the old `aggregated_verdict_reason` phrasing ("reverted
pose holds by majority ... median ra_eff X Å") lacks the new
"seeds individually pose_holds" clause, so you can tell by string
match which version ran.

## Confirmed final state: three runs, zero passing

With all five bugs (A/B/C/D/E) fixed and verified by re-postprocess,
the final tally across all three live shakedown runs:

| Run | Config | Unique seqs | Genuine pose_holds |
|-----|--------|-------------|--------------------|
| v5  | 3 muts / pool 12 / 30 designs | 30 | 0 |
| v6  | 4 muts / pool 8 / 50 designs  | 50 | 0 |
| v7  | 6 muts / pool 10 / 50 designs | 50 | 0 |

130 unique steered sequences tried across three parameter regimes;
zero produced a reverted pose that held by genuine per-seed
majority.  The pipeline is now considered correct: the zero-survivor
result is the target's actual signal, not an artifact of a bug.

## 100-seed variance baseline on the unmodified sequence (NEW)

Ran `multiseed_variance_test.py` on the unmodified (no-substitutions)
baseline sequence for design_3 with 100 seeds.  The result is
striking and reshapes the outlook below:

| Metric | Distribution |
|--------|--------------|
| ra_eff_vs_truth | min 28.43 Å, max 31.91 Å, median 29.64 Å, std 0.72 Å |
| Seeds with ra_eff < 5 Å | 0 / 100 |
| Seeds with ra_eff < 15 Å | 0 / 100 |
| Seeds with receptor intact | 100 / 100 |
| iptm | mean 0.89, min 0.81, max 0.93 |

Interpretation:

1. Boltz2 does NOT stochastically explore the pose space on this
   scaffold.  All 100 seeds converge to essentially the same
   off-target pose with very low variance (std 0.72 Å over a 3.5 Å
   range).
2. There is NO lucky-seed phenomenon to exploit.  Zero seeds in 100
   land anywhere near the correct interface.
3. iptm is consistently high — the predictor is confident about the
   wrong answer.

Implications for the three shakedown runs:

- The 25 "pose_collapses" rows in v6 make sense: steering pulls the
  prediction to ra_eff < 5, removing substitutions snaps it back to
  the ~29 Å off-target mode — which is the only mode the predictor
  knows on this scaffold.
- Any steered prediction at ra_eff < 5 is ENTIRELY driven by the
  substitutions.  There's no underlying landscape the substitutions
  are "unmasking."
- The v5 diagnostic priority of "variance baseline" is DONE and
  answered.  No need to re-run.

## Summary of files changed in this addendum

| File | What changed |
|------|-------------|
| `boltz2_iterate_steering.py` | Bugs A, B, C, D, E all fixed.  `_populate_reverted_mutations` now populates reverted contact data for all reverted rows.  `_per_seed_verdict_breakdown` prefers the authoritative verdict column.  `cmd_compute_final_metrics` bypasses the reversion_dropped gate when `--populate-all` is set and reads the post-rename `steered_mutations_chimerax` column.  `_classify_aggregated_verdict` adds a per-seed-majority gate for pose_holds. |
| `submit_postprocess.sh` | New `--force-recompute` flag that passes `--no-skip-existing` to compute-final-metrics.  Default behaviour unchanged. |

No other files changed in this addendum.

## Outlook and next direction (revised post-variance-baseline)

Three runs (v5, v6, v7) across 30-50 designs and 3-6 substitution
counts produced zero genuine pose_holds survivors after the bugs
were fixed.  Combined with the 100-seed variance baseline result,
the picture is now clearer: **Boltz2 has a strong off-target
attractor at ra_eff ≈ 29 Å on this scaffold and no access to the
correct interface via random sampling**.  Substitutions can push
the prediction away from the attractor transiently, but removing
the substitutions releases the prediction back to it.

The next chat should skip the variance-baseline step (already done)
and prioritize differently:

1. **Forced-anchor sanity check** — HIGHEST PRIORITY.  Run Boltz on
   the unmodified sequence with a template/anchor forcing the
   partner chain at the true interface.  Does it stay there, or
   drift back to the ~29 Å attractor?  This distinguishes two cases
   very cleanly:
   - If it STAYS: the correct pose exists in Boltz's landscape but
     is inaccessible from the baseline.  Smarter steering (better
     anchor positions, different substitution patterns) could
     plausibly find it.
   - If it DRIFTS: the predictor cannot represent the target state
     on this scaffold at all.  No steering workflow will succeed
     here; move to a different candidate design.
2. **Structural inspection in ChimeraX** of v6 near-hits (design_00,
   17, 25, 36).  At steered ra_eff < 5 Å, is the partner chain
   actually at the correct interface, or at a third site
   coincidentally within 5 Å?  These were the closest the workflow
   has come to passing; understanding their structure is the
   strongest diagnostic for whether the substitutions are doing
   something real.

If the forced-anchor check PASSES (Boltz can hold the correct pose
when shown it):

- **Exhaustive small-substitution search** instead of random sampling.
  For max-subs=2 over pool=8 there are only 28 combinations — fully
  tractable and covers the space.
- **Position stratification** of the candidate pool: bias away from
  the recurring-contamination positions (22, 26, 32) and toward
  positions that produced v6 near-hits.
- **Different steering modes** beyond "mild": conservative or alanine
  probe different failure modes.
- **Relaxed threshold temporarily** (`--rmsd-threshold 8.0`) so
  near-misses reach reversion — confirms whether the issue is pose
  acquisition or pose stability.

If the forced-anchor check FAILS: switch to a different candidate
design whose baseline variance distribution actually places some
seeds near the correct interface, or whose scaffold is more
amenable to Boltz prediction.  No amount of steering will fix a
scaffold the predictor cannot represent.

**Longer-term orthogonal approach**: train a fine-tuned ESM-2 on
in-planta assay data as an independent scoring layer.  The current
pipeline uses Boltz for both pose geometry AND ranking; decoupling
structure from function scoring could surface designs that Boltz
mis-poses but that are functionally plausible.  Requires assay data
accumulation before feasible.

**`--n-cycles 2` is NOT on the table.**  Cycle 1 refines cycle-0
survivors; with zero survivors it produces zero output.  Only
useful after cycle 0 yields something.
