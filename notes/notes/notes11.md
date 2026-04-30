# Notes 11 — Gated `new_contamination`, ranking method, passing_summary trim, and the v11 wet-lab triage

## Context

Notes 10 closed with v10 in flight: 19 per-seed `pose_holds`, 6 aggregated
promoted to wet-lab triage. The strict-protection patch (v11) was staged
but not yet validated. v11 was expected to remove design-region
contamination entirely — and it did — but it surfaced a deeper bug
in the post-reversion contamination check that was hiding the majority of
the pipeline's actual successes.

This session covered:

1. v11 deployment and validation of the strict-protection patch.
2. Discovery and fix of the **gated `new_contamination` bug** — the
   contamination check on the reverted prediction was not applying the
   same gating policy as the pre-reversion check, falsely flagging 61 of
   167 reversion runs as contamination at off-interface positions.
3. Establishing the **tier-then-composite ranking method** for wet-lab
   selection.
4. Trimming `passing_summary.csv` from 44 → 43 columns: removed 8
   redundant/tautological columns, added 7 mutation-set columns
   (steered, reverted majority, reverted union).
5. End-to-end documentation of all output CSVs in `columns.md`.

## v11 strict-protection results — the patch worked

The v11 plan.json confirmed strict protection was active:

```
protected_set_source: design_region_union
max_mutations_effective: 5
max_mutations_reduced_reason: "pool_size<6 (strict protection preserved)"
candidate_pool: {5, 7, 22, 26, 32}
protected_set_idx: 26 positions (true interface ∪ design region)
```

Design-region positions 36, 44, 46 were NOT in the candidate pool, NOT
mutated, NOT contaminating any prediction. The pre-reversion side of the
contamination story was finally clean.

Per-seed verdict distribution (before the gated-classifier patch):

| Verdict           | Count |
|-------------------|-------|
| pose_holds        | 29    |
| pose_collapses    | 77    |
| new_contamination | 61    |

A 50% jump in pose_holds vs v10 (29 vs 19), but the new_contamination
count seemed high — and on inspection every single one was at off-
interface positions.

## The gated `new_contamination` bug

### What was happening

`classify_reversion_verdict` in `reversion.py` was firing
`new_contamination` whenever the reverted prediction had any mutated
residue contacting the effector — regardless of whether the position
was inside the design region or true interface. This was inconsistent
with the pre-reversion contamination policy: that side already excluded
contacts at off-interface positions (the "reaching bit" pattern at
positions 5, 7, 22, 26, 32).

The asymmetry meant a design could:

1. Have the pre-reversion check accept off-interface contacts at, say,
   positions 22/26 as valid steering work (don't revert those).
2. Then have the post-reversion check fire `new_contamination` on the
   exact same off-interface positions, dropping the design.

In v11, **all 61 of the 61 `new_contamination` seeds** were flagged at
off-interface positions only. Position breakdown:

```
pos 5:  14 seed-hits  [off-iface]
pos 7:  13 seed-hits  [off-iface]
pos 22: 33 seed-hits  [off-iface]
pos 26: 23 seed-hits  [off-iface]
pos 32:  8 seed-hits  [off-iface]
```

Zero design-region or true-interface hits across the entire 61-seed set.
The `new_contamination` verdict was firing exclusively on contacts the
pipeline had already decided to ignore.

### The fix

Two-file patch:

- `reversion.py`: `classify_reversion_verdict` now accepts
  `contamination_gating_positions: Optional[Set[int]] = None`. When
  provided, mutated contacts at positions outside the gating set are
  reported in the verdict reason (for transparency) but don't change the
  verdict from `pose_holds`. When `None`, legacy ungated behaviour
  prevails.
- `boltz2_iterate_steering.py`:
  - `cmd_harvest_reversions` builds the gating set from `plan.json`'s
    `design_region_idx ∪ true_interface_idx` (0-based → 1-based) and
    passes it to `classify_reversion_verdict`.
  - `_classify_aggregated_verdict` accepts the same gating set and
    filters `reverted_mutated_contact_positions_majority` against it
    before firing the aggregated `new_contamination` verdict. The caller
    in `cmd_aggregate_per_sequence` loads the gating set from
    `cycle_0/plan.json`.

### Result

Per-seed verdict distribution after the patch:

| Verdict           | Count | Δ from before |
|-------------------|-------|---------------|
| pose_holds        | 90    | +61           |
| pose_collapses    | 77    | 0             |
| new_contamination | 0     | −61           |

Aggregated verdict distribution:

| Verdict           | Count |
|-------------------|-------|
| pose_holds        | 34    |
| no_reversion      | 41    |
| pose_collapses    | 25    |
| singleton         | 1     |

`passing_summary.csv` went from 18 rows (10 pose_holds + 8 no_reversion)
to **42 rows (34 pose_holds + 8 no_reversion)**. Three pose_holds
designs were under composite rank 5; the strongest are tied across true
interface overlap and confidence.

### Why this was easy to miss

The pre-reversion gating was implemented as an effector-atom filter
(restricting which effector atoms count as "interface"), and the post-
reversion check inherited that filter via the same `compute_metrics.py`
invocation — so the filter WAS applied, but at the **atom** level. The
position-level contamination decision (`reverted_mutated_contact_positions
= remaining_muts ∩ contact_residues`) didn't apply any further gating
on which receptor positions counted. The atom filter restricts which
effector atoms register as contacts in the first place, but does not
restrict which receptor positions are considered. So contacts to
ungated receptor positions (5, 7, 22, 26, 32) involving filtered
effector atoms still came through as `mutated_contact_positions`, and
the unconditional `if mutated_contacts: → new_contamination` on line
929 of `reversion.py` then fired.

The fix gates at the receptor-position level, mirroring the pre-
reversion policy.

## Ranking method — tier then composite

The pipeline emits two rank columns (`rank_by_composite_score`,
`rank_by_ra_eff`) but neither incorporates per-seed agreement. With
v11's 90 pose_holds seeds spread across 34 sequences, many designs are
3/3, 2/3, or 1/3 pose_holds. Reproducibility across independent Boltz
seeds is a stronger transferability signal than a small interface-
overlap difference.

Adopted policy:

| Tier | Criterion                                                              |
|------|------------------------------------------------------------------------|
| A    | `aggregated_verdict == no_reversion` OR `n_seeds_pose_holds == n_seeds` (3/3) |
| B    | `1 < n_seeds_pose_holds < n_seeds` (2/3 of 3 seeds)                    |
| C    | `n_seeds_pose_holds == 1` (1/3)                                        |

Within each tier, rank by `rank_by_composite_score`. Take the best
non-empty tier. Same rule applied per MPNN sequence when comparing
across the 100 candidates.

Confidence metrics (`iptm_median`, `actifptm_median`, `ipsae_min_median`)
are tie-breakers for near-identical composite scores within a tier —
not the primary ranker, because prior analysis (notes 10 / metric
reliability memo) showed iPTM and ipSAE have anti-informative AUROC for
correct-vs-wrong discrimination on this system.

This rule is recomputable from the columns already in
`passing_summary.csv` (`n_seeds_pose_holds`, `n_seeds`,
`aggregated_verdict`, `rank_by_composite_score`); it's NOT written
back as a `tier` column to avoid baking the policy into the data
files. If the policy evolves the data files don't need to be
regenerated.

## v11 top candidates

Under tier-then-composite (top three of tier A):

**design_94 (no_reversion, 0/3 pose_holds — interpreted as 3/3 clean
steered)**: composite rank 1. Steered ra_eff 2.71 Å (median across 3
seeds), MAD 0.29. iPTM 0.892, actifPTM 0.758, ipSAE 0.754,
true_jaccard 0.720. Mutation at position 7 contacts the effector in all
three steered seeds — load-bearing steering at a non-design-region
position. The wet-lab construct is the unmutated MPNN sequence.

**design_32 (pose_holds, 3/3)**: composite rank 3. Reverted ra_eff
3.53 Å, true_jaccard 0.739 (highest in the set). iPTM 0.845, actifPTM
0.533, ipSAE 0.526. All three seeds reverted to the same set of
remaining mutations.

**design_68 (pose_holds, 3/3)**: composite rank 4. Reverted ra_eff
3.28 Å, MAD 0.01 (tightest of the top tier). iPTM 0.782, actifPTM 0.441,
ipSAE 0.429, true_jaccard 0.708.

design_15 (composite rank 2 by composite alone) drops to tier B because
it's 2/3 pose_holds — the third seed produced an off-interface-only
contamination flag pre-patch, which under the new gating reclassifies
to pose_holds in 2/3 of seeds.

## passing_summary.csv trim and mutation columns

### Removed (8 columns)

These were redundant or tautological:

- `ra_eff_vs_truth_mad`, `ra_eff_vs_truth_n_used` (MAD in
  aggregated_results)
- `true_jaccard_mad`, `true_jaccard_n_used` (MAD in aggregated_results)
- `receptor_intact_majority`, `receptor_intact_n_positive`,
  `receptor_intact_n_used` (tautologically 1 for any row that reaches
  passing_summary)
- `reverted` (determinable from `aggregated_verdict`)

### Added (7 mutation columns)

The wet-lab construct never carries any steering mutations, but the
scored prediction does — the user wanted to see exactly which mutations
are in each scored sequence. Per-seed mutation strings were already in
`raw_per_seed_results.csv`, but `passing_summary.csv` only had a
median count.

For each unique sequence's three seeds:

- `steered_mutations_chimerax` / `steered_mutations_aa` — full set,
  identical across seeds of one sequence_group (same sequence by
  definition).
- `reverted_mutations_majority_chimerax` / `_aa` — positions remaining
  in **more than half** the seeds (for n=3, ≥2 seeds).
- `reverted_mutations_all_chimerax` / `_aa` — **union** of remaining
  positions across all seeds.
- `reverted_mutations_seeds_identical` — `1`/`0`/blank flag indicating
  whether all seeds had the same reverted set.

`extract_passing.py` now requires `--per-seed
raw_per_seed_results.csv` (defaults to same dir as `--input`) to
populate these.

### Mutation column AA-lookup bug fix

First implementation of the AA helper looked up identities from the
**first seed** that had non-empty data. For the union-across-seeds
case, this missed positions retained only in non-first seeds. Fixed
by aggregating AA tokens across all seeds into a position-to-token
map, with steered_mutations_aa as a fallback (since steered always
contains every mutated position).

## A side observation: pool uniformity in v11

Every one of the 100 designs uses the exact same set of mutated
positions: `{5, 7, 22, 26, 32}`. This is a consequence of the strict-
protection patch reducing `max_mutations_effective` to 5 (the pool
size). With pool size = max mutations, every design HAS to mutate the
full pool — only AA-identity choice differs between designs.

For future runs, if positional diversity across designs matters
(e.g. for the assay's ability to sample sequence-binding mappings),
either the candidate pool needs to be larger (loosen protection
beyond design-region union) or `max_mutations` needs to be set
strictly less than the pool size. The latter would relax steering
pressure proportionally — worth thinking about whether the resulting
designs are still steered enough.

## Output file final state

| File                                  | Cols | Role                                                         |
|---------------------------------------|------|--------------------------------------------------------------|
| `all_results_multicycle.csv`          | ~68  | Per-seed, written by `aggregate`. Confidence cols blank.     |
| `all_results_multicycle_with_metrics.csv` | ~68  | Same, with confidence cols populated by `compute-final-metrics`. |
| `raw_per_seed_results.csv`            | ~68  | Renamed copy by `aggregate-per-sequence`.                    |
| `aggregated_results.csv`              | ~230 | Per unique sequence: median + MAD + n_used per metric, position majorities. |
| `passing_summary.csv`                 | 43   | Wet-lab triage: filtered + ranked + mutation summaries.      |

Full schema reference is in `columns.md`.

## Files patched / created this session

- `reversion.py` — gated `classify_reversion_verdict`.
- `boltz2_iterate_steering.py` — gated `cmd_harvest_reversions` and
  `_classify_aggregated_verdict`, with caller updates in
  `cmd_aggregate_per_sequence`.
- `extract_passing.py` — column trim, per-seed input, mutation
  summarisation helpers.
- `columns.md` — full schema documentation for all five output CSVs,
  ranking method, gating semantics, per-seed-vs-aggregated verdict
  reconciliation, n_used redundancy note.
- `notes11.md` — this file.

## On the horizon (carried from notes 10)

1. **Multi-seed validation of passing candidates** — promoted
   sequences through additional independent seeds before wet-lab.
   Now more pressing because the v11 `pose_holds` set has 34 entries,
   and many are 2/3 or 1/3 — a confirmation pass would tighten the
   ranking signal.

2. **Reversion deduplication across designs** — when many designs
   revert to the same residual mutation set, the corresponding Boltz
   reversion predictions are redundant. notes 8 / 9 partially addressed
   this for steered sequences; the same logic on the reverted side
   would save GPU.

3. **Generative ESM-2 for steering proposals** — fine-tune ESM-2 (LoRA
   or full fine-tune, depending on cost tolerance) on production
   NextFlow data to generate candidate steering mutation sets given a
   (receptor, effector) pair. Boltz remains the gate; the model is a
   pre-filter to avoid wasting GPU on bad mutation sets. Training data
   accumulates as a free byproduct of production runs once the
   classifier semantics are stable. Discussed under "fool Boltz" framing
   in this session — well-suited as a generator since the goal is to
   produce mutation sets that pass Boltz, not to predict in-planta
   binding directly.

4. **Orthogonal scoring layer (separate ESM-2)** — once PCP assay data
   is available, train a second model to predict in-planta binding from
   (receptor sequence, effector sequence) and use it to rank Boltz-
   passing candidates. Different model, different data, different
   purpose from the steering generator above.

5. **Pool diversity for the next run** — see the v11 pool-uniformity
   observation. If the assay benefits from sequence-binding diversity
   across designs at this single MPNN scaffold, expand the pool or
   reduce max_mutations.
