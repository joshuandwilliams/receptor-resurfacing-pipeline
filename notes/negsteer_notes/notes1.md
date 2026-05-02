# Negative steering — project notes

Living document. Update as experiments run. Last updated: 2026-04-10.

## What this project does

A multi-cycle pipeline for **negative steering of Boltz2 protein–effector
predictions**. Goal: when Boltz2 predicts the wrong binding site for a
known complex, identify the wrong-interface receptor residues, mutate
them disruptively, re-predict, and iterate. Each cycle treats the
previous (still-wrong) prediction as a new starting point and tries to
push the effector somewhere new. Targeting sub-3 Å (definitely
sub-5 Å) `receptor_aligned_effector_rmsd` vs ground truth.

Built on top of `boltz2_negative_steering.py` (single-cycle core,
~1900 lines, stable) with `boltz2_iterate_steering.py` as the
multi-cycle orchestrator. SLURM-dispatched, container-isolated GPU
prediction stages, host-side CPU collect/sbatch stages.

## Current pipeline state (as of v5)

**Filter chain:** `generated → intact → novel → selected`

1. **Intact:** `independent_receptor_rmsd ≤ 3.0` AND
   `independent_effector_rmsd ≤ 3.0`. Drops candidates where
   either chain has folded badly enough that downstream metrics
   are unreliable.
2. **Novel:** `min_dist_to_upstream ≥ --novelty-cutoff` (default
   10 Å). Distance is computed against *every* prior prediction
   in the pathway including the wild-type. Forces each cycle to
   find a binding mode that hasn't been seen before in the chain.
3. **Selected:** maximin by `min_dist_to_upstream`, with
   `receptor_aligned_effector_rmsd_vs_truth` as a tie-breaker
   only. Up to `--max-passing` survivors per parent per cycle.

Mutation accumulation: cycle N+1's candidate pool excludes positions
already mutated in any prior cycle in the pathway. The protected set
(true binding-site overlap) carries forward unchanged across cycles.
Effector template is always pinned to ground-truth.

**Outputs per run:**
- `params.yml` — full snapshot at submission time (see below)
- `submission.log` — captured banner + sbatch JIDs + footer
- `cycle_statistics.csv` — appended incrementally as each pathway
  finishes its collect stage; one row per (cycle, pathway) with
  full distance/ra_eff lists for bell-curve analysis
- `all_results_multicycle.csv` — every prediction across the tree,
  full per-design metrics, leaf pathway labels
- `pathways.json` — tree structure with parent linkage

## Key files

| file | role |
|---|---|
| `boltz2_negative_steering.py` | single-cycle core: plan → predict-one array → collect. Unchanged across all multi-cycle work. |
| `boltz2_iterate_steering.py` | multi-cycle orchestrator. Subcommands: `iterate-plan`, `compute-distances`, `iterate-collect`, `kickoff-distances`, `kickoff`, `aggregate`. |
| `submit_boltz2_negative_steering.sh` | top-level submitter. Single-cycle if `--n-cycles 1`, else writes `params.yml`/`submission.log` and chains a kickoff job after cycle 0. |
| `submit_boltz2_iterate_steering.sh` | per-cycle submitter for cycles 1+. Self-resubmitted by `iterate-collect` for surviving pathways. |
| `boltz2_verify_binding.py`, `boltz2_filter.py` | validation/filtering helpers used inside the container. |

## Architecture invariants worth knowing

- **Phase split:** every pathway's collect job runs as `singularity exec
  ... compute-distances && python3 ... iterate-collect`. Phase 1 needs
  numpy/Bio.Align (in container); phase 2 needs sbatch (host only).
  The container lacks SLURM client tools, so sbatch must run on the
  host side.
- **Pathway directories are named after parent labels.** Cycle-N
  pathway directory is `cycle_N/pathway_<parent_label>/`. The N designs
  inside it have *leaf labels* like `<parent>.cNd00`, `<parent>.cNd01`,
  etc. Leaf labels are design entries, NOT directories. This bit me at
  cycle 2 in v2 — see "Bug history" below.
- **Cycle-0 lives at `cycle_0/`** — no pathway prefix, because the
  wild-type is the only "parent."
- **`distances.json` uses `None`, not `NaN`** for missing values, so
  the host-side phase 2 doesn't need numpy to parse it.
- **Aggregator has self-healing** for old runs that lack the embedded
  `result` field in passing.json — it falls back to reading
  `result.json` directly.
- **`exec > >(tee -a "$log") 2>&1`** plus a `sleep 0.1` in
  `cleanup_on_failure` is what makes `submission.log` capture the full
  banner without truncating the tail when the script exits.

## Bug history

1. **Aggregate row labels were parent labels not leaf labels.** Fixed
   by computing `make_pathway_label(parent, cycle, design_idx)` per
   row with a separate `parent_pathway` column.
2. **`cmd_aggregate` glob `cycle_*` matched `cycle_statistics.csv`.**
   Fixed with `is_dir()` check.
3. **Cycle 2+ chain-walker bug (broke v2, also broke v3 before fix).**
   `read_upstream_chain` was calling `pathway_workdir(leaf_label)`,
   which built a path like `cycle_1/pathway_c0d00.c1d01/` — a directory
   that doesn't exist. Cycle-1 pathway directories are named
   `cycle_1/pathway_c0d00/` (after the *parent* label), and `c1d01` is
   a design entry inside it, not a directory. Fixed by deriving the
   parent label inline and using `cycle_N/pathway_<parent_label>/` as
   the workdir for cycle N>0 entries. The bug only manifested at
   cycle 2 because cycle 1 happens to work by accident (single-segment
   parent label = single-segment directory name).
4. **kickoff calling sbatch from inside the container.** The container
   lacks SLURM client tools. Fixed by splitting kickoff into
   `kickoff-distances` (in container) and `kickoff` (on host).
5. **The quality filter itself.** See "Lessons learned" below — added
   in the gap between v1 and v2, removed entirely after v5 confirmed
   it was based on a wrong mental model.

## Experiment log (7QPX)

Target: 7QPX, receptor chain A (73 residues), effector chain C
(81 residues). Wild-type Boltz2 prediction is ~31.47 Å from truth —
i.e., Boltz2 confidently predicts the wrong binding site for this
complex. Cycle 0 with `--mode mild` reliably gets to ~14.46 Å on a
*different* wrong site. The objective is to push past 14.46 Å toward
sub-5 Å.

| run | cycles | designs/cyc | filter on? | best ra_eff (intact) | tree size | notes |
|---|---|---|---|---|---|---|
| v1 | 1 | 6 | n/a | 14.46 Å (cycle 0) | 6 | baseline single-cycle |
| v2 | 3 | 6 | yes | n/a | crashed | killed by chain-walker bug at cycle 2 |
| v3 | 3 | 6 | yes | 14.46 Å (cycle 0 still best) | ~50 | bug fixed; best tJac=0.29 in cycle 2 |
| v4 | 5 (planned) | 6 | yes | 14.46 Å | 49 | tree died at cycle 3, only 1 survivor; quality filter killed legitimately novel modes |
| v5 | 5 | 6 | NO | **14.46 Å** still | 331 | ~7× wider tree, still no improvement; best tJac=**0.31** at cycle 3 |

**The 14.46 Å plateau is structural to this target.** Two different
selection regimes (with filter and without) explored very different
amounts of pose space and both bottomed out at the same number.

### Best results in v5 (intact only)

| pathway | cycle | ra_eff | tJac | wJac |
|---|---|---|---|---|
| `cycle_0/design_01` | 0 | 14.46 | 0.18 | 0.36 |
| `c0d01.c1d03` | 1 | 14.60 | 0.22 | 0.82 |
| `c0d01.c1d05.c2d00.c3d04` | 3 | 14.66 | 0.20 | 0.03 |
| `c0d02.c1d05.c2d00.c3d04.c4d04` | 4 | 14.95 | 0.20 | 0.16 |
| `c0d00.c1d05.c2d03.c3d03` | 3 | 27.26 | **0.31** | 0.46 |

Note that the best ra_eff (14.46–14.95) and the best true_jaccard
(0.31) come from completely different pathways. **The two metrics
disagree on this target** and ra_eff is misleading — it's possible to
get a "good" ra_eff by binding a *third* wrong site that happens to be
the same distance from truth as the original wrong site. true_jaccard
is the metric showing genuine progress: it climbed from 0.18 (cycle 0)
→ 0.22 (cycle 1) → 0.29 (cycle 2) → 0.31 (cycle 3) across cycles.
That's the only signal that iteration is doing anything useful.

## Lessons learned

### The quality filter was based on a wrong mental model

Between v1 and v2, the orchestrator gained `--min-improvement` (drop
candidates worse than wild-type by more than this) and
`--quality-weight` (weight on the improvement term in maximin scoring).
The intuition was hill-climbing: "designs that move closer to truth
are better candidates for the next cycle than designs that don't."

This is wrong for the kind of pipeline we're running. Each Boltz2
prediction is a fresh sample from the model's pose distribution
given the perturbed sequence — not a gradient step from the parent.
Cycle N's pose doesn't predict cycle N+1's pose. The only thing
that actually carries forward between cycles is the *mutation set*.
A design at ra_eff=45 Å (very wrong) is just as likely to produce a
sub-5 Å child as a design at ra_eff=20 Å, because the child is
sampled fresh from a different perturbed receptor.

What the quality filter actually did: it killed novel binding modes
that happened to be in worse-than-wild-type parts of pose space, even
though those were exactly the directions the iteration was supposed
to explore. v4 cycle 3 had `c0d00.c1d01.c2d03.c3d02` (intact, novel,
ra_eff=45.85, min_dist=31.65) — the most novel design of the cycle —
killed by the filter and the pathway terminated. Without that filter
(v5), the tree grew to 331 designs across 56 pathways. Removed
entirely after v5; the maximin selector now uses pure novelty with
ra_eff as a tie-breaker only.

### Novelty filter is doing its job

In v4 design_01 of cycle 3 was an intact design with the same wrong
site Boltz2 keeps snapping back to (min_dist 7.35 Å from one upstream
pose). Killed by the novelty filter. **This is the correct behaviour**
— the whole point of negative steering is to *not* keep iterating
around the same wrong site after we've already disrupted the residues
that were producing it. If 4 cycles of mutation accumulation can't
push Boltz2 off a pose, more cycles of variants of the same pose
won't help either.

### ra_eff is a poor optimisation target on this kind of run

When you're trying to find a binding site, "distance from truth" is
ambiguous. A receptor surface 14 Å from the true site is at the same
ra_eff as the true site rotated 180°. Two designs can have identical
ra_eff and be in completely different parts of pose space. On 7QPX,
v5 found multiple designs at ~14–15 Å scattered across different
wrong surfaces — none of them on the right one.

**`true_jaccard` (overlap of design's contact residues with the true
binding-site residues) is the metric with signal on this target.**
It climbed monotonically across cycles (0.18 → 0.22 → 0.29 → 0.31)
even while ra_eff stayed flat. By cycle 3 we have a design where
~30% of the receptor's contact residues are on the *correct* binding
surface, even though the effector pose is geometrically wrong. This
suggests the iteration is working at the contact level — finding
receptor surfaces that look like they want to bind there — but the
pose itself doesn't follow because Boltz2's prior over poses is too
strong.

### Targets that are too well-predicted by Boltz2 are bad targets

7QPX wild-type prediction is at 31 Å, and cycle 0 with mild
mutations gets to 14 Å. That's a big single-shot improvement, but
it leaves very little headroom for cycles to do work. **The targets
where iteration would actually help are ones where cycle 0 only
gets to ~25 Å.** Look for those next.

## Open work items

1. **`--target-true-jaccard` benchmark mode.** Now genuinely
   actionable based on v5 results. Add an alternative scoring
   objective alongside (not replacing) the current novelty-based
   maximin: when ground truth is known, score survivors by
   `true_jaccard - parent_true_jaccard` so the selector rewards
   designs moving the contact pattern toward the truth. Suggested
   flag: `--scoring novelty|jaccard` (default `novelty`). Only
   useful for benchmark runs against solved structures.
2. **Try `--mode alanine` on 7QPX.** Most disruptive mutation
   strategy. If anything is going to dislodge Boltz2 from the
   ~14 Å wrong site on this target, it's this. One run, ~3-4
   hours, definitive answer for whether negative steering can
   crack 7QPX with current code.
3. **Find or pick a better benchmark target.** Criteria: Boltz2
   wild-type prediction should be ≥25 Å from truth, cycle 0 with
   `mild` mode should plateau in the 18–25 Å range (not 14 — that's
   too easy and leaves no headroom). Iteration only buys you
   something on targets where one-shot mutation doesn't already
   half-solve the problem.
4. **Consider adding a per-cycle novelty cutoff schedule.** Right
   now `--novelty-cutoff 10.0` is fixed across all cycles. By cycle
   3+, the upstream pose set has grown to 4-5 reference poses, and
   requiring 10 Å from each is an increasingly large region to
   exclude. Could be worth allowing the cutoff to relax with cycle
   depth, e.g. `10 / (1 + 0.2*cycle)`. Speculative — only useful if
   future runs hit "exhausted from over-strict novelty" rather than
   the 14 Å plateau we keep seeing.

## Things to remember when working on this code

- **Minimal targeted changes only.** Catches over-refactoring early.
- **Run submit scripts via `./script.sh` not `sh script.sh`** — dash
  doesn't have bash arrays and the scripts use them.
- **All numpy/gemmi/BioPython only inside container.** Phase 2 must
  be pure stdlib.
- **Never call sbatch from inside singularity** — no SLURM client
  tools in the image. Use the phase split.
- **Cycle-0 plan.json `designs` list must be dense** — all indices
  0..N-1 present, even for designs that fail.
- **`make_pathway_label(parent, cycle, design_idx)` is the only
  way to construct leaf labels.** Don't string-concatenate.
- **If `read_upstream_chain` is showing weird burnt-position sets,
  check that it's using `cycle_N/pathway_<parent_label>/` and not
  `cycle_N/pathway_<leaf_label>/`.** That bug broke v2 and would
  break again the same way if someone "simplifies" it.
