# 15 — Discovery run path coverage and fixture candidates

Companion to `14_phase_2_revision_per_module_tests.md` §2.4 / §2.5. Inventories
which designs and sequences from the discovery run (`params_full_test.yml`,
results at `tests/full_test_run/results/` on HPC) cover which paths through
each stage, and proposes minimum-coverage candidates for the per-module test
fixtures. Verification commands at the end of each section run on the HPC
login node from inside `tests/full_test_run/results/`.

The taxonomies in this report are derived from the codebase's own vocabulary,
not pre-imposed. Where notes documents disagreed I cite the producer code as
the tie-breaker and flag the disagreement under "Open questions".

The candidate file sets in this report are based on five helper artefacts the
human pulled to Mac at `tests/full_test_run/reference_data_helpers/` —
`rfdiffusion_metrics.json`, `rosetta_filter_metrics.json`, `top_metadata.csv`,
`cross_sequence_summary_with_interface_metrics.csv`,
`survivors_with_orthogonal_metrics.csv`. Per-design and per-seed schema
details inside the per-sequence workdirs were inferred from the producer
scripts and from the archived supervisor-demo layout at
`tests/full_test_run/example_output_files/negative_steering/runs/<seq>/`.

## Summary

Class distribution observed for the 122 cross-sequence rows the discovery run
produced (120 steered + 2 controls), in the producer code's vocabulary:

| Aggregated class (representative-row level) | Count | Tier dist | Sub-case |
|---|---:|---|---|
| `no_reversion` + `representative_design == "initial"` | 13 | 13 A | `cold_start_all_clean` (notes12 / patch b) |
| `no_reversion` + `representative_design != "initial"` (tiered) | 5 | 5 A/B/C | `steered_clean`-tiered: at least one design had `steered_ra_eff < 5.0` AND intact, makes it into passing_summary |
| `no_reversion` + `representative_design != "initial"` (tier-none, all_zeros) | 63 | 63 none | `steered_clean`-not-tiered: every design `no_reversion` but no design's steered ra_eff < 5.0 → empty passing_summary, tier-none fallback fires |
| `no_reversion` representative + sequence has ≥1 `pose_collapses` design | 6 | 6 none | mixed: most designs `no_reversion`-but-bad-pose; at least one design had reversion run that aggregated to `pose_collapses`. passing_summary still empty |
| `pose_holds` | 9 | 3 A, 5 B, 1 C | reversion ran, aggregator passed |
| `pose_collapses` representative | 26 | 26 none | reversion ran, aggregator failed for every design that had reversion |
| `new_contamination` | 0 | — | **NOT OBSERVED** at the aggregated level |
| `singleton` (per-row, not per-sequence) | 122 (one per sequence) | n/a | every sequence has exactly one singleton row — see Axis 3 below |

Candidate counts proposed per stage: RFDiffusion 1 (no fail observed),
Rosetta 4 (1 marginal pass + 1 strong pass + 2 fails), ProteinMPNN 4 (2
designs feeding 4 sequences across negsteer outcomes), Negative steering 7
(2 cold-start, 1 steered-clean, 1 pose_holds-3/3, 1 pose_holds-2/3, 1
pose_holds-1/3, 1 pose_collapses), Orthogonal metrics 2 (one tier-A
representative + one near-passing). Several taxonomy classes are uncovered
or thinly covered — see "Coverage gaps and risks".

## Per-stage analysis

### RFDiffusion

#### Paths

`bin/rfdiffusion_filter.py` writes one `passes_filter` boolean per design
into `rfdiffusion/rfdiffusion_metrics.json["designs"][i].passes_filter`.
The two paths are:

- **`passes_filter == true`** — design has enough Cα contacts in the
  designed region (and meets `min_hotspot_frac` if hotspots were supplied).
  Only these are split into two-chain PDBs under `rfdiffusion/split/` and
  written to `rfdiffusion/passing_designs.txt`.
- **`passes_filter == false`** — design dropped before the split step;
  appears in the metrics JSON but not in `split/` or `passing_designs.txt`.

For the per-module RFDiffusion test, the path classification only matters
downstream: the test exercises RFDIFFUSION → RFDIFFUSION_FILTER →
RFDIFFUSION_PLOTS on a fixed input PDB; what we need from the discovery
run is the input PDB and contigs (already on disk under
`tests/full_test_run/af3_pikp1_native_avrpikf_complex.pdb`) plus a sense of
the metric ranges the filter sees.

#### Discriminating signals

| Path | File | Field |
|---|---|---|
| pass | `rfdiffusion/rfdiffusion_metrics.json` | `designs[i].passes_filter == true` |
| pass | same | `designs[i].n_contact_pairs` (≥1), `frac_contacts_in_design`, `design_region_coverage` |
| pass | `rfdiffusion/passing_designs.txt` | design stem present |
| pass | `rfdiffusion/split/<stem>.pdb` | file exists |
| fail | `rfdiffusion_metrics.json` | `passes_filter == false`, typically `n_contact_pairs == 0` |
| fail | `passing_designs.txt` | absent |
| fail | `rfdiffusion/split/` | no file |

#### Proposed candidates

The discovery run's filter config (`min_hotspot_frac=0.0` and `hotspot=""`
in `params_full_test.yml`) means the contact filter degenerates to "any
contact at all". From `rfdiffusion_metrics.json` (helper): **64/64 designs
pass** — the fail path was not exercised. So the only candidates are
pass-path designs:

- **design_0** — pass, `n_contact_pairs=24`, `frac_contacts_in_design=0.46`,
  `design_region_coverage=0.35`. Median-ish metrics; representative.
- (no fail candidate available; see Coverage gaps)

For the per-module fixture, the test workflow consumes only
`params.pdb_file` (the input complex) and the contigs string — both
already in `tests/rfdiffusion/data/` per the decoupling work. The
discovery-run output is **not needed as a fixture input** for this test;
it just provides reference outputs. The reference set should mirror what
this test produces from its own input.

#### Verification commands

```bash
# Confirm all 64 designs pass; expect "n_total=64 n_pass=64 n_fail=0"
python3 -c '
import json
d=json.load(open("rfdiffusion/rfdiffusion_metrics.json"))
ds=d["designs"]
print("n_total=%d n_pass=%d n_fail=%d" % (len(ds),
  sum(1 for x in ds if x["passes_filter"]),
  sum(1 for x in ds if not x["passes_filter"])))'

# Confirm passing_designs.txt count matches
wc -l rfdiffusion/passing_designs.txt
# Should print "64 ..."

# Spot-check that design_0 is split out
ls rfdiffusion/split/design_0.pdb
# Should print the path; absent file means RFDiffusion was not split at all

# Stage-output sizes for the eventual fixture-tarball footprint
du -sh rfdiffusion/ rfdiffusion/split/ rfdiffusion/traj/ 2>/dev/null
```

### Rosetta filtering

#### Paths

`bin/rosetta_filter_collect.py` writes one `passes_filter` per design into
`rosetta_filtering/rosetta_filter_metrics.json["designs"]`. Threshold is
`sc_value >= sc_threshold` (default 0.5; `sc_threshold` is recorded in the
JSON). The two paths:

- **pass** — `sc_value >= sc_threshold`; design's split PDB is copied into
  `rosetta_filtering/passing/`.
- **fail** — `sc_value < sc_threshold`; design appears in the metrics JSON
  but not in `passing/`.

The per-module test (`test_rosetta_filtering.nf`) takes split PDBs from
`data/rfdiffusion_split/design_*.pdb`. To exercise both filter outcomes,
the fixture needs PDBs for both passing and failing designs.

#### Discriminating signals

| Path | File | Field |
|---|---|---|
| pass | `rosetta_filtering/rosetta_filter_metrics.json` | `designs[i].passes_filter == true`, `sc_value >= 0.5` |
| pass | `rosetta_filtering/passing/<stem>.pdb` | file exists |
| pass | `rosetta_filtering/rosetta_filter_summary.csv` | row present, `passes_filter == True` |
| fail | metrics JSON | `passes_filter == false`, `sc_value < 0.5` |
| fail | `rosetta_filtering/passing/` | no file |

#### Proposed candidates

From the helper `rosetta_filter_metrics.json` (61 pass, 3 fail):

- **design_0** — strong pass, `sc=0.569`. (Note: this design also has
  Rosetta metrics in the helper's design[0] entry, sc 0.569; cross-check
  with `sc=0.504` for design_1 below — they're different designs.)
- **design_42** — top pass, `sc=0.633`. Best-case for the pass path.
- **design_1** — marginal pass, `sc=0.504`. Boundary case for the gate.
- **design_14** — fail, `sc=0.488` (just below threshold). Boundary fail.
- **design_59** — fail, `sc=0.405`. Clear fail.

Recommendation: curate **design_0, design_1, design_14, design_59** into
`tests/rosetta_filtering/data/rfdiffusion_split/`. Four PDBs covers both
gate outcomes plus boundary cases either side of `sc_threshold=0.5`.

Files needed per design:
- `rfdiffusion/split/<stem>.pdb` — the two-chain (A=receptor, B=effector)
  split-PDB produced by `rfdiffusion_filter.py:write_split_pdb`. Copy
  these into `tests/rosetta_filtering/data/rfdiffusion_split/`.

No other inputs needed; the test reads only the glob.

#### Verification commands

```bash
# Confirm fail count and sc values; expect 3 fails: design_14, design_25, design_59
python3 -c '
import json
d=json.load(open("rosetta_filtering/rosetta_filter_metrics.json"))
fails=[x for x in d["designs"] if not x["passes_filter"]]
print("threshold:", d["sc_threshold"])
print("n_fail:", len(fails))
for x in fails: print(" ", x["design_stem"], "sc=", x["sc_value"])'

# Confirm marginal pass design_1 sc value; expect ~0.504
python3 -c '
import json
d=json.load(open("rosetta_filtering/rosetta_filter_metrics.json"))
for x in d["designs"]:
    if x["design_stem"] in ("design_0","design_1","design_42"):
        print(x["design_stem"], "sc=", x["sc_value"])'

# Confirm split-PDB existence for the four candidates
for s in design_0 design_1 design_14 design_59; do
  ls -la rfdiffusion/split/${s}.pdb 2>&1 | head -1
done
# Each should print a real file path; "No such file" means the split did
# not happen for that design (would block the candidate).

# Stage-output sizes
du -sh rosetta_filtering/ rosetta_filtering/passing/
```

### ProteinMPNN

#### Paths

ProteinMPNN itself doesn't have a binary pass/fail per design; it samples
`num_seqs` sequences per input PDB and the downstream MPNN_SELECT_TOP
process picks the top-N by composite score across all sequences. The
discovery-run partition surfaces at the **negsteer stage**:

- A given (design, seq) pair downstream lands in one of the negsteer
  classes (cold_start_all_clean / steered-clean / pose_holds /
  pose_collapses / new_contamination / no_reversion-tier-none).
- For ProteinMPNN's per-module test, what we want is fixture inputs
  (RFDiffusion split-PDBs) such that running MPNN + downstream-equivalent
  scoring produces both the "produces a viable cold-start input" and
  "produces a sequence that fails downstream" partitions.

But the per-module MPNN test as written (`test_proteinmpnn.nf`) does not
itself run negsteer; it stops at MPNN_SELECT_TOP. So the partition the
test exercises is "appears in `mpnn_top_n` selection" vs "filtered out by
QC / clustering". The discriminating signals at the test's own output
boundary are:

- **selected for top-N** — design appears in
  `sequences/top_metadata.csv`; FASTA exists at
  `sequences/top_fastas/design_<N>_seq_<S>.fasta`
- **dropped by QC** — design's MPNN sequence has `qc_pass == False` in
  `sequences/qc_metadata.csv` (or is absent from
  `sequences/qc_fastas/`)
- **clustered out** — sequence assigned to a non-canonical cluster
  representative; visible in MPNN_CLUSTER's output under the per-design
  `mpnn/design_<N>/` subtree

`bin/mpnn_select_top.py` writes `top_metadata.csv` to the dir passed
via `--output-dir`, and `MPNN_SELECT_TOP` in `modules/proteinmpnn.nf`
publishes that dir to `${params.outdir}/sequences/`. Per-design MPNN
output lives at `${params.outdir}/mpnn/<pdb_baseName>/` (line 60 of
proteinmpnn.nf). The discovery run's `sequences/top_metadata.csv`
(helper) has 120 rows across 61 unique designs (61 × 2 sequences each,
with one design having only one). Every passing-Rosetta design appears
to feed at least one MPNN sequence.

#### Discriminating signals

| Path | File | Field |
|---|---|---|
| selected for top-N | `sequences/top_fastas/design_<N>_seq_<S>.fasta` | exists |
| selected | `sequences/top_metadata.csv` | row with `design=N seq=S` |
| selected | `mpnn/design_<N>/mpnn_results.json` | sequence present |
| QC-dropped | `sequences/qc_fastas/design_<N>_seq_<S>.fasta` | absent |
| QC-dropped | look at `qc_metadata.csv` `qc_pass` column | `False` |

#### Proposed candidates

Pick 2 designs from `top_metadata.csv` whose MPNN sequences span the
downstream partition (one becomes cold_start_all_clean, one fails to
reach tier A/B/C). Cross-referencing with the cross_sequence_summary
helper:

- **design_28** — both sequences make it through; `design_28_seq_1`
  becomes the top-ranked tier-A `cold_start_all_clean` representative
  (cross_rank 1, ra_eff 2.72). Strong "feeds clean cold-start" example.
- **design_0** — `design_0_seq_0` becomes tier=none/no_reversion (steered
  but no design passed); `design_0_seq_1` becomes tier-A/no_reversion
  cold-start. Within-design split.

Recommendation: curate **design_0, design_28** into
`tests/proteinmpnn/data/rosetta_passing/` as the fixture input PDBs (these
go through MPNN), plus the existing `data/input_complex.pdb` (already on
disk from the per-module decoupling work).

Files needed per design:
- `rosetta_filtering/passing/<stem>.pdb` — the post-Rosetta two-chain
  PDB. Copy into `tests/proteinmpnn/data/rosetta_passing/`.

#### Verification commands

```bash
# Confirm both candidate designs are in MPNN top selection; expect both stems present
python3 -c '
import csv
with open("sequences/top_metadata.csv") as f:
    rows=list(csv.DictReader(f))
for d in ("0","28"):
    seqs=[r["seq"] for r in rows if r["design"]==d]
    print("design_%s: %d seqs (%s)" % (d, len(seqs), ",".join(seqs)))'
# Should print "design_0: 2 seqs (0,1)" and "design_28: 2 seqs (0,1)"

# Confirm the rosetta_passing PDBs exist for both
ls -la rosetta_filtering/passing/design_0.pdb rosetta_filtering/passing/design_28.pdb 2>&1

# Confirm per-design MPNN output dirs exist
ls -d mpnn/design_0 mpnn/design_28 2>&1

# Stage-output sizes
du -sh mpnn/ sequences/top_fastas/ sequences/qc_fastas/
```

(The first command in the previous version of this report pointed at
`mpnn/top_metadata.csv`, which doesn't exist — `MPNN_SELECT_TOP`
publishes to `${params.outdir}/sequences/`, not `mpnn/`. The corrected
path is verified against `modules/proteinmpnn.nf:315`.)

### Negative steering

#### Pathway taxonomy

Derived from `bin/boltz2_negative_steering.py::cmd_plan` (cold-start
outcome), `bin/reversion.py::classify_reversion_verdict` (per-seed
post-reversion verdict), `bin/boltz2_iterate_steering.py::_per_seed_verdict_breakdown`
and `_classify_aggregated_verdict` (per-seed counts and aggregation rule),
and `bin/cross_sequence_summary.py::_tier_for_row` (cohort tier).
Pre-`14`-era assistant prompts that talked about "P1 cold-start fail / P2
clean steering / P3 reversion succeeds / P4 reversion fails" do not match
this vocabulary and should not be used.

The taxonomy has **three orthogonal axes** that the codebase tracks
separately:

##### Axis 1 — Cold-start outcome (per `boltz2_negative_steering.py:1937–2127`)

Per cold-start seed (default `num_seeds=3`, plan stage runs all of them):

- **clean** — `ra_eff <= rmsd_threshold` AND `ind_rec <= 5.0` AND
  `ind_eff <= 5.0` (see `boltz2_negative_steering.py:2018-2030`,
  `RECEPTOR_INTACT_CUTOFF` / `EFFECTOR_INTACT_CUTOFF`).
- **dirty** — any of the three checks fail.

Per-sequence aggregation:
- `cold_start_all_clean = (n_clean == n_total_cold and n_total_cold > 0)`
  (line 2036). When true, `cmd_plan` sets `plan["skip_steering"] = True`
  and writes per-seed cold-start rows via `_write_initial_multiseed_csv`.
  Steering is skipped entirely; the cold-start predictions ARE the
  candidate set.

There is **no separate "cold-start fail" path**. A sequence whose
cold-start fails just enters the steered-design pipeline normally; the
cold-start results are screening-only artefacts.

##### Axis 2 — Per-seed verdict for steered designs (per `_per_seed_verdict_breakdown`)

Each (negsteer-design, seed) pair emerges from the steered + reversion
chain with one of five per-seed states (`bin/boltz2_iterate_steering.py:3689–3722`):

- **`clean_steered`** — `steered_n_contacts_on_mutated_positions == 0`;
  no reversion needed; the steered structure is the final per-seed result.
- **`pose_holds`** — reversion ran; `classify_reversion_verdict` returned
  `pose_holds` (reverted prediction structurally OK, no gated-position
  contamination).
- **`pose_collapses`** — reversion ran; reverted prediction failed
  intact filter or `reverted_ra_eff >= 5.0`.
- **`new_contamination`** — reversion ran; reverted prediction has at
  least one mutated residue contacting the effector at a position inside
  the **gated set** (design region ∪ true interface). The gating was
  added in notes11 to mirror the pre-reversion build-contaminated
  policy; before notes11, *any* mutated contact would fire
  new_contamination (`bin/reversion.py:1062–1144`).
- **`no_data`** — seed prediction errored or has no usable metrics.

##### Axis 3 — Aggregated `aggregated_verdict` (per `_classify_aggregated_verdict`)

The per-(negsteer-design)-row aggregator in
`boltz2_iterate_steering.py:3514–3686` re-derives a verdict from the
median + position-majority columns. Possible values:

- **`no_reversion`** — `any_reverted` was False over the design's seeds
  (every seed was `clean_steered` OR no seed went through reversion at
  all). All `n_seeds_*` columns blank.
- **`pose_holds`** — aggregated structural filter passes
  (`reverted_receptor_intact_majority == 1` AND
  `reverted_ra_eff_vs_truth_median < 5.0` AND no gated-position
  contamination by majority) AND `n_pass_eff = n_seeds_pose_holds +
  n_seeds_clean_steered >= 1` AND failure modes don't dominate
  (`n_pass_eff + n_seeds_no_data >= max(n_pose_collapses, n_new_contamination)`).
- **`pose_collapses`** — aggregator's structural filter or any-seed gate
  failed; downgrade. Overrides on majority-failure too.
- **`new_contamination`** — aggregated reverted contact set has a
  position-majority position inside the gated set.
- **`singleton`** — emitted by the aggregator at
  `boltz2_iterate_steering.py:4596–4604` for any per-design row whose
  `sequence_group` is blank (`""` or `None`). Such a row is treated as
  its own group of size 1 with no aggregation, gets `n_seeds=1`, and is
  written into `aggregated_results.csv` with `aggregated_verdict =
  "singleton"`. **The cycle-0 cold-start baseline row written by
  `cmd_collect` has `sequence_group=""` by construction**, so every
  per-sequence `aggregated_results.csv` carries exactly one `singleton`
  row regardless of `num_seeds`. Verified locally against
  `tests/full_test_run/example_output_files/negative_steering/runs/design_0_seq_0/aggregated_results.csv`
  (and design_13_seq_2): each shows one row with `sg='', design='initial',
  verdict='singleton', n_seeds=1`. The g-deep verification on HPC
  found exactly 69 singleton rows across 69 Class 3 sequences (1 per
  sequence) and the design_27_seq_0 verdict-set check returned
  `{'no_reversion', 'pose_collapses', 'singleton'}`, both consistent
  with this rule. The previous version of this report incorrectly
  stated singletons were "not produced when num_seeds > 1"; that
  claim was wrong.

Per-seed columns (`n_seeds_pose_holds`, `n_seeds_pose_collapses`,
`n_seeds_new_contamination`, `n_seeds_no_data`,
`n_seeds_clean_steered`) are populated **only when reversion ran**
for that design (i.e. `aggregated_verdict ∈ {pose_holds,
pose_collapses, new_contamination}`). For `no_reversion` and
`singleton` rows the aggregator writes blank strings to all five
columns (`boltz2_iterate_steering.py:4705–4709`). Numeric coercion
of those blanks to integers yields 0; a row that reads as
`n_seeds_clean_steered=0, n_seeds_no_data=0, …` after coercion is
indistinguishable from one where steering ran with all-failing seeds.
Discriminating the two requires reading the `aggregated_verdict`
column directly, not the per-seed counts.

The "any-seed gate" (notes10 Bug-E fix) replaces an earlier
ceil(num_seeds/2) majority gate; with `num_seeds=3` and a per-seed
pose_holds rate of ~10–15%, majority-of-3 was architecturally
unreachable, so the aggregator was relaxed.

##### Axis 4 — Cohort tier (per `cross_sequence_summary.py::_tier_for_row`)

Computed downstream of the aggregator from the columns above; **not**
written into per-sequence aggregated_results.csv but rather into
`cross_sequence_summary.csv["cross_tier"]`:

- `n_pass = n_seeds_pose_holds + n_seeds_clean_steered`
- **A**: `aggregated_verdict == "no_reversion"` OR `n_pass == n_seeds`
- **B**: `1 < n_pass < n_seeds`  (e.g. 2/3)
- **C**: `n_pass == 1` (e.g. 1/3)
- **none**: otherwise (zero pass-equivalent seeds, or pose_collapses /
  new_contamination aggregated, or insufficient data).

##### Per-MPNN-sequence outcome classes actually observable

A representative is picked per MPNN sequence (best-non-empty-tier wins;
within tier, lowest `rank_by_composite_score`). Cross-tabbing
`(cross_tier, aggregated_verdict, representative_design)` over the
helper file gives the classes the discovery run produced:

| # | Class (producer vocabulary) | Predicate | Discovery count |
|---|---|---|---:|
| 1 | `cold_start_all_clean` | tier=A, verdict=`no_reversion`, rep_design=`initial` | 13 |
| 2 | `steered_clean` (tiered) | tier∈{A,B,C}, verdict=`no_reversion`, rep_design≠`initial` — at least one design hit `steered_ra_eff < 5.0` AND intact, made it into passing_summary | 5 |
| 3a | `steered_clean` (tier-none, all-`no_reversion`) | tier=`none`, verdict=`no_reversion`, rep_design≠`initial` AND aggregated_results.csv contains zero `pose_collapses` / `pose_holds` / `new_contamination` rows | 63 |
| 3b | `steered_clean` (tier-none, mixed with `pose_collapses`) | tier=`none`, verdict=`no_reversion` representative, but at least one design in aggregated_results.csv had `aggregated_verdict == "pose_collapses"` | 6 |
| 4 | `pose_holds` aggregated, all-3-seed | tier=A, verdict=`pose_holds`, n_pass=3/3 | 3 |
| 5 | `pose_holds` aggregated, majority | tier=B, verdict=`pose_holds`, n_pass=2/3 | 5 |
| 6 | `pose_holds` aggregated, minority | tier=C, verdict=`pose_holds`, n_pass=1/3 | 1 |
| 7 | `pose_collapses` aggregated representative | tier=`none`, verdict=`pose_collapses` | 26 |
| 8 | `new_contamination` aggregated | tier=`none`, verdict=`new_contamination` | **0 (NOT OBSERVED at the aggregated level; per-seed level also 0 — see Coverage gap 2)** |

The Class 3 split was identified by the g-deep verification command on
HPC. The shape of the underlying `aggregated_results.csv` files
(verified locally for design_0_seq_0 and design_13_seq_2) confirms
**both 3a and 3b sequences have an empty `passing_summary.csv`**
because none of their design-aggregates pass the `_agg_row_is_clean_steered`
gate (`boltz2_iterate_steering.py:4749` — requires `steered_ra_eff <
5.0` AND `steered_receptor_intact_majority == 1`). The `cross_tier`
of `none` therefore arises from the cross_summary's
**tier-none-fallback** code path (`cross_sequence_summary.py:341–436`)
which picks a representative directly from `aggregated_results.csv`
when `passing_summary.csv` is empty. This is structurally distinct
from Class 2's path, which goes through `extract_passing.extract_row`
on a populated passing_summary. The human's reading is correct:
Class 3 = "steering ran but the steered designs landed in
confidently-misplaced poses (steered ra_eff ≥ 5.0)" — a meaningfully
different design content from Class 2's "steering ran and landed
correctly".

Class 3b's distinguishing feature is that one design in the sequence
reached the steered structural filter, was found contaminated by
`build-contaminated`, went through reversion, and the aggregator
classified the reversion as `pose_collapses`. Other designs in the
same sequence stayed `no_reversion` (clean steered but bad pose).
Cross_summary's representative pick for these 6 sequences is one of
the `no_reversion` rows (because the tier-none fallback picks by the
same composite-score logic, and a `pose_collapses` row's reverted
metrics often score worse than a `no_reversion`'s steered metrics).

The two **negative controls** (`input_control_polyA`,
`input_control_scrambled`) are generated by the test workflow itself
from `params.input_pdb` via `DERIVE_INPUT_INDICES` + `NEGSTEER_CONTROLS`
and are **not curated as fixtures** — they will appear automatically when
`test_negative_steering.nf` runs.

#### Discriminating signals

| Class | File | Field |
|---|---|---|
| `cold_start_all_clean` | `runs/<seq>/cycle_0/plan.json` | `cold_start_all_clean == true` AND `skip_steering == true` |
| `cold_start_all_clean` | same | `cold_start_n_clean == num_seeds` |
| `cold_start_all_clean` | `cross_sequence_summary.csv` | `representative_design == "initial"` AND `representative_aggregated_verdict == "no_reversion"` AND `cross_tier == "A"` |
| `steered_clean`, Class 2 / Class 3a (strict) | `runs/<seq>/aggregated_results.csv` | every non-singleton row has `aggregated_verdict == "no_reversion"`; no row has `reversion_verdict` populated |
| `steered_clean`, Class 3b (mixed) | `runs/<seq>/aggregated_results.csv` | most non-singleton rows are `no_reversion`, but at least one row has `aggregated_verdict == "pose_collapses"`; cross_summary still picks a `no_reversion` row as the representative |
| `steered_clean` (tiered) | `cross_sequence_summary.csv` | `representative_design != "initial"` AND `representative_aggregated_verdict == "no_reversion"` AND `cross_tier ∈ {A,B,C}` |
| `pose_holds` | `runs/<seq>/aggregated_results.csv` | at least one row with `aggregated_verdict == "pose_holds"` |
| `pose_holds` | same row | `n_seeds_pose_holds`, `n_seeds_pose_collapses`, `n_seeds_new_contamination`, `n_seeds_no_data`, `n_seeds_clean_steered` |
| `pose_holds` | `runs/<seq>/cycle_0/reversion_results.json` | per-(steered-label) entries with `verdict == "pose_holds"` |
| `pose_collapses` aggregated | `cross_sequence_summary.csv` | `representative_aggregated_verdict == "pose_collapses"` |
| `pose_collapses` aggregated | `runs/<seq>/cycle_0/reversion_results.json` | majority of entries `verdict == "pose_collapses"` |
| reversion was attempted at all | `runs/<seq>/cycle_0/contaminated.json` | `n_contaminated > 0` |
| reversion was attempted, what was reverted | `runs/<seq>/cycle_0/reversion_plan.json` | per-(label) `positions_to_revert` |
| reversion per-seed outcomes | `runs/<seq>/cycle_0/reversion_results_per_seed.json` | `{label: {seed_idx: result}}` |
| `new_contamination` | `runs/<seq>/cycle_0/reversion_results.json` | per-entry `verdict == "new_contamination"`, `reason` mentions "gated position(s)" |
| Class 3a vs Class 3b | `runs/<seq>/aggregated_results.csv` | 3b: at least one row has `aggregated_verdict == "pose_collapses"`; 3a: every non-singleton row is `no_reversion` |
| Class 3 vs Class 2 (passing_summary structural shape) | `runs/<seq>/passing_summary.csv` | Class 3: header-only CSV (zero data rows). Class 2: ≥1 data row. |
| singleton row presence | `runs/<seq>/aggregated_results.csv` | exactly one row with `sequence_group == ""` AND `design == "initial"` AND `aggregated_verdict == "singleton"` |

The `aggregated_verdict_reason` column on `aggregated_results.csv`
records the exact rule that fired, which is the most direct human-readable
signal for which downgrade path applied at the aggregator.

#### Proposed candidates

Aim for one fixture sequence per observable class (1–7 above), plus
boundary cases at the tier-A / tier-B and tier-B / tier-C transitions.
Class 8 (`new_contamination`) is uncovered — see Coverage gaps. Where the
class is large, pick the highest-cross-rank representative so the
fixture's outputs include the most-informative metric values.

| # | Candidate | Class claim | n_pass / n_seeds | Notes |
|---|---|---|---|---|
| 1 | **design_28_seq_1** | Class 1: `cold_start_all_clean` | n/a (no reversion) | Cross-rank 1, top representative. HPC verification (a) confirmed `cold_start_all_clean=True, skip_steering=True, cold_start_n_clean=3 of 3`. Workdir size 21 MB (smallest). |
| 2 | **design_5_seq_0** | Class 1 (second example) | n/a | Cross-rank 3. Diversity check: different parent design. *Optional* — may drop if footprint is tight. |
| 3 | **design_62_seq_0** | Class 2: `steered_clean` tiered (rep_design ≠ initial) | n/a (no reversion) | Tier A. Representative PDB at `cycle_0/steered/design_02_s0/prediction.pdb` (verified by HPC command (c)). Closest-to-passing on orthogonal metrics (af3_ra=14.53). Doubles as the Class-2 fixture and as one of the orthogonal-metrics survivors below. Workdir size 104 MB. |
| 4 | **design_0_seq_0** | Class 3a: `steered_clean` (tier-none, all-`no_reversion`) — empty passing_summary, tier-none fallback fires | n/a | g-deep verified: 4 design rows all `no_reversion`, 1 singleton, all per-seed columns blank. Parent `design_0` already curated for the Rosetta + MPNN fixtures, so no extra parent-design fixture needed. Exercises the empty-passing_summary branch of `cross_sequence_summary.py:341–436`. *Recommended add* — see decision note below. |
| 5 | **design_44_seq_1** | Class 4: `pose_holds` 3/3 (tier A) | 3/3 | HPC verification resolved Open Q3: representative is `design_00` with `n_pose_holds=1, n_clean_steered=2`, so `n_pass=3=n_seeds` → tier A. The helper file showed only one of those two contributing columns, which is what made the row look like `pose_holds 1/3` initially. Other designs in this sequence are `pose_holds 2/1/0` over their seeds. Mixed within-sequence example. |
| 6 | **design_42_seq_0** | Class 4: `pose_holds` 3/3 (tier A) | 3/3 | HPC verification (b) confirmed `verdict=pose_holds, n_pose_holds=3, n_pose_collapses=0, n_no_data=0, n_clean_steered=0, design=design_03`. Cleanest 3/3 example (no clean_steered admixture). Workdir size 147 MB. |
| 7 | **design_3_seq_1** | Class 5: `pose_holds` 2/3 (tier B) | 2/3 | HPC verification (c) confirmed `design=design_00, n_pose_holds=2, n_pose_collapses=1`. |
| 8 | **design_55_seq_1** | Class 6: `pose_holds` 1/3 (tier C) | 1/3 | Only tier-C row in the run. HPC verification confirmed the pose_holds row: `design=design_02, n_pose_holds=1, n_pose_collapses=0, n_no_data=2, n_clean_steered=0`. The sequence's `aggregated_results.csv` carries all four verdict types — `{'no_reversion', 'pose_collapses', 'pose_holds', 'singleton'}` — making it the most heterogeneous fixture in the set and the highest-value tier-C example. |
| 9 | **design_27_seq_0** | Class 7: `pose_collapses` aggregated, tier none | 0/3 | HPC verification (d) confirmed verdict-set `{'no_reversion', 'pose_collapses', 'singleton'}`; (e) confirmed `n_checked=4, n_contaminated=4`. Workdir size 142 MB. |
| 10 | **design_19_seq_0** | Class 7 boundary: `pose_collapses` aggregated despite 1 individual pose_holds seed | 1/3 in pose_holds, 1/3 collapse, 1/3 no_data | Bug-E-style downgrade case. *Optional*. |
| 11 | **(none available)** | Class 8: `new_contamination` aggregated | n/a | Not produced. See Coverage gap 2 — standalone task for tomorrow. |

**Class 3 fixture decision.** Recommend adding **design_0_seq_0** (Class 3a)
as candidate #4. Two reasons:

1. The per-module test exercises a **structurally different** code path
   for Class 3 vs Class 2 inputs — Class 3's `passing_summary.csv` is
   header-only (zero data rows), which routes the cross_sequence
   aggregator through `cross_sequence_summary.py:341–436`'s tier-none
   fallback (calling `extract_passing.extract_row` on a row sourced
   from `aggregated_results.csv` rather than from a populated
   `passing_summary.csv`). Class 2 routes through the normal
   passing_summary read path. A characterization fixture covering only
   Class 2 will not pin the fallback path's output shape.
2. The fixture cost is negligible: parent `design_0` is already curated
   for the Rosetta and MPNN fixtures, so adding `design_0_seq_0` adds
   only the per-sequence workdir under
   `tests/negative_steering/data/runs/design_0_seq_0/` — and that
   workdir is the small "no reversion ran" shape (no `cycle_0/steered/`
   subtree, no `reversion_*.json`, no `contamination_scratch/`),
   comparable in size to the Class-1 candidates (~21 MB observed for
   design_28_seq_1; expect similar).

**Class 3b skipped.** The 6 `has_pose_collapses` Class 3 sequences are
covered for path-shape purposes by candidate #9 (`design_27_seq_0`,
straight Class 7) plus candidate #4 (Class 3a). Adding a third would
test the same combination of code paths — a sequence with a mix of
`no_reversion` and `pose_collapses` rows where the cross_summary
representative-selection logic happens to pick a `no_reversion` row.
That is a representative-selection artefact, not an additional
production code path. Skip unless the human disagrees.

Plus, for Classes 3a and 3b, the `cross_summary` representative
columns alone don't distinguish them from each other — the
discrimination requires reading the full per-sequence
`aggregated_results.csv` and counting `pose_collapses` rows. The
discriminator table above includes the verification-command form.

For each candidate, the smallest set of files needed to seed the
`tests/negative_steering/data/` fixture set is the same for all of them,
because the test workflow consumes by glob:

- `data/fastas/<seq_name>.fasta` — the MPNN-corrected receptor + effector
  FASTA. Source: `runs/<seq_name>/inputs/{receptor,effector}.fasta`
  combined back into one file, OR sourced from the `sequences/top_fastas/`
  output of MPNN.
- `data/design_pdbs/<design_stem>.pdb` — the parent RFDiffusion design
  PDB (chain A = Cα-only designed receptor, chain B = native effector).
  Source: `rfdiffusion/passing/<design_stem>.pdb`. Multiple sequences
  share a parent design, so dedupe.
- `data/rfdiffusion_metrics.json` — the metrics blob with per-design
  contact_residues / position_order / per_design_design_residues.
  Source: `rfdiffusion/rfdiffusion_metrics.json`. One file shared by
  all candidates.
- `data/input_complex.pdb` — the receptor+effector complex used for
  controls. Already on disk at `tests/full_test_run/af3_pikp1_native_avrpikf_complex.pdb`.

The proposed 8-sequence fixture (candidates 1, 3, 4, 5, 6, 7, 8, 9
above — i.e. all required candidates including the new Class 3 and
Class 2 additions) covers parent designs `design_28, design_62,
design_0, design_15, design_42, design_3, design_55, design_27`
(8 distinct parent designs). Dropping the optional candidates 2 and 10
keeps the count at 8 across 8 parents — one example per observable
class plus one Class 3a fixture for the empty-passing_summary fallback
path.

#### Verification commands

```bash
# (a) Confirm cold_start_all_clean fired for design_28_seq_1; expect "true"
python3 -c '
import json
p=json.load(open("negative_steering/runs/design_28_seq_1/cycle_0/plan.json"))
print("cold_start_all_clean=", p.get("cold_start_all_clean"))
print("skip_steering=", p.get("skip_steering"))
print("cold_start_n_clean=", p.get("cold_start_n_clean"), "of", p.get("num_seeds"))'

# (b) Confirm pose_holds 3/3 for design_42_seq_0; expect 3 pose_holds, 0 collapses
python3 -c '
import csv
for r in csv.DictReader(open("negative_steering/runs/design_42_seq_0/aggregated_results.csv")):
    if r.get("aggregated_verdict")=="pose_holds":
        print("verdict=", r["aggregated_verdict"],
              "n_pose_holds=", r["n_seeds_pose_holds"],
              "n_pose_collapses=", r["n_seeds_pose_collapses"],
              "n_no_data=", r["n_seeds_no_data"],
              "n_clean_steered=", r["n_seeds_clean_steered"],
              "design=", r.get("design",""))'

# (c) Confirm pose_holds 2/3 for design_3_seq_1; expect n_pose_holds=2 n_pose_collapses=1
python3 -c '
import csv
for r in csv.DictReader(open("negative_steering/runs/design_3_seq_1/aggregated_results.csv")):
    if r.get("aggregated_verdict")=="pose_holds":
        print("design=", r.get("design"), "n_pose_holds=", r["n_seeds_pose_holds"],
              "n_pose_collapses=", r["n_seeds_pose_collapses"])'

# (d) Confirm pose_collapses aggregated for design_27_seq_0; expect a row with verdict pose_collapses
python3 -c '
import csv
verdicts=set()
for r in csv.DictReader(open("negative_steering/runs/design_27_seq_0/aggregated_results.csv")):
    verdicts.add(r.get("aggregated_verdict",""))
print("design_27_seq_0 verdicts:", sorted(verdicts))'

# (e) Confirm reversion was actually attempted (n_contaminated>0) for the pose_holds and pose_collapses candidates
for s in design_42_seq_0 design_3_seq_1 design_27_seq_0; do
  python3 -c "
import json
d=json.load(open('negative_steering/runs/$s/cycle_0/contaminated.json'))
print('$s', 'n_checked=', d.get('n_checked'), 'n_contaminated=', d.get('n_contaminated'))"
done
# Each should print n_contaminated > 0; n_contaminated == 0 means reversion was correctly skipped
# (would refute a pose_holds / pose_collapses claim).

# (f) Search the whole cohort for any aggregated new_contamination verdict; expect zero hits
python3 -c '
import csv, glob
hits=[]
for f in glob.glob("negative_steering/runs/*/aggregated_results.csv"):
    for r in csv.DictReader(open(f)):
        if r.get("aggregated_verdict")=="new_contamination":
            hits.append((f.split("/")[-2], r.get("design",""), r.get("aggregated_verdict_reason","")[:80]))
print("aggregated new_contamination occurrences:", len(hits))
for h in hits[:5]: print(" ", h)'
# Expect "0". A non-zero count means we have a missed candidate for class 8 — surface it.

# (g) For class 3 disambiguation: list per-design n_seeds_no_data dominance for tier-none/no_reversion sequences
python3 -c '
import csv, glob, sys
# Load cross_summary to get tier-none/no_reversion sequences
cross=list(csv.DictReader(open("negative_steering/cross_sequence_summary.csv")))
target=[r["mpnn_sequence"] for r in cross
        if r["cross_tier"]=="none" and r["representative_aggregated_verdict"]=="no_reversion"
        and r["representative_design"]!="initial"]
print(f"tier=none/no_reversion/non-initial sequences: {len(target)}")
print()
print("Per-sequence: total_seeds_clean_steered vs total_seeds_no_data across all designs")
for seq in target[:6]:
    p="negative_steering/runs/"+seq+"/aggregated_results.csv"
    n_cs=n_nd=n_designs=0
    for r in csv.DictReader(open(p)):
        ncs=int(r.get("n_seeds_clean_steered","0") or 0)
        nnd=int(r.get("n_seeds_no_data","0") or 0)
        if ncs+nnd>0: n_designs+=1
        n_cs+=ncs
        n_nd+=nnd
    print(f"  {seq}: n_designs_with_data={n_designs} sum_clean_steered={n_cs} sum_no_data={n_nd}")'
# Sequences with sum_clean_steered >> sum_no_data are "every seed clean_steered but ra_eff too high"
# Sequences with sum_no_data >> sum_clean_steered are "cold-start cascade collapse"

# (h) Stage-output sizes
du -sh negative_steering/ negative_steering/runs/ negative_steering/indices/
# Estimate per-sequence workdir size:
du -sh negative_steering/runs/design_28_seq_1/ negative_steering/runs/design_42_seq_0/ negative_steering/runs/design_27_seq_0/
```

### Orthogonal metrics

#### Paths

The orthogonal-metrics test fans out three per-survivor streams (AF3-no-MSA,
biophysical, Rosetta) and then merges them into
`survivors_with_orthogonal_metrics.csv`. `bin/merge_orthogonal_metrics.py`
applies four threshold gates and writes `passes_orthogonal_filters` (1 or 0)
plus a comma-separated `orthogonal_flags` list naming each gate that
failed:

- `af3_nomsa_ra_eff_too_high` — `af3_nomsa_best_ra_eff > orthogonal_filter_af3_ra_max` (default 5.0)
- `interface_plddt_too_low` — `interface_plddt < orthogonal_filter_plddt_min` (default 0.75)
- `sc_too_low` — `sc < orthogonal_filter_sc_min` (default 0.55)
- `bsa_too_low` — `bsa < orthogonal_filter_bsa_min` (default 600)
- `af3_nomsa_missing` — AF3 stream had no rows for the survivor

A "survivor" in the manifest sense (per `bin/extract_survivor_manifest.py`)
is **any** cross_summary row with a non-empty `representative_canonical_pdb`
and a workdir containing `plan.json` + a valid `ground_truth` PDB +
effector_template_cif. **The manifest does NOT gate on `cross_tier`** — so
tier-none rows are also fed into the orthogonal pipeline. In the discovery
run, all 122 cross_summary rows became survivors and went through AF3 +
biophys + rosetta.

The two paths the test workflow exercises are:

- **survivor present in manifest, all four metrics computed** — appears in
  `survivors_with_orthogonal_metrics.csv` with non-blank
  `af3_nomsa_*`, `bsa`, `sc`, `interface_plddt` columns.
- **survivor passes all gates** — `passes_orthogonal_filters == 1`, empty
  `orthogonal_flags`. **NOT OBSERVED in this discovery run** (0/122).
- **survivor fails one or more gates** — `passes_orthogonal_filters == 0`,
  `orthogonal_flags` lists the failing gates.

#### Discriminating signals

| Path | File | Field |
|---|---|---|
| in manifest | `orthogonal_metrics/survivor_manifest.csv` | row with `seq_name == <seq>` |
| AF3 stream computed | `orthogonal_metrics/af3_nomsa_summary.csv` | row present |
| biophys stream computed | `orthogonal_metrics/biophysical_summary.csv` | row present |
| rosetta stream computed | `orthogonal_metrics/rosetta_summary.csv` | row present |
| passes all gates | `orthogonal_metrics/survivors_with_orthogonal_metrics.csv` | `passes_orthogonal_filters == "1"` AND `orthogonal_flags == ""` |
| fails AF3 ra | same | `orthogonal_flags` contains `af3_nomsa_ra_eff_too_high` |
| fails plddt | same | contains `interface_plddt_too_low` |
| fails sc | same | contains `sc_too_low` |
| AF3 absent | same | contains `af3_nomsa_missing` |

#### Proposed candidates

The orthogonal-metrics test takes a **completed negsteer run as input**
(`data/negsteer_run/cross_sequence_summary.csv` + per-sequence
`runs/<seq>/`). Curating a fixture means selecting a small set of
sequences whose per-sequence workdirs we keep. Since the orthogonal
streams are per-survivor, the fixture should include:

- **design_28_seq_1** — top tier-A; representative-canonical-PDB present
  and workdir intact. Will exercise the AF3 + biophys + rosetta + merge
  cascade end-to-end with high-quality inputs. Expected output:
  `passes_orthogonal_filters == 0`, `orthogonal_flags ==
  "af3_nomsa_ra_eff_too_high:32.86"`.
- **design_62_seq_0** — closest-to-passing in the discovery run
  (`af3_ra=14.53`, only one flag); useful as the "most-likely-to-pass"
  example. May be worth re-running with adjusted thresholds in a future
  test to exercise the pass path; for now, exercises the
  near-boundary case of `af3_nomsa_ra_eff_too_high`.

Two survivors is enough to exercise the full stream cascade and the merge
gate logic. Adding a third (e.g., `design_42_seq_0`) would let the test
cover the `pose_holds` representative path through the orthogonal layer
as well.

For each survivor, the per-sequence workdir bundle the test consumes is
substantial:

- `runs/<seq>/cycle_0/initial_prediction.pdb` (the canonical PDB for
  cold-start representatives) OR `runs/<seq>/cycle_0/steered/<label>/prediction.pdb`
  (for steered/pose_holds representatives)
- `runs/<seq>/cycle_0/effector_template.cif`
- `runs/<seq>/cycle_0/plan.json` (provides ground_truth path,
  effector_template_cif path, etc.)
- `runs/<seq>/aggregated_results.csv` and per-cycle CSVs (consumed by
  EXTRACT_SURVIVOR_MANIFEST + NEGSTEER_INTERFACE_METRICS)
- the input cross_summary CSV, copied to
  `data/negsteer_run/cross_sequence_summary.csv`

Plus, the test reads `cross_sequence_summary.csv` itself — but only
rows whose workdir is present become survivors. So the curated fixture
needs the cross_summary trimmed to only the sequences whose workdirs
are present.

#### Verification commands

```bash
# (a) Survivor count (should equal cross_summary rows whose workdir + plan.json + ground_truth all exist)
wc -l orthogonal_metrics/survivor_manifest.csv
# Subtract 1 for header. In this run, expect ~122 (all cross_summary rows became survivors).

# (b) Confirm zero pass; print the flag distribution
python3 -c '
import csv
from collections import Counter
rows=list(csv.DictReader(open("orthogonal_metrics/survivors_with_orthogonal_metrics.csv")))
print("total survivors:", len(rows))
print("passes_orthogonal_filters dist:", Counter(r["passes_orthogonal_filters"] for r in rows))
flag_count=Counter()
for r in rows:
    for f in r.get("orthogonal_flags","").split(","):
        f=f.split(":")[0].strip()
        if f: flag_count[f]+=1
print("flag occurrences:", dict(flag_count))'
# Expect "passes:0=122" and four gate names with counts.

# (c) Confirm canonical PDB exists for the two top candidates
for s in design_28_seq_1 design_62_seq_0; do
  pdb=$(python3 -c "
import csv
for r in csv.DictReader(open('negative_steering/cross_sequence_summary.csv')):
    if r['mpnn_sequence']=='$s':
        print(r['representative_canonical_pdb']); break")
  echo "$s -> $pdb"
  ls -la "$pdb" 2>&1 | head -1
done

# (d) For design_28_seq_1, confirm only one orthogonal flag set; expect af3_nomsa_ra_eff_too_high:32.86
python3 -c '
import csv
for r in csv.DictReader(open("orthogonal_metrics/survivors_with_orthogonal_metrics.csv")):
    if r["mpnn_sequence"]=="design_28_seq_1":
        print("flags:", r["orthogonal_flags"])
        print("af3_ra:", r["af3_nomsa_best_ra_eff"], "sc:", r["sc"], "bsa:", r["bsa"], "iface_plddt:", r.get("interface_plddt_median",""))'

# (e) Stage-output sizes
du -sh orthogonal_metrics/
# And per-survivor workdir, since the orthogonal-metrics test consumes them:
du -sh negative_steering/runs/design_28_seq_1/ negative_steering/runs/design_62_seq_0/
```

## Coverage gaps and risks

1. **RFDiffusion `passes_filter == false` is not exercised** by the
   discovery run, and the human has confirmed it has not been observed
   in any pipeline run, real or test. The Cα-contact filter is a
   defensive code path against a degenerate input the upstream diffusion
   does not produce. **Status: not a characterization gap.** This
   branch will be closed by a small unit-style test added during Phase
   2.9 (characterization-test restructuring), not by a fixture from
   any pipeline run. No re-run of RFDiffusion is needed for fixtures.

2. **Aggregated AND per-seed `new_contamination` verdicts both not
   produced.** Verification command (f) on HPC found zero sequences
   with any per-seed `n_seeds_new_contamination > 0`, in addition to
   the zero aggregated `new_contamination` rows already known. This
   means the discovery run produced no examples of the
   `new_contamination` code path at either resolution level. The
   notes11 gating change (narrowing `new_contamination` to "majority
   contact at gated positions" only) appears to have moved this from
   "rare" to "absent" on this scaffold. Per-module fixtures from this
   discovery run will not cover the `new_contamination` code path at
   all. **Status: deferred to a standalone task the human will pick
   up tomorrow morning.** The eventual fix will likely be a small
   manufactured fixture or a targeted unit test against
   `bin/reversion.py::classify_reversion_verdict` and
   `bin/boltz2_iterate_steering.py::_classify_aggregated_verdict`,
   not a re-run of the full discovery pipeline.

3. **`passes_orthogonal_filters == 1` not produced** (0/122 survivors
   passed). The flag-frequency breakdown from HPC verification (b):
   `af3_nomsa_ra_eff_too_high` 120, `interface_plddt_too_low` 70,
   `sc_too_low` 4, `af3_nomsa_missing` 2. The `af3_nomsa_ra_eff` gate
   in particular is **warning-style**, not a hard reject — AF3-no-MSA
   places effectors imprecisely in general, so flagging at the 5 Å
   threshold reflects a property of AF3-no-MSA's accuracy floor, not
   a defect in the candidate sequences. The 0/122 pass rate is
   expected behaviour. **Status: not a characterization gap.** The
   pass branch of `merge_orthogonal_metrics.py` remains exercisable
   in the per-module test by overriding the threshold params (e.g.
   `orthogonal_filter_af3_ra_max=50.0`) without changing input
   fixtures; flag this in the test's README so a future reader
   doesn't misinterpret the 0-pass result.

4. **Class 3 split confirmed by g-deep verification: 63 / 6.**
   Verification of all 69 Class 3 sequences showed the underlying
   `aggregated_results.csv` files split cleanly into 63 sequences
   with all `no_reversion` design rows (Class 3a) and 6 sequences
   with at least one `pose_collapses` design row (Class 3b). The
   "all_zeros" appearance in per-seed columns is an artefact of the
   producer blanking those columns when verdict is `no_reversion`
   (`boltz2_iterate_steering.py:4705–4709`); it does NOT mean
   "steering produced no data". The shape was previously called a
   coverage risk because the cross_summary representative didn't
   reveal the underlying split — but the candidate-list update above
   adds a Class 3a fixture (`design_0_seq_0`) that exercises the
   empty-`passing_summary.csv` + tier-none-fallback code path the
   other classes don't. Class 3b is left uncovered as a
   representative-selection artefact (the same code paths fire as
   for Class 7 + Class 3a together). **Status: addressed by adding
   `design_0_seq_0` to the candidate list.**

5. **`representative_design == "initial"` only fires in the
   `cold_start_all_clean` path** in this discovery run. The
   `_agg_row_is_clean_steered` function in `boltz2_iterate_steering.py:4749`
   will admit a single-row "initial" baseline when `aggregated_verdict ==
   "no_reversion"` and the legacy ra_eff cap holds — but no row in the
   discovery run takes that path independently of `cold_start_all_clean`.
   Not a coverage gap for the fixtures (the path is exercised), just
   worth noting that the two code paths share a representative shape.

6. **The 13 `cold_start_all_clean` sequences all skipped reversion
   entirely** — their per-sequence workdirs will lack
   `cycle_0/contaminated.json`, `cycle_0/reversion_*.json`, the
   `cycle_0/steered/` subtree, and the contamination_scratch outputs.
   The fixture file list for Class-1 candidates is therefore noticeably
   smaller than for pose_holds / pose_collapses candidates. Worth
   flagging in the curation prompt so the tarball construction logic
   doesn't expect those files for cold_start_all_clean sequences.

## Open questions for the human

1. **Tier rule disagreement: notes11 vs pipeline_notes3.**
   `notes/negsteer_notes/notes11.md` defines Tier A as
   `aggregated_verdict == "no_reversion" OR n_seeds_pose_holds == n_seeds`
   and the docstring at the top of `bin/cross_sequence_summary.py`
   reproduces that wording verbatim. But the **implementation** at
   `cross_sequence_summary.py::_tier_for_row` (line 200) uses
   `n_pass = n_seeds_pose_holds + n_seeds_clean_steered` as the
   pass-equivalent count. `pipeline_notes/pipeline_notes3.md` agrees with
   the implementation (clean_steered counted as a pass). So the
   docstring at the top of `cross_sequence_summary.py` is stale relative
   to the implementation — it says "n_seeds_pose_holds == n_seeds"
   instead of "n_pass == n_seeds". Worth a docstring fix in Phase 3,
   but **for taxonomy purposes the implementation wins** and
   `clean_steered` is counted toward Tier A. This report uses the
   implementation's definition.

2. **Notes11 on the `_classify_aggregated_verdict` "any-seed gate".**
   Notes 6–8 describe a strict majority-of-N gate; notes 10 documents
   the Bug-E relaxation to "any-seed-pass + failure-count guard"; notes
   11 confirms the gate as it currently stands. The producer
   implementation matches notes 11. No outstanding contradiction, but
   if the human reads notes 6/7 in isolation they'll get the wrong gate
   semantics — flagged here so it doesn't surprise.

3. **RESOLVED — `design_44_seq_1` is tier A via clean_steered, not via
   pose_holds.** HPC verification of `design_44_seq_1`'s
   `aggregated_results.csv` returned three pose_holds rows for that
   sequence:
   - `design_00`: `n_pose_holds=1, n_clean_steered=2` → `n_pass=3=n_seeds` → tier A
   - `design_02`: `n_pose_holds=2, n_pose_collapses=1, n_clean_steered=0` → tier B
   - `design_03`: `n_pose_holds=2, n_pose_collapses=1, n_clean_steered=0` → tier B

   The representative pick (best non-empty tier; within tier, lowest
   composite-score rank) lands on `design_00` because it's the only
   tier-A row in the sequence. The helper file's
   `representative_n_seeds_pose_holds=1` reflects design_00's value;
   `n_seeds_clean_steered=2` is the column it didn't show. Candidate
   stands.

4. **No `new_contamination` aggregated OR per-seed rows in the discovery
   cohort.** Verification command (f) returned 0 per-seed occurrences
   across all 122 sequences. The notes11 gating to design-region ∪
   true-interface positions has moved this from "rare" to "absent" on
   this scaffold. The human's plan: tackle this as a standalone task
   tomorrow (small manufactured fixture or unit test against
   `classify_reversion_verdict` + `_classify_aggregated_verdict`),
   not by re-running the discovery pipeline. See Coverage gap 2.

5. **Survivor manifest does not gate on `cross_tier`.** This means
   the orthogonal-metrics pipeline computes AF3 + biophys + rosetta
   for tier-none representatives too. In the discovery run that's 95
   AF3 invocations on tier-none sequences — fully expected to fail
   the `af3_nomsa_ra_eff_too_high` gate. This is not a bug per se,
   but worth confirming whether the human wants the production survivor
   manifest to filter on tier (saving GPU minutes) or wants to keep
   the diagnostic-completeness behaviour.

6. **`representative_canonical_pdb` paths in `cross_sequence_summary.csv`
   were rewritten by `_rewrite_workdir_path_to_published`** to point at
   the publishDir path under `negative_steering/runs/<seq>/...`. This
   means the path values in the helper file refer to the HPC
   publishDir layout, not the workdir paths. When pulling fixtures
   to Mac, the curation prompt will need to either preserve that
   directory layout or rewrite the paths. Flag in the curation prompt.

## Verification results applied

The verification commands a–h listed in the per-stage sections were
run on HPC by the human (output pasted into the session, not into the
report). The following table logs which results changed the report and
how. "VC=Verification command".

| Result | Surprised? | Report change |
|---|---|---|
| RFDiffusion VC: 64/64 pass, 0 fail | No (expected) | Coverage gap 1 reframed: human confirmed no real run has produced a fail; gap closed by Phase 2.9 unit test, not fixture |
| Rosetta VC: 3 fails (design_14/25/59), all four candidate PDBs exist | No | none |
| ProteinMPNN VC: `mpnn/top_metadata.csv` not found | Yes | Path corrected to `sequences/top_metadata.csv` throughout the MPNN section; verification command updated; root cause noted (`MPNN_SELECT_TOP` publishes to `sequences/`, not `mpnn/`) |
| Negsteer VC (a): `cold_start_all_clean=True, skip_steering=True, 3 of 3` for design_28_seq_1 | No | none — Class 1 candidate confirmed |
| Negsteer VC (b): design_42_seq_0 returns clean 3/0/0/0/0 → tier A | No | none — Class 4 candidate confirmed |
| Negsteer VC (b extra): design_44_seq_1 returns three pose_holds rows; representative is design_00 with 1 pose_holds + 2 clean_steered (n_pass=3) | Yes — resolved Open Q3 | Open Q3 marked resolved with the explanation that the helper file showed only one of the two contributing per-seed columns |
| Negsteer VC (c): design_3_seq_1 design_00 row n_pose_holds=2, n_pose_collapses=1 | No | none |
| Negsteer VC (d): design_27_seq_0 verdicts are `{'no_reversion', 'pose_collapses', 'singleton'}` | Yes — singleton was claimed absent | Axis 3 rewritten with the correct rule (singleton fires per cycle-0 baseline row, always present); observable-classes table updated; new singleton-presence row added to the discriminator table |
| Negsteer VC (e): n_contaminated > 0 for design_42_seq_0 (3), design_3_seq_1 (12), design_27_seq_0 (4) | No | none |
| Negsteer VC follow-up on design_55_seq_1: pose_holds row `design=design_02, n_pose_holds=1, n_pose_collapses=0, n_no_data=2, n_clean_steered=0` → tier C; verdict-set `{'no_reversion', 'pose_collapses', 'pose_holds', 'singleton'}` (most heterogeneous in the cohort) | Useful — Class 6 candidate carries all four verdict types | Candidate #8 note rewritten to record the verdict-set heterogeneity as the candidate's value-add, replacing the prior "Verify on HPC before curating" caveat |
| Negsteer VC (f): zero per-seed `new_contamination` occurrences cohort-wide | Yes — was hoped "may still be > 0" | Coverage gap 2 reframed from "partial" to "total at both resolutions"; deferred to standalone task |
| Negsteer VC (g) → g-deep: 69 Class 3 sequences split 63 (`all_zeros`) / 6 (`has_pose_collapses`); aggregate verdicts 264 no_reversion + 69 singleton + 6 pose_collapses | Yes — Class 3 was treated as monolithic | Observable-classes table split Class 3 into 3a (63) and 3b (6); new candidate `design_0_seq_0` added to the candidate list as the Class 3a fixture; reasoning recorded under "Class 3 fixture decision" |
| Negsteer VC (h): per-sequence workdir sizes 21M (cold-start), 142M (pose_collapses), 147M (pose_holds); total `negative_steering/` is 13G | No | candidate notes annotated with workdir sizes so the curation prompt can budget the tarball footprint |
| Orthogonal VC (a): 122 survivors | No | none |
| Orthogonal VC (b): 0/122 pass; flag dist `af3_ra_too_high:120, plddt_too_low:70, sc_too_low:4, af3_missing:2` | Yes — was treated as a coverage gap | Coverage gap 3 reframed: af3_nomsa_ra_eff is a warning-style gate, 0/122 is expected, gap removed |
| Orthogonal VC (c): design_62_seq_0's representative_canonical_pdb is at `cycle_0/steered/design_02_s0/prediction.pdb` (not `initial_prediction.pdb`) | Useful — confirms a Class 2 fixture exists | New candidate `design_62_seq_0` added to the candidate list as both the Class 2 fixture and the orthogonal-metrics survivor (single workdir doubles for both tests) |
| Orthogonal VC (d): design_28_seq_1 flags = `af3_nomsa_ra_eff_too_high:32.86`, sc=0.732, bsa=913.64, iface_plddt blank | Useful (iface_plddt source col is named differently) | none — proceed; the `interface_plddt` column populated for the gate may have a different name; future verification can pin down |
