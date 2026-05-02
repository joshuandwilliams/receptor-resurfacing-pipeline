# Negative Steering — Reversion Validation Pass — Session Summary

## Context

Entering this session, the pipeline had been validated on four native benchmark
complexes (notes 3) and two RFDiffusion design runs (notes 4–5). The open question
from notes 5 was whether the high-confidence steered predictions for Pikp1 design_3
were trustworthy as wet-lab candidates, given that the steering mutations were
landing inside or adjacent to the protected design region and forming contacts in
the rescued pose. The mutation-reliance check added at the end of notes 5 flagged
51 of 57 passing designs as contaminated — but that check was using mutation-position
intersection rather than design-region intersection, which was too broad.

This session had two goals: (1) build the full reversion validation pass — a
per-cycle mechanism that, for any steered design whose rescue depends on mutations
at protected-set positions, reverts those mutations and re-runs Boltz to determine
whether the pose holds without them, and (2) run it end-to-end on a fresh design_3
experiment to get the first real data on how many steered candidates actually survive.

## What we built

### 1. Refined mutation-reliance check: protected-set intersection

The notes-5 check flagged any design whose contact residues overlapped the mutated
positions. This over-flagged receptor-resurfacing designs because the steering
deliberately picks mutations adjacent to the true interface.

The corrected check asks whether any mutated contact residue falls inside the
**protected set**, defined as:

```
protected_set = design_region_positions ∪ true_interface_positions
```

`design_region_positions` comes from `rfdiffusion_metrics.json`'s
`per_design_design_residues` field, mapped through `receptor_position_order` to
1-based positional indices via the new `derive_design_region.py` helper.
`true_interface_positions` comes from `cycle_0/plan.json`'s `true_interface_idx`
field (0-based, converted to 1-based at the boundary).

The union is resolved inside `cmd_compute_final_metrics`, so the user supplies only
the design-region file and the orchestrator unions it with plan.json's true interface
automatically.

New columns emitted by `compute_metrics.py`:
- `n_contact_on_protected_positions`
- `n_mutated_contact_on_protected_positions`
- `mutated_contact_protected_positions` (comma list)
- `construct_reliance_flag` (`clean` / `contaminated` / `""`)

The old `mutation_reliance_flag` remains for backward compat but is now a noisier
diagnostic; `construct_reliance_flag` is the operationally meaningful one.

The key insight from applying the corrected check to the previous design_3 CSV: the
distinction between position 44 (inside design region 33-46 → protected → flag)
and position 26 (outside design region, outside true interface → not protected →
don't flag) is what separates the c0d40 family from the cycle-0 d40 parent. Under
the corrected check, cycle-0 design_40 (A26R Q46D) is clean; all its c0d40.c1dXX
children are contaminated because Y44X contacts appear in the rescued poses.

### 2. `derive_design_region.py`

New standalone stdlib-only script. Takes `rfdiffusion_metrics.json` + a design name,
derives 1-based positional indices of the de novo region, writes them as a collapsed
range string (e.g. `33-46,74-79`), and optionally cross-checks against
`cycle_0/plan.json`'s `true_interface_idx` to confirm coordinate-system consistency.

Usage:
```bash
python3 derive_design_region.py \
    --metrics-json resurface_pipeline_test/rfdiffusion_metrics.json \
    --design design_3 \
    --output resurface_pipeline_test/design_3_design_region.txt \
    --cross-check-plan resurface_pipeline_test/runs/<run>/cycle_0/plan.json
```

For Pikp1 design_3 the output is `33-46,74-79` (two de novo blocks, 20 positions
total). The cross-check WARNING is expected and correct: design_3's true interface
includes 6 native-scaffold anchor residues (3, 48, 70-73) that RFDiffusion kept
outside the de novo region. Those anchors are covered by the union in the orchestrator
even though they're not in the design-region file.

### 3. `--design-region-indices-file` flag in `cmd_plan`

Added to `boltz2_negative_steering.py`. Accepts a file in the same range/list format
that `derive_design_region.py` produces. Converts from 1-based (file convention) to
0-based and stores `design_region_idx` in `cycle_0/plan.json` alongside the existing
`true_interface_idx`. Future runs will read `design_region_idx` from plan.json
automatically; existing runs can use the `--design-region-positions-file` override in
`compute-final-metrics`.

### 4. Full reversion validation pass

The reversion pass runs at the end of every cycle (cycle 0 and cycle 1+) as a chain
of new SLURM stages inserted between the existing collect stage and the
novelty/maximin selection. The full per-cycle flow is now:

```
compute-distances (container)
  → iterate-collect-prefilter (host): intact filter → writes prefilter.json
  → build-contaminated (container): contact analysis → writes contaminated.json
  → plan-reversions (container): stages reverted YAMLs → writes reversion_plan.json
  → predict-reversion-one × N (GPU array): Boltz on reverted sequences
  → harvest-reversions (container): metrics on reverted predictions
  → iterate-collect-finalize (host): verdict routing → novelty → maximin → passing.json
```

Cycle 0 uses `kickoff-prefilter` / `kickoff-finalize` in the same pattern.

**Verdict logic in finalize:**

| condition | verdict | outcome |
|---|---|---|
| Reverted pose passes structural filter AND `construct_reliance_flag == clean` | `pose_holds` | Design kept; `receptor.fasta` overwritten with reverted sequence so cycle N+1 starts from the validated sequence |
| Reverted ra_eff ≥ threshold OR receptor not intact | `pose_collapses` | Design dropped from `passing.json`; kept in `all_results_multicycle.csv` with verdict |
| Reverted pose passes structural filter BUT `construct_reliance_flag == contaminated` | `new_contamination` | Design dropped; kept in CSV |
| Design not flagged contaminated by `build-contaminated` | (no verdict) | Passes through unchanged |

The ordering is: structural filter → reversion → repeat structural filter on reverted
poses → novelty → maximin. This ensures the novelty selector operates on the
validated-sequence versions of surviving designs, not on the steered sequences.

**Key design decision — destructive FASTA rewrite:** For `pose_holds` designs, the
design's `receptor.fasta` is overwritten with the reverted sequence. The original
steered sequence is preserved implicitly in the reversion workdir
(`cycle_N/pathway_<label>/reversions/<label>/receptor.fasta`). The existing
`mutations.tsv` is left intact (it records that the position was touched, which is
correct for the burnt-positions rule — cycle N+1 should not re-mutate a position
that was already edited and reverted, since that would just recreate the contamination).

### 5. New subcommands

| subcommand | runs where | purpose |
|---|---|---|
| `iterate-collect-prefilter` | host | intact filter; writes `prefilter.json` |
| `build-contaminated` | container | contact analysis per intact candidate; writes `contaminated.json` |
| `plan-reversions` | container | stages reverted Boltz YAMLs; writes `reversion_plan.json` |
| `predict-reversion-one` | container (GPU) | Boltz on one reversion YAML by array index |
| `harvest-reversions` | container | metrics on reverted predictions; writes `reversion_results.json` |
| `iterate-collect-finalize` | host | verdict routing + novelty + maximin + resubmit |
| `kickoff-prefilter` | host | cycle-0 equivalent of `iterate-collect-prefilter` |
| `kickoff-finalize` | host | cycle-0 equivalent of `iterate-collect-finalize` |

### 6. `reversion.py`

New module with four pure-logic functions:

- `build_reverted_sequence(steered_fasta_path, cumulative_mutations, positions_to_revert) → str`
  Reverts exactly the specified positions to their earliest wild-type residue from the
  mutation history. Includes a sanity check that verifies the current residue in the
  FASTA matches the most-recent edit in the history before reverting. This check caught
  the off-by-one label bug (see Bug history) correctly — it hard-failed rather than
  silently reverting the wrong positions.
- `write_reversion_plan(workdir, contaminated_designs, plan_meta) → Path`
- `harvest_reversion_results(workdir, ...) → Dict`
- `classify_reversion_verdict(steered_metrics, reverted_metrics, ...) → str`

### 7. Submit script restructuring

Both `submit_boltz2_negative_steering.sh` (cycle 0) and
`submit_boltz2_iterate_steering.sh` (cycle 1+) now chain six SLURM jobs where there
were previously three:

Old: plan → predict-array → collect (compute-distances + iterate-collect in one wrap)

New: plan → predict-array → collect-prefilter → reversion-prep (build-contaminated +
plan-reversions) → reversion-predict-array → finalize (harvest-reversions +
iterate-collect-finalize)

The reversion array is sized to `ARRAY_MAX` (same as the main predict array); tasks
whose index exceeds the number of contaminated designs exit 0 immediately. If no
designs are contaminated, the reversion prep produces an empty `reversion_plan.json`
and all array tasks are no-ops.

### 8. `cmd_aggregate` reversion columns

`all_results_multicycle.csv` now carries per-row reversion data:

- `reversion_verdict` — `pose_holds` / `pose_collapses` / `new_contamination` / blank
- `reversion_dropped` — 1 if design was excluded from passing.json, else 0
- `reverted_sequence_used_as_parent` — 1 if receptor.fasta was rewritten, else 0
- `steered_ra_eff_vs_truth` — original pre-reversion metric (preserved even for kept designs)
- `steered_independent_receptor_rmsd`
- `steered_independent_effector_rmsd`
- `reverted_ra_eff_vs_truth`
- `reverted_independent_receptor_rmsd`
- `reverted_independent_effector_rmsd`
- `reverted_construct_reliance_flag`

Both steered and reverted metrics appear on every row that went through the verdict
phase, including dropped designs, so the CSV is a complete audit trail.

## Experiments and results

### design_3 reversion test — first successful end-to-end run

**Configuration:** `--mode mild --max-mutations 6 --candidate-pool-size 6
--n-designs 50 --n-cycles 2 --max-passing 3`

Ground truth: `design_3.pdb` (RFDiffusion Cα-only backbone, chains A/B)
Receptor FASTA: `design_3_seq_3_receptor.fasta` (79 residues, ProteinMPNN seq 3)
Protected set: 26 positions — design region (33-46, 74-79) ∪ native anchors (3, 48, 70-79)

**Cycle 0 results:**

| metric | value |
|---|---|
| generated | 50 |
| intact | 36 |
| contaminated | 0 |
| pose_holds | 0 (no reversion needed) |
| selected for cycle 1 | 5 (c0d00, c0d03, c0d09, c0d18, c0d36) |
| best ra_eff (design_38) | 1.15 Å |

0 of 36 intact designs were contaminated. All 5 selected cycle-0 survivors passed
through without reversion.

**Cycle 1 results (all 5 pathways, 125 total rows in aggregate):**

| verdict | count |
|---|---|
| pose_holds | 39 |
| pose_collapses | 30 |
| new_contamination | 2 |
| no verdict (cycle-0 rows) | 54 |

**Top wet-lab candidates (pose_holds, sub-2 Å reverted ra_eff):**

| pathway | ra_eff (steered) | verdict | notes |
|---|---|---|---|
| cycle_0 design_38 | 1.15 Å | — (c0, no reversion needed) | Primary candidate |
| c0d03.c1d41 | 1.24 Å | pose_holds | Best reversion survivor |
| c0d03.c1d21 | 1.50 Å | pose_holds | |
| c0d03.c1d31 | 1.80 Å | pose_holds | |

All three cycle-1 survivors come from the c0d03 parent. The c0d00 pathway produced
zero pose_holds survivors — all 15 intact cycle-1 children of c0d00 collapsed when
reverted, with median reverted ra_eff ~28.4 Å (essentially returning to cold-start).

### The c0d00 vs c0d03 divergence

This is the central finding of the session and the clearest validation that the
reversion pass is working as intended.

Both c0d00 and c0d03 were selected as cycle-0 survivors with comparable steered
ra_eff. Their cycle-1 children have ra_eff in the 3–5 Å band and both look
structurally promising. The reversion pass cleanly separates them:

**c0d00 pathway:** 15 intact cycle-1 children, 15/15 contaminated, 15/15
pose_collapses. Steered ra_eff median ~3.8 Å, reverted ra_eff median ~28.4 Å.
The cycle-1 mutations on this pathway were doing all the work of rescuing the pose.
Without those mutations, Boltz returns to its cold-start prior. Without the reversion
pass, all 15 of these would have been presented as promising wet-lab candidates.

**c0d03 pathway:** multiple pose_holds survivors. Cycle-1 mutations on c0d03's
descendants are non-load-bearing at protected positions — when reverted, Boltz still
finds the correct binding mode because the receptor's inherent sequence (the actual
ProteinMPNN design) is already sufficient to direct the prediction to the right site.

The distinction matters enormously: sending c0d00 descendants to the wet lab would
have been a waste of synthesis budget on sequences that Boltz does not predict will
bind when expressed as designed.

## Bug history

### Off-by-one in cycle-1+ label construction

`cmd_build_contaminated` and `cmd_iterate_collect_finalize._label_for` were calling
`make_pathway_label(parent_label, cycle - 1, design_idx)` instead of
`make_pathway_label(parent_label, cycle, design_idx)`. For cycle 1, this produced
labels like `c0d00.c0d04` (two cycle-0 segments) instead of the correct `c0d00.c1d04`.

When `read_cumulative_mutations` was called with the bogus label, it parsed two
cycle-0 entries and walked back through two unrelated cycle-0 designs' mutation
files. The cumulative history it returned was design_00's mutations combined with
design_04's mutations — completely unrelated to what cycle-1 design_04 actually had
in its FASTA.

`build_reverted_sequence`'s sanity check correctly caught this and hard-failed:
`steered sequence mismatch at position 36: history says 'R' but FASTA has 'L'`.
This caused plan-reversions to skip every contaminated design and produce an empty
`reversion_plan.json`, which cascaded to empty reversion_results.json and incorrect
pass-through of all contaminated designs in finalize.

Fixed by dropping the `- 1` from both label constructions.

**Note:** the sanity check in `build_reverted_sequence` working correctly here is
important. If it had been permissive, the bug would have silently built garbage
reverted sequences and Boltz would have predicted nonsense — probably all appearing
as pose_collapses by coincidence, which would have been a hard-to-diagnose false
negative. The hard-fail was the right behaviour.

### Cycle-0 `plan.json` schema mismatch

`cmd_build_contaminated` and `cmd_harvest_reversions` were written assuming cycle-1+
plan.json schema (which carries `experiment_root`). Cycle-0 plan.json (written by
`boltz2_negative_steering.py::cmd_plan`) does not have that field.

Fixed by detecting `plan.get("cycle", 0)`. For cycle 0, `experiment_root` is derived
from `workdir.parent` and `cycle0_plan` is the workdir's own plan.json. For cycle
1+, the existing path is unchanged.

### `build-contaminated` silently skips mutation check for cycle-0 designs

In `cmd_build_contaminated`'s cycle-0 branch, `label = label_suffix` produces
`"design_11"`. This gets passed to `read_cumulative_mutations(experiment_root,
"design_11")`, which calls `parse_pathway_label("design_11")`, which raises
`ValueError` (expects `cNdMM` format). The `except Exception` around the call
silently catches it, sets `cumulative = []`, gives `mutated_cli = ""`, and omits
`--mutated-positions` from the `compute_metrics.py` call. Without that argument,
`construct_reliance_flag` is never computed — every cycle-0 design comes back clean.

Result: the entire cycle-0 reversion pass was a no-op. design_11 (Y44E contacting
effector at protected position 44) and design_38 (L36R at position 36) were never
reverted. This was caught by visually inspecting the `cm_construct_reliance_flag`
column in `compute-final-metrics` output and noticing it showed contaminated when the
reversion pass had reported none.

Fix: read `mutations.tsv` directly from `pdb.parent` for cycle-0 designs rather than
going through `read_cumulative_mutations`.

### Reverted confidence metrics not forwarded to aggregate CSV

`harvest-reversions` was computing `reverted_ipsae_min`, `reverted_iptm`,
`reverted_complex_plddt` etc. and storing them in `reversion_results.json`. But the
finalize functions only copied four structural fields from `rev` into the candidate
dict. The confidence metrics never reached `passing.json` or the aggregate CSV.

This is critical because `reverted_ipsae_min` is the correct ranking signal for
wet-lab triage — not `cm_ipsae_min` (steered ipSAE), which includes contributions
from mutations that will not exist in the expressed construct.

Fix: added `_REVERTED_CONFIDENCE_FIELDS` and `_copy_reverted_confidence` helper.
Extended `reversion.py::harvest_reversion_results` to emit 11 reverted confidence
fields. Wired `_copy_reverted_confidence` into all six verdict sites across both
finalize functions. Updated `cmd_aggregate` to forward all 11 fields into the CSV.

### Submitting with `--n-cycles 1` silently skips reversion

The kickoff block in `submit_boltz2_negative_steering.sh` is gated by
`[[ "$N_CYCLES" -gt 1 ]]`. Submitting with `--n-cycles 1` produces a single-cycle
run but never fires the kickoff chain, so the reversion pass is silently skipped.
For reversion validation to run, always use `--n-cycles 2` or higher.

### Wrong ground-truth PDB (glycine placeholder sequences)

The first run used `af3_pikp1_native_avrpikf_complex.pdb` as ground truth and did
not pass `--receptor-fasta`. The plan stage extracted the receptor sequence from the
PDB's chain A, which is the RFDiffusion Cα-only backbone with glycine placeholders
in the de novo region (positions 33-40, 74-79 all G). Boltz then tried to fold a
poly-glycine stretch and produced a receptor backbone at 14.5 Å from the reference.

Fix: always pass `--receptor-fasta design_3_seq_3_receptor.fasta` for RFDiffusion
designs. The ground-truth PDB is used only for structural reference (RMSD measurement
and chain alignment), not as a sequence source. The MPNN FASTA is what Boltz folds.

## Key findings

### 1. The reversion pass works and is essential for receptor-resurfacing designs

39 pose_holds and 30 pose_collapses out of 71 cycle-1 designs that went through
verdict routing. 42% of structurally-plausible steered predictions collapse when
the load-bearing protected-position mutations are removed. Without reversion, all of
those would have been promoted to wet-lab candidates.

The c0d00 vs c0d03 divergence (15/15 collapse vs multiple survivors) shows that the
pass is discriminating between superficially similar pathways in a mechanistically
interpretable way.

### 2. Cycle-0 contamination rate was zero on this design

None of cycle 0's 36 intact designs were flagged contaminated. The cycle-0 candidate
pool (positions 26, 32, 36, 44, 46, 62) has 3 positions inside the design region
(36, 44, 46) and 3 outside (26, 32, 62). In practice, cycle-0 mutations either
landed outside the design region entirely, or landed inside it at positions that
didn't end up contacting the effector in the rescued pose. This meant the reversion
pass had nothing to do at cycle 0, and the 5 cycle-0 survivors passed through as-is.

This is consistent with the notes-3 observation that cycle-0 single-shot steering
is often sufficient to produce structurally correct poses — the contamination problem
accumulates as later cycles build on progressively more-mutated parents.

### 3. The `new_contamination` rate is low (2/41 post-reversion survivors)

Only 2 designs introduced new protected-position contacts after reversion that
weren't there in the original steered pose. Both were dropped. The 5% rate suggests
this edge case is rare, not systematic. Not enough data to draw conclusions about
what triggers it.

### 4. All three pose_holds sub-2 Å survivors come from the c0d03 parent

This suggests c0d03's cycle-0 mutations established a receptor configuration that
inherently supports the designed interface even without further steering. What makes
c0d03's receptor more self-sufficient than c0d00's is not yet clear — the mutations
are drawn from the same pool of 6 candidate positions. Worth checking what c0d03's
3 cycle-0 mutations actually were vs c0d00's.

## Run commands for this design class

Standard recipe for a Pikp1 RFDiffusion design:

```bash
# Step 1: derive design region file
python3 derive_design_region.py \
    --metrics-json resurface_pipeline_test/rfdiffusion_metrics.json \
    --design design_N \
    --output resurface_pipeline_test/design_N_design_region.txt

# Step 2: submit
./submit_boltz2_negative_steering.sh \
    --workdir         resurface_pipeline_test/runs/design_N_seq_M_v1 \
    --ground-truth    resurface_pipeline_test/design_N.pdb \
    --receptor        A \
    --effector        B \
    --receptor-fasta  resurface_pipeline_test/design_N_seq_M_receptor.fasta \
    --mode            mild \
    --max-mutations   6 \
    --candidate-pool-size 6 \
    --n-designs       50 \
    --n-cycles        2 \
    --true-interface-indices-file resurface_pipeline_test/design_N_true_interface.txt \
    --design-region-indices-file  resurface_pipeline_test/design_N_design_region.txt

# Step 3: aggregate (after all jobs complete)
python3 boltz2_iterate_steering.py aggregate \
    --experiment-root resurface_pipeline_test/runs/design_N_seq_M_v1

# Step 4: inspect results
python3 -c "
import json
p = json.load(open('resurface_pipeline_test/runs/design_N_seq_M_v1/cycle_0/passing.json'))
print(f'pose_holds={p[\"pose_holds_count\"]} collapses={p[\"pose_collapses_count\"]} new_contam={p[\"new_contamination_count\"]}')
"
```

**Critical:** always pass `--n-cycles 2` (not 1) or the kickoff chain never fires
and the reversion pass is silently skipped.

**Critical:** always pass `--receptor-fasta` pointing at the ProteinMPNN-designed
FASTA. The RFDiffusion PDB has glycine placeholders in the de novo region and is not
a valid sequence source.

## Files produced this session

- `boltz2_negative_steering.py` — `--design-region-indices-file` flag; writes
  `design_region_idx` into `plan.json`
- `boltz2_iterate_steering.py` — 8 new subcommands; `iterate-collect` and
  `kickoff` split into prefilter/finalize phases; `cmd_aggregate` reversion columns;
  cycle-0 schema-awareness in `cmd_build_contaminated` and `cmd_harvest_reversions`
- `reversion.py` — new module with `build_reverted_sequence`,
  `write_reversion_plan`, `harvest_reversion_results`, `classify_reversion_verdict`
- `compute_metrics.py` — `construct_reliance_flag` and protected-set intersection
  columns; `--protected-positions` CLI flag
- `derive_design_region.py` — new standalone helper
- `submit_boltz2_negative_steering.sh` — 6-stage chain (was 3)
- `submit_boltz2_iterate_steering.sh` — 6-stage chain (was 3)

## Additional bugs found in subsequent run

### `build-contaminated` silently skips mutation check for cycle-0 designs

After the first successful run appeared to show 0/36 cycle-0 designs contaminated,
`compute-final-metrics` (run post-hoc) flagged design_11 and design_38 as
contaminated. The two checks disagreed — same protected-set logic, different answers.

Root cause: in `cmd_build_contaminated`'s cycle-0 branch, `label = label_suffix`
produces `"design_11"` (the raw design name). This gets passed to
`read_cumulative_mutations(experiment_root, "design_11")`, which calls
`parse_pathway_label("design_11")`, which raises `ValueError` (expects `cNdMM`
format). The exception is caught silently, `cumulative` becomes `[]`,
`mutated_cli = ""`, and `--mutated-positions` is never passed to `compute_metrics.py`.
Without that argument, `construct_reliance_flag` can't be computed, so every cycle-0
design came back clean.

Fix: for cycle 0, read `mutations.tsv` directly from `pdb.parent` (the design
workdir is already in scope) rather than going through `read_cumulative_mutations`.
The `mutations.tsv` path is `pdb.parent / "mutations.tsv"` — straightforward stdlib
read, no pathway label needed.

Consequence for the run where this was found: the cycle-0 reversion pass was a
complete no-op. design_11 (Y44E contacting effector at protected position 44) and
design_38 (L36R contacting effector at protected position 36) were never reverted.
All cycle-1 designs were unaffected — their labels (`c0d03.c1d32` etc.) parse
correctly through `read_cumulative_mutations`, so contamination was properly detected
and reverted for cycle 1.

### `reverted_ipsae_min` and confidence suite not forwarded to aggregate CSV

After the reversion pass ran correctly for cycle 1, `harvest-reversions` was already
computing `reverted_ipsae_min`, `reverted_iptm`, `reverted_complex_plddt` and storing
them in `reversion_results.json`. However, `cmd_kickoff_finalize` and
`cmd_iterate_collect_finalize` only copied four structural fields from `rev` into the
candidate dict — `reverted_ra_eff_vs_truth`, `reverted_independent_receptor_rmsd`,
`reverted_independent_effector_rmsd`, `reverted_construct_reliance_flag`. The
confidence metrics never made it into `passing.json` or the aggregate CSV.

This was the critical gap: `reverted_ipsae_min` is the correct ranking signal for
wet-lab triage — Boltz's confidence in the interface using only the residues that
will actually exist in the expressed construct. Using `cm_ipsae_min` (steered
prediction) for ranking is systematically biased upward for contaminated designs,
because the steered prediction includes the contribution of mutations that will not
be present in the construct.

Fix: added `_REVERTED_CONFIDENCE_FIELDS` list and `_copy_reverted_confidence(c, rev)`
helper. Extended `reversion.py::harvest_reversion_results` to emit the full suite
(`reverted_ipsae_min/ab/ba`, `reverted_iptm`, `reverted_ptm`, `reverted_actifptm`,
`reverted_complex_plddt`, `reverted_avg_plddt`, `reverted_ipae`, `reverted_pae_mean`,
`reverted_pae_pass_frac`). Wired `_copy_reverted_confidence` into all six verdict
sites (pose_holds / pose_collapses / new_contamination × 2 finalize functions).
Updated `cmd_aggregate` to forward all 11 fields into the CSV for both cycle-0 and
cycle-1+ rows.

**Important note on `ipTM`:** `reverted_iptm` is the Boltz-native scalar from the
confidence JSON of the reverted prediction, not a recomputed value. `ipSAE` is the
more principled ranking signal because it derives from the PAE matrix over actual
interface pairs — but `reverted_iptm` is also available and moves consistently with
`reverted_ipsae_min` across candidates.

**Why not mask mutated positions from `cm_ipsae_min`?** A construct-aware ipSAE
computed by masking mutated rows from the steered prediction's PAE matrix would still
be contaminated — Boltz folded the structure with the mutations in place, so every
residue's PAE reflects a backbone that was conditioned on the mutations being there.
Masking rows from an already-mutation-conditioned PAE matrix gives a score for "the
non-mutated residues in a structure as predicted with mutations present" — which is
not the same as "confidence in binding if only the construct sequence existed."
The only clean answer is to run Boltz on the reverted sequence, which is exactly
what the reversion pass does. `compute_metrics.py` was therefore kept unchanged.

## Final results (third run: all bugs fixed, full reverted confidence suite)

**Cycle-0 contamination:** 32 of 36 intact designs were contaminated (previously 0
due to the `build-contaminated` label bug). Design_38 pose_collapses (reverted
ra_eff 5.03 Å — L36R was load-bearing, confirmed false positive). Design_11
pose_holds (reverted ra_eff 1.97 Å, reverted ipSAE 0.751 — Y44E contact was
incidental, pose moved only 0.5 Å on reversion).

**Final shortlist ranked by `reverted_ipsae_min`:**

| rank | pathway | cycle | steered ra_eff | reverted ra_eff | rev ipSAE | rev ipTM |
|------|---------|-------|---------------|----------------|-----------|---------|
| 1 | c0d03.c1d35 | 1 | 3.78 Å | 3.78 Å | 0.884 | 0.911 |
| 2 | c0d03.c1d33 | 1 | 2.14 Å | 2.14 Å | 0.862 | 0.905 |
| 3 | c0d03.c1d42 | 1 | 1.81 Å | 1.81 Å | 0.844 | 0.895 |
| 4 | c0d18.c1d21 | 1 | 3.40 Å | 3.40 Å | 0.838 | 0.899 |
| 5 | c0d03.c1d46 | 1 | 1.99 Å | 1.99 Å | 0.815 | 0.886 |
| 6 | c0d03.c1d45 | 1 | 2.24 Å | 2.24 Å | 0.805 | 0.886 |
| 7 | c0d03.c1d32 | 1 | 2.54 Å | 2.54 Å | 0.799 | — |
| 10 | cycle_0 design_11 | 0 | 1.44 Å | 1.97 Å | 0.751 | 0.866 |
| 11 | c0d03.c1d41 | 1 | 1.24 Å | 1.24 Å | 0.732 | — |

**Design_38 confirmed false positive.** The `low_pass_frac` flag on `cm_ipsae_min`
(0.016) was already catching the problem — steered confidence was near-zero despite
the 1.15 Å ra_eff. The reversion pass confirms it: L36R was the only thing holding
Boltz's effector at the correct site.

**c0d03.c1d35 is the top candidate by reverted ipSAE (0.884) despite the largest
ra_eff (3.78 Å) on the shortlist.** Reverted ra_eff exactly matches steered ra_eff,
meaning the pose does not move at all when contaminating mutations are removed. Boltz
is confident in the interface on the construct's own sequence. Worth visual inspection
in ChimeraX to confirm this is a valid alternate binding geometry rather than an
artefact.

**The ranking split between ra_eff and reverted ipSAE is the intended pipeline
output.** c0d03.c1d41 has the best structural score (1.24 Å) but ranks 11th on
reverted ipSAE (0.732). c0d03.c1d35 has the best confidence but a middling structural
score. Neither axis alone is sufficient; both are needed.

**For wet-lab triage:** rank by `reverted_ipsae_min` for contaminated designs,
fall back to `cm_ipsae_min` for uncontaminated designs (where the two are equivalent
— no mutations contacted the effector, so steered and reverted predictions are
identical).

**For wet-lab synthesis:** the construct sequence is the steered FASTA (not the
reverted one). `pose_holds` means "the pose doesn't depend on the load-bearing
mutations" — but the wet-lab construct carries the full steered sequence. The FASTA
to order is at `cycle_N/pathway_<label>/steered/design_XX/receptor.fasta`.

**ChimeraX top-5 paths for visual QC:**
```
# Reference
resurface_pipeline_test/design_3.pdb  (chains A=receptor, B=effector)

# Ranked by reverted_ipsae_min
cycle_1/pathway_c0d03/steered/design_35/prediction.pdb  (rev ipSAE 0.884)
cycle_1/pathway_c0d03/steered/design_33/prediction.pdb  (rev ipSAE 0.862)
cycle_1/pathway_c0d03/steered/design_42/prediction.pdb  (rev ipSAE 0.844)
cycle_1/pathway_c0d03/steered/design_46/prediction.pdb  (rev ipSAE 0.815)
cycle_1/pathway_c0d03/steered/design_32/prediction.pdb  (rev ipSAE 0.799)
# cycle-0 candidate if diversity wanted
cycle_0/steered/design_11/prediction.pdb                (rev ipSAE 0.751)
```

## Open work items

1. **ChimeraX visual QC on the final top-5.** Confirm c0d03.c1d35's 3.78 Å ra_eff
   is a valid alternate binding geometry and not a wrong pose at a coincidentally high
   confidence. Check reverted predictions alongside steered counterparts to confirm
   pose stability on reversion.

2. **Understand why c0d00 collapsed universally and c0d03 did not.** Read the
   cycle-0 mutations.tsv files for both designs and compare what positions they
   mutated and what residues they chose. The answer will inform whether candidate
   selection in the plan stage can be improved to prefer cycle-0 mutations that
   set up the receptor for self-sufficient binding.

3. **Run more designs.** design_3 is one RFDiffusion design from one run. The
   pipeline is now validated for this design class. Run it on design_1 (with the
   corrected MPNN sequence — see notes 4) and on higher-scoring designs from
   `scored_metadata.csv` to get a broader picture.

4. **Integrate reversion pass into the main Nextflow pipeline.** Currently the
   reversion validation is a standalone negative-steering step. Once the pipeline
   is mature enough to move into Nextflow DSL2, the per-cycle reversion chain needs
   to be wired in as native Nextflow processes. Deferred until the codebase is
   restructured.

5. **Production tier update.** The notes-3 tier (`ra_eff < 3.0 AND tJac >= 0.7
   AND intact`) predates the reversion pass. With reversion available, the criterion
   for sending a design to the wet lab should be:
   `reversion_verdict == 'pose_holds' OR (cycle == 0 AND intact AND ra_eff < 3.0)`
   The confidence-aware tier from notes 5 (`ra_eff < 5 AND tJac >= 0.7 AND
   ipSAE >= 0.5`) also still needs formalising as a separate RFDiffusion-specific
   tier in `compare_versions.py`.

6. **iRMSD / DockQ metric.** Carried from notes 3. `ra_eff` is still the primary
   filter and it's still a whole-effector Cα RMSD with no interface restriction.
   The 7B1I flexible-termini problem from notes 3 hasn't been addressed. iRMSD is
   the easiest drop-in; DockQ is the most standard.

7. **`boltz2_iterate_steering.py` is now 4175 lines (2868 actual code).** Natural
   split points: cycle-0/cycle-1 phase handlers, reversion subcommands, and
   aggregate/compute-final-metrics. Target ~3 files of 800–1000 lines each. Deferred
   until Nextflow migration restructures the codebase anyway.

8. **Silent pass-through of unvalidated contaminated designs.** If
   `plan-reversions` crashes or stages zero entries despite non-empty
   `contaminated.json`, `iterate-collect-finalize` currently treats those designs
   as non-contaminated and keeps them. They should be treated as drops. The correct
   fix is to compare `contaminated.json`'s `n_contaminated` against the number of
   entries in `reversion_results.json` at finalize time, and drop anything that
   never got a verdict. Filed but not fixed this session.

## Things to remember when working on this code

- **Always pass `--n-cycles 2` for reversion validation.** `--n-cycles 1` silently
  skips the entire kickoff chain.
- **Always pass `--receptor-fasta` for RFDiffusion designs.** Never let the plan
  stage extract the receptor sequence from the Cα-only ground-truth PDB.
- **The `build_reverted_sequence` sanity check is intentionally strict.** It will
  hard-fail if the cumulative mutation history is inconsistent with the FASTA. This
  is a feature. If you see "steered sequence mismatch at position X", the most likely
  cause is a wrong label being passed to `read_cumulative_mutations` — check that
  `make_pathway_label(parent_label, cycle, design_idx)` is using the child's cycle,
  not `cycle - 1`.
- **Cycle-0 plan.json does not have `experiment_root` or `cycle`.** Any new
  subcommand that reads plan.json must handle this gracefully. Use
  `plan.get("cycle", 0)` and derive experiment_root from `workdir.parent` when
  cycle is 0.
- **`derive_design_region.py` uses `--metrics-json`, not `--rfdiffusion-metrics-json`.**
- **`submit_boltz2_negative_steering.sh` uses `--workdir`, `--receptor`,
  `--effector`** — not `--experiment-root`, `--receptor-chain`, `--effector-chain`.
- **The protected set = design_region ∪ true_interface.** Neither alone is
  sufficient. design_region misses the native-scaffold anchor residues. true_interface
  alone would protect all interface positions including ones the steering is supposed
  to mutate.
- **`receptor.fasta` is overwritten destructively for `pose_holds` designs.**
  The pre-reversion steered sequence is preserved only in
  `cycle_N/pathway_<label>/reversions/<label>/receptor.fasta`. The `mutations.tsv`
  is intentionally left as the steered history to keep burnt-positions correct.
