# Negative Steering Diagnostic — Session Summary

## The problem

The negative-steering pipeline converged cleanly on 6G10 (best effector
RMSD ~1.8 Å, `true_jaccard` up to 0.88) but failed on 7QPX across 1309
designs spanning 7 modes and 5 cycles — a hard floor at ~12.5 Å
effector RMSD with max `true_jaccard` only 0.31. The two complexes
had nearly identical sequences (4 interior point edits plus terminal
trims), which framed the initial question as: *which residues, when
introduced into the working background, flip Boltz from the correct
pose to the 7QPX-style wrong pose?*

## What we built

### 1. Sequence-override mechanism in the core pipeline

Added two flags to `boltz2_negative_steering.py plan`:

- `--receptor-fasta PATH` — override the receptor sequence Boltz folds
  while still using the ground truth for contact filtering,
  true-site protection, and RMSD ranking.
- `--effector-fasta PATH` — same for the effector.
- `--skip-steering` — run only the initial cold-start prediction and
  stop; useful for cheap sequence-walk diagnostics.

The override loader (`_load_seq_override` in `cmd_plan`) validates
length against the ground-truth chain, reports all differing
positions, and persists the override paths in `plan.json` so the
iterate driver can pick them up. `rec_seq_truth` and `eff_seq_truth`
are kept separate from `rec_seq` / `eff_seq` so that
`find_contact_residues_heavy` on the ground truth uses the truth
sequence (satisfying its `expected_rec_seq` assertion) while the
call on the prediction uses the override (because that's what was
actually folded). Without overrides, behaviour is fully backwards
compatible.

### 2. Effector-override propagation in the iterate driver

`boltz2_iterate_steering.py cmd_iterate_plan` was silently
re-extracting the effector from the ground truth every cycle at
`eff_seq = bns.get_chain_sequence(ground_truth, truth_eff_chain)`.
Patched to read `effector_fasta_override` from `cycle_0/plan.json`
if present and fall back to the old behaviour otherwise. The
receptor side propagates by itself because each cycle's `parent_seq`
is read from the parent's `receptor.fasta`, which already carries
the override forward.

### 3. Singularity bind-path fix in the submission script

`submit_boltz2_negative_steering.sh` now intercepts
`--receptor-fasta` / `--effector-fasta` during arg parsing, validates
the files exist on the login node (fast-fail instead of after sbatch
has queued jobs), and adds the parent directories to the Singularity
bind list so the overrides are visible inside the container. The
flags still pass through the existing `EXTRA_PLAN_ARGS` mechanism to
reach `plan`.

### 4. Waypoint generators

- **`generate_waypoints.py`** — emits the original 8 7QPX-anchored +
  2 6G10-anchored chimera waypoints for the first diagnostic pass.
- **`generate_waypoints_7QPX_chainB.py`** — regenerated the
  7QPX-anchored waypoints at 73-aa chain-B length after we discovered
  the chain assignment was wrong.

Both produce per-waypoint FASTA pairs and a ready-to-run submission
shell script with narrow-steering settings
(`--mode mild --max-mutations 4 --candidate-pool-size 8 --n-designs 4 --n-cycles 2 --no-effector-template`),
about 8 Boltz calls per waypoint.

### 5. Aggregation and combine scripts

- **`aggregate_waypoints.sh`** — SLURM job that runs `aggregate` per
  waypoint inside the boltz container (which has gemmi, numpy,
  Bio.Align) and then attempts a combined pandas summary inside the
  TensorFlow container (which, despite the `.def` file, turned out
  not to have pandas).
- **`combine_waypoints.py`** — stdlib-only combiner that reads the
  per-waypoint `all_results_multicycle.csv` files directly, produces
  `runs/waypoint_combined_results.csv` (every design, tagged with
  `anchor` and `waypoint`) and `runs/waypoint_summary.csv` (one row
  per waypoint with initial / best / max statistics). Works on the
  login node with no container.

## The diagnostic walk and what it revealed

Ran 10 waypoints (8 anchored to 7QPX chain A, 2 chimera controls
anchored to 6G10) against the 7QPX chain-A receptor + chain-C
effector. Results:

| anchor | waypoint | best eff Å | max true_J | reads as |
|---|---|---|---|---|
| walk_6G10 | chimera_recEdits_only | 2.56 | 0.70 | rescued |
| walk_6G10 | chimera_effEdits_only | 3.01 | 0.74 | rescued |
| walk_7QPX | rev_recE72S | 14.08 | 0.21 | floor |
| walk_7QPX | rev_effN16H_D37A | 14.36 | 0.33 | floor |
| walk_7QPX | baseline_7QPX | 14.46 | 0.22 | floor |
| walk_7QPX | rev_recE72S_effN16H | 14.83 | 0.24 | floor |
| walk_7QPX | rev_effD37A | 14.89 | 0.27 | floor |
| walk_7QPX | rev_recE72S_effD37A | 16.27 | 0.23 | floor |
| walk_7QPX | rev_effN16H | 16.31 | 0.33 | floor |
| walk_7QPX | rev_all3 | 16.65 | 0.25 | floor |

Sharp gap — 3 Å vs 14 Å, nothing in between — with the gap aligning
perfectly with the *anchor*, not with which edits were applied.
`rev_all3` (every 7QPX→6G10 reversion expressible within the 7QPX
scaffold) sat at the same floor as `baseline_7QPX`. Every sequence
hypothesis was ruled out.

## The actual cause

You inspected 7QPX in ChimeraX and found that the interacting pair is
**chain B + chain C**, not chain A + chain C. The original pipeline
had been pointed at the wrong receptor chain, so
`find_contact_residues_heavy(7QPX.pdb, "A", "C", ...)` was reading
the "true site" off an unrelated chain pair and Boltz was being
measured against an interaction that doesn't exist in the relevant
biological assembly.

Chain comparison from the 7QPX PDB:

| chain | length | notes |
|---|---|---|
| A | 72 aa | receptor, not in contact with chain C |
| B | 73 aa | receptor, actual chain-C partner. More N-terminally deleted than A, has 3 extra C-terminal residues `AKE` |
| C | 81 aa | effector, unchanged |
| D, E | 70, 72 aa | additional receptor copies |
| F | 82 aa | additional effector copy |

The 6G10→chain-B edit set is different from the chain-A version: 3
receptor substitutions (S72E at chain-B pos 69, N75K at 72, K76E at
73) plus the 2 effector substitutions (H16N, D35A).

## Chain-B rerun (diagnostic)

Regenerated the 7QPX-anchored waypoints against chain B (73 aa) and
reran with narrow steering:

| anchor | waypoint | best eff Å | max true_J |
|---|---|---|---|
| walk_7QPX_chainB | rev_recE73K | **2.23** | 0.87 |
| walk_7QPX_chainB | baseline_7QPX | **2.41** | 0.67 |
| walk_7QPX_chainB | rev_effN14H | **2.78** | 0.83 |
| walk_7QPX_chainB | rev_effD35A | **2.92** | 0.91 |
| walk_7QPX_chainB | rev_rec_all | 4.50 | 0.63 |
| walk_7QPX_chainB | rev_all | 4.90 | 0.67 |
| walk_7QPX_chainB | rev_recE69S | 6.34 | 0.62 |
| walk_7QPX_chainB | rev_recK72N | 6.80 | 0.64 |

Every chain-B waypoint converged. `baseline_7QPX` went from 14.46 Å
(chain A) to 2.41 Å (chain B) with no sequence editing at all. The
sequence walk was never the right question; the chain assignment
was the only thing wrong. The elaborate hypotheses built around the
chain-A results (the `del_E12` deletion theory, the "Boltz
fundamentally can't predict 7QPX" hypothesis, the alanine-padding
diagnostic) can all be discarded.

## Full `mild_default` reruns with corrected chains

With the chain fix confirmed, reran the original `mild_default`
experiment (`--mode mild --max-mutations 6 --candidate-pool-size 6
--n-designs 10 --n-cycles 3`) on the three affected complexes:

| complex | chains | cold start eff RMSD | best after steering | true_jaccard | regime |
|---|---|---|---|---|---|
| 7B1I | B / C | 5.47 Å | n/a (no steering needed) | — | already correct |
| 7QPX | B / C | 24.6 Å | **1.04 Å** | ≈ | rescued in 1 cycle, no further headroom |
| 7QZD | B / C | 24.6 Å | **1.07 Å** | 0.92 | rescued + refined across 3 cycles |

**7B1I** converged on the cold start alone. Boltz placed the effector
at 5.47 Å with both chains folding correctly (rec 2.55 Å, eff 0.51
Å), and all 19 residues at the "wrong" interface were also on the
true interface — no surface to push against. The plan stage
correctly wrote `skip_steering=true` and exited cleanly. This is
success, not failure; negative steering is the wrong tool for
complexes Boltz already gets right.

**7QPX** converged in a single cycle. Cold start was 24.6 Å; the two
surviving cycle-0 designs landed at 1.04 Å and 1.35 Å effector RMSD.
Cycle 1 exhausted immediately because every remaining wrong-interface
residue now belonged to the true interface — same "nowhere left to
push" condition as 7B1I, but after steering rather than before.

**7QZD** showed the canonical negative-steering pattern: wildly
wrong cold start (24.6 Å), dramatic rescue by cycle 0 (2.06 Å),
continued refinement through cycle 1 (1.71 Å, `true_jaccard` 0.92,
22 of 24 true-interface residues correctly identified) and cycle 2
(1.07 Å). The single cleanest positive example across the three.

## Files produced

Code (drop into `~/receptor_design/negative_steering/`, overwriting
the first three):

- `boltz2_negative_steering.py` — override flags + skip-steering
- `boltz2_iterate_steering.py` — effector override propagation
- `submit_boltz2_negative_steering.sh` — bind-path + fast-fail
- `generate_waypoints.py` — initial 10-waypoint diagnostic (chain A)
- `generate_waypoints_7QPX_chainB.py` — 8-waypoint rerun (chain B)
- `combine_waypoints.py` — stdlib-only results combiner
- `aggregate_waypoints.sh` — SLURM aggregate driver

## Lessons for the benchmark at large

All three of the "failing" complexes (7QPX, 7B1I, 7QZD) turned out to
need chains B and C rather than A and C. That's not a coincidence —
it's consistent with a crystal-form or deposition convention where
chain A is a monomer or a crystallographic partner and the
biological interaction sits on later chains. A small audit script
that, for each benchmark entry, counts the number of heavy-atom
contacts between every chain pair and flags entries where the
user-specified pair has noticeably fewer contacts than the best pair
would catch this across your whole benchmark in a single pass.
Worth a few hours of work.

The override infrastructure, the waypoint generators, the iterate
driver fix, and the combine script all worked correctly. They were
exercising the right machinery against the wrong target. They are
useful infrastructure for any future experiment where you legitimately
need to walk a sequence between two complexes — worth keeping even
though the specific diagnostic they were built for turned out to be
unnecessary.
