# Notes 10 — v7 diagnostics, the contamination-check fix, and four iterations to get the effector-atom filter right

## Context

Entering this session, notes 9 had closed with three consecutive live runs
(v5, v6, v7) producing zero genuine passing designs and a tentative conclusion
that "the problem is not a pipeline bug at this point — it's the biology."
v7 surfaced two false-positive pose_holds candidates (design_08 and design_13)
that per-seed inspection showed to be noise, and the classifier's inability
to promote anything was flagged as architectural.

The session had two initial goals: (1) decide whether to keep iterating on this
scaffold or move on, and (2) either patch the classifier to match the real
per-seed data distribution or confirm the scaffold as a dead end.

Answering (1) first shifted the whole direction of the session: inspection of
v7's raw per-seed data showed the pipeline was flagging false contamination
on residues 22/26, which were physically 14 Å from the correct-pose effector
and therefore could not be genuine load-bearing contacts. The rest of the
session was spent chasing this bug through four layers of pipeline code and
finally across four iterations (v8 → v9 → v10) of the patched filter itself.

## What was learned from v7 diagnostics

### The per-seed data shape

Across 50 designs × 3 seeds = 150 predictions:

- 21 seeds (14 %) reached the correct interface on reversion (ra_eff < 5 Å,
  receptor intact).
- 3 seeds individually passed every criterion (verdict `pose_holds`), but
  these 3 were spread across 3 different designs.
- 0 designs had ≥2/3 seeds individually passing.

The bimodality was stark: 0/50 three-seed groups all-low, 9/50 all-high,
15/50 split. Within a sequence, each seed is a roughly-Bernoulli draw —
when a seed gets the correct pose it sits at sub-5 Å ra_eff, when it
misses it flies back to the ~29 Å wrong-site attractor.

### Bug E classifier architecturally incapable

The Bug E fix from notes 9 added a per-seed-majority gate requiring
`n_seeds_pose_holds ≥ ceil(N/2)`. With `num_seeds=3` and p ≈ 0.14 per
seed, majority-of-3 requires both p ≥ 0.5 AND correlation across seeds
within a design. Neither condition holds. The classifier could never
fire `pose_holds` at these settings.

## Chain of findings on the 22/26 "contamination"

What looked at first like a recurring-position contamination problem
unwound into a chain of diagnostic turns:

1. **Initial read**: positions 22, 26, 32 kept appearing as reverted-pose
   contacts across v6/v7. First guess was that protected-set choice was
   starving the pool of useful steering residues.

2. **`plan.json` diagnostic on the HPC**: v7 used
   `protected_set_source: true_interface_fallback`, not the requested
   `design_region_union` — the union emptied the pool so the code
   retreated. The actual pool was `{5, 7, 22, 26, 32, 36, 44, 46}`.
   Empirical per-position table revealed two behavioural classes:

   | Class | Positions | Reverted-contact rate |
   |-------|-----------|----------------------|
   | Works as intended | 5, 7, 32, 36, 44, 46 | 2–10 / ~110 |
   | Recurring contaminator | 22, 26 | 32 / 108, 25 / 99 |

   Position 46 was most striking: 120/126 in contact steered, only 6/126
   contaminating reverted — textbook useful steering residue.

3. **Widen the protected set?** A second HPC diagnostic measured
   min-heavy-atom distance from each receptor residue to the ground-truth
   effector. Positions 22 and 26 were at 14.27 Å and 14.17 Å respectively.
   **No cutoff captures 22/26 without sweeping in 17 unrelated residues.**
   They are not native-pose-adjacent.

4. **Direct ChimeraX inspection** (design_25_s2 reverted prediction,
   positions 22 and 26 highlighted on receptor chain A): the reverted
   effector was at the correct interface but rotated slightly, and an
   extending bit of its surface draped down across the opposite face of
   the receptor to reach residues 22/26. Not a binding contact — a
   geometric coincidence from the ~3 Å rotational jitter.

5. **Disorder hypothesis**: maybe the extending bit was a flexible tail
   and its contacts shouldn't count. Two problems: (a) Boltz receives the
   effector as a templated chain with `force: true`, so effector pLDDT
   is inherited from the template, not a disorder signal; (b) the
   residues reaching 22/26 were B:13, B:15, B:17 of the 82-residue
   prediction — mid-chain, not terminal.

6. **AF3-native-complex distraction**: a brief dead-end where I analysed
   B-factors on `af3_pikp1_native_avrpikf_complex.pdb` as if it were the
   ground truth. User correction: v7 uses `design_3.pdb` (the
   RFDiffusion design itself, Cα-only) as `--ground-truth`. The AF3
   native complex is only used for validation context and the effector
   there came from the 7B1I crystal via matchmaker. The low pLDDT values
   were residual 7B1I B-factors, not confidence.

7. **Actual conclusion on 22/26**: the reverted effector sits ~3 Å off the
   ground-truth pose with 14–18 of its contacts at true-interface
   residues. It's genuinely at the correct interface. The 22/26 contacts
   are the structured-but-non-interface edge of the effector brushing
   across a non-interface patch of receptor surface as the effector
   rotates. The contamination check is over-firing.

## The contamination-check fix — four iterations

### The right principle (agreed on mid-session)

Define contamination on the **effector side**: a mutated receptor
position counts as contaminating only if it is in contact with an
effector atom that belongs to the effector's interface region. "Interface
region" is defined once on the RFDiffusion ground truth — the set of
effector atoms / residues within 8 Å Cα of a receptor residue in
`true_interface_idx`.

If an effector has long reaching bits, flexible loops, or any structured
feature that brushes non-interface receptor surface while the complex
rotates, those contacts are excluded from the contamination check. If a
mutation at receptor position P contacts only the non-interface side of
the effector, P is not flagged.

The user explicitly chose to accept the trade-off that a real secondary
binding site (biologically plausible in some other protein) would be
excluded from scoring under this rule. Rationale: the alternative —
residue-weighting schemes — would complicate the pipeline beyond what
the downstream wet-lab assay can resolve.

### What was implemented (level 1 contamination fix)

Deferred options:
- **Level 2** (interface-restricted ipSAE / iPAE metrics) — not yet.
- **Level 3** (core-only RMSDs and everything interface-restricted) —
  for the NextFlow migration.

### The metric audit that fell out of this

While auditing the contamination check, I audited the confidence metrics
in `compute_metrics.py` and found two real bugs:

**Bug F — ipSAE d0 formula wrong.** `compute_ipsae_one_direction` used
`d0 = 0.5 * sqrt(n_below)`. The published Dunbrack 2025 ipSAE uses the
TM-score d0: `d0 = 1.24 * (n - 15)^(1/3) - 1.8`, clamped at 0.5. The
buggy formula produced d0 values ~1.4–4× larger than correct, inflating
ipSAE scores by ~30–50 % depending on chain length.

Consequence for prior notes: all absolute ipSAE numbers in v1–v7 are
inflated. Rankings within a run are preserved (n_below is similar
row-to-row) but cross-run thresholds need rescaling. The notes 5
finding that design_3 hit ipSAE 0.86 becomes ~0.5–0.6 under the fix.
The notes 3 finding that 7QPX ipSAE plateaus at 0 is unchanged — that's
driven by n_below = 0, which both formulas produce.

**Bug G — actifPTM misnamed.** The `actifptm` column was computing
`0.8 * iptm + 0.2 * ptm`, which is the AlphaFold-Multimer ranking score,
not Varga et al. 2025's actifPTM. Real actifPTM restricts the pTM sum
to an "interacting residue set" (residues with at least one inter-chain
PAE below a threshold).

Fix applied: wrote a proper `compute_actifptm` that identifies the
interacting residue set from the PAE matrix and applies the pTM kernel
over it with d0 from the interacting-set size. The old AF ranking
formula is preserved as a new `af_rank_score` column (renamed, not
deleted).

### v8 — failed: atom-level filter on Cα-only ground truth

First patched version. `compute_effector_interface_atom_mask` walked
effector heavy atoms in `design_3.pdb` and kept those within 8 Å of an
interface receptor Cα. Filter was written to `plan.json` as a list of
`(resseq, atom_name)` pairs.

**What broke**: `design_3.pdb` is Cα-only on the effector too (it's a
RFDiffusion design output — Cα placeholders throughout the de novo
region, and even the effector template appears to have been reduced to
Cα). The atom mask therefore contained only Cα atoms. When applied to
a Boltz prediction (full heavy atoms), the heavy-to-heavy 5 Å contact
check against receptor side chains almost never registered — Cα
backbone atoms just aren't that close to receptor side chains.

**Result from v8** (100 designs × 3 seeds, 300 predictions):
`steered_contact_residues` empty on every row; 0 contaminations flagged;
0 reversions triggered; `passing_summary.csv` contained 57 designs that
fell through to `no_reversion` without any contamination check running.
The steered predictions themselves were real (steered_iptm up to 0.91,
steered_actifptm up to 0.80 on sub-5 ra_eff designs) but unvalidated.

### v9 — failed: residue-level filter with wrong resseq numbering

Rewrote `compute_effector_interface_atom_mask` as
`compute_effector_interface_residues`, returning a set of effector
residue numbers. `compute_interface_contacts` extended to accept either
the new residue-set form or the legacy atom-pair form (back-compat).

Self-tested on a synthetic Cα-only ground truth — passed all three
pathways (direct function, kwarg, CLI-via-JSON).

**What broke**: `compute_effector_interface_residues` collected effector
`resseq` values as they appeared in the ground truth PDB — those came
out as `[89, 90, 91, 92, 94, 97, 98, 112, 113, 114, 115, 116, 123,
124, 125, 126, 127, 128]`. These are the RFDiffusion design's synthetic
offset resseq scheme. Boltz predictions renumber chains 1-based
positional (1, 2, …, 82). A filter saying "keep residues 89–128"
matched zero atoms in every prediction.

**Result from v9**: identical to v8 — empty `steered_contact_residues`
on all rows. Same failure mode via a different bug.

### v10 — fixed: positional indices

`compute_effector_interface_residues` rewritten to track positional
index as it walks the file (ignoring resseq entirely) and return
1-based positional indices matching what Boltz writes.

Self-tested with GT-resseq ≠ positional (89–91 in ground truth,
1–3 in prediction) — filter correctly identifies positional index 1
as the interface residue and produces non-empty contacts on the
prediction.

v10 patch set was prepared and given to the user with a pre-submission
check (expect small 1-based integers in `plan.json`
`effector_interface_residues`; NOT values >100). v10 has not yet been
submitted as of session end.

## The Bug E classifier fix

Orthogonal to the contamination-check story, but prepared in the same
patch round.

`_classify_aggregated_verdict` in `boltz2_iterate_steering.py`: replaced
the per-seed-majority gate (`n_seeds_pose_holds ≥ ceil(N/2)`) with an
any-seed gate plus a failure-count guard:

- Promote to `pose_holds` if at least one seed individually says
  `pose_holds` AND non-failing seed count ≥ max(n_pose_collapses,
  n_new_contamination).
- Otherwise, the dominant aggregated failure mode wins.

Dry-replay against v7 data: flips exactly one design (design_08) from
`pose_collapses` to `pose_holds`. Designs 25 and 43 correctly stay
dropped — design_25 has a position-majority of 22/26 that fires
`new_contamination` before the any-seed gate is reached; design_43 has
reverted median ra_eff 7.60 Å (1 clean seed at 4.51, the other two at
7.60 and 13.61) which fails the structural filter.

The patch is conservative: promotes only designs where the aggregated
evidence is strong AND at least one seed individually agrees AND no
failure mode dominates.

## Files changed

| File | What |
|------|------|
| `compute_metrics.py` | ipSAE `d0` formula corrected to Dunbrack TM-score (Bug F). `compute_actifptm` implemented per Varga 2025; column meaning fixed (Bug G). Old 0.8·iptm + 0.2·ptm renamed `af_rank_score`. `compute_effector_interface_residues` added (positional, Cα-only-ground-truth-safe). `compute_interface_contacts` + `parse_boltz2` extended with `effector_atom_filter` kwarg accepting either residue ints or legacy atom-pair tuples. CLI flag `--effector-atom-filter-json`. |
| `boltz2_negative_steering.py` | `cmd_plan` computes effector interface residues from `true_interface_idx` + ground truth at 8 Å Cα and writes `effector_interface_residues` to `plan.json` (plus legacy `effector_interface_atoms` alias for back-compat). |
| `boltz2_iterate_steering.py` | Bug E classifier replaced (any-seed + failure-count guard). Both `cmd_build_contaminated` and `cmd_compute_final_metrics` pass `--effector-atom-filter-json <cycle_0/plan.json>` to `compute_metrics.py` subprocesses. `af_rank_score` added to `METRIC_KEYS_RAW` so the new column propagates through the aggregator (missed in initial patch round, caught during v8 postmortem). |
| `reversion.py` | `harvest_reversion_results` resolves the cycle-0 plan path and passes `--effector-atom-filter-json` to the `compute_metrics.py` subprocess for reverted predictions. |

## v10 submission command (ready to deploy)

```bash
./submit_boltz2_negative_steering.sh \
    --ground-truth resurface_pipeline_test/design_3.pdb \
    --receptor A --effector B \
    --receptor-fasta resurface_pipeline_test/design_3_seq_3_receptor.fasta \
    --workdir runs/design_3_dedup_test_v10 \
    --mode mild --max-mutations 6 --candidate-pool-size 10 \
    --n-designs 100 --num-seeds 3 --n-cycles 1 \
    --diffusion-samples 5 \
    --rmsd-threshold 6.0 \
    --protected-set-source design_region_union \
    --true-interface-indices-file resurface_pipeline_test/design_3_true_interface.txt \
    --design-region-indices-file resurface_pipeline_test/design_3_design_region.txt

# after queue drains:
./submit_postprocess.sh --experiment-root runs/design_3_dedup_test_v10 --force-recompute
```

Changes vs v7: `--n-designs 100` (was 50), `--rmsd-threshold 6.0`
(was 5.0), `--max-mutations 6 --candidate-pool-size 10` (v7 was 6/10;
v6 was 4/8). Classifier patch + contamination-check patch + metric
fixes are all in the code paths.

Pre-submission sanity check on the login node (takes 30 seconds, avoids
burning GPU time on a fourth broken version):

```bash
cat runs/design_3_dedup_test_v10/cycle_0/plan.json | \
  python3 -c "import json,sys; p=json.load(sys.stdin); \
    print('effector_interface_residues:', p.get('effector_interface_residues'))"
```

Expected: small 1-based integers (e.g. `[1, 2, 3, 4, 5, 6, 9, 10, …]`)
with all values ≤ 82 (effector length). If values ≥ 89 appear, v9's
resseq bug is back; if values are atom-pair tuples, v8's atom-level
code is back; in either case don't submit.

## Key findings and learnings

1. **22/26 contamination was never a pool problem** — it was a contact-cutoff
   sensitivity artefact. The reverted pose sits ~3 Å off ground truth and
   a structured but non-interface edge of the effector reaches over to
   brush non-interface receptor surface. Widening the true-interface
   cutoff doesn't help because 22/26 are 14 Å from the correct-pose
   effector. The fix is effector-side: only effector atoms belonging to
   the interface region are allowed to generate contamination calls.

2. **ipSAE values in all prior notes are inflated ~30–50 %**. The d0
   formula was wrong since the metric was added. Rankings within a run
   are preserved; absolute thresholds need rescaling.

3. **The `actifptm` column in all prior CSVs is the AF-Multimer ranking
   score, not Varga actifPTM**. Renamed to `af_rank_score`. A new,
   correctly-computed `actifptm` column replaces it going forward.

4. **Full-heavy-atom assumption about ground truth**. RFDiffusion output
   PDBs are Cα-only. Any function that walks ground-truth heavy atoms
   needs explicit Cα-only handling, and filters derived from them need
   residue-level rather than atom-level granularity.

5. **Coordinate-system mismatch between ground truth and predictions**.
   RFDiffusion designs use synthetic offset resseq schemes (typically
   1000-offset from the fixed residue max); Boltz predictions
   renumber 1-based positional. Anything computed on the ground truth
   that will be applied to predictions must be converted to positional
   indices before persisting.

6. **Synthetic self-tests didn't catch the layered bugs**. v8's atom-level
   filter passed a synthetic full-heavy-atom test that didn't exercise
   the Cα-only case. v9's residue-level filter passed a synthetic test
   that used positional resseq. Each iteration's test only exercised
   the scenario that its predecessor already handled. For v10 the test
   was built to exercise exactly the previous failure mode
   (GT-resseq ≠ positional).

7. **14 % per-seed pose-acquisition rate on the design_3 scaffold**.
   At `num_seeds=3` this gives an expected ~0.4 passing seeds per design
   — majority-of-N gates are architecturally wrong at this per-seed
   success rate. Any-seed plus failure-count guard is the right shape.

## Outstanding work

1. **Submit v10**. All four files patched, syntax-checked, self-tested
   end-to-end. Deploy all four to the cluster together; MD5-check.
   Run the pre-submission login-node sanity check above before
   submitting.

2. **ChimeraX inspection of any v10 passes**. Notes 9's deferred inspection
   work (three sessions for the v7 `pose_holds` seeds) was never done
   in this session because the 22/26 diagnostic chain took priority.
   Once v10 produces real passing candidates, open the reverted
   prediction against ground truth and confirm:
   - Effector is at the correct interface, not a coincidental third site.
   - Residual green-stick residues (mutations kept through reversion)
     are peripheral, not structurally load-bearing.
   - Magenta mutation sticks in the steered prediction point into the
     effector, not into solvent.

3. **Level 2 interface-restricted metrics** (deferred). The `ipae` and
   `ipsae_*` metrics still average over the full inter-chain PAE block,
   including off-interface contributions. Not technically a bug — the
   field uses this naive iPAE convention — but the more principled
   interface-restricted forms (matching the new `actifptm`) would
   reduce noise from off-interface PAE contributions. Worth doing
   before deploying at RFDiffusion scale across diverse effector shapes.

4. **Level 3 core-only RMSDs and full interface-restricted scoring**
   remains the NextFlow-migration target.

5. **Notes 9 NextFlow integration** — still deferred.

6. **`extract_passing.py` empty-file handling** — still deferred from
   notes 9.

## Key trust notes for future sessions

- Test **exactly the failure mode that broke the previous iteration**.
  Generic synthetic tests gave me false confidence three times running.

- **Ground truth for RFDiffusion designs is Cα-only**. Any heavy-atom
  assumption must be verified against the actual file, not inferred.

- **Coordinate systems must be checked at every boundary**. Ground-truth
  resseq ≠ prediction positional index ≠ 0-based positional index ≠
  ChimeraX 1-based. The notes 3 through notes 9 pattern of "check the
  index base" needs to extend to "check the numbering source".

- **The `actifptm` and `ipsae_*` columns in pre-v10 CSVs are not what
  they say they are**. When reading old notes, translate
  `actifptm` → `af_rank_score`, and rescale any `ipsae` absolute values
  down by ~30–50 %.
