# Experiment Notes 03 — Strategy 4e expansion and PWT7 v1_4e post-mortem

Third session built around two threads: (a) extending Strategy 4e
("Josh Intuition" — α1/β2-face binding hypothesis) to the remaining
campaigns, culminating in a PWT3-specific assembly using the AFDB
model in place of the in-house AF3 monomer; and (b) diagnosing the
completed PWT7 v1_4e run, which failed with 0/64 designs passing
the Rosetta SC filter and exposed a pipeline bug in
`NEGSTEER_WITHIN_SEQUENCE_PLOTS`.

## Headline

By session close:

- Strategy 4e runs exist for all five campaigns. PBY2 and PWT7 4e
  were committed between sessions; PWT3 4e was assembled and committed
  this session using a different structural strategy (AFDB model
  aligned to AVR-Pia rather than an in-house AF3 monomer).
- PWT7 v1_4e has completed and failed: 0/64 designs passed Rosetta SC
  (best SC = 0.370, threshold 0.5). The failure is informative — the
  α1/β2-face hypothesis for PWT7 may not be viable, or the starting
  pose is too far from a designable geometry.
- A pipeline bug identified: `NEGSTEER_WITHIN_SEQUENCE_PLOTS` fails
  with a missing-output error when no design sequences make it through
  the Rosetta gate. The script exits 0 and prints "Nothing to plot"
  but Nextflow expects `negsteer_*.png` output files regardless. Needs
  a fix before the next run with a potentially failing SC gate.
- ChimeraX syntax should always be looked up via WebFetch rather than
  recalled from memory. Two wrong commands were given before the docs
  were consulted.

## What landed between sessions (commits after d261310)

Two commits between the session 2 notes and the start of this session:

- **`7ce321c` — negsteer: treat skip_steering plan exit as soft
  failure, not hard abort.** Pipeline fix. A `skip_steering` exit from
  the negative-steering planner was being treated as a hard abort,
  killing the entire negsteer arm. Changed to a soft failure so the
  run completes with whatever sequences made it through, rather than
  dying entirely. Relevant to the PWT7 4e result — the run did
  complete and produce control data.

- **`0e2f663` — campaigns: strategy 4e for pby2 and pwt7.** Added
  `v1_4e/params.yml` for both `pikp1_pby2` and `pikp1_pwt7`, and
  assembled the corresponding `_4e_complex.pdb` input files. Both
  use the α1/β2-face contig `A1-13/22-36/A43-78 B` (matching
  `pikp1_avrpia v1_4b`). The PWT7 4e complex was assembled by
  double-anchoring onto 6Q76: Pikp AF3 superposed onto 6Q76/A, PWT7
  (9TFP chain A) superposed onto 6Q76/B (AVR-Pia position).

## Strategy 4e — the hypothesis and the contig

Strategy 4e is the hypothesis that a given effector is a MAX effector
and binds Pikp-HMA on the same α1/β2 face as AVR-Pia. The template
for both the complex assembly and the contig is the native
Pikp-HMA/AVR-Pia complex (6Q76, 1.90 Å, Varden 2019).

Assembly recipe for 4e: superpose the AF3 Pikp-HMA onto 6Q76/A
(receptor anchor), superpose the effector onto 6Q76/B (AVR-Pia
position, fold-only alignment), delete the template, combine. The
effector sits in the AVR-Pia pocket by construction; the quality
check is the alignment RMSD — < 3.5 Å over ≥ 40 Cα pairs supports
the MAX hypothesis, > 5 Å is a poor match.

Contig (all 4e runs): `A1-13/22-36/A43-78 B`

Design region: chain A residues 14–42 (α1/β2 face, 29 residues in
the native scaffold). Anchors: A1-13 (β1 + N-terminal tail) and
A43-78 (β3 through α2 + C-tail). The 22-36 range allows RFDiffusion
to use between 22 and 36 residues to span the interface — 7 below
and 7 above native length.

## PWT3 4e — assembling from the AFDB model

The existing PWT3 campaign (`v1_4b`) uses an in-house AF3 monomer
(`AF3/Pwt3.pdb`) that models only the mature chain (residues 19–141
of UniProt A0A223ZP76, numbered 1–123). For the 4e run, the user
provided `Pwt3_AFDB.pdb` downloaded directly from the AlphaFold EBI
database — the full precursor with UniProt numbering (signal peptide
~1–18, mature domain ~19–141).

### Why a different structural input for 4e

The 4e assembly strategy replaces the double-anchor-onto-6G10
approach used in `v1_4b` (which aligns PWT3 to AVR-PikD) with a
direct fold-only alignment of PWT3 onto AVR-Pia (chain B of the
assembled avrpia complex). This tests a specific variant of the MAX
hypothesis: that PWT3 binds the α1/β2 face like AVR-Pia rather than
the AVR-PikD face used in the v1_4b template.

The AFDB model was used rather than the in-house AF3 monomer because
it was the available structure at session start. The numbering
difference (full precursor vs mature-only) is noted in the v1_4e
params.yml header and does not affect the fold-only alignment or the
contig, which operates on chain A (Pikp-HMA) only.

### Assembly commands

The full ChimeraX workflow for the PWT3 4e complex:

```chimerax
# Open Pwt3_AFDB and the assembled avrpia complex
open .../experiments/inputs/structures/AF3/Pwt3_AFDB.pdb
open .../experiments/campaigns/pikp1_avrpia/inputs/pikp1_avrpia_complex.pdb

# Fold-only alignment of Pwt3 onto AVR-Pia (chain B)
matchmaker #1/A to #2/B alg sw ssFraction 1.0

# Drop AVR-Pia to avoid chain B collision in combine
delete #2/B

# Relabel Pwt3 A→B for contig convention
changechains #1/A B

# Combine: Pikp-1 first (stays chain A), Pwt3 second (chain B)
combine #2 #1 modelId 4 name pikp1_pwt3_4e_complex close true

# Save
save .../experiments/campaigns/pikp1_pwt3/inputs/pikp1_pwt3_4e_complex.pdb #4
```

The key insight for chain ordering: in ChimeraX's `combine`, the
output chain IDs are inherited from the source models in the order
listed. Listing Pikp-1 (`#2`, chain A) first and Pwt3 (`#1`, chain B,
relabelled) second produces chain A = Pikp-1, chain B = Pwt3 — the
convention every other campaign uses.

### Commits this session

- **`f2e171e`** — `Pwt3_AFDB.pdb` and `pikp1_pwt3/runs/v1_4e/params.yml`
- **`f7bac84`** — `pikp1_pwt3_4e_complex.pdb`

The complex PDB was committed separately because ChimeraX had not yet
been run when the params were committed; the save command was provided
for the user to run interactively, and the file was committed in the
following message.

## PWT7 v1_4e run — failure analysis

Run completed in 1h 27m. **Success: false.**

### Root cause: 0/64 designs passed Rosetta SC

All 64 RFDiffusion backbones failed the shape complementarity filter:

| Metric | Range |
|---|---|
| SC value | 0.204 – 0.370 (best: design_29) |
| packstat | **0.0 for every design** |
| dG_separated | +4,495 to +5,251 REU (all positive) |
| sc_threshold | 0.5 (none pass) |

Because nothing passed Rosetta, MPNN never ran, `NEGSTEER_RUN_ONE`
never ran, and the negsteer arm operated on controls only (polyA +
scrambled).

The SC scores (0.20–0.37) are uniformly poor. A well-packed
protein–protein interface typically scores > 0.55. Even the best
design at 0.370 would be considered marginal. The scores are not
clustered just below threshold — they are genuinely low, suggesting
the 4e complex geometry does not provide a surface that RFDiffusion
can design tightly against.

The uniformly positive dG_separated values reinforce this: the
backbone is not producing a favourable interaction with PWT7 in the
superposed pose. This is consistent with the α1/β2-face hypothesis
being wrong for PWT7 (PWT7 may not bind Pikp-HMA in the AVR-Pia
pose at all), or with the fold-only alignment placing PWT7 in a
non-interacting orientation.

**Packstat = 0.0 for all 64 designs** is suspicious enough to
warrant a separate check. This could indicate a genuine packing
failure or a Rosetta scoring issue with this complex. It does not
change the outcome (the SC scores independently fail all designs)
but should be investigated before drawing strong conclusions.

### Pipeline bug: NEGSTEER_WITHIN_SEQUENCE_PLOTS

When no design sequences reach the negsteer arm, the
`NEGSTEER_WITHIN_SEQUENCE_PLOTS` process runs with only controls and
prints `Nothing to plot` to stdout (exit status 0). However, Nextflow
declares the process as failed because it expects `negsteer_*.png`
output files that are never created. The process retried once and
failed again with the same result, crashing the run.

The fix: the script should either write a placeholder PNG in the
no-data case, or the Nextflow process definition should not declare
`negsteer_*.png` as a required output when the input summary contains
zero sequences with rep-sg data. This needs to be fixed in the
pipeline codebase before the next run where the SC gate may eliminate
all designs.

### Ranker miscalibration warning

The scrambled receptor control scored `ipSAE_min = 16.49`, well above
the `controls_warning_ipsae_max = 0.5` threshold. Boltz2 predicted
the scrambled control as a strong binder. This would be a serious
concern if any designs had made it through — it would suggest the
ranker cannot discriminate genuine binding from noise for this target.
Since no designs passed Rosetta it is moot for this run, but the
miscalibration flag should be carried forward if a future PWT7 run
produces survivors.

### Best design inspection

design_29 (best SC = 0.370) was inspected in ChimeraX alongside the
input complex. The designed backbone was nearly identical to the input
— RFDiffusion essentially recovered the native Pikp-HMA scaffold
rather than remodelling toward a tighter interface. This is
characteristic of weak effector surface signal: when the target
surface does not present a complementary geometry, diffusion collapses
back toward the prior (the native scaffold).

Contact check via `contacts #2/A:14-42 restrict #2/B` returned 139
contacts between the design region and PWT7, indicating the effector
is geometrically proximal. The question of whether the contacts are
productive (i.e., whether the effector side chains are oriented toward
the Pikp-HMA surface) was raised but not fully resolved — to be
discussed at lab meeting.

### Is more design region flexibility the answer?

The hypothesis that increasing the upper bound of the design region
(e.g. `22-36` → `22-50`) would help was considered and judged
unlikely to resolve the core issue. The SC failures are not borderline
(best 0.370); they indicate the posed complex doesn't provide a
surface the design region can pack against. More chain slack might
produce longer loops, but if the effector surface is unfavourable or
the pose is wrong, extra residues will not improve packing. The more
informative next step is visual inspection of the side-chain
orientations at the interface before deciding whether to rerun with
different parameters.

## ChimeraX command syntax — operating constraint updated

Three incorrect ChimeraX commands were given during this session
before the documentation was consulted:

1. `select sel extend true` — returned "Expected a keyword."
   `extend` is not a parameter of `select`.
2. `zone #2/A:14-42 range 4.5 models #2/B extend true` — did not
   work. `zone` in ChimeraX controls display of atomic detail around
   a focused residue; it is not a proximity selection tool and does
   not have the semantics assumed.
3. The correct command for selecting full residues within a distance
   (found via WebFetch of the ChimeraX docs):

```chimerax
select zone <ref-spec> <cutoff> <other-spec> residues true
```

The `residues true` parameter is what expands the selection to full
residues rather than only the specific atoms within range. For the
session's use case:

```chimerax
hide #2/B atoms
select zone #2/A:14-42 4.5 #2/B residues true
show sel atoms
style sel stick
```

**Operating constraint**: ChimeraX command syntax must be looked up
via WebFetch against the official docs before being given to the user.
Recall from training is unreliable for specific parameter names and
argument order. This applies especially to `select`, `zone`, `show`,
`changechains`, and `matchmaker` which have non-obvious argument
conventions. The session 2 notes already documented that `changechains`
syntax was verified via WebFetch; this pattern should be the default
for any ChimeraX command, not just the ones that have already failed.

## What's running where

| Campaign | Run | Status |
|---|---|---|
| `pikp1_avrpia` | `v1_4b` | Complete (results not reviewed this session) |
| `pikp1_avrpia` | `v1_4a` | Params ready, not submitted |
| `pikp1_avrpikf` | `v1_4a` | Complete (from session 1) |
| `pikp1_avrpikf` | `v2_4b` | Pending re-run (carry-over from session 1) |
| `pikp1_pby2` | `v1_4a` | Params ready, not submitted |
| `pikp1_pby2` | `v1_4b` | Params ready, not submitted |
| `pikp1_pby2` | `v1_4e` | Params ready, not submitted |
| `pikp1_pwt3` | `v1_4b` | Params ready, gated on fold investigation |
| `pikp1_pwt3` | `v1_4e` | **Params and complex committed this session; ready to submit** |
| `pikp1_pwt7` | `v1_4a` | Params ready, not submitted |
| `pikp1_pwt7` | `v1_4b` | Params ready, not submitted |
| `pikp1_pwt7` | `v1_4e` | **Complete — failed (0/64 Rosetta SC)** |

## Tasks for next session

1. **Lab meeting discussion — PWT7 4e.** 139 contacts between the
   design region and PWT7 exist in the 4e pose, but all designs
   recovered the native scaffold. Is the effector side-chain
   orientation at the interface the issue? Discuss with Daniel whether
   the α1/β2-face hypothesis for PWT7 is worth pursuing further or
   should be abandoned.

2. **Fix `NEGSTEER_WITHIN_SEQUENCE_PLOTS` pipeline bug.** The process
   must not declare `negsteer_*.png` as required output when there
   are no design sequences. Produces a spurious pipeline crash that
   obscures the real failure mode (Rosetta SC). Fix before next run.

3. **Investigate packstat = 0.0 for all PWT7 4e designs.** Verify
   whether this is a Rosetta scoring issue with this complex type or
   a genuine signal about the pose.

4. **Review `pikp1_avrpia v1_4b` results.** These are complete but
   were not inspected this session. This is the first native-complex
   4b run and should set the baseline for what a well-calibrated run
   looks like.

5. **PWT3 4e submission decision.** The complex is assembled and
   params are ready. Before submitting, verify the matchmaker RMSD
   for the PWT3 AFDB → AVR-Pia alignment to confirm the MAX-fold
   hypothesis holds (< 3.5 Å over ≥ 40 Cα pairs). If the RMSD is
   poor, reconsider whether to submit.

6. **Submit `pikp1_pby2 v1_4a` and `pikp1_pwt7 v1_4a`** (carried
   from session 2). These are the geometrically-derived runs for the
   cross-system campaigns.

7. **Resubmit `pikp1_avrpikf v2_4b`** (carried from sessions 1 and 2).

8. **Revise `pikp1_pby2` and `pikp1_pwt7` literature sweeps** (carried
   from session 2) — update §4 and §6 to reflect measured residue
   ranges and explicitly document the literature-insufficiency finding.
