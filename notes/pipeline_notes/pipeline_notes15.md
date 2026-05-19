# Pipeline Notes 15 — pikp1_avrpia pose_solved campaign: analysis, diagnostics, and next steps

## Context

Two runs — `pose_solved_3a` (13 de novo residues, 3Å contact cutoff) and `pose_solved_5a` (25 de
novo residues, 5Å contact cutoff) — were submitted based on a constraint-solved complex pose from
`tests/haddock/pose_solver/`. pose_solved_5a achieved several Tier-A candidates (best: Ra_eff 3.58
Å, iPTM 0.654) while pose_solved_3a found only one Tier-B candidate. This document records the
diagnostic discussion and the design of the next iteration (`pose_solved_5a_2`).

---

## Why negsteering struggled: five root causes

### 1. Structural monoculture — hotspot over-constrains RFDiffusion

The hotspot (`B22,B24,B31,B33`) gave RFDiffusion a strong directional contact signal. Every
independent run converged to the same anti-parallel β-strand geometry, confirmed by the clustering
heatmaps showing nearly all 64 designs within a few Å of each other. Since negsteering only
modifies sequence (not backbone), working on 64 structurally identical scaffolds is equivalent to
working on 1.

**Fix for 5a_2:** Remove the hotspot entirely. The pose solver already positioned the effector
correctly; RFDiffusion should use that geometry plus the fixed-region context to guide inpainting
without the additional rigidifying signal.

### 2. Input complex has too many clashes, compressing inpainting geometry

With 13 heavy-atom clashes (< 2 Å) in the most recent solved pose, the effector sits very close to
— or partially within — the volume RFDiffusion must inpaint. The fixed flanking backbone and the
nearby effector surface together constrain the inpainted backbone to a narrow geometric channel,
reinforcing the monoculture. Note: the clashing sidechains themselves are irrelevant (RFDiffusion
strips the design-region backbone), but the *position* of the effector (set by the clashes) is what
matters. Fewer clashes → effector slightly further away → more geometric freedom for inpainting.

**Fix:** Re-run pose_solver with higher `W_INTERP` (global mode), targeting < 5 clashes. Accept
looser pair constraint satisfaction as the cost. See section "Pose solver re-run (high W_INTERP)"
below.

### 3. All four pair constraint distances violated in input (4.9–5.3 Å)

The solved pose used as RFDiffusion input had all four pair contacts above the 4 Å threshold.
RFDiffusion therefore inpainted a backbone designed for a slightly offset geometry. The 5–10 Å
Ra_eff band in pose_solved_5a (24/104 sequences) represents designs where Boltz-2 found the
effector in the right region but with wrong contact specificity — confirmed by manual inspection of
a 7 Å Ra_eff structure that failed to reproduce the pair contacts.

**Fix:** Use a solved pose where the key pair contacts (A73-B31, A71-B33) are tightly satisfied
(< 4 Å), even at the cost of some other properties. The "best" backup pose (loss ~12,991, A71-B33
at 4.015 Å) may have been a better RFDiffusion input than the most recent pose.

### 4. Rosetta Sc threshold (0.5) filtered correct scaffolds before MPNN

pose_solved_3a filtered 29/64 designs (45%) at Sc ≥ 0.5, vs 13/64 (20%) for 5a. The Rosetta Sc is
computed on RFDiffusion-placed sidechains, not MPNN-optimised ones. Scaffolds with excellent
backbone geometry but unfavourable RFDiffusion rotamers fail this filter. Note: the comparison of
"design_47" between runs was incorrectly cited — designs numbered identically across runs are
entirely unrelated. The key finding is only that the pass rate difference (55% vs 80%) reduced the
effective pool entering negsteering for 3a.

**Fix for 5a_2:** Lower `sc_threshold` to 0.3 to allow more scaffolds through. Accept that some
lower-quality designs will enter MPNN; downstream filtering (negsteering, orthogonal metrics)
handles further quality control.

### 5. Negsteer pool depleted upstream, leaving insufficient structural diversity

After Rosetta filtering, the effective structural diversity entering negsteering is much lower than
the nominal sequence count suggests. With monoculture backbones and 29/64 filtered, the pool may
represent only 2–3 truly distinct structural families. The 97/104 `no_reversion` outcomes confirm
negsteering had almost nothing to work with. Note: the protected set MUST remain as
`design_region_union` — this is the entire point of negsteering, assessing the RFDiffusion-designed
interface unchanged.

**Fix for 5a_2:** Increase `num_designs` from 64 to 128 to increase the probability of sampling a
correct-geometry family. Combined with hotspot removal (Reason 1) and lower Sc threshold
(Reason 4), this should substantially increase the structural diversity reaching negsteering.

---

## Proposed experiment: standalone Boltz-2 on the input complex

**Idea:** Before submitting any new RFDiffusion run, run Boltz-2 directly on the native Pikp-1_HMA
+ AVR-Pia sequence pair (from the solved pose) to ask: can Boltz-2 recognise this interface at all?

**Rationale:** Since RFDiffusion reliably reproduces the same beta-strand geometry, and negsteering
evaluates Boltz-2 predictions of that geometry, a prerequisite question is whether Boltz-2 can find
the correct binding mode for the sequences involved. If Boltz-2 consistently achieves Ra_eff < 5 Å
on the native sequence pair with the solved-pose geometry as the "truth", that validates that the
geometry is credible and Boltz-2 should in principle be able to evaluate designs at this interface.
If Ra_eff > 15 Å even for the native, the geometry itself may not be recognisable to Boltz-2 —
which would explain all negsteering failures regardless of sequence quality.

**Caveat:** Pikp-1_HMA does not naturally engage AVR-Pia at this interface via the HMA domain in
the same way a designed binder would. Poor Ra_eff could reflect biology, not geometry. Still worth
running as a quick (minutes) sanity check before days of RFDiffusion compute.

**Implementation:** Run Boltz-2 directly on the input complex sequences (native, no RFDiffusion) in
"prediction" mode with AVR-Pia as the effector. Measure Ra_eff relative to the solved pose.
Multiple seeds (3–5) to assess consistency.

---

## New run: pose_solved_5a_2

Based on the above, `pose_solved_5a_2` makes the following changes vs `pose_solved_5a`:

| Parameter | pose_solved_5a | pose_solved_5a_2 | Reason |
|---|---|---|---|
| `hotspot` | `B22,B24,B31,B33` | (removed) | Allow backbone diversity |
| `sc_threshold` | 0.5 | 0.3 | Let more scaffolds through to MPNN |
| `num_designs` | 64 | 128 | More independent structural samples |
| Everything else | same | same | — |

Input complex, contigs, weights, negsteer params all unchanged.

---

## Pose solver re-run (high W_INTERP)

Re-running the pose solver with `--global-interp` and `W_INTERP` increased from 10 to 100, keeping
all other parameters the same (4 constraints: A73-B31, A71-B33, A8-B22, A7-B24; exclusion
A8-B33@4.5; VALID_TOL=-1.0; 1000 restarts). Goal: reduce clashes to < 5 while accepting looser
pair constraint satisfaction. Results reported below once run completes.

---

## Files changed

- `experiments/campaigns/pikp1_avrpia/runs/pose_solved_5a_2/params.yml` — new run
- `notes/pipeline_notes/pipeline_notes15.md` — this document
