# HADDOCK Protocol Experiments (GOHREP)

A running record of all ideas tried and proposed for the HADDOCK module,
formatted as **Goal / Hypothesis / Rationale / Experimental Plan** (+
Results/Conclusion for completed work).  Companion to
`notes/design_audit.md` Session 7 (which captured the architectural
decisions) and `notes/remediation_state.md` (which tracks branch state).
Last updated: 2026-05-18.

---

## Background: HADDOCK's role in this pipeline

HADDOCK is used to place two monomer PDBs (receptor HMA domain + effector
target) in a compact, biologically plausible relative orientation to serve
as input to RFDiffusion.  Its job is **geometric placement** — not affinity
prediction, not energetic optimisation.

The receptor and target both come from crystal structures.  They are
biologically correct.  Neither should be structurally modified.  Sidechain
rotamer optimisation in HADDOCK is irrelevant because RFDiffusion +
ProteinMPNN replace the design-region sequences, Boltz-2 reproduces the
structure, and Rosetta FastRelax re-relaxes for the orthogonal metrics —
all downstream.

---

## Standing clarifications

These represent decisions or observations locked in over the course of this
work.  Do not relitigate without new evidence.

1. **Rigid-body equivalence**: In pure rigid-body docking, "freeze the
   receptor and move the effector" is equivalent to "let both move."  The
   relative configuration space is identical.  The only real difference is
   refinement scope (flexref/emref on both vs one body).

2. **HADDOCK never produces output clashes**: Empirically verified across
   multiple production runs.  flexref/emref refines clashes away.  The
   absence of clashes in HADDOCK output is not a constraint violation — it
   is the pipeline working as designed.

3. **Both proteins are rigid throughout**: Sidechain and backbone
   refinement in HADDOCK is wasted computation because all design-region
   residues are replaced downstream.  Fixed-region sidechains in their
   crystal rotamers are correct by definition.

4. **Hint pose quality is a manual responsibility**: Without HADDOCK
   performing the initial placement (which it cannot do in the presence of
   clashes from regions about to be redesigned), the only starting point is
   a ChimeraX-aligned complex.  The user accepts this.

5. **HADDOCK as placement optimiser, not energy minimiser**: The right
   framing is "find the best effector position given the anchored receptor
   and the specified restraints," not "minimise HADDOCK score."

---

## Literature context

### The core analogue: partial diffusion in RFdiffusion

The closest published precedent for pose-hint → adaptive refinement is
**RFdiffusion partial diffusion** (Watson et al. 2023, *Nature*; Glögl
et al. 2024, *Science*).  Partial diffusion takes a roughly-placed
structure, adds a controlled amount of noise, and lets RFdiffusion denoise
back — sampling near the hint rather than from pure random noise.  The
Glögl et al. TNFR paper demonstrated that partial diffusion for affinity
maturation improved hit rates from 7–35% to 30–46% for challenging targets
including GPCRs.  Picomolar affinity was achieved via a partial diffusion
step applied after initial nanomolar designs.

**Key difference from this pipeline**: In RFdiffusion partial diffusion,
the design region (motif vs designed) is still pre-specified by the user.
The partial diffusion step adjusts backbone geometry, not the boundary
between "keep" and "redesign."  The automatic derivation of the design
region from clashes is not a feature of this workflow.

### Antibody design as the closest paradigm

**RFantibody** (Bennett, Watson et al. 2025, *Nature*, "Atomically accurate
de novo design of antibodies with RFdiffusion") is the most directly
analogous published workflow.  The antibody framework (analogous to the HMA
domain) is provided as a template encoding pairwise geometry but NOT the
rigid-body relationship to the target — so RFdiffusion samples diverse
docking orientations.  CDRs (analogous to the design regions) are the only
parts designed de novo.  ProteinMPNN designs CDR sequences; the framework
sequence is frozen.

**The parallel to this pipeline:**
- Antibody framework → Pikp-1 HMA fixed regions (A1-32 + A50-68)
- CDR loops → HMA design regions (33-49 + 69-78)
- Target antigen → AvrPia effector
- Hotspot residues on target → B20-26 hotspots

**Key difference**: RFantibody doesn't start from a user-specified docked
pose or use HADDOCK.  It samples the dock directly via RFdiffusion.  The
design region boundaries come from antibody anatomy (CDR annotations) not
from clashes at a user pose.

### BindCraft: One-shot from hotspots

**BindCraft** (Pacesa et al. 2025, *Nature*, "One-shot design of functional
protein binders") is a fully automated pipeline: specify target hotspots
→ AF2 backpropagation generates binder backbone → ProteinMPNN designs
sequences.  If hotspots are not specified, AF2 selects the binding site
automatically.

**Relevance**: BindCraft can in principle design a novel binding loop
against a specified target face — which is what this pipeline attempts via
RFDiffusion.  The pipeline described here differs in that it engineers an
EXISTING receptor rather than designing a binder de novo.

### PPDiff: joint sequence-structure design of binder complexes

**PPDiff** (Song et al. ICML 2025) jointly diffuses the sequence and
backbone of protein-protein complexes, achieving 50% success on pretraining
tasks and 16-23% on antibody/mini-binder applications.  Does not take a
user-specified docked pose; samples from the complex distribution
conditioned on target structure.

### ProtFill: interface redesign via inpainting

**ProtFill** (Kozlova et al. 2023/24, NeurIPS MLSB) inpaints protein
structure and sequence simultaneously using SE(3) equivariant diffusion.
Demonstrated for antibody interface redesign.  Requires the design region
to be pre-specified; does not derive it from clash analysis.

### The Rosetta era: AnchoredDesign and loop remodeling

**AnchoredDesign** (Fleishman et al. 2011) grafts a key interaction from a
natural binder onto a scaffold surface loop, then redesigns the surface
with backbone flexibility.  Represents the pre-deep-learning version of the
user's workflow: scaffold (analogous to HMA fixed regions) + designed loop
(analogous to design regions) + fixed target.  Still requires the user to
specify which loops to redesign.

### The gap: clash-driven automatic design region selection

No published method identifies, from a user-provided docked pose, which
receptor residues are clashing and automatically designates those as the
redesign region.  This specific formulation — clash detection at a manually
placed complex → adaptive design region boundary → downstream generative
redesign of only the carved residues — appears to be methodologically
novel.

The closest are:
- Antibody humanisation tools that flag CDR residues clashing with a human
  framework as candidates for back-mutation (but use CDR anatomy, not
  structural clashes, as the primary selector).
- Rosetta RosettaRemodel protocols that remodel loops to resolve clashes
  (but require the user to specify which loops).

---

## Implemented experiments

---

### E1 — Original HADDOCK module restructure

**Status**: Implemented, in production.
**Commits**: 1d0d542 (Phase 4 types), a40f7cf (Nextflow rewiring), 71abda3
(EXTRACT_HOTSPOTS deleted), 0e8a68a (contig-derived design region), 5fc900c
(numpy fix), eb4ec37 (plot polish), 375b369 (sidechain strip + DR1 active).

**Goal**: Replace the pre-Session-7 HADDOCK module, which derived AIR
active residues from the contig's de novo regions (phantom residues that
don't exist on the input PDB at HADDOCK time), with a user-driven restraint
system and a BSA-primary cluster ranking.

**Hypothesis**: User-specified pair contacts + effector hotspot AIRs +
BSA-driven cluster ranking + manual checkpoint will produce biologically
better docked poses than the original contig-derived AIRs and HADDOCK-score
ranking.

**Rationale**: The original design derived AIRs from contig de novo
segments, but those residues don't exist on the native receptor PDB at
HADDOCK time.  AIRs referenced phantom residues; the auto-pick used
HADDOCK score (misaligned with the geometric-placement goal); no mechanism
existed for the user to inspect poses before committing to RFDiffusion.

**Experimental plan**:
1. New user-facing restraint params: `haddock_contact_pairs` (hard CA-CA
   pins), `haddock_receptor_active_residues` and
   `haddock_effector_active_residues` (soft AIR sides).
2. Cluster ranking by pair contact fraction (primary) + BSA (secondary).
3. `stop_after_haddock` gate + `haddock_chosen_cluster` resume.
4. HADDOCK3_PREPARE relabels input PDBs to chain A/B.
5. EXTRACT_HOTSPOTS deleted: auto-derived hotspots biased toward the
   existing docked interface and reduced design diversity.
6. 8 new diagnostic plots: score-vs-BSA (coloured by cluster), cluster
   overview heatmap, cluster ranking scatter, pair satisfaction, hotspot
   placement, clash breakdown, cluster sizes, interface heatmap.

**Test config (short run)**:
- `haddock_sampling = 100`, `haddock_seletop = 20`,
  `haddock_min_cluster_size = 2`
- `haddock_contact_pairs = "A73-B31 A72-B32 A71-B33"`
- `haddock_effector_active_residues = "20-26"`

**Results**: Pipeline runs end-to-end.  Pair pins satisfied (3/3).  Hotspot
B20-26 markers all 8-15 Å from receptor — none in contact range (<5 Å).
Effector beta strand bent down toward B20-26 instead of receptor rotating.

**Conclusion**: Functional but does not recover the intended pose.  The
effector deforms rather than the receptor rotating.  Restraint-following
works; pose recovery for the intended interface does not.

---

### E2 — Sidechain stripping of design region

**Status**: Implemented, in production.
**Commits**: 375b369 (initial), 0e8a68a (contig-derived design region).

**Goal**: Remove the receptor's design-region sidechain steric profile
during HADDOCK so the effector can approach the receptor backbone without
being deflected by sidechains that will be replaced by RFDiffusion anyway.

**Hypothesis**: The design-region sidechains are accidentally constraining
HADDOCK's pose search by acting as immovable steric walls, causing the
effector to bend rather than the receptor to rotate.

**Rationale**: HADDOCK never produces output clashes — it resolves them by
deforming the more flexible body (typically the smaller/softer one).  The
receptor's design-region sidechains prevent the effector from approaching
the design face, so the effector's own interface region bends down to meet
the effector hotspot from below.

**Experimental plan**:
1. `haddock_strip_design_sidechains: true`.
2. All atoms in contig de novo residues (33-49 + 69-78) stripped to N/CA/C/O
   and renamed GLY in `receptor_haddock.pdb`.  HADDOCK-only — RFD replaces
   these residues entirely.
3. Run at production-length (`haddock_sampling = 10000`).

**Results**: Correct — 27 residues stripped.  Run completes.  Hotspot B20-26
markers still >5 Å from receptor.  One cluster rotated the receptor in the
OPPOSITE direction, a different loop (not DR1) lying parallel to the hotspot
strand.  Sidechain stripping alone insufficient.

**Conclusion**: Necessary but not sufficient.  The design-region BACKBONE
(not just sidechains) is also constraining the docked pose.  The correct
form of the intervention is backbone removal.

---

### E3 — DR1 receptor active residues

**Status**: Implemented, in production.
**Commits**: 375b369 (alongside sidechain stripping).

**Goal**: Provide a soft directional pull forcing design region 1 (residues
33-49) of the receptor toward the effector hotspot patch B20-26, in addition
to the hard pair pins that anchor design region 2 (69-78) at the existing
AvrPikF-like interface.

**Hypothesis**: Adding DR1 to the receptor active list creates a two-sided
AIR ("each DR1 residue must contact one of B20-26") that, combined with the
pair pins on DR2, geometrically requires the receptor to rotate so both
design regions face the effector.

**Rationale**: The user confirmed DR1 faces the effector in their ChimeraX
placement.  A two-sided AIR makes the rotation a soft requirement.

**Experimental plan**:
1. `haddock_receptor_active_residues = "33-49"`.
2. Combined with sidechain stripping (E2).
3. Production-length run (10000/400/4).

**Results**: Run completed (37 min).  Auto-picked cluster: pair pins 3/3,
hotspot B20-26 all >5 Å.  One cluster rotated far in the WRONG direction —
a different loop lay parallel to the hotspot strand.  The receptor's
design-region BACKBONE prevented the intended rotation.

**Conclusion**: Soft AIR + sidechain stripping insufficient.  The design-region
backbone must be physically removed to allow the rotation.  Proceed to E4.

---

## Active proposals

---

### E4 — Carve-and-Place (backbone removal + frozen receptor + rigid-body placement)

**Status**: Proposed, not yet implemented.

**Goal**: Allow HADDOCK to find a biologically optimal effector placement at
the truncated receptor without any steric obstruction from design-region
backbone or sidechains.  The user's ChimeraX pose provides the starting hint.
Both proteins are rigid throughout (no flexref/emref).

**Hypothesis**: The design-region backbone is the steric barrier preventing
HADDOCK from reaching the user's intended interface.  Removing all design-region
atoms from the receptor PDB (not just sidechains) and locking the truncated
receptor at the user's pose orientation will allow HADDOCK to optimally place
the effector via rigid-body sampling, guided by the restraint set.

**Rationale**:

- E2 (sidechain stripping) was insufficient because HADDOCK still
  sees backbone atoms in the design region that the effector cannot
  penetrate, even at very tight sampling.
- After stripping, HADDOCK found alternative poses (e.g., opposite receptor
  rotation) because the backbone still prevented the desired rotation.
- Removing backbone entirely creates space for the effector to reach the
  intended receptor face.
- Freezing the remaining receptor is not restrictive: relative motion between
  frozen receptor and mobile effector spans the same configuration space as
  letting both move (rigid-body equivalence).  The restriction is to what
  the user *wants* — the user's pose orientation for the fixed regions.
- HADDOCK never produces output clashes, so the absence of design-region
  backbone means no backbone-mediated steric barriers exist.
- Both proteins rigid: no flexref/emref.  Their crystal-structure
  coordinates are correct; sidechain optimisation is irrelevant given
  downstream redesign.

**Experimental plan**:

1. **Clash detection at the hint pose**: Load the user's ChimeraX-aligned
   complex.  Compute heavy-atom pairwise distances (receptor vs effector).
   Identify all receptor residues with any atom within a threshold (start:
   2.0 Å hard clashes + 3.5 Å near-contacts) of any effector heavy atom.
   These form the initial carve set.

2. **Expand carve set**: Union of (a) identified clashing residues + (b) the
   existing contig de novo regions (33-49 + 69-78), since those will be
   redesigned regardless.  The carve set may be a subset or superset of the
   contig design region.

3. **Build truncated receptor PDB**: Remove ALL atoms (backbone + sidechain)
   for carve-set residues.  The remaining receptor contains only the fixed
   contig segments (A1-32 + A50-68).  These are kept at the user's hint
   coordinates.

4. **Handle chain breaks**: The two fixed segments are no longer covalently
   connected.  Options: (a) treat as two separate molecules in HADDOCK3
   `molecules = [...]`; (b) keep as one chain with an explicit chain-break
   record.  Option (a) is more robust — HADDOCK3 supports multi-molecule
   complexes natively.

5. **Pair pins**: If pair-pin residues are in the carve set (e.g., A71-73
   which are in DR2 / residues 69-78), convert them to distance restraints
   to virtual Cα positions derived from the user's hint pose.  CNS supports
   distance-to-coordinate restraints.  If pair pins reference FIXED residues
   (A1-32 or A50-68), they work unchanged.

6. **Configure HADDOCK3**:
   - Reduced sampling (200–500, not 10000) to exploit the hint-initialisation.
   - Set rigid-body initial coordinates from the hint complex.
   - Skip or stub-out flexref/emref (mol1_fix=true + mol2_fix=true, or
     remove those stages from docking.cfg).
   - Keep clustfcc + seletopclusts + caprieval as normal.

7. **Pair distance correction**: Change pair-pin target from `2,2,4` (CA-CA
   distance window 0-6 Å) to `5,2,3` (window 3-8 Å), matching biologically
   realistic Cα-Cα contact distances for strand-mediated hydrogen bonding.

8. **Post-docking splice**: Before BUILD_CONTIGS, splice placeholder residues
   back into the receptor PDB at carve-set positions (Cα-only, with the
   correct residue count for contig compatibility), so RFDiffusion gets the
   correct protein length and contig string.

**Predicted outcome**: HADDOCK finds effector poses where B20-26 lands on
the DR1 side of the receptor (the only side with actual backbone to bind
against), pair pins satisfy at biologically realistic distances.

**Failure modes**:

- The two fixed-segment "molecules" may be too small / sparse for HADDOCK
  to compute meaningful BSA, leading to degenerate clustering.
- If only a handful of pair-pin-compatible poses exist, the cluster may have
  very few members.  Not a scientific problem but may confuse the
  cluster-size diagnostics.
- If the user's hint is ≥30° off the true orientation, HADDOCK's reduced
  sampling may not escape the local basin.  Mitigation: increase sampling
  or run blind first (see E5).

---

### E5 — Iterative Carve-and-Search (fallback for E4)

**Status**: Proposed, not yet implemented.

**Goal**: Discover the minimal carve set that allows HADDOCK to reach the
user's intended pose, when the initial single-pass carve (E4) fails.

**Hypothesis**: The E4 hint-pose clash detection may miss residues that are
not clashing at the hint but become clashing when the effector rotates to
its optimal position against the truncated receptor.  A second pass would
detect and remove these.

**Rationale**: HADDOCK never produces output clashes (clarification #2), so
any "new" clashes after E4 won't appear in the output — but they might be
biologically real if the optimal effector placement puts the effector face
against a receptor region that was not clashing in the raw hint pose.

**Experimental plan**:

1. Run E4.
2. Compute the RMSD of the E4 auto-picked cluster to the user's hint complex
   (receptor fixed-region RMSD + effector RMSD).
3. If RMSD < threshold (say 3 Å) on both: accept E4 result.
4. If RMSD > threshold: inspect which receptor residues are within 3.5 Å of
   the effector in the E4 output that were NOT in the original carve set.
   Add these to the carve set and repeat from E4 step 3.
5. Maximum 2-3 iterations.

**Note**: Per the user's empirical observation that HADDOCK never produces
output clashes, the iterative expansion may rarely be triggered.  Include
as a defined fallback, not the primary path.

---

### E6 — Pose-Driven Design Region (eliminate manual contig string)

**Status**: Proposed, under discussion.

**Goal**: Remove the requirement for the user to manually specify a contig
string defining the design region.  Instead, the design region is
automatically derived from the clash set at the user's input pose.

**Hypothesis**: The residues that need to be redesigned are exactly those
that clash with the effector in the user's intended pose.  If the user can
see and approve the auto-derived design region in ChimeraX before submitting,
this is a sufficient specification.

**Rationale**: The current pipeline requires the user to provide two
consistent inputs: (1) a ChimeraX-aligned complex encoding the intended
pose, and (2) a contig string encoding the design region.  These should be
the same decision — "which receptor residues need to change to accommodate
this pose?" — but are expressed separately.  Auto-derivation from clashes
unifies them into a single input decision.

**Potential workflow**:

1. User places both proteins in ChimeraX.
2. Pipeline (or a ChimeraX plugin) identifies clashing receptor residues in
   real time.
3. User sees the clashing residues highlighted as the proposed design region.
4. User adjusts the placement until the design region matches their biological
   intent (e.g., "I want the de novo loop to be roughly 15 residues").
5. User submits — pose + auto-derived design region.
6. Contig string is AUTO-GENERATED from the carve set.

**PhD supervisor critique — drawbacks**:

1. **Loss of biological control over the design region.** Not all redesign
   decisions are clash-driven.  A user might want to redesign a
   surface-exposed loop for stability, a conserved patch for specificity, or
   avoid redesigning a structurally important motif even if it clashes.
   Constraining the design region to clashes alone removes this agency.

2. **Unpredictable design region size.** RFDiffusion runtime, success rate,
   and convergence scale strongly with design region size.  A 1 Å nudge in
   ChimeraX can add or remove 5 clashing residues.  The user loses the ability
   to say "I have GPU budget for a 15-residue redesign."

3. **Multi-region topology is lost.** The current contig syntax encodes two
   separate de novo regions of specified lengths.  Auto-detection from clashes
   produces a set of residues, not a structured multi-region topology.

4. **Threshold ambiguity.** What counts as a clash: 2 Å? 3 Å? 4 Å?
   Replacing "explicit residue list" with "clash threshold" doesn't reduce
   decisions — it relocates one in a less biologically interpretable form.

5. **Coupling pose quality to design region quality.** A poor ChimeraX
   placement produces both a poor pose AND a poorly-derived design region.
   With explicit contigs these are separable; the user can fix the pose
   without changing the design region.

6. **Provenance and reproducibility.** The design region would be an emergent
   property of a PDB file + a threshold.  Re-execution with a slightly
   different input PDB produces a different design region.  Explicit contigs
   are version-controlled and self-documenting.

7. **Non-interface redesign intent is inaccessible.** Redesign for binder
   stability, reduced immunogenicity, specificity against a panel of effectors
   — none of these are derivable from clashes at a single pose.

**Mitigations that make this worth pursuing**:

- **Auto-detect as default, manual override available**: a YAML param
  `haddock_design_region_overrides: {add: [90,95], remove: [45]}` lets the
  user adjust the auto-detected set.
- **Pre-submission visualisation** (the ChimeraX plugin idea): the user SEES
  the carved region before submitting.  Catches "oh, that's only 8 residues,
  I expected 15."  This is the most important mitigation.
- **Two thresholds**: 2 Å (hard clash, definitely must carve) + 3.5 Å (soft
  near-contact, propose-but-require-user-confirmation).
- **Contig string preserved as canonical output**: the auto-detection populates
  a contig string that is saved in params.yml for provenance.  Auto-detection
  populates it; the user can still inspect and edit.

**Verdict**: Pursue as a longer-term UX improvement, built on top of E4.
Not a replacement for explicit contigs in the near term.  The ChimeraX plugin
for real-time clash visualisation is a worthwhile standalone tool even if the
fully-automated contig derivation is never adopted.

---

## Discarded ideas

### D1 — Dense pseudo-pair restraints

**Why**: HADDOCK overrides restraints to avoid clashes.  Providing many
restraints derived from the intended pose does not force HADDOCK into that
pose — it will find the best satisfied-restraint configuration that avoids
sterically impossible overlaps.  The effector will bend or rotate to the
nearest restraint-compatible clash-free pose, which may not be the intended
one.  Rejected based on user's empirical observation.

### D2 — Drop HADDOCK entirely (pose-first workflow)

**Why**: The user's ChimeraX placement is a rough hint, not the final answer.
Orientation and placement may both be suboptimal; contact distances may be
biologically implausible.  Taking the hint directly into RFDiffusion without
any refinement step gives up the scientific basis the user explicitly wants.

### D3 — Receptor sidechain refinement / flexref

**Why**: Sidechain rotamers in the design region will be completely replaced
by RFDiffusion + ProteinMPNN.  Sidechain rotamers in the fixed region come
from the crystal structure and are correct by definition.  HADDOCK's flexref/
emref refinement of either is wasted computation.  Explicitly rejected by
the user.

### D4 — Effector-only mobile with random initialisation

**Why**: Equivalent to "both mobile" in rigid-body docking (clarification #1).
The "freeze the receptor" decision is about limiting sampling to a hint-initialised
basin and skipping receptor refinement — not about changing the configuration
space.  A hint-initialised run achieves the same relative pose quality.

---

## Literature context summary

| Method | Pose-hint input | Auto design region | Backbone removal | Both proteins rigid |
|---|---|---|---|---|
| RFDiffusion partial diffusion | ✓ (rough scaffold) | ✗ (user-specified) | ✗ | ✗ (diffusion) |
| RFantibody | ✗ (full dock sampled) | Partial (CDR anatomy) | ✗ | ✗ (framework fixed) |
| BindCraft | ✗ (hotspots only) | ✗ (target-driven) | ✗ | ✗ |
| PPDiff | ✗ | ✗ | ✗ | ✗ |
| AnchoredDesign (Rosetta) | Partial (anchor graft) | ✗ (user-specified loop) | ✗ | ✗ (flexbb) |
| **E4 (this work)** | **✓** | **✓ (from clashes)** | **✓** | **✓** |

**E4 occupies a combination of features that no single published method
addresses.**  The closest published analogue is RFdiffusion partial diffusion
(pose-hint input + refinement) but it differs in that the design region
boundary is pre-specified not auto-derived, and the proteins are not rigid
(diffusion resamples backbone).  The clash-driven carve is the novel
methodological contribution if this approach is published.

---

## Decision log

| # | Decision | Status | Source |
|---|---|---|---|
| 1 | HADDOCK's role = geometric placement, not affinity prediction | Locked | Session 7 A129 |
| 2 | Cluster ranking primary key = pair_contact_fraction, secondary = BSA | Locked | Session 7 A154 |
| 3 | Contig string no longer drives HADDOCK AIR derivation | Locked | Session 7 A131 |
| 4 | Both proteins rigid throughout HADDOCK (no flexref/emref) | Locked | This document |
| 5 | HADDOCK never produces output clashes (empirical) | Locked | This document |
| 6 | Rigid-body freezing is equivalent to letting both move (rigid-body equivalence) | Locked | This document |
| 7 | Backbone removal from design region is required (sidechain stripping insufficient) | Locked | E2/E3 results |
| 8 | Pair pin distances should be biologically realistic (5 Å target, not 2 Å) | Proposed | This document |
| 9 | ChimeraX plugin for real-time clash visualisation is the right UX approach | Proposed | This document |

---

*See also*:
- `notes/design_audit.md` Session 7 (Q129–Q157) for the HADDOCK module architectural decisions
- `notes/SESSION_HANDOFF.md` for current branch state and immediate next steps
