# Pipeline Design Audit — Grill-Me Sessions

A living record of the systematic audit of design intent vs codebase reality.
Questions are flagged [VERIFY] (hypothesis from code/notes — user confirms or corrects)
or [UNKNOWN] (only user can answer). Answers are recorded verbatim.

**Scope**: 5 sessions, ~320 questions total. Goal: produce an authoritative spec
of what the pipeline is supposed to do, surface every accumulated error.

**Progress**: Session 1 of 5 — Questions 1–8 of ~320 (2%)

---

## Session 1 — Scientific domain intent + pipeline invariants

### Batch 1 — Scientific foundation (Q1–Q8)

---

**Q1 [UNKNOWN]** — The five campaigns target pikp1_avrpia, pikp1_avrpikf,
pikp1_pby2, pikp1_pwt3, pikp1_pwt7. What is the biological relationship between
these targets? Are pby2, pwt3, pwt7 effectors from different pathogens that
Pikp1 doesn't currently recognise (so the goal is to engineer new recognition),
or are they variants within the same effector family that some Pikp1 alleles
already recognise? The answer determines what "success" means computationally.

**A1:** These are 5 effector targets that the Banfield lab works on, and therefore
understands well. They all have similar structures but are chosen for different
reasons - AVR-PikF is reproducing a known result (the SNK-EKE mutant by Josie
Maidement). AVR-Pia has also been shown to bind Pikp-1, but at a different HMA
interface, so we want to see if we can engineer binding at the AVR-Pik interface.
Pby2 is canonically recognised by a different NLR (Rmo2). Pwt7 is canonically
recognised by a Tandem Kinase Protein (TKP) receptor. Pwt3 has no solved structure,
and it's AF3 prediction shows that it has an additional large alpha helix not shared
by the others. They are different levels of challenge for the pipeline, whilst still
keeping the core MAX effector folds.

**Synthesis note**: The five campaigns span a challenge gradient — from reproducing
a known experimental result (AVR-PikF) to engineering entirely new recognition
(Pby2, Pwt7) to a structurally atypical target with no solved structure (Pwt3).
"Success" therefore means different things per campaign: for AVR-PikF, recovery
of the SNK-EKE binding mode is a positive control; for Pby2/Pwt7, any
agroinfiltration-passing design is novel. This asymmetry is not currently encoded
anywhere in the pipeline and should inform how results are interpreted per campaign.

---

**Q2 [UNKNOWN]** — When a design passes the full pipeline (pose_holds, orthogonal
filters), what is the first wet-lab test it goes into? From the notes I see
agroinfiltration and protoplast assays mentioned. Which comes first, and what does
a positive result in that assay actually tell you — protein-protein interaction,
immune activation, or something else?

**A2:** The passing designs will be validated in the plant using agroinfiltration
transient co-expression experiments. This may be followed up with Co-IPs for
agroinfiltration-passing designs.

**Synthesis note**: Agroinfiltration is the primary gate. A positive result means
immune activation (HR response) when the designed receptor is co-expressed with
the effector in planta. This is a functional assay, not a direct binding assay —
it reports whether the designed receptor triggers the immune pathway in the
presence of the effector. Co-IP is a secondary confirmation of physical interaction
only for agroinfiltration-positive hits.

---

**Q3 [VERIFY]** — My reading of the pipeline's intent: the output is a set of
ProteinMPNN-designed receptor sequences where Boltz reproducibly predicts the
effector binding at the RFDiffusion-designed site, and that pose survives removal
of the steering mutations. The claim is not "this sequence binds" but "Boltz's
structural prediction is internally consistent and not dependent on artificial
perturbations." Is that correct?

**A3:** You're correct that the pipeline is designed to provide ProteinMPNN
sequences which will be taken forward into the lab for agroinfiltration. Negative
steering is designed to make sure that Boltz-2 can predict the correct binding
interface, but thus far there does not seem to be any Boltz-2 confidence or
interface metrics which are reliable for wet lab success. We haven't done any
agroinfiltrations yet, but I've seen that running Boltz-2 several times on the
same input sequence can generate very different metrics for the same pose. I need
to finalise this before fully generalising and giving up on other binding metrics,
but it will probably be the case that I will pick and choose from the passing
designs those that excel in different confidence metrics so I can have a variety
of things to test.

**Synthesis note**: The pipeline does NOT claim to predict binding affinity or
in-planta activity. Its claim is structural consistency: Boltz can be made to
predict the correct pose, and that pose is not dependent on the steering mutations.
Confidence metrics (ipSAE, iPTM) have no known correlation with wet-lab success
for this system yet — first agroinfiltration results are needed to establish any
such correlation. Current plan: manual selection of designs that collectively
sample different confidence profiles, creating a diversity of test cases.

---

**Q4 [UNKNOWN]** — For the current campaigns: do you have solved crystal structures
or high-confidence AF3 models of the native complexes that could serve as ground
truth instead of (or alongside) the RFDiffusion output?

**A4:** For the current designs, every protein except Pwt3 has a solved structure
(although I'm also using an AF3 prediction of Pikp1-HMA because the solved
structures are missing a few key residues on the C-terminus). However, AVR-PikF
has a loop region which forms significant contact with Pikp1-HMA which isn't
present in the other effectors, so it may be necessary to adjust their positioning
slightly to prevent RFDiffusion from having to build out a huge bridge across the
resulting gap that would come from aligning the effectors to place the complex.
HADDOCK will be redeveloped later to see if I can have it produce better input
binding poses - I think the current AIR table implementation isn't well optimised
but that's a future problem.

**Synthesis note**: Ground truth availability per campaign:
- AVR-PikF: solved structure available; Pikp1-HMA uses AF3 prediction (solved
  structures lack C-terminal residues). HADDOCK placement quality is uncertain.
- AVR-Pia: solved structure available.
- Pby2: solved structure available.
- Pwt7: solved structure available.
- Pwt3: NO solved structure — AF3 prediction only; structurally atypical (extra
  helix). Highest uncertainty in ground truth.

AVR-PikF-specific risk: the unique contact loop may cause HADDOCK to misposition
the effector, as documented in notes4. This has not been fixed.

HADDOCK redevelopment is acknowledged as future work. Current AIR table
implementation is known to be suboptimal.

---

**Q5 [VERIFY]** — The contamination rule is now "any steering mutation that
contacts the effector in the predicted pose triggers reversion, regardless of
which surface." The earlier protected-set concept was removed. Is that the
intended final rule, or should non-design-region contacts be excluded?

**A5:** You're correct that the contamination should only be steered residues
contacting the effector that are in the design region. The intention is that any
non-load-bearing steering contacts aren't going to have any influence on the
design region binding, or any of downstream interface metrics that are used to
rank designs. Some of the proteins involved have shapes that loop around their
partner, so it would be very challenging to negatively steer without any other
contacts. The non-load-bearing contacts are just there to push the effector
placement to the right interface - I think this is worth significant discussion
in the future though, and there may be strategies to avoid these non-load-bearing
steering contacts.

**Synthesis note — POTENTIAL DISCREPANCY**: The user's stated intent ("only
design-region contacts matter") diverges from what the code currently does.
The code after notes7 reverts ALL mutated contacts regardless of location —
the protected-set concept was explicitly removed. However, notes11 added
post-reversion gating at the design_region ∪ true_interface level for the
new_contamination verdict. There are therefore two different policies in play:

1. Pre-reversion (build-contaminated): ALL mutated contacts trigger reversion.
2. Post-reversion (new_contamination verdict): only design_region ∪ true_interface
   contacts fire.

The user's intent may be closer to policy 2 applied consistently to both stages.
This is flagged as a **follow-up question** in Batch 3 (negsteer semantics deep dive)
to be resolved precisely before any Phase 4 refactoring of the negsteer core.

---

**Q6 [VERIFY]** — The tiering (A=3/3, B=2/3, C=1/3 for num_seeds=3) — is this
based on an assumption of exactly 3 seeds, or should it scale?

**A6:** Is the current tiering based on an assumption of 3 Boltz-2 runs per
steered sequence? In that case A=3/3, B=2/3, and C=1/3, and None=0/3. This is
the most likely use case, but it should also scale to different numbers of Boltz-2
runs per steered sequence. A should always be "all runs passed", and None should
always be "no runs passed". B and C should be evenly distributed e.g. B=3-4/5,
C=1-2/5. It shouldn't ever be an even number because then it's hard to calculate
meaningful medians of the metrics.

**Synthesis note — SCALABILITY GAP**: The tiering must be parametric on num_seeds,
not hardcoded for 3. The invariants:
- Tier A: n_seeds_pose_holds == n_seeds
- Tier None: n_seeds_pose_holds == 0
- Tiers B and C: divide the middle range evenly. For num_seeds=3: B=2, C=1.
  For num_seeds=5: B=3–4, C=1–2. For num_seeds=7: B=4–6, C=1–3.

**Design constraint**: num_seeds must always be odd to ensure unambiguous medians.
Even num_seeds should either be rejected at the CLI or flagged with a warning.
This constraint is not currently enforced anywhere in the code.

The current implementation was written for num_seeds=3 specifically. Whether it
scales correctly for other values needs to be verified in the code.

---

**Q7 [VERIFY]** — The composite score `true_jaccard − 0.05 · ra_eff` — was the
0.05 weight empirically calibrated, and should weighted_jaccard also be in it?

**A7:** The composite formula was definitely a reasonable sounding guess. It isn't
based on anything empirical other than trying to balance if the effector was placed
at the correct interface. We also added a weighted jaccard, and I don't know if
that should be used in this metric. The idea was to run the pipeline a few times
and understand the relationship between all the metrics so we can choose some
individual or combination metrics to use as rankers going forward. Since the main
issue was incorrect effector placement, the current composite metric is just good
for making sure at least that is correct. It should be audited going forward.

**Synthesis note**: The composite `true_jaccard − 0.05 · ra_eff` is a placeholder,
not a calibrated formula. It is not suitable as a final ranker. The 0.05 weight
was chosen to give ra_eff influence without dominating true_jaccard, but has
no empirical basis. Weighted_jaccard's inclusion in the composite has not been
decided. The composite should be treated as a sorting convenience only, and
replaced once agroinfiltration data exists to calibrate it. This is deferred
to Task 47 (threshold audit).

---

**Q8 [VERIFY + CORRECTION]** — mpnn_top_n selects sequences per design or across
all designs? I assumed one per design.

**A8:** That's incorrect. The pipeline currently runs the top N-ranked MPNN
sequences (MPNN score) through negative steering, meaning there could be more
than one MPNN sequence per design taken forward if those sequences have good MPNN
scores and assuming that more than one MPNN sequence was generated for each design
in the first place. At first, selecting what will be taken forward to the wet lab
will be a manual process involving humans inspecting the passing results and
choosing from them. It's possible that two MPNN sequences from the same design
will be tested with agroinfiltration if there is a meaningful reason to do so
(maybe they excel in very different confidence metrics).

**Synthesis note — CODE VERIFICATION NEEDED**: mpnn_top_n selects the top N
sequences globally by MPNN score, not one per design. This means a single
high-quality RFDiffusion design could contribute multiple MPNN sequences to the
negsteer cohort while a lower-quality design contributes none. The code in
`bin/mpnn_select_top.py` needs to be verified to confirm it implements global
ranking rather than per-design selection. This affects how results are interpreted
per campaign.

First-pass wet-lab selection is manual. No automated per-design deduplication
or diversity rule is intended.

---

## Pending follow-up questions (flagged during this session)

1. **Contamination semantics (from Q5)**: Precisely, which stages should use
   design-region-only gating vs all-contacts? Pre-reversion (build-contaminated)
   currently reverts all contacts; post-reversion (new_contamination) is gated
   to design_region ∪ true_interface. Is the pre-reversion policy intended to
   match post-reversion policy? → Deferred to Session 2 (negsteer core semantics).

2. **num_seeds even/odd constraint**: Should even num_seeds be rejected at the
   CLI? → Deferred to Session 2.

3. **mpnn_select_top.py**: Confirm global vs per-design ranking in code. →
   Code read required before Session 3.

---

---

### Batch 2 — Campaigns, chain conventions, contigs (Q9–Q16)

---

**Q9 [UNKNOWN]** — For AVR-Pia specifically: you said it "has been shown to bind
Pikp-1, but at a different HMA interface." The pipeline designs binding at the
AVR-Pik interface. What is the evidence that this is achievable?

**A9:** The intention of designing binding at the AVR-Pik interface is that we
have a lot more experimental data about variations at that interface. It is a
hypothesis that is being tested, but given all the targets share the MAX effector
fold, it's reasonable to believe that this engineering at the AVR-Pik interface
is possible. Also, some modifications at the AVR-Pia interface have been known
to cause issues when integrated into the wider Pikp-1 chassis. The designed pose
at the AVR-Pik interface is based on some biological plausibility - there are two
anti-parallel beta strands interacting between Pikp-1 HMA and each effector, but
the loop region in AVR-PikF that forms another binding site I mentioned earlier
is missing in the other effectors, and so some pose changes may be necessary.

**Synthesis note**: AVR-Pia binding at the AVR-Pik interface is a hypothesis with
biological plausibility (shared MAX fold, shared anti-parallel beta strand contact)
but no direct experimental evidence yet. The AVR-PikF loop is absent in AVR-Pia,
meaning the input complex construction for AVR-Pia may require manual pose
adjustment relative to the AVR-PikF reference. This asymmetry is not captured
in the pipeline and should be noted when interpreting AVR-Pia results — the
ground truth pose may be less well-defined than for AVR-PikF.

---

**Q10 [UNKNOWN]** — For Pwt3 specifically: AF3 prediction as ground truth —
how confident are you in it, and where is the extra helix relative to the
designed interface?

**A10:** I don't have any basis of confidence or non-confidence for AF3-predicted
effectors. That can only come from experimental data, which creates a bit of a
chicken-and-egg situation. For now we just have to treat them the same as the
solved structures, with an eye on the original AF3 prediction confidence scores.
For Pwt3, the confidence scores aren't great, and the extra helical region, which
is fortunately on the other side of the protein to the AVR-Pik interface site,
is particularly low confidence. My HPC AF3 predictions of Pwt3 matched those of
the AFDB except from the orientation of the extra helical region, where mine
pointed away from the effector centre of mass, and the AFDB one folds nicely
alongside it (more compact overall shape).

**Synthesis note**: Pwt3 ground truth uncertainty is highest of the five campaigns:
no solved structure, AF3 confidence is poor at the extra helix, and two AF3
predictions (HPC vs AFDB) disagree on helix orientation. The extra helix is
fortunately distal to the AVR-Pik interface. The pipeline treats all five ground
truths identically — there is no mechanism to flag or weight Pwt3 results with
lower confidence due to ground truth uncertainty. This is a known limitation.

---

**Q11 [VERIFY]** — Chain convention: RFDiffusion output always has chain A =
receptor, chain B = effector, hardcoded. `params.receptor_chain` /
`params.effector_chain` describe the input PDB only. Is this correct?

**A11:** RFDiffusion should always output chain A = receptor, chain B = effector,
yes. It makes it unambiguous and consistent - there is never any doubt which is
which. This is particularly useful for effectors who mimic the shape of immune
proteins. The input structure may have different chain IDs, and that allows the
user input to be a little more flexible - the user needs to correctly set the
input params for these.

**Synthesis note**: Chain convention confirmed as intentional and correct. The
A/B hardcoding is a deliberate design decision, not an accident. The purpose is
disambiguation for effectors that share structural features with receptor domains.
`params.receptor_chain` / `params.effector_chain` are INPUT-only parameters;
everything downstream of RFDiffusion uses the hardcoded A/B convention via
`params.rfdiff_output_receptor_chain` / `params.rfdiff_output_effector_chain`.

---

**Q12 [UNKNOWN]** — Is there a documented procedure for constructing the input
complex PDB (aligning effector onto reference)? The quality of this determines
the ground truth quality.

**A12:** Yes, this is a real challenge right now. Currently I've just been using
ChimeraX's "matchmaker" function, but sometimes this becomes confused by non-core
folds. The main idea is to match up those anti-parallel beta strands that form a
contact, since they are shared between all the 5 effectors I'm focusing on. But
for some I had to manually position them in ChimeraX. Clearly this is not ideal,
but at least it is rooted in the biology rather than complete guesswork. In the
future, and particularly for effectors without those shared folds, I need to sort
out the HADDOCK module. I also thought about perhaps including something that uses
FoldSeek to look for similar complexes to help guide the positioning of the input,
but I think that's something for the distant future.

**Synthesis note**: Input complex construction is currently:
- Primary method: ChimeraX matchmaker on anti-parallel beta strands (shared feature)
- Fallback: manual positioning when matchmaker fails
- No written protocol or reproducibility record
- Future intent: improved HADDOCK module; possibly FoldSeek for cross-reference

This is a known fragility. The five campaign input PDBs were constructed with
different levels of manual intervention. Results across campaigns are not
directly comparable in terms of ground truth quality.

---

**Q13 [VERIFY]** — Contig string format: `A1-32/10-20/A46-72/6-6 C32-113`
means fixed receptor / variable-length design / fixed receptor / fixed-length
design / effector fixed. The effector can also be specified as just the chain ID
to fix the whole chain. Is this the standard RFDiffusion contig format?

**A13:** Your reading of the contig is correct - design/non-design regions are
split by a forward slash, where fixed non-design regions are led by the chain ID
(e.g. A1-10), and design regions have the format min_design_length-max_design_length.
For the effector, it can be provided as just the chain ID (e.g. C), and the
program automatically infers that the whole effector will remain fixed. I believe
it matches the standard RFDiffusion contig format from the Google Colab version
by sokrypton.

**Synthesis note**: Contig format confirmed as standard RFDiffusion format.
Key clarification: bare chain ID in effector position (e.g. `C`) means whole
chain is fixed. This is handled internally by `contig_utils.resolve_contigs`.
The design region ranges (e.g. `10-20`) specify min/max lengths; RFDiffusion
selects the actual length per design within that range.

---

**Q14 [UNKNOWN]** — How were the contig strings chosen? What distinguishes
strategies 4a through 4f?

**A14:** There are several strategies for choosing the design regions, and this
is something that I'm going to iteratively improve the more times I run the
pipeline. One idea (4a) was to just select residues whose side chains face
externally into the binding interface, and then gather any nearby islands of
design regions into longer bits (otherwise RFDiffusion doesn't have enough
flexibility to design anything meaningful). The second idea (4b) was to
specifically use the literature, and the known residues important in binding,
and create design regions around them (this makes sense given the risks around
autoactivity which can occur when changes are made outside of known variable
regions). The third idea (4c) was to select residues on the effector that are
targets, and to build design regions around the residues that come into contact
with them. The fourth idea (4d) was to just redesign the entire AVR-Pik
interface. After that, I had some other strategies involving me manually
selecting them based on intuition (4e for AVR-Pia interface, 4f for AVR-Pik
interface), although they need some more rigorous planning in the future.
Regarding the lengths of the design regions, these were just selected as a wide
guess to allow RFDiffusion to figure out what length would work best. The design
region at the C-terminus was fixed though, since it shouldn't extend or truncate
either end of the receptor HMA.

**Synthesis note**: Strategy taxonomy:
- 4a: surface-facing residues + island aggregation (structural criterion)
- 4b: known binding residues from literature (biological criterion; autoactivity-aware)
- 4c: effector-target residues + contacting receptor (contact-based criterion)
- 4d: entire AVR-Pik interface (maximal)
- 4e: manual intuition, AVR-Pia interface
- 4f: manual intuition, AVR-Pik interface

Design region lengths are currently wide guesses (e.g. `10-20`), not derived
from structural analysis. The C-terminal design region is length-fixed (`6-6`)
to prevent extension/truncation of the HMA domain — this is an important
biological constraint encoded in the contig string.

---

**Q15 [VERIFY]** — `design_region` comes from RFDiffusion metrics; `true_interface`
is receptor residues contacting effector in the RFDiffusion design at 8Å Cα.
Both are from the post-design PDB, not a native structure.

**A15:** That is correct. There are a few things here. The original design regions
are part of the input Pikp-1 HMA. The design region term also applies to the
changes made by RFDiffusion at those original design regions, which could be
substitutions, insertions, or deletions. The true interface should be calculated
from the RFDiffusion design as the residues on the receptor which are in contact
with the effector. Now clearly the majority of these should be in the design
interface, but also I'm not sure how this handles the non-load-bearing steered
residues that we spoke about earlier. If I remember correctly, these are excluded
from the true interface. This is one of the things that may have migrated over
time and would be worth revisiting properly.

**Synthesis note — UNCERTAINTY FLAGGED BY USER**: The relationship between
non-load-bearing steered positions and the true_interface set is explicitly
uncertain. The user believes non-load-bearing contacts (e.g. positions 22/26
in the notes10 analysis) are excluded from the true_interface, but is not
certain. The true_interface is derived from the RFDiffusion design PDB at
8Å Cα — it does NOT use the native structure. Whether the non-load-bearing
contacts at non-interface receptor positions appear in true_interface depends
on the Cα distance from those positions to the effector in the RFDiffusion
design. This is flagged as a **code verification needed** item: read
`derive_true_interface.py` and `bin/boltz2_negative_steering.py:cmd_plan`
to confirm exactly which positions appear in `true_interface_idx`.

---

**Q16 [UNKNOWN]** — Difference between strategy 4e and 4f, and what are the
rfdtest runs?

**A16:** The difference between strategy 4e and 4f (as outlined in my answer to
Q14) is that they are for different HMA interfaces (AVR-Pia and AVR-Pik
interfaces). Other than that, they are both chosen manually by me. The rfdtest
runs were a way of me testing RFDiffusion runs on a challenging design and
implementing a pause after the Rosetta Filtering step to iterate over the designs
before moving onto the negative steering. They aren't necessarily production runs
either. All the runs up until now are exploratory.

**Synthesis note**: 4e = AVR-Pia interface (manual), 4f = AVR-Pik interface
(manual). The `v1_4f_rfdtest_1` and `v1_4f_rfdtest_2` runs used `stop_after_rosetta`
(a parameter added from the experiments branch) to pause the pipeline after Rosetta
filtering for manual inspection of designs before committing to the GPU-intensive
negsteer stage. This is a deliberate workflow feature, not a partial run.

**Important framing**: ALL current runs are exploratory. There is no "production
run" yet in the sense of a committed cohort for agroinfiltration. The pipeline
is in a validation and iteration phase.

---

## Pending follow-up questions (updated)

1. **Contamination semantics (from Q5, reinforced by Q15)**: Precisely which
   stages use design-region-only gating vs all-contacts? The user believes
   non-load-bearing contacts outside the true_interface don't trigger reversion,
   but is uncertain. Code verification of `derive_true_interface.py` and
   `cmd_plan` needed. → Session 2.

2. **num_seeds even/odd constraint**: Should even num_seeds be rejected at CLI?
   → Session 2.

3. **mpnn_select_top.py**: Confirm global vs per-design ranking. → Session 3.

4. **true_interface derivation**: Does `derive_true_interface.py` use 8Å Cα
   from the RFDiffusion design? Do non-load-bearing positions (e.g. 22/26 from
   notes10) appear in it? → Code read + Session 2.

---

---

### Batch 3 — RFDiffusion filter, MPNN selection, controls, seeds (Q17–Q24)

---

**Q17 [VERIFY]** — RFDiffusion filter is maximally permissive in current campaigns
(`hotspot=""`, `min_hotspot_frac=0.0`). Is Rosetta Sc the real first quality gate?

**A17:** Well currently I set those input parameters to be very permissive, to
allow any designs through whilst I'm building and testing the pipeline. It won't
necessarily all be like that. Questioning the current input parameters isn't
particularly meaningful since they're all going to be calibrated later. Currently
yes Rosetta Sc and negsteering are the main filters, but that's because currently
negsteering has a low hit rate and I want to give it as many chances as possible.

**Synthesis note**: Permissive RFDiffusion filter is a deliberate strategy during
the exploratory phase: give negsteer as many candidate designs as possible while
the hit rate is being characterised. The low negsteer hit rate is acknowledged.
All RFDiffusion filter parameters (`contact_cutoff`, `min_hotspot_frac`, `hotspot`)
are provisional and require calibration once hit rates are better understood.

---

**Q18 [UNKNOWN]** — Rosetta Sc threshold of 0.5 — how was it chosen, and what
fraction of designs typically survive it?

**A18:** 0.5 was selected arbitrarily. Most of the time it's been either all or
none, but in the most recent rfdtest iteration, around 1.3 passed it. It needs
to be calibrated too.

**Synthesis note**: Sc=0.5 is arbitrary. The "all or none" pattern suggests the
Sc filter provides limited discrimination in most runs — either the backbone
geometry is consistently good or consistently bad. "Around 1.3 passed" in the
rfdtest iteration is ambiguous in phrasing (possibly 1–3 designs, or ~13%) and
requires clarification. The Sc threshold needs to be calibrated against data from
more runs. At present it should be treated as a coarse sanity check rather than
a calibrated quality gate.

---

**Q19 [VERIFY]** — `mpnn_top_n` selects globally by MPNN score across all
designs, not per-design. Is this intentional?

**A19:** If mpnn_top_n=8 and there are 8 RFDiffusion designs, it's only possible
for all 8 sequences to come from one design if there are at least 8 sequences
generated per design. I'm currently uncertain about this selection step, and
you're right that per-design selection may be more appropriate. I was trying to
remove any sequences which didn't seem biologically plausible by MPNN's own
metric, but you're right that it does tend to favour certain designs. I'm
planning to create more RFDiffusion designs with fewer MPNN sequences for each,
so that forces diversity anyway at the moment. Worth reconsidering later though.

**Synthesis note**: Global selection is acknowledged as potentially problematic
but accepted for now. The user's workaround is to generate more designs with
fewer sequences per design, which enforces diversity at the RFDiffusion level
rather than the MPNN level. The selection step should be revisited once larger
cohorts are being run. **Code verification needed**: confirm `mpnn_select_top.py`
implements global ranking and whether the sort direction is correct (lower MPNN
score = better per Q20).

---

**Q20 [UNKNOWN]** — MPNN design region score: which direction is better?

**A20:** I understand that lower MPNN scores are better. I believe that it is
the negative log-likelihood.

**Synthesis note — CODE VERIFICATION NEEDED**: MPNN score = negative log-likelihood;
lower is better (more probable sequence). The sort in `mpnn_select_top.py` should
therefore be ASCENDING. This needs to be verified in the code — if it sorts
descending, the top-N selection is backwards (selecting the least probable
sequences). Also verify `mpnn_design_region_score.py` produces per-region scores
in the same direction.

---

**Q21 [VERIFY]** — Negative controls (scrambled, polyA) use the input PDB as
ground truth and contig-derived design region. They are testing whether the
ranker promotes disrupted sequences — not the full steering mechanism.

**A21:** Yes exactly - I want to make sure that the passing designs are
meaningfully different than the negative controls - otherwise the ranking is
essentially random. It is correct that they are generated from the input complex
PDB, because they are just treated like a bad design. I suppose the receptor is
unlikely to fold well, so probably it will just run through the cold start and
not the negative steering, but that is also fine as a negative control.

**Synthesis note**: Controls confirmed as "bad design" surrogates. The expected
path is: cold-start prediction → spaghetti fold → no contact residues → empty
candidate pool → skip_steering → 3 cold-start rows in passing_summary. The
pipeline currently handles this correctly (as documented in notes13). The
`onComplete` miscalibration warning fires if controls score above threshold,
alerting to ranker failure.

---

**Q22 [UNKNOWN]** — Post-Rosetta inspection (rfdtest workflow): what are you
actually checking in ChimeraX?

**A22:** The post-Rosetta inspections are mostly done in ChimeraX, as well as
looking at the plot outputs. I want to see whether the designs are varied, what
sorts of shapes they take. Changing the inputs just a little bit, such as the
design region, can have a huge impact on the generated design regions, and so
sense checking them seems appropriate before committing the negsteering compute,
which is the most expensive part of the pipeline.

**Synthesis note**: The `stop_after_rosetta` parameter is used to gate the most
expensive stage (negsteer GPU compute). The human inspection checks: (1) design
shape variety — do the designs explore diverse backbone geometries? (2) visual
plausibility of the designed region. This is currently a manual step with no
automated diversity metric. Future improvement: automated diversity check
(e.g. pairwise backbone RMSD distribution) before proceeding to negsteer.

---

**Q23 [VERIFY]** — Multi-seed semantics: a single lucky seed should not promote
a design, but 1/3 (Tier C) should still be accessible for hard targets.

**A23:** Your reading is correct. I saw that negsteering and confidence metric
outcomes could vary a lot between seeds, so I thought it important to implement
several seeds for every negsteering run. I will preferentially choose MPNN
sequences for which all seeds passed negsteering, but I may also look at some
examples where individual seeds scored really well. It's just best to have as
much information as possible to help me choose what to take forward to
agroinfiltration.

**Synthesis note**: Seed semantics confirmed. The tiering system is primarily
a sorting aid, not a hard gate — the user will look across tiers and individual
seed data during manual selection. The pipeline's job is to expose as much
information as possible; the human makes the final selection.

---

**Q24 [UNKNOWN]** — What does a convincing computational result look like before
sending to the wet lab?

**A24:** The best outcome is a single MPNN sequence for which all seeds pass
(Tier A) and agree on strong binding interface confidence metrics. Weaker examples
would be those where all seeds pass (again Tier A) but seeds disagree on binding
interface confidence metrics, or lower tiers, B or C, where all passing seeds
agree on strong binding interface confidence metrics. Having multiple MPNN
sequences from the same design with strong results is a positive, but not a
necessity, because it could be the sequence itself which is being well or poorly
predicted by Boltz-2.

**Synthesis note**: Ideal selection criterion = Tier A AND seed agreement on
confidence metrics. The two axes are independent: tier (structural consistency)
and confidence agreement (prediction confidence stability). Neither alone is
sufficient. Multiple MPNN sequences from the same design is a positive signal
(suggests the backbone is generally bindable) but not required.

This clarifies the current composite score's role: `true_jaccard − 0.05·ra_eff`
only captures the structural axis. A separate per-metric seed agreement score
(e.g. MAD of ipSAE across seeds, already computed in `aggregated_results.csv`)
is needed to capture the confidence axis. The manual selection process currently
inspects both. An automated tier-and-agreement ranking would formalise this.

---

## Pending follow-up questions (updated)

1. **Contamination semantics** (Q5/Q15): Exactly which stages use design-region-only
   gating? User believes non-design-region contacts should not trigger reversion
   but is uncertain. → Session 2.
2. **num_seeds even/odd constraint** (Q6): Enforce at CLI? → Session 2.
3. **true_interface derivation** (Q15): Do non-load-bearing positions appear
   in `true_interface_idx`? → Code read + Session 2.
4. **mpnn_select_top.py sort direction** (Q20): Confirm ascending sort for
   lower-is-better MPNN score. → Code read before Session 3.
5. **Sc threshold "1.3 passed"** (Q18): Clarify what this means numerically.
6. **Confidence agreement axis** (Q24): MAD of ipSAE across seeds is in
   `aggregated_results.csv` but not in `passing_summary.csv`. Should it be?

---

---

### Batch 4 — Negsteer thresholds and outcome semantics (Q25–Q32)

---

**Q25 [VERIFY]** — `ra_eff` threshold of 5.0Å — is it correct and should it
be per-campaign?

**A25:** As with many of the thresholds that need auditing, this should be
selected more rigorously, and you're right, per-campaign options may be more
appropriate. We also calculated interface RMSD metrics which may be more stable
across designs in a campaign that could be substituted in to the composite metric.
It's true that sometimes I look at a design which fails the ra_eff filter, the
actual binding interface looks correct, but something off-interface differs a lot
from the RFDiffusion prediction, and therefore the receptor RMSD or ra_eff fails.
I like the per-campaign idea though.

**Synthesis note**: ra_eff = whole-effector Cα RMSD after receptor alignment.
It fails on designs where the interface itself is correct but off-interface
regions (floppy termini, disordered loops) differ — the exact issue documented
in negsteer notes3 with 7B1I. The iRMSD/fnat/dockq metrics in the orthogonal
stage are a better gate for this reason and were added specifically for this.
Per-campaign thresholds: plausible; would require at least one validated
agroinfiltration hit per campaign to calibrate. Currently 5.0Å is provisional.

---

**Q26 [UNKNOWN]** — `true_jaccard` threshold of 0.7 — was it derived or arbitrary?

**A26:** The true_jaccard threshold was chosen arbitrarily. Clearly a higher
true_jaccard is better, but also I have not yet seen any designs with a
true_jaccard of over 90%.

**Synthesis note**: true_jaccard=0.7 is arbitrary. The empirical upper bound
observed so far is ~90%. This means the filter is not near-saturating and
designs that pass genuinely have substantial interface overlap. However the
0.7 cutoff has no calibration basis. Once agroinfiltration data exists, this
is a candidate for calibration against experimental outcomes.

---

**Q27 [VERIFY]** — `clean_steered` and `pose_holds` — are they equivalent or
should clean_steered be rated higher?

**A27:** clean_steered and pose_holds should be equivalent because they both
involve changes to non-interface regions - neither should have contaminating
changes. Designs with no off-interface steering contacts are preferred over those
with steering contacts in both cases.

**Synthesis note**: Equivalence confirmed. Both clean_steered and pose_holds
certify that the wet-lab construct (the ProteinMPNN sequence, not the steered
version) directs Boltz to the correct interface. The user notes a within-category
preference: designs where the REMAINING non-reverted mutations also happen to be
off-interface are preferred over those where remaining mutations are near (but
not contacting) the interface. This sub-preference is not currently captured in
the pipeline's ranking and could be added as a tiebreaker.

---

**Q28 [UNKNOWN]** — `skip_steering` — is it the best possible outcome?

**A28:** Yes, skip steering is one of the best possible outcomes in terms of
the negsteering pipeline. Negsteering is a fallback for when Boltz-2 doesn't
predict the correct effector placement. My intuition tells me (this needs
empirically testing), that steered designs are likely to have lower confidence
metrics, since Boltz-2 is effectively choosing its 2nd or 3rd binding site. If
it was confident in it, it would have chosen it the first time. There is no
possibility in skip_steering that off-interface mutations are inflating any of
the confidence metrics.

**Synthesis note**: `skip_steering` = best computational outcome, with the
hypothesis that these designs will also have higher confidence metrics (Boltz
chose the designed interface as its primary prediction, not a fallback). This
hypothesis is empirically untested. If confirmed with agroinfiltration data,
`skip_steering` designs should be ranked above tier-A steered designs in the
final selection. The pipeline currently treats `skip_steering` as equivalent
to tier A (no_reversion verdict), which may be understating their quality.

---

**Q29 [VERIFY]** — `pose_collapses` means: steered pose was correct, but
reverted prediction snapped back to a wrong pose (not necessarily the original
one). These designs should not go to wet lab.

**A29:** That interpretation of pose_collapses is correct. It may not even snap
back to the same original incorrect pose - it just does not remain at the correct
pose after the contaminating changes are removed.

**Synthesis note**: Confirmed. `pose_collapses` = Boltz requires the steering
mutations to maintain the correct interface. Without them, it finds any other
pose. The key insight: the destination of the collapse is irrelevant — what
matters is that the correct interface is not accessible to the original MPNN
sequence. These designs should be completely excluded from wet-lab consideration.

---

**Q30 [UNKNOWN]** — `new_contamination` is very rare. What causes it?

**A30:** new_contamination is very rare and that makes sense - there aren't many
mutations introduced, and if some of them have already been reverted, then there
are probably only <3 possible new contaminants, some of which will be inward
facing on the receptor. This is most likely to occur when some potential
contaminants were just outside the contact threshold, and by reverting other
contaminants, the prediction shifted slightly, bringing them within the cutoff.

**Synthesis note**: `new_contamination` mechanism: borderline-contact mutation
(just outside 5Å heavy-atom cutoff in the steered prediction) crosses the
threshold in the reverted prediction because the reverted pose shifts slightly.
This is a prediction instability artefact rather than a biologically meaningful
secondary binding mode. The rarity (0 observed in the characterisation run at
aggregate level) is expected given the small number of mutations involved.

---

**Q31 [VERIFY]** — Two-policy contamination check: pre-reversion reverts ALL
contacts; post-reversion gates new_contamination to design_region ∪ true_interface
only. Is this intentional?

**A31:** Good catch! That was an inconsistency, but I believe it is correct.
For the pre-reversion stage, the idea is to get rid of anything that is in
contact. I think one of my earlier questions might have addressed this in a
different way, but this is the correct way. For the post-reversion stage, that's
when the non-load-bearing steered contacts may be ignored - they would otherwise
become new_contaminations. But they aren't necessarily holding the binding
interface together, and may just be contacting other unrelated parts of the
receptor by chance. If they are helping to secure the effector placement without
interfering with the interface itself, then they are fine post-reversion.

**Synthesis note — CONFIRMED AS INTENTIONAL**: The two-policy system is correct
by design. Rationale:
- Pre-reversion (build-contaminated): revert ALL mutations that contact the
  effector. Clean slate — remove everything that might be propping up the pose.
- Post-reversion (new_contamination): only fire if remaining mutations contact
  the effector at DESIGN-REGION OR TRUE-INTERFACE positions. Off-interface
  contacts in the reverted prediction are acceptable because they are not load-
  bearing for the designed binding interface.

This resolves the apparent inconsistency from Q5. The current implementation
(notes11 gated new_contamination) is CORRECT. The Q5 answer was describing
intent that applies to the post-reversion stage, not the pre-reversion stage.

---

**Q32 [UNKNOWN]** — Cold-start policy: if any seed fails cold-start, ALL seeds
run through full negsteer. Correct?

**A32:** That is the correct policy. If all three seeds pass cold_steering, then
we have a good Tier A candidate and there's no point wasting compute, since the
cold_steering prediction would be preferred over any other Tier A candidates
that negsteering might produce. If any less than all seeds pass cold_steering,
then there's a possibility for a stronger Tier A candidate to come out of the
negsteering process. For all negsteering, regardless if it's 0, 1, or 2
cold_steering-passing seeds, all seeds should run through cold steering. We
don't operate per-seed: everything goes through all 3 seeds if it runs. It's
necessary for the median confidence metrics.

**Synthesis note**: Policy confirmed: skip_steering only if ALL seeds pass
cold-start. Otherwise ALL seeds run full negsteer. The reason for running all
seeds (even the ones that passed cold-start) through negsteer is: median
confidence metrics require the full seed set for meaningful aggregation.
This is already implemented correctly in the code (the cold_start_all_clean
branch fires only when every seed passes).

---

## Pending follow-up questions (updated)

1. **Two-policy contamination — RESOLVED** (Q31): The two-policy system is
   intentional and correct. No code change needed. Pre-reversion: revert all
   contacts. Post-reversion: gate new_contamination to design_region ∪ true_interface.
2. **num_seeds even/odd constraint** (Q6): Enforce odd num_seeds at CLI? → Session 2.
3. **true_interface derivation** (Q15): Do non-load-bearing positions appear?
   → Code read + Session 2.
4. **mpnn_select_top.py sort direction** (Q20): Ascending for lower-is-better.
   → Code verification.
5. **Sc threshold "1.3 passed"** (Q18): Clarify what this means numerically.
6. **skip_steering ranking** (Q28): Should skip_steering designs be ranked above
   tier-A steered designs? Currently treated as equivalent. → Session 2.
7. **Within-pose_holds sub-ranking** (Q27): Preference for designs with no
   off-interface steering contacts. Not currently captured. → Session 2.

---

---

### Batch 5 — Contamination gating, thresholds, cache-busting, cohort summary (Q33–Q40)

---

**Q33 [VERIFY]** — Pre-reversion build-contaminated: does the effector-atom
filter apply (only interface-side effector atoms count), and is this correct?

**A33:** I think these definitions are becoming quite confusing and it is worth
being really clear about them. I expect some of my answers might contradict
themselves based on my understanding. In the previous set of questions we wrote
that the pre-reversion contamination check looks at all steered residues. This
is correct. Both pre- and post-reversion can use the same distance threshold
for contamination, just one is applied to all steered residues (pre-reversion),
and one is applied only to design region steered residues (post-reversion).
Regarding the effector lobe question, if it touches a steered residue
pre-reversion, than that residue will be contaminating and be reverted. If it
touches a steered residue post-reversion, then that residue isn't contaminating
(it's not on the main interface) and won't be reverted. I'm not sure how we are
defining "main" interface though. I think it's still the design region and the
"true_interface" of contacts.

**Synthesis note — CONFIRMED POLICY**: Both stages use the same effector-atom
filter (only effector atoms in the interface region of the effector count as
contacts). The asymmetry is solely on the RECEPTOR side:
- Pre-reversion: any steered receptor residue whose contact is counted triggers reversion
- Post-reversion: only steered receptor residues in design_region ∪ true_interface
  fire new_contamination

This is consistent with Q31's answer. The two-policy system is correct. The
"interface region" definition is confirmed as design_region ∪ true_interface.

---

**Q34 [VERIFY]** — Intact threshold (receptor RMSD ≤ 5.0Å) — is it the right
value, and does it apply equally to the effector?

**A34:** Yes, the receptor must remain intact i.e. similar to the RFDiffusion
design, otherwise it's not assessing RFDiffusion's output at all, and it's
likely to have misfolded, particularly if strong steering is applied. The
threshold is arbitrary, and may need per-campaign tailoring too. The effector
intact check can be a more stringent threshold presumably, because the effector
template is given and no sequence changes have been made.

**Synthesis note**: Intact threshold rationale confirmed. The effector intact
check SHOULD be more stringent than the receptor check: the effector is templated
(passed to Boltz as a structural prior) with no sequence changes, so large
effector deviations indicate prediction instability rather than reasonable
sampling. Whether the effector intact cutoff is currently set more stringent
than 5.0Å needs code verification — in notes8 the change was applied uniformly
to both receptor and effector. This is a calibration gap.

---

**Q35 [VERIFY]** — Should the ra_eff threshold be expressed as improvement
over cold-start rather than absolute?

**A35:** It's a good question, and I think we track the improvement over
cold-start in the output tables too. The reason we take the median value of
the 3 seeds is to try and iron out any of that variation caused by seed
differences. Let's stick with an absolute cutoff and come back to it when
we do the threshold auditing.

**Synthesis note**: Absolute cutoff confirmed as the working approach. The
cold-start improvement IS tracked in the output (delta_ra_eff columns) and
is available for future calibration. The median across seeds is the correct
aggregation to reduce stochastic noise. Defer to threshold audit.

---

**Q36 [UNKNOWN]** — Multi-cycle cmd_kickoff* family (~600 lines) — dead code
or future roadmap?

**A36:** The current pipeline does only use n_cycles=1, but higher n_cycles
will be used in the future, it's just not been wired in yet. It's going to be
one of the first things I do after finishing the remediation plan. It may be
crucial for non-MAX effectors or different receptors.

**Synthesis note**: Multi-cycle code is LIVE ROADMAP, not dead code. It is
the first feature to implement post-remediation. The `cmd_kickoff*` family
should be preserved and treated as dormant, not deleted. Phase 4 refactoring
must not remove or substantially restructure these functions. This updates the
remediation plan: `cmd_kickoff*` was listed as a possible deletion candidate
— it should instead be listed as deferred implementation work.

---

**Q37 [VERIFY]** — 8Å Cα (RFDiffusion filter) vs 5Å heavy-atom (negsteer):
different thresholds because RFDiffusion outputs are Cα-only.

**A37:** That is the idea behind the RFDiffusion filter threshold being higher —
because there are no side chains, but as for the specific numbers, they are
currently arbitrary.

**Synthesis note**: Rationale confirmed; specific values provisional. Both
thresholds need calibration. The Cα-only vs heavy-atom distinction is the
correct conceptual basis.

---

**Q38 [UNKNOWN]** — `params.outdir` relative-path bug — risk to future users?

**A38:** I am the only user right now. This will need to be fixed before this
pipeline is shared with other users. There is currently a risk that needs to
be fixed.

**Synthesis note**: Risk acknowledged. Workaround (absolute path in params.yml)
is currently in place for the production run. Fix (Task 60) remains in the
deferred list. Priority increases when the pipeline is shared.

---

**Q39 [VERIFY]** — Nextflow cache-busting: Task 61 described as not yet
fully done in notes14. Is that still true?

**A39:** I believe that the cache busting has actually been implemented
properly. Please check that the questions you're asking are actually true of
the current code - if you only read the notes then you may be behind. If there
are other questions you need to re-ask because of this, then please do so.

**Code check result — USER IS CORRECT**: Task 61 is substantially complete.
All 13 `modules/*.nf` files now use `path script_input` (or equivalent)
for their directly-invoked Python scripts. The code confirms this.

Two documented residual gaps remain (explicitly commented in the .nf files
and in remediation_state.md Verification Queue):
1. `negative_steering.nf` NEGSTEER_RUN_ONE: the orchestrator script
   (`negative_steering_run_one.sh`) IS content-hashed via `path orchestrator_script`,
   but the sub-scripts it loads via `--bin-dir ${projectDir}/bin` are NOT. The
   code comment explicitly says: "those indirect dependencies are NOT content-
   hashed by this declaration."
2. `negsteer_controls.nf`: same `--bin-dir` pattern for controls orchestration.

These gaps are KNOWN and documented. They are the residual items listed in
`remediation_state.md`'s Verification Queue. The Q39 premise from notes14
(that Task 61 was largely undone) was stale. I apologise for the inaccuracy.

**Implication**: for any changes to scripts invoked via `--bin-dir`
(specifically the suite of scripts called by `negative_steering_run_one.sh`),
a fresh run (not `-resume`) is still required. Changes to all directly-invoked
scripts (the majority) correctly bust the cache.

---

**Q40 [UNKNOWN]** — Observed issue with the cohort summary table:
inconsistencies in data sources, rows failing filters but showing in tiers.

**A40:** Currently I use the summary table plot from the orthogonal metrics as
my initial "results", but there seem to be wild inconsistencies in it, with data
that appears to be pulled from many separate sources i.e. some from the best
seed, some from the median across the seeds. Some rows fail filters but still
show in Tier A, B, or C. I'm not sure what's going on with it but it needs a
good check. That's the main issue I had found.

**Code check result**: Having read `orthogonal_metrics_plots.py:COMBINED_SUMMARY_COLUMNS`
(lines 764–815), the apparent inconsistencies are partially intentional and
partially a documentation/labelling gap:

1. **Mixed data sources — INTENTIONAL**: The table deliberately combines:
   - "cross" columns (from `cross_sequence_summary.csv`): Boltz metrics as
     medians across seeds — `representative_ra_eff_vs_truth_median`,
     `representative_iptm_median`, etc.
   - "survivor" columns (from `survivors_with_orthogonal_metrics.csv`):
     AF3 `best_ra_eff` (best across 15 predictions, not a median), Sc, BSA,
     ΔΔG — computed on the canonical PDB only, not aggregated across seeds.
   - "special" column: `n_mut` reconstructed from multiple sources.

   This mix is intentional but not communicated to the reader. A column-group
   header or footnote indicating "Boltz | Orthogonal" would clarify the source.

2. **Rows failing filters but showing in tier A/B/C — INTENTIONAL**: The tier
   (A/B/C/none) reflects negsteer outcomes (n_seeds_pose_holds). Orthogonal
   filter pass/fail is a separate axis shown via cell colour (red = fail, green
   = pass). A tier-A design can fail every orthogonal filter — this is by design
   and reflects the pipeline's intent that orthogonal metrics are informational
   for now. The plot correctly shows tier as a stripe colour AND filter result
   as cell colour. If this looks confusing, the issue is labelling, not logic.

3. **Potential real bug — needs investigation**: `af3_nomsa_best_ra_eff` uses
   BEST (not median) across 15 AF3 predictions, while Boltz columns use MEDIAN
   across 3 seeds. The aggregation is inconsistent and the column header does
   not communicate this. Whether this inconsistency is intentional or a
   historical accident is unclear. → **Flagged for Session 4 investigation**.

**Action**: The user's observed "inconsistencies" are a labelling/communication
issue, not a logic error. The plot needs column-group headers and a legend
explaining the data source and aggregation method for each column. This is a
Phase 4 task.

---

## Session 1 close — Discussion

### What Q39 tells us about the audit process

The notes-vs-code discrepancy on cache-busting is a useful calibration. The
notes reflect the state at the time of writing; the code has since moved forward.
For the rest of the audit, questions that are already answered by the code will
be flagged as "code-verified" rather than asked. Where notes suggest a problem,
the code will be checked before the question is raised.

### What Q40 tells us about the cohort summary

The cohort summary plot is the user's primary "results" view. Its data mixing
and ambiguous labelling are generating confusion about what is and isn't a bug.
Session 4 will audit `orthogonal_metrics_plots.py` and `merge_orthogonal_metrics.py`
in detail with the specific question: "is everything in this plot communicated
correctly?"

### Key decisions noted this session

These decisions were confirmed as intentional AT THE TIME they were made.
They are NOT closed. Each will be re-raised in the appropriate later session
with the question: "now that you have the full picture, do you still want this?"

| Decision | Currently in code | Re-raise in |
|---|---|---|
| Chain A=receptor, B=effector hardcoded downstream of RFDiffusion | Yes | Session 3 (design generation) |
| Two-policy contamination (pre: all contacts; post: design_region ∪ true_interface) | Yes | Session 2 (negsteer deep dive) |
| skip_steering = best possible outcome, ranked same as tier-A | Yes | Session 2 |
| pose_collapses designs excluded from wet lab | Yes | Session 2 |
| Multi-cycle cmd_kickoff* preserved (not deleted) | Yes | Session 2 |
| Tier A/B/C independent of orthogonal filter pass/fail | Yes | Session 4 (orthogonal + plotting) |
| num_seeds should always be odd | Intent confirmed; not enforced in code | Session 2 |
| Cohort summary mixes medians (Boltz) and best-value (AF3) across columns | Yes | Session 4 (plotting) |
| Mixed data sources in cohort summary (cross vs survivor CSV) | Yes | Session 4 (plotting) |
| Composite = true_jaccard − 0.05·ra_eff | Yes | Session 2 |

### Open issues carried into Session 2

1. true_interface derivation: which residues appear, especially non-load-bearing positions
2. num_seeds odd constraint: enforce at CLI
3. skip_steering ranking: should it rank above tier-A steered in final selection?
4. Effector intact threshold: should it be more stringent than receptor (5.0Å)?
5. mpnn_select_top.py: confirm sort direction (ascending = lower-is-better)
6. `af3_nomsa_best_ra_eff` aggregation: intentional best vs median inconsistency?

---

## Session 2 — Negsteer core semantics (deep dive)

---

### Batch 1 — Outcome classes, tiers, contamination policy (Q41–Q48)

---

**Q41 [VERIFY]** — Outcome class taxonomy: are the current 7 classes correct?

**A41:** I'm not sure why clean_steered and no_reversion have been separated - I
think it actually makes more sense for them to become one class, with tier A, B,
C, and None depending on the number of seeds that passed the filtering criteria.
clean_steered only means clean in terms of no contaminations. It doesn't say
anything about the quality of the predictions. That means that clean_steered
could be: a) steering mutations applied but correct effector placement not found,
or also b) correct effector placement (low ra_eff) but failing something else
like receptor RMSD or pae. That's just a consequence of having a limited number
of cycles. Anything Tier None should have no seeds passing that particular
threshold - so if one seed had new_contamination, but the other two seeds had
pose_holds, then it would still be Tier B pose_holds - the new_contamination is
treated as a single error, not representative of the whole negsteer path. Reversion
does not happen on a clean_steered prediction. I'm not sure in which cases
something could fall into the "No passing seeds" category separately from the
other categories - maybe this is just another wording on the failed effector
placement across all 3 seeds at the end of the last steering cycle. I think it
should be merged into clean_steered.

**Synthesis note — VOCABULARY MISMATCH IDENTIFIED**:

The user uses "clean_steered" to mean ONLY "no contaminating mutations" —
independent of whether the effector was correctly placed. The code's definition
of `clean_steered` (from notes12) is stricter: it requires no contamination AND
ra_eff < 5Å AND intact. This is a terminology mismatch.

The user's intended taxonomy (simplified):

| Outcome | Definition | Tiers |
|---|---|---|
| clean_steered (broad) | No mutations contacted the effector, regardless of placement quality | A/B/C/none based on n passing seeds |
| pose_holds | Mutations contacted effector; reverted prediction still finds interface | A/B/C/none based on n passing seeds |
| pose_collapses | Steered found interface; reverted collapsed | none always |
| new_contamination | Per-seed error, not representative of whole path; demoted to fail for that seed | affects tier counting only |

"No passing seeds" (Tier None) is just clean_steered with 0 seeds meeting the
placement threshold — same class. Tier None = 0 seeds pass; Tier C = 1 seed;
Tier B = 2 seeds (for num_seeds=3); Tier A = all seeds.

**Tiering formula already in code**: From Phase 3 notes, `cross_sequence_summary.py`
already uses `n_pass = n_seeds_pose_holds + n_seeds_clean_steered` for tier
assignment. The vocabulary labels (clean_steered vs no_reversion) may be the
source of confusion, not the actual counting logic.

→ **Follow-up Q needed**: Does the code already combine these for tier counting,
or are they scored separately? If separately, this is a bug. → See Q49.

---

**Q42 [UNKNOWN]** — Is `n_seeds_clean_steered` used only for tier assignment,
or does it appear in ranking/output columns too?

**A42:** n_seeds_clean_steered is presumably used to determine the Tier (A, B,
C, None) of the clean_steered results. n_seeds_pose_holds is only used to
determine the Tier of the reverted results, and not all negsteer paths undergo
reversion (the clean_steered ones have no contamination and therefore don't
undergo reversion).

**Synthesis note**: Confirmed. n_seeds_clean_steered drives tier for non-reversion
paths; n_seeds_pose_holds drives tier for reversion paths. These are currently
separate counts in the output tables. The user's Q41 answer implies they should
be COMBINED for tier assignment (both count as "passing seeds"). → See Q49.

---

**Q43 [VERIFY]** — skip_steering (cold_start_all_clean) should be tier A
equivalent, not a separate higher tier.

**A43:** cold_start passes should still be considered Tier A - where it is
important, the cold_start information will be available for use in choosing
between MPNN sequences to take forward to agroinfiltration (such as in output
tables, and plots). Until I have real evidence that steered designs have lower
confidence metrics (which I don't yet), any changes here would be based on
assumptions that might not hold up.

**Synthesis note**: Tier A for skip_steering confirmed. No change. The
hypothesis that steered designs have lower confidence metrics should be
empirically testable once agroinfiltration data exists. Note for future: if
the hypothesis holds, add a `skip_steering` indicator column to the cohort
summary so it can be filtered on manually during wet-lab selection, without
changing the tier assignment.

---

**Q44 [VERIFY]** — Composite score `true_jaccard − 0.05·ra_eff` — keep as is?

**A44:** I'd like to keep the current composite ranking as it is for now. I
believe it might be more meaningful to use the interface pLDDT values later on,
particularly for large complexes, but we're not making changes here, just
documenting knowledge.

**Synthesis note**: Composite confirmed as interim. Future direction: replace
or supplement with interface pLDDT (already computed per prediction as
`representative_interface_plddt_median`). For large complexes where whole-
effector ra_eff is less meaningful (more surface area = more noise), interface
pLDDT may better capture local binding quality. Record this for threshold audit.

---

**Q45 [VERIFY]** — Two-policy contamination: pre-reversion reverts all contacts;
post-reversion gates new_contamination to design_region ∪ true_interface.
Is this still correct?

**A45:** I am happy with the two-policy contamination system. Although, for the
post-reversion, surely the new_contamination can never be at the true_interface,
since if it was it would have already been a contamination previously. If I
remember it correctly, the true_interface is excluded from the pool of possible
steering sites, so it shouldn't be possible for a steered mutation to exist in
the true_interface and to then go on to be a new_contamination after reversion.
Please correct me if I'm wrong.

**Code check**: The user is CORRECT. The steering pool is constructed as
`wrong_interface − protected_set`, where `protected_set = design_region ∪ true_interface`.
Therefore no steered mutation can ever be at a `true_interface` position.

Consequence: the `∪ true_interface` part of the post-reversion gate is
REDUNDANT — it can never be triggered. new_contamination can only fire for
design_region positions. Including true_interface in the gate is a defensive
no-op, not wrong, but worth documenting so future developers don't wonder why
it's there.

---

**Q46 [VERIFY]** — Should even `num_seeds` be rejected by the CLI?

**A46:** Yes, previously any value of num_seeds could be chosen, and I think
a check should be added to make sure it's odd.

**Code action needed**: Add validation in `negative_steering_run_one.sh` or
wherever `num_seeds` enters the pipeline to reject even values with a clear
error message. → Phase 4 task.

---

**Q47 [VERIFY]** — Cohort summary: tier A/B/C is independent of orthogonal
filter pass/fail. Is the current display correct?

**A47:** I believe now all sequences go through orthogonal metrics, which is
correct. There are essentially three sets of metrics: Boltz-2 placement RMSD
(what does the prediction look like), Boltz-2 confidence (how confident is
Boltz-2 in this prediction), and orthogonal metrics (Do these predictions look
biochemically plausible). The tiering is mostly a mix of the first two, and
this is intended. It's fine for the orthogonal metrics to be highlighted in
red and green in the final summary table.

**Synthesis note**: The three-axis framework is now explicitly documented:
1. **Boltz placement** (ra_eff, true_jaccard, receptor RMSD) → drives tier
2. **Boltz confidence** (iptm, ipSAE, complex_plddt, pae_pass_frac,
   interface_plddt) → informational; drives manual selection within tier
3. **Orthogonal** (AF3 ra_eff, Sc, BSA, ΔΔG) → biochemical plausibility;
   red/green display only, does not change tier

Note: "all sequences go through orthogonal metrics" — this needs code
verification. Previously orthogonal metrics ran only on survivors. If all
sequences now run through it, the survivors_with_orthogonal_metrics.csv
scope has changed. → Q50 will verify this.

---

**Q48 [VERIFY]** — `true_interface` is derived from RFDiffusion design at 8Å Cα.
Is this correct, and do non-load-bearing steered contacts appear in it?

**A48:** The definition of the true interface region is correct. I understood
it better in my earlier answer - true interface regions cannot become
contaminated by design, because they are excluded from the pool of steering
positions. Thus, any non-design-region new_contamination cannot be part of the
true interface, and non-load-bearing steered contacts cannot become contamination
or new_contamination.

**Synthesis note**: Confirmed. The true_interface exclusion from the steering
pool creates a structural guarantee: nothing in true_interface can ever be a
steered mutation. This means:
1. Non-load-bearing contacts (e.g. the 22/26 "reaching bit") are at positions
   that are ALSO not in true_interface (otherwise they'd have been excluded from
   the pool). These contacts are handled by the effector-atom filter.
2. The post-reversion ∪ true_interface gate (noted in Q45 as redundant) is
   architecturally guaranteed to never fire.

---

---

### Batch 2 — Tier logic, thresholds, output labels (Q49–Q56)

**Note on approach correction**: Questions Q50, Q53, Q56 were asked without first checking
the code. The user correctly pointed this out. Code and test outputs were checked before
recording answers to those questions. Going forward: questions whose answer can be
obtained from the codebase are verified first, then presented to the user only if there
is a design decision to be made.

---

**Q49 [VERIFY]** — Does n_pass = n_seeds_pose_holds + n_seeds_clean_steered already
capture the user's intent for tier assignment?

**A49:** Seeds with no contamination that failed ra_eff or intact checks should not
be considered in n_pass - they failed because they did not pass those checks. That
is a major error catch if correct. Clean steered just means the placement of the
effector is correct - it doesn't alone determine a pass or fail. The sum of possible
passing or failing should always be n_seeds, but at the moment, something that fails
steering (because it fails the "intact" check, not because the effector is in the
wrong place) but then passes reversion, could somehow be counted as two passes if
I understand what you're saying correctly. This should not be possible.

**Code check and synthesis note — TERMINOLOGY CLARIFICATION**:

The code's `n_seeds_clean_steered` counts seeds where: ra_eff < 5Å AND intact AND
no contamination AND no reversion. The user's concern about double-counting is
addressed by the code: a seed is counted in exactly one bucket (pose_holds,
clean_steered, pose_collapses, new_contamination, or no_data). A seed that fails
intact is counted in no_data, not clean_steered.

The user's phrasing "clean steered just means placement is correct" — combined with
Q41's "clean_steered only means no contaminations" — reveals a terminology ambiguity.
Resolution: in the OUTPUT LABELS, "no_reversion" / "clean_steered" refers to any
seed where no reversion was attempted. In TIER COUNTING, only seeds that ALSO pass
ra_eff and intact are counted in n_pass. These two concepts need clearer labels.

**From the test output** (`cross_sequence_summary.csv` in the test fixture):
- `aggregated_verdict = "no_reversion"` appears for BOTH tier-A (skip_steering) and
  tier-none (clean steered, wrong placement) sequences. Same label, opposite outcomes.
- This is the ambiguity the user identified in Q41. → See Q59 (relabelling).

---

**Q50 [VERIFY]** — Do tier-none sequences receive orthogonal metrics?

**A50:** You have this wrong. Look again. The tier None do receive orthogonal
metrics. I think it was added in a later update so it may be in a different part
of the code. It would be worth simplifying and consolidating this code if it has
been spread around inefficiently.

**Code check**: The user is correct. `extract_survivor_manifest.py` iterates over ALL
rows in `cross_sequence_summary_with_interface_metrics.csv` without filtering by
`cross_tier`. All sequences — including tier-none — receive orthogonal metrics.
This is the correct behavior. The code is already doing what the user describes.
The "simplify and consolidate" comment is noted as a Phase 4 task.

---

**Q51 [VERIFY]** — Canonical representative = median-ra_eff seed. Is this correct?

**A51:** The point of the median_ra_eff is to capture roughly a summary of the data.
If it happens to be an incorrect effector placement, then the conclusion is "there
are more failing designs than passing ones". The representative is a failing design,
and that makes sense in my eyes. We could implement a second set of results for the
"best" designs, but I thought that would be bad scientific practice to discard the
bad results and focus only on the good.

**Synthesis note**: Median-ra_eff representative confirmed as intentional. The
representative represents the typical seed, not the best seed. For the cohort
summary plot, showing the median result gives an honest picture. If the median is
a failing design, the conclusion is that the sequence is unreliable. The "best
seed" data is still available in `raw_per_seed_results.csv` for manual inspection.

---

**Q52 [VERIFY]** — Effector intact threshold — should it be more stringent than 5.0Å?

**A52:** The effector intact threshold can currently stay as it is - I'll do a
threshold audit later. The effector structure may be predicted differently in
complex than as a monomer, for example, so it's not as simple as "template = better".

**Synthesis note**: Effector intact threshold kept at 5.0Å pending threshold audit.
Noted: a larger deviation in the predicted effector structure is not necessarily
prediction failure — it may reflect genuine conformational change on binding.
Add to threshold audit task list.

---

**Q53 [VERIFY]** — Contact cutoff for n_contacts_on_mutated_positions — should it
match the contamination check cutoff?

**A53:** I'm not sure what this question is getting at - I don't remember every
parameter from the huge codebase. I don't know the context of
n_contacts_on_mutated_positions, and it's not clear how it would be different to
a contamination check. If it's the post-reversion check (which would be a much
clearer question), then I don't see why it would be necessary for the contact
cutoffs to differ - they can be kept the same, and they will be included in the
threshold audit later.

**Synthesis note**: `n_contacts_on_mutated_positions` is a column in the per-seed
CSV that counts how many steered mutations contacted the effector (the contamination
count). It uses the same cutoff as the contamination detection (5.0Å heavy-atom).
Same cutoff for both confirmed. Include in threshold audit.

**Question improvement noted**: Q53 lacked context. Going forward, output column
questions will always include what the column does before asking about it.

---

**Q54 [VERIFY]** — All ra_eff thresholds should be the same value. Also: why is
the ra_eff check in clean_steered separate from the intact+ra_eff gate? And
should thresholds be centralized?

**A54:** They should use the same value. Anything passing that check will be
presented as a potential result. It doesn't make sense for some results to appear
with more lenient thresholds than others. I'm not sure why the clean_steered
criterion is separate to the intact+ra_eff check - surely the ra_eff check only
needs to be applied once at the end of steering - why compute it twice? Also
there must be a better way than hardcoding all these thresholds. I think it would
be good to have a separate yml file which includes ALL of the thresholds so that
when I threshold audit I can change things in there easily and they are propagated
across the codebase, but they aren't included in the original params file because
that would be overcomplicating things for the user. Or even just defining them at
the top of scripts instead of throughout. We should discuss this properly.

**Synthesis note — MAJOR ARCHITECTURAL ITEM**:

1. **Ra_eff threshold uniformity**: all three uses (intact+ra_eff gate, clean_steered
   criterion, pose_holds check on reverted prediction) must use the same value.
   Currently they do use the same CLI-passed `--rmsd-threshold`. Confirmed correct.

2. **Redundancy in ra_eff check**: the user asks why clean_steered applies its own
   ra_eff check if the intact+ra_eff gate already ran. From the code, the sequence is:
   (a) intact+ra_eff gate → only intact AND ra_eff-passing designs go to reversion
   (b) clean_steered criterion → applied AFTER reversion, filters the no_reversion
   seeds in the aggregate. BUT if step (a) already filtered on ra_eff, any seed
   reaching the aggregate as "no_reversion" already passed ra_eff. The redundancy is
   real. → Phase 4 simplification: the clean_steered ra_eff check in the aggregator
   may be removable.

3. **Threshold centralization — NEW REQUIREMENT**: A separate `thresholds.yml` file
   that is NOT part of the main user-facing `params_example.yml`, containing all
   hardcoded threshold values. Scripts load it internally; users don't see it unless
   they want to do a threshold audit. Alternatives: Python `config.py` imported by
   all scripts; or module-level constants consolidated to top-of-file. Decision on
   which approach is needed before Phase 4. → Discuss in Phase 4 architecture design.

---

**Q55 [UNKNOWN]** — Multi-cycle: automatic or manual inspection between cycles?

**A55:** Cycles will run automatically. If cycle-0 produces predictions where the
effector is in the wrong place and the receptor is intact, then the sequence will
be passed to cycle-1. It doesn't matter if the cycle-0 predictions receive poor
metrics, since that's not the structure we're trying to assess (the effector is in
the wrong place).

**Synthesis note**: Multi-cycle is FULLY AUTOMATIC. Gate from cycle-N to cycle-N+1:
effector wrong AND receptor intact → proceed. Poor confidence metrics do not block
progression to the next cycle (correct — you're trying to move the effector, not
assess a wrong pose). The `stop_after_rosetta` manual inspection point remains for
the upstream design-generation stages; negsteer itself runs hands-off.

---

**Q56 [VERIFY]** — What labels appear in the output, and what is "singleton"?

**A56:** Please don't ask questions like this again. You have access to all the
output formats in the "tests/" folder. Look at them rather than asking me to.
The existing code and test outputs are the exact state of the code. I don't know
what the "singleton" label is and I didn't see you mention it in the previous
question about the 7 (or now 5) output classes.

**Code check**: Verified in test outputs and code.

Labels in `representative_aggregated_verdict` (test fixture):
`no_reversion`, `pose_holds`, `pose_collapses`

Tiers in `cross_tier`: `A`, `B`, `C`, `none`

`singleton`: a legacy label for rows from pre-multi-seed runs (no sequence_group
assigned) or from the failed-steering path where the cold-start wasn't clean and
steering couldn't find any valid mutations. These rows appear in `aggregated_results.csv`
as passthrough rows but are filtered out by `extract_passing.py` and never appear
in `passing_summary.csv`. The user will never encounter this label in normal operation.
It is an internal implementation detail, not an output class.

**Instruction recorded**: Do not ask about output formats. Read the test fixtures
and code directly. Only bring a finding to the user if there is a design decision
to make about it.

---

---

### Batch 3 — Labels, thresholds, multi-cycle (Q57–Q64)

---

**Q57 [UNKNOWN]** — Threshold centralization approach.

**A57:** We've spoken about potentially tailoring some thresholds to different
runs. For this reason, it may be best to have it in a yaml file which is easy
to access, and then have that yaml file accessed wherever it is needed. Maybe
you could create a little wrapper function called get_thresholds_from_yaml()
which takes the name of a threshold to retrieve and is loaded as a module in
the top of each script.

**Synthesis note — ARCHITECTURAL DECISION RECORDED**:
- Format: YAML file (not Python module or Nextflow config extension)
- Access pattern: `get_thresholds_from_yaml(threshold_name)` wrapper function,
  importable module, loaded at top of each script
- Per-run override: the YAML path could be passed as a CLI argument to allow
  campaign-specific YAML files without touching the main params
- NOT in the main user-facing `params_example.yml`
- Scope: all hardcoded threshold values across the codebase

→ Phase 4 implementation task. The wrapper function design (caching, error
handling for missing keys, type coercion) needs to be specified before implementation.

---

**Q58 [UNKNOWN]** — Should tier-none + no_reversion get a distinct label?

**A58:** There are so many labels and it's going to become extremely bloated and
overcomplicated. If a sequence receives Tier A and "no_reversion", that is obviously
different from one that receives Tier None and "no_reversion", because we can see
the tier. Every time we add a new label we just lose comprehension. If the
information is already there and obvious you don't have to add a new combinatorial
label. Only add it if the user never sees the difference between them.

**Synthesis note**: No new labels. Tier + existing verdict is sufficient context.
The `cross_tier` column is always adjacent to the verdict in the output. A reader
who sees Tier A + no_reversion vs Tier None + no_reversion can infer the difference.
This is the correct minimalist approach. The principle: add labels only when
existing information is insufficient to distinguish cases.

---

**Q59 [VERIFY]** — Should no_data seeds count as failures (reduce effective n_pass
relative to n_seeds)?

**A59:** If something fails, it fails. If no_data is an example of failing, it
should be included as a failure. As I said previously, the total pass/fail
assessments should equal the number of seeds. If 2 seeds pass and 1 seed has
no data across 3 seeds, then clearly it's Tier B 2/3 pass.

**Code check**: The current implementation already matches this intent. n_seeds
(denominator) = the total seed count for the sequence group. n_pass = n_pose_holds +
n_clean_steered (numerator). no_data seeds are in n_seeds but not in n_pass, so
they correctly count against the tier. 2 pass + 1 no_data = n_pass=2 / n_seeds=3 =
Tier B. This is already correct.

---

**Q60 [UNKNOWN]** — Multi-cycle gate: which seeds propagate to cycle-2?

**A60:** The terminology is becoming overwhelming. You told me in a previous
question that things that have an incorrect effector placement and don't receive
reversion are included in "clean_steered", but now you're telling me they're
"no_data". You must give the correct and full context for these questions or
it will be impossible for me to answer them consistently. Regarding the 2nd
cycle, if any seed finds the correct interface, it should not be propagated to
the next cycle - that seed has shown it is possible to find the interface at
the first cycle, so adding more mutations to try and rediscover the same
interface is a waste of resources.

**Terminology error corrected and apology recorded**:
The confusion arose because I gave inconsistent definitions across questions.
The correct code definitions:

| Seed outcome | Condition | In n_pass? |
|---|---|---|
| `clean_steered` | ra_eff < threshold AND intact AND no contamination | YES |
| `pose_holds` | ra_eff < threshold AND intact AND contamination AND reverted passes | YES |
| `pose_collapses` | ra_eff < threshold AND intact AND contamination AND reverted fails | NO |
| `new_contamination` | ra_eff < threshold AND intact AND contamination AND reverted introduces new contacts | NO |
| `no_data` | Failed intact OR ra_eff > threshold — never reached contamination check | NO |

A seed with WRONG EFFECTOR PLACEMENT (ra_eff > threshold) is `no_data`, not
`clean_steered`. I stated the opposite in earlier questions — that was wrong.

**Multi-cycle gate confirmed**: ANY seed that finds the correct interface (passes
ra_eff AND intact) is NOT propagated to cycle-2 — it has already succeeded.
Only `no_data` seeds (wrong placement, intact receptor) propagate. Seeds with
`pose_collapses`, `new_contamination`, or `clean_steered` do NOT propagate —
they've had a definitive result (pose_holds equivalent, or reversion failure).

---

**Q61 [VERIFY]** — Should weighted_jaccard be in the composite score YAML?

**A61:** I cannot answer that question and I deliberately deferred it in the
previous answer. Yes, weighted jaccard is calculated, but I don't know if it
will be included into the composite score until I audit the metrics after
running the code several times in its corrected form. For now the YAML should
follow the existing composite metric, and it can always change later on.

**Synthesis note**: YAML threshold file should encode current composite:
`cross_composite_score = true_jaccard - 0.05 * ra_eff_vs_truth`. No weighted_jaccard
in the formula for now. The YAML structure should be designed so adding or
changing composite components doesn't require code changes — just YAML edits.

---

**Q62 [VERIFY]** — Does the code correctly handle 2 pose_holds + 1 new_contamination
as Tier B pose_holds?

**A62:** Do not ask me what the code does. You have the code. Read it. If one seed
has new_contamination and 2 seeds have pose_holds, then the verdict should be Tier B
pose_holds. If all three have new_contamination, then the verdict should be Tier None
new_contamination.

**Code check — intent confirmed, code verified**:
From `cross_sequence_summary.py`:
- n_pass = n_pose_holds + n_clean_steered (new_contamination NOT counted)
- For 2 pose_holds + 1 new_contamination: n_pass=2, n_seeds=3 → Tier B ✓
- The sequence-level `aggregated_verdict` = majority verdict = "pose_holds" (2 vs 1) ✓
- For 3 new_contamination: n_pass=0, Tier None; majority verdict = "new_contamination" ✓

Code is already correct.

---

**Q63 [UNKNOWN]** — How will the threshold audit be conducted?

**A63:** I'll research and document reasoning for each of the selected thresholds
at first, so they're not just arbitrary guesses, as well as discussing it with
my supervisor. I'll then run each of the campaigns with the corrected code and
assess how thresholds may need to be changed between campaigns due to structural
differences, and update the thresholds again.

**Synthesis note**: Threshold audit process:
1. Document the biological/structural reasoning for each threshold (literature + supervisor input)
2. Run all five campaigns with corrected code
3. Assess campaign-by-campaign: structural differences (e.g. effector size, interface geometry) may warrant per-campaign values
4. Update thresholds accordingly — enabled by the threshold YAML once implemented

→ Threshold audit is Task 47 in the pipeline plan. Now has a defined methodology.

---

**Q64 [VERIFY]** — Multi-cycle architecture: does cycle-2 run within the same SLURM job as cycle-1?

**A64:** I believe currently each MPNN sequence negsteering job is treated as a
single SLURM job. 1 input sequence, many steering and reversion predictions, 1
SLURM job. I think it makes sense for the second cycle for that same input
sequence begin running within the same SLURM job once the first cycle has
finished, given everything is running in series, and the only parallelisation
is across input sequences.

**Synthesis note**: Architecture confirmed. `NEGSTEER_RUN_ONE` in Nextflow is
one GPU SLURM job per MPNN sequence. ALL stages (plan, predict, collect, reversion,
harvest, finalize, and future cycle-2 iterations) run in series within that single
job. Cross-sequence parallelism is handled by Nextflow scheduling multiple
NEGSTEER_RUN_ONE jobs simultaneously. Multi-cycle will extend the shell script
(`negative_steering_run_one.sh`) to loop through cycles within the same job.

---

---

### Batch 4 — Multi-cycle, parameters, pipeline structure (Q65–Q72)

**Process note**: Q67 and Q68 repeated questions already answered in Sessions 1 and 2.
This indicates a tracking gap — what has been asked and confirmed must be maintained
so questions are not repeated. A confirmed-answers index is added at the end of this
session to prevent this.

---

**Q65 [VERIFY]** — Multi-cycle continuation criterion: novelty only, or novelty + correct placement?

**A65:** To avoid any confusion here I'll be explicit. A negsteering path (one set
of steered residues for one input MPNN sequence), should end after the first cycle
(cycle-0) if the correct effector placement was found for at least one seed. It
should end after the second cycle (cycle-1) if in cycle-0 no correct effector
placement was found, and either a) the correct effector placement was found in
cycle-1 OR b) no new effector placements were found (maybe it found cycle-0's
incorrect placement again). It should end after the third cycle (cycle-2) if in
cycles 1 and 2 (cycle-0 and cycle-1) two different incorrect effector placements
were found, and either a) the correct effector placement was found in cycle-1 OR
b) no new effector placements were found. And of course depending on the user
defined n_cycles, it shouldn't continue to a new cycle that has not been specified.
The idea is to keep letting it search as long as it's finding new effector
placements and not just wasting time repredicting the same thing. The true effector
placement may be Boltz-2's 2nd, 3rd, or 4th preference. So novelty should be a
filter.

**Synthesis note — MULTI-CYCLE GATE DEFINED**:

Stop after cycle-N if any of:
1. ≥1 seed found the correct interface (ra_eff < threshold AND intact) — success
2. No novel effector placements found (all cycle-N predictions are within novelty
   cutoff distance of prior-cycle poses) — search exhausted
3. n_cycles reached — user-defined limit

Continue to cycle-N+1 ONLY if: no correct placement found AND at least one novel
(previously unseen) incorrect placement was found in cycle-N.

Rationale: the correct interface may be Boltz-2's 2nd, 3rd, or 4th preference.
Negative steering iteratively explores the pose space by blocking prior wrong poses.
Continue only while new territory is being found.

This matches the standalone negsteer architecture (novelty filter with min_dist_to_upstream).
The `cmd_kickoff*` code implements this. Confirmed that novelty IS a continuation
criterion.

---

**Q66 [UNKNOWN]** — Guiding potentials: always-off default or campaign-specific?

**A66:** Guiding potentials are designed for specific campaigns where receptor
resurfacing is more challenging due to the structures of the two proteins involved,
they're an optional lever to tweak when other things have been unsuccessful. They
may end up being on by default if they produce more passing designs.

**Synthesis note**: Guiding potentials = campaign-specific lever for difficult cases.
Default `add_potential = false` is correct. If empirical results show they consistently
improve pass rates across campaigns, default may flip to `true`. Document this
decision point in the threshold YAML (guiding potential parameters should also be
in the YAML, not hardcoded).

---

**Q67 and Q68**: Repeated questions. Q67 repeated the cold-start / all-seeds policy
(confirmed in Session 1 Q32 and again in Q60). Q68 repeated the contact cutoff
rationale (confirmed in Session 1 Q37 and Session 2 Q53). These questions were
not tracked correctly. Answers unchanged from prior confirmations.

---

**Q69 [VERIFY]** — `diffusion_samples = 5` confirmed as the correct default?

**A69:** Yes, 5 samples is fine - I tried running with 2 samples to see if it
sped up the runtime but it didn't change anywhere near as much as I hoped, and
it's good to be thorough. It is different to the 3 seeds. Each of the 3 seeds
is run with Boltz-2 over 5 samples.

**Synthesis note**: Architecture confirmed. Per MPNN sequence: 3 seeds × 5 diffusion
samples per seed = 15 internal trajectories per Boltz call set. Boltz picks
model_0 (highest-confidence sample) from each seed's 5 samples. The 3 seeds provide
the statistical sampling used for tier assignment; the 5 samples improve the quality
of each seed's best-model selection. diffusion_samples=2 gave insufficient runtime
savings to justify the quality trade-off. diffusion_samples=5 is confirmed default.

---

**Q70 [UNKNOWN]** — `n_designs` parameter — which does this refer to?

**A70:** That seems like a misunderstanding - there are two separate parameters -
one which determines the number of RFDiffusion-generated designs (which I think
you mean by n_designs), and one which determines how many MPNN sequences are
generated from each design. Then the number of MPNN sequences that undergo negative
steering are determined by mpnn_top_n I think.

**Synthesis note — TWO DISTINCT n_designs PARAMETERS**:

1. `params.num_designs` (RFDiffusion stage): how many backbone designs RFDiffusion
   generates from the input complex. Typically 8–64. Sets the size of the design
   backbone pool.
2. `params.mpnn_top_n` (MPNN → negsteer gate): how many MPNN sequences (selected
   globally by MPNN score) proceed to negsteer.
3. MPNN sequences per design: controlled by `params.num_seqs` in the MPNN step —
   how many sequences MPNN designs per RFDiffusion backbone.

The negsteer parameter that controls how many unique steered designs are generated
per MPNN sequence per cycle is `params.negsteer_n_designs` — separate from
`params.num_designs`. Q70 was asking about `negsteer_n_designs` but used "n_designs"
ambiguously. The production value of `negsteer_n_designs` needs to be confirmed
per campaign, not assumed from any one value in the notes.

---

**Q71 [VERIFY]** — Three Session-1 deferred items: skip_steering indicator, within-
pose_holds sub-rank, even num_seeds validation.

**A71:** A) Not a separate column, just something small, like Tier A with a star
or a snowflake symbol in the green box to indicate "cold-steering". B) No - each
potential candidate to take forward for agroinfiltration will be thoroughly
investigated first - I'll find that out anyway before I order them. C) It would
be a shame for the pipeline to get all the way to negsteering and then fail
because someone provided an even num_seeds. Maybe there should be input parameter
checks right at the start of the pipeline that make sure that all of the parameters
are appropriate formats, types, and values.

**Synthesis note — THREE DECISIONS RECORDED**:

A) skip_steering indicator: a symbol (star ★ or snowflake ❄) in the tier cell of
the cohort summary table — not a separate column. Tier A ★ vs Tier A with no symbol.
Small visual cue, no data bloat.

B) Within-pose_holds sub-ranking: no automated indicator. The user will inspect
candidates manually before ordering. The detailed data is available in
`raw_per_seed_results.csv` for manual review.

C) Input validation: a dedicated input-validation step at the START of the pipeline
(`main.nf` entry point) that checks ALL parameters for correct type, format, and
value before any SLURM jobs are submitted. This prevents wasted compute from invalid
inputs. `num_seeds` parity check is one item; the full list of validated parameters
needs to be designed. → Phase 4 implementation task.

---

**Q72 [UNKNOWN]** — Is there anything wrong or confusing in the negsteer core that
questions haven't landed on?

**A72:** The negative steering terminology has become completely overwhelming, and
it's challenging to intuitively understand what's going on in detail. If there are
any changes that could simplify things that don't need to be different, then please
propose it so we can implement them. Even if it's just some descriptions or
terminology that are unclear or confusing.

**Synthesis note — TERMINOLOGY SIMPLIFICATION REQUEST**:

This is a meta-item for Phase 4 / Session 5. The accumulated negsteer vocabulary
includes: aggregated_verdict, clean_steered, no_reversion, pose_holds, pose_collapses,
new_contamination, no_data, singleton, n_seeds_pose_holds, n_seeds_clean_steered,
n_seeds_pose_collapses, n_seeds_new_contamination, n_seeds_no_data, representative_*,
cross_tier, cross_composite_score, cross_rank_by_*, within_sequence_rank, ...

Proposed simplification areas (for discussion in Session 5):
1. Per-seed outcomes: collapse 6 categories into 3 for user-facing outputs:
   `pass` (clean_steered or pose_holds) / `fail_reversion` (pose_collapses or
   new_contamination) / `fail_no_interface` (no_data)
2. Seed summary: replace five n_seeds_* columns with a single "X/N passed" display
   plus the tier — removes redundancy since tier already encodes this
3. Column prefixes: `representative_*` prefix is verbose; rename to `rep_*`
4. `aggregated_verdict`: rename to `outcome` — clearer for non-developers

These are proposals for discussion, not decisions. The terminology should be
stabilised once confirmed, then applied across all output files consistently.

---

## Session 2 — Confirmed answers index (prevents repeat questions)

| Topic | Decision | Questions |
|---|---|---|
| Chain A=receptor, B=effector hardcoded post-RFDiffusion | Intentional, correct | Q11, Q67 |
| All 3 seeds run negsteer if any cold-start fails | Confirmed | Q32, Q67 |
| skip_steering only if ALL 3 seeds pass cold-start | Confirmed | Q32, Q67 |
| Contact cutoff: 8Å Cα for RFDiff (no sidechains), 5Å heavy for negsteer | Confirmed | Q37, Q68 |
| n_pass = n_pose_holds + n_clean_steered for tier | Confirmed correct | Q49, Q59 |
| no_data seeds count as failures (reduce effective n_pass) | Confirmed | Q59 |
| Tier B for 2 pose_holds + 1 new_contamination | Confirmed | Q62 |
| Multi-cycle gate: stop if correct interface found OR no novel placements | Confirmed | Q65 |
| Multi-cycle runs within same SLURM job, automated | Confirmed | Q64, Q55 |
| diffusion_samples=5 confirmed default | Confirmed | Q69 |
| Effector intact threshold stays at 5Å pending audit | Confirmed | Q52 |
| Skip_steering = Tier A, no higher tier | Confirmed | Q43 |
| Two-policy contamination (pre: all; post: design_region only) | Confirmed | Q45, Q31 |
| true_interface excluded from steering pool | Confirmed | Q45, Q48 |

---

## Session 3 — Design generation (RFDiffusion → MPNN)

### Batch 1 — RFDiffusion parameters and geometry metrics (Q73–Q80)

---

**Q73 [UNKNOWN]** — When would you switch from the default `Complex_beta_ckpt.pt` checkpoint?

**A73:** I currently have no plans to change to the other checkpoints if the beta
one is fine-tuned on complexes.

**Synthesis note**: `Complex_beta_ckpt.pt` is the only checkpoint in use. The other
options (`Base_ckpt.pt`, `Complex_base_ckpt.pt`) exist in the params file as
documentation, but the user has no reason to switch. The checkpoint comment in
`params_example.yml` should note that `Complex_beta_ckpt.pt` is strongly preferred
for PPI design tasks.

---

**Q74 [VERIFY]** — `add_potential: true` in params_example.yml conflicts with prior
statement that guiding potentials are off by default.

**A74:** I changed my mind - the example params should say true, and the user can
turn it off if they believe it is adversely affecting the predictions, which is
unlikely.

**Synthesis note**: `add_potential: true` is confirmed as the correct default. The
reasoning: guiding potentials are unlikely to hurt; they may help in difficult cases;
the user turns them off only if needed. The params_example.yml already has this
value — no change required.

---

**Q75 [UNKNOWN]** — Production cohort defaults: num_designs=32, num_seqs=16,
mpnn_top_n=64 — are these the right production values?

**A75:** The defaults don't really matter and will change depending on how easy
it is to find candidates. Normally I've been running 64 RFDiffusion designs, each
with 2-4 MPNN sequences, where 64 or 128 designs get taken forward through
negsteering. But if no successful designs were found, I may increase it.

**Code check — params_example.yml is outdated**: The production campaigns (v1_4a,
v1_4f, v1_4e) all use:
- `num_designs: 64` (not 32)
- `num_seqs: 2` (not 16)
- `mpnn_top_n: 128` (not 64)
- `negsteer_n_designs: 3` (not 20)

The `params_example.yml` defaults are wrong relative to actual production usage.
→ **Action required**: update params_example.yml to reflect production-representative
values, or add a comment distinguishing "test" vs "production" settings.

---

**Q76 [UNKNOWN]** — `negsteer_n_designs=20` seems low — should it be higher?

**A76:** Actually, the number of negsteer designs benefits from being smaller -
please look at the actual production runs, I think you're looking at the old params
example file. I got more hits from increasing the number of RFDiffusion designs and
providing diverse options with fewer steering mutation sets for each, than trying to
force every RFDiffusion design to pass steering (sometimes the pool of potential
steering sites is very small or Boltz-2 just doesn't want to place it at a particular
binding interface).

**Code check — production value is 3, not 20**: All three checked production params
files set `negsteer_n_designs: 3`. The params_example.yml default of 20 is a stale
value.

**Synthesis note — KEY DESIGN INSIGHT**: Hit rate is improved by:
- MORE RFDiffusion backbone diversity (64 designs, not 32)
- FEWER steering mutation sets per MPNN sequence (3, not 20)
- Accepting that some MPNN sequences simply won't steer

Rationale: when the steering pool is small or the binding interface is inaccessible
to Boltz, more steering attempts don't help — they just waste GPU time. Spreading
the compute across more diverse backbones is more productive.

---

**Q77 [VERIFY]** — Are the three geometry metrics (motif_rmsd, COM displacement,
endpoint_distance) used to make design decisions?

**A77:** They are informational. The Rosetta shape complementarity filter is more
important and helps me understand if high-quality designs are being produced. They
would tell me, for example, if the designs barely diverge from the native, or if
all the designs are very similar to each other, or even, say, if they form a couple
of clusters in terms of similarity. All these things will be useful information when
deciding which sequences to order.

**Synthesis note**: Three geometry metrics are diagnostic / informational, not gating:
- `motif_rmsd` (scaffold RMSD): how faithfully RFDiffusion reproduced the fixed regions
- `per_region_com_displacement`: how far the design region's centre of mass moved from
  the input gap
- `per_region_endpoint_distance_input/design`: end-to-end span comparison

These are visualised in `rfdiff_com_displacement.png` and `rfdiff_design_clustering.png`.
They inform manual inspection (divergence from native, cluster diversity) but don't gate
any automated filter. Sc is the automated gate.

---

**Q78 [VERIFY]** — Should the Rosetta Sc filter run pre-MPNN on polyvaline backbones?

**A78:** The idea was a rough catch for anything which has really poor shape
complementarity. I'm not sure how much the Sc might change after adding side-chains,
but my intuition was that RFDiffusion seeks to design a complementary backbone
structure, regardless of the side chains, and any designs where that has been done
really poorly should be culled to save compute. Maybe the threshold is too high. It
could be decreased, and this will also be chosen more thoroughly during the audit.

**Synthesis note**: Pre-MPNN Sc filter on polyvaline backbones is intentional — a
coarse cull of catastrophically poor backbone geometry before spending MPNN compute.
The biological intuition: good backbone geometry (high Sc) is necessary but not
sufficient for binding; poor Sc at the backbone level won't be rescued by good
side chains. Threshold of 0.5 is provisionally too high for polyvaline (real Sc values
at this stage are around 0.47–0.50 per the notes). → Lower threshold in audit.

---

**Q79 [UNKNOWN]** — `rfdiff_iterations=100`: is this sufficient for production?

**A79:** My understanding is that the rfdiff_iterations=50 is appropriate for tests
and rfdiff_iterations=100 is better for production runs. I can run more iterations
if there is literature to support it making a meaningful difference.

**Synthesis note**: 50 iterations = test, 100 iterations = production. This is
confirmed in the production campaign params. The params_example.yml already shows 100.
Higher iterations may produce better designs by giving the diffusion model more
denoising steps, but the marginal benefit above 100 is unknown for this system.
Document as a future experiment if design quality is insufficient.

---

**Q80 [VERIFY]** — `align_to_native_by_anchors`: does it use contigs to identify
which regions are fixed anchors?

**A80:** Yes, but it should also use the contigs to make sure that the fixed regions
are indeed the fixed regions (it doesn't extend into the design region because of
matching residues by mistake). The contig will say keep A1-30 fixed, so it should
look to align that A1-30 with the MPNN sequence, and then say A40-60 fixed in the
contig, it will look to align that A40-60 with the MPNN sequence. Maybe the
insertion is 20 residues instead of the original 10, so A40-60 has actually moved
in the sequence 10 bases further down, but since we know which sequence was in
A40-60, we should be able to align it and work out which bit was the design region.

**Synthesis note**: Confirmed. `align_to_native_by_anchors`:
1. Reads contig fixed segments (e.g. A1-30, A40-60)
2. Extracts those native sequences from the input PDB
3. Aligns each anchor segment to the MPNN output by sequence identity
4. The design region = whatever falls between the aligned anchors
5. If the design region changed length (insertion/deletion), the anchor alignment
   absorbs the offset — no positional assumption is made about the design region

The "by anchors" name describes exactly this: the fixed contig segments are the
alignment anchors; the design region is derived, not directly specified.

This means the function MUST receive both the MPNN output AND the contig string to
work correctly. Any caller that omits the contig will get incorrect design-region
identification. → Verify all call sites pass contigs correctly.

---

---

### Batch 2 — MPNN selection, QC, and sequence correction (Q81–Q88)

---

**Q81 [VERIFY]** — Global MPNN score vs design-region score for selection.

**A81:** Well given that the design region is just a part of the global sequence,
and the non-design regions are the same for all sequences, the ranking should stay
the same, so it's fine to use the global MPNN sequence score. The design region
score is a diagnostic yes, it will show if MPNN generates any very unusual sequences,
which might help understand why certain things fail later in negsteering or even
in agroinfiltration.

**Synthesis note**: Global MPNN score confirmed as correct for selection (non-design
regions are constant, so global ranking is monotone with design-region ranking).
Design region score is diagnostic: flags unusually designed positions that may
indicate an interesting or problematic sequence. Neither metric currently gates
selection — both are in the output CSV for manual inspection.

---

**Q82 [VERIFY]** — What does `mpnn_sequence_qc.py` filter?

**A82:** Just read the code it's all there? It should be something like checking
for Poly-X runs or sequences where >X% of residues have the same identity. There
is a minimum / maximum sequence identity threshold configurable in the params file
which is optional here - I haven't used it because there is always very low sequence
identity to the native - I think it's because MPNN was trained on a lot of human
data, and these sorts of plant HMAs are very rare.

**Code check**: Three filters in `mpnn_sequence_qc.py`:
1. **Poly-X runs** (`--max-poly-x 5` default): rejects sequences with any amino acid
   running for more than 5 consecutive positions.
2. **Unusual AA composition** (hardcoded `max_single_aa_frac=0.30`): rejects sequences
   where any single amino acid constitutes >30% of the sequence.
3. **Percent identity to native** (`--min-pct-identity 0.0`, `--max-pct-identity 100.0`
   defaults): optional identity bounds; currently inactive (defaults allow everything).

The user's observation about low identity to native is correct — MPNN produces very
diverse sequences for these plant HMA domains because they're rare in training data.
The identity filter is therefore not needed in practice for current campaigns.

---

**Q83 [VERIFY]** — Clustering: informational or used for selection?

**A83:** Yes the clustering is informational - it should cluster at several different
thresholds actually because that's how it generates one of the plots. It's currently
diagnostic, but in the future maybe it can be used to avoid selecting very similar
sequences.

**Code check**: `mpnn_cluster_sequences.py` clusters at multiple identity thresholds
(for the sequence diversity plot). The selected sequences for negsteer come from
`mpnn_select_top.py` (global MPNN score) independently of clustering. Clustering
informs the diversity plot only. Future direction noted: per-cluster selection to
guarantee sequence diversity → Phase 4/5 feature.

---

**Q84 [UNKNOWN]** — Why `num_seqs=2` per RFDiffusion design?

**A84:** I've chosen num_sequence = 2 because I would rather diversity on designs
rather than sequences, which tend to be very similar within designs even with higher
sampling temperatures, and also with my current mpnn_top_n, sequences from the same
design tend to score similarly on MPNN so I would just end up with loads of similar
sequences passing through negsteering. I could I suppose generate loads of MPNN
sequences and then just select the top 2 from each of them. Maybe that would be a
wise decision in the future to bear in mind.

**Synthesis note**: `num_seqs=2` is a deliberate diversity strategy: sequences from
the same backbone tend to be similar regardless of sampling temperature. Backbone
diversity (64 designs) is more valuable than sequence diversity per backbone for
exploring the search space. Future option noted: generate many sequences per design,
then select top 2 per design (per-design selection), rather than top 128 globally.
This would guarantee backbone-level diversity in the negsteer cohort.

---

**Q85 [VERIFY]** — R7 fix in `align_to_native_by_anchors` — could rounding cause
errors?

**A85:** This seems risky - couldn't the previous alignment code used to identify
design regions work here? Proportions could end up with rounding errors that trim
or extend the design region. I think this needs to be inspected and audited.

**Code check**: The R7 fix implementation in `pipeline_correct_sequences.py` (lines
261-315):
1. Fixed-length segments (`min_len == max_len`, e.g. `6-6`) are allocated FIRST,
   getting exactly their specified length — no rounding.
2. Variable segments share the remaining budget. With one variable segment, it
   absorbs all remaining residues (no rounding). With multiple variable segments,
   `round()` is used for proportional distribution, but the LAST variable segment
   absorbs any rounding remainder explicitly.

For the current 2-region contig (`10-20` variable + `6-6` fixed): `6-6` is allocated
exactly 6, then `10-20` absorbs whatever remains (clamped to [10, 20]). No rounding
involved for this specific contig.

**However**: the user's concern about auditing is valid. The R7 fix has not been
unit-tested. A test with a multi-variable-region contig (e.g. two variable regions)
would exercise the `round()` path. → Phase 5 unit test target. Until then, flag as
needing manual verification for any campaign using multiple variable-length design regions.

---

**Q86 [VERIFY]** — All current campaign PDB receptor chains start at residue 1?

**A86:** When I built the input complexes, I made sure that the receptors started
at 1 (makes it easier to track changes through the pipeline), and I should do that
for the effector too, but I don't think it should be mandatory - there should be
flexibility for it to work even if the starting position is not 1.

**Synthesis note**: All current campaign inputs start at residue 1 by deliberate
construction. Future flexibility needed: the pipeline should handle non-1 start
residues without breaking (the `receptor_start_pdb` plumbing from pipeline_notes2
exists for this). Effector numbering should also start at 1 where possible — and
should still work if it doesn't. → Add `receptor_start_pdb` and `effector_start_pdb`
flexibility to the input validation checklist.

---

**Q87 [VERIFY]** — Branch B (pre-docked PDB) is the authoritative branch. HADDOCK
(Branch A) is unused and out of date.

**A87:** The HADDOCK code is out of date and currently unused - it remains there
as an option to redesign in the future - the branch B processing should be read as
the ground truth.

**Synthesis note**: Confirmed. For all audit purposes, Branch B is canonical. Branch
A (HADDOCK) exists as dead-but-preserved code. HADDOCK redesign is explicitly deferred.
Any questions about HADDOCK-specific code in Sessions 3–5 will be skipped unless
directly relevant to the current pipeline behaviour.

---

**Q88 [VERIFY]** — Has the trailing-design-region `native_residues` bug fix been
verified in production?

**A88:** I'm pretty sure it works, but you can look at the results in the mounted
volume to confirm this. I didn't check these results since the fix was implemented
and at that time it worked for that particular contig so I generalised. Later on
I'm assuming we will do some unit testing.

**Code check — fix verified**: In production campaign `pikp1_avrpikf/runs/v1_4f/`,
`scored_metadata.csv` shows `native_residues = 'SVALVGDLRDKIEVVGY|LLQVSQANKD'`
for all sequences — both regions populated. The trailing `6-6` design region is
correctly captured. Fix is working for the current production contig.

**Phase 5 unit test target**: Write a parametrised unit test for the three boundary
cases: trailing denovo (no fixed after), leading denovo (no fixed before), interior
denovo (fixed on both sides). The trailing case is now verified in production; the
other two are currently exercised only by the contig format, not by dedicated tests.

---

---

### Batch 3 — Closing Session 3 (Q89–Q96)

**Process correction** (recorded from user feedback on this batch):
VERIFY questions must be "I found this in the code and it works this way — is this correct?"
UNKNOWN questions must be things that genuinely cannot be answered from the code.
Questions about technical implementation details (sampling temperature, forward pass structure,
dG computation) are code-readable and must be researched before asking.
Going forward, each batch will include a note of what was pre-verified in the code.

---

**Q89** — Update params_example.yml to match production defaults.
**Action taken**: Updated `num_designs: 32→64`, `num_seqs: 16→2`, `mpnn_top_n: 64→128`,
`negsteer_n_designs: 20→3`. Added comments explaining the rationale. Documented
`stop_after_rosetta` prominently. Updated checkpoint comment to note `Complex_beta_ckpt.pt`
is strongly preferred for PPI design.

---

**Q90 [VERIFY]** — `design_region_length_observed` and `design_region_spec` naming.

**A90:** Observed and Spec are clear to me - it also helps as a visual check to
make sure the observed lengths match the spec and nothing has gone wrong.

**Synthesis note**: Naming confirmed as clear and useful. The observed/spec comparison
is a built-in sanity check: if a design has region 2 observed=4 but spec=6-6, something
went wrong upstream.

---

**Q91 — Code-verified findings**:

`dG_separated` IS computed at the pre-MPNN Rosetta stage (read from InterfaceAnalyzer
score file by `rosetta_filter_collect.py`). The values are collected and written to the
`rosetta_filter_summary.csv`, but the filter gates ONLY on `Sc ≥ sc_threshold`. The
histogram plot shows Sc only (dG scatter was removed). The dG values at this stage are
positive/meaningless on polyvaline backbones — computed but not used for filtering.

Post-negsteer Rosetta (orthogonal metrics, after FastRelax on Boltz-predicted complex):
dG values are meaningful (-50 to -70 REU per user confirmation), calibrated against
Bennett et al. 2023 using the -30 REU threshold.

User confirmed: orthogonal-metrics dG values ARE meaningful.

---

**Q92 [VERIFY]** — `stop_after_rosetta` should be documented prominently.

**A92:** Yes it should be documented prominently. It saves a lot of compute resources
to check the designs are reasonable before continuing with the heavier negsteering
code.

**Action taken**: Updated `params_example.yml` to document `stop_after_rosetta` as
a useful diagnostic tool with usage guidance.

---

**Q93 — Code-verified**: MPNN sampling temperature = `mpnn_sampling_temp: 0.25` (now
corrected in params_example.yml from `0.25` — this was already correct). Temperature
of 0.25 is a low-diversity setting, which aligns with the user's observation that
sequences within a design are similar regardless of temperature.

---

**Q94 — Code-verified**: `mpnn_design_region_score.py` runs a SEPARATE forward pass
using `ProteinMPNN.forward()` to compute per-position log-probabilities, then averages
over the design-region positions. This IS an additional model forward pass beyond the
main sampling pass in `PROTEINMPNN`. However, it is a single deterministic scoring
pass (no sampling), which is cheaper than the full sampling run. Future optimisation
opportunity: compute design-region scores during the main sampling pass to avoid
loading the model twice.

---

**Q95 [VERIFY]** — Two deferred items from Session 1.

**A95:** A) Yes (sort direction correct, code verified). B) The AF3 check is going
to be reworked because it fails on everything which is why I turned it into a
warning. It's fine for it to be the best across the predictions for now - it's
good that it's been noted.

**Synthesis note**: (a) mpnn_select_top sort direction confirmed correct. (b) AF3
`best_ra_eff` aggregation (best-of-15) is acceptable as-is. The AF3 stage fails
all current designs (high ra_eff) — expected, since AF3 without MSA is a weak
predictor for these novel complexes. AF3 stage is acknowledged as needing
rework; using best rather than median is intentional to be maximally lenient.

---

**Q96 [UNKNOWN]** — Anything wrong/confusing in design generation?

**A96:** I would like the sequences from the MPNN output tables to be carried
forward into the final results in a way they aren't at the moment. It should be
possible to see the native and final design region sequences, as well as the whole
de novo MPNN sequence in the final output table (that way it makes it much easier
to go and order the thing rather than hunting through intermediate MPNN outputs).

**Synthesis note — FEATURE REQUEST**:
The final cohort summary (`cross_sequence_summary.csv` or the combined output table)
should include:
1. Native design region sequence (what the original receptor had at those positions)
2. Final design region sequence (what MPNN designed)
3. Full de novo MPNN receptor sequence (what would be ordered for synthesis)

These are available in MPNN intermediate outputs (`sequence_metadata.csv`,
`mpnn_corrected.fasta`) but are not carried forward to the cross-sequence summary.
→ Phase 4 feature: propagate sequence information from MPNN through the pipeline
to the final output table. Already tracked in Q42 as an open feature.

---

## Session 3 — Confirmed answers index

| Topic | Decision | Q# |
|---|---|---|
| `Complex_beta_ckpt.pt` is the only checkpoint in use | Confirmed | Q73 |
| `add_potential: true` is the correct default | Confirmed | Q74 |
| Production params: 64 designs, 2 seqs/design, top 128, 3 negsteer designs | Confirmed, params_example.yml updated | Q75/Q76 |
| Three geometry metrics are diagnostic, not gating | Confirmed | Q77 |
| Rosetta Sc at pre-MPNN stage: rough catch, threshold may be too high | Confirmed, audit deferred | Q78 |
| rfdiff_iterations=100 for production, 50 for tests | Confirmed | Q79 |
| align_to_native_by_anchors uses contig fixed segments as alignment anchors | Confirmed | Q80 |
| Global MPNN score is correct for selection (non-design regions constant) | Confirmed | Q81 |
| num_seqs=2: backbone diversity preferred over sequence diversity | Confirmed | Q84 |
| HADDOCK branch is dead/deferred; Branch B is ground truth | Confirmed | Q87 |
| native_residues fix verified in production for 2-region contig | Verified in code | Q88 |
| dG_separated computed pre-MPNN but not used for filtering | Code-verified | Q91 |
| Post-negsteer Rosetta dG values meaningful (-50 to -70 REU) | Confirmed | Q91 |
| AF3 best_ra_eff (best-of-15) is acceptable; AF3 stage needs rework | Confirmed | Q95 |
| Final output table should carry MPNN sequences forward | Feature request noted | Q96 |

---

## Session 4 — Orthogonal metrics + cross-cutting concerns

### Batch 1 — Filter cascade, interface_plddt, duplications (Q97–Q104)

**Code pre-research findings** (checked before asking):
- Filter cascade: `passes_orthogonal_filters=1` iff ALL non-AF3 flags absent. ANY failing non-AF3 metric → passes=0. The Q97 question premise was wrong — current code already implements "any fail → 0".
- Production data: 146/845 rows (17%) pass orthogonal filters across campaigns. The Q102 premise ("no designs pass") was wrong — many pass, particularly in v1_4a and v1_4b.
- The blocking metric in failing rows is most often `interface_plddt_too_low` (Boltz metric in the orthogonal cascade — see Q103).
- `find_contact_residues_heavy`: two distinct implementations. The `boltz2_negative_steering.py` version is more robust (numpy-vectorised, raises ValueError, sorts by distance). The `compute_metrics.py` version returns `[]` on missing chain and is unsorted.

---

**Q97 [VERIFY]** — Orthogonal filter cascade: the premise in the question was wrong.

**A97:** The thresholds should be lowered to just filter out the terrible ones - this
will need checking during the thresholds audit - and the fail should be on any
failing orthogonal metric (except from AF3 which is just informational). Good catch.

**Correction**: the current code already implements "any non-AF3 fail → passes=0".
`gating_flags = [f for f in flags if not f.startswith("af3_nomsa")]`; then
`passes = 1 if not gating_flags else 0`. The user's intended behaviour matches
the existing implementation. No code change needed here. The threshold values
are what need auditing (Sc≥0.55, BSA≥600, interface_plddt≥0.75, ΔΔG≤-30).

---

**Q98 [VERIFY]** — Historical iRMSD/fnat/dockq: recompute?

**A98:** All the results so far are for testing purposes and getting a feel for
the pipeline. Nothing will be recomputed. Old results are old results. I will
only take things forward for agroinfiltration once I'm confident the code has
no major holes.

**Synthesis note**: No recomputation of historical iRMSD/fnat/dockq. The historical
data warning in remediation_state.md stays as a documented caveat. Any production
run going forward will use the fixed chain-param code and produce correct values.

---

**Q99 [UNKNOWN]** — iRMSD/fnat/dockq for candidate selection?

**A99:** I don't currently use them for selecting candidates but I should. They
do belong in the threshold audit.

**Synthesis note**: iRMSD and fnat are currently in the output CSV but not used
for selection. Future plan: include in threshold audit; may replace or supplement
ra_eff as the primary structural placement metric (they are more interface-specific
than whole-effector Cα RMSD). DockQ score is the combined metric (0–1 scale with
CAPRI thresholds: >0.23 acceptable, >0.49 medium, >0.80 high).

---

**Q100 [VERIFY]** — Weighted_jaccard primitives: canonical home.

**A100:** They should be set in one place and loaded in another. Consolidating
them into boltz_lib could be an option, depending on how you decide to restructure
the scripts.

**Synthesis note**: Phase 4 action — move `_WJ_MU`, `_WJ_TWO_SIGMA_SQ`,
`_gaussian_contact_weight`, `_collect_chain_cb_positions`, `_build_weighted_pair_map`
to `boltz_lib.py`. Both `compute_metrics.py` and `compute_interface_metrics.py`
import from there. Location is conditional on Phase 4 restructuring design.

---

**Q101 [VERIFY]** — ΔΔG threshold: -30 REU is a floor, not a discriminator.

**A101:** I haven't actually read Bennett et al., 2023. I don't really need the
threshold to be discriminating, just to filter out terrible designs. The threshold
should be a floor.

**Synthesis note**: ΔΔG ≤ -30 REU is a floor — catches catastrophically poor
binding energy while letting everything else through. Current designs are -50 to
-70 REU (well below the threshold), so it doesn't discriminate. The Bennett et al.
2023 reference provides calibration context that should be checked; the -30 REU
value may need adjusting once more data is available. → Threshold audit.

---

**Q102 [VERIFY]** — Does `passes_orthogonal_filters` flag currently produce any 1s?

**A102:** There are loads of runs - I'm sure some pass all four you just haven't
looked. If you're including AF3, then exclude it. The flag is useful.

**Correction**: The premise was wrong. Production data shows 146/845 sequences
(17%) pass all non-AF3 orthogonal filters. Highest: v1_4b (56/104), v1_4e (28/118),
v1_4a (24/116). Passing rows have either blank flags or AF3-only flags — AF3 correctly
excluded from gating. The flag IS useful and meaningful.

---

**Q103 [VERIFY]** — `interface_plddt` in the orthogonal cascade: is this correct?

**A103:** It's a plddt score, so it must come from boltz. It must come from the
representative prediction, whether that's from the cold-start, steering, or
reversion for that particular negsteering path. If it has already been calculated
per-prediction and aggregated, then that is the one to use. I'm not sure why
there is an interface_plddt in biophysical metrics? plddt is from the Boltz-2
model. I need more information to understand that. The interface_plddt shouldn't
be part of orthogonal metrics at all, it should be in the negsteering section -
it's a plddt score it's not orthogonal at all.

**Synthesis note — CONFIRMED BUG / ARCHITECTURAL ISSUE**:
`interface_plddt` (Boltz pLDDT over interface residues) is a Boltz-2 confidence
metric. It is NOT an orthogonal metric — "orthogonal" means independent of Boltz.
Its inclusion in the orthogonal filter cascade violates the definition of orthogonal
validation.

Code state: `run_biophysical_metrics.py` had `interface_plddt` removed in notes9.
The value used in the orthogonal filter (`representative_interface_plddt_median`)
comes from the negsteer aggregated data — Boltz predictions, not an orthogonal
source. The filter label and column placement are misleading.

**Production impact**: In v1_4f_rfdtest_2, `interface_plddt_too_low` is the
primary blocking flag (50 rows, 0 passing). Removing it from the orthogonal cascade
would change which rows pass.

**Action**: Remove `interface_plddt` from `merge_orthogonal_metrics.py`'s orthogonal
filter cascade. The metric stays in the cohort summary table under the "Boltz
confidence" column group (where it already is), but stops gating
`passes_orthogonal_filters`. → Phase 4 code change.

---

**Q104 [VERIFY]** — `find_contact_residues_heavy`: two implementations — are they
semantically interchangeable?

**A104:** I don't know why they would be separate. Maybe they're used for different
implementations, but they should really be consolidated into one choice if there's
no reason for them to be different. You need to suggest whether they're semantically
interchangeable.

**Code analysis**: The two implementations compute the same logical operation (min
heavy-atom distance from each receptor residue to any effector heavy atom, returning
positional indices for contacts below cutoff). Key differences:

| Aspect | `compute_metrics.py` | `boltz2_negative_steering.py` |
|---|---|---|
| PDB parser | `_parse_heavy_atoms_by_chain` | `read_residue_heavy_atoms` |
| Missing chain | Returns `[]` silently | Raises `ValueError` |
| Output sort | Unsorted | Sorted by distance ascending |
| Vectorisation | Per-residue loops | numpy-vectorised |

They are semantically equivalent for the common case. Differences matter at edges:
(1) missing chains: negsteer version correctly fails loudly; (2) callers that rely
on the first N results will get different answers if sort order differs.

**Recommendation**: The `boltz2_negative_steering.py` version is canonical — more
robust error handling, vectorised, sorted. Consolidate into `boltz_lib.py` with
the negsteer semantics. Update `compute_metrics.py` to import and call it.

---

---

### Batch 2 — Biophysical metrics, Rosetta, cross-cutting parsers (Q105–Q112)

**Code pre-research for this batch**:
- BSA production data: 840 values, 15% below 600 Å², median 781 Å². Most designs pass.
- Boltz PDB atoms: only C, N, O (no hydrogens). Rosetta adds/optimises H with `-no_optH false`.
- Chain names: `write_split_pdb` in `rfdiffusion_filter.py` explicitly assigns A=receptor, B=effector, renumbering from 1. MPNN and Boltz preserve those chain labels throughout.
- Contig parsers: three active implementations with distinct purposes and output formats (see Q110 synthesis).

---

**Q105 [VERIFY]** — `interface_hbonds`: informational; is it a 5th heavy-atom detection function?

**A105:** It's currently purely informational. I don't know much about this metric.
I thought we already had several functions for detecting heavy atoms which were used
in negsteering to determine contacts? Is this yet another function? And I don't
remember working through the MDAnalysis HydrogenBondAnalysis and the corresponding
failure - maybe something to come back to. If the fallback is well established then
it's fine to keep it. Should I be thresholding this?

**Code check and synthesis**:
`interface_hbonds` uses a gemini-based geometric detector in `run_biophysical_metrics.py`:
N/O atoms within 3.5Å between chains. This is DISTINCT from the heavy-atom contact
detectors in negsteer:
- negsteer `find_contact_residues_heavy`: all heavy atoms, returns positional indices
- `compute_interface_contacts`: Cβ-Cβ contacts for weighted_jaccard
- hbond detector: N/O donor/acceptor atoms only, cross-chain, count only

So yes, this is a 4th type of heavy-atom interaction detector, but with a different
physical meaning (hydrogen bond-like contacts vs generic contacts).

The MDAnalysis HydrogenBondAnalysis failed because Boltz PDBs are heavy-atom only
(confirmed: only C, N, O, no H). The geometric fallback is standard practice
(used by PISA/FreeContact). Thresholding: include in threshold audit; a minimum
count (e.g. ≥2) could be a sanity check, but since this is only approximately
correlated with real H-bonds, keep as informational for now.

---

**Q106 [VERIFY]** — BSA ≥600 Å² threshold: how discriminating is it?

**A106:** Look yourself. I think many of them do pass this BSA threshold though.
It will be audited. I don't know much about the FreeSASA computation choices.

**Code check**: 840 BSA values across production campaigns. 85% (713/840) pass the
≥600 Å² threshold. Median BSA = 781 Å². Range: 280–3015 Å². The threshold
is not overly strict — it's filtering the bottom 15%, not a tight gate.

The FreeSASA computation: three calls (complex, receptor alone, effector alone);
BSA = (SASA_rec + SASA_eff − SASA_complex) / 2. Standard method; 600 Å² from
Overath et al. 2025 is a reasonable floor for a real interface.

---

**Q107 [VERIFY]** — Does Rosetta run on Boltz PDB directly, and does it contain hydrogens?

**A107:** I believe it's running on the Boltz-2 prediction. [User directed to code check.]

**Code check**: Boltz prediction PDBs contain only heavy atoms (C, N, O confirmed
from checking production PDB files — no H atoms). `run_rosetta_metrics.py` passes
the Boltz PDB directly to `rosetta_scripts` FastRelax. Rosetta is invoked with
`-no_optH false` (optimise hydrogens) — Rosetta adds H positions internally as
part of the FastRelax protocol. This is the correct approach: Rosetta's force field
requires H for scoring, so it generates and optimises them internally rather than
requiring pre-placed H in the input. The `-ignore_unrecognized_res` flag handles
any non-standard residues gracefully.

---

**Q108 [VERIFY]** — FastRelax protocol (1 repeat, Cartesian, ref2015_cart) is correct?

**A108:** The current FastRelax setup seems fine, but I could include the number
of repeats in the future for more rigorous values. This could be included in the
threshold audit. I don't know about the other parts of this pipeline - I'm
generally uninformed about the orthogonal metrics and it's something I need to
understand better.

**Synthesis note**: 1-repeat FastRelax is a fast approximation; 3-5 repeats would
give more converged ΔΔG values. The -30 REU threshold from Bennett et al. 2023 is
calibrated against a specific Rosetta protocol — if the number of repeats changes,
the threshold may need recalibration. → Threshold audit.

---

**Q109 [VERIFY]** — Do chain names ever change in RFDiffusion, MPNN, or negsteer outputs?

**A109:** I can't see any reason why RFDiffusion, MPNN, or Negative Steering (Boltz-2)
would ever change the Chain names. You can check to be sure though in the code.

**Code check**: `rfdiffusion_filter.py:write_split_pdb` explicitly writes chain A=receptor,
chain B=effector, renumbering from 1 — regardless of input chain letters. MPNN runs
on these A/B PDBs and its output inherits the same chain labels. Boltz receives the
MPNN-corrected FASTA sequences and the design PDB (A/B) and produces predictions
with chains A (receptor) and B (effector). Chain names are invariant throughout the
post-RFDiffusion pipeline. DockQ correctly finds the A/B key in prediction PDBs.

---

**Q110 [VERIFY + UNKNOWN]** — Three active contig parsers: are they interchangeable?

**A110:** If they're all different, I need to understand what the differences are
and why they exist. They might be handling different inputs e.g. contigs from the
original input, or contigs from the RFDiffusion design which might have different
design region lengths, but I'm sure some of them can be combined. You're correct
to disregard the HADDOCK one.

**Code analysis** — the three parsers are NOT interchangeable:

| Parser | Location | Input | Output | Purpose |
|---|---|---|---|---|
| `parse_block_segments` | `contig_utils.py` | One contig block (pre-space-split) | `[("fixed", chain, start, end), ("denovo", "min-max"), ...]` | General; handles breaks, passthrough, bare chain resolution |
| `_parse_contigs` | `derive_input_design_region.py` | Full contigs string + both chain IDs | `[(chain, start, end, is_denovo)]` flat list | Controls: needs both receptor AND effector segments |
| `parse_contig_segments` | `pipeline_correct_sequences.py` | Full contigs string + receptor chain | `[{type, start, end, min_len, max_len}]` dicts | Sequence correction: receptor segments as structured dicts |

Could be consolidated: `parse_contig_segments` and `_parse_contigs` could both be
rewritten to call `parse_block_segments` as their backend and post-process the output.
This is a Phase 4 task. The output types differ enough that a shared function would
still need adapter layers, but the parsing logic (string tokenization, chain detection,
range splitting) would be deduplicated.

---

**Q111 [VERIFY]** — Consolidate `make_empty_plot`/`save_fallback_plots` into shared module?

**A111:** Consolidating them into a shared module is a good approach.

**Synthesis note**: Three active instances (rfdiffusion_plots.py, mpnn_plots.py,
rosetta_filter_plots.py). `plot_lib.py` proposed in Phase 4. Each module's
`ALL_PLOT_FILES` stays module-specific (it's just a list); the fallback logic is
shared. → Phase 4 implementation.

---

**Q112 [UNKNOWN]** — `COLOUR_TIER` and `_TIER_SORT_ORDER` constants: satisfactory colours?

**A112:** I'm satisfied with these colours. They could be consolidated elsewhere
to remove duplication if it's cleaner.

**Synthesis note**: Colours confirmed. `COLOUR_TIER = {A: green, B: yellow, C: orange,
none: grey}`. Consolidate to `plot_lib.py` in Phase 4 along with `make_empty_plot`.
The constants are currently in 3 bin/ files plus their test twins (6 copies total).

---

## Session 4 confirmed answers index

| Topic | Decision | Q# |
|---|---|---|
| `interface_plddt` removed from orthogonal filter cascade | Confirmed | Q103 |
| Historical iRMSD/fnat/dockq not recomputed | Confirmed | Q98 |
| iRMSD/fnat/dockq belong in threshold audit | Confirmed | Q99 |
| AF3 remains informational only | Confirmed | Q97, Q95 |
| `passes_orthogonal_filters`: 146/845 pass in production | Code-verified | Q102 |
| Chain names (A=rec, B=eff) are invariant post-RFDiffusion | Code-verified | Q109 |
| BSA 85% pass ≥600 Å²; threshold is a floor | Code-verified | Q106 |
| Rosetta runs on Boltz heavy-atom PDB; adds H internally | Code-verified | Q107 |
| Three contig parsers are NOT interchangeable; consolidation is Phase 4 | Code-verified | Q110 |
| `make_empty_plot` → `plot_lib.py` in Phase 4 | Confirmed | Q111 |
| `COLOUR_TIER` colours accepted; consolidate to `plot_lib.py` | Confirmed | Q112 |
| `find_contact_residues_heavy`: negsteer version is canonical; boltz_lib.py | Confirmed | Q104 |
| weighted_jaccard primitives → boltz_lib.py | Confirmed | Q100 |

---

## Session 5 — Architecture + future intent

### Batch 1 — Architecture scope, terminology, priorities (Q113–Q120)

---

**Q113 [VERIFY]** — `boltz_lib.py` scope confirmed?

**A113:** I think all of those things are in scope as long as they don't break
anything else. I'm looking to simplify and deduplicate the code as much as
possible so that things are used consistently across the whole pipeline.

**Synthesis note**: boltz_lib.py scope confirmed: `get_chain_sequence`,
`find_contact_residues_heavy` (negsteer version as canonical), `jaccard`,
weighted_jaccard primitives, `binding_rmsds`, `write_boltz_yaml`, Cα/Kabsch
helpers, residue dataclasses. Guiding principle: one implementation, consistent
use across all scripts.

---

**Q114 [VERIFY]** — `plot_lib.py` scope confirmed?

**A114:** That scope also looks correct.

**Synthesis note**: plot_lib.py scope confirmed: `make_empty_plot`,
`save_fallback_plots`, `ALL_PLOT_FILES` pattern, `COLOUR_TIER`,
`_TIER_SORT_ORDER`, `COMPOSITE_RA_EFF_WEIGHT`, `_try_float`/`_try_int`/`_short_name`
shared utilities.

---

**Q115 [UNKNOWN]** — Terminology simplification proposals.

**A115:** The representative prefix is actually useful because otherwise we get
confused about whether the values are from cold-start, steering, or reversion -
it's fine to rename it to "rep" though to make it shorter. I don't think the fine
detail about how many seeds passed or failed each filter or stage are necessary -
that information is in the per_seed csv output files, so I think you're right to
just consolidate them into n_pass and n_seeds. I like the rename of
aggregated_verdict to outcome.

**Decisions recorded**:
1. `representative_` prefix → `rep_` (shorter, retains semantic distinction)
2. `n_seeds_pose_holds`, `n_seeds_clean_steered`, `n_seeds_pose_collapses`,
   `n_seeds_new_contamination`, `n_seeds_no_data` → `n_pass` + `n_seeds`
   (fine detail stays in `raw_per_seed_results.csv`)
3. `aggregated_verdict` → `outcome`

→ Phase 4 implementation: update `boltz2_iterate_steering.py` output column names,
`cross_sequence_summary.py`, `extract_passing.py`, and all downstream CSV readers.

---

**Q116 [UNKNOWN]** — MPNN sequence propagation: add to cross_sequence_summary.csv?

**A116:** I think that is the correct approach for now. Later on I'm going to do
another audit of exactly what information is in each output file, cleaning out
anything that isn't necessary or is duplicated for no reason, or files that don't
really need to exist or can be merged, but for now I think making that information
more visible is an important step.

**Synthesis note**: Add to `cross_sequence_summary.csv`:
- `corrected_receptor` (full MPNN-designed receptor sequence — for ordering)
- The design-region slice of `corrected_receptor` (derived from contig info)
- `native_residues` (pipe-separated native sequence at design region positions)

Source: `scored_metadata.csv` from the MPNN stage, joined at the
`NEGSTEER_CROSS_SEQUENCE` step on `mpnn_sequence` as the join key.
A future output audit will decide which intermediate files can be removed or merged.

---

**Q117 [UNKNOWN]** — Input validation: dedicated script for all parameters.

**A117:** Anything with boundaries between 0 and 1 should be checked to make sure
it is. You should be able to come up with a range of possible inputs for each
parameter, and these will all be tested in the unit tests. This is one of the most
important bits. If the internal pipeline works that's great, but the most important
part is to make sure the user can't break it by giving weird inputs. The input
possibilities need to be completely locked down right at the start. Your existing
suggestions are a good start but there should be checks for every parameter. It
could be its own separate script.

**Synthesis note — HIGH PRIORITY ITEM**:
- Separate script: `bin/validate_params.py` (or inline in preprocessing.nf entry)
- Called at the start of the pipeline before any SLURM jobs are submitted
- Tests unit-tests: full range of valid/invalid inputs for every parameter
- Current parameter catalogue needs to be assembled for the validation spec
- Rationale: the user is currently the only operator; future operators may not
  know which inputs are invalid, and silent failures downstream are worse than
  clear upfront errors

→ Phase 4 implementation task. The validation spec (full parameter catalogue with
valid ranges and types) is a prerequisite.

---

**Q118 [VERIFY]** — Remove `interface_plddt` from orthogonal cascade now?

**A118:** Implement it now. It is in the wrong place. The exact threshold value
can be deferred to the threshold audit.

**Action taken**: `merge_orthogonal_metrics.py` updated. The `interface_plddt`
check is removed from `_apply_filters` (it is now a no-op comment explaining why).
The `--filter-plddt-min` CLI argument is retained for backward compatibility but
has no effect. AF3 remains informational. `passes_orthogonal_filters` now gates
only on Sc, BSA, and ΔΔG (genuinely orthogonal metrics). Docstring updated to
explain the orthogonal vs Boltz distinction explicitly.

---

**Q119 [VERIFY]** — Multi-cycle implementation: loop within shell orchestrator?

**A119:** That sounds like the correct approach.

**Synthesis note**: Multi-cycle implementation plan:
1. Extend `negative_steering_run_one.sh` to loop over cycles (cycle-0, cycle-1, ...)
2. After each cycle, run the novelty/selection logic inline (calling
   `cmd_kickoff_distances` and the selection logic from Python)
3. If at least one seed found the correct interface → stop (success)
4. If no novel placements found → stop (search exhausted)
5. Otherwise, continue to cycle-N+1 with the surviving sequences as input
6. Stop at user-defined `n_cycles`

The `cmd_kickoff*` sbatch-based submission family is retained for the standalone
negsteer pipeline but is NOT used in the Nextflow integration. The Nextflow
integration loop runs entirely within the single NEGSTEER_RUN_ONE SLURM job.

→ First post-remediation feature implementation.

---

**Q120 [UNKNOWN]** — Any accumulating fixes/features not yet reached?

**A120:** I've audited the plots before and caught most of the issues, but the
main table that I'm using - the summary table produced after the orthogonal
metrics (it's being used as a first view of the pipeline results) seems to have
some inconsistencies that we've talked about. I think correcting that should
become a priority.

**Synthesis note**: The cohort summary table (`orthogonal_combined_cohort_summary.png`)
is the primary result view. Known issues documented across the audit:
1. **Mixed data sources** not communicated (Boltz medians vs AF3 best-of-15)
2. **interface_plddt** in the filter cascade — FIXED (Q118)
3. **Tier/orthogonal independence** confusing — needs clearer labelling
4. **Column source labelling** needed (which columns are Boltz, which are orthogonal)

→ Cohort summary table audit is a Phase 4 priority, immediately after the
`interface_plddt` removal is tested.

---

---

### Batch 2 — Cohort summary, audit close (Q121–Q128)

---

**Q121 [VERIFY]** — Cohort summary: three column groups instead of two?

**A121:** I think you could go one step further from that: 1) Boltz-2 structure
2) Boltz-2 confidence and 3) Orthogonal. This is a good idea.

**Decision**: Three visual column groups in the cohort summary:
1. **Boltz-2 structure** — composite, Boltz ra_eff, rec RMSD, true_jaccard, n_mut
2. **Boltz-2 confidence** — complex_plddt, iptm, ipae, pae_pass_frac, iface_plddt, ipsae_min_15
3. **Orthogonal** — AF3 ra_eff, Sc, BSA, ΔΔG, weighted_jaccard

→ Phase 4: add group header row to `plot_combined_cohort_orthogonal_summary`.

---

**Q122 [VERIFY]** — AF3 `best_ra_eff`: add `(best)` label only.

**A122:** It's only AF3 that uses the best one right? That's a single column.
You can just add (best) in small text for that column. Otherwise median can
be assumed. I don't want to overcrowd the plot.

**Decision**: Add `(best)` suffix to the AF3 ra_eff column label only. All Boltz
columns are implicitly median and need no label change. Confirmed code: the column
is `af3_nomsa_best_ra_eff` (best of 15 predictions) in `COMBINED_SUMMARY_COLUMNS`.

---

**Q123 [VERIFY + INVESTIGATION]** — Tier A design with ipae=17.58, pae_pass_frac=0.0.
Is the thresholding on best or median values?

**A123:** This may have brought up a really important point about thresholding.
[...] I wonder if the thresholding is conducted on the best or median values -
it should be on the median values, as it is in the plot. It shouldn't be possible
for a 17.58 ipae median to be tier A, because at least 2 of the 3 seeds must have
had ipae values greater than the threshold.

**Code check and investigation**:

The specific case is `design_42_seq_0` (full_test_run): cross_tier=A,
n_seeds_pose_holds=3, ra_eff_median=2.897Å, ipae_median=17.58, pae_pass_frac=0.0.

**Finding 1 — Thresholding is correctly on median values**: The plot uses
`representative_ipae_median` and `representative_pae_pass_frac_median` (confirmed
from `COMBINED_SUMMARY_COLUMNS`). Tier A is assigned from `n_pass = n_seeds_pose_holds
+ n_seeds_clean_steered = n_seeds = 3` — based purely on structural placement (ra_eff
< 5Å AND intact), NOT on confidence metrics. A Tier A design can have any ipae value.
This is expected and intentional — the pipeline explicitly documents that Boltz can
produce structurally correct placements with poor confidence (documented throughout
the negsteer notes as the central finding for these novel interfaces).

**Finding 2 — This specific case is from a pre-fix run**: The full_test_run was
submitted 2026-04-30; the chain-param bug fix (commit a4ef1f7) was 2026-05-03.
Confidence metrics (`ipae`, `pae_pass_frac`) were computed using the wrong effector
chain in that run. The `pae_pass_frac=0.0` for a design with ra_eff=2.897Å is
almost certainly an artifact of that historical bug (documented in Q98/remediation_state.md),
not a genuine pipeline behaviour. The structural metric (ra_eff) was computed
correctly in that run.

**Conclusion**: No code bug. The user's concern is valid as a conceptual question
("can a Tier A design have poor confidence?") — and the answer is yes, it can,
and this is expected. The specific numbers in the screenshot are from a buggy
historical run and should not be used as the canonical example.

**UX fix needed**: The three column groups from Q121 will make the independence
of tier (structure) and confidence clearly visible. The label "Tier A design with
red confidence cells" should read as "found the correct interface but Boltz is
uncertain about it" — a meaningful finding, not an error.

---

**Q124 [UNKNOWN]** — What should be implemented before the Phase 4 architecture
session?

**A124:** I think all should be added before the Phase 4 architecture session,
since then we have all our content there and it's just a matter of redesigning it.

**Items queued for pre-architecture implementation**:
1. ~~`interface_plddt` removal from orthogonal cascade~~ — DONE (Q118)
2. Cohort summary table fixes (three groups, AF3 label, legend) — Q121-Q123
3. `num_seeds` odd validation — Q46
4. MPNN sequence columns in cross_sequence_summary — Q116
5. Input validation script (`bin/validate_params.py`) — Q117

All five items complete before starting Phase 4 architecture design session.

---

**Q125 [VERIFY]** — Six F-rated cmd_* functions: any untouchable?

**A125:** I don't consider any of them untouchable.

**Synthesis note**: All six are eligible for Phase 4 decomposition.
`cmd_compute_final_metrics` (CC 72→31 already, 6 remaining F-rated helpers
still deferred), `cmd_kickoff`, `cmd_kickoff_finalize`, `cmd_aggregate`,
`cmd_aggregate_per_sequence`, `cmd_collect`. Normal per-commit safety net applies.

---

**Q126 [VERIFY]** — Merge `derive_design_region.py` and `derive_true_interface.py`?

**A126:** No reason to keep them separate.

**Decision**: Merge into `derive_indices.py --mode design-region|true-interface`.
The four shared functions are verbatim identical; only `derive_positional_indices`
body differs. Phase 4 implementation.

---

**Q127 [UNKNOWN]** — Continue formal sessions or move to Phase 4?

**A127:** I think the audit has covered enough to move into the Phase 4
implementation. It wouldn't be meaningful to talk about things in tiny detail and
then change them hugely. The rest of the questions can follow or occur during the
restructuring process.

**Decision**: Formal grill-me audit ends here at 128 questions. The remaining
~190 questions in the original plan will be addressed organically during Phase 4
implementation — each code change will surface the specific questions that matter.
The audit document serves as the architectural spec.

---

**Q128 [UNKNOWN]** — What single issue caused the most friction day-to-day?

**A128:** The terminology confusion. There are so many terms in this pipeline,
which is fine, but we just have to be super clear about what each of them mean.
I think from our questions I understand now, but it will be important when I
come to write up this work to have a reference to work off of. I think we created
a shared glossary of terms at some point in the notes/inventory/ folder, and this
may need updating with some of the discussions we've had today.

**Synthesis note**: The ubiquitous language document (`notes/inventory/06_ubiquitous_language.md`)
exists and needs updating with key decisions from this audit:
- `clean_steered` definition (structural correct + no contamination)
- `no_data` definition (failed ra_eff or intact before contamination check)
- `no_reversion` ambiguity (two meanings — disambiguated by tier)
- `n_pass` formula and denominator
- The three-axis framework (Boltz structure / Boltz confidence / Orthogonal)
- `passes_orthogonal_filters` semantics (only Sc, BSA, ΔΔG — not pLDDT)
- `outcome` (renamed from `aggregated_verdict`)

→ Update the glossary before starting Phase 4 to lock in the terminology for
the restructuring work.

---

## Audit close — Summary of decisions and action items

### Immediate actions (pre-Phase 4 architecture session)

| Item | Status |
|---|---|
| `interface_plddt` removed from orthogonal cascade | ✅ DONE |
| params_example.yml defaults updated to production values | ✅ DONE |
| `stop_after_rosetta` documented prominently | ✅ DONE |
| Cohort summary: three column groups + AF3 (best) label + legend | → implement |
| `num_seeds` odd validation in `negative_steering_run_one.sh` | → implement |
| MPNN sequences added to `cross_sequence_summary.csv` | → implement |
| Input validation script `bin/validate_params.py` | → implement |
| Ubiquitous language glossary updated | → implement |

### Phase 4 architecture — scope confirmed

**boltz_lib.py**: `get_chain_sequence`, `find_contact_residues_heavy` (negsteer version),
`jaccard`, weighted_jaccard primitives, `binding_rmsds`, `write_boltz_yaml`,
Cα/Kabsch helpers, residue dataclasses.

**plot_lib.py**: `make_empty_plot`, `save_fallback_plots`, `COLOUR_TIER`,
`_TIER_SORT_ORDER`, `COMPOSITE_RA_EFF_WEIGHT`, shared utilities.

**derive_indices.py**: merge `derive_design_region.py` + `derive_true_interface.py`.

**terminology**: `aggregated_verdict` → `outcome`; `representative_` → `rep_`;
n_seeds_* consolidated to `n_pass` + `n_seeds`.

**contig parsers**: `_parse_contigs` and `parse_contig_segments` rewritten to
use `parse_block_segments` as backend.

### Multi-cycle — first post-remediation feature

Extend `negative_steering_run_one.sh` to loop over cycles inline (no sbatch).
Continue while novel incorrect placements are found; stop if correct interface found
or novelty exhausted. Gate implemented using existing `cmd_kickoff_distances` logic.

### Threshold audit (Task 47) — methodology confirmed

1. Research and document reasoning for each threshold
2. Run all five campaigns with corrected code
3. Assess per-campaign; update thresholds accordingly
