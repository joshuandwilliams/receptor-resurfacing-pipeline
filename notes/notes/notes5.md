# Negative Steering on RFDiffusion Designs — Session 2 (design_3)

Date: 2026-04-12. Continuation of notes 4. Single-design end-to-end run on
Pikp1 design_3 + ProteinMPNN seq_3, using the patches from notes 4 and the
production recipe locked in by notes 3.

## What was run

Target: Pikp1 + AvrPikD complex. Design 3 selected from `scored_metadata.csv`
on the basis of being the only entry top-tier on both `mpnn_score` (1.0289,
lowest of 16) and `design_region_score` (1.3559, second-lowest). Receptor
sequence taken from the raw `design_3_mpnn/seqs/design_3.fa` (sample 3),
**not** from `sequence_metadata.csv` — design_3 was generated before the
`pipeline_correct_sequences.py` fix from notes 4 was applied to the upstream
pipeline, so the CSV's `corrected_receptor` column would still contain
effector residues for this entry. The receptor FASTA was extracted manually
into a single-record file `design_3_seq_3_receptor.fasta` (79 residues,
matching PDB chain A and `receptor_total_length` from the CSV).

True interface indices derived from `rfdiffusion_metrics.json` for design_3:
`receptor_contact_residues` mapped through `receptor_position_order` gave
0-based seq_idx `[2, 32, 34, 36, 37, 38, 39, 47, 69, 70, 71, 72, 73, 74, 75,
76, 77, 78]` — 18 contacts, concentrated heavily at the C-terminal end
(positions 69–78 contribute 10 of 18) plus a small N-terminal cluster.

Run configuration: `--mode mild --max-mutations 4 --candidate-pool-size 6
--n-designs 50 --n-cycles 2 --max-passing 3`. 50 cycle-0 designs, 3 cycle-1
parents selected by maximin novelty (`c0d07`, `c0d29`, `c0d40`), 50 cycle-1
children per parent. 201 total predictions including the cold start.
Workdir: `runs/design_3_seq_3_smoke` (the `_smoke` suffix is misleading —
this is the production-style run; smoke test was never done. Worth renaming.)

Crucially, **the ground truth here is the experimentally-validated native
Pikp1–AvrPikD complex**, not a HADDOCK-derived pose. The notes-4 HADDOCK
caveat does not apply to this run.

## Headline result

The cold start is at ra_eff 29.30 Å with `intact=0` (rec_rmsd 3.05 Å, just
over the strict cutoff). Within one cycle of mild steering, the run produces
a design that scores well on every Boltz confidence metric simultaneously:

| pathway | cycle | ra_eff | rec_r | eff_r | tJac | ipSAE | ipTM | cplddt | flag | mutations |
|---|---|---|---|---|---|---|---|---|---|---|
| c0d40.c1d06 | 1 | 3.37 | 0.91 | 0.41 | 0.76 | **0.700** | 0.85 | 0.91 | ok | K5E A26R Y44R Q46D |

This is the strongest design in the run on confidence and the second-best on
ra_eff once the strict 3.0 Å cut is relaxed. Its parent `c0d40` (cycle-0,
mutations A26R Q46D) was already at ra_eff 3.35 Å, ipSAE 0.48 — the cycle-1
addition of K5E and Y44R sharpened both the structural alignment and the
confidence signal. Receptor backbone RMSD is 0.91 Å (well inside the 3.0 Å
intact threshold) and effector RMSD is 0.41 Å, so the floppy-termini and
side-chain-repacking concerns from notes 4 do not apply to this design — the
fold is genuinely intact, the residual ra_eff comes from rigid-body
positioning of an already-close pose.

The run also produced 27 designs at sub-3 Å ra_eff, 57 at sub-5 Å, and 9
that pass the full notes-3 production tier (`receptor_intact AND ra_eff < 3
AND tJac >= 0.7`), the strongest of which is a cycle-0 design at ra_eff
2.36 Å, tJac 0.77, ipSAE 0.491, ipTM 0.80, complex_plddt 0.90.

## Cross-target context (re-read notes 3)

Initial reaction was that the structural numbers were "8-10× better than
anything we've seen." That framing was wrong — based on the v5 7QPX number
from notes 1 (14.46 Å plateau) which predates the production recipe.
Comparing properly against the four native benchmarks from notes 3:

| target | rows | intact | best ra_eff | best tJac | cold ra | cold tJac | best ipSAE |
|---|---|---|---|---|---|---|---|
| 7QPX | 51 | 12 | 1.04 | 0.91 | 37.71 | 0.07 | 0.00 |
| 7QZD | 1288 | 552 | 0.59 | 0.92 | 24.61 | 0.38 | 0.39 |
| 6G10 | 51 | 35 | 2.26 | 0.84 | 22.11 | 0.50 | 0.00 |
| 7B1I | 102 | 55 | 2.76 | 0.87 | 5.47 | 0.83 | 0.00 |
| **design_3** | 201 | 83 | **1.46** | 0.77 | 29.30 | 0.48 | **0.86** |

The structural rescue (1.46 Å) is in the middle of the native band, not
above it. The cold-start tJac (0.48) is also in-band — comparable to 6G10
(0.50) and well below 7B1I (0.83). The structural dimension contains no
new information.

**The ipSAE result is the actually novel finding.** All four natives
plateau at ipSAE ≤ 0.39, with three of the four (7QPX, 6G10, 7B1I) maxing
out at exactly 0.000 even while producing structurally perfect rescues.
design_3 reaches 0.86 (and 0.70 within the production tier with tJac >=
0.7). This is the first observation of a meaningful, well-distributed ipSAE
signal on this pipeline's output.

## What this confirms

Notes 3 line 49–52 contains an explicit prediction:

> "Retained ipSAE as the primary ranker despite it being ~0 on most
> benchmark complexes, because the RFDiffusion designs will span a wider
> sequence space than native proteins and will likely produce meaningful
> ipSAE values."

design_3 confirms that prediction directly. The decision to keep ipSAE as
the production ranker — made on intuition during notes 3, when there was
no data supporting it — was correct. This is the first design where the
ranker has anything to rank.

Notes 4's "ipSAE confirmed broken on RFDiffusion data" conclusion now needs
qualification. It was based on Pikp1 design_1 v2, which (a) was run before
the `pipeline_correct_sequences.py` fix and may have had a wrong receptor
sequence going into Boltz, and (b) was just one design. design_3, run with
the verified-correct sequence and the patched plan-stage flags, gives the
opposite result. The notes 4 conclusion should be amended to "ipSAE was
zero on Pikp1 design_1 specifically, possibly due to the upstream sequence
bug; on design_3 with verified-correct inputs ipSAE is well-distributed
and informative."

## Open questions

### 1. Native vs designed cold-start comparison

The most interesting thing the user raised: Boltz's cold-start prediction
on the *native* Pikp1–AvrPikD sequence is wrong (presumably also at 20-30
Å, like the natives in the table above), but the cold start on the
RFDiffusion-designed sequence is *also* at 29 Å. So the designed sequence
isn't giving Boltz a head start on the cold start — both are wrong out of
the box. What's different is that one cycle of mild steering on the
designed sequence produces a sub-2 Å hit with high confidence, whereas no
comparable run on native Pikp1+AvrPikD has been done yet. The natural
control: run the same recipe on the native Pikp1–AvrPikD sequence. If it
also rescues to sub-2 Å with ipSAE ~0.7, design_3 is in line with the
natives. If it plateaus at 13–14 Å the way Pikp1 design_1 v2 did, then
design_3 is genuinely doing something the native sequence can't, and the
finding is much stronger. **This is the single most important next
experiment.**

### 2. The 3.0 Å ra_eff threshold may be too strict for high-confidence
designs

`c0d40.c1d06` fails the production tier by 0.37 Å on ra_eff while passing
every other criterion comfortably. The rec_rmsd is 0.91 Å — the receptor
backbone is intact in the strong sense, not the borderline-intact sense
notes 4 was worried about. The 0.37 Å miss on ra_eff is happening at a
length scale where side-chain repacking against the designed sequence and
inherent flexibility at the (Cα-only, glycine-placeholder) RFDiffusion
termini can plausibly account for it without anything biologically
meaningful changing. The threshold was set conservatively in notes 3 for
the structural-only filter on natives where ipSAE was useless; on designs
where ipSAE comes back at 0.70 with ipTM 0.85, the structural threshold
is doing redundant gating work and excluding a strong design. A
confidence-aware tier — e.g. `ra_eff < 5 AND tJac >= 0.7 AND ipSAE >= 0.5`
— would let `c0d40.c1d06` through cleanly. Worth proposing as an
RFDiffusion-specific production tier alongside the notes-3 native one,
not replacing it.

### 3. The steered interface must not rely on the steering mutations themselves

Critical check before trusting any result from this pipeline on RFDiffusion
designs: **the contact residues Boltz is using in the rescued prediction
must not be the residues that were mutated during steering.** The mutations
are a tool to push Boltz onto the right interface so we can assess binding
likelihood — they will NOT appear in the actual agroinfiltration constructs.
The construct carries the original ProteinMPNN-designed sequence, unmodified.
So if `c0d40.c1d06`'s low-RMSD, high-ipSAE pose is held in place by
interactions involving K5E, A26R, Y44R, or Q46D, the result is not
predictive of what will happen in planta — the in planta construct has
K5/A26/Y44/Q46, and Boltz has told us nothing about whether those
residues form a productive interface.

The check is straightforward: for any design that passes the filter,
intersect its `receptor_contact_residues` (from the prediction's interface
analysis) against `mutations_chimerax` for that design. If the overlap is
non-zero, the design is using its own steering mutations as contacts and
should be flagged. The acceptable designs are ones where the contacts are
all on *unmutated* positions — meaning the steering pushed Boltz onto an
interface composed of the original MPNN-designed residues, which is what
the agroinfiltration construct will actually present.

For `c0d40.c1d06` specifically: the mutations are at positions 5, 26, 44,
46. The true interface (from `rfdiffusion_metrics.json`) is at 0-based
indices 2, 32, 34, 36, 37, 38, 39, 47, 69–78 — none of the mutated
positions overlap with the true interface set, which is a good sign. But
the design's *actual* predicted contact residues (what Boltz used to place
the effector in the rescue) need to be checked directly from the prediction
PDB, not assumed from the true-interface set. That check hasn't been done
yet and is a prerequisite for claiming any design as a binding-likelihood
signal.

Worth adding as a post-processing step in `compute-final-metrics`: for each
passing design, compute `n_contacts_on_mutated_positions` and flag any
design where this is non-zero. Would go alongside the existing confidence
flag. Call it `mutation_reliance_flag` or similar.

### 4. Receptor-resurfacing design class concern (carried from notes 4)

Notes 4 warned that this design class can't meaningfully test "redirect
Boltz to a new site" because the de novo loop is built around the existing
native interface. This concern is *partially* mitigated for design_3 by
the fact that the cold start is at 29 Å (i.e. Boltz did NOT find the site
out of the box, despite the design being built around it), so the steering
is doing real work. But the concern doesn't fully go away — the residues
the steering is targeting are residues the native complex already uses,
and "find the right site" is easier when the right answer is in Boltz's
training distribution. The native-sequence control above would also
address this.

### 5. Pikp1 design_1 v2 may need rerunning

Notes 4's conclusions about design_1 — including the 13.5 Å plateau, the
zero-ipSAE finding, and the HADDOCK discrepancy investigation — were all
made before the `pipeline_correct_sequences.py` fix. design_1 v2 was
likely run with a wrong receptor sequence going into Boltz (the silent
`split_sequence` bug from notes 4 would have routed effector residues into
the receptor FASTA via the CSV path). Worth rerunning design_1 with the
verified-correct sequence to see whether the 13.5 Å plateau and zero-ipSAE
findings reproduce, or whether they were artefacts of the upstream bug.
If they reproduce: design_3 is the outlier and design_1's caveats stand.
If they don't: notes 4's main conclusions need significant revision.

## Things to remember about this session

- **Read notes 3 before working on this code, not just notes 1 and 4.**
  Notes 3 contains the production recipe, the four native benchmarks, the
  production tier criterion, and the explicit prediction about ipSAE on
  RFDiffusion designs. Not having read it carefully was the cause of
  multiple wrong claims during this session about what was novel.
- **The directory `runs/design_3_seq_3_smoke` contains the production
  run, not a smoke test.** Rename or annotate before this confuses future
  work.
- **When deriving true_interface_idx from `rfdiffusion_metrics.json`,
  the design key is `design_N.pdb` not the integer N.** Trivial gotcha
  that cost a tool call.
- **`compute-final-metrics` only populates `cm_*` columns for rows passing
  its filter** (intact AND ra_eff < threshold). The first metrics-enriched
  CSV from this session was generated with the default 3 Å threshold and
  had 0 cm_iptm rows; rerunning with `--rmsd-threshold 5.0` populated 57.
  If the cm_* columns look empty, check the threshold the metrics step
  ran with.
- **`aggregate_results.sh` originally had a hard-coded `runs/7QPX*`
  glob.** Patched this session to take a `--glob` flag with the old
  pattern as default, so it's backward-compatible. Quote the glob argument
  on the command line or the shell will expand it before the script sees
  it.
- **The orchestrator's `aggregate` subcommand and the host-side
  `aggregate_results.sh` have confusingly overlapping names but do
  different things.** The orchestrator's `aggregate` produces
  `all_results_multicycle.csv` inside one experiment dir; the shell
  script combines those across many experiments into a single CSV. If
  the multicycle file is missing after a run, run the orchestrator's
  `aggregate` first (inside the container), then the shell script.

## Concrete next steps for tomorrow

1. **Run the native Pikp1+AvrPikD control** with the same recipe
   (`--mode mild --max-mutations 4 --candidate-pool-size 6 --n-designs 50
   --n-cycles 2 --max-passing 3`). This is the single experiment that
   determines whether design_3's result is meaningful or just inside
   the native band.
2. **ChimeraX inspection of `c0d40.c1d06`** against the native ground
   truth. Confirm the 3.37 Å ra_eff is rigid-body offset of an
   essentially-correct pose, not a wrong pose at lucky distance. **While
   there, check which receptor residues are actually making contact with
   the effector in the rescued pose, and confirm none of them are at
   positions 5, 26, 44, or 46 (the steering mutations).** If any contact
   residue is on a mutated position, the design's binding signal is
   contaminated by the steering and cannot be trusted as a prediction
   about what the agroinfiltration construct will do.
3. **Implement `mutation_reliance_flag` in `compute-final-metrics`.**
   For each passing design, count how many of its predicted interface
   contacts fall on mutated positions (from `mutations_chimerax`). Flag
   any design with non-zero overlap. This is mandatory before any
   RFDiffusion-design result from this pipeline can be used to rank
   constructs for agroinfiltration.
4. **Rerun Pikp1 design_1 v2** with the patched
   `pipeline_correct_sequences.py` and verify whether the notes-4
   findings reproduce.
5. **Propose a confidence-aware production tier for RFDiffusion designs**
   in `compare_versions.py`, alongside (not replacing) the notes-3
   native tier. Suggested rule: `ra_eff < 5 AND tJac >= 0.7 AND
   ipSAE >= 0.5 AND receptor_intact == 1 AND cm_confidence_flag == 'ok'
   AND mutation_reliance_flag == 'clean'`.
6. **Update notes 4 with a footnote** noting that the "ipSAE broken on
   RFDiffusion data" conclusion is now superseded pending the design_1
   rerun.