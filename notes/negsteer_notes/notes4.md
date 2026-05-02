# Negative Steering on RFDiffusion Designs — Session Summary

## Context

Entering this session, the negative-steering pipeline had been validated on
four native benchmark complexes (notes 3). The open question for production
deployment was whether the same workflow could be run on **RFDiffusion-designed
binders** with **ProteinMPNN-designed sequences**, where the "ground truth"
is a Cα-only design backbone rather than a solved native complex.

Three structural mismatches against the existing pipeline made this
non-trivial:

1. RFDiffusion outputs are Cα-only — no side chains. The negative-steering
   pipeline's `find_contact_residues_heavy` requires heavy atoms to compute
   the protected interface set on the ground truth.
2. The de novo region of an RFDiffusion design has glycine placeholders in
   the PDB but a real designed sequence in the corresponding ProteinMPNN
   FASTA. The two have to be reconciled before Boltz can fold anything.
3. The receptor sequence Boltz folds (the ProteinMPNN design) is not the
   same as the sequence the ground-truth PDB encodes (glycines in the
   placeholder region). The existing `--receptor-fasta` override mechanism
   was added in the notes 2 session for sequence-walk diagnostics; it
   needed to be re-purposed here.

The session worked through these mismatches one by one, hit and fixed two
pipeline bugs along the way, ran an end-to-end validation, and surfaced an
unexpected finding about HADDOCK's role upstream of negative steering.

## What we built

### 1. `cmd_plan` patch: pre-computed true-interface indices

Added two new flags to `boltz2_negative_steering.py plan`:

- `--true-interface-indices "i1,i2,..."` — comma-separated list of 0-based
  receptor seq_idx values that constitute the protected (true) interface
- `--true-interface-indices-file PATH` — same data, read from a file with
  permissive whitespace/comma/newline separation and `#` comments

When either flag is set, `cmd_plan` skips the `find_contact_residues_heavy`
call on the ground truth entirely and constructs `true_contacts` from the
pre-computed list with `NaN` distances. The plan stage records
`true_interface_source: "pre_computed"` (or `"heavy_atom_contact"` on the
default path) so any future diagnostic run can tell which path was taken.

A new helper `_load_true_interface_indices` validates every index against
the receptor length and rejects out-of-range or non-integer input loudly.
Patched 11 unit tests covering CLI parsing, file parsing, deduplication,
and error cases all pass.

The override-fasta mechanism from notes 2 handles the sequence side: the
patched workflow uses `--receptor-fasta` to make Boltz fold the ProteinMPNN
sequence while the ground-truth PDB (with its glycine de novo region) is
used only for the Cα Kabsch alignment in `binding_rmsds`. The two are kept
length-consistent by `_load_seq_override`.

### 2. `split_sequence` bug fix in the receptor-resurfacing pipeline

While building the input files for the test run, we discovered that
`pipeline_correct_sequences.py` in the upstream RFDiffusion → ProteinMPNN
pipeline was producing systematically wrong `corrected_receptor` sequences
in `sequence_metadata.csv` and `af2_fastas/`. The investigation tracked it
down to `split_sequence`, which had a hard-coded assumption that the chain
order in MPNN's output was `effector / receptor`. The actual order in this
pipeline (where the filter writes split PDBs as receptor=A, effector=B) is
the opposite: `receptor / effector`.

The result was that every "corrected receptor" sequence written to the CSV
was actually the *effector* residues, chopped to fit the receptor contig
structure. For Pikp1+AvrPikF this meant every CSV row had
`denovo_length=16` and a designed region of `QDAGTCYSSWKMDKKV` (residues
33–48 of AvrPikF) regardless of which RFDiffusion design it was for. The
failure was completely silent — every downstream consumer (the AF2 FASTAs,
the QC step, the design-region scoring) was operating on garbage receptor
sequences.

The fix: rewrite `split_sequence` to pick whichever half of the
slash-separated MPNN output is closest in length to the known
`effector_len`, and assign that as the effector. Robust against either
chain ordering. The redundant "swap if lengths look wrong" block in
`correct_sequences` (which used a `> 5` threshold and missed cases where
the length difference was 3) was removed since the new `split_sequence`
handles swapping internally. Patched file produced as a drop-in
replacement; behavioural tests pass on receptor-first slash, effector-first
slash, colon separator, and equal-length ambiguous inputs.

The bug had been silently affecting every previous receptor-resurfacing run.
The actual MPNN-designed receptor sequences exist in the raw
`design_N_mpnn/seqs/design_N.fa` files and were never propagated to the
CSV correctly until the fix.

### 3. Deriving `true_interface_idx` from `rfdiffusion_metrics.json`

The RFDiffusion filter (`rfdiffusion_filter.py`) computes per-design Cα–Cα
contacts at 8 Å between the receptor and effector and writes them to the
metrics JSON as `receptor_contact_residues`. The list uses the filter's own
synthetic resnum scheme, which differs from the contiguous 1..N numbering
in the split PDB. Both schemes describe the same residues in the same
file order — the difference is purely cosmetic labelling.

The mapping from contact resnums to 0-based seq_idx values is one
dictionary lookup against `receptor_position_order` from the same JSON
entry. No PDB parsing needed. The derivation is self-contained in the JSON.

For Pikp1 design_1 in the dr1 run, this produced
`true_interface_idx = [34, 35, 36, 37, 38, 39, 40, 41, 81, 82, 83]` —
11 contacts: 8 in the de novo loop and 3 in the C-terminal end of the
fixed region.

### 4. The 8 Å Cα cutoff caveat

The RFDiffusion filter uses 8 Å Cα–Cα for contact detection (because no
side chains exist). The negative-steering pipeline normally uses 4.5 Å
heavy-atom. These are different criteria and the resulting interface
residue sets are not directly comparable to the notes 3 native benchmark.
The 8 Å Cα set is approximately a superset of the 4.5 Å heavy-atom set,
which means using it as the protected set is *more* conservative — Claude
notes 3's 0.7 true_jaccard threshold may run high on RFDiffusion data, but
in the v2 run this turned out to be moot for other reasons (see below).

## Experiments and results

### Smoke test — Pikp1 design_1, n_designs=4, n_cycles=1

Successful end-to-end run. The patched `cmd_plan` accepted the pre-computed
indices, stored them in `plan.json`, the override FASTA loaded correctly
(84 residues matching the PDB chain A length), and the array job ran 4
Boltz predictions to completion.

Cold start: ra_eff 32.54 Å, intact 0.
Best steered (design_02): ra_eff 17.13 Å, rec_rmsd 3.89, eff_rmsd 1.42, intact 0.
All 4 designs failed `receptor_intact == 1` (rec_rmsd 2.82–6.12 Å).
All 4 reduced ra_eff vs cold start (Δ between 1.98 and 15.42 Å).

ChimeraX inspection of design_02 showed the receptor core β-sheet and
α-helix were preserved; deviations were concentrated at the disordered
termini and immediately around the 6 mutated positions. Not a global fold
failure — the 3.0 Å intact threshold is overly strict for RFDiffusion
designs whose termini are inherently floppy.

The cold-start `independent_effector_rmsd` of 8.80 Å was traced (via
ChimeraX) to a disordered C-terminal tail in the effector that adopted
different conformations between Boltz's prediction and the RFDiffusion
ground truth. The β-sheet core overlaid cleanly. Same threshold-strictness
issue as the receptor termini, on the effector side.

### v2 — Pikp1 design_1, n_designs=50, n_cycles=2, max_passing=3

Configuration: `--mode mild --max-mutations 4 --candidate-pool-size 6
--n-designs 50 --n-cycles 2 --max-passing 3`. Reduced max_mutations from 6
to 4 to allow positional variation across designs (with 6 candidates and 6
mutations all designs hit identical positions; with 4 they cover 15
possible position subsets).

Cycle 0: 50 designs, ra_eff range 14.6–41.0 Å, median 25.2.
Cycle 1: 150 designs across 3 pathways (parents c0d00, c0d20, c0d33),
ra_eff range 13.5–45.4 Å, median 25.7.

Best result across the whole run: `c0d00.c1d18` at ra_eff 13.49 Å,
rec_rmsd 2.96, eff_rmsd 0.72, intact 1, true_jaccard 0.39.

**Hard floor at ~13.5 Å.** Across 200 mutated predictions and 2 cycles, no
design got within 10 Å of the HADDOCK ground truth. Cycle 1 improved the
floor by exactly 1 Å over cycle 0 (14.6 → 13.5). The candidate pool was
fully exhausted in cycle 0, and cycle 1's selection of diverse parents
(c0d00, c0d20, c0d33 — *not* the best ra_eff designs) by maximin distance
explored alternative starting points but converged on the same general
ra_eff band.

### The HADDOCK discrepancy

The headline finding emerged when the user noticed that the
negative-steered predictions, despite having ra_eff ~14 Å against the
HADDOCK pose used as ground truth, appeared visually closer to the **native
6G10 effector pose** than to HADDOCK's. ChimeraX measurements:

| comparison | distance |
|---|---|
| HADDOCK ground truth ↔ Boltz cold start | 32.5 Å |
| HADDOCK ground truth ↔ best steered prediction | 13.5 Å |
| Native 6G10 ↔ best steered prediction | **10.2 Å** |
| Native 6G10 ↔ Boltz cold start | **26.4 Å** |

The four numbers form a triangle. Boltz's cold start is 26 Å from native
and 33 Å from HADDOCK — far from both, in a third location. The steering
pulled it toward HADDOCK and got partway (13.5 Å), and that trajectory
incidentally passed through territory closer to native (10.2 Å) than to
the actual HADDOCK target.

#### Interpretation

The 26.4 Å cold start vs native rules out the hypothesis "Boltz natively
finds the native pose." Boltz fails the cold start on this complex the
same way it fails the cold start on every notes-3 native benchmark
complex — at distances of 24–32 Å. There's nothing special about Boltz's
preferences here.

The user pointed out the more fundamental issue: in this design class
(receptor resurfacing of HMA-effector binding), the RFDiffusion design
region is built **around the existing native interface**. The de novo loop
reshapes residues that were already at the binding site. AvrPikF and the
original Pik effector are nearly identical, so the interface isn't being
moved — it's being decorated. Negative steering can't meaningfully test
"can mutations push Boltz away from the original site toward a designed
site" because the original site and the designed site are the same place.

The 10.2 Å distance from native to the steered prediction is geometrically
incidental: HADDOCK's pose, the native pose, and the steered cluster are
all on the same general face of the receptor, so any trajectory between
them passes through nearby territory.

#### The actual root cause: HADDOCK

The HADDOCK step in the upstream receptor-resurfacing pipeline placed the
effector in a wrong rigid-body orientation despite being constrained to
keep it at the canonical binding site. The pipeline then built the entire
RFDiffusion design around that wrong placement. The negative-steering run
was therefore being asked to optimise Boltz toward a target that wasn't
biologically meaningful in the first place.

**Success in the negative-steering metric, on this design, would have
corresponded to a worse biological prediction.** A "production-tier
rescue" with ra_eff < 3 Å against the HADDOCK pose would have meant
Boltz had been successfully pushed *away* from the native-like region
toward the HADDOCK error. This is a critical caveat for any future
deployment of the workflow on receptor-resurfacing outputs.

## Key findings

### 1. The workflow is validated for RFDiffusion designs

The patched negative-steering pipeline runs end-to-end on a Cα-only
RFDiffusion ground truth with a ProteinMPNN sequence. The patches work,
the pre-computed indices flow through the plan and predict-one stages
correctly, and `aggregate` + `compute-final-metrics` produce the same
output shape as on native complexes. The infrastructure question that
opened the session is closed.

### 2. ipSAE is uniformly zero on RFDiffusion data too

Notes 3 left open the question of whether RFDiffusion designs (which span
wider sequence space than native variants) might produce more meaningful
ipSAE values than the native benchmark complexes. A subset of v2 designs
had `cm_iptm` populated by an early `compute-final-metrics` pass. All of
them — including the best-ranked designs by ra_eff — had `cm_ipsae_min =
0.000`, `cm_iptm` 0.08–0.10, and `cm_pae_pass_frac = 0.0`.

This **disconfirms** the notes 3 hypothesis. RFDiffusion designs don't
unlock confidence signal. The production tier must remain
**structural-only**, not confidence-gated. ipSAE is not a portable ranker
for any of the regimes the pipeline has now been tested on.

### 3. The negative-steering validator can't test receptor-resurfacing designs

Where the de novo region rebuilds an existing interface rather than
creating a new one, the conserved residues that dominate Boltz's binding
prediction are unchanged. The mutated positions are decorations on a fixed
binding mode. Negative steering can't push the effector to a "new"
binding site because there isn't one — there's just the existing site,
maybe with a small rigid-body offset.

This is a property of the design class, not a failure of the workflow.
The interesting test for the validator is a design where the de novo
region creates a binding site spatially distinct from any binding mode
Boltz would natively predict. Receptor-resurfacing designs by construction
don't satisfy that criterion.

### 4. HADDOCK is the weak link in the upstream pipeline

For Pikp1+AvrPikF specifically, HADDOCK appears to have placed the effector
in a wrong orientation despite explicit constraints. This propagated
forward through RFDiffusion, ProteinMPNN, and into the negative-steering
ground truth, where it caused 200 Boltz predictions to chase an incorrect
target. The receptor-resurfacing pipeline cannot produce trustworthy
designs without first establishing that its docking step produces
biologically sensible starting poses. This is a separate workstream from
the negative-steering pipeline itself.

### 5. The 3 Å receptor intact threshold is too strict for RFDiffusion designs

Notes 3 already flagged that ra_eff 3.0 Å might need loosening for
RFDiffusion data. The same is true of `RECEPTOR_INTACT_CUTOFF`. RFDiffusion
designs have inherently floppy termini (no anchoring residues outside the
modeled chain) and the 3.0 Å threshold rejects designs whose β-sheet cores
are perfectly preserved but whose termini have wandered. A core-only or
secondary-structure-only RMSD would be more meaningful.

## File deliverables from this session

- `boltz2_negative_steering.py` — patched with `--true-interface-indices`
  and `--true-interface-indices-file` flags, the `_load_true_interface_indices`
  helper, and a `true_interface_source` field in `plan.json`
- `pipeline_correct_sequences.py` — drop-in replacement for the
  receptor-resurfacing pipeline's `bin/pipeline_correct_sequences.py`,
  with the `split_sequence` chain-order bug fixed and the redundant
  swap-block in `correct_sequences` removed
- `design_1_seq_4_receptor.fasta` — the manually-extracted MPNN sequence
  for Pikp1 design_1 sample 4 (mpnn_score 0.9484, 84 residues), used as
  the `--receptor-fasta` for the validation runs
- Two complete v2 experiment outputs in
  `runs/Pikp1_design_1_test/` (4-design smoke test) and
  `runs/Pikp1_design_1_v2/` (50-design × 2-cycle production-style run)

## Status summary

- **Negative-steering workflow validated for RFDiffusion + ProteinMPNN
  inputs** — patches work, output is sensible, end-to-end pipeline runs
  cleanly
- **ipSAE confirmed broken on RFDiffusion data** — production filter must
  remain structural-only
- **Receptor-resurfacing designs are not the right test case** for the
  validator's "redirect to new binding site" capability — this is a
  design-class limitation, not a workflow bug
- **HADDOCK in the upstream pipeline is the immediate blocker** for
  meaningful negative-steering results on the current Pikp1+AvrPikF
  design — needs investigation tomorrow
- **Two RMSD threshold issues** (effector tail, receptor termini) deferred
  for future work as Option C

## Next steps

### Tomorrow — HADDOCK day

1. **HADDOCK positive controls.** Run the pipeline's HADDOCK step on a
   complex with a known answer (e.g. 6G10 itself, with effector and
   receptor stripped to sequences) and measure how close the docked output
   is to the deposited structure. Tells us whether HADDOCK is broken in
   general or only on Pikp1+AvrPikF.

2. **Inspect HADDOCK constraint generation.** Pull the AIR file the
   pipeline produces and check whether the active/passive residue
   assignments correspond to known HMA-effector interface residues. The
   most likely "just a bug" candidate.

3. **Cluster re-ranking by homology agreement.** HADDOCK outputs multiple
   clustered solutions ranked by HADDOCK score. Compute Cα RMSD of every
   cluster's effector against an external reference (6G10 for the current
   design) and pick the one closest to reference rather than the
   top-scored one. Cheapest fix for the immediate Pikp1+AvrPikF case.

4. **Tighten constraint specificity.** If the AIRs are loose (large
   passive lists, weak distance bounds), tighten them and re-run. HADDOCK
   gives the effector less freedom when constraints are stricter.

5. **HADDOCK protocol parameter sweep** — increased rigid-body sampling,
   more refinement iterations, lower temperature schedule. Default protocols
   favour speed.

### Negative-steering pipeline improvements (deferred)

6. **Option C: mask disordered tails in RMSD.** Both effector
   (`independent_effector_rmsd`) and receptor (`independent_receptor_rmsd`
   and the intact threshold). Either core-only RMSD against a manually-
   specified residue range, or automated tail-detection from B-factor or
   secondary-structure annotation. Should bring the intact flag's behaviour
   into agreement with visual fold-quality assessment on RFDiffusion designs.

7. **Switch ra_eff to an interface-geometry metric** (carried from notes
   3). DockQ, lDDT-PPI, fnat, or interface-only RMSD. Notes 3 identified
   iRMSD as the easiest drop-in. Still the right priority.

8. **Loosen `RECEPTOR_INTACT_CUTOFF` for RFDiffusion runs**, or replace it
   with a core-only RMSD. The 3.0 Å threshold rejects fold-preserved
   designs whose only deviation is at the floppy termini.

### Design-class question

9. **Find an RFDiffusion design where the validator can actually be
   tested.** Specifically, a design where the de novo region creates a
   binding site spatially distinct from any binding mode Boltz would
   natively predict. Receptor-resurfacing designs (where the existing
   interface is decorated rather than redirected) don't satisfy this.
   Possibilities: hallucinated binders against effectors with no known
   receptors, or scaffolded mini-binders against new targets. This is the
   right next experiment once HADDOCK is sorted.

### Carried forward from notes 3

10. **RFDiffusion production validator wrapper script** — iterate over a
    directory of designs, run the full pipeline per design, produce one
    summary CSV row per design. Becomes a real need once HADDOCK is fixed
    and we want to run the validator at scale.

11. **Second-shell pool tuning on 60–80% overlap complexes** — never
    tested empirically. Open since notes 3.

12. **Pose clustering within production tier** — for collapsing
    near-duplicate rescues into representatives. Open since notes 3.

13. **Benchmark chain-pair audit script** (carried from notes 2) — count
    heavy-atom contacts between every chain pair across the benchmark and
    flag entries where the user-specified pair has fewer contacts than the
    best pair. A few hours of work, never done.

### Filed but not fixed

14. **Region 2 length bug in RFDiffusion contig handling.** The
    `denovo_length=4 or 5 instead of 6` issue from the original
    7QPX/design_21 session. Affects designs with multiple variable de novo
    regions. The single-region contig used in dr1 sidesteps it. Worth
    filing upstream once HADDOCK is sorted.

### Deprecated or moot

15. Mixed-sequence 6G10/7B1I experiment (deprecated in notes 3, still
    deprecated).
16. Cross-version confidence gating on borderline tier (moot — ipSAE
    confirmed uniformly zero on RFDiffusion data, see Key Finding 2).
