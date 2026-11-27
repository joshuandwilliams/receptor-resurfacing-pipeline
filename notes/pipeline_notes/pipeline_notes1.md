# Receptor-Resurfacing Pipeline Cleanup Pass — Session Summary

## Context

Entering this session, the receptor-resurfacing Nextflow pipeline had been
running for some time and had accumulated cruft from earlier development
phases. An unused AF2-monomer / FASTA-input branch survived from before
HADDOCK was wired in. Container path strings were duplicated across the main
pipeline and five per-module test scripts. Several Python scripts carried
copy-paste comments that no longer matched the code. At least one silent bug
class around design-region length reporting was already known, and the user
had been quietly working around it.

The plan was a per-block cleanup pass: read each block's `.nf` module, audit
its Python scripts, fix anything broken or misleading, and run the per-module
test scripts to validate. Strict ground rule: **minimal targeted changes
only — no functional changes, just cleanup**. Earlier interactions had
established that overengineered "while we're here" refactors get rejected
hard, so the calibration was diff-minimisation as a first principle.

The session worked through preprocessing, HADDOCK, RFDiffusion, and Rosetta
filtering in order. RFDiffusion became a much deeper audit than expected
when a multi-denovo test exposed three real bugs the simpler tests had hidden.
Rosetta was small and clean. ProteinMPNN, MSA, Boltz2, and the aggregate block
remain for a future session, with a substantial deferred-work list to carry
forward.

## What we built

### 1. AF2 monomer / FASTA input path removed

The pipeline originally supported a third entry point alongside the
HADDOCK-from-two-PDBs path and the pre-docked-complex-PDB path: a "give me
two FASTAs and I'll AF2-fold them" branch. It hadn't been used in months
and was the source of several stale code paths and config keys. Archived
verbatim to `archive/af2_monomer_fasta_input/`. The archive holds the
entire `af2_monomer.nf` plus a `removed_blocks.md` listing every removed
snippet from `main.nf`, `preprocessing.nf`, the slurm wrapper, and
`params_example.yml`. A README explains the removal rationale and the
resurrection path. The pipeline is now strictly PDB-only, either two
PDBs into HADDOCK or one complex PDB straight to RFDiffusion.

### 2. Container path consolidation

Four container paths (`rfdiff_container`, `rosetta_container`, `boltz2_container`,
`colabfold_container`) were previously declared in `main.nf` AND in all five
`test_*.nf` scripts. At least three of the test scripts pointed at stale
subdirectory layouts. They used `RFDiffusion/Rosetta.img` instead of
`HADDOCK_RFDiffusion_ProteinMPNN_MMseqs2/Rosetta.img`, and `Boltz2/boltz2.img`
which didn't exist on disk. All four params now live exclusively in
`nextflow.config` as the single source of truth, inherited automatically
by every test. Touched: `main.nf`, `nextflow.config`, all five
`test_*.nf` files.

### 3. Preprocessing block fixes (3 bugs)

- **Bug 1: dummy mapping byte-count off-by-one.** `WRITE_DUMMY_MAPPING`
  was using `echo '{}'`, which writes 3 bytes (`{`, `}`, `\n`). Changed
  to `printf '{}'`, which writes 2 bytes. The downstream consumer had a
  size check that misbehaved on the 3-byte version on certain filesystems.
- **Bug 2(a): `receptor_start_pdb` plumbing into sequences.json.**
  `EXTRACT_SEQUENCES` now writes `receptor_start_pdb` to its JSON output,
  matching what `BUILD_CONTIGS` already does for `updated_params.json`.
  The field was previously absent and downstream consumers fell back to
  a default of 1, silently breaking PDBs that started at non-1 residue
  numbers. (The MPNN-side consumers still default to 1, see deferred
  Bug 2(b) below.)
- **Bug 3: stale docstring.** `EXTRACT_SEQUENCES`'s docstring still
  described the old AF2-monomer behaviour. Updated to match what it
  actually does.

### 4. HADDOCK block — 12 fixes applied, test validated

The HADDOCK audit produced a long fix list. The headline items:

- **H6 (high severity, recurring across the codebase):** the substring
  pattern `chain.upper() in block.upper()` matched a contig block to a
  chain ID. It gives false positives when chain letters are multi-character
  or share characters with other tokens in the block. Replaced with a
  per-segment chain-prefix check. The same pattern surfaced again later in
  RFDiffusion (R1, R2) and was fixed identically.
- **H7:** `HADDOCK3_PLOTS` test was being called with 5 arguments where
  the process expected 6. Caught by reading the test file against the
  process declaration. It would have failed at runtime on first invocation.
- **H10 (high severity, hard to find):** chain disambiguation in
  `extract_hotspots.py` was based purely on chain length, which fails
  when receptor and effector lengths are similar. Added a sequence-identity
  fallback using a stdlib sliding-window helper `_best_identity` and two
  new CLI flags (`--receptor-seq`, `--effector-seq`). This is the same
  class of bug as the `split_sequence` chain-order bug from the previous
  session. Length-based disambiguation is fragile when the two chains
  are within ~5 residues of each other.
- **H15 (documentation only):** the hardcoded `min_population=4`
  threshold appears in four locations across the HADDOCK block
  (`haddock.nf`, `collect_haddock3_dock.py`, `haddock3_plots.py`,
  `haddock_utils.py`). Added cross-referencing comments at all four
  locations so a future change in one place flags the others. Proper
  consolidation deferred to ProteinMPNN block where the param-passing
  pattern is clearer.
- Plus seven smaller fixes covering docstrings, error messages, test
  arity, and a dropped unused emit.

The HADDOCK test ran cleanly after the fixes: 7/7 processes succeeded,
13m 27s wall, against PikP1HMA + AVRPikF.

### 5. RFDiffusion block — 9 fixes including three high-severity bugs surfaced by multi-denovo testing

This block was the deepest audit of the session. The first test run
(single denovo region, contig `B A1-390/20-40/A421-438`) passed cleanly
with 4 designs producing sensible metrics. The user then ran a
**multi-denovo test with a ChimeraX-built complex PDB** (chains B and C
instead of A and B, contig `B187-208/10-20/B223-262/6-6 C`). Three new
bugs surfaced in rapid succession:

**R9 (defensive fix, surfaced by user test):** the user's first multi-denovo
test had the wrong chain letters in the contig (referencing `A` against a
PDB with chains `B` and `C`). RFDiffusion's parser failed 5 minutes into
GPU spin-up with a confusing assertion `('A', np.int64(1)) is not in pdb
file!`. Added a guard in `contig_utils.remap_segments_to_pdb` for fixed
segments that reference a chain with no atoms in the PDB. It fails loudly
with a clear error naming the chain and the referenced segments. The
re-run with the same broken contig now fails in seconds with a useful
message instead of waiting for RFDiffusion to load and crash.

**R10 (high severity, latent for months):** with the chain letters fixed,
the next run produced wildly broken metrics. There were 779 contact pairs
(impossible for an 81-residue interface), 105 Å design-region RMSDs, and
three of four designs reporting identical contact counts. The cause was in
`rfdiffusion_filter.py:main()`, where the per-design analysis loop had a
chain-ID-based shortcut that assumed the user's input chain letters
(`receptor_chain="B"`, `effector_chain="C"`) survived RFDiffusion's
processing. They don't. RFDiffusion always renames output chains to A/B
in the order they appear in the contig, regardless of input. So
`args.receptor_chain="B"` looked up RFDiffusion's chain B, which was the
*effector*. And `args.effector_chain="C"` wasn't in design_chains at all,
so it fell through to `design_chains[1]="B"`. Both pointed at the same
chain, and `find_contacts(rec_atoms, eff_atoms)` was computing contacts of
a chain with itself.

The bug had been latent forever because the user's earlier successful
tests all happened to use `receptor_chain="A"` and `effector_chain="B"`,
which match RFDiffusion's output convention by coincidence. The fix:
delete the chain-ID shortcut entirely, always use
`split_at_chain_break + identify_segments`, which works on chain-break
detection and length matching rather than chain letters. Validated on
re-run: receptor sizes 78–86 residues (matching contig spec), contacts
35–39 per design, designs genuinely distinct.

**R7 (medium severity, multi-denovo only):** with chain identification
correct, the per-region length reporting was still wrong. Designs showed
`per_region_lengths = [18, 5]`, `[13, 4]`, `[12, 4]`, `[18, 6]` for a
contig spec of `10-20` and `6-6`. Three of four designs reported region
2 as fewer than 6 residues, violating the contig's hard `6-6` constraint.
The bug was in `build_receptor_resnum_map`'s allocation logic for the
multi-denovo branch. It used a purely proportional formula
(`alloc = round(total_denovo * desc[2] / total_max)`) clamped to
`[min_len, max_len]`, then re-clamped against the remaining budget. That
second clamp could drop a segment below its minimum. The final segment
got whatever was left over with no bounds check at all.

Fix strategy: any segment with `min_len == max_len` (a fixed-length
denovo region like `6-6`) gets exactly that length first. The leftover
budget is distributed to variable-length segments only. If there are
multiple variable-length segments, the leftover is distributed
proportionally with min/max clamping, with the final segment absorbing
the rounding remainder. Validated on re-run: every design now reports
region 2 as exactly 6 residues, region 1 in `[10, 20]`, math adds up.

The corresponding pattern in `pipeline_correct_sequences.py` (~lines
280–320 in `align_to_native_by_anchors`) has the same bug and needs the
same fix back-applied — deferred to MPNN block.

**R6 (medium severity, latent footgun):** `plot_design_lengths` had a
fallback path for designs whose metrics lacked `per_region_lengths`. It
called `_split_design_into_regions`, which uses gap-based sorting on
residue numbers. But the de novo residues in the metrics JSON have
*synthetic* resnums starting at `max_fixed_resnum + 1000` to avoid
collisions with real fixed residues. Sorting them puts all de novo
residues into one giant block at offset +1000, collapsing N regions into
1 with wildly wrong sizes. The primary path uses `per_region_lengths`
directly, and it is correct and is what gets used today. But if any
future caller forgot to write `per_region_lengths`, the fallback would
silently produce nonsense plots. The fallback now refuses to plot and
emits a clear warning instead, so the next person to hit this sees it
rather than a silently misleading plot.

**R12 (deferred):** the `design_region_rmsd` field uses positional
alignment without Kabsch superposition. For the multi-denovo test it
reports values around 124 Å, which is consistent with "two coordinate
sets dropped into the same frame with no rotation/translation". Confirmed
this is *not* a filter gate. The filter passes/fails on
`frac_contacts_in_design`, not RMSD, so it doesn't affect pipeline
behaviour, only the reported RMSD value. Proper fix requires Kabsch
superposition on the fixed segments before RMSD-ing the design region,
which is ~30 lines of careful code. Deferred to ProteinMPNN block where
there'll be more context to do it properly.

Plus four smaller R1–R5/R8 fixes. They cover the H6-class chain-substring
pattern (R1, R2), a dropped misleading comment (R3), dead code removal
(R8), and a documentation note about RFDiffusion's two coexisting
residue numbering schemes. The metrics JSON uses contig-space numbers
matching the input PDB. The split PDBs use 1-based output-space numbers
that always start from 1, which was undocumented previously and is
exactly the kind of thing that bites later.

### 6. Rosetta filtering block — 5 fixes, small and self-contained

This block was clean by comparison. Total 715 lines across all five
files. Five fixes:

- **Ros1:** `rosetta_filter_collect.py` had cryptic NaN checks using the
  `value != value` trick repeated 7 times. Replaced with `math.isnan()`
  and two helper functions (`_round_or_none`, `_int_or_none`).
- **Ros2 (real bug):** `rosetta_filter_plots.py:plot_sc_vs_dg` looked up
  the best-model name with `d.get("design_name", d.get("model", f"design_{i}"))`.
  The JSON written by the collect script uses field names `"design"`
  and `"design_stem"` instead. The fallback always fired, so the "best
  model" star marker was always labelled with a positional index instead
  of the actual filename. Fixed.
- **Ros3:** the `.nf` module's input comment said the channel held
  "flattened tuples" consumed positionally, but the Python script
  actually globs the work dir for files. Comment-only fix.
- **Ros4:** `axvspan` was being called *before* the scatter, so
  `ax.get_xlim()` returned matplotlib's default `(0, 1)` instead of
  data-driven limits. Worked by coincidence for typical Sc data
  (0.4–0.8) but brittle. Moved `axvspan` after the scatter.
- **Ros5:** function docstring claimed `"sc_value (y) vs dG_separated (x)"`
  but the code (and axis labels) put Sc on x and dG on y. Fixed the
  docstring to match.

The Rosetta test ran cleanly on real RFDiffusion split PDBs from the
previous block's output. Sc values were 0.474–0.501, realistic for
unrelaxed polyvaline backbones, and 1/4 designs passed at the Sc≥0.5
threshold. `nres_int` varied per design (67–70), confirming R10's chain
identification fix propagates through to Rosetta cleanly.

The metrics had positive `dG_separated` values (+62 to +85 REU) and
`packstat=0.0` for every design. This initially looked wrong but makes
sense. RFDiffusion outputs are polyvaline backbones with no real
sidechains, so Rosetta's energy and packing terms are dominated by
valine-clash artefacts at this stage. The pipeline correctly uses
`Sc` (purely geometric) as the filter gate rather than dG, because dG
isn't physically meaningful here. Added a docstring caveat to
`rosetta_filter_collect.py` explaining this so future readers don't
misinterpret the energy terms.

## Experiments and results

### HADDOCK module test

7/7 processes succeeded, 13m 27s wall, real test against PikP1HMA +
AVRPikF. Validated H6 (chain matching), H7 (test arity), H8 (channel
wiring), H10 (sequence-identity disambiguation), and the container
consolidation. Single test run covered all 12 fixes.

### RFDiffusion module test #1 — single denovo region

Contig `B A1-390/20-40/A421-438`, 4 designs. Clean run, ~13 min wall.
Per-region lengths 26, 27, 25, 35 (within the 20–40 spec). Receptor
totals 434/435/433/443 = 408 fixed + denovo length, math correct.
RMSDs 2.2–7.9 Å (the 7.9 outlier is the longest design region, which
is biologically plausible — more flexibility, more divergence).

### RFDiffusion module test #2 — multi-denovo, ChimeraX input

Contig `A1-32/10-20/A46-72/6-6 B` initially (wrong chain letters for
the PDB which had B and C). **R9 caught it in seconds** with a clear
error message after the fix was applied.

Re-ran with `B187-208/10-20/B223-262/6-6 C`. First post-R9 run produced
the broken metrics that exposed R10 (779 contacts, identical metrics
across designs, 105 Å RMSDs). Second post-R10 run produced sensible
contact counts but per_region_lengths still wrong (R7). Third post-R7
run was clean. Every design had `per_region_lengths` region 2 = exactly
6 and region 1 in [10, 20], with the math adding up against
`n_receptor_residues`. Contacts were in the 35–39 range varying per
design, and frac_contacts_in_design varied meaningfully. Total wall
~2 min per re-run on a small (4-design) test.

### Rosetta filtering test

4 designs from the multi-denovo RFDiffusion test, single SLURM job
per design (parallelised by Nextflow), wall time a couple of minutes.
Sc 0.474–0.501, 1/4 passes at threshold 0.5, all 7 metrics present,
no NaN/None values, varying nres_int and dSASA_int across designs
confirming chain identification correct.

## Key findings and lessons

### 1. Length-based chain disambiguation is fragile when chains are similar in size

This class of bug surfaced three times in this session. The instances
were H10 in `extract_hotspots.py`, R10 in `rfdiffusion_filter.py`, and
the `split_sequence` bug from the previous session that's referenced in
notes 4. For PikP1HMA + AVRPikF specifically, the receptor is 76
residues (HMA domain) and the effector is 81 residues. They are within 5
residues of each other, well inside the noise margin of any
length-based heuristic. The right fix is sequence-identity-based
disambiguation (like H10's `_best_identity` helper) or chain-order-based
(like R10's "always use split + identify_segments"). Length alone is
not enough.

### 2. Multi-region contigs are the stress test that surfaces a whole class of bugs

The single-denovo test passed cleanly and would have let R7, R10, and
R6 all stay latent indefinitely. The user's decision to test with a
two-region contig (`10-20/.../6-6`) and a non-A/B chain PDB (chains B
and C from ChimeraX) exposed all three. **Any future code touching
`build_receptor_resnum_map`-equivalents needs to be tested with at
least one multi-denovo contig**. Ideally that contig has both a variable
region and a fixed-length region (`6-6` style), to exercise the
allocation logic in both branches.

### 3. The "skim plotting code" calibration was wrong for the filter

Earlier in the session I made a calibration call to skim
`rfdiffusion_filter.py` (650 lines) and `rfdiffusion_plots.py` (835
lines). I used function inventories and grep patterns instead of reading
them in full, on the basis that "plotting code is rarely
correctness-critical". The user pushed back, pointing out that the
filter is the gate that decides which designs get forwarded to MPNN,
not plotting at all. The deeper read found R6, R7, R10, and R12, all
real bugs that the smoke tests had cleared. **For the MPNN block, read
every Python file in full.** Function inventories and greps are good
for surface-area checks but miss logic bugs by definition.

### 4. Fast-failing on bad inputs is high-leverage

R9 cost about 30 minutes to design and implement. It paid for itself
on the first user test that hit it. It saved ~5 minutes of GPU spin-up
time. It also converted a confusing assertion error from RFDiffusion's
internals into a clear "your contig references chain X but the PDB
has chains Y and Z" message. For a pipeline that gets run repeatedly
by humans who occasionally typo chain letters, this kind of defensive
fix has high real-world value.

### 5. The "design region length" concept is genuinely ambiguous in user-facing outputs

Throughout the session the user kept seeing CSV columns and JSON fields
labelled "design region length" with values that didn't match what
they'd asked for. The investigation revealed three different concepts
all being reported under similar names. (a) The contig spec range,
e.g. "20-40", is what was *asked for*, and the honest answer is a
range. (b) The per-design observed length is what RFDiffusion
*actually produced*. (c) Gap arithmetic on the input PDB gives the size
of the gap *in the input scaffold*, which is only meaningful for
fixed-length design tasks. The proposed rename is `design_region_spec`
for (a) and `design_region_length_observed` for (b), with (c) dropped
entirely. It is deferred to the MPNN block because it touches consumers
across multiple files and should be done atomically.

## File deliverables from this session

- `preprocessing.nf`, `main.nf`, `nextflow.config`, `params_example.yml`,
  `run_pipeline_slurm.sh` — preprocessing fixes plus AF2-monomer removal
  plus container path consolidation
- `haddock.nf`, `haddock3_prepare.py`, `haddock_utils.py`,
  `collect_haddock3_dock.py`, `haddock3_plots.py`, `extract_hotspots.py`,
  `test_haddock.nf` — 12 HADDOCK fixes
- `rfdiffusion_filter.py` — R1, R7, R8, R10 plus dual-numbering docstring
- `rfdiffusion_plots.py` — R6
- `contig_utils.py` — R2, R9
- `test_rfdiffusion.nf` — R3
- `rosetta_filtering.nf` — Ros3
- `rosetta_filter_collect.py` — Ros1 plus polyvaline-caveat docstring
- `rosetta_filter_plots.py` — Ros2, Ros4, Ros5
- `test_proteinmpnn.nf`, `test_boltz2.nf`, `test_rosetta_filtering.nf`
  — container path consolidation only
- `pipeline_correct_sequences.py` — TODO marker only, fix deferred
- `archive/af2_monomer_fasta_input/` — verbatim archive of removed code
  with README explaining resurrection path
- `deferred_to_proteinmpnn.md` — handover document for the next session

## Status summary

- **Steps 1–5 of the cleanup pass complete** (preprocessing, HADDOCK,
  RFDiffusion, Rosetta filtering, container consolidation)
- **HADDOCK test validated** end-to-end (7/7 processes, 13m 27s)
- **RFDiffusion tests validated** for both single-denovo and multi-denovo
  contigs after R9, R10, R7 fixes
- **Rosetta test validated** on real RFDiffusion split PDBs from previous
  block's output
- **AF2-monomer / FASTA-input branch retired** and archived
- **Container paths consolidated** into `nextflow.config` as single
  source of truth
- **ProteinMPNN block remains** (step 6) — substantial deferred-work
  list ready, will be done in a fresh conversation
- **MSA, Boltz2, aggregate blocks remain** (steps 7–9) — not yet audited
- **Final whole-pipeline pass on `main.nf`** remains

## Next steps

### Tomorrow / next session — ProteinMPNN block

Start in a fresh conversation with the deferred handover doc
(`deferred_to_proteinmpnn.md`) pasted in. The MPNN block is where
several deferred items from this session land:

1. **`receptor_start_pdb` plumbing into MPNN.** TODO marker in
   `main.nf:352`. `BUILD_CONTIGS` and `EXTRACT_SEQUENCES` both write
   the real value to JSON now, but `proteinmpnn.nf:41,131` still uses
   `params.receptor_start_pdb ?: 1` and silently defaults to 1. Fix:
   parse from JSON in the seqs_ch map, expose as a value channel,
   pass through `MPNN_FIXED_POSITIONS` and `PROTEINMPNN` as a process
   input.

2. **`denovo_length` rename and Branch B correctness in
   `pipeline_correct_sequences.py`.** TODO marker in place. Two issues:
   (a) Rename `denovo_length` → `design_region_length_observed`. Add a
   sibling `design_region_spec` column carrying the contig spec range.
   Drop the misleading gap-arithmetic `design_residues` field from
   `rfdiffusion_metrics.json`. Apply all of this consistently across
   CSV/JSON/plots. (b) The calculation assumes the contig string is
   PDB-resolved, which is true for Branch A (HADDOCK + BUILD_CONTIGS)
   but not Branch B (pre-docked complex PDB). Either call
   `resolve_contigs` inside the script when unresolved input is
   detected, or always run resolution in `main.nf` Branch B.

3. **R7 cross-pollination.** The same proportional-allocation pattern
   exists in `pipeline_correct_sequences.py:280-320` (`align_to_native_by_anchors`)
   and has the same bounds-violation bug. The R7 fix in
   `rfdiffusion_filter.py:build_receptor_resnum_map` is the reference
   implementation. Copy the structure (fixed-length segments first,
   variable-length segments share leftover with proportional split and
   per-segment clamping).

4. **R12 design_region_rmsd.** Needs Kabsch superposition on the fixed
   segments before RMSD-ing the design region. Currently produces
   meaningless 100+ Å values. Not a filter gate, so not blocking, but
   the reported field is currently noise. ~30 lines of careful numpy
   in `calc_design_region_rmsd`. Or: report fixed-segment-only RMSD
   (how much the scaffold moved during diffusion) and skip the design
   region for variable-length tasks where per-residue correspondence
   isn't even well-defined.

5. **H15 `min_population` consolidation follow-up.** If MPNN's Nextflow
   param-passing pattern is clean enough to use as a model, back-apply
   to the four hardcoded `=4` constants in the HADDOCK block. Otherwise
   leave the cross-referencing comments in place.

### After ProteinMPNN — remaining cleanup blocks

6. **Step 7: MSA block.** Not yet looked at. Probably small —
   ColabFold MSA generation with caching.

7. **Step 8: Boltz2 block.** The biggest remaining block. Multiple
   processes (`BOLTZ2_PREPARE`, `BOLTZ2_PREDICT`, `BOLTZ2_VERIFY_BINDING`,
   `BOLTZ2_FILTER_AND_RANK`, `BOLTZ2_PLOTS`). Check whether
   `boltz2_prepare.py:extract_design_residues_from_contigs` has the
   same H6-class chain-substring bug pattern surfaced repeatedly
   elsewhere.

8. **Step 9: Aggregate block.** `AGGREGATE_RESULTS`. Probably small
   if it's just collecting per-design rows into a final summary.

9. **Final pass on `main.nf`** as a whole-pipeline review once every
   block is individually clean. Check that channel wiring between
   blocks is consistent and that no orphan emit channels exist. Check
   that the Branch A vs Branch B logic is symmetrical, and that
   `main.nf`'s consumers have picked up all the per-block fixes.

### Future work (not urgent)

10. **RFDiffusion parallelisation across SLURM submissions.** Currently
    in Claude's persistent memory. RFDiffusion runs as a single GPU job
    that produces all `num_designs` designs in one go. The MPNN block
    uses a `max_af2_parallel`-style fan-out across many concurrent GPU
    jobs. The same pattern should be applied to RFDiffusion. Best done
    after the MPNN audit when the fan-out pattern is fresh.

### Filed but not fixed

11. **Two-singularity-exec waste in `ROSETTA_SC`.** The process does
    one `singularity exec` to detect the binary suffix, then another
    to actually run InterfaceAnalyzer. Could be combined into one
    shell session. Not a bug, just wasteful — left alone in
    cleanup-only mode.

12. **`-no_optH false` on polyvaline backbones in `ROSETTA_SC`.** Runs
    hydrogen optimisation on structures that have no real sidechains
    to optimise. Wasted compute, not a bug.

### Carried forward from previous sessions

13. **The `split_sequence` chain-order bug** was fixed in notes 4's
    session and is no longer an open issue. But the *class* of bug has
    now been independently re-discovered three times in this pipeline
    (notes 4, plus H10 and R10 from this session). That class is
    length-based chain disambiguation failing on similar-length chains.
    Worth keeping it on the watch list for any future code that needs
    to distinguish receptor from effector by anything other than chain
    ID or sequence identity.
