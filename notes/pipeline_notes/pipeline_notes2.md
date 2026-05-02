# Receptor-Resurfacing Pipeline Cleanup Pass — ProteinMPNN Block Session

## Context

This session picked up the deferred-work list left at the end of
`pipeline_notes1.md`. The previous session had cleared steps 1–5 of the
cleanup pass (preprocessing, HADDOCK, RFDiffusion, Rosetta filtering,
container consolidation) and left a five-item handover specifically tied
to the ProteinMPNN block:

1. `receptor_start_pdb` plumbing through MPNN.
2. `denovo_length` rename and Branch B contig-resolution correctness.
3. R7 cross-pollination (the proportional-allocation bug previously
   fixed in `rfdiffusion_filter.build_receptor_resnum_map`, still
   present in `pipeline_correct_sequences.align_to_native_by_anchors`).
4. R12 design-region RMSD (currently producing 100+ Å nonsense values
   because there's no Kabsch superposition before the RMSD calculation).
5. H15 `min_population` consolidation across the four sites in the
   HADDOCK block.

The plan was to do all five, run the four available test suites (HADDOCK,
RFDiffusion, Rosetta filtering, ProteinMPNN), and stop. Same ground rule
as last session — minimal targeted changes, no scope creep — but with
one explicit exception: the user indicated upfront that "if it's broken,
fix it now" was preferred over "leave a TODO marker for next time."

What actually happened: items 1, 2, 3, 5 went roughly as planned. Item 4
expanded into a much larger investigation when the test results made it
clear that the original RMSD function had an atom-correspondence bug
inherited from before this session, and that "fix the RMSD" was the wrong
goal in the first place — what the user actually wanted to measure was
"how much did RFDiffusion change the structure," for which RMSD is a poor
choice. We replaced the broken RMSD function with three new
geometry-based metrics, then audited and re-fixed each of them as new
edge cases surfaced. Then a follow-on observation about the
`test_proteinmpnn.nf` test file led to a small but important architectural
fix: the test had its own inline duplicate of `EXTRACT_SEQUENCES`
diverging from production, which we removed and replaced with the
production module's actual process.

## What we built

### 1. `receptor_start_pdb` plumbing through MPNN (Item 1)

The previous session had fixed `EXTRACT_SEQUENCES` and `BUILD_CONTIGS` to
write `receptor_start_pdb` to their JSON outputs, but the MPNN-side
consumers (`MPNN_FIXED_POSITIONS` and `SEQUENCE_CORRECTION` in
`proteinmpnn.nf`) were still using `${params.receptor_start_pdb ?: 1}` as
a fallback that silently defaulted to 1 even when the producer had
computed the real value. This silently mis-aligned MPNN fixed positions
for any input PDB whose receptor chain didn't start at residue 1.

The fix is straightforward plumbing: parse `receptor_start_pdb` out of
the JSON in `main.nf`'s `seqs_ch` map alongside the existing receptor
and effector sequences, expose it as a value channel
(`receptor_start_pdb_ch`), and pass it as a new `val receptor_start_pdb`
input to both `MPNN_FIXED_POSITIONS` and `SEQUENCE_CORRECTION`. Each
process's script block then references the input variable instead of the
params fallback.

Touched: `main.nf`, `proteinmpnn.nf`, `test_proteinmpnn.nf`. The test
file's inline `EXTRACT_SEQUENCES` (which we ended up deleting later in
the session — see Item 6 below) was also updated to write
`receptor_start_pdb` to match production at the time.

While reading `main.nf` we also found dead code from the old
user-provided-sequences entry point (`if (params.receptor_seq && params.effector_seq)`
branches in the seqs_ch resolution and Branch B's `EXTRACT_SEQUENCES`
conditional). This branch was retired with the FASTA-input archive in
the previous session but the surrounding `if/else` shells remained.
Removed both, plus a `Channel.empty()` dead else-branch in Branch B that
was unreachable once the params were always null.

Two `params.receptor_seq` / `params.effector_seq` references were left
in place: their `null` defaults at the top of `main.nf` and the
`params.receptor_seq ?: ""` arguments passed to `EXTRACT_HOTSPOTS` as
the H10 chain-disambiguation fallback hint. Both effectively dead but
removing them would require touching `extract_hotspots.py`'s CLI, which
was out of scope.

### 2(a). `denovo_length` → `design_region_length_observed` rename

The CSV column `denovo_length` in `sequence_metadata.csv` was
ambiguous — it could mean the contig spec range ("what was asked for"),
the per-design observed length ("what RFDiffusion actually produced"),
or gap arithmetic on the input PDB ("what the original gap was"). The
three concepts had been silently conflated, which led to misleading
plots and confused analysis.

Rename: `denovo_length` → `design_region_length_observed`. Added a new
sibling column `design_region_spec` carrying the per-segment "min-max"
range as written in the contig (e.g. `"20-40"` for a single variable
region, `"5-5|10-30"` for a two-region contig). Both columns are
pipe-joined for multi-region contigs so each region's observed length
and spec sit at the same index in their respective columns.

The rename propagated through `pipeline_correct_sequences.py`: the
function return value, the CSV row dict, the FASTA header in
`mpnn_corrected.fasta`, and the per-design progress print line. Also
threaded through `sequence_metadata.csv` → `qc_metadata.csv` →
`scored_metadata.csv` automatically because the downstream consumers
read the dict from the upstream CSV unchanged.

Touched: `pipeline_correct_sequences.py` only.

Also dropped the misleading top-level `design_residues` field from
`rfdiffusion_metrics.json` (written by `rfdiffusion_filter.py` lines 489
and 674). This field was gap arithmetic from `parse_design_region` on the
input PDB, which only matched reality for fixed-length contigs. The
authoritative per-design source `per_design_design_residues` (computed
via the R7-fixed `build_receptor_resnum_map`) is still present, so
consumers should use that. `rfdiffusion_plots.py` still references the
top-level field at line 834 with a `.get(..., [])` default, which means
plots silently use an empty fallback set rather than erroring — verified
this didn't break the plot test, the fallback path is not actually used
because `_get_per_design_residues` prefers the per-design field.

Touched: `rfdiffusion_filter.py` (two write-site edits). `rfdiffusion_plots.py`
deliberately not updated yet — we did rewrite it later in the session,
but for a different reason (the new COM displacement plot from Item 4).

### 2(b). Branch B contig resolution via new `RESOLVE_CONTIGS` process

Branch A (HADDOCK input) resolves the user contig string to PDB
coordinates inside `BUILD_CONTIGS`, which also handles HADDOCK-driven
residue remapping. Branch B (pre-docked complex PDB) was just passing
`params.contigs` through unresolved as a value channel, meaning every
downstream consumer had to either re-resolve internally or assume
resolved input. `pipeline_correct_sequences.py` assumed resolved input
and silently mis-counted residues for Branch B inputs.

The fix is structural: the invariant should be **`contigs_ch` is always
PDB-resolved by the time it leaves the branch block**, regardless of
which branch produced it. New process `RESOLVE_CONTIGS` in
`preprocessing.nf` wraps the existing `bin/rfdiffusion_contigs.py` CLI
(which already exposed `contig_utils.resolve_contigs` as a standalone
script). Branch B in `main.nf` now calls it on the input PDB and
produces a resolved `contigs_ch` value channel mirroring Branch A's
shape (`.map { it.text.trim() }.first()`).

Branch A's `BUILD_CONTIGS` is left alone — it does more than just
resolution (HADDOCK coordinate remapping, sequences JSON), so unifying
the two would have been overengineering.

Touched: `preprocessing.nf` (new process), `main.nf` (include + Branch B
wiring).

### 3. R7 cross-pollination in `align_to_native_by_anchors`

The previous session's R7 fix in `rfdiffusion_filter.build_receptor_resnum_map`
restructured the de novo segment allocation as: fixed-length segments
first (exact assignment), then variable-length segments share leftover
budget with proportional split and per-segment clamping. The same bug
existed in `pipeline_correct_sequences.align_to_native_by_anchors`
(lines 273–320) and had the same shape — purely proportional allocation
that violated fixed-length specs and didn't properly clamp the final
segment.

Ported the R7 structure verbatim, adapted to the local dict-keyed
segment representation (`s['min_len']`, `s['max_len']` instead of
tuples). Replaced the inline allocation-during-walk loop with a
precomputed `seg_index_to_alloc` dict; the walk loop now just looks up
the allocation per de novo segment by index. Three classes of bug
fixed:

- `min_len == max_len` segments (e.g. `"6-6"`) are now respected exactly
  instead of being proportionally split.
- Final variable segment is now bounds-clamped via the variable-length
  share rather than swallowing all remainder unconditionally.
- Clamping order corrected so allocations can never fall below `min_len`.

Also added a length-mismatch warning when the contig accounting doesn't
match the MPNN receptor length (`total_fixed + sum(allocs) != len(mpnn_receptor)`).
This caught the edge case where a fixed-length de novo segment's
declared length disagrees with the actual design — previously silent.

Touched: `pipeline_correct_sequences.py` only.

### 4. R12 design-region RMSD — replaced wholesale

Started as "fix the RMSD function." Ended as "replace it with metrics
that actually answer the question."

The original `calc_design_region_rmsd` did a positional-index comparison
(`input_coords[i]` vs `design_coords[i]`) without any superposition,
producing 100+ Å values dominated by global rigid-body motion. First
attempt was to add Kabsch superposition on the fixed scaffold residues
and compute RMSDs on the aligned coordinates — which fixed the
superposition issue but still had a hidden bug: the positional-index
correspondence between input and design atom lists is only valid when
the design region has the same length as the input gap, which it
doesn't for variable-length contigs (the common case).

The deeper problem was that **per-residue RMSD is the wrong metric** for
"how much did RFDiffusion change the structure." For variable-length
designs there's no per-residue correspondence between the new region
and what was there before; for fixed-length designs there is but you'd
have to commit to a strategy that handles both. After a back-and-forth
discussion about what we actually wanted to measure, we replaced the
single-RMSD function with three new metrics:

- **`scaffold_rmsd`** — single float, RMSD over the fixed-segment Cα
  atoms after Kabsch superposition on those same atoms. Sanity check
  on RFDiffusion's faithfulness to the contig — should be small
  (typically <1 Å) for a well-behaved design.
- **`endpoint_distance_input`** / **`endpoint_distance_design`** —
  per-design-region lists. End-to-end Cα distance of the input gap
  residues (input) and the design region's residues (design). Compares
  whether RFDiffusion built something spanning the same physical reach
  as the original gap.
- **`design_region_com_displacement`** — per-design-region list.
  Distance from the centroid of the design region's Cα atoms to the
  centroid of the input gap's Cα atoms, in the scaffold-aligned frame.
  Length-independent, so it works for both fixed and variable-length
  contigs.

The function is now `calc_scaffold_and_region_metrics`. Critically, it
uses **resnum-keyed atom correspondence** rather than positional
indexing — each design atom is matched to its input atom via
`resnum_map`'s recorded input residue number, not by list position.
This was the real underlying bug: positional indexing fails as soon as
the design region has a different atom count than the input gap,
shifting every fixed atom that follows.

Three subtleties fell out of testing:

- **Trailing-anchor-only design regions** (a `6-6` C-terminal extension
  beyond the last fixed segment, for example) need special handling.
  Endpoint distance is genuinely undefined (only one anchor), but COM
  displacement is well-defined as long as the input PDB contains the
  residues immediately following the preceding fixed segment. Added
  leading-anchor-only / trailing-anchor-only / both / neither branching
  in the function. For one-sided cases the comparison uses
  `n_region` input residues immediately adjacent to the available
  anchor — same atom count on both sides, symmetrical comparison.
- **The first test fixture didn't reach far enough into the receptor
  sequence** to actually exercise the leading-only branch with real
  data. The contig asked for 6 residues past the input PDB's
  C-terminus, the residues didn't exist, and the function correctly
  returned NaN. Confused this for a function bug initially. Resolved
  by switching the test fixture to an AF3-derived complex of the full
  Pikp-1 HMA domain (78 residues, 1–78) with a contig
  (`A1-32/10-20/A46-72/6-6 C`) that places the `6-6` region inside the
  natural domain so all six residues exist for comparison.
- **The original fixture had been built from the wrong PDB.** The user
  noticed mid-debugging that they'd built the input from 6G10 instead
  of 7QZD and had been missing the C-terminal residues entirely. This
  was a fixture issue, not a code issue, but it sent us through a
  longer-than-expected loop of "is the function wrong or is the input
  wrong?"

Also added a new diagnostic plot to `rfdiffusion_plots.py`:
`rfdiff_com_displacement.png`. Grouped bar chart, one slot per design,
sorted left-to-right by maximum per-region displacement, viridis
colours per region with a legend for multi-region contigs. Defensive
about None entries and missing data. Verified that the four existing
plots in `rfdiffusion_plots.py` still run cleanly after the
`design_residues` JSON-field drop from Item 2(a) — none of them
actually used the field, just the `_get_per_design_residues` helper
which prefers the still-present per-design source.

Touched: `rfdiffusion_filter.py` (function rewrite, call site update,
metrics dict update, progress print update, new
`_parse_receptor_segments` precompute outside the design loop),
`rfdiffusion_plots.py` (new plot function, register in `ALL_PLOT_FILES`,
add to `main()`).

### 5. H15 `min_population` consolidation

Four hardcoded `min_population = 4` / `MIN_CLUSTER_SIZE = 4` constants
across the HADDOCK block (`haddock.nf` heredoc, `collect_haddock3_dock.py`,
`haddock3_plots.py`, `haddock_utils.py` docstring). Previous session
left cross-referencing comments at all four sites with the consolidation
deferred to the MPNN block "where the param-passing pattern is clearer."

Single source of truth now lives in `nextflow.config` as
`params.haddock_min_cluster_size = 4`. The value is substituted into
the HADDOCK config heredoc in `haddock.nf` line 127
(`min_population = ${params.haddock_min_cluster_size}`), and passed via
`--min-cluster-size` CLI flag to both Python script invocations
(`collect_haddock3_dock.py` line 147, `haddock3_plots.py` line 184).

The Python scripts use a small idiom to keep the diff minimal: the
module-level `MIN_CLUSTER_SIZE` constant is preserved, the new
`--min-cluster-size` CLI argument defaults to it, and `main()`
re-binds the constant from the CLI arg via `global MIN_CLUSTER_SIZE`
before any other code runs. This means all the existing internal
references (function bodies, `parse_clustfcc_tsv` calls, error messages)
keep working unchanged — no function signature plumbing required.

`haddock_utils.py` (docstring reference only) left as-is — out of scope
for the consolidation since it's not a real call site.

Touched: `nextflow.config`, `haddock.nf`, `collect_haddock3_dock.py`,
`haddock3_plots.py`.

### 6. `test_proteinmpnn.nf` test/production divergence — fixed

Discovered late in the session, after the user pushed back on a
chain-renaming workaround for what looked like a chain-layout problem.
The actual problem was that `test_proteinmpnn.nf` had its own inline
copy of `EXTRACT_SEQUENCES` (a simplified Python heredoc inside the
test file's process block) that read chain letters from `params.input_pdb`
in a way that conflicted with how the MPNN block's downstream processes
treated the same chain-letter params for the split design PDBs. The
test was effectively reimplementing the input-PDB contig and chain
plumbing in a buggier way than `main.nf` Branch B did it.

Surveyed all five test files: only `test_proteinmpnn.nf` had this
problem. The other four (`test_haddock.nf`, `test_rfdiffusion.nf`,
`test_rosetta_filtering.nf`, `test_boltz2.nf`) all already imported
their helper processes from production modules without inline
duplicates.

The fix was to delete the inline `EXTRACT_SEQUENCES` and rewrite the
test workflow to mirror `main.nf` Branch B end-to-end up to the MPNN
block: load the input PDB, call the production `RESOLVE_CONTIGS` to
resolve the contigs, call the production `EXTRACT_SEQUENCES` to extract
sequences and `receptor_start_pdb`, parse the JSON into value channels,
feed everything into the MPNN processes. The test now exercises the
exact same code paths as production, so chain-letter handling, contig
resolution, residue numbering, and `receptor_start_pdb` plumbing all
behave identically here and in `main.nf`.

The chain-layout confusion that triggered this work resolved itself
once production `EXTRACT_SEQUENCES` was in use — `params.receptor_chain`
in production is the chain letter the user specifies for *the input
PDB* (and for the contig string, which uses the same convention),
and the MPNN block handles split-design-PDB chain layout internally
via length-based identification in `pipeline_correct_sequences.py`.
The two roles of `params.receptor_chain` were never actually conflated
in production; only the test file conflated them via its inline copy.

Touched: `test_proteinmpnn.nf` (substantial rewrite, smaller and
cleaner than the original).

## Test results

All four available tests pass end-to-end on the new AF3-derived
PikP1/AVR-PikF fixture with the literature-derived contig
`A1-32/10-20/A46-72/6-6 C`:

- **HADDOCK** ✅ — H15 plumbing verified, cluster threshold of 4
  honoured end-to-end, top model selected from cluster 1 (4 members),
  capri scores match expected values.
- **RFDiffusion** ✅ — Item 4 fix verified with sub-Å scaffold RMSDs
  (0.11–0.14 Å) across all four designs, both endpoint distances
  reported correctly per region, COM displacement reported correctly
  for both regions including the trailing `6-6` via the leading-only
  branch (3.4–4.3 Å for region 2). Item 2(a) verified with
  `design_residues` top-level field absent from JSON. The new
  `rfdiff_com_displacement.png` plot rendered successfully.
- **Rosetta filtering** ✅ — untouched by this session, regression-tested
  as a side effect, 4/4 passing on the new fixture.
- **ProteinMPNN** ✅ — exercises Items 1, 2(a), 2(b), 3, and 6 all at
  once. `processed_contigs.txt` shows correct Branch B resolution
  (`A1-32/10-20/A46-72/6-6 C32-113`). `sequence_metadata.csv` headers
  show both `design_region_length_observed` and `design_region_spec`
  columns populated correctly. Per-design observed lengths (14, 18,
  19, 20 for region 1; 6 for region 2) match RFDiffusion's output and
  the contig spec. Zero `denovo` references anywhere in the metadata
  pipeline. No length-mismatch warnings, meaning the Item 3 allocation
  logic produces self-consistent accounting on the multi-region contig.
  16 sequences passed QC, 16 sequences cluster into 7 groups at 30%
  identity with sensible monotonic increase to 16 unique groups at
  100%. End-to-end runtime ~4 minutes.

## Test fixture updates

The test fixture for the second half of the session was rebuilt to
match the actual biology of the design problem. Original fixture used
a 6G10-derived complex that was missing the C-terminal residues of the
Pikp-1 HMA domain. New fixture is `pikp1_avrpikf_complex.pdb`,
constructed in ChimeraX from:

- **Receptor**: AlphaFold3 prediction of the full Pikp-1 HMA domain
  (residues 186–263 of wild-type Pikp-1, renumbered to 1–78 in the
  PDB), aligned by sequence-based superposition onto 7QZD's chain B.
- **Effector**: AVR-PikF chain C from 7QZD (residues 32–113), kept in
  its crystal coordinates from the alignment.

The complex is saved with chain A (receptor, 1–78) and chain C
(effector, 32–113), no other chains present. Contig
`A1-32/10-20/A46-72/6-6 C` corresponds to (in wild-type numbering):
fix residues 186–217, design 10–20 residues, fix residues 231–257,
design exactly 6 C-terminal residues (residues 258–263). All design
regions fall inside the natural domain so every metric has meaningful
input residues to compare against.

## File deliverables from this session

- `main.nf` — Items 1, 2(b) (RESOLVE_CONTIGS wiring + dead-branch removal)
- `proteinmpnn.nf` — Item 1
- `test_proteinmpnn.nf` — Items 1, 6 (rewritten to mirror Branch B)
- `preprocessing.nf` — Item 2(b) (new RESOLVE_CONTIGS process)
- `pipeline_correct_sequences.py` — Items 2(a) part 1, 3
- `rfdiffusion_filter.py` — Item 2(a) part 2, Item 4 (function replaced
  with `calc_scaffold_and_region_metrics`, three new metrics, resnum-keyed
  correspondence)
- `rfdiffusion_plots.py` — Item 4 follow-on (new
  `plot_com_displacement` function, registered in `ALL_PLOT_FILES`)
- `nextflow.config` — Item 5 (new `haddock_min_cluster_size` param)
- `haddock.nf` — Item 5 (heredoc substitution + CLI flag plumbing)
- `collect_haddock3_dock.py` — Item 5 (CLI flag + global re-bind idiom)
- `haddock3_plots.py` — Item 5 (CLI flag + global re-bind idiom)

## Status summary

- **All five deferred items from `pipeline_notes1.md` complete** and
  verified by the per-module test suites.
- **Item 4 audit and fix complete** — the original RMSD bug was deeper
  than expected (atom-correspondence, not just superposition), and the
  metric was the wrong measurement in the first place. Replaced with
  three geometry-based metrics that are length-independent and answer
  the user's actual question.
- **Test/production divergence in `test_proteinmpnn.nf` fixed** —
  inline `EXTRACT_SEQUENCES` duplicate deleted, test now mirrors
  `main.nf` Branch B end-to-end via production module imports.
- **All four available tests pass** on the new AF3-derived fixture.
- **MSA, Boltz2, and aggregate blocks remain** (steps 7–9) — never
  audited, no work done this session.
- **Final whole-pipeline pass on `main.nf`** remains.

## Next steps

### Tomorrow / next session — MSA, Boltz2, aggregate blocks

These three blocks have never been audited. Pick them up in order:

1. **MSA block.** Probably small — ColabFold MSA generation with
   caching. Should be fast to read and audit.

2. **Boltz2 block.** The biggest remaining block by far. Multiple
   processes (`BOLTZ2_PREPARE`, `BOLTZ2_PREDICT`, `BOLTZ2_VERIFY_BINDING`,
   `BOLTZ2_FILTER_AND_RANK`, `BOLTZ2_PLOTS`). The handover doc from
   `pipeline_notes1.md` flagged a possible H6-class chain-substring
   bug in `boltz2_prepare.py:extract_design_residues_from_contigs`
   that nobody has confirmed or fixed.

3. **Aggregate block.** `AGGREGATE_RESULTS` only. Probably small.

### Deferred follow-ups from this session

1. **`rfdiffusion_plots.py` cleanup pass.** The four legacy plots
   (`plot_contact_map`, `plot_design_clustering`,
   `plot_specificity_coverage`, `plot_design_lengths`) still
   read the now-removed top-level `design_residues` field with a
   `.get(..., [])` default. Verified harmless because the
   `_get_per_design_residues` helper prefers the per-design source,
   but worth removing the dead fallback path for clarity.

2. **`haddock_utils.py` `min_population` docstring reference.** Cosmetic
   only; the docstring still mentions the constant by name but it's
   not a real call site. Update the docstring to point to
   `params.haddock_min_cluster_size` in `nextflow.config` as the
   single source of truth.

3. **HADDOCK test threshold-fragility.** The test fixture produces
   exactly 4 cluster members at the `min_population = 4` threshold.
   Any stochastic variation that produces one fewer member will fail
   the test with "No cluster has >= 4 models." Either lower the
   test-specific `params.haddock_min_cluster_size` to 3 or pick a
   fixture that produces more obvious clustering.

4. **`params.receptor_seq` / `params.effector_seq` declarations in
   `main.nf`.** Lines 48–49 still declare these as `null` defaults,
   and lines 288–289 still pass them to `EXTRACT_HOTSPOTS` as the H10
   chain-disambiguation fallback hint. Both effectively dead. Removing
   them would require touching `extract_hotspots.py`'s CLI to drop the
   two flags. Out of scope for this session, worth a small follow-up.

5. **Final whole-pipeline pass on `main.nf`.** After Boltz2 and
   aggregate are done, one read-through of `main.nf` end-to-end to
   catch anything that drifted during all the per-block work.
