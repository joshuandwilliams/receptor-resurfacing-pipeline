# Pipeline Notes 11 — Session 27 Apr 2026

Plot iteration: RFDiffusion, Rosetta filtering, ProteinMPNN. Step 3 of the 1–15 plan is partially closed.

## Headline

Three module-plot scripts iterated to v2 and integrated to production:

- `rfdiffusion_plots.py` — five plots tuned across many micro-fixes (calibrated dendrogram placement, hatch-line visibility, x-tick collision fix, per-region clustering calibration, consistent physical bar widths in design-lengths). Production replacement verified end-to-end on synthetic data.
- `rosetta_filter_plots.py` — replaced `sc_vs_dg` scatter with a 1D Sc histogram. dG_separated dropped because at this pre-MPNN stage it's dominated by valine-clash artefacts and is not physically meaningful (documented in notes 1). Best-model star marker dropped — redundant in 1D and was the gnarliest code in the file. Filtered region now flush against the left axis spine.
- `mpnn_plots.py` — score panels flipped (design-region top, primary; global bottom, sanity check) and re-ranked by design-region median; AA composition gets a native-residue overlay (red diamonds); new physicochem plot (mean Kyte–Doolittle hydrophobicity + net charge per residue, with native horizontal reference); sequence diversity left as-is at user request.

One real upstream bug found and fixed in `pipeline_correct_sequences.py` along the way — dropped `native_residues` for trailing denovo regions. Discovered because the MPNN AA-composition native overlay looked too sparse on real data (only one region's worth of residues).

Step 3 is now closed for the upstream three modules. Tasks 43 (negsteer plots) and 44 (orthogonal-metrics plots) deferred to a fresh conversation — distinct enough in scope that mixing them with this session's work would have been counter-productive.

## Workflow that emerged this session

For each module, I built two test-harness files in `tests/<module>/`:

- `test_<module>_plots.py` — standalone iteration script that reads cached test outputs (the metrics JSON/CSV from a previous test run) and re-emits the plots without re-running the module. No GPU, no Nextflow.
- `run_test_<module>_plots_slurm.sh` — SLURM wrapper that invokes the test script inside the right container with the right paths.

This let me iterate plot-only changes on the order of seconds-per-revision instead of minutes (or hours, for MPNN). Once the user signed off on the test outputs, the corresponding production `bin/<module>_plots.py` was patched in-place with the same logic and re-verified end-to-end against the same synthetic data the test script used.

Key environment knowledge that re-surfaced (carried over from earlier sessions):

- `nextflow.config` overrides `params.outdir`, so the real test results dir is `tests/<module>/receptor_resurfacing_results/...`, NOT `tests/<module>/results/...`.
- LRR_Pipeline container at `/hpc-home/jowillia/singularity/HADDOCK_RFDiffusion_ProteinMPNN_MMseqs2/LRR_Pipeline.img` — must use `python` (not `python3`), since `python3` is a system stub without numpy.
- Rosetta container at `/hpc-home/jowillia/singularity/Rosetta/Rosetta.img` for the rosetta plots test.
- SLURM wrappers should write log files named `test_plots_%j.out`/`.err`, partition `jic-medium`, mail to `jowillia@nbi.ac.uk`.
- The production `publishDir "${params.outdir}/plots", mode: 'copy'` glob `<module>_*.png` already matches any new filename starting with the module prefix — no `.nf` changes needed when the only thing that changes is the plot's filename suffix.

## RFDiffusion plot fixes

Five plots were already in place. None needed adding; all five were touched.

### Plot 1 — contact map (`rfdiff_contact_map.png`)

The "calibrated G" layout from the previous session held up. Fixes this session were narrower:

- **X-tick endpoint collision (the `82/83` and `131/132` overlap).** Original code blindly appended the final residue index to `xtick_positions` if not already present. With a typical step of 3 across ~83 positions, the natural ticks ended at index 81 ("82") and the appended endpoint at index 82 ("83") — two labels one position apart, visibly overlapping. Fix: only append the endpoint if it's at least one full tick step from the previous tick; otherwise replace the previous-to-last tick with the endpoint. Same fix applied to the effector-frequency bar's x-ticks (it had the same pattern).
- **Beyond-protein-end hatching invisible.** The grey "beyond protein end" rectangles in the legend showed diagonal hatching, but the actual rectangles in the plot rendered as solid grey with no diagonals. Root cause: matplotlib renders hatching as ornament on the rectangle's edge, and `edgecolor="none"` suppresses it. Fix: `edgecolor="#888888", linewidth=0` — invisible border but visible hatch lines.

### Plot 2 — design clustering (`rfdiff_design_clustering.png`)

Per-region pairwise RMSD heatmaps with dendrogram insets, one panel per design region. This took the most iteration. Final state:

- **Per-region dendrogram calibration.** Heatmap left edge stays fixed across panels (anchored on the widest panel's longest-label width) so the heatmap columns line up vertically across regions for cross-region comparison. But the dendrogram inset is positioned per-panel — its right edge sits `LABEL_DENDRO_GAP_IN` (0.20") to the left of *that panel's own* longest label. So a panel with shorter labels gets its dendrogram drifted rightward, keeping the visual label-to-dendro gap consistent everywhere instead of leaving big variable white-space gaps when label widths differ.
- **Square heatmap cells** (`aspect="equal"`) — clearer per-cell value reading than aspect="auto".
- **Dedicated colorbar axis** to the right of each heatmap, with a `FormatStrFormatter("%.1f")` so tick labels read as `0.0, 0.5, 1.0` instead of `0, 1, 2`.
- **Pass/Filtered legend repositioned** above the heatmap, right-aligned (matches the contact map's legend style — was previously inside the heatmap upper-right corner, looked cramped).
- **Vertical spacing tightened** after the user pointed out that there was too much whitespace between panels. Final values: `panel_gap=0.2"`, `xtick_room=0.4"`, `title_room=0.25"`, `legend_room=0.25"`. Top/bottom margins 0.3" each (overridden locally; module-level constants are 0.7" tuned for the contact map).

The label-width measurement helper `_measure_label_width_inches` is the calibration primitive the layout depends on — it renders the longest label to a throwaway figure at known DPI and reads the bounding box width. This is what makes the per-panel calibration possible without hard-coded magic numbers.

### Plot 3 — specificity coverage (`rfdiff_specificity_coverage.png`)

Untouched. User signed off in an earlier session.

### Plot 4 — design lengths (`rfdiff_design_lengths.png`)

Most-iterated plot in the session. Discrete integer bars, one panel per design region. Three rounds of fixes:

1. **Per-panel sizing.** Original used `subplots(n_regions, 1)` with each panel filling the figure width regardless of how many unique values it had. A region with one unique length (e.g. region 2 in the test data, with all 8 designs landing on length=6) got a single bar stretched across the full panel width — looked like one giant bar. Fix: each panel's axis width scales with its number of unique length values, using `fig.add_axes` for absolute control.
2. **Black bar borders** (`edgecolor="black", linewidth=0.5`) to make each bar visually distinct.
3. **Consistent physical bar width across panels.** First attempt computed `target_data_range = bar_width_data * ax_w_in / target_bar_w_in` from a pre-set `ax_w_in`, but for a 6-bar region this gave a data range slightly smaller than the bar span, so the first and last bars rendered half-clipped by xlim. Reversed the dependency: compute `target_data_range = (n_bars - 1) + bar_width + 2 * edge_pad` first, derive `ax_w_in` from that. Bars now render flush. **And then** the user pointed out that with the `min_panel_w_in=3.0` floor (needed so the xlabel "Design-region length (residues)" doesn't crop on narrow panels), single-bar panels had axes wider than `natural_ax_w_in`, which inflated the bar's physical width back to ~1.5". Final fix: when `ax_w_in > natural`, scale `target_data_range` up proportionally, so the bar stays at 0.6" physical regardless of how many bars share the panel.
4. **X-ticks.** User asked for ticks between bars but not extending past the data. Implementation: `set_xticks(range(min_unique, max_unique + 1))` — every integer between the lowest and highest unique value, no padding. For dense data (no gaps in unique values) this matches the unique values exactly; for sparse data (e.g. unique = [13, 15, 17]) it adds ticks at 14 and 16 too, but never below 13 or above 17. Crucially, removed the `MaxNLocator(integer=True)` call on the x-axis — it was over-riding the explicit `set_xticks` and producing too many ticks.

### Plot 5 — COM displacement (`rfdiff_com_displacement.png`)

Per-design centre-of-mass displacement of each design region from the input gap residues, scaffold-aligned. Two small fixes:

- **Title removed** — y-axis label now carries the meaning instead.
- **Y-axis label expanded** to two lines: `"Design region COM displacement\nfrom input gap, scaffold-aligned (Å)"`.

Bar borders (`edgecolor="black", linewidth=0.5`) added to match the design-lengths plot's idiom.

### What's in the file now

- New module-level constants `DENDRO_W_INCHES`, `LABEL_DENDRO_GAP_IN`, `RIGHT_MARGIN_IN`, `TOP_MARGIN_IN`, `BOTTOM_MARGIN_IN` for layout tunables.
- New helpers `_label_with_region_len`, `_per_region_clustering_data`, `_n_regions`, `_build_contact_map_data`, `_build_effector_freq`, `_build_design_labels`, `_draw_contact_heatmap_content`, `_set_contact_xticks`, `_draw_effector_freq_bar`, `_contact_legend_handles`, `_draw_clustering_heatmap_content`, `_measure_label_width_inches`. Most were promoted from inline code inside the original plot functions to enable reuse and per-panel calibration.
- `FormatStrFormatter` added to the matplotlib import block.
- All five `plot_*` function bodies replaced (specificity_coverage left unchanged in the diff but the surrounding section structure was reformatted).

Touched: `rfdiffusion_plots.py` (954 → 1248 lines).

## Rosetta filtering plot fixes

Single plot. Substantial change in scope.

### From `rosetta_sc_vs_dg.png` to `rosetta_sc_histogram.png`

The previous plot was a 2D scatter of Sc (x) vs dG_separated (y), with the filtered region shaded to the left of the Sc threshold and a star marker on the best-Sc design.

Notes 1 had documented that pre-MPNN PDBs are RFDiffusion polyvaline backbones with no real sidechains, so `dG_separated`, `packstat`, and `delta_unsatHbonds` are all dominated by valine-clash artefacts and are not physically meaningful at this stage. The pipeline correctly uses `Sc` (purely geometric) as the filter gate and ignores the energy terms here. So putting dG on the y-axis was actively misleading — implying the dG values were informative when they're not.

Three changes:

- **Replaced the scatter with a 1D Sc histogram.** Bin count by Freedman–Diaconis when n ≥ 12, otherwise a sensible fixed count.
- **Dropped the best-model star marker.** With only one axis, "best model" is just the rightmost histogram bar — redundant with the histogram itself. Removing it killed the gnarliest code in the original file: the post-render annotation bbox measurement and conditional x-axis-extension to fit the label text. ~30 lines of brittle code gone.
- **Filtered region flush against the left axis spine.** Original bug: `ax.axvspan` was being called before xlim was finalised, so the span's left edge sat at whatever xlim happened to be at that moment, and matplotlib's later autoscaling could shift the spine without shifting the span — leaving a visible white gap between the y-axis and the start of the red region. Fix: `ax.set_xlim(xlim_left, xlim_right)` *first* with explicitly computed bounds, then `ax.axvspan(xlim_left, threshold, ...)`. The span's left edge now coincides with the spine deterministically.

A vertical dashed line at the threshold marks the boundary between filtered and passing crisply (in case the pink fill is too faint on a particular display). Pass/fail counts go in the legend as a single line — `Sc < 0.62 (filtered, n=40 of 45)` — rather than separate "all designs" and "filtered" entries (the user's call: avoids an orphan blue patch in the legend that doesn't add information).

### Filename change

`rosetta_sc_vs_dg.png` → `rosetta_sc_histogram.png`. The Nextflow `ROSETTA_FILTER_PLOTS` process publishes with `path "rosetta_*.png"` so the new filename is picked up automatically — no `.nf` change needed.

`ALL_PLOT_FILES` list added to the file (single source of truth for the no-mpl fallback path), matching the pattern in `rfdiffusion_plots.py`. Removed the unused `mpatches` import (no star marker = no patch handles needed). Added `MaxNLocator` for integer y-tick locator on the count axis.

Touched: `rosetta_filter_plots.py` (196 → 209 lines).

## ProteinMPNN plot fixes

Three plots existed; this session ended with four. Big design discussion before any code was written, summarised here so the rationale survives the conversation.

### Design discussion

The user asked four questions before we started:

1. *Should the AA composition plot include the native sequence as a comparison?* Answer: yes, native (gap residues from the input PDB). Most informative for receptor resurfacing where the native interface was under selection. Comparing to non-designed regions of the same binder was the alternative considered and rejected — those regions were MPNN's *input*, so comparing to them is comparing to MPNN's prior, not to anything biologically meaningful. UniProt baseline considered as a third option but skipped to keep the plot scoped to the local biology question.

2. *Should the score-distribution panels be flipped, with design-region on top and ranking by design-region median?* Answer: yes. Design-region score is what we're actually evaluating (the newly designed bit). Global score is dominated by the long fixed regions and barely moves between designs — a sanity check, not a primary signal.

3. *How to improve the sequence-diversity plot?* User's view (correct, after I initially proposed a heatmap): keep the line chart. At production scale (hundreds-to-thousands of sequences) a pairwise distance heatmap is unreadable and expensive to compute. The line chart trades per-pair detail for scalability and stays.

4. *What other plots might be appropriate, given variable-length design regions?* Three candidates discussed: hydrophobicity/charge per design (chosen — sidesteps positional analysis by reducing to one number per axis); recovery-vs-score scatter (skipped — `mpnn_sequence_qc.py` already filters on `pct_identity`, so recovery is bounded to a known range and a scatter would be redundant); MPNN score per design region (skipped — covered by the score-distribution plot once it's flipped).

### Plot 1 — score distribution

Panels flipped (design-region top, global bottom). Ranked by median design-region score. X-axis label updated to `"Design (ranked by median design-region score)"`. Falls back to global-only / global-ranked when no design-region scores are present.

### Plot 2 — AA composition with native overlay

Heatmap + boxplot layout unchanged. New: red diamond markers per AA on the boxplot showing the native gap-residue composition. Legend in the upper-right reads `"Native"`.

`_aa_composition_pct` helper added — iterates a string (skipping `|` separators and any non-AA characters) and returns a 20-element AA frequency array in percent. Used for both the design rows and the native overlay.

### Plot 3 — sequence diversity

Unchanged. Production line chart kept.

### Plot 4 (NEW) — physicochem

Per-design boxplots of mean Kyte–Doolittle hydrophobicity (top panel) and net charge per residue (bottom panel), with native value drawn as a horizontal red dashed reference line on each panel. Charge is normalised by sequence length so designs with longer regions aren't artificially penalised. K, R = +1; D, E = −1; H skipped (≈ +0.1 at neutral pH, conventionally omitted).

`_mean_hydrophobicity` and `_net_charge` helpers added. `KD_HYDROPHOBICITY` and `CHARGE` constants added at module level.

### Bug found via the AA-composition overlay

After the user dropped the test scripts on real data, the native overlay looked far too sparse — 7 AAs at ~7.5% and 3 at ~15%, rest at zero. Total looked like ~10 distinct AAs across ~13 positions. The user's expected design region across both regions was substantially larger.

Investigation steps:

1. Re-checked the plot code — `_aa_composition_pct` correctly skips `|` separators and counts both regions in a multi-region native string. Plot code was fine.
2. User ran `awk 'NR==1 {for(i=1;i<=NF;i++) if($i=="native_residues") col=i; next} {print $col}' scored_metadata.csv | sort -u`. Result: one unique value, `SVALVGDLRDKIE|`. Region 1 is `SVALVGDLRDKIE` (13 residues). Region 2 is empty (just the trailing `|`).
3. Traced upstream to `pipeline_correct_sequences.py:382-411` — the `native_at_design_parts` construction. It only appends a non-empty native slice for denovo regions that have flanking fixed segments **on both sides**. For a contig like `B187-208/13-13/B223-262/6-6 C` (`fixed/denovo/fixed/denovo`), the trailing denovo has `prev_end` (the second fixed segment) but no `next_start` — so the function appends `""` for that region.
4. Result: `native_residues = "SVALVGDLRDKIE" + "|" + "" = "SVALVGDLRDKIE|"`. Region 2's native gap residues silently dropped from the metadata.

**Fix in `pipeline_correct_sequences.py`** (lines 382-411): handle three boundary cases that the original code missed.

- *Trailing denovo* (the user's bug, no fixed after): use `gap_end_0 = len(native_receptor)`.
- *Leading denovo* (no fixed before): use `gap_start_0 = 0`.
- *Both* (denovo as the entire receptor — pathological but representable): both clauses fire, slice is the whole native.

Tightened the bounds check from `gap_start_0 < gap_end_0 and gap_start_0 >= 0` to `0 <= gap_start_0 < gap_end_0 <= len(native_receptor)` to catch any out-of-range cases the new boundary logic might produce on weird inputs.

After re-running the ProteinMPNN test, the awk command showed `SVALVGDLRDKIE|<region2>` — both regions present. The plot's native overlay then matched the expected residue count.

Two reflection points worth recording:

- The plot was correct all along. The bug was upstream in metadata construction. The plot was faithfully showing what the metadata contained — it just contained too little. Initial diagnostic instinct was to keep editing the plot ("am I aggregating correctly?"); the right move was to inspect the data the plot was reading, which the user prompted by running the awk command. **For sparse-looking plots on real data, always check the input first.**
- The defensive collection of unique `native_residues` strings (gather to a set, pick the longest if multiple, warn) caught nothing here because the metadata had exactly one unique value — `SVALVGDLRDKIE|`. The defence isn't useless: it'd catch concatenated metadata from multiple runs. But the bug in this case was orthogonal to that defence.

### What's in the file now

- New constants: `KD_HYDROPHOBICITY`, `CHARGE`, `COLOUR_NATIVE`, `COLOUR_REGION`. `COLOUR_DESIGN` named in place of inline `"#4C72B0"`.
- `ALL_PLOT_FILES` includes `mpnn_physicochem.png`.
- `plot_score_distribution` body replaced — flipped, re-ranked, fallback handling.
- `plot_aa_composition` body replaced — defensive native collection, red-diamond overlay.
- `plot_sequence_diversity` unchanged.
- New `plot_physicochem` function.
- `main()` calls `plot_physicochem(rows)` after the existing three plots.

CLI shape unchanged: `--receptor-seq` and `--contigs` still accepted (still unused — kept for `.nf` compatibility). The `MPNN_PLOTS` Nextflow process publishes with `path "mpnn_*.png"` so the new physicochem plot is published automatically.

Touched: `mpnn_plots.py` (346 → 550 lines), `pipeline_correct_sequences.py` (748 → 770 lines).

## Lessons re-confirmed this session

### When the plot looks wrong, check the data the plot is reading

The MPNN AA composition story: I went through three rounds of "fixing" the native overlay logic before the user prompted the awk inspection that exposed the upstream bug. The plot was right; the metadata was missing a region. Cost ~30 minutes of plot-code iteration before the right move surfaced.

The general pattern: when a plot looks sparse, off-scale, or otherwise wrong, inspect the *input* with a one-line shell command before editing the plot. The plot's job is to render what it's given; if what it's given is wrong, no amount of plot-code fiddling will make it right.

### Test scaffold + production split worked well

Building `test_<module>_plots.py` + `run_test_<module>_plots_slurm.sh` as standalone iteration tools paid for itself many times over this session. Each plot iteration cycle took seconds (not minutes) because the test script reads cached metrics JSON/CSV instead of re-running the module. Once the user signed off, porting back to production was mechanical: same logic, same imports, plus production conventions (`make_empty_plot`, `save_fallback_plots`, `print(f"Saved ...")`).

Worth doing the same for negsteer (Task 43) and orthogonal-metrics (Task 44) plots — same scaffold, different metrics CSV. Both have csvs already produced by step 1+2 runs.

### Bar-width-vs-axis-width invariants

The design-lengths plot took three rounds because I kept setting one quantity and deriving the other in the wrong direction. The invariant: **one of `{ax_w_in, target_data_range, target_bar_w_in}` is the free variable; the other two are determined by it**. Pick the right free variable depending on what you're trying to keep constant.

In the design-lengths case: `target_bar_w_in = 0.6"` was the user-set constant. `target_data_range` should be derived from it given the data span, and `ax_w_in` then derived from that. When a min-axis-width floor kicks in (so the xlabel doesn't crop), it inflates `ax_w_in` past the natural value — and `target_data_range` has to be re-inflated proportionally to maintain the bar-width invariant. Otherwise the bar gets stretched.

This is a re-occurring pattern in matplotlib layouts where you want consistent physical sizing across panels with different data ranges. Worth keeping in mind for negsteer/orthogonal plot work.

### Variable-length design regions are a real constraint

Several plot ideas were considered and rejected because they required positional alignment across designs (sequence logos, position-by-position AA frequency heatmaps, MSA-style displays). Variable-length regions block all of these. The two plots that ended up in MPNN that *do* work for variable-length data:

- **Composition plots** (AA composition, hydrophobicity, charge) — reduce each design to summary statistics (frequencies, means, sums per residue). Length cancels out in normalisation.
- **Distance-based clustering** (sequence diversity, structural clustering in rfdiff) — pairwise distances handle length differences naturally (Levenshtein for sequences, RMSD-after-trim for coordinates).

Useful framing for thinking about negsteer plot design too — if a metric requires positional alignment, it's probably wrong for this pipeline.

## Outputs of the session

Five files staged at `/mnt/user-data/outputs/` for cluster deployment:

| File | Destination |
|---|---|
| `rfdiffusion_plots.py` | `bin/rfdiffusion_plots.py` |
| `rosetta_filter_plots.py` | `bin/rosetta_filter_plots.py` |
| `mpnn_plots.py` | `bin/mpnn_plots.py` |
| `pipeline_correct_sequences.py` | `bin/pipeline_correct_sequences.py` |
| `test_*_plots.py` + `run_test_*_plots_slurm.sh` (six files total) | `tests/<module>/` |

After deploying:

1. **Re-run `test_proteinmpnn`** so `scored_metadata.csv` regenerates with the fixed `native_residues` column. The existing CSV has the bug baked in and won't be fixed by re-running plots alone. *(User confirmed this step was completed before integration.)*
2. Diagnostic awk on the new CSV should show `<region1>|<region2>` (both regions present). *(User confirmed.)*
3. The MPNN AA composition and physicochem plots should then show native overlays/reference-lines reflecting all design regions, not just region 1.

## Open at session end

### Still to do under "Step 3 — plots refresh"

- **Task 43** P1 — Negative-steering plots from scratch. Inputs: `cross_sequence_summary_with_interface_metrics.csv` (~68 columns). Plot ideas listed in notes 10: tier distribution, composite-score histogram, ra_eff vs ipSAE scatter, contamination rates, control-row diagnostics. **Scoped for a separate conversation** — the existing rfdiff/rosetta/mpnn plots all had previous-session structure to iterate on top of, but negsteer plots are net new and the metric semantics deserve a clean-slate discussion.
- **Task 44** P1 — Orthogonal-metrics plots from scratch. Inputs: `survivors_with_orthogonal_metrics.csv` (~83 columns). Plot ideas in notes 10: AF3 vs Boltz ra_eff scatter, Sc/BSA/interface pLDDT/ΔΔG histograms, weighted_jaccard. Same logic — separate conversation.
- **HADDOCK plots** — explicitly excluded ("HADDOCK work is going to come later on" per user). Not in any current task list.

### Status against the 1–15 plan

| Step | Task | Status |
|---|---|---|
| 1 | Individual section tests | DONE in earlier sessions |
| 2 | Full pipeline run (small) | DONE 26 Apr (notes 10) |
| **3** | **Plots refresh** | **PARTIAL — rfdiff/rosetta/mpnn done this session; negsteer (Task 43) and orthogonal (Task 44) outstanding** |
| 4+ | RFDiffusion parallelisation, multi-cycle steering, ... | Not started |

## Environmental gotchas (re-statement, unchanged from notes 10)

- Login node has stdlib only. Anything importing numpy, gemmi, freesasa, MDAnalysis MUST run inside a container.
- LRR_Pipeline container has TWO pythons: `/usr/bin/python3` (system stub, no scientific stack) and `/usr/local/bin/python3.11` (full stack — use this). Or just call `python` (no version suffix), which is symlinked to the right one.
- Production outdir: `${projectDir}/${params.project_name}_results` — NOT `${params.outdir}` from `params.yml` (config overrides).
- Test results live at `tests/<module>/receptor_resurfacing_results/<module>/`, not `tests/<module>/results/<module>/`.
- For SLURM wrappers, log output names matter: production tests use `slurm_%j.out`, plot iteration tests use `test_plots_%j.out` (avoids collision when both run side-by-side).

## Closing note

This session was scoped narrowly to plot iteration on three already-existing modules. The cumulative diff is ~600 lines added across four production files plus six new test-harness files. The one substantive bug fix (`pipeline_correct_sequences.py` boundary cases) was discovered through plot work but is independent of plotting — it would have surfaced eventually wherever else `native_residues` got consumed.

Negsteer and orthogonal-metrics plots intentionally not started in this session. The user's instinct that they need a fresh conversation is right — both have larger scope (5+ plots each, brand-new code, novel metric semantics to think about) and would have over-stuffed this session if mixed in. Notes 12 will pick up there.
