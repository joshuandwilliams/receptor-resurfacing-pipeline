# Prompt 3 — HPC-tier characterization tests

Run from the repo root, with the `receptor-tests` mamba environment active. This produces every Wave 1 hpc-tier characterization test in a single pass.

The plan document at `notes/inventory/11_phase_2_plan.md` is the authoritative reference. Read §2 (scope of comparison) and §3 (test framework architecture) before writing anything.

---

```
Continuing Phase 2 Wave 1. Prompts 1 and 2 produced the framework
foundation and the comparator unit tests. This prompt produces the
hpc-tier characterization tests — the actual safety net.

Read these in order before writing:
  1. notes/inventory/11_phase_2_plan.md §2 (file-by-file scope) and
     §3 (framework architecture).
  2. tests/characterization/README.md for the framework conventions.
  3. tests/characterization/TRACEABILITY.md for the schema you must
     populate.
  4. tests/full_test_run/example_output_files/ — the actual reference
     data on disk. The structure mirrors a fresh pipeline run's
     results/ directory exactly. Use `ls` and `find` to confirm
     which files exist before writing tests for them. Do NOT pin
     files that aren't in the reference set.

GENERAL CONVENTIONS:

- Tier all tests with @pytest.mark.hpc.
- Add @pytest.mark.wave1 alongside (selection marker, orthogonal to
  tier).
- Test naming: `test_<short_output_description>_from_<producer>`.
  Producer is the .py or .nf file/function that emits the output.
  Examples:
    test_processed_contigs_from_rfdiffusion_contigs_py
    test_cross_sequence_summary_from_aggregate
- Each test has a one-line docstring noting the producer's full
  path and the comparison strategy.
- Resolve paths via reference_root and output_root fixtures, never
  hardcoded.
- Use compare_csv_exact, compare_csv_struct, compare_json_deep,
  compare_json_modulo_paths, compare_png_perceptual, or
  compare_png_exists from tests.characterization.helpers as
  appropriate.
- Call .assert_passed() to raise on mismatch; pytest will show the
  comparator's formatted failure message.
- Default tolerances are 1e-6 absolute / 1e-9 relative — passed
  implicitly via the comparator defaults. Don't override unless
  you have a specific reason.
- For per-sequence tests, use @pytest.mark.parametrize over the
  trio (input_control_polyA, design_0_seq_0, design_13_seq_2). Make
  the list a module-level constant SEQUENCE_TRIO so it's shared
  across tests and visible at the top of each per-sequence file.

WHAT TO PIN — by stage. Eight test files, one per pipeline stage.
For each stage I list the files to pin, the comparator, and any
notes. Verify each file exists in the reference set before adding
its test; if a file isn't there, skip that test and note it in the
final summary.

==========================================================================
1. tests/characterization/test_preprocessing.py
==========================================================================

Files in tests/full_test_run/example_output_files/preprocessing/:
  - processed_contigs.txt → CSV-EXACT not applicable (plain text);
    use a simple compare via the EXACT-equality file comparison
    (read both files, assert text equality). The processed_contigs
    output is not strictly CSV but its content is line-deterministic.
    If comparator helpers don't have a plain-text comparator, write
    a one-off helper inline in this file:

      def _compare_text_exact(ref, act, strategy="TEXT-EXACT"):
          # returns ComparisonResult with text-byte equality
          ...

    Producer: bin/rfdiffusion_contigs.py (RESOLVE_CONTIGS in
    modules/preprocessing.nf)

  - sequences.json → JSON-DEEP
    Producer: EXTRACT_SEQUENCES inline Python in modules/preprocessing.nf

==========================================================================
2. tests/characterization/test_rfdiffusion.py
==========================================================================

Files in tests/full_test_run/example_output_files/rfdiffusion/:
  - rfdiffusion_metrics.json → JSON-DEEP
    Producer: bin/rfdiffusion_filter.py (RFDIFFUSION_FILTER in
    modules/rfdiffusion.nf)
    Note: contains design_region_coords[][3] — RFDiffusion is
    non-deterministic across runs unless seeded. If first HPC
    round-trip fails on this file, the team will demote to
    JSON-STRUCT (keys/types/ranges only). For this prompt: pin
    JSON-DEEP and let the round-trip tell us.

  - filter_summary.json → JSON-DEEP
    Same producer.

  - passing_designs.txt → text-exact
    Same producer.

==========================================================================
3. tests/characterization/test_rosetta_filtering.py
==========================================================================

Files in tests/full_test_run/example_output_files/rosetta_filtering/:
  - rosetta_filter_metrics.json → JSON-DEEP
    Producer: bin/rosetta_filter.py (ROSETTA_FILTER in
    modules/rosetta_filtering.nf)
  - rosetta_filter_summary.json → JSON-DEEP
    Same producer.
  - rosetta_passing_designs.txt → text-exact
    Same producer.

==========================================================================
4. tests/characterization/test_mpnn.py
==========================================================================

Files in tests/full_test_run/example_output_files/mpnn/:
  Walk the directory. Each design_X subdirectory should contain
  design-level outputs (JSON metrics, etc.). For the cohort-level
  Wave 1 we don't pin per-design subtrees individually — instead,
  pin any cohort-level files at mpnn/ root (if present), and rely
  on the sequences/ directory (Stage 5) for the consolidated per-
  design outputs.

  Files at sequences/ root:
  - mpnn_cluster_counts.csv → CSV-EXACT
    Producer: bin/cluster_sequences.py (MPNN_CLUSTER in
    modules/proteinmpnn.nf)
  - mpnn_corrected.fasta → text-exact
    Producer: bin/correct_mpnn_sequences.py (SEQUENCE_CORRECTION)
  - qc_metadata.csv → CSV-EXACT
    Producer: bin/sequence_qc.py (SEQUENCE_QC)
  - qc_report.txt → text-exact
    Same producer.
  - scored_metadata.csv → CSV-EXACT
    Producer: bin/score_design_region.py (MPNN_DESIGN_REGION_SCORE)
  - sequence_metadata.csv → CSV-EXACT
    Producer: bin/sequence_qc.py
  - top_metadata.csv → CSV-EXACT
    Producer: bin/select_top_sequences.py (MPNN_SELECT_TOP)

  Verify the design_*.fasta files in sequences/ — these are the
  per-design top sequences. Pin one (design_0_seq_0.fasta) as a
  representative text-exact comparison; spot-checking any one is
  sufficient since the full set would be redundant with the
  metadata CSVs.

==========================================================================
5. tests/characterization/test_negsteer_cohort.py
==========================================================================

Files in tests/full_test_run/example_output_files/negative_steering/
(top level — NOT inside runs/):

  - cross_sequence_summary.csv → CSV-EXACT
    Producer: bin/cross_sequence_summary.py::aggregate
    (NEGSTEER_CROSS_SEQUENCE in modules/negative_steering.nf)
  - cross_sequence_summary_with_interface_metrics.csv → CSV-EXACT
    Producer: bin/compute_interface_metrics.py
    (NEGSTEER_INTERFACE_METRICS in modules/negsteer_interface_metrics.nf)

  Files in negative_steering/controls_inputs/ and
  negative_steering/indices/:
  - Walk these directories and pin any small text/JSON files
    present. These are the input-derivation outputs. Use
    text-exact for .txt, JSON-DEEP for .json.
    Producer for indices/: bin/derive_input_indices.py
    Producer for controls_inputs/: bin/derive_input_indices.py
    (DERIVE_INPUT_INDICES in modules/negsteer_controls.nf)

==========================================================================
6. tests/characterization/test_negsteer_per_sequence.py
==========================================================================

Use @pytest.mark.parametrize over SEQUENCE_TRIO at module top:

  SEQUENCE_TRIO = ["input_control_polyA", "design_0_seq_0", "design_13_seq_2"]

For each sequence, the per-sequence directory is
  negative_steering/runs/<sequence>/

Pin the following files per sequence. Each is a separate
parametrized test:

Per-sequence root files (top-level of the sequence dir):
  - aggregated_results.csv → CSV-EXACT
  - all_results_multicycle.csv → CSV-EXACT
  - all_results_multicycle_with_metrics.csv → CSV-EXACT
  - cycle_statistics.csv → CSV-EXACT
  - passing_summary.csv → CSV-EXACT (note: input_control_polyA may
    have an empty body — header-only — verify your comparator
    handles that case correctly; it should, per the unit tests)
  - pathways.json → JSON-DEEP (note: contains workdir paths; use
    JSON-MODULO-PATHS instead with the path_normalizer fixture)
  - raw_per_seed_results.csv → CSV-EXACT
  - run_one_runtime_sec.txt → SKIP (timing-dependent; will not
    reproduce. Do not pin.)
  - inputs/effector.fasta → text-exact
  - inputs/receptor.fasta → text-exact

Per-sequence cycle_0/ root files:
  - cycle_0/contaminated.json → JSON-MODULO-PATHS (contains workdir
    paths in design_workdir field; verified for design_13_seq_2)
  - cycle_0/effector_template.cif → text-exact
  - cycle_0/kickoff_distances.json → JSON-DEEP
  - cycle_0/passing.json → JSON-DEEP
  - cycle_0/plan.json → JSON-MODULO-PATHS (likely contains paths)
  - cycle_0/prefilter.json → JSON-DEEP
  - cycle_0/reversion_plan.json → JSON-MODULO-PATHS (paths)
  - cycle_0/reversion_results_per_seed.json → JSON-MODULO-PATHS (paths)
  - cycle_0/reversion_results.json → JSON-MODULO-PATHS (paths)
  - cycle_0/steered_results.csv → CSV-EXACT
  - cycle_0/steered_results_aggregate.csv → CSV-EXACT
  - cycle_0/summary.txt → text-exact
  - cycle_0/true_interface_residues.txt → text-exact
  - cycle_0/wrong_interface_residues.txt → text-exact

For each parametrized test, if the expected file does NOT exist
for a particular sequence (e.g., the control may not have all of
these), use pytest.skip() with a message explaining which sequence
lacks which file. Do not let a missing-file case fail the whole
test; the pinning of "this file doesn't exist for this sequence"
is itself useful information.

For JSON-MODULO-PATHS comparisons: depend on the path_normalizer
fixture and pass it as the second positional arg to
compare_json_modulo_paths.

Producer for the negative-steering core: bin/negative_steering_run_one.sh
orchestrating multiple bin/*.py scripts (boltz2_negative_steering.py,
boltz2_iterate_steering.py, reversion.py, etc.) — see plan §2.5
and 04_functional_categorization.md for the full producer list.

==========================================================================
7. tests/characterization/test_orthogonal_metrics.py
==========================================================================

Files in tests/full_test_run/example_output_files/orthogonal_metrics/:
  - af3_nomsa_summary.csv → CSV-EXACT
    Producer: bin/parse_af3_output.py (AF3_PARSE_OUTPUT in
    modules/negsteer_af3_nomsa.nf)
  - biophysical_summary.csv → CSV-EXACT
    Producer: bin/run_biophysical_metrics.py
    (NEGSTEER_BIOPHYSICAL_METRICS in modules/negsteer_biophysical_metrics.nf)
  - rosetta_summary.csv → CSV-EXACT
    Producer: bin/run_rosetta_metrics.py
    (NEGSTEER_ROSETTA_METRICS in modules/negsteer_rosetta_metrics.nf)
  - survivor_manifest.csv → CSV-EXACT
    Producer: bin/extract_survivor_manifest.py
    (EXTRACT_SURVIVOR_MANIFEST in modules/negsteer_manifest.nf)
    Note: contains paths; use compare_csv_exact's string_columns
    parameter to flag the path columns as strings, OR if there are
    workdir-hash paths in cells, switch to a custom comparison —
    inspect the file before writing the test.
  - survivors_with_orthogonal_metrics.csv → CSV-EXACT
    Producer: bin/merge_orthogonal_metrics.py
    (NEGSTEER_ORTHOGONAL_METRICS in
    modules/negsteer_orthogonal_metrics.nf)

==========================================================================
8. tests/characterization/test_plots.py
==========================================================================

Files in tests/full_test_run/example_output_files/plots/:

  All plots get PNG-PERCEPTUAL with default ssim_threshold=0.95
  and size_tolerance_pct=5.0. Real matplotlib outputs are typically
  large enough that the size tolerance isn't an issue (the
  comparator's docstring covers this).

  Walk the plots/ directory and write one test per .png file
  found. Producers vary by plot category — broadly:

    rfdiff_*.png         → RFDIFFUSION_PLOTS / bin/rfdiffusion_plots.py
    rosetta_*.png        → ROSETTA_FILTER_PLOTS / bin/rosetta_filter_plots.py
    mpnn_*.png           → MPNN_PLOTS / bin/mpnn_plots.py
    negsteer_*.png       → NEGSTEER_PLOTS or NEGSTEER_WITHIN_SEQUENCE_PLOTS
    orthogonal_*.png     → ORTHOG_PLOTS / bin/orthogonal_plots.py

  Don't try to pin per-plot SSIM thresholds yet. The default 0.95
  is the starting calibration; the first HPC round-trip will tell
  us which (if any) plots need per-plot overrides.

==========================================================================
SHARED HELPERS
==========================================================================

If a plain-text-exact comparator doesn't already exist in
helpers/, add one to a new helpers/text_compare.py:

  def compare_text_exact(reference: Path, actual: Path) -> ComparisonResult:
      """Strategy: TEXT-EXACT. Byte-equal comparison of two text files."""
      ...

Add corresponding unit tests to
helpers/tests/test_text_compare.py covering: identity-passes,
single-byte-difference-fails, missing-file-fails. Mark these
@pytest.mark.local_unit. They should be small (~5 tests).

==========================================================================
TRACEABILITY.md
==========================================================================

For every test you add (each parametrize case counts separately),
add a row to tests/characterization/TRACEABILITY.md per the
schema in that file. The Notes column is for things like "JSON-
MODULO-PATHS because of workdir hashes" or "skipped for control
sequence — file not present."

==========================================================================
ACCEPTANCE CRITERIA
==========================================================================

  - All test files exist under tests/characterization/ with the
    naming above.
  - `pytest -m local_unit` still passes (no regressions in the
    comparator unit tests, including new text_compare tests if
    added).
  - `pytest --collect-only -m hpc` succeeds and reports the
    expected test count. RECEPTOR_OUTPUT_ROOT is unset, so all
    hpc-tier tests will be skipped — that's fine; we're checking
    collection, not execution.
  - `pytest -m hpc` (with RECEPTOR_OUTPUT_ROOT unset) reports all
    hpc tests as skipped, with the skip message clearly indicating
    why (env var not set).
  - TRACEABILITY.md has one row per parametrized test instance.
  - Each test file has a module-level docstring naming the stage
    it covers and the producer modules involved.

DO NOT run git add, git commit, or git push. Staging is the human's
job.

In the final summary, list:
  - The number of test functions and the total parametrized test
    count.
  - Any files in the plan that you skipped because they weren't on
    disk (so the human can decide whether the absence is expected).
  - Any tests where the comparison strategy was non-obvious and
    you made a judgment call (so the human can review).
```
