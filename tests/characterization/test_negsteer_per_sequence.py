"""Stage 5 (per-sequence) — Negative steering per-run outputs.

Each parametrized case pins the per-sequence subtree
``negative_steering/runs/<sequence>/`` for one of the curated cohort
sequences. The cohort spans the seven observable per-MPNN-sequence outcome
classes the discovery run produced, plus the two synthesised controls that
every per-module test of negative_steering generates from
``params.input_pdb``.

Cohort (see notes/inventory/15_discovery_run_path_coverage.md
§Negative steering "Per-MPNN-sequence outcome classes actually observable"
for the path-coverage rationale):

- ``design_28_seq_1`` — Class 1: cold_start_all_clean (steering skipped)
- ``design_62_seq_0`` — Class 2: steered_clean tiered (rep_design ≠ initial)
- ``design_0_seq_0`` — Class 3a: steered_clean tier-none (empty passing_summary)
- ``design_42_seq_0`` — Class 4: pose_holds 3/3 (tier A)
- ``design_3_seq_1``  — Class 5: pose_holds 2/3 (tier B)
- ``design_55_seq_1`` — Class 6: pose_holds 1/3 (tier C)
- ``design_44_seq_1`` — Class 4 mixed (pose_holds + clean_steered)
- ``design_27_seq_0`` — Class 7: pose_collapses aggregated
- ``input_control_polyA``     — synthesised control (cold-start path only)
- ``input_control_scrambled`` — synthesised control (cold-start path only)

Class 1 sequences (e.g. ``design_28_seq_1``) skip steering entirely, so
``cycle_0/steered/``, ``cycle_0/passing.json``, ``cycle_0/summary.txt``,
``cycle_0/steered_results*.csv``, and ``cycle_statistics.csv`` are all
absent. The ``_maybe_skip`` helper handles the per-sequence absence.

Producers (per ``modules/negative_steering.nf`` and the bin/ scripts it
orchestrates via ``bin/negative_steering_run_one.sh``):
- ``bin/boltz2_negative_steering.py``   (cmd_plan, cmd_collect)
- ``bin/boltz2_iterate_steering.py``    (per-cycle orchestration)
- ``bin/reversion.py``                  (contamination + reversion harvest)
- ``bin/extract_passing.py``            (passing_summary.csv)
- ``bin/cross_sequence_summary.py``     (aggregated_results.csv at per-sequence scope)

Files marked JSON-MODULO-PATHS contain Nextflow workdir hashes that vary
between runs; values are normalised via the ``path_normalizer`` fixture
(``canonicalize_workdir_paths``).
"""
from __future__ import annotations

import pytest

from tests.characterization.helpers.csv_compare import compare_csv_exact
from tests.characterization.helpers.json_compare import compare_json_modulo_paths
from tests.characterization.helpers.text_compare import compare_text_exact

SEQUENCE_COHORT = [
    "design_28_seq_1",
    "design_62_seq_0",
    "design_0_seq_0",
    "design_42_seq_0",
    "design_3_seq_1",
    "design_55_seq_1",
    "design_44_seq_1",
    "design_27_seq_0",
    "input_control_polyA",
    "input_control_scrambled",
]


@pytest.fixture
def ref(stage_reference_root):
    return stage_reference_root("negative_steering")


@pytest.fixture
def out(stage_output_root):
    return stage_output_root("negative_steering")


def _per_sequence_dir(ref_root, out_root, sequence: str):
    return (
        ref_root / "negative_steering" / "runs" / sequence,
        out_root / "negative_steering" / "runs" / sequence,
    )


def _maybe_skip(ref_path, sequence: str, label: str) -> None:
    if not ref_path.exists():
        pytest.skip(
            f"{label} not present for {sequence!r} in the per-module reference set "
            f"(missing: {ref_path}). Pinning the absence is itself useful information."
        )


# ---------------------------------------------------------------------------
# Per-sequence root files
# ---------------------------------------------------------------------------


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("sequence", SEQUENCE_COHORT)
def test_aggregated_results_from_boltz2_iterate_steering_py(ref, out, sequence):
    """bin/boltz2_iterate_steering.py — CSV-EXACT (per-sequence aggregated_results.csv)."""
    ref_dir, act_dir = _per_sequence_dir(ref, out, sequence)
    target = ref_dir / "aggregated_results.csv"
    _maybe_skip(target, sequence, "aggregated_results.csv")
    compare_csv_exact(target, act_dir / "aggregated_results.csv").assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("sequence", SEQUENCE_COHORT)
def test_all_results_multicycle_from_boltz2_iterate_steering_py(ref, out, sequence):
    """bin/boltz2_iterate_steering.py — CSV-EXACT (all_results_multicycle.csv)."""
    ref_dir, act_dir = _per_sequence_dir(ref, out, sequence)
    target = ref_dir / "all_results_multicycle.csv"
    _maybe_skip(target, sequence, "all_results_multicycle.csv")
    compare_csv_exact(target, act_dir / "all_results_multicycle.csv").assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("sequence", SEQUENCE_COHORT)
def test_all_results_multicycle_with_metrics_from_boltz2_iterate_steering_py(
    ref, out, sequence,
):
    """bin/boltz2_iterate_steering.py + interface metrics — CSV-EXACT."""
    ref_dir, act_dir = _per_sequence_dir(ref, out, sequence)
    target = ref_dir / "all_results_multicycle_with_metrics.csv"
    _maybe_skip(target, sequence, "all_results_multicycle_with_metrics.csv")
    compare_csv_exact(target, act_dir / "all_results_multicycle_with_metrics.csv").assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("sequence", SEQUENCE_COHORT)
def test_cycle_statistics_from_boltz2_iterate_steering_py(ref, out, sequence):
    """bin/boltz2_iterate_steering.py — CSV-EXACT (cycle_statistics.csv).

    Class 1 sequences (e.g. ``design_28_seq_1``) skip steering and emit no
    cycle_statistics.csv; ``_maybe_skip`` handles the absence.
    """
    ref_dir, act_dir = _per_sequence_dir(ref, out, sequence)
    target = ref_dir / "cycle_statistics.csv"
    _maybe_skip(target, sequence, "cycle_statistics.csv")
    compare_csv_exact(target, act_dir / "cycle_statistics.csv").assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("sequence", SEQUENCE_COHORT)
def test_passing_summary_from_extract_passing_py(ref, out, sequence):
    """bin/extract_passing.py — CSV-EXACT.

    ``input_control_polyA``'s passing_summary.csv may be header-only
    (zero rows). The CSV-EXACT comparator handles that case; see the
    ``test_compare_csv_exact_empty_body_passes`` unit test.
    """
    ref_dir, act_dir = _per_sequence_dir(ref, out, sequence)
    target = ref_dir / "passing_summary.csv"
    _maybe_skip(target, sequence, "passing_summary.csv")
    compare_csv_exact(target, act_dir / "passing_summary.csv").assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("sequence", SEQUENCE_COHORT)
def test_pathways_from_boltz2_iterate_steering_py(ref, out, path_normalizer, sequence):
    """bin/boltz2_iterate_steering.py — JSON-MODULO-PATHS (workdir paths in cells)."""
    ref_dir, act_dir = _per_sequence_dir(ref, out, sequence)
    target = ref_dir / "pathways.json"
    _maybe_skip(target, sequence, "pathways.json")
    compare_json_modulo_paths(target, act_dir / "pathways.json", path_normalizer).assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("sequence", SEQUENCE_COHORT)
def test_raw_per_seed_results_from_reversion_harvest(ref, out, sequence):
    """bin/reversion.py (harvest) — CSV-EXACT (raw_per_seed_results.csv)."""
    ref_dir, act_dir = _per_sequence_dir(ref, out, sequence)
    target = ref_dir / "raw_per_seed_results.csv"
    _maybe_skip(target, sequence, "raw_per_seed_results.csv")
    compare_csv_exact(target, act_dir / "raw_per_seed_results.csv").assert_passed()


# ---------------------------------------------------------------------------
# Per-sequence inputs/
# ---------------------------------------------------------------------------


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("sequence", SEQUENCE_COHORT)
def test_inputs_effector_fasta_from_negative_steering_run_one_sh(ref, out, sequence):
    """bin/negative_steering_run_one.sh staging — TEXT-EXACT (inputs/effector.fasta)."""
    ref_dir, act_dir = _per_sequence_dir(ref, out, sequence)
    target = ref_dir / "inputs" / "effector.fasta"
    _maybe_skip(target, sequence, "inputs/effector.fasta")
    compare_text_exact(target, act_dir / "inputs" / "effector.fasta").assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("sequence", SEQUENCE_COHORT)
def test_inputs_receptor_fasta_from_negative_steering_run_one_sh(ref, out, sequence):
    """bin/negative_steering_run_one.sh staging — TEXT-EXACT (inputs/receptor.fasta).

    Both controls (``input_control_polyA``, ``input_control_scrambled``) lack
    ``inputs/receptor.fasta`` — they instead carry
    ``control_polyA_receptor.fasta`` and ``control_scrambled_receptor.fasta``,
    pinned separately by ``test_control_inputs_fasta_from_negative_steering_run_one_sh``.
    The skip on those sequences pins the receptor.fasta absence.
    """
    ref_dir, act_dir = _per_sequence_dir(ref, out, sequence)
    target = ref_dir / "inputs" / "receptor.fasta"
    _maybe_skip(target, sequence, "inputs/receptor.fasta")
    compare_text_exact(target, act_dir / "inputs" / "receptor.fasta").assert_passed()


# Control-specific FASTAs. Both ``input_control_polyA`` and
# ``input_control_scrambled`` runs stage the per-control receptor variants
# and ``source_effector.fasta`` under their respective ``inputs/`` dirs.
# We pin the polyA copy as the canonical reference; the scrambled copy is
# byte-identical for the cohort-shared FASTAs and not separately pinned.
CONTROL_INPUT_FASTAS = [
    "control_polyA_receptor.fasta",
    "control_scrambled_receptor.fasta",
    "source_effector.fasta",
]


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("filename", CONTROL_INPUT_FASTAS)
def test_control_inputs_fasta_from_negative_steering_run_one_sh(ref, out, filename):
    """bin/negative_steering_run_one.sh staging — TEXT-EXACT (control-specific FASTA).

    Pins the control-only FASTAs under
    ``negative_steering/runs/input_control_polyA/inputs/``: per-control
    receptor variants (``control_polyA_receptor.fasta``,
    ``control_scrambled_receptor.fasta``) and ``source_effector.fasta``.
    These do not exist for the design sequences.
    """
    rel = f"negative_steering/runs/input_control_polyA/inputs/{filename}"
    target = ref / rel
    if not target.exists():
        pytest.skip(f"{filename} not present in control inputs reference set: {target}")
    compare_text_exact(target, out / rel).assert_passed()


# ---------------------------------------------------------------------------
# Per-sequence cycle_0/ files
# ---------------------------------------------------------------------------


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("sequence", SEQUENCE_COHORT)
def test_cycle_0_contaminated_from_reversion_py(ref, out, path_normalizer, sequence):
    """bin/reversion.py — JSON-MODULO-PATHS (cycle_0/contaminated.json)."""
    ref_dir, act_dir = _per_sequence_dir(ref, out, sequence)
    target = ref_dir / "cycle_0" / "contaminated.json"
    _maybe_skip(target, sequence, "cycle_0/contaminated.json")
    compare_json_modulo_paths(target, act_dir / "cycle_0" / "contaminated.json", path_normalizer).assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("sequence", SEQUENCE_COHORT)
def test_cycle_0_effector_template_from_extract_effector_template_cif(
    ref, out, sequence,
):
    """bin/boltz2_negative_steering.py::extract_effector_template_cif — TEXT-EXACT."""
    ref_dir, act_dir = _per_sequence_dir(ref, out, sequence)
    target = ref_dir / "cycle_0" / "effector_template.cif"
    _maybe_skip(target, sequence, "cycle_0/effector_template.cif")
    compare_text_exact(target, act_dir / "cycle_0" / "effector_template.cif").assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("sequence", SEQUENCE_COHORT)
def test_cycle_0_kickoff_distances_from_boltz2_negative_steering_py(
    ref, out, path_normalizer, sequence,
):
    """bin/boltz2_negative_steering.py — JSON-MODULO-PATHS.

    Plan §2.5 listed this as JSON-DEEP, but inspection of the reference file
    shows ~120 occurrences of workdir paths in candidate-PDB references; a
    plain JSON-DEEP would fail every HPC round-trip. Demoted to
    JSON-MODULO-PATHS — judgment call, flagged for human review.
    """
    ref_dir, act_dir = _per_sequence_dir(ref, out, sequence)
    target = ref_dir / "cycle_0" / "kickoff_distances.json"
    _maybe_skip(target, sequence, "cycle_0/kickoff_distances.json")
    compare_json_modulo_paths(target, act_dir / "cycle_0" / "kickoff_distances.json", path_normalizer).assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("sequence", SEQUENCE_COHORT)
def test_cycle_0_passing_from_boltz2_negative_steering_py(
    ref, out, path_normalizer, sequence,
):
    """bin/boltz2_negative_steering.py — JSON-MODULO-PATHS.

    Plan §2.5 listed JSON-DEEP, but reference contains workdir paths in
    candidate references (~94 occurrences). Demoted — judgment call.

    Class 1 sequences (cold_start_all_clean) skip steering — passing.json
    is absent. ``_maybe_skip`` handles the absence.
    """
    ref_dir, act_dir = _per_sequence_dir(ref, out, sequence)
    target = ref_dir / "cycle_0" / "passing.json"
    _maybe_skip(target, sequence, "cycle_0/passing.json")
    compare_json_modulo_paths(target, act_dir / "cycle_0" / "passing.json", path_normalizer).assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("sequence", SEQUENCE_COHORT)
def test_cycle_0_plan_from_boltz2_negative_steering_py(
    ref, out, path_normalizer, sequence,
):
    """bin/boltz2_negative_steering.py::cmd_plan — JSON-MODULO-PATHS."""
    ref_dir, act_dir = _per_sequence_dir(ref, out, sequence)
    target = ref_dir / "cycle_0" / "plan.json"
    _maybe_skip(target, sequence, "cycle_0/plan.json")
    compare_json_modulo_paths(target, act_dir / "cycle_0" / "plan.json", path_normalizer).assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("sequence", SEQUENCE_COHORT)
def test_cycle_0_prefilter_from_boltz2_negative_steering_py(
    ref, out, path_normalizer, sequence,
):
    """bin/boltz2_negative_steering.py — JSON-MODULO-PATHS.

    Plan §2.5 listed JSON-DEEP; reference contains ~84 workdir paths.
    Demoted — judgment call.
    """
    ref_dir, act_dir = _per_sequence_dir(ref, out, sequence)
    target = ref_dir / "cycle_0" / "prefilter.json"
    _maybe_skip(target, sequence, "cycle_0/prefilter.json")
    compare_json_modulo_paths(target, act_dir / "cycle_0" / "prefilter.json", path_normalizer).assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("sequence", SEQUENCE_COHORT)
def test_cycle_0_reversion_plan_from_reversion_py(
    ref, out, path_normalizer, sequence,
):
    """bin/reversion.py — JSON-MODULO-PATHS (cycle_0/reversion_plan.json)."""
    ref_dir, act_dir = _per_sequence_dir(ref, out, sequence)
    target = ref_dir / "cycle_0" / "reversion_plan.json"
    _maybe_skip(target, sequence, "cycle_0/reversion_plan.json")
    compare_json_modulo_paths(target, act_dir / "cycle_0" / "reversion_plan.json", path_normalizer).assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("sequence", SEQUENCE_COHORT)
def test_cycle_0_reversion_results_per_seed_from_reversion_py(
    ref, out, path_normalizer, sequence,
):
    """bin/reversion.py — JSON-MODULO-PATHS (per-seed harvest)."""
    ref_dir, act_dir = _per_sequence_dir(ref, out, sequence)
    target = ref_dir / "cycle_0" / "reversion_results_per_seed.json"
    _maybe_skip(target, sequence, "cycle_0/reversion_results_per_seed.json")
    compare_json_modulo_paths(
        target, act_dir / "cycle_0" / "reversion_results_per_seed.json", path_normalizer,
    ).assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("sequence", SEQUENCE_COHORT)
def test_cycle_0_reversion_results_from_reversion_py(
    ref, out, path_normalizer, sequence,
):
    """bin/reversion.py — JSON-MODULO-PATHS (consolidated harvest)."""
    ref_dir, act_dir = _per_sequence_dir(ref, out, sequence)
    target = ref_dir / "cycle_0" / "reversion_results.json"
    _maybe_skip(target, sequence, "cycle_0/reversion_results.json")
    compare_json_modulo_paths(
        target, act_dir / "cycle_0" / "reversion_results.json", path_normalizer,
    ).assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("sequence", SEQUENCE_COHORT)
def test_cycle_0_steered_results_from_boltz2_negative_steering_py(
    ref, out, sequence,
):
    """bin/boltz2_negative_steering.py — CSV-EXACT (steered_results.csv).

    Class 1 sequences skip steering; the file exists for them in some
    layouts as a header-only CSV but not in others. ``_maybe_skip``
    handles either shape — present for design_28_seq_1 in the per-module
    fixture; pinned identically.
    """
    ref_dir, act_dir = _per_sequence_dir(ref, out, sequence)
    target = ref_dir / "cycle_0" / "steered_results.csv"
    _maybe_skip(target, sequence, "cycle_0/steered_results.csv")
    compare_csv_exact(target, act_dir / "cycle_0" / "steered_results.csv").assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("sequence", SEQUENCE_COHORT)
def test_cycle_0_steered_results_aggregate_from_boltz2_negative_steering_py(
    ref, out, sequence,
):
    """bin/boltz2_negative_steering.py — CSV-EXACT (steered_results_aggregate.csv).

    Class 1 sequences skip steering and produce no aggregate CSV;
    ``_maybe_skip`` handles the absence.
    """
    ref_dir, act_dir = _per_sequence_dir(ref, out, sequence)
    target = ref_dir / "cycle_0" / "steered_results_aggregate.csv"
    _maybe_skip(target, sequence, "cycle_0/steered_results_aggregate.csv")
    compare_csv_exact(target, act_dir / "cycle_0" / "steered_results_aggregate.csv").assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("sequence", SEQUENCE_COHORT)
def test_cycle_0_summary_from_boltz2_iterate_steering_py(ref, out, sequence):
    """bin/boltz2_iterate_steering.py — TEXT-EXACT (cycle_0/summary.txt).

    Class 1 sequences skip steering and produce no summary.txt;
    ``_maybe_skip`` handles the absence.
    """
    ref_dir, act_dir = _per_sequence_dir(ref, out, sequence)
    target = ref_dir / "cycle_0" / "summary.txt"
    _maybe_skip(target, sequence, "cycle_0/summary.txt")
    compare_text_exact(target, act_dir / "cycle_0" / "summary.txt").assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("sequence", SEQUENCE_COHORT)
def test_cycle_0_true_interface_residues_from_derive_input_indices_py(
    ref, out, sequence,
):
    """bin/derive_input_indices.py — TEXT-EXACT (cycle_0/true_interface_residues.txt).

    Class 1 sequences skip steering; this file is staged only when the
    steering branch runs. ``_maybe_skip`` handles the absence.
    """
    ref_dir, act_dir = _per_sequence_dir(ref, out, sequence)
    target = ref_dir / "cycle_0" / "true_interface_residues.txt"
    _maybe_skip(target, sequence, "cycle_0/true_interface_residues.txt")
    compare_text_exact(target, act_dir / "cycle_0" / "true_interface_residues.txt").assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("sequence", SEQUENCE_COHORT)
def test_cycle_0_wrong_interface_residues_from_derive_input_indices_py(
    ref, out, sequence,
):
    """bin/derive_input_indices.py — TEXT-EXACT (cycle_0/wrong_interface_residues.txt).

    Class 1 sequences skip steering; this file is staged only when the
    steering branch runs. ``_maybe_skip`` handles the absence.
    """
    ref_dir, act_dir = _per_sequence_dir(ref, out, sequence)
    target = ref_dir / "cycle_0" / "wrong_interface_residues.txt"
    _maybe_skip(target, sequence, "cycle_0/wrong_interface_residues.txt")
    compare_text_exact(target, act_dir / "cycle_0" / "wrong_interface_residues.txt").assert_passed()
