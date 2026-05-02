"""Stage 5 (per-sequence) — Negative steering per-run outputs.

Each parametrized case pins the per-sequence subtree
``negative_steering/runs/<sequence>/`` for one of three representative
sequences:

- ``input_control_polyA``   — control (no contact residues / cold-start path)
- ``design_0_seq_0``        — design with contaminated.json populated
                              but reversion attempted
- ``design_13_seq_2``       — design exercising the full reversion harvest
                              path (populated reversion subtree)

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
from tests.characterization.helpers.json_compare import (
    compare_json_deep,
    compare_json_modulo_paths,
)
from tests.characterization.helpers.text_compare import compare_text_exact

SEQUENCE_TRIO = ["input_control_polyA", "design_0_seq_0", "design_13_seq_2"]


def _per_sequence_dir(reference_root, output_root, sequence: str):
    return (
        reference_root / "negative_steering" / "runs" / sequence,
        output_root / "negative_steering" / "runs" / sequence,
    )


def _maybe_skip(ref_path, sequence: str, label: str) -> None:
    if not ref_path.exists():
        pytest.skip(
            f"{label} not present for {sequence!r} in the reference set "
            f"(missing: {ref_path}). Pinning the absence is itself useful information."
        )


# ---------------------------------------------------------------------------
# Per-sequence root files
# ---------------------------------------------------------------------------


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("sequence", SEQUENCE_TRIO)
def test_aggregated_results_from_boltz2_iterate_steering_py(reference_root, output_root, sequence):
    """bin/boltz2_iterate_steering.py — CSV-EXACT (per-sequence aggregated_results.csv)."""
    ref_dir, act_dir = _per_sequence_dir(reference_root, output_root, sequence)
    ref = ref_dir / "aggregated_results.csv"
    _maybe_skip(ref, sequence, "aggregated_results.csv")
    compare_csv_exact(ref, act_dir / "aggregated_results.csv").assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("sequence", SEQUENCE_TRIO)
def test_all_results_multicycle_from_boltz2_iterate_steering_py(reference_root, output_root, sequence):
    """bin/boltz2_iterate_steering.py — CSV-EXACT (all_results_multicycle.csv)."""
    ref_dir, act_dir = _per_sequence_dir(reference_root, output_root, sequence)
    ref = ref_dir / "all_results_multicycle.csv"
    _maybe_skip(ref, sequence, "all_results_multicycle.csv")
    compare_csv_exact(ref, act_dir / "all_results_multicycle.csv").assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("sequence", SEQUENCE_TRIO)
def test_all_results_multicycle_with_metrics_from_boltz2_iterate_steering_py(
    reference_root, output_root, sequence,
):
    """bin/boltz2_iterate_steering.py + interface metrics — CSV-EXACT."""
    ref_dir, act_dir = _per_sequence_dir(reference_root, output_root, sequence)
    ref = ref_dir / "all_results_multicycle_with_metrics.csv"
    _maybe_skip(ref, sequence, "all_results_multicycle_with_metrics.csv")
    compare_csv_exact(ref, act_dir / "all_results_multicycle_with_metrics.csv").assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("sequence", SEQUENCE_TRIO)
def test_cycle_statistics_from_boltz2_iterate_steering_py(reference_root, output_root, sequence):
    """bin/boltz2_iterate_steering.py — CSV-EXACT (cycle_statistics.csv)."""
    ref_dir, act_dir = _per_sequence_dir(reference_root, output_root, sequence)
    ref = ref_dir / "cycle_statistics.csv"
    _maybe_skip(ref, sequence, "cycle_statistics.csv")
    compare_csv_exact(ref, act_dir / "cycle_statistics.csv").assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("sequence", SEQUENCE_TRIO)
def test_passing_summary_from_extract_passing_py(reference_root, output_root, sequence):
    """bin/extract_passing.py — CSV-EXACT.

    Note: ``input_control_polyA``'s passing_summary.csv may be header-only
    (zero rows). The CSV-EXACT comparator handles that case; see the
    `test_compare_csv_exact_empty_body_passes` unit test.
    """
    ref_dir, act_dir = _per_sequence_dir(reference_root, output_root, sequence)
    ref = ref_dir / "passing_summary.csv"
    _maybe_skip(ref, sequence, "passing_summary.csv")
    compare_csv_exact(ref, act_dir / "passing_summary.csv").assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("sequence", SEQUENCE_TRIO)
def test_pathways_from_boltz2_iterate_steering_py(reference_root, output_root, path_normalizer, sequence):
    """bin/boltz2_iterate_steering.py — JSON-MODULO-PATHS (workdir paths in cells)."""
    ref_dir, act_dir = _per_sequence_dir(reference_root, output_root, sequence)
    ref = ref_dir / "pathways.json"
    _maybe_skip(ref, sequence, "pathways.json")
    compare_json_modulo_paths(ref, act_dir / "pathways.json", path_normalizer).assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("sequence", SEQUENCE_TRIO)
def test_raw_per_seed_results_from_reversion_harvest(reference_root, output_root, sequence):
    """bin/reversion.py (harvest) — CSV-EXACT (raw_per_seed_results.csv)."""
    ref_dir, act_dir = _per_sequence_dir(reference_root, output_root, sequence)
    ref = ref_dir / "raw_per_seed_results.csv"
    _maybe_skip(ref, sequence, "raw_per_seed_results.csv")
    compare_csv_exact(ref, act_dir / "raw_per_seed_results.csv").assert_passed()


# ---------------------------------------------------------------------------
# Per-sequence inputs/
# ---------------------------------------------------------------------------


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("sequence", SEQUENCE_TRIO)
def test_inputs_effector_fasta_from_negative_steering_run_one_sh(reference_root, output_root, sequence):
    """bin/negative_steering_run_one.sh staging — TEXT-EXACT (inputs/effector.fasta)."""
    ref_dir, act_dir = _per_sequence_dir(reference_root, output_root, sequence)
    ref = ref_dir / "inputs" / "effector.fasta"
    _maybe_skip(ref, sequence, "inputs/effector.fasta")
    compare_text_exact(ref, act_dir / "inputs" / "effector.fasta").assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("sequence", SEQUENCE_TRIO)
def test_inputs_receptor_fasta_from_negative_steering_run_one_sh(reference_root, output_root, sequence):
    """bin/negative_steering_run_one.sh staging — TEXT-EXACT (inputs/receptor.fasta).

    ``input_control_polyA`` does not have ``inputs/receptor.fasta`` — it
    instead has ``control_polyA_receptor.fasta`` and
    ``control_scrambled_receptor.fasta``, pinned separately by
    :func:`test_control_inputs_fasta_from_negative_steering_run_one_sh`.
    This test pins the file's absence for the control sequence by skipping.
    """
    ref_dir, act_dir = _per_sequence_dir(reference_root, output_root, sequence)
    ref = ref_dir / "inputs" / "receptor.fasta"
    _maybe_skip(ref, sequence, "inputs/receptor.fasta")
    compare_text_exact(ref, act_dir / "inputs" / "receptor.fasta").assert_passed()


# Control-specific FASTAs. The ``input_control_polyA`` run stages additional
# FASTAs that the design sequences do not (the per-control receptor variants
# and the source effector). Pinning them here closes the gap left by the
# control's skip on inputs/receptor.fasta above.
CONTROL_INPUT_FASTAS = [
    "control_polyA_receptor.fasta",
    "control_scrambled_receptor.fasta",
    "source_effector.fasta",
]


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("filename", CONTROL_INPUT_FASTAS)
def test_control_inputs_fasta_from_negative_steering_run_one_sh(
    reference_root, output_root, filename,
):
    """bin/negative_steering_run_one.sh staging — TEXT-EXACT (control-specific FASTA).

    Pins the control-only FASTAs under
    ``negative_steering/runs/input_control_polyA/inputs/``: per-control
    receptor variants (``control_polyA_receptor.fasta``,
    ``control_scrambled_receptor.fasta``) and ``source_effector.fasta``.
    These do not exist for the design sequences.
    """
    rel = f"negative_steering/runs/input_control_polyA/inputs/{filename}"
    ref = reference_root / rel
    if not ref.exists():
        pytest.skip(f"{filename} not present in control inputs reference set: {ref}")
    compare_text_exact(ref, output_root / rel).assert_passed()


# ---------------------------------------------------------------------------
# Per-sequence cycle_0/ files
# ---------------------------------------------------------------------------


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("sequence", SEQUENCE_TRIO)
def test_cycle_0_contaminated_from_reversion_py(reference_root, output_root, path_normalizer, sequence):
    """bin/reversion.py — JSON-MODULO-PATHS (cycle_0/contaminated.json)."""
    ref_dir, act_dir = _per_sequence_dir(reference_root, output_root, sequence)
    ref = ref_dir / "cycle_0" / "contaminated.json"
    _maybe_skip(ref, sequence, "cycle_0/contaminated.json")
    compare_json_modulo_paths(ref, act_dir / "cycle_0" / "contaminated.json", path_normalizer).assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("sequence", SEQUENCE_TRIO)
def test_cycle_0_effector_template_from_extract_effector_template_cif(
    reference_root, output_root, sequence,
):
    """bin/boltz2_negative_steering.py::extract_effector_template_cif — TEXT-EXACT."""
    ref_dir, act_dir = _per_sequence_dir(reference_root, output_root, sequence)
    ref = ref_dir / "cycle_0" / "effector_template.cif"
    _maybe_skip(ref, sequence, "cycle_0/effector_template.cif")
    compare_text_exact(ref, act_dir / "cycle_0" / "effector_template.cif").assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("sequence", SEQUENCE_TRIO)
def test_cycle_0_kickoff_distances_from_boltz2_negative_steering_py(
    reference_root, output_root, path_normalizer, sequence,
):
    """bin/boltz2_negative_steering.py — JSON-MODULO-PATHS.

    Plan §2.5 listed this as JSON-DEEP, but inspection of the reference file
    shows ~120 occurrences of workdir paths in candidate-PDB references; a
    plain JSON-DEEP would fail every HPC round-trip. Demoted to
    JSON-MODULO-PATHS — judgment call, flagged for human review.
    """
    ref_dir, act_dir = _per_sequence_dir(reference_root, output_root, sequence)
    ref = ref_dir / "cycle_0" / "kickoff_distances.json"
    _maybe_skip(ref, sequence, "cycle_0/kickoff_distances.json")
    compare_json_modulo_paths(ref, act_dir / "cycle_0" / "kickoff_distances.json", path_normalizer).assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("sequence", SEQUENCE_TRIO)
def test_cycle_0_passing_from_boltz2_negative_steering_py(
    reference_root, output_root, path_normalizer, sequence,
):
    """bin/boltz2_negative_steering.py — JSON-MODULO-PATHS.

    Plan §2.5 listed JSON-DEEP, but reference contains workdir paths in
    candidate references (~94 occurrences). Demoted — judgment call.
    """
    ref_dir, act_dir = _per_sequence_dir(reference_root, output_root, sequence)
    ref = ref_dir / "cycle_0" / "passing.json"
    _maybe_skip(ref, sequence, "cycle_0/passing.json")
    compare_json_modulo_paths(ref, act_dir / "cycle_0" / "passing.json", path_normalizer).assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("sequence", SEQUENCE_TRIO)
def test_cycle_0_plan_from_boltz2_negative_steering_py(
    reference_root, output_root, path_normalizer, sequence,
):
    """bin/boltz2_negative_steering.py::cmd_plan — JSON-MODULO-PATHS."""
    ref_dir, act_dir = _per_sequence_dir(reference_root, output_root, sequence)
    ref = ref_dir / "cycle_0" / "plan.json"
    _maybe_skip(ref, sequence, "cycle_0/plan.json")
    compare_json_modulo_paths(ref, act_dir / "cycle_0" / "plan.json", path_normalizer).assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("sequence", SEQUENCE_TRIO)
def test_cycle_0_prefilter_from_boltz2_negative_steering_py(
    reference_root, output_root, path_normalizer, sequence,
):
    """bin/boltz2_negative_steering.py — JSON-MODULO-PATHS.

    Plan §2.5 listed JSON-DEEP; reference contains ~84 workdir paths.
    Demoted — judgment call.
    """
    ref_dir, act_dir = _per_sequence_dir(reference_root, output_root, sequence)
    ref = ref_dir / "cycle_0" / "prefilter.json"
    _maybe_skip(ref, sequence, "cycle_0/prefilter.json")
    compare_json_modulo_paths(ref, act_dir / "cycle_0" / "prefilter.json", path_normalizer).assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("sequence", SEQUENCE_TRIO)
def test_cycle_0_reversion_plan_from_reversion_py(
    reference_root, output_root, path_normalizer, sequence,
):
    """bin/reversion.py — JSON-MODULO-PATHS (cycle_0/reversion_plan.json)."""
    ref_dir, act_dir = _per_sequence_dir(reference_root, output_root, sequence)
    ref = ref_dir / "cycle_0" / "reversion_plan.json"
    _maybe_skip(ref, sequence, "cycle_0/reversion_plan.json")
    compare_json_modulo_paths(ref, act_dir / "cycle_0" / "reversion_plan.json", path_normalizer).assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("sequence", SEQUENCE_TRIO)
def test_cycle_0_reversion_results_per_seed_from_reversion_py(
    reference_root, output_root, path_normalizer, sequence,
):
    """bin/reversion.py — JSON-MODULO-PATHS (per-seed harvest)."""
    ref_dir, act_dir = _per_sequence_dir(reference_root, output_root, sequence)
    ref = ref_dir / "cycle_0" / "reversion_results_per_seed.json"
    _maybe_skip(ref, sequence, "cycle_0/reversion_results_per_seed.json")
    compare_json_modulo_paths(
        ref, act_dir / "cycle_0" / "reversion_results_per_seed.json", path_normalizer,
    ).assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("sequence", SEQUENCE_TRIO)
def test_cycle_0_reversion_results_from_reversion_py(
    reference_root, output_root, path_normalizer, sequence,
):
    """bin/reversion.py — JSON-MODULO-PATHS (consolidated harvest)."""
    ref_dir, act_dir = _per_sequence_dir(reference_root, output_root, sequence)
    ref = ref_dir / "cycle_0" / "reversion_results.json"
    _maybe_skip(ref, sequence, "cycle_0/reversion_results.json")
    compare_json_modulo_paths(
        ref, act_dir / "cycle_0" / "reversion_results.json", path_normalizer,
    ).assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("sequence", SEQUENCE_TRIO)
def test_cycle_0_steered_results_from_boltz2_negative_steering_py(
    reference_root, output_root, sequence,
):
    """bin/boltz2_negative_steering.py — CSV-EXACT (steered_results.csv)."""
    ref_dir, act_dir = _per_sequence_dir(reference_root, output_root, sequence)
    ref = ref_dir / "cycle_0" / "steered_results.csv"
    _maybe_skip(ref, sequence, "cycle_0/steered_results.csv")
    compare_csv_exact(ref, act_dir / "cycle_0" / "steered_results.csv").assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("sequence", SEQUENCE_TRIO)
def test_cycle_0_steered_results_aggregate_from_boltz2_negative_steering_py(
    reference_root, output_root, sequence,
):
    """bin/boltz2_negative_steering.py — CSV-EXACT (steered_results_aggregate.csv)."""
    ref_dir, act_dir = _per_sequence_dir(reference_root, output_root, sequence)
    ref = ref_dir / "cycle_0" / "steered_results_aggregate.csv"
    _maybe_skip(ref, sequence, "cycle_0/steered_results_aggregate.csv")
    compare_csv_exact(ref, act_dir / "cycle_0" / "steered_results_aggregate.csv").assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("sequence", SEQUENCE_TRIO)
def test_cycle_0_summary_from_boltz2_iterate_steering_py(reference_root, output_root, sequence):
    """bin/boltz2_iterate_steering.py — TEXT-EXACT (cycle_0/summary.txt)."""
    ref_dir, act_dir = _per_sequence_dir(reference_root, output_root, sequence)
    ref = ref_dir / "cycle_0" / "summary.txt"
    _maybe_skip(ref, sequence, "cycle_0/summary.txt")
    compare_text_exact(ref, act_dir / "cycle_0" / "summary.txt").assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("sequence", SEQUENCE_TRIO)
def test_cycle_0_true_interface_residues_from_derive_input_indices_py(
    reference_root, output_root, sequence,
):
    """bin/derive_input_indices.py — TEXT-EXACT (cycle_0/true_interface_residues.txt)."""
    ref_dir, act_dir = _per_sequence_dir(reference_root, output_root, sequence)
    ref = ref_dir / "cycle_0" / "true_interface_residues.txt"
    _maybe_skip(ref, sequence, "cycle_0/true_interface_residues.txt")
    compare_text_exact(ref, act_dir / "cycle_0" / "true_interface_residues.txt").assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
@pytest.mark.parametrize("sequence", SEQUENCE_TRIO)
def test_cycle_0_wrong_interface_residues_from_derive_input_indices_py(
    reference_root, output_root, sequence,
):
    """bin/derive_input_indices.py — TEXT-EXACT (cycle_0/wrong_interface_residues.txt)."""
    ref_dir, act_dir = _per_sequence_dir(reference_root, output_root, sequence)
    ref = ref_dir / "cycle_0" / "wrong_interface_residues.txt"
    _maybe_skip(ref, sequence, "cycle_0/wrong_interface_residues.txt")
    compare_text_exact(ref, act_dir / "cycle_0" / "wrong_interface_residues.txt").assert_passed()
