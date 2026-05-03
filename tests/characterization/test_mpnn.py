"""Stage 4 — ProteinMPNN cohort outputs.

Pins the consolidated outputs of ``modules/proteinmpnn.nf`` and the
sequence-QC/scoring chain that feeds Stage 5. The per-design ``mpnn/`` tree
contains intermediate per-design fixed-positions / FASTA artefacts; for
Wave 1 we pin the cohort-level files at ``sequences/`` (the consolidated
view), plus one representative per-design top-FASTA as a spot check.

The per-module fixture has two parent designs (design_0, design_28). Only
design_0's sequences pass QC and reach ``sequences/qc_fastas/`` and
``sequences/top_fastas/``; design_28's MPNN sequences are present in the
per-design ``mpnn/design_28/`` tree but are absent from the QC/top sets
(they did not survive ``bin/sequence_qc.py``). The representative spot
check therefore stays on ``design_0_seq_0``.

Producers:
- ``MPNN_CLUSTER`` → ``bin/cluster_sequences.py`` (mpnn_cluster_counts.csv)
- ``SEQUENCE_CORRECTION`` → ``bin/correct_mpnn_sequences.py`` (mpnn_corrected.fasta)
- ``SEQUENCE_QC`` → ``bin/sequence_qc.py`` (qc_metadata.csv, qc_report.txt, sequence_metadata.csv)
- ``MPNN_DESIGN_REGION_SCORE`` → ``bin/score_design_region.py`` (scored_metadata.csv)
- ``MPNN_SELECT_TOP`` → ``bin/select_top_sequences.py`` (top_metadata.csv, top_fastas/*.fasta)
"""
from __future__ import annotations

import pytest

from tests.characterization.helpers.csv_compare import compare_csv_exact
from tests.characterization.helpers.text_compare import compare_text_exact


@pytest.fixture
def ref(stage_reference_root):
    return stage_reference_root("proteinmpnn")


@pytest.fixture
def out(stage_output_root):
    return stage_output_root("proteinmpnn")


@pytest.mark.hpc
@pytest.mark.wave1
def test_mpnn_cluster_counts_from_cluster_sequences_py(ref, out):
    """bin/cluster_sequences.py — CSV-EXACT."""
    rel = "sequences/mpnn_cluster_counts.csv"
    compare_csv_exact(ref / rel, out / rel).assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
def test_mpnn_corrected_fasta_from_correct_mpnn_sequences_py(ref, out):
    """bin/correct_mpnn_sequences.py — TEXT-EXACT."""
    rel = "sequences/mpnn_corrected.fasta"
    compare_text_exact(ref / rel, out / rel).assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
def test_qc_metadata_from_sequence_qc_py(ref, out):
    """bin/sequence_qc.py — CSV-EXACT."""
    rel = "sequences/qc_metadata.csv"
    compare_csv_exact(ref / rel, out / rel).assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
def test_qc_report_from_sequence_qc_py(ref, out):
    """bin/sequence_qc.py — TEXT-EXACT."""
    rel = "sequences/qc_report.txt"
    compare_text_exact(ref / rel, out / rel).assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
def test_scored_metadata_from_score_design_region_py(ref, out):
    """bin/score_design_region.py — CSV-EXACT."""
    rel = "sequences/scored_metadata.csv"
    compare_csv_exact(ref / rel, out / rel).assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
def test_sequence_metadata_from_sequence_qc_py(ref, out):
    """bin/sequence_qc.py — CSV-EXACT."""
    rel = "sequences/sequence_metadata.csv"
    compare_csv_exact(ref / rel, out / rel).assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
def test_top_metadata_from_select_top_sequences_py(ref, out):
    """bin/select_top_sequences.py — CSV-EXACT."""
    rel = "sequences/top_metadata.csv"
    compare_csv_exact(ref / rel, out / rel).assert_passed()


@pytest.mark.hpc
@pytest.mark.wave1
def test_design_0_seq_0_top_fasta_from_select_top_sequences_py(ref, out):
    """bin/select_top_sequences.py — TEXT-EXACT (representative spot-check of top_fastas/)."""
    rel = "sequences/top_fastas/design_0_seq_0.fasta"
    compare_text_exact(ref / rel, out / rel).assert_passed()
