"""Inputs for the three non-plot modules of the ProteinMPNN sequence stage.

`bin/mpnn_sequence_qc.py`, `bin/mpnn_design_region_score.py` and
`bin/mpnn_cluster_sequences.py` run in that order inside the SEQUENCE_QC,
MPNN_DESIGN_REGION_SCORE and MPNN_CLUSTER processes of
`modules/proteinmpnn.nf`. Each one reads the previous one's CSV plus some part
of the MPNN run tree, so they share one set of builders rather than three.

Every builder takes the values that decide the answer instead of hiding them
in a default, the same way `fixtures_plots.py` does. A test that wants a
sequence to be 50% identical to its native says so in the call.

Column and key names come from the producers and from the committed reference
set, never from a docstring:

    sequence_metadata.csv   bin/pipeline_correct_sequences.py, and
                            tests/proteinmpnn/example_output_files/sequences/
                            sequence_metadata.csv
    fixed_positions JSONL   tests/proteinmpnn/example_output_files/mpnn/
                            design_0/fixed_positions_design_0.jsonl
    MPNN seqs FASTA         tests/proteinmpnn/example_output_files/mpnn/
                            design_0/design_0_mpnn/seqs/design_0.fa
    af2 FASTA               tests/proteinmpnn/example_output_files/sequences/
                            qc_fastas/design_0_seq_0.fasta

Two of the three modules reach out to something this repository does not
install: ProteinMPNN's model, and the MMseqs2 binary. Both are reached through
an indirection the modules already have — an import off sys.path, and the
`--mmseqs-bin` flag — so `mpnn_module_with_stub` and `fake_mmseqs` supply a
stand-in on disk and let the module's own path find it. Neither patches the
module under test.
"""
from __future__ import annotations

import csv
import importlib.util
import json
import pickle
import stat
import sys
from pathlib import Path
from typing import Mapping, Sequence

REPO = Path(__file__).resolve().parents[2]
BIN = REPO / "bin"
STUB_DIR = Path(__file__).resolve().parent / "stubs"


# ══════════════════════════════════════════════════════════════════════
# sequence_metadata.csv — what SEQUENCE_QC reads
# ══════════════════════════════════════════════════════════════════════

# The header bin/pipeline_correct_sequences.py writes. It is
# fixtures_plots.MPNN_METADATA_COLUMNS without design_region_score, which the
# later MPNN_DESIGN_REGION_SCORE stage appends.
SEQUENCE_METADATA_COLUMNS = [
    "design", "seq", "mpnn_score", "native_residues", "designed_residues",
    "num_changes", "design_region_length_observed", "design_region_spec",
    "corrected_receptor", "receptor_total_length",
]


def sequence_metadata_csv(path: Path, rows: Sequence[Mapping],
                          columns: Sequence[str] | None = None) -> Path:
    """Write a metadata CSV, header-only when `rows` is empty.

    `columns` defaults to the SEQUENCE_QC input header. A test driving
    MPNN_DESIGN_REGION_SCORE passes the qc_metadata.csv header instead, which
    is the same thing — SEQUENCE_QC copies its input columns through.
    """
    fieldnames = list(columns) if columns is not None else SEQUENCE_METADATA_COLUMNS
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames,
                                extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return path


def qc_row(
    *,
    design: str = "0",
    seq: str = "0",
    mpnn_score: str | float = 1.0,
    corrected_receptor: str = "MKTAYIAKQRQISFVKSHFSRQ",
    designed_residues: str = "AAKK",
    native_residues: str = "AADD",
) -> dict:
    """One row of sequence_metadata.csv.

    `corrected_receptor` is the whole receptor chain, which is what the
    poly-X and composition predicates run on. Its default clears every
    predicate: no run longer than two, no residue above 18% of the chain,
    nothing outside the twenty canonical amino acids.

    `designed_residues` and `native_residues` are the pipe-joined per-region
    strings, and only their pairwise identity matters to QC, so a call that
    wants 50% identity states two four-residue strings differing in two
    places. That is the default, deliberately off both ends of the default
    0-100 identity window: an all-alanine design against an all-alanine
    native would be 100% identical and so rejected by the default
    --max-pct-identity.
    """
    designed_flat = designed_residues.replace("|", "")
    native_flat = native_residues.replace("|", "")
    return {
        "design": design,
        "seq": seq,
        "mpnn_score": str(mpnn_score),
        "native_residues": native_residues,
        "designed_residues": designed_residues,
        "num_changes": str(sum(1 for a, b in zip(designed_flat, native_flat)
                               if a != b)),
        "design_region_length_observed": str(len(designed_flat)),
        "design_region_spec": f"{len(designed_flat)}-{len(designed_flat)}",
        "corrected_receptor": corrected_receptor,
        "receptor_total_length": str(len(corrected_receptor)),
    }


def af2_fasta(path: Path, *, receptor: str | None,
              effector: str = "AIDLSRERDPNFFDNADIPV") -> Path:
    """The two-record FASTA that AF2_FASTAS writes, named design_D_seq_S.fasta.

    `receptor=None` writes the effector record only. That is the input where
    `mpnn_sequence_qc.main` finds no header containing "receptor" and falls
    back to the metadata CSV's corrected_receptor column, so the fallback is
    reachable without hand-writing a malformed file.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    records = []
    if receptor is not None:
        records.append(f">receptor\n{receptor}\n")
    records.append(f">effector\n{effector}\n")
    path.write_text("".join(records))
    return path


# ══════════════════════════════════════════════════════════════════════
# The MPNN run tree — what MPNN_DESIGN_REGION_SCORE reads
# ══════════════════════════════════════════════════════════════════════

def mpnn_seqs_fasta(root: Path, design: str, *, native: str,
                    designed: Sequence[str]) -> Path:
    """ProteinMPNN's own output FASTA at <root>/design_D_mpnn/seqs/design_D.fa.

    The first record is the native sequence and carries the run's provenance
    header; the sampled sequences follow. `score_sequences_for_design` drops
    record 0 and numbers the rest from 0, so `designed[0]` is the sequence
    metadata row seq=0 refers to.

    Chains are "/"-separated, as ProteinMPNN writes them.
    """
    name = f"design_{design}"
    path = root / f"{name}_mpnn" / "seqs" / f"{name}.fa"
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        f">{name}, score=3.8529, global_score=2.0076, fixed_chains=[], "
        "designed_chains=['A', 'B'], model_name=v_48_020, seed=42",
        native,
    ]
    for index, sequence in enumerate(designed, start=1):
        lines.append(f">T=0.25, sample={index}, score=1.3460, "
                     "global_score=1.6773, seq_recovery=0.1500")
        lines.append(sequence)
    path.write_text("\n".join(lines) + "\n")
    return path


def fixed_positions_jsonl(root: Path, design: str,
                          fixed: Mapping[str, Sequence[int]]) -> Path:
    """<root>/fixed_positions_design_D.jsonl, one JSON object on one line.

    Positions are 1-based within their own chain, which is what ProteinMPNN's
    driver writes and what `load_fixed_positions` hands straight to
    `tied_featurize`. Every position NOT listed here is a free, designable
    position, and the design-region score is the mean over exactly those.
    """
    path = root / f"fixed_positions_design_{design}.jsonl"
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {f"design_{design}": {chain: list(positions)
                                    for chain, positions in fixed.items()}}
    path.write_text(json.dumps(payload) + "\n")
    return path


def mpnn_checkpoint(path: Path, *, noise_level: float | None = 0.2,
                    num_edges: int | None = 48) -> Path:
    """A checkpoint carrying the three keys `load_model` reads.

    Written with pickle rather than `torch.save`, because the stub torch
    beside the ProteinMPNN stub reads it with `pickle.load`. `load_model` only
    ever indexes the mapping that comes back, so the archive format is not
    under test.

    Passing None for either hyperparameter leaves the key out, which is the
    input where `load_model`'s `.get` default applies instead.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    payload: dict = {"model_state_dict": {"stub_layer.weight": [0.0]}}
    if noise_level is not None:
        payload["noise_level"] = noise_level
    if num_edges is not None:
        payload["num_edges"] = num_edges
    with open(path, "wb") as handle:
        pickle.dump(payload, handle)
    return path


def mpnn_module_with_stub():
    """Import bin/mpnn_design_region_score.py with HAS_MPNN True.

    The module's `import torch` and `from protein_mpnn_utils import ...` are
    plain top-level imports off sys.path, and off the cluster at least one of
    them raises: CI installs neither, and a laptop that has torch has a torch
    whose `Tensor.numpy()` does not work against numpy 2.x. HAS_MPNN then goes
    False and `main()` exits 1 before any arithmetic runs.

    Putting tests/characterization/stubs/ at the front of sys.path for the
    duration of one `exec_module` makes both imports resolve to the stubs
    there, with no patching of the module under test. Any real torch already
    imported is set aside first and restored afterwards, so the shim is used
    whether or not torch is installed and the result does not depend on what
    an earlier test imported.

    The module is loaded under its own name each call, so a test's
    monkeypatching of MPNN_DIR or `find_checkpoint` cannot leak into the next
    one. Coverage still attributes the lines to bin/mpnn_design_region_score.py,
    which is keyed on the file rather than the module name. The stubs are
    dropped from sys.modules afterwards so no other test can import them by
    accident; the loaded module keeps its own references, so it goes on
    working.
    """
    shadowed = [name for name in list(sys.modules)
                if name == "torch" or name.startswith("torch.")
                or name == "protein_mpnn_utils"]
    saved = {name: sys.modules.pop(name) for name in shadowed}
    sys.path.insert(0, str(STUB_DIR))
    try:
        spec = importlib.util.spec_from_file_location(
            "mpnn_design_region_score__stubbed",
            BIN / "mpnn_design_region_score.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    finally:
        sys.path.remove(str(STUB_DIR))
        sys.modules.pop("protein_mpnn_utils", None)
        sys.modules.pop("mpnn_design_region_score__stubbed", None)
        for name in [n for n in list(sys.modules)
                     if n == "torch" or n.startswith("torch.")]:
            del sys.modules[name]
        sys.modules.update(saved)
    return module


# ══════════════════════════════════════════════════════════════════════
# A stand-in MMseqs2 — what MPNN_CLUSTER shells out to
# ══════════════════════════════════════════════════════════════════════

_FAKE_MMSEQS_BODY = '''\
"""Stand-in for `mmseqs easy-cluster`, driven by a plan file.

Writes <prefix>_cluster.tsv in MMseqs2's own two-column rep/member format and
records every invocation, so a test can assert both the cluster count the
module derives and the command line it built.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PLAN = json.load(open(os.path.join(HERE, "plan.json")))

argv = sys.argv[1:]
fasta, prefix = argv[1], argv[2]
min_seq_id = argv[argv.index("--min-seq-id") + 1]

names = [line[1:].strip() for line in open(fasta) if line.startswith(">")]
with open(os.path.join(HERE, "calls.jsonl"), "a") as log:
    log.write(json.dumps({"argv": argv, "min_seq_id": min_seq_id,
                          "names": names}) + "\\n")

action = PLAN.get(min_seq_id, len(names))
if action == "fail":
    sys.stderr.write("stub mmseqs: refusing threshold %s\\n" % min_seq_id)
    sys.exit(1)
if action == "no_tsv":
    sys.exit(0)

# `action` representatives, with every remaining sequence a member of the
# first one.  Representatives therefore repeat across lines, which is what
# makes counting DISTINCT first columns the thing under test.
n_reps = min(int(action), len(names))
with open(prefix + "_cluster.tsv", "w") as out:
    for index, name in enumerate(names):
        rep = name if index < n_reps else names[0]
        out.write("%s\\t%s\\n" % (rep, name))
'''


def fake_mmseqs(root: Path, plan: Mapping[str, object] | None = None) -> Path:
    """Write an executable stand-in for the mmseqs binary; return its path.

    `plan` maps the `--min-seq-id` string the module formats ("0.50", "1.00")
    to what should happen at that threshold: an int number of clusters to
    emit, "fail" to exit non-zero the way a real mmseqs does on bad input, or
    "no_tsv" to exit 0 having written nothing. An unlisted threshold emits one
    cluster per sequence.

    The result is a /bin/sh wrapper around a Python script, invoked through
    the interpreter running the tests, so the fixture does not depend on
    `python3` being on PATH.
    """
    root.mkdir(parents=True, exist_ok=True)
    script = root / "fake_mmseqs.py"
    script.write_text(_FAKE_MMSEQS_BODY)
    (root / "plan.json").write_text(json.dumps(dict(plan or {})))

    wrapper = root / "mmseqs"
    wrapper.write_text(f'#!/bin/sh\nexec "{sys.executable}" "{script}" "$@"\n')
    wrapper.chmod(wrapper.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP
                  | stat.S_IXOTH)
    return wrapper


def mmseqs_calls(root: Path) -> list[dict]:
    """Every invocation `fake_mmseqs` recorded, in order."""
    log = root / "calls.jsonl"
    if not log.exists():
        return []
    return [json.loads(line) for line in log.read_text().splitlines() if line]
