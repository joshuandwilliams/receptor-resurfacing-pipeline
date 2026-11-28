"""Build run trees and structures whose answers are known before the code runs.

Testing against a real pipeline output only proves the code still agrees with
its own last output. These builders state the answer first, then construct the
input that must produce it, so a wrong result is a failure rather than a new
blessed reference.

Three kinds of builder:

`complex_pdb` places a receptor and an effector at a geometry you specify, so
a contact count, a chain split and a centre-of-mass displacement are all known
in advance.

`rfdiffusion_tree` writes the directory layout an RFDiffusion stage produces,
with the metrics JSON and the design PDBs agreeing with each other, so the
readers and the filters can be driven end to end without a cluster.

`negsteer_workdir` and `boltz_prediction_dir` write the per-sequence workdir
the negative-steering engine leaves behind. The engine itself is an installed
Singularity container and is not part of this repository, so its output tree
on disk is the only stand-in `NegativeSteeringRun.from_workdir` needs.

Layouts are taken from the committed reference sets under
`tests/<module>/example_output_files/`, not from docstrings.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Mapping, Sequence

# One CA per residue is enough for every metric in rfdiffusion_filter, which
# reads only CA records. Writing full backbones would make the fixtures larger
# without changing a single assertion.
_ATOM = ("ATOM  {serial:5d}  CA  {res} {chain}{resnum:4d}    "
         "{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C")


def ca_line(serial: int, chain: str, resnum: int, x: float, y: float,
            z: float, res: str = "ALA") -> str:
    return _ATOM.format(serial=serial, res=res, chain=chain, resnum=resnum,
                        x=x, y=y, z=z)


def complex_pdb(
    path: Path,
    *,
    n_receptor: int = 20,
    n_effector: int = 10,
    receptor_chain: str = "A",
    effector_chain: str = "B",
    spacing: float = 3.8,
    separation: float = 4.0,
    receptor_start: int = 1,
    effector_start: int = 1,
    single_chain: bool = False,
    contact_residues: tuple[int, ...] | None = None,
) -> Path:
    """Two chains on parallel lines, `separation` angstroms apart in z.

    The receptor runs along x at z = 0. The effector runs along x at
    z = `separation`, directly above receptor residues 0..n_effector-1 unless
    `contact_residues` says otherwise.

    So every effector residue is exactly `separation` from the receptor
    residue beneath it, and sqrt(separation^2 + (k*spacing)^2) from the one k
    places away. That makes the contact count at a given cutoff arithmetic
    rather than something to measure.

    `single_chain=True` writes both segments under the receptor chain id with
    continuous numbering, which is what an unsplit RFDiffusion output looks
    like. The geometric gap at the join is what the break detector must find.

    `contact_residues` places the effector above those receptor indices only,
    with the rest of it moved far away in y, for building a specific interface.
    """
    lines, serial = [], 1
    for i in range(n_receptor):
        lines.append(ca_line(serial, receptor_chain, receptor_start + i,
                             i * spacing, 0.0, 0.0))
        serial += 1

    eff_chain = receptor_chain if single_chain else effector_chain
    eff_first = (receptor_start + n_receptor) if single_chain else effector_start
    for j in range(n_effector):
        if contact_residues is None:
            x, y = j * spacing, 0.0
        elif j < len(contact_residues):
            x, y = contact_residues[j] * spacing, 0.0
        else:
            x, y = j * spacing, 500.0          # far enough to contact nothing
        lines.append(ca_line(serial, eff_chain, eff_first + j, x, y, separation))
        serial += 1

    lines.append("END")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")
    return path


def lattice_pdb(
    path: Path,
    *,
    chain: str = "A",
    origin: tuple[float, float, float] = (0.0, 0.0, 0.0),
    shape: tuple[int, int, int] = (3, 3, 3),
    step: float = 4.0,
    start: int = 1,
    res: str = "ALA",
    append: bool = False,
) -> Path:
    """One CA per point of an axis-aligned cuboid lattice.

    `complex_pdb` places each chain on a straight line, which is the right
    shape for a contact count and the wrong shape for anything that takes a
    convex hull. Qhull rejects a collinear point set, so the pose solver's
    hull depth, validity and interpenetration terms cannot be driven by it.

    This builder gives a non-degenerate point cloud whose geometry is still
    arithmetic. Point (i, j, k) sits at `origin + (i, j, k) * step` and
    carries residue number `start + i*ny*nz + j*nz + k`, so a residue number
    and a coordinate convert into each other without reading the file.

    The hull of a full lattice is the cuboid
    `[origin, origin + (nx-1, ny-1, nz-1) * step]`, and interior points are
    inside it by a known depth.

    `append=True` adds a second chain to an existing file, which is how a
    two-chain complex with two hulls is built.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    kept: list[str] = []
    if append and path.exists():
        kept = [ln for ln in path.read_text().splitlines()
                if ln.startswith("ATOM")]

    lines, serial, resnum = [], len(kept) + 1, start
    for i in range(shape[0]):
        for j in range(shape[1]):
            for k in range(shape[2]):
                lines.append(ca_line(serial, chain, resnum,
                                     origin[0] + i * step,
                                     origin[1] + j * step,
                                     origin[2] + k * step, res=res))
                serial += 1
                resnum += 1

    path.write_text("\n".join(kept + lines) + "\nEND\n")
    return path


def lattice_resnum(i: int, j: int, k: int, *,
                   shape: tuple[int, int, int] = (3, 3, 3),
                   start: int = 1) -> int:
    """The residue number `lattice_pdb` gives to lattice point (i, j, k)."""
    _, ny, nz = shape
    return start + i * ny * nz + j * nz + k


def rfdiffusion_tree(
    root: Path,
    *,
    designs: dict[int, dict] | None = None,
    n_receptor: int = 20,
    n_effector: int = 10,
    contigs: str = "A1-5/10-10/A16-20 B",
) -> Path:
    """An RFDiffusion stage output directory.

    Layout copied from tests/rfdiffusion/example_output_files/:

        <root>/rfdiffusion/design_<i>.pdb
        <root>/rfdiffusion/rfdiffusion_metrics.json
        <root>/rfdiffusion/filter_summary.json
        <root>/rfdiffusion/passing_designs.txt

    `designs` maps a design index to the keyword arguments for `complex_pdb`,
    so each design's geometry, and therefore each design's expected metrics,
    is stated at the call site.
    """
    designs = designs or {0: {}, 1: {"separation": 20.0}}
    d = root / "rfdiffusion"
    d.mkdir(parents=True, exist_ok=True)

    entries, passing = [], []
    for idx, kwargs in sorted(designs.items()):
        kw = {"n_receptor": n_receptor, "n_effector": n_effector, **kwargs}
        complex_pdb(d / f"design_{idx}.pdb", **kw)
        # The position order is the receptor's residue numbers in file order,
        # which is what every positional index downstream is counted against.
        start = kw.get("receptor_start", 1)
        order = list(range(start, start + kw["n_receptor"]))
        entries.append({
            "design": f"design_{idx}.pdb",
            "receptor_position_order": order,
            "per_design_design_residues": order[5:15],
            "receptor_contact_residues": order[5:10],
            "n_contact_pairs": kwargs.get("n_contact_pairs", 10),
        })
        if kwargs.get("passes", True):
            passing.append(f"design_{idx}")

    (d / "rfdiffusion_metrics.json").write_text(
        json.dumps({"contigs": contigs, "designs": entries}, indent=2))
    (d / "filter_summary.json").write_text(json.dumps({
        "n_designs": len(entries),
        "n_passing": len(passing),
        "contigs": contigs,
    }, indent=2))
    (d / "passing_designs.txt").write_text("\n".join(passing) + "\n")
    return root


# ══════════════════════════════════════════════════════════════════════
# The negative-steering engine's per-sequence workdir
# ══════════════════════════════════════════════════════════════════════
#
# Layout copied from
# tests/negative_steering/example_output_files/negative_steering/runs/
#     design_42_seq_0/cycle_0/...
# and from the layout `NegativeSteeringRun.from_workdir` documents:
#
#     <workdir>/row_type.txt
#     <workdir>/run_one_runtime_sec.txt
#     <workdir>/cycle_<C>/initial[_sN]/                    cold-start seeds
#     <workdir>/cycle_<C>/steered/design_NN_sS/            steered
#     <workdir>/cycle_<C>/reversions/rev_design_NN_sS_sR/  reversion

# The fields a confidence sidecar carries in the committed reference set. The
# values are arbitrary but well-formed, so a test that cares about a number
# passes its own payload instead.
BOLTZ_CONFIDENCE_FIELDS = {
    "confidence_score": 0.6,
    "ptm": 0.5,
    "iptm": 0.4,
    "complex_plddt": 0.7,
}


def boltz_prediction_dir(
    parent: Path,
    *,
    model_sample: int = 0,
    write_pdb: bool = True,
    write_confidence: bool = True,
    flat_pdb: bool = False,
    confidence: Mapping[str, object] | None = None,
    mutations: Sequence[tuple[int, str, str]] | str | None = None,
    n_receptor: int = 6,
    n_effector: int = 4,
) -> Path:
    """One Boltz prediction directory, as the engine writes it.

        <parent>/boltz_results_input/predictions/input/pdb/input_model_<M>.pdb
        <parent>/boltz_results_input/predictions/input/confidence_input_model_<M>.json
        <parent>/mutations.tsv

    `_find_boltz_outputs` needs the PDB and its confidence sidecar together
    and reports None for either that is absent, so `write_pdb=False` and
    `write_confidence=False` build the two half-written directories a killed
    GPU job leaves behind.

    `flat_pdb=True` writes the PDB beside the sidecar instead of under `pdb/`,
    which is the older engine layout `_find_boltz_outputs` falls back to.

    `mutations` writes mutations.tsv, either as (pos1, wt, mut) triples or as
    a string written verbatim for the malformed-file cases. pos1 is 1-based on
    the designed-frame receptor. It is NOT the engine plan.json's
    `mutated_positions`, which are 0-based.
    """
    pred_root = parent / "boltz_results_input" / "predictions" / "input"
    pred_root.mkdir(parents=True, exist_ok=True)

    if write_pdb:
        name = f"input_model_{model_sample}.pdb"
        dest = pred_root / name if flat_pdb else pred_root / "pdb" / name
        complex_pdb(dest, n_receptor=n_receptor, n_effector=n_effector)

    if write_confidence:
        payload = dict(BOLTZ_CONFIDENCE_FIELDS if confidence is None
                       else confidence)
        (pred_root / f"confidence_input_model_{model_sample}.json").write_text(
            json.dumps(payload))

    if mutations is not None:
        if isinstance(mutations, str):
            text = mutations
        else:
            rows = "".join(f"{pos}\t{wt}\t{mut}\n" for pos, wt, mut in mutations)
            text = "pos1\twt\tmut\n" + rows
        (parent / "mutations.tsv").write_text(text)

    return parent


def negsteer_workdir(
    root: Path,
    *,
    cycle: int = 0,
    cold_start: Sequence[int] | Mapping[int, Mapping] = (0, 1, 2),
    steered: Mapping[str, Mapping] | None = None,
    reversions: Mapping[str, Mapping] | None = None,
    row_type: str | None = "steered",
    runtime_sec: object | None = None,
) -> Path:
    """One MPNN sequence's finished negsteer workdir.

    `cold_start` names the seed indices to write. Seed 0 goes in `initial/`
    and seed N in `initial_sN/`, which is the engine's naming. Pass a mapping
    instead to give a seed its own `boltz_prediction_dir` keywords.

    `steered` and `reversions` are keyed by the subdirectory name, so the
    engine's `design_NN_sS` and `rev_design_NN_sS_sR` conventions stay visible
    at the call site rather than being assembled here. Each value is the
    keywords for that directory's `boltz_prediction_dir`.

    `row_type=None` omits row_type.txt, which is the sidecar-absent case.
    `runtime_sec` is written verbatim, so a test can pass a string that does
    not parse as a float.
    """
    cycle_dir = root / f"cycle_{cycle}"
    cycle_dir.mkdir(parents=True, exist_ok=True)

    seeds = (cold_start if isinstance(cold_start, Mapping)
             else {s: {} for s in cold_start})
    for seed_idx, kwargs in seeds.items():
        sub = "initial" if seed_idx == 0 else f"initial_s{seed_idx}"
        boltz_prediction_dir(cycle_dir / sub, **dict(kwargs))

    for group, entries in (("steered", steered), ("reversions", reversions)):
        if entries is None:
            continue
        (cycle_dir / group).mkdir(parents=True, exist_ok=True)
        for name, kwargs in entries.items():
            boltz_prediction_dir(cycle_dir / group / name, **dict(kwargs))

    if row_type is not None:
        (root / "row_type.txt").write_text(row_type)
    if runtime_sec is not None:
        (root / "run_one_runtime_sec.txt").write_text(str(runtime_sec))
    return root
