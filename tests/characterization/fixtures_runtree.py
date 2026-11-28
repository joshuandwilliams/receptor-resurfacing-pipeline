"""Build run trees and structures whose answers are known before the code runs.

Testing against a real pipeline output only proves the code still agrees with
its own last output. These builders state the answer first, then construct the
input that must produce it, so a wrong result is a failure rather than a new
blessed reference.

Two kinds of builder:

`complex_pdb` places a receptor and an effector at a geometry you specify, so
a contact count, a chain split and a centre-of-mass displacement are all known
in advance.

`rfdiffusion_tree` writes the directory layout an RFDiffusion stage produces,
with the metrics JSON and the design PDBs agreeing with each other, so the
readers and the filters can be driven end to end without a cluster.

`af3_output_tree` writes the directory layout an AlphaFold3-no-MSA run
produces, with a ground truth beside it, so each prediction's receptor-aligned
effector RMSD is a number stated at the call site.

Layouts are taken from the committed reference sets under
`tests/<module>/example_output_files/`, not from docstrings.
"""
from __future__ import annotations

import json
from pathlib import Path

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


def complex_cif(path: Path, **kwargs) -> Path:
    """`complex_pdb`'s geometry written as a real mmCIF.

    The AF3 walk globs for `.cif`, and gemmi chooses its parser from the
    extension, so a `.cif` file holding PDB text raises instead of loading.
    Going through gemmi rather than renaming the file keeps the fixture honest
    about that.

    gemmi is imported here rather than at module scope so the PDB builders
    stay usable in an environment without it.
    """
    import gemmi

    scratch = path.parent / f".{path.stem}.pdb"
    complex_pdb(scratch, **kwargs)
    st = gemmi.read_pdb_string(scratch.read_text())
    scratch.unlink()
    st.setup_entities()
    st.make_mmcif_document().write_file(str(path))
    return path


def af3_output_tree(
    path: Path,
    *,
    predictions,
    truth_separation: float = 4.0,
    n_receptor: int = 8,
    n_effector: int = 5,
    seed: int = 1,
) -> tuple[Path, Path]:
    """The AF3-no-MSA layout `af3_confidence` walks, with known ra_eff values.

    Returns `(ground_truth_pdb, output_dir)`:

        <path>/ground_truth.pdb
        <path>/output/seed-<seed>_sample-<m>/model.cif
        <path>/output/seed-<seed>_sample-<m>/summary_confidences.json

    `predictions` is a sequence of `(separation, iptm)` pairs, one per
    seed x diffusion sample. The ground truth is `complex_pdb` at
    `truth_separation` and each prediction is the same complex with the
    effector moved to its own separation along z, the receptor untouched. A
    rigid displacement of the effector alone means the receptor alignment is
    the identity, so that prediction's receptor-aligned effector RMSD is
    exactly `abs(separation - truth_separation)` angstroms.

    `iptm=None` omits that prediction's summary_confidences.json, which is how
    a prediction whose confidence write never landed looks on disk.
    """
    path.mkdir(parents=True, exist_ok=True)
    truth = complex_pdb(path / "ground_truth.pdb", n_receptor=n_receptor,
                        n_effector=n_effector, separation=truth_separation)
    out = path / "output"
    for sample, (separation, iptm) in enumerate(predictions):
        d = out / f"seed-{seed}_sample-{sample}"
        complex_cif(d / "model.cif", n_receptor=n_receptor,
                    n_effector=n_effector, separation=separation)
        if iptm is not None:
            (d / "summary_confidences.json").write_text(json.dumps({"iptm": iptm}))
    return truth, out


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
