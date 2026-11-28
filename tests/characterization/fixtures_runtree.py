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
