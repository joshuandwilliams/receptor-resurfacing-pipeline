"""Per-residue secondary structure for the design backbones, via ChimeraX DSSP.

The RFdiffusion outputs carry N, CA, C and O, which is everything DSSP needs to
find backbone hydrogen bonds, so the strand assignment is computed rather than
inferred from geometry by hand. mkdssp is not installed locally; ChimeraX has
the same algorithm built in and is already used elsewhere in this chapter for
rendering, so it is one fewer dependency.

Writes design, residue number, and whether the residue is strand, helix or coil.

Usage:
    ChimeraX --exit --nogui --script dssp_extract.py
"""

import csv
import glob
import os

from chimerax.core.commands import run

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
OUT = os.path.join(HERE, "design_secondary_structure.csv")

# The wild type goes through the same DSSP call as the designs, so its strand
# assignment is made the same way rather than quoted from the literature.
paths = sorted(glob.glob(os.path.join(DATA, "rfdiffusion", "design_*.pdb")))
paths.append(os.path.join(DATA, "wt_reference.pdb"))

rows = []
for path in paths:
    design = os.path.basename(path).replace(".pdb", "")
    if design == "wt_reference":
        design = "wild-type"
    structures = run(session, f"open {path}")  # noqa: F821
    st = structures[0]
    run(session, "dssp")  # noqa: F821
    for r in st.residues:
        if r.chain_id != "A":
            continue
        if r.is_strand:
            ss = "strand"
        elif r.is_helix:
            ss = "helix"
        else:
            ss = "coil"
        rows.append({"design": design, "resnum": r.number, "ss": ss})
    run(session, "close")  # noqa: F821

with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["design", "resnum", "ss"])
    w.writeheader()
    w.writerows(rows)

print(f"wrote {OUT} ({len(rows)} residues, "
      f"{len(set(r['design'] for r in rows))} designs)")
