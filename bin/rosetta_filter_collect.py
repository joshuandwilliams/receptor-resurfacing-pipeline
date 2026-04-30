#!/usr/bin/env python3
"""
rosetta_filter_collect.py
-------------------------
Collect Rosetta InterfaceAnalyzer score files from parallel ROSETTA_SC runs,
apply shape complementarity (Sc) threshold, and produce:

    rosetta_filter_metrics.json    Per-design metrics (all designs)
    rosetta_passing_designs.txt    PDB filenames passing the Sc filter
    rosetta_filter_summary.json    Overall filter statistics
    passing/                       Copies of PDBs that pass the filter

Score files are expected to be named interface_scores_<design_stem>.sc
and located in --score-dir.  Corresponding PDBs are <design_stem>.pdb
in --pdb-dir.

Metrics extracted from InterfaceAnalyzer:
    sc_value            Shape complementarity (0–1; >0.65 is good)
    dG_separated        Binding energy in REU (more negative = stronger).
                        NOTE: at this pre-MPNN stage the input PDBs are
                        polyvaline backbones with no real sidechains, so
                        dG_separated, packstat, and delta_unsatHbonds are
                        dominated by valine-clash artefacts and should
                        NOT be interpreted as physically meaningful
                        binding energies.  Only sc_value (purely geometric)
                        is used as a filter gate at this stage.  The
                        energy terms become meaningful after ProteinMPNN
                        sequence design and downstream relaxation.
    dSASA_int           Buried surface area at the interface (Å²)
    dG_separated/dSASAx100  Binding energy density (same caveat as above)
    packstat            Interface packing quality (0–1, see caveat above)
    delta_unsatHbonds   Buried unsatisfied H-bonds (see caveat above)
    nres_int            Number of interface residues
"""

import argparse
import glob
import json
import math
import os
import re
import shutil
import sys


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--score-dir", default=".",
                        help="Directory containing interface_scores_*.sc files")
    parser.add_argument("--pdb-dir", default=".",
                        help="Directory containing the corresponding design PDBs")
    parser.add_argument("--sc-threshold", type=float, default=0.62,
                        help="Minimum sc_value to pass (default: 0.62)")
    return parser.parse_args()


def parse_rosetta_scorefile(path):
    """
    Parse a Rosetta score file into a list of dicts.

    Rosetta score files have the format:
        SEQUENCE: ...
        SCORE: <header fields...>  description
        SCORE: <data fields...>    <pdb_name>
    """
    results = []
    header = None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("SEQUENCE"):
                continue
            if line.startswith("SCORE:") and header is None:
                header = line.split()[1:]
                continue
            if line.startswith("SCORE:") and header is not None:
                fields = line.split()[1:]
                # Rosetta sometimes has fewer fields than headers; zip handles it
                data = dict(zip(header, fields))
                results.append(data)
    return header, results


def safe_float(val, default=float("nan")):
    """Convert string to float, returning default on failure."""
    try:
        return float(val)
    except (ValueError, TypeError):
        return default


def main():
    args = parse_args()

    # ── Discover score files ──────────────────────────────────────────────
    score_files = sorted(glob.glob(
        os.path.join(args.score_dir, "interface_scores_*.sc")
    ))

    if not score_files:
        print("WARNING: No InterfaceAnalyzer score files found")
        # Write empty outputs so Nextflow doesn't fail
        with open("rosetta_filter_metrics.json", "w") as f:
            json.dump({"designs": [], "sc_threshold": args.sc_threshold}, f)
        with open("rosetta_passing_designs.txt", "w") as f:
            pass
        with open("rosetta_filter_summary.json", "w") as f:
            json.dump({"n_total": 0, "n_passing": 0, "n_filtered": 0,
                        "sc_threshold": args.sc_threshold}, f, indent=2)
        os.makedirs("passing", exist_ok=True)
        sys.exit(0)

    print(f"Found {len(score_files)} score files")

    # ── Key metrics to extract ────────────────────────────────────────────
    metric_keys = [
        "sc_value",
        "dG_separated",
        "dSASA_int",
        "dG_separated/dSASAx100",
        "packstat",
        "delta_unsatHbonds",
        "nres_int",
    ]

    # ── Parse all score files ─────────────────────────────────────────────
    all_designs = []
    passing = []
    filtered = []

    for sc_file in score_files:
        # Extract design stem from filename: interface_scores_<stem>.sc
        basename = os.path.basename(sc_file)
        match = re.match(r"interface_scores_(.+)\.sc$", basename)
        if not match:
            print(f"  WARNING: Could not parse design name from {basename}, skipping")
            continue
        design_stem = match.group(1)
        pdb_filename = f"{design_stem}.pdb"

        header, results = parse_rosetta_scorefile(sc_file)

        if not results:
            print(f"  WARNING: No scores in {basename}")
            continue

        # InterfaceAnalyzer with -out:no_nstruct_label produces one row
        row = results[0]

        # Extract metrics
        metrics = {}
        for key in metric_keys:
            metrics[key] = safe_float(row.get(key))

        sc_val = metrics["sc_value"]
        sc_is_nan = math.isnan(sc_val)
        passes = (not sc_is_nan) and sc_val >= args.sc_threshold

        def _round_or_none(v, ndigits):
            return None if math.isnan(v) else round(v, ndigits)

        def _int_or_none(v):
            return None if math.isnan(v) else int(v)

        design_entry = {
            "design": pdb_filename,
            "design_stem": design_stem,
            "sc_value": _round_or_none(sc_val, 4),
            "dG_separated": _round_or_none(metrics["dG_separated"], 3),
            "dSASA_int": _round_or_none(metrics["dSASA_int"], 1),
            "dG_dSASA_density": _round_or_none(metrics["dG_separated/dSASAx100"], 4),
            "packstat": _round_or_none(metrics["packstat"], 4),
            "delta_unsatHbonds": _round_or_none(metrics["delta_unsatHbonds"], 1),
            "nres_int": _int_or_none(metrics["nres_int"]),
            "passes_filter": passes,
        }
        all_designs.append(design_entry)

        status = "PASS" if passes else "FAIL"
        sc_str = f"{sc_val:.3f}" if not sc_is_nan else "N/A"
        dg_str = (f"{metrics['dG_separated']:.1f}"
                  if not math.isnan(metrics["dG_separated"]) else "N/A")
        dsasa_str = (f"{metrics['dSASA_int']:.0f}"
                     if not math.isnan(metrics["dSASA_int"]) else "N/A")
        print(f"  {pdb_filename}: Sc={sc_str}  dG={dg_str}  "
              f"dSASA={dsasa_str}  [{status}]")

        (passing if passes else filtered).append(pdb_filename)

    # ── Copy passing PDBs ─────────────────────────────────────────────────
    os.makedirs("passing", exist_ok=True)
    for pdb_name in passing:
        src = os.path.join(args.pdb_dir, pdb_name)
        if os.path.exists(src):
            shutil.copy(src, os.path.join("passing", pdb_name))
        else:
            print(f"  WARNING: PDB {src} not found for passing design")

    # ── Write outputs ─────────────────────────────────────────────────────
    with open("rosetta_filter_metrics.json", "w") as f:
        json.dump({
            "designs": all_designs,
            "sc_threshold": args.sc_threshold,
        }, f, indent=2)

    with open("rosetta_passing_designs.txt", "w") as f:
        for name in passing:
            f.write(name + "\n")

    with open("rosetta_filter_summary.json", "w") as f:
        json.dump({
            "n_total": len(all_designs),
            "n_passing": len(passing),
            "n_filtered": len(filtered),
            "sc_threshold": args.sc_threshold,
            "filtered_designs": filtered,
        }, f, indent=2)

    print(f"\nRosetta filter summary: {len(passing)}/{len(all_designs)} designs pass "
          f"(Sc threshold: {args.sc_threshold})")

    if not passing:
        print("WARNING: No designs passed the Sc filter!", file=sys.stderr)


if __name__ == "__main__":
    main()