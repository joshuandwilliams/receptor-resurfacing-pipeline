#!/usr/bin/env python3
"""
run_rosetta_metrics.py
----------------------
The Rosetta arm of the orthogonal-metrics streams, run on one cohort
survivor.  Two stages:

  1. FastRelax (1 repeat, Cartesian, ref2015_cart) on the bound
     complex.  Removes clashes from the predicted (Boltz) structure
     and brings it onto Rosetta's energy surface.  Run through
     RosettaScripts with the pinned XML at bin/fastrelax_for_ia.xml.

  2. InterfaceAnalyzer on the relaxed PDB.  Computes Lawrence-Colman
     shape complementarity (Sc) and the separated ΔG (rosetta_ddg)
     with sidechain repacking on the bound and unbound states.

Relax is the standard pre-treatment for ΔΔG scoring of predicted
structures.  Without it, ΔΔG is inflated 10-20 REU by clashes and
idealised-bond-geometry artefacts that Rosetta penalises but a real
folded complex would not have.  1 repeat matches Bennett et al. 2023's
de novo binder calibration, so the −30 REU threshold
(pipeline_thresholds.orthogonal_ddg_max) is directly comparable.

Emits a one-row CSV with columns seq_name, sc, rosetta_ddg,
rosetta_failures.

Runtime: ~3-6 min per survivor.  Relax dominates, the IA call is ~15 s.

Assumes Rosetta.img has both `rosetta_scripts.{flavour}.linuxgccrelease`
and `InterfaceAnalyzer.{flavour}.linuxgccrelease` on PATH.
"""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import List, Optional, Tuple

# ── Binary discovery ─────────────────────────────────────────────────────
# Rosetta ships several build flavours.  Try each in priority order,
# mpi first because it is the most common in HPC singularity images.
ROSETTA_FLAVOURS = ("mpi", "default", "static", "linuxgccrelease")


def _find_rosetta_binary(stem: str) -> Optional[str]:
    """Locate a Rosetta binary by stem (e.g. 'InterfaceAnalyzer',
    'rosetta_scripts').  Returns the first existing variant in the
    flavour priority order, or None if none are on PATH."""
    candidates = [
        f"{stem}.{flavour}.linuxgccrelease"
        for flavour in ROSETTA_FLAVOURS
        if flavour != "linuxgccrelease"
    ]
    candidates.append(f"{stem}.linuxgccrelease")
    for b in candidates:
        path = shutil.which(b)
        if path:
            return path
    return None


# ── Score-file parser ────────────────────────────────────────────────────

def _parse_ia_scorefile(
    score_path: Path,
) -> Tuple[Optional[float], Optional[float], Optional[str]]:
    """Parse the score.sc file InterfaceAnalyzer produces.

    The file starts with SEQUENCE lines, then a SCORE: header line of
    column names, then SCORE: rows of values.  Two columns are read:

      sc_value      → Lawrence-Colman shape complementarity
      dG_separated  → ΔG on chain separation (our rosetta_ddg)

    Returns (sc, ddg, error_string).  error_string is None on full
    success, otherwise a comma-separated list of token-style failure
    reasons such as 'sc_column_missing,ddg_not_numeric'.
    """
    try:
        with open(score_path) as f:
            lines = [ln.strip() for ln in f if ln.strip()]
    except OSError as e:
        return None, None, f"scorefile_read_failed:{type(e).__name__}"

    header: Optional[List[str]] = None
    data_row: Optional[List[str]] = None
    for ln in lines:
        if not ln.startswith("SCORE:"):
            continue
        parts = ln.split()
        if header is None and "total_score" in parts:
            header = parts
        elif header is not None and parts[0] == "SCORE:":
            data_row = parts
            break

    if header is None or data_row is None:
        return None, None, "scorefile_no_data_row"

    def get(col_name: str) -> Optional[str]:
        try:
            idx = header.index(col_name)
            return data_row[idx]
        except (ValueError, IndexError):
            return None

    sc_str = get("sc_value")
    ddg_str = get("dG_separated")

    sc_f: Optional[float] = None
    ddg_f: Optional[float] = None
    err_tokens: List[str] = []
    if sc_str is not None:
        try:
            sc_f = float(sc_str)
        except ValueError:
            err_tokens.append("sc_not_numeric")
    else:
        err_tokens.append("sc_column_missing")
    if ddg_str is not None:
        try:
            ddg_f = float(ddg_str)
        except ValueError:
            err_tokens.append("ddg_not_numeric")
    else:
        err_tokens.append("ddg_column_missing")

    err = ",".join(err_tokens) if err_tokens else None
    return sc_f, ddg_f, err


# ── Stage 1: FastRelax ───────────────────────────────────────────────────

def _run_fast_relax(
    rosetta_scripts_bin: str,
    xml_path: Path,
    input_pdb: Path,
    workdir: Path,
) -> Tuple[Optional[Path], List[str]]:
    """Run RosettaScripts FastRelax on `input_pdb`, write the relaxed
    PDB into `workdir`, and return its path (or None on failure).

    The XML at `xml_path` defines one Cartesian FastRelax mover scored
    with ref2015_cart.  The repeat count is NOT in the XML — this build
    rejects `default_repeats` as an attribute, so it is passed below as
    the `-relax:default_repeats 1` command-line option.  1 repeat
    matches Bennett 2023's de novo binder calibration, keeping the ΔΔG
    values comparable to its −30 REU threshold.

    Failures (non-zero exit, no output PDB, timeout) come back as a
    list of token-style strings the caller appends to
    `rosetta_failures` in the emitted CSV.
    """
    failures: List[str] = []

    # rosetta_scripts appends a per-nstruct numeric suffix by default.
    # -out:no_nstruct_label suppresses that and -out:suffix adds a fixed
    # _0001, so the output filename is deterministic.
    out_pdb = workdir / f"{input_pdb.stem}_0001.pdb"
    cmd = [
        rosetta_scripts_bin,
        "-database", "/opt/rosetta/main/database",
        "-parser:protocol", str(xml_path),
        "-s", str(input_pdb),
        "-out:path:all", str(workdir),
        "-out:suffix", "_0001",
        "-out:no_nstruct_label",
        "-nstruct", "1",
        "-ignore_unrecognized_res",
        # FastRelax repeat count.  Rosetta 2025.37's schema rejects
        # `default_repeats` as an XML mover attribute, so the only
        # place to set it is here as a top-level option.
        "-relax:default_repeats", "1",
        "-mute", "all",
    ]

    try:
        # 1-repeat FastRelax on a ~150-residue complex takes 2-4 min.
        # Cold-start CPU contention on the cluster can push it past 10,
        # hence the generous 15-minute cap.
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=900,
        )
    except subprocess.TimeoutExpired:
        failures.append("relax_timed_out")
        return None, failures

    if result.returncode != 0:
        # Rosetta writes diagnostics to STDOUT, not stderr.  Even
        # -mute all leaves the [ ERROR ] block on stdout.  An earlier
        # version echoed stderr only, lost the whole error message, and
        # left ROSETTA_CRASH.log as the only clue.  Echo the tail of
        # both streams so failures are diagnosable from .command.err.
        for ln in (result.stdout or "").splitlines()[-25:]:
            print(f"[relax-stdout] {ln}", file=sys.stderr)
        for ln in (result.stderr or "").splitlines()[-10:]:
            print(f"[relax-stderr] {ln}", file=sys.stderr)
        failures.append(f"relax_nonzero_exit:{result.returncode}")
        return None, failures

    if not out_pdb.is_file():
        # rosetta_scripts sometimes writes to a slightly different
        # name depending on suffix handling; fall back to the first
        # .pdb in workdir that isn't the input itself.
        candidates = [
            p for p in workdir.glob("*.pdb")
            if p.resolve() != input_pdb.resolve()
        ]
        if candidates:
            out_pdb = candidates[0]
        else:
            failures.append("relax_no_output_pdb")
            return None, failures

    return out_pdb, failures


# ── Stage 2: InterfaceAnalyzer ───────────────────────────────────────────

def _run_interface_analyzer(
    ia_bin: str,
    relaxed_pdb: Path,
    workdir: Path,
) -> Tuple[Optional[float], Optional[float], List[str]]:
    """Run InterfaceAnalyzer on the relaxed PDB and return
    (sc, rosetta_ddg, failures).

    The relaxed input is already on Rosetta's energy surface, so the
    bound-state score is meaningful.  `pack_separated` and
    `pack_input` pack the unbound and bound sidechains before scoring.
    Both are defensive.  They should be no-ops on a relaxed PDB, but
    they stop the numbers drifting if an input ever bypasses relax.
    """
    failures: List[str] = []
    score_path = workdir / "score.sc"

    # InterfaceAnalyzer in this Rosetta build (2025.37) does NOT
    # recognise -compute_separated_dG or -compute_packstat as
    # top-level flags — those metrics are computed by default.
    # -compute_interface_sc IS still a valid flag.  -database is
    # required to avoid an init-time crash looking for default
    # database paths that don't exist in this container.
    cmd = [
        ia_bin,
        "-database", "/opt/rosetta/main/database",
        "-s", str(relaxed_pdb),
        "-no_optH", "false",
        "-ignore_unrecognized_res",
        "-pack_separated",
        "-pack_input",
        "-compute_interface_sc", "true",
        "-add_regular_scores_to_scorefile",
        "-use_input_sc",
        "-out:file:score_only", str(score_path),
        "-out:no_nstruct_label",
        "-mute", "all",
    ]

    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=300,
        )
    except subprocess.TimeoutExpired:
        failures.append("ia_timed_out")
        return None, None, failures

    if result.returncode != 0:
        # Same rationale as _run_fast_relax: Rosetta logs to stdout,
        # not stderr.  Echo both so future failures are diagnosable
        # from .command.err alone.
        for ln in (result.stdout or "").splitlines()[-25:]:
            print(f"[ia-stdout] {ln}", file=sys.stderr)
        for ln in (result.stderr or "").splitlines()[-10:]:
            print(f"[ia-stderr] {ln}", file=sys.stderr)
        failures.append(f"ia_nonzero_exit:{result.returncode}")

    if not score_path.is_file():
        failures.append("ia_no_scorefile_produced")
        return None, None, failures

    sc, ddg, parse_err = _parse_ia_scorefile(score_path)
    if parse_err:
        failures.append(parse_err)
    return sc, ddg, failures


# ── Top-level orchestration ──────────────────────────────────────────────

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seq-name", required=True)
    ap.add_argument("--canonical-pdb", required=True, type=Path)
    ap.add_argument("--receptor-chain", default="A")
    ap.add_argument("--effector-chain", default="B")
    ap.add_argument(
        "--fast-relax-xml", type=Path, default=None,
        help="Path to the FastRelax RosettaScripts XML.  Defaults to "
             "<this script's dir>/fastrelax_for_ia.xml.",
    )
    ap.add_argument("--output-csv", required=True, type=Path)
    args = ap.parse_args()

    if args.fast_relax_xml is None:
        args.fast_relax_xml = (
            Path(__file__).resolve().parent / "fastrelax_for_ia.xml"
        )
    if not args.fast_relax_xml.is_file():
        print(
            f"ERROR: FastRelax XML not found at {args.fast_relax_xml}",
            file=sys.stderr,
        )
        # Still emit a row so the merge step doesn't lose this survivor;
        # all metrics blank, single failure tag.
        _emit_failure_row(
            args.output_csv, args.seq_name,
            ["fastrelax_xml_not_found"],
        )
        return 1

    failures: List[str] = []
    sc_val: Optional[float] = None
    ddg_val: Optional[float] = None

    rs_bin = _find_rosetta_binary("rosetta_scripts")
    ia_bin = _find_rosetta_binary("InterfaceAnalyzer")

    if rs_bin is None:
        failures.append("rosetta_scripts_binary_not_found")
    if ia_bin is None:
        failures.append("interface_analyzer_binary_not_found")

    if rs_bin is not None and ia_bin is not None:
        with tempfile.TemporaryDirectory(prefix="rosetta_run_") as tmpdir:
            tmp = Path(tmpdir)
            relax_dir = tmp / "relax"
            ia_dir = tmp / "ia"
            relax_dir.mkdir()
            ia_dir.mkdir()

            relaxed_pdb, relax_failures = _run_fast_relax(
                rs_bin, args.fast_relax_xml, args.canonical_pdb, relax_dir,
            )
            failures.extend(relax_failures)

            if relaxed_pdb is not None:
                sc_val, ddg_val, ia_failures = _run_interface_analyzer(
                    ia_bin, relaxed_pdb, ia_dir,
                )
                failures.extend(ia_failures)

    _emit_row(
        args.output_csv, args.seq_name,
        sc_val, ddg_val, failures,
    )
    return 0


def _emit_row(
    output_csv: Path,
    seq_name: str,
    sc: Optional[float],
    ddg: Optional[float],
    failures: List[str],
) -> None:
    row = {
        "seq_name": seq_name,
        "sc": f"{sc:.4f}" if sc is not None else "",
        "rosetta_ddg": f"{ddg:.3f}" if ddg is not None else "",
        "rosetta_failures": ",".join(failures),
    }
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(output_csv, "w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=["seq_name", "sc", "rosetta_ddg", "rosetta_failures"],
            extrasaction="ignore",
        )
        w.writeheader()
        w.writerow(row)
    print(
        f"[{seq_name}] sc={row['sc']}  ddg={row['rosetta_ddg']}  "
        f"flags={row['rosetta_failures']}"
    )


def _emit_failure_row(
    output_csv: Path, seq_name: str, failures: List[str],
) -> None:
    """Convenience: emit an all-blank row with only failure tags
    populated.  Used when we can't even start the run (missing XML,
    missing binaries) so the downstream merge step still has a row."""
    _emit_row(output_csv, seq_name, None, None, failures)


if __name__ == "__main__":
    raise SystemExit(main())
