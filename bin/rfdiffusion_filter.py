#!/usr/bin/env python3
"""
rfdiffusion_filter.py
---------------------
Analyse and filter RFDiffusion designs by interface contact overlap with
the design (de novo) region.

For each design PDB:
  - Split single-chain output into receptor/effector at the Cα chain break
  - Map receptor residues to original contig numbering
  - Find receptor–effector Cα contacts within a distance cutoff
  - Compute design-region RMSD relative to the input PDB
  - Compute fraction of contacts inside vs outside the design region

Outputs:
    rfdiffusion_metrics.json  - per-design metrics for plotting
    passing_designs.txt       - PDB filenames that pass the filter
    filter_summary.json       - overall filter statistics
    split/design_*.pdb        - two-chain PDBs (receptor=A, effector=B)

NUMBERING CONVENTIONS — TWO DIFFERENT SPACES IN PLAY:
  - rfdiffusion_metrics.json uses "contig-space" residue numbers — these
    match the input PDB's original numbering (e.g. B187-262 receptor,
    C33-113 effector).  De novo positions get synthetic numbers
    >= max_fixed_resnum + 1000 to avoid collisions.
  - split/design_*.pdb uses "1-based output space" — every chain is
    renumbered from 1, regardless of the input PDB.  This matches what
    RFDiffusion's run_inference.py emits.
  Downstream code that joins metrics-JSON values to split-PDB residues
  must translate between the two spaces.
"""

import argparse
import glob
import json
import math
import os
import sys

import numpy as np

# Local helper for centroid/COM math.
_BIN_DIR = os.path.dirname(os.path.abspath(__file__))
if _BIN_DIR not in sys.path:
    sys.path.insert(0, _BIN_DIR)
import structure_metrics  # noqa: E402

from contig_utils import (
    resolve_contigs, get_expected_chain_lengths, parse_design_region,
)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-pdb", required=True, help="Input PDB for RFDiffusion")
    parser.add_argument("--design-dir", default=".", help="Directory with design_*.pdb")
    parser.add_argument("--contigs", required=True, help="Contig string")
    parser.add_argument("--hotspot", default="", help="Hotspot residues (e.g. 'B24,B25')")
    parser.add_argument("--receptor-chain", default="A")
    parser.add_argument("--effector-chain", default="B")
    parser.add_argument("--rfdiff-contact-cutoff", type=float, default=8.0,
                        help="Cα–Cα distance cutoff in Å (default: 8.0)")
    parser.add_argument("--min-hotspot-frac", type=float, default=0.0,
                        help="Min fraction of contacts in design region (default: 0.0)")
    return parser.parse_args()


# ═══════════════════════════════════════════════════════════════════════════
# PDB parsing
# ═══════════════════════════════════════════════════════════════════════════

def read_ca_atoms(pdb_path, chain=None):
    """Extract Cα atoms as list of dicts: [{chain, resnum, x, y, z}, ...].

    See also ``boltz2_negative_steering.read_ca_atoms`` — same name,
    different shape: that one returns ``List[CAEntry]`` with numpy xyz
    plus per-chain seq_index and is altloc/icode aware.  This one is
    a plain coord-dict reader optimised for RFDiffusion output PDBs
    (clean, no altlocs, no insertion codes).
    """
    atoms = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                ch = line[21]
                if chain is not None and ch != chain:
                    continue
                atoms.append({
                    "chain": ch,
                    "resnum": int(line[22:26].strip()),
                    "x": float(line[30:38]),
                    "y": float(line[38:46]),
                    "z": float(line[46:54]),
                })
    return atoms


def ca_coords_array(atoms):
    """Convert atom dicts to (N, 3) numpy array."""
    if not atoms:
        return np.zeros((0, 3))
    return np.array([[a["x"], a["y"], a["z"]] for a in atoms])


def count_chains(pdb_path):
    """Return sorted list of unique chain IDs."""
    chains = set()
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM"):
                chains.add(line[21])
    return sorted(chains)


# ═══════════════════════════════════════════════════════════════════════════
# Chain break detection
# ═══════════════════════════════════════════════════════════════════════════

def _find_break_index(coords_or_atoms):
    """
    Find the index of the largest Cα–Cα gap in sequential residues.
    Accepts either a list of atom dicts or a list of (x,y,z) tuples.
    Returns (break_idx, break_dist).
    """
    n = len(coords_or_atoms)
    if n < 2:
        return 0, 0.0

    def _xyz(item):
        if isinstance(item, dict):
            return item["x"], item["y"], item["z"]
        return item  # tuple

    max_dist = 0.0
    max_idx = 0
    for i in range(n - 1):
        x1, y1, z1 = _xyz(coords_or_atoms[i])
        x2, y2, z2 = _xyz(coords_or_atoms[i + 1])
        dist = math.sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)
        if dist > max_dist:
            max_dist = dist
            max_idx = i
    return max_idx, max_dist


def split_at_chain_break(ca_atoms):
    """Split Cα atoms at the chain break. Returns (seg1, seg2, break_dist)."""
    idx, dist = _find_break_index(ca_atoms)
    return ca_atoms[:idx + 1], ca_atoms[idx + 1:], dist


def identify_segments(seg1_atoms, seg2_atoms, expected_rec_len, expected_eff_len):
    """
    Determine which segment is receptor vs effector by matching to expected
    lengths. Effector length is always fixed, so match against that.
    Returns (rec_atoms, eff_atoms, rec_is_seg1).
    """
    s1, s2 = len(seg1_atoms), len(seg2_atoms)

    # Exact match
    if s1 == expected_rec_len and s2 == expected_eff_len:
        return seg1_atoms, seg2_atoms, True
    if s1 == expected_eff_len and s2 == expected_rec_len:
        return seg2_atoms, seg1_atoms, False

    # Fuzzy match by effector length
    if abs(s1 - expected_eff_len) < abs(s2 - expected_eff_len):
        print(f"  Fuzzy match: seg1 ({s1}) ≈ effector ({expected_eff_len}), "
              f"seg2 ({s2}) = receptor")
        return seg2_atoms, seg1_atoms, False
    else:
        print(f"  Fuzzy match: seg1 ({s1}) = receptor, "
              f"seg2 ({s2}) ≈ effector ({expected_eff_len})")
        return seg1_atoms, seg2_atoms, True


def write_split_pdb(input_pdb_path, output_path, expected_rec_len, expected_eff_len):
    """
    Split a single-chain RFDiffusion PDB at the chain break, assign chain A
    to receptor and B to effector, renumber from 1, write two-chain PDB.
    """
    with open(input_pdb_path) as f:
        all_lines = f.readlines()

    # Collect residue order and per-residue lines in a single pass
    residue_order = []
    residue_set = set()
    ca_coords = {}
    residue_lines = {}

    for line in all_lines:
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        resnum = int(line[22:26].strip())
        if resnum not in residue_set:
            residue_set.add(resnum)
            residue_order.append(resnum)
            residue_lines[resnum] = []
        residue_lines[resnum].append(line)
        if line[12:16].strip() == "CA":
            ca_coords[resnum] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))

    # Find chain break using shared helper
    ordered_ca = [ca_coords[r] for r in residue_order if r in ca_coords]
    # Map back: we need the break in residue_order space, not ca space
    # Since every residue should have a CA, these lists align
    ca_residues = [r for r in residue_order if r in ca_coords]
    break_idx, _ = _find_break_index(ordered_ca)

    seg1_residues = ca_residues[:break_idx + 1]
    seg2_residues = ca_residues[break_idx + 1:]

    # Assign receptor vs effector by effector length
    if abs(len(seg1_residues) - expected_eff_len) < abs(len(seg2_residues) - expected_eff_len):
        eff_residues, rec_residues = seg1_residues, seg2_residues
    else:
        rec_residues, eff_residues = seg1_residues, seg2_residues

    # Write two-chain PDB
    with open(output_path, "w") as f:
        for new_resnum, resnum in enumerate(rec_residues, start=1):
            for line in residue_lines.get(resnum, []):
                f.write(line[:21] + "A" + f"{new_resnum:4d}" + line[26:])
        f.write("TER\n")
        for new_resnum, resnum in enumerate(eff_residues, start=1):
            for line in residue_lines.get(resnum, []):
                f.write(line[:21] + "B" + f"{new_resnum:4d}" + line[26:])
        f.write("TER\nEND\n")


# ═══════════════════════════════════════════════════════════════════════════
# Residue mapping
# ═══════════════════════════════════════════════════════════════════════════

def _parse_receptor_segments(contigs, rec_chain):
    """Parse the receptor block from contigs into segment descriptors.

    Returns list of ('fixed', pdb_start, pdb_end) or ('denovo', min_len, max_len).
    """
    blocks = contigs.split()
    rec_block = None
    for block in blocks:
        # Match by checking whether any segment in this block starts with
        # the receptor chain letter, rather than substring matching the
        # whole block (which gives false positives for multi-letter chain
        # IDs or chain letters that happen to appear elsewhere).
        for seg in block.split("/"):
            seg = seg.strip()
            if seg and seg[0].upper() == rec_chain.upper():
                rec_block = block
                break
        if rec_block is not None:
            break
    if rec_block is None:
        return []

    seg_descs = []
    for seg in rec_block.split("/"):
        seg = seg.strip()
        if not seg or seg == "0":
            continue
        if seg[0].isalpha() and seg[0].upper() == rec_chain.upper():
            rest = seg[1:]
            if "-" in rest:
                parts = rest.split("-")
                seg_descs.append(("fixed", int(parts[0]), int(parts[1])))
            else:
                v = int(rest)
                seg_descs.append(("fixed", v, v))
        elif seg[0].isdigit():
            if "-" in seg:
                parts = seg.split("-")
                seg_descs.append(("denovo", int(parts[0]), int(parts[1])))
            else:
                v = int(seg)
                seg_descs.append(("denovo", v, v))
    return seg_descs


def build_receptor_resnum_map(rec_atoms, contigs, rec_chain):
    """
    Map each receptor Cα atom (by index) to a residue number for
    contact classification (design vs fixed).

    The total receptor length is known from the design PDB.  The total
    fixed length is known from the contig.  The difference is the total
    de novo length.  We walk the contig segments:

      - Fixed segments: consume exactly (end - start + 1) atoms and
        assign original PDB residue numbers.
      - De novo segments: we know the total de novo atom budget
        (n_atoms - total_fixed).  For multiple de novo segments, the
        last one gets the remainder after earlier ones are allocated
        proportionally by max_len.  Atoms are assigned unique residue
        numbers that don't collide with any fixed residue, starting
        from a synthetic base above all fixed residue numbers.

    Returns (mapping, segment_info) where:
      mapping: {positional_index: residue_number}
      segment_info: dict with:
        - per_region_lengths: [int, ...] one per de novo segment
        - per_region_atom_ranges: [(start_idx, end_idx), ...] atom index
          ranges for each de novo region (inclusive start, exclusive end)
        - receptor_position_order: [resnum, ...] all receptor resnums in
          contig segment walk order (correct positional ordering)

    Design atoms get numbers >= synthetic_base (guaranteed not in any
    fixed segment).
    """
    seg_descs = _parse_receptor_segments(contigs, rec_chain)
    if not seg_descs:
        return {}, {"per_region_lengths": [], "per_region_atom_ranges": [],
                     "receptor_position_order": []}

    n_atoms = len(rec_atoms)
    total_fixed = sum(d[2] - d[1] + 1 for d in seg_descs if d[0] == "fixed")
    total_denovo = max(0, n_atoms - total_fixed)

    # Synthetic base for de novo residue numbers: above all fixed residues
    # so they never collide.  These numbers will be in the design_residues set.
    max_fixed_resnum = max(
        (d[2] for d in seg_descs if d[0] == "fixed"), default=0
    )
    synthetic_base = max_fixed_resnum + 1000  # large gap to avoid any overlap

    # Compute actual length for each de novo segment.
    #
    # Strategy: any segment with min_len == max_len has a known fixed
    # length and gets exactly that.  All remaining (variable-length)
    # segments share the leftover budget.  If there are multiple
    # variable-length segments, the leftover is distributed proportionally
    # by max_len, with min/max clamping per segment and the final segment
    # absorbing any rounding remainder.
    #
    # This is more correct than the previous purely-proportional approach
    # because it can't violate fixed-length constraints (e.g. "6-6" must
    # be exactly 6), and it respects the user's contig spec instead of
    # imposing a heuristic split.
    denovo_segs = [(i, d) for i, d in enumerate(seg_descs) if d[0] == "denovo"]
    denovo_actual_lens = {}

    if len(denovo_segs) == 1:
        denovo_actual_lens[denovo_segs[0][0]] = total_denovo
    elif len(denovo_segs) > 1:
        # Step 1: assign fixed-length segments their exact length.
        fixed_len_segs = [(i, d) for i, d in denovo_segs if d[1] == d[2]]
        var_len_segs = [(i, d) for i, d in denovo_segs if d[1] != d[2]]
        for seg_idx, desc in fixed_len_segs:
            denovo_actual_lens[seg_idx] = desc[1]
        used_by_fixed = sum(d[1] for _, d in fixed_len_segs)
        remaining_for_var = max(0, total_denovo - used_by_fixed)

        # Step 2: distribute the leftover to variable-length segments.
        if len(var_len_segs) == 0:
            pass  # all segments are fixed-length, nothing to do
        elif len(var_len_segs) == 1:
            seg_idx, desc = var_len_segs[0]
            # Clamp to [min_len, max_len] to honour the spec, but warn
            # if the actual receptor length doesn't permit it.
            alloc = max(desc[1], min(desc[2], remaining_for_var))
            if alloc != remaining_for_var:
                print(f"WARNING: variable denovo segment {seg_idx} got {alloc} "
                      f"residues after clamping to [{desc[1]}, {desc[2]}], "
                      f"but receptor length suggests {remaining_for_var}.  "
                      f"Per-region split may not match RFDiffusion's actual "
                      f"output for this design.", file=sys.stderr)
            denovo_actual_lens[seg_idx] = alloc
        else:
            # Multiple variable-length segments: distribute proportionally
            # by max_len.  Final segment absorbs remainder.
            total_max = sum(d[2] for _, d in var_len_segs)
            allocated = 0
            for k, (seg_idx, desc) in enumerate(var_len_segs):
                if k < len(var_len_segs) - 1:
                    if total_max > 0:
                        alloc = round(remaining_for_var * desc[2] / total_max)
                    else:
                        alloc = remaining_for_var // len(var_len_segs)
                    alloc = max(desc[1], min(desc[2], alloc))
                    alloc = min(alloc, remaining_for_var - allocated)
                    denovo_actual_lens[seg_idx] = alloc
                    allocated += alloc
                else:
                    denovo_actual_lens[seg_idx] = remaining_for_var - allocated

    # Walk segments, build mapping + per-region metadata
    mapping = {}
    atom_idx = 0
    synthetic_counter = synthetic_base
    receptor_position_order = []
    per_region_lengths = []
    per_region_atom_ranges = []

    for s_idx, desc in enumerate(seg_descs):
        if desc[0] == "fixed":
            for r in range(desc[1], desc[2] + 1):
                if atom_idx < n_atoms:
                    mapping[atom_idx] = r
                    receptor_position_order.append(r)
                    atom_idx += 1

        elif desc[0] == "denovo":
            actual_len = denovo_actual_lens.get(s_idx, 0)
            region_start_idx = atom_idx
            for k in range(actual_len):
                if atom_idx < n_atoms:
                    mapping[atom_idx] = synthetic_counter
                    receptor_position_order.append(synthetic_counter)
                    synthetic_counter += 1
                    atom_idx += 1
            per_region_lengths.append(atom_idx - region_start_idx)
            per_region_atom_ranges.append((region_start_idx, atom_idx))

    # Leftover atoms (appended as extra de novo)
    if atom_idx < n_atoms:
        leftover_start = atom_idx
        while atom_idx < n_atoms:
            mapping[atom_idx] = synthetic_counter
            receptor_position_order.append(synthetic_counter)
            synthetic_counter += 1
            atom_idx += 1
        # Add leftover to last region or create a new one
        if per_region_atom_ranges:
            last_start, _ = per_region_atom_ranges[-1]
            per_region_atom_ranges[-1] = (last_start, atom_idx)
            per_region_lengths[-1] = atom_idx - last_start
        else:
            per_region_lengths.append(atom_idx - leftover_start)
            per_region_atom_ranges.append((leftover_start, atom_idx))

    segment_info = {
        "per_region_lengths": per_region_lengths,
        "per_region_atom_ranges": per_region_atom_ranges,
        "receptor_position_order": receptor_position_order,
    }

    return mapping, segment_info


# ═══════════════════════════════════════════════════════════════════════════
# Interface analysis (NumPy-accelerated)
# ═══════════════════════════════════════════════════════════════════════════

def find_contacts(rec_atoms, eff_atoms, cutoff):
    """
    Find Cα–Cα contacts within cutoff using NumPy broadcasting.
    Returns list of (rec_idx, eff_idx, distance) tuples.
    """
    if not rec_atoms or not eff_atoms:
        return []

    rec_coords = ca_coords_array(rec_atoms)  # (R, 3)
    eff_coords = ca_coords_array(eff_atoms)  # (E, 3)

    # Broadcast distance matrix: (R, E)
    diff = rec_coords[:, np.newaxis, :] - eff_coords[np.newaxis, :, :]
    dist_matrix = np.sqrt((diff ** 2).sum(axis=2))

    # Extract contacts below cutoff
    rec_idx, eff_idx = np.where(dist_matrix < cutoff)
    contacts = [
        (int(r), int(e), float(dist_matrix[r, e]))
        for r, e in zip(rec_idx, eff_idx)
    ]
    return contacts


def calc_motif_and_region_metrics(
    input_rec_atoms, design_rec_atoms, resnum_map, fixed_residues,
    seg_descs, per_region_atom_ranges,
):
    """
    Compute motif and per-region structural change metrics for one design.

    Uses **resnum-keyed correspondence**: each design atom is matched to its
    input atom via resnum_map (which records, per design positional index,
    the input PDB residue number assigned by build_receptor_resnum_map).
    Positional indexing across input and design lists is unsafe — RFDiffusion
    can change atom counts in the design region, which shifts every fixed
    atom that follows the gap by an unpredictable offset.  Looking atoms up
    by resnum sidesteps that entirely.

    Metrics returned (all in Å):

      motif_rmsd: float
        RMSD over the fixed-segment (motif) Cα atoms after Kabsch
        superposition on those same atoms.  Sanity check on RFDiffusion's
        faithfulness to the contig — should be small (typically <1 Å) for
        a well-behaved design where the motif residues were genuinely held
        in place.

      per_region_endpoint_distance_input: list[float]
      per_region_endpoint_distance_design: list[float]
        For each design region, the Cα–Cα distance between the first and
        last residue OF THAT REGION, measured in the input PDB (over the
        input gap residues between the flanking fixed segments) and in the
        design (over RFDiffusion's built region) respectively.  Different
        lengths are fine — only the endpoints matter.  Comparing these
        tells you whether RFDiffusion built something that spans the same
        physical reach as the original gap, or pulled the geometry tighter
        / stretched it further.  The input-side value is NaN whenever the
        design region is one-sided (terminal extension at the N or C end
        of the receptor): you cannot define a span without two anchors,
        even if there are input residues to compare to for COM purposes.

      per_region_com_displacement: list[float]
        For each design region, the distance between the centroid of the
        design region's Cα atoms and the centroid of the corresponding
        input residues, measured in the scaffold-aligned frame (post-
        Kabsch).  "Corresponding input residues" depends on the contig
        anchoring of the design region:

          * Both anchors present (most common): all input residues in
            the inclusive gap between the flanking fixed segments.
          * Leading anchor only (C-terminal extension): the N input
            residues immediately after the preceding fixed segment,
            where N is the design region's actual produced atom count.
          * Trailing anchor only (N-terminal extension): the N input
            residues immediately before the following fixed segment.
          * Neither anchor (free-floating chain, very rare): NaN.

        Length-independent: works for fixed and variable-length contigs.
        A value of 0 means RFDiffusion built something occupying the same
        average position as the corresponding input residues; large
        values mean the new region's centre of mass is somewhere
        different in space.  NaN only when there are no corresponding
        input residues at all.

    Returns NaN motif_rmsd (and empty per-region lists) if there are
    no fixed residues to align on at all.
    """
    if not input_rec_atoms or not design_rec_atoms:
        return float("nan"), [], [], []

    # ── Build input resnum → coord lookup (one entry per Cα). ───────────
    input_by_resnum = {a["resnum"]: np.array([a["x"], a["y"], a["z"]])
                       for a in input_rec_atoms}

    design_coords_all = ca_coords_array(design_rec_atoms)
    n_design = len(design_coords_all)

    # ── Gather (design_idx, input_resnum) pairs for fixed residues ──────
    matched_design_idx = []
    matched_input_coords = []
    for i in range(n_design):
        rn = resnum_map.get(i)
        if rn is None or rn not in fixed_residues:
            continue
        if rn not in input_by_resnum:
            # Fixed residue from contig spec not found in input PDB —
            # contig and PDB are inconsistent.  Skip silently; the caller
            # already reports such mismatches elsewhere.
            continue
        matched_design_idx.append(i)
        matched_input_coords.append(input_by_resnum[rn])

    if not matched_design_idx:
        return float("nan"), [], [], []

    P = design_coords_all[matched_design_idx]   # mobile (design)
    Q = np.array(matched_input_coords)          # reference (input)

    # ── Kabsch on the fixed residues ────────────────────────────────────
    # Centroid helper from bin/structure_metrics.py.
    P_centroid = structure_metrics.centroid(P)
    Q_centroid = structure_metrics.centroid(Q)
    P_centred = P - P_centroid
    Q_centred = Q - Q_centroid

    H = P_centred.T @ Q_centred
    U, _S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T

    # Transform the entire design coord set into the input frame.
    design_aligned = (design_coords_all - P_centroid) @ R.T + Q_centroid

    # Motif RMSD: how well the fixed residues line up after Kabsch.
    motif_diffs = design_aligned[matched_design_idx] - Q
    motif_rmsd = float(np.sqrt(np.mean(np.sum(motif_diffs ** 2, axis=1))))

    # ── Per-region endpoint distances and COM displacement ──────────────
    # For each design region we identify two pieces of input information:
    #
    #   * The flanking fixed segment anchors (prev_end_resnum,
    #     next_start_resnum), either of which may be None when the design
    #     region is at the N or C terminus of the receptor.
    #
    #   * The set of input residues this design region "replaces" for
    #     the purpose of COM comparison.  See the per-region loop for
    #     the per-anchor-case rules.
    #
    # Endpoint-distance and COM-displacement use these differently:
    # endpoint distance needs both anchors (you can't span without two
    # ends); COM displacement only needs *some* corresponding input
    # residues, which is well-defined even with a one-sided anchor.

    denovo_anchors = []  # list of (prev_end, next_start) per denovo segment
    for s_idx, desc in enumerate(seg_descs):
        if desc[0] != "denovo":
            continue
        prev_end = None
        for j in range(s_idx - 1, -1, -1):
            if seg_descs[j][0] == "fixed":
                prev_end = seg_descs[j][2]
                break
        next_start = None
        for j in range(s_idx + 1, len(seg_descs)):
            if seg_descs[j][0] == "fixed":
                next_start = seg_descs[j][1]
                break
        denovo_anchors.append((prev_end, next_start))

    per_region_endpoint_distance_input = []
    per_region_endpoint_distance_design = []
    per_region_com_displacement = []

    for region_idx, (atom_start, atom_end) in enumerate(per_region_atom_ranges):
        # ── Design side ─────────────────────────────────────────────────
        if atom_end <= atom_start or atom_end > n_design:
            per_region_endpoint_distance_input.append(float("nan"))
            per_region_endpoint_distance_design.append(float("nan"))
            per_region_com_displacement.append(float("nan"))
            continue
        design_region_coords = design_aligned[atom_start:atom_end]
        n_region = len(design_region_coords)
        if n_region >= 2:
            ep_design = float(np.linalg.norm(
                design_region_coords[-1] - design_region_coords[0]))
        else:
            ep_design = float("nan")
        design_com = structure_metrics.centroid(design_region_coords)

        # ── Input side: pick the residues this region corresponds to ───
        if region_idx >= len(denovo_anchors):
            per_region_endpoint_distance_input.append(float("nan"))
            per_region_endpoint_distance_design.append(ep_design)
            per_region_com_displacement.append(float("nan"))
            continue
        prev_end, next_start = denovo_anchors[region_idx]

        # COM-comparison residues:
        #   both anchors:    inclusive gap between them
        #   leading-only:    N residues immediately after prev_end
        #   trailing-only:   N residues immediately before next_start
        #   neither:         no comparison
        # (N here is n_region — the design region's actual atom count, so
        # we always compare equally-sized atom sets when we can.)
        if prev_end is not None and next_start is not None:
            gap_resnums = list(range(prev_end + 1, next_start))
        elif prev_end is not None:
            gap_resnums = list(range(prev_end + 1, prev_end + 1 + n_region))
        elif next_start is not None:
            gap_resnums = list(range(next_start - n_region, next_start))
        else:
            gap_resnums = []

        input_gap_coords = [input_by_resnum[r] for r in gap_resnums
                            if r in input_by_resnum]

        # Endpoint distance (input): only meaningful when both anchors
        # exist AND the gap actually has ≥2 residues in the input PDB.
        # A leading- or trailing-only anchor can't define an endpoint
        # span, even if we picked some residues for the COM comparison.
        if (prev_end is not None and next_start is not None
                and len(input_gap_coords) >= 2):
            input_gap_arr = np.array(input_gap_coords)
            ep_input = float(np.linalg.norm(input_gap_arr[-1] - input_gap_arr[0]))
        else:
            ep_input = float("nan")

        # COM displacement: works as long as we found *any* corresponding
        # input residues at all.
        if input_gap_coords:
            input_com = structure_metrics.centroid(np.array(input_gap_coords))
            com_disp = float(np.linalg.norm(design_com - input_com))
        else:
            com_disp = float("nan")

        per_region_endpoint_distance_input.append(ep_input)
        per_region_endpoint_distance_design.append(ep_design)
        per_region_com_displacement.append(com_disp)

    return (motif_rmsd,
            per_region_endpoint_distance_input,
            per_region_endpoint_distance_design,
            per_region_com_displacement)


# ═══════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════

def main():
    args = parse_args()

    design_files = sorted(glob.glob(os.path.join(args.design_dir, "design_*.pdb")))
    if not design_files:
        print("WARNING: No design PDBs found")
        with open("rfdiffusion_metrics.json", "w") as f:
            json.dump({"designs": [], "fixed_residues": []}, f)
        with open("passing_designs.txt", "w") as f:
            pass
        with open("filter_summary.json", "w") as f:
            json.dump({"n_total": 0, "n_passing": 0, "n_filtered": 0,
                        "min_hotspot_frac": args.min_hotspot_frac}, f, indent=2)
        sys.exit(0)

    # ── Resolve contigs ──────────────────────────────────────────────────
    resolved_contigs = resolve_contigs(args.contigs, args.input_pdb)
    print(f"Raw contigs:      {args.contigs}")
    print(f"Resolved contigs: {resolved_contigs}")

    expected_rec_len, expected_eff_len = get_expected_chain_lengths(
        resolved_contigs, args.receptor_chain, args.effector_chain
    )
    print(f"Expected receptor length: {expected_rec_len}")
    print(f"Expected effector length: {expected_eff_len}")

    design_residues, fixed_residues = parse_design_region(
        resolved_contigs, args.receptor_chain
    )
    print(f"Design region residues: {sorted(design_residues)} ({len(design_residues)} total)")
    print(f"Fixed receptor residues: {len(fixed_residues)}")

    # ── Parse hotspot residues ───────────────────────────────────────────
    hotspot_residues = set()
    if args.hotspot:
        for token in args.hotspot.split(","):
            token = token.strip()
            if not token:
                continue
            num_part = token[1:] if token[0].isalpha() else token
            try:
                hotspot_residues.add(int(num_part))
            except ValueError:
                pass
    print(f"Hotspot (effector) residues: {len(hotspot_residues)}")

    # ── Load input PDB receptor atoms (RMSD reference) ───────────────────
    input_chains = count_chains(args.input_pdb)
    if len(input_chains) >= 2:
        input_rec_chain = (args.receptor_chain
                           if args.receptor_chain in input_chains
                           else input_chains[0])
        input_rec_atoms = read_ca_atoms(args.input_pdb, chain=input_rec_chain)
    else:
        input_all = read_ca_atoms(args.input_pdb)
        seg1, seg2, _ = split_at_chain_break(input_all)
        input_rec_atoms, _, _ = identify_segments(
            seg1, seg2, expected_rec_len, expected_eff_len
        )
    print(f"Input receptor Cα atoms: {len(input_rec_atoms)}")

    # Pre-parse the contig segment descriptors once.  Same input for every
    # design, used by calc_motif_and_region_metrics to find input gaps.
    seg_descs_for_metrics = _parse_receptor_segments(resolved_contigs, args.receptor_chain)

    os.makedirs("split", exist_ok=True)

    # ── Analyse each design ──────────────────────────────────────────────
    metrics = []
    passing = []
    filtered = []

    for dpdb in design_files:
        name = os.path.basename(dpdb)
        all_atoms = read_ca_atoms(dpdb)

        # Always use chain-break-split + length-based identify_segments,
        # regardless of how many chains the design PDB has.  RFDiffusion
        # renumbers output chains alphabetically (A, B, ...) in the order
        # they appear in the contig — these letters do NOT match the
        # user's input chain letters (e.g. receptor_chain="B" in the
        # input PDB becomes chain "A" in RFDiffusion output if it's the
        # first chain in the contig).  The previous chain-ID lookup
        # silently produced wrong results when input chain letters
        # weren't already A/B in the right order.
        seg1, seg2, break_dist = split_at_chain_break(all_atoms)
        rec_atoms, eff_atoms, _ = identify_segments(
            seg1, seg2, expected_rec_len, expected_eff_len
        )
        print(f"  {name}: chain break at {break_dist:.1f} Å, "
              f"seg1={len(seg1)} seg2={len(seg2)}")

        print(f"  {name}: {len(all_atoms)} total Cα, rec={len(rec_atoms)}, eff={len(eff_atoms)}")

        resnum_map, segment_info = build_receptor_resnum_map(rec_atoms, resolved_contigs, args.receptor_chain)
        contacts = find_contacts(rec_atoms, eff_atoms, args.rfdiff_contact_cutoff)
        n_contacts = len(contacts)

        rec_contact_resnums = sorted(set(
            r for r in (resnum_map.get(c[0], -1) for c in contacts) if r >= 0
        ))
        eff_contact_resnums = sorted(set(eff_atoms[c[1]]["resnum"] for c in contacts))

        # Paired contacts: unique (receptor_resnum, effector_resnum) tuples
        contact_pairs = sorted(set(
            (resnum_map.get(c[0], -1), eff_atoms[c[1]]["resnum"])
            for c in contacts
            if resnum_map.get(c[0], -1) >= 0
        ))

        # Classify contacts: "in design" = mapped residue NOT in fixed set
        n_in = sum(1 for c in contacts
                   if resnum_map.get(c[0], -1) >= 0
                   and resnum_map[c[0]] not in fixed_residues)
        n_out = n_contacts - n_in
        frac_in = n_in / n_contacts if n_contacts > 0 else 0.0

        # Design-region coverage: fraction of this design's de novo atoms
        # that have >= 1 contact
        per_design_design_indices = [
            idx for idx in range(len(rec_atoms))
            if resnum_map.get(idx) is not None
            and resnum_map[idx] not in fixed_residues
        ]
        n_design_this = len(per_design_design_indices)
        design_contacted_indices = set(
            c[0] for c in contacts
            if resnum_map.get(c[0]) is not None
            and resnum_map[c[0]] not in fixed_residues
        )
        design_coverage = (len(design_contacted_indices) / n_design_this
                           if n_design_this > 0 else 0.0)

        (motif_rmsd,
         endpoint_dist_input,
         endpoint_dist_design,
         com_displacement) = calc_motif_and_region_metrics(
            input_rec_atoms, rec_atoms, resnum_map, fixed_residues,
            seg_descs_for_metrics, segment_info["per_region_atom_ranges"],
        )
        passes = frac_in >= args.min_hotspot_frac if n_contacts > 0 else False

        # Extract design-region Cα coords for inter-design clustering
        design_region_coords = []
        per_design_design_residues = []
        for idx in per_design_design_indices:
            a = rec_atoms[idx]
            design_region_coords.append([round(a["x"], 3),
                                         round(a["y"], 3),
                                         round(a["z"], 3)])
            per_design_design_residues.append(resnum_map[idx])  # synthetic resnum

        # Per-region coords: list of coord lists, one per de novo region
        per_region_coords = []
        for rstart, rend in segment_info["per_region_atom_ranges"]:
            region_coords = []
            for idx in range(rstart, rend):
                if idx < len(rec_atoms):
                    a = rec_atoms[idx]
                    region_coords.append([round(a["x"], 3),
                                          round(a["y"], 3),
                                          round(a["z"], 3)])
            per_region_coords.append(region_coords)

        write_split_pdb(dpdb, os.path.join("split", name),
                        expected_rec_len, expected_eff_len)

        metrics.append({
            "design": name,
            "n_receptor_residues": len(rec_atoms),
            "n_contact_pairs": n_contacts,
            "n_unique_receptor_residues": len(rec_contact_resnums),
            "n_unique_effector_residues": len(eff_contact_resnums),
            "n_contacts_in_design": n_in,
            "n_contacts_outside_design": n_out,
            "frac_contacts_in_design": round(frac_in, 4),
            "design_region_coverage": round(design_coverage, 4),
            "motif_rmsd": round(motif_rmsd, 3) if not np.isnan(motif_rmsd) else None,
            "endpoint_distance_input": [
                round(v, 3) if not np.isnan(v) else None for v in endpoint_dist_input
            ],
            "endpoint_distance_design": [
                round(v, 3) if not np.isnan(v) else None for v in endpoint_dist_design
            ],
            "design_region_com_displacement": [
                round(v, 3) if not np.isnan(v) else None for v in com_displacement
            ],
            "receptor_contact_residues": rec_contact_resnums,
            "effector_contact_residues": eff_contact_resnums,
            "contact_pairs": contact_pairs,
            "design_region_coords": design_region_coords,
            "per_design_design_residues": per_design_design_residues,
            "per_region_lengths": segment_info["per_region_lengths"],
            "per_region_coords": per_region_coords,
            "receptor_position_order": segment_info["receptor_position_order"],
            "passes_filter": passes,
        })

        (passing if passes else filtered).append(name)

        motif_str = f"{motif_rmsd:.2f}" if not np.isnan(motif_rmsd) else "N/A"
        com_summary = ",".join(
            f"{v:.1f}" if not np.isnan(v) else "N/A" for v in com_displacement
        ) or "N/A"
        status = "PASS" if passes else "FAIL"
        print(f"    contacts={n_contacts} (in={n_in}, out={n_out}, "
              f"frac={frac_in:.2f}) motif_rmsd={motif_str}Å "
              f"com_disp=[{com_summary}]Å [{status}]")

    # Write outputs.  Note: per-design design residues live inside each
    # entry of `designs` (as `per_design_design_residues`, computed via
    # build_receptor_resnum_map) and are the authoritative source.  The
    # old top-level `design_residues` field — gap arithmetic from
    # parse_design_region on the input PDB — was misleading because it
    # only matched reality for fixed-length contigs and didn't reflect
    # what RFDiffusion actually produced.  Removed; consumers should
    # use the per-design field.
    with open("rfdiffusion_metrics.json", "w") as f:
        json.dump({
            "designs": metrics,
            "fixed_residues": sorted(fixed_residues),
            "hotspot_residues": sorted(hotspot_residues),
            "rfdiff_contact_cutoff": args.rfdiff_contact_cutoff,
            "min_hotspot_frac": args.min_hotspot_frac,
        }, f, indent=2)

    with open("passing_designs.txt", "w") as f:
        for name in passing:
            f.write(name + "\n")

    with open("filter_summary.json", "w") as f:
        json.dump({
            "n_total": len(design_files),
            "n_passing": len(passing),
            "n_filtered": len(filtered),
            "min_hotspot_frac": args.min_hotspot_frac,
            "filtered_designs": filtered,
        }, f, indent=2)

    print(f"\nFilter summary: {len(passing)}/{len(design_files)} designs pass "
          f"(threshold: {args.min_hotspot_frac:.0%} contacts in design region)")

    if not passing:
        print("WARNING: No designs passed the filter!", file=sys.stderr)


if __name__ == "__main__":
    main()