#!/usr/bin/env python
"""
Sequence correction for ProteinMPNN outputs.

Handles the case where RFDiffusion merges both chains into a single chain (A),
so ProteinMPNN outputs one continuous sequence without chain separators.

Supports variable-length inpaintings (e.g. "20-30" de novo segments) by using
the contig structure to identify fixed anchor regions, then aligning the MPNN
output against native anchors rather than assuming fixed positions.

Modes:
  gen_fixed_positions: Generate fixed_positions JSONL for MPNN (PDB-aware)
  correct:            Process MPNN FASTAs, restore native at non-designed sites, write AF2 FASTAs
"""

import argparse
import csv
import json
import os


# ---------------------------------------------------------------------------
# Contig parsing
# ---------------------------------------------------------------------------

def parse_contig_segments(contigs, receptor_chain="A"):
    """Parse the contig string into an ordered list of segments for the receptor.

    Returns a list of dicts, each with:
      type:  'fixed' or 'denovo'
      start: PDB start residue (fixed only)
      end:   PDB end residue (fixed only)
      min_len: minimum de novo length (denovo only)
      max_len: maximum de novo length (denovo only)

    Example for "A1-400/20-30/A426-435":
      [{'type':'fixed', 'start':1, 'end':400},
       {'type':'denovo', 'min_len':20, 'max_len':30},
       {'type':'fixed', 'start':426, 'end':435}]
    """
    parts = contigs.split()
    receptor_part = None
    for part in parts:
        # Find block that contains receptor chain segments
        if receptor_chain.upper() in part.upper():
            receptor_part = part
            break
    if not receptor_part:
        return []

    segments = []
    for seg in receptor_part.split("/"):
        seg = seg.strip()
        if not seg or seg == "0":
            continue

        if seg[0].isalpha() and seg[0].upper() == receptor_chain.upper():
            # Fixed receptor segment: A1-400, A403, A426-435
            rest = seg[1:]
            if "-" in rest:
                parts_r = rest.split("-")
                start, end = int(parts_r[0]), int(parts_r[1])
            else:
                start = end = int(rest)
            segments.append({'type': 'fixed', 'start': start, 'end': end})
        elif seg[0].isdigit():
            # De novo segment: 1-1, 20-30, 5-5
            if "-" in seg:
                parts_r = seg.split("-")
                min_len, max_len = int(parts_r[0]), int(parts_r[1])
            else:
                min_len = max_len = int(seg)
            segments.append({'type': 'denovo', 'min_len': min_len, 'max_len': max_len})

    return segments


def get_fixed_residue_set(segments, receptor_start_pdb=1):
    """Get the set of 0-indexed receptor positions that are fixed (not de novo).

    Only works for fixed-length contigs (all de novo segments have min_len == max_len).
    For variable-length contigs, returns None to signal that PDB-based detection is needed.
    """
    has_variable = any(
        s['type'] == 'denovo' and s['min_len'] != s['max_len']
        for s in segments
    )
    if has_variable:
        return None

    fixed_positions = set()
    for seg in segments:
        if seg['type'] == 'fixed':
            for r in range(seg['start'], seg['end'] + 1):
                fixed_positions.add(r - receptor_start_pdb)
    return fixed_positions


def get_native_anchor_regions(segments):
    """Return list of (native_start_0idx, native_end_0idx) for each fixed segment.

    These are 0-indexed positions within the native receptor sequence corresponding
    to each fixed segment (derived from PDB numbering minus receptor_start_pdb).
    """
    anchors = []
    for seg in segments:
        if seg['type'] == 'fixed':
            anchors.append((seg['start'], seg['end']))
    return anchors


# ---------------------------------------------------------------------------
# PDB utilities
# ---------------------------------------------------------------------------

def get_pdb_chain_residues(pdb_path):
    """Get residue numbers per chain from PDB file."""
    chains = {}
    seen = set()
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM"):
                chain = line[21]
                resnum = int(line[22:26].strip())
                key = (chain, resnum)
                if key not in seen:
                    seen.add(key)
                    chains.setdefault(chain, []).append(resnum)
    return chains


def get_pdb_sequence(pdb_path, chain_id):
    """Extract amino acid sequence from a PDB file for a given chain."""
    THREE_TO_ONE = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
        'MSE': 'M', 'SEC': 'U', 'PYL': 'O',
    }
    residues = {}
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") and line[21] == chain_id:
                resnum = int(line[22:26].strip())
                resname = line[17:20].strip()
                if resnum not in residues:
                    residues[resnum] = resname
    resnums = sorted(residues.keys())
    seq = ''.join(THREE_TO_ONE.get(residues[r], 'X') for r in resnums)
    return seq, resnums


# ---------------------------------------------------------------------------
# FASTA utilities
# ---------------------------------------------------------------------------

def parse_mpnn_fasta(fasta_path):
    """Parse a multi-entry FASTA file."""
    entries = []
    header, seq_lines = None, []
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    entries.append((header, "".join(seq_lines)))
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)
    if header is not None:
        entries.append((header, "".join(seq_lines)))
    return entries


def extract_mpnn_score(header):
    """Extract MPNN score from FASTA header."""
    for part in header.split(","):
        part = part.strip()
        if part.startswith("score="):
            try:
                return float(part.split("=")[1])
            except (ValueError, IndexError):
                pass
    return 0.0


def split_sequence(seq, effector_len):
    """Split a concatenated MPNN sequence into (effector, receptor) parts.

    ProteinMPNN concatenates its per-chain output using '/' or ':' as
    the inter-chain separator.  Which chain comes first depends on how
    the PDB was written.  In the receptor-resurfacing pipeline the
    split PDB has receptor=A and effector=B, so the order in MPNN's
    output is receptor then effector — the opposite of what this
    function originally assumed.  A previous implementation hard-coded
    "effector first" and silently produced corrected_receptor strings
    built from the *effector* residues, which then fed downstream as
    if they were MPNN's receptor design.

    To stay robust against chain ordering, when a separator is
    present we assign the part whose length is closest to
    `effector_len` as the effector and the other as the receptor.
    The no-separator fallback keeps the legacy positional split for
    compatibility; it is only reached if MPNN stops emitting
    separators, which the current pipeline does not do.

    Returns (effector_seq, receptor_seq) in that order.
    """
    for sep in ["/", ":"]:
        if sep in seq:
            parts = seq.split(sep)
            if len(parts) >= 2:
                cand_eff = parts[0]
                cand_rec = sep.join(parts[1:])
                # Pick whichever half is closer to effector_len as the
                # effector.  This makes the function correct regardless
                # of which chain MPNN wrote first.
                if abs(len(cand_rec) - effector_len) < abs(len(cand_eff) - effector_len):
                    cand_eff, cand_rec = cand_rec, cand_eff
                return cand_eff, cand_rec

    # No separator -- legacy positional split (effector first).
    # Not reached in the current pipeline; MPNN always emits '/'.
    return seq[:effector_len], seq[effector_len:]


# ---------------------------------------------------------------------------
# Anchor-based sequence alignment
# ---------------------------------------------------------------------------

def align_to_native_by_anchors(mpnn_receptor, native_receptor, segments, receptor_start_pdb):
    """Align MPNN receptor output to native using fixed segments as anchors.

    For each fixed segment in the contig, the native sequence at those positions
    is known. We walk through the contig segments in order, consuming the
    appropriate number of residues from the MPNN output for fixed segments
    (matching native length) and whatever remains for de novo segments.

    Returns:
      corrected_receptor: str - native at fixed positions, MPNN at de novo positions
      designed_residues:  str - just the residues at de novo positions
      native_at_design:   str - native residues at positions replaced by de novo
                                (empty string if insertion is longer than original gap)
      design_region_length_observed:
                          str - per-segment "|"-joined actual lengths of
                                the de novo region(s) in this design
    """
    # Calculate expected lengths for each segment
    # Fixed segments: length = end - start + 1 (from native)
    # De novo segments: length = whatever is left after accounting for fixed segments
    total_fixed = sum(s['end'] - s['start'] + 1 for s in segments if s['type'] == 'fixed')
    total_denovo = len(mpnn_receptor) - total_fixed

    if total_denovo < 0:
        print(f"    WARNING: MPNN receptor ({len(mpnn_receptor)}) shorter than "
              f"total fixed ({total_fixed}). Falling back to length-based split.")
        total_denovo = 0

    # Count de novo segments to distribute length
    denovo_segments = [s for s in segments if s['type'] == 'denovo']
    n_denovo = len(denovo_segments)

    # Compute the per-segment de novo allocation up front, mirroring the
    # R7 fix in rfdiffusion_filter.build_receptor_resnum_map.  Strategy:
    #
    #   1. Fixed-length segments (min_len == max_len) get exactly their
    #      declared length.  These are spec-mandated and must not be
    #      proportionally split.
    #   2. Variable-length segments share the leftover budget.  With one
    #      var segment, it absorbs the remainder clamped to its bounds.
    #      With several, the leftover is distributed proportionally by
    #      max_len, every allocation is clamped to [min_len, max_len],
    #      and the final var segment absorbs any rounding remainder.
    #
    # This is more correct than the previous purely-proportional approach
    # because it can't violate fixed-length specs (e.g. "6-6" must be 6)
    # and every segment is bounds-checked individually.
    seg_index_to_alloc = {}
    if n_denovo > 0:
        denovo_idx_pairs = [
            (i, s) for i, s in enumerate(segments) if s['type'] == 'denovo'
        ]
        fixed_len_pairs = [
            (i, s) for i, s in denovo_idx_pairs if s['min_len'] == s['max_len']
        ]
        var_len_pairs = [
            (i, s) for i, s in denovo_idx_pairs if s['min_len'] != s['max_len']
        ]

        for i, s in fixed_len_pairs:
            seg_index_to_alloc[i] = s['min_len']
        used_by_fixed = sum(s['min_len'] for _, s in fixed_len_pairs)
        remaining_for_var = max(0, total_denovo - used_by_fixed)

        if len(var_len_pairs) == 0:
            pass  # all de novo segments are fixed-length, nothing left to do
        elif len(var_len_pairs) == 1:
            i, s = var_len_pairs[0]
            alloc = max(s['min_len'], min(s['max_len'], remaining_for_var))
            if alloc != remaining_for_var:
                print(f"    WARNING: variable denovo segment {i} got {alloc} "
                      f"residues after clamping to [{s['min_len']}, {s['max_len']}], "
                      f"but receptor length suggests {remaining_for_var}. "
                      f"Per-region split may not match RFDiffusion's actual "
                      f"output for this design.")
            seg_index_to_alloc[i] = alloc
        else:
            total_max = sum(s['max_len'] for _, s in var_len_pairs)
            allocated = 0
            for k, (i, s) in enumerate(var_len_pairs):
                if k < len(var_len_pairs) - 1:
                    if total_max > 0:
                        alloc = round(remaining_for_var * s['max_len'] / total_max)
                    else:
                        alloc = remaining_for_var // len(var_len_pairs)
                    alloc = max(s['min_len'], min(s['max_len'], alloc))
                    alloc = min(alloc, remaining_for_var - allocated)
                    seg_index_to_alloc[i] = alloc
                    allocated += alloc
                else:
                    seg_index_to_alloc[i] = remaining_for_var - allocated

    # Sanity check: the sum of fixed-segment lengths and de novo allocations
    # should equal the MPNN receptor length.  A mismatch means the contig
    # spec can't account for what RFDiffusion actually produced — usually a
    # fixed-length de novo segment whose declared length doesn't match the
    # design's true length, or a contig that's structurally wrong for this
    # PDB.  We warn loudly but don't raise: the walk below will silently
    # drop or short-read the unaccounted residues, and the caller should
    # treat the resulting corrected sequence as suspect.
    total_alloc = sum(seg_index_to_alloc.values())
    expected_total = total_fixed + total_alloc
    if expected_total != len(mpnn_receptor):
        diff = len(mpnn_receptor) - expected_total
        print(f"    WARNING: MPNN receptor length ({len(mpnn_receptor)}) does "
              f"not match contig accounting (fixed={total_fixed} + "
              f"denovo={total_alloc} = {expected_total}, diff={diff:+d}). "
              f"This usually means a fixed-length denovo segment's declared "
              f"length disagrees with the actual design, or the contig spec "
              f"doesn't match this PDB. Corrected sequence may be truncated "
              f"or have unaccounted residues silently dropped.")

    # Walk through segments, consuming from mpnn_receptor
    corrected_parts = []
    designed_parts = []
    native_design_parts = []
    mpnn_pos = 0

    for seg_idx, seg in enumerate(segments):
        if seg['type'] == 'fixed':
            seg_len = seg['end'] - seg['start'] + 1
            # Use native sequence for fixed positions
            native_start_0 = seg['start'] - receptor_start_pdb
            native_chunk = native_receptor[native_start_0:native_start_0 + seg_len]

            # Consume the corresponding positions from MPNN output but replace
            # with native (this is the "correction" — fixing any MPNN drift)
            corrected_parts.append(native_chunk)
            mpnn_pos += seg_len

        elif seg['type'] == 'denovo':
            alloc = seg_index_to_alloc.get(seg_idx, 0)
            denovo_chunk = mpnn_receptor[mpnn_pos:mpnn_pos + alloc]
            corrected_parts.append(denovo_chunk)
            designed_parts.append(denovo_chunk)

            # What was native at these positions (if they existed)?
            # For variable-length inpaintings, there may be no native equivalent
            # (e.g. native gap from 401-425 didn't exist, we're inserting new residues)
            # We report what the native had between the flanking fixed segments
            native_design_parts.append("")  # will be filled below

            mpnn_pos += alloc

    corrected_receptor = "".join(corrected_parts)
    designed_residues = "|".join(designed_parts)

    # Figure out the native residues at the de novo positions
    # These are the residues in the native between consecutive fixed segments
    #
    # Boundary handling:
    # - "interior" denovo: flanked by fixed segments on both sides → take
    #   the slice between them.
    # - "leading" denovo: no fixed segment before it → slice from the
    #   start of the native receptor (gap_start_0 = 0) up to the first
    #   following fixed segment.
    # - "trailing" denovo: no fixed segment after it → slice from the
    #   previous fixed segment's end to the end of the native receptor.
    # - denovo as the entire receptor (no fixed at all): use the whole
    #   native sequence.  Pathological but representable.
    #
    # Without these boundary cases, a contig like "B187-208/13-13/B223-262/6-6 C"
    # (fixed/denovo/fixed/denovo) silently drops the trailing denovo
    # region from `native_at_design`, producing a string like
    # "<region1>|" with an empty region 2.  Downstream plots that rely on
    # this column then under-report the native composition.
    native_at_design_parts = []
    for i, seg in enumerate(segments):
        if seg['type'] != 'denovo':
            continue
        # Find flanking fixed segments
        prev_end = None
        next_start = None
        for j in range(i - 1, -1, -1):
            if segments[j]['type'] == 'fixed':
                prev_end = segments[j]['end']
                break
        for j in range(i + 1, len(segments)):
            if segments[j]['type'] == 'fixed':
                next_start = segments[j]['start']
                break

        # Resolve gap boundaries with the boundary cases above.
        if prev_end is not None:
            gap_start_0 = prev_end + 1 - receptor_start_pdb
        else:
            gap_start_0 = 0   # leading denovo

        if next_start is not None:
            gap_end_0 = next_start - receptor_start_pdb
        else:
            gap_end_0 = len(native_receptor)   # trailing denovo

        if 0 <= gap_start_0 < gap_end_0 <= len(native_receptor):
            native_at_design_parts.append(native_receptor[gap_start_0:gap_end_0])
        else:
            native_at_design_parts.append("")

    native_at_design = "|".join(native_at_design_parts)

    # Per-segment observed lengths of the de novo region(s) in this design,
    # joined with "|".  This is the per-design ACTUAL length produced by
    # RFDiffusion + MPNN, which may differ from the contig spec range
    # (e.g. spec "20-40" but this design produced 27).  See callers, which
    # also write a sibling design_region_spec column carrying the spec.
    design_region_length_observed = "|".join(str(len(p)) for p in designed_parts)

    return corrected_receptor, designed_residues, native_at_design, design_region_length_observed


# ---------------------------------------------------------------------------
# Mode: Generate fixed positions JSONL for ProteinMPNN
# ---------------------------------------------------------------------------

def generate_fixed_positions(args):
    """Generate fixed_positions JSONL.

    RFDiffusion outputs a single chain A containing [effector | receptor].
    We fix all effector residues and all non-designed receptor residues.

    For variable-length inpaintings, we read the actual PDB to determine
    which residues are in the de novo region (residues between the fixed
    anchor endpoints that don't appear in the native numbering).
    """
    segments = parse_contig_segments(args.contigs, args.receptor_chain)
    effector_len = len(args.effector_seq)
    receptor_len = len(args.receptor_seq)

    print(f"Contig segments: {segments}")
    print(f"Expected chain layout: effector ({effector_len}) + receptor (variable)")

    # Build list of PDB files to process
    if args.pdb_file:
        pdb_files = [args.pdb_file]
    else:
        pdb_files = []
        for m in range(args.num_designs):
            p = os.path.join(args.rfdiff_dir, f"design_{m}.pdb")
            if os.path.exists(p):
                pdb_files.append(p)

    # Build set of fixed PDB residue numbers from the contig segments
    # These are the residue numbers that MUST be kept as native
    fixed_pdb_resnums = set()
    for seg in segments:
        if seg['type'] == 'fixed':
            for r in range(seg['start'], seg['end'] + 1):
                fixed_pdb_resnums.add(r)

    for pdb_path in pdb_files:
        chain_residues = get_pdb_chain_residues(pdb_path)
        chains = sorted(chain_residues.keys())
        pdb_stem = os.path.basename(pdb_path).replace(".pdb", "")

        if len(chains) == 1:
            # Single chain -- RFDiffusion merged everything into chain A
            chain_id = chains[0]
            all_resnums = sorted(chain_residues[chain_id])
            total_res = len(all_resnums)

            # The PDB has: [effector residues 1..eff_len] [receptor residues]
            # Receptor residues in the PDB are numbered sequentially, but
            # RFDiffusion may renumber them.  We identify fixed residues by
            # matching against the contig structure:
            #   - First effector_len residues → all fixed
            #   - Remaining residues → walk the contig segments to assign

            # Strategy: walk the PDB residues after the effector, walk the
            # contig segments, and match fixed segments by their expected length
            receptor_resnums = all_resnums[effector_len:]

            # Walk contig segments and consume receptor residues
            fixed_resnums = list(all_resnums[:effector_len])  # all effector fixed
            pos = 0
            for seg in segments:
                if seg['type'] == 'fixed':
                    seg_len = seg['end'] - seg['start'] + 1
                    # These residues are fixed
                    for j in range(seg_len):
                        if pos + j < len(receptor_resnums):
                            fixed_resnums.append(receptor_resnums[pos + j])
                    pos += seg_len
                elif seg['type'] == 'denovo':
                    # Figure out how many residues RFDiffusion actually placed here
                    # = total receptor residues minus all fixed segment lengths
                    # For single de novo region this is straightforward
                    total_fixed_receptor = sum(
                        s['end'] - s['start'] + 1 for s in segments if s['type'] == 'fixed'
                    )
                    total_denovo_actual = len(receptor_resnums) - total_fixed_receptor

                    # For multiple de novo regions, distribute proportionally
                    denovo_segs = [s for s in segments if s['type'] == 'denovo']
                    if len(denovo_segs) == 1:
                        actual_len = total_denovo_actual
                    else:
                        total_max = sum(s['max_len'] for s in denovo_segs)
                        if total_max > 0:
                            actual_len = round(total_denovo_actual * seg['max_len'] / total_max)
                        else:
                            actual_len = total_denovo_actual // len(denovo_segs)
                        actual_len = max(seg['min_len'], min(seg['max_len'], actual_len))

                    # De novo residues are NOT fixed (skip them)
                    pos += actual_len

            n_free = total_res - len(fixed_resnums)
            jsonl_path = os.path.join(args.output_dir, f"fixed_positions_{pdb_stem}.jsonl")
            with open(jsonl_path, "w") as f:
                f.write(json.dumps({pdb_stem: {chain_id: fixed_resnums}}) + "\n")

            print(f"  {pdb_stem}: {n_free} free / {len(fixed_resnums)} fixed "
                  f"(single chain {chain_id}, {total_res} total)")

        elif len(chains) >= 2:
            # Multi-chain case
            # The effector length is always known exactly (fixed), but the
            # receptor length varies due to de novo insertions.  We cannot
            # compare against native receptor_len because the actual receptor
            # chain may be much longer.  Instead, assign the chain whose
            # length is closest to effector_len as the effector, and the
            # other as the receptor.  For >2 chains, pick the best effector
            # match and assign the largest remaining chain as receptor.
            chain_lengths = {c: len(chain_residues[c]) for c in chains}

            # Find the chain closest to effector_len
            effector_chain_id = min(
                chains, key=lambda c: abs(chain_lengths[c] - effector_len)
            )
            # Receptor is the other chain (for 2-chain case) or the largest
            # remaining chain (for >2 chains)
            remaining = [c for c in chains if c != effector_chain_id]
            if remaining:
                receptor_chain_id = max(remaining, key=lambda c: chain_lengths[c])
            else:
                receptor_chain_id = None

            if receptor_chain_id is None:
                print(f"  WARNING: Could not identify chains for {pdb_stem}")
                continue

            fixed = {effector_chain_id: chain_residues[effector_chain_id]}

            # Walk contig segments to determine fixed receptor residues
            receptor_resnums = sorted(chain_residues[receptor_chain_id])
            fixed_receptor = []
            pos = 0
            total_fixed_receptor = sum(
                s['end'] - s['start'] + 1 for s in segments if s['type'] == 'fixed'
            )
            total_denovo_actual = len(receptor_resnums) - total_fixed_receptor

            for seg in segments:
                if seg['type'] == 'fixed':
                    seg_len = seg['end'] - seg['start'] + 1
                    for j in range(seg_len):
                        if pos + j < len(receptor_resnums):
                            fixed_receptor.append(receptor_resnums[pos + j])
                    pos += seg_len
                elif seg['type'] == 'denovo':
                    denovo_segs = [s for s in segments if s['type'] == 'denovo']
                    if len(denovo_segs) == 1:
                        actual_len = total_denovo_actual
                    else:
                        total_max = sum(s['max_len'] for s in denovo_segs)
                        if total_max > 0:
                            actual_len = round(total_denovo_actual * seg['max_len'] / total_max)
                        else:
                            actual_len = total_denovo_actual // len(denovo_segs)
                        actual_len = max(seg['min_len'], min(seg['max_len'], actual_len))
                    pos += actual_len

            fixed[receptor_chain_id] = fixed_receptor

            jsonl_path = os.path.join(args.output_dir, f"fixed_positions_{pdb_stem}.jsonl")
            with open(jsonl_path, "w") as f:
                f.write(json.dumps({pdb_stem: fixed}) + "\n")

            n_free = len(receptor_resnums) - len(fixed_receptor)
            print(f"  {pdb_stem}: {n_free} free / {len(fixed_receptor)} fixed receptor / "
                  f"{len(chain_residues[effector_chain_id])} fixed effector")


# ---------------------------------------------------------------------------
# Mode: Correct sequences and write AF2 FASTAs
# ---------------------------------------------------------------------------

def correct_sequences(args):
    native_lrr = args.receptor_seq
    native_effector = args.effector_seq
    receptor_len = len(native_lrr)
    effector_len = len(native_effector)

    segments = parse_contig_segments(args.contigs, args.receptor_chain)
    receptor_start_pdb = args.receptor_start_pdb

    # Build the contig spec string for the de novo region(s).  This is the
    # per-segment "min-max" range as written in the contig (e.g. "20-40"
    # for a single variable region, "5-5|10-30" for two regions), and is
    # invariant across all designs in this run.  Reported alongside the
    # per-design observed length so consumers can distinguish "what was
    # asked for" from "what RFDiffusion actually produced".
    design_region_spec = "|".join(
        f"{s['min_len']}-{s['max_len']}"
        for s in segments if s['type'] == 'denovo'
    )

    # For reporting: what were the native residues at the design positions?
    # (only meaningful for fixed-length de novo regions that replace existing residues)
    fixed_seg_lengths = sum(s['end'] - s['start'] + 1 for s in segments if s['type'] == 'fixed')

    print(f"Contig segments: {segments}")
    print(f"Fixed receptor residues: {fixed_seg_lengths}")
    print(f"Expected: effector ({effector_len}) + receptor (variable, native={receptor_len})")

    fasta_dir = os.path.join(args.output_dir, "af2_fastas")
    os.makedirs(fasta_dir, exist_ok=True)

    corrected_fasta_lines = []
    metadata_rows = []

    for m in range(args.num_designs):
        mpnn_fasta = os.path.join(args.mpnn_dir, f"design_{m}", "seqs", f"design_{m}.fa")
        if not os.path.exists(mpnn_fasta):
            print(f"  WARNING: {mpnn_fasta} not found, skipping")
            continue

        entries = parse_mpnn_fasta(mpnn_fasta)
        designed_entries = entries[1:]  # skip native/score=0 entry

        for n, (header, seq) in enumerate(designed_entries):
            if n >= args.num_seqs:
                break

            mpnn_score = extract_mpnn_score(header)

            # Split into (effector, receptor).  split_sequence picks
            # the assignment by length regardless of chain order in the
            # raw MPNN output, so no further swap is needed here.
            effector_pred, lrr_pred = split_sequence(seq, effector_len)

            # Align MPNN receptor to native using anchor-based approach
            corrected_lrr, designed_residues, native_at_design, design_region_length_observed = \
                align_to_native_by_anchors(lrr_pred, native_lrr, segments, receptor_start_pdb)

            # Count mutations at design positions — per region
            designed_parts = designed_residues.split("|")
            native_parts = native_at_design.split("|")
            changes_per_region = []
            for dp, np_ in zip(designed_parts, native_parts):
                if dp and np_:
                    min_len_cmp = min(len(dp), len(np_))
                    ch = sum(1 for d, na in zip(dp[:min_len_cmp], np_[:min_len_cmp])
                             if d != na)
                    ch += abs(len(dp) - len(np_))
                    changes_per_region.append(ch)
                elif dp:
                    changes_per_region.append(len(dp))
                else:
                    changes_per_region.append(0)
            # Any extra designed regions without native counterparts
            for dp in designed_parts[len(native_parts):]:
                changes_per_region.append(len(dp))
            num_changes = "|".join(str(c) for c in changes_per_region)

            fasta_name = f"design_{m}_seq_{n}"
            fasta_path = os.path.join(fasta_dir, f"{fasta_name}.fasta")
            with open(fasta_path, "w") as f:
                f.write(f">receptor\n{corrected_lrr}\n>effector\n{native_effector}\n")

            corrected_fasta_lines.append(
                f">{fasta_name}|mpnn={mpnn_score:.3f}|designed={designed_residues}"
                f"|design_region_length_observed={design_region_length_observed}"
                f"|design_region_spec={design_region_spec}|changes={num_changes}"
            )
            corrected_fasta_lines.append(f"{corrected_lrr}:{native_effector}")

            metadata_rows.append({
                "design": m, "seq": n, "mpnn_score": mpnn_score,
                "native_residues": native_at_design,
                "designed_residues": designed_residues,
                "num_changes": num_changes,
                "design_region_length_observed": design_region_length_observed,
                "design_region_spec": design_region_spec,
                "corrected_receptor": corrected_lrr,
                "receptor_total_length": len(corrected_lrr),
            })
            print(f"  design_{m}_seq_{n}: changes={num_changes}, "
                  f"design_region_length_observed={design_region_length_observed}, "
                  f"mpnn_score={mpnn_score:.3f}, "
                  f"total_receptor_len={len(corrected_lrr)}")

    with open(os.path.join(args.output_dir, "mpnn_corrected.fasta"), "w") as f:
        f.write("\n".join(corrected_fasta_lines) + "\n")

    if metadata_rows:
        csv_path = os.path.join(args.output_dir, "sequence_metadata.csv")
        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=metadata_rows[0].keys())
            writer.writeheader()
            writer.writerows(metadata_rows)

    print(f"\nDone: {len(metadata_rows)} corrected sequences -> {fasta_dir}/")


# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Sequence correction for ProteinMPNN outputs in receptor resurfacing pipeline"
    )
    parser.add_argument("--mode", default="correct",
                        choices=["correct", "gen_fixed_positions"])
    parser.add_argument("--mpnn_dir")
    parser.add_argument("--rfdiff_dir")
    parser.add_argument("--pdb_file", help="Single PDB file (bypasses design_N.pdb loop)")
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--receptor_seq", required=True)
    parser.add_argument("--effector_seq", required=True)
    parser.add_argument("--contigs", required=True)
    parser.add_argument("--receptor_chain", default="A")
    parser.add_argument("--receptor_start_pdb", type=int, default=1)
    parser.add_argument("--num_designs", type=int, required=True)
    parser.add_argument("--num_seqs", type=int, default=2)
    args = parser.parse_args()

    if args.mode == "gen_fixed_positions":
        assert args.rfdiff_dir or args.pdb_file, \
            "--rfdiff_dir or --pdb_file required for gen_fixed_positions"
        generate_fixed_positions(args)
    else:
        assert args.mpnn_dir, "--mpnn_dir required for correct mode"
        correct_sequences(args)


if __name__ == "__main__":
    main()