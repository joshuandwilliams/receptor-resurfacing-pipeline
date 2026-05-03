"""Strategy 4a contig derivation — interface-facing side chains, geometric.

Selects binder residues that face outward (Cβ-vector heuristic) and
contact the target (heavy-atom distance). Interface-facing residues
become the *designable* region (handed to RFDiffusion); non-interface
residues become *anchors* (kept native).

Short internal anchor blocks (≤ --merge-threshold residues) that sit
between two design regions are absorbed into the surrounding design
("bridging"), so a periodic helix-face contact pattern collapses into
one contiguous designable stretch instead of dozens of single-residue
fragments. Bridging only applies to anchors with a design region on
both sides — a short anchor at either terminus is preserved.

Bridging here is a fixed-threshold approach: any internal anchor
shorter than the threshold is absorbed. A future enhancement could
instead use a density-based approach — slide a window along the
binder and mark high-design-residue-density windows as design
regardless of the lengths of internal anchor runs they contain. Not
implemented now; flagged as a forward pointer.

See experiments/README.md, Stage 4 for the strategy taxonomy.
"""

from __future__ import annotations

import experiments._path_setup  # noqa: F401  - adds bin/ to sys.path

import argparse
import shlex
import sys
from pathlib import Path

try:
    from Bio.PDB import PDBParser
except ModuleNotFoundError:
    sys.exit(
        "BioPython is required (Bio.PDB.PDBParser).  "
        "Activate the LRR_Resurface mamba env, or install biopython."
    )

from boltz2_negative_steering import (
    find_contact_residues_heavy,
    is_surface_exposed,
    read_residue_heavy_atoms,
)


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Derive a Strategy 4a (interface-facing side chains) "
                    "RFDiffusion contig from an assembled complex PDB."
    )
    p.add_argument("--complex-pdb", required=True, type=Path,
                   help="Path to assembled complex PDB.")
    p.add_argument("--binder-chain", required=True,
                   help="Chain ID of the binder.")
    p.add_argument("--target-chain", required=True,
                   help="Chain ID of the target.")
    p.add_argument("--output", required=True, type=Path,
                   help="Output path for the contig file.")
    p.add_argument("--cutoff", type=float, default=5.0,
                   help="Heavy-atom distance cutoff in angstrom (default: 5.0).")
    p.add_argument("--min-multiplier", type=float, default=0.7,
                   help="Lower bound of design region as multiple of native "
                        "design length (default: 0.7).")
    p.add_argument("--max-multiplier", type=float, default=1.5,
                   help="Upper bound of design region as multiple of native "
                        "design length (default: 1.5).")
    p.add_argument("--merge-threshold", type=int, default=3,
                   help="Maximum length of an internal anchor block (one "
                        "bordered by design regions on both sides) that "
                        "will be absorbed into the surrounding design "
                        "(default: 3).")
    return p.parse_args(argv)


def list_chains(pdb_path: Path) -> list[str]:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", str(pdb_path))
    chains: list[str] = []
    for model in structure:
        for chain in model:
            cid = chain.id
            if cid and cid not in chains:
                chains.append(cid)
        break
    return chains


def design_token(native_length: int, min_mult: float, max_mult: float) -> str:
    if native_length == 1:
        return "1-1"
    min_len = max(1, round(native_length * min_mult))
    max_len = round(native_length * max_mult)
    if max_len < min_len:
        max_len = min_len
    return f"{min_len}-{max_len}"


def segment_binder(
    binder_resnums_sorted: list[int],
    interface_set: set[int],
) -> list[tuple[str, int, int]]:
    """Walk the binder N-to-C and group residues into (kind, start, end)
    segments where kind is 'design' (interface-facing) or 'anchor'
    (everything else).  Numbering gaps in the binder break a segment
    just like a kind change does.
    """
    segments: list[tuple[str, int, int]] = []
    if not binder_resnums_sorted:
        return segments

    cur_kind: str | None = None
    cur_start: int | None = None
    prev: int | None = None
    for r in binder_resnums_sorted:
        kind = "design" if r in interface_set else "anchor"
        if cur_kind is None:
            cur_kind, cur_start, prev = kind, r, r
            continue
        if kind != cur_kind or r != prev + 1:
            segments.append((cur_kind, cur_start, prev))
            cur_kind, cur_start = kind, r
        prev = r
    segments.append((cur_kind, cur_start, prev))
    return segments


def bridge_segments(
    segments: list[tuple[str, int, int]],
    merge_threshold: int,
) -> tuple[list[tuple[str, int, int]], list[dict]]:
    """Iteratively absorb short internal anchor segments into surrounding
    design segments.  An anchor is "internal" iff it has a design segment
    on BOTH sides; terminus anchors are never absorbed.

    Each pass scans left-to-right and absorbs the first eligible anchor;
    the loop repeats until a full pass finds nothing to absorb. Returns
    (new_segments, bridges), where each bridge dict records the absorbed
    anchor's range and the merged design's final range.
    """
    bridges: list[dict] = []
    segs = list(segments)

    while True:
        merged_this_pass = False
        for i in range(1, len(segs) - 1):
            kind, start, end = segs[i]
            if kind != "anchor":
                continue
            if (end - start + 1) > merge_threshold:
                continue
            left = segs[i - 1]
            right = segs[i + 1]
            if left[0] != "design" or right[0] != "design":
                continue

            merged = ("design", left[1], right[2])
            bridges.append({
                "absorbed_start": start,
                "absorbed_end": end,
                "absorbed_length": end - start + 1,
                "merged_start": merged[1],
                "merged_end": merged[2],
            })
            segs = segs[:i - 1] + [merged] + segs[i + 2:]
            merged_this_pass = True
            break

        if not merged_this_pass:
            break

    return segs, bridges


def main() -> None:
    args = parse_args()

    if not args.complex_pdb.is_file():
        sys.exit(f"Complex PDB not found: {args.complex_pdb}")

    chains = list_chains(args.complex_pdb)
    missing = [c for c in (args.binder_chain, args.target_chain) if c not in chains]
    if missing:
        sys.exit(
            f"Chain(s) not found in {args.complex_pdb}: {', '.join(missing)}.\n"
            f"Chains present: {', '.join(chains) or '(none)'}.\n"
            "To list chains in a PDB:\n"
            f"  grep '^ATOM' '{args.complex_pdb}' | awk '{{print $5}}' | sort -u"
        )

    # Use the negsteer reader so seq_index values line up with what
    # find_contact_residues_heavy returns.
    residues = read_residue_heavy_atoms(args.complex_pdb)
    binder_residues = [r for r in residues if r.chain == args.binder_chain]
    target_residues = [r for r in residues if r.chain == args.target_chain]
    n_binder = len(binder_residues)
    n_target = len(target_residues)

    seq_to_resnum = {r.seq_index: r.resnum for r in binder_residues}
    seq_to_residue = {r.seq_index: r for r in binder_residues}

    contacts = find_contact_residues_heavy(
        args.complex_pdb,
        args.binder_chain,
        args.target_chain,
        args.cutoff,
    )

    interface_resnums: list[int] = []
    for seq_idx, _dist in contacts:
        residue = seq_to_residue[seq_idx]
        if is_surface_exposed(residue, binder_residues):
            interface_resnums.append(seq_to_resnum[seq_idx])
    interface_resnums = sorted(set(interface_resnums))

    if not interface_resnums:
        sys.exit(
            "No interface-facing residues found on the binder. "
            "Try a larger --cutoff or check that the chain assignments are correct."
        )

    binder_resnums_sorted = sorted(r.resnum for r in binder_residues)
    interface_set = set(interface_resnums)

    raw_segments = segment_binder(binder_resnums_sorted, interface_set)
    final_segments, bridges = bridge_segments(raw_segments, args.merge_threshold)

    tokens: list[str] = []
    for kind, start, end in final_segments:
        if kind == "anchor":
            tokens.append(f"{args.binder_chain}{start}-{end}")
        else:
            native_len = end - start + 1
            tokens.append(design_token(native_len, args.min_multiplier, args.max_multiplier))

    binder_block = "/".join(tokens)
    contig = f"{binder_block} {args.target_chain}"

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(contig + "\n")

    notes_path = args.output.with_suffix(".notes")
    invocation = " ".join(shlex.quote(a) for a in sys.argv)

    final_anchors = [s for s in final_segments if s[0] == "anchor"]
    final_designs = [s for s in final_segments if s[0] == "design"]

    notes: list[str] = [
        "# Strategy 4a contig derivation",
        "",
        f"Invocation: {invocation}",
        f"Complex PDB: {args.complex_pdb}",
        f"Binder chain: {args.binder_chain}",
        f"Target chain: {args.target_chain}",
        f"Cutoff (A): {args.cutoff}",
        f"Min multiplier: {args.min_multiplier}",
        f"Max multiplier: {args.max_multiplier}",
        f"Merge threshold: {args.merge_threshold}",
        "",
        "## Binder summary",
        f"Total binder residues: {n_binder} "
        f"(resnums {binder_resnums_sorted[0]}-{binder_resnums_sorted[-1]})",
        f"Interface-facing residues identified (raw, pre-bridging): "
        f"{len(interface_resnums)}",
        f"Anchor blocks absorbed during bridging: {len(bridges)}",
        f"Final design regions count: {len(final_designs)}",
        f"Final anchor blocks count: {len(final_anchors)}",
        "",
        "## Target summary",
        f"Total target residues: {n_target}",
        "",
        "## Anchor blocks (non-interface, kept native)",
    ]
    if not final_anchors:
        notes.append("  (none)")
    else:
        for _kind, start, end in final_anchors:
            notes.append(
                f"  {args.binder_chain}{start}-{end}  (length {end - start + 1})"
            )

    notes.append("")
    notes.append("## Design regions (interface-facing, handed to RFDiffusion)")
    if not final_designs:
        notes.append("  (none)")
    else:
        for _kind, start, end in final_designs:
            nlen = end - start + 1
            tok = design_token(nlen, args.min_multiplier, args.max_multiplier)
            notes.append(
                f"  resnums {start}-{end}  native_length={nlen}  -> {tok}"
            )

    notes.append("")
    notes.append("## Bridging")
    if not bridges:
        notes.append("No anchor blocks absorbed.")
    else:
        for b in bridges:
            notes.append(
                f"Absorbed {args.binder_chain}{b['absorbed_start']}-"
                f"{b['absorbed_end']} (native length {b['absorbed_length']}) "
                f"into surrounding design spanning "
                f"{args.binder_chain}{b['merged_start']}-{b['merged_end']}"
            )

    notes += [
        "",
        "## Final contig",
        contig,
        "",
        "## Notes on helpers",
        "Outward-facing test: is_surface_exposed (bin/boltz2_negative_steering.py).",
        "  Cbeta-vector heuristic: (Calpha->Cbeta).(Calpha->neighbour-centroid) < 0.",
        "  Glycine (no Cbeta) is excluded.",
        "Contact test: find_contact_residues_heavy (bin/boltz2_negative_steering.py).",
        "  Compares all heavy atoms (backbone + side chain) on both sides;",
        "  the 4a strategy doc specifies side-chain heavy on the binder,",
        "  but we reuse the helper rather than reimplementing.",
    ]

    notes_path.write_text("\n".join(notes) + "\n")

    print(f"Contig written to {args.output}: {contig}")
    print(f"Notes written to {notes_path}")


if __name__ == "__main__":
    main()
