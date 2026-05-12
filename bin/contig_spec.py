"""
contig_spec.py
--------------
ContigSpec — RFDiffusion contig strings parsed into structured per-chain
segments.  Canonical input form: ``A1-10/5/A15-20 B`` (see memory
project_contig_string_format.md).

Tier 0 type — no upstream dependencies on other Phase 4 types.

ContigSpec represents the *resolved* form: after RFDiffusion has sampled
a specific length for each de novo segment, the contig becomes a
concrete layout.  The constraint form (with length ranges like ``5-7``)
is an input string to RFDiffusion and is NOT a runtime type — it's
plumbed through as text only.

This module supersedes the 4 existing contig parsers (consolidation
deferred to higher tiers):

  - bin/contig_utils.py:parse_block_segments / resolve_contigs /
    parse_design_region — the closest existing parser; ContigSpec
    absorbs its functionality here.  contig_utils.py keeps working
    until callers are migrated.
  - bin/derive_input_design_region.py:_parse_contigs
  - bin/haddock3_prepare.py:parse_contig_segments
  - bin/pipeline_correct_sequences.py:parse_contig_segments
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple


# ═════════════════════════════════════════════════════════════════════════
# Segment types
# ═════════════════════════════════════════════════════════════════════════


@dataclass(frozen=True)
class FixedSegment:
    """A fixed (anchored) segment: chain + native residue range
    (inclusive on both ends)."""
    start: int
    end: int

    def __post_init__(self) -> None:
        if self.start > self.end:
            raise ValueError(
                f"FixedSegment start={self.start} > end={self.end}"
            )

    @property
    def length(self) -> int:
        """Number of native residues in this fixed segment."""
        return self.end - self.start + 1


@dataclass(frozen=True)
class DeNovoSegment:
    """A de novo (design) segment: specific resolved length, no chain
    prefix (RFDiffusion fills these in).  Carries no native residue
    references — these positions have no native counterpart by
    definition."""
    length: int

    def __post_init__(self) -> None:
        if self.length < 0:
            raise ValueError(f"DeNovoSegment length={self.length} must be >= 0")


# ═════════════════════════════════════════════════════════════════════════
# ContigChain
# ═════════════════════════════════════════════════════════════════════════


@dataclass(frozen=True)
class ContigChain:
    """The contig segments for a single chain, in order."""
    chain_id: str
    segments: Tuple  # Tuple[Union[FixedSegment, DeNovoSegment], ...]

    def __post_init__(self) -> None:
        # Coerce a list to a tuple so frozen+hashable holds.
        if not isinstance(self.segments, tuple):
            object.__setattr__(self, "segments", tuple(self.segments))

    @property
    def total_length(self) -> int:
        """Total residues in the resolved chain
        (sum of fixed and de novo lengths)."""
        return sum(s.length for s in self.segments)

    @property
    def fixed_segments(self) -> List[FixedSegment]:
        return [s for s in self.segments if isinstance(s, FixedSegment)]

    @property
    def denovo_segments(self) -> List[DeNovoSegment]:
        return [s for s in self.segments if isinstance(s, DeNovoSegment)]

    def designed_position_to_native(self, designed_pos: int) -> Optional[int]:
        """For a 1-based designed-frame position on this chain, return
        the corresponding native residue number, or None if the position
        falls in a de novo segment.

        Walks the segments in order, tracking the cumulative designed
        position; when the running cursor crosses ``designed_pos``, the
        position is in the current segment.
        """
        if designed_pos < 1:
            return None
        cursor = 0
        for seg in self.segments:
            seg_start = cursor + 1
            seg_end = cursor + seg.length
            if seg_start <= designed_pos <= seg_end:
                if isinstance(seg, FixedSegment):
                    offset = designed_pos - seg_start
                    return seg.start + offset
                return None  # de novo — no native counterpart
            cursor = seg_end
        return None

    def native_position_to_designed(self, native_pos: int) -> Optional[int]:
        """For a native residue number on this chain, return the
        corresponding 1-based designed-frame position, or None if the
        native position is not anchored by any fixed segment (i.e. it
        falls in or between the de novo / non-anchored regions).
        """
        cursor = 0
        for seg in self.segments:
            if isinstance(seg, FixedSegment):
                if seg.start <= native_pos <= seg.end:
                    offset = native_pos - seg.start
                    return cursor + 1 + offset
            cursor += seg.length
        return None

    def is_design_region_position(self, designed_pos: int) -> bool:
        """True iff designed_pos is in a de novo segment on this chain."""
        cursor = 0
        for seg in self.segments:
            seg_start = cursor + 1
            seg_end = cursor + seg.length
            if seg_start <= designed_pos <= seg_end:
                return isinstance(seg, DeNovoSegment)
            cursor = seg_end
        return False

    def design_region_positions(self) -> List[int]:
        """The 1-based designed-frame positions that fall in de novo
        segments on this chain."""
        out: List[int] = []
        cursor = 0
        for seg in self.segments:
            if isinstance(seg, DeNovoSegment):
                out.extend(range(cursor + 1, cursor + seg.length + 1))
            cursor += seg.length
        return out

    def fixed_anchor_positions(self) -> List[int]:
        """The 1-based designed-frame positions that DO have native
        counterparts (i.e. fall in fixed segments)."""
        out: List[int] = []
        cursor = 0
        for seg in self.segments:
            if isinstance(seg, FixedSegment):
                out.extend(range(cursor + 1, cursor + seg.length + 1))
            cursor += seg.length
        return out


# ═════════════════════════════════════════════════════════════════════════
# ContigSpec
# ═════════════════════════════════════════════════════════════════════════


# Bare chain identifier (e.g. ` B` after a space) — entire chain
# included with no design regions.  Used for the effector chain
# typically.  Encoded as a ContigChain with one "PassthroughSegment"
# (modelled as a FixedSegment of unknown length but with placeholder
# start/end pending PDB lookup).  Most consumers don't need to know
# the actual residue range of a passthrough chain at the spec level —
# RFDiffusion expands it from the input PDB at run time.  For Tier 0
# the simplest model: a passthrough chain is a ContigChain with an
# empty segments tuple; downstream code that needs the actual range
# reads it from the input PDB directly.

# Match a chain segment block (slash-separated) for one chain.
_SEG_PATTERN = re.compile(r"([A-Za-z]\d+(?:-\d+)?|\d+(?:-\d+)?)")


@dataclass(frozen=True)
class ContigSpec:
    """A resolved RFDiffusion contig: per-chain structured segments.

    Construct via :meth:`from_resolved_string`.  Direct construction is
    supported for tests and for programmatic generation.
    """
    chains: Tuple  # Tuple[ContigChain, ...]
    contig_string: str

    def __post_init__(self) -> None:
        if not isinstance(self.chains, tuple):
            object.__setattr__(self, "chains", tuple(self.chains))

    # ── Construction ─────────────────────────────────────────────────

    @classmethod
    def from_resolved_string(cls, s: str) -> "ContigSpec":
        """Parse a resolved contig string (specific de novo lengths).

        Format (see memory project_contig_string_format.md):
        - Space-separated chain blocks.
        - Within a block, segments are slash-separated.
        - Fixed segments are chain-prefixed native residue ranges:
          ``A1-10`` = chain A residues 1..10.
        - De novo segments are bare length integers: ``5`` = a de novo
          region of length 5.  Ranges (``5-7``) are NOT accepted here —
          ContigSpec is the resolved form.  Use the upstream
          length-resolver if you have an unresolved contig.
        - A bare chain identifier (e.g. ``B`` alone) = entire chain
          included with no design regions; encoded as a ContigChain with
          no segments (length-unknown until PDB lookup).
        """
        raw = s.strip()
        if not raw:
            raise ValueError("ContigSpec.from_resolved_string: empty input")
        blocks = raw.split()
        chains: List[ContigChain] = []
        for block in blocks:
            chain = _parse_block(block)
            chains.append(chain)
        return cls(chains=tuple(chains), contig_string=raw)

    # ── Per-chain accessors ──────────────────────────────────────────

    def chain(self, chain_id: str) -> ContigChain:
        for c in self.chains:
            if c.chain_id == chain_id:
                return c
        raise KeyError(f"chain {chain_id!r} not in contig (have: {self.chain_ids})")

    @property
    def chain_ids(self) -> List[str]:
        return [c.chain_id for c in self.chains]

    def chain_length(self, chain_id: str) -> int:
        return self.chain(chain_id).total_length

    def fixed_segments(self, chain_id: str) -> List[FixedSegment]:
        return self.chain(chain_id).fixed_segments

    def denovo_segments(self, chain_id: str) -> List[DeNovoSegment]:
        return self.chain(chain_id).denovo_segments

    # ── Frame conversion (per Q7 — lives here, called by PositionSet) ──

    def designed_to_native(
        self, pos: int, chain: str
    ) -> Optional[int]:
        """Map a designed-frame position to its native counterpart, or
        None if the position falls in a de novo segment.
        """
        return self.chain(chain).designed_position_to_native(pos)

    def native_to_designed(
        self, pos: int, chain: str
    ) -> Optional[int]:
        """Map a native residue number to its designed-frame position,
        or None if the native position is not anchored by any fixed
        segment.
        """
        return self.chain(chain).native_position_to_designed(pos)

    def prediction_to_designed(self, pos: int, chain: str) -> int:
        """Boltz / AF3 emit per-chain 1-based numbering.  For the
        receptor chain the prediction frame matches the designed frame
        (the receptor PDB went straight into Boltz).  Documented as a
        no-op for explicitness; included in the API so callers don't
        special-case the prediction frame elsewhere."""
        return pos

    # ── Design region access ─────────────────────────────────────────

    def design_region_positions(self, chain: str) -> List[int]:
        return self.chain(chain).design_region_positions()

    def fixed_anchor_positions(self, chain: str) -> List[int]:
        return self.chain(chain).fixed_anchor_positions()

    def is_design_region_position(self, pos: int, chain: str) -> bool:
        return self.chain(chain).is_design_region_position(pos)


# ═════════════════════════════════════════════════════════════════════════
# Parsing helpers
# ═════════════════════════════════════════════════════════════════════════


def _parse_block(block: str) -> ContigChain:
    """Parse one space-separated chain block, e.g. ``A1-10/5/A15-20`` or
    ``B`` (bare chain).
    """
    block = block.strip()
    if not block:
        raise ValueError("empty chain block")

    # Bare chain identifier: a single letter.
    if len(block) == 1 and block.isalpha():
        return ContigChain(chain_id=block, segments=())

    # Multi-segment block.  Determine the chain ID from the first segment
    # that has a chain prefix.  All fixed segments in a block must share
    # the same chain.
    parts = block.split("/")
    segments: List = []
    chain_id: Optional[str] = None
    for raw_seg in parts:
        seg = raw_seg.strip()
        if not seg:
            continue
        if seg[0].isalpha():
            # Fixed segment: chain-prefixed range
            chain_letter = seg[0]
            rest = seg[1:]
            if not rest:
                raise ValueError(
                    f"fixed segment {seg!r} has chain prefix but no range"
                )
            if "-" in rest:
                lo_s, hi_s = rest.split("-", 1)
                try:
                    lo, hi = int(lo_s), int(hi_s)
                except ValueError:
                    raise ValueError(
                        f"fixed segment {seg!r}: residue range must be int-int"
                    )
            else:
                try:
                    lo = hi = int(rest)
                except ValueError:
                    raise ValueError(
                        f"fixed segment {seg!r}: residue ref must be an int"
                    )
            if chain_id is None:
                chain_id = chain_letter
            elif chain_letter != chain_id:
                raise ValueError(
                    f"block {block!r} mixes chains: {chain_id!r} and "
                    f"{chain_letter!r}.  Use space-separated blocks for "
                    f"different chains."
                )
            segments.append(FixedSegment(start=lo, end=hi))
        elif seg[0].isdigit():
            # De novo segment: must be a single resolved length, not a range
            if "-" in seg:
                raise ValueError(
                    f"de novo segment {seg!r}: ContigSpec accepts resolved "
                    f"contigs only (length must be a single integer, not a "
                    f"range).  Resolve the range upstream before constructing "
                    f"the ContigSpec."
                )
            try:
                length = int(seg)
            except ValueError:
                raise ValueError(
                    f"de novo segment {seg!r}: length must be an integer"
                )
            segments.append(DeNovoSegment(length=length))
        else:
            raise ValueError(
                f"segment {seg!r}: unrecognised first character"
            )

    if chain_id is None:
        # The block was entirely de novo segments with no chain anchor —
        # unusual.  Refuse rather than guess.
        raise ValueError(
            f"block {block!r} has no chain-prefixed fixed segment to "
            f"determine the chain identifier"
        )

    return ContigChain(chain_id=chain_id, segments=tuple(segments))
