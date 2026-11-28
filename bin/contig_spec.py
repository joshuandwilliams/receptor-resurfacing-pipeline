"""
contig_spec.py
--------------
ContigSpec — RFDiffusion contig strings parsed into structured per-chain
segments.  Canonical input form: ``A1-10/5/A15-20 B``.  Chain blocks are
space-separated, segments within a block are slash-separated, and commas
are never used.

This module imports none of the other type modules in bin/, so it can be
imported on its own.

The resolved form is the main case.  Once RFDiffusion has sampled a
specific length for each de novo segment, the contig is a concrete
layout.  The constraint form carries length ranges like ``5-7``.
``from_string`` accepts both forms.  ``from_resolved_string`` rejects
ranges.

The three older contig parsers are now thin adapters that delegate here:

  - bin/contig_utils.py:parse_block_segments / resolve_contigs /
    parse_design_region
  - bin/derive_input_design_region.py:_parse_contigs
  - bin/pipeline_correct_sequences.py:parse_contig_segments
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import List, Optional, Tuple

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
class BreakSegment:
    """A chain-break marker (literal ``0`` token in RFDiffusion contig
    grammar).  Carries no length and no positions.  It only tells
    RFDiffusion that the flanking residues are not contiguous.
    Position-math methods on ContigChain ignore breaks.
    """

    def __post_init__(self) -> None:
        return None


@dataclass(frozen=True)
class PassthroughSegment:
    """A bare-chain-letter segment whose residue range cannot be
    resolved without a reference PDB.  Carried verbatim through the
    contig pipeline; RFDiffusion expands it from the input PDB at run
    time.  Like BreakSegment, has no length at the spec level.
    """
    chain: str

    def __post_init__(self) -> None:
        if not isinstance(self.chain, str) or not self.chain:
            raise ValueError("PassthroughSegment.chain must be non-empty")


@dataclass(frozen=True)
class DeNovoSegment:
    """A de novo (design) segment.  Carries a length RANGE (min_len,
    max_len) for both the constraint form (``5-7``) and the resolved
    form (``5`` is stored as ``min_len == max_len == 5``).  One type
    covers both, so there is no separate resolved segment class.

    Most downstream code only sees resolved segments, because
    RFDiffusion has already sampled a concrete length.  ``length``
    returns that integer when ``is_resolved`` is True and raises
    otherwise.
    """
    min_len: int
    max_len: int

    def __post_init__(self) -> None:
        if self.min_len < 0:
            raise ValueError(f"DeNovoSegment min_len={self.min_len} must be >= 0")
        if self.max_len < self.min_len:
            raise ValueError(
                f"DeNovoSegment max_len={self.max_len} < min_len={self.min_len}"
            )

    @property
    def is_resolved(self) -> bool:
        """True iff the segment has a single concrete length (min == max)."""
        return self.min_len == self.max_len

    @property
    def length(self) -> int:
        """Resolved length.  Raises ValueError if the segment is still a
        range (callers wanting the bounds should read min_len/max_len)."""
        if not self.is_resolved:
            raise ValueError(
                f"DeNovoSegment is unresolved ({self.min_len}-{self.max_len}); "
                f"no single length.  Use min_len / max_len explicitly."
            )
        return self.min_len

    @classmethod
    def resolved(cls, length: int) -> "DeNovoSegment":
        """Convenience: a fixed-length de novo segment.  Equivalent to
        ``DeNovoSegment(length, length)``."""
        return cls(min_len=length, max_len=length)

    def to_string(self) -> str:
        """Serialize: ``5`` if resolved, else ``5-7``."""
        if self.is_resolved:
            return str(self.min_len)
        return f"{self.min_len}-{self.max_len}"


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
    def is_resolved(self) -> bool:
        """True iff every de novo segment in this chain has a single
        concrete length AND the chain carries no Passthrough segments
        (whose length is PDB-dependent and unknown at spec level)."""
        for s in self.segments:
            if isinstance(s, DeNovoSegment) and not s.is_resolved:
                return False
            if isinstance(s, PassthroughSegment):
                return False
        return True

    @property
    def _length_bearing(self) -> List:
        """Length-bearing segments only — Fixed + DeNovo.  Break and
        Passthrough segments are excluded.  Position-math methods
        iterate this list rather than ``self.segments``."""
        return [
            s for s in self.segments
            if isinstance(s, (FixedSegment, DeNovoSegment))
        ]

    @property
    def total_length(self) -> int:
        """Total residues in the resolved chain (sum of fixed and de
        novo lengths).  Raises ValueError if any de novo segment is
        still a range.  For the constraint form use min_total_length /
        max_total_length instead."""
        if not self.is_resolved:
            raise ValueError(
                "ContigChain is unresolved; total_length is undefined.  "
                "Resolve de novo segment lengths or use min_total_length / "
                "max_total_length."
            )
        return sum(s.length for s in self._length_bearing)

    @property
    def min_total_length(self) -> int:
        """Minimum total residues across the chain (constraint form).
        Passthrough segments contribute 0 (unknown length until PDB)."""
        return sum(
            s.length if isinstance(s, FixedSegment) else s.min_len
            for s in self._length_bearing
        )

    @property
    def max_total_length(self) -> int:
        """Maximum total residues across the chain (constraint form).
        Passthrough segments contribute 0."""
        return sum(
            s.length if isinstance(s, FixedSegment) else s.max_len
            for s in self._length_bearing
        )

    @property
    def fixed_segments(self) -> List[FixedSegment]:
        return [s for s in self.segments if isinstance(s, FixedSegment)]

    @property
    def denovo_segments(self) -> List[DeNovoSegment]:
        return [s for s in self.segments if isinstance(s, DeNovoSegment)]

    @property
    def has_breaks(self) -> bool:
        return any(isinstance(s, BreakSegment) for s in self.segments)

    @property
    def passthrough_segments(self) -> List[PassthroughSegment]:
        return [s for s in self.segments if isinstance(s, PassthroughSegment)]

    def designed_position_to_native(self, designed_pos: int) -> Optional[int]:
        """For a 1-based designed-frame position on this chain, return
        the corresponding native residue number, or None if the position
        falls in a de novo segment.
        """
        if designed_pos < 1:
            return None
        cursor = 0
        for seg in self._length_bearing:
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
        corresponding 1-based designed-frame position, or None if no
        fixed segment anchors that native position.
        """
        cursor = 0
        for seg in self._length_bearing:
            if isinstance(seg, FixedSegment):
                if seg.start <= native_pos <= seg.end:
                    offset = native_pos - seg.start
                    return cursor + 1 + offset
            cursor += seg.length
        return None

    def is_design_region_position(self, designed_pos: int) -> bool:
        """True iff designed_pos is in a de novo segment on this chain."""
        cursor = 0
        for seg in self._length_bearing:
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
        for seg in self._length_bearing:
            if isinstance(seg, DeNovoSegment):
                out.extend(range(cursor + 1, cursor + seg.length + 1))
            cursor += seg.length
        return out

    def fixed_anchor_positions(self) -> List[int]:
        """The 1-based designed-frame positions that DO have native
        counterparts (i.e. fall in fixed segments)."""
        out: List[int] = []
        cursor = 0
        for seg in self._length_bearing:
            if isinstance(seg, FixedSegment):
                out.extend(range(cursor + 1, cursor + seg.length + 1))
            cursor += seg.length
        return out


# ═════════════════════════════════════════════════════════════════════════
# ContigSpec
# ═════════════════════════════════════════════════════════════════════════


# A block that is a bare chain identifier (e.g. ``B`` in ``A1-10/5 B``)
# means the entire chain, with no design regions.  Typically the
# effector chain.  _parse_block encodes it as a ContigChain with an
# empty segments tuple, because the residue range is unknown without
# the input PDB.  RFDiffusion expands it from the PDB at run time, and
# downstream code that needs the real range reads the PDB directly.
# A bare chain letter appearing as one segment *inside* a block becomes
# a PassthroughSegment instead.

# Segment token: chain-prefixed residue ref, or a bare length / range.
# Currently unused — _parse_block splits on "/" and inspects each token.
_SEG_PATTERN = re.compile(r"([A-Za-z]\d+(?:-\d+)?|\d+(?:-\d+)?)")


@dataclass(frozen=True)
class ContigSpec:
    """An RFDiffusion contig as per-chain structured segments.

    Construct via :meth:`from_string`, or :meth:`from_resolved_string`
    to reject unresolved de novo ranges.  Direct construction is
    supported for tests and for programmatic generation.
    """
    chains: Tuple  # Tuple[ContigChain, ...]
    contig_string: str

    def __post_init__(self) -> None:
        if not isinstance(self.chains, tuple):
            object.__setattr__(self, "chains", tuple(self.chains))

    # ── Construction ─────────────────────────────────────────────────

    @classmethod
    def from_string(cls, s: str) -> "ContigSpec":
        """Parse a contig string in either form.

        Format:
        - Space-separated chain blocks.
        - Within a block, segments are slash-separated.
        - Fixed segments are chain-prefixed native residue ranges:
          ``A1-10`` = chain A residues 1..10.
        - De novo segments are bare lengths.  ``5`` is resolved
          (min_len == max_len == 5).  ``5-7`` is constraint form
          (RFDiffusion samples a length in [5, 7]).  Both round-trip
          through DeNovoSegment(min_len, max_len).
        - A bare chain identifier (e.g. ``B`` alone) means the entire
          chain with no design regions.  It is encoded as a ContigChain
          with no segments, since the length is unknown until PDB
          lookup.
        - ``0`` is RFDiffusion's chain-break marker.
        """
        raw = s.strip()
        if not raw:
            raise ValueError("ContigSpec.from_string: empty input")
        blocks = raw.split()
        chains: List[ContigChain] = []
        for block in blocks:
            chain = _parse_block(block)
            chains.append(chain)
        return cls(chains=tuple(chains), contig_string=raw)

    @classmethod
    def from_resolved_string(cls, s: str) -> "ContigSpec":
        """Parse a contig string and assert every de novo segment is
        resolved.  Raises ValueError if any segment is a range.  Use
        ``from_string`` to tolerate the constraint form."""
        spec = cls.from_string(s)
        if not spec.is_resolved:
            raise ValueError(
                f"ContigSpec.from_resolved_string({s!r}): contains "
                f"unresolved de novo segment(s).  Resolve lengths upstream "
                f"or call from_string instead."
            )
        return spec

    @property
    def is_resolved(self) -> bool:
        """True iff every de novo segment across every chain is a single
        concrete length (min_len == max_len)."""
        return all(c.is_resolved for c in self.chains)

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

    # ── Frame conversion (lives here, called by PositionSet) ─────────

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
        receptor chain the prediction frame already matches the designed
        frame, because the receptor PDB went straight into Boltz.  This
        is a deliberate no-op, kept in the API so callers need no
        special case for the prediction frame."""
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
        if seg == "0":
            # RFDiffusion chain-break marker.  Carries no length.
            segments.append(BreakSegment())
            continue
        if seg[0].isalpha():
            chain_letter = seg[0]
            rest = seg[1:]
            if not rest:
                # Bare chain letter within a block — passthrough.
                if chain_id is None:
                    chain_id = chain_letter
                segments.append(PassthroughSegment(chain=chain_letter))
                continue
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
            # De novo segment: bare length `5` (resolved) or range `5-7`
            # (constraint).  Both are stored as (min_len, max_len) with
            # min == max for the resolved case.
            if "-" in seg:
                lo_s, hi_s = seg.split("-", 1)
                try:
                    lo_n, hi_n = int(lo_s), int(hi_s)
                except ValueError:
                    raise ValueError(
                        f"de novo segment {seg!r}: range must be int-int"
                    )
                segments.append(DeNovoSegment(min_len=lo_n, max_len=hi_n))
            else:
                try:
                    length = int(seg)
                except ValueError:
                    raise ValueError(
                        f"de novo segment {seg!r}: length must be an integer"
                    )
                segments.append(DeNovoSegment(min_len=length, max_len=length))
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
