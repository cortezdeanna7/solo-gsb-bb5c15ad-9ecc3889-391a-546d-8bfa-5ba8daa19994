"""Core restriction-digest logic (no web framework dependencies)."""

from __future__ import annotations

from dataclasses import dataclass

# IUPAC degenerate base codes -> the concrete bases they match.
IUPAC_MATCH: dict[str, frozenset[str]] = {
    "A": frozenset("A"),
    "C": frozenset("C"),
    "G": frozenset("G"),
    "T": frozenset("T"),
    "R": frozenset("AG"),
    "Y": frozenset("CT"),
    "S": frozenset("GC"),
    "W": frozenset("AT"),
    "K": frozenset("GT"),
    "M": frozenset("AC"),
    "B": frozenset("CGT"),
    "D": frozenset("AGT"),
    "H": frozenset("ACT"),
    "V": frozenset("ACG"),
    "N": frozenset("ACGT"),
}

# Complement of each IUPAC code (used to build the reverse complement site).
IUPAC_COMPLEMENT: dict[str, str] = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "R": "Y",
    "Y": "R",
    "S": "S",
    "W": "W",
    "K": "M",
    "M": "K",
    "B": "V",
    "V": "B",
    "D": "H",
    "H": "D",
    "N": "N",
}


@dataclass(frozen=True)
class Fragment:
    """A digest product: zero-based start on the original sequence, length,
    and the bases in the original (5'->3') direction."""

    start: int
    length: int
    sequence: str


def reverse_complement(site: str) -> str:
    """Reverse complement of an IUPAC recognition string."""
    return "".join(IUPAC_COMPLEMENT[base] for base in reversed(site))


def _match_starts(sequence: str, pattern: str, circular: bool) -> list[int]:
    """All start positions where ``pattern`` matches ``sequence``.

    Overlapping matches are kept.  For circular molecules a match may span
    the origin, so starts run over every position and the tail of the
    sequence is prepended conceptually by scanning an extended string.
    """
    n = len(sequence)
    m = len(pattern)
    if circular:
        haystack = sequence + sequence[: m - 1]
        limit = n
    else:
        haystack = sequence
        limit = n - m + 1

    starts: list[int] = []
    for s in range(limit):
        for i, code in enumerate(pattern):
            if haystack[s + i] not in IUPAC_MATCH[code]:
                break
        else:
            starts.append(s)
    return starts


def enzyme_cut_positions(
    sequence: str, site: str, cut_offset: int, circular: bool
) -> list[int]:
    """Sorted, de-duplicated zero-based cut positions for one enzyme.

    The recognition string and its reverse complement are both searched:
    a forward-strand match at start ``s`` cuts at ``s + k``; a
    reverse-complement match cuts at ``s + len(site) - k``.  On circular
    molecules positions are taken modulo the sequence length.
    """
    n = len(sequence)
    m = len(site)
    cuts: set[int] = set()

    for s in _match_starts(sequence, site, circular):
        cuts.add((s + cut_offset) % n if circular else s + cut_offset)

    rc_site = reverse_complement(site)
    for s in _match_starts(sequence, rc_site, circular):
        cuts.add((s + m - cut_offset) % n if circular else s + m - cut_offset)

    return sorted(cuts)


def digest(sequence: str, cut_positions: list[int], circular: bool) -> list[Fragment]:
    """Cut ``sequence`` at all ``cut_positions`` simultaneously.

    Fragments are returned sorted by start position, each carrying its
    bases in the original sequence direction.  With no cut positions the
    original sequence is returned as a single fragment; empty fragments
    are never produced.
    """
    n = len(sequence)
    cuts = sorted(set(cut_positions))

    if circular:
        cuts = sorted({c % n for c in cuts})
        if not cuts:
            return [Fragment(start=0, length=n, sequence=sequence)]
        fragments: list[Fragment] = []
        for i, c in enumerate(cuts):
            nxt = cuts[(i + 1) % len(cuts)]
            length = (nxt - c) % n
            if length == 0:
                # A single cut linearises the circle at that point.
                length = n
            if nxt <= c:  # wraps around the origin
                bases = sequence[c:] + sequence[:nxt]
            else:
                bases = sequence[c:nxt]
            fragments.append(Fragment(start=c, length=length, sequence=bases))
        return fragments

    boundaries = sorted({0, n, *cuts})
    return [
        Fragment(start=a, length=b - a, sequence=sequence[a:b])
        for a, b in zip(boundaries, boundaries[1:])
        if b > a
    ]
