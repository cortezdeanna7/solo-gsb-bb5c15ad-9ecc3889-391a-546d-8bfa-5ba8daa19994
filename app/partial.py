"""Exact fragment-length expectations for partial (probabilistic) digests.

Each enzyme cuts at each of its cut positions independently with its own
probability ``p``; a nick forms at a position when at least one enzyme cuts
there, so a position shared by several enzymes has nick probability
``1 - prod(1 - p_i)``.  Because every (enzyme, position) cut event is
independent, nick events at distinct positions are independent as well, so
the exact fragment-length distribution falls out of products over the
sorted positions in O(K^2) time (K = number of effective cut positions) --
no sampling and no enumeration of nick subsets.
"""

from __future__ import annotations

from dataclasses import dataclass

from .digest import enzyme_cut_positions


@dataclass(frozen=True)
class PartialEnzyme:
    """One restriction enzyme with a per-position cut probability."""

    name: str
    site: str
    cut_offset: int
    p: float


def effective_cut_positions(
    sequence: str, enzymes: list[PartialEnzyme], circular: bool
) -> tuple[list[int], list[float]]:
    """Sorted effective cut positions and their nick probabilities.

    Reuses the deterministic matching and de-duplication rules per enzyme,
    then merges enzymes sharing a position into the nick probability
    ``1 - prod(1 - p_i)``.  Positions where every enzyme has ``p == 0``
    never nick and are dropped (they are not effective cut positions).  On
    linear molecules the two ends (0 and ``len(sequence)``) are not random
    nicks and are excluded.
    """
    n = len(sequence)
    stay_intact: dict[int, float] = {}  # position -> prod(1 - p_i) so far
    for enzyme in enzymes:
        if enzyme.p == 0.0:
            continue
        for pos in enzyme_cut_positions(
            sequence, enzyme.site, enzyme.cut_offset, circular
        ):
            if not circular and (pos == 0 or pos == n):
                continue
            stay_intact[pos] = stay_intact.get(pos, 1.0) * (1.0 - enzyme.p)
    ordered = sorted(
        (pos, 1.0 - intact)
        for pos, intact in stay_intact.items()
        if intact < 1.0
    )
    return [pos for pos, _ in ordered], [q for _, q in ordered]


def linear_expectations(
    n: int, positions: list[int], probs: list[float]
) -> tuple[dict[int, float], float]:
    """Exact expected fragment counts by length for a linear molecule.

    ``positions`` (sorted, strictly inside (0, n)) and ``probs`` describe
    the independent nick events.  A fragment (a, b) appears exactly when a
    and b are nicked (or are the molecule ends) and every cut position
    between them stays intact; independence makes its probability a
    product, accumulated here over all O(K^2) intervals with running
    products (no division, so probabilities 0 and 1 stay exact).  Returns
    the per-length expectations and the probability of no nick at all.
    """
    expected: dict[int, float] = {}

    def add(length: int, prob: float) -> None:
        if prob:
            expected[length] = expected.get(length, 0.0) + prob

    k = len(positions)
    # Fragments whose left end is the molecule's left end (position 0).
    intact = 1.0
    for j in range(k):
        add(positions[j], intact * probs[j])
        intact *= 1.0 - probs[j]
    add(n, intact)  # no nick anywhere -> one full-length fragment
    no_cut = intact

    # Fragments whose left end is each cut position.
    for i in range(k):
        run = probs[i]
        for j in range(i + 1, k):
            add(positions[j] - positions[i], run * probs[j])
            run *= 1.0 - probs[j]
        add(n - positions[i], run)
    return expected, no_cut


def circular_expectations(
    n: int, positions: list[int], probs: list[float]
) -> tuple[dict[int, float], float]:
    """Exact expected fragment counts by length for a circular molecule.

    With two or more nicks the fragments are the arcs between consecutive
    nicked positions; with zero or exactly one nick the molecule yields a
    single full-length fragment.  The returned no-cut probability counts
    only the zero-nick case.  Runs in O(K^2) over the sorted positions.
    """
    expected: dict[int, float] = {}

    def add(length: int, prob: float) -> None:
        if prob:
            expected[length] = expected.get(length, 0.0) + prob

    k = len(positions)
    # Prefix/suffix products of the stay-intact probabilities.
    pref = [1.0] * (k + 1)
    for i in range(k):
        pref[i + 1] = pref[i] * (1.0 - probs[i])
    suff = [1.0] * (k + 1)
    for i in range(k - 1, -1, -1):
        suff[i] = suff[i + 1] * (1.0 - probs[i])

    no_cut = pref[k]
    exactly_one = sum(probs[i] * pref[i] * suff[i + 1] for i in range(k))
    add(n, no_cut + exactly_one)

    # Arcs: for each nicked start i, extend the nicked end j forward.
    for i in range(k):
        run = probs[i]
        for step in range(1, k):
            j = (i + step) % k
            add((positions[j] - positions[i]) % n, run * probs[j])
            run *= 1.0 - probs[j]
    return expected, no_cut
