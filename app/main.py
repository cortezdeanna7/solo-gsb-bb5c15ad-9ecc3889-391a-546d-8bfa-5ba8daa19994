"""FastAPI application exposing the restriction-digest simulation."""

from __future__ import annotations

from fastapi import FastAPI

from .digest import digest, enzyme_cut_positions
from .panel import PanelEnzyme, PanelSample, design_panel
from .partial import (
    PartialEnzyme,
    circular_expectations,
    effective_cut_positions,
    linear_expectations,
)
from .schemas import (
    DesignPanelRequest,
    DesignPanelResponse,
    DigestRequest,
    DigestResponse,
    EnzymeCuts,
    FragmentOut,
    PairDistinction,
    PartialDigestRequest,
    PartialDigestResponse,
    PartialFragmentExpectation,
    SampleLaneBands,
    SamplePair,
)

app = FastAPI(
    title="DNA Restriction Digest Simulator",
    version="1.0.0",
    description=(
        "Simulates restriction-enzyme digestion of linear or circular DNA. "
        "Recognition sites (IUPAC) and their reverse complements are searched "
        "with overlaps kept; all enzymes cut simultaneously. "
        "Also computes exact fragment-length expectations for partial "
        "(probabilistic) digests, and designs minimal enzyme panels that "
        "distinguish DNA samples."
    ),
)


@app.get("/health", tags=["meta"])
def health() -> dict[str, str]:
    return {"status": "ok"}


@app.post("/digest", response_model=DigestResponse, tags=["digest"])
def run_digest(request: DigestRequest) -> DigestResponse:
    circular = request.topology == "circular"

    per_enzyme: list[EnzymeCuts] = []
    all_cuts: list[int] = []
    for enzyme in request.enzymes:
        positions = enzyme_cut_positions(
            request.sequence, enzyme.site, enzyme.cut_offset, circular
        )
        per_enzyme.append(EnzymeCuts(name=enzyme.name, cuts=positions))
        all_cuts.extend(positions)

    fragments = digest(request.sequence, all_cuts, circular)

    return DigestResponse(
        cuts=per_enzyme,
        fragments=[FragmentOut(start=f.start, length=f.length, sequence=f.sequence) for f in fragments],
    )


@app.post("/partial-digest", response_model=PartialDigestResponse, tags=["partial-digest"])
def run_partial_digest(request: PartialDigestRequest) -> PartialDigestResponse:
    circular = request.topology == "circular"
    positions, probs = effective_cut_positions(
        request.sequence,
        [
            PartialEnzyme(name=e.name, site=e.site, cut_offset=e.cut_offset, p=e.p)
            for e in request.enzymes
        ],
        circular,
    )

    n = len(request.sequence)
    if circular:
        expected, no_cut = circular_expectations(n, positions, probs)
    else:
        expected, no_cut = linear_expectations(n, positions, probs)

    return PartialDigestResponse(
        fragments=[
            PartialFragmentExpectation(length=length, expected_count=count)
            for length, count in sorted(expected.items())
            if count != 0.0
        ],
        no_cut_probability=no_cut,
    )


@app.post("/design-panel", response_model=DesignPanelResponse, tags=["design-panel"])
def run_design_panel(request: DesignPanelRequest) -> DesignPanelResponse:
    result = design_panel(
        samples=[
            PanelSample(
                name=s.name,
                sequence=s.sequence,
                circular=s.topology == "circular",
            )
            for s in request.samples
        ],
        enzymes=[
            PanelEnzyme(name=e.name, site=e.site, cut_offset=e.cut_offset)
            for e in request.enzymes
        ],
        tolerance_bp=request.tolerance_bp,
    )

    return DesignPanelResponse(
        enzymes=result.enzymes,
        lanes=[
            SampleLaneBands(sample=name, bands=bands)
            for name, bands in result.lanes.items()
        ],
        pairs=[
            PairDistinction(samples=list(pair), enzymes=enzymes)
            for pair, enzymes in result.pair_enzymes
        ],
        undistinguishable_pairs=[
            SamplePair(samples=list(pair)) for pair in result.undistinguishable_pairs
        ],
    )
