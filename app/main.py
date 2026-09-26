"""FastAPI application exposing the restriction-digest simulation."""

from __future__ import annotations

from fastapi import FastAPI

from .digest import digest, enzyme_cut_positions
from .panel import PanelEnzyme, PanelSample, design_panel
from .schemas import (
    DesignPanelRequest,
    DesignPanelResponse,
    DigestRequest,
    DigestResponse,
    EnzymeCuts,
    FragmentOut,
    PairDistinction,
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
        "Also designs minimal enzyme panels that distinguish DNA samples."
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
