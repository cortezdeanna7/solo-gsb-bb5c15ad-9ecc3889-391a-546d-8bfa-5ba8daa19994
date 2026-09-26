"""Pydantic request/response models for the digest API."""

from __future__ import annotations

import re
from typing import Literal

from pydantic import BaseModel, Field, field_validator, model_validator

from .digest import IUPAC_MATCH

_SEQUENCE_RE = re.compile(r"^[ACGT]+$")
_SITE_RE = re.compile(r"^[ACGTRYSWKMBDHVN]+$")


class EnzymeRule(BaseModel):
    """One restriction enzyme: unique name, IUPAC recognition string and the
    in-string cut offset k (0 <= k <= len(site))."""

    name: str = Field(min_length=1, examples=["EcoRI"])
    site: str = Field(min_length=1, examples=["GAATTC"])
    cut_offset: int = Field(ge=0, examples=[1])

    @field_validator("site")
    @classmethod
    def site_must_be_iupac(cls, value: str) -> str:
        value = value.upper()
        if not _SITE_RE.fullmatch(value):
            raise ValueError(
                "site must contain only IUPAC codes: " + "".join(sorted(IUPAC_MATCH))
            )
        return value

    @model_validator(mode="after")
    def offset_within_site(self) -> "EnzymeRule":
        if self.cut_offset > len(self.site):
            raise ValueError(
                f"cut_offset ({self.cut_offset}) must be <= len(site) ({len(self.site)})"
            )
        return self


class DigestRequest(BaseModel):
    sequence: str = Field(min_length=1, examples=["ACGTGAATTCGGGAATTCAA"])
    topology: Literal["linear", "circular"] = "linear"
    enzymes: list[EnzymeRule] = Field(default_factory=list)

    @field_validator("sequence")
    @classmethod
    def sequence_must_be_acgt(cls, value: str) -> str:
        value = value.upper()
        if not _SEQUENCE_RE.fullmatch(value):
            raise ValueError("sequence must contain only A, C, G, T")
        return value

    @model_validator(mode="after")
    def check_enzymes(self) -> "DigestRequest":
        names = [e.name for e in self.enzymes]
        if len(names) != len(set(names)):
            raise ValueError("enzyme names must be unique")
        for e in self.enzymes:
            if len(e.site) > len(self.sequence):
                raise ValueError(
                    f"recognition site of enzyme {e.name!r} "
                    f"({len(e.site)} bp) is longer than the sequence "
                    f"({len(self.sequence)} bp)"
                )
        return self


class PartialEnzymeRule(EnzymeRule):
    """Enzyme rule plus the probability p (0..1) of cutting at each site."""

    p: float = Field(ge=0.0, le=1.0, examples=[0.5])


class PartialDigestRequest(DigestRequest):
    """Partial-digest input: sequence (<= 2000 bp), topology and enzymes
    each carrying a cut probability p.  Inherits the /digest validations
    (ACGT sequence, unique enzyme names, site not longer than sequence)."""

    sequence: str = Field(
        min_length=1, max_length=2000, examples=["ACGTGAATTCGGGAATTCAA"]
    )
    enzymes: list[PartialEnzymeRule] = Field(default_factory=list)


class PartialFragmentExpectation(BaseModel):
    """Expected number of fragments of one length (zero entries omitted)."""

    length: int
    expected_count: float


class PartialDigestResponse(BaseModel):
    fragments: list[PartialFragmentExpectation]
    no_cut_probability: float


class EnzymeCuts(BaseModel):
    name: str
    cuts: list[int]


class FragmentOut(BaseModel):
    start: int
    length: int
    sequence: str


class DigestResponse(BaseModel):
    cuts: list[EnzymeCuts]
    fragments: list[FragmentOut]


class PanelSampleIn(BaseModel):
    """One DNA sample: unique name, ACGT sequence and topology."""

    name: str = Field(min_length=1, examples=["wild-type"])
    sequence: str = Field(min_length=1, examples=["AAGAATTCTT"])
    topology: Literal["linear", "circular"] = "linear"

    @field_validator("sequence")
    @classmethod
    def sequence_must_be_acgt(cls, value: str) -> str:
        value = value.upper()
        if not _SEQUENCE_RE.fullmatch(value):
            raise ValueError("sequence must contain only A, C, G, T")
        return value


class DesignPanelRequest(BaseModel):
    """Panel design input: 2-12 uniquely named samples, 1-12 candidate
    enzymes (same rules as /digest) and a non-negative band tolerance."""

    samples: list[PanelSampleIn] = Field(min_length=2, max_length=12)
    enzymes: list[EnzymeRule] = Field(min_length=1, max_length=12)
    tolerance_bp: int = Field(ge=0, examples=[0])

    @model_validator(mode="after")
    def check_panel(self) -> "DesignPanelRequest":
        names = [s.name for s in self.samples]
        if len(names) != len(set(names)):
            raise ValueError("sample names must be unique")
        enzyme_names = [e.name for e in self.enzymes]
        if len(enzyme_names) != len(set(enzyme_names)):
            raise ValueError("enzyme names must be unique")
        for s in self.samples:
            for e in self.enzymes:
                if len(e.site) > len(s.sequence):
                    raise ValueError(
                        f"recognition site of enzyme {e.name!r} "
                        f"({len(e.site)} bp) is longer than the sequence "
                        f"of sample {s.name!r} ({len(s.sequence)} bp)"
                    )
        return self


class SampleLaneBands(BaseModel):
    """One sample's bands in each selected lane (enzyme)."""

    sample: str
    bands: dict[str, list[int]]


class SamplePair(BaseModel):
    samples: list[str] = Field(min_length=2, max_length=2)


class PairDistinction(SamplePair):
    """The selected enzymes able to tell this sample pair apart."""

    enzymes: list[str]


class DesignPanelResponse(BaseModel):
    enzymes: list[str]
    lanes: list[SampleLaneBands]
    pairs: list[PairDistinction]
    undistinguishable_pairs: list[SamplePair]
