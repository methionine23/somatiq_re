"""Locus catalog: load, validate, and resolve per-sample ploidy."""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class Interruption:
    motifs: tuple[str, ...] = ()
    position: str | None = None
    pure_tract: str | None = None
    notes: str = ""


@dataclass(frozen=True)
class HipstrRegion:
    start: int  # 0-based
    end: int  # 1-based inclusive end (BED-style half-open with 0-based start)
    period: int
    ref_copies: int
    name: str


@dataclass(frozen=True)
class Locus:
    locus_id: str
    disease: str
    chrom: str
    ref_start: int  # 1-based inclusive
    ref_end: int  # 1-based inclusive
    ref_motif: str
    disease_motif: str
    period: int
    ploidy: str  # "diploid" | "haploid" | "sex_dependent"
    eh_variant_id: str
    hipstr_region: HipstrRegion
    interruption: Interruption = field(default_factory=Interruption)

    def resolve_ploidy(self, sex: str | None) -> int:
        """Return effective copy number of chromosomes at this locus.

        sex is "M"/"F"/None. sex_dependent (AR on chrX) is 1 for males, 2 for
        females. Explicit haploid/diploid override sex.
        """
        if self.ploidy == "haploid":
            return 1
        if self.ploidy == "diploid":
            return 2
        if self.ploidy == "sex_dependent":
            if sex is None:
                raise ValueError(
                    f"locus {self.locus_id} is sex_dependent; --sex is required"
                )
            s = sex.strip().upper()
            if s in ("M", "MALE", "1"):
                return 1
            if s in ("F", "FEMALE", "2"):
                return 2
            raise ValueError(f"unrecognized sex {sex!r} (use M/F)")
        raise ValueError(f"unknown ploidy {self.ploidy!r} for {self.locus_id}")


@dataclass(frozen=True)
class Catalog:
    genome_build: str
    chrom_style: str
    loci: tuple[Locus, ...]

    def __iter__(self):
        return iter(self.loci)

    def __len__(self):
        return len(self.loci)

    def get(self, locus_id: str) -> Locus:
        for locus in self.loci:
            if locus.locus_id == locus_id:
                return locus
        raise KeyError(locus_id)

    def with_chrom_style(self, style: str) -> "Catalog":
        """Return a catalog whose chrom names match the reference style.

        style is "chr" (chr18) or "nochr" (18). No-op if already matching.
        """
        if style not in ("chr", "nochr"):
            raise ValueError("chrom style must be 'chr' or 'nochr'")

        def conv(chrom: str) -> str:
            has = chrom.startswith("chr")
            if style == "chr":
                return chrom if has else "chr" + chrom
            return chrom[3:] if has else chrom

        new = tuple(
            Locus(**{**locus.__dict__, "chrom": conv(locus.chrom)}) for locus in self.loci
        )
        return Catalog(self.genome_build, style, new)


_REQUIRED_LOCUS_KEYS = {
    "locus_id",
    "disease",
    "chrom",
    "ref_region",
    "ref_motif",
    "disease_motif",
    "period",
    "ploidy",
    "eh_variant_id",
    "hipstr_region",
}


def _parse_locus(d: dict[str, Any]) -> Locus:
    missing = _REQUIRED_LOCUS_KEYS - d.keys()
    if missing:
        raise ValueError(f"locus {d.get('locus_id', '?')} missing keys: {sorted(missing)}")
    rr = d["ref_region"]
    if rr.get("base", 1) != 1:
        raise ValueError("ref_region.base must be 1 (1-based inclusive)")
    hr = d["hipstr_region"]
    itr = d.get("interruption") or {}
    return Locus(
        locus_id=d["locus_id"],
        disease=d["disease"],
        chrom=d["chrom"],
        ref_start=int(rr["start"]),
        ref_end=int(rr["end"]),
        ref_motif=d["ref_motif"].upper(),
        disease_motif=d["disease_motif"].upper(),
        period=int(d["period"]),
        ploidy=d["ploidy"],
        eh_variant_id=d["eh_variant_id"],
        hipstr_region=HipstrRegion(
            start=int(hr["start"]),
            end=int(hr["end"]),
            period=int(hr["period"]),
            ref_copies=int(hr["ref_copies"]),
            name=hr["name"],
        ),
        interruption=Interruption(
            motifs=tuple(m.upper() for m in itr.get("motifs", [])),
            position=itr.get("position"),
            pure_tract=(itr.get("pure_tract") or d["ref_motif"]).upper(),
            notes=itr.get("notes", ""),
        ),
    )


def load_catalog(path: str | Path) -> Catalog:
    data = json.loads(Path(path).read_text())
    loci = tuple(_parse_locus(x) for x in data["loci"])
    if not loci:
        raise ValueError("catalog has no loci")
    ids = [locus.locus_id for locus in loci]
    if len(set(ids)) != len(ids):
        raise ValueError("duplicate locus_id in catalog")
    cat = Catalog(
        genome_build=data.get("genome_build", "GRCh38"),
        chrom_style=data.get("chrom_style", "chr"),
        loci=loci,
    )
    validate_catalog(cat)
    return cat


def validate_catalog(cat: Catalog) -> list[str]:
    """Return a list of non-fatal warnings; raise on fatal inconsistencies."""
    warnings: list[str] = []
    for locus in cat.loci:
        if len(locus.ref_motif) != locus.period:
            raise ValueError(
                f"{locus.locus_id}: ref_motif len != period ({locus.ref_motif}/{locus.period})"
            )
        if locus.ref_end < locus.ref_start:
            raise ValueError(f"{locus.locus_id}: ref_region end < start")
        region_len = locus.ref_end - locus.ref_start + 1
        derived = region_len / locus.period
        if abs(derived - locus.hipstr_region.ref_copies) > 1.0:
            warnings.append(
                f"{locus.locus_id}: hipstr ref_copies={locus.hipstr_region.ref_copies} "
                f"differs from region_len/period={derived:.2f}; validate region"
            )
    return warnings
