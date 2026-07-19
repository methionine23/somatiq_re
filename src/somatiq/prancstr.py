"""prancSTR wrapper: command construction + .tab output parsing.

prancSTR (TRTools) is HipSTR-only. It reads the HipSTR VCF (MALLREADS + the
INFRAME_* stutter params) and writes '<out>.tab' with columns:
  sample chrom pos locus motif A B C f pval reads mosaic_support
  'stutter parameter u' 'stutter paramter d' 'stutter paramter rho'
  'quality factor' 'read depth'
(A,B = germline diploid alleles; C = mosaic allele; f = mosaic fraction).
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class PrancResult:
    sample: str
    chrom: str
    pos: int
    locus: str
    motif: str
    allele_a: float | None  # germline allele A
    allele_b: float | None  # germline allele B
    mosaic_c: float | None  # C: mosaic allele
    f: float | None  # mosaic fraction
    pval: float | None
    reads: int | None
    mosaic_support: int | None
    read_depth: int | None
    raw: dict[str, str]


def build_prancstr_command(
    vcf: str,
    out_prefix: str,
    *,
    vcftype: str = "hipstr",
    samples: list[str] | None = None,
    region: str | None = None,
    readfield: str = "MALLREADS",
    only_passing: bool = False,
    output_all: bool = False,
    extra: Iterable[str] = (),
) -> list[str]:
    """Construct the prancSTR command (verified against TRTools source)."""
    cmd = [
        "prancSTR",
        "--vcf",
        vcf,
        "--out",
        out_prefix,
        "--vcftype",
        vcftype,
        "--readfield",
        readfield,
    ]
    if samples:
        cmd += ["--samples", ",".join(samples)]
    if region:
        cmd += ["--region", region]
    if only_passing:
        cmd.append("--only-passing")
    if output_all:
        cmd.append("--output-all")
    cmd += list(extra)
    return cmd


# Header names in the prancSTR .tab (some carry the upstream spelling quirks).
_COL = {
    "sample": "sample",
    "chrom": "chrom",
    "pos": "pos",
    "locus": "locus",
    "motif": "motif",
    "a": "A",
    "b": "B",
    "c": "C",
    "f": "f",
    "pval": "pval",
    "reads": "reads",
    "mosaic_support": "mosaic_support",
    "read_depth": "read depth",
}


def _num(row: dict[str, str], key: str, cast):
    v = row.get(key)
    if v is None or v.strip() in (".", "", "NA", "nan"):
        return None
    try:
        return cast(float(v)) if cast is int else cast(v)
    except (ValueError, TypeError):
        return None


def parse_prancstr_tab(path: str | Path) -> list[PrancResult]:
    lines = Path(path).read_text().splitlines()
    return parse_prancstr_lines(lines)


def parse_prancstr_lines(lines: list[str]) -> list[PrancResult]:
    rows: list[PrancResult] = []
    header: list[str] | None = None
    for line in lines:
        if not line.strip():
            continue
        fields = line.rstrip("\n").split("\t")
        if header is None:
            header = fields
            continue
        row = dict(zip(header, fields))
        rows.append(
            PrancResult(
                sample=row.get(_COL["sample"], ""),
                chrom=row.get(_COL["chrom"], ""),
                pos=_num(row, _COL["pos"], int) or 0,
                locus=row.get(_COL["locus"], ""),
                motif=row.get(_COL["motif"], ""),
                allele_a=_num(row, _COL["a"], float),
                allele_b=_num(row, _COL["b"], float),
                mosaic_c=_num(row, _COL["c"], float),
                f=_num(row, _COL["f"], float),
                pval=_num(row, _COL["pval"], float),
                reads=_num(row, _COL["reads"], int),
                mosaic_support=_num(row, _COL["mosaic_support"], int),
                read_depth=_num(row, _COL["read_depth"], int),
                raw=row,
            )
        )
    return rows
