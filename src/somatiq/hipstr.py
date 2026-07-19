"""HipSTR wrapper: command construction + VCF parsing.

We run HipSTR only on the catalog regions (targeted; never a full-CRAM scan) so it
learns the per-locus stutter model prancSTR needs and emits MALLREADS. The parser is
stdlib-only and unit-tested; execution is a thin subprocess call.
"""
from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from .catalog import Catalog, Locus


@dataclass(frozen=True)
class HipstrCall:
    chrom: str
    pos: int
    ref_seq: str
    alt_seqs: tuple[str, ...]
    gt: tuple[int, ...]  # allele indices (0=ref), () if no-call
    gb_bp: tuple[int, ...]  # per-called-allele repeat bp (from GB), () if absent
    q: float | None
    dp: int | None
    dstutter: int | None
    dflankindel: int | None
    mallreads: str
    allreads: str
    inframe_up: float | None
    inframe_down: float | None
    inframe_pgeom: float | None

    def called_allele_seqs(self) -> list[str]:
        seqs = [self.ref_seq, *self.alt_seqs]
        out = []
        for i in self.gt:
            if 0 <= i < len(seqs):
                out.append(seqs[i])
        return out

    def allele_copies(self, period: int) -> list[float]:
        return [len(s) / period for s in self.called_allele_seqs()]

    @property
    def stutter_learned(self) -> bool:
        return self.inframe_pgeom is not None


def build_hipstr_command(
    bams: list[str],
    fasta: str,
    regions_bed: str,
    out_vcf: str,
    *,
    def_stutter_model: bool = True,
    viz_out: str | None = None,
    extra: Iterable[str] = (),
) -> list[str]:
    """Construct the HipSTR command line (flags confirmed vs installed version at
    build step 1, docs/P0_SPEC.md 1)."""
    cmd = [
        "HipSTR",
        "--bams",
        ",".join(bams),
        "--fasta",
        fasta,
        "--regions",
        regions_bed,
        "--str-vcf",
        out_vcf,
    ]
    if def_stutter_model:
        cmd.append("--def-stutter-model")
    if viz_out:
        cmd += ["--viz-out", viz_out]
    cmd += list(extra)
    return cmd


def write_hipstr_regions(cat: Catalog, path: str | Path) -> Path:
    """Write the HipSTR region file: chrom start end period ref_copies name."""
    path = Path(path)
    lines = []
    for locus in cat.loci:
        r = locus.hipstr_region
        lines.append(
            f"{locus.chrom}\t{r.start}\t{r.end}\t{r.period}\t{r.ref_copies}\t{r.name}"
        )
    path.write_text("\n".join(lines) + "\n")
    return path


def _open_text(path: str | Path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def _info_dict(info: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for field in info.split(";"):
        if "=" in field:
            k, _, v = field.partition("=")
            out[k] = v
        elif field:
            out[field] = "true"
    return out


def _to_float(x: str | None) -> float | None:
    if x is None or x in (".", ""):
        return None
    try:
        return float(x)
    except ValueError:
        return None


def _to_int(x: str | None) -> int | None:
    f = _to_float(x)
    return int(f) if f is not None else None


def _parse_gt(gt: str) -> tuple[int, ...]:
    if gt in (".", "./.", ".|.", ""):
        return ()
    sep = "|" if "|" in gt else "/"
    out = []
    for a in gt.split(sep):
        if a == ".":
            continue
        try:
            out.append(int(a))
        except ValueError:
            continue
    return tuple(out)


def parse_hipstr_record(line: str, sample_cols: list[str], sample: str) -> HipstrCall | None:
    """Parse one VCF data line for one sample into a HipstrCall (None if no data)."""
    f = line.rstrip("\n").split("\t")
    if len(f) < 10:
        return None
    chrom, pos, _id, ref, alt, _qual, _filt, info, fmt = f[:9]
    try:
        col = 9 + sample_cols.index(sample)
    except ValueError:
        raise KeyError(f"sample {sample!r} not in VCF")
    sample_field = f[col]
    keys = fmt.split(":")
    vals = sample_field.split(":")
    kv = dict(zip(keys, vals))
    info_d = _info_dict(info)

    gb_bp: tuple[int, ...] = ()
    if kv.get("GB", ".") not in (".", ""):
        gb_bp = tuple(
            v for v in (_to_int(x) for x in kv["GB"].split("|")) if v is not None
        )

    return HipstrCall(
        chrom=chrom,
        pos=int(pos),
        ref_seq=ref,
        alt_seqs=tuple(a for a in alt.split(",") if a != "."),
        gt=_parse_gt(kv.get("GT", ".")),
        gb_bp=gb_bp,
        q=_to_float(kv.get("Q")),
        dp=_to_int(kv.get("DP")),
        dstutter=_to_int(kv.get("DSTUTTER")),
        dflankindel=_to_int(kv.get("DFLANKINDEL")),
        mallreads=kv.get("MALLREADS", "."),
        allreads=kv.get("ALLREADS", "."),
        inframe_up=_to_float(info_d.get("INFRAME_UP")),
        inframe_down=_to_float(info_d.get("INFRAME_DOWN")),
        inframe_pgeom=_to_float(info_d.get("INFRAME_PGEOM")),
    )


def parse_hipstr_vcf(path: str | Path, sample: str) -> list[HipstrCall]:
    """Parse all records for one sample."""
    calls: list[HipstrCall] = []
    sample_cols: list[str] = []
    with _open_text(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                sample_cols = line.rstrip("\n").split("\t")[9:]
                continue
            if not line.strip():
                continue
            rec = parse_hipstr_record(line, sample_cols, sample)
            if rec is not None:
                calls.append(rec)
    return calls


def match_call_to_locus(calls: list[HipstrCall], locus: Locus, pad: int = 10) -> HipstrCall | None:
    """Pick the HipSTR record overlapping a locus (by chrom + position window)."""
    for c in calls:
        if c.chrom != locus.chrom:
            continue
        if (locus.ref_start - pad) <= c.pos <= (locus.ref_end + pad):
            return c
    return None
