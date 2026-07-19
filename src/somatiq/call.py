"""Orchestration: HipSTR -> prancSTR -> scored rows for one sample.

The heavy calls (HipSTR, prancSTR, pysam) run where the data lives. The pure merge
logic lives in schema.assemble_score_row and is unit-tested; this module wires I/O.
"""
from __future__ import annotations

import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

from .catalog import Catalog, Locus
from .eh import parse_eh_vcf
from .hipstr import (
    HipstrCall,
    build_hipstr_command,
    match_call_to_locus,
    parse_hipstr_vcf,
    write_hipstr_regions,
)
from .prancstr import PrancResult, build_prancstr_command, parse_prancstr_tab
from .qc import QCConfig
from .readlen import detect_read_length
from .schema import ScoreRow, assemble_score_row


def match_pranc_to_locus(
    rows: list[PrancResult], locus: Locus, pad: int = 10
) -> PrancResult | None:
    """Match a prancSTR row to a locus by region name, else by chrom+position."""
    for r in rows:
        if r.locus and r.locus == locus.hipstr_region.name:
            return r
    for r in rows:
        if r.chrom == locus.chrom and (locus.ref_start - pad) <= r.pos <= (locus.ref_end + pad):
            return r
    return None


@dataclass
class RunConfig:
    data_mode: str = "wgs"
    def_stutter_model: bool = True
    eh_annotate: bool = False
    eh_vcf: str | None = None
    eh_thresholds: dict[str, float] | None = None
    qc: QCConfig | None = None
    keep_intermediate: bool = False


def run_sample(
    *,
    alignment: str,
    reference: str,
    catalog: Catalog,
    sample: str,
    sex: str | None,
    config: RunConfig | None = None,
    workdir: str | None = None,
    runner=subprocess.run,
) -> list[ScoreRow]:
    cfg = config or RunConfig()
    tmp = Path(workdir) if workdir else Path(tempfile.mkdtemp(prefix="somatiq_"))
    tmp.mkdir(parents=True, exist_ok=True)

    regions = write_hipstr_regions(catalog, tmp / "regions.bed")
    hipstr_vcf = str(tmp / f"{sample}.hipstr.vcf.gz")
    runner(
        build_hipstr_command(
            [alignment], reference, str(regions), hipstr_vcf,
            def_stutter_model=cfg.def_stutter_model,
        ),
        check=True,
    )

    pranc_prefix = str(tmp / f"{sample}.prancstr")
    runner(
        build_prancstr_command(hipstr_vcf, pranc_prefix, samples=[sample]),
        check=True,
    )

    calls = parse_hipstr_vcf(hipstr_vcf, sample)
    pranc_rows = parse_prancstr_tab(pranc_prefix + ".tab")

    # read length from the first locus region only (targeted; no full scan)
    first = catalog.loci[0]
    read_length = detect_read_length(
        alignment, reference, (first.chrom, first.ref_start, first.ref_end)
    )

    eh_ann = {}
    if cfg.eh_annotate and cfg.eh_vcf:
        eh_ann = parse_eh_vcf(cfg.eh_vcf, sample, cfg.eh_thresholds)

    rows = assemble_rows(
        sample=sample, sex=sex, catalog=catalog, calls=calls, pranc_rows=pranc_rows,
        read_length=read_length, data_mode=cfg.data_mode, eh_ann=eh_ann, qc=cfg.qc,
    )
    return rows


def assemble_rows(
    *,
    sample: str,
    sex: str | None,
    catalog: Catalog,
    calls: list[HipstrCall],
    pranc_rows: list[PrancResult],
    read_length: int | None,
    data_mode: str,
    eh_ann: dict | None = None,
    qc: QCConfig | None = None,
) -> list[ScoreRow]:
    """Pure assembly across loci (unit-tested)."""
    eh_ann = eh_ann or {}
    out: list[ScoreRow] = []
    for locus in catalog.loci:
        hip = match_call_to_locus(calls, locus)
        pr = match_pranc_to_locus(pranc_rows, locus)
        out.append(
            assemble_score_row(
                sample=sample,
                locus=locus,
                sex=sex,
                data_mode=data_mode,
                hipstr=hip,
                pranc=pr,
                read_length=read_length,
                eh_annotation=eh_ann.get(locus.eh_variant_id) or eh_ann.get(locus.locus_id),
                qc_config=qc,
            )
        )
    return out
