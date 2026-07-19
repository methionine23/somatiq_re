"""Per-sample x locus score row: assembly and output."""
from __future__ import annotations

import csv
from dataclasses import asdict, dataclass, fields
from pathlib import Path

from .catalog import Locus
from .hipstr import HipstrCall
from .interruptions import summarize_alleles
from .mallreads import parse_read_field
from .prancstr import PrancResult
from .qc import QCConfig, QCResult, level1_qc


@dataclass
class ScoreRow:
    sample_id: str
    locus_id: str
    disease: str
    sex: str | None
    ploidy_used: int | None
    data_mode: str
    # germline (HipSTR)
    germline_gt_copies: str | None  # e.g. "20|20"
    germline_gt_bp: str | None
    germline_q: float | None
    dp: int | None
    dstutter: int | None
    dflankindel: int | None
    stutter_up: float | None
    stutter_down: float | None
    stutter_pgeom: float | None
    stutter_model_source: str  # learned | default
    # somatic (prancSTR) -- the P0 score
    prancstr_f: float | None
    prancstr_mosaic_allele: float | None
    prancstr_pval: float | None
    mosaic_support: int | None
    # power / evidence
    informative_reads: int | None
    detectable_f_floor: float | None
    read_length: int | None
    allele_exceeds_readlen: bool
    # interruptions (annotation)
    interruption_flag: bool
    interruption_seq: str | None
    # optional EH annotation
    eh_carrier: bool | None
    eh_gt_short: float | None
    eh_gt_long: float | None
    # QC
    qc_level1: str
    qc_reasons: str

    @classmethod
    def columns(cls) -> list[str]:
        return [f.name for f in fields(cls)]


def _copies_str(vals) -> str | None:
    if not vals:
        return None
    return "|".join(f"{v:g}" for v in vals)


def assemble_score_row(
    *,
    sample: str,
    locus: Locus,
    sex: str | None,
    data_mode: str,
    hipstr: HipstrCall | None,
    pranc: PrancResult | None,
    read_length: int | None,
    eh_annotation: dict | None = None,
    qc_config: QCConfig | None = None,
) -> ScoreRow:
    """Merge parsed HipSTR + prancSTR (+ optional EH) into one scored row.

    Pure function: all I/O already done by the caller, so this is unit-tested.
    """
    ploidy = None
    try:
        ploidy = locus.resolve_ploidy(sex)
    except ValueError:
        ploidy = None

    # read-length evidence from MALLREADS
    hist = None
    informative = None
    max_copies = None
    if hipstr is not None:
        hist = parse_read_field(hipstr.mallreads, locus.period, locus.hipstr_region.ref_copies)
        informative = hist.informative_reads
        max_copies = hist.max_copies()

    # germline
    gt_copies = gt_bp = None
    germ_q = dp = dstutter = dflankindel = None
    up = down = pgeom = None
    stutter_source = "default"
    interruption_flag = False
    interruption_seq = None
    if hipstr is not None:
        gt_copies = _copies_str(hipstr.allele_copies(locus.period))
        gt_bp = _copies_str(hipstr.gb_bp) if hipstr.gb_bp else None
        germ_q, dp = hipstr.q, hipstr.dp
        dstutter, dflankindel = hipstr.dstutter, hipstr.dflankindel
        up, down, pgeom = hipstr.inframe_up, hipstr.inframe_down, hipstr.inframe_pgeom
        stutter_source = "learned" if hipstr.stutter_learned else "default"
        itr = summarize_alleles(
            hipstr.called_allele_seqs(), locus.interruption.pure_tract or locus.ref_motif,
            locus.interruption.motifs,
        )
        interruption_flag = itr.flag
        interruption_seq = ",".join(itr.impure_units) or None

    # ceiling: does any observed allele (germline or mosaic) exceed the read length?
    allele_exceeds = False
    if read_length is not None:
        candidate_bp = []
        if max_copies is not None:
            candidate_bp.append(max_copies * locus.period)
        if pranc is not None and pranc.mosaic_c is not None:
            candidate_bp.append(pranc.mosaic_c * locus.period)
        if hipstr is not None:
            candidate_bp.extend(len(s) for s in hipstr.called_allele_seqs())
        if candidate_bp and max(candidate_bp) > read_length:
            allele_exceeds = True

    # somatic score
    pf = pmos = ppval = msup = None
    if pranc is not None:
        pf, pmos, ppval, msup = pranc.f, pranc.mosaic_c, pranc.pval, pranc.mosaic_support

    qc = level1_qc(
        informative_reads=informative,
        stutter_learned=(stutter_source == "learned"),
        germline_q=germ_q,
        allele_exceeds_readlen=allele_exceeds,
        interruption_flag=interruption_flag,
        config=qc_config,
    )

    eh = eh_annotation or {}
    return ScoreRow(
        sample_id=sample,
        locus_id=locus.locus_id,
        disease=locus.disease,
        sex=sex,
        ploidy_used=ploidy,
        data_mode=data_mode,
        germline_gt_copies=gt_copies,
        germline_gt_bp=gt_bp,
        germline_q=germ_q,
        dp=dp,
        dstutter=dstutter,
        dflankindel=dflankindel,
        stutter_up=up,
        stutter_down=down,
        stutter_pgeom=pgeom,
        stutter_model_source=stutter_source,
        prancstr_f=pf,
        prancstr_mosaic_allele=pmos,
        prancstr_pval=ppval,
        mosaic_support=msup,
        informative_reads=informative,
        detectable_f_floor=qc.detectable_f_floor,
        read_length=read_length,
        allele_exceeds_readlen=allele_exceeds,
        interruption_flag=interruption_flag,
        interruption_seq=interruption_seq,
        eh_carrier=eh.get("carrier"),
        eh_gt_short=eh.get("gt_short"),
        eh_gt_long=eh.get("gt_long"),
        qc_level1=qc.flag,
        qc_reasons=";".join(qc.reasons),
    )


def _fmt(v) -> str:
    if v is None:
        return ""
    if isinstance(v, bool):
        return "1" if v else "0"
    if isinstance(v, float):
        return f"{v:g}"
    return str(v)


def write_tsv(rows: list[ScoreRow], path: str | Path) -> Path:
    path = Path(path)
    cols = ScoreRow.columns()
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(cols)
        for r in rows:
            d = asdict(r)
            w.writerow([_fmt(d[c]) for c in cols])
    return path


def write_parquet(rows: list[ScoreRow], path: str | Path) -> Path:
    """Optional Parquet output (requires the 'parquet' extra)."""
    import pandas as pd  # noqa: PLC0415 (optional dep)

    df = pd.DataFrame([asdict(r) for r in rows], columns=ScoreRow.columns())
    df.to_parquet(path, index=False)
    return Path(path)
