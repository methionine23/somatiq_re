"""Optional annotation from PRE-COMPUTED ExpansionHunter results.

We do NOT run ExpansionHunter — it is the owner's upstream sample selector. When
--eh-annotate is set we only parse a provided EH VCF to attach short/long allele
sizes (and, if a per-locus threshold is supplied, a carrier flag). Default: off.
"""
from __future__ import annotations

import gzip
from pathlib import Path


def _open_text(path: str | Path):
    path = str(path)
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def _info_dict(info: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for f in info.split(";"):
        if "=" in f:
            k, _, v = f.partition("=")
            out[k] = v
    return out


def parse_eh_vcf(
    path: str | Path,
    sample: str,
    thresholds: dict[str, float] | None = None,
) -> dict[str, dict]:
    """Return {locus_id: {gt_short, gt_long, carrier}} from an EH VCF.

    EH encodes the locus in INFO REPID/VARID and repeat copy numbers in the FORMAT
    REPCN field ('n1/n2'). carrier is set only if a threshold for that locus is given.
    """
    thresholds = thresholds or {}
    out: dict[str, dict] = {}
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
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            info = _info_dict(f[7])
            locus_id = info.get("REPID") or info.get("VARID")
            if not locus_id:
                continue
            try:
                col = 9 + sample_cols.index(sample)
            except ValueError:
                continue
            kv = dict(zip(f[8].split(":"), f[col].split(":")))
            repcn = kv.get("REPCN", ".")
            if repcn in (".", ""):
                continue
            copies = []
            for a in repcn.replace("|", "/").split("/"):
                try:
                    copies.append(float(a))
                except ValueError:
                    pass
            if not copies:
                continue
            short, long = min(copies), max(copies)
            carrier = None
            if locus_id in thresholds:
                carrier = long >= thresholds[locus_id]
            out[locus_id] = {"gt_short": short, "gt_long": long, "carrier": carrier}
    return out
