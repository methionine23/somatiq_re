"""Per-sample read length from the TARGET-LOCUS reads only (no full-CRAM scan).

We reuse the same targeted region access the pipeline already needs, so this adds
no extra egress. pysam is imported lazily so the rest of the package works without
the CRAM stack installed.
"""
from __future__ import annotations

from collections import Counter
from typing import Iterable


def mode_read_length(lengths: Iterable[int]) -> int | None:
    """Most common read length (ties -> the larger, to be conservative on ceiling)."""
    counts = Counter(int(x) for x in lengths if x and x > 0)
    if not counts:
        return None
    best = max(counts.items(), key=lambda kv: (kv[1], kv[0]))
    return best[0]


def detect_read_length(
    alignment_path: str,
    reference: str | None,
    region: tuple[str, int, int],
    max_reads: int = 200,
) -> int | None:
    """Fetch up to max_reads at one small region and return the modal read length.

    region is (chrom, start_1based, end_1based). Requires a co-located index (.crai/
    .bai) and, for CRAM, the reference for decode.
    """
    try:
        import pysam  # noqa: PLC0415 (optional dep)
    except ImportError as e:  # pragma: no cover
        raise RuntimeError("pysam is required for read-length detection; pip install '.[cram]'") from e

    chrom, start, end = region
    lengths: list[int] = []
    kwargs = {"reference_filename": reference} if reference else {}
    with pysam.AlignmentFile(alignment_path, **kwargs) as af:
        for i, read in enumerate(af.fetch(chrom, max(0, start - 1), end)):
            if i >= max_reads:
                break
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            ln = read.infer_read_length() or read.query_length
            if ln:
                lengths.append(ln)
    return mode_read_length(lengths)
