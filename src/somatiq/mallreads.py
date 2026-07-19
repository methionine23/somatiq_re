"""Parse HipSTR MALLREADS / ALLREADS per-read length evidence.

HipSTR encodes per-read support as a ';'-separated list of 'bpdiff|count', where
bpdiff is the base-pair difference of the read's (maximum-likelihood) alignment from
the reference allele. Repeat copy number = ref_copies + bpdiff/period.
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ReadHist:
    """Histogram of per-read length evidence at one locus for one sample."""

    # bp difference from reference -> read count
    by_bpdiff: dict[int, int]
    period: int
    ref_copies: int

    @property
    def informative_reads(self) -> int:
        return sum(self.by_bpdiff.values())

    def by_copies(self) -> dict[float, int]:
        """Map repeat-copy-number -> read count (copies may be fractional if a
        bpdiff is not a clean multiple of the period, e.g. an indel/interruption)."""
        out: dict[float, int] = {}
        for bpdiff, count in self.by_bpdiff.items():
            copies = self.ref_copies + bpdiff / self.period
            out[copies] = out.get(copies, 0) + count
        return out

    def modal_bpdiff(self) -> int | None:
        if not self.by_bpdiff:
            return None
        return max(self.by_bpdiff.items(), key=lambda kv: (kv[1], -abs(kv[0])))[0]

    def max_copies(self) -> float | None:
        cps = self.by_copies()
        return max(cps) if cps else None

    def detectable_f_floor(self) -> float | None:
        n = self.informative_reads
        return (1.0 / n) if n > 0 else None


def parse_read_field(value: str, period: int, ref_copies: int) -> ReadHist:
    """Parse a MALLREADS/ALLREADS FORMAT value into a ReadHist.

    Accepts '.', empty, or 'bpdiff|count;bpdiff|count'. Malformed entries are
    skipped rather than raising, so one bad token does not drop a whole sample.
    """
    by: dict[int, int] = {}
    if value and value not in (".", ""):
        for token in value.split(";"):
            token = token.strip()
            if not token or "|" not in token:
                continue
            bp_s, _, ct_s = token.partition("|")
            try:
                bpdiff = int(bp_s)
                count = int(ct_s)
            except ValueError:
                continue
            if count <= 0:
                continue
            by[bpdiff] = by.get(bpdiff, 0) + count
    return ReadHist(by_bpdiff=by, period=period, ref_copies=ref_copies)
