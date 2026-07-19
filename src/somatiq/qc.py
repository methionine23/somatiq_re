"""Level-1 (per-sample, data-based) QC for the P0 score row."""
from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class QCConfig:
    min_informative_reads: int = 20
    min_q_warn: float = 0.9
    min_q_fail: float = 0.5


@dataclass(frozen=True)
class QCResult:
    flag: str  # PASS | WARN | FAIL
    reasons: tuple[str, ...] = ()
    detectable_f_floor: float | None = None


def level1_qc(
    *,
    informative_reads: int | None,
    stutter_learned: bool,
    germline_q: float | None,
    allele_exceeds_readlen: bool,
    interruption_flag: bool,
    config: QCConfig | None = None,
) -> QCResult:
    cfg = config or QCConfig()
    reasons: list[str] = []
    fail = False
    warn = False

    n = informative_reads or 0
    floor = (1.0 / n) if n > 0 else None
    if n < cfg.min_informative_reads:
        fail = True
        reasons.append(f"low_depth:{n}<{cfg.min_informative_reads}")

    if germline_q is not None:
        if germline_q < cfg.min_q_fail:
            fail = True
            reasons.append(f"germline_q_fail:{germline_q:.2f}")
        elif germline_q < cfg.min_q_warn:
            warn = True
            reasons.append(f"germline_q_low:{germline_q:.2f}")

    if not stutter_learned:
        warn = True
        reasons.append("stutter_default_model")

    if allele_exceeds_readlen:
        warn = True
        reasons.append("allele_exceeds_readlen")

    if interruption_flag:
        reasons.append("interruption_present")  # annotation, not a downgrade

    flag = "FAIL" if fail else ("WARN" if warn else "PASS")
    return QCResult(flag=flag, reasons=tuple(reasons), detectable_f_floor=floor)
