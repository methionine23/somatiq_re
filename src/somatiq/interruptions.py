"""Detect repeat interruptions in a HipSTR allele sequence (P0: annotation only).

An interruption is any period-sized unit inside the tract that differs from the
pure (reference-strand) motif. We report whether any were found, the count, and the
non-pure units observed (rotations of the pure motif are treated as pure, since the
tract phase is arbitrary).
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class InterruptionCall:
    flag: bool
    n_impure_units: int
    impure_units: tuple[str, ...]  # distinct non-pure units, in order of first sight
    pure_fraction: float | None  # pure units / total units


def _rotations(motif: str) -> set[str]:
    return {motif[i:] + motif[:i] for i in range(len(motif))}


def detect_interruptions(
    allele_seq: str,
    pure_motif: str,
    known_motifs: tuple[str, ...] = (),
) -> InterruptionCall:
    seq = (allele_seq or "").upper()
    period = len(pure_motif)
    if period == 0 or len(seq) < period:
        return InterruptionCall(False, 0, (), None)
    pure_set = _rotations(pure_motif.upper())
    n_units = len(seq) // period
    impure: list[str] = []
    n_impure = 0
    seen: list[str] = []
    for i in range(n_units):
        unit = seq[i * period : (i + 1) * period]
        if unit in pure_set:
            continue
        n_impure += 1
        if unit not in seen:
            seen.append(unit)
            impure.append(unit)
    pure_fraction = (n_units - n_impure) / n_units if n_units else None
    return InterruptionCall(
        flag=n_impure > 0,
        n_impure_units=n_impure,
        impure_units=tuple(impure),
        pure_fraction=pure_fraction,
    )


def summarize_alleles(
    allele_seqs: list[str], pure_motif: str, known_motifs: tuple[str, ...] = ()
) -> InterruptionCall:
    """Aggregate across a sample's called alleles: flag set if any allele impure."""
    calls = [detect_interruptions(s, pure_motif, known_motifs) for s in allele_seqs]
    if not calls:
        return InterruptionCall(False, 0, (), None)
    flag = any(c.flag for c in calls)
    n = sum(c.n_impure_units for c in calls)
    units: list[str] = []
    for c in calls:
        for u in c.impure_units:
            if u not in units:
                units.append(u)
    fracs = [c.pure_fraction for c in calls if c.pure_fraction is not None]
    pure_fraction = min(fracs) if fracs else None
    return InterruptionCall(flag, n, tuple(units), pure_fraction)
