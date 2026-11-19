"""Confidence evaluation.

pLDDT is AlphaFold's per-residue *prediction confidence*. It is not a measured
accuracy and it is not ground truth. This module

1. validates pLDDT values (count, range, residue numbering),
2. cross-checks the confidence JSON against the B-factor column of the model
   file (AlphaFold stores pLDDT there),
3. relates pLDDT bins to the experimentally referenced ``lDDT-Cα`` so the report
   can show where confidence did and did not correspond to measured agreement.
"""

from __future__ import annotations

import numpy as np

from nexusforge.errors import ConfidenceError, ErrorCode
from nexusforge.models import ConfidenceEvaluation, ConfidenceScores, PlddtBin
from nexusforge.validation.structure import CaTrace

PLDDT_BINS: tuple[tuple[str, float, float], ...] = (
    ("Very low", 0.0, 50.0),
    ("Low", 50.0, 70.0),
    ("Confident", 70.0, 90.0),
    ("Very high", 90.0, 100.0000001),
)


def validate_confidence(
    scores: ConfidenceScores,
    *,
    sequence_length: int,
    sequence_start: int | None = None,
    trace_label: str = "prediction",
) -> dict[int, float]:
    """Validate pLDDT coverage against the modelled sequence.

    Returns a mapping ``residue_number -> pLDDT``. A mismatch in count, or a
    numbering scheme that does not cover the modelled window exactly, raises
    ``CONFIDENCE_LENGTH_MISMATCH`` instead of being silently truncated/padded.
    """
    start = sequence_start or 1
    expected_numbers = list(range(start, start + sequence_length))
    if len(scores.scores) != sequence_length:
        raise ConfidenceError(
            ErrorCode.CONFIDENCE_LENGTH_MISMATCH,
            f"{trace_label} has {sequence_length} modelled residue(s) but "
            f"{len(scores.scores)} pLDDT value(s)",
            details={
                "sequence_length": sequence_length,
                "confidence_count": len(scores.scores),
                "sequence_start": start,
            },
        )
    actual = scores.residue_numbers
    if actual and actual != expected_numbers:
        if len(actual) == len(expected_numbers) and actual[0] != expected_numbers[0]:
            raise ConfidenceError(
                ErrorCode.CONFIDENCE_LENGTH_MISMATCH,
                f"{trace_label} pLDDT numbering starts at {actual[0]} but the modelled "
                f"sequence starts at {expected_numbers[0]}",
                details={
                    "confidence_first": actual[0],
                    "expected_first": expected_numbers[0],
                    "count": len(actual),
                },
            )
        raise ConfidenceError(
            ErrorCode.CONFIDENCE_LENGTH_MISMATCH,
            f"{trace_label} pLDDT residue numbering is not the expected contiguous "
            f"range {expected_numbers[0]}..{expected_numbers[-1]}",
            details={"confidence_first": actual[0], "confidence_last": actual[-1]},
        )
    values = np.asarray(scores.scores, dtype=float)
    if not np.isfinite(values).all():
        raise ConfidenceError(
            ErrorCode.INVALID_API_RESPONSE, "pLDDT array contains non-finite values"
        )
    if values.min() < 0.0 or values.max() > 100.0:
        raise ConfidenceError(
            ErrorCode.INVALID_API_RESPONSE,
            f"pLDDT values outside [0, 100] (min {values.min():.3f}, max {values.max():.3f})",
        )
    return dict(zip(actual, scores.scores, strict=True))


def plddt_from_trace(trace: CaTrace, *, source: str = "alphafold_cif_bfactor") -> ConfidenceScores:
    """Read per-residue pLDDT from a model file's B-factor column.

    Only valid for AlphaFold models, where the B-factor column stores pLDDT.
    """
    values = np.asarray(trace.bfactors, dtype=float)
    if values.size == 0:
        raise ConfidenceError(
            ErrorCode.CONFIDENCE_LENGTH_MISMATCH, f"{trace.label} trace carries no B-factors"
        )
    if values.min() < 0.0 or values.max() > 100.0:
        raise ConfidenceError(
            ErrorCode.INVALID_API_RESPONSE,
            f"{trace.label} B-factors fall outside the pLDDT range [0, 100] "
            f"(min {values.min():.2f}, max {values.max():.2f}); this does not look like an "
            "AlphaFold model file",
        )
    return ConfidenceScores(
        source=source,  # type: ignore[arg-type]
        residue_numbers=list(trace.residue_numbers),
        scores=[float(v) for v in values],
    )


def cross_check_confidence(primary: dict[int, float], secondary: dict[int, float]) -> float | None:
    """Largest absolute difference on shared residue numbers, or ``None``."""
    shared = sorted(set(primary) & set(secondary))
    if not shared:
        return None
    return max(abs(primary[number] - secondary[number]) for number in shared)


def _average_ranks(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=float)
    sorted_values = values[order]
    index = 0
    while index < len(values):
        end = index
        while end + 1 < len(values) and sorted_values[end + 1] == sorted_values[index]:
            end += 1
        ranks[order[index : end + 1]] = (index + end) / 2.0 + 1.0
        index = end + 1
    return ranks


def _pearson(x: np.ndarray, y: np.ndarray) -> float | None:
    if len(x) < 3:
        return None
    x_centered = x - x.mean()
    y_centered = y - y.mean()
    denominator = float(np.sqrt(np.sum(x_centered**2) * np.sum(y_centered**2)))
    if denominator == 0.0:
        return None
    return float(np.sum(x_centered * y_centered) / denominator)


def summarise_confidence(
    plddt: dict[int, float],
    lddt_ca: dict[int, float],
    *,
    source: str,
    cross_check_max_abs_delta: float | None = None,
) -> ConfidenceEvaluation:
    """Bin pLDDT and relate each bin to measured per-residue lDDT-Cα."""
    if not plddt:
        raise ConfidenceError(
            ErrorCode.CONFIDENCE_LENGTH_MISMATCH, "no pLDDT values available to summarise"
        )
    positions = sorted(plddt)
    values = np.asarray([plddt[position] for position in positions], dtype=float)
    total = len(values)

    bins: list[PlddtBin] = []
    for label, lower, upper in PLDDT_BINS:
        members = [p for p in positions if lower <= plddt[p] < upper]
        correlated = [p for p in members if p in lddt_ca and np.isfinite(lddt_ca[p])]
        bins.append(
            PlddtBin(
                label=label,
                lower=lower,
                upper=min(upper, 100.0),
                residue_count=len(members),
                fraction=len(members) / total,
                mean_lddt_ca=(
                    float(np.mean([lddt_ca[p] for p in correlated])) if correlated else None
                ),
            )
        )

    shared = [p for p in positions if p in lddt_ca and np.isfinite(lddt_ca[p])]
    pearson = spearman = None
    if len(shared) >= 3:
        x = np.asarray([plddt[p] for p in shared], dtype=float)
        y = np.asarray([lddt_ca[p] for p in shared], dtype=float)
        pearson = _pearson(x, y)
        spearman = _pearson(_average_ranks(x), _average_ranks(y))

    return ConfidenceEvaluation(
        source=source,
        residue_count=total,
        min_plddt=float(values.min()),
        max_plddt=float(values.max()),
        mean_plddt=float(values.mean()),
        median_plddt=float(np.median(values)),
        bins=bins,
        cross_check_max_abs_delta=cross_check_max_abs_delta,
        pearson_plddt_lddt=pearson,
        spearman_plddt_lddt=spearman,
        correlated_residue_count=len(shared),
    )


__all__ = [
    "PLDDT_BINS",
    "cross_check_confidence",
    "plddt_from_trace",
    "summarise_confidence",
    "validate_confidence",
]
