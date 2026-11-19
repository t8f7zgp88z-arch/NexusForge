"""Sequence-level validation.

This module answers one question before any structural number is computed:
*does this chain actually correspond to this UniProt sequence, and over which
region?* When the answer is ambiguous the pipeline stops instead of producing a
plausible-looking but meaningless metric.
"""

from __future__ import annotations

import warnings
from typing import Any

from Bio.Align import PairwiseAligner

from nexusforge.errors import ErrorCode, ValidationError
from nexusforge.models import AlignmentResult, SequenceStats

STANDARD_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")
AMBIGUOUS_AMINO_ACIDS = frozenset("BJOUXZ")
#: Anything else (including gap characters, '*' and whitespace) is invalid in a
#: raw sequence record.
_ALLOWED = STANDARD_AMINO_ACIDS | AMBIGUOUS_AMINO_ACIDS

#: Alignment scoring tuned for near-identical protein chains: a single-residue
#: terminal truncation must remain one gap, not be spread over several columns.
MATCH_SCORE = 2.0
MISMATCH_SCORE = -3.0
GAP_OPEN_SCORE = -8.0
GAP_EXTEND_SCORE = -0.5

#: Counting optimal alignments is exponential-ish on long, divergent pairs; past
#: this many cells the count is skipped and reported as 1.
_ALIGNMENT_COUNT_LIMIT = 2_000_000


def validate_sequence(sequence: str, *, label: str = "sequence") -> SequenceStats:
    """Validate a raw amino-acid sequence, raising ``INVALID_SEQUENCE`` if unusable."""
    normalized = (sequence or "").strip().upper()
    if not normalized:
        raise ValidationError(
            ErrorCode.INVALID_SEQUENCE,
            f"{label} is empty",
            details={"label": label},
        )
    invalid = sorted({char for char in normalized if char not in _ALLOWED})
    if invalid:
        raise ValidationError(
            ErrorCode.INVALID_SEQUENCE,
            f"{label} contains {len(invalid)} non-amino-acid character(s): {''.join(invalid)}",
            details={"label": label, "invalid_characters": invalid},
        )
    ambiguous = sorted({char for char in normalized if char in AMBIGUOUS_AMINO_ACIDS})
    return SequenceStats(
        length=len(normalized), invalid_characters=[], ambiguous_characters=ambiguous
    )


def build_aligner() -> PairwiseAligner:
    """Global protein aligner with affine gaps (deterministic scoring)."""
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = MATCH_SCORE
    aligner.mismatch_score = MISMATCH_SCORE
    aligner.open_gap_score = GAP_OPEN_SCORE
    aligner.extend_gap_score = GAP_EXTEND_SCORE
    return aligner


def _count_optimal(alignment: Any, query_length: int, target_length: int) -> int:
    if query_length * target_length > _ALIGNMENT_COUNT_LIMIT:
        return 1
    try:
        return max(1, len(alignment))
    except (MemoryError, OverflowError):  # pragma: no cover - defensive
        return 1


def align_sequences(
    query: str,
    target: str,
    *,
    query_id: str,
    target_id: str,
    min_identity: float | None = None,
    min_coverage: float | None = None,
    preview_length: int = 80,
) -> AlignmentResult:
    """Globally align two protein sequences and summarise the relationship.

    ``mapping`` maps 0-based query indices to 0-based target indices for aligned
    (non-gap) residue pairs, which is what downstream residue matching uses.
    """
    query = query.strip().upper()
    target = target.strip().upper()
    if not query or not target:
        raise ValidationError(
            ErrorCode.INVALID_SEQUENCE,
            "cannot align an empty sequence",
            details={"query_id": query_id, "target_id": target_id},
        )

    aligner = build_aligner()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        alignments = aligner.align(query, target)
    n_optimal = _count_optimal(alignments, len(query), len(target))
    alignment = alignments[0]
    gapped_query = str(alignment[0])
    gapped_target = str(alignment[1])

    matches = mismatches = gap_positions = gap_openings = 0
    mapping: dict[int, int] = {}
    query_index = target_index = 0
    in_gap = False
    for q_char, t_char in zip(gapped_query, gapped_target, strict=True):
        if q_char == "-":
            gap_positions += 1
            target_index += 1
            in_gap = True
            continue
        if t_char == "-":
            gap_positions += 1
            query_index += 1
            in_gap = True
            continue
        if in_gap:
            gap_openings += 1
            in_gap = False
        mapping[query_index] = target_index
        if q_char == t_char:
            matches += 1
        else:
            mismatches += 1
        query_index += 1
        target_index += 1
    if in_gap:
        gap_openings += 1

    aligned_columns = matches + mismatches
    identity = matches / aligned_columns if aligned_columns else 0.0
    coverage_query = len(mapping) / len(query)
    coverage_target = len(mapping) / len(target)

    result = AlignmentResult(
        query_id=query_id,
        target_id=target_id,
        query_length=len(query),
        target_length=len(target),
        aligned_length=len(gapped_query),
        matches=matches,
        mismatches=mismatches,
        gap_openings=gap_openings,
        gap_positions=gap_positions,
        identity=round(identity, 6),
        coverage_query=round(coverage_query, 6),
        coverage_target=round(coverage_target, 6),
        n_optimal_alignments=n_optimal,
        mapping=mapping,
        aligned_query_preview=gapped_query[:preview_length],
        aligned_target_preview=gapped_target[:preview_length],
    )
    if min_identity is not None or min_coverage is not None:
        require_compatible(
            result,
            min_identity=min_identity if min_identity is not None else 0.0,
            min_coverage=min_coverage if min_coverage is not None else 0.0,
        )
    return result


def require_compatible(
    alignment: AlignmentResult,
    *,
    min_identity: float,
    min_coverage: float,
) -> AlignmentResult:
    """Raise ``SEQUENCE_MISMATCH`` unless identity and coverage clear the floors."""
    if alignment.identity < min_identity:
        raise ValidationError(
            ErrorCode.SEQUENCE_MISMATCH,
            f"{alignment.query_id} and {alignment.target_id} share only "
            f"{alignment.identity:.1%} identity over aligned residues "
            f"(minimum {min_identity:.0%}); refusing to treat them as the same protein",
            details={
                "identity": alignment.identity,
                "min_identity": min_identity,
                "query_id": alignment.query_id,
                "target_id": alignment.target_id,
            },
        )
    if alignment.coverage_query < min_coverage:
        raise ValidationError(
            ErrorCode.SEQUENCE_MISMATCH,
            f"only {alignment.coverage_query:.1%} of {alignment.query_id} is covered by "
            f"aligned residues (minimum {min_coverage:.0%})",
            details={
                "coverage_query": alignment.coverage_query,
                "min_coverage": min_coverage,
                "query_id": alignment.query_id,
                "target_id": alignment.target_id,
            },
        )
    return alignment


def describe_gaps(alignment: AlignmentResult) -> str:
    """Human-readable gap summary used in reports."""
    if alignment.gap_openings == 0:
        return "no gaps"
    return (
        f"{alignment.gap_openings} gap(s), {alignment.gap_positions} residue position(s) "
        f"absent from one side"
    )


__all__ = [
    "AMBIGUOUS_AMINO_ACIDS",
    "STANDARD_AMINO_ACIDS",
    "align_sequences",
    "build_aligner",
    "describe_gaps",
    "require_compatible",
    "validate_sequence",
]
