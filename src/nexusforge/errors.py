"""Structured error taxonomy.

Every failure in NexusForge is reported as a typed error carrying:

* a stable machine-readable :class:`ErrorCode`,
* the :class:`PipelineStage` where the failure happened,
* the data source involved (when applicable),
* optional structured details for debugging / reproduction,
* whether the failure is retryable.

This is what makes a failed run diagnosable: retrieval, sequence mapping,
structure parsing and evaluation failures are never collapsed into one opaque
"something went wrong" message.
"""

from __future__ import annotations

from enum import StrEnum
from typing import Any


class ErrorCode(StrEnum):
    """Stable error codes. Values are part of the public report contract."""

    # --- input ---
    INVALID_INPUT = "INVALID_INPUT"

    # --- retrieval ---
    SOURCE_NOT_FOUND = "SOURCE_NOT_FOUND"
    SOURCE_UNAVAILABLE = "SOURCE_UNAVAILABLE"
    NETWORK_ERROR = "NETWORK_ERROR"
    RATE_LIMITED = "RATE_LIMITED"
    INVALID_API_RESPONSE = "INVALID_API_RESPONSE"

    # --- scientific validation ---
    INVALID_SEQUENCE = "INVALID_SEQUENCE"
    SEQUENCE_MISMATCH = "SEQUENCE_MISMATCH"
    CHAIN_MAPPING_AMBIGUOUS = "CHAIN_MAPPING_AMBIGUOUS"
    NO_BASELINE = "NO_BASELINE"
    INVALID_COORDINATES = "INVALID_COORDINATES"
    CONFIDENCE_LENGTH_MISMATCH = "CONFIDENCE_LENGTH_MISMATCH"

    # --- reporting / runtime ---
    REPORT_WRITE_FAILED = "REPORT_WRITE_FAILED"
    INTERNAL_ERROR = "INTERNAL_ERROR"


class PipelineStage(StrEnum):
    """Where in the pipeline a failure occurred."""

    INPUT_VALIDATION = "input_validation"
    DATA_RETRIEVAL = "data_retrieval"
    SEQUENCE_VALIDATION = "sequence_validation"
    STRUCTURE_PARSING = "structure_parsing"
    STRUCTURE_EVALUATION = "structure_evaluation"
    CONFIDENCE_EVALUATION = "confidence_evaluation"
    REPORTING = "reporting"
    INTERNAL = "internal"


_CODE_STAGES: dict[ErrorCode, PipelineStage] = {
    ErrorCode.INVALID_INPUT: PipelineStage.INPUT_VALIDATION,
    ErrorCode.SOURCE_NOT_FOUND: PipelineStage.DATA_RETRIEVAL,
    ErrorCode.SOURCE_UNAVAILABLE: PipelineStage.DATA_RETRIEVAL,
    ErrorCode.NETWORK_ERROR: PipelineStage.DATA_RETRIEVAL,
    ErrorCode.RATE_LIMITED: PipelineStage.DATA_RETRIEVAL,
    ErrorCode.INVALID_API_RESPONSE: PipelineStage.DATA_RETRIEVAL,
    ErrorCode.INVALID_SEQUENCE: PipelineStage.SEQUENCE_VALIDATION,
    ErrorCode.SEQUENCE_MISMATCH: PipelineStage.SEQUENCE_VALIDATION,
    ErrorCode.CHAIN_MAPPING_AMBIGUOUS: PipelineStage.SEQUENCE_VALIDATION,
    ErrorCode.INVALID_COORDINATES: PipelineStage.STRUCTURE_PARSING,
    ErrorCode.NO_BASELINE: PipelineStage.STRUCTURE_PARSING,
    ErrorCode.CONFIDENCE_LENGTH_MISMATCH: PipelineStage.CONFIDENCE_EVALUATION,
    ErrorCode.REPORT_WRITE_FAILED: PipelineStage.REPORTING,
    ErrorCode.INTERNAL_ERROR: PipelineStage.INTERNAL,
}

_RETRYABLE_CODES = frozenset(
    {
        ErrorCode.NETWORK_ERROR,
        ErrorCode.RATE_LIMITED,
        ErrorCode.SOURCE_UNAVAILABLE,
    }
)

#: Process exit code per error family. Documented in the README so scripts can
#: branch on the failure kind instead of parsing log text.
_EXIT_CODES = {
    ErrorCode.INVALID_INPUT: 2,
    ErrorCode.SOURCE_NOT_FOUND: 3,
    ErrorCode.SOURCE_UNAVAILABLE: 4,
    ErrorCode.NETWORK_ERROR: 4,
    ErrorCode.RATE_LIMITED: 4,
    ErrorCode.INVALID_API_RESPONSE: 3,
    ErrorCode.INVALID_SEQUENCE: 3,
    ErrorCode.SEQUENCE_MISMATCH: 3,
    ErrorCode.CHAIN_MAPPING_AMBIGUOUS: 3,
    ErrorCode.NO_BASELINE: 3,
    ErrorCode.INVALID_COORDINATES: 3,
    ErrorCode.CONFIDENCE_LENGTH_MISMATCH: 3,
    ErrorCode.REPORT_WRITE_FAILED: 5,
    ErrorCode.INTERNAL_ERROR: 5,
}
EXIT_OK = 0
EXIT_UNEXPECTED = 5


class NexusForgeError(Exception):
    """Base class for every structured NexusForge failure."""

    default_stage: PipelineStage = PipelineStage.INTERNAL

    def __init__(
        self,
        code: ErrorCode,
        message: str,
        *,
        stage: PipelineStage | None = None,
        source: str | None = None,
        details: dict[str, Any] | None = None,
        http_status: int | None = None,
    ) -> None:
        self.code = code
        if stage is not None:
            self.stage = stage
        elif type(self).default_stage is not PipelineStage.INTERNAL:
            self.stage = type(self).default_stage
        else:
            self.stage = _CODE_STAGES.get(code, PipelineStage.INTERNAL)
        self.message = message
        self.source = source
        self.details: dict[str, Any] = details or {}
        self.http_status = http_status
        self.retryable = code in _RETRYABLE_CODES
        super().__init__(self.__str__())

    @property
    def exit_code(self) -> int:
        """Process exit code a CLI invocation should return for this error."""
        return _EXIT_CODES.get(self.code, EXIT_UNEXPECTED)

    def to_dict(self) -> dict[str, Any]:
        payload: dict[str, Any] = {
            "code": self.code.value,
            "stage": self.stage.value,
            "message": self.message,
            "retryable": self.retryable,
        }
        if self.source is not None:
            payload["source"] = self.source
        if self.http_status is not None:
            payload["http_status"] = self.http_status
        if self.details:
            payload["details"] = self.details
        return payload

    def __str__(self) -> str:
        where = f" [{self.stage.value}]"
        origin = f" source={self.source}" if self.source else ""
        return f"{self.code.value}{where}{origin}: {self.message}"

    def __repr__(self) -> str:
        return f"{type(self).__name__}(code={self.code.value!r}, stage={self.stage.value!r})"


class InputError(NexusForgeError):
    """A user-supplied identifier or option is not acceptable."""

    default_stage = PipelineStage.INPUT_VALIDATION


class SourceError(NexusForgeError):
    """A remote data source could not deliver the requested data."""

    default_stage = PipelineStage.DATA_RETRIEVAL


class ValidationError(NexusForgeError):
    """Data was retrieved but does not satisfy scientific preconditions."""

    default_stage = PipelineStage.SEQUENCE_VALIDATION


class StructureError(ValidationError):
    """Structure parsing or coordinate handling failed."""

    default_stage = PipelineStage.STRUCTURE_PARSING


class ConfidenceError(ValidationError):
    """Per-residue confidence data is missing or inconsistent."""

    default_stage = PipelineStage.CONFIDENCE_EVALUATION


class ReportingError(NexusForgeError):
    """Run artifacts could not be written."""

    default_stage = PipelineStage.REPORTING


__all__ = [
    "EXIT_OK",
    "EXIT_UNEXPECTED",
    "ConfidenceError",
    "ErrorCode",
    "InputError",
    "NexusForgeError",
    "PipelineStage",
    "ReportingError",
    "SourceError",
    "StructureError",
    "ValidationError",
]
