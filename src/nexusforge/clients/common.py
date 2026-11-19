"""Helpers shared by the data-source clients."""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from datetime import UTC, datetime
from typing import Any

from nexusforge.errors import ErrorCode, SourceError


def parse_json_object(payload: bytes, *, source: str, url: str) -> Any:
    """Decode a JSON payload, turning decode errors into a structured failure."""
    try:
        return json.loads(payload.decode("utf-8"))
    except UnicodeDecodeError as exc:
        raise SourceError(
            ErrorCode.INVALID_API_RESPONSE,
            f"{source} returned non-UTF-8 content",
            source=source,
            details={"url": url, "cause": str(exc)},
        ) from exc
    except json.JSONDecodeError as exc:
        preview = payload[:200].decode("utf-8", errors="replace")
        raise SourceError(
            ErrorCode.INVALID_API_RESPONSE,
            f"{source} returned invalid JSON ({exc.msg} at char {exc.pos})",
            source=source,
            details={"url": url, "body_preview": preview},
        ) from exc


def pick(
    entry: Mapping[str, Any],
    names: Sequence[str],
    *,
    source: str,
    url: str,
    required: bool = True,
    notes: list[str] | None = None,
) -> Any | None:
    """Read the first present field among ``names``.

    ``names[0]`` is the field documented by the current API schema; the rest are
    tolerated fallbacks. When a fallback is used the fact is recorded in
    ``notes`` so the run report shows that the source schema differs from the
    documented one — instead of silently guessing.
    """
    for index, name in enumerate(names):
        if name in entry and entry[name] is not None:
            if index > 0 and notes is not None:
                notes.append(
                    f"field '{names[0]}' missing from {source}; used fallback field '{name}'"
                )
            return entry[name]
    if required:
        raise SourceError(
            ErrorCode.INVALID_API_RESPONSE,
            f"{source} response is missing the required field '{names[0]}'",
            source=source,
            details={"url": url, "available_fields": sorted(entry.keys())[:40]},
        )
    return None


def require_mapping(value: Any, *, field: str, source: str, url: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise SourceError(
            ErrorCode.INVALID_API_RESPONSE,
            f"{source} field '{field}' is not an object",
            source=source,
            details={"url": url, "actual_type": type(value).__name__},
        )
    return value


def parse_datetime(value: Any) -> datetime | None:
    """Parse ISO-8601 timestamps seen in UniProt/AlphaFold/RCSB payloads."""
    if value is None:
        return None
    if isinstance(value, datetime):
        return value if value.tzinfo else value.replace(tzinfo=UTC)
    if not isinstance(value, str) or not value.strip():
        return None
    raw = value.strip()
    if raw.endswith("Z"):
        raw = raw[:-1] + "+00:00"
    try:
        parsed = datetime.fromisoformat(raw)
    except ValueError:
        return None
    return parsed if parsed.tzinfo else parsed.replace(tzinfo=UTC)


def require_float(value: Any, *, field: str, source: str) -> float | None:
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise SourceError(
            ErrorCode.INVALID_API_RESPONSE,
            f"{source} field '{field}' is not numeric",
            source=source,
            details={"actual_type": type(value).__name__},
        )
    return float(value)


__all__ = [
    "parse_datetime",
    "parse_json_object",
    "pick",
    "require_float",
    "require_mapping",
]
