"""HTTP retrieval with bounded retries, rate-limit handling and provenance.

Design rules:

* **Bounded** exponential backoff — never retry forever.
* **Retry-After is honoured** when a source sends it, but capped so a run
  cannot hang for minutes because a server asked for a long wait.
* **4xx fails fast** (except 408/429) — a wrong accession must not be retried.
* **Every byte is hashed** with SHA-256 and timestamped, so reports can be
  reproduced and verified later.
* No jitter: retries are deterministic, which keeps the failure-injection tests
  exact and run logs comparable between runs.
"""

from __future__ import annotations

import hashlib
import logging
import time
from collections.abc import Callable
from dataclasses import dataclass, field
from datetime import UTC, datetime
from email.utils import parsedate_to_datetime
from pathlib import Path

import httpx

from nexusforge.config import Settings
from nexusforge.errors import ErrorCode, SourceError
from nexusforge.models import ArtifactRecord, Provenance

logger = logging.getLogger("nexusforge.http")

RETRYABLE_STATUS = frozenset({408, 425, 429, 500, 502, 503, 504})


def utcnow() -> datetime:
    return datetime.now(UTC)


def sha256_hex(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


@dataclass(frozen=True)
class AttemptRecord:
    """One HTTP attempt, kept for the run log and failure diagnosis."""

    attempt: int
    url: str
    status: int | None
    outcome: str
    delay_before_seconds: float = 0.0
    retry_after_seconds: float | None = None
    detail: str | None = None

    def as_dict(self) -> dict[str, object]:
        payload: dict[str, object] = {
            "attempt": self.attempt,
            "url": self.url,
            "status": self.status,
            "outcome": self.outcome,
            "delay_before_seconds": round(self.delay_before_seconds, 3),
        }
        if self.retry_after_seconds is not None:
            payload["retry_after_seconds"] = round(self.retry_after_seconds, 3)
        if self.detail:
            payload["detail"] = self.detail
        return payload


@dataclass
class FetchResult:
    """Raw bytes plus everything needed to cite and verify them."""

    content: bytes
    provenance: Provenance
    attempts: list[AttemptRecord] = field(default_factory=list)


def parse_retry_after(value: str | None, *, now: datetime | None = None) -> float | None:
    """Parse an HTTP ``Retry-After`` header (delta-seconds or HTTP-date)."""
    if not value:
        return None
    raw = value.strip()
    if not raw:
        return None
    try:
        seconds = float(raw)
    except ValueError:
        try:
            when = parsedate_to_datetime(raw)
        except (TypeError, ValueError):
            return None
        if when.tzinfo is None:
            when = when.replace(tzinfo=UTC)
        reference = now or utcnow()
        return max(0.0, (when - reference).total_seconds())
    return max(0.0, seconds)


class HttpClient:
    """Thin, testable HTTPX wrapper with a documented retry policy."""

    def __init__(
        self,
        settings: Settings,
        *,
        client: httpx.Client | None = None,
        sleep: Callable[[float], None] = time.sleep,
        clock: Callable[[], datetime] = utcnow,
    ) -> None:
        self._settings = settings
        self._sleep = sleep
        self._clock = clock
        self._owns_client = client is None
        self._client = client or httpx.Client(
            timeout=httpx.Timeout(settings.timeout_seconds),
            follow_redirects=True,
            headers={"User-Agent": settings.effective_user_agent},
        )
        self.attempts: list[AttemptRecord] = []
        #: Every successful response, in fetch order. The pipeline archives these
        #: bytes under ``raw/`` so a run keeps its exact inputs.
        self.responses: list[FetchResult] = []

    # -- lifecycle ---------------------------------------------------------
    def close(self) -> None:
        if self._owns_client:
            self._client.close()

    def __enter__(self) -> HttpClient:
        return self

    def __exit__(self, *exc_info: object) -> None:
        self.close()

    # -- internals ---------------------------------------------------------
    def _backoff(self, retry_number: int) -> float:
        """Delay before retry ``retry_number`` (1-based), capped and deterministic."""
        raw = self._settings.retry_base_delay * (2 ** (retry_number - 1))
        return min(raw, self._settings.retry_max_delay)

    def _resolve_delay(self, retry_number: int, retry_after: float | None) -> float:
        backoff = self._backoff(retry_number)
        if retry_after is None:
            return backoff
        # Honour the server, but never wait longer than the configured ceiling.
        return min(max(retry_after, 0.0), max(backoff, self._settings.retry_max_delay))

    @staticmethod
    def _headers(accept: str | None) -> dict[str, str]:
        return {"Accept": accept} if accept else {}

    # -- primary entry point ----------------------------------------------
    def get(
        self,
        url: str,
        *,
        source: str,
        accept: str | None = "application/json",
    ) -> FetchResult:
        """GET ``url``, retrying only the failures that are worth retrying."""
        attempts: list[AttemptRecord] = []
        last_error: SourceError | None = None

        for attempt in range(1, self._settings.max_attempts + 1):
            delay = 0.0
            try:
                response = self._client.get(url, headers=self._headers(accept))
            except httpx.TimeoutException as exc:
                record = AttemptRecord(
                    attempt=attempt,
                    url=url,
                    status=None,
                    outcome="timeout",
                    delay_before_seconds=0.0,
                    detail=f"{type(exc).__name__}: {exc}",
                )
                attempts.append(record)
                self.attempts.append(record)
                last_error = SourceError(
                    ErrorCode.NETWORK_ERROR,
                    f"request to {url} timed out after {self._settings.timeout_seconds}s",
                    source=source,
                    details={"url": url, "attempt": attempt, "cause": type(exc).__name__},
                )
                logger.warning(
                    "timeout (attempt %d/%d) %s",
                    attempt,
                    self._settings.max_attempts,
                    url,
                )
                if attempt < self._settings.max_attempts:
                    delay = self._backoff(attempt)
                    self._sleep(delay)
                    continue
                raise last_error from exc
            except httpx.HTTPError as exc:
                record = AttemptRecord(
                    attempt=attempt,
                    url=url,
                    status=None,
                    outcome="transport_error",
                    detail=f"{type(exc).__name__}: {exc}",
                )
                attempts.append(record)
                self.attempts.append(record)
                last_error = SourceError(
                    ErrorCode.NETWORK_ERROR,
                    f"transport error talking to {url}: {type(exc).__name__}",
                    source=source,
                    details={"url": url, "attempt": attempt, "cause": str(exc)},
                )
                logger.warning("transport error (attempt %d) %s: %s", attempt, url, exc)
                if attempt < self._settings.max_attempts:
                    delay = self._backoff(attempt)
                    self._sleep(delay)
                    continue
                raise last_error from exc

            status = response.status_code
            if 200 <= status < 300:
                content = response.content
                record = AttemptRecord(attempt=attempt, url=url, status=status, outcome="success")
                attempts.append(record)
                self.attempts.append(record)
                if not content:
                    raise SourceError(
                        ErrorCode.INVALID_API_RESPONSE,
                        f"{source} returned an empty body (HTTP {status}) for {url}",
                        source=source,
                        http_status=status,
                        details={"url": url},
                    )
                provenance = Provenance(
                    source=source,
                    url=str(response.url),
                    fetched_at=self._clock(),
                    http_status=status,
                    sha256=sha256_hex(content),
                    size_bytes=len(content),
                    content_type=response.headers.get("content-type"),
                    attempts=attempt,
                )
                logger.info(
                    "fetched %s (%d bytes, sha256=%s…)",
                    url,
                    len(content),
                    provenance.sha256[:12],
                )
                result = FetchResult(content=content, provenance=provenance, attempts=attempts)
                self.responses.append(result)
                return result

            retry_after = parse_retry_after(response.headers.get("retry-after"), now=self._clock())

            if status in RETRYABLE_STATUS:
                outcome = "rate_limited" if status == 429 else "server_error"
                record = AttemptRecord(
                    attempt=attempt,
                    url=url,
                    status=status,
                    outcome=outcome,
                    retry_after_seconds=retry_after,
                )
                attempts.append(record)
                self.attempts.append(record)
                if status == 429:
                    last_error = SourceError(
                        ErrorCode.RATE_LIMITED,
                        f"{source} rate-limited the request (HTTP 429) after {attempt} attempt(s)",
                        source=source,
                        http_status=status,
                        details={"url": url, "retry_after": retry_after},
                    )
                else:
                    last_error = SourceError(
                        ErrorCode.SOURCE_UNAVAILABLE,
                        f"{source} returned HTTP {status} after {attempt} attempt(s)",
                        source=source,
                        http_status=status,
                        details={"url": url},
                    )
                logger.warning(
                    "HTTP %d (attempt %d/%d) %s",
                    status,
                    attempt,
                    self._settings.max_attempts,
                    url,
                )
                if attempt < self._settings.max_attempts:
                    delay = self._resolve_delay(attempt, retry_after)
                    self._sleep(delay)
                    continue
                raise last_error

            # Non-retryable client errors: fail fast, and say why.
            code = self._client_error_code(status)
            message = self._client_error_message(source, status, url)
            record = AttemptRecord(
                attempt=attempt, url=url, status=status, outcome="client_error", detail=message
            )
            attempts.append(record)
            self.attempts.append(record)
            raise SourceError(
                code,
                message,
                source=source,
                http_status=status,
                details={"url": url, "body_preview": self._preview(response)},
            )

        # Defensive: the loop always either returns or raises.
        raise last_error or SourceError(
            ErrorCode.INTERNAL_ERROR, f"retry loop exhausted for {url}", source=source
        )

    @staticmethod
    def _preview(response: httpx.Response, limit: int = 300) -> str:
        try:
            text = response.text
        except Exception:  # pragma: no cover - defensive
            return "<unreadable body>"
        return text[:limit]

    @staticmethod
    def _client_error_code(status: int) -> ErrorCode:
        if status in {404, 410}:
            return ErrorCode.SOURCE_NOT_FOUND
        if status in {401, 403}:
            return ErrorCode.SOURCE_UNAVAILABLE
        return ErrorCode.INVALID_API_RESPONSE

    @staticmethod
    def _client_error_message(source: str, status: int, url: str) -> str:
        if status in {404, 410}:
            return f"{source} has no record at {url} (HTTP {status})"
        if status in {401, 403}:
            return (
                f"{source} refused the request (HTTP {status}); the endpoint may require "
                "a different client identity or the resource may be restricted"
            )
        return f"{source} rejected the request (HTTP {status}) for {url}"

    # -- file downloads ----------------------------------------------------
    def download(
        self,
        url: str,
        destination: Path,
        *,
        source: str,
        kind: str,
        accept: str | None = None,
    ) -> tuple[ArtifactRecord, FetchResult]:
        """Fetch ``url`` and write it atomically to ``destination``."""
        result = self.get(url, source=source, accept=accept)
        destination.parent.mkdir(parents=True, exist_ok=True)
        temporary = destination.with_name(destination.name + ".partial")
        temporary.write_bytes(result.content)
        temporary.replace(destination)
        artifact = ArtifactRecord(
            kind=kind,
            path=destination.name if destination.parent else str(destination),
            url=result.provenance.url,
            sha256=result.provenance.sha256,
            size_bytes=result.provenance.size_bytes,
        )
        return artifact, result


__all__ = [
    "AttemptRecord",
    "FetchResult",
    "HttpClient",
    "parse_retry_after",
    "sha256_hex",
    "utcnow",
]
