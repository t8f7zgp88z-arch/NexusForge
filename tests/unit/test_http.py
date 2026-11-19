"""HTTP retry policy and failure handling."""

from __future__ import annotations

from datetime import UTC, datetime, timedelta

import httpx
import pytest

from nexusforge.clients.common import parse_json_object
from nexusforge.config import load_settings
from nexusforge.errors import ErrorCode, SourceError
from nexusforge.http import HttpClient, parse_retry_after

from ..conftest import fixture_bytes


def _http(handler, *, max_attempts: int = 3, retry_max_delay: float = 0.0):
    settings = load_settings(
        env={},
        retry_base_delay=0.25,
        retry_max_delay=retry_max_delay,
        max_attempts=max_attempts,
        timeout_seconds=1.0,
    )
    sleeps: list[float] = []
    client = httpx.Client(transport=httpx.MockTransport(handler))
    http = HttpClient(settings, client=client, sleep=sleeps.append)
    return http, sleeps, client


def test_success_records_provenance() -> None:
    http, sleeps, client = _http(lambda request: httpx.Response(200, content=b'{"a": 1}'))
    result = http.get("https://example.org/x", source="example")
    assert result.content == b'{"a": 1}'
    assert result.provenance.sha256 == __import__("hashlib").sha256(b'{"a": 1}').hexdigest()
    assert result.provenance.attempts == 1
    assert result.provenance.http_status == 200
    assert sleeps == []
    client.close()


def test_429_is_retried_and_retry_after_is_honoured() -> None:
    calls = {"n": 0}

    def handler(request: httpx.Request) -> httpx.Response:
        calls["n"] += 1
        if calls["n"] == 1:
            return httpx.Response(429, headers={"retry-after": "0"})
        return httpx.Response(200, content=b"ok")

    http, sleeps, client = _http(handler)
    result = http.get("https://example.org/limited", source="example")
    assert result.content == b"ok"
    assert calls["n"] == 2
    assert sleeps == [0.0]
    assert result.provenance.attempts == 2
    client.close()


def test_retry_after_is_capped_by_the_configured_ceiling() -> None:
    calls = {"n": 0}

    def handler(request: httpx.Request) -> httpx.Response:
        calls["n"] += 1
        if calls["n"] == 1:
            return httpx.Response(429, headers={"retry-after": "3600"})
        return httpx.Response(200, content=b"ok")

    http, sleeps, client = _http(handler, retry_max_delay=2.0)
    http.get("https://example.org/limited", source="example")
    assert sleeps == [2.0]
    client.close()


def test_persistent_429_raises_rate_limited() -> None:
    http, _sleeps, client = _http(lambda request: httpx.Response(429), max_attempts=3)
    with pytest.raises(SourceError) as excinfo:
        http.get("https://example.org/limited", source="example")
    assert excinfo.value.code is ErrorCode.RATE_LIMITED
    assert excinfo.value.retryable is True
    assert excinfo.value.exit_code == 4
    client.close()


def test_5xx_is_retried_then_succeeds() -> None:
    calls = {"n": 0}

    def handler(request: httpx.Request) -> httpx.Response:
        calls["n"] += 1
        if calls["n"] < 3:
            return httpx.Response(503)
        return httpx.Response(200, content=b"recovered")

    http, sleeps, client = _http(handler, max_attempts=3, retry_max_delay=5.0)
    result = http.get("https://example.org/flaky", source="example")
    assert result.content == b"recovered"
    assert calls["n"] == 3
    assert len(sleeps) == 2
    assert sleeps[0] == pytest.approx(0.25)
    assert sleeps[1] == pytest.approx(0.5)
    client.close()


def test_persistent_5xx_raises_source_unavailable() -> None:
    http, _sleeps, client = _http(lambda request: httpx.Response(500), max_attempts=2)
    with pytest.raises(SourceError) as excinfo:
        http.get("https://example.org/down", source="example")
    assert excinfo.value.code is ErrorCode.SOURCE_UNAVAILABLE
    assert excinfo.value.http_status == 500
    client.close()


def test_404_fails_fast_without_retry() -> None:
    calls = {"n": 0}

    def handler(request: httpx.Request) -> httpx.Response:
        calls["n"] += 1
        return httpx.Response(404)

    http, sleeps, client = _http(handler)
    with pytest.raises(SourceError) as excinfo:
        http.get("https://example.org/missing", source="example")
    assert excinfo.value.code is ErrorCode.SOURCE_NOT_FOUND
    assert excinfo.value.exit_code == 3
    assert calls["n"] == 1
    assert sleeps == []
    client.close()


def test_timeout_is_retried_then_reported() -> None:
    calls = {"n": 0}

    def handler(request: httpx.Request) -> httpx.Response:
        calls["n"] += 1
        raise httpx.ReadTimeout("simulated timeout", request=request)

    http, sleeps, client = _http(handler, max_attempts=2)
    with pytest.raises(SourceError) as excinfo:
        http.get("https://example.org/slow", source="example")
    assert excinfo.value.code is ErrorCode.NETWORK_ERROR
    assert excinfo.value.retryable is True
    assert calls["n"] == 2
    assert len(sleeps) == 1
    client.close()


def test_timeout_then_success() -> None:
    calls = {"n": 0}

    def handler(request: httpx.Request) -> httpx.Response:
        calls["n"] += 1
        if calls["n"] == 1:
            raise httpx.ConnectTimeout("simulated", request=request)
        return httpx.Response(200, content=b"ok")

    http, _sleeps, client = _http(handler, max_attempts=3)
    assert http.get("https://example.org/slow", source="example").content == b"ok"
    client.close()


def test_empty_body_is_an_invalid_response() -> None:
    http, _sleeps, client = _http(lambda request: httpx.Response(200, content=b""))
    with pytest.raises(SourceError) as excinfo:
        http.get("https://example.org/empty", source="example")
    assert excinfo.value.code is ErrorCode.INVALID_API_RESPONSE
    assert "empty body" in excinfo.value.message
    client.close()


def test_invalid_json_is_reported_with_context() -> None:
    payload = b"<html>not json</html>"
    with pytest.raises(SourceError) as excinfo:
        parse_json_object(payload, source="example", url="https://example.org/x")
    assert excinfo.value.code is ErrorCode.INVALID_API_RESPONSE
    assert excinfo.value.details["body_preview"].startswith("<html>")


def test_403_is_not_retried_and_says_so(offline_http) -> None:
    """AlphaFold's file endpoint returned 403 without a client User-Agent."""
    http, _sleeps, client = _http(lambda request: httpx.Response(403))
    with pytest.raises(SourceError) as excinfo:
        http.get("https://example.org/forbidden", source="example")
    assert excinfo.value.code is ErrorCode.SOURCE_UNAVAILABLE
    assert "refused" in excinfo.value.message
    client.close()


def test_retry_after_parsing() -> None:
    assert parse_retry_after("5") == pytest.approx(5.0)
    assert parse_retry_after("0") == pytest.approx(0.0)
    assert parse_retry_after(None) is None
    assert parse_retry_after("not-a-date") is None
    now = datetime(2026, 11, 19, 12, 0, 0, tzinfo=UTC)
    header = (now + timedelta(seconds=30)).strftime("%a, %d %b %Y %H:%M:%S GMT")
    assert parse_retry_after(header, now=now) == pytest.approx(30.0)


def test_download_writes_bytes_and_hashes(tmp_path) -> None:
    content = fixture_bytes("alphafold_P69905.json")
    http, _sleeps, client = _http(lambda request: httpx.Response(200, content=content))
    destination = tmp_path / "prediction.json"
    artifact, result = http.download(
        "https://example.org/prediction", destination, source="example", kind="prediction"
    )
    assert destination.read_bytes() == content
    assert artifact.sha256 == result.provenance.sha256
    assert artifact.size_bytes == len(content)
    client.close()
