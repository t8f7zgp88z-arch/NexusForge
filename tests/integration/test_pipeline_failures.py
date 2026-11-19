"""Failure injection: every injected fault must produce a typed, staged failure.

None of these runs may be reported as SUCCESS.
"""

from __future__ import annotations

import json
from pathlib import Path

import httpx
import pytest

from nexusforge.config import load_settings
from nexusforge.errors import ErrorCode
from nexusforge.http import HttpClient
from nexusforge.models import RunStatus
from nexusforge.pipeline.protein_smoke import ProteinSmokePipeline, SmokeRequest

from ..conftest import FIXTURE_ROUTES, fixture_bytes, offline_handler

CONFIDENCE_URL = "https://alphafold.ebi.ac.uk/files/AF-P69905-F1-confidence_v6.json"


def _run(tmp_path: Path, overrides=None):
    settings = load_settings(
        env={}, retry_base_delay=0.0, retry_max_delay=0.0, max_attempts=2, timeout_seconds=5.0
    )
    client = httpx.Client(transport=httpx.MockTransport(offline_handler(overrides)))
    pipeline = ProteinSmokePipeline(settings=settings, http=HttpClient(settings, client=client))
    try:
        return pipeline.run(SmokeRequest("P69905", "4HHB", "A", tmp_path / "run"))
    finally:
        pipeline.close()
        client.close()


def test_missing_uniprot_entry_fails_at_retrieval(tmp_path: Path) -> None:
    outcome = _run(
        tmp_path,
        {
            "https://rest.uniprot.org/uniprotkb/P69905.json": (
                lambda request: httpx.Response(404, text="not found")
            )
        },
    )
    assert outcome.error is not None
    assert outcome.manifest.status is RunStatus.FAILED
    assert outcome.error.code is ErrorCode.SOURCE_NOT_FOUND
    assert outcome.error.stage.value == "data_retrieval"
    assert outcome.exit_code == 3
    report = json.loads((tmp_path / "run" / "report.json").read_text())
    assert report["status"] == "FAILED"
    assert report["error"]["code"] == "SOURCE_NOT_FOUND"
    assert "structure" not in report["metrics"]
    assert report["metrics"].get("sequence") is None


def test_invalid_json_response_fails_with_context(tmp_path: Path) -> None:
    outcome = _run(
        tmp_path,
        {
            "https://rest.uniprot.org/uniprotkb/P69905.json": (
                lambda request: httpx.Response(200, text="<html>maintenance</html>")
            )
        },
    )
    assert outcome.error is not None
    assert outcome.error.code is ErrorCode.INVALID_API_RESPONSE
    assert outcome.error.stage.value == "data_retrieval"


def test_empty_response_fails_at_retrieval(tmp_path: Path) -> None:
    outcome = _run(
        tmp_path,
        {
            "https://rest.uniprot.org/uniprotkb/P69905.json": (
                lambda request: httpx.Response(200, content=b"")
            )
        },
    )
    assert outcome.error is not None
    assert outcome.error.code is ErrorCode.INVALID_API_RESPONSE
    assert "empty body" in outcome.error.message


def test_chain_belongs_to_another_protein(tmp_path: Path) -> None:
    """4HHB chain B is beta-globin (P68871); requesting it as P69905 must fail."""
    settings = load_settings(env={}, retry_base_delay=0.0, retry_max_delay=0.0, max_attempts=2)
    client = httpx.Client(transport=httpx.MockTransport(offline_handler()))
    pipeline = ProteinSmokePipeline(settings=settings, http=HttpClient(settings, client=client))
    try:
        outcome = pipeline.run(SmokeRequest("P69905", "4HHB", "B", tmp_path / "run"))
    finally:
        pipeline.close()
        client.close()
    assert outcome.error is not None
    assert outcome.manifest.status is RunStatus.FAILED
    assert outcome.error.code is ErrorCode.SEQUENCE_MISMATCH
    assert outcome.error.stage.value == "sequence_validation"
    assert "P68871" in outcome.error.message


def test_unknown_chain_is_not_silently_mapped(tmp_path: Path) -> None:
    settings = load_settings(env={}, retry_base_delay=0.0, retry_max_delay=0.0, max_attempts=2)
    client = httpx.Client(transport=httpx.MockTransport(offline_handler()))
    pipeline = ProteinSmokePipeline(settings=settings, http=HttpClient(settings, client=client))
    try:
        outcome = pipeline.run(SmokeRequest("P69905", "4HHB", "Z", tmp_path / "run"))
    finally:
        pipeline.close()
        client.close()
    assert outcome.error is not None
    assert outcome.error.code is ErrorCode.SEQUENCE_MISMATCH
    assert "polymer chains" in outcome.error.message


def test_confidence_length_mismatch_blocks_evaluation(tmp_path: Path) -> None:
    real = json.loads(fixture_bytes("alphafold_plddt_P69905.json"))
    truncated = {
        "residueNumber": real["residueNumber"][:-5],
        "confidenceScore": real["confidenceScore"][:-5],
        "confidenceCategory": real["confidenceCategory"][:-5],
    }
    outcome = _run(
        tmp_path,
        {CONFIDENCE_URL: (lambda request: httpx.Response(200, json=truncated))},
    )
    assert outcome.error is not None
    assert outcome.error.code is ErrorCode.CONFIDENCE_LENGTH_MISMATCH
    assert outcome.error.stage.value == "confidence_evaluation"
    assert outcome.manifest.status is RunStatus.FAILED


def test_out_of_range_plddt_is_rejected(tmp_path: Path) -> None:
    real = json.loads(fixture_bytes("alphafold_plddt_P69905.json"))
    real["confidenceScore"][3] = 140.0
    outcome = _run(
        tmp_path,
        {CONFIDENCE_URL: (lambda request: httpx.Response(200, json=real))},
    )
    assert outcome.error is not None
    assert outcome.error.code is ErrorCode.INVALID_API_RESPONSE
    assert "outside [0, 100]" in outcome.error.message


def test_failure_keeps_retrieved_artifacts_for_diagnosis(tmp_path: Path) -> None:
    """A late failure must not discard data already fetched."""
    real = json.loads(fixture_bytes("alphafold_plddt_P69905.json"))
    real["confidenceScore"] = real["confidenceScore"][:-3]
    real["residueNumber"] = real["residueNumber"][:-3]
    outcome = _run(
        tmp_path,
        {CONFIDENCE_URL: (lambda request: httpx.Response(200, json=real))},
    )
    run_dir = tmp_path / "run"
    assert outcome.error is not None
    assert (run_dir / "raw" / "uniprot__P69905.json").exists()
    assert (run_dir / "manifest.json").exists()
    manifest = json.loads((run_dir / "manifest.json").read_text())
    assert manifest["status"] == "FAILED"
    assert manifest["error"]["code"] == "CONFIDENCE_LENGTH_MISMATCH"
    assert manifest["sources"], "sources fetched before the failure must be recorded"


def test_unknown_url_in_fixture_map_is_a_404(tmp_path: Path) -> None:
    """Guard the fixture map itself: unknown URLs must fail, never silently pass."""
    with pytest.raises(KeyError):
        FIXTURE_ROUTES["https://example.org/not-in-map"]
