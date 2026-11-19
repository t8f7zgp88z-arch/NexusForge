"""Shared fixtures: real captured API responses served by an offline transport."""

from __future__ import annotations

import gzip
import json
from collections.abc import Callable
from pathlib import Path
from typing import Any

import httpx
import pytest

from nexusforge.config import Settings, load_settings
from nexusforge.http import HttpClient

FIXTURES = Path(__file__).parent / "fixtures"

#: Exact live URLs used by the v0.1 case, mapped to captured response files.
FIXTURE_ROUTES: dict[str, str] = {
    "https://rest.uniprot.org/uniprotkb/P69905.json": "uniprot_P69905.json",
    "https://alphafold.ebi.ac.uk/api/prediction/P69905": "alphafold_P69905.json",
    "https://alphafold.ebi.ac.uk/files/AF-P69905-F1-model_v6.cif": "AF-P69905-F1-model_v6.cif",
    "https://alphafold.ebi.ac.uk/files/AF-P69905-F1-confidence_v6.json": (
        "alphafold_plddt_P69905.json"
    ),
    "https://alphafold.ebi.ac.uk/files/AF-P69905-F1-predicted_aligned_error_v6.json": (
        "alphafold_pae_P69905.json"
    ),
    "https://data.rcsb.org/rest/v1/core/entry/4HHB": "rcsb_4HHB_entry.json",
    "https://data.rcsb.org/rest/v1/core/polymer_entity/4HHB/1": "rcsb_4HHB_entity_1.json",
    "https://data.rcsb.org/rest/v1/core/polymer_entity/4HHB/2": "rcsb_4HHB_entity_2.json",
    "https://files.rcsb.org/download/4HHB.cif": "4HHB.cif",
}


def fixture_bytes(name: str) -> bytes:
    """Decompress a captured response (fixtures are stored gzipped)."""
    return gzip.decompress((FIXTURES / f"{name}.gz").read_bytes())


def fixture_json(name: str) -> Any:
    return json.loads(fixture_bytes(name))


def offline_handler(
    overrides: dict[str, Callable[[httpx.Request], httpx.Response]] | None = None,
) -> Callable[[httpx.Request], httpx.Response]:
    """Serve captured responses; unknown URLs 404 like the real APIs would."""
    overrides = overrides or {}

    def handler(request: httpx.Request) -> httpx.Response:
        url = str(request.url)
        if url in overrides:
            return overrides[url](request)
        name = FIXTURE_ROUTES.get(url)
        if name is None:
            return httpx.Response(404, json={"error": f"no fixture for {url}"})
        return httpx.Response(
            200,
            content=fixture_bytes(name),
            headers={"content-type": "application/json"},
        )

    return handler


@pytest.fixture
def settings() -> Settings:
    return load_settings(
        env={}, retry_base_delay=0.0, retry_max_delay=0.0, max_attempts=2, timeout_seconds=5.0
    )


@pytest.fixture
def offline_http(settings: Settings) -> HttpClient:
    client = httpx.Client(transport=httpx.MockTransport(offline_handler()))
    http = HttpClient(settings, client=client)
    yield http
    client.close()
