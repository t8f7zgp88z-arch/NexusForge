"""RCSB client: entry/entity parsing and chain↔accession enforcement."""

from __future__ import annotations

import copy

import httpx
import pytest

from nexusforge.clients.rcsb import RcsbClient, validate_chain_id, validate_pdb_id
from nexusforge.errors import ErrorCode, InputError, SourceError, ValidationError
from nexusforge.http import HttpClient

from ..conftest import fixture_bytes, fixture_json

ENTRY = fixture_json("rcsb_4HHB_entry.json")
ENTITY_1 = fixture_json("rcsb_4HHB_entity_1.json")
ENTITY_2 = fixture_json("rcsb_4HHB_entity_2.json")


def _client(handler, settings) -> RcsbClient:
    http = HttpClient(settings, client=httpx.Client(transport=httpx.MockTransport(handler)))
    return RcsbClient(http)


def test_parse_real_entry() -> None:
    entry = RcsbClient.parse_entry(ENTRY, requested_pdb_id="4HHB")
    assert entry.pdb_id == "4HHB"
    assert "DEOXYHAEMOGLOBIN" in entry.title.upper()
    assert entry.experimental_methods == ["X-RAY DIFFRACTION"]
    assert entry.resolution_angstrom == pytest.approx(1.74)
    assert entry.polymer_entity_ids == ["1", "2"]
    assert entry.deposited_atom_count == 4779
    assert entry.released_at is not None


def test_parse_entry_rejects_wrong_pdb_id() -> None:
    with pytest.raises(SourceError) as excinfo:
        RcsbClient.parse_entry(ENTRY, requested_pdb_id="1A2B")
    assert excinfo.value.code is ErrorCode.INVALID_API_RESPONSE


def test_parse_entity_maps_chain_and_uniprot() -> None:
    entity = RcsbClient.parse_polymer_entity(ENTITY_1, requested_pdb_id="4HHB", entity_id="1")
    assert entity.entity_id == "1"
    assert entity.auth_asym_ids == ["A", "C"]
    assert entity.uniprot_accessions == ["P69905"]
    assert entity.sequence.startswith("VLSPADKTNV")
    assert len(entity.sequence) == 141
    assert entity.coverage_of_reference == pytest.approx(0.9929577464788732)


def test_resolve_chain_success(settings) -> None:
    client = _client(
        lambda request: httpx.Response(
            200,
            json=(
                ENTITY_1
                if str(request.url).endswith("/1")
                else ENTITY_2
                if str(request.url).endswith("/2")
                else ENTRY
            ),
        ),
        settings,
    )
    mapping, provenance = client.resolve_chain("4HHB", "A", "P69905")
    assert mapping.auth_chain_id == "A"
    assert mapping.entity_id == "1"
    assert mapping.uniprot_accessions == ["P69905"]
    # entry metadata plus every polymer entity inspected (1 and 2 for 4HHB)
    assert len(provenance) == 3


def test_resolve_chain_rejects_wrong_chain_mapped_to_other_protein(settings) -> None:
    """Chain B of 4HHB is the beta chain (P68871), not alpha (P69905)."""

    def handler(request: httpx.Request) -> httpx.Response:
        url = str(request.url)
        if "/entry/" in url:
            return httpx.Response(200, json=ENTRY)
        if url.endswith("/1"):
            return httpx.Response(200, json=ENTITY_1)
        return httpx.Response(200, json=ENTITY_2)

    client = _client(handler, settings)
    with pytest.raises(ValidationError) as excinfo:
        client.resolve_chain("4HHB", "B", "P69905")
    assert excinfo.value.code is ErrorCode.SEQUENCE_MISMATCH
    assert "P68871" in excinfo.value.message


def test_resolve_chain_unknown_chain_is_explicit(settings) -> None:
    client = _client(
        lambda request: httpx.Response(
            200,
            json=(
                ENTITY_1
                if str(request.url).endswith("/1")
                else ENTITY_2
                if str(request.url).endswith("/2")
                else ENTRY
            ),
        ),
        settings,
    )
    with pytest.raises(ValidationError) as excinfo:
        client.resolve_chain("4HHB", "Z", "P69905")
    assert excinfo.value.code is ErrorCode.SEQUENCE_MISMATCH
    assert "polymer chains" in excinfo.value.message


def test_resolve_chain_refuses_ambiguous_entity(settings) -> None:
    ambiguous_entry = copy.deepcopy(ENTRY)
    ambiguous_entity = copy.deepcopy(ENTITY_1)
    ambiguous_entity["rcsb_polymer_entity_container_identifiers"]["entity_id"] = "2"

    client = _client(
        lambda request: httpx.Response(
            200, json=ambiguous_entity if "polymer_entity" in str(request.url) else ambiguous_entry
        ),
        settings,
    )
    with pytest.raises(ValidationError) as excinfo:
        client.resolve_chain("4HHB", "A", "P69905")
    assert excinfo.value.code is ErrorCode.CHAIN_MAPPING_AMBIGUOUS


def test_download_structure_writes_cif(settings, tmp_path) -> None:
    cif = fixture_bytes("4HHB.cif")

    def handler(request: httpx.Request) -> httpx.Response:
        return httpx.Response(200, content=cif, headers={"content-type": "chemical/x-mmcif"})

    client = _client(handler, settings)
    destination = tmp_path / "4HHB.cif"
    artifact, result = client.download_structure("4HHB", destination)
    assert destination.exists()
    assert destination.read_bytes() == cif
    assert artifact.sha256 == result.provenance.sha256
    assert len(artifact.sha256) == 64


def test_identifier_validation() -> None:
    assert validate_pdb_id("4hhb") == "4HHB"
    with pytest.raises(InputError):
        validate_pdb_id("HHB4")
    assert validate_chain_id("a") == "A"
    with pytest.raises(InputError):
        validate_chain_id("chain")
