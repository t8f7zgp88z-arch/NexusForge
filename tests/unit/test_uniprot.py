"""UniProt client parsing and validation."""

from __future__ import annotations

import pytest

from nexusforge.clients.uniprot import UniProtClient, validate_accession
from nexusforge.errors import ErrorCode, InputError, SourceError

from ..conftest import fixture_json

PAYLOAD = fixture_json("uniprot_P69905.json")


def test_parse_real_uniprot_entry() -> None:
    record = UniProtClient.parse(PAYLOAD, requested_accession="P69905")
    assert record.accession == "P69905"
    assert record.entry_name == "HBA_HUMAN"
    assert record.reviewed is True
    assert record.protein_name == "Hemoglobin subunit alpha"
    assert record.organism == "Homo sapiens"
    assert record.taxon_id == 9606
    assert record.sequence_length == 142
    assert len(record.sequence) == 142
    assert record.sequence.startswith("MVLSPADKTNVKAAWGKVGA")
    assert record.sequence.endswith("VSTVLTSKYR")
    assert "HBA1" in record.gene_names


@pytest.mark.parametrize(
    "accession",
    ["P69905", "p69905", " A0A123B4C5 ", "Q8N158"],
)
def test_validate_accession_accepts_realistic_ids(accession: str) -> None:
    assert validate_accession(accession) == accession.strip().upper()


@pytest.mark.parametrize("accession", ["", "P6990", "not-an-accession", "P6990!", "12345"])
def test_validate_accession_rejects_bad_ids(accession: str) -> None:
    with pytest.raises(InputError) as excinfo:
        validate_accession(accession)
    assert excinfo.value.code is ErrorCode.INVALID_INPUT
    assert excinfo.value.stage.value == "input_validation"


def test_parse_rejects_entry_for_another_accession() -> None:
    with pytest.raises(SourceError) as excinfo:
        UniProtClient.parse(PAYLOAD, requested_accession="Q8N158")
    assert excinfo.value.code is ErrorCode.INVALID_API_RESPONSE
    assert "Q8N158" in excinfo.value.message


def test_parse_rejects_search_payload() -> None:
    with pytest.raises(SourceError) as excinfo:
        UniProtClient.parse({"results": []}, requested_accession="P69905")
    assert excinfo.value.code is ErrorCode.INVALID_API_RESPONSE


def test_parse_rejects_sequence_length_disagreement() -> None:
    broken = {
        "primaryAccession": "P69905",
        "entryType": "UniProtKB reviewed (Swiss-Prot)",
        "uniProtkbId": "HBA_HUMAN",
        "organism": {"scientificName": "Homo sapiens"},
        "sequence": {"value": "MVLSPADKTNV", "length": 999},
    }
    with pytest.raises(SourceError) as excinfo:
        UniProtClient.parse(broken, requested_accession="P69905")
    assert excinfo.value.code is ErrorCode.INVALID_API_RESPONSE
    assert "disagrees" in excinfo.value.message


def test_fetch_uses_the_real_endpoint(offline_http) -> None:
    client = UniProtClient(offline_http)
    record, provenance = client.fetch("p69905")
    assert record.accession == "P69905"
    assert provenance.source == "uniprot"
    assert provenance.http_status == 200
    assert provenance.url.endswith("/uniprotkb/P69905.json")
    assert len(provenance.sha256) == 64


def test_fetch_unknown_accession_returns_source_not_found(offline_http) -> None:
    client = UniProtClient(offline_http)
    with pytest.raises(SourceError) as excinfo:
        client.fetch("Q9NQ11")
    assert excinfo.value.code is ErrorCode.SOURCE_NOT_FOUND
    assert excinfo.value.http_status == 404
