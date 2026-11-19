"""AlphaFold DB client: schema handling, selection and confidence validation."""

from __future__ import annotations

import copy

import pytest

from nexusforge.clients.alphafold import AlphaFoldClient
from nexusforge.errors import ConfidenceError, ErrorCode, SourceError

from ..conftest import fixture_json

PAYLOAD = fixture_json("alphafold_P69905.json")
CONFIDENCE = fixture_json("alphafold_plddt_P69905.json")
PAE = fixture_json("alphafold_pae_P69905.json")


def test_parse_real_prediction() -> None:
    prediction = AlphaFoldClient.parse_prediction(PAYLOAD, accession="P69905")
    assert prediction.entry_id == "AF-P69905-F1"
    assert prediction.uniprot_accession == "P69905"
    assert prediction.latest_version == 6
    assert prediction.tool_used == "AlphaFold Monomer v2.0 pipeline"
    assert len(prediction.sequence) == 142
    assert prediction.sequence_start == 1
    assert prediction.sequence_end == 142
    assert prediction.chain_id == "A"
    assert prediction.global_plddt == pytest.approx(98.06)
    assert prediction.files.cif_url and prediction.files.cif_url.endswith(".cif")
    assert prediction.files.plddt_doc_url
    assert prediction.files.pae_doc_url
    # No fallback was needed with the live schema.
    assert prediction.schema_notes == []


def test_parse_uses_current_fields_not_hardcoded_urls() -> None:
    payload = copy.deepcopy(PAYLOAD)
    payload[0]["cifUrl"] = "https://example.org/custom/model.cif"
    payload[0]["latestVersion"] = 7
    prediction = AlphaFoldClient.parse_prediction(payload, accession="P69905")
    assert prediction.files.cif_url == "https://example.org/custom/model.cif"
    assert prediction.latest_version == 7


def test_parse_records_fallback_field_names() -> None:
    legacy = copy.deepcopy(PAYLOAD)
    entry = legacy[0]
    entry.pop("modelEntityId")
    entry["entryId"] = "AF-P69905-F1"
    entry.pop("sequence")
    entry["uniprotSequence"] = "M" * 142
    prediction = AlphaFoldClient.parse_prediction(legacy, accession="P69905")
    assert prediction.entry_id == "AF-P69905-F1"
    assert len(prediction.sequence) == 142
    assert any("entryId" in note for note in prediction.schema_notes)
    assert any("uniprotSequence" in note for note in prediction.schema_notes)


def test_select_entry_prefers_highest_version() -> None:
    payload = copy.deepcopy(PAYLOAD)
    older = copy.deepcopy(PAYLOAD[0])
    older["latestVersion"] = 3
    older["modelEntityId"] = "AF-P69905-F1-v3"
    payload.append(older)
    prediction = AlphaFoldClient.parse_prediction(payload, accession="P69905")
    assert prediction.latest_version == 6


def test_select_entry_rejects_mismatched_accession() -> None:
    with pytest.raises(SourceError) as excinfo:
        AlphaFoldClient.select_entry([{"uniprotAccession": "Q8N158"}], accession="P69905")
    assert excinfo.value.code is ErrorCode.SOURCE_NOT_FOUND


def test_empty_prediction_list_is_not_found() -> None:
    with pytest.raises(SourceError) as excinfo:
        AlphaFoldClient.select_entry([], accession="P69905")
    assert excinfo.value.code is ErrorCode.SOURCE_NOT_FOUND


def test_prediction_without_structure_file_is_rejected() -> None:
    payload = copy.deepcopy(PAYLOAD)
    payload[0].pop("cifUrl")
    payload[0].pop("pdbUrl")
    with pytest.raises(SourceError) as excinfo:
        AlphaFoldClient.parse_prediction(payload, accession="P69905")
    assert excinfo.value.code is ErrorCode.INVALID_API_RESPONSE
    assert "structure file" in excinfo.value.message


def test_confidence_document_parses_real_data() -> None:
    scores = AlphaFoldClient.parse_confidence(CONFIDENCE)
    assert scores.source == "alphafold_confidence_json"
    assert len(scores.scores) == 142
    assert scores.residue_numbers[0] == 1
    assert min(scores.scores) >= 0.0
    assert max(scores.scores) <= 100.0


def test_confidence_length_mismatch_is_typed() -> None:
    broken = {
        "residueNumber": CONFIDENCE["residueNumber"],
        "confidenceScore": CONFIDENCE["confidenceScore"][:-1],
    }
    with pytest.raises(ConfidenceError) as excinfo:
        AlphaFoldClient.parse_confidence(broken)
    assert excinfo.value.code is ErrorCode.CONFIDENCE_LENGTH_MISMATCH


def test_confidence_out_of_range_is_typed() -> None:
    broken = copy.deepcopy(CONFIDENCE)
    broken["confidenceScore"][5] = 123.0
    with pytest.raises(SourceError) as excinfo:
        AlphaFoldClient.parse_confidence(broken)
    assert excinfo.value.code is ErrorCode.INVALID_API_RESPONSE
    assert "outside [0, 100]" in excinfo.value.message


def test_confidence_non_monotonic_residues_rejected() -> None:
    broken = {
        "residueNumber": [1, 3, 2],
        "confidenceScore": [90.0, 90.0, 90.0],
    }
    with pytest.raises(SourceError) as excinfo:
        AlphaFoldClient.parse_confidence(broken)
    assert "strictly increasing" in excinfo.value.message


def test_pae_summary_from_real_document() -> None:
    summary = AlphaFoldClient.parse_pae(PAE)
    assert summary.shape == (142, 142)
    assert summary.mean_predicted_aligned_error > 0
    assert summary.max_predicted_aligned_error == pytest.approx(31.75)
    assert 0 < summary.p90_predicted_aligned_error <= summary.max_predicted_aligned_error


def test_pae_non_square_matrix_rejected() -> None:
    with pytest.raises(SourceError) as excinfo:
        AlphaFoldClient.parse_pae([{"predicted_aligned_error": [[1.0, 2.0]]}])
    assert "square" in excinfo.value.message


def test_fetch_pae_404_fails_fast(offline_http) -> None:
    client = AlphaFoldClient(offline_http)
    with pytest.raises(SourceError) as excinfo:
        client.fetch_pae("https://alphafold.ebi.ac.uk/files/missing-pae.json")
    assert excinfo.value.code is ErrorCode.SOURCE_NOT_FOUND
