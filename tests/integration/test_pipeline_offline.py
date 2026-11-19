"""Offline end-to-end integration tests driven by captured API responses."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import httpx
import pytest

from nexusforge.config import load_settings
from nexusforge.errors import ErrorCode, InputError
from nexusforge.http import HttpClient
from nexusforge.models import RunStatus
from nexusforge.pipeline.protein_smoke import ProteinSmokePipeline, SmokeRequest

from ..conftest import offline_handler


@pytest.fixture
def pipeline():
    settings = load_settings(
        env={}, retry_base_delay=0.0, retry_max_delay=0.0, max_attempts=2, timeout_seconds=5.0
    )
    client = httpx.Client(transport=httpx.MockTransport(offline_handler()))
    http = HttpClient(settings, client=client)
    instance = ProteinSmokePipeline(settings=settings, http=http)
    yield instance, client
    instance.close()
    client.close()


def test_full_pipeline_succeeds_offline(pipeline, tmp_path: Path) -> None:
    instance, _client = pipeline
    run_dir = tmp_path / "P69905-4HHB-A"
    outcome = instance.run(SmokeRequest("P69905", "4HHB", "A", run_dir))

    assert outcome.error is None
    assert outcome.exit_code == 0
    assert outcome.manifest.status is RunStatus.SUCCESS
    assert outcome.manifest.reproducible is True

    metrics = outcome.manifest.metrics
    sequence = metrics["sequence"]
    assert sequence["alphafold_identity"] == 1.0
    assert sequence["pdb_entity_identity"] == 1.0
    assert sequence["pdb_entity_coverage_query"] == pytest.approx(141 / 142)
    assert sequence["pdb_entity_gap_openings"] == 1

    structure = metrics["structure"]
    assert structure["mapped_residue_count"] == 141
    assert 0.0 < structure["rmsd_angstrom"] < 2.0
    assert 0.85 <= structure["lddt_ca"] <= 1.0
    assert structure["lddt_ca_pairs"] > 1000
    assert structure["lddt_ca_thresholds"] == [0.5, 1.0, 2.0, 4.0]

    confidence = metrics["confidence"]
    assert confidence["residue_count"] == 142
    assert 0.0 <= confidence["min_plddt"] <= 100.0
    assert confidence["source"] == "alphafold_confidence_json"

    # every request is recorded as a source with a digest (9 fetches for this case)
    assert len(outcome.manifest.sources) == 9
    assert all(source.url.startswith("https://") for source in outcome.manifest.sources)
    assert all(len(source.sha256) == 64 for source in outcome.manifest.sources)

    assert "pae" in metrics
    assert metrics["pae"]["shape"] == [142, 142]
    assert any(item["metric"] == "TM-score" for item in metrics["not_evaluated"])

    # artifacts exist and every recorded digest matches the bytes on disk
    assert (run_dir / "manifest.json").exists()
    assert (run_dir / "report.json").exists()
    assert (run_dir / "report.md").exists()
    assert (run_dir / "metrics.json").exists()
    assert (run_dir / "logs" / "run.log").exists()
    for artifact in outcome.manifest.artifacts:
        path = run_dir / artifact.path
        assert path.exists(), artifact.path
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        assert digest == artifact.sha256, artifact.path

    report = json.loads((run_dir / "report.json").read_text())
    assert report["status"] == "SUCCESS"
    assert report["metrics"]["structure"]["lddt_ca"] == structure["lddt_ca"]
    assert (
        report["manifest_sha256"]
        == hashlib.sha256((run_dir / "manifest.json").read_bytes()).hexdigest()
    )

    markdown = (run_dir / "report.md").read_text()
    assert "lDDT-Cα" in markdown
    assert "P69905" in markdown and "4HHB" in markdown
    assert "not a full-atom lDDT" in markdown


def test_second_run_refuses_to_overwrite(pipeline, tmp_path: Path) -> None:
    instance, _client = pipeline
    run_dir = tmp_path / "run"
    instance.run(SmokeRequest("P69905", "4HHB", "A", run_dir))
    with pytest.raises(InputError) as excinfo:
        instance.run(SmokeRequest("P69905", "4HHB", "A", run_dir))
    assert excinfo.value.code is ErrorCode.INVALID_INPUT
    assert "--force" in excinfo.value.message


def test_force_overwrites_existing_directory(pipeline, tmp_path: Path) -> None:
    instance, _client = pipeline
    run_dir = tmp_path / "run"
    first = instance.run(SmokeRequest("P69905", "4HHB", "A", run_dir))
    second = instance.run(SmokeRequest("P69905", "4HHB", "A", run_dir, overwrite=True))
    assert first.error is None and second.error is None
    assert (run_dir / "manifest.json").exists()


def test_missing_pae_is_reported_not_invented(tmp_path: Path) -> None:
    settings = load_settings(
        env={}, retry_base_delay=0.0, retry_max_delay=0.0, max_attempts=2, timeout_seconds=5.0
    )
    overrides = {
        "https://alphafold.ebi.ac.uk/files/AF-P69905-F1-predicted_aligned_error_v6.json": (
            lambda request: httpx.Response(404)
        )
    }
    client = httpx.Client(transport=httpx.MockTransport(offline_handler(overrides)))
    instance = ProteinSmokePipeline(settings=settings, http=HttpClient(settings, client=client))
    try:
        outcome = instance.run(SmokeRequest("P69905", "4HHB", "A", tmp_path / "run"))
    finally:
        instance.close()
        client.close()

    assert outcome.error is None
    assert "pae" not in outcome.manifest.metrics
    not_evaluated = outcome.manifest.metrics["not_evaluated"]
    assert any(item["metric"] == "PAE analysis" for item in not_evaluated)
    assert any("PAE data was advertised" in warning for warning in outcome.manifest.warnings)
