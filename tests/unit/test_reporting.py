"""Manifest consistency and report rendering."""

from __future__ import annotations

from datetime import UTC, datetime

import pytest

from nexusforge.config import load_settings
from nexusforge.errors import ErrorCode, NexusForgeError, ReportingError
from nexusforge.models import CheckResult, CheckStatus, RequestSpec, RunStatus
from nexusforge.reporting.json_report import build_report
from nexusforge.reporting.manifest import build_manifest, ensure_consistent, fail_manifest
from nexusforge.reporting.markdown_report import render


def _manifest():
    return build_manifest(
        request=RequestSpec(uniprot_accession="P69905", pdb_id="4HHB", chain_id="A"),
        settings=load_settings(env={}),
        run_id="P69905-4HHB-A",
        generated_at=datetime(2026, 11, 19, 12, 0, 0, tzinfo=UTC),
    )


def test_failed_manifest_requires_an_error() -> None:
    with pytest.raises(ReportingError):
        ensure_consistent(_manifest())


def test_success_manifest_must_not_carry_an_error() -> None:
    manifest = _manifest()
    manifest.status = RunStatus.SUCCESS
    manifest.reproducible = True
    manifest.metrics = {"sequence": {}, "structure": {}}
    manifest.error = NexusForgeError(ErrorCode.INVALID_SEQUENCE, "boom").to_dict()
    with pytest.raises(ReportingError):
        ensure_consistent(manifest)


def test_success_manifest_requires_metrics() -> None:
    manifest = _manifest()
    manifest.status = RunStatus.SUCCESS
    manifest.reproducible = True
    with pytest.raises(ReportingError):
        ensure_consistent(manifest)


def test_success_manifest_requires_reproducible_flag() -> None:
    manifest = _manifest()
    manifest.status = RunStatus.SUCCESS
    manifest.metrics = {"sequence": {}, "structure": {}}
    with pytest.raises(ReportingError):
        ensure_consistent(manifest)


def test_failed_manifest_renders_failure_and_no_metrics() -> None:
    manifest = _manifest()
    error = NexusForgeError(
        ErrorCode.SEQUENCE_MISMATCH,
        "chain B of 4HHB maps to UniProt P68871, not P69905",
        source="rcsb-entity",
        details={"chain": "B"},
    )
    fail_manifest(manifest, error)
    markdown = render(manifest)
    assert "Status: FAILED" in markdown
    assert "SEQUENCE_MISMATCH" in markdown
    assert "sequence_validation" in markdown
    assert "No scientific conclusion" in markdown
    assert "lDDT-Cα" not in markdown.split("## Failure")[0]


def test_report_for_failed_run_draws_no_conclusion() -> None:
    manifest = _manifest()
    fail_manifest(manifest, NexusForgeError(ErrorCode.SOURCE_NOT_FOUND, "no entry"))
    report = build_report(manifest)
    assert report["status"] == "FAILED"
    assert report["error"]["code"] == "SOURCE_NOT_FOUND"
    assert "No scientific conclusion" in report["conclusion"]
    assert report["metrics"] == {}


def test_success_report_contains_provenance_and_checks() -> None:
    manifest = _manifest()
    manifest.status = RunStatus.SUCCESS
    manifest.reproducible = True
    manifest.metrics = {
        "sequence": {"pdb_entity_identity": 1.0},
        "structure": {"rmsd_angstrom": 0.33, "lddt_ca": 0.99},
        "confidence": {"mean_plddt": 98.0},
        "limitations": ["pLDDT is not ground truth"],
        "not_evaluated": [{"metric": "TM-score", "reason": "not implemented"}],
    }
    manifest.checks = [
        CheckResult(
            check_id="uniprot_sequence_valid",
            title="UniProt sequence is valid",
            status=CheckStatus.PASSED,
            detail="142 residues",
        )
    ]
    report = build_report(manifest)
    assert report["status"] == "SUCCESS"
    assert report["error"] is None
    assert report["checks"][0]["check_id"] == "uniprot_sequence_valid"
    assert report["not_evaluated"][0]["metric"] == "TM-score"
    assert report["reproduction_command"].startswith("nexusforge protein-smoke --uniprot P69905")


def test_markdown_marks_scope_of_lddt() -> None:
    manifest = _manifest()
    manifest.status = RunStatus.SUCCESS
    manifest.reproducible = True
    manifest.metrics = {
        "structure": {
            "mapped_residue_count": 141,
            "rmsd_angstrom": 0.335,
            "lddt_ca": 0.9921,
            "lddt_ca_pairs": 3078,
            "lddt_ca_residues_scored": 141,
            "lddt_ca_residues_unscored": 0,
            "lddt_ca_inclusion_radius": 15.0,
            "lddt_ca_per_threshold": {"0.5": 0.97, "1": 1.0},
        },
        "confidence": {
            "residue_count": 142,
            "mean_plddt": 98.06,
            "median_plddt": 98.56,
            "min_plddt": 65.38,
            "max_plddt": 98.88,
            "bins": [],
        },
    }
    markdown = render(manifest)
    assert "lDDT-Cα" in markdown
    assert "not a full-atom lDDT" in markdown
    assert "pLDDT is a prediction-confidence estimate, not ground truth" in markdown
