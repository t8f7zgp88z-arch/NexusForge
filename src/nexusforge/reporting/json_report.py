"""Machine-readable JSON report.

The report is deliberately separate from the manifest: the manifest is the audit
trail (provenance, hashes, environment), the report is the scientific result
what was measured, what was not, and why.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from nexusforge.errors import ErrorCode, ReportingError
from nexusforge.models import ArtifactRecord, RunManifest, RunStatus
from nexusforge.reporting.manifest import (
    manifest_sha256,
    reproduction_command,
    serialize,
    write_bytes,
)

REPORT_SCHEMA_VERSION = "1.0"


def build_report(manifest: RunManifest) -> dict[str, Any]:
    """Assemble the report payload from a fully populated manifest."""
    if manifest.status is RunStatus.FAILED and manifest.error is None:
        raise ReportingError(
            ErrorCode.REPORT_WRITE_FAILED, "cannot build a report for an anonymised failure"
        )
    report: dict[str, Any] = {
        "schema_version": REPORT_SCHEMA_VERSION,
        "run_id": manifest.run_id,
        "generated_at": manifest.generated_at.isoformat(),
        "status": manifest.status.value,
        "request": manifest.request.model_dump(mode="json"),
        "software": manifest.software.model_dump(mode="json"),
        "reproducible": manifest.reproducible,
        "reproduction_command": reproduction_command(manifest),
        "manifest_sha256": manifest_sha256(manifest),
        "data": manifest.data,
        "sources": [
            {
                "source": provenance.source,
                "url": provenance.url,
                "fetched_at": provenance.fetched_at.isoformat(),
                "http_status": provenance.http_status,
                "sha256": provenance.sha256,
                "size_bytes": provenance.size_bytes,
            }
            for provenance in manifest.sources
        ],
        "alignment": (
            manifest.alignment.model_dump(mode="json") if manifest.alignment is not None else None
        ),
        "metrics": manifest.metrics,
        "checks": [check.model_dump(mode="json") for check in manifest.checks],
        "warnings": manifest.warnings,
        "not_evaluated": manifest.metrics.get("not_evaluated", []),
        "limitations": manifest.metrics.get("limitations", []),
        "error": manifest.error,
    }
    if manifest.status is RunStatus.FAILED:
        report["conclusion"] = (
            "Run failed. No scientific conclusion is drawn and no metric from this run may be "
            "cited as a result."
        )
    else:
        report["conclusion"] = (
            "Run completed. Metrics below describe the agreement between an AlphaFold DB "
            "prediction and one experimental structure, restricted to common Cα positions."
        )
    return report


def write_report(run_dir: Path, manifest: RunManifest) -> tuple[ArtifactRecord, dict[str, Any]]:
    report = build_report(manifest)
    artifact = write_bytes(run_dir / "report.json", serialize(report))
    return artifact, report


def write_metrics(run_dir: Path, manifest: RunManifest) -> ArtifactRecord:
    payload = {
        "run_id": manifest.run_id,
        "status": manifest.status.value,
        "metrics": manifest.metrics,
        "checks": [check.model_dump(mode="json") for check in manifest.checks],
    }
    return write_bytes(run_dir / "metrics.json", serialize(payload))


__all__ = ["REPORT_SCHEMA_VERSION", "build_report", "write_metrics", "write_report"]
