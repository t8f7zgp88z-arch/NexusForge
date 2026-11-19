"""Run manifest construction and writing.

The manifest is the audit trail: what was requested, what was fetched (with
SHA-256), which software versions ran, which checks executed, and whether the
run succeeded. :func:`ensure_consistent` makes the "failed runs are never
reported as SUCCESS" rule executable rather than aspirational.
"""

from __future__ import annotations

import hashlib
import json
import platform
import sys
from datetime import UTC, datetime
from importlib import metadata
from pathlib import Path
from typing import Any

from nexusforge import __version__
from nexusforge.config import Settings
from nexusforge.errors import ErrorCode, NexusForgeError, ReportingError
from nexusforge.models import (
    ArtifactRecord,
    RequestSpec,
    RunManifest,
    RunStatus,
    SoftwareInfo,
)

MANIFEST_SCHEMA_VERSION = "1.0"
_TRACKED_DEPENDENCIES = ("biopython", "httpx", "numpy", "pydantic", "typer")


def utcnow() -> datetime:
    return datetime.now(UTC)


def dependency_versions() -> dict[str, str]:
    versions: dict[str, str] = {}
    for name in _TRACKED_DEPENDENCIES:
        try:
            versions[name] = metadata.version(name)
        except metadata.PackageNotFoundError:  # pragma: no cover - environment dependent
            versions[name] = "not-installed"
    return versions


def build_software_info() -> SoftwareInfo:
    return SoftwareInfo(
        name="nexusforge",
        version=__version__,
        python_version=sys.version.split()[0],
        platform=f"{platform.system()} {platform.release()} ({platform.machine()})",
        dependencies=dependency_versions(),
    )


def build_manifest(
    *,
    request: RequestSpec,
    settings: Settings,
    run_id: str,
    generated_at: datetime | None = None,
) -> RunManifest:
    return RunManifest(
        schema_version=MANIFEST_SCHEMA_VERSION,
        run_id=run_id,
        generated_at=generated_at or utcnow(),
        status=RunStatus.FAILED,
        request=request,
        software=build_software_info(),
        settings=settings.public_dict(),
        reproducible=False,
    )


def ensure_consistent(manifest: RunManifest) -> None:
    """Reject internally contradictory manifests before they reach disk."""
    if manifest.status is RunStatus.SUCCESS:
        if manifest.error is not None:
            raise ReportingError(
                ErrorCode.REPORT_WRITE_FAILED,
                "manifest claims SUCCESS but carries an error",
            )
        missing = [key for key in ("sequence", "structure") if key not in manifest.metrics]
        if missing:
            raise ReportingError(
                ErrorCode.REPORT_WRITE_FAILED,
                f"manifest claims SUCCESS but is missing metrics: {', '.join(missing)}",
            )
        if not manifest.reproducible:
            raise ReportingError(
                ErrorCode.REPORT_WRITE_FAILED,
                "manifest claims SUCCESS but is not marked reproducible",
            )
    elif manifest.error is None:
        raise ReportingError(
            ErrorCode.REPORT_WRITE_FAILED,
            "manifest is FAILED but carries no structured error",
        )


def serialize(model: Any) -> bytes:
    """Canonical JSON bytes for a pydantic model or plain payload."""
    payload = model.model_dump(mode="json") if hasattr(model, "model_dump") else model
    text = json.dumps(payload, indent=2, ensure_ascii=False, sort_keys=False)
    return (text + "\n").encode("utf-8")


def write_bytes(path: Path, content: bytes) -> ArtifactRecord:
    """Write ``content`` atomically and return a hashed artifact record."""
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_name(path.name + ".partial")
        temporary.write_bytes(content)
        temporary.replace(path)
    except OSError as exc:
        raise ReportingError(
            ErrorCode.REPORT_WRITE_FAILED,
            f"could not write {path}: {exc}",
            details={"path": str(path)},
        ) from exc
    return ArtifactRecord(
        kind=path.name,
        path=path.name,
        sha256=hashlib.sha256(content).hexdigest(),
        size_bytes=len(content),
    )


def artifact_for_file(path: Path, run_dir: Path, *, kind: str) -> ArtifactRecord:
    """Describe an already-written file as a hashed artifact."""
    content = path.read_bytes()
    return ArtifactRecord(
        kind=kind,
        path=str(path.relative_to(run_dir)),
        sha256=hashlib.sha256(content).hexdigest(),
        size_bytes=len(content),
    )


def write_manifest(run_dir: Path, manifest: RunManifest) -> ArtifactRecord:
    ensure_consistent(manifest)
    return write_bytes(run_dir / "manifest.json", serialize(manifest))


def manifest_sha256(manifest: RunManifest) -> str:
    return hashlib.sha256(serialize(manifest)).hexdigest()


def reproduction_command(manifest: RunManifest) -> str:
    return (
        "nexusforge protein-smoke "
        f"--uniprot {manifest.request.uniprot_accession} "
        f"--pdb {manifest.request.pdb_id} "
        f"--chain {manifest.request.chain_id} "
        f"--output runs/{manifest.run_id}"
    )


def fail_manifest(manifest: RunManifest, error: NexusForgeError) -> RunManifest:
    """Record a failure without erasing what was already retrieved."""
    manifest.status = RunStatus.FAILED
    manifest.error = error.to_dict()
    manifest.reproducible = False
    return manifest


__all__ = [
    "MANIFEST_SCHEMA_VERSION",
    "artifact_for_file",
    "build_manifest",
    "build_software_info",
    "dependency_versions",
    "ensure_consistent",
    "fail_manifest",
    "manifest_sha256",
    "reproduction_command",
    "serialize",
    "utcnow",
    "write_bytes",
    "write_manifest",
]
