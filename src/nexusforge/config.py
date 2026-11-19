"""Runtime configuration.

Settings are explicit, inspectable and serialised into every run manifest so a
run can be reproduced later. Environment variables use the ``NEXUSFORGE_``
prefix (see :func:`load_settings`).
"""

from __future__ import annotations

import os
from typing import Any

from pydantic import BaseModel, ConfigDict, Field, field_validator

from nexusforge import __version__

DEFAULT_USER_AGENT = (
    f"NexusForge/{__version__} (reproducible protein structure validation; "
    "+https://github.com/NexusForge/NexusForge)"
)


class Settings(BaseModel):
    """Immutable-ish configuration object for one run."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    # --- HTTP / provenance ---
    user_agent: str = DEFAULT_USER_AGENT
    contact: str | None = None
    timeout_seconds: float = Field(default=30.0, gt=0)
    max_attempts: int = Field(default=4, ge=1, le=10)
    retry_base_delay: float = Field(default=0.5, ge=0.0)
    retry_max_delay: float = Field(default=8.0, ge=0.0)

    # --- artefacts ---
    keep_raw_downloads: bool = True

    # --- scientific thresholds ---
    min_sequence_identity: float = Field(default=0.90, ge=0.0, le=1.0)
    min_sequence_coverage: float = Field(default=0.50, ge=0.0, le=1.0)
    lddt_inclusion_radius: float = Field(default=15.0, gt=0)
    lddt_thresholds: tuple[float, ...] = (0.5, 1.0, 2.0, 4.0)

    # --- optional resources ---
    fetch_pae: bool = True

    @field_validator("lddt_thresholds")
    @classmethod
    def _check_thresholds(cls, value: tuple[float, ...]) -> tuple[float, ...]:
        if not value:
            raise ValueError("lddt_thresholds must contain at least one distance threshold")
        if any(t <= 0 for t in value):
            raise ValueError("lddt_thresholds must be positive distances in angstroms")
        return tuple(sorted(value))

    @property
    def effective_user_agent(self) -> str:
        """User agent string, optionally carrying a contact address.

        AlphaFold DB serves its file endpoint only to non-generic clients, and
        every source in this project asks for a crawler-identifying agent.
        """
        if self.contact:
            return f"{self.user_agent} contact:{self.contact}"
        return self.user_agent

    def public_dict(self) -> dict[str, Any]:
        """Serialisable settings for the run manifest (no secrets are stored)."""
        data = self.model_dump(mode="json")
        data["lddt_thresholds"] = list(self.lddt_thresholds)
        return data


_ENV_FIELDS: dict[str, tuple[str, type]] = {
    "NEXUSFORGE_USER_AGENT": ("user_agent", str),
    "NEXUSFORGE_CONTACT": ("contact", str),
    "NEXUSFORGE_TIMEOUT_SECONDS": ("timeout_seconds", float),
    "NEXUSFORGE_MAX_ATTEMPTS": ("max_attempts", int),
    "NEXUSFORGE_RETRY_BASE_DELAY": ("retry_base_delay", float),
    "NEXUSFORGE_RETRY_MAX_DELAY": ("retry_max_delay", float),
    "NEXUSFORGE_KEEP_RAW_DOWNLOADS": ("keep_raw_downloads", bool),
    "NEXUSFORGE_MIN_SEQUENCE_IDENTITY": ("min_sequence_identity", float),
    "NEXUSFORGE_MIN_SEQUENCE_COVERAGE": ("min_sequence_coverage", float),
    "NEXUSFORGE_FETCH_PAE": ("fetch_pae", bool),
}


def _coerce(raw: str, kind: type) -> Any:
    if kind is bool:
        lowered = raw.strip().lower()
        if lowered in {"1", "true", "yes", "on"}:
            return True
        if lowered in {"0", "false", "no", "off"}:
            return False
        raise ValueError(f"cannot interpret {raw!r} as a boolean")
    return kind(raw)


def load_settings(env: dict[str, str] | None = None, **overrides: Any) -> Settings:
    """Build :class:`Settings` from environment variables plus explicit overrides.

    Explicit keyword overrides win over environment variables, which win over
    defaults. Unknown ``NEXUSFORGE_*`` variables are ignored so the config stays
    forward compatible.
    """
    source = os.environ if env is None else env
    values: dict[str, Any] = {}
    for env_name, (field, kind) in _ENV_FIELDS.items():
        if env_name in source and source[env_name] != "":
            values[field] = _coerce(source[env_name], kind)
    values.update({k: v for k, v in overrides.items() if v is not None})
    return Settings(**values)


__all__ = ["DEFAULT_USER_AGENT", "Settings", "load_settings"]
