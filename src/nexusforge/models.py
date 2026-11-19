"""Serialisable domain contract shared by clients, validation and reporting.

Everything that leaves the pipeline (reports, manifests, metrics) is built from
these models, so the report schema is explicit and versioned rather than a pile
of ad-hoc dictionaries.
"""

from __future__ import annotations

from datetime import datetime
from enum import StrEnum
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator


class RunStatus(StrEnum):
    """Overall outcome of a run. A failed run is never reported as SUCCESS."""

    SUCCESS = "SUCCESS"
    FAILED = "FAILED"


class CheckStatus(StrEnum):
    PASSED = "passed"
    FAILED = "failed"
    SKIPPED = "skipped"


class Base(BaseModel):
    model_config = ConfigDict(extra="forbid", validate_assignment=True)


class Provenance(Base):
    """Where bytes came from, when, and how to verify them again."""

    source: str = Field(description="Logical data source, e.g. 'uniprot'")
    url: str
    fetched_at: datetime = Field(description="UTC timestamp of the successful fetch")
    http_status: int
    sha256: str
    size_bytes: int = Field(ge=0)
    content_type: str | None = None
    attempts: int = Field(default=1, ge=1)

    @field_validator("sha256")
    @classmethod
    def _sha256_shape(cls, value: str) -> str:
        if len(value) != 64 or any(c not in "0123456789abcdef" for c in value):
            raise ValueError("sha256 must be 64 lowercase hex characters")
        return value


class ArtifactRecord(Base):
    """A file written into the run directory (raw download or report)."""

    kind: str
    path: str = Field(description="Path relative to the run directory")
    url: str | None = None
    sha256: str
    size_bytes: int = Field(ge=0)


class SoftwareInfo(Base):
    name: str = "nexusforge"
    version: str
    python_version: str
    platform: str
    dependencies: dict[str, str] = Field(default_factory=dict)


class RequestSpec(Base):
    """The exact request a run answers."""

    uniprot_accession: str
    pdb_id: str
    chain_id: str
    comparison: str = "alphafold_prediction_vs_experimental_structure"


# --------------------------------------------------------------------------- #
# UniProt
# --------------------------------------------------------------------------- #
class UniProtRecord(Base):
    accession: str
    entry_name: str
    reviewed: bool
    protein_name: str | None = None
    alternative_names: list[str] = Field(default_factory=list)
    gene_names: list[str] = Field(default_factory=list)
    organism: str
    taxon_id: int | None = None
    sequence: str
    sequence_length: int = Field(gt=0)
    sequence_md5: str | None = None
    sequence_crc64: str | None = None

    @field_validator("sequence")
    @classmethod
    def _upper(cls, value: str) -> str:
        return value.strip().upper()


# --------------------------------------------------------------------------- #
# AlphaFold DB
# --------------------------------------------------------------------------- #
class AlphaFoldFiles(Base):
    """Model file URLs exactly as advertised by the prediction record."""

    cif_url: str | None = None
    pdb_url: str | None = None
    plddt_doc_url: str | None = None
    pae_doc_url: str | None = None
    pae_image_url: str | None = None
    msa_url: str | None = None


class AlphaFoldPrediction(Base):
    entry_id: str
    model_entity_id: str | None = None
    uniprot_accession: str
    tool_used: str | None = None
    latest_version: int | None = None
    model_created_date: datetime | None = None
    global_plddt: float | None = Field(default=None, ge=0, le=100)
    fraction_plddt_very_low: float | None = None
    fraction_plddt_low: float | None = None
    fraction_plddt_confident: float | None = None
    fraction_plddt_very_high: float | None = None
    sequence: str
    sequence_start: int | None = None
    sequence_end: int | None = None
    chain_id: str | None = None
    is_complex: bool | None = None
    files: AlphaFoldFiles
    schema_notes: list[str] = Field(
        default_factory=list,
        description="Field-name fallbacks or missing optional fields observed while parsing",
    )


class ConfidenceScores(Base):
    """Per-residue pLDDT (a prediction-confidence value, never ground truth)."""

    source: Literal["alphafold_confidence_json", "alphafold_cif_bfactor"]
    residue_numbers: list[int]
    scores: list[float]

    @field_validator("scores")
    @classmethod
    def _range(cls, value: list[float]) -> list[float]:
        if any(not (0.0 <= s <= 100.0) for s in value):
            raise ValueError("pLDDT values must lie in [0, 100]")
        return value

    def as_mapping(self) -> dict[int, float]:
        return dict(zip(self.residue_numbers, self.scores, strict=True))


class PaeSummary(Base):
    """Aggregate statistics only: PAE matrices are large and stay in raw files."""

    shape: tuple[int, int]
    max_predicted_aligned_error: float | None = None
    mean_predicted_aligned_error: float
    p90_predicted_aligned_error: float | None = None


# --------------------------------------------------------------------------- #
# RCSB PDB
# --------------------------------------------------------------------------- #
class RcsbEntry(Base):
    pdb_id: str
    title: str
    experimental_methods: list[str] = Field(default_factory=list)
    resolution_angstrom: float | None = Field(default=None, gt=0)
    deposited_at: datetime | None = None
    released_at: datetime | None = None
    last_revision_at: datetime | None = None
    polymer_entity_ids: list[str] = Field(default_factory=list)
    deposited_atom_count: int | None = None
    authors: list[str] = Field(default_factory=list)


class RcsbPolymerEntity(Base):
    """A polymer entity as described by the RCSB entity service."""

    pdb_id: str
    entity_id: str
    asym_ids: list[str] = Field(default_factory=list)
    auth_asym_ids: list[str] = Field(default_factory=list)
    uniprot_accessions: list[str] = Field(default_factory=list)
    entity_name: str | None = None
    sequence: str
    coverage_of_reference: float | None = Field(default=None, ge=0, le=1)


class RcsbChainMapping(Base):
    """How an author chain in a PDB entry relates to a UniProt accession."""

    pdb_id: str
    entity_id: str
    auth_chain_id: str
    asym_ids: list[str] = Field(default_factory=list)
    uniprot_accessions: list[str] = Field(default_factory=list)
    entity_name: str | None = None
    sequence: str
    coverage_of_reference: float | None = Field(default=None, ge=0, le=1)

    @field_validator("sequence")
    @classmethod
    def _upper(cls, value: str) -> str:
        return value.strip().upper()


class CaChainSummary(Base):
    """Summary of extracted Cα traces, enough to audit without shipping coordinates."""

    label: str
    chain_id: str | None
    residue_count: int = Field(ge=0)
    first_residue_number: int | None = None
    last_residue_number: int | None = None
    non_standard_residues: list[str] = Field(default_factory=list)
    duplicate_residue_numbers: list[int] = Field(default_factory=list)


# --------------------------------------------------------------------------- #
# Sequence validation
# --------------------------------------------------------------------------- #
class SequenceStats(Base):
    length: int = Field(gt=0)
    invalid_characters: list[str] = Field(default_factory=list)
    ambiguous_characters: list[str] = Field(default_factory=list)


class AlignmentResult(Base):
    query_id: str
    target_id: str
    query_length: int
    target_length: int
    aligned_length: int
    matches: int
    mismatches: int
    gap_openings: int
    gap_positions: int
    identity: float = Field(ge=0, le=1)
    coverage_query: float = Field(ge=0, le=1)
    coverage_target: float = Field(ge=0, le=1)
    n_optimal_alignments: int = Field(default=1, ge=1)
    mapping: dict[int, int] = Field(
        default_factory=dict,
        description="0-based query index -> 0-based target index for aligned residues",
    )
    aligned_query_preview: str | None = None
    aligned_target_preview: str | None = None


# --------------------------------------------------------------------------- #
# Structure evaluation
# --------------------------------------------------------------------------- #
class ResidueScore(Base):
    uniprot_position: int
    experimental_residue_number: int
    predicted_sequence_position: int
    residue_name: str
    plddt: float | None = None
    lddt_ca: float | None = None


class StructureEvaluation(Base):
    mapped_residue_count: int
    experimental_residue_count: int
    predicted_residue_count: int
    experimental_only_residues: list[int] = Field(default_factory=list)
    predicted_only_residues: list[int] = Field(default_factory=list)
    rmsd_angstrom: float
    lddt_ca: float = Field(ge=0, le=1)
    lddt_ca_per_threshold: dict[str, float] = Field(default_factory=dict)
    lddt_ca_pairs: int
    lddt_ca_residues_scored: int
    lddt_ca_residues_unscored: int
    lddt_ca_inclusion_radius: float
    lddt_ca_thresholds: list[float]
    per_residue: list[ResidueScore] = Field(default_factory=list)
    notes: list[str] = Field(default_factory=list)


# --------------------------------------------------------------------------- #
# Confidence evaluation
# --------------------------------------------------------------------------- #
class PlddtBin(Base):
    label: str
    lower: float
    upper: float
    residue_count: int
    fraction: float = Field(ge=0, le=1)
    mean_lddt_ca: float | None = None


class ConfidenceEvaluation(Base):
    source: str
    residue_count: int
    min_plddt: float
    max_plddt: float
    mean_plddt: float
    median_plddt: float
    bins: list[PlddtBin] = Field(default_factory=list)
    cross_check_max_abs_delta: float | None = Field(
        default=None,
        description="Largest |Δ| between confidence-JSON pLDDT and CIF B-factor pLDDT",
    )
    pearson_plddt_lddt: float | None = None
    spearman_plddt_lddt: float | None = None
    correlated_residue_count: int = 0
    interpretation: str = (
        "pLDDT is an AlphaFold per-residue confidence estimate. It is compared with, "
        "not equated to, the experimentally referenced lDDT-Cα."
    )


# --------------------------------------------------------------------------- #
# Checks / run bookkeeping
# --------------------------------------------------------------------------- #
class CheckResult(Base):
    check_id: str
    title: str
    status: CheckStatus
    detail: str | None = None
    metrics: dict[str, Any] = Field(default_factory=dict)


class AlignmentSummary(Base):
    uniprot_vs_alphafold: AlignmentResult
    uniprot_vs_pdb_entity: AlignmentResult


class RunManifest(Base):
    schema_version: str = "1.0"
    run_id: str
    generated_at: datetime
    status: RunStatus
    request: RequestSpec
    software: SoftwareInfo
    settings: dict[str, Any] = Field(default_factory=dict)
    sources: list[Provenance] = Field(default_factory=list)
    artifacts: list[ArtifactRecord] = Field(default_factory=list)
    data: dict[str, Any] = Field(
        default_factory=dict, description="Resolved identifiers, versions and model metadata"
    )
    alignment: AlignmentSummary | None = None
    checks: list[CheckResult] = Field(default_factory=list)
    metrics: dict[str, Any] = Field(default_factory=dict)
    warnings: list[str] = Field(default_factory=list)
    error: dict[str, Any] | None = None
    reproducible: bool = False
    reproduction_command: str | None = None


__all__ = [
    "AlignmentResult",
    "AlignmentSummary",
    "AlphaFoldFiles",
    "AlphaFoldPrediction",
    "ArtifactRecord",
    "CaChainSummary",
    "CheckResult",
    "CheckStatus",
    "ConfidenceEvaluation",
    "ConfidenceScores",
    "PaeSummary",
    "PlddtBin",
    "Provenance",
    "RcsbChainMapping",
    "RcsbEntry",
    "RcsbPolymerEntity",
    "RequestSpec",
    "ResidueScore",
    "RunManifest",
    "RunStatus",
    "SequenceStats",
    "SoftwareInfo",
    "StructureEvaluation",
    "UniProtRecord",
]
