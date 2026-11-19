"""End-to-end protein structure smoke pipeline.

Stages, each with its own failure domain:

1. ``data_retrieval``      — UniProt, AlphaFold DB, RCSB (metadata + files)
2. ``sequence_validation`` — is this chain really this UniProt protein?
3. ``structure_parsing``   — Cα traces, coordinate sanity
4. ``structure_evaluation``— superposition, RMSD, lDDT-Cα
5. ``confidence_evaluation``— pLDDT validation, binning, correlation with lDDT-Cα
6. ``reporting``           — manifest, JSON/Markdown reports, logs

A failure in any stage stops the run, is recorded with its stage, and produces a
FAILED manifest. Artifacts retrieved before the failure are kept so the run can
be diagnosed.
"""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from urllib.parse import urlsplit

from nexusforge.clients.alphafold import AlphaFoldClient
from nexusforge.clients.rcsb import RcsbClient, validate_chain_id, validate_pdb_id
from nexusforge.clients.uniprot import UniProtClient, validate_accession
from nexusforge.config import Settings, load_settings
from nexusforge.errors import (
    ConfidenceError,
    ErrorCode,
    InputError,
    NexusForgeError,
    ReportingError,
    StructureError,
    ValidationError,
)
from nexusforge.http import HttpClient
from nexusforge.models import (
    AlignmentSummary,
    ArtifactRecord,
    CheckResult,
    CheckStatus,
    PaeSummary,
    RequestSpec,
    ResidueScore,
    RunManifest,
    RunStatus,
    StructureEvaluation,
)
from nexusforge.reporting import markdown_report
from nexusforge.reporting.json_report import write_metrics, write_report
from nexusforge.reporting.manifest import (
    artifact_for_file,
    build_manifest,
    fail_manifest,
    reproduction_command,
    write_bytes,
    write_manifest,
)
from nexusforge.validation.confidence import (
    cross_check_confidence,
    plddt_from_trace,
    summarise_confidence,
    validate_confidence,
)
from nexusforge.validation.sequence import align_sequences, describe_gaps, validate_sequence
from nexusforge.validation.structure import (
    CaTrace,
    evaluate_structure,
    extract_ca_trace,
    match_residues,
)

logger = logging.getLogger("nexusforge.pipeline")

LIMITATIONS = [
    "AlphaFold DB predictions are consumed as published; this project does not run the "
    "AlphaFold model itself.",
    "The comparison uses one experimental chain. Quaternary context (e.g. haemoglobin "
    "tetramer assembly) is not modelled or scored.",
    "Only positions resolved in both structures are compared; unobserved residues are "
    "reported, never imputed.",
    "lDDT-Cα is restricted to Cα atoms and is not comparable with a full-atom lDDT value.",
    "pLDDT is a model confidence estimate; the lDDT-Cα column is measured against one "
    "experimental structure, so agreement or disagreement between them is an observation, "
    "not a calibrated error.",
]


@dataclass(frozen=True)
class SmokeRequest:
    uniprot_accession: str
    pdb_id: str
    chain_id: str
    output_dir: Path
    overwrite: bool = False


@dataclass
class RunOutcome:
    run_id: str
    run_dir: Path
    manifest: RunManifest
    report: dict[str, Any] | None
    error: NexusForgeError | None

    @property
    def succeeded(self) -> bool:
        return self.error is None and self.manifest.status is RunStatus.SUCCESS

    @property
    def exit_code(self) -> int:
        return self.error.exit_code if self.error is not None else 0


def _invert(mapping: dict[int, int]) -> dict[int, int]:
    return {target: source for source, target in mapping.items()}


def _compose(first: dict[int, int], second: dict[int, int]) -> dict[int, int]:
    return {source: second[mid] for source, mid in first.items() if mid in second}


def _url_filename(url: str, fallback: str) -> str:
    name = Path(urlsplit(url).path).name
    return name or fallback


def _finite(value: float) -> bool:
    return value == value and value not in (float("inf"), float("-inf"))


class ProteinSmokePipeline:
    """Retrieve → validate → evaluate → report, with full provenance."""

    def __init__(
        self,
        settings: Settings | None = None,
        *,
        http: HttpClient | None = None,
    ) -> None:
        self._settings = settings or load_settings()
        self._owns_http = http is None
        self._http = http or HttpClient(self._settings)
        self._uniprot = UniProtClient(self._http)
        self._alphafold = AlphaFoldClient(self._http)
        self._rcsb = RcsbClient(self._http)
        self._archived_urls: set[str] = set()
        self._artifacts: list[ArtifactRecord] = []
        self._log_handler: logging.Handler | None = None

    # -- lifecycle ---------------------------------------------------------
    def close(self) -> None:
        self._detach_logging()
        if self._owns_http:
            self._http.close()

    def __enter__(self) -> ProteinSmokePipeline:
        return self

    def __exit__(self, *exc_info: object) -> None:
        self.close()

    # -- public API --------------------------------------------------------
    def run(self, request: SmokeRequest) -> RunOutcome:
        accession = validate_accession(request.uniprot_accession)
        pdb_id = validate_pdb_id(request.pdb_id)
        chain_id = validate_chain_id(request.chain_id)
        run_dir = Path(request.output_dir)
        run_id = f"{accession}-{pdb_id}-{chain_id}"

        self._guard_output_dir(run_dir, request.overwrite)
        run_dir.mkdir(parents=True, exist_ok=True)
        (run_dir / "raw").mkdir(exist_ok=True)
        (run_dir / "logs").mkdir(exist_ok=True)

        manifest = build_manifest(
            request=RequestSpec(
                uniprot_accession=accession,
                pdb_id=pdb_id,
                chain_id=chain_id,
            ),
            settings=self._settings,
            run_id=run_id,
        )
        self._attach_logging(run_dir)
        started = time.monotonic()
        error: NexusForgeError | None = None
        try:
            self._execute(manifest, run_dir)
            manifest.status = RunStatus.SUCCESS
            manifest.reproducible = True
            manifest.reproduction_command = reproduction_command(manifest)
        except NexusForgeError as exc:
            error = exc
            logger.error("run failed: %s", exc)
            fail_manifest(manifest, exc)
        except Exception as exc:
            import traceback

            logger.exception("unexpected failure")
            error = NexusForgeError(
                ErrorCode.INTERNAL_ERROR,
                f"unexpected {type(exc).__name__}: {exc}",
                details={"traceback": traceback.format_exc()},
            )
            fail_manifest(manifest, error)
        finally:
            manifest.metrics["timings_seconds"] = {"total": round(time.monotonic() - started, 3)}
            self._detach_logging()
            report = self._finalize(manifest, run_dir)
        return RunOutcome(
            run_id=run_id, run_dir=run_dir, manifest=manifest, report=report, error=error
        )

    # -- setup helpers -----------------------------------------------------
    @staticmethod
    def _guard_output_dir(run_dir: Path, overwrite: bool) -> None:
        if run_dir.exists() and any(run_dir.iterdir()) and not overwrite:
            raise InputError(
                ErrorCode.INVALID_INPUT,
                f"output directory {run_dir} already contains files; pass --force to overwrite",
                details={"output_dir": str(run_dir)},
            )

    def _attach_logging(self, run_dir: Path) -> None:
        log_path = run_dir / "logs" / "run.log"
        handler = logging.FileHandler(log_path, mode="w", encoding="utf-8")
        handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)-8s %(name)s: %(message)s"))
        handler.setLevel(logging.DEBUG)
        logging.getLogger("nexusforge").addHandler(handler)
        logging.getLogger("nexusforge").setLevel(logging.DEBUG)
        self._log_handler = handler

    def _detach_logging(self) -> None:
        if self._log_handler is not None:
            root = logging.getLogger("nexusforge")
            root.removeHandler(self._log_handler)
            self._log_handler.close()
            self._log_handler = None

    # -- main flow ---------------------------------------------------------
    def _execute(self, manifest: RunManifest, run_dir: Path) -> None:
        request = manifest.request
        accession = request.uniprot_accession
        pdb_id = request.pdb_id
        chain_id = request.chain_id
        settings = self._settings
        checks = manifest.checks
        warnings = manifest.warnings

        # --- 1. UniProt ----------------------------------------------------
        uniprot, uniprot_provenance = self._uniprot.fetch(accession)
        manifest.sources.append(uniprot_provenance)
        sequence_stats = validate_sequence(uniprot.sequence, label="uniprot sequence")
        manifest.data["uniprot"] = {
            "accession": uniprot.accession,
            "entry_name": uniprot.entry_name,
            "protein_name": uniprot.protein_name,
            "alternative_names": uniprot.alternative_names,
            "gene_names": uniprot.gene_names,
            "organism": uniprot.organism,
            "taxon_id": uniprot.taxon_id,
            "reviewed": uniprot.reviewed,
            "sequence_length": uniprot.sequence_length,
            "sequence_md5": uniprot.sequence_md5,
            "ambiguous_residues": sequence_stats.ambiguous_characters,
        }
        checks.append(
            CheckResult(
                check_id="uniprot_sequence_valid",
                title="UniProt sequence uses valid amino-acid characters",
                status=CheckStatus.PASSED,
                detail=f"{sequence_stats.length} residues, "
                f"{len(sequence_stats.ambiguous_characters)} ambiguous code(s)",
                metrics={"length": sequence_stats.length},
            )
        )
        self._archive_responses(run_dir)

        # --- 2. AlphaFold metadata ----------------------------------------
        prediction, alphafold_provenance = self._alphafold.fetch_prediction(accession)
        manifest.sources.append(alphafold_provenance)
        if prediction.uniprot_accession.upper() != accession:
            raise ValidationError(
                ErrorCode.SEQUENCE_MISMATCH,
                f"AlphaFold entry {prediction.entry_id} belongs to "
                f"{prediction.uniprot_accession}, not {accession}",
            )
        manifest.data["alphafold"] = {
            "entry_id": prediction.entry_id,
            "model_entity_id": prediction.model_entity_id,
            "tool_used": prediction.tool_used,
            "latest_version": prediction.latest_version,
            "model_created_date": (
                prediction.model_created_date.isoformat() if prediction.model_created_date else None
            ),
            "global_plddt": prediction.global_plddt,
            "sequence_length": len(prediction.sequence),
            "sequence_start": prediction.sequence_start,
            "sequence_end": prediction.sequence_end,
            "chain_id": prediction.chain_id,
            "is_complex": prediction.is_complex,
            "files": prediction.files.model_dump(mode="json"),
            "pae_available": bool(prediction.files.pae_doc_url),
        }
        for note in prediction.schema_notes:
            warnings.append(f"AlphaFold schema note: {note}")
        self._archive_responses(run_dir)

        # --- 3. UniProt vs AlphaFold sequence ------------------------------
        alphafold_alignment = align_sequences(
            uniprot.sequence,
            prediction.sequence,
            query_id=f"uniprot:{accession}",
            target_id=f"alphafold:{prediction.entry_id}",
            min_identity=settings.min_sequence_identity,
            min_coverage=settings.min_sequence_coverage,
        )
        checks.append(
            CheckResult(
                check_id="alphafold_sequence_matches_uniprot",
                title="AlphaFold model sequence corresponds to the UniProt entry",
                status=CheckStatus.PASSED,
                detail=f"identity {alphafold_alignment.identity:.2%}, coverage "
                f"{alphafold_alignment.coverage_query:.2%}",
                metrics={
                    "identity": alphafold_alignment.identity,
                    "coverage": alphafold_alignment.coverage_query,
                },
            )
        )

        if prediction.files.cif_url is None:
            raise ValidationError(
                ErrorCode.INVALID_API_RESPONSE,
                "AlphaFold prediction advertises no mmCIF structure file",
                details={"entry_id": prediction.entry_id},
            )

        # --- 4. AlphaFold model file + confidence --------------------------
        prediction_cif = (
            run_dir / "raw" / _url_filename(prediction.files.cif_url, f"{prediction.entry_id}.cif")
        )
        artifact, download_result = self._http.download(
            prediction.files.cif_url,
            prediction_cif,
            source="alphafold-file",
            kind="predicted_structure_cif",
            accept="chemical/x-mmcif,text/plain,*/*",
        )
        manifest.sources.append(download_result.provenance)
        self._record_artifact(manifest, artifact, run_dir, prediction_cif)
        self._archived_urls.add(prediction.files.cif_url)
        self._archive_responses(run_dir)

        if prediction.files.plddt_doc_url is None:
            raise ConfidenceError(
                ErrorCode.CONFIDENCE_LENGTH_MISMATCH,
                "AlphaFold prediction advertises no per-residue confidence document",
                details={"entry_id": prediction.entry_id},
            )
        confidence_scores, confidence_provenance = self._alphafold.fetch_confidence(
            prediction.files.plddt_doc_url
        )
        manifest.sources.append(confidence_provenance)
        self._archive_responses(run_dir)

        validated_scores = validate_confidence(
            confidence_scores,
            sequence_length=len(prediction.sequence),
            sequence_start=prediction.sequence_start,
            trace_label=f"alphafold:{prediction.entry_id}",
        )
        checks.append(
            CheckResult(
                check_id="plddt_count_matches_sequence",
                title="pLDDT vector length equals the modelled sequence length",
                status=CheckStatus.PASSED,
                detail=f"{len(validated_scores)} values for "
                f"{len(prediction.sequence)} modelled residues",
                metrics={"values": len(validated_scores)},
            )
        )

        pae_summary: PaeSummary | None = None
        if settings.fetch_pae and prediction.files.pae_doc_url:
            try:
                pae_summary, pae_provenance = self._alphafold.fetch_pae(
                    prediction.files.pae_doc_url
                )
                manifest.sources.append(pae_provenance)
                self._archive_responses(run_dir)
            except NexusForgeError as exc:
                warnings.append(
                    f"PAE data was advertised but could not be used ({exc.code.value}); "
                    "PAE statistics are omitted"
                )
                logger.warning("PAE fetch failed: %s", exc)

        # --- 5. RCSB metadata + chain resolution ---------------------------
        entry, entry_provenance = self._rcsb.fetch_entry(pdb_id)
        manifest.sources.append(entry_provenance)
        self._archive_responses(run_dir)
        mapping, mapping_provenance = self._rcsb.resolve_chain(
            pdb_id, chain_id, accession, entry=entry
        )
        manifest.sources.extend(mapping_provenance)
        self._archive_responses(run_dir)
        manifest.data["experimental"] = {
            "pdb_id": entry.pdb_id,
            "title": entry.title,
            "experimental_methods": entry.experimental_methods,
            "resolution_angstrom": entry.resolution_angstrom,
            "deposited_at": entry.deposited_at.isoformat() if entry.deposited_at else None,
            "released_at": entry.released_at.isoformat() if entry.released_at else None,
            "last_revision_at": (
                entry.last_revision_at.isoformat() if entry.last_revision_at else None
            ),
            "deposited_atom_count": entry.deposited_atom_count,
            "entity_id": mapping.entity_id,
            "entity_name": mapping.entity_name,
            "uniprot_accessions": mapping.uniprot_accessions,
            "coverage_of_reference": mapping.coverage_of_reference,
            "chain_id": mapping.auth_chain_id,
        }
        checks.append(
            CheckResult(
                check_id="pdb_chain_maps_to_uniprot",
                title="Experimental chain is cross-referenced to the UniProt accession",
                status=CheckStatus.PASSED,
                detail=f"entity {mapping.entity_id} of {pdb_id} lists "
                f"{', '.join(mapping.uniprot_accessions)}",
                metrics={
                    "entity_id": mapping.entity_id,
                    "uniprot_accessions": mapping.uniprot_accessions,
                },
            )
        )

        # --- 6. Experimental coordinates -----------------------------------
        experimental_cif = run_dir / "raw" / f"{pdb_id}.cif"
        artifact, download_result = self._rcsb.download_structure(pdb_id, experimental_cif)
        manifest.sources.append(download_result.provenance)
        self._record_artifact(manifest, artifact, run_dir, experimental_cif)
        self._archived_urls.add(self._rcsb.structure_url(pdb_id))

        # --- 7. UniProt vs PDB entity sequence -----------------------------
        entity_alignment = align_sequences(
            uniprot.sequence,
            mapping.sequence,
            query_id=f"uniprot:{accession}",
            target_id=f"pdb:{pdb_id}/entity:{mapping.entity_id}",
            min_identity=settings.min_sequence_identity,
            min_coverage=settings.min_sequence_coverage,
        )
        checks.append(
            CheckResult(
                check_id="pdb_entity_sequence_matches_uniprot",
                title="Experimental entity sequence corresponds to the UniProt entry",
                status=CheckStatus.PASSED,
                detail=f"identity {entity_alignment.identity:.2%}, coverage "
                f"{entity_alignment.coverage_query:.2%}, {describe_gaps(entity_alignment)}",
                metrics={
                    "identity": entity_alignment.identity,
                    "coverage_query": entity_alignment.coverage_query,
                    "gap_openings": entity_alignment.gap_openings,
                },
            )
        )
        manifest.alignment = AlignmentSummary(
            uniprot_vs_alphafold=alphafold_alignment,
            uniprot_vs_pdb_entity=entity_alignment,
        )

        # --- 8. Cα traces ---------------------------------------------------
        experimental_trace = extract_ca_trace(
            experimental_cif, label=f"pdb:{pdb_id}", chain_id=chain_id
        )
        predicted_trace = extract_ca_trace(
            prediction_cif, label=f"alphafold:{prediction.entry_id}", chain_id=prediction.chain_id
        )
        if predicted_trace.non_standard:
            warnings.append(
                f"AlphaFold model contains non-standard residues: "
                f"{', '.join(predicted_trace.non_standard)}"
            )
        if experimental_trace.non_standard:
            warnings.append(
                f"Experimental chain contains non-standard residues: "
                f"{', '.join(experimental_trace.non_standard)}"
            )
        if experimental_trace.missing_ca:
            warnings.append(
                f"{len(experimental_trace.missing_ca)} experimental residue(s) have no Cα atom "
                "and were excluded from the comparison"
            )
        if predicted_trace.missing_ca:
            warnings.append(
                f"{len(predicted_trace.missing_ca)} predicted residue(s) have no Cα atom"
            )
        checks.append(
            CheckResult(
                check_id="coordinates_parsed",
                title="Both structures parsed into Cα traces",
                status=CheckStatus.PASSED,
                detail=f"experimental {len(experimental_trace)} Cα, prediction "
                f"{len(predicted_trace)} Cα",
                metrics={
                    "experimental_residues": len(experimental_trace),
                    "predicted_residues": len(predicted_trace),
                },
            )
        )
        manifest.data["traces"] = {
            "experimental": experimental_trace.summary().model_dump(mode="json"),
            "prediction": predicted_trace.summary().model_dump(mode="json"),
        }

        # --- 9. Residue-level reference mapping ----------------------------
        experimental_to_entity = align_sequences(
            experimental_trace.sequence,
            mapping.sequence,
            query_id=f"pdb:{pdb_id}/chain:{chain_id} (observed Cα)",
            target_id=f"pdb:{pdb_id}/entity:{mapping.entity_id}",
            min_identity=0.98,
        )
        if experimental_to_entity.coverage_target < 0.5:
            raise ValidationError(
                ErrorCode.SEQUENCE_MISMATCH,
                f"observed chain {chain_id} covers only "
                f"{experimental_to_entity.coverage_target:.1%} of the entity sequence; "
                "the coordinate set does not represent the entity",
            )
        predicted_to_alphafold = align_sequences(
            predicted_trace.sequence,
            prediction.sequence,
            query_id=f"alphafold:{prediction.entry_id} (observed Cα)",
            target_id=f"alphafold:{prediction.entry_id} (sequence)",
            min_identity=0.98,
        )
        experimental_to_uniprot = _compose(
            experimental_to_entity.mapping, _invert(entity_alignment.mapping)
        )
        predicted_to_uniprot = _compose(
            predicted_to_alphafold.mapping, _invert(alphafold_alignment.mapping)
        )
        if not experimental_to_uniprot or not predicted_to_uniprot:
            raise StructureError(
                ErrorCode.NO_BASELINE,
                "residue mapping to the UniProt reference is empty",
                details={
                    "experimental_mapped": len(experimental_to_uniprot),
                    "predicted_mapped": len(predicted_to_uniprot),
                },
            )

        matching = match_residues(
            experimental_trace,
            predicted_trace,
            experimental_trace_to_reference=experimental_to_uniprot,
            predicted_trace_to_reference=predicted_to_uniprot,
        )

        # --- 10. Structure evaluation --------------------------------------
        structure_metrics = evaluate_structure(
            matching,
            experimental_trace,
            predicted_trace,
            inclusion_radius=settings.lddt_inclusion_radius,
            thresholds=settings.lddt_thresholds,
        )
        lddt = structure_metrics.lddt
        checks.append(
            CheckResult(
                check_id="lddt_ca_defined",
                title="lDDT-Cα has defined residue pairs and scored residues",
                status=CheckStatus.PASSED,
                detail=f"{lddt.pairs} pairs within {lddt.radius} Å, "
                f"{lddt.scored_residues} residues scored",
                metrics={"pairs": lddt.pairs, "scored_residues": lddt.scored_residues},
            )
        )
        checks.append(
            CheckResult(
                check_id="rmsd_computed",
                title="Cα RMSD computed after optimal superposition",
                status=CheckStatus.PASSED,
                detail=f"{structure_metrics.rmsd_angstrom:.3f} Å over "
                f"{len(matching.uniprot_positions)} residues",
                metrics={"rmsd_angstrom": structure_metrics.rmsd_angstrom},
            )
        )

        # --- 11. Confidence evaluation -------------------------------------
        start = prediction.sequence_start or 1
        plddt_by_uniprot: dict[int, float] = {}
        for uniprot_index, alphafold_index in alphafold_alignment.mapping.items():
            number = alphafold_index + start
            if number in validated_scores:
                plddt_by_uniprot[uniprot_index + 1] = validated_scores[number]
        if len(plddt_by_uniprot) != len(validated_scores):
            raise ConfidenceError(
                ErrorCode.CONFIDENCE_LENGTH_MISMATCH,
                f"only {len(plddt_by_uniprot)} of {len(validated_scores)} pLDDT values map "
                "onto the UniProt reference sequence",
                details={
                    "mapped": len(plddt_by_uniprot),
                    "available": len(validated_scores),
                },
            )

        trace_plddt = plddt_from_trace(predicted_trace)
        cross_check = cross_check_confidence(
            validated_scores,
            dict(zip(trace_plddt.residue_numbers, trace_plddt.scores, strict=True)),
        )
        if cross_check is not None and cross_check > 0.5:
            raise ConfidenceError(
                ErrorCode.INVALID_API_RESPONSE,
                f"pLDDT from the confidence document disagrees with the model file B-factor "
                f"column by {cross_check:.3f}; refusing to mix two different confidence sources",
                details={"max_abs_delta": cross_check},
            )
        if cross_check is not None and cross_check > 0.01:
            warnings.append(
                f"pLDDT from the confidence document differs from the model B-factor column by "
                f"up to {cross_check:.3f} (rounding expected)"
            )
        checks.append(
            CheckResult(
                check_id="plddt_matches_model_bfactors",
                title="Confidence JSON agrees with the pLDDT stored in the model file",
                status=CheckStatus.SKIPPED if cross_check is None else CheckStatus.PASSED,
                detail=(
                    "no shared residues to compare"
                    if cross_check is None
                    else f"largest |Δ| = {cross_check:.4f}"
                ),
                metrics={"max_abs_delta": cross_check},
            )
        )

        lddt_by_position = {
            position: float(lddt.per_residue[index])
            for index, position in enumerate(matching.uniprot_positions)
            if _finite(float(lddt.per_residue[index]))
        }
        confidence_evaluation = summarise_confidence(
            plddt_by_uniprot,
            lddt_by_position,
            source="alphafold_confidence_json",
            cross_check_max_abs_delta=cross_check,
        )
        checks.append(
            CheckResult(
                check_id="plddt_range_valid",
                title="pLDDT values lie in [0, 100] and cover the model",
                status=CheckStatus.PASSED,
                detail=f"min {confidence_evaluation.min_plddt:.2f}, "
                f"max {confidence_evaluation.max_plddt:.2f}, "
                f"{confidence_evaluation.residue_count} residues",
                metrics={
                    "min_plddt": confidence_evaluation.min_plddt,
                    "max_plddt": confidence_evaluation.max_plddt,
                },
            )
        )

        # --- 12. Assemble metrics ------------------------------------------
        structure_evaluation = self._build_structure_evaluation(
            matching, experimental_trace, predicted_trace, structure_metrics, plddt_by_uniprot
        )
        manifest.metrics.update(
            {
                "sequence": {
                    "uniprot_length": len(uniprot.sequence),
                    "alphafold_identity": alphafold_alignment.identity,
                    "alphafold_coverage": alphafold_alignment.coverage_query,
                    "pdb_entity_identity": entity_alignment.identity,
                    "pdb_entity_coverage_query": entity_alignment.coverage_query,
                    "pdb_entity_coverage_target": entity_alignment.coverage_target,
                    "pdb_entity_gap_openings": entity_alignment.gap_openings,
                    "pdb_entity_mismatches": entity_alignment.mismatches,
                },
                "structure": structure_evaluation.model_dump(mode="json"),
                "confidence": confidence_evaluation.model_dump(mode="json"),
                "limitations": LIMITATIONS,
                "not_evaluated": self._not_evaluated(pae_summary),
            }
        )
        if pae_summary is not None:
            manifest.metrics["pae"] = pae_summary.model_dump(mode="json")
        checks.append(
            CheckResult(
                check_id="pae_statistics",
                title="Predicted aligned error statistics extracted",
                status=CheckStatus.PASSED if pae_summary is not None else CheckStatus.SKIPPED,
                detail=(
                    f"mean {pae_summary.mean_predicted_aligned_error:.3f} Å"
                    if pae_summary is not None
                    else "PAE not available from this prediction record"
                ),
            )
        )
        checks.append(
            CheckResult(
                check_id="tm_score",
                title="TM-score comparison",
                status=CheckStatus.SKIPPED,
                detail=(
                    "not implemented in v0.1: no validated structure-alignment program is bundled"
                ),
            )
        )

    # -- helpers -----------------------------------------------------------
    def _build_structure_evaluation(
        self,
        matching: Any,
        experimental_trace: CaTrace,
        predicted_trace: CaTrace,
        metrics: Any,
        plddt_by_uniprot: dict[int, float],
    ) -> StructureEvaluation:
        lddt = metrics.lddt
        per_residue: list[ResidueScore] = []
        for index, (position, experimental_index, predicted_index) in enumerate(
            zip(
                matching.uniprot_positions,
                matching.experimental_indices,
                matching.predicted_indices,
                strict=True,
            )
        ):
            score = float(lddt.per_residue[index])
            per_residue.append(
                ResidueScore(
                    uniprot_position=position,
                    experimental_residue_number=experimental_trace.residue_numbers[
                        experimental_index
                    ],
                    predicted_sequence_position=predicted_trace.residue_numbers[predicted_index],
                    residue_name=experimental_trace.residue_names[experimental_index],
                    plddt=plddt_by_uniprot.get(position),
                    lddt_ca=score if _finite(score) else None,
                )
            )
        return StructureEvaluation(
            mapped_residue_count=len(matching.uniprot_positions),
            experimental_residue_count=len(experimental_trace),
            predicted_residue_count=len(predicted_trace),
            experimental_only_residues=matching.experimental_only_positions,
            predicted_only_residues=matching.predicted_only_positions,
            rmsd_angstrom=round(metrics.rmsd_angstrom, 6),
            lddt_ca=round(lddt.score, 6),
            lddt_ca_per_threshold={k: round(v, 6) for k, v in lddt.per_threshold.items()},
            lddt_ca_pairs=lddt.pairs,
            lddt_ca_residues_scored=lddt.scored_residues,
            lddt_ca_residues_unscored=lddt.unscored_residues,
            lddt_ca_inclusion_radius=lddt.radius,
            lddt_ca_thresholds=list(lddt.thresholds),
            per_residue=per_residue,
            notes=[
                "lDDT-Cα uses Cα atoms only; it is not a full-atom lDDT.",
                "RMSD uses the same Cα correspondence after Kabsch superposition.",
                "lDDT-Cα is superposition-free by definition; the superposition is only for RMSD.",
            ],
        )

    @staticmethod
    def _not_evaluated(pae: PaeSummary | None) -> list[dict[str, str]]:
        items = [
            {
                "metric": "full-atom lDDT",
                "reason": (
                    "v0.1 computes the Cα-only variant (lDDT-Cα); all-atom scoring is planned"
                ),
            },
            {
                "metric": "TM-score",
                "reason": "no validated structure-alignment program is bundled in v0.1",
            },
            {
                "metric": "AlphaFold inference",
                "reason": "predictions are read from AlphaFold DB; no model is run by this project",
            },
            {
                "metric": "oligomeric assembly assessment",
                "reason": "only the single requested chain is compared",
            },
        ]
        if pae is None:
            items.append(
                {"metric": "PAE analysis", "reason": "the prediction record provided no PAE data"}
            )
        else:
            items.append(
                {
                    "metric": "per-domain PAE interpretation",
                    "reason": "only aggregate PAE statistics are computed in v0.1",
                }
            )
        return items

    def _record_artifact(
        self,
        manifest: RunManifest,
        artifact: ArtifactRecord,
        run_dir: Path,
        path: Path,
    ) -> None:
        artifact.path = str(path.relative_to(run_dir))
        manifest.artifacts.append(artifact)
        self._artifacts.append(artifact)

    def _archive_responses(self, run_dir: Path) -> None:
        """Persist every fetched response byte-for-byte under ``raw/``."""
        if not self._settings.keep_raw_downloads:
            return
        for result in self._http.responses:
            url = result.provenance.url
            if url in self._archived_urls:
                continue
            basename = _url_filename(url, "response.bin")
            if not Path(basename).suffix:
                content_type = result.provenance.content_type or ""
                basename = f"{basename}.json" if "json" in content_type else basename
            filename = f"{result.provenance.source}__{basename}"
            artifact = write_bytes(run_dir / "raw" / filename, result.content)
            artifact.url = result.provenance.url
            artifact.kind = f"{result.provenance.source}_response"
            artifact.path = str((run_dir / "raw" / filename).relative_to(run_dir))
            self._artifacts.append(artifact)
            self._archived_urls.add(url)

    # -- finalisation ------------------------------------------------------
    def _finalize(self, manifest: RunManifest, run_dir: Path) -> dict[str, Any] | None:
        """Write metrics, Markdown, report.json and manifest.json.

        Ordering matters: everything the manifest hashes is written first, then
        report.json (which embeds the manifest digest), then the manifest itself
        so the digest in the report is the digest of the manifest on disk.
        """
        try:
            self._detach_logging()
            for artifact in list(self._artifacts):
                if artifact not in manifest.artifacts:
                    manifest.artifacts.append(artifact)

            metrics_artifact = write_metrics(run_dir, manifest)
            manifest.artifacts.append(metrics_artifact)

            markdown_artifact = markdown_report.write_markdown(run_dir, manifest)
            manifest.artifacts.append(markdown_artifact)

            log_path = run_dir / "logs" / "run.log"
            if log_path.exists():
                manifest.artifacts.append(artifact_for_file(log_path, run_dir, kind="run_log"))

            report_artifact, report = write_report(run_dir, manifest)
            manifest_artifact = write_manifest(run_dir, manifest)
            logger.info(
                "wrote %s (%d bytes), %s (%d bytes)",
                report_artifact.path,
                report_artifact.size_bytes,
                manifest_artifact.path,
                manifest_artifact.size_bytes,
            )
            return report
        except ReportingError as exc:
            logger.error("failed to write run artifacts: %s", exc)
            return None


def run_protein_smoke(
    *,
    uniprot: str,
    pdb_id: str,
    chain_id: str,
    output_dir: Path | str,
    settings: Settings | None = None,
    overwrite: bool = False,
    http: HttpClient | None = None,
) -> RunOutcome:
    """Convenience wrapper around :class:`ProteinSmokePipeline`."""
    pipeline = ProteinSmokePipeline(settings=settings, http=http)
    try:
        return pipeline.run(
            SmokeRequest(
                uniprot_accession=uniprot,
                pdb_id=pdb_id,
                chain_id=chain_id,
                output_dir=Path(output_dir),
                overwrite=overwrite,
            )
        )
    finally:
        pipeline.close()


__all__ = ["ProteinSmokePipeline", "RunOutcome", "SmokeRequest", "run_protein_smoke"]
