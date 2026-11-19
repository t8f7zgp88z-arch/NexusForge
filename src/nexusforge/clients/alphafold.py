"""AlphaFold Protein Structure Database client.

Important scientific boundary: this project **reads existing AlphaFold DB
predictions**. It does not run the AlphaFold model itself, and no code here
should be described as "running AlphaFold".

Schema discipline: the parser reads the field names the live API returns today
(``modelEntityId``, ``sequence``, ``cifUrl``, ``plddtDocUrl``, ``paeDocUrl``,
``latestVersion``, …). Retired-looking aliases are accepted only as ordered
fallbacks, and every fallback is recorded in ``schema_notes`` so a schema
migration shows up in the run report instead of silently breaking the run. No
file name or model version is hardcoded: URLs always come from the API record.
"""

from __future__ import annotations

from itertools import pairwise
from typing import Any

import numpy as np

from nexusforge.clients.common import (
    parse_datetime,
    parse_json_object,
    pick,
    require_float,
)
from nexusforge.errors import ConfidenceError, ErrorCode, SourceError
from nexusforge.http import HttpClient
from nexusforge.models import (
    AlphaFoldFiles,
    AlphaFoldPrediction,
    ConfidenceScores,
    PaeSummary,
    Provenance,
)

BASE_URL = "https://alphafold.ebi.ac.uk"

#: Field name candidates, most current first. Later entries are tolerated
#: fallbacks and are always reported in ``AlphaFoldPrediction.schema_notes``.
_FIELDS: dict[str, tuple[str, ...]] = {
    "entry_id": ("modelEntityId", "entryId"),
    "sequence": ("sequence", "uniprotSequence"),
    "sequence_start": ("sequenceStart", "uniprotStart"),
    "sequence_end": ("sequenceEnd", "uniprotEnd"),
    "reviewed": ("isUniProtReviewed", "isReviewed"),
    "cif_url": ("cifUrl",),
    "pdb_url": ("pdbUrl",),
    "plddt_doc_url": ("plddtDocUrl",),
    "pae_doc_url": ("paeDocUrl",),
    "pae_image_url": ("paeImageUrl",),
    "msa_url": ("msaUrl",),
    "global_plddt": ("globalMetricValue",),
    "latest_version": ("latestVersion",),
    "tool_used": ("toolUsed",),
    "created": ("modelCreatedDate",),
    "chain_id": ("chainId",),
    "is_complex": ("isComplex",),
}


class AlphaFoldClient:
    """Typed accessor for AlphaFold DB prediction records and their files."""

    def __init__(self, http: HttpClient, base_url: str = BASE_URL) -> None:
        self._http = http
        self._base_url = base_url.rstrip("/")

    @property
    def base_url(self) -> str:
        return self._base_url

    def prediction_url(self, accession: str) -> str:
        return f"{self._base_url}/api/prediction/{accession}"

    # -- retrieval ---------------------------------------------------------
    def fetch_prediction(self, accession: str) -> tuple[AlphaFoldPrediction, Provenance]:
        url = self.prediction_url(accession)
        result = self._http.get(url, source="alphafold")
        payload = parse_json_object(result.content, source="alphafold", url=url)
        prediction = self.parse_prediction(payload, accession=accession, url=url)
        return prediction, result.provenance

    def fetch_confidence(self, url: str) -> tuple[ConfidenceScores, Provenance]:
        result = self._http.get(url, source="alphafold-confidence")
        payload = parse_json_object(result.content, source="alphafold-confidence", url=url)
        scores = self.parse_confidence(payload, url=url)
        return scores, result.provenance

    def fetch_pae(self, url: str) -> tuple[PaeSummary, Provenance]:
        result = self._http.get(url, source="alphafold-pae")
        payload = parse_json_object(result.content, source="alphafold-pae", url=url)
        summary = self.parse_pae(payload, url=url)
        return summary, result.provenance

    # -- parsing -----------------------------------------------------------
    @staticmethod
    def select_entry(payload: Any, *, accession: str, url: str = "<payload>") -> dict[str, Any]:
        """Choose the monomer prediction that answers this accession."""
        if not isinstance(payload, list):
            raise SourceError(
                ErrorCode.INVALID_API_RESPONSE,
                "alphafold prediction endpoint did not return a JSON list",
                source="alphafold",
                details={"url": url, "actual_type": type(payload).__name__},
            )
        entries = [e for e in payload if isinstance(e, dict)]
        if not entries:
            raise SourceError(
                ErrorCode.SOURCE_NOT_FOUND,
                f"AlphaFold DB has no prediction entry for {accession}",
                source="alphafold",
                details={"url": url},
            )
        matching = [
            e for e in entries if str(e.get("uniprotAccession", "")).upper() == accession.upper()
        ]
        if not matching:
            raise SourceError(
                ErrorCode.SOURCE_NOT_FOUND,
                f"AlphaFold DB returned {len(entries)} entry/entries but none for {accession}",
                source="alphafold",
                details={
                    "url": url,
                    "returned_accessions": sorted(
                        {str(e.get("uniprotAccession")) for e in entries}
                    ),
                },
            )
        monomers = [e for e in matching if e.get("isComplex") is not True]
        candidates = monomers or matching
        # Highest model version wins; creation date breaks ties deterministically.
        candidates.sort(
            key=lambda e: (
                int(e.get("latestVersion") or 0),
                str(e.get("modelCreatedDate") or ""),
            ),
            reverse=True,
        )
        return candidates[0]

    @classmethod
    def parse_prediction(
        cls, payload: Any, *, accession: str, url: str = "<payload>"
    ) -> AlphaFoldPrediction:
        entry = cls.select_entry(payload, accession=accession, url=url)
        notes: list[str] = []
        entry_id = pick(entry, _FIELDS["entry_id"], source="alphafold", url=url, notes=notes)
        sequence = pick(entry, _FIELDS["sequence"], source="alphafold", url=url, notes=notes)
        if not isinstance(sequence, str) or not sequence.strip():
            raise SourceError(
                ErrorCode.INVALID_API_RESPONSE,
                "alphafold entry has an empty or non-string sequence",
                source="alphafold",
                details={"url": url},
            )
        sequence = sequence.strip().upper()

        files = AlphaFoldFiles(
            cif_url=pick(entry, _FIELDS["cif_url"], source="alphafold", url=url, required=False),
            pdb_url=pick(entry, _FIELDS["pdb_url"], source="alphafold", url=url, required=False),
            plddt_doc_url=pick(
                entry,
                _FIELDS["plddt_doc_url"],
                source="alphafold",
                url=url,
                required=False,
            ),
            pae_doc_url=pick(
                entry, _FIELDS["pae_doc_url"], source="alphafold", url=url, required=False
            ),
            pae_image_url=pick(
                entry, _FIELDS["pae_image_url"], source="alphafold", url=url, required=False
            ),
            msa_url=pick(entry, _FIELDS["msa_url"], source="alphafold", url=url, required=False),
        )
        if not files.cif_url and not files.pdb_url:
            raise SourceError(
                ErrorCode.INVALID_API_RESPONSE,
                "alphafold entry advertises neither a mmCIF nor a PDB structure file",
                source="alphafold",
                details={"url": url, "fields": sorted(entry.keys())},
            )

        def _int(name: str) -> int | None:
            value = entry.get(name)
            return int(value) if isinstance(value, int) and not isinstance(value, bool) else None

        def _frac(name: str) -> float | None:
            value = entry.get(name)
            if isinstance(value, (int, float)) and not isinstance(value, bool):
                return float(value)
            return None

        return AlphaFoldPrediction(
            entry_id=str(entry_id),
            model_entity_id=entry.get("modelEntityId")
            if isinstance(entry.get("modelEntityId"), str)
            else None,
            uniprot_accession=str(entry.get("uniprotAccession", accession)).upper(),
            tool_used=pick(
                entry, _FIELDS["tool_used"], source="alphafold", url=url, required=False
            ),
            latest_version=_int("latestVersion"),
            model_created_date=parse_datetime(
                pick(entry, _FIELDS["created"], source="alphafold", url=url, required=False)
            ),
            global_plddt=require_float(
                pick(entry, _FIELDS["global_plddt"], source="alphafold", url=url, required=False),
                field="globalMetricValue",
                source="alphafold",
            ),
            fraction_plddt_very_low=_frac("fractionPlddtVeryLow"),
            fraction_plddt_low=_frac("fractionPlddtLow"),
            fraction_plddt_confident=_frac("fractionPlddtConfident"),
            fraction_plddt_very_high=_frac("fractionPlddtVeryHigh"),
            sequence=sequence,
            sequence_start=_int("sequenceStart"),
            sequence_end=_int("sequenceEnd"),
            chain_id=pick(entry, _FIELDS["chain_id"], source="alphafold", url=url, required=False),
            is_complex=pick(
                entry, _FIELDS["is_complex"], source="alphafold", url=url, required=False
            ),
            files=files,
            schema_notes=notes,
        )

    @staticmethod
    def parse_confidence(payload: Any, *, url: str = "<payload>") -> ConfidenceScores:
        """Parse the per-residue confidence document into validated pLDDT scores."""
        residue_numbers: list[int]
        scores: list[float]

        if isinstance(payload, dict):
            raw_numbers = payload.get("residueNumber")
            raw_scores = payload.get("confidenceScore")
            if raw_numbers is None or raw_scores is None:
                raise SourceError(
                    ErrorCode.INVALID_API_RESPONSE,
                    "alphafold confidence document lacks residueNumber/confidenceScore",
                    source="alphafold-confidence",
                    details={"url": url, "keys": sorted(payload.keys())},
                )
            if not isinstance(raw_numbers, list) or not isinstance(raw_scores, list):
                raise SourceError(
                    ErrorCode.INVALID_API_RESPONSE,
                    "alphafold confidence fields must be JSON arrays",
                    source="alphafold-confidence",
                    details={"url": url},
                )
            residue_numbers = [int(n) for n in raw_numbers]
            scores = [float(s) for s in raw_scores]
        elif isinstance(payload, list) and payload and isinstance(payload[0], dict):
            residue_numbers = [int(item["residueNumber"]) for item in payload]
            scores = [float(item["confidenceScore"]) for item in payload]
        else:
            raise SourceError(
                ErrorCode.INVALID_API_RESPONSE,
                "alphafold confidence document has an unexpected shape",
                source="alphafold-confidence",
                details={"url": url, "actual_type": type(payload).__name__},
            )

        if len(residue_numbers) != len(scores):
            raise ConfidenceError(
                ErrorCode.CONFIDENCE_LENGTH_MISMATCH,
                f"confidence arrays disagree: {len(residue_numbers)} residue numbers vs "
                f"{len(scores)} scores",
                details={"url": url},
            )
        if not scores:
            raise ConfidenceError(
                ErrorCode.CONFIDENCE_LENGTH_MISMATCH,
                "confidence document contains no residues",
                details={"url": url},
            )
        if any(not np.isfinite(s) for s in scores):
            raise SourceError(
                ErrorCode.INVALID_API_RESPONSE,
                "confidence document contains non-finite pLDDT values",
                source="alphafold-confidence",
                details={"url": url},
            )
        out_of_range = [s for s in scores if s < 0.0 or s > 100.0]
        if out_of_range:
            raise SourceError(
                ErrorCode.INVALID_API_RESPONSE,
                f"{len(out_of_range)} pLDDT value(s) fall outside [0, 100]; "
                f"first offending value {out_of_range[0]}",
                source="alphafold-confidence",
                details={"url": url},
            )
        if any(b <= a for a, b in pairwise(residue_numbers)):
            raise SourceError(
                ErrorCode.INVALID_API_RESPONSE,
                "confidence residue numbers are not strictly increasing",
                source="alphafold-confidence",
                details={"url": url},
            )
        return ConfidenceScores(
            source="alphafold_confidence_json", residue_numbers=residue_numbers, scores=scores
        )

    @staticmethod
    def parse_pae(payload: Any, *, url: str = "<payload>") -> PaeSummary:
        entry: Any = payload
        if isinstance(payload, list):
            if not payload:
                raise SourceError(
                    ErrorCode.INVALID_API_RESPONSE,
                    "alphafold PAE document is an empty list",
                    source="alphafold-pae",
                    details={"url": url},
                )
            entry = payload[0]
        if not isinstance(entry, dict):
            raise SourceError(
                ErrorCode.INVALID_API_RESPONSE,
                "alphafold PAE document has an unexpected shape",
                source="alphafold-pae",
                details={"url": url},
            )
        matrix = entry.get("predicted_aligned_error")
        if not isinstance(matrix, list):
            raise SourceError(
                ErrorCode.INVALID_API_RESPONSE,
                "alphafold PAE document lacks a predicted_aligned_error matrix",
                source="alphafold-pae",
                details={"url": url},
            )
        array: np.ndarray = np.asarray(matrix, dtype=float)
        shape = tuple(int(dimension) for dimension in array.shape)
        if len(shape) != 2 or shape[0] != shape[1]:
            raise SourceError(
                ErrorCode.INVALID_API_RESPONSE,
                f"PAE matrix is not square (shape {shape})",
                source="alphafold-pae",
                details={"url": url},
            )
        if not np.isfinite(array).all():
            raise SourceError(
                ErrorCode.INVALID_API_RESPONSE,
                "PAE matrix contains non-finite values",
                source="alphafold-pae",
                details={"url": url},
            )
        max_error = require_float(
            entry.get("max_predicted_aligned_error"),
            field="max_predicted_aligned_error",
            source="alphafold-pae",
        )
        return PaeSummary(
            shape=(shape[0], shape[1]),
            max_predicted_aligned_error=max_error,
            mean_predicted_aligned_error=float(array.mean()),
            p90_predicted_aligned_error=float(np.percentile(array, 90)),
        )


__all__ = ["BASE_URL", "AlphaFoldClient"]
