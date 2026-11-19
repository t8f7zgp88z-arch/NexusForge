"""UniProt REST client (https://rest.uniprot.org).

Only the UniProtKB JSON endpoint is used. The accession is validated locally
first so a typo fails as an input error instead of a confusing HTTP 400.
"""

from __future__ import annotations

import re
from typing import Any

from nexusforge.clients.common import (
    parse_json_object,
    pick,
    require_mapping,
)
from nexusforge.errors import ErrorCode, InputError, SourceError
from nexusforge.http import HttpClient
from nexusforge.models import Provenance, UniProtRecord

BASE_URL = "https://rest.uniprot.org"

#: Canonical UniProtKB accession shapes (6 or 10 characters).
ACCESSION_PATTERN = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$"
)


def validate_accession(accession: str) -> str:
    """Normalise and validate a UniProt accession, or raise :class:`InputError`."""
    candidate = (accession or "").strip().upper()
    if not ACCESSION_PATTERN.fullmatch(candidate):
        raise InputError(
            ErrorCode.INVALID_INPUT,
            f"{accession!r} is not a valid UniProtKB accession "
            "(expected e.g. P69905 or A0A123B4C5)",
        )
    return candidate


def _full_name(block: Any) -> str | None:
    if not isinstance(block, dict):
        return None
    full = block.get("fullName")
    if isinstance(full, dict) and isinstance(full.get("value"), str):
        return full["value"]
    return None


class UniProtClient:
    """Minimal, typed accessor for a single UniProtKB entry."""

    def __init__(self, http: HttpClient, base_url: str = BASE_URL) -> None:
        self._http = http
        self._base_url = base_url.rstrip("/")

    @property
    def base_url(self) -> str:
        return self._base_url

    def entry_url(self, accession: str) -> str:
        return f"{self._base_url}/uniprotkb/{accession}.json"

    def fetch(self, accession: str) -> tuple[UniProtRecord, Provenance]:
        accession = validate_accession(accession)
        url = self.entry_url(accession)
        result = self._http.get(url, source="uniprot")
        payload = parse_json_object(result.content, source="uniprot", url=url)
        record = self.parse(payload, requested_accession=accession, url=url)
        return record, result.provenance

    @staticmethod
    def parse(payload: Any, *, requested_accession: str, url: str = "<payload>") -> UniProtRecord:
        """Parse a UniProtKB entry payload into :class:`UniProtRecord`."""
        if not isinstance(payload, dict):
            raise SourceError(
                ErrorCode.INVALID_API_RESPONSE,
                "uniprot returned a JSON document that is not an entry object",
                source="uniprot",
                details={"url": url, "actual_type": type(payload).__name__},
            )
        if "results" in payload:
            # A search-style payload arriving where an entry was expected.
            raise SourceError(
                ErrorCode.INVALID_API_RESPONSE,
                "uniprot returned a search result list instead of a single entry",
                source="uniprot",
                details={"url": url},
            )

        notes: list[str] = []
        accession = pick(payload, ["primaryAccession"], source="uniprot", url=url)
        if not isinstance(accession, str) or accession.upper() != requested_accession:
            raise SourceError(
                ErrorCode.INVALID_API_RESPONSE,
                f"uniprot returned entry {accession!r} for requested {requested_accession!r}",
                source="uniprot",
                details={"url": url},
            )

        sequence_block = require_mapping(
            pick(payload, ["sequence"], source="uniprot", url=url),
            field="sequence",
            source="uniprot",
            url=url,
        )
        sequence = pick(sequence_block, ["value"], source="uniprot", url=url)
        length = pick(sequence_block, ["length"], source="uniprot", url=url)
        if not isinstance(sequence, str) or not isinstance(length, int):
            raise SourceError(
                ErrorCode.INVALID_API_RESPONSE,
                "uniprot sequence block is missing a string value or integer length",
                source="uniprot",
                details={"url": url},
            )
        if len(sequence) != length:
            raise SourceError(
                ErrorCode.INVALID_API_RESPONSE,
                f"uniprot sequence length {len(sequence)} disagrees with declared length {length}",
                source="uniprot",
                details={"url": url},
            )

        protein_description = payload.get("proteinDescription")
        protein_name = None
        alternative_names: list[str] = []
        if isinstance(protein_description, dict):
            recommended = protein_description.get("recommendedName")
            protein_name = _full_name(recommended)
            if protein_name is None:
                submission_names = protein_description.get("submissionNames")
                if isinstance(submission_names, list) and submission_names:
                    protein_name = _full_name(submission_names[0])
            for alternative in protein_description.get("alternativeNames") or []:
                name = _full_name(alternative)
                if name:
                    alternative_names.append(name)

        gene_names: list[str] = []
        genes = payload.get("genes")
        if isinstance(genes, list):
            for gene in genes:
                if not isinstance(gene, dict):
                    continue
                gene_name = gene.get("geneName")
                if isinstance(gene_name, dict) and isinstance(gene_name.get("value"), str):
                    value = gene_name["value"]
                    if value not in gene_names:
                        gene_names.append(value)

        organism_block = require_mapping(
            pick(payload, ["organism"], source="uniprot", url=url),
            field="organism",
            source="uniprot",
            url=url,
        )
        organism = pick(organism_block, ["scientificName"], source="uniprot", url=url)
        taxon_id = organism_block.get("taxonId")
        entry_type = pick(payload, ["entryType"], source="uniprot", url=url)
        entry_name = pick(payload, ["uniProtkbId"], source="uniprot", url=url, notes=notes)

        return UniProtRecord(
            accession=accession.upper(),
            entry_name=str(entry_name),
            reviewed="reviewed" in str(entry_type).lower(),
            protein_name=protein_name,
            alternative_names=alternative_names,
            gene_names=gene_names,
            organism=str(organism),
            taxon_id=int(taxon_id) if isinstance(taxon_id, int) else None,
            sequence=sequence,
            sequence_length=length,
            sequence_md5=(
                sequence_block.get("md5") if isinstance(sequence_block.get("md5"), str) else None
            ),
            sequence_crc64=(
                sequence_block.get("crc64")
                if isinstance(sequence_block.get("crc64"), str)
                else None
            ),
        )


__all__ = ["ACCESSION_PATTERN", "BASE_URL", "UniProtClient", "validate_accession"]
