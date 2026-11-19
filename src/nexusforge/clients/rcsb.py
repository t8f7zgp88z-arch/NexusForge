"""RCSB Protein Data Bank client.

Two services are used, both public and documented:

* ``https://data.rcsb.org/rest/v1/core/...`` — entry and polymer-entity metadata.
* ``https://files.rcsb.org/download/<PDB>.cif`` — the experimental coordinates.

The experimental chain is *verified* against the UniProt accession through the
entity's UniProt cross-references. Nothing is assumed to match just because a
configuration file says so: a mismatch, a non-polymer chain, or an ambiguous
chain-to-entity relationship is a hard, typed failure.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

from nexusforge.clients.common import (
    parse_datetime,
    parse_json_object,
    pick,
    require_float,
    require_mapping,
)
from nexusforge.errors import ErrorCode, InputError, SourceError, ValidationError
from nexusforge.http import FetchResult, HttpClient
from nexusforge.models import (
    ArtifactRecord,
    Provenance,
    RcsbChainMapping,
    RcsbEntry,
    RcsbPolymerEntity,
)

DATA_BASE_URL = "https://data.rcsb.org/rest/v1/core"
FILES_BASE_URL = "https://files.rcsb.org/download"

PDB_ID_PATTERN = re.compile(r"^[0-9][A-Za-z0-9]{3}$")
CHAIN_ID_PATTERN = re.compile(r"^[A-Za-z0-9]{1,2}$")


def validate_pdb_id(pdb_id: str) -> str:
    candidate = (pdb_id or "").strip().upper()
    if not PDB_ID_PATTERN.fullmatch(candidate):
        raise InputError(
            ErrorCode.INVALID_INPUT,
            f"{pdb_id!r} is not a valid PDB identifier (expected four characters, e.g. 4HHB)",
        )
    return candidate


def validate_chain_id(chain_id: str) -> str:
    candidate = (chain_id or "").strip()
    if not CHAIN_ID_PATTERN.fullmatch(candidate):
        raise InputError(
            ErrorCode.INVALID_INPUT,
            f"{chain_id!r} is not a valid PDB author chain identifier",
        )
    return candidate.upper()


class RcsbClient:
    """Accessor for experimental entry metadata, chain mapping and coordinates."""

    def __init__(
        self,
        http: HttpClient,
        data_base_url: str = DATA_BASE_URL,
        files_base_url: str = FILES_BASE_URL,
    ) -> None:
        self._http = http
        self._data_base_url = data_base_url.rstrip("/")
        self._files_base_url = files_base_url.rstrip("/")

    # -- URLs --------------------------------------------------------------
    def entry_url(self, pdb_id: str) -> str:
        return f"{self._data_base_url}/entry/{pdb_id}"

    def entity_url(self, pdb_id: str, entity_id: str) -> str:
        return f"{self._data_base_url}/polymer_entity/{pdb_id}/{entity_id}"

    def structure_url(self, pdb_id: str) -> str:
        return f"{self._files_base_url}/{pdb_id}.cif"

    # -- retrieval ---------------------------------------------------------
    def fetch_entry(self, pdb_id: str) -> tuple[RcsbEntry, Provenance]:
        pdb_id = validate_pdb_id(pdb_id)
        url = self.entry_url(pdb_id)
        result = self._http.get(url, source="rcsb-entry")
        payload = parse_json_object(result.content, source="rcsb-entry", url=url)
        return self.parse_entry(payload, requested_pdb_id=pdb_id, url=url), result.provenance

    def fetch_polymer_entity(
        self, pdb_id: str, entity_id: str
    ) -> tuple[RcsbPolymerEntity, Provenance]:
        pdb_id = validate_pdb_id(pdb_id)
        url = self.entity_url(pdb_id, entity_id)
        result = self._http.get(url, source="rcsb-entity")
        payload = parse_json_object(result.content, source="rcsb-entity", url=url)
        return (
            self.parse_polymer_entity(
                payload, requested_pdb_id=pdb_id, entity_id=entity_id, url=url
            ),
            result.provenance,
        )

    def resolve_chain(
        self,
        pdb_id: str,
        chain_id: str,
        uniprot_accession: str,
        *,
        entry: RcsbEntry | None = None,
    ) -> tuple[RcsbChainMapping, list[Provenance]]:
        """Bind an author chain to the UniProt accession, or fail explicitly."""
        pdb_id = validate_pdb_id(pdb_id)
        chain_id = validate_chain_id(chain_id)
        accession = uniprot_accession.strip().upper()

        provenance: list[Provenance] = []
        if entry is None:
            entry, entry_provenance = self.fetch_entry(pdb_id)
            provenance.append(entry_provenance)
        if entry.pdb_id.upper() != pdb_id:
            raise SourceError(
                ErrorCode.INVALID_API_RESPONSE,
                f"entry metadata is for {entry.pdb_id}, not {pdb_id}",
                source="rcsb-entry",
            )
        if not entry.polymer_entity_ids:
            raise SourceError(
                ErrorCode.INVALID_API_RESPONSE,
                f"RCSB entry {pdb_id} lists no polymer entities",
                source="rcsb-entry",
                details={"url": self.entry_url(pdb_id)},
            )

        entities: list[RcsbPolymerEntity] = []
        for entity_id in entry.polymer_entity_ids:
            entity, entity_provenance = self.fetch_polymer_entity(pdb_id, entity_id)
            provenance.append(entity_provenance)
            entities.append(entity)

        all_polymer_chains = sorted({chain for e in entities for chain in e.auth_asym_ids})
        matches = [e for e in entities if chain_id in e.auth_asym_ids]

        if not matches:
            raise ValidationError(
                ErrorCode.SEQUENCE_MISMATCH,
                f"author chain {chain_id!r} is not a polymer chain of {pdb_id}; "
                f"polymer chains are {all_polymer_chains or 'unknown'}",
                details={
                    "pdb_id": pdb_id,
                    "chain_id": chain_id,
                    "polymer_chains": all_polymer_chains,
                },
            )

        if len(matches) > 1:
            raise ValidationError(
                ErrorCode.CHAIN_MAPPING_AMBIGUOUS,
                f"author chain {chain_id!r} of {pdb_id} appears in "
                f"{len(matches)} polymer entities ({[m.entity_id for m in matches]}); "
                "refusing to guess which entity is the intended target",
                details={
                    "pdb_id": pdb_id,
                    "chain_id": chain_id,
                    "entity_ids": [m.entity_id for m in matches],
                },
            )

        entity = matches[0]
        if accession not in [a.upper() for a in entity.uniprot_accessions]:
            raise ValidationError(
                ErrorCode.SEQUENCE_MISMATCH,
                f"chain {chain_id!r} of {pdb_id} maps to UniProt "
                f"{entity.uniprot_accessions or 'no accession'}, not {accession}",
                details={
                    "pdb_id": pdb_id,
                    "chain_id": chain_id,
                    "entity_id": entity.entity_id,
                    "uniprot_accessions": entity.uniprot_accessions,
                    "requested_accession": accession,
                },
            )

        mapping = RcsbChainMapping(
            pdb_id=pdb_id,
            entity_id=entity.entity_id,
            auth_chain_id=chain_id,
            asym_ids=entity.asym_ids,
            uniprot_accessions=[a.upper() for a in entity.uniprot_accessions],
            entity_name=entity.entity_name,
            sequence=entity.sequence,
            coverage_of_reference=entity.coverage_of_reference,
        )
        return mapping, provenance

    def download_structure(
        self, pdb_id: str, destination: Path
    ) -> tuple[ArtifactRecord, FetchResult]:
        pdb_id = validate_pdb_id(pdb_id)
        url = self.structure_url(pdb_id)
        return self._http.download(
            url,
            destination,
            source="rcsb-file",
            kind="experimental_structure_cif",
            accept="chemical/x-mmcif,text/plain,*/*",
        )

    # -- parsing -----------------------------------------------------------
    @staticmethod
    def parse_entry(payload: Any, *, requested_pdb_id: str, url: str = "<payload>") -> RcsbEntry:
        if not isinstance(payload, dict):
            raise SourceError(
                ErrorCode.INVALID_API_RESPONSE,
                "rcsb entry response is not a JSON object",
                source="rcsb-entry",
                details={"url": url, "actual_type": type(payload).__name__},
            )
        pdb_id = payload.get("rcsb_id")
        if not isinstance(pdb_id, str) or pdb_id.upper() != requested_pdb_id:
            raise SourceError(
                ErrorCode.INVALID_API_RESPONSE,
                f"rcsb returned entry {pdb_id!r} for requested {requested_pdb_id!r}",
                source="rcsb-entry",
                details={"url": url},
            )

        struct = require_mapping(
            pick(payload, ["struct"], source="rcsb-entry", url=url),
            field="struct",
            source="rcsb-entry",
            url=url,
        )
        title = str(
            pick(struct, ["title"], source="rcsb-entry", url=url, required=False) or ""
        ).strip()

        methods: list[str] = []
        for experiment in payload.get("exptl") or []:
            if isinstance(experiment, dict) and isinstance(experiment.get("method"), str):
                methods.append(experiment["method"])

        resolution: float | None = None
        refine = payload.get("refine")
        if isinstance(refine, list):
            for block in refine:
                if isinstance(block, dict):
                    resolution = require_float(
                        block.get("ls_d_res_high"), field="ls_d_res_high", source="rcsb-entry"
                    )
                    if resolution is not None:
                        break
        if resolution is None:
            entry_info = payload.get("rcsb_entry_info")
            if isinstance(entry_info, dict):
                combined = entry_info.get("resolution_combined")
                if isinstance(combined, list) and combined:
                    resolution = require_float(
                        combined[0], field="resolution_combined", source="rcsb-entry"
                    )

        accession_info = payload.get("rcsb_accession_info")
        deposited = released = revision = None
        if isinstance(accession_info, dict):
            deposited = parse_datetime(accession_info.get("deposit_date"))
            released = parse_datetime(accession_info.get("initial_release_date"))
            revision = parse_datetime(accession_info.get("revision_date"))

        container = payload.get("rcsb_entry_container_identifiers")
        entity_ids: list[str] = []
        if isinstance(container, dict):
            raw_ids = container.get("polymer_entity_ids")
            if isinstance(raw_ids, list):
                entity_ids = [str(e) for e in raw_ids]

        entry_info = payload.get("rcsb_entry_info")
        atom_count = None
        if isinstance(entry_info, dict):
            raw_count = entry_info.get("deposited_atom_count")
            if isinstance(raw_count, int) and not isinstance(raw_count, bool):
                atom_count = raw_count

        authors: list[str] = []
        for author in payload.get("audit_author") or []:
            if isinstance(author, dict) and isinstance(author.get("name"), str):
                authors.append(author["name"])

        return RcsbEntry(
            pdb_id=pdb_id.upper(),
            title=title,
            experimental_methods=methods,
            resolution_angstrom=resolution,
            deposited_at=deposited,
            released_at=released,
            last_revision_at=revision,
            polymer_entity_ids=entity_ids,
            deposited_atom_count=atom_count,
            authors=authors,
        )

    @staticmethod
    def parse_polymer_entity(
        payload: Any, *, requested_pdb_id: str, entity_id: str, url: str = "<payload>"
    ) -> RcsbPolymerEntity:
        if not isinstance(payload, dict):
            raise SourceError(
                ErrorCode.INVALID_API_RESPONSE,
                "rcsb polymer entity response is not a JSON object",
                source="rcsb-entity",
                details={"url": url, "actual_type": type(payload).__name__},
            )
        entity_poly = require_mapping(
            pick(payload, ["entity_poly"], source="rcsb-entity", url=url),
            field="entity_poly",
            source="rcsb-entity",
            url=url,
        )
        sequence = pick(
            entity_poly,
            ["pdbx_seq_one_letter_code_can", "pdbx_seq_one_letter_code"],
            source="rcsb-entity",
            url=url,
        )
        if not isinstance(sequence, str):
            raise SourceError(
                ErrorCode.INVALID_API_RESPONSE,
                "rcsb entity one-letter sequence is not a string",
                source="rcsb-entity",
                details={"url": url},
            )

        container = require_mapping(
            pick(
                payload,
                ["rcsb_polymer_entity_container_identifiers"],
                source="rcsb-entity",
                url=url,
            ),
            field="rcsb_polymer_entity_container_identifiers",
            source="rcsb-entity",
            url=url,
        )
        auth_asym_ids = [str(c) for c in container.get("auth_asym_ids") or []]
        asym_ids = [str(c) for c in container.get("asym_ids") or []]
        uniprot_ids = [str(u).upper() for u in container.get("uniprot_ids") or []]

        entity_name = None
        entity_block = payload.get("rcsb_polymer_entity")
        if isinstance(entity_block, dict) and isinstance(entity_block.get("pdbx_description"), str):
            entity_name = entity_block["pdbx_description"]

        coverage: float | None = None
        for reference in container.get("reference_sequence_identifiers") or []:
            if not isinstance(reference, dict):
                continue
            if str(reference.get("database_name", "")).lower() != "uniprot":
                continue
            accessions = {
                str(reference.get("database_accession", "")).upper(),
            }
            if accessions & set(uniprot_ids):
                coverage = require_float(
                    reference.get("reference_sequence_coverage"),
                    field="reference_sequence_coverage",
                    source="rcsb-entity",
                )
                break

        return RcsbPolymerEntity(
            pdb_id=requested_pdb_id.upper(),
            entity_id=str(container.get("entity_id", entity_id)),
            asym_ids=asym_ids,
            auth_asym_ids=auth_asym_ids,
            uniprot_accessions=uniprot_ids,
            entity_name=entity_name,
            sequence=sequence,
            coverage_of_reference=coverage,
        )


__all__ = [
    "CHAIN_ID_PATTERN",
    "DATA_BASE_URL",
    "FILES_BASE_URL",
    "PDB_ID_PATTERN",
    "RcsbClient",
    "validate_chain_id",
    "validate_pdb_id",
]
