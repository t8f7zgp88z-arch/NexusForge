# Data sources — what was actually observed

All statements below were verified against the live public APIs on
**2026-11-19 (UTC)** while building v0.1. They are recorded so a future schema
change can be diagnosed quickly.

## UniProt REST

* Endpoint: `GET https://rest.uniprot.org/uniprotkb/P69905.json` → `200`,
  241,776 bytes.
* Fields used: `primaryAccession`, `entryType`, `uniProtkbId`,
  `proteinDescription.recommendedName.fullName.value`, `organism.scientificName`,
  `organism.taxonId`, `genes[].geneName.value`, `sequence.{value,length,md5,crc64}`.
* A malformed accession is rejected locally (`INVALID_INPUT`) before any request,
  so a typo never becomes an ambiguous HTTP 400.
* `404` means "no such entry" (`SOURCE_NOT_FOUND`).

## AlphaFold DB

* Prediction endpoint:
  `GET https://alphafold.ebi.ac.uk/api/prediction/P69905` → `200`, a JSON **list**
  with one monomer record (1,987 bytes for P69905).
* File endpoint (mmCIF, confidence JSON, PAE JSON):
  `https://alphafold.ebi.ac.uk/files/...` → **HTTP 403 for a generic client
  identity**. A descriptive `User-Agent` (default
  `NexusForge/<version> (...)`) returns `200`. This is why the HTTP client always
  sends one.
* Observed record fields (P69905, model v6):
  `modelEntityId` = `entryId` = `AF-P69905-F1`, `toolUsed` =
  "AlphaFold Monomer v2.0 pipeline", `latestVersion` 6,
  `modelCreatedDate` 2025-08-01, `globalMetricValue` 98.06,
  `fractionPlddtVeryLow/Low/Confident/VeryHigh`, `sequence` (142 aa),
  `sequenceStart` 1, `sequenceEnd` 142, `chainId` "A", `isComplex` false,
  `cifUrl`, `pdbUrl`, `plddtDocUrl`, `paeDocUrl`, `msaUrl`.
* **Current vs legacy fields.** The live record still contains `entryId`,
  `uniprotSequence`, `isReviewed` and `paeImageUrl` alongside the newer
  `modelEntityId`/`sequence`/`isUniProtReviewed`. NexusForge prefers the current
  field names and treats the others as fallbacks; if a fallback is used it is
  written into `schema_notes` and surfaced as a run warning. No URL or model
  version is hard-coded, so `v6 → v7` requires no code change.
* Confidence document shape:
  `{"residueNumber": [...], "confidenceScore": [...], "confidenceCategory": [...]}`
  — 142 values in `[0, 100]` for P69905.
* PAE document shape: `[{"predicted_aligned_error": [[...]], "max_predicted_aligned_error": 31.75}]`
  — a 142×142 matrix; v0.1 reports aggregate statistics only.

## RCSB PDB

* Entry metadata: `GET https://data.rcsb.org/rest/v1/core/entry/4HHB` → `200`,
  26,055 bytes. Used: `struct.title`, `exptl[].method`,
  `refine[].ls_d_res_high` (1.74), `rcsb_entry_info.deposited_atom_count`,
  `rcsb_entry_container_identifiers.polymer_entity_ids` (`["1","2"]`),
  `rcsb_accession_info` dates.
* Polymer entity: `GET .../polymer_entity/4HHB/1` → entity `1`,
  `auth_asym_ids` `["A","C"]`, `uniprot_ids` `["P69905"]`, sequence length 141
  (initiator methionine not present), SIFTS `reference_sequence_coverage`
  0.99296.
* Coordinates: `GET https://files.rcsb.org/download/4HHB.cif` → `200`,
  772,198 bytes. Chain A contains 141 amino-acid residues plus HEM/water
  heteroatoms; only standard polymer residues with a Cα are used.
* Chain B belongs to entity 2 (`P68871`, beta-globin). Requesting
  `P69905 + 4HHB + B` therefore fails with `SEQUENCE_MISMATCH` — this is used as
  a real negative test.

## Reproducibility contract

For every response the run records: URL, HTTP status, UTC fetch time, byte size
and SHA-256; the bytes themselves are stored under `raw/`. Upstream files can
change, so reproduction means *re-fetch and compare digests*, not assuming
immutability.

The offline test fixtures under `tests/fixtures/` are gzipped copies of exactly
these responses; `tests/fixtures/SHA256SUMS.txt` records the digest of each
uncompressed payload and its capture date.
