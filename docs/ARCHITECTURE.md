# Architecture

## Design goals

1. **Verified inputs before metrics.** No structural number is produced until the
   chain, the model sequence and the UniProt entry are provably related.
2. **Typed, staged failures.** Every failure names *what* broke and *where*.
3. **Reproducibility as an artifact.** A run leaves hashed inputs plus the
   software versions needed to re-derive its numbers.
4. **No silent substitution.** Missing data is reported, never zero-filled.

## Data flow

```
CLI (cli.py)
  │  SmokeRequest(accession, pdb_id, chain_id, output_dir)
  ▼
ProteinSmokePipeline.run()                      pipeline/protein_smoke.py
  ├─ 1 data_retrieval
  │    UniProtClient.fetch                 → UniProtRecord + Provenance
  │    AlphaFoldClient.fetch_prediction    → AlphaFoldPrediction + Provenance
  │    AlphaFoldClient.fetch_confidence    → ConfidenceScores + Provenance
  │    AlphaFoldClient.fetch_pae           → PaeSummary + Provenance   (optional)
  │    RcsbClient.fetch_entry/resolve_chain→ RcsbEntry, RcsbChainMapping
  │    RcsbClient.download_structure       → mmCIF file (+ SHA-256)
  ├─ 2 sequence_validation
  │    align(uniprot, alphafold) / align(uniprot, pdb_entity)
  │    require_compatible(identity, coverage)
  ├─ 3 structure_parsing
  │    extract_ca_trace(experimental, chain) / extract_ca_trace(prediction)
  │    match_residues via UniProt reference frame
  ├─ 4 structure_evaluation
  │    kabsch_superpose → RMSD
  │    compute_lddt_ca  → lDDT-Cα (superposition-free)
  ├─ 5 confidence_evaluation
  │    validate_confidence → pLDDT range/count/numbering
  │    plddt_from_trace   → B-factor cross-check
  │    summarise_confidence → bins + Pearson/Spearman vs lDDT-Cα
  └─ 6 reporting
       metrics.json, report.md, report.json, manifest.json, logs/run.log
```

Failures short-circuit to `FAILED`; everything already fetched is still written.

## Module responsibilities

| Module | Owns | Must not |
| --- | --- | --- |
| `http.py` | retry policy, timeouts, provenance, atomic downloads | interpret domain payloads |
| `clients/common.py` | JSON decode errors, ordered field fallbacks, datetime/float parsing | know biology |
| `clients/uniprot.py` | accession validation, entry parsing | compute metrics |
| `clients/alphafold.py` | prediction selection, confidence/PAE validation | invent URLs or versions |
| `clients/rcsb.py` | entry/entity parsing, chain ↔ UniProt enforcement | align sequences |
| `validation/sequence.py` | sequence validity, global alignment, identity/coverage gates | parse structures |
| `validation/structure.py` | Cα extraction, residue matching, superposition, lDDT-Cα | read JSON APIs |
| `validation/confidence.py` | pLDDT validation, binning, correlation | call lDDT "pLDDT" |
| `pipeline/protein_smoke.py` | ordering, checks, metrics, failure attribution | silently fill gaps |
| `reporting/*` | manifest/report rendering and hashing | invent measurements |

## Key decisions

**Current-schema-first parsing.** `pick()` takes an ordered list of field names.
The first name is the field the live API returns today; later names are tolerated
fallbacks. Whenever a fallback is used, the fact is recorded in
`schema_notes`/`warnings`, so a schema migration is visible in the report instead
of being masked by compatibility code.

**Residue correspondence through a reference frame.** The experimental trace and
the predicted trace are each aligned to their own sequence, and those sequences
are aligned to UniProt. Composing the maps gives a residue-level correspondence
that tolerates unresolved residues on either side without assuming that residue
numbering is interchangeable between PDB and AlphaFold files.

**Deterministic retries.** Backoff is `base * 2^(n-1)` capped at `max_delay`,
without jitter. `Retry-After` is honoured up to the ceiling. Determinism keeps
failure-injection tests exact and makes logs comparable across runs.

**Manifest/report split.** The manifest is the audit trail (sources + digests +
environment); the report is the scientific result (metrics, checks,
`not_evaluated`). `report.json` embeds the manifest digest, so a report can be
tied to the exact manifest that produced it.

**Failure cannot claim success.** `ensure_consistent()` refuses to write a
manifest that says `SUCCESS` without metrics or with an error attached, and
refuses `FAILED` without a structured error.

## Extension points

* `PipelineStage` and `ErrorCode` are stable enums: new stages (e.g.
  `literature_retrieval`) extend the same vocabulary.
* `Settings` is frozen and serialised into every manifest, so new thresholds are
  reproducible by construction.
* Validation functions are pure (`numpy` in, model out) and independently
  testable; a future evaluation agent can reuse them without the CLI.
* The report schema is versioned (`schema_version`), allowing multi-agent
  orchestration on top without breaking consumers.
