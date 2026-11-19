# Validation and failure semantics

## The two confidence concepts

| Quantity | What it is | What it is not |
| --- | --- | --- |
| **pLDDT** | AlphaFold's per-residue confidence estimate, 0–100, stored in the model's B-factor column and in the confidence JSON | not an experimental measurement, not lDDT, not ground truth |
| **lDDT-Cα** | fraction of Cα–Cα distances conserved relative to an experimental structure, at 0.5/1/2/4 Å thresholds, within a 15 Å inclusion radius | not full-atom lDDT; not computed by AlphaFold |

The report never converts one into the other. It shows pLDDT bins alongside the
measured mean lDDT-Cα of the residues in each bin, plus Pearson/Spearman
correlation over residues where both exist. A low correlation is a result, not an
error.

## Sequence validation gates

| Gate | Default | Failure |
| --- | --- | --- |
| Valid amino-acid alphabet | required | `INVALID_SEQUENCE` |
| Identity vs UniProt (AF model, PDB entity) | ≥ 90 % | `SEQUENCE_MISMATCH` |
| UniProt coverage | ≥ 50 % | `SEQUENCE_MISMATCH` |
| Observability of the trace against its sequence | ≥ 98 % identity | `SEQUENCE_MISMATCH` |
| Chain → polymer entity resolution | exactly one | `CHAIN_MAPPING_AMBIGUOUS` |
| Entity → UniProt cross-reference | must contain the accession | `SEQUENCE_MISMATCH` |
| Shared residues for comparison | ≥ 3 | `NO_BASELINE` |

Thresholds are configurable (`Settings`) and are written into every manifest, so
a result can be read together with the gates that produced it.

## Structural metrics

* **RMSD** — Cα RMSD after optimal rigid-body superposition (Kabsch, with
  reflection correction). Sensitive to a single large displacement; reported in Å.
* **lDDT-Cα** — superposition-free local distance test. Per residue: mean over
  four distance thresholds of the conserved-pair fraction. Global: mean over
  residues that have at least one neighbour within the inclusion radius.
  Residues without neighbours are counted as `unscored`, not as zero.
* **Pairs** — the number of Cα pairs inside the 15 Å inclusion radius in the
  experimental structure; reported so the score's support is visible.

## Failure semantics

1. Every failure raises `NexusForgeError` with a stable `ErrorCode` and a
   `PipelineStage`.
2. The pipeline catches it, sets `status = FAILED`, records the error in the
   manifest, and still writes `manifest.json`, `report.json`, `report.md`,
   `metrics.json`, the run log and any raw data fetched before the failure.
3. `ensure_consistent()` rejects a `SUCCESS` manifest that lacks metrics or
   carries an error, and a `FAILED` manifest without an error — so the contract
   cannot be violated by a future code change silently.
4. The CLI exit code maps to the error family (2 input, 3 data/validation,
   4 network/source, 5 internal/reporting).

### Injected failures covered by tests

| Injection | Expected outcome |
| --- | --- |
| HTTP 404 on UniProt | `SOURCE_NOT_FOUND`, stage `data_retrieval`, exit 3 |
| HTTP 500 always | `SOURCE_UNAVAILABLE` after bounded retries |
| HTTP 500 twice then 200 | success, 2 sleeps, attempt count 3 |
| HTTP 429 with `Retry-After` | retry honours header, capped by ceiling |
| Read timeout | `NETWORK_ERROR`, retryable, bounded attempts |
| Empty 200 body | `INVALID_API_RESPONSE` ("empty body") |
| HTML instead of JSON | `INVALID_API_RESPONSE` with body preview |
| Chain B of 4HHB requested as P69905 | `SEQUENCE_MISMATCH` (P68871 ≠ P69905) |
| Unknown chain Z | `SEQUENCE_MISMATCH` listing polymer chains |
| Truncated pLDDT vector | `CONFIDENCE_LENGTH_MISMATCH`, stage `confidence_evaluation` |
| pLDDT value 140 | `INVALID_API_RESPONSE` (outside `[0, 100]`) |
| Missing PAE document | run succeeds, PAE omitted, listed under `not_evaluated` |

No injected failure is ever reported as `SUCCESS`; each assertion is encoded in
`tests/integration/test_pipeline_failures.py` and `tests/unit/test_http.py`.
