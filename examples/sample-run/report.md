# Protein Structure Confidence Report — P69905 / 4HHB chain A

**Status: SUCCESS.** The pipeline retrieved real data, mapped the chains and computed the metrics listed below.

| Field | Value |
| --- | --- |
| Run id | `P69905-4HHB-A` |
| Generated at (UTC) | 2026-11-19T15:21:13.913594+00:00 |
| UniProt accession | P69905 |
| PDB entry | 4HHB |
| Author chain | A |
| NexusForge version | 0.1.0 |
| Reproducible | yes |

## What data was used

| Item | Value |
| --- | --- |
| UniProt entry | P69905 — Hemoglobin subunit alpha (Homo sapiens), reviewed: yes |
| AlphaFold model | AF-P69905-F1 (v6), created 2025-08-01T00:00:00+00:00, pipeline AlphaFold Monomer v2.0 pipeline |
| AlphaFold global pLDDT | 98.06 (model-level, reported by AlphaFold DB) |
| PAE data | available |
| Experimental entry | 4HHB — THE CRYSTAL STRUCTURE OF HUMAN DEOXYHAEMOGLOBIN AT 1.74 ANGSTROMS RESOLUTION |
| Experimental method | X-RAY DIFFRACTION |
| Resolution | 1.74 Å |
| Chain → entity | entity 1 (Hemoglobin subunit alpha), UniProt P69905 |

## Where it came from

| Source | URL | Fetched (UTC) | SHA-256 | Bytes |
| --- | --- | --- | --- | --- |
| uniprot | [https://rest.uniprot.org/uniprotkb/P69905.json](https://rest.uniprot.org/uniprotkb/P69905.json) | 2026-11-19T15:21:15+00:00 | `78a781ffcbc7ee16…` | 241776 |
| alphafold | [https://alphafold.ebi.ac.uk/api/prediction/P69905](https://alphafold.ebi.ac.uk/api/prediction/P69905) | 2026-11-19T15:21:15+00:00 | `af907bb73c411790…` | 1987 |
| alphafold-file | [https://alphafold.ebi.ac.uk/files/AF-P69905-F1-model_v6.cif](https://alphafold.ebi.ac.uk/files/AF-P69905-F1-model_v6.cif) | 2026-11-19T15:21:16+00:00 | `eb7adbdd79775b0d…` | 137559 |
| alphafold-confidence | [https://alphafold.ebi.ac.uk/files/AF-P69905-F1-confidence_v6.json](https://alphafold.ebi.ac.uk/files/AF-P69905-F1-confidence_v6.json) | 2026-11-19T15:21:17+00:00 | `523095e99862ad9f…` | 1930 |
| alphafold-pae | [https://alphafold.ebi.ac.uk/files/AF-P69905-F1-predicted_aligned_error_v6.json](https://alphafold.ebi.ac.uk/files/AF-P69905-F1-predicted_aligned_error_v6.json) | 2026-11-19T15:21:17+00:00 | `9d6bfee2b49f0391…` | 40977 |
| rcsb-entry | [https://data.rcsb.org/rest/v1/core/entry/4HHB](https://data.rcsb.org/rest/v1/core/entry/4HHB) | 2026-11-19T15:21:18+00:00 | `a95a79421b902116…` | 26055 |
| rcsb-entity | [https://data.rcsb.org/rest/v1/core/polymer_entity/4HHB/1](https://data.rcsb.org/rest/v1/core/polymer_entity/4HHB/1) | 2026-11-19T15:21:18+00:00 | `f7585e98ac51bca4…` | 53683 |
| rcsb-entity | [https://data.rcsb.org/rest/v1/core/polymer_entity/4HHB/2](https://data.rcsb.org/rest/v1/core/polymer_entity/4HHB/2) | 2026-11-19T15:21:18+00:00 | `b0fda3e54a7b4fbe…` | 49486 |
| rcsb-file | [https://files.rcsb.org/download/4HHB.cif](https://files.rcsb.org/download/4HHB.cif) | 2026-11-19T15:21:19+00:00 | `2d3ca1bf21ddff35…` | 772198 |

## Sequence validation

| Comparison | Lengths | Identity | Coverage | Gaps | Mismatches |
| --- | --- | --- | --- | --- | --- |
| UniProt → AlphaFold | 142 / 142 | 100.00% | 100.00% | 0 | 0 |
| UniProt → PDB entity | 142 / 141 | 100.00% | 99.30% | 1 | 0 |

Identity is computed over aligned columns only; coverage is the fraction of the query sequence that participates in the alignment. A missing N-terminal methionine therefore appears as one gap and reduced coverage, not as a residue mismatch.

## Structure evaluation

| Metric | Value | Meaning |
| --- | --- | --- |
| Mapped Cα residues | 141 | residues resolved in both structures |
| RMSD | 0.335 Å | Cα RMSD after optimal rigid-body superposition |
| lDDT-Cα | 0.9921 | Cα-only local distance difference test vs the experimental structure |
| lDDT-Cα pairs | 3078 | pairs within 15.0 Å in the experimental structure |
| lDDT-Cα residues scored | 141 (0 with no neighbour) | residues contributing to the mean |

### lDDT-Cα by distance threshold

| Threshold (Å) | Conserved fraction |
| --- | --- |
| 0.5 | 0.9685 |
| 1 | 1.0000 |
| 2 | 1.0000 |
| 4 | 1.0000 |

**Scope:** lDDT-Cα is restricted to Cα atoms. It is not a full-atom lDDT and is not comparable to a full-atom lDDT value reported elsewhere.

## Confidence evaluation

| Statistic | pLDDT |
| --- | --- |
| Residues | 142 |
| Mean | 98.06 |
| Median | 98.56 |
| Min | 65.38 |
| Max | 98.88 |

### pLDDT bins vs measured lDDT-Cα

| pLDDT bin | Range | Residues | Share | Mean lDDT-Cα |
| --- | --- | --- | --- | --- |
| Very low | 0.0–50.0 | 0 | 0.0% | not available |
| Low | 50.0–70.0 | 1 | 0.7% | not available |
| Confident | 70.0–90.0 | 0 | 0.0% | not available |
| Very high | 90.0–100.0 | 141 | 99.3% | 0.9921 |

| Relationship | Value |
| --- | --- |
| Residues with both pLDDT and lDDT-Cα | 141 |
| Pearson r (pLDDT vs lDDT-Cα) | 0.030 |
| Spearman ρ (pLDDT vs lDDT-Cα) | 0.192 |
| pLDDT cross-check (confidence JSON vs model B-factors) | 0.0000 |

**pLDDT is a prediction-confidence estimate, not ground truth.** The lDDT-Cα column is measured against one experimental structure; disagreements between the two are informative, not errors in either quantity.

## Predicted aligned error (PAE)

| Statistic | Value |
| --- | --- |
| Matrix shape | (142, 142) |
| Mean PAE | 2.062 Å |
| 90th percentile PAE | 3.000 Å |
| Source-declared maximum | 31.750 |

## Validation checks

| Check | What it verifies | Status | Detail |
| --- | --- | --- | --- |
| `uniprot_sequence_valid` | UniProt sequence uses valid amino-acid characters | PASSED | 142 residues, 0 ambiguous code(s) |
| `alphafold_sequence_matches_uniprot` | AlphaFold model sequence corresponds to the UniProt entry | PASSED | identity 100.00%, coverage 100.00% |
| `plddt_count_matches_sequence` | pLDDT vector length equals the modelled sequence length | PASSED | 142 values for 142 modelled residues |
| `pdb_chain_maps_to_uniprot` | Experimental chain is cross-referenced to the UniProt accession | PASSED | entity 1 of 4HHB lists P69905 |
| `pdb_entity_sequence_matches_uniprot` | Experimental entity sequence corresponds to the UniProt entry | PASSED | identity 100.00%, coverage 99.30%, 1 gap(s), 1 residue position(s) absent from one side |
| `coordinates_parsed` | Both structures parsed into Cα traces | PASSED | experimental 141 Cα, prediction 142 Cα |
| `lddt_ca_defined` | lDDT-Cα has defined residue pairs and scored residues | PASSED | 3078 pairs within 15.0 Å, 141 residues scored |
| `rmsd_computed` | Cα RMSD computed after optimal superposition | PASSED | 0.335 Å over 141 residues |
| `plddt_matches_model_bfactors` | Confidence JSON agrees with the pLDDT stored in the model file | PASSED | largest \|Δ\| = 0.0000 |
| `plddt_range_valid` | pLDDT values lie in [0, 100] and cover the model | PASSED | min 65.38, max 98.88, 142 residues |
| `pae_statistics` | Predicted aligned error statistics extracted | PASSED | mean 2.062 Å |
| `tm_score` | TM-score comparison | SKIPPED | not implemented in v0.1: no validated structure-alignment program is bundled |

## Not evaluated in this run

- **full-atom lDDT** — v0.1 computes the Cα-only variant (lDDT-Cα); all-atom scoring is planned
- **TM-score** — no validated structure-alignment program is bundled in v0.1
- **AlphaFold inference** — predictions are read from AlphaFold DB; no model is run by this project
- **oligomeric assembly assessment** — only the single requested chain is compared
- **per-domain PAE interpretation** — only aggregate PAE statistics are computed in v0.1

## Limitations

- AlphaFold DB predictions are consumed as published; this project does not run the AlphaFold model itself.
- The comparison uses one experimental chain. Quaternary context (e.g. haemoglobin tetramer assembly) is not modelled or scored.
- Only positions resolved in both structures are compared; unobserved residues are reported, never imputed.
- lDDT-Cα is restricted to Cα atoms and is not comparable with a full-atom lDDT value.
- pLDDT is a model confidence estimate; the lDDT-Cα column is measured against one experimental structure, so agreement or disagreement between them is an observation, not a calibrated error.

## Warnings from this run

- Experimental chain contains non-standard residues: HEM

## Reproduction

```bash
nexusforge protein-smoke --uniprot P69905 --pdb 4HHB --chain A --output runs/P69905-4HHB-A
```

Every request URL, HTTP status, fetch timestamp and SHA-256 digest is recorded in `manifest.json`; raw source responses are kept under `raw/`. Re-running the command re-fetches the sources, so byte-identical inputs require comparing those digests rather than assuming the upstream files are immutable.
