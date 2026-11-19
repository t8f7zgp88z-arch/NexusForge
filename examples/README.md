# Example run

`sample-run/` is a trimmed copy of a real execution:

```bash
nexusforge protein-smoke --uniprot P69905 --pdb 4HHB --chain A --output runs/P69905-4HHB-A
```

Executed 2026-11-19 against the live UniProt / AlphaFold DB / RCSB APIs.

| File | Contents |
| --- | --- |
| `report.md` | Human-readable report for P69905 vs 4HHB chain A |
| `metrics.json` | Compact metric extract (sequence, structure, confidence) |
| `manifest.json` | Full audit trail: sources, SHA-256 digests, software versions, checks |

Not copied here (git-ignored because of size): `raw/` with the byte-identical
API responses (mmCIF structures ≈ 0.9 MB) and `logs/run.log`. Reproduce them by
running the command above; the digests in `manifest.json` let you verify that the
upstream files are unchanged.

## Headline result

| Metric | Value |
| --- | --- |
| Sequence identity (UniProt → PDB entity) | 100.00 % over 141 aligned residues |
| Coverage of UniProt sequence | 99.30 % (initiator Met not resolved in the crystal) |
| Mapped Cα residues | 141 |
| Cα RMSD | 0.335 Å |
| lDDT-Cα | 0.9921 |
| pLDDT | mean 98.06, min 65.38, max 98.88 |

Interpretation: the AlphaFold v6 model of human α-globin reproduces the
experimental Cα geometry of 4HHB chain A closely at 1.74 Å resolution. lDDT-Cα is
Cα-only; do not compare it with full-atom lDDT values.
