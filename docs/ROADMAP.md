# Roadmap

## Now — v0.1 (implemented)

One real workflow: retrieve → verify → evaluate → report for a
(UniProt accession, PDB entry, chain) triple, with typed failures and hashed
provenance.

## Next — v0.2 (in progress / planned)

| Item | Status | Notes |
| --- | --- | --- |
| All-atom lDDT | planned | side-chain aware scoring; interface will reuse `evaluate_structure` |
| TM-score | planned | only via a validated external aligner with version pinning; not approximated in-house |
| Multi-model comparison | planned | compare AlphaFold model versions for the same entry |
| Batch runs | planned | one manifest per input, shared settings snapshot |
| Optional caching layer | planned | content-addressed by SHA-256, offline replay from cache |

## Later — the research factory (planned, not implemented)

The platform thesis is that autonomous research agents need the same contract
this v0.1 enforces for one workflow: verified inputs, typed failures, hashed
provenance, and explicit `not_evaluated` lists.

Planned multi-agent surface:

* **Literature mining agent** — retrieve and cite primary sources for a claim.
* **Experiment design agent** — propose a reproducible experiment plan with
  pre-registered success/failure criteria.
* **Code generation agent** — implement the plan behind the same validation gates.
* **Evaluation agent** — run the plan, compare against the reference, refuse to
  report unverified numbers.
* **Reporting agent** — assemble a manuscript-grade artifact with provenance.

None of these exist in this repository today. They are listed so the direction is
explicit, not to imply maturity.
