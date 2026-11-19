"""Human-readable Markdown report.

Every number rendered here comes from the manifest. When a metric was not
computed the report says so explicitly; it never fills a gap with a plausible
default.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from nexusforge.models import ArtifactRecord, RunManifest, RunStatus
from nexusforge.reporting.manifest import reproduction_command, write_bytes


def _fmt(value: Any, digits: int = 3, suffix: str = "") -> str:
    if value is None:
        return "not available"
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, (int,)):
        return f"{value}{suffix}"
    if isinstance(value, float):
        return f"{value:.{digits}f}{suffix}"
    return f"{value}{suffix}"


def _table(headers: list[str], rows: list[list[str]]) -> str:
    lines = ["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |"]
    for row in rows:
        lines.append("| " + " | ".join(row) + " |")
    return "\n".join(lines)


def render(manifest: RunManifest) -> str:
    """Render the Markdown report for a finished (successful or failed) run."""
    status = manifest.status
    lines: list[str] = []
    lines.append(
        f"# Protein Structure Confidence Report — {manifest.request.uniprot_accession} / "
        f"{manifest.request.pdb_id} chain {manifest.request.chain_id}"
    )
    lines.append("")
    verdict = (
        "**Status: SUCCESS.** The pipeline retrieved real data, mapped the chains and computed "
        "the metrics listed below."
        if status is RunStatus.SUCCESS
        else "**Status: FAILED.** No scientific conclusion is drawn from this run."
    )
    lines.append(verdict)
    lines.append("")
    lines.append(
        _table(
            ["Field", "Value"],
            [
                ["Run id", f"`{manifest.run_id}`"],
                ["Generated at (UTC)", manifest.generated_at.isoformat()],
                ["UniProt accession", manifest.request.uniprot_accession],
                ["PDB entry", manifest.request.pdb_id],
                ["Author chain", manifest.request.chain_id],
                ["NexusForge version", manifest.software.version],
                ["Reproducible", _fmt(manifest.reproducible)],
            ],
        )
    )
    lines.append("")

    if status is RunStatus.FAILED and manifest.error:
        error = manifest.error
        lines.append("## Failure")
        lines.append("")
        lines.append(
            _table(
                ["Field", "Value"],
                [
                    ["Error code", f"`{error.get('code')}`"],
                    ["Stage", f"`{error.get('stage')}`"],
                    ["Source", str(error.get("source", "—"))],
                    ["Retryable", _fmt(bool(error.get("retryable")))],
                    ["Message", str(error.get("message", "")).replace("|", "\\|")],
                ],
            )
        )
        details = error.get("details")
        if details:
            lines.append("")
            lines.append("Structured details:")
            lines.append("")
            lines.append("```json")
            import json

            lines.append(json.dumps(details, indent=2, ensure_ascii=False))
            lines.append("```")
        lines.append("")
        lines.append(
            "The manifest, sources and any downloaded artefacts are still written to the run "
            "directory so the failure can be reproduced and diagnosed."
        )
        lines.append("")

    # -- data provenance ---------------------------------------------------
    if manifest.data:
        lines.append("## What data was used")
        lines.append("")
        data = manifest.data
        rows: list[list[str]] = []
        if "uniprot" in data:
            uniprot = data["uniprot"]
            rows.append(
                [
                    "UniProt entry",
                    f"{uniprot.get('accession')} — {uniprot.get('protein_name') or 'unnamed'} "
                    f"({uniprot.get('organism')}), reviewed: {_fmt(uniprot.get('reviewed'))}",
                ]
            )
        if "alphafold" in data:
            af = data["alphafold"]
            rows.append(
                [
                    "AlphaFold model",
                    f"{af.get('entry_id')} (v{af.get('latest_version')}), created "
                    f"{af.get('model_created_date')}, pipeline "
                    f"{af.get('tool_used') or 'not reported'}",
                ]
            )
            rows.append(
                [
                    "AlphaFold global pLDDT",
                    _fmt(af.get("global_plddt"), 2) + " (model-level, reported by AlphaFold DB)",
                ]
            )
            if af.get("pae_available") is not None:
                rows.append(
                    ["PAE data", "available" if af["pae_available"] else "not provided by source"]
                )
        if "experimental" in data:
            exp = data["experimental"]
            rows.append(
                [
                    "Experimental entry",
                    f"{exp.get('pdb_id')} — {exp.get('title')}",
                ]
            )
            rows.append(
                [
                    "Experimental method",
                    ", ".join(exp.get("experimental_methods") or []) or "not reported",
                ]
            )
            rows.append(
                ["Resolution", _fmt(exp.get("resolution_angstrom"), 2) + " Å"]
                if exp.get("resolution_angstrom")
                else ["Resolution", "not reported (non-diffraction entry)"]
            )
            rows.append(
                [
                    "Chain → entity",
                    f"entity {exp.get('entity_id')} ({exp.get('entity_name') or 'unnamed'}), "
                    f"UniProt {', '.join(exp.get('uniprot_accessions') or []) or 'none'}",
                ]
            )
        if rows:
            lines.append(_table(["Item", "Value"], rows))
            lines.append("")

    if manifest.sources:
        lines.append("## Where it came from")
        lines.append("")
        rows = []
        for provenance in manifest.sources:
            rows.append(
                [
                    provenance.source,
                    f"[{provenance.url}]({provenance.url})",
                    provenance.fetched_at.isoformat(timespec="seconds"),
                    f"`{provenance.sha256[:16]}…`",
                    str(provenance.size_bytes),
                ]
            )
        lines.append(_table(["Source", "URL", "Fetched (UTC)", "SHA-256", "Bytes"], rows))
        lines.append("")

    # -- sequence validation ----------------------------------------------
    if manifest.alignment is not None:
        lines.append("## Sequence validation")
        lines.append("")
        rows = []
        for label, alignment in (
            ("UniProt → AlphaFold", manifest.alignment.uniprot_vs_alphafold),
            ("UniProt → PDB entity", manifest.alignment.uniprot_vs_pdb_entity),
        ):
            rows.append(
                [
                    label,
                    f"{alignment.query_length} / {alignment.target_length}",
                    f"{alignment.identity:.2%}",
                    f"{alignment.coverage_query:.2%}",
                    str(alignment.gap_openings),
                    str(alignment.mismatches),
                ]
            )
        lines.append(
            _table(["Comparison", "Lengths", "Identity", "Coverage", "Gaps", "Mismatches"], rows)
        )
        lines.append("")
        lines.append(
            "Identity is computed over aligned columns only; coverage is the fraction of the "
            "query sequence that participates in the alignment. A missing N-terminal methionine "
            "therefore appears as one gap and reduced coverage, not as a residue mismatch."
        )
        lines.append("")

    # -- structural metrics ------------------------------------------------
    structure = manifest.metrics.get("structure")
    if structure:
        lines.append("## Structure evaluation")
        lines.append("")
        lines.append(
            _table(
                ["Metric", "Value", "Meaning"],
                [
                    [
                        "Mapped Cα residues",
                        str(structure.get("mapped_residue_count")),
                        "residues resolved in both structures",
                    ],
                    [
                        "RMSD",
                        _fmt(structure.get("rmsd_angstrom"), 3) + " Å",
                        "Cα RMSD after optimal rigid-body superposition",
                    ],
                    [
                        "lDDT-Cα",
                        _fmt(structure.get("lddt_ca"), 4),
                        "Cα-only local distance difference test vs the experimental structure",
                    ],
                    [
                        "lDDT-Cα pairs",
                        str(structure.get("lddt_ca_pairs")),
                        f"pairs within {structure.get('lddt_ca_inclusion_radius')} Å in the "
                        "experimental structure",
                    ],
                    [
                        "lDDT-Cα residues scored",
                        f"{structure.get('lddt_ca_residues_scored')} "
                        f"({structure.get('lddt_ca_residues_unscored')} with no neighbour)",
                        "residues contributing to the mean",
                    ],
                ],
            )
        )
        lines.append("")
        per_threshold = structure.get("lddt_ca_per_threshold") or {}
        if per_threshold:
            lines.append("### lDDT-Cα by distance threshold")
            lines.append("")
            lines.append(
                _table(
                    ["Threshold (Å)", "Conserved fraction"],
                    [[key, f"{value:.4f}"] for key, value in per_threshold.items()],
                )
            )
            lines.append("")
        lines.append(
            "**Scope:** lDDT-Cα is restricted to Cα atoms. It is not a full-atom lDDT and is "
            "not comparable to a full-atom lDDT value reported elsewhere."
        )
        lines.append("")

    # -- confidence --------------------------------------------------------
    confidence = manifest.metrics.get("confidence")
    if confidence:
        lines.append("## Confidence evaluation")
        lines.append("")
        lines.append(
            _table(
                ["Statistic", "pLDDT"],
                [
                    ["Residues", str(confidence.get("residue_count"))],
                    ["Mean", _fmt(confidence.get("mean_plddt"), 2)],
                    ["Median", _fmt(confidence.get("median_plddt"), 2)],
                    ["Min", _fmt(confidence.get("min_plddt"), 2)],
                    ["Max", _fmt(confidence.get("max_plddt"), 2)],
                ],
            )
        )
        lines.append("")
        bins = confidence.get("bins") or []
        if bins:
            lines.append("### pLDDT bins vs measured lDDT-Cα")
            lines.append("")
            rows = []
            for entry in bins:
                rows.append(
                    [
                        entry.get("label", ""),
                        f"{entry.get('lower')}–{min(entry.get('upper', 100), 100)}",
                        str(entry.get("residue_count")),
                        f"{entry.get('fraction', 0):.1%}",
                        _fmt(entry.get("mean_lddt_ca"), 4),
                    ]
                )
            lines.append(_table(["pLDDT bin", "Range", "Residues", "Share", "Mean lDDT-Cα"], rows))
            lines.append("")
        lines.append(
            _table(
                ["Relationship", "Value"],
                [
                    [
                        "Residues with both pLDDT and lDDT-Cα",
                        str(confidence.get("correlated_residue_count", 0)),
                    ],
                    ["Pearson r (pLDDT vs lDDT-Cα)", _fmt(confidence.get("pearson_plddt_lddt"))],
                    ["Spearman ρ (pLDDT vs lDDT-Cα)", _fmt(confidence.get("spearman_plddt_lddt"))],
                    [
                        "pLDDT cross-check (confidence JSON vs model B-factors)",
                        _fmt(confidence.get("cross_check_max_abs_delta"), 4),
                    ],
                ],
            )
        )
        lines.append("")
        lines.append(
            "**pLDDT is a prediction-confidence estimate, not ground truth.** The lDDT-Cα column "
            "is measured against one experimental structure; disagreements between the two are "
            "informative, not errors in either quantity."
        )
        lines.append("")

    pae = manifest.metrics.get("pae")
    if pae:
        lines.append("## Predicted aligned error (PAE)")
        lines.append("")
        lines.append(
            _table(
                ["Statistic", "Value"],
                [
                    ["Matrix shape", str(tuple(pae.get("shape", [])))],
                    ["Mean PAE", _fmt(pae.get("mean_predicted_aligned_error"), 3) + " Å"],
                    ["90th percentile PAE", _fmt(pae.get("p90_predicted_aligned_error"), 3) + " Å"],
                    ["Source-declared maximum", _fmt(pae.get("max_predicted_aligned_error"), 3)],
                ],
            )
        )
        lines.append("")

    # -- checks ------------------------------------------------------------
    if manifest.checks:
        lines.append("## Validation checks")
        lines.append("")
        rows = []
        for check in manifest.checks:
            rows.append(
                [
                    f"`{check.check_id}`",
                    check.title,
                    check.status.value.upper(),
                    (check.detail or "").replace("|", "\\|"),
                ]
            )
        lines.append(_table(["Check", "What it verifies", "Status", "Detail"], rows))
        lines.append("")

    # -- not evaluated -----------------------------------------------------
    not_evaluated = manifest.metrics.get("not_evaluated") or []
    lines.append("## Not evaluated in this run")
    lines.append("")
    if not_evaluated:
        for item in not_evaluated:
            lines.append(f"- **{item.get('metric')}** — {item.get('reason')}")
    else:
        lines.append("- none")
    lines.append("")

    limitations = manifest.metrics.get("limitations") or []
    if limitations:
        lines.append("## Limitations")
        lines.append("")
        for item in limitations:
            lines.append(f"- {item}")
        lines.append("")

    if manifest.warnings:
        lines.append("## Warnings from this run")
        lines.append("")
        for warning in manifest.warnings:
            lines.append(f"- {warning}")
        lines.append("")

    # -- reproduction ------------------------------------------------------
    lines.append("## Reproduction")
    lines.append("")
    lines.append("```bash")
    lines.append(reproduction_command(manifest))
    lines.append("```")
    lines.append("")
    lines.append(
        "Every request URL, HTTP status, fetch timestamp and SHA-256 digest is recorded in "
        "`manifest.json`; raw source responses are kept under `raw/`. Re-running the command "
        "re-fetches the sources, so byte-identical inputs require comparing those digests rather "
        "than assuming the upstream files are immutable."
    )
    lines.append("")
    return "\n".join(lines)


def write_markdown(run_dir: Path, manifest: RunManifest) -> ArtifactRecord:
    content = render(manifest).encode("utf-8")
    return write_bytes(run_dir / "report.md", content)


__all__ = ["render", "write_markdown"]
