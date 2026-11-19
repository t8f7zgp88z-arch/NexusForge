"""Command-line interface.

Exit codes are part of the contract:

``0`` success · ``2`` bad input · ``3`` data/validation failure · ``4`` network
or source failure · ``5`` internal/reporting failure.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import typer

from nexusforge import __version__
from nexusforge.clients.alphafold import BASE_URL as ALPHAFOLD_BASE
from nexusforge.clients.rcsb import DATA_BASE_URL, FILES_BASE_URL
from nexusforge.clients.uniprot import BASE_URL as UNIPROT_BASE
from nexusforge.config import load_settings
from nexusforge.errors import EXIT_OK, NexusForgeError
from nexusforge.pipeline.protein_smoke import ProteinSmokePipeline, SmokeRequest
from nexusforge.reporting.manifest import dependency_versions

app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help=(
        "NexusForge — evidence-first research engineering.\n\n"
        "v0.1 product: Protein Structure Confidence Explorer. Retrieves a UniProt "
        "entry, its AlphaFold DB prediction and an experimental PDB chain, verifies "
        "that the sequences genuinely correspond, then reports Cα RMSD, lDDT-Cα and "
        "pLDDT statistics with full provenance."
    ),
)


def _print_version(value: bool) -> None:
    if value:
        typer.echo(f"nexusforge {__version__}")
        raise typer.Exit(EXIT_OK)


@app.callback()
def main(
    version: bool = typer.Option(
        False,
        "--version",
        "-V",
        help="Show the NexusForge version and exit.",
        callback=_print_version,
        is_eager=True,
    ),
) -> None:
    """NexusForge command line."""


@app.command("version")
def version_command() -> None:
    """Print the version and the versions of tracked dependencies."""
    typer.echo(f"nexusforge {__version__}")
    typer.echo(f"python {sys.version.split()[0]}")
    for name, value in dependency_versions().items():
        typer.echo(f"{name} {value}")


@app.command("protein-smoke")
def protein_smoke(
    uniprot: str = typer.Option(..., "--uniprot", help="UniProt accession, e.g. P69905"),
    pdb: str = typer.Option(..., "--pdb", help="PDB entry id, e.g. 4HHB"),
    chain: str = typer.Option(..., "--chain", help="Author chain id, e.g. A"),
    output: Path = typer.Option(
        Path("runs"),
        "--output",
        help="Run directory to create (manifest, reports, raw sources).",
    ),
    force: bool = typer.Option(False, "--force", help="Write into a non-empty output directory."),
    timeout: float | None = typer.Option(
        None, "--timeout", help="Per-request HTTP timeout in seconds."
    ),
    no_pae: bool = typer.Option(False, "--no-pae", help="Skip PAE retrieval."),
    quiet: bool = typer.Option(False, "--quiet", help="Only print the run directory and status."),
) -> None:
    """Run the protein structure confidence pipeline end to end."""
    settings = load_settings(
        timeout_seconds=timeout,
        fetch_pae=False if no_pae else None,
    )
    pipeline = ProteinSmokePipeline(settings=settings)
    try:
        outcome = pipeline.run(
            SmokeRequest(
                uniprot_accession=uniprot,
                pdb_id=pdb,
                chain_id=chain,
                output_dir=output,
                overwrite=force,
            )
        )
    except NexusForgeError as exc:
        typer.secho(
            f"error [{exc.stage.value}] {exc.code.value}: {exc.message}", err=True, fg="red"
        )
        raise typer.Exit(exc.exit_code) from exc
    finally:
        pipeline.close()

    if outcome.error is not None:
        error = outcome.error
        typer.secho(
            f"run failed [{error.stage.value}] {error.code.value}: {error.message}",
            err=True,
            fg="red",
        )
        typer.echo(f"artifacts: {outcome.run_dir}")
        raise typer.Exit(error.exit_code)

    manifest = outcome.manifest
    metrics = manifest.metrics
    if not quiet:
        typer.echo(f"status:     {manifest.status.value}")
        typer.echo(f"run id:     {manifest.run_id}")
        typer.echo(f"directory:  {outcome.run_dir}")
        if manifest.data.get("uniprot"):
            uniprot_data = manifest.data["uniprot"]
            typer.echo(
                f"uniprot:    {uniprot_data['accession']} — {uniprot_data['protein_name']} "
                f"({uniprot_data['organism']})"
            )
        if manifest.data.get("alphafold"):
            alphafold = manifest.data["alphafold"]
            typer.echo(
                f"alphafold:  {alphafold['entry_id']} (v{alphafold['latest_version']}, "
                f"global pLDDT {alphafold['global_plddt']})"
            )
        if manifest.data.get("experimental"):
            experimental = manifest.data["experimental"]
            typer.echo(
                f"pdb:        {experimental['pdb_id']} — entity {experimental['entity_id']}, "
                f"chain {experimental['chain_id']}"
            )
        sequence = metrics.get("sequence", {})
        structure = metrics.get("structure", {})
        confidence = metrics.get("confidence", {})
        if sequence:
            typer.echo(
                f"identity:   {sequence['pdb_entity_identity']:.2%} over aligned residues, "
                f"coverage {sequence['pdb_entity_coverage_query']:.2%}"
            )
        if structure:
            typer.echo(
                f"rmsd:       {structure['rmsd_angstrom']:.3f} Å over "
                f"{structure['mapped_residue_count']} mapped Cα"
            )
            typer.echo(f"lddt-ca:    {structure['lddt_ca']:.4f}")
        if confidence:
            typer.echo(
                f"plddt:      mean {confidence['mean_plddt']:.2f}, "
                f"median {confidence['median_plddt']:.2f}"
            )
        typer.echo(f"report:     {outcome.run_dir / 'report.md'}")
    else:
        typer.echo(f"{manifest.status.value} {outcome.run_dir}")
    raise typer.Exit(EXIT_OK)


@app.command("doctor")
def doctor(
    timeout: float = typer.Option(10.0, "--timeout", help="Per-request timeout in seconds."),
) -> None:
    """Check that the three public data sources are reachable right now."""
    import httpx

    settings = load_settings()
    targets = [
        ("uniprot", f"{UNIPROT_BASE}/uniprotkb/P69905.json"),
        ("alphafold", f"{ALPHAFOLD_BASE}/api/prediction/P69905"),
        ("rcsb-metadata", f"{DATA_BASE_URL}/entry/4HHB"),
        ("rcsb-files", f"{FILES_BASE_URL}/4HHB.cif"),
    ]
    failures = 0
    for name, url in targets:
        started = time.monotonic()
        try:
            response = httpx.get(
                url,
                timeout=timeout,
                follow_redirects=True,
                headers={"User-Agent": settings.effective_user_agent},
            )
            elapsed_ms = (time.monotonic() - started) * 1000
            ok = response.status_code < 400
            status = f"HTTP {response.status_code}"
        except httpx.HTTPError as exc:
            elapsed_ms = (time.monotonic() - started) * 1000
            ok = False
            status = f"{type(exc).__name__}"
        marker = "ok  " if ok else "FAIL"
        typer.echo(f"{marker} {name:14s} {status:16s} {elapsed_ms:7.0f} ms  {url}")
        if not ok:
            failures += 1
    if failures:
        typer.secho(f"{failures} target(s) unreachable", err=True, fg="red")
        raise typer.Exit(4)
    typer.echo("all data sources reachable")


if __name__ == "__main__":  # pragma: no cover
    app()
