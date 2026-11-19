"""Live network smoke test — opt-in, never part of CI.

Run explicitly:

    NEXUSFORGE_LIVE=1 pytest -q -m network tests/integration/test_live_smoke.py

It hits UniProt, AlphaFold DB and RCSB for real. Third-party outages must not
fail the default test suite, which is why this module is skipped unless asked
for.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from nexusforge.config import load_settings
from nexusforge.models import RunStatus
from nexusforge.pipeline.protein_smoke import ProteinSmokePipeline, SmokeRequest

pytestmark = [
    pytest.mark.network,
    pytest.mark.skipif(
        os.environ.get("NEXUSFORGE_LIVE") != "1",
        reason="live API test disabled; set NEXUSFORGE_LIVE=1 to enable",
    ),
]


def test_live_p69905_against_4hhb_chain_a(tmp_path: Path) -> None:
    settings = load_settings(env={}, timeout_seconds=60.0, max_attempts=3)
    pipeline = ProteinSmokePipeline(settings=settings)
    try:
        outcome = pipeline.run(SmokeRequest("P69905", "4HHB", "A", tmp_path / "live"))
    finally:
        pipeline.close()

    assert outcome.error is None, outcome.error
    assert outcome.manifest.status is RunStatus.SUCCESS
    metrics = outcome.manifest.metrics
    assert metrics["sequence"]["pdb_entity_identity"] == pytest.approx(1.0)
    assert 0.0 < metrics["structure"]["rmsd_angstrom"] < 2.0
    assert 0.8 <= metrics["structure"]["lddt_ca"] <= 1.0
    assert metrics["confidence"]["residue_count"] == 142
    assert (tmp_path / "live" / "report.json").exists()
