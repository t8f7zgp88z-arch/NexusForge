"""Scientific validation: sequence alignment, geometry, lDDT-Cα, pLDDT handling."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from nexusforge.errors import ConfidenceError, ErrorCode, StructureError, ValidationError
from nexusforge.models import ConfidenceScores
from nexusforge.validation.confidence import (
    cross_check_confidence,
    plddt_from_trace,
    summarise_confidence,
    validate_confidence,
)
from nexusforge.validation.sequence import (
    align_sequences,
    require_compatible,
    validate_sequence,
)
from nexusforge.validation.structure import (
    compute_lddt_ca,
    evaluate_structure,
    extract_ca_trace,
    kabsch_superpose,
    match_residues,
    rmsd,
)

FIXTURES = Path(__file__).parent.parent / "fixtures"


# --------------------------------------------------------------------------- #
# Sequence validation
# --------------------------------------------------------------------------- #
def test_validate_sequence_accepts_real_protein() -> None:
    stats = validate_sequence("MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHG")
    assert stats.length == 60
    assert stats.invalid_characters == []


@pytest.mark.parametrize("bad", ["", "   ", "MVLS-PADK", "MVLS*PADK", "MVLS123"])
def test_validate_sequence_rejects_invalid(bad: str) -> None:
    with pytest.raises(ValidationError) as excinfo:
        validate_sequence(bad)
    assert excinfo.value.code is ErrorCode.INVALID_SEQUENCE


def test_alignment_recovers_identity_and_terminal_gap() -> None:
    query = "MVLSPADKTNVKAAWGKVGA"
    target = "VLSPADKTNVKAAWGKVGA"  # missing initiator methionine
    result = align_sequences(query, target, query_id="uniprot", target_id="pdb")
    assert result.matches == 19
    assert result.mismatches == 0
    assert result.gap_openings == 1
    assert result.identity == 1.0
    assert result.coverage_query == pytest.approx(19 / 20)
    assert result.coverage_target == 1.0
    assert result.mapping[1] == 0  # second query residue maps to first target residue


def test_alignment_detects_mismatch() -> None:
    result = align_sequences("AAAAAAAAAA", "AAAAAKAAAA", query_id="a", target_id="b")
    assert result.mismatches == 1
    assert result.identity == pytest.approx(0.9)


def test_require_compatible_rejects_unrelated_sequences() -> None:
    result = align_sequences("AAAAKKKK", "WWWWWWWW", query_id="a", target_id="b")
    with pytest.raises(ValidationError) as excinfo:
        require_compatible(result, min_identity=0.9, min_coverage=0.5)
    assert excinfo.value.code is ErrorCode.SEQUENCE_MISMATCH


# --------------------------------------------------------------------------- #
# Geometry
# --------------------------------------------------------------------------- #
def _random_ca(n: int = 30, seed: int = 7) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return np.cumsum(rng.normal(scale=2.0, size=(n, 3)), axis=0)


def test_kabsch_recovers_known_rotation() -> None:
    reference = _random_ca()
    angle = np.deg2rad(37.0)
    rotation = np.array(
        [
            [np.cos(angle), -np.sin(angle), 0.0],
            [np.sin(angle), np.cos(angle), 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    translation = np.array([4.0, -2.0, 9.0])
    mobile = (reference - translation) @ rotation.T

    recovered_rotation, recovered_translation = kabsch_superpose(mobile, reference)
    value = rmsd(mobile, reference, recovered_rotation, recovered_translation)
    assert value == pytest.approx(0.0, abs=1e-6)
    assert np.linalg.det(recovered_rotation) == pytest.approx(1.0)


def test_kabsch_rejects_too_few_points() -> None:
    with pytest.raises(StructureError):
        kabsch_superpose(np.zeros((2, 3)), np.zeros((2, 3)))


def test_rmsd_is_zero_for_identical_coordinates() -> None:
    coords = _random_ca(12)
    rotation, translation = kabsch_superpose(coords, coords)
    assert rmsd(coords, coords, rotation, translation) == pytest.approx(0.0, abs=1e-9)


def test_rmsd_grows_with_perturbation() -> None:
    reference = _random_ca(20, seed=3)
    mobile = reference.copy()
    mobile[5] += np.array([3.0, 0.0, 0.0])
    rotation, translation = kabsch_superpose(mobile, reference)
    assert rmsd(mobile, reference, rotation, translation) > 0.1


# --------------------------------------------------------------------------- #
# lDDT-Cα
# --------------------------------------------------------------------------- #
def test_lddt_ca_is_one_for_identical_structures() -> None:
    coords = _random_ca(25)
    result = compute_lddt_ca(coords, coords)
    assert result.score == pytest.approx(1.0)
    assert result.per_threshold == {"0.5": 1.0, "1": 1.0, "2": 1.0, "4": 1.0}
    assert result.pairs > 0
    assert result.unscored_residues == 0


def test_lddt_ca_drops_when_one_residue_moves() -> None:
    coords = _random_ca(25, seed=11)
    moved = coords.copy()
    moved[12] += np.array([6.0, 0.0, 0.0])
    result = compute_lddt_ca(coords, moved)
    assert result.score < 1.0
    assert result.per_residue[12] < 1.0


def test_lddt_ca_excludes_self_pairs() -> None:
    coords = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [6.0, 0.0, 0.0]])
    result = compute_lddt_ca(coords, coords, inclusion_radius=5.0)
    assert result.pairs == 2  # (0,1) and (1,2); (0,2) is 6 Å away


def test_lddt_ca_is_superposition_free() -> None:
    coords = _random_ca(20)
    rotated = coords @ np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]], dtype=float) + 50.0
    assert compute_lddt_ca(coords, rotated).score == pytest.approx(1.0)


def test_lddt_ca_requires_pairs_within_radius() -> None:
    coords = np.array([[0.0, 0.0, 0.0], [100.0, 0.0, 0.0]])
    with pytest.raises(StructureError) as excinfo:
        compute_lddt_ca(coords, coords, inclusion_radius=15.0)
    assert excinfo.value.code is ErrorCode.NO_BASELINE


def test_lddt_ca_rejects_mismatched_arrays() -> None:
    with pytest.raises(StructureError):
        compute_lddt_ca(np.zeros((5, 3)), np.zeros((6, 3)))


# --------------------------------------------------------------------------- #
# Traces from the real fixture structures
# --------------------------------------------------------------------------- #
def _write(tmp_path, name: str) -> Path:
    import gzip

    path = tmp_path / name
    path.write_bytes(gzip.decompress((FIXTURES / f"{name}.gz").read_bytes()))
    return path


def test_extract_ca_trace_from_prediction(tmp_path) -> None:
    path = _write(tmp_path, "AF-P69905-F1-model_v6.cif")
    trace = extract_ca_trace(path, label="alphafold", chain_id="A")
    assert len(trace) == 142
    assert trace.sequence.startswith("MVLS")
    assert trace.residue_numbers[0] == 1
    assert trace.bfactors.min() >= 0.0


def test_extract_ca_trace_from_experimental_structure(tmp_path) -> None:
    path = _write(tmp_path, "4HHB.cif")
    trace = extract_ca_trace(path, label="pdb", chain_id="A")
    assert len(trace) == 141
    assert trace.sequence.startswith("VLSP")  # initiator methionine not resolved
    assert trace.chain_id == "A"


def test_extract_ca_trace_missing_chain(tmp_path) -> None:
    path = _write(tmp_path, "4HHB.cif")
    with pytest.raises(StructureError) as excinfo:
        extract_ca_trace(path, label="pdb", chain_id="Z")
    assert excinfo.value.code is ErrorCode.INVALID_COORDINATES
    assert "not present" in excinfo.value.message


def test_extract_ca_trace_rejects_garbage(tmp_path) -> None:
    path = tmp_path / "broken.cif"
    path.write_text("data_broken\n#\n")
    with pytest.raises(StructureError):
        extract_ca_trace(path, label="broken")


# --------------------------------------------------------------------------- #
# End-to-end structural comparison on the real fixtures
# --------------------------------------------------------------------------- #
def test_real_structure_comparison_produces_expected_ranges(tmp_path) -> None:
    predicted_path = _write(tmp_path, "AF-P69905-F1-model_v6.cif")
    experimental_path = _write(tmp_path, "4HHB.cif")
    experimental = extract_ca_trace(experimental_path, label="pdb", chain_id="A")
    predicted = extract_ca_trace(predicted_path, label="alphafold", chain_id="A")

    # chain A of 4HHB is UniProt P69905 with the initiator methionine removed
    trace_to_reference = {index: index + 1 for index in range(len(experimental))}
    predicted_to_reference = {index: index for index in range(len(predicted))}

    matching = match_residues(
        experimental,
        predicted,
        experimental_trace_to_reference=trace_to_reference,
        predicted_trace_to_reference=predicted_to_reference,
    )
    assert len(matching.uniprot_positions) == 141
    assert matching.predicted_only_positions == [1]

    metrics = evaluate_structure(matching, experimental, predicted)
    # A high-confidence AlphaFold model against a 1.74 Å crystal structure.
    assert 0.0 <= metrics.rmsd_angstrom < 2.0
    assert 0.8 <= metrics.lddt.score <= 1.0
    assert metrics.lddt.pairs > 1000


# --------------------------------------------------------------------------- #
# Confidence
# --------------------------------------------------------------------------- #
def test_validate_confidence_accepts_real_length() -> None:
    scores = ConfidenceScores(
        source="alphafold_confidence_json",
        residue_numbers=list(range(1, 143)),
        scores=[90.0] * 142,
    )
    mapping = validate_confidence(scores, sequence_length=142, sequence_start=1)
    assert len(mapping) == 142
    assert mapping[1] == 90.0


def test_validate_confidence_length_mismatch() -> None:
    scores = ConfidenceScores(
        source="alphafold_confidence_json",
        residue_numbers=list(range(1, 142)),
        scores=[90.0] * 141,
    )
    with pytest.raises(ConfidenceError) as excinfo:
        validate_confidence(scores, sequence_length=142)
    assert excinfo.value.code is ErrorCode.CONFIDENCE_LENGTH_MISMATCH


def test_validate_confidence_offset_numbering() -> None:
    scores = ConfidenceScores(
        source="alphafold_confidence_json",
        residue_numbers=list(range(2, 12)),
        scores=[80.0] * 10,
    )
    mapping = validate_confidence(scores, sequence_length=10, sequence_start=2)
    assert sorted(mapping) == list(range(2, 12))


def test_summarise_confidence_bins_and_correlation() -> None:
    plddt = {position: float(position) * 10 for position in range(1, 11)}
    lddt = {position: float(position) / 10 for position in range(1, 11)}
    summary = summarise_confidence(plddt, lddt, source="unit-test")
    assert summary.residue_count == 10
    assert sum(entry.residue_count for entry in summary.bins) == 10
    assert all(0.0 <= entry.fraction <= 1.0 for entry in summary.bins)
    assert summary.pearson_plddt_lddt == pytest.approx(1.0)
    assert summary.spearman_plddt_lddt == pytest.approx(1.0)


def test_cross_check_confidence_reports_max_delta() -> None:
    primary = {1: 90.0, 2: 80.0, 3: 70.0}
    secondary = {1: 90.5, 2: 79.0, 4: 60.0}
    assert cross_check_confidence(primary, secondary) == pytest.approx(1.0)
    assert cross_check_confidence({1: 1.0}, {}) is None


def test_plddt_from_trace_uses_bfactors(tmp_path) -> None:
    path = _write(tmp_path, "AF-P69905-F1-model_v6.cif")
    trace = extract_ca_trace(path, label="alphafold", chain_id="A")
    scores = plddt_from_trace(trace)
    assert len(scores.scores) == 142
    assert scores.source == "alphafold_cif_bfactor"
    assert 60.0 < min(scores.scores) < 100.0
