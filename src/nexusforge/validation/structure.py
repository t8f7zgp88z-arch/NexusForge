"""Structural validation: Cα extraction, rigid-body superposition, lDDT-Cα.

Scientific boundaries enforced here:

* ``pLDDT`` (AlphaFold's per-residue confidence) and ``lDDT`` (a distance test
  against an experimental reference) are different quantities. This module only
  computes **lDDT-Cα**, and every label says so.
* The score is Cα-only, so it is *not* a full-atom lDDT. Per-atom scoring is not
  implemented in v0.1 and is reported as not evaluated.
* Superposition is used for RMSD only. lDDT is superposition-free by definition.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from Bio.Data.IUPACData import protein_letters_3to1_extended
from Bio.PDB import MMCIFParser
from Bio.PDB.PDBExceptions import PDBConstructionException
from Bio.PDB.Polypeptide import is_aa

from nexusforge.errors import ErrorCode, StructureError
from nexusforge.models import CaChainSummary

_WATER = frozenset({"HOH", "WAT", "DOD", "H2O"})
#: Biopython keys the table by title case (``Ala``); PDB files use upper case.
_RESIDUE_ONE_LETTER = {
    name.upper(): letter for name, letter in protein_letters_3to1_extended.items()
}


@dataclass(frozen=True)
class CaTrace:
    """An ordered Cα trace with its per-residue metadata."""

    label: str
    chain_id: str | None
    residue_numbers: list[int]
    residue_names: list[str]
    one_letter: list[str]
    coords: np.ndarray
    bfactors: np.ndarray
    non_standard: tuple[str, ...] = ()
    missing_ca: tuple[int, ...] = ()
    duplicate_numbers: tuple[int, ...] = ()

    @property
    def sequence(self) -> str:
        return "".join(self.one_letter)

    def __len__(self) -> int:
        return len(self.one_letter)

    def summary(self) -> CaChainSummary:
        return CaChainSummary(
            label=self.label,
            chain_id=self.chain_id,
            residue_count=len(self),
            first_residue_number=self.residue_numbers[0] if self.residue_numbers else None,
            last_residue_number=self.residue_numbers[-1] if self.residue_numbers else None,
            non_standard_residues=list(self.non_standard),
            duplicate_residue_numbers=list(self.duplicate_numbers),
        )


@dataclass(frozen=True)
class ResidueMatching:
    """Residues that exist in both structures, expressed on the UniProt frame."""

    uniprot_positions: list[int]
    experimental_indices: list[int]
    predicted_indices: list[int]
    experimental_only_positions: list[int]
    predicted_only_positions: list[int]


@dataclass(frozen=True)
class LddtResult:
    score: float
    per_threshold: dict[str, float]
    per_residue: np.ndarray
    pairs: int
    scored_residues: int
    unscored_residues: int
    radius: float
    thresholds: tuple[float, ...]


@dataclass(frozen=True)
class StructureMetrics:
    rmsd_angstrom: float
    rotation: np.ndarray
    translation: np.ndarray
    lddt: LddtResult


# --------------------------------------------------------------------------- #
# Extraction
# --------------------------------------------------------------------------- #
def extract_ca_trace(
    path: Path,
    *,
    label: str,
    chain_id: str | None = None,
) -> CaTrace:
    """Parse an mmCIF file and extract one chain's Cα trace.

    Residues without a Cα atom are skipped and reported, never imputed.
    """
    parser = MMCIFParser(QUIET=True)
    try:
        structure = parser.get_structure(label, str(path))
    except (PDBConstructionException, ValueError, OSError, IndexError, KeyError) as exc:
        raise StructureError(
            ErrorCode.INVALID_COORDINATES,
            f"could not parse {label} structure file as mmCIF: {type(exc).__name__}",
            details={"path": str(path), "cause": str(exc)},
        ) from exc

    models = list(structure)
    if not models:
        raise StructureError(
            ErrorCode.INVALID_COORDINATES,
            f"{label} structure file contains no model",
            details={"path": str(path)},
        )
    model = models[0]

    chains = list(model)
    if not chains:
        raise StructureError(
            ErrorCode.INVALID_COORDINATES,
            f"{label} structure file contains no chains",
            details={"path": str(path)},
        )

    if chain_id is not None:
        if chain_id not in model:
            raise StructureError(
                ErrorCode.INVALID_COORDINATES,
                f"chain {chain_id!r} is not present in the {label} structure file "
                f"(chains: {sorted(c.id for c in chains)})",
                details={"path": str(path), "chain_id": chain_id},
            )
        chain = model[chain_id]
        resolved_chain_id = chain_id
    else:
        # Deterministic choice: the chain with the most amino-acid residues.
        chain = max(chains, key=lambda c: sum(1 for r in c if is_aa(r, standard=False)))
        resolved_chain_id = chain.id

    residue_numbers: list[int] = []
    residue_names: list[str] = []
    one_letter: list[str] = []
    coordinates: list[list[float]] = []
    bfactors: list[float] = []
    non_standard: list[str] = []
    missing_ca: list[int] = []
    seen_numbers: set[int] = set()
    duplicates: list[int] = []

    for residue in chain:
        resname = residue.get_resname().strip().upper()
        if resname in _WATER or resname not in _RESIDUE_ONE_LETTER:
            if resname not in _WATER:
                non_standard.append(resname)
            continue
        residue_number = int(residue.id[1])
        if "CA" not in residue:
            missing_ca.append(residue_number)
            continue
        calcium = residue["CA"]
        coord = np.asarray(calcium.get_coord(), dtype=float)
        if coord.shape != (3,) or not np.isfinite(coord).all():
            raise StructureError(
                ErrorCode.INVALID_COORDINATES,
                f"{label} residue {residue_number} has a non-finite Cα coordinate",
                details={"path": str(path), "residue_number": residue_number},
            )
        if float(np.max(np.abs(coord))) > 1.0e6:
            raise StructureError(
                ErrorCode.INVALID_COORDINATES,
                f"{label} residue {residue_number} has an implausible Cα coordinate "
                f"({coord.tolist()})",
                details={"path": str(path), "residue_number": residue_number},
            )
        if residue_number in seen_numbers:
            duplicates.append(residue_number)
        seen_numbers.add(residue_number)
        residue_numbers.append(residue_number)
        residue_names.append(resname)
        one_letter.append(_RESIDUE_ONE_LETTER[resname])
        coordinates.append(coord.tolist())
        bfactors.append(float(calcium.get_bfactor()))

    if not coordinates:
        raise StructureError(
            ErrorCode.INVALID_COORDINATES,
            f"{label} chain {resolved_chain_id!r} has no Cα atoms",
            details={"path": str(path), "chain_id": resolved_chain_id},
        )

    coords = np.asarray(coordinates, dtype=float)
    if coords.ndim != 2 or coords.shape[1] != 3 or not np.isfinite(coords).all():
        raise StructureError(
            ErrorCode.INVALID_COORDINATES,
            f"{label} Cα coordinate array is malformed (shape {coords.shape})",
            details={"path": str(path)},
        )
    if len(coords) < 3:
        raise StructureError(
            ErrorCode.INVALID_COORDINATES,
            f"{label} chain {resolved_chain_id!r} has only {len(coords)} Cα atom(s); "
            "at least three are required for a structural comparison",
            details={"path": str(path), "chain_id": resolved_chain_id},
        )

    return CaTrace(
        label=label,
        chain_id=resolved_chain_id,
        residue_numbers=residue_numbers,
        residue_names=residue_names,
        one_letter=one_letter,
        coords=coords,
        bfactors=np.asarray(bfactors, dtype=float),
        non_standard=tuple(sorted(set(non_standard))),
        missing_ca=tuple(missing_ca),
        duplicate_numbers=tuple(duplicates),
    )


# --------------------------------------------------------------------------- #
# Residue matching
# --------------------------------------------------------------------------- #
def match_residues(
    experimental: CaTrace,
    predicted: CaTrace,
    *,
    experimental_trace_to_reference: dict[int, int],
    predicted_trace_to_reference: dict[int, int],
) -> ResidueMatching:
    """Intersect two trace→reference index maps on shared reference positions.

    Reference positions are 0-based UniProt sequence indices. Both maps are
    produced by :func:`align_sequences` against the same reference sequence, so
    the intersection is meaningful even when one structure has missing residues.
    """
    exp_by_reference = {
        ref: index
        for index, ref in experimental_trace_to_reference.items()
        if 0 <= index < len(experimental)
    }
    pred_by_reference = {
        ref: index
        for index, ref in predicted_trace_to_reference.items()
        if 0 <= index < len(predicted)
    }
    shared = sorted(set(exp_by_reference) & set(pred_by_reference))
    if len(shared) < 3:
        raise StructureError(
            ErrorCode.NO_BASELINE,
            f"only {len(shared)} residue position(s) are resolved in both the experimental "
            "chain and the prediction; a structural comparison is not meaningful",
            details={
                "shared_residues": len(shared),
                "experimental_residues": len(experimental),
                "predicted_residues": len(predicted),
            },
        )
    return ResidueMatching(
        uniprot_positions=[ref + 1 for ref in shared],
        experimental_indices=[exp_by_reference[ref] for ref in shared],
        predicted_indices=[pred_by_reference[ref] for ref in shared],
        experimental_only_positions=sorted(
            ref + 1 for ref in set(exp_by_reference) - set(pred_by_reference)
        ),
        predicted_only_positions=sorted(
            ref + 1 for ref in set(pred_by_reference) - set(exp_by_reference)
        ),
    )


# --------------------------------------------------------------------------- #
# Geometry
# --------------------------------------------------------------------------- #
def kabsch_superpose(mobile: np.ndarray, reference: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Rotation ``R`` and translation ``t`` minimising ``||R·mobile + t − reference||``."""
    mobile = np.asarray(mobile, dtype=float)
    reference = np.asarray(reference, dtype=float)
    if mobile.shape != reference.shape or mobile.ndim != 2 or mobile.shape[1] != 3:
        raise StructureError(
            ErrorCode.INVALID_COORDINATES,
            f"superposition needs matching (N, 3) arrays, got {mobile.shape} and {reference.shape}",
        )
    if len(mobile) < 3:
        raise StructureError(
            ErrorCode.INVALID_COORDINATES,
            "at least three atom pairs are required for a rigid-body superposition",
        )
    mobile_centroid = mobile.mean(axis=0)
    reference_centroid = reference.mean(axis=0)
    mobile_centered = mobile - mobile_centroid
    reference_centered = reference - reference_centroid
    covariance = mobile_centered.T @ reference_centered
    u_matrix, _singular, vt_matrix = np.linalg.svd(covariance)
    determinant = np.sign(np.linalg.det(vt_matrix.T @ u_matrix.T))
    correction = np.diag([1.0, 1.0, determinant])
    rotation = vt_matrix.T @ correction @ u_matrix.T
    translation = reference_centroid - rotation @ mobile_centroid
    return rotation, translation


def rmsd(
    mobile: np.ndarray,
    reference: np.ndarray,
    rotation: np.ndarray,
    translation: np.ndarray,
) -> float:
    """RMSD after applying ``rotation``/``translation`` to ``mobile``."""
    transformed = (np.asarray(mobile, dtype=float) @ rotation.T) + translation
    difference = transformed - np.asarray(reference, dtype=float)
    return float(np.sqrt(np.mean(np.sum(difference * difference, axis=1))))


def compute_lddt_ca(
    experimental_coords: np.ndarray,
    predicted_coords: np.ndarray,
    *,
    inclusion_radius: float = 15.0,
    thresholds: tuple[float, ...] = (0.5, 1.0, 2.0, 4.0),
) -> LddtResult:
    """Cα-based local Distance Difference Test.

    For every residue pair closer than ``inclusion_radius`` **in the
    experimental structure**, the absolute difference between the experimental
    and predicted Cα–Cα distances is compared with each threshold. A residue's
    score is the mean fraction of conserved pairs across thresholds; the global
    score is the mean over scored residues. This is the standard lDDT definition
    restricted to Cα atoms, hence the name ``lDDT-Cα``.
    """
    experimental_coords = np.asarray(experimental_coords, dtype=float)
    predicted_coords = np.asarray(predicted_coords, dtype=float)
    if experimental_coords.shape != predicted_coords.shape:
        raise StructureError(
            ErrorCode.INVALID_COORDINATES,
            f"lDDT-Cα needs paired coordinate arrays, got {experimental_coords.shape} and "
            f"{predicted_coords.shape}",
        )
    if experimental_coords.ndim != 2 or experimental_coords.shape[1] != 3:
        raise StructureError(
            ErrorCode.INVALID_COORDINATES,
            f"coordinate arrays must be (N, 3), got {experimental_coords.shape}",
        )
    if not np.isfinite(experimental_coords).all() or not np.isfinite(predicted_coords).all():
        raise StructureError(
            ErrorCode.INVALID_COORDINATES, "coordinate arrays contain non-finite values"
        )
    thresholds = tuple(sorted(float(t) for t in thresholds))
    if not thresholds or any(t <= 0 for t in thresholds):
        raise StructureError(
            ErrorCode.INVALID_COORDINATES, "lDDT thresholds must be positive distances"
        )
    n_residues = len(experimental_coords)
    if n_residues < 2:
        raise StructureError(
            ErrorCode.NO_BASELINE,
            "at least two residues are required to evaluate lDDT-Cα",
        )

    difference_exp = experimental_coords[:, None, :] - experimental_coords[None, :, :]
    distances_exp = np.sqrt(np.sum(difference_exp * difference_exp, axis=-1))
    difference_pred = predicted_coords[:, None, :] - predicted_coords[None, :, :]
    distances_pred = np.sqrt(np.sum(difference_pred * difference_pred, axis=-1))

    within_radius = distances_exp <= inclusion_radius
    np.fill_diagonal(within_radius, False)
    pair_count = int(within_radius.sum() // 2)
    if pair_count == 0:
        raise StructureError(
            ErrorCode.NO_BASELINE,
            f"no Cα pairs fall within {inclusion_radius} Å in the experimental structure; "
            "lDDT-Cα is not defined for this pair of structures",
            details={"inclusion_radius": inclusion_radius, "residues": n_residues},
        )

    distance_error = np.abs(distances_pred - distances_exp)
    per_threshold_conserved = np.stack(
        [(distance_error < threshold) & within_radius for threshold in thresholds]
    )  # (T, N, N)

    per_residue = np.full(n_residues, np.nan, dtype=float)
    for index in range(n_residues):
        partners = within_radius[index]
        if not partners.any():
            continue
        per_threshold_scores = per_threshold_conserved[:, index, partners].mean(axis=1)
        per_residue[index] = float(per_threshold_scores.mean())

    scored = int(np.count_nonzero(~np.isnan(per_residue)))
    unscored = n_residues - scored
    if scored == 0:
        raise StructureError(
            ErrorCode.NO_BASELINE,
            "no residue has a neighbour within the inclusion radius; lDDT-Cα is undefined",
        )

    per_threshold: dict[str, float] = {}
    for threshold_index, threshold in enumerate(thresholds):
        values = np.full(n_residues, np.nan, dtype=float)
        for index in range(n_residues):
            partners = within_radius[index]
            if not partners.any():
                continue
            values[index] = float(per_threshold_conserved[threshold_index, index, partners].mean())
        per_threshold[f"{threshold:g}"] = float(np.nanmean(values))

    return LddtResult(
        score=float(np.nanmean(per_residue)),
        per_threshold=per_threshold,
        per_residue=per_residue,
        pairs=pair_count,
        scored_residues=scored,
        unscored_residues=unscored,
        radius=inclusion_radius,
        thresholds=thresholds,
    )


def evaluate_structure(
    matching: ResidueMatching,
    experimental: CaTrace,
    predicted: CaTrace,
    *,
    inclusion_radius: float = 15.0,
    thresholds: tuple[float, ...] = (0.5, 1.0, 2.0, 4.0),
) -> StructureMetrics:
    """Superpose on shared Cα atoms, then score local distance agreement."""
    experimental_coords = experimental.coords[matching.experimental_indices]
    predicted_coords = predicted.coords[matching.predicted_indices]
    rotation, translation = kabsch_superpose(predicted_coords, experimental_coords)
    rmsd_value = rmsd(predicted_coords, experimental_coords, rotation, translation)
    lddt_result = compute_lddt_ca(
        experimental_coords,
        predicted_coords,
        inclusion_radius=inclusion_radius,
        thresholds=thresholds,
    )
    return StructureMetrics(
        rmsd_angstrom=rmsd_value,
        rotation=rotation,
        translation=translation,
        lddt=lddt_result,
    )


__all__ = [
    "CaTrace",
    "LddtResult",
    "ResidueMatching",
    "StructureMetrics",
    "compute_lddt_ca",
    "evaluate_structure",
    "extract_ca_trace",
    "kabsch_superpose",
    "match_residues",
    "rmsd",
]
