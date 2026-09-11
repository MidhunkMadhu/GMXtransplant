"""Experimental-cholesterol restoration and composition correction."""

from __future__ import annotations

import os
import random
import tempfile
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from types import SimpleNamespace
from typing import Dict, List, Optional, Tuple

import numpy as np

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import MDAnalysis as mda
    from MDAnalysis.lib.distances import distance_array

from align import AlignmentResult, align_replacement_to_box
from charmm36_cholesterol import run as run_charmm36_conversion
from clashes import ClashResult, detect_and_remove_clashes
from classify import classify_resnames
from config import CholesterolSpec, ClashDetectionSpec
from masks import resolve_mask
from output import canonicalize_system, write_pdb
from box_validation import BoxValidationError, validate_box, read_environment


class CholesterolError(RuntimeError):
    pass


@dataclass
class CompositionResult:
    midplane_z: float
    box_volume_liters: float
    lipid_targets: Dict[str, Dict[str, int]]
    lipid_counts_before: Dict[str, Dict[str, int]]
    lipid_counts_after: Dict[str, Dict[str, int]]
    lipid_residues_removed: List[dict]
    reference_ion_counts: Dict[str, int]
    ion_counts_before: Dict[str, int]
    ion_target_counts: Dict[str, int]
    ion_counts_after: Dict[str, int]
    ions_removed: List[dict]
    neutralizing_excess: int
    target_concentration_molar: float
    reference_salt_pairs: int
    final_salt_pairs: int
    estimated_final_concentration_molar: Optional[float]
    concentration_deviation_fraction: Optional[float]
    concentration_tolerance_fraction: float
    concentration_within_tolerance: Optional[bool]
    output_ion_counts: Dict[str, Dict[str, int]] = field(default_factory=dict)
    warnings: List[str] = field(default_factory=list)


@dataclass
class CholesterolResult:
    converted_structure_path: Optional[str]
    n_experimental_cholesterols: int
    n_experimental_cholesterol_atoms: int
    alignment: AlignmentResult
    clash_result: ClashResult
    composition: CompositionResult
    merged_pdb_path: Optional[str]
    diagnostic_paths: List[str]
    heavy_atom_audit: List[dict] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)


def _ensure_parent(path: str) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)


def convert_experimental_cholesterols(spec: CholesterolSpec):
    """Load the experimental system after optional CHARMM36 conversion.

    Conversion products are retained only when write_diagnostics is enabled.
    Otherwise they live in a temporary directory and the returned Universe is
    detached from those files before the directory is removed.
    """
    if not spec.convert_to_charmm36:
        return mda.Universe(spec.experimental_structure_path), None, [], []

    def run_conversion(raw_path, fixed_path, aligned_path, placed_path):
        args = SimpleNamespace(
            input=Path(spec.experimental_structure_path),
            reference=(Path(spec.charmm36_reference_path)
                       if spec.charmm36_reference_path else None),
            raw_out=Path(raw_path),
            fixed_out=Path(fixed_path),
            aligned_reference_out=Path(aligned_path),
            reference_placed_out=Path(placed_path),
            chol_resnames=",".join(spec.cholesterol_resnames),
            obabel=spec.obabel_command,
            verify_tolerance=spec.verify_tolerance,
            reconstruct_missing_heavy_atoms=spec.reconstruct_missing_heavy_atoms,
            max_missing_heavy_atoms=spec.max_missing_heavy_atoms,
            max_heavy_atom_fit_rmsd=spec.max_heavy_atom_fit_rmsd,
        )
        try:
            summary = run_charmm36_conversion(args)
        except Exception as exc:
            raise CholesterolError(
                f"CHARMM36 cholesterol conversion failed: {exc}"
            ) from exc
        loaded = mda.Universe(str(fixed_path))
        # Merge copies topology and coordinates into an in-memory trajectory,
        # allowing temporary conversion files to be removed immediately.
        detached = mda.Merge(loaded.atoms)
        return detached, summary.get("heavy_atom_audit", [])

    if spec.write_diagnostics:
        paths = (
            spec.raw_protonated_pdb_path,
            spec.converted_pdb_path,
            spec.aligned_reference_pdb_path,
            spec.reference_placed_pdb_path,
        )
        for path in paths:
            _ensure_parent(path)
        experimental, audit = run_conversion(*paths)
        diagnostics = [paths[0], paths[2], paths[3]]
        return experimental, spec.converted_pdb_path, diagnostics, audit

    with tempfile.TemporaryDirectory(prefix="cholesterol_conversion_") as directory:
        scratch = Path(directory)
        paths = (
            scratch / "protonated_raw.pdb",
            scratch / "charmm36_converted.pdb",
            scratch / "aligned_reference.pdb",
            scratch / "reference_placed.pdb",
        )
        experimental, audit = run_conversion(*paths)
    return experimental, None, [], audit


def _load_target_with_box(spec: CholesterolSpec, validation_spec=None):
    try:
        u = read_environment(spec.target_system_path, spec.box_dimensions, validation_spec,
                             format=spec.target_system_format)
    except BoxValidationError as exc:
        raise CholesterolError(str(exc)) from exc
    return u, list(u.box_preflight["notes"])


def _write_alignment_inspection(target_receptor, aligned_receptor, cholesterols, path: str) -> None:
    _ensure_parent(path)
    inspection = mda.Merge(target_receptor, aligned_receptor, cholesterols)
    chain_ids = (
        ["X"] * len(target_receptor)
        + ["Y"] * len(aligned_receptor)
        + ["Z"] * len(cholesterols)
    )
    inspection.add_TopologyAttr("chainIDs", chain_ids)
    inspection.atoms.write(path)


def _residue_reference_z(residue, charmm_resname: str) -> float:
    preferred = ("O3",) if residue.resname == charmm_resname else ("P", "P1")
    for atom in residue.atoms:
        if atom.name in preferred:
            return float(atom.position[2])
    return float(residue.atoms.positions[:, 2].mean())


def _lipid_midplane(universe) -> float:
    classes = classify_resnames(universe.residues.resnames)
    lipid_resindices = {
        res.resindex for res in universe.residues if classes[res.resname] == "lipid"
    }
    if not lipid_resindices:
        raise CholesterolError("No membrane lipids were recognized in the merged system.")
    ref_atoms = universe.atoms[
        np.isin(universe.atoms.resindices, list(lipid_resindices))
        & np.isin(universe.atoms.names, ["P", "P1"])
    ]
    if len(ref_atoms):
        return float(ref_atoms.positions[:, 2].mean())
    lipid_atoms = universe.atoms[np.isin(universe.atoms.resindices, list(lipid_resindices))]
    return float(lipid_atoms.positions[:, 2].mean())


def _leaflet_counts(universe, midplane_z: float, charmm_resname: str) -> Dict[str, Dict[str, int]]:
    classes = classify_resnames(universe.residues.resnames)
    counts: Dict[str, Dict[str, int]] = {}
    for res in universe.residues:
        if classes[res.resname] != "lipid":
            continue
        leaf = "upper" if _residue_reference_z(res, charmm_resname) > midplane_z else "lower"
        counts.setdefault(res.resname, {"upper": 0, "lower": 0})[leaf] += 1
    return counts


def _protein_atoms(universe):
    classes = classify_resnames(universe.residues.resnames)
    protein_resindices = [
        res.resindex for res in universe.residues if classes[res.resname] == "protein"
    ]
    return universe.atoms[np.isin(universe.atoms.resindices, protein_resindices)]


def _eligible_far_from_protein(residues, protein_atoms, cutoff: float):
    if len(protein_atoms) == 0:
        raise CholesterolError("No protein residues were recognized for distance filtering.")
    eligible = []
    box = protein_atoms.universe.dimensions
    if box is None or len(box) != 6:
        raise CholesterolError(
            "PBC-aware protein-distance filtering requires six valid box dimensions."
        )
    for res in residues:
        distances = distance_array(
            res.atoms.positions, protein_atoms.positions, box=box
        )
        if float(distances.min()) > cutoff:
            eligible.append(res)
    return eligible


def _box_volume_liters(dimensions) -> float:
    try:
        validation = validate_box(dimensions, label="cholesterol target box")
    except BoxValidationError as exc:
        raise CholesterolError(str(exc)) from exc
    return validation.volume * 1e-27


def count_salt_ions(universe, cation: str, anion: str) -> Dict[str, int]:
    """Count the configured monatomic salt species by residue."""
    cation = cation.upper()
    anion = anion.upper()
    return {
        cation: sum(str(res.resname).upper() == cation for res in universe.residues),
        anion: sum(str(res.resname).upper() == anion for res in universe.residues),
    }


def compare_salt_concentration(
    reference_counts: Dict[str, int],
    final_counts: Dict[str, int],
    cation: str,
    anion: str,
    reference_concentration_molar: float,
    tolerance_fraction: float,
) -> Dict[str, object]:
    """Infer final molarity from salt-pair retention relative to a reference."""
    cation = cation.upper()
    anion = anion.upper()
    reference_pairs = min(reference_counts.get(cation, 0), reference_counts.get(anion, 0))
    final_pairs = min(final_counts.get(cation, 0), final_counts.get(anion, 0))
    if reference_pairs <= 0:
        return {
            "reference_salt_pairs": reference_pairs,
            "final_salt_pairs": final_pairs,
            "estimated_final_concentration_molar": None,
            "concentration_deviation_fraction": None,
            "concentration_within_tolerance": None,
        }
    retained_fraction = final_pairs / reference_pairs
    deviation = abs(retained_fraction - 1.0)
    return {
        "reference_salt_pairs": reference_pairs,
        "final_salt_pairs": final_pairs,
        "estimated_final_concentration_molar": (
            reference_concentration_molar * retained_fraction
        ),
        "concentration_deviation_fraction": deviation,
        "concentration_within_tolerance": deviation <= tolerance_fraction,
    }


def adjust_composition(
    universe, reference_universe, spec: CholesterolSpec,
    protein_chain_atom_counts=None,
):
    comp = spec.composition
    rng = random.Random(comp.random_seed)
    midplane_z = _lipid_midplane(universe)
    before_lipids = _leaflet_counts(universe, midplane_z, spec.charmm_resname)
    protein_atoms = _protein_atoms(universe)
    classes = classify_resnames(universe.residues.resnames)
    remove_resindices = set()
    lipid_removed: List[dict] = []
    notes: List[str] = []

    for resname, leaves in sorted(comp.lipid_targets.items()):
        for leaf in ("upper", "lower"):
            target = int(leaves[leaf])
            candidates = [
                res for res in universe.residues
                if classes[res.resname] == "lipid"
                and res.resname.upper() == resname.upper()
                and ("upper" if _residue_reference_z(res, spec.charmm_resname) > midplane_z else "lower") == leaf
            ]
            surplus = len(candidates) - target
            if surplus < 0:
                notes.append(
                    f"{resname} {leaf} is {abs(surplus)} residue(s) below target; "
                    "the mode does not add lipids."
                )
                continue
            if surplus == 0:
                continue
            eligible = _eligible_far_from_protein(
                candidates, protein_atoms, comp.distance_from_protein
            )
            n_pick = min(surplus, len(eligible))
            if n_pick < surplus:
                notes.append(
                    f"{resname} {leaf}: {surplus} removal(s) were requested, but only "
                    f"{len(eligible)} residue(s) were farther than "
                    f"{comp.distance_from_protein:.1f} A from the protein."
                )
            for res in rng.sample(eligible, n_pick):
                remove_resindices.add(res.resindex)
                lipid_removed.append({
                    "resindex": int(res.resindex),
                    "resid": int(res.resid),
                    "resname": str(res.resname),
                    "leaflet": leaf,
                })

    cation = comp.salt_cation.upper()
    anion = comp.salt_anion.upper()
    ref_counts = count_salt_ions(reference_universe, cation, anion)
    neutralizing_excess = ref_counts[cation] - ref_counts[anion]
    volume_l = _box_volume_liters(universe.dimensions)
    ion_before = count_salt_ions(universe, cation, anion)
    ions_removed: List[dict] = []

    keep_mask = ~np.isin(universe.atoms.resindices, list(remove_resindices))
    final_universe = mda.Merge(universe.atoms[keep_mask])
    final_universe = canonicalize_system(
        final_universe, universe.dimensions, protein_chain_atom_counts
    )
    after_lipids = _leaflet_counts(final_universe, midplane_z, spec.charmm_resname)
    ion_after = count_salt_ions(final_universe, cation, anion)
    concentration = compare_salt_concentration(
        ref_counts,
        ion_after,
        cation,
        anion,
        comp.target_concentration,
        comp.concentration_tolerance_fraction,
    )
    if concentration["concentration_within_tolerance"] is None:
        notes.append(
            "Salt concentration could not be compared because the reference "
            f"contains no complete {cation}/{anion} salt pairs."
        )
    elif not concentration["concentration_within_tolerance"]:
        notes.append(
            "Final salt concentration differs from the reference by "
            f"{concentration['concentration_deviation_fraction']:.1%}, exceeding "
            f"the accepted {comp.concentration_tolerance_fraction:.1%} tolerance "
            f"({concentration['reference_salt_pairs']} reference pairs versus "
            f"{concentration['final_salt_pairs']} final pairs; estimated "
            f"{concentration['estimated_final_concentration_molar']:.4f} M versus "
            f"the assumed {comp.target_concentration:.4f} M reference)."
        )

    result = CompositionResult(
        midplane_z=midplane_z,
        box_volume_liters=volume_l,
        lipid_targets=comp.lipid_targets,
        lipid_counts_before=before_lipids,
        lipid_counts_after=after_lipids,
        lipid_residues_removed=lipid_removed,
        reference_ion_counts=ref_counts,
        ion_counts_before=ion_before,
        # Retained for report/API compatibility. The authoritative target ion
        # composition is now the reference structure itself.
        ion_target_counts=dict(ref_counts),
        ion_counts_after=ion_after,
        ions_removed=ions_removed,
        neutralizing_excess=neutralizing_excess,
        target_concentration_molar=comp.target_concentration,
        reference_salt_pairs=concentration["reference_salt_pairs"],
        final_salt_pairs=concentration["final_salt_pairs"],
        estimated_final_concentration_molar=concentration[
            "estimated_final_concentration_molar"
        ],
        concentration_deviation_fraction=concentration[
            "concentration_deviation_fraction"
        ],
        concentration_tolerance_fraction=comp.concentration_tolerance_fraction,
        concentration_within_tolerance=concentration[
            "concentration_within_tolerance"
        ],
        warnings=notes,
    )
    return final_universe, result


def restore_cholesterols(cfg, workdir: str = ".", receptor_layout=None):
    """Run cholesterol preparation through final composition correction."""
    spec = cfg.cholesterol
    # Audit/repair and map cholesterol before touching the target environment.
    # This keeps malformed experimental cholesterol as the true first-stage
    # failure and makes the workflow order match its user-facing description.
    experimental, converted_path, diagnostics, heavy_atom_audit = (
        convert_experimental_cholesterols(spec)
    )
    target, notes = _load_target_with_box(spec, cfg.box_validation)

    align_res = align_replacement_to_box(
        target,
        experimental,
        spec.target_receptor_mask,
        spec.experimental_receptor_mask,
        cfg.alignment,
        workdir=workdir,
    )

    original_positions = experimental.atoms.positions.astype(np.float64).copy()
    experimental.atoms.positions = (
        (align_res.rotation @ original_positions.T).T + align_res.translation
    )
    aligned_receptor = resolve_mask(experimental, spec.experimental_receptor_mask)
    aligned_receptor.positions = align_res.new_positions

    recognized = set(spec.cholesterol_resnames) | {spec.charmm_resname}
    chol_residues = experimental.residues[
        np.isin(np.char.upper(experimental.residues.resnames.astype(str)), list(recognized))
    ]
    if len(chol_residues) == 0:
        cholesterol_source = (
            converted_path or f"converted {spec.experimental_structure_path}"
        )
        raise CholesterolError(
            f"No cholesterol residues with names {sorted(recognized)} were found "
            f"in {cholesterol_source}."
        )
    chol_u = mda.Merge(chol_residues.atoms)
    chol_u.residues.resnames = [spec.charmm_resname] * len(chol_u.residues)

    target_receptor = resolve_mask(target, spec.target_receptor_mask)
    if spec.write_diagnostics:
        _write_alignment_inspection(
            target_receptor, aligned_receptor, chol_u.atoms,
            cfg.output.inspection_pdb_path,
        )

    from dataclasses import replace as dataclass_replace
    protected_classes = sorted(set(cfg.clash_detection.keep_classes) | {"protein", "other"})
    clash_cfg: ClashDetectionSpec = dataclass_replace(
        cfg.clash_detection, keep_classes=protected_classes
    )
    clash_res = detect_and_remove_clashes(
        target.atoms, chol_u.atoms, target.dimensions, clash_cfg
    )

    merged = mda.Merge(clash_res.kept_ag, chol_u.atoms)
    chain_counts = receptor_layout.atom_counts if receptor_layout else None
    merged = canonicalize_system(merged, target.dimensions, chain_counts)
    if spec.write_diagnostics:
        write_pdb(merged, spec.merged_pdb_path, chain_counts)

    ref_path = spec.composition.reference_system_path or spec.target_system_path
    reference = target if os.path.abspath(ref_path) == os.path.abspath(spec.target_system_path) \
        else mda.Universe(ref_path)
    final_universe, composition = adjust_composition(
        merged, reference, spec, chain_counts
    )

    result = CholesterolResult(
        converted_structure_path=converted_path,
        n_experimental_cholesterols=len(chol_u.residues),
        n_experimental_cholesterol_atoms=len(chol_u.atoms),
        alignment=align_res,
        clash_result=clash_res,
        composition=composition,
        merged_pdb_path=(spec.merged_pdb_path if spec.write_diagnostics else None),
        diagnostic_paths=diagnostics,
        heavy_atom_audit=heavy_atom_audit,
        notes=notes,
    )
    return target, final_universe, result
