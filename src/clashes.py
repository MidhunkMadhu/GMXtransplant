"""
PBC-aware steric clash detection and whole-molecule removal.

Flags any environment residue (lipid / water / ion / other) that has at
least one atom within a configurable, per-class distance threshold of the
inserted protein+ligand block, then removes that residue in its
entirety -- never a partial residue -- using minimum-image-convention
distances so clashes across periodic boundaries are not missed.
"""

from __future__ import annotations

import warnings
import os
import tempfile
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import MDAnalysis as mda
    from MDAnalysis.lib.distances import capped_distance

from classify import classify_resnames, is_hydrogen_name
from clash_logic import closest_clash_per_residue
from config import ClashDetectionSpec


def _serialized_coordinates(atomgroup, box_dimensions, suffix):
    """Use the real writer/parser to reproduce coordinate and cell rounding."""
    copy = mda.Merge(atomgroup)
    copy.dimensions = box_dimensions
    with tempfile.TemporaryDirectory(prefix="clash-rounding-") as directory:
        path = os.path.join(directory, "coordinates." + suffix)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            copy.atoms.write(path)
            loaded = mda.Universe(path)
        return loaded.atoms.positions.copy(), loaded.dimensions


@dataclass
class RemovedMolecule:
    resindex: int
    resname: str
    resid: int
    molclass: str
    min_distance: float
    n_atoms: int
    cutoff_used: float = 0.0
    # The specific atom pair that produced min_distance -- env side (e.g.
    # the lipid atom) and the inserted-block side (e.g. the protein/ligand
    # atom), so a flagged contact is actually actionable, not just "some
    # atom of this residue is close." Empty string if unavailable.
    contact_env_atom: str = ""
    contact_partner_resname: str = ""
    contact_partner_resid: int = -1
    contact_partner_atom: str = ""
    contact_context: str = "inserted_block"
    removal_decision: str = "remove"


@dataclass
class ClashResult:
    kept_ag: "mda.core.groups.AtomGroup"
    removed_ag: "mda.core.groups.AtomGroup"
    removed_molecules: List[RemovedMolecule]
    counts_by_class: Dict[str, int]
    counts_by_resname: Dict[str, int]
    threshold_used: Dict[str, float]
    n_lipids_removed: int
    n_lipids_total: int
    lipid_removal_fraction: float
    lipid_removal_flagged: bool
    other_class_resnames: List[str]
    # Residues that were within threshold for a class listed in
    # clash_cfg.keep_classes (e.g. "lipid") -- detected and reported, but
    # deliberately NOT removed. Same shape as removed_molecules.
    flagged_but_kept: List[RemovedMolecule] = field(default_factory=list)
    flagged_but_kept_counts_by_resname: Dict[str, int] = field(default_factory=dict)
    notes: List[str] = field(default_factory=list)


def detect_and_remove_clashes(
    environment_ag,
    inserted_ag,
    box_dimensions: np.ndarray,
    clash_cfg: ClashDetectionSpec,
    classification_overrides: Optional[Dict[str, str]] = None,
) -> ClashResult:
    notes: List[str] = []

    env_resname_class = classify_resnames(environment_ag.resnames, classification_overrides)
    other_resnames = sorted(rn for rn, cls in env_resname_class.items() if cls == "other")
    if other_resnames:
        notes.append(
            f"Environment resnames classified as 'other' (unrecognized lipid/water/ion): "
            f"{other_resnames}. They use the 'other' clash threshold "
            f"({clash_cfg.threshold_for('other')} A) unless overridden."
        )

    env_atom_class = np.array([env_resname_class[rn] for rn in environment_ag.resnames])

    if clash_cfg.heavy_atoms_only:
        env_heavy_mask = np.array([not is_hydrogen_name(n) for n in environment_ag.names])
        ins_heavy_mask = np.array([not is_hydrogen_name(n) for n in inserted_ag.names])
    else:
        env_heavy_mask = np.ones(len(environment_ag), dtype=bool)
        ins_heavy_mask = np.ones(len(inserted_ag), dtype=bool)

    env_query_ag = environment_ag[env_heavy_mask]
    ins_query_ag = inserted_ag[ins_heavy_mask]

    max_threshold = max(
        clash_cfg.threshold,
        *(clash_cfg.thresholds.values() or [clash_cfg.threshold]),
    )

    box = box_dimensions if clash_cfg.use_pbc else None

    if len(env_query_ag) == 0 or len(ins_query_ag) == 0:
        raise ValueError(
            "Clash detection requires at least one selected environment atom and "
            "one selected inserted-block atom after heavy-atom filtering."
        )

    keep_classes = set(clash_cfg.keep_classes)
    if keep_classes:
        notes.append(
            f"clash_detection.keep_classes = {sorted(keep_classes)}: residues of these "
            f"classes are still checked and reported when within threshold, but are "
            f"deliberately NOT removed."
        )

    min_dist_per_resindex, min_pair_per_resindex, contact_formats = {}, {}, {}
    for coordinate_format in ("original", "pdb", "gro"):
        if coordinate_format == "original":
            env_positions, ins_positions, query_box = env_query_ag.positions, ins_query_ag.positions, box
        else:
            env_positions, rounded_box = _serialized_coordinates(env_query_ag, box_dimensions, coordinate_format)
            ins_positions, _ = _serialized_coordinates(ins_query_ag, box_dimensions, coordinate_format)
            query_box = rounded_box if clash_cfg.use_pbc else None
        pairs, distances = capped_distance(
            env_positions, ins_positions, max_cutoff=max_threshold,
            box=query_box, return_distances=True,
        )
        minima, closest_pairs = closest_clash_per_residue(
            pairs, distances, env_atom_class[env_heavy_mask], env_query_ag.resindices,
            clash_cfg.threshold, clash_cfg.thresholds,
        )
        for ri, distance in minima.items():
            if ri not in min_dist_per_resindex or distance < min_dist_per_resindex[ri]:
                min_dist_per_resindex[ri] = distance
                min_pair_per_resindex[ri] = closest_pairs[ri]
                contact_formats[ri] = coordinate_format
    notes.append("Clashes checked before and after actual PDB/GRO serialization, including "
                 "cell rounding. Minimum reported distances identify the coordinate format; "
                 "whole-residue removal and keep_classes apply before charge/topology assembly.")

    resindex_to_resname = dict(zip(environment_ag.resindices, environment_ag.resnames))
    resindex_to_resid = dict(zip(environment_ag.resindices, environment_ag.resids))

    # Preserve the exact closest atom pair. Every removal candidate comes only
    # from a contact with the newly inserted protein/ligand block. Contacts
    # among molecules already present in the equilibrated environment are not
    # clash-removal candidates.
    closest_contact = {}
    for ri, distance in min_dist_per_resindex.items():
        env_q_idx, ins_q_idx = min_pair_per_resindex[ri]
        closest_contact[ri] = (
            float(distance), env_query_ag[env_q_idx], ins_query_ag[ins_q_idx],
            "inserted_block:" + contact_formats[ri],
        )

    flagged_resindices = set(closest_contact)

    def _make_record(ri: int, decision: str) -> RemovedMolecule:
        rn = resindex_to_resname[ri]
        cls = env_resname_class[rn]
        n_atoms = int(np.sum(environment_ag.resindices == ri))
        min_distance, env_atom, ins_atom, context = closest_contact[ri]
        return RemovedMolecule(
            resindex=ri, resname=rn, resid=int(resindex_to_resid[ri]),
            molclass=cls, min_distance=min_distance, n_atoms=n_atoms,
            cutoff_used=float(clash_cfg.threshold_for(cls)),
            contact_env_atom=str(env_atom.name),
            contact_partner_resname=str(ins_atom.resname),
            contact_partner_resid=int(ins_atom.resid),
            contact_partner_atom=str(ins_atom.name),
            contact_context=context,
            removal_decision=decision,
        )

    removal_resindices = set()
    keep_resindices = set()
    for ri in flagged_resindices:
        cls = env_resname_class[resindex_to_resname[ri]]
        (keep_resindices if cls in keep_classes else removal_resindices).add(ri)

    remove_atom_mask = np.isin(environment_ag.resindices, list(removal_resindices))
    removed_ag = environment_ag[remove_atom_mask]
    kept_ag = environment_ag[~remove_atom_mask]

    removed_molecules: List[RemovedMolecule] = []
    counts_by_class: Dict[str, int] = {}
    counts_by_resname: Dict[str, int] = {}
    for ri in sorted(removal_resindices):
        rec = _make_record(ri, "remove")
        removed_molecules.append(rec)
        counts_by_class[rec.molclass] = counts_by_class.get(rec.molclass, 0) + 1
        counts_by_resname[rec.resname] = counts_by_resname.get(rec.resname, 0) + 1

    flagged_but_kept: List[RemovedMolecule] = [
        _make_record(ri, "keep_configured_class") for ri in sorted(keep_resindices)
    ]
    flagged_but_kept_counts_by_resname: Dict[str, int] = {}
    for rec in flagged_but_kept:
        flagged_but_kept_counts_by_resname[rec.resname] = (
            flagged_but_kept_counts_by_resname.get(rec.resname, 0) + 1
        )
    if flagged_but_kept:
        by_class: Dict[str, int] = {}
        for rec in flagged_but_kept:
            by_class[rec.molclass] = by_class.get(rec.molclass, 0) + 1
        notes.append(
            f"{len(flagged_but_kept)} residue(s) within threshold but kept "
            f"(clash_detection.keep_classes): {by_class}. See flagged_but_kept for the "
            f"specific contact atoms -- closest distances range from "
            f"{min(r.min_distance for r in flagged_but_kept):.2f} to "
            f"{max(r.min_distance for r in flagged_but_kept):.2f} A."
        )

    n_lipids_total = int(np.sum(np.array([env_resname_class[rn] for rn in environment_ag.residues.resnames]) == "lipid"))
    n_lipids_removed = counts_by_class.get("lipid", 0)
    lipid_removal_fraction = (n_lipids_removed / n_lipids_total) if n_lipids_total else 0.0
    lipid_flagged = lipid_removal_fraction > clash_cfg.flag_lipid_removal_fraction
    if lipid_flagged:
        notes.append(
            f"FLAG: {n_lipids_removed}/{n_lipids_total} lipids removed "
            f"({lipid_removal_fraction:.1%}), exceeding the configured "
            f"flag_lipid_removal_fraction ({clash_cfg.flag_lipid_removal_fraction:.1%}). "
            f"This may leave a large membrane cavity -- recommend manual inspection."
        )

    threshold_used = {
        cls: clash_cfg.threshold_for(cls)
        for cls in set(list(env_resname_class.values()) + ["lipid", "water", "ion", "other"])
    }

    return ClashResult(
        kept_ag=kept_ag,
        removed_ag=removed_ag,
        removed_molecules=removed_molecules,
        counts_by_class=counts_by_class,
        counts_by_resname=counts_by_resname,
        threshold_used=threshold_used,
        n_lipids_removed=n_lipids_removed,
        n_lipids_total=n_lipids_total,
        lipid_removal_fraction=lipid_removal_fraction,
        lipid_removal_flagged=lipid_flagged,
        other_class_resnames=other_resnames,
        flagged_but_kept=flagged_but_kept,
        flagged_but_kept_counts_by_resname=flagged_but_kept_counts_by_resname,
        notes=notes,
    )
