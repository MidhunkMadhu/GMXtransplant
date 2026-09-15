"""
Ligand-only replacement mode.

A special case of the whole workflow: the protein, membrane, water, and
ions do NOT change at all. Only one ligand inside an existing, otherwise
untouched structure is swapped for a different one, supplied as its own
separate coordinate file (optionally with its own topology .itp).

Three ways to position the incoming ligand:
  - "pairfit"  -- rigid-body superpose the new ligand onto the OLD ligand's
                  pose, using an explicit, user-supplied list of matched
                  atom names (one list per side, paired positionally). This
                  is a plain Kabsch/SVD best-fit rotation+translation over
                  just those atoms, then applied to every atom of the new
                  ligand -- no sequence alignment, no PyMOL involved (that
                  machinery in align.py is for the protein's much larger,
                  sequence-based alignment problem; a handful of named atom
                  pairs on a small molecule doesn't need it).
  - "autofit"  -- automatically match every heavy atom by its unique atom
                  name and calculate a rigid-body Kabsch/SVD best fit. The
                  old and new ligands must have identical heavy-atom name
                  sets; file order and hydrogen names/counts may differ.
  - "nofit"    -- use the new ligand's own coordinates exactly as given, no
                  transformation at all. For when you've already positioned
                  it yourself (e.g. docked directly into this structure).

After the swap, the SAME downstream steps as the protein-replacement flow
apply unchanged: PBC-aware clash detection (with the untouched protein
itself automatically protected -- see run_pipeline.py), charge accounting
and neutralization, PDB/GRO output, and (if enabled) topology/index
generation.
"""

from __future__ import annotations

import warnings
import re
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import numpy as np

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import MDAnalysis as mda

from config import LigandReplaceSpec
from classify import is_hydrogen_name
from itp import ITPParseError, parse_itp
from kabsch import kabsch_fit as _kabsch_fit


class LigandReplaceError(RuntimeError):
    pass


def select_by_resname(universe: "mda.Universe", selector: str, label: str):
    """Select exactly one ligand residue by ``RES`` or ``RES:N``.

    ``N`` is the 1-based positional occurrence among residues named ``RES``
    in coordinate-file order. A bare name is accepted only when exactly one
    residue has that name.
    """
    match = re.fullmatch(r"([A-Za-z0-9_+\-]+)(?::([1-9]\d*))?", selector.strip())
    if not match:
        raise LigandReplaceError(
            f"{label} '{selector}' must use RES or RES:N syntax, for example "
            "LDP or LDP:2. N is a 1-based positional occurrence."
        )
    resname, occurrence_text = match.groups()
    matches = [res for res in universe.residues if str(res.resname) == resname]
    if not matches:
        present = sorted(set(universe.atoms.resnames))
        raise LigandReplaceError(
            f"{label} resname '{resname}' matched zero atoms. Resnames present in this "
            f"file: {present}"
        )
    if occurrence_text is None:
        if len(matches) != 1:
            choices = ", ".join(f"{resname}:{i}" for i in range(1, len(matches) + 1))
            raise LigandReplaceError(
                f"{label} '{resname}' matches {len(matches)} residues. Select exactly "
                f"one positional occurrence: {choices}"
            )
        selected = matches[0]
    else:
        occurrence = int(occurrence_text)
        if occurrence > len(matches):
            raise LigandReplaceError(
                f"{label} '{selector}' requests occurrence {occurrence}, but only "
                f"{len(matches)} {resname} residue(s) exist."
            )
        selected = matches[occurrence - 1]
    return selected.atoms


def _apply_itp_atom_names(
    universe: "mda.Universe", itp_path: str, resname: str
) -> None:
    """Apply the matching ITP's atom names to a coordinate file in atom order.

    Some MOL2 producers write only element symbols in the atom-name column
    (``C``, ``N``, ...). The ITP is the authoritative final atom ordering for
    this workflow and restores names such as ``C1`` and ``N2`` before fitting.
    """
    try:
        definitions = parse_itp(itp_path)
    except ITPParseError as exc:
        raise LigandReplaceError(str(exc)) from exc

    candidates = [
        definition
        for definition in definitions.values()
        if definition.atom_count == len(universe.atoms)
        and (
            not resname
            or resname in definition.resnames
            or definition.name == resname
        )
    ]
    if len(candidates) != 1:
        candidate_names = [definition.name for definition in candidates]
        raise LigandReplaceError(
            f"new_ligand.itp_path '{itp_path}' must contain exactly one "
            f"moleculetype matching resname '{resname}' and the coordinate "
            f"atom count ({len(universe.atoms)}); matched {candidate_names or 'none'}."
        )

    definition = candidates[0]
    current_hydrogen_flags = [
        is_hydrogen_name(str(name)) for name in universe.atoms.names
    ]
    itp_hydrogen_flags = [
        is_hydrogen_name(str(name)) for name in definition.atom_names
    ]
    if current_hydrogen_flags != itp_hydrogen_flags:
        mismatches = [
            index
            for index, (current, expected) in enumerate(
                zip(current_hydrogen_flags, itp_hydrogen_flags), start=1
            )
            if current != expected
        ]
        raise LigandReplaceError(
            f"new_ligand.coord_path and new_ligand.itp_path disagree on whether "
            f"atom(s) {mismatches[:10]} are hydrogen. ITP atom names can only be "
            "applied when both files use the same atom order."
        )
    universe.atoms.names = list(definition.atom_names)


def load_new_ligand(
    coord_path: str,
    format: str,
    resname_override: str,
    itp_path: str = "",
) -> "mda.Universe":
    fmt = None if format in (None, "", "auto") else format.lower()
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            u = mda.Universe(coord_path, format=fmt) if fmt else mda.Universe(coord_path)
    except Exception as exc:
        raise LigandReplaceError(
            f"Cannot read new_ligand.coord_path '{coord_path}'"
            + (f" as format '{format}'" if fmt else "")
            + f": {exc}"
        ) from exc
    if len(u.atoms) == 0:
        raise LigandReplaceError(f"new_ligand.coord_path '{coord_path}' contains zero atoms.")
    if resname_override:
        u.residues.resnames = [resname_override] * len(u.residues)
    if itp_path:
        _apply_itp_atom_names(u, itp_path, resname_override)
    return u


def _resolve_fit_coords(ag, atom_names: List[str], label: str) -> np.ndarray:
    """Return an (N,3) coordinate array, one row per name in atom_names, in
    that exact order -- so old/new fit arrays line up positionally.
    Fail-stop on a missing or ambiguous (duplicate) atom name."""
    coords = []
    names_present = list(ag.names)
    for name in atom_names:
        idx = [i for i, n in enumerate(names_present) if n == name]
        if not idx:
            raise LigandReplaceError(
                f"{label}: atom name '{name}' not found. Atom names present: {sorted(set(names_present))}"
            )
        if len(idx) > 1:
            raise LigandReplaceError(
                f"{label}: atom name '{name}' matches {len(idx)} atoms (ambiguous) -- "
                f"fit atom names must be unique within the ligand."
            )
        coords.append(ag.positions[idx[0]])
    return np.asarray(coords, dtype=np.float64)


def _unique_heavy_atom_indices(ag, label: str) -> dict[str, int]:
    """Return ``{atom_name: atom_index}`` for an unambiguous heavy-atom set."""
    heavy_names = [
        str(name) for name in ag.names if not is_hydrogen_name(str(name))
    ]
    duplicates = sorted({name for name in heavy_names if heavy_names.count(name) > 1})
    if duplicates:
        raise LigandReplaceError(
            f"{label} has duplicate heavy-atom name(s) {duplicates}; autofit pairs "
            "heavy atoms by exact unique name, so the correspondence is ambiguous. "
            "Use fit.method: pairfit with explicit atom lists instead."
        )
    if not heavy_names:
        raise LigandReplaceError(f"{label} contains no heavy atoms for autofit.")
    return {
        str(name): index
        for index, name in enumerate(ag.names)
        if not is_hydrogen_name(str(name))
    }


def _resolve_autofit_coords(old_ag, new_ag) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """Match all heavy atoms by name and return old/new coordinates in one order."""
    old_by_name = _unique_heavy_atom_indices(old_ag, "original ligand")
    new_by_name = _unique_heavy_atom_indices(new_ag, "new ligand")
    old_names = set(old_by_name)
    new_names = set(new_by_name)
    if old_names != new_names:
        missing_from_new = sorted(old_names - new_names)
        only_in_new = sorted(new_names - old_names)
        raise LigandReplaceError(
            "autofit requires the original and new ligands to have identical, "
            "uniquely named heavy-atom sets so every fitted atom has an unambiguous "
            f"correspondence. Missing from new ligand: {missing_from_new or 'none'}; "
            f"only in new ligand: {only_in_new or 'none'}. Use fit.method: pairfit "
            "to define a shared core explicitly, or fit.method: nofit for an already "
            "positioned ligand."
        )

    # Preserve the original ligand's heavy-atom order, but look up the incoming
    # coordinates by name so the coordinate files need not use the same order.
    ordered_names = list(old_by_name)
    old_coords = np.asarray(
        [old_ag.positions[old_by_name[name]] for name in ordered_names],
        dtype=np.float64,
    )
    new_coords = np.asarray(
        [new_ag.positions[new_by_name[name]] for name in ordered_names],
        dtype=np.float64,
    )
    return old_coords, new_coords, ordered_names


def _require_non_collinear(coords: np.ndarray, label: str) -> None:
    if len(coords) < 3:
        raise LigandReplaceError(
            f"{label} requires at least three atom pairs to determine a unique "
            "three-dimensional rigid transform."
        )
    centered = coords - coords.mean(axis=0)
    if np.linalg.matrix_rank(centered, tol=1e-8) < 2:
        raise LigandReplaceError(
            f"{label} fit atoms are collinear; choose at least three non-collinear atoms."
        )


@dataclass
class LigandReplaceResult:
    merged_universe: "mda.Universe"
    environment_ag_size: int
    inserted_ag_size: int
    n_original_ligand_atoms_removed: int
    n_new_ligand_atoms_inserted: int
    old_ligand_resname: str
    new_ligand_resname: str
    fit_method: str
    n_fit_atoms: int
    fit_rmsd: Optional[float]
    fit_atom_pairs: List[Tuple[str, str]] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)
    # Kept for an optional before/after inspection overlay (old ligand pose
    # vs. the new ligand's final, already-positioned pose) -- both AtomGroups
    # are independent of merged_universe (the old ligand isn't part of it at
    # all, since it was removed) so they have to be carried separately.
    old_ligand_ag: Optional["mda.core.groups.AtomGroup"] = None
    new_ligand_ag_positioned: Optional["mda.core.groups.AtomGroup"] = None


def assemble_ligand_replacement(
    structure_universe: "mda.Universe",
    ligand_replace_cfg: LigandReplaceSpec,
) -> LigandReplaceResult:
    old_selector = ligand_replace_cfg.original_ligand.resname
    old_ligand_ag = select_by_resname(structure_universe, old_selector, "ligand_replace.original_ligand.resname")
    old_resname = str(old_ligand_ag.residues[0].resname)

    all_indices = np.arange(len(structure_universe.atoms))
    old_ligand_indices = set(old_ligand_ag.indices.tolist())
    environment_mask = ~np.isin(all_indices, list(old_ligand_indices))
    environment_ag = structure_universe.atoms[environment_mask]

    new_u = load_new_ligand(
        ligand_replace_cfg.new_ligand.coord_path,
        ligand_replace_cfg.new_ligand.format,
        ligand_replace_cfg.new_ligand.resname,
        ligand_replace_cfg.new_ligand.itp_path,
    )
    new_ligand_ag = new_u.atoms
    new_resname = ligand_replace_cfg.new_ligand.resname or (
        new_ligand_ag.resnames[0] if len(new_ligand_ag.resnames) else old_resname
    )

    notes: List[str] = []
    fit = ligand_replace_cfg.fit
    fit_rmsd: Optional[float] = None
    fit_atom_pairs: List[Tuple[str, str]] = []

    if fit.method == "pairfit":
        old_fit_coords = _resolve_fit_coords(old_ligand_ag, fit.old_ligand_fit_atoms, "ligand_replace.fit.old_ligand_fit_atoms")
        new_fit_coords = _resolve_fit_coords(new_ligand_ag, fit.new_ligand_fit_atoms, "ligand_replace.fit.new_ligand_fit_atoms")
        _require_non_collinear(old_fit_coords, "ligand_replace.fit.old_ligand_fit_atoms")
        _require_non_collinear(new_fit_coords, "ligand_replace.fit.new_ligand_fit_atoms")
        R, mobile_centroid, ref_centroid, fit_rmsd = _kabsch_fit(new_fit_coords, old_fit_coords)
        all_new_coords = new_ligand_ag.positions.astype(np.float64)
        transformed = (R @ (all_new_coords - mobile_centroid).T).T + ref_centroid
        new_ligand_ag.positions = transformed.astype(np.float32)
        fit_atom_pairs = list(zip(fit.old_ligand_fit_atoms, fit.new_ligand_fit_atoms))
        notes.append(
            f"pairfit: superposed {len(fit_atom_pairs)} matched atom pair(s) onto the original "
            f"'{old_resname}' pose, RMSD {fit_rmsd:.3f} A over the fit atoms; the transform was "
            f"then applied to all {len(new_ligand_ag)} atoms of the new ligand."
        )
    elif fit.method == "autofit":
        old_fit_coords, new_fit_coords, fit_names = _resolve_autofit_coords(
            old_ligand_ag, new_ligand_ag
        )
        _require_non_collinear(old_fit_coords, "original ligand autofit heavy atoms")
        _require_non_collinear(new_fit_coords, "new ligand autofit heavy atoms")
        R, mobile_centroid, ref_centroid, fit_rmsd = _kabsch_fit(
            new_fit_coords, old_fit_coords
        )
        all_new_coords = new_ligand_ag.positions.astype(np.float64)
        transformed = (R @ (all_new_coords - mobile_centroid).T).T + ref_centroid
        new_ligand_ag.positions = transformed.astype(np.float32)
        fit_atom_pairs = [(name, name) for name in fit_names]
        notes.append(
            f"autofit: automatically matched all {len(fit_names)} uniquely named heavy "
            f"atoms and superposed the new ligand onto the original '{old_resname}' pose, "
            f"RMSD {fit_rmsd:.3f} A over the heavy atoms; the transform was then applied "
            f"to all {len(new_ligand_ag)} atoms of the new ligand."
        )
    elif fit.method == "nofit":
        notes.append(
            "nofit: the new ligand's coordinates were used exactly as given in "
            f"'{ligand_replace_cfg.new_ligand.coord_path}', with no transformation. Make sure "
            "they are already in the same reference frame as ligand_replace.structure_path "
            "(e.g. docked directly into this structure)."
        )
    else:
        raise LigandReplaceError(f"Unknown ligand_replace.fit.method '{fit.method}'")

    merged = mda.Merge(environment_ag, new_ligand_ag)
    merged.dimensions = structure_universe.dimensions.copy()

    return LigandReplaceResult(
        merged_universe=merged,
        environment_ag_size=len(environment_ag),
        inserted_ag_size=len(new_ligand_ag),
        n_original_ligand_atoms_removed=len(old_ligand_ag),
        n_new_ligand_atoms_inserted=len(new_ligand_ag),
        old_ligand_resname=old_resname,
        new_ligand_resname=new_resname,
        fit_method=fit.method,
        n_fit_atoms=len(fit_atom_pairs),
        fit_rmsd=fit_rmsd,
        fit_atom_pairs=fit_atom_pairs,
        notes=notes,
        old_ligand_ag=old_ligand_ag,
        new_ligand_ag_positioned=new_ligand_ag,
    )
