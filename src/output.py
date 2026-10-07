"""
Final system assembly and PDB/GRO output.
"""

from __future__ import annotations

import os
import tempfile
import warnings
from typing import List

import numpy as np

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import MDAnalysis as mda

from chainids import apply_chain_ids


class OutputError(RuntimeError):
    pass


def _promote_staged_outputs(staged) -> None:
    """Promote a set of staged files and restore every old destination on error."""
    backups = []
    promoted = []
    try:
        for _temporary, destination in staged:
            if os.path.isdir(destination) and not os.path.islink(destination):
                raise OutputError(
                    f"Coordinate output destination is a directory: '{destination}'"
                )
            backup = None
            if os.path.lexists(destination):
                fd, backup = tempfile.mkstemp(
                    prefix=f".{os.path.basename(destination)}.backup.",
                    dir=os.path.dirname(destination),
                )
                os.close(fd)
                os.unlink(backup)
                os.replace(destination, backup)
            backups.append((destination, backup))

        for temporary, destination in staged:
            os.replace(temporary, destination)
            promoted.append(destination)
    except (OSError, OutputError) as exc:
        rollback_errors = []
        for destination in reversed(promoted):
            try:
                os.unlink(destination)
            except FileNotFoundError:
                pass
            except OSError as rollback_exc:
                rollback_errors.append(str(rollback_exc))
        for destination, backup in reversed(backups):
            if backup is None:
                continue
            try:
                os.replace(backup, destination)
            except OSError as rollback_exc:
                rollback_errors.append(str(rollback_exc))
        if rollback_errors:
            raise OutputError(
                "Coordinate promotion failed and prior outputs could not be fully "
                f"restored: promotion={exc}; rollback={rollback_errors}"
            ) from exc
        raise OutputError(
            f"Coordinate promotion failed; prior outputs were restored: {exc}"
        ) from exc
    else:
        for _destination, backup in backups:
            if backup is not None:
                try:
                    os.unlink(backup)
                except OSError:
                    # The validated outputs are already in place. A stale hidden
                    # backup is safer than reporting a failed coordinate write.
                    pass

# Canonical whole-system ordering, enforced unconditionally on every run --
# not a config option. classify_resnames()'s "other" bucket (anything not
# recognized as protein/water/ion/lipid -- i.e. every ligand/cofactor) is
# folded into the "ligand" bucket here, matching the exact same convention
# topology.py's own determine_final_molecule_sequence() already uses to
# read this order back out when building [ molecules ]. Keep these two in
# sync if either ever changes.
_CANONICAL_CLASS_BUCKET = {"protein": "protein", "ligand": "ligand", "other": "ligand", "lipid": "lipid", "ion": "ion", "water": "water"}


def _reorder_canonical(universe, protein_chain_atom_counts=None):
    """Return a new Universe with residues reordered into the fixed
    protein -> ligand -> lipid -> ion -> water sequence. Residues within
    each class keep their original relative file order (e.g. multiple
    protein chains stay in their existing chain order, multiple lipid
    species keep whatever order they were already in) -- only the order
    *between* classes changes.

    Note: this deliberately does NOT do `mda.Merge(universe.residues[perm].atoms)`
    with a single already-permuted AtomGroup. mda.Merge() rebuilds each
    ResidueAttr (resname, resid, ...) from `ag.residues`, and
    AtomGroup.residues always comes back in ascending resindex order
    regardless of the order atoms/residues were selected in -- so a single
    combined AtomGroup silently gets its residue TABLE reset to the
    original order even though the returned Universe's raw atom/coordinate
    array is correctly permuted. That mismatch would leave
    step5_input.pdb/gro correctly reordered while topol.top (built
    from final_universe.residues, via topology.py's
    determine_final_molecule_sequence) and chain-ID assignment (chainids.py,
    also residues-table based) silently stayed in the OLD order. Passing
    one AtomGroup per class as separate positional args avoids this: each
    per-class group's own residues are already in ascending order (no
    within-class reordering needed), so Merge's residues-table rebuild is
    a no-op for each one, and concatenating the args in canonical class
    order is what actually sets the final between-class order correctly
    for atoms, residues, and every attribute alike."""
    from classify import classify_resnames

    class_ags = []
    if protein_chain_atom_counts is not None:
        protein_atom_count = sum(int(count) for count in protein_chain_atom_counts)
        if protein_atom_count > len(universe.atoms):
            raise RuntimeError(
                f"Topology describes {protein_atom_count} protein atoms, but the "
                f"assembled system contains only {len(universe.atoms)} atoms."
            )
        if protein_atom_count:
            class_ags.append(universe.atoms[:protein_atom_count])
        remaining_residues = universe.atoms[protein_atom_count:].residues
        if len(remaining_residues):
            rn_class = classify_resnames(set(remaining_residues.resnames))
            resolved_class = np.array([
                _CANONICAL_CLASS_BUCKET[rn_class[rn]]
                for rn in remaining_residues.resnames
            ])
        else:
            resolved_class = np.array([], dtype=str)
        class_order = ("ligand", "lipid", "ion", "water")
        residue_source = remaining_residues
    else:
        class_order = ("protein", "ligand", "lipid", "ion", "water")
        residue_source = universe.residues
        rn_class = classify_resnames(set(residue_source.resnames))
        resolved_class = np.array([
            _CANONICAL_CLASS_BUCKET[rn_class[rn]] for rn in residue_source.resnames
        ])

    for cls in class_order:
        mask = resolved_class == cls
        if mask.any():
            class_ags.append(residue_source[mask].atoms)

    if not class_ags:
        return universe

    reordered = mda.Merge(*class_ags)
    return reordered


def canonicalize_system(universe, box_dimensions=None, protein_chain_atom_counts=None):
    """Return a canonically ordered copy with the requested box dimensions."""
    reordered = _reorder_canonical(universe, protein_chain_atom_counts)
    dims = box_dimensions if box_dimensions is not None else universe.dimensions
    if dims is not None:
        reordered.dimensions = np.asarray(dims, dtype=np.float32).copy()
    return reordered


def write_pdb(universe, path: str, protein_chain_atom_counts=None):
    """Write one inspection or intermediate PDB with meaningful chain IDs."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    apply_chain_ids(universe, protein_chain_atom_counts)
    universe.atoms.write(path)


def finalize_system(
    kept_environment_ag,
    inserted_ag,
    box_dimensions,
    ions_removed: List[dict],
    protein_chain_atom_counts=None,
    inserted_contains_protein: bool = True,
):
    """Remove neutralization-selected ions from kept_environment_ag, merge
    with inserted_ag, and reorder the result into the fixed canonical
    sequence (protein, then ligand, then lipid, then ion, then water) that
    every output -- step5_input.pdb/gro AND topol.top's [ molecules ]
    section, since that's derived by walking this same final_universe in
    file order -- must always follow."""
    if ions_removed:
        remove_resindices = {r["resindex"] for r in ions_removed}
        keep_mask = ~np.isin(kept_environment_ag.resindices, list(remove_resindices))
        final_env_ag = kept_environment_ag[keep_mask]
    else:
        final_env_ag = kept_environment_ag

    if inserted_contains_protein:
        final_universe = mda.Merge(inserted_ag, final_env_ag)
    else:
        final_universe = mda.Merge(final_env_ag, inserted_ag)
    return canonicalize_system(
        final_universe, box_dimensions, protein_chain_atom_counts
    )


def write_outputs(
    final_universe, pdb_path: str, gro_path: str,
    protein_chain_atom_counts=None,
    protected_input_paths=None,
    staged_validator=None,
):
    pdb_abs = os.path.abspath(pdb_path)
    gro_abs = os.path.abspath(gro_path)
    if os.path.realpath(pdb_abs) == os.path.realpath(gro_abs):
        raise OutputError("output.pdb_path and output.gro_path must be different files")
    os.makedirs(os.path.dirname(pdb_abs), exist_ok=True)
    os.makedirs(os.path.dirname(gro_abs), exist_ok=True)
    # GRO has no chain-ID column and doesn't need this; the PDB writer
    # would otherwise stamp a meaningless placeholder chain ('X') on
    # every atom, so assign real per-chain letters first (see chainids.py).
    apply_chain_ids(final_universe, protein_chain_atom_counts)

    protected = {
        os.path.realpath(os.path.abspath(path))
        for path in (protected_input_paths or []) if path
    }
    collisions = sorted(
        path for path in (pdb_abs, gro_abs) if os.path.realpath(path) in protected
    )

    staged = []
    try:
        for destination in (pdb_abs, gro_abs):
            suffix = os.path.splitext(destination)[1]
            fd, temporary = tempfile.mkstemp(
                prefix=f".{os.path.basename(destination)}.",
                suffix=suffix,
                dir=os.path.dirname(destination),
            )
            os.close(fd)
            staged.append((temporary, destination))

        final_universe.atoms.write(staged[0][0])
        final_universe.atoms.write(staged[1][0])

        expected_names = [str(name) for name in final_universe.atoms.names]
        expected_resnames = [str(name) for name in final_universe.atoms.resnames]
        for temporary, destination in staged:
            try:
                check = mda.Universe(temporary)
            except Exception as exc:
                raise OutputError(
                    f"Staged coordinate output for '{destination}' cannot be read back: {exc}"
                ) from exc
            if len(check.atoms) != len(final_universe.atoms):
                raise OutputError(
                    f"Staged coordinate output for '{destination}' has {len(check.atoms)} "
                    f"atoms; expected {len(final_universe.atoms)}"
                )
            if [str(name) for name in check.atoms.names] != expected_names:
                raise OutputError(
                    f"Staged coordinate output for '{destination}' changed atom names or order"
                )
            if [str(name) for name in check.atoms.resnames] != expected_resnames:
                raise OutputError(
                    f"Staged coordinate output for '{destination}' changed residue names or order"
                )
            # Both writers round by format (PDB 0.001 Å, GRO 0.01 Å).
            # Validate units and the cell as well as names before promotion.
            expected_positions = final_universe.atoms.positions.astype(np.float64)
            storage_tolerance = 2 * float(np.max(np.spacing(np.abs(final_universe.atoms.positions))))
            coordinate_tolerance = (0.00051 if destination == pdb_abs else 0.0051) + storage_tolerance
            if not np.allclose(check.atoms.positions, expected_positions,
                               atol=coordinate_tolerance, rtol=0):
                raise OutputError(
                    f"Staged coordinate output for '{destination}' changed coordinates "
                    f"beyond format rounding ({coordinate_tolerance:g} Å)"
                )
            if final_universe.dimensions is not None:
                # PDB lengths: 0.001 Å; angles: 0.01 degrees. GRO stores vectors
                # at 0.0001 Å precision, with angles reconstructed on reading.
                if check.dimensions is None or not np.allclose(
                    check.dimensions, final_universe.dimensions,
                    atol=np.array([0.001, 0.001, 0.001, 0.02, 0.02, 0.02]), rtol=0,
                ):
                    raise OutputError(
                        f"Staged coordinate output for '{destination}' lost or changed "
                        "the cell beyond PDB/GRO precision (Å and degrees)"
                    )

        additional_validation = None
        if staged_validator is not None:
            additional_validation = staged_validator(staged[0][0], staged[1][0])

        # Promotion happens only after both complete files have been read back
        # and any caller-supplied topology/charge audit has succeeded. Existing
        # destinations are backed up on their own filesystems and restored if
        # either rename fails, including when an output is also an input.
        _promote_staged_outputs(staged)
        staged.clear()
    finally:
        for temporary, _destination in staged:
            try:
                os.unlink(temporary)
            except FileNotFoundError:
                pass

    return {
        "pdb_path": pdb_abs,
        "gro_path": gro_abs,
        "atomic_write": True,
        "input_path_collisions_safely_replaced": collisions,
        "staged_validation": additional_validation,
    }


def refgro_completeness_check(final_universe, refgro_path: str) -> dict:
    """Compare the final system's resnames against a reference GRO's own
    non-protein resnames -- every non-protein species refgro contains is
    expected to also appear somewhere in the final output. Informational
    only (never fail-stop): refgro is a shape/completeness reference, not
    necessarily this run's own donor file, so a genuinely apo run or a
    protein-only environment swap can legitimately be missing some of
    refgro's species."""
    from classify import classify_resnames

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ref_u = mda.Universe(refgro_path)
    ref_resnames = sorted(set(ref_u.residues.resnames))
    ref_class = classify_resnames(ref_resnames)
    ref_nonprotein = sorted(rn for rn in ref_resnames if ref_class[rn] != "protein")

    final_resnames = set(final_universe.residues.resnames)
    missing = sorted(rn for rn in ref_nonprotein if rn not in final_resnames)

    return {
        "refgro_path": refgro_path,
        "refgro_nonprotein_resnames": ref_nonprotein,
        "missing_from_final": missing,
    }


def write_inspection_pdb(target_full_ag, aligned_mobile_ag, path: str):
    """Write a combined inspection structure: the box's original protein
    (chain X) overlaid with the newly-aligned replacement protein
    (chain Y), so alignment quality can be checked visually before the
    original protein is removed."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    merged = mda.Merge(target_full_ag, aligned_mobile_ag)
    n1 = len(target_full_ag)
    chain_ids = ["X"] * n1 + ["Y"] * (len(merged.atoms) - n1)
    merged.add_TopologyAttr("chainIDs", chain_ids)
    merged.atoms.write(path)
