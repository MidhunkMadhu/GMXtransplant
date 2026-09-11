"""
PyMOL-based automatic structural alignment.

Uses PyMOL's built-in ``align`` (sequence alignment via Needleman-Wunsch +
iterative outlier-rejection superposition) or ``cealign`` (structure-based,
sequence-independent) to compute the rotation/translation that best
superposes the replacement receptor block onto the receptor already
present in the original simulation box, and applies that transform to
every atom of the replacement block (not just the atoms used for the fit).
"""

from __future__ import annotations

import os
import tempfile
import warnings
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import MDAnalysis as mda

from masks import atom_names_in_mask, resolve_mask, residue_cpptraj_numbers
from config import AlignmentSpec
from kabsch import kabsch_fit


class AlignmentError(RuntimeError):
    pass


@dataclass
class AlignmentResult:
    rmsd_before: float
    rmsd_after: float
    n_atom_pairs_before: int
    n_atom_pairs_after: int
    n_cycles_executed: int
    fit_atom_selection: str
    mobile_block_global_resnums: List[int]
    target_block_global_resnums: List[int]
    excluded_as_outliers: int
    new_positions: np.ndarray            # (n_atoms_in_mobile_block, 3), same order as mobile_ag
    rotation: np.ndarray                 # 3x3
    translation: np.ndarray              # 3,
    notes: List[str] = field(default_factory=list)


def _compress_ranges(nums: List[int]) -> str:
    """Compress a sorted list of ints into 'a-b,c,d-e' form."""
    if not nums:
        return ""
    nums = sorted(nums)
    ranges = []
    start = prev = nums[0]
    for n in nums[1:]:
        if n == prev + 1:
            prev = n
            continue
        ranges.append((start, prev))
        start = prev = n
    ranges.append((start, prev))
    parts = [f"{a}-{b}" if a != b else f"{a}" for a, b in ranges]
    return ",".join(parts)


def _legacy_atom_names(atoms_mode: Optional[str]) -> Optional[List[str]]:
    if atoms_mode in (None, "", "all"):
        return None
    if atoms_mode == "CA":
        return ["CA"]
    if atoms_mode == "backbone":
        return ["C", "CA", "N", "O"]
    raise AlignmentError(f"Unknown legacy alignment.atoms mode '{atoms_mode}'")


def _atom_filter_clause(atom_names: Optional[List[str]]) -> str:
    if not atom_names:
        return ""
    return " and name " + "+".join(atom_names)


def _filter_by_atom_names(ag, atom_names: Optional[List[str]]):
    """Order-preserving fallback for the legacy alignment.atoms option."""
    if not atom_names:
        return ag
    return ag[np.isin(ag.names, atom_names)]


def _local_residue_ordinals(ag) -> np.ndarray:
    order = {}
    result = []
    for resindex in ag.resindices:
        key = int(resindex)
        if key not in order:
            order[key] = len(order)
        result.append(order[key])
    return np.asarray(result, dtype=int)


def _write_block_pdb(universe: "mda.Universe", block_ag, global_nums: List[int], path: str, chain: str = "A"):
    """Write block_ag to a PDB with a single chain and sequential resSeq
    numbers matching position-in-block (1..N), so region sub-selections
    (expressed in the *original* global cpptraj numbering) can be mapped
    to local PyMOL 'resi' selections unambiguously.

    global_nums must be the sorted list of global cpptraj residue numbers
    corresponding, in order, to block_ag.residues (i.e. the output of
    residue_cpptraj_numbers(block_ag, universe)).
    """
    merged = mda.Merge(block_ag)
    # merged.residues is in the same order as block_ag.residues, so the
    # local resSeq for residue i (0-based) is simply i+1, and global_nums[i]
    # is its corresponding global cpptraj number.
    merged.residues.resids = np.arange(1, len(merged.residues) + 1)
    merged.add_TopologyAttr("chainIDs", [chain] * len(merged.atoms))
    merged.atoms.positions = block_ag.positions
    merged.atoms.write(path)


def align_replacement_to_box(
    original_universe: "mda.Universe",
    replacement_universe: "mda.Universe",
    original_receptor_mask: str,
    replacement_receptor_mask: str,
    alignment: AlignmentSpec,
    workdir: str,
) -> AlignmentResult:
    os.makedirs(workdir, exist_ok=True)

    target_full_ag = resolve_mask(original_universe, original_receptor_mask)
    mobile_full_ag = resolve_mask(replacement_universe, replacement_receptor_mask)

    target_global_nums = residue_cpptraj_numbers(target_full_ag, original_universe).tolist()
    mobile_global_nums = residue_cpptraj_numbers(mobile_full_ag, replacement_universe).tolist()

    region_orig_mask = alignment.region_mask_original or original_receptor_mask
    region_repl_mask = alignment.region_mask_replacement or replacement_receptor_mask

    region_orig_ag = resolve_mask(original_universe, region_orig_mask)
    region_repl_ag = resolve_mask(replacement_universe, region_repl_mask)
    target_atom_names = atom_names_in_mask(region_orig_mask)
    mobile_atom_names = atom_names_in_mask(region_repl_mask)
    legacy_atom_names = _legacy_atom_names(alignment.atoms)
    if target_atom_names is None and legacy_atom_names:
        region_orig_ag = _filter_by_atom_names(region_orig_ag, legacy_atom_names)
        target_atom_names = legacy_atom_names
        notes_legacy_target = True
    else:
        notes_legacy_target = False
    if mobile_atom_names is None and legacy_atom_names:
        region_repl_ag = _filter_by_atom_names(region_repl_ag, legacy_atom_names)
        mobile_atom_names = legacy_atom_names
        notes_legacy_mobile = True
    else:
        notes_legacy_mobile = False
    region_orig_global = set(residue_cpptraj_numbers(region_orig_ag, original_universe).tolist())
    region_repl_global = set(residue_cpptraj_numbers(region_repl_ag, replacement_universe).tolist())

    notes: List[str] = []
    if notes_legacy_target or notes_legacy_mobile:
        notes.append(
            "Deprecated alignment.atoms was used because one or both region masks "
            "did not contain an @atom-name clause. Put the atom selection directly "
            "in both alignment region masks instead."
        )
    missing_orig = region_orig_global - set(target_global_nums)
    missing_repl = region_repl_global - set(mobile_global_nums)
    if missing_orig:
        notes.append(
            f"alignment.region_mask_original selects {len(missing_orig)} residue(s) outside "
            f"target_box.receptor_mask; they will be ignored for the fit."
        )
    if missing_repl:
        notes.append(
            f"alignment.region_mask_replacement selects {len(missing_repl)} residue(s) outside "
            f"replacement_structure.receptor_mask; they will be ignored for the fit."
        )
    region_orig_global &= set(target_global_nums)
    region_repl_global &= set(mobile_global_nums)
    if not region_orig_global or not region_repl_global:
        raise AlignmentError("Alignment region resolved to zero residues on one side; check masks.")

    # Atoms outside the receptor_mask block are dropped from the fit
    # selection (same policy the pymol_* methods below apply via their own
    # local-resi intersection; the notes above already warn if this
    # happened).
    region_orig_ag = region_orig_ag[np.isin(region_orig_ag.resindices, target_full_ag.residues.resindices)]
    region_repl_ag = region_repl_ag[np.isin(region_repl_ag.resindices, mobile_full_ag.residues.resindices)]

    if target_atom_names is not None and mobile_atom_names is not None:
        if target_atom_names != mobile_atom_names:
            raise AlignmentError(
                "Alignment fit masks specify different ordered atom-name lists: "
                f"target {target_atom_names} vs replacement {mobile_atom_names}."
            )
    for label, atomgroup in (
        ("target", region_orig_ag), ("replacement", region_repl_ag),
    ):
        if len(atomgroup) < 3:
            raise AlignmentError(
                f"Alignment {label} fit selection contains {len(atomgroup)} point(s); "
                "at least three non-collinear fit points are required."
            )
        centered = atomgroup.positions.astype(np.float64)
        centered -= centered.mean(axis=0)
        if np.linalg.matrix_rank(centered, tol=1e-8) < 2:
            raise AlignmentError(
                f"Alignment {label} fit selection is collinear; choose at least "
                "three non-collinear fit points."
            )

    if alignment.method == "mask_fit":
        # cpptraj-'rms'-style direct fit: two masks, matched positionally
        # (mask order == file order here), no sequence alignment, no
        # iterative outlier rejection. A plain Kabsch/SVD best fit is
        # computed over the atom selections encoded directly in the masks and
        # applied rigidly to every atom of the WHOLE replacement
        # receptor_mask block -- same "fit on a subset, move the whole
        # block" behavior as the pymol_* methods, just without PyMOL.
        fit_target_ag = region_orig_ag
        fit_mobile_ag = region_repl_ag

        if len(fit_target_ag) != len(fit_mobile_ag):
            raise AlignmentError(
                f"alignment.method 'mask_fit' requires the two fit selections to have exactly "
                f"the same number of atoms -- they are matched positionally (cpptraj-style: no "
                f"sequence alignment, no outlier rejection), so a count mismatch means the "
                f"pairing would be wrong even if it didn't error. Got {len(fit_target_ag)} "
                f"atom(s) from target_box (region_mask_original, or target_box.receptor_mask if "
                "unset) vs "
                f"{len(fit_mobile_ag)} atom(s) from replacement_structure (region_mask_replacement, "
                "or replacement_structure.receptor_mask if unset). Adjust the masks "
                f"so both sides select the same residues, in the same order."
            )
        if len(fit_target_ag) == 0:
            raise AlignmentError("alignment.method 'mask_fit': fit selection is empty; check your masks.")

        target_names = np.asarray(fit_target_ag.names, dtype=str)
        mobile_names = np.asarray(fit_mobile_ag.names, dtype=str)
        if not np.array_equal(target_names, mobile_names):
            mismatch = int(np.flatnonzero(target_names != mobile_names)[0])
            raise AlignmentError(
                "alignment.method 'mask_fit' requires element-wise atom-name "
                f"correspondence. Pair {mismatch + 1} is target "
                f"'{target_names[mismatch]}' vs replacement "
                f"'{mobile_names[mismatch]}'. Equal atom counts alone are not safe."
            )
        target_ordinals = _local_residue_ordinals(fit_target_ag)
        mobile_ordinals = _local_residue_ordinals(fit_mobile_ag)
        if not np.array_equal(target_ordinals, mobile_ordinals):
            mismatch = int(np.flatnonzero(target_ordinals != mobile_ordinals)[0])
            raise AlignmentError(
                "alignment.method 'mask_fit' has different per-residue atom "
                f"grouping at pair {mismatch + 1}: target local residue "
                f"{target_ordinals[mismatch] + 1}, replacement local residue "
                f"{mobile_ordinals[mismatch] + 1}."
            )
        target_resnames = [str(res.resname) for res in fit_target_ag.residues]
        mobile_resnames = [str(res.resname) for res in fit_mobile_ag.residues]
        resname_mismatches = [
            (i + 1, left, right)
            for i, (left, right) in enumerate(zip(target_resnames, mobile_resnames))
            if left != right
        ]
        if len(target_resnames) != len(mobile_resnames):
            raise AlignmentError(
                "alignment.method 'mask_fit' selected different numbers of residues: "
                f"{len(target_resnames)} target vs {len(mobile_resnames)} replacement."
            )
        if resname_mismatches:
            notes.append(
                "mask_fit residue-name differences were retained as intentional "
                f"variants after atom correspondence validation: {resname_mismatches}"
            )

        ref_coords = fit_target_ag.positions.astype(np.float64)
        mobile_coords = fit_mobile_ag.positions.astype(np.float64)
        try:
            R, mobile_centroid, ref_centroid, rmsd_after = kabsch_fit(
                mobile_coords, ref_coords
            )
        except ValueError as exc:
            raise AlignmentError(str(exc)) from exc

        old_full = mobile_full_ag.positions.astype(np.float64)
        new_full = (R @ (old_full - mobile_centroid).T).T + ref_centroid

        notes.append(
            f"mask_fit: direct Kabsch/SVD best fit over {len(fit_target_ag)} positionally-matched "
            f"atom pair(s) selected directly by the region masks, RMSD {rmsd_after:.3f} A; no "
            f"sequence alignment, no outlier rejection (equivalent to cpptraj's "
            f"'rms <ref_mask> <mask>' given two explicit masks). Applied rigidly to all "
            f"{len(mobile_full_ag)} atoms of the replacement receptor_mask block."
        )
        return AlignmentResult(
            rmsd_before=float("nan"),  # not meaningful here: the two files' raw coordinate
                                        # frames are independent before any fit, unlike PyMOL's
                                        # own before/after outlier-rejection comparison.
            rmsd_after=rmsd_after,
            n_atom_pairs_before=len(fit_target_ag),
            n_atom_pairs_after=len(fit_target_ag),
            n_cycles_executed=0,
            fit_atom_selection=(
                f"mask_fit: target '{region_orig_mask}' ({len(fit_target_ag)} atoms) "
                f"vs replacement '{region_repl_mask}' ({len(fit_mobile_ag)} atoms), matched positionally"
            ),
            mobile_block_global_resnums=mobile_global_nums,
            target_block_global_resnums=target_global_nums,
            excluded_as_outliers=0,
            new_positions=new_full,
            rotation=R,
            translation=ref_centroid - R @ mobile_centroid,
            notes=notes,
        )

    target_pdb = os.path.join(workdir, "_align_target.pdb")
    mobile_pdb = os.path.join(workdir, "_align_mobile.pdb")
    _write_block_pdb(original_universe, target_full_ag, target_global_nums, target_pdb, chain="A")
    _write_block_pdb(replacement_universe, mobile_full_ag, mobile_global_nums, mobile_pdb, chain="A")

    target_local_to_global = {i + 1: g for i, g in enumerate(target_global_nums)}
    mobile_local_to_global = {i + 1: g for i, g in enumerate(mobile_global_nums)}
    target_global_to_local = {g: i for i, g in target_local_to_global.items()}
    mobile_global_to_local = {g: i for i, g in mobile_local_to_global.items()}

    target_region_local = sorted(target_global_to_local[g] for g in region_orig_global)
    mobile_region_local = sorted(mobile_global_to_local[g] for g in region_repl_global)

    target_atom_clause = _atom_filter_clause(target_atom_names)
    mobile_atom_clause = _atom_filter_clause(mobile_atom_names)
    target_resi_clause = f" and resi {_compress_ranges(target_region_local)}"
    mobile_resi_clause = f" and resi {_compress_ranges(mobile_region_local)}"

    try:
        import pymol2
    except ImportError as exc:
        raise AlignmentError(
            f"alignment.method '{alignment.method}' requires the pymol2 module. "
            "Install PyMOL or use alignment.method: mask_fit."
        ) from exc

    p = pymol2.PyMOL()
    p.start()
    cmd = p.cmd
    try:
        cmd.load(mobile_pdb, "mobile")
        cmd.load(target_pdb, "target")

        n_mobile_pre = cmd.count_atoms("mobile")
        n_target_pre = cmd.count_atoms("target")

        mobile_sel = f"mobile{mobile_resi_clause}{mobile_atom_clause}"
        target_sel = f"target{target_resi_clause}{target_atom_clause}"

        n_mobile_fit = cmd.count_atoms(mobile_sel)
        n_target_fit = cmd.count_atoms(target_sel)
        if n_mobile_fit == 0 or n_target_fit == 0:
            raise AlignmentError(
                f"Alignment fit selection is empty (mobile: {n_mobile_fit} atoms, "
                f"target: {n_target_fit} atoms). Check the alignment region masks."
            )

        if alignment.method == "pymol_align":
            result = cmd.align(
                mobile_sel, target_sel,
                cycles=alignment.cycles, cutoff=alignment.cutoff,
                object="aln",
            )
            # PyMOL align() returns:
            # (RMSD_after, n_atoms_after, n_cycles, RMSD_before,
            #  n_atoms_before, raw_score, n_residues_aligned)
            rmsd_after, n_after, n_cycles, rmsd_before, n_before, raw_score, n_res_aligned = result
        elif alignment.method == "pymol_cealign":
            result = cmd.cealign(target_sel, mobile_sel)
            rmsd_after = result["RMSD"]
            n_after = result["alignment_length"]
            rmsd_before = float("nan")
            n_before = n_after
            n_cycles = 0
        else:
            raise AlignmentError(f"Unknown alignment method '{alignment.method}'")

        if int(n_after) < 3:
            raise AlignmentError(
                f"PyMOL retained only {int(n_after)} aligned atom pair(s); at least "
                "three non-collinear pairs are required."
            )

        # IMPORTANT: PyMOL does NOT preserve file order internally -- on
        # load it reorders atoms within each residue by its own canonical
        # priority (verified empirically: N/CA/C/O backbone first, then
        # sidechain, independent of the PDB file's atom order). So we
        # cannot assume cmd.get_coords("mobile")[i] corresponds to
        # mobile_full_ag.positions[i]. What IS preserved is each atom's
        # PDB serial number (exposed as Atom.id in get_model()), and the
        # PDB we wrote via _write_block_pdb has serials 1..N in exactly
        # mobile_full_ag's order. We use that as an explicit, verified
        # mapping back to our original atom order rather than relying on
        # positional correspondence.
        model_after = cmd.get_model("mobile")
        if len(model_after.atom) != n_mobile_pre:
            raise AlignmentError("Could not retrieve transformed mobile coordinates from PyMOL.")
        if len(mobile_full_ag) != n_mobile_pre:
            raise AlignmentError(
                f"Atom-count mismatch between MDAnalysis extraction ({len(mobile_full_ag)}) "
                f"and PyMOL-loaded object ({n_mobile_pre}); refusing to proceed."
            )
        coord_by_serial = {a.id: a.coord for a in model_after.atom}
        if set(coord_by_serial.keys()) != set(range(1, n_mobile_pre + 1)):
            raise AlignmentError(
                "PyMOL atom serial numbers for 'mobile' are not a clean 1..N range; "
                "cannot safely map coordinates back to the original atom order."
            )
        new = np.array([coord_by_serial[i + 1] for i in range(n_mobile_pre)], dtype=np.float64)

        # Recover rotation+translation via Kabsch between original mobile
        # positions and the new ones, now that ordering is verified correct
        # (this is exact since PyMOL applies one rigid-body transform to
        # the whole object; residual is a floating-point-level sanity check).
        old = mobile_full_ag.positions.astype(np.float64)
        old_c = old.mean(axis=0)
        new_c = new.mean(axis=0)
        H = (old - old_c).T @ (new - new_c)
        U, S, Vt = np.linalg.svd(H)
        d = np.sign(np.linalg.det(Vt.T @ U.T))
        D = np.diag([1, 1, d])
        R = Vt.T @ D @ U.T
        t = new_c - R @ old_c
        reconstructed = (R @ old.T).T + t
        max_err = np.max(np.linalg.norm(reconstructed - new, axis=1))
        if max_err > 1e-2:
            notes.append(
                f"WARNING: reconstructed rigid-body transform has max per-atom error "
                f"{max_err:.4f} Angstrom vs PyMOL's actual output; the workflow uses "
                f"PyMOL's own coordinates directly for the insertion (this warning is "
                f"informational only and would indicate a non-rigid discrepancy)."
            )

        excluded = n_mobile_fit - n_after if isinstance(n_after, int) else 0

        return AlignmentResult(
            rmsd_before=float(rmsd_before) if rmsd_before == rmsd_before else float("nan"),
            rmsd_after=float(rmsd_after),
            n_atom_pairs_before=int(n_before),
            n_atom_pairs_after=int(n_after),
            n_cycles_executed=int(n_cycles) if isinstance(n_cycles, (int, float)) else 0,
            fit_atom_selection=(
                f"mobile '{region_repl_mask}' vs target '{region_orig_mask}'"
            ),
            mobile_block_global_resnums=mobile_global_nums,
            target_block_global_resnums=target_global_nums,
            excluded_as_outliers=int(excluded),
            new_positions=new,  # in mobile_full_ag order
            rotation=R,
            translation=t,
            notes=notes,
        )
    finally:
        p.stop()
