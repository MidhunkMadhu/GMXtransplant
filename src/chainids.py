"""
Meaningful chain-ID assignment for the final PDB output.

Neither of the input coordinate files carries real chain identifiers
(GRO has none at all; the source PDB's chain column is blank), so
without this step MDAnalysis's PDB writer would fall back to stamping
every single atom with the same placeholder chain 'X' -- technically
harmless (residue/atom order is what actually matters for re-simulation,
and the GRO output is unaffected since GRO has no chain column), but
unhelpful for visual inspection in PyMOL/VMD. When a protein topology is
available, protein chains come from the ordered protein moleculetypes and
their exact ITP atom counts. Each ligand residue gets its own chain, all
lipids share one membrane chain, and ions and water share one solvent chain.
A residue-gap heuristic remains only as a compatibility fallback when no
protein topology was configured.
"""

from __future__ import annotations

import string
from typing import Optional

import numpy as np

from classify import classify_resnames, AMINO_ACID_RESNAMES

_ALPHABET = list(string.ascii_uppercase) + list(string.ascii_lowercase) + list(string.digits)


def assign_chain_ids(residue_resnames: np.ndarray, residue_resids: np.ndarray) -> np.ndarray:
    """Given per-residue resname/resid arrays (in file order), return a
    per-residue chain-letter array. This is the compatibility fallback used
    only without topology-derived protein boundaries."""
    n = len(residue_resnames)
    resname_class = classify_resnames(residue_resnames)
    letters = np.empty(n, dtype="<U1")

    shared_letter = {}  # membrane or solvent group -> letter
    next_letter_idx = 0

    def take_letter():
        nonlocal next_letter_idx
        if next_letter_idx >= len(_ALPHABET):
            # Extremely unlikely (>62 distinct chain groups); wrap rather
            # than crash, chain ID becomes cosmetic only at that point.
            next_letter_idx = 0
        letter = _ALPHABET[next_letter_idx]
        next_letter_idx += 1
        return letter

    prev_class = None
    prev_resid = None
    current_letter = None

    for i in range(n):
        rn = residue_resnames[i]
        rid = residue_resids[i]
        cls = resname_class[rn]

        if cls in ("lipid", "water", "ion"):
            group = "lipid" if cls == "lipid" else "solvent"
            if group not in shared_letter:
                shared_letter[group] = take_letter()
            letters[i] = shared_letter[group]
            prev_class = cls
            prev_resid = rid
            continue

        if cls == "protein":
            new_chain = (
                prev_class != "protein"
                or current_letter is None
                or (prev_resid is not None and rid != prev_resid + 1 and rid != prev_resid)
            )
            if new_chain:
                current_letter = take_letter()
            letters[i] = current_letter
            prev_class = "protein"
            prev_resid = rid
            continue

        # 'other' (ligand/cofactor): each residue is one ligand chain.
        current_letter = take_letter()
        letters[i] = current_letter
        prev_class = "other"
        prev_resid = rid

    return letters


def _topology_atom_chain_ids(universe, protein_chain_atom_counts) -> np.ndarray:
    """Assign protein boundaries from ITP atom counts, never resid gaps."""
    atom_letters = np.empty(len(universe.atoms), dtype="<U1")
    next_letter_idx = 0

    def take_letter():
        nonlocal next_letter_idx
        if next_letter_idx >= len(_ALPHABET):
            next_letter_idx = 0
        letter = _ALPHABET[next_letter_idx]
        next_letter_idx += 1
        return letter

    offset = 0
    for chain_index, count in enumerate(protein_chain_atom_counts, 1):
        count = int(count)
        if count <= 0 or offset + count > len(universe.atoms):
            raise RuntimeError(
                f"Invalid topology-derived atom count {count} for protein chain "
                f"{chain_index} at atom offset {offset}."
            )
        stop = offset + count
        if stop < len(universe.atoms) and (
            universe.atoms[stop - 1].resindex == universe.atoms[stop].resindex
        ):
            left = universe.atoms[stop - 1]
            right = universe.atoms[stop]
            raise RuntimeError(
                f"Topology-derived protein chain {chain_index} boundary at atom "
                f"{stop} splits coordinate residue {left.resname}{left.resid}; "
                f"boundary atoms are {left.name}/{right.name}."
            )
        atom_letters[offset:stop] = take_letter()
        offset = stop

    shared = {}
    remaining_residues = universe.atoms[offset:].residues
    classes = classify_resnames(remaining_residues.resnames)
    for residue in remaining_residues:
        cls = classes[residue.resname]
        if cls == "protein":
            raise RuntimeError(
                "Protein residue outside topology-derived protein chains: "
                f"{residue.resname}{residue.resid} at residue index "
                f"{residue.resindex + 1}."
            )
        if cls in ("lipid", "ion", "water"):
            # All lipid types share one membrane chain. Ions and waters
            # together share one solvent chain, as requested for final PDBs.
            group = "lipid" if cls == "lipid" else "solvent"
            if group not in shared:
                shared[group] = take_letter()
            letter = shared[group]
        else:
            # Each ligand residue is its own chain.
            letter = take_letter()
        atom_letters[residue.atoms.indices] = letter
    return atom_letters


def apply_chain_ids(universe, protein_chain_atom_counts=None) -> None:
    """Compute and set a per-atom chainIDs TopologyAttr on `universe`
    in place, based on its residues in file order."""
    if protein_chain_atom_counts is not None:
        atom_letters = _topology_atom_chain_ids(
            universe, protein_chain_atom_counts
        )
    else:
        res_letters = assign_chain_ids(universe.residues.resnames, universe.residues.resids)
        counts = np.array([r.atoms.n_atoms for r in universe.residues])
        atom_letters = np.repeat(res_letters, counts)
    if universe.atoms.n_atoms != len(atom_letters):
        raise RuntimeError("Chain-ID assignment size mismatch; refusing to write inconsistent chain IDs.")
    if hasattr(universe, "add_TopologyAttr"):
        try:
            universe.del_TopologyAttr("chainIDs")
        except Exception:
            pass
        universe.add_TopologyAttr("chainIDs", list(atom_letters))
