"""GROMACS index generator for generic and CHARMM-GUI-style groups."""

from __future__ import annotations

import os
from typing import Dict, Tuple

import numpy as np

from classify import classify_resnames
from masks import resolve_mask, MaskError


def _mask_for(resname_class: Dict[str, str], resnames, *classes) -> np.ndarray:
    wanted = set(classes)
    return np.array([resname_class[rn] in wanted for rn in resnames])


def build_groups(universe, style: str = "both") -> Dict[str, np.ndarray]:
    """Return one-based atom indices for the requested index-group style."""
    resnames = universe.residues.resnames
    resname_class = classify_resnames(resnames)

    def atoms_for(*classes) -> np.ndarray:
        res_mask = _mask_for(resname_class, resnames, *classes)
        return universe.residues[res_mask].atoms.ix + 1  # 1-based

    generic = {
        "Protein": atoms_for("protein"),
        "Ligand": atoms_for("other"),
        "Protein_Ligand": atoms_for("protein", "other"),
        "Lipid": atoms_for("lipid"),
        "Water": atoms_for("water"),
        "Ion": atoms_for("ion"),
        "Solvent": atoms_for("water", "ion"),
        "Membrane": atoms_for("protein", "other", "lipid"),
        "System": np.arange(1, universe.atoms.n_atoms + 1),
    }
    charmm_gui = {
        "SOLU": atoms_for("protein", "other"),
        "MEMB": atoms_for("lipid"),
        "SOLV": atoms_for("water", "ion"),
        "SOLU_MEMB": atoms_for("protein", "other", "lipid"),
        "SYSTEM": np.arange(1, universe.atoms.n_atoms + 1),
    }
    if style == "generic":
        return generic
    if style == "charmm_gui":
        return charmm_gui
    if style == "both":
        return {**generic, **charmm_gui}
    raise ValueError(f"Unknown ndx.style '{style}'")


def write_ndx(universe, ndx_cfg) -> Tuple[str, Dict[str, int], list]:
    groups = build_groups(universe, ndx_cfg.style)
    notes = []

    for name, mask_str in (ndx_cfg.extra_groups or {}).items():
        try:
            ag = resolve_mask(universe, mask_str)
        except MaskError as e:
            raise MaskError(f"ndx.extra_groups['{name}']: {e}")
        groups[name] = ag.ix + 1

    out_path = ndx_cfg.output_path
    os.makedirs(os.path.dirname(os.path.abspath(out_path)) or ".", exist_ok=True)

    counts = {}
    with open(out_path, "w") as fh:
        for name, indices in groups.items():
            if len(indices) == 0 and name != "System":
                notes.append(f"Group '{name}' is empty in the final system -- not written to index.ndx.")
                continue
            counts[name] = len(indices)
            fh.write(f"[ {name} ]\n")
            for i in range(0, len(indices), 15):
                fh.write(" ".join(str(int(x)) for x in indices[i:i + 15]) + "\n")
            fh.write("\n")

    return out_path, counts, notes
