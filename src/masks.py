"""
cpptraj-style residue/atom mask parsing.

cpptraj numbers *every* residue in a topology sequentially starting at 1,
in file order, regardless of the PDB/GRO ``resid`` field. This matters a
lot for CHARMM-GUI-style systems (like the ones this workflow targets),
where the resid field commonly *resets to 1* at the start of every
segment (protein chain A, chain B, ligand, lipid, ion, water, ...). Using
the raw resid field to define "residues 1-963" would therefore silently
grab the wrong atoms. This module always resolves masks against the
sequential (file-order) residue index, exactly like cpptraj's default
numbering.

Supported syntax (a practical subset of full cpptraj masks, sufficient
for the residue-range + atom-name selections this workflow needs):

    ":1-963"              residues 1 through 963 (inclusive), sequential
    ":1,5,10-20"          comma-separated list of residues/ranges
    ":1-963@CA"           residues 1-963, atom name CA only
    ":1-963@CA,C,N,O"     residues 1-963, backbone atom names
    ":1-963&!:962-963"    residues 1-963 excluding 962-963 (set difference)

A bare atom-name clause with no leading ":residue" part (e.g. "@CA") is
also accepted and applies to the whole universe.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Optional

import numpy as np

try:
    import MDAnalysis as mda
    from MDAnalysis.core.groups import AtomGroup
except ImportError:  # pragma: no cover
    mda = None
    AtomGroup = None


class MaskError(ValueError):
    pass


@dataclass
class ParsedMask:
    raw: str
    residue_indices: np.ndarray  # 0-based sequential residue indices (cpptraj_res - 1)
    atom_names: Optional[list]   # None means "all atoms in these residues"
    negate: bool = False


def _parse_residue_ranges(spec: str) -> np.ndarray:
    """Parse '1-963,970,975-980' (1-based, inclusive) into a sorted
    0-based numpy index array."""
    spec = spec.strip()
    if not spec:
        return np.array([], dtype=int)
    indices = []
    seen = set()
    tokens = spec.split(",")
    if any(not token.strip() for token in tokens):
        raise MaskError(f"Residue selection '{spec}' contains an empty comma-separated token")
    for token in tokens:
        token = token.strip()
        m = re.fullmatch(r"(\d+)\s*-\s*(\d+)", token)
        if m:
            lo, hi = int(m.group(1)), int(m.group(2))
            if lo < 1 or hi < lo:
                raise MaskError(f"Invalid residue range '{token}' in mask")
            expanded = list(range(lo - 1, hi))  # convert to 0-based
        elif re.fullmatch(r"\d+", token):
            v = int(token)
            if v < 1:
                raise MaskError(f"Invalid residue index '{token}' in mask")
            expanded = [v - 1]
        else:
            raise MaskError(f"Could not parse residue token '{token}' in mask '{spec}'")
        repeated = sorted(seen.intersection(expanded))
        if repeated:
            display = ", ".join(str(value + 1) for value in repeated[:10])
            raise MaskError(
                f"Residue selection '{spec}' repeats residue position(s) {display}; "
                "remove duplicate or overlapping terms"
            )
        seen.update(expanded)
        indices.extend(expanded)
    return np.array(sorted(indices), dtype=int)


def parse_mask(mask_str: str) -> ParsedMask:
    """Parse a single cpptraj-style mask clause (no '&' combination)."""
    original = mask_str
    mask_str = mask_str.strip()
    negate = False
    if mask_str.startswith("!"):
        negate = True
        mask_str = mask_str[1:].strip()

    residue_part = None
    atom_part = None

    if mask_str.count("@") > 1:
        raise MaskError(
            f"Mask clause '{original}' contains more than one atom-name clause; "
            "only one @name1,name2 clause is supported"
        )

    if mask_str.startswith(":"):
        rest = mask_str[1:]
        if "@" in rest:
            residue_part, atom_part = rest.split("@", 1)
            if not residue_part.strip():
                raise MaskError(
                    f"Mask clause '{original}' has an empty residue clause; use "
                    "a bare '@name' mask for all residues"
                )
        else:
            residue_part = rest
    elif mask_str.startswith("@"):
        atom_part = mask_str[1:]
    else:
        raise MaskError(
            f"Mask clause '{original}' must start with ':' (residue) or '@' (atom name)"
        )

    if residue_part is not None and not residue_part.strip():
        raise MaskError(f"Mask clause '{original}' has an empty residue selection")
    residue_indices = _parse_residue_ranges(residue_part) if residue_part else np.array([], dtype=int)
    atom_names = None
    if atom_part is not None:
        raw_atom_names = atom_part.split(",")
        if any(not name.strip() for name in raw_atom_names):
            raise MaskError(f"Mask clause '{original}' has an empty atom-name selection")
        atom_names = [a.strip() for a in raw_atom_names]
        duplicates = sorted({name for name in atom_names if atom_names.count(name) > 1})
        if duplicates:
            raise MaskError(
                f"Mask clause '{original}' repeats atom name(s) {duplicates}; "
                "repeated atom clauses are not silently discarded"
            )

    if negate and atom_names is not None:
        raise MaskError(
            f"Negated atom-name clause '{original}' is not implemented. "
            "Use one positive @name1,name2 clause or change the residue mask."
        )

    return ParsedMask(raw=original, residue_indices=residue_indices, atom_names=atom_names, negate=negate)


def resolve_mask(universe: "mda.Universe", mask_str: str) -> "AtomGroup":
    """Resolve a cpptraj-style mask string into an MDAnalysis AtomGroup,
    using sequential (file-order) residue numbering starting at 1.

    Supports '&' to combine clauses (set intersection when both are
    residue-only; here we mainly support the common pattern of
    ':range' optionally combined with '&!:sub-range' to exclude a
    sub-range, or ':range@names' to restrict to specific atom names).
    """
    if mda is None:
        raise ImportError("MDAnalysis is required for mask resolution")

    raw_clauses = mask_str.split("&")
    if any(not clause.strip() for clause in raw_clauses):
        raise MaskError(f"Mask '{mask_str}' contains an empty '&' clause")
    clauses = [c.strip() for c in raw_clauses]
    if not clauses:
        raise MaskError(f"Empty mask string: '{mask_str}'")

    n_res = len(universe.residues)
    selected_res_mask = None  # boolean array over residues
    atom_name_filter = None
    atom_clause = None

    for clause in clauses:
        pm = parse_mask(clause)
        if pm.residue_indices.size or (not clause.lstrip("!").startswith("@")):
            bad = pm.residue_indices[(pm.residue_indices < 0) | (pm.residue_indices >= n_res)]
            if bad.size:
                raise MaskError(
                    f"Mask '{mask_str}' references residue(s) beyond the system "
                    f"(system has {n_res} residues, cpptraj numbering 1-{n_res})"
                )
            clause_mask = np.zeros(n_res, dtype=bool)
            clause_mask[pm.residue_indices] = True
            if pm.negate:
                clause_mask = ~clause_mask
            selected_res_mask = clause_mask if selected_res_mask is None else (selected_res_mask & clause_mask)
        if pm.atom_names is not None:
            if atom_name_filter is not None:
                raise MaskError(
                    f"Mask '{mask_str}' contains multiple atom-name clauses "
                    f"('{atom_clause}' and '{clause}'); their combination is not implemented"
                )
            atom_name_filter = pm.atom_names
            atom_clause = clause

    if selected_res_mask is None:
        selected_res_mask = np.ones(n_res, dtype=bool)

    residue_index_list = np.nonzero(selected_res_mask)[0]
    if residue_index_list.size == 0:
        raise MaskError(f"Mask '{mask_str}' resolved to zero residues")

    residues = universe.residues[residue_index_list]
    ag = residues.atoms
    if atom_name_filter:
        ag = ag[np.isin(ag.names, atom_name_filter)]
        if len(ag) == 0:
            raise MaskError(
                f"Mask '{mask_str}' selected residues but atom name(s) "
                f"{atom_name_filter} matched zero atoms"
            )
    return ag


def atom_names_in_mask(mask_str: str) -> Optional[list]:
    """Return the one explicit atom-name clause in a validated mask."""
    raw_clauses = mask_str.split("&")
    if any(not clause.strip() for clause in raw_clauses):
        raise MaskError(f"Mask '{mask_str}' contains an empty '&' clause")
    found = []
    source = []
    for clause in raw_clauses:
        parsed = parse_mask(clause.strip())
        if parsed.atom_names is not None:
            found.append(parsed.atom_names)
            source.append(clause.strip())
    if len(found) > 1:
        raise MaskError(
            f"Mask '{mask_str}' contains multiple atom-name clauses {source}; "
            "their combination is not implemented"
        )
    return found[0] if found else None


def residue_cpptraj_numbers(atomgroup: "AtomGroup", universe: "mda.Universe") -> np.ndarray:
    """Return the 1-based cpptraj sequential residue numbers for every
    residue represented in atomgroup (unique, sorted)."""
    all_resindices = universe.residues.resindices
    resindex_to_seq = {ri: i + 1 for i, ri in enumerate(all_resindices)}
    unique_resix = np.unique(atomgroup.resindices)
    return np.array(sorted(resindex_to_seq[ri] for ri in unique_resix))
