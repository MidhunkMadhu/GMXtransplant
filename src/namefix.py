"""Restore nonprotein residue names truncated by a PDB file.

Two fail-safe methods are available. ``atom_signature`` compares residue
atom-name sets with an untruncated reference structure. ``itp_atom_count``
uses a PDB-alias-to-full-name mapping, validates those full names against a
reference topol.top, and distinguishes shared aliases using per-molecule atom
counts from the environment ITP files. Reference topology molecule counts and
order are deliberately ignored. Restoration happens before clash detection,
charge accounting, and final topology assembly.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from typing import Dict, FrozenSet, List, Optional, Sequence, Set, Tuple

import numpy as np

from classify import COMMON_LIPID_RESNAMES
from topology import (
    parse_top_molecules,
    scan_toppar_atom_counts,
    scan_toppar_dir,
)


class NameFixError(RuntimeError):
    pass


@dataclass
class NameFixResult:
    renamed_residue_count: int
    renames_by_pair: Dict[Tuple[str, str], int]   # (old, new) -> count
    unresolved_resnames: Dict[str, int]            # resname -> residue count left unresolved
    candidate_atom_counts: Dict[str, int]          # full resname -> atoms per molecule
    candidate_itp_paths: Dict[str, str] = field(default_factory=dict)
    notes: List[str] = field(default_factory=list)
    method: str = "atom_signature"
    source_path: Optional[str] = None
    reference_lipid_types: List[str] = field(default_factory=list)
    toppar_lipid_types: List[str] = field(default_factory=list)
    reference_lipid_atom_counts: Dict[str, int] = field(default_factory=dict)
    reference_lipid_itp_paths: Dict[str, str] = field(default_factory=dict)


@dataclass
class ITPNameReference:
    valid_full_names: Set[str]
    reference_lipid_types: List[str]
    toppar_lipid_types: List[str]
    atom_counts: Dict[str, int]
    itp_paths: Dict[str, str]
    notes: List[str] = field(default_factory=list)


def build_resname_templates(reference_ag) -> Dict[str, FrozenSet[str]]:
    """One atom-name-set template per unique resname found in reference_ag
    (take the first residue seen of each resname as representative)."""
    templates: Dict[str, FrozenSet[str]] = {}
    for res in reference_ag.residues:
        if res.resname not in templates:
            templates[res.resname] = frozenset(res.atoms.names.tolist())
    return templates


def restore_full_resnames(
    target_universe,
    templates: Dict[str, FrozenSet[str]],
    already_full_resnames: Optional[set] = None,
    min_jaccard: float = 0.95,
) -> NameFixResult:
    """Rename residues in target_universe (in place) whose resname does not
    already appear as a template key, by matching their atom-name set
    against the templates. Exact set equality is preferred; a close
    (>= min_jaccard) and unambiguous match is accepted as a fallback.
    Residues whose resname already matches a template key are left alone
    (nothing to fix). Returns a report; never renames on an ambiguous or
    low-confidence match."""
    already_full_resnames = already_full_resnames or set(templates.keys())

    resnames = target_universe.residues.resnames
    # Use object dtype so assigning a four-character name such as POPC to
    # a PDB-derived three-character array cannot silently truncate it back
    # to POP before MDAnalysis receives the updated values.
    new_resnames = np.asarray([str(name) for name in resnames], dtype=object)

    renamed_count = 0
    renames_by_pair: Dict[Tuple[str, str], int] = {}
    unresolved: Dict[str, int] = {}
    notes: List[str] = []

    # Group residues by (resname, frozenset(atom names)) so we only do the
    # matching work once per distinct signature, not once per residue.
    residues = target_universe.residues
    sig_cache: Dict[Tuple[str, FrozenSet[str]], Optional[str]] = {}

    for i, res in enumerate(residues):
        rn = res.resname
        if rn in already_full_resnames:
            continue  # already a recognized, untruncated name

        atom_names = frozenset(res.atoms.names.tolist())
        key = (rn, atom_names)
        if key not in sig_cache:
            sig_cache[key] = _best_match(atom_names, templates, min_jaccard)
        matched = sig_cache[key]

        if matched is None:
            unresolved[rn] = unresolved.get(rn, 0) + 1
            continue

        new_resnames[i] = matched
        renamed_count += 1
        pair = (rn, matched)
        renames_by_pair[pair] = renames_by_pair.get(pair, 0) + 1

    if renamed_count:
        target_universe.residues.resnames = new_resnames

    if unresolved:
        notes.append(
            f"Could not confidently restore a full resname for: {unresolved} "
            f"(left unchanged -- these will be classified as 'other' downstream "
            f"unless you add them to config classification overrides)."
        )

    return NameFixResult(
        renamed_residue_count=renamed_count,
        renames_by_pair=renames_by_pair,
        unresolved_resnames=unresolved,
        candidate_atom_counts={k: len(v) for k, v in templates.items()},
        notes=notes,
        method="atom_signature",
    )


def build_itp_name_reference(
    reference_topol: str,
    environment_toppar_dir: str,
    pdb_to_full_resname: Dict[str, List[str]],
) -> ITPNameReference:
    """Validate the name map and collect atom counts from environment ITPs."""
    topology_molecules = parse_top_molecules(reference_topol)
    topology_names = {str(name) for name, _count in topology_molecules}
    resnames_by_mtype, _definition_paths = scan_toppar_dir(
        environment_toppar_dir
    )
    counts_by_mtype, count_paths = scan_toppar_atom_counts(
        environment_toppar_dir
    )

    valid_full_names = set(topology_names)
    for mtype in topology_names:
        valid_full_names.update(resnames_by_mtype.get(mtype, set()))

    mapped_full_names = {
        full_name
        for candidates in pdb_to_full_resname.values()
        for full_name in candidates
    }
    extra_mapping_names = sorted(mapped_full_names - valid_full_names)
    if extra_mapping_names:
        raise NameFixError(
            "pdb_to_full_resname contains full residue name(s) that are not "
            "present in reference_topol or represented by one of its molecule "
            f"ITPs: {extra_mapping_names}"
        )

    toppar_full_names = set(resnames_by_mtype)
    for resnames in resnames_by_mtype.values():
        toppar_full_names.update(resnames)

    reference_lipids = set(topology_names) & COMMON_LIPID_RESNAMES
    for mtype in topology_names:
        reference_lipids.update(
            resnames_by_mtype.get(mtype, set()) & COMMON_LIPID_RESNAMES
        )
    toppar_lipids = toppar_full_names & COMMON_LIPID_RESNAMES

    missing_lipid_itps = sorted(reference_lipids - toppar_lipids)
    if missing_lipid_itps:
        raise NameFixError(
            "reference_topol contains lipid type(s) with no matching molecule "
            "definition in topology.environment_toppar_dir: "
            f"{missing_lipid_itps}"
        )

    notes = []
    extra_toppar_lipids = sorted(toppar_lipids - reference_lipids)
    if extra_toppar_lipids:
        notes.append(
            "The environment toppar directory also contains unused common lipid "
            f"ITPs: {extra_toppar_lipids}. They are allowed and are not added to "
            "the output topology unless present in the final coordinates."
        )

    atom_counts: Dict[str, int] = {}
    itp_paths: Dict[str, str] = {}
    names_needing_counts = mapped_full_names | reference_lipids
    for full_name in names_needing_counts:
        matching_mtypes = sorted(
            mtype
            for mtype in topology_names
            if (
                mtype == full_name
                or full_name in resnames_by_mtype.get(mtype, set())
            )
            and mtype in counts_by_mtype
        )
        if len(matching_mtypes) == 1:
            mtype = matching_mtypes[0]
            atom_counts[full_name] = counts_by_mtype[mtype]
            itp_paths[full_name] = count_paths[mtype]
        elif len(matching_mtypes) > 1:
            raise NameFixError(
                f"Full residue name '{full_name}' is defined by multiple molecule "
                f"types permitted by reference_topol: {matching_mtypes}. "
                "The ITP atom-count source is ambiguous."
            )

    missing_lipid_atom_counts = sorted(reference_lipids - set(atom_counts))
    if missing_lipid_atom_counts:
        raise NameFixError(
            "Could not obtain a populated ITP [ atoms ] count for reference "
            f"lipid type(s): {missing_lipid_atom_counts}"
        )

    return ITPNameReference(
        valid_full_names=valid_full_names,
        reference_lipid_types=sorted(reference_lipids),
        toppar_lipid_types=sorted(toppar_lipids),
        atom_counts=atom_counts,
        itp_paths=itp_paths,
        notes=notes,
    )


def itp_atom_count_resname_plan(
    observed_residues: Sequence[Tuple[str, int]],
    pdb_to_full_resname: Dict[str, List[str]],
    valid_full_names: Set[str],
    itp_atom_counts: Dict[str, int],
) -> Tuple[List[str], Dict[Tuple[str, str], int], Dict[str, int]]:
    """Plan PDB alias restoration using exact ITP molecule atom counts."""
    mapping = {
        str(alias).strip(): [str(name).strip() for name in candidates]
        for alias, candidates in pdb_to_full_resname.items()
    }
    if not mapping:
        raise NameFixError("pdb_to_full_resname is empty")
    if any(not alias or not candidates for alias, candidates in mapping.items()):
        raise NameFixError(
            "pdb_to_full_resname cannot contain an empty alias or candidate list"
        )
    if any(not name for candidates in mapping.values() for name in candidates):
        raise NameFixError("pdb_to_full_resname cannot contain empty full names")
    if any(len(candidates) != len(set(candidates)) for candidates in mapping.values()):
        raise NameFixError("pdb_to_full_resname candidate lists cannot contain duplicates")

    mapped_full_names = {
        name for candidates in mapping.values() for name in candidates
    }
    extra_mapping_names = sorted(mapped_full_names - set(valid_full_names))
    if extra_mapping_names:
        raise NameFixError(
            "pdb_to_full_resname contains full residue name(s) not supported by "
            f"reference_topol: {extra_mapping_names}"
        )

    planned = [str(name) for name, _atom_count in observed_residues]
    renames_by_pair: Counter = Counter()
    used_atom_counts: Dict[str, int] = {}
    for residue_index, (current_name_raw, observed_atom_count) in enumerate(observed_residues):
        current_name = str(current_name_raw)
        if current_name in mapped_full_names:
            continue
        candidates = mapping.get(current_name)
        if not candidates:
            continue
        if len(candidates) == 1:
            full_name = candidates[0]
        else:
            missing_counts = sorted(set(candidates) - set(itp_atom_counts))
            if missing_counts:
                raise NameFixError(
                    f"Cannot resolve shared PDB alias '{current_name}': no populated "
                    f"[ atoms ] definition was found for {missing_counts} in "
                    "topology.environment_toppar_dir"
                )
            for name in candidates:
                used_atom_counts[name] = itp_atom_counts[name]
            matching_names = [
                name
                for name in candidates
                if itp_atom_counts[name] == observed_atom_count
            ]
            if len(matching_names) != 1:
                expected = ", ".join(
                    f"{name}={itp_atom_counts[name]} atoms" for name in candidates
                )
                raise NameFixError(
                    f"Cannot uniquely restore target environment residue "
                    f"{residue_index + 1} named '{current_name}' with "
                    f"{observed_atom_count} atoms. Candidate ITP counts: {expected}."
                )
            full_name = matching_names[0]
        planned[residue_index] = full_name
        if current_name != full_name:
            renames_by_pair[(current_name, full_name)] += 1

    return planned, dict(renames_by_pair), dict(sorted(used_atom_counts.items()))


def restore_full_resnames_from_itp_atom_counts(
    target_environment_ag,
    pdb_to_full_resname: Dict[str, List[str]],
    reference: ITPNameReference,
    reference_topol: str,
) -> NameFixResult:
    """Restore mapped target-environment names using ITP atom counts."""
    residues = target_environment_ag.residues
    observed_residues = [
        (str(res.resname), len(res.atoms)) for res in residues
    ]
    planned, renames_by_pair, atom_counts_used = itp_atom_count_resname_plan(
        observed_residues,
        pdb_to_full_resname,
        reference.valid_full_names,
        reference.atom_counts,
    )
    residues.resnames = np.asarray(planned, dtype=object)

    shared_aliases = sorted(
        alias
        for alias, candidates in pdb_to_full_resname.items()
        if len(candidates) > 1
    )
    notes = [
        "Restored mapped nonprotein residue names before clash detection, "
        "charge accounting, and topology assembly.",
        "Reference topol.top molecule counts and order were intentionally ignored.",
    ]
    notes.extend(reference.notes)
    if shared_aliases:
        notes.append(
            f"Shared PDB alias(es) {shared_aliases} were disambiguated using "
            "exact per-molecule atom counts from environment ITP files."
        )

    return NameFixResult(
        renamed_residue_count=sum(renames_by_pair.values()),
        renames_by_pair=renames_by_pair,
        unresolved_resnames={},
        candidate_atom_counts=atom_counts_used,
        candidate_itp_paths={
            name: reference.itp_paths[name] for name in atom_counts_used
        },
        notes=notes,
        method="itp_atom_count",
        source_path=reference_topol,
        reference_lipid_types=reference.reference_lipid_types,
        toppar_lipid_types=reference.toppar_lipid_types,
        reference_lipid_atom_counts={
            name: reference.atom_counts[name]
            for name in reference.reference_lipid_types
        },
        reference_lipid_itp_paths={
            name: reference.itp_paths[name]
            for name in reference.reference_lipid_types
        },
    )


def _best_match(atom_names: FrozenSet[str], templates: Dict[str, FrozenSet[str]], min_jaccard: float) -> Optional[str]:
    exact = [rn for rn, tmpl in templates.items() if tmpl == atom_names]
    if len(exact) == 1:
        return exact[0]
    if len(exact) > 1:
        return None  # ambiguous even with exact match (shouldn't happen with real data)

    scored = []
    for rn, tmpl in templates.items():
        union = len(atom_names | tmpl)
        if union == 0:
            continue
        jac = len(atom_names & tmpl) / union
        if jac >= min_jaccard:
            scored.append((jac, rn))
    if not scored:
        return None
    scored.sort(reverse=True)
    if len(scored) > 1 and scored[0][0] == scored[1][0]:
        return None  # tie, ambiguous
    return scored[0][1]
