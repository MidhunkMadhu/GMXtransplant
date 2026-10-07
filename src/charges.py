"""
Topology-aware net-charge accounting and counterion-based neutralization.

When topology inputs are configured, exact molecule charges are summed from
each ITP ``[ atoms ]`` section, including protein chains, lipids, ligands,
water, and ions. The residue-level table is only a fallback. Lipids are not
assumed neutral because the lipid registry includes anionic species. Any
unresolved final-system residue stops the workflow before neutralization.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import MDAnalysis as mda
    from MDAnalysis.lib.distances import capped_distance, distance_array

from classify import (
    classify_resnames, DEFAULT_WATER_RESNAMES,
    AMINO_ACID_RESNAMES, is_hydrogen_name,
)
from config import ChargeSpec
from ion_selection import select_counterions

# Standard amino-acid formal charges by resname (protonation-state aware).
# Everything not listed here is assumed neutral (0).
AMINO_ACID_CHARGES = {
    "ASP": -1.0, "GLU": -1.0,
    "ASH": 0.0, "GLH": 0.0,          # protonated (neutral) Asp/Glu
    "LYS": 1.0, "ARG": 1.0,
    "LYN": 0.0,                       # neutral Lys
    "HSP": 1.0, "HIP": 1.0,           # doubly-protonated (charged) His
    "HSD": 0.0, "HSE": 0.0, "HIS": 0.0, "HID": 0.0, "HIE": 0.0,
    "CYX": 0.0,
}

ION_CHARGES = {
    "SOD": 1.0, "NA": 1.0, "NA+": 1.0,
    "CLA": -1.0, "CL": -1.0, "CL-": -1.0,
    "POT": 1.0, "K": 1.0,
    "CAL": 2.0, "CA2": 2.0, "CA": 2.0,
    "MG": 2.0, "ZN": 2.0, "ZN2": 2.0,
    "CES": 1.0, "CS": 1.0,
    "LIT": 1.0, "LI": 1.0,
    "RUB": 1.0, "RB": 1.0,
}


class ChargeError(RuntimeError):
    pass


def build_charge_table(overrides: Dict[str, Optional[float]]) -> Dict[str, float]:
    table = {}
    # All standard amino acids default to neutral; AMINO_ACID_CHARGES then
    # overrides the charged/protonation-state-dependent subset.
    for rn in AMINO_ACID_RESNAMES:
        table[rn] = 0.0
    table.update(AMINO_ACID_CHARGES)
    # Lipids are deliberately not assigned a generic zero charge. Neutral
    # and anionic lipid names share the classification registry, so their
    # charges must come from an ITP or an explicit override.
    table.update(ION_CHARGES)
    for wn in DEFAULT_WATER_RESNAMES:
        table[wn] = 0.0
    for rn, val in (overrides or {}).items():
        if val is not None:
            table[rn] = float(val)
    return table


def residue_charge_map(resnames: np.ndarray, charge_table: Dict[str, float]) -> Tuple[np.ndarray, List[str]]:
    """Return (per-atom charge array is NOT what we want -- charges are
    per-RESIDUE for this formal-charge approach; caller should aggregate
    at residue granularity). Here we return, for the unique resnames
    present, which ones are unknown (missing from charge_table and not
    implicitly zero via lipid/other default-neutral assumption)."""
    unique = sorted(set(resnames))
    unknown = [rn for rn in unique if rn not in charge_table]
    return unique, unknown


def charge_of_atomgroup(ag, charge_table: Dict[str, float], default_unknown_is_zero: bool = False) -> Tuple[float, List[str]]:
    """Sum formal charge over all residues in ag. Unknown resnames (not in
    charge_table) are reported; if default_unknown_is_zero is False, an
    unknown resname makes the total charge NaN (unreliable) rather than
    silently assuming zero."""
    total = 0.0
    unknown = set()
    for res in ag.residues:
        rn = res.resname
        if rn in charge_table:
            total += charge_table[rn]
        else:
            unknown.add(rn)
            if not default_unknown_is_zero:
                continue  # will report as NaN below
    if unknown and not default_unknown_is_zero:
        return float("nan"), sorted(unknown)
    return total, sorted(unknown)


def charge_of_atomgroup_from_itp(ag, model, include_protein_prefix: bool = False):
    """Sum exact ITP molecule charges for a complete-molecule AtomGroup."""
    total = 0.0
    unknown = set()
    atom_offset = 0
    if include_protein_prefix:
        layout = model.layout
        if len(ag) < layout.total_atoms:
            raise ChargeError(
                f"Coordinate block has {len(ag)} atoms, fewer than the "
                f"{layout.total_atoms} protein atoms defined by "
                f"{layout.moleculetypes}."
            )
        for mtype, definition in zip(layout.moleculetypes, layout.definitions):
            stop = atom_offset + definition.atom_count
            actual_names = [str(name) for name in ag[atom_offset:stop].names]
            expected_names = list(definition.atom_names)
            if actual_names != expected_names:
                mismatch = next(
                    i for i, (got, expected) in enumerate(zip(actual_names, expected_names))
                    if got != expected
                )
                raise ChargeError(
                    f"Protein moleculetype {mtype} does not match coordinate atom "
                    f"order at molecule atom {mismatch + 1}: coordinates have "
                    f"'{actual_names[mismatch]}', ITP expects '{expected_names[mismatch]}'."
                )
            total += definition.charge
            atom_offset = stop
        if atom_offset < len(ag) and ag[atom_offset - 1].resindex == ag[atom_offset].resindex:
            raise ChargeError(
                "Topology-derived protein boundary falls inside a coordinate residue."
            )

    remaining = ag[atom_offset:]
    for residue in remaining.residues:
        resname = str(residue.resname)
        mtype = model.resname_to_mtype.get(resname)
        if mtype is None:
            unknown.add(resname)
            continue
        definition = model.definitions[mtype]
        if len(residue.atoms) != definition.atom_count:
            raise ChargeError(
                f"Coordinate residue {resname}{residue.resid} has "
                f"{len(residue.atoms)} atoms, but moleculetype '{mtype}' in "
                f"'{model.files[mtype]}' has {definition.atom_count}."
            )
        actual_names = [str(name) for name in residue.atoms.names]
        expected_names = list(definition.atom_names)
        if actual_names != expected_names:
            mismatch = next(
                i for i, (got, expected) in enumerate(zip(actual_names, expected_names))
                if got != expected
            )
            raise ChargeError(
                f"Coordinate residue {resname}{residue.resid} does not match "
                f"moleculetype '{mtype}' atom order at atom {mismatch + 1}: "
                f"coordinates have '{actual_names[mismatch]}', ITP expects "
                f"'{expected_names[mismatch]}'."
            )
        total += definition.charge
    if unknown:
        return float("nan"), sorted(unknown)
    return total, []


@dataclass
class ChargeReport:
    charge_table_used: Dict[str, float]
    unknown_resnames: List[str]
    charge_original_protein: float
    charge_replacement_protein: float
    charge_diff_from_replacement: float
    charge_environment_untouched: float
    net_charge_before_clash_removal: float
    charge_removed_by_clashes: float
    charge_removed_by_class: Dict[str, float]
    charge_removed_by_resname: Dict[str, float]
    net_charge_after_clash_removal: float
    target_net_charge: float
    charge_excess: float
    neutralization_performed: bool
    ions_removed: List[dict]
    final_net_charge: Optional[float]
    ion_selection_summary: Dict[str, object] = field(default_factory=dict)
    stopped_reason: Optional[str] = None
    notes: List[str] = field(default_factory=list)


def compute_pre_neutralization_report(
    removed_original_ag,
    inserted_ag,
    environment_untouched_ag,
    clash_removed_ag,
    charge_cfg: ChargeSpec,
    topology_charge_model=None,
    inserted_has_protein: bool = False,
    environment_has_protein: bool = False,
    original_charge_from_itp: bool = True,
    resolve_original_charge: bool = True,
) -> ChargeReport:
    table = build_charge_table(charge_cfg.charge_table_overrides)
    notes = []

    if topology_charge_model is not None:
        for resname, itp_charge in topology_charge_model.residue_charge_table.items():
            override = charge_cfg.charge_table_overrides.get(resname)
            if override is not None and abs(float(override) - itp_charge) > charge_cfg.tolerance:
                raise ChargeError(
                    f"charge.charge_table_overrides gives {resname} charge "
                    f"{float(override):+.6f}, but its selected ITP moleculetype "
                    f"sums to {itp_charge:+.6f}."
                )
            table[resname] = itp_charge
        charge_repl, unk_repl = charge_of_atomgroup_from_itp(
            inserted_ag, topology_charge_model, inserted_has_protein
        )
        charge_env, unk_env = charge_of_atomgroup_from_itp(
            environment_untouched_ag,
            topology_charge_model,
            environment_has_protein,
        )
        charge_clash_removed, unk_clash = charge_of_atomgroup_from_itp(
            clash_removed_ag, topology_charge_model, False
        )
        if not resolve_original_charge:
            charge_orig, unk_orig = float("nan"), []
            notes.append(
                "The removed original ligand was intentionally excluded from "
                "topology and charge resolution. Its charge is not used for "
                "final-system charge or neutralization."
            )
        elif original_charge_from_itp:
            charge_orig, unk_orig = charge_of_atomgroup_from_itp(
                removed_original_ag, topology_charge_model, False
            )
        else:
            charge_orig, unk_orig = charge_of_atomgroup(removed_original_ag, table)
            notes.append(
                "The removed original block charge is a residue-table diagnostic "
                "because no matching original-block ITP was supplied. It is not "
                "used for final-system charge or neutralization."
            )
        notes.append(
            "Final-system molecule charges were summed from the selected ITP "
            "moleculetype [ atoms ] sections."
        )
    else:
        if resolve_original_charge:
            charge_orig, unk_orig = charge_of_atomgroup(removed_original_ag, table)
        else:
            charge_orig, unk_orig = float("nan"), []
            notes.append(
                "The removed original ligand was intentionally excluded from "
                "charge resolution. Its charge is not used for final-system "
                "charge or neutralization."
            )
        charge_repl, unk_repl = charge_of_atomgroup(inserted_ag, table)
        charge_env, unk_env = charge_of_atomgroup(environment_untouched_ag, table)
        charge_clash_removed, unk_clash = charge_of_atomgroup(clash_removed_ag, table)

    # Unknown charges in the removed original block do not affect the final
    # system. Only final or clash-removed components can stop neutralization.
    unknown_all = sorted(set(unk_repl) | set(unk_env) | set(unk_clash))
    if unk_orig:
        notes.append(f"Removed original block has unresolved diagnostic charges: {unk_orig}")

    # Per-class breakdown of charge removed by clash detection (for report).
    removed_class = classify_resnames(clash_removed_ag.resnames)
    charge_removed_by_class: Dict[str, float] = {}
    charge_removed_by_resname: Dict[str, float] = {}
    for res in clash_removed_ag.residues:
        cls = removed_class[res.resname]
        c = table.get(res.resname, 0.0)
        charge_removed_by_class[cls] = charge_removed_by_class.get(cls, 0.0) + c
        charge_removed_by_resname[res.resname] = (
            charge_removed_by_resname.get(res.resname, 0.0) + c
        )

    net_before = charge_env + charge_repl if not (np.isnan(charge_env) or np.isnan(charge_repl)) else float("nan")
    net_after = net_before - charge_clash_removed if not np.isnan(net_before) and not np.isnan(charge_clash_removed) else float("nan")

    return ChargeReport(
        charge_table_used=table,
        unknown_resnames=unknown_all,
        charge_original_protein=charge_orig,
        charge_replacement_protein=charge_repl,
        charge_diff_from_replacement=(charge_repl - charge_orig) if not (np.isnan(charge_repl) or np.isnan(charge_orig)) else float("nan"),
        charge_environment_untouched=charge_env,
        net_charge_before_clash_removal=net_before,
        charge_removed_by_clashes=charge_clash_removed,
        charge_removed_by_class=charge_removed_by_class,
        charge_removed_by_resname=charge_removed_by_resname,
        net_charge_after_clash_removal=net_after,
        target_net_charge=charge_cfg.target_net_charge,
        charge_excess=(net_after - charge_cfg.target_net_charge) if not np.isnan(net_after) else float("nan"),
        neutralization_performed=False,
        ions_removed=[],
        final_net_charge=None,
        notes=notes,
    )


def neutralize(
    kept_environment_ag,
    inserted_ag,
    charge_report: ChargeReport,
    charge_cfg: ChargeSpec,
    universe_for_write,
) -> ChargeReport:
    """Select and remove bulk counterions to bring the system to
    charge_cfg.target_net_charge, respecting exclusion distances from the
    protein, ligand, and lipid heavy atoms. Returns an updated
    ChargeReport; the caller is responsible for actually deleting the
    chosen ions from the working AtomGroup (this function reports the
    residue indices to remove via charge_report.ions_removed, keyed by
    resindex within kept_environment_ag's parent universe)."""
    if charge_report.unknown_resnames:
        charge_report.stopped_reason = (
            "Cannot compute a reliable net charge: the following residue(s) have no known "
            "formal charge and none was supplied via replacement_ligands or "
            f"charge.charge_table_overrides: {charge_report.unknown_resnames}. Add their net "
            "charge and ITP metadata to the config and re-run. "
            "The workflow stops here rather than silently assuming a value."
        )
        return charge_report

    if not charge_cfg.neutralize:
        charge_report.stopped_reason = None
        charge_report.final_net_charge = charge_report.net_charge_after_clash_removal
        return charge_report

    excess = charge_report.charge_excess
    if abs(excess) <= charge_cfg.tolerance:
        charge_report.neutralization_performed = True
        charge_report.final_net_charge = charge_report.net_charge_after_clash_removal
        return charge_report

    resname_class = classify_resnames(kept_environment_ag.resnames)
    ion_mask = np.array([resname_class[rn] == "ion" for rn in kept_environment_ag.resnames])
    ion_ag = kept_environment_ag[ion_mask]
    if len(ion_ag) == 0:
        charge_report.stopped_reason = (
            f"Net charge is off target by {excess:+.3f} e and there are no ions left in the "
            "bulk solvent to remove. The workflow stops rather than adding new ions "
            "automatically -- add ions manually and re-run, or adjust target_net_charge."
        )
        return charge_report

    protein_atom_mask = np.isin(inserted_ag.resnames, list(AMINO_ACID_RESNAMES))
    protein_only_ag = inserted_ag[protein_atom_mask]
    ligand_only_ag = inserted_ag[~protein_atom_mask]
    # Any protein atoms already present in the KEPT environment (this
    # happens in ligand_replace mode, where the protein itself is
    # untouched and therefore lives in the environment rather than the
    # inserted block) get the same protein exclusion-distance protection.
    # In the normal protein-swap flow the environment never contains
    # protein atoms (they were stripped via protein_mask), so this is a
    # no-op there.
    env_protein_mask = np.isin(kept_environment_ag.resnames, list(AMINO_ACID_RESNAMES))
    env_protein_ag = kept_environment_ag[env_protein_mask]
    if len(env_protein_ag) > 0:
        protein_only_ag = (protein_only_ag + env_protein_ag) if len(protein_only_ag) > 0 else env_protein_ag

    ion_pos = ion_ag.positions
    box = universe_for_write.dimensions
    if box is None or len(box) != 6:
        raise ChargeError(
            "PBC-aware ion exclusion requires six valid box dimensions."
        )
    d_to_protein = distance_array(
        ion_pos, protein_only_ag.positions, box=box
    ).min(axis=1) if len(protein_only_ag) else np.full(len(ion_ag), np.inf)
    d_to_ligand = distance_array(
        ion_pos, ligand_only_ag.positions, box=box
    ).min(axis=1) if len(ligand_only_ag) else np.full(len(ion_ag), np.inf)

    # A global lipid min-z/max-z slab is invalid for wrapped, unwrapped,
    # tilted, or displaced membranes. Protect membrane-associated ions using
    # direct minimum-image distance to lipid heavy atoms instead. This is
    # independent of the membrane orientation and remains valid when molecule
    # coordinates lie outside the primary unit cell.
    near_lipid = np.zeros(len(ion_ag), dtype=bool)
    lipid_cutoff = float(charge_cfg.exclusion_distance_from_lipid)
    if charge_cfg.exclude_membrane_interior and lipid_cutoff > 0:
        lipid_mask_env = np.array([
            resname_class[rn] == "lipid" for rn in kept_environment_ag.resnames
        ])
        lipid_ag = kept_environment_ag[lipid_mask_env]
        if len(lipid_ag):
            lipid_heavy = lipid_ag[np.array([
                not is_hydrogen_name(name) for name in lipid_ag.names
            ])]
            if len(lipid_heavy):
                pairs = capped_distance(
                    ion_pos,
                    lipid_heavy.positions,
                    max_cutoff=lipid_cutoff,
                    box=box,
                    return_distances=False,
                )
                if len(pairs):
                    near_lipid[np.unique(pairs[:, 0]).astype(int)] = True
        charge_report.notes.append(
            "Neutralization excluded ions within "
            f"{lipid_cutoff:.3f} A of any lipid heavy atom using minimum-image "
            "distances; no absolute-z membrane slab was used."
        )
    elif not charge_cfg.exclude_membrane_interior:
        charge_report.notes.append(
            "Lipid-proximity exclusion disabled; protein and ligand distance "
            "protections remain active."
        )

    # Expand every atom-level exclusion to its complete residue. This matters
    # for any multi-atom ion model and guarantees molecule-level selection.
    def _expand_to_residues(atom_mask):
        excluded = np.unique(ion_ag.resindices[np.asarray(atom_mask, dtype=bool)])
        return np.isin(ion_ag.resindices, excluded)

    near_protein = _expand_to_residues(
        d_to_protein < charge_cfg.exclusion_distance_from_protein
    )
    near_ligand = _expand_to_residues(
        d_to_ligand < charge_cfg.exclusion_distance_from_ligand
    )
    near_lipid = _expand_to_residues(near_lipid)
    eligible_mask = ~(near_protein | near_ligand | near_lipid)
    eligible_ion_ag = ion_ag[eligible_mask]

    table = charge_report.charge_table_used
    needed_ion_sign = 1.0 if excess > 0 else -1.0
    candidates = []
    for res in eligible_ion_ag.residues:
        c = table.get(res.resname, 0.0)
        if c * needed_ion_sign > 0:
            candidates.append((res.resindex, res.resname, c))

    def _counts(atomgroup):
        counts: Dict[str, int] = {}
        for resname in atomgroup.residues.resnames:
            counts[str(resname)] = counts.get(str(resname), 0) + 1
        return counts

    charge_report.ion_selection_summary = {
        "remaining_ions_total": len(ion_ag.residues),
        "remaining_ions_by_resname": _counts(ion_ag),
        "excluded_near_protein": int(len(np.unique(ion_ag.resindices[near_protein]))),
        "excluded_near_ligand": int(len(np.unique(ion_ag.resindices[near_ligand]))),
        "excluded_near_lipid": int(len(np.unique(ion_ag.resindices[near_lipid]))),
        "eligible_ions_total": len(eligible_ion_ag.residues),
        "eligible_ions_by_resname": _counts(eligible_ion_ag),
        "needed_sign_candidates": len(candidates),
    }

    removed, remaining = select_counterions(
        candidates,
        excess_charge=excess,
        tolerance=charge_cfg.tolerance,
        random_seed=charge_cfg.random_seed,
    )

    if abs(remaining) > charge_cfg.tolerance:
        charge_report.stopped_reason = (
            f"Could not fully neutralize the system by removing eligible bulk ions: "
            f"{abs(remaining):.3f} e of net charge remains uncorrected after removing "
            f"{len(removed)} ion(s) from {len(candidates)} eligible ion(s) of the needed sign. "
            f"The available discrete charges could not reach the target within "
            f"tolerance={charge_cfg.tolerance}. Review ion_selection_summary in the report; "
            "adjust the geometric exclusions, add suitable ions, or change the target charge."
        )
        charge_report.ions_removed = removed
        return charge_report

    charge_report.neutralization_performed = True
    charge_report.ions_removed = removed
    charge_report.final_net_charge = charge_report.net_charge_after_clash_removal - sum(r["charge"] for r in removed)
    return charge_report
