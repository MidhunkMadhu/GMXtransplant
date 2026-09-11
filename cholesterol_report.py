"""Text and JSON reports for cholesterol-restoration mode."""

from __future__ import annotations

import json
import os
from datetime import datetime, timezone
from box_validation import render_preflight


def build_report(
    cfg,
    result,
    final_universe,
    target_universe,
    topology_result=None,
    topology_error=None,
    ndx_counts=None,
    ndx_notes=None,
    ndx_error=None,
    refgro_check=None,
):
    clash = result.clash_result
    comp = result.composition
    align = result.alignment
    return {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "box_validation": getattr(target_universe, "box_preflight", None),
        "mode": "chl",
        "inputs": {
            "experimental_structure": cfg.cholesterol.experimental_structure_path,
            "converted_structure": result.converted_structure_path,
            "target_system": cfg.cholesterol.target_system_path,
            "composition_reference_system": (
                cfg.cholesterol.composition.reference_system_path
                or cfg.cholesterol.target_system_path
            ),
            "target_receptor_mask": cfg.cholesterol.target_receptor_mask,
            "experimental_receptor_mask": cfg.cholesterol.experimental_receptor_mask,
        },
        "conversion": {
            "enabled": cfg.cholesterol.convert_to_charmm36,
            "diagnostics_retained": cfg.cholesterol.write_diagnostics,
            "charmm_resname": cfg.cholesterol.charmm_resname,
            "n_cholesterol_residues": result.n_experimental_cholesterols,
            "n_cholesterol_atoms": result.n_experimental_cholesterol_atoms,
            "diagnostic_paths": result.diagnostic_paths,
            "heavy_atom_audit": result.heavy_atom_audit,
        },
        "alignment": {
            "method": cfg.alignment.method,
            "region_mask_original": cfg.alignment.region_mask_original,
            "region_mask_replacement": cfg.alignment.region_mask_replacement,
            "fit_selection": align.fit_atom_selection,
            "n_atom_pairs_before": align.n_atom_pairs_before,
            "n_atom_pairs_after": align.n_atom_pairs_after,
            "rmsd_before_angstrom": align.rmsd_before,
            "rmsd_after_angstrom": align.rmsd_after,
            "excluded_as_outliers": align.excluded_as_outliers,
            "notes": align.notes,
        },
        "clash_detection": {
            "thresholds_angstrom": clash.threshold_used,
            "use_pbc": cfg.clash_detection.use_pbc,
            "heavy_atoms_only": cfg.clash_detection.heavy_atoms_only,
            "counts_by_class": clash.counts_by_class,
            "n_removed": len(clash.removed_molecules),
            "removed_molecules": [
                {
                    "resname": x.resname,
                    "resid": x.resid,
                    "class": x.molclass,
                    "min_distance_angstrom": x.min_distance,
                    "cutoff_angstrom": x.cutoff_used,
                    "contact_env_atom": x.contact_env_atom,
                    "contact_cholesterol_resname": x.contact_partner_resname,
                    "contact_cholesterol_resid": x.contact_partner_resid,
                    "contact_cholesterol_atom": x.contact_partner_atom,
                    "contact_context": x.contact_context,
                    "removal_decision": x.removal_decision,
                }
                for x in clash.removed_molecules
            ],
            "flagged_but_kept": [
                {
                    "resname": x.resname,
                    "resid": x.resid,
                    "class": x.molclass,
                    "min_distance_angstrom": x.min_distance,
                    "cutoff_angstrom": x.cutoff_used,
                    "contact_env_atom": x.contact_env_atom,
                    "contact_cholesterol_resname": x.contact_partner_resname,
                    "contact_cholesterol_resid": x.contact_partner_resid,
                    "contact_cholesterol_atom": x.contact_partner_atom,
                    "contact_context": x.contact_context,
                    "removal_decision": x.removal_decision,
                }
                for x in clash.flagged_but_kept
            ],
            "notes": clash.notes,
        },
        "composition": {
            "midplane_z_angstrom": comp.midplane_z,
            "distance_from_protein_angstrom": (
                cfg.cholesterol.composition.distance_from_protein
            ),
            "random_seed": cfg.cholesterol.composition.random_seed,
            "lipid_targets": comp.lipid_targets,
            "lipid_counts_before": comp.lipid_counts_before,
            "lipid_counts_after": comp.lipid_counts_after,
            "lipid_residues_removed": comp.lipid_residues_removed,
            "box_volume_liters": comp.box_volume_liters,
            "assumed_reference_concentration_molar": comp.target_concentration_molar,
            # Legacy alias retained for consumers of existing JSON reports.
            "target_concentration_molar": comp.target_concentration_molar,
            "concentration_tolerance_fraction": comp.concentration_tolerance_fraction,
            "reference_salt_pairs": comp.reference_salt_pairs,
            "final_salt_pairs": comp.final_salt_pairs,
            "estimated_final_concentration_molar": (
                comp.estimated_final_concentration_molar
            ),
            "concentration_deviation_fraction": (
                comp.concentration_deviation_fraction
            ),
            "concentration_within_tolerance": comp.concentration_within_tolerance,
            "reference_ion_counts": comp.reference_ion_counts,
            "neutralizing_excess": comp.neutralizing_excess,
            "ion_counts_before": comp.ion_counts_before,
            "ion_target_counts": comp.ion_target_counts,
            "ion_counts_after": comp.ion_counts_after,
            "ions_removed": comp.ions_removed,
            "output_ion_counts": comp.output_ion_counts,
            "warnings": comp.warnings,
        },
        "final_system": {
            "n_atoms": len(final_universe.atoms),
            "n_residues": len(final_universe.residues),
            "box_dimensions": list(map(float, final_universe.dimensions)),
            "target_box_dimensions": list(map(float, target_universe.dimensions)),
            "box_unchanged": bool(
                (final_universe.dimensions == target_universe.dimensions).all()
            ),
        },
        "topology": None if topology_result is None and topology_error is None else {
            "top_path": None if topology_result is None else topology_result.top_path,
            "toppar_dir": None if topology_result is None else topology_result.toppar_dir,
            "molecules_written": None if topology_result is None else topology_result.molecules_written,
            "includes_written": None if topology_result is None else topology_result.includes_written,
            "charge_audit": None if topology_result is None else topology_result.charge_audit,
            "notes": [] if topology_result is None else topology_result.notes,
            "error": topology_error,
        },
        "ndx": None if ndx_counts is None and ndx_error is None else {
            "counts": ndx_counts,
            "notes": ndx_notes or [],
            "error": ndx_error,
        },
        "refgro_check": refgro_check,
        "notes": result.notes,
    }


def render_text(data) -> str:
    lines = []
    add = lines.append
    add("=" * 78)
    add("CHOLESTEROL RESTORATION WORKFLOW: VALIDATION REPORT")
    add(f"Generated: {data['generated_at_utc']}")
    add("=" * 78)
    inp = data["inputs"]
    add("\nINPUTS")
    add(f"Experimental structure: {inp['experimental_structure']}")
    converted_label = (
        inp["converted_structure"] or "temporary conversion (not retained)"
    )
    add(f"Converted structure:    {converted_label}")
    add(f"Target system:          {inp['target_system']}")
    add(f"Composition reference:  {inp['composition_reference_system']}")

    conversion = data["conversion"]
    add("\n1. CHOLESTEROL PREPARATION")
    add(f"CHARMM36 conversion enabled: {conversion['enabled']}")
    add(f"Diagnostic PDBs retained:     {conversion['diagnostics_retained']}")
    add(f"Inserted cholesterol residues: {conversion['n_cholesterol_residues']}")
    add(f"Inserted cholesterol atoms:    {conversion['n_cholesterol_atoms']}")
    for item in conversion["heavy_atom_audit"]:
        modeled = item["modeled_heavy_atom_names"]
        add(
            f"  {item['resname']}{item['resid']}: "
            f"{item['observed_heavy_atom_count']}/"
            f"{item['expected_heavy_atom_count']} observed heavy atoms; "
            f"status={item['status']}"
            + (f"; modeled={modeled}" if modeled else "")
            + (
                f"; whole-fit RMSD="
                f"{item['reference_fit_rmsd_angstrom']:.3f} A; "
                f"local-fit RMSDs="
                f"{item['modeled_atom_local_fit_rmsd_angstrom']}"
                if modeled else ""
            )
        )

    alignment = data["alignment"]
    add("\n2. RECEPTOR-BASED ALIGNMENT")
    add(f"Method: {alignment['method']}")
    add(f"Fit selection:    {alignment['fit_selection']}")
    add(f"Final RMSD:       {alignment['rmsd_after_angstrom']:.3f} A")

    clash = data["clash_detection"]
    add("\n3. CLASH REMOVAL")
    add(f"Removed molecules: {clash['n_removed']} {clash['counts_by_class']}")
    for item in clash["removed_molecules"]:
        add(
            f"  {item['class']:7s} {item['resname']}{item['resid']}: "
            f"{item['min_distance_angstrom']:.3f} A; "
            f"{item['contact_env_atom']} to "
            f"{item['contact_cholesterol_resname']}{item['contact_cholesterol_resid']}:"
            f"{item['contact_cholesterol_atom']}"
        )

    comp = data["composition"]
    add("\n4. COMPOSITION CORRECTION")
    add(f"Bilayer midplane z: {comp['midplane_z_angstrom']:.3f} A")
    add(f"Lipid targets:      {comp['lipid_targets']}")
    add(f"Counts before:      {comp['lipid_counts_before']}")
    add(f"Counts after:       {comp['lipid_counts_after']}")
    add(f"Lipids removed:     {len(comp['lipid_residues_removed'])}")
    add(
        f"Reference salt:     {comp['assumed_reference_concentration_molar']:.4f} M; "
        f"ion counts {comp['reference_ion_counts']}"
    )
    add(f"Ion counts before:  {comp['ion_counts_before']}")
    add(f"Ion counts after:   {comp['ion_counts_after']}")
    add(f"Output ion counts:  {comp['output_ion_counts']}")
    if comp["estimated_final_concentration_molar"] is None:
        add("Final salt estimate: unavailable (no complete reference salt pairs)")
    else:
        add(
            f"Final salt estimate: {comp['estimated_final_concentration_molar']:.4f} M; "
            f"deviation {comp['concentration_deviation_fraction']:.2%}; "
            f"tolerance {comp['concentration_tolerance_fraction']:.2%}; "
            f"within tolerance: {comp['concentration_within_tolerance']}"
        )
    for warning in comp["warnings"]:
        add(f"WARNING: {warning}")

    final = data["final_system"]
    add("\n5. FINAL SYSTEM")
    add(f"Atoms:          {final['n_atoms']}")
    add(f"Residues:       {final['n_residues']}")
    add(f"Box dimensions: {final['box_dimensions']}")
    add(f"Box unchanged:  {final['box_unchanged']}")
    if data.get("topology"):
        add(f"Topology:       {data['topology']}")
    if data.get("ndx"):
        add(f"Index groups:   {data['ndx']}")
    if data.get("refgro_check"):
        add(f"Reference check:{data['refgro_check']}")
    for note in data.get("notes", []):
        add(f"NOTE: {note}")
    add("\n" + "=" * 78)
    if data.get("box_validation"):
        lines.append(render_preflight(data["box_validation"]))
    return "\n".join(lines) + "\n"


def write_report(data, report_path_no_ext: str):
    os.makedirs(os.path.dirname(os.path.abspath(report_path_no_ext)) or ".", exist_ok=True)
    text_path = report_path_no_ext + ".txt"
    json_path = report_path_no_ext + ".json"
    with open(text_path, "w") as handle:
        handle.write(render_text(data))
    with open(json_path, "w") as handle:
        json.dump(data, handle, indent=2, default=str)
    return text_path, json_path
