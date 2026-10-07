"""Text and JSON reports for cholesterol-restoration mode."""

from __future__ import annotations

import json
import os
from datetime import datetime, timezone
from collections import Counter
from box_validation import render_preflight


def _ranges(numbers):
    """[1, 2, 3, 7] -> '1-3, 7'."""
    numbers, parts = sorted(numbers), []
    for n in numbers:
        if parts and n == parts[-1][1] + 1:
            parts[-1][1] = n
        else:
            parts.append([n, n])
    return ", ".join(str(a) if a == b else f"{a}-{b}" for a, b in parts)


def restoration_summary(final_universe, result, charmm_resname, distance_from_protein):
    """Where the restored cholesterols are in the final system, and what was removed for them.

    Final residue numbers are sequential positions in the written GRO/PDB and
    topol.top (1-based); the GRO's own residue IDs keep the experimental
    numbers, so they can repeat across the system.
    """
    from repack import restored_residue_indices
    inserted = result.visualization_inserted
    matched = restored_residue_indices(final_universe, inserted, charmm_resname)
    chl_order = [int(r.resindex) for r in final_universe.residues
                 if str(r.resname).upper() == charmm_resname.upper()]
    position = {resindex: number + 1 for number, resindex in enumerate(chl_order)}
    nearby_by_chl = {}
    for item in result.clash_result.removed_molecules:
        nearby_by_chl.setdefault(int(item.contact_partner_resid), []).append(item)
    restored = []
    for source, resindex in zip(inserted.residues, matched):
        residue = final_universe.residues[resindex]
        removed = sorted(nearby_by_chl.get(int(source.resid), []), key=lambda x: x.min_distance)
        restored.append({
            "experimental": f"{source.resname} {int(source.resid)}",
            "final_residue_number": resindex + 1,
            "resid_in_files": int(residue.resid),
            "atoms": [int(residue.atoms[0].index) + 1, int(residue.atoms[-1].index) + 1],
            "cholesterol_number": position[resindex],
            "cholesterols_in_system": len(chl_order),
            "removed_nearby": [{"resname": x.resname, "resid": int(x.resid), "class": x.molclass,
                                "closest_angstrom": round(float(x.min_distance), 3)} for x in removed],
        })
    nearby = result.clash_result.removed_molecules
    far = result.composition.lipid_residues_removed
    far_distances = [x["distance_to_protein_angstrom"] for x in far if "distance_to_protein_angstrom" in x]
    return {
        "restored": restored,
        "removed_nearby": {
            "count": len(nearby),
            "by_resname": dict(Counter(x.resname for x in nearby)),
            "closest_angstrom": ([round(min(x.min_distance for x in nearby), 3),
                                  round(max(x.min_distance for x in nearby), 3)] if nearby else None),
        },
        "removed_far": {
            "count": len(far),
            "by_resname": dict(Counter(x["resname"] for x in far)),
            "by_leaflet": dict(Counter(f"{x['resname']} {x['leaflet']}" for x in far)),
            "only_beyond_angstrom": distance_from_protein,
            "distance_to_protein_angstrom": ([min(far_distances), max(far_distances)]
                                             if far_distances else None),
            "ions": len(result.composition.ions_removed),
        },
    }


def summary_lines(summary):
    """Short console lines, also read by the GUI's 'Progress at a glance'."""
    restored = summary["restored"]
    lines = []
    if restored:
        total = restored[0]["cholesterols_in_system"]
        lines.append(
            f"Restored cholesterols: {len(restored)}, final residues "
            f"{_ranges(r['final_residue_number'] for r in restored)} (cholesterol "
            f"{_ranges(r['cholesterol_number'] for r in restored)} of {total} in topol.top)")
    near = summary["removed_nearby"]
    if near["count"]:
        low, high = near["closest_angstrom"]
        lines.append(f"Removed next to them: {near['count']} molecule(s) ({_counts_dict(near['by_resname'])}), "
                     f"{low:.2f}-{high:.2f} A from a restored cholesterol")
    else:
        lines.append("Removed next to them: none")
    far = summary["removed_far"]
    if far["count"]:
        low, high = far["distance_to_protein_angstrom"]
        lines.append(f"Removed far away: {far['count']} molecule(s) ({_counts_dict(far['by_resname'])}), "
                     f"{low:.1f}-{high:.1f} A from the protein (only beyond {far['only_beyond_angstrom']:g} A)")
    else:
        lines.append("Removed far away: none")
    return lines


def _counts_dict(counts):
    return ", ".join(f"{name} {count}" for name, count in sorted(counts.items(), key=lambda kv: -kv[1]))


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
    output_validation=None,
    coordinates_written=False,
    repack=None,
    summary=None,
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
            "target_protein_mask": cfg.cholesterol.target_protein_mask,
            "experimental_protein_mask": cfg.cholesterol.experimental_protein_mask,
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
            "coordinates_written": coordinates_written,
            "coordinate_output_validation": output_validation,
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
        "repack": repack,
        "restoration_summary": summary,
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

    summary = data.get("restoration_summary")
    if summary:
        add("\nSUMMARY")
        for line in summary_lines(summary):
            add(line)
        add("\nRestored cholesterols in the final system (residue and atom numbers count from 1 in")
        add("step5_input.gro/pdb and topol.top; the GRO residue ID keeps the experimental number):")
        for item in summary["restored"]:
            add(f"  {item['experimental']:>10s} -> residue {item['final_residue_number']}, atoms "
                f"{item['atoms'][0]}-{item['atoms'][1]}, cholesterol {item['cholesterol_number']} of "
                f"{item['cholesterols_in_system']} (GRO resid {item['resid_in_files']})")
            removed = item["removed_nearby"]
            add("            removed next to it: " + (", ".join(
                f"{x['resname']} {x['resid']} ({x['closest_angstrom']:.2f} A)" for x in removed) or "none"))
        far = summary["removed_far"]
        if far["count"]:
            add(f"Removed far from the protein to meet the leaflet targets (only molecules beyond "
                f"{far['only_beyond_angstrom']:g} A of the protein are eligible):")
            add("  " + ", ".join(f"{k} {v}" for k, v in sorted(far["by_leaflet"].items())))
            low, high = far["distance_to_protein_angstrom"]
            add(f"  closest atom to the protein: {low:.1f}-{high:.1f} A; full list in section 4")
        if far["ions"]:
            add(f"Ions removed: {far['ions']}")

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
    add("\n2. PROTEIN-BASED ALIGNMENT")
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
    for item in sorted(comp["lipid_residues_removed"],
                       key=lambda x: x.get("distance_to_protein_angstrom", 0.0)):
        distance = item.get("distance_to_protein_angstrom")
        add(f"  {item['resname']} {item['resid']} ({item['leaflet']} leaflet)"
            + (f": {distance:.1f} A from the protein" if distance is not None else ""))
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
    add(f"Coordinates written: {final.get('coordinates_written', False)}")
    add(f"Coordinate output validation: {final.get('coordinate_output_validation')}")
    if data.get("topology"):
        add(f"Topology:       {data['topology']}")
    if data.get("ndx"):
        add(f"Index groups:   {data['ndx']}")
    repack = data.get("repack")
    if repack and repack.get("error"):
        add(f"Repacking bundle: NOT written: {repack['error']}")
    elif repack:
        add(f"Repacking bundle: {repack['directory']} ({repack['restored_cholesterols']} restored "
            f"cholesterol(s) as {repack['restored_moleculetype']}, {repack['md_length_ps'] / 1000:g} ns of MD; "
            "topol.top, toppar/, GRO, index and mdp files; a sample run script is also saved)")
        for warning in repack["warnings"]:
            add(f"  WARNING: {warning}")
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
