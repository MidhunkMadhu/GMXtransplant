"""
Human-readable + JSON validation report covering every number the project
spec asks for: alignment quality, clash detection/removal, charge
accounting, and final system validation.
"""

from __future__ import annotations

import json
import os
from dataclasses import asdict
from datetime import datetime, timezone
from typing import Any, Dict
from box_validation import render_preflight


def _f(x, nd=3):
    try:
        if x is None:
            return "n/a"
        if isinstance(x, float) and (x != x):  # NaN
            return "n/a (unresolved -- see unknown_resnames)"
        return f"{x:.{nd}f}"
    except Exception:
        return str(x)


def build_report_dict(
    config,
    align_res,
    rep_res,
    clash_res,
    charge_report,
    final_universe,
    original_universe,
    namefix_summary=None,
    topology_result=None,
    topology_error=None,
    ndx_counts=None,
    ndx_notes=None,
    ndx_error=None,
    refgro_check=None,
    ligand_replace_result=None,
    resolved_replacement_ligands=None,
    output_validation=None,
) -> Dict[str, Any]:
    if ligand_replace_result is not None:
        lrr = ligand_replace_result
        inputs = {
            "structure_path": config.ligand_replace.structure_path,
            "original_ligand_resname": lrr.old_ligand_resname,
            "new_ligand_coord_path": config.ligand_replace.new_ligand.coord_path,
            "new_ligand_resname": lrr.new_ligand_resname,
        }
        alignment = None
        replacement = None
        ligand_replacement = {
            "fit_method": lrr.fit_method,
            "old_ligand_resname": lrr.old_ligand_resname,
            "new_ligand_resname": lrr.new_ligand_resname,
            "n_fit_atoms": lrr.n_fit_atoms,
            "fit_atom_pairs": lrr.fit_atom_pairs,
            "fit_rmsd_angstrom": lrr.fit_rmsd,
            "n_original_ligand_atoms_removed": lrr.n_original_ligand_atoms_removed,
            "n_new_ligand_atoms_inserted": lrr.n_new_ligand_atoms_inserted,
            "notes": lrr.notes,
        }
    else:
        inputs = {
            "target_box": config.target_box.path,
            "target_box_receptor_mask": config.target_box.receptor_mask,
            "replacement_structure": config.replacement_structure.path,
            "replacement_structure_receptor_mask": config.replacement_structure.receptor_mask,
        }
        alignment = {
            "method": config.alignment.method,
            "region_mask_original": config.alignment.region_mask_original,
            "region_mask_replacement": config.alignment.region_mask_replacement,
            "fit_selection": align_res.fit_atom_selection,
            "n_atom_pairs_before_refinement": align_res.n_atom_pairs_before,
            "n_atom_pairs_after_refinement": align_res.n_atom_pairs_after,
            "excluded_as_outliers": align_res.excluded_as_outliers,
            "cycles_executed": align_res.n_cycles_executed,
            "rmsd_before_angstrom": align_res.rmsd_before,
            "rmsd_after_angstrom": align_res.rmsd_after,
            "n_residues_in_replacement_block": len(align_res.mobile_block_global_resnums),
            "n_residues_in_original_receptor_block": len(align_res.target_block_global_resnums),
            "notes": align_res.notes,
        }
        replacement = {
            "selection_policy": (
                "insert every residue selected by replacement_structure.receptor_mask"
            ),
            "n_original_receptor_atoms_removed": rep_res.n_original_receptor_atoms_removed,
            "n_replacement_atoms_inserted": rep_res.n_replacement_atoms_inserted,
            "inserted_residue_counts": rep_res.inserted_residue_counts,
            "replacement_ligands": [
                {
                    "selector": ligand.selector,
                    "resname": ligand.resname,
                    "moleculetype": ligand.moleculetype,
                    "net_charge_per_residue": ligand.net_charge,
                    "itp_charge_sum": ligand.itp_charge,
                    "charge_source": ligand.charge_source,
                    "itp_path": ligand.itp_path,
                }
                for ligand in (resolved_replacement_ligands or [])
            ],
            "environment_atoms_untouched": rep_res.environment_ag_size,
            "notes": rep_res.notes,
        }
        ligand_replacement = None

    return {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "box_validation": getattr(original_universe, "box_preflight", None),
        "mode": "ligand_replace" if ligand_replace_result is not None else "receptor_replace",
        "inputs": inputs,
        "name_restoration": None if not namefix_summary else {
            label: {
                "renamed_residue_count": r.renamed_residue_count,
                "renames_by_pair": {f"{a}->{b}": c for (a, b), c in r.renames_by_pair.items()},
                "unresolved_resnames": r.unresolved_resnames,
                "method": getattr(r, "method", "atom_signature"),
                "source_path": getattr(r, "source_path", None),
                "candidate_atom_counts": r.candidate_atom_counts,
                "candidate_itp_paths": getattr(r, "candidate_itp_paths", {}),
                "reference_topology_lipid_types": getattr(
                    r, "reference_lipid_types", []
                ),
                "environment_toppar_lipid_types": getattr(
                    r, "toppar_lipid_types", []
                ),
                "reference_lipid_atom_counts": getattr(
                    r, "reference_lipid_atom_counts", {}
                ),
                "reference_lipid_itp_paths": getattr(
                    r, "reference_lipid_itp_paths", {}
                ),
                "notes": r.notes,
            }
            for label, r in namefix_summary.items()
        },
        "alignment": alignment,
        "replacement": replacement,
        "ligand_replacement": ligand_replacement,
        "clash_detection": {
            "threshold_angstrom_by_class": clash_res.threshold_used,
            "use_pbc": config.clash_detection.use_pbc,
            "heavy_atoms_only": config.clash_detection.heavy_atoms_only,
            "counts_by_class": clash_res.counts_by_class,
            "counts_by_resname": clash_res.counts_by_resname,
            "n_lipids_removed": clash_res.n_lipids_removed,
            "n_lipids_total_before_removal": clash_res.n_lipids_total,
            "lipid_removal_fraction": clash_res.lipid_removal_fraction,
            "lipid_removal_flagged": clash_res.lipid_removal_flagged,
            "other_class_resnames_seen": clash_res.other_class_resnames,
            "n_molecules_removed_total": len(clash_res.removed_molecules),
            "removed_molecules": [
                {
                    "resname": m.resname, "resid": m.resid, "class": m.molclass,
                    "min_distance_angstrom": round(m.min_distance, 3), "n_atoms": m.n_atoms,
                    "cutoff_angstrom": m.cutoff_used,
                    "contact": f"{m.contact_env_atom} -- {m.contact_partner_resname}{m.contact_partner_resid}:{m.contact_partner_atom}",
                    "contact_context": m.contact_context,
                    "removal_decision": m.removal_decision,
                }
                for m in clash_res.removed_molecules
            ],
            "keep_classes": config.clash_detection.keep_classes,
            "n_flagged_but_kept": len(clash_res.flagged_but_kept),
            "flagged_but_kept": [
                {
                    "resname": m.resname, "resid": m.resid, "class": m.molclass,
                    "min_distance_angstrom": round(m.min_distance, 3), "n_atoms": m.n_atoms,
                    "cutoff_angstrom": m.cutoff_used,
                    "contact": f"{m.contact_env_atom} -- {m.contact_partner_resname}{m.contact_partner_resid}:{m.contact_partner_atom}",
                    "contact_context": m.contact_context,
                    "removal_decision": m.removal_decision,
                }
                for m in clash_res.flagged_but_kept
            ],
            "flagged_but_kept_counts_by_resname": (
                clash_res.flagged_but_kept_counts_by_resname
            ),
            "notes": clash_res.notes,
        },
        "charge": {
            "unknown_resnames": charge_report.unknown_resnames,
            "charge_original_receptor": charge_report.charge_original_receptor,
            "charge_replacement_receptor": charge_report.charge_replacement_receptor,
            "charge_diff_from_replacement": charge_report.charge_diff_from_replacement,
            "charge_removed_by_clashes": charge_report.charge_removed_by_clashes,
            "charge_removed_by_class": charge_report.charge_removed_by_class,
            "charge_removed_by_resname": charge_report.charge_removed_by_resname,
            "net_charge_before_clash_removal": charge_report.net_charge_before_clash_removal,
            "net_charge_after_clash_removal": charge_report.net_charge_after_clash_removal,
            "target_net_charge": charge_report.target_net_charge,
            "charge_excess_before_neutralization": charge_report.charge_excess,
            "neutralization_performed": charge_report.neutralization_performed,
            "n_ions_removed_for_neutralization": len(charge_report.ions_removed),
            "ions_removed_for_neutralization": charge_report.ions_removed,
            "ion_selection_summary": charge_report.ion_selection_summary,
            "final_net_charge": charge_report.final_net_charge,
            "stopped_reason": charge_report.stopped_reason,
            "notes": charge_report.notes,
        },
        "final_system": None if final_universe is None else {
            "n_atoms": len(final_universe.atoms),
            "n_residues": len(final_universe.residues),
            "box_dimensions": list(map(float, final_universe.dimensions)),
            "target_box_dimensions": list(map(float, original_universe.dimensions)),
            "box_unchanged": bool(
                (final_universe.dimensions == original_universe.dimensions).all()
            ),
            "coordinate_output_validation": output_validation,
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
    }


def render_text_report(d: Dict[str, Any]) -> str:
    lines = []
    A = lines.append
    lr_mode = d.get("mode") == "ligand_replace"
    A("=" * 78)
    A(
        "LIGAND REPLACEMENT WORKFLOW - VALIDATION REPORT"
        if lr_mode else
        "RECEPTOR REPLACEMENT WORKFLOW - VALIDATION REPORT"
    )
    A(f"Generated: {d['generated_at_utc']}")
    A("=" * 78)

    A("\n--- INPUTS ---")
    if lr_mode:
        A(f"Structure (unchanged except one ligand): {d['inputs']['structure_path']}")
        A(f"  original ligand removed: {d['inputs']['original_ligand_resname']}")
        A(f"New ligand:              {d['inputs']['new_ligand_coord_path']}")
        A(f"  inserted as resname:   {d['inputs']['new_ligand_resname']}")
    else:
        A(f"Target box:            {d['inputs']['target_box']}")
        A(f"  receptor mask removed: {d['inputs']['target_box_receptor_mask']}")
        A(f"Replacement structure:   {d['inputs']['replacement_structure']}")
        A(f"  receptor mask used:    {d['inputs']['replacement_structure_receptor_mask']}")

    nrfix = d.get("name_restoration")
    if nrfix:
        A("\n--- 0. RESNAME RESTORATION (PDB truncation fix) ---")
        for label, r in nrfix.items():
            A(f"{label}: method {r['method']}")
            if r["source_path"]:
                A(f"{label}: reference topology {r['source_path']}")
            if r["reference_topology_lipid_types"]:
                A(
                    f"{label}: reference topology lipid types "
                    f"{r['reference_topology_lipid_types']}"
                )
            if r["environment_toppar_lipid_types"]:
                A(
                    f"{label}: environment toppar lipid types "
                    f"{r['environment_toppar_lipid_types']}"
                )
            if r["reference_lipid_atom_counts"]:
                A(
                    f"{label}: reference lipid ITP atom counts "
                    f"{r['reference_lipid_atom_counts']}"
                )
                A(
                    f"{label}: reference lipid ITP files "
                    f"{r['reference_lipid_itp_paths']}"
                )
            A(f"{label}: candidate atom counts {r['candidate_atom_counts']}")
            if r["candidate_itp_paths"]:
                A(f"{label}: candidate ITP files {r['candidate_itp_paths']}")
            A(f"{label}: renamed {r['renamed_residue_count']} residues: {r['renames_by_pair']}")
            if r["unresolved_resnames"]:
                A(f"{label}: could not confidently restore: {r['unresolved_resnames']}")
            for note in r["notes"]:
                A(f"NOTE: {note}")

    if lr_mode:
        lrp = d["ligand_replacement"]
        A("\n--- 1. LIGAND REPLACEMENT (receptor/membrane/solvent/ions untouched) ---")
        A(f"Fit method:                    {lrp['fit_method']}")
        A(f"Original ligand resname:       {lrp['old_ligand_resname']}")
        A(f"New ligand resname:            {lrp['new_ligand_resname']}")
        if lrp["fit_method"] in ("pairfit", "autofit"):
            A(f"Fit atoms used ({lrp['n_fit_atoms']}, old->new): {lrp['fit_atom_pairs']}")
            A(f"Fit RMSD (over fit atoms):     {_f(lrp['fit_rmsd_angstrom'])} A")
        A(f"Original ligand atoms removed: {lrp['n_original_ligand_atoms_removed']}")
        A(f"New ligand atoms inserted:     {lrp['n_new_ligand_atoms_inserted']}")
        for n in lrp["notes"]:
            A(f"NOTE: {n}")
    else:
        al = d["alignment"]
        if al["method"] == "mask_fit":
            # No PyMOL, no sequence alignment, no outlier rejection for this
            # method -- a dedicated, honestly-labeled block instead of the
            # PyMOL-flavored one below (which would otherwise print
            # "RMSD after outlier rejection" for a fit that never rejected
            # any outliers, and a meaningless "Refinement cycles executed: 0").
            A("\n--- 1. RECEPTOR ALIGNMENT (mask_fit: direct Kabsch/SVD fit, no outlier rejection) ---")
            A(f"Method:                        {al['method']}")
            A(f"Fit selection:                 {al['fit_selection']}")
            A(f"Residues in replacement block: {al['n_residues_in_replacement_block']}")
            A(f"Residues in box receptor block:{al['n_residues_in_original_receptor_block']}")
            A(f"Atom pairs fit (positional):   {al['n_atom_pairs_after_refinement']}")
            A(f"RMSD (positionally-matched fit, no outlier rejection): {_f(al['rmsd_after_angstrom'])} A")
        else:
            A("\n--- 1. RECEPTOR ALIGNMENT (PyMOL automatic align) ---")
            A(f"Method:                        {al['method']}")
            A(f"Fit selection:                 {al['fit_selection']}")
            A(f"Residues in replacement block: {al['n_residues_in_replacement_block']}")
            A(f"Residues in box receptor block:{al['n_residues_in_original_receptor_block']}")
            A(f"Atom pairs used (initial fit): {al['n_atom_pairs_before_refinement']}")
            A(f"Atom pairs used (final fit):   {al['n_atom_pairs_after_refinement']}")
            A(f"Excluded as outliers:          {al['excluded_as_outliers']}")
            A(f"Refinement cycles executed:    {al['cycles_executed']}")
            A(f"RMSD before outlier rejection: {_f(al['rmsd_before_angstrom'])} A")
            A(f"RMSD after outlier rejection:  {_f(al['rmsd_after_angstrom'])} A")
        for n in al["notes"]:
            A(f"NOTE: {n}")

        rp = d["replacement"]
        A("\n--- 2. RECEPTOR BLOCK REPLACEMENT ---")
        A(f"Selection policy:                  {rp['selection_policy']}")
        A(f"Original receptor atoms removed:   {rp['n_original_receptor_atoms_removed']}")
        A(f"Replacement atoms inserted:        {rp['n_replacement_atoms_inserted']}")
        A(f"Inserted residue counts:           {rp['inserted_residue_counts']}")
        if rp["replacement_ligands"]:
            A("Incoming ligand charge/topology metadata:")
            for ligand in rp["replacement_ligands"]:
                A(
                    f"  {ligand['resname']}: charge {ligand['net_charge_per_residue']:+.6f} e "
                    f"per residue ({ligand['charge_source']}), moleculetype "
                    f"{ligand['moleculetype']}, ITP {ligand['itp_path']}"
                )
        A(f"Environment atoms left untouched:  {rp['environment_atoms_untouched']}")
        for n in rp["notes"]:
            A(f"NOTE: {n}")

    cd = d["clash_detection"]
    A("\n--- 3. STERIC CLASH DETECTION & REMOVAL ---")
    A(f"Clash thresholds used (A): {cd['threshold_angstrom_by_class']}")
    A(f"PBC-aware: {cd['use_pbc']}   Heavy-atoms-only: {cd['heavy_atoms_only']}")
    A(f"Molecules removed by class: {cd['counts_by_class']}")
    A(f"Molecules removed by resname: {cd['counts_by_resname']}")
    A(f"Total molecules removed:    {cd['n_molecules_removed_total']}")
    A(f"Lipids removed:             {cd['n_lipids_removed']} / {cd['n_lipids_total_before_removal']} "
      f"({cd['lipid_removal_fraction']:.1%})"
      f"{'  <-- FLAGGED for manual inspection' if cd['lipid_removal_flagged'] else ''}")
    if cd["other_class_resnames_seen"]:
        A(f"Unrecognized environment resnames (classified 'other'): {cd['other_class_resnames_seen']}")
    if cd.get("keep_classes"):
        A(f"Classes checked but never removed (clash_detection.keep_classes): {cd['keep_classes']}")
        A(f"Flagged but kept (within threshold, not removed): {cd['n_flagged_but_kept']}")
        for m in cd["flagged_but_kept"]:
            A(f"  {m['class']:7s} {m['resname']}{m['resid']:<6d} "
              f"{_f(m['min_distance_angstrom'])} A   contact: {m['contact']} "
              f"({m['contact_context']}); cutoff={m['cutoff_angstrom']} A; "
              f"decision={m['removal_decision']}")
    if cd["removed_molecules"]:
        A("Removed molecule clash details:")
        for m in cd["removed_molecules"]:
            A(f"  {m['class']:7s} {m['resname']}{m['resid']:<6d} "
              f"{_f(m['min_distance_angstrom'])} A   contact: {m['contact']} "
              f"({m['contact_context']}); cutoff={m['cutoff_angstrom']} A; "
              f"decision={m['removal_decision']}")
    for n in cd["notes"]:
        A(f"NOTE: {n}")

    ch = d["charge"]
    A("\n--- 4. CHARGE CALCULATION & NEUTRALIZATION ---")
    if ch["unknown_resnames"]:
        A(f"UNKNOWN RESIDUE CHARGES (must be supplied in config): {ch['unknown_resnames']}")
        A(f"STOPPED: {ch['stopped_reason']}")
    else:
        label_orig = "original ligand" if lr_mode else "original receptor block"
        label_new = "new ligand" if lr_mode else "replacement receptor block"
        A(f"Charge of {label_orig}:".ljust(41) + f"{_f(ch['charge_original_receptor'])} e")
        A(f"Charge of {label_new}:".ljust(41) + f"{_f(ch['charge_replacement_receptor'])} e")
        A("Charge diff (new - original):".ljust(41) + f"{_f(ch['charge_diff_from_replacement'])} e")
        A(f"Charge removed by clash deletion:        {_f(ch['charge_removed_by_clashes'])} e "
          f"{ch['charge_removed_by_class']}")
        A(f"Charge removed by resname:                {ch['charge_removed_by_resname']}")
        A(f"Net charge BEFORE clash removal:         {_f(ch['net_charge_before_clash_removal'])} e")
        A(f"Net charge AFTER clash removal:          {_f(ch['net_charge_after_clash_removal'])} e")
        A(f"Target net charge:                       {_f(ch['target_net_charge'])} e")
        if ch["neutralization_performed"]:
            A(f"Ions removed for neutralization:         {ch['n_ions_removed_for_neutralization']}")
            from collections import Counter
            resn_counts = Counter(r["resname"] for r in ch["ions_removed_for_neutralization"])
            A(f"  by resname: {dict(resn_counts)}")
            A(f"FINAL NET CHARGE:                        {_f(ch['final_net_charge'])} e")
        if ch.get("ion_selection_summary"):
            A(f"Ion selection summary:                   {ch['ion_selection_summary']}")
        for n in ch.get("notes", []):
            A(f"NOTE: {n}")
        if ch["stopped_reason"]:
            A(f"STOPPED: {ch['stopped_reason']}")

    fs = d["final_system"]
    A("\n--- 5. FINAL SYSTEM VALIDATION ---")
    if fs is None:
        A("Final system was NOT written (pipeline stopped before output; see above).")
    else:
        A(f"Final atom count:     {fs['n_atoms']}")
        A(f"Final residue count:  {fs['n_residues']}")
        A(f"Box dimensions:       {fs['box_dimensions']}")
        A(f"Target box dims:    {fs['target_box_dimensions']}")
        A(f"Box unchanged:        {fs['box_unchanged']}")
        A(f"Coordinate output validation: {fs.get('coordinate_output_validation')}")

    tp = d.get("topology")
    if tp is not None:
        A("\n--- 6. TOPOLOGY ASSEMBLY (topol.top + toppar/) ---")
        if tp["error"]:
            A(f"NOT WRITTEN -- {tp['error']}")
        else:
            A(f"Wrote: {tp['top_path']}")
            A(f"toppar/ directory: {tp['toppar_dir']}")
            A(f"[ molecules ] written: {tp['molecules_written']}")
            A(f"#include files copied: {len(tp['includes_written'])}")
            if tp.get("charge_audit"):
                A(f"Final topology/coordinate charge audit: {tp['charge_audit']}")
        for n in tp["notes"]:
            A(f"NOTE: {n}")

    ix = d.get("ndx")
    if ix is not None:
        A("\n--- 7. GROMACS INDEX (index.ndx) ---")
        if ix["error"]:
            A(f"NOT WRITTEN -- {ix['error']}")
        else:
            A(f"Group sizes (atoms): {ix['counts']}")
        for n in ix["notes"]:
            A(f"NOTE: {n}")

    rg = d.get("refgro_check")
    if rg is not None:
        A("\n--- 8. REFERENCE-GRO COMPLETENESS CHECK ---")
        A(f"refgro: {rg['refgro_path']}")
        A(f"Non-protein resnames in refgro: {rg['refgro_nonprotein_resnames']}")
        if rg["missing_from_final"]:
            A(f"NOTE: present in refgro but ABSENT from the final system: {rg['missing_from_final']}")
        else:
            A("All of refgro's non-protein resnames are present in the final system.")

    A("\n" + "=" * 78)
    if d.get("box_validation"):
        lines.append(render_preflight(d["box_validation"]))
    return "\n".join(lines)


def write_report(d: Dict[str, Any], report_path_no_ext: str):
    os.makedirs(os.path.dirname(os.path.abspath(report_path_no_ext)), exist_ok=True)
    text = render_text_report(d)
    with open(report_path_no_ext + ".txt", "w") as fh:
        fh.write(text)
    with open(report_path_no_ext + ".json", "w") as fh:
        json.dump(d, fh, indent=2, default=str)
    return report_path_no_ext + ".txt", report_path_no_ext + ".json"
