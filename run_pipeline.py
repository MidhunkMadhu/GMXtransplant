#!/usr/bin/env python3
"""Entry point for receptor, ligand, and cholesterol modification modes.

Usage:
    python3 run_pipeline.py --mode receptor -i receptor_replace.yaml
    python3 run_pipeline.py --mode lig -i ligand_replace.yaml
    python3 run_pipeline.py --mode chl -i cholesterol_restore.yaml

RECEPTOR MODE:
implements the project's core workflow:
  1-3.  read target box, identify its receptor, read replacement structure
  4-6.  identify shared atoms, align (PyMOL), apply transform to full block
  7-9.  remove target box's selected block, insert the full aligned replacement block

LIGAND MODE:
the receptor/membrane/solvent/ions do NOT change at all and only a bound
ligand is swapped for a different one (see ligand_replace.py's module
docstring):
  1-3.  read the one structure, remove the named old ligand, load + position
        (pairfit, autofit, or nofit) the new ligand from its own file

Both replacement modes then share the same downstream steps:
  10-11. detect steric clashes (PBC-aware), remove complete clashing molecules
  12-13. calculate net charge, restore intended charge via counterion removal
  14.   write PDB/GRO + validate + report
  15.   (optional) assemble topol.top, toppar/, and index.ndx for the final system.
  16.   (optional) refgro completeness check from the active mode YAML.

Fails loudly (with a clear message, no partial/garbage output) rather than
silently producing a system with an unintended composition or net charge.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
import sysconfig
import traceback
from dataclasses import replace as dc_replace

import numpy as np
import MDAnalysis as mda

from config import load_config, ConfigError
from align import align_replacement_to_box, AlignmentError
from replace import assemble_replacement, build_environment, ReplaceError
from replacement_ligands import resolve_replacement_ligands, ReplacementLigandError
from ligand_replace import assemble_ligand_replacement, LigandReplaceError
from clashes import detect_and_remove_clashes
from charges import compute_pre_neutralization_report, neutralize, ChargeError
from output import (
    finalize_system, write_outputs, write_inspection_pdb,
    refgro_completeness_check, OutputError,
)
from report import build_report_dict, render_text_report, write_report
from masks import resolve_mask, MaskError
from namefix import (
    build_resname_templates,
    build_itp_name_reference,
    restore_full_resnames,
    restore_full_resnames_from_itp_atom_counts,
    NameFixError,
)
from classify import AMINO_ACID_RESNAMES
from topology import (
    assemble_topology,
    build_receptor_chain_layout,
    build_topology_charge_model,
    audit_final_topology,
    TopologyError,
)
from ndx import write_ndx
from cholesterol import restore_cholesterols, CholesterolError
from box_validation import (
    BoxValidationError, validate_box, read_environment,
    validate_staged_hard_clashes,
)
from cholesterol_report import (
    build_report as build_cholesterol_report,
    write_report as write_cholesterol_report,
)

__version__ = "0.1.0"

_EXAMPLE_FILES = {
    "receptor": "receptor_replace.yaml",
    "lig": "ligand_replace.yaml",
    "chl": "cholesterol_restore.yaml",
}


def _read_example(mode: str) -> str:
    filename = _EXAMPLE_FILES[mode]
    candidates = (
        Path(__file__).resolve().parent / filename,
        Path(sysconfig.get_path("data"))
        / "share"
        / "gmxtransplant"
        / "examples"
        / filename,
    )
    for path in candidates:
        if path.is_file():
            return path.read_text(encoding="utf-8")
    raise ConfigError(
        f"Bundled example '{filename}' was not found. Reinstall GMXtransplant."
    )


class _StageFailure(Exception):
    """Internal control-flow only: a prepare-stage helper raises this after
    printing its own error message, so main() can just return the exit code."""
    def __init__(self, exit_code: int):
        self.exit_code = exit_code


def _read_box_with_override(
    label: str,
    path: str,
    box_dimensions,
    fatal_if_missing: bool,
    format: str = "auto",
) -> "mda.Universe":
    explicit_format = None if format in (None, "", "auto") else format.upper()
    try:
        u = mda.Universe(path, format=explicit_format) if explicit_format else mda.Universe(path)
    except Exception as exc:
        print(f"[INPUT ERROR] Cannot read {label} ('{path}'): {exc}", file=sys.stderr)
        raise _StageFailure(2) from exc
    if box_dimensions is not None:
        u.dimensions = np.asarray(box_dimensions, dtype=np.float32)
    elif u.dimensions is None:
        msg = (
            f"{label} ('{path}') has no readable box dimensions "
            f"(for example, a PDB without CRYST1). Use an authoritative box-bearing "
            f"input or set matching box_dimensions: [a, b, c, alpha, beta, gamma] "
            f"in Å and degrees if this structure needs periodic processing."
        )
        if fatal_if_missing:
            print(f"[CONFIG ERROR] {msg}", file=sys.stderr)
            raise _StageFailure(2)
        else:
            print(f"      WARNING: {msg}")
    if u.dimensions is not None:
        try:
            validation = validate_box(
                u.dimensions, label=f"{label} box"
            )
        except BoxValidationError as exc:
            print(f"[CONFIG ERROR] {exc}", file=sys.stderr)
            raise _StageFailure(2)
        for note in validation.warnings:
            print(f"      WARNING: {note}")
    return u


def _prepare_receptor_replace(cfg, resolved_replacement_ligands):
    """Front half of receptor-replace mode. Returns a dict of everything
    the shared downstream tail needs."""
    print(f"[1-3] Reading target box ({cfg.target_box.path}) and "
          f"replacement structure ({cfg.replacement_structure.path}) ...")
    try:
        u_orig = read_environment(cfg.target_box.path, cfg.target_box.box_dimensions,
                                  cfg.box_validation, format=cfg.target_box.format)
    except BoxValidationError as exc:
        raise _StageFailure(2) from exc
    u_repl = _read_box_with_override(
        "replacement_structure.path",
        cfg.replacement_structure.path,
        cfg.replacement_structure.box_dimensions,
        fatal_if_missing=False,
        format=cfg.replacement_structure.format,
    )

    print(f"      target box:   {len(u_orig.atoms)} atoms, {len(u_orig.residues)} residues, "
          f"box (lengths Å, angles degrees)={None if u_orig.dimensions is None else list(u_orig.dimensions)}")
    print(f"      replacement:  {len(u_repl.atoms)} atoms, {len(u_repl.residues)} residues, "
          f"box (lengths Å, angles degrees)={None if u_repl.dimensions is None else list(u_repl.dimensions)}")
    _target_receptor_ag = resolve_mask(u_orig, cfg.target_box.receptor_mask)
    _mobile_receptor_ag = resolve_mask(u_repl, cfg.replacement_structure.receptor_mask)
    print(f"      target_box.receptor_mask '{cfg.target_box.receptor_mask}' resolves to "
          f"{len(_target_receptor_ag)} atoms / {len(_target_receptor_ag.residues)} residues")
    print(f"      replacement_structure.receptor_mask '{cfg.replacement_structure.receptor_mask}' resolves to "
          f"{len(_mobile_receptor_ag)} atoms / {len(_mobile_receptor_ag.residues)} residues")

    selected_resindices = set(int(i) for i in _mobile_receptor_ag.residues.resindices)
    for ligand in resolved_replacement_ligands:
        matches = [
            residue for residue in u_repl.residues
            if str(residue.resname).upper() == ligand.resname
        ]
        if ligand.occurrence is None and len(matches) > 1:
            choices = ", ".join(
                f"{ligand.resname}:{i}" for i in range(1, len(matches) + 1)
            )
            print(
                f"[CONFIG ERROR] replacement_ligands selector '{ligand.selector}' "
                f"matches {len(matches)} residues. Use one of: {choices}",
                file=sys.stderr,
            )
            raise _StageFailure(2)
        if ligand.occurrence is not None and ligand.occurrence > len(matches):
            print(
                f"[CONFIG ERROR] replacement_ligands selector '{ligand.selector}' "
                f"requests occurrence {ligand.occurrence}, but only {len(matches)} "
                f"{ligand.resname} residue(s) exist.",
                file=sys.stderr,
            )
            raise _StageFailure(2)
        chosen = (
            matches[ligand.occurrence - 1]
            if ligand.occurrence is not None and matches else
            (matches[0] if matches else None)
        )
        if chosen is None or int(chosen.resindex) not in selected_resindices:
            print(
                f"[CONFIG ERROR] replacement_ligands declares {ligand.selector}, but that "
                "specific residue is not selected by "
                f"replacement_structure.receptor_mask '{cfg.replacement_structure.receptor_mask}'.",
                file=sys.stderr,
            )
            raise _StageFailure(2)
        print(
            f"      incoming ligand {ligand.selector}: selected residue position "
            f"{int(chosen.resindex) + 1}, "
            f"charge {ligand.net_charge:+.6f} e per residue, moleculetype "
            f"'{ligand.moleculetype}'"
        )

    namefix_summary = None
    if cfg.name_restoration.enabled:
        nr = cfg.name_restoration
        try:
            if nr.method == "itp_atom_count":
                print(
                    "[0]   Restoring target-environment PDB resnames using "
                    "reference topology names and environment ITP atom counts ..."
                )
                name_reference = build_itp_name_reference(
                    nr.reference_topol,
                    cfg.topology.environment_toppar_dir,
                    nr.pdb_to_full_resname,
                )
                target_environment_ag, _ = build_environment(
                    u_orig, cfg.target_box.receptor_mask
                )
                res = restore_full_resnames_from_itp_atom_counts(
                    target_environment_ag,
                    nr.pdb_to_full_resname,
                    name_reference,
                    nr.reference_topol,
                )
                namefix_summary = {"target_box": res}
                print(
                    "      reference_topol lipid types: "
                    f"{name_reference.reference_lipid_types}"
                )
                print(
                    "      environment toppar lipid types: "
                    f"{name_reference.toppar_lipid_types}"
                )
                if name_reference.reference_lipid_types:
                    print("      reference lipid ITP atom counts:")
                    for name in name_reference.reference_lipid_types:
                        print(
                            f"        {name}: {name_reference.atom_counts[name]} atoms, "
                            f"{name_reference.itp_paths[name]}"
                        )
                if res.candidate_atom_counts:
                    print("      ambiguous-candidate ITP definitions:")
                    for name, count in res.candidate_atom_counts.items():
                        print(
                            f"        {name}: {count} atoms, "
                            f"{res.candidate_itp_paths[name]}"
                        )
                print(
                    f"      target_box: renamed {res.renamed_residue_count} residues "
                    f"({dict(res.renames_by_pair)})"
                )
                for note in res.notes:
                    print(f"      NOTE: {note}")
            else:
                print(
                    "[0]   Restoring truncated PDB resnames using atom-name "
                    "signatures from the reference structure ..."
                )
                ref_u, ref_mask = (
                    (u_orig, cfg.target_box.receptor_mask)
                    if nr.reference == "target_box"
                    else (u_repl, cfg.replacement_structure.receptor_mask)
                )
                ref_env_ag, _ = build_environment(ref_u, ref_mask)
                templates = build_resname_templates(ref_env_ag)
                already_full = set(templates.keys()) | AMINO_ACID_RESNAMES
                targets = []
                if nr.apply_to in ("target_box", "both"):
                    targets.append(("target_box", u_orig))
                if nr.apply_to in ("replacement_structure", "both"):
                    targets.append(("replacement_structure", u_repl))
                namefix_summary = {}
                for label, u in targets:
                    res = restore_full_resnames(
                        u,
                        templates,
                        already_full_resnames=already_full,
                        min_jaccard=nr.min_jaccard,
                    )
                    namefix_summary[label] = res
                    print(
                        f"      {label}: renamed {res.renamed_residue_count} residues "
                        f"({dict(res.renames_by_pair)})"
                    )
                    if res.unresolved_resnames:
                        print(
                            f"      {label}: could not confidently restore: "
                            f"{res.unresolved_resnames}"
                        )
        except (NameFixError, TopologyError) as exc:
            print(f"[NAME RESTORATION ERROR] {exc}", file=sys.stderr)
            raise _StageFailure(2)

    print(f"[4-6] Aligning replacement receptor onto box receptor (alignment.method: {cfg.alignment.method}) ...")
    try:
        align_res = align_replacement_to_box(
            u_orig, u_repl,
            cfg.target_box.receptor_mask, cfg.replacement_structure.receptor_mask,
            cfg.alignment, workdir=".",
        )
    except AlignmentError as e:
        print(f"[ALIGNMENT ERROR] {e}", file=sys.stderr)
        raise _StageFailure(3)
    if cfg.alignment.method == "mask_fit":
        print(f"      RMSD (positionally-matched, no outlier rejection): "
              f"{align_res.rmsd_after:.3f} A ({align_res.n_atom_pairs_after} atom pairs)")
    else:
        print(f"      RMSD before/after outlier rejection: "
              f"{align_res.rmsd_before:.3f} / {align_res.rmsd_after:.3f} A "
              f"({align_res.n_atom_pairs_after}/{align_res.n_atom_pairs_before} atom pairs kept)")

    # Inspection file: box's original receptor overlaid with the newly
    # aligned replacement receptor, written BEFORE the original is removed.
    target_full_ag = resolve_mask(u_orig, cfg.target_box.receptor_mask)
    mobile_full_ag = resolve_mask(u_repl, cfg.replacement_structure.receptor_mask)
    mobile_full_ag.positions = align_res.new_positions
    write_inspection_pdb(target_full_ag, mobile_full_ag, cfg.output.inspection_pdb_path)
    print(f"      Inspection structure written: {cfg.output.inspection_pdb_path}")

    print("[7-9] Removing the target block and inserting the full aligned replacement block ...")
    try:
        rep_res = assemble_replacement(
            u_orig, u_repl,
            cfg.target_box.receptor_mask, cfg.replacement_structure.receptor_mask,
            align_res,
        )
    except ReplaceError as e:
        print(f"[REPLACEMENT ERROR] {e}", file=sys.stderr)
        raise _StageFailure(4)
    print(f"      Removed {rep_res.n_original_receptor_atoms_removed} atoms, "
          f"inserted {rep_res.n_replacement_atoms_inserted} atoms")

    merged = rep_res.merged_universe
    env_ag = merged.atoms[: rep_res.environment_ag_size]
    ins_ag = merged.atoms[rep_res.environment_ag_size:]
    _, removed_ag_orig = build_environment(u_orig, cfg.target_box.receptor_mask)
    environment_ag_orig, _ = build_environment(u_orig, cfg.target_box.receptor_mask)

    return {
        "merged": merged,
        "env_ag": env_ag,
        "ins_ag": ins_ag,
        "removed_ag_orig": removed_ag_orig,
        "environment_ag_orig": environment_ag_orig,
        "clash_cfg": cfg.clash_detection,
        "original_universe_for_report": u_orig,
        "report_kwargs": {
            "align_res": align_res,
            "rep_res": rep_res,
            "namefix_summary": namefix_summary,
            "resolved_replacement_ligands": resolved_replacement_ligands,
        },
    }


def _prepare_ligand_replace(cfg):
    """Front half of ligand-only replace mode. Returns the same shape of
    dict as _prepare_receptor_replace so the downstream tail is identical."""
    lr = cfg.ligand_replace
    print(f"[1-3] Reading structure ({lr.structure_path}) ...")
    try:
        u = read_environment(lr.structure_path, lr.box_dimensions, cfg.box_validation, format=lr.format)
    except BoxValidationError as exc:
        raise _StageFailure(2) from exc
    print(f"      structure:  {len(u.atoms)} atoms, {len(u.residues)} residues, "
          f"box (lengths Å, angles degrees)={None if u.dimensions is None else list(u.dimensions)}")

    print(f"[4-6] Removing '{lr.original_ligand.resname}' and positioning the new ligand "
          f"({lr.fit.method}) from '{lr.new_ligand.coord_path}' ...")
    try:
        lrr = assemble_ligand_replacement(u, lr)
    except LigandReplaceError as e:
        print(f"[LIGAND REPLACE ERROR] {e}", file=sys.stderr)
        raise _StageFailure(4)
    print(f"      Removed {lrr.n_original_ligand_atoms_removed} atom(s) of '{lrr.old_ligand_resname}', "
          f"inserted {lrr.n_new_ligand_atoms_inserted} atom(s) of '{lrr.new_ligand_resname}'")
    if lrr.fit_method in ("pairfit", "autofit"):
        atom_description = "heavy atom(s)" if lrr.fit_method == "autofit" else "matched atom(s)"
        print(f"      Fit RMSD over {lrr.n_fit_atoms} {atom_description}: {lrr.fit_rmsd:.3f} A")
    for n in lrr.notes:
        print(f"      NOTE: {n}")

    # Inspection file: old ligand pose (chain X) overlaid with the new
    # ligand's final positioned pose (chain Y) -- same visual-QC idea as
    # the receptor-mode alignment overlay, reusing the same writer.
    write_inspection_pdb(lrr.old_ligand_ag, lrr.new_ligand_ag_positioned, cfg.output.inspection_pdb_path)
    print(f"      Inspection structure written: {cfg.output.inspection_pdb_path}")

    merged = lrr.merged_universe
    env_ag = merged.atoms[: lrr.environment_ag_size]
    ins_ag = merged.atoms[lrr.environment_ag_size:]

    # The environment here still contains the (untouched) receptor protein
    # -- unlike receptor-replace mode, where build_environment() already
    # stripped it out. Force "protein" into keep_classes for clash
    # detection so any protein/new-ligand contact is measured and reported
    # (useful QC on the fit/placement) but the receptor itself is never a
    # candidate for removal.
    keep_classes = sorted(set(cfg.clash_detection.keep_classes) | {"protein"})
    clash_cfg = dc_replace(cfg.clash_detection, keep_classes=keep_classes)
    print(f"      (ligand_replace mode: 'protein' added to clash_detection.keep_classes automatically "
          f"-- the untouched receptor is checked/reported for contacts with the new ligand, never removed)")

    return {
        "merged": merged,
        "env_ag": env_ag,
        "ins_ag": ins_ag,
        "removed_ag_orig": lrr.old_ligand_ag,
        "environment_ag_orig": env_ag,
        "clash_cfg": clash_cfg,
        "original_universe_for_report": u,
        "report_kwargs": {"ligand_replace_result": lrr},
    }


def _run_cholesterol(cfg) -> int:
    """Execute the dedicated cholesterol-restoration mode."""
    cs = cfg.cholesterol
    print("=" * 78)
    print("CHOLESTEROL RESTORATION WORKFLOW")
    print("=" * 78)
    print(f"[config] cholesterol.experimental_structure_path = {cs.experimental_structure_path}")
    print(f"[config] cholesterol.target_system_path          = {cs.target_system_path}")
    print(f"[config] alignment.method                        = {cfg.alignment.method}")
    print(
        f"[config] alignment region masks                   = "
        f"{cfg.alignment.region_mask_original} / "
        f"{cfg.alignment.region_mask_replacement}"
    )
    print(f"[config] cholesterol.convert_to_charmm36         = {cs.convert_to_charmm36}")
    print(
        f"[config] cholesterol.missing_heavy_reconstruction = "
        f"{cs.reconstruct_missing_heavy_atoms} (max {cs.max_missing_heavy_atoms}, "
        f"fit RMSD <= {cs.max_heavy_atom_fit_rmsd:.3f} A)"
    )
    print(f"[config] cholesterol.write_diagnostics           = {cs.write_diagnostics}")
    print(f"[config] cholesterol.composition.lipid_targets   = {cs.composition.lipid_targets}")
    print(
        f"[config] cholesterol.composition.reference_conc  = "
        f"{cs.composition.target_concentration} M (assumed from reference; "
        f"tolerance {cs.composition.concentration_tolerance_fraction:.1%})"
    )
    print(f"[config] output.pdb_path / gro_path               = {cfg.output.pdb_path} / {cfg.output.gro_path}")
    print("-" * 78)

    receptor_layout = None
    try:
        if cfg.topology.receptor_toppar_dir:
            receptor_layout = build_receptor_chain_layout(cfg.topology)
        print("[1] Preparing CHARMM36 cholesterol coordinates ...")
        print("[2] Aligning the experimental receptor to the target receptor ...")
        print("[3] Inserting cholesterols and removing clashing molecules ...")
        print("[4] Correcting leaflet lipids and comparing salt with the reference ...")
        target_universe, final_universe, result = restore_cholesterols(
            cfg, workdir=".", receptor_layout=receptor_layout
        )
    except (CholesterolError, AlignmentError, MaskError, TopologyError) as exc:
        print(f"[CHOLESTEROL ERROR] {exc}", file=sys.stderr)
        return 4

    print(
        f"      Inserted {result.n_experimental_cholesterols} cholesterol residue(s) "
        f"({result.n_experimental_cholesterol_atoms} atoms)."
    )
    modeled_heavy = sum(
        len(item["modeled_heavy_atom_names"]) for item in result.heavy_atom_audit
    )
    print(
        f"      Heavy-atom audit: {len(result.heavy_atom_audit)} cholesterol(s) "
        f"checked; {modeled_heavy} missing heavy atom(s) modeled."
    )
    print(
        f"      Alignment RMSD: {result.alignment.rmsd_after:.3f} A; "
        f"clash removals: {result.clash_result.counts_by_class}"
    )
    print(
        f"      Composition removals: {len(result.composition.lipid_residues_removed)} lipid(s), "
        f"{len(result.composition.ions_removed)} ion(s)."
    )
    salt = result.composition
    if salt.estimated_final_concentration_molar is None:
        print("      Salt concentration comparison unavailable: no reference salt pairs.")
    else:
        print(
            f"      Salt pairs reference/final: {salt.reference_salt_pairs}/"
            f"{salt.final_salt_pairs}; estimated final concentration "
            f"{salt.estimated_final_concentration_molar:.4f} M "
            f"({salt.concentration_deviation_fraction:.2%} deviation)."
        )

    topology_result = None
    topology_error = None
    ndx_counts = None
    ndx_notes = None
    ndx_error = None
    refgro_check = None
    partial_failure = False

    if cfg.topology.enabled:
        print("[5] Assembling and auditing GROMACS topology files ...")
        try:
            topology_result = assemble_topology(final_universe, cfg.topology)
            audit = audit_final_topology(final_universe, topology_result)
            print(f"      Wrote {topology_result.top_path}")
            print(f"      Audited topology/coordinate atom order and count "
                  f"({audit['coordinate_atoms']} atoms; net charge "
                  f"{audit['net_charge']:+.6f} e).")
        except TopologyError as exc:
            topology_error = str(exc)
            partial_failure = True
            print(f"[TOPOLOGY ERROR] {exc}", file=sys.stderr)

    print("[6] Writing and validating final coordinates ...")
    try:
        def validate_cholesterol_outputs(staged_pdb, staged_gro):
            geometry = validate_staged_hard_clashes(
                staged_pdb, staged_gro, cfg.box_validation
            )
            cation = cs.composition.salt_cation.upper()
            anion = cs.composition.salt_anion.upper()
            output_counts = {}
            for label, path in (("pdb", staged_pdb), ("gro", staged_gro)):
                written = mda.Universe(path)
                counts = {
                    cation: sum(
                        str(res.resname).upper() == cation
                        for res in written.residues
                    ),
                    anion: sum(
                        str(res.resname).upper() == anion
                        for res in written.residues
                    ),
                }
                if counts != result.composition.ion_counts_after:
                    raise OutputError(
                        f"Staged {label.upper()} changed salt-ion counts: "
                        f"expected {result.composition.ion_counts_after}, found {counts}"
                    )
                output_counts[label] = counts
            result.composition.output_ion_counts = output_counts
            topology_charges = {}
            if topology_result is not None and topology_error is None:
                for label, path in (("pdb", staged_pdb), ("gro", staged_gro)):
                    written = mda.Universe(path)
                    topology_charges[label] = audit_final_topology(
                        written, topology_result
                    )["net_charge"]
                topology_result.charge_audit["coordinate_files_checked"] = {
                    "pdb": cfg.output.pdb_path,
                    "gro": cfg.output.gro_path,
                }
                topology_result.charge_audit[
                    "coordinate_file_net_charges"
                ] = topology_charges
            return {
                "periodic_box_and_hard_clashes": geometry,
                "salt_ion_counts": output_counts,
                "topology_net_charges": topology_charges,
            }

        write_outputs(
            final_universe,
            cfg.output.pdb_path,
            cfg.output.gro_path,
            protein_chain_atom_counts=(
                receptor_layout.atom_counts if receptor_layout else None
            ),
            protected_input_paths=[
                cs.experimental_structure_path,
                cs.target_system_path,
                cs.composition.reference_system_path,
            ],
            staged_validator=validate_cholesterol_outputs,
        )
    except (OutputError, BoxValidationError, TopologyError, OSError, ValueError) as exc:
        print(f"[OUTPUT VALIDATION ERROR] {exc}", file=sys.stderr)
        return 6
    print(f"      Wrote {cfg.output.pdb_path}")
    print(f"      Wrote {cfg.output.gro_path}")

    if cfg.ndx.enabled:
        print("[7] Writing GROMACS index groups ...")
        try:
            ndx_path, ndx_counts, ndx_notes = write_ndx(final_universe, cfg.ndx)
            print(f"      Wrote {ndx_path} ({ndx_counts})")
        except (MaskError, ValueError) as exc:
            ndx_error = str(exc)
            partial_failure = True
            print(f"[NDX ERROR] {exc}", file=sys.stderr)

    if cfg.refgro:
        print("[8] Checking final species against the reference GRO ...")
        refgro_check = refgro_completeness_check(final_universe, cfg.refgro)
        if refgro_check["missing_from_final"]:
            print(
                "      Species present in refgro but absent from the final system: "
                f"{refgro_check['missing_from_final']}"
            )

    report_data = build_cholesterol_report(
        cfg,
        result,
        final_universe,
        target_universe,
        topology_result=topology_result,
        topology_error=topology_error,
        ndx_counts=ndx_counts,
        ndx_notes=ndx_notes,
        ndx_error=ndx_error,
        refgro_check=refgro_check,
    )
    txt_path, json_path = write_cholesterol_report(
        report_data, cfg.output.report_path
    )
    print(f"      Wrote {txt_path}")
    print(f"      Wrote {json_path}")

    print("\n" + "=" * 78)
    print("Final files:")
    print(f"  {cfg.output.pdb_path}")
    print(f"  {cfg.output.gro_path}")
    print(f"  {txt_path}")
    print(f"  {json_path}")
    if cfg.cholesterol.write_diagnostics:
        print(f"  {cfg.output.inspection_pdb_path}")
        if result.merged_pdb_path:
            print(f"  {result.merged_pdb_path}")
    for path in result.diagnostic_paths:
        print(f"  {path}")
    if topology_result is not None:
        print(f"  {topology_result.top_path}")
        print(f"  {topology_result.toppar_dir}/")
    if ndx_counts is not None:
        print(f"  {cfg.ndx.output_path}")
    print("=" * 78)
    if result.composition.warnings:
        print("FINAL WARNINGS:")
        for warning in result.composition.warnings:
            print(f"  WARNING: {warning}")
    return 6 if partial_failure else 0


def main(config_path: str, mode: str, dry_run: bool = False) -> int:
    try:
        cfg = load_config(config_path, mode=mode, check_paths=not dry_run)
    except ConfigError as e:
        print(f"[CONFIG ERROR] {e}", file=sys.stderr)
        return 2

    if dry_run:
        print("Configuration dry-run passed (schema and values only; input paths were not opened).")
        print(f"  mode: {mode}")
        print(f"  config: {config_path}")
        print(f"  final PDB/GRO: {cfg.output.pdb_path} / {cfg.output.gro_path}")
        print(f"  topology: {cfg.topology.output_dir}/topol.top")
        print(f"  index: {cfg.ndx.output_path}")
        print(f"  reports: {cfg.output.report_path}.txt / {cfg.output.report_path}.json")
        return 0

    if mode == "chl":
        return _run_cholesterol(cfg)

    lr_mode = mode == "lig"
    resolved_replacement_ligands = []
    if not lr_mode:
        try:
            resolved_replacement_ligands = resolve_replacement_ligands(cfg)
        except ReplacementLigandError as exc:
            print(f"[REPLACEMENT LIGAND ERROR] {exc}", file=sys.stderr)
            return 2

    print("=" * 78)
    print("LIGAND-ONLY REPLACEMENT WORKFLOW" if lr_mode else "RECEPTOR REPLACEMENT WORKFLOW")
    print("=" * 78)
    if lr_mode:
        print(f"[config] ligand_replace.structure_path      = {cfg.ligand_replace.structure_path}")
        print(f"[config] ligand_replace.original_ligand      = {cfg.ligand_replace.original_ligand.resname}")
        print(f"[config] ligand_replace.new_ligand.coord_path = {cfg.ligand_replace.new_ligand.coord_path}")
        print(f"[config] ligand_replace.new_ligand.resname    = {cfg.ligand_replace.new_ligand.resname or '(kept from coord_path)'}")
        print(f"[config] ligand_replace.fit.method            = {cfg.ligand_replace.fit.method}")
        if cfg.ligand_replace.fit.method == "pairfit":
            print(f"[config] ligand_replace.fit atom pairs        = "
                  f"{list(zip(cfg.ligand_replace.fit.old_ligand_fit_atoms, cfg.ligand_replace.fit.new_ligand_fit_atoms))}")
        elif cfg.ligand_replace.fit.method == "autofit":
            print("[config] ligand_replace.fit atoms             = all uniquely named heavy atoms (matched by name)")
    else:
        print(f"[config] target_box.path              = {cfg.target_box.path}")
        print(f"[config] target_box.receptor_mask      = {cfg.target_box.receptor_mask}")
        print(f"[config] replacement_structure.path    = {cfg.replacement_structure.path}")
        print(f"[config] replacement_structure.mask    = {cfg.replacement_structure.receptor_mask}")
        print(f"[config] alignment.method              = {cfg.alignment.method}")
        print(
            f"[config] alignment region masks         = "
            f"{cfg.alignment.region_mask_original} / "
            f"{cfg.alignment.region_mask_replacement}"
        )
        print(
            f"[config] name_restoration               = "
            f"{cfg.name_restoration.enabled} ({cfg.name_restoration.method})"
        )
        if (
            cfg.name_restoration.enabled
            and cfg.name_restoration.method == "itp_atom_count"
        ):
            print(
                f"[config] name_restoration.reference_topol = "
                f"{cfg.name_restoration.reference_topol}"
            )
            print(
                "[config] pdb_to_full_resname             = "
                f"{cfg.name_restoration.pdb_to_full_resname}"
            )
            print(
                "[config] reference topology counts/order = ignored"
            )
        if resolved_replacement_ligands:
            print("[config] replacement_ligands:")
            for ligand in resolved_replacement_ligands:
                print(
                    f"           {ligand.selector}: charge={ligand.net_charge:+.6f} e "
                    f"({ligand.charge_source}), moleculetype={ligand.moleculetype}, "
                    f"itp={ligand.itp_path}"
                )
        else:
            print("[config] replacement_ligands           = []")
    print(f"[config] clash_detection.thresholds    = {{'default': {cfg.clash_detection.threshold}, "
          f"**{cfg.clash_detection.thresholds}}}")
    if cfg.clash_detection.keep_classes:
        print(f"[config] clash_detection.keep_classes  = {cfg.clash_detection.keep_classes} "
              f"(reported, never removed)")
    print(f"[config] charge.target_net_charge      = {cfg.charge.target_net_charge}   "
          f"neutralize={cfg.charge.neutralize}")
    if cfg.charge.exclude_membrane_interior:
        print(
            "[config] charge lipid protection        = PBC distance to lipid heavy atoms "
            f">= {cfg.charge.exclusion_distance_from_lipid} A (no z slab)"
        )
    if cfg.charge.charge_table_overrides:
        print(f"[config] charge.charge_table_overrides = {cfg.charge.charge_table_overrides}")
    print(f"[config] topology.enabled              = {cfg.topology.enabled}"
          + (f"  -> output_dir={cfg.topology.output_dir}" if cfg.topology.enabled else ""))
    print(f"[config] ndx.enabled                   = {cfg.ndx.enabled}"
          + (f"  -> output_path={cfg.ndx.output_path}" if cfg.ndx.enabled else ""))
    print(f"[config] output.pdb_path / gro_path    = {cfg.output.pdb_path} / {cfg.output.gro_path}")
    print(f"[config] output.report_path            = {cfg.output.report_path} (.txt + .json)")
    print("-" * 78)

    replacement_ligand_itp_paths = []

    # topology.ligand_itp_paths automatically picks up ligand_replace's own
    # new_ligand.itp_path, if set, so it doesn't need to be listed twice.
    if lr_mode and cfg.ligand_replace.new_ligand.itp_path:
        replacement_ligand_itp_paths.append(
            cfg.ligand_replace.new_ligand.itp_path
        )
        if cfg.ligand_replace.new_ligand.itp_path not in cfg.topology.ligand_itp_paths:
            cfg.topology.ligand_itp_paths.append(cfg.ligand_replace.new_ligand.itp_path)
            print(f"      (added ligand_replace.new_ligand.itp_path to topology.ligand_itp_paths automatically)")

    try:
        if lr_mode:
            prep = _prepare_ligand_replace(cfg)
        else:
            prep = _prepare_receptor_replace(cfg, resolved_replacement_ligands)
    except _StageFailure as sf:
        return sf.exit_code

    merged = prep["merged"]
    env_ag = prep["env_ag"]
    ins_ag = prep["ins_ag"]
    removed_ag_orig = prep["removed_ag_orig"]
    environment_ag_orig = prep["environment_ag_orig"]
    clash_cfg = prep["clash_cfg"]
    original_universe_for_report = prep["original_universe_for_report"]
    report_kwargs = prep["report_kwargs"]

    print("[10-11] Detecting and removing steric clashes (PBC-aware) ...")
    try:
        clash_res = detect_and_remove_clashes(
            env_ag, ins_ag, merged.dimensions, clash_cfg
        )
    except (ValueError, RuntimeError) as exc:
        print(f"[CLASH ERROR] {exc}", file=sys.stderr)
        return 4
    print(f"      Removed by class: {clash_res.counts_by_class}  "
          f"(lipids {clash_res.n_lipids_removed}/{clash_res.n_lipids_total}, "
          f"{clash_res.lipid_removal_fraction:.1%})")
    print(f"      Removed by resname: {clash_res.counts_by_resname}")
    if clash_res.lipid_removal_flagged:
        print("      *** FLAGGED: unusually large lipid removal -- recommend manual inspection ***")

    print("[12-13] Computing net charge and correcting via counterion removal ...")
    topology_charge_model = None
    receptor_layout = None
    if cfg.topology.receptor_toppar_dir and cfg.topology.environment_toppar_dir:
        try:
            receptor_layout = build_receptor_chain_layout(cfg.topology)
            if lr_mode:
                env_nonprotein = environment_ag_orig[receptor_layout.total_atoms:]
                ins_nonprotein = ins_ag
            else:
                env_nonprotein = environment_ag_orig
                if len(ins_ag) < receptor_layout.total_atoms:
                    raise TopologyError(
                        "The inserted receptor block is shorter than the receptor "
                        "described by the reference topology."
                    )
                ins_nonprotein = ins_ag[receptor_layout.total_atoms:]
            nonprotein_resnames = set(env_nonprotein.residues.resnames) | set(
                ins_nonprotein.residues.resnames
            )
            topology_charge_model = build_topology_charge_model(
                cfg.topology,
                nonprotein_resnames,
                replacement_ligand_itp_paths=replacement_ligand_itp_paths,
            )
            print(
                "      ITP receptor chains: "
                + ", ".join(
                    f"{mtype}={charge:+.3f} e"
                    for mtype, charge in zip(
                        receptor_layout.moleculetypes, receptor_layout.charges
                    )
                )
            )
        except TopologyError as exc:
            print(f"[CHARGE TOPOLOGY ERROR] {exc}", file=sys.stderr)
            return 5
    try:
        charge_report = compute_pre_neutralization_report(
            removed_ag_orig,
            ins_ag,
            environment_ag_orig,
            clash_res.removed_ag,
            cfg.charge,
            topology_charge_model=topology_charge_model,
            inserted_has_receptor=not lr_mode,
            environment_has_receptor=lr_mode,
            original_charge_from_itp=lr_mode,
            resolve_original_charge=not lr_mode,
        )
        charge_report = neutralize(
            clash_res.kept_ag, ins_ag, charge_report, cfg.charge, merged
        )
    except ChargeError as exc:
        print(f"[CHARGE ERROR] {exc}", file=sys.stderr)
        return 5

    final_universe = None
    topology_result = None
    topology_error = None
    ndx_counts = None
    ndx_notes = None
    ndx_error = None
    refgro_check = None
    output_validation = None
    partial_failure = False
    coordinates_written = False
    if charge_report.stopped_reason:
        print(f"[STOPPED] {charge_report.stopped_reason}", file=sys.stderr)
    else:
        print(f"      Net charge before/after clash removal: "
              f"{charge_report.net_charge_before_clash_removal:.3f} / "
              f"{charge_report.net_charge_after_clash_removal:.3f} e")
        if charge_report.neutralization_performed and charge_report.ions_removed:
            print(f"      Removed {len(charge_report.ions_removed)} counterion(s) -> "
                  f"final net charge {charge_report.final_net_charge:.3f} e")
        if charge_report.ion_selection_summary:
            print(
                "      Ion selection: "
                f"{charge_report.ion_selection_summary}"
            )

        print("[14] Building and validating the final molecular composition ...")
        final_universe = finalize_system(
            clash_res.kept_ag,
            ins_ag,
            merged.dimensions,
            charge_report.ions_removed,
            protein_chain_atom_counts=(
                receptor_layout.atom_counts if receptor_layout else None
            ),
            inserted_contains_receptor=not lr_mode,
        )
        if cfg.topology.enabled:
            print("[15] Assembling and auditing GROMACS topol.top + toppar/ ...")
            try:
                topology_result = assemble_topology(
                    final_universe,
                    cfg.topology,
                    replacement_ligand_itp_paths=replacement_ligand_itp_paths,
                )
                audit = audit_final_topology(final_universe, topology_result)
                if charge_report.final_net_charge is None or abs(
                    audit["net_charge"] - charge_report.final_net_charge
                ) > cfg.charge.tolerance:
                    raise TopologyError(
                        "Final topology charge does not match pre-output charge "
                        f"accounting: generated topology={audit['net_charge']:+.6f} e, "
                        f"accounting={charge_report.final_net_charge!r} e."
                    )
                charge_report.notes.append(
                    "Final charge was independently confirmed from generated "
                    "topol.top molecule counts and copied ITP [ atoms ] sums."
                )
                print(f"      Wrote {topology_result.top_path}")
                print(f"      Built {topology_result.toppar_dir}/ with "
                      f"{len(topology_result.includes_written)} file(s):")
                for inc in topology_result.includes_written:
                    print(f"        {inc}")
                print(f"      [ molecules ]:")
                for mtype, count in topology_result.molecules_written:
                    print(f"        {mtype:<16s} {count}")
                for note in topology_result.notes:
                    print(f"      NOTE: {note}")
                print(f"      Audited final topology charge: {audit['net_charge']:+.6f} e")
            except TopologyError as e:
                topology_error = str(e)
                partial_failure = True
                print(f"[TOPOLOGY ERROR] {e}", file=sys.stderr)

        if topology_error is None:
            protected_inputs = (
                [
                    cfg.ligand_replace.structure_path,
                    cfg.ligand_replace.new_ligand.coord_path,
                ]
                if lr_mode else
                [cfg.target_box.path, cfg.replacement_structure.path]
            )
            print("[16] Writing final PDB/GRO transactionally and reading both back ...")
            try:
                def validate_staged_coordinates(staged_pdb, staged_gro):
                    staged_geometry = validate_staged_hard_clashes(
                        staged_pdb, staged_gro, cfg.box_validation, clash_cfg
                    )
                    file_charges = {}
                    if topology_result is not None:
                        for label, path in (
                            ("final_pdb", staged_pdb),
                            ("final_gro", staged_gro),
                        ):
                            written = mda.Universe(path)
                            file_charges[label] = audit_final_topology(
                                written, topology_result
                            )["net_charge"]
                        topology_result.charge_audit["coordinate_files_checked"] = {
                            "pdb": cfg.output.pdb_path,
                            "gro": cfg.output.gro_path,
                        }
                        topology_result.charge_audit["coordinate_file_net_charges"] = file_charges
                    return {
                        "periodic_box_and_hard_clashes": staged_geometry,
                        "topology_net_charges": file_charges,
                    }

                output_validation = write_outputs(
                    final_universe,
                    cfg.output.pdb_path,
                    cfg.output.gro_path,
                    protein_chain_atom_counts=(
                        receptor_layout.atom_counts if receptor_layout else None
                    ),
                    protected_input_paths=protected_inputs,
                    staged_validator=validate_staged_coordinates,
                )
                coordinates_written = True
                print(f"      Wrote {cfg.output.pdb_path}")
                print(f"      Wrote {cfg.output.gro_path}")
                if output_validation["input_path_collisions_safely_replaced"]:
                    print(
                        "      Input/output collision handled through validated temporary "
                        "files and atomic rename: "
                        f"{output_validation['input_path_collisions_safely_replaced']}"
                    )
            except (OutputError, TopologyError, OSError, ValueError) as e:
                partial_failure = True
                print(f"[OUTPUT VALIDATION ERROR] {e}", file=sys.stderr)

        if coordinates_written and cfg.ndx.enabled:
            print("[15] Writing index.ndx ...")
            try:
                ndx_path, ndx_counts, ndx_notes = write_ndx(final_universe, cfg.ndx)
                print(f"      Wrote {ndx_path}  ({ndx_counts})")
            except MaskError as e:
                ndx_error = str(e)
                partial_failure = True
                print(f"[NDX ERROR] {e}", file=sys.stderr)

        if coordinates_written and cfg.refgro:
            print("[16] Checking final system against refgro (non-protein completeness) ...")
            refgro_check = refgro_completeness_check(final_universe, cfg.refgro)
            if refgro_check["missing_from_final"]:
                print(f"      NOTE: present in refgro but absent from final system: "
                      f"{refgro_check['missing_from_final']}")
            else:
                print("      All of refgro's non-protein resnames are present in the final system.")

    report_dict = build_report_dict(
        cfg, report_kwargs.get("align_res"), report_kwargs.get("rep_res"),
        clash_res, charge_report,
        final_universe if coordinates_written else None,
        original_universe_for_report,
        namefix_summary=report_kwargs.get("namefix_summary"),
        topology_result=topology_result, topology_error=topology_error,
        ndx_counts=ndx_counts, ndx_notes=ndx_notes, ndx_error=ndx_error,
        refgro_check=refgro_check,
        ligand_replace_result=report_kwargs.get("ligand_replace_result"),
        resolved_replacement_ligands=report_kwargs.get("resolved_replacement_ligands"),
        output_validation=output_validation,
    )
    txt_path, json_path = write_report(report_dict, cfg.output.report_path)
    print(f"      Wrote {txt_path}")
    print(f"      Wrote {json_path}")

    if charge_report.stopped_reason:
        print("\nPipeline stopped before producing final output (see [STOPPED] message above "
              "and the report for full context). No PDB/GRO was written.")
        return 5

    print("\n" + "=" * 78)
    print("Final files:")
    if coordinates_written:
        print(f"  {cfg.output.pdb_path}")
        print(f"  {cfg.output.gro_path}")
    print(f"  {cfg.output.inspection_pdb_path}")
    print(f"  {txt_path}")
    print(f"  {json_path}")
    if topology_result is not None:
        print(f"  {topology_result.top_path}")
        print(f"  {topology_result.toppar_dir}/  ({len(topology_result.includes_written)} files)")
    if ndx_counts is not None:
        print(f"  {cfg.ndx.output_path}")
    print("=" * 78)

    if partial_failure:
        if coordinates_written:
            print("\nCoordinate output was validated, but a later optional output failed; "
                  "see the reported error above.")
        else:
            print("\nFinal coordinate output was not promoted because topology or output "
                  "validation failed; see the reported error above.")
        return 6

    print("\nDone.")
    return 0


def _build_cli_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="gmxtransplant",
        description="Modify a prepared molecular-simulation system."
    )
    parser.add_argument(
        "--mode", choices=("receptor", "lig", "chl"),
        help="receptor replacement, ligand replacement, or cholesterol restoration",
    )
    parser.add_argument(
        "-i", "--input", metavar="CONFIG.yaml",
        help="mode-specific YAML configuration",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="validate configuration schema/values without opening configured input paths",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument(
        "--show-example",
        choices=("receptor", "lig", "chl"),
        metavar="MODE",
        help="print a bundled mode-specific YAML example and exit",
    )
    return parser


def cli(argv=None) -> int:
    """Console-script entry point used by both Python and installed packages."""
    parser = _build_cli_parser()
    cli_args = parser.parse_args(argv)
    if cli_args.show_example:
        try:
            sys.stdout.write(_read_example(cli_args.show_example))
        except ConfigError as exc:
            parser.error(str(exc))
        return 0
    if not cli_args.mode or not cli_args.input:
        parser.error("--mode and --input are required unless --show-example is used")
    try:
        return main(cli_args.input, cli_args.mode, dry_run=cli_args.dry_run)
    except Exception:
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(cli())
