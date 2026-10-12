"""Mode-specific YAML configuration for the system-modification pipeline.

The definitions live in config_base (errors, YAML, paths) and config_schema (the section
dataclasses); this module loads and validates a YAML file and re-exports both, so
`from config import X` works for every name.
"""

from __future__ import annotations

import os
import math
import re
import sys
import sysconfig
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Union

import yaml
from yaml.nodes import MappingNode

from mode_registry import EDIT_MODES, path_kind as edit_path_kind

from config_base import *  # noqa: F401,F403
from config_base import _require_mapping, _UNANCHORED_OUTPUTS, _UniqueKeySafeLoader, _construct_unique_mapping, _is_finite_number, _package_roots, _resolve_path_references  # noqa: F401
from config_schema import *  # noqa: F401,F403
from config_schema import _default_repack_stages, _repack_errors  # noqa: F401


def _get(d, key, default=None):
    return d.get(key, default) if d else default


_COMMON_SECTIONS = {
    "paths", "alignment", "clash_detection", "topology", "ndx", "refgro",
    "output", "box_validation", "minimization"
}
_MODE_SECTIONS = {
    "protein": _COMMON_SECTIONS | {
        "target_box", "replacement_structure", "name_restoration",
        "replacement_ligands", "charge",
    },
    "lig": _COMMON_SECTIONS | {"ligand_replace", "charge"},
    "chl": _COMMON_SECTIONS | {"cholesterol"},
}

_SECTION_KEYS = {
    "box_validation": set(BoxValidationSpec.__dataclass_fields__),
    "target_box": {"path", "format", "protein_mask", "box_dimensions"},
    "replacement_structure": {"path", "format", "protein_mask", "box_dimensions"},
    # atoms is retained only as a deprecated compatibility fallback. New
    # configurations put @atom names directly in each region mask.
    "alignment": {
        "method", "atoms", "region_mask_original", "region_mask_replacement",
        "cycles", "cutoff",
    },
    "clash_detection": {
        "threshold", "thresholds", "use_pbc", "heavy_atoms_only",
        "flag_lipid_removal_fraction", "keep_classes",
    },
    "charge": {
        "target_net_charge", "neutralize", "tolerance",
        "exclusion_distance_from_protein", "exclusion_distance_from_ligand",
        "exclusion_distance_from_lipid", "membrane_z_margin", "random_seed", "charge_table_overrides",
        "exclude_membrane_interior", "membrane_z_override",
    },
    "output": {"pdb_path", "gro_path", "report_path", "inspection_pdb_path"},
    "minimization": set(MinimizationSpec.__dataclass_fields__),
    "name_restoration": {
        "enabled", "method", "reference", "apply_to", "min_jaccard",
        "reference_topol", "pdb_to_full_resname",
        # Accepted legacy names, converted explicitly below.
        "topology_path", "topology_to_pdb_resname",
    },
    "topology": {
        "enabled", "protein_toppar_dir", "protein_template_top",
        "environment_toppar_dir", "ligand_itp_paths",
        "moleculetype_overrides", "output_dir",
    },
    "ndx": {"enabled", "output_path", "style", "extra_groups"},
    "ligand_replace": {
        "structure_path", "format", "box_dimensions", "original_ligand",
        "new_ligand", "fit",
    },
    "cholesterol": {
        "experimental_structure_path", "target_system_path",
        "target_system_format", "box_dimensions", "target_protein_mask",
        "experimental_protein_mask", "cholesterol_resnames", "charmm_resname",
        "convert_to_charmm36", "reconstruct_missing_heavy_atoms",
        "max_missing_heavy_atoms", "max_heavy_atom_fit_rmsd", "write_diagnostics",
        "charmm36_reference_path", "obabel_command",
        "verify_tolerance", "raw_protonated_pdb_path", "converted_pdb_path",
        "aligned_reference_pdb_path", "reference_placed_pdb_path",
        "merged_pdb_path", "composition", "repack",
    },
}



def _reject_unknown_keys(mapping: dict, label: str, allowed) -> None:
    unexpected = sorted(set(mapping) - set(allowed))
    if unexpected:
        raise ConfigError(
            f"Unknown key(s) in {label}: {', '.join(unexpected)}"
        )


def _validate_nested_keys(raw: dict, mode: str) -> None:
    for section in sorted(_MODE_SECTIONS[mode] - {"paths", "refgro", "replacement_ligands"}):
        if section not in raw:
            continue
        mapping = _require_mapping(raw[section], section)
        _reject_unknown_keys(mapping, section, _SECTION_KEYS[section])

    if "paths" in raw:
        paths = _require_mapping(raw["paths"], "paths")
        for key, value in paths.items():
            if not isinstance(key, str) or not re.fullmatch(r"[a-z][a-z0-9_]*", key):
                raise ConfigError(
                    "paths keys must use lower-case snake_case names"
                )
            if value is not None and not isinstance(value, str):
                raise ConfigError(f"paths.{key} must be a path string or null")

    mapping_fields = (
        ("clash_detection", "thresholds"),
        ("charge", "charge_table_overrides"),
        ("topology", "moleculetype_overrides"),
        ("ndx", "extra_groups"),
        ("minimization", "defines"),
    )
    for section, key in mapping_fields:
        value = (raw.get(section) or {}).get(key)
        if value is not None and not isinstance(value, dict):
            raise ConfigError(f"{section}.{key} must be a mapping")
    list_fields = (
        ("clash_detection", "keep_classes"),
        ("topology", "ligand_itp_paths"),
    )
    for section, key in list_fields:
        value = (raw.get(section) or {}).get(key)
        if value is not None and not isinstance(value, list):
            raise ConfigError(f"{section}.{key} must be a list")

    raw_ligands = raw.get("replacement_ligands", []) or []
    if not isinstance(raw_ligands, list):
        raise ConfigError("replacement_ligands must be a list")
    for index, item in enumerate(raw_ligands, 1):
        mapping = _require_mapping(item, f"replacement_ligands[{index}]")
        _reject_unknown_keys(
            mapping, f"replacement_ligands[{index}]",
            {"resname", "itp_path", "charge", "moleculetype", "forcefield_path"},
        )

    if "ligand_replace" in raw:
        lr = _require_mapping(raw["ligand_replace"], "ligand_replace")
        for key, allowed in {
            "original_ligand": {"resname"},
            "new_ligand": {"coord_path", "format", "resname", "itp_path", "forcefield_path"},
            "fit": {"method", "old_ligand_fit_atoms", "new_ligand_fit_atoms"},
        }.items():
            if key in lr:
                nested = _require_mapping(lr[key], f"ligand_replace.{key}")
                _reject_unknown_keys(nested, f"ligand_replace.{key}", allowed)
        fit = lr.get("fit") or {}
        for key in ("old_ligand_fit_atoms", "new_ligand_fit_atoms"):
            value = fit.get(key)
            if value is not None and not isinstance(value, list):
                raise ConfigError(f"ligand_replace.fit.{key} must be a list")

    if "cholesterol" in raw:
        cholesterol = _require_mapping(raw["cholesterol"], "cholesterol")
        value = cholesterol.get("cholesterol_resnames")
        if value is not None and not isinstance(value, list):
            raise ConfigError("cholesterol.cholesterol_resnames must be a list")
        if "repack" in cholesterol:
            repack = _require_mapping(cholesterol["repack"], "cholesterol.repack")
            _reject_unknown_keys(repack, "cholesterol.repack", set(RepackSpec.__dataclass_fields__))
            stages = repack.get("stages")
            if stages is not None:
                if not isinstance(stages, list) or not stages:
                    raise ConfigError("cholesterol.repack.stages must be a nonempty list")
                for number, stage in enumerate(stages):
                    _reject_unknown_keys(_require_mapping(stage, f"cholesterol.repack.stages[{number}]"),
                                         f"cholesterol.repack.stages[{number}]",
                                         set(RepackStageSpec.__dataclass_fields__))
        if "composition" in cholesterol:
            composition = _require_mapping(
                cholesterol["composition"], "cholesterol.composition"
            )
            _reject_unknown_keys(
                composition,
                "cholesterol.composition",
                {
                    "reference_system_path", "distance_from_protein",
                    "target_concentration", "concentration_tolerance_fraction",
                    "salt_cation", "salt_anion",
                    "random_seed", "lipid_targets",
                },
            )
            targets = composition.get("lipid_targets", {}) or {}
            if not isinstance(targets, dict):
                raise ConfigError(
                    "cholesterol.composition.lipid_targets must be a mapping"
                )
            for resname, leaves in targets.items():
                leaf_map = _require_mapping(
                    leaves,
                    f"cholesterol.composition.lipid_targets.{resname}",
                )
                _reject_unknown_keys(
                    leaf_map,
                    f"cholesterol.composition.lipid_targets.{resname}",
                    {"upper", "lower"},
                )


def load_config(path: str, mode: str = "protein", check_paths: bool = True,
                output_root=None) -> Config:
    if mode == "charmprot":
        from charmprot import load_charmprot_config
        return load_charmprot_config(path, check_paths=check_paths, output_root=output_root)
    if mode not in _MODE_SECTIONS:
        raise ConfigError("Mode must be one of: protein, lig, chl")
    try:
        with open(path, encoding="utf-8") as fh:
            raw = yaml.load(fh, Loader=_UniqueKeySafeLoader)
    except ConfigError:
        raise
    except yaml.YAMLError as exc:
        raise ConfigError(f"Cannot parse YAML config '{path}': {exc}") from exc
    except OSError as exc:
        raise ConfigError(f"Cannot read config file '{path}': {exc}") from exc
    if not raw:
        raise ConfigError(f"Config file '{path}' is empty")
    if not isinstance(raw, dict):
        raise ConfigError(f"Config file '{path}' must contain a top-level mapping")
    unexpected = sorted(set(raw) - _MODE_SECTIONS[mode])
    if unexpected:
        raise ConfigError(
            f"The following top-level section(s) do not belong in {mode} mode: "
            + ", ".join(unexpected)
        )

    raw = _resolve_path_references(raw)
    raw = anchor_config_paths(raw, config_path=path, output_root=output_root)
    _validate_nested_keys(raw, mode)

    ob = raw.get("target_box", {})
    rs = raw.get("replacement_structure", {})
    al = raw.get("alignment", {})
    raw_replacement_ligands = raw.get("replacement_ligands", []) or []
    if not isinstance(raw_replacement_ligands, list):
        raise ConfigError("replacement_ligands must be a list")
    if any(not isinstance(item, dict) for item in raw_replacement_ligands):
        raise ConfigError("every replacement_ligands entry must be a mapping")
    cd = raw.get("clash_detection", {})
    ch = raw.get("charge", {})
    out = raw.get("output", {})
    nr = raw.get("name_restoration", {})
    raw_namefix_method = str(_get(nr, "method", "atom_signature")).strip()
    # Backward compatibility for the two earlier names of the topology-based
    # restoration method.
    namefix_method = (
        "itp_atom_count"
        if raw_namefix_method in ("topology_order", "topology_alias")
        else raw_namefix_method
    )
    raw_pdb_to_full = _get(nr, "pdb_to_full_resname", {}) or {}
    raw_topology_to_pdb = _get(nr, "topology_to_pdb_resname", {}) or {}
    if raw_pdb_to_full and raw_topology_to_pdb:
        raise ConfigError(
            "Use name_restoration.pdb_to_full_resname only; do not also provide "
            "the legacy topology_to_pdb_resname mapping"
        )
    if not isinstance(raw_pdb_to_full, dict):
        raise ConfigError("name_restoration.pdb_to_full_resname must be a mapping")
    if not isinstance(raw_topology_to_pdb, dict):
        raise ConfigError("name_restoration.topology_to_pdb_resname must be a mapping")

    pdb_to_full_resname: Dict[str, List[str]] = {}
    if raw_pdb_to_full:
        for pdb_name, raw_candidates in raw_pdb_to_full.items():
            if isinstance(raw_candidates, str):
                candidates = [raw_candidates]
            elif isinstance(raw_candidates, list):
                candidates = raw_candidates
            else:
                raise ConfigError(
                    "Each name_restoration.pdb_to_full_resname value must be a "
                    "full-name string or a list of full-name strings"
                )
            if any(not isinstance(name, str) for name in candidates):
                raise ConfigError(
                    "Every name_restoration.pdb_to_full_resname candidate must be a string"
                )
            pdb_to_full_resname[str(pdb_name).strip()] = [
                name.strip() for name in candidates
            ]
    else:
        # Convert the legacy full-name-to-PDB-alias direction into the new
        # PDB-alias-to-full-name list direction.
        for full_name, pdb_name in raw_topology_to_pdb.items():
            alias = str(pdb_name).strip()
            pdb_to_full_resname.setdefault(alias, []).append(str(full_name).strip())

    new_reference_topol = str(_get(nr, "reference_topol", "")).strip()
    legacy_topology_path = str(_get(nr, "topology_path", "")).strip()
    if new_reference_topol and legacy_topology_path and new_reference_topol != legacy_topology_path:
        raise ConfigError(
            "name_restoration.reference_topol and legacy topology_path refer to "
            "different files; keep only reference_topol"
        )
    reference_topol = new_reference_topol or legacy_topology_path
    tp = raw.get("topology", {})
    if not reference_topol and tp.get("environment_toppar_dir"):
        reference_topol = os.path.join(
            os.path.dirname(os.path.normpath(tp["environment_toppar_dir"])), "topol.top"
        )
    ix = raw.get("ndx", {})
    mn = raw.get("minimization", {})
    lr = raw.get("ligand_replace", {})
    lr_ol = _get(lr, "original_ligand", {}) or {}
    lr_nl = _get(lr, "new_ligand", {}) or {}
    lr_fit = _get(lr, "fit", {}) or {}
    cs = raw.get("cholesterol", {})
    cs_comp = _get(cs, "composition", {}) or {}
    cs_repack = _get(cs, "repack", {}) or {}
    repack_defaults = RepackSpec()
    raw_stages = cs_repack.get("stages")
    repack_stages = (repack_defaults.stages if raw_stages is None else
                     [RepackStageSpec(**stage) for stage in raw_stages])
    raw_targets = _get(cs_comp, "lipid_targets", {}) or {}
    if not isinstance(raw_targets, dict):
        raise ConfigError("cholesterol.composition.lipid_targets must be a mapping")
    lipid_targets = {
        str(resname).upper(): {
            str(leaf).lower(): count for leaf, count in (leaves or {}).items()
        }
        for resname, leaves in raw_targets.items()
        if isinstance(leaves, dict)
    }

    cfg = Config(
        mode=mode,
        box_validation=BoxValidationSpec(**(raw.get("box_validation") or {})),
        target_box=StructureSpec(
            path=_get(ob, "path", ""),
            format=_get(ob, "format", "auto"),
            protein_mask=_get(ob, "protein_mask", ""),
            box_dimensions=_get(ob, "box_dimensions", None),
        ),
        replacement_structure=StructureSpec(
            path=_get(rs, "path", ""),
            format=_get(rs, "format", "auto"),
            protein_mask=_get(rs, "protein_mask", ""),
            box_dimensions=_get(rs, "box_dimensions", None),
        ),
        alignment=AlignmentSpec(
            method=_get(al, "method", "pymol_align"),
            atoms=_get(al, "atoms", None),
            region_mask_original=_get(al, "region_mask_original", None),
            region_mask_replacement=_get(al, "region_mask_replacement", None),
            cycles=_get(al, "cycles", 5),
            cutoff=_get(al, "cutoff", 2.0),
        ),
        replacement_ligands=[
            ReplacementLigandSpec(
                resname=str(_get(item, "resname", "")).strip().upper(),
                itp_path=str(_get(item, "itp_path", "")).strip(),
                charge=_get(item, "charge", "from_itp"),
                moleculetype=str(_get(item, "moleculetype", "")).strip(),
                forcefield_path=str(_get(item, "forcefield_path", "") or "").strip(),
            )
            for item in raw_replacement_ligands
        ],
        clash_detection=ClashDetectionSpec(
            threshold=_get(cd, "threshold", 3.0),
            thresholds=_get(cd, "thresholds", {}) or {},
            use_pbc=_get(cd, "use_pbc", True),
            heavy_atoms_only=_get(cd, "heavy_atoms_only", True),
            flag_lipid_removal_fraction=_get(cd, "flag_lipid_removal_fraction", 0.05),
            keep_classes=_get(cd, "keep_classes", []) or [],
        ),
        charge=ChargeSpec(
            target_net_charge=_get(ch, "target_net_charge", 0.0),
            neutralize=_get(ch, "neutralize", True),
            tolerance=_get(ch, "tolerance", 0.01),
            exclusion_distance_from_protein=_get(ch, "exclusion_distance_from_protein", 10.0),
            exclusion_distance_from_ligand=_get(ch, "exclusion_distance_from_ligand", 10.0),
            exclusion_distance_from_lipid=_get(
                ch, "exclusion_distance_from_lipid", _get(ch, "membrane_z_margin", 5.0)
            ),
            membrane_z_margin=_get(ch, "membrane_z_margin", None),
            random_seed=_get(ch, "random_seed", 42),
            charge_table_overrides=_get(ch, "charge_table_overrides", {}) or {},
            exclude_membrane_interior=_get(ch, "exclude_membrane_interior", True),
            membrane_z_override=_get(ch, "membrane_z_override", None),
        ),
        output=OutputSpec(
            pdb_path=_get(out, "pdb_path", "step5_input.pdb"),
            gro_path=_get(out, "gro_path", "step5_input.gro"),
            report_path=_get(out, "report_path", "replacement_report"),
            inspection_pdb_path=_get(out, "inspection_pdb_path", "aligned_replacement_inspection.pdb"),
        ),
        minimization=MinimizationSpec(
            enabled=_get(mn, "enabled", False),
            output_dir=_get(mn, "output_dir", "openmm_minimization"),
            resources=_get(mn, "resources", {}),
            restraint_force_constant_kj_mol_nm2=_get(mn, "restraint_force_constant_kj_mol_nm2", 1000.0),
            restraint_residue_classes=_get(mn, "restraint_residue_classes", {}),
            coordinates_path=_get(mn, "coordinates_path", ""),
            topology_path=_get(mn, "topology_path", ""),
            include_dir=_get(mn, "include_dir", ""),
            output_gro_path=_get(mn, "output_gro_path", "minimized.gro"),
            report_path=_get(mn, "report_path", "minimization_report"),
            algorithm=_get(mn, "algorithm", "lbfgs"),
            tolerance_kj_mol_nm=_get(mn, "tolerance_kj_mol_nm", 10.0),
            max_iterations=_get(mn, "max_iterations", 5000),
            constraint_tolerance=_get(mn, "constraint_tolerance", 1.0e-5),
            nonbonded_method=_get(mn, "nonbonded_method", "PME"),
            nonbonded_cutoff_nm=_get(mn, "nonbonded_cutoff_nm", 1.2),
            switch_distance_nm=_get(mn, "switch_distance_nm", 1.0),
            constraints=_get(mn, "constraints", "h_bonds"),
            rigid_water=_get(mn, "rigid_water", True),
            ewald_error_tolerance=_get(mn, "ewald_error_tolerance", 5.0e-4),
            use_dispersion_correction=_get(mn, "use_dispersion_correction", False),
            platform=_get(mn, "platform", "auto"),
            precision=_get(mn, "precision", "mixed"),
            device_index=str(_get(mn, "device_index", "")),
            defines=_get(mn, "defines", {}) or {},
        ),
        name_restoration=NameRestorationSpec(
            enabled=_get(nr, "enabled", False),
            method=namefix_method,
            reference=_get(nr, "reference", "replacement_structure"),
            apply_to=_get(nr, "apply_to", "target_box"),
            min_jaccard=_get(nr, "min_jaccard", 0.95),
            reference_topol=reference_topol,
            pdb_to_full_resname=pdb_to_full_resname,
        ),
        topology=TopologySpec(
            enabled=_get(tp, "enabled", False),
            protein_toppar_dir=_get(tp, "protein_toppar_dir", ""),
            protein_template_top=_get(tp, "protein_template_top", ""),
            environment_toppar_dir=_get(tp, "environment_toppar_dir", ""),
            ligand_itp_paths=_get(tp, "ligand_itp_paths", []) or [],
            moleculetype_overrides=_get(tp, "moleculetype_overrides", {}) or {},
            output_dir=_get(tp, "output_dir", "."),
        ),
        ndx=NdxSpec(
            enabled=_get(ix, "enabled", True),
            output_path=_get(ix, "output_path", "index.ndx"),
            style=_get(ix, "style", "both"),
            extra_groups=_get(ix, "extra_groups", {}) or {},
        ),
        refgro=raw.get("refgro", None),
        ligand_replace=LigandReplaceSpec(
            enabled=(mode == "lig"),
            structure_path=_get(lr, "structure_path", ""),
            format=_get(lr, "format", "auto"),
            box_dimensions=_get(lr, "box_dimensions", None),
            original_ligand=OriginalLigandSpec(
                resname=_get(lr_ol, "resname", ""),
            ),
            new_ligand=NewLigandSpec(
                coord_path=_get(lr_nl, "coord_path", ""),
                format=_get(lr_nl, "format", "auto"),
                resname=_get(lr_nl, "resname", ""),
                itp_path=_get(lr_nl, "itp_path", ""),
                forcefield_path=_get(lr_nl, "forcefield_path", "") or "",
            ),
            fit=LigandFitSpec(
                method=_get(lr_fit, "method", "autofit"),
                old_ligand_fit_atoms=_get(lr_fit, "old_ligand_fit_atoms", []) or [],
                new_ligand_fit_atoms=_get(lr_fit, "new_ligand_fit_atoms", []) or [],
            ),
        ),
        cholesterol=CholesterolSpec(
            experimental_structure_path=_get(cs, "experimental_structure_path", ""),
            target_system_path=_get(cs, "target_system_path", ""),
            target_system_format=_get(cs, "target_system_format", "auto"),
            box_dimensions=_get(cs, "box_dimensions", None),
            target_protein_mask=_get(cs, "target_protein_mask", ""),
            experimental_protein_mask=_get(cs, "experimental_protein_mask", ""),
            cholesterol_resnames=[
                str(x).upper() for x in (
                    _get(cs, "cholesterol_resnames", ["CHL1", "CHL", "CHOL", "CLR"])
                    or []
                )
            ],
            charmm_resname=str(_get(cs, "charmm_resname", "CHL1")).upper(),
            convert_to_charmm36=_get(cs, "convert_to_charmm36", True),
            reconstruct_missing_heavy_atoms=_get(
                cs, "reconstruct_missing_heavy_atoms", True
            ),
            max_missing_heavy_atoms=_get(cs, "max_missing_heavy_atoms", 2),
            max_heavy_atom_fit_rmsd=_get(cs, "max_heavy_atom_fit_rmsd", 1.5),
            write_diagnostics=_get(cs, "write_diagnostics", False),
            charmm36_reference_path=_get(cs, "charmm36_reference_path", ""),
            obabel_command=_get(cs, "obabel_command", "obabel"),
            verify_tolerance=_get(cs, "verify_tolerance", 0.40),
            raw_protonated_pdb_path=_get(
                cs, "raw_protonated_pdb_path", "cholesterol_protonated_raw.pdb"
            ),
            converted_pdb_path=_get(
                cs, "converted_pdb_path", "cholesterol_charmm36.pdb"
            ),
            aligned_reference_pdb_path=_get(
                cs, "aligned_reference_pdb_path", "cholesterol_reference_aligned.pdb"
            ),
            reference_placed_pdb_path=_get(
                cs, "reference_placed_pdb_path", "cholesterol_reference_placed.pdb"
            ),
            merged_pdb_path=_get(cs, "merged_pdb_path", "merged_system.pdb"),
            composition=CholesterolCompositionSpec(
                reference_system_path=_get(cs_comp, "reference_system_path", ""),
                distance_from_protein=_get(cs_comp, "distance_from_protein", 20.0),
                target_concentration=_get(cs_comp, "target_concentration", 0.15),
                concentration_tolerance_fraction=_get(
                    cs_comp, "concentration_tolerance_fraction", 0.05
                ),
                salt_cation=str(_get(cs_comp, "salt_cation", "SOD")).upper(),
                salt_anion=str(_get(cs_comp, "salt_anion", "CLA")).upper(),
                random_seed=_get(cs_comp, "random_seed", 42),
                lipid_targets=lipid_targets,
            ),
            repack=RepackSpec(
                enabled=cs_repack.get("enabled", repack_defaults.enabled),
                output_dir=cs_repack.get("output_dir", repack_defaults.output_dir),
                moleculetype=str(cs_repack.get("moleculetype", repack_defaults.moleculetype)),
                lipid_restraint=cs_repack.get("lipid_restraint", repack_defaults.lipid_restraint),
                temperature=cs_repack.get("temperature", repack_defaults.temperature),
                stages=repack_stages,
            ),
        ),
    )
    anchor_output_defaults(cfg, output_root)
    cfg.validate(check_paths=check_paths)
    return cfg


def load_minimization_config(
    path: str, check_paths: bool = True, output_root=None, bundle: bool = False
) -> MinimizationSpec:
    """Load a standalone minimization-only YAML file.

    `bundle` is set only by the portable runner inside a prepared minimization
    folder, where the configuration and its inputs travel together: there the
    folder holding the configuration is the explicit root, so the bundle stays
    relocatable.
    """
    try:
        with open(path, encoding="utf-8") as fh:
            raw = yaml.load(fh, Loader=_UniqueKeySafeLoader)
    except ConfigError:
        raise
    except yaml.YAMLError as exc:
        raise ConfigError(f"Cannot parse YAML config '{path}': {exc}") from exc
    except OSError as exc:
        raise ConfigError(f"Cannot read config file '{path}': {exc}") from exc
    if not raw:
        raise ConfigError(f"Config file '{path}' is empty")
    if not isinstance(raw, dict):
        raise ConfigError(f"Config file '{path}' must contain a top-level mapping")
    unexpected = sorted(set(raw) - {"paths", "minimization"})
    if unexpected:
        raise ConfigError(
            "Standalone minimization accepts only paths and minimization; found: "
            + ", ".join(unexpected)
        )
    raw = _resolve_path_references(raw)
    raw = anchor_config_paths(raw, config_path=path, output_root=output_root,
                              input_root=Path(path).resolve().parent if bundle else None)
    paths = _require_mapping(raw.get("paths"), "paths")
    for key, value in paths.items():
        if not isinstance(key, str) or not re.fullmatch(r"[a-z][a-z0-9_]*", key):
            raise ConfigError("paths keys must use lower-case snake_case names")
        if value is not None and not isinstance(value, str):
            raise ConfigError(f"paths.{key} must be a path string or null")
    values = _require_mapping(raw.get("minimization"), "minimization")
    _reject_unknown_keys(values, "minimization", _SECTION_KEYS["minimization"])
    if "defines" in values and not isinstance(values["defines"], dict):
        raise ConfigError("minimization.defines must be a mapping")
    try:
        spec = MinimizationSpec(**values)
    except TypeError as exc:
        raise ConfigError(f"Invalid minimization configuration: {exc}") from exc
    if not Path(spec.output_dir).expanduser().is_absolute():
        spec.output_dir = str(resolve_output_root(output_root) / spec.output_dir)
    spec.validate(check_paths=check_paths, require_inputs=True)
    return spec
