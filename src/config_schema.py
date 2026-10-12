"""The dataclasses of the pipeline YAML sections (protein, ligand, cholesterol modes)."""

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


@dataclass
class StructureSpec:
    path: str = ""
    format: str = "auto"          # auto | gro | pdb
    protein_mask: str = ""       # cpptraj-style mask, sequential numbering
    # Explicit [a, b, c, alpha, beta, gamma] in Å/degrees when the input
    # has no usable box. A different file's box is never inferred implicitly.
    # Leave null to use the frame's own dimensions.
    box_dimensions: Optional[List[float]] = None

    def resolved_format(self) -> str:
        if self.format != "auto":
            return self.format
        ext = os.path.splitext(self.path)[1].lower().lstrip(".")
        if ext not in ("gro", "pdb"):
            raise ConfigError(f"Cannot auto-detect format for '{self.path}' (extension '.{ext}')")
        return ext


@dataclass
class AlignmentSpec:
    # pymol_align   -- PyMOL cmd.align(): sequence alignment (Needleman-
    #                  Wunsch) + iterative outlier-rejection superposition.
    #                  Tolerant of the two sides not being 1:1 the same
    #                  length/order -- PyMOL works that out itself.
    # pymol_cealign -- PyMOL cmd.cealign(): structure-based (no sequence
    #                  needed), also its own internal correspondence-finding.
    # mask_fit      -- cpptraj-'rms'-style: region_mask_original/
    #                  region_mask_replacement (or the full protein_mask on
    #                  either side if left unset) are matched POSITIONALLY,
    #                  one atom to the next in mask/file order -- no
    #                  sequence alignment, no outlier rejection. Plain
    #                  Kabsch/SVD best fit (see kabsch.py), same as
    #                  cpptraj's `rms <ref_mask> <mask>` given two explicit
    #                  masks. Requires both fit selections to resolve to
    #                  exactly the same atom count (fail-stop if not) --
    #                  it's on you to make sure the two masks correspond
    #                  residue-for-residue, same as in cpptraj. cycles/
    #                  cutoff below are ignored for this method (no PyMOL,
    #                  no outlier-rejection cycles).
    method: str = "pymol_align"          # pymol_align | pymol_cealign | mask_fit
    # Deprecated compatibility fallback. New configurations put @CA or
    # @C,CA,N,O directly in each region mask.
    atoms: Optional[str] = None           # legacy: CA | backbone | all
    region_mask_original: Optional[str] = None
    region_mask_replacement: Optional[str] = None
    cycles: int = 5
    cutoff: float = 2.0


@dataclass
class ReplacementLigandSpec:
    """Charge and topology metadata for a ligand inside the replacement mask."""

    resname: str = ""
    itp_path: str = ""
    # Give a numeric net charge or use "from_itp" to sum the [ atoms ] charges.
    charge: Union[float, str] = "from_itp"
    # Optional when the ITP contains one moleculetype or its name matches resname.
    moleculetype: str = ""
    # Force field with this ligand's bonded parameters (the forcefield.itp
    # CHARMM-GUI generated with the ITP). Empty: use a forcefield.itp next to
    # itp_path if present, else the environment force field.
    forcefield_path: str = ""


@dataclass
class ClashDetectionSpec:
    # One removal cutoff for every mode and every molecule class, matching the
    # CHARMM-GUI transplant mode. No class -- cholesterol included -- gets its
    # own implicit cutoff; override per class only via `thresholds`.
    threshold: float = 0.80
    thresholds: Dict[str, float] = field(default_factory=dict)
    use_pbc: bool = True
    heavy_atoms_only: bool = True
    flag_lipid_removal_fraction: float = 0.05
    # Classes to check and REPORT but never remove (e.g. ["lipid"] to keep
    # the membrane fully intact regardless of protein overlap, while still
    # surfacing exactly which contacts exceeded the threshold -- see the
    # report's "flagged_but_kept" list, which names the specific contacting
    # atom on each side, not just the residue).
    keep_classes: List[str] = field(default_factory=list)

    def threshold_for(self, molclass: str) -> float:
        return self.thresholds.get(molclass, self.threshold)


@dataclass
class ChargeSpec:
    target_net_charge: float = 0.0
    neutralize: bool = True
    tolerance: float = 0.01
    exclusion_distance_from_protein: float = 10.0
    exclusion_distance_from_ligand: float = 10.0
    # Minimum PBC distance from any lipid heavy atom required for an ion to
    # be eligible for neutralization. This replaces the unsafe global lipid
    # min-z/max-z slab used by older versions.
    exclusion_distance_from_lipid: float = 5.0
    # Deprecated input alias for exclusion_distance_from_lipid. It remains in
    # the dataclass so old YAML files fail safely and predictably.
    membrane_z_margin: Optional[float] = None
    random_seed: int = 42
    charge_table_overrides: Dict[str, Optional[float]] = field(default_factory=dict)
    # When enabled, neutralization candidates must also be at least
    # exclusion_distance_from_lipid from every lipid heavy atom under PBC.
    # The historical name is retained for YAML compatibility; no z slab is
    # constructed.
    exclude_membrane_interior: bool = True
    # Deprecated because a fixed absolute-z slab is unsafe for wrapped,
    # unwrapped, tilted, or displaced membranes. A non-null value is rejected.
    membrane_z_override: Optional[List[float]] = None


@dataclass
class NameRestorationSpec:
    # Fixes nonprotein resnames truncated by a strict PDB field, such as
    # TIP3->TIP, CHL1->CHL, and both POPC/POPE->POP.
    enabled: bool = False
    # atom_signature uses an untruncated coordinate reference. itp_atom_count
    # validates names against reference_topol and resolves shared aliases by
    # comparing target-residue atom counts with environment ITP [ atoms ] rows.
    method: str = "atom_signature"       # atom_signature | itp_atom_count
    reference: str = "replacement_structure"   # "target_box" | "replacement_structure"
    apply_to: str = "target_box"             # "target_box" | "replacement_structure" | "both"
    min_jaccard: float = 0.95
    reference_topol: str = ""
    # Observed target-PDB alias -> one or more permitted full residue names.
    # Example: {POP: [POPC, POPE], TIP: [TIP3], CHL: [CHL1]}.
    pdb_to_full_resname: Dict[str, List[str]] = field(default_factory=dict)


@dataclass
class OriginalLigandSpec:
    # RES selects the only residue with that name. RES:N selects the Nth
    # occurrence in file order when multiple residues share the name.
    resname: str = ""


@dataclass
class NewLigandSpec:
    coord_path: str = ""     # coordinate file for the incoming ligand (any MDAnalysis-readable format)
    format: str = "auto"
    # Resname to assign to every atom of the loaded ligand in the final
    # system. Leave "" to keep whatever resname is already in coord_path.
    resname: str = ""
    # Standalone .itp for this ligand -- automatically added to
    # topology.ligand_itp_paths at run time if not already present there
    # (only matters when topology.enabled is true), so you don't have to
    # list it in two places.
    itp_path: str = ""
    # Force field with this ligand's bonded parameters (the forcefield.itp
    # CHARMM-GUI generated with the ITP). Empty: use a forcefield.itp next to
    # itp_path if present, else the environment force field.
    forcefield_path: str = ""


@dataclass
class LigandFitSpec:
    method: str = "autofit"   # "pairfit" | "autofit" | "mcsfit" | "nofit"
    # Atom NAMES used for the fit, matched positionally: old_ligand_fit_atoms[i]
    # in the structure's own original ligand corresponds to
    # new_ligand_fit_atoms[i] in the incoming ligand. Only used when
    # method: pairfit. Ignored (may be left empty) for autofit and nofit.
    old_ligand_fit_atoms: List[str] = field(default_factory=list)
    new_ligand_fit_atoms: List[str] = field(default_factory=list)


@dataclass
class LigandReplaceSpec:
    # In lig mode, the protein, membrane, solvent, and ions stay in place.
    # Only one ligand is swapped inside the supplied complete structure.
    enabled: bool = False
    structure_path: str = ""   # the ONE file: protein+membrane+solvent+ions+old ligand, all kept except the old ligand
    format: str = "auto"
    # Same purpose as StructureSpec.box_dimensions above -- only needed if
    # structure_path's own box info can't be read.
    box_dimensions: Optional[List[float]] = None
    original_ligand: OriginalLigandSpec = field(default_factory=OriginalLigandSpec)
    new_ligand: NewLigandSpec = field(default_factory=NewLigandSpec)
    fit: LigandFitSpec = field(default_factory=LigandFitSpec)


@dataclass
class CholesterolCompositionSpec:
    """Final per-leaflet lipid targets and reference salt baseline."""

    reference_system_path: str = ""
    distance_from_protein: float = 20.0
    # The reference system is assumed to have this concentration.  The final
    # concentration is inferred from its retained salt-pair fraction rather
    # than recalculated from the box volume.
    target_concentration: float = 0.15
    concentration_tolerance_fraction: float = 0.05
    salt_cation: str = "SOD"
    salt_anion: str = "CLA"
    random_seed: int = 42
    # Shape: {"CHL1": {"upper": 40, "lower": 40}, ...}
    lipid_targets: Dict[str, Dict[str, int]] = field(default_factory=dict)


@dataclass
class RepackStageSpec:
    """One stage of the repacking equilibration."""

    name: str = ""
    integrator: str = "npt"  # steep (minimization) | nvt | npt
    nsteps: int = 0
    dt: float = 0.002        # ps; not used by steep
    # kJ/mol/nm^2 on every heavy atom of the protein, ligands and restored cholesterols.
    restraint: float = 1000.0


def _default_repack_stages():
    # Protein, ligands and restored cholesterols: 4000 -> 2000 -> 1000, then held
    # at 1000 to the end; about 1.4 ns of MD, the last stage 1 ns.
    return [
        RepackStageSpec("step6.0_minimization", "steep", 5000, 0.0, 4000.0),
        RepackStageSpec("step6.1_equilibration", "nvt", 125000, 0.001, 4000.0),
        RepackStageSpec("step6.2_equilibration", "nvt", 125000, 0.001, 2000.0),
        RepackStageSpec("step6.3_equilibration", "npt", 125000, 0.001, 1000.0),
        RepackStageSpec("step6.4_equilibration", "npt", 500000, 0.002, 1000.0),
    ]


@dataclass
class RepackSpec:
    """GROMACS repacking equilibration written next to the restored system (repack.py)."""

    enabled: bool = True
    output_dir: str = "repack"          # folder next to the final GRO
    moleculetype: str = "CHLR"          # restrained copy of the cholesterol moleculetype
    lipid_restraint: float = 0.0        # POSRES_FC_LIPID and DIHRES_FC for every other lipid
    temperature: float = 310.0          # K
    stages: List[RepackStageSpec] = field(default_factory=_default_repack_stages)


@dataclass
class CholesterolSpec:
    """Inputs and conversion settings for experimental cholesterol restoration."""

    experimental_structure_path: str = ""
    target_system_path: str = ""
    target_system_format: str = "auto"
    box_dimensions: Optional[List[float]] = None
    target_protein_mask: str = ""
    experimental_protein_mask: str = ""
    cholesterol_resnames: List[str] = field(
        default_factory=lambda: ["CHL1", "CHL", "CHOL", "CLR"]
    )
    charmm_resname: str = "CHL1"
    convert_to_charmm36: bool = True
    # Missing heavy atoms may be placed from the fitted CHARMM36 template only
    # when the surviving reference names make the correspondence unambiguous.
    reconstruct_missing_heavy_atoms: bool = True
    max_missing_heavy_atoms: int = 2
    max_heavy_atom_fit_rmsd: float = 1.5
    # Conversion intermediates and visual inspection PDBs are opt-in. The
    # final PDB/GRO, reports, topology, and index are unaffected by this flag.
    write_diagnostics: bool = False
    charmm36_reference_path: str = ""
    obabel_command: str = "obabel"
    verify_tolerance: float = 0.40
    raw_protonated_pdb_path: str = "cholesterol_protonated_raw.pdb"
    converted_pdb_path: str = "cholesterol_charmm36.pdb"
    aligned_reference_pdb_path: str = "cholesterol_reference_aligned.pdb"
    reference_placed_pdb_path: str = "cholesterol_reference_placed.pdb"
    merged_pdb_path: str = "merged_system.pdb"
    composition: CholesterolCompositionSpec = field(
        default_factory=CholesterolCompositionSpec
    )
    repack: RepackSpec = field(default_factory=RepackSpec)


@dataclass
class OutputSpec:
    # Flat by default -- everything lands in the current working directory
    # (the directory from which run_pipeline.py is launched),
    # same as the old merge_pdb.sh: no separate output/ or reports/
    # subfolders unless you put a "/" in one of these paths yourself.
    pdb_path: str = "step5_input.pdb"
    gro_path: str = "step5_input.gro"
    report_path: str = "replacement_report"
    inspection_pdb_path: str = "aligned_replacement_inspection.pdb"


@dataclass
class MinimizationSpec:
    """Preparation settings plus physical settings for the exported OpenMM runner."""

    enabled: bool = False
    output_dir: str = "openmm_minimization"
    resources: dict = field(default_factory=dict)
    restraint_force_constant_kj_mol_nm2: float = 1000.0
    restraint_residue_classes: dict = field(default_factory=dict)
    # Standalone preparation inputs. For post-pipeline preparation these
    # are deliberately replaced with the freshly written output.gro_path and
    # topology.output_dir/topol.top, preventing an older system from being
    # exported by mistake.
    coordinates_path: str = ""
    topology_path: str = ""
    include_dir: str = ""
    output_gro_path: str = "minimized.gro"
    report_path: str = "minimization_report"

    # OpenMM's supported, constrained local minimizer is L-BFGS.
    algorithm: str = "lbfgs"
    tolerance_kj_mol_nm: float = 10.0
    max_iterations: int = 5000
    constraint_tolerance: float = 1.0e-5

    # Force calculation choices that would otherwise live in a GROMACS MDP.
    nonbonded_method: str = "PME"
    nonbonded_cutoff_nm: float = 1.2
    # null or 0 disables LJ switching.
    switch_distance_nm: Optional[float] = 1.0
    constraints: str = "h_bonds"  # none | h_bonds | all_bonds | h_angles
    rigid_water: bool = True
    ewald_error_tolerance: float = 5.0e-4
    use_dispersion_correction: bool = False

    platform: str = "auto"  # auto | CUDA | HIP | OpenCL | CPU | Reference
    precision: str = "mixed"  # single | mixed | double (when supported)
    device_index: str = ""
    # Topology preprocessor definitions; positional restraints must be inactive.
    defines: Dict[str, Union[str, int, float]] = field(default_factory=dict)

    def validate(self, check_paths: bool = True, require_inputs: bool = False):
        errors = []
        if not isinstance(self.output_dir, str) or not self.output_dir.strip():
            errors.append("minimization.output_dir must be a non-empty path")
        if not isinstance(self.resources, dict):
            errors.append("minimization.resources must be a mapping")
        if (not _is_finite_number(self.restraint_force_constant_kj_mol_nm2)
                or self.restraint_force_constant_kj_mol_nm2 <= 0):
            errors.append("minimization.restraint_force_constant_kj_mol_nm2 must be positive")
        if not isinstance(self.restraint_residue_classes, dict):
            errors.append("minimization.restraint_residue_classes must be a mapping")
        else:
            for name, role in self.restraint_residue_classes.items():
                if (not isinstance(name, str) or not name.strip()
                        or role not in ("protein", "ligand", "lipid", "water", "ion")):
                    errors.append("restraint_residue_classes must map residue names to protein, ligand, lipid, water or ion")
        if type(self.enabled) is not bool:
            errors.append("minimization.enabled must be boolean")
        if self.algorithm != "lbfgs":
            errors.append("minimization.algorithm must be lbfgs")
        if not _is_finite_number(self.tolerance_kj_mol_nm) or self.tolerance_kj_mol_nm <= 0:
            errors.append("minimization.tolerance_kj_mol_nm must be > 0")
        if type(self.max_iterations) is not int or self.max_iterations < 0:
            errors.append("minimization.max_iterations must be an integer >= 0")
        if not _is_finite_number(self.constraint_tolerance) or not 0 < self.constraint_tolerance < 1:
            errors.append("minimization.constraint_tolerance must be between 0 and 1")
        methods = {"NoCutoff", "CutoffNonPeriodic", "CutoffPeriodic", "Ewald", "PME", "LJPME"}
        if self.nonbonded_method not in methods:
            errors.append(
                "minimization.nonbonded_method must be one of: "
                + ", ".join(sorted(methods))
            )
        if not _is_finite_number(self.nonbonded_cutoff_nm) or self.nonbonded_cutoff_nm <= 0:
            errors.append("minimization.nonbonded_cutoff_nm must be > 0")
        switch = self.switch_distance_nm
        if switch is not None:
            if not _is_finite_number(switch) or switch < 0:
                errors.append("minimization.switch_distance_nm must be null or >= 0")
            elif switch > 0 and _is_finite_number(self.nonbonded_cutoff_nm) and switch >= self.nonbonded_cutoff_nm:
                errors.append(
                    "minimization.switch_distance_nm must be smaller than "
                    "minimization.nonbonded_cutoff_nm"
                )
        if self.constraints not in {"none", "h_bonds", "all_bonds", "h_angles"}:
            errors.append(
                "minimization.constraints must be one of: none, h_bonds, "
                "all_bonds, h_angles"
            )
        for name in ("rigid_water", "use_dispersion_correction"):
            if type(getattr(self, name)) is not bool:
                errors.append(f"minimization.{name} must be boolean")
        if (
            not _is_finite_number(self.ewald_error_tolerance)
            or not 0 < self.ewald_error_tolerance < 1
        ):
            errors.append("minimization.ewald_error_tolerance must be between 0 and 1")
        if self.platform not in {"auto", "CUDA", "HIP", "OpenCL", "CPU", "Reference"}:
            errors.append(
                "minimization.platform must be one of: auto, CUDA, HIP, "
                "OpenCL, CPU, Reference"
            )
        if self.precision not in {"single", "mixed", "double"}:
            errors.append("minimization.precision must be single, mixed, or double")
        if not isinstance(self.device_index, str):
            errors.append("minimization.device_index must be a string")
        if not isinstance(self.defines, dict):
            errors.append("minimization.defines must be a mapping")
        else:
            for key, value in self.defines.items():
                if not isinstance(key, str) or not key.strip():
                    errors.append("minimization.defines keys must be non-empty strings")
                if isinstance(value, (dict, list, tuple, set)) or value is None:
                    errors.append(
                        f"minimization.defines.{key} must be a scalar string or number"
                    )
        for name in ("output_gro_path", "report_path"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip():
                errors.append(f"minimization.{name} must be a non-empty path string")
        for name in ("coordinates_path", "topology_path", "include_dir"):
            if not isinstance(getattr(self, name), str):
                errors.append(f"minimization.{name} must be a path string")
        if require_inputs:
            for name in ("coordinates_path", "topology_path"):
                value = getattr(self, name)
                if not isinstance(value, str) or not value.strip():
                    errors.append(f"minimization.{name} is required")
            if (
                isinstance(self.coordinates_path, str)
                and self.coordinates_path
                and os.path.splitext(self.coordinates_path)[1].lower() != ".gro"
            ):
                errors.append("minimization.coordinates_path must be a .gro file")
            if check_paths and self.coordinates_path and not os.path.isfile(self.coordinates_path):
                errors.append(
                    f"minimization.coordinates_path does not exist: {self.coordinates_path}"
                )
            if check_paths and self.topology_path and not os.path.isfile(self.topology_path):
                errors.append(
                    f"minimization.topology_path does not exist: {self.topology_path}"
                )
        if self.include_dir and check_paths and not os.path.isdir(self.include_dir):
            errors.append(f"minimization.include_dir does not exist: {self.include_dir}")

        paths = {
            "coordinates_path": self.coordinates_path,
            "topology_path": self.topology_path,
            "output_gro_path": self.output_gro_path,
            "report_path (.txt)": self.report_path + ".txt",
            "report_path (.json)": self.report_path + ".json",
        }
        by_realpath = {}
        for label, value in paths.items():
            if isinstance(value, str) and value:
                by_realpath.setdefault(os.path.realpath(value), []).append(label)
        for labels in by_realpath.values():
            if len(labels) > 1:
                errors.append(
                    "Minimization input/output paths collide: " + ", ".join(labels)
                )
        if errors:
            raise ConfigError("Invalid configuration:\n  - " + "\n  - ".join(errors))


@dataclass
class TopologySpec:
    # Builds a full GROMACS input set (topol.top + a self-contained
    # toppar/ directory) describing the final system. Generation is
    # disabled by default. See topology.py's
    # module docstring for the full explanation of how moleculetype names
    # get resolved; this is a system-agnostic mechanism, not tied to any
    # particular protein/ligand -- every name below is a path you supply.
    enabled: bool = False
    # toppar/ folder belonging to whichever structure the NEW protein (and
    # usually its ligand) came from -- e.g. replacement_structure's own
    # ".../gromacs/toppar/". Must sit next to that structure's own
    # topol.top (standard CHARMM-GUI layout), unless you override that via
    # protein_template_top below.
    protein_toppar_dir: str = ""
    # Optional override if protein_toppar_dir's sibling topol.top isn't
    # at the standard "<parent of toppar_dir>/topol.top" location.
    protein_template_top: str = ""
    # toppar/ folder with correct parameters for the kept environment's
    # lipid/water/ion species -- can be from any CHARMM-GUI system with a
    # matching composition (its own topol.top is not read, only its .itp
    # files).
    environment_toppar_dir: str = ""
    # Extra standalone ligand .itp file(s), only needed if the ligand's
    # parameters aren't already inside protein_toppar_dir.
    ligand_itp_paths: List[str] = field(default_factory=list)
    # resname -> moleculetype name overrides for anything that doesn't
    # follow the "moleculetype name == resname" convention (water's
    # TIP3->"SOL" fallback is automatic; add others here as needed).
    moleculetype_overrides: Dict[str, str] = field(default_factory=dict)
    # "." = current working directory -- topol.top and toppar/ land right
    # next to everything else, exactly where merge_pdb.sh used to put them.
    # Set to a subfolder name (e.g. "gromacs_input") if you'd rather keep
    # them separate.
    output_dir: str = "."


@dataclass
class NdxSpec:
    # index.ndx is always written with the system. `enabled: false` is still
    # accepted so older configurations keep working; the GUI never offers it.
    enabled: bool = True
    output_path: str = "index.ndx"
    # generic | charmm_gui | both
    style: str = "both"
    # Extra/custom groups: {group_name: cpptraj-style mask}, resolved
    # against the FINAL assembled system. Use this to add your own group
    # names/aliases alongside the built-in generic set (Protein / Ligand /
    # Lipid / Water / Ion / Protein_Ligand / Solvent / Membrane / System).
    extra_groups: Dict[str, str] = field(default_factory=dict)


@dataclass
class BoxValidationSpec:
    mode: str = "strict"
    hard_clash_cutoff: float = 1.2  # Å; nonbonded pairs only
    molecule_protrusion_allowance: float = 5.0  # Å; diagnostic, not a rejection rule
    search_sample_size: int = 20000
    random_seed: int = 42
    wrap_after_validation: bool = False
    search_starts: int = 12
    search_max_iterations: int = 150
    max_candidate_verifications: int = 12
    max_search_seconds: float = 45.0

    def validate(self):
        errors = []
        if self.mode not in ("strict", "repair", "off"):
            errors.append("box_validation.mode must be strict, repair, or off")
        for name, minimum in (("hard_clash_cutoff", 0), ("max_search_seconds", 0),
                              ("molecule_protrusion_allowance", -1)):
            value = getattr(self, name)
            if (isinstance(value, bool) or not isinstance(value, (int, float))
                    or not math.isfinite(value) or value <= minimum
                    or (name == "molecule_protrusion_allowance" and value < 0)):
                errors.append(f"box_validation.{name} must be finite and {'>= 0' if minimum == -1 else '> 0'}")
        for name, minimum in (("search_sample_size", 8), ("random_seed", 0),
                              ("search_starts", 1), ("search_max_iterations", 1),
                              ("max_candidate_verifications", 1)):
            value = getattr(self, name)
            if type(value) is not int or value < minimum:
                errors.append(f"box_validation.{name} must be an integer >= {minimum}")
        if type(self.wrap_after_validation) is not bool:
            errors.append("box_validation.wrap_after_validation must be boolean")
        if errors:
            raise ConfigError("; ".join(errors))


@dataclass
class Config:
    mode: str = "protein"
    box_validation: BoxValidationSpec = field(default_factory=BoxValidationSpec)
    target_box: StructureSpec = field(default_factory=StructureSpec)
    replacement_structure: StructureSpec = field(default_factory=StructureSpec)
    alignment: AlignmentSpec = field(default_factory=AlignmentSpec)
    replacement_ligands: List[ReplacementLigandSpec] = field(default_factory=list)
    clash_detection: ClashDetectionSpec = field(default_factory=ClashDetectionSpec)
    charge: ChargeSpec = field(default_factory=ChargeSpec)
    output: OutputSpec = field(default_factory=OutputSpec)
    minimization: MinimizationSpec = field(default_factory=MinimizationSpec)
    name_restoration: NameRestorationSpec = field(default_factory=NameRestorationSpec)
    topology: TopologySpec = field(default_factory=TopologySpec)
    ndx: NdxSpec = field(default_factory=NdxSpec)
    # Ligand-only settings used when mode is "lig".
    ligand_replace: LigandReplaceSpec = field(default_factory=LigandReplaceSpec)
    cholesterol: CholesterolSpec = field(default_factory=CholesterolSpec)
    # A reference GRO used two ways: (1) sanity-checks the final GRO
    # writer's >99999-atom index wraparound (confirmed: GROMACS/MDAnalysis
    # wrap 99999 -> 0 -> 1 ..., not restarting at 1) against a real
    # large-atom-count example; (2) every non-protein resname it contains
    # is expected to also appear somewhere in the final output -- if one
    # is entirely missing, the report notes it (not a hard failure --
    # refgro is a shape reference, not necessarily this run's own donor
    # file). Optional; leave unset to skip both checks.
    refgro: Optional[str] = None

    def validate(self, check_paths: bool = True):
        errors: List[str] = []
        self.box_validation.validate()
        self.minimization.validate(check_paths=check_paths, require_inputs=False)
        structure_boxes = ({"target_box": self.target_box,
                            "replacement_structure": self.replacement_structure}
                           if self.mode == "protein" else
                           {"ligand_replace": self.ligand_replace}
                           if self.mode == "lig" else {"cholesterol": self.cholesterol})
        for label, structure in structure_boxes.items():
            dimensions = structure.box_dimensions
            if dimensions is None:
                continue
            if (not isinstance(dimensions, (list, tuple)) or len(dimensions) != 6
                    or any(isinstance(x, bool) or not isinstance(x, (int, float)) for x in dimensions)):
                raise ConfigError(f"{label}.box_dimensions must be six numeric values [a,b,c,alpha,beta,gamma] in Å/degrees")
            from box_validation import validate_box, BoxValidationError
            try:
                validate_box(dimensions)
            except BoxValidationError as exc:
                raise ConfigError(f"{label}.box_dimensions: {exc}") from exc
        if type(self.clash_detection.use_pbc) is not bool:
            errors.append("clash_detection.use_pbc must be boolean")
        elif not self.clash_detection.use_pbc and self.box_validation.mode != "off":
            errors.append("clash_detection.use_pbc: false requires explicit box_validation.mode: off")
        if self.mode not in ("protein", "lig", "chl"):
            errors.append("mode must be one of: protein, lig, chl")
        if self.mode == "lig":
            lr = self.ligand_replace
            if lr.format not in ("auto", "pdb", "gro"):
                errors.append("ligand_replace.format must be one of: auto, pdb, gro")
            if not isinstance(lr.structure_path, str) or not lr.structure_path.strip():
                errors.append("ligand_replace.structure_path is required in lig mode")
            elif check_paths and not os.path.exists(lr.structure_path):
                errors.append(f"ligand_replace.structure_path does not exist: {lr.structure_path}")
            if lr.box_dimensions is not None and len(lr.box_dimensions) != 6:
                errors.append("ligand_replace.box_dimensions must have exactly 6 values [x, y, z, alpha, beta, gamma]")
            if not isinstance(lr.original_ligand.resname, str) or not lr.original_ligand.resname.strip():
                errors.append("ligand_replace.original_ligand.resname is required (the ligand resname to remove)")
            elif not re.fullmatch(
                r"[A-Za-z0-9_+\-]+(?::[1-9]\d*)?", lr.original_ligand.resname.strip()
            ):
                errors.append(
                    "ligand_replace.original_ligand.resname must use RES or RES:N syntax"
                )
            if not isinstance(lr.new_ligand.coord_path, str) or not lr.new_ligand.coord_path.strip():
                errors.append("ligand_replace.new_ligand.coord_path is required")
            elif check_paths and not os.path.exists(lr.new_ligand.coord_path):
                errors.append(f"ligand_replace.new_ligand.coord_path does not exist: {lr.new_ligand.coord_path}")
            if lr.new_ligand.forcefield_path and (
                    not isinstance(lr.new_ligand.forcefield_path, str)
                    or (check_paths and not os.path.isfile(lr.new_ligand.forcefield_path))):
                errors.append("ligand_replace.new_ligand.forcefield_path does not exist or is not a file: "
                              f"{lr.new_ligand.forcefield_path}")
            if not isinstance(lr.new_ligand.format, str) or not lr.new_ligand.format.strip():
                errors.append("ligand_replace.new_ligand.format must be a non-empty format name")
            if lr.new_ligand.itp_path and (
                not isinstance(lr.new_ligand.itp_path, str)
                or (check_paths and not os.path.isfile(lr.new_ligand.itp_path))
            ):
                errors.append(
                    "ligand_replace.new_ligand.itp_path does not exist or is not a file: "
                    f"{lr.new_ligand.itp_path}"
                )
            if lr.fit.method not in ("pairfit", "autofit", "mcsfit", "nofit"):
                errors.append(
                    "ligand_replace.fit.method must be one of: pairfit, autofit, mcsfit, nofit"
                )
            elif lr.fit.method == "pairfit":
                if not lr.fit.old_ligand_fit_atoms or not lr.fit.new_ligand_fit_atoms:
                    errors.append(
                        "ligand_replace.fit.method is 'pairfit' but old_ligand_fit_atoms/"
                        "new_ligand_fit_atoms is empty -- list the atom names to fit on, or "
                        "switch to fit.method: autofit or nofit"
                    )
                elif len(lr.fit.old_ligand_fit_atoms) != len(lr.fit.new_ligand_fit_atoms):
                    errors.append(
                        f"ligand_replace.fit.old_ligand_fit_atoms ({len(lr.fit.old_ligand_fit_atoms)} atoms) "
                        f"and new_ligand_fit_atoms ({len(lr.fit.new_ligand_fit_atoms)} atoms) must be the same "
                        f"length -- they're matched positionally, one pair per fit atom"
                    )
                elif len(lr.fit.old_ligand_fit_atoms) < 3:
                    errors.append(
                        "ligand_replace.fit.pairfit requires at least three atom pairs"
                    )
                if len(lr.fit.old_ligand_fit_atoms) != len(set(lr.fit.old_ligand_fit_atoms)):
                    errors.append("ligand_replace.fit.old_ligand_fit_atoms contains duplicate names")
                if len(lr.fit.new_ligand_fit_atoms) != len(set(lr.fit.new_ligand_fit_atoms)):
                    errors.append("ligand_replace.fit.new_ligand_fit_atoms contains duplicate names")
                if any(not isinstance(name, str) or not name for name in lr.fit.old_ligand_fit_atoms):
                    errors.append("ligand_replace.fit.old_ligand_fit_atoms must contain non-empty strings")
                if any(not isinstance(name, str) or not name for name in lr.fit.new_ligand_fit_atoms):
                    errors.append("ligand_replace.fit.new_ligand_fit_atoms must contain non-empty strings")
        elif self.mode == "protein":
            for label, spec in (("target_box", self.target_box),
                                 ("replacement_structure", self.replacement_structure)):
                if not isinstance(spec.path, str) or not spec.path.strip():
                    errors.append(f"{label}.path is required")
                elif check_paths and not os.path.exists(spec.path):
                    errors.append(f"{label}.path does not exist: {spec.path}")
                if spec.format not in ("auto", "pdb", "gro"):
                    errors.append(f"{label}.format must be one of: auto, pdb, gro")
                if not spec.protein_mask:
                    errors.append(f"{label}.protein_mask is required (cpptraj-style, e.g. ':1-963')")
                if spec.box_dimensions is not None and len(spec.box_dimensions) != 6:
                    errors.append(f"{label}.box_dimensions must have exactly 6 values [x, y, z, alpha, beta, gamma]")
            if self.alignment.atoms not in (None, "", "CA", "backbone", "all"):
                errors.append("legacy alignment.atoms must be one of: CA, backbone, all")
            if self.alignment.method not in ("pymol_align", "pymol_cealign", "mask_fit"):
                errors.append("alignment.method must be one of: pymol_align, pymol_cealign, mask_fit")
            seen_ligands = set()
            for i, ligand in enumerate(self.replacement_ligands, 1):
                prefix = f"replacement_ligands[{i}]"
                if not ligand.resname:
                    errors.append(f"{prefix}.resname is required")
                elif not re.fullmatch(
                    r"[A-Za-z0-9_+\-]+(?::[1-9]\d*)?", ligand.resname
                ):
                    errors.append(f"{prefix}.resname must use RES or RES:N syntax")
                elif ligand.resname in seen_ligands:
                    errors.append(f"replacement_ligands contains duplicate resname: {ligand.resname}")
                seen_ligands.add(ligand.resname)
                if not ligand.itp_path:
                    errors.append(f"{prefix}.itp_path is required")
                elif check_paths and not os.path.isfile(ligand.itp_path):
                    errors.append(f"{prefix}.itp_path does not exist: {ligand.itp_path}")
                if ligand.forcefield_path and (
                        not isinstance(ligand.forcefield_path, str)
                        or (check_paths and not os.path.isfile(ligand.forcefield_path))):
                    errors.append(f"{prefix}.forcefield_path does not exist: {ligand.forcefield_path}")
                if isinstance(ligand.charge, str):
                    if ligand.charge.lower() != "from_itp":
                        errors.append(f"{prefix}.charge must be numeric or 'from_itp'")
                elif not _is_finite_number(ligand.charge):
                    errors.append(f"{prefix}.charge must be numeric or 'from_itp'")
            if type(self.name_restoration.enabled) is not bool:
                errors.append("name_restoration.enabled must be boolean")
            if self.name_restoration.reference not in ("target_box", "replacement_structure"):
                errors.append("name_restoration.reference must be one of: target_box, replacement_structure")
            if self.name_restoration.apply_to not in ("target_box", "replacement_structure", "both"):
                errors.append("name_restoration.apply_to must be one of: target_box, replacement_structure, both")
            if self.name_restoration.method not in ("atom_signature", "itp_atom_count"):
                errors.append(
                    "name_restoration.method must be one of: atom_signature, itp_atom_count"
                )
            if (
                not _is_finite_number(self.name_restoration.min_jaccard)
                or not 0.0 <= self.name_restoration.min_jaccard <= 1.0
            ):
                errors.append("name_restoration.min_jaccard must be between 0 and 1")
            if self.name_restoration.enabled and self.name_restoration.method == "itp_atom_count":
                if self.name_restoration.apply_to != "target_box":
                    errors.append(
                        "name_restoration.method itp_atom_count requires apply_to: target_box"
                    )
                if (
                    not isinstance(self.name_restoration.reference_topol, str)
                    or not self.name_restoration.reference_topol.strip()
                ):
                    errors.append(
                        "name_restoration.reference_topol is required for method: itp_atom_count"
                    )
                elif check_paths and not os.path.isfile(self.name_restoration.reference_topol):
                    errors.append(
                        "name_restoration.reference_topol does not exist: "
                        f"{self.name_restoration.reference_topol}"
                    )
                if not self.name_restoration.pdb_to_full_resname:
                    errors.append(
                        "name_restoration.pdb_to_full_resname cannot be empty for "
                        "method: itp_atom_count"
                    )
                all_full_names = []
                for pdb_name, full_names in self.name_restoration.pdb_to_full_resname.items():
                    if not pdb_name or not full_names or any(not name for name in full_names):
                        errors.append(
                            "name_restoration.pdb_to_full_resname cannot contain empty names "
                            "or empty candidate lists"
                        )
                    if len(full_names) != len(set(full_names)):
                        errors.append(
                            f"name_restoration.pdb_to_full_resname.{pdb_name} contains "
                            "duplicate full names"
                        )
                    all_full_names.extend(full_names)
                repeated_full_names = sorted(
                    name for name in set(all_full_names) if all_full_names.count(name) > 1
                )
                if repeated_full_names:
                    errors.append(
                        "Each full residue name may occur under only one PDB alias; "
                        f"repeated: {repeated_full_names}"
                    )
                if not self.topology.environment_toppar_dir:
                    errors.append(
                        "topology.environment_toppar_dir is required for "
                        "name_restoration.method: itp_atom_count"
                    )
                elif check_paths and not os.path.isdir(self.topology.environment_toppar_dir):
                    errors.append(
                        "topology.environment_toppar_dir does not exist: "
                        f"{self.topology.environment_toppar_dir}"
                    )
        elif self.mode == "chl":
            ch = self.cholesterol
            if ch.target_system_format not in ("auto", "pdb", "gro"):
                errors.append(
                    "cholesterol.target_system_format must be one of: auto, pdb, gro"
                )
            for label, value in (
                ("cholesterol.experimental_structure_path", ch.experimental_structure_path),
                ("cholesterol.target_system_path", ch.target_system_path),
            ):
                if not isinstance(value, str) or not value.strip():
                    errors.append(f"{label} is required in chl mode")
                elif check_paths and not os.path.exists(value):
                    errors.append(f"{label} does not exist: {value}")
            if ch.box_dimensions is not None and len(ch.box_dimensions) != 6:
                errors.append(
                    "cholesterol.box_dimensions must have exactly 6 values "
                    "[x, y, z, alpha, beta, gamma]"
                )
            if not ch.target_protein_mask:
                errors.append("cholesterol.target_protein_mask is required in chl mode")
            if not ch.experimental_protein_mask:
                errors.append("cholesterol.experimental_protein_mask is required in chl mode")
            if not ch.cholesterol_resnames:
                errors.append("cholesterol.cholesterol_resnames cannot be empty")
            elif any(not isinstance(name, str) or not name for name in ch.cholesterol_resnames):
                errors.append("cholesterol.cholesterol_resnames must contain non-empty strings")
            if not ch.charmm_resname:
                errors.append("cholesterol.charmm_resname cannot be empty")
            if type(ch.reconstruct_missing_heavy_atoms) is not bool:
                errors.append(
                    "cholesterol.reconstruct_missing_heavy_atoms must be boolean"
                )
            if type(ch.write_diagnostics) is not bool:
                errors.append("cholesterol.write_diagnostics must be boolean")
            if type(ch.convert_to_charmm36) is not bool:
                errors.append("cholesterol.convert_to_charmm36 must be boolean")
            if type(ch.max_missing_heavy_atoms) is not int or ch.max_missing_heavy_atoms < 0:
                errors.append(
                    "cholesterol.max_missing_heavy_atoms must be an integer >= 0"
                )
            if not _is_finite_number(ch.max_heavy_atom_fit_rmsd) or ch.max_heavy_atom_fit_rmsd <= 0:
                errors.append("cholesterol.max_heavy_atom_fit_rmsd must be > 0")
            if not _is_finite_number(ch.verify_tolerance) or ch.verify_tolerance < 0:
                errors.append("cholesterol.verify_tolerance must be >= 0")
            if ch.charmm36_reference_path and not isinstance(ch.charmm36_reference_path, str):
                errors.append("cholesterol.charmm36_reference_path must be a path string")
            elif check_paths and ch.charmm36_reference_path and not os.path.exists(ch.charmm36_reference_path):
                errors.append(
                    "cholesterol.charmm36_reference_path does not exist: "
                    f"{ch.charmm36_reference_path}"
                )
            if ch.convert_to_charmm36 and not ch.obabel_command:
                errors.append(
                    "cholesterol.obabel_command is required when convert_to_charmm36 is true"
                )
            comp = ch.composition
            if type(comp.random_seed) is not int or comp.random_seed < 0:
                errors.append(
                    "cholesterol.composition.random_seed must be an integer >= 0"
                )
            if check_paths and comp.reference_system_path and not os.path.exists(comp.reference_system_path):
                errors.append(
                    "cholesterol.composition.reference_system_path does not exist: "
                    f"{comp.reference_system_path}"
                )
            if not _is_finite_number(comp.distance_from_protein) or comp.distance_from_protein <= 0:
                errors.append("cholesterol.composition.distance_from_protein must be > 0")
            if not _is_finite_number(comp.target_concentration) or comp.target_concentration <= 0:
                errors.append("cholesterol.composition.target_concentration must be > 0")
            if (
                not _is_finite_number(comp.concentration_tolerance_fraction)
                or not 0.0 <= comp.concentration_tolerance_fraction <= 1.0
            ):
                errors.append(
                    "cholesterol.composition.concentration_tolerance_fraction "
                    "must be between 0 and 1"
                )
            if not comp.lipid_targets:
                errors.append("cholesterol.composition.lipid_targets cannot be empty")
            for resname, leaves in comp.lipid_targets.items():
                if not isinstance(leaves, dict):
                    errors.append(
                        f"cholesterol.composition.lipid_targets.{resname} must map upper/lower to counts"
                    )
                    continue
                missing = {"upper", "lower"} - {str(k).lower() for k in leaves}
                if missing:
                    errors.append(
                        f"cholesterol.composition.lipid_targets.{resname} is missing: "
                        + ", ".join(sorted(missing))
                    )
                for leaf, count in leaves.items():
                    if str(leaf).lower() not in ("upper", "lower"):
                        errors.append(
                            f"cholesterol.composition.lipid_targets.{resname} has invalid leaflet '{leaf}'"
                        )
                    if type(count) is not int or count < 0:
                        errors.append(
                            f"cholesterol.composition.lipid_targets.{resname}.{leaf} must be a non-negative integer"
                        )
            errors.extend(_repack_errors(ch.repack))
            if self.alignment.atoms not in (None, "", "CA", "backbone", "all"):
                errors.append("legacy alignment.atoms must be one of: CA, backbone, all")
            if self.alignment.method not in ("pymol_align", "pymol_cealign", "mask_fit"):
                errors.append("alignment.method must be one of: pymol_align, pymol_cealign, mask_fit")
        if type(self.alignment.cycles) is not int or self.alignment.cycles < 0:
            errors.append("alignment.cycles must be an integer >= 0")
        if not _is_finite_number(self.alignment.cutoff) or self.alignment.cutoff <= 0:
            errors.append("alignment.cutoff must be > 0")
        if not _is_finite_number(self.clash_detection.threshold) or self.clash_detection.threshold <= 0:
            errors.append("clash_detection.threshold must be > 0")
        for molclass, threshold in self.clash_detection.thresholds.items():
            if molclass not in {"protein", "water", "ion", "lipid", "other"}:
                errors.append(
                    f"clash_detection.thresholds has unrecognized class '{molclass}'"
                )
            if (
                not _is_finite_number(threshold)
                or threshold <= 0
            ):
                errors.append(
                    f"clash_detection.thresholds.{molclass} must be > 0"
                )
        if type(self.clash_detection.heavy_atoms_only) is not bool:
            errors.append("clash_detection.heavy_atoms_only must be boolean")
        if (
            not _is_finite_number(self.clash_detection.flag_lipid_removal_fraction)
            or not 0.0 <= self.clash_detection.flag_lipid_removal_fraction <= 1.0
        ):
            errors.append(
                "clash_detection.flag_lipid_removal_fraction must be between 0 and 1"
            )
        valid_classes = {"protein", "water", "ion", "lipid", "other"}
        nonstring_classes = [
            value for value in self.clash_detection.keep_classes
            if not isinstance(value, str)
        ]
        bad_classes = (
            set(self.clash_detection.keep_classes) - valid_classes
            if not nonstring_classes else set()
        )
        if nonstring_classes:
            errors.append("clash_detection.keep_classes must contain only class-name strings")
        elif bad_classes:
            errors.append(
                f"clash_detection.keep_classes has unrecognized class name(s) {sorted(bad_classes)}; "
                f"valid values are {sorted(valid_classes)}"
            )
        for label, value in (
            ("charge.exclusion_distance_from_protein", self.charge.exclusion_distance_from_protein),
            ("charge.exclusion_distance_from_ligand", self.charge.exclusion_distance_from_ligand),
            ("charge.exclusion_distance_from_lipid", self.charge.exclusion_distance_from_lipid),
        ):
            if (
                not _is_finite_number(value)
                or value < 0
            ):
                errors.append(f"{label} must be >= 0")
        if (
            not _is_finite_number(self.charge.tolerance)
            or self.charge.tolerance <= 0
        ):
            errors.append("charge.tolerance must be > 0")
        if not _is_finite_number(self.charge.target_net_charge):
            errors.append("charge.target_net_charge must be a finite number")
        if type(self.charge.neutralize) is not bool:
            errors.append("charge.neutralize must be boolean")
        if type(self.charge.exclude_membrane_interior) is not bool:
            errors.append("charge.exclude_membrane_interior must be boolean")
        if type(self.charge.random_seed) is not int or self.charge.random_seed < 0:
            errors.append("charge.random_seed must be an integer >= 0")
        if self.charge.membrane_z_override is not None:
            errors.append(
                "charge.membrane_z_override is no longer supported because absolute-z "
                "membrane slabs are unsafe for periodic coordinates; use "
                "charge.exclusion_distance_from_lipid instead"
            )
        if type(self.topology.enabled) is not bool:
            errors.append("topology.enabled must be boolean")
        elif self.topology.enabled:
            if (
                not isinstance(self.topology.protein_toppar_dir, str)
                or not self.topology.protein_toppar_dir.strip()
            ):
                errors.append("topology.protein_toppar_dir is required when topology.enabled is true")
            elif check_paths and not os.path.isdir(self.topology.protein_toppar_dir):
                errors.append(f"topology.protein_toppar_dir does not exist: {self.topology.protein_toppar_dir}")
            if (
                not isinstance(self.topology.environment_toppar_dir, str)
                or not self.topology.environment_toppar_dir.strip()
            ):
                errors.append("topology.environment_toppar_dir is required when topology.enabled is true")
            elif check_paths and not os.path.isdir(self.topology.environment_toppar_dir):
                errors.append(f"topology.environment_toppar_dir does not exist: {self.topology.environment_toppar_dir}")
            if not isinstance(self.topology.output_dir, str) or not self.topology.output_dir.strip():
                errors.append("topology.output_dir must be a non-empty path string")
            if self.topology.protein_template_top and (
                not isinstance(self.topology.protein_template_top, str)
                or (check_paths and not os.path.isfile(self.topology.protein_template_top))
            ):
                errors.append(
                    "topology.protein_template_top does not exist or is not a file: "
                    f"{self.topology.protein_template_top}"
                )
        if self.minimization.enabled and not self.topology.enabled:
            errors.append(
                "minimization.enabled requires topology.enabled so the freshly "
                "assembled coordinates have a matching topol.top"
            )
        for index, path in enumerate(self.topology.ligand_itp_paths, 1):
            if not isinstance(path, str) or not path.strip():
                errors.append(f"topology.ligand_itp_paths[{index}] must be a non-empty path string")
            elif check_paths and not os.path.isfile(path):
                errors.append(f"topology.ligand_itp_paths[{index}] does not exist: {path}")
        for resname, mtype in self.topology.moleculetype_overrides.items():
            if not isinstance(resname, str) or not isinstance(mtype, str) or not resname or not mtype:
                errors.append("topology.moleculetype_overrides must map non-empty strings to non-empty strings")
        for resname, charge in self.charge.charge_table_overrides.items():
            if not isinstance(resname, str) or not resname:
                errors.append("charge.charge_table_overrides keys must be non-empty residue names")
            if charge is not None and not _is_finite_number(charge):
                errors.append(
                    f"charge.charge_table_overrides.{resname} must be a finite number or null"
                )
        output_values = {
            "output.pdb_path": self.output.pdb_path,
            "output.gro_path": self.output.gro_path,
            "output.inspection_pdb_path": self.output.inspection_pdb_path,
            "output.report_path": self.output.report_path,
        }
        for label, value in output_values.items():
            if not isinstance(value, str) or not value.strip():
                errors.append(f"{label} must be a non-empty path string")
        concrete_outputs = {}
        if all(isinstance(value, str) and value.strip() for value in output_values.values()):
            concrete_outputs = {
                "output.pdb_path": self.output.pdb_path,
                "output.gro_path": self.output.gro_path,
                "output.inspection_pdb_path": self.output.inspection_pdb_path,
                "output.report_path (.txt)": self.output.report_path + ".txt",
                "output.report_path (.json)": self.output.report_path + ".json",
            }
        if type(self.ndx.enabled) is not bool:
            errors.append("ndx.enabled must be boolean")
        if self.ndx.style not in ("generic", "charmm_gui", "both"):
            errors.append("ndx.style must be one of: generic, charmm_gui, both")
        if self.ndx.enabled and (
            not isinstance(self.ndx.output_path, str) or not self.ndx.output_path.strip()
        ):
            errors.append("ndx.output_path must be a non-empty path string when enabled")
        elif self.ndx.enabled and concrete_outputs:
            concrete_outputs["ndx.output_path"] = self.ndx.output_path
        for group_name, mask in self.ndx.extra_groups.items():
            if (
                not isinstance(group_name, str) or not group_name.strip()
                or not isinstance(mask, str) or not mask.strip()
            ):
                errors.append("ndx.extra_groups must map non-empty names to non-empty masks")
        if (
            self.topology.enabled
            and isinstance(self.topology.output_dir, str)
            and self.topology.output_dir.strip()
            and concrete_outputs
        ):
            concrete_outputs["topology topol.top"] = os.path.join(
                self.topology.output_dir, "topol.top"
            )
        if self.mode == "chl" and self.cholesterol.write_diagnostics:
            diagnostics = {
                "cholesterol.raw_protonated_pdb_path": self.cholesterol.raw_protonated_pdb_path,
                "cholesterol.converted_pdb_path": self.cholesterol.converted_pdb_path,
                "cholesterol.aligned_reference_pdb_path": self.cholesterol.aligned_reference_pdb_path,
                "cholesterol.reference_placed_pdb_path": self.cholesterol.reference_placed_pdb_path,
                "cholesterol.merged_pdb_path": self.cholesterol.merged_pdb_path,
            }
            for label, value in diagnostics.items():
                if not isinstance(value, str) or not value.strip():
                    errors.append(f"{label} must be a non-empty path string")
                elif concrete_outputs:
                    concrete_outputs[label] = value
        if concrete_outputs:
            by_realpath = {}
            for label, value in concrete_outputs.items():
                by_realpath.setdefault(os.path.realpath(value), []).append(label)
            for labels in by_realpath.values():
                if len(labels) > 1:
                    errors.append("Configured output paths collide: " + ", ".join(labels))
        input_paths = []
        if self.mode == "protein":
            input_paths = [self.target_box.path, self.replacement_structure.path]
        elif self.mode == "lig":
            input_paths = [
                self.ligand_replace.structure_path,
                self.ligand_replace.new_ligand.coord_path,
            ]
        elif self.mode == "chl":
            input_paths = [
                self.cholesterol.experimental_structure_path,
                self.cholesterol.target_system_path,
                self.cholesterol.composition.reference_system_path,
            ]
        safe_input_overwrites = {"output.pdb_path", "output.gro_path"}
        for label, output_path in concrete_outputs.items():
            if label in safe_input_overwrites:
                continue
            output_real = os.path.realpath(output_path)
            for value in input_paths:
                if isinstance(value, str) and value and output_real == os.path.realpath(value):
                    errors.append(
                        f"{label} must not overwrite an input file"
                    )
        if self.refgro is not None and not isinstance(self.refgro, str):
            errors.append("refgro must be a path string or null")
        elif check_paths and self.refgro and not os.path.exists(self.refgro):
            errors.append(f"refgro does not exist: {self.refgro}")
        if errors:
            raise ConfigError("Invalid configuration:\n  - " + "\n  - ".join(errors))


def _repack_errors(repack) -> List[str]:
    errors = []
    label = "cholesterol.repack"
    if type(repack.enabled) is not bool:
        errors.append(f"{label}.enabled must be true or false")
    if not isinstance(repack.output_dir, str) or not re.fullmatch(r"[A-Za-z0-9_.-]+", repack.output_dir or ""):
        errors.append(f"{label}.output_dir must be a simple folder name")
    if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]{0,15}", repack.moleculetype or ""):
        errors.append(f"{label}.moleculetype must be a short name (letters, digits, _)")
    for key in ("lipid_restraint", "temperature"):
        value = getattr(repack, key)
        if isinstance(value, bool) or not isinstance(value, (int, float)) or value < 0:
            errors.append(f"{label}.{key} must be a number >= 0")
    names = []
    for number, stage in enumerate(repack.stages):
        where = f"{label}.stages[{number}]"
        if not isinstance(stage.name, str) or not re.fullmatch(r"[A-Za-z0-9_.-]+", stage.name):
            errors.append(f"{where}.name must be a simple file name")
        names.append(stage.name)
        if stage.integrator not in ("steep", "nvt", "npt"):
            errors.append(f"{where}.integrator must be steep, nvt or npt")
        if type(stage.nsteps) is not int or stage.nsteps <= 0:
            errors.append(f"{where}.nsteps must be a whole number > 0")
        for key in ("dt", "restraint"):
            value = getattr(stage, key)
            if isinstance(value, bool) or not isinstance(value, (int, float)) or value < 0:
                errors.append(f"{where}.{key} must be a number >= 0")
        if stage.integrator != "steep" and not (isinstance(stage.dt, (int, float)) and 0 < stage.dt <= 0.004):
            errors.append(f"{where}.dt must be between 0 and 0.004 ps")
    if len(set(names)) != len(names):
        errors.append(f"{label}.stages names must be unique")
    return errors
