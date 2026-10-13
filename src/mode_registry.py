"""Metadata of the modes that edit one prepared system (no heavy imports).

The command line, the configuration path rules and the desktop application read
this table, so a system-edit mode is declared once: its YAML section, module,
path fields, example and GUI texts.
"""
from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class EditMode:
    key: str                  # --mode value and YAML section
    title: str                # mode name in the GUI
    module: str               # module with load_config(path, check_paths, output_root) and run(path, dry_run, output_root)
    spec: str                 # dataclass name of the section in `module`
    summary: str              # one line for --help and the README
    description: str          # GUI mode description
    example_dir: str          # examples/<example_dir>/
    example_yaml: str
    example_title: str
    example_description: str
    progress: dict = field(default_factory=dict)   # "[tag]" -> GUI progress label
    directories: tuple = ("host",)                 # input folders
    files: tuple = ("host_gro",)                   # input files
    file_lists: tuple = ()                         # lists of input files
    required: tuple = ("host",)
    primary: tuple = ("host",)                     # fields shown first in the editor
    labels: dict = field(default_factory=dict)     # field -> GUI label
    help: dict = field(default_factory=dict)       # field -> GUI tooltip
    choices: dict = field(default_factory=dict)    # field path tuple (after the section) -> options
    section_title: str = "Host system"
    section_help: str = ("Choose the prepared membrane system: a CHARMM-GUI-style folder with topol.top, "
                         "toppar/ and the step5 GRO.")


_HOST_LABELS = {"host": "Host system folder", "host_gro": "Host coordinates (GRO)",
                "environment_overrides": "Molecule classification overrides",
                "equilibration": "Equilibration (GROMACS)",
                "equilibration.enabled": "Write the equilibration folder",
                "equilibration.output_dir": "Folder name", "equilibration.temperature": "Temperature (K)",
                "equilibration.stages": "Stages"}
_HOST_HELP = {
    "host": "The prepared system, e.g. a CHARMM-GUI build: a folder with topol.top, toppar/ and step5_input.gro.",
    "host_gro": "Tick to use a coordinate file other than the host folder's own step5 GRO, e.g. an MD frame.",
    "environment_overrides": "Classify molecule types the automatic detector does not recognise, for example "
                             "a custom lipid MYLIP → lipid.",
    "equilibration": "A GROMACS folder with the CHARMM-GUI step6 schedule, a sample run script and an analysis "
                     "script. Nothing is run here.",
}


def _labels(**extra):
    return {**_HOST_LABELS, **extra}


def _help(**extra):
    return {**_HOST_HELP, **extra}


EDIT_MODES = {
    "salt": EditMode(
        key="salt", title="Salt exchange", module="salt", spec="SaltSpec",
        summary="Set the salt type and concentration.",
        description="Set the salt of a prepared system: remove the bulk ions and add the chosen cation and anion "
                    "(NaCl, KCl, CaCl2, MgCl2 ...) at a concentration in mol/L, replacing random bulk waters "
                    "away from the protein and membrane, then neutralize. Missing ion parameters come from the "
                    "bundled CHARMM36 set.",
        example_dir="salt", example_yaml="salt.yaml",
        example_title="Salt exchange: D1R",
        example_description="150 mM KCl instead of NaCl.",
        progress={"[salt]": "Placing ions"},
        files=("host_gro",), file_lists=("ion_itps",),
        primary=("host", "concentration", "cation", "anion", "host_gro"),
        labels=_labels(concentration="Salt concentration (mol/L)", cation="Cation", anion="Anion",
                       replace_existing="Remove the existing bulk ions first",
                       keep_bound_distance="Keep ions this close to the protein (Å)",
                       ion_itps="Extra ion ITPs", solute_distance="Smallest ion distance to protein/ligands (Å)",
                       membrane_distance="Smallest ion distance to lipids (Å)",
                       ion_spacing="Smallest distance between ions (Å)"),
        help=_help(cation="Ion moleculetype: SOD, POT, CAL, MG, LIT, CES, ZN2 ... (CHARMM36 names).",
                   anion="Ion moleculetype: CLA (CHARMM36 chloride).",
                   keep_bound_distance="Ions with an atom this close to a protein heavy atom are kept "
                                       "(e.g. a structural sodium); 0 removes all bulk ions.")),
    "protonate": EditMode(
        key="protonate", title="Protonation states", module="protonate", spec="ProtonateSpec",
        summary="Change residue protonation states.",
        description="Change the protonation state of chosen residues of a prepared system: histidine HSD/HSE/HSP, "
                    "protonated ASP/GLU, neutral LYS. The protein ITP gets the CHARMM36 atoms, charges and bonded "
                    "terms, hydrogens are placed from the CHARMM internal coordinates, and counterions keep the "
                    "system neutral.",
        example_dir="protonate", example_yaml="protonate.yaml",
        example_title="Protonation states: D1R",
        example_description="D2.50 charged, H164 doubly protonated.",
        progress={"[protonate]": "Changing protonation states"},
        primary=("host", "residues", "host_gro"),
        labels=_labels(residues="Residues", **{"residues.*.molecule": "Protein moleculetype",
                                               "residues.*.resid": "Residue number",
                                               "residues.*.state": "New state"}),
        help=_help(residues="Each residue: its protein moleculetype (e.g. PROA), residue number as in the ITP, "
                            "and the new state."),
        choices={("residues", "*", "state"): ["HSD", "HSE", "HSP", "ASP", "ASPP", "GLU", "GLUP", "LYS", "LSN"]},
        section_title="Host system and residues"),
    "mutate": EditMode(
        key="mutate", title="Point mutations", module="mutate", spec="MutateSpec",
        summary="Mutate residues; side chains built in place.",
        description="Mutate residues of a prepared system. The side chain is built from the CHARMM36 internal "
                    "coordinates, its rotamer chosen to avoid clashes, the protein ITP gets the new residue's "
                    "atoms and bonded terms, overlapping water and ions are removed and the charge rebalanced.",
        example_dir="mutate", example_yaml="mutate.yaml",
        example_title="Point mutations: D1R",
        example_description="D3.32A and S5.42A.",
        progress={"[mutate]": "Building mutated residues"},
        primary=("host", "mutations", "host_gro"),
        labels=_labels(mutations="Mutations", **{"mutations.*.molecule": "Protein moleculetype",
                                                 "mutations.*.resid": "Residue number",
                                                 "mutations.*.to": "New residue"}),
        help=_help(mutations="Each mutation: protein moleculetype (e.g. PROA), residue number as in the ITP, and "
                             "the new residue (CHARMM36 name: ALA, ARG ... HSD, HSE, HSP, ASPP ...)."),
        choices={("mutations", "*", "to"): ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HSD",
                                            "HSE", "HSP", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
                                            "THR", "TRP", "TYR", "VAL"]},
        section_title="Host system and mutations"),
    "hmr": EditMode(
        key="hmr", title="Hydrogen mass repartitioning", module="hmr", spec="HmrSpec",
        summary="Heavier hydrogens for 4 fs time steps.",
        description="Repartition hydrogen masses in every molecule of a prepared system (water excluded): each "
                    "hydrogen gets the chosen mass, taken from the heavy atom it is bonded to, so the total mass "
                    "is unchanged. Allows a 4 fs time step with h-bond constraints.",
        example_dir="hmr", example_yaml="hmr.yaml",
        example_title="Hydrogen mass repartitioning: D1R",
        example_description="3.024 u hydrogens for 4 fs steps.",
        progress={"[hmr]": "Repartitioning hydrogen masses"},
        primary=("host", "hydrogen_mass", "host_gro"),
        labels=_labels(hydrogen_mass="Hydrogen mass (u)", skip_water="Leave water unchanged",
                       time_step="Time step written to the mdp files (ps)"),
        help=_help(hydrogen_mass="3.024 u is three times the hydrogen mass (CHARMM-GUI's HMR choice).")),
    "resize": EditMode(
        key="resize", title="Water layer extension", module="resize", spec="ResizeSpec",
        summary="Add bulk water to make the box taller.",
        description="Make the water layer of a prepared membrane system thicker: slabs of the host's own bulk "
                    "water and ions are copied into the new space, so the box grows along z without changing "
                    "the membrane or the salt concentration.",
        example_dir="resize", example_yaml="resize.yaml",
        example_title="Water layer extension: D1R",
        example_description="20 Å more water.",
        progress={"[resize]": "Adding water slabs"},
        primary=("host", "add_water", "host_gro"),
        labels=_labels(add_water="Water added along z (Å)",
                       slab_margin="Copied slab's distance from protein and lipids (Å)",
                       water_clash_distance="Water overlap distance at the seams (Å)"),
        help=_help(add_water="The box grows by this much along z. The water above and below the membrane is "
                             "one layer through the periodic boundary, so it all goes in the middle of that layer.")),
    "flood": EditMode(
        key="flood", title="Ligand flooding", module="flood", spec="FloodSpec",
        summary="Add many copies of a small molecule to the water.",
        description="Place many copies of a small molecule (a ligand or cosolvent probe) at random positions and "
                    "orientations in the water of a prepared system, away from the protein, the membrane and "
                    "each other. Overlapping water is removed and the charge rebalanced. Used for unbiased "
                    "binding (ligand flooding) and mixed-solvent mapping.",
        example_dir="flood", example_yaml="flood.yaml",
        example_title="Ligand flooding: D1R",
        example_description="10 dopamine copies in the water.",
        progress={"[flood]": "Placing ligand copies"},
        files=("host_gro", "molecule_coordinates", "molecule_itp", "molecule_forcefield"),
        primary=("host", "molecule_coordinates", "molecule_itp", "molecule_forcefield", "copies",
                 "concentration", "host_gro"),
        required=("host", "molecule_coordinates", "molecule_itp"),
        labels=_labels(molecule_coordinates="Molecule coordinates (PDB, GRO or MOL2)",
                       molecule_itp="Molecule ITP", molecule_forcefield="Molecule force field (forcefield.itp)",
                       copies="Number of copies", concentration="Or concentration (mol/L)",
                       min_solute_distance="Smallest distance to protein and ligands (Å)",
                       min_membrane_distance="Smallest distance to lipids (Å)",
                       min_copy_distance="Smallest distance between copies (Å)", seed="Random seed"),
        help=_help(copies="0 uses the concentration instead.")),
    "oligomer": EditMode(
        key="oligomer", title="Protein copies", module="oligomer", spec="OligomerSpec",
        summary="Add copies of the membrane protein.",
        description="Add copies of the membrane protein (with its ligands) to the same membrane: each copy is "
                    "shifted in the membrane plane and rotated about the membrane normal. Overlapping lipids "
                    "are removed whole, water and ions as needed, and the charge rebalanced.",
        example_dir="oligomer", example_yaml="oligomer.yaml",
        example_title="Receptor dimer: D1R",
        example_description="A second receptor, turned 180°.",
        progress={"[oligomer]": "Placing protein copies"},
        primary=("host", "copies", "host_gro"),
        labels=_labels(copies="Copies", **{"copies.*.shift": "Shift in the membrane plane x y (Å)",
                                           "copies.*.rotation": "Rotation about the membrane normal (degrees)",
                                           "copies.*.flip": "Flip (antiparallel copy)"},
                       lipid_clash_distance="Lipid removal distance (Å)",
                       min_protein_gap="Smallest gap between protein copies and images (Å)"),
        help=_help(copies="Each copy: shift x y in Å from the original, rotation about the membrane normal "
                          "through the protein centre.")),
}


def path_kind(section, key):
    mode = EDIT_MODES.get(section)
    if mode is None:
        return None
    if key in mode.directories or key == "output_dir":
        return "directory"
    if key in mode.files or key in mode.file_lists:
        return "file"
    return None


def section_title(key):
    return EDIT_MODES[key].section_title
