"""Change residue protonation states in a prepared system.

Histidine HSD/HSE/HSP, aspartate ASP/ASPP, glutamate GLU/GLUP, lysine LYS/LSN
and cysteine CYS/CYM. The protein ITP gets the CHARMM36 atoms, types and
charges of the new state with its bonded terms regenerated (see protein_edit);
new hydrogens are placed from the CHARMM internal coordinates and removed ones
dropped; every other atom keeps its coordinates. The state of each residue is
read from its atoms (CHARMM-GUI keeps the name ASP for a protonated ASPP).
"""
from __future__ import annotations

from dataclasses import dataclass, field

from config import ConfigError
from equilibration import EquilibrationSpec
from protein_edit import VARIANTS, replace_residue, residue_state
from residue_edit import report_lines, run_changes
from system_edit import load_section, positive, run_edit

STATES = ("HSD", "HSE", "HSP", "ASP", "ASPP", "GLU", "GLUP", "LYS", "LSN", "CYS", "CYM")


def _no_equilibration():
    return EquilibrationSpec(enabled=False)


@dataclass
class ResidueChange:
    molecule: str = "PROA"     # protein moleculetype
    resid: int = 1             # residue number as in the ITP
    state: str = "HSP"


def _default_residues():
    return [ResidueChange()]


@dataclass
class ProtonateSpec:
    host: str
    host_gro: str | None = None
    residues: list[ResidueChange] = field(default_factory=_default_residues)
    water_clash_distance: float = 2.0    # A, new atoms to water/ion heavy atoms: those are removed
    neutralize: bool = True
    target_net_charge: float = 0.0
    seed: int = 42
    environment_overrides: dict = field(default_factory=dict)
    output_dir: str = "protonate_output"
    equilibration: EquilibrationSpec = field(default_factory=_no_equilibration)


def load_config(path, check_paths=True, output_root=None):
    cfg = load_section(path, "protonate", ProtonateSpec, check_paths, output_root)
    positive(cfg, "protonate", "water_clash_distance")
    if not cfg.residues:
        raise ConfigError("protonate.residues must list at least one residue")
    for number, change in enumerate(cfg.residues):
        if change.state not in STATES:
            raise ConfigError(f"protonate.residues[{number}].state must be one of {', '.join(STATES)}")
    return cfg


def run(config_path, dry_run=False, output_root=None):
    cfg = load_config(config_path, check_paths=not dry_run, output_root=output_root)

    def apply(itp, positions, resid, target, environment, mol):
        names = [a.name for a in itp.atoms if a.resnr == resid]
        resname = next(a.resname for a in itp.atoms if a.resnr == resid)
        state = residue_state(resname, names)
        if state is None or target not in VARIANTS.get(state, ()):
            raise ConfigError(f"{mol.name} {resname} {resid} cannot become {target}: a protonation change keeps the "
                              "residue (use the mutate mode to change it)")
        new_positions, change = replace_residue(itp, resid, target, positions,
                                                host_parameter_files=[str(p) for p in host_files])
        change["to_label"] = "" if state != target else " (unchanged)"
        return new_positions, change

    host_files = []

    def build(stage, host, report):
        host_files.extend(host.files)
        changes = [(c.molecule, c.resid, c.state) for c in cfg.residues]
        run_changes(stage, host, report, cfg, "protonate", changes, apply,
                    purpose="Equilibration with the new protonation states")

    return run_edit("protonate", cfg, "CHANGE PROTONATION STATES", build, "change protonation states",
                    report_lines("protonate", "PROTONATION STATES"), dry_run=dry_run)
