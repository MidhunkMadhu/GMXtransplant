"""Flood the water of a prepared system with copies of a small molecule.

Copies of one molecule (a ligand for unbiased binding simulations, or a
cosolvent probe for mixed-solvent mapping) are placed at random positions and
orientations in the water, at least `min_solute_distance` from the protein and
ligands, `min_membrane_distance` from the lipids and `min_copy_distance` from
each other, through every periodic image. Water and ions overlapping a copy are
removed, and counterions restore the net charge.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.spatial.transform import Rotation

from charmprot import _heavy, _indices
from config import ConfigError
from equilibration import EquilibrationSpec
from system_edit import (MEMBRANE_CATEGORIES, SOLVENT_CATEGORIES, WATER_MOLARITY, _near, balance_charge, ion_block,
                         load_section, molecules_near, monovalent_ions, net_charge, order, place, positive, run_edit,
                         salt_summary, write_system)
from topology import TopologyError, check_parameters


def _no_equilibration():
    return EquilibrationSpec(enabled=False)


@dataclass
class FloodSpec:
    host: str
    molecule_coordinates: str = ""
    molecule_itp: str = ""
    molecule_forcefield: str | None = None
    host_gro: str | None = None
    copies: int = 10                     # 0: use concentration
    concentration: float = 0.0           # mol/L in the water (copies / waters x 55.51)
    side: str = "both"                   # both | upper | lower water layer
    min_solute_distance: float = 10.0    # A, copy heavy atoms to protein/ligand heavy atoms
    min_membrane_distance: float = 5.0   # A, copy heavy atoms to lipid heavy atoms
    min_copy_distance: float = 8.0       # A between the heavy atoms of two copies
    water_clash_distance: float = 2.5    # A, copy atoms to water/ion heavy atoms: those are removed
    attempts: int = 20000
    neutralize: bool = True
    target_net_charge: float = 0.0
    seed: int = 42
    environment_overrides: dict = field(default_factory=dict)
    output_dir: str = "flood_output"
    equilibration: EquilibrationSpec = field(default_factory=_no_equilibration)


def load_config(path, check_paths=True, output_root=None):
    cfg = load_section(path, "flood", FloodSpec, check_paths, output_root,
                       optional_paths=("host_gro", "molecule_coordinates", "molecule_itp", "molecule_forcefield"))
    for key in ("molecule_coordinates", "molecule_itp"):
        if not getattr(cfg, key):
            raise ConfigError(f"flood.{key} must be a file path")
    positive(cfg, "flood", "min_copy_distance", "water_clash_distance")
    positive(cfg, "flood", "min_solute_distance", "min_membrane_distance", "concentration", allow_zero=True)
    if cfg.copies < 0 or cfg.attempts <= 0:
        raise ConfigError("flood.copies must be 0 or more and flood.attempts more than 0")
    if cfg.copies == 0 and cfg.concentration <= 0:
        raise ConfigError("Set flood.copies, or flood.concentration with copies: 0")
    if cfg.side not in ("both", "upper", "lower"):
        raise ConfigError("flood.side must be both, upper or lower")
    return cfg


def _random_points(rng, dims, count):
    from MDAnalysis.lib.mdamath import triclinic_vectors
    vectors = triclinic_vectors(dims)
    return rng.random((count, 3)) @ vectors


def run(config_path, dry_run=False, output_root=None):
    cfg = load_config(config_path, check_paths=not dry_run, output_root=output_root)

    def build(stage, host, report):
        from addbinder import load_binder, membrane_frame
        work = stage.parent / "work"
        work.mkdir()
        molecule = load_binder(SimpleNamespace(binder_itp=[cfg.molecule_itp],
                                               binder_coordinates=cfg.molecule_coordinates,
                                               binder_forcefield=cfg.molecule_forcefield), work, host)
        definition = molecule.definitions[0]
        forcefields = [str(p) for p in molecule.files[:-1]] + [str(p) for p in host.files]
        try:
            check_parameters([definition.source_path], forcefields,
                             "Set molecule_forcefield to the forcefield.itp generated with the molecule's ITP.")
        except TopologyError as exc:
            raise ConfigError(str(exc)) from exc
        print(f"[flood] molecule {definition.name}: {definition.atom_count} atoms, charge {definition.charge:+.3f} e")
        universe = host.universe
        box = universe.dimensions
        protected = [m for m in host.molecules if m.category in ("protein", "ligand")]
        membrane = [m for m in host.molecules if m.category in MEMBRANE_CATEGORIES]
        solvent = [m for m in host.molecules if m.category in SOLVENT_CATEGORIES]
        waters = sum(m.category == "water" for m in solvent)
        count = cfg.copies or int(round(cfg.concentration * waters / WATER_MOLARITY))
        if count <= 0:
            raise ConfigError(f"{cfg.concentration:g} M in {waters} waters is less than one copy")
        protected_heavy = _heavy(universe.atoms[_indices(protected)]).positions
        membrane_heavy = _heavy(universe.atoms[_indices(membrane)]).positions if membrane else np.zeros((0, 3))
        heavy_mask = np.array([not a.atom_type.upper().startswith("H") for a in definition.atoms])
        local = molecule.positions - molecule.positions[heavy_mask].mean(axis=0)
        frame = membrane_frame(universe, host.molecules) if membrane else None
        rng = np.random.default_rng(cfg.seed)
        placed, placed_heavy, tries = [], np.zeros((0, 3)), 0
        while len(placed) < count and tries < cfg.attempts:
            tries += 1
            centre = _random_points(rng, box, 1)[0]
            if frame is not None and cfg.side != "both":
                z = np.mod(centre[2] + frame.shift, frame.lengths[2])
                if (cfg.side == "upper") != (z > frame.lengths[2] / 2):
                    continue
            coords = Rotation.random(random_state=rng).apply(local) + centre
            heavy = coords[heavy_mask]
            if _near(heavy, protected_heavy, cfg.min_solute_distance, box):
                continue
            if len(membrane_heavy) and _near(heavy, membrane_heavy, cfg.min_membrane_distance, box):
                continue
            if len(placed_heavy) and _near(heavy, placed_heavy, cfg.min_copy_distance, box):
                continue
            placed.append(coords)
            placed_heavy = np.concatenate([placed_heavy, heavy])
        if len(placed) < count:
            raise ConfigError(f"Placed only {len(placed)} of {count} copies in {cfg.attempts} attempts; lower the "
                              "distances or the number of copies, or raise attempts")
        all_copy_atoms = np.concatenate(placed)
        clashing = molecules_near(host, solvent, all_copy_atoms, cfg.water_clash_distance, box)
        removed = [m for i, m in enumerate(solvent) if i in clashing]
        solvent = [m for i, m in enumerate(solvent) if i not in clashing]
        print(f"[flood] {len(placed)} copies placed in {tries} attempts; removed "
              + ", ".join(f"{n} {k}" for k, n in Counter(m.name for m in removed).items()))
        block = place([(definition, "ligand", p) for p in placed], work, label="flood")
        new = [(block, m) for m in block.molecules]
        net = net_charge([(host, m) for m in protected + membrane + solvent]) + len(placed) * definition.charge
        added = []
        report["charge"] = {"before": round(net, 6), "after": round(net, 6), "removed_ions": [], "added_ions": {}}
        if cfg.neutralize:
            cation, anion = monovalent_ions(host)
            if cation is None or anion is None:
                raise ConfigError("The host has no monovalent cation and anion to neutralize with")
            avoid = [(protected_heavy, 5.0), (placed_heavy, 5.0)]
            solvent, added, report["charge"] = balance_charge(host, solvent, net, cfg.target_net_charge, cation,
                                                              anion, avoid, box, seed=cfg.seed + 3)
            print(f"[charge] {report['charge']['before']:+.3f} e -> {report['charge']['after']:+.3f} e")
        new += ion_block(added, work)
        ordered = order(host, protected + membrane + solvent, new=new)
        final_waters = sum(m.category == "water" for _, m in ordered)
        report["flood"] = {"moleculetype": definition.name, "copies": len(placed), "attempts": tries,
                           "molecule_charge": definition.charge, "renamed_from": molecule.report["renamed_from"],
                           "concentration_molar": round(len(placed) / final_waters * WATER_MOLARITY, 4),
                           "removed_solvent": dict(Counter(m.name for m in removed))}
        report["salt"] = {"host": salt_summary([(host, m) for m in host.molecules]), "final": salt_summary(ordered)}
        flood_ids = {id(m) for m in block.molecules}
        report.update(write_system(
            stage, host, ordered, "flood", parameter_files=[Path(p) for p in molecule.files[:-1]],
            equilibration=cfg.equilibration, purpose=f"Equilibration of the system with {len(placed)} "
                                                     f"{definition.name} copies",
            extra_index={f"Flood_{definition.name}": flood_ids},
            scene=({"inserted_flood": flood_ids},
                   {"inserted_flood": (f"{definition.name} copies", [.95, .45, .1], True)},
                   "Host coordinates unchanged; the copies are new."),
            restrained=[m.definition for m in protected]))

    def lines(report):
        info = report.get("flood")
        if not info:
            return []
        return ["", f"FLOOD: {info['copies']} {info['moleculetype']} copies ({info['concentration_molar']} M in the "
                    f"water), placed in {info['attempts']} attempts",
                "  solvent removed: " + (", ".join(f"{k} {n}" for k, n in info["removed_solvent"].items()) or "none"),
                f"  index group Flood_{info['moleculetype']} holds every copy"]

    return run_edit("flood", cfg, "FLOOD WITH LIGAND COPIES", build, "flood with ligand copies", lines,
                    dry_run=dry_run)
