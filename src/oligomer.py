"""Add copies of the membrane protein to the same membrane (dimers, oligomers).

Each copy of the host's proteins (and, by default, their ligands) is rotated
about the membrane normal through the protein centre, optionally flipped
(antiparallel) about an in-plane axis through the membrane midplane, and
shifted in the membrane plane. Lipids overlapping a copy are removed whole,
water and ions as needed, and counterions restore the charge. Every copy must
stay `min_protein_gap` from every other copy and from all periodic images.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field

import numpy as np

from charmprot import _heavy, _indices
from config import ConfigError
from equilibration import EquilibrationSpec
from system_edit import (MEMBRANE_CATEGORIES, SOLVENT_CATEGORIES, balance_charge, ion_block, load_section,
                         molecules_near, monovalent_ions, net_charge, order, place, positive, run_edit,
                         salt_summary, unwrap, write_system)


def _no_equilibration():
    return EquilibrationSpec(enabled=False)


@dataclass
class CopySpec:
    shift: list[float] | None = None     # x y (A) in the membrane plane; empty: half the box along x
    rotation: float = 0.0                # degrees about the membrane normal through the protein centre
    flip: bool = False                   # turn upside down (antiparallel), about the x axis


def _default_copies():
    return [CopySpec(rotation=180.0)]


@dataclass
class OligomerSpec:
    host: str
    host_gro: str | None = None
    copies: list[CopySpec] = field(default_factory=_default_copies)
    include_ligands: bool = True
    min_protein_gap: float = 10.0        # A between heavy atoms of different copies, periodic images included
    lipid_clash_distance: float = 3.0    # A; lipids closer to a copy are removed whole
    water_clash_distance: float = 2.5    # A; water/ions closer to a copy are removed
    neutralize: bool = True
    target_net_charge: float = 0.0
    seed: int = 42
    environment_overrides: dict = field(default_factory=dict)
    output_dir: str = "oligomer_output"
    equilibration: EquilibrationSpec = field(default_factory=_no_equilibration)


def load_config(path, check_paths=True, output_root=None):
    cfg = load_section(path, "oligomer", OligomerSpec, check_paths, output_root)
    positive(cfg, "oligomer", "min_protein_gap", "lipid_clash_distance", "water_clash_distance")
    if not cfg.copies:
        raise ConfigError("oligomer.copies must list at least one copy")
    for number, copy in enumerate(cfg.copies):
        if copy.shift is not None and len(copy.shift) != 2:
            raise ConfigError(f"oligomer.copies[{number}].shift must be two numbers: x y in A")
    return cfg


def _rotation_z(degrees):
    a = np.radians(degrees)
    return np.array([[np.cos(a), -np.sin(a), 0.0], [np.sin(a), np.cos(a), 0.0], [0.0, 0.0, 1.0]])


def _min_gap(a, b, box, cutoff):
    """Smallest distance between two point sets under PBC, or None when above cutoff."""
    from MDAnalysis.lib.distances import capped_distance
    _, distances = capped_distance(a, b, cutoff, box=box)
    return float(distances.min()) if len(distances) else None


def run(config_path, dry_run=False, output_root=None):
    cfg = load_config(config_path, check_paths=not dry_run, output_root=output_root)

    def build(stage, host, report):
        from addbinder import membrane_frame
        work = stage.parent / "work"
        work.mkdir()
        universe = host.universe
        box = universe.dimensions
        lengths = box[:3].astype(float)
        frame = membrane_frame(universe, host.molecules)
        midplane = float(np.mod(lengths[2] / 2 - frame.shift, lengths[2]))
        categories = ("protein", "ligand") if cfg.include_ligands else ("protein",)
        unit = [m for m in host.molecules if m.category in categories]
        # The unit made whole as one body: each molecule whole, then moved next to the first protein.
        blocks, anchor = [], None
        for mol in unit:
            block = universe.atoms[mol.start:mol.stop].positions.astype(float)
            block = unwrap(block, lengths) if len(block) > 1 else block
            centre = block.mean(axis=0)
            if anchor is None:
                anchor = centre
            delta = centre - anchor
            block = block - lengths * np.round(delta / lengths)
            blocks.append(block)
        protein_heavy_mask = [np.array([not a.atom_name.upper().startswith(("H", "LP")) for a in m.definition.atoms])
                              for m in unit]
        unit_heavy = np.concatenate([b[mask] for b, mask in zip(blocks, protein_heavy_mask)])
        pivot = unit_heavy.mean(axis=0)
        pivot[2] = midplane
        new_units, unit_reports = [], []
        heavy_sets = [unit_heavy]
        for number, copy in enumerate(cfg.copies):
            shift = np.array([*(copy.shift if copy.shift is not None else (lengths[0] / 2, 0.0)), 0.0])
            matrix = _rotation_z(copy.rotation)
            if copy.flip:
                matrix = matrix @ np.diag([1.0, -1.0, -1.0])
            moved = [(b - pivot) @ matrix.T + pivot + shift for b in blocks]
            heavy = np.concatenate([b[mask] for b, mask in zip(moved, protein_heavy_mask)])
            gaps = [g for other in heavy_sets
                    for g in [_min_gap(heavy, other, box, cfg.min_protein_gap)] if g is not None]
            if gaps:
                raise ConfigError(f"Copy {number + 1} comes {min(gaps):.1f} A from another copy or a periodic "
                                  f"image (min_protein_gap {cfg.min_protein_gap:g} A); change its shift or rotation")
            heavy_sets.append(heavy)
            new_units.append(moved)
            unit_reports.append({"shift_angstrom": [round(float(v), 3) for v in shift[:2]],
                                 "rotation_degrees": copy.rotation, "flip": copy.flip})
        # Closest approach between copies, including periodic images (for the report).
        closest = None
        for i in range(len(heavy_sets)):
            for j in range(i + 1, len(heavy_sets)):
                g = _min_gap(heavy_sets[i], heavy_sets[j], box, 60.0)
                if g is not None:
                    closest = g if closest is None else min(closest, g)
        new_heavy = np.concatenate(heavy_sets[1:])
        new_all = np.concatenate([b for moved in new_units for b in moved])
        membrane = [m for m in host.molecules if m.category in MEMBRANE_CATEGORIES]
        clashing = molecules_near(host, membrane, new_heavy, cfg.lipid_clash_distance, box)
        removed_lipids = [m for i, m in enumerate(membrane) if i in clashing]
        membrane = [m for i, m in enumerate(membrane) if i not in clashing]
        solvent = [m for m in host.molecules if m.category in SOLVENT_CATEGORIES]
        clashing = molecules_near(host, solvent, new_all, cfg.water_clash_distance, box)
        removed_solvent = [m for i, m in enumerate(solvent) if i in clashing]
        solvent = [m for i, m in enumerate(solvent) if i not in clashing]
        others = [m for m in host.molecules if m.category not in MEMBRANE_CATEGORIES + SOLVENT_CATEGORIES]
        items = [(m.definition, m.category, b) for moved in new_units for m, b in zip(unit, moved)]
        block = place(items, work, label="oligomer")
        new = [(block, m) for m in block.molecules]
        print(f"[oligomer] {len(cfg.copies)} cop{'y' if len(cfg.copies) == 1 else 'ies'} of "
              + "+".join(m.name for m in unit) + (f"; closest approach {closest:.1f} A" if closest else "") + "; removed "
              + ", ".join(f"{n} {k}" for k, n in Counter(m.name for m in removed_lipids + removed_solvent).items()))
        net = net_charge([(host, m) for m in others + membrane + solvent]) + sum(d.charge for d, _, _ in items)
        report["charge"] = {"before": round(net, 6), "after": round(net, 6), "removed_ions": [], "added_ions": {}}
        added = []
        if cfg.neutralize:
            cation, anion = monovalent_ions(host)
            if cation is None or anion is None:
                raise ConfigError("The host has no monovalent cation and anion to neutralize with")
            protein_heavy = np.concatenate([_heavy(universe.atoms[_indices(unit)]).positions, new_heavy])
            solvent, added, report["charge"] = balance_charge(host, solvent, net, cfg.target_net_charge, cation,
                                                              anion, [(protein_heavy, 8.0)], box, seed=cfg.seed)
            print(f"[charge] {report['charge']['before']:+.3f} e -> {report['charge']['after']:+.3f} e")
        new += ion_block(added, work)
        ordered = order(host, others + membrane + solvent, new=new)
        leaflets = Counter()
        for m in removed_lipids:
            z = np.mod(universe.atoms[m.start:m.stop].positions[:, 2].mean() + frame.shift, lengths[2])
            leaflets[("upper" if z > lengths[2] / 2 else "lower", m.name)] += 1
        report["oligomer"] = {"unit": [m.name for m in unit], "copies": unit_reports,
                              "closest_approach_angstrom": round(closest, 3) if closest is not None else None,
                              "removed_lipids": {f"{side} {name}": n for (side, name), n in sorted(leaflets.items())},
                              "removed_solvent": dict(Counter(m.name for m in removed_solvent))}
        report["salt"] = {"host": salt_summary([(host, m) for m in host.molecules]), "final": salt_summary(ordered)}
        copy_ids = {id(m) for m in block.molecules}
        report.update(write_system(
            stage, host, ordered, "oligomer", equilibration=cfg.equilibration,
            purpose="Equilibration of the oligomer: the lipids close the gaps around the new copies",
            extra_index={"Copies": copy_ids, "Original": {id(m) for m in unit}},
            scene=({"inserted_copy": copy_ids}, {"inserted_copy": ("Added copies", [.85, .35, .75], True)},
                   "The original protein keeps its coordinates; the copies are new.")))

    def lines(report):
        info = report.get("oligomer")
        if not info:
            return []
        out = ["", f"OLIGOMER: unit {'+'.join(info['unit'])}; closest approach between copies (periodic images "
                   f"included) {info['closest_approach_angstrom']} A"]
        for number, copy in enumerate(info["copies"], 1):
            out.append(f"  copy {number}: shift {copy['shift_angstrom']} A, rotation {copy['rotation_degrees']:g} deg"
                       + (", flipped" if copy["flip"] else ""))
        out.append("  lipids removed: " + (", ".join(f"{k} {n}" for k, n in info["removed_lipids"].items()) or "none"))
        out.append("  solvent removed: " + (", ".join(f"{k} {n}" for k, n in info["removed_solvent"].items()) or "none"))
        return out

    return run_edit("oligomer", cfg, "ADD PROTEIN COPIES TO THE MEMBRANE", build, "add protein copies", lines,
                    dry_run=dry_run)
