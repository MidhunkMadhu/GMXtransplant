"""Make the water layer of a prepared membrane system thicker.

The box grows along z by `add_water`. In a periodic membrane box the water above
and below the bilayer is one layer joined through the z boundary, so the new
water goes there, where the layer is farthest from the membrane: every molecule
moves so the bilayer centre sits mid-box, and the gap left at the boundary is
filled with copies of the host's own bulk water slab (water and ions, already
equilibrated at the host's density and salt). Copied waters that overlap at the
seams are removed, and copied ions are dropped as needed to keep the charge.
Lateral box vectors are unchanged.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field

import numpy as np

from config import ConfigError
from equilibration import EquilibrationSpec
from system_edit import (SOLVENT_CATEGORIES, _near, load_section, net_charge, order, place, positive, run_edit,
                         salt_summary, unwrap, write_system)

WATER_DENSITY = 0.0334   # waters per A^3 at 300 K


def _no_equilibration():
    return EquilibrationSpec(enabled=False)


@dataclass
class ResizeSpec:
    host: str
    host_gro: str | None = None
    add_water: float = 20.0              # A added to the box along z
    slab_margin: float = 4.0             # A between the copied slab and protein/lipid atoms
    water_clash_distance: float = 2.4    # A, heavy atoms of a copied water/ion to any kept heavy atom
    target_net_charge: float = 0.0
    seed: int = 42
    environment_overrides: dict = field(default_factory=dict)
    output_dir: str = "resize_output"
    equilibration: EquilibrationSpec = field(default_factory=_no_equilibration)


def load_config(path, check_paths=True, output_root=None):
    cfg = load_section(path, "resize", ResizeSpec, check_paths, output_root)
    positive(cfg, "resize", "add_water", "water_clash_distance")
    positive(cfg, "resize", "slab_margin", allow_zero=True)
    if cfg.add_water > 500:
        raise ConfigError("resize.add_water must be at most 500 A")
    return cfg


def run(config_path, dry_run=False, output_root=None):
    cfg = load_config(config_path, check_paths=not dry_run, output_root=output_root)

    def build(stage, host, report):
        from addbinder import membrane_frame
        universe = host.universe
        dims = universe.dimensions.astype(float).copy()
        if abs(dims[3] - 90) > 1e-3 or abs(dims[4] - 90) > 1e-3:
            raise ConfigError("resize needs a box whose third vector is along z (alpha = beta = 90 degrees)")
        lengths, height = dims[:3], float(dims[2])
        frame = membrane_frame(universe, host.molecules)
        added = float(cfg.add_water)
        new_height = height + added
        positions = universe.atoms.positions.astype(float)
        heavy_mask = np.array([not str(n).upper().startswith(("H", "LP", "MW")) for n in universe.atoms.names])
        # Each molecule whole, its reference (first heavy) atom at z' in [0, height): the centred
        # frame, membrane centre at height / 2.
        reference = {}
        for mol in host.molecules:
            block = positions[mol.start:mol.stop]
            if mol.stop - mol.start > 1:
                block = unwrap(block, lengths)
            first = int(np.argmax(heavy_mask[mol.start:mol.stop]))
            z = np.mod(block[first, 2] + frame.shift, height)
            block[:, 2] += z - block[first, 2]
            positions[mol.start:mol.stop] = block
            reference[id(mol)] = (first, z)
        solid = [m for m in host.molecules if m.category not in SOLVENT_CATEGORIES]
        solid_z = np.concatenate([positions[m.start:m.stop, 2] for m in solid])
        top, bottom = float(solid_z.max()), float(solid_z.min())
        free = (height - top) + bottom - 2 * cfg.slab_margin
        if free < 6.0:
            raise ConfigError(f"Only {max(free + 2 * cfg.slab_margin, 0):.1f} A of bulk water lies between the "
                              f"periodic images of the protein and membrane (margin {cfg.slab_margin:g} A each "
                              "side); at least 6 A beyond the margins is needed to copy")
        copies = int(np.ceil(added / free))
        slab = added / copies
        # The gap opens in the middle of the bulk water band, between `top` and `bottom + height`.
        gap = np.mod((top + bottom + height) / 2, height)

        def placed(block, ref_index, target):
            out = block.copy()
            out[:, 2] += target - out[ref_index, 2]
            return out

        # w: distance above the gap (mod height). Originals go to w + added, shifted so the
        # membrane is centred in the taller box.
        moved = positions.copy()
        source = []
        for mol in host.molecules:
            first, z = reference[id(mol)]
            w = np.mod(z - gap, height)
            moved[mol.start:mol.stop] = placed(positions[mol.start:mol.stop], first,
                                               np.mod(w + added / 2, new_height))
            if mol.category in ("water", "ion"):
                u = w if w < height / 2 else w - height
                if -slab / 2 <= u < slab / 2:
                    source.append((mol, first, u))
        new_box = dims.copy()
        new_box[2] = new_height
        candidates = []
        for j in range(copies):
            for mol, first, u in source:
                target = np.mod(j * slab + slab / 2 + u - added / 2, new_height)
                candidates.append((mol.definition, mol.category, placed(positions[mol.start:mol.stop], first, target)))
        # Seams: copies within water_clash_distance of a kept heavy atom or an earlier copy are dropped.
        kept_heavy = moved[heavy_mask]
        accepted, accepted_heavy = [], np.zeros((0, 3))
        for item in candidates:
            definition, _, block = item
            heavy = block[[not a.atom_name.upper().startswith(("H", "LP", "MW")) for a in definition.atoms]]
            if _near(heavy, kept_heavy, cfg.water_clash_distance, new_box) or \
                    _near(heavy, accepted_heavy, cfg.water_clash_distance, new_box):
                continue
            accepted.append(item)
            accepted_heavy = np.concatenate([accepted_heavy, heavy])
        # Charge: drop copied ions of the excess sign.
        before = net_charge([(host, m) for m in host.molecules])
        excess = before + sum(d.charge for d, _, _ in accepted) - cfg.target_net_charge
        dropped = Counter()
        if abs(excess) > 0.01:
            for n in np.random.default_rng(cfg.seed).permutation(len(accepted)):
                definition, category, _ = accepted[n]
                if category == "ion" and definition.charge * excess > 0 and abs(definition.charge) <= abs(excess) + 0.01:
                    excess -= definition.charge
                    dropped[definition.name] += 1
                    accepted[n] = None
                if abs(excess) <= 0.01:
                    break
            accepted = [item for item in accepted if item is not None]
            if abs(excess) > 0.01:
                raise ConfigError(f"Cannot restore the net charge with the copied ions ({excess:+.3f} e left)")
        universe.atoms.positions = moved
        universe.dimensions = new_box
        new = []
        for name in dict.fromkeys(d.name for d, _, _ in accepted):
            block = place([item for item in accepted if item[0].name == name], stage, label=f"copied_{name}")
            new += [(block, m) for m in block.molecules]
        ordered = order(host, host.molecules, new=new)
        waters_added = sum(1 for _, c, _ in accepted if c == "water")
        area = abs(float(np.linalg.det(np.array([[dims[0], 0.0], [dims[1] * np.cos(np.radians(dims[5])),
                                                                   dims[1] * np.sin(np.radians(dims[5]))]]))))
        expected = round(added * area * WATER_DENSITY)
        report["resize"] = {
            "box_before_angstrom": [round(float(v), 3) for v in dims[:3]],
            "box_after_angstrom": [round(float(v), 3) for v in new_box[:3]],
            "added_angstrom": added, "copies": copies, "slab_thickness_angstrom": round(slab, 3),
            "bulk_water_between_images_angstrom": round(free + 2 * cfg.slab_margin, 3),
            "source_slab_molecules": dict(Counter(m.name for m, _, _ in source)),
            "added_molecules": dict(Counter(d.name for d, _, _ in accepted)),
            "removed_at_seams": len(candidates) - len(accepted) - sum(dropped.values()),
            "ions_dropped_for_charge": dict(dropped), "expected_waters_at_bulk_density": expected}
        report["salt"] = {"host": salt_summary([(host, m) for m in host.molecules]), "final": salt_summary(ordered)}
        report["charge"] = {"before": round(before, 6), "after": round(net_charge(ordered), 6),
                            "removed_ions": [], "added_ions": {}}
        print(f"[resize] box z {height:.2f} -> {new_height:.2f} A: {copies} cop{'y' if copies == 1 else 'ies'} of a "
              f"{slab:.1f} A bulk slab; {waters_added} waters (bulk density: {expected}), "
              f"{sum(1 for _, c, _ in accepted if c == 'ion')} ions added")
        report.update(write_system(stage, host, ordered, "resize", equilibration=cfg.equilibration,
                                   purpose="Equilibration of the taller box"))

    def lines(report):
        info = report.get("resize")
        if not info:
            return []
        before, after = info["box_before_angstrom"], info["box_after_angstrom"]
        out = ["", f"BOX: {before[0]:.2f} x {before[1]:.2f} x {before[2]:.2f} A -> "
                   f"{after[0]:.2f} x {after[1]:.2f} x {after[2]:.2f} A",
               f"  {info['copies']} cop{'y' if info['copies'] == 1 else 'ies'} of a {info['slab_thickness_angstrom']} A "
               f"slab of bulk water ({info['bulk_water_between_images_angstrom']} A of bulk water between the "
               "periodic images of the protein and membrane)",
               "  added: " + ", ".join(f"{k} {n}" for k, n in info["added_molecules"].items()),
               f"  copies removed at the seams: {info['removed_at_seams']}"]
        if info.get("expected_waters_at_bulk_density"):
            out.append(f"  waters expected at bulk density: {info['expected_waters_at_bulk_density']}; the NPT "
                       "equilibration settles the density at the seams")
        return out

    return run_edit("resize", cfg, "EXTEND THE WATER LAYER", build, "extend the water layer", lines,
                    dry_run=dry_run)
