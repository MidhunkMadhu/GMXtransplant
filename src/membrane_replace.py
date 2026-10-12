"""Replace the membrane of a prepared system with any equilibrated bilayer.

The generalisation of the popc mode: the reference may hold several lipid types
and sterols (e.g. POPC/cholesterol 70:30, POPC/POPE/POPS, asymmetric leaflets).
Every host lipid and sterol is removed; the reference bilayer is repeated
periodically, cut to the host box with its midplane on the host's, and its
lipids replace them. The protein and ligands keep their coordinates (the OPM
placement). Reference lipids overlapping the protein or ligands, or a periodic
copy at the box edge, are removed whole; the reference hydration water replaces
host water near the new headgroups; counterions and salt are restored.

Each reference residue is matched to a moleculetype: the host's own (e.g.
POPC.itp, CHL1.itp) or one from lipid_itps; resname_map renames residues first
(e.g. CLOL -> CHL1). Atom names must match the ITP; a different atom order is
fixed by name.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
import warnings

import MDAnalysis as mda
import numpy as np

from charmm_environment import environment_name
from charmprot import _heavy, _indices, _included_files
from config import ConfigError
from equilibration import EquilibrationSpec
from itp import parse_itp
from popc_membrane import (_anchors, _cell, _heavy_mask, _midplane_shift, _plane, _reference_waters,
                           _resolve_contacts, _tail_end, _touching, _unwrap_sequential, place_water,
                           restore_salt, tile)
from system_edit import (MEMBRANE_CATEGORIES, SOLVENT_CATEGORIES, balance_charge, ion_block, load_section,
                         monovalent_ions, net_charge, order, place, positive, run_edit, salt_summary, write_system)
from topology import TopologyError, check_parameters

LIPID_CATEGORIES = ("lipid", "sterol")


@dataclass
class MembraneSpec:
    host: str
    reference_gro: str = ""
    lipid_itps: list[str] = field(default_factory=list)     # ITPs of reference lipids the host lacks
    lipid_forcefield: str | None = None   # parameters for lipid_itps the host forcefield.itp lacks
    resname_map: dict = field(default_factory=dict)          # reference residue name -> moleculetype
    host_gro: str | None = None
    membrane_center_z: float | None = None
    max_tilt: float = 5.0
    anchor_count: int = 3
    anchor_separation: float = 20.0
    protein_clash_distance: float = 3.0
    lateral_search: int = 16
    water_layer: float = 10.0
    protein_hydration_distance: float = 4.0
    reference_contact_distance: float = 2.0
    water_clash_distance: float = 2.5
    neutralize: bool = True
    target_net_charge: float = 0.0
    restore_salt: bool = True
    ion_exclusion_distance: float = 10.0
    lipid_ion_exclusion_distance: float = 5.0
    environment_overrides: dict = field(default_factory=dict)
    output_dir: str = "membrane_output"
    equilibration: EquilibrationSpec = field(default_factory=EquilibrationSpec)


def load_config(path, check_paths=True, output_root=None):
    cfg = load_section(path, "membrane", MembraneSpec, check_paths, output_root,
                       optional_paths=("host_gro", "reference_gro", "lipid_itps", "lipid_forcefield"))
    if not cfg.reference_gro:
        raise ConfigError("membrane.reference_gro must be the equilibrated bilayer GRO")
    positive(cfg, "membrane", "max_tilt", "anchor_separation", "protein_clash_distance", "water_clash_distance",
             "reference_contact_distance", "protein_hydration_distance", "ion_exclusion_distance",
             "lipid_ion_exclusion_distance")
    positive(cfg, "membrane", "water_layer", allow_zero=True)
    if cfg.max_tilt >= 45:
        raise ConfigError("membrane.max_tilt must be below 45 degrees")
    if cfg.anchor_count < 3 or not 1 <= cfg.lateral_search <= 400:
        raise ConfigError("membrane.anchor_count must be 3 or more and lateral_search 1 to 400")
    if any(not isinstance(k, str) or not isinstance(v, str) for k, v in cfg.resname_map.items()):
        raise ConfigError("membrane.resname_map must map residue names to moleculetype names")
    return cfg


# ---------------------------------------------------------------------------
# Host and reference
# ---------------------------------------------------------------------------

def _head_index(definition):
    """The headgroup atom: P for phospholipids, else the first oxygen (sterol O3), else atom 0."""
    names = [a.atom_name for a in definition.atoms]
    for candidate in ("P",):
        if candidate in names:
            return names.index(candidate)
    return next((i for i, a in enumerate(definition.atoms) if a.atom_name.upper().startswith("O")), 0)


@dataclass
class Frame:
    shift: float
    height: float
    midplane: float
    planes: dict
    counts: dict


def measure_host(host, cfg):
    """Midplane, leaflet planes and lipid counts of the host bilayer (any phospholipids)."""
    universe = host.universe
    cell, height = _cell(universe.dimensions)
    lengths = universe.dimensions[:3].astype(float)
    membrane = [m for m in host.molecules if m.category in MEMBRANE_CATEGORIES]
    heads, tails, others = [], [], []
    for mol in membrane:
        names = [a.atom_name for a in mol.definition.atoms]
        if "P" in names:
            whole = _unwrap_sequential(universe.atoms[mol.start:mol.stop].positions.astype(float), lengths)
            p = names.index("P")
            heads.append((mol, universe.atoms[mol.start + p].position.astype(float)))
            tails.append(_tail_end(whole, p) - whole[p, 2] + heads[-1][1][2])
        else:
            others.append(mol)
    if len(heads) < 6:
        raise ConfigError("Fewer than six phospholipids in the host; cannot measure the bilayer")
    positions = np.array([p for _, p in heads])
    shift = float(_midplane_shift(np.array(tails), height))
    zc = np.mod(positions[:, 2] + shift, height)
    leaflet = np.where(zc > height / 2, "upper", "lower")
    planes, counts = {}, {}
    for side in ("upper", "lower"):
        on = leaflet == side
        if on.sum() < cfg.anchor_count:
            raise ConfigError(f"Found {int(on.sum())} phosphorus atoms in the {side} leaflet; the membrane normal "
                              "must lie along the box z axis")
        _, tilt, spread = _plane(positions[on, :2], zc[on])
        if tilt > cfg.max_tilt:
            raise ConfigError(f"The {side} leaflet's P plane is tilted {tilt:.1f} degrees from the xy plane "
                              f"(max_tilt {cfg.max_tilt:g})")
        chosen, spacing = _anchors(positions[on, :2], cfg.anchor_count, cfg.anchor_separation, cell)
        if spacing + 1e-6 < cfg.anchor_separation:
            raise ConfigError(f"The {side} leaflet has no {cfg.anchor_count} phospholipids whose P atoms are "
                              f"{cfg.anchor_separation:g} A apart (best {spacing:.1f} A)")
        planes[side] = {"mean_z": round(float(np.mean(zc[on]) - shift), 3), "tilt_degrees": round(tilt, 3),
                        "spread_angstrom": round(spread, 3), "phosphorus_atoms": int(on.sum()),
                        "anchor_separation_angstrom": round(spacing, 3)}
        counts[side] = Counter(heads[i][0].name for i in np.flatnonzero(on))
    for mol in others:
        whole = _unwrap_sequential(universe.atoms[mol.start:mol.stop].positions.astype(float), lengths)
        centre = np.mod(whole[:, 2].mean() + shift, height)
        counts["upper" if centre > height / 2 else "lower"][mol.name] += 1
    thickness = planes["upper"]["mean_z"] - planes["lower"]["mean_z"]
    if not 20.0 <= thickness <= 70.0:
        raise ConfigError(f"Headgroup planes are {thickness:.1f} A apart; expected a flat bilayer normal to z")
    centre = (planes["upper"]["mean_z"] + planes["lower"]["mean_z"]) / 2 + shift
    if cfg.membrane_center_z is not None:
        centre = cfg.membrane_center_z + shift
    return Frame(shift, height, float(np.mod(centre - shift, height)), planes, counts)


@dataclass
class Reference:
    path: str
    lipids: list      # per lipid: (definition, whole positions, midplane at z = 0, ITP order)
    leaflets: list
    heads: np.ndarray
    waters: np.ndarray
    cell: np.ndarray
    height: float
    report: dict


def _definitions(host, cfg):
    """Moleculetypes the reference lipids may use: the host's lipids and sterols plus lipid_itps."""
    found = {name: d for name, d in host.definitions.items()}
    extra = {}
    for path in cfg.lipid_itps:
        for name, definition in parse_itp(path).items():
            if name in found and Path(found[name].source_path).read_bytes() != Path(path).read_bytes():
                raise ConfigError(f"lipid_itps defines {name}, which the host defines differently; rename one")
            extra[name] = definition
    return {**found, **extra}, extra


def load_reference(path, definitions, cfg):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        universe = mda.Universe(str(path))
    if universe.dimensions is None or np.any(np.abs(universe.dimensions[3:] - 90) > 1e-3):
        raise ConfigError("The reference bilayer needs a rectangular box")
    cell, height = _cell(universe.dimensions)
    lengths = universe.dimensions[:3].astype(float)
    lipids, skipped, reordered = [], Counter(), 0
    for residue in universe.residues:
        resname = str(residue.resname)
        name = cfg.resname_map.get(resname, resname)
        known = environment_name(name) or environment_name(resname)
        category = known[0] if known else ("lipid" if name in definitions and resname in cfg.resname_map else None)
        if category not in LIPID_CATEGORIES:
            skipped[resname] += 1      # water, ions and anything else: not part of the bilayer
            continue
        if name not in definitions:
            raise ConfigError(f"Reference lipid {resname} has no moleculetype {name}: add its ITP to lipid_itps "
                              "or map it with resname_map")
        definition = definitions[name]
        expected = [a.atom_name for a in definition.atoms]
        names = [str(n) for n in residue.atoms.names]
        order_ = slice(None)
        if names != expected:
            if sorted(names) != sorted(expected) or len(set(names)) != len(names):
                raise ConfigError(f"Reference {resname} {residue.resid} atoms do not match the {name} ITP")
            order_ = [names.index(n) for n in expected]
            reordered += 1
        lipids.append((definition, _unwrap_sequential(residue.atoms.positions.astype(float)[order_], lengths)))
    if len(lipids) < 20:
        raise ConfigError(f"The reference {path} holds {len(lipids)} lipids; expected a bilayer")
    phospholipids = [(d, x) for d, x in lipids if "P" in [a.atom_name for a in d.atoms]]
    if len(phospholipids) < 6:
        raise ConfigError("The reference has fewer than six phospholipids; cannot find its midplane")
    shift = _midplane_shift(np.array([_tail_end(x, _head_index(d)) for d, x in phospholipids]), height)
    inverse = np.linalg.inv(cell)
    out, leaflets = [], []
    for definition, whole in lipids:
        h = _head_index(definition)
        whole = whole.copy()
        whole[:, 2] += shift - height * np.floor((whole[h, 2] + shift) / height) - height / 2
        whole[:, :2] -= np.floor(whole[h, :2] @ inverse) @ cell
        out.append((definition, whole))
        leaflets.append("upper" if whole[h, 2] > 0 else "lower")
    heads = np.array([x[_head_index(d)] for d, x in out])
    planes = {}
    for side in ("upper", "lower"):
        on = [i for i, (d, _) in enumerate(out) if leaflets[i] == side and "P" in [a.atom_name for a in d.atoms]]
        _, tilt, spread = _plane(heads[on, :2], heads[on, 2])
        if tilt > 5.0:
            raise ConfigError(f"The reference {side} leaflet is tilted {tilt:.1f} degrees; expected a flat bilayer")
        composition = Counter(d.name for i, (d, _) in enumerate(out) if leaflets[i] == side)
        planes[side] = {"lipids": int(sum(composition.values())), "composition": dict(composition),
                        "mean_z": round(float(heads[on, 2].mean()), 3), "spread_angstrom": round(spread, 3)}
    waters = _reference_waters(universe, lengths, shift, cell, height)
    area = abs(np.linalg.det(cell))
    report = {"path": str(path), "lipids": len(out), "waters": len(waters),
              "box_angstrom": [round(float(v), 3) for v in lengths],
              "phosphate_thickness_angstrom": round(planes["upper"]["mean_z"] - planes["lower"]["mean_z"], 3),
              "area_per_lipid_angstrom2": {s: round(area / planes[s]["lipids"], 2) for s in planes},
              "leaflets": planes, "residues_reordered_to_itp": reordered, "ignored_residues": dict(skipped)}
    return Reference(str(path), out, leaflets, heads, waters, cell, height, report)


# ---------------------------------------------------------------------------
# Placement
# ---------------------------------------------------------------------------

def _lipids_at(cfg, reference, host_cell, frame, origin, protected_heavy, box):
    candidates = []
    for k, offset in tile(reference.heads[:, :2], reference.cell, host_cell, origin):
        definition, whole = reference.lipids[k]
        coords = whole.copy()
        coords[:, :2] += offset
        coords[:, 2] += frame.midplane
        candidates.append((k, reference.leaflets[k], coords, _heavy_mask(definition)))
    clash = _touching([c[2][c[3]] for c in candidates], protected_heavy, cfg.protein_clash_distance, box)
    keep = [c for i, c in enumerate(candidates) if i not in clash]
    edge = _resolve_contacts([c[2][c[3]] for c in keep], cfg.reference_contact_distance, box)
    kept = [c for i, c in enumerate(keep) if i not in edge]

    def by_type(items):
        return dict(Counter(f"{side} {reference.lipids[k][0].name}" for k, side, _, _ in items))
    report = {"tiled": by_type(candidates), "removed_protein_overlap": by_type([candidates[i] for i in clash]),
              "removed_box_edge_overlap": by_type([keep[i] for i in edge])}
    return kept, report


def place_lipids(cfg, host, frame, reference, protected):
    box = host.universe.dimensions
    host_cell, _ = _cell(box)
    protected_heavy = _heavy(host.universe.atoms[_indices(protected)]).positions
    side = int(np.ceil(np.sqrt(cfg.lateral_search)))
    origins = [(i / side) * reference.cell[0] + (j / side) * reference.cell[1]
               for i in range(side) for j in range(side)][:cfg.lateral_search]
    best = None
    for number, origin in enumerate(origins):
        kept, report = _lipids_at(cfg, reference, host_cell, frame, origin, protected_heavy, box)
        if best is None or len(kept) > len(best[0]):
            best = (kept, report, origin, number)
    kept, report, origin, number = best
    kept.sort(key=lambda item: (reference.lipids[item[0]][0].name, item[1] != "upper", item[0]))
    report["placed"] = dict(Counter(f"{s} {reference.lipids[k][0].name}" for k, s, _, _ in kept))
    report["lateral_search"] = {"registrations_tried": len(origins), "chosen": number,
                                "origin_angstrom": [round(float(v), 3) for v in origin]}
    return kept, origin, report


def run(config_path, dry_run=False, output_root=None):
    cfg = load_config(config_path, check_paths=not dry_run, output_root=output_root)

    def build(stage, host, report):
        work = stage.parent / "work"
        work.mkdir()
        frame = measure_host(host, cfg)
        removed = [m for m in host.molecules if m.category in MEMBRANE_CATEGORIES]
        report["host_membrane"] = {
            "midplane_z": round(frame.midplane, 3),
            "phosphate_thickness_angstrom": round(frame.planes["upper"]["mean_z"] - frame.planes["lower"]["mean_z"], 3),
            "planes": frame.planes, "lipids_removed": {s: dict(frame.counts[s]) for s in ("upper", "lower")}}
        print(f"[host] midplane z {frame.midplane:.2f} A; removing "
              + ", ".join(f"{n} {k}" for k, n in Counter(m.name for m in removed).items()))
        definitions, extra = _definitions(host, cfg)
        if extra:
            forcefields = [str(p) for p in host.files]
            if cfg.lipid_forcefield:
                forcefields += [str(p) for p in _included_files(Path(cfg.lipid_forcefield).resolve(),
                                                                Path(cfg.lipid_forcefield).resolve().parent)]
            try:
                check_parameters([d.source_path for d in extra.values()], forcefields,
                                 "Set lipid_forcefield to a file with the missing lipid parameters.")
            except TopologyError as exc:
                raise ConfigError(str(exc)) from exc
        reference = load_reference(cfg.reference_gro, definitions, cfg)
        report["reference"] = reference.report
        print(f"[reference] {reference.report['lipids']} lipids: upper "
              + ", ".join(f"{n} {k}" for k, n in reference.report["leaflets"]["upper"]["composition"].items())
              + "; lower " + ", ".join(f"{n} {k}" for k, n in reference.report["leaflets"]["lower"]["composition"].items())
              + f"; P-P {reference.report['phosphate_thickness_angstrom']:.1f} A (host "
              f"{report['host_membrane']['phosphate_thickness_angstrom']:.1f} A)")
        difference = reference.report["phosphate_thickness_angstrom"] - report["host_membrane"]["phosphate_thickness_angstrom"]
        if abs(difference) > 3.0:
            report["warnings"].append(f"The reference bilayer is {difference:+.1f} A thicker than the host's; the water "
                                      "layers adjust during the NPT equilibration")
        protected = [m for m in host.molecules if m.category not in MEMBRANE_CATEGORIES + SOLVENT_CATEGORIES]
        lipids, origin, report["placement"] = place_lipids(cfg, host, frame, reference, protected)
        print("[membrane] placed " + ", ".join(f"{k} {n}" for k, n in report["placement"]["placed"].items())
              + f" (best of {report['placement']['lateral_search']['registrations_tried']} lateral registrations)")
        # popc's water step: reference hydration water and the host solvent kept around it.
        waters, water_definition, kept_solvent, removed_solvent, report["water"] = place_water(
            cfg, host, frame, _AsPopcReference(reference), lipids, None, protected, origin)
        report["solvent_removed"] = {reason: dict(Counter(m.name for m in items))
                                     for reason, items in removed_solvent.items()}
        if report["water"].get("reference_waters"):
            print(f"[water] {report['water']['reference_waters']} reference waters hydrate the bilayer")
        elif report["water"].get("note"):
            report["warnings"].append(report["water"]["note"])
        box = host.universe.dimensions
        lipid_heavy = np.concatenate([c[2][c[3]] for c in lipids]) if lipids else np.zeros((0, 3))
        protected_heavy = _heavy(host.universe.atoms[_indices(protected)]).positions
        net = net_charge([(host, m) for m in protected + kept_solvent]) \
            + sum(reference.lipids[k][0].charge for k, _, _, _ in lipids)
        added = []
        report["charge"] = {"before": round(net, 6), "after": round(net, 6), "removed_ions": [], "added_ions": {}}
        if cfg.neutralize:
            cation, anion = monovalent_ions(host)
            if cation is None or anion is None:
                raise ConfigError("The host has no monovalent cation and anion to neutralize with")
            avoid = [(protected_heavy, cfg.ion_exclusion_distance), (lipid_heavy, cfg.lipid_ion_exclusion_distance)]
            kept_solvent, added, report["charge"] = balance_charge(host, kept_solvent, net, cfg.target_net_charge,
                                                                   cation, anion, avoid, box)
            print(f"[charge] {report['charge']['before']:+.3f} e -> {report['charge']['after']:+.3f} e")
        half = report["water"]["slab_half_thickness_angstrom"]
        kept_solvent, salt_ions, salt_report = restore_salt(cfg, host, kept_solvent, len(waters), frame, half,
                                                            protected)
        if salt_report.get("note"):
            report["warnings"].append(salt_report["note"])
        added += salt_ions
        # Assemble: new lipids by type at the host membrane's place, reference water after the host water.
        def category(definition):
            known = environment_name(definition.name)
            return known[0] if known and known[0] in LIPID_CATEGORIES else "lipid"
        items = [(reference.lipids[k][0], category(reference.lipids[k][0]), coords) for k, _, coords, _ in lipids]
        lipid_block = place(items, work, label="lipids")
        leaflet_of = {id(m): side for m, (_, side, _, _) in zip(lipid_block.molecules, lipids)}
        new = [(lipid_block, m) for m in lipid_block.molecules]
        if len(waters):
            water_block = place([(water_definition, "water", w) for w in waters], work, label="water")
            new += [(water_block, m) for m in water_block.molecules]
        new += ion_block(added, work)
        ordered = order(host, protected + kept_solvent, new=new)
        report["salt"] = {"host": salt_summary([(host, m) for m in host.molecules]), "final": salt_summary(ordered),
                          "restore": salt_report}
        parameter_files = []
        if cfg.lipid_forcefield:
            parameter_files = list(_included_files(Path(cfg.lipid_forcefield).resolve(),
                                                   Path(cfg.lipid_forcefield).resolve().parent))
        report.update(write_system(
            stage, host, ordered, "membrane", parameter_files=parameter_files, equilibration=cfg.equilibration,
            purpose="Equilibration of the new bilayer",
            extra_index={"Upper_leaflet": {k for k, v in leaflet_of.items() if v == "upper"},
                         "Lower_leaflet": {k for k, v in leaflet_of.items() if v == "lower"}},
            scene=({"inserted_lipids": set(leaflet_of)},
                   {"inserted_lipids": ("Equilibrated bilayer", [.95, .68, .12], True)},
                   "Protein and ligands keep their input coordinates; the bilayer is new."),
            restrained=[m.definition for m in protected]))

    def lines(report):
        out = []
        host = report.get("host_membrane")
        if host:
            out += ["", "HOST MEMBRANE (protein and ligands are not moved: OPM placement kept)",
                    f"  midplane z {host['midplane_z']:.2f} A; P-P thickness {host['phosphate_thickness_angstrom']:.2f} A",
                    "  removed upper: " + ", ".join(f"{k} {n}" for k, n in host["lipids_removed"]["upper"].items()),
                    "  removed lower: " + ", ".join(f"{k} {n}" for k, n in host["lipids_removed"]["lower"].items())]
        ref = report.get("reference")
        if ref:
            out += ["", "REFERENCE BILAYER", f"  {ref['path']}",
                    f"  {ref['lipids']} lipids, box {ref['box_angstrom'][0]:.2f} x {ref['box_angstrom'][1]:.2f} A, "
                    f"P-P thickness {ref['phosphate_thickness_angstrom']:.2f} A, area per lipid "
                    f"{ref['area_per_lipid_angstrom2']['upper']:.1f} / {ref['area_per_lipid_angstrom2']['lower']:.1f} A^2"]
            for side in ("upper", "lower"):
                out.append(f"  {side}: " + ", ".join(f"{k} {n}" for k, n in ref["leaflets"][side]["composition"].items()))
        placement = report.get("placement")
        if placement:
            out += ["", "PLACED: " + ", ".join(f"{k} {n}" for k, n in placement["placed"].items()),
                    "  removed for protein overlap: " + (", ".join(
                        f"{k} {n}" for k, n in placement["removed_protein_overlap"].items()) or "none"),
                    "  removed at the box edge: " + (", ".join(
                        f"{k} {n}" for k, n in placement["removed_box_edge_overlap"].items()) or "none")]
        water = report.get("water")
        if water and water.get("reference_waters"):
            out += [f"WATER: {water['reference_waters']} reference waters hydrate the new bilayer"]
        return out

    return run_edit("membrane", cfg, "REPLACE THE MEMBRANE WITH AN EQUILIBRATED BILAYER", build,
                    "replace the membrane with an equilibrated bilayer", lines, dry_run=dry_run)


class _AsPopcReference:
    """The parts of popc's ReferenceBilayer that its water step reads."""

    def __init__(self, reference):
        self.waters, self.cell = reference.waters, reference.cell
        self.report = {"leaflets": {s: {"mean_z": reference.report["leaflets"][s]["mean_z"]}
                                    for s in ("upper", "lower")}}

