"""Shared machinery of the modes that edit one prepared system.

A host is a CHARMM-GUI-style folder (topol.top, toppar/, step5 GRO). A mode
keeps, removes, replaces and adds molecules; this module turns the resulting
ordered molecule list into a run-ready folder: coordinates, a merged topology,
CHARMM-GUI index groups, an optional GROMACS equilibration folder, viewer
scenes and a report. Output is staged in a temporary folder and published only
when everything succeeded.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import asdict, fields, is_dataclass
from datetime import datetime, timezone
from math import gcd
from pathlib import Path
import shutil
from typing import Any, Union, get_args, get_origin, get_type_hints
import json
import tempfile
import types
import warnings

import MDAnalysis as mda
from MDAnalysis.lib.distances import capped_distance
import numpy as np
import yaml

from charmprot import Molecule, inspect_system, write_topology, _heavy, _indices, _publish_directory
from config import ConfigError, _UniqueKeySafeLoader, anchor_config_paths
from equilibration import EquilibrationSpec, parse_equilibration, write_equilibration
from ion_selection import select_counterions
from itp import parse_itp
from output import write_outputs
from topology import audit_final_topology, TopologyError

WATER_MOLARITY = 55.51          # mol/L of pure water: concentration = units / waters * 55.51
SOLVENT_CATEGORIES = ("water", "ion", "solvent")
MEMBRANE_CATEGORIES = ("lipid", "sterol", "detergent")
CATEGORY_ORDER = ("protein", "ligand", "lipid", "sterol", "detergent", "solvent", "ion", "water")
ENVIRONMENT_CATEGORIES = ("water", "ion", "lipid", "sterol", "detergent", "solvent")
COMMON_GENERATED = {"step5_input.pdb", "step5_input.gro", "topol.top", "toppar", "index.ndx", "equilibration",
                    "repack", "host_protocol",
                    "visualization", "view.pml", "view.vmd", "view.pse"}


def generated(mode):
    """Files a mode owns in its output folder (replaced on a rerun)."""
    return COMMON_GENERATED | {f"{mode}_report.json", f"{mode}_report.txt"}


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

def _alternatives(annotation):
    return get_args(annotation) if get_origin(annotation) in (Union, types.UnionType) else (annotation,)


def check_value(value, annotation, where):
    """Validate one YAML value against a dataclass field annotation; return the coerced value."""
    options = _alternatives(annotation)
    if value is None:
        if type(None) in options:
            return None
        raise ConfigError(f"{where} must be set")
    annotation = next((t for t in options if t is not type(None)), Any)
    origin = get_origin(annotation)
    if annotation is Any:
        return value
    if annotation is bool:
        if type(value) is not bool:
            raise ConfigError(f"{where} must be true or false")
        return value
    if annotation is int:
        if type(value) is not int:
            raise ConfigError(f"{where} must be a whole number")
        return value
    if annotation is float:
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not np.isfinite(value):
            raise ConfigError(f"{where} must be a number")
        return float(value)
    if annotation is str:
        if not isinstance(value, str):
            raise ConfigError(f"{where} must be text")
        return value
    if origin is list or annotation is list:
        if not isinstance(value, list):
            raise ConfigError(f"{where} must be a list")
        item = (get_args(annotation) or (Any,))[0]
        return [check_value(v, item, f"{where}[{i}]") for i, v in enumerate(value)]
    if origin is dict or annotation is dict:
        if not isinstance(value, dict):
            raise ConfigError(f"{where} must be a mapping")
        return value
    if is_dataclass(annotation):
        if isinstance(value, annotation):
            value = asdict(value)
        if not isinstance(value, dict):
            raise ConfigError(f"{where} must be a mapping")
        return parse_dataclass(annotation, value, where)
    return value


def parse_dataclass(cls, raw, where):
    unknown = set(raw) - {f.name for f in fields(cls)}
    if unknown:
        raise ConfigError(f"Unknown {where} options: {sorted(unknown)}")
    hints = get_type_hints(cls)
    values = {}
    for f in fields(cls):
        if f.name in raw:
            values[f.name] = check_value(raw[f.name], hints[f.name], f"{where}.{f.name}")
    return cls(**values)


def load_section(path, section, spec_class, check_paths=True, output_root=None, optional_paths=("host_gro",)):
    """Read and validate the `section:` mapping of a system-edit YAML file.

    Every field is checked against `spec_class`'s annotations; `equilibration` is
    parsed with the shared schedule rules. Paths are anchored like every other mode:
    inputs must be absolute unless the file is a bundled example.
    """
    try:
        raw = yaml.load(Path(path).read_text(), Loader=_UniqueKeySafeLoader)
    except (OSError, yaml.YAMLError) as exc:
        raise ConfigError(str(exc)) from exc
    if not isinstance(raw, dict) or not isinstance(raw.get(section), dict):
        raise ConfigError(f"{section} YAML must contain a {section}: mapping")
    from config import _resolve_path_references
    raw = _resolve_path_references(raw)
    raw.pop("paths", None)
    if set(raw) - {section}:
        raise ConfigError(f"Unknown {section} sections: {sorted(set(raw) - {section})}")
    raw = dict(raw[section])
    unknown = set(raw) - {f.name for f in fields(spec_class)}
    if unknown:
        raise ConfigError(f"Unknown {section} options: {sorted(unknown)}")
    if not isinstance(raw.get("host"), str) or not raw["host"].strip():
        raise ConfigError(f"{section}.host must be a nonempty folder path")
    hints = get_type_hints(spec_class)
    for key, value in list(raw.items()):
        if value == "" and type(None) in _alternatives(hints[key]):
            raw[key] = None  # an empty optional field in the GUI means "not set"
    for key in optional_paths:
        value = raw.get(key)
        if value is not None and (not isinstance(value, str) or not value.strip()):
            raise ConfigError(f"{section}.{key} must be a file path")
    raw.setdefault("output_dir", next(f.default for f in fields(spec_class) if f.name == "output_dir"))
    schedules = [f for f in fields(spec_class) if hints[f.name] is EquilibrationSpec]
    given = {f.name: raw.pop(f.name, None) for f in schedules}
    raw = anchor_config_paths({section: raw}, config_path=path, output_root=output_root)[section]
    cfg = parse_dataclass(spec_class, raw, section)
    for f in schedules:
        default = f.default_factory()
        setattr(cfg, f.name, parse_equilibration(given[f.name], section, enabled=default.enabled, name=f.name,
                                                 default=default))
    overrides = getattr(cfg, "environment_overrides", {}) or {}
    if any(not isinstance(k, str) or v not in ENVIRONMENT_CATEGORIES for k, v in overrides.items()):
        raise ConfigError(f"{section}.environment_overrides must map molecule names to environment categories")
    output, host = Path(cfg.output_dir), Path(cfg.host)
    if output == host or host in output.parents or output in host.parents:
        raise ConfigError(f"{section}.output_dir must be separate from the host folder")
    if check_paths:
        if not host.is_dir():
            raise ConfigError(f"host folder does not exist: {host}")
        if not (host / "topol.top").is_file() or not (host / "toppar").is_dir():
            raise ConfigError(f"host must contain topol.top and toppar/: {host}")
        for key in optional_paths:
            value = getattr(cfg, key, None)
            for item in (value if isinstance(value, list) else [value]):
                if item and not Path(item).is_file():
                    raise ConfigError(f"{section}.{key} does not exist: {item}")
    return cfg


def positive(cfg, section, *keys, allow_zero=False):
    for key in keys:
        value = getattr(cfg, key)
        if value < 0 or (value == 0 and not allow_zero):
            raise ConfigError(f"{section}.{key} must be {'0 or more' if allow_zero else 'greater than 0'}")


# ---------------------------------------------------------------------------
# Molecules created by a mode
# ---------------------------------------------------------------------------

class ParameterSource:
    """Extra force-field files whose parameter sections join the merged forcefield.itp."""

    def __init__(self, files):
        self.files = [Path(f) for f in files]


class Placed:
    """New molecules: a universe and its Molecule entries, like a BuildSystem."""

    def __init__(self, universe, molecules):
        self.universe, self.molecules = universe, molecules


def place(items, workdir, label="placed"):
    """Molecules from (definition, category, positions) items, as one Placed block.

    The block is written to GRO and read back, so it merges with the GRO-read host.
    """
    items = list(items)
    if not items:
        return Placed(None, [])
    resindex, residues, names, previous = [], [], [], None
    for number, (definition, _, _) in enumerate(items):
        for atom in definition.atoms:
            key = (number, atom.resnr)
            if key != previous:
                residues.append((atom.resnr, atom.resname))
                previous = key
            resindex.append(len(residues) - 1)
            names.append(atom.atom_name)
    universe = mda.Universe.empty(len(resindex), n_residues=len(residues), atom_resindex=resindex,
                                  trajectory=True)
    universe.add_TopologyAttr("names", names)
    universe.add_TopologyAttr("resnames", [r[1] for r in residues])
    universe.add_TopologyAttr("resids", [r[0] for r in residues])
    universe.atoms.positions = np.concatenate([np.asarray(p, dtype=float).reshape(-1, 3) for _, _, p in items])
    universe.dimensions = [100.0, 100.0, 100.0, 90.0, 90.0, 90.0]
    path = Path(workdir) / f".{label}_{id(items)}.gro"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        universe.atoms.write(str(path))
        reread = mda.Universe(str(path))
    path.unlink()
    if len(reread.atoms) != len(resindex):
        raise ConfigError(f"Writing the new {label} molecules lost atoms")
    molecules, offset, counts = [], 0, Counter()
    for definition, category, _ in items:
        counts[definition.name] += 1
        molecules.append(Molecule(definition.name, counts[definition.name], offset,
                                  offset + definition.atom_count, category, True, definition))
        offset += definition.atom_count
    return Placed(reread, molecules)


def positions_of(system, mol):
    return system.universe.atoms[mol.start:mol.stop].positions.astype(float)


def unwrap(positions, lengths):
    """Make a molecule whole by chaining each atom to the previous one (orthorhombic images)."""
    from addbinder import _unwrap_sequential
    return _unwrap_sequential(np.asarray(positions, dtype=float), np.asarray(lengths, dtype=float))


def order(host, kept, new=(), replacements=None):
    """The final molecule order: host order of the kept molecules, new molecules inserted.

    kept: host molecules to keep; replacements: {id(host molecule): (system, molecule)}
    put in that molecule's place; new: (system, molecule) pairs, each placed right after the
    last molecule of its moleculetype, else after the last molecule of an earlier category
    (protein, ligand, lipid, ..., ion, water), else at the end.
    """
    replacements = replacements or {}
    keep = {id(m) for m in kept}
    ordered = []
    for mol in host.molecules:
        if id(mol) in replacements:
            ordered.append(replacements[id(mol)])
        elif id(mol) in keep:
            ordered.append((host, mol))
    for system, mol in new:
        rank = CATEGORY_ORDER.index(mol.category) if mol.category in CATEGORY_ORDER else len(CATEGORY_ORDER)
        position = None
        for i in range(len(ordered) - 1, -1, -1):
            if ordered[i][1].name == mol.name:
                position = i + 1
                break
        if position is None:
            for i in range(len(ordered) - 1, -1, -1):
                other = ordered[i][1].category
                if (CATEGORY_ORDER.index(other) if other in CATEGORY_ORDER else len(CATEGORY_ORDER)) <= rank:
                    position = i + 1
                    break
        ordered.insert(len(ordered) if position is None else position, (system, mol))
    return ordered


def merge(ordered, dimensions):
    groups, current, block = [], None, []
    for system, mol in ordered:
        if current is not None and system is not current:
            groups.append(current.universe.atoms[_indices(block)])
            block = []
        current = system
        block.append(mol)
    if block:
        groups.append(current.universe.atoms[_indices(block)])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        final = mda.Merge(*groups)
    final.dimensions = np.array(dimensions, dtype=float).copy()
    return final


def net_charge(ordered):
    return float(sum(mol.definition.charge for _, mol in ordered))


def write_index(path, ordered, extra=None):
    """CHARMM-GUI index groups (SOLU, MEMB, SOLV ...) plus readable ones; extra: name -> {id(mol)}."""
    names = ("System", "Protein", "Ligand", "Water", "Ion", "Membrane", "Environment",
             "SOLU", "MEMB", "SOLV", "SOLU_MEMB", "SYSTEM")
    extra = extra or {}
    groups = {name: [] for name in (*names, *extra)}
    charmm_gui = {"protein": ("SOLU", "SOLU_MEMB"), "ligand": ("SOLU", "SOLU_MEMB"),
                  "lipid": ("MEMB", "SOLU_MEMB"), "sterol": ("MEMB", "SOLU_MEMB"),
                  "detergent": ("MEMB", "SOLU_MEMB"), "water": ("SOLV",), "ion": ("SOLV",), "solvent": ("SOLV",)}
    labels = {"protein": "Protein", "ligand": "Ligand", "water": "Water", "ion": "Ion", "lipid": "Membrane",
              "sterol": "Membrane", "detergent": "Membrane"}
    offset = 0
    for _, mol in ordered:
        ids = list(range(offset + 1, offset + 1 + mol.stop - mol.start))
        groups["System"].extend(ids)
        groups["SYSTEM"].extend(ids)
        for name in charmm_gui.get(mol.category, ("SOLU", "SOLU_MEMB")):
            groups[name].extend(ids)
        if mol.category in labels:
            groups[labels[mol.category]].extend(ids)
        if mol.category not in ("protein", "ligand"):
            groups["Environment"].extend(ids)
        for name, members in extra.items():
            if id(mol) in members:
                groups[name].extend(ids)
        offset += len(ids)
    Path(path).write_text("\n".join(
        f"[ {name} ]\n" + "\n".join(" ".join(map(str, ids[i:i + 15])) for i in range(0, len(ids), 15))
        for name, ids in groups.items() if ids) + "\n")
    return {name: ids for name, ids in groups.items() if ids}


# ---------------------------------------------------------------------------
# Water and ions
# ---------------------------------------------------------------------------

def heavy_positions(system, molecules):
    if not molecules:
        return np.zeros((0, 3))
    return _heavy(system.universe.atoms[_indices(molecules)]).positions.astype(float)


def _near(points, others, cutoff, box):
    """Indices of `points` within cutoff of any of `others`."""
    if not len(points) or not len(others) or cutoff <= 0:
        return set()
    pairs = capped_distance(np.asarray(points, dtype=float), np.asarray(others, dtype=float), cutoff, box=box,
                            return_distances=False)
    return {int(a) for a, _ in pairs}


def molecules_near(system, molecules, others, cutoff, box):
    """Molecules with a heavy atom (any atom for ions) within cutoff of `others`."""
    if not molecules or not len(others):
        return set()
    owner, points = [], []
    for i, mol in enumerate(molecules):
        atoms = system.universe.atoms[mol.start:mol.stop]
        heavy = _heavy(atoms) if len(atoms) > 1 else atoms
        points.append(heavy.positions.astype(float))
        owner.extend([i] * len(heavy))
    owner = np.array(owner)
    return {int(owner[i]) for i in _near(np.concatenate(points), others, cutoff, box)}


def salt_units(cation_charge, anion_charge):
    """(cations, anions) per neutral formula unit: (1, 1) NaCl, (1, 2) CaCl2."""
    zc, za = int(round(cation_charge)), int(round(-anion_charge))
    if zc <= 0 or za <= 0 or abs(cation_charge - zc) > 1e-3 or abs(-anion_charge - za) > 1e-3:
        raise ConfigError("Salt ions need whole charges: a positive cation and a negative anion")
    g = gcd(zc, za)
    return za // g, zc // g


def concentration(waters, cations, anions, units):
    """Salt concentration (mol/L) from whole formula units: min(cations/a, anions/b) / waters * 55.51."""
    if not waters:
        return None
    count = min(cations // units[0], anions // units[1])
    return round(count / waters * WATER_MOLARITY, 4)


def water_oxygens(system, waters):
    return np.array([_heavy(system.universe.atoms[m.start:m.stop]).positions[0] for m in waters], dtype=float) \
        if waters else np.zeros((0, 3))


def replace_waters(system, solvent, species, avoid, box, exclusion=5.0, spacing=5.0, seed=42, occupied=None,
                   allowed=None):
    """Turn bulk waters into ions.

    solvent: the kept solvent molecules of `system`; species: [(ion definition, count)];
    avoid: [(positions, distance)] to stay away from. Waters are drawn at random (seeded), at
    least `spacing` from every ion. allowed: optional predicate on a water oxygen position.
    Returns (kept solvent, [(definition, position)], {name: shortfall}).
    """
    waters = [(i, m) for i, m in enumerate(solvent) if m.category == "water"]
    oxygens = water_oxygens(system, [m for _, m in waters])
    blocked = set()
    for positions, distance in avoid:
        blocked |= _near(oxygens, positions, distance, box)
    candidates = [n for n in range(len(waters)) if n not in blocked
                  and (allowed is None or allowed(oxygens[n]))]
    ions = [m for m in solvent if m.category == "ion"]
    taken = [p for p in (occupied if occupied is not None else
                         [system.universe.atoms[m.start].position.astype(float) for m in ions])]
    rng = np.random.default_rng(seed)
    order_ = rng.permutation(len(candidates))
    # Interleave the species so a shortfall hits them evenly.
    queue = []
    remaining = {d.name: c for d, c in species}
    while any(remaining.values()):
        for definition, _ in species:
            if remaining[definition.name]:
                queue.append(definition)
                remaining[definition.name] -= 1
    added, chosen = [], set()
    cursor = 0
    for n in order_:
        if cursor == len(queue):
            break
        oxygen = oxygens[candidates[n]]
        if taken and _near(oxygen[None, :], np.array(taken), spacing, box):
            continue
        added.append((queue[cursor], oxygen))
        chosen.add(waters[candidates[n]][0])
        taken.append(oxygen)
        cursor += 1
    shortfall = Counter(d.name for d in queue[cursor:])
    return [m for i, m in enumerate(solvent) if i not in chosen], added, dict(shortfall)


def balance_charge(system, solvent, net, target, cation, anion, avoid, box, spacing=5.0, seed=42,
                   tolerance=0.01, allow_removal=True, occupied=None):
    """Bring the net charge to `target`: remove counterions far from `avoid`, then add ions.

    cation/anion: ion definitions to add (monovalent or not; the remainder must be a multiple
    of their charge). Returns (kept solvent, added [(definition, position)], report).
    """
    report = {"before": round(net, 6), "target": target, "removed_ions": [], "added_ions": {}}
    excess = net - target
    if abs(excess) <= tolerance:
        report["after"] = round(net, 6)
        return solvent, [], report
    if allow_removal:
        blocked = set()
        ions = [(i, m) for i, m in enumerate(solvent) if m.category == "ion" and m.definition.charge * excess > 0]
        points = np.array([system.universe.atoms[m.start].position for _, m in ions], dtype=float) \
            if ions else np.zeros((0, 3))
        for positions, distance in avoid:
            blocked |= _near(points, positions, distance, box)
        candidates = [(i, m.name, m.definition.charge) for n, (i, m) in enumerate(ions) if n not in blocked]
        # Remove what can be removed exactly; whatever is left is added below.
        removable = sum(c for _, _, c in candidates)
        goal = excess if abs(removable) >= abs(excess) else removable
        removed, _ = select_counterions(candidates, goal, tolerance, seed) if abs(goal) > tolerance else ([], 0)
        ids = {item["resindex"] for item in removed}
        report["removed_ions"] = [{"moleculetype": solvent[i].name, "occurrence": solvent[i].occurrence,
                                   "charge": solvent[i].definition.charge} for i in sorted(ids)]
        solvent = [m for i, m in enumerate(solvent) if i not in ids]
        excess -= sum(item["charge"] for item in removed)
    added = []
    if abs(excess) > tolerance:
        ion = anion if excess > 0 else cation
        count = abs(excess) / abs(ion.charge)
        if abs(count - round(count)) > 0.02:
            raise ConfigError(f"Cannot neutralize a remaining charge of {excess:+.3f} e with {ion.name} "
                              f"({ion.charge:+g} e); choose a monovalent counterion")
        if occupied is not None:
            occupied = [p for p in occupied] + [system.universe.atoms[m.start].position.astype(float)
                                                for m in solvent if m.category == "ion"]
        solvent, added, shortfall = replace_waters(system, solvent, [(ion, int(round(count)))], avoid, box,
                                                   spacing=spacing, seed=seed + 1, occupied=occupied)
        if shortfall:
            raise ConfigError(f"Found room for only {len(added)} of {int(round(count))} {ion.name} counterions in "
                              "the bulk water")
        excess += sum(d.charge for d, _ in added)
        report["added_ions"] = dict(Counter(d.name for d, _ in added))
    report["after"] = round(target + excess, 6)
    return solvent, added, report


def ion_block(added, workdir):
    """New ions as Placed blocks per moleculetype: [(Placed, Molecule)]."""
    new = []
    for name in dict.fromkeys(d.name for d, _ in added):
        items = [(d, "ion", p[None, :]) for d, p in added if d.name == name]
        block = place(items, workdir, label=f"ions_{name}")
        new += [(block, m) for m in block.molecules]
    return new


def salt_summary(ordered):
    """Waters, ions per type and the monovalent-pair concentration of a molecule list."""
    waters = sum(m.category == "water" for _, m in ordered)
    ions = Counter(m.name for _, m in ordered if m.category == "ion")
    charge = {m.name: m.definition.charge for _, m in ordered if m.category == "ion"}
    cations = sum(n for name, n in ions.items() if charge[name] > 0)
    anions = sum(n for name, n in ions.items() if charge[name] < 0)
    return {"waters": waters, "ions": dict(ions),
            "concentration_molar": round(min(cations, anions) / waters * WATER_MOLARITY, 4) if waters else None}


def monovalent_ions(host):
    """The host's most common monovalent cation and anion definitions (or None)."""
    ions = Counter(m.name for m in host.molecules if m.category == "ion")
    definitions = host.definitions
    cations = [n for n in ions if abs(definitions[n].charge - 1) < 1e-3]
    anions = [n for n in ions if abs(definitions[n].charge + 1) < 1e-3]
    pick = (lambda names: definitions[max(names, key=lambda n: ions[n])] if names else None)
    return pick(cations), pick(anions)


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

CLASH_DISTANCE = 0.2   # A; atoms of different molecules closer than this are reported (minimization fails)


def find_clashes(positions, sizes, box, cutoff=CLASH_DISTANCE, limit=20):
    """Atom pairs of different molecules closer than cutoff (hydrogens included).
    `sizes` are the molecules' atom counts in order. Returns (count, closest examples)."""
    owner = np.repeat(np.arange(len(sizes)), sizes)
    pairs, distances = capped_distance(positions, positions, cutoff, box=box)
    across = owner[pairs[:, 0]] < owner[pairs[:, 1]]
    pairs, distances = pairs[across], distances[across]
    order_ = np.argsort(distances)[:limit]
    return int(len(pairs)), [(int(pairs[i, 0]), int(pairs[i, 1]), float(distances[i])) for i in order_]


def write_system(stage, host, ordered, mode, parameter_files=(), equilibration=None, purpose="",
                 extra_index=None, scene=None, extra_mdp="", restrained=None, more_equilibrations=()):
    """Write the run-ready folder into `stage`; returns a report fragment.

    parameter_files: extra force-field files (atom types, bonded types) for new molecules.
    scene: (roles {role: [id(mol)]}, palette, alignment text) for the PyMOL/VMD views.
    restrained: definitions whose POSRES coverage the equilibration report checks
    (default: every protein and ligand).
    """
    final = merge(ordered, host.universe.dimensions)
    donor = ParameterSource(list(host.files) + [Path(p) for p in parameter_files])
    topology = write_topology(stage, host, donor, ordered)
    top = Path(topology.top_path)
    top.write_text(top.read_text().replace("CHARMM-GUI protein transplant", f"GMXtransplant {mode}"))
    audit = audit_final_topology(final, topology)
    expected = net_charge(ordered)
    if abs(audit["net_charge"] - expected) > 0.01:
        raise ConfigError("Final topology charge disagrees with molecule charge accounting")
    protein_counts = [m.stop - m.start for _, m in ordered if m.category == "protein"]

    def validate_written(_pdb, gro):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return audit_final_topology(mda.Universe(gro), topology)

    write_outputs(final, str(stage / "step5_input.pdb"), str(stage / "step5_input.gro"),
                  protein_chain_atom_counts=protein_counts, staged_validator=validate_written)
    index_groups = write_index(stage / "index.ndx", ordered, extra_index)
    result = {"final_molecules": topology.molecules_written, "final_net_charge": round(audit["net_charge"], 6),
              "atoms": len(final.atoms), "box_angstrom": [round(float(v), 3) for v in final.dimensions[:3]]}
    # The finished system, every atom: anything this close will not minimize.
    count, examples = find_clashes(final.atoms.positions, [m.stop - m.start for _, m in ordered], final.dimensions)
    result["final_clash_check"] = {"distance_angstrom": CLASH_DISTANCE, "pairs": count, "closest": [
        f"{final.atoms[a].resname}{final.atoms[a].resid}:{final.atoms[a].name} (atom {a + 1}) - "
        f"{final.atoms[b].resname}{final.atoms[b].resid}:{final.atoms[b].name} (atom {b + 1}) {d:.2f} A"
        for a, b, d in examples]}
    if count:
        print(f"[clash check] WARNING: {count} atom pair(s) of different molecules closer than {CLASH_DISTANCE:g} A; "
              f"closest: {result['final_clash_check']['closest'][0]}. Minimization may fail; see the report.")
    else:
        print(f"[clash check] no atoms of different molecules closer than {CLASH_DISTANCE:g} A")
    if restrained is None:
        restrained = {m.name: m.definition for _, m in ordered if m.category in ("protein", "ligand")}.values()
    # The main equilibration, then any further sets (key, spec, purpose), e.g. the membrane repack.
    for key, spec, text in [("equilibration", equilibration, purpose), *more_equilibrations]:
        if spec is None or not spec.enabled:
            continue
        result[key] = write_equilibration(spec, stage, topology, list(restrained), index_groups,
                                          purpose=text, extra_mdp=extra_mdp)
        print(f"[{key}] {spec.output_dir}/: {result[key]['md_length_ps'] / 1000:g} ns "
              f"in {len(spec.stages)} stages; run sample_run.sh")
        for warning in result[key]["warnings"]:
            print(f"      WARNING: {warning}")
    if scene is not None:
        _write_scene(stage, final, ordered, mode, *scene)
    return result


def write_host_protocol(stage, host, folder="host_protocol"):
    """The host's own run inputs on the new system: every mdp file of the host folder
    (mdout.mdp excepted), its shell scripts and README, with the new system's topology,
    index and coordinates under the GRO name the host used (e.g. system.gro), so the
    host's protocol and scripts run unchanged. A -D define for a macro that no file of
    the new topology uses (e.g. the restraint of a removed molecule) is dropped, since
    grompp refuses unused macros. Returns a report, or None when the host has no mdp files.
    """
    import re
    source = Path(host.folder)
    mdps = [p for p in sorted(source.glob("*.mdp")) if p.name != "mdout.mdp"]
    if not mdps:
        return None
    bundle = stage / folder
    bundle.mkdir()
    shutil.copytree(stage / "toppar", bundle / "toppar")
    for name in ("topol.top", "index.ndx"):
        shutil.copy2(stage / name, bundle / name)
    gro_name = Path(host.gro).name
    shutil.copy2(stage / "step5_input.gro", bundle / gro_name)
    topology_text = (bundle / "topol.top").read_text() + "".join(
        p.read_text(errors="replace") for p in sorted((bundle / "toppar").rglob("*")) if p.is_file())
    dropped = {}
    for path in mdps:
        lines = path.read_text().splitlines()
        for number, line in enumerate(lines):
            match = re.match(r"(\s*define\s*=\s*)([^;]*)(.*)", line, re.I)
            if not match:
                continue
            keep, gone = [], []
            for token in match[2].split():
                name = token[2:].split("=", 1)[0] if token.startswith("-D") else None
                (gone if name and not re.search(rf"\b{re.escape(name)}\b", topology_text) else keep).append(token)
            if gone:
                lines[number] = f"{match[1]}{' '.join(keep)}{match[3]}"
                dropped[path.name] = gone
        header = []
        if path.name in dropped:
            header = [f"; GMXtransplant: dropped {' '.join(dropped[path.name])} (not used by the new topology)"]
        (bundle / path.name).write_text("\n".join(header + lines) + "\n")
    extras = [p for p in sorted(source.iterdir())
              if p.is_file() and (p.suffix == ".sh" or p.name.upper().startswith("README"))]
    for path in extras:
        shutil.copy2(path, bundle / path.name)
    return {"directory": folder, "coordinates": gro_name, "mdp_files": [p.name for p in mdps],
            "other_files": [p.name for p in extras], "dropped_defines": dropped}


def _write_scene(stage, final, ordered, mode, roles, palette, alignment):
    from visualization import safe_scene
    role_of = {member: role for role, members in roles.items() for member in members}
    groups, offset = {}, 0
    for _, mol in ordered:
        role = role_of.get(id(mol)) or ("retained_" + {"solvent": "other", "detergent": "lipid",
                                                        "sterol": "lipid"}.get(mol.category, mol.category))
        groups.setdefault(role, []).extend(range(offset, offset + mol.stop - mol.start))
        offset += mol.stop - mol.start
    safe_scene(stage, {role: final.atoms[ids] for role, ids in groups.items()}, mode, palette=palette,
               alignment=alignment)


def write_report(folder, mode, report, heading, lines):
    """<mode>_report.json and a readable <mode>_report.txt (heading, status, mode lines, common facts)."""
    folder = Path(folder)
    (folder / f"{mode}_report.json").write_text(json.dumps(report, indent=2, default=str) + "\n")
    text = [f"GMXtransplant: {heading}", f"Status: {report['status']}"]
    if "error" in report:
        text += ["", f"Error: {report['error']}"]
    text += list(lines)
    if "charge" in report:
        charge = report["charge"]
        text += ["", f"CHARGE: {charge['before']:+.3f} e -> {charge['after']:+.3f} e; "
                     f"{len(charge.get('removed_ions', []))} counterion(s) removed"
                 + (", added " + ", ".join(f"{n} {k}" for k, n in charge["added_ions"].items())
                    if charge.get("added_ions") else "")]
    if "salt" in report:
        text += [f"SALT: host {report['salt']['host']['concentration_molar']} M, "
                 f"final {report['salt']['final']['concentration_molar']} M"]
    if "final_molecules" in report:
        text += ["", "FINAL MOLECULES"] + [f"  {name:<10} {count}" for name, count in report["final_molecules"]]
    clashes = report.get("clashes")
    if clashes and clashes.get("pairs"):
        text += ["", f"CLASHES: {clashes['pairs']} atom pair(s) of different molecules closer than "
                     f"{clashes['distance_angstrom']:g} A (hydrogens included); "
                 + ("removed " + ", ".join(f"{k} {n}" for k, n in clashes["removed"].items())
                    if clashes.get("removed") else "kept (remove_clashes: false)")]
        text += [f"  {e}" for e in clashes["examples"][:10]]
    final_check = report.get("final_clash_check")
    if final_check is not None:
        text += ["", f"FINAL CLASH CHECK: {final_check['pairs']} atom pair(s) of different molecules closer "
                     f"than {final_check['distance_angstrom']:g} A"
                 + ("" if final_check["pairs"] else " (none: ready for minimization)")]
        text += [f"  {e}" for e in final_check["closest"][:10]]
    for key, title in (("repack", "REPACK"), ("equilibration", "EQUILIBRATION")):
        eq = report.get(key)
        if eq:
            text += ["", f"{title}: {eq['directory']}/ ({eq['md_length_ps'] / 1000:g} ns of MD, "
                         f"{len(eq['stages'])} stages); run sample_run.sh"]
            text += [f"  WARNING: {w}" for w in eq["warnings"]]
    protocol = report.get("host_protocol")
    if protocol:
        text += ["", f"HOST PROTOCOL: {protocol['directory']}/ (the host's {len(protocol['mdp_files'])} mdp file(s) "
                     f"and {', '.join(protocol['other_files']) or 'no scripts'} on {protocol['coordinates']})"]
        text += [f"  {name}: dropped {' '.join(gone)} (not in the new topology)"
                 for name, gone in protocol["dropped_defines"].items()]
    if report.get("warnings"):
        text += ["", "WARNINGS"] + [f"  {w}" for w in report["warnings"]]
    (folder / f"{mode}_report.txt").write_text("\n".join(text) + "\n")


def run_edit(mode, cfg, banner, build, heading, report_lines, dry_run=False, dry_run_note=""):
    """Stage, build and publish one system-edit run.

    build(stage, host, report) fills `report` and writes the system into `stage`.
    report_lines(report) gives the mode's lines of the text report.
    """
    if dry_run:
        print(f"{mode} configuration validated; no inputs opened. {dry_run_note}output: {cfg.output_dir}")
        return 0
    destination = Path(cfg.output_dir)
    if destination.exists() and not destination.is_dir():
        raise ConfigError(f"output_dir exists and is not a directory: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    print("=" * 78)
    print(banner)
    print("=" * 78)
    settings = {k: v for k, v in asdict(cfg).items() if k != "equilibration"}
    report = {"mode": mode, "status": "running", "generated_at": datetime.now(timezone.utc).isoformat(),
              "host": cfg.host, "settings": settings, "warnings": []}
    owned = generated(mode)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        with tempfile.TemporaryDirectory(prefix=f".{mode}-", dir=destination.parent) as temporary:
            stage = Path(temporary) / "result"
            stage.mkdir()
            try:
                host = inspect_system(cfg.host, gro=getattr(cfg, "host_gro", None),
                                      overrides=getattr(cfg, "environment_overrides", None))
                print(f"[host] {host.gro}: proteins={host.report['selected_proteins']}, "
                      f"ligands={host.report['selected_ligands']}")
                build(stage, host, report)
                if (stage / "view.pml").is_file():
                    pml = stage / "view.pml"
                    pml.write_text(pml.read_text().replace(repr(str(stage)), repr(str(destination))))
                report["status"] = "complete"
                write_report(stage, mode, report, heading, report_lines(report))
                _publish_directory(stage, destination, owned)
            except (ConfigError, TopologyError) as exc:
                report["status"], report["error"] = "failed", str(exc)
                failure = Path(temporary) / "failure"
                failure.mkdir()
                write_report(failure, mode, report, heading, report_lines(report))
                _publish_directory(failure, destination, owned)
                raise ConfigError(str(exc)) from exc
    for warning in report["warnings"]:
        print(f"WARNING: {warning}")
    print(f"[complete] {destination}")
    from run_pipeline import _announce_output_location
    _announce_output_location(destination)
    return 0


def rewritten_definition(text, workdir, name):
    """Write an edited molecule ITP into the stage and parse it back."""
    path = Path(workdir) / f"{name}.itp"
    path.write_text(text)
    definitions = parse_itp(str(path))
    if name not in definitions:
        raise ConfigError(f"The rewritten ITP for {name} does not define {name}")
    return definitions[name]


__all__ = ["EquilibrationSpec", "Molecule", "Placed", "ParameterSource", "load_section", "place", "order",
           "merge", "write_system", "run_edit", "balance_charge", "replace_waters", "ion_block"]
