"""Topology-driven CHARMM-GUI protein/ligand transplantation from two folders.

Unlike residue-name selections, topology molecule boundaries also support
multi-residue lipids, covalent cofactors and repeated molecule instances.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
import json
import re
import shutil
import tempfile
import warnings

import MDAnalysis as mda
from MDAnalysis.lib.distances import capped_distance
import numpy as np
import yaml

from charmm_environment import CATALOGUE, environment_name
from classify import AMINO_ACID_RESNAMES, is_hydrogen_name
from config import ConfigError, _UniqueKeySafeLoader, anchor_config_paths
from itp import parse_itp, ITPMoleculeType
from topology import parse_top_molecules, TopologyError, TopologyResult, audit_final_topology
from box_validation import validate_box
from kabsch import kabsch_fit
from ion_selection import select_counterions
from output import write_outputs, write_inspection_pdb


@dataclass
class CharmProtSpec:
    reference: str
    transplant: str
    reference_proteins: list[str] | None = None
    transplant_proteins: list[str] | None = None
    reference_ligands: list[str] | None = None
    transplant_ligands: list[str] | None = None
    ligands: str = "auto"
    reference_gro: str | None = None
    transplant_gro: str | None = None
    output_dir: str = "charmprot_output"
    clash_distance: float = 0.80
    # Keep every clashing lipid and sterol; only other environment molecules
    # (water, ions, detergents, solvents) are removed. Minimize afterwards.
    keep_lipids: bool = False
    neutralize: bool = True
    target_net_charge: float = 0.0
    ion_exclusion_distance: float = 10.0
    lipid_ion_exclusion_distance: float = 5.0
    environment_overrides: dict[str, str] = field(default_factory=dict)


def load_charmprot_config(path, check_paths=True, output_root=None):
    try:
        raw = yaml.load(Path(path).read_text(), Loader=_UniqueKeySafeLoader)
    except (OSError, yaml.YAMLError) as exc:
        raise ConfigError(str(exc)) from exc
    if not isinstance(raw, dict):
        raise ConfigError("charmprot YAML must be a mapping with reference and transplant folders")
    from config import _resolve_path_references
    raw = _resolve_path_references(raw)
    raw.pop("paths", None)
    minimization = raw.pop("minimization", {})
    if "charmprot" in raw:
        unknown_sections = set(raw) - {"charmprot"}
        if unknown_sections:
            raise ConfigError(f"Unknown charmprot sections: {sorted(unknown_sections)}")
        raw = raw["charmprot"]
        if not isinstance(raw, dict):
            raise ConfigError("charmprot must be a mapping")
    unknown = set(raw) - set(CharmProtSpec.__dataclass_fields__)
    if unknown:
        raise ConfigError(f"Unknown charmprot options: {sorted(unknown)}")
    for key in ("reference", "transplant", "output_dir"):
        if key not in raw and key == "output_dir":
            continue
        if not isinstance(raw.get(key), str) or not raw[key].strip():
            raise ConfigError(f"{key} must be a nonempty folder path")
    for key in ("reference_gro", "transplant_gro"):
        value = raw.get(key)
        if value is not None and (not isinstance(value, str) or not value.strip()):
            raise ConfigError(f"{key} must be a path to a GRO file")
    # Input folders and coordinate overrides are explicit paths; only the
    # generated output folder may be a bare name, resolved against --output.
    raw.setdefault("output_dir", CharmProtSpec.__dataclass_fields__["output_dir"].default)
    anchored = anchor_config_paths({"charmprot": raw, "minimization": minimization},
                                   config_path=path, output_root=output_root)
    raw, minimization = anchored["charmprot"], anchored["minimization"]
    cfg = CharmProtSpec(**raw)
    for key in ("reference_proteins", "transplant_proteins", "reference_ligands", "transplant_ligands"):
        value = getattr(cfg, key)
        if value is not None and (not isinstance(value, list) or
                any(not isinstance(x, str) or not x.strip() for x in value) or
                len(set(value)) != len(value) or (key.endswith("proteins") and not value)):
            raise ConfigError(f"{key} must be a list of unique topology molecule names")
    if cfg.ligands not in ("auto", "ignore"):
        raise ConfigError("ligands must be auto or ignore")
    # ignore applies to the transplant only: its ligands are left out, while the
    # reference's proteins and ligands are still removed as with auto.
    if cfg.ligands == "ignore" and cfg.transplant_ligands is not None:
        raise ConfigError("ligands: ignore leaves out every transplant ligand, so it cannot be "
                          "combined with transplant_ligands")
    for key in ("clash_distance", "ion_exclusion_distance", "lipid_ion_exclusion_distance", "target_net_charge"):
        value = getattr(cfg, key)
        if isinstance(value, bool) or not isinstance(value, (float, int)) or not np.isfinite(value):
            raise ConfigError(f"{key} must be a finite number")
        if key != "target_net_charge" and value <= 0:
            raise ConfigError(f"{key} must be positive (Angstrom)")
    for key in ("neutralize", "keep_lipids"):
        if type(getattr(cfg, key)) is not bool:
            raise ConfigError(f"{key} must be boolean")
    if not isinstance(cfg.environment_overrides, dict) or any(
        not isinstance(k, str) or v not in ("water", "ion", "lipid", "sterol", "detergent", "solvent")
        for k, v in cfg.environment_overrides.items()
    ):
        raise ConfigError("environment_overrides must map molecule names to environment categories")
    from config import MinimizationSpec, _SECTION_KEYS
    if not isinstance(minimization, dict) or set(minimization) - _SECTION_KEYS['minimization']:
        raise ConfigError('minimization must contain supported minimization settings')
    cfg.minimization = MinimizationSpec(**minimization)
    cfg.minimization.coordinates_path = str(Path(cfg.output_dir) / 'step5_input.gro')
    cfg.minimization.topology_path = str(Path(cfg.output_dir) / 'topol.top')
    cfg.minimization.include_dir = cfg.output_dir
    cfg.minimization.output_dir = str(Path(
        minimization.get('output_dir') or (Path(cfg.output_dir) / 'openmm_minimization')).resolve())
    cfg.minimization.validate(check_paths=False, require_inputs=False)
    output = Path(cfg.output_dir)
    for key in ("reference", "transplant"):
        folder = Path(getattr(cfg, key))
        if output == folder or folder in output.parents or output in folder.parents:
            raise ConfigError("output_dir must be separate from both input folders")
        if check_paths and not folder.is_dir():
            raise ConfigError(f"{key} folder does not exist: {folder}")
        override = getattr(cfg, key + '_gro')
        if override:
            coordinate = Path(override).resolve()
            if coordinate == output or output in coordinate.parents:
                raise ConfigError(f"{key}_gro must not be inside output_dir")
            if check_paths and not coordinate.is_file():
                raise ConfigError(f"{key}_gro does not exist: {coordinate}")
        if check_paths:
            if not (folder / 'topol.top').is_file() or not (folder / 'toppar').is_dir():
                raise ConfigError(f"{key} must contain topol.top and toppar/: {folder}")
            if not override and not (folder / 'step5_input.gro').is_file() and len(list(folder.glob('step5*.gro'))) != 1:
                raise ConfigError(f"{key} needs step5_input.gro or one step5*.gro; set its GRO override")
    return cfg


def _included_files(top, boundary):
    """Read the actual include graph, with no guessed inactive molecule branches."""
    visited, active, result = set(), set(), []
    def visit(path):
        path = path.resolve()
        if path in active:
            raise ConfigError(f"Cyclic topology include: {path}")
        if path in visited:
            return
        if not path.is_relative_to(boundary):
            raise ConfigError(f"Topology include outside system folder: {path}")
        if not path.is_file():
            raise ConfigError(f"Missing topology include: {path}")
        active.add(path)
        depth = 0
        for raw in path.read_text().splitlines():
            line = raw.split(";", 1)[0].strip()
            if re.match(r"#\s*(if|ifdef|ifndef)\b", line):
                depth += 1
            elif re.match(r"#\s*endif\b", line):
                depth -= 1
            match = re.match(r'#\s*include\s+["<]([^">]+)[">]', line)
            if match:
                if depth:
                    raise ConfigError(f"Conditional include in {path}; supply a preprocessed topology")
                visit(path.parent / match[1])
        active.remove(path)
        visited.add(path)
        result.append(path)
    visit(top)
    return result


@dataclass
class Molecule:
    name: str
    occurrence: int
    start: int
    stop: int
    category: str
    selected: bool
    definition: ITPMoleculeType


@dataclass
class BuildSystem:
    folder: Path
    gro: Path
    universe: object
    definitions: dict
    files: list
    molecules: list[Molecule]
    report: dict


def inspect_system(folder, proteins=None, ligands=None, ligand_policy="auto", gro=None, overrides=None):
    folder = Path(folder).resolve()
    if not (folder / "toppar").is_dir():
        raise ConfigError(f"Missing toppar directory in {folder}")
    if gro:
        coordinate = (folder / gro).resolve()
    elif (folder / "step5_input.gro").is_file():
        coordinate = folder / "step5_input.gro"
    else:
        candidates = sorted(folder.glob("step5*.gro"))
        if len(candidates) != 1:
            raise ConfigError(f"Expected one step5*.gro in {folder}; found {len(candidates)}. Set the GRO override.")
        coordinate = candidates[0]
    files = _included_files(folder / "topol.top", folder)
    definitions = {}
    for path in files:
        for name, definition in parse_itp(str(path)).items():
            if name in definitions:
                raise ConfigError(f"Duplicate active moleculetype {name} in {folder}")
            definitions[name] = definition
    rows = parse_top_molecules(str(folder / "topol.top"))
    names = {name for name, count in rows if count > 0}
    for label, selection in (("proteins", proteins), ("ligands", ligands)):
        if selection is not None and set(selection) - names:
            raise ConfigError(f"Unknown {label} in {folder}: {sorted(set(selection) - names)}")
    if set(proteins or ()) & set(ligands or ()):
        raise ConfigError("A molecule cannot be selected as both protein and ligand")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        universe = mda.Universe(str(coordinate))
    validate_box(universe.dimensions, universe.atoms.positions, label=str(coordinate))
    categories, evidence = {}, {}
    for name in names:
        if name not in definitions:
            raise ConfigError(f"No included ITP defines {name} in {folder}")
        definition = definitions[name]
        if proteins is not None and name in proteins:
            category, why = "protein", {"method": "explicit_yaml"}
        elif ligands is not None and name in ligands:
            category, why = "ligand", {"method": "explicit_yaml"}
        elif name in (overrides or {}):
            category, why = overrides[name], {"method": "environment_override"}
        elif re.fullmatch(r"PRO[A-Z0-9_]+", name.upper()) or (
                len({a.resnr for a in definition.atoms}) > 1 and
                sum(a.atom_name == "CA" and a.resname in AMINO_ACID_RESNAMES for a in definition.atoms) >= 2):
            category, why = "protein", {"method": "protein_name_or_backbone"}
        else:
            known = environment_name(name)
            if known is None:
                matches = [environment_name(rn) for rn in definition.resnames]
                if matches and all(matches) and len({m[0] for m in matches}) == 1:
                    known = matches[0]
            if known:
                category, why = known
            elif definition.atom_count == 1 and abs(definition.charge) > 0.01:
                category, why = "ion", {"method": "single_charged_atom"}
            else:
                category, why = "ligand", {"method": "unmatched_nonprotein", "inferred": True}
        categories[name], evidence[name] = category, why
    protein_names = proteins if proteins is not None else [n for n, _ in rows if categories.get(n) == "protein"]
    if not protein_names:
        raise ConfigError(f"No proteins detected in {folder}; set protein names explicitly")
    ligand_names = ([] if ligand_policy == "ignore" else ligands if ligands is not None else
                    [n for n, _ in rows if categories.get(n) == "ligand"])
    chosen = set(protein_names) | set(ligand_names)
    molecules, offset, counts = [], 0, Counter()
    entries, seen_environment = [], False
    seen_protein = False
    for name, count in rows:
        if count < 0:
            raise ConfigError(f"Negative molecule count for {name}")
        if count == 0:
            continue
        definition = definitions[name]
        category = categories[name]
        entries.append({"moleculetype": name, "count": count, "resnames": sorted(definition.resnames),
                        "category": category, "environment": category not in ("protein", "ligand"),
                        "selected": name in chosen, "evidence": evidence[name],
                        "after_protein_before_environment": seen_protein and not seen_environment,
                        "charge_per_molecule": definition.charge,
                        "atoms_per_molecule": definition.atom_count})
        seen_protein |= category == "protein"
        seen_environment |= category not in ("protein", "ligand")
        expected = [(a.resname[:5], a.atom_name[:5]) for a in definition.atoms]
        for _ in range(count):
            stop = offset + definition.atom_count
            actual = universe.atoms[offset:stop]
            observed = list(zip(actual.resnames, actual.names))
            if observed != expected:
                raise ConfigError(f"GRO/ITP atom-order mismatch in {folder}: {name}, starting at atom {offset + 1}")
            if offset and universe.atoms[offset - 1].resindex == universe.atoms[offset].resindex:
                raise ConfigError(f"GRO merges residues across molecule boundary at atom {offset + 1}")
            counts[name] += 1
            molecules.append(Molecule(name, counts[name], offset, stop, category, name in chosen, definition))
            offset = stop
    if offset != len(universe.atoms):
        raise ConfigError(f"Topology describes {offset} atoms but {coordinate} contains {len(universe.atoms)}")
    return BuildSystem(folder, coordinate, universe, definitions, files, molecules,
                       {"folder": str(folder), "coordinates": str(coordinate), "molecules": entries,
                        "selected_proteins": list(dict.fromkeys(protein_names)),
                        "protein_order_explicit": proteins is not None,
                        "selected_ligands": list(dict.fromkeys(ligand_names)),
                        "ligand_policy": ligand_policy,
                        "unknown_names_inferred_as_ligands": sorted(n for n in names if evidence[n].get("inferred"))})


def _ca_atoms(system, molecule):
    atoms = system.universe.atoms[molecule.start:molecule.stop]
    return atoms[np.array([a.atom_name == "CA" for a in molecule.definition.atoms])]


def _sequence_pairs(left, right):
    """Global residue alignment, +2 match, -1 substitution, -2 gap."""
    n, m = len(left), len(right)
    if not n or not m:
        raise ConfigError("Each selected protein chain needs CA atoms for alignment")
    if n * m > 25_000_000:
        raise ConfigError("Protein chain is too large for automatic sequence alignment; split into chain molecule types")
    score = np.empty((n + 1, m + 1), dtype=np.int32)
    trace = np.zeros((n + 1, m + 1), dtype=np.uint8)
    score[:, 0] = -2 * np.arange(n + 1)
    score[0, :] = -2 * np.arange(m + 1)
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            values = (score[i-1, j-1] + (2 if left[i-1] == right[j-1] else -1),
                      score[i-1, j] - 2, score[i, j-1] - 2)
            move = max(range(3), key=values.__getitem__)
            score[i, j], trace[i, j] = values[move], move
    pairs, i, j = [], n, m
    while i or j:
        move = trace[i, j] if i and j else 1 if i else 2
        if move == 0:
            pairs.append((i-1, j-1))
            i, j = i-1, j-1
        elif move == 1:
            i -= 1
        else:
            j -= 1
    return pairs[::-1]


def _protein_order(system):
    if not system.report["protein_order_explicit"]:
        return [mol for mol in system.molecules if mol.selected and mol.category == "protein"]
    return [mol for name in system.report["selected_proteins"] for mol in system.molecules
            if mol.name == name and mol.selected]


def align_systems(reference, donor):
    refs, mobiles = _protein_order(reference), _protein_order(donor)
    if len(refs) != len(mobiles):
        raise ConfigError("Different selected protein chain counts; specify corresponding protein lists")
    ref_indices, mobile_indices, chains = [], [], []
    aliases = {"HSD": "HIS", "HSE": "HIS", "HSP": "HIS", "HID": "HIS", "HIE": "HIS",
               "HIP": "HIS", "CYX": "CYS", "ASH": "ASP", "GLH": "GLU", "LYN": "LYS"}
    for ref, mobile in zip(refs, mobiles):
        ra, ma = _ca_atoms(reference, ref), _ca_atoms(donor, mobile)
        for system, atoms, mol in ((reference, ra, ref), (donor, ma, mobile)):
            from MDAnalysis.lib.distances import minimize_vectors
            steps = np.diff(atoms.positions.astype(float), axis=0)
            nearest = minimize_vectors(steps, system.universe.dimensions)
            if np.any(np.linalg.norm(steps - nearest, axis=1) > 0.1):
                raise ConfigError(f"Protein {mol.name} in {system.folder} crosses periodic boundaries; "
                                  "make the protein whole before alignment")
        rs = [aliases.get(str(x), str(x)) for x in ra.resnames]
        ms = [aliases.get(str(x), str(x)) for x in ma.resnames]
        pairs = _sequence_pairs(rs, ms)
        identity = sum(rs[i] == ms[j] for i, j in pairs) / max(len(pairs), 1)
        if identity < 0.3:
            raise ConfigError(f"Low sequence identity for {ref.name}/{mobile.name}: {identity:.1%}; check chain pairing")
        ref_indices.extend(int(ra.indices[i]) for i, _ in pairs)
        mobile_indices.extend(int(ma.indices[j]) for _, j in pairs)
        chains.append({"reference": ref.name, "transplant": mobile.name, "matched_CA": len(pairs),
                       "sequence_identity": identity, "reference_CA": len(ra), "transplant_CA": len(ma)})
    R, mobile_center, ref_center, rmsd = kabsch_fit(
        donor.universe.atoms[mobile_indices].positions.astype(float),
        reference.universe.atoms[ref_indices].positions.astype(float))
    donor.universe.atoms.positions = (donor.universe.atoms.positions - mobile_center) @ R.T + ref_center
    return {"method": "global_sequence_CA_kabsch", "chains": chains, "rmsd_angstrom": rmsd,
            "fit_atoms": len(ref_indices), "rotation": R.tolist(),
            "mobile_center": mobile_center.tolist(), "reference_center": ref_center.tolist()}


def _indices(molecules):
    return np.concatenate([np.arange(m.start, m.stop) for m in molecules]) if molecules else np.array([], dtype=int)


def _heavy(atoms):
    return atoms[np.array([not is_hydrogen_name(str(n)) and not str(n).upper().startswith(("LP", "MW"))
                           for n in atoms.names], dtype=bool)]


def remove_clashes(reference, incoming, cfg):
    remaining = [m for m in reference.molecules if not m.selected]
    env = _heavy(reference.universe.atoms[_indices(remaining)])
    block = _heavy(incoming)
    # 0.02 A margin covers GRO's 0.01 A coordinate rounding at the cutoff.
    pairs, distances = capped_distance(env.positions, block.positions, cfg.clash_distance + 0.02,
                                       box=reference.universe.dimensions, return_distances=True)
    atom_to_mol = np.full(len(reference.universe.atoms), -1, dtype=int)
    for i, mol in enumerate(remaining):
        atom_to_mol[mol.start:mol.stop] = i
    close = {}
    for pair, distance in zip(pairs, distances):
        i = int(atom_to_mol[env.indices[pair[0]]])
        close[i] = min(close.get(i, float("inf")), float(distance))
    removed, kept, retained = [], [], []
    for i, mol in enumerate(remaining):
        if i not in close:
            kept.append(mol)
        elif cfg.keep_lipids and mol.category in ("lipid", "sterol"):
            kept.append(mol)
            retained.append({"moleculetype": mol.name, "occurrence": mol.occurrence,
                             "category": mol.category, "minimum_distance_angstrom": close[i]})
        elif mol.category in ("protein", "ligand"):
            # Explicitly omitted/ignored biological molecules must not vanish
            # as a side effect of environment clash cleanup.
            raise ConfigError(f"Retained {mol.category} {mol.name}:{mol.occurrence} clashes with inserted block "
                              f"at {close[i]:.3f} A; include it in the replacement or resolve the overlap")
        else:
            removed.append({"moleculetype": mol.name, "occurrence": mol.occurrence,
                            "category": mol.category, "atoms": mol.stop-mol.start,
                            "minimum_distance_angstrom": close[i], "charge": mol.definition.charge})
    if retained:
        closest = min(entry["minimum_distance_angstrom"] for entry in retained)
        print(f"WARNING: keep_lipids retained {len(retained)} clashing lipid/sterol molecule(s) "
              f"(closest {closest:.2f} A); minimize before simulation.")
    return kept, removed, retained


def neutralize_system(reference, donor_molecules, kept, incoming, cfg):
    net = sum(m.definition.charge for m in donor_molecules + kept)
    report = {"before": net, "target": cfg.target_net_charge, "enabled": cfg.neutralize, "removed_ions": []}
    if not cfg.neutralize or abs(net - cfg.target_net_charge) <= 0.01:
        report["after"] = net
        return kept, report
    protected = _heavy(incoming)
    bio = [m for m in kept if m.category in ("protein", "ligand")]
    lipid = [m for m in kept if m.category in ("lipid", "sterol", "detergent")]
    groups = [(protected.positions, cfg.ion_exclusion_distance)]
    for molecules, distance in ((bio, cfg.ion_exclusion_distance), (lipid, cfg.lipid_ion_exclusion_distance)):
        if molecules:
            groups.append((_heavy(reference.universe.atoms[_indices(molecules)]).positions, distance))
    candidates = []
    for i, mol in enumerate(kept):
        if mol.category != "ion" or mol.definition.charge * (net - cfg.target_net_charge) <= 0:
            continue
        atoms = reference.universe.atoms[mol.start:mol.stop]
        if all(not len(capped_distance(atoms.positions, positions, distance,
                                       box=reference.universe.dimensions, return_distances=False))
               for positions, distance in groups):
            candidates.append((i, mol.name, mol.definition.charge))
    removed, excess = select_counterions(candidates, net - cfg.target_net_charge, 0.01, 42)
    if abs(excess) > 0.01:
        raise ConfigError(f"Cannot reach charge {cfg.target_net_charge} by removing eligible counterions; "
                          f"net charge {net:.6f}, residual {excess:.6f}")
    ids = {item["resindex"] for item in removed}
    report["removed_ions"] = [{"moleculetype": kept[i].name, "occurrence": kept[i].occurrence,
                               "charge": kept[i].definition.charge} for i in sorted(ids)]
    report["after"] = net - sum(item["charge"] for item in removed)
    return [mol for i, mol in enumerate(kept) if i not in ids], report


def _parameter_lines(system):
    """Collect active parameter sections, excluding molecule definitions."""
    sections = {}
    allowed = {"defaults", "atomtypes", "bondtypes", "constrainttypes", "pairtypes", "angletypes",
               "dihedraltypes", "nonbond_params", "cmaptypes", "implicit_genborn_params"}
    for path in system.files:
        section, continuation = None, ""
        for raw in path.read_text().splitlines():
            line = raw.split(";", 1)[0].strip()
            if not line:
                continue
            match = re.match(r"\[\s*([^]]+?)\s*\]", line)
            if match:
                section = match[1].lower()
                continue
            if line.startswith("#"):
                if section in allowed and not line.startswith("#include"):
                    raise ConfigError(f"Conditional/macro parameter section in {path}; preprocess the topology first")
                continue
            if section in allowed:
                continuation += " " + line.rstrip("\\").strip()
                if not line.endswith("\\"):
                    sections.setdefault(section, []).append(tuple(continuation.split()))
                    continuation = ""
        if continuation:
            raise ConfigError(f"Incomplete parameter continuation in {path}")
    return sections


def write_topology(stage, reference, donor, ordered):
    """Merge compatible flat CHARMM parameters and copy exact molecule sources.

    ITP atom signatures alone do not prove force-field compatibility. Parameter
    keys with different values are rejected rather than silently overridden.
    """
    from decimal import Decimal, InvalidOperation
    def normalized(tokens):
        result = []
        for token in tokens:
            try:
                result.append(Decimal(token))
            except InvalidOperation:
                result.append(token)
        return tuple(result)
    key_lengths = {"defaults": 0, "atomtypes": 1, "bondtypes": 3, "constrainttypes": 3,
                   "pairtypes": 3, "angletypes": 4, "dihedraltypes": 5,
                   "nonbond_params": 3, "cmaptypes": 6, "implicit_genborn_params": 1}
    merged, seen = {}, {}
    for system in (donor, reference):
        for section, rows in _parameter_lines(system).items():
            for tokens in rows:
                width = key_lengths[section]
                key = tokens[:width]
                symmetric = {"bondtypes": 2, "constrainttypes": 2, "pairtypes": 2,
                             "nonbond_params": 2, "angletypes": 3}
                if section in symmetric:
                    ntypes = symmetric[section]
                    types = key[:ntypes]
                    canonical = min(types, types[::-1])
                    key = canonical + key[ntypes:]
                    tokens = canonical + tokens[ntypes:]
                if section == "dihedraltypes" and len(tokens) >= 8 and tokens[4] in ("1", "4", "9"):
                    key += (tokens[7],)  # distinct periodic multiplicities coexist
                identity = (section, key)
                if identity in seen:
                    if normalized(tokens) != normalized(seen[identity]):
                        raise ConfigError(f"Incompatible force-field parameters: [{section}] {' '.join(key)}")
                    continue
                seen[identity] = tokens
                merged.setdefault(section, []).append(tokens)
    if "defaults" not in merged or "atomtypes" not in merged:
        raise ConfigError("Both systems must supply explicit GROMACS defaults and atom types")
    toppar = stage / "toppar"
    toppar.mkdir()
    parameter_order = list(key_lengths)
    def format_row(section, row):
        if section != "cmaptypes" or len(row) <= 8:
            return " ".join(row)
        # GROMACS reads at most 4095 characters per line; a CMAP grid joined
        # onto one line exceeds that, so continue it ten values per line.
        head, grid = " ".join(row[:8]), row[8:]
        lines = [head] + [" ".join(grid[i:i + 10]) for i in range(0, len(grid), 10)]
        return "\\\n".join(lines)
    (toppar / "forcefield.itp").write_text("\n".join(
        f"[ {section} ]\n" + "\n".join(format_row(section, row) for row in merged[section]) + "\n"
        for section in parameter_order if section in merged))
    selected_definitions, molecules = {}, []
    for system, mol in ordered:
        definition = mol.definition
        if mol.name in selected_definitions and selected_definitions[mol.name].source_path != definition.source_path:
            old = selected_definitions[mol.name]
            if Path(old.source_path).read_bytes() != Path(definition.source_path).read_bytes():
                raise ConfigError(f"Retained and inserted molecules share incompatible type {mol.name}; rename one type")
        selected_definitions[mol.name] = definition
        if molecules and molecules[-1][0] == mol.name:
            molecules[-1] = (mol.name, molecules[-1][1] + 1)
        else:
            molecules.append((mol.name, 1))
    # Keep molecule-local restraint conditionals, stripping global parameter
    # sections already merged above. Reject nested includes here: flattening
    # them without the active defines could silently change a molecule.
    includes, copied = ["toppar/forcefield.itp"], set()
    for name, definition in selected_definitions.items():
        source = Path(definition.source_path)
        if source in copied:
            continue
        copied.add(source)
        text = source.read_text()
        if re.search(r'^\s*#\s*include\b', text, re.M):
            raise ConfigError(f"Molecule ITP has nested includes: {source}; preprocess first")
        match = re.search(r'^\s*\[\s*moleculetype\s*\]', text, re.M | re.I)
        if not match:
            raise ConfigError(f"Missing molecule definition in {source}")
        target = f"molecule_{len(copied):04d}.itp"
        body = text[match.start():]
        end = re.search(r'^\s*\[\s*(system|molecules)\s*\]', body, re.M | re.I)
        if end:
            body = body[:end.start()]
        # A preprocessor wrapper before the first molecule would otherwise
        # lose its opener when extracting the molecule body.
        prefix_depth = 0
        for line in text[:match.start()].splitlines():
            if re.match(r"\s*#\s*(if|ifdef|ifndef)\b", line):
                prefix_depth += 1
            elif re.match(r"\s*#\s*endif\b", line):
                prefix_depth -= 1
        if prefix_depth:
            raise ConfigError(f"Conditional molecule definition in {source}; preprocess first")
        (toppar / target).write_text(body)
        includes.append("toppar/" + target)
    (stage / "topol.top").write_text("\n".join(f'#include "{inc}"' for inc in includes) +
        "\n\n[ system ]\nCHARMM-GUI protein transplant\n\n[ molecules ]\n" +
        "\n".join(f"{name:<20} {count}" for name, count in molecules) + "\n")
    return TopologyResult(str(stage / "topol.top"), str(toppar), molecules, includes, [])


def _write_report(folder, report):
    (folder / "charmprot_report.json").write_text(json.dumps(report, indent=2) + "\n")
    lines = ["CHARMM-GUI protein transplant", f"Status: {report['status']}",
             "Environment = water, ions, lipids, sterols, detergents, solvents.",
             "Ligands = explicit selections or non-protein molecules unmatched by environment detection.",
             f"Catalogue: {len(CATALOGUE['names'])} names, retrieved {CATALOGUE['retrieved_at']}"]
    for side in ("reference", "transplant"):
        if side not in report:
            continue
        lines += ["", side.upper(), report[side]["folder"]]
        for entry in report[side]["molecules"]:
            lines.append(f"  {entry['moleculetype']} x{entry['count']}: {entry['category']}; "
                         f"selected={entry['selected']}; {entry['evidence']['method']}")
        lines.append("  Unmatched names inferred as ligands: " +
                     ", ".join(report[side]["unknown_names_inferred_as_ligands"]))
    if "alignment" in report:
        lines.append(f"\nAlignment CA RMSD: {report['alignment']['rmsd_angstrom']:.4f} A")
    for key in ("charge", "clash_removals", "retained_clashing_lipids", "final_molecules", "error"):
        if key in report:
            lines += ["", key + ":", json.dumps(report[key], indent=2)]
    (folder / "charmprot_report.txt").write_text("\n".join(lines) + "\n")


GENERATED = {"step5_input.pdb", "step5_input.gro", "topol.top", "toppar", "index.ndx",
             "aligned_inspection.pdb", "charmprot_report.json", "charmprot_report.txt",
             "visualization", "view.pml", "view.vmd", "view.pse"}


def _publish_directory(source, destination):
    """Transactionally replace owned artifacts while preserving GUI locks and user files."""
    destination.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".charmprot-backup-", dir=destination.parent) as tmp:
        backup = Path(tmp)
        moved, installed = [], []
        try:
            for name in GENERATED:
                old = destination / name
                if old.exists() or old.is_symlink():
                    old.rename(backup / name)
                    moved.append(name)
            for item in source.iterdir():
                item.rename(destination / item.name)
                installed.append(item.name)
        except Exception:
            for name in installed:
                item = destination / name
                if item.is_dir() and not item.is_symlink():
                    shutil.rmtree(item)
                else:
                    item.unlink()
            for name in moved:
                (backup / name).rename(destination / name)
            raise


def run_charmprot(config_path, dry_run=False, prepare=False, output_root=None):
    cfg = load_charmprot_config(config_path, check_paths=not dry_run, output_root=output_root)
    prepare = prepare or cfg.minimization.enabled
    if prepare:
        from minimization_bundle import validate_preparation, check_destination
        validate_preparation(cfg.minimization)
        if not dry_run:
            check_destination(cfg.minimization.output_dir)
    if dry_run:
        print(f"charmprot configuration validated; no inputs opened. Output: {cfg.output_dir}")
        return 0
    destination = Path(cfg.output_dir)
    if destination.exists() and not destination.is_dir():
        raise ConfigError(f"output_dir exists and is not a directory: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    report = {"mode": "charmprot", "status": "detecting", "generated_at": datetime.now(timezone.utc).isoformat(),
              "catalogue": {"names": len(CATALOGUE["names"]), "retrieved_at": CATALOGUE["retrieved_at"],
                            "sources": CATALOGUE["sources"]}}
    with tempfile.TemporaryDirectory(prefix=".charmprot-", dir=destination.parent) as temporary:
        stage = Path(temporary) / "result"
        stage.mkdir()
        try:
            systems = []
            for side in ("reference", "transplant"):
                policy = cfg.ligands if side == "transplant" else "auto"
                system = inspect_system(getattr(cfg, side), getattr(cfg, side + "_proteins"),
                    getattr(cfg, side + "_ligands"), policy, getattr(cfg, side + "_gro"), cfg.environment_overrides)
                report[side] = system.report
                systems.append(system)
                print(f"[{side}] proteins={system.report['selected_proteins']}, ligands={system.report['selected_ligands']}")
            reference, donor = systems
            report["alignment"] = align_systems(reference, donor)
            print(f"[alignment] {report['alignment']['fit_atoms']} CA pairs; RMSD {report['alignment']['rmsd_angstrom']:.4f} A")
            selected = _protein_order(donor) + [m for m in donor.molecules if m.selected and m.category != "protein"]
            incoming = donor.universe.atoms[_indices(selected)]
            write_inspection_pdb(reference.universe.atoms[_indices([m for m in reference.molecules if m.selected])],
                                 incoming, str(stage / "aligned_inspection.pdb"))
            print("[clashes] Checking overlaps with the retained environment")
            kept, report["clash_removals"], report["retained_clashing_lipids"] = remove_clashes(
                reference, incoming, cfg)
            print("[charge] Checking charge and adjusting counterions")
            kept, report["charge"] = neutralize_system(reference, selected, kept, incoming, cfg)
            from collections import Counter
            print(f"Removed by class: {dict(Counter(m['category'] for m in report['clash_removals']))}")
            # All protein molecules form the prefix required by chain-ID output.
            ordered = ([(donor, m) for m in selected if m.category == "protein"] +
                       [(reference, m) for m in kept if m.category == "protein"] +
                       [(donor, m) for m in selected if m.category != "protein"] +
                       [(reference, m) for m in kept if m.category != "protein"])
            # One merge group per source block keeps the cost linear for water.
            groups = []
            current_system, blocks = None, []
            for system, mol in ordered:
                if current_system is not None and system is not current_system:
                    groups.append(current_system.universe.atoms[_indices(blocks)])
                    blocks = []
                current_system = system
                blocks.append(mol)
            if blocks:
                groups.append(current_system.universe.atoms[_indices(blocks)])
            final = mda.Merge(*groups)
            final.dimensions = reference.universe.dimensions.copy()
            print("[topology] Building and validating topology and coordinates")
            result = write_topology(stage, reference, donor, ordered)
            report["topology_audit"] = audit_final_topology(final, result)
            if abs(report["topology_audit"]["net_charge"] - report["charge"]["after"]) > 0.01:
                raise ConfigError("Final topology charge disagrees with molecule charge accounting")
            print(f"Audited final topology charge: {report['topology_audit']['net_charge']:+.6f} e")
            protein_counts = [m.stop-m.start for _, m in ordered if m.category == "protein"]
            def validate_written(_pdb, gro):
                loaded = mda.Universe(gro)
                return audit_final_topology(loaded, result)
            write_outputs(final, str(stage / "step5_input.pdb"), str(stage / "step5_input.gro"),
                          protein_chain_atom_counts=protein_counts, staged_validator=validate_written)
            # Index groups are based on molecule identity, including overridden ligands.
            index_groups = {name: [] for name in ("System", "Protein", "Ligand", "Water", "Ion", "Membrane", "Environment")}
            offset = 0
            for _, mol in ordered:
                ids = list(range(offset+1, offset+1+mol.stop-mol.start))
                index_groups["System"].extend(ids)
                category = {"protein": "Protein", "ligand": "Ligand", "water": "Water", "ion": "Ion",
                            "lipid": "Membrane", "sterol": "Membrane", "detergent": "Membrane"}.get(mol.category)
                if category:
                    index_groups[category].extend(ids)
                if mol.category not in ("protein", "ligand"):
                    index_groups["Environment"].extend(ids)
                offset += len(ids)
            (stage / "index.ndx").write_text("\n".join(
                f"[ {name} ]\n" + "\n".join(" ".join(map(str, ids[i:i+15])) for i in range(0, len(ids), 15))
                for name, ids in index_groups.items()))
            report["final_molecules"] = result.molecules_written
            report["status"] = "complete"
            # Do not expose temporary paths that disappear after publication.
            report["topology_audit"]["definition_files"] = {
                name: str(destination / Path(value).relative_to(stage))
                for name, value in report["topology_audit"]["definition_files"].items()}
            try:
                from visualization import comparison_groups, safe_scene
                replaced = reference.universe.atoms[_indices([m for m in reference.molecules if m.selected])]
                scene = comparison_groups(final, reference.universe, incoming, replaced, donor.universe)
                # Exact topology categories override residue-name heuristics for this mode.
                for role in list(scene):
                    if role.startswith(('inserted_', 'retained_')):
                        del scene[role]
                roles, offset = {}, 0
                for system, mol in ordered:
                    category = {'sterol': 'lipid', 'detergent': 'lipid', 'solvent': 'other'}.get(mol.category, mol.category)
                    role = ('inserted_' if system is donor else 'retained_') + category
                    roles.setdefault(role, []).extend(range(offset, offset + mol.stop - mol.start))
                    offset += mol.stop - mol.start
                scene.update({role: final.atoms[ids] for role, ids in roles.items()})
                safe_scene(stage, scene, 'charmprot')
            except Exception as exc:
                print(f'[VISUALIZATION WARNING] Scientific outputs are available; comparison views failed: {exc}')
            _write_report(stage, report)
            if (stage / 'view.pml').is_file():
                pml = stage / 'view.pml'
                pml.write_text(pml.read_text().replace(repr(str(stage)), repr(str(destination))))
            _publish_directory(stage, destination)
        except Exception as exc:
            report["status"], report["error"] = "failed", str(exc)
            failure = Path(temporary) / "failure"
            failure.mkdir()
            _write_report(failure, report)
            _publish_directory(failure, destination)
            raise
    print(f"[complete] Coordinates, topology, index and classification reports: {destination}")
    if prepare:
        from minimization_bundle import prepare_minimization
        prepare_minimization(cfg.minimization, assembly_mode="charmprot")
    from run_pipeline import _announce_output_location
    _announce_output_location(destination)
    return 0
