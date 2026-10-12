"""Shared run of the protonate and mutate modes: residue changes in protein ITPs.

Every copy of an edited protein moleculetype gets the same new ITP; each copy's
coordinates are rebuilt. Missing CHARMM36 parameters are added from the bundled
library, water and ions overlapping new atoms are removed (lipids only when asked),
and counterions restore the net charge.
"""
from __future__ import annotations

from collections import Counter
from pathlib import Path

import numpy as np

from charmprot import _heavy, _indices
from config import ConfigError
from protein_edit import ProteinItp, residue_state
from system_edit import (MEMBRANE_CATEGORIES, SOLVENT_CATEGORIES, balance_charge, ion_block, molecules_near,
                         monovalent_ions, net_charge, order, place, rewritten_definition, salt_summary, write_system)
import charmm36_library as c36


def check_changes(host, changes, label):
    """Group (molecule, resid, target) changes by protein moleculetype; validate the targets exist."""
    grouped = {}
    for number, (molecule, resid, target) in enumerate(changes):
        where = f"{label}[{number}]"
        if molecule not in host.definitions or not any(m.name == molecule for m in host.molecules):
            raise ConfigError(f"{where}: no moleculetype {molecule} in the host")
        if next(m for m in host.molecules if m.name == molecule).category != "protein":
            raise ConfigError(f"{where}: {molecule} is not a protein")
        definition = host.definitions[molecule]
        names = [a.atom_name for a in definition.atoms if a.resnr == resid]
        if not names:
            raise ConfigError(f"{where}: {molecule} has no residue {resid} (numbers as in its ITP)")
        resname = next(a.resname for a in definition.atoms if a.resnr == resid)
        grouped.setdefault(molecule, []).append((resid, target, resname, residue_state(resname, names)))
        if [r for r, *_ in grouped[molecule]].count(resid) > 1:
            raise ConfigError(f"{where}: residue {resid} of {molecule} is listed twice")
    return grouped


def run_changes(stage, host, report, cfg, mode, changes, apply, purpose, remove_lipids=False):
    """apply(itp, positions, resid, target, environment) -> (positions, change report) for one residue."""
    work = stage.parent / "work"
    work.mkdir()
    box = host.universe.dimensions
    grouped = check_changes(host, changes, f"{mode}.{'residues' if mode == 'protonate' else 'mutations'}")
    replacements, new_definitions, texts, change_reports = {}, {}, [], []
    new_atom_positions = []
    for molecule, items in grouped.items():
        definition = host.definitions[molecule]
        source = Path(definition.source_path).read_text()
        instances = [m for m in host.molecules if m.name == molecule]
        results = []
        for number, mol in enumerate(instances):
            itp = ProteinItp(source)
            positions = host.universe.atoms[mol.start:mol.stop].positions.astype(float)
            others = [m for m in host.molecules if m is not mol and m.category not in SOLVENT_CATEGORIES]
            environment = _heavy(host.universe.atoms[_indices(others)]).positions.astype(float) if others \
                else np.zeros((0, 3))
            reports = []
            for resid, target, resname, state in items:
                before = {a.name for a in itp.atoms if a.resnr == resid}
                positions, change = apply(itp, positions, resid, target, environment, mol)
                change["moleculetype"], change["copy"] = molecule, number + 1
                reports.append(change)
                for i, atom in enumerate(itp.atoms):
                    if atom.resnr == resid and (atom.name not in before or atom.name in change.get("atoms_rebuilt", ())):
                        new_atom_positions.append(positions[i])
            results.append((itp.text(), positions, reports))
        if len({text for text, _, _ in results}) != 1:
            raise ConfigError(f"The copies of {molecule} gave different topologies")
        text = results[0][0]
        new_definition = rewritten_definition(text, work, molecule)
        new_definitions[molecule] = new_definition
        texts.append(text)
        for mol, (_, positions, reports) in zip(instances, results):
            block = place([(new_definition, "protein", positions)], work, label=f"{molecule}_{mol.occurrence}")
            replacements[id(mol)] = (block, block.molecules[0])
            change_reports += reports
    for change in change_reports:
        if change["copy"] == 1:
            print(f"[{mode}] {change['moleculetype']} {change['from']} {change['resnr']} -> {change['to']}{change['to_label']}"
                  f": +{len(change['atoms_added'])} / -{len(change['atoms_removed'])} atoms, charge "
                  f"{change['charge_change']:+g} e")
    # Parameters the edited residues need and the host lacks.
    types = {a.atom_type for d in host.definitions.values() for a in d.atoms}
    supplement_text, supplement = c36.supplement(host.files, texts, types)
    parameter_files = []
    if supplement:
        path = work / "residue_parameters.itp"
        path.write_text(supplement_text)
        parameter_files.append(path)
        print("[parameters] added from the CHARMM36 library: " + ", ".join(f"{k} {n}" for k, n in supplement.items()))
    # Environment overlapping the new atoms.
    new_atoms = np.array(new_atom_positions) if new_atom_positions else np.zeros((0, 3))
    solvent = [m for m in host.molecules if m.category in SOLVENT_CATEGORIES]
    clashing = molecules_near(host, solvent, new_atoms, cfg.water_clash_distance, box) if len(new_atoms) else set()
    removed = [m for i, m in enumerate(solvent) if i in clashing]
    solvent = [m for i, m in enumerate(solvent) if i not in clashing]
    membrane = [m for m in host.molecules if m.category in MEMBRANE_CATEGORIES]
    lipid_clash = molecules_near(host, membrane, new_atoms, getattr(cfg, "lipid_clash_distance", 0.0), box) \
        if len(new_atoms) and getattr(cfg, "lipid_clash_distance", 0.0) > 0 else set()
    if lipid_clash and not remove_lipids:
        names = Counter(membrane[i].name for i in lipid_clash)
        report["warnings"].append("New atoms lie close to " + ", ".join(f"{n} {k}" for k, n in names.items())
                                  + "; minimize before MD")
        lipid_clash = set()
    removed += [m for i, m in enumerate(membrane) if i in lipid_clash]
    membrane = [m for i, m in enumerate(membrane) if i not in lipid_clash]
    others = [m for m in host.molecules if m.category not in SOLVENT_CATEGORIES + MEMBRANE_CATEGORIES]
    kept = others + membrane + solvent
    net = net_charge([replacements.get(id(m), (host, m)) for m in others]) \
        + net_charge([(host, m) for m in membrane + solvent])
    report["charge"] = {"before": round(net, 6), "after": round(net, 6), "removed_ions": [], "added_ions": {}}
    added = []
    if cfg.neutralize:
        cation, anion = monovalent_ions(host)
        if cation is None or anion is None:
            raise ConfigError("The host has no monovalent cation and anion to neutralize with")
        protein_heavy = _heavy(host.universe.atoms[_indices(others)]).positions
        solvent, added, report["charge"] = balance_charge(host, solvent, net, cfg.target_net_charge, cation, anion,
                                                          [(protein_heavy, 8.0)], box, seed=cfg.seed)
        kept = others + membrane + solvent
        print(f"[charge] {report['charge']['before']:+.3f} e -> {report['charge']['after']:+.3f} e")
    ordered = order(host, kept, new=ion_block(added, work), replacements=replacements)
    report[mode] = {"changes": change_reports, "parameters_added": supplement,
                    "removed": dict(Counter(m.name for m in removed))}
    report["salt"] = {"host": salt_summary([(host, m) for m in host.molecules]), "final": salt_summary(ordered)}
    edited = {id(block_mol) for block, block_mol in replacements.values()}
    report.update(write_system(stage, host, ordered, mode, parameter_files=parameter_files,
                               equilibration=cfg.equilibration, purpose=purpose,
                               scene=({"inserted_edited_protein": edited},
                                      {"inserted_edited_protein": ("Edited protein", [.3, .75, .95], True)},
                                      "Coordinates unchanged except the edited residues."),
                               restrained=[d for d in new_definitions.values()]
                               + [m.definition for m in others if m.name not in new_definitions]))


def report_lines(mode, heading):
    def lines(report):
        info = report.get(mode)
        if not info:
            return []
        out = ["", heading]
        for change in info["changes"]:
            if change["copy"] != 1:
                continue
            out.append(f"  {change['moleculetype']} {change['from']} {change['resnr']} -> {change['to']}"
                       f"{change['to_label']}: atoms added {', '.join(change['atoms_added']) or 'none'}; removed "
                       f"{', '.join(change['atoms_removed']) or 'none'}; charge {change['charge_change']:+g} e")
            if change.get("rotamer"):
                out.append(f"    side chain: chi {change['rotamer']['chi']} (clash score "
                           f"{change['rotamer']['score']:.2f}, {change['rotamer']['tried']} rotamers tried)")
        copies = {c["copy"] for c in info["changes"]}
        if len(copies) > 1:
            out.append(f"  applied to {len(copies)} copies of each edited protein")
        if info.get("parameters_added"):
            out.append("  CHARMM36 parameters added: " + ", ".join(f"{k} {n}" for k, n in info["parameters_added"].items()))
        out.append("  removed: " + (", ".join(f"{k} {n}" for k, n in info["removed"].items()) or "none"))
        return out
    return lines
