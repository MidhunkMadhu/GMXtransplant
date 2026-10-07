"""Protein/ligand classification and harmonic heavy-atom restraints.

Restrained: protein and ligand heavy atoms, and lipid head groups (the polar
heavy atoms of lipids and sterols: phosphorus, oxygen, nitrogen, sulfur).
Lipid tail carbons, hydrogens, water and ions stay free to relax.
"""
from collections import Counter

import numpy as np

from classify import classify_resnames


def molecule_atom_roles(resnames, overrides):
    """Classify a topology molecule, including protein caps and modifications.

    Unrecognized molecules are ligands. Known environment molecules are left
    unrestrained. Explicit residue overrides always take precedence.
    """
    classes = classify_resnames(resnames, overrides)
    known = set(classes.values()) - {"other"}
    if "protein" in known:
        default = "protein"
    elif "ligand" in known:
        default = "ligand"
    elif known:
        default = sorted(known)[0]
    else:
        default = "ligand"
    return [overrides.get(name, default) for name in resnames]


def add_heavy_atom_restraints(openmm, unit, top, system, positions, spec):
    roles = []
    inferred = set()
    for name, count in top._molecules:
        names = [atom[3] for atom in top._moleculeTypes[name].atoms]
        molecule_roles = molecule_atom_roles(names, spec.restraint_residue_classes)
        classes = classify_resnames(names, spec.restraint_residue_classes)
        inferred.update(rn for rn, role in zip(names, molecule_roles)
                        if role == "ligand" and classes[rn] == "other")
        roles.extend(molecule_roles * count)
    atoms = list(top.topology.atoms())
    if len(roles) != len(atoms) or len(atoms) != system.getNumParticles():
        raise ValueError("Cannot map restraint selection to topology particles")
    periodic = system.usesPeriodicBoundaryConditions()
    distance = ("periodicdistance(x,y,z,x0,y0,z0)^2" if periodic else
                "((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    force = openmm.CustomExternalForce("0.5*restraint_k*" + distance)
    force.setName("Protein/ligand heavy-atom and lipid head-group positional restraints")
    force.addGlobalParameter("restraint_k", spec.restraint_force_constant_kj_mol_nm2)
    for parameter in ("x0", "y0", "z0"):
        force.addPerParticleParameter(parameter)
    reference = np.asarray(positions.value_in_unit(unit.nanometer), dtype=float)
    selected, counts = [], Counter()
    for index, (atom, role) in enumerate(zip(atoms, roles)):
        if role not in {"protein", "ligand", "lipid"}:
            continue
        if system.getParticleMass(index).value_in_unit(unit.dalton) == 0:
            continue  # Virtual/dummy sites are not independent physical atoms.
        if atom.element is None:
            raise ValueError(f"Unknown element for restrained atom {index + 1} ({atom.name})")
        if atom.element.atomic_number == 1:
            continue
        if role == "lipid" and atom.element.atomic_number == 6:
            continue  # Carbon chains (tails, ring bodies) relax; polar head groups stay.
        force.addParticle(index, reference[index].tolist())
        selected.append(index)
        counts[role] += 1
    if not selected:
        raise ValueError("No protein, ligand or lipid head-group atoms found for positional restraints")
    system.addForce(force)
    return {"force_constant_kj_mol_nm2": spec.restraint_force_constant_kj_mol_nm2,
            "reference": "input coordinates", "periodic": periodic,
            "heavy_atom_count": len(selected),
            "protein_heavy_atoms": counts["protein"],
            "ligand_heavy_atoms": counts["ligand"],
            "lipid_headgroup_atoms": counts["lipid"],
            "inferred_ligand_resnames": sorted(inferred),
            "residue_class_overrides": spec.restraint_residue_classes}, selected
