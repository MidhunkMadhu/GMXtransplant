"""Maximum common substructure (MCS) matching for ligand placement.

Used by ligand_replace ``fit.method: mcsfit``. The heavy atoms of the original
and incoming ligands are turned into bond graphs; the largest connected set of
atoms the two share -- same element on each matched atom, identical bonds
between every matched pair -- is found by exhaustive search. Symmetric
alternatives (a ring matched clockwise or anticlockwise, equivalent substituent
positions) are all kept, and the one with the lowest rigid-fit RMSD onto the
original ligand's pose is used.

Bonds come from each ligand's ITP when one is available, otherwise from
interatomic distances with covalent radii. Elements come from ITP masses or,
without an ITP, from atom names.
"""

from __future__ import annotations

import glob
import os
import re
from typing import Dict, List, Optional, Sequence, Set, Tuple

import numpy as np

from kabsch import kabsch_fit

# Nearest-mass element assignment for ITP [ atoms ] masses.
_MASSES = {"H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998, "P": 30.974,
           "S": 32.06, "CL": 35.45, "BR": 79.904, "I": 126.904, "B": 10.81, "SE": 78.97}
# Covalent radii (Å) for distance-based bond detection.
_RADII = {"H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "P": 1.07, "S": 1.05,
          "CL": 1.02, "BR": 1.20, "I": 1.39, "B": 0.84, "SE": 1.20}
SEARCH_BUDGET = 2_000_000
MAX_MAPPINGS = 5000


class SubstructureError(RuntimeError):
    pass


class LigandGraph:
    """Heavy-atom bond graph of one ligand, keyed by atom name."""

    def __init__(self, names: Sequence[str], elements: Sequence[str], bonds: Set[Tuple[int, int]], source: str):
        self.names, self.elements, self.source = list(names), list(elements), source
        self.neighbors: List[Set[int]] = [set() for _ in self.names]
        for i, j in bonds:
            self.neighbors[i].add(j)
            self.neighbors[j].add(i)

    def __len__(self):
        return len(self.names)


def _element_from_mass(mass: float) -> str:
    return min(_MASSES, key=lambda element: abs(_MASSES[element] - mass))


def _element_from_name(name: str) -> str:
    letters = re.sub(r"[^A-Za-z]", "", name).upper()
    if letters[:2] in ("CL", "BR", "SE"):
        return letters[:2]
    return letters[:1] or "X"


def _read_itp_molecules(path: str):
    """{moleculetype: (atom names, masses, bonds as 0-based index pairs)} from one ITP."""
    molecules, current, section = {}, None, None
    with open(path, encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.split(";", 1)[0].strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith("["):
                section = line.strip("[] ").lower()
                continue
            fields = line.split()
            if section == "moleculetype":
                current = fields[0]
                molecules[current] = ([], [], [])
            elif current is None:
                continue
            elif section == "atoms" and len(fields) >= 8:
                molecules[current][0].append(fields[4])
                molecules[current][1].append(float(fields[7]))
            elif section == "bonds" and len(fields) >= 2:
                molecules[current][2].append((int(fields[0]) - 1, int(fields[1]) - 1))
    return molecules


def _graph_from_itp(atom_names: Sequence[str], path: str, molecule: Tuple) -> LigandGraph:
    names, masses, bonds = molecule
    heavy = [i for i, mass in enumerate(masses) if _element_from_mass(mass) != "H"]
    position = {old: new for new, old in enumerate(heavy)}
    heavy_bonds = {(position[i], position[j]) for i, j in bonds if i in position and j in position}
    return LigandGraph([names[i] for i in heavy], [_element_from_mass(masses[i]) for i in heavy],
                       heavy_bonds, f"bonds from {path}")


def _find_itp(atom_names: Sequence[str], candidates: Sequence[str]):
    """The ITP moleculetype whose atom names equal these coordinates' names, in order."""
    for path in candidates:
        try:
            molecules = _read_itp_molecules(path)
        except (OSError, ValueError, IndexError):
            continue
        for molecule in molecules.values():
            if list(molecule[0]) == list(atom_names):
                return path, molecule
    return None, None


def ligand_graph(atom_group, itp_paths: Sequence[str] = (), search_dirs: Sequence[str] = ()) -> LigandGraph:
    names = [str(n) for n in atom_group.names]
    candidates = [p for p in itp_paths if p and os.path.isfile(p)]
    for folder in search_dirs:
        if folder and os.path.isdir(folder):
            candidates += sorted(glob.glob(os.path.join(folder, "**", "*.itp"), recursive=True))
    path, molecule = _find_itp(names, candidates)
    if molecule is not None:
        return _graph_from_itp(names, path, molecule)
    # No matching ITP: detect bonds from the coordinates.
    elements = [_element_from_name(n) for n in names]
    heavy = [i for i, e in enumerate(elements) if e != "H"]
    positions = atom_group.positions.astype(float)[heavy]
    bonds = set()
    for a in range(len(heavy)):
        for b in range(a + 1, len(heavy)):
            ra = _RADII.get(elements[heavy[a]], 0.8)
            rb = _RADII.get(elements[heavy[b]], 0.8)
            if np.linalg.norm(positions[a] - positions[b]) <= ra + rb + 0.45:
                bonds.add((a, b))
    return LigandGraph([names[i] for i in heavy], [elements[i] for i in heavy], bonds,
                       "bonds detected from interatomic distances")


def maximum_common_substructures(g1: LigandGraph, g2: LigandGraph, budget: int = SEARCH_BUDGET):
    """All largest connected, bond-preserving element matches between g1 and g2.

    Returns a list of mappings, each a tuple of (g1 index, g2 index) pairs.
    """
    best_size, best = 0, {}
    nodes = 0

    def compatible(a, b, mapping):
        for x, y in mapping.items():
            if (x in g1.neighbors[a]) != (y in g2.neighbors[b]):
                return False
        return True

    def extend(mapping, used, excluded):
        nonlocal best_size, nodes
        nodes += 1
        if nodes > budget:
            raise SubstructureError(
                "The substructure search between these ligands is too large to finish. "
                "Use fit.method: pairfit with explicit atom lists instead.")
        frontier = sorted({n for x in mapping for n in g1.neighbors[x]} - set(mapping) - excluded)
        # Upper bound: every atom still reachable could at most be matched once.
        if len(mapping) + len(set(range(len(g1))) - set(mapping) - excluded) < best_size:
            return
        if not frontier:
            size = len(mapping)
            if size > best_size:
                best_size = size
                best.clear()
            if size == best_size and len(best) < MAX_MAPPINGS:
                best[frozenset(mapping.items())] = None
            return
        a = frontier[0]
        options = [b for b in {n for y in mapping.values() for n in g2.neighbors[y]} - used
                   if g2.elements[b] == g1.elements[a] and compatible(a, b, mapping)]
        for b in sorted(options):
            mapping[a] = b
            used.add(b)
            extend(mapping, used, excluded)
            del mapping[a]
            used.discard(b)
        excluded.add(a)          # or leave atom a out of the common substructure
        extend(mapping, used, excluded)
        excluded.discard(a)

    for a in range(len(g1)):
        for b in range(len(g2)):
            if g1.elements[a] == g2.elements[b]:
                # Seeds below `a` were already explored as seeds themselves.
                extend({a: b}, {b}, set(range(a)))
    return [tuple(sorted(m)) for m in best]


def best_substructure_fit(old_ag, new_ag, old_graph: LigandGraph, new_graph: LigandGraph):
    """Pick the largest common substructure with the lowest fit RMSD onto the old pose.

    Returns (old names, new names, rmsd, number of equivalent mappings tried).
    """
    mappings = maximum_common_substructures(old_graph, new_graph)
    if not mappings or len(mappings[0]) < 3:
        raise SubstructureError(
            "The original and new ligands share fewer than three bonded heavy atoms of the same "
            "element, so no common substructure can place the new ligand. Use fit.method: pairfit "
            "or nofit.")
    old_positions = {str(n): p for n, p in zip(old_ag.names, old_ag.positions.astype(float))}
    new_positions = {str(n): p for n, p in zip(new_ag.names, new_ag.positions.astype(float))}
    best = None
    for mapping in mappings:
        old_names = [old_graph.names[i] for i, _ in mapping]
        new_names = [new_graph.names[j] for _, j in mapping]
        A = np.array([new_positions[n] for n in new_names])
        B = np.array([old_positions[n] for n in old_names])
        if np.linalg.matrix_rank(A - A.mean(0), tol=1e-6) < 2:
            continue
        rmsd = kabsch_fit(A, B)[3]
        if best is None or rmsd < best[2] - 1e-9:
            best = (old_names, new_names, float(rmsd))
    if best is None:
        raise SubstructureError("Every common substructure found is collinear; use fit.method: pairfit.")
    return best + (len(mappings),)
