"""Replace residues inside a protein ITP with CHARMM36 templates.

Used by the protonate and mutate modes. A residue R is rebuilt from a template
of the bundled CHARMM36 library: its atoms get the template's types and charges,
and every bonded term that touches R beyond its backbone N and C (the atoms the
neighbouring residues bond to) is regenerated the way pdb2gmx does for CHARMM36:
bonds and impropers from the template, angles, proper dihedrals and 1-4 pairs
from the bond graph, CMAP from the template. Terms that touch R only through N
or C belong to the neighbours and are kept unchanged. Position restraints follow
CHARMM-GUI: POSRES_FC_BB on N, CA, C, O and POSRES_FC_SC on the other heavy atoms.

New atom coordinates come from the template's CHARMM internal coordinates.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
import re

import numpy as np

from config import ConfigError
import charmm36_library as c36

INDEX_COLUMNS = {"bonds": 2, "pairs": 2, "angles": 3, "dihedrals": 4, "cmap": 5, "position_restraints": 1,
                 "constraints": 2, "dihedral_restraints": 4, "distance_restraints": 2, "settles": 1}
BACKBONE = ("N", "CA", "C", "O")
# Residue names CHARMM-GUI may keep for another protonation state (e.g. ASP with HD2 = ASPP).
VARIANTS = {"ASP": ("ASP", "ASPP"), "ASPP": ("ASP", "ASPP"), "GLU": ("GLU", "GLUP"), "GLUP": ("GLU", "GLUP"),
            "LYS": ("LYS", "LSN"), "LSN": ("LYS", "LSN"), "HSD": ("HSD", "HSE", "HSP"),
            "HSE": ("HSD", "HSE", "HSP"), "HSP": ("HSD", "HSE", "HSP"), "HIS": ("HSD", "HSE", "HSP"),
            "CYS": ("CYS", "CYM"), "CYM": ("CYS", "CYM")}


def residue_state(resname, names):
    """The CHARMM36 template a residue's atoms actually match (ASP with HD2 is ASPP)."""
    names = set(names)
    for candidate in VARIANTS.get(resname, (resname,)):
        try:
            template = c36.residue(candidate)
        except ConfigError:
            continue
        if {a[0] for a in template["atoms"]} == names:
            return candidate
    return None
ANCHORS = ("N", "C")


@dataclass
class Atom:
    nr: int
    type: str
    resnr: int
    resname: str
    name: str
    cgnr: int
    charge: float
    mass: float


@dataclass
class Item:
    kind: str                 # "text" (kept verbatim), "section", "row", "directive"
    text: str = ""
    section: str = ""
    tokens: list = field(default_factory=list)
    comment: str = ""


class ProteinItp:
    """One-moleculetype ITP kept as a sequence of items so it can be written back faithfully."""

    def __init__(self, text):
        self.items, self.atoms = [], []
        section, seen_molecule = None, False
        for raw in text.splitlines():
            stripped = raw.strip()
            body, _, comment = raw.partition(";")
            match = re.match(r"\s*\[\s*([^]]+?)\s*\]", body)
            if match:
                section = match[1].lower()
                if section == "moleculetype":
                    if seen_molecule:
                        raise ConfigError("The protein ITP defines more than one moleculetype")
                    seen_molecule = True
                self.items.append(Item("section", raw, section))
            elif stripped.startswith("#"):
                self.items.append(Item("directive", raw, section))
            elif not body.strip():
                self.items.append(Item("text", raw, section))
            elif section == "atoms":
                f = body.split()
                self.atoms.append(Atom(int(f[0]), f[1], int(f[2]), f[3], f[4], int(f[5]), float(f[6]), float(f[7])))
                self.items.append(Item("atom", "", section))
            else:
                self.items.append(Item("row", raw, section, body.split(), comment))
        if not self.atoms:
            raise ConfigError("The protein ITP has no [ atoms ]")

    # -- queries ------------------------------------------------------------
    def residues(self):
        """Ordered residues: resnr -> list of atom indices (0-based)."""
        out = defaultdict(list)
        for i, atom in enumerate(self.atoms):
            out[atom.resnr].append(i)
        return dict(out)

    def bonds(self):
        return [(int(it.tokens[0]) - 1, int(it.tokens[1]) - 1) for it in self.items
                if it.kind == "row" and it.section == "bonds"]

    # -- writing ------------------------------------------------------------
    def text(self):
        out, qtot, written = [], 0.0, False
        for item in self.items:
            if item.kind == "atom":
                # The atom list may have grown or shrunk: write all of it at the first [ atoms ] row.
                if written:
                    continue
                written = True
                for atom in self.atoms:
                    qtot += atom.charge
                    out.append(f"{atom.nr:6d} {atom.type:>10} {atom.resnr:6d} {atom.resname:>8} {atom.name:>6} "
                               f"{atom.cgnr:6d} {atom.charge:12.6f} {atom.mass:10.4f}   ; qtot {qtot:6.3f}")
            elif item.kind == "row" and item.tokens:
                out.append(" ".join(f"{t:>5}" for t in item.tokens) + (f" ;{item.comment}" if item.comment else ""))
            else:
                out.append(item.text)
        return "\n".join(out) + "\n"


def _index_columns(item):
    if item.section == "dihedrals" or item.section in INDEX_COLUMNS:
        return INDEX_COLUMNS.get(item.section, 0)
    return 0


def _paths(graph, atoms):
    """Angles, proper dihedrals and 1-4 pairs of the bond graph that touch `atoms`."""
    angles, dihedrals = set(), set()
    for j in graph:
        for i in graph[j]:
            for k in graph[j]:
                if i < k:
                    triple = (i, j, k)
                    if atoms & set(triple):
                        angles.add(triple)
    for j in graph:
        for k in graph[j]:
            if j >= k:
                continue
            for i in graph[j]:
                if i == k:
                    continue
                for l in graph[k]:
                    if l in (i, j):
                        continue
                    quad = (i, j, k, l)
                    if atoms & set(quad):
                        dihedrals.add(quad if quad < quad[::-1] else quad[::-1])
    close = {(a, b) for a in graph for b in graph[a]}
    close |= {(i, k) for i, _, k in angles} | {(k, i) for i, _, k in angles}
    # 1-3 pairs anywhere in the graph (rings: a 1-4 pair can also be a 1-3 pair)
    for j in graph:
        for i in graph[j]:
            for k in graph[j]:
                if i != k:
                    close.add((i, k))
    pairs = set()
    for i, _, _, l in dihedrals:
        # A pair belongs to the residue when one of its two atoms does (removal uses the same rule).
        if (i, l) in close or i == l or not ({i, l} & atoms):
            continue
        pairs.add((min(i, l), max(i, l)))
    return angles, dihedrals, pairs


def replace_residue(itp, resnr, template_name, positions, keep=None, host_parameter_files=(), build=True,
                    ic_table=None, coordinates_only=False):
    """Rebuild residue `resnr` of `itp` (a ProteinItp) from a CHARMM36 template.

    positions: (n_atoms, 3) coordinates in ITP order. keep: names whose coordinates are kept
    (default: every name the old and new residue share). ic_table: the template's IC table
    with changed values (a rotamer). Returns (new positions, report); atoms are renumbered in
    place and `itp` is modified. With coordinates_only, returns ({name: xyz} of the new
    residue, {old atom indices of the residue and its peptide neighbours}) and leaves `itp` alone.
    """
    template = c36.residue(template_name)
    residues = itp.residues()
    if resnr not in residues:
        raise ConfigError(f"Residue {resnr} is not in the protein ITP")
    old = residues[resnr]
    old_names = [itp.atoms[i].name for i in old]
    old_resname = itp.atoms[old[0]].resname
    old_charge = sum(itp.atoms[i].charge for i in old)
    old_state = residue_state(old_resname, old_names)
    if old_state is None:
        raise ConfigError(f"Residue {old_resname} {resnr} ({', '.join(old_names)}) matches no CHARMM36 residue: a "
                          "terminal, capped or patched residue cannot be changed")
    bonds = itp.bonds()
    index_of = {itp.atoms[i].name: i for i in old}
    old_set = set(old)
    # Neighbours through the peptide bonds.
    prev_c = next_n = None
    for a, b in bonds:
        for x, y in ((a, b), (b, a)):
            if x == index_of.get("N") and y not in old_set and itp.atoms[y].name == "C":
                prev_c = y
            if x == index_of.get("C") and y not in old_set and itp.atoms[y].name == "N":
                next_n = y
    if prev_c is None or next_n is None:
        raise ConfigError(f"Residue {old_resname} {resnr} is at a chain end; terminal residues cannot be changed")
    anchors = {index_of[n] for n in ANCHORS}
    for a, b in bonds:
        for x, y in ((a, b), (b, a)):
            if x in old_set and y not in old_set and x not in anchors:
                raise ConfigError(f"Atom {itp.atoms[x].name} of {old_resname} {resnr} is bonded to another residue "
                                  f"({itp.atoms[y].resname} {itp.atoms[y].resnr} {itp.atoms[y].name}): a disulfide or "
                                  "link cannot be changed")
    new_atoms = template["atoms"]
    new_names = [a[0] for a in new_atoms]
    for required in BACKBONE:
        if required not in new_names:
            raise ConfigError(f"Template {template_name} lacks backbone atom {required}")
    keep = set(old_names) & set(new_names) if keep is None else set(keep) & set(old_names) & set(new_names)
    keep |= set(ANCHORS)

    # -- coordinates ---------------------------------------------------------
    coords = {n: positions[index_of[n]] for n in keep}
    neighbours = {"-C": positions[prev_c], "+N": positions[next_n]}
    for a, b in bonds:
        for x, y in ((a, b), (b, a)):
            if x == prev_c and itp.atoms[y].name == "CA" and itp.atoms[y].resnr != resnr:
                neighbours["-CA"] = positions[y]
            if x == next_n and itp.atoms[y].name == "CA":
                neighbours["+CA"] = positions[y]
    types = {a[0]: a[1] for a in new_atoms}
    built = {}
    missing = [n for n in new_names if n not in coords]
    if missing and build:
        built = c36.build_atoms({**coords, **neighbours}, ic_table or template["ic"], missing, types,
                                host_parameter_files)
    if coordinates_only:
        near = set(old) | {prev_c, next_n} | {y for a, b in bonds for x, y in ((a, b), (b, a))
                                                 if x in (prev_c, next_n)}
        return {**coords, **built}, near

    # -- atoms ---------------------------------------------------------------
    first, last = min(old), max(old)
    mapping = {}
    new_list = []
    for i, atom in enumerate(itp.atoms[:first]):
        mapping[i] = len(new_list)
        new_list.append(atom)
    uniform_cg = all(a.cgnr == a.nr for a in itp.atoms)
    start = len(new_list)
    residue_atoms = []
    for name, atype, charge, group in new_atoms:
        mass = c36.atomtype_mass(atype, host_parameter_files)
        atom = Atom(0, atype, resnr, template_name, name, group, charge, mass)
        if name in keep:
            mapping[index_of[name]] = len(new_list)
        residue_atoms.append(atom)
        new_list.append(atom)
    for i in range(last + 1, len(itp.atoms)):
        mapping[i] = len(new_list)
        new_list.append(itp.atoms[i])
    if old != list(range(first, last + 1)):
        raise ConfigError(f"The atoms of residue {resnr} are not consecutive in the ITP")
    # renumber; charge groups: one per atom (CHARMM-GUI) or shifted template groups
    group_base = itp.atoms[first - 1].cgnr if first else 0
    for n, atom in enumerate(new_list):
        if uniform_cg:
            atom.cgnr = n + 1
        elif atom in residue_atoms:
            atom.cgnr = group_base + atom.cgnr
        atom.nr = n + 1
    if not uniform_cg:
        last_group = max(a.cgnr for a in residue_atoms)
        shift = last_group - itp.atoms[last].cgnr
        for atom in new_list[start + len(residue_atoms):]:
            atom.cgnr += shift
    edited = {start + k for k, (name, *_rest) in enumerate(new_atoms) if name not in ANCHORS}
    new_index = {name: start + k for k, (name, *_rest) in enumerate(new_atoms)}
    removed_atoms = [n for n in old_names if n not in new_names]

    # -- bonded terms --------------------------------------------------------
    kept_rows = []
    old_edited = {index_of[n] for n in old_names if n not in ANCHORS}
    positions_new = np.zeros((len(new_list), 3))
    for i, j in mapping.items():
        positions_new[j] = positions[i]
    for name, xyz in built.items():
        positions_new[new_index[name]] = xyz
    for item in itp.items:
        if item.kind != "row":
            kept_rows.append(item)
            continue
        width = _index_columns(item)
        if not width:
            kept_rows.append(item)
            continue
        indices = [int(t) - 1 for t in item.tokens[:width]]
        if set(indices) & old_edited or (item.section == "position_restraints" and set(indices) & old_set):
            continue  # regenerated below
        if any(i not in mapping for i in indices):
            continue
        item.tokens = [str(mapping[i] + 1) for i in indices] + item.tokens[width:]
        kept_rows.append(item)
    itp.items = kept_rows
    itp.atoms = new_list

    def resolve(name):
        """Template atom names, with -C and +N for the neighbours' peptide atoms."""
        if name == "-C":
            return mapping[prev_c]
        if name == "+N":
            return mapping[next_n]
        return None if name[0] in "+-" else new_index.get(name)

    graph = defaultdict(set)
    for item in itp.items:
        if item.kind == "row" and item.section == "bonds":
            a, b = int(item.tokens[0]) - 1, int(item.tokens[1]) - 1
            graph[a].add(b)
            graph[b].add(a)
    new_bonds = []
    for a, b in template["bonds"]:
        ia, ib = resolve(a), resolve(b)
        if ia is None or ib is None:
            continue
        if not ({ia, ib} & edited):
            continue
        new_bonds.append((ia, ib))
        graph[ia].add(ib)
        graph[ib].add(ia)
    angles, dihedrals, pairs = _paths(graph, edited)
    impropers = []
    for quad in template["impropers"]:
        ids = [resolve(n) for n in quad]
        if any(i is None for i in ids):
            raise ConfigError(f"Cannot resolve improper {quad} of {template_name}")
        if set(ids) & edited:
            impropers.append(ids)
    cmaps = []
    for five in template["cmap"]:
        ids = [resolve(n) for n in five]
        if None not in ids and set(ids) & edited:
            cmaps.append(ids)
    posres = []
    for name, atype, *_ in new_atoms:
        if atype.upper().startswith("H"):
            continue
        posres.append((new_index[name], "POSRES_FC_BB" if name in BACKBONE else "POSRES_FC_SC"))

    def append(section, rows, predicate=lambda item: True):
        """Add rows after the last row of the first matching section block."""
        target = None
        for n, item in enumerate(itp.items):
            if item.kind == "section" and item.section == section:
                block = n
                last_row = n
                for m in range(n + 1, len(itp.items)):
                    other = itp.items[m]
                    if other.kind == "section":
                        break
                    if other.kind == "row":
                        last_row = m
                rows_in = [itp.items[m] for m in range(n + 1, last_row + 1) if itp.items[m].kind == "row"]
                if not rows_in or predicate(rows_in[0]):
                    target = (block, last_row)
                    break
        if target is None:
            itp.items.append(Item("section", f"[ {section} ]", section))
            target = (len(itp.items) - 1, len(itp.items) - 1)
        itp.items[target[1] + 1:target[1] + 1] = [Item("row", "", section, list(map(str, r))) for r in rows]

    append("bonds", [(a + 1, b + 1, 1) for a, b in sorted(new_bonds)])
    append("pairs", [(a + 1, b + 1, 1) for a, b in sorted(pairs)])
    append("angles", [(a + 1, b + 1, c + 1, 5) for a, b, c in sorted(angles)])
    append("dihedrals", [(a + 1, b + 1, c + 1, d + 1, 9) for a, b, c, d in sorted(dihedrals)],
           lambda item: len(item.tokens) > 4 and item.tokens[4] == "9")
    append("dihedrals", [(a + 1, b + 1, c + 1, d + 1, 2) for a, b, c, d in impropers],
           lambda item: len(item.tokens) > 4 and item.tokens[4] == "2")
    append("cmap", [(*[i + 1 for i in ids], 1) for ids in cmaps])
    if any(it.kind == "row" and it.section == "position_restraints" for it in itp.items):
        append("position_restraints", [(i + 1, 1, macro, macro, macro) for i, macro in posres])
    report = {"resnr": resnr, "from": old_state, "residue_name_in_input": old_resname, "to": template_name, "atoms_removed": removed_atoms,
              "atoms_added": sorted(set(new_names) - set(old_names)),
              "atoms_rebuilt": sorted(n for n in set(new_names) & set(old_names) if n not in keep),
              "charge_change": round(sum(a[2] for a in new_atoms) - old_charge, 4),
              "regenerated_terms": {"bonds": len(new_bonds), "pairs": len(pairs), "angles": len(angles),
                                    "dihedrals": len(dihedrals), "impropers": len(impropers), "cmap": len(cmaps)}}
    return positions_new, report

