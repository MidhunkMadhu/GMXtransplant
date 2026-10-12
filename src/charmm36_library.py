"""The bundled CHARMM36 protein and ion library, and the parameters a host is missing.

CHARMM-GUI writes a forcefield.itp that holds only the parameters its own system
uses. A new ion type, a protonation state or a mutation can need atom types or
bonded parameters that file lacks. This module finds the rows a set of molecule
ITPs needs, the way grompp would look them up (exact type match first, then the
CHARMM wildcards), and returns only those the host does not define. Rows the host
defines are never replaced, so every molecule already in the host keeps its
parameters.

Data: charmm36_protein_ions.json.gz, built by tools/build_charmm36_library.py from
the July 2022 CHARMM36 GROMACS port and top_all36_prot.rtf (MacKerell lab).
"""
from __future__ import annotations

from functools import lru_cache
from importlib.resources import files
import gzip
import json
import re

import numpy as np

from config import ConfigError

LIBRARY = "charmm36_protein_ions.json.gz"
SECTION_ORDER = ("atomtypes", "nonbond_params", "pairtypes", "bondtypes", "constrainttypes", "angletypes",
                 "dihedraltypes", "cmaptypes")
TYPE_COLUMNS = {"atomtypes": 1, "pairtypes": 2, "bondtypes": 2, "constrainttypes": 2, "angletypes": 3,
                "dihedraltypes": 4, "cmaptypes": 5, "nonbond_params": 2}


@lru_cache(None)
def library():
    with gzip.open(str(files("gmxtransplant").joinpath("data", LIBRARY)), "rt", encoding="utf-8") as handle:
        return json.load(handle)


def residue(name):
    residues = library()["residues"]
    if name not in residues:
        raise ConfigError(f"No CHARMM36 template for residue {name}; available: {', '.join(sorted(residues))}")
    return residues[name]


def ion(name):
    ions = library()["ions"]
    if name not in ions:
        raise ConfigError(f"No CHARMM36 ion {name} in the bundled library; available: {', '.join(sorted(ions))}")
    return ions[name]


def source():
    return library()["source"]


# ---------------------------------------------------------------------------
# Parameter rows
# ---------------------------------------------------------------------------

def parameter_rows(paths):
    """Parameter-section rows of force-field files: {section: [tokens]}."""
    rows = {}
    for path in paths:
        section, continuation = None, ""
        for raw in open(path, encoding="utf-8", errors="replace"):
            line = raw.split(";", 1)[0].strip()
            if not line or line.startswith("#"):
                continue
            match = re.fullmatch(r"\[\s*([a-z_]+)\s*\]", line)
            if match:
                section = match[1]
                continue
            if section not in TYPE_COLUMNS:
                continue
            continuation += " " + line.rstrip("\\").strip()
            if line.endswith("\\"):
                continue
            rows.setdefault(section, []).append(continuation.split())
            continuation = ""
    return rows


def _key(section, row):
    """Rows defining the same parameter share a key; dihedrals also by function."""
    width = TYPE_COLUMNS[section]
    types = tuple(row[:width])
    if section in ("pairtypes", "bondtypes", "constrainttypes", "nonbond_params", "angletypes", "dihedraltypes"):
        types = min(types, types[::-1])
    if section == "dihedraltypes":
        return types + (row[4],)
    return types


def _score(pattern, types):
    """Specific (non-X) matches of a wildcard pattern to atom types, or -1 if it does not match."""
    if len(pattern) != len(types) or any(p != "X" and p != t for p, t in zip(pattern, types)):
        return -1
    return sum(p != "X" for p in pattern)


class _Index:
    def __init__(self, rows):
        self.rows = rows
        self.by_key = {}
        for section, items in rows.items():
            for row in items:
                self.by_key.setdefault((section, _key(section, row)), []).append(row)
        self.dihedral_keys = {}
        for section, key in self.by_key:
            if section == "dihedraltypes":
                self.dihedral_keys.setdefault(key[4], []).append(key[:4])

    def has(self, section, key):
        return (section, key) in self.by_key

    def best_dihedral(self, types, function):
        """The key grompp would use for a dihedral of these types: most specific match, either direction."""
        best, score = None, -1
        for pattern in self.dihedral_keys.get(function, []):
            for candidate in (types, types[::-1]):
                s = _score(pattern, candidate)
                if s > score:
                    best, score = pattern, s
        return None if best is None else best + (function,)


@lru_cache(None)
def _library_index():
    return _Index(library()["parameters"])


def molecule_terms(itp_text):
    """Atom types and type tuples of every bonded term in one molecule ITP text."""
    types, section = {}, None
    terms = {"atomtypes": set(), "bonds": set(), "pairs": set(), "angles": set(), "dihedrals9": set(),
             "dihedrals2": set(), "cmap": set()}
    for raw in itp_text.splitlines():
        line = raw.split(";", 1)[0].strip()
        if not line or line.startswith("#"):
            continue
        match = re.fullmatch(r"\[\s*([a-z_]+)\s*\]", line)
        if match:
            section = match[1]
            continue
        fields = line.split()
        if section == "atoms":
            types[fields[0]] = fields[1]
            terms["atomtypes"].add(fields[1])
        elif section in ("bonds", "pairs") and len(fields) >= 2:
            terms[section].add(tuple(types[i] for i in fields[:2]))
        elif section == "angles" and len(fields) >= 3:
            terms["angles"].add(tuple(types[i] for i in fields[:3]))
        elif section == "dihedrals" and len(fields) >= 5:
            terms["dihedrals" + fields[4]] = terms.get("dihedrals" + fields[4], set())
            terms["dihedrals" + fields[4]].add(tuple(types[i] for i in fields[:4]))
        elif section == "cmap" and len(fields) >= 5:
            terms["cmap"].add(tuple(types[i] for i in fields[:5]))
    return terms


def supplement(host_parameter_files, itp_texts, system_types=()):
    """Parameter rows the new molecules need and the host lacks.

    host_parameter_files: the host's force-field files; itp_texts: the new or edited molecule
    ITPs; system_types: every atom type in the final system (for NBFIX pairs of new types).
    Returns (text of a parameter-only ITP, report). Raises ConfigError when the bundled
    library lacks a needed parameter too.
    """
    host = _Index(parameter_rows(host_parameter_files))
    lib = _library_index()
    needed, missing = {}, []

    def add(section, key):
        rows = lib.by_key.get((section, key))
        if rows:
            needed.setdefault(section, {})[key] = rows
            return True
        return False

    terms = {}
    for text in itp_texts:
        for name, values in molecule_terms(text).items():
            terms.setdefault(name, set()).update(values)
    host_types = {k[0] for (s, k) in host.by_key if s == "atomtypes"}
    new_types = set()
    for t in sorted(terms.get("atomtypes", ())):
        if t not in host_types:
            if not add("atomtypes", (t,)):
                missing.append(f"atom type {t}")
            new_types.add(t)
    for section, term in (("bondtypes", "bonds"), ("angletypes", "angles")):
        for types in sorted(terms.get(term, ())):
            key = min(types, types[::-1])
            if not host.has(section, key) and not add(section, key):
                missing.append(f"{term} {'-'.join(types)}")
    for function in ("9", "2"):
        for types in sorted(terms.get("dihedrals" + function, ())):
            best = lib.best_dihedral(types, function)
            if best is None:
                if host.best_dihedral(types, function) is None:
                    missing.append(f"dihedral ({function}) {'-'.join(types)}")
                continue
            if not host.has("dihedraltypes", min(best[:4], best[:4][::-1]) + (function,)):
                add("dihedraltypes", min(best[:4], best[:4][::-1]) + (function,))
    for types in sorted(terms.get("cmap", ())):
        if not host.has("cmaptypes", types) and not add("cmaptypes", types):
            missing.append(f"cmap {'-'.join(types)}")
    # 1-4 pairs listed in the new molecules that use CHARMM's special parameters; NBFIX pairs of new types.
    for types in sorted(terms.get("pairs", ())):
        key = min(types, types[::-1])
        if not host.has("pairtypes", key):
            add("pairtypes", key)
    everything = set(system_types) | host_types | set(terms.get("atomtypes", ()))
    for t in new_types:
        for other in everything:
            key = min((t, other), (other, t))
            if not host.has("nonbond_params", key):
                add("nonbond_params", key)
    # A pair type equal to the combination rule (what gen-pairs makes from the atom types) adds nothing.
    lj = {}
    for rows in (lib.rows.get("atomtypes", []), host.rows.get("atomtypes", [])):
        for row in rows:
            lj[row[0]] = (float(row[5]), float(row[6]))
    for key in list(needed.get("pairtypes", {})):
        row = needed["pairtypes"][key][0]
        if key[0] in lj and key[1] in lj:
            sigma = (lj[key[0]][0] + lj[key[1]][0]) / 2
            epsilon = (lj[key[0]][1] * lj[key[1]][1]) ** 0.5
            if abs(float(row[3]) - sigma) <= 1e-6 * sigma and abs(float(row[4]) - epsilon) <= 1e-6 * max(epsilon, 1e-12):
                del needed["pairtypes"][key]
    if missing:
        raise ConfigError("Neither the host forcefield.itp nor the bundled CHARMM36 library has parameters for: "
                          + ", ".join(missing[:8]) + (f" (+{len(missing) - 8} more)" if len(missing) > 8 else ""))
    lines = ["; Parameters added by GMXtransplant from the bundled CHARMM36 library",
             f"; ({source()})", "; only rows the host forcefield.itp does not define"]
    report = {}
    for section in SECTION_ORDER:
        rows = [row for key in sorted(needed.get(section, {})) for row in needed[section][key]]
        if not rows:
            continue
        report[section] = len(rows)
        lines.append(f"\n[ {section} ]")
        for row in rows:
            if section == "cmaptypes" and len(row) > 8:
                head, grid = " ".join(row[:8]), row[8:]
                lines.append("\\\n".join([head] + [" ".join(grid[i:i + 10]) for i in range(0, len(grid), 10)]))
            else:
                lines.append(" ".join(row))
    return "\n".join(lines) + "\n", report


def atomtype_mass(name, host_parameter_files=()):
    """Mass of an atom type, from the host force field or the library."""
    for path_rows in (parameter_rows(host_parameter_files).get("atomtypes", []),
                      library()["parameters"]["atomtypes"]):
        for row in path_rows:
            if row[0] == name:
                return float(row[2])
    raise ConfigError(f"Atom type {name} has no mass in the host force field or the CHARMM36 library")


def ion_itp(name):
    """A one-atom moleculetype ITP for a library ion."""
    item = ion(name)
    return ("; CHARMM36 ion from the bundled library (GROMACS port of July 2022)\n"
            f"[ moleculetype ]\n; name nrexcl\n{name} 1\n\n[ atoms ]\n"
            "; nr type resnr residu atom cgnr charge mass\n"
            f"     1 {item['type']:>6}      1 {item['resname']:>6} {item['atom']:>6}      1 "
            f"{item['charge']:>10.6f} {item['mass']:>10.4f}\n")


# ---------------------------------------------------------------------------
# Internal coordinates
# ---------------------------------------------------------------------------

def _place(a, b, c, bond, angle, torsion):
    """Position d with |cd| = bond, angle(b, c, d) = angle and dihedral(a, b, c, d) = torsion (degrees)."""
    angle, torsion = np.radians(angle), np.radians(torsion)
    bc = c - b
    bc /= np.linalg.norm(bc)
    n = np.cross(b - a, bc)
    norm = np.linalg.norm(n)
    if norm < 1e-8:
        raise ConfigError("Internal-coordinate reference atoms are collinear")
    n /= norm
    m = np.cross(n, bc)
    d2 = np.array([-bond * np.cos(angle), bond * np.sin(angle) * np.cos(torsion),
                   bond * np.sin(angle) * np.sin(torsion)])
    return c + d2[0] * bc + d2[1] * m + d2[2] * n


def dihedral(a, b, c, d):
    b0, b1, b2 = a - b, c - b, d - c
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w))))


def _filled(ic, types, host_parameter_files):
    """CHARMM leaves some IC values 0 (patches): take them from the bond and angle parameters."""
    i, j, k, l = (n.lstrip("*") for n in ic[:4])
    improper = ic[2].startswith("*")
    r_ij, t_ijk, phi, t_jkl, r_kl = ic[4:9]
    rows = None

    def bond(x, y):
        nonlocal rows
        rows = rows or {**library()["parameters"], **parameter_rows(host_parameter_files)}
        tx, ty = types.get(x), types.get(y)
        for row in rows.get("bondtypes", []):
            if {row[0], row[1]} == {tx, ty} or (row[0], row[1]) in ((tx, ty), (ty, tx)):
                return float(row[3]) * 10.0
        raise ConfigError(f"No bond parameter for {x}-{y} to complete an internal coordinate")

    def angle(x, y, z):
        nonlocal rows
        rows = rows or {**library()["parameters"], **parameter_rows(host_parameter_files)}
        tx, ty, tz = types.get(x), types.get(y), types.get(z)
        for row in rows.get("angletypes", []):
            if row[1] == ty and {row[0], row[2]} == {tx, tz}:
                return float(row[4])
        raise ConfigError(f"No angle parameter for {x}-{y}-{z} to complete an internal coordinate")

    if r_ij == 0:
        r_ij = bond(i, k) if improper else bond(i, j)
    if t_ijk == 0:
        t_ijk = angle(i, k, j) if improper else angle(i, j, k)
    if t_jkl == 0:
        t_jkl = angle(j, k, l)
    if r_kl == 0:
        r_kl = bond(k, l)
    return r_ij, t_ijk, phi, t_jkl, r_kl


def build_atoms(known, ic_table, wanted, types, host_parameter_files=()):
    """Place the `wanted` atoms from a residue's IC table.

    known: {name: xyz} for atoms with coordinates, neighbours as '-C', '+N', '+CA'.
    types: {name: CHARMM atom type}, used to fill zero IC values from parameters.
    Returns {name: xyz} for the placed atoms; raises if one cannot be placed.
    """
    known = {k: np.asarray(v, dtype=float) for k, v in known.items()}
    wanted = set(wanted)
    placed = {}
    progress = True
    while progress and wanted - set(placed):
        progress = False
        for ic in ic_table:
            names = [n.lstrip("*") for n in ic[:4]]
            i, j, k, l = names
            improper = ic[2].startswith("*")
            have = {n: (n in known or n in placed) for n in names}
            pos = {**known, **placed}
            if l in wanted and l not in placed and have[i] and have[j] and have[k]:
                r_ij, t_ijk, phi, t_jkl, r_kl = _filled(ic, types, host_parameter_files)
                placed[l] = _place(pos[i], pos[j], pos[k], r_kl, t_jkl, phi)
                progress = True
            elif i in wanted and i not in placed and have[j] and have[k] and have[l]:
                r_ij, t_ijk, phi, t_jkl, r_kl = _filled(ic, types, host_parameter_files)
                if improper:
                    # I is bonded to K: angle I-K-J and the torsion I-J-K-L about the J-K axis.
                    candidate = _place(pos[l], pos[j], pos[k], r_ij, t_ijk, -phi)
                    if abs(((dihedral(candidate, pos[j], pos[k], pos[l]) - phi) + 180) % 360 - 180) > 1e-3:
                        candidate = _place(pos[l], pos[j], pos[k], r_ij, t_ijk, phi)
                    placed[i] = candidate
                else:
                    placed[i] = _place(pos[l], pos[k], pos[j], r_ij, t_ijk, phi)
                progress = True
    left = wanted - set(placed)
    if left:
        raise ConfigError(f"Could not place atoms {sorted(left)} from the CHARMM36 internal coordinates")
    return placed
