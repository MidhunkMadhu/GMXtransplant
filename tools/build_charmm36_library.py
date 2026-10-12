#!/usr/bin/env python3
"""Build src/gmxtransplant/data/charmm36_protein_ions.json.gz from the official CHARMM36 files.

Inputs (from the MacKerell lab, https://mackerell.umaryland.edu/charmm_ff.shtml):
  charmm36-jul2022.ff/            the GROMACS port (aminoacids.rtp, ffnonbonded.itp,
                                  ffbonded.itp, cmap.itp, nbfix.itp, ions.itp)
  toppar/top_all36_prot.rtf       from toppar_c36_jul24.tgz, for the internal coordinates (IC)

    python3 tools/build_charmm36_library.py charmm36-jul2022.ff toppar/top_all36_prot.rtf

The library holds the standard amino acids and their protonation variants (atoms,
types, charges, bonds, impropers, CMAP, IC tables), the CHARMM36 ions, and every
parameter row whose atom types all belong to proteins, ions or TIP3P water.
"""
import gzip
import json
from pathlib import Path
import re
import sys

RESIDUES = ["ALA", "ARG", "ASN", "ASP", "ASPP", "CYS", "CYM", "GLN", "GLU", "GLUP", "GLY", "HSD", "HSE", "HSP",
            "ILE", "LEU", "LYS", "LSN", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
# Patches of top_all36_prot.rtf applied on a base residue's IC table.
PATCHED = {"ASPP": "ASP", "GLUP": "GLU", "LSN": "LYS"}
IONS = ["SOD", "POT", "CLA", "CAL", "MG", "LIT", "CES", "RUB", "BAR", "ZN2", "CD2"]
WATER_TYPES = {"OT", "HT"}
SECTIONS = ("atomtypes", "pairtypes", "bondtypes", "constrainttypes", "angletypes", "dihedraltypes",
            "cmaptypes", "nonbond_params")
# Number of atom-type columns that key each section's rows.
TYPE_COLUMNS = {"atomtypes": 1, "pairtypes": 2, "bondtypes": 2, "constrainttypes": 2, "angletypes": 3,
                "dihedraltypes": 4, "cmaptypes": 5, "nonbond_params": 2}


def parse_rtp(path):
    residues, current, section = {}, None, None
    for raw in Path(path).read_text().splitlines():
        line = raw.split(";", 1)[0].strip()
        if not line:
            continue
        match = re.fullmatch(r"\[\s*(\S+)\s*\]", line)
        if match:
            name = match[1]
            if name in ("atoms", "bonds", "impropers", "cmap", "dihedrals", "exclusions"):
                section = name
            else:
                current = residues.setdefault(name, {"atoms": [], "bonds": [], "impropers": [], "cmap": []})
                section = None
            continue
        if current is None or section is None:
            continue
        fields = line.split()
        if section == "atoms":
            current["atoms"].append([fields[0], fields[1], float(fields[2]), int(fields[3])])
        elif section in ("bonds", "impropers", "cmap"):
            current[section].append(fields)
    return residues


def parse_rtf_ic(path):
    """IC lines per RESI/PRES: [I, J, K, L, R(IJ), T(IJK), PHI, T(JKL), R(KL)]; '*K' marks an improper IC."""
    ics, deleted, current = {}, {}, None
    for raw in Path(path).read_text().splitlines():
        line = raw.split("!", 1)[0].strip()
        fields = line.split()
        if not fields:
            continue
        keyword = fields[0].upper()
        if keyword in ("RESI", "PRES"):
            current = fields[1]
            ics[current], deleted[current] = [], []
        elif keyword == "END":
            current = None
        elif current and keyword == "IC" and len(fields) >= 10:
            ics[current].append(fields[1:5] + [float(v) for v in fields[5:10]])
        elif current and keyword == "DELETE" and len(fields) >= 3 and fields[1].upper() == "ATOM":
            deleted[current].append(fields[2])
    return ics, deleted


def parse_parameter_files(paths):
    rows = {section: [] for section in SECTIONS}
    defaults = None
    for path in paths:
        section, continuation = None, ""
        for raw in Path(path).read_text().splitlines():
            line = raw.split(";", 1)[0].strip()
            if not line or line.startswith("#"):
                continue
            match = re.fullmatch(r"\[\s*([a-z_]+)\s*\]", line)
            if match:
                section = match[1]
                continue
            if section == "defaults":
                defaults = line.split()
                continue
            if section not in rows:
                continue
            continuation += " " + line.rstrip("\\").strip()
            if line.endswith("\\"):
                continue
            rows[section].append(continuation.split())
            continuation = ""
    return defaults, rows


def parse_ions(path):
    ions, name, section = {}, None, None
    for raw in Path(path).read_text().splitlines():
        line = raw.split(";", 1)[0].strip()
        if not line:
            continue
        match = re.fullmatch(r"\[\s*(\S+)\s*\]", line)
        if match:
            section = match[1]
            continue
        fields = line.split()
        if section == "moleculetype":
            name = fields[0]
        elif section == "atoms" and name:
            ions[name] = {"type": fields[1], "charge": float(fields[6]), "mass": float(fields[7]),
                          "atom": fields[4], "resname": fields[3]}
    return ions


def main(ff_dir, rtf_path, output):
    ff = Path(ff_dir)
    rtp = parse_rtp(ff / "aminoacids.rtp")
    ics, deleted = parse_rtf_ic(rtf_path)
    residues = {}
    for name in RESIDUES:
        entry = dict(rtp[name])
        base = PATCHED.get(name)
        table = list(ics[base]) if base else list(ics[name])
        if base:
            removed = set(deleted.get(name, []))
            table = [ic for ic in table if not removed & set(n.lstrip("*") for n in ic[:4])] + ics[name]
        entry["ic"] = table
        residues[name] = entry
    ions = {name: item for name, item in parse_ions(ff / "ions.itp").items() if name in IONS}
    types = {atom[1] for r in residues.values() for atom in r["atoms"]} | {i["type"] for i in ions.values()} \
        | WATER_TYPES
    defaults, rows = parse_parameter_files([ff / "ffnonbonded.itp", ff / "ffbonded.itp", ff / "cmap.itp",
                                            ff / "nbfix.itp"])
    kept = {}
    for section, items in rows.items():
        width = TYPE_COLUMNS[section]
        kept[section] = [row for row in items if all(t == "X" or t in types for t in row[:width])]
    document = {
        "source": ("CHARMM36 (July 2022 GROMACS port, charmm2gmx, https://mackerell.umaryland.edu/"
                   "charmm_ff.shtml) for residues, ions and parameters; internal coordinates from "
                   "top_all36_prot.rtf of toppar_c36_jul24. MacKerell lab, University of Maryland."),
        "defaults": defaults, "residues": residues, "ions": ions, "parameters": kept,
    }
    with gzip.open(output, "wt", encoding="utf-8") as handle:
        json.dump(document, handle, separators=(",", ":"))
    print(f"{output}: {len(residues)} residues, {len(ions)} ions, "
          + ", ".join(f"{s} {len(r)}" for s, r in kept.items()))


if __name__ == "__main__":
    if len(sys.argv) not in (3, 4):
        sys.exit(__doc__)
    target = sys.argv[3] if len(sys.argv) == 4 else str(Path(__file__).resolve().parent.parent / "src" /
                                                          "gmxtransplant" / "data" / "charmm36_protein_ions.json.gz")
    main(sys.argv[1], sys.argv[2], target)
