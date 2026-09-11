#!/usr/bin/env python3
"""Convert cholesterol residues in a full PDB to CHARMM36 atom format.

Run ``python3 charmm36_cholesterol.py --help`` for the command-line interface.
"""

from __future__ import annotations

import argparse
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from itertools import permutations
from typing import Iterable, Sequence

import numpy as np


# This is the 74-atom CHARMM36 cholesterol reference used when -d/--reference
# is omitted. Atom names, atom order, hydrogen count, and coordinates are all
# embedded so the default does not depend on an external reference PDB.
DEFAULT_REFERENCE_PDB = """\
ATOM      1  C3  CHL1    1     105.630  53.526 100.147  1.00  0.00      MEMB C
ATOM      2  H3  CHL1    1     106.475  53.055  99.588  1.00  0.00      MEMB H
ATOM      3  O3  CHL1    1     105.854  53.341 101.534  1.00  0.00      MEMB O
ATOM      4  H3' CHL1    1     105.168  53.848 101.969  1.00  0.00      MEMB H
ATOM      5  C4  CHL1    1     105.579  55.036  99.850  1.00  0.00      MEMB C
ATOM      6  H4A CHL1    1     104.942  55.574 100.576  1.00  0.00      MEMB H
ATOM      7  H4B CHL1    1     106.632  55.289 100.063  1.00  0.00      MEMB H
ATOM      8  C5  CHL1    1     105.224  55.427  98.428  1.00  0.00      MEMB C
ATOM      9  C6  CHL1    1     105.883  56.390  97.773  1.00  0.00      MEMB C
ATOM     10  H6  CHL1    1     106.534  57.068  98.320  1.00  0.00      MEMB H
ATOM     11  C7  CHL1    1     105.816  56.685  96.316  1.00  0.00      MEMB C
ATOM     12  H7A CHL1    1     105.901  57.766  96.085  1.00  0.00      MEMB H
ATOM     13  H7B CHL1    1     106.680  56.204  95.816  1.00  0.00      MEMB H
ATOM     14  C8  CHL1    1     104.573  56.084  95.685  1.00  0.00      MEMB C
ATOM     15  H8  CHL1    1     103.733  56.767  95.902  1.00  0.00      MEMB H
ATOM     16  C14 CHL1    1     104.791  56.002  94.201  1.00  0.00      MEMB C
ATOM     17  H14 CHL1    1     105.552  55.193  94.042  1.00  0.00      MEMB H
ATOM     18  C15 CHL1    1     105.387  57.224  93.493  1.00  0.00      MEMB C
ATOM     19 H15A CHL1    1     104.909  58.164  93.842  1.00  0.00      MEMB H
ATOM     20 H15B CHL1    1     106.481  57.318  93.646  1.00  0.00      MEMB H
ATOM     21  C16 CHL1    1     105.006  56.938  92.041  1.00  0.00      MEMB C
ATOM     22 H16A CHL1    1     104.390  57.783  91.687  1.00  0.00      MEMB H
ATOM     23 H16B CHL1    1     105.875  56.803  91.389  1.00  0.00      MEMB H
ATOM     24  C17 CHL1    1     104.138  55.653  92.007  1.00  0.00      MEMB C
ATOM     25  H17 CHL1    1     104.842  54.790  91.912  1.00  0.00      MEMB H
ATOM     26  C13 CHL1    1     103.557  55.581  93.421  1.00  0.00      MEMB C
ATOM     27  C18 CHL1    1     102.443  56.655  93.569  1.00  0.00      MEMB C
ATOM     28 H18A CHL1    1     101.570  56.411  92.936  1.00  0.00      MEMB H
ATOM     29 H18B CHL1    1     102.090  56.725  94.615  1.00  0.00      MEMB H
ATOM     30 H18C CHL1    1     102.797  57.664  93.274  1.00  0.00      MEMB H
ATOM     31  C12 CHL1    1     103.156  54.186  93.937  1.00  0.00      MEMB C
ATOM     32 H12A CHL1    1     102.205  53.846  93.487  1.00  0.00      MEMB H
ATOM     33 H12B CHL1    1     103.959  53.472  93.646  1.00  0.00      MEMB H
ATOM     34  C11 CHL1    1     103.004  54.141  95.467  1.00  0.00      MEMB C
ATOM     35 H11A CHL1    1     102.122  54.762  95.730  1.00  0.00      MEMB H
ATOM     36 H11B CHL1    1     102.792  53.099  95.778  1.00  0.00      MEMB H
ATOM     37  C9  CHL1    1     104.254  54.676  96.205  1.00  0.00      MEMB C
ATOM     38  H9  CHL1    1     105.121  54.037  95.903  1.00  0.00      MEMB H
ATOM     39  C10 CHL1    1     104.171  54.560  97.746  1.00  0.00      MEMB C
ATOM     40  C19 CHL1    1     102.850  55.209  98.204  1.00  0.00      MEMB C
ATOM     41 H19A CHL1    1     102.705  56.142  97.636  1.00  0.00      MEMB H
ATOM     42 H19B CHL1    1     101.983  54.558  97.987  1.00  0.00      MEMB H
ATOM     43 H19C CHL1    1     102.834  55.443  99.289  1.00  0.00      MEMB H
ATOM     44  C1  CHL1    1     104.183  53.071  98.195  1.00  0.00      MEMB C
ATOM     45  H1A CHL1    1     103.240  52.592  97.858  1.00  0.00      MEMB H
ATOM     46  H1B CHL1    1     105.029  52.556  97.688  1.00  0.00      MEMB H
ATOM     47  C2  CHL1    1     104.328  52.873  99.713  1.00  0.00      MEMB C
ATOM     48  H2A CHL1    1     104.332  51.792  99.968  1.00  0.00      MEMB H
ATOM     49  H2B CHL1    1     103.480  53.355 100.248  1.00  0.00      MEMB H
ATOM     50  C20 CHL1    1     103.223  55.617  90.781  1.00  0.00      MEMB C
ATOM     51  H20 CHL1    1     102.626  56.558  90.745  1.00  0.00      MEMB H
ATOM     52  C21 CHL1    1     102.250  54.427  90.881  1.00  0.00      MEMB C
ATOM     53 H21A CHL1    1     101.521  54.598  91.697  1.00  0.00      MEMB H
ATOM     54 H21B CHL1    1     102.808  53.488  91.077  1.00  0.00      MEMB H
ATOM     55 H21C CHL1    1     101.656  54.302  89.957  1.00  0.00      MEMB H
ATOM     56  C22 CHL1    1     104.060  55.543  89.476  1.00  0.00      MEMB C
ATOM     57 H22A CHL1    1     104.784  56.386  89.442  1.00  0.00      MEMB H
ATOM     58 H22B CHL1    1     104.652  54.602  89.480  1.00  0.00      MEMB H
ATOM     59  C23 CHL1    1     103.231  55.614  88.187  1.00  0.00      MEMB C
ATOM     60 H23A CHL1    1     102.474  54.803  88.165  1.00  0.00      MEMB H
ATOM     61 H23B CHL1    1     102.690  56.585  88.159  1.00  0.00      MEMB H
ATOM     62  C24 CHL1    1     104.097  55.468  86.931  1.00  0.00      MEMB C
ATOM     63 H24A CHL1    1     104.721  56.379  86.800  1.00  0.00      MEMB H
ATOM     64 H24B CHL1    1     104.789  54.607  87.072  1.00  0.00      MEMB H
ATOM     65  C25 CHL1    1     103.275  55.206  85.654  1.00  0.00      MEMB C
ATOM     66  H25 CHL1    1     102.758  54.224  85.770  1.00  0.00      MEMB H
ATOM     67  C26 CHL1    1     102.216  56.290  85.396  1.00  0.00      MEMB C
ATOM     68 H26A CHL1    1     102.698  57.289  85.336  1.00  0.00      MEMB H
ATOM     69 H26B CHL1    1     101.699  56.094  84.433  1.00  0.00      MEMB H
ATOM     70 H26C CHL1    1     101.457  56.311  86.205  1.00  0.00      MEMB H
ATOM     71  C27 CHL1    1     104.214  55.118  84.452  1.00  0.00      MEMB C
ATOM     72 H27A CHL1    1     104.687  56.111  84.311  1.00  0.00      MEMB H
ATOM     73 H27B CHL1    1     105.014  54.367  84.620  1.00  0.00      MEMB H
ATOM     74 H27C CHL1    1     103.658  54.857  83.527  1.00  0.00      MEMB H
END
"""

# A displayed token such as ``CHL1R`` in a whitespace-split PDB line normally
# means residue name CHL1 followed immediately by chain R. The fixed-width PDB
# residue-name field itself is only four characters, so CHL1R is not included.
DEFAULT_CHOL_RESNAMES = ("CHL1", "CHL", "CHOL", "CLR")
H_BOND_CUTOFF = 1.25
HEAVY_BOND_CUTOFF = 1.85
DEFAULT_VERIFY_TOLERANCE = 0.40

COVALENT_RADII = {
    "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
    "P": 1.07, "S": 1.05, "CL": 1.02, "BR": 1.20, "I": 1.39,
    "NA": 1.66, "MG": 1.41, "K": 2.03, "CA": 1.76, "ZN": 1.22,
    "FE": 1.32, "MN": 1.39, "CU": 1.32,
}


class PipelineError(RuntimeError):
    """A user-facing conversion error."""


def build_parser() -> argparse.ArgumentParser:
    description = """\
Convert cholesterol residues in a protein or full-system PDB to the atom
naming, atom order, hydrogen count, connectivity, and equivalent-hydrogen
label orientation of CHARMM36 cholesterol.

The program finds all recognized cholesterol residues, removes their existing
hydrogens, uses Open Babel to regenerate hydrogens from the 3D structure, maps
each residue onto a CHARMM36 cholesterol topology, removes excess hydrogens,
renames and reorders all cholesterol atoms, assigns geometrically equivalent
hydrogens such as H4A/H4B by comparison with the aligned reference, and
independently verifies the result. Protein, water, ions, lipids, and other
non-cholesterol atom records are carried into the outputs unchanged, apart
from global atom serial renumbering. Stale CONECT records are removed because
their serial numbers are no longer valid after rebuilding and reordering
cholesterol.

If -d/--reference is omitted, the script uses its built-in 74-atom CHARMM36
cholesterol reference. The built-in atom names, order, hydrogen count, and
coordinates are identical to the default reference used to develop this tool.
"""
    epilog = """\
OUTPUTS
  Unless overridden, every output is placed beside INPUT_PDB and uses the
  input filename stem:

    <stem>_protonated_raw.pdb
        Intermediate full-system PDB after Open Babel protonation, before
        CHARMM36 pruning, naming, and ordering.

    <stem>_fixed.pdb
        Main result. Cholesterol residues use CHARMM36 atom names, order,
        hydrogen count, connectivity, and equivalent-hydrogen label
        orientation.

    <stem>_reference_aligned_to_first_chol.pdb
        One copy of the reference cholesterol rigidly aligned to the first
        cholesterol residue found in the input.

    <stem>_reference_placed.pdb
        Full-system PDB in which every detected cholesterol is replaced by a
        rigidly aligned copy of the reference cholesterol.

EXAMPLES
  Use the built-in CHARMM36 cholesterol reference:
    python3 charmm36_cholesterol.py -f receptor_with_cholesterol.pdb

  Use an external reference:
    python3 charmm36_cholesterol.py -f system.pdb -d chl.pdb

  Override the main and intermediate output paths:
    python3 charmm36_cholesterol.py -f system.pdb \\
        -o system_fixed.pdb -r system_raw.pdb

  Override every generated path:
    python3 charmm36_cholesterol.py -f system.pdb -d chl.pdb \\
        -r raw.pdb -o fixed.pdb -a aligned_reference.pdb \\
        -p reference_placed.pdb

DEPENDENCIES
  Python 3, NumPy, and Open Babel's 'obabel' executable are required.
  Open Babel can commonly be installed with:
    conda install -c conda-forge openbabel
  or, on Debian/Ubuntu:
    sudo apt-get install openbabel
"""
    parser = argparse.ArgumentParser(
        prog="charmm36_cholesterol.py",
        description=description,
        epilog=epilog,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    paths = parser.add_argument_group("input and reference")
    paths.add_argument(
        "-f", "--input", required=True, metavar="INPUT_PDB", type=Path,
        help="protein or full-system PDB containing one or more cholesterols",
    )
    paths.add_argument(
        "-d", "--reference", metavar="REFERENCE_PDB", type=Path,
        help=("single-residue CHARMM36 cholesterol reference PDB; if omitted, "
              "the built-in reference is used"),
    )

    outputs = parser.add_argument_group("output paths")
    outputs.add_argument(
        "-r", "--raw-out", metavar="PDB", type=Path,
        help="raw protonated intermediate; default: INPUT_STEM_protonated_raw.pdb",
    )
    outputs.add_argument(
        "-o", "--fixed-out", metavar="PDB", type=Path,
        help="main CHARMM36-formatted output; default: INPUT_STEM_fixed.pdb",
    )
    outputs.add_argument(
        "-a", "--aligned-reference-out", metavar="PDB", type=Path,
        help=("reference aligned to the first input cholesterol; default: "
              "INPUT_STEM_reference_aligned_to_first_chol.pdb"),
    )
    outputs.add_argument(
        "-p", "--reference-placed-out", metavar="PDB", type=Path,
        help=("full system with aligned reference cholesterols; default: "
              "INPUT_STEM_reference_placed.pdb"),
    )

    advanced = parser.add_argument_group("advanced options")
    advanced.add_argument(
        "--chol-resnames", metavar="NAMES",
        default=",".join(DEFAULT_CHOL_RESNAMES),
        help=("comma-separated cholesterol residue names to detect; default: "
              + ",".join(DEFAULT_CHOL_RESNAMES)),
    )
    advanced.add_argument(
        "--obabel", default="obabel", metavar="COMMAND",
        help="Open Babel executable name or path; default: obabel",
    )
    advanced.add_argument(
        "--verify-tolerance", type=float, default=DEFAULT_VERIFY_TOLERANCE,
        metavar="ANGSTROM",
        help="extra distance tolerance used for final bond verification; default: 0.40",
    )
    advanced.add_argument(
        "--no-reconstruct-missing-heavy-atoms",
        dest="reconstruct_missing_heavy_atoms",
        action="store_false",
        help="fail instead of fitting the reference to reconstruct missing heavy atoms",
    )
    advanced.set_defaults(reconstruct_missing_heavy_atoms=True)
    advanced.add_argument(
        "--max-missing-heavy-atoms", type=int, default=2, metavar="N",
        help="maximum heavy atoms reconstructed in one cholesterol; default: 2",
    )
    advanced.add_argument(
        "--max-heavy-atom-fit-rmsd", type=float, default=1.5,
        metavar="ANGSTROM",
        help="maximum observed-heavy-atom RMSD for reconstruction; default: 1.5",
    )
    return parser


def derive_output_paths(args: argparse.Namespace) -> None:
    source = args.input.expanduser()
    base = source.with_suffix("") if source.suffix.lower() == ".pdb" else source
    args.input = source
    args.reference = args.reference.expanduser() if args.reference else None
    args.raw_out = (args.raw_out or Path(f"{base}_protonated_raw.pdb")).expanduser()
    args.fixed_out = (args.fixed_out or Path(f"{base}_fixed.pdb")).expanduser()
    args.aligned_reference_out = (
        args.aligned_reference_out
        or Path(f"{base}_reference_aligned_to_first_chol.pdb")
    ).expanduser()
    args.reference_placed_out = (
        args.reference_placed_out
        or Path(f"{base}_reference_placed.pdb")
    ).expanduser()


def parse_resnames(value: str) -> set[str]:
    names = {item.strip().upper() for item in value.split(",") if item.strip()}
    if not names:
        raise PipelineError("--chol-resnames must contain at least one residue name")
    too_long = sorted(name for name in names if len(name) > 4)
    if too_long:
        raise PipelineError(
            "PDB residue names can contain at most four characters; invalid: "
            + ", ".join(too_long)
        )
    return names


def validate_paths(args: argparse.Namespace) -> None:
    if not args.input.is_file():
        raise PipelineError(f"input PDB does not exist or is not a file: {args.input}")
    if args.reference and not args.reference.is_file():
        raise PipelineError(
            f"reference PDB does not exist or is not a file: {args.reference}"
        )
    if args.verify_tolerance < 0:
        raise PipelineError("--verify-tolerance cannot be negative")
    if args.max_missing_heavy_atoms < 0:
        raise PipelineError("--max-missing-heavy-atoms cannot be negative")
    if not math.isfinite(args.max_heavy_atom_fit_rmsd) or args.max_heavy_atom_fit_rmsd <= 0:
        raise PipelineError("--max-heavy-atom-fit-rmsd must be finite and > 0")

    outputs = [
        args.raw_out,
        args.fixed_out,
        args.aligned_reference_out,
        args.reference_placed_out,
    ]
    normalized = [os.path.abspath(path) for path in outputs]
    if len(set(normalized)) != len(normalized):
        raise PipelineError("each output option must point to a different file")

    protected = {os.path.abspath(args.input)}
    if args.reference:
        protected.add(os.path.abspath(args.reference))
    collisions = [path for path in outputs if os.path.abspath(path) in protected]
    if collisions:
        raise PipelineError(
            "an output path would overwrite an input/reference file: "
            + ", ".join(map(str, collisions))
        )
    for path in outputs:
        path.parent.mkdir(parents=True, exist_ok=True)


def is_hydrogen(name: str, element: str) -> bool:
    return element.strip().upper() == "H" or name.strip().upper().startswith("H")


def atom_element(atom: dict) -> str:
    element = atom["element"].strip().upper()
    if element:
        return element
    return "H" if atom["is_h"] else atom["name"].strip()[0].upper()


def distance(a: dict, b: dict) -> float:
    return math.sqrt(
        (a["x"] - b["x"]) ** 2
        + (a["y"] - b["y"]) ** 2
        + (a["z"] - b["z"]) ** 2
    )


def read_pdb_records(path: Path) -> list[dict]:
    records: list[dict] = []
    try:
        handle = path.open()
    except OSError as exc:
        raise PipelineError(f"cannot open PDB file {path}: {exc}") from exc

    with handle:
        for index, line in enumerate(handle, start=1):
            record: dict = {
                "idx": index,
                "line": line,
                "is_atom": line.startswith(("ATOM  ", "HETATM")),
            }
            if record["is_atom"]:
                try:
                    name = line[12:16].strip()
                    resname = line[17:21].strip()
                    element = line[76:78].strip() if len(line) >= 78 else ""
                    record.update(
                        serial=int(line[6:11]),
                        name=name,
                        resname=resname,
                        chain=line[21].strip(),
                        resid=int(line[22:26]),
                        icode=line[26].strip() if len(line) > 26 else "",
                        x=float(line[30:38]),
                        y=float(line[38:46]),
                        z=float(line[46:54]),
                        element=element,
                        is_h=is_hydrogen(name, element),
                    )
                except (ValueError, IndexError) as exc:
                    raise PipelineError(
                        f"cannot parse ATOM/HETATM record at {path}:{index}: "
                        f"{line.rstrip()}"
                    ) from exc
            records.append(record)
    return records


def residue_key(atom: dict) -> tuple[str, str, int, str]:
    return atom["chain"], atom["resname"], atom["resid"], atom["icode"]


def format_atom_line(
    serial: int,
    name: str,
    resname: str,
    chain: str,
    resid: int,
    icode: str,
    x: float,
    y: float,
    z: float,
    element: str,
    is_het: bool = True,
) -> str:
    line = list(" " * 80)
    line[0:6] = list("HETATM" if is_het else "ATOM  ")
    line[6:11] = list(f"{serial:5d}")
    atom_name = name[:4] if len(name) >= 4 else f" {name:<3}"[:4]
    line[12:16] = list(atom_name)
    line[17:21] = list(f"{resname:<4}"[:4])
    line[21] = (chain or " ")[:1]
    line[22:26] = list(f"{resid:4d}")
    line[26] = (icode or " ")[:1]
    line[30:38] = list(f"{x:8.3f}")
    line[38:46] = list(f"{y:8.3f}")
    line[46:54] = list(f"{z:8.3f}")
    line[54:60] = list(f"{1.00:6.2f}")
    line[60:66] = list(f"{0.00:6.2f}")
    line[72:76] = list("MEMB")
    line[76:78] = list(f"{element:>2}"[:2])
    return "".join(line) + "\n"


def check_obabel(command: str) -> str:
    resolved = shutil.which(command)
    if resolved is None and Path(command).is_file():
        resolved = command
    if resolved is None:
        raise PipelineError(
            f"Open Babel executable not found: {command!r}\n"
            "Install it with 'conda install -c conda-forge openbabel' or "
            "'sudo apt-get install openbabel', or pass --obabel PATH."
        )
    return resolved


def write_embedded_reference(path: Path) -> None:
    path.write_text(DEFAULT_REFERENCE_PDB)


def group_selected_residues(
    records: Iterable[dict], resnames: set[str]
) -> dict[tuple[str, str, int, str], list[dict]]:
    groups: dict[tuple[str, str, int, str], list[dict]] = defaultdict(list)
    for record in records:
        if record["is_atom"] and record["resname"].upper() in resnames:
            groups[residue_key(record)].append(record)
    return groups


def reference_atoms(records: Sequence[dict], path: Path) -> list[dict]:
    groups: dict[tuple[str, str, int, str], list[dict]] = defaultdict(list)
    for record in records:
        if record["is_atom"]:
            groups[residue_key(record)].append(record)
    if len(groups) != 1:
        raise PipelineError(
            f"reference {path} must contain exactly one residue; found "
            f"{len(groups)}: {list(groups)}"
        )
    atoms = next(iter(groups.values()))
    names = [atom["name"] for atom in atoms]
    if len(names) != len(set(names)):
        duplicates = sorted({name for name in names if names.count(name) > 1})
        raise PipelineError(
            f"reference {path} has duplicate atom names: {', '.join(duplicates)}"
        )
    return atoms


def protonate_heavy_atoms(
    heavy_atoms: Sequence[dict], tag: str, workdir: Path, obabel: str
) -> tuple[list[dict], list[dict]]:
    input_path = workdir / f"{tag}_heavy_in.pdb"
    output_path = workdir / f"{tag}_h_out.pdb"
    with input_path.open("w") as handle:
        for serial, atom in enumerate(heavy_atoms, start=1):
            handle.write(
                format_atom_line(
                    serial, atom["name"], "LIG", "", 1, "",
                    atom["x"], atom["y"], atom["z"], atom_element(atom),
                )
            )
        handle.write("END\n")

    result = subprocess.run(
        [obabel, str(input_path), "-O", str(output_path), "-h"],
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0 or not output_path.exists():
        raise PipelineError(
            f"Open Babel failed for {tag} (exit {result.returncode}):\n"
            f"{result.stderr.strip()}"
        )

    output_atoms = [record for record in read_pdb_records(output_path) if record["is_atom"]]
    output_heavy = [atom for atom in output_atoms if not atom["is_h"]]
    output_hydrogens = [atom for atom in output_atoms if atom["is_h"]]

    by_name: dict[str, list[dict]] = defaultdict(list)
    for atom in output_heavy:
        by_name[atom["name"]].append(atom)
    input_names = [atom["name"] for atom in heavy_atoms]
    output_names = [atom["name"] for atom in output_heavy]
    if sorted(input_names) != sorted(output_names):
        missing = sorted(set(input_names) - set(output_names))
        extra = sorted(set(output_names) - set(input_names))
        raise PipelineError(
            f"Open Babel changed the heavy-atom name set for {tag}; "
            f"missing={missing}, unexpected={extra}"
        )
    duplicates = sorted(name for name, atoms in by_name.items() if len(atoms) > 1)
    if duplicates:
        raise PipelineError(
            f"Open Babel produced duplicate heavy-atom names for {tag}: {duplicates}"
        )
    ordered_heavy = [by_name[atom["name"]][0] for atom in heavy_atoms]
    return ordered_heavy, output_hydrogens


def audit_and_reconstruct_heavy_atoms(
    heavy_atoms: Sequence[dict],
    ref_atoms: Sequence[dict],
    key: tuple[str, str, int, str],
    reconstruct_missing: bool,
    max_missing: int,
    max_fit_rmsd: float,
) -> tuple[list[dict], dict]:
    """Validate heavy atoms and conservatively place named missing atoms."""
    reference_heavy = [atom for atom in ref_atoms if not atom["is_h"]]
    expected_count = len(reference_heavy)
    observed_count = len(heavy_atoms)
    audit = {
        "chain": key[0],
        "resname": key[1],
        "resid": key[2],
        "icode": key[3],
        "expected_heavy_atom_count": expected_count,
        "observed_heavy_atom_count": observed_count,
        "missing_heavy_atom_names": [],
        "modeled_heavy_atom_names": [],
        "reference_fit_rmsd_angstrom": None,
        "modeled_atom_local_fit_rmsd_angstrom": {},
        "status": "complete",
    }
    if observed_count > expected_count:
        raise PipelineError(
            f"cholesterol {key} has {observed_count} heavy atoms; the reference "
            f"has {expected_count}. Remove alternate/duplicate or non-cholesterol atoms."
        )

    observed_names = [atom["name"] for atom in heavy_atoms]
    duplicates = sorted({name for name in observed_names if observed_names.count(name) > 1})
    if duplicates:
        raise PipelineError(
            f"cholesterol {key} has duplicate heavy-atom names: {duplicates}. "
            "Resolve alternate locations before conversion."
        )

    reference_by_name = {atom["name"]: atom for atom in reference_heavy}
    observed_by_name = {atom["name"]: atom for atom in heavy_atoms}
    reference_names = set(reference_by_name)
    observed_name_set = set(observed_by_name)

    # A complete residue can still use a different naming convention; the
    # topology mapper below will identify and rename it by graph structure.
    if observed_count == expected_count:
        if observed_name_set == reference_names:
            ref_xyz = coordinates([reference_by_name[name] for name in observed_names])
            target_xyz = coordinates(heavy_atoms)
            rotation, mobile_centroid, target_centroid = kabsch_transform(
                ref_xyz, target_xyz
            )
            fitted = apply_transform(
                ref_xyz, rotation, mobile_centroid, target_centroid
            )
            audit["reference_fit_rmsd_angstrom"] = float(
                np.sqrt(np.mean(np.sum((fitted - target_xyz) ** 2, axis=1)))
            )
        return list(heavy_atoms), audit

    missing_count = expected_count - observed_count
    if not reconstruct_missing:
        raise PipelineError(
            f"cholesterol {key} has {observed_count}/{expected_count} heavy atoms. "
            "Missing-heavy-atom reconstruction is disabled."
        )
    if missing_count > max_missing:
        raise PipelineError(
            f"cholesterol {key} is missing {missing_count} heavy atoms, exceeding "
            f"max_missing_heavy_atoms={max_missing}."
        )
    unexpected = sorted(observed_name_set - reference_names)
    if unexpected:
        raise PipelineError(
            f"cholesterol {key} is incomplete and uses heavy-atom names not present "
            f"in the reference: {unexpected}. Partial reconstruction requires an "
            "unambiguous reference-name subset."
        )
    missing_names = [
        atom["name"] for atom in reference_heavy if atom["name"] not in observed_name_set
    ]
    if len(heavy_atoms) < 3:
        raise PipelineError(
            f"cholesterol {key} has fewer than three heavy-atom anchors; cannot fit "
            "the reference safely."
        )

    ref_anchor_xyz = coordinates(
        [reference_by_name[atom["name"]] for atom in heavy_atoms]
    )
    target_anchor_xyz = coordinates(heavy_atoms)
    if (
        np.linalg.matrix_rank(ref_anchor_xyz - ref_anchor_xyz.mean(axis=0)) < 2
        or np.linalg.matrix_rank(target_anchor_xyz - target_anchor_xyz.mean(axis=0)) < 2
    ):
        raise PipelineError(
            f"cholesterol {key} does not provide three non-collinear heavy-atom "
            "anchors for reconstruction."
        )
    rotation, mobile_centroid, target_centroid = kabsch_transform(
        ref_anchor_xyz, target_anchor_xyz
    )
    fitted_anchors = apply_transform(
        ref_anchor_xyz, rotation, mobile_centroid, target_centroid
    )
    fit_rmsd = float(
        np.sqrt(np.mean(np.sum((fitted_anchors - target_anchor_xyz) ** 2, axis=1)))
    )
    if fit_rmsd > max_fit_rmsd:
        raise PipelineError(
            f"cholesterol {key} is missing {missing_names}, but the surviving heavy "
            f"atoms fit the reference at {fit_rmsd:.3f} A RMSD, exceeding "
            f"max_heavy_atom_fit_rmsd={max_fit_rmsd:.3f} A."
        )

    # Confirm the observed subgraph is compatible before adding coordinates.
    ref_bonds = build_bonds(reference_heavy)
    target_bonds = build_bonds(heavy_atoms)
    ref_index = {atom["name"]: i for i, atom in enumerate(reference_heavy)}
    for i, first in enumerate(heavy_atoms):
        for j, second in enumerate(heavy_atoms[i + 1:], start=i + 1):
            ref_edge = ref_index[second["name"]] in ref_bonds[ref_index[first["name"]]]
            target_edge = j in target_bonds[i]
            if ref_edge != target_edge:
                raise PipelineError(
                    f"cholesterol {key} observed heavy-atom connectivity conflicts "
                    f"with the reference between {first['name']} and {second['name']}; "
                    "missing atoms were not reconstructed."
                )

    reconstructed = list(heavy_atoms)
    chain, resname, resid, icode = key
    local_fit_rmsds: dict[str, float] = {}
    for ref_index_value, ref_atom in enumerate(reference_heavy):
        if ref_atom["name"] not in missing_names:
            continue

        # Cholesterol's hydrocarbon tail is flexible, so a whole-molecule
        # transform is suitable for validating identity but not necessarily
        # for placing an absent tail atom. Expand outward through the reference
        # graph until the nearest observed atoms provide a local 3D frame.
        graph_distance = {ref_index_value: 0}
        frontier = [ref_index_value]
        while frontier:
            current = frontier.pop(0)
            for neighbor in ref_bonds[current]:
                if neighbor not in graph_distance:
                    graph_distance[neighbor] = graph_distance[current] + 1
                    frontier.append(neighbor)
        observed_ref_indices = [
            ref_index[name]
            for name in observed_names
            if name in reference_by_name
        ]
        local_anchor_indices = None
        max_graph_distance = max(
            (graph_distance.get(index, 0) for index in observed_ref_indices),
            default=0,
        )
        for radius in range(1, max_graph_distance + 1):
            candidates = [
                index for index in observed_ref_indices
                if graph_distance.get(index, math.inf) <= radius
            ]
            if len(candidates) < 3:
                continue
            local_ref_xyz = coordinates(
                [reference_heavy[index] for index in candidates]
            )
            local_target_xyz = coordinates(
                [observed_by_name[reference_heavy[index]["name"]]
                 for index in candidates]
            )
            if (
                np.linalg.matrix_rank(local_ref_xyz - local_ref_xyz.mean(axis=0)) >= 2
                and np.linalg.matrix_rank(
                    local_target_xyz - local_target_xyz.mean(axis=0)
                ) >= 2
            ):
                local_anchor_indices = candidates
                break
        if local_anchor_indices is None:
            raise PipelineError(
                f"cholesterol {key} has no non-collinear local reference frame "
                f"for missing heavy atom {ref_atom['name']}."
            )

        local_ref_xyz = coordinates(
            [reference_heavy[index] for index in local_anchor_indices]
        )
        local_target_xyz = coordinates(
            [observed_by_name[reference_heavy[index]["name"]]
             for index in local_anchor_indices]
        )
        local_rotation, local_mobile_centroid, local_target_centroid = (
            kabsch_transform(local_ref_xyz, local_target_xyz)
        )
        fitted_local_anchors = apply_transform(
            local_ref_xyz,
            local_rotation,
            local_mobile_centroid,
            local_target_centroid,
        )
        local_fit_rmsd = float(
            np.sqrt(np.mean(np.sum(
                (fitted_local_anchors - local_target_xyz) ** 2, axis=1
            )))
        )
        if local_fit_rmsd > max_fit_rmsd:
            anchor_names = [
                reference_heavy[index]["name"] for index in local_anchor_indices
            ]
            raise PipelineError(
                f"cholesterol {key} local anchors {anchor_names} for missing "
                f"heavy atom {ref_atom['name']} fit at {local_fit_rmsd:.3f} A "
                f"RMSD, exceeding max_heavy_atom_fit_rmsd="
                f"{max_fit_rmsd:.3f} A."
            )
        fitted_missing = apply_transform(
            coordinates([ref_atom]),
            local_rotation,
            local_mobile_centroid,
            local_target_centroid,
        )[0]
        x, y, z = map(float, fitted_missing)
        modeled = dict(ref_atom)
        modeled.update(
            serial=0,
            resname=resname,
            chain=chain,
            resid=resid,
            icode=icode,
            x=x,
            y=y,
            z=z,
            is_h=False,
        )
        modeled["line"] = format_atom_line(
            0, modeled["name"], resname, chain, resid, icode,
            x, y, z, atom_element(modeled),
        )
        reconstructed.append(modeled)
        local_fit_rmsds[ref_atom["name"]] = local_fit_rmsd

    # A guessed coordinate is accepted only if the completed heavy-atom graph
    # exactly matches the reference graph. This catches a plausible-looking
    # coordinate attached to the wrong branch or creating a spurious bond.
    reconstructed_bonds = build_bonds(reconstructed)
    reconstructed_index = {
        atom["name"]: index for index, atom in enumerate(reconstructed)
    }
    graph_conflicts = []
    for first_index, first_atom in enumerate(reference_heavy):
        for second_index in range(first_index + 1, len(reference_heavy)):
            second_atom = reference_heavy[second_index]
            expected_edge = second_index in ref_bonds[first_index]
            observed_edge = (
                reconstructed_index[second_atom["name"]]
                in reconstructed_bonds[reconstructed_index[first_atom["name"]]]
            )
            if expected_edge != observed_edge:
                graph_conflicts.append(
                    f"{first_atom['name']}-{second_atom['name']} "
                    f"({'missing' if expected_edge else 'unexpected'} bond)"
                )
    if graph_conflicts:
        raise PipelineError(
            f"cholesterol {key} reconstructed heavy-atom connectivity does not "
            f"match the reference: {graph_conflicts[:8]}."
        )

    audit.update(
        missing_heavy_atom_names=missing_names,
        modeled_heavy_atom_names=missing_names,
        reference_fit_rmsd_angstrom=fit_rmsd,
        modeled_atom_local_fit_rmsd_angstrom=local_fit_rmsds,
        status="reconstructed",
    )
    return reconstructed, audit


def detect_and_protonate(
    source_path: Path,
    raw_output: Path,
    resnames: set[str],
    obabel: str,
    ref_atoms: Sequence[dict],
    reconstruct_missing_heavy_atoms: bool,
    max_missing_heavy_atoms: int,
    max_heavy_atom_fit_rmsd: float,
) -> tuple[set[str], list[tuple[str, str, int, str]], list[dict]]:
    print("=" * 72)
    print("STEP 1: DETECT AND PROTONATE CHOLESTEROL")
    print("=" * 72)
    records = read_pdb_records(source_path)
    residues = group_selected_residues(records, resnames)
    if not residues:
        counts: dict[str, int] = defaultdict(int)
        for record in records:
            if record["is_atom"]:
                counts[record["resname"]] += 1
        present = ", ".join(f"{name} ({count})" for name, count in sorted(counts.items()))
        raise PipelineError(
            f"no residues named {sorted(resnames)} were found in {source_path}.\n"
            f"Residue names present, with atom counts: {present or 'none'}\n"
            "Use --chol-resnames to supply the cholesterol residue name."
        )

    print(f"Input: {source_path}")
    print(f"Detected {len(residues)} cholesterol residue(s):")
    for key, atoms in sorted(residues.items()):
        n_h = sum(atom["is_h"] for atom in atoms)
        print(f"  {key}: {len(atoms)} atoms, including {n_h} hydrogen(s)")

    protonated: dict[tuple[str, str, int, str], list[dict]] = {}
    heavy_atom_audit: list[dict] = []
    with tempfile.TemporaryDirectory(prefix="charmm36_chol_obabel_") as directory:
        workdir = Path(directory)
        for number, (key, atoms) in enumerate(sorted(residues.items()), start=1):
            heavy = [atom for atom in atoms if not atom["is_h"]]
            stripped = len(atoms) - len(heavy)
            heavy, audit = audit_and_reconstruct_heavy_atoms(
                heavy,
                ref_atoms,
                key,
                reconstruct_missing_heavy_atoms,
                max_missing_heavy_atoms,
                max_heavy_atom_fit_rmsd,
            )
            heavy_atom_audit.append(audit)
            tag = f"chol_{number}_{key[1]}_{key[2]}"
            heavy_out, hydrogens = protonate_heavy_atoms(heavy, tag, workdir, obabel)
            protonated[key] = heavy_out + hydrogens
            print(
                f"  {key}: retained {len(heavy)} heavy atoms, stripped {stripped} "
                f"old H, added {len(hydrogens)} H"
            )
            if audit["modeled_heavy_atom_names"]:
                print(
                    f"      MODELED missing heavy atoms from fitted reference: "
                    f"{audit['modeled_heavy_atom_names']} "
                    f"(whole-molecule anchor RMSD "
                    f"{audit['reference_fit_rmsd_angstrom']:.3f} A; local fits "
                    f"{audit['modeled_atom_local_fit_rmsd_angstrom']})"
                )

    output_lines: list[str] = []
    written: set[tuple[str, str, int, str]] = set()
    for record in records:
        if not record["is_atom"]:
            output_lines.append(record["line"])
            continue
        key = residue_key(record)
        if key not in residues:
            output_lines.append(record["line"])
            continue
        if key in written:
            continue
        written.add(key)
        chain, resname, resid, icode = key
        for atom in protonated[key]:
            output_lines.append(
                format_atom_line(
                    0, atom["name"], resname, chain, resid, icode,
                    atom["x"], atom["y"], atom["z"], atom_element(atom),
                )
            )

    write_renumbered(raw_output, output_lines, strip_conect=True)
    print(f"Wrote raw protonated intermediate: {raw_output}")
    found_names = {key[1].upper() for key in residues}
    return found_names, sorted(residues), heavy_atom_audit


def build_bonds(atoms: Sequence[dict]) -> dict[int, set[int]]:
    bonds: dict[int, set[int]] = defaultdict(set)
    heavy_positions = [i for i, atom in enumerate(atoms) if not atom["is_h"]]
    hydrogen_positions = [i for i, atom in enumerate(atoms) if atom["is_h"]]

    # Heavy-atom connectivity is perceived from the geometry.
    for offset, i in enumerate(heavy_positions):
        for j in heavy_positions[offset + 1:]:
            if distance(atoms[i], atoms[j]) <= HEAVY_BOND_CUTOFF:
                bonds[i].add(j)
                bonds[j].add(i)

    # Assign every hydrogen to at most one heavy atom. Using every pair under
    # the cutoff can spuriously attach a hydrogen to two adjacent carbons.
    for hydrogen_i in hydrogen_positions:
        candidates = [
            (distance(atoms[hydrogen_i], atoms[heavy_i]), heavy_i)
            for heavy_i in heavy_positions
            if distance(atoms[hydrogen_i], atoms[heavy_i]) <= H_BOND_CUTOFF
        ]
        if candidates:
            _bond_length, heavy_i = min(candidates)
            bonds[hydrogen_i].add(heavy_i)
            bonds[heavy_i].add(hydrogen_i)
    return bonds


def atom_signature(
    index: int, atoms: Sequence[dict], bonds: dict[int, set[int]]
) -> tuple[str, bool, int, int]:
    neighbors = [atoms[j] for j in bonds[index]]
    atom = atoms[index]
    return (
        atom_element(atom),
        atom["is_h"],
        sum(neighbor["is_h"] for neighbor in neighbors),
        sum(not neighbor["is_h"] for neighbor in neighbors),
    )


def find_mapping(ref_atoms: Sequence[dict], target_atoms: Sequence[dict]) -> dict[int, int]:
    if len(ref_atoms) != len(target_atoms):
        raise PipelineError(
            f"cannot map residues with different atom counts: reference "
            f"{len(ref_atoms)}, target {len(target_atoms)}"
        )
    ref_bonds = build_bonds(ref_atoms)
    target_bonds = build_bonds(target_atoms)
    ref_signatures = {
        i: atom_signature(i, ref_atoms, ref_bonds) for i in range(len(ref_atoms))
    }
    target_signatures = {
        i: atom_signature(i, target_atoms, target_bonds)
        for i in range(len(target_atoms))
    }
    candidates: dict[int, list[int]] = {}
    for target_i in range(len(target_atoms)):
        choices = [
            ref_i
            for ref_i in range(len(ref_atoms))
            if ref_signatures[ref_i] == target_signatures[target_i]
        ]
        choices.sort(
            key=lambda ref_i: (
                target_atoms[target_i]["name"] != ref_atoms[ref_i]["name"],
                ref_i,
            )
        )
        candidates[target_i] = choices
    impossible = [i for i, choices in candidates.items() if not choices]
    if impossible:
        names = [target_atoms[i]["name"] for i in impossible]
        raise PipelineError(
            f"no reference mapping candidates for target atoms: {names}"
        )

    target_order = sorted(candidates, key=lambda i: (len(candidates[i]), i))
    mapping: dict[int, int] = {}
    used_ref: set[int] = set()

    def compatible(target_i: int, ref_i: int) -> bool:
        return all(
            ((other_target in target_bonds[target_i])
             == (other_ref in ref_bonds[ref_i]))
            for other_target, other_ref in mapping.items()
        )

    def backtrack(position: int) -> bool:
        if position == len(target_order):
            return True
        target_i = target_order[position]
        for ref_i in candidates[target_i]:
            if ref_i in used_ref or not compatible(target_i, ref_i):
                continue
            mapping[target_i] = ref_i
            used_ref.add(ref_i)
            if backtrack(position + 1):
                return True
            used_ref.remove(ref_i)
            del mapping[target_i]
        return False

    if not backtrack(0):
        raise PipelineError("could not map target cholesterol topology to reference")
    return mapping


def bonded_hydrogens(atoms: Sequence[dict]) -> dict[int, list[int]]:
    result: dict[int, list[int]] = defaultdict(list)
    bonds = build_bonds(atoms)
    for heavy_i, heavy in enumerate(atoms):
        if not heavy["is_h"]:
            result[heavy_i] = sorted(
                (i for i in bonds[heavy_i] if atoms[i]["is_h"]),
                key=lambda i: distance(heavy, atoms[i]),
            )
    return result


def prune_extra_hydrogens(
    target_atoms: Sequence[dict], ref_atoms: Sequence[dict]
) -> tuple[list[dict], set[int]]:
    target_heavy = [atom for atom in target_atoms if not atom["is_h"]]
    ref_heavy = [atom for atom in ref_atoms if not atom["is_h"]]
    heavy_mapping = find_mapping(ref_heavy, target_heavy)

    ref_h_connections = bonded_hydrogens(ref_atoms)
    target_h_connections = bonded_hydrogens(target_atoms)
    ref_heavy_positions = [i for i, atom in enumerate(ref_atoms) if not atom["is_h"]]
    target_heavy_positions = [i for i, atom in enumerate(target_atoms) if not atom["is_h"]]
    remove_positions: set[int] = set()

    for target_heavy_i, ref_heavy_i in heavy_mapping.items():
        target_full_i = target_heavy_positions[target_heavy_i]
        ref_full_i = ref_heavy_positions[ref_heavy_i]
        expected = len(ref_h_connections.get(ref_full_i, []))
        observed_positions = target_h_connections.get(target_full_i, [])
        if len(observed_positions) > expected:
            remove_positions.update(observed_positions[expected:])

    pruned = [atom for i, atom in enumerate(target_atoms) if i not in remove_positions]
    return pruned, remove_positions


def format_atom_name(name: str, element: str) -> str:
    if len(name) == 4:
        return name
    if len(element.strip()) == 1:
        return f" {name:<3}"
    return f"{name:<4}"


def replace_atom_name(line: str, new_name: str, element: str) -> str:
    return line[:12] + format_atom_name(new_name, element) + line[16:]


def replace_xyz(line: str, xyz: Sequence[float]) -> str:
    return line[:30] + f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}" + line[54:]


def replace_residue_info(line: str, target_atom: dict) -> str:
    result = line[:17] + f"{target_atom['resname']:<4}"[:4] + line[21:]
    result = result[:21] + f"{target_atom['chain'][:1]:1}" + result[22:]
    result = result[:22] + f"{target_atom['resid']:4d}" + result[26:]
    result = result[:26] + f"{target_atom['icode'][:1]:1}" + result[27:]
    return result


def update_serial(line: str, serial: int) -> str:
    return line[:6] + f"{serial:5d}" + line[11:]


def write_renumbered(path: Path, lines: Iterable[str], strip_conect: bool) -> None:
    serial = 1
    with path.open("w") as handle:
        for line in lines:
            if strip_conect and line.startswith("CONECT"):
                continue
            if line.startswith(("ATOM  ", "HETATM")):
                line = update_serial(line, serial)
                serial += 1
            handle.write(line)


def coordinates(atoms: Sequence[dict]) -> np.ndarray:
    return np.array([[a["x"], a["y"], a["z"]] for a in atoms], dtype=float)


def kabsch_transform(
    mobile: np.ndarray, target: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mobile_centroid = mobile.mean(axis=0)
    target_centroid = target.mean(axis=0)
    covariance = (mobile - mobile_centroid).T @ (target - target_centroid)
    u, _singular_values, vt = np.linalg.svd(covariance)
    rotation = vt.T @ u.T
    if np.linalg.det(rotation) < 0:
        vt[-1, :] *= -1
        rotation = vt.T @ u.T
    return rotation, mobile_centroid, target_centroid


def apply_transform(
    xyz: np.ndarray,
    rotation: np.ndarray,
    mobile_centroid: np.ndarray,
    target_centroid: np.ndarray,
) -> np.ndarray:
    return (xyz - mobile_centroid) @ rotation.T + target_centroid


def reference_to_target_transform(
    ref_atoms: Sequence[dict], ordered_target_atoms: Sequence[dict]
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Align the reference to a mapped target using heavy atoms only.

    Equivalent hydrogens can initially be mapped in either order. Excluding
    them from the Kabsch fit prevents those arbitrary assignments from
    influencing the reference orientation used to resolve their labels.
    """
    heavy_positions = [i for i, atom in enumerate(ref_atoms) if not atom["is_h"]]
    if len(heavy_positions) < 3:
        raise PipelineError("at least three reference heavy atoms are required for alignment")
    ref_heavy_xyz = coordinates([ref_atoms[i] for i in heavy_positions])
    target_heavy_xyz = coordinates([ordered_target_atoms[i] for i in heavy_positions])
    return kabsch_transform(ref_heavy_xyz, target_heavy_xyz)


def atom_to_xyz_distance(atom: dict, xyz: Sequence[float]) -> float:
    return math.sqrt(
        (atom["x"] - xyz[0]) ** 2
        + (atom["y"] - xyz[1]) ** 2
        + (atom["z"] - xyz[2]) ** 2
    )


def disambiguate_equivalent_hydrogens(
    ref_atoms: Sequence[dict],
    ref_bonds: dict[int, set[int]],
    ordered_target_atoms: Sequence[dict],
    placed_ref_xyz: np.ndarray,
) -> tuple[list[dict], list[tuple[str, tuple[str, ...], float]]]:
    """Resolve topologically equivalent hydrogen labels using geometry.

    A topology mapping cannot distinguish hydrogens attached to the same
    heavy atom. For each CH2 or CH3 group, this function tries every possible
    assignment of the mapped target hydrogen coordinates to the reference
    hydrogen-name slots and retains the assignment with the lowest total
    distance to the aligned reference positions. Coordinates are not moved;
    only the atom-to-reference-name assignment is changed.
    """
    corrected = list(ordered_target_atoms)
    changes: list[tuple[str, tuple[str, ...], float]] = []

    for heavy_i, heavy in enumerate(ref_atoms):
        if heavy["is_h"]:
            continue
        hydrogen_slots = sorted(
            i for i in ref_bonds[heavy_i] if ref_atoms[i]["is_h"]
        )
        if len(hydrogen_slots) < 2:
            continue

        current_atoms = [corrected[i] for i in hydrogen_slots]
        identity = tuple(range(len(hydrogen_slots)))

        def assignment_cost(permutation: tuple[int, ...]) -> float:
            return sum(
                atom_to_xyz_distance(
                    current_atoms[permutation[position]],
                    placed_ref_xyz[hydrogen_slots[position]],
                )
                for position in range(len(hydrogen_slots))
            )

        identity_cost = assignment_cost(identity)
        best_permutation = identity
        best_cost = identity_cost
        for candidate in permutations(identity):
            candidate_cost = assignment_cost(candidate)
            if candidate_cost < best_cost - 1.0e-9:
                best_permutation = candidate
                best_cost = candidate_cost

        if best_permutation != identity:
            for position, ref_i in enumerate(hydrogen_slots):
                corrected[ref_i] = current_atoms[best_permutation[position]]
            changes.append(
                (
                    heavy["name"],
                    tuple(ref_atoms[i]["name"] for i in hydrogen_slots),
                    identity_cost - best_cost,
                )
            )

    return corrected, changes


def write_replaced_system(
    output_path: Path,
    target_records: Sequence[dict],
    residues: dict[tuple[str, str, int, str], list[dict]],
    replacement_lines: dict[tuple[str, str, int, str], list[str]],
) -> None:
    lines: list[str] = []
    written: set[tuple[str, str, int, str]] = set()
    for record in target_records:
        if not record["is_atom"]:
            lines.append(record["line"])
            continue
        key = residue_key(record)
        if key not in residues:
            lines.append(record["line"])
            continue
        if key in written:
            continue
        lines.extend(replacement_lines[key])
        written.add(key)
    write_renumbered(output_path, lines, strip_conect=True)


def build_outputs(
    reference_path: Path,
    raw_path: Path,
    fixed_path: Path,
    aligned_reference_path: Path,
    reference_placed_path: Path,
    found_resnames: set[str],
) -> None:
    print("\n" + "=" * 72)
    print("STEP 2: MAP TO CHARMM36 REFERENCE AND BUILD OUTPUTS")
    print("=" * 72)
    ref_records = read_pdb_records(reference_path)
    ref_atoms = reference_atoms(ref_records, reference_path)
    target_records = read_pdb_records(raw_path)
    residues = group_selected_residues(target_records, found_resnames)
    if not residues:
        raise PipelineError(f"no cholesterol residues found in intermediate {raw_path}")

    ref_xyz = coordinates(ref_atoms)
    ref_bonds = build_bonds(ref_atoms)
    ref_lines = [atom["line"] for atom in ref_atoms]
    fixed_lines_by_residue: dict[tuple[str, str, int, str], list[str]] = {}
    placed_lines_by_residue: dict[tuple[str, str, int, str], list[str]] = {}
    first_key = sorted(residues)[0]
    first_ordered_atoms: list[dict] | None = None
    total_removed = 0
    total_orientation_corrections = 0

    for key, atoms in sorted(residues.items()):
        pruned, removed_positions = prune_extra_hydrogens(atoms, ref_atoms)
        total_removed += len(removed_positions)
        if len(pruned) != len(ref_atoms):
            raise PipelineError(
                f"after removing excess hydrogens, cholesterol {key} has "
                f"{len(pruned)} atoms but reference has {len(ref_atoms)}. "
                "Open Babel may have under-protonated the residue or the input "
                "may not have cholesterol connectivity."
            )
        mapping = find_mapping(ref_atoms, pruned)
        ref_to_target = {ref_i: pruned[target_i] for target_i, ref_i in mapping.items()}
        ordered = [ref_to_target[i] for i in range(len(ref_atoms))]

        rotation, mobile_centroid, target_centroid = reference_to_target_transform(
            ref_atoms, ordered
        )
        placed_xyz = apply_transform(
            ref_xyz, rotation, mobile_centroid, target_centroid
        )
        ordered, orientation_changes = disambiguate_equivalent_hydrogens(
            ref_atoms, ref_bonds, ordered, placed_xyz
        )
        total_orientation_corrections += len(orientation_changes)
        for heavy_name, hydrogen_names, improvement in orientation_changes:
            print(
                f"  {key}: reassigned {', '.join(hydrogen_names)} on "
                f"{heavy_name} using reference geometry "
                f"(distance improved by {improvement:.3f} A)"
            )

        fixed_lines_by_residue[key] = [
            replace_atom_name(target["line"], ref["name"], atom_element(target))
            for ref, target in zip(ref_atoms, ordered)
        ]
        if key == first_key:
            first_ordered_atoms = ordered
        placed_lines_by_residue[key] = [
            replace_xyz(replace_residue_info(line, ordered[i]), placed_xyz[i])
            for i, line in enumerate(ref_lines)
        ]
        print(
            f"  {key}: removed {len(removed_positions)} excess H; "
            "renamed, reordered, and aligned"
        )

    write_replaced_system(
        fixed_path, target_records, residues, fixed_lines_by_residue
    )
    write_replaced_system(
        reference_placed_path, target_records, residues, placed_lines_by_residue
    )

    assert first_ordered_atoms is not None
    rotation, mobile_centroid, target_centroid = reference_to_target_transform(
        ref_atoms, first_ordered_atoms
    )
    aligned_xyz = apply_transform(
        ref_xyz, rotation, mobile_centroid, target_centroid
    )
    aligned_lines: list[str] = []
    atom_index = 0
    for record in ref_records:
        line = record["line"]
        if record["is_atom"]:
            line = replace_xyz(line, aligned_xyz[atom_index])
            atom_index += 1
        aligned_lines.append(line)
    write_renumbered(aligned_reference_path, aligned_lines, strip_conect=True)

    print(f"Removed {total_removed} excess hydrogen(s) in total")
    print(
        "Corrected equivalent-hydrogen orientation for "
        f"{total_orientation_corrections} heavy-atom group(s)"
    )
    print(f"Wrote main CHARMM36-formatted PDB: {fixed_path}")
    print(f"Wrote reference aligned to first cholesterol: {aligned_reference_path}")
    print(f"Wrote full system with placed references: {reference_placed_path}")
    print(f"First cholesterol used for the single-reference alignment: {first_key}")


def guess_element(atom_name: str, element_field: str) -> str | None:
    element = element_field.strip().upper()
    if element in COVALENT_RADII:
        return element
    name = "".join(char for char in atom_name.strip() if not char.isdigit())
    name = name.replace("'", "")
    if len(name) >= 2 and name[:2].upper() in COVALENT_RADII:
        return name[:2].upper()
    if name and name[0].upper() in COVALENT_RADII:
        return name[0].upper()
    return None


class VerificationAtom:
    __slots__ = ("name", "resname", "chain", "resid", "icode", "x", "y", "z", "element")

    def __init__(
        self, name: str, resname: str, chain: str, resid: int, icode: str,
        x: float, y: float, z: float, element: str | None,
    ) -> None:
        self.name = name
        self.resname = resname
        self.chain = chain
        self.resid = resid
        self.icode = icode
        self.x = x
        self.y = y
        self.z = z
        self.element = element


def verification_parse(path: Path) -> list[VerificationAtom]:
    atoms: list[VerificationAtom] = []
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            try:
                name = line[12:16].strip()
                element_field = line[76:78] if len(line) >= 78 else ""
                atoms.append(
                    VerificationAtom(
                        name=name,
                        resname=line[17:21].strip(),
                        chain=line[21].strip(),
                        resid=int(line[22:26]),
                        icode=line[26].strip(),
                        x=float(line[30:38]),
                        y=float(line[38:46]),
                        z=float(line[46:54]),
                        element=guess_element(name, element_field),
                    )
                )
            except (ValueError, IndexError) as exc:
                raise PipelineError(
                    f"independent verifier could not parse {path}:{line_number}"
                ) from exc
    return atoms


def group_verification_atoms(
    atoms: Iterable[VerificationAtom], resnames: set[str] | None = None
) -> dict[tuple[str, str, int, str], list[VerificationAtom]]:
    groups: dict[tuple[str, str, int, str], list[VerificationAtom]] = defaultdict(list)
    for atom in atoms:
        if resnames is None or atom.resname.upper() in resnames:
            groups[(atom.chain, atom.resname, atom.resid, atom.icode)].append(atom)
    return groups


def verification_bonds(
    atoms: Sequence[VerificationAtom], tolerance: float
) -> set[frozenset[str]]:
    edges: set[frozenset[str]] = set()
    heavy_positions = [i for i, atom in enumerate(atoms) if atom.element != "H"]
    hydrogen_positions = [i for i, atom in enumerate(atoms) if atom.element == "H"]

    for offset, i in enumerate(heavy_positions):
        atom = atoms[i]
        radius = COVALENT_RADII.get(atom.element)
        for j in heavy_positions[offset + 1:]:
            other = atoms[j]
            other_radius = COVALENT_RADII.get(other.element)
            if radius is None or other_radius is None:
                continue
            squared_distance = (
                (atom.x - other.x) ** 2
                + (atom.y - other.y) ** 2
                + (atom.z - other.z) ** 2
            )
            cutoff = radius + other_radius + tolerance
            if squared_distance <= cutoff * cutoff:
                edges.add(frozenset((atom.name, other.name)))

    for hydrogen_i in hydrogen_positions:
        hydrogen = atoms[hydrogen_i]
        hydrogen_radius = COVALENT_RADII["H"]
        candidates: list[tuple[float, int]] = []
        for heavy_i in heavy_positions:
            heavy = atoms[heavy_i]
            heavy_radius = COVALENT_RADII.get(heavy.element)
            if heavy_radius is None:
                continue
            squared_distance = (
                (hydrogen.x - heavy.x) ** 2
                + (hydrogen.y - heavy.y) ** 2
                + (hydrogen.z - heavy.z) ** 2
            )
            cutoff = hydrogen_radius + heavy_radius + tolerance
            if squared_distance <= cutoff * cutoff:
                candidates.append((squared_distance, heavy_i))
        if candidates:
            _squared_distance, heavy_i = min(candidates)
            edges.add(frozenset((hydrogen.name, atoms[heavy_i].name)))
    return edges


def edge_string(edge: frozenset[str]) -> str:
    return "-".join(sorted(edge))


def verify_residue(
    ref_atoms: Sequence[VerificationAtom],
    target_atoms: Sequence[VerificationAtom],
    tolerance: float,
) -> list[str]:
    ref_names = [atom.name for atom in ref_atoms]
    target_names = [atom.name for atom in target_atoms]
    ref_set, target_set = set(ref_names), set(target_names)
    problems: list[str] = []
    if ref_set != target_set:
        problems.append(
            f"atom-name mismatch; only in reference: {sorted(ref_set - target_set)}, "
            f"only in target: {sorted(target_set - ref_set)}"
        )
    if ref_names != target_names:
        limit = min(len(ref_names), len(target_names))
        first_difference = next(
            (i for i in range(limit) if ref_names[i] != target_names[i]), limit
        )
        problems.append(
            f"atom-order mismatch; first difference at position {first_difference + 1}"
        )

    common = ref_set & target_set
    ref_edges = {
        edge for edge in verification_bonds(ref_atoms, tolerance) if set(edge) <= common
    }
    target_edges = {
        edge
        for edge in verification_bonds(target_atoms, tolerance)
        if set(edge) <= common
    }
    if ref_edges != target_edges:
        problems.append(
            "connectivity mismatch; only in reference: "
            f"{sorted(edge_string(edge) for edge in ref_edges - target_edges)}, "
            "only in target: "
            f"{sorted(edge_string(edge) for edge in target_edges - ref_edges)}"
        )
    return problems


def verification_distance(
    first: VerificationAtom, second: VerificationAtom
) -> float:
    return math.sqrt(
        (first.x - second.x) ** 2
        + (first.y - second.y) ** 2
        + (first.z - second.z) ** 2
    )


def verify_equivalent_hydrogen_orientation(
    placed_reference_atoms: Sequence[VerificationAtom],
    target_atoms: Sequence[VerificationAtom],
    tolerance: float,
) -> list[str]:
    """Check that equivalent H labels give the closest reference assignment."""
    placed_by_name = {atom.name: atom for atom in placed_reference_atoms}
    target_by_name = {atom.name: atom for atom in target_atoms}
    placed_edges = verification_bonds(placed_reference_atoms, tolerance)
    heavy_to_hydrogens: dict[str, set[str]] = defaultdict(set)

    for edge in placed_edges:
        first_name, second_name = tuple(edge)
        first = placed_by_name.get(first_name)
        second = placed_by_name.get(second_name)
        if first is None or second is None:
            continue
        if first.element == "H" and second.element != "H":
            heavy_to_hydrogens[second.name].add(first.name)
        elif second.element == "H" and first.element != "H":
            heavy_to_hydrogens[first.name].add(second.name)

    problems: list[str] = []
    reference_name_order = [atom.name for atom in placed_reference_atoms]
    for heavy_name, hydrogen_name_set in sorted(heavy_to_hydrogens.items()):
        if len(hydrogen_name_set) < 2:
            continue
        hydrogen_names = [
            name for name in reference_name_order if name in hydrogen_name_set
        ]
        missing = [name for name in hydrogen_names if name not in target_by_name]
        if missing:
            problems.append(
                f"cannot verify equivalent hydrogens on {heavy_name}; "
                f"target is missing {missing}"
            )
            continue

        current_atoms = [target_by_name[name] for name in hydrogen_names]
        identity = tuple(range(len(hydrogen_names)))

        def assignment_cost(permutation: tuple[int, ...]) -> float:
            return sum(
                verification_distance(
                    current_atoms[permutation[position]],
                    placed_by_name[hydrogen_names[position]],
                )
                for position in range(len(hydrogen_names))
            )

        identity_cost = assignment_cost(identity)
        best_permutation = identity
        best_cost = identity_cost
        for candidate in permutations(identity):
            candidate_cost = assignment_cost(candidate)
            if candidate_cost < best_cost - 1.0e-4:
                best_permutation = candidate
                best_cost = candidate_cost

        if best_permutation != identity:
            reassignment = ", ".join(
                f"{hydrogen_names[position]}<-{hydrogen_names[best_permutation[position]]}"
                for position in range(len(hydrogen_names))
            )
            problems.append(
                f"equivalent-hydrogen orientation mismatch on {heavy_name} "
                f"({reassignment}; distance would improve by "
                f"{identity_cost - best_cost:.3f} A)"
            )

    return problems


def verify_output(
    reference_path: Path,
    fixed_path: Path,
    reference_placed_path: Path,
    found_resnames: set[str],
    tolerance: float,
) -> None:
    print("\n" + "=" * 72)
    print("STEP 3: INDEPENDENT VERIFICATION")
    print("=" * 72)
    ref_groups = group_verification_atoms(verification_parse(reference_path))
    if len(ref_groups) != 1:
        raise PipelineError(
            f"verifier expected one reference residue; found {len(ref_groups)}"
        )
    ref_atoms = next(iter(ref_groups.values()))
    target_groups = group_verification_atoms(
        verification_parse(fixed_path), found_resnames
    )
    if not target_groups:
        raise PipelineError(f"verifier found no cholesterol in {fixed_path}")
    placed_groups = group_verification_atoms(
        verification_parse(reference_placed_path), found_resnames
    )
    if not placed_groups:
        raise PipelineError(
            f"verifier found no placed reference cholesterol in {reference_placed_path}"
        )

    passed = 0
    failed = 0
    for key, atoms in sorted(target_groups.items()):
        problems = verify_residue(ref_atoms, atoms, tolerance)
        if key not in placed_groups:
            problems.append(f"no matching placed reference cholesterol found for {key}")
        else:
            problems.extend(
                verify_equivalent_hydrogen_orientation(
                    placed_groups[key], atoms, tolerance
                )
            )
        if problems:
            failed += 1
            print(f"[FAIL] {key} ({len(atoms)} atoms)")
            for problem in problems:
                print(f"  {problem}")
        else:
            passed += 1
            print(
                f"[PASS] {key} ({len(atoms)} atoms): names, order, and "
                "connectivity match the reference; equivalent-hydrogen "
                "labels match the aligned reference geometry"
            )
    print(f"Summary: {passed} PASS, {failed} FAIL, {passed + failed} checked")
    if failed:
        raise PipelineError(f"verification failed for {failed} cholesterol residue(s)")


def print_configuration(args: argparse.Namespace, reference_label: str) -> None:
    print("Configuration")
    print(f"  Input PDB:                 {args.input}")
    print(f"  Reference:                 {reference_label}")
    print(f"  Raw intermediate:          {args.raw_out}")
    print(f"  Main fixed output:         {args.fixed_out}")
    print(f"  Aligned reference output:  {args.aligned_reference_out}")
    print(f"  Reference-placed output:   {args.reference_placed_out}")
    print(f"  Cholesterol residue names: {', '.join(sorted(args.resnames))}")
    print(
        "  Missing-heavy reconstruction: "
        f"{args.reconstruct_missing_heavy_atoms} "
        f"(max {args.max_missing_heavy_atoms}, fit RMSD <= "
        f"{args.max_heavy_atom_fit_rmsd:.3f} A)"
    )
    print(f"  Open Babel:                {args.obabel_resolved}")
    print()


def run(args: argparse.Namespace) -> dict:
    derive_output_paths(args)
    args.resnames = parse_resnames(args.chol_resnames)
    validate_paths(args)
    args.obabel_resolved = check_obabel(args.obabel)

    with tempfile.TemporaryDirectory(prefix="charmm36_chol_reference_") as directory:
        if args.reference:
            reference_path = args.reference
            reference_label = str(reference_path)
        else:
            reference_path = Path(directory) / "builtin_charmm36_cholesterol.pdb"
            write_embedded_reference(reference_path)
            reference_label = "built-in CHARMM36 cholesterol (74 atoms)"

        ref_atoms = reference_atoms(read_pdb_records(reference_path), reference_path)
        if len(ref_atoms) != 74:
            raise PipelineError(
                f"reference contains {len(ref_atoms)} atoms; CHARMM36 cholesterol "
                "is expected to contain 74"
            )
        print_configuration(args, reference_label)
        found_resnames, _residue_keys, heavy_atom_audit = detect_and_protonate(
            args.input,
            args.raw_out,
            args.resnames,
            args.obabel_resolved,
            ref_atoms,
            args.reconstruct_missing_heavy_atoms,
            args.max_missing_heavy_atoms,
            args.max_heavy_atom_fit_rmsd,
        )
        build_outputs(
            reference_path,
            args.raw_out,
            args.fixed_out,
            args.aligned_reference_out,
            args.reference_placed_out,
            found_resnames,
        )
        verify_output(
            reference_path,
            args.fixed_out,
            args.reference_placed_out,
            found_resnames,
            args.verify_tolerance,
        )
        return {"heavy_atom_audit": heavy_atom_audit}


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        run(args)
    except (PipelineError, OSError) as exc:
        parser.exit(1, f"ERROR: {exc}\n")
    print("\nConversion and verification completed successfully.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
