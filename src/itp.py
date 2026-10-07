"""Strict, shared parser for GROMACS ITP molecule definitions.

Every pipeline feature that needs moleculetype names, residue names, atom
counts, atom names, or charge sums uses this module.  Keeping one parser
prevents topology assembly, name restoration, and charge accounting from
silently interpreting the same ITP differently.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from typing import Dict, Tuple


class ITPParseError(RuntimeError):
    pass


_SECTION_RE = re.compile(r"^\s*\[\s*([^]]+)\s*\]")
_CONDITIONAL_DIRECTIVES = {
    "#if", "#ifdef", "#ifndef", "#elif", "#else", "#endif",
}


@dataclass(frozen=True)
class ITPAtom:
    nr: int
    atom_type: str
    resnr: int
    resname: str
    atom_name: str
    cgnr: int
    charge: Decimal


@dataclass(frozen=True)
class ITPMoleculeType:
    name: str
    atoms: Tuple[ITPAtom, ...]
    source_path: str

    @property
    def atom_count(self) -> int:
        return len(self.atoms)

    @property
    def charge(self) -> float:
        return float(sum((atom.charge for atom in self.atoms), Decimal("0")))

    @property
    def resnames(self):
        return {atom.resname for atom in self.atoms}

    @property
    def atom_names(self) -> Tuple[str, ...]:
        return tuple(atom.atom_name for atom in self.atoms)

    def content_signature(self):
        """Definition identity independent of source filename and comments."""
        return tuple(
            (
                atom.nr, atom.atom_type, atom.resnr, atom.resname,
                atom.atom_name, atom.cgnr, atom.charge,
            )
            for atom in self.atoms
        )


def _clean(raw: str) -> str:
    return raw.split(";", 1)[0].strip()


def parse_itp(path: str) -> Dict[str, ITPMoleculeType]:
    """Parse every populated moleculetype in one ITP.

    Parameter-only ITP files return an empty mapping. Preprocessor directives
    outside ``[ moleculetype ]`` and ``[ atoms ]`` are ignored. Conditional
    directives inside either of those sections are rejected because choosing
    an active branch without GROMACS preprocessor definitions would be a guess.
    """
    section = ""
    current_name = None
    need_name = False
    atoms_by_name: Dict[str, list] = {}

    try:
        handle = open(path, encoding="utf-8", errors="strict")
    except OSError as exc:
        raise ITPParseError(f"Cannot read ITP '{path}': {exc}") from exc

    try:
        with handle:
            for line_number, raw in enumerate(handle, 1):
                line = _clean(raw)
                if not line:
                    continue

                if line.startswith("#"):
                    directive = line.split(None, 1)[0].lower()
                    if section in {"moleculetype", "atoms"}:
                        qualifier = (
                            "conditional preprocessor directive"
                            if directive in _CONDITIONAL_DIRECTIVES
                            else "preprocessor directive"
                        )
                        raise ITPParseError(
                            f"{path}:{line_number}: {qualifier} "
                            f"'{directive}' occurs inside [ {section} ]. Supply an ITP "
                            "whose [ atoms ] definition is fully expanded, or preprocess "
                            "the topology with the active GROMACS defines. Silently "
                            "skipping a directive here would change atom counts or charge."
                        )
                    continue

                match = _SECTION_RE.match(line)
                if match:
                    section = match.group(1).strip().lower()
                    if section == "moleculetype":
                        current_name = None
                        need_name = True
                    continue

                fields = line.split()
                if section == "moleculetype" and need_name:
                    current_name = fields[0]
                    if current_name in atoms_by_name:
                        raise ITPParseError(
                            f"{path}:{line_number}: duplicate moleculetype "
                            f"'{current_name}' in the same ITP"
                        )
                    atoms_by_name[current_name] = []
                    need_name = False
                    continue

                if section != "atoms":
                    continue
                if current_name is None:
                    raise ITPParseError(
                        f"{path}:{line_number}: found an [ atoms ] row before a "
                        "[ moleculetype ] name"
                    )
                if len(fields) < 7:
                    raise ITPParseError(
                        f"{path}:{line_number}: [ atoms ] row has fewer than seven "
                        f"columns: '{line}'"
                    )
                try:
                    atom = ITPAtom(
                        nr=int(fields[0]),
                        atom_type=fields[1],
                        resnr=int(fields[2]),
                        resname=fields[3],
                        atom_name=fields[4],
                        cgnr=int(fields[5]),
                        charge=Decimal(fields[6]),
                    )
                except (ValueError, InvalidOperation) as exc:
                    raise ITPParseError(
                        f"{path}:{line_number}: cannot parse [ atoms ] row: '{line}'"
                    ) from exc
                atoms_by_name[current_name].append(atom)
    except UnicodeDecodeError as exc:
        raise ITPParseError(f"ITP '{path}' is not valid UTF-8: {exc}") from exc

    result = {}
    for name, atoms in atoms_by_name.items():
        if not atoms:
            continue
        expected = list(range(1, len(atoms) + 1))
        observed = [atom.nr for atom in atoms]
        if observed != expected:
            raise ITPParseError(
                f"ITP '{path}' moleculetype '{name}' has non-sequential [ atoms ] "
                f"numbers; expected 1..{len(atoms)}"
            )
        result[name] = ITPMoleculeType(name, tuple(atoms), path)
    return result


def molecule_charges(path: str) -> Dict[str, float]:
    return {name: definition.charge for name, definition in parse_itp(path).items()}


def molecule_atom_counts(path: str) -> Dict[str, int]:
    return {name: definition.atom_count for name, definition in parse_itp(path).items()}


def molecule_resnames(path: str):
    return {name: definition.resnames for name, definition in parse_itp(path).items()}
