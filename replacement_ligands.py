"""Resolve charges and topology files for ligands inside a receptor-replacement mask."""

from __future__ import annotations

from dataclasses import dataclass
import re
from typing import Dict, List

from itp import ITPParseError, molecule_charges, parse_itp


class ReplacementLigandError(RuntimeError):
    pass


@dataclass
class ResolvedReplacementLigand:
    selector: str
    resname: str
    occurrence: int | None
    moleculetype: str
    net_charge: float
    itp_charge: float
    itp_path: str
    charge_source: str


def read_itp_molecule_charges(path: str) -> Dict[str, float]:
    """Return the summed [ atoms ] charge for every moleculetype in an ITP."""
    try:
        populated = molecule_charges(path)
    except ITPParseError as exc:
        raise ReplacementLigandError(str(exc)) from exc
    if not populated:
        raise ReplacementLigandError(
            f"No populated [ atoms ] section was found in ligand ITP '{path}'"
        )
    return populated


def _parse_selector(selector: str):
    match = re.fullmatch(r"([A-Za-z0-9_+\-]+)(?::([1-9]\d*))?", selector.strip())
    if not match:
        raise ReplacementLigandError(
            f"Replacement ligand selector '{selector}' must use RES or RES:N syntax, "
            "for example LDP or LDP:2."
        )
    name, occurrence = match.groups()
    return name, int(occurrence) if occurrence else None


def _select_moleculetype(spec, resname: str, definitions) -> str:
    if spec.moleculetype:
        if spec.moleculetype not in definitions:
            raise ReplacementLigandError(
                f"Ligand {spec.resname} requests moleculetype '{spec.moleculetype}', but "
                f"{spec.itp_path} contains {sorted(definitions)}"
            )
        if resname not in definitions[spec.moleculetype].resnames:
            raise ReplacementLigandError(
                f"Ligand {spec.resname} requests moleculetype '{spec.moleculetype}', "
                f"but its [ atoms ] rows use residue name(s) "
                f"{sorted(definitions[spec.moleculetype].resnames)}."
            )
        return spec.moleculetype
    candidates = sorted(
        name for name, definition in definitions.items()
        if resname in definition.resnames
    )
    if len(candidates) == 1:
        return candidates[0]
    if len(candidates) > 1:
        raise ReplacementLigandError(
            f"Ligand {spec.resname}: multiple moleculetype definitions in "
            f"{spec.itp_path} use that residue name: {candidates}. Set "
            "replacement_ligands.moleculetype explicitly."
        )
    raise ReplacementLigandError(
        f"Ligand {spec.resname}: no moleculetype in {spec.itp_path} has [ atoms ] "
        f"rows with that residue name (definitions: {sorted(definitions)})."
    )


def resolve_replacement_ligands(cfg, tolerance: float = 0.01) -> List[ResolvedReplacementLigand]:
    """Populate charge overrides and topology inputs from replacement_ligands."""
    resolved = []
    for spec in cfg.replacement_ligands:
        resname, occurrence = _parse_selector(spec.resname)
        try:
            definitions = parse_itp(spec.itp_path)
        except ITPParseError as exc:
            raise ReplacementLigandError(str(exc)) from exc
        if not definitions:
            raise ReplacementLigandError(
                f"No populated [ atoms ] section was found in ligand ITP '{spec.itp_path}'"
            )
        moleculetype = _select_moleculetype(spec, resname, definitions)
        itp_charge = definitions[moleculetype].charge

        if isinstance(spec.charge, str) and spec.charge.lower() == "from_itp":
            net_charge = itp_charge
            source = "ITP [ atoms ]"
        else:
            net_charge = float(spec.charge)
            source = "YAML"
            if abs(net_charge - itp_charge) > tolerance:
                raise ReplacementLigandError(
                    f"Ligand {spec.resname}: YAML charge {net_charge:+.6f} does not match "
                    f"the charge sum {itp_charge:+.6f} in {spec.itp_path} "
                    f"(tolerance {tolerance})."
                )

        existing = cfg.charge.charge_table_overrides.get(resname)
        if existing is not None and abs(float(existing) - net_charge) > tolerance:
            raise ReplacementLigandError(
                f"Ligand {spec.resname}: replacement_ligands resolves to {net_charge:+.6f}, "
                f"but charge.charge_table_overrides contains {float(existing):+.6f}."
            )
        cfg.charge.charge_table_overrides[resname] = net_charge

        if spec.itp_path not in cfg.topology.ligand_itp_paths:
            cfg.topology.ligand_itp_paths.append(spec.itp_path)
        if moleculetype != resname:
            existing_type = cfg.topology.moleculetype_overrides.get(resname)
            if existing_type and existing_type != moleculetype:
                raise ReplacementLigandError(
                    f"Ligand {spec.resname}: ITP resolves to moleculetype '{moleculetype}', "
                    f"but topology.moleculetype_overrides specifies '{existing_type}'."
                )
            cfg.topology.moleculetype_overrides[resname] = moleculetype

        resolved.append(
            ResolvedReplacementLigand(
                selector=spec.resname,
                resname=resname,
                occurrence=occurrence,
                moleculetype=moleculetype,
                net_charge=net_charge,
                itp_charge=itp_charge,
                itp_path=spec.itp_path,
                charge_source=source,
            )
        )
    return resolved
