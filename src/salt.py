"""Set the salt of a prepared system: ion types and concentration.

The bulk ions are removed (ions bound to the protein can be kept), then formula
units of the chosen salt (NaCl, KCl, CaCl2, MgCl2 ...) replace random bulk waters
away from the protein, ligands and membrane until the concentration is reached:
units / waters x 55.51 mol/L, counting the waters left after replacement. The
system is then neutralized. Ion types the host does not have come from an ITP
given in ion_itps or from the bundled CHARMM36 ions, with their parameters (atom
type, NBFIX pairs) added to the force field only where the host lacks them.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

from charmprot import _heavy, _indices
from config import ConfigError
from equilibration import EquilibrationSpec
from itp import parse_itp
from system_edit import (MEMBRANE_CATEGORIES, SOLVENT_CATEGORIES, WATER_MOLARITY, balance_charge, concentration,
                         ion_block, load_section, molecules_near, net_charge, order, positive, replace_waters,
                         run_edit, salt_summary, salt_units, write_system)
import charmm36_library as c36


def _no_equilibration():
    return EquilibrationSpec(enabled=False)


@dataclass
class SaltSpec:
    host: str
    host_gro: str | None = None
    concentration: float = 0.15          # mol/L of formula units (CaCl2 counts once)
    cation: str = "SOD"                  # moleculetype names (CHARMM36: SOD POT LIT CES RUB CAL MG BAR ZN2 CD2)
    anion: str = "CLA"
    replace_existing: bool = True        # remove the host's bulk ions first
    keep_bound_distance: float = 0.0     # A; ions this close to a protein heavy atom are kept
    ion_itps: list[str] = field(default_factory=list)   # ITPs of ion types the host lacks
    solute_distance: float = 5.0         # A, new ion to protein/ligand heavy atoms
    membrane_distance: float = 5.0       # A, new ion to lipid heavy atoms
    ion_spacing: float = 5.0             # A between ions
    neutralize: bool = True
    target_net_charge: float = 0.0
    seed: int = 42
    environment_overrides: dict = field(default_factory=dict)
    output_dir: str = "salt_output"
    equilibration: EquilibrationSpec = field(default_factory=_no_equilibration)


def load_config(path, check_paths=True, output_root=None):
    cfg = load_section(path, "salt", SaltSpec, check_paths, output_root, optional_paths=("host_gro", "ion_itps"))
    positive(cfg, "salt", "concentration", "keep_bound_distance", "solute_distance", "membrane_distance",
             "ion_spacing", allow_zero=True)
    if cfg.concentration >= 5:
        raise ConfigError("salt.concentration must be below 5 mol/L")
    for key in ("cation", "anion"):
        if not isinstance(getattr(cfg, key), str) or not getattr(cfg, key).strip():
            raise ConfigError(f"salt.{key} must be an ion moleculetype name")
    return cfg


def _ion_definition(name, host, cfg, work, created):
    """The ion's definition: host, ion_itps, or the bundled CHARMM36 library (written into work/)."""
    if name in host.definitions:
        return host.definitions[name]
    for path in cfg.ion_itps:
        definitions = parse_itp(path)
        if name in definitions:
            return definitions[name]
    target = Path(work) / f"{name}.itp"
    target.write_text(c36.ion_itp(name))
    created.append(target)
    return parse_itp(str(target))[name]


def run(config_path, dry_run=False, output_root=None):
    cfg = load_config(config_path, check_paths=not dry_run, output_root=output_root)

    def build(stage, host, report):
        work = stage.parent / "work"
        work.mkdir()
        created = []
        cation = _ion_definition(cfg.cation, host, cfg, work, created)
        anion = _ion_definition(cfg.anion, host, cfg, work, created)
        for definition, sign in ((cation, 1), (anion, -1)):
            if definition.atom_count != 1 or definition.charge * sign <= 0:
                raise ConfigError(f"{definition.name} must be a single-atom {'cation' if sign > 0 else 'anion'}")
        units = salt_units(cation.charge, anion.charge)
        box = host.universe.dimensions
        protected = [m for m in host.molecules if m.category in ("protein", "ligand")]
        membrane = [m for m in host.molecules if m.category in MEMBRANE_CATEGORIES]
        solvent = [m for m in host.molecules if m.category in SOLVENT_CATEGORIES]
        protected_heavy = _heavy(host.universe.atoms[_indices(protected)]).positions
        membrane_heavy = _heavy(host.universe.atoms[_indices(membrane)]).positions if membrane else []
        host_ions = [m for m in solvent if m.category == "ion"]
        removed = []
        if cfg.replace_existing:
            bound = set()
            if cfg.keep_bound_distance > 0:
                bound = molecules_near(host, host_ions, protected_heavy, cfg.keep_bound_distance, box)
            removed = [m for i, m in enumerate(host_ions) if i not in bound]
            gone = {id(m) for m in removed}
            solvent = [m for m in solvent if id(m) not in gone]
            print(f"[salt] removed {len(removed)} bulk ion(s): "
                  + (", ".join(f"{n} {k}" for k, n in Counter(m.name for m in removed).items()) or "none")
                  + (f"; kept {len(bound)} bound to the protein" if bound else ""))
        waters = sum(m.category == "water" for m in solvent)
        have = Counter(m.name for m in solvent if m.category == "ion")
        existing = min(have[cation.name] // units[0], have[anion.name] // units[1])
        per_unit = sum(units)
        # units + n = c (W - n per_unit) / 55.51, each unit replacing per_unit waters
        wanted = cfg.concentration * waters / WATER_MOLARITY
        add = int(round((wanted - existing) / (1 + cfg.concentration * per_unit / WATER_MOLARITY)))
        avoid = [(protected_heavy, cfg.solute_distance)]
        if len(membrane_heavy):
            avoid.append((membrane_heavy, cfg.membrane_distance))
        added = []
        if add > 0:
            solvent, added, shortfall = replace_waters(host, solvent, [(cation, add * units[0]), (anion, add * units[1])],
                                                       avoid, box, spacing=cfg.ion_spacing, seed=cfg.seed)
            if shortfall:
                report["warnings"].append("Not every ion found room in the bulk water: short "
                                          + ", ".join(f"{n} {k}" for k, n in shortfall.items()))
        elif add < 0:
            surplus = [m for m in solvent if m.name in (cation.name, anion.name)]
            drop, counts = set(), Counter()
            for mol in surplus:
                need = -add * (units[0] if mol.name == cation.name else units[1])
                if counts[mol.name] < need:
                    drop.add(id(mol))
                    counts[mol.name] += 1
            removed += [m for m in solvent if id(m) in drop]
            solvent = [m for m in solvent if id(m) not in drop]
        print(f"[salt] {cfg.concentration:g} M {cation.name}/{anion.name} ({units[0]}:{units[1]}): "
              f"{max(add, 0)} formula unit(s) added")
        report["salt"] = {"host": salt_summary([(host, m) for m in host.molecules])}
        kept = protected + membrane + solvent
        ordered_charge = net_charge([(host, m) for m in kept]) + sum(d.charge for d, _ in added)
        neutralizers = (cation if abs(cation.charge - 1) < 1e-3 else _ion_definition("SOD", host, cfg, work, created),
                        anion if abs(anion.charge + 1) < 1e-3 else _ion_definition("CLA", host, cfg, work, created))
        if cfg.neutralize:
            solvent, extra, report["charge"] = balance_charge(
                host, solvent, ordered_charge, cfg.target_net_charge, neutralizers[0], neutralizers[1], avoid, box,
                spacing=cfg.ion_spacing, seed=cfg.seed + 7, occupied=[p for _, p in added])
            added += extra
            print(f"[charge] {report['charge']['before']:+.3f} e -> {report['charge']['after']:+.3f} e")
        new = ion_block(added, work)
        ordered = order(host, protected + membrane + solvent, new=new)
        final_waters = sum(m.category == "water" for _, m in ordered)
        final_ions = Counter(m.name for _, m in ordered if m.category == "ion")
        report["salt"]["final"] = {"concentration_molar": concentration(final_waters, final_ions[cation.name],
                                                                        final_ions[anion.name], units),
                                   "waters": final_waters, "ions": dict(final_ions)}
        report["salt"].update({"cation": cation.name, "anion": anion.name, "formula_units": list(units),
                               "target_molar": cfg.concentration,
                               "removed": dict(Counter(m.name for m in removed)),
                               "added": dict(Counter(d.name for d, _ in added)),
                               "ion_sources": {d.name: d.source_path for d in (cation, anion)}})
        print(f"[salt] final {report['salt']['final']['concentration_molar']} M; ions "
              + ", ".join(f"{k} {n}" for k, n in final_ions.items()))
        parameter_files = []
        if created or cfg.ion_itps:
            types = {a.atom_type for d in host.definitions.values() for a in d.atoms}
            texts = [Path(d.source_path).read_text() for d in {cation.name: cation, anion.name: anion,
                                                                **{n.name: n for n in neutralizers}}.values()
                     if d.name not in host.definitions]
            text, supplement = c36.supplement(host.files, texts, types)
            if supplement:
                path = work / "ion_parameters.itp"
                path.write_text(text)
                parameter_files.append(path)
                report["salt"]["parameters_added"] = supplement
        report.update(write_system(stage, host, ordered, "salt", parameter_files=parameter_files,
                                   equilibration=cfg.equilibration, purpose="Equilibration with the new salt"))

    def lines(report):
        salt = report.get("salt")
        if not salt or "final" not in salt:
            return []
        out = ["", f"SALT: {salt['target_molar']:g} M {salt['cation']} / {salt['anion']} "
                   f"(formula {salt['formula_units'][0]}:{salt['formula_units'][1]})",
               "  removed: " + (", ".join(f"{k} {n}" for k, n in salt["removed"].items()) or "none"),
               "  added:   " + (", ".join(f"{k} {n}" for k, n in salt["added"].items()) or "none"),
               f"  final:   {salt['final']['concentration_molar']} M; "
               + ", ".join(f"{k} {n}" for k, n in salt["final"]["ions"].items())
               + f"; {salt['final']['waters']} waters"]
        if salt.get("parameters_added"):
            out.append("  CHARMM36 parameters added: "
                       + ", ".join(f"{k} {n}" for k, n in salt["parameters_added"].items()))
        return out

    return run_edit("salt", cfg, "SET THE SALT", build, "set the salt", lines, dry_run=dry_run)
