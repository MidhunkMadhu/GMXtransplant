"""Hydrogen mass repartitioning (HMR) of a prepared system's topology.

Every hydrogen bonded to a heavy atom gets `hydrogen_mass` (3.024 u by default,
three times the hydrogen mass, as CHARMM-GUI does); the mass is taken from that
heavy atom, so each molecule keeps its total mass. Water is left unchanged by
default (it is rigid under SETTLE). Coordinates are not touched. With h-bond
constraints this allows a 4 fs time step; the host's mdp files are copied with
the production time step changed.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace
from pathlib import Path
import re

from config import ConfigError
from equilibration import EquilibrationSpec
from system_edit import load_section, order, positive, run_edit, write_system

MIN_HEAVY_MASS = 4.0   # u; a heavy atom lighter than this after repartitioning is refused


def _no_equilibration():
    return EquilibrationSpec(enabled=False)


@dataclass
class HmrSpec:
    host: str
    host_gro: str | None = None
    hydrogen_mass: float = 3.024         # u, the mass every bonded hydrogen gets
    skip_water: bool = True              # leave rigid water models unchanged
    time_step: float = 0.004             # ps, written to the copied production mdp
    copy_run_inputs: bool = True         # copy the host's mdp files and README, production dt changed
    environment_overrides: dict = field(default_factory=dict)
    output_dir: str = "hmr_output"
    equilibration: EquilibrationSpec = field(default_factory=_no_equilibration)


def load_config(path, check_paths=True, output_root=None):
    cfg = load_section(path, "hmr", HmrSpec, check_paths, output_root)
    positive(cfg, "hmr", "hydrogen_mass", "time_step")
    if not 1.0 < cfg.hydrogen_mass <= 4.1:
        raise ConfigError("hmr.hydrogen_mass must be between 1 and 4.1 u")
    if cfg.time_step > 0.005:
        raise ConfigError("hmr.time_step must be at most 0.005 ps")
    return cfg


def _is_hydrogen(atom_type, mass):
    return atom_type.upper().startswith("H") and 0.9 < mass < 1.2


def repartition_text(text, hydrogen_mass, skip=frozenset()):
    """Rewrite the masses of every moleculetype in one ITP text; returns (text, report per moleculetype)."""
    lines = text.splitlines(keepends=True)
    molecules, section, current, pending = {}, None, None, False
    for number, raw in enumerate(lines):
        line = raw.split(";", 1)[0].strip()
        if not line or line.startswith("#"):
            continue
        match = re.match(r"\[\s*([^]]+?)\s*\]", line)
        if match:
            section = match[1].lower()
            pending = section == "moleculetype"
            continue
        fields = line.split()
        if pending:
            current, pending = fields[0], False
            molecules[current] = {"atoms": {}, "bonds": [], "rigid": False}
            continue
        if current is None:
            continue
        if section == "atoms":
            if len(fields) < 8:
                raise ConfigError(f"[ atoms ] row without a mass in moleculetype {current}: {line}")
            molecules[current]["atoms"][int(fields[0])] = (number, fields[1], float(fields[7]))
        elif section == "bonds" and len(fields) >= 2:
            molecules[current]["bonds"].append((int(fields[0]), int(fields[1])))
        elif section == "settles":
            molecules[current]["rigid"] = True
    new_mass, report = {}, {}
    for name, mol in molecules.items():
        if name in skip or mol["rigid"]:
            report[name] = {"changed_hydrogens": 0, "note": "left unchanged (water)" if mol["rigid"] else "skipped"}
            continue
        masses = {nr: mass for nr, (_, _, mass) in mol["atoms"].items()}
        before = sum(masses.values())
        hydrogens = 0
        for a, b in mol["bonds"]:
            for h, heavy in ((a, b), (b, a)):
                _, h_type, h_mass = mol["atoms"][h]
                _, heavy_type, heavy_mass = mol["atoms"][heavy]
                if _is_hydrogen(h_type, h_mass) and not _is_hydrogen(heavy_type, heavy_mass):
                    delta = hydrogen_mass - masses[h]
                    masses[h] += delta
                    masses[heavy] -= delta
                    hydrogens += 1
        lightest = min((m for nr, m in masses.items() if not _is_hydrogen(mol["atoms"][nr][1], mol["atoms"][nr][2])
                        and mol["atoms"][nr][2] > 0), default=None)
        if lightest is not None and lightest < MIN_HEAVY_MASS:
            raise ConfigError(f"{name}: a heavy atom would weigh {lightest:.3f} u after repartitioning; "
                              "use a smaller hydrogen_mass")
        if abs(sum(masses.values()) - before) > 1e-6:
            raise ConfigError(f"{name}: repartitioning changed the total mass")
        for nr, (line_number, _, mass) in mol["atoms"].items():
            if abs(masses[nr] - mass) > 1e-9:
                new_mass[line_number] = masses[nr]
        report[name] = {"changed_hydrogens": hydrogens, "total_mass": round(before, 4),
                        "lightest_heavy_atom": round(lightest, 4) if lightest is not None else None}
    out = []
    for number, raw in enumerate(lines):
        if number in new_mass:
            body, sep, comment = raw.rstrip("\n").partition(";")
            fields = body.split()
            fields[7] = f"{new_mass[number]:.4f}"
            raw = (f"{fields[0]:>6} {fields[1]:>10} {fields[2]:>6} {fields[3]:>8} {fields[4]:>6} {fields[5]:>6} "
                   f"{float(fields[6]):>12.6f} {fields[7]:>10}" + "".join(f" {f}" for f in fields[8:])
                   + (f"   ;{comment}" if sep else "") + "\n")
        out.append(raw)
    return "".join(out), report


def _copy_run_inputs(host_folder, stage, time_step):
    """The host's mdp files and README; the production stage gets the HMR time step, same length."""
    copied = []
    for path in sorted(Path(host_folder).iterdir()):
        if not path.is_file() or not (path.suffix == ".mdp" or path.name == "README"):
            continue
        text = path.read_text()
        if path.name.startswith("step7") and path.suffix == ".mdp":
            dt = re.search(r"^\s*dt\s*=\s*([0-9.eE+-]+)", text, re.M)
            steps = re.search(r"^\s*nsteps\s*=\s*([0-9]+)", text, re.M)
            if dt:
                old = float(dt[1])
                text = re.sub(r"^(\s*dt\s*=\s*)[0-9.eE+-]+", lambda m: f"{m[1]}{time_step:g}", text, count=1,
                              flags=re.M)
                if steps:
                    length = int(steps[1]) * old
                    text = re.sub(r"^(\s*nsteps\s*=\s*)[0-9]+",
                                  lambda m: f"{m[1]}{int(round(length / time_step))}", text, count=1, flags=re.M)
                text = f"; GMXtransplant hmr: dt {old:g} -> {time_step:g} ps (hydrogen mass repartitioning)\n" + text
        (stage / path.name).write_text(text)
        copied.append(path.name)
    return copied


def run(config_path, dry_run=False, output_root=None):
    cfg = load_config(config_path, check_paths=not dry_run, output_root=output_root)
    copied_inputs = []

    def build(stage, host, report):
        work = stage.parent / "work"
        work.mkdir()
        skip = {m.name for m in host.molecules if m.category == "water"} if cfg.skip_water else set()
        rewritten, summary = {}, {}
        for path in host.files:
            names = [n for n, d in host.definitions.items() if Path(d.source_path).resolve() == path.resolve()]
            if not names:
                continue
            text, info = repartition_text(path.read_text(), cfg.hydrogen_mass, skip)
            summary.update({k: v for k, v in info.items() if k in names})
            if any(info.get(n, {}).get("changed_hydrogens") for n in names):
                target = work / path.name
                target.write_text(text)
                from itp import parse_itp
                rewritten.update(parse_itp(str(target)))
        replacements = {id(m): (host, replace(m, definition=rewritten[m.name]))
                        for m in host.molecules if m.name in rewritten}
        ordered = order(host, host.molecules, replacements=replacements)
        changed = sum(v.get("changed_hydrogens", 0) * sum(m.name == n for m in host.molecules)
                      for n, v in summary.items())
        print(f"[hmr] hydrogen mass {cfg.hydrogen_mass:g} u in {len(rewritten)} moleculetype(s); "
              f"{changed} hydrogens repartitioned")
        report["hmr"] = {"hydrogen_mass": cfg.hydrogen_mass, "moleculetypes": summary,
                         "hydrogens_repartitioned": changed}
        report.update(write_system(stage, host, ordered, "hmr", equilibration=cfg.equilibration,
                                   purpose="Equilibration of the repartitioned system"))
        if cfg.copy_run_inputs:
            copied_inputs.extend(_copy_run_inputs(host.folder, stage, cfg.time_step))
            report["hmr"]["copied_run_inputs"] = list(copied_inputs)
        report["warnings"].append(
            f"Use dt = {cfg.time_step:g} ps only with constraints = h-bonds; the equilibration stages keep "
            "their own time steps")

    def lines(report):
        info = report.get("hmr")
        if not info:
            return []
        out = ["", f"HYDROGEN MASS REPARTITIONING: {info['hydrogen_mass']:g} u, "
                   f"{info['hydrogens_repartitioned']} hydrogens in the system"]
        for name, item in info["moleculetypes"].items():
            if item.get("changed_hydrogens"):
                out.append(f"  {name:<10} {item['changed_hydrogens']} H per molecule; lightest heavy atom "
                           f"{item['lightest_heavy_atom']:.3f} u")
            else:
                out.append(f"  {name:<10} {item.get('note', 'no bonded hydrogens')}")
        if info.get("copied_run_inputs"):
            out.append(f"  mdp files copied, production dt {cfg.time_step:g} ps: "
                       + ", ".join(info["copied_run_inputs"]))
        return out

    return run_edit("hmr", cfg, "HYDROGEN MASS REPARTITIONING", build, "hydrogen mass repartitioning", lines,
                    dry_run=dry_run)

