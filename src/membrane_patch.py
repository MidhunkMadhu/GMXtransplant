"""A reference bilayer cut to the host's own box, set up for re-equilibration.

The membrane mode tiles the reference bilayer onto the host box. When the two
boxes differ, lipids that cross the host box edge meet the periodic images of
lipids from the other side; those are removed, so the new membrane can hold
noticeably fewer lipids than fit (about 20 per leaflet for a 12.85 nm patch cut
to 10.6 nm). This module cuts the reference to the host box at the reference's
own density, with equal leaflets, and writes a GROMACS folder that minimizes it
and re-equilibrates it at that fixed area. Its last frame, used as the membrane
mode's reference_gro, fits the host box with nothing lost at the edges:

    gmxtransplant membrane-patch -i membrane.yaml -o patch/
    (run patch/sample_run.sh or patch/submit.sh)
    membrane: reference_gro: patch/patch3_npt.gro
"""
from __future__ import annotations

from collections import Counter
from pathlib import Path
import json
import os
import shutil

import numpy as np

from charmprot import inspect_system
from config import ConfigError
from membrane import (REMOVED_CATEGORIES, _cell, _definitions, default_reference_path, load_membrane_config,
                      load_reference, tile)

STAGES = ("patch0_minimization", "patch1_nvt", "patch2_npt_restrained", "patch3_npt")
FINAL = STAGES[-1] + ".gro"


def _origins(cell, steps=24):
    return [(i / steps) * cell[0] + (j / steps) * cell[1] for i in range(steps) for j in range(steps)]


def choose_cut(reference, host_cell):
    """The lateral registration of the reference whose copy in the host cell has equal
    leaflets, the reference's density and its composition. Returns (origin, placed, report)."""
    host_area = abs(np.linalg.det(host_cell))
    reference_area = abs(np.linalg.det(reference.cell))
    composition = {s: Counter(reference.lipids[i][0].name for i in range(len(reference.lipids))
                              if reference.leaflets[i] == s) for s in ("upper", "lower")}
    target = {s: sum(composition[s].values()) * host_area / reference_area for s in composition}
    best = None
    for origin in _origins(reference.cell):
        placed = tile(reference.heads[:, :2], reference.cell, host_cell, origin)
        counts = {s: Counter() for s in ("upper", "lower")}
        for k, _ in placed:
            counts[reference.leaflets[k]][reference.lipids[k][0].name] += 1
        totals = {s: sum(counts[s].values()) for s in counts}
        drift = sum(abs(counts[s][name] / max(totals[s], 1) - n / sum(composition[s].values()))
                    for s in counts for name, n in composition[s].items())
        score = (abs(totals["upper"] - totals["lower"]),
                 round(sum(abs(totals[s] - target[s]) for s in totals)), round(drift, 3))
        if best is None or score < best[0]:
            best = (score, origin, placed, counts, totals)
    _, origin, placed, counts, totals = best
    report = {"host_area_nm2": round(host_area / 100, 3),
              "reference_area_nm2": round(reference_area / 100, 3),
              "target_per_leaflet": {s: round(v, 1) for s, v in target.items()},
              "lipids": {s: dict(counts[s]) for s in counts},
              "area_per_lipid_angstrom2": {s: round(host_area / totals[s], 2) for s in totals},
              "origin_angstrom": [round(float(v), 3) for v in origin]}
    return origin, placed, report


def _gro_box(cell, height):
    """GRO box line (nm): rectangular, or the nine triclinic components."""
    a, b = cell / 10.0
    c = height / 10.0
    if abs(a[1]) < 1e-6 and abs(b[0]) < 1e-6:
        return f"{a[0]:10.5f}{b[1]:10.5f}{c:10.5f}\n"
    return (f"{a[0]:10.5f}{b[1]:10.5f}{c:10.5f}{a[1]:10.5f}{0.0:10.5f}{b[0]:10.5f}"
            f"{0.0:10.5f}{0.0:10.5f}{0.0:10.5f}\n")


def _write_gro(path, title, molecules, cell, height):
    lines, atom = [], 0
    for residue, (definition, coords) in enumerate(molecules, 1):
        for item, xyz in zip(definition.atoms, coords / 10.0):
            atom += 1
            lines.append(f"{residue % 100000:5d}{item.resname[:5]:<5s}{item.atom_name[:5]:>5s}{atom % 100000:5d}"
                         f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}")
    Path(path).write_text(f"{title}\n{atom}\n" + "\n".join(lines) + "\n" + _gro_box(cell, height))
    return atom


def _ndx(groups):
    text = ""
    for name, (first, last) in groups.items():
        ids = list(range(first, last + 1))
        text += f"[ {name} ]\n" + "".join(" ".join(map(str, ids[i:i + 15])) + "\n" for i in range(0, len(ids), 15))
    return text


COMMON = """cutoff-scheme           = Verlet
nstlist                 = 20
rlist                   = 1.2
vdwtype                 = Cut-off
vdw-modifier            = Force-switch
rvdw_switch             = 1.0
rvdw                    = 1.2
coulombtype             = PME
rcoulomb                = 1.2
"""


def _md(temperature):
    return f"""nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
nstcalcenergy           = 100
nstenergy               = 1000
nstlog                  = 1000
nstxout-compressed      = 50000
tcoupl                  = v-rescale
tc_grps                 = MEMB SOLV
tau_t                   = 1.0 1.0
ref_t                   = {temperature:g} {temperature:g}
constraints             = h-bonds
constraint_algorithm    = LINCS
nstcomm                 = 100
comm_mode               = linear
comm_grps               = MEMB SOLV
"""


NPT = """pcoupl                  = C-rescale
pcoupltype              = semiisotropic
tau_p                   = 5.0
compressibility         = 0.0  4.5e-5      ; lateral area fixed at the host box, z free
ref_p                   = 1.0  1.0
refcoord_scaling        = com
"""


def write_mdps(folder, temperature, nanoseconds, restraints):
    """Minimization, 250 ps NVT and 1 ns NPT with the headgroups held in z (when the
    lipid ITPs offer CHARMM-GUI's POSRES_FC_LIPID/DIHRES_FC), then unrestrained NPT."""
    def define(lipid, dihedral):
        parts = []
        if "posres" in restraints:
            parts += ["-DPOSRES", f"-DPOSRES_FC_LIPID={lipid:.1f}"]
        if "dihres" in restraints:
            parts += ["-DDIHRES", f"-DDIHRES_FC={dihedral:.1f}"]
        return f"define                  = {' '.join(parts)}\n" if parts else ""
    steps = int(round(nanoseconds * 500000))
    texts = {
        STAGES[0]: "; Remove the contacts where the reference was cut to the host box\n"
                   + define(1000.0, 1000.0) + "integrator              = steep\nemtol                   = 500.0\n"
                   "nsteps                  = 50000\n" + COMMON + "constraints             = none\n",
        STAGES[1]: "; NVT, 250 ps at 1 fs, new velocities, headgroups held in z\n" + define(1000.0, 1000.0)
                   + "integrator              = md\ndt                      = 0.001\nnsteps                  = 250000\n"
                   + COMMON + _md(temperature)
                   + f"gen-vel                 = yes\ngen-temp                = {temperature:g}\ngen-seed                = -1\n",
        STAGES[2]: "; NPT, 1 ns at 2 fs, restraints lowered, lateral area fixed\n" + define(400.0, 200.0)
                   + "integrator              = md\ndt                      = 0.002\nnsteps                  = 500000\n"
                   + COMMON + _md(temperature) + NPT + "continuation            = yes\n",
        STAGES[3]: f"; NPT, {nanoseconds:g} ns at 2 fs, no restraints, lateral area fixed at the host box\n"
                   + "integrator              = md\ndt                      = 0.002\n"
                   + f"nsteps                  = {steps}\n" + COMMON + _md(temperature) + NPT
                   + "continuation            = yes\n",
    }
    for name, text in texts.items():
        (Path(folder) / f"{name}.mdp").write_text(text)


def write_run_script(folder, nanoseconds):
    lines = ["#!/bin/bash",
             "# Sample run script written by GMXtransplant (membrane-patch). It minimizes the cut",
             f"# bilayer and re-equilibrates it at the host's lateral area ({nanoseconds:g} ns unrestrained NPT).",
             "# On a cluster, submit it through your scheduler (gmxtransplant slurm-script).",
             f"# Then use {FINAL} as membrane.reference_gro.",
             "set -e", 'cd "$(dirname "$0")"', 'GMX=${GMX:-gmx}', 'MDRUN=${MDRUN:-"$GMX mdrun"}',
             "previous=patch_start.gro", "checkpoint="]
    for name in STAGES:
        lines += [f'echo "==> {name}"',
                  f"$GMX grompp -f {name}.mdp -o {name}.tpr -c $previous -r patch_start.gro -p topol.top "
                  "-n index.ndx -maxwarn 0 $checkpoint",
                  f"$MDRUN -v -deffnm {name}", f"previous={name}.gro"]
        if name != STAGES[0]:
            lines.append(f'checkpoint="-t {name}.cpt"')
    (Path(folder) / "sample_run.sh").write_text("\n".join(lines) + "\n")
    os.chmod(Path(folder) / "sample_run.sh", 0o755)


def build_patch(config_path, output_dir, nanoseconds=50.0):
    """Write the patch folder; returns its report."""
    cfg = load_membrane_config(config_path)
    if cfg.lipid_forcefield:
        raise ConfigError("membrane-patch does not take lipid_forcefield yet; build the patch for lipids whose "
                          "parameters are in the host forcefield.itp")
    if not 1 <= nanoseconds <= 1000:
        raise ConfigError("--ns must be between 1 and 1000")
    output = Path(output_dir).resolve()
    if output.exists() and any(output.iterdir()):
        raise ConfigError(f"{output} is not empty; choose a new folder")
    reference_path = cfg.reference_gro or str(default_reference_path())
    host = inspect_system(cfg.host, gro=cfg.host_gro, overrides=cfg.environment_overrides)
    definitions, extra = _definitions(host, cfg)
    membrane = sorted({m.name for m in host.molecules if m.category in REMOVED_CATEGORIES} | set(extra))
    reference = load_reference(reference_path, definitions, cfg, membrane)
    water = [m.definition for m in host.molecules if m.category == "water"]
    if not water:
        raise ConfigError("The host has no water moleculetype to hydrate the patch with")
    water = water[0]
    if len(reference.waters) == 0:
        raise ConfigError(f"The reference {reference_path} has no three-site water; the patch needs its hydration")
    host_cell, _ = _cell(host.universe.dimensions)
    origin, placed, report = choose_cut(reference, host_cell)
    height = reference.height
    lipids = []
    for k, offset in placed:
        definition, whole = reference.lipids[k]
        coords = whole.copy()
        coords[:, :2] += offset
        coords[:, 2] += height / 2
        lipids.append((definition, coords))
    order = {name: i for i, name in enumerate(dict.fromkeys(d.name for d, _ in reference.lipids))}
    lipids.sort(key=lambda item: order[item[0].name])
    waters = []
    for k, offset in tile(reference.waters[:, 0, :2], reference.cell, host_cell, origin):
        coords = reference.waters[k].copy()
        coords[:, :2] += offset
        coords[:, 2] += height / 2
        waters.append((water, coords))
    output.mkdir(parents=True, exist_ok=True)
    toppar = output / "toppar"
    toppar.mkdir()
    used = {d.name for d, _ in lipids} | {water.name}
    includes = []
    top_path = (Path(host.folder) / "topol.top").resolve()
    for path in host.files:
        if Path(path).resolve() == top_path:
            continue
        names = {d.name for d in definitions.values() if Path(d.source_path).resolve() == Path(path).resolve()}
        if names and not names & used:
            continue
        shutil.copy2(path, toppar / Path(path).name)
        includes.append(Path(path).name)
    for name in sorted(used & set(extra)):
        shutil.copy2(extra[name].source_path, toppar / Path(extra[name].source_path).name)
        includes.append(Path(extra[name].source_path).name)
    counts = Counter(d.name for d, _ in lipids)
    top = ["; Bilayer patch cut to the host box by GMXtransplant (membrane-patch)"]
    top += [f'#include "toppar/{name}"' for name in includes]
    top += ["", "[ system ]", "Bilayer patch in the host box", "", "[ molecules ]"]
    top += [f"{name:<8s}{n}" for name, n in counts.items()] + [f"{water.name:<8s}{len(waters)}"]
    (output / "topol.top").write_text("\n".join(top) + "\n")
    lipid_atoms = sum(len(c) for _, c in lipids)
    total = _write_gro(output / "patch_start.gro",
                       f"Bilayer patch cut from {Path(reference_path).name} to the host box",
                       lipids + waters, host_cell, height)
    (output / "index.ndx").write_text(_ndx({"System": (1, total), "MEMB": (1, lipid_atoms),
                                            "SOLV": (lipid_atoms + 1, total)}))
    texts = "".join(Path(d.source_path).read_text() for d in {d.name: d for d, _ in lipids}.values())
    restraints = {key for key, macro in (("posres", "POSRES_FC_LIPID"), ("dihres", "DIHRES_FC")) if macro in texts}
    temperature = cfg.equilibration.temperature
    write_mdps(output, temperature, nanoseconds, restraints)
    write_run_script(output, nanoseconds)
    report.update({"reference": reference_path, "host": str(host.gro), "waters": len(waters), "atoms": total,
                   "temperature": temperature, "unrestrained_ns": nanoseconds, "stages": list(STAGES),
                   "final": FINAL, "headgroup_restraints": sorted(restraints)})
    (output / "patch_report.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def main(argv=None):
    import argparse
    parser = argparse.ArgumentParser(prog="gmxtransplant membrane-patch",
                                     description="Cut the membrane mode's reference bilayer to the host box and "
                                                 "write a GROMACS folder that re-equilibrates it at that area.")
    parser.add_argument("-i", "--input", required=True, help="membrane YAML (host, reference_gro, lipid_itps ...)")
    parser.add_argument("-o", "--output", required=True, help="new folder for the patch")
    parser.add_argument("--ns", type=float, default=50.0, help="unrestrained NPT length in ns (default 50)")
    args = parser.parse_args(argv)
    try:
        report = build_patch(args.input, args.output, args.ns)
    except ConfigError as exc:
        print(f"[CONFIG ERROR] {exc}")
        return 2
    lipids = report["lipids"]
    print("[patch] " + "; ".join(f"{side} " + ", ".join(f"{n} {k}" for k, n in lipids[side].items())
                                 for side in ("upper", "lower"))
          + f"; area per lipid {report['area_per_lipid_angstrom2']['upper']:.1f} / "
            f"{report['area_per_lipid_angstrom2']['lower']:.1f} A^2; {report['waters']} waters")
    out = Path(args.output).resolve()
    print(f"[patch] wrote {out}: run sample_run.sh (or gmxtransplant slurm-script {out}), then set "
          f"membrane.reference_gro to {out / FINAL}")
    return 0
