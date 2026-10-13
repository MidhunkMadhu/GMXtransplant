"""GROMACS equilibration folders written next to a built system.

Every mode that builds a system can write `equilibration/`: the CHARMM-GUI
Membrane Builder step6 schedule (minimization, then NVT and NPT stages whose
position restraints are released step by step), a copy of the topology and
coordinates, a sample run script and an analysis script. Nothing is run here.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
import os
import shutil

import numpy as np

from config import ConfigError


@dataclass
class EquilibrationStage:
    """One stage; restraints in kJ/mol/nm^2 (POSRES_FC_BB, POSRES_FC_SC, POSRES_FC_LIPID/DIHRES_FC)."""

    name: str = ""
    integrator: str = "npt"      # steep (minimization) | nvt | npt
    nsteps: int = 0
    dt: float = 0.002            # ps; not used by steep
    backbone: float = 0.0
    sidechain: float = 0.0
    lipid: float = 0.0


def default_stages():
    # The CHARMM-GUI Membrane Builder schedule: the protein is released slowly, so it
    # keeps its OPM placement while the lipids pack around it; 1.875 ns of MD.
    return [
        EquilibrationStage("step6.0_minimization", "steep", 5000, 0.0, 4000.0, 2000.0, 1000.0),
        EquilibrationStage("step6.1_equilibration", "nvt", 125000, 0.001, 4000.0, 2000.0, 1000.0),
        EquilibrationStage("step6.2_equilibration", "nvt", 125000, 0.001, 2000.0, 1000.0, 400.0),
        EquilibrationStage("step6.3_equilibration", "npt", 125000, 0.001, 1000.0, 500.0, 400.0),
        EquilibrationStage("step6.4_equilibration", "npt", 250000, 0.002, 500.0, 200.0, 200.0),
        EquilibrationStage("step6.5_equilibration", "npt", 250000, 0.002, 200.0, 50.0, 40.0),
        EquilibrationStage("step6.6_equilibration", "npt", 250000, 0.002, 50.0, 0.0, 0.0),
    ]


@dataclass
class EquilibrationSpec:
    enabled: bool = True
    output_dir: str = "equilibration"    # folder inside the output folder
    temperature: float = 310.0           # K
    stages: list = field(default_factory=default_stages)


def _number(value, key, positive=True):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not np.isfinite(value):
        raise ConfigError(f"{key} must be a number")
    if positive and value <= 0:
        raise ConfigError(f"{key} must be greater than 0")
    return float(value)


def parse_equilibration(raw, section="membrane", enabled=True):
    """Validate an `equilibration:` mapping of the `section` YAML section."""
    where = f"{section}.equilibration"
    if raw is None:
        return EquilibrationSpec(enabled=enabled)
    if isinstance(raw, EquilibrationSpec):
        raw = asdict(raw)
    if not isinstance(raw, dict):
        raise ConfigError(f"{where} must be a mapping")
    unknown = set(raw) - set(EquilibrationSpec.__dataclass_fields__)
    if unknown:
        raise ConfigError(f"Unknown {where} options: {sorted(unknown)}")
    stages = raw.get("stages")
    spec = EquilibrationSpec(**{"enabled": enabled, **{k: v for k, v in raw.items() if k != "stages"}})
    if stages is not None:
        if not isinstance(stages, list) or not stages:
            raise ConfigError(f"{where}.stages must be a nonempty list")
        parsed = []
        for number, stage in enumerate(stages):
            label = f"{where}.stages[{number}]"
            if isinstance(stage, EquilibrationStage):
                stage = asdict(stage)
            if not isinstance(stage, dict):
                raise ConfigError(f"{label} must be a mapping")
            unknown = set(stage) - set(EquilibrationStage.__dataclass_fields__)
            if unknown:
                raise ConfigError(f"Unknown {label} options: {sorted(unknown)}")
            parsed.append(EquilibrationStage(**stage))
        spec.stages = parsed
    if type(spec.enabled) is not bool:
        raise ConfigError(f"{where}.enabled must be true or false")
    if (not isinstance(spec.output_dir, str) or not spec.output_dir
            or Path(spec.output_dir).name != spec.output_dir or spec.output_dir in (".", "..")):
        raise ConfigError(f"{where}.output_dir must be a simple folder name")
    spec.temperature = _number(spec.temperature, f"{where}.temperature")
    names = set()
    for number, stage in enumerate(spec.stages):
        label = f"{where}.stages[{number}]"
        if not isinstance(stage.name, str) or not stage.name or Path(stage.name).name != stage.name:
            raise ConfigError(f"{label}.name must be a simple file name")
        if stage.name in names:
            raise ConfigError(f"{label}.name {stage.name} is used twice")
        names.add(stage.name)
        if stage.integrator not in ("steep", "nvt", "npt"):
            raise ConfigError(f"{label}.integrator must be steep, nvt or npt")
        if type(stage.nsteps) is not int or stage.nsteps <= 0:
            raise ConfigError(f"{label}.nsteps must be a whole number > 0")
        for key in ("backbone", "sidechain", "lipid"):
            setattr(stage, key, _number(getattr(stage, key), f"{label}.{key}", positive=False))
            if getattr(stage, key) < 0:
                raise ConfigError(f"{label}.{key} must not be negative")
        stage.dt = _number(stage.dt, f"{label}.dt", positive=False)
        if stage.integrator != "steep" and not 0 < stage.dt <= 0.004:
            raise ConfigError(f"{label}.dt must be between 0 and 0.004 ps")
    return spec


def write_equilibration(spec, stage_dir, topology, solute_definitions, index_groups,
                        purpose="Equilibration of the built system", extra_mdp=""):
    """A self-contained GROMACS equilibration folder; returns its report.

    `extra_mdp` is appended to every MD stage (e.g. a pull-code restraint).
    """
    from repack import _mdp, restraint_coverage
    stage_dir = Path(stage_dir)
    bundle = stage_dir / spec.output_dir
    bundle.mkdir()
    shutil.copytree(topology.toppar_dir, bundle / "toppar")
    shutil.copy2(topology.top_path, bundle / "topol.top")
    shutil.copy2(stage_dir / "step5_input.gro", bundle / "step5_input.gro")
    shutil.copy2(stage_dir / "index.ndx", bundle / "index.ndx")
    warnings_, coverage = [], {}
    for definition in solute_definitions:
        info = restraint_coverage(definition)
        coverage[definition.name] = info
        if info["heavy_atoms"] and info["missing_heavy_atoms"]:
            warnings_.append(f"{definition.name}: {len(info['missing_heavy_atoms'])} of {info['heavy_atoms']} heavy "
                             "atoms have no POSRES restraint and move freely during equilibration")
    groups = (("SOLU_MEMB SOLV", "SOLU_MEMB SOLV") if {"SOLU_MEMB", "SOLV"} <= set(index_groups)
              else ("System", "System"))
    first_md = next((i for i, s in enumerate(spec.stages) if s.integrator != "steep"), None)
    for number, stage in enumerate(spec.stages):
        define = (f"-DPOSRES -DPOSRES_FC_BB={stage.backbone:.1f} -DPOSRES_FC_SC={stage.sidechain:.1f} "
                  f"-DPOSRES_FC_LIPID={stage.lipid:.1f} -DDIHRES -DDIHRES_FC={stage.lipid:.1f}")
        text = _mdp(stage, spec, number == first_md, groups, define=define)
        if extra_mdp and stage.integrator != "steep":
            text += extra_mdp
        (bundle / f"{stage.name}.mdp").write_text(text)
    md_ps = sum(s.nsteps * s.dt for s in spec.stages if s.integrator != "steep")
    script = ["#!/bin/bash",
              "# Sample run script written by GMXtransplant. Adapt it to your machine, or submit",
              "# the same grompp/mdrun lines through your scheduler (gmxtransplant slurm-script).",
              f"# {purpose}, {md_ps / 1000:g} ns of MD. The protein is held and released step by",
              "# step (CHARMM-GUI schedule) while the rest of the system relaxes around it.",
              "# Restraints (kJ/mol/nm^2): backbone / side chains / lipid P z + dihedrals"]
    script += [f"#   {s.name}: {s.backbone:g} / {s.sidechain:g} / {s.lipid:g}" for s in spec.stages]
    script += ["# Afterwards, check the stages with analyze_equilibration.py and continue production",
               "# from the last GRO with the same topol.top.",
               "set -e", 'cd "$(dirname "$0")"', 'GMX=${GMX:-gmx}', 'MDRUN=${MDRUN:-"$GMX mdrun"}',
               "previous=step5_input.gro"]
    for stage in spec.stages:
        script += [f'echo "==> {stage.name}"',
                   f"$GMX grompp -f {stage.name}.mdp -o {stage.name}.tpr -c $previous -r step5_input.gro "
                   "-p topol.top -n index.ndx -maxwarn 0",
                   f"$MDRUN -v -deffnm {stage.name}", f"previous={stage.name}.gro"]
    (bundle / "sample_run.sh").write_text("\n".join(script) + "\n")
    os.chmod(bundle / "sample_run.sh", 0o755)
    analysis = bundle / "analyze_equilibration.py"
    analysis.write_text(ANALYSIS_SCRIPT)
    os.chmod(analysis, 0o755)
    return {"directory": spec.output_dir, "md_length_ps": md_ps, "temperature": spec.temperature,
            "coupling_groups": list(groups), "stages": [asdict(s) for s in spec.stages],
            "restraint_coverage": coverage, "warnings": warnings_}


# A standalone script (numpy only), so it runs on a cluster next to the GROMACS output.
ANALYSIS_SCRIPT = r'''#!/usr/bin/env python3
"""Check an equilibration written by GMXtransplant: box, membrane and protein per stage.

Run it in the equilibration folder after sample_run.sh (or on any finished stages):
    python3 analyze_equilibration.py            # table on screen
    python3 analyze_equilibration.py --json     # the same as JSON
It reads step5_input.gro, index.ndx and every stepN*.gro present. Per stage:
  box           x, y, z in nm
  APL           projected area per lipid of each leaflet, box xy area / phospholipids
                (the protein's cross-section is included, so this is an upper bound)
  P-P           mean z distance between the leaflets' phosphorus planes, nm
  RMSD          protein C-alpha RMSD from step5_input.gro after a best fit, nm, and the
                drift of the protein centre in z (membrane placement), nm
Missing stages and stages with no membrane or no protein are reported as such.
"""
import glob
import json
import os
import re
import sys

import numpy as np


def read_gro(path):
    with open(path) as handle:
        lines = handle.read().splitlines()
    n = int(lines[1])
    rows = lines[2:2 + n]
    names = [r[10:15].strip() for r in rows]
    resids = [int(r[0:5]) for r in rows]
    xyz = np.array([[float(r[20:28]), float(r[28:36]), float(r[36:44])] for r in rows])
    box = [float(v) for v in lines[2 + n].split()]
    return names, resids, xyz, box


def read_ndx(path):
    groups, current = {}, None
    if not os.path.isfile(path):
        return groups
    for line in open(path):
        match = re.match(r"\s*\[\s*(.+?)\s*\]", line)
        if match:
            current = groups.setdefault(match.group(1), [])
        elif current is not None:
            current.extend(int(v) - 1 for v in line.split())
    return groups


def kabsch_rmsd(a, b):
    a = a - a.mean(0)
    b = b - b.mean(0)
    u, s, vt = np.linalg.svd(a.T @ b)
    d = np.sign(np.linalg.det(u @ vt))
    s[-1] *= d
    return float(np.sqrt(max((a ** 2).sum() + (b ** 2).sum() - 2 * s.sum(), 0.0) / len(a)))


def stage_key(path):
    match = re.search(r"step(\d+)(?:\.(\d+))?", os.path.basename(path))
    return (int(match.group(1)), int(match.group(2) or 0)) if match else (999, 0)


def main():
    groups = read_ndx("index.ndx")
    names0, _, ref, _ = read_gro("step5_input.gro")
    protein = groups.get("Protein", [])
    ca = [i for i in protein if names0[i] == "CA"]
    membrane = groups.get("Membrane") or groups.get("MEMB") or []
    phosphorus = [i for i in membrane if names0[i] == "P"]
    files = ["step5_input.gro"] + sorted((p for p in glob.glob("step6*.gro") + glob.glob("step7*.gro")),
                                         key=stage_key)
    rows = []
    for path in files:
        _, _, xyz, box = read_gro(path)
        row = {"stage": os.path.splitext(path)[0], "box_nm": [round(v, 3) for v in box[:3]]}
        if len(phosphorus) >= 6:
            z = xyz[phosphorus, 2]
            # Leaflets from the circular mean around the membrane centre.
            angle = 2 * np.pi * z / box[2]
            centre = np.arctan2(np.sin(angle).mean(), np.cos(angle).mean()) * box[2] / (2 * np.pi)
            dz = z - centre
            dz -= box[2] * np.round(dz / box[2])
            upper, lower = dz > 0, dz <= 0
            area = box[0] * box[1]
            row["apl_nm2"] = {"upper": round(area / max(upper.sum(), 1), 4),
                              "lower": round(area / max(lower.sum(), 1), 4)}
            row["pp_thickness_nm"] = round(float(dz[upper].mean() - dz[lower].mean()), 3)
        if len(ca) >= 3:
            # Minimum-image displacements from the input, so a protein crossing the box edge
            # between frames does not count as moving a box length.
            lengths = np.array(box[:3])
            moved = xyz[ca] - ref[ca]
            moved -= lengths * np.round(moved / lengths)
            current = ref[ca] + moved
            row["ca_rmsd_nm"] = round(kabsch_rmsd(current, ref[ca]), 4)
            row["protein_z_drift_nm"] = round(float(current[:, 2].mean() - ref[ca, 2].mean()), 4)
        rows.append(row)
    if "--json" in sys.argv:
        print(json.dumps(rows, indent=2))
        return 0
    print(f"{'stage':<24} {'box x/y/z (nm)':<22} {'APL up/low (nm2)':<18} {'P-P (nm)':<9} "
          f"{'CA RMSD (nm)':<13} {'z drift (nm)'}")
    for row in rows:
        box = "/".join(f"{v:.2f}" for v in row["box_nm"])
        apl = (f"{row['apl_nm2']['upper']:.3f}/{row['apl_nm2']['lower']:.3f}" if "apl_nm2" in row else "-")
        pp = f"{row['pp_thickness_nm']:.2f}" if "pp_thickness_nm" in row else "-"
        rmsd = f"{row['ca_rmsd_nm']:.3f}" if "ca_rmsd_nm" in row else "-"
        drift = f"{row['protein_z_drift_nm']:+.3f}" if "protein_z_drift_nm" in row else "-"
        print(f"{row['stage']:<24} {box:<22} {apl:<18} {pp:<9} {rmsd:<13} {drift}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
'''
