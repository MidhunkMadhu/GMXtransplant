"""GROMACS repacking bundle for cholesterol restoration.

Inserting experimental cholesterols removes the lipids that overlapped them
and then trims the leaflet counts, so the restored cholesterols start with
under-packed surroundings. Released straight into MD they can leave their
experimental poses before the lipids close the gaps. The bundle written here
is a short GROMACS equilibration (about 1.4 ns by default) that holds every
heavy atom of the protein, the ligands and the restored cholesterols while the
other lipids, water and ions move freely and repack around them.

Written to <output>/repack/, a self-contained folder:
  toppar/          the run's toppar plus CHLR.itp: the cholesterol ITP with only
                   the name and the position restraints changed (every heavy
                   atom in x, y, z at POSRES_FC_CHLR)
  topol.top        the run's topology with the restored cholesterols as CHLR
  step5_input.gro  the restored system (also the -r reference)
  index.ndx        the run's index
  step6.*.mdp      one file per stage
  sample_run.sh    a sample run script; how to run or submit it is up to the user

Protein and ligands are held through their own CHARMM-GUI POSRES blocks
(POSRES_FC_BB / POSRES_FC_SC); the bundle checks that those blocks cover every
heavy atom and reports any that are not.
"""
from __future__ import annotations

import os
import re
import shutil
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

from classify import classify_resnames
from itp import parse_itp
from topology import scan_toppar_definitions


class RepackError(RuntimeError):
    pass


_SECTION = re.compile(r"^\s*\[\s*([A-Za-z_]+)\s*\]")


def _is_heavy(atom) -> bool:
    # CHARMM/CGenFF hydrogen types all start with H; no heavy-atom type does.
    # CGenFF halogen lone pairs (LPH ...) are massless virtual sites, not atoms to restrain.
    kind = atom.atom_type.upper()
    return not (kind.startswith("H") or kind.startswith("LP"))


def restored_residue_indices(final_universe, inserted_atoms, resname: str) -> List[int]:
    """Residue indices in final_universe of the inserted cholesterols, matched by coordinates,
    in the order of inserted_atoms.residues."""
    box = final_universe.dimensions[:3].astype(float) if final_universe.dimensions is not None else None
    wanted = [res.atoms.positions.astype(float) for res in inserted_atoms.residues]
    found = []
    candidates = [r for r in final_universe.residues if str(r.resname).upper() == resname.upper()]
    for positions in wanted:
        match = None
        for residue in candidates:
            if len(residue.atoms) != len(positions):
                continue
            delta = residue.atoms.positions.astype(float) - positions
            if box is not None and np.all(box > 0):
                delta -= box * np.round(delta / box)
            if np.abs(delta).max() < 0.01:
                match = residue
                break
        if match is None:
            raise RepackError("A restored cholesterol could not be found in the final system")
        found.append(int(match.resindex))
    return found  # in the order of inserted_atoms.residues


def _posres_rows(source_path: str, mtype: str):
    """(atom index, function, parameters) of the [ position_restraints ] rows under #ifdef POSRES."""
    rows, current, section, stack = [], None, None, []
    pending_name = False
    for raw in Path(source_path).read_text().splitlines():
        line = raw.split(";", 1)[0].strip()
        if not line:
            continue
        if line.startswith("#ifdef") or line.startswith("#ifndef"):
            stack.append(line.split()[1] if len(line.split()) > 1 else "")
            continue
        if line.startswith("#endif"):
            if stack:
                stack.pop()
            continue
        if line.startswith("#"):
            continue
        header = _SECTION.match(line)
        if header:
            section = header.group(1).lower()
            pending_name = section == "moleculetype"
            continue
        if pending_name:
            current, pending_name = line.split()[0], False
            continue
        if current == mtype and section == "position_restraints" and "POSRES" in stack:
            fields = line.split()
            rows.append((int(fields[0]), fields[1], fields[2:]))
    return rows


def restraint_coverage(definition) -> dict:
    """How many of a moleculetype's heavy atoms its POSRES block restrains, and with which values."""
    heavy = {a.nr for a in definition.atoms if _is_heavy(a)}
    rows = _posres_rows(definition.source_path, definition.name)
    restrained = {nr for nr, _, _ in rows}
    values = sorted({p for _, _, params in rows for p in params if not re.fullmatch(r"[-+0-9.eE]+", p)})
    return {"heavy_atoms": len(heavy), "restrained_heavy_atoms": len(heavy & restrained),
            "missing_heavy_atoms": sorted(heavy - restrained), "restraint_values": values}


def make_restrained_cholesterol_itp(source_path: str, mtype: str, new_name: str) -> str:
    """The cholesterol ITP renamed, with its POSRES block replaced by x/y/z restraints on every heavy atom."""
    definition = parse_itp(source_path)[mtype]
    text = Path(source_path).read_text()
    if len(parse_itp(source_path)) != 1:
        raise RepackError(f"{source_path} defines more than one moleculetype; expected only {mtype}")
    text, count = re.subn(r"(\[\s*moleculetype\s*\][^\n]*\n(?:\s*;[^\n]*\n)*\s*)" + re.escape(mtype) + r"\b",
                          lambda m: m[1] + new_name, text, count=1, flags=re.I)
    if count != 1:
        raise RepackError(f"Could not rename moleculetype {mtype} in {source_path}")
    # Drop the existing POSRES block (CHARMM-GUI: a z restraint on O3 at POSRES_FC_LIPID).
    text = re.sub(r"#ifdef\s+POSRES\s*\n\s*\[\s*position_restraints\s*\].*?#endif[^\n]*\n?", "", text,
                  flags=re.S | re.I)
    macro = f"POSRES_FC_{new_name}"
    rows = "\n".join(f"{a.nr:5d}     1    {macro}    {macro}    {macro}"
                     for a in definition.atoms if _is_heavy(a))
    header = (f"; {new_name} = {mtype} with identical parameters, used for the cholesterols restored from the\n"
              f"; experimental structure. Only the moleculetype name and the position restraints differ:\n"
              f"; every heavy atom is restrained in x, y and z at {macro}.\n")
    return (header + text.rstrip("\n") + f"\n\n#ifdef POSRES\n[ position_restraints ]\n{rows}\n#endif\n")


def _mdp(stage, cfg, first_dynamics: bool, groups) -> str:
    restraint = float(stage.restraint)
    lipid = float(cfg.lipid_restraint)
    define = (f"-DPOSRES -DPOSRES_FC_BB={restraint:.1f} -DPOSRES_FC_SC={restraint:.1f} "
              f"-DPOSRES_FC_LIPID={lipid:.1f} -DDIHRES -DDIHRES_FC={lipid:.1f} "
              f"-DPOSRES_FC_{cfg.moleculetype}={restraint:.1f}")
    common = ["cutoff-scheme           = Verlet", "rlist                   = 1.2",
              "vdwtype                 = Cut-off", "vdw-modifier            = Force-switch",
              "rvdw_switch             = 1.0", "rvdw                    = 1.2",
              "coulombtype             = PME", "rcoulomb                = 1.2"]
    if stage.integrator == "steep":
        lines = [f"define                  = {define}", "integrator              = steep",
                 "emtol                   = 1000.0", f"nsteps                  = {stage.nsteps}",
                 "nstlist                 = 10", *common, ";",
                 "constraints             = h-bonds", "constraint_algorithm    = LINCS"]
        return "\n".join(lines) + "\n"
    tc, comm = groups
    n = len(tc.split())
    temperature = f"{cfg.temperature:g}"
    lines = [f"define                  = {define}", "integrator              = md",
             f"dt                      = {stage.dt:g}",
             f"nsteps                  = {stage.nsteps}      ; {stage.nsteps * stage.dt:g} ps",
             "nstxout-compressed      = 5000", "nstxout                 = 0", "nstvout                 = 0",
             "nstfout                 = 0", "nstcalcenergy           = 100", "nstenergy               = 1000",
             "nstlog                  = 1000", ";", "nstlist                 = 20", *common, ";",
             "tcoupl                  = v-rescale", f"tc_grps                 = {tc}",
             f"tau_t                   = {' '.join(['1.0'] * n)}",
             f"ref_t                   = {' '.join([temperature] * n)}", ";"]
    if stage.integrator == "npt":
        lines += ["pcoupl                  = C-rescale", "pcoupltype              = semiisotropic",
                  "tau_p                   = 5.0", "compressibility         = 4.5e-5  4.5e-5",
                  "ref_p                   = 1.0     1.0", "refcoord_scaling        = com", ";"]
    lines += ["constraints             = h-bonds", "constraint_algorithm    = LINCS"]
    if first_dynamics:
        lines += ["continuation            = no", ";", "gen-vel                 = yes",
                  f"gen-temp                = {temperature}", "gen-seed                = -1"]
    else:
        lines += ["continuation            = yes"]
    lines += [";", "nstcomm                 = 100", "comm_mode               = linear",
              f"comm_grps               = {comm}"]
    return "\n".join(lines) + "\n"


def _top_with(main_top: str, itp_include: str, molecules) -> str:
    """The run's topol.top with the restrained cholesterol ITP included and [ molecules ] replaced."""
    out, section, placed = [], None, False
    for raw in Path(main_top).read_text().splitlines():
        header = _SECTION.match(raw)
        if header:
            section = header.group(1).lower()
            if section == "molecules":
                out.append(raw)
                out.append("; Compound        #mols")
                out += [f"{name:<16} {count}" for name, count in molecules]
                continue
        if section == "molecules":
            continue
        out.append(raw)
        include = re.match(r'^\s*#include\s+"([^"]+)"', raw)
        if include and not placed and "forcefield" in Path(include.group(1)).name.lower():
            out.append(f'#include "{itp_include}"')
            placed = True
    if not placed:
        raise RepackError("Could not find the force-field include in topol.top to place the restrained "
                          "cholesterol ITP after it")
    return "\n".join(out) + "\n"


def write_repack_bundle(cfg, final_universe, inserted_atoms, topology_result, gro_path: str,
                        ndx_path: Optional[str], ndx_groups, resname: str) -> dict:
    """Write the repack bundle next to the GRO; return a report dictionary."""
    output_dir = Path(gro_path).resolve().parent
    bundle = output_dir / cfg.output_dir
    restored = set(restored_residue_indices(final_universe, inserted_atoms, resname))
    definitions, files = scan_toppar_definitions(topology_result.toppar_dir)

    # Expand [ molecules ] into one entry per molecule, with its first atom index.
    renamed, offset, starts = [], 0, {}
    first_atom_of = {int(final_universe.residues[i].atoms[0].index): i for i in restored}
    for mtype, count in topology_result.molecules_written:
        if mtype not in definitions:
            raise RepackError(f"Moleculetype {mtype} from topol.top was not found in {topology_result.toppar_dir}")
        size = definitions[mtype].atom_count
        for _ in range(count):
            name = mtype
            if offset in first_atom_of:
                if mtype != resname and definitions[mtype].atoms[0].resname.upper() != resname.upper():
                    raise RepackError(f"Restored cholesterol at atom {offset + 1} is in moleculetype {mtype}")
                starts[offset] = mtype
                name = cfg.moleculetype
            if renamed and renamed[-1][0] == name:
                renamed[-1] = (name, renamed[-1][1] + 1)
            else:
                renamed.append((name, 1))
            offset += size
    if len(starts) != len(restored):
        raise RepackError(f"Found {len(starts)} of {len(restored)} restored cholesterols in the topology")
    if cfg.moleculetype in definitions:
        raise RepackError(f"repack.moleculetype {cfg.moleculetype} is already a moleculetype in the run")
    source_mtype = next(iter(set(starts.values())))

    coverage, warnings = {}, []
    used = {name for name, _ in topology_result.molecules_written}
    classes = classify_resnames({a.resname for m in used for a in definitions[m].atoms[:1]})
    for mtype in sorted(used):
        definition = definitions[mtype]
        if classes.get(definition.atoms[0].resname) not in ("protein", "other"):
            continue
        info = restraint_coverage(definition)
        coverage[mtype] = info
        if info["missing_heavy_atoms"]:
            warnings.append(f"{mtype}: {len(info['missing_heavy_atoms'])} of {info['heavy_atoms']} heavy atoms "
                            "have no POSRES restraint and will move freely")
        other = [v for v in info["restraint_values"] if v not in ("POSRES_FC_BB", "POSRES_FC_SC")]
        if other:
            warnings.append(f"{mtype}: restrained with {other}, not the repack values")

    groups = ("System", "System")
    if ndx_groups and {"SOLU", "MEMB", "SOLV", "SOLU_MEMB"} <= set(ndx_groups):
        groups = ("SOLU MEMB SOLV", "SOLU_MEMB SOLV")

    # A self-contained folder: its own toppar/, topol.top, GRO, index and mdp files.
    bundle.mkdir(parents=True, exist_ok=True)
    toppar = bundle / "toppar"
    if toppar.exists():
        shutil.rmtree(toppar)
    shutil.copytree(topology_result.toppar_dir, toppar)
    itp_name = f"{cfg.moleculetype}.itp"
    (toppar / itp_name).write_text(make_restrained_cholesterol_itp(files[source_mtype], source_mtype,
                                                                   cfg.moleculetype))
    (bundle / "topol.top").write_text(_top_with(topology_result.top_path, f"toppar/{itp_name}", renamed))
    gro_name = Path(gro_path).name
    shutil.copy2(gro_path, bundle / gro_name)
    index_arg = ""
    if ndx_path:
        shutil.copy2(ndx_path, bundle / "index.ndx")
        index_arg = " -n index.ndx"
    first_md = next((i for i, s in enumerate(cfg.stages) if s.integrator != "steep"), None)
    for number, stage in enumerate(cfg.stages):
        (bundle / f"{stage.name}.mdp").write_text(_mdp(stage, cfg, number == first_md, groups))
    for old in ("run_repack.sh", "README.md", itp_name):  # earlier layout
        if (bundle / old).is_file():
            (bundle / old).unlink()

    md_ps = sum(s.nsteps * s.dt for s in cfg.stages if s.integrator != "steep")
    script = ["#!/bin/bash",
              "# Sample run script written by GMXtransplant. Adapt it to your machine, or submit",
              "# the same grompp/mdrun lines through your scheduler.",
              f"# Repacking equilibration, {md_ps / 1000:g} ns of MD: every heavy atom of the protein, ligands",
              f"# and the {len(restored)} restored cholesterol(s) ({cfg.moleculetype}) is restrained; the other "
              "lipids, water",
              "# and ions are free to fill the gaps left around the restored cholesterols.",
              "# Restraint (kJ/mol/nm^2) on protein, ligands and restored cholesterols per stage:"]
    script += [f"#   {s.name}: {s.restraint:g}" for s in cfg.stages]
    script += ["# Afterwards, continue from the last GRO with the run's own topol.top (same atoms).",
               "set -e", 'cd "$(dirname "$0")"', 'GMX=${GMX:-gmx}', 'MDRUN=${MDRUN:-"$GMX mdrun"}',
               f"previous={gro_name}"]
    for stage in cfg.stages:
        script += [f'echo "==> {stage.name}"',
                   f"$GMX grompp -f {stage.name}.mdp -o {stage.name}.tpr -c $previous -r {gro_name} "
                   f"-p topol.top{index_arg}",
                   f"$MDRUN -v -deffnm {stage.name}", f"previous={stage.name}.gro"]
    (bundle / "sample_run.sh").write_text("\n".join(script) + "\n")
    os.chmod(bundle / "sample_run.sh", 0o755)
    return {"directory": str(bundle), "restored_cholesterols": len(restored),
            "restored_moleculetype": cfg.moleculetype, "molecules": renamed,
            "stages": [{"name": s.name, "integrator": s.integrator, "nsteps": s.nsteps, "dt": s.dt,
                        "restraint": s.restraint} for s in cfg.stages],
            "md_length_ps": md_ps, "lipid_restraint": cfg.lipid_restraint,
            "temperature": cfg.temperature, "coupling_groups": list(groups),
            "solute_restraint_coverage": coverage, "warnings": warnings}
