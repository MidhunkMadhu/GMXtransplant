"""Point mutations in a prepared system.

Each mutated residue keeps its backbone (N, HN, CA, HA, C, O) and, when both
residues have one, its CB. The new side chain is built from the CHARMM36 internal
coordinates; its rotamer is chosen from staggered (and, for planar groups,
perpendicular and in-plane) chi angles by the smallest steric clash with the
rest of the system, water and ions excluded. The protein ITP gets the new
residue's atoms and bonded terms (see protein_edit), water and ions overlapping
the side chain are removed, lipids too when asked, and counterions restore the
net charge. Terminal, capped and disulfide-bonded residues cannot be mutated.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from itertools import product

import numpy as np
from MDAnalysis.lib.distances import capped_distance

from config import ConfigError
from equilibration import EquilibrationSpec
from protein_edit import replace_residue, residue_state
from residue_edit import report_lines, run_changes
from system_edit import load_section, positive, run_edit
import charmm36_library as c36

TARGETS = ("ALA", "ARG", "ASN", "ASP", "ASPP", "CYS", "GLN", "GLU", "GLUP", "GLY", "HSD", "HSE", "HSP", "ILE",
           "LEU", "LYS", "LSN", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL")
STAGGERED = (-60.0, 180.0, 60.0)
PLANAR = (90.0, -90.0, 0.0, 180.0)
SP2_TYPES = {"CA", "CPH1", "CPH2", "CPT", "CY", "CC", "CD", "C", "NC2", "CAI", "CN"}
MAX_ROTAMERS = 400
CLASH = 3.3   # A; heavy-atom pairs closer than this are penalised (squared overlap)


def _no_equilibration():
    return EquilibrationSpec(enabled=False)


@dataclass
class Mutation:
    molecule: str = "PROA"     # protein moleculetype
    resid: int = 1             # residue number as in the ITP
    to: str = "ALA"            # CHARMM36 residue name


def _default_mutations():
    return [Mutation()]


@dataclass
class MutateSpec:
    host: str
    host_gro: str | None = None
    mutations: list[Mutation] = field(default_factory=_default_mutations)
    water_clash_distance: float = 2.5    # A, new atoms to water/ion heavy atoms: those are removed
    lipid_clash_distance: float = 2.0    # A, new heavy atoms to lipid heavy atoms
    remove_clashing_lipids: bool = True  # remove those lipids whole (else only warn)
    neutralize: bool = True
    target_net_charge: float = 0.0
    seed: int = 42
    environment_overrides: dict = field(default_factory=dict)
    output_dir: str = "mutate_output"
    equilibration: EquilibrationSpec = field(default_factory=_no_equilibration)


def load_config(path, check_paths=True, output_root=None):
    cfg = load_section(path, "mutate", MutateSpec, check_paths, output_root)
    positive(cfg, "mutate", "water_clash_distance")
    positive(cfg, "mutate", "lipid_clash_distance", allow_zero=True)
    if not cfg.mutations:
        raise ConfigError("mutate.mutations must list at least one mutation")
    for number, mutation in enumerate(cfg.mutations):
        if mutation.to not in TARGETS:
            raise ConfigError(f"mutate.mutations[{number}].to must be one of {', '.join(TARGETS)}")
    return cfg


def _ring_bonds(template):
    """Bonds of the template's heavy-atom graph that lie in a ring."""
    graph = defaultdict(set)
    heavy = {a[0] for a in template["atoms"] if not a[1].upper().startswith("H")}
    for a, b in template["bonds"]:
        if a in heavy and b in heavy:
            graph[a].add(b)
            graph[b].add(a)
    rings = set()
    for a in graph:
        for b in graph[a]:
            # (a, b) is in a ring when b is reachable from a without that bond.
            seen, stack = {a}, [n for n in graph[a] if n != b]
            while stack:
                node = stack.pop()
                if node == b:
                    rings.add(frozenset((a, b)))
                    break
                if node not in seen:
                    seen.add(node)
                    stack.extend(graph[node] - seen)
    return rings


def chi_options(name):
    """{IC line index: candidate dihedrals} for the side-chain chi angles of a template."""
    if name in ("PRO", "GLY", "ALA"):
        return {}
    template = c36.residue(name)
    types = {a[0]: a[1] for a in template["atoms"]}
    rings = _ring_bonds(template)
    backbone = {"N", "CA", "C", "O"}
    placed, options = set(), {}
    for n, ic in enumerate(template["ic"]):
        i, j, k, l = ic[:4]
        if k.startswith("*") or any(x.startswith(("H", "+", "-")) for x in (i, j, k, l)):
            continue
        # A chi: proper dihedral about a side-chain bond (CA-CB, CB-CG, ...), outside rings.
        if k in backbone or j in ("N", "C", "O") or l in backbone or l in placed:
            continue
        if frozenset((j, k)) in rings or (types.get(j) in SP2_TYPES and types.get(k) in SP2_TYPES):
            continue  # ring bonds and conjugated (partial double) bonds stay as built
        placed.add(l)
        options[n] = PLANAR if types.get(k) in SP2_TYPES else STAGGERED
    return options


def best_rotamer(itp, resnr, target, positions, keep, environment, box, host_files):
    """Try chi combinations; return (ic table, report)."""
    template = c36.residue(target)
    options = chi_options(target)
    if not options:
        return None, None
    keys = sorted(options)
    combos = list(product(*(options[k] for k in keys)))
    if len(combos) > MAX_ROTAMERS:
        combos = combos[:MAX_ROTAMERS]
    side_heavy = [a[0] for a in template["atoms"]
                  if not a[1].upper().startswith("H") and a[0] not in ("N", "CA", "C", "O", "CB")]
    best = None
    for combo in combos:
        table = [list(ic) for ic in template["ic"]]
        for key, value in zip(keys, combo):
            table[key][6] = value
        coords, near = replace_residue(itp, resnr, target, positions, keep=keep, host_parameter_files=host_files,
                                       ic_table=table, coordinates_only=True)
        if best is None:
            own = np.array([i for i in range(len(itp.atoms))
                            if i not in near and not itp.atoms[i].type.upper().startswith("H")])
            context = np.concatenate([positions[own], environment]) if len(environment) else positions[own]
        points = np.array([coords[n] for n in side_heavy])
        _, distances = capped_distance(points, context, CLASH, box=box)
        score = float(np.sum((CLASH - distances) ** 2))
        if best is None or score < best[0] - 1e-9:
            best = (score, combo, table)
    score, combo, table = best
    return table, {"chi": [round(v) for v in combo], "score": round(score, 3), "tried": len(combos)}


def run(config_path, dry_run=False, output_root=None):
    cfg = load_config(config_path, check_paths=not dry_run, output_root=output_root)
    host_files = []
    box_holder = {}

    def apply(itp, positions, resid, target, environment, mol):
        names = [a.name for a in itp.atoms if a.resnr == resid]
        resname = next(a.resname for a in itp.atoms if a.resnr == resid)
        state = residue_state(resname, names)
        if state is None:
            raise ConfigError(f"{mol.name} {resname} {resid} matches no CHARMM36 residue (terminal, capped or "
                              "patched); it cannot be mutated")
        new_names = {a[0] for a in c36.residue(target)["atoms"]}
        keep = {"N", "HN", "CA", "C", "O"} | ({"HA"} & new_names) | ({"CB"} & set(names) & new_names)
        table, rotamer = best_rotamer(itp, resid, target, positions, keep, environment, box_holder["box"],
                                      [str(p) for p in host_files])
        new_positions, change = replace_residue(itp, resid, target, positions, keep=keep,
                                                host_parameter_files=[str(p) for p in host_files], ic_table=table)
        change["to_label"] = ""
        if rotamer:
            change["rotamer"] = rotamer
        return new_positions, change

    def build(stage, host, report):
        host_files.extend(host.files)
        box_holder["box"] = host.universe.dimensions
        changes = [(m.molecule, m.resid, m.to) for m in cfg.mutations]
        run_changes(stage, host, report, cfg, "mutate", changes, apply,
                    purpose="Equilibration of the mutated protein", remove_lipids=cfg.remove_clashing_lipids)

    return run_edit("mutate", cfg, "POINT MUTATIONS", build, "point mutations",
                    report_lines("mutate", "MUTATIONS"), dry_run=dry_run)
