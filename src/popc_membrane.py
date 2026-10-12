"""Replace the membrane of a prepared system with a well-equilibrated pure POPC bilayer.

The host is one prepared CHARMM-GUI-style folder (topol.top, toppar/, GRO) with
POPC among its lipids. Its proteins and ligands never move, so the membrane
placement CHARMM-GUI derived from OPM is kept. The host bilayer is measured
from its lipid phosphorus atoms: each leaflet must hold at least three POPC
whose P atoms lie at least 20 A apart laterally, and the plane fitted through
each leaflet's P atoms must be flat and normal to z.

Every host lipid and sterol is then removed and replaced by POPC from an
equilibrated reference bilayer (by default the 500 ns frame of the CHARMM36
POPC/150 mM NaCl simulation of Zenodo 10.5281/zenodo.1010083). The reference
is repeated periodically and cut to the host box, its midplane on the host
midplane. Reference lipids that overlap the protein or ligands, or a periodic
copy of another reference lipid at the box edge, are removed whole. Host water
and ions that overlap the new lipids are removed, counterions restore the net
charge, and a GROMACS equilibration folder lets the lipids close the gaps
around the protein before production.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from importlib.resources import files
from pathlib import Path
import json
import os
import shutil
import tempfile
import warnings

import MDAnalysis as mda
from MDAnalysis.lib.distances import capped_distance, self_capped_distance
from MDAnalysis.lib.mdamath import triclinic_vectors
import numpy as np
import yaml

from addbinder import _unwrap_sequential
from charmprot import (Molecule, inspect_system, write_topology, _heavy, _indices,
                       _publish_directory)
from config import ConfigError, _UniqueKeySafeLoader, anchor_config_paths
from ion_selection import select_counterions
from output import write_outputs
from topology import audit_final_topology, TopologyError

LIPID = "POPC"
REMOVED_CATEGORIES = ("lipid", "sterol", "detergent")
SOLVENT_CATEGORIES = ("water", "ion", "solvent")
WATER_MOLARITY = 55.51   # mol/L of pure water: salt concentration = pairs / waters * 55.51
DEFAULT_REFERENCE = "popc512_charmm36_nacl150mM_310K_500ns.gro.gz"
REFERENCE_SOURCE = ("POPC bilayer, 512 lipids, CHARMM36, 150 mM NaCl, 310 K, frame at 500 ns. "
                    "H. Martinez-Seara, Zenodo, doi:10.5281/zenodo.1010083 (CC BY 4.0); "
                    "lipids only, water and ions removed")
GENERATED = {"step5_input.pdb", "step5_input.gro", "topol.top", "toppar", "index.ndx",
             "popc_report.json", "popc_report.txt", "equilibration",
             "visualization", "view.pml", "view.vmd", "view.pse"}


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

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


def _default_stages():
    # The CHARMM-GUI Membrane Builder schedule: the protein is released slowly, so it
    # keeps its OPM placement while the new lipids pack around it; 1.875 ns of MD.
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
    stages: list = field(default_factory=_default_stages)


@dataclass
class PopcSpec:
    host: str
    host_gro: str | None = None
    # Equilibrated bilayer to take the POPC from; empty: the bundled Zenodo 1010083 frame.
    reference_gro: str | None = None
    # Midplane z (A, host GRO frame); empty: midway between the host's P planes.
    membrane_center_z: float | None = None
    anchor_count: int = 3                # POPC per leaflet whose P atoms are ...
    anchor_separation: float = 20.0      # ... at least this far apart laterally (A)
    max_tilt: float = 5.0                # degrees between a leaflet's P plane and the xy plane
    protein_clash_distance: float = 3.0  # A, POPC heavy atom to protein/ligand heavy atom
    # Lateral registrations of the reference tried; the one keeping the most POPC is used.
    lateral_search: int = 16
    # Reference water hydrates the new bilayer up to this far beyond each P plane (A)
    # and replaces the host water there, except host water near the protein.
    water_layer: float = 10.0
    protein_hydration_distance: float = 4.0  # A; host water/ions this close to the protein stay
    reference_contact_distance: float = 2.0  # A; reference copies meeting at the box edge closer than this
    water_clash_distance: float = 2.5    # A, host water/ion heavy atom to POPC or reference water
    neutralize: bool = True
    target_net_charge: float = 0.0
    # Bulk waters outside the new hydration layer become salt pairs until the host's
    # salt concentration is reached again (the reference water brings no ions).
    restore_salt: bool = True
    ion_exclusion_distance: float = 10.0
    lipid_ion_exclusion_distance: float = 5.0
    environment_overrides: dict[str, str] = field(default_factory=dict)
    output_dir: str = "popc_output"
    equilibration: EquilibrationSpec = field(default_factory=EquilibrationSpec)


def _number(value, key, positive=True):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not np.isfinite(value):
        raise ConfigError(f"popc.{key} must be a number")
    if positive and value <= 0:
        raise ConfigError(f"popc.{key} must be greater than 0")
    return float(value)


def _parse_equilibration(raw):
    if raw is None:
        return EquilibrationSpec()
    if not isinstance(raw, dict):
        raise ConfigError("popc.equilibration must be a mapping")
    unknown = set(raw) - set(EquilibrationSpec.__dataclass_fields__)
    if unknown:
        raise ConfigError(f"Unknown popc.equilibration options: {sorted(unknown)}")
    stages = raw.get("stages")
    spec = EquilibrationSpec(**{k: v for k, v in raw.items() if k != "stages"})
    if stages is not None:
        if not isinstance(stages, list) or not stages:
            raise ConfigError("popc.equilibration.stages must be a nonempty list")
        parsed = []
        for number, stage in enumerate(stages):
            where = f"popc.equilibration.stages[{number}]"
            if isinstance(stage, EquilibrationStage):
                stage = asdict(stage)
            if not isinstance(stage, dict):
                raise ConfigError(f"{where} must be a mapping")
            unknown = set(stage) - set(EquilibrationStage.__dataclass_fields__)
            if unknown:
                raise ConfigError(f"Unknown {where} options: {sorted(unknown)}")
            parsed.append(EquilibrationStage(**stage))
        spec.stages = parsed
    if type(spec.enabled) is not bool:
        raise ConfigError("popc.equilibration.enabled must be true or false")
    if (not isinstance(spec.output_dir, str) or not spec.output_dir
            or Path(spec.output_dir).name != spec.output_dir or spec.output_dir in (".", "..")):
        raise ConfigError("popc.equilibration.output_dir must be a simple folder name")
    spec.temperature = _number(spec.temperature, "equilibration.temperature")
    names = set()
    for number, stage in enumerate(spec.stages):
        where = f"equilibration.stages[{number}]"
        if not isinstance(stage.name, str) or not stage.name or Path(stage.name).name != stage.name:
            raise ConfigError(f"popc.{where}.name must be a simple file name")
        if stage.name in names:
            raise ConfigError(f"popc.{where}.name {stage.name} is used twice")
        names.add(stage.name)
        if stage.integrator not in ("steep", "nvt", "npt"):
            raise ConfigError(f"popc.{where}.integrator must be steep, nvt or npt")
        if type(stage.nsteps) is not int or stage.nsteps <= 0:
            raise ConfigError(f"popc.{where}.nsteps must be a whole number > 0")
        for key in ("backbone", "sidechain", "lipid"):
            setattr(stage, key, _number(getattr(stage, key), f"{where}.{key}", positive=False))
            if getattr(stage, key) < 0:
                raise ConfigError(f"popc.{where}.{key} must not be negative")
        stage.dt = _number(stage.dt, f"{where}.dt", positive=False)
        if stage.integrator != "steep" and not 0 < stage.dt <= 0.004:
            raise ConfigError(f"popc.{where}.dt must be between 0 and 0.004 ps")
    return spec


def load_popc_config(path, check_paths=True, output_root=None):
    try:
        raw = yaml.load(Path(path).read_text(), Loader=_UniqueKeySafeLoader)
    except (OSError, yaml.YAMLError) as exc:
        raise ConfigError(str(exc)) from exc
    if not isinstance(raw, dict) or not isinstance(raw.get("popc"), dict):
        raise ConfigError("popc YAML must contain a popc: mapping")
    from config import _resolve_path_references
    raw = _resolve_path_references(raw)
    raw.pop("paths", None)
    if set(raw) - {"popc"}:
        raise ConfigError(f"Unknown popc sections: {sorted(set(raw) - {'popc'})}")
    raw = dict(raw["popc"])
    unknown = set(raw) - set(PopcSpec.__dataclass_fields__)
    if unknown:
        raise ConfigError(f"Unknown popc options: {sorted(unknown)}")
    if not isinstance(raw.get("host"), str) or not raw["host"].strip():
        raise ConfigError("popc.host must be a nonempty folder path")
    for key in ("host_gro", "reference_gro", "membrane_center_z"):
        if raw.get(key) == "":
            raw[key] = None  # an empty optional field in the GUI means "not set"
    for key in ("host_gro", "reference_gro"):
        value = raw.get(key)
        if value is not None and (not isinstance(value, str) or not value.strip()):
            raise ConfigError(f"popc.{key} must be a file path")
    raw.setdefault("output_dir", PopcSpec.__dataclass_fields__["output_dir"].default)
    equilibration = raw.pop("equilibration", None)
    raw = anchor_config_paths({"popc": raw}, config_path=path, output_root=output_root)["popc"]
    cfg = PopcSpec(**raw)
    cfg.equilibration = _parse_equilibration(equilibration)
    if cfg.membrane_center_z is not None:
        cfg.membrane_center_z = _number(cfg.membrane_center_z, "membrane_center_z", positive=False)
    if type(cfg.anchor_count) is not int or cfg.anchor_count < 3:
        raise ConfigError("popc.anchor_count must be a whole number, 3 or more (three points define a plane)")
    if type(cfg.lateral_search) is not int or not 1 <= cfg.lateral_search <= 400:
        raise ConfigError("popc.lateral_search must be a whole number from 1 to 400")
    for key in ("anchor_separation", "max_tilt", "protein_clash_distance", "reference_contact_distance",
                "water_layer", "protein_hydration_distance",
                "water_clash_distance", "ion_exclusion_distance", "lipid_ion_exclusion_distance"):
        setattr(cfg, key, _number(getattr(cfg, key), key))
    if cfg.max_tilt >= 45:
        raise ConfigError("popc.max_tilt must be below 45 degrees")
    cfg.target_net_charge = _number(cfg.target_net_charge, "target_net_charge", positive=False)
    for key in ("neutralize", "restore_salt"):
        if type(getattr(cfg, key)) is not bool:
            raise ConfigError(f"popc.{key} must be true or false")
    if not isinstance(cfg.environment_overrides, dict) or any(
            not isinstance(k, str) or v not in ("water", "ion", "lipid", "sterol", "detergent", "solvent")
            for k, v in cfg.environment_overrides.items()):
        raise ConfigError("popc.environment_overrides must map molecule names to environment categories")
    output, host = Path(cfg.output_dir), Path(cfg.host)
    if output == host or host in output.parents or output in host.parents:
        raise ConfigError("popc.output_dir must be separate from the host folder")
    if check_paths:
        if not host.is_dir():
            raise ConfigError(f"host folder does not exist: {host}")
        if not (host / "topol.top").is_file() or not (host / "toppar").is_dir():
            raise ConfigError(f"host must contain topol.top and toppar/: {host}")
        for key in ("host_gro", "reference_gro"):
            value = getattr(cfg, key)
            if value and not Path(value).is_file():
                raise ConfigError(f"popc.{key} does not exist: {value}")
    return cfg


def default_reference_path():
    return Path(str(files("gmxtransplant").joinpath("data", DEFAULT_REFERENCE)))


# ---------------------------------------------------------------------------
# Membrane geometry
# ---------------------------------------------------------------------------

def _cell(dimensions):
    """xy cell vectors (2x2, rows a and b) and the box height; c must lie along z."""
    vectors = triclinic_vectors(dimensions)
    if abs(vectors[2, 0]) > 1e-3 or abs(vectors[2, 1]) > 1e-3:
        raise ConfigError("The box's third vector is not along z; a membrane box needs alpha = beta = 90 degrees")
    return vectors[:2, :2].astype(float), float(vectors[2, 2])


def _midplane_shift(z, height):
    """Shift that puts the bilayer centre at height / 2, from the z of the lipid tail ends.

    The tail ends gather at the bilayer centre, so their circular mean is defined however
    thin the water layer is; the mean of the two P planes is not when they lie half a box apart.
    """
    angle = 2 * np.pi * z / height
    centre = np.mod(np.arctan2(np.sin(angle).mean(), np.cos(angle).mean()) * height / (2 * np.pi), height)
    return height / 2 - centre


def _min_image_xy(delta, cell):
    """Minimum-image xy displacements in a (possibly hexagonal) cell."""
    inverse = np.linalg.inv(cell)
    fractional = delta @ inverse
    fractional -= np.round(fractional)
    out = fractional @ cell
    # A hexagonal cell can need a neighbouring image as well.
    best = out.copy()
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            trial = out + i * cell[0] + j * cell[1]
            shorter = np.linalg.norm(trial, axis=-1) < np.linalg.norm(best, axis=-1)
            best[shorter] = trial[shorter]
    return best


def _tail_end(whole, p_index):
    """z of the lipid atom farthest from its P atom: a tail end, near the bilayer centre."""
    return whole[np.argmax(np.linalg.norm(whole - whole[p_index], axis=1)), 2]


def _plane(xy, z):
    """Least-squares plane z = a x + b y + c: (a, b, c), tilt in degrees, residual spread."""
    design = np.column_stack([xy, np.ones(len(xy))])
    coefficients, *_ = np.linalg.lstsq(design, z, rcond=None)
    residual = z - design @ coefficients
    tilt = float(np.degrees(np.arctan(np.hypot(coefficients[0], coefficients[1]))))
    return coefficients, tilt, float(residual.std())


def _anchors(xy, count, separation, cell):
    """Farthest-point choice of `count` P atoms; returns (indices, smallest pairwise lateral distance)."""
    centroid = xy.mean(axis=0)
    chosen = [int(np.argmin(np.linalg.norm(_min_image_xy(xy - centroid, cell), axis=1)))]
    nearest = np.linalg.norm(_min_image_xy(xy - xy[chosen[0]], cell), axis=1)
    while len(chosen) < count:
        pick = int(np.argmax(nearest))
        if nearest[pick] <= 0:
            break
        chosen.append(pick)
        nearest = np.minimum(nearest, np.linalg.norm(_min_image_xy(xy - xy[pick], cell), axis=1))
    pairs = [np.linalg.norm(_min_image_xy(xy[i] - xy[j], cell))
             for n, i in enumerate(chosen) for j in chosen[n + 1:]]
    return chosen, (float(min(pairs)) if pairs else 0.0)


@dataclass
class HostMembrane:
    shift: float          # added to z to put the midplane at height / 2
    height: float
    midplane: float       # GRO frame, A
    planes: dict          # leaflet -> report
    anchors: dict         # leaflet -> list of anchor reports
    leaflet_counts: dict  # leaflet -> Counter of host lipid moleculetypes


def measure_host(host, cfg):
    universe = host.universe
    cell, height = _cell(universe.dimensions)
    membrane = [m for m in host.molecules if m.category in REMOVED_CATEGORIES]
    popc = [m for m in membrane if m.name == LIPID]
    if not popc:
        raise ConfigError(f"The host has no {LIPID} moleculetype; this mode replaces a membrane that contains "
                          f"{LIPID} (lipid moleculetypes found: {sorted({m.name for m in membrane}) or 'none'})")
    if not any(a.atom_name == "P" for a in popc[0].definition.atoms):
        raise ConfigError(f"The host {LIPID} has no P atom to measure the headgroup planes")
    lengths = universe.dimensions[:3].astype(float)
    heads, others, tails = [], [], []   # (molecule, P position); molecules without P (sterols ...)
    for mol in membrane:
        names = [a.atom_name for a in mol.definition.atoms]
        if "P" in names:
            whole = _unwrap_sequential(universe.atoms[mol.start:mol.stop].positions.astype(float), lengths)
            p_index = names.index("P")
            heads.append((mol, universe.atoms[mol.start + p_index].position.astype(float)))
            tails.append(_tail_end(whole, p_index) - whole[p_index, 2] + heads[-1][1][2])
        else:
            others.append(mol)
    if len(heads) < 6:
        raise ConfigError("Fewer than six phospholipids in the host; cannot measure the bilayer")
    positions = np.array([p for _, p in heads])
    shift = float(_midplane_shift(np.array(tails), height))
    zc = np.mod(positions[:, 2] + shift, height)
    leaflet = np.where(zc > height / 2, "upper", "lower")
    planes, anchors, counts = {}, {}, {}
    for side in ("upper", "lower"):
        on_side = leaflet == side
        if on_side.sum() < 3:
            raise ConfigError(f"Found {int(on_side.sum())} phosphorus atom(s) in the {side} leaflet; "
                              "the membrane normal must lie along the box z axis")
        (a, b, c), tilt, spread = _plane(positions[on_side, :2], zc[on_side])
        if tilt > cfg.max_tilt:
            raise ConfigError(f"The {side} leaflet's P plane is tilted {tilt:.1f} degrees from the xy plane "
                              f"(max_tilt {cfg.max_tilt:g}); the membrane must be flat and normal to z")
        planes[side] = {"mean_z": round(float(np.mean(zc[on_side]) - shift), 3), "tilt_degrees": round(tilt, 3),
                        "spread_angstrom": round(spread, 3), "phosphorus_atoms": int(on_side.sum())}
        indices = [i for i in np.flatnonzero(on_side) if heads[i][0].name == LIPID]
        if len(indices) < cfg.anchor_count:
            raise ConfigError(f"The {side} leaflet holds {len(indices)} {LIPID}; at least {cfg.anchor_count} "
                              f"are needed, {cfg.anchor_separation:g} A apart")
        chosen, spacing = _anchors(positions[indices, :2], cfg.anchor_count, cfg.anchor_separation, cell)
        if spacing + 1e-6 < cfg.anchor_separation:
            raise ConfigError(f"The {side} leaflet has no {cfg.anchor_count} {LIPID} whose P atoms are "
                              f"{cfg.anchor_separation:g} A apart (best {spacing:.1f} A)")
        triangle = positions[[indices[i] for i in chosen[:3]], :2]
        legs = _min_image_xy(triangle[1:] - triangle[0], cell)
        area = 0.5 * abs(legs[0, 0] * legs[1, 1] - legs[0, 1] * legs[1, 0])
        if area < 0.1 * cfg.anchor_separation ** 2:
            raise ConfigError(f"The {side} leaflet's {LIPID} anchors lie on one line; no plane is defined")
        anchors[side] = [{"moleculetype": heads[indices[i]][0].name, "occurrence": heads[indices[i]][0].occurrence,
                          "resid": int(universe.atoms[heads[indices[i]][0].start].resid),
                          "p_position_angstrom": [round(float(v), 3) for v in positions[indices[i]]]}
                         for i in chosen]
        anchors[side + "_min_separation_angstrom"] = round(spacing, 3)
        counts[side] = Counter(heads[i][0].name for i in np.flatnonzero(on_side))
    for mol in others:
        whole = _unwrap_sequential(universe.atoms[mol.start:mol.stop].positions.astype(float),
                                   universe.dimensions[:3].astype(float))
        centre = np.mod(whole[:, 2].mean() + shift, height)
        counts["upper" if centre > height / 2 else "lower"][mol.name] += 1
    midplane_centred = (planes["upper"]["mean_z"] + planes["lower"]["mean_z"]) / 2 + shift
    if cfg.membrane_center_z is not None:
        midplane_centred = cfg.membrane_center_z + shift
    thickness = planes["upper"]["mean_z"] - planes["lower"]["mean_z"]
    if not 20.0 <= thickness <= 70.0:
        raise ConfigError(f"Headgroup planes are {thickness:.1f} A apart; expected a flat bilayer normal to z")
    return HostMembrane(shift, height, float(np.mod(midplane_centred - shift, height)), planes, anchors, counts)


# ---------------------------------------------------------------------------
# Reference bilayer
# ---------------------------------------------------------------------------

@dataclass
class ReferenceBilayer:
    path: str
    positions: list       # per lipid: whole positions (n_atoms, 3), midplane at z = 0, ITP atom order
    leaflets: list        # per lipid: "upper" | "lower"
    heads: np.ndarray     # per lipid: P position, inside the reference cell
    waters: np.ndarray    # (n, 3, 3) whole waters, O first, O inside the reference cell
    cell: np.ndarray
    height: float
    report: dict


def _reference_waters(universe, lengths, shift, cell, height):
    """Three-site waters, O then H H, whole, wrapped by their O like the lipids."""
    from charmm_environment import environment_name
    waters = []
    for residue in universe.residues:
        known = environment_name(str(residue.resname))
        if not known or known[0] != "water":
            continue
        names = [str(n) for n in residue.atoms.names]
        oxygen = [i for i, n in enumerate(names) if n.upper().startswith("O")]
        hydrogen = [i for i, n in enumerate(names) if n.upper().startswith("H")]
        if len(names) != 3 or len(oxygen) != 1 or len(hydrogen) != 2:
            raise ConfigError(f"Reference water {residue.resname} {residue.resid} is not a three-site water")
        waters.append(_unwrap_sequential(residue.atoms.positions.astype(float)[oxygen + hydrogen], lengths))
    if not waters:
        return np.zeros((0, 3, 3))
    waters = np.array(waters)
    waters[:, :, 2] += shift
    waters[:, :, 2] -= (height * np.floor(waters[:, 0, 2] / height) + height / 2)[:, None]
    whole = np.floor(waters[:, 0, :2] @ np.linalg.inv(cell))
    waters[:, :, :2] -= (whole @ cell)[:, None, :]
    return waters


def load_reference(path, definition):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        universe = mda.Universe(str(path))
    if universe.dimensions is None:
        raise ConfigError(f"The reference {path} has no box")
    if np.any(np.abs(universe.dimensions[3:] - 90) > 1e-3):
        raise ConfigError("The reference bilayer needs a rectangular box")
    cell, height = _cell(universe.dimensions)
    lengths = universe.dimensions[:3].astype(float)
    expected = [a.atom_name for a in definition.atoms]
    residues = [r for r in universe.residues if str(r.resname).upper() == LIPID]
    if len(residues) < 20:
        raise ConfigError(f"The reference {path} holds {len(residues)} {LIPID} residues; expected a bilayer")
    lipids, reordered = [], 0
    for residue in residues:
        names = [str(n) for n in residue.atoms.names]
        order = slice(None)
        if names != expected:
            if sorted(names) != sorted(expected) or len(set(names)) != len(names):
                raise ConfigError(f"Reference {LIPID} {residue.resid} atoms do not match the host {LIPID}.itp")
            order = [names.index(n) for n in expected]
            reordered += 1
        lipids.append(_unwrap_sequential(residue.atoms.positions.astype(float)[order], lengths))
    p_index = expected.index("P")
    shift = _midplane_shift(np.array([_tail_end(lipid, p_index) for lipid in lipids]), height)
    # Whole lipids move by box vectors only: P z into [-height/2, height/2) about the
    # midplane, P xy into the reference cell.
    out, leaflets = [], []
    inverse = np.linalg.inv(cell)
    for lipid in lipids:
        lipid = lipid.copy()
        lipid[:, 2] += shift - height * np.floor((lipid[p_index, 2] + shift) / height) - height / 2
        lipid[:, :2] -= np.floor(lipid[p_index, :2] @ inverse) @ cell
        out.append(lipid)
        leaflets.append("upper" if lipid[p_index, 2] > 0 else "lower")
    heads = np.array([lipid[p_index] for lipid in out])
    planes = {}
    for side in ("upper", "lower"):
        mask = np.array(leaflets) == side
        _, tilt, spread = _plane(heads[mask, :2], heads[mask, 2])
        if tilt > 5.0:
            raise ConfigError(f"The reference {side} leaflet is tilted {tilt:.1f} degrees; expected a flat bilayer")
        planes[side] = {"lipids": int(mask.sum()), "mean_z": round(float(heads[mask, 2].mean()), 3),
                        "spread_angstrom": round(spread, 3)}
    waters = _reference_waters(universe, lengths, shift, cell, height)
    area = abs(np.linalg.det(cell))
    report = {"path": str(path), "lipids": len(out), "waters": len(waters),
              "box_angstrom": [round(float(v), 3) for v in lengths],
              "phosphate_thickness_angstrom": round(planes["upper"]["mean_z"] - planes["lower"]["mean_z"], 3),
              "area_per_lipid_angstrom2": {side: round(area / planes[side]["lipids"], 2) for side in planes},
              "leaflets": planes, "residues_reordered_to_itp": reordered}
    return ReferenceBilayer(str(path), out, leaflets, heads, waters, cell, height, report)


def tile(points, ref_cell, host_cell, origin):
    """Periodic copies of reference points that fall inside the host cell: [(index, xy shift)]."""
    inverse = np.linalg.inv(host_cell)
    corners = np.array([[0, 0], host_cell[0], host_cell[1], host_cell[0] + host_cell[1]]) - origin
    span = corners @ np.linalg.inv(ref_cell)
    placed = []
    for i in range(int(np.floor(span[:, 0].min())) - 1, int(np.ceil(span[:, 0].max())) + 1):
        for j in range(int(np.floor(span[:, 1].min())) - 1, int(np.ceil(span[:, 1].max())) + 1):
            offset = i * ref_cell[0] + j * ref_cell[1] + origin
            fractional = (points + offset) @ inverse
            inside = np.all((fractional >= 0) & (fractional < 1), axis=1)
            placed += [(int(k), offset) for k in np.flatnonzero(inside)]
    return placed


# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------

def _heavy_mask(definition):
    # CHARMM hydrogen types start with H; lone pairs (LP*) are massless sites.
    return np.array([not a.atom_type.upper().startswith(("H", "LP")) for a in definition.atoms])


def _owner_map(blocks):
    """Concatenate per-molecule position blocks; return positions and the molecule of each row."""
    if not blocks:
        return np.zeros((0, 3)), np.zeros(0, dtype=int)
    return (np.concatenate(blocks),
            np.concatenate([np.full(len(b), i, dtype=int) for i, b in enumerate(blocks)]))


def _touching(blocks, others, cutoff, box):
    """Indices of the blocks with an atom within cutoff of any of `others`."""
    positions, owner = _owner_map(blocks)
    if not len(positions) or not len(others):
        return set()
    pairs = capped_distance(positions, others, cutoff, box=box, return_distances=False)
    return {int(owner[a]) for a, _ in pairs}


def _resolve_contacts(blocks, cutoff, box):
    """Fewest blocks to drop so no two blocks come within cutoff (copies meeting at the box edge)."""
    positions, owner = _owner_map(blocks)
    if not len(positions):
        return set()
    edges = {}
    for a, b in self_capped_distance(positions, cutoff, box=box, return_distances=False):
        i, j = int(owner[a]), int(owner[b])
        if i != j:
            edges.setdefault(i, set()).add(j)
            edges.setdefault(j, set()).add(i)
    dropped = set()
    while any(edges.values()):
        worst = max(edges, key=lambda i: (len(edges[i]), i))
        dropped.add(worst)
        for other in edges.pop(worst):
            edges[other].discard(worst)
    return dropped


def _lipids_at(cfg, reference, host_cell, frame, origin, heavy, protected_heavy, box):
    candidates = []
    for k, offset in tile(reference.heads[:, :2], reference.cell, host_cell, origin):
        lipid = reference.positions[k].copy()
        lipid[:, :2] += offset
        lipid[:, 2] += frame.midplane
        candidates.append((k, reference.leaflets[k], lipid))
    clash = _touching([lipid[heavy] for _, _, lipid in candidates], protected_heavy,
                      cfg.protein_clash_distance, box)
    keep = [c for i, c in enumerate(candidates) if i not in clash]
    edge = _resolve_contacts([lipid[heavy] for _, _, lipid in keep], cfg.reference_contact_distance, box)
    kept = [c for i, c in enumerate(keep) if i not in edge]
    report = {"tiled": dict(Counter(side for _, side, _ in candidates)),
              "removed_protein_overlap": dict(Counter(candidates[i][1] for i in clash)),
              "removed_box_edge_overlap": dict(Counter(keep[i][1] for i in edge))}
    return kept, report


def place_popc(cfg, host, frame, reference, definition, protected):
    """Reference POPC on the host midplane, cut to the host box; drop overlaps. Returns lipids and a report."""
    box = host.universe.dimensions
    host_cell, _ = _cell(box)
    heavy = _heavy_mask(definition)
    protected_heavy = _heavy(host.universe.atoms[_indices(protected)]).positions
    # Registrations on a square grid spanning one reference cell; ties keep the earliest.
    side = int(np.ceil(np.sqrt(cfg.lateral_search)))
    origins = [(i / side) * reference.cell[0] + (j / side) * reference.cell[1]
               for i in range(side) for j in range(side)][:cfg.lateral_search]
    best = None
    for number, origin in enumerate(origins):
        kept, report = _lipids_at(cfg, reference, host_cell, frame, origin, heavy, protected_heavy, box)
        if best is None or len(kept) > len(best[0]):
            best = (kept, report, origin, number)
    kept, report, origin, number = best
    # Upper leaflet first, each in reference order: a tidy, reproducible POPC block.
    kept.sort(key=lambda item: (item[1] != "upper", item[0]))
    report["placed"] = dict(Counter(side for _, side, _ in kept))
    report["lateral_search"] = {"registrations_tried": len(origins), "chosen": number,
                                "origin_angstrom": [round(float(v), 3) for v in origin]}
    return kept, origin, report


def _wrapped_dz(z, midplane, height):
    dz = z - midplane
    return dz - height * np.round(dz / height)


def place_water(cfg, host, frame, reference, lipids, definition, protected, origin):
    """Reference water on the new bilayer and the host solvent kept around it.

    Returns (reference waters [(n, 3, 3) positions], kept host solvent, removed host solvent by reason, report).
    """
    universe = host.universe
    box = universe.dimensions
    host_cell, height = _cell(box)
    heavy = _heavy_mask(definition)
    lipid_heavy = np.concatenate([lipid[heavy] for _, _, lipid in lipids]) if lipids else np.zeros((0, 3))
    protected_heavy = _heavy(universe.atoms[_indices(protected)]).positions
    solvent = [m for m in host.molecules if m.category in SOLVENT_CATEGORIES]
    solvent_heavy = [_heavy(universe.atoms[m.start:m.stop]).positions.astype(float) for m in solvent]
    removed = {"replaced_by_reference_water": [], "overlap_with_popc": [], "overlap_with_reference_water": []}

    reach = float(np.mean(np.abs([reference.report["leaflets"][s]["mean_z"] for s in ("upper", "lower")])))
    half = min(reach + cfg.water_layer, height / 2 - 3.0)
    water_definition = next((m.definition for m in solvent if m.category == "water"), None)
    use_water = len(reference.waters) and water_definition is not None and cfg.water_layer > 0
    report = {"slab_half_thickness_angstrom": round(half, 3), "reference_waters": 0}
    if use_water:
        types = [a.atom_type.upper() for a in water_definition.atoms]
        if len(types) != 3 or types[0].startswith("H") or not all(t.startswith("H") for t in types[1:]):
            use_water = False
            report["note"] = (f"Host water {water_definition.name} is not a three-site O H H water; "
                              "reference water was not used")
    if not use_water and "note" not in report:
        report["note"] = "No reference water (none in the reference, none in the host, or water_layer 0)"

    near_protein = _touching(solvent_heavy, protected_heavy, cfg.protein_hydration_distance, box)
    in_slab = set()
    if use_water:
        for i, positions in enumerate(solvent_heavy):
            if np.abs(_wrapped_dz(positions[:, 2], frame.midplane, height)).min() < half:
                in_slab.add(i)
        for i in sorted(in_slab - near_protein):
            removed["replaced_by_reference_water"].append(i)
    gone = set(removed["replaced_by_reference_water"])
    remaining = [i for i in range(len(solvent)) if i not in gone]
    for j in _touching([solvent_heavy[i] for i in remaining], lipid_heavy, cfg.water_clash_distance, box):
        removed["overlap_with_popc"].append(remaining[j])
    gone |= set(removed["overlap_with_popc"])

    waters = np.zeros((0, 3, 3))
    if use_water:
        placed = []
        for k, offset in tile(reference.waters[:, 0, :2], reference.cell, host_cell, origin):
            if abs(reference.waters[k, 0, 2]) < half:
                water = reference.waters[k].copy()
                water[:, :2] += offset
                water[:, 2] += frame.midplane
                placed.append(water)
        report["tiled_waters"] = len(placed)
        oxygens = [w[:1] for w in placed]
        drop = _touching(oxygens, protected_heavy, cfg.water_clash_distance, box)
        drop |= _touching(oxygens, lipid_heavy, cfg.reference_contact_distance, box)
        survivors = [i for i in range(len(placed)) if i not in drop]
        edge = _resolve_contacts([oxygens[i] for i in survivors], cfg.reference_contact_distance, box)
        survivors = [i for n, i in enumerate(survivors) if n not in edge]
        # Host solvent near the protein keeps its place; elsewhere the reference water does.
        kept_host = [i for i in range(len(solvent)) if i not in gone]
        host_positions, host_owner = _owner_map([solvent_heavy[i] for i in kept_host])
        ref_positions = np.concatenate([oxygens[i] for i in survivors]) if survivors else np.zeros((0, 3))
        losing_ref = set()
        if len(host_positions) and len(ref_positions):
            for a, b in capped_distance(host_positions, ref_positions, cfg.water_clash_distance, box=box,
                                        return_distances=False):
                host_index = kept_host[int(host_owner[a])]
                if host_index in in_slab:
                    losing_ref.add(int(b))
                elif host_index not in gone:
                    removed["overlap_with_reference_water"].append(host_index)
                    gone.add(host_index)
        survivors = [i for n, i in enumerate(survivors) if n not in losing_ref]
        waters = np.array([placed[i] for i in survivors]) if survivors else waters
        report.update({"reference_waters": len(waters), "removed_protein_overlap": len(drop),
                       "removed_box_edge_overlap": len(edge), "yielded_to_host_water": len(losing_ref),
                       "host_solvent_kept_near_protein_in_slab": len(in_slab & near_protein)})
    kept = [m for i, m in enumerate(solvent) if i not in gone]
    removed = {reason: [solvent[i] for i in sorted(set(items))] for reason, items in removed.items()}
    return waters, water_definition, kept, removed, report


def neutralize(cfg, host, kept_solvent, protected, lipids, definition):
    universe = host.universe
    net = sum(m.definition.charge for m in protected + kept_solvent) + len(lipids) * definition.charge
    report = {"before": round(net, 6), "target": cfg.target_net_charge, "enabled": cfg.neutralize,
              "removed_ions": []}
    if not cfg.neutralize or abs(net - cfg.target_net_charge) <= 0.01:
        report["after"] = round(net, 6)
        return kept_solvent, report
    excess = net - cfg.target_net_charge
    groups = []
    if protected:
        groups.append((_heavy(universe.atoms[_indices(protected)]).positions, cfg.ion_exclusion_distance))
    if lipids:
        heavy = _heavy_mask(definition)
        groups.append((np.concatenate([lipid[heavy] for _, _, lipid in lipids]), cfg.lipid_ion_exclusion_distance))
    candidates = []
    for i, mol in enumerate(kept_solvent):
        if mol.category != "ion" or mol.definition.charge * excess <= 0:
            continue
        atoms = universe.atoms[mol.start:mol.stop]
        if all(not len(capped_distance(atoms.positions, positions, distance, box=universe.dimensions,
                                       return_distances=False)) for positions, distance in groups):
            candidates.append((i, mol.name, mol.definition.charge))
    removed, residual = select_counterions(candidates, excess, 0.01, 42)
    if abs(residual) > 0.01:
        raise ConfigError(f"Cannot reach net charge {cfg.target_net_charge:g} by removing counterions far from "
                          f"the protein and lipids: net charge {net:+.3f}, residual {residual:+.3f}")
    ids = {item["resindex"] for item in removed}
    report["removed_ions"] = [{"moleculetype": kept_solvent[i].name, "occurrence": kept_solvent[i].occurrence,
                               "charge": kept_solvent[i].definition.charge} for i in sorted(ids)]
    report["after"] = round(net - sum(item["charge"] for item in removed), 6)
    return [m for i, m in enumerate(kept_solvent) if i not in ids], report


ION_MIN_SPACING = 5.0    # A between an added ion and any other ion


def restore_salt(cfg, host, kept_solvent, n_reference_waters, frame, half, protected):
    """Turn bulk host waters into salt pairs until the host concentration is restored.

    Returns (kept solvent, [(ion definition, position)], report).
    """
    universe = host.universe
    solvent = [m for m in host.molecules if m.category in SOLVENT_CATEGORIES]
    ions = [m for m in solvent if m.category == "ion"]
    counts = Counter(m.name for m in ions)
    cations = [n for n in counts if next(m for m in ions if m.name == n).definition.charge > 0.99
               and next(m for m in ions if m.name == n).definition.charge < 1.01]
    anions = [n for n in counts if -1.01 < next(m for m in ions if m.name == n).definition.charge < -0.99]
    report = {"enabled": cfg.restore_salt, "added_pairs": 0}
    if not cfg.restore_salt or not cations or not anions:
        if cfg.restore_salt:
            report["note"] = "The host has no monovalent salt to restore"
        return kept_solvent, [], report
    cation = max(cations, key=lambda n: counts[n])
    anion = max(anions, key=lambda n: counts[n])
    definitions = {m.name: m.definition for m in ions}
    host_waters = sum(m.category == "water" for m in solvent)
    target = min(counts[cation], counts[anion]) / host_waters * WATER_MOLARITY if host_waters else 0.0
    kept_counts = Counter(m.name for m in kept_solvent)
    pairs = min(kept_counts[cation], kept_counts[anion])
    waters = sum(m.category == "water" for m in kept_solvent) + n_reference_waters
    # pairs + n = target (waters - 2 n) / 55.51: each added pair replaces two waters.
    needed = int(round((target * waters / WATER_MOLARITY - pairs) / (1 + 2 * target / WATER_MOLARITY)))
    report.update({"target_concentration_molar": round(target, 4), "cation": cation, "anion": anion})
    if needed <= 0:
        return kept_solvent, [], report
    box = universe.dimensions
    _, height = _cell(box)
    protected_heavy = _heavy(universe.atoms[_indices(protected)]).positions
    candidates = []
    for i, mol in enumerate(kept_solvent):
        if mol.category != "water":
            continue
        oxygen = _heavy(universe.atoms[mol.start:mol.stop]).positions[0].astype(float)
        if abs(_wrapped_dz(oxygen[2], frame.midplane, height)) >= half:
            candidates.append((i, oxygen))
    if candidates and len(protected_heavy):
        far = capped_distance(np.array([o for _, o in candidates]), protected_heavy, cfg.ion_exclusion_distance,
                              box=box, return_distances=False)
        close = {int(a) for a, _ in far}
        candidates = [c for n, c in enumerate(candidates) if n not in close]
    occupied = [universe.atoms[m.start].position.astype(float) for m in kept_solvent if m.category == "ion"]
    rng = np.random.default_rng(42)
    order = rng.permutation(len(candidates))
    picked = []   # (kept_solvent index, ion definition, position), cation and anion alternating
    for n in order:
        if len(picked) == 2 * needed:
            break
        i, oxygen = candidates[n]
        if occupied and len(capped_distance(oxygen[None, :], np.array(occupied), ION_MIN_SPACING, box=box,
                                            return_distances=False)):
            continue
        picked.append((i, definitions[cation if len(picked) % 2 == 0 else anion], oxygen))
        occupied.append(oxygen)
    picked = picked[:len(picked) // 2 * 2]   # whole pairs only
    added = [(definition, position) for _, definition, position in picked]
    chosen = {i for i, _, _ in picked}
    report["added_pairs"] = len(added) // 2
    if report["added_pairs"] < needed:
        report["note"] = f"Only {report['added_pairs']} of {needed} salt pairs found room in the bulk water"
    return [m for i, m in enumerate(kept_solvent) if i not in chosen], added, report


def _placed_universe(definition, blocks, workdir):
    """New molecules of one moleculetype as a GRO-read universe, one residue each."""
    n_atoms, count = definition.atom_count, len(blocks)
    universe = mda.Universe.empty(n_atoms * count, n_residues=count,
                                  atom_resindex=np.repeat(np.arange(count), n_atoms), trajectory=True)
    universe.add_TopologyAttr("names", [a.atom_name for a in definition.atoms] * count)
    universe.add_TopologyAttr("resnames", [definition.atoms[0].resname] * count)
    universe.add_TopologyAttr("resids", np.arange(1, count + 1))
    universe.atoms.positions = np.concatenate(blocks)
    universe.dimensions = [100.0, 100.0, 100.0, 90.0, 90.0, 90.0]
    path = Path(workdir) / f"placed_{definition.name}.gro"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        universe.atoms.write(str(path))
        reread = mda.Universe(str(path))
    path.unlink()
    return reread


class _Placed:
    """The placed POPC block: a universe and its molecules, like a BuildSystem."""

    def __init__(self, universe, molecules):
        self.universe, self.molecules = universe, molecules


def _salt(molecules):
    waters = sum(m.category == "water" for m in molecules)
    ions = Counter(m.name for m in molecules if m.category == "ion")
    cations = sum(n for name, n in ions.items()
                  if next(m.definition.charge for m in molecules if m.name == name) > 0)
    anions = sum(ions.values()) - cations
    pairs = min(cations, anions)
    return {"waters": waters, "ions": dict(ions),
            "concentration_molar": round(pairs / waters * WATER_MOLARITY, 4) if waters else None}


def _write_index(path, ordered, leaflet_of):
    groups = {name: [] for name in ("System", "Protein", "Ligand", "Water", "Ion", "Membrane",
                                    "POPC_upper", "POPC_lower", "Environment",
                                    "SOLU", "MEMB", "SOLV", "SOLU_MEMB", "SYSTEM")}
    # CHARMM-GUI groups, so CHARMM-GUI mdp files run unchanged.
    charmm_gui = {"protein": ("SOLU", "SOLU_MEMB"), "ligand": ("SOLU", "SOLU_MEMB"),
                  "lipid": ("MEMB", "SOLU_MEMB"), "water": ("SOLV",), "ion": ("SOLV",), "solvent": ("SOLV",)}
    labels = {"protein": "Protein", "ligand": "Ligand", "water": "Water", "ion": "Ion", "lipid": "Membrane"}
    offset = 0
    for _, mol in ordered:
        ids = list(range(offset + 1, offset + 1 + mol.stop - mol.start))
        groups["System"].extend(ids)
        groups["SYSTEM"].extend(ids)
        for name in charmm_gui.get(mol.category, ("SOLU", "SOLU_MEMB")):
            groups[name].extend(ids)
        if mol.category in labels:
            groups[labels[mol.category]].extend(ids)
        if id(mol) in leaflet_of:
            groups[f"POPC_{leaflet_of[id(mol)]}"].extend(ids)
        if mol.category not in ("protein", "ligand"):
            groups["Environment"].extend(ids)
        offset += len(ids)
    Path(path).write_text("\n".join(
        f"[ {name} ]\n" + "\n".join(" ".join(map(str, ids[i:i + 15])) for i in range(0, len(ids), 15))
        for name, ids in groups.items() if ids) + "\n")
    return groups


# ---------------------------------------------------------------------------
# Equilibration
# ---------------------------------------------------------------------------

def write_equilibration(spec, stage_dir, topology, solute_definitions, index_groups):
    """A self-contained GROMACS equilibration folder; returns its report."""
    from repack import _mdp, restraint_coverage
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
    groups = ("SOLU_MEMB SOLV", "SOLU_MEMB SOLV") if {"SOLU_MEMB", "SOLV"} <= set(index_groups) else ("System", "System")
    first_md = next((i for i, s in enumerate(spec.stages) if s.integrator != "steep"), None)
    for number, stage in enumerate(spec.stages):
        define = (f"-DPOSRES -DPOSRES_FC_BB={stage.backbone:.1f} -DPOSRES_FC_SC={stage.sidechain:.1f} "
                  f"-DPOSRES_FC_LIPID={stage.lipid:.1f} -DDIHRES -DDIHRES_FC={stage.lipid:.1f}")
        (bundle / f"{stage.name}.mdp").write_text(_mdp(stage, spec, number == first_md, groups, define=define))
    md_ps = sum(s.nsteps * s.dt for s in spec.stages if s.integrator != "steep")
    script = ["#!/bin/bash",
              "# Sample run script written by GMXtransplant. Adapt it to your machine, or submit",
              "# the same grompp/mdrun lines through your scheduler.",
              f"# Equilibration of the new POPC bilayer, {md_ps / 1000:g} ns of MD. The protein is held and",
              "# released step by step (CHARMM-GUI schedule), so it keeps its OPM placement while the",
              "# POPC packs around it. Restraints (kJ/mol/nm^2): backbone / side chains / lipid P z + dihedrals",]
    script += [f"#   {s.name}: {s.backbone:g} / {s.sidechain:g} / {s.lipid:g}" for s in spec.stages]
    script += ["# Afterwards, continue production from the last GRO with the same topol.top.",
               "set -e", 'cd "$(dirname "$0")"', 'GMX=${GMX:-gmx}', 'MDRUN=${MDRUN:-"$GMX mdrun"}',
               "previous=step5_input.gro"]
    for stage in spec.stages:
        script += [f'echo "==> {stage.name}"',
                   f"$GMX grompp -f {stage.name}.mdp -o {stage.name}.tpr -c $previous -r step5_input.gro "
                   "-p topol.top -n index.ndx -maxwarn 0",
                   f"$MDRUN -v -deffnm {stage.name}", f"previous={stage.name}.gro"]
    (bundle / "sample_run.sh").write_text("\n".join(script) + "\n")
    os.chmod(bundle / "sample_run.sh", 0o755)
    return {"directory": spec.output_dir, "md_length_ps": md_ps, "temperature": spec.temperature,
            "coupling_groups": list(groups), "stages": [asdict(s) for s in spec.stages],
            "restraint_coverage": coverage, "warnings": warnings_}


# ---------------------------------------------------------------------------
# Report and run
# ---------------------------------------------------------------------------

def _write_report(folder, report):
    (folder / "popc_report.json").write_text(json.dumps(report, indent=2) + "\n")
    lines = ["GMXtransplant: replace the membrane with equilibrated POPC", f"Status: {report['status']}"]
    if "error" in report:
        lines += ["", f"Error: {report['error']}"]
    host = report.get("host_membrane")
    if host:
        lines += ["", "HOST MEMBRANE (protein and ligands are not moved: OPM placement kept)",
                  f"  midplane z {host['midplane_z']:.2f} A ({host['midplane_source']}); "
                  f"P-P thickness {host['phosphate_thickness_angstrom']:.2f} A"]
        for side in ("upper", "lower"):
            plane = host["planes"][side]
            lines.append(f"  {side}: P plane z {plane['mean_z']:.2f} A, tilt {plane['tilt_degrees']:.2f} deg, "
                         f"spread {plane['spread_angstrom']:.2f} A; lipids {host['lipids_removed'][side]}")
            anchors = host["anchors"][side]
            lines.append(f"    {len(anchors)} POPC anchors at least "
                         f"{host['anchors'][side + '_min_separation_angstrom']:.1f} A apart: "
                         + ", ".join(f"POPC {a['resid']}" for a in anchors))
    ref = report.get("reference")
    if ref:
        lines += ["", "REFERENCE BILAYER", f"  {ref['path']}", f"  {ref.get('source', 'user supplied')}",
                  f"  {ref['lipids']} POPC, box {ref['box_angstrom'][0]:.2f} x {ref['box_angstrom'][1]:.2f} A, "
                  f"P-P thickness {ref['phosphate_thickness_angstrom']:.2f} A, area per lipid "
                  f"{ref['area_per_lipid_angstrom2']['upper']:.1f} / {ref['area_per_lipid_angstrom2']['lower']:.1f} A^2"]
    placement = report.get("placement")
    if placement:
        lines += ["", "POPC PLACEMENT (upper / lower)"]
        for key, label in (("tiled", "cut to the host box"), ("removed_protein_overlap", "removed, protein/ligand overlap"),
                           ("removed_box_edge_overlap", "removed, overlap at the box edge"), ("placed", "placed")):
            values = placement[key]
            lines.append(f"  {label:<36} {values.get('upper', 0):>5} / {values.get('lower', 0):<5}")
    water = report.get("water")
    if water:
        lines += ["", "WATER"]
        if water.get("reference_waters"):
            lines.append(f"  {water['reference_waters']} reference waters within "
                         f"{water['slab_half_thickness_angstrom']:.1f} A of the midplane hydrate the new bilayer; "
                         f"{water['host_solvent_kept_near_protein_in_slab']} host water/ion molecules near the "
                         "protein kept there")
        if water.get("note"):
            lines.append(f"  {water['note']}")
    for reason, counts in report.get("solvent_removed", {}).items():
        lines.append(f"  host solvent removed, {reason.replace('_', ' ')}: "
                     + (", ".join(f"{k} {v}" for k, v in counts.items()) or "none"))
    if "charge" in report:
        charge = report["charge"]
        lines += ["", f"CHARGE: {charge['before']:+.3f} e -> {charge['after']:+.3f} e; "
                      f"{len(charge['removed_ions'])} counterion(s) removed"]
    if "salt" in report:
        lines += [f"SALT: host {report['salt']['host']['concentration_molar']} M, "
                  f"final {report['salt']['final']['concentration_molar']} M"]
    if "final_molecules" in report:
        lines += ["", "FINAL MOLECULES"] + [f"  {name:<10} {count}" for name, count in report["final_molecules"]]
    eq = report.get("equilibration")
    if eq:
        lines += ["", f"EQUILIBRATION: {eq['directory']}/ ({eq['md_length_ps'] / 1000:g} ns of MD, "
                      f"{len(eq['stages'])} stages); run sample_run.sh"]
        lines += [f"  WARNING: {w}" for w in eq["warnings"]]
    if report.get("warnings"):
        lines += ["", "WARNINGS"] + [f"  {w}" for w in report["warnings"]]
    (folder / "popc_report.txt").write_text("\n".join(lines) + "\n")


def _write_scene(stage, final, ordered):
    from visualization import safe_scene
    roles, offset = {}, 0
    for system, mol in ordered:
        if isinstance(system, _Placed):
            role = "inserted_popc"
        else:
            role = "retained_" + {"solvent": "other", "detergent": "lipid"}.get(mol.category, mol.category)
        roles.setdefault(role, []).extend(range(offset, offset + mol.stop - mol.start))
        offset += mol.stop - mol.start
    groups = {role: final.atoms[ids] for role, ids in roles.items()}
    safe_scene(stage, groups, "popc", palette={"inserted_popc": ("Equilibrated POPC", [.95, .68, .12], True)},
               alignment="Protein and ligands keep their input coordinates; the POPC bilayer is new.")


def run_popc(config_path, dry_run=False, output_root=None):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        return _run_popc(config_path, dry_run, output_root)


def _run_popc(config_path, dry_run=False, output_root=None):
    cfg = load_popc_config(config_path, check_paths=not dry_run, output_root=output_root)
    reference_path = cfg.reference_gro or str(default_reference_path())
    if dry_run:
        print(f"popc configuration validated; no inputs opened. Reference: {reference_path}; "
              f"output: {cfg.output_dir}")
        return 0
    destination = Path(cfg.output_dir)
    if destination.exists() and not destination.is_dir():
        raise ConfigError(f"output_dir exists and is not a directory: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    print("=" * 78)
    print("REPLACE THE MEMBRANE WITH EQUILIBRATED POPC")
    print("=" * 78)
    report = {"mode": "popc", "status": "running", "generated_at": datetime.now(timezone.utc).isoformat(),
              "host": cfg.host, "settings": {k: v for k, v in asdict(cfg).items() if k != "equilibration"},
              "warnings": []}
    with tempfile.TemporaryDirectory(prefix=".popc-", dir=destination.parent) as temporary:
        stage = Path(temporary) / "result"
        stage.mkdir()
        try:
            host = inspect_system(cfg.host, gro=cfg.host_gro, overrides=cfg.environment_overrides)
            print(f"[host] {host.gro}: proteins={host.report['selected_proteins']}, "
                  f"ligands={host.report['selected_ligands']}")
            frame = measure_host(host, cfg)
            removed_lipids = [m for m in host.molecules if m.category in REMOVED_CATEGORIES]
            report["host_membrane"] = {
                "midplane_z": round(frame.midplane, 3),
                "midplane_source": "membrane_center_z" if cfg.membrane_center_z is not None
                                   else "midway between the P planes",
                "phosphate_thickness_angstrom": round(frame.planes["upper"]["mean_z"]
                                                      - frame.planes["lower"]["mean_z"], 3),
                "planes": frame.planes, "anchors": frame.anchors,
                "lipids_removed": {side: dict(frame.leaflet_counts[side]) for side in ("upper", "lower")},
                "molecules_removed": dict(Counter(m.name for m in removed_lipids))}
            for side in ("upper", "lower"):
                plane = frame.planes[side]
                print(f"[host] {side} leaflet: P plane z {plane['mean_z']:.1f} A, tilt {plane['tilt_degrees']:.1f} deg; "
                      f"{cfg.anchor_count} POPC anchors >= {frame.anchors[side + '_min_separation_angstrom']:.1f} A apart")
            print(f"[host] midplane z {frame.midplane:.2f} A; removing "
                  + ", ".join(f"{n} {name}" for name, n in report["host_membrane"]["molecules_removed"].items()))
            definition = host.definitions[LIPID]
            reference = load_reference(reference_path, definition)
            report["reference"] = dict(reference.report)
            if not cfg.reference_gro:
                report["reference"]["source"] = REFERENCE_SOURCE
            print(f"[reference] {reference.report['lipids']} POPC, P-P "
                  f"{reference.report['phosphate_thickness_angstrom']:.1f} A (host "
                  f"{report['host_membrane']['phosphate_thickness_angstrom']:.1f} A)")
            difference = (reference.report["phosphate_thickness_angstrom"]
                          - report["host_membrane"]["phosphate_thickness_angstrom"])
            if abs(difference) > 3.0:
                report["warnings"].append(
                    f"The reference bilayer is {difference:+.1f} A thicker than the host's; the water layers "
                    "adjust during the NPT equilibration")

            protected = [m for m in host.molecules if m.category not in REMOVED_CATEGORIES + SOLVENT_CATEGORIES]
            lipids, origin, report["placement"] = place_popc(cfg, host, frame, reference, definition, protected)
            search = report["placement"]["lateral_search"]
            print(f"[popc] {report['placement']['placed'].get('upper', 0)} upper / "
                  f"{report['placement']['placed'].get('lower', 0)} lower placed (best of "
                  f"{search['registrations_tried']} lateral registrations); removed for protein overlap "
                  f"{sum(report['placement']['removed_protein_overlap'].values())}, at the box edge "
                  f"{sum(report['placement']['removed_box_edge_overlap'].values())}")
            waters, water_definition, kept_solvent, removed, report["water"] = place_water(
                cfg, host, frame, reference, lipids, definition, protected, origin)
            report["solvent_removed"] = {reason: dict(Counter(m.name for m in items))
                                         for reason, items in removed.items()}
            if report["water"]["reference_waters"]:
                print(f"[water] {report['water']['reference_waters']} reference waters hydrate the bilayer "
                      f"(within {report['water']['slab_half_thickness_angstrom']:.1f} A of the midplane); host "
                      f"solvent near the protein kept: {report['water']['host_solvent_kept_near_protein_in_slab']}")
            elif "note" in report["water"]:
                report["warnings"].append(report["water"]["note"])
            for reason, counts in report["solvent_removed"].items():
                if counts:
                    print(f"[water] host solvent removed, {reason.replace('_', ' ')}: "
                          + ", ".join(f"{n} {name}" for name, n in counts.items()))
            kept_solvent, report["charge"] = neutralize(cfg, host, kept_solvent, protected, lipids, definition)
            print(f"[charge] {report['charge']['before']:+.3f} e -> {report['charge']['after']:+.3f} e; "
                  f"{len(report['charge']['removed_ions'])} counterion(s) removed")
            kept_solvent, added_ions, salt_report = restore_salt(
                cfg, host, kept_solvent, len(waters), frame, report["water"]["slab_half_thickness_angstrom"],
                protected)
            final_ions = Counter(m.name for m in kept_solvent if m.category == "ion")
            final_ions.update(d.name for d, _ in added_ions)
            final_waters = sum(m.category == "water" for m in kept_solvent) + len(waters)
            pairs = min(final_ions[salt_report["cation"]], final_ions[salt_report["anion"]]) \
                if "cation" in salt_report else 0
            report["salt"] = {"host": _salt([m for m in host.molecules if m.category in SOLVENT_CATEGORIES]),
                              "final": {"waters": final_waters, "ions": dict(final_ions),
                                        "concentration_molar": round(pairs / final_waters * WATER_MOLARITY, 4)
                                        if final_waters else None},
                              "restore": salt_report}
            print(f"[salt] host {report['salt']['host']['concentration_molar']} M -> final "
                  f"{report['salt']['final']['concentration_molar']} M"
                  + (f"; {salt_report['added_pairs']} pair(s) added in the bulk water"
                     if salt_report["added_pairs"] else ""))
            if salt_report.get("note"):
                report["warnings"].append(salt_report["note"])

            size = definition.atom_count
            placed_molecules = [Molecule(LIPID, n + 1, n * size, (n + 1) * size, "lipid", True, definition)
                                for n in range(len(lipids))]
            placed = _Placed(_placed_universe(definition, [lipid for _, _, lipid in lipids], stage),
                             placed_molecules)
            leaflet_of = {id(mol): side for mol, (_, side, _) in zip(placed_molecules, lipids)}
            water_block = None
            if len(waters):
                water_molecules = [Molecule(water_definition.name, n + 1, 3 * n, 3 * n + 3, "water", True,
                                            water_definition) for n in range(len(waters))]
                water_block = _Placed(_placed_universe(water_definition, list(waters), stage), water_molecules)
            # New molecules join the end of their type's block in the host order: added
            # ions after the host's last kept ion of that type, reference water after its water.
            appended = {}
            if water_block is not None:
                appended[water_definition.name] = [water_block]
            for ion_definition in {d.name: d for d, _ in added_ions}.values():
                positions = [p[None, :] for d, p in added_ions if d.name == ion_definition.name]
                block = _Placed(_placed_universe(ion_definition, positions, stage),
                                [Molecule(ion_definition.name, n + 1, n, n + 1, "ion", True, ion_definition)
                                 for n in range(len(positions))])
                appended.setdefault(ion_definition.name, []).append(block)
            kept_ids = {id(m) for m in protected + kept_solvent}
            last = {}
            for i, mol in enumerate(host.molecules):
                if id(mol) in kept_ids:
                    last[mol.name] = i
            ordered, inserted = [], False
            for i, mol in enumerate(host.molecules):
                if mol.category in REMOVED_CATEGORIES:
                    if not inserted:
                        ordered += [(placed, m) for m in placed_molecules]
                        inserted = True
                elif id(mol) in kept_ids:
                    ordered.append((host, mol))
                if last.get(mol.name) == i:
                    for block in appended.pop(mol.name, []):
                        ordered += [(block, m) for m in block.molecules]
            for blocks in appended.values():
                for block in blocks:
                    ordered += [(block, m) for m in block.molecules]
            groups, current, blocks = [], None, []
            for system, mol in ordered:
                if current is not None and system is not current:
                    groups.append(current.universe.atoms[_indices(blocks)])
                    blocks = []
                current = system
                blocks.append(mol)
            groups.append(current.universe.atoms[_indices(blocks)])
            final = mda.Merge(*groups)
            final.dimensions = host.universe.dimensions.copy()

            print("[topology] Writing and auditing topol.top, toppar/ and the coordinates")
            topology = write_topology(stage, host, host, ordered)
            top = Path(topology.top_path)
            top.write_text(top.read_text().replace("CHARMM-GUI protein transplant", "GMXtransplant popc"))
            audit = audit_final_topology(final, topology)
            if abs(audit["net_charge"] - report["charge"]["after"]) > 0.01:
                raise ConfigError("Final topology charge disagrees with molecule charge accounting")
            protein_counts = [m.stop - m.start for _, m in ordered if m.category == "protein"]

            def validate_written(_pdb, gro):
                return audit_final_topology(mda.Universe(gro), topology)

            write_outputs(final, str(stage / "step5_input.pdb"), str(stage / "step5_input.gro"),
                          protein_chain_atom_counts=protein_counts, staged_validator=validate_written)
            index_groups = _write_index(stage / "index.ndx", ordered, leaflet_of)
            report["final_molecules"] = topology.molecules_written
            report["final_net_charge"] = round(audit["net_charge"], 6)
            if cfg.equilibration.enabled:
                solute = {mol.name: mol.definition for mol in protected}
                report["equilibration"] = write_equilibration(cfg.equilibration, stage, topology,
                                                              list(solute.values()), index_groups)
                print(f"[equilibration] {cfg.equilibration.output_dir}/: "
                      f"{report['equilibration']['md_length_ps'] / 1000:g} ns in {len(cfg.equilibration.stages)} "
                      "stages; run sample_run.sh")
                for warning in report["equilibration"]["warnings"]:
                    print(f"      WARNING: {warning}")
            _write_scene(stage, final, ordered)
            if (stage / "view.pml").is_file():
                pml = stage / "view.pml"
                pml.write_text(pml.read_text().replace(repr(str(stage)), repr(str(destination))))
            report["status"] = "complete"
            _write_report(stage, report)
            _publish_directory(stage, destination, GENERATED)
        except (ConfigError, TopologyError) as exc:
            report["status"], report["error"] = "failed", str(exc)
            failure = Path(temporary) / "failure"
            failure.mkdir()
            _write_report(failure, report)
            _publish_directory(failure, destination, GENERATED)
            raise ConfigError(str(exc)) from exc
    for warning in report["warnings"]:
        print(f"WARNING: {warning}")
    print(f"[complete] {destination}")
    from run_pipeline import _announce_output_location
    _announce_output_location(destination)
    return 0
