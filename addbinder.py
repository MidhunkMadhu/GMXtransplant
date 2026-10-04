"""Place an incoming binder (ligand or protein) in the water above or below a membrane protein.

The host is one prepared CHARMM-GUI-style folder (topol.top, toppar/, GRO) that
is kept unchanged except for water and ions. The binder comes as coordinates
plus its ITP (and, for CHARMM-GUI ligands, the forcefield.itp generated with
it). A binder made of several molecules (e.g. a G-protein heterotrimer) gives
one ITP per molecule; the coordinates hold them in that order and the set is
placed as one rigid body. Each requested pose places the binder ``distance`` Angstrom beyond the
host's extreme heavy atom (the tip) on the chosen membrane side, checks that it
clears the host, the membrane and every periodic image, and writes a complete
GRO/PDB + topol.top + toppar/ + index.ndx system per accepted pose.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace
from datetime import datetime, timezone
from itertools import product
from pathlib import Path
from typing import Dict, List, Optional, Union
import json
import re
import zlib
import shutil
import tempfile
import warnings

import MDAnalysis as mda
from MDAnalysis.lib.distances import capped_distance
import numpy as np
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation
import yaml

from charmprot import (inspect_system, Molecule, write_topology, _included_files, _indices,
                       _heavy)
from classify import AMINO_ACID_RESNAMES
from config import ConfigError, _UniqueKeySafeLoader, anchor_config_paths
from ion_selection import select_counterions
from itp import parse_itp
from output import write_outputs
from topology import audit_final_topology, check_ligand_parameters, TopologyError


ORIENTATIONS = ("auto", "flat", "end_on", "edge", "as_is", "euler")
_POSE_KEYS = {"name", "distance", "lateral_offset", "orientation", "angles", "spin", "flip", "approach", "centroid", "from_input_position",
              "distance_to", "reduce_distance_by", "min_distance"}
RANDOM_ATTEMPTS = 2000
CENTROID_MIN_GAP = 3.0  # A; a centroid pose closer than this to the host overlaps it
# auto poses take these (orientation, flip) in turn, skipping any an explicit
# pose already uses; past the end the cycle repeats turned 45 degrees further.
AUTO_ORIENTATIONS = [("flat", False), ("end_on", False), ("end_on", True),
                     ("edge", False), ("flat", True), ("edge", True)]


@dataclass
class PoseSpec:
    name: str = ""                    # empty: binderpose1, binderpose2, ...
    distance: Optional[float] = None  # None: the top-level distance
    lateral_offset: List[float] = field(default_factory=lambda: [0.0, 0.0])
    # flat: thinnest binder axis along the membrane normal, longest along x.
    # end_on: longest axis along the normal, its far-reaching end towards the
    # host. edge: middle axis along the normal. as_is: keep the input
    # orientation. euler: ZYZ angles (degrees) about the binder centre. auto:
    # resolved at load time to a different orientation for each pose.
    # flip turns the binder upside down (180 degrees about x); spin (degrees,
    # about the normal) is applied last.
    orientation: str = "auto"
    angles: List[float] = field(default_factory=lambda: [0.0, 0.0, 0.0])
    spin: float = 0.0
    flip: bool = False
    # [tilt, azimuth] in degrees: come in along this direction instead of straight
    # down the membrane normal, stopping at distance from the nearest host heavy atom.
    approach: Optional[List[float]] = None
    # [x, y, z] in Angstrom in the host GRO's frame (GRO nm x 10): put the binder's
    # heavy-atom centre exactly here; distance is then only reported.
    centroid: Optional[List[float]] = None
    # Start from where the binder is in its input file (e.g. its bound pose, already
    # in the host's frame) and move it straight out along the membrane normal (or
    # along approach) until its nearest heavy atom is distance from the host.
    from_input_position: bool = False
    # nearest_atom: distance is the gap between the nearest heavy atoms. centroid:
    # distance is from the tip (or host centre, see lateral_reference) to the
    # binder's heavy-atom centre; with approach [tilt, azimuth] this is the
    # binder's position in spherical coordinates (r, theta, phi).
    distance_to: str = "nearest_atom"
    # When the pose fails the box checks, retry at distance - reduce_distance_by,
    # and so on, down to min_distance. 0 keeps the one distance.
    reduce_distance_by: float = 0.0
    min_distance: float = 10.0


@dataclass
class AddBinderSpec:
    host: str
    binder_coordinates: str
    # One ITP path, or a list of ITPs (one molecule each) for a multi-molecule binder.
    binder_itp: Union[str, List[str]]
    binder_forcefield: str | None = None
    host_gro: str | None = None
    side: str = "upper"
    distance: float = 10.0
    min_image_gap: float = 10.0
    min_membrane_gap: float = 5.0
    lateral_reference: str = "tip"
    water_clash_distance: float = 2.4
    neutralize: bool = True
    target_net_charge: float = 0.0
    same_orientation: bool = False   # every auto pose gets the first auto orientation
    copy_run_inputs: bool = True     # copy the host's *.mdp and README into each pose folder
    # Extra poses from random directions, distances and orientations, numbered on after the listed ones
    random_poses: int = 3
    random_distance: Optional[List[float]] = None   # [min, max] A; None: [distance, distance + 5]
    random_max_tilt: float = 60.0                   # degrees from the membrane normal
    random_min_angle: float = 25.0                  # degrees between any two approach directions
    random_min_separation: float = 10.0             # A between binder centres of any two poses
    random_seed: int = 1
    # How many random poses also get a random rotation (the last ones); the
    # others keep the input orientation. None: all of them.
    random_rotated: Optional[int] = None
    random_from_input_position: bool = False  # true: move out from the input position (e.g. a bound pose)
    ion_exclusion_distance: float = 10.0
    lipid_ion_exclusion_distance: float = 5.0
    ion_side: str = "any"            # any: the whole box; binder_side: the binder's side first
    # Salt after the binder displaces water and ions: "host" keeps the host's
    # concentration (pairs x 55.51 / waters), or a number in mol/L.
    salt_concentration: Union[str, float] = "host"
    output_dir: str = "addbinder_output"
    environment_overrides: Dict[str, str] = field(default_factory=dict)
    poses: List[PoseSpec] = field(default_factory=list)


def _number(value, key, positive=True):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not np.isfinite(value):
        raise ConfigError(f"{key} must be a finite number")
    if positive and value <= 0:
        raise ConfigError(f"{key} must be positive (Angstrom)")
    return float(value)


def _parse_pose(raw, index, default_distance):
    if not isinstance(raw, dict):
        raise ConfigError(f"poses[{index}] must be a mapping")
    unknown = set(raw) - _POSE_KEYS
    if unknown:
        raise ConfigError(f"Unknown options in poses[{index}]: {sorted(unknown)}")
    name = raw.get("name") or f"binderpose{index + 1}"
    if not isinstance(name, str) or not re.fullmatch(r"[A-Za-z0-9_.-]+", name):
        raise ConfigError(f"poses[{index}].name must be a simple folder name")
    offset = raw.get("lateral_offset", [0.0, 0.0])
    if not isinstance(offset, (list, tuple)) or len(offset) != 2:
        raise ConfigError(f"poses[{index}].lateral_offset must be [dx, dy] in Angstrom")
    angles = raw.get("angles", [0.0, 0.0, 0.0])
    if not isinstance(angles, (list, tuple)) or len(angles) != 3:
        raise ConfigError(f"poses[{index}].angles must be three ZYZ angles in degrees")
    orientation = raw.get("orientation", "auto")
    if orientation not in ORIENTATIONS:
        raise ConfigError(f"poses[{index}].orientation must be one of {ORIENTATIONS}")
    if orientation != "euler" and any(v != 0 for v in angles):
        raise ConfigError(f"poses[{index}].angles requires orientation: euler")
    flip = raw.get("flip", False)
    if type(flip) is not bool:
        raise ConfigError(f"poses[{index}].flip must be true or false")
    if flip and orientation == "auto":
        raise ConfigError(f"poses[{index}].flip needs an explicit orientation; auto picks its own")
    approach = raw.get("approach")
    if approach is not None:
        if not isinstance(approach, (list, tuple)) or len(approach) != 2:
            raise ConfigError(f"poses[{index}].approach must be [tilt, azimuth] in degrees")
        approach = [_number(v, f"poses[{index}].approach", False) for v in approach]
        if not 0 <= approach[0] < 90:
            raise ConfigError(f"poses[{index}].approach tilt must be 0-90 degrees from the membrane normal")
    from_input = raw.get("from_input_position", False)
    if type(from_input) is not bool:
        raise ConfigError(f"poses[{index}].from_input_position must be true or false")
    if from_input and (raw.get("centroid") is not None or any(v != 0 for v in offset)):
        raise ConfigError(f"poses[{index}].from_input_position starts from the input position; leave out "
                          "centroid and lateral_offset")
    if from_input and orientation == "auto":
        orientation = "as_is"  # the input pose is the point of starting from the input position
    distance_to = raw.get("distance_to", "nearest_atom")
    if distance_to not in ("nearest_atom", "centroid"):
        raise ConfigError(f"poses[{index}].distance_to must be nearest_atom or centroid")
    if distance_to == "centroid" and from_input:
        raise ConfigError(f"poses[{index}].from_input_position measures distance to the nearest atom; "
                          "leave out distance_to: centroid")
    reduce_by = _number(raw.get("reduce_distance_by", 0.0), f"poses[{index}].reduce_distance_by", False)
    min_distance = _number(raw.get("min_distance", 10.0), f"poses[{index}].min_distance")
    if reduce_by < 0:
        raise ConfigError(f"poses[{index}].reduce_distance_by must not be negative")
    centroid = raw.get("centroid")
    if centroid is not None and (distance_to == "centroid" or reduce_by):
        raise ConfigError(f"poses[{index}].centroid fixes the position; leave out distance_to and "
                          "reduce_distance_by")
    if centroid is not None:
        if not isinstance(centroid, (list, tuple)) or len(centroid) != 3:
            raise ConfigError(f"poses[{index}].centroid must be [x, y, z] in Angstrom")
        centroid = [_number(v, f"poses[{index}].centroid", False) for v in centroid]
        if approach is not None or any(v != 0 for v in offset) or raw.get("distance") is not None:
            raise ConfigError(f"poses[{index}].centroid fixes the position; leave out approach, "
                              "lateral_offset and distance")
    return PoseSpec(
        name=name,
        distance=_number(default_distance if raw.get("distance") is None else raw["distance"],
                         f"poses[{index}].distance"),
        lateral_offset=[_number(v, f"poses[{index}].lateral_offset", False) for v in offset],
        orientation=orientation,
        angles=[_number(v, f"poses[{index}].angles", False) for v in angles],
        spin=_number(raw.get("spin", 0.0), f"poses[{index}].spin", False),
        flip=flip,
        approach=approach,
        centroid=centroid,
        from_input_position=from_input,
        distance_to=distance_to,
        reduce_distance_by=reduce_by,
        min_distance=min_distance,
    )


def _resolve_auto_orientations(poses, same):
    """Give each auto pose its own orientation, or the first one when same is set."""
    taken = {(p.orientation, p.flip) for p in poses if p.orientation != "auto" and p.spin == 0}
    free = [o for o in AUTO_ORIENTATIONS if o not in taken] or AUTO_ORIENTATIONS
    for number, pose in enumerate(p for p in poses if p.orientation == "auto"):
        if same:
            number = 0
        pose.orientation, pose.flip = free[number % len(free)]
        pose.spin += 45.0 * (number // len(free))


def load_addbinder_config(path, check_paths=True, output_root=None):
    try:
        raw = yaml.load(Path(path).read_text(), Loader=_UniqueKeySafeLoader)
    except (OSError, yaml.YAMLError) as exc:
        raise ConfigError(str(exc)) from exc
    if not isinstance(raw, dict) or not isinstance(raw.get("addbinder"), dict):
        raise ConfigError("addbinder YAML must contain an addbinder: mapping")
    from config import _resolve_path_references
    raw = _resolve_path_references(raw)
    raw.pop("paths", None)
    if set(raw) - {"addbinder"}:
        raise ConfigError(f"Unknown addbinder sections: {sorted(set(raw) - {'addbinder'})}")
    raw = raw["addbinder"]
    unknown = set(raw) - set(AddBinderSpec.__dataclass_fields__)
    if unknown:
        raise ConfigError(f"Unknown addbinder options: {sorted(unknown)}")
    for key in ("host", "binder_coordinates"):
        if not isinstance(raw.get(key), str) or not raw[key].strip():
            raise ConfigError(f"{key} must be a nonempty path")
    itps = raw.get("binder_itp")
    if isinstance(itps, str):
        itps = [itps]
    if (not isinstance(itps, list) or not itps
            or any(not isinstance(v, str) or not v.strip() for v in itps)):
        raise ConfigError("binder_itp must be a nonempty path or a nonempty list of paths")
    if len(set(itps)) != len(itps):
        raise ConfigError("binder_itp lists the same file twice")
    raw["binder_itp"] = itps
    for key in ("binder_forcefield", "host_gro"):
        if raw.get(key) == "":
            raw[key] = None  # an empty optional field in the GUI means "not set"
        value = raw.get(key)
        if value is not None and (not isinstance(value, str) or not value.strip()):
            raise ConfigError(f"{key} must be a file path")
    raw.setdefault("output_dir", AddBinderSpec.__dataclass_fields__["output_dir"].default)
    raw = anchor_config_paths({"addbinder": raw}, config_path=path,
                              output_root=output_root)["addbinder"]
    poses = raw.pop("poses", None)
    cfg = AddBinderSpec(**raw)
    if cfg.side not in ("upper", "lower"):
        raise ConfigError("side must be upper (+z) or lower (-z)")
    if cfg.lateral_reference not in ("tip", "host_center"):
        raise ConfigError("lateral_reference must be tip or host_center")
    if cfg.ion_side not in ("binder_side", "any"):
        raise ConfigError("ion_side must be binder_side or any")
    if isinstance(cfg.salt_concentration, str) and cfg.salt_concentration.strip() != "host":
        try:  # the GUI's text field
            cfg.salt_concentration = float(cfg.salt_concentration)
        except ValueError:
            raise ConfigError("salt_concentration must be host or a concentration in mol/L") from None
    elif isinstance(cfg.salt_concentration, str):
        cfg.salt_concentration = "host"
    if cfg.salt_concentration != "host" and (
            isinstance(cfg.salt_concentration, bool) or not isinstance(cfg.salt_concentration, (int, float))
            or not 0 <= cfg.salt_concentration < 5):
        raise ConfigError("salt_concentration must be host or a concentration in mol/L (0-5)")
    for key in ("distance", "min_image_gap", "min_membrane_gap", "water_clash_distance",
                "ion_exclusion_distance", "lipid_ion_exclusion_distance"):
        setattr(cfg, key, _number(getattr(cfg, key), key))
    cfg.target_net_charge = _number(cfg.target_net_charge, "target_net_charge", False)
    for key in ("random_poses", "random_seed"):
        if type(getattr(cfg, key)) is not int or getattr(cfg, key) < 0:
            raise ConfigError(f"{key} must be a whole number, 0 or more")
    if cfg.random_rotated is not None and (type(cfg.random_rotated) is not int
                                           or not 0 <= cfg.random_rotated <= cfg.random_poses):
        raise ConfigError("random_rotated must be a whole number from 0 to random_poses (empty: all)")
    if cfg.random_distance is None:
        cfg.random_distance = [cfg.distance, cfg.distance + 5.0]
    if not isinstance(cfg.random_distance, (list, tuple)) or len(cfg.random_distance) != 2:
        raise ConfigError("random_distance must be [min, max] in Angstrom")
    cfg.random_distance = [_number(v, "random_distance") for v in cfg.random_distance]
    if cfg.random_distance[0] > cfg.random_distance[1]:
        raise ConfigError("random_distance must be [min, max] with min <= max")
    cfg.random_max_tilt = _number(cfg.random_max_tilt, "random_max_tilt", False)
    if not 0 <= cfg.random_max_tilt < 90:
        raise ConfigError("random_max_tilt must be 0-90 degrees")
    cfg.random_min_angle = _number(cfg.random_min_angle, "random_min_angle", False)
    cfg.random_min_separation = _number(cfg.random_min_separation, "random_min_separation", False)
    if cfg.random_min_angle < 0 or cfg.random_min_separation < 0:
        raise ConfigError("random_min_angle and random_min_separation must not be negative")
    for key in ("neutralize", "same_orientation", "copy_run_inputs", "random_from_input_position"):
        if type(getattr(cfg, key)) is not bool:
            raise ConfigError(f"{key} must be true or false")
    if not isinstance(cfg.environment_overrides, dict) or any(
            not isinstance(k, str) or v not in ("water", "ion", "lipid", "sterol", "detergent", "solvent")
            for k, v in cfg.environment_overrides.items()):
        raise ConfigError("environment_overrides must map molecule names to environment categories")
    if poses is None:
        poses = [{}]
    if not isinstance(poses, list) or not poses:
        raise ConfigError("poses must be a nonempty list")
    cfg.poses = [_parse_pose(p, i, cfg.distance) for i, p in enumerate(poses)]
    _resolve_auto_orientations(cfg.poses, cfg.same_orientation)
    if len({p.name for p in cfg.poses}) != len(cfg.poses):
        raise ConfigError("pose names must be unique")
    taken = {p.name for p in cfg.poses} & set(random_pose_names(cfg))
    if taken:
        raise ConfigError(f"pose names {sorted(taken)} are used by the random poses; rename them")
    output = Path(cfg.output_dir)
    host = Path(cfg.host)
    if output == host or host in output.parents or output in host.parents:
        raise ConfigError("output_dir must be separate from the host folder")
    if check_paths:
        if not host.is_dir():
            raise ConfigError(f"host folder does not exist: {host}")
        if not (host / "topol.top").is_file() or not (host / "toppar").is_dir():
            raise ConfigError(f"host must contain topol.top and toppar/: {host}")
        for key in ("binder_coordinates", "binder_forcefield", "host_gro"):
            value = getattr(cfg, key)
            if value and not Path(value).is_file():
                raise ConfigError(f"{key} does not exist: {value}")
        for value in cfg.binder_itp:
            if not Path(value).is_file():
                raise ConfigError(f"binder_itp does not exist: {value}")
    return cfg


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------

def _unwrap_sequential(positions, lengths):
    """Make one molecule whole by stepping atom to atom with minimum image.

    Consecutive ITP atoms are bonded or near-bonded, so each step is far below
    half a box length even for proteins that cross a periodic boundary.
    """
    out = positions.copy()
    for i in range(1, len(out)):
        delta = out[i] - out[i - 1]
        out[i] = out[i - 1] + delta - lengths * np.round(delta / lengths)
    return out


@dataclass
class MembraneFrame:
    lengths: np.ndarray
    shift: float              # added to z to put the membrane centre at Lz/2
    upper_plane: float        # centred-frame z of the upper headgroup plane
    lower_plane: float
    headgroup_atoms: str
    leaflet_spread: tuple

    def centred(self, positions):
        out = positions.copy()
        out[:, 2] = np.mod(out[:, 2] + self.shift, self.lengths[2])
        return out


def membrane_frame(universe, molecules):
    lengths = universe.dimensions[:3].astype(float)
    lipid = [m for m in molecules if m.category in ("lipid", "sterol")]
    if not lipid:
        raise ConfigError("No lipids found in the host; addbinder needs a membrane to define the sides")
    atoms = universe.atoms[_indices([m for m in lipid if m.category == "lipid"] or lipid)]
    heads = atoms[atoms.names == "P"]
    label = "lipid P atoms"
    if len(heads) < 10:
        heads, label = _heavy(atoms), "all lipid heavy atoms (no phosphorus found)"
    z = heads.positions[:, 2].astype(float)
    # Circular mean along z finds the bilayer centre even when it wraps.
    angle = 2 * np.pi * z / lengths[2]
    centre = np.mod(np.arctan2(np.sin(angle).mean(), np.cos(angle).mean()) * lengths[2] / (2 * np.pi),
                    lengths[2])
    shift = lengths[2] / 2 - centre
    zc = np.mod(z + shift, lengths[2])
    upper, lower = zc[zc > lengths[2] / 2], zc[zc <= lengths[2] / 2]
    if len(upper) < 5 or len(lower) < 5:
        raise ConfigError("Could not find two lipid leaflets along z; addbinder expects the membrane "
                          "normal along the box z axis")
    frame = MembraneFrame(lengths, float(shift), float(np.median(upper)), float(np.median(lower)),
                          label, (float(np.std(upper)), float(np.std(lower))))
    thickness = frame.upper_plane - frame.lower_plane
    if not 20.0 <= thickness <= 70.0:
        raise ConfigError(f"Headgroup planes are {thickness:.1f} A apart; expected a flat bilayer "
                          "normal to z (20-70 A)")
    return frame


def _orientation_matrix(pose, heavy_centred, side="upper"):
    if pose.orientation == "as_is":
        rotation = np.eye(3)
    elif pose.orientation == "euler":
        rotation = Rotation.from_euler("zyz", pose.angles, degrees=True).as_matrix()
    else:
        _, vectors = np.linalg.eigh(np.cov(heavy_centred.T))
        longest, middle, thinnest = vectors[:, 2], vectors[:, 1], vectors[:, 0]
        if pose.orientation == "flat":
            # Rows map binder axes onto the box: longest -> x, middle -> y, thinnest -> z.
            rotation = np.array([longest, middle, thinnest])
            if np.linalg.det(rotation) < 0:
                rotation[2] *= -1
        else:
            # Point the longest axis at the end reaching furthest from the centre
            # (positive third moment), e.g. a ligand's tail rather than its ring.
            if np.sum((heavy_centred @ longest) ** 3) < 0:
                longest = -longest
            if pose.orientation == "end_on":
                towards_host = -1.0 if side == "upper" else 1.0
                x, z = middle, towards_host * longest
            else:  # edge
                x, z = longest, middle
            rotation = np.array([x, np.cross(z, x), z])
    if pose.flip:
        rotation = np.diag([1.0, -1.0, -1.0]) @ rotation
    spin = Rotation.from_euler("z", pose.spin, degrees=True).as_matrix()
    return spin @ rotation


def _image_shifts(lengths):
    return [np.array(s) * lengths for s in product((-1, 0, 1), repeat=3) if any(s)]


def approach_direction(pose, side):
    """Unit vector from the tip out into the water: the membrane normal, or tilted by approach."""
    sign = 1.0 if side == "upper" else -1.0
    if pose.approach is None:
        return np.array([0.0, 0.0, sign])
    tilt, azimuth = np.radians(pose.approach)
    return np.array([np.sin(tilt) * np.cos(azimuth), np.sin(tilt) * np.sin(azimuth), sign * np.cos(tilt)])


def _slide_out(placed, heavy_mask, start, direction, host_heavy, distance):
    """Move the binder from start along direction until its nearest heavy atom is distance from the host."""
    tree = cKDTree(host_heavy)
    heavy = placed[heavy_mask]

    def gap(step):
        return float(tree.query(heavy + start + step * direction)[0].min())
    low, high = 0.0, 0.5
    while gap(high) < distance:
        low, high = high, high + 0.5
        if high > 1000:
            raise ConfigError("Could not move the binder clear of the host along the approach direction")
    for _ in range(40):
        middle = (low + high) / 2
        low, high = (middle, high) if gap(middle) < distance else (low, middle)
    return placed + start + high * direction


def place_binder(pose, cfg, frame, host_heavy, binder_positions, binder_heavy_mask, tip):
    """Return binder coordinates in the host's centred frame for one pose."""
    heavy = binder_positions[binder_heavy_mask]
    centre = heavy.mean(axis=0)
    rotation = _orientation_matrix(pose, heavy - centre, cfg.side)
    placed = (binder_positions - centre) @ rotation.T
    if pose.centroid is not None:
        target = np.array(pose.centroid) + [0.0, 0.0, frame.shift]
        return placed + target - placed[binder_heavy_mask].mean(axis=0), rotation
    if pose.from_input_position:
        # The input coordinates are in the host GRO frame; bring their centre into the
        # membrane-centred frame, at the periodic image nearest the host.
        start = centre + [0.0, 0.0, frame.shift]
        delta = start - host_heavy.mean(axis=0)
        start = start - frame.lengths * np.round(delta / frame.lengths)
        return _slide_out(placed, binder_heavy_mask, start, approach_direction(pose, cfg.side),
                          host_heavy, pose.distance), rotation
    if pose.distance_to == "centroid":
        # Spherical coordinates: the binder's centre at distance along the approach direction.
        origin = tip["position"] if cfg.lateral_reference == "tip" else host_heavy.mean(axis=0)
        target = origin + np.array([*pose.lateral_offset, 0.0]) + pose.distance * approach_direction(pose, cfg.side)
        return placed + target - placed[binder_heavy_mask].mean(axis=0), rotation
    if pose.approach is not None:
        start = tip["position"] + np.array([*pose.lateral_offset, 0.0])
        return _slide_out(placed, binder_heavy_mask, start, approach_direction(pose, cfg.side),
                          host_heavy, pose.distance), rotation
    placed_heavy = placed[binder_heavy_mask]
    reference = tip["position"][:2] if cfg.lateral_reference == "tip" else host_heavy[:, :2].mean(axis=0)
    target_xy = reference + np.asarray(pose.lateral_offset)
    placed[:, :2] += target_xy - placed_heavy[:, :2].mean(axis=0)
    if cfg.side == "upper":
        placed[:, 2] += tip["position"][2] + pose.distance - placed_heavy[:, 2].min()
    else:
        placed[:, 2] += tip["position"][2] - pose.distance - placed_heavy[:, 2].max()
    return placed, rotation


def check_pose(pose, cfg, frame, host_heavy, lipid_heavy, placed_heavy):
    """Distances to the host, its periodic images, the membrane and the binder's own images."""
    lengths = frame.lengths
    tree = cKDTree(host_heavy)
    primary = float(tree.query(placed_heavy)[0].min())
    image_gap, image_shift = np.inf, None
    for shift in _image_shifts(lengths):
        d = float(tree.query(placed_heavy - shift)[0].min())
        if d < image_gap:
            image_gap, image_shift = d, shift
    own = cKDTree(placed_heavy)
    self_gap = min(float(own.query(placed_heavy - s)[0].min()) for s in _image_shifts(lengths))
    membrane_gap = float("inf")
    if len(lipid_heavy):
        membrane_gap = float(cKDTree(lipid_heavy).query(placed_heavy)[0].min())
        for shift in _image_shifts(lengths):
            membrane_gap = min(membrane_gap, float(cKDTree(lipid_heavy + shift).query(placed_heavy)[0].min()))
    solute_z = np.concatenate([host_heavy[:, 2], lipid_heavy[:, 2]]) if len(lipid_heavy) else host_heavy[:, 2]
    free_column = float(lengths[2] - np.ptp(solute_z))
    height = float(np.ptp(placed_heavy[:, 2]))
    extent = np.ptp(placed_heavy, axis=0)
    checks = {
        "gap_to_host_angstrom": round(primary, 3),
        "gap_to_host_image_angstrom": round(image_gap, 3),
        "closest_host_image_cell": [int(round(v)) for v in image_shift / lengths],
        "gap_to_own_image_angstrom": round(self_gap, 3),
        "gap_to_membrane_angstrom": round(membrane_gap, 3),
        "binder_extent_angstrom": [round(float(v), 2) for v in extent],
        "binder_height_along_normal_angstrom": round(height, 2),
        "free_water_column_angstrom": round(free_column, 2),
        "column_room_for_binder_angstrom": round(free_column - pose.distance - cfg.min_image_gap, 2),
    }
    failures = []
    # GRO stores 0.01 A precision; allow that rounding at the requested limits.
    if pose.centroid is not None or pose.distance_to == "centroid":
        if primary < CENTROID_MIN_GAP:
            failures.append(f"binder is {primary:.2f} A from the host at the given centroid; it overlaps the "
                            f"host (under {CENTROID_MIN_GAP:.0f} A)")
    elif primary < pose.distance - 0.02:
        failures.append(f"binder is {primary:.2f} A from the host, closer than distance {pose.distance:.2f} A")
    if image_gap < cfg.min_image_gap - 0.02:
        cell = checks["closest_host_image_cell"]
        message = (f"binder is {image_gap:.2f} A from a periodic image of the host (cell {cell}); "
                   f"min_image_gap is {cfg.min_image_gap:.2f} A")
        if cell[2] != 0 and cell[0] == 0 and cell[1] == 0:
            # The binder must fit in the free column, not just clear the image by the missing gap.
            shortfall = max(cfg.min_image_gap - image_gap, height - checks["column_room_for_binder_angstrom"])
            message += (f". About {shortfall:.1f} A more box height (water) is "
                        "needed, or a smaller distance / flatter orientation")
        else:
            message += ". Move the binder laterally or use a larger box in x/y"
        failures.append(message)
    if self_gap < cfg.min_image_gap - 0.02:
        failures.append(f"binder is {self_gap:.2f} A from its own periodic image; min_image_gap is "
                        f"{cfg.min_image_gap:.2f} A (box too small for this binder)")
    if membrane_gap < cfg.min_membrane_gap - 0.02:
        failures.append(f"binder is {membrane_gap:.2f} A from lipid heavy atoms; min_membrane_gap is "
                        f"{cfg.min_membrane_gap:.2f} A")
    return checks, failures


# ---------------------------------------------------------------------------
# Binder input
# ---------------------------------------------------------------------------

def random_pose_names(cfg):
    return [f"binderpose{len(cfg.poses) + number + 1}" for number in range(cfg.random_poses)]


def _random_rotation_angles(rng):
    quaternion = rng.normal(size=4)
    return Rotation.from_quat(quaternion / np.linalg.norm(quaternion)).as_euler("zyz", degrees=True)


def random_poses(cfg, frame, host_heavy, lipid_heavy, binder_positions, heavy_mask, tip):
    """Draw cfg.random_poses extra poses: random approach direction, distance and orientation.

    A draw is redrawn when it fails the box/membrane checks, comes from within
    random_min_angle of another pose's direction, or puts the binder centre
    within random_min_separation of another pose's. Returns (pose, failures)
    pairs; failures is empty for a pose that was found.
    """
    rng = np.random.default_rng(cfg.random_seed)
    taken = []  # (direction, binder heavy-atom centre) of every pose so far
    for pose in cfg.poses:
        placed, _ = place_binder(pose, cfg, frame, host_heavy, binder_positions, heavy_mask, tip)
        taken.append((approach_direction(pose, cfg.side), placed[heavy_mask].mean(axis=0)))
    low, high = cfg.random_distance
    min_cos = np.cos(np.radians(cfg.random_max_tilt))
    result = []
    for number, name in enumerate(random_pose_names(cfg)):
        rotated = cfg.random_rotated is None or number >= cfg.random_poses - cfg.random_rotated
        reasons = {"box or membrane checks": 0, "too close in angle": 0, "too close in position": 0}
        for attempt in range(1, RANDOM_ATTEMPTS + 1):
            tilt = float(np.degrees(np.arccos(rng.uniform(min_cos, 1.0))))  # uniform over the cap
            angles = [round(float(v), 3) for v in _random_rotation_angles(rng)]
            pose = PoseSpec(name=name, distance=round(float(rng.uniform(low, high)), 2),
                            orientation="euler" if rotated else "as_is",
                            angles=angles if rotated else [0.0, 0.0, 0.0],
                            approach=[round(tilt, 3), round(float(rng.uniform(0.0, 360.0)), 3)],
                            from_input_position=cfg.random_from_input_position)
            placed, _ = place_binder(pose, cfg, frame, host_heavy, binder_positions, heavy_mask, tip)
            direction, centre = approach_direction(pose, cfg.side), placed[heavy_mask].mean(axis=0)
            if any(np.degrees(np.arccos(np.clip(direction @ d, -1, 1))) < cfg.random_min_angle
                   for d, _ in taken):
                reasons["too close in angle"] += 1
                continue
            if any(np.linalg.norm(centre - c) < cfg.random_min_separation for _, c in taken):
                reasons["too close in position"] += 1
                continue
            if check_pose(pose, cfg, frame, host_heavy, lipid_heavy, placed[heavy_mask])[1]:
                reasons["box or membrane checks"] += 1
                continue
            taken.append((direction, centre))
            result.append((pose, [], attempt))
            break
        else:
            counts = ", ".join(f"{k} {v}" for k, v in reasons.items() if v)
            result.append((PoseSpec(name=name, distance=low), [
                f"no random placement passed after {RANDOM_ATTEMPTS} draws ({counts}); widen "
                "random_distance or random_max_tilt, or lower random_min_angle / random_min_separation"],
                RANDOM_ATTEMPTS))
    return result


@dataclass
class BinderSystem:
    files: list
    universe: object
    molecules: list
    category: str
    report: dict


@dataclass
class BinderInput:
    """The binder as read from its inputs, before any pose is applied."""
    names: list           # moleculetype names, in coordinate order
    definitions: list     # ITPMoleculeType per molecule
    categories: list      # "protein" or "ligand" per molecule
    positions: np.ndarray
    files: list           # force-field includes followed by each molecule's ITP
    report: dict

    @property
    def atoms(self):
        return [a for d in self.definitions for a in d.atoms]

    @property
    def charge(self):
        return sum(d.charge for d in self.definitions)

    @property
    def category(self):
        return "protein" if "protein" in self.categories else "ligand"


def _rename_moleculetype(source, name, renamed, workdir):
    text = Path(source).read_text()
    text, count = re.subn(r"(\[\s*moleculetype\s*\][^\n]*\n(?:\s*;[^\n]*\n)*\s*)" + re.escape(name) + r"\b",
                          lambda m: m[1] + renamed, text, count=1, flags=re.I)
    if count != 1:
        raise ConfigError(f"Could not rename binder moleculetype {name}")
    target = Path(workdir) / f"{renamed}.itp"
    target.write_text(text)
    return parse_itp(str(target))[renamed]


def _binder_definitions(cfg, workdir, host):
    """One (name, definition, renamed_from) per binder ITP, in the order given."""
    result = []
    for path in cfg.binder_itp:
        definitions = parse_itp(path)
        if len(definitions) != 1:
            raise ConfigError(f"Each binder_itp must define exactly one moleculetype; {path} defines "
                              f"{sorted(definitions)}")
        name, definition = next(iter(definitions.items()))
        host_definition = host.definitions.get(name)
        renamed = None
        if host_definition is not None and (
                Path(host_definition.source_path).read_bytes() != Path(path).read_bytes()):
            # CHARMM-GUI names every protein chain PROA, PROB ...; keep both by renaming the binder's type.
            renamed = f"{name}_BND"
            definition = _rename_moleculetype(path, name, renamed, workdir)
            name, renamed = renamed, name
        result.append((name, definition, renamed))
    names = [name for name, _, _ in result]
    duplicates = sorted({n for n in names if names.count(n) > 1})
    if duplicates:
        raise ConfigError(f"binder_itp defines moleculetype {duplicates} more than once; give each "
                          "binder molecule its own type name")
    return result


def load_binder(cfg, workdir, host):
    entries = _binder_definitions(cfg, workdir, host)
    definitions = [d for _, d, _ in entries]
    label = "+".join(name for name, _, _ in entries)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        source = mda.Universe(cfg.binder_coordinates)
    atoms = source.atoms
    expected = [a.atom_name for d in definitions for a in d.atoms]
    observed = [str(n) for n in atoms.names]
    notes = []
    if len(observed) != len(expected):
        sizes = ", ".join(f"{name} {d.atom_count}" for name, d, _ in entries)
        raise ConfigError(f"binder_coordinates has {len(observed)} atoms but binder_itp defines "
                          f"{len(expected)} ({sizes}); supply each binder molecule once, in binder_itp "
                          "order, hydrogens included")
    if observed != expected:
        if len(entries) == 1 and sorted(observed) == sorted(expected) and len(set(expected)) == len(expected):
            order = {n: i for i, n in enumerate(observed)}
            atoms = atoms[[order[n] for n in expected]]
            notes.append("binder atoms reordered by name to match the ITP")
        else:
            mismatch = next(i for i, (a, b) in enumerate(zip(observed, expected)) if a != b)
            offset, owner = mismatch, entries[0][0]
            for name, d, _ in entries:
                if offset < d.atom_count:
                    owner = name
                    break
                offset -= d.atom_count
            raise ConfigError(f"binder atom {mismatch + 1} is {observed[mismatch]} in the coordinates but "
                              f"{expected[mismatch]} in the ITP ({owner} atom {offset + 1}); atom names "
                              "and order must match binder_itp")
    positions = atoms.positions.astype(float)
    if np.ptp(positions, axis=0).max() > 250:
        raise ConfigError("binder coordinates span more than 250 A; make the binder whole first")
    categories = []
    for name, d, _ in entries:
        is_protein = (re.fullmatch(r"PRO[A-Z0-9_]+", name.upper()) is not None or sum(
            a.atom_name == "CA" and a.resname in AMINO_ACID_RESNAMES for a in d.atoms) >= 2)
        categories.append("protein" if is_protein else "ligand")
    files = []
    if cfg.binder_forcefield:
        forcefield = Path(cfg.binder_forcefield).resolve()
        files.extend(_included_files(forcefield, forcefield.parent))
    files.extend(Path(d.source_path) for d in definitions)
    renamed = {name: old for name, _, old in entries if old}
    report = {"moleculetype": label, "renamed_from": renamed or None,
              "category": "protein" if "protein" in categories else "ligand",
              "molecules": [{"moleculetype": name, "category": c, "atoms": d.atom_count,
                             "net_charge": d.charge, "itp": itp}
                            for (name, d, _), c, itp in zip(entries, categories, cfg.binder_itp)],
              "atoms": sum(d.atom_count for d in definitions),
              "net_charge": round(sum(d.charge for d in definitions), 6),
              "coordinates": cfg.binder_coordinates, "itp": cfg.binder_itp,
              "forcefield": cfg.binder_forcefield, "notes": notes}
    return BinderInput([n for n, _, _ in entries], definitions, categories, positions, files, report)


def _binder_universe(definitions, positions, workdir):
    resindex, residues, previous = [], [], None
    for number, definition in enumerate(definitions):
        for atom in definition.atoms:
            key = (number, atom.resnr)
            if key != previous:
                residues.append((atom.resnr, atom.resname))
                previous = key
            resindex.append(len(residues) - 1)
    universe = mda.Universe.empty(len(resindex), n_residues=len(residues),
                                  atom_resindex=resindex, trajectory=True)
    universe.add_TopologyAttr("names", [a.atom_name for d in definitions for a in d.atoms])
    universe.add_TopologyAttr("resnames", [r[1] for r in residues])
    universe.add_TopologyAttr("resids", [r[0] for r in residues])
    universe.atoms.positions = positions
    universe.dimensions = [100.0, 100.0, 100.0, 90.0, 90.0, 90.0]
    path = Path(workdir) / "binder_placed.gro"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        universe.atoms.write(str(path))
        # Re-read through the GRO parser so the binder merges with the GRO host.
        reread = mda.Universe(str(path))
    if len(reread.residues) != len(residues):
        raise ConfigError("Two consecutive binder molecules meet at residues with the same number and "
                          "name, which GRO cannot keep apart; renumber one of them")
    return reread


# ---------------------------------------------------------------------------
# Environment cleanup and charge
# ---------------------------------------------------------------------------

def _remove_overlapping_solvent(host, environment, binder_positions, cfg):
    universe = host.universe
    env = _heavy(universe.atoms[_indices(environment)])
    pairs, distances = capped_distance(binder_positions, env.positions, cfg.water_clash_distance,
                                       box=universe.dimensions, return_distances=True)
    owner = np.full(len(universe.atoms), -1, dtype=int)
    for i, mol in enumerate(environment):
        owner[mol.start:mol.stop] = i
    close = {}
    for (_, j), distance in zip(pairs, distances):
        i = int(owner[env.indices[j]])
        close[i] = min(close.get(i, np.inf), float(distance))
    kept, removed = [], []
    for i, mol in enumerate(environment):
        if i not in close:
            kept.append(mol)
        elif mol.category in ("water", "ion", "solvent"):
            removed.append({"moleculetype": mol.name, "occurrence": mol.occurrence,
                            "category": mol.category, "charge": mol.definition.charge,
                            "minimum_distance_angstrom": round(close[i], 3)})
        else:
            raise ConfigError(f"{mol.category} {mol.name}:{mol.occurrence} lies {close[i]:.2f} A from the "
                              "binder; only water and ions may be removed")
    return kept, removed


WATER_MOLARITY = 55.51   # mol/L of pure water: salt concentration = pairs / waters * 55.51
ION_MIN_SPACING = 5.0    # A between an added ion and any other ion


def _salt_species(molecules):
    """The main monovalent cation and anion type of a system, by count."""
    counts = {}
    for mol in molecules:
        if mol.category == "ion" and abs(abs(mol.definition.charge) - 1) < 0.01:
            key = (mol.name, round(mol.definition.charge))
            counts[key] = counts.get(key, 0) + 1
    cation = max((k for k in counts if k[1] > 0), key=counts.get, default=None)
    anion = max((k for k in counts if k[1] < 0), key=counts.get, default=None)
    return (cation and cation[0]), (anion and anion[0])


def host_salt(host):
    """Salt pairs, waters and concentration (M) of the host as built."""
    cation, anion = _salt_species(host.molecules)
    waters = sum(m.category == "water" for m in host.molecules)
    n_cat = sum(m.name == cation for m in host.molecules)
    n_an = sum(m.name == anion for m in host.molecules)
    pairs = min(n_cat, n_an)
    return {"cation": cation, "anion": anion, "cations": n_cat, "anions": n_an, "waters": waters,
            "pairs": pairs, "concentration_molar": (pairs * WATER_MOLARITY / waters) if waters else 0.0}


def _ion_universe(definitions, positions, workdir):
    """Single-atom ions as their own small universe, one residue each, read back through GRO."""
    n = len(definitions)
    universe = mda.Universe.empty(n, n_residues=n, atom_resindex=np.arange(n), trajectory=True)
    universe.add_TopologyAttr("names", [d.atoms[0].atom_name for d in definitions])
    universe.add_TopologyAttr("resnames", [d.atoms[0].resname for d in definitions])
    universe.add_TopologyAttr("resids", np.arange(1, n + 1))
    universe.atoms.positions = positions
    path = Path(workdir) / "added_ions.gro"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        universe.atoms.write(str(path))
        reread = mda.Universe(str(path))
    path.unlink()
    return reread


def _adjust_salt(host, kept, solute, binder_charge, binder_universe, frame, cfg, stage, seed):
    """Set the final counter- and co-ion counts: the host's salt concentration on the remaining
    water, plus the counterions that bring the system to target_net_charge.

    Extra ions are removed and missing ones replace bulk waters, at random among
    molecules clear of the protein, the binder and the lipids. Returns the kept
    environment, the added ions as (molecule, system) pairs, and the report.
    """
    universe = host.universe
    salt = host_salt(host)
    cation, anion = salt["cation"], salt["anion"]
    if cfg.salt_concentration == "host":
        concentration = salt["concentration_molar"]
    else:
        concentration = float(cfg.salt_concentration)
    definitions = {m.name: m.definition for m in host.molecules if m.name in (cation, anion)}
    is_salt = lambda m: m.name in (cation, anion)
    other_charge = (sum(m.definition.charge for m in solute) + binder_charge
                    + sum(m.definition.charge for m in kept if not is_salt(m)))
    n_cat = sum(m.name == cation for m in kept)
    n_an = sum(m.name == anion for m in kept)
    before = other_charge + n_cat - n_an
    if cfg.neutralize:
        excess = cfg.target_net_charge - other_charge
        if abs(excess - round(excess)) > 0.01:
            raise ConfigError(f"The system without salt ions carries {other_charge:.3f} e; a whole number of "
                              f"monovalent ions cannot reach the target charge {cfg.target_net_charge}")
        excess = int(round(excess))  # cations minus anions needed
    else:
        excess = n_cat - n_an
    report = {"species": [cation, anion], "host": salt, "target_concentration_molar": round(concentration, 5),
              "concentration_source": "host" if cfg.salt_concentration == "host" else "salt_concentration",
              "before": round(before, 6), "target": cfg.target_net_charge, "enabled": cfg.neutralize,
              "ion_side": cfg.ion_side, "removed_ions": [], "added_ions": []}

    waters = sum(m.category == "water" for m in kept)
    added_total = 0
    for _ in range(3):  # added ions replace waters, so settle the pair count on the final water count
        pairs = int(round(concentration * (waters - added_total) / WATER_MOLARITY))
        targets = {cation: pairs + max(excess, 0), anion: pairs + max(-excess, 0)}
        added_total = max(targets[cation] - n_cat, 0) + max(targets[anion] - n_an, 0)
    report.update({"waters_before_adjustment": waters, "target_pairs": pairs,
                   "target_counts": dict(targets), "counts_before_adjustment": {cation: n_cat, anion: n_an}})
    if (cation is None or anion is None) and (pairs or excess):
        raise ConfigError("The host has no monovalent cation or anion to copy; set salt_concentration: 0 and "
                          "neutralize: false, or use a host built with salt")

    rng = np.random.default_rng(seed)
    box = universe.dimensions
    solute_heavy = np.concatenate([_heavy(universe.atoms[_indices(solute)]).positions,
                                   _heavy(binder_universe.atoms).positions])
    lipid = [m for m in kept if m.category in ("lipid", "sterol", "detergent")]
    lipid_heavy = _heavy(universe.atoms[_indices(lipid)]).positions if lipid else np.zeros((0, 3))

    def clear(points, extra=()):
        """Mask of points clear of the solute, the lipids and any extra (positions, cutoff) groups."""
        ok = np.ones(len(points), bool)
        for group, cutoff in ((solute_heavy, cfg.ion_exclusion_distance),
                              (lipid_heavy, cfg.lipid_ion_exclusion_distance), *extra):
            if len(group) and len(points):
                pairs_found = capped_distance(points, group, cutoff, box=box, return_distances=False)
                ok[np.unique(pairs_found[:, 0])] = False
        return ok

    def on_binder_side(points):
        upper = frame.centred(points)[:, 2] > frame.lengths[2] / 2
        return upper == (cfg.side == "upper")

    def preference(indices, points):
        """Random order; binder-side molecules first when ion_side is binder_side."""
        order = rng.permutation(len(indices))
        if cfg.ion_side == "binder_side":
            side = on_binder_side(points)
            order = np.concatenate([order[side[order]], order[~side[order]]])
        return [indices[i] for i in order]

    # Remove extra ions.
    drop = set()
    for name, target in targets.items():
        current = sum(m.name == name for m in kept)
        if current <= target:
            continue
        candidates = [i for i, m in enumerate(kept) if m.name == name]
        points = universe.atoms[[kept[i].start for i in candidates]].positions
        eligible = [c for c, ok in zip(candidates, clear(points)) if ok]
        if len(eligible) < current - target:
            raise ConfigError(f"Only {len(eligible)} {name} ions are clear of the protein and lipids; "
                              f"{current - target} must be removed")
        points = universe.atoms[[kept[i].start for i in eligible]].positions
        for i in preference(eligible, points)[:current - target]:
            drop.add(i)
            report["removed_ions"].append({"moleculetype": name, "occurrence": kept[i].occurrence,
                                           "charge": kept[i].definition.charge,
                                           "side": "binder_side" if on_binder_side(
                                               universe.atoms[[kept[i].start]].positions)[0] else "other_side"})

    # Add missing ions in place of bulk waters.
    to_add = [(name, target - sum(m.name == name for m in kept)) for name, target in targets.items()]
    to_add = [(name, n) for name, n in to_add if n > 0]
    sites, site_names = [], []
    if to_add:
        water_ids = [i for i, m in enumerate(kept) if m.category == "water" and i not in drop]
        oxygen = universe.atoms[[kept[i].start for i in water_ids]].positions
        ion_points = universe.atoms[[m.start for i, m in enumerate(kept)
                                     if m.category == "ion" and i not in drop]].positions
        ok = clear(oxygen, [(ion_points, ION_MIN_SPACING)])
        candidates = [w for w, good in zip(water_ids, ok) if good]
        points = {w: p for w, p in zip(water_ids, oxygen)}
        needed = [name for name, n in to_add for _ in range(n)]
        rng.shuffle(needed)
        for w in preference(candidates, np.array([points[w] for w in candidates]) if candidates else
                            np.zeros((0, 3))):
            if len(sites) == len(needed):
                break
            p = points[w]
            if sites and len(capped_distance(p[None], np.array(sites), ION_MIN_SPACING, box=box,
                                             return_distances=False)):
                continue
            sites.append(p)
            site_names.append(needed[len(sites) - 1])
            drop.add(w)
            report["added_ions"].append({"moleculetype": site_names[-1], "replaced_water": kept[w].occurrence,
                                         "side": "binder_side" if on_binder_side(p[None])[0] else "other_side"})
        if len(sites) < len(needed):
            raise ConfigError(f"Only {len(sites)} of {len(needed)} ions could be placed in bulk water clear of "
                              "the protein, lipids and other ions")

    kept = [m for i, m in enumerate(kept) if i not in drop]
    added = []
    if sites:
        ion_universe = _ion_universe([definitions[n] for n in site_names], np.array(sites), stage)
        system = BuildSystemStub(ion_universe)
        added = [(system, Molecule(n, k + 1, k, k + 1, "ion", True, definitions[n]))
                 for k, n in enumerate(site_names)]
    final_cat = sum(m.name == cation for m in kept) + site_names.count(cation)
    final_an = sum(m.name == anion for m in kept) + site_names.count(anion)
    final_waters = sum(m.category == "water" for m in kept)
    final_pairs = min(final_cat, final_an)
    report.update({"final_counts": {cation: final_cat, anion: final_an}, "final_waters": final_waters,
                   "final_concentration_molar": round(final_pairs * WATER_MOLARITY / final_waters, 5)
                   if final_waters else 0.0,
                   "after": round(other_charge + final_cat - final_an, 6)})
    return kept, added, report


@dataclass
class BuildSystemStub:
    """Holds the coordinates of molecules addbinder creates (added ions)."""
    universe: object


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

GENERATED = {"step5_input.pdb", "step5_input.gro", "topol.top", "toppar", "index.ndx",
             "addbinder_report.json", "addbinder_report.txt", "README"}


def _run_inputs(host_folder):
    """The host's CHARMM-GUI run inputs: step6/step7 mdp files and the README run script."""
    folder = Path(host_folder)
    return sorted(p for p in folder.iterdir()
                  if p.is_file() and (p.suffix == ".mdp" or p.name == "README"))


def _publish(stage, destination):
    """Replace a pose folder only when it is empty or was written by addbinder."""
    if destination.exists():
        if not destination.is_dir():
            raise ConfigError(f"{destination} exists and is not a folder")
        foreign = {p.name for p in destination.iterdir() if p.suffix != ".mdp"} - GENERATED
        if foreign:
            raise ConfigError(f"{destination} contains files addbinder did not write: {sorted(foreign)}")
        shutil.rmtree(destination)
    stage.rename(destination)


def _describe_pose(pose):
    text = pose["orientation"] + (" flipped" if pose["flip"] else "")
    if pose["orientation"] == "euler":
        text += f" {pose['angles']}"
    if pose["spin"]:
        text += f", spin {pose['spin']:g}"
    if pose.get("approach"):
        text += f", approach tilt {pose['approach'][0]:g} azimuth {pose['approach'][1]:g}"
    if pose.get("centroid"):
        text += f", centroid {pose['centroid']} A"
    if pose.get("from_input_position"):
        text += ", moved out from its input position"
    return text


def _write_report(folder, report):
    (folder / "addbinder_report.json").write_text(json.dumps(report, indent=2) + "\n")
    lines = [f"GMXtransplant addbinder: pose {report['pose']['name']}", f"Status: {report['status']}"]
    if report.get("failures"):
        lines += ["", "Rejected because:"] + [f"  - {f}" for f in report["failures"]]
    if "frame" in report:
        f = report["frame"]
        lines += ["", f"Membrane: upper/lower headgroup planes {f['upper_plane_z']:.2f} / {f['lower_plane_z']:.2f} A "
                      f"(original frame), thickness {f['thickness_angstrom']:.2f} A, from {f['headgroup_atoms']}"]
    if "tip" in report:
        t = report["tip"]
        lines.append(f"Tip ({report['side']} side): {t['moleculetype']} {t['resname']}{t['resid']} {t['atom']} "
                     f"at z = {t['z_original_frame']:.2f} A, {t['beyond_headgroup_plane_angstrom']:.2f} A beyond "
                     "the headgroup plane")
    if "binder" in report:
        b = report["binder"]
        lines.append(f"Binder: {b['moleculetype']} ({b['category']}), {b['atoms']} atoms, net charge "
                     f"{b['net_charge']:+.3f} e")
    pose = report["pose"]
    if pose.get("centroid"):
        lines.append(f"Pose: {_describe_pose(pose)} (host GRO frame)")
    else:
        if pose.get("distance_to") == "centroid":
            where = (f"from the {'tip' if report['settings']['lateral_reference'] == 'tip' else 'host centre'} "
                     "to the binder's centre")
        elif pose.get("approach") or pose.get("from_input_position"):
            where = "from the nearest host heavy atom"
        else:
            where = "beyond the tip"
        lines.append(f"Pose: {_describe_pose(pose)}, {pose['distance']:.2f} A {where}, lateral offset "
                     f"{pose['lateral_offset']}")
        if report.get("distance_attempts"):
            lines.append("Distance reduced until the pose fit: " + "; ".join(
                f"{a['distance']:g} A " + ("fits" if not a["failures"] else "does not fit")
                for a in report["distance_attempts"]))
    if report.get("random"):
        r = report["random"]
        lines.append(f"Random pose: seed {r['seed']}, accepted on draw {r['draw']}")
    if "checks" in report:
        lines += ["", "Distances (heavy atoms, Angstrom):"]
        lines += [f"  {k}: {v}" for k, v in report["checks"].items()]
    if "solvent_removed" in report:
        from collections import Counter
        counts = Counter(item["moleculetype"] for item in report["solvent_removed"])
        lines += ["", "Removed for overlap with the binder (within water_clash_distance): " +
                  (", ".join(f"{n} x{c}" for n, c in counts.items()) or "none")]
    salt = report.get("charge")
    if salt and "final_counts" in salt:
        from collections import Counter
        removed = Counter(x["moleculetype"] for x in salt["removed_ions"])
        added = Counter(x["moleculetype"] for x in salt["added_ions"])
        lines += ["", f"Salt: {salt['target_concentration_molar']:.4f} M ({salt['concentration_source']}; host "
                      f"{salt['host']['concentration_molar']:.4f} M), final {salt['final_concentration_molar']:.4f} M "
                      f"on {salt['final_waters']} waters; {salt['species'][0]}/{salt['species'][1]} "
                      f"{salt['counts_before_adjustment']} -> {salt['final_counts']}",
                  "  removed: " + (", ".join(f"{n} x{c}" for n, c in removed.items()) or "none")
                  + "; added in place of bulk water: " + (", ".join(f"{n} x{c}" for n, c in added.items()) or "none")
                  + f"; net charge {salt['before']:+.3f} -> {salt['after']:+.3f} e"]
    for key in ("charge", "final_molecules"):
        if key in report:
            lines += ["", f"{key}:", json.dumps(report[key], indent=2)]
    if report.get("run_inputs_copied"):
        lines += ["", "Copied from the host to run this pose: " + ", ".join(report["run_inputs_copied"])]
    if "error" in report:
        lines += ["", f"Error: {report['error']}"]
    (folder / "addbinder_report.txt").write_text("\n".join(lines) + "\n")


def _write_index(path, ordered, binder_system, tip_index):
    groups = {name: [] for name in ("System", "Protein", "Ligand", "Host", "Binder", "Host_tip",
                                    "Solute", "Water", "Ion", "Membrane", "Environment",
                                    "SOLU", "MEMB", "SOLV", "SOLU_MEMB", "SYSTEM")}
    # CHARMM-GUI groups, so the host's step6/step7 mdp files run unchanged.
    charmm_gui = {"protein": ("SOLU", "SOLU_MEMB"), "ligand": ("SOLU", "SOLU_MEMB"),
                  "lipid": ("MEMB", "SOLU_MEMB"), "sterol": ("MEMB", "SOLU_MEMB"),
                  "detergent": ("MEMB", "SOLU_MEMB"), "water": ("SOLV",), "ion": ("SOLV",),
                  "solvent": ("SOLV",)}
    offset = 0
    for system, mol in ordered:
        ids = list(range(offset + 1, offset + 1 + mol.stop - mol.start))
        groups["System"].extend(ids)
        groups["SYSTEM"].extend(ids)
        for name in charmm_gui.get(mol.category, ()):
            groups[name].extend(ids)
        label = {"protein": "Protein", "ligand": "Ligand", "water": "Water", "ion": "Ion",
                 "lipid": "Membrane", "sterol": "Membrane", "detergent": "Membrane"}.get(mol.category)
        if label:
            groups[label].extend(ids)
        if mol.category in ("protein", "ligand"):
            groups["Solute"].extend(ids)
            groups["Binder" if system is binder_system else "Host"].extend(ids)
            if system is not binder_system and mol.start <= tip_index < mol.stop:
                groups["Host_tip"].append(offset + 1 + tip_index - mol.start)
        else:
            groups["Environment"].extend(ids)
        offset += len(ids)
    Path(path).write_text("\n".join(
        f"[ {name} ]\n" + "\n".join(" ".join(map(str, ids[i:i + 15])) for i in range(0, len(ids), 15))
        for name, ids in groups.items()) + "\n")
    return groups


def _build_pose(pose, cfg, host, frame, solute, environment, host_heavy,
                lipid_heavy, tip, binder, stage, report):
    binder_report = binder.report
    heavy_mask = np.array([_is_heavy(a) for a in binder.atoms])
    distances = [pose.distance]
    if pose.reduce_distance_by > 0:
        while distances[-1] - pose.reduce_distance_by >= pose.min_distance - 1e-9:
            distances.append(round(distances[-1] - pose.reduce_distance_by, 3))
    attempts = []
    for distance in distances:
        trial = replace(pose, distance=distance)
        placed, rotation = place_binder(trial, cfg, frame, host_heavy, binder.positions, heavy_mask, tip)
        checks, failures = check_pose(trial, cfg, frame, host_heavy, lipid_heavy, placed[heavy_mask])
        attempts.append({"distance": distance, "failures": failures})
        if not failures:
            break
    if len(distances) > 1:
        report["distance_attempts"] = attempts
        report["distance_used"] = distance
        pose.distance = distance
        report["pose"]["distance"] = distance
    report["checks"] = checks
    report["rotation_matrix"] = np.round(rotation, 6).tolist()
    if failures:
        report["status"], report["failures"] = "rejected", failures
        return None
    # Back to the host's original frame; keep the binder whole, centre inside the box.
    placed[:, 2] -= frame.shift
    centre = placed.mean(axis=0)
    placed -= frame.lengths * np.floor(centre / frame.lengths)
    binder_universe = _binder_universe(binder.definitions, placed, stage)
    (stage / "binder_placed.gro").unlink()
    binder_molecules, offset = [], 0
    for name, definition, category in zip(binder.names, binder.definitions, binder.categories):
        binder_molecules.append(Molecule(name, 1, offset, offset + definition.atom_count, category, True,
                                         definition))
        offset += definition.atom_count
    binder_system = BinderSystem(binder.files, binder_universe, binder_molecules, binder.category, binder_report)

    kept, report["solvent_removed"] = _remove_overlapping_solvent(host, environment, placed, cfg)
    seed = cfg.random_seed * 1000003 + zlib.crc32(pose.name.encode())
    kept, added, report["charge"] = _adjust_salt(host, kept, solute, binder.charge, binder_universe, frame,
                                                 cfg, stage, seed)
    environment_order = [(host, m) for m in kept]
    for system, mol in added:  # each added ion joins the end of its type's block
        last = max((i for i, (_, m) in enumerate(environment_order) if m.name == mol.name), default=None)
        if last is None:
            last = max((i for i, (_, m) in enumerate(environment_order) if m.category == "ion"),
                       default=len(environment_order) - 1)
        environment_order.insert(last + 1, (system, mol))

    # Host proteins, binder proteins, host ligands, binder ligands, then environment.
    ordered = ([(host, m) for m in solute if m.category == "protein"]
               + [(binder_system, m) for m in binder_molecules if m.category == "protein"]
               + [(host, m) for m in solute if m.category != "protein"]
               + [(binder_system, m) for m in binder_molecules if m.category != "protein"]
               + environment_order)
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

    result = write_topology(stage, host, binder_system, ordered)
    top = Path(result.top_path)
    top.write_text(top.read_text().replace("CHARMM-GUI protein transplant", "GMXtransplant addbinder"))
    audit = audit_final_topology(final, result)
    if abs(audit["net_charge"] - report["charge"]["after"]) > 0.01:
        raise ConfigError("Final topology charge disagrees with molecule charge accounting")
    protein_counts = [m.stop - m.start for _, m in ordered if m.category == "protein"]

    def validate_written(_pdb, gro):
        return audit_final_topology(mda.Universe(gro), result)

    write_outputs(final, str(stage / "step5_input.pdb"), str(stage / "step5_input.gro"),
                  protein_chain_atom_counts=protein_counts, staged_validator=validate_written)
    index_groups = _write_index(stage / "index.ndx", ordered, binder_system, tip["index"])
    if cfg.copy_run_inputs:
        for source in _run_inputs(cfg.host):
            shutil.copy2(source, stage / source.name)
        report["run_inputs_copied"] = [p.name for p in _run_inputs(cfg.host)]
    report["final_molecules"] = result.molecules_written
    report["final_net_charge"] = round(audit["net_charge"], 6)
    report["status"] = "accepted"
    return final, index_groups


def _is_heavy(atom):
    # Every CHARMM/CGenFF hydrogen type starts with H (HA, HGA1, HN1, HT ...);
    # no heavy-atom type does.
    return not atom.atom_type.upper().startswith("H")


# One colour per pose; the host keeps the shared palette colours.
POSE_COLOURS = [[.84, .22, .67], [.95, .55, .10], [.20, .70, .30], [.55, .30, .85], [.90, .20, .20],
                [.10, .60, .80], [.75, .65, .10], [.40, .40, .40]]


def _write_scene(destination, built):
    """One PyMOL/VMD scene: host and membrane from the first pose, every pose's binder."""
    from visualization import safe_scene
    first_universe, first_groups = built[0][1]

    def atoms(universe, ids):
        return universe.atoms[np.array(ids, dtype=int) - 1] if ids else universe.atoms[[]]
    groups = {"retained_protein": atoms(first_universe, [i for i in first_groups["Host"]
                                                          if i in set(first_groups["Protein"])]),
              "retained_ligand": atoms(first_universe, [i for i in first_groups["Host"]
                                                         if i in set(first_groups["Ligand"])]),
              "retained_lipid": atoms(first_universe, first_groups["Membrane"]),
              "retained_water": atoms(first_universe, first_groups["Water"]),
              "retained_ion": atoms(first_universe, first_groups["Ion"]),
              "host_tip": atoms(first_universe, first_groups["Host_tip"])}
    palette = {"host_tip": ("Protein tip atom", [.98, .80, .10], True)}
    for number, (name, (universe, index_groups)) in enumerate(built):
        role = f"binder_{name}"
        groups[role] = atoms(universe, index_groups["Binder"])
        palette[role] = (f"Binder, pose {name}", POSE_COLOURS[number % len(POSE_COLOURS)], True)
    safe_scene(destination, groups, "addbinder", palette=palette,
               alignment="All poses share the host frame; each binder is a separate object.")


class _Skip(Exception):
    """A pose with nothing to build (no random placement was found)."""


def run_addbinder(config_path, dry_run=False, output_root=None):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        return _run_addbinder(config_path, dry_run, output_root)


def _run_addbinder(config_path, dry_run=False, output_root=None):
    cfg = load_addbinder_config(config_path, check_paths=not dry_run, output_root=output_root)
    if dry_run:
        print(f"addbinder configuration validated; no inputs opened. {len(cfg.poses)} pose(s) "
              f"+ {cfg.random_poses} random; "
              f"output: {cfg.output_dir}")
        return 0
    destination = Path(cfg.output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    print("=" * 78)
    print("ADD BINDER WORKFLOW")
    print("=" * 78)
    host = inspect_system(cfg.host, gro=cfg.host_gro, overrides=cfg.environment_overrides)
    solute = [m for m in host.molecules if m.category in ("protein", "ligand")]
    environment = [m for m in host.molecules if m.category not in ("protein", "ligand")]
    print(f"[host] {host.gro}: proteins={host.report['selected_proteins']}, "
          f"ligands={host.report['selected_ligands']}")
    frame = membrane_frame(host.universe, host.molecules)
    universe = host.universe
    pieces = []
    for mol in solute:
        atoms = universe.atoms[mol.start:mol.stop]
        whole = _unwrap_sequential(frame.centred(atoms.positions.astype(float)), frame.lengths)
        pieces.append(whole)
    solute_positions = np.concatenate(pieces)
    solute_atoms = universe.atoms[_indices(solute)]
    heavy_mask = np.array([_is_heavy(a) for mol in solute for a in mol.definition.atoms])
    host_heavy = solute_positions[heavy_mask]
    host_heavy_index = solute_atoms.indices[heavy_mask]
    lipid_atoms = _heavy(universe.atoms[_indices([m for m in environment
                                                   if m.category in ("lipid", "sterol", "detergent")])])
    lipid_heavy = frame.centred(lipid_atoms.positions.astype(float)) if len(lipid_atoms) else np.zeros((0, 3))
    pick = np.argmax(host_heavy[:, 2]) if cfg.side == "upper" else np.argmin(host_heavy[:, 2])
    tip_atom = universe.atoms[int(host_heavy_index[pick])]
    tip_mol = next(m for m in solute if m.start <= tip_atom.index < m.stop)
    plane = frame.upper_plane if cfg.side == "upper" else frame.lower_plane
    beyond = (host_heavy[pick, 2] - plane) * (1 if cfg.side == "upper" else -1)
    tip = {"index": int(tip_atom.index), "position": host_heavy[pick].copy()}
    tip_report = {"moleculetype": tip_mol.name, "resname": str(tip_atom.resname), "resid": int(tip_atom.resid),
                  "atom": str(tip_atom.name), "atom_number": int(tip_atom.index) + 1,
                  "z_original_frame": round(float(host_heavy[pick, 2] - frame.shift), 3),
                  "beyond_headgroup_plane_angstrom": round(float(beyond), 3)}
    frame_report = {"upper_plane_z": round(frame.upper_plane - frame.shift, 3),
                    "lower_plane_z": round(frame.lower_plane - frame.shift, 3),
                    "thickness_angstrom": round(frame.upper_plane - frame.lower_plane, 3),
                    "headgroup_atoms": frame.headgroup_atoms,
                    "leaflet_z_spread_angstrom": [round(v, 3) for v in frame.leaflet_spread],
                    "box_angstrom": [round(float(v), 3) for v in frame.lengths]}
    print(f"[frame] headgroup planes {frame_report['upper_plane_z']:.1f} / {frame_report['lower_plane_z']:.1f} A; "
          f"tip {tip_report['resname']}{tip_report['resid']}:{tip_report['atom']} "
          f"{tip_report['beyond_headgroup_plane_angstrom']:.1f} A beyond the {cfg.side} plane")
    if beyond <= 0:
        print("WARNING: the tip does not reach beyond the headgroup plane on this side")

    summary, built = [], []
    with tempfile.TemporaryDirectory(prefix=".addbinder-", dir=destination) as temporary:
        binder = load_binder(cfg, temporary, host)
        binder_report = binder.report
        print(f"[binder] {binder_report['moleculetype']} ({binder_report['category']}), "
              f"{binder_report['atoms']} atoms, charge {binder_report['net_charge']:+.3f} e"
              + (f"; renamed " + ", ".join(f"{old} -> {new}" for new, old in binder_report["renamed_from"].items())
                 if binder_report["renamed_from"] else ""))
        ligand_itps = [d.source_path for d, c in zip(binder.definitions, binder.categories) if c == "ligand"]
        if ligand_itps:
            forcefields = ([str(p) for p in binder.files[:-len(binder.definitions)]]
                           + [str(p) for p in host.files])
            try:
                check_ligand_parameters(ligand_itps, forcefields)
            except TopologyError as exc:
                raise ConfigError(str(exc)) from exc
        heavy_mask = np.array([_is_heavy(a) for a in binder.atoms])
        drawn = random_poses(cfg, frame, host_heavy, lipid_heavy, binder.positions, heavy_mask, tip)
        if drawn:
            print(f"[random] {len(drawn)} pose(s) from random directions, seed {cfg.random_seed}")
        plan = [(pose, [], None) for pose in cfg.poses] + drawn
        for pose, unplaced, draw in plan:
            stage = Path(temporary) / pose.name
            stage.mkdir()
            report = {"mode": "addbinder", "generated_at": datetime.now(timezone.utc).isoformat(),
                      "host": cfg.host, "side": cfg.side, "pose": dict(pose.__dict__),
                      "settings": {"min_image_gap": cfg.min_image_gap, "min_membrane_gap": cfg.min_membrane_gap,
                                   "lateral_reference": cfg.lateral_reference,
                                   "water_clash_distance": cfg.water_clash_distance},
                      "frame": frame_report, "tip": tip_report, "binder": binder_report}
            if draw is not None:
                report["random"] = {"seed": cfg.random_seed, "draw": draw,
                                    "distance_range": cfg.random_distance, "max_tilt": cfg.random_max_tilt,
                                    "min_angle": cfg.random_min_angle,
                                    "min_separation": cfg.random_min_separation}
            try:
                if unplaced:
                    report["status"], report["failures"] = "rejected", unplaced
                    raise _Skip
                result = _build_pose(pose, cfg, host, frame, solute, environment, host_heavy, lipid_heavy,
                                     tip, binder, stage, report)
                if result is not None:
                    built.append((pose.name, result))
            except _Skip:
                pass
            except (ConfigError, TopologyError) as exc:
                report["status"], report["error"] = "failed", str(exc)
                for item in list(stage.iterdir()):
                    shutil.rmtree(item) if item.is_dir() else item.unlink()
            _write_report(stage, report)
            _publish(stage, destination / pose.name)
            checks = report.get("checks", {})
            summary.append({"pose": pose.name, "status": report["status"],
                            "gap_to_host": checks.get("gap_to_host_angstrom"),
                            "gap_to_host_image": checks.get("gap_to_host_image_angstrom"),
                            "gap_to_own_image": checks.get("gap_to_own_image_angstrom"),
                            "gap_to_membrane": checks.get("gap_to_membrane_angstrom"),
                            "binder_height": checks.get("binder_height_along_normal_angstrom"),
                            "waters_ions_removed": len(report.get("solvent_removed", [])),
                            "ions_removed": len(report.get("charge", {}).get("removed_ions", [])),
                            "ions_added": len(report.get("charge", {}).get("added_ions", [])),
                            "reasons": report.get("failures") or ([report["error"]] if "error" in report else [])})
            where = "" if pose.centroid else f" at {pose.distance:g} A"
            line = f"[{pose.name}] {_describe_pose(report['pose'])}{where}: {report['status']}"
            if checks:
                line += (f": host {checks['gap_to_host_angstrom']:.1f} A, host image "
                         f"{checks['gap_to_host_image_angstrom']:.1f} A, own image "
                         f"{checks['gap_to_own_image_angstrom']:.1f} A, membrane "
                         f"{checks['gap_to_membrane_angstrom']:.1f} A")
            print(line)
            if report.get("distance_attempts") and len(report["distance_attempts"]) > 1:
                print("      distance reduced: " + ", ".join(f"{a['distance']:g} A" for a in report["distance_attempts"])
                      + (" (fits)" if report["status"] == "accepted" else " (none fit)"))
            salt = report.get("charge", {})
            if "final_counts" in salt:
                print(f"      salt {salt['final_concentration_molar']:.3f} M (host "
                      f"{salt['host']['concentration_molar']:.3f} M): {len(salt['removed_ions'])} ion(s) removed, "
                      f"{len(salt['added_ions'])} added; net charge {salt['after']:+.0f} e")
            for reason in summary[-1]["reasons"]:
                print(f"      {reason}")
    (destination / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    width = max(14, *(len(s["pose"]) + 2 for s in summary))
    columns = (("gap_to_host", "host", 7), ("gap_to_host_image", "h-image", 9),
               ("gap_to_own_image", "own-img", 9), ("gap_to_membrane", "membr", 8),
               ("binder_height", "height", 8), ("waters_ions_removed", "W/ion", 7),
               ("ions_removed", "ions-", 7), ("ions_added", "ions+", 7))
    rows = [f"{'pose':<{width}}{'status':<10}" + "".join(f"{label:>{w}}" for _, label, w in columns)]
    for s in summary:
        cells = []
        for key, _, w in columns:
            v = s[key]
            cells.append(f"{v:>{w}.1f}" if isinstance(v, float) else f"{'-' if v is None else v:>{w}}")
        rows.append(f"{s['pose']:<{width}}{s['status']:<10}" + "".join(cells))
        rows += [f"    {reason}" for reason in s["reasons"]]
    (destination / "summary.txt").write_text("Distances in Angstrom (heavy atoms)\n" + "\n".join(rows) + "\n")
    if built:
        _write_scene(destination, built)
    accepted = sum(s["status"] == "accepted" for s in summary)
    print(f"[complete] {accepted}/{len(summary)} pose(s) built; see {destination / 'summary.txt'}")
    from run_pipeline import _announce_output_location
    _announce_output_location(destination)
    return 0 if accepted else 3
