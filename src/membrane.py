"""Replace the membrane of a prepared system with an equilibrated bilayer: POPC by default.

The host is one prepared CHARMM-GUI-style folder (topol.top, toppar/, GRO). Its
proteins and ligands never move, so the membrane placement CHARMM-GUI derived
from OPM is kept. Every host lipid and sterol is removed. The host midplane is
found from its lipids (the tail ends of the phospholipids), and an equilibrated
reference bilayer is repeated periodically and cut to the host box with its own
midplane on the host's: by default the bundled 500 ns CHARMM36 POPC bilayer
(Zenodo 10.5281/zenodo.1010083), or any bilayer GRO given as reference_gro
(mixed lipids, cholesterol, asymmetric leaflets).

Each reference residue becomes a host moleculetype: the one of the same name,
else the host lipid or sterol (or a lipid_itps entry) with exactly the same atom
names, so CLOL becomes CHL1 without being told; resname_map overrides this.
Reference lipids that overlap the protein or ligands, or a periodic copy of
another at the box edge, are removed whole. The reference hydration water
replaces host water near the new headgroups; host water and ions that overlap
are removed; counterions and the host salt concentration are restored; a
GROMACS equilibration folder lets the lipids close the gaps around the protein.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from importlib.resources import files
from pathlib import Path
import warnings

import MDAnalysis as mda
from MDAnalysis.lib.distances import capped_distance, self_capped_distance
from MDAnalysis.lib.mdamath import triclinic_vectors
import numpy as np

from addbinder import _unwrap_sequential
from charmm_environment import environment_name
from charmprot import _heavy, _indices, _included_files
from config import ConfigError
from equilibration import EquilibrationSpec, repack_stages
from itp import parse_itp
from system_edit import (write_host_protocol, balance_charge, generated, ion_block, load_section, monovalent_ions, net_charge, order,
                         place, positive, run_edit, salt_summary, write_system)
from topology import TopologyError, check_parameters

LIPID = "POPC"
REMOVED_CATEGORIES = ("lipid", "sterol", "detergent")
LIPID_CATEGORIES = ("lipid", "sterol")
SOLVENT_CATEGORIES = ("water", "ion", "solvent")
WATER_MOLARITY = 55.51   # mol/L of pure water: salt concentration = pairs / waters * 55.51
ION_MIN_SPACING = 5.0    # A between an added ion and any other ion
DEFAULT_REFERENCE = "popc512_charmm36_nacl150mM_310K_500ns.gro.gz"
REFERENCE_SOURCE = ("POPC bilayer, 512 lipids, CHARMM36, 150 mM NaCl, 310 K, frame at 500 ns. "
                    "H. Martinez-Seara, Zenodo, doi:10.5281/zenodo.1010083 (CC BY 4.0); "
                    "lipids only, water and ions removed")
GENERATED = generated("membrane")
# Settings of earlier versions that are no longer used (the anchor rule is gone).
RETIRED = ("anchor_count", "anchor_separation")


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

@dataclass
class MembraneSpec:
    host: str
    host_gro: str | None = None
    # Equilibrated bilayer to take the lipids from; empty: the bundled 500 ns POPC bilayer.
    reference_gro: str | None = None
    lipid_itps: list[str] = field(default_factory=list)   # ITPs of reference lipids the host does not have
    lipid_forcefield: str | None = None   # their parameters, when the host forcefield.itp lacks them
    resname_map: dict = field(default_factory=dict)        # reference residue -> moleculetype (overrides)
    # Midplane z (A, host GRO frame); empty: found from the host lipids.
    membrane_center_z: float | None = None
    max_tilt: float = 5.0                # degrees between a leaflet's P plane and the xy plane
    protein_clash_distance: float = 2.0  # A, lipid heavy atom to protein/ligand heavy atom: lipid removed
    # Lateral registrations of the reference tried; the one keeping the most lipids is used.
    lateral_search: int = 16
    # Reference water hydrates the new bilayer up to this far beyond each P plane (A)
    # and replaces the host water there, except host water near the protein.
    water_layer: float = 10.0
    protein_hydration_distance: float = 4.0  # A; host water/ions this close to the protein stay
    # A; lipids meeting their neighbours across the box edge closer than this are removed (and
    # refilled, below). Minimization resolves contacts down to 1 A (tested).
    reference_contact_distance: float = 1.0
    # The box is never changed. A lipid removed at the box edge is replaced by another
    # equilibrated lipid of the same leaflet at the same head position, turned until
    # its heavy atoms stay this far (A) from every lipid; minimization and repack/
    # relax what is left.
    refill_edge: bool = True
    edge_refill_distance: float = 1.0
    # Final clash check, every atom (hydrogens included) of different molecules: pairs
    # closer than clash_distance (A) are reported; with remove_clashes the new lipid or
    # the water/ion of each pair is removed (never the protein or ligands).
    clash_distance: float = 0.2
    remove_clashes: bool = True
    water_clash_distance: float = 2.0    # A, host water/ion heavy atom to lipid or reference water: removed
    neutralize: bool = True
    target_net_charge: float = 0.0
    # Bulk waters outside the new hydration layer become salt pairs until the host's
    # salt concentration is reached again (the reference water brings no ions).
    restore_salt: bool = True
    ion_exclusion_distance: float = 10.0
    lipid_ion_exclusion_distance: float = 5.0
    environment_overrides: dict = field(default_factory=dict)
    output_dir: str = "membrane_output"
    equilibration: EquilibrationSpec = field(default_factory=EquilibrationSpec)
    # A second set: relaxation of the new membrane with the protein and ligands held
    # (like the cholesterol mode's repack/); written next to equilibration/.
    # A third set: the host folder's own mdp files, scripts and README on the new
    # system (written to host_protocol/ when the host folder has mdp files).
    copy_host_protocol: bool = True
    repack: EquilibrationSpec = field(
        default_factory=lambda: EquilibrationSpec(output_dir="repack", stages=repack_stages()))


def load_membrane_config(path, check_paths=True, output_root=None):
    import yaml
    from config import _UniqueKeySafeLoader
    retired, legacy = [], False
    try:
        raw = yaml.load(Path(path).read_text(), Loader=_UniqueKeySafeLoader)
        if isinstance(raw, dict) and "membrane" not in raw and isinstance(raw.get("popc"), dict):
            raw["membrane"] = raw.pop("popc")  # files written for the former name of this mode
            legacy = True
        if isinstance(raw, dict) and isinstance(raw.get("membrane"), dict):
            retired = [k for k in RETIRED if k in raw["membrane"]]
    except (OSError, yaml.YAMLError):
        pass
    optional = ("host_gro", "reference_gro", "lipid_itps", "lipid_forcefield")
    if retired or legacy:
        # Older files: rewrite them into a temporary copy (next to the original, so
        # relative paths still resolve) before the strict check.
        import tempfile
        raw["membrane"] = {k: v for k, v in raw["membrane"].items() if k not in RETIRED}
        handle = tempfile.NamedTemporaryFile("w", suffix=".yaml", dir=Path(path).parent, delete=False)
        with handle:
            yaml.safe_dump(raw, handle)
        try:
            cfg = load_section(handle.name, "membrane", MembraneSpec, check_paths, output_root,
                               optional_paths=optional)
        finally:
            Path(handle.name).unlink()
    else:
        cfg = load_section(path, "membrane", MembraneSpec, check_paths, output_root, optional_paths=optional)
    positive(cfg, "membrane", "max_tilt", "protein_clash_distance", "reference_contact_distance", "edge_refill_distance", "clash_distance",
             "protein_hydration_distance", "water_clash_distance", "ion_exclusion_distance",
             "lipid_ion_exclusion_distance")
    positive(cfg, "membrane", "water_layer", allow_zero=True)
    if cfg.max_tilt >= 45:
        raise ConfigError("membrane.max_tilt must be below 45 degrees")
    if not 1 <= cfg.lateral_search <= 400:
        raise ConfigError("membrane.lateral_search must be a whole number from 1 to 400")
    if any(not isinstance(k, str) or not isinstance(v, str) for k, v in cfg.resname_map.items()):
        raise ConfigError("membrane.resname_map must map residue names to moleculetype names")
    cfg.retired_settings = retired
    return cfg


load_config = load_membrane_config


# ---------------------------------------------------------------------------
# Geometry (membrane's tested helpers)
# ---------------------------------------------------------------------------

def default_reference_path():
    return Path(str(files("gmxtransplant").joinpath("data", DEFAULT_REFERENCE)))


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


def _wrapped_dz(z, midplane, height):
    dz = z - midplane
    return dz - height * np.round(dz / height)


def _lipid_heavy(lipids, definition):
    """Heavy-atom positions of placed lipids: (k, side, positions) items of one moleculetype,
    or (k, side, positions, heavy_mask) items of several."""
    if not lipids:
        return np.zeros((0, 3))
    heavy = _heavy_mask(definition) if definition is not None else None
    return np.concatenate([item[2][item[3] if len(item) > 3 else heavy] for item in lipids])


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


def place_water(cfg, host, frame, reference, lipids, definition, protected, origin):
    """Reference water on the new bilayer and the host solvent kept around it.

    Returns (reference waters [(n, 3, 3) positions], kept host solvent, removed host solvent by reason, report).
    """
    universe = host.universe
    box = universe.dimensions
    host_cell, height = _cell(box)
    lipid_heavy = _lipid_heavy(lipids, definition)
    protected_heavy = _heavy(universe.atoms[_indices(protected)]).positions
    solvent = [m for m in host.molecules if m.category in SOLVENT_CATEGORIES]
    solvent_heavy = [_heavy(universe.atoms[m.start:m.stop]).positions.astype(float) for m in solvent]
    removed = {"replaced_by_reference_water": [], "overlap_with_new_lipids": [], "overlap_with_reference_water": []}

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
        removed["overlap_with_new_lipids"].append(remaining[j])
    gone |= set(removed["overlap_with_new_lipids"])

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
        drop |= _touching(oxygens, lipid_heavy, cfg.water_clash_distance, box)
        survivors = [i for i in range(len(placed)) if i not in drop]
        edge = _resolve_contacts([oxygens[i] for i in survivors], cfg.water_clash_distance, box)
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

# ---------------------------------------------------------------------------
# Host midplane and reference bilayer
# ---------------------------------------------------------------------------

def _head_index(definition):
    """The headgroup atom: P for phospholipids, else the first oxygen (sterol O3), else atom 0."""
    names = [a.atom_name for a in definition.atoms]
    if "P" in names:
        return names.index("P")
    return next((i for i, a in enumerate(definition.atoms) if a.atom_name.upper().startswith("O")), 0)


@dataclass
class HostMembrane:
    shift: float          # added to z to put the midplane at height / 2
    height: float
    midplane: float       # GRO frame, A
    planes: dict          # leaflet -> report (phosphorus planes), or {} without phospholipids
    counts: dict          # leaflet -> Counter of host lipid moleculetypes
    source: str


def measure_host(host, cfg):
    """The host bilayer's midplane (tail ends of the phospholipids) and its leaflets."""
    universe = host.universe
    _, height = _cell(universe.dimensions)
    lengths = universe.dimensions[:3].astype(float)
    membrane = [m for m in host.molecules if m.category in REMOVED_CATEGORIES]
    if not membrane:
        raise ConfigError("The host has no lipids or sterols: there is no membrane to replace")
    heads, tails, others = [], [], []
    for mol in membrane:
        names = [a.atom_name for a in mol.definition.atoms]
        whole = _unwrap_sequential(universe.atoms[mol.start:mol.stop].positions.astype(float), lengths)
        if "P" in names:
            p = names.index("P")
            heads.append((mol, whole[p]))
            tails.append(_tail_end(whole, p))
        else:
            others.append((mol, whole))
    if len(heads) >= 6:
        shift = float(_midplane_shift(np.array(tails), height))
        source = "tail ends of the phospholipids"
    else:
        # No phospholipids to speak of: the circular mean of every membrane heavy atom.
        z = np.concatenate([_heavy(universe.atoms[m.start:m.stop]).positions[:, 2] for m in membrane])
        shift = float(_midplane_shift(z.astype(float), height))
        source = "membrane heavy atoms (fewer than six phospholipids)"
    planes, counts = {}, {"upper": Counter(), "lower": Counter()}
    if heads:
        positions = np.array([p for _, p in heads])
        zc = np.mod(positions[:, 2] + shift, height)
        for side, on in (("upper", zc > height / 2), ("lower", zc <= height / 2)):
            counts[side].update(heads[i][0].name for i in np.flatnonzero(on))
            if on.sum() >= 3:
                _, tilt, spread = _plane(positions[on, :2], zc[on])
                if tilt > cfg.max_tilt:
                    raise ConfigError(f"The {side} leaflet's P plane is tilted {tilt:.1f} degrees from the xy plane "
                                      f"(max_tilt {cfg.max_tilt:g}); the membrane must be flat and normal to z")
                planes[side] = {"mean_z": round(float(np.mean(zc[on]) - shift), 3), "tilt_degrees": round(tilt, 3),
                                "spread_angstrom": round(spread, 3), "phosphorus_atoms": int(on.sum())}
    for mol, whole in others:
        centre = np.mod(whole[:, 2].mean() + shift, height)
        counts["upper" if centre > height / 2 else "lower"][mol.name] += 1
    if len(planes) == 2:
        thickness = planes["upper"]["mean_z"] - planes["lower"]["mean_z"]
        if not 20.0 <= thickness <= 70.0:
            raise ConfigError(f"Headgroup planes are {thickness:.1f} A apart; expected a flat bilayer normal to z")
    midplane = float(np.mod(height / 2 - shift, height))
    if cfg.membrane_center_z is not None:
        midplane, source = float(cfg.membrane_center_z), "membrane_center_z"
    return HostMembrane(shift, height, midplane, planes, counts, source)


def _definitions(host, cfg):
    """Moleculetypes the reference may use: every host definition plus lipid_itps."""
    found = dict(host.definitions)
    extra = {}
    for path in cfg.lipid_itps:
        for name, definition in parse_itp(path).items():
            if name in found and Path(found[name].source_path).read_bytes() != Path(path).read_bytes():
                raise ConfigError(f"lipid_itps defines {name}, which the host defines differently; rename one")
            extra[name] = definition
    return {**found, **extra}, extra


def _match(resname, names, definitions, membrane_names, cfg):
    """The moleculetype for a reference residue: resname_map, same name, or a membrane
    moleculetype (host lipid or sterol, or lipid_itps) with exactly the same atom names."""
    if resname in cfg.resname_map:
        return cfg.resname_map[resname], "resname_map"
    if resname in definitions:
        return resname, "same name"
    candidates = [name for name in membrane_names
                  if sorted(a.atom_name for a in definitions[name].atoms) == sorted(names)]
    if len(candidates) == 1:
        return candidates[0], "same atom names"
    return None, None


@dataclass
class ReferenceBilayer:
    path: str
    lipids: list      # per lipid: (definition, whole positions, midplane at z = 0, ITP order)
    leaflets: list
    heads: np.ndarray
    waters: np.ndarray
    cell: np.ndarray
    height: float
    report: dict


def load_reference(path, definitions, cfg, membrane_names=()):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        universe = mda.Universe(str(path))
    if universe.dimensions is None or np.any(np.abs(universe.dimensions[3:] - 90) > 1e-3):
        raise ConfigError(f"The reference {path} needs a rectangular box")
    cell, height = _cell(universe.dimensions)
    lengths = universe.dimensions[:3].astype(float)
    lipids, ignored, reordered, mapped = [], Counter(), 0, {}
    for residue in universe.residues:
        resname = str(residue.resname)
        names = [str(n) for n in residue.atoms.names]
        known = environment_name(resname)
        if known and known[0] in ("water", "ion", "solvent"):
            ignored[resname] += 1      # water is read separately; ions and solvent are not used
            continue
        name, how = _match(resname, names, definitions, membrane_names, cfg)
        if name is None:
            if known and known[0] in LIPID_CATEGORIES or len(names) > 20:
                raise ConfigError(f"Reference residue {resname} matches no moleculetype: add its ITP to lipid_itps "
                                  "or map it with resname_map")
            ignored[resname] += 1
            continue
        definition = definitions[name]
        if how != "same name":
            mapped[resname] = f"{name} ({how})"
        expected = [a.atom_name for a in definition.atoms]
        order_ = slice(None)
        if names != expected:
            if sorted(names) != sorted(expected) or len(set(names)) != len(names):
                raise ConfigError(f"Reference {resname} {residue.resid} atoms do not match the {name} ITP")
            order_ = [names.index(n) for n in expected]
            reordered += 1
        lipids.append((definition, _unwrap_sequential(residue.atoms.positions.astype(float)[order_], lengths)))
    if len(lipids) < 20:
        raise ConfigError(f"The reference {path} holds {len(lipids)} lipids; expected a bilayer")
    phospholipids = [(d, x) for d, x in lipids if "P" in [a.atom_name for a in d.atoms]]
    if len(phospholipids) < 6:
        raise ConfigError("The reference has fewer than six phospholipids; cannot find its midplane")
    shift = _midplane_shift(np.array([_tail_end(x, _head_index(d)) for d, x in phospholipids]), height)
    inverse = np.linalg.inv(cell)
    out, leaflets = [], []
    for definition, whole in lipids:
        h = _head_index(definition)
        whole = whole.copy()
        whole[:, 2] += shift - height * np.floor((whole[h, 2] + shift) / height) - height / 2
        whole[:, :2] -= np.floor(whole[h, :2] @ inverse) @ cell
        out.append((definition, whole))
        leaflets.append("upper" if whole[h, 2] > 0 else "lower")
    heads = np.array([x[_head_index(d)] for d, x in out])
    planes = {}
    for side in ("upper", "lower"):
        on = [i for i, (d, _) in enumerate(out) if leaflets[i] == side and "P" in [a.atom_name for a in d.atoms]]
        _, tilt, spread = _plane(heads[on, :2], heads[on, 2])
        if tilt > 5.0:
            raise ConfigError(f"The reference {side} leaflet is tilted {tilt:.1f} degrees; expected a flat bilayer")
        composition = Counter(d.name for i, (d, _) in enumerate(out) if leaflets[i] == side)
        planes[side] = {"lipids": int(sum(composition.values())), "composition": dict(composition),
                        "mean_z": round(float(heads[on, 2].mean()), 3), "spread_angstrom": round(spread, 3)}
    waters = _reference_waters(universe, lengths, shift, cell, height)
    area = abs(np.linalg.det(cell))
    report = {"path": str(path), "lipids": len(out), "waters": len(waters),
              "box_angstrom": [round(float(v), 3) for v in lengths],
              "phosphate_thickness_angstrom": round(planes["upper"]["mean_z"] - planes["lower"]["mean_z"], 3),
              "area_per_lipid_angstrom2": {s: round(area / planes[s]["lipids"], 2) for s in planes},
              "leaflets": planes, "residues_reordered_to_itp": reordered, "residues_matched": mapped,
              "ignored_residues": dict(ignored)}
    return ReferenceBilayer(str(path), out, leaflets, heads, waters, cell, height, report)


# ---------------------------------------------------------------------------
# Placement
# ---------------------------------------------------------------------------

def _lipids_at(cfg, reference, host_cell, frame, origin, protected_heavy, box):
    candidates = []
    for k, offset in tile(reference.heads[:, :2], reference.cell, host_cell, origin):
        definition, whole = reference.lipids[k]
        coords = whole.copy()
        coords[:, :2] += offset
        coords[:, 2] += frame.midplane
        candidates.append((k, reference.leaflets[k], coords, _heavy_mask(definition)))
    clash = _touching([c[2][c[3]] for c in candidates], protected_heavy, cfg.protein_clash_distance, box)
    keep = [c for i, c in enumerate(candidates) if i not in clash]
    edge = _resolve_contacts([c[2][c[3]] for c in keep], cfg.reference_contact_distance, box)
    kept = [c for i, c in enumerate(keep) if i not in edge]

    def by_type(items):
        return dict(Counter(f"{side} {reference.lipids[k][0].name}" for k, side, _, _ in items))
    report = {"tiled": by_type(candidates), "removed_protein_overlap": by_type([candidates[i] for i in clash]),
              "removed_box_edge_overlap": by_type([keep[i] for i in edge])}
    return kept, report, [keep[i] for i in edge]


REFILL_ROTATIONS = 24    # orientations about z tried per reference lipid
REFILL_CANDIDATES = 256  # reference lipids (of the same leaflet) tried per position


def refill_edge(cfg, reference, kept, dropped, protected_all, protected_heavy, box, host_cell, seed=7):
    """Put a lipid back at each head position emptied at the box edge.

    A lipid removed where the reference was cut to the host box leaves a hole of one
    lipid. Other equilibrated lipids of the same leaflet are tried there, head on the
    same spot, turned about z. The first is used whose heavy atoms stay
    edge_refill_distance from the other lipids' and protein_clash_distance from the
    protein's, and whose every atom (hydrogens included: a refilled lipid sits in a
    tight spot) stays clash_distance from every atom of the lipids, protein, ligands
    and the host water kept near them. Minimization and repack/ relax the rest.
    """
    rng = np.random.default_rng(seed)
    added = []
    for k_old, side, coords_old, _ in dropped:
        head = coords_old[_head_index(reference.lipids[k_old][0])]
        current = np.concatenate([item[2] for item in kept + added])
        current_heavy = np.concatenate([item[2][item[3]] for item in kept + added])

        def nearby(points):
            return points[np.linalg.norm(_min_image_xy(points[:, :2] - head[:2], host_cell), axis=1) < 30.0]
        near_lipids, near_lipids_heavy = nearby(current), nearby(current_heavy)
        near_protein, near_protein_heavy = nearby(protected_all), nearby(protected_heavy)
        pool = [k for k in range(len(reference.lipids)) if reference.leaflets[k] == side]
        found = None
        for k in rng.choice(pool, size=min(REFILL_CANDIDATES, len(pool)), replace=False):
            definition, whole = reference.lipids[k]
            mask = _heavy_mask(definition)
            pivot = whole[_head_index(definition)]
            for angle in np.linspace(0, 2 * np.pi, REFILL_ROTATIONS, endpoint=False):
                c, s = np.cos(angle), np.sin(angle)
                coords = whole.copy()
                coords[:, :2] = (coords[:, :2] - pivot[:2]) @ np.array([[c, s], [-s, c]]) + head[:2]
                coords[:, 2] += head[2] - pivot[2]
                if len(near_lipids_heavy) and len(capped_distance(coords[mask], near_lipids_heavy,
                                                                  cfg.edge_refill_distance, box=box,
                                                                  return_distances=False)):
                    continue
                if len(near_lipids) and len(capped_distance(coords, near_lipids, cfg.clash_distance, box=box,
                                                            return_distances=False)):
                    continue
                if len(near_protein) and len(capped_distance(coords, near_protein, cfg.clash_distance,
                                                             box=box, return_distances=False)):
                    continue
                if len(near_protein_heavy) and len(capped_distance(coords[mask], near_protein_heavy,
                                                                   cfg.protein_clash_distance, box=box,
                                                                   return_distances=False)):
                    continue
                found = (int(k), side, coords, mask)
                break
            if found:
                break
        if found:
            added.append(found)
    return added


def place_lipids(cfg, host, frame, reference, protected):
    """Reference lipids on the host midplane, cut to the host box; the best of lateral_search registrations."""
    box = host.universe.dimensions
    host_cell, _ = _cell(box)
    protected_heavy = _heavy(host.universe.atoms[_indices(protected)]).positions
    side = int(np.ceil(np.sqrt(cfg.lateral_search)))
    origins = [(i / side) * reference.cell[0] + (j / side) * reference.cell[1]
               for i in range(side) for j in range(side)][:cfg.lateral_search]
    best = None
    for number, origin in enumerate(origins):
        kept, report, dropped = _lipids_at(cfg, reference, host_cell, frame, origin, protected_heavy, box)
        sides = Counter(item[1] for item in kept + dropped)
        score = (abs(sides["upper"] - sides["lower"]) > 2, -len(kept) - len(dropped), len(dropped))
        if best is None or score < best[4]:
            best = (kept, report, origin, number, score, dropped)
    kept, report, origin, number, _, dropped = best
    report["refilled_box_edge"] = {}
    if cfg.refill_edge and dropped:
        # Everything a refilled lipid must not touch: protein and ligands, and the host
        # water and ions kept near them (protein_hydration_distance).
        universe = host.universe
        solvent = [m for m in host.molecules if m.category in SOLVENT_CATEGORIES]
        near = _touching([_heavy(universe.atoms[m.start:m.stop]).positions for m in solvent], protected_heavy,
                         cfg.protein_hydration_distance, box)
        guarded = [universe.atoms[_indices(protected)].positions] + \
                  [universe.atoms[solvent[i].start:solvent[i].stop].positions for i in sorted(near)]
        added = refill_edge(cfg, reference, kept, dropped, np.concatenate(guarded), protected_heavy, box,
                            host_cell)
        kept = kept + added
        report["refilled_box_edge"] = dict(Counter(f"{s} {reference.lipids[k][0].name}" for k, s, _, _ in added))
    kept.sort(key=lambda item: (reference.lipids[item[0]][0].name, item[1] != "upper", item[0]))
    report["placed"] = dict(Counter(f"{s} {reference.lipids[k][0].name}" for k, s, _, _ in kept))
    report["lateral_search"] = {"registrations_tried": len(origins), "chosen": number,
                                "origin_angstrom": [round(float(v), 3) for v in origin]}
    return kept, origin, report


def resolve_clashes(cfg, host, protected, lipids, waters, kept_solvent, reference):
    """Every atom, hydrogens included, of the new lipids, the reference water and the
    kept host solvent against each other and the protein/ligands: pairs closer than
    clash_distance. With remove_clashes, the water or ion of a pair goes first, then the
    lipid; the protein and ligands never move. Returns the kept items and a report."""
    universe, box = host.universe, host.universe.dimensions
    blocks = [("fixed", None, universe.atoms[_indices(protected)].positions)] if protected else []
    blocks += [("lipid", i, item[2]) for i, item in enumerate(lipids)]
    blocks += [("water", i, w) for i, w in enumerate(waters)]
    blocks += [("solvent", i, universe.atoms[m.start:m.stop].positions) for i, m in enumerate(kept_solvent)]
    positions = np.concatenate([b[2] for b in blocks])
    owner = np.concatenate([np.full(len(b[2]), n) for n, b in enumerate(blocks)])
    # A small margin: the written GRO rounds coordinates to 0.01 A, which can move a pair
    # just under the cutoff that the final check of the written system applies.
    pairs, distances = capped_distance(positions, positions, cfg.clash_distance + 0.02, box=box)
    across = owner[pairs[:, 0]] < owner[pairs[:, 1]]
    pairs, distances = pairs[across], distances[across]

    def label(n, atom):
        kind, i, _ = blocks[n]
        if kind == "fixed":
            a = universe.atoms[_indices(protected)][atom]
            return f"{a.resname}{a.resid}:{a.name}"
        if kind == "lipid":
            definition = reference.lipids[lipids[i][0]][0]
            return f"new {definition.name}:{definition.atoms[atom].atom_name}"
        if kind == "water":
            return f"reference water:{('OH2', 'H1', 'H2')[atom]}"
        a = universe.atoms[kept_solvent[i].start + atom]
        return f"host {kept_solvent[i].name}{a.resid}:{a.name}"
    starts = np.concatenate([[0], np.cumsum([len(b[2]) for b in blocks])])
    examples = [f"{label(owner[a], a - starts[owner[a]])} - {label(owner[b], b - starts[owner[b]])} {d:.2f} A"
                for (a, b), d in sorted(zip(pairs.tolist(), distances), key=lambda x: x[1])[:20]]
    report = {"distance_angstrom": cfg.clash_distance, "pairs": int(len(pairs)), "examples": examples,
              "removed": {}}
    if not len(pairs) or not cfg.remove_clashes:
        return lipids, waters, kept_solvent, report
    priority = {"water": 0, "solvent": 0, "lipid": 1, "fixed": 9}
    gone = set()
    for a, b in sorted({(int(owner[a]), int(owner[b])) for a, b in pairs}):
        if a in gone or b in gone:
            continue
        drop = min((a, b), key=lambda n: (priority[blocks[n][0]], -n))
        if blocks[drop][0] == "fixed":
            continue
        gone.add(drop)
    removed = Counter()
    for n in gone:
        kind, i, _ = blocks[n]
        removed[f"new {reference.lipids[lipids[i][0]][0].name}" if kind == "lipid" else
                "reference water" if kind == "water" else f"host {kept_solvent[i].name}"] += 1
    report["removed"] = dict(removed)
    lipids = [x for i, x in enumerate(lipids) if ("lipid", i) not in {blocks[n][:2] for n in gone}]
    keep_w = [i for i in range(len(waters)) if ("water", i) not in {blocks[n][:2] for n in gone}]
    waters = waters[keep_w] if len(waters) else waters
    kept_solvent = [m for i, m in enumerate(kept_solvent) if ("solvent", i) not in {blocks[n][:2] for n in gone}]
    return lipids, waters, kept_solvent, report


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

def run_membrane(config_path, dry_run=False, output_root=None):
    cfg = load_membrane_config(config_path, check_paths=not dry_run, output_root=output_root)
    reference_path = cfg.reference_gro or str(default_reference_path())

    def build(stage, host, report):
        work = stage.parent / "work"
        work.mkdir()
        for key in cfg.retired_settings:
            report["warnings"].append(f"membrane.{key} is no longer used and was ignored (the midplane is found from "
                                      "the host lipids)")
        frame = measure_host(host, cfg)
        removed = [m for m in host.molecules if m.category in REMOVED_CATEGORIES]
        thickness = (round(frame.planes["upper"]["mean_z"] - frame.planes["lower"]["mean_z"], 3)
                     if len(frame.planes) == 2 else None)
        report["host_membrane"] = {"midplane_z": round(frame.midplane, 3), "midplane_source": frame.source,
                                   "phosphate_thickness_angstrom": thickness, "planes": frame.planes,
                                   "lipids_removed": {s: dict(frame.counts[s]) for s in ("upper", "lower")},
                                   "molecules_removed": dict(Counter(m.name for m in removed))}
        print(f"[host] midplane z {frame.midplane:.2f} A; removing "
              + ", ".join(f"{n} {k}" for k, n in Counter(m.name for m in removed).items()))
        definitions, extra = _definitions(host, cfg)
        if extra:
            forcefields = [str(p) for p in host.files]
            if cfg.lipid_forcefield:
                forcefields += [str(p) for p in _included_files(Path(cfg.lipid_forcefield).resolve(),
                                                                Path(cfg.lipid_forcefield).resolve().parent)]
            try:
                check_parameters([d.source_path for d in extra.values()], forcefields,
                                 "Set lipid_forcefield to a file with the missing lipid parameters.")
            except TopologyError as exc:
                raise ConfigError(str(exc)) from exc
        membrane_names = sorted({m.name for m in removed} | set(extra))
        reference = load_reference(reference_path, definitions, cfg, membrane_names)
        report["reference"] = dict(reference.report)
        if not cfg.reference_gro:
            report["reference"]["source"] = REFERENCE_SOURCE
        print(f"[reference] {reference.report['lipids']} lipids: "
              + "; ".join(f"{side} " + ", ".join(f"{n} {k}" for k, n in reference.report["leaflets"][side]["composition"].items())
                          for side in ("upper", "lower"))
              + f"; P-P {reference.report['phosphate_thickness_angstrom']:.1f} A"
              + (f" (host {thickness:.1f} A)" if thickness else ""))
        for resname, how in reference.report["residues_matched"].items():
            print(f"[reference] {resname} -> {how}")
        difference = reference.report["phosphate_thickness_angstrom"] - thickness if thickness else 0.0
        if abs(difference) > 3.0:
            report["warnings"].append(
                f"The reference bilayer is {abs(difference):.1f} A {'thicker' if difference > 0 else 'thinner'} "
                "than the host's; the water layers adjust during the NPT equilibration")
        protected = [m for m in host.molecules if m.category not in REMOVED_CATEGORIES + SOLVENT_CATEGORIES]
        lipids, origin, report["placement"] = place_lipids(cfg, host, frame, reference, protected)
        edge = sum(report["placement"]["removed_box_edge_overlap"].values())
        refilled = sum(report["placement"]["refilled_box_edge"].values())
        if refilled:
            print(f"[membrane] box edge: {edge} overlapping lipids removed, {refilled} put back as other "
                  "equilibrated lipids")
        if edge - refilled > 0.02 * max(len(lipids), 1):
            report["warnings"].append(
                f"{edge - refilled} lipids removed where the reference was cut to the host box could not be put back, so the new "
                "membrane holds fewer lipids than fit; `gmxtransplant membrane-patch` cuts the reference to this "
                "box and re-equilibrates it, after which nothing is lost (see the membrane guide)")
        print("[membrane] placed " + ", ".join(f"{k} {n}" for k, n in report["placement"]["placed"].items())
              + f" (best of {report['placement']['lateral_search']['registrations_tried']} lateral registrations)")
        waters, water_definition, kept_solvent, removed_solvent, report["water"] = place_water(
            cfg, host, frame, reference, lipids, None, protected, origin)
        lipids, waters, kept_solvent, report["clashes"] = resolve_clashes(cfg, host, protected, lipids, waters,
                                                                           kept_solvent, reference)
        if report["clashes"]["pairs"]:
            text = (f"{report['clashes']['pairs']} atom pair(s) of different molecules closer than "
                    f"{cfg.clash_distance:g} A")
            if cfg.remove_clashes:
                print(f"[clashes] {text}: removed " + ", ".join(f"{k} {n}" for k, n in
                                                               report["clashes"]["removed"].items()))
                report["placement"]["placed"] = dict(Counter(f"{s} {reference.lipids[k][0].name}"
                                                             for k, s, _, _ in lipids))
                print("[membrane] final " + ", ".join(f"{k} {n}" for k, n in report["placement"]["placed"].items()))
            else:
                report["warnings"].append(f"{text}; they stay (remove_clashes: false): " +
                                          "; ".join(report["clashes"]["examples"][:5]))
        report["solvent_removed"] = {reason: dict(Counter(m.name for m in items))
                                     for reason, items in removed_solvent.items()}
        if report["water"].get("reference_waters"):
            print(f"[water] {report['water']['reference_waters']} reference waters hydrate the bilayer")
        elif report["water"].get("note"):
            report["warnings"].append(report["water"]["note"])
        box = host.universe.dimensions
        lipid_heavy = np.concatenate([c[2][c[3]] for c in lipids]) if lipids else np.zeros((0, 3))
        protected_heavy = _heavy(host.universe.atoms[_indices(protected)]).positions
        net = net_charge([(host, m) for m in protected + kept_solvent]) \
            + sum(reference.lipids[k][0].charge for k, _, _, _ in lipids)
        added = []
        report["charge"] = {"before": round(net, 6), "after": round(net, 6), "removed_ions": [], "added_ions": {}}
        if cfg.neutralize:
            cation, anion = monovalent_ions(host)
            if cation is None or anion is None:
                raise ConfigError("The host has no monovalent cation and anion to neutralize with")
            avoid = [(protected_heavy, cfg.ion_exclusion_distance), (lipid_heavy, cfg.lipid_ion_exclusion_distance)]
            kept_solvent, added, report["charge"] = balance_charge(host, kept_solvent, net, cfg.target_net_charge,
                                                                   cation, anion, avoid, box)
            print(f"[charge] {report['charge']['before']:+.3f} e -> {report['charge']['after']:+.3f} e")
        half = report["water"]["slab_half_thickness_angstrom"]
        kept_solvent, salt_ions, salt_report = restore_salt(cfg, host, kept_solvent, len(waters), frame, half,
                                                            protected)
        if salt_report.get("note"):
            report["warnings"].append(salt_report["note"])
        added += salt_ions

        def category(definition):
            known = environment_name(definition.name)
            return known[0] if known and known[0] in LIPID_CATEGORIES else "lipid"
        items = [(reference.lipids[k][0], category(reference.lipids[k][0]), coords) for k, _, coords, _ in lipids]
        lipid_block = place(items, work, label="lipids")
        leaflet_of = {id(m): (side, m.name) for m, (_, side, _, _) in zip(lipid_block.molecules, lipids)}
        new = [(lipid_block, m) for m in lipid_block.molecules]
        if len(waters):
            water_block = place([(water_definition, "water", w) for w in waters], work, label="water")
            new += [(water_block, m) for m in water_block.molecules]
        new += ion_block(added, work)
        ordered = order(host, protected + kept_solvent, new=new)
        report["salt"] = {"host": salt_summary([(host, m) for m in host.molecules]), "final": salt_summary(ordered),
                          "restore": salt_report}
        print(f"[salt] host {report['salt']['host']['concentration_molar']} M -> final "
              f"{report['salt']['final']['concentration_molar']} M")
        groups = {"Upper_leaflet": set(), "Lower_leaflet": set()}
        for key, (side, name) in leaflet_of.items():
            groups["Upper_leaflet" if side == "upper" else "Lower_leaflet"].add(key)
            groups.setdefault(f"{name}_{side}", set()).add(key)
        parameter_files = []
        if cfg.lipid_forcefield:
            parameter_files = list(_included_files(Path(cfg.lipid_forcefield).resolve(),
                                                   Path(cfg.lipid_forcefield).resolve().parent))
        report.update(write_system(
            stage, host, ordered, "membrane", parameter_files=parameter_files, equilibration=cfg.equilibration,
            purpose="Equilibration of the new bilayer", extra_index=groups,
            more_equilibrations=[("repack", cfg.repack, "Repacking of the new membrane around the held protein "
                                  "and ligands (lipids unrestrained)")],
            scene=({"inserted_membrane": set(leaflet_of)},
                   {"inserted_membrane": ("Equilibrated bilayer", [.95, .68, .12], True)},
                   "Protein and ligands keep their input coordinates; the bilayer is new."),
            restrained=[m.definition for m in protected]))
        if cfg.copy_host_protocol:
            protocol = write_host_protocol(stage, host)
            if protocol:
                report["host_protocol"] = protocol
                print(f"[host_protocol] host_protocol/: the host's {len(protocol['mdp_files'])} mdp file(s) and "
                      f"scripts on the new system ({protocol['coordinates']})")
                for name, gone in protocol["dropped_defines"].items():
                    print(f"      {name}: dropped {' '.join(gone)} (not in the new topology)")

    def lines(report):
        out = []
        host = report.get("host_membrane")
        if host:
            out += ["", "HOST MEMBRANE (protein and ligands are not moved: OPM placement kept)",
                    f"  midplane z {host['midplane_z']:.2f} A, from the {host['midplane_source']}"
                    + (f"; P-P thickness {host['phosphate_thickness_angstrom']:.2f} A"
                       if host["phosphate_thickness_angstrom"] else ""),
                    "  removed: " + ", ".join(f"{k} {n}" for k, n in host["molecules_removed"].items())]
        ref = report.get("reference")
        if ref:
            out += ["", "REFERENCE BILAYER", f"  {ref['path']}", f"  {ref.get('source', 'user supplied')}",
                    f"  {ref['lipids']} lipids, box {ref['box_angstrom'][0]:.2f} x {ref['box_angstrom'][1]:.2f} A, "
                    f"P-P thickness {ref['phosphate_thickness_angstrom']:.2f} A, area per lipid "
                    f"{ref['area_per_lipid_angstrom2']['upper']:.1f} / {ref['area_per_lipid_angstrom2']['lower']:.1f} A^2"]
            for side in ("upper", "lower"):
                out.append(f"  {side}: " + ", ".join(f"{k} {n}" for k, n in ref["leaflets"][side]["composition"].items()))
            for resname, how in ref.get("residues_matched", {}).items():
                out.append(f"  {resname} -> {how}")
        placement = report.get("placement")
        if placement:
            out += ["", "PLACED: " + ", ".join(f"{k} {n}" for k, n in placement["placed"].items()),
                    "  removed for protein overlap: " + (", ".join(
                        f"{k} {n}" for k, n in placement["removed_protein_overlap"].items()) or "none"),
                    "  removed at the box edge: " + (", ".join(
                        f"{k} {n}" for k, n in placement["removed_box_edge_overlap"].items()) or "none")]
        water = report.get("water")
        if water and water.get("reference_waters"):
            out += [f"WATER: {water['reference_waters']} reference waters hydrate the new bilayer"]
        for reason, counts in report.get("solvent_removed", {}).items():
            out.append(f"  host solvent removed, {reason.replace('_', ' ')}: "
                       + (", ".join(f"{k} {v}" for k, v in counts.items()) or "none"))
        return out

    return run_edit("membrane", cfg, "REPLACE THE MEMBRANE WITH AN EQUILIBRATED BILAYER", build,
                    "replace the membrane with an equilibrated bilayer", lines, dry_run=dry_run,
                    dry_run_note=f"Reference: {reference_path}; ")


run = run_membrane
