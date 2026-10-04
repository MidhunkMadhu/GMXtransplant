"""Validation for periodic boxes used by distance and concentration code."""

from __future__ import annotations

from dataclasses import dataclass, field
import itertools
from pathlib import Path
import sys
import time
from typing import List

import numpy as np
from MDAnalysis.exceptions import NoDataError
from scipy.optimize import minimize
from scipy.spatial import ConvexHull, QhullError, cKDTree
from scipy.spatial.transform import Rotation


class BoxValidationError(ValueError):
    """Preflight refusal, optionally carrying the complete JSON-ready diagnostic."""

    def __init__(self, message, result=None):
        super().__init__(message)
        self.result = result


@dataclass
class BoxValidationResult:
    dimensions: np.ndarray
    volume: float
    coordinate_span: np.ndarray
    cell_cartesian_span: np.ndarray
    warnings: List[str] = field(default_factory=list)


def _cell_vectors(dimensions: np.ndarray) -> np.ndarray:
    a, b, c, alpha_deg, beta_deg, gamma_deg = dimensions
    if alpha_deg == beta_deg == gamma_deg == 90.0:
        return np.diag([a, b, c])
    alpha, beta, gamma = np.deg2rad([alpha_deg, beta_deg, gamma_deg])
    sin_gamma = np.sin(gamma)
    if abs(sin_gamma) < 1e-10:
        raise BoxValidationError("Box gamma angle produces a singular unit cell")
    avec = np.array([a, 0.0, 0.0])
    bvec = np.array([b * np.cos(gamma), b * sin_gamma, 0.0])
    cx = c * np.cos(beta)
    cy = c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / sin_gamma
    cz2 = c * c - cx * cx - cy * cy
    if cz2 <= 0.0:
        raise BoxValidationError("Box lengths and angles do not define a positive-volume cell")
    cvec = np.array([cx, cy, np.sqrt(cz2)])
    return np.vstack([avec, bvec, cvec])


def validate_box(
    dimensions,
    positions=None,
    *,
    label: str = "periodic box",
    fatal_span_ratio: float = 3.0,
    warning_span_ratio: float = 1.25,
) -> BoxValidationResult:
    dims = np.asarray(dimensions, dtype=float)
    if dims.shape != (6,):
        raise BoxValidationError(
            f"{label} must contain exactly [a, b, c, alpha, beta, gamma]"
        )
    if not np.all(np.isfinite(dims)):
        raise BoxValidationError(f"{label} contains non-finite values: {dims.tolist()}")
    if np.any(dims[:3] <= 0.0):
        raise BoxValidationError(f"{label} lengths must be positive: {dims[:3].tolist()}")
    if np.any(dims[3:] <= 0.0) or np.any(dims[3:] >= 180.0):
        raise BoxValidationError(
            f"{label} angles must be strictly between 0 and 180 degrees: "
            f"{dims[3:].tolist()}"
        )

    vectors = _cell_vectors(dims)
    volume = abs(float(np.linalg.det(vectors)))
    if not np.isfinite(volume) or volume <= 0.0:
        raise BoxValidationError(f"{label} has non-positive volume")

    cell_span = np.sum(np.abs(vectors), axis=0)
    coordinate_span = np.zeros(3, dtype=float)
    notes: List[str] = []
    if positions is not None:
        coords = np.asarray(positions, dtype=float)
        if coords.ndim != 2 or coords.shape[1] != 3 or len(coords) == 0:
            raise BoxValidationError(
                f"{label} coordinate array must be a non-empty N x 3 array"
            )
        if not np.all(np.isfinite(coords)):
            raise BoxValidationError(f"{label} coordinates contain non-finite values")
        coordinate_span = np.ptp(coords, axis=0)
        ratios = np.divide(
            coordinate_span,
            cell_span,
            out=np.zeros(3, dtype=float),
            where=cell_span > 0,
        )
        if np.any(ratios > fatal_span_ratio):
            raise BoxValidationError(
                f"{label} is incompatible with the coordinate bounding box: coordinate "
                f"span {coordinate_span.round(3).tolist()} A versus cell Cartesian span "
                f"{cell_span.round(3).tolist()} A (maximum ratio {ratios.max():.2f})."
            )
        if np.any(ratios > warning_span_ratio):
            notes.append(
                f"{label}: coordinate span {coordinate_span.round(3).tolist()} A exceeds "
                f"the cell Cartesian span {cell_span.round(3).tolist()} A in at least one "
                "direction. This may be valid for unwrapped coordinates; verify the box."
            )
        reverse = np.divide(
            cell_span,
            coordinate_span,
            out=np.full(3, np.inf, dtype=float),
            where=coordinate_span > 0,
        )
        if np.any(reverse > fatal_span_ratio):
            notes.append(
                f"{label}: the cell is more than {fatal_span_ratio:g} times the coordinate "
                "span in at least one direction. This can be valid for a sparse system, "
                "but a solvated-system box may be stale."
            )

    return BoxValidationResult(
        dimensions=dims,
        volume=volume,
        coordinate_span=coordinate_span,
        cell_cartesian_span=cell_span,
        warnings=notes,
    )


IMAGING_NOTE = (
    "The periodic origin is arbitrary: a membrane split across cell edges with "
    "solvent in the middle is physically equivalent. Atom-wise wrapping can "
    "split molecules; this validator translates whole connected molecules only. "
    "For conventional visualization use centering and molecule imaging, e.g. "
    "gmx trjconv -s matching.tpr -f input.gro -o centered.gro -pbc mol -center "
    "with appropriate interactive groups and a matching topology. No external "
    "command was run."
)

# GRO writes coordinates to 0.001 nm (0.01 Å). Two independently rounded
# atoms can change a 3-D separation by at most sqrt(3)*0.01 Å. Repaired
# orientations use a slightly larger verification cutoff so a passing
# in-memory frame cannot acquire a sub-cutoff contact during GRO output.
GRO_ROUNDING_GUARD_ANGSTROM = float(np.sqrt(3) * 0.01 + 0.0005)


def proper_axis_rotations():
    """24 right-handed signed permutations (axis assignments, not all symmetries
    of an unequal cuboid). Column-coordinate convention: x_new = R @ x_old.
    """
    matrices = []
    for perm in itertools.permutations(range(3)):
        for signs in itertools.product((-1, 1), repeat=3):
            matrix = np.eye(3)[list(perm)] * np.array(signs)[:, None]
            if np.linalg.det(matrix) > 0:
                matrices.append(matrix)
    return matrices


class ContactCounter:
    """Exact unique i<j counts with a bounded diagnostic sample.

    Without bonds, same-residue candidates in a covalent-distance window are
    inferred once from direct coordinates, never from a suspect periodic image.
    Ambiguous intramolecular contacts need a matching topology for certainty.
    Residues are identified by resindex, never by possibly repeated residue ID.
    """

    def __init__(self, universe):
        self.u = universe
        self.n = len(universe.atoms)
        self.resindices = np.asarray(universe.atoms.resindices)
        elements = np.array([str(n).lstrip("0123456789")[:1].upper() for n in universe.atoms.names])
        name_elements = elements.copy()
        try:
            provided = np.char.upper(np.asarray(universe.atoms.elements, dtype=str))
            known = provided != ""
            elements = elements.astype("U3")
            elements[known] = provided[known]
        except (AttributeError, NoDataError):
            pass
        try:
            bonds = np.sort(universe.bonds.indices, axis=1)
        except (AttributeError, NoDataError):
            bonds = np.empty((0, 2), dtype=int)
        self.bond_keys = bonds[:, 0].astype(np.int64) * self.n + bonds[:, 1]
        self.policy = "explicit bonds excluded; contacts below 0.1 Å always counted"
        covered = np.zeros(self.n, dtype=bool)
        covered[bonds.ravel()] = True
        if not np.all(covered):
            # Approximate covalent radii (Å); unsupported element/name candidates
            # are not excluded. Any disagreement is reported without renaming.
            radii = {"H": 0.31, "D": 0.31, "T": 0.31, "C": 0.76,
                     "N": 0.71, "O": 0.66, "F": 0.57, "P": 1.07, "S": 1.05}
            radius = np.array([radii.get(e, 0) for e in elements])
            name_radius = np.array([radii.get(e, 0) for e in name_elements])
            coords = np.asarray(universe.atoms.positions, dtype=float)
            tree = cKDTree(coords)
            if (int(tree.count_neighbors(tree, 1.25)) - self.n) // 2 > 5_000_000:
                raise BoxValidationError("Too many bond candidates for bounded-memory preflight")
            pairs = tree.query_pairs(1.25, output_type="ndarray")
            i, j = pairs.T
            d = np.linalg.norm(coords[i] - coords[j], axis=1)
            candidate = np.zeros(len(pairs), dtype=bool)
            for radii_source in (radius, name_radius):
                total_radius = radii_source[i] + radii_source[j]
                candidate |= ((self.resindices[i] == self.resindices[j])
                              & ~(covered[i] & covered[j])
                              & (radii_source[i] > 0) & (radii_source[j] > 0)
                              & (d >= np.maximum(0.65, 0.75 * total_radius))
                              & (d <= total_radius + 0.15))
            self.bond_keys = np.concatenate((self.bond_keys,
                               i[candidate].astype(np.int64) * self.n + j[candidate]))
            self.policy = ("explicit bonds where present; for atoms without connectivity, "
                           "same-residue direct-distance bond candidates "
                           "excluded at max(0.65 Å, 0.75*(r_i+r_j)) <= d <= "
                           "min(1.25 Å, r_i+r_j+0.15 Å), using approximate covalent radii; "
                           "element or atom-name candidates. Ambiguous intramolecular "
                           "contacts require topology; duplicates below 0.1 Å always counted")
            disagreements = int(np.sum(elements != name_elements))
            if disagreements:
                self.policy += (f"; {disagreements} element/name-prefix disagreements "
                                "(some may be valid two-letter elements); metadata was not changed")

    def allowed(self, pairs, distances):
        i, j = pairs.T
        excluded = np.isin(i.astype(np.int64) * self.n + j, self.bond_keys)
        return ~excluded | (distances < 0.1)

    def pairs(self, positions, cutoff, lengths=None):
        coords = np.asarray(positions, dtype=np.float64)
        if lengths is not None:
            coords = np.mod(coords, lengths)
            # Rounding of very small negative values may give exactly L.
            coords = np.minimum(coords, np.nextafter(lengths, 0))
        tree = cKDTree(coords, boxsize=lengths)
        # Refuse pathological quadratic outputs before allocating pair arrays.
        raw_count = (int(tree.count_neighbors(tree, cutoff)) - len(coords)) // 2
        if raw_count > 5_000_000:
            raise BoxValidationError(
                "Contact search exceeds 5,000,000 candidate pairs; severe crowding "
                "or an excessive cutoff prevents bounded-memory verification."
            )
        pairs = tree.query_pairs(cutoff, output_type="ndarray")
        delta = coords[pairs[:, 0]] - coords[pairs[:, 1]]
        if lengths is not None:
            delta -= np.rint(delta / lengths) * lengths
        distances = np.linalg.norm(delta, axis=1)
        keep = (distances < cutoff) & self.allowed(pairs, distances)
        return pairs[keep], distances[keep]

    def count(self, positions, cutoff, lengths=None, pair_filter=None):
        pairs, distances = self.pairs(positions, cutoff, lengths)
        if pair_filter is not None:
            selected = pair_filter(pairs, distances)
            pairs, distances = pairs[selected], distances[selected]
        examples = []
        # Stable, actionable sample; count still includes every pair.
        order = np.lexsort((pairs[:, 1], pairs[:, 0]))[:20]
        for k in order:
            atoms = []
            for index in pairs[k]:
                atom = self.u.atoms[int(index)]
                atoms.append(dict(index=int(index) + 1, name=str(atom.name),
                                  resname=str(atom.resname), resid=int(atom.resid),
                                  resindex=int(atom.resindex) + 1,
                                  chain=str(getattr(atom, "chainID", "")),
                                  segment=str(getattr(atom, "segid", ""))))
            examples.append(dict(atoms=atoms, distance_angstrom=float(distances[k])))
        self_images = 0
        if lengths is not None and min(lengths) < cutoff:
            # cKDTree query_pairs omits i==j, but a too-small cell can put an
            # atom's own nonzero periodic image inside the hard cutoff.
            self_images = self.n
            for index in range(min(self.n, 20 - len(examples))):
                atom = self.u.atoms[index]
                label = dict(index=index + 1, name=str(atom.name), resname=str(atom.resname),
                             resid=int(atom.resid), resindex=int(atom.resindex) + 1,
                             chain=str(getattr(atom, "chainID", "")), segment=str(getattr(atom, "segid", "")))
                examples.append(dict(atoms=[label, label.copy()], distance_angstrom=float(min(lengths)),
                                     self_periodic_image=True))
        return dict(count=len(pairs) + self_images, examples=examples,
                    self_image_atom_count=self_images)


def _coordinate_stats(positions, vectors, allowance):
    fractional = positions @ np.linalg.inv(vectors)
    outside = (fractional < 0) | (fractional >= 1)
    span = np.ptp(positions, axis=0)
    cell_span = np.sum(np.abs(vectors), axis=0)
    return dict(min_angstrom=positions.min(axis=0).tolist(),
                max_angstrom=positions.max(axis=0).tolist(),
                extent_angstrom=span.tolist(), cell_span_angstrom=cell_span.tolist(),
                outside_count=outside.sum(axis=0).tolist(),
                outside_fraction=outside.mean(axis=0).tolist(),
                outside_any_fraction=float(outside.any(axis=1).mean()),
                extent_flag=(span > cell_span + allowance).tolist(),
                large_outside_fraction=bool(outside.any(axis=1).mean() > 0.25))


def render_preflight(result):
    """One renderer shared by console and persisted success reports."""
    lines = ["Periodic-box preflight", f"Classification: {result['classification']} "
             f"({result['status']}); mode={result['mode']}", result.get("reason", ""),
             f"Box source: {result.get('box_source', 'unavailable')}"]
    if "cell_matrix_angstrom" in result:
        lines += [f"Cell matrix (Å, row vectors): {result['cell_matrix_angstrom']}",
                  f"Lengths (Å): {result['lengths_angstrom']}; angles (degrees): {result['angles_degrees']}",
                  f"Off-diagonal entries (Å): {result['off_diagonal_angstrom']}",
                  f"Volume (Å³): {result['volume_angstrom3']}; orthorhombic: {result['orthorhombic']}"]
    for key in ("coordinates_before", "coordinates_after"):
        if key in result:
            s = result[key]
            lines += [f"{key}: min/max/extents (Å) = {s['min_angstrom']} / "
                      f"{s['max_angstrom']} / {s['extent_angstrom']}",
                      f"Cell Cartesian spans (Å): {s['cell_span_angstrom']}; "
                      f"extent flags: {s['extent_flag']}",
                      f"Outside per axis: {s['outside_count']} atoms; fractions {s['outside_fraction']}; "
                      f"outside any face fraction: {s['outside_any_fraction']:.6f}"]
    if "contact_policy" in result:
        lines.append(f"Hard-clash cutoff (Å): {result['hard_clash_cutoff_angstrom']}; {result['contact_policy']}")
    for key in ("nonperiodic_before", "periodic_before", "periodic_after"):
        if key in result:
            lines.append(f"{key}: {result[key]['count']} nonbonded hard-clash pairs")
            for pair in result[key]["examples"]:
                labels = [f"atom {a['index']} {a['name']} residue {a['resname']} {a['resid']} "
                          f"(residue index {a['resindex']}, chain {a['chain']!r}, segment {a['segment']!r})"
                          for a in pair['atoms']]
                lines.append("  " + " -- ".join(labels) + f": {pair['distance_angstrom']:.6f} Å")
    if "rotation_matrix" in result:
        lines += [f"Rotation (dimensionless, x_new = R @ x_old + t): {result['rotation_matrix']}",
                  f"Translation (Å): {result['translation_angstrom']}",
                  f"Search: {result.get('search', {})}",
                  f"Geometry verification: {result.get('invariants', {})}"]
    lines.extend(result.get("notes", []))
    return "\n".join(lines)


class _SearchDeadline(Exception):
    pass


def _place_densest_primary_image(positions, lengths):
    """Translate one rigid frame so the primary cell contains most atoms.

    Molecular-whole coordinates can have an extent larger than a cell, so no
    single translation can put every atom inside.  Per axis, choose the
    length-L half-open slab containing the most coordinates and map that slab
    to [0,L).  Ties prefer a slab whose occupied midpoint is nearest L/2 after
    translation, then the lowest start, making the result deterministic.
    """
    translation = np.zeros(3, dtype=float)
    axes = []
    for axis, length in enumerate(lengths):
        values = np.sort(np.asarray(positions[:, axis], dtype=float))
        right = 0
        candidates = []
        best_count = -1
        for left in range(len(values)):
            if right < left:
                right = left
            while right < len(values) and values[right] < values[left] + length:
                right += 1
            count = right - left
            if count > best_count:
                best_count, candidates = count, [(left, right)]
            elif count == best_count:
                candidates.append((left, right))
        def tie_key(pair):
            left, right = pair
            occupied_midpoint = (values[left] + values[right - 1]) / 2
            # Mapping start=values[left] to zero gives this midpoint.
            return (abs((occupied_midpoint - values[left]) - length / 2), values[left])
        left, right = min(candidates, key=tie_key)
        start = values[left]
        translation[axis] = -start
        axes.append(dict(axis="xyz"[axis], slab_start_angstrom=float(start),
                         slab_end_angstrom=float(start + length),
                         atoms_inside=int(best_count),
                         fraction_inside=float(best_count / len(values))))
    return translation, dict(method="densest half-open cell-width slab per axis",
                             axes=axes)


def search_orientation(positions, lengths, counter, spec):
    """Global hull/random seeds, all axis assignments, then 3-D refinement.

    Trimmed enclosing widths have a useful slope far from the correct frame;
    local contact refinement is used only after global shape fitting. This is
    a bounded search, not a proof that a failed box cannot be repaired.
    """
    start = time.monotonic()
    rng = np.random.default_rng(spec.random_seed)
    center = positions.mean(axis=0)
    centered = positions - center
    sample = centered[np.sort(rng.choice(len(positions),
                       min(len(positions), spec.search_sample_size), replace=False))]
    # A smaller global sample keeps many orientations cheap; full sample refines.
    coarse = sample[np.linspace(0, len(sample) - 1, min(1500, len(sample)), dtype=int)]
    stats = dict(sample_atoms=len(sample), coarse_sample_atoms=len(coarse),
                 objective_evaluations=0, candidate_verifications=0,
                 axis_assignments=24, deadline_reached=False, random_seed=spec.random_seed,
                 max_search_seconds=spec.max_search_seconds, starts=spec.search_starts,
                 max_iterations=spec.search_max_iterations,
                 output_rounding_guard_angstrom=GRO_ROUNDING_GUARD_ANGSTROM,
                 guarded_verification_cutoff_angstrom=(
                     spec.hard_clash_cutoff + GRO_ROUNDING_GUARD_ANGSTROM
                 ))
    guarded_cutoff = spec.hard_clash_cutoff + GRO_ROUNDING_GUARD_ANGSTROM

    def deadline():
        if time.monotonic() - start > spec.max_search_seconds:
            raise _SearchDeadline

    def objective(matrix, points):
        deadline()
        stats["objective_evaluations"] += 1
        projected = points @ matrix.T
        # Trim a small fraction of protruding molecular atoms. Minimize widths,
        # not just zero excess, so sparse samples don't give a flat plateau.
        lo, hi = np.quantile(projected, [0.005, 0.995], axis=0)
        ratio = (hi - lo) / lengths
        return float(np.sum(ratio * ratio) + 4 * np.sum(np.maximum(ratio - 1, 0)**2))

    candidates = []
    try:
        seeds = [np.eye(3)]
        seeds.extend(Rotation.random(max(24, spec.search_starts * 4), random_state=rng).as_matrix())
        # Large hull faces reveal cell axes without depending on solute inertia.
        try:
            hull = ConvexHull(sample)
            triangles = sample[hull.simplices]
            areas = np.linalg.norm(np.cross(triangles[:, 1] - triangles[:, 0],
                                            triangles[:, 2] - triangles[:, 0]), axis=1)
            normals = hull.equations[np.argsort(-areas)[:32], :3]
            for i, a in enumerate(normals):
                for b in normals[i + 1:]:
                    if abs(a @ b) > 0.3:
                        continue
                    b = b - a * (a @ b)
                    b /= np.linalg.norm(b)
                    seeds.append(np.vstack((a, b, np.cross(a, b))))
        except QhullError:
            pass
        # Axis assignments are evaluated explicitly. Keep only diverse seeds.
        ranked = []
        for seed in seeds:
            for assignment in proper_axis_rotations():
                matrix = assignment @ seed
                ranked.append((objective(matrix, coarse), matrix))
        ranked.sort(key=lambda x: x[0])
        chosen = []
        for score, matrix in ranked:
            # Equivalent sign flips need no duplicate local optimization.
            if any(np.max(np.abs(np.abs(matrix @ other.T) - np.eye(3))) < 0.08
                   for other in chosen):
                continue
            chosen.append(matrix)
            if len(chosen) >= spec.search_starts:
                break
        for matrix in chosen:
            candidates.append((objective(matrix, sample), matrix))
            fitted = minimize(lambda v: objective(Rotation.from_rotvec(v).as_matrix() @ matrix, sample),
                              np.zeros(3), method="Powell",
                              options={"maxiter": spec.search_max_iterations, "xtol": 1e-7, "ftol": 1e-9})
            fitted_matrix = Rotation.from_rotvec(fitted.x).as_matrix() @ matrix
            candidates.append((objective(fitted_matrix, sample), fitted_matrix))
    except _SearchDeadline:
        stats["deadline_reached"] = True

    candidates.sort(key=lambda x: x[0])
    # Leave no partial transform behind: all verification operates on copies.
    for score, matrix in candidates[:spec.max_candidate_verifications]:
        transformed = centered @ matrix.T
        stats["candidate_verifications"] += 1
        guarded_count = counter.count(transformed, guarded_cutoff, lengths)
        if guarded_count["count"] == 0:
            placement, placement_stats = _place_densest_primary_image(transformed, lengths)
            stats.update(objective_value=score, elapsed_seconds=time.monotonic() - start)
            stats["primary_cell_placement"] = placement_stats
            return (matrix, -matrix @ center + placement,
                    counter.count(transformed, spec.hard_clash_cutoff, lengths), stats)
        # Refine a near solution using periodic contact residuals. The variables
        # remain exactly three global rotation parameters, never atomic positions.
        if not stats["deadline_reached"]:
            try:
                from scipy.optimize import least_squares
                for _ in range(3):
                    deadline()
                    pairs, _ = counter.pairs(transformed, guarded_cutoff + 0.4, lengths)
                    if not len(pairs):
                        break
                    deltas = centered[pairs[:, 0]] - centered[pairs[:, 1]]
                    def residual(v):
                        deadline()
                        moved = deltas @ (Rotation.from_rotvec(v).as_matrix() @ matrix).T
                        moved -= np.rint(moved / lengths) * lengths
                        return np.maximum(guarded_cutoff + 0.002 - np.linalg.norm(moved, axis=1), 0)
                    fitted = least_squares(residual, np.zeros(3), max_nfev=40,
                                           bounds=(-0.035, 0.035), ftol=1e-10, gtol=1e-10, xtol=1e-10)
                    matrix = Rotation.from_rotvec(fitted.x).as_matrix() @ matrix
                    transformed = centered @ matrix.T
                    stats["candidate_verifications"] += 1
                    guarded_count = counter.count(transformed, guarded_cutoff, lengths)
                    if guarded_count["count"] == 0:
                        placement, placement_stats = _place_densest_primary_image(transformed, lengths)
                        stats.update(objective_value=objective(matrix, sample),
                                     elapsed_seconds=time.monotonic() - start)
                        stats["primary_cell_placement"] = placement_stats
                        return (matrix, -matrix @ center + placement,
                                counter.count(transformed, spec.hard_clash_cutoff, lengths), stats)
            except _SearchDeadline:
                stats["deadline_reached"] = True
    stats["elapsed_seconds"] = time.monotonic() - start
    return None, None, None, stats


def _wrap_molecules(universe, coords, lengths):
    try:
        if len(universe.bonds) == 0:
            raise ValueError("empty connectivity")
        groups = universe.atoms.fragindices
        # Incomplete CONECT records often cover only a ligand. Do not treat
        # unconnected atoms elsewhere as independently wrappable molecules.
        resindices = universe.atoms.resindices
        starts = np.unique(resindices, return_index=True)[1]
        first_fragments = groups[starts]
        if np.any(groups != first_fragments[resindices]):
            raise ValueError("connectivity splits at least one residue")
    except (AttributeError, ValueError, NoDataError):
        # A residue is not a molecule: declining optional wrapping is safer than
        # translating parts of a protein independently. Rotation needs no bonds.
        return coords, np.zeros_like(coords), "Wrapping skipped: no reliable connectivity; whole-system imaging retained."
    _, inverse = np.unique(groups, return_inverse=True)
    counts = np.bincount(inverse)
    centroids = np.column_stack([np.bincount(inverse, weights=coords[:, a]) / counts for a in range(3)])
    shifts = -np.floor(centroids / lengths)[inverse] * lengths
    return coords + shifts, shifts, "Wrapped connected fragments by integer cell vectors; no molecules were made whole atom by atom."


def preflight_environment(universe, spec=None, *, box_source="file", cell_vectors=None,
                          parser_notes=(), stream=None):
    """Validate/repair in memory; return a JSON-ready report or raise with it.

    Success attaches the report to the original Universe for pipeline reporting.
    Strict mode can search diagnostically but never applies a recovered rotation.
    """
    from config import BoxValidationSpec
    spec = spec or BoxValidationSpec()
    spec.validate()
    result = dict(mode=spec.mode, classification="(e) no_usable_box", status="refused",
                  box_source=box_source, notes=list(parser_notes) + [IMAGING_NOTE],
                  hard_clash_cutoff_angstrom=spec.hard_clash_cutoff)

    def refuse(label, reason):
        result.update(classification=label, status="refused", reason=reason)
        print(render_preflight(result), file=sys.stderr if stream is None else stream)
        raise BoxValidationError(f"{label}: {reason}", result)

    if spec.mode == "off":
        result.update(classification="skipped", status="unchecked", reason=
                      "EXPLICIT OPT-OUT: environment periodicity was not validated; downstream PBC remains unchecked.")
        universe.box_preflight = result
        print(render_preflight(result), file=stream)
        return result
    try:
        validation = validate_box(universe.dimensions)
    except (ValueError, TypeError) as exc:
        refuse("(e) no_usable_box", str(exc))
    dims = validation.dimensions
    vectors = _cell_vectors(dims) if cell_vectors is None else np.asarray(cell_vectors, dtype=float)
    orthorhombic = bool(np.allclose(vectors, np.diag(np.diag(vectors)), atol=1e-6, rtol=0)
                        and np.all(np.diag(vectors) > 0))
    result.update(cell_matrix_angstrom=vectors.tolist(), lengths_angstrom=dims[:3].tolist(),
                  angles_degrees=dims[3:].tolist(), volume_angstrom3=validation.volume,
                  off_diagonal_angstrom=vectors[~np.eye(3, dtype=bool)].tolist(), orthorhombic=orthorhombic)
    positions = np.asarray(universe.atoms.positions, dtype=float).copy()
    if not len(positions) or not np.isfinite(positions).all():
        refuse("(d) geometry_broken", "Coordinates must be nonempty and finite (Å).")
    result["coordinates_before"] = _coordinate_stats(positions, vectors, spec.molecule_protrusion_allowance)
    if not orthorhombic:
        refuse("(e) no_usable_box", "Unsupported triclinic or non-axis-aligned cell: full matrix reported; no diagonal approximation was used.")
    try:
        counter = ContactCounter(universe)
    except BoxValidationError as exc:
        refuse("indeterminate", str(exc))
    result["contact_policy"] = counter.policy
    lengths = dims[:3]
    try:
        result["nonperiodic_before"] = counter.count(positions, spec.hard_clash_cutoff)
        result["periodic_before"] = counter.count(positions, spec.hard_clash_cutoff, lengths)
    except BoxValidationError as exc:
        refuse("indeterminate", str(exc))
    if result["nonperiodic_before"]["count"]:
        refuse("(d) geometry_broken", "Nonbonded hard clashes exist without PBC; a cell rotation cannot remove them.")
    if result["periodic_before"]["self_image_atom_count"]:
        refuse("(c) box_wrong", "A cell vector is shorter than the hard-clash cutoff (Å); each atom clashes with its own periodic image regardless of orientation.")
    matrix, translation = np.eye(3), np.zeros(3)
    result["classification"] = "(a) consistent"
    result["reason"] = "No nonbonded hard clashes under PBC; outside-cell coordinates alone do not indicate failure."
    if result["periodic_before"]["count"]:
        try:
            matrix, translation, after, stats = search_orientation(positions, lengths, counter, spec)
        except BoxValidationError as exc:
            refuse("indeterminate", str(exc))
        result["search"] = stats
        if matrix is None:
            refuse("(c) box_wrong", "Suspected stale/mismatched box or unresolved orientation: bounded orientation search found no zero-clash solution. This is not proof that no rotation exists; supply an authoritative unaligned frame/box.")
        result.update(classification="(b) orientation_broken", rotation_matrix=matrix.tolist(),
                      translation_angstrom=translation.tolist(), periodic_after=after,
                      reason="A single global rigid rotation yields zero full-system periodic hard clashes; orientation repair is supported by this test, not uniquely inferred from coordinates.")
        if spec.mode == "strict":
            refuse("(b) orientation_broken", "Strict mode refuses a frame requiring rotation. A zero-clash candidate was found but not applied; use repair mode to permit it.")
    transformed = positions.copy() if result["classification"] == "(a) consistent" else positions @ matrix.T + translation
    shifts = np.zeros_like(transformed)
    if spec.wrap_after_validation:
        transformed, shifts, note = _wrap_molecules(universe, transformed, lengths)
        result["notes"].append(note)
    lattice_units = shifts / lengths
    unique_shifts, shift_counts = np.unique(np.rint(lattice_units).astype(np.int64), axis=0, return_counts=True)
    result["wrapping"] = dict(requested=spec.wrap_after_validation,
                              translated_atoms=int(np.any(shifts != 0, axis=1).sum()),
                              lattice_shift_atom_counts=[dict(cell_vector_multiples=s.tolist(), atoms=int(n))
                                                         for s, n in zip(unique_shifts, shift_counts)])
    # Check transform identity for every atom, plus deterministic sampled pair
    # distances, after undoing allowed integer lattice translations.
    unshifted = transformed - shifts
    expected = positions @ matrix.T + translation
    residual = float(np.max(np.abs(unshifted - expected)))
    indices = np.linspace(0, len(positions) - 1, min(4096, len(positions)), dtype=int)
    partners = np.roll(indices, 1)
    pair_error = float(np.max(np.abs(np.linalg.norm(unshifted[indices] - unshifted[partners], axis=1)
                                    - np.linalg.norm(positions[indices] - positions[partners], axis=1))))
    if (not np.allclose(matrix @ matrix.T, np.eye(3), atol=1e-10, rtol=0)
            or abs(np.linalg.det(matrix) - 1) > 1e-10 or residual > 1e-8 or pair_error > 1e-8
            or not np.allclose(lattice_units, np.rint(lattice_units), atol=1e-10, rtol=0)):
        refuse("indeterminate", "Rigid-transform geometry invariant failed (Å tolerance 1e-8).")
    # MDAnalysis stores float32 coordinates: verify that actual stored precision
    # also passes, before assigning anything to the caller's Universe.
    stored = transformed.astype(universe.atoms.positions.dtype).astype(float)
    result["invariants"] = dict(transform_residual_angstrom=residual,
                                sampled_pair_error_angstrom=pair_error,
                                storage_rounding_max_angstrom=float(np.max(np.abs(stored - transformed))),
                                determinant=float(np.linalg.det(matrix)), cell_changed=False)
    try:
        result["periodic_after"] = counter.count(stored, spec.hard_clash_cutoff, lengths)
    except BoxValidationError as exc:
        refuse("indeterminate", str(exc))
    if result["periodic_after"]["count"]:
        refuse(result["classification"], "Candidate failed zero-clash verification at coordinate storage precision; nothing applied.")
    result["coordinates_after"] = _coordinate_stats(stored, vectors, spec.molecule_protrusion_allowance)
    result["status"] = "repaired" if result["classification"] == "(b) orientation_broken" else "passed"
    result["rotation_matrix"] = matrix.tolist()
    result["translation_angstrom"] = translation.tolist()
    if not np.array_equal(stored, positions):
        universe.atoms.positions = stored
    universe.box_preflight = result
    print(render_preflight(result), file=stream)
    return result


def read_environment(path, box_dimensions=None, spec=None, *, format="auto"):
    """PDB/GRO reader with explicit box provenance and no sibling-box fallback."""
    import MDAnalysis as mda
    from config import BoxValidationSpec
    spec = spec or BoxValidationSpec()
    notes, raw_vectors = [], None
    path = Path(path)
    file_format = path.suffix.lower().lstrip(".") if format == "auto" else format.lower()
    try:
        if file_format == "pdb":
            raw_names = []
            with path.open() as handle:
                for line in handle:
                    if line.startswith(("ATOM  ", "HETATM")):
                        raw_names.append(line[17:21].strip())
                    if line.startswith("ENDMDL"):
                        break
            # No guessed restoration: a generic note covers unknowable aliases.
            if any(len(name) == 3 for name in raw_names):
                notes.append("Three-character PDB residue names may be legitimate or already truncated. Ambiguous names require topology-assisted name_restoration; no names were guessed.")
        elif file_format == "gro":
            with path.open() as handle:
                next(handle)
                n_atoms = int(next(handle))
                for _ in range(n_atoms):
                    next(handle)
                values = np.array([float(x) for x in next(handle).split()]) * 10
            if len(values) not in (3, 9):
                raise ValueError("GRO box line must have three or nine values (nm in file)")
            raw_vectors = np.diag(values[:3])
            if len(values) == 9:
                raw_vectors = np.array([[values[0], values[3], values[4]],
                                        [values[5], values[1], values[6]],
                                        [values[7], values[8], values[2]]])
        else:
            raise ValueError("Environment input must be PDB or GRO")
        u = mda.Universe(str(path), format=file_format.upper(), topology_format=file_format.upper())
        if file_format == "pdb":
            if len(raw_names) != len(u.atoms):
                raise ValueError("PDB atom/name count mismatch")
            raw_names = np.asarray(raw_names, dtype=str)
            resindices = u.atoms.resindices
            starts = np.unique(resindices, return_index=True)[1]
            names = raw_names[starts]
            if np.any(raw_names != names[resindices]):
                raise ValueError("PDB residue contains conflicting four-character names")
            u.residues.resnames = names
        source = "file"
        if box_dimensions is not None:
            validate_box(box_dimensions)
            if u.dimensions is not None and not np.allclose(u.dimensions, box_dimensions, atol=1e-4, rtol=0):
                raise ValueError("Explicit box differs from the file box. Changing dimensions is not a repair; provide an authoritative corrected input instead.")
            if u.dimensions is None:
                u.dimensions = np.asarray(box_dimensions, dtype=np.float32)
                source = "explicit YAML override"
                if raw_vectors is None or not np.any(raw_vectors):
                    raw_vectors = _cell_vectors(np.asarray(box_dimensions, dtype=float))
            else:
                source = "file (matching explicit YAML box confirmed; file values retained)"
    except (ValueError, TypeError, OSError, StopIteration) as exc:
        result = dict(mode=spec.mode, classification="(e) no_usable_box", status="refused",
                      box_source="unreadable input or override", reason=str(exc), notes=notes)
        print(render_preflight(result), file=sys.stderr)
        raise BoxValidationError(f"(e) no_usable_box: {exc}", result) from exc
    preflight_environment(u, spec, box_source=source, cell_vectors=raw_vectors, parser_notes=notes)
    return u


def validate_staged_hard_clashes(pdb_path, gro_path, spec, clash_spec=None):
    """Re-read both staged formats and enforce the final all-atom PBC pass."""
    import MDAnalysis as mda
    if spec.mode == "off":
        return dict(status="unchecked", reason="box_validation.mode is off")
    checked = {}
    for label, path in (("pdb", pdb_path), ("gro", gro_path)):
        universe = mda.Universe(path)
        validation = validate_box(universe.dimensions)
        vectors = _cell_vectors(validation.dimensions)
        if not np.allclose(vectors, np.diag(np.diag(vectors)), atol=1e-6, rtol=0):
            raise BoxValidationError(
                f"Final staged {label.upper()} has an unsupported non-orthorhombic cell"
            )
        positions = universe.atoms.positions.astype(float)
        counter = ContactCounter(universe)
        contacts = counter.count(
            positions, spec.hard_clash_cutoff, validation.dimensions[:3]
        )
        allowed_contacts = None
        if clash_spec is not None:
            from classify import classify_resnames
            mapping = classify_resnames(universe.atoms.resnames)
            classes = np.array([mapping[rn] for rn in universe.atoms.resnames])
            cutoffs = np.array([min(spec.hard_clash_cutoff, clash_spec.threshold_for(c))
                                for c in classes])
            protected = np.isin(classes, clash_spec.keep_classes)

            def permitted(pairs, distances):
                i, j = pairs.T
                return (protected[i] | protected[j]
                        | (distances >= np.minimum(cutoffs[i], cutoffs[j])))

            allowed_contacts = counter.count(
                positions, spec.hard_clash_cutoff, validation.dimensions[:3], permitted
            )
            contacts = counter.count(
                positions, spec.hard_clash_cutoff, validation.dimensions[:3],
                lambda pairs, distances: ~permitted(pairs, distances),
            )
            # Self-image contacts indicate an invalid cell, never a class exception.
            allowed_contacts['count'] -= allowed_contacts['self_image_atom_count']
            allowed_contacts['self_image_atom_count'] = 0
            allowed_contacts['examples'] = [e for e in allowed_contacts['examples']
                                            if not e.get('self_periodic_image')]
            if allowed_contacts['count']:
                print(f"      WARNING: {label.upper()}: {allowed_contacts['count']} contacts "
                      f"below {spec.hard_clash_cutoff:g} Å accepted by per-class cutoffs "
                      "or keep_classes; coordinate output permitted (see report).")
        wrapped_copy = np.mod(positions, validation.dimensions[:3])
        periodic_equivalent_outside = np.any(
            (wrapped_copy < 0) | (wrapped_copy >= validation.dimensions[:3]), axis=1
        )
        checked[label] = dict(
            hard_clashes=contacts,
            accepted_contacts=allowed_contacts,
            coordinates=_coordinate_stats(
                positions, vectors, spec.molecule_protrusion_allowance
            ),
            periodic_equivalent_outside_count=int(periodic_equivalent_outside.sum()),
            box_lengths_angstrom=validation.dimensions[:3].tolist(),
        )
        if contacts["count"]:
            example = contacts["examples"][0]
            indices = [atom["index"] for atom in example["atoms"]]
            raise BoxValidationError(
                f"Final staged {label.upper()} has {contacts['count']} nonbonded "
                f"PBC hard clash(es) below {spec.hard_clash_cutoff:g} Å; "
                f"first pair uses 1-based atom indices {indices}. No coordinate "
                "output was promoted."
            )
    accepted = any((v.get('accepted_contacts') or {}).get('count', 0) for v in checked.values())
    return dict(status="passed_with_exceptions" if accepted else "passed",
                keep_classes=list(clash_spec.keep_classes) if clash_spec else [],
                class_thresholds=dict(clash_spec.thresholds) if clash_spec else {},
                hard_clash_cutoff_angstrom=spec.hard_clash_cutoff,
                formats=checked)
