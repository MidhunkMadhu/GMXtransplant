"""Runtime support for the exported standalone OpenMM runner, not assembly."""

from __future__ import annotations

import json
import math
import os
from pathlib import Path
import tempfile
import time
import warnings

import numpy as np

from config import MinimizationSpec
from output import _promote_staged_outputs
from minimization_restraints import add_heavy_atom_restraints


class MinimizationError(RuntimeError):
    pass


def _load_openmm():
    try:
        import openmm
        from openmm import app, unit
    except (ImportError, ModuleNotFoundError) as exc:
        raise MinimizationError(
            "OpenMM is required to execute this standalone runner. Install the optional dependency "
            "from the repository with `python -m pip install '.[minimize]'`, "
            "install it directly with `python -m pip install openmm`, or use "
            "`conda install -c conda-forge openmm`."
        ) from exc
    return openmm, app, unit


def _force_statistics(state, unit) -> tuple[float, float]:
    values = np.asarray(
        state.getForces(asNumpy=True).value_in_unit(
            unit.kilojoule_per_mole / unit.nanometer
        ),
        dtype=np.float64,
    )
    if values.ndim != 2 or values.shape[1] != 3 or not np.all(np.isfinite(values)):
        raise MinimizationError("OpenMM returned non-finite or malformed forces")
    rms_components = float(np.sqrt(np.mean(values * values)))
    max_atom_norm = float(np.max(np.linalg.norm(values, axis=1))) if len(values) else 0.0
    return rms_components, max_atom_norm


def _energy_kj_mol(state, unit) -> float:
    value = float(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
    if not math.isfinite(value):
        raise MinimizationError(
            "OpenMM returned a non-finite potential energy; the coordinates or "
            "topology likely contain an unresolved overlap or invalid parameter"
        )
    return value


def _platform_candidates(openmm, requested: str):
    available = {
        openmm.Platform.getPlatform(index).getName():
        openmm.Platform.getPlatform(index)
        for index in range(openmm.Platform.getNumPlatforms())
    }
    if requested != "auto":
        if requested not in available:
            raise MinimizationError(
                f"Requested OpenMM platform '{requested}' is unavailable; "
                f"available platforms: {', '.join(available) or '(none)'}"
            )
        return [(requested, available[requested])]
    priority = ("CUDA", "HIP", "OpenCL", "CPU", "Reference")
    return [(name, available[name]) for name in priority if name in available]


def _create_simulation(openmm, app, unit, topology, system, spec: MinimizationSpec):
    failures = []
    for name, platform in _platform_candidates(openmm, spec.platform):
        integrator = openmm.VerletIntegrator(0.001 * unit.picoseconds)
        integrator.setConstraintTolerance(spec.constraint_tolerance)
        property_names = set(platform.getPropertyNames())
        properties = {}
        if "Precision" in property_names:
            properties["Precision"] = spec.precision
        if name == "CPU" and "Threads" in property_names and os.environ.get("OPENMM_CPU_THREADS"):
            properties["Threads"] = os.environ["OPENMM_CPU_THREADS"]
        if spec.device_index:
            if "DeviceIndex" not in property_names:
                if spec.platform != "auto":
                    raise MinimizationError(
                        f"OpenMM platform '{name}' does not support device_index"
                    )
            else:
                properties["DeviceIndex"] = spec.device_index
        try:
            simulation = app.Simulation(
                topology, system, integrator, platform, properties
            )
        except Exception as exc:
            failures.append(f"{name}: {exc}")
            del integrator
            if spec.platform != "auto":
                break
            continue
        return simulation, integrator, name, properties, failures
    details = "; ".join(failures) or "no OpenMM platforms were registered"
    raise MinimizationError(f"Could not create an OpenMM execution context: {details}")


def _make_progress_reporter(openmm):
    class ProgressReporter(openmm.MinimizationReporter):
        def __init__(self):
            super().__init__()
            self.callback_count = 0
            self.restraint_stages = 1
            self.last_iteration = None
            self.last_objective_rms = None
            self.last_statistics = {}

        def report(self, iteration, x, grad, args):
            if self.last_iteration is not None and iteration <= self.last_iteration:
                self.restraint_stages += 1
            self.last_iteration = int(iteration)
            self.callback_count += 1
            gradient = np.asarray(grad, dtype=np.float64)
            self.last_objective_rms = float(np.sqrt(np.mean(gradient * gradient)))
            self.last_statistics = {str(key): float(value) for key, value in args.items()}
            return False

    return ProgressReporter()


def _write_minimized_gro(input_gro: str, output_path: str, positions_angstrom) -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        import MDAnalysis as mda

        universe = mda.Universe(input_gro)
    positions = np.asarray(positions_angstrom, dtype=np.float64)
    if positions.shape != (len(universe.atoms), 3):
        raise MinimizationError(
            f"Minimized coordinate shape {positions.shape} does not match the "
            f"input GRO atom count {len(universe.atoms)}"
        )
    if not np.all(np.isfinite(positions)):
        raise MinimizationError("Minimized coordinates contain NaN or infinity")
    universe.atoms.positions = positions.astype(np.float32)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        universe.atoms.write(output_path)
        check = mda.Universe(output_path)
    if len(check.atoms) != len(universe.atoms):
        raise MinimizationError(
            f"Staged minimized GRO has {len(check.atoms)} atoms; expected "
            f"{len(universe.atoms)}"
        )
    if list(check.atoms.names) != list(universe.atoms.names):
        raise MinimizationError("Staged minimized GRO changed atom names or order")
    if list(check.atoms.resnames) != list(universe.atoms.resnames):
        raise MinimizationError("Staged minimized GRO changed residue names or order")
    if not np.allclose(check.atoms.positions, positions, atol=0.0052, rtol=0):
        raise MinimizationError(
            "Staged minimized GRO changed coordinates beyond GRO rounding precision"
        )
    if universe.dimensions is not None and (
        check.dimensions is None
        or not np.allclose(
            check.dimensions,
            universe.dimensions,
            atol=np.array([0.001, 0.001, 0.001, 0.02, 0.02, 0.02]),
            rtol=0,
        )
    ):
        raise MinimizationError("Staged minimized GRO lost or changed the periodic box")


def _validate_coordinate_topology_identity(input_gro: str, topology) -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        import MDAnalysis as mda

        universe = mda.Universe(input_gro)
    coordinate_names = [str(name) for name in universe.atoms.names]
    coordinate_resnames = [str(name) for name in universe.atoms.resnames]
    # Public OpenMM atom names are normalized (for example HN -> H). Audit
    # against the source records actually parsed, without renaming coordinates.
    try:
        records = [atom for name, count in topology._molecules
                   for _ in range(count) for atom in topology._moleculeTypes[name].atoms]
        topology_names = [atom[4] for atom in records]
        topology_resnames = [atom[3] for atom in records]
    except (AttributeError, KeyError, IndexError, TypeError) as exc:
        raise MinimizationError("OpenMM source atom records unavailable; incompatible parser version") from exc
    if len(records) != len(coordinate_names):
        raise MinimizationError("GRO/topology source atom counts differ")
    for index, (coord_name, top_name, coord_resname, top_resname) in enumerate(
        zip(
            coordinate_names,
            topology_names,
            coordinate_resnames,
            topology_resnames,
        ),
        1,
    ):
        if coord_name != top_name or coord_resname != top_resname:
            raise MinimizationError(
                "GRO/topology atom-order mismatch at atom "
                f"{index}: GRO={coord_resname}:{coord_name}, "
                f"topology={top_resname}:{top_name}"
            )


def _render_text_report(report: dict) -> str:
    before = report["before"]
    after = report["after"]
    settings = report["settings"]
    progress = report["minimizer_progress"]
    lines = [
        "GMXtransplant OpenMM minimization report",
        "=" * 46,
        f"Status: {report['status']}",
        f"Coordinates: {report['inputs']['coordinates']}",
        f"Topology: {report['inputs']['topology']}",
        f"Output GRO: {report['outputs']['gro']}",
        f"OpenMM version: {report['openmm']['version']}",
        f"Platform: {report['openmm']['platform']}",
        f"Platform properties: {report['openmm']['properties']}",
        f"Positional restraints: {report['positional_restraints']}",
        f"Energy definition: {report['energy_definition']}",
        f"Algorithm: {settings['algorithm']}",
        f"Nonbonded method: {settings['nonbonded_method']}",
        f"Cutoff / switch (nm): {settings['nonbonded_cutoff_nm']} / "
        f"{settings['switch_distance_nm']}",
        f"Constraints / rigid water: {settings['constraints']} / "
        f"{settings['rigid_water']}",
        f"Tolerance (kJ mol^-1 nm^-1): {settings['tolerance_kj_mol_nm']}",
        f"Maximum iterations: {settings['max_iterations']}",
        f"Minimizer callbacks / constraint stages: "
        f"{progress['reporter_callbacks']} / {progress['constraint_restraint_stages']}",
        f"Last objective RMS gradient (kJ mol^-1 nm^-1): "
        f"{progress['last_objective_rms_kj_mol_nm']}",
        f"Potential energy before (kJ/mol): {before['potential_energy_kj_mol']:.8g}",
        f"Potential energy after  (kJ/mol): {after['potential_energy_kj_mol']:.8g}",
        f"RMS force before (kJ mol^-1 nm^-1): {before['rms_force_kj_mol_nm']:.8g}",
        f"RMS force after  (kJ mol^-1 nm^-1): {after['rms_force_kj_mol_nm']:.8g}",
        f"Maximum atom-force norm after (kJ mol^-1 nm^-1): "
        f"{after['max_force_kj_mol_nm']:.8g}",
        f"Elapsed seconds: {report['elapsed_seconds']:.3f}",
    ]
    for warning in report.get("warnings", []):
        lines.append(f"WARNING: {warning}")
    return "\n".join(lines) + "\n"


def run_minimization(spec: MinimizationSpec) -> dict:
    """Minimize one matching GRO/topology pair and atomically write results."""
    spec.validate(check_paths=True, require_inputs=True)
    output_gro = os.path.abspath(spec.output_gro_path)
    report_txt = os.path.abspath(spec.report_path + ".txt")
    report_json = os.path.abspath(spec.report_path + ".json")
    input_gro = os.path.abspath(spec.coordinates_path)
    topology_path = os.path.abspath(spec.topology_path)
    for destination in (output_gro, report_txt, report_json):
        if os.path.realpath(destination) in {
            os.path.realpath(input_gro), os.path.realpath(topology_path)
        }:
            raise MinimizationError(
                f"Minimization output must not overwrite an input: {destination}"
            )

    openmm, app, unit = _load_openmm()
    start = time.monotonic()
    try:
        gro = app.GromacsGroFile(input_gro)
        box_vectors = gro.getPeriodicBoxVectors()
        include_dir = os.path.abspath(spec.include_dir) if spec.include_dir else str(
            Path(topology_path).parent
        )
        # OpenMM 8.6 currently leaves its internal topology file handle for
        # garbage collection, which can emit a harmless ResourceWarning under
        # unittest's warning settings. It does not affect parsed data.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ResourceWarning)
            top = app.GromacsTopFile(
                topology_path,
                periodicBoxVectors=box_vectors,
                includeDir=include_dir,
                defines=spec.defines,
            )
        coordinate_count = len(gro.getPositions())
        topology_count = top.topology.getNumAtoms()
        if coordinate_count != topology_count:
            raise MinimizationError(
                f"GRO/topology atom-count mismatch: coordinates={coordinate_count}, "
                f"topology={topology_count}"
            )
        _validate_coordinate_topology_identity(input_gro, top)

        method = getattr(app, spec.nonbonded_method)
        constraints = {
            "none": None,
            "h_bonds": app.HBonds,
            "all_bonds": app.AllBonds,
            "h_angles": app.HAngles,
        }[spec.constraints]
        switch_distance = (
            None
            if spec.switch_distance_nm in (None, 0)
            else spec.switch_distance_nm * unit.nanometer
        )
        system = top.createSystem(
            nonbondedMethod=method,
            nonbondedCutoff=spec.nonbonded_cutoff_nm * unit.nanometer,
            constraints=constraints,
            rigidWater=spec.rigid_water,
            ewaldErrorTolerance=spec.ewald_error_tolerance,
            switchDistance=switch_distance,
            useDispersionCorrection=spec.use_dispersion_correction,
        )
        restraints, restrained_indices = add_heavy_atom_restraints(
            openmm, unit, top, system, gro.getPositions(), spec)
        simulation, integrator, platform_name, properties, platform_failures = (
            _create_simulation(openmm, app, unit, top.topology, system, spec)
        )
        simulation.context.setPeriodicBoxVectors(*box_vectors)
        simulation.context.setPositions(gro.getPositions())
        initial_state = simulation.context.getState(getEnergy=True, getForces=True)
        initial_energy = _energy_kj_mol(initial_state, unit)
        initial_rms, initial_max = _force_statistics(initial_state, unit)

        progress = _make_progress_reporter(openmm)
        simulation.minimizeEnergy(
            tolerance=(
                spec.tolerance_kj_mol_nm
                * unit.kilojoule_per_mole
                / unit.nanometer
            ),
            maxIterations=spec.max_iterations,
            reporter=progress,
        )
        final_state = simulation.context.getState(
            getEnergy=True, getForces=True, getPositions=True
        )
        final_energy = _energy_kj_mol(final_state, unit)
        final_rms, final_max = _force_statistics(final_state, unit)
        positions = final_state.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
        displacement = (positions - np.asarray(gro.getPositions().value_in_unit(unit.angstrom)))[restrained_indices]
        if restraints["periodic"]:
            box = np.asarray(box_vectors.value_in_unit(unit.angstrom))
            fractional = displacement @ np.linalg.inv(box)
            displacement = (fractional - np.round(fractional)) @ box
        distances = np.linalg.norm(displacement, axis=1)
        restraints["rms_displacement_angstrom"] = float(np.sqrt(np.mean(distances ** 2)))
        restraints["max_displacement_angstrom"] = float(np.max(distances))
    except MinimizationError:
        raise
    except Exception as exc:
        raise MinimizationError(f"OpenMM minimization failed: {exc}") from exc

    # The reporter's gradient is for OpenMM's actual minimization objective,
    # including temporary forces used to satisfy constraints. Raw physical
    # forces from State are valuable diagnostics but are not a valid
    # convergence test along constrained degrees of freedom.
    objective_rms = progress.last_objective_rms
    converged = (
        spec.max_iterations == 0
        or (
            objective_rms is not None
            and objective_rms <= spec.tolerance_kj_mol_nm * (1.0 + 1.0e-6)
        )
    )
    report = {
        "schema_version": 1,
        "status": "converged" if converged else "not_converged",
        "converged": converged,
        "inputs": {
            "coordinates": input_gro,
            "topology": topology_path,
            "include_dir": include_dir,
        },
        "outputs": {
            "gro": output_gro,
            "text_report": report_txt,
            "json_report": report_json,
        },
        "system": {"particles": topology_count},
        "positional_restraints": restraints,
        "energy_definition": "force-field potential plus protein/ligand and lipid head-group positional restraints",
        "minimizer_progress": {
            "reporter_callbacks": progress.callback_count,
            "constraint_restraint_stages": progress.restraint_stages,
            "last_iteration_in_stage": progress.last_iteration,
            "last_objective_rms_kj_mol_nm": objective_rms,
            "last_statistics": progress.last_statistics,
        },
        "openmm": {
            "version": getattr(openmm, "__version__", "unknown"),
            "platform": platform_name,
            "properties": properties,
            "failed_auto_platforms": platform_failures,
        },
        "settings": {
            "algorithm": spec.algorithm,
            "tolerance_kj_mol_nm": spec.tolerance_kj_mol_nm,
            "max_iterations": spec.max_iterations,
            "constraint_tolerance": spec.constraint_tolerance,
            "nonbonded_method": spec.nonbonded_method,
            "nonbonded_cutoff_nm": spec.nonbonded_cutoff_nm,
            "switch_distance_nm": spec.switch_distance_nm,
            "constraints": spec.constraints,
            "rigid_water": spec.rigid_water,
            "ewald_error_tolerance": spec.ewald_error_tolerance,
            "use_dispersion_correction": spec.use_dispersion_correction,
            "requested_platform": spec.platform,
            "requested_precision": spec.precision,
            "device_index": spec.device_index,
            "defines": spec.defines,
        },
        "before": {
            "potential_energy_kj_mol": initial_energy,
            "rms_force_kj_mol_nm": initial_rms,
            "max_force_kj_mol_nm": initial_max,
        },
        "after": {
            "potential_energy_kj_mol": final_energy,
            "rms_force_kj_mol_nm": final_rms,
            "max_force_kj_mol_nm": final_max,
        },
        "elapsed_seconds": time.monotonic() - start,
        "warnings": [],
    }
    if final_energy > initial_energy:
        report["warnings"].append(
            "The potential energy including positional restraints increased. OpenMM minimizes a combined "
            "potential-plus-constraint-restraint objective while satisfying constraints."
        )
    if not converged:
        report["warnings"].append(
            "The minimizer objective gradient remained above the requested "
            "tolerance. Increase max_iterations or inspect the system before "
            "using these coordinates."
        )

    destinations = (output_gro, report_txt, report_json)
    for destination in destinations:
        os.makedirs(os.path.dirname(destination), exist_ok=True)
    staged = []
    try:
        for destination in destinations:
            suffix = os.path.splitext(destination)[1]
            fd, temporary = tempfile.mkstemp(
                prefix=f".{os.path.basename(destination)}.",
                suffix=suffix,
                dir=os.path.dirname(destination),
            )
            os.close(fd)
            staged.append((temporary, destination))
        _write_minimized_gro(input_gro, staged[0][0], positions)
        Path(staged[1][0]).write_text(_render_text_report(report), encoding="utf-8")
        with open(staged[2][0], "w", encoding="utf-8") as handle:
            json.dump(report, handle, indent=2, sort_keys=True, allow_nan=False)
            handle.write("\n")
        _promote_staged_outputs(staged)
        staged.clear()
    except MinimizationError:
        raise
    except Exception as exc:
        raise MinimizationError(
            f"Could not write minimized outputs transactionally: {exc}"
        ) from exc
    finally:
        for temporary, _destination in staged:
            try:
                os.unlink(temporary)
            except FileNotFoundError:
                pass

    return report
