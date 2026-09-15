"""Prepare portable minimization inputs. Never import or launch either engine."""
from __future__ import annotations

from dataclasses import asdict, replace
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import tempfile

import MDAnalysis as mda
import numpy as np
import yaml

from config import ConfigError, MinimizationSpec
from itp import parse_itp
from topology import parse_top_molecules
from gmxtransplant import __version__


class BundleError(ValueError):
    pass


GROMACS_DEFAULTS = {
    "integrator": "steep", "emtol": 100.0, "emstep": 0.01, "nsteps": 5000,
    "constraints": "h-bonds", "coulombtype": "PME", "rcoulomb": 1.2,
    "vdwtype": "Cut-off", "vdw_modifier": "Force-switch",
    "rvdw_switch": 1.0, "rvdw": 1.2, "disp_corr": "no",
    "ewald_rtol": 5e-4, "defines": {}, "restrained_defines": {},
}
RESOURCE_DEFAULTS = {
    "cpus_per_task": 8, "gromacs_mpi_tasks": 1,
    "time": "00:30:00", "memory": "8G",
}
INCLUDE = re.compile(r'^\s*#\s*include\s+["<]([^">]+)[">]\s*(?:;.*|//.*)?$')
IDENTIFIER = re.compile(r"^[A-Za-z_][A-Za-z_0-9]*$")


def _defines(values, label):
    if not isinstance(values, dict):
        raise BundleError(f"{label} must be a mapping")
    for name, value in values.items():
        if not isinstance(name, str) or not IDENTIFIER.fullmatch(name):
            raise BundleError(f"Invalid preprocessor name in {label}: {name!r}")
        if not isinstance(value, (str, int, float)) or isinstance(value, bool):
            raise BundleError(f"{label}.{name} must be a string or number")
        if not re.fullmatch(r"[A-Za-z0-9_.+\-]+", str(value)):
            raise BundleError(f"Unsafe preprocessor value in {label}.{name}")


def validate_preparation(spec: MinimizationSpec):
    """Validate without reading molecular input files or loading an engine."""
    spec.validate(check_paths=False, require_inputs=True)
    if set(spec.gromacs) - GROMACS_DEFAULTS.keys():
        raise BundleError("Unknown minimization.gromacs keys: " +
                          ", ".join(sorted(set(spec.gromacs) - GROMACS_DEFAULTS.keys())))
    if set(spec.resources) - RESOURCE_DEFAULTS.keys():
        raise BundleError("Unknown minimization.resources keys: " +
                          ", ".join(sorted(set(spec.resources) - RESOURCE_DEFAULTS.keys())))
    gmx = GROMACS_DEFAULTS | spec.gromacs
    resources = RESOURCE_DEFAULTS | spec.resources
    for name in ("cpus_per_task", "gromacs_mpi_tasks"):
        if type(resources[name]) is not int or resources[name] < 1:
            raise BundleError(f"resources.{name} must be a positive integer")
    if not isinstance(resources["time"], str) or not re.fullmatch(
            r"(?:\d+-)?\d{2,3}:[0-5]\d:[0-5]\d", resources["time"]):
        raise BundleError("resources.time must be HH:MM:SS or D-HH:MM:SS")
    if not isinstance(resources["memory"], str) or not re.fullmatch(
            r"[1-9]\d*[KMGT]", resources["memory"]):
        raise BundleError("resources.memory must be a positive size, for example 8G")
    choices = {"integrator": {"steep"},
               "constraints": {"none", "h-bonds", "all-bonds", "h-angles"},
               "coulombtype": {"PME", "Cut-off"}, "vdwtype": {"Cut-off"},
               "vdw_modifier": {"Force-switch", "Potential-switch", "Potential-shift", "None"},
               "disp_corr": {"no", "EnerPres", "Ener"}}
    for name, allowed in choices.items():
        if not isinstance(gmx[name], str) or gmx[name] not in allowed:
            raise BundleError(f"gromacs.{name} must be one of {sorted(allowed)}")
    for name in ("emtol", "emstep", "rcoulomb", "rvdw", "ewald_rtol"):
        if (type(gmx[name]) not in (int, float)
                or not np.isfinite(gmx[name]) or gmx[name] <= 0):
            raise BundleError(f"gromacs.{name} must be finite and positive")
    if type(gmx["nsteps"]) is not int or gmx["nsteps"] < 1:
        raise BundleError("gromacs.nsteps must be a positive integer")
    switch = gmx["rvdw_switch"]
    if (type(switch) not in (int, float) or not np.isfinite(switch)
            or not 0 <= switch < gmx["rvdw"]):
        raise BundleError("gromacs.rvdw_switch must be >= 0 and below rvdw")
    if gmx["ewald_rtol"] >= 1:
        raise BundleError("gromacs.ewald_rtol must be less than 1")
    _defines(spec.defines, "minimization.defines")
    _defines(gmx["defines"], "gromacs.defines")
    _defines(gmx["restrained_defines"], "gromacs.restrained_defines")
    if gmx["restrained_defines"] and gmx["defines"]:
        raise BundleError("Two-stage GROMACS preparation requires empty final-stage defines")
    if "openmm" in spec.engines and spec.nonbonded_method not in {"PME", "CutoffPeriodic"}:
        raise BundleError("Portable OpenMM preparation currently supports PME or CutoffPeriodic")
    if spec.max_iterations == 0:
        raise BundleError("Prepared runs require a finite positive max_iterations")
    return gmx, resources


def check_destination(destination):
    destination = Path(destination).absolute()
    if destination.exists() or destination.is_symlink():
        raise BundleError(f"Refusing to overwrite existing bundle: {destination}")
    return destination


def _copy_topology(topology, include_dir, destination):
    """Copy all branches of a literal include graph, rewriting local paths.

    Unique dependency names prevent basename collisions. Inactive dependencies
    must also resolve, so users can safely change defines after copying a folder.
    """
    root = Path(topology).resolve(strict=True)
    include_root = Path(include_dir).resolve() if include_dir else root.parent
    mapping, visiting, originals = {}, set(), {}

    def visit(source):
        source = source.resolve(strict=True)
        if source in visiting:
            raise BundleError(f"Cyclic topology includes: {source.name}; preprocess first")
        if source in mapping:
            return mapping[source]
        relative = (Path("topol.top") if not mapping else
                    Path("toppar") / f"{len(mapping):04d}_{source.name}")
        mapping[source] = relative
        visiting.add(source)
        raw = source.read_text(encoding="utf-8")
        originals[relative.as_posix()] = hashlib.sha256(source.read_bytes()).hexdigest()
        lines = []
        for line in raw.splitlines():
            if re.match(r"^\s*#\s*include\b", line):
                match = INCLUDE.fullmatch(line)
                if match is None:
                    raise BundleError(f"Unsupported include in {source.name}: {line}")
                target = match.group(1)
                candidates = [source.parent / target, root.parent / target, include_root / target]
                found = next((p for p in candidates if p.is_file()), None)
                if found is None:
                    raise BundleError(f"Unresolved include {target!r} in {source.name}")
                child = visit(found)
                local = os.path.relpath(child, relative.parent).replace(os.sep, "/")
                line = f'#include "{local}"'
            lines.append(line)
        out = destination / relative
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text("\n".join(lines) + "\n", encoding="utf-8")
        visiting.remove(source)
        return relative

    visit(root)
    (destination / "toppar").mkdir(exist_ok=True)
    return originals


def _preprocess(path, defines):
    """Strict subset of the GROMACS preprocessor for atom-identity validation.

    Unsupported constructs are rejected, never silently ignored. This audit is
    not an engine/force-field compatibility check; grompp remains mandatory.
    """
    macros = {k: str(v) for k, v in defines.items()}

    def expand(file, stack):
        file = file.resolve()
        if file in stack:
            raise BundleError("Cyclic include during topology audit")
        active, conditions, output = True, [], []
        for raw in file.read_text(encoding="utf-8").splitlines():
            line = raw.split(";", 1)[0].strip()
            if not line:
                continue
            if line.startswith("#"):
                parts = line[1:].strip().split(None, 1)
                directive, argument = parts[0], parts[1] if len(parts) > 1 else ""
                if directive in ("ifdef", "ifndef"):
                    if not IDENTIFIER.fullmatch(argument):
                        raise BundleError(f"Malformed #{directive} in {file.name}")
                    condition = argument in macros
                    if directive == "ifndef":
                        condition = not condition
                    conditions.append((active, condition, False))
                    active = active and condition
                elif directive == "else":
                    if not conditions or conditions[-1][2]:
                        raise BundleError(f"Unmatched/repeated #else in {file.name}")
                    parent, condition, _ = conditions[-1]
                    active = parent and not condition
                    conditions[-1] = parent, condition, True
                elif directive == "endif":
                    if not conditions:
                        raise BundleError(f"Unmatched #endif in {file.name}")
                    active = conditions.pop()[0]
                elif directive in ("if", "elif"):
                    raise BundleError(f"Unsupported #{directive}; supply a preprocessed topology")
                elif active:
                    if directive == "include":
                        match = INCLUDE.fullmatch(line)
                        if not match:
                            raise BundleError("Only literal includes are supported")
                        output.extend(expand(file.parent / match.group(1), stack + [file]))
                    elif directive == "define":
                        pair = argument.split(None, 1)
                        if not pair or not IDENTIFIER.fullmatch(pair[0]):
                            raise BundleError("Function-like macros are unsupported; preprocess first")
                        macros[pair[0]] = pair[1] if len(pair) > 1 else ""
                    elif directive == "undef":
                        macros.pop(argument, None)
                    else:
                        raise BundleError(f"Unsupported #{directive}; preprocess topology first")
                continue
            if active:
                # Substitute complete tokens, not substrings of atom names.
                output.append(" ".join(macros.get(token, token) for token in line.split()))
        if conditions:
            raise BundleError(f"Unclosed conditional in {file.name}")
        return output

    return "\n".join(expand(Path(path), [])) + "\n"


def _audit(folder, defines, label):
    expanded = _preprocess(folder / "topol.top", defines)
    audit_path = folder / f"audit_{label}.top"
    audit_path.write_text(expanded, encoding="utf-8")
    definitions = parse_itp(str(audit_path))
    molecules = parse_top_molecules(str(audit_path))
    universe = mda.Universe(str(folder / "input.gro"))
    if not np.all(np.isfinite(universe.atoms.positions)):
        raise BundleError("Coordinates contain NaN or infinity")
    box = universe.dimensions
    if (box is None or not np.all(np.isfinite(box)) or np.any(box[:3] <= 0)
            or not np.allclose(box[3:], 90, atol=0.001)):
        raise BundleError("Preparation requires a finite orthorhombic periodic box")
    offset, charge = 0, 0.0
    for name, count in molecules:
        if name not in definitions or count < 0:
            raise BundleError(f"Missing definition or invalid count for {name}")
        definition = definitions[name]
        for _ in range(count):
            for atom in definition.atoms:
                if offset >= len(universe.atoms):
                    raise BundleError("Topology has more atoms than coordinates")
                observed = universe.atoms[offset]
                if (str(observed.resname), str(observed.name)) != (atom.resname, atom.atom_name):
                    raise BundleError(f"Coordinate/topology identity mismatch at atom {offset + 1}")
                offset += 1
        charge += count * definition.charge
    if offset != len(universe.atoms):
        raise BundleError("Coordinate/topology atom counts differ")
    sections = re.findall(r"\[\s*([^]]+?)\s*\]", expanded)
    has_restraints = "position_restraints" in [s.strip() for s in sections]
    return {"atom_count": offset, "net_charge": charge,
            "box_angstrom_degrees": box.tolist(), "position_restraints": has_restraints,
            "defines": defines, "preprocessed_topology": audit_path.name}


def _mdp(settings, defines):
    values = {key.replace("_", "-"): value for key, value in settings.items()
              if key not in ("defines", "restrained_defines")}
    values["DispCorr"] = values.pop("disp-corr")
    values.update({"cutoff-scheme": "Verlet", "nstlist": 10,
                   "pbc": "xyz", "tcoupl": "no", "pcoupl": "no",
                   "define": " ".join(f"-D{k}={v}" for k, v in defines.items())})
    return ("; Review for the source force field before running.\n"
            "; emtol is maximum force, not OpenMM RMS objective tolerance.\n" +
            "\n".join(f"{k} = {v}" for k, v in values.items()) + "\n")


def _write(folder, name, text, executable=False):
    path = folder / name
    path.write_text(text, encoding="utf-8")
    if executable:
        path.chmod(0o755)


def _slurm(engine, resources, gpu=False):
    tasks = resources["gromacs_mpi_tasks"] if engine == "gromacs" else 1
    text = f'''#!/usr/bin/env bash
#SBATCH --job-name={engine}_minimize
#SBATCH --nodes=1
#SBATCH --ntasks={tasks}
#SBATCH --cpus-per-task={resources['cpus_per_task']}
#SBATCH --time={resources['time']}
#SBATCH --mem={resources['memory']}
#SBATCH --output=slurm-%j.out
'''
    if gpu:
        text += "#SBATCH --gpus=1\n"
    text += '''set -euo pipefail
# Submit from this engine directory. Load your software environment here.
cd "${SLURM_SUBMIT_DIR:?Submit with sbatch from the engine folder}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
'''
    if engine == "gromacs":
        text += 'export GMX="${GMX:-gmx_mpi}"\nexport GMX_MPI=1\nbash run.sh\n'
    else:
        platform = "CUDA" if gpu else "CPU"
        text += ('export OPENMM_CPU_THREADS="${SLURM_CPUS_PER_TASK}"\n'
                 f'srun --ntasks=1 bash run.sh --platform {platform}\n')
    return text


def _gromacs_runner(two_stage, cpus):
    stages = "restrained minimization" if two_stage else "minimization"
    return r'''#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
GMX="${GMX:-gmx}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-${OMP_NUM_THREADS:-__DEFAULT_THREADS__}}"
[[ "$OMP_NUM_THREADS" =~ ^[1-9][0-9]*$ ]] || { echo "Invalid thread count" >&2; exit 2; }
command -v "$GMX" >/dev/null
command -v python3 >/dev/null
# mkdir is intentionally exclusive: never overwrite a previous run.
mkdir results
"$GMX" --version > results/engine_version.txt 2>&1
input=input.gro
for stage in STAGES; do
    "$GMX" grompp -f "${stage}.mdp" -c "$input" -r input.gro \
        -p topol.top -o "results/${stage}.tpr" \
        -po "results/${stage}_resolved.mdp" -pp "results/${stage}_processed.top" \
        -maxwarn 0 > "results/${stage}_grompp.log" 2>&1
    if [[ "${GMX_MPI:-0}" == 1 ]]; then
        srun --ntasks="${SLURM_NTASKS:-1}" "$GMX" mdrun \
            -deffnm "results/${stage}" -ntomp "$OMP_NUM_THREADS" \
            > "results/${stage}_console.log" 2>&1
    else
        "$GMX" mdrun -deffnm "results/${stage}" -ntmpi 1 \
            -ntomp "$OMP_NUM_THREADS" > "results/${stage}_console.log" 2>&1
    fi
    python3 summarize_gromacs.py "results/${stage}"
    input="results/${stage}.gro"
done
'''.replace("__DEFAULT_THREADS__", str(cpus)).replace("STAGES", stages)


GROMACS_SUMMARY = r'''"""Report actual convergence; a zero mdrun exit code is insufficient."""
import json
from pathlib import Path
import re
import sys

prefix = Path(sys.argv[1])
log = prefix.with_suffix(".log").read_text()
converged = bool(re.search(r"converged to Fmax\s*<", log))
summary = {"converged": converged, "status": "converged" if converged else "not_converged",
           "coordinate_file": str(prefix.with_suffix(".gro")), "log": str(prefix.with_suffix(".log"))}
for key, label in (("potential_energy_kj_mol", "Potential Energy"),
                   ("max_force_kj_mol_nm", "Maximum force")):
    matches = re.findall(label + r"\s*=\s*([+\-0-9.eE]+)", log)
    if matches:
        summary[key] = float(matches[-1])
prefix.with_suffix(".json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
print(json.dumps(summary))
sys.exit(0 if converged else 7)
'''


OPENMM_RUNNER = '''#!/usr/bin/env python3
"""Standalone OpenMM runner; no installed GMXtransplant package required."""
import argparse
from dataclasses import replace
import json
from pathlib import Path
import sys
import os

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "runner_support"))
from config import load_minimization_config
from minimization import run_minimization

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--platform", choices=["CUDA", "HIP", "OpenCL", "CPU", "Reference"])
    args = parser.parse_args()
    os.chdir(ROOT)
    spec = load_minimization_config("minimization.yaml")
    if args.platform:
        spec = replace(spec, platform=args.platform)
    # No silent Reference/CPU fallback when GPU execution was requested.
    Path("results").mkdir(exist_ok=False)
    try:
        result = run_minimization(spec)
    except Exception as exc:
        Path("results/failure.json").write_text(json.dumps({"status": "failed", "error": str(exc)}) + "\\n")
        print(str(exc), file=sys.stderr)
        return 7
    print(json.dumps(result, indent=2))
    return 0 if result["converged"] else 7

if __name__ == "__main__":
    sys.exit(main())
'''


def prepare_minimization(spec: MinimizationSpec, *, assembly_mode=None):
    gmx, resources = validate_preparation(spec)
    spec.validate(check_paths=True, require_inputs=True)
    destination = check_destination(spec.output_dir)
    destination.parent.mkdir(parents=True, exist_ok=True)
    lock = destination.parent / f".{destination.name}.lock"
    try:
        fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError as exc:
        raise BundleError(f"Another preparation may be using {destination}") from exc
    os.close(fd)
    try:
        check_destination(destination)
        with tempfile.TemporaryDirectory(prefix=f".{destination.name}.", dir=destination.parent) as tmp:
            staging = Path(tmp) / "bundle"
            staging.mkdir()
            manifest = {"schema_version": 1, "gmxtransplant_version": __version__,
                        "status": "prepared_not_run", "assembly_mode": assembly_mode,
                        "engines": {}, "resources": resources,
                        "notes": ["Input identity audit is not an engine compatibility test.",
                                  "GROMACS emtol uses maximum force; OpenMM tolerance uses RMS objective gradient.",
                                  "Force-switch and OpenMM potential switching are not numerically equivalent."]}
            coordinates = Path(spec.coordinates_path).read_bytes()
            for engine in spec.engines:
                folder = staging / engine
                folder.mkdir()
                (folder / "input.gro").write_bytes(coordinates)
                originals = _copy_topology(spec.topology_path, spec.include_dir, folder)
                engine_data = {"source_sha256": originals, "audits": {}}
                manifest["engines"][engine] = engine_data
                if engine == "gromacs":
                    stages = [("minimization", gmx["defines"])]
                    if gmx["restrained_defines"]:
                        stages.insert(0, ("restrained", gmx["restrained_defines"]))
                    for stage, defines in stages:
                        audit = _audit(folder, defines, stage)
                        if stage == "restrained" and not audit["position_restraints"]:
                            raise BundleError("Restrained stage requested but no active position_restraints found")
                        if stage == "minimization" and len(stages) > 1 and audit["position_restraints"]:
                            raise BundleError("Final unrestrained stage still has active position restraints")
                        engine_data["audits"][stage] = audit
                        _write(folder, f"{stage}.mdp", _mdp(gmx, defines))
                    shortest_nm = min(audit["box_angstrom_degrees"][:3]) / 10
                    if 2 * max(gmx["rcoulomb"], gmx["rvdw"]) >= shortest_nm:
                        raise BundleError("GROMACS cutoff must be below half the shortest box side")
                    engine_data["settings"] = gmx
                    _write(folder, "run.sh", _gromacs_runner(len(stages) > 1, resources["cpus_per_task"]), True)
                    _write(folder, "summarize_gromacs.py", GROMACS_SUMMARY)
                else:
                    # GromacsTopFile implicitly enables FLEXIBLE to obtain water
                    # bonds/angles before createSystem applies rigidWater.
                    audit = _audit(folder, {"FLEXIBLE": 1} | spec.defines, "openmm")
                    if audit["position_restraints"]:
                        raise BundleError("OpenMM export does not translate GROMACS position_restraints; use GROMACS for this protocol")
                    if 2 * spec.nonbonded_cutoff_nm >= min(audit["box_angstrom_degrees"][:3]) / 10:
                        raise BundleError("OpenMM cutoff must be below half the shortest box side")
                    engine_data["audits"]["openmm"] = audit
                    # Use expanded topology so the audited active definitions are exactly those executed.
                    settings = asdict(replace(spec, enabled=False, coordinates_path="input.gro",
                        topology_path="audit_openmm.top", include_dir=".",
                        output_gro_path="results/minimized.gro", report_path="results/minimization_report",
                        platform="CPU" if spec.platform == "auto" else spec.platform,
                        defines={}, output_dir="unused"))
                    engine_data["settings"] = settings
                    _write(folder, "minimization.yaml", yaml.safe_dump({"minimization": settings}, sort_keys=False))
                    support = folder / "runner_support"
                    support.mkdir()
                    for name in ("minimization.py", "config.py", "output.py", "chainids.py", "classify.py"):
                        shutil.copyfile(Path(__file__).parent / name, support / name)
                    _write(folder, "minimize_openmm.py", OPENMM_RUNNER, True)
                    _write(folder, "requirements.txt", "numpy>=1.23\nPyYAML>=6\nMDAnalysis>=2.4\nopenmm>=8.2\n")
                    _write(folder, "run.sh", '''#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
export OPENMM_CPU_THREADS="${SLURM_CPUS_PER_TASK:-${OPENMM_CPU_THREADS:-''' + str(resources["cpus_per_task"]) + '''}}"
exec python3 minimize_openmm.py "$@"
''', True)
                    _write(folder, "submit_gpu.slurm", _slurm(engine, resources, True), True)
                _write(folder, "submit_cpu.slurm", _slurm(engine, resources), True)
                _write(folder, "README.md", _engine_readme(engine, bool(gmx["restrained_defines"])))
            manifest["coordinate_source_sha256"] = hashlib.sha256(coordinates).hexdigest()
            _write(staging, "README.md", "# Minimization inputs\n\nPrepared only; nothing has been run or submitted.\n"
                   "Each engine folder is independent and uses the same input coordinates.\n"
                   "Review its README and scientific settings before running.\n"
                   "Checksums and input audits are in manifest.json. No original files are overwritten.\n")
            manifest["files_sha256"] = {p.relative_to(staging).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
                                        for p in sorted(staging.rglob("*")) if p.is_file()}
            _write(staging, "manifest.json", json.dumps(manifest, indent=2, allow_nan=False) + "\n")
            check_destination(destination)
            os.rename(staging, destination)
        return {"output_dir": str(destination), "engines": spec.engines, "status": "prepared_not_run"}
    finally:
        lock.unlink()


def _engine_readme(engine, restrained):
    common = ("\nLoad/activate your software environment before running. Submit from this folder; "
              "the generic SLURM templates intentionally contain no site-specific configuration.\n\n"
              "Local: `bash run.sh`. Scheduled CPU: `sbatch submit_cpu.slurm`. "
              "Edit resource directives as needed. Scripts create an exclusive results/ directory. "
              "For a rerun, preserve the previous results by renaming that directory yourself.\n\n"
              "Coordinate atom identity and charge were audited; engine parsing and force-field "
              "compatibility are only checked at execution. Inspect convergence and the final structure. "
              "Minimization is not equilibration.\n")
    if engine == "gromacs":
        return ("# GROMACS minimization\n\nRequires GROMACS, Bash, and Python 3 (standard library only). "
                "Local execution uses thread-MPI gmx with one rank and multiple OpenMP threads. "
                "The CPU SLURM template uses external-MPI gmx_mpi and srun. "
                "Override GMX with an executable path, not a shell command; do not put spaces/arguments in it.\n\n"
                "Review minimization.mdp: the default is steepest descent with PME and LJ force switching, "
                "not a universal force-field protocol. GROMACS applies explicit topology constraints/SETTLE "
                "even when constraints=none; review source water definitions. grompp uses -maxwarn 0. "
                "No TPR is generated until run.sh executes with your installed GROMACS.\n\n"
                + ("Two stages: restrained.mdp then minimization.mdp. The original input.gro is the "
                   "restraint reference. A nonconverged stage stops the sequence.\n" if restrained else "")
                + "Outputs: results/minimization.gro, .log, .edr, .tpr, and .json summary, plus preprocessing logs. "
                "Exit 7 means convergence was not confirmed; coordinates may still exist.\n" + common)
    return ("# OpenMM minimization\n\nInstall with `python3 -m pip install -r requirements.txt`. "
            "The bundled runner_support/ is a source snapshot; no GMXtransplant installation is required.\n\n"
            "GPU: `sbatch submit_gpu.slurm` or `bash run.sh --platform CUDA`. "
            "The GPU template requests one GPU. HIP/OpenCL can be selected explicitly in a suitably edited template. "
            "Device visibility is left to the scheduler. CPU threads follow OPENMM_CPU_THREADS. "
            "Requested platforms fail explicitly rather than silently falling back to Reference.\n\n"
            "The YAML uses audit_openmm.top, expanded with the requested defines during preparation. "
            "The audit includes OpenMM's implicit FLEXIBLE definition for water bonds/angles. "
            "Regenerate the bundle if changing topology defines. The original include tree is also provided. "
            "No GROMACS position restraints are silently translated. OpenMM potential switching differs from "
            "GROMACS force switching; adjust both engine protocols deliberately. RMS objective tolerance is "
            "not GROMACS maximum-force emtol.\n\n"
            "Outputs: results/minimized.gro and minimization_report.{txt,json}. "
            "The report includes platform, energy, force, elapsed time and convergence. "
            "Exit 7 means failure or unconfirmed convergence.\n" + common)
