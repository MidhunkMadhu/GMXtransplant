"""Prepare portable minimization inputs. Do not import or launch OpenMM."""
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


RESOURCE_DEFAULTS = {
    "cpus_per_task": 8,
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
    if set(spec.resources) - RESOURCE_DEFAULTS.keys():
        raise BundleError("Unknown minimization.resources keys: " +
                          ", ".join(sorted(set(spec.resources) - RESOURCE_DEFAULTS.keys())))
    resources = RESOURCE_DEFAULTS | spec.resources
    if type(resources["cpus_per_task"]) is not int or resources["cpus_per_task"] < 1:
        raise BundleError("resources.cpus_per_task must be a positive integer")
    if not isinstance(resources["time"], str) or not re.fullmatch(
            r"(?:\d+-)?\d{2,3}:[0-5]\d:[0-5]\d", resources["time"]):
        raise BundleError("resources.time must be HH:MM:SS or D-HH:MM:SS")
    if not isinstance(resources["memory"], str) or not re.fullmatch(
            r"[1-9]\d*[KMGT]", resources["memory"]):
        raise BundleError("resources.memory must be a positive size, for example 8G")
    _defines(spec.defines, "minimization.defines")
    if spec.nonbonded_method not in {"PME", "CutoffPeriodic"}:
        raise BundleError("Portable OpenMM preparation supports PME or CutoffPeriodic")
    if spec.max_iterations == 0:
        raise BundleError("Prepared runs require a finite positive max_iterations")
    return resources


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
    not an OpenMM/force-field compatibility check; execution checks compatibility.
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


def _write(folder, name, text, executable=False):
    path = folder / name
    path.write_text(text, encoding="utf-8")
    if executable:
        path.chmod(0o755)


def _slurm(resources, gpu=False):
    text = f'''#!/usr/bin/env bash
#SBATCH --job-name=openmm_minimize
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={resources['cpus_per_task']}
#SBATCH --time={resources['time']}
#SBATCH --mem={resources['memory']}
#SBATCH --output=slurm-%j.out
'''
    if gpu:
        text += "#SBATCH --gpus=1\n"
    text += '''set -euo pipefail
# Submit from openmm_minimization. Activate your software environment here.
cd "${SLURM_SUBMIT_DIR:?Submit from the minimization folder}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENMM_CPU_THREADS="${SLURM_CPUS_PER_TASK}"
'''
    platform = "CUDA" if gpu else "CPU"
    return text + f"srun --ntasks=1 bash run.sh --platform {platform}\n"


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
    resources = validate_preparation(spec)
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
            folder = Path(tmp) / "bundle"
            folder.mkdir()
            coordinates = Path(spec.coordinates_path).read_bytes()
            (folder / "input.gro").write_bytes(coordinates)
            originals = _copy_topology(spec.topology_path, spec.include_dir, folder)
            # GromacsTopFile implicitly enables FLEXIBLE for water bonds/angles.
            audit = _audit(folder, {"FLEXIBLE": 1} | spec.defines, "openmm")
            if audit["position_restraints"]:
                raise BundleError("Deactivate topology position_restraints before preparation; "
                                  "the OpenMM runner adds protein/ligand heavy-atom restraints")
            if 2 * spec.nonbonded_cutoff_nm >= min(audit["box_angstrom_degrees"][:3]) / 10:
                raise BundleError("OpenMM cutoff must be below half the shortest box side")
            settings = asdict(replace(spec, enabled=False, coordinates_path="input.gro",
                topology_path="audit_openmm.top", include_dir=".",
                output_gro_path="results/minimized.gro", report_path="results/minimization_report",
                platform="CPU" if spec.platform == "auto" else spec.platform,
                defines={}, output_dir="unused"))
            _write(folder, "minimization.yaml", yaml.safe_dump({"minimization": settings}, sort_keys=False))
            support = folder / "runner_support"
            support.mkdir()
            for name in ("minimization.py", "minimization_restraints.py", "config.py",
                         "output.py", "chainids.py", "classify.py"):
                shutil.copyfile(Path(__file__).parent / name, support / name)
            _write(folder, "minimize_openmm.py", OPENMM_RUNNER, True)
            _write(folder, "requirements.txt", "numpy>=1.23\nPyYAML>=6\nMDAnalysis>=2.4\nopenmm>=8.2\n")
            _write(folder, "run.sh", '''#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
export OPENMM_CPU_THREADS="${SLURM_CPUS_PER_TASK:-${OPENMM_CPU_THREADS:-''' + str(resources["cpus_per_task"]) + '''}}"
exec python3 minimize_openmm.py "$@"
''', True)
            _write(folder, "submit_gpu.slurm", _slurm(resources, True), True)
            _write(folder, "submit_cpu.slurm", _slurm(resources), True)
            _write(folder, "README.md", BUNDLE_README)
            manifest = {"schema_version": 2, "gmxtransplant_version": __version__,
                        "status": "prepared_not_run", "assembly_mode": assembly_mode,
                        "engine": "openmm", "resources": resources, "audit": audit,
                        "settings": settings, "source_sha256": originals,
                        "coordinate_source_sha256": hashlib.sha256(coordinates).hexdigest()}
            manifest["files_sha256"] = {
                p.relative_to(folder).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
                for p in sorted(folder.rglob("*")) if p.is_file()}
            _write(folder, "manifest.json", json.dumps(manifest, indent=2, allow_nan=False) + "\n")
            check_destination(destination)
            os.rename(folder, destination)
        return {"output_dir": str(destination), "engine": "openmm", "status": "prepared_not_run"}
    finally:
        lock.unlink()


BUNDLE_README = """# OpenMM minimization

This folder contains portable inputs for restrained clash relaxation. Preparation
does not run or submit a job. Copy the whole folder to your preferred machine.

## Run

Activate your Python environment and install the runner dependencies:

```bash
python3 -m pip install -r requirements.txt
bash run.sh
```

The default is a multicore CPU run. Set OPENMM_CPU_THREADS to change the local
thread count. For a GPU, use `bash run.sh --platform CUDA` (or HIP/OpenCL).
Requested platforms fail explicitly if unavailable.

For SLURM, edit the generic resource directives and software setup, then submit
from this folder with `sbatch submit_cpu.slurm` or `sbatch submit_gpu.slurm`.
The CPU template uses one process with multiple threads; the GPU template requests
one GPU. Neither template contains a cluster name or account.

## Settings and restraints

Review minimization.yaml before running. Harmonic positional restraints keep
protein and ligand heavy atoms near input.gro, with a default force constant
of 1000 kJ/mol/nm². Hydrogens, water, ions and recognized membrane lipids are
unrestrained. Unrecognized molecules are treated as ligands; use
restraint_residue_classes to classify custom lipids or other unusual residues.

The runner uses L-BFGS, a finite iteration limit, and an RMS objective-gradient
tolerance. Review nonbonded settings for your force field. OpenMM potential
switching is not the same as LJ force switching in the source force field.

audit_openmm.top is the expanded topology used for execution; topol.top and
toppar/ preserve its source dependencies. Regenerate the folder if changing
topology defines. Active topology position_restraints are rejected because the
runner supplies its own heavy-atom restraints.

## Results

results/minimized.gro contains the relaxed coordinates.
results/minimization_report.json and .txt record convergence, energies, forces,
restraint selection and heavy-atom displacement. Energies include positional
restraints. The runner exits with code 7 on failure or unconfirmed convergence.

Existing results are never overwritten. Rename a previous results/ directory
before rerunning. Inspect the structure and convergence report before continuing;
this clash-relaxation step is not equilibration.
"""
