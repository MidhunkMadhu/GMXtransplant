"""Write a Slurm batch script that runs an equilibration folder's sample_run.sh.

    gmxtransplant slurm-script popc_output/equilibration --profile dardel --account naiss2025-3-21
    sbatch popc_output/equilibration/submit.sh

A profile fills in the cluster's partition, modules and GROMACS commands; every
value can be given or overridden on the command line. Nothing is submitted.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import os
import shlex

PROFILES = {
    "generic": {"partition": None, "modules": [], "gmx": "gmx", "mdrun": "gmx mdrun", "nodes": 1,
                "ntasks_per_node": 1, "cpus_per_task": 32, "time": "08:00:00", "gpus": 0},
    "dardel": {"partition": "main", "modules": ["PDC/26.03", "gromacs/2026.1-cpeGNU-26.03"],
               "gmx": "srun -n 1 gmx_mpi", "mdrun": "srun gmx_mpi mdrun -pin on", "nodes": 1,
               "ntasks_per_node": 128, "cpus_per_task": 1, "time": "04:00:00", "gpus": 0},
    "gpu": {"partition": None, "modules": [], "gmx": "gmx", "mdrun": "gmx mdrun -nb gpu -pme gpu -update gpu",
            "nodes": 1, "ntasks_per_node": 1, "cpus_per_task": 16, "time": "08:00:00", "gpus": 1},
}


def script(folder, profile="generic", account=None, partition=None, modules=None, gmx=None, mdrun=None,
           nodes=None, ntasks_per_node=None, cpus_per_task=None, time=None, gpus=None, job_name=None):
    values = dict(PROFILES[profile])
    for key, value in (("partition", partition), ("modules", modules), ("gmx", gmx), ("mdrun", mdrun),
                       ("nodes", nodes), ("ntasks_per_node", ntasks_per_node), ("cpus_per_task", cpus_per_task),
                       ("time", time), ("gpus", gpus)):
        if value is not None:
            values[key] = value
    folder = Path(folder).resolve()
    if not (folder / "sample_run.sh").is_file():
        raise SystemExit(f"{folder} has no sample_run.sh: point at an equilibration folder written by GMXtransplant")
    name = job_name or f"gmxt-{folder.parent.name}"[:40]
    lines = ["#!/bin/bash", f"#SBATCH -J {name}"]
    if account:
        lines.append(f"#SBATCH -A {account}")
    if values["partition"]:
        lines.append(f"#SBATCH -p {values['partition']}")
    lines += [f"#SBATCH --nodes={values['nodes']}", f"#SBATCH --ntasks-per-node={values['ntasks_per_node']}",
              f"#SBATCH --cpus-per-task={values['cpus_per_task']}", f"#SBATCH -t {values['time']}",
              f"#SBATCH -o {name}.%j.out", f"#SBATCH -e {name}.%j.err"]
    if values["gpus"]:
        lines.append(f"#SBATCH --gpus={values['gpus']}")
    lines += ["# Written by GMXtransplant (gmxtransplant slurm-script). Submit with: sbatch submit.sh", "set -e",
              f"cd {shlex.quote(str(folder))}"]
    if values["modules"]:
        lines.append("module load " + " ".join(values["modules"]))
    lines += [f"export OMP_NUM_THREADS={values['cpus_per_task']}",
              f"export GMX={shlex.quote(values['gmx'])}", f"export MDRUN={shlex.quote(values['mdrun'])}",
              "./sample_run.sh",
              "if command -v python3 >/dev/null; then python3 analyze_equilibration.py | tee analysis.txt; fi"]
    target = folder / "submit.sh"
    target.write_text("\n".join(lines) + "\n")
    os.chmod(target, 0o755)
    return target


def main(argv=None):
    parser = argparse.ArgumentParser(prog="gmxtransplant slurm-script", description=__doc__.split("\n\n")[0])
    parser.add_argument("folder", help="an equilibration folder (with sample_run.sh)")
    parser.add_argument("--profile", choices=sorted(PROFILES), default="generic")
    parser.add_argument("--account", "-A")
    parser.add_argument("--partition", "-p")
    parser.add_argument("--module", dest="modules", action="append", help="module to load (repeat)")
    parser.add_argument("--gmx", help="command for grompp, e.g. 'srun -n 1 gmx_mpi'")
    parser.add_argument("--mdrun", help="command for mdrun, e.g. 'srun gmx_mpi mdrun -pin on'")
    parser.add_argument("--nodes", type=int)
    parser.add_argument("--ntasks-per-node", type=int)
    parser.add_argument("--cpus-per-task", type=int)
    parser.add_argument("--time", "-t")
    parser.add_argument("--gpus", type=int)
    parser.add_argument("--job-name", "-J")
    args = parser.parse_args(argv)
    target = script(args.folder, args.profile, args.account, args.partition, args.modules, args.gmx, args.mdrun,
                    args.nodes, args.ntasks_per_node, args.cpus_per_task, args.time, args.gpus, args.job_name)
    print(f"Wrote {target}; submit it with: sbatch {target}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
