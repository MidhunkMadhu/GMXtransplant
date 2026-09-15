# Portable usage

Install from the source directory:

```bash
python3 -m pip install .
gmxtransplant --help
```

Create and edit an example in a separate job directory:

```bash
gmxtransplant --show-example receptor > receptor_replace.yaml
gmxtransplant --mode receptor -i receptor_replace.yaml --dry-run
gmxtransplant --mode receptor -i receptor_replace.yaml --prepare-minimization
```

Modes `lig` and `chl` support the same preparation option. Nothing is minimized
or submitted during assembly. Version 0.3 removes `--minimize`; YAML
`minimization.enabled: true` now requests preparation only.

To prepare an existing system:

```bash
gmxtransplant prepare-minimization --coordinates step5_input.gro --topology topol.top
```

The new `minimization_inputs/` directory contains independent GROMACS and OpenMM
folders with copied GRO/topology inputs, complete topology dependencies, runners,
generic multicore SLURM templates and instructions. An existing destination is
never overwritten. No OpenMM or GROMACS installation is needed for preparation.

Move either engine folder to the execution machine, activate its software
environment, review settings, and run `bash run.sh` or `sbatch submit_cpu.slurm`
from that folder. OpenMM also has `submit_gpu.slurm` and `requirements.txt`.
Read the generated README before execution. All runtime results go under a new
`results/` directory. The full configuration and limitations are in [README.md](README.md).
