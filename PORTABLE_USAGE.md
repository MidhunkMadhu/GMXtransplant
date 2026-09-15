# Portable usage

Install from the source directory:

```bash
python3 -m pip install .
```

Create and edit a configuration in your job directory:

```bash
gmxtransplant --show-example protein > protein_replace.yaml
gmxtransplant --mode protein -i protein_replace.yaml --dry-run
gmxtransplant --mode protein -i protein_replace.yaml --prepare-minimization
```

Use modes `lig` or `chl` for ligand replacement or cholesterol restoration.
Omit the preparation flag and leave `minimization.enabled: false` for assembly alone.

To prepare an existing matching system:

```bash
gmxtransplant prepare-minimization --coordinates step5_input.gro --topology topol.top
```

Both workflows create one `openmm_minimization/` folder. Preparation copies and
audits the inputs; it does not minimize or submit a job. The folder includes all
topology dependencies, the standalone runner, settings, and generic CPU/GPU
SLURM templates.

Move the whole folder to your execution machine, activate your Python environment,
and review the settings:

```bash
cd openmm_minimization
python3 -m pip install -r requirements.txt
bash run.sh
```

The runner relaxes clashes with harmonic restraints on protein and ligand heavy
atoms. Recognized membrane lipids, water, ions, and hydrogens remain unrestrained.
Classify custom lipids with `minimization.restraint_residue_classes` so they are
not treated as ligands.

For SLURM, edit and submit either `submit_cpu.slurm` or `submit_gpu.slurm` from
this folder. Inspect `results/minimized.gro` and the convergence reports before
continuing. Existing inputs and results are never overwritten.

See [README.md](README.md) for all modes, settings, and limitations.
