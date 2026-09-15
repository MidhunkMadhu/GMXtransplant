# GMXtransplant

GMXtransplant builds a coordinate/topology bundle by inserting molecular components into an existing prepared environment. Its main use is placing a receptor, a receptor-ligand complex, a replacement ligand, or experimentally resolved cholesterol into an equilibrated membrane system while preserving as much of that environment as possible.

The pipeline removes complete environment molecules that clash with the inserted atoms, corrects charge when requested, regenerates the GROMACS molecule list, writes an index, and produces detailed text and JSON reports. It can also prepare independent GROMACS and OpenMM minimization folders containing copied inputs, running scripts, and generic multicore SLURM submission templates.

This is a structure-assembly tool and expects every molecular component to be parameterized already. After assembly, the user can inspect the outputs and use them directly or prepare portable minimization inputs as described below.

## What it can do

| Mode or command | Purpose | Configuration |
|---|---|---|
| `--mode receptor` | Insert a receptor or receptor-ligand block into a prepared membrane and solvent environment | `receptor_replace.yaml` |
| `--mode lig` | Replace a bound ligand while retaining the receptor | `ligand_replace.yaml` |
| `--mode chl` | Restore experimental cholesterol and correct leaflet composition | `cholesterol_restore.yaml` |
| `--prepare-minimization` | Prepare minimization inputs immediately after an assembly run | `minimization` section in the assembly YAML |
| `prepare-minimization` | Prepare minimization inputs from an existing matching GRO/topology bundle | Command-line inputs or `minimization.yaml` |

In receptor mode, ligands selected together with the incoming receptor undergo the same rigid transformation and are inserted as one block.

## Input requirements

The environment coordinates may come from:

- CHARMM-GUI GROMACS output, such as `step5_input.gro`.
- Another molecular-system builder.
- A frame exported from an MD trajectory.
- Another prepared PDB/GRO structure with matching topology files.

For an MD frame, use the topology belonging to that trajectory and preserve atom order, residue names, whole molecules, and box vectors. GRO is preferred because it preserves the periodic box and four-character GROMACS residue names.

The environment must have:

- A readable `topol.top` with a final `[ molecules ]` section.
- All referenced force-field and molecular topology files.
- An unambiguous `[ moleculetype ]` and `[ atoms ]` definition for every retained component.
- Coordinate atom names and ordering that match the topology.

The incoming receptor or ligand must also have matching GROMACS parameters.

Validated periodic processing requires an axis-aligned orthorhombic box. Supply a box-bearing GRO/PDB file or specify `[a, b, c, alpha, beta, gamma]` in Å and degrees in the YAML.

## Installation

Python 3.10 or newer is required. Installation includes NumPy, SciPy, PyYAML, and MDAnalysis.

From a downloaded or cloned source directory:

```bash
cd /path/to/GMXtransplant
python3 -m pip install .
gmxtransplant --help
```

For an editable installation:

```bash
python3 -m pip install -e .
```

You can also run directly from the source checkout:

```bash
python3 run_pipeline.py --help
```

### Optional dependencies

| Dependency | Required for |
|---|---|
| Open Babel (`obabel`) | Cholesterol conversion when `cholesterol.convert_to_charmm36: true` |
| PyMOL with `pymol2` | `pymol_align` and `pymol_cealign` alignment |
| GROMACS | Running the generated GROMACS minimization workflow |
| OpenMM | Running the generated OpenMM minimization workflow |

The default `mask_fit` alignment does not require PyMOL. Neither GROMACS nor OpenMM is required to prepare minimization folders.

To install OpenMM support in the current environment:

```bash
python3 -m pip install '.[minimize]'
python3 -m openmm.testInstallation
```

Here, `minimize` is an optional dependency group used by `pip`; this command only installs dependencies.

A Conda environment can also be used:

```bash
conda create -n gmxtransplant -c conda-forge \
  python=3.11 pip numpy scipy pyyaml mdanalysis openbabel

conda activate gmxtransplant
python -m pip install /path/to/GMXtransplant
```

Install OpenMM in that environment if you intend to run the OpenMM workflow:

```bash
conda install -c conda-forge openmm
```

## Quick start

Create a separate working directory and generate an example configuration:

```bash
mkdir my_gmxtransplant_job
cd my_gmxtransplant_job

gmxtransplant --show-example receptor > receptor_replace.yaml
```

Edit the input paths, selections, topology sources, and settings for your system.

Validate the configuration:

```bash
gmxtransplant --mode receptor -i receptor_replace.yaml --dry-run
```

Run assembly:

```bash
gmxtransplant --mode receptor -i receptor_replace.yaml
```

To prepare minimization inputs immediately after assembly:

```bash
gmxtransplant --mode receptor -i receptor_replace.yaml --prepare-minimization
```

The same option works with `--mode lig` and `--mode chl`.

To assemble without preparing minimization inputs, omit `--prepare-minimization` and leave `minimization.enabled: false` in the YAML.

## Configuration files

Generate the example corresponding to your task:

```bash
gmxtransplant --show-example receptor > receptor_replace.yaml
gmxtransplant --show-example lig > ligand_replace.yaml
gmxtransplant --show-example chl > cholesterol_restore.yaml
gmxtransplant --show-example minimize > minimization.yaml
```

`--show-example` only prints a configuration file. It does not assemble or minimize a system.

Examples include explanatory comments by default. To omit comments:

```bash
gmxtransplant --show-example receptor --no-comments > receptor_replace_simple.yaml
```

To include comments explicitly:

```bash
gmxtransplant --show-example receptor --comments > receptor_replace.yaml
```

Both forms contain the same active settings.

### Paths

Assembly examples group file paths under `paths:`. Define a path once and reference it elsewhere:

```yaml
paths:
  environment_gro: environment.gro

target_box:
  path: "${environment_gro}"
```

References must occupy the entire value. References inside the `paths:` block and environment-variable expansion are not supported.

Relative paths are resolved from the directory where the command is run. Absolute paths can be used from any directory.

Unknown and duplicate YAML keys are rejected.

### Receptor or receptor-ligand insertion

The main settings are:

- `target_box`: the prepared system supplying the retained environment.
- `replacement_structure`: the incoming receptor and any selected bound ligands.
- `alignment`: the fitting method and corresponding atom selections.
- `replacement_ligands`: charge and ITP metadata for incoming ligands.
- `topology`: topology sources for the receptor and environment.

The selection in `target_box.receptor_mask` is removed. Complete residues selected by `replacement_structure.receptor_mask` are transformed and inserted.

To insert a receptor-ligand complex, include its bound ligands in the replacement mask and list their metadata under `replacement_ligands`.

The receptor example contains protein positions 1 through 961, LIG1 at 962, and LIG2 at 963. Replace these example selections with your system's sequential residue positions.

For one ligand, remove the LIG2 metadata and path, and adjust the incoming mask accordingly. Set the target removal mask independently. Alignment masks should select corresponding protein atoms.

`mask_fit` requires:

- Equal atom counts.
- Matching atom-name sequences.
- Matching per-residue grouping.
- At least three non-collinear points.

`pymol_align` and `pymol_cealign` can determine correspondence through PyMOL.

For truncated PDB residue names such as `POP`, `TIP`, or `CHL`, enable `name_restoration`. A unique candidate is renamed directly; multiple candidates are distinguished using ITP atom counts. Ambiguous assignments stop the run.

### Ligand replacement

`ligand_replace.structure_path` specifies the prepared system containing the original ligand.

`original_ligand.resname` selects the residue to replace. If the name occurs more than once, specify a one-based occurrence, such as `UJU:2`.

Available placement methods:

| Method | Behavior |
|---|---|
| `autofit` | Match uniquely named heavy atoms and perform a Kabsch fit |
| `pairfit` | Fit explicitly paired atom-name lists |
| `nofit` | Preserve the incoming coordinates because they are already in the target coordinate frame |

The receptor is retained. Protein-ligand contacts are reported, while protein residues are protected from clash-based removal.

### Cholesterol restoration

The experimental receptor is aligned to the receptor in the target environment. Only the experimental cholesterol coordinates are transferred.

With `convert_to_charmm36: true`, existing cholesterol hydrogen atoms are removed, hydrogens are regenerated with Open Babel, and the structure is mapped and checked against the embedded 74-atom CHARMM36 cholesterol reference.

Missing heavy-atom reconstruction is controlled by:

- `max_missing_heavy_atoms`.
- `max_heavy_atom_fit_rmsd`.

Set `write_diagnostics: true` to retain conversion and placement intermediates.

`composition.lipid_targets` specifies desired upper and lower leaflet counts. Surplus molecules can be removed outside `distance_from_protein`; missing lipids are not created.

Salt is compared with `composition.reference_system_path`, which supplies the reference concentration.

## Minimization

GMXtransplant prepares minimization inputs for GROMACS and OpenMM. You then run one of the generated scripts locally or submit it through SLURM.

There are two ways to prepare the same folder structure:

| Option | When to use it | Inputs used |
|---|---|---|
| `--prepare-minimization` | During an assembly run | The GRO and topology produced by that run |
| `prepare-minimization` | When matching coordinates and topology already exist | Files specified on the command line or in `minimization.yaml` |

The leading `--` makes `--prepare-minimization` an option added to an assembly command. Without the leading `--`, `prepare-minimization` is a separate command for existing files.

**Both prepare inputs. Neither executes minimization.**

### Prepare inputs after assembly

Add the option to any assembly mode:

```bash
gmxtransplant --mode receptor -i receptor_replace.yaml --prepare-minimization
```

Configure the `minimization` section in that mode's YAML. Set `topology.enabled: true` so assembly produces the matching topology.

Alternatively, enable preparation in the YAML:

```yaml
minimization:
  enabled: true
  output_dir: minimization_inputs
  engines: [gromacs, openmm]
```

Then run assembly without the extra flag:

```bash
gmxtransplant --mode receptor -i receptor_replace.yaml
```

Preparation uses the freshly generated GRO and `topol.top`. The standalone `coordinates_path`, `topology_path`, and `include_dir` fields do not select different inputs for this workflow.

### Prepare inputs for an existing system

Use the separate command when assembly is complete or the prepared system comes from another source:

```bash
gmxtransplant prepare-minimization \
  --coordinates step5_input.gro \
  --topology topol.top \
  --output minimization_inputs
```

This does not repeat receptor insertion, ligand replacement, or cholesterol restoration.

To configure physical settings and computing resources, generate and edit the standalone example:

```bash
gmxtransplant --show-example minimize > minimization.yaml
```

Here, `minimize` selects the example configuration to print. After editing its paths and settings:

```bash
gmxtransplant prepare-minimization -i minimization.yaml --dry-run
gmxtransplant prepare-minimization -i minimization.yaml
```

`--dry-run` checks the configuration without reading coordinate or topology files. The normal command copies and validates those files.

Available standalone options:

| Option | Purpose |
|---|---|
| `--coordinates FILE.gro` | Select the input coordinates |
| `--topology FILE.top` | Select the matching topology |
| `--include-dir DIRECTORY` | Add a search directory for topology includes |
| `--output DIRECTORY` | Choose a new output directory |
| `--engines gromacs` | Prepare only the GROMACS folder |
| `--engines openmm` | Prepare only the OpenMM folder |

Command-line values override the corresponding YAML settings. Both engines are prepared by default.

For preparation after assembly, choose the engines through `minimization.engines` in the assembly YAML.

### Generated folders

The default output directory is `minimization_inputs/`.

| Location | Contents |
|---|---|
| Root directory | README and a manifest containing settings, checksums, and validation results |
| `gromacs/` | `input.gro`, `topol.top`, `toppar/`, audited topology, `minimization.mdp`, `run.sh`, `submit_cpu.slurm`, a result-summary script, and README |
| `openmm/` | `input.gro`, `topol.top`, `toppar/`, audited topology, `minimization.yaml`, `minimize_openmm.py`, supporting Python files, `requirements.txt`, `run.sh`, CPU/GPU submission scripts, and README |

Each engine folder contains its own copies of the inputs and can be moved independently. GMXtransplant is not required on the execution machine.

Preparation checks atom names, order, counts, charge, finite coordinates, and the periodic box. It copies topology dependencies and adjusts their include paths.

An existing output directory is not overwritten. Choose a new destination when preparing another set of inputs.

### Run GROMACS minimization

Activate your GROMACS environment and enter the generated folder:

```bash
cd minimization_inputs/gromacs
```

Run locally:

```bash
bash run.sh
```

Or submit a CPU job from that directory:

```bash
sbatch submit_cpu.slurm
```

The local runner uses `gmx` with one thread-MPI rank and multiple OpenMP threads. The SLURM script uses `gmx_mpi` through `srun`.

Set `GMX` to the executable path if needed. The runner creates the TPR with `grompp`, then runs minimization.

### Run OpenMM minimization

Enter the generated folder and install its dependencies:

```bash
cd minimization_inputs/openmm
python3 -m pip install -r requirements.txt
```

Run locally on the CPU:

```bash
bash run.sh --platform CPU
```

Or run locally on a compatible GPU:

```bash
bash run.sh --platform CUDA
```

For SLURM execution, submit one of the supplied scripts from that directory:

```bash
sbatch submit_cpu.slurm
```

Alternatively:

```bash
sbatch submit_gpu.slurm
```

The GPU template requests one GPU and selects CUDA. The runner reports an error if the requested platform cannot initialize.

### CPU and submission settings

Configure resources in the preparation YAML:

```yaml
minimization:
  resources:
    cpus_per_task: 8
    gromacs_mpi_tasks: 1
    time: "00:30:00"
    memory: 8G
```

The generic SLURM templates request one node. GROMACS can use multiple MPI tasks, each with the configured CPU threads. OpenMM uses one process with multiple CPU threads or one GPU.

Review the generated submission script, adjust resource requests, and add any software-environment setup needed before submitting. Submit from the engine folder.

Local CPU threading uses:

- `OMP_NUM_THREADS` for GROMACS.
- `OPENMM_CPU_THREADS` for OpenMM.

In the supplied SLURM scripts, thread counts follow `SLURM_CPUS_PER_TASK`. The scheduler controls GPU visibility.

### Scientific settings

Review minimization settings against the force field and protocol used to prepare your system.

| Engine | Settings in the preparation YAML | Generated settings file |
|---|---|---|
| GROMACS | `minimization.gromacs` | `gromacs/minimization.mdp` |
| OpenMM | Physical fields directly under `minimization`, such as `tolerance_kj_mol_nm` and `nonbonded_cutoff_nm` | `openmm/minimization.yaml` |

GROMACS defaults to steepest descent, PME, 1.2 nm cutoffs, and LJ force switching from 1.0 nm. OpenMM uses L-BFGS. These are configurable starting settings.

The convergence criteria differ:

- GROMACS `emtol` measures maximum force.
- OpenMM's tolerance uses the RMS minimizer objective gradient.

OpenMM potential switching also differs from GROMACS force switching. Matching numerical values does not guarantee equivalent settings or results.

The OpenMM example selects CPU. A platform passed to `run.sh` overrides the YAML setting. If `platform: auto` is specified during preparation, the exported YAML selects CPU; the GPU submission script explicitly selects CUDA.

### Optional restrained GROMACS stage

To minimize with position restraints first, followed by an unrestrained stage:

```yaml
minimization:
  gromacs:
    restrained_defines: {POSRES: 1}
    defines: {}
```

This requires compatible position-restraint definitions already present in the topology. Preparation checks that restraints are active in the first stage and inactive in the final stage.

The runner executes `restrained.mdp` first, using `input.gro` as the restraint reference. After confirmed convergence, it uses the resulting coordinates for `minimization.mdp`.

With empty `restrained_defines`, only the single minimization stage is generated.

OpenMM preparation rejects active GROMACS `position_restraints`. Use the GROMACS workflow for that protocol.

### Minimization results

Each runner creates a new `results/` directory inside its engine folder.

| Engine | Main results |
|---|---|
| GROMACS | `results/minimization.gro`, log and energy files, TPR, resolved MDP/topology, preprocessing logs, and a JSON convergence summary |
| OpenMM | `results/minimized.gro`, `results/minimization_report.txt`, and `results/minimization_report.json` |

Existing results are not overwritten. Preserve them by renaming the `results/` directory before another run.

Inspect the convergence report and final structure before proceeding. Status 7 indicates failure or unconfirmed convergence; output coordinates may still exist.

Minimization does not replace equilibration.

### Topology requirements for minimization preparation

All literal includes must resolve, including files in inactive conditional branches.

Supported preprocessing includes:

- `#ifdef`, `#ifndef`, `#else`, and `#endif`.
- Object macros.
- Literal includes.

Expressions such as `#if`, function-like macros, and cyclic include graphs require a preprocessed topology.

OpenMM runs the expanded `audit_openmm.top` created during preparation. Regenerate the inputs if changing its topology defines. The expansion includes OpenMM's `FLEXIBLE` definition for water bonds and angles; the runner then applies the configured rigid-water constraints.

Preparation validates input identity and consistency. Engine compatibility is checked at execution, including GROMACS preprocessing with `-maxwarn 0`.

If preparation fails after assembly, the validated assembly outputs remain available.

Engine documentation:

- [GROMACS MDP options](https://manual.gromacs.org/current/user-guide/mdp-options.html)
- [OpenMM platform properties](https://docs.openmm.org/latest/userguide/library/04_platform_specifics.html)

## Residue masks

Masks use one-based sequential residue positions in coordinate-file order. They do not use the residue numbers printed in PDB/GRO files.

| Mask | Meaning |
|---|---|
| `:1-298` | Residues 1 through 298 |
| `:1,5,10-20` | Residues 1, 5, and 10 through 20 |
| `:1-298@CA` | CA atoms in residues 1 through 298 |
| `:1-298@C,CA,N,O` | Backbone atoms in residues 1 through 298 |
| `:1-298&!:297-298` | Residues 1 through 296 |
| `@CA` | Every CA atom |

Unsupported cpptraj syntax is rejected.

## Clash, box, and charge validation

Environment preflight runs before alignment or coordinate output.

| Setting | Behavior |
|---|---|
| `strict` | Validate the existing coordinates and cell |
| `repair` | Attempt a validated whole-system rotation and translation |
| `off` | Disable preflight validation; required for nonperiodic clash detection |

The validator checks finite coordinates, the cell, direct and periodic hard contacts, and output-format rounding. Repair does not rotate molecules independently.

Clash detection compares environment residues with the inserted block using class-specific cutoffs and minimum-image distances. If an atom pair crosses the cutoff, the complete environment residue is removed.

Existing environment-environment contacts are left unchanged. Classes listed in `keep_classes` are retained, with their contacts included in the reports.

Charge accounting uses ITP `[ atoms ]` charges where available. Explicit overrides are for species without a usable ITP and cannot contradict an available ITP charge.

Neutralization removes eligible counterions only when this moves the system toward `target_net_charge`. Ion candidates can be excluded near receptor, ligand, and lipid heavy atoms.

## Assembly outputs

With the corresponding output options enabled, a successful run produces:

- `step5_input.gro` and `step5_input.pdb`.
- An alignment or placement inspection PDB.
- Text and JSON reports.
- `topol.top` and a self-contained `toppar/` directory.
- `index.ndx`.

Requested minimization preparation additionally creates `minimization_inputs/`. Minimized coordinates appear only after a generated runner is executed.

The final atom order is protein, ligand/cofactor, lipid, ion, then water. Written coordinate files are checked for atom names, order, coordinates, box dimensions, and remaining clashes. When topology output is enabled, topology atom order, counts, and charge are also checked.

Before using the outputs, review:

- Alignment RMSD and the inspection PDB.
- Removed molecules and `flagged_but_kept` contacts.
- Lipid-removal fractions and final leaflet counts.
- Unknown residue names or topology ambiguity.
- Final net charge.
- Coordinate and topology validation results.

## Common problems

**No usable box**

Use a GRO containing the correct periodic box or supply the six box dimensions from an authoritative source.

**Unknown or ambiguous molecule type**

Supply the correct topology files, add the ligand ITP, restore truncated PDB residue names, or provide a specific `moleculetype_overrides` entry.

**Alignment selections differ**

Compare the atom counts and atom-name order in the two masks. Use `pairfit` for an explicitly selected shared ligand core or PyMOL alignment for non-identical receptor sequences.

**Open Babel not found**

Install it through Conda:

```bash
conda install -c conda-forge openbabel
```

Alternatively, set `cholesterol.obabel_command` to the executable path.

**Preparation refuses an existing directory**

Choose a new `--output` directory for standalone preparation, or change `minimization.output_dir` in the assembly YAML.

**OpenMM is unavailable**

Install the dependencies from the generated OpenMM folder:

```bash
python3 -m pip install -r requirements.txt
```

**Minimization exits with status 7**

Inspect the result reports and logs. Check parameters, overlaps, and nonbonded settings before increasing iteration limits.

**Large lipid-removal warning**

Inspect the inserted pose, periodic box, and selections before proceeding.

## Tests and build

Run the tests from the repository root:

```bash
python3 -m unittest discover -v
```

Build the source archive and wheel:

```bash
python3 -m pip install build
python3 -m build
```

Build outputs are written under `dist/`.

GitHub Actions tests Python 3.10, 3.11, and 3.12, builds the distributions, and runs GROMACS and OpenMM smoke tests in environments with the corresponding engine installed.

## Limitations

- Validated periodic processing supports PDB/GRO coordinates with orthorhombic cells.
- PDB has finite-width atom and residue fields; use GRO as the authoritative GROMACS coordinate output.
- Assembly topology parsing requires preprocessed definitions when conditionals occur inside `[ moleculetype ]`, `[ atoms ]`, or `[ molecules ]`.
- Cholesterol composition correction removes surplus molecules but does not add missing lipids or solvent.
- Prepared minimization supports GROMACS steepest descent and OpenMM L-BFGS.
- The two engines' force and convergence settings are not assumed to be equivalent.
