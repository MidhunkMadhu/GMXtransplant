# GMXtransplant

GMXtransplant assembles membrane-protein systems using an existing prepared membrane, solvent, and ion environment. You can replace a protein or protein-ligand complex, replace a bound ligand, or restore experimentally resolved cholesterol.

This is a structure-assembly tool and expects every molecular component to be parameterized already. After assembly, inspect the outputs and use them directly, or prepare portable OpenMM inputs for restrained clash relaxation.

## Choose a task

| Task | Command option | Example configuration |
|---|---|---|
| Insert a protein or protein-ligand complex into a prepared environment | `--mode protein` | `protein_replace.yaml` |
| Replace a bound ligand while retaining the protein | `--mode lig` | `ligand_replace.yaml` |
| Restore experimental cholesterol and adjust leaflet composition | `--mode chl` | `cholesterol_restore.yaml` |

Protein mode supports membrane proteins generally, including receptors, channels, and transporters. Any ligands included in the incoming selection move with the protein during alignment.

## Install

Use Python 3.10 or newer. From the downloaded or cloned source directory:

```bash
python3 -m pip install .
gmxtransplant --help
```

For an editable installation, use `python3 -m pip install -e .`. To run from source without installing the package, use `python3 run_pipeline.py` in place of `gmxtransplant`.

NumPy, SciPy, PyYAML, and MDAnalysis are installed with the package. Optional dependencies are:

- **OpenMM:** needed only when running the generated minimization script.
- **PyMOL with pymol2:** needed for `pymol_align` or `pymol_cealign`. The default `mask_fit` does not need PyMOL.
- **Open Babel:** needed when converting experimental cholesterol with `convert_to_charmm36: true`.

## Quick start

Work in a separate directory and generate the example for your task:

```bash
mkdir my_gmxtransplant_job
cd my_gmxtransplant_job
gmxtransplant --show-example protein > protein_replace.yaml
```

Edit the paths, residue selections, topology sources, and scientific settings. Then validate the configuration and run assembly:

```bash
gmxtransplant --mode protein -i protein_replace.yaml --dry-run
gmxtransplant --mode protein -i protein_replace.yaml
```

`--dry-run` checks the configuration without opening molecular input files. Input identity, topology, and geometry checks happen during the normal run.

Use `--show-example lig` or `--show-example chl` for the other modes. Add `--no-comments` to print a compact example with the same active settings.

## Prepare your inputs

The target environment can come from CHARMM-GUI, another system builder, or an MD trajectory frame. Supply matching coordinate and GROMACS topology files, including all force-field and molecule definitions. Incoming proteins and ligands also need matching parameters from a compatible force field.

GRO is preferred: it preserves the periodic box and four-character residue names such as POPC and TIP3. If using an MD frame, preserve atom order, whole molecules, and the correct box. Validated periodic processing requires an axis-aligned orthorhombic cell.

Example files group paths under `paths:` and reuse them as `"${name}"`. Relative paths are resolved from the directory where you run the command. Unknown or duplicate YAML keys are rejected.

### Protein insertion

Set these together in `protein_replace.yaml`:

- `target_box.protein_mask`: the complete old protein/ligand block to remove.
- `replacement_structure.protein_mask`: the complete incoming block to insert, including any bound ligands.
- `alignment`: corresponding protein atoms used for fitting.
- `replacement_ligands`: residue names, ITP files, and charges for incoming ligands. Use `[]` when there are none.
- `topology.protein_toppar_dir` and `topology.environment_toppar_dir`: the protein and retained environment parameter sources.

The alignment transform is applied to the whole incoming block. Fit a stable protein region so that a changed ligand pose does not bias placement. For example, a receptor complex can be aligned using corresponding receptor backbone atoms.

`mask_fit` requires matching atom-name order and residue grouping, equal atom counts, and at least three non-collinear points. Use PyMOL alignment when correspondence must be established between different protein sequences.

### Ligand replacement

In `ligand_replace.yaml`, select the prepared system and the original ligand residue. Use `RESNAME:N` to select a specific occurrence when its name appears more than once.

| Placement | Behavior |
|---|---|
| `autofit` | Fit matching, uniquely named heavy atoms |
| `pairfit` | Fit explicit pairs of atom names |
| `nofit` | Keep incoming coordinates that are already in the target frame |

The protein is retained and protected from clash-based removal. Protein-ligand contacts are reported for inspection.

### Cholesterol restoration

In `cholesterol_restore.yaml`, the experimental protein is aligned to the target protein and the associated cholesterol coordinates are transferred.

Optional CHARMM36 conversion regenerates hydrogens with Open Babel and validates the result against the included cholesterol reference. The example explains the limits for reconstructing missing heavy atoms.

Leaflet composition settings can remove surplus molecules outside the protected region around the protein. They do not create missing lipids or solvent.

**This restoration option supports cholesterol only. There is currently no corresponding option to preserve and transfer other experimentally resolved bound lipids.** Existing lipids in the target membrane are retained according to the clash and composition settings.

### Residue selections

Masks use **one-based sequential residue positions in coordinate-file order**, not the residue numbers printed in PDB or GRO files.

| Mask | Selection |
|---|---|
| `:1-298` | All atoms in residues 1 through 298 |
| `:1,5,10-20` | Selected residues or ranges |
| `:1-298@CA` | CA atoms in those residues |
| `:1-298@C,CA,N,O` | Backbone atoms in those residues |
| `:1-298&!:297-298` | Exclude residues 297 and 298 |

## Optional OpenMM minimization

Prepare one `openmm_minimization/` folder for a short clash-relaxation step. **Preparation does not run or submit minimization.** Execute the generated runner separately when ready.

Harmonic positional restraints keep protein and ligand heavy atoms near their starting coordinates. Hydrogens, recognized membrane lipids, water, and ions remain unrestrained. These restraints allow small movements; they do not freeze the protein or ligand.

### Prepare the folder

| Use case | Command |
|---|---|
| Assemble and then prepare inputs from that run | `gmxtransplant --mode protein -i protein_replace.yaml --prepare-minimization` |
| Prepare inputs from an existing matching system | `gmxtransplant prepare-minimization --coordinates step5_input.gro --topology topol.top` |

`--prepare-minimization` is an assembly option and works with all three modes. It requires `topology.enabled: true` and uses the freshly generated GRO and topology. Setting `minimization.enabled: true` in the assembly YAML has the same effect.

`prepare-minimization` is a separate command for existing files and does not repeat assembly. Use `--include-dir` for an additional topology search directory and `--output` to choose a different destination.

For adjustable standalone settings:

```bash
gmxtransplant --show-example minimize > minimization.yaml
# Edit the input paths and settings.
gmxtransplant prepare-minimization -i minimization.yaml
```

The folder contains coordinate and topology copies, all topology dependencies, an expanded topology for execution, a settings YAML, the Python runner and its support files, `run.sh`, generic multicore CPU and GPU SLURM templates, a README, and an input audit with checksums. Move the complete folder to your execution machine; GMXtransplant itself is not required there.

### Settings

The assembly and standalone YAML examples use the same `minimization` settings:

```yaml
minimization:
  enabled: false
  output_dir: openmm_minimization
  restraint_force_constant_kj_mol_nm2: 1000.0
  restraint_residue_classes: {}
  tolerance_kj_mol_nm: 10.0
  max_iterations: 5000
  resources:
    cpus_per_task: 8
    time: "00:30:00"
    memory: 8G
```

Restraints use the input coordinates as references and a default force constant of 1000 kJ/mol/nm². Protein molecules are recognized from their amino-acid residues, including attached caps and modifications. Unrecognized molecules are treated as ligands. Classify custom lipids and unusual residues explicitly, for example `restraint_residue_classes: {MYLIP: lipid}`. Supported classes are `protein`, `ligand`, `lipid`, `water`, and `ion`.

The runner uses OpenMM L-BFGS with a finite iteration limit and an RMS objective-gradient tolerance. Review the example's nonbonded settings for your force field. OpenMM potential switching differs from LJ force switching; this is a configurable clash-relaxation step, not a universal simulation protocol.

Topology `position_restraints` must be inactive because the runner adds its own protein/ligand restraints. Preparation rejects active definitions instead of ignoring them.

### Run separately

From the generated folder:

```bash
cd openmm_minimization
python3 -m pip install -r requirements.txt
bash run.sh
```

The default run uses multiple CPU threads. Set `OPENMM_CPU_THREADS` to adjust the local thread count, or select a GPU with `bash run.sh --platform CUDA`. HIP and OpenCL can also be selected when available.

For SLURM, review the resource directives and software setup, then submit **one** template from this folder:

```bash
sbatch submit_cpu.slurm
# Alternatively:
sbatch submit_gpu.slurm
```

The CPU template uses one process with multiple threads. The GPU template requests one GPU and selects CUDA. An unavailable requested platform produces an error.

### Inspect results

The runner writes:

- `results/minimized.gro`
- `results/minimization_report.txt`
- `results/minimization_report.json`

Reports include convergence, energies, forces, selected restraint counts, inferred ligand names, and heavy-atom displacement. Energies include the positional restraints. Exit code 7 means failure or unconfirmed convergence; inspect the report before using any output coordinates.

Existing preparation and results directories are never overwritten. Choose a new destination for preparation, or preserve an earlier `results/` directory by renaming it before rerunning. Minimization does not replace equilibration.

## Assembly checks and outputs

Assembly checks the periodic environment, aligns the incoming component, handles clashes, optionally corrects charge, and rebuilds the topology molecule list.

Clash cutoffs act on complete environment residues. Classes in `clash_detection.keep_classes` are retained and their contacts are reported. Neutralization removes eligible existing counterions; it does not add ions or set a new salt concentration.

With the corresponding outputs enabled, a successful run produces GRO/PDB coordinates, an inspection PDB, text/JSON reports, `topol.top`, `toppar/`, and `index.ndx`. Review the alignment, removed molecules, retained clashes, lipid counts, net charge, and coordinate/topology checks before continuing.

The commented example files describe box validation, name restoration for truncated PDB names, clash thresholds, ion protection, and output options.

## Limitations

- All components must already have compatible parameters. GMXtransplant does not parameterize new molecules.
- Validated periodic processing supports orthorhombic boxes. Use GRO as the authoritative simulation coordinate output.
- Experimental bound-lipid restoration supports cholesterol only; other bound lipids have no dedicated preservation/restoration option.
- Composition correction removes surplus lipids but does not build missing lipids or solvent.
- Conditional topology definitions may need preprocessing. Minimization preparation supports literal includes, object macros, and `#ifdef`/`#ifndef` branches; unsupported expressions and cyclic includes are rejected.
- Preparation checks input consistency without importing OpenMM. OpenMM parsing and force-field compatibility are checked when the runner executes.
- Exact overlaps, poor placement, or incorrect parameters may need correction before minimization can succeed.

## Development

```bash
python3 -m unittest discover -v
python3 -m pip install build
python3 -m build
```

GitHub Actions tests Python 3.10, 3.11, and 3.12, builds the package, and runs OpenMM tests with the optional dependency installed.
