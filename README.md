# GMXtransplant

**New first mode: two-folder CHARMM-GUI protein transplantation.** Supply only reference and transplant folders; proteins and ligands are detected and replaced automatically. The other assembly modes remain available. See [the mode guide](CHARMPROT.md), [supplied example](examples/charmprot/README.md), and [PyMOL/VMD comparisons](VISUALIZATION.md).


GMXtransplant assembles membrane-protein systems using an existing prepared membrane, solvent, and ion environment. You can replace a protein or protein-ligand complex, replace a bound ligand, or restore experimentally resolved cholesterol.

**Two ways to work, one pipeline.** The desktop application is the main way to use GMXtransplant: browse for inputs, adjust settings in a guided form, run the bundled examples in one click, and follow progress and results in one window. For scripts, clusters and remote servers, every mode is also available as a complete command-line tool (`gmxtransplant`). Both run the same pipeline and use the same YAML configurations: every application run saves its setup as `run.yaml` in the output folder, which runs unchanged with `gmxtransplant --mode <mode> -i run.yaml`.

This is a structure-assembly tool and expects every molecular component to be parameterized already. After assembly, inspect the outputs and use them directly, or prepare portable OpenMM inputs for restrained clash relaxation.

## Choose a task

| Task | Command option | Example configuration |
|---|---|---|
| Replace proteins and ligands using only two CHARMM-GUI folders | `--mode charmprot` | `charmprot.yaml` |
| Insert a protein or protein-ligand complex into a prepared environment | `--mode protein` | `protein_replace.yaml` |
| Replace a bound ligand while retaining the protein | `--mode lig` | `ligand_replace.yaml` |
| Restore experimental cholesterol and adjust leaflet composition | `--mode chl` | `cholesterol_restore.yaml` |
| Add a ligand or protein binder in the water above/below a membrane protein ([guide](ADDBINDER.md), first version) | `--mode addbinder` | `addbinder.yaml` |

Protein mode supports membrane proteins generally, including receptors, channels, and transporters. Any ligands included in the incoming selection move with the protein during alignment.

## Install

Use Python 3.10 or newer on Linux or macOS. From the downloaded or cloned source
directory, one command installs the desktop application, the command-line tool
and everything they need:

```bash
python3 -m pip install .
gmxtransplant-gui        # desktop application
gmxtransplant --help     # command-line tool
```

A separate environment keeps GMXtransplant's dependencies apart from other
software: `python3 -m venv .venv && source .venv/bin/activate` before installing
(or `conda create -n gmxtransplant python=3.12 && conda activate gmxtransplant`).

### What is installed, and why

| Package | Why GMXtransplant needs it |
|---|---|
| NumPy, SciPy | Coordinate arithmetic, rigid-body fitting and distance searches |
| PyYAML | Reading and writing the YAML configurations |
| MDAnalysis | Reading and writing GRO/PDB/MOL2 structures and periodic-box geometry |
| PySide6 (Qt) | The desktop application's windows. The command-line tool never loads it, so it also works on servers without a display |
| openbabel-wheel (Open Babel) | The `obabel` command that adds hydrogens and CHARMM36 names to experimental cholesterol (cholesterol restoration with `convert_to_charmm36: true`). Installed automatically on Linux and on macOS 14 or newer |

### Optional, installed separately

| Package | Needed for | How to install |
|---|---|---|
| OpenMM | Running a prepared minimization folder | `python3 -m pip install ".[minimize]"` or `conda install -c conda-forge openmm` |
| PyMOL (with pymol2) | The `pymol_align` / `pymol_cealign` alignment methods; the default `mask_fit` does not need it | `conda install -c conda-forge pymol-open-source` |
| VMD | The **Open in VMD** comparison view | From the VMD website |

On macOS, PyMOL and VMD installed as applications in `/Applications` are found
automatically. For a copy elsewhere, set `GMXTRANSPLANT_PYMOL` or `GMXTRANSPLANT_VMD`
to its full path.

### If a dependency cannot be installed

- **Open Babel on macOS 12 or 13:** no ready-built pip package exists, so it is
  skipped. Install it with `brew install open-babel` or
  `conda install -c conda-forge openbabel`. On other systems where the pip package
  fails, install GMXtransplant with `GMXTRANSPLANT_NO_OPENBABEL=1 python3 -m pip install .`
  and use one of those commands (or `sudo apt-get install openbabel`). Only
  cholesterol restoration needs it; check with `obabel -V`.
- **Qt on older Linux systems or some clusters:** install just the command-line
  tool with `GMXTRANSPLANT_CLI_ONLY=1 python3 -m pip install .`. Every mode, example
  and option remains available through `gmxtransplant`; the application can be
  added later with `python3 -m pip install "PySide6>=6.6,<7"`.

For an editable installation, use `python3 -m pip install -e .`. To run the command line from source without installing the package, use `python3 run_pipeline.py` in place of `gmxtransplant`.

## Desktop application

Launch the application from any directory:

```bash
gmxtransplant-gui
# Equivalent module entry point:
python -m gmxtransplant.gui
```

The desktop application supports macOS, Linux, and WSL 2 with WSLg. It includes an input
editor for all configuration sections, file and folder browsing, configuration import/export,
one-click examples, validation, readable progress summaries, cancellation, and result browsing.
Detailed logs are collapsed until **View Command Progress** is clicked. YAML is
generated automatically in the background; there is no YAML editing tab.
Select an output folder first. Runs write directly to `charmprot/`, `protein/`, `ligand/`, or
`cholesterol/`; example runs use `examples/<mode>/`. Repeating a run replaces its
generated products, including `toppar/`, without timestamped run directories.

Complete example inputs are included in the package for offline use. This adds
approximately 58 MB of installed data (the full wheel is about 12 MB compressed).
Read the [GUI guide](GUI_USAGE.md) for path handling, output replacement, and WSL setup.

On a MacBook Air, run these commands in Bash (or zsh) from the repository root:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install .
gmxtransplant-gui
```

Use a native Python installation matching your Mac's architecture. See the
[macOS instructions](GUI_USAGE.md#macos-and-macbook-air) for prerequisites
(including Open Babel on macOS 12 and 13).

## Quick start

### In the desktop application

1. Run `gmxtransplant-gui` and choose an output folder at the top of the window.
2. Open the **Examples** tab and click **Run example** on any card to see a complete
   run, or choose a mode on the **Configuration** tab and browse for your own inputs.
3. Click **Check configuration**, then **Run pipeline**, and follow the run on
   **Run & results**.

### On the command line

Work in a separate directory and generate the example for your task:

```bash
mkdir my_gmxtransplant_job
cd my_gmxtransplant_job
gmxtransplant --show-example protein > protein_replace.yaml
```

Edit the paths, residue selections, topology sources, and scientific settings. Then validate the configuration and run assembly:

```bash
gmxtransplant --mode protein -i protein_replace.yaml --dry-run
gmxtransplant --mode protein -i protein_replace.yaml --output results/
```

`--dry-run` checks the configuration without opening molecular input files. Input identity, topology, and geometry checks happen during the normal run.

`--output DIR` chooses where the generated files go. Without it they are written to the directory you ran the command from. Either way, the run's closing line names the absolute folder it wrote to.

Use `--show-example lig` or `--show-example chl` for the other modes. Add `--no-comments` to print a compact example with the same active settings.

Source-checkout example configurations are available for
[protein insertion](examples/protein_insertion/README.md),
[ligand replacement](examples/ligand_replacement/README.md), and
[cholesterol restoration](examples/cholesterol_restoration/README.md).
The ligand and cholesterol examples contain input folders for your own
structures. Their READMEs list the exact filenames and settings to supply.

## Prepare your inputs

The target environment can come from CHARMM-GUI, another system builder, or an MD trajectory frame. Supply matching coordinate and GROMACS topology files, including all force-field and molecule definitions. Incoming proteins and ligands also need matching parameters from a compatible force field.

GRO is preferred: it preserves the periodic box and four-character residue names such as POPC and TIP3. If using an MD frame, preserve atom order, whole molecules, and the correct box. Validated periodic processing requires an axis-aligned orthorhombic cell.

Example files group paths under `paths:` and reuse them as `"${name}"`. Unknown or duplicate YAML keys are rejected.

Every input path must be written out in full. No input is resolved against the working directory, against the configuration file's folder, or against any shared input root, and a relative input path is rejected with a message naming the offending key. The one exception is the examples shipped with the package: their paths resolve against the `examples/` folder installed alongside the code, so they run anywhere without referring to your filesystem.

Generated files are the mirror image: an output name such as `step5_input.gro` is a plain filename that lands in the run's output folder. That folder is `--output DIR` when you pass it, and otherwise the directory you ran the command from. Every run ends by printing the absolute path it wrote to.

### Protein insertion

Set these together in `protein_replace.yaml`:

- `target_box.protein_mask`: the complete old protein/ligand block to remove.
- `replacement_structure.protein_mask`: the complete incoming block to insert, including any bound ligands.
- `alignment`: corresponding protein atoms used for fitting.
- `replacement_ligands`: coordinate residue names and ITP files for incoming ligands; charge and molecule type are read automatically from the ITP. Use `[]` when there are none. Each ligand may also give `forcefield_path` (see below).
- `topology.protein_toppar_dir` and `topology.environment_toppar_dir`: the protein and retained environment parameter sources.

The alignment transform is applied to the whole incoming block. Fit a stable protein region so that a changed ligand pose does not bias placement. For example, a receptor complex can be aligned using corresponding receptor backbone atoms.

`mask_fit` requires matching atom-name order and residue grouping, equal atom counts, and at least three non-collinear points. Use PyMOL alignment when correspondence must be established between different protein sequences.

### Ligand replacement

In `ligand_replace.yaml`, select the prepared system and the original ligand residue. Use `RESNAME:N` to select a specific occurrence when its name appears more than once.

The incoming ligand needs three files: its coordinates (`new_ligand.coord_path`,
e.g. a `.mol2`), its ITP (`new_ligand.itp_path`), and the force field holding its
bond, angle and dihedral parameters (`new_ligand.forcefield_path`). A CHARMM-GUI
ligand ITP only lists which atoms are bonded; the parameters are in the
`forcefield.itp` CHARMM-GUI generated together with it. When `forcefield_path` is
empty, a `forcefield.itp` next to the ITP is used, otherwise the environment's.
Before writing the topology the run checks that every bonded term of the new
ligand has parameters, and stops naming any that are missing.

| Placement | Behavior |
|---|---|
| `autofit` | Fit matching, uniquely named heavy atoms |
| `pairfit` | Fit explicit pairs of atom names |
| `mcsfit` | Detect the largest shared substructure automatically (same elements and bonds), fit on its closest symmetric match, and report the pairs |
| `nofit` | Keep incoming coordinates that are already in the target frame |

The protein is retained and protected from clash-based removal. Protein-ligand contacts are reported for inspection.

### Cholesterol restoration

In `cholesterol_restore.yaml`, the experimental protein is aligned to the target protein and the associated cholesterol coordinates are transferred.

Optional CHARMM36 conversion regenerates hydrogens with Open Babel and validates the result against the included cholesterol reference. The example explains the limits for reconstructing missing heavy atoms.

Leaflet composition settings can remove surplus molecules outside the protected region around the protein. They do not create missing lipids or solvent.

#### Repacking before MD (`repack/`)

**Issue:** inserting the experimental cholesterols removes the lipids that
overlap them, and the composition step then trims lipids far from the protein.
The restored cholesterols start with few lipid neighbours (0–1 lipid heavy
atoms within 5 Å instead of about 45). Released straight into MD they left
their experimental poses within 25–50 ps in a D1R test, while the membrane
shrank by 20%.

**What is written:** every run also writes `repack/`, a system set up for a
short GROMACS equilibration (about 1.4 ns) that holds every heavy atom of the
protein, the ligands and the restored cholesterols, and leaves all other
lipids, water and ions free to fill the gaps. It has its own `topol.top`,
`toppar/`, `step5_input.gro`, `index.ndx` and mdp files. A sample run script
(`sample_run.sh`) is also saved. Afterwards, continue from the last GRO with
the main `topol.top`.

| Stage | Length | Protein, ligands, restored cholesterols (kJ/mol/nm²) | Other lipids |
|---|---|---|---|
| step6.0 minimization | 5000 steps | 4000 | free |
| step6.1 NVT | 125 ps | 4000 | free |
| step6.2 NVT | 125 ps | 2000 | free |
| step6.3 NPT | 125 ps | 1000 | free |
| step6.4 NPT | 1 ns | 1000 | free |

The restored cholesterols are their own moleculetype there (`toppar/CHLR.itp`:
CHL1 with all heavy atoms restrained in x, y, z); protein and ligands use their CHARMM-GUI
`POSRES` blocks, and the run warns if a heavy atom has none. In a D1R test
with the protein and restored cholesterols held for 1.4 ns, the lipid heavy
atoms within 5 Å of each restored cholesterol rose from 0–1 to 35–56 (normal:
59–73), the gaps closed within about 1 ns, and the poses held. The membrane
area keeps settling for 10–15 ns afterwards, so release the remaining restraints gradually. Change the stages,
restraints or temperature (310 K) under `cholesterol.repack`, or set
`enabled: false`.

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

Harmonic positional restraints keep protein and ligand heavy atoms and lipid head groups (the polar P, O, N and S atoms of lipids and sterols) near their starting coordinates. Lipid tails, hydrogens, water, and ions relax freely. These restraints allow small movements; they do not freeze the protein, ligand or membrane.

### Prepare the folder

| Use case | Command |
|---|---|
| Assemble and then prepare inputs from that run | `gmxtransplant --mode protein -i protein_replace.yaml --prepare-minimization` |
| Prepare inputs from an existing matching system | `gmxtransplant prepare-minimization --coordinates step5_input.gro --topology topol.top` |

`--prepare-minimization` is an assembly option and works with all four assembly modes. It requires `topology.enabled: true` and uses the freshly generated GRO and topology. Setting `minimization.enabled: true` in the assembly YAML has the same effect.

`prepare-minimization` is a separate command for existing files and does not repeat assembly. Use `--include-dir` for an additional topology search directory and `--output` to choose a different destination.

For adjustable standalone settings:

```bash
gmxtransplant --show-example minimize > minimization.yaml
# Edit the input paths and settings.
gmxtransplant prepare-minimization -i minimization.yaml
```

The folder contains coordinate and topology copies, all topology dependencies, an expanded topology for execution, a settings YAML, the Python runner and its support files, `run.sh`, generic multicore CPU and GPU SLURM templates, a README, and an input audit with checksums. Move the complete folder to your execution machine; GMXtransplant itself is not required there.

### Settings

Every setting has a default suited to CHARMM36 membrane systems, so an assembly
configuration needs only `minimization: {enabled: true}`. The settings you are
most likely to change are:

```yaml
minimization:
  enabled: true
  restraint_force_constant_kj_mol_nm2: 1000.0
  max_iterations: 5000
  platform: auto            # or CPU, CUDA, HIP, OpenCL
  restraint_residue_classes: {}   # e.g. {MYLIP: lipid} for a custom lipid
```

Nonbonded, constraint and SLURM `resources` settings remain available in YAML
for special cases; the GUI shows only the settings above.

Restraints use the input coordinates as references and a default force constant of 1000 kJ/mol/nm². Protein molecules are recognized from their amino-acid residues, including attached caps and modifications. Unrecognized molecules are treated as ligands. Classify custom lipids and unusual residues explicitly, for example `restraint_residue_classes: {MYLIP: lipid}`. Supported classes are `protein`, `ligand`, `lipid`, `water`, and `ion`.

The runner uses OpenMM L-BFGS with a finite iteration limit and an RMS objective-gradient tolerance. Review the example's nonbonded settings for your force field. OpenMM potential switching differs from LJ force switching; this is a configurable clash-relaxation step, not a universal simulation protocol.

Topology `position_restraints` must be inactive because the runner adds its own restraints. Preparation rejects active definitions instead of ignoring them.

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

Clash cutoffs act on complete environment residues. Every mode uses the same `clash_detection.threshold` of 0.80 Å by default, and no molecule class carries an implicit cutoff of its own: cholesterol is judged exactly like any other lipid. Override a single class through `clash_detection.thresholds` when a system needs it. Classes in `clash_detection.keep_classes` are retained and their contacts are reported. Neutralization removes eligible existing counterions; it does not add ions or set a new salt concentration.

Cholesterol is treated as a distinct molecule only in cholesterol-restoration mode, which is the mode that reasons about it. In protein, ligand and CHARMM-GUI transplant mode it is classified, coloured and reported as an ordinary membrane lipid.

A successful run produces GRO/PDB coordinates, an inspection PDB, text/JSON reports and `index.ndx`, plus `topol.top` and `toppar/` when topology writing is on (as in the examples). Review the alignment, removed molecules, retained clashes, lipid counts, net charge, and coordinate/topology checks before continuing.

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

## Runnable protein insertion example

The repository includes a complete input system in
[`examples/protein_insertion`](examples/protein_insertion/README.md).
Its relative paths are resolved against the installed example folder, so it runs
from anywhere:

```bash
python -m pip install .
gmxtransplant --mode protein \
  -i "$(python3 -c 'from gmxtransplant.gui.model import example_path; print(example_path("protein"))')" \
  --output ~/gmxtransplant-example
```

`environment/` supplies the target box, membrane, solvent and ions.
`replacement/` supplies the selected incoming protein and bound ligands only.
Protein topology definitions come from the replacement; retained environment
molecules use environment definitions. Explicit incoming ligand ITPs supersede
old copies in either directory. Conflicts within a source remain errors.

In the generated examples, incoming ligand names and ITP paths are together
near the top. `charge: from_itp` and `moleculetype` need not be repeated.
Coordinate `resname` and ITP `[ moleculetype ]` names need not be identical:
the unique matching `[ atoms ]` residue name determines the type. Ambiguous
ITPs still require an explicit molecule type. Old detailed YAMLs remain supported.

For protein-mode PDB name restoration, the environment's `topol.top` is found
beside `environment_toppar_dir`; no separate reference topology is needed.
A nonstandard layout can still set `name_restoration.reference_topol` explicitly.
The example enables name restoration so truncated POP/TIP/CHL names are
classified correctly before clash removal. Scientific choices (masks, cutoffs,
protection and composition) remain explicit in the YAMLs; inactive/default
plumbing is omitted. Standalone minimization settings remain in minimization.yaml.
