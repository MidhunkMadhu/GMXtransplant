# GMXtransplant

GMXtransplant builds a new simulation-ready coordinate/topology bundle by
inserting molecular components into an existing prepared environment. Its main
use is placing a receptor, a receptor–ligand complex, a replacement ligand, or
experimentally resolved cholesterol into an equilibrated membrane system while
preserving as much of that environment as possible.

The pipeline removes complete environment molecules that clash with the inserted
atoms, corrects charge when requested, regenerates the GROMACS molecule list,
writes an index, and produces detailed text and JSON reports. It can optionally
energy-minimize the assembled system with the free, open-source OpenMM Python
API. It validates inputs and reads written coordinate files back before accepting
them.

This is a structure-assembly tool and expects every molecular component to be
parameterized already. After assembly, the user can inspect and use the generated
bundle directly, request optional OpenMM minimization in the same command, or run
minimization later as a standalone step before continuing with any preferred
downstream protocol.

## What it can do

| Mode | Use case | Example configuration |
|---|---|---|
| `receptor` | Insert a receptor or receptor–ligand block into a prepared membrane/solvent environment | `receptor_replace.yaml` |
| `lig` | Replace one bound ligand while keeping the receptor and environment | `ligand_replace.yaml` |
| `chl` | Restore experimental cholesterol and correct leaflet composition | `cholesterol_restore.yaml` |
| `--minimize` | Minimize an assembled or existing GROMACS bundle with OpenMM L-BFGS | `minimization.yaml` |

In receptor mode, ligands selected together with the incoming receptor are moved
by the same rigid transform and inserted as one block. This is the intended route
for inserting an intact receptor–ligand complex into an equilibrated membrane.

## What can be used as the environment?

The environment coordinates do not have to come from CHARMM-GUI. They may be:

- a CHARMM-GUI GROMACS output such as `step5_input.gro`;
- a prepared system from another builder;
- a frame exported from an MD trajectory; or
- any other PDB/GRO structure that meets the requirements below.

An MD-generated frame is therefore a valid environment and is often useful when
the goal is to insert a receptor–ligand system into an already equilibrated
membrane. Export the frame using the topology that belongs to that trajectory so
atom order, residue names, whole molecules, and box vectors remain consistent.
A GRO frame is preferred because it carries the periodic box and four-character
GROMACS residue names reliably.

Regardless of origin, the environment must have matching GROMACS topology data
in the format expected by this project:

- a readable `topol.top` with a final `[ molecules ]` section;
- the referenced force-field and molecule `.itp` files, normally under
  `toppar/`;
- one unambiguous `[ moleculetype ]`/`[ atoms ]` definition for every retained
  receptor chain, lipid, ion, water, ligand, and cofactor; and
- coordinate atom names and atom order matching those definitions.

The incoming receptor or ligand must likewise have its matching GROMACS
parameters. The generated bundle is only as physically valid as those supplied
parameters.

For periodic processing, the environment needs an authoritative unit cell. Use a
box-bearing GRO/PDB or provide exactly `[a, b, c, alpha, beta, gamma]` in
Å/degrees in the YAML. A box is never guessed from a nearby file. Only
axis-aligned orthorhombic cells are currently accepted by the preflight validator.

## Installation

Python 3.10 or newer is required. Installing the project installs NumPy, SciPy,
PyYAML, and MDAnalysis automatically.

### Downloaded ZIP or source directory

```bash
cd /path/to/GMXtransplant
python3 -m pip install .
gmxtransplant --help
```

To enable energy minimization, install the optional OpenMM dependency:

```bash
python3 -m pip install '.[minimize]'
python3 -m openmm.testInstallation
```

For an editable developer installation:

```bash
python3 -m pip install -e .
```

You can also install directly from a Git repository after replacing the example
URL with the repository’s real URL:

```bash
python3 -m pip install 'git+https://github.com/OWNER/REPOSITORY.git'
```

After any installation, create a local editable configuration without finding
the installed data directory:

```bash
gmxtransplant --show-example receptor > receptor_replace.yaml
# Alternatives: --show-example lig, --show-example chl, or --show-example minimize
```

If `pip` is not available in the active Python environment, create a clean
environment first. Conda is recommended when cholesterol conversion or PyMOL
alignment is needed because those tools include non-Python components:

```bash
conda create -n gmxtransplant -c conda-forge \
  python=3.11 pip numpy scipy pyyaml mdanalysis openbabel
conda activate gmxtransplant
python -m pip install /path/to/GMXtransplant
conda install -c conda-forge openmm  # Only when --minimize is wanted.
```

Optional tools:

- **Open Babel** (`obabel`) is required only when
  `cholesterol.convert_to_charmm36: true`.
- **PyMOL** with the `pymol2` Python module is required only for
  `pymol_align` and `pymol_cealign`. The default examples use `mask_fit` and do
  not need PyMOL.
- **OpenMM** is required only for `--minimize`. It is an optional dependency,
  not imported by assembly-only runs.

The old script form remains available from a source checkout:

```bash
python3 run_pipeline.py --help
```

## Quick start

Work in a separate job directory so generated files cannot clutter the source
checkout:

```bash
mkdir my_gmxtransplant_job
cd my_gmxtransplant_job
gmxtransplant --show-example receptor > receptor_replace.yaml
```

Open the copied YAML. All file paths are grouped under `paths:` at the top.
Write each path once as `name: path`; use `"${name}"` in the settings below
to reference it. GMXtransplant resolves these references from `paths:` (including
`null` values). References must occupy the whole value; environment-variable
expansion and references inside the paths block are not supported. Existing YAML
anchors and aliases remain supported. Then adjust masks, residue names, and the
scientific choices for the new system.

Generate either version with these options (comments are included by default):

```bash
gmxtransplant --show-example receptor --comments > receptor_replace.yaml
gmxtransplant --show-example receptor --no-comments > receptor_replace_simple.yaml
```

Both versions contain the same active settings. These options also work with
`lig`, `chl`, and `minimize` examples.

The receptor example includes two ligands: protein positions 1–961, LIG1 at 962,
and LIG2 at 963. Replace these placeholders with your system's sequential residue
positions. For one ligand, remove the LIG2 metadata and path and change the incoming
mask to `":1-962"`. To add more, add one metadata entry and ITP path per ligand and
extend the incoming mask. Set the target mask independently to the old protein and
ligands you want removed; keep the alignment masks on corresponding protein atoms.

Validate the YAML schema and values without opening large coordinate/topology
files:

```bash
gmxtransplant --mode receptor -i receptor_replace.yaml --dry-run
```

Run the pipeline:

```bash
gmxtransplant --mode receptor -i receptor_replace.yaml
```

Run the assembly and then minimize its validated GRO/topology outputs:

```bash
gmxtransplant --mode receptor -i receptor_replace.yaml --minimize
```

The same flag works with `lig` and `chl`. Alternatively set
`minimization.enabled: true` in that mode's YAML.

Use `--mode lig` with `ligand_replace.yaml` or `--mode chl` with
`cholesterol_restore.yaml`. The mode and file must agree. Unknown and duplicate
YAML keys are rejected.

Relative paths are resolved from the directory where the command is run, not
from the installed package. Running from the job directory is therefore the most
predictable approach. Absolute paths work from any directory.

## Configurations

The included YAML files preserve the current working options. Alternative valid
choices are shown as comments beside the active values.

### Receptor/receptor–ligand insertion

The important inputs are:

- `target_box`: complete prepared system that supplies the retained receptor
  environment;
- `replacement_structure`: coordinates containing the incoming receptor block
  and any bound ligands selected with it;
- `alignment`: fit method and corresponding masks;
- `replacement_ligands`: charge and ITP metadata for ligands already inside the
  incoming selection; and
- `topology`: receptor and environment GROMACS topology sources.

`target_box.receptor_mask` is removed from the target. Every complete residue in
`replacement_structure.receptor_mask` is transformed and inserted. To include a
bound ligand, include its residue in that replacement mask and list it under
`replacement_ligands`. To omit it, omit it from the mask.

`mask_fit` requires equal atom counts, the same atom-name sequence, the same
per-residue grouping, and at least three non-collinear points. `pymol_align` and
`pymol_cealign` can determine correspondence through PyMOL.

When PDB input has truncated names such as `POP`, `TIP`, or `CHL`, enable
`name_restoration`. A single candidate is renamed directly; multiple candidates
are distinguished by exact ITP `[ atoms ]` counts. Ambiguity stops the run.

### Ligand replacement

`ligand_replace.structure_path` is one complete prepared system containing the
old ligand. `original_ligand.resname` selects the residue to remove. If the name
occurs more than once, use a one-based occurrence such as `UJU:2`.

The placement methods are:

- `autofit`: match every uniquely named heavy atom and perform a Kabsch fit;
- `pairfit`: fit explicitly paired atom-name lists; or
- `nofit`: preserve incoming coordinates exactly because the ligand is already
  in the target coordinate frame.

The receptor is never removed in this mode. Protein–ligand contacts are reported,
but protein residues are protected automatically.

### Cholesterol restoration

The experimental receptor is aligned to the receptor in the target environment;
only the experimental cholesterol coordinates are transferred. With
`convert_to_charmm36: true`, existing cholesterol hydrogen atoms are removed,
hydrogens are regenerated with Open Babel, and the result is mapped and checked
against the embedded 74-atom CHARMM36 cholesterol reference.

Missing heavy-atom reconstruction is conservative and limited by
`max_missing_heavy_atoms` and `max_heavy_atom_fit_rmsd`. Set
`write_diagnostics: true` to retain the conversion and placement intermediates.

`composition.lipid_targets` gives desired upper/lower leaflet counts. The mode
can remove surplus molecules outside `distance_from_protein`; it does not create
missing lipids. Salt is compared with `composition.reference_system_path`, which
is treated as the authoritative concentration baseline.

### Optional OpenMM energy minimization

Post-assembly minimization always uses the GRO just written by that run and the
generated `topol.top`; values in `minimization.coordinates_path` and
`minimization.topology_path` are ignored in this combined workflow. Consequently,
`topology.enabled` must be `true`. This prevents accidentally minimizing a new
coordinate file against an older topology.

An existing matching GROMACS bundle can instead be minimized independently:

```bash
gmxtransplant --show-example minimize > minimization.yaml
gmxtransplant --minimize -i minimization.yaml --dry-run
gmxtransplant --minimize -i minimization.yaml
```

With no YAML, standalone mode deliberately uses conservative filename defaults
in the current directory:

```bash
gmxtransplant --minimize
# Inputs: step5_input.gro and topol.top
# Outputs: minimized.gro and minimization_report.{txt,json}
```

Atom identity is checked against the original, preprocessed GROMACS `[ atoms ]`
names in `[ molecules ]` order. OpenMM's public topology normalizes some names
(for example `SER:HN` becomes `SER:H`, and water names can also change). Comparing
GRO names directly with that normalized topology caused false atom-order errors
in earlier versions. No coordinate or ITP renaming is needed for this case;
genuine source-name/count/order mismatches still stop minimization. The check
uses OpenMM's parsed molecule records and reports an explicit compatibility
error if a future OpenMM version no longer exposes them. After updating the
package, an already assembled bundle can be retried with the standalone command
above without repeating receptor replacement.

OpenMM uses constrained L-BFGS. `tolerance_kj_mol_nm` is the RMS objective-force
tolerance; `max_iterations: 0` asks OpenMM to continue until convergence. The
finite default avoids an accidentally unbounded command and reports a nonzero
exit status if the objective remains above tolerance.

A GROMACS topology does not contain every run choice normally held in an MDP.
The YAML therefore makes the nonbonded method, cutoff, LJ switch distance,
constraints, rigid-water choice, Ewald tolerance, and dispersion correction
explicit. Match these values to the force-field protocol that produced the
system. OpenMM's switching function should not be assumed numerically identical
to every GROMACS `vdw-modifier` variant.

`platform: auto` tries CUDA, HIP, OpenCL, CPU, then Reference, falling back when
a registered accelerator cannot create a context. A selected GPU can be set with
`device_index`. `defines: {POSRES: 1}` activates compatible conditional position
restraints already present in the supplied topology; it does not create new
restraints.

The original assembled GRO is never overwritten. The minimized GRO and both
reports are staged, read back, and promoted together. Reports include initial
and final potential energy, physical-force diagnostics, the minimizer objective,
convergence status, OpenMM platform/settings, and any accelerator fallbacks.

## Residue masks

Masks use one-based sequential residue positions in coordinate-file order. They
do not use the residue number printed in PDB/GRO.

| Mask | Meaning |
|---|---|
| `:1-298` | residues 1 through 298 |
| `:1,5,10-20` | residues 1, 5, and 10 through 20 |
| `:1-298@CA` | CA atoms in residues 1 through 298 |
| `:1-298@C,CA,N,O` | backbone atoms in residues 1 through 298 |
| `:1-298&!:297-298` | residues 1 through 296 |
| `@CA` | every CA atom |

Unsupported cpptraj syntax is rejected rather than simplified silently.

## Clash, box, and charge safety

Environment preflight runs before alignment or coordinate output:

- `strict` validates the existing coordinate/cell frame;
- `repair` may apply one validated whole-system rotation and translation; and
- `off` is an explicit opt-out required for nonperiodic clash detection.

The validator checks finite coordinates, the full cell, direct and periodic hard
contacts, and output-format rounding. Repair never rotates molecules separately.
A failed bounded orientation search means the frame remains unresolved; it does
not prove that the supplied box is historically wrong.

Replacement clash detection compares each environment residue with the inserted
block using the configured class cutoff and minimum-image distance. When one atom
pair crosses the cutoff, the complete environment residue is removed. Existing
environment–environment contacts are left unchanged. Classes in `keep_classes`
are retained but their contacts remain visible in the reports.

Charge accounting prefers ITP `[ atoms ]` charges. Explicit overrides are for
species without a usable ITP and cannot contradict one. Neutralization removes
eligible counterions only when doing so moves the system toward
`target_net_charge`; candidates can be excluded near receptor, ligand, and lipid
heavy atoms.

## Outputs

With the enabled options in the examples, a successful run produces:

- `step5_input.gro` and `step5_input.pdb`;
- an alignment/placement inspection PDB;
- `<report>.txt` and `<report>.json`;
- `topol.top` plus a self-contained `toppar/`; and
- `index.ndx`.

When minimization is requested, it additionally produces `minimized.gro`,
`minimization_report.txt`, and `minimization_report.json`.

The final order is protein, ligand/cofactor, lipid, ion, then water. PDB and GRO
are written to temporary files, read back, and checked for atom names/order,
coordinates, cell, final clashes, and—when topology output is enabled—exact
topology atom order/count and charge. Topology is built and audited in temporary
storage before replacing an existing `topol.top`/`toppar` bundle.

Review at least:

- alignment RMSD and the inspection PDB;
- every removed molecule and `flagged_but_kept` contact;
- lipid-removal fraction and final leaflet counts;
- unknown residue names or topology ambiguity;
- final net charge; and
- topology/coordinate audit status.

## Common problems

**“No usable box”** — use a GRO containing the correct frame box or set the six
explicit dimensions from an authoritative source. Do not copy a box from an
unrelated frame.

**Unknown or ambiguous molecule type** — supply the correct environment/receptor
`toppar/`, add a ligand ITP, restore truncated PDB residue names, or use a
specific `moleculetype_overrides` entry.

**Alignment selections differ** — compare the atom count and atom-name order of
the two masks. Use `pairfit` for a shared ligand core or PyMOL alignment for
non-identical receptor sequences.

**Open Babel not found** — install it with Conda (`conda install -c conda-forge
openbabel`) or set `cholesterol.obabel_command` to its executable path.

**OpenMM is required for --minimize** — install with `python -m pip install
'.[minimize]'` from a downloaded checkout, or `conda install -c
conda-forge openmm`, then run `python -m openmm.testInstallation`.

**OpenMM cannot parse the topology** — make sure every quoted include is inside
the generated bundle. Set `minimization.include_dir` only for an additional
include root, and provide required preprocessor values through
`minimization.defines`. Advanced GROMACS-specific topology constructs may still
require conversion or simplification before OpenMM can use them.

**Minimization exits with status 7** — the assembly outputs remain valid and
untouched. Read `minimization_report.json` if it was written; increase
`max_iterations` only after checking for bad parameters, severe overlaps, and
incorrect nonbonded settings.

**Large lipid-removal warning** — inspect the inserted pose, cell, and mask before
continuing. Keeping all lipids can hide serious overlaps; it is a deliberate
expert choice, not an automatic repair.

## Tests and release build

Run all tests from the repository root:

```bash
python3 -m unittest discover -v
```

Build the installable source archive and wheel:

```bash
python3 -m pip install build
python3 -m build
```

Artifacts are written under `dist/`. GitHub Actions runs the core suite on
Python 3.10–3.12, builds both distributions, and runs the real OpenMM
minimization smoke test in a separate optional-dependency job. Before publishing
the repository, add the owner’s chosen `LICENSE` file and update any repository
URL used in installation examples.

At handoff, users can keep the validated assembly unchanged, use the optional
OpenMM-minimized coordinates, apply a different compatible minimizer, or take
the standard GROMACS-format bundle into their chosen downstream protocol.

## Current limitations

- Environment coordinates are limited to PDB/GRO and orthorhombic cells for
  validated periodic processing.
- PDB has finite-width atom/residue fields; use GRO as the authoritative GROMACS
  coordinate output.
- Topology conditionals inside `[ moleculetype ]`, `[ atoms ]`, or `[ molecules ]`
  require a preprocessed topology.
- Cholesterol composition correction removes surplus molecules but cannot add
  missing lipids or solvent.
- Optional energy minimization uses OpenMM L-BFGS; users who prefer a different
  minimizer can take the validated GRO/topology bundle into another compatible
  workflow.
