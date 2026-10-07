# Cholesterol restoration

Place your files in this layout (folders are already created):

```text
cholesterol_restoration/
├── cholesterol_restore.yaml
├── environment/
│   ├── step5_input.gro
│   ├── topol.top
│   └── toppar/         # all matching force-field and molecule ITPs
└── experimental/
    └── D1R_chl_fixed.pdb
```

The environment is your complete prepared protein, membrane, water and ions.
The experimental PDB contains the protein and its resolved cholesterol.
Preserve all topology include dependencies and include matching CHL1 parameters.

Both protein masks and alignment masks select the matching 328-residue D1R
receptor in the supplied files. The target also contains the G-protein complex;
it is retained but does not define the alignment. The experimental file contains
five CHL1 residues. Masks count sequential positions, not printed residue IDs.

The YAML uses the fresh GRO directly, including its stored periodic box.
No environment PDB or explicit box override is required. This input contains
223,627 atoms, correctly grouped waters, and CHL1/POPE/POPC counts 83/130/210.

The updated source topology and ITPs match the fresh GRO, including HD2 in
protein residue 100 (ASP), 151 CLA and 52,519 TIP3. All 15,160 protein atoms
match the ITP atom order. Keep these updated parameters with this input.

There is only one parameter path, `paths.system_toppar`, used for both protein
and environment roles. This directory also supplies CHL1 parameters.
`composition.reference_system_path` is omitted: the target itself supplies
the baseline salt counts. A separate reference is useful only when deliberately
comparing against a different starting system. Leaflet targets are explicit
desired counts, not automatically taken from that reference. Review those
targets and the 0.8 Å removal cutoff for your intended composition.

Open Babel is required for the enabled CHARMM36 conversion. It is installed with
GMXtransplant on Linux and macOS 14+; elsewhere see the README's install section.
Use `convert_to_charmm36: false` only for complete, already compatible CHL1
coordinates. Composition correction removes surplus lipids but cannot add
missing ones. Salt concentration is compared against the environment,
assumed to be 0.15 M in this configuration.

Choose one of the following run options. Run from this example directory
because the YAML uses relative input and output paths. With the installed
package, you can also copy this directory and its inputs to another location
and run there; the source checkout is not required.

**Installed command** (after `python -m pip install .` or
`python -m pip install -e .` from the repository root containing `pyproject.toml`):

```bash
gmxtransplant --mode chl -i cholesterol_restore.yaml --dry-run
gmxtransplant --mode chl -i cholesterol_restore.yaml
```

**Installed Python module** (uses the package installed in the selected Python environment):

```bash
python -B -m gmxtransplant --mode chl -i cholesterol_restore.yaml --dry-run
python -B -m gmxtransplant --mode chl -i cholesterol_restore.yaml
```

**Source checkout** (run from this example directory inside the checkout,
with the dependencies installed):

```bash
python -B ../../src/run_pipeline.py --mode chl -i cholesterol_restore.yaml --dry-run
python -B ../../src/run_pipeline.py --mode chl -i cholesterol_restore.yaml
```

Dry-run checks configuration only; full validation requires your files.
Outputs are written here, like the protein example: `step5_input.pdb`,
`step5_input.gro`, `topol.top`, `toppar/`, `index.ndx` and reports.
The run also writes `repack/`: `topol.top`, `toppar/`, `step5_input.gro`,
`index.ndx` and mdp files for a 1.4 ns equilibration that holds the protein,
ligands and the five restored cholesterols while the other lipids repack around
them. A sample run script is also saved (see the main README, "Repacking
before MD").
The report opens with a summary: where each restored cholesterol is in the
final files (residue and atom numbers, its place among the CHL1 in
`topol.top`), the molecules removed next to it, and the molecules removed far
from the protein with their distances. The GUI shows the same in short.
Inspection/conversion files are written only with `write_diagnostics: true`.
Inputs remain in their separate folders. Generated outputs are ignored by
Git. Add `--prepare-minimization` only when you want a portable minimization
folder; it does not run minimization.

`-B` prevents Python from creating `__pycache__` bytecode caches.

For the fresh GRO, all five cholesterol conversions pass and 328 CA pairs
align at 0.153 Å RMSD. Insertion removes 7 clashing lipids and composition
adjustment removes another 24. The full run passes: 220,422 final atoms,
audited net charge 0 e, and zero final hard clashes in both PDB and GRO.
All 143 SOD and 151 CLA remain. Protein coordinates are preserved within
output precision. Final leaflet counts are CHL1 40/40, POPE 60/59 and
POPC 100/98. No minimization or simulation was run, and source molecular
inputs remain unchanged.
