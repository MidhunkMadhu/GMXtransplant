# Ligand replacement

Place your files in this layout (folders are already created):

```text
ligand_replacement/
├── ligand_replace.yaml
├── environment/
│   ├── step5_input.gro
│   ├── topol.top
│   └── toppar/         # all matching force-field and molecule ITPs
└── replacement/
    ├── LI2.mol2
    ├── LI2.itp
    └── forcefield.itp  # bonded parameters generated with LI2.itp
```

The environment is your complete prepared protein, old ligand, membrane,
water and ions. The replacement contains only the incoming ligand and its
matching parameters. Preserve all topology include dependencies.

This example replaces one bound ligand with another. The ligands are given the
generic residue names **LI1** (the ligand already in the prepared system) and
**LI2** (the incoming ligand, here phenethylamine); substitute your own names in
the YAML. The supplied system has one LI1 ligand after 298 protein residues. The
incoming MOL2 has generic atom labels; matching names and order are supplied by
`replacement/LI2.itp`, and the output residue is named LI2.

LI1 and LI2 are different molecules, and the LI2 coordinates start away from
the binding site. LI2 is therefore placed with explicit `pairfit` on the 9
heavy atoms the two share: LI1 contains the same amine–ethyl–benzene motif as
LI2. The pairs follow the bond graphs (amine N, the two chain carbons, then the
ring). Of the two ring orientations, the closer fit was selected (about 0.40 Å
RMSD). Name-based `autofit` is inappropriate because equal names in LI1 and LI2
do not denote corresponding atoms. Review the inspection PDB.

The pairs do not have to be worked out by hand. `fit.method: mcsfit` detects the
largest substructure the two ligands share automatically: same element on each
matched atom and identical bonds, read from the ITPs (or from interatomic
distances when no ITP is available). It compares every symmetric equivalent,
such as a ring matched in either direction, fits on the closest, and prints the
chosen pairs in the `old_ligand_fit_atoms`/`new_ligand_fit_atoms` format. For
this example it finds the same 9 pairs as the YAML (0.40 Å). The example keeps
the explicit `pairfit` lists so the placement is fixed and visible in the file;
to try the automatic detection, change `method: pairfit` to `method: mcsfit`.
A common substructure places the shared core only; for very different ligands,
inspect the pose.
Keep the supplied replacement ITP: its atom order and parameters differ from
the environment's LI1 ITP. It is included automatically in the output, and
LI1.itp is dropped because LI1 no longer appears in the system. The bond, angle
and dihedral parameters of LI2 come from `replacement/forcefield.itp`, selected
by `new_ligand.forcefield_path`.

`paths.system_toppar` is shared by the retained protein and environment.
It includes cholesterol parameters too; there is no need to duplicate it.

Choose one of the following run options. Run from this example directory
because the YAML uses relative input and output paths. With the installed
package, you can also copy this directory and its inputs to another location
and run there; the source checkout is not required.

**Installed command** (after `python -m pip install .` or
`python -m pip install -e .` from the repository root containing `pyproject.toml`):

```bash
gmxtransplant --mode lig -i ligand_replace.yaml --dry-run
gmxtransplant --mode lig -i ligand_replace.yaml
```

**Installed Python module** (uses the package installed in the selected Python environment):

```bash
python -B -m gmxtransplant --mode lig -i ligand_replace.yaml --dry-run
python -B -m gmxtransplant --mode lig -i ligand_replace.yaml
```

**Source checkout** (run from this example directory inside the checkout,
with the dependencies installed):

```bash
python -B ../../src/run_pipeline.py --mode lig -i ligand_replace.yaml --dry-run
python -B ../../src/run_pipeline.py --mode lig -i ligand_replace.yaml
```

Dry-run checks configuration only; full validation requires your files.
Outputs are written here, like the protein example: `step5_input.pdb`,
`step5_input.gro`, `topol.top`, `toppar/`, `index.ndx`, inspection PDB and
reports. Inputs remain in their separate folders. Generated outputs are
ignored by Git. Add `--prepare-minimization` only when you want a portable
minimization folder; it does not run minimization.

`-B` prevents Python from creating `__pycache__` bytecode caches.

Validated with the supplied files: 58,217 final atoms, 0.399 Å pairfit RMSD,
no clash removals or ion removals, and audited net charge 0 e. Final PDB and
GRO topology/atom-order and hard-clash checks pass, and `gmx grompp` accepts the
output topology. The retained protein
coordinates are unchanged within GRO rounding precision. Inspect the pose
before subsequent minimization or simulation; neither was run by this test.
