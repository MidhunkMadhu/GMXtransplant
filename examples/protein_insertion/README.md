# Protein and bound-ligand insertion

Choose one of the following run options. Run from this example directory
because the YAML uses relative input and output paths. With the installed
package, you can also copy this directory and its inputs to another location
and run there; the source checkout is not required.

**Installed command** (after `python -m pip install .` or
`python -m pip install -e .` from the repository root containing `pyproject.toml`):

```bash
gmxtransplant --mode protein -i protein_replace.yaml --dry-run
gmxtransplant --mode protein -i protein_replace.yaml
```

**Installed Python module** (uses the package installed in the selected Python environment):

```bash
python -B -m gmxtransplant --mode protein -i protein_replace.yaml --dry-run
python -B -m gmxtransplant --mode protein -i protein_replace.yaml
```

**Source checkout** (run from this example directory inside the checkout,
with the dependencies installed):

```bash
python -B ../../src/run_pipeline.py --mode protein -i protein_replace.yaml --dry-run
python -B ../../src/run_pipeline.py --mode protein -i protein_replace.yaml
```

Dry-run checks configuration only; full validation requires your files.
`-B` prevents Python from creating `__pycache__` bytecode caches.

- `environment/frames_100ns.pdb`: target environment; its protein/ligands at
  sequential residue positions 1–963 are removed. Its box, membrane, water and ions survive.
- `replacement/step5_input.gro`: incoming protein and LDP/G4C at positions 1–963.
  Its surrounding environment is not imported. Alignment uses protein backbone
  atoms at positions 1–961 in both files.
- Each directory includes its matching `topol.top` and `toppar/` input files.
- `protein_replace.yaml`: runnable configuration, retaining the supplied masks,
  alignment, cutoffs and relative input paths, with simplified metadata and
  PDB name restoration enabled. Original folders were renamed to environment/replacement.

PDB POP is resolved by environment topology and atom counts, not guessed.
Without restoration, these lipids would be classified as “other” and lose lipid protection.

Validation of the bundled system: 220,306 final atoms, all 397 lipids retained,
33 clashing waters removed and one counterion removed; final audited charge 0 e.
The final PDB/GRO checks report three close contacts allowed by the configured
lipid protection/class cutoffs. Inspect the report before subsequent simulation;
this assembly check is not an energy minimization or dynamics validation.

Outputs are written here and ignored by Git. The committed example consists of
input files and configuration, so another user can reproduce the run.
