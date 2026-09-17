# Protein and bound-ligand insertion

Run from this directory after installing the repository with `python -m pip install -e .`:

```bash
gmxtransplant --mode protein -i protein_replace.yaml
```

- `environment/frames_100ns.pdb`: target environment; its protein/ligands at
  sequential residue positions 1–963 are removed. Its box, membrane, water and ions survive.
- `replacement/step5_input.gro`: incoming protein and LDP/G4C at positions 1–963.
  Its surrounding environment is not imported. Alignment uses protein backbone
  atoms at positions 1–961 in both files.
- Each directory includes its matching `topol.top` and `toppar/` input files.
- `protein_replace.yaml`: runnable configuration, retaining the supplied masks,
  alignment, cutoffs and relative input paths, with simplified metadata and
  PDB name restoration enabled. Original folders were renamed to environment/replacement.
- `protein_replace.original.yaml`: exact original user configuration, for reference;
  its old folder names apply to the original job directory, not this copy.

PDB POP is resolved by environment topology and atom counts, not guessed.
Without restoration, these lipids would be classified as “other” and lose lipid protection.

Validation of the bundled system: 220,306 final atoms, all 397 lipids retained,
33 clashing waters removed and one counterion removed; final audited charge 0 e.
The final PDB/GRO checks report three close contacts allowed by the configured
lipid protection/class cutoffs. Inspect the report before subsequent simulation;
this assembly check is not an energy minimization or dynamics validation.

Outputs are written here and ignored by Git. The committed example consists of
input files and configuration, so another user can reproduce the run.
