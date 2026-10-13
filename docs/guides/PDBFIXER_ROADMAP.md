# PDBFixer pipelines: roadmap

Status (2026-10-13): planned, not implemented. PDBFixer (OpenMM project, MIT
licence, `conda install -c conda-forge pdbfixer`) would be an optional
dependency, like OpenMM for minimization: `pip install ".[pdbfixer]"`. Every
mode keeps working without it.

## Why

GMXtransplant edits systems that are already built (CHARMM-GUI folders). PDBFixer
repairs and changes *structures*: missing atoms and loops, nonstandard
residues, mutations, hydrogens at a given pH. Used before or inside our modes,
it removes the main reasons a structure cannot go through them.

## The rule for every pipeline

PDBFixer works in PDB/OpenMM naming. GMXtransplant keeps its own topology step
(`protein_edit.py`, the bundled CHARMM36 library), so a PDBFixer result is only
a source of **heavy-atom coordinates**:

1. PDBFixer builds or changes the heavy atoms.
2. Names are mapped to CHARMM36 (HIS → HSD/HSE/HSP, ILE CD1 → CD, terminal and
   cap names); PDBFixer's hydrogens are dropped.
3. Hydrogens are rebuilt from the CHARMM internal coordinates, and the ITP from
   the CHARMM36 templates, exactly as now.
4. grompp checks the result (`--gromacs-check`).

## Pipelines, in order

### 1. mutate: PDBFixer side chains (`mutate.side_chains: pdbfixer`)
- `PDBFixer.applyMutations()` on the protein extracted from the host, then
  `addMissingAtoms()` builds the new side chain from its residue templates.
- Keeps the current clash scan as a check and as the fallback when PDBFixer is
  not installed (`side_chains: ic_scan`, the default today).
- Allows insertions/deletions later (see 4).

### 2. protonate: states from pH (`protonate.ph: 7.0`)
- `Modeller.addHydrogens(pH=...)` picks HIS/ASP/GLU/LYS/CYS variants from
  residue pKa values and hydrogen-bond patterns; map them to HSD/HSE/HSP, ASPP,
  GLUP, LSN, CYM and apply them with the existing protonate machinery.
- Report every residue whose state changes; explicit `residues:` entries win.
- Later: PROPKA for environment-dependent pKa values.

### 3. prepare: fix an incoming structure before protein insertion or addbinder
- New mode or a step of `protein`/`addbinder`: `findMissingResidues()`,
  `findNonstandardResidues()` / `replaceNonstandardResidues()` (MSE → MET …),
  `findMissingAtoms()` / `addMissingAtoms()`, `removeHeterogens(keepWater=False)`.
- Short missing loops are built; long ones only on request (`max_loop_length`),
  with a warning that built loops need equilibration.
- The output is a cleaned PDB that feeds CHARMM-GUI or the protein mode.

### 4. Loops and termini in mutate/protonate
- Build missing loop residues in a host protein (topology through the
  templates, with peptide bonds and CMAP).
- Change terminal residues and caps (NTER/CTER, ACE/CT3 patches), which are
  refused today.

### 5. Comparison against the current builders
- For every standard residue: PDBFixer side chain vs internal-coordinate build
  vs crystal structures (heavy-atom RMSD, clashes after minimization).
- Keep the bundled topology test: the ITP must still match CHARMM-GUI exactly.

## Not planned

- `Modeller.addMembrane()` / `addSolvent()`: CHARMM-GUI and the membrane, popc and
  resize modes already cover this with GROMACS topologies.
- Using OpenMM force-field XML topologies in place of the CHARMM36 ITPs.
