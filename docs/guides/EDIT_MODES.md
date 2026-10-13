# System-edit modes: salt, protonate, mutate, hmr, resize, flood, oligomer

These seven modes each change one prepared system (the host) and write a
run-ready folder: `step5_input.gro/.pdb`, `topol.top`, `toppar/`, `index.ndx`
with the CHARMM-GUI groups, a report, PyMOL/VMD scenes and, when asked, a
GROMACS equilibration folder. The host is a CHARMM-GUI-style folder
(`topol.top`, `toppar/`, `step5_input.gro`; `host_gro` picks another GRO, such
as an MD frame). Nothing outside the output folder is changed.

In the desktop application these modes are under **More modes** on the
Configuration page, and their examples under **More examples** on the Examples
page. Every mode has a D1 receptor example in `examples/<mode>/`, built on the
apo D1R host of the addbinder example.

| Mode | What it does | Example |
|---|---|---|
| `salt` | Sets the salt type and concentration (NaCl, KCl, CaCl₂, MgCl₂ …) | 150 mM KCl instead of NaCl |
| `protonate` | Changes protonation states (HSD/HSE/HSP, ASP/ASPP, GLU/GLUP, LYS/LSN, CYS/CYM) | D2.50 charged, H164 → HSP |
| `mutate` | Point mutations with a CHARMM36 topology and a clash-free rotamer | D3.32A, S5.42A |
| `hmr` | Hydrogen mass repartitioning for 4 fs time steps | 3.024 u hydrogens |
| `resize` | Makes the water layer thicker with copies of the host's bulk water | +20 Å |
| `flood` | Adds many copies of a small molecule to the water (ligand flooding, cosolvent probes) | 10 dopamine |
| `oligomer` | Adds copies of the membrane protein to the same membrane | D1R dimer lattice |

```bash
gmxtransplant --show-example salt > salt.yaml      # set host (and the mode's settings)
gmxtransplant --mode salt -i salt.yaml --dry-run
gmxtransplant --mode salt -i salt.yaml --output results/
```

Every input path must be given in full (the bundled examples are the
exception). The output goes to `<output>/<mode>_output/` unless `output_dir`
is set.

## Settings every mode shares

| Key | Meaning |
|---|---|
| `host` | The prepared system folder (required). |
| `host_gro` | Optional coordinates instead of the folder's `step5_input.gro`. |
| `environment_overrides` | Molecule types the automatic classifier does not know, e.g. `{MYLIP: lipid}`. |
| `neutralize`, `target_net_charge` | After the change, counterions far from the protein are removed or added (replacing bulk waters) to reach the target charge (default 0). The host's own monovalent ions are used. |
| `equilibration` | `enabled: true` writes `equilibration/` (below). Off by default; the examples turn it on. |
| `output_dir` | Output folder name or path. |

To replace the membrane (pure POPC, or any equilibrated bilayer from a GRO
file), use the [popc mode](POPC.md).

## salt

`concentration` (mol/L of formula units), `cation` and `anion` (moleculetype
names: CHARMM36 `SOD POT LIT CES RUB CAL MG BAR ZN2 CD2`, `CLA`). With
`replace_existing: true` the host's bulk ions are removed first;
`keep_bound_distance` (Å) keeps ions bound to the protein, e.g. a structural
Na⁺. New ions replace random bulk waters at least `solute_distance` from the
protein and ligands, `membrane_distance` from the lipids and `ion_spacing`
apart. The concentration is units / waters × 55.51, counting the waters left
after replacement; CaCl₂ counts once per Ca²⁺. An ion type the host lacks comes
from `ion_itps` or the bundled CHARMM36 ions, with its atom type and NBFIX
pairs added to forcefield.itp; nothing the host defines is changed.

## protonate

```yaml
protonate:
  host: /path/to/system
  residues:
    - {molecule: PROA, resid: 164, state: HSP}
    - {molecule: PROA, resid: 70, state: ASP}
```

`resid` is the residue number in the protein ITP. A residue's present state is
read from its atoms, not its name (CHARMM-GUI names a protonated aspartate ASP,
with HD2). States: `HSD HSE HSP ASP ASPP GLU GLUP LYS LSN CYS CYM`; a change
must stay within one residue (use mutate otherwise).

The residue's atoms get the CHARMM36 types and charges of the new state; every
bonded term that touches the residue beyond its backbone N and C is regenerated
the way pdb2gmx does for CHARMM36: bonds, impropers and CMAP from the residue
template, angles, proper dihedrals and 1-4 pairs from the bond graph. This
reproduces CHARMM-GUI's own ITP exactly for every standard residue (tested).
New hydrogens are placed from the CHARMM internal coordinates (within 0.1 Å of
CHARMM-built positions); position restraints follow CHARMM-GUI (POSRES_FC_BB on
N, CA, C, O; POSRES_FC_SC on other heavy atoms). Parameters the new state needs
and the host forcefield.itp lacks (e.g. the HSP atom types) are added from the
bundled CHARMM36 library. Every copy of an edited moleculetype is changed.
Terminal, capped and disulfide-bonded residues are refused. Water within
`water_clash_distance` (2 Å) of a new atom is removed.

## mutate

```yaml
mutate:
  host: /path/to/system
  mutations:
    - {molecule: PROA, resid: 103, to: ALA}
```

`to` is a CHARMM36 residue (`ALA ARG ASN ASP ASPP CYS GLN GLU GLUP GLY HSD HSE
HSP ILE LEU LYS LSN MET PHE PRO SER THR TRP TYR VAL`). The backbone (N, HN,
CA, HA, C, O) keeps its coordinates, and CB when both residues have one. The
side chain is built from the CHARMM internal coordinates; its rotamer is chosen
from staggered chi angles (−60°, 180°, 60°; for planar groups ±90°, 0°, 180°) by
the smallest steric overlap with the rest of the system, water and ions
excluded. The report gives the chi angles and the score. The topology is
rebuilt as for protonate. Water and ions within `water_clash_distance` (2.5 Å)
of a new atom are removed; lipids within `lipid_clash_distance` (2 Å) are
removed whole (`remove_clashing_lipids: false` only warns). Minimize before MD.
PDBFixer side chains are planned ([roadmap](PDBFIXER_ROADMAP.md)).

## hmr

Every hydrogen bonded to a heavy atom gets `hydrogen_mass` (3.024 u), taken
from that heavy atom, so each molecule keeps its mass; rigid water is left
unchanged (`skip_water`). With `copy_run_inputs` the host's mdp files are copied
and the production stage (step7) gets `time_step` (0.004 ps) with the same
simulated length. Use 4 fs only with `constraints = h-bonds`.

## resize

`add_water` (Å) is added to the box along z. The water above and below a
membrane is one layer joined through the periodic boundary, so the new water
goes in its middle: the molecules move so the bilayer is centred in the taller
box, and the gap is filled with copies of the host's own bulk water slab (water
and ions, at the host's density and salt), taken at least `slab_margin` (4 Å)
from the protein and lipids. Copies within `water_clash_distance` (2.4 Å) of
another molecule are dropped at the seams, so the new layer starts about 10%
below bulk density; the NPT equilibration closes the gap. Copied ions are
dropped as needed to keep the charge.

## flood

`molecule_coordinates`, `molecule_itp` and `molecule_forcefield` (the
forcefield.itp generated with a CGenFF ligand) describe the molecule; `copies`
(or `concentration` in mol/L of the water) how many. Each copy gets a random
position and orientation in the water (`side: both | upper | lower`), at least
`min_solute_distance` (10 Å) from the protein and ligands, `min_membrane_distance`
(5 Å) from the lipids and `min_copy_distance` (8 Å) from other copies, periodic
images included. Water and ions within `water_clash_distance` of a copy are
removed. The index group `Flood_<name>` holds every copy.

## oligomer

```yaml
oligomer:
  host: /path/to/system
  copies:
    - {shift: [64.76, 64.76], rotation: 180.0, flip: false}
```

Each copy of the host's proteins (and their ligands, `include_ligands`) is
turned by `rotation` about the membrane normal through the protein centre,
optionally flipped upside down about an in-plane axis through the membrane
midplane, and shifted by `shift` (x, y in Å). Every copy must stay
`min_protein_gap` (10 Å) from every other copy and all periodic images, or the
run stops with the distance found. Lipids within `lipid_clash_distance` (3 Å) of
a copy are removed whole, water and ions within `water_clash_distance`. The
index groups `Original` and `Copies` separate them. In a small box, try a
diagonal shift (half the box in x and y) before a shift along one axis.

## The equilibration folder

`equilibration/` holds the system, `topol.top`, `toppar/`, `index.ndx`, one mdp
per stage and `sample_run.sh`. The default stages are the CHARMM-GUI Membrane
Builder step6 schedule (minimization, two NVT and four NPT stages, 1.875 ns),
releasing the backbone, side-chain and lipid restraints step by step; stages
can be listed in `equilibration.stages`. Run it with

```bash
cd results/salt_output/equilibration
GMX="gmx" ./sample_run.sh
python3 analyze_equilibration.py
```

`analyze_equilibration.py` (numpy only) reads every finished stage and prints
the box, the projected area per phospholipid of each leaflet, the P–P
thickness, the protein C-alpha RMSD from the input (best fit, periodic images
taken into account) and the drift of the protein centre along z.

On a cluster, `gmxtransplant slurm-script` writes `submit.sh` next to
`sample_run.sh`:

```bash
gmxtransplant slurm-script results/salt_output/equilibration --profile dardel --account naiss2025-3-21
sbatch results/salt_output/equilibration/submit.sh
```

Profiles: `generic`, `dardel` (PDC modules, `srun gmx_mpi`, 128 ranks) and
`gpu`; `--partition`, `--module`, `--gmx`, `--mdrun`, `--nodes`,
`--ntasks-per-node`, `--cpus-per-task`, `--time` and `--gpus` override them.
Nothing is submitted.

## The GROMACS check

After every successful run, every mode checks what it wrote with
`gmx grompp -maxwarn 0` (a plain minimization mdp, no restraints) when GROMACS
is found: `GMXTRANSPLANT_GMX`, `GMX`, then `gmx`, `gmx_mpi` on PATH.
`--gromacs-check always` requires it, `never` skips it. A rejected system keeps
its files and the run exits with code 8. Any folder can be checked later:

```bash
gmxtransplant check-gromacs results/salt_output results/mutate_output
```

## Validation

All seven examples (and both popc examples) were built on Dardel (GROMACS 2026.1) with
`--gromacs-check always`: grompp passed each with 0 warnings and 0 notes; the
step6.0 minimization and 5 ps of the first restrained NVT stage ran without
LINCS warnings, and the HMR system ran 10 ps at 4 fs without LINCS warnings.
