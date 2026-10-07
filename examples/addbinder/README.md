# addbinder examples: the apo dopamine D1 receptor

Two examples share one host, the apo D1R membrane build:

| Config | Binder | Side |
|---|---|---|
| `addbinder.yaml` | dopamine (ligand) | extracellular (`upper`) |
| `addbinder_gprotein.yaml` | Gs heterotrimer (protein) | intracellular (`lower`) |

From the source checkout root, open the example directory and run either
configuration:

```bash
cd examples/addbinder
python -B ../../run_pipeline.py --mode addbinder -i addbinder.yaml --output results/
python -B ../../run_pipeline.py --mode addbinder -i addbinder_gprotein.yaml --output results_gprotein/
```

The installed package also includes these inputs. In the desktop application,
both runs are on the Examples tab.
For your own system, use the [addbinder mode guide](../../ADDBINDER.md) to set
the host, binder coordinates, ITP files, force field, side, and poses.

## Inputs

| Folder | Contents | Where it came from |
|---|---|---|
| `host/` | `step5_input.gro`, `topol.top`, `toppar/`, `step6.*`/`step7` mdp files, `README` (run script) | CHARMM-GUI Membrane Builder: apo D1R (PROA) in POPC/POPE/cholesterol with Na⁺/Cl⁻. Its `forcefield.itp` is protein-only, so it has no CGenFF types. |
| `dopamine/` | `dopamine.pdb`, `LDP.itp`, `forcefield.itp` | Dopamine (CGenFF residue `LDP`, +1 e) from `examples/cholesterol_restoration/environment`. The PDB is the bound pose from that system. Its `forcefield.itp` holds the CGenFF parameters. |

The installed package leaves out `host/step5_input.pdb` (30 MB), which
addbinder does not read.

| `gprotein/` | `gs_trimer.pdb`, `PROB.itp`, `PROC.itp`, `PROD.itp`, `forcefield.itp` | The Gs heterotrimer (Gα PROB, Gβ PROC, Gγ PROD; 9,954 atoms, −7 e) from the D1R–Gs complex in `examples/cholesterol_restoration/environment`, placed in this host's frame in its receptor-bound arrangement (see below). |

You do not edit any `forcefield.itp` yourself. addbinder merges the host and
binder force fields into each pose's `toppar/forcefield.itp`, and it stops with
an error if the two files give one parameter different values. For this example
they merge cleanly: the CGenFF types dopamine needs (`CG2R61`, `NG3P3`, `OG311`,
`CG324`, `HGP2`, ...) are added to the receptor's CHARMM36m parameters.

# Example 1: dopamine, extracellular

## Poses

The D1R N-terminus points to +z, so `side: upper` is the extracellular side.
The tip is THR283 OG1, 15.3 Å above the upper headgroup plane. The free water
column is about 158 Å, so there is plenty of room.

| Pose | Placement (seed 1) | Distance |
|---|---|---|
| `binderpose1` | straight above the tip, flat | 20 Å beyond the tip |
| `binderpose2` | tilted 41° from the normal, azimuth 298°, random orientation | 24.75 Å from the receptor |
| `binderpose3` | tilted 50°, azimuth 100°, random orientation | 22.43 Å |
| `binderpose4` | tilted 55°, azimuth 221°, random orientation | 24.85 Å |

The random poses differ from each other and from `binderpose1` by at least 25° in
direction and 10 Å in position. If a draw comes closer than that, or fails a
box check, it is drawn again. Change `random_seed` for a different set, or
`random_poses` for more or fewer.

Other things you can set per pose in `addbinder.yaml`:

```yaml
  poses:
    - distance: 20.0                   # binderpose1, straight above
    - distance: 15.0                   # binderpose2, its own distance
    - approach: [30.0, 90.0]           # tilt from the normal, azimuth (degrees)
      distance: 22.0
    - centroid: [70.0, 50.0, 185.0]    # dopamine centre in A, host GRO frame (nm x 10)
      orientation: as_is
```

## Output: four folders, each ready to run

Each pose folder holds `step5_input.gro/.pdb`, `topol.top`, `toppar/`,
`index.ndx`, the host's mdp files and its CHARMM-GUI `README` run script. The
index has the CHARMM-GUI groups `SOLU MEMB SOLV SOLU_MEMB` that the mdp files
use, plus `Binder` and `Host_tip`. To start a run:

```bash
cd results/addbinder_output/binderpose1
csh README          # or run its grompp/mdrun lines by hand
```

During step6 (`-DPOSRES`), dopamine is position-restrained like the protein,
because `LDP.itp` has its own `POSRES` block. In step7 it is free.

To see all four poses at once, open `view.pse` or `view.pml` in PyMOL, or
`view.vmd` in VMD, or use the PyMOL/VMD buttons in the GUI. The tip atom is the
yellow sphere.

# Example 2: Gs heterotrimer, intracellular

The D1R C-terminus points to −z, so `side: lower` is the intracellular side.
**Force field.** CHARMM-GUI writes into `forcefield.itp` only the parameters
its own system uses. The D1R-only host therefore lacks some bonded terms that
Gs needs, e.g. the `CC-CT1` bond and angles of Gα's C-terminal carboxylate
(LEU 246) and some proline/aspartate backbone dihedrals. addbinder checks this
before building and would stop with those terms listed. The example therefore
sets `binder_forcefield: gprotein/forcefield.itp`, the force field of the
D1R–Gs complex the trimer came from; it merges with the host's without a
single conflicting value.

**Bound arrangement.** `gs_trimer.pdb` holds the trimer where it sits on the
receptor in the D1R–Gs complex, moved into this host's frame:

- the complex was rotated about z only, fitted on the 143 transmembrane-helix
  Cα atoms of the two D1R copies (2.55 Å RMSD: apo vs. G-protein-bound
  receptor);
- it was moved in z so that the two membranes' midplanes coincide.

So the G protein keeps its exact orientation relative to the membrane. In this
arrangement it slightly overlaps the apo receptor, whose intracellular cavity
is closed (0.2 Å closest gap).

**Poses.** All three start from the bound arrangement
(`from_input_position: true`) and move the trimer out until its nearest heavy
atom is the distance from the receptor:

| Pose | Direction | Orientation | Distance | Gap to membrane |
|---|---|---|---|---|
| `binderpose1` | straight down the membrane normal | bound | 20 Å | 43.9 Å |
| `binderpose2` | random: tilted 31°, azimuth 298° | bound | 15 Å | 34.3 Å |
| `binderpose3` | random: tilted 38°, azimuth 100° | random rotation | 15 Å | 15.2 Å |

- `binderpose1` asks for 20 Å with `reduce_distance_by: 2` and
  `min_distance: 10`: if 20 Å did not fit the box (10 Å to every periodic
  image, 5 Å to lipids), 18, 16 ... 10 Å would be tried. Here 20 Å fits.
- `binderpose2` and `binderpose3` are random poses
  (`random_from_input_position: true`): random side and tilt (up to 45°), at
  least 25° apart in direction and 10 Å in position from each other and from
  `binderpose1`. With `random_rotated: 1` only the last one, `binderpose3`, is
  also randomly rotated; `binderpose2` keeps the bound orientation.
- The trimer moves further than the distance (43 Å for 20 Å) because Gα's α5
  helix sits deep in the receptor in the bound arrangement.

**Water and ions.** The trimer replaces about 3,900 waters and the ions in
that space, and brings −7 e. addbinder then sets the ion counts so the salt
stays at the host's 0.1537 M on the remaining water and the system is neutral
(D1R +14, Gs −7: seven more Cl⁻ than Na⁺):

| Pose | Na⁺ / Cl⁻ after overlap removal | Removed | Added (in place of bulk water) | Final | Salt |
|---|---|---|---|---|---|
| `binderpose1` | 286 / 312 | 16 Cl⁻ | 3 Na⁺ | 289 / 296 | 0.1536 M |
| `binderpose2` | 293 / 306 | 4 Na⁺, 10 Cl⁻ | none | 289 / 296 | 0.1536 M |
| `binderpose3` | 289 / 303 | 7 Cl⁻ | none | 289 / 296 | 0.1536 M |

Ions are taken from, and added to, bulk water anywhere in the box (at least
10 Å from protein, 5 Å from lipids and other ions).
