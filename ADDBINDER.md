# addbinder mode: add a binder outside a membrane protein

`--mode addbinder` takes a prepared membrane system and adds one incoming
molecule (a ligand or a protein) in the water above or below the membrane
protein, at a set distance from the protein's outermost heavy atom. Typical
uses:
- an apo receptor with its ligand placed outside the pocket, ready for
  unbiased binding simulations
- a partner protein placed near the receptor's extracellular or intracellular face

```bash
gmxtransplant --show-example addbinder > addbinder.yaml   # edit the paths
gmxtransplant --mode addbinder -i addbinder.yaml --dry-run
gmxtransplant --mode addbinder -i addbinder.yaml --output results/
```

To try the included dopamine and G-protein examples, see the
[addbinder example guide](examples/addbinder/README.md). The desktop application
also has a card for each example on its Examples tab.

## Inputs

| Key | What to give |
|---|---|
| `host` | Folder holding the system the binder joins: `topol.top`, `toppar/` and `step5_input.gro` (or one `step5*.gro`; `host_gro` overrides it). A CHARMM-GUI membrane build fits directly. The GRO must match the topology atom for atom. |
| `binder_coordinates` | PDB, GRO or MOL2 of the binder, hydrogens included. Atom names must match the ITP(s); for a single molecule with unique names the order may differ. For a multi-molecule binder, the molecules appear in `binder_itp` order. Position and orientation do not matter except for `orientation: as_is`. |
| `binder_itp` | The binder's ITP with exactly one `[ moleculetype ]`, or a **list of ITPs, one molecule each**, for a binder made of several molecules (e.g. a G-protein heterotrimer `[PROB.itp, PROC.itp, PROD.itp]`). The molecules are placed as one rigid body. A name that clashes with a different host molecule type (CHARMM-GUI calls protein chains `PROA`, `PROB`, ...) is renamed `<name>_BND`. |
| `binder_forcefield` | The force field holding the binder's parameters, e.g. the `forcefield.itp` CHARMM-GUI Ligand Reader generated with the ITP. Parameters are merged with the host's; the same key with different values stops the run. Before building, every binder molecule (protein or ligand) is checked: each atom type must be in `[ atomtypes ]`, and each bond, angle, dihedral and CMAP term must have parameters, in the host's force field plus this one. A gap stops the run and names the terms and residues (e.g. a non-standard residue). Note that a CHARMM-GUI `forcefield.itp` holds only what its own system uses, so a protein from a different build usually needs its own. |

For your own configuration, give input paths in full. The supplied examples
resolve their relative input paths against the example folder.

## How placement works

1. **Membrane frame.** The two headgroup planes are the median z of the lipid
   phosphorus atoms in each leaflet. The bilayer must be flat and normal to z.
   `side: upper` is +z and `side: lower` is −z. A coordinate file does not say
   which side is extracellular, so you choose.
2. **Tip.** The host protein and its bound ligands are made whole across the
   periodic boundary. The tip is their outermost heavy atom on the chosen side.
   The report names it, e.g. `PROA HSD159 NE2`.
3. **Pose.** The binder is rotated (see below), centred over the tip in x/y
   (`lateral_reference: host_center` uses the protein's centre instead), shifted
   by `lateral_offset`, and moved along the normal until its nearest heavy atom
   is `distance` Å beyond the tip plane.
4. **Checks** (heavy atoms, all 26 periodic images). Any failure rejects the
   pose and explains why:
   - distance to the host ≥ `distance`
   - distance to any periodic image of the host ≥ `min_image_gap`; when the
     closest image is along z, the message says how much more box height is needed
   - distance to the binder's own images ≥ `min_image_gap`
   - distance to lipid heavy atoms ≥ `min_membrane_gap`
5. **Water and ions.** Waters and ions with any heavy atom within
   `water_clash_distance` (2.4 Å) of any binder atom are removed. Protein,
   ligand and lipid molecules are never removed.
6. **Salt and charge.** The binder displaces water and the ions in it, and
   brings its own charge. The final ion counts are therefore set directly:
   - the host's salt concentration is measured as salt pairs × 55.51 / waters,
     using its main monovalent cation and anion (e.g. SOD/CLA);
   - on the remaining water, pairs = concentration × waters / 55.51 (rounded);
   - the counterion species gets pairs + |net charge of everything else|, so
     the system reaches `target_net_charge`; the other species gets pairs.

   Extra ions are removed and missing ones replace bulk waters, chosen at
   random among molecules at least `ion_exclusion_distance` (10 Å) from the
   protein and binder, `lipid_ion_exclusion_distance` (5 Å) from lipids, and
   (for added ions) 5 Å from any other ion. By default this uses the whole box
   (`ion_side: any`): in a single-bilayer periodic box the two water layers are
   one connected layer, so the ions redistribute during equilibration.
   `ion_side: binder_side` uses the binder's side first. Set
   `salt_concentration` to a number (mol/L) to use a fixed concentration
   instead of the host's.
7. **Output.** Coordinates in topology order: host proteins, binder proteins,
   host ligands, binder ligands, environment. The binder ends up inside the box
   and whole. The written topology is re-read and audited against the
   coordinates, including the charge.

### Orientation

| `orientation` | Meaning |
|---|---|
| `auto` (default) | A different orientation for each pose, whether or not the poses share a distance. |
| `flat` | Principal axes: longest along x, thinnest along the membrane normal. This uses the least box height. |
| `end_on` | Longest axis along the normal. The end that reaches furthest from the centre (e.g. a ligand's tail rather than its ring) faces the protein. |
| `edge` | Middle axis along the normal, longest along x. |
| `as_is` | Keep the input orientation, e.g. a docked or pocket pose. |
| `euler` | Rotate by `angles: [z, y, z]` (degrees) about the binder centre. |

`flip: true` turns the binder upside down (180° about x). `spin` (degrees) then
turns it about the membrane normal.

`auto` poses take, in turn: flat, end_on, end_on flipped, edge, flat flipped,
edge flipped. It skips any of these that a pose already uses explicitly. After
six it starts again, turned 45° further about the normal. With
`same_orientation: true`, every `auto` pose gets the first one (flat). Give a
pose its own `distance` to vary the distance too.

### Where each pose goes

| Pose setting | Placement | `distance` means |
|---|---|---|
| (default) | straight out along the membrane normal, above (`upper`) or below (`lower`) the tip, shifted by `lateral_offset` | gap along the normal from the tip atom to the nearest binder heavy atom |
| `approach: [tilt, azimuth]` | out from the tip along a direction tilted `tilt`° from the normal towards `azimuth`° (0 = +x, 90 = +y) | gap from the nearest host heavy atom to the nearest binder heavy atom |
| `centroid: [x, y, z]` | binder heavy-atom centre exactly at x, y, z (Å, in the host GRO's coordinates, i.e. GRO nm × 10) | not used; the pose is rejected only if it overlaps the host (under 3 Å) or fails the box checks |
| `distance_to: centroid` (with optional `approach`) | the binder's heavy-atom centre at `distance` from the tip (or the host centre, `lateral_reference: host_center`) along the normal or the `approach` direction: spherical coordinates (r = `distance`, θ, φ = `approach`). Rotate with `orientation: euler` and `angles: [z, y, z]` | centre-to-reference distance; rejected only if it overlaps the host (under 3 Å) or fails the box checks |
| `from_input_position: true` | starts where the binder is in its input file (coordinates already in the host's frame, e.g. a G protein in its receptor-bound arrangement) and moves it straight out along the normal, or along `approach` | gap from the nearest host heavy atom to the nearest binder heavy atom |

Every pose can have its own `distance`. As many poses as you like can be
listed. With `reduce_distance_by` (e.g. 2) a pose that does not fit the box is
retried at a shorter distance, down to `min_distance` (default 10 Å); the
report lists every distance tried.

### Random poses

After the listed poses, `random_poses` (default 3) more are drawn, named
numbered on from the listed ones (`binderpose2`, `binderpose3`, ...). Each one
gets:

- a random approach direction, up to `random_max_tilt` (60°) from the normal
- a random distance in `random_distance` (default: `distance` to `distance + 5`)
- a random orientation

A draw is thrown away and drawn again (up to 2000 times) when:

- it fails the box or membrane checks
- its direction is within `random_min_angle` (25°) of another pose's
- its binder centre is within `random_min_separation` (10 Å) of another pose's

The same `random_seed` gives the same poses. Each report records the values
drawn, so a pose you like can be copied into `poses` as `distance`,
`orientation: euler`, `angles` and `approach`. Set `random_poses: 0` to place
only the listed poses.

For a binder whose input is already a meaningful pose (e.g. a G protein in its
receptor-bound arrangement), `random_from_input_position: true` moves the
binder out from its input position along each random direction, and
`random_rotated: N` gives only the last N random poses a random rotation; the
others keep the input orientation, so only their side and tilt are random.

### Choosing distances

The box repeats along z. Above the extracellular tip lies the intracellular side
of the next periodic copy. The binder therefore needs:

```
distance + binder height along z + min_image_gap ≤ free water column
```

The free water column is the box height minus the solute's z extent. Each
report gives `free_water_column_angstrom` and `column_room_for_binder_angstrom`.
The real 3D check is usually more forgiving, because the two tips are rarely
directly above each other. In the bundled D1R example the column is about
158 Å, so dopamine at 20 Å fits easily and so would a G protein.

## Output

```
<output>/addbinder_output/
  summary.txt, summary.json      one row per pose: status, gaps, removed molecules
  <pose>/step5_input.gro/.pdb    full system with the binder
  <pose>/topol.top, toppar/      merged force field + one ITP per molecule type
  <pose>/index.ndx               System Protein Ligand Host Binder Host_tip Solute Water Ion Membrane Environment
                                 + CHARMM-GUI SOLU MEMB SOLV SOLU_MEMB SYSTEM
  <pose>/*.mdp, README           the host's CHARMM-GUI run inputs (copy_run_inputs)
  <pose>/addbinder_report.txt/.json
  view.pml, view.pse, view.vmd   one PyMOL/VMD scene: host, tip atom, every pose's binder
```

A rejected pose writes only its report. A pose folder is replaced on reruns
only if it holds nothing but addbinder's own files. `Host_tip` and `Binder` are
ready-made groups for a distance restraint or pull coordinate during
equilibration. Without one, a free binder diffuses away.

## All options

| Option | Default | Meaning |
|---|---|---|
| `side` | `upper` | `upper` (+z) or `lower` (−z) |
| `distance` | 10.0 | Å, default for poses without their own |
| `min_image_gap` | 10.0 | Å, to periodic images of host and binder |
| `min_membrane_gap` | 5.0 | Å, to lipid heavy atoms |
| `lateral_reference` | `tip` | `tip` or `host_center` |
| `water_clash_distance` | 2.4 | Å, binder atom ↔ water/ion heavy atom |
| `same_orientation` | false | true: every `auto` pose gets the same orientation |
| `copy_run_inputs` | true | copy the host's `*.mdp` and `README` into each pose folder |
| `random_poses` | 3 | extra poses from random directions (0 for none) |
| `random_distance` | [distance, distance + 5] | Å, range for the random poses |
| `random_max_tilt` / `random_min_angle` | 60 / 25 | degrees |
| `random_min_separation` | 10 | Å between binder centres of any two poses |
| `random_seed` | 1 | same seed, same random poses |
| `random_rotated` / `random_from_input_position` | all / false | how many random poses are also rotated; start from the input position |
| `salt_concentration` | `host` | `host` (keep the host's) or mol/L |
| `neutralize` / `target_net_charge` | true / 0 | add counterions to reach the target charge |
| `ion_side` | `any` | `any` (whole box) or `binder_side` (binder's side first) |
| `ion_exclusion_distance` / `lipid_ion_exclusion_distance` | 10 / 5 | Å, protected ions |
| `environment_overrides` | {} | molecule type → water/ion/lipid/sterol/detergent/solvent |
| `output_dir` | `addbinder_output` | inside `--output` |
| `poses[]` | one pose straight above the tip | `name`, `distance`, `lateral_offset`, `orientation`, `angles`, `flip`, `spin`, `approach`, `centroid`, `from_input_position`, `distance_to`, `reduce_distance_by`, `min_distance` |

## Desktop application

In the desktop application, choose **Add binder outside the membrane protein** on
the Configuration tab, or run one of its two cards on the Examples tab: dopamine
(extracellular) or the Gs trimer (intracellular). Results go to `addbinder/`
(examples: `examples/addbinder/` and `examples/addbinder_gprotein/`) in the
output folder.

## Limitations

These are not supported yet:

- `--prepare-minimization`
- automatic box extension
- protein binders build (the Gs example) but have not yet been simulated
