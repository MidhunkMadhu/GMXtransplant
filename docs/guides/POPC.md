# popc mode: replace the membrane with equilibrated POPC

`--mode popc` replaces the membrane of a prepared system with pure POPC taken
from a well-equilibrated bilayer. The protein and its ligands are not moved,
so the OPM placement that CHARMM-GUI built around is kept.

Use it to:
- start from lipids with equilibrated packing, tail order and headgroup
  hydration instead of a freshly built membrane
- turn a mixed membrane (e.g. POPC/POPE/cholesterol) into pure POPC

## Quick start

```bash
gmxtransplant --show-example popc > popc.yaml   # set host to your system folder
gmxtransplant --mode popc -i popc.yaml --dry-run
gmxtransplant --mode popc -i popc.yaml --output results/
```

The only required setting is `host`. To try the included D1 receptor system
first, see the [popc example](../../examples/popc/README.md). It is also a
card on the Examples tab of the desktop application.

## Inputs

| Key | What to give |
|---|---|
| `host` | Your system folder: `topol.top`, `toppar/` and `step5_input.gro` (or one `step5*.gro`; `host_gro` picks another GRO, e.g. an MD frame). A CHARMM-GUI Membrane Builder folder fits directly. The membrane must contain POPC: the host's own `POPC.itp` is used for the new lipids. |
| `reference_gro` | Optional. The equilibrated bilayer to take the POPC from. Leave it out to use the bundled bilayer (below). Your own bilayer works if its POPC atom names match the host `POPC.itp`; the atom order may differ. Three-site water in it is used as hydration water. |
| `membrane_center_z` | Optional. The membrane midplane in Å, in the host GRO frame. Leave it out to use the point midway between the host's phosphate planes, the centre the host was built around. |

### The bundled bilayer

The default reference is the last frame of a 500 ns simulation: 512 POPC and
20,480 waters, CHARMM36, 150 mM NaCl, 310 K. It comes from H. Martinez-Seara,
Zenodo [10.5281/zenodo.1010083](https://doi.org/10.5281/zenodo.1010083),
licensed CC BY 4.0. Only its ions were removed. It has a P–P thickness of
39.0 Å and 64.5 Å² per lipid. Cite the Zenodo record when you publish with it.

## What happens

1. **Measure the host membrane.** The midplane is found from the lipid tail
   ends, and each lipid is assigned to the upper or lower leaflet. In each
   leaflet:
   - at least 3 POPC must have P atoms 20 Å or more apart laterally, and not
     on one line. These are the anchors, which define the plane.
   - the plane fitted through the P atoms must be within 5° of the xy plane.
     The membrane must be flat and normal to z.
2. **Remove the host membrane.** Every lipid, sterol and detergent is
   removed. The protein, ligands, water and ions stay.
3. **Place the new bilayer.** The reference is repeated periodically and cut
   to the host box. Its midplane goes on the host midplane, and it keeps its
   own equilibrated thickness. Any box size works, including hexagonal boxes.
4. **Remove overlapping POPC, whole lipids only.**
   - Any POPC with a heavy atom within 3 Å of a protein or ligand heavy atom
     is removed. At 3 Å, no lipid tail can thread through an aromatic ring.
   - Where copies of the reference meet at the box edge, the fewest are
     removed so that none overlap.
   - The bilayer is tried at 16 lateral positions. The one that keeps the most
     POPC is used.
5. **Hydrate the new bilayer.** The reference water hydrates the new
   headgroups, up to 10 Å beyond each phosphate plane, and replaces the host
   water there. Host water and ions close to the protein are kept, including
   internal waters; the reference water yields to them. Outside that layer,
   host water and ions are removed only where they overlap the new POPC or
   the reference water.
6. **Fix charge and salt.** Counterions far from the protein and lipids are
   removed to reach the target charge, for example after anionic lipids were
   replaced. The reference water brings no ions. Bulk waters away from the
   protein therefore become salt pairs, using the host's own ion types, until
   the host's salt concentration is reached again.
7. **Write and check.** The topology is checked against the coordinates atom
   by atom, including the net charge. A failed run leaves only its report,
   with the reason.

## Output

In `<output>/popc_output/`:

| File | Contents |
|---|---|
| `step5_input.gro`, `step5_input.pdb` | The new system |
| `topol.top`, `toppar/` | Its topology |
| `index.ndx` | CHARMM-GUI groups (`SOLU`, `MEMB`, `SOLV`, `SOLU_MEMB`) plus `Protein`, `Membrane`, `POPC_upper`, `POPC_lower`, ... |
| `popc_report.txt`, `.json` | Host planes and anchors; lipids removed per leaflet; the reference used; POPC placed per leaflet; water and ions removed and why; charge and salt; final molecules |
| `equilibration/` | The GROMACS equilibration (next section) |
| `view.pml`, `view.vmd`, `view.pse` | PyMOL/VMD scenes, with the new POPC in its own colour |

## Equilibration

The removed lipids leave gaps around the protein, and the edges of the
hydration layer meet the host water. Run `equilibration/` before production.
It is a self-contained GROMACS folder (`topol.top`, `toppar/`, `step5_input.gro`,
`index.ndx`, one mdp per stage and `sample_run.sh`). It follows the
CHARMM-GUI Membrane Builder step6 schedule, about 1.9 ns in total:

| Stage | Type | Length | Restraints: backbone / side chains / lipid (kJ/mol/nm²) |
|---|---|---|---|
| step6.0_minimization | steepest descent | 5000 steps | 4000 / 2000 / 1000 |
| step6.1_equilibration | NVT | 125 ps (1 fs) | 4000 / 2000 / 1000 |
| step6.2_equilibration | NVT | 125 ps (1 fs) | 2000 / 1000 / 400 |
| step6.3_equilibration | NPT | 125 ps (1 fs) | 1000 / 500 / 400 |
| step6.4_equilibration | NPT | 500 ps (2 fs) | 500 / 200 / 200 |
| step6.5_equilibration | NPT | 500 ps (2 fs) | 200 / 50 / 40 |
| step6.6_equilibration | NPT | 500 ps (2 fs) | 50 / 0 / 0 |

The protein is released step by step, so it keeps its OPM placement while
the POPC packs around it. The restraints use the `POSRES_FC_BB` and
`POSRES_FC_SC` blocks in the protein and ligand ITPs, and the `POSRES_FC_LIPID`
(P in z) and `DIHRES_FC` blocks in the POPC ITP. NPT stages use semi-isotropic
C-rescale coupling at 1 bar and 310 K. The report warns if a protein or
ligand ITP leaves heavy atoms without a restraint.

Run it on your workstation:

```bash
cd results/popc_output/equilibration
./sample_run.sh
```

On a cluster, set the GROMACS commands. For example, in a Slurm job with an
MPI build:

```bash
GMX="srun -n 1 gmx_mpi" MDRUN="srun gmx_mpi mdrun" ./sample_run.sh
```

Then check the stages with `python3 analyze_equilibration.py` (box, area per
phospholipid, P–P thickness, protein C-alpha RMSD and z drift per stage), and
continue production from `step6.6_equilibration.gro` with the same `topol.top`.
`gmxtransplant slurm-script results/popc_output/equilibration --profile dardel
--account <project>` writes a Slurm job for the whole schedule.

On Dardel (GROMACS 2026.1, one node, 53 ns/day) the D1R example ran the full
schedule with no LINCS warnings: the box shrank laterally from 12.95 to 12.56 nm
as the gaps closed, the P–P thickness settled at 3.83 nm, and the receptor
C-alpha RMSD was 0.09 nm after the restraints were released, its centre
moving less than 0.1 nm along z. To change a stage, e.g. to lengthen the last one, edit
`equilibration.stages` in the YAML.

## All settings

| Key | Default | Meaning |
|---|---|---|
| `host` | (required) | System folder |
| `host_gro` | step5 GRO | Other coordinates from the host folder |
| `reference_gro` | bundled bilayer | Equilibrated POPC bilayer |
| `membrane_center_z` | between the P planes | Midplane z (Å, host GRO frame) |
| `anchor_count` | 3 | POPC anchors needed per leaflet |
| `anchor_separation` | 20.0 | Smallest lateral distance between anchors (Å) |
| `max_tilt` | 5.0 | Largest leaflet tilt from the xy plane (degrees) |
| `protein_clash_distance` | 3.0 | POPC heavy atom to protein/ligand heavy atom (Å) |
| `lateral_search` | 16 | Lateral positions of the bilayer tried |
| `water_layer` | 10.0 | Reference hydration beyond each P plane (Å); 0 keeps only host water |
| `protein_hydration_distance` | 4.0 | Host water/ions this close to the protein are kept (Å) |
| `reference_contact_distance` | 2.0 | Allowed contact between reference copies at the box edge (Å) |
| `water_clash_distance` | 2.5 | Host water/ion to new POPC or reference water (Å) |
| `neutralize`, `target_net_charge` | true, 0 | Remove counterions to reach this charge |
| `ion_exclusion_distance` | 10.0 | Ions this close to the protein are never removed (Å) |
| `lipid_ion_exclusion_distance` | 5.0 | Ions this close to lipids are never removed (Å) |
| `restore_salt` | true | Restore the host salt concentration |
| `environment_overrides` | {} | Classify unknown molecule types, e.g. `MYLIP: lipid` |
| `equilibration.enabled` | true | Write `equilibration/` |
| `equilibration.temperature` | 310.0 | K |
| `equilibration.stages` | table above | Each stage: `name`, `integrator` (steep/nvt/npt), `nsteps`, `dt`, `backbone`, `sidechain`, `lipid` |

## Limits

- The host membrane must contain POPC. The result is pure POPC: other lipids
  and sterols are not kept.
- The membrane normal must lie along z, the box's third vector.
- The new bilayer keeps the reference's thickness. If it differs from the
  host's by more than 3 Å, the report warns; NPT equilibration adjusts the
  water layers.
