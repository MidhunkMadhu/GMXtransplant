# addbinder mode — roadmap

Status (2026-10-04): ligand and protein binders build end to end (examples:
dopamine extracellular, Gs heterotrimer intracellular); neither has been
simulated yet. The user guide is
[ADDBINDER.md](ADDBINDER.md). It covers every option, how poses are placed, and
the output layout. This file lists only what is still open.

addbinder takes a prepared membrane system (CHARMM-GUI folder) and adds one
binder, a ligand or a protein (one or several molecules), in the water above or
below the membrane protein. It writes one run-ready folder per pose. Only water
and ions are removed. Code: `addbinder.py`; tests: `test_addbinder.py`; example:
`examples/addbinder/` (apo D1R host; dopamine above, 4 poses; Gs trimer below,
3 poses).

## Open tasks (in order)

### 1. `face_mask` orientation
Point chosen binder residues (e.g. Gα α5, a nanobody CDR) at the host, then
apply `spin`. This is needed for protein–protein setups.

### 2. Tip control
- `tip_mask` restricts the tip search to a region.
- `tip_method: percentile` stops a floppy terminus or a single side chain from
  setting the tip.
- `side_check` stops the run if a given residue is on the wrong side.

### 3. Restraint output
- Write a GROMACS pull-code or flat-bottom `.mdp` fragment between `Host_tip`
  and `Binder`, so a free binder does not drift away.
- Add matching restraints to the OpenMM bundle.
- `--prepare-minimization` is refused for addbinder until this is done.

### 4. Smaller items
- `run.yaml` per pose, so each pose reruns on its own.
- A distance line (tip → binder) in the PyMOL/VMD scene.
- A test with a protein binder whose ITP is also named `PROA`, to exercise the
  `_BND` rename and chain IDs.

### 5. Later (v2)
- Extend the box by replicating a bulk-water slab when it is too short (now
  the run stops and reports the height needed).
- Nudge the binder out on a clash.
- Swap water for ions when the binder's side lacks counterions.
- Several binder copies.

## Known limitations
- **Ion balance:** in a single-bilayer periodic box, the two water layers
  connect through the z boundary. `ion_side: binder_side` balances the ion
  counts at setup, but it cannot keep the compartments separate during MD.
- **Distance definition:**
  - Straight-above poses measure `distance` from the tip plane along the
    normal, so the true 3D gap is ≥ `distance`.
  - `approach` and random poses measure it from the nearest host heavy atom.
  - Both values are in each report.

## Done
- 2026-10-01: first version.
  - Membrane frame from lipid P atoms; tip = outermost heavy atom.
  - Orientations `flat`, `as_is`, `euler`, plus `spin` and `lateral_offset`.
  - Checks against the host, the membrane and all 26 periodic images.
  - Water/ion overlap removal and neutralization from the binder's side.
  - `_BND` rename on name clashes, ITP atom-name check, CGenFF coverage check.
  - Per-pose folders and a summary.
  - `gmx grompp` (0 warnings) and 500-step minimization checked on the old
    LI2 example.
- 2026-10-02: PyMOL/VMD scene and GUI mode with an example card.
- 2026-10-03: multi-molecule binders (`binder_itp` list, placed as one rigid
  body).
- 2026-10-04:
  - The example is now D1R + dopamine; the dopamine and receptor force fields
    merge without conflicts.
  - `auto` orientations (a different one per pose), plus `end_on`, `edge` and
    `flip`.
  - CHARMM-GUI index groups, and the host's mdp/README copied into each pose.
  - Random poses (direction, distance, orientation; redrawn when too similar).
  - Per-pose `approach` and `centroid`.
  - Poses are named `binderpose1`, `binderpose2`, …
  - Salt: the host's concentration is kept on the remaining water; extra ions
    removed, missing ones replace bulk waters; whole box by default.
  - Gs heterotrimer example below D1R (bound arrangement, `from_input_position`).
  - `reduce_distance_by` / `min_distance`, `distance_to: centroid` (spherical
    placement), `random_from_input_position`, `random_rotated`.
