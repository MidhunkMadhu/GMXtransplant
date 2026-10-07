# Two-folder CHARMM-GUI protein transplantation

Run the installed package:

```bash
gmxtransplant --mode charmprot -i charmprot.yaml --dry-run
gmxtransplant --mode charmprot -i charmprot.yaml
```

Or from the source checkout:

```bash
python -B run_pipeline.py --mode charmprot -i charmprot.yaml --dry-run
python -B run_pipeline.py --mode charmprot -i charmprot.yaml
```

The entire required YAML is:

```yaml
charmprot:
  reference: /path/to/reference/gromacs
  transplant: /path/to/donor/gromacs
```

Each folder must contain `topol.top`, `toppar/`, and `step5_input.gro` (or exactly
one `step5*.gro`). Both folder paths, and any `reference_gro`/`transplant_gro`
override, must be written out in full: nothing is resolved against the YAML's
directory or the working directory. Inputs are read-only. `output_dir` defaults
to the name `charmprot_output`, placed in the `--output` directory when you pass
one and otherwise in the directory you ran the command from. Rewriting an existing `output_dir` is
the norm: generated artifacts are replaced transactionally, with rollback if publication fails. Unrelated files and GUI logs/locks are preserved. A failed assembly publishes failure reports and removes stale generated results.
`--dry-run` validates the YAML without opening the input
folders.

## Detection and selection

The program reads the actual topology include graph and expands `[ molecules ]`
using ITP atom counts. Every coordinate atom must agree with its ITP residue and
atom name in file order. GRO serial-number wrapping does not affect selection.
Unreferenced ITP files do not affect classification.

`PROA`, `PROB`, and other `PRO...` topology molecule types identify proteins.
An additional backbone check recognizes protein types with other names. Every
detected protein is selected by default. A covalent cofactor inside the same
protein moleculetype travels with that molecule; `ligands: ignore` only applies
to separate ligand molecules of the transplant.

Environment categories are water, ion, lipid, sterol, detergent and solvent.
The bundled catalogue is a dated snapshot of **943 distinct identifiers from
19 primary-source pages**, supplemented with exported water/ion/lipid aliases.
It includes every download identifier on the retrieved Individual Lipid Library
page and every `RESI` identifier in the retrieved lipid, water/ion and noble-gas
topology files linked from the CHARMM-GUI topology index. It does not treat
patch identifiers (`PRES`) as complete molecules.

Primary sources:

- [CHARMM-GUI Individual Lipid Library](https://www.charmm-gui.org/?doc=archive&lib=lipid)
- [CHARMM-GUI topology index](https://www.charmm-gui.org/?doc=toppar)
- [Water and ions topology](https://www.charmm-gui.org/?doc=open_toppar&filename=toppar_water_ions.str)
- [Extended ion topology](https://www.charmm-gui.org/?doc=open_toppar&filename=stream/misc/toppar_ions_won.str)
- [CHARMM-GUI water-model FAQ](https://www.charmm-gui.org/?doc=faq)
- [OpenMM CHARMM water distribution](https://github.com/openmm/openmmforcefields/blob/main/charmm/files/waters.yaml)

The catalogue records a retrieval date, source hashes and per-name source
references in `gmxtransplant/charmm_environment.json`. `charmm_environment.py`
contains export aliases. Ordinary runs work entirely offline.

Names are matched exactly, case-insensitively. Both molecule names and the ITP
residue names are examined. Thus all catalogued `POP...` lipids are covered
without an unrestricted `POP*` rule. Single charged atoms also qualify as ions
even when their name is custom. Unknown non-protein molecules become **inferred
ligands** and are named explicitly in the report. Topology position after the
protein and before the environment is recorded as supporting evidence, not
proof of chemical identity.

No fixed name list can guarantee every future, custom or renamed CHARMM-GUI
component. The snapshot's coverage is stated precisely above; unfamiliar names
remain visible in the report. Use an explicit ligand list or environment
override when a name has a different biological role (e.g. a bound metal or
functional cholesterol molecule).

## Optional YAML

```yaml
charmprot:
  reference: /path/to/reference/gromacs
  transplant: /path/to/donor/gromacs
  reference_proteins: [PROA, PROB]
  transplant_proteins: [PROB, PROA]
  reference_ligands: [LI1]
  transplant_ligands: [LI2]
  output_dir: transplant_result
```

Protein lists specify which complete molecule types to replace and the chain
pairing order. Multiple instances of a selected type are all selected, in file
order. Omitted protein types in the reference remain in place. Ligand lists
take precedence over the environment catalogue; `[]` explicitly selects no
ligands on that side. Omission means automatic detection on that side.

`ligands` chooses what is inserted from the transplant:

- `auto` (default): the reference's proteins and ligands are removed, and the
  transplant's proteins and ligands are inserted.
- `ignore`: the reference's proteins and ligands are removed in the same way,
  but only the transplant's proteins are inserted; its ligands are left out.
  `reference_ligands` may still name which reference ligands to remove;
  `transplant_ligands` cannot be combined with `ignore`.

If a retained protein or ligand (one not selected for removal) clashes with the
incoming block, the run stops and reports the conflict.

Additional options:

| Option | Default | Meaning |
|---|---|---|
| `reference_gro`, `transplant_gro` | automatic | Full path to a coordinate file, when the folder's own `step5*.gro` is not the one to use |
| `environment_overrides` | `{}` | Explicit molecule-name to category mapping, e.g. `{MYLIP: lipid}` |
| `clash_distance` | `0.80` | Heavy-atom environment clash distance in Å; removes complete molecules |
| `keep_lipids` | `false` | Keep every lipid and sterol, even when it clashes; only other clashing molecules (water, ions, solvent) are removed. Retained clashes are listed in the report; minimize before simulating |
| `neutralize` | `true` | Remove eligible counterions to restore target charge |
| `target_net_charge` | `0.0` | Desired charge in elementary-charge units |
| `ion_exclusion_distance` | `10.0` | Minimum ion distance from inserted/retained biological molecules in Å |
| `lipid_ion_exclusion_distance` | `5.0` | Minimum ion distance from retained membrane heavy atoms in Å |

## Alignment, topology and outputs

Corresponding protein chains are globally aligned by residue sequence, followed
by a single Cα Kabsch rigid fit over all matched chains. This supports missing
residues and substitutions without PyMOL. Chain counts must agree; low-identity
pairings stop with a clear error. Symmetric chain permutations are not searched
automatically: give corresponding protein lists if file order differs.
All selected donor atoms, including ligands, receive exactly the same transform.
The donor protein should be whole under periodic boundaries before alignment.

Environment clashes use minimum-image distances with the reference box. Whole
topology molecules are removed, including multi-residue lipids. A 0.02 Å margin
allows for GRO rounding. Exact molecule charges come from ITP `[ atoms ]`
records. Charge balancing removes distant counterions with the needed charge;
it stops if the requested charge cannot be reached. It does not add ions.

The parameter sets are combined, identical definitions are deduplicated, and
conflicting definitions stop the run. Protein and ligand molecule definitions
come from the donor; retained molecules keep reference definitions. Conditional
include graphs, parameter macros and nested molecule includes need preprocessing;
unsupported ambiguity is rejected. Standard flat CHARMM-GUI toppar inputs are
the intended format.

On success the output folder contains:

- `step5_input.gro`, `step5_input.pdb`: combined coordinates in the reference box
- `topol.top`, `toppar/`: assembled parameters and molecule definitions
- `index.ndx`: topology-based protein, ligand, environment and component groups
- `aligned_inspection.pdb`: reference and transformed donor selection overlay
- `charmprot_report.txt`, `charmprot_report.json`: both systems' classifications,
  selection decisions, unknown names, alignment, molecule removals and charge audit

The final topology is reparsed independently and checked against the written GRO
atom order and charge. Final files are published together from temporary storage.
On failure only diagnostic reports are published, with `status: failed`.
The reports are the record of exactly what was considered environment and ligand.

`--prepare-minimization` also exports a portable minimization bundle after a
successful transplant; it does not run minimization.

## GUI, examples and visualization

In the desktop application, select **CHARMM-GUI protein transplant** and choose the Reference and Transplant folders. Select the output folder at the top of the window. Other settings are under Advanced settings. The Examples card runs the supplied datasets in `examples/charmprot/`.

Canonical configurations use a `charmprot:` section and may use the shared `paths:` registry. The original flat two-folder YAML is still accepted. A root `minimization:` section uses the same settings as other modes; `--prepare-minimization` is also supported.

All successful modes produce `view.pml`, `view.vmd`, and a `visualization/` asset directory. A `view.pse` session is also created when PyMOL is installed. See [Visualization](VISUALIZATION.md) for colors and input overlays.
