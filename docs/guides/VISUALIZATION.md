# Molecular comparison views

Each successful assembly exports `view.pml`, `view.vmd`, and `visualization/` alongside its coordinates. A real `view.pse` session is generated automatically when PyMOL is installed; PyMOL is found on the PATH, in the running Python environment, or on macOS in `/Applications`, and `GMXTRANSPLANT_PYMOL` can point to any other copy. PyMOL and VMD are optional external viewers; neither becomes a CLI or GUI dependency. Minimization preparation shows the prepared input, not an unperformed minimization result.

Use **Open in PyMOL** or **Open in VMD** in GUI results, or run from the output folder:

```bash
pymol view.pml
vmd -e view.vmd
```

Every mode that produces a visualization writes `view.vmd` and `view.pse` into
its output folder, alongside `view.pml` and the `visualization/` assets, so the
assembled system can be reopened later without rerunning the pipeline. `view.pse`
is written during the run when PyMOL is on `PATH`; when it is not, the run says
so and opening `view.pml` saves `view.pse`. The PSE is self-contained. When moving a scripted scene, copy the entire output folder, including `visualization/`, and launch from that folder. VMD locates its assets relative to the `.vmd` file. For Python-script loading in PyMOL, use `run visualization/load.py, module`.

Every view also contains **what the run changed**, as objects that start
hidden: switch one on to see what was removed or replaced. Nothing hidden is
part of the final system.

| Component | Color | Initially visible |
|---|---|---|
| Incoming protein | Blue | Yes |
| Retained protein | Slate | Yes |
| Incoming ligand | Magenta | Yes |
| Retained ligand (topology-classified transplant mode) | Purple | Yes |
| Incoming cholesterol | Gold | Yes |
| Retained cholesterol | Teal | Yes |
| Retained membrane | Pale gray | Yes |
| New molecules of the mode (new bilayer, flood copies, protein copies, edited protein) | Mode colour | Yes |
| Water / ions | Pale blue | No |
| `removed_protein`: original protein replaced, or the original atoms of edited residues (mutate, protonate) | Orange | No |
| `removed_ligand`: original ligand replaced | Orange | No |
| `removed_membrane`: original lipids and sterols removed or replaced (e.g. the host membrane in the membrane mode, lipids under oligomer copies) | Red-orange | No |
| `removed_water`, `removed_ion`: water and ions removed (salt, flood, membrane, addbinder's first pose, ...) | Light red / red | No |
| `removed_environment`: other removed molecules | Red | No |
| `host_input`: the whole input, when the mode moves molecules (resize) | Gray | No |
| Reference input | Amber | No |
| Aligned transplant input | Cyan | No |

How to show them: in PyMOL, click an object's name in the panel on the right
(or `enable removed_membrane`); in VMD, double-click the **D** column of the
molecule in the Main window (or `mol on <id>`). In PyMOL they are drawn semi-transparent,
so the final system stays readable. `visualization/LEGEND.txt` and `scene.json`
list the objects and atom counts actually present in that run.

Source overlays use the same rigid transformation as assembly, in the final protein's coordinate frame. Ligand mode shows the positioned incoming ligand; cholesterol mode shows the converted, aligned experimental structure. Box repair, when enabled, is already reflected in the reference overlay. These views do not apply a second fit, wrap molecules, or change scientific coordinates.

Inserted membership comes from the actual inserted atom block, rather than guessing by residue name: an old and a new ligand with the same residue name remain distinct. Two-folder mode additionally uses exact topology molecule categories. Other modes use their standard residue classification for display; unknown retained components are labeled “Other retained molecules.” Viewer files are for comparison, not simulation inputs; use the validated final GRO and topology for simulation.

Viewer export failures are reported without invalidating successfully written scientific outputs. Details remain in View Command Progress. The generated scripts work without network access on supported viewer installations on Linux, WSL, and macOS.
