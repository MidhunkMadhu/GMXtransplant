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

| Component | Color | Initially visible |
|---|---|---|
| Incoming protein | Blue | Yes |
| Retained protein | Slate | Yes |
| Incoming ligand | Magenta | Yes |
| Retained ligand (topology-classified transplant mode) | Purple | Yes |
| Incoming cholesterol | Gold | Yes |
| Retained cholesterol | Teal | Yes |
| Retained membrane | Pale gray | Yes |
| Original replaced protein | Orange | No |
| Original replaced ligand | Orange | Yes |
| Removed environment | Red | No |
| Reference input | Amber | No |
| Aligned transplant input | Cyan | No |
| Water / ions | Pale blue | No |

Toggle named objects in PyMOL or molecules in VMD. The old ligand is a comparison overlay, **not part of the final molecular system**. Full input overlays and removed environment are initially hidden to avoid obscuring the result. `visualization/LEGEND.txt` and `scene.json` list the roles and atom counts actually present in that run.

Source overlays use the same rigid transformation as assembly, in the final protein's coordinate frame. Ligand mode shows the positioned incoming ligand; cholesterol mode shows the converted, aligned experimental structure. Box repair, when enabled, is already reflected in the reference overlay. These views do not apply a second fit, wrap molecules, or change scientific coordinates.

Inserted membership comes from the actual inserted atom block, rather than guessing by residue name: an old and a new ligand with the same residue name remain distinct. Two-folder mode additionally uses exact topology molecule categories. Other modes use their standard residue classification for display; unknown retained components are labeled “Other retained molecules.” Viewer files are for comparison, not simulation inputs; use the validated final GRO and topology for simulation.

Viewer export failures are reported without invalidating successfully written scientific outputs. Details remain in View Command Progress. The generated scripts work without network access on supported viewer installations on Linux, WSL, and macOS.
