# Two-folder protein transplant example

From this directory, run either the installed package:

```bash
gmxtransplant --mode charmprot -i charmprot.yaml --dry-run
gmxtransplant --mode charmprot -i charmprot.yaml
```

or the source checkout:

```bash
python -B ../../run_pipeline.py --mode charmprot -i charmprot.yaml --dry-run
python -B ../../run_pipeline.py --mode charmprot -i charmprot.yaml
```

The reference and transplant datasets were supplied with the new modality. Input folders are read only. Results replace generated files in `charmprot_output/`. In the GUI, use the first card in Examples; results go to `<chosen output>/examples/charmprot/`.
