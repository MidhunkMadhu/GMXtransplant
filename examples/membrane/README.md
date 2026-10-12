# membrane example: an equilibrated POPC/cholesterol membrane

The POPC/POPE/cholesterol membrane of the apo D1R build (`../addbinder/host`) is
replaced by the lipids of an equilibrated POPC/cholesterol 70:30 bilayer
(`reference/`, 100 ns of CHARMM36 at 303 K, Zenodo 10.5281/zenodo.209362, CC BY 4.0;
see `reference/README.md`). The receptor keeps its OPM placement.

```bash
cd examples/membrane
python -B ../../src/run_pipeline.py --mode membrane -i membrane.yaml --output results/
```

Expect 399 POPC and 199 CHL1 placed (best of 16 lateral registrations), 11,166
reference waters, salt 0.154 M as in the host, net charge 0. The reference is
7.9 Å thicker than the host membrane (cholesterol orders the tails); the NPT
equilibration adjusts the water layers. Guide: [EDIT_MODES.md](../../docs/guides/EDIT_MODES.md).
