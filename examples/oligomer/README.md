# oligomer example: two D1 receptors in one membrane

A copy of the apo D1R is placed half the box away in x and y (64.76 Å, 64.76 Å)
and turned 180° about the membrane normal: a checkerboard lattice in which each
receptor's nearest neighbours are copies of the other. The receptor is 66 Å wide
in x and 54 Å in y, so a shift along one axis alone would bring the copy within
10 Å of the original or its periodic image.

```bash
cd examples/oligomer
python -B ../../src/run_pipeline.py --mode oligomer -i oligomer.yaml --output results/
```

Expect 2 PROA, closest approach 27.7 Å (periodic images included), 57 lipids
removed whole (25 POPC, 23 POPE, 9 CHL1), 12 Na⁺ removed to neutralize.
Guide: [EDIT_MODES.md](../../docs/guides/EDIT_MODES.md).
