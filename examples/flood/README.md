# flood example: dopamine flooding of the apo D1 receptor

Ten dopamine molecules (CGenFF residue LDP, +1 e each, from `../addbinder/dopamine`)
are placed at random positions and orientations in the water around the apo D1R
build, at least 10 Å from the receptor, 5 Å from the membrane and 8 Å from each
other: the start of an unbiased binding (ligand flooding) simulation.

```bash
cd examples/flood
python -B ../../src/run_pipeline.py --mode flood -i flood.yaml --output results/
```

Expect 10 LDP, 143 waters and 1 Cl⁻ removed where they overlap, 11 Na⁺
removed to neutralize, index group `Flood_LDP`. Guide: [EDIT_MODES.md](../../docs/guides/EDIT_MODES.md).
