# resize example: a thicker water layer

20 Å of water is added along z to the apo D1R build (`../addbinder/host`): the
box grows from 245.9 to 265.9 Å. Copies of the host's own bulk water slab
(water and ions) fill the new space in the middle of the water layer, so the
salt stays at the host's concentration.

```bash
cd examples/resize
python -B ../../src/run_pipeline.py --mode resize -i resize.yaml --output results/
```

Expect 9,886 waters and 30 Na⁺/30 Cl⁻ added (bulk density would give about
11,200: overlapping copies are dropped at the two seams, and the NPT
equilibration closes the gap), salt 0.155 M, net charge 0.
Guide: [EDIT_MODES.md](../../docs/guides/EDIT_MODES.md).
