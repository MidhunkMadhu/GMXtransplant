# salt example: 150 mM KCl instead of NaCl

The 300 Na⁺ and 314 Cl⁻ of the apo D1R build (`../addbinder/host`) are removed and
150 mM KCl is added in the bulk water, at least 5 Å from the receptor and the
lipids and 5 Å apart. POT is not in the host's toppar: its ITP and parameters
(atom type, NBFIX pairs) come from the bundled CHARMM36 ions.

```bash
cd examples/salt
python -B ../../src/run_pipeline.py --mode salt -i salt.yaml --output results/
```

Expect 291 K⁺ and 305 Cl⁻ (14 extra Cl⁻ neutralize the +14 e receptor),
0.150 M on 107,765 waters, net charge 0. Guide: [EDIT_MODES.md](../../docs/guides/EDIT_MODES.md).
