# hmr example: hydrogen mass repartitioning of the apo D1R build

The host is the apo D1R build from the addbinder example (`../addbinder/host`):
PROA in a POPC/POPE/cholesterol membrane with 150 mM NaCl. Every hydrogen bonded
to a heavy atom gets 3.024 u, taken from that heavy atom; water (rigid, SETTLE)
is unchanged. The host's mdp files are copied; the production stage (step7) gets
dt = 0.004 ps with the same simulated length.

```bash
cd examples/hmr
python -B ../../src/run_pipeline.py --mode hmr -i hmr.yaml --output results/
```

Expect `results/hmr_output/`: 40,718 hydrogens repartitioned (PROA, CHL1, POPC,
POPE), lightest heavy atom 5.96 u, coordinates unchanged, `equilibration/`.
Checked on Dardel (GROMACS 2026.1): grompp 0 warnings; minimization, 5 ps NVT and
10 ps at 4 fs ran with no LINCS warnings. Guide: [EDIT_MODES.md](../../docs/guides/EDIT_MODES.md).
