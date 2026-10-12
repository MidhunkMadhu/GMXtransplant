# popc example: the apo D1 receptor in equilibrated POPC

The host is the apo D1R build from the addbinder example (`../addbinder/host`).
CHARMM-GUI Membrane Builder made it: PROA in a POPC/POPE/cholesterol membrane
with Na⁺/Cl⁻, placed by OPM. This run replaces the membrane with pure POPC
from the bundled 500 ns CHARMM36 bilayer. The receptor keeps its coordinates.

## Run it

From the source checkout root:

```bash
cd examples/popc
python -B ../../src/run_pipeline.py --mode popc -i popc.yaml --output results/
```

In the desktop application, use the **Equilibrated POPC membrane: D1R** card
on the Examples tab. For your own system, see the
[popc mode guide](../../docs/guides/POPC.md).

## What to expect

The run takes about 20 seconds and writes `results/popc_output/`:

| | |
|---|---|
| Host membrane removed | 261 POPC, 156 POPE, 105 CHL1 |
| Anchors | 3 POPC per leaflet, over 60 Å apart; both leaflets tilted 0.1° |
| P–P thickness, host / reference | 37.3 Å / 39.0 Å |
| POPC placed, upper / lower | 228 / 229 (55 removed for overlap with the receptor) |
| Hydration water | 10,501 reference waters; 185 host water/ion molecules kept at the receptor |
| Charge, salt | 0 e; 0.154 M, as in the host |

Then run the equilibration (about 1.9 ns):

```bash
cd results/popc_output/equilibration
./sample_run.sh
```

The system has about 380,000 atoms, so run the equilibration on a compute
node rather than a laptop. The guide shows the settings for a Slurm job.
