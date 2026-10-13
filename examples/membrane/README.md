# membrane example: the apo D1 receptor in equilibrated POPC

The host is the apo D1R build from the addbinder example (`../addbinder/host`).
CHARMM-GUI Membrane Builder made it: PROA in a POPC/POPE/cholesterol membrane
with Na⁺/Cl⁻, placed by OPM. This run replaces the membrane with pure POPC
from the bundled 500 ns CHARMM36 bilayer. The receptor keeps its coordinates.

## Run it

From the source checkout root:

```bash
cd examples/membrane
python -B ../../src/run_pipeline.py --mode membrane -i membrane.yaml --output results/
```

In the desktop application, use the **POPC membrane: D1R** card
on the Examples tab. For your own system, see the
[membrane mode guide](../../docs/guides/MEMBRANE.md).

## What to expect

The run takes about 20 seconds and writes `results/membrane_output/`:

| | |
|---|---|
| Host membrane removed | 261 POPC, 156 POPE, 105 CHL1 |
| Host midplane | z = 123.0 Å, from the lipid tail ends; both leaflets tilted 0.1° |
| P–P thickness, host / reference | 37.3 Å / 39.0 Å |
| POPC placed, upper / lower | 232 / 231 (49 removed within 2 Å of the receptor) |
| Hydration water | 10,684 reference waters; host water and ions near the receptor kept |
| Charge, salt | 0 e; 0.154 M, as in the host |

Then run the equilibration (about 1.9 ns):

```bash
cd results/membrane_output/equilibration
./sample_run.sh
```

The system has about 380,000 atoms, so run the equilibration on a compute
node rather than a laptop. The guide shows the settings for a Slurm job.

## Your own bilayer: POPC/cholesterol 70:30

`popc_cholesterol.yaml` gives a bilayer GRO in place of the bundled POPC:
`reference/popc_chol_70-30_charmm36_303K_100ns.gro.gz`, 100 ns of CHARMM36 at
303 K (J. Melcr, Zenodo 10.5281/zenodo.209362, CC BY 4.0; see
`reference/README.md`). Its cholesterol is named CLOL; the run matches it to the
host's CHL1 by its atom names, so the YAML holds only `host` and `reference_gro`.

```bash
python -B ../../src/run_pipeline.py --mode membrane -i popc_cholesterol.yaml --output results/
```

Expect about 400 POPC and 200 CHL1 placed, reference hydration water, the
host's salt concentration and a net charge of 0. The reference is about 8 Å
thicker than the host membrane (cholesterol orders the tails); the NPT
equilibration adjusts the water layers. In the desktop application this is the
**POPC/cholesterol membrane: D1R** card under **More examples**.
