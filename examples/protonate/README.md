# protonate example: protonation states in the apo D1R build

CHARMM-GUI built D2.50 (ASP70) and D3.49 (ASP120, DRY motif) protonated: their
atoms are those of ASPP although the residue name is ASP. This example makes
D2.50 charged again (ASPP → ASP, HD2 removed) and doubly protonates histidine 164
(HSD → HSP, HE2 added from the CHARMM internal coordinates). Residue numbers
are those of `PROA.itp`.

```bash
cd examples/protonate
python -B ../../src/run_pipeline.py --mode protonate -i protonate.yaml --output results/
```

Expect the receptor's charge unchanged (−1 + 1), the HSP atom types (NR3, HR1)
and their bonded parameters added from the bundled CHARMM36 library, net charge 0.
Guide: [EDIT_MODES.md](../../docs/guides/EDIT_MODES.md).
