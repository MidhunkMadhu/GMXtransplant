# mutate example: D3.32A and S5.42A in the apo D1R build

Two binding-site residues of D1R become alanine: D3.32 (ASP103), whose
carboxylate binds the dopamine amine, and S5.42 (SER198). The backbone and CB
keep their coordinates; the alanine hydrogens are built from the CHARMM
internal coordinates; the protein ITP gets the ALA atoms and bonded terms.
Residue numbers are those of `PROA.itp`.

```bash
cd examples/mutate
python -B ../../src/run_pipeline.py --mode mutate -i mutate.yaml --output results/
```

Expect the receptor's charge +1 e higher (one aspartate gone), one Na⁺ removed to
neutralize, net charge 0. For a larger side chain (e.g. `to: PHE`) the report
gives the chosen rotamer and its clash score. Guide: [EDIT_MODES.md](../../docs/guides/EDIT_MODES.md).
