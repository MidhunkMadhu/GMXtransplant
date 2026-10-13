# Bundled reference bilayer

`popc512_charmm36_nacl150mM_310K_500ns.gro.gz` is the last frame (500 ns) of
the simulation `POPC_512_NaCl_150mM_310K` from

> H. Martinez-Seara, *Simulations POPC bilayers (512 lipids) using charmm36 ff in gromacs*,
> Zenodo (2017), doi:[10.5281/zenodo.1010083](https://doi.org/10.5281/zenodo.1010083).
> Licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

System: 512 POPC, 20,480 water, 55 Na⁺ and 55 Cl⁻ (150 mM NaCl), CHARMM36,
310 K, GROMACS 5.0.4, semi-isotropic pressure coupling, 500 ns.

Change made for GMXtransplant: the Na⁺ and Cl⁻ ions were removed. The POPC
and water coordinates and the box are unchanged from
`POPC_512_NaCl_150mM_310K_md_500ns.gro`. The `membrane` mode uses this frame by
default; give `membrane.reference_gro` to use another equilibrated bilayer.

# CHARMM36 protein and ion library

`charmm36_protein_ions.json.gz` holds the standard amino acids and their
protonation variants (atoms, types, charges, bonds, impropers, CMAP and the
CHARMM internal-coordinate tables), the CHARMM36 ions, and every parameter row
whose atom types belong to proteins, ions or TIP3P water. The protonate, mutate
and salt modes add from it only the parameters a host's forcefield.itp lacks.

It is built by `tools/build_charmm36_library.py` from the MacKerell lab's
CHARMM36 force field (https://mackerell.umaryland.edu/charmm_ff.shtml):

- `charmm36-jul2022.ff` (the GROMACS port made with charmm2gmx by A. Wacha,
  after the port by E. P. Raman, J. A. Lemkul, R. Best and A. D. MacKerell Jr.):
  `aminoacids.rtp`, `ions.itp`, `ffnonbonded.itp`, `ffbonded.itp`, `cmap.itp`, `nbfix.itp`;
- `top_all36_prot.rtf` from `toppar_c36_jul24.tgz`, for the internal coordinates.

References: J. Huang et al., *Nat. Methods* 14, 71 (2017) (CHARMM36m);
R. B. Best et al., *J. Chem. Theory Comput.* 8, 3257 (2012).
