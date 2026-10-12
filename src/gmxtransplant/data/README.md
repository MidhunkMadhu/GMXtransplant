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
`POPC_512_NaCl_150mM_310K_md_500ns.gro`. The `popc` mode uses this frame by
default; give `popc.reference_gro` to use another equilibrated bilayer.
