# Reference bilayer: POPC/cholesterol 70:30

`popc_chol_70-30_charmm36_303K_100ns.gro.gz` is `confout.gro`, the frame at the
end of the 100 ns simulation of

> J. Melcr, *POPC/Cholesterol (70:30) lipid membrane, 303K, Charmm36 force field,
> simulation files and 100 ns trajectory for Gromacs simulation engine v5.1.2*,
> Zenodo (2016), doi:[10.5281/zenodo.209362](https://doi.org/10.5281/zenodo.209362).
> Licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

System: 84 POPC and 36 cholesterol (42 POPC + 18 cholesterol per leaflet),
4,800 TIP3P waters, CHARMM36, 303 K, made for the NMRlipids project. Unchanged
apart from gzip compression; cholesterol is named CLOL (atom names as CHARMM's
CHL1, which `resname_map: {CLOL: CHL1}` in `membrane.yaml` matches).
