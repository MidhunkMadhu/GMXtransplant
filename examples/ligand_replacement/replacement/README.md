The incoming ligand, named LI2 in this example (phenethylamine): LI2.mol2
(coordinates), its matching LI2.itp, and forcefield.itp (its bond, angle and
dihedral parameters). The coordinate atom order must match the ITP. Edit the old
and new residue names independently in ligand_replace.yaml. This example uses
pairfit on the atoms LI1 and LI2 share; the mapping is in the YAML.

forcefield.itp is the CHARMM-GUI force field generated together with LI2.itp.
CHARMM-GUI ligand ITPs list which atoms are bonded, but the bond, angle and
dihedral parameters are in forcefield.itp. It is selected with
new_ligand.forcefield_path; without it, a forcefield.itp next to the ITP is
used, and if the parameters are missing the run stops and names them.
