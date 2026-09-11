"""
Residue -> molecule-class classification (lipid / water / ion / ligand /
other), used for per-class clash thresholds, reporting, and eligible-ion
selection during charge neutralization.
"""

from __future__ import annotations

from typing import Dict, Set

# Reasonably comprehensive defaults for CHARMM-GUI / CHARMM & AMBER force
# field naming conventions. Extend via the config's classification
# overrides for anything project-specific (custom lipids, non-standard
# ion names, etc).
DEFAULT_WATER_RESNAMES = {
    "TIP3", "TIP", "TIP4", "TIP5", "SPC", "SPCE", "WAT", "HOH",
    "T3P", "T4P", "SOL",
}

DEFAULT_ION_RESNAMES = {
    "SOD", "NA", "NA+",
    "CLA", "CL", "CL-",
    "POT", "K",
    "CAL", "CA2", "CA",
    "MG", "ZN", "ZN2",
    "CES", "CS",
    "LIT", "LI",
    "RUB", "RB",
}

# Common membrane residue and moleculetype names used by CHARMM-GUI,
# CHARMM36, AMBER lipid force fields, and related GROMACS systems. This set is
# also used to compare a reference topol.top with the available environment
# ITP files before restoring PDB-truncated names.
COMMON_LIPID_RESNAMES = {
    # Phosphatidylcholine
    "POPC", "DOPC", "DPPC", "DMPC", "DLPC", "SOPC", "DSPC",
    "PSPC", "OSPC", "PLPC", "OLPC", "SAPC",
    # Phosphatidylethanolamine
    "POPE", "DOPE", "DPPE", "DMPE", "DLPE", "SOPE", "DSPE",
    "PSPE", "OSPE", "PLPE", "OLPE", "SAPE",
    # Phosphatidylserine
    "POPS", "DOPS", "DPPS", "DMPS", "DLPS", "SOPS", "DSPS",
    "PLPS", "SAPS",
    # Phosphatidylglycerol
    "POPG", "DOPG", "DPPG", "DMPG", "DLPG", "SOPG", "DSPG",
    "PLPG", "SAPG",
    # Phosphatidic acid
    "POPA", "DOPA", "DPPA", "DMPA", "DLPA", "SOPA", "DSPA",
    "PLPA", "SAPA",
    # Phosphatidylinositol and phosphoinositides
    "POPI", "DOPI", "DPPI", "DMPI", "DLPI", "SOPI", "DSPI",
    "PLPI", "SAPI", "PIP2", "PIP3",
    # Sphingomyelin and ceramide families
    "PSM", "SSM", "DPSM", "BSM", "CER", "CER1", "CER2",
    # Cardiolipins
    "TOCL", "TLCL2", "CDL2",
    # Cholesterol names
    "CHL1", "CHOL", "CHL",
}

DEFAULT_LIPID_RESNAMES = COMMON_LIPID_RESNAMES

AMINO_ACID_RESNAMES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "HSD", "HSE", "HSP", "HID", "HIE", "HIP", "CYX", "ASH", "GLH", "LYN",
}


def is_hydrogen_name(atom_name: str) -> bool:
    """True if this looks like a hydrogen atom in standard protein/lipid/
    water/ion force-field naming.

    An earlier version of this function excluded names starting with
    "HE"/"HG"/"HF"/"HO"/"HS" as a defensive guard against the real
    elements He/Hg/Hf/Ho/Hs. That guard was wrong for this domain: "HE"
    and "HG" in particular are extremely common HYDROGEN name prefixes in
    standard amino acid nomenclature (His HE1, Trp HE1/HE3, Arg HE, Tyr
    HE1/HE2, Gln HE21/HE22, Glu HE2 on the protonated form; Ile HG1/HG2/
    HG11/HG12/HG13, Val/Leu HG-series, Ser/Cys/Thr hydroxyl/thiol HG),
    so that guard was silently letting real sidechain hydrogens slip
    through heavy-atom-only clash filtering. None of Hg/Hf/Ho/Hs/He ever
    appear as real atoms in a standard protein/lipid/water/ion force
    field, so there is no ambiguity worth guarding against here: any
    atom name starting with 'H' is treated as hydrogen.
    """
    name = atom_name.strip().upper()
    if not name:
        return False
    return name[0] == "H"


def classify_resnames(resnames, overrides: Dict[str, str] = None) -> Dict[str, str]:
    """Return {resname: class} for every unique resname in `resnames`.
    class in {'protein', 'water', 'ion', 'lipid', 'ligand', 'other'}.
    `overrides` (resname -> class) takes precedence over the defaults."""
    overrides = overrides or {}
    result = {}
    for rn in set(resnames):
        if rn in overrides:
            result[rn] = overrides[rn]
        elif rn in AMINO_ACID_RESNAMES:
            result[rn] = "protein"
        elif rn in DEFAULT_WATER_RESNAMES:
            result[rn] = "water"
        elif rn in DEFAULT_ION_RESNAMES:
            result[rn] = "ion"
        elif rn in DEFAULT_LIPID_RESNAMES:
            result[rn] = "lipid"
        else:
            result[rn] = "other"
    return result
