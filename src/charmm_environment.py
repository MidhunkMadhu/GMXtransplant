"""Versioned CHARMM-GUI environment names and GROMACS naming aliases.

The JSON snapshot contains every archive download ID and RESI identifier from
the sources listed in it. These names are facts, not copied parameter tables.
No finite catalogue covers future releases or user-defined molecules; unknown
names remain visible as inferred ligands. Never use an unrestricted POP* match.
"""
import json
from importlib.resources import files

CATALOGUE = json.loads(files("gmxtransplant").joinpath("charmm_environment.json").read_text())

# Water model labels and common exported coordinate/topology aliases. The GUI
# FAQ documents TIP3P, TIP4P-EW and OPC; the CHARMM/OpenMM water distribution
# also documents SPC/E, TIP4P/2005, TIP5P and Drude waters. A model label is not
# necessarily the exported residue name (usually TIP3, SOL or WAT).
ALIASES = {
    "water": "TIP TIP3 TIP3P TIPS3P T3P TP3M TIP4 TIP4P TIP4PEW TIP4P-EW T4P "
             "TIP4P2005 TIP4P/2005 TIP4P-FB TIP4PFB TIP5 TIP5P TIP5PE TIP5PEW "
             "SPC SPCE SPC/E SPC216 SOL WAT HOH H2O OH2 OPC OPC3 SWM4 SWM4DP "
             "SWM4NDP SWM4-NDP SWM6 W PW WF",
    "ion": "SOD NA NA+ POT K K+ CLA CL CL- LIT LI LI+ RUB RB RB+ CES CS CS+ "
           "CAL CA CA2 CA2+ MG MG2 MG2+ ZN ZN2 ZN2+ CD CD2 CD2+ BAR BA BA2 "
           "BA2+ SR SR2 SR2+ F F- BR BR- I I- CU CU1 CU2 FE FE2 FE3 MN MN2 "
           "CO CO2 NI NI2 AL AL3 AG AG+ AU AU1 HG HG2 PB PB2 CSOD CPOT CCLA",
    "sterol": "CHL CHL1 CHOL CLR ERG ERGO SITO STIG",
    "lipid": "POP PIP2 PIP3 CER CER1 CER2 BSM CDL2 TOCL",
}
ALIAS_SOURCES = [
    "https://www.charmm-gui.org/?doc=faq",
    "https://github.com/openmm/openmmforcefields/blob/main/charmm/files/waters.yaml",
    "https://manual.gromacs.org/current/user-guide/force-fields.html",
]


def environment_name(name):
    """Return (category, evidence) for an exact case-insensitive name, or None."""
    name = str(name).strip().upper()
    if name in CATALOGUE["names"]:
        item = CATALOGUE["names"][name]
        return item["category"], {
            "method": "official_catalogue", "matched_name": name,
            "sources": [CATALOGUE["sources"][i]["url"] for i in item["sources"]],
        }
    for category, words in ALIASES.items():
        if name in words.split():
            return category, {"method": "export_alias", "matched_name": name,
                              "sources": ALIAS_SOURCES}
    return None
