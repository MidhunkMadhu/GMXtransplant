"""
Pure-Python GROMACS topology (.top/.itp) assembly.

Builds a final topol.top and self-contained toppar/ directory describing
the final system after alignment, clash removal, and neutralization.

-------------------------------------------------------------------------
How this works, and the assumptions it makes
-------------------------------------------------------------------------
You point it at two existing, already-correct CHARMM-GUI GROMACS
directories:

  protein_toppar_dir     -- the toppar/ folder belonging to whichever
                              structure the NEW protein (and, usually,
                              its ligand) came from. Its OWN sibling
                              topol.top (one directory up) is read
                              automatically -- not configured separately
                              -- purely to learn which moleculetype name
                              corresponds to each protein chain, in the
                              correct order and exact ITP atom counts. This
                              supplies reliable chain boundaries even when
                              GRO has no chain IDs or residue numbers reset.

  environment_toppar_dir  -- a toppar/ folder (from any CHARMM-GUI system)
                              containing correct parameters for the kept
                              environment's lipid/water/ion species. Its
                              own topol.top is NOT read -- only the .itp
                              files themselves, each parsed for its own
                              [ moleculetype ] name and the resname(s) its
                              [ atoms ] section uses.

Everything else -- the actual final composition (how many of each protein
chain, how many lipids/waters/ions survived clash removal and
neutralization, whether a ligand is present) -- comes fresh from the
final assembled Universe, never from either template's own counts.

Moleculetype-name resolution:
  - A protein chain's moleculetype name is read positionally from
    protein_toppar_dir's sibling topol.top's own [ molecules ] section
    (filtered to entries whose ITP contains at least one standard amino
    acid). This includes protein moleculetypes with caps, patched termini,
    or nonstandard covalent residues. Each molecule instance is one chain.
  - A ligand/cofactor resname's moleculetype name is found by scanning
    every .itp under protein_toppar_dir, `ligand_itp_paths`, AND
    environment_toppar_dir for a moleculetype whose atoms use that
    resname. In protein mode, each `replacement_ligands[].itp_path` is
    added to `ligand_itp_paths` automatically, so it need not be repeated.
  - A lipid/ion resname is matched to the ITP whose [ atoms ] rows use that
    resname. Zero matches or multiple matches are errors unless an explicit
    `topology.moleculetype_overrides` entry resolves the case.
  - Water is the one well-known CHARMM-GUI naming exception: the
    coordinate file's resname (TIP3, TIP4, ...) usually does NOT match
    its own moleculetype name (commonly "SOL"). If no direct itp match is
    found for a water-classified resname, "SOL" is used as a fallback.
    Override anything here via `topology.moleculetype_overrides`.

Anything that doesn't resolve this way raises TopologyError with a
specific, actionable message -- this module never silently guesses a
mapping it isn't confident about.
-------------------------------------------------------------------------
"""

from __future__ import annotations

import os
import re
import shutil
import hashlib
import tempfile
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple

from classify import classify_resnames, AMINO_ACID_RESNAMES, DEFAULT_WATER_RESNAMES
from itp import ITPMoleculeType, ITPParseError, parse_itp


class TopologyError(Exception):
    pass


_INCLUDE_RE = re.compile(r'^\s*#include\s+"([^"]+)"')
_SECTION_RE = re.compile(r'^\s*\[\s*([A-Za-z_]+)\s*\]')


def _strip_comment(line: str) -> str:
    idx = line.find(";")
    return line if idx == -1 else line[:idx]


# ---------------------------------------------------------------------------
# .top / .itp parsing
# ---------------------------------------------------------------------------

def parse_top_includes(top_path: str) -> List[str]:
    """Return every #include target string, in file order, exactly as
    written (e.g. "toppar/charmm36.ff/forcefield.itp")."""
    includes = []
    with open(top_path) as fh:
        for line in fh:
            m = _INCLUDE_RE.match(line)
            if m:
                includes.append(m.group(1))
    return includes


def parse_top_molecules(top_path: str) -> List[Tuple[str, int]]:
    """Return the [ molecules ] section as an ordered list of
    (moleculetype_name, count)."""
    molecules = []
    in_section = False
    with open(top_path) as fh:
        for raw in fh:
            line = _strip_comment(raw).strip()
            if not line:
                continue
            if line.startswith("#"):
                directive = line.split(None, 1)[0].lower()
                if in_section and directive in {
                    "#if", "#ifdef", "#ifndef", "#elif", "#else", "#endif",
                }:
                    raise TopologyError(
                        f"Conditional preprocessor directive '{directive}' occurs "
                        f"inside [ molecules ] in '{top_path}'. Supply a topology "
                        "preprocessed with the active GROMACS defines."
                    )
                continue
            m = _SECTION_RE.match(line)
            if m:
                in_section = m.group(1).lower() == "molecules"
                continue
            if in_section:
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        count = int(parts[1])
                    except ValueError:
                        continue
                    molecules.append((parts[0], count))
    if not molecules:
        raise TopologyError(f"No [ molecules ] entries found in '{top_path}'")
    return molecules


def parse_itp_moleculetypes(itp_path: str) -> Dict[str, Set[str]]:
    """Parse one .itp file (it may define more than one [ moleculetype ]
    block) and return {moleculetype_name: {resnames used in its
    [ atoms ] section}}. Files with no [ moleculetype ] block (bonded /
    nonbonded parameter files, cmap tables, ...) return {}."""
    try:
        return {name: definition.resnames for name, definition in parse_itp(itp_path).items()}
    except ITPParseError as exc:
        raise TopologyError(str(exc)) from exc


def parse_itp_molecule_atom_counts(itp_path: str) -> Dict[str, int]:
    """Return the number of [ atoms ] rows for each moleculetype in an ITP.

    The count represents the number of atoms expected in one complete
    molecule. Parameter-only ITP files return an empty dictionary.
    """
    try:
        return {name: definition.atom_count for name, definition in parse_itp(itp_path).items()}
    except ITPParseError as exc:
        raise TopologyError(str(exc)) from exc


def scan_toppar_definitions(toppar_dir: str) -> Tuple[Dict[str, ITPMoleculeType], Dict[str, str]]:
    """Parse every molecule ITP under a directory in deterministic order.

    Repeated, byte-independent definitions of the same moleculetype are
    accepted. Conflicting definitions stop the run instead of allowing the
    filesystem traversal order to choose one.
    """
    definitions: Dict[str, ITPMoleculeType] = {}
    file_by_mtype: Dict[str, str] = {}
    for root, dirs, files in os.walk(toppar_dir):
        dirs.sort()
        for fn in sorted(files):
            if not fn.lower().endswith(".itp"):
                continue
            path = os.path.normpath(os.path.join(root, fn))
            try:
                parsed = parse_itp(path)
            except ITPParseError as exc:
                raise TopologyError(str(exc)) from exc
            for mtype, definition in parsed.items():
                if mtype in definitions:
                    previous = definitions[mtype]
                    if previous.content_signature() != definition.content_signature():
                        raise TopologyError(
                            f"Moleculetype '{mtype}' has conflicting definitions in "
                            f"'{file_by_mtype[mtype]}' and '{path}'"
                        )
                    continue
                definitions[mtype] = definition
                file_by_mtype[mtype] = path
    return definitions, file_by_mtype


def scan_toppar_dir(toppar_dir: str) -> Tuple[Dict[str, Set[str]], Dict[str, str]]:
    definitions, file_by_mtype = scan_toppar_definitions(toppar_dir)
    return ({name: definition.resnames for name, definition in definitions.items()},
            file_by_mtype)


def scan_toppar_atom_counts(toppar_dir: str) -> Tuple[Dict[str, int], Dict[str, str]]:
    """Recursively read molecule atom counts from every ITP in toppar_dir."""
    definitions, file_by_mtype = scan_toppar_definitions(toppar_dir)
    return ({name: definition.atom_count for name, definition in definitions.items()},
            file_by_mtype)


def _parse_replacement_itp_definitions(replacement_ligand_itp_paths):
    """Parse ligand ITPs that intentionally replace template definitions.

    Replacement ITPs are authoritative for every moleculetype they define.
    They are parsed as a separate pool so that an old ligand definition found
    while scanning the protein/environment templates cannot win merely due to
    traversal order. Conflicts between replacement ITPs remain errors.
    """
    definitions: Dict[str, ITPMoleculeType] = {}
    files: Dict[str, str] = {}
    realpaths: Set[str] = set()
    for path in replacement_ligand_itp_paths or ():
        normalized = os.path.realpath(path)
        if normalized in realpaths:
            continue
        realpaths.add(normalized)
        try:
            parsed = parse_itp(path)
        except ITPParseError as exc:
            raise TopologyError(str(exc)) from exc
        for mtype, definition in parsed.items():
            if mtype in definitions:
                if definitions[mtype].content_signature() != definition.content_signature():
                    raise TopologyError(
                        f"Replacement ligand moleculetype '{mtype}' has conflicting "
                        f"definitions in '{files[mtype]}' and '{path}'."
                    )
                continue
            definitions[mtype] = definition
            files[mtype] = path
    return definitions, files, realpaths


def collect_topology_definitions(
    topology_cfg, replacement_ligand_itp_paths=()
):
    """Collect all configured molecule definitions with conflict detection.

    A ligand ITP explicitly identified as a replacement shadows a definition
    with the same moleculetype name in either template directory. All other
    conflicting definitions still stop the run.
    """
    definitions: Dict[str, ITPMoleculeType] = {}
    files: Dict[str, str] = {}
    replacement_defs, replacement_files, replacement_realpaths = (
        _parse_replacement_itp_definitions(replacement_ligand_itp_paths)
    )
    replaced_mtypes = set(replacement_defs)

    def merge(new_definitions, new_files):
        for mtype, definition in new_definitions.items():
            if mtype in replaced_mtypes:
                continue
            if mtype in definitions:
                if definitions[mtype].content_signature() != definition.content_signature():
                    raise TopologyError(
                        f"Moleculetype '{mtype}' has conflicting definitions in "
                        f"'{files[mtype]}' and '{new_files[mtype]}'."
                    )
                continue
            definitions[mtype] = definition
            files[mtype] = new_files[mtype]

    seen_dirs = set()
    for directory in (
        topology_cfg.protein_toppar_dir,
        topology_cfg.environment_toppar_dir,
    ):
        normalized = os.path.realpath(directory)
        if normalized in seen_dirs:
            continue
        seen_dirs.add(normalized)
        new_definitions, new_files = scan_toppar_definitions(directory)
        merge(new_definitions, new_files)

    for path in topology_cfg.ligand_itp_paths:
        if os.path.realpath(path) in replacement_realpaths:
            continue
        try:
            parsed = parse_itp(path)
        except ITPParseError as exc:
            raise TopologyError(str(exc)) from exc
        merge(parsed, {name: path for name in parsed})
    definitions.update(replacement_defs)
    files.update(replacement_files)
    return definitions, files


@dataclass
class TopologyChargeModel:
    layout: "ProteinChainLayout"
    definitions: Dict[str, ITPMoleculeType]
    files: Dict[str, str]
    resname_to_mtype: Dict[str, str]
    unresolved_resnames: List[str]

    @property
    def residue_charge_table(self) -> Dict[str, float]:
        return {
            resname: self.definitions[mtype].charge
            for resname, mtype in self.resname_to_mtype.items()
        }


def build_topology_charge_model(
    topology_cfg, nonprotein_resnames, replacement_ligand_itp_paths=()
) -> TopologyChargeModel:
    """Resolve exact ITP charge sources for every requested nonprotein name."""
    layout = build_protein_chain_layout(topology_cfg)
    definitions, files = collect_topology_definitions(
        topology_cfg, replacement_ligand_itp_paths
    )
    protein_mtypes = set(layout.moleculetypes)
    resolved: Dict[str, str] = {}
    unresolved: List[str] = []
    for resname in sorted(set(str(name) for name in nonprotein_resnames)):
        override = topology_cfg.moleculetype_overrides.get(resname)
        if override:
            definition = definitions.get(override)
            if definition is None:
                raise TopologyError(
                    f"topology.moleculetype_overrides maps {resname} to '{override}', "
                    "but no configured ITP defines that moleculetype."
                )
            if resname not in definition.resnames:
                raise TopologyError(
                    f"topology.moleculetype_overrides maps {resname} to '{override}', "
                    f"but its ITP uses residue names {sorted(definition.resnames)}."
                )
            resolved[resname] = override
            continue

        candidates = sorted(
            mtype for mtype, definition in definitions.items()
            if mtype not in protein_mtypes and resname in definition.resnames
        )
        if not candidates and resname in DEFAULT_WATER_RESNAMES and "SOL" in definitions:
            candidates = ["SOL"]
        if len(candidates) > 1:
            raise TopologyError(
                f"Residue name '{resname}' occurs in multiple nonprotein "
                f"moleculetypes: {candidates}. Set topology.moleculetype_overrides "
                "to choose explicitly."
            )
        if not candidates:
            unresolved.append(resname)
        else:
            resolved[resname] = candidates[0]
    return TopologyChargeModel(layout, definitions, files, resolved, unresolved)


def find_forcefield_include(top_path: str) -> str:
    """Return the #include target string (as literally written) whose
    basename is 'forcefield.itp' -- the standard CHARMM-GUI convention for
    the shared force-field parameter entry point."""
    for inc in parse_top_includes(top_path):
        if os.path.basename(inc).lower() == "forcefield.itp":
            return inc
    raise TopologyError(
        f"Could not find a '#include \".../forcefield.itp\"' line in '{top_path}' -- "
        f"is this a real CHARMM-GUI GROMACS topol.top?"
    )


# ---------------------------------------------------------------------------
# Final-system molecule sequence (from the assembled Universe, not from
# either template's own counts)
# ---------------------------------------------------------------------------

@dataclass
class MoleculeEntry:
    kind: str            # "protein" | "ligand" | "environment"
    resname: str          # representative resname (chain: first residue's resname)
    count: int            # number of molecule instances in this run
    chain_index: Optional[int] = None   # only set for kind == "protein"


@dataclass
class ProteinChainLayout:
    """Protein molecule instances from the protein reference topology."""

    template_top: str
    moleculetypes: List[str]
    atom_counts: List[int]
    charges: List[float]
    definitions: List[ITPMoleculeType]

    @property
    def total_atoms(self) -> int:
        return sum(self.atom_counts)


def protein_template_top_path(topology_cfg) -> str:
    if topology_cfg.protein_template_top:
        return topology_cfg.protein_template_top
    return os.path.join(
        os.path.dirname(os.path.normpath(topology_cfg.protein_toppar_dir)),
        "topol.top",
    )


def build_protein_chain_layout(topology_cfg) -> ProteinChainLayout:
    """Read protein molecule boundaries from topol.top and protein ITPs.

    A moleculetype is a protein when its ITP contains at least one
    recognized amino-acid residue. This deliberately allows the same molecule
    to contain caps, patched termini, or non-standard covalent residues.
    """
    top_path = protein_template_top_path(topology_cfg)
    if not os.path.isfile(top_path):
        raise TopologyError(f"Protein reference topology does not exist: '{top_path}'")
    definitions, _files = scan_toppar_definitions(topology_cfg.protein_toppar_dir)
    moleculetypes: List[str] = []
    atom_counts: List[int] = []
    charges: List[float] = []
    chain_definitions: List[ITPMoleculeType] = []
    for mtype, count in parse_top_molecules(top_path):
        definition = definitions.get(mtype)
        if definition is None or not (definition.resnames & AMINO_ACID_RESNAMES):
            continue
        for _ in range(count):
            moleculetypes.append(mtype)
            atom_counts.append(definition.atom_count)
            charges.append(definition.charge)
            chain_definitions.append(definition)
    if not moleculetypes:
        raise TopologyError(
            f"No protein moleculetype was found by comparing '{top_path}' "
            f"with ITP files under '{topology_cfg.protein_toppar_dir}'."
        )
    return ProteinChainLayout(
        template_top=top_path,
        moleculetypes=moleculetypes,
        atom_counts=atom_counts,
        charges=charges,
        definitions=chain_definitions,
    )


def _validate_protein_prefix(final_universe, layout: ProteinChainLayout) -> int:
    """Validate topology-derived protein chains and return residue count used."""
    if len(final_universe.atoms) < layout.total_atoms:
        raise TopologyError(
            f"Final coordinates contain {len(final_universe.atoms)} atoms, fewer than "
            f"the {layout.total_atoms} protein atoms described by "
            f"{layout.moleculetypes}."
        )
    atom_offset = 0
    for chain_index, (mtype, definition) in enumerate(
        zip(layout.moleculetypes, layout.definitions), 1
    ):
        stop = atom_offset + definition.atom_count
        actual = final_universe.atoms[atom_offset:stop]
        actual_names = [str(name) for name in actual.names]
        expected_names = list(definition.atom_names)
        if actual_names != expected_names:
            mismatch = next(
                i for i, (got, expected) in enumerate(zip(actual_names, expected_names))
                if got != expected
            )
            raise TopologyError(
                f"Protein chain {chain_index} ({mtype}) does not match its ITP atom "
                f"order at molecule atom {mismatch + 1}: coordinates have "
                f"'{actual_names[mismatch]}', ITP expects '{expected_names[mismatch]}'. "
                "The final protein cannot be assigned topology-derived chain boundaries."
            )
        if stop < len(final_universe.atoms):
            if final_universe.atoms[stop - 1].resindex == final_universe.atoms[stop].resindex:
                raise TopologyError(
                    f"Topology-derived boundary after protein chain {chain_index} "
                    f"({mtype}, atom {stop}) falls inside coordinate residue "
                    f"{final_universe.atoms[stop].resname} "
                    f"{final_universe.atoms[stop].resid}."
                )
        atom_offset = stop

    if atom_offset == 0:
        return 0
    protein_resindices = set(final_universe.atoms[:atom_offset].resindices.tolist())
    return len(protein_resindices)


def determine_final_molecule_sequence(
    final_universe, layout: ProteinChainLayout
) -> Tuple[List[MoleculeEntry], int]:
    """Walk final_universe's residues in file order and group them into
    molecule instances: one instance per contiguous protein chain, one
    instance per contiguous same-residue ligand/cofactor run, one instance
    per contiguous same-resname lipid/water/ion run. Returns (entries,
    n_protein_chains)."""
    n_protein_residues = _validate_protein_prefix(final_universe, layout)
    resnames = final_universe.residues.resnames
    resname_class = classify_resnames(resnames)

    n = len(resnames)
    entries: List[MoleculeEntry] = [
        MoleculeEntry("protein", definition.atoms[0].resname, 1, chain_index=i)
        for i, definition in enumerate(layout.definitions)
    ]
    i = n_protein_residues
    while i < n:
        cls = resname_class[resnames[i]]
        if cls == "protein":
            previous = final_universe.residues[i - 1] if i else None
            current = final_universe.residues[i]
            raise TopologyError(
                "A protein-classified residue occurs outside the topology-derived "
                f"protein prefix at final residue position {i + 1}: "
                f"{current.resname}{current.resid}. The preceding boundary is "
                + (f"after {previous.resname}{previous.resid}." if previous else "at atom 0.")
            )
        if cls == "other":
            entries.append(MoleculeEntry("ligand", resnames[i], 1))
            i += 1
            continue
        # lipid / water / ion
        j = i
        while j < n and resnames[j] == resnames[i] and resname_class[resnames[j]] == cls:
            j += 1
        entries.append(MoleculeEntry("environment", resnames[i], j - i))
        i = j
    return entries, len(layout.moleculetypes)


# ---------------------------------------------------------------------------
# Moleculetype resolution
# ---------------------------------------------------------------------------

@dataclass
class ResolvedTopology:
    protein_chain_mtypes: List[str]          # in chain-index order
    ligand_mtype_for_resname: Dict[str, str]
    environment_mtype_for_resname: Dict[str, str]
    mtype_file: Dict[str, str]                # moleculetype -> absolute itp path
    forcefield_dir: str                       # absolute path of the forcefield.itp's own directory
    forcefield_itp_relname: str                # e.g. "charmm36.ff/forcefield.itp" -- path *within* forcefield_dir's parent
    forcefield_is_flat: bool                   # copy entry point/dependencies, not its whole parent directory
    notes: List[str] = field(default_factory=list)


def resolve_topology(
    protein_toppar_dir: str,
    environment_toppar_dir: str,
    ligand_itp_paths: List[str],
    final_ligand_resnames: Set[str],
    n_protein_chains_final: int,
    moleculetype_overrides: Dict[str, str],
    protein_template_top: Optional[str] = None,
    replacement_ligand_itp_paths=(),
) -> ResolvedTopology:
    if not protein_template_top:
        parent = os.path.dirname(os.path.normpath(protein_toppar_dir))
        protein_template_top = os.path.join(parent, "topol.top")
    if not os.path.isfile(protein_template_top):
        raise TopologyError(
            f"Could not find a topol.top next to protein_toppar_dir "
            f"(looked for '{protein_template_top}'). Standard CHARMM-GUI layout puts "
            f"topol.top one directory above its own toppar/ folder -- if yours differs, "
            f"set topology.protein_template_top explicitly."
        )

    replacement_defs, replacement_files, replacement_realpaths = (
        _parse_replacement_itp_definitions(replacement_ligand_itp_paths)
    )
    replaced_mtypes = set(replacement_defs)

    recv_defs, recv_file_by_mtype = scan_toppar_definitions(protein_toppar_dir)
    for itp_path in ligand_itp_paths:
        if os.path.realpath(itp_path) in replacement_realpaths:
            continue
        if not os.path.isfile(itp_path):
            raise TopologyError(f"topology.ligand_itp_paths entry not found: '{itp_path}'")
        try:
            extra_defs = parse_itp(itp_path)
        except ITPParseError as exc:
            raise TopologyError(str(exc)) from exc
        for mtype, definition in extra_defs.items():
            if mtype in replaced_mtypes:
                continue
            if mtype in recv_defs and (
                recv_defs[mtype].content_signature() != definition.content_signature()
            ):
                raise TopologyError(
                    f"Moleculetype '{mtype}' has conflicting definitions in "
                    f"'{recv_file_by_mtype[mtype]}' and '{itp_path}'"
                )
            recv_defs.setdefault(mtype, definition)
            recv_file_by_mtype.setdefault(mtype, itp_path)

    env_defs, env_file_by_mtype = scan_toppar_definitions(environment_toppar_dir)
    notes: List[str] = []
    for mtype in sorted(replaced_mtypes):
        old_protein_definition = recv_defs.get(mtype)
        if old_protein_definition is not None and (
            old_protein_definition.resnames & AMINO_ACID_RESNAMES
        ):
            raise TopologyError(
                f"Replacement ligand ITP '{replacement_files[mtype]}' attempts to "
                f"shadow protein moleculetype '{mtype}' from "
                f"'{recv_file_by_mtype[mtype]}'."
            )

        superseded_paths = {
            source_files[mtype]
            for source_defs, source_files in (
                (recv_defs, recv_file_by_mtype),
                (env_defs, env_file_by_mtype),
            )
            if mtype in source_defs
            and os.path.realpath(source_files[mtype])
            != os.path.realpath(replacement_files[mtype])
        }
        recv_defs.pop(mtype, None)
        recv_file_by_mtype.pop(mtype, None)
        env_defs.pop(mtype, None)
        env_file_by_mtype.pop(mtype, None)
        recv_defs[mtype] = replacement_defs[mtype]
        recv_file_by_mtype[mtype] = replacement_files[mtype]
        if superseded_paths:
            notes.append(
                f"Replacement ligand moleculetype '{mtype}' from "
                f"'{replacement_files[mtype]}' superseded template definition(s): "
                + ", ".join(f"'{path}'" for path in sorted(superseded_paths))
                + "."
            )

    for mtype in sorted(set(recv_defs) & set(env_defs)):
        if recv_defs[mtype].content_signature() != env_defs[mtype].content_signature():
            raise TopologyError(
                f"Moleculetype '{mtype}' conflicts between protein ITP "
                f"'{recv_file_by_mtype[mtype]}' and environment ITP "
                f"'{env_file_by_mtype[mtype]}'."
            )
    recv_resnames_by_mtype = {
        name: definition.resnames for name, definition in recv_defs.items()
    }
    env_resnames_by_mtype = {
        name: definition.resnames for name, definition in env_defs.items()
    }

    # --- protein chain moleculetypes, positional, from the protein's own
    #     topol.top [ molecules ] section -----------------------------------
    protein_molecules = parse_top_molecules(protein_template_top)
    protein_chain_mtypes: List[str] = []
    for mtype, count in protein_molecules:
        resnames_here = recv_resnames_by_mtype.get(mtype)
        if resnames_here and (resnames_here & AMINO_ACID_RESNAMES):
            protein_chain_mtypes.extend([mtype] * count)
    if len(protein_chain_mtypes) != n_protein_chains_final:
        raise TopologyError(
            f"Protein chain count mismatch: the final assembled system has "
            f"{n_protein_chains_final} contiguous protein chain(s), but "
            f"'{protein_template_top}' [ molecules ] lists {len(protein_chain_mtypes)} "
            f"protein moleculetype instance(s) ({protein_chain_mtypes}) found in "
            f"'{protein_toppar_dir}'. This usually means protein_toppar_dir doesn't "
            "belong to the protein in the final system. Detected topology boundaries: "
            + "; ".join(
                f"chain {i + 1} {mtype} ({recv_defs[mtype].atom_count} atoms)"
                for i, mtype in enumerate(protein_chain_mtypes)
            )
        )

    # --- ligand moleculetype(s), by resname, searched across BOTH pools ----
    combined_resnames_by_mtype: Dict[str, Set[str]] = {}
    combined_file_by_mtype: Dict[str, str] = {}
    for d in (recv_resnames_by_mtype, env_resnames_by_mtype):
        for mtype, resnames in d.items():
            combined_resnames_by_mtype.setdefault(mtype, set()).update(resnames)
    for d in (recv_file_by_mtype, env_file_by_mtype):
        combined_file_by_mtype.update(d)

    ligand_mtype_for_resname: Dict[str, str] = {}
    for rn in sorted(final_ligand_resnames):
        if rn in moleculetype_overrides:
            override = moleculetype_overrides[rn]
            if override not in combined_resnames_by_mtype:
                raise TopologyError(
                    f"topology.moleculetype_overrides maps {rn} to '{override}', "
                    "but no configured ITP defines that moleculetype."
                )
            if rn not in combined_resnames_by_mtype[override]:
                raise TopologyError(
                    f"topology.moleculetype_overrides maps {rn} to '{override}', "
                    f"but its ITP uses residue names "
                    f"{sorted(combined_resnames_by_mtype[override])}."
                )
            ligand_mtype_for_resname[rn] = override
            continue
        candidates = []
        for mtype, resnames in combined_resnames_by_mtype.items():
            if mtype in protein_chain_mtypes:
                continue
            if rn in resnames:
                candidates.append(mtype)
        if not candidates:
            raise TopologyError(
                f"Could not find a moleculetype for ligand/cofactor resname '{rn}' in "
                f"protein_toppar_dir, topology.ligand_itp_paths, or environment_toppar_dir. "
                f"Declare it under replacement_ligands (protein mode), add its .itp file "
                f"to topology.ligand_itp_paths, or add an explicit "
                f"topology.moleculetype_overrides: {{{rn}: <moleculetype_name>}}."
            )
        if len(candidates) > 1:
            raise TopologyError(
                f"Ligand/cofactor resname '{rn}' occurs in multiple moleculetype "
                f"definitions: {sorted(candidates)}. Set "
                f"topology.moleculetype_overrides: {{{rn}: <moleculetype_name>}} "
                "instead of relying on file traversal order."
            )
        ligand_mtype_for_resname[rn] = candidates[0]

    # --- environment (lipid/water/ion) moleculetypes, by resname -----------
    environment_mtype_for_resname: Dict[str, str] = {}

    def resolve_environment_resname(rn: str, molclass: str) -> str:
        if rn in moleculetype_overrides:
            override = moleculetype_overrides[rn]
            if override not in env_resnames_by_mtype:
                raise TopologyError(
                    f"topology.moleculetype_overrides maps environment residue "
                    f"{rn} to '{override}', but environment_toppar_dir does not "
                    "define that moleculetype."
                )
            if rn not in env_resnames_by_mtype[override]:
                raise TopologyError(
                    f"topology.moleculetype_overrides maps {rn} to '{override}', "
                    f"but its ITP uses residue names "
                    f"{sorted(env_resnames_by_mtype[override])}."
                )
            return override
        candidates = [
            mtype for mtype, resnames in env_resnames_by_mtype.items()
            if rn in resnames
        ]
        if len(candidates) == 1:
            return candidates[0]
        if len(candidates) > 1:
            raise TopologyError(
                f"Environment resname '{rn}' occurs in multiple moleculetype "
                f"definitions: {sorted(candidates)}. Set "
                f"topology.moleculetype_overrides: {{{rn}: <moleculetype_name>}}."
            )
        if molclass == "water" or rn in DEFAULT_WATER_RESNAMES:
            if "SOL" in env_resnames_by_mtype:
                notes.append(
                    f"Water resname '{rn}' had no direct itp match in environment_toppar_dir; "
                    f"used the standard CHARMM-GUI fallback moleculetype 'SOL'."
                )
                return "SOL"
        raise TopologyError(
            f"Could not find a moleculetype for environment resname '{rn}' anywhere under "
            f"'{environment_toppar_dir}'. Add an explicit "
            f"topology.moleculetype_overrides: {{{rn}: <moleculetype_name>}}."
        )

    # --- forcefield entry point ---------------------------------------------
    ff_include = find_forcefield_include(protein_template_top)
    ff_abs = os.path.normpath(os.path.join(os.path.dirname(protein_template_top), ff_include))
    if not os.path.isfile(ff_abs):
        raise TopologyError(
            f"'{protein_template_top}' references forcefield '{ff_include}' but "
            f"'{ff_abs}' does not exist on disk."
        )
    forcefield_is_flat = (
        os.path.realpath(os.path.dirname(ff_abs))
        == os.path.realpath(protein_toppar_dir)
    )

    # CHARMM-GUI folds ligand-specific atom and bonded parameters into the
    # toppar/forcefield.itp file rather than the ligand molecule ITP. When a
    # replacement ITP has a sibling forcefield.itp, the two files are one
    # parameter set and must be selected together. Standalone replacement
    # ITPs without such a sibling continue to use the protein force field.
    replacement_forcefields: Dict[str, str] = {}
    for replacement_path in sorted(set(replacement_files.values())):
        candidate = os.path.join(os.path.dirname(replacement_path), "forcefield.itp")
        if not os.path.isfile(candidate):
            continue
        digest = _file_sha256(candidate)
        replacement_forcefields.setdefault(digest, candidate)
    if len(replacement_forcefields) > 1:
        raise TopologyError(
            "Replacement ligand ITPs have different sibling forcefield.itp "
            "files. A generated topology can use only one force-field parameter "
            "set: "
            + ", ".join(
                f"'{path}'" for path in sorted(replacement_forcefields.values())
            )
            + "."
        )
    if replacement_forcefields:
        replacement_ff_abs = next(iter(replacement_forcefields.values()))
        if os.path.realpath(replacement_ff_abs) != os.path.realpath(ff_abs):
            notes.append(
                f"Selected replacement ligand force field '{replacement_ff_abs}' "
                f"instead of template force field '{ff_abs}'."
            )
        ff_abs = replacement_ff_abs
        forcefield_is_flat = True

    resolved = ResolvedTopology(
        protein_chain_mtypes=protein_chain_mtypes,
        ligand_mtype_for_resname=ligand_mtype_for_resname,
        environment_mtype_for_resname=environment_mtype_for_resname,
        mtype_file=combined_file_by_mtype,
        forcefield_dir=os.path.dirname(ff_abs),
        forcefield_itp_relname=os.path.basename(ff_abs),
        forcefield_is_flat=forcefield_is_flat,
        notes=notes,
    )
    # stash the resolver closure's results lazily -- environment resnames are
    # resolved on demand in assemble_topology() once we know exactly which
    # ones appear in the final system; expose the helper via the object.
    resolved._resolve_environment_resname = resolve_environment_resname  # type: ignore[attr-defined]
    resolved._recv_file_by_mtype = recv_file_by_mtype                    # type: ignore[attr-defined]
    resolved._env_file_by_mtype = env_file_by_mtype                      # type: ignore[attr-defined]
    return resolved


# ---------------------------------------------------------------------------
# Assembly + writing
# ---------------------------------------------------------------------------

@dataclass
class TopologyResult:
    top_path: str
    toppar_dir: str
    molecules_written: List[Tuple[str, int]]
    includes_written: List[str]
    notes: List[str] = field(default_factory=list)
    charge_audit: Dict[str, object] = field(default_factory=dict)


def _collect_forcefield_dependencies(ff_abs_path: str, boundary_dir: str) -> List[str]:
    """Recursively resolve every #include forcefield.itp itself needs
    (e.g. ffbonded.itp/ffnonbonded.itp) for the case where forcefield.itp
    sits directly in the same toppar/ folder as the per-molecule .itp
    files rather than its own separate subfolder -- there's no single
    self-contained directory to copy wholesale in that layout, so this
    walks forcefield.itp's own #include chain instead and returns the
    absolute paths to copy (forcefield.itp first, then whatever it needs,
    each listed once). Stays within boundary_dir -- an #include that
    resolves outside of it, or to a file that doesn't exist, is skipped
    rather than followed (a moleculetype .itp is never legitimately
    #included BY forcefield.itp in standard GROMACS layouts, so this
    should only ever pick up genuine parameter-only files)."""
    boundary_norm = os.path.normpath(boundary_dir) + os.sep
    seen: List[str] = []
    seen_set = set()

    def _walk(path: str) -> None:
        path = os.path.normpath(path)
        if path in seen_set or not os.path.isfile(path):
            return
        seen_set.add(path)
        seen.append(path)
        for inc in parse_top_includes(path):
            inc_abs = os.path.normpath(os.path.join(os.path.dirname(path), inc))
            if not inc_abs.startswith(boundary_norm):
                continue
            _walk(inc_abs)

    _walk(ff_abs_path)
    return seen


def _copy_file_flat(src: str, dest_dir: str, notes: List[str]) -> str:
    """Copy src into dest_dir (flat), returning the basename used. Adds a
    disambiguating prefix on a genuine name collision (same basename,
    different content already copied from a different source)."""
    os.makedirs(dest_dir, exist_ok=True)
    base = os.path.basename(src)
    dest = os.path.join(dest_dir, base)
    if os.path.exists(dest):
        src_digest = _file_sha256(src)
        if _file_sha256(dest) == src_digest:
            return base  # already copied (e.g. same file needed by two moleculetypes)
        stem, ext = os.path.splitext(base)
        digest_length = 12
        while True:
            base = f"{stem}__{src_digest[:digest_length]}{ext}"
            dest = os.path.join(dest_dir, base)
            if not os.path.exists(dest):
                break
            if _file_sha256(dest) == src_digest:
                return base
            digest_length += 4
            if digest_length > len(src_digest):
                raise TopologyError(
                    f"Could not create a collision-free deterministic filename for '{src}'."
                )
        notes.append(f"Renamed '{os.path.basename(src)}' -> '{base}' in toppar/ to avoid a filename collision.")
    shutil.copy2(src, dest)
    return base


def _file_sha256(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _remove_path(path: str) -> None:
    if os.path.isdir(path) and not os.path.islink(path):
        shutil.rmtree(path)
    else:
        try:
            os.unlink(path)
        except FileNotFoundError:
            pass


def _promote_topology_bundle(
    staged_top: str,
    staged_toppar: str,
    final_top: str,
    final_toppar: str,
    staging_root: str,
) -> None:
    """Replace topol.top and toppar/ together, restoring old files on error."""
    backup_top = os.path.join(staging_root, "previous.topol.top")
    backup_toppar = os.path.join(staging_root, "previous.toppar")
    top_backed_up = False
    toppar_backed_up = False
    top_promoted = False
    toppar_promoted = False
    try:
        if os.path.lexists(final_top):
            os.replace(final_top, backup_top)
            top_backed_up = True
        if os.path.lexists(final_toppar):
            os.replace(final_toppar, backup_toppar)
            toppar_backed_up = True
        os.replace(staged_toppar, final_toppar)
        toppar_promoted = True
        os.replace(staged_top, final_top)
        top_promoted = True
    except OSError as exc:
        try:
            if top_promoted:
                _remove_path(final_top)
            if toppar_promoted:
                _remove_path(final_toppar)
            if toppar_backed_up:
                os.replace(backup_toppar, final_toppar)
            if top_backed_up:
                os.replace(backup_top, final_top)
        except OSError as rollback_exc:
            raise TopologyError(
                "Topology promotion failed and the previous output could not be "
                f"fully restored: promotion={exc}; rollback={rollback_exc}"
            ) from rollback_exc
        raise TopologyError(
            f"Topology promotion failed; previous topol.top/toppar were restored: {exc}"
        ) from exc


def _stage_topology_files(resolved, used_mtypes, molecules, notes, staging_root):
    """Build a complete topology bundle below a temporary staging directory."""
    toppar_out = os.path.join(staging_root, "toppar")
    os.makedirs(toppar_out, exist_ok=True)
    includes_written: List[str] = []
    seen_src_paths: Dict[str, str] = {}

    ff_dir_norm = os.path.normpath(resolved.forcefield_dir)
    if resolved.forcefield_is_flat:
        ff_abs = os.path.join(
            resolved.forcefield_dir, resolved.forcefield_itp_relname
        )
        for dep in _collect_forcefield_dependencies(
            ff_abs, resolved.forcefield_dir
        ):
            if dep in seen_src_paths:
                continue
            base = _copy_file_flat(dep, toppar_out, notes)
            rel = f"toppar/{base}"
            seen_src_paths[dep] = rel
            includes_written.append(rel)
    else:
        ff_dest_dirname = os.path.basename(ff_dir_norm)
        ff_dest_dir = os.path.join(toppar_out, ff_dest_dirname)
        shutil.copytree(resolved.forcefield_dir, ff_dest_dir)
        includes_written.append(
            f"toppar/{ff_dest_dirname}/{resolved.forcefield_itp_relname}"
        )

    for mtype in used_mtypes:
        src = resolved.mtype_file.get(mtype)
        if not src:
            raise TopologyError(
                f"Internal error: moleculetype '{mtype}' was selected for "
                "[ molecules ] but no source .itp file is on record for it."
            )
        if src in seen_src_paths:
            continue
        base = _copy_file_flat(src, toppar_out, notes)
        rel = f"toppar/{base}"
        seen_src_paths[src] = rel
        includes_written.append(rel)

    top_path = os.path.join(staging_root, "topol.top")
    with open(top_path, "w", encoding="utf-8") as fh:
        for inc in includes_written:
            fh.write(f'#include "{inc}"\n')
        fh.write(
            "\n[ system ]\n; Name\nGenerated by the molecular-system "
            "modification pipeline\n\n"
        )
        fh.write("[ molecules ]\n; Compound        #mols\n")
        for mtype, count in molecules:
            fh.write(f"{mtype:<16s} {count}\n")
    return top_path, toppar_out, includes_written


def assemble_topology(
    final_universe, topology_cfg, replacement_ligand_itp_paths=()
) -> TopologyResult:
    layout = build_protein_chain_layout(topology_cfg)
    entries, n_chains = determine_final_molecule_sequence(final_universe, layout)
    final_ligand_resnames = {e.resname for e in entries if e.kind == "ligand"}

    resolved = resolve_topology(
        protein_toppar_dir=topology_cfg.protein_toppar_dir,
        environment_toppar_dir=topology_cfg.environment_toppar_dir,
        ligand_itp_paths=topology_cfg.ligand_itp_paths,
        final_ligand_resnames=final_ligand_resnames,
        n_protein_chains_final=n_chains,
        moleculetype_overrides=topology_cfg.moleculetype_overrides,
        protein_template_top=topology_cfg.protein_template_top or None,
        replacement_ligand_itp_paths=replacement_ligand_itp_paths,
    )
    notes = list(resolved.notes)

    # Resolve environment resnames on demand, now that we know exactly which
    # ones appear in the final system.
    env_resname_class = classify_resnames({e.resname for e in entries if e.kind == "environment"})
    for e in entries:
        if e.kind == "environment" and e.resname not in resolved.environment_mtype_for_resname:
            resolved.environment_mtype_for_resname[e.resname] = resolved._resolve_environment_resname(  # type: ignore[attr-defined]
                e.resname, env_resname_class[e.resname]
            )

    # --- build [ molecules ], merging adjacent identical moleculetypes -----
    molecules: List[Tuple[str, int]] = []
    for e in entries:
        if e.kind == "protein":
            mtype = resolved.protein_chain_mtypes[e.chain_index]
        elif e.kind == "ligand":
            mtype = resolved.ligand_mtype_for_resname[e.resname]
        else:
            mtype = resolved.environment_mtype_for_resname[e.resname]
        if molecules and molecules[-1][0] == mtype:
            molecules[-1] = (mtype, molecules[-1][1] + e.count)
        else:
            molecules.append((mtype, e.count))

    used_mtypes: List[str] = []
    for mtype, _count in molecules:
        if mtype not in used_mtypes:
            used_mtypes.append(mtype)

    # --- copy files, build #include list ------------------------------------
    output_dir = os.path.abspath(topology_cfg.output_dir)
    final_toppar = os.path.join(output_dir, "toppar")
    final_top = os.path.join(output_dir, "topol.top")
    toppar_out_real = os.path.realpath(final_toppar)
    input_toppar_dirs = {
        os.path.realpath(topology_cfg.protein_toppar_dir),
        os.path.realpath(topology_cfg.environment_toppar_dir),
    }
    if toppar_out_real in input_toppar_dirs:
        raise TopologyError(
            f"Refusing to rebuild output toppar directory '{final_toppar}' because "
            "it is also configured as an input topology directory. Choose a "
            "different topology.output_dir."
        )
    template_top = protein_template_top_path(topology_cfg)
    if os.path.realpath(final_top) == os.path.realpath(template_top):
        raise TopologyError(
            f"Refusing to overwrite protein template topology '{template_top}'. "
            "Choose a different topology.output_dir."
        )
    os.makedirs(output_dir, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix=".gmxtransplant-topology-", dir=output_dir
    ) as staging_root:
        top_path, toppar_out, includes_written = _stage_topology_files(
            resolved, used_mtypes, molecules, notes, staging_root
        )
        staged_result = TopologyResult(
            top_path=top_path,
            toppar_dir=toppar_out,
            molecules_written=molecules,
            includes_written=includes_written,
            notes=notes,
        )
        audit_final_topology(final_universe, staged_result)
        _promote_topology_bundle(
            top_path, toppar_out, final_top, final_toppar, staging_root
        )
    notes.append(
        "Topology bundle was built and audited in temporary storage, then "
        "transactionally promoted; an existing bundle is restored on promotion failure."
    )
    return TopologyResult(
        top_path=final_top,
        toppar_dir=final_toppar,
        molecules_written=molecules,
        includes_written=includes_written,
        notes=notes,
    )


def audit_final_topology(final_universe, topology_result: TopologyResult) -> Dict[str, object]:
    """Verify that generated topology molecules match coordinate atom order.

    The audit deliberately reparses the written ``topol.top`` and copied ITP
    files.  This makes it independent of the in-memory resolver used during
    assembly and ensures the files handed to GROMACS describe exactly the
    coordinate sequence that will be written.
    """
    molecules = parse_top_molecules(topology_result.top_path)
    definitions, definition_paths = scan_toppar_definitions(
        topology_result.toppar_dir
    )

    expected_atoms = []
    molecule_charges: Dict[str, float] = {}
    net_charge = 0.0
    for mtype, count in molecules:
        if count < 1:
            raise TopologyError(
                f"Generated topology has invalid count {count} for "
                f"moleculetype '{mtype}'."
            )
        definition = definitions.get(mtype)
        if definition is None:
            raise TopologyError(
                f"Generated topology references moleculetype '{mtype}', but no "
                f"matching definition exists under '{topology_result.toppar_dir}'."
            )
        molecule_charges[mtype] = definition.charge
        net_charge += definition.charge * count
        atom_signature = [
            (atom.resname, atom.atom_name) for atom in definition.atoms
        ]
        for _ in range(count):
            expected_atoms.extend(atom_signature)

    coordinate_atoms = [
        (str(atom.resname), str(atom.name)) for atom in final_universe.atoms
    ]
    if len(coordinate_atoms) != len(expected_atoms):
        raise TopologyError(
            "Generated topology does not match the final coordinates: "
            f"topology describes {len(expected_atoms)} atoms, coordinates contain "
            f"{len(coordinate_atoms)}."
        )

    for index, (observed, expected) in enumerate(
        zip(coordinate_atoms, expected_atoms), 1
    ):
        if observed != expected:
            raise TopologyError(
                "Generated topology atom order does not match the final "
                f"coordinates at atom {index}: coordinates have "
                f"{observed[0]}:{observed[1]}, topology expects "
                f"{expected[0]}:{expected[1]}."
            )

    audit: Dict[str, object] = {
        "net_charge": net_charge,
        "coordinate_atoms": len(coordinate_atoms),
        "topology_atoms": len(expected_atoms),
        "molecule_charges": molecule_charges,
        "definition_files": {
            mtype: definition_paths[mtype]
            for mtype, _count in molecules
        },
    }
    topology_result.charge_audit = audit
    return audit
