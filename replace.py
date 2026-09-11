"""Remove the target receptor block and insert the full aligned replacement block.

Both blocks are controlled only by their configured masks. Every residue selected
by ``replacement_structure.receptor_mask``, including any bound ligands or
cofactors, is inserted after the one rigid alignment transform. The untouched
membrane/solvent/ion environment and the target box vectors are retained.
"""

from __future__ import annotations

import warnings
from collections import Counter
from dataclasses import dataclass, field
from typing import Dict, List

import numpy as np

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import MDAnalysis as mda

from masks import resolve_mask
from align import AlignmentResult


class ReplaceError(RuntimeError):
    pass


@dataclass
class ReplacementResult:
    merged_universe: "mda.Universe"
    environment_ag_size: int
    inserted_ag_size: int
    n_original_receptor_atoms_removed: int
    n_replacement_atoms_inserted: int
    inserted_residue_counts: Dict[str, int] = field(default_factory=dict)
    notes: List[str] = field(default_factory=list)


def build_environment(original_universe: "mda.Universe", original_receptor_mask: str):
    """Return every target atom outside ``target_box.receptor_mask``."""
    removed_ag = resolve_mask(original_universe, original_receptor_mask)
    removed_indices = set(removed_ag.indices.tolist())
    all_indices = np.arange(len(original_universe.atoms))
    keep_mask = ~np.isin(all_indices, list(removed_indices))
    environment_ag = original_universe.atoms[keep_mask]
    return environment_ag, removed_ag


def assemble_replacement(
    original_universe: "mda.Universe",
    replacement_universe: "mda.Universe",
    original_receptor_mask: str,
    replacement_receptor_mask: str,
    alignment_result: AlignmentResult,
) -> ReplacementResult:
    environment_ag, removed_ag = build_environment(original_universe, original_receptor_mask)

    mobile_full_ag = resolve_mask(replacement_universe, replacement_receptor_mask)
    if len(mobile_full_ag) != alignment_result.new_positions.shape[0]:
        raise ReplaceError(
            "Alignment result atom count does not match replacement_structure.receptor_mask "
            f"({alignment_result.new_positions.shape[0]} vs {len(mobile_full_ag)}); "
            "was the alignment computed against a different mask?"
        )
    # Write the aligned coordinates into the replacement universe in place
    # so downstream selections/merges see the post-alignment geometry.
    mobile_full_ag.positions = alignment_result.new_positions

    inserted_ag = mobile_full_ag
    inserted_residue_counts = dict(
        sorted(Counter(str(name) for name in inserted_ag.residues.resnames).items())
    )

    merged = mda.Merge(environment_ag, inserted_ag)
    merged.dimensions = original_universe.dimensions.copy()

    return ReplacementResult(
        merged_universe=merged,
        environment_ag_size=len(environment_ag),
        inserted_ag_size=len(inserted_ag),
        n_original_receptor_atoms_removed=len(removed_ag),
        n_replacement_atoms_inserted=len(inserted_ag),
        inserted_residue_counts=inserted_residue_counts,
        notes=[
            "Inserted the complete block selected by "
            "replacement_structure.receptor_mask; no ligand-specific stripping was applied."
        ],
    )
