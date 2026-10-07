"""MDAnalysis-independent reduction of atom contacts to molecule clashes."""

from __future__ import annotations

from typing import Dict, Iterable, Mapping, Sequence, Tuple


def closest_clash_per_residue(
    pairs: Iterable[Sequence[int]],
    distances: Sequence[float],
    environment_atom_classes: Sequence[str],
    environment_atom_resindices: Sequence[int],
    default_threshold: float,
    thresholds: Mapping[str, float],
) -> Tuple[Dict[int, float], Dict[int, Tuple[int, int]]]:
    """Reduce candidate atom pairs to the closest crossing per residue.

    ``pairs`` contains indices into the queried environment and inserted atom
    arrays. The class-specific cutoff is applied inclusively. The returned atom
    indices identify the exact closest contact for reporting.
    """
    minimum_distance: Dict[int, float] = {}
    minimum_pair: Dict[int, Tuple[int, int]] = {}
    for pair, distance in zip(pairs, distances):
        env_index = int(pair[0])
        inserted_index = int(pair[1])
        molclass = str(environment_atom_classes[env_index])
        cutoff = float(thresholds.get(molclass, default_threshold))
        distance = float(distance)
        if distance > cutoff:
            continue
        resindex = int(environment_atom_resindices[env_index])
        if resindex not in minimum_distance or distance < minimum_distance[resindex]:
            minimum_distance[resindex] = distance
            minimum_pair[resindex] = (env_index, inserted_index)
    return minimum_distance, minimum_pair
