"""Deterministic exact selection of counterions by discrete formal charge."""

from __future__ import annotations

import random
from typing import Dict, List, Sequence, Tuple


IonCandidate = Tuple[int, str, float]


def select_counterions(
    candidates: Sequence[IonCandidate],
    excess_charge: float,
    tolerance: float,
    random_seed: int,
) -> Tuple[List[dict], float]:
    """Return a subset whose removable charge matches ``excess_charge``.

    Candidates must already have the sign needed to correct the excess. A
    bounded subset-sum search is used instead of greedy removal, so mixed
    monovalent and multivalent ions cannot overshoot or miss a valid exact
    combination. Randomization affects which equivalent ions are chosen but
    remains reproducible for a fixed seed.
    """
    if tolerance <= 0:
        raise ValueError("tolerance must be > 0")
    if abs(excess_charge) <= tolerance:
        return [], excess_charge

    needed_sign = 1.0 if excess_charge > 0 else -1.0
    filtered = [item for item in candidates if item[2] * needed_sign > 0]
    rng = random.Random(random_seed)
    rng.shuffle(filtered)

    quantum = max(float(tolerance), 1.0e-6)
    target_units = int(round(abs(excess_charge) / quantum))
    tolerance_units = max(1, int(round(tolerance / quantum)))

    # sum_units -> tuple of candidate indices. Keeping the first path to each
    # sum makes equivalent choices follow the seeded shuffle deterministically.
    reachable: Dict[int, Tuple[int, ...]] = {0: ()}
    for index, (_, _, charge) in enumerate(filtered):
        units = int(round(abs(charge) / quantum))
        if units <= 0 or units > target_units + tolerance_units:
            continue
        additions: Dict[int, Tuple[int, ...]] = {}
        for current, chosen in list(reachable.items()):
            new_sum = current + units
            if new_sum > target_units + tolerance_units:
                continue
            if new_sum not in reachable and new_sum not in additions:
                additions[new_sum] = chosen + (index,)
        reachable.update(additions)

    acceptable = [
        value for value in reachable
        if abs(value - target_units) <= tolerance_units
    ]
    if not acceptable:
        return [], excess_charge

    best = min(acceptable, key=lambda value: (abs(value - target_units), value))
    removed = [
        {
            "resindex": int(filtered[index][0]),
            "resname": str(filtered[index][1]),
            "charge": float(filtered[index][2]),
        }
        for index in reachable[best]
    ]
    remaining = excess_charge - sum(item["charge"] for item in removed)
    return removed, remaining
