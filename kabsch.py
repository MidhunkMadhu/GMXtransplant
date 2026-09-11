"""
Shared plain Kabsch/SVD rigid-body best-fit superposition.

Used wherever the pipeline needs "take these N matched points and rotate +
translate them onto these other N matched points, minimizing RMSD" without
any sequence alignment or outlier rejection: ligand_replace.py's `pairfit`
and `autofit` modes, and align.py's `mask_fit` alignment method (the cpptraj-`rms`-style
"give two masks, fit atom-by-atom, no outlier rejection" option).
"""

from __future__ import annotations

import numpy as np


def kabsch_fit(mobile_coords: np.ndarray, ref_coords: np.ndarray):
    """Best-fit rigid-body rotation+translation taking mobile_coords onto
    ref_coords (classic Kabsch/SVD superposition, proper-rotation-only via
    the determinant-sign correction).

    mobile_coords and ref_coords must be (N, 3) arrays with matched rows --
    mobile_coords[i] is assumed to correspond to ref_coords[i]. This
    function does no correspondence-finding of its own (no sequence
    alignment, no nearest-neighbor matching, no outlier rejection); the
    caller is responsible for the two arrays already being in 1:1 order,
    exactly like cpptraj's `rms <ref_mask> <mask>` with two explicit masks.

    Returns (R, mobile_centroid, ref_centroid, rmsd). To apply the same
    transform to any other point set X (e.g. every atom of a ligand/
    receptor, not just the fit subset):
        X_fitted = (R @ (X - mobile_centroid).T).T + ref_centroid
    """
    if mobile_coords.shape != ref_coords.shape:
        raise ValueError(
            f"kabsch_fit: mobile_coords {mobile_coords.shape} and ref_coords "
            f"{ref_coords.shape} must have the same shape (matched N x 3 arrays)."
        )
    if mobile_coords.ndim != 2 or mobile_coords.shape[1] != 3:
        raise ValueError("kabsch_fit: coordinates must have shape N x 3")
    if len(mobile_coords) < 3:
        raise ValueError(
            "kabsch_fit: at least three point pairs are required for a unique "
            "three-dimensional rigid transform"
        )
    for label, coords in (("mobile", mobile_coords), ("reference", ref_coords)):
        centered = coords - coords.mean(axis=0)
        if np.linalg.matrix_rank(centered, tol=1e-8) < 2:
            raise ValueError(
                f"kabsch_fit: {label} points are collinear; at least three "
                "non-collinear pairs are required"
            )
    mobile_centroid = mobile_coords.mean(axis=0)
    ref_centroid = ref_coords.mean(axis=0)
    mobile_c = mobile_coords - mobile_centroid
    ref_c = ref_coords - ref_centroid
    H = mobile_c.T @ ref_c
    U, _S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    fitted = (R @ mobile_c.T).T + ref_centroid
    rmsd = float(np.sqrt(np.mean(np.sum((fitted - ref_coords) ** 2, axis=1))))
    return R, mobile_centroid, ref_centroid, rmsd
