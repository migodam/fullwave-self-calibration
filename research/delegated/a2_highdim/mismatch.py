"""Controlled model-mismatch scenes for the a2_highdim comparison.

All mismatches are known ahead of the runs and are never estimated by an
augmented nuisance model in this worker, so they test misspecification of
the coherent geometry-only parameterisation rather than a geometry proof.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from common import (
    BG_ALPHA,
    KS,
    MISMATCH_SEED,
    N_ALPHA,
    N_DATA,
    SQRT2,
    _complex_to_json,
    complex_from_json,
    gaussian_basis_for_grid,
    make_grid_points,
    make_model,
    make_scene,
    nominal_sigma,
)

CLOCK_OFFSET_M = 0.04
COUPLING_STRENGTH = 0.05


def _rows_per_freq_block(y: np.ndarray) -> list[np.ndarray]:
    y = np.asarray(y, dtype=np.complex128).reshape(-1)
    n = y.size // 4
    return [y[i * n : (i + 1) * n] for i in range(4)]


def apply_clock_phase(y: np.ndarray) -> np.ndarray:
    blocks = _rows_per_freq_block(y)
    out = [
        block * np.exp(1j * float(KS[i]) * CLOCK_OFFSET_M)
        for i, block in enumerate(blocks)
    ]
    return np.concatenate(out)


def _receiver_coupling_matrix() -> np.ndarray:
    R = 12
    C = np.eye(R, dtype=complex)
    P = np.zeros((R, R), dtype=complex)
    for i in range(R):
        P[i, (i + 1) % R] = 1.0
        P[i, (i - 1) % R] = 1.0
    C = C + COUPLING_STRENGTH * P
    return C


def apply_receiver_coupling(y: np.ndarray) -> np.ndarray:
    C = _receiver_coupling_matrix()
    blocks = _rows_per_freq_block(y)
    out = []
    for block in blocks:
        n_sub = block.size // 12
        parts = [
            C @ block[i * 12 : (i + 1) * 12] for i in range(n_sub)
        ]
        out.append(np.concatenate(parts))
    return np.concatenate(out)


def _outside_basis_truth_map() -> np.ndarray:
    """N32 pixel-wise positive contrast: the declared 49-basis background
    plus a narrow off-grid Gaussian inclusion that the 49 smooth basis
    cannot represent."""
    points, _ = make_grid_points(N_DATA)
    Phi = gaussian_basis_for_grid(points)

    background = Phi @ np.full(N_ALPHA, BG_ALPHA)
    centre = np.array([0.185, -0.285], dtype=float)
    width = 0.025
    d = points - centre[None, :]
    inclusion = 0.75 * np.exp(
        -np.einsum("ij,ij->i", d, d) / (2.0 * width**2)
    )
    return background + inclusion
def make_mismatch_scene(kind: str) -> dict[str, Any]:
    """Deterministic mismatch scene at seed 881 with 30 dB parent noise
    derived from the same nominal N32 reference as the in-basis scenes."""
    if kind not in ("clock_phase", "receiver_coupling", "outside_basis"):
        raise ValueError(f"unknown mismatch kind {kind!r}")
    sigma = nominal_sigma()
    base = make_scene(MISMATCH_SEED)
    alpha_true = np.asarray(base["alpha_true"], dtype=float)
    y_true_base = complex_from_json(base["y_true"])
    x_true = np.zeros(3)

    if kind in ("clock_phase", "receiver_coupling"):
        if kind == "clock_phase":
            y_true = apply_clock_phase(y_true_base)
            corruption = "shared_frequency_dependent_clock_phase_exp(i*k*0.04m)"
        else:
            y_true = apply_receiver_coupling(y_true_base)
            corruption = "nearest_neighbour_receiver_coupling_0.05"
        truth_spatial = None
        alpha_out = alpha_true
    else:
        model_pixel = make_model(N_DATA, "pixel32")
        chi_map = _outside_basis_truth_map()
        fw = model_pixel.forward(chi_map, x_true, jacobian=False)
        y_true = fw["total"]
        corruption = "outside_49_gaussian_basis_narrow_offgrid_inclusion"
        truth_spatial = chi_map
        alpha_out = None

    rng = np.random.default_rng(MISMATCH_SEED)
    noise = (sigma / SQRT2) * (
        rng.standard_normal(y_true.shape)
        + 1j * rng.standard_normal(y_true.shape)
    )
    rng_val = np.random.default_rng(10_000 + MISMATCH_SEED)
    noise_val = (sigma / SQRT2) * (
        rng_val.standard_normal(y_true.shape)
        + 1j * rng_val.standard_normal(y_true.shape)
    )
    return {
        "seed": MISMATCH_SEED,
        "mismatch_kind": kind,
        "corruption": corruption,
        "alpha_true": None if alpha_out is None else alpha_out.tolist(),
        "truth_spatial_chi_n32": (
            None if truth_spatial is None else truth_spatial.tolist()
        ),
        "x_true": x_true.tolist(),
        "y_true": _complex_to_json(y_true),
        "y": _complex_to_json(y_true + noise),
        "y_val": _complex_to_json(y_true + noise_val),
        "sigma": float(sigma),
    }


def mismatch_scene_arrays(scene: dict[str, Any]) -> dict[str, Any]:
    out = {
        "kind": str(scene["mismatch_kind"]),
        "alpha_true": (
            None
            if scene["alpha_true"] is None
            else np.asarray(scene["alpha_true"], dtype=float)
        ),
        "x_true": np.asarray(scene["x_true"], dtype=float),
        "y_true": np.asarray(scene["y_true"], dtype=complex),
        "y": np.asarray(scene["y"], dtype=complex),
        "y_val": np.asarray(scene["y_val"], dtype=complex),
        "sigma": float(scene["sigma"]),
        "truth_spatial": (
            None
            if scene["truth_spatial_chi_n32"] is None
            else np.asarray(scene["truth_spatial_chi_n32"], dtype=float)
        ),
    }
    return out
