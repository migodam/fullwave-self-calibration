"""Predeclared reference design metrics and A2 margins (no test data)."""
from __future__ import annotations

import numpy as np

from common import CostLedger, SQRT2, write_json


def realify(z: np.ndarray, sigma: float) -> np.ndarray:
    return (SQRT2 / sigma) * np.concatenate([np.real(z), np.imag(z)])


def lowpass_task_operator(model32, ell_res: float) -> np.ndarray:
    """9x9 low-pass task operator for the fixed Gaussian material family.

    Q maps alpha onto the resolution-limited cell map, then projects back
    onto the same basis.  Width is the declared resolution length.
    """
    pts = model32.points
    n = model32.n_cells
    h = model32.h
    Phi = model32.basis
    # Width is explicitly one standard deviation, not a relabelled 4-sigma
    # support. Periodic boundary and coefficient projection are declared.
    sigma = float(ell_res)
    dx = np.abs(pts[:, 0, None] - pts[None, :, 0])
    dy = np.abs(pts[:, 1, None] - pts[None, :, 1])
    # periodic wrap on the world-fixed square
    dx = np.minimum(dx, 1.0 - dx)
    dy = np.minimum(dy, 1.0 - dy)
    K = np.exp(-(dx**2 + dy**2) / (2.0 * sigma**2))
    K = K / K.sum(axis=1, keepdims=True)
    smoothed = K @ Phi  # cells x 9
    gram = Phi.T @ Phi
    proj = np.linalg.solve(gram + 1e-12 * np.eye(9), Phi.T @ smoothed)
    return proj


def reference_metrics(
    model16,
    model32,
    sigma: float,
    aperture: str,
    lambda_min: float = 0.5,
    lever: float = 1.5,
) -> dict:
    """Nominal reference covariance at alpha=0.5, x=0 (N16 design)."""
    ref_ledger = CostLedger()
    alpha_ref = np.full(model16.n_alpha, 0.5)
    x_ref = np.zeros(3)
    fw = model16.forward(alpha_ref, x_ref, jacobian=True)
    ref_ledger.charge_model_work(fw["work"])
    A = fw["A"]
    B = fw["B"]
    AB = np.hstack([A, B])
    JR = realify(AB, sigma)
    gram = JR.T @ JR
    s = np.linalg.svd(gram, compute_uv=False)
    rank = int(np.count_nonzero(s > max(1e-10 * s[0], 1e-12)))
    full_rank = bool(rank == 12 and float(s[-1]) > 1e-8)
    if full_rank:
        C = np.linalg.inv(gram)
    else:
        C = np.linalg.pinv(gram, rcond=1e-10)
    C_alpha = C[:9, :9]
    C_x = C[9:, 9:]
    # lever-arm canonical coordinates u=(t, R*theta)
    Tc = np.array([1.0, 1.0, 1.0 / lever])
    Jx = JR[:, 9:]
    Ju = Jx * Tc[None, :]
    # Unknown-map pose covariance: use the joint inverse, not the known-map
    # inverse J_u^T J_u. The latter understated the noninferiority scale.
    C_u = np.diag([1,1,lever]) @ C_x @ np.diag([1,1,lever])
    # grid spacing of the inverse grid, resolution length
    grid_h = float(model16.h)
    na = 1.0 if aperture == "full" else np.sin(np.pi / 4.0)
    ell_res = max(2.0 * grid_h, lambda_min / (2.0 * na))
    Q = lowpass_task_operator(model32, ell_res)
    task_cov = Q @ C_alpha @ Q.T
    d_x = float(min(lambda_min / 16.0, float(ell_res) / 4.0))
    dx_ni = float(min(lambda_min / 32.0, float(np.sqrt(np.trace(C_u) / 3.0))))
    d_chi = float(2.0 * np.sqrt(np.trace(task_cov) / 9.0))
    dchi_ni = float(d_chi / 2.0)
    return {
        "aperture": aperture,
        "lambda_min": lambda_min,
        "lever": lever,
        "sigma": sigma,
        "ell_res": float(ell_res),
        "NA": float(na),
        "gram_rank": rank,
        "full_rank": full_rank,
        "design_singular_values": [float(v) for v in s],
        "condition_number": float(s[0] / s[-1]) if full_rank else None,
        "reference_flag": "ok" if full_rank else "parent_review_required",
        "tr_C_u_per_3": float(np.trace(C_u) / 3.0),
        "d_x_pose_m": float(d_x),
        "Delta_x_noninf_m": float(dx_ni),
        "d_chi_task": float(d_chi),
        "Delta_chi_noninf": float(dchi_ni),
        "q_T": 9,
        "reference_ledger": ref_ledger.snapshot(),
        "Q_res": Q.tolist(),
    }


def predeclare_reference(
    model16_full,
    model16_limited,
    model32_full,
    model32_limited,
    sigma_by_aperture: dict[str, float],
    out_path,
) -> dict:
    out = {}
    for ap, m16, m32 in (
        ("full", model16_full, model32_full),
        ("limited", model16_limited, model32_limited),
    ):
        out[ap] = reference_metrics(
            m16, m32, sigma_by_aperture[ap], ap
        )
    write_json(out_path, out)
    return out
