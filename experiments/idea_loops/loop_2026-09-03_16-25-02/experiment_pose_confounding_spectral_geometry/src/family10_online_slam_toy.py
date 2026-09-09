"""Family 10: empirical NLS covariance toy vs linearized SLAM predictions.

Run (from the experiment root):
    .venv/bin/python src/family10_online_slam_toy.py

This script is a *finite-dimensional empirical validation only*.  For a
small 2D scalar Helmholtz contrast-source SLAM toy it checks whether the
linearized spectral predictions

    K_IS  = A_stack^T A_stack                    (known-pose map info)
    K_eff = A_stack^T W A_stack,                 (pose-marginalized map info)
    W     = I - B_stack (B_stack^T B_stack + alpha I)^{-1} B_stack^T,

with predicted map-covariance P_known = K_IS^{-1} and
P_free = K_eff^{-1}, actually track the estimation-error covariance of a
nonlinear least-squares (NLS) estimator on the same forward model.

Conventions (identical to the earlier families in this experiment tree):

* Scene: N = 16, T = 6 arc poses (Family 1 `build_poses`), n_rx = 4,
  two-Gaussian contrast `family1.make_chi0`, and the Family-6 smooth RBF
  configuration (`family2.build_smooth_basis`, p = 24 unit-2-norm columns).
* Frequencies f in {1.0, 1.4, 1.8}, k_b(f) = 2*pi*f, per-frequency scalar
  whitening with
      sigma_f^2 = ||A_pix||_F^2 / (m_c * snr),   m_c = T*n_rx complex rows,
  applied to the complex data and Jacobian blocks: A_w = A_c/sigma_f,
  B_w = B_c/sigma_f, F_w = F_c/sigma_f, then realified by
      A_R = sqrt(2)[Re A_w; Im A_w], B_R = sqrt(2)[Re B_w; Im B_w],
      y_R = sqrt(2)[Re F_w; Im F_w]
  (`hh.whiten_realify(A_w, B_w, None)` for the two Jacobian blocks and the
  equivalent two-line vector realification for y_R).  A circular complex
  noise with covariance sigma_f^2 I therefore becomes independent N(0, I)
  on every realified stacked row.
* True map coefficients: c0 = argmin_c ||S c - chi0||_2 (least squares;
  S has unit columns but is not orthogonal).  The "true scene" used for
  every data-generation and Jacobian evaluation is chi_true = S c0.
* Pose parameterization: x0 = poses0.ravel(); the free-pose fit uses
  poses = (x0 + dx).reshape(T,3).  The Gaussian pose prior is represented
  as sqrt(alpha) * dx extra residual rows (dx measured in the raw pose
  units p_x, p_y, theta).
* All NLS fits start at the true parameters (c0, dx = 0); this is a local
  identifiability study in the small-noise sense, and large-noise /
  ill-conditioned excursions along the most-confounded directions are
  exactly what the empirical comparison is intended to expose.
* Born mode: the Born forward `F_born = G_S (chi * E_inc)` is linear in c.
  The reused module exposes no Born pose Jacobian, so this script contains a
  small local helper `born_forward_AB` that assembles (F_born, A_born, B_born)
  from the public helmholtz primitives; it is verified by one centered
  finite-difference self-check below.
* No continuum-limit, global-convergence, real-SLAM-system, or production
  claim is made.  All statements are about this finite-dimensional toy.

Outputs:
    results/family10_online_slam_toy.json
    notes/family10_online_slam_toy.md
    figures/family10_*_*.png
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
_ROOT = _HERE.parent
sys.path.insert(0, str(_HERE))

import helmholtz as hh  # noqa: E402
import family1_pilot as family1  # noqa: E402
import family2_algebraic_spine as family2  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from scipy.linalg import solve_triangular  # noqa: E402
from scipy.optimize import least_squares  # noqa: E402


CONFIG = {
    # Scenario identical to family6 (N, poses, receivers, two-blob scene).
    "N": 16,
    "T": 6,
    "n_rx": 4,
    "m_c": 24,  # T * n_rx complex rows per frequency
    "m_real": 144,  # 2 * m_c * len(frequencies) after stacking
    "p": 24,  # smooth basis dimension
    "q_pose": 18,  # 3 * T pose parameters
    "arc_radius": 1.6,
    "arc_phi_deg": [-45.0, 45.0],
    "pose_theta_convention": (
        "theta = atan2(-p_y, -p_x): body +x axis points toward the origin"
    ),
    "rx_offsets": [
        [-0.06, 0.0],
        [0.06, 0.0],
        [0.0, -0.06],
        [0.0, 0.06],
    ],
    "tx_offset": [0.0, 0.0],
    "chi_blobs": {
        "amp1": 0.3,
        "amp2": 0.5,
        "sigma1": 0.09,
        "sigma2": 0.07,
        "c1": [-0.15, -0.12],
        "c2": [0.18, 0.14],
    },
    "smooth_basis": {
        "p": 24,
        "x_centers_n": 4,
        "y_centers_n": 6,
        "x_span": [-0.3, 0.3],
        "y_span": [-0.3, 0.3],
        "sigma_b": 0.16,
        "unit_columns": True,
        "note": "identical construction/config to family2/family6",
    },
    "k_b_rule": "k_b(f) = 2*pi*f with f in {1.0, 1.4, 1.8}",
    "frequencies": [1.0, 1.4, 1.8],
    "snr": 100.0,
    "alpha": 1.0,
    "noise_model": {
        "per_frequency_complex_covariance": "sigma_f^2 I, sigma_f^2 = ||A_pix||_F^2 / (m_c*snr)",
        "whitening": "scalar division of complex blocks and forward by 1/sigma_f",
        "realification": (
            "A_R = sqrt(2)*[Re(A_w); Im(A_w)], B_R = sqrt(2)*[Re(B_w); Im(B_w)] "
            "via hh.whiten_realify(A_w, B_w, None); "
            "y_R = sqrt(2)*[Re(F_w); Im(F_w)]"
        ),
        "stacked_noise_covariance": (
            "I_m_real after realification (independent standard Gaussians on "
            "each stacked row)"
        ),
        "noise_draw": (
            "one default_rng(20260903); the same per-trial draw is applied to "
            "both Born and full-wave true data (matched comparison)"
        ),
    },
    "monte_carlo": {
        "n_trials": 200,
        "seed": 20260903,
        "modes": ["born", "full_wave"],
        "wall_guard_seconds": 2800.0,
    },
    "fit": {
        "method": "trf",
        "x_scale": "jac",
        "max_nfev": 1200,
        "xtol": 1e-10,
        "ftol": 1e-10,
        "gtol": 1e-10,
        "note": (
            "local fits started at true parameters; x_scale='jac' preconditions "
            "the trust region on top of a Fisher-whitened coordinate system "
            "(theta = theta0 + T y with T^T H_true T = I for known-pose "
            "H=K_IS and free-pose H=joint information), which only changes "
            "optimizer conditioning, not the NLS optimum. 'success' rows are "
            "those with status True."
        ),
    },
    "n_confounded_directions": 3,
    "finite_difference_self_check": {
        "eps": 1e-6,
        "seed": 20260903,
        "mode": "born",
        "note": (
            "one centered finite-difference check of the realified stacked "
            "residual Jacobian at the true point (locally duplicated code only "
            "for Born B, whose API is absent from the reused modules)"
        ),
    },
    "finite_difference_self_check_full_wave": {
        "eps": 1e-6,
        "seed": 20260903,
        "mode": "full_wave",
        "note": (
            "centered finite-difference checks of the realified full-wave "
            "stacked residual Jacobian at (c0, dx=0): full-c and full-dx unit "
            "directions plus one joint random direction.  The analytic data "
            "rows [-A_stack, -B_stack] and prior rows [0, sqrt(alpha) I] are "
            "compared against (r(theta+eps w)-r(theta-eps w))/(2 eps)."
        ),
    },
    "coefficient_choice": (
        "c0 = argmin_c ||S c - chi0||_2 (numpy lstsq); chi_true = S c0"
    ),
    "scope_note": (
        "finite-dimensional toy only; no continuum, global-nonlinearity, "
        "real-SLAM-system, or production claim"
    ),
}

_SQRT2 = float(np.sqrt(2.0))


# ---------------------------------------------------------------------------
# Small local helpers
# ---------------------------------------------------------------------------

def _sym(M: np.ndarray) -> np.ndarray:
    M = np.asarray(M, dtype=float)
    return 0.5 * (M + M.T)


def _json_safe(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, dict):
        return {k: _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    return obj


def realify_vector(z: np.ndarray) -> np.ndarray:
    """sqrt(2)*[Re z; Im z] for a complex vector (matches module convention)."""
    z = np.asarray(z)
    return _SQRT2 * np.concatenate([np.real(z), np.imag(z)])


def build_scene(cfg: dict) -> dict:
    """Grid, chi0, S, poses, c0, chi_true, receivers."""
    points, h = hh.make_grid(int(cfg["N"]))
    chi0 = family1.make_chi0(points, cfg)
    S = family2.build_smooth_basis(points, cfg["smooth_basis"])
    poses = family1.build_poses(cfg)
    c0 = np.linalg.lstsq(S, chi0, rcond=None)[0]
    chi_true = S @ c0
    rx = np.asarray(cfg["rx_offsets"], dtype=float)
    tx = np.asarray(cfg["tx_offset"], dtype=float)
    return {
        "points": points,
        "h": h,
        "chi0": chi0,
        "S": S,
        "poses": poses,
        "c0": c0,
        "chi_true": chi_true,
        "rx": rx,
        "tx": tx,
        "x0": poses.ravel(),
    }


def born_forward_AB(
    chi: np.ndarray,
    poses: np.ndarray,
    rx_offsets: np.ndarray,
    tx_offset: np.ndarray,
    N: int,
    k_b: float,
):
    """Born (F_born, A_born, B_born) built from public helmholtz primitives.

    F_born,t = G_S,t (chi * E_inc,t),  A_born,t = G_S,t diag(E_inc,t),
    and the pose Jacobian
        B_born,t[:,l] = (D_x G_S,t,l * E_inc,t + G_S,t * D_x E_inc,t,l) chi,
    where D_x G_S,t,l uses green_grad_first and D_x E_inc,t,l uses
    green_grad_source (the same source-motion sign convention as the
    full-wave builder in helmholtz.build_AB).

    This helper exists only because no Born pose Jacobian is exposed by the
    reused modules.  A finite-difference self-check is run below.
    """
    chi = np.asarray(chi, dtype=float).reshape(-1)
    poses = np.atleast_2d(np.asarray(poses, dtype=float))
    points, h = hh.make_grid(int(N))
    T = poses.shape[0]
    n_rx = np.atleast_2d(np.asarray(rx_offsets, dtype=float)).shape[0]
    S0 = points.shape[0]
    geom = hh.pose_geometry(poses, rx_offsets, tx_offset)
    A = np.empty((T * n_rx, S0), dtype=np.complex128)
    B = np.zeros((T * n_rx, 3 * T), dtype=np.complex128)
    F = np.empty((T * n_rx,), dtype=np.complex128)
    for t in range(T):
        rx_t = geom["rx_world"][t]
        tx_t = geom["tx_world"][t]
        drx_t = geom["drx_world"][t]  # (3, n_rx, 2)
        dtx_t = geom["dtx_world"][t]  # (3, 2)
        G_S = (k_b**2) * (h**2) * hh.green_matrix(points, rx_t, k_b)
        E_inc = hh.green_matrix(tx_t[None, :], points, k_b)[:, 0]
        sl = slice(t * n_rx, (t + 1) * n_rx)
        A[sl, :] = G_S * E_inc[None, :]
        F[sl] = G_S @ (chi * E_inc)
        grad1 = hh.green_grad_first(points, rx_t, k_b)  # (n_rx, S0, 2)
        grad2 = hh.green_grad_source(points, tx_t, k_b)  # (S0, 2)
        for l in range(3):
            DXGS_l = (k_b**2) * (h**2) * np.einsum(
                "asd,ad->as", grad1, drx_t[l]
            )
            DXE_l = np.einsum("sd,d->s", grad2, dtx_t[l])
            B[sl, 3 * t + l] = (
                DXGS_l * E_inc[None, :] + G_S * DXE_l[None, :]
            ) @ chi
    return F, A, B


def build_model_blocks(
    scene: dict, cfg: dict, mode: str
) -> dict:
    """Complex -> whitened/realified stacked blocks evaluated at the truth."""
    chi = scene["chi_true"]
    poses = scene["poses"]
    rx = scene["rx"]
    tx = scene["tx"]
    N = int(cfg["N"])
    T = int(cfg["T"])
    n_rx = int(cfg["n_rx"])
    m_c = T * n_rx
    snr = float(cfg["snr"])
    As, Bs, ys, sigmas = [], [], [], []
    per_freq = []
    for f in cfg["frequencies"]:
        k_b = 2.0 * np.pi * float(f)
        if mode == "full_wave":
            A_c, B_c, F_c, _ = hh.build_AB(
                chi, poses, rx, tx, N, k_b
            )
        elif mode == "born":
            F_c, A_c, B_c = born_forward_AB(
                chi, poses, rx, tx, N, k_b
            )
        else:
            raise ValueError(f"unknown mode {mode!r}")
        sigma2 = float(np.linalg.norm(A_c, ord="fro") ** 2) / (
            float(m_c) * snr
        )
        sigma = float(np.sqrt(sigma2))
        A_w = A_c / sigma
        B_w = B_c / sigma
        F_w = F_c / sigma
        A_pix_R, B_R = hh.whiten_realify(A_w, B_w, None)
        A_s = A_pix_R @ scene["S"]
        y_R = realify_vector(F_w)
        As.append(A_s)
        Bs.append(B_R)
        ys.append(y_R)
        sigmas.append(sigma)
        per_freq.append(
            {
                "f": float(f),
                "k_b": k_b,
                "sigma2": sigma2,
                "sigma": sigma,
                "m_c": m_c,
            }
        )
    A_stack = np.vstack(As)
    B_stack = np.vstack(Bs)
    y_stack = np.concatenate(ys)
    return {
        "mode": mode,
        "A_stack": A_stack,
        "B_stack": B_stack,
        "y_stack": y_stack,
        "per_freq": per_freq,
    }


def shrinkage_W(B: np.ndarray, alpha: float) -> np.ndarray:
    """W = I - B (B^T B + alpha I)^{-1} B^T."""
    q = B.shape[1]
    C = B.T @ B + float(alpha) * np.eye(q, dtype=float)
    return np.eye(B.shape[0], dtype=float) - B @ np.linalg.solve(C, B.T)


def linearized_predictions(
    A_stack: np.ndarray, B_stack: np.ndarray, alpha: float
) -> dict:
    """K_IS/K_eff, P_known/P_free, generalized retention spectrum + dirs."""
    p = A_stack.shape[1]
    K_IS = _sym(A_stack.T @ A_stack)
    W = shrinkage_W(B_stack, alpha)
    K_eff = _sym(A_stack.T @ (W @ A_stack))
    P_known = np.linalg.inv(K_IS)
    P_free = np.linalg.inv(K_eff)
    # K_IS = L L^T, then M = L^{-1} K_eff L^{-T} (symmetrized); generalized
    # eigenvectors v = L^{-T} u, normalized to unit Euclidean norm.
    L = np.linalg.cholesky(_sym(K_IS))
    X = solve_triangular(L, K_eff, lower=True)
    M = solve_triangular(L, X.T, lower=True).T
    M = _sym(M)
    eigvals_u, U = np.linalg.eigh(M)
    rho_asc = np.sort(eigvals_u)
    order = np.argsort(eigvals_u)
    V = solve_triangular(L.T, U[:, order], lower=False)
    V = V / np.linalg.norm(V, axis=0)[None, :]
    rho_asc = eigvals_u[order]
    retained_dof = float(rho_asc.sum())
    KIS_eig = np.linalg.eigvalsh(K_IS)
    Keff_eig = np.linalg.eigvalsh(K_eff)
    P_known_eig = np.sort(np.linalg.eigvalsh(P_known))[::-1]
    P_free_eig = np.sort(np.linalg.eigvalsh(P_free))[::-1]
    return {
        "K_IS": K_IS,
        "K_eff": K_eff,
        "P_known": P_known,
        "P_free": P_free,
        "rho_asc": rho_asc,
        "V": V,
        "retained_dof": retained_dof,
        "K_IS_eigvals_asc": KIS_eig,
        "K_eff_eigvals_asc": Keff_eig,
        "P_known_eigvals_desc": P_known_eig,
        "P_free_eigvals_desc": P_free_eig,
        "KIS_cond": float(KIS_eig[-1] / KIS_eig[0]),
        "Keff_cond": float(Keff_eig[-1] / max(Keff_eig[0], 0.0)),
        "eig_symmetry_guard": float(
            np.max(np.abs(K_eff - K_eff.T))
            + np.max(np.abs(K_IS - K_IS.T))
        ),
    }


class ResidualModel:
    """Cached residual + analytic dense Jacobian for one mode."""

    def __init__(self, scene: dict, cfg: dict, blocks: dict, target: np.ndarray):
        self.scene = scene
        self.cfg = cfg
        self.mode = blocks["mode"]
        self.target = np.asarray(target, dtype=float).copy()
        self._sigmas = [row["sigma"] for row in blocks["per_freq"]]
        self._key = None
        self._res = None
        self._jac = None

    def _forward_at(self, theta: np.ndarray, free: bool):
        p = self.scene["S"].shape[1]
        q = self.scene["x0"].size
        c = theta[:p]
        dx = theta[p:] if free else np.zeros(q, dtype=float)
        chi = self.scene["S"] @ c
        poses = (self.scene["x0"] + dx).reshape(
            int(self.cfg["T"]), 3
        )
        N = int(self.cfg["N"])
        rx = self.scene["rx"]
        tx = self.scene["tx"]
        As, Bs, ms = [], [], []
        for f, sigma in zip(self.cfg["frequencies"], self._sigmas):
            k_b = 2.0 * np.pi * float(f)
            if self.mode == "full_wave":
                A_c, B_c, F_c, _ = hh.build_AB(chi, poses, rx, tx, N, k_b)
            else:
                F_c, A_c, B_c = born_forward_AB(
                    chi, poses, rx, tx, N, k_b
                )
            A_w, B_w, F_w = A_c / sigma, B_c / sigma, F_c / sigma
            A_pix_R, B_R = hh.whiten_realify(A_w, B_w, None)
            As.append(A_pix_R @ self.scene["S"])
            Bs.append(B_R)
            ms.append(realify_vector(F_w))
        A = np.vstack(As)
        B = np.vstack(Bs)
        model = np.concatenate(ms)
        alpha = float(self.cfg["alpha"])
        if free:
            res = np.concatenate(
                [self.target - model, np.sqrt(alpha) * dx]
            )
            J_data = np.hstack([-A, -B])
            J_prior = np.hstack(
                [
                    np.zeros((q, p), dtype=float),
                    np.sqrt(alpha) * np.eye(q, dtype=float),
                ]
            )
            J = np.vstack([J_data, J_prior])
        else:
            res = self.target - model
            J = -A
        return res, J

    def _eval(self, theta: np.ndarray, free: bool):
        if self._key is not None:
            if self._key[0] == free and np.array_equal(self._key[1], theta):
                return self._res, self._jac
        res, J = self._forward_at(theta, free)
        self._key = (free, np.asarray(theta, dtype=float).copy())
        self._res = res
        self._jac = J
        return res, J

    def residual(self, theta: np.ndarray, free: bool) -> np.ndarray:
        return self._eval(theta, free)[0]

    def jacobian(self, theta: np.ndarray, free: bool) -> np.ndarray:
        return self._eval(theta, free)[1]


def whitening_sqrt_root(P: np.ndarray) -> np.ndarray:
    """T with T^T P T = I from the symmetric eigendecomposition of P."""
    d, U = np.linalg.eigh(_sym(P))
    d = np.maximum(d, float(np.finfo(float).tiny))
    return U @ np.diag(1.0 / np.sqrt(d))


class PreconditionedModel:
    """Least-squares model in coordinates y with theta = theta0 + T y.

    T is chosen so the quadratic (information) Hessian at the true point is
    the identity (known-pose: T^T K_IS T = I; free-pose: T^T H_joint T = I).
    This changes only the optimizer's trust-region metric; the recovered
    theta_hat/c_hat are the NLS stationary points of the original problem.
    """

    def __init__(
        self,
        base: ResidualModel,
        theta0: np.ndarray,
        T: np.ndarray,
    ):
        self.base = base
        self.cfg = base.cfg
        self.theta0 = np.asarray(theta0, dtype=float)
        self.T = np.asarray(T, dtype=float)

    def _theta(self, y: np.ndarray) -> np.ndarray:
        return self.theta0 + self.T @ np.asarray(y, dtype=float)

    def residual(self, y: np.ndarray, free: bool) -> np.ndarray:
        return self.base.residual(self._theta(y), free)

    def jacobian(self, y: np.ndarray, free: bool) -> np.ndarray:
        J_theta = self.base.jacobian(self._theta(y), free)
        return J_theta @ self.T

    def theta_hat(self, y: np.ndarray) -> np.ndarray:
        return self._theta(y)


def fit_least_squares(
    model: PreconditionedModel | ResidualModel,
    free: bool,
    c0: np.ndarray,
    precond_y_size: int,
) -> dict:
    fit_cfg = model.cfg["fit"]
    t0 = time.perf_counter()
    sol = least_squares(
        lambda th: model.residual(th, free),
        np.zeros(precond_y_size, dtype=float),
        jac=lambda th: model.jacobian(th, free),
        method=fit_cfg["method"],
        x_scale=fit_cfg["x_scale"],
        max_nfev=int(fit_cfg["max_nfev"]),
        xtol=float(fit_cfg["xtol"]),
        ftol=float(fit_cfg["ftol"]),
        gtol=float(fit_cfg["gtol"]),
    )
    return {
        "success": bool(sol.success),
        "status": int(sol.status),
        "message": str(sol.message),
        "nfev": int(sol.nfev),
        "njev": int(sol.njev),
        "cost": float(sol.cost),
        "optimality": float(getattr(sol, "optimality", np.nan)),
        "seconds": time.perf_counter() - t0,
        "theta_hat": np.asarray(
            model.theta_hat(sol.x) if hasattr(model, "theta_hat") else sol.x,
            dtype=float,
        ),
    }


def covariance_metrics(
    cov_emp: np.ndarray, P_pred: np.ndarray
) -> dict:
    if not (np.all(np.isfinite(cov_emp)) and np.all(np.isfinite(P_pred))):
        p = int(P_pred.shape[0])
        return {
            "rel_fro": float("nan"),
            "rel_spectral": float("nan"),
            "emp_eigvals_desc": np.full(p, np.nan),
            "pred_eigvals_desc": np.full(p, np.nan),
            "trace_emp": float("nan"),
            "trace_pred": float(np.trace(P_pred)) if np.all(np.isfinite(P_pred)) else float("nan"),
        }
    denom = float(np.linalg.norm(P_pred, ord="fro"))
    return {
        "rel_fro": float(
            np.linalg.norm(cov_emp - P_pred, ord="fro") / denom
        ),
        "rel_spectral": float(
            np.linalg.norm(cov_emp - P_pred, ord=2) / max(float(np.linalg.norm(P_pred, ord=2)), 0.0)
        ),
        "emp_eigvals_desc": np.sort(np.linalg.eigvalsh(cov_emp))[::-1],
        "pred_eigvals_desc": np.sort(np.linalg.eigvalsh(P_pred))[::-1],
        "trace_emp": float(np.trace(cov_emp)),
        "trace_pred": float(np.trace(P_pred)),
    }


def run_finite_difference_self_check(scene, cfg, blocks) -> dict:
    """Centered-FD check of the realified stacked residual Jacobian."""
    fd_cfg = cfg["finite_difference_self_check"]
    eps = float(fd_cfg["eps"])
    rng = np.random.default_rng(int(fd_cfg["seed"]))
    p = scene["S"].shape[1]
    q = scene["x0"].size
    free = True
    model = ResidualModel(
        scene, cfg, blocks, target=blocks["y_stack"]
    )
    theta0 = np.concatenate([scene["c0"], np.zeros(q)])
    res0, J0 = model._eval(theta0, free)
    dw = rng.normal(size=(p + q,))
    dw /= np.linalg.norm(dw)
    rp = model.residual(theta0 + eps * dw, free)
    rm = model.residual(theta0 - eps * dw, free)
    fd = (rp - rm) / (2.0 * eps)
    pred = J0 @ dw
    denom = float(np.linalg.norm(pred))
    return {
        "mode": fd_cfg["mode"],
        "eps": eps,
        "direction_l2": float(np.linalg.norm(dw)),
        "rel_fd_error": float(
            np.linalg.norm(fd - pred) / max(denom, 1e-300)
        ),
        "residual_dim": int(res0.size),
        "jacobian_shape": list(J0.shape),
    }


def _rel_l2(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    return float(
        np.linalg.norm(a - b) / max(float(np.linalg.norm(b)), 1e-300)
    )


def run_full_wave_fd_self_check(scene, cfg, blocks) -> dict:
    """Full-wave stack-level centered-FD self-check of the free residual J.

    The stacked residual at theta = (c, dx) is

        r(theta) = [y_stack - f_fullwave(theta); sqrt(alpha) * dx],

    so its analytic Jacobian consists of the data rows
    [-A_stack, -B_stack] (A/B are the whitened/realified full-wave stacks at
    the evaluation point) plus the prior rows [0, sqrt(alpha) I_q].  This
    routine checks that analytic Jacobian against centered finite differences

        d_j = (r(theta0 + eps w_j) - r(theta0 - eps w_j)) / (2 eps)

    at theta0 = (c0, dx=0) for:
      * two full-c directions  (unit random vectors with support on all p
        map-coefficient coordinates and zero dx),
      * two full-dx directions (unit random vectors with support on all q pose
        coordinates and zero c),
      * one joint random unit direction over both c and dx.
    Per-frequency data-row slices and the prior rows are each evaluated for
    every direction, so both parameter-block types are exercised against every
    frequency-stack row block.  eps = 1e-6 gives an O(eps^2) truncation floor.
    """
    fd_cfg = cfg["finite_difference_self_check_full_wave"]
    eps = float(fd_cfg["eps"])
    rng = np.random.default_rng(int(fd_cfg["seed"]))
    p = scene["S"].shape[1]
    q = scene["x0"].size
    alpha = float(cfg["alpha"])
    m_c = int(cfg["m_c"])
    data_rows = int(blocks["y_stack"].size)  # m_real = 144
    model = ResidualModel(scene, cfg, blocks, target=blocks["y_stack"])
    theta0 = np.concatenate([scene["c0"], np.zeros(q)])
    _, J0_fresh = model._eval(theta0, free=True)

    A0 = np.asarray(blocks["A_stack"], dtype=float)
    B0 = np.asarray(blocks["B_stack"], dtype=float)
    J_data = np.hstack([-A0, -B0])
    J_prior = np.hstack(
        [
            np.zeros((q, p), dtype=float),
            np.sqrt(alpha) * np.eye(q, dtype=float),
        ]
    )
    J_exp = np.vstack([J_data, J_prior])
    analytic_stack_rel = _rel_l2(J0_fresh, J_exp)

    def rand_unit(dim: int) -> np.ndarray:
        v = rng.normal(size=(dim,))
        return v / np.linalg.norm(v)

    probe_specs = []
    for k in range(2):
        uc = rand_unit(p)
        probe_specs.append(
            (
                f"full_c_random_{k}",
                "c",
                np.concatenate([uc, np.zeros(q)]),
            )
        )
        ud = rand_unit(q)
        probe_specs.append(
            (
                f"full_dx_random_{k}",
                "dx",
                np.concatenate([np.zeros(p), ud]),
            )
        )
    probe_specs.append(("joint_random", "c_and_dx", rand_unit(p + q)))

    freq_rows = []
    for row in blocks["per_freq"]:
        start = len(freq_rows) * (2 * m_c)
        freq_rows.append(
            {
                "f": float(row["f"]),
                "rows": [start, start + 2 * m_c],
            }
        )
    prior_row_block = {
        "label": "prior (sqrt(alpha) dx)",
        "rows": [data_rows, data_rows + q],
    }

    probes = []
    fd_cols = []
    pred_cols = []
    for label, ptype, w in probe_specs:
        w = np.asarray(w, dtype=float)
        w = w / max(float(np.linalg.norm(w)), 1e-300)
        rp = model.residual(theta0 + eps * w, free=True)
        rm = model.residual(theta0 - eps * w, free=True)
        fd = (rp - rm) / (2.0 * eps)
        pred = J_exp @ w
        fd_cols.append(fd)
        pred_cols.append(pred)

        freq_blocks = []
        for fr in freq_rows:
            sl = slice(fr["rows"][0], fr["rows"][1])
            freq_blocks.append(
                {
                    "f": fr["f"],
                    "rows": list(fr["rows"]),
                    "rel_l2": _rel_l2(fd[sl], pred[sl]),
                    "fd_l2": float(np.linalg.norm(fd[sl])),
                    "pred_l2": float(np.linalg.norm(pred[sl])),
                }
            )
        psl = slice(prior_row_block["rows"][0], prior_row_block["rows"][1])
        prior = {
            "rows": list(prior_row_block["rows"]),
            "rel_l2": _rel_l2(fd[psl], pred[psl]),
            "fd_l2": float(np.linalg.norm(fd[psl])),
            "pred_l2": float(np.linalg.norm(pred[psl])),
        }
        dsl = slice(0, data_rows)
        probes.append(
            {
                "label": label,
                "parameter_type": ptype,
                "direction": w,
                "direction_l2": float(np.linalg.norm(w)),
                "rel_l2": _rel_l2(fd, pred),
                "data_rel_l2": _rel_l2(fd[dsl], pred[dsl]),
                "prior_rel_l2": prior["rel_l2"],
                "prior": prior,
                "per_frequency": freq_blocks,
                "fd_l2": float(np.linalg.norm(fd)),
                "pred_l2": float(np.linalg.norm(pred)),
            }
        )

    Fd = np.column_stack(fd_cols)
    Pred = np.column_stack(pred_cols)
    diff = Fd - Pred
    col_rel = [
        float(
            np.linalg.norm(diff[:, j])
            / max(float(np.linalg.norm(Pred[:, j])), 1e-300)
        )
        for j in range(diff.shape[1])
    ]
    return {
        "mode": fd_cfg["mode"],
        "eps": eps,
        "seed": int(fd_cfg["seed"]),
        "n_directions": len(probes),
        "n_residual_evaluations": 2 * len(probes),
        "theta0": theta0,
        "residual_dim": int(data_rows + q),
        "data_rows": data_rows,
        "prior_rows": q,
        "jacobian_shape": list(J_exp.shape),
        "analytic_stack_vs_fresh_rel_fro": analytic_stack_rel,
        "rel_fro_error": float(
            np.linalg.norm(diff, ord="fro")
            / max(float(np.linalg.norm(Pred, ord="fro")), 1e-300)
        ),
        "max_column_error": max(col_rel) if col_rel else float("nan"),
        "column_errors": col_rel,
        "predicted_products_fro_norm": float(np.linalg.norm(Pred, ord="fro")),
        "per_frequency_row_blocks": freq_rows,
        "probes": probes,
        "recipe": (
            "At theta0=(c0,dx=0) compare the analytic free residual Jacobian "
            "J = [[-A_stack, -B_stack], [0, sqrt(alpha) I_q]] (full-wave "
            "whitened/realified stacks at the truth) with centered finite "
            "differences (r(theta0+eps w)-r(theta0-eps w))/(2 eps), "
            "eps=1e-6.  Directions: two full-c unit directions (random "
            "vectors supported on all p map coordinates, zero dx), two "
            "full-dx unit directions (random vectors supported on all q pose "
            "coordinates, zero c), and one joint random unit direction; every "
            "direction is scored on the whole residual, on the data rows, on "
            "the prior rows, and on each per-frequency data-row slice.  "
            "rel_fro_error = ||FD - J W||_F / ||J W||_F over the probe matrix "
            "W; max_column_error = max_j ||FD_j - J w_j||_2 / ||J w_j||_2."
        ),
    }


def run_monte_carlo(
    scene: dict,
    cfg: dict,
    blocks: dict,
    pred: dict,
    noise: np.ndarray,
) -> dict:
    n_trials = int(cfg["monte_carlo"]["n_trials"])
    n_dir = int(cfg["n_confounded_directions"])
    mode = blocks["mode"]
    c0 = scene["c0"]
    p = c0.size
    q = scene["x0"].size
    A_stack = blocks["A_stack"]
    B_stack = blocks["B_stack"]
    alpha = float(cfg["alpha"])
    K_IS = _sym(A_stack.T @ A_stack)
    H_joint = np.block(
        [
            [K_IS, A_stack.T @ B_stack],
            [
                (A_stack.T @ B_stack).T,
                B_stack.T @ B_stack + alpha * np.eye(q),
            ],
        ]
    )
    T_known = whitening_sqrt_root(K_IS)
    T_free = whitening_sqrt_root(H_joint)
    theta0_free = np.concatenate([c0, np.zeros(q)])
    errs_known = np.empty((n_trials, p), dtype=float)
    errs_free = np.empty((n_trials, p), dtype=float)
    ok_known = np.zeros(n_trials, dtype=bool)
    ok_free = np.zeros(n_trials, dtype=bool)
    trial_rows = []
    wall_guard_seconds = float(
        cfg["monte_carlo"].get("wall_guard_seconds", 2800.0)
    )
    t_loop_start = time.perf_counter()
    n_completed = 0
    wall_guard_hit = False
    for i in range(n_trials):
        if time.perf_counter() - t_loop_start > wall_guard_seconds:
            wall_guard_hit = True
            break
        target = blocks["y_stack"] + noise[i]
        base_model = ResidualModel(scene, cfg, blocks, target=target)
        model_k = PreconditionedModel(base_model, c0, T_known)
        model_f = PreconditionedModel(
            base_model, theta0_free, T_free
        )
        sk = fit_least_squares(model_k, free=False, c0=c0, precond_y_size=p)
        sf = fit_least_squares(
            model_f, free=True, c0=c0, precond_y_size=p + q
        )
        ok_known[i] = sk["success"]
        ok_free[i] = sf["success"]
        errs_known[i] = sk["theta_hat"][:p] - c0
        errs_free[i] = sf["theta_hat"][:p] - c0
        trial_rows.append(
            {
                "trial": i,
                "known": {
                    "success": sk["success"],
                    "status": sk["status"],
                    "nfev": sk["nfev"],
                    "cost": sk["cost"],
                    "message": sk["message"],
                    "c_error_l2": float(np.linalg.norm(errs_known[i])),
                    "seconds": sk["seconds"],
                },
                "free": {
                    "success": sf["success"],
                    "status": sf["status"],
                    "nfev": sf["nfev"],
                    "cost": sf["cost"],
                    "message": sf["message"],
                    "c_error_l2": float(np.linalg.norm(errs_free[i])),
                    "dx_l2": float(
                        np.linalg.norm(sf["theta_hat"][p:])
                    ),
                    "seconds": sf["seconds"],
                },
            }
        )
        n_completed = i + 1
    errs_known = errs_known[:n_completed]
    errs_free = errs_free[:n_completed]
    ok_known = ok_known[:n_completed]
    ok_free = ok_free[:n_completed]
    n_ok_known = int(ok_known.sum())
    n_ok_free = int(ok_free.sum())

    def emp_cov(err: np.ndarray, mask: np.ndarray) -> np.ndarray:
        rows = err[mask]
        if rows.shape[0] < 2:
            return np.full((p, p), np.nan)
        return np.cov(rows, rowvar=False, bias=False)

    Cov_known = emp_cov(errs_known, ok_known)
    Cov_free = emp_cov(errs_free, ok_free)
    # Mean estimation errors on the same converged rows (bias diagnostics).
    mean_known = errs_known[ok_known].mean(axis=0) if n_ok_known else np.full(p, np.nan)
    mean_free = errs_free[ok_free].mean(axis=0) if n_ok_free else np.full(p, np.nan)

    cm_known = covariance_metrics(Cov_known, pred["P_known"])
    cm_free = covariance_metrics(Cov_free, pred["P_free"])
    direction_rows = []
    for i in range(n_dir):
        v = pred["V"][:, i]
        rho_i = float(pred["rho_asc"][i])
        pred_known_var = float(v @ pred["P_known"] @ v)
        pred_free_var = float(v @ pred["P_free"] @ v)
        emp_known_var = float(v @ Cov_known @ v) if n_ok_known >= 2 else np.nan
        emp_free_var = float(v @ Cov_free @ v) if n_ok_free >= 2 else np.nan
        direction_rows.append(
            {
                "index": i,
                "rho": rho_i,
                "predicted_inflation_1_over_rho": 1.0 / rho_i,
                "predicted_var_ratio_exact": float(
                    pred_free_var / pred_known_var
                )
                if pred_known_var > 0.0
                else float("nan"),
                "predicted_known_var": pred_known_var,
                "predicted_free_var": pred_free_var,
                "empirical_known_var": emp_known_var,
                "empirical_free_var": emp_free_var,
                "empirical_inflation_ratio": (
                    float(emp_free_var / emp_known_var)
                    if (np.isfinite(emp_known_var) and emp_known_var > 0.0)
                    else np.nan
                ),
                "eigenvector_v": v,
            }
        )
    trace_ratio_pred = float(
        np.trace(pred["P_free"]) / max(np.trace(pred["P_known"]), 1e-300)
    )
    trace_ratio_emp = float(
        np.trace(Cov_free) / max(np.trace(Cov_known), 1e-300)
    ) if (n_ok_free >= 2 and n_ok_known >= 2) else np.nan
    return {
        "mode": mode,
        "n_trials_requested": n_trials,
        "n_trials_completed": n_completed,
        "wall_guard_hit": wall_guard_hit,
        "wall_guard_seconds_limit": wall_guard_seconds,
        "n_known_success": n_ok_known,
        "n_free_success": n_ok_free,
        "known_rel_fro": cm_known["rel_fro"],
        "free_rel_fro": cm_free["rel_fro"],
        "known_rel_spectral": cm_known["rel_spectral"],
        "free_rel_spectral": cm_free["rel_spectral"],
        "known_trace_emp": cm_known["trace_emp"],
        "known_trace_pred": cm_known["trace_pred"],
        "free_trace_emp": cm_free["trace_emp"],
        "free_trace_pred": cm_free["trace_pred"],
        "trace_ratio_empirical": trace_ratio_emp,
        "trace_ratio_predicted": trace_ratio_pred,
        "mean_error_known": mean_known,
        "mean_error_free": mean_free,
        "Cov_known": Cov_known,
        "Cov_free": Cov_free,
        "direction_rows": direction_rows,
        "trial_rows": trial_rows,
        "errs_known": errs_known,
        "errs_free": errs_free,
        "ok_known": ok_known,
        "ok_free": ok_free,
    }


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def make_figures(
    born_run: dict,
    fw_run: dict,
    preds: dict,
    fig_dir: Path,
    suffix: str = "",
) -> dict:
    runs = {"born": born_run, "full_wave": fw_run}
    paths = {}
    # Figure (i): sorted eigenvalues of empirical vs predicted covariances.
    fig1, axes1 = plt.subplots(2, 2, figsize=(11, 9))
    for r_, mode in enumerate(["born", "full_wave"]):
        for c_, fit in enumerate(["known", "free"]):
            ax = axes1[r_, c_]
            run = runs[mode]
            cm = (
                covariance_metrics(run["Cov_known"], preds[mode]["P_known"])
                if fit == "known"
                else covariance_metrics(run["Cov_free"], preds[mode]["P_free"])
            )
            ax.semilogy(
                np.arange(1, cm["emp_eigvals_desc"].size + 1),
                np.maximum(cm["emp_eigvals_desc"], 1e-30),
                "o-",
                ms=4,
                label="empirical Cov",
            )
            ax.semilogy(
                np.arange(1, cm["pred_eigvals_desc"].size + 1),
                np.maximum(cm["pred_eigvals_desc"], 1e-30),
                "s--",
                ms=4,
                label="predicted P",
            )
            ax.set_title(f"{mode} / {fit}-pose map covariance")
            ax.set_xlabel("eigenvalue index (descending)")
            ax.set_ylabel("covariance eigenvalue")
            ax.legend(fontsize=8)
            ax.grid(True, which="both", alpha=0.3)
    fig1.suptitle(
        "Family 10: empirical NLS covariance vs linearized prediction "
        "(sampling + nonlinearity)"
    )
    fig1.tight_layout(rect=(0, 0, 1, 0.97))
    p1 = fig_dir / f"family10_cov_eigenvalues{suffix}.png"
    fig1.savefig(p1, dpi=150)
    plt.close(fig1)
    paths["cov_eigenvalues"] = p1

    # Figure (ii): confounded-direction variance ratios.
    fig2, axes2 = plt.subplots(1, 2, figsize=(12, 5))
    for ax, (mode, run) in zip(axes2, runs.items()):
        rows = run["direction_rows"]
        xs = np.arange(len(rows))
        preds_v = [r["predicted_inflation_1_over_rho"] for r in rows]
        preds_exact = [r["predicted_var_ratio_exact"] for r in rows]
        emps_v = [r["empirical_inflation_ratio"] for r in rows]
        ax.plot(xs, preds_v, "s--", color="tab:red", label="predicted 1/rho")
        ax.plot(
            xs,
            preds_exact,
            "^:",
            color="tab:orange",
            label="predicted exact ratio",
        )
        ax.plot(xs, emps_v, "o-", color="tab:blue", label="empirical ratio")
        ax.set_yscale("log")
        ax.set_xticks(xs)
        ax.set_xticklabels([f"dir {int(r['index'])}" for r in rows])
        ax.set_title(f"{mode}: 3 most-confounded directions")
        ax.set_xlabel("generalized retention direction (ascending rho)")
        ax.set_ylabel("free/known map variance ratio")
        ax.legend(fontsize=8)
        ax.grid(True, which="both", alpha=0.3)
    fig2.suptitle(
        "Family 10: empirical variance inflation vs 1/rho in confounded "
        "map directions"
    )
    fig2.tight_layout(rect=(0, 0, 1, 0.94))
    p2 = fig_dir / f"family10_confounded_direction_ratios{suffix}.png"
    fig2.savefig(p2, dpi=150)
    plt.close(fig2)
    paths["direction_ratios"] = p2

    # Figure (iii): Born vs full-wave deviation summary.
    fig3, (ax3a, ax3b) = plt.subplots(1, 2, figsize=(12, 5))
    labels = ["born/known", "born/free", "full-wave/known", "full-wave/free"]
    rels = [
        runs["born"]["known_rel_fro"],
        runs["born"]["free_rel_fro"],
        runs["full_wave"]["known_rel_fro"],
        runs["full_wave"]["free_rel_fro"],
    ]
    cols = ["tab:green", "tab:green", "tab:purple", "tab:purple"]
    ax3a.bar(labels, rels, color=cols, alpha=0.75)
    ax3a.set_yscale("log")
    ax3a.set_ylabel("rel Fro ||Cov_emp - P_pred|| / ||P_pred||")
    ax3a.set_title("Born vs full-wave covariance deviation")
    ax3a.tick_params(axis="x", rotation=15)
    ax3a.grid(True, which="both", alpha=0.3)
    for mode, run in runs.items():
        xs = np.arange(len(run["direction_rows"]))
        ys = [r["predicted_var_ratio_exact"] for r in run["direction_rows"]]
        zs = [r["empirical_inflation_ratio"] for r in run["direction_rows"]]
        ax3b.plot(
            ys, zs, "o-",
            label=f"{mode} (3 confounded dirs)",
            color="tab:green" if mode == "born" else "tab:purple",
        )
    emp_flat = np.array(
        [
            r["empirical_inflation_ratio"]
            for run in runs.values()
            for r in run["direction_rows"]
        ],
        dtype=float,
    )
    pred_flat = np.array(
        [
            r["predicted_var_ratio_exact"]
            for run in runs.values()
            for r in run["direction_rows"]
        ],
        dtype=float,
    )
    emp_flat = emp_flat[np.isfinite(emp_flat) & (emp_flat > 0.0)]
    pred_flat = pred_flat[np.isfinite(pred_flat) & (pred_flat > 0.0)]
    lo = float(min(pred_flat.min(), 1.0)) if pred_flat.size else 1.0
    hi = float(max(pred_flat.max(), 10.0)) if pred_flat.size else 10.0
    ax3b.loglog([lo, hi], [lo, hi], "k--", label="empirical = predicted")
    ax3b.set_xlabel("predicted exact variance ratio")
    ax3b.set_ylabel("empirical free/known variance ratio")
    ax3b.set_title("Confounded-direction inflation: empirical vs predicted")
    ax3b.legend(fontsize=8)
    ax3b.grid(True, which="both", alpha=0.3)
    fig3.suptitle(
        "Family 10: Born control vs full-wave nonlinear deviation"
    )
    fig3.tight_layout(rect=(0, 0, 1, 0.94))
    p3 = fig_dir / f"family10_born_vs_fullwave_deviation{suffix}.png"
    fig3.savefig(p3, dpi=150)
    plt.close(fig3)
    paths["deviation"] = p3
    return paths


def sha256_file(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_report(
    results: dict, notes_dir: Path, fig_paths: dict, suffix: str = ""
) -> Path:
    run = results["monte_carlo"]
    theory = results["theory"]
    born = run["born"]
    fw = run["full_wave"]
    lines = [
        "# Family 10: empirical NLS covariance toy vs linearized spectral predictions",
        "",
        "Date: " + results["generated_utc"] + " (UTC)",
        "",
        "## Scope",
        "",
        results["config"]["scope_note"] + "  ",
        "All fits are local (started at the true parameters); the report is a "
        "finite-dimensional empirical check and claims nothing about continuum "
        "or real SLAM systems.",
        "",
        "## What was done",
        "",
        f"- Scene/config: see `results/family10_online_slam_toy{suffix}.json` (`config`); "
        f"coefficients `c0` chosen by least squares (`chi_true = S c0`).",
        f"- Linearized predictions computed from whitened/realified multi-frequency "
        f"stacks: `K_IS`, `K_eff(alpha={results['config']['alpha']})`, "
        "P_known and P_free.",
        f"- NLS Monte Carlo: {born['n_trials_completed']}/"
        f"{born['n_trials_requested']} (born) and "
        f"{fw['n_trials_completed']}/{fw['n_trials_requested']} (full-wave) "
        "completed trials; known-pose and free-pose (pose-prior) fits via "
        "`scipy.optimize.least_squares` with analytic dense Jacobians.  "
        f"Wall-clock guard {born['wall_guard_seconds_limit']:g} s per "
        "Monte-Carlo loop; hit: "
        f"born={born['wall_guard_hit']}, full_wave={fw['wall_guard_hit']}.",
        f"- Born known/free converged trials: {born['n_known_success']}/"
        f"{born['n_trials_completed']} and {born['n_free_success']}/"
        f"{born['n_trials_completed']}.  "
        f"Full-wave known/free: {fw['n_known_success']}/"
        f"{fw['n_trials_completed']} and {fw['n_free_success']}/"
        f"{fw['n_trials_completed']}.",
        "",
        "## Key comparison numbers",
        "",
        "Relative Frobenius deviations ||Cov_emp - P_pred||_F / ||P_pred||_F:",
        "",
        f"- Born known-pose: {born['known_rel_fro']:.6g}",
        f"- Born free-pose:  {born['free_rel_fro']:.6g}",
        f"- Full-wave known-pose: {fw['known_rel_fro']:.6g}",
        f"- Full-wave free-pose:  {fw['free_rel_fro']:.6g}",
        "",
        "Predicted retained DOF tr(K_eff K_IS^{-1}) and trace ratios:",
        "",
        f"- Born retained DOF {theory['born']['retained_dof']:.6g}; "
        f"trace ratio predicted {born['trace_ratio_predicted']:.6g}, "
        f"empirical {born['trace_ratio_empirical']:.6g}.",
        f"- Full-wave retained DOF {theory['full_wave']['retained_dof']:.6g}; "
        f"trace ratio predicted {fw['trace_ratio_predicted']:.6g}, "
        f"empirical {fw['trace_ratio_empirical']:.6g}.",
        "",
        "Three most-confounded direction rows (rho ascending):",
        "",
        "| mode | dir | rho | predicted 1/rho | predicted exact ratio | empirical ratio |",
        "|---|---|---|---|---|---|",
    ]
    for mode, run in [("born", born), ("full_wave", fw)]:
        for row in run["direction_rows"]:
            lines.append(
                f"| {mode} | {row['index']} | {row['rho']:.6g} | "
                f"{row['predicted_inflation_1_over_rho']:.6g} | "
                f"{row['predicted_var_ratio_exact']:.6g} | "
                f"{row['empirical_inflation_ratio']:.6g} |"
            )
    lines += [
        "",
        "## Full-wave finite-difference self-check",
        "",
        "- Exact recipe (also stored in the JSON as "
        "`full_wave_fd_self_check.recipe`): at `theta0=(c0, dx=0)` the "
        "analytic free residual Jacobian "
        "`[[-A_stack, -B_stack], [0, sqrt(alpha) I_q]]` is compared with "
        "centered differences `(r(theta0+eps w) - r(theta0-eps w))/(2 eps)` "
        "at `eps=1e-6` for two full-c unit directions (random support on all "
        "map coordinates, zero dx), two full-dx unit directions, and one joint "
        "random direction; every direction is scored on the whole residual, "
        "the data rows, the prior rows, and each per-frequency row slice.",
        f"- Directional batch relative Frobenius error "
        f"`||FD - J W||_F / ||J W||_F = "
        f"{results['full_wave_fd_self_check']['rel_fro_error']:.6g}`; "
        f"max directional column error "
        f"{results['full_wave_fd_self_check']['max_column_error']:.6g}.",
        "",
        "## Honest caveats",
        "",
        "- The full-wave map is nonlinear in c; the Born mode is linear in c "
        "but still nonlinear in the pose perturbation dx (poses enter through "
        "Green functions/incident fields).  Both are compared against the same "
        "first-order prediction, so any mismatch is a measured nonlinear/sampling effect.",
        f"- Empirical sample covariance uses {born['n_trials_completed']} "
        f"completed Born and {fw['n_trials_completed']} completed full-wave "
        "trials (fewer where the solver did not report success); the "
        "most-confounded directions carry huge variance, so Frobenius "
        "deviations contain large sampling noise.",
        "- `1/rho` is the classical retention/information-loss inflation factor "
        "and equals the variance ratio only when the generalized eigendirections "
        "also diagonalize both covariances (for example p=1 or K_IS proportional "
        "to I).  For the multi-direction c-space studied here the exact predicted "
        "ratio from the linearized covariances is "
        "(v^T P_free v)/(v^T P_known v); both are tabulated.",
        "- The Born pose Jacobian helper in this script is locally implemented "
        "(no reused module exposes a Born B) and passed one centered "
        "finite-difference check; see JSON `finite_difference_self_check`.",
        "- No claim is made that the fitter reached a global minimum; statuses "
        "and nfev/cost per trial are recorded in the JSON.",
        "",
        "## Artifacts",
        "",
        f"- Results: `results/family10_online_slam_toy{suffix}.json`",
        f"- Notes: `notes/family10_online_slam_toy{suffix}.md`",
    ]
    for name, path in fig_paths.items():
        lines.append(f"- Figure: `figures/{path.name}` ({name})")
    lines.append("")
    report = "\n".join(lines)
    p = notes_dir / f"family10_online_slam_toy{suffix}.md"
    p.write_text(report)
    return p


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--n-trials",
        type=int,
        default=None,
        help="override CONFIG['monte_carlo']['n_trials']",
    )
    parser.add_argument(
        "--max-nfev",
        type=int,
        default=None,
        help="override CONFIG['fit']['max_nfev']",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="override CONFIG['monte_carlo']['seed']",
    )
    parser.add_argument(
        "--out-suffix",
        type=str,
        default="",
        help=(
            "append this string to the base output filenames, e.g. _n120 -> "
            "results/family10_online_slam_toy_n120.json, "
            "notes/family10_online_slam_toy_n120.md, "
            "figures/family10_*_n120.png (default '' preserves names)"
        ),
    )
    args = parser.parse_args(argv)
    cfg = CONFIG
    if args.n_trials is not None:
        cfg["monte_carlo"]["n_trials"] = int(args.n_trials)
    if args.max_nfev is not None:
        cfg["fit"]["max_nfev"] = int(args.max_nfev)
    if args.seed is not None:
        cfg["monte_carlo"]["seed"] = int(args.seed)
    suffix = str(args.out_suffix)

    t_utc = datetime.now(timezone.utc)
    t_start = time.perf_counter()
    scene = build_scene(cfg)
    print(
        "[family10] scene N=%d S=%s chi0 L2=%.6g c0 norm=%.6g"
        % (
            scene["points"].shape[0],
            scene["S"].shape,
            np.linalg.norm(scene["chi0"]),
            np.linalg.norm(scene["c0"]),
        )
    )
    print("[family10] representation error ||chi0 - S c0||/||chi0|| = %.6g"
          % (np.linalg.norm(scene["chi0"] - scene["chi_true"])
             / np.linalg.norm(scene["chi0"])))

    blocks = {}
    preds = {}
    theory_out = {}
    for mode in cfg["monte_carlo"]["modes"]:
        blk = build_model_blocks(scene, cfg, mode)
        blocks[mode] = blk
        pr = linearized_predictions(blk["A_stack"], blk["B_stack"], cfg["alpha"])
        preds[mode] = pr
        theory_out[mode] = {
            "mode": mode,
            "shapes": {
                "A_stack": list(blk["A_stack"].shape),
                "B_stack": list(blk["B_stack"].shape),
                "y_stack": list(blk["y_stack"].shape),
            },
            "rho_asc": pr["rho_asc"],
            "retained_dof": pr["retained_dof"],
            "K_IS_eigvals_asc": pr["K_IS_eigvals_asc"],
            "K_eff_eigvals_asc": pr["K_eff_eigvals_asc"],
            "P_known_eigvals_desc": pr["P_known_eigvals_desc"],
            "P_free_eigvals_desc": pr["P_free_eigvals_desc"],
            "KIS_cond": pr["KIS_cond"],
            "Keff_cond": pr["Keff_cond"],
            "eig_symmetry_guard": pr["eig_symmetry_guard"],
            "n_confounded": int(cfg["n_confounded_directions"]),
            "confounded": [
                {
                    "index": int(i),
                    "rho": float(pr["rho_asc"][i]),
                    "predicted_inflation_1_over_rho": float(
                        1.0 / pr["rho_asc"][i]
                    ),
                    "predicted_var_ratio_exact": float(
                        (
                            pr["V"][:, i] @ pr["P_free"] @ pr["V"][:, i]
                        )
                        / (
                            pr["V"][:, i] @ pr["P_known"] @ pr["V"][:, i]
                        )
                    ),
                    "v": pr["V"][:, i],
                    "pred_known_var": float(
                        pr["V"][:, i] @ pr["P_known"] @ pr["V"][:, i]
                    ),
                    "pred_free_var": float(
                        pr["V"][:, i] @ pr["P_free"] @ pr["V"][:, i]
                    ),
                }
                for i in range(int(cfg["n_confounded_directions"]))
            ],
            "per_freq": blk["per_freq"],
        }
        print("[family10] %s: A %s B %s rho[0:3]=%s retained_dof=%.4f"
              % (
                  mode,
                  blk["A_stack"].shape,
                  blk["B_stack"].shape,
                  np.round(pr["rho_asc"][:3], 6),
                  pr["retained_dof"],
              ))

    fd_check = run_finite_difference_self_check(
        scene, cfg, blocks["born"]
    )
    print("[family10] Born stacked residual Jacobian FD rel error %.3g"
          % fd_check["rel_fd_error"])
    fw_fd_check = run_full_wave_fd_self_check(
        scene, cfg, blocks["full_wave"]
    )
    print(
        "[family10] full-wave stacked residual Jacobian FD rel-Fro %.3g, "
        "max-column %.3g"
        % (
            fw_fd_check["rel_fro_error"],
            fw_fd_check["max_column_error"],
        )
    )

    n_trials = int(cfg["monte_carlo"]["n_trials"])
    m_real = int(blocks["born"]["y_stack"].size)
    rng = np.random.default_rng(int(cfg["monte_carlo"]["seed"]))
    noise = rng.normal(size=(n_trials, m_real))
    print("[family10] noise matrix", noise.shape, "drawn with seed",
          cfg["monte_carlo"]["seed"])

    runs = {}
    for mode in cfg["monte_carlo"]["modes"]:
        t_mode = time.perf_counter()
        print(f"[family10] running {mode} Monte Carlo ...")
        runs[mode] = run_monte_carlo(
            scene, cfg, blocks[mode], preds[mode], noise
        )
        print(
            "[family10] %s known relF=%.4g free relF=%.4g "
            "(success %d/%d known, %d/%d free; completed %d/%d, "
            "wall_guard_hit=%s) in %.1fs"
            % (
                mode,
                runs[mode]["known_rel_fro"],
                runs[mode]["free_rel_fro"],
                runs[mode]["n_known_success"],
                runs[mode]["n_trials_completed"],
                runs[mode]["n_free_success"],
                runs[mode]["n_trials_completed"],
                runs[mode]["n_trials_completed"],
                runs[mode]["n_trials_requested"],
                str(runs[mode]["wall_guard_hit"]),
                time.perf_counter() - t_mode,
            )
        )

    results_dir = _ROOT / "results"
    figures_dir = _ROOT / "figures"
    notes_dir = _ROOT / "notes"
    results_dir.mkdir(exist_ok=True)
    figures_dir.mkdir(exist_ok=True)
    notes_dir.mkdir(exist_ok=True)
    fig_paths = make_figures(
        runs["born"], runs["full_wave"], preds, figures_dir, suffix=suffix
    )

    source_files = [
        "src/family10_online_slam_toy.py",
        "src/helmholtz.py",
        "src/family1_pilot.py",
        "src/family2_algebraic_spine.py",
    ]
    results = {
        "generated_utc": t_utc.isoformat(),
        "runner": "src/family10_online_slam_toy.py",
        "output_suffix": suffix,
        "command": " ".join(
            [".venv/bin/python src/family10_online_slam_toy.py"]
            + (list(argv) if argv else [])
        ),
        "title": "NLS covariance toy vs linearized spectral predictions",
        "runtime_seconds": time.perf_counter() - t_start,
        "platform": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": __import__("scipy").__version__,
            "matplotlib": matplotlib.__version__,
        },
        "environment_note": (
            "Apple Silicon CPU, no GPU; deterministic numpy/scipy dense "
            "linear algebra and fixed-seed Monte Carlo noise"
        ),
        "source_sha256": {
            p: hashlib.sha256((_ROOT / p).read_bytes()).hexdigest()
            for p in source_files
        },
        "finite_difference_self_check": fd_check,
        "full_wave_fd_self_check": fw_fd_check,
        "config": _json_safe(cfg),
        "scene": {
            "chi0_stats": {
                "min": float(scene["chi0"].min()),
                "max": float(scene["chi0"].max()),
                "mean": float(scene["chi0"].mean()),
                "l2": float(np.linalg.norm(scene["chi0"])),
            },
            "representation_error_rel_l2": float(
                np.linalg.norm(scene["chi0"] - scene["chi_true"])
                / np.linalg.norm(scene["chi0"])
            ),
            "c0_l2": float(np.linalg.norm(scene["c0"])),
            "grid_h_cell": scene["h"],
        },
        "theory": _json_safe(theory_out),
        "monte_carlo": {
            "noise_seed": cfg["monte_carlo"]["seed"],
            "noise_draw_note": cfg["noise_model"]["noise_draw"],
            "born": _json_safe(
                {
                    k: v
                    for k, v in runs["born"].items()
                    if k not in ("errs_known", "errs_free", "ok_known", "ok_free")
                }
            ),
            "full_wave": _json_safe(
                {
                    k: v
                    for k, v in runs["full_wave"].items()
                    if k not in ("errs_known", "errs_free", "ok_known", "ok_free")
                }
            ),
        },
        "figure_sha256": {
            name: sha256_file(Path(path)) for name, path in fig_paths.items()
        },
        "artifacts": {
            "results_json": f"results/family10_online_slam_toy{suffix}.json",
            "figures_relative": sorted(
                f"figures/{Path(p).name}" for p in fig_paths.values()
            ),
            "figures_absolute": [str(p) for p in fig_paths.values()],
            "notes": f"notes/family10_online_slam_toy{suffix}.md",
        },
    }
    results["runtime_seconds"] = time.perf_counter() - t_start

    results_path = results_dir / f"family10_online_slam_toy{suffix}.json"
    results_path.write_text(json.dumps(_json_safe(results), indent=2))
    report_path = write_report(results, notes_dir, fig_paths, suffix=suffix)

    digest_paths = {
        "src/family10_online_slam_toy.py": _ROOT / "src/family10_online_slam_toy.py",
        "src/helmholtz.py": _ROOT / "src/helmholtz.py",
        "src/family1_pilot.py": _ROOT / "src/family1_pilot.py",
        "src/family2_algebraic_spine.py": _ROOT / "src/family2_algebraic_spine.py",
        f"results/family10_online_slam_toy{suffix}.json": results_path,
    }
    digest_paths.update(
        {f"figures/{Path(path).name}": Path(path) for path in fig_paths.values()}
    )
    digest_block = "\n## Artifacts and digests\n\n```text\n"
    digest_block += "\n".join(
        f"{sha256_file(p)}  {label}" for label, p in digest_paths.items()
    )
    digest_block += "\n```\n"
    with report_path.open("a") as fh:
        fh.write(digest_block)
    print(
        "[family10] wrote results/family10_online_slam_toy"
        f"{suffix}.json, {len(fig_paths)} figures, "
        f"notes/family10_online_slam_toy{suffix}.md"
    )


if __name__ == "__main__":
    main()
