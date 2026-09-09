#!/usr/bin/env python
"""Independent finite-difference/adjoint checks for the fixed-chart ROM.

These checks validate the implemented algebra only (a reusable numerical
foundation).  They are not physical-model or scientific-acceptance evidence.

Run:
    experiments/.../.venv/bin/python research/delegated/a3_rom/test_rom.py
"""

from __future__ import annotations

import json
import platform
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np
from scipy.linalg import qr, solve_triangular

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from common import (  # noqa: E402
    A2,
    FREQ_IDS,
    POSE_SCALE,
    WorkLedger,
    config_json,
    make_model,
    params_from_scaled,
    scaled_from_params,
    select_rows,
    write_json,
)
from basis import sensing_chart, twofold_chart  # noqa: E402
from wave import (  # noqa: E402
    chi_and_Tmat,
    fixed_chart_evaluate,
    solve_normal_via_qr,
    state_geometry,
)

ALPHA0 = np.array(
    [0.42, -0.18, 0.27, 0.11, 0.35, -0.05, -0.12, 0.21, 0.33]
)
X0 = np.array([0.02, -0.015, 0.06])
FD_H_ALPHA = 2e-5
FD_H_POSE = 2e-5


def _chart_for_test(model):
    alpha = ALPHA0
    x = X0
    geom = model._effective_geometry(x)
    ledger = WorkLedger()
    fi = 0
    U = twofold_chart(model, alpha, geom, fi, 48, ledger)
    return {fi: U}, ledger


def _fd_deriv(
    f, n_alpha: int, n_pose: int, alpha: np.ndarray, x: np.ndarray
) -> np.ndarray:
    """Centered finite-difference Jacobian d f / d[alpha, x]."""
    p = n_alpha + n_pose
    cols = []
    for j in range(p):
        if j < n_alpha:
            ap = alpha.copy()
            am = alpha.copy()
            ap[j] += FD_H_ALPHA
            am[j] -= FD_H_ALPHA
            fp = f(ap, x)
            fm = f(am, x)
        else:
            xp = x.copy()
            xm = x.copy()
            xp[j - n_alpha] += FD_H_POSE
            xm[j - n_alpha] -= FD_H_POSE
            fp = f(alpha, xp)
            fm = f(alpha, xm)
        cols.append((fp - fm) / (2.0 * FD_H_POSE if j >= n_alpha else 2.0 * FD_H_ALPHA))
    return np.column_stack(cols)


def _state_fd_one(model, chart, alpha, x):
    fi = 0
    geom = model._effective_geometry(x)
    chi, Tmat = chi_and_Tmat(model, alpha, fi)
    g = state_geometry(model, geom, fi, 0, 0)
    U = chart[fi]
    D = model._domain_operator(fi)
    DU = D @ U
    C = U - chi[:, None] * DU
    Qred, Rtri = qr(C, mode="economic", check_finite=False)
    QH = Qred.conj().T
    b = chi * g["E_inc"]
    c = solve_triangular(Rtri, QH @ b, lower=False, check_finite=False)
    return U @ c


def check_proxy_derivative_has_residual_term() -> dict:
    model = make_model(n=8)
    chart, _ = _chart_for_test(model)
    alpha = ALPHA0
    x = X0
    out = fixed_chart_evaluate(
        model, alpha, x, chart, np.zeros(288), 1.0, (0,), True
    )
    rec = out["states"][0]
    C = rec["C"]
    DU = rec["DU"]
    Rtri = rec["Rtri"]
    c = rec["c"]
    chi, Tmat = chi_and_Tmat(model, alpha, 0)
    Etot = rec["E_total_approx"]
    state_res = rec["state_residual"]

    # coefficient Jacobians from the full formula and the omitted term.
    rhs_with = (
        C.conj().T @ (Tmat * Etot[:, None])
        - DU.conj().T @ (np.conj(Tmat) * state_res[:, None])
    )
    c_a_with = solve_normal_via_qr(Rtri, rhs_with, WorkLedger(), 9)
    c_a_omitted = solve_normal_via_qr(
        Rtri, C.conj().T @ (Tmat * Etot[:, None]), WorkLedger(), 9
    )

    fd_j = np.empty((model.n_cells, model.n_alpha), dtype=np.complex128)
    for j in range(model.n_alpha):
        ap = alpha.copy()
        am = alpha.copy()
        ap[j] += FD_H_ALPHA
        am[j] -= FD_H_ALPHA
        fd_j[:, j] = (
            _state_fd_one(model, chart, ap, x) - _state_fd_one(model, chart, am, x)
        ) / (2.0 * FD_H_ALPHA)
    err_with = float(np.linalg.norm((rec["U"] @ c_a_with) - fd_j) / np.linalg.norm(fd_j))
    err_omit = float(
        np.linalg.norm((rec["U"] @ c_a_omitted) - fd_j) / np.linalg.norm(fd_j)
    )
    if err_with > 5e-5:
        raise AssertionError(
            f"full fixed-chart derivative formula failed FD {err_with:.3e}"
        )
    return {
        "rank": int(rec["U"].shape[1]),
        "err_with_residual_term": err_with,
        "err_omitted_residual_term": err_omit,
        "residual_abs": float(np.linalg.norm(state_res)),
    }


def check_proxy_output_jacobian_fd() -> dict:
    model = make_model(n=8)
    chart, _ = _chart_for_test(model)

    def f_total(alpha, x):
        return fixed_chart_evaluate(
            model,
            np.asarray(alpha, dtype=float),
            np.asarray(x, dtype=float),
            chart,
            np.zeros(288),
            1.0,
            (0,),
            False,
        )["total"]

    out = fixed_chart_evaluate(
        model, ALPHA0, X0, chart, np.zeros(288), 1.0, (0,), True
    )
    Jc = np.hstack([out["A"], out["B"]])
    Jfd = _fd_deriv(f_total, 9, 3, ALPHA0, X0)
    err = float(np.linalg.norm(Jc - Jfd) / max(np.linalg.norm(Jfd), 1e-300))
    if err > 5e-5:
        raise AssertionError(f"fixed-chart output Jacobian FD mismatch {err:.3e}")
    return {"proxy_output_jacobian_rel_error": err}


def check_reduced_loss_gradient_fd() -> dict:
    model = make_model(n=8)
    chart, _ = _chart_for_test(model)
    y = np.linspace(0.2, 0.9, 288) + 1j * np.linspace(-0.4, 0.4, 288)

    def f_loss(z):
        alpha, x = params_from_scaled(z)
        return fixed_chart_evaluate(
            model, alpha, x, chart, y, 1.0, (0,), False
        )["loss"]

    z0 = scaled_from_params(ALPHA0, X0)
    out = fixed_chart_evaluate(
        model, ALPHA0, X0, chart, y, 1.0, (0,), True
    )
    J = np.hstack([out["A"], out["B"] / POSE_SCALE[None, :]])
    res = out["total"] - select_rows(y, (0,))
    from common import realify

    grad = (realify(J, 1.0).T @ realify(res, 1.0)).real
    grads_fd = []
    for j in range(12):
        h = 2e-5
        zp = z0.copy()
        zm = z0.copy()
        zp[j] += h
        zm[j] -= h
        grads_fd.append((f_loss(zp) - f_loss(zm)) / (2.0 * h))
    err = float(
        np.linalg.norm(grad - np.asarray(grads_fd))
        / max(np.linalg.norm(grads_fd), 1e-300)
    )
    if err > 5e-5:
        raise AssertionError(f"reduced adjoint gradient FD mismatch {err:.3e}")
    return {"reduced_gradient_rel_error": err, "grad_norm": float(np.linalg.norm(grad))}


def check_fixed_chart_rank_trend() -> dict:
    """Smoke trend: increasing sensing/twofold rank lowers offline error."""
    model = make_model(n=8)
    from audit import compare_to_exact

    results = []
    for method in ("sensing", "twofold"):
        for r in (4, 12, 32, 60, 64):
            geom = model._effective_geometry(X0)
            led = WorkLedger()
            if method == "sensing":
                U = sensing_chart(model, geom, 3, r, led)
            else:
                U = twofold_chart(model, ALPHA0, geom, 3, r, led)
            red = fixed_chart_evaluate(
                model, ALPHA0, X0, {3: U}, np.zeros(288), 1.0, (3,), True
            )
            err = compare_to_exact(model, ALPHA0, X0, red, (3,))
            results.append(
                {
                    "method": method,
                    "rank": int(U.shape[1]),
                    "state_rel": err["state_rel_max"],
                    "output_rel": err["output_rel"],
                    "map_jac_rel": err["map_jac_rel"],
                    "pose_jac_rel": err["pose_jac_rel"],
                }
            )
    state_ok = results[-1]["state_rel"] < 1e-8 and results[-1]["output_rel"] < 1e-8
    if not state_ok:
        raise AssertionError("highest-rank twofold smoke audit did not improve")
    return {"rows": results}


def main() -> None:
    checks = {
        "proxy_derivative_residual_term": check_proxy_derivative_has_residual_term(),
        "proxy_output_jacobian_fd": check_proxy_output_jacobian_fd(),
        "reduced_loss_gradient_fd": check_reduced_loss_gradient_fd(),
        "fixed_chart_rank_trend_smoke": check_fixed_chart_rank_trend(),
    }
    env = {
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "scipy": __import__("scipy").__version__,
        "platform": platform.platform(),
    }
    payload = {
        "timestamp": datetime.now().astimezone().isoformat(timespec="seconds"),
        "status": "PASS",
        "checks": checks,
        "environment": env,
        "config": config_json(n=8),
        "notes": [
            "grid8 smoke only; finite differences are algebra checks, not physical acceptance",
            "no final E4 seeds and no nonlinear runs were executed by this test",
        ],
    }
    write_json(HERE / "checks.json", payload)
    print("A3 ROM algebra checks: PASS")
    for k, v in checks.items():
        print(f"- {k}: {json.dumps(v, sort_keys=True)}")
    (HERE / "environment.json").write_text(
        json.dumps(env, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
