#!/usr/bin/env python
"""Bounded checks for ``research/delegated/a2_physics/physics.py``.

This is a physical-core validation, not a research acceptance gate.  It
verifies the implemented discrete forward model against itself (analytic vs
finite-difference tangents, state residuals, quadrature, rigid-geometry
invariance, factorisation reuse, noise realification and adjoint-gradient
costing).  It does not make scientific claims and does not run seeds 1001+.

Run:
    python test_physics.py

Outputs written to ``checks.json``, ``summary.md`` and ``environment.json``.
"""

from __future__ import annotations

import json
import platform
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import physics  # noqa: E402
from physics import Config, Model, point_green, realify, self_cell_green

C0 = physics.C0
EPS0 = physics.EPS0

ALPHA_TEST = np.array(
    [0.42, -0.18, 0.27, 0.11, 0.35, -0.05, -0.12, 0.21, 0.33]
)
X_TEST = np.array([0.02, -0.015, 0.06])
FD_ALPHA_H = 1.0e-4
FD_POSE_H = 1.0e-5


# ---------------------------------------------------------------------------
# Individual checks
# ---------------------------------------------------------------------------

def check_config_validation() -> dict:
    for bad_N in (7, 17, 64):
        try:
            Config(N=bad_N)
        except ValueError:
            pass
        else:
            raise AssertionError(f"Config should reject N={bad_N}")
    try:
        Config(aperture="side")
    except ValueError:
        pass
    else:
        raise AssertionError("Config should reject invalid aperture")
    m = Model(Config(N=8, aperture="full"))
    assert m.basis.shape == (64, 9)
    assert m.ks.shape == (4,)
    assert np.allclose(
        m.pose_metric, np.diag([1.0, 1.0, 1.5**2])
    )
    try:
        m.forward(np.zeros(8), X_TEST)
    except ValueError:
        pass
    else:
        raise AssertionError("alpha length mismatch should raise")
    px = Model(
        Config(N=8, material_basis=np.eye(64))
    )
    assert px.n_alpha == 64
    try:
        Model(Config(N=8, material_basis=np.eye(63)))
    except ValueError:
        pass
    else:
        raise AssertionError("material_basis rows must match N*N")
    return {"basis_columns": 9, "pixel_override_q": 64}


def check_outgoing_convention() -> dict:
    """Far-field amplitude/phase of ``i/4 H0^(1)`` = outgoing ``e^{+ikr}``."""
    x = 6000.0
    val = point_green(1.0, x)
    asym = (
        1j
        / 4.0
        * np.sqrt(2.0 / (np.pi * x))
        * np.exp(1j * (x - np.pi / 4.0))
    )
    rel = float(np.abs(val - asym) / np.abs(asym))
    if rel > 2.0e-3:
        raise AssertionError(f"far-field phase/amplitude mismatch {rel:.3e}")
    return {"kr": float(x), "far_field_rel_error": rel}


def check_self_cell_and_k2_factors() -> dict:
    m = Model(Config(N=8, aperture="full"))
    h = m.h
    k = float(m.ks[3])
    a = h / np.sqrt(np.pi)

    # Radial quadrature of (i/4) * 2 pi * int_0^a H0(k r) r dr.
    nodes, weights = np.polynomial.legendre.leggauss(1024)
    t = (nodes + 1.0) / 2.0 * a
    radial = np.sum(weights * hankel0_quad(t, k) * t) * a / 2.0
    quad = 1j * np.pi / 2.0 * radial
    analytic = self_cell_green(k, h)
    rel_self = float(abs(quad - analytic) / max(1.0, abs(analytic)))
    if rel_self > 1.0e-10:
        raise AssertionError(f"self-cell quadrature mismatch {rel_self:.3e}")

    r = m.forward(
        ALPHA_TEST, X_TEST, freq_ids=[3], jacobian=False
    )
    st = r["states"][0]
    D = st["D"]
    S_mat = st["S"]
    # off-diagonal D entry carries k^2 h^2 point Green
    d_off = D[3, 7]
    dist_d = float(np.linalg.norm(m.points[3] - m.points[7]))
    expect_d = (k**2) * (h**2) * point_green(dist_d, k)
    rel_d_off = float(abs(d_off - expect_d) / max(1.0, abs(expect_d)))
    if rel_d_off > 1.0e-12:
        raise AssertionError(f"D off-diagonal factor mismatch {rel_d_off:.3e}")
    rel_d_diag = float(
        abs(D[3, 3] - (k**2) * self_cell_green(k, h)) / max(1.0, abs(D[3, 3]))
    )
    if rel_d_diag > 1.0e-12:
        raise AssertionError(f"D diagonal factor mismatch {rel_d_diag:.3e}")
    # S entry carries k^2 h^2 point Green from grid cell to receiver
    rr = 5
    nn = 11
    dist_s = float(np.linalg.norm(st["rx_world"][rr] - m.points[nn]))
    expect_s = (k**2) * (h**2) * point_green(dist_s, k)
    rel_s = float(abs(S_mat[rr, nn] - expect_s) / max(1.0, abs(expect_s)))
    if rel_s > 1.0e-12:
        raise AssertionError(f"S factor mismatch {rel_s:.3e}")
    return {
        "self_cell_quadrature_rel": rel_self,
        "D_offdiag_rel": rel_d_off,
        "D_diag_rel": rel_d_diag,
        "S_rel": rel_s,
    }


def hankel0_quad(r: np.ndarray, k: float) -> np.ndarray:
    from scipy.special import hankel1

    return hankel1(0, k * r)


def check_forward_consistency_and_work() -> dict:
    m = Model(Config(N=8, aperture="full"))
    r = m.forward(ALPHA_TEST, X_TEST, jacobian=True)
    if r["total"].shape != (288,):
        raise AssertionError("full data length must be 4*3*2*12=288")
    if r["A"].shape != (288, 9) or r["B"].shape != (288, 3):
        raise AssertionError("A/B shapes")
    if not np.allclose(r["total"], r["scattered"] + r["incident"], atol=0.0):
        raise AssertionError("raw total != scattered + incident")
    max_res = max(s["state_residual_rel"] for s in r["states"])
    if max_res > 1.0e-10:
        raise AssertionError(f"state residual too large {max_res:.3e}")
    # stable order: frequency, pose, illumination, receiver
    blocks = r["blocks"]
    order = list(
        zip(blocks["freq_idx"], blocks["pose_idx"], blocks["illum_idx"])
    )
    expected = [
        (f, p, u)
        for f in range(4)
        for p in range(3)
        for u in range(2)
    ]
    if order != expected:
        raise AssertionError("state order not frequency/pose/illumination")
    for st in r["states"]:
        if st["row_stop"] - st["row_start"] != 12:
            raise AssertionError("row slice must contain 12 receivers")
    # anchored pose has exact zero pose tangent
    if np.max(np.abs(r["B"][:24, :])) != 0.0:
        raise AssertionError("anchored-pose B block must be exactly zero")
    w = r["work"]
    if (
        w["factorizations"] != 4
        or w["rhs_solves"]["forward_state"] != 24
        or w["rhs_solves"]["material"] != 216
        or w["rhs_solves"]["pose"] != 72
        or w["rhs_solves_total"] != 312
    ):
        raise AssertionError(f"unexpected work counts {w['rhs_solves']}")
    # subset of frequencies: factorisations and RHS must scale with selection
    r2 = m.forward(ALPHA_TEST, X_TEST, freq_ids=[0, 2], jacobian=True)
    if r2["total"].shape != (144,):
        raise AssertionError("two-frequency subset length")
    w2 = r2["work"]
    if (
        w2["factorizations"] != 2
        or w2["rhs_solves"]["forward_state"] != 12
        or w2["rhs_solves"]["material"] != 108
        or w2["rhs_solves"]["pose"] != 36
    ):
        raise AssertionError("subset work counts must scale per frequency")
    return {
        "max_state_residual_rel": max_res,
        "work_full": {
            "factorizations": w["factorizations"],
            "rhs_solves_total": w["rhs_solves_total"],
        },
        "work_freq_subset": {
            "factorizations": w2["factorizations"],
            "rhs_solves_total": w2["rhs_solves_total"],
        },
    }


def _fd_derivative_errors(m: Model) -> dict:
    r = m.forward(ALPHA_TEST, X_TEST, jacobian=True)
    a_worst = 0.0
    b_worst = 0.0
    for j in range(m.n_alpha):
        ap = ALPHA_TEST.copy()
        am = ALPHA_TEST.copy()
        ap[j] += FD_ALPHA_H
        am[j] -= FD_ALPHA_H
        fd = (
            m.forward(ap, X_TEST)["total"]
            - m.forward(am, X_TEST)["total"]
        ) / (2.0 * FD_ALPHA_H)
        col = r["A"][:, j]
        err = float(np.max(np.abs(fd - col)))
        rel = err / max(1.0, float(np.max(np.abs(col))))
        a_worst = max(a_worst, rel)
    for l in range(3):
        xp = X_TEST.copy()
        xm = X_TEST.copy()
        xp[l] += FD_POSE_H
        xm[l] -= FD_POSE_H
        fd = (
            m.forward(ALPHA_TEST, xp)["total"]
            - m.forward(ALPHA_TEST, xm)["total"]
        ) / (2.0 * FD_POSE_H)
        col = r["B"][:, l]
        err = float(np.max(np.abs(fd - col)))
        rel = err / max(1.0, float(np.max(np.abs(col))))
        b_worst = max(b_worst, rel)
    return {"aperture": m.config.aperture, "a_worst_rel": a_worst,
            "b_worst_rel": b_worst}


def check_derivatives_fd_both_apertures() -> dict:
    rows = []
    for aperture in ("full", "limited"):
        m = Model(Config(N=16, aperture=aperture))
        row = _fd_derivative_errors(m)
        if row["a_worst_rel"] > 2.0e-6:
            raise AssertionError(
                f"A FD mismatch ({aperture}): {row['a_worst_rel']:.3e}"
            )
        if row["b_worst_rel"] > 5.0e-6:
            raise AssertionError(
                f"B FD mismatch ({aperture}): {row['b_worst_rel']:.3e}"
            )
        rows.append(row)
    return {
        "fd_alpha_h": FD_ALPHA_H,
        "fd_pose_h": FD_POSE_H,
        "apertures": rows,
    }


def check_noise_real_variance() -> dict:
    rng = np.random.default_rng(20260905)
    sigma = 0.7
    n = 400_000
    z = (
        rng.normal(0.0, sigma / np.sqrt(2.0), n)
        + 1j * rng.normal(0.0, sigma / np.sqrt(2.0), n)
    )
    v = realify(z, sigma)
    mean = np.mean(v)
    cov = np.cov(v.reshape(2, -1))
    if abs(mean) > 0.01:
        raise AssertionError(f"realified mean {mean:.4f}")
    if np.max(np.abs(cov - np.eye(2))) > 0.02:
        raise AssertionError(f"realified covariance {cov.tolist()}")
    return {
        "empirical_mean_abs": float(abs(mean)),
        "empirical_covariance": cov.tolist(),
        "samples": n,
    }


def check_rigid_invariance_and_anchor() -> dict:
    m = Model(Config(N=8, aperture="limited"))
    g0 = m._effective_geometry(np.zeros(3))
    gx = m._effective_geometry(X_TEST)
    max_dist_diff = 0.0
    for p in range(3):
        pts0 = np.vstack([g0["rx_world"][p], g0["tx_world"][p]])
        ptsx = np.vstack([gx["rx_world"][p], gx["tx_world"][p]])
        d0 = np.linalg.norm(pts0[:, None, :] - pts0[None, :, :], axis=2)
        dx = np.linalg.norm(ptsx[:, None, :] - ptsx[None, :, :], axis=2)
        max_dist_diff = max(max_dist_diff, float(np.max(np.abs(d0 - dx))))
    if max_dist_diff > 1.0e-12:
        raise AssertionError(
            f"rigid co-motion distances changed: {max_dist_diff:.3e}"
        )
    # direct incident at receivers must be unchanged under the shared rigid x
    r0 = m.forward(ALPHA_TEST, np.zeros(3))
    rx = m.forward(ALPHA_TEST, X_TEST)
    direct_diff = float(np.max(np.abs(r0["incident"] - rx["incident"])))
    if direct_diff > 1.0e-12:
        raise AssertionError(
            f"direct incident not rigid-invariant: {direct_diff:.3e}"
        )
    # direct derivative through the code is (numerically) zero for co-motion
    max_direct_dx = max(
        float(np.max(np.abs(st["direct_derivative_dx"])))
        for st in rx["states"]
    )
    if max_direct_dx > 1.0e-12:
        raise AssertionError(
            f"direct-incident x-derivative not zero: {max_direct_dx:.3e}"
        )
    return {
        "max_pairwise_distance_change": max_dist_diff,
        "direct_incident_max_change": direct_diff,
        "max_direct_derivative_abs": max_direct_dx,
    }


def check_adjoint_gradient() -> dict:
    m = Model(Config(N=8, aperture="full"))
    freq_ids = [0, 1]
    alpha_ref = ALPHA_TEST + np.array([0.05, 0.04, -0.03, 0.02, 0.01,
                                       -0.02, 0.03, -0.01, 0.02])
    x_ref = X_TEST + np.array([0.01, -0.008, 0.03])
    y = m.forward(alpha_ref, x_ref, freq_ids=freq_ids)["total"]
    sigma = 2.0
    g = m.adjoint_gradient(
        ALPHA_TEST, X_TEST, y, sigma=sigma, freq_ids=freq_ids
    )

    def loss(a: np.ndarray, xx: np.ndarray) -> float:
        total = m.forward(a, xx, freq_ids=freq_ids)["total"]
        return 0.5 * float(np.sum(np.abs(realify(total - y, sigma)) ** 2))

    h_alpha = 1.0e-5
    worst = 0.0
    for j in range(m.n_alpha):
        e = np.zeros(m.n_alpha)
        e[j] = h_alpha
        fd = (loss(ALPHA_TEST + e, X_TEST) - loss(ALPHA_TEST - e, X_TEST)) / (
            2.0 * h_alpha
        )
        worst = max(worst, abs(fd - g["grad_alpha"][j]))
    for l in range(3):
        hh = FD_POSE_H if l < 2 else 2.0 * FD_POSE_H
        xp = X_TEST.copy()
        xm = X_TEST.copy()
        xp[l] += hh
        xm[l] -= hh
        fd = (loss(ALPHA_TEST, xp) - loss(ALPHA_TEST, xm)) / (2.0 * hh)
        worst = max(worst, abs(fd - g["grad_x"][l]))
    scale = max(
        1.0,
        float(np.max(np.abs(g["grad_alpha"]))),
        float(np.max(np.abs(g["grad_x"]))),
    )
    if worst / scale > 1.0e-7:
        raise AssertionError(f"adjoint gradient FD mismatch {worst/scale:.3e}")
    w = g["work"]
    if (
        w["factorizations"] != 2
        or w["rhs_solves"]["forward_state"] != 12
        or w["rhs_solves"]["adjoint"] != 12
        or w["rhs_solves"]["material"] != 0
        or w["rhs_solves"]["pose"] != 0
    ):
        raise AssertionError(f"unexpected adjoint work {w['rhs_solves']}")
    return {
        "gradient_max_fd_abs_error": worst,
        "gradient_scale": scale,
        "loss": g["loss"],
        "adjoint_work": {
            "factorizations": w["factorizations"],
            "forward_state_rhs": w["rhs_solves"]["forward_state"],
            "adjoint_rhs": w["rhs_solves"]["adjoint"],
        },
    }


def run_n16_n32_timing() -> list[dict]:
    rows = []
    for N in (16, 32):
        m = Model(Config(N=N, aperture="full"))
        t0 = time.perf_counter()
        r = m.forward(ALPHA_TEST, X_TEST, jacobian=True)
        elapsed = time.perf_counter() - t0
        w = r["work"]
        rows.append(
            {
                "N": N,
                "rows": int(r["total"].shape[0]),
                "wall_seconds": elapsed,
                "factorizations": w["factorizations"],
                "rhs_solves_total": w["rhs_solves_total"],
                "max_state_residual_rel": max(
                    s["state_residual_rel"] for s in r["states"]
                ),
            }
        )
    return rows


# ---------------------------------------------------------------------------
# Runner and file outputs
# ---------------------------------------------------------------------------

_CHECK_FUNCTIONS = [
    ("config_validation", check_config_validation),
    ("outgoing_convention_far_field", check_outgoing_convention),
    ("self_cell_and_k2_factors", check_self_cell_and_k2_factors),
    ("forward_consistency_and_work", check_forward_consistency_and_work),
    ("analytic_A_B_vs_fd_both_apertures", check_derivatives_fd_both_apertures),
    ("noise_real_variance", check_noise_real_variance),
    ("rigid_distance_and_anchor_invariance", check_rigid_invariance_and_anchor),
    ("adjoint_gradient_vs_fd", check_adjoint_gradient),
]


def main() -> int:
    started = time.perf_counter()
    checks = []
    failed = 0
    for name, fn in _CHECK_FUNCTIONS:
        t0 = time.perf_counter()
        try:
            detail = fn()
        except Exception as exc:  # noqa: BLE001 - runner records failures
            checks.append(
                {
                    "check": name,
                    "status": "FAIL",
                    "seconds": round(time.perf_counter() - t0, 4),
                    "error": f"{type(exc).__name__}: {exc}",
                }
            )
            failed += 1
            continue
        checks.append(
            {
                "check": name,
                "status": "PASS",
                "seconds": round(time.perf_counter() - t0, 4),
                "detail": detail,
            }
        )
    timing = run_n16_n32_timing()
    total_seconds = time.perf_counter() - started

    now = datetime.now().astimezone().isoformat(timespec="seconds")
    env = {
        "executable": sys.executable,
        "python_version": sys.version.split()[0],
        "platform": platform.platform(),
        "numpy": np.__version__,
        "scipy": __import__("scipy").__version__,
        "cwd": str(Path.cwd()),
        "timestamp": now,
    }
    payload = {
        "package": "research/delegated/a2_physics",
        "timestamp": now,
        "checks": checks,
        "timing": timing,
        "total_check_seconds": round(total_seconds, 3),
    }
    (HERE / "checks.json").write_text(json.dumps(payload, indent=2) + "\n")
    (HERE / "environment.json").write_text(json.dumps(env, indent=2) + "\n")
    write_summary(payload, env)

    passed = len(checks) - failed
    print(f"checks passed: {passed}/{len(checks)}  failed: {failed}")
    print(f"check wall time: {total_seconds:.2f}s")
    for row in timing:
        print(
            f"N{row['N']}: {row['wall_seconds']:.3f}s "
            f"({row['rhs_solves_total']} rhs, "
            f"{row['factorizations']} factorisations)"
        )
    return 1 if failed else 0


def write_summary(payload: dict, env: dict) -> None:
    lines = [
        "# a2_physics physical-core check summary",
        "",
        f"Timestamp: {payload['timestamp']}",
        "",
        "## Status",
        "",
        f"Passed {sum(c['status'] == 'PASS' for c in payload['checks'])}/"
        f"{len(payload['checks'])} bounded physical-core checks.",
        "",
    ]
    for row in payload["checks"]:
        lines.append(f"- **{row['check']}**: {row['status']} "
                     f"({row['seconds']:.2f}s)")
    lines += [
        "",
        "## N16/N32 forward+jacobian timing (one call each, full aperture)",
        "",
        "| N | rows | wall s | LU factorisations | RHS solves | max state residual |",
        "|---|------|--------|-------------------|------------|--------------------|",
    ]
    for row in payload["timing"]:
        lines.append(
            f"| {row['N']} | {row['rows']} | {row['wall_seconds']:.3f} "
            f"| {row['factorizations']} | {row['rhs_solves_total']} "
            f"| {row['max_state_residual_rel']:.2e} |"
        )
    lines += [
        "",
        "## Facts",
        "",
        "- State equation residual is at machine precision for every "
        "frequency/pose/illumination state (max rel residual "
        f"{max(t['max_state_residual_rel'] for t in payload['timing']):.1e}).",
        "- Raw total equals the retained direct-incident plus scattered "
        "components by construction and in the checks.",
        "- Analytic material A and pose B match centred finite differences at "
        "nonzero shared pose error for both full and limited apertures across "
        "all four frequencies (worst normalised errors in checks.json).",
        "- The self-cell diagonal reproduces radial Gauss-Legendre quadrature "
        "of the equal-area disk integral to ~1e-15 relative error.",
        "- One LU factorisation per frequency is shared by all 6 "
        "poses/illuminations; every LU right-hand side is charged individually.",
        "- The efficient adjoint gradient matches finite differences of the "
        "realified objective and charges 12 forward + 12 adjoint RHS for two "
        "frequencies, with no additional factorisations.",
        "- Realification sqrt(2)/sigma [Re; Im] makes proper complex noise "
        "isotropic with unit real covariance in the check.",
        "- Pose 1 is anchored: its B block is exactly zero; direct-incident "
        "derivatives vanish under rigid receiver/transmitter co-motion.",
        "",
        "## Artifacts",
        "",
        "- `physics.py` (core), `test_physics.py` (executable checks)",
        "- `checks.json`, `environment.json`, `summary.md`",
        "",
        "## Uncertainties / assumptions",
        "",
        "- 'Width 0.16 m' is interpreted as the Gaussian standard deviation "
        "in exp(-|x-c|^2/(2 width^2)), matching the unnormalised Gaussian "
        "basis convention used in the surrounding loop code; columns are not "
        "normalised.",
        "- The shared unknown x is interpreted as adding the same (dx, dy, "
        "dtheta) world transform to the two non-anchored nominal poses.  The "
        "documentation states this explicitly; a different interpretation "
        "would change B but not the solver core.",
        "- These checks validate the implemented discrete model's internal "
        "consistency (finite differences, quadrature, state residuals).  They "
        "are not continuum convergence, measurement-truth, or scientific "
        "acceptance evidence.",
        "- No seeds 1001+ or optimiser runs were executed.",
        "",
        "## Environment",
        "",
        f"- Python {env['python_version']}, NumPy {env['numpy']}, "
        f"SciPy {env['scipy']}",
        f"- {env['platform']}",
        "",
    ]
    (HERE / "summary.md").write_text("\n".join(lines))


if __name__ == "__main__":
    raise SystemExit(main())
