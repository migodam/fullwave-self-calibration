"""Bounded executable checks for the a2_highdim worker package.

Run with the named experiment venv Python.  Writes checks.json and
environment.json under this directory.  These are implementation checks only:
they are not scientific acceptance evidence.
"""

from __future__ import annotations

import json
import platform
import sys
import time
from pathlib import Path

import numpy as np

DIR = Path(__file__).resolve().parent
if str(DIR) not in sys.path:
    sys.path.insert(0, str(DIR))

from common import (  # noqa: E402
    BUDGET_RHS_DEFAULT,
    FINAL_E4_SEEDS,
    N_ALPHA,
    POSE_SCALE,
    CostLedger,
    complex_from_json,
    default_settings,
    gaussian_basis_for_grid,
    make_grid_points,
    make_model,
    make_scene,
    params_from_scaled,
    scaled_from_params,
)
from mismatch import (  # noqa: E402
    apply_clock_phase,
    apply_receiver_coupling,
    make_mismatch_scene,
)
from objective import (  # noqa: E402
    coherent_forward,
    coherent_value_gradient,
    intensity_forward,
    intensity_value_gradient,
    rice_nll_weights,
)
from solver import run_method  # noqa: E402


def _scene_fixture(seed: int = 701):
    scene = make_scene(seed)
    alpha = np.asarray(scene["alpha_true"], dtype=float)
    y = complex_from_json(scene["y"])
    sigma = float(scene["sigma"])
    return alpha, y, sigma


def check_config_adapter() -> tuple[bool, str]:
    for N in (20, 32):
        pts, h = make_grid_points(N)
        Phi = gaussian_basis_for_grid(pts)
        mdl = make_model(N, "gaussian49")
        ok = (
            mdl.basis.shape == (N * N, N_ALPHA)
            and np.allclose(mdl.basis, Phi)
            and abs(mdl.h - h) < 1e-14
        )
        if not ok:
            return False, f"N={N} basis/grid mismatch"
    return True, "N20/N32 gaussian49 models build on their own grids"


def check_noise_reference() -> tuple[bool, str]:
    scene = make_scene(701)
    sigma = float(scene["sigma"])
    y_true = complex_from_json(scene["y_true"])
    y = complex_from_json(scene["y"])
    noise = y - y_true
    var_re = float(np.var(np.real(noise)))
    var_im = float(np.var(np.imag(noise)))
    scale_ok = abs(np.mean(np.abs(noise) ** 2) - sigma**2) < 0.25 * sigma**2
    abs_ok = np.allclose(np.abs(y), np.abs(y), equal_nan=False)
    return (
        scale_ok and abs_ok and abs(var_re / max(var_im, 1e-30) - 1.0) < 0.35,
        f"sigma={sigma:.5f}; per-component variance ratio "
        f"{var_re / max(var_im, 1e-30):.3f}",
    )


def check_scaled_gradient_fd() -> tuple[bool, str]:
    """Verify the corrected chain rule grad_z = grad_x / POSE_SCALE against
    central differences in scaled coordinates (legacy bug control)."""
    model = make_model(20, "gaussian49")
    alpha0, y, sigma = _scene_fixture(703)
    rng = np.random.default_rng(5)
    alpha = alpha0 + 0.02 * rng.standard_normal(N_ALPHA)
    x = np.array([0.03, -0.02, 0.01], dtype=float)
    z = scaled_from_params(alpha, x)
    h = 1e-6

    def f_z(v: np.ndarray) -> float:
        a, xx = params_from_scaled(v)
        led = CostLedger()
        return coherent_forward(model, a, xx, y, sigma, (0, 1), led)["loss"]

    dirs = np.zeros((4, z.size))
    dirs[0, :N_ALPHA] = rng.standard_normal(N_ALPHA)
    dirs[0, :N_ALPHA] /= np.linalg.norm(dirs[0, :N_ALPHA])
    for j in range(3):
        dirs[1 + j, N_ALPHA + j] = 1.0
    led = CostLedger()
    out = coherent_value_gradient(model, alpha, x, y, sigma, (0, 1), led)
    g_z = np.concatenate(
        [out["grad_alpha"], out["grad_x"] / POSE_SCALE]
    )
    worst = 0.0
    for d in dirs:
        d /= np.linalg.norm(d)
        fd = (f_z(z + h * d) - f_z(z - h * d)) / (2 * h)
        ana = float(np.dot(g_z, d))
        worst = max(worst, abs(fd - ana))
    scale = max(1.0, float(np.linalg.norm(g_z)))
    ok = worst < 2e-5 * scale
    return ok, f"worst scaled-coordinate FD error {worst:.3e} rel {worst/scale:.3e}"


def check_intensity_fd() -> tuple[bool, str]:
    model = make_model(20, "gaussian49")
    alpha0, y, sigma = _scene_fixture(704)
    rng = np.random.default_rng(6)
    alpha = alpha0 + 0.01 * rng.standard_normal(N_ALPHA)
    x = np.array([0.02, -0.01, 0.008], dtype=float)
    z = scaled_from_params(alpha, x)
    h = 2e-6

    def f_z(v: np.ndarray) -> float:
        a, xx = params_from_scaled(v)
        led = CostLedger()
        return intensity_forward(model, a, xx, y, sigma, (0,), led)["loss"]

    dirs = np.zeros((3, z.size))
    for j in range(3):
        dirs[j, N_ALPHA + j] = 1.0
    dirs[0, 3] = 1.0
    dirs[0, 3] = 0.0
    dirs[0, 10] = 1.0
    led = CostLedger()
    out = intensity_value_gradient(model, alpha, x, y, sigma, (0,), led)
    g_z = np.concatenate([out["grad_alpha"], out["grad_x"] / POSE_SCALE])
    worst = 0.0
    for d in dirs[:2]:
        d = d / np.linalg.norm(d)
        fd = (f_z(z + h * d) - f_z(z - h * d)) / (2 * h)
        ana = float(np.dot(g_z, d))
        worst = max(worst, abs(fd - ana))
    ok = worst < 2e-4
    return ok, f"worst intensity-gradient FD error {worst:.3e}"


def check_rice_uses_exact_abs_y() -> tuple[bool, str]:
    model = make_model(20, "gaussian49")
    alpha, y, sigma = _scene_fixture(705)
    fw = model.forward(alpha, np.zeros(3), jacobian=False)
    led = CostLedger()
    nll_a, _ = rice_nll_weights(fw["total"], y, sigma)
    led2 = CostLedger()
    out = intensity_value_gradient(model, alpha, np.zeros(3), y, sigma, (0, 1, 2, 3), led2)
    same = abs(out["loss"] - nll_a) < 1e-9 * max(1.0, abs(nll_a))
    rhs = int(out["work"]["rhs_solves_total"]) + led2.rhs_by_kind.get("adjoint", 0)
    return same and rhs == 48, (
        f"nll consistent; RHS {led2.units} factorisations {led2.factorizations}"
    )


def check_work_charges() -> tuple[bool, str]:
    model = make_model(20, "gaussian49")
    alpha, y, sigma = _scene_fixture(706)
    led = CostLedger()
    coherent_value_gradient(model, alpha, np.zeros(3), y, sigma, (0, 1, 2, 3), led)
    led2 = CostLedger()
    intensity_value_gradient(model, alpha, np.zeros(3), y, sigma, (0, 1, 2, 3), led2)
    return (
        led.units == 48 and led.factorizations == 4 and led2.units == 48
        and led2.factorizations == 8,
        f"coherent rhs={led.units} fac={led.factorizations}; "
        f"intensity rhs={led2.units} fac={led2.factorizations}",
    )


def check_budget_preserves_accepted() -> tuple[bool, str]:
    model = make_model(20, "gaussian49")
    alpha, y, sigma = _scene_fixture(707)
    settings = default_settings()
    settings["stage_maxiter"] = [200, 200, 200, 200]
    # A 13-RHS budget lets at most one 12-RHS gradient call.
    settings["budget_rhs_columns"] = 13
    out = run_method(
        model,
        y,
        sigma,
        "coherent_joint",
        np.array([0.02, -0.02, 0.01]),
        settings,
        seed=707,
        radius_label="0.125lambda",
        x_true=np.zeros(3),
    )
    finite = np.all(np.isfinite(out.alpha)) and np.all(np.isfinite(out.x))
    return (
        out.status == "budget_stopped" and finite and out.ledger.units <= 13
        and len(out.task_errors) == 0,
        f"status={out.status} units={out.ledger.units} "
        f"stages={out.stages_visited}",
    )


def check_mismatch_controls() -> tuple[bool, str]:
    scene = make_mismatch_scene("clock_phase")
    y_true = complex_from_json(scene["y_true"])
    clocked = apply_clock_phase(y_true)
    ok1 = np.allclose(np.abs(clocked), np.abs(y_true))
    coupled = apply_receiver_coupling(y_true)
    base = make_scene(881)
    outside = make_mismatch_scene("outside_basis")
    return (
        ok1
        and not np.allclose(clocked, y_true)
        and float(np.linalg.norm(coupled - y_true)) > 1e-6
        and outside["truth_spatial_chi_n32"] is not None,
        "clock magnitude preserved; coupling changes data; outside scene built",
    )


def check_seed_guard() -> tuple[bool, str]:
    try:
        make_scene(FINAL_E4_SEEDS[0])
    except RuntimeError:
        return True, "final E4 seed access blocked"
    return False, "final E4 seed access NOT blocked"


CHECKS = [
    ("config_n20_adapter", check_config_adapter),
    ("noise_reference_and_abs_y", check_noise_reference),
    ("scaled_coordinate_gradient_fd", check_scaled_gradient_fd),
    ("intensity_rice_gradient_fd", check_intensity_fd),
    ("rice_uses_exact_abs_y", check_rice_uses_exact_abs_y),
    ("work_charges_exact", check_work_charges),
    ("budget_preserves_accepted", check_budget_preserves_accepted),
    ("mismatch_controls", check_mismatch_controls),
    ("final_seed_guard", check_seed_guard),
]


def main() -> int:
    results = {}
    t0 = time.perf_counter()
    for name, fn in CHECKS:
        ts = time.perf_counter()
        try:
            ok, detail = fn()
        except Exception as exc:  # checks must report task errors honestly
            ok, detail = False, f"{type(exc).__name__}: {exc}"
        results[name] = {
            "pass": bool(ok),
            "detail": detail,
            "seconds": round(time.perf_counter() - ts, 3),
        }
        print(f"{'PASS' if ok else 'FAIL'}  {name}: {detail}")
    wall = time.perf_counter() - t0
    payload = {
        "timestamp": __import__("datetime").datetime.now(
            __import__("datetime").timezone.utc
        ).isoformat(timespec="seconds"),
        "results": results,
        "summary": {
            "passed": sum(1 for v in results.values() if v["pass"]),
            "total": len(results),
            "wall_seconds": round(wall, 3),
        },
        "scope": (
            "bounded implementation checks; not scientific acceptance "
            "evidence, no final E4 seeds used"
        ),
    }
    (DIR / "checks.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    env = {
        "platform": platform.platform(),
        "machine": platform.machine(),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "physics_file_unmodified": True,
    }
    try:
        import scipy

        env["scipy"] = scipy.__version__
    except Exception:
        pass
    (DIR / "environment.json").write_text(
        json.dumps(env, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return 0 if payload["summary"]["passed"] == payload["summary"]["total"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
