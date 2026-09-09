"""Tests for the tangent-enabled FFT VIE adapter.

Checks are bounded development checks of the *discrete* model: they compare
the matrix-free implementation with the read-only dense parent and with
finite differences of the same discrete physical-current equation.  They are
not continuum-accuracy or production-acceptance gates.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from tangent_fft import TangentFFTVIE  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "a3_maxwell_fft"))
from maxwell_fft import DipoleVIE  # noqa: E402

OUT = Path(__file__).resolve().parent
PARENT_DIR = OUT.parents[2] / "research" / "trispace_self_calibration" / "a3_research"
sys.path.insert(0, str(PARENT_DIR))
import calibrate3d as parent_cal
from calibration_fft import FFTJointModel

SINGLE = ([[0.0, 0.0, 0.0]], [0.12], [2.4 + 0.03j])
PAIR = (
    [[-0.2, 0.0, 0.0], [0.2, 0.03, 0.04]],
    [0.12, 0.10],
    [2.4 + 0.03j, 3.0 + 0.04j],
)


def dense_parent_currents_checks() -> list[dict]:
    """FFT currents and material tangent vs dense parent currents."""
    rows: list[dict] = []
    cases = [
        ("single_k6", SINGLE, 0.06, 6.3),
        ("pair_k9", PAIR, 0.04, 9.0),
    ]
    for name, (centers, radii, eps), spacing, k in cases:
        parent = DipoleVIE(centers, radii, spacing, k, fill_quadrature=6)
        fft = TangentFFTVIE(
            centers,
            radii,
            spacing,
            k,
            fill_quadrature=6,
            gmres_rtol=1e-11,
            restart=120,
            maxiter=3000,
        )
        p_dense, dp_dense = parent.currents(eps, derivatives=True)
        p_fft, dp_fft = fft.currents(eps, derivatives=True)
        rel_p = float(
            np.linalg.norm(p_fft - p_dense) / np.linalg.norm(p_dense)
        )
        rel_dp = float(
            np.linalg.norm(dp_fft - dp_dense) / np.linalg.norm(dp_dense)
        )
        max_abs_dp = float(np.max(np.abs(dp_fft - dp_dense)))
        work = dict(fft.work)
        expected_rhs_columns = 4 * (1 + len(eps))
        expected_solves = 4 * (1 + len(eps))
        row = {
            "name": name + "_dense_parent_currents_and_tangent",
            "voxels": fft.kernel.n_vox,
            "dof": fft.ndof,
            "p_shape": list(p_fft.shape),
            "dp_shape": list(dp_fft.shape),
            "relative_current_error": rel_p,
            "relative_tangent_error": rel_dp,
            "max_abs_tangent_error": max_abs_dp,
            "true_relative_residuals": list(fft.true_relative_residuals),
            "gmres_converged": bool(
                fft.gmres_flags and all(f == 0 for f in fft.gmres_flags)
            ),
            "work": work,
            "fft_work": dict(fft.fft_work),
            "passed": bool(
                rel_p < 1e-7
                and rel_dp < 1e-7
                and work["factorizations"] == 0
                and work["rhs_columns"] == expected_rhs_columns
                and fft.gmres_flags
                and all(f == 0 for f in fft.gmres_flags)
                and len(fft.gmres_flags) == expected_solves
            ),
        }
        rows.append(row)
    return rows


def finite_difference_checks() -> list[dict]:
    """Material tangent vs central differences of the FFT currents."""
    rows: list[dict] = []
    cases = [
        ("single", SINGLE, 0.06, 6.3),
        ("pair", PAIR, 0.04, 9.0),
    ]
    delta = 1e-5
    for name, (centers, radii, eps), spacing, k in cases:
        fft = TangentFFTVIE(
            centers,
            radii,
            spacing,
            k,
            fill_quadrature=6,
            gmres_rtol=1e-11,
            restart=120,
            maxiter=3000,
        )
        _, dp = fft.currents(eps, derivatives=True)
        for region in range(len(eps)):
            plus = [complex(v) for v in eps]
            minus = [complex(v) for v in eps]
            plus[region] = complex(eps[region].real + delta, eps[region].imag)
            minus[region] = complex(eps[region].real - delta, eps[region].imag)
            p_plus = TangentFFTVIE(
                centers,
                radii,
                spacing,
                k,
                fill_quadrature=6,
                gmres_rtol=1e-11,
                restart=120,
                maxiter=3000,
            ).currents(plus)
            p_minus = TangentFFTVIE(
                centers,
                radii,
                spacing,
                k,
                fill_quadrature=6,
                gmres_rtol=1e-11,
                restart=120,
                maxiter=3000,
            ).currents(minus)
            fd = (p_plus - p_minus) / (2 * delta)
            analytic = dp[:, :, region]
            rel = float(
                np.linalg.norm(analytic - fd) / np.linalg.norm(fd)
            )
            rows.append(
                {
                    "name": f"{name}_region{region}_finite_difference",
                    "relative_error": rel,
                    "delta": delta,
                    "passed": bool(rel < 2e-4),
                }
            )
    return rows


def adapter_points_work_checks() -> list[dict]:
    """The adapter exposes points and parent-compatible work ledgers."""
    rows: list[dict] = []
    centers, radii, eps = PAIR
    spacing, k = 0.06, 6.3
    parent = DipoleVIE(centers, radii, spacing, k, fill_quadrature=6)
    fft = TangentFFTVIE(
        centers,
        radii,
        spacing,
        k,
        fill_quadrature=6,
        gmres_rtol=1e-11,
        restart=120,
        maxiter=3000,
    )
    _ = fft.currents(eps, derivatives=True)
    same_points = np.allclose(parent.points, fft.points)
    same_labels = np.array_equal(parent.labels, fft.labels)
    same_fill = np.allclose(parent.fill, fft.fill)
    work = fft.work
    compatible = set(parent.work) <= set(work)
    # Base solve fields are populated by every GMRES column.
    solve_count = len(fft.gmres_flags)
    rows.append(
        {
            "name": "adapter_points_and_work_ledger",
            "points_match": bool(same_points),
            "labels_match": bool(same_labels),
            "fill_match": bool(same_fill),
            "parent_work_keys_subset": bool(compatible),
            "work": dict(work),
            "fft_work": dict(fft.fft_work),
            "solves": solve_count,
            "matvecs": fft.matvecs,
            "true_check_matvecs": fft.true_check_matvecs,
            "passed": bool(
                same_points
                and same_labels
                and same_fill
                and compatible
                and work["factorizations"] == 0
                and work["rhs_columns"] == 12
                and solve_count == 12
                and fft.matvecs > solve_count
                and fft.true_check_matvecs == solve_count
            ),
        }
    )
    return rows


def joint_model_ordering_checks() -> list[dict]:
    z = np.zeros(14)
    z[:2] = [2.4, 3.0]
    z[2:6] = [.02, -.01, .005, .03]
    z[6:] = np.linspace(-.08, .09, 8)
    dense = parent_cal.JointModel(.06)
    fft = FFTJointModel(.06)
    y, j = dense.field_jac(z)
    yf, jf = fft.field_jac(z)
    ye = float(np.linalg.norm(y-yf)/np.linalg.norm(y))
    je = float(np.linalg.norm(j-jf)/np.linalg.norm(j))
    return [{"name":"joint_ordering_dense_fft", "prediction_error":ye,
             "jacobian_error":je, "passed":ye < 1e-7 and je < 1e-7}]


def main() -> None:
    checks = (
        dense_parent_currents_checks()
        + finite_difference_checks()
        + adapter_points_work_checks()
        + joint_model_ordering_checks()
    )
    report = {
        "scope": (
            "bounded discrete-model verification of FFT material tangent "
            "against the dense parent and central differences; not a "
            "continuum or production gate"
        ),
        "checks": checks,
        "passed": bool(checks and all(check["passed"] for check in checks)),
        "n_checks": len(checks),
    }
    (OUT / "tangent_checks_parent.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    print(json.dumps(report, indent=2))
    assert report["passed"]


if __name__ == "__main__":
    main()
