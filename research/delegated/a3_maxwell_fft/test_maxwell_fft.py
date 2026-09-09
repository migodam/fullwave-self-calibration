"""Test suite for the matrix-free FFT dyadic-convolution VIE solver.

Every test compares against the read-only parent dense implementation
``research/trispace_self_calibration/a3_research/maxwell3d.py``.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from maxwell_fft import (  # noqa: E402
    FFTVIE,
    DipoleVIE,
    dipole_kernel,
    receivers,
    relative_field_error,
    treams_field,
)

OUT = Path(__file__).resolve().parent

SINGLE = ([[0.0, 0.0, 0.0]], [0.12], [2.4 + 0.03j])
PAIR = (
    [[-0.2, 0.0, 0.0], [0.2, 0.03, 0.04]],
    [0.12, 0.10],
    [2.4 + 0.03j, 3.0 + 0.04j],
)
SHIFTED = (
    [[0.113, -0.071, 0.131]],
    [0.085],
    [2.7 + 0.05j],
)


def dense_fft_checks(seed: int = 9101) -> list[dict]:
    """FFT/dense kernel equality on full and occupied grids, plus reciprocity."""
    rng = np.random.default_rng(seed)
    rows: list[dict] = []
    configs = [
        ("single_origin", SINGLE[:2], SINGLE[2], 0.06, 7.3),
        ("pair_parent", PAIR[:2], PAIR[2], 0.06, 9.0),
        ("shifted_off_origin", SHIFTED[:2], SHIFTED[2], 0.045, 5.1),
    ]
    for name, (centers, radii), eps, spacing, k in configs:
        model = FFTVIE(centers, radii, spacing, k, fill_quadrature=6)
        n_rect = model.kernel.n_rect
        full_points = model.kernel.full_points
        source_full = (
            rng.normal(size=(n_rect, 3)) + 1j * rng.normal(size=(n_rect, 3))
        )
        source_grid = source_full.reshape(model.rect_shape + (3,))
        dense_matrix = dipole_kernel(full_points, full_points, model.k)
        dense_out = dense_matrix @ source_full.reshape(-1)
        fft_out = model.kernel.convolve(source_grid).reshape(-1)
        rel_full = float(
            np.linalg.norm(fft_out - dense_out) / np.linalg.norm(dense_out)
        )
        rows.append(
            {
                "name": name + "_fft_vs_dense_full",
                "rect_shape": list(model.rect_shape),
                "padded_shape": list(model.pad_shape),
                "full_cells": n_rect,
                "voxels": model.kernel.n_vox,
                "relative_error": rel_full,
                "passed": rel_full < 1e-11,
            }
        )

        # Occupied-grid operator equality.
        dense_occ = dipole_kernel(model.points, model.points, model.k)
        x_occ = source_full[model.kernel.active_flat]
        dense_occ_out = dense_occ @ x_occ.reshape(-1)
        fft_occ_out = model.kernel.from_full_grid(
            model.kernel.convolve(model.kernel.to_full_grid(x_occ))
        ).reshape(-1)
        rel_occ = float(
            np.linalg.norm(fft_occ_out - dense_occ_out)
            / np.linalg.norm(dense_occ_out)
        )
        rows.append(
            {
                "name": name + "_fft_vs_dense_occupied",
                "relative_error": rel_occ,
                "passed": rel_occ < 1e-11,
            }
        )

        # Reciprocity identities: u^T K v = v^T K u (ordinary transpose, as in
        # the parent dyadic reciprocity test).
        other_full = (
            rng.normal(size=(n_rect, 3)) + 1j * rng.normal(size=(n_rect, 3))
        )
        kv = model.kernel.convolve(
            other_full.reshape(model.rect_shape + (3,))
        ).reshape(-1)
        ku = model.kernel.convolve(source_grid).reshape(-1)
        left = np.dot(source_full.reshape(-1), kv)
        right = np.dot(other_full.reshape(-1), ku)
        rec_full = abs(left - right) / max(abs(left), abs(right), 1e-30)
        rows.append(
            {
                "name": name + "_reciprocity_full",
                "relative_error": float(rec_full),
                "passed": float(rec_full) < 1e-10,
            }
        )

        a_occ = source_full[model.kernel.active_flat].reshape(-1)
        b_occ = other_full[model.kernel.active_flat].reshape(-1)
        kv_occ = model.kernel.from_full_grid(
            model.kernel.convolve(model.kernel.to_full_grid(b_occ.reshape(-1, 3)))
        ).reshape(-1)
        ku_occ = model.kernel.from_full_grid(
            model.kernel.convolve(model.kernel.to_full_grid(a_occ.reshape(-1, 3)))
        ).reshape(-1)
        left = np.dot(a_occ, kv_occ)
        right = np.dot(b_occ, ku_occ)
        rec_occ = abs(left - right) / max(abs(left), abs(right), 1e-30)
        rows.append(
            {
                "name": name + "_reciprocity_occupied",
                "relative_error": float(rec_occ),
                "passed": float(rec_occ) < 1e-10,
            }
        )

        # Component ordering of the FFT block operator.
        # The operator includes per-component alpha scaling; compare it with
        # the dense operator so ordering is exercised component by component.
        model.set_eps(eps)
        probe = np.zeros(model.ndof, dtype=complex)
        probe[0] = 1.0
        column = model.apply_operator(probe, count=False)
        alpha_rep = np.repeat(model._alpha_values, 3)
        dense_eye = np.eye(model.ndof, dtype=complex) - alpha_rep[:, None] * dense_occ
        dense_column = dense_eye @ probe
        rel_col = float(
            np.linalg.norm(column - dense_column) / np.linalg.norm(dense_column)
        )
        rows.append(
            {
                "name": name + "_operator_column_ordering",
                "relative_error": rel_col,
                "passed": rel_col < 1e-10,
            }
        )
    return rows


def alpha_incident_checks() -> list[dict]:
    rows: list[dict] = []
    for name, (centers, radii, eps) in [
        ("single", SINGLE),
        ("pair", PAIR),
        ("shifted", SHIFTED),
    ]:
        spacing, k = 0.06, 6.3
        parent = DipoleVIE(centers, radii, spacing, k, fill_quadrature=6)
        fft_model = FFTVIE(centers, radii, spacing, k, fill_quadrature=6)
        same_grid = np.allclose(parent.points, fft_model.points, atol=0.0)
        same_labels = np.array_equal(parent.labels, fft_model.labels)
        same_fill = np.allclose(parent.fill, fft_model.fill)
        same_incident = np.allclose(parent.incident, fft_model.incident)
        # alpha from the same closed-form definition used by parent currents().
        e = np.asarray(eps)[fft_model.labels]
        a0 = 3 * spacing**3 * fft_model.fill * (e - 1) / (e + 2)
        alpha = a0 / (1 - 1j * k**3 * a0 / (6 * np.pi))
        fft_model.set_eps(eps)
        same_alpha = np.allclose(alpha, fft_model._alpha_values)
        passed = all(
            [same_grid, same_labels, same_fill, same_incident, same_alpha]
        )
        rows.append(
            {
                "name": name + "_embedding_polarizability_incident_match",
                "grid_match": bool(same_grid),
                "labels_match": bool(same_labels),
                "fill_match": bool(same_fill),
                "incident_match": float(
                    np.linalg.norm(parent.incident - fft_model.incident)
                    / np.linalg.norm(parent.incident)
                )
                if same_incident
                else None,
                "alpha_match": bool(same_alpha),
                "voxels": parent.points.shape[0],
                "passed": passed,
            }
        )
    return rows


def gmres_dense_checks() -> list[dict]:
    rows: list[dict] = []
    cases = [
        ("single_k9", SINGLE, 0.06, 9.0),
        ("pair_k6", PAIR, 0.06, 6.3),
    ]
    rx = receivers()
    for name, (centers, radii, eps), spacing, k in cases:
        parent = DipoleVIE(centers, radii, spacing, k, fill_quadrature=6)
        fft_model = FFTVIE(
            centers,
            radii,
            spacing,
            k,
            fill_quadrature=6,
            gmres_rtol=1e-10,
            restart=120,
            maxiter=1000,
        )
        p_dense = parent.currents(eps)
        p_fft = fft_model.solve(eps)
        rel_current = float(
            np.linalg.norm(p_fft - p_dense) / np.linalg.norm(p_dense)
        )
        stats = fft_model.solve_stats_row(eps)
        rows.append(
            {
                "name": name + "_gmres_vs_dense_currents",
                "relative_current_error": rel_current,
                "true_relative_residuals": stats["true_relative_residuals"],
                "iterations": stats["iterations"],
                "matvecs": stats["matvecs"],
                "solve_wall_seconds": stats["solve_wall_seconds"],
                "memory_estimate_mib": stats["memory_estimate"][
                    "estimated_peak_bytes"
                ]
                / 2**20,
                "passed": bool(
                    rel_current < 1e-7
                    and all(r < 1e-8 for r in stats["true_relative_residuals"])
                    and stats["gmres_converged"]
                ),
            }
        )
        field_parent = parent.field(eps, rx)
        # Rebuild a fresh model to keep solve counts in one place.
        fft_field_model = FFTVIE(
            centers,
            radii,
            spacing,
            k,
            fill_quadrature=6,
            gmres_rtol=1e-10,
            restart=120,
            maxiter=1000,
        )
        field_fft = fft_field_model.field(eps, rx)
        rel_field = relative_field_error(field_fft, field_parent)
        rows.append(
            {
                "name": name + "_gmres_vs_dense_field",
                "relative_field_error": rel_field,
                "passed": bool(rel_field < 1e-7),
            }
        )
    return rows


def parent_fill_consistency_checks() -> list[dict]:
    """The GMRES route reproduces the parent dense fill-study errors at k=9."""
    parent_rows = {
        (row["scene"], row["spacing"]): row
        for row in json.loads(
            (
                OUT.parents[2]
                / "research"
                / "trispace_self_calibration"
                / "a3_research"
                / "results"
                / "maxwell3d_fill_development.json"
            ).read_text()
        )
    }
    rows: list[dict] = []
    rx = receivers()
    scenes = [("single", SINGLE), ("pair", PAIR)]
    spacing, k = 0.03, 9.0
    for name, (centers, radii, eps) in scenes:
        reference, _ = treams_field(centers, radii, eps, k, rx, 5)
        model = FFTVIE(
            centers,
            radii,
            spacing,
            k,
            fill_quadrature=6,
            gmres_rtol=1e-10,
            restart=150,
            maxiter=1200,
        )
        predicted = model.field(eps, rx)
        error = relative_field_error(predicted, reference)
        expected = parent_rows[(name, spacing)]["relative_error"]
        rows.append(
            {
                "name": name + "_k9_h003_parent_fill_consistency",
                "relative_field_error": error,
                "parent_relative_field_error": expected,
                "absolute_difference": float(abs(error - expected)),
                "passed": bool(abs(error - expected) < 2e-6),
            }
        )
    return rows


def main() -> None:
    checks = (
        dense_fft_checks()
        + alpha_incident_checks()
        + gmres_dense_checks()
        + parent_fill_consistency_checks()
    )
    report = {
        "scope": (
            "development matrix-free FFT dyadic convolution + GMRES; "
            "bounded comparisons against read-only parent dense VIE, not a "
            "continuum-convergence or production acceptance gate"
        ),
        "checks": checks,
        "passed": bool(checks and all(check["passed"] for check in checks)),
        "n_checks": len(checks),
    }
    (OUT / "checks.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    assert report["passed"]


if __name__ == "__main__":
    main()
