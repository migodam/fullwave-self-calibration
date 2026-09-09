"""Reproduction script for the matrix-free FFT VIE solver and Treams audit.

Usage (from the repository root):

    research/trispace_self_calibration/a3_research/.venv3d/bin/python \
        research/delegated/a3_maxwell_fft/run_maxwell_fft.py

Steps:
1. independent validation checks against the read-only parent dense VIE;
2. short independent Treams sphere/cluster convergence audit at k = 9 and 18
   for grid spacings .03, .02, .015;
3. preservation of raw JSON (validation + audit), checks.json, environment
   and parent-file hashes.

The audit stops (preserving already-written rows) if a single GMRES solve
exceeds 120 s or the documented peak-memory estimate exceeds 4 GiB.
"""

from __future__ import annotations

import hashlib
import json
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(OUT))

from maxwell_fft import (  # noqa: E402
    FFTVIE,
    pointwise_errors,
    receivers,
    relative_field_error,
    treams_field,
)
import test_maxwell_fft as tests  # noqa: E402

PARENT_A3 = OUT.parents[2] / "research" / "trispace_self_calibration" / "a3_research"

SCENES = [
    (
        "single",
        [[0.0, 0.0, 0.0]],
        [0.12],
        [2.4 + 0.03j],
    ),
    (
        "pair",
        [[-0.2, 0.0, 0.0], [0.2, 0.03, 0.04]],
        [0.12, 0.10],
        [2.4 + 0.03j, 3.0 + 0.04j],
    ),
]

# Same geometry/material states as the parent fill study; reference order pairs
# are lmax5 vs 8 at k=9 (low) and lmax8 vs 10 at k=18 (high).
REFERENCE_LMAX = {9.0: (5, 8), 18.0: (8, 10)}
SPACINGS = [0.03, 0.02, 0.015]
FILL_QUADRATURE = 6
MEMORY_LIMIT_BYTES = 4 * 2**30
SOLVE_LIMIT_SECONDS = 120.0


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def environment_record() -> dict[str, Any]:
    return {
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "scipy": __import__("scipy").__version__,
        "cwd": str(Path.cwd()),
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
    }


def parent_hashes() -> dict[str, str]:
    names = [
        "maxwell3d.py",
        "validate_maxwell3d_fill.py",
        "calibrate3d.py",
        "results/maxwell3d_fill_development.json",
    ]
    return {name: sha256(PARENT_A3 / name) for name in names}


def validation_report() -> dict[str, Any]:
    checks = (
        tests.dense_fft_checks()
        + tests.alpha_incident_checks()
        + tests.gmres_dense_checks()
        + tests.parent_fill_consistency_checks()
    )
    return {
        "scope": (
            "independent bounded validation of matrix-free FFT dyadic "
            "convolution + GMRES against the read-only parent dense VIE; "
            "small grids only"
        ),
        "checks": checks,
        "passed": bool(checks and all(check["passed"] for check in checks)),
        "n_checks": len(checks),
    }


def audit_case(
    scene: tuple[str, list[list[float]], list[float], list[complex]],
    k: float,
    spacing: float,
    rx: np.ndarray,
) -> tuple[dict[str, Any], bool]:
    """One audit case; returns ``(row, stop_audit)``."""
    name, centers, radii, eps = scene
    low_lmax, high_lmax = REFERENCE_LMAX[k]
    row: dict[str, Any] = {
        "scene": name,
        "k": k,
        "resolution_m": spacing,
        "spacing": spacing,
        "fill_quadrature": FILL_QUADRATURE,
        "reference_lmax_low": low_lmax,
        "reference_lmax_high": high_lmax,
        "status": "ok",
    }
    case_start = time.perf_counter()

    model = FFTVIE(
        centers,
        radii,
        spacing,
        k,
        fill_quadrature=FILL_QUADRATURE,
        gmres_rtol=1e-10,
        gmres_atol=0.0,
        restart=80,
        maxiter=3000,
        wall_limit_seconds=SOLVE_LIMIT_SECONDS,
    )
    memory = model.memory_estimate()
    row.update(
        voxels=model.kernel.n_vox,
        dof=model.ndof,
        rect_shape=list(model.rect_shape),
        padded_shape=list(model.pad_shape),
        padded_cells=model.kernel.n_pad,
        fill_volume=float(model.fill.sum() * spacing**3),
        memory_estimate_mib=memory["estimated_peak_bytes"] / 2**20,
        memory_components={
            key: float(value) / 2**20 for key, value in memory.items()
        },
        volume_relative_error=float(
            abs(
                model.fill.sum() * spacing**3
                - sum(4 * np.pi * np.asarray(radii) ** 3 / 3)
            )
            / sum(4 * np.pi * np.asarray(radii) ** 3 / 3)
        ),
    )
    if memory["estimated_peak_bytes"] > MEMORY_LIMIT_BYTES:
        row["status"] = "stopped_memory_guard"
        row["elapsed_seconds"] = time.perf_counter() - case_start
        return row, True

    references: dict[str, np.ndarray] = {}
    xs_by_lmax: dict[int, list[list[float]]] = {}
    for label, lmax in [("low", low_lmax), ("high", high_lmax)]:
        start = time.perf_counter()
        ref, xs = treams_field(centers, radii, eps, k, rx, lmax)
        references[label] = np.asarray(ref)
        xs_by_lmax[lmax] = [[float(v) for v in pair] for pair in xs]
        row["reference_" + label + "_wall_seconds"] = (
            time.perf_counter() - start
        )
    row["xs_sca_ext_low"] = xs_by_lmax[low_lmax]
    row["xs_sca_ext_high"] = xs_by_lmax[high_lmax]
    passive_ok = all(
        sca <= ext + 1e-12
        for sca, ext in xs_by_lmax[high_lmax]
    )
    row["passive_extinction_ok_high"] = passive_ok

    solve_start = time.perf_counter()
    try:
        currents = model.solve(eps)
    except RuntimeError as exc:  # wall-limit guard
        row["status"] = "stopped_solve_time_guard"
        row["error"] = str(exc)
        row["elapsed_seconds"] = time.perf_counter() - case_start
        return row, True
    solve_wall = time.perf_counter() - solve_start
    stats = model.solve_stats_row(eps)
    stats.pop("scene_eps", None)
    row.update(
        currents_shape=list(currents.shape),
        solve_wall_seconds=solve_wall,
        gmres=stats,
    )
    predicted = (
        tests.dipole_kernel(rx, model.points, k) @ currents
    ).reshape(len(rx), 3, -1)

    row["reference_relative_convergence"] = relative_field_error(
        references["low"], references["high"]
    )
    row["relative_field_error_low"] = relative_field_error(
        predicted, references["low"]
    )
    row["relative_field_error_high"] = relative_field_error(
        predicted, references["high"]
    )
    row["pointwise_error_low"] = pointwise_errors(
        predicted, references["low"]
    )
    row["pointwise_error_high"] = pointwise_errors(
        predicted, references["high"]
    )
    row["reference_shape"] = list(references["high"].shape)
    row["predicted_shape"] = list(predicted.shape)
    row["elapsed_seconds"] = time.perf_counter() - case_start
    return row, False


def run_audit() -> list[dict[str, Any]]:
    rx = receivers()
    rows: list[dict[str, Any]] = []
    stop = False
    for scene in SCENES:
        for k in sorted(REFERENCE_LMAX):
            for spacing in SPACINGS:
                if stop:
                    row = {
                        "scene": scene[0],
                        "k": k,
                        "spacing": spacing,
                        "status": "skipped_after_guard",
                    }
                    rows.append(row)
                    continue
                row, stop = audit_case(scene, k, spacing, rx)
                rows.append(row)
                (OUT / "treams_audit_raw.json").write_text(
                    json.dumps(rows, indent=2) + "\n"
                )
                print(json.dumps(row), flush=True)
    return rows


def main() -> None:
    t0 = time.perf_counter()
    validation = validation_report()
    (OUT / "validation_raw.json").write_text(
        json.dumps(validation, indent=2) + "\n"
    )
    (OUT / "checks.json").write_text(json.dumps(validation, indent=2) + "\n")
    if not validation["passed"]:
        raise SystemExit("validation checks failed; no audit run")
    print("validation passed", validation["n_checks"], "checks")

    audit = run_audit()
    (OUT / "treams_audit_raw.json").write_text(
        json.dumps(audit, indent=2) + "\n"
    )
    env = environment_record()
    env.update(
        parent_file_sha256=parent_hashes(),
        reproduction_command=(
            "research/trispace_self_calibration/a3_research/.venv3d/bin/python "
            "research/delegated/a3_maxwell_fft/run_maxwell_fft.py"
        ),
        total_wall_seconds=time.perf_counter() - t0,
        n_audit_cases=len(audit),
        audit_statuses={
            status: sum(row["status"] == status for row in audit)
            for status in sorted({row["status"] for row in audit})
        },
    )
    (OUT / "environment.json").write_text(json.dumps(env, indent=2) + "\n")
    print("environment", json.dumps(env, indent=2))
    print("audit cases", len(audit))


if __name__ == "__main__":
    main()
