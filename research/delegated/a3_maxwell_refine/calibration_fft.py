"""Independent parent-calibration driver with the FFT tangent adapter.

The parent ``calibrate3d`` estimator is imported read-only and driven as a
callable (``parent_cal.fit``) against a JointModel-compatible wrapper whose
models are ``TangentFFTVIE`` matrix-free adapters.  No parent source is
modified; the wrapper reproduces the parent ``JointModel.field_jac`` ordering
so that predictions, tangent columns, fitted indices and work-ledger keys are
identical.

Scope remains the parent's stated DEVELOPMENT scope: known two-region support,
world-anchored illumination and Rx-translation calibration, not general 3D
imaging or hardware.  No true gains are given to the estimator; the
``unified_reference`` method receives the same independent noisy electronics
observation that the parent reference-study driver constructs.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence

import numpy as np

OUT = Path(__file__).resolve().parent
_PARENT_DIR = OUT.parents[2] / "research" / "trispace_self_calibration" / "a3_research"
sys.path.insert(0, str(_PARENT_DIR))
sys.path.insert(0, str(OUT))

import calibrate3d as parent_cal  # noqa: E402

from tangent_fft import TangentFFTVIE, dipole_kernel, receivers  # noqa: E402

CENTERS = parent_cal.CENTERS
RADII = parent_cal.RADII
LOSS = parent_cal.LOSS


class FFTJointModel:
    """JointModel-compatible wrapper over the FFT tangent VIE adapters."""

    def __init__(
        self,
        spacing: float = 0.01,
        frequencies: Sequence[float] | None = None,
        gmres_rtol: float = 1e-11,
        restart: int = 120,
        maxiter: int | None = 5000,
        wall_limit_seconds: float | None = 120.0,
    ) -> None:
        self.frequencies = np.asarray(
            [3.0, 6.0, 9.0] if frequencies is None else frequencies,
            dtype=float,
        )
        self.spacing = float(spacing)
        self.models = [
            TangentFFTVIE(
                CENTERS,
                RADII,
                self.spacing,
                float(k),
                fill_quadrature=6,
                gmres_rtol=gmres_rtol,
                restart=restart,
                maxiter=maxiter,
                wall_limit_seconds=wall_limit_seconds,
            )
            for k in self.frequencies
        ]
        self.rx = receivers()
        self.cached = None
        self.n_predict = 0

    def field_jac(self, z: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Prediction and exact Jacobian, byte-for-byte the parent ordering."""
        material = tuple(z[:2])
        if self.cached is None or self.cached[0] != material:
            self.cached = (
                material,
                [
                    m.currents(np.array(material) + 1j * LOSS, True)
                    for m in self.models
                ],
            )
        yy: list[np.ndarray] = []
        jj: list[np.ndarray] = []
        for k, m, (p, dp) in zip(
            self.frequencies, self.models, self.cached[1]
        ):
            pts = self.rx + z[2:5]
            op = dipole_kernel(pts, m.points, k)
            raw = (op @ p).reshape(len(pts), 3, 4)
            jac = np.zeros(raw.shape + (14,), complex)
            jac[..., :2] = (op @ dp.reshape(len(p), -1)).reshape(
                raw.shape + (2,)
            )
            for d in range(3):
                shift = np.eye(3)[d] * 1e-5
                dop = (
                    dipole_kernel(pts + shift, m.points, k)
                    - dipole_kernel(pts - shift, m.points, k)
                ) / (2e-5)
                jac[..., 2 + d] = (dop @ p).reshape(raw.shape)
            factor = np.exp(z[6:10] + 1j * z[10:14] + 1j * k * z[5])
            y = raw * factor
            jac = jac * factor[None, None, :, None]
            jac[..., 5] = 1j * k * y
            for t in range(4):
                jac[:, :, t, 6 + t] = y[:, :, t]
                jac[:, :, t, 10 + t] = 1j * y[:, :, t]
            yy.append(y)
            jj.append(jac)
        self.n_predict += 1
        return np.array(yy), np.array(jj)


def noisy_reference(
    z: np.ndarray, frequencies: np.ndarray, seed: int
) -> np.ndarray:
    """Independent noisy electronics observation, parent main's construction."""
    rng = np.random.default_rng(seed + 10000)
    reference = np.exp(
        z[6:10] + 1j * z[10:14] + 1j * frequencies[:, None] * z[5]
    )
    reference += 0.02 / np.sqrt(2.0) * (
        rng.normal(size=reference.shape)
        + 1j * rng.normal(size=reference.shape)
    )
    return reference


def write_rows(rows: list[dict[str, Any]]) -> None:
    """Append-style raw JSON preservation (rewrites after each row)."""
    (OUT / "calibration_fft_raw.json").write_text(
        json.dumps(rows, indent=2) + "\n"
    )


def run_one(
    rows: list[dict[str, Any]],
    seed: int,
    spacing: float,
    frequencies: list[float],
    method: str,
    max_nfev: int,
) -> dict[str, Any]:
    """Run one DEVELOPMENT fit; raises on wall/memory/resource guards."""
    old_ks = parent_cal.KS
    parent_cal.KS = np.asarray(frequencies, dtype=float)
    try:
        z, y, y0, sigma = parent_cal.truth(seed)
        reference = noisy_reference(z, parent_cal.KS, seed)
        model = FFTJointModel(
            spacing=spacing,
            frequencies=frequencies,
            wall_limit_seconds=120.0,
        )
        # Memory guard (documented estimate, per process).
        largest_mib = max(
            m.memory_estimate()["estimated_peak_bytes"] / 2**20
            for m in model.models
        )
        if largest_mib > 4096.0:
            row = {
                "status": "stopped_memory_guard",
                "seed": seed,
                "spacing": spacing,
                "frequencies": frequencies,
                "method": method,
                "largest_memory_estimate_mib": largest_mib,
                "scope": parent_cal.fit.__doc__ or "",
            }
            rows.append(row)
            write_rows(rows)
            return row
        before_work = [dict(m.work) for m in model.models]
        before_fft = [dict(m.fft_work) for m in model.models]
        start = time.perf_counter()
        row = parent_cal.fit(
            model,
            y,
            sigma,
            z,
            method,
            max_nfev=max_nfev,
            reference=reference if method == "unified_reference" else None,
        )
        # Cross-check the parent-fit ledger against our own deltas.
        manual_work = {
            key: float(
                sum(
                    m.work[key] - b[key]
                    for m, b in zip(model.models, before_work)
                )
            )
            for key in before_work[0]
        }
        fft_deltas = [
            {
                key: float(m.fft_work[key] - b[key])
                for key in b
            }
            for m, b in zip(model.models, before_fft)
        ]
        row.update(
            optimizer_status=row.get("status"),
            seed=seed,
            spacing=spacing,
            frequencies=list(frequencies),
            true=z.tolist(),
            noise_sigma=float(sigma),
            solver="fft_gmres_tangent",
            status="ok",
            total_elapsed_seconds=time.perf_counter() - start,
            work_ledger_match=bool(
                all(
                    abs(row["work"][key] - manual_work[key])
                    < 1e-9 * max(1.0, abs(manual_work[key]))
                    for key in manual_work
                )
            ),
            manual_work=manual_work,
            fft_work_deltas=fft_deltas,
            largest_memory_estimate_mib=largest_mib,
            per_frequency={
                "frequencies": list(frequencies),
                "voxels": [m.kernel.n_vox for m in model.models],
                "dof": [m.ndof for m in model.models],
                "rect_shape": [list(m.rect_shape) for m in model.models],
                "pad_shape": [list(m.pad_shape) for m in model.models],
                "memory_estimate_mib": [
                    m.memory_estimate()["estimated_peak_bytes"] / 2**20
                    for m in model.models
                ],
                "work": [dict(m.work) for m in model.models],
                "fft_work": [dict(m.fft_work) for m in model.models],
                "solves": [len(m.gmres_flags) for m in model.models],
                "matvecs": [m.matvecs for m in model.models],
                "true_check_matvecs": [
                    m.true_check_matvecs for m in model.models
                ],
                "true_relative_residuals": [
                    list(m.true_relative_residuals) for m in model.models
                ],
                "iterations": [list(m.iterations) for m in model.models],
                "gmres_flags": [list(m.gmres_flags) for m in model.models],
            },
        )
    except RuntimeError as exc:
        row = {
            "status": "stopped_resource_guard",
            "error": str(exc),
            "seed": seed,
            "spacing": spacing,
            "frequencies": list(frequencies),
            "method": method,
        }
        rows.append(row)
        write_rows(rows)
        return row
    finally:
        parent_cal.KS = old_ks
    rows.append(row)
    write_rows(rows)
    print(json.dumps(row), flush=True)
    return row


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--spacing", type=float, default=0.01)
    parser.add_argument("--max-nfev", type=int, default=45)
    parser.add_argument(
        "--frequency-sets",
        action="append",
        nargs="+",
        type=float,
        default=None,
        help="repeatable; default [[3,6,9],[3,9,18]]",
    )
    parser.add_argument("--seeds", nargs="+", type=int, default=[6101, 6102])
    parser.add_argument(
        "--methods",
        nargs="+",
        default=["unified", "unified_reference"],
    )
    parser.add_argument(
        "--frequency-sets-json",
        help="JSON file with [[...],[...]] frequency sets",
    )
    args = parser.parse_args()
    if args.frequency_sets_json:
        frequency_sets = json.loads(Path(args.frequency_sets_json).read_text())
    elif args.frequency_sets is not None:
        frequency_sets = [list(group) for group in args.frequency_sets]
    else:
        frequency_sets = [[3.0, 6.0, 9.0], [3.0, 9.0, 18.0]]

    output = OUT / "calibration_fft_raw.json"
    rows: list[dict[str, Any]] = json.loads(output.read_text()) if output.exists() else []
    for frequencies in frequency_sets:
        for seed in args.seeds:
            for method in args.methods:
                if any(r.get("seed") == seed and r.get("spacing") == args.spacing
                       and r.get("frequencies") == frequencies and r.get("method") == method
                       and r.get("status") == "ok" for r in rows):
                    continue
                row = run_one(
                    rows,
                    seed,
                    args.spacing,
                    [float(k) for k in frequencies],
                    method,
                    args.max_nfev,
                )
                if row.get("status", "ok").startswith("stopped_"):
                    print(
                        "guard checkpoint saved; not faking completion",
                        flush=True,
                    )
                    return
    summary = {
        "n_rows": len(rows),
        "statuses": {
            status: sum(r.get("status") == status for r in rows)
            for status in sorted({r.get("status") for r in rows})
        },
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "note": (
            "DEVELOPMENT calibration rows; raw JSON written row-by-row so "
            "guard interruptions leave an exact checkpoint."
        ),
    }
    (OUT / "calibration_fft_status.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
