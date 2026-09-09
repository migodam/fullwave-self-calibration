"""Repeatable N16 timing microbenchmark for reduced-work calibration.

Runs only on tuning seeds / before final freeze; saves raw samples and
conservative factors in ``calibration_tune.json``.
"""
from __future__ import annotations

import time
from typing import Any

import numpy as np
from scipy.linalg import qr, solve_triangular, svd

from common import Calibration, write_json


def _batch_timer(fn, reps: int):
    # warm up
    fn()
    ts = []
    for _ in range(int(reps)):
        t0 = time.perf_counter()
        fn()
        ts.append(time.perf_counter() - t0)
    return ts


def run_microbenchmark(model: Any, reps: int = 5, out_path=None) -> Calibration:
    n = model.n_cells
    rng = np.random.default_rng(12345)
    cal = Calibration(nominal_n=int(n))
    raw: dict[str, Any] = {}

    # full RHS solves and factorisations as timed by the physical core
    rhs_times: list[float] = []
    lu_times: list[float] = []
    for _ in range(reps):
        fw = model.forward(
            np.full(model.n_alpha, 0.5), np.zeros(3), jacobian=True
        )
        w = fw["work"]
        rhs_times.append(float(w["rhs_solve_seconds"]) / max(1, int(w["rhs_solves_total"])))
        lu_times.append(
            float(w["factorization_seconds"]) / max(1, int(w["factorizations"]))
        )
    cal.full_rhs_seconds = float(np.median(rhs_times))
    cal.lu_factor_seconds = float(np.median(lu_times))
    raw["full_rhs_per_solve"] = rhs_times
    raw["lu_per_factor"] = lu_times

    # dense n x n matvec (operator product estimate)
    D = model._domain_operator(3).copy()
    v = rng.standard_normal(n) + 1j * rng.standard_normal(n)

    def matvec():
        D @ v

    matvec()
    mv = _batch_timer(matvec, max(3, reps * 3))
    cal.dense_matvec_seconds = float(np.median(mv))
    raw["dense_matvec"] = mv

    # reduced QR of n x r and triangular solve r x r for r in ladder
    ranks = sorted({8, 16, 24, 36, 64, 96, 128, 192, 256} & set(range(1, n + 1)))
    for r in ranks:
        U = rng.standard_normal((n, r)) + 1j * rng.standard_normal((n, r))
        U, _ = qr(U, mode="economic")
        R = U - (0.5 * rng.standard_normal(n))[:, None] * (
            rng.standard_normal((n, r)) + 1j * rng.standard_normal((n, r))
        )
        rhs = rng.standard_normal((n, 12)) + 1j * rng.standard_normal((n, 12))

        def qr_fn():
            return qr(R, mode="economic")

        qr_fn()
        Q, Rr = qr(R, mode="economic")
        qts = _batch_timer(qr_fn, max(2, reps))
        cal.reduced_qr_samples[r] = float(np.median(qts))
        raw[f"qr_r{r}"] = qts

        def tri_fn():
            return solve_triangular(Rr, Q.conj().T @ rhs, lower=False)

        tri_fn()
        sts = _batch_timer(tri_fn, max(3, reps * 2))
        # 12 RHS through the r x r factor
        cal.reduced_solve_samples[r] = float(np.median(sts)) / 12.0
        raw[f"tri_r{r}"] = sts

    # sensing-stack economy SVD for realistic row counts
    for m in (36, 72, 108, 144):
        A = rng.standard_normal((m, n)) + 1j * rng.standard_normal((m, n))

        def svd_fn():
            return svd(A, full_matrices=False)

        svd_fn()
        s_ts = _batch_timer(svd_fn, max(2, reps))
        cal.svd_samples[m] = float(np.median(s_ts))
        raw[f"svd_m{m}"] = s_ts

    # RRQR of projected domain vectors
    for k in (16, 36, 64):
        A = rng.standard_normal((n, k)) + 1j * rng.standard_normal((n, k))

        def piv_fn():
            return qr(A, mode="economic", pivoting=True)

        piv_fn()
        p_ts = _batch_timer(piv_fn, max(2, reps))
        cal.rrqr_samples[k] = float(np.median(p_ts))
        raw[f"rrqr_k{k}"] = p_ts

    cal.raw_samples = raw
    if out_path is not None:
        write_json(out_path, cal.to_dict())
    return cal


if __name__ == "__main__":
    import sys

    sys.path.insert(0, "research/delegated/a2_physics")
    from physics import Config, Model

    out = sys.argv[1] if len(sys.argv) > 1 else "research/delegated/a2_solver/calibration_tune.json"
    reps = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    cal = run_microbenchmark(Model(Config(N=16, aperture="full")), reps=reps, out_path=out)
    print(cal.to_dict())
