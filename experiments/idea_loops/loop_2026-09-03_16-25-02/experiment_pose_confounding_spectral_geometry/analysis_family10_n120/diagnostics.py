#!/usr/bin/env python
"""Read-only diagnostics for results/family10_online_slam_toy_n120.json.

Does not run the experiment script; operates only on the stored JSON.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
JSON_PATH = ROOT / "results" / "family10_online_slam_toy_n120.json"
OUT = Path(__file__).resolve().parent / "diagnostics_summary.txt"


def fmt(x: float, nd: int = 4) -> str:
    return f"{x:.{nd}g}"


def med_mean_min_max(x: np.ndarray) -> tuple[str, str, str, str]:
    return (
        fmt(float(np.min(x))),
        fmt(float(np.median(x))),
        fmt(float(np.mean(x))),
        fmt(float(np.max(x))),
    )


def main() -> None:
    data = json.loads(JSON_PATH.read_text())
    mc = data["monte_carlo"]

    lines: list[str] = []
    def p(*args) -> None:
        s = " ".join(str(a) for a in args)
        print(s)
        lines.append(s)

    p("=" * 90)
    p("Family10 n120 trial diagnostics (stored JSON only)")
    p(f"file: {JSON_PATH}")
    p(f"generated_utc: {data['generated_utc']}")
    p("=" * 90)

    # ------------------------------------------------------------------
    # 0. Structure of stored per-trial records
    # ------------------------------------------------------------------
    p("\n[0] Stored fields per trial row (both modes):")
    p("  common: trial, known{success,status,nfev,cost,message,c_error_l2,seconds},")
    p("          free{success,status,nfev,cost,message,c_error_l2,dx_l2,seconds}")
    p("  c_error_l2 = ||c_hat - c0||_2  (map-estimation error norm, p=24)")
    p("  dx_l2      = ||dx_hat||_2      (free-pose MAP displacement norm, q=18;")
    p("              pose truth is dx=0, so this is the pose error norm)")
    p("  Per-trial map error VECTORS and dx_hat VECTORS are NOT stored;")
    p("  serialization drops errs_known/errs_free/ok_* (source line ~1523).")

    summary: dict[str, dict] = {}
    modes = ["born", "full_wave"]
    kinds = ["known", "free"]
    for mode in modes:
        rows = mc[mode]["trial_rows"]
        for kind in kinds:
            key = f"{mode}_{kind}"
            summary[key] = {
                "n_completed": len(rows),
                "n_success": sum(1 for r in rows if r[kind]["success"]),
                "success": np.array([r[kind]["success"] for r in rows], dtype=bool),
                "e2": np.array([r[kind]["c_error_l2"] ** 2 for r in rows], dtype=float),
                "e": np.array([r[kind]["c_error_l2"] for r in rows], dtype=float),
                "status": np.array([r[kind]["status"] for r in rows], dtype=int),
                "cost": np.array([r[kind]["cost"] for r in rows], dtype=float),
                "nfev": np.array([r[kind]["nfev"] for r in rows], dtype=int),
                "dx": np.array(
                    [r[kind].get("dx_l2", np.nan) for r in rows], dtype=float
                ),
                "rows": rows,
            }

    # ------------------------------------------------------------------
    # 1. Per-mode summary tables
    # ------------------------------------------------------------------
    p("\n[1] Per-mode trial/success and error-norm summary (successful trials only)")
    p(f"{'mode':<22}{'ntrials':>8}{'nsucc':>7}{'succ%':>7}"
      f"{'|e|^2 min':>14}{'median':>12}{'mean':>12}{'max':>12}")
    e2_rows = {}
    for key in [f"{m}_{k}" for m in modes for k in kinds]:
        s = summary[key]
        ok = s["success"]
        e2 = s["e2"][ok]
        e2_rows[key] = e2
        mn, md, me, mx = med_mean_min_max(e2)
        p(f"{key:<22}{s['n_completed']:>8}{s['n_success']:>7}"
          f"{100*s['n_success']/s['n_completed']:>7.1f}{mn:>14}{md:>12}{me:>12}{mx:>12}")

    p("\nSame, |e| (norm as stored)  and free-mode pose-error norm ||dx_hat||:")
    p(f"{'mode':<22}{'|e| min':>12}{'median':>12}{'mean':>12}{'max':>12}"
      f"{'|dx| min':>12}{'median':>12}{'mean':>12}{'max':>12}")
    for key in [f"{m}_{k}" for m in modes for k in kinds]:
        s = summary[key]
        ok = s["success"]
        e = s["e"][ok]
        dx = s["dx"][ok]
        dxvals = med_mean_min_max(dx) if np.isfinite(dx).all() else ("n/a",) * 4
        p(f"{key:<22}" + "".join(f"{v:>12}" for v in med_mean_min_max(e) + dxvals))
    p("  known-pose fits fix dx=0, so their pose error is identically zero (not stored).")

    # Failure rows (all trials)
    p("\nNon-success rows (all 120 trials used for n_completed):")
    for mode in modes:
        for kind in kinds:
            s = summary[f"{mode}_{kind}"]
            bad = np.where(~s["success"])[0]
            if bad.size:
                p(f"  {mode}_{kind}: {bad.size} failures, trials={bad.tolist()}")
                for i in bad:
                    r = s["rows"][i][kind]
                    p(f"    trial {i}: status={r['status']} nfev={r['nfev']} "
                      f"cost={r['cost']:.4g} e2={s['e2'][i]:.4g} "
                      f"dx_l2={r.get('dx_l2', np.nan):.4g} msg={r['message'][:60]}")
            else:
                p(f"  {mode}_{kind}: none")

    # ------------------------------------------------------------------
    # 2. Outliers (>20x median ||e_free||^2, successful trials)
    # ------------------------------------------------------------------
    p("\n[2] Free-pose outliers: ||e_free||^2 > 20 * median(successful ||e_free||^2)")
    outlier_info = {}
    for mode in modes:
        key = f"{mode}_free"
        known_key = f"{mode}_known"
        sf = summary[key]
        sk = summary[known_key]
        okf = sf["success"]
        okk = sk["success"]
        e2f = sf["e2"][okf]
        med = float(np.median(e2f))
        thr = 20.0 * med
        inds = np.where(okf & (sf["e2"] > thr))[0]
        outlier_info[key] = {"median_e2": med, "threshold": thr, "trials": inds.tolist()}
        p(f"\n  {mode} free: median e2={med:.5g}, threshold={thr:.5g}, "
          f"outliers={inds.size} of {okf.sum()} successful trials")
        if inds.size == 0:
            continue
        # How many of these are also huge in the known-pose fit?
        also_20x_known = 0
        p(f"  {'trial':>6}{'free e2':>13}{'known e2':>13}{'free/known':>11}"
          f"{'free |dx|':>11}{'free status':>12}{'known status':>13}"
          f"{'free nfev':>11}{'free cost':>12}")
        for i in inds:
            rf = sf["rows"][i]["free"]
            rk = sk["rows"][i]["known"]
            kmed = float(np.median(sk["e2"][okk]))
            also = sk["e2"][i] > 20.0 * kmed
            also_20x_known += int(also)
            ratio = sf["e2"][i] / sk["e2"][i] if sk["e2"][i] > 0 else np.inf
            p(f"  {i:>6}{sf['e2'][i]:>13.5g}{sk['e2'][i]:>13.5g}{ratio:>11.4g}"
              f"{rf['dx_l2']:>11.4g}{rf['status']:>12}{rk['status']:>13}"
              f"{rf['nfev']:>11}{rf['cost']:>12.5g}")
        p(f"  -> {also_20x_known} of {inds.size} outlier trials are also >20x "
          f"the known-pose median; the rest are free-pose-specific.")

    # Known-free inflation across all successful trials (helps interpretation)
    p("\n  Known-vs-free squared-error inflation on successful trials "
      "(per-trial free e2 / known e2):")
    for mode in modes:
        sf = summary[f"{mode}_free"]
        sk = summary[f"{mode}_known"]
        both = sf["success"] & sk["success"]
        ratio = sf["e2"][both] / np.maximum(sk["e2"][both], 1e-300)
        p(f"  {mode:<10} n={both.sum()} median={np.median(ratio):.4g} "
          f"mean={np.mean(ratio):.4g} max={np.max(ratio):.4g}")

    # ------------------------------------------------------------------
    # 3. Trace ratio excluding outliers
    # ------------------------------------------------------------------
    p("\n[3] Empirical trace ratio on non-outlier trials (exactness caveat)")
    for mode in modes:
        key = f"{mode}_free"
        info = outlier_info[key]
        s = summary[key]
        ok = s["success"]
        inds = info["trials"]
        keep = np.ones(ok.size, dtype=bool)
        keep[inds] = False
        e2_full_ok = s["e2"][ok]
        e2_clean = s["e2"][ok & keep]
        med = info["median_e2"]
        p(f"\n  {mode} free:")
        p(f"    full sample: n={ok.sum()}, sum|e|^2={e2_full_ok.sum():.6g}, "
          f"median|e|^2={med:.5g}")
        p(f"    clean (below 20x median): n={e2_clean.size}, "
          f"sum|e|^2={e2_clean.sum():.6g}, median|e|^2={np.median(e2_clean):.5g}")

        # Exact full-sample trace via identity: tr(Cov) =
        # (sum_i ||e_i||^2 - n ||mean||^2)/(n-1), checked against stored trace.
        mean_full = np.asarray(mc[mode]["mean_error_free"])
        n_full = ok.sum()
        trace_from_scalars = (e2_full_ok.sum() - n_full * (mean_full @ mean_full)) / (
            n_full - 1
        )
        stored_trace = mc[mode]["free_trace_emp"]
        p(f"    identity check tr(Cov_full): from stored norms+mean = {trace_from_scalars:.7g}, "
          f"stored = {stored_trace:.7g}")

        # For the clean subset the exact centered trace is NOT recoverable from
        # the stored JSON because the clean-subset mean vector is not stored.
        # Two usable estimates:
        n_clean = int(e2_clean.size)
        # (a) exact uncentered second-moment version (omit mean subtraction),
        # (b) subtract full-sample mean-vector squared norm (approximation).
        prox_centered = (e2_clean.sum() - n_clean * (mean_full @ mean_full)) / (
            n_clean - 1
        )
        p(f"    clean tr approx (full-sample mean used): {prox_centered:.6g}")
        p(f"    clean tr upper bound-ish uncentered: "
          f"{e2_clean.sum() / (n_clean - 1):.6g}")

        denom_trace = mc[mode]["known_trace_emp"]
        p(f"    stored full trace ratio={mc[mode]['trace_ratio_empirical']:.5g} "
          f"(predicted {mc[mode]['trace_ratio_predicted']:.5g})")
        p(f"    clean approx trace ratio (approx numerator / full known trace)="
          f"{prox_centered / denom_trace:.5g}")

    # ------------------------------------------------------------------
    # 4. Pose (dx) diagnostics, free modes
    # ------------------------------------------------------------------
    p("\n[4] Pose MAP displacement (dx_hat) statistics, free successful trials")
    p("  Stored per-trial pose data is scalar dx_l2 = ||dx_hat|| only; component-")
    p("  wise |dx| or dx_hat vectors are not stored.")
    for mode in modes:
        s = summary[f"{mode}_free"]
        dx = s["dx"][s["success"]]
        p(f"  {mode} free: n={dx.size}, mean||dx||={np.mean(dx):.5g}, "
          f"median={np.median(dx):.5g}, max||dx||={np.max(dx):.5g}")
        big = np.where(dx > np.median(dx) * 5)[0]
        if big.size:
            p(f"    trials with ||dx|| > 5x median: {big.tolist()}")

    # ------------------------------------------------------------------
    # 5. Born free: predicted vs empirical covariance scale
    # ------------------------------------------------------------------
    p("\n[5] Born free: predicted vs empirical covariance scale")
    pred_eig = np.asarray(data["theory"]["born"]["P_free_eigvals_desc"])
    pred_tr = float(pred_eig.sum())
    emp = np.asarray(mc["born"]["Cov_free"])
    emp_tr = float(np.trace(emp))
    emp_diag = np.diag(emp)
    p(f"  predicted P_free eigenvalues (desc): min={pred_eig.min():.5g} "
      f"median={np.median(pred_eig):.5g} max={pred_eig.max():.5g}")
    p(f"  predicted P_free trace = sum(eig) = {pred_tr:.7g}")
    p(f"  predicted P_free per-element DIAGONAL is not stored in the JSON "
      f"(only eigenvalues/trace), so it cannot be reported without re-running "
      f"the linearized algebra.")
    p(f"  empirical Cov_free trace = {emp_tr:.7g}")
    p(f"  empirical Cov_free diagonal: min={emp_diag.min():.5g} "
      f"median={np.median(emp_diag):.5g} mean={np.mean(emp_diag):.5g} "
      f"max={emp_diag.max():.5g}")
    p(f"  understatement factor (emp trace / pred trace) = {emp_tr / pred_tr:.4g}")
    p(f"  stored per-mode traces: born free_emp={mc['born']['free_trace_emp']:.7g} "
      f"pred={mc['born']['free_trace_pred']:.7g}; "
      f"fw free_emp={mc['full_wave']['free_trace_emp']:.7g} "
      f"pred={mc['full_wave']['free_trace_pred']:.7g}")

    OUT.write_text("\n".join(lines) + "\n")
    p(f"\n[written] {OUT}")


if __name__ == "__main__":
    main()
