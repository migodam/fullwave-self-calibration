"""Bounded full-wave Gaussian parameter-tangent SOM/TSOMG development probe.

Route A only: the current is eliminated by the full VIE at every iterate.
The second operator ranks parameter directions; it is never added as data.
Run with ``Gaussian/.venv_nn/bin/python Gaussian/code/a2/som.py``.
"""
from __future__ import annotations

import os
for _name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_name] = "2"

import argparse
import hashlib
import json
import platform
import sys
import time
from pathlib import Path
from dataclasses import dataclass

import numpy as np
import scipy
import scipy.linalg as la

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "code"))
from a2.physics import Geometry, VIE
from a2.ports import GaussianComponent, render

FREQUENCIES_HZ = (1.50e9, 2.75e9)
SCALES = np.tile(np.array([0.25, 0.04, 0.04, 0.35, 0.35, 1.0]), 2)
LOWER = np.tile(np.array([0.02, -0.135, -0.135, np.log(.012), np.log(.012), -np.pi]), 2)
UPPER = np.tile(np.array([1.20, .135, .135, np.log(.085), np.log(.085), np.pi]), 2)


@dataclass(frozen=True)
class FitConfig:
    """Reusable A2 parameter-tangent fit controls (all values are recorded)."""
    max_iter: int = 5
    first_rank: int = 6
    second_rank: int = 2
    rtol: float = 2e-9
    max_trial_steps: int = 5
    scaled_trust_radius: float = 1.0
    discrepancy_mse: float = 1.25
    gradient_inf_tolerance: float = 1e-5
    adaptive_rank: bool = False
    adaptive_gradient_deficit: float = .15
    adaptive_first_cap: int = 10
    adaptive_second_cap: int = 4


def components(theta: np.ndarray, frequency_hz: float) -> list[GaussianComponent]:
    """Known fixed loss law; theta itself is a shared real material vector."""
    eta = 1.0 + .02j * 2.25e9 / frequency_hz
    answer = []
    for a, cx, cy, lsx, lsy, angle in np.asarray(theta).reshape(-1, 6):
        answer.append(GaussianComponent(complex(a) * eta, (float(cx), float(cy)),
                                        (float(np.exp(lsx)), float(np.exp(lsy))), float(angle)))
    return answer


def chi_and_partials(theta: np.ndarray, vie: VIE, frequency_hz: float) -> tuple[np.ndarray, np.ndarray]:
    """Material and analytic dchi/dq in the declared physical metric."""
    theta = np.asarray(theta, float)
    chi = np.zeros(vie.points.shape[0], complex)
    partial = np.empty((theta.size, chi.size), complex)
    eta = 1.0 + .02j * 2.25e9 / frequency_hz
    x, y = vie.points.T
    for k, (a, cx, cy, lsx, lsy, angle) in enumerate(theta.reshape(-1, 6)):
        sx, sy = np.exp(lsx), np.exp(lsy); c, s = np.cos(angle), np.sin(angle)
        dx, dy = x-cx, y-cy; u, v = c*dx+s*dy, -s*dx+c*dy
        term = eta * a * np.exp(-.5*((u/sx)**2+(v/sy)**2)); chi += term
        raw = (term/a, term*(c*u/sx**2-s*v/sy**2), term*(s*u/sx**2+c*v/sy**2), term*u**2/sx**2, term*v**2/sy**2, term*u*v*(1/sy**2-1/sx**2))
        partial[6*k:6*k+6] = np.asarray(raw) * SCALES[6*k:6*k+6,None]
    return chi, partial


def realify(x: np.ndarray) -> np.ndarray:
    return np.concatenate((np.asarray(x).real, np.asarray(x).imag), axis=0)


def solve_stack(theta: np.ndarray, vies: list[VIE], *, rtol: float) -> tuple[list[np.ndarray], list[dict], list[np.ndarray]]:
    chi, fw, deriv = [], [], []
    for frequency, vie in zip(FREQUENCIES_HZ, vies):
        c, d = chi_and_partials(theta, vie, frequency)
        chi.append(c); fw.append(vie.forward(c, rtol=rtol)); deriv.append(d)
    return chi, fw, deriv


def tangent_operators(theta: np.ndarray, vies: list[VIE], train_rx: np.ndarray, noise_sigma: float, *, rtol: float):
    """Return full-state external A and internal B at one real parameter point."""
    chi, fw, dchi = solve_stack(theta, vies, rtol=rtol)
    acol, bcol, tangent_solves = [], [], 0
    norm_b = vies[0].h / np.sqrt(len(vies) * vies[0].E.shape[1] * vies[0].D.npix)
    for fi, vie in enumerate(vies):
        for p in range(theta.size):
            tan = vie.material_tangent(chi[fi], dchi[fi][p], forward=fw[fi], rtol=rtol)
            tangent_solves += vie.E.shape[1]
            acol.append((tan["scattered"][train_rx] / noise_sigma).ravel())
            bcol.append((norm_b * vie.D.matmat(tan["current"])).ravel())
    # Each frequency parameter block was appended sequentially; collect columns.
    # Recompute column order as p-major across all frequencies.
    a = np.empty((2 * len(vies) * int(np.count_nonzero(train_rx)) * vies[0].E.shape[1], theta.size))
    b = np.empty((2 * len(vies) * vies[0].D.npix * vies[0].E.shape[1], theta.size))
    for p in range(theta.size):
        aa = [acol[fi * theta.size + p] for fi in range(len(vies))]
        bb = [bcol[fi * theta.size + p] for fi in range(len(vies))]
        a[:, p] = realify(np.concatenate(aa))
        b[:, p] = realify(np.concatenate(bb))
    return a, b, chi, fw, {"full_forward_rhs": len(vies) * vies[0].E.shape[1], "tangent_rhs": tangent_solves}


def folds(a: np.ndarray, b: np.ndarray, first_target: int, second_target: int, relative_cutoff: float = 1e-8, absolute_data_tol: float = 1e-10, absolute_internal_tol: float = 1e-12):
    _, sa, vh = la.svd(a, full_matrices=False)
    v = vh.T
    numerical_rank = int(np.count_nonzero(sa > max(absolute_data_tol, (sa[0] if len(sa) else 0.) * relative_cutoff)))
    r1 = min(first_target, numerical_rank)
    v1 = v[:, :r1]
    vd = v[:, r1:]
    if vd.shape[1] == 0:
        return v1, np.empty((a.shape[1], 0)), {"first_rank": r1, "second_rank": 0, "data_singular_values": sa.tolist(), "internal_weak_singular_values": []}
    _, sb, zt = la.svd(b @ vd, full_matrices=False)
    r2 = min(second_target, int(np.count_nonzero(sb > max(absolute_internal_tol, (sb[0] if len(sb) else 0.) * relative_cutoff))))
    v2 = vd @ zt.T[:, :r2]
    return v1, v2, {"first_rank": r1, "second_rank": r2, "relative_singular_cutoff": relative_cutoff, "absolute_data_tolerance": absolute_data_tol, "absolute_internal_tolerance": absolute_internal_tol, "data_singular_values": sa.tolist(), "internal_weak_singular_values": sb.tolist(), "twofold_internal_energy_fraction": float(np.sum(sb[:r2]**2) / max(np.sum(sb**2), 1e-30))}


def field_stack(fw: list[dict], train_rx: np.ndarray) -> np.ndarray:
    return np.concatenate([q["scattered"][train_rx].ravel() for q in fw])


def state_residual(vies: list[VIE], chi: list[np.ndarray], fw: list[dict]) -> float:
    vals = []
    for vie, c, out in zip(vies, chi, fw):
        r = out["current"] - c[:, None] * (vie.E + vie.D.matmat(out["current"]))
        vals.append(np.linalg.norm(r) / max(np.linalg.norm(c[:, None] * vie.E), 1e-30))
    return float(max(vals))


def fit_parameter_tangent_som(name: str, theta0: np.ndarray, truth: np.ndarray, observed: list[np.ndarray], clean: list[np.ndarray], vies: list[VIE], train_rx: np.ndarray, held_rx: np.ndarray, noise_sigma: float, config: FitConfig = FitConfig()):
    """Fit one paired case with ordinary LM, first-fold SOM, or twofold TSOMG.

    This is the reusable inverse API for a frozen campaign.  It returns every
    accepted/capped iteration and actual VIE right-hand-side calls; callers own
    data generation and any independent-grid validation.
    """
    theta = theta0.copy(); damping = 1e-2; history = []; started = time.perf_counter(); cost = {"full_forward_rhs": 0, "tangent_rhs": 0, "attempts": 0}
    last_meta = None
    adaptive_first, adaptive_second = config.first_rank, config.second_rank
    for iteration in range(config.max_iter):
        a, b, chi, fw, ledger = tangent_operators(theta, vies, train_rx, noise_sigma, rtol=config.rtol)
        for key, value in ledger.items(): cost[key] += value
        target = np.concatenate([observed[i][train_rx].ravel() - fw[i]["scattered"][train_rx].ravel() for i in range(len(vies))]) / noise_sigma
        residual = realify(target)
        v1, v2, meta = folds(a, b, first_target=adaptive_first, second_target=adaptive_second)
        if name == "ordinary_lm": basis = np.eye(theta.size)
        elif name == "singlefold_som": basis = v1
        elif name == "twofold_tsomg": basis = np.column_stack((v1, v2))
        elif name in ("random_complement_tsomg", "gsvd_ratio_tsomg"):
            _, _, vh = la.svd(a, full_matrices=False); vd = vh.T[:, meta["first_rank"]:]
            r2 = meta["second_rank"]
            if name == "random_complement_tsomg":
                # Fixed seed by iterate: a comparable equal-dimensional null-complement control.
                q, _ = la.qr(np.random.default_rng(20260911 + iteration).normal(size=(vd.shape[1], r2)), mode="economic")
                extra = vd @ q
            else:
                av, bv = a @ vd, b @ vd; ridge = 1e-10 * max(float(np.linalg.norm(av, 2)**2), 1.)
                _, vec = la.eigh(bv.T @ bv, av.T @ av + ridge * np.eye(vd.shape[1]))
                extra = vd @ vec[:, -r2:]
                meta["gsvd_ratio_ridge"] = ridge
            basis = np.column_stack((v1, extra))
        elif name == "adaptive_tsomg":
            basis = np.column_stack((v1, v2)); full_gradient = a.T @ residual
            deficit = 1. - float(np.linalg.norm(basis.T @ full_gradient)**2 / max(np.linalg.norm(full_gradient)**2, 1e-30))
            change = "none"
            if deficit > config.adaptive_gradient_deficit:
                if adaptive_first < config.adaptive_first_cap and adaptive_first < theta.size:
                    adaptive_first += 1; change = "increase_first"
                elif adaptive_second < config.adaptive_second_cap:
                    adaptive_second += 1; change = "increase_second"
                if change != "none":
                    v1, v2, meta = folds(a, b, first_target=adaptive_first, second_target=adaptive_second)
                    basis = np.column_stack((v1, v2))
            meta.update({"adaptive_gradient_deficit": deficit, "adaptive_rank_change": change, "adaptive_requested_first": adaptive_first, "adaptive_requested_second": adaptive_second})
        else: raise ValueError(name)
        old_loss = float(np.mean(np.abs(target)**2)); gradient_inf = float(np.linalg.norm((a @ basis).T @ residual, ord=np.inf)); accepted = False; selected_fw = fw
        if old_loss <= config.discrepancy_mse or gradient_inf <= config.gradient_inf_tolerance:
            history.append({"iteration": iteration, "train_weighted_mse": old_loss, "accepted": False, "stop_reason": "discrepancy" if old_loss <= config.discrepancy_mse else "gradient", "gradient_inf": gradient_inf, "damping": damping, **meta})
            last_meta = meta; break
        for _ in range(config.max_trial_steps):
            z = la.solve((a @ basis).T @ (a @ basis) + damping * np.eye(basis.shape[1]), (a @ basis).T @ residual, assume_a="pos")
            dq = basis @ z; dq *= min(1., config.scaled_trust_radius / max(np.linalg.norm(dq), 1e-30))
            proposal = np.clip(theta + SCALES * dq, LOWER, UPPER)
            cnew, fnew, _ = solve_stack(proposal, vies, rtol=config.rtol); cost["full_forward_rhs"] += len(vies) * vies[0].E.shape[1]; cost["attempts"] += 1
            new_target = np.concatenate([observed[i][train_rx].ravel() - fnew[i]["scattered"][train_rx].ravel() for i in range(len(vies))]) / noise_sigma
            new_loss = float(np.mean(np.abs(new_target)**2))
            if new_loss < old_loss:
                theta = proposal; selected_fw = fnew; damping = max(damping / 3., 1e-8); accepted = True
                break
            damping = min(damping * 10., 1e8)
        held_num = np.sqrt(sum(np.linalg.norm(selected_fw[i]["scattered"][held_rx] - clean[i]["scattered"][held_rx])**2 for i in range(len(vies))))
        held_den = np.sqrt(sum(np.linalg.norm(clean[i]["scattered"][held_rx])**2 for i in range(len(vies))))
        history.append({"iteration": iteration, "train_weighted_mse": new_loss if accepted else old_loss, "accepted": accepted, "stop_reason": "accepted" if accepted else "rejected_step", "gradient_inf": gradient_inf, "damping": damping, "held_scattered_relative": float(held_num / held_den), **meta})
        last_meta = meta
        if not accepted: break
    chi, fw, _ = solve_stack(theta, vies, rtol=config.rtol); cost["full_forward_rhs"] += len(vies) * vies[0].E.shape[1]
    material = render(components(theta, FREQUENCIES_HZ[0]), vies[0].points).real
    true_material = render(components(truth, FREQUENCIES_HZ[0]), vies[0].points).real
    held_num = np.sqrt(sum(np.linalg.norm(fw[i]["scattered"][held_rx] - clean[i]["scattered"][held_rx])**2 for i in range(len(vies))))
    held_den = np.sqrt(sum(np.linalg.norm(clean[i]["scattered"][held_rx])**2 for i in range(len(vies))))
    terminal = history[-1]["stop_reason"] if history else "no_iteration"
    return {"method": name, "fit_config": config.__dict__, "parameters": theta.tolist(), "parameter_scaled_error": float(np.linalg.norm((theta-truth)/SCALES) / np.sqrt(theta.size)), "material_relative_l2": float(np.linalg.norm(material-true_material)/np.linalg.norm(true_material)), "held_scattered_relative": float(held_num / held_den), "state_relative_residual": state_residual(vies, chi, fw), "history": history, "termination": "iteration_cap" if len(history) == config.max_iter and history[-1]["accepted"] else terminal, "final_fold": last_meta, "cost": cost, "wall_seconds": time.perf_counter()-started}


def run_method(*args, **kwargs):
    """Compatibility alias; new campaign callers should use the named API."""
    return fit_parameter_tangent_som(*args, **kwargs)


def selection_diagnostics(theta: np.ndarray, vies: list[VIE], train_rx: np.ndarray, noise_sigma: float) -> dict:
    """Equal-dimensional random and generalized-ratio diagnostics at one point.

    They diagnose whether V2 carries more B-energy than a seeded random
    complement and how it compares with an A/B generalized ratio.  Neither is
    a reconstruction baseline or an identifiability claim.
    """
    a, b, _, _, _ = tangent_operators(theta, vies, train_rx, noise_sigma, rtol=2e-9)
    v1, v2, meta = folds(a, b, first_target=6, second_target=2)
    _, _, vh = la.svd(a, full_matrices=True); vd = vh.T[:, meta["first_rank"]:]
    rng = np.random.default_rng(20260911)
    q, _ = la.qr(rng.normal(size=(vd.shape[1], max(meta["second_rank"], 1))), mode="economic")
    vrand = vd @ q[:, :meta["second_rank"]]
    # Stable generalized B/A ratio in the same data-weak complement.
    av, bv = a @ vd, b @ vd
    ridge = 1e-10 * max(float(np.linalg.norm(av, 2)**2), 1.)
    vals, vecs = la.eigh(bv.T @ bv, av.T @ av + ridge * np.eye(vd.shape[1]))
    vgsvd = vd @ vecs[:, -meta["second_rank"]:]
    energy = lambda v: float(np.linalg.norm(b @ v, "fro")**2)
    data = lambda v: float(np.linalg.norm(a @ v, "fro")**2)
    return {"scope": "one-point equal-dimension selector diagnostic; not a GSVD reconstruction baseline", "seed": 20260911, "twofold": {"internal_energy": energy(v2), "data_energy": data(v2)}, "random_complement": {"internal_energy": energy(vrand), "data_energy": data(vrand)}, "generalized_internal_over_data": {"internal_energy": energy(vgsvd), "data_energy": data(vgsvd), "ridge": ridge, "generalized_eigenvalues_desc": vals[::-1].tolist()}, "fold": meta}


def frozen_campaign_call_ledger() -> dict:
    """Pre-execution ledger for the requested 3x100x3 paired-noise campaign."""
    objects, noises, methods = 3 * 100, 3, ("ordinary_lm", "singlefold_som", "twofold_tsomg")
    cfg = FitConfig()
    rhs_per_forward = len(FREQUENCIES_HZ) * 6
    tangent_rhs_per_outer = 12 * rhs_per_forward
    maximum_outer_rhs = rhs_per_forward + tangent_rhs_per_outer + cfg.max_trial_steps * rhs_per_forward
    return {"status": "prepared_not_executed", "campaign": {"families": 3, "objects_per_family": 100, "paired_noises_per_object": noises, "paired_cases": objects * noises, "methods": methods, "inverse_grid": "N64 requested", "data_grid": "N128 requested independently; awaiting frozen quadrature/data API from A2 ports work"}, "fit_config": cfg.__dict__, "accounting_unit": "full VIE right-hand-side solve; a multi-RHS GMRES call reports each RHS", "per_method_case_upper_bound": {"outer_iterations": cfg.max_iter, "forward_rhs_per_outer_for_current_state": rhs_per_forward, "tangent_rhs_per_outer": tangent_rhs_per_outer, "trial_forward_rhs_per_outer_max": cfg.max_trial_steps * rhs_per_forward, "final_forward_rhs": rhs_per_forward, "upper_bound_rhs": cfg.max_iter * maximum_outer_rhs + rhs_per_forward}, "whole_campaign_upper_bound_rhs": int(objects * noises * len(methods) * (cfg.max_iter * maximum_outer_rhs + rhs_per_forward)), "actual_cost_rule": "record accepted/rejected trials, converged or capped termination, forward_rhs and tangent_rhs separately for every paired case; never substitute this upper bound for observed cost", "guard": "do not launch until root freezes protocol, data quadrature, material families, split, stopping/cap and independent N128 validation."}


def main() -> None:
    parser = argparse.ArgumentParser(); parser.add_argument("--max-iter", type=int, default=5); parser.add_argument("--prepare-frozen", action="store_true"); args = parser.parse_args()
    out = ROOT / "runs/a2/som"; out.mkdir(parents=True, exist_ok=True)
    if args.prepare_frozen:
        (out / "frozen_campaign_call_ledger.json").write_text(json.dumps(frozen_campaign_call_ledger(), indent=2))
        print(json.dumps(frozen_campaign_call_ledger(), indent=2)); return
    # Compact CPU development geometry; no held-out choices are tuned after observing its result.
    geom = Geometry(n=32, n_tx=6, n_rx=24, aperture="half")
    vies = [VIE(geom, f) for f in FREQUENCIES_HZ]
    truth = np.array([.58, -.042, -.015, np.log(.031), np.log(.024), .18, .46, .046, .018, np.log(.026), np.log(.036), -.31])
    theta0 = np.array([.42, -.030, -.004, np.log(.041), np.log(.031), .05, .34, .030, .007, np.log(.036), np.log(.046), -.10])
    chis, clean, _ = solve_stack(truth, vies, rtol=2e-9)
    rng = np.random.default_rng(20260911); all_sca = np.concatenate([x["scattered"].ravel() for x in clean]); noise_sigma = .01 * np.linalg.norm(all_sca) / np.sqrt(all_sca.size)
    observed = [x["scattered"] + noise_sigma / np.sqrt(2) * (rng.normal(size=x["scattered"].shape) + 1j*rng.normal(size=x["scattered"].shape)) for x in clean]
    train_rx = np.arange(geom.n_rx) % 2 == 0; held_rx = ~train_rx
    cfg = FitConfig(max_iter=args.max_iter)
    diagnostic = selection_diagnostics(theta0, vies, train_rx, noise_sigma)
    methods = [fit_parameter_tangent_som(name, theta0, truth, observed, clean, vies, train_rx, held_rx, noise_sigma, cfg) for name in ("ordinary_lm", "singlefold_som", "twofold_tsomg")]
    payload = {"scope": "one deterministic bounded development case; not a statistical comparison, certification, current-space TSOM implementation, or novelty result", "protocol": "Gaussian/delegated/a2_som/METHODS_AND_PROTOCOL.md", "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), "environment": {"python": sys.version, "platform": platform.platform(), "numpy": np.__version__, "scipy": scipy.__version__, "thread_cap": 2}, "geometry": {"n": geom.n, "side_m": geom.side, "n_tx": geom.n_tx, "n_rx": geom.n_rx, "aperture": geom.aperture, "frequencies_hz": FREQUENCIES_HZ}, "parameter_metric_scales": SCALES.tolist(), "noise_sigma_complex_component": float(noise_sigma), "split": {"train_receiver_indices": np.flatnonzero(train_rx).tolist(), "held_receiver_indices": np.flatnonzero(held_rx).tolist()}, "truth": truth.tolist(), "initial": theta0.tolist(), "twofold_definition": "V1 from realified whitened full-wave receiver tangent A; V2 from SVD of B Vd, with B normalized D*K internal coupling; B is selection only", "selector_diagnostics": diagnostic, "explicit_tsom_current_variant": "not implemented: needs separately constructed current-space weak/strong basis and state-residual penalty continuation", "methods": methods}
    (out / "results.json").write_text(json.dumps(payload, indent=2))
    print(json.dumps({m["method"]: {"material_relative_l2": m["material_relative_l2"], "parameter_scaled_error": m["parameter_scaled_error"], "wall_seconds": m["wall_seconds"]} for m in methods}, indent=2))


if __name__ == "__main__":
    main()
