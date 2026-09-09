#!/usr/bin/env python
"""Threshold and negative-control ablations for the A2 replication.

Bounded supplement only.  The script reuses the deterministic fixtures and
linear-algebra primitives from ``a2val`` (E1/E2/E3/E5) and the read-only
shared physical core.  It does not invent new theorems, does not modify any
original package or physics file, and records raw data for every threshold
value (no selection).

Axes
----
A1  physical reduced-state vs declared free-current envelope:
    per physical seed 201-212 the exact map/pose tangent has no free-current
    block (C0 empty); a conservative envelope C_w uses declared left singular
    directions of the realified tangent.  For each envelope width in
    {2,4,8} and rank thresholds in {1e-12,1e-10,1e-8} we record pose Gram,
    visible residual ranks and the envelope-saturation flag.  A zero envelope
    score is recorded as an envelope failure only (G0-D), never as physical
    nonidentifiability.
A2  relative-only vs absolute + bias gates:
    reuse ``e3.b_zero_control`` and ``e3.confounding_sweep`` exact examples
    and sweep a relative map-retention tolerance, an absolute pose-information
    gate, and a residual/bias gate over three values each.
A3  shared-map vs independent-map:
    reuse the complementary-pair and random 3-frame examples; add three
    per-frame pose-noise scales and sweep the Loewner PSD gate and an
    independent-map visibility gate.
A4  greedy vs pair-lookahead:
    reuse the fixed E5 8-candidate fixture, sweeping the regularisation gate,
    pose-noise scale and the greedy-certification gap tolerance (three values
    each).  Raw policy logdets and gaps are stored without pruning.

Usage (from the replication root):
    python scripts/run_ablations.py

Outputs (relative to the working directory):
    results/ablations_results.json, results/ablations_summary.csv,
    figures/ablation_2x2.png (optional compact figure)
"""

from __future__ import annotations

import csv
import itertools
import json
import os
import random as _random
import sys
import time

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

PHYSICS_DIR = (
    "/Volumes/migodam's-external-brain/Research/Inv_SLAM/"
    "research/delegated/a2_physics"
)
if PHYSICS_DIR not in sys.path:
    sys.path.insert(0, PHYSICS_DIR)
import physics  # noqa: E402  (read-only shared core)

from a2val import common, e2, e3, e5, physical_runner as pr


EPS = float(np.finfo(float).eps)
PHYS_SEEDS = list(range(201, 213))


def _sym(X):
    X = np.asarray(X, dtype=float)
    return 0.5 * (X + X.T)


def _fro(X):
    return float(np.linalg.norm(np.asarray(X, dtype=float), "fro"))


def _min_eig(X):
    return float(np.linalg.eigvalsh(_sym(X))[0])


def _eig_asc(X):
    return np.sort(np.linalg.eigvalsh(_sym(X))).tolist()


def _jsonable(obj):
    if isinstance(obj, dict):
        return {k: _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, np.generic):
        return obj.item()
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, complex):
        return {"re": obj.real, "im": obj.imag}
    if isinstance(obj, float) and (np.isnan(obj) or np.isinf(obj)):
        return None
    return obj


def _reference(seed: int):
    """Deterministic reference rule used by the physical runners."""
    rng = np.random.default_rng(int(seed))
    alpha = 0.2 + 0.05 * rng.standard_normal(9)
    alpha = np.clip(alpha, -0.5, 0.5)
    x = 0.02 * rng.standard_normal(3)
    x = np.clip(x, -0.05, 0.05)
    return alpha, x, rng


def _physical_tangent(seed: int):
    """N=8 full-aperture realified map/pose tangents at the declared seed."""
    cfg = physics.Config(N=8, aperture="full")
    model = physics.Model(cfg)
    alpha, x, _rng = _reference(seed)
    out = model.forward(alpha, x, jacobian=True)
    A_r = physics.realify_jacobian(out["A"])
    B_r = physics.realify_jacobian(out["B"])
    state_res = [st["state_residual_rel"] for st in out["states"]]
    return {
        "seed": int(seed),
        "alpha": alpha,
        "x": x,
        "A_r": A_r,
        "B_r": B_r,
        "out": out,
        "state_residual_rel": state_res,
        "work": out["work"],
    }


# ---------------------------------------------------------------------------
# A1 physical vs free-current envelope
# ---------------------------------------------------------------------------


def _rank_counts(sv, thresholds):
    """Absolute and relative numerical-rank counts at three thresholds."""
    sv = np.asarray(sv, dtype=float)
    sv = sv[sv > 0.0]
    out = []
    for t in thresholds:
        if sv.size == 0:
            out.append({"threshold": t, "rank_rel": 0, "rank_abs": 0})
            continue
        out.append({
            "threshold": t,
            "rank_rel": int(np.count_nonzero(sv > t * sv[0])),
            "rank_abs": int(np.count_nonzero(sv > t)),
        })
    return out


def run_axis_a1():
    widths = (2, 4, 8)
    rank_ratios = (1e-12, 1e-10, 1e-8)
    support_tols = rank_ratios
    rows = []
    for seed in PHYS_SEEDS:
        pt = _physical_tangent(seed)
        A, B = pt["A_r"], pt["B_r"]
        T = np.hstack([A, B])
        U, sT, _ = np.linalg.svd(T, full_matrices=False)
        C0 = np.zeros((A.shape[0], 0), dtype=float)
        spaces = [("physical_r_free_0", C0, 0)]
        for w in widths:
            spaces.append(("envelope_C%d" % w, U[:, :w], int(w)))
        for label, C, w in spaces:
            Jx = e2.Jx(A, B, C)
            Bv = e2.Bv(A, B, C)
            sv = np.linalg.svd(Bv, compute_uv=False)
            abs_tol = pr.backward_rank_tol(B)
            for rt in rank_ratios:
                rho_rec = e2.rho_values(A, B, C, support_tol_rel=rt)
                rank_rec = _rank_counts(sv, (rt, rt, abs_tol))
                row = {
                    "axis": "A1",
                    "seed": int(seed),
                    "space": label,
                    "envelope_width": w,
                    "support_tol_rel_rho": rt,
                    "tr_Jx": float(np.trace(Jx)),
                    "min_eig_Jx": _min_eig(Jx),
                    "Jx_eig_asc": _eig_asc(Jx),
                    "Jx_fro": _fro(Jx),
                    "Bv_sv_abs": np.sort(sv)[::-1].tolist(),
                    "rank_rel_at_ratio": rank_rec[0],
                    "rank_abs_at_same_absolute_value": rank_rec[1],
                    "rank_abs_at_backward_tol": rank_rec[2],
                    "backward_rank_tol": abs_tol,
                    "map_support_dim": int(rho_rec["support_dim"]),
                    "map_rho_values": rho_rec["rho"],
                    "k0_singular_values": rho_rec["k0_singular_values"],
                    "envelope_saturated": bool(
                        _min_eig(Jx) <= max(1e-12, 2.0 * abs_tol ** 2)
                    ),
                    "state_residual_rel_max": float(np.max(
                        pt["state_residual_rel"]
                    )),
                    "n_physical_rows": int(A.shape[0]),
                    "n_pose_cols": int(B.shape[1]),
                    "n_map_cols": int(A.shape[1]),
                    "work_rhs_solves": int(pt["work"]["rhs_solves_total"]),
                }
                rows.append(row)
    # Envelope-width T5 controls at each declared width (one representative
    # threshold record per seed/adjacent pair).
    t5_rows = []
    for seed in PHYS_SEEDS:
        pt = _physical_tangent(seed)
        A, B = pt["A_r"], pt["B_r"]
        T = np.hstack([A, B])
        U, _sT, _ = np.linalg.svd(T, full_matrices=False)
        for lo, hi, name in (
            (0, 2, "C0->C2"),
            (2, 4, "C2->C4"),
            (4, 8, "C4->C8"),
        ):
            Cl = np.zeros((A.shape[0], 0)) if lo == 0 else U[:, :lo]
            Ch = U[:, :hi]
            rec = pr.t5_adjacent_pair(A, B, Cl, Ch, name)
            t5_rows.append({
                "axis": "A1_t5",
                "seed": int(seed),
                "pair": name,
                "declared_E_dim": rec["declared_dims"]["E_dim"],
                "loss_residual_fro": rec["loss_residual_fro"],
                "backward_scaled_residual": rec["backward_scaled_residual"],
                "rank_drop_actual": rec["rank_drop_actual"],
                "rank_drop_predicted": rec["rank_drop_predicted"],
                "t5b_matches": rec["t5b_matches"],
            })
    summary = {
        "physical_seeds": PHYS_SEEDS,
        "envelope_widths": list(widths),
        "rank_ratio_thresholds": list(rank_ratios),
        "support_tols": list(support_tols),
        "n_rows": len(rows),
        "n_t5_rows": len(t5_rows),
        "physical_info_positive": {
            "count": sum(
                1 for r in rows
                if r["space"] == "physical_r_free_0" and r["min_eig_Jx"] > 0.0
            ),
            "min_eig_Jx_over_seeds": min(
                r["min_eig_Jx"] for r in rows
                if r["space"] == "physical_r_free_0"
            ),
        },
        "envelope_saturation_examples": [
            {"seed": r["seed"], "width": r["envelope_width"],
             "min_eig_Jx": r["min_eig_Jx"], "space": r["space"]}
            for r in rows
            if r["envelope_saturated"]
        ],
        "note": (
            "Zero envelope visibility is an envelope diagnostic only; it is "
            "never labelled physical nonidentifiability (G0-D)."
        ),
    }
    return {"rows": rows, "t5_rows": t5_rows, "summary": summary}


# ---------------------------------------------------------------------------
# A2 relative-only vs absolute + bias gates
# ---------------------------------------------------------------------------


def run_axis_a2():
    rel_tols = (1e-10, 1e-8, 1e-6)
    abs_gates = (1e-12, 1e-8, 1e-4)
    res_gates = (1e-6, 1e-3, 1e-1)
    bzero = e3.b_zero_control()
    sweep = e3.confounding_sweep()
    slb = sweep["small_residual_large_bias"]
    conf_rows = sweep["records"]
    rows = []
    # b_zero fixture: full relative map retention, zero absolute pose info.
    map_spectrum = bzero["map_spectrum"]
    jx_sv = bzero["J_x_singular_values"]
    min_rho = float(min(map_spectrum)) if map_spectrum else 0.0
    for rt in rel_tols:
        rel_pass = bool(min_rho >= 1.0 - rt)
        for ag in abs_gates:
            abs_pass = bool(jx_sv and float(min(jx_sv)) > ag)
            rows.append({
                "axis": "A2",
                "fixture": "b_zero_control",
                "map_retention_min_rho": min_rho,
                "relative_gate_tol": rt,
                "relative_only_gate_pass": rel_pass,
                "absolute_pose_gate_threshold": ag,
                "absolute_pose_gate_pass": abs_pass,
                "pose_Jx_min_singular_value": (
                    float(min(jx_sv)) if jx_sv else None
                ),
                "false_assurance_relative_only": bool(rel_pass and not abs_pass),
            })
    # Confounding eps sweep: relative/absolute/bias-risk comparison.
    for cr in conf_rows:
        for rt in rel_tols:
            # relative map retention is 1 by construction in this 2-obs model;
            # the misleading gate is the residual scale, so relative gates are
            # not informative here.
            for rg in res_gates:
                rows.append({
                    "axis": "A2",
                    "fixture": "confounding_eps",
                    "eps_param": cr["eps_param"],
                    "full_ls_risk": cr["full_ls_risk"],
                    "truncated_risk": cr["truncated_risk"],
                    "truncation_better_analytic": cr["truncation_better"],
                    "relative_map_gate_tol": rt,
                    "residual_gate_threshold": rg,
                    "relative_or_residual_gate_not_computed": True,
                    "bias_aware_criterion": cr["criterion"],
                })
    # Small-residual/large-bias fixture: residual-only gate vs bias-aware.
    res_rel = slb["residual_y2_relative_to_pose_error"]
    for rg in res_gates:
        res_pass = bool(res_rel <= rg)
        for ag in abs_gates:
            # pose error norm 1.0; absolute error gate is meaningful only if
            # the declared precision permits error <= threshold.
            rows.append({
                "axis": "A2",
                "fixture": "small_residual_large_bias",
                "residual_relative_to_pose_error": res_rel,
                "pose_bias_c_star": slb["truncation_bias"],
                "residual_gate_threshold": rg,
                "residual_only_gate_pass": res_pass,
                "absolute_error_gate_threshold": ag,
                "absolute_error_gate_pass": bool(slb["truncation_bias"] <= ag),
                "false_assurance_residual_only": bool(
                    res_pass and slb["truncation_bias"] > ag
                ),
            })
    summary = {
        "relative_only_tols": list(rel_tols),
        "absolute_gates": list(abs_gates),
        "residual_gates": list(res_gates),
        "relative_only_false_assurance_count": int(sum(
            1 for r in rows
            if r.get("false_assurance_relative_only")
        )),
        "residual_only_false_assurance_count": int(sum(
            1 for r in rows
            if r.get("false_assurance_residual_only")
        )),
        "note": (
            "Relative-only and residual-only gates pass by construction on "
            "b_zero and small-residual-large-bias when their thresholds are "
            "loose; absolute/bias-aware gates are the intended controls."
        ),
    }
    return {"rows": rows, "summary": summary}


# ---------------------------------------------------------------------------
# A3 shared-map vs independent-map
# ---------------------------------------------------------------------------


def _perturb_frames(frames, noise_scale, seed):
    """Add deterministic per-frame pose noise (raw rows are recorded)."""
    rng = np.random.default_rng(int(seed))
    out = []
    for a, b in frames:
        a = np.asarray(a, dtype=float)
        b = np.asarray(b, dtype=float)
        if noise_scale > 0.0:
            b = b + noise_scale * rng.standard_normal(b.shape)
        out.append((a.copy(), b.copy()))
    return out


def _frame_stats(frames):
    Jper = e5._per_frame_infos(frames)
    A, B = e5._stack_AB(frames)
    Jstack = e5._pose_info(A, B)
    Jsum = sum(Jper)
    gap = _sym(Jstack - Jsum)
    return {
        "J_per_frame_fro": [_fro(J) for J in Jper],
        "J_per_frame_min_eig": [_min_eig(J) for J in Jper],
        "sum_J_min_eig": _min_eig(Jsum),
        "sum_J_fro": _fro(Jsum),
        "J_stack_min_eig": _min_eig(Jstack),
        "J_stack_fro": _fro(Jstack),
        "min_eig_stack_minus_sum": _min_eig(gap),
        "fro_gap": _fro(gap),
        "J_stack_eig": _eig_asc(Jstack),
        "J_sum_eig": _eig_asc(Jsum),
    }


def run_axis_a3():
    noise_scales = (0.0, 0.01, 0.1)
    psd_gates = (1e-12, 1e-8, 1e-6)
    vis_gates = (1e-12, 1e-8, 1e-4)
    rows = []

    # Complementary pair (deterministic).
    comp = [
        (np.array([[1.0]]), np.array([[1.0]])),
        (np.array([[1.0]]), np.array([[-1.0]])),
    ]
    for nscale in noise_scales:
        fr = _perturb_frames(comp, nscale, seed=41_000 + int(nscale * 1000))
        st = _frame_stats(fr)
        for pg in psd_gates:
            for vg in vis_gates:
                rows.append({
                    "axis": "A3",
                    "fixture": "complementary_pair",
                    "noise_scale": nscale,
                    "psd_gate": pg,
                    "shared_map_stack_ge_sum_gate": bool(
                        st["min_eig_stack_minus_sum"] >= -pg
                    ),
                    "independent_map_visible_gate": bool(
                        st["sum_J_min_eig"] > vg
                    ),
                    "false_assurance_independent_map": bool(
                        st["J_stack_min_eig"] > vg
                        and st["sum_J_min_eig"] <= vg
                    ),
                    **{k: v for k, v in st.items()},
                })

    # Random 3-frame (seeded like e5.stack_vs_sum, 3/4/5 rows).
    base = e5._random_frames(2750, (3, 4, 5), 2, 3)
    for nscale in noise_scales:
        fr = _perturb_frames(base, nscale, seed=41_300 + int(nscale * 1000))
        st = _frame_stats(fr)
        for pg in psd_gates:
            for vg in vis_gates:
                rows.append({
                    "axis": "A3",
                    "fixture": "random_three_frame",
                    "noise_scale": nscale,
                    "psd_gate": pg,
                    "shared_map_stack_ge_sum_gate": bool(
                        st["min_eig_stack_minus_sum"] >= -pg
                    ),
                    "independent_map_visible_gate": bool(
                        st["sum_J_min_eig"] > vg
                    ),
                    "false_assurance_independent_map": bool(
                        st["J_stack_min_eig"] > vg
                        and st["sum_J_min_eig"] <= vg
                    ),
                    **{k: v for k, v in st.items()},
                })

    # Hidden common-compensator example (three frames with identical map
    # directions): independent map looks visible only under perturbation.
    c_hidden = [np.array([[0.7], [-1.3]])] * 3

    def _hidden_frames(rng0):
        frames = []
        for _c in c_hidden:
            a = rng0.standard_normal((3, 2))
            frames.append((a, a @ _c))
        return frames

    rng0 = np.random.default_rng(2700)
    hidden_base = _hidden_frames(rng0)
    for nscale in noise_scales:
        rng1 = np.random.default_rng(27_000 + int(nscale * 1000))
        fr = []
        for a, b in hidden_base:
            if nscale > 0.0:
                b = b + nscale * rng1.standard_normal(b.shape)
            fr.append((a.copy(), b.copy()))
        st = _frame_stats(fr)
        for pg in psd_gates:
            for vg in vis_gates:
                rows.append({
                    "axis": "A3",
                    "fixture": "hidden_common_compensator",
                    "noise_scale": nscale,
                    "psd_gate": pg,
                    "shared_map_stack_ge_sum_gate": bool(
                        st["min_eig_stack_minus_sum"] >= -pg
                    ),
                    "independent_map_visible_gate": bool(
                        st["sum_J_min_eig"] > vg
                    ),
                    "false_assurance_independent_map": bool(
                        st["J_stack_min_eig"] > vg
                        and st["sum_J_min_eig"] <= vg
                    ),
                    **{k: v for k, v in st.items()},
                })

    exact = {
        "minimal_pair": e5.minimal_pair(),
        "compensation_criterion": e5.compensation_criterion(),
        "stack_vs_sum": e5.stack_vs_sum(),
    }
    summary = {
        "noise_scales": list(noise_scales),
        "psd_gates": list(psd_gates),
        "visibility_gates": list(vis_gates),
        "exact_examples_preserved": {
            k: {"pass": bool(exact[k].get("pass"))} for k in exact
        },
        "complementary_stack_minus_sum_gate_pass_counts": {
            str(g): sum(
                1 for r in rows
                if r["fixture"] == "complementary_pair"
                and r["psd_gate"] == g
                and r["shared_map_stack_ge_sum_gate"]
            )
            for g in psd_gates
        },
        "independent_map_false_assurance_count": int(sum(
            1 for r in rows if r.get("false_assurance_independent_map")
        )),
        "note": (
            "The shared-map stack (not the per-frame sum) is the exact "
            "criterion; independent-map visibility is the mismatch control."
        ),
    }
    return {"rows": rows, "exact": exact, "summary": summary}


# ---------------------------------------------------------------------------
# A4 greedy vs pair-lookahead
# ---------------------------------------------------------------------------


def _fixture_logdet(A0, B0, frames, subset, Lam):
    if subset:
        A = np.vstack([A0] + [frames[i][0] for i in subset])
        B = np.vstack([B0] + [frames[i][1] for i in subset])
    else:
        A, B = A0, B0
    J = e5._pose_info(A, B)
    return e5._logdet(J + Lam)


def _one_policy_run(A0, B0, frames, Lam, noise_scale, seed_base, row_idx):
    rng = np.random.default_rng(seed_base + row_idx)
    B0n = B0.copy()
    if noise_scale > 0.0:
        B0n = B0 + noise_scale * rng.standard_normal(B0.shape)
    fn = []
    for a, b in frames:
        bn = b.copy()
        if noise_scale > 0.0:
            bn = bn + noise_scale * rng.standard_normal(bn.shape)
        fn.append((a.copy(), bn.copy()))
    frames = fn
    K = 3
    n = 8
    # Greedy
    greedy_set = []
    remaining = set(range(n))
    greedy_evals = 0
    while len(greedy_set) < K:
        best = None
        best_val = None
        for i in remaining:
            greedy_evals += 1
            val = _fixture_logdet(A0, B0n, frames, greedy_set + [i], Lam)
            if best is None or val > best_val + 1e-15:
                best = i
                best_val = val
        greedy_set.append(best)
        remaining.remove(best)
    # Pair look-ahead (best pair + one greedy finish)
    pair_evals = 0
    best_pair = None
    best_pair_val = None
    for pair in itertools.combinations(range(n), 2):
        pair_evals += 1
        val = _fixture_logdet(A0, B0n, frames, list(pair), Lam)
        if best_pair is None or val > best_pair_val + 1e-15:
            best_pair = pair
            best_pair_val = val
    best_pair = list(best_pair)
    finish_evals = 0
    finish = None
    finish_val = None
    for i in range(n):
        if i in best_pair:
            continue
        finish_evals += 1
        val = _fixture_logdet(A0, B0n, frames, best_pair + [i], Lam)
        if finish is None or val > finish_val + 1e-15:
            finish = i
            finish_val = val
    pair_set = best_pair + [finish]
    # Exhaustive
    ex_evals = 0
    best_comb = None
    best_ex = None
    for comb in itertools.combinations(range(n), K):
        ex_evals += 1
        val = _fixture_logdet(A0, B0n, frames, list(comb), Lam)
        if best_ex is None or val > best_ex + 1e-15:
            best_ex = val
            best_comb = list(comb)
    # Random selection, same 12 seeds as the synthetic policy control
    random_rows = []
    rnd_evals = 0
    for seed_r in range(401, 413):
        rnd = _random.Random(seed_r)
        comb = tuple(sorted(rnd.sample(range(n), K)))
        rnd_evals += 1
        val = _fixture_logdet(A0, B0n, frames, list(comb), Lam)
        random_rows.append({"seed": seed_r, "subset": list(comb), "logdet": val})
    random_best = max(random_rows, key=lambda t: t["logdet"])
    greedy_val = _fixture_logdet(A0, B0n, frames, greedy_set, Lam)
    pair_val = _fixture_logdet(A0, B0n, frames, pair_set, Lam)
    return {
        "greedy": {"subset": greedy_set, "logdet": greedy_val},
        "pair_lookahead": {"subset": pair_set, "logdet": pair_val},
        "exhaustive": {"best_subset": best_comb, "best_logdet": best_ex},
        "random_best_of_12": random_best,
        "random_rows": random_rows,
        "evaluation_counts": {
            "greedy": greedy_evals,
            "pair_lookahead": pair_evals + finish_evals,
            "exhaustive": ex_evals,
            "random_best_of_12": rnd_evals,
        },
    }


def run_axis_a4():
    gates = (0.01, 0.1, 1.0)
    noise_scales = (0.0, 0.02, 0.2)
    cert_tols = (1e-9, 1e-6, 1e-3)
    A0, B0, frames = e5._policy_fixture(seed=3032)
    p = B0.shape[1]
    rows = []
    row_idx = 0
    for lam in gates:
        Lam = lam * np.eye(p)
        for noise in noise_scales:
            run = _one_policy_run(
                A0, B0, frames, Lam, noise, seed_base=60_000, row_idx=row_idx
            )
            row_idx += 1
            greedy_gap = float(run["exhaustive"]["best_logdet"]
                               - run["greedy"]["logdet"])
            pair_gap = float(run["exhaustive"]["best_logdet"]
                             - run["pair_lookahead"]["logdet"])
            for ct in cert_tols:
                rows.append({
                    "axis": "A4",
                    "lambda_x": lam,
                    "noise_scale": noise,
                    "certification_gap_tol": ct,
                    "greedy_logdet": float(run["greedy"]["logdet"]),
                    "pair_logdet": float(run["pair_lookahead"]["logdet"]),
                    "exhaustive_logdet": float(run["exhaustive"]["best_logdet"]),
                    "random_best_logdet": float(
                        run["random_best_of_12"]["logdet"]
                    ),
                    "greedy_subset": run["greedy"]["subset"],
                    "pair_subset": run["pair_lookahead"]["subset"],
                    "exhaustive_subset": run["exhaustive"]["best_subset"],
                    "random_best_subset": run["random_best_of_12"]["subset"],
                    "greedy_optimality_gap": greedy_gap,
                    "pair_optimality_gap": pair_gap,
                    "greedy_certified_optimal": bool(greedy_gap <= ct),
                    "pair_certified_optimal": bool(pair_gap <= ct),
                    "evaluation_counts": run["evaluation_counts"],
                    "raw_random_rows": run["random_rows"],
                })
    summary = {
        "lambda_x_values": list(gates),
        "noise_scales": list(noise_scales),
        "certification_tols": list(cert_tols),
        "greedy_gap_min": min(r["greedy_optimality_gap"] for r in rows),
        "greedy_gap_max": max(r["greedy_optimality_gap"] for r in rows),
        "exhaustive_strictly_better_count": sum(
            1 for r in rows if r["exhaustive_logdet"]
            > r["greedy_logdet"] + 1e-12
        ),
        "pair_recovered_exhaustive_count": sum(
            1 for r in rows if r["pair_logdet"]
            >= r["exhaustive_logdet"] - 1e-9
        ),
        "note": (
            "Raw policy logdets are recorded for every (lambda, noise) "
            "combination; certification tolerances only label the gap, they "
            "do not change the algorithm or the objective."
        ),
    }
    return {"rows": rows, "summary": summary}


# ---------------------------------------------------------------------------
# Figure and output
# ---------------------------------------------------------------------------


def _compact_figure(a1, a2, a3, a4):
    """2x2 compact summary; purely descriptive, no fitted model."""
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.6))
    # A1 panel
    ax = axes[0, 0]
    phys_min = min(r["min_eig_Jx"] for r in a1["rows"]
                   if r["space"] == "physical_r_free_0")
    xw = sorted({r["envelope_width"] for r in a1["rows"]
                 if r["envelope_width"]})
    env_min = {w: min(r["min_eig_Jx"] for r in a1["rows"]
                      if r["envelope_width"] == w) for w in xw}
    ax.axhline(phys_min, color="#1b6ca8", ls="-", lw=2,
               label="physical C0 min eig(Jx)")
    for w in xw:
        ax.plot([w], [env_min[w]], "o", color="#b03a2e", ms=7,
                label="envelope C%d min eig(Jx)" % w)
    ax.set_yscale("log")
    ax.set_xticks(list(xw))
    ax.set_xlabel("declared envelope width")
    ax.set_ylabel("min eig(Jx)")
    ax.set_title("A1 physical vs envelope\n(seeds 201-212, worst min eig)")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7)
    # A2 panel
    ax = axes[0, 1]
    rf = [r for r in a2["rows"] if r["fixture"] == "b_zero_control"]
    x = [r["relative_gate_tol"] for r in rf]
    y = [1.0 if r["relative_only_gate_pass"] else 0.0 for r in rf]
    ax.plot(x, y, "o-", color="#1b6ca8", ms=5, label="relative-only gate")
    y2 = [1.0 if r["absolute_pose_gate_pass"] else 0.0 for r in rf]
    ax.plot(x, y2, "s-", color="#b03a2e", ms=5, label="absolute Jx gate")
    ax.set_xscale("log")
    ax.set_ylim(-0.1, 1.1)
    ax.set_yticks([0, 1])
    ax.set_xlabel("relative retention tolerance")
    ax.set_title("A2 b_zero gates\n(relative passes, absolute fails)")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7)
    # A3 panel
    ax = axes[1, 0]
    for fixture, col in (("complementary_pair", "#1b6ca8"),
                         ("random_three_frame", "#b03a2e"),
                         ("hidden_common_compensator", "#3c8c40")):
        xs = sorted({r["noise_scale"] for r in a3["rows"]
                     if r["fixture"] == fixture})
        ys = [min(r["min_eig_stack_minus_sum"] for r in a3["rows"]
                  if r["fixture"] == fixture and r["noise_scale"] == xx)
              for xx in xs]
        ax.plot(xs, ys, "o-", ms=5, color=col, label=fixture)
    ax.axhline(0.0, color="0.5", lw=0.8)
    ax.set_xscale("log")
    ax.set_xlabel("per-frame pose noise scale")
    ax.set_title("A3 min eig(J_stack - sum J_l)")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7)
    # A4 panel
    ax = axes[1, 1]
    for lam in sorted({r["lambda_x"] for r in a4["rows"]}):
        rows_lam = [r for r in a4["rows"] if r["lambda_x"] == lam]
        xs = sorted({r["noise_scale"] for r in rows_lam})
        ys = [next(r["greedy_optimality_gap"] for r in rows_lam
                   if r["noise_scale"] == xx and r["certification_gap_tol"]
                   == 1e-9) for xx in xs]
        ax.plot(xs, ys, "o-", ms=5, label="lambda=%.2g" % lam)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("pose-noise scale")
    ax.set_ylabel("greedy-exhaustive gap")
    ax.set_title("A4 greedy suboptimality\n(8 candidates, K=3)")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7)
    fig.suptitle("A2 replication negative controls (threshold sweep)", y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    os.makedirs("figures", exist_ok=True)
    fig.savefig("figures/ablation_2x2.png", dpi=150)
    plt.close(fig)


def _csv_dump(records, path):
    if not records:
        return
    keys = []
    for rec in records:
        for k in rec:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=keys,
                                extrasaction="ignore")
        writer.writeheader()
        for rec in records:
            writer.writerow({k: _jsonable(v) for k, v in rec.items()})


def main():
    t_start = time.perf_counter()
    os.makedirs("results", exist_ok=True)
    os.makedirs("figures", exist_ok=True)
    print("A1 physical vs free-current envelope ...")
    a1 = run_axis_a1()
    print("A2 relative-only vs absolute+bias gates ...")
    a2 = run_axis_a2()
    print("A3 shared-map vs independent-map ...")
    a3 = run_axis_a3()
    print("A4 greedy vs pair-lookahead ...")
    a4 = run_axis_a4()
    doc = {
        "title": "A2 replication threshold and negative-control ablations",
        "date": "2026-09-06",
        "package": "a2val (replication copy) + shared physics core",
        "environment": {
            "interpreter": sys.executable,
            "physics_module": physics.__file__,
        },
        "settings": {
            "A1": {
                "physical_seeds": PHYS_SEEDS,
                "envelope_widths": [2, 4, 8],
                "rank_thresholds": [1e-12, 1e-10, 1e-8],
                "space_label": (
                    "physical_r_free_0 = exact reduced-state map/pose tangent "
                    "with no free-current block; envelope_Cw = declared left "
                    "singular directions of the realified data tangent"
                ),
            },
            "A2": {
                "fixtures": ["b_zero_control", "confounding_eps",
                             "small_residual_large_bias"],
                "relative_tols": [1e-10, 1e-8, 1e-6],
                "absolute_gates": [1e-12, 1e-8, 1e-4],
                "residual_gates": [1e-6, 1e-3, 1e-1],
            },
            "A3": {
                "fixtures": ["complementary_pair", "random_three_frame",
                             "hidden_common_compensator"],
                "noise_scales": [0.0, 0.01, 0.1],
                "psd_gates": [1e-12, 1e-8, 1e-6],
                "visibility_gates": [1e-12, 1e-8, 1e-4],
            },
            "A4": {
                "fixture_seed": 3032,
                "n_candidates": 8,
                "K": 3,
                "lambda_x_values": [0.01, 0.1, 1.0],
                "noise_scales": [0.0, 0.02, 0.2],
                "certification_tols": [1e-9, 1e-6, 1e-3],
                "random_seeds": list(range(401, 413)),
            },
        },
        "axes": {
            "A1": a1,
            "A2": a2,
            "A3": a3,
            "A4": a4,
        },
        "runtime_wall_seconds": time.perf_counter() - t_start,
    }
    with open("results/ablations_results.json", "w") as fh:
        json.dump(_jsonable(doc), fh, indent=2)
        fh.write("\n")
    flat = a1["rows"] + a1["t5_rows"] + a2["rows"] + a3["rows"] + a4["rows"]
    _csv_dump(flat, "results/ablations_summary.csv")
    _compact_figure(a1, a2, a3, a4)
    print("Saved: results/ablations_results.json, "
          "results/ablations_summary.csv, figures/ablation_2x2.png")
    print("Row counts:", {
        "A1": len(a1["rows"]),
        "A1_t5": len(a1["t5_rows"]),
        "A2": len(a2["rows"]),
        "A3": len(a3["rows"]),
        "A4": len(a4["rows"]),
    })
    print("A1 physical positive over seeds:",
          a1["summary"]["physical_info_positive"])
    print("A4 greedy gap range: [%g, %g]" % (
        a4["summary"]["greedy_gap_min"], a4["summary"]["greedy_gap_max"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
