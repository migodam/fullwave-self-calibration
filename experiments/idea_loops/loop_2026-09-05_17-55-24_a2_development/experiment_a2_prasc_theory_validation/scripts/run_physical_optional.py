#!/usr/bin/env python
"""OPTIONAL second-stage physical-tangent runner for E2/E3/E5 (N=8 smoke).

Usage (from this experiment directory, exact interpreter required):

    python scripts/run_physical_optional.py

The shared physical core is imported read-only from
``research/delegated/a2_physics``; nothing under that directory is modified.
Every result is labelled a physical tangent smoke with declared (no-oracle)
free-current cutoffs; no physical rank, envelope nonidentifiability, nonlinear
basin or method-superiority claim is made.
"""

from __future__ import annotations

import csv
import json
import os
import sys
import time
import traceback

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

PHYSICS_DIR = (
    "/Volumes/migodam's-external-brain/Research/Inv_SLAM/"
    "research/delegated/a2_physics"
)
if PHYSICS_DIR not in sys.path:
    sys.path.insert(0, PHYSICS_DIR)
import physics  # noqa: E402  (read-only shared core)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from a2val import physical_runner as pr  # noqa: E402

EPS = float(np.finfo(float).eps)


def _jsonable(obj):
    if isinstance(obj, dict):
        return {k: _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, np.generic):
        return obj.item()
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, float) and (np.isnan(obj) or np.isinf(obj)):
        return None
    return obj


def _fmt(v, nd=6):
    if v is None:
        return ""
    return f"{v:.{nd}g}"


def _reference(seed: int):
    """Deterministic non-oracle reference, clamped as specified."""
    rng = np.random.default_rng(seed)
    alpha = 0.2 + 0.05 * rng.standard_normal(9)
    alpha = np.clip(alpha, -0.5, 0.5)
    x = 0.02 * rng.standard_normal(3)
    x = np.clip(x, -0.05, 0.05)
    return alpha, x, rng


def _forward_realified(model, seed: int):
    alpha, x, rng = _reference(seed)
    out = model.forward(alpha, x, jacobian=True)
    A_r = physics.realify_jacobian(out["A"])
    B_r = physics.realify_jacobian(out["B"])
    return alpha, x, rng, out, A_r, B_r


def _run_phase0(model) -> dict:
    rec = {
        "status": "ok",
        "label": pr.declared_note(),
        "config": {"N": 8, "aperture": "full"},
        "import_physics": {"module": physics.__name__, "file": physics.__file__},
    }
    try:
        alpha0 = np.full(9, 0.2)
        x0 = np.zeros(3)
        t0 = time.perf_counter()
        out = model.forward(alpha0, x0, jacobian=True)
        wall = time.perf_counter() - t0
        A_r = physics.realify_jacobian(out["A"])
        B_r = physics.realify_jacobian(out["B"])
        rec["forward"] = {
            "baseline_alpha": alpha0.tolist(),
            "baseline_x": x0.tolist(),
            "A_complex_shape": list(out["A"].shape),
            "B_complex_shape": list(out["B"].shape),
            "A_realified_shape": list(A_r.shape),
            "B_realified_shape": list(B_r.shape),
            "wall_seconds": wall,
            "work": {
                k: out["work"][k]
                for k in (
                    "wall_seconds",
                    "factorizations",
                    "factorization_seconds",
                    "rhs_solves_total",
                    "operator_products",
                )
            },
        }
        freq_idx = list(out["blocks"]["freq_idx"])
        rec["row_structure"] = {
            "stable_order": "frequency -> pose -> illumination -> receiver",
            "n_complex_rows": int(out["A"].shape[0]),
            "n_realified_rows": int(A_r.shape[0]),
            "frequencies": sorted(set(freq_idx)),
            "rows_per_frequency": {
                int(fi): sum(
                    1
                    for f, rs, re in zip(
                        out["blocks"]["freq_idx"],
                        out["blocks"]["row_start"],
                        out["blocks"]["row_stop"],
                    )
                    if f == fi
                    for _ in range(int(rs), int(re))
                )
                for fi in sorted(set(freq_idx))
            },
            "blocks_checked": {
                "phase0_shapes_ok": bool(
                    A_r.shape == (576, 9) and B_r.shape == (576, 3)
                ),
            },
        }
    except Exception as exc:  # noqa: BLE001 - record instead of crashing
        rec["status"] = "failed"
        rec["error"] = {"type": type(exc).__name__, "message": str(exc)}
        rec["error_traceback"] = traceback.format_exc()
    return rec


def _run_phase1(model, representative_seed: int = 201) -> dict:
    """E2/E3 physical tangent checks on seeds 201-212 (no oracle rank)."""
    seeds = list(range(201, 213))
    records = []
    errors = []
    rep = None
    for seed in seeds:
        try:
            alpha, x, rng, out, A_r, B_r = _forward_realified(model, seed)
            T = np.hstack([A_r, B_r])
            U, sT, _ = np.linalg.svd(T, full_matrices=False)
            C0 = np.zeros((A_r.shape[0], 0), dtype=float)
            C1 = U[:, :4]
            C2 = U[:, :8]
            pairs = [
                ("C0->C1", C0, C1),
                ("C1->C2", C1, C2),
            ]
            t5 = [
                pr.t5_adjacent_pair(A_r, B_r, Clo, Chi, label)
                for label, Clo, Chi in pairs
            ]
            e3_thm3 = pr.theorem3_physical(A_r, B_r, C1)
            prior = pr.prior_sandwich_physical(
                A_r, B_r, C1, 0.1 * np.eye(3)
            )
            thm7 = pr.theorem7_physical(A_r, B_r, C1, rng)
            rec = {
                "seed": int(seed),
                "alpha": alpha.tolist(),
                "x": x.tolist(),
                "work_wall_seconds": out["work"]["wall_seconds"],
                "declared": {
                    "note": pr.declared_note(),
                    "C1_width": 4,
                    "C2_width": 8,
                    "C_space_source": (
                        "leading left singular vectors of the realified full "
                        "data tangent T=hstack([A_r,B_r]); declared candidate "
                        "current directions, NOT a physical rank claim"
                    ),
                    "T_singular_values": sT.tolist(),
                },
                "t5_pairs": t5,
                "e3_theorem3": e3_thm3,
                "e3_prior": prior,
                "e3_theorem7": thm7,
            }
            if seed == representative_seed:
                rec["representative_tangents"] = {
                    "seed": int(seed),
                    "A_r_shape": list(A_r.shape),
                    "B_r_shape": list(B_r.shape),
                    "A_r": A_r.tolist(),
                    "B_r": B_r.tolist(),
                }
                rep = {
                    "seed": int(seed),
                    "sv_A_r": np.linalg.svd(A_r, compute_uv=False).tolist(),
                    "sv_B_r": np.linalg.svd(B_r, compute_uv=False).tolist(),
                    "sv_T": sT.tolist(),
                    "C1_width": 4,
                    "C2_width": 8,
                }
            records.append(rec)
        except Exception as exc:  # noqa: BLE001
            errors.append({
                "seed": int(seed),
                "type": type(exc).__name__,
                "message": str(exc),
            })

    def _mx(key, default=0.0):
        vals = [r[key] for r in records if key in r]
        return max(vals, default=default)

    def _max_t5(key):
        vals = [
            float(rr[key])
            for r in records
            for rr in r["t5_pairs"]
            if key in rr
        ]
        return max(vals, default=0.0)

    def _min_prior(key):
        vals = [float(r["e3_prior"][key]) for r in records]
        return min(vals, default=0.0)

    thm3_residuals = [
        max(
            r["e3_theorem3"]["map_residual_max_abs"],
            r["e3_theorem3"]["pose_residual_max_abs"],
        )
        for r in records
    ]
    t5b_count = sum(1 for r in records for rr in r["t5_pairs"] if rr["t5b_matches"])
    aggregate = {
        "label": pr.declared_note(),
        "seeds": seeds,
        "n_seeds": len(seeds),
        "n_success": len(records),
        "n_checks_t5": int(len(records) * 2),
        "max_t5a_residual_fro": _max_t5("loss_residual_fro"),
        "max_t5a_backward_scaled_residual": _max_t5(
            "backward_scaled_residual"
        ),
        "t5b_matches_count": int(t5b_count),
        "max_theorem3_residual_max_abs": max(thm3_residuals, default=0.0),
        "min_prior_min_eig_K_eL_minus_K_e": _min_prior(
            "min_eig_K_eL_minus_K_e"
        ),
        "min_prior_min_eig_K0_minus_K_eL": _min_prior(
            "min_eig_K0_minus_K_eL"
        ),
        "prior_sandwich_ok_all": bool(
            records and all(r["e3_prior"]["loewner_sandwich_ok"] for r in records)
        ),
        "theorem7_records": [
            {"seed": r["seed"], **r["e3_theorem7"]} for r in records
        ],
        "theorem7_analytic_evaluated_count": int(
            sum(1 for r in records if r["e3_theorem7"].get("full_rank_B_v"))
        ),
        "theorem7_rank_deficient_skipped_count": int(
            sum(1 for r in records if not r["e3_theorem7"].get("full_rank_B_v"))
        ),
        "theorem7_note": (
            "with declared C1 the residual visible pose matrix B_v = "
            "(I-P_[C1,A])B is numerically zero on every physical seed, so "
            "Theorem 7 (which assumes full column rank of B_v) is not "
            "evaluable; recorded as rank-deficient smoke result, not a "
            "finite-risk claim."
        ),
        "tolerances": {
            "t5a_scaled_pass_threshold": 100.0,
            "theorem3_max_abs_pass_threshold": 1e-6,
            "prior_relative_loewner_tol": 1e-8,
        },
    }
    status = "ok"
    if errors:
        status = "failed" if not records else "partial"
    pass_flag = bool(
        status == "ok"
        and records
        and aggregate["max_t5a_backward_scaled_residual"] < 100.0
        and t5b_count == int(len(records) * 2)
        and aggregate["max_theorem3_residual_max_abs"] <= 1e-6
        and aggregate["prior_sandwich_ok_all"]
    )
    return {
        "status": status,
        "errors": errors,
        "aggregate": aggregate,
        "per_seed_records": records,
        "pass": pass_flag,
        "representative": rep,
    }


def _run_phase2(model) -> dict:
    """E5 shared-map multifrequency physical checks on seeds 401-412."""
    seeds = list(range(401, 413))
    records = []
    errors = []
    for seed in seeds:
        try:
            alpha, x, rng, out, A_r, B_r = _forward_realified(model, seed)
            frames = pr.split_frequency_frames(A_r, B_r, out["blocks"])
            if len(frames) != 4:
                raise RuntimeError(
                    f"expected 4 frequency frames, got {len(frames)}"
                )
            stack = pr.stack_vs_sum_physical(frames)
            seq = pr.sequential_innovation(frames)
            budget = pr.rank_budget_physical(
                frames[0][0],
                frames[0][1],
                frames[1][0],
                frames[1][1],
                seed=int(seed),
                label="physical_freq0_to_freq0+1_C4_to_C8",
            )
            records.append({
                "seed": int(seed),
                "alpha": alpha.tolist(),
                "x": x.tolist(),
                "work_wall_seconds": out["work"]["wall_seconds"],
                "label": pr.declared_note(),
                "frames": stack,
                "stack_vs_sum": stack,
                "sequential_innovation": seq,
                "rank_budget": budget,
            })
        except Exception as exc:  # noqa: BLE001
            errors.append({
                "seed": int(seed),
                "type": type(exc).__name__,
                "message": str(exc),
            })

    stack_psd = [r["stack_vs_sum"] for r in records]
    seq_recs = [
        rec
        for r in records
        for rec in r["sequential_innovation"]["records"]
    ]
    budget_recs = [r["rank_budget"] for r in records]
    innovation_residuals = [
        float(rec["formula_direct_fro_residual"]) for rec in seq_recs
    ]
    innovation_rel = [
        float(rec["formula_direct_rel_residual"]) for rec in seq_recs
    ]
    budget_residuals = [float(rec["identity_fro_residual"]) for rec in budget_recs]
    budget_rel = [float(rec["identity_rel_residual"]) for rec in budget_recs]
    aggregate = {
        "label": pr.declared_note(),
        "seeds": seeds,
        "n_seeds": len(seeds),
        "n_success": len(records),
        "n_sequential_steps": len(seq_recs),
        "n_sequential_singular_steps": int(
            sum(r["sequential_innovation"]["singular_steps"] for r in records)
        ),
        "n_rank_budget_singular_fallbacks": int(
            sum(1 for r in budget_recs if r["singular_fallback"])
        ),
        "max_innovation_residual_fro": max(innovation_residuals, default=0.0),
        "max_innovation_rel_residual": max(innovation_rel, default=0.0),
        "max_budget_identity_residual_fro": max(
            budget_residuals, default=0.0
        ),
        "max_budget_identity_rel_residual": max(budget_rel, default=0.0),
        "stack_ge_sum_count": int(
            sum(1 for s in stack_psd if s["stack_ge_sum"])
        ),
        "rank_budget_loewner_I_ge_L_count": int(
            sum(1 for r in budget_recs if r["loewner_I_ge_L"])
        ),
        "min_eig_J_stack_minus_sum_over_seeds": min(
            (s["min_eig_J_stack_minus_sum"] for s in stack_psd), default=0.0
        ),
        "tolerances": {
            "innovation_rel_pass_threshold": 1e-6,
            "budget_rel_pass_threshold": 1e-6,
            "psd_tolerance_note": (
                "per-check absolute tolerance = "
                "max(1e-10, 20*eps*max(1,||B||_F^2)*n); raw min eigenvalues "
                "are recorded"
            ),
        },
    }
    status = "ok"
    if errors:
        status = "failed" if not records else "partial"
    pass_flag = bool(
        status == "ok"
        and records
        and aggregate["max_innovation_rel_residual"] <= 1e-6
        and aggregate["max_budget_identity_rel_residual"] <= 1e-6
        and aggregate["stack_ge_sum_count"] == len(records)
    )
    return {
        "status": status,
        "errors": errors,
        "aggregate": aggregate,
        "per_seed_records": records,
        "pass": pass_flag,
    }


def _make_tangent_figure(rep: dict) -> None:
    """Declared-cutoff singular-value figure for one representative seed."""
    xA = np.arange(1, len(rep["sv_A_r"]) + 1)
    xB = np.arange(1, len(rep["sv_B_r"]) + 1)
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    ax.plot(
        xA, rep["sv_A_r"], marker="o", ms=4.5, ls="", color="#1f4e79",
        label=r"$\sigma_i(A_r)$ (realified map tangent)",
    )
    ax.plot(
        xB, rep["sv_B_r"], marker="s", ms=4.5, ls="", color="#b03a2e",
        label=r"$\sigma_i(B_r)$ (realified pose tangent)",
    )
    for cut, name in ((rep["C1_width"], "C1=4"), (rep["C2_width"], "C2=8")):
        ax.axvline(cut + 0.5, color="#6da8d4", ls="--", lw=1.2,
                   label=f"declared cutoff {name}")
    ax.set_yscale("log")
    ax.set_xlabel("singular-value index")
    ax.set_ylabel("absolute singular value")
    ax.set_title(
        "Physical tangents (seed %d): A_r/B_r singular values\n"
        "vertical lines = DECLARED free-current cutoffs, NOT oracle ranks"
        % rep["seed"]
    )
    ax.legend(fontsize=8)
    ax.grid(alpha=0.25, which="both")
    fig.tight_layout()
    fig.savefig("figures/physical_tangents.png", dpi=160)
    plt.close(fig)


def _failed_phase(phase: str, exc: Exception) -> dict:
    return {
        "status": "failed",
        "label": pr.declared_note(),
        "error": {"type": type(exc).__name__, "message": str(exc)},
        "error_traceback": traceback.format_exc(),
        "pass": False,
    }


def main() -> int:
    os.makedirs("results", exist_ok=True)
    os.makedirs("figures", exist_ok=True)

    print("Phase 0: N=8 full-aperture sanity forward")
    try:
        cfg = physics.Config(N=8, aperture="full")
        model = physics.Model(cfg)
        phase0 = _run_phase0(model)
    except Exception as exc:  # noqa: BLE001 - record instead of crashing
        phase0 = _failed_phase("phase0", exc)
        model = None
    print("  status:", phase0.get("status"), "shapes:",
          phase0.get("forward", {}).get("A_realified_shape"),
          phase0.get("forward", {}).get("B_realified_shape"))

    print("Phase 1: E2/E3 physical tangents, seeds 201-212")
    if model is None:
        phase1 = _failed_phase("phase1", RuntimeError("model unavailable"))
    else:
        try:
            phase1 = _run_phase1(model)
        except Exception as exc:  # noqa: BLE001
            phase1 = _failed_phase("phase1", exc)
    agg1 = phase1.get("aggregate", {})
    print("  status:", phase1["status"], "| max T5a scaled:",
          _fmt(agg1.get("max_t5a_backward_scaled_residual", None), 3),
          "| T5b matches:",
          f"{agg1.get('t5b_matches_count')}/{agg1.get('n_checks_t5')}",
          "| max Thm3 resid:",
          _fmt(agg1.get("max_theorem3_residual_max_abs", None), 3))

    print("Phase 2: E5 physical shared-map frames, seeds 401-412")
    if model is None:
        phase2 = _failed_phase("phase2", RuntimeError("model unavailable"))
    else:
        try:
            phase2 = _run_phase2(model)
        except Exception as exc:  # noqa: BLE001
            phase2 = _failed_phase("phase2", exc)
    agg2 = phase2.get("aggregate", {})
    print("  status:", phase2["status"], "| max innovation rel:",
          _fmt(agg2.get("max_innovation_rel_residual", None), 3),
          "| max budget rel:",
          _fmt(agg2.get("max_budget_identity_rel_residual", None), 3),
          "| stack>=sum:",
          f"{agg2.get('stack_ge_sum_count')}/{agg2.get('n_seeds')}")

    if phase1.get("representative"):
        _make_tangent_figure(phase1["representative"])
        print("Saved: figures/physical_tangents.png")

    settings = {
        "label": pr.declared_note(),
        "date": "2026-09-05",
        "interpreter": (
            "/Volumes/migodam's-external-brain/Research/Inv_SLAM/experiments/"
            "idea_loops/loop_2026-09-04_02-58-31/"
            "experiment_geometry_lifted_trispace_som/.venv/bin/python"
        ),
        "physics": {
            "path": PHYSICS_DIR,
            "config": "Config(N=8, aperture='full')",
            "shapes": "A (288,9) complex, B (288,3) complex; "
                      "realified A_r (576,9), B_r (576,3)",
            "realification": (
                "sqrt(2)*[Re(J); Im(J)] via physics.realify_jacobian "
                "(unit-noise whitened); row order frequency->pose->"
                "illumination->receiver"
            ),
            "read_only_boundary": (
                "nothing under research/delegated/a2_physics/ is modified"
            ),
        },
        "phase1_seeds": list(range(201, 213)),
        "phase2_seeds": list(range(401, 413)),
        "declared_C": {
            "C0_dim": 0,
            "C1_dim": 4,
            "C2_dim": 8,
            "source": (
                "leading left singular vectors of realified full data tangent; "
                "declared candidate current directions, no oracle rank"
            ),
        },
        "phase2_frames": (
            "4 acquisition frames = the 4 frequencies of the full-wave model "
            "(shared map/pose parameters, per-frequency current)"
        ),
        "tolerances": {
            "rank": "absolute 10*eps*max(1,||B||_F)*n",
            "t5a_scaled_threshold": 100.0,
            "theorem3_max_abs": 1e-6,
            "prior_loewner_relative": 1e-8,
            "innovation_rel": 1e-6,
            "budget_rel": 1e-6,
            "psd": "max(1e-10, 20*eps*max(1,||B||_F^2)*n) per check",
            "G_pd_min_eig": 1e-10,
        },
        "theorem7_metric": (
            "declared pose error metric M_x = I_3 (unit metric on the three "
            "pose coordinates); bias directions from residual left singular "
            "vector or random unit fallback"
        ),
        "claims_boundary": (
            "No physical-rank, envelope-nonidentifiability, nonlinear-basin, "
            "or method-superiority claim is made by this smoke runner."
        ),
    }

    verdict = {
        "phase0": bool(
            phase0["status"] == "ok"
            and phase0.get("row_structure", {}).get(
                "blocks_checked", {}
            ).get("phase0_shapes_ok")
        ),
        "phase1": bool(phase1["pass"]),
        "phase2": bool(phase2["pass"]),
        "statuses": {
            "phase0": phase0["status"],
            "phase1": phase1["status"],
            "phase2": phase2["status"],
        },
    }

    doc = {
        "title": (
            "OPTIONAL physical second-stage tangent runner "
            "(E2 Theorem 5, E3 Theorems 3/7, E5 Theorem 9/Corollary 10)"
        ),
        "label": pr.declared_note(),
        "package": "a2val",
        "settings": settings,
        "results": {
            "phase0": phase0,
            "phase1": phase1,
            "phase2": phase2,
        },
        "verdict": verdict,
    }
    with open("results/physical_optional_results.json", "w") as fh:
        json.dump(_jsonable(doc), fh, indent=2)
        fh.write("\n")

    rows = []

    def _add(phase, status, key_values, passed, note):
        rows.append({
            "phase": phase,
            "status": status,
            "key_values": key_values,
            "pass": passed,
            "note": note,
        })

    f0 = phase0.get("forward", {})
    _add(
        "phase0",
        phase0["status"],
        (
            f"A_r={f0.get('A_realified_shape')}, "
            f"B_r={f0.get('B_realified_shape')}, "
            f"wall={_fmt(f0.get('wall_seconds'), 3)}s"
        ),
        bool(verdict["phase0"]),
        "import physics; N=8 full-aperture forward; realified shape sanity",
    )
    _add(
        "phase1",
        phase1["status"],
        (
            f"maxT5aScaled={_fmt(agg1.get('max_t5a_backward_scaled_residual'), 3)}, "
            f"T5b={agg1.get('t5b_matches_count')}/{agg1.get('n_checks_t5')}, "
            f"maxThm3={_fmt(agg1.get('max_theorem3_residual_max_abs'), 3)}, "
            f"priorMinEig=[{_fmt(agg1.get('min_prior_min_eig_K_eL_minus_K_e'), 3)}, "
            f"{_fmt(agg1.get('min_prior_min_eig_K0_minus_K_eL'), 3)}]"
        ),
        bool(phase1["pass"]),
        "E2 T5a/T5b + E3 Theorem3/prior/Theorem7 on physical tangents; "
        "declared cutoffs, no oracle rank",
    )
    _add(
        "phase2",
        phase2["status"],
        (
            f"maxInnovRel={_fmt(agg2.get('max_innovation_rel_residual'), 3)}, "
            f"maxBudgetRel={_fmt(agg2.get('max_budget_identity_rel_residual'), 3)}, "
            f"stack>=sum={agg2.get('stack_ge_sum_count')}/{agg2.get('n_seeds')}"
        ),
        bool(phase2["pass"]),
        "E5 shared-map innovation + Corollary-10 rank budget on 4-frequency "
        "physical frames; declared C4/C8 enlargement",
    )

    with open("results/physical_optional_summary.csv", "w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=["phase", "status", "key_values", "pass", "note"]
        )
        writer.writeheader()
        for r in rows:
            writer.writerow(r)

    commands = [
        "cd " + ROOT,
        (
            "/Volumes/migodam's-external-brain/Research/Inv_SLAM/experiments/"
            "idea_loops/loop_2026-09-04_02-58-31/"
            "experiment_geometry_lifted_trispace_som/.venv/bin/python "
            "tests/test_physical_runner.py"
        ),
        (
            "/Volumes/migodam's-external-brain/Research/Inv_SLAM/experiments/"
            "idea_loops/loop_2026-09-04_02-58-31/"
            "experiment_geometry_lifted_trispace_som/.venv/bin/python "
            "scripts/run_physical_optional.py"
        ),
    ]
    with open("results/physical_optional_commands.txt", "w") as fh:
        fh.write("Exact commands run (physical optional second stage)\n")
        fh.write("\n".join(commands))
        fh.write("\n")

    print("Saved: results/physical_optional_results.json, "
          "results/physical_optional_summary.csv, "
          "results/physical_optional_commands.txt")
    print("Verdict:", verdict)
    return 0 if all(
        (verdict["phase0"], verdict["phase1"], verdict["phase2"])
    ) else 1


if __name__ == "__main__":
    raise SystemExit(main())
