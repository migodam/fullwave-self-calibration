#!/usr/bin/env python
"""E5 physical acquisition extension: 8 receiver-configuration candidates.

Bounded tangent-level supplement on the N=8 scalar full-wave core.  The
decision criterion is the exact shared-map pose-information stack
``J(A,B) = B^T (I - P_Ran(A)) B`` over fixed per-frequency tangent blocks
(frequency rows = frames, shared map/pose parameters).  All candidate
acquisition decisions use ONLY the fixed declared reference state
``seed 501``.  Truth states (seeds 1001 and 1002) are used only to evaluate
the predeclared policies; no truth tangent, oracle pose score or nonlinear
reconstruction is used in the decision.

Eight candidate actions are 12-channel receiver configurations supported
directly by ``physics.Config(receiver_angles_full=...)`` with N=8 and
aperture='full'.  The old model is the full-ring configuration at the two
lowest frequencies; each candidate action adds the same new frequency with a
different receiver configuration.

Usage (from the replication root):
    python scripts/run_e5_physical_acquisition.py

Outputs:
    results/e5_physical_acquisition_results.json
    results/e5_physical_acquisition_summary.csv
    figures/e5_physical_acquisition_2x2.png
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

from a2val import e5, physical_runner as pr


DECISION_SEED = 501
TRUTH_SEEDS = (1001, 1002)
OLD_FREQS = (0, 1)
NEW_FREQ = 2
K = 3
LAMBDA = 0.1
RANDOM_SEEDS = list(range(401, 413))


def _ang(*deg):
    return tuple(np.radians(float(d)) for d in deg)


def _lin_deg(a, b, n=12):
    return tuple(np.linspace(np.radians(float(a)), np.radians(float(b)), n))


def _concat(*groups):
    return tuple(np.concatenate([np.asarray(g, dtype=float) for g in groups]))


RX_CONFIGS = {
    "ring12_full": _ang(*[30.0 * m for m in range(12)]),
    "ring12_rot15": _ang(*[30.0 * m + 15.0 for m in range(12)]),
    "arc12_front": _lin_deg(-45.0, 45.0),
    "arc12_right": _lin_deg(45.0, 135.0),
    "arc12_back": _lin_deg(135.0, 225.0),
    "arc12_left": _lin_deg(225.0, 315.0),
    "arc12_top": _lin_deg(-90.0, 90.0),
    "twin_arcs_fb": _concat(_lin_deg(-45.0, 45.0, 6),
                            _lin_deg(135.0, 225.0, 6)),
}


def _sym(X):
    X = np.asarray(X, dtype=float)
    return 0.5 * (X + X.T)


def _fro(X):
    return float(np.linalg.norm(np.asarray(X, dtype=float), "fro"))


def _min_eig(X):
    return float(np.linalg.eigvalsh(_sym(X))[0])


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


def _reference(seed: int):
    """Declared reference rule identical to the physical runners."""
    rng = np.random.default_rng(int(seed))
    alpha = np.clip(0.2 + 0.05 * rng.standard_normal(9), -0.5, 0.5)
    x = np.clip(0.02 * rng.standard_normal(3), -0.05, 0.05)
    return alpha, x


def _make_models():
    models = {}
    for name, angles in RX_CONFIGS.items():
        cfg = physics.Config(
            N=8, aperture="full", receiver_angles_full=np.asarray(angles)
        )
        models[name] = physics.Model(cfg)
    return models


def _tangent_block(model, alpha, x, freq_ids):
    """Realified (A,B) tangent and per-block metadata for selected freqs."""
    out = model.forward(alpha, x, freq_ids=list(freq_ids), jacobian=True)
    A_r = physics.realify_jacobian(out["A"])
    B_r = physics.realify_jacobian(out["B"])
    frames = pr.split_frequency_frames(A_r, B_r, out["blocks"])
    return {
        "out": out,
        "A_r": A_r,
        "B_r": B_r,
        "frames": frames,
        "state_residual_rel": [st["state_residual_rel"]
                               for st in out["states"]],
        "work": out["work"],
    }


def _block_state(models, alpha, x):
    """Old (full ring, freqs 0-1) and eight candidate (new freq) blocks."""
    old = _tangent_block(models["ring12_full"], alpha, x, OLD_FREQS)
    cand = {}
    for name, model in models.items():
        if name == "ring12_full":
            # Candidate action with the default ring is the same receiver
            # geometry at a NEW frequency (separate forward, freq 2 only).
            pass
        blk = _tangent_block(model, alpha, x, [NEW_FREQ])
        cand[name] = blk
    return old, cand


def _old_frames(blk):
    """Two baseline frequency frames from the full-ring block."""
    if isinstance(blk, (list, tuple)) and blk and isinstance(
        blk[0], tuple
    ):
        return list(blk)
    if len(blk["frames"]) != len(OLD_FREQS):
        raise RuntimeError("expected one frame per old frequency")
    return list(blk["frames"])


def _candidate_frame(blk):
    """Single candidate frequency frame (realified rows)."""
    if isinstance(blk, tuple) and len(blk) == 2 and isinstance(
        blk[0], np.ndarray
    ):
        return blk
    if len(blk["frames"]) != 1:
        raise RuntimeError("expected a single candidate frequency frame")
    return blk["frames"][0]


def _exact_logdet(old, cand, cand_subset, Lam):
    """Shared-map exact criterion over the old stack + selected candidates."""
    frames = _old_frames(old)
    for name in cand_subset:
        frames.append(_candidate_frame(cand[name]))
    A, B = e5._stack_AB(frames)
    J = e5._pose_info(A, B)
    return float(e5._logdet(_sym(J + Lam))), J


def _independent_logdet(old, cand, cand_subset, Lam):
    """Incorrectly independent per-frame map criterion (sum of J_l)."""
    frames_old = _old_frames(old)
    A_old, B_old = e5._stack_AB(frames_old)
    J_base = e5._pose_info(A_old, B_old)
    total = J_base.copy()
    for name in cand_subset:
        a, b = _candidate_frame(cand[name])
        total = total + e5._pose_info(a, b)
    return float(e5._logdet(_sym(total + Lam))), total


def _all_logdets(old, cand, Lam):
    rows = []
    for comb in itertools.combinations(range(8), K):
        subset = [list(RX_CONFIGS.keys())[i] for i in comb]
        ex, Jex = _exact_logdet(old, cand, subset, Lam)
        ind, Jind = _independent_logdet(old, cand, subset, Lam)
        rows.append({
            "subset": subset,
            "indices": list(comb),
            "exact_logdet": ex,
            "independent_logdet": ind,
            "J_exact_eig": np.linalg.eigvalsh(_sym(Jex)).tolist(),
            "J_ind_eig": np.linalg.eigvalsh(_sym(Jind)).tolist(),
        })
    return rows


def _policy_rows(all_rows, old, cand, Lam, n_candidates=8):
    """Rank-1/stack policy records; scoring uses precomputed tangent blocks."""
    exact_by_set = {frozenset(r["subset"]): r["exact_logdet"]
                    for r in all_rows}
    ind_by_set = {frozenset(r["subset"]): r["independent_logdet"]
                  for r in all_rows}
    names = list(RX_CONFIGS.keys())

    def score_exact(subset):
        key = frozenset(subset)
        if key in exact_by_set:
            return exact_by_set[key]
        return float(_exact_logdet(old, cand, list(subset), Lam)[0])

    def score_ind(subset):
        key = frozenset(subset)
        if key in ind_by_set:
            return ind_by_set[key]
        return float(_independent_logdet(old, cand, list(subset), Lam)[0])

    def greedy(score):
        chosen = []
        remaining = set(range(n_candidates))
        while len(chosen) < K:
            best = None
            bestv = None
            for i in remaining:
                val = score([names[j] for j in chosen + [i]])
                if best is None or val > bestv + 1e-15:
                    best = i
                    bestv = val
            chosen.append(best)
            remaining.remove(best)
        return [names[j] for j in chosen]

    def pair_lookahead(score):
        best_pair = None
        bestv = None
        for pair in itertools.combinations(range(n_candidates), 2):
            val = score([names[j] for j in pair])
            if best_pair is None or val > bestv + 1e-15:
                best_pair = pair
                bestv = val
        chosen = list(best_pair)
        finish = None
        finishv = None
        for i in range(n_candidates):
            if i in chosen:
                continue
            val = score([names[j] for j in chosen + [i]])
            if finish is None or val > finishv + 1e-15:
                finish = i
                finishv = val
        return [names[j] for j in chosen + [finish]]

    def exhaustive(score):
        best = None
        bestv = None
        for r in all_rows:
            val = score(r["subset"])
            if best is None or val > bestv + 1e-15:
                best = r["subset"]
                bestv = val
        return best

    greedy_ex = greedy(score_exact)
    pair_ex = pair_lookahead(score_exact)
    ex_ex = exhaustive(score_exact)
    greedy_ind = greedy(score_ind)
    ex_ind = exhaustive(score_ind)
    random_rows = []
    for seed_r in RANDOM_SEEDS:
        rnd = _random.Random(seed_r)
        comb = tuple(sorted(rnd.sample(range(n_candidates), K)))
        subset = [names[j] for j in comb]
        random_rows.append({
            "seed": seed_r,
            "subset": subset,
            "exact_logdet": score_exact(subset),
            "independent_logdet": score_ind(subset),
        })
    random_best = max(random_rows, key=lambda t: t["exact_logdet"])
    return {
        "greedy_exact": {"subset": greedy_ex,
                         "logdet": score_exact(greedy_ex)},
        "pair_exact": {"subset": pair_ex, "logdet": score_exact(pair_ex)},
        "exhaustive_exact": {"subset": ex_ex,
                             "logdet": score_exact(ex_ex)},
        "greedy_independent": {"subset": greedy_ind,
                               "logdet": score_exact(greedy_ind),
                               "independent_logdet": score_ind(greedy_ind)},
        "exhaustive_independent": {"subset": ex_ind,
                                   "logdet": score_exact(ex_ind),
                                   "independent_logdet": score_ind(ex_ind)},
        "random_best_of_12_exact": random_best,
        "random_rows": random_rows,
        "greedy_gap_exact": float(score_exact(ex_ex) - score_exact(greedy_ex)),
        "pair_gap_exact": float(score_exact(ex_ex) - score_exact(pair_ex)),
        "independent_gap_exact": float(
            score_exact(ex_ex) - score_exact(ex_ind)
        ),
        "all_rows": all_rows,
    }


def _truth_evaluation(old_t, cand_t, policies_ref, Lam):
    """Truth-only evaluation of the reference policy subsets."""
    names = list(RX_CONFIGS.keys())
    all_rows = _all_logdets(old_t, cand_t, Lam)
    by_set = {frozenset(r["subset"]): r for r in all_rows}
    ref_ex = by_set[frozenset(policies_ref["exhaustive_exact"]["subset"])]
    keys = ("greedy_exact", "pair_exact", "random_best_of_12_exact",
            "exhaustive_independent")
    out = {
        "exhaustive_reference_subset_truth": ref_ex,
        "policy_truth_rows": {},
    }
    for key in keys:
        subset = policies_ref[key]["subset"]
        out["policy_truth_rows"][key] = by_set[frozenset(subset)]
    # Reference-truth rank agreement over all 56 triples.
    order_ref = [r["exact_logdet"] for r in policies_ref["all_rows"]]
    order_tr = [r["exact_logdet"] for r in all_rows]
    # Pairwise index orderings are identical because all_rows iterates the
    # same itertools.combinations order for a fixed n_candidates=8.
    rk_ref = np.argsort(np.argsort(order_ref))
    rk_tr = np.argsort(np.argsort(order_tr))
    out["spearman_exact_truth"] = float(np.corrcoef(rk_ref, rk_tr)[0, 1])
    out["max_abs_logdet_drift"] = float(np.max(np.abs(
        np.asarray(order_tr) - np.asarray(order_ref)
    )))
    out["n_triples"] = len(all_rows)
    out["all_rows"] = all_rows
    out["names"] = names
    return out


def _fig(ref_pol, evals):
    names = list(RX_CONFIGS.keys())
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.2))
    ax = axes[0]
    labels = ["greedy\nexact", "pair\nexact", "exhaustive\nexact",
              "greedy\nindependent", "random\nbest-of-12"]
    subsets = [
        ref_pol["greedy_exact"]["subset"],
        ref_pol["pair_exact"]["subset"],
        ref_pol["exhaustive_exact"]["subset"],
        ref_pol["greedy_independent"]["subset"],
        ref_pol["random_best_of_12_exact"]["subset"],
    ]
    x = np.arange(len(labels))
    vals_ref = [ref_pol["greedy_exact"]["logdet"],
                ref_pol["pair_exact"]["logdet"],
                ref_pol["exhaustive_exact"]["logdet"],
                ref_pol["greedy_independent"]["logdet"],
                ref_pol["random_best_of_12_exact"]["exact_logdet"]]
    ax.bar(x, vals_ref, width=0.58, color="#1b6ca8", alpha=0.85)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=7)
    ax.set_ylabel("exact shared-map logdet (reference seed 501)")
    ax.set_title("Reference policy objective (8 receiver configs, K=3)")
    ax.grid(alpha=0.25, axis="y")
    for xi, ss in zip(x, subsets):
        ax.annotate(",".join(str(names.index(s)) for s in ss),
                    (xi, vals_ref[xi]), textcoords="offset points",
                    xytext=(0, 4), ha="center", fontsize=6)
    ax = axes[1]
    # truth-policy panel (reference subsets at truth seeds)
    x2 = np.arange(len(labels) - 1)  # omit duplicate exhaustive-independent
    width = 0.36
    keys = ("greedy_exact", "pair_exact", "random_best_of_12_exact",
            "exhaustive_independent")
    for j, seed in enumerate(evals):
        ev = evals[seed]
        truth_vals = [ev["policy_truth_rows"][k]["exact_logdet"]
                      for k in keys]
        ax.bar(x2 + (j - 0.5) * width, truth_vals, width=width,
               label="truth seed %d" % seed)
    ax.set_xticks(x2)
    ax.set_xticklabels(keys, fontsize=6.5)
    ax.set_ylabel("exact shared-map logdet at truth")
    ax.set_title("Truth evaluation of reference subsets\n"
                 "(truth never used in selection)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.25, axis="y")
    fig.tight_layout()
    os.makedirs("figures", exist_ok=True)
    fig.savefig("figures/e5_physical_acquisition_2x2.png", dpi=150)
    plt.close(fig)


def main():
    t_start = time.perf_counter()
    os.makedirs("results", exist_ok=True)
    os.makedirs("figures", exist_ok=True)
    models = _make_models()
    Lam = LAMBDA * np.eye(3)

    # ---- Decision stage (reference seed only) -----------------------------
    alpha_r, x_r = _reference(DECISION_SEED)
    old_r, cand_r = _block_state(models, alpha_r, x_r)
    old_frames_r = _old_frames(old_r)
    A_base, B_base = e5._stack_AB(old_frames_r)
    G_base = A_base.T @ A_base
    policies = _policy_rows(
        _all_logdets(old_r, cand_r, Lam), old_r, cand_r, Lam
    )
    ref_forward_work = {
        "baseline_forward": dict(old_r["work"]),
        "candidate_forwards": {name: dict(blk["work"])
                               for name, blk in cand_r.items()},
    }
    ref_state_res = {
        "baseline_max": float(np.max(old_r["state_residual_rel"])),
        "candidate_max": {
            name: float(np.max(blk["state_residual_rel"]))
            for name, blk in cand_r.items()
        },
    }
    candidate_metrics = {}
    for name, blk in cand_r.items():
        a, b = _candidate_frame(blk)
        J = e5._pose_info(a, b)
        candidate_metrics[name] = {
            "n_receivers": int(a.shape[0] // 12),  # real rows = 12*R
            "rows_real": int(a.shape[0]),
            "J_eig": np.linalg.eigvalsh(_sym(J)).tolist(),
            "min_eig_J": _min_eig(J),
        }

    # ---- Truth evaluation (seeds 1001/1002 only) --------------------------
    evals = {}
    total_work = {
        "rhs_solves_total": ref_forward_work["baseline_forward"][
            "rhs_solves_total"],
        "operator_products": ref_forward_work["baseline_forward"][
            "operator_products"],
    }
    for name, w in ref_forward_work["candidate_forwards"].items():
        total_work["rhs_solves_total"] += w["rhs_solves_total"]
        total_work["operator_products"] += w["operator_products"]
    for seed in TRUTH_SEEDS:
        alpha_t, x_t = _reference(seed)
        old_t, cand_t = _block_state(models, alpha_t, x_t)
        ev = _truth_evaluation(old_t, cand_t, policies, Lam)
        ev["forward_work"] = {
            "baseline": dict(old_t["work"]),
            "candidates": {name: dict(blk["work"])
                           for name, blk in cand_t.items()},
        }
        ev["state_residual_rel_max"] = float(np.max(
            old_t["state_residual_rel"]
        ))
        evals[seed] = ev
        for name, w in ev["forward_work"]["candidates"].items():
            total_work["rhs_solves_total"] += w["rhs_solves_total"]
            total_work["operator_products"] += w["operator_products"]
        total_work["rhs_solves_total"] += (
            ev["forward_work"]["baseline"]["rhs_solves_total"]
        )
        total_work["operator_products"] += (
            ev["forward_work"]["baseline"]["operator_products"]
        )

    # Evaluation-only logdet counts (no extra forwards).
    total_work["offline_exact_logdet_evaluations"] = (
        len(policies["all_rows"]) * (1 + len(TRUTH_SEEDS))
    )

    _fig(policies, evals)

    doc = {
        "title": (
            "E5 physical acquisition extension: 8 receiver configurations, "
            "reference-decision and separate truth evaluation"
        ),
        "date": "2026-09-06",
        "package": "a2val (replication copy) + shared physics core",
        "settings": {
            "physics": {"module": physics.__name__, "file": physics.__file__},
            "config": "physics.Config(N=8, aperture='full', "
                      "receiver_angles_full=action-specific 12 angles)",
            "decision_seed": DECISION_SEED,
            "truth_seeds": list(TRUTH_SEEDS),
            "old_model": {
                "receiver_config": "ring12_full",
                "frequencies": list(OLD_FREQS),
                "frames_per_frequency": 2,
                "shared_map_pose": True,
            },
            "candidate_action": {
                "receiver_config": "one of the 8 listed configs",
                "frequency": NEW_FREQ,
                "candidates": list(RX_CONFIGS.keys()),
            },
            "budget_K": K,
            "lambda_x": LAMBDA,
            "objective_note": (
                "g(S)=logdet(J_stack(old + candidates(S)) + 0.1*I3), "
                "exact shared-map criterion; no forward/adjoint is required "
                "for triple scoring after each candidate tangent block is "
                "precomputed"
            ),
            "reference_state_note": (
                "alpha=clip(0.2+0.05*N(0,1),[-0.5,0.5]), "
                "x=clip(0.02*N(0,1),[-0.05,0.05]) from default_rng(seed); "
                "no nonlinear inversion or oracle pose score"
            ),
            "truth_boundary": (
                "truth seeds used only for evaluating predeclared subsets; "
                "no truth tangent enters any decision"
            ),
            "claim_boundary": (
                "tangent-level physical acquisition diagnostic on one "
                "reference and two truth states; no nonlinear reconstruction, "
                "no method-superiority or universal frame-count claim"
            ),
        },
        "candidate_receiver_configs": {
            name: {"angles_deg": [float(np.degrees(a)) for a in angles],
                   "n_receivers": len(angles)}
            for name, angles in RX_CONFIGS.items()
        },
        "reference_state": {
            "seed": DECISION_SEED,
            "alpha": alpha_r.tolist(),
            "x": x_r.tolist(),
        },
        "reference": {
            "candidate_metrics": candidate_metrics,
            "old_base_J_eig": np.linalg.eigvalsh(_sym(
                e5._pose_info(A_base, B_base)
            )).tolist(),
            "old_base_G_min_eig": _min_eig(G_base),
            "policies": policies,
            "forward_work": ref_forward_work,
            "state_residual_rel_max": ref_state_res,
        },
        "truth_evaluations": evals,
        "charged_work": total_work,
        "runtime_wall_seconds": time.perf_counter() - t_start,
    }
    with open("results/e5_physical_acquisition_results.json", "w") as fh:
        json.dump(_jsonable(doc), fh, indent=2)
        fh.write("\n")

    csv_rows = []
    for r in policies["all_rows"]:
        csv_rows.append({
            "stage": "reference",
            "seed": DECISION_SEED,
            "subset": "+".join(r["subset"]),
            "exact_logdet": r["exact_logdet"],
            "independent_logdet": r["independent_logdet"],
        })
    for seed, ev in evals.items():
        for r in ev["all_rows"]:
            csv_rows.append({
                "stage": "truth",
                "seed": seed,
                "subset": "+".join(r["subset"]),
                "exact_logdet": r["exact_logdet"],
                "independent_logdet": r["independent_logdet"],
            })
    keys = ["stage", "seed", "subset", "exact_logdet", "independent_logdet"]
    with open("results/e5_physical_acquisition_summary.csv", "w",
              newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=keys)
        writer.writeheader()
        for r in csv_rows:
            writer.writerow({k: _fmtv(r[k]) for k in keys})

    print("E5 physical acquisition extension")
    print("  greedy exact:", policies["greedy_exact"])
    print("  pair exact:", policies["pair_exact"])
    print("  exhaustive exact:", policies["exhaustive_exact"])
    print("  independent exhaustive:", policies["exhaustive_independent"])
    print("  greedy gap:", policies["greedy_gap_exact"])
    print("  independent-map gap:", policies["independent_gap_exact"])
    for seed, ev in evals.items():
        print("  truth %d: spearman exact=%.4f drift=%.4g" % (
            seed, ev["spearman_exact_truth"], ev["max_abs_logdet_drift"]
        ))
    print("  charged RHS solves:", total_work["rhs_solves_total"],
          "| operator products:", total_work["operator_products"])
    print("Saved: results/e5_physical_acquisition_results.json, "
          "results/e5_physical_acquisition_summary.csv, "
          "figures/e5_physical_acquisition_2x2.png")
    return 0


def _fmtv(v):
    if isinstance(v, float):
        return f"{v:.10g}"
    if isinstance(v, (list, tuple, dict)):
        return json.dumps(_jsonable(v))
    return v


if __name__ == "__main__":
    raise SystemExit(main())
