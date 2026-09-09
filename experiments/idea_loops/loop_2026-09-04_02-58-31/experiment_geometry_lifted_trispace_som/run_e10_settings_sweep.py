"""E10 -- settings sweep of the retained-rank SOM lift and reduced self-calibration.

E10 broadens the E5/E9 core claims across the receiver count M, the
wavenumber k, and limited-aperture receiver arcs, using the canonical E5/E9
base scene (N=16, K=3 Gaussian basis, alpha_true=[1.5,2.0,0.0],
alpha_init=[1.0,1.0,0.0], p_true=[0.08,-0.06,0.05], p_init=[0,0,0], 30 dB).

Settings:
    M8_k12_full            M=8  k=12 full circle
    M12_k12_full           M=12 k=12 full circle
    M16_k12_full           M=16 k=12 full circle
    M12_k8_full            M=12 k=8  full circle
    M12_k16_full           M=12 k=16 full circle
    M12_k12_half           M=12 k=12 half-circle arc (pi aperture at 0)
    M16_k12_threequarter   M=16 k=12 3/4-circle arc (3*pi/2 aperture at 0)

Per setting (deterministic, at (alpha_init, p_init)):
  * rank checks of G_s (singular values, gaps, numerical rank),
  * unrestricted (full-current) vs retained-rank (r=4,6) data lifts of the
    canonical receiver-side pose perturbation b_pert = B_r[:,0],
  * visible-pose subspace sweep over retained rank r in {2,...,8} (r <= M).

Nonlinear runs (seeds 0,1,2; E8/E9 noise: ||noise||/||d_true|| = 10^-1.5):
  wrongpose  (alpha only, p fixed at p_init=0)
  known_pose (alpha only, p fixed at p_true -- oracle contrast baseline)
  known_alpha(p only, alpha fixed at alpha_true -- pose-only oracle)
  direct     ([alpha; p])
  reduced_r4 ([alpha; q], p = p_init + V_vis(r4) @ q)
  reduced_r6 ([alpha; q], p = p_init + V_vis(r6) @ q)
plus, for M12_k12_full only, reduced_r{r} for every r in {2,...,8}.

Pose errors are always decomposed on the fixed r=6 V_vis at (alpha_init,
p_init) so hidden/visible components are comparable across methods; for
reduced_r6 that basis is also the solver's own basis.

CPU only, deterministic, reuses the E5/E9 forward/jacobian/subspace code.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import least_squares

from geom_som_core import Params, pixel_grid, green_domain_matrix
from run_e5 import (
    basis_matrix,
    forward,
    jacobians,
    state_witness,
    _realified_residual,
)
from run_e5_final import visible_pose_subspace, LS_KWARGS


HERE = Path(__file__).resolve().parent

SNR_DB = 30.0
NOISE_RATIO = float(10.0 ** (-SNR_DB / 20.0))  # 10^-1.5
NOISE_SEEDS = (0, 1, 2)
SUCCESS_OPTIMALITY = 1e-7
SV_RANK_REL = 1e-10  # numerical-rank threshold for G_s singular values

BASE_N = 16
ALPHA_TRUE = np.array([1.5, 2.0, 0.0], dtype=float)
ALPHA_INIT = np.array([1.0, 1.0, 0.0], dtype=float)
P_TRUE = np.array([0.08, -0.06, 0.05], dtype=float)
P_INIT = np.zeros(3)

BASE_P = dict(
    N=BASE_N,
    R_r=1.6,
    R_t=2.0,
    phi0=0.7,
    s=0.12,
)


def _angles(spec: str, M: int) -> tuple | None:
    """Receiver-angle spec -> tuple of body-frame angles (or None = full)."""
    if spec == "full":
        return None
    if spec == "half":
        return tuple(float(a) for a in np.linspace(-np.pi / 2.0, np.pi / 2.0, M))
    if spec == "threequarter":
        return tuple(
            float(a) for a in np.linspace(-3.0 * np.pi / 4.0, 3.0 * np.pi / 4.0, M)
        )
    raise ValueError(f"unknown aperture spec {spec!r}")


SETTINGS = (
    dict(label="M8_k12_full", M=8, k=12.0, aperture="full"),
    dict(label="M12_k12_full", M=12, k=12.0, aperture="full"),
    dict(label="M16_k12_full", M=16, k=12.0, aperture="full"),
    dict(label="M12_k8_full", M=12, k=8.0, aperture="full"),
    dict(label="M12_k16_full", M=12, k=16.0, aperture="full"),
    dict(label="M12_k12_half", M=12, k=12.0, aperture="half"),
    dict(label="M16_k12_threequarter", M=16, k=12.0, aperture="threequarter"),
)

SETTING_LABELS = [s["label"] for s in SETTINGS]
FULL_CIRCLE_M12 = "M12_k12_full"
RANK_SWEEP_RS = (2, 3, 4, 5, 6, 7, 8)


def make_params(setting: dict) -> Params:
    angles = _angles(setting["aperture"], int(setting["M"]))
    if angles is not None:
        assert len(angles) == int(setting["M"])
    return Params(
        N=BASE_N,
        k=float(setting["k"]),
        M=int(setting["M"]),
        R_r=BASE_P["R_r"],
        R_t=BASE_P["R_t"],
        phi0=BASE_P["phi0"],
        s=BASE_P["s"],
        receiver_angles=angles,
    )


def aperture_span(setting: dict) -> float:
    """Aperture arc length in radians (2*pi for full circle)."""
    if setting["aperture"] == "full":
        return 2.0 * np.pi
    if setting["aperture"] == "half":
        return np.pi
    if setting["aperture"] == "threequarter":
        return 1.5 * np.pi
    raise ValueError(setting["aperture"])


# ---------------------------------------------------------------------------
# Noise (exactly the E8/E9 convention: ||noise||/||d_true|| = 10^-1.5)
# ---------------------------------------------------------------------------
def make_noise(seed: int, M: int) -> np.ndarray:
    rng = np.random.default_rng(int(seed))
    return (
        rng.standard_normal(M) + 1j * rng.standard_normal(M)
    ) / np.sqrt(2.0)


def scaled_observation(
    d_true: np.ndarray, raw: np.ndarray, snr_db: float = SNR_DB
) -> tuple[np.ndarray, float]:
    ratio = float(10.0 ** (-snr_db / 20.0))
    scale = ratio * float(np.linalg.norm(d_true)) / float(np.linalg.norm(raw))
    return d_true + scale * raw, ratio


# ---------------------------------------------------------------------------
# Deterministic linear / rank / lift diagnostics
# ---------------------------------------------------------------------------
def rank_check(jb: dict, M: int) -> dict:
    G_s = jb["G_s"]
    sv = np.asarray(np.linalg.svd(G_s, compute_uv=False), dtype=float)
    assert sv.shape[0] == M
    gaps = [float(sv[j] - sv[j + 1]) for j in range(len(sv) - 1)]
    num_rank = int(np.count_nonzero(sv > SV_RANK_REL * sv[0]))
    return {
        "singular_values": [float(v) for v in sv],
        "gaps_consecutive": gaps,
        "numerical_rank": num_rank,
        "rank_tol_relative": SV_RANK_REL,
        "note": "rank = count(sv > 1e-10 * sv[0]); gaps are sv[j]-sv[j+1]",
    }


def lift_check(jb: dict, r: int) -> dict:
    """Retained-rank lift of b_pert = B_r[:,0] through V_r = first r of G_s."""
    G_s = jb["G_s"]
    b = np.asarray(jb["B_r"][:, 0], dtype=complex).copy()
    norm_b = float(np.linalg.norm(b))
    _, _, Vh = np.linalg.svd(G_s, full_matrices=False)
    V_r = Vh.conj().T[:, :r]
    Q = G_s @ V_r
    c = np.linalg.pinv(Q) @ b
    J_ret = V_r @ c
    R = Q @ c - b
    QH_R = Q.conj().T @ R
    return {
        "retained_rank": int(r),
        "norm_b_pert": norm_b,
        "residual_ratio_ret": float(np.linalg.norm(R) / norm_b),
        "orthogonality_QH_R_max_abs": float(np.max(np.abs(QH_R))),
        "orthogonality_QH_R_norm_over_R_norm": float(
            np.linalg.norm(QH_R) / np.linalg.norm(R)
        ),
        "norm_J_ret": float(np.linalg.norm(J_ret)),
    }


def unrestricted_lift_check(jb: dict) -> dict:
    G_s = jb["G_s"]
    b = np.asarray(jb["B_r"][:, 0], dtype=complex).copy()
    norm_b = float(np.linalg.norm(b))
    J_full = np.linalg.pinv(G_s) @ b
    R = G_s @ J_full - b
    return {
        "retained_rank": None,
        "norm_b_pert": norm_b,
        "residual_ratio_full": float(np.linalg.norm(R) / norm_b),
        "norm_J_full": float(np.linalg.norm(J_full)),
    }


def subspace_sweep_json(
    sp: dict, label: str
) -> dict:
    return {
        "setting": label,
        "rank": int(sp["rank"]),
        "n_vis": int(sp["n_vis"]),
        "hidden_rank": int(sp["hidden_rank"]),
        "hid_sv": [float(v) for v in sp["hid_sv"]],
        "threshold": float(sp["threshold"]),
        "hidden_directions": sp["hidden_directions"],
        "V_vis_columns": np.asarray(sp["V_vis"]).tolist(),
    }


# ---------------------------------------------------------------------------
# One nonlinear scipy TRF run (E5-final fixed-V_vis reduced pattern)
# ---------------------------------------------------------------------------
def run_case(
    method: str,
    reduced_r: int | None,
    seed: int,
    setting: dict,
    P: Params,
    d_obs: np.ndarray,
    xs: np.ndarray,
    h: float,
    k: float,
    G_D: np.ndarray,
    Phi: np.ndarray,
    solver_subspace: dict | None,
    decomp_subspace: dict,
) -> dict:
    """Run one least-squares case; method in the E10 five-way method set."""
    K = Phi.shape[1]
    p_init = P_INIT
    V_dec = decomp_subspace["V_vis"]  # r=6 basis, used for decomposition
    V_sol = None if solver_subspace is None else solver_subspace["V_vis"]

    def xy(theta: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """(alpha, p) for the method's parameterization."""
        if method == "wrongpose":
            return np.asarray(theta, dtype=float), P_INIT
        if method == "known_pose":
            return np.asarray(theta, dtype=float), P_TRUE
        if method == "known_alpha":
            return ALPHA_TRUE, np.asarray(theta, dtype=float)
        if method == "direct":
            return np.asarray(theta[:K], float), np.asarray(theta[K:], float)
        # reduced
        alpha = np.asarray(theta[:K], dtype=float)
        q = np.asarray(theta[K:], dtype=float)
        return alpha, p_init + V_sol @ q

    def residual(theta: np.ndarray) -> np.ndarray:
        alpha, p = xy(theta)
        return _realified_residual(
            forward(alpha, p, xs, h, k, G_D, Phi, P), d_obs
        )

    def jacobian(theta: np.ndarray) -> np.ndarray:
        alpha, p = xy(theta)
        jb = jacobians(alpha, p, xs, h, k, G_D, Phi, P)
        if method in ("wrongpose", "known_pose"):
            return jb["J_alpha_real"]
        if method == "known_alpha":
            return jb["B_real"]
        if method == "direct":
            return np.hstack([jb["J_alpha_real"], jb["B_real"]])
        return np.hstack(
            [jb["J_alpha_real"], jb["B_real"] @ V_sol]
        )

    if method in ("wrongpose", "known_pose"):
        x0 = np.asarray(ALPHA_INIT, dtype=float)
    elif method == "known_alpha":
        x0 = np.asarray(p_init, dtype=float)
    elif method == "direct":
        x0 = np.concatenate(
            [np.asarray(ALPHA_INIT, dtype=float), np.asarray(p_init, dtype=float)]
        )
    else:
        n_sol = int(solver_subspace["n_vis"])
        x0 = np.concatenate(
            [
                np.asarray(ALPHA_INIT, dtype=float),
                np.zeros(n_sol, dtype=float),
            ]
        )

    result = least_squares(residual, x0, jac=jacobian, **LS_KWARGS)
    status = int(result.status)
    optimality = float(result.optimality)
    success = bool((status > 0) and (optimality < SUCCESS_OPTIMALITY))
    njev = getattr(result, "njev", None)
    njev = int(njev) if njev is not None else None

    alpha_est, p_est = xy(np.asarray(result.x, dtype=float))
    alpha_est = np.asarray(alpha_est, dtype=float)
    p_est = np.asarray(p_est, dtype=float)

    ref_norm = float(np.linalg.norm(d_obs))
    final_residual = float(
        np.linalg.norm(
            _realified_residual(
                forward(alpha_est, p_est, xs, h, k, G_D, Phi, P), d_obs
            )
        )
        / ref_norm
    )
    chi_true = Phi @ ALPHA_TRUE
    chi_est = Phi @ alpha_est
    map_error = float(
        np.linalg.norm(chi_est - chi_true) / np.linalg.norm(chi_true)
    )
    delta_p = p_est - P_TRUE
    pose_error = float(np.linalg.norm(delta_p))
    visible_proj = V_dec @ (V_dec.T @ delta_p)
    pose_error_visible = float(np.linalg.norm(visible_proj))
    pose_error_hidden = float(np.linalg.norm(delta_p - visible_proj))

    t_u = state_witness(
        alpha_est, p_est, 6, xs, h, k, G_D, Phi, P
    )

    return {
        "setting": str(setting["label"]),
        "M": int(P.M),
        "k": float(P.k),
        "aperture": str(setting["aperture"]),
        "aperture_span_rad": float(aperture_span(setting)),
        "method": method if method != "reduced" else f"reduced_r{reduced_r}",
        "r": int(reduced_r) if reduced_r is not None else None,
        "noise_seed": int(seed),
        "success": success,
        "status": status,
        "optimality": optimality,
        "nfev": int(result.nfev),
        "njev": njev,
        "final_residual_ratio": final_residual,
        "alpha_est": [float(v) for v in alpha_est],
        "p_est": [float(v) for v in p_est],
        "map_error": map_error,
        "pose_error": pose_error,
        "pose_error_visible": pose_error_visible,
        "pose_error_hidden": pose_error_hidden,
        "decomposition_basis": "r6",
        "n_vis": int(decomp_subspace["n_vis"]),
        "hidden_rank": int(decomp_subspace["hidden_rank"]),
        "solver_n_vis": (
            int(solver_subspace["n_vis"]) if solver_subspace is not None else None
        ),
        "solver_hidden_rank": (
            int(solver_subspace["hidden_rank"])
            if solver_subspace is not None
            else None
        ),
        "T_U": float(t_u),
        "T_U_rank": 6,
    }


# ---------------------------------------------------------------------------
# Results assembly
# ---------------------------------------------------------------------------
def build_settings_results() -> tuple[dict, dict]:
    """Return (json-like results, memory subspaces per label/r)."""
    xs, h = pixel_grid(BASE_N)
    Phi = basis_matrix(xs, Params(N=BASE_N, s=BASE_P["s"]))
    K = Phi.shape[1]

    gd_cache: dict[float, np.ndarray] = {}
    params: dict[str, Params] = {}
    geometry: dict[str, dict] = {}
    for setting in SETTINGS:
        P = make_params(setting)
        params[setting["label"]] = P
        k = float(P.k)
        if k not in gd_cache:
            gd_cache[k] = green_domain_matrix(xs, h, k)
        geometry[setting["label"]] = dict(
            xs=xs, h=h, k=k, G_D=gd_cache[k], Phi=Phi, P=P
        )

    results: dict = {
        "experiment": "run_e10_settings_sweep",
        "parameters": {
            "N": BASE_N,
            "K": int(K),
            "R_r": BASE_P["R_r"],
            "R_t": BASE_P["R_t"],
            "phi0": BASE_P["phi0"],
            "s": BASE_P["s"],
            "alpha_true": [float(v) for v in ALPHA_TRUE],
            "alpha_init": [float(v) for v in ALPHA_INIT],
            "p_true": [float(v) for v in P_TRUE],
            "p_init": [float(v) for v in P_INIT],
            "SNR_dB": SNR_DB,
            "noise_ratio": NOISE_RATIO,
            "noise_seeds": [int(s) for s in NOISE_SEEDS],
            "success_optimality_max": SUCCESS_OPTIMALITY,
            "settings": [
                {
                    "label": s["label"],
                    "M": int(s["M"]),
                    "k": float(s["k"]),
                    "aperture": s["aperture"],
                    "aperture_span_rad": float(aperture_span(s)),
                    "receiver_angles": (
                        None
                        if _angles(s["aperture"], int(s["M"])) is None
                        else [float(a) for a in _angles(s["aperture"], int(s["M"]))]
                    ),
                }
                for s in SETTINGS
            ],
            "least_squares": dict(LS_KWARGS),
        },
        "settings": {},
    }

    # In-memory full subspace dicts (contain V_vis) per setting/rank.
    spaces: dict[str, dict[int, dict]] = {}

    for setting in SETTINGS:
        label = setting["label"]
        M = int(setting["M"])
        k = geometry[label]["k"]
        xs = geometry[label]["xs"]
        h = geometry[label]["h"]
        G_D = geometry[label]["G_D"]
        Phi = geometry[label]["Phi"]
        P = params[label]

        jb_init = jacobians(ALPHA_INIT, P_INIT, xs, h, k, G_D, Phi, P)
        rc = rank_check(jb_init, M)
        full_lift = unrestricted_lift_check(jb_init)
        ret_lifts = [lift_check(jb_init, r) for r in (4, 6)]
        lift_checks = {
            "unrestricted": full_lift,
            "retained": ret_lifts,
            "b_pert_source": (
                "B_r[:,0] from jacobians(alpha_init, p_init, ...); "
                "receiver-tx complex data-space pose perturbation"
            ),
        }

        sweep_json = []
        sweep_mem = {}
        for r in RANK_SWEEP_RS:
            if r > M:
                break
            sp = visible_pose_subspace(
                ALPHA_INIT, P_INIT, int(r), xs, h, k, G_D, Phi, P
            )
            sweep_mem[int(r)] = sp
            sweep_json.append(subspace_sweep_json(sp, label))
        spaces[label] = sweep_mem

        # d_true / observations for the setting.
        d_true = forward(ALPHA_TRUE, P_TRUE, xs, h, k, G_D, Phi, P)
        obs = {}
        for seed in NOISE_SEEDS:
            raw = make_noise(int(seed), M)
            obs[int(seed)], _ = scaled_observation(d_true, raw)

        runs = []
        for method, rr in (
            ("wrongpose", None),
            ("known_pose", None),
            ("known_alpha", None),
            ("direct", None),
            ("reduced", 4),
            ("reduced", 6),
        ):
            solver_sub = None if rr is None else sweep_mem[int(rr)]
            decomp_sub = sweep_mem[6]
            for seed in NOISE_SEEDS:
                runs.append(
                    run_case(
                        method,
                        rr,
                        int(seed),
                        setting,
                        P,
                        obs[int(seed)],
                        xs,
                        h,
                        k,
                        G_D,
                        Phi,
                        solver_sub,
                        decomp_sub,
                    )
                )
                print(
                    f"[E10] {label:>20s} {runs[-1]['method']:>12s} "
                    f"seed={int(seed)} status={runs[-1]['status']} "
                    f"pose={runs[-1]['pose_error']:.3e} "
                    f"map={runs[-1]['map_error']:.3e}"
                )

        results["settings"][label] = {
            "M": M,
            "k": float(k),
            "aperture": setting["aperture"],
            "aperture_span_rad": float(aperture_span(setting)),
            "receiver_angles": (
                None
                if _angles(setting["aperture"], M) is None
                else [float(a) for a in _angles(setting["aperture"], M)]
            ),
            "rank_checks": rc,
            "lift_checks": lift_checks,
            "subspace_sweep": sweep_json,
            "runs": runs,
        }

    # ------------------------------------------------------------------
    # STEP 5: retained-rank nonlinear ablation on M12_k12_full.
    # ------------------------------------------------------------------
    label = FULL_CIRCLE_M12
    setting = next(s for s in SETTINGS if s["label"] == label)
    M = int(setting["M"])
    k = geometry[label]["k"]
    xs = geometry[label]["xs"]
    h = geometry[label]["h"]
    G_D = geometry[label]["G_D"]
    Phi = geometry[label]["Phi"]
    P = params[label]
    d_true = forward(ALPHA_TRUE, P_TRUE, xs, h, k, G_D, Phi, P)
    obs = {}
    for seed in NOISE_SEEDS:
        raw = make_noise(int(seed), M)
        obs[int(seed)], _ = scaled_observation(d_true, raw)

    rank_sweep_runs = []
    for r in RANK_SWEEP_RS:
        for seed in NOISE_SEEDS:
            rec = run_case(
                "reduced",
                int(r),
                int(seed),
                setting,
                P,
                obs[int(seed)],
                xs,
                h,
                k,
                G_D,
                Phi,
                spaces[label][int(r)],
                spaces[label][6],
            )
            rank_sweep_runs.append(rec)
            print(
                f"[E10 rank-sweep] {rec['method']:>12s} seed={int(seed)} "
                f"status={rec['status']} pose={rec['pose_error']:.3e} "
                f"map={rec['map_error']:.3e}"
            )
    results["rank_sweep_setting"] = label
    results["rank_sweep_runs"] = rank_sweep_runs

    # Counts / seed policy.
    n_total = sum(len(results["settings"][lab]["runs"]) for lab in SETTING_LABELS)
    n_ok = sum(
        1
        for lab in SETTING_LABELS
        for rec in results["settings"][lab]["runs"]
        if rec["success"]
    )
    results["counts"] = {
        "n_setting_method_seed_runs": n_total,
        "n_converged": n_ok,
        "n_rank_sweep_runs": len(rank_sweep_runs),
        "n_rank_sweep_converged": sum(
            1 for rec in rank_sweep_runs if rec["success"]
        ),
    }
    results["note"] = (
        "Reduced runs optimize (alpha, q) with p = p_init + V_vis @ q using "
        "the fixed V_vis computed at (alpha_init, p_init) from the "
        "geometry-lift hiding condition (run_e5_final.visible_pose_subspace). "
        "All pose errors (including direct/known_* and the M12 rank sweep) "
        "are decomposed on the r=6 V_vis at (alpha_init, p_init) so visible "
        "vs hidden components are comparable across methods. T_U is "
        "run_e5.state_witness at retained rank 6 at the final estimate. "
        "Noise: raw=(N(0,1)+1j N(0,1))/sqrt(2), scaled so "
        "||noise||/||d_true||=10^(-SNR_dB/20) exactly (E8/E9 convention). "
        "runs whose subspace hidden_rank>0 cannot move pose along hidden "
        "directions; runs with hidden_rank=0 reduce to the direct pose "
        "parameterization and are expected to match direct."
    )
    return results, spaces


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def _median(records: list[dict], field: str) -> float:
    vals = np.array([float(rec[field]) for rec in records], dtype=float)
    if vals.size == 0:
        return float("nan")
    return float(np.median(vals))


def make_plot(results: dict, png_path: Path) -> None:
    labels = SETTING_LABELS
    full_m = ("M8_k12_full", "M12_k12_full", "M16_k12_full")
    full_colors = {"M8_k12_full": "#1f77b4", "M12_k12_full": "#2ca02c", "M16_k12_full": "#d62728"}

    fig, axes = plt.subplots(3, 2, figsize=(17.0, 17.0))
    plt.rcParams.update({"font.size": 10})

    # (a) G_s singular-value spectra for full-circle M values.
    ax = axes[0, 0]
    for lab in full_m:
        sv = results["settings"][lab]["rank_checks"]["singular_values"]
        ax.semilogy(
            np.arange(1, len(sv) + 1),
            sv,
            "o-",
            ms=4,
            color=full_colors[lab],
            label=f"{lab} (M={results['settings'][lab]['M']})",
        )
    ax.set_xlabel("singular-value index")
    ax.set_ylabel(r"$\sigma_j(G_s)$")
    ax.set_title("(a) Receiver data-operator spectrum (k=12, full circle)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=9)

    # (b) n_vis / hidden_rank vs retained rank r (full-circle M values).
    ax = axes[0, 1]
    for lab in full_m:
        rows = results["settings"][lab]["subspace_sweep"]
        rs = [r["rank"] for r in rows]
        nv = [r["n_vis"] for r in rows]
        hk = [r["hidden_rank"] for r in rows]
        ax.plot(rs, nv, "o-", color=full_colors[lab], label=f"n_vis {lab}")
        ax.plot(
            rs,
            hk,
            "s--",
            ms=4,
            color=full_colors[lab],
            label=f"hidden_rank {lab}",
        )
    ax.set_xticks(list(RANK_SWEEP_RS))
    ax.set_xlabel("retained rank r")
    ax.set_ylabel("pose DOF")
    ax.set_ylim(-0.2, 3.2)
    ax.set_title("(b) Visible / hidden pose DOF vs r (full circle)")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, ncol=2)

    # (c) median pose error: direct vs reduced_r4 vs reduced_r6 per setting.
    ax = axes[1, 0]
    methods_c = ("direct", "reduced_r4", "reduced_r6")
    bar_colors = {"direct": "#1f77b4", "reduced_r4": "#2ca02c", "reduced_r6": "#9467bd"}
    x = np.arange(len(labels))
    width = 0.26
    for j, method in enumerate(methods_c):
        vals = []
        frozen = []
        for lab in labels:
            recs = [
                r
                for r in results["settings"][lab]["runs"]
                if r["method"] == method
            ]
            vals.append(_median(recs, "pose_error"))
            sub_hidden = recs[0]["solver_hidden_rank"]
            frozen.append(bool(sub_hidden and sub_hidden > 0))
        off = (j - 1) * width
        bars = ax.bar(
            x + off,
            np.maximum(vals, 1e-16),
            width,
            color=bar_colors[method],
            label=method,
        )
        for xi, v, fr in zip(x, vals, frozen):
            if fr:
                ax.plot(
                    xi + off,
                    max(v, 1e-16) * 4.0,
                    "*",
                    color="black",
                    ms=13,
                    zorder=5,
                )
    ax.set_yscale("log")
    ax.set_ylim(1e-9, 1e1)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=28, ha="right", fontsize=9)
    ax.set_ylabel("median pose error")
    ax.set_title("(c) Pose error by setting (3 seeds; * = solver subspace hides pose)")
    ax.grid(axis="y", alpha=0.3)
    ax.legend(fontsize=9, loc="upper left")

    # (d) median map error: contrast methods per setting.
    ax = axes[1, 1]
    methods_d = ("wrongpose", "direct", "reduced_r4", "reduced_r6", "known_pose")
    colors_d = {
        "wrongpose": "#d62728",
        "direct": "#1f77b4",
        "reduced_r4": "#2ca02c",
        "reduced_r6": "#9467bd",
        "known_pose": "#ff7f0e",
    }
    for j, method in enumerate(methods_d):
        vals = [
            _median(
                [
                    r
                    for r in results["settings"][lab]["runs"]
                    if r["method"] == method
                ],
                "map_error",
            )
            for lab in labels
        ]
        ax.bar(
            x + (j - 2) * 0.165,
            np.maximum(vals, 1e-16),
            0.155,
            color=colors_d[method],
            label=method,
        )
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=28, ha="right", fontsize=9)
    ax.set_ylabel("median map error")
    ax.set_title("(d) Contrast reconstruction error by setting (3 seeds)")
    ax.grid(axis="y", alpha=0.3)
    ax.legend(fontsize=8, ncol=2)

    # (e) M12 rank sweep.
    ax = axes[2, 0]
    rs_recs = results["rank_sweep_runs"]
    direct_recs = [
        r
        for r in results["settings"][FULL_CIRCLE_M12]["runs"]
        if r["method"] == "direct"
    ]
    direct_med = _median(direct_recs, "pose_error")
    rvals = list(RANK_SWEEP_RS)
    pose_med = [
        _median(
            [r for r in rs_recs if r["method"] == f"reduced_r{r}"], "pose_error"
        )
        for r in rvals
    ]
    vis_med = [
        _median(
            [r for r in rs_recs if r["method"] == f"reduced_r{r}"],
            "pose_error_visible",
        )
        for r in rvals
    ]
    hid_med = [
        _median(
            [r for r in rs_recs if r["method"] == f"reduced_r{r}"],
            "pose_error_hidden",
        )
        for r in rvals
    ]
    ax.axhline(
        direct_med,
        color="#1f77b4",
        ls="--",
        lw=1.4,
        label=f"direct median ({direct_med:.2e})",
    )
    ax.plot(rvals, np.maximum(pose_med, 1e-15), "o-", color="#111111", label="reduced_r pose")
    ax.plot(rvals, np.maximum(vis_med, 1e-15), "s-", color="#2ca02c", label="pose visible (r6 basis)")
    ax.plot(rvals, np.maximum(hid_med, 1e-15), "^--", color="#d62728", label="pose hidden (r6 basis)")
    ax.set_yscale("log")
    ax.set_ylim(1e-15, 1e1)
    ax.set_xticks(rvals)
    ax.set_xlabel("retained rank r (reduced_r solver)")
    ax.set_ylabel("median pose error")
    ax.set_title("(e) M12 full-circle retained-rank sweep (3 seeds)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=9)

    # (f) lift residual ratios across settings.
    ax = axes[2, 1]
    series = {
        "full (unrestricted)": [
            results["settings"][lab]["lift_checks"]["unrestricted"][
                "residual_ratio_full"
            ]
            for lab in labels
        ],
        "retained r4": [
            results["settings"][lab]["lift_checks"]["retained"][0][
                "residual_ratio_ret"
            ]
            for lab in labels
        ],
        "retained r6": [
            results["settings"][lab]["lift_checks"]["retained"][1][
                "residual_ratio_ret"
            ]
            for lab in labels
        ],
    }
    for name, vals in series.items():
        ax.semilogy(
            np.arange(len(labels)),
            np.maximum(vals, 1e-18),
            "o-",
            ms=5,
            label=name,
        )
    ax.set_ylim(1e-18, 1e1)
    ax.set_xticks(np.arange(len(labels)))
    ax.set_xticklabels(labels, rotation=28, ha="right", fontsize=9)
    ax.set_ylabel(r"$\|G_s J - b_{\rm pert}\|/\|b_{\rm pert}\|$")
    ax.set_title("(f) Lift residual of receiver-tx pose perturbation")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=9)

    fig.tight_layout()
    fig.savefig(png_path, dpi=160)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Console summary
# ---------------------------------------------------------------------------
def print_summary(results: dict) -> None:
    print("\n--- E10 settings sweep summary ---")
    for lab in SETTING_LABELS:
        s = results["settings"][lab]
        rk = s["rank_checks"]
        subs = {r["rank"]: r for r in s["subspace_sweep"]}
        nv4 = subs[4]["n_vis"] if 4 in subs else None
        nv6 = subs[6]["n_vis"] if 6 in subs else None
        hk6 = subs[6]["hidden_rank"] if 6 in subs else None
        print(
            f"\n{lab}: M={s['M']} k={s['k']:g} aperture={s['aperture']} "
            f"rank(G_s)={rk['numerical_rank']}/{s['M']} "
            f"n_vis(r4)={nv4} n_vis(r6)={nv6} hidden_rank(r6)={hk6}"
        )
        for method in ("direct", "reduced_r4", "reduced_r6", "known_pose", "wrongpose"):
            recs = [r for r in s["runs"] if r["method"] == method]
            if not recs:
                continue
            ok = sum(int(r["success"]) for r in recs)
            print(
                f"  {method:>12s}: ok={ok}/{len(recs)} "
                f"med_pose={_median(recs, 'pose_error'):.3e} "
                f"med_map={_median(recs, 'map_error'):.3e}"
            )
    rs = results["rank_sweep_runs"]
    print(f"\nM12 rank sweep ({len(rs)} runs):")
    for r in RANK_SWEEP_RS:
        recs = [x for x in rs if x["method"] == f"reduced_r{r}"]
        ok = sum(int(x["success"]) for x in recs)
        print(
            f"  reduced_r{r}: ok={ok}/{len(recs)} "
            f"med_pose={_median(recs, 'pose_error'):.3e} "
            f"med_hidden={_median(recs, 'pose_error_hidden'):.3e} "
            f"med_map={_median(recs, 'map_error'):.3e}"
        )
    counts = results["counts"]
    print(
        f"\nConverged {counts['n_converged']}/{counts['n_setting_method_seed_runs']} "
        f"setting runs and {counts['n_rank_sweep_converged']}/"
        f"{counts['n_rank_sweep_runs']} rank-sweep runs "
        f"(success = status>0 and optimality<1e-7)."
    )


def main() -> dict:
    results, _spaces = build_settings_results()
    json_path = HERE / "results_e10_settings_sweep.json"
    png_path = HERE / "plot_e10_settings_sweep.png"
    with open(json_path, "w") as fh:
        json.dump(results, fh, indent=2)
    make_plot(results, png_path)

    # Verification: JSON round-trip, run counts, PNG loads.
    with open(json_path) as fh:
        loaded = json.load(fh)
    expected_runs = len(SETTINGS) * 6 * len(NOISE_SEEDS)
    assert (
        sum(len(loaded["settings"][lab]["runs"]) for lab in SETTING_LABELS)
        == expected_runs
    ), expected_runs
    assert len(loaded["rank_sweep_runs"]) == len(RANK_SWEEP_RS) * len(NOISE_SEEDS)
    assert loaded["settings"] == results["settings"]
    assert loaded["rank_sweep_runs"] == results["rank_sweep_runs"]
    assert png_path.exists() and png_path.stat().st_size > 0
    try:
        import PIL.Image

        im = PIL.Image.open(png_path)
        im.load()
    except Exception:
        pass  # matplotlib savefig already confirmed a non-empty file

    print_summary(results)
    print("Artifacts:", json_path, png_path)
    return results


if __name__ == "__main__":
    main()
