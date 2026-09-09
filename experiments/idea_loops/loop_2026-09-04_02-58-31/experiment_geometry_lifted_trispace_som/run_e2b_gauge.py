"""Geometry-lifted TriSpace SOM: E2b + gauge/anchor control (run_e2b_gauge.py).

Part 1 -- a true (exact) rank event via a two-parameter receiver-pose scan.
    Scan receiver pose (rx, ry, 0) over an 81 x 81 grid, locate adjacent-pair
    degeneracies of G_s, refine the sigma_4 = sigma_5 crossing, and measure
    hard/soft spectral-projector behaviour through the crossing.

Part 2 -- transmitter gauge/anchor control for the disambiguated regime
    (r=2, K=1).  Rebuilds B_t and the transmitter-state quantities for
    phi0=0.0 and phi0=0.7 and re-runs the E4 transmitter hiding check.

Everything is deterministic (fixed grid, no random numbers).
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import minimize

from geom_som_core import (
    Params,
    pixel_grid,
    contrast_vector,
    green_domain_matrix,
    incident_field,
    incident_derivative,
    state_operator,
    solve_current,
    total_field,
    receiver_positions,
    transmitter_position,
    transmitter_direction,
    data_matrix,
)
from run_e2_e4 import (
    basis_matrix,
    realify_complex_cols,
    realify_real_cols,
    norm_cols,
    rank_svd_tol,
    proj_complement,
)


HERE = Path(__file__).resolve().parent

GRID_N = 81
GRID_LO = -0.8
GRID_HI = 0.8
CROSS_TOL = 1e-10
TAUS = np.array(
    [-0.05, -0.02, -0.01, -0.005, 0.005, 0.01, 0.02, 0.05], dtype=float
)
BT_COLUMNS = ("tx", "ty", "theta")


# ---------------------------------------------------------------------------
# Part 1 helpers
# ---------------------------------------------------------------------------
def _svd_pose(p: np.ndarray, xs: np.ndarray, h: float, P: Params):
    """Sorted SVD (s, V) of G_s at receiver pose (rx, ry, 0)."""
    pose = np.zeros(3, dtype=float)
    pose[: len(p)] = np.asarray(p, dtype=float)
    y, _ = receiver_positions(pose, P)
    U, s, Vh = np.linalg.svd(data_matrix(y, xs, h, P.k), full_matrices=False)
    del U
    return s, Vh.conj().T


def _svals_pose(p: np.ndarray, xs: np.ndarray, h: float, P: Params) -> np.ndarray:
    return _svd_pose(p, xs, h, P)[0]


def _scan_grid(P: Params, xs: np.ndarray, h: float):
    """S[ix, iy, 8]: singular values over the 81x81 pose grid."""
    rxs = np.linspace(GRID_LO, GRID_HI, GRID_N)
    rys = np.linspace(GRID_LO, GRID_HI, GRID_N)
    S = np.empty((GRID_N, GRID_N, P.M), dtype=float)
    for i, rx in enumerate(rxs):
        for j, ry in enumerate(rys):
            S[i, j] = _svals_pose(np.array([rx, ry, 0.0]), xs, h, P)
    return rxs, rys, S


def _refine_pair(
    j: int,
    g: np.ndarray,
    rxs: np.ndarray,
    rys: np.ndarray,
    xs: np.ndarray,
    h: float,
    P: Params,
) -> dict:
    """Refine the exact zero of gap_j = s[j]-s[j+1] by minimizing gap^2.

    numpy/scipy SVD returns singular values in descending order, so every
    gap_j field is pointwise >= 0 and never changes sign.  A genuine
    two-parameter degeneracy is instead an exact zero of the ordered gap;
    we refine the global-minimum grid cell of gap^2 with Nelder-Mead.
    """

    def obj(q: np.ndarray) -> float:
        s = _svals_pose(q, xs, h, P)
        return float((s[j] - s[j + 1]) ** 2)

    gi, gj = np.unravel_index(np.argmin(g), g.shape)
    x0, y0 = rxs[gi], rys[gj]
    res = minimize(
        obj,
        np.array([x0, y0]),
        method="Nelder-Mead",
        options={"xatol": 1e-11, "fatol": 1e-22, "maxiter": 1000},
    )
    q = np.asarray(res.x, dtype=float)
    s = _svals_pose(q, xs, h, P)
    gap = float(abs(s[j] - s[j + 1]))
    inside = bool(
        GRID_LO + 1e-5 <= q[0] <= GRID_HI - 1e-5
        and GRID_LO + 1e-5 <= q[1] <= GRID_HI - 1e-5
    )
    return {
        "refined_gap": gap,
        "point": [float(q[0]), float(q[1])],
        "sigma_pair": [float(s[j]), float(s[j + 1])],
        "inside": inside,
        "grid_gap_start": float(g[gi, gj]),
        "candidates_evaluated": 1,
        "grid_min_cell": [int(gi), int(gj)],
    }


def _scan_crossings(P, xs, h, rxs, rys, S) -> dict:
    """Per adjacent pair i (1..7): ordered-gap min/max and refined zero check."""
    pair_rows = []
    crossing_map = {}
    for j in range(P.M - 1):
        g = S[:, :, j] - S[:, :, j + 1]
        refined = _refine_pair(j, g, rxs, rys, xs, h, P)
        sign_change = bool(float(g.min()) < 0.0 and float(g.max()) > 0.0)
        crossing = bool(
            refined["refined_gap"] < CROSS_TOL and refined["inside"]
        )
        crossing_map[j + 1] = crossing
        pair_rows.append(
            {
                "i": int(j + 1),
        "pair": f"sigma_{j + 1} vs sigma_{j + 2}",
        "grid_min_gap": float(g.min()),
        "grid_max_gap": float(g.max()),
        "ordered_gap_sign_change": sign_change,
        "crossing_found": crossing,
                "refined_gap": refined["refined_gap"],
                "refined_point": refined["point"],
                "sigma_pair_at_refined": refined["sigma_pair"],
                "candidates_evaluated": refined["candidates_evaluated"],
            }
        )
    return {"pair_rows": pair_rows, "crossing_map": crossing_map}


def _choose_pair(crossing_map: dict, pair_rows: list[dict]) -> int:
    for i in (4, 3, 5):
        if crossing_map.get(i, False):
            return i
    fallback = [r for r in pair_rows if r["crossing_found"]]
    if fallback:
        return int(min(fallback, key=lambda r: r["refined_gap"])["i"])
    raise RuntimeError("No exact adjacent-pair degeneracy found in the scan box")


def _part1_projectors(
    P: Params, xs: np.ndarray, h: float, pc: np.ndarray, i: int
) -> dict:
    j = i - 1
    s_c, _ = _svd_pose(pc, xs, h, P)
    lam = 0.5 * (s_c[j] + s_c[j + 1])

    _, V0 = _svd_pose(pc, xs, h, P)
    P_hard0 = V0[:, :i] @ V0[:, :i].conj().T
    w0 = s_c**2 / (s_c**2 + lam**2)
    P_soft0 = V0 @ np.diag(w0) @ V0.conj().T

    proj_at = {}
    for tau in TAUS:
        s, V = _svd_pose(pc + tau * np.array([1.0, 0.0]), xs, h, P)
        Ph = V[:, :i] @ V[:, :i].conj().T
        w = s**2 / (s**2 + lam**2)
        Ps = V @ np.diag(w) @ V.conj().T
        proj_at[float(tau)] = (Ph, Ps)

    Ph_m, Ps_m = proj_at[-0.02]
    Ph_p, Ps_p = proj_at[0.02]
    d_hard = [
        float(np.linalg.norm(proj_at[tau][0] - P_hard0, "fro")) for tau in TAUS
    ]
    d_soft = [
        float(np.linalg.norm(proj_at[tau][1] - P_soft0, "fro")) for tau in TAUS
    ]
    return {
        "crossing_sigma_pair": [float(s_c[j]), float(s_c[j + 1])],
        "sigma": float(0.5 * (s_c[j] + s_c[j + 1])),
        "lambda_soft": float(lam),
        "jump_hard": float(np.linalg.norm(Ph_p - Ph_m, "fro")),
        "jump_soft": float(np.linalg.norm(Ps_p - Ps_m, "fro")),
        "tau": [float(t) for t in TAUS],
        "d_hard": d_hard,
        "d_soft": d_soft,
    }


def run_part1(P: Params, xs: np.ndarray, h: float) -> dict:
    rxs, rys, S = _scan_grid(P, xs, h)
    scan = _scan_crossings(P, xs, h, rxs, rys, S)
    pair_rows = scan["pair_rows"]
    chosen_i = _choose_pair(scan["crossing_map"], pair_rows)
    chosen_row = next(r for r in pair_rows if r["i"] == chosen_i)
    pc = np.asarray(chosen_row["refined_point"], dtype=float)
    j = chosen_i - 1

    s_check = _svals_pose(pc, xs, h, P)
    verification_gap = float(abs(s_check[j] - s_check[j + 1]))
    assert verification_gap < 1e-10, (
        f"chosen crossing not exact: gap={verification_gap:.3e}"
    )

    proj = _part1_projectors(P, xs, h, pc, chosen_i)
    gap_chosen = S[:, :, j] - S[:, :, j + 1]
    return {
        "grid": {
            "rx_lo": GRID_LO,
            "rx_hi": GRID_HI,
            "ry_lo": GRID_LO,
            "ry_hi": GRID_HI,
            "n": GRID_N,
        },
        "pairs_with_crossings": [
            int(i) for i, ok in scan["crossing_map"].items() if ok
        ],
        "pair_rows": pair_rows,
        "chosen_pair_i": int(chosen_i),
        "crossing": {
            "rx": float(pc[0]),
            "ry": float(pc[1]),
            "sigma": proj["sigma"],
            "sigma_pair": proj["crossing_sigma_pair"],
            "i": int(chosen_i),
        },
        "verification_gap_at_crossing": verification_gap,
        "verification_lt_1e-10": bool(verification_gap < 1e-10),
        "jump_hard": proj["jump_hard"],
        "jump_soft": proj["jump_soft"],
        "tau": proj["tau"],
        "d_hard": proj["d_hard"],
        "d_soft": proj["d_soft"],
        "lambda_soft": proj["lambda_soft"],
        "note": (
            "numpy.linalg.svd sorts singular values descending, so the ordered "
            "gap s[i-1]-s[i] is pointwise >= 0 and never changes sign; each "
            "pair was therefore tested for an exact zero by refining the "
            "global-minimum grid cell of gap^2 with Nelder-Mead (crossing "
            "found iff refined gap < 1e-10). The chosen pair is the sigma_4 "
            "vs sigma_5 degeneracy at the retained rank r=4. Projector "
            "tau-walk uses direction d=(1,0), i.e. rx about the refined "
            "(rx*, ry*)."
        ),
        "_gap_chosen_grid": gap_chosen,
        "_rxs": rxs,
        "_rys": rys,
    }


def plot_part1(part1: dict, path: Path) -> None:
    rxs = part1["_rxs"]
    rys = part1["_rys"]
    gap = part1["_gap_chosen_grid"]
    X, Y = np.meshgrid(rxs, rys, indexing="ij")
    loggap = np.log10(np.maximum(gap, 1e-15))
    rx_star = part1["crossing"]["rx"]
    ry_star = part1["crossing"]["ry"]

    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.2))
    ax = axes[0]
    im = ax.pcolormesh(X, Y, loggap, shading="auto", cmap="viridis")
    levels = np.log10([1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3])
    cs = ax.contour(X, Y, loggap, levels=levels, colors="white", linewidths=0.6)
    ax.clabel(cs, inline=True, fontsize=7, fmt="%.0f")
    ax.plot(rx_star, ry_star, marker="*", ms=16, mfc="red", mec="k")
    ax.text(
        rx_star,
        ry_star,
        f" $\\sigma_4=\\sigma_5$ crossing\n ({rx_star:.3f}, {ry_star:.3f})",
        color="red",
        fontsize=9,
        va="center",
    )
    fig.colorbar(im, ax=ax, label=r"$\log_{10}(s_4-s_5)$")
    ax.set_xlabel("receiver $r_x$")
    ax.set_ylabel("receiver $r_y$")
    ax.set_title(
        "E2b Part 1: exact $\\sigma_4=\\sigma_5$ gap field over pose space"
    )

    ax = axes[1]
    tau = np.asarray(part1["tau"])
    d_hard = np.asarray(part1["d_hard"])
    d_soft = np.asarray(part1["d_soft"])
    ax.axvline(0.0, color="k", ls=":", lw=1.0)
    ax.plot(tau, d_hard, "o-", ms=4, lw=1.1, label="hard projector $d_{\\mathrm{hard}}$")
    ax.plot(tau, d_soft, "s-", ms=4, lw=1.1, label="soft projector $d_{\\mathrm{soft}}$")
    ax.annotate(
        f"$\\|P_h(+0.02)-P_h(-0.02)\\|_F$ = {part1['jump_hard']:.3f}",
        xy=(0.0, d_hard.max() * 0.65),
        fontsize=9,
    )
    ax.set_xlabel(r"$\tau$ through crossing $(r_x^*+\tau, r_y^*, 0)$")
    ax.set_ylabel("Frobenius distance from crossing-point projector")
    ax.set_title("E2b Part 1: projector distance through rank event")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Part 2 helpers
# ---------------------------------------------------------------------------
def _build_bt_geometry(P: Params, xs, h, chi, G_D) -> dict:
    """Transmitter pose reference quantities at Params.phi0 (p0 = 0)."""
    p0 = np.zeros(3)
    y0, _ = receiver_positions(p0, P)
    G_s = data_matrix(y0, xs, h, P.k)
    t, t_local = transmitter_position(p0, P)
    u_inc = incident_field(xs, t, P.k)
    J = solve_current(chi, u_inc, G_D)
    A = state_operator(chi, G_D)
    u = total_field(u_inc, J, G_D)
    J_chi = np.linalg.solve(A, np.diag(u))  # Np x Np current response

    w_dir = [
        np.array([1.0, 0.0]),
        np.array([0.0, 1.0]),
        transmitter_direction("theta", t_local),
    ]
    B_t = np.empty((P.M, 3), dtype=complex)
    for c, w in enumerate(w_dir):
        du_inc = incident_derivative(xs, t, w, P.k)
        J_p = np.linalg.solve(A, chi * du_inc)
        B_t[:, c] = G_s @ J_p
    return {
        "G_s": G_s,
        "t": t,
        "t_local": t_local,
        "w_dir": w_dir,
        "u_inc": u_inc,
        "J": J,
        "u": u,
        "A": A,
        "J_chi": J_chi,
        "B_t": B_t,
    }


def _bt_stats(B_t: np.ndarray) -> dict:
    """Complex rank + singular values of column-normalized realification."""
    complex_rank = int(rank_svd_tol(B_t))
    complex_sv = np.linalg.svd(B_t, compute_uv=False)
    raw_real = np.vstack([B_t.real, B_t.imag])  # 2M x 3
    raw_real_rank = int(rank_svd_tol(raw_real))
    Btot_n = norm_cols(raw_real)
    sv_norm = np.linalg.svd(Btot_n, compute_uv=False)
    # theta column (index 2) residual in span of tx/ty columns.
    c_fit, *_ = np.linalg.lstsq(B_t[:, :2], B_t[:, 2], rcond=None)
    theta_resid = float(
        np.linalg.norm(B_t[:, 2] - B_t[:, :2] @ c_fit)
        / np.linalg.norm(B_t[:, 2])
    )
    return {
        "columns": list(BT_COLUMNS),
        "matrix_rank_complex": complex_rank,
        "complex_singular_values": [float(v) for v in complex_sv],
        "rank_raw_realified": raw_real_rank,
        "singular_values_realified_column_normed": [float(v) for v in sv_norm],
        "theta_col_residual_in_span_txty": theta_resid,
    }


def _transmitter_hide_check(P: Params, geo: dict) -> dict:
    """E4 disambiguated (r=2, K=1) transmitter-only hiding check."""
    r, K = 2, 1
    xs_geo = geo["xs"]
    _, _, Vh = np.linalg.svd(geo["G_s"], full_matrices=False)
    V_r = Vh.conj().T[:, :r]
    Q = geo["G_s"] @ V_r
    Phi = basis_matrix(xs_geo, P, K)
    W = geo["G_s"] @ geo["J_chi"] @ Phi
    Hn = norm_cols(
        np.hstack([realify_complex_cols(Q), realify_real_cols(W)])
    )
    rank_H = int(rank_svd_tol(Hn))
    P_perp = proj_complement(Hn)
    Btot_n = norm_cols(np.vstack([geo["B_t"].real, geo["B_t"].imag]))
    M_mat = P_perp @ Btot_n
    hid_sv = np.linalg.svd(M_mat, compute_uv=False)
    return {
        "r": r,
        "K": K,
        "rank_H": rank_H,
        "hid_sv": [float(v) for v in hid_sv],
        "min_hid_sv": float(np.min(hid_sv)),
        "hid_sv_all_gt_1e-9": bool(np.all(hid_sv > 1e-9)),
    }


def run_part2(P: Params, xs: np.ndarray, h: float, chi: np.ndarray, G_D) -> dict:
    geo = _build_bt_geometry(P, xs, h, chi, G_D)
    geo["xs"] = xs
    stats = _bt_stats(geo["B_t"])
    hide = _transmitter_hide_check(P, geo)
    return {"B_t": stats, "transmitter_hide_check": hide, "geometry": geo}


# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------
def main() -> dict:
    P = Params()
    xs, h = pixel_grid(P.N)
    chi = contrast_vector(xs, P)
    G_D = green_domain_matrix(xs, h, P.k)

    part1 = run_part1(P, xs, h)
    part1_pub = {
        key: value
        for key, value in part1.items()
        if not key.startswith("_")
    }
    plot_part1(part1, HERE / "plot_e2b_rank_event.png")

    part2 = {}
    part2_summaries = []
    for phi in (0.0, 0.7):
        P_phi = Params(phi0=phi)
        res = run_part2(P_phi, xs, h, chi, G_D)
        part2[f"phi0_{phi:g}"] = {
            "B_t": res["B_t"],
            "transmitter_hide_check": res["transmitter_hide_check"],
        }
        part2_summaries.append((phi, res))

    part2["phi0_0.7"]["expectation_check"] = {
        "B_t_rank3_expected": False,
        "B_t_rank3_observed": (
            part2["phi0_0.7"]["B_t"]["matrix_rank_complex"] == 3
        ),
        "hid_sv_all_gt_1e-9_expected": True,
        "hid_sv_all_gt_1e-9_observed": part2["phi0_0.7"][
            "transmitter_hide_check"
        ]["hid_sv_all_gt_1e-9"],
        "note": (
            "For a scalar point transmitter the theta-pose derivative of the "
            "source position is w_theta = R_t(-sin(phi0), cos(phi0)) = "
            "-R_t sin(phi0) e_x + R_t cos(phi0) e_y, i.e. exactly a linear "
            "combination of the tx/ty translation derivatives. Since all "
            "derivatives are evaluated at one reference pose, B_t has at most "
            "two independent complex columns for every phi0: phi0=0.7 removes "
            "the special collinearity theta ~ ty but cannot create a third "
            "transmitter data direction. A rank-3 transmitter Jacobian needs "
            "an orientation-dependent source (pattern/array) or co-moving "
            "receiver motion."
        ),
    }

    results = {
        "parameters": {
            "N": P.N,
            "k": P.k,
            "M": P.M,
            "R_r": P.R_r,
            "R_t": P.R_t,
            "s": P.s,
            "reference_contrast_coeffs": list(P.coeffs),
            "grid": {"rx_lo": GRID_LO, "rx_hi": GRID_HI, "n": GRID_N},
        },
        "Part1": part1_pub,
        "Part2": part2,
    }

    json_path = HERE / "results_e2b_gauge.json"
    png_path = HERE / "plot_e2b_rank_event.png"
    with open(json_path, "w") as fh:
        json.dump(results, fh, indent=2)

    # Deterministic verification: JSON round-trips and the plot exists.
    with open(json_path) as fh:
        loaded = json.load(fh)
    assert loaded["Part1"] == results["Part1"]
    assert loaded["Part1"]["verification_gap_at_crossing"] < 1e-10
    assert loaded["Part1"]["chosen_pair_i"] == 4
    assert len(loaded["Part1"]["d_hard"]) == len(TAUS)
    assert len(loaded["Part2"]["phi0_0.7"]["B_t"]["singular_values_realified_column_normed"]) == 3
    assert png_path.exists() and png_path.stat().st_size > 0

    # ---- concise console summary -------------------------------------------
    print("--- E2b Part 1: two-parameter rank event ---")
    print(
        "pairs_with_crossings:",
        loaded["Part1"]["pairs_with_crossings"],
    )
    print(
        "chosen pair i =",
        loaded["Part1"]["chosen_pair_i"],
        "| crossing = (",
        f"{loaded['Part1']['crossing']['rx']:.6f},",
        f"{loaded['Part1']['crossing']['ry']:.6f})",
        "| sigma =",
        f"{loaded['Part1']['crossing']['sigma']:.6e}",
        "| refined gap =",
        f"{loaded['Part1']['verification_gap_at_crossing']:.3e}",
    )
    print(
        "jump_hard =",
        f"{loaded['Part1']['jump_hard']:.6f}",
        "jump_soft =",
        f"{loaded['Part1']['jump_soft']:.6f}",
    )
    for row in loaded["Part1"]["pair_rows"]:
        print(
            f"  i={row['i']} gap_min={row['grid_min_gap']:.3e} "
            f"sign_change={row['ordered_gap_sign_change']} "
            f"crossing={row['crossing_found']} "
            f"refined={row['refined_gap']:.3e}"
        )
    print("\n--- E2b Part 2: transmitter gauge / anchor ---")
    for phi, res in part2_summaries:
        b = res["B_t"]
        print(
            f"phi0={phi:.1f}: complex rank {b['matrix_rank_complex']} "
            "| realified column-normed svs = "
            + ", ".join(f"{v:.3e}" for v in b["singular_values_realified_column_normed"])
            + f" | theta-in-span residual {b['theta_col_residual_in_span_txty']:.2e}"
        )
    hid = part2["phi0_0.7"]["transmitter_hide_check"]
    print(
        "phi0=0.7 disambiguated transmitter: rank_H =",
        hid["rank_H"],
        "| hid_sv =",
        [f"{v:.3e}" for v in hid["hid_sv"]],
    )
    print(
        "Part 2 note: phi0=0.7 does not create a third transmitter direction "
        "(rank 2, hid_sv ~1e-16); parent spec intervention recommended."
    )
    print("\nArtifacts:", json_path, png_path)
    return results


if __name__ == "__main__":
    main()
