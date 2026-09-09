"""Geometry-lifted TriSpace SOM: experiments E2 and E4 (run_e2_e4.py).

E2 -- spectral-coordinate drift of the receiver data operator and a deliberate
rank event.  Receiver-only pose rx in [-0.6, 0.6]; hard rank-r truncation
projectors are compared with fixed-threshold soft spectral projectors.

E4 -- exact finite-dimensional hiding after realification.  In a 16-dimensional
realified data space we test whether receiver/transmitter/comoving pose data
perturbations can be explained by the retained SOM current columns plus real
contrast modes (confounded r=6,K=3 vs disambiguated r=2,K=1).

Everything is deterministic (no nonlinear optimization, no random numbers).
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from geom_som_core import (
    Params,
    pixel_grid,
    contrast_vector,
    green_domain_matrix,
    incident_field,
    state_operator,
    solve_current,
    total_field,
    receiver_positions,
    receiver_directions,
    transmitter_position,
    transmitter_direction,
    data_matrix,
    receiver_data_derivative,
    incident_derivative,
)


HERE = Path(__file__).resolve().parent
DIRS = ("rx", "ry", "theta")
TOL = 1e-10
HIDDEN_THRESHOLD = 1e-9


# ---------------------------------------------------------------------------
# Shared helpers (E4)
# ---------------------------------------------------------------------------
def basis_matrix(xs: np.ndarray, P: Params, K: int) -> np.ndarray:
    """Np x K real Gaussian basis.  K=1/2/3 as in the experiment spec."""
    centers = ((-0.15, 0.10), (0.20, -0.10), (0.00, 0.00))[:K]
    Phi = np.empty((xs.shape[0], K), dtype=float)
    for j, center in enumerate(centers):
        delta = xs - np.asarray(center, dtype=float)
        Phi[:, j] = np.exp(
            -np.einsum("ij,ij->i", delta, delta) / (2.0 * P.s**2)
        )
    return Phi


def realify_complex_cols(Q: np.ndarray) -> np.ndarray:
    """2M x 2r real matrix [[Re Q, -Im Q], [Im Q, Re Q]]."""
    return np.block([[Q.real, -Q.imag], [Q.imag, Q.real]])


def realify_real_cols(W: np.ndarray) -> np.ndarray:
    """2M x K real matrix [[Re W], [Im W]] (real coefficients)."""
    return np.vstack([W.real, W.imag])


def norm_cols(M: np.ndarray) -> np.ndarray:
    """Column-normalize with unit 2-norm; leave exactly-zero columns as zero."""
    out = np.array(M, dtype=float, copy=True)
    norms = np.linalg.norm(out, axis=0)
    nz = norms > 0.0
    out[:, nz] /= norms[nz]
    return out


def rank_svd_tol(X: np.ndarray, tol: float = TOL) -> int:
    """Rank by the number of singular values above tol."""
    return int(np.count_nonzero(np.linalg.svd(X, compute_uv=False) > tol))


def proj_complement(Xn: np.ndarray, tol: float = TOL) -> np.ndarray:
    """P_perp = I - U_r U_r^H for column-normalized real Xn (2M x cols)."""
    sv = np.linalg.svd(Xn, compute_uv=False)
    keep = int(np.count_nonzero(sv > tol))
    U = np.linalg.svd(Xn, full_matrices=True)[0]
    U_r = U[:, :keep]
    return np.eye(Xn.shape[0]) - U_r @ U_r.T


def realified_unit(b: np.ndarray) -> np.ndarray:
    """Unit-norm [Re b; Im b] for a complex M-vector b."""
    x = np.concatenate([b.real, b.imag])
    n = np.linalg.norm(x)
    return x / n if n > 0.0 else x


# ---------------------------------------------------------------------------
# E2
# ---------------------------------------------------------------------------
def _receiver_svd_at(rx: float, xs: np.ndarray, h: float, P: Params):
    y, _ = receiver_positions(np.array([rx, 0.0, 0.0]), P)
    G_s = data_matrix(y, xs, h, P.k)
    U, s, Vh = np.linalg.svd(G_s, full_matrices=False)
    V = Vh.conj().T
    return G_s, U, s, V


def _projectors_at(rx: float, xs, h, P: Params, r: int, lam: float):
    """Return hard/soft projectors P_hard, P_soft for receiver pose rx."""
    _, _, s, V = _receiver_svd_at(rx, xs, h, P)
    P_hard = V[:, :r] @ V[:, :r].conj().T
    w = s**2 / (s**2 + lam**2)
    P_soft = V @ np.diag(w) @ V.conj().T
    return P_hard, P_soft


def run_e2(P: Params) -> dict:
    xs, h = pixel_grid(P.N)
    rxs = np.linspace(-0.6, 0.6, 121)
    S_all = np.empty((len(rxs), P.M), dtype=float)
    Vs = []
    for i, rx in enumerate(rxs):
        _, _, s, V = _receiver_svd_at(rx, xs, h, P)
        S_all[i] = s
        Vs.append(V)

    # ---- rank-event scan: strict sign change of s[r-1] - s[r] ---------------
    scan = {}
    crossing_r = None
    crossing_rx = None
    for r_try in (4, 3, 5):
        d = S_all[:, r_try - 1] - S_all[:, r_try]
        hits = []
        for i in range(len(rxs) - 1):
            if d[i] * d[i + 1] < 0.0:
                rx_i = rxs[i]
                rx_j = rxs[i + 1]
                # linear zero interpolation of s[r-1](rx) - s[r](rx)
                denom = d[i + 1] - d[i]
                rx_c = rx_i - d[i] * (rx_j - rx_i) / denom
                hits.append(float(rx_c))
        scan[f"r{r_try}_crossings"] = hits
        if crossing_r is None and hits:
            crossing_r = r_try
            crossing_rx = hits[0]

    # No strict crossing in the requested window -> closest avoided crossing.
    r_use = crossing_r if crossing_r is not None else 4
    gap = np.abs(S_all[:, r_use - 1] - S_all[:, r_use])
    min_idx = int(np.argmin(gap))
    min_gap = float(gap[min_idx])
    min_gap_rx = float(rxs[min_idx])
    event_rx = crossing_rx if crossing_rx is not None else min_gap_rx
    event_kind = (
        "sign_change_crossing"
        if crossing_rx is not None
        else "min_gap_avoided_crossing_no_sign_change"
    )

    # ---- projectors and drift metrics for the selected retained rank --------
    _, _, s0, V0 = _receiver_svd_at(0.0, xs, h, P)
    lam = 0.5 * (s0[r_use - 1] + s0[r_use])
    P_hard0 = V0[:, :r_use] @ V0[:, :r_use].conj().T
    w0 = s0**2 / (s0**2 + lam**2)
    P_soft0 = V0 @ np.diag(w0) @ V0.conj().T

    d_hard = np.empty_like(rxs)
    d_soft = np.empty_like(rxs)
    leakage_hard = np.empty_like(rxs)
    leakage_soft = np.empty_like(rxs)
    for i, rx in enumerate(rxs):
        V = Vs[i]
        s = S_all[i]
        P_hard = V[:, :r_use] @ V[:, :r_use].conj().T
        w = s**2 / (s**2 + lam**2)
        P_soft = V @ np.diag(w) @ V.conj().T
        d_hard[i] = np.linalg.norm(P_hard - P_hard0, "fro")
        d_soft[i] = np.linalg.norm(P_soft - P_soft0, "fro")
        leakage_hard[i] = np.linalg.norm(
            (np.eye(P_hard0.shape[0]) - P_hard0) @ P_hard, "fro"
        )
        leakage_soft[i] = np.linalg.norm(
            (np.eye(P_soft0.shape[0]) - P_soft0) @ P_soft, "fro"
        )

    # ---- jumps evaluated around the event point (offset +- 0.02) ------------
    P_hard_m, P_soft_m = _projectors_at(event_rx - 0.02, xs, h, P, r_use, lam)
    P_hard_p, P_soft_p = _projectors_at(event_rx + 0.02, xs, h, P, r_use, lam)
    _, _, s_m, _ = _receiver_svd_at(event_rx - 0.02, xs, h, P)
    _, _, s_p, _ = _receiver_svd_at(event_rx + 0.02, xs, h, P)
    jump_hard = float(np.linalg.norm(P_hard_p - P_hard_m, "fro"))
    jump_soft = float(np.linalg.norm(P_soft_p - P_soft_m, "fro"))

    # ---- downsample 121 -> 31 points ---------------------------------------
    ds = np.linspace(0, len(rxs) - 1, 31).astype(int)
    result = {
        "r_used": int(r_use),
        "r_scan": [4, 3, 5],
        "rank_event_scan": scan,
        "crossing_detected": crossing_rx is not None,
        "crossing_r": int(crossing_r) if crossing_r is not None else None,
        "crossing_rx": float(crossing_rx) if crossing_rx is not None else None,
        "rank_event_kind": event_kind,
        "rank_event_rx": float(event_rx),
        "min_gap": min_gap,
        "min_gap_rx": min_gap_rx,
        "singular_values_rx_c_minus_0p02": [float(x) for x in s_m],
        "singular_values_rx_c_plus_0p02": [float(x) for x in s_p],
        "singular_values_at_reference_rx0": [float(x) for x in s0],
        "lambda_soft_fixed_at_rx0": float(lam),
        "jump_hard": jump_hard,
        "jump_soft": jump_soft,
        "note": (
            "No strict sign-change crossing of the retained-rank gap was found "
            "in [-0.6, 0.6] for r=4, r=3, or r=5 (see rank_event_scan). "
            "rank_event_rx is therefore the minimum-gap avoided-crossing point "
            "for r=4, and the jump metrics evaluate the projector discontinuity "
            "around that point."
        ),
        "downsampled": {
            "rx": [float(x) for x in rxs[ds]],
            "d_hard": [float(x) for x in d_hard[ds]],
            "d_soft": [float(x) for x in d_soft[ds]],
            "leakage_hard": [float(x) for x in leakage_hard[ds]],
            "leakage_soft": [float(x) for x in leakage_soft[ds]],
            "sigma_3": [float(x) for x in S_all[ds, 2]],
            "sigma_4": [float(x) for x in S_all[ds, 3]],
            "sigma_5": [float(x) for x in S_all[ds, 4]],
            "sigma_6": [float(x) for x in S_all[ds, 5]],
        },
    }
    return result, rxs, S_all, d_hard, d_soft, event_rx


def plot_e2(
    rxs: np.ndarray,
    S_all: np.ndarray,
    d_hard: np.ndarray,
    d_soft: np.ndarray,
    event_rx: float,
    path: Path,
) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(8.2, 8.0), sharex=True)
    for j, label in ((2, "$\\sigma_3$"), (3, "$\\sigma_4$"), (4, "$\\sigma_5$"), (5, "$\\sigma_6$")):
        axes[0].plot(rxs, S_all[:, j], lw=1.3, label=label)
    axes[0].axvline(event_rx, color="k", ls=":", lw=1.0)
    axes[0].text(
        event_rx,
        axes[0].get_ylim()[1] * 0.985,
        f"rank event rx={event_rx:.3f}",
        ha="center",
        va="top",
        fontsize=9,
    )
    axes[0].set_ylabel("singular value")
    axes[0].set_title("E2: receiver data-operator singular values vs rx")
    axes[0].legend(loc="best", fontsize=9, ncol=2)
    axes[0].grid(alpha=0.3)

    axes[1].plot(rxs, d_hard, "o-", ms=3, lw=0.9, label="hard rank projector $d_{\\mathrm{hard}}$")
    axes[1].plot(rxs, d_soft, "s-", ms=3, lw=0.9, label="soft filter projector $d_{\\mathrm{soft}}$")
    axes[1].axvline(event_rx, color="k", ls=":", lw=1.0)
    axes[1].set_xlabel("receiver pose $r_x$")
    axes[1].set_ylabel("Frobenius distance from rx=0 projector")
    axes[1].set_title("E2: projector drift from reference pose")
    axes[1].legend(fontsize=9)
    axes[1].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


# ---------------------------------------------------------------------------
# E4
# ---------------------------------------------------------------------------
def _build_e4_geometry(P: Params):
    """All fixed reference-pose quantities used by both regimes."""
    xs, h = pixel_grid(P.N)
    G_D = green_domain_matrix(xs, h, P.k)
    chi = contrast_vector(xs, P)
    p0 = np.zeros(3)
    y0, y_local = receiver_positions(p0, P)
    t0, t_local = transmitter_position(p0, P)
    G_s = data_matrix(y0, xs, h, P.k)
    u_inc = incident_field(xs, t0, P.k)
    J = solve_current(chi, u_inc, G_D)
    u = total_field(u_inc, J, G_D)
    A = state_operator(chi, G_D)
    J_chi = np.linalg.solve(A, np.diag(u))  # Np x Np current response

    B_r = np.empty((P.M, 3), dtype=complex)
    B_t = np.empty((P.M, 3), dtype=complex)
    for j, direction in enumerate(DIRS):
        v = receiver_directions(direction, y_local)
        B_r[:, j] = receiver_data_derivative(y0, v, xs, h, P.k) @ J
        w = transmitter_direction(direction, t_local)
        du_inc_dc = incident_derivative(xs, t0, w, P.k)
        J_p = np.linalg.solve(A, chi * du_inc_dc)
        B_t[:, j] = G_s @ J_p
    return dict(
        xs=xs, h=h, chi=chi, G_D=G_D, A=A, u=u, J=J, G_s=G_s,
        J_chi=J_chi, B_r=B_r, B_t=B_t,
    )


def _regime_blocks(P, geo, r, K):
    """Return Q, V_r, W, Hn, P_perp, W_real, P_perp_W for a regime."""
    _, _, Vh = np.linalg.svd(geo["G_s"], full_matrices=False)
    V_r = Vh.conj().T[:, :r]
    Q = geo["G_s"] @ V_r
    Phi = basis_matrix(geo["xs"], P, K)
    W = geo["G_s"] @ geo["J_chi"] @ Phi
    Hn = norm_cols(np.hstack([realify_complex_cols(Q), realify_real_cols(W)]))
    rank_H = rank_svd_tol(Hn)
    P_perp = proj_complement(Hn)
    W_real = norm_cols(realify_real_cols(W))
    P_perp_W = proj_complement(W_real)
    return {
        "V_r": V_r,
        "Q": Q,
        "W": W,
        "Hn": Hn,
        "rank_H": rank_H,
        "P_perp": P_perp,
        "W_real": W_real,
        "P_perp_W": P_perp_W,
    }


def _run_e4_regime(P, geo, blocks, r, K, regime_name):
    B_r = geo["B_r"]
    B_t = geo["B_t"]
    motions = {"receiver": B_r, "transmitter": B_t, "comoving": B_r + B_t}
    per_motion = {}
    for motion_name, B_total in motions.items():
        Btot_raw = np.vstack([B_total.real, B_total.imag])  # 2M x 3
        Btot_n = norm_cols(Btot_raw)
        colnorms = np.linalg.norm(Btot_raw, axis=0)
        M_mat = blocks["P_perp"] @ Btot_n
        _, hid_sv, Vh_m = np.linalg.svd(M_mat, full_matrices=False)

        direction_rows = []
        for c, direction in enumerate(DIRS):
            b = B_total[:, c]
            bn = realified_unit(b)
            residual_data_hidden = float(np.linalg.norm(blocks["P_perp"] @ bn))
            rank_inc = int(
                rank_svd_tol(np.hstack([blocks["Hn"], bn[:, None]]))
                - blocks["rank_H"]
            )
            hidden_flag = bool(residual_data_hidden < HIDDEN_THRESHOLD)

            # Receiver-part lift witnesses (computed from B_r[:, c]).
            b_r_c = B_r[:, c]
            c_r = np.linalg.pinv(blocks["Q"]) @ b_r_c
            dJ_r = blocks["V_r"] @ c_r
            R_U = b_r_c - blocks["Q"] @ c_r
            r_data = float(
                np.linalg.norm(R_U) / np.linalg.norm(b_r_c)
                if np.linalg.norm(b_r_c) > 0.0
                else np.nan
            )
            delta_chi_implied = (geo["A"] @ dJ_r) / geo["u"]
            T_U = float(
                np.linalg.norm(delta_chi_implied.imag)
                / np.linalg.norm(delta_chi_implied)
            )
            resid_state_only = float(
                np.linalg.norm(blocks["P_perp_W"] @ bn)
            )
            direction_rows.append(
                {
                    "direction": direction,
                    "residual_data_hidden": residual_data_hidden,
                    "rank_inc": rank_inc,
                    "hidden_flag": hidden_flag,
                    "r_data": r_data,
                    "T_U": T_U,
                    "resid_state_only": resid_state_only,
                }
            )

        per_motion[motion_name] = {
            "rank_H": int(blocks["rank_H"]),
            "hid_sv": [float(x) for x in hid_sv],
            "directions": direction_rows,
            "min_hid_sv": float(np.min(hid_sv)),
            "B_column_norms": [float(x) for x in colnorms],
        }
    return per_motion


def _demo_directions(P, geo, blocks):
    """Confounded receiver: visible/hidden direction demonstration.

    M_mat = P_perp @ Btot_n is formed with column-normalized motion columns,
    so its right singular vectors live in the normalized-coordinate basis.
    Mapping back to physical pose coefficients requires division by the column
    norms; that is the direction for which the residual is genuinely ~0.
    Both variants are stored below.
    """
    B_total = geo["B_r"]
    Btot_raw = np.vstack([B_total.real, B_total.imag])
    Btot_n = norm_cols(Btot_raw)
    colnorms = np.linalg.norm(Btot_raw, axis=0)
    M_mat = blocks["P_perp"] @ Btot_n
    _, hid_sv, Vh_m = np.linalg.svd(M_mat, full_matrices=False)
    v_visible = Vh_m[0]
    v_hidden = Vh_m[-1]

    rows = {}
    for label, v in (("visible", v_visible), ("hidden", v_hidden)):
        # Physical-coefficient direction consistent with the normalized SVD.
        delta_p = v / colnorms
        delta_p = delta_p / np.linalg.norm(delta_p)
        b = B_total @ delta_p
        bn = realified_unit(b)
        resid_vs_Sdata = float(np.linalg.norm(blocks["P_perp"] @ bn))
        resid_vs_Sstate = float(np.linalg.norm(blocks["P_perp_W"] @ bn))
        c_r = np.linalg.pinv(blocks["Q"]) @ (B_total @ delta_p)
        dJ_r = blocks["V_r"] @ c_r
        b_norm = np.linalg.norm(B_total @ delta_p)
        r_data = float(
            np.linalg.norm(B_total @ delta_p - blocks["Q"] @ c_r) / b_norm
            if b_norm > 0.0
            else np.nan
        )
        delta_chi_implied = (geo["A"] @ dJ_r) / geo["u"]
        T_U = float(
            np.linalg.norm(delta_chi_implied.imag)
            / np.linalg.norm(delta_chi_implied)
        )
        rows[label] = {
            "kind": f"confounded_receiver_{label}_direction",
            "delta_p_physical": [float(x) for x in delta_p],
            "delta_p_normalized_svd": [float(x) for x in v],
            "resid_vs_Sdata": resid_vs_Sdata,
            "resid_vs_Sstate": resid_vs_Sstate,
            "r_data": r_data,
            "T_U": T_U,
            "M_mat_singular_value": float(hid_sv[0] if label == "visible" else hid_sv[-1]),
        }
    return rows


def run_e4(P: Params) -> tuple[dict, dict]:
    geo = _build_e4_geometry(P)
    regimes = (("confounded", 6, 3), ("disambiguated", 2, 1))
    out = {}
    demos = {}
    for regime_name, r, K in regimes:
        blocks = _regime_blocks(P, geo, r, K)
        motion_data = _run_e4_regime(P, geo, blocks, r, K, regime_name)
        out[regime_name] = {
            "r": r,
            "K": K,
            "rank_H": blocks["rank_H"],
            "motion": motion_data,
        }
        if regime_name == "confounded":
            demos = _demo_directions(P, geo, blocks)
    return out, demos


def plot_e4(out: dict, demos: dict, path: Path) -> None:
    # Panel 1: per-coordinate residuals grouped by regime and motion type.
    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.6))
    ax = axes[0]
    group_keys = [
        ("confounded", "receiver"),
        ("disambiguated", "receiver"),
        ("confounded", "transmitter"),
        ("disambiguated", "transmitter"),
        ("confounded", "comoving"),
        ("disambiguated", "comoving"),
    ]
    x = np.arange(len(group_keys))
    offsets = (-0.27, 0.0, 0.27)
    colors = ("#1f77b4", "#ff7f0e", "#2ca02c")
    for gi, (regime, motion) in enumerate(group_keys):
        rows = out[regime]["motion"][motion]["directions"]
        for ci, (off, col, row) in enumerate(zip(offsets, colors, rows)):
            val = max(float(row["residual_data_hidden"]), 1e-16)
            ax.bar(
                gi + off, val, width=0.25, color=col,
                label=(DIRS[ci] if gi == 0 else None),
            )
    ax.axhline(HIDDEN_THRESHOLD, color="k", ls="--", lw=1.0)
    ax.text(
        5.45, HIDDEN_THRESHOLD, "1e-9", ha="right", va="bottom", fontsize=8
    )
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(
        ["conf\nrecv", "dis\nrecv", "conf\ntrans", "dis\ntrans", "conf\ncomov", "dis\ncomov"],
        fontsize=8,
    )
    ax.set_ylabel("residual_data_hidden")
    ax.set_title("E4: per-coordinate pose data residual vs col([Q, W])")
    ax.grid(axis="y", alpha=0.3)
    ax.legend(fontsize=8, ncol=3, loc="upper left")

    # Panel 2: confounded receiver hidden/visible directions.
    ax = axes[1]
    labels = ["hidden", "visible"]
    xpos = np.arange(len(labels))
    for j, metric in enumerate(("resid_vs_Sdata", "resid_vs_Sstate")):
        vals = [
            max(demos[label][metric], 1e-16) for label in labels
        ]
        ax.bar(
            xpos + (j - 0.5) * 0.34,
            vals,
            0.32,
            label="data-only $\\mathrm{col}(Q)$" if metric == "resid_vs_Sdata" else "joint data+state $\\mathrm{col}([Q,W])$",
        )
    ax.set_yscale("log")
    ax.set_xticks(xpos)
    ax.set_xticklabels(labels)
    ax.set_ylabel("residual")
    ax.set_ylim(1e-16, 2.0)
    for label in labels:
        tu = demos[label]["T_U"]
        rd = demos[label]["r_data"]
        ax.text(
            xpos[labels.index(label)],
            0.95,
            f"$T_U$={tu:.3f}\n$r_{{\\mathrm{{data}}}}$={rd:.3f}",
            ha="center",
            va="top",
            fontsize=8,
        )
    ax.set_title("E4: hidden/visible receiver directions\n(residuals; T_U, r_data annotated)")
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> dict:
    P = Params()
    e2_result, rxs, S_all, d_hard, d_soft, event_rx = run_e2(P)
    e4_result, e4_demos = run_e4(P)
    e4_result["note"] = (
        "r_data/T_U are receiver-part witnesses computed from B_r[:, c] in "
        "every motion-type entry (they are receiver-lift statements, not "
        "transmitter-lift statements). With phi0=0 the single transmitter has "
        "B_t[:, ry] proportional to B_t[:, theta], so transmitter-motion M_mat "
        "is exactly column-degenerate; its zero hid_sv is a parameter-space "
        "degeneracy, not an additional hiding direction."
    )

    results = {
        "parameters": {
            "N": P.N,
            "k": P.k,
            "M": P.M,
            "R_r": P.R_r,
            "R_t": P.R_t,
            "phi0": P.phi0,
            "s": P.s,
            "reference_contrast_coeffs": list(P.coeffs),
        },
        "E2": e2_result,
        "E4": e4_result,
        "E4_demonstration": {
            "confounded_receiver_hidden_direction": e4_demos["hidden"],
            "confounded_receiver_visible_direction": e4_demos["visible"],
            "note": (
                "Demo directions are the right singular vectors of "
                "M_mat=P_perp@Btot_n mapped back to physical pose coefficients "
                "(division by Btot_n column norms), so the hidden residual is "
                "the exact ~0 physical-direction residual. Unscaled SVD vectors "
                "are stored as delta_p_normalized_svd."
            ),
        },
    }

    json_path = HERE / "results_e2_e4.json"
    e2_png = HERE / "plot_e2_rank_event.png"
    e4_png = HERE / "plot_e4_hiding.png"
    with open(json_path, "w") as fh:
        json.dump(results, fh, indent=2)
    plot_e2(rxs, S_all, d_hard, d_soft, event_rx, e2_png)
    plot_e4(e4_result, e4_demos, e4_png)

    # Deterministic verification: reload and check sizes/shape invariants.
    with open(json_path) as fh:
        loaded = json.load(fh)
    assert loaded["E2"] == results["E2"]
    assert len(loaded["E2"]["downsampled"]["rx"]) == 31
    assert loaded["E2"]["rank_event_rx"] == results["E2"]["rank_event_rx"]
    assert {"confounded", "disambiguated"} <= set(loaded["E4"])
    for regime in ("confounded", "disambiguated"):
        for motion in loaded["E4"][regime]["motion"]:
            entry = loaded["E4"][regime]["motion"][motion]
            assert len(entry["hid_sv"]) == 3
            assert len(entry["directions"]) == 3
    for png in (e2_png, e4_png):
        assert png.exists() and png.stat().st_size > 0

    print("--- E2 ---")
    e2 = loaded["E2"]
    print(
        f"r_used={e2['r_used']} crossing_detected={e2['crossing_detected']} "
        f"event_kind={e2['rank_event_kind']} rank_event_rx={e2['rank_event_rx']:.6f}"
    )
    print(
        f"min_gap={e2['min_gap']:.3e} at rx={e2['min_gap_rx']:.6f} | "
        f"jump_hard={e2['jump_hard']:.6f} jump_soft={e2['jump_soft']:.6f}"
    )
    print("rank_event_scan:", e2["rank_event_scan"])

    print("\n--- E4 ---")
    for regime in ("confounded", "disambiguated"):
        rinfo = loaded["E4"][regime]
        print(f"{regime}: r={rinfo['r']} K={rinfo['K']} rank_H={rinfo['rank_H']}")
        for motion in ("receiver", "transmitter", "comoving"):
            entry = rinfo["motion"][motion]
            print(
                f"  {motion:>11s}: min hid_sv = {min(entry['hid_sv']):.3e} "
                f"hid_sv = {[f'{v:.2e}' for v in entry['hid_sv']]}"
            )
    dem = loaded["E4_demonstration"]
    for key in ("confounded_receiver_hidden_direction", "confounded_receiver_visible_direction"):
        d = dem[key]
        print(
            f"{key}: resid_vs_Sdata={d['resid_vs_Sdata']:.3e} "
            f"resid_vs_Sstate={d['resid_vs_Sstate']:.3e} "
            f"r_data={d['r_data']:.3e} T_U={d['T_U']:.3f}"
        )
    print("\nArtifacts:", json_path, e2_png, e4_png)
    return results


if __name__ == "__main__":
    main()
