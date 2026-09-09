"""Geometry-lifted TriSpace SOM: 2D scalar Helmholtz numerical core (E1 & E3).

Self-contained implementation of the forward model and analytic derivatives
specified for the experiment "Geometry-Lifted TriSpace SOM for Self-Calibrating
Full-Wave Inverse Scattering".

Convention note (state equation ordering)
-----------------------------------------
With A(chi) = I - diag(chi) @ G_D, the internally consistent induced/contrast
current is the Lippmann-Schwinger contrast-source solution
        J = A^{-1} (chi .* u_inc),    i.e.  A J = chi .* u_inc,
which is algebraically J = chi .* u with the physical total field
u = u_inc + G_D J (so J = chi (u_inc + G_D J) <=> A J = chi u_inc).
Under this convention the analytic transmitter-pose derivative stated in the
spec, J_p = solve(A, chi .* du_inc_dc), is exactly the derivative of the
current used in the finite-difference checks, so E3 converges at O(eps^2).
(Using instead u = solve(A, u_inc) with J = chi*u would make the transmitter /
co-moving FD checks inconsistent with J_p at the O(chi G) level.)
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy.special import hankel1

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception as exc:  # pragma: no cover - plotting is optional
    plt = None


HERE = Path(__file__).resolve().parent


# ---------------------------------------------------------------------------
# Model parameters
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Params:
    N: int = 24
    k: float = 12.0
    M: int = 8
    R_r: float = 1.6
    R_t: float = 2.0
    phi0: float = 0.0
    r: int = 4  # retained SOM rank
    s: float = 0.12  # Gaussian basis width
    centers: tuple = field(default_factory=lambda: ((-0.15, 0.10), (0.20, -0.10)))
    coeffs: tuple = field(default_factory=lambda: (1.5, 2.0))
    eps_list: tuple = field(default_factory=lambda: (1e-2, 1e-3, 1e-4, 1e-5))
    receiver_angles: tuple | None = None

    @property
    def h(self) -> float:
        return 1.0 / self.N

    @property
    def Np(self) -> int:
        return self.N * self.N


def _rotation(theta: float) -> np.ndarray:
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, -s], [s, c]])


# ---------------------------------------------------------------------------
# Grid, contrast, Green operators
# ---------------------------------------------------------------------------
def pixel_grid(N: int) -> tuple[np.ndarray, float]:
    """Return (centers, h): N*N x 2 pixel centers and side h in [-0.5, 0.5]^2."""
    h = 1.0 / N
    one_d = -0.5 + (np.arange(N) + 0.5) * h
    X, Y = np.meshgrid(one_d, one_d, indexing="xy")
    centers = np.stack((X.ravel(), Y.ravel()), axis=1)  # row-major: i*N + j
    return centers, h


def contrast_vector(xs: np.ndarray, P: Params) -> np.ndarray:
    """chi(x) = sum_j a_j exp(-|x-c_j|^2 / (2 s^2)), real length-Np array."""
    chi = np.zeros(xs.shape[0], dtype=float)
    for coeff, center in zip(P.coeffs, P.centers):
        delta = xs - np.asarray(center, dtype=float)
        chi += coeff * np.exp(-np.einsum("ij,ij->i", delta, delta) / (2.0 * P.s**2))
    return chi


def _green_kernel(r: np.ndarray, k: float, order: int = 0) -> np.ndarray:
    return hankel1(order, k * r)


def green_domain_matrix(xs: np.ndarray, h: float, k: float) -> np.ndarray:
    """G_D (Np x Np): 2D free-space Helmholtz Green on the pixel grid."""
    Np = xs.shape[0]
    dx = xs[:, None, :] - xs[None, :, :]  # (Np, Np, 2)
    r = np.sqrt(np.einsum("ijk,ijk->ij", dx, dx))
    G = np.zeros((Np, Np), dtype=complex)
    off = ~np.eye(Np, dtype=bool)
    G[off] = h**2 * (1j / 4.0) * _green_kernel(r[off], k, 0)

    a = h / np.sqrt(np.pi)  # equivalent disk radius (same area as a pixel)
    diag = (1j * np.pi * a / (2.0 * k)) * hankel1(1, k * a) - 1.0 / k**2
    G[np.diag_indices(Np)] = diag
    return G


def diag_selfcell_diff(G_D: np.ndarray, h: float) -> float:
    """Max |diag(G_D) - 1j*h^2/4| (the stated quadrature self-cell check)."""
    diag = np.diag(G_D)
    return float(np.max(np.abs(diag - 1j * h**2 / 4.0)))


def incident_field(xs: np.ndarray, t: np.ndarray, k: float) -> np.ndarray:
    """u_inc at pixels from a 2D point source at t."""
    delta = xs - t[None, :]
    R = np.sqrt(np.einsum("ij,ij->i", delta, delta))
    return (1j / 4.0) * _green_kernel(R, k, 0)


def state_operator(chi: np.ndarray, G_D: np.ndarray) -> np.ndarray:
    return np.eye(chi.size) - np.diag(chi) @ G_D


def solve_current(chi: np.ndarray, u_inc: np.ndarray, G_D: np.ndarray) -> np.ndarray:
    """J = A^{-1} (chi .* u_inc), A = I - diag(chi) G_D (contrast current)."""
    A = state_operator(chi, G_D)
    return np.linalg.solve(A, chi * u_inc)


def total_field(u_inc: np.ndarray, J: np.ndarray, G_D: np.ndarray) -> np.ndarray:
    """Physical total field u = u_inc + G_D J  (J = chi*u by the state equation)."""
    return u_inc + G_D @ J


# ---------------------------------------------------------------------------
# Geometry / pose
# ---------------------------------------------------------------------------
def receiver_positions(p_r: np.ndarray, P: Params) -> tuple[np.ndarray, np.ndarray]:
    """Receiver points y_m (Lrx2) and their local base points y_local (Lrx2).

    If P.receiver_angles is not None, the Lr = len(P.receiver_angles)
    body-frame angles given there are used (e.g. limited-aperture arcs).
    Otherwise the historical full-circle grid 2*pi*arange(P.M)/P.M is kept
    exactly (Lr = P.M).
    """
    if P.receiver_angles is not None:
        th0 = np.asarray(P.receiver_angles, dtype=float)
    else:
        th0 = 2.0 * np.pi * np.arange(P.M) / P.M
    y_local = P.R_r * np.stack((np.cos(th0), np.sin(th0)), axis=1)
    R = _rotation(p_r[2])
    y = y_local @ R.T + p_r[:2]
    return y, y_local


def receiver_directions(direction: str, y_local: np.ndarray) -> np.ndarray:
    """Receiver derivative vectors v_m (Mx2) for a pose coordinate."""
    direction = direction.lower()
    if direction in ("rx", "x", "tx"):
        return np.tile(np.array([1.0, 0.0]), (y_local.shape[0], 1))
    if direction in ("ry", "y"):
        return np.tile(np.array([0.0, 1.0]), (y_local.shape[0], 1))
    if direction in ("theta", "th"):
        return np.stack((-y_local[:, 1], y_local[:, 0]), axis=1)
    raise ValueError(f"unknown receiver direction {direction!r}")


def transmitter_position(p_t: np.ndarray, P: Params) -> tuple[np.ndarray, np.ndarray]:
    """Transmitter point t and its local base position t_local."""
    t_local = P.R_t * np.array([np.cos(P.phi0), np.sin(P.phi0)])
    t = _rotation(p_t[2]) @ t_local + p_t[:2]
    return t, t_local


def transmitter_direction(direction: str, t_local: np.ndarray) -> np.ndarray:
    """Transmitter derivative vector w_t (2,) for a pose coordinate."""
    direction = direction.lower()
    if direction in ("rx", "x", "tx"):
        return np.array([1.0, 0.0])
    if direction in ("ry", "y"):
        return np.array([0.0, 1.0])
    if direction in ("theta", "th"):
        return np.array([-t_local[1], t_local[0]])
    raise ValueError(f"unknown transmitter direction {direction!r}")


# ---------------------------------------------------------------------------
# Data operators and analytic derivatives
# ---------------------------------------------------------------------------
def data_matrix(y: np.ndarray, xs: np.ndarray, h: float, k: float) -> np.ndarray:
    """G_s (M x Np): receiver Green/data matrix (no self terms)."""
    delta = y[:, None, :] - xs[None, :, :]  # (M, Np, 2)
    r = np.sqrt(np.einsum("ijk,ijk->ij", delta, delta))
    return h**2 * (1j / 4.0) * _green_kernel(r, k, 0)


def receiver_data_derivative(
    y: np.ndarray, v: np.ndarray, xs: np.ndarray, h: float, k: float
) -> np.ndarray:
    """dGs/dc (M x Np) for receiver direction vectors v (M x 2).

    Entry: h^2 (i/4) (-H1(k r)) k ((y_m - x_n).v_m)/r  with r=|y_m - x_n|.
    """
    delta = y[:, None, :] - xs[None, :, :]  # (M, Np, 2)
    r = np.sqrt(np.einsum("ijk,ijk->ij", delta, delta))
    assert np.all(r > 0.0), "receiver pixel distance must be positive"
    dot = np.einsum("ijk,ik->ij", delta, v)
    return (
        h**2
        * (1j / 4.0)
        * (-_green_kernel(r, k, 1))
        * k
        * dot
        / r
    )


def incident_derivative(
    xs: np.ndarray, t: np.ndarray, w: np.ndarray, k: float
) -> np.ndarray:
    """du_inc/dc at pixels for transmitter motion direction w.

    Entry: (i/4) H1(k R) k ((x_n - t).w)/R  with R = |x_n - t|.
    """
    delta = xs - t[None, :]
    R = np.sqrt(np.einsum("ij,ij->i", delta, delta))
    assert np.all(R > 0.0), "transmitter pixel distance must be positive"
    dot = delta @ w
    return (1j / 4.0) * _green_kernel(R, k, 1) * k * dot / R


# ---------------------------------------------------------------------------
# SOM retained-space lift
# ---------------------------------------------------------------------------
def svd_data_operator(G_s: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compact SVD G_s = U s Vh; returns U, s (descending), V (columns)."""
    U, s, Vh = np.linalg.svd(G_s, full_matrices=False)
    V = Vh.conj().T
    return U, s, V


def som_lift(
    G_s: np.ndarray, V_r: np.ndarray, b: np.ndarray
) -> dict:
    """Lift data vector b into span(V_r): c=pinv(Q)b, dJ=V_r c, D=Q c, R=b-D."""
    Q = G_s @ V_r
    c_r = np.linalg.pinv(Q) @ b
    dJ_r = V_r @ c_r
    D_U = Q @ c_r
    R_U = b - D_U
    return dict(Q=Q, c_r=c_r, dJ_r=dJ_r, D_U=D_U, R_U=R_U)


def unrestricted_lift(G_s: np.ndarray, b: np.ndarray) -> dict:
    c_full = np.linalg.pinv(G_s) @ b
    dJ_full = c_full
    D_full = G_s @ c_full
    R_full = b - D_full
    return dict(c_full=c_full, dJ_full=dJ_full, D_full=D_full, R_full=R_full)


# ---------------------------------------------------------------------------
# Experiment E1
# ---------------------------------------------------------------------------
def run_E1(P: Params) -> dict:
    xs, h = pixel_grid(P.N)
    chi = contrast_vector(xs, P)
    G_D = green_domain_matrix(xs, h, P.k)
    t0, _ = transmitter_position(np.zeros(3), P)
    u_inc = incident_field(xs, t0, P.k)
    J = solve_current(chi, u_inc, G_D)

    y0, y_local0 = receiver_positions(np.zeros(3), P)
    G_s = data_matrix(y0, xs, h, P.k)

    rank_Gs = int(np.linalg.matrix_rank(G_s))
    U, s, V = svd_data_operator(G_s)
    del U, V  # singular values only needed for E1 report

    v_rx = receiver_directions("rx", y_local0)
    dGs_drx = receiver_data_derivative(y0, v_rx, xs, h, P.k)
    b = dGs_drx @ J  # data-space derivative for receiver x-translation

    unf = unrestricted_lift(G_s, b)
    V_r = svd_data_operator(G_s)[2][:, : P.r]
    ret = som_lift(G_s, V_r, b)

    diag_diff = diag_selfcell_diff(G_D, h)
    return {
        "parameters": dict(
            N=P.N,
            k=P.k,
            M=P.M,
            R_r=P.R_r,
            R_t=P.R_t,
            r=P.r,
            singular_values_s=s.real.tolist(),
            rank_G_s=rank_Gs,
            diag_selfcell_max_abs_diff_vs_ih2o4=diag_diff,
        ),
        "E1": dict(
            unrestricted_norm_R_full_over_norm_b=float(
                np.linalg.norm(unf["R_full"]) / np.linalg.norm(b)
            ),
            unrestricted_norm_dJ_full=float(np.linalg.norm(unf["dJ_full"])),
            retained_r4_norm_R_U_over_norm_b=float(
                np.linalg.norm(ret["R_U"]) / np.linalg.norm(b)
            ),
            retained_r4_norm_dJ_r=float(np.linalg.norm(ret["dJ_r"])),
            retained_r4_orthogonality_residual=float(
                np.linalg.norm(ret["Q"].conj().T @ ret["R_U"])
                / np.linalg.norm(ret["R_U"])
            ),
        ),
    }


# ---------------------------------------------------------------------------
# Experiment E3
# ---------------------------------------------------------------------------
def _state_current_at_transmitter(p_t: np.ndarray, P: Params, xs, G_D, chi):
    t, _ = transmitter_position(p_t, P)
    return solve_current(chi, incident_field(xs, t, P.k), G_D)


def _receiver_fd(p_r_plus, p_r_minus, J, xs, h, k, P, eps):
    yp, _ = receiver_positions(p_r_plus, P)
    ym, _ = receiver_positions(p_r_minus, P)
    Gp, Gm = data_matrix(yp, xs, h, k), data_matrix(ym, xs, h, k)
    return (Gp @ J - Gm @ J) / (2.0 * eps)


def _transmitter_fd(
    p_t_plus, p_t_minus, xs, h, k, G_D, chi, P, eps
) -> tuple[np.ndarray, np.ndarray]:
    Jp = _state_current_at_transmitter(p_t_plus, P, xs, G_D, chi)
    Jm = _state_current_at_transmitter(p_t_minus, P, xs, G_D, chi)
    Jp_fd = (Jp - Jm) / (2.0 * eps)
    y0, _ = receiver_positions(np.zeros(3), P)
    G_s = data_matrix(y0, xs, h, k)
    return Jp_fd, G_s @ Jp_fd


def run_E3(P: Params) -> list[dict]:
    xs, h = pixel_grid(P.N)
    chi = contrast_vector(xs, P)
    G_D = green_domain_matrix(xs, h, P.k)
    A = state_operator(chi, G_D)

    p0 = np.zeros(3)
    t0, t_local0 = transmitter_position(p0, P)
    u_inc0 = incident_field(xs, t0, P.k)
    J0 = solve_current(chi, u_inc0, G_D)
    y0, y_local0 = receiver_positions(p0, P)
    G_s0 = data_matrix(y0, xs, h, P.k)

    rows: list[dict] = []
    for direction in ("rx", "theta"):
        v = receiver_directions(direction, y_local0)
        w = transmitter_direction(direction, t_local0)

        # Analytic pieces at the reference pose.
        B_r = receiver_data_derivative(y0, v, xs, h, P.k) @ J0
        du_inc = incident_derivative(xs, t0, w, P.k)
        J_p = np.linalg.solve(A, chi * du_inc)
        B_t = G_s0 @ J_p
        B_total = B_r + B_t

        for eps in P.eps_list:
            coord = {"rx": 0, "theta": 2}[direction]
            e = np.zeros(3)
            e[coord] = eps

            # Receiver-only FD (transmitter and J fixed).
            B_r_fd = _receiver_fd(p0 + e, p0 - e, J0, xs, h, P.k, P, eps)
            rel_rec = float(np.linalg.norm(B_r_fd - B_r) / np.linalg.norm(B_r))

            # Transmitter-only FD (receivers fixed, full current re-solve).
            J_p_fd, B_t_fd_data = _transmitter_fd(
                p0 + e, p0 - e, xs, h, P.k, G_D, chi, P, eps
            )
            rel_tra = float(np.linalg.norm(J_p_fd - J_p) / np.linalg.norm(J_p))
            rel_tra_data = float(
                np.linalg.norm(B_t_fd_data - B_t) / np.linalg.norm(B_t)
            )

            # Co-moving FD (receivers and transmitter moved together).
            yp, _ = receiver_positions(p0 + e, P)
            ym, _ = receiver_positions(p0 - e, P)
            tp, _ = transmitter_position(p0 + e, P)
            tm, _ = transmitter_position(p0 - e, P)
            Jp = solve_current(chi, incident_field(xs, tp, P.k), G_D)
            Jm = solve_current(chi, incident_field(xs, tm, P.k), G_D)
            d_p = data_matrix(yp, xs, h, P.k) @ Jp
            d_m = data_matrix(ym, xs, h, P.k) @ Jm
            d_total_fd = (d_p - d_m) / (2.0 * eps)
            rel_com = float(
                np.linalg.norm(d_total_fd - B_total) / np.linalg.norm(B_total)
            )

            rows.append(
                {
                    "direction": direction,
                    "type": "receiver",
                    "eps": eps,
                    "rel_err": rel_rec,
                }
            )
            rows.append(
                {
                    "direction": direction,
                    "type": "transmitter",
                    "eps": eps,
                    "rel_err": rel_tra,
                }
            )
            rows.append(
                {
                    "direction": direction,
                    "type": "comoving",
                    "eps": eps,
                    "rel_err": rel_com,
                }
            )

            print(
                f"[E3] dir={direction:>5s} eps={eps:.0e} | "
                f"recv={rel_rec:.3e} tra(J)={rel_tra:.3e} "
                f"tra(data)={rel_tra_data:.3e} comov={rel_com:.3e}"
            )
    return rows


def make_fd_plot(P: Params, rows: list[dict], path: Path) -> None:
    if plt is None:  # pragma: no cover
        raise RuntimeError("matplotlib not available")
    fig, ax = plt.subplots(figsize=(7.0, 5.2))
    for kind, style in (
        ("receiver", "o-"),
        ("transmitter", "s-"),
        ("comoving", "^-"),
    ):
        pts = [
            (r["eps"], r["rel_err"])
            for r in rows
            if r["type"] == kind and r["direction"] == "rx"
        ]
        eps = np.array([p[0] for p in pts])
        err = np.array([p[1] for p in pts])
        ax.loglog(eps, err, style, label=kind)
    # Reference O(eps^2) line anchored at the eps=1e-2 error of receiver FD.
    ref = next(r["rel_err"] for r in rows if r["type"] == "receiver" and r["eps"] == 1e-2)
    eps_ref = np.array([1e-2, 1e-4])
    ax.loglog(eps_ref, ref * (eps_ref / 1e-2) ** 2, "k--", lw=0.8, label=r"$O(\varepsilon^2)$")
    ax.set_xlabel("step $\\varepsilon$")
    ax.set_ylabel("relative error vs analytic derivative")
    ax.set_title("E3: centered finite-difference vs analytic derivatives (rx)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------
def main() -> dict:
    P = Params()
    out_path = HERE / "results_e1_e3.json"
    png_path = HERE / "plot_fd_errors.png"

    e1 = run_E1(P)
    e3_rows = run_E3(P)

    results = {
        "parameters": e1["parameters"],
        "E1": e1["E1"],
        "E3": e3_rows,
    }
    with open(out_path, "w") as fh:
        json.dump(results, fh, indent=2)

    make_fd_plot(P, e3_rows, png_path)

    # Verification: reload JSON and confirm plot exists.
    with open(out_path) as fh:
        loaded = json.load(fh)
    assert loaded["E1"] == results["E1"], "JSON reload mismatch"
    assert len(loaded["E3"]) == 2 * 3 * len(P.eps_list)
    assert png_path.exists() and png_path.stat().st_size > 0

    print("\n--- Summary ---")
    print("rank(G_s) =", loaded["parameters"]["rank_G_s"])
    print(
        "max |diag(G_D) - i h^2/4| =",
        loaded["parameters"]["diag_selfcell_max_abs_diff_vs_ih2o4"],
    )
    print("singular values:", np.round(loaded["parameters"]["singular_values_s"], 6))
    print("E1 R_full/b:", loaded["E1"]["unrestricted_norm_R_full_over_norm_b"])
    print("E1 R_U/b:", loaded["E1"]["retained_r4_norm_R_U_over_norm_b"])
    print(
        "E1 orthogonality residual:",
        loaded["E1"]["retained_r4_orthogonality_residual"],
    )
    print("Artifacts:", out_path, png_path)
    return results


if __name__ == "__main__":
    main()
