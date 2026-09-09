"""Independent audit of the retained-range pose-rank transition.

This script does not modify the locked Agentic-AI-Scientist experiment.  It
recomputes the pose-visible rank with a projector-codimension check that is
well defined when the nuisance span saturates the realified data space.

The historical experiment used only ``s_i > 1e-8*s_1`` on the singular
values of ``B_red = P_perp B``.  If ``P_perp`` is numerically zero, all three
singular values are roundoff-sized but have ordinary relative ratios, so that
rule incorrectly returns rank three.  Here a full nuisance-space rank implies
``P_perp = 0`` and hence visible rank zero exactly.  Away from saturation we
combine a relative cutoff with a backward-error-scaled absolute cutoff.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
EXP = (
    ROOT
    / "experiments/idea_loops/loop_2026-09-04_02-58-31"
    / "experiment_geometry_lifted_trispace_som"
)
sys.path.insert(0, str(EXP))

from geom_som_core import Params, green_domain_matrix, pixel_grid  # noqa: E402
from run_e2_e4 import (  # noqa: E402
    norm_cols,
    realify_complex_cols,
    realify_real_cols,
)
from run_e5 import basis_matrix, forward, jacobians  # noqa: E402
from run_e10_settings_sweep import (  # noqa: E402
    ALPHA_INIT,
    ALPHA_TRUE,
    P_INIT,
    P_TRUE,
    make_noise,
    scaled_observation,
)
from run_e10b_m12_m16_highrank import run_case  # noqa: E402


H_RANK_TOL = 1e-10
B_REL_TOL = 1e-8
B_BACKWARD_FACTOR = 100.0
SEEDS = (0, 1, 2)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def params_for(M: int) -> Params:
    return Params(
        N=16,
        k=12.0,
        M=M,
        R_r=1.6,
        R_t=2.0,
        phi0=0.7,
        s=0.12,
    )


def corrected_subspace(
    alpha: np.ndarray,
    pose: np.ndarray,
    r: int,
    xs: np.ndarray,
    h: float,
    G_D: np.ndarray,
    Phi: np.ndarray,
    P: Params,
) -> dict:
    """Return a rank-safe visible pose basis and complete diagnostics."""
    jb = jacobians(alpha, pose, xs, h, P.k, G_D, Phi, P)
    _, _, Vh = np.linalg.svd(jb["G_s"], full_matrices=False)
    V_r = Vh.conj().T[:, :r]
    Q = jb["G_s"] @ V_r
    Hn = norm_cols(
        np.hstack(
            [
                realify_complex_cols(Q),
                realify_real_cols(jb["J_alpha"]),
            ]
        )
    )

    U_h, s_h, _ = np.linalg.svd(Hn, full_matrices=True)
    rank_h = int(np.count_nonzero(s_h > H_RANK_TOL))
    data_dim = int(Hn.shape[0])
    codim = data_dim - rank_h
    U_range = U_h[:, :rank_h]
    P_perp = np.eye(data_dim) - U_range @ U_range.T
    B = jb["B_real"]
    B_red = P_perp @ B
    _, s_b, Vh_b = np.linalg.svd(B_red, full_matrices=False)

    b_scale = float(np.linalg.norm(B, ord=2))
    backward_tol = (
        B_BACKWARD_FACTOR
        * np.finfo(float).eps
        * max(B_red.shape)
        * max(b_scale, np.finfo(float).tiny)
    )
    relative_tol = B_REL_TOL * float(s_b[0]) if s_b.size else 0.0
    rank_tol = max(relative_tol, backward_tol)

    # If the nuisance range is the whole data space, its orthogonal
    # complement is exactly {0}; roundoff in P_perp cannot create visibility.
    if codim == 0:
        n_vis = 0
    else:
        n_vis = min(
            codim,
            int(np.count_nonzero(np.asarray(s_b) > rank_tol)),
            B.shape[1],
        )
    V_vis = Vh_b[:n_vis].T
    hidden_rank = int(B.shape[1] - n_vis)

    return {
        "rank": int(r),
        "rank_H": rank_h,
        "data_dim": data_dim,
        "nuisance_codimension": codim,
        "nuisance_column_count": int(Hn.shape[1]),
        "n_vis": int(n_vis),
        "hidden_rank": hidden_rank,
        "hid_sv": [float(v) for v in s_b],
        "rank_tol": float(rank_tol),
        "relative_tol": float(relative_tol),
        "backward_tol": float(backward_tol),
        "projector_norm": float(np.linalg.norm(P_perp, ord=2)),
        "B_norm": b_scale,
        "V_vis": V_vis,
        "hidden_directions": [
            [float(x) for x in Vh_b[i]] for i in range(n_vis, B.shape[1])
        ],
    }


def summarize_cases(records: list[dict]) -> list[dict]:
    groups: dict[tuple, list[dict]] = {}
    for row in records:
        key = (row["setting"], row["method"], row["regime"])
        groups.setdefault(key, []).append(row)
    out = []
    for (setting, method, regime), rows in sorted(groups.items()):
        out.append(
            {
                "setting": setting,
                "method": method,
                "regime": regime,
                "n": len(rows),
                "successes": int(sum(bool(x["success"]) for x in rows)),
                "median_pose_error": float(
                    np.median([x["pose_error"] for x in rows])
                ),
                "median_map_error": float(
                    np.median([x["map_error"] for x in rows])
                ),
                "median_residual": float(
                    np.median([x["final_residual_ratio"] for x in rows])
                ),
            }
        )
    return out


def write_markdown(payload: dict) -> None:
    lines = [
        "# Scale-aware audit of the retained-range pose rank",
        "",
        "This is an independent parent-level recomputation. It preserves the "
        "locked auto-research outputs and corrects only their scientific "
        "interpretation at nuisance-range saturation.",
        "",
        "## Rank table",
        "",
        "| M | r | rank(H) / 2M | codim | visible pose rank | hidden pose rank | ||P_perp||2 | max sv(P_perp B) |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in payload["rank_rows"]:
        lines.append(
            "| {M} | {r} | {rank_H}/{data_dim} | {nuisance_codimension} | "
            "{n_vis} | {hidden_rank} | {projector_norm:.3e} | {max_sv:.3e} |".format(
                max_sv=max(row["hid_sv"]), **row
            )
        )
    lines.extend(
        [
            "",
            "The transition follows a dimension obstruction, not a special "
            "Fourier/Hankel theorem. With p=3 real pose variables and K=3 "
            "real contrast variables, the generic nuisance column count is "
            "2r+K in a 2M-dimensional real data space. At r=M-2 the "
            "orthogonal complement has dimension at most one, so at least "
            "two pose directions are hidden. At r>=M-1 a full-rank nuisance "
            "span fills the data space, so all three pose directions are "
            "hidden; the old relative-only rule mistook roundoff for three "
            "visible directions.",
            "",
            "## Corrected nonlinear control",
            "",
            "| setting | method | regime | success | median pose error | median map error | median residual |",
            "|---|---|---|---:|---:|---:|---:|",
        ]
    )
    for row in payload["case_summary"]:
        lines.append(
            "| {setting} | {method} | {regime} | {successes}/{n} | "
            "{median_pose_error:.4e} | {median_map_error:.4e} | "
            "{median_residual:.4e} |".format(**row)
        )
    lines.extend(
        [
            "",
            "A reduced solver with zero visible pose directions freezes pose "
            "at the nominal value. Therefore equality with the direct solver "
            "at r=M-1 in the locked report was a numerical-rank artifact, not "
            "evidence of lossless reduction.",
            "",
            "## Scientific consequence",
            "",
            "The reliable result is a pose-hiding diagnostic and a dimension "
            "bound. It is not a hidden-direction recovery algorithm. Full "
            "pose reparameterization (visible rank three) is algebraically "
            "equivalent to direct joint optimization, while truncating to the "
            "visible row space intentionally removes hidden directions.",
        ]
    )
    (HERE / "rank_saturation_summary.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def make_plot(payload: dict) -> None:
    """Create a compact corrected figure for the manuscript."""
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 4.9))

    ax = axes[0]
    colors = {8: "#1f77b4", 12: "#2ca02c", 16: "#d62728"}
    for M in (8, 12, 16):
        rows = [x for x in payload["rank_rows"] if x["M"] == M]
        rows.sort(key=lambda x: x["r"])
        rs = [x["r"] for x in rows]
        q = [x["n_vis"] for x in rows]
        h = [x["hidden_rank"] for x in rows]
        ax.plot(rs, q, "o-", color=colors[M], label=f"visible, M={M}")
        ax.plot(
            rs,
            h,
            "s--",
            color=colors[M],
            alpha=0.72,
            label=f"hidden, M={M}",
        )
    ax.set_xlabel("retained complex-current rank r")
    ax.set_ylabel("real pose dimension")
    ax.set_yticks((0, 1, 2, 3))
    ax.set_title("Corrected pose rank after nuisance projection")
    ax.grid(alpha=0.25)
    ax.legend(ncol=2, fontsize=8)

    ax = axes[1]
    summary = payload["case_summary"]
    labels = []
    pose = []
    mape = []
    residual = []
    for M in (12, 16):
        for regime, short in (
            ("direct_control", "direct"),
            ("codim_one", "codim=1"),
            ("nuisance_saturated", "codim=0"),
        ):
            row = next(
                x
                for x in summary
                if x["setting"].startswith(f"M{M}_")
                and x["regime"] == regime
            )
            labels.append(f"M{M}\n{short}")
            pose.append(row["median_pose_error"])
            mape.append(row["median_map_error"])
            residual.append(row["median_residual"])
    x = np.arange(len(labels))
    width = 0.24
    ax.bar(x - width, pose, width, label="pose error", color="#1f77b4")
    ax.bar(x, mape, width, label="map error", color="#ff7f0e")
    ax.bar(x + width, residual, width, label="data residual", color="#2ca02c")
    ax.axhline(10.0 ** (-30.0 / 20.0), color="k", ls=":", lw=1.1, label="30 dB ratio")
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("median over three noise seeds")
    ax.set_title("Direct control vs rank-safe reduced coordinates")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(fontsize=8, ncol=2)

    fig.suptitle(
        "Nuisance-space saturation hides pose; relative-only rank creates a false recovery",
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(HERE / "rank_saturation_corrected.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    rank_rows: list[dict] = []
    cases: list[dict] = []
    for M in (8, 12, 16):
        P = params_for(M)
        xs, h = pixel_grid(P.N)
        G_D = green_domain_matrix(xs, h, P.k)
        Phi = basis_matrix(xs, P)
        subspaces: dict[int, dict] = {}
        for r in range(2, M + 1):
            sp = corrected_subspace(
                ALPHA_INIT, P_INIT, r, xs, h, G_D, Phi, P
            )
            subspaces[r] = sp
            rank_rows.append(
                {
                    "M": M,
                    "r": r,
                    **{
                        k: v
                        for k, v in sp.items()
                        if k not in ("rank", "V_vis")
                    },
                }
            )

        # Nonlinear controls are most informative for M=12 and M=16, for
        # which the locked E10b report claimed a trailing lossless state.
        if M not in (12, 16):
            continue
        setting = {
            "label": f"M{M}_k12_full_corrected",
            "M": M,
            "k": 12.0,
            "aperture": "full",
        }
        d_true = forward(ALPHA_TRUE, P_TRUE, xs, h, P.k, G_D, Phi, P)
        transition_sp = subspaces[M - 2]
        saturated_sp = subspaces[M - 1]
        assert transition_sp["n_vis"] == 1
        assert saturated_sp["n_vis"] == 0
        for seed in SEEDS:
            d_obs, _ = scaled_observation(d_true, make_noise(seed, M))
            cases.append(
                run_case(
                    "direct",
                    None,
                    seed,
                    setting,
                    P,
                    d_obs,
                    xs,
                    h,
                    P.k,
                    G_D,
                    Phi,
                    None,
                    transition_sp,
                    "direct_control",
                )
            )
            cases.append(
                run_case(
                    "reduced",
                    M - 2,
                    seed,
                    setting,
                    P,
                    d_obs,
                    xs,
                    h,
                    P.k,
                    G_D,
                    Phi,
                    transition_sp,
                    transition_sp,
                    "codim_one",
                )
            )
            cases.append(
                run_case(
                    "reduced",
                    M - 1,
                    seed,
                    setting,
                    P,
                    d_obs,
                    xs,
                    h,
                    P.k,
                    G_D,
                    Phi,
                    saturated_sp,
                    transition_sp,
                    "nuisance_saturated",
                )
            )

    payload = {
        "schema": "trispace_rank_saturation_audit_v1",
        "historical_outputs_modified": False,
        "rank_rule": {
            "H_rank_tol_after_column_normalization": H_RANK_TOL,
            "B_relative_tol": B_REL_TOL,
            "B_backward_factor": B_BACKWARD_FACTOR,
            "saturation_rule": "rank(H)=data_dim implies n_vis=0 exactly",
        },
        "dimension_bound": (
            "hidden_dim >= max(0, p + dim(N_U) - data_dim), with "
            "p=3 and data_dim=2M"
        ),
        "rank_rows": rank_rows,
        "cases": cases,
        "case_summary": summarize_cases(cases),
        "source_sha256": {
            "this_script": sha256(Path(__file__)),
            "run_e5_final.py": sha256(EXP / "run_e5_final.py"),
            "run_e10b_m12_m16_highrank.py": sha256(
                EXP / "run_e10b_m12_m16_highrank.py"
            ),
        },
    }
    out = HERE / "rank_saturation_results.json"
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    write_markdown(payload)
    make_plot(payload)
    print(json.dumps(payload["case_summary"], indent=2))


if __name__ == "__main__":
    main()
