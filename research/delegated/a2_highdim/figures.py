"""PNG figures for the a2_highdim exploratory comparison."""

from __future__ import annotations

from typing import Any

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from common import DIR, make_model, make_scene  # noqa: E402

FIGS = DIR / "figures"

METHOD_LABELS = {
    "coherent_fixedpose": "coherent, fixed pose",
    "coherent_joint": "coherent joint map+pose",
    "intensity_joint": "matched intensity joint",
    "oracle": "oracle (known true pose)",
}


def _records_of(records: list[dict[str, Any]], phase: str) -> list[dict[str, Any]]:
    return [r for r in records if r.get("phase") == phase]


def reconstruction_figure(records: list[dict[str, Any]]) -> None:
    test = _records_of(records, "test")
    picks: list[dict[str, Any]] = []
    for method in ("coherent_fixedpose", "coherent_joint", "intensity_joint", "oracle"):
        match = [
            r for r in test if r["seed"] == 801
            and r["radius_index"] == 0 and r["method"] == method
        ]
        if match:
            picks.append(match[0])
    if not picks:
        print("[figures] no reconstruction records for seed 801/.125lambda")
        return
    model_data = make_model(32, "gaussian49")
    seed = picks[0]["seed"]
    scene = make_scene(seed)
    alpha_true = np.asarray(scene["alpha_true"], dtype=float)
    chi_true = model_data.basis @ alpha_true
    vmin, vmax = float(chi_true.min()), float(chi_true.max())
    maps = [("truth", chi_true, None, False)]
    for r in picks:
        chi = model_data.basis @ np.asarray(r["alpha_est"], dtype=float)
        rmse = float(np.sqrt(np.mean((chi - chi_true) ** 2)))
        maps.append((r["method"], chi, rmse, bool(r["is_oracle"])))

    n = len(maps)
    fig, axes = plt.subplots(
        1, n, figsize=(3.1 * n, 3.0), squeeze=False, layout="constrained"
    )
    fig.suptitle(
        f"Exploratory 2D reconstruction, seed {seed}, .125 lambda pose error "
        "(spatial epsilon_r - 1 on N32 grid)"
    )
    for ax, (method, chi, rmse, oracle) in zip(axes[0], maps):
        im = ax.imshow(
            chi.reshape(32, 32).T,
            origin="lower",
            extent=(-0.5, 0.5, -0.5, 0.5),
            vmin=vmin,
            vmax=vmax,
        )
        label = "truth" if method == "truth" else METHOD_LABELS.get(method, method)
        if oracle:
            label += "\n(explicit oracle upper reference)"
        title = label
        if rmse is not None:
            title += f"\nspatial RMSE {rmse:.4f}"
        ax.set_title(title, fontsize=9)
        ax.set_xlabel("x (m)")
        if ax is axes[0][0]:
            ax.set_ylabel("y (m)")
        else:
            ax.set_yticks([])
    fig.colorbar(
        im, ax=axes[0].tolist(), shrink=0.85, pad=0.01, label=r"$\epsilon_r-1$"
    )
    path = FIGS / "reconstruction_map.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)
    print(f"[figures] wrote {path}")


def _test_metric_arrays(test: list[dict[str, Any]], key: str, method: str):
    vals = [
        r["metrics"].get(key)
        for r in test
        if r["method"] == method and r["metrics"].get(key) is not None
    ]
    return np.asarray(vals, dtype=float)


def error_distributions(records: list[dict[str, Any]]) -> None:
    test = _records_of(records, "test")
    methods = ["coherent_fixedpose", "coherent_joint", "intensity_joint", "oracle"]
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.4))
    titles = [
        ("spatial_map_rmse", "Spatial map RMSE (epsilon_r-1)"),
        ("coefficient_rmse", "49-coefficient RMSE"),
        ("pose_error_lever_m", "Pose lever-arm error (m)"),
    ]
    for ax, (key, title) in zip(axes, titles):
        data = [_test_metric_arrays(test, key, m) for m in methods]
        colors = ["#7f7f7f", "#1f77b4", "#d62728", "#2ca02c"]
        bp = ax.boxplot(
            data,
            tick_labels=[METHOD_LABELS[m] for m in methods],
            patch_artist=True,
            showfliers=False,
            widths=0.45,
        )
        for patch, c in zip(bp["boxes"], colors):
            patch.set_facecolor(c)
            patch.set_alpha(0.5)
        for idx, (d, c) in enumerate(zip(data, colors), start=1):
            jitter = np.random.default_rng(12345 + idx).uniform(-0.13, 0.13, d.size)
            ax.scatter(idx + jitter, d, color=c, s=12, alpha=0.8, zorder=3)
        ax.set_title(title, fontsize=10)
        ax.tick_params(axis="x", labelsize=7.5)
        ax.set_ylabel("unconditional error (all 12 final runs per method)")
    fig.tight_layout()
    path = FIGS / "error_distributions.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)
    print(f"[figures] wrote {path}")


def phase_residual_figure(records: list[dict[str, Any]]) -> None:
    test = _records_of(records, "test")
    methods = ["coherent_fixedpose", "coherent_joint", "intensity_joint", "oracle"]
    fig, ax = plt.subplots(figsize=(7.6, 4.4))
    bins = np.linspace(-np.pi, np.pi, 81)
    colors = ["#7f7f7f", "#1f77b4", "#d62728", "#2ca02c"]
    for method, c in zip(methods, colors):
        samples = [
            v
            for r in test
            if r["method"] == method
            for v in r.get("phase_residual_masked_rad", [])
        ]
        ax.hist(
            samples,
            bins=bins,
            density=True,
            alpha=0.38,
            color=c,
            label=f"{METHOD_LABELS[method]} (n={len(samples):,})",
        )
    ax.set_xlim(-np.pi, np.pi)
    ax.set_xticks([-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi])
    ax.set_xticklabels([r"$-\pi$", r"$-\pi/2$", "0", r"$\pi/2$", r"$\pi$"])
    ax.axvline(0, color="k", lw=0.7, ls=":")
    ax.set_xlabel("wrapped forward-phase error (rad), vs clean truth total")
    ax.set_ylabel("density")
    ax.set_title(
        "Masked channels only: |truth total| > 3 sigma "
        "(evaluation-only mask)"
    )
    ax.legend(fontsize=8)
    fig.tight_layout()
    path = FIGS / "phase_residuals.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)
    print(f"[figures] wrote {path}")


def mismatch_figure(records: list[dict[str, Any]]) -> None:
    mm = _records_of(records, "mismatch")
    if not mm:
        print("[figures] no mismatch records")
        return
    kinds = ["clock_phase", "receiver_coupling", "outside_basis"]
    methods = ["coherent_fixedpose", "coherent_joint"]
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2))
    xpos = np.arange(len(kinds))
    width = 0.34
    for ai, (ax, key) in enumerate(
        zip(axes, ("validation_chi2_p", "fit_chi2_p"))
    ):
        for mi, method in enumerate(methods):
            vals = [
                next(
                    r["metrics"][key]
                    for r in mm
                    if r["mismatch_kind"] == k and r["method"] == method
                )
                for k in kinds
            ]
            ax.bar(
                xpos + (mi - 0.5) * width,
                vals,
                width=width,
                label=METHOD_LABELS[method],
                color=["#1f77b4", "#d62728"][mi],
                alpha=0.75,
            )
        ax.axhline(0.01, color="k", lw=0.8, ls="--", label="1% reference line")
        ax.set_xticks(xpos)
        ax.set_xticklabels(["clock phase", "receiver\ncoupling", "outside 49\nbasis"],
                           fontsize=8)
        ax.set_yscale("log")
        ax.set_title("whitened residual chi-square p" if "val" in key
                     else "in-sample whitened chi-square p", fontsize=9)
        ax.legend(fontsize=7)
    fig.suptitle(
        "Controlled mismatch controls (seed 881): residual diagnostics, "
        "no calibration guarantee", y=1.02)
    fig.tight_layout()
    path = FIGS / "mismatch_residuals.png"
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"[figures] wrote {path}")


def make_all_figures(records: list[dict[str, Any]]) -> None:
    FIGS.mkdir(exist_ok=True)
    reconstruction_figure(records)
    error_distributions(records)
    phase_residual_figure(records)
    mismatch_figure(records)
