"""Model-fit figures for the bounded 2001 Fresnel pilot."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import model as M

ROOT = Path(__file__).resolve().parent
PLOTS = ROOT / "plots"


def main() -> int:
    PLOTS.mkdir(exist_ok=True)
    res = json.load(open(ROOT / "results_model.json"))
    view, rec, freq, total, incident = M.load_data()
    scattered = total - incident

    # Fig 1: measured vs modeled scattered field, one train view, 2 and 8 GHz
    joint = res["fits"]["joint_3param"]["params"]
    nominal = res["fits"]["nominal"]["params"]
    v = 2
    fig, axes = plt.subplots(2, 2, figsize=(11, 7))
    for col, fq in enumerate([2, 8]):
        sel = (view == v) & (freq == fq)
        rs = rec[sel]
        s_meas = scattered[sel]
        # gains profiled on train views for the joint params (identical nuisance)
        key, _ = M.build_index(view, rec, freq, scattered)
        train = np.arange(2, 37, 2)
        _, meta = M.gains_and_residual(
            (joint["cx_m"], joint["cy_m"], joint["eps_r"]),
            joint["eps_r"] + 0j,
            train,
            np.arange(1, 9),
            scattered,
            key,
        )
        g = meta["gain"][int(fq)]
        Mj = M.forward_scattered(fq, v, (joint["cx_m"], joint["cy_m"]), joint["eps_r"] + 0j)
        selr = M.rec_idx_for_view(v)
        s_mod = g * Mj[selr]
        axes[0, col].plot(rs, np.abs(s_meas), "o-", ms=3, lw=1, label="measured")
        axes[0, col].plot(rs, np.abs(s_mod), "s-", ms=3, lw=1, label="model (joint fit)")
        axes[0, col].set_title(f"{fq} GHz, view {v}")
        axes[0, col].set_ylabel("|scattered field| (arbitrary units)")
        axes[0, col].legend(fontsize=8)
        axes[0, col].grid(alpha=0.3)
        axes[1, col].plot(rs, np.unwrap(np.angle(s_meas)), "o-", ms=3, lw=1)
        axes[1, col].plot(rs, np.unwrap(np.angle(s_mod)), "s-", ms=3, lw=1)
        axes[1, col].set_ylabel("arg scattered (rad)")
        axes[1, col].set_xlabel("receiver index (5 deg steps)")
        axes[1, col].grid(alpha=0.3)
    fig.suptitle(
        "Measured vs modeled scattered field (per-frequency gain profiled on train views)",
        fontsize=10,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(PLOTS / "fit_scattered_comparison.png", dpi=150)
    plt.close(fig)

    # Fig 2: residual by frequency (nominal vs joint) and band-stability of eps
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    fq_list = [str(i) for i in range(1, 9)]
    nom_tr = [res["fits"]["nominal"]["train_per_freq"][k] for k in fq_list]
    jnt_tr = [res["fits"]["joint_3param"]["train_per_freq"][k] for k in fq_list]
    jnt_te = [res["fits"]["joint_3param"]["test_per_freq"][k] for k in fq_list]
    x = np.arange(8)
    w = 0.27
    axes[0].bar(x - w, nom_tr, w, label="nominal geometry (train)")
    axes[0].bar(x, jnt_tr, w, label="joint fit (train)")
    axes[0].bar(x + w, jnt_te, w, label="joint fit (held-out)")
    axes[0].set_xticks(x)
    axes[0].set_xticklabels([f"{i} GHz" for i in range(1, 9)])
    axes[0].set_ylabel("relative residual ||g M - s||^2 / ||s||^2")
    axes[0].legend(fontsize=8)
    axes[0].grid(alpha=0.3, axis="y")

    diag = res["low_frequency_diagnostic"]
    bands = ["up_to_2GHz", "up_to_3GHz", "up_to_4GHz", "up_to_8GHz"]
    eps_by_band = [
        diag["up_to_2GHz"]["eps_r"],
        diag["up_to_3GHz"]["eps_r"],
        diag["up_to_4GHz"]["eps_r"],
        res["fits"]["joint_3param"]["params"]["eps_r"],
    ]
    rad_by_band = [
        diag["up_to_2GHz"]["center_radius_mm"],
        diag["up_to_3GHz"]["center_radius_mm"],
        diag["up_to_4GHz"]["center_radius_mm"],
        res["fits"]["joint_3param"]["center_radius_mm"],
    ]
    axes[1].plot(bands, eps_by_band, "o-", label="fitted eps_r")
    axes[1].axhspan(2.7, 3.3, color="green", alpha=0.15, label="published eps_r = 3 +/- 0.3")
    axes[1].set_ylabel("fitted eps_r")
    ax2 = axes[1].twinx()
    ax2.plot(bands, rad_by_band, "s--", color="red", label="center radius (mm)")
    ax2.axhline(30, color="red", alpha=0.4, ls=":")
    ax2.set_ylabel("center radius (mm)")
    axes[1].set_title("band stability: center stable, eps_r not identifiable")
    axes[1].legend(loc="upper left", fontsize=8)
    ax2.legend(loc="upper right", fontsize=8)
    fig.suptitle("dielTM_dec8f: fit quality and parameter stability", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(PLOTS / "fit_residual_and_stability.png", dpi=150)
    plt.close(fig)

    # Fig 3: robustness pilot - residual landscape over (cx, cy) at fixed eps
    inj = res["robustness_pilot"]["injected"]
    recv = res["robustness_pilot"]["recovered"]
    # regenerate the synthetic data exactly as the pilot does
    rng = np.random.default_rng(7)
    c_inj = np.array([inj["cx_m"], inj["cy_m"]])
    eps_inj = inj["eps_r"] + 0j
    train_views = np.arange(2, 37, 2)
    freqs = np.arange(1, 9)
    key_syn = {}
    for fq in freqs:
        batch = M.forward_view_batch(int(fq), range(1, 37), c_inj, eps_inj)
        syn = np.concatenate([batch[v][M.rec_idx_for_view(v)] for v in range(1, 37)])
        sigma = 0.05 * np.sqrt(np.mean(np.abs(syn) ** 2))
        syn = syn + sigma * (
            rng.standard_normal(syn.size) + 1j * rng.standard_normal(syn.size)
        )
        for j, vv in enumerate(range(1, 37)):
            if vv in train_views:
                key_syn[(int(vv), int(fq))] = syn[j * 49 : (j + 1) * 49]

    def resid(cx, cy):
        tot_num = tot_den = 0.0
        for fq in freqs:
            batch = M.forward_view_batch(int(fq), train_views, (cx, cy), eps_inj)
            ml, sl = [], []
            for vv in train_views:
                ml.append(batch[int(vv)][M.rec_idx_for_view(vv)])
                sl.append(key_syn[(int(vv), int(fq))])
            Mv, sv = np.concatenate(ml), np.concatenate(sl)
            g = np.sum(np.conj(Mv) * sv) / np.sum(np.abs(Mv) ** 2)
            tot_num += np.sum(np.abs(g * Mv - sv) ** 2)
            tot_den += np.sum(np.abs(sv) ** 2)
        return tot_num / tot_den

    xs = np.linspace(0.012, 0.036, 25)
    ys = np.linspace(-0.02, 0.004, 25)
    XX, YY = np.meshgrid(xs, ys)
    ZZ = np.empty_like(XX)
    for i in range(XX.shape[0]):
        for j in range(XX.shape[1]):
            ZZ[i, j] = resid(XX[i, j], YY[i, j])
    fig, ax = plt.subplots(figsize=(7, 5.5))
    cs = ax.contourf(XX * 1e3, YY * 1e3, ZZ, levels=30, cmap="viridis")
    fig.colorbar(cs, label="relative residual")
    ax.plot(inj["cx_m"] * 1e3, inj["cy_m"] * 1e3, "w*", ms=16, label="injected (truth)")
    ax.plot(recv["cx_m"] * 1e3, recv["cy_m"] * 1e3, "r.", ms=14, label="recovered")
    ax.set_xlabel("cx (mm)")
    ax.set_ylabel("cy (mm)")
    ax.set_title("Robustness pilot: residual basin at fixed eps_r = 2.7")
    ax.legend()
    fig.tight_layout()
    fig.savefig(PLOTS / "robustness_basin.png", dpi=150)
    plt.close(fig)

    print(f"wrote model plots to {PLOTS}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
