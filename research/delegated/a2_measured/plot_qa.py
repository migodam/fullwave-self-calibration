"""QA figures for the 2001 Fresnel dielTM_dec8f ingestion pass."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from ingestion import read_rows

ROOT = Path(__file__).resolve().parent
PLOTS = ROOT / "plots"


def main() -> int:
    PLOTS.mkdir(exist_ok=True)
    arr, _ = read_rows()
    view = arr[:, 0].astype(int)
    rec = arr[:, 1].astype(int)
    freq = arr[:, 2].astype(int)
    total = arr[:, 3] + 1j * arr[:, 4]
    incident = arr[:, 5] + 1j * arr[:, 6]

    # Fig 1: measured amplitude / phase summary with resolved column labels
    fig, axes = plt.subplots(3, 2, figsize=(11, 9), sharex=True)
    for row, fq in enumerate([1, 4, 8]):
        sel = (view == 1) & (freq == fq)
        rs = rec[sel]
        th = np.deg2rad((rs - 1) * 5.0)
        t, i = total[sel], incident[sel]
        s = t - i
        axes[row, 0].plot(np.rad2deg(th), np.abs(t), label="total field E", lw=1.5)
        axes[row, 0].plot(
            np.rad2deg(th), np.abs(i), label="incident field E (target removed)", lw=1.5
        )
        axes[row, 0].plot(np.rad2deg(th), np.abs(s), label="scattered = total - incident", lw=1.5)
        axes[row, 0].set_ylabel(f"{fq} GHz\n|E| (arbitrary units)")
        axes[row, 0].grid(alpha=0.3)
        axes[row, 1].plot(np.rad2deg(th), np.unwrap(np.angle(t)), lw=1.0)
        axes[row, 1].plot(np.rad2deg(th), np.unwrap(np.angle(i)), lw=1.0)
        axes[row, 1].set_ylabel("arg E (rad)")
        axes[row, 1].grid(alpha=0.3)
    axes[0, 0].set_title("amplitude |E| vs receiver angle (view 1)")
    axes[0, 1].set_title("unwrapped phase arg E vs receiver angle (view 1)")
    axes[0, 0].legend(fontsize=8)
    for ax in axes[-1]:
        ax.set_xlabel("receiver angle (deg); receiver index = 1 + angle/5")
    fig.suptitle(
        "dielTM_dec8f.exp: TM total vs incident electric field (exp(+i w t), "
        "E parallel to cylinder axis)",
        fontsize=10,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(PLOTS / "qa_amplitude_phase.png", dpi=150)
    plt.close(fig)

    # Fig 2: incident-field rotation consistency across views (emitter pattern)
    fig, axes = plt.subplots(2, 2, figsize=(11, 7))
    for ax, fq in zip(axes.flat, [1, 3, 5, 8]):
        s1 = (view == 1) & (freq == fq)
        e1 = incident[s1]
        r1 = rec[s1]
        amp_v1 = np.abs(e1)
        # view 2 receiver window is +2 indices; compare overlapping receivers
        s2 = (view == 2) & (freq == fq)
        e2, r2 = incident[s2], rec[s2]
        m = {rr: v for rr, v in zip(r1, amp_v1)}
        shared = [rr for rr in r2 if (rr - 2) in m]
        ratio = [np.abs(e2[list(r2).index(rr)]) / m[rr - 2] for rr in shared]
        ax.plot([int(rr) for rr in shared], ratio, "o", ms=3)
        ax.axhline(1.0, color="k", lw=0.8, ls="--")
        ax.set_title(f"{fq} GHz")
        ax.set_ylabel("|E_inc(view2, rec-2)| / |E_inc(view1, rec)|")
        ax.set_xlabel("receiver index (view 2)")
        ax.set_ylim(0.95, 1.05)
        ax.grid(alpha=0.3)
    fig.suptitle(
        "Incident-field rotational consistency: view v pattern == view 1 pattern "
        "rotated by the source step",
        fontsize=10,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(PLOTS / "qa_incident_rotation_consistency.png", dpi=150)
    plt.close(fig)

    print(f"wrote plots to {PLOTS}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
