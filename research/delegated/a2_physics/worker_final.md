Done. All bounded checks pass (8/8, ~5 s), and the N16/N32 timing runs complete cleanly.

Artifacts in [research/delegated/a2_physics](/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/a2_physics):
- [physics.py](/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/a2_physics/physics.py) — `Config`/`Model`, full-wave forward, analytic A/B, realify, material-basis override, and an adjoint `gradient` with honest per-RHS/factorisation accounting.
- [test_physics.py](/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/a2_physics/test_physics.py) — executable check suite; rerun with the named venv python.
- [checks.json](/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/a2_physics/checks.json), [summary.md](/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/a2_physics/summary.md), [environment.json](/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/a2_physics/environment.json).

Facts
- State residuals are machine precision for every frequency/pose/illumination (max rel ~1.8e-15 at N=32); raw total equals retained direct + scattered components.
- Analytic material A and pose B match centred finite differences at nonzero pose for both apertures across all four frequencies: worst normalised errors 2.5e-10 (A) and 2.7e-10 (B).
- Corrected self-cell disk integral reproduces radial quadrature to ~1e-15; D and S entries verified to carry the specified k² and cell-integration factors.
- Adjoint gradient matches finite differences of the realified objective (~2e-13 abs) and charges 12 forward + 12 adjoint RHS over two frequencies with no extra factorisations.
- One LU factorisation per frequency is shared by all 6 pose/illumination states; all RHS columns are charged individually (312 total, 4 factorisations per full forward).
- N16 forward+jacobian: 0.167 s; N32: 2.34 s (same 288 rows, 312 RHS).
- Pose 1 is anchored (B block exactly zero); direct-incident derivatives are numerically zero under rigid co-motion as specified.

Uncertainties
- “Width 0.16 m” was implemented as the standard deviation in exp(−|x−c|²/(2·0.16²)), unnormalised, consistent with the surrounding loop-code convention; and the shared x is interpreted as adding the same (dx, dy, dθ) transform to both non-anchored nominal poses. Both are documented in `physics.py` and `summary.md`.
- These are internal-consistency checks of the implemented discrete core, not continuum convergence or production acceptance evidence. No seeds 1001+ or optimizer runs were executed.

No intervention needed from Codex for this bounded package.