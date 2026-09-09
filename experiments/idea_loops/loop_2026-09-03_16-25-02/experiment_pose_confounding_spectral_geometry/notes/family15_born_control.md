# Family 15: Born analytic / literature control

Date: 2026-09-04 (SGT); UTC stamp in `results/family15_born_control.json`.
Experiment: `experiment_pose_confounding_spectral_geometry`.

This family is a formula-level/internal analytic control that extends, and does not duplicate, the earlier Born checks: Family 3/3b verified `A(chi=0) == A_born` and diagnosed the FD slope gates, and Family 7 reported Born-vs-full-wave A discrepancy plus full-wave no-prior retention over a coarse contrast grid.  Family 15 compares the Born operator, `K_IS`, and `K_eff` against the full-wave local quantities in the smooth p=24 coefficient space, under one shared full-wave pose confounder block, and records observed convergence orders and the overlap of the most-confounded subspaces.

## Exact command, runtime, platform

```bash
.venv/bin/python src/family15_born_control.py
```

Wall runtime: 0.770 s (also in the JSON as `wall_runtime_seconds`).
Platform: Apple Silicon CPU (`arm64`), macOS-26.6.2-arm64-arm-64bit.
Python 3.12.13, numpy 2.5.2, scipy 1.18.1, matplotlib 3.11.1.  Deterministic dense linear algebra; no RNG used.

## Scene and configuration

* N = 16 (`h_cell = 1/16`), k_b = 2*pi*f with f = 1.0, T = 6, n_rx = 4.
* family1 90-degree arc, radius 1.6, phi in [-45, 45] deg, theta = atan2(-p_y, -p_x); rx offsets [[-0.06, 0.0], [0.06, 0.0], [0.0, -0.06], [0.0, 0.06]], tx offset [0.0, 0.0].
* `chi_s = s * chi0`, `chi0 = family1.make_chi0(points, cfg)` (two_blob); scales [0.001, 0.01, 0.03, 0.1, 0.3, 1.0].  s=0.1 is the fixed named low-contrast row (`0.1 * two_blob`).
* Smooth p=24 unit-2-norm Gaussian RBF basis `family2.build_smooth_basis` with the family-2/4 config (4 x 6 centres, sigma_b = 0.16); `A_s = A_pix_R @ S` in every row.
* Whitening: identity noise, `whiten_realify(A_c, B_c, None)` (sqrt(2) Re/Im stacking); real data rows 48.
* Finite prior alpha = 1.0.
* Rank tolerance: tol(M) = max(M.shape) * eps_machine * sigma_1(M).

## Born composition identity

`hh.born_forward(chi_s, poses, rx_offsets, tx_offset, N, k_b)` returns `(F_born, A_born)` with complex rows grouped by pose; each block is `A_born,t = G_S_t diag(E_inc_t)` (equivalently the full-wave Jacobian at chi = 0).  The identity is verified twice: once against an explicit `G_S diag(E_inc)` stack built from `hh.build_operators`, and once in the whitened/projected smooth space `A_s_born` used below.

| s | rel Fro explicit-vs-born (complex) | rel Fro explicit-vs-born (real, smooth) |
| --- | --- | --- |
| 0.001 | 0 | 0 |
| 0.01 | 0 | 0 |
| 0.03 | 0 | 0 |
| 0.1 | 0 | 0 |
| 0.3 | 0 | 0 |
| 1 | 0 | 0 |

The residuals are exactly 0.0 (bitwise-identical blockwise assembly): `hh.born_forward` and the explicit `G_S diag(E_inc)` stack evaluate the same deterministic green-matrix formula with the same inputs.  This confirms, at the code-path level, that the Born forward map used in the FIM comparison is the self-derived `A_born = G_S diag(E_inc)` map; the independent numerical identity `A(chi=0) == A_born` (7e-17 relative Frobenius) is Family 3b check D.

## No-prior retention sweep (P_perp = I - Z Z^T)

`Z = family2.range_basis(B_R)` with the full-wave `B_R` at chi_s; the same `P_perp` is applied to full-wave and Born rows (see JSON `config.retention.same_B_R_note`).

| s | rdof full | rdof Born | log-vol full | log-vol Born | rho_min full | rho_min Born | theta_min full (deg) | theta_min Born (deg) | rel A | rel K_IS | rel K_eff | d_rdof (full-born) | d_rho_min (full-born) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0.001 | 9.15161 | 9.15117 | -161.631 | -161.632 | 8.358360e-11 | 8.362399e-11 | 5.238213e-04 | 5.239479e-04 | 0.000499888 | 0.000701769 | 0.000701769 | 0.000438214 | -4.03947e-14 |
| 0.01 | 9.15633 | 9.1517 | -161.635 | -161.645 | 8.339258e-11 | 8.379677e-11 | 5.232226e-04 | 5.244888e-04 | 0.00500362 | 0.0070262 | 0.0070262 | 0.00463626 | -4.04192e-13 |
| 0.03 | 9.16842 | 9.15287 | -161.651 | -161.673 | 8.297610e-11 | 8.418974e-11 | 5.219141e-04 | 5.257173e-04 | 0.0150417 | 0.0211333 | 0.0211333 | 0.0155581 | -1.21364e-12 |
| 0.1 | 9.22456 | 9.15703 | -161.779 | -161.771 | 8.159503e-11 | 8.563664e-11 | 5.175530e-04 | 5.302154e-04 | 0.0504733 | 0.0710059 | 0.0710046 | 0.0675323 | -4.04161e-12 |
| 0.3 | 9.35787 | 9.16949 | -162.844 | -162.048 | 7.838932e-11 | 9.027141e-11 | 5.072838e-04 | 5.443749e-04 | 0.153554 | 0.215586 | 0.21555 | 0.18838 | -1.18821e-11 |
| 1 | 9.40716 | 9.21936 | -167.746 | -163.055 | 5.510027e-11 | 1.047412e-10 | 4.253044e-04 | 5.863834e-04 | 0.504938 | 0.652688 | 0.651679 | 0.187799 | -4.9641e-11 |

## Observed convergence order (log10 error vs log10 s, full grid)

| quantity | slope | intercept | r^2 | n |
| --- | --- | --- | --- | --- |
| A operator | 1.00261 | -0.294199 | 0.999994 | 6 |
| K_IS | 0.994365 | -0.16368 | 0.999841 | 6 |
| K_eff | 0.994199 | -0.164027 | 0.999837 | 6 |

Reported as observed, not forced.  A slope near 1 is expected if the relative Born error is dominated by the O(s) first correction; deviation at the largest scales (and any imperfect log-log linearity) is left visible in the JSON and figure.

## Finite-prior table (alpha = 1.0)

`W = I - B_R (B_R^T B_R + alpha I)^{-1} B_R^T` and `K_eff = A_s^T W A_s` for both full-wave and Born A_s.

| s | rdof eff full | rdof eff Born | rho_min eff full | rho_min eff Born | theta_min full (deg) | theta_min Born (deg) | rel K_eff |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0.001 | 24 | 24 | 1.000000e+00 | 1.000000e+00 | 8.999439e+01 | 8.999439e+01 | 0.000701769 |
| 0.01 | 24 | 24 | 9.999990e-01 | 9.999990e-01 | 8.994392e+01 | 8.994392e+01 | 0.0070262 |
| 0.03 | 24 | 24 | 9.999914e-01 | 9.999914e-01 | 8.983191e+01 | 8.983190e+01 | 0.0211333 |
| 0.1 | 23.9996 | 23.9996 | 9.999050e-01 | 9.999050e-01 | 8.944144e+01 | 8.944142e+01 | 0.0710046 |
| 0.3 | 23.9966 | 23.9966 | 9.991564e-01 | 9.991565e-01 | 8.833565e+01 | 8.833574e+01 | 0.21555 |
| 1 | 23.9602 | 23.9602 | 9.895111e-01 | 9.895107e-01 | 8.412172e+01 | 8.412161e+01 | 0.651679 |

## Most-confounded subspace comparison (full-wave vs Born)

The three most pose-confounded directions (`u0/u1/u2`, ascending generalized retention) were computed with `family4.generalized_eigen_directions(A_s, P_perp)`; the Born directions use the same full-wave `P_perp`.

| s | diag sq overlap u0 | u1 | u2 | best-pairing sum | principal angle 1 (deg) | angle 2 (deg) | angle 3 (deg) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0.001 | 1 | 0.999965 | 1 | 2.99997 | 9.723736e-04 | 1.812114e-02 | 4.486077e-02 |
| 0.01 | 1 | 0.996494 | 0.999993 | 2.99649 | 9.565515e-03 | 1.779941e-01 | 4.454117e-01 |
| 0.03 | 0.999998 | 0.96811 | 0.999948 | 2.96806 | 2.758982e-02 | 5.124475e-01 | 1.319206e+00 |
| 0.1 | 0.99998 | 0.704295 | 0.999642 | 2.70392 | 7.668802e-02 | 1.463618e+00 | 4.325013e+00 |
| 0.3 | 0.999731 | 0.215607 | 0.99087 | 2.20621 | 2.221044e-01 | 3.806400e+00 | 1.665510e+01 |
| 1 | 0.953359 | 0.0271854 | 0.958752 | 2.12697 | 9.713510e-01 | 1.532509e+01 | 6.485976e+01 |

The full 3x3 normalized squared-overlap matrix, top-3 direction coefficient vectors, and the generalized spectra are stored in the JSON for audit.

## Literature field and scope

**Citation (verbatim):** M. L. Diong, A. Roueff, P. Lasaygues, A. Litman, "Impact of the Born approximation on the estimation error in 2D inverse scattering", Inverse Problems, 2016.

The paper quantifies the Born approximation's effect on estimation error by comparing the linear Born MLE variance with the full nonlinear Cramer-Rao bound (CRB).  Family 15 is a formula-level/internal analytic control with the same conceptual structure: the Born FIM `K_IS_born = A_born^T A_born` is the inverse covariance of the linear Born Gaussian estimator, and the full-wave `K_IS`/`K_eff` is the local FIM of the nonlinear full-wave map.

**No quantitative external benchmark was possible because the published setup's geometry, normalization, noise model, parameterization, and units were not reproduced here; this is not independent validation.**

## Scope note

All results are finite-dimensional N=16 toy-model numbers only (m=48 real data rows, p=24 smooth coefficients, q_pose=18).  No continuum-limit, universal-trajectory, or external validation claim is made; the Born/FIM differences and slopes are internal analytic controls for this discrete experiment.

## Artifacts and digests

```text
ddf45031bd4b1077d4cf031b450d0e07f5b3978975c67a2320b4ff3c35fac2e1  [family15_born_control.py](src/family15_born_control.py)
645024f64fcc058b5920b2fab1da88663aa4d3b4c614dfcf45bc2d43d17995fa  [family15_born_control.json](/Volumes/migodam's-external-brain/Research/Inv_SLAM/experiments/idea_loops/loop_2026-09-03_16-25-02/experiment_pose_confounding_spectral_geometry/results/family15_born_control.json)
9909ca6273235d3090f9aee320bebaafb80c322976317b648d2a3062f7a3373c  [family15_born_control.png](/Volumes/migodam's-external-brain/Research/Inv_SLAM/experiments/idea_loops/loop_2026-09-03_16-25-02/experiment_pose_confounding_spectral_geometry/figures/family15_born_control.png)
ee3b4ac45ad296dba533a8ccefd107a6cc53e4357b7a768d07f7db90f3ee6585  src/helmholtz.py
762af3bdd33af1ab80563c7257b3805e7972f786fa88af1bf5e1b51a768526d7  src/family1_pilot.py
173d2b5d9f66be274c9a8422f4ac2c381bb91a60132b37c734748414cf3eb235  src/family2_algebraic_spine.py
ee7cac8b6233bcb0dbbb215df02452746dabf50774760910be5b16a0c0a40242  src/family4_frequency_trajectory.py
```

This report and the new result/figure files are additive; no existing source, result, figure, or note was modified.
