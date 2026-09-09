# Family 10: empirical NLS covariance toy vs linearized spectral predictions

Date: 2026-09-03T13:50:52.216913+00:00 (UTC)

## Scope

finite-dimensional toy only; no continuum, global-nonlinearity, real-SLAM-system, or production claim  
All fits are local (started at the true parameters); the report is a finite-dimensional empirical check and claims nothing about continuum or real SLAM systems.

## What was done

- Scene/config: see `results/family10_online_slam_toy_n120.json` (`config`); coefficients `c0` chosen by least squares (`chi_true = S c0`).
- Linearized predictions computed from whitened/realified multi-frequency stacks: `K_IS`, `K_eff(alpha=1.0)`, P_known and P_free.
- NLS Monte Carlo: 120/120 (born) and 120/120 (full-wave) completed trials; known-pose and free-pose (pose-prior) fits via `scipy.optimize.least_squares` with analytic dense Jacobians.  Wall-clock guard 2800 s per Monte-Carlo loop; hit: born=False, full_wave=False.
- Born known/free converged trials: 120/120 and 115/120.  Full-wave known/free: 120/120 and 119/120.

## Key comparison numbers

Relative Frobenius deviations ||Cov_emp - P_pred||_F / ||P_pred||_F:

- Born known-pose: 0.0778738
- Born free-pose:  13.4882
- Full-wave known-pose: 0.474486
- Full-wave free-pose:  4.04332

Predicted retained DOF tr(K_eff K_IS^{-1}) and trace ratios:

- Born retained DOF 14.0265; trace ratio predicted 6.06832, empirical 82.9206.
- Full-wave retained DOF 14.3655; trace ratio predicted 3.62905, empirical 11.0187.

Three most-confounded direction rows (rho ascending):

| mode | dir | rho | predicted 1/rho | predicted exact ratio | empirical ratio |
|---|---|---|---|---|---|
| born | 0 | 0.000764505 | 1308.04 | 19.1762 | 121.827 |
| born | 1 | 0.00389207 | 256.933 | 7.11919 | 82.3399 |
| born | 2 | 0.0207158 | 48.2724 | 8.69062 | 257.094 |
| full_wave | 0 | 0.000497834 | 2008.7 | 7.5052 | 4.46572 |
| full_wave | 1 | 0.00792522 | 126.179 | 2.83403 | 7.55445 |
| full_wave | 2 | 0.0273153 | 36.6095 | 2.53753 | 11.7863 |

## Full-wave finite-difference self-check

- Exact recipe (also stored in the JSON as `full_wave_fd_self_check.recipe`): at `theta0=(c0, dx=0)` the analytic free residual Jacobian `[[-A_stack, -B_stack], [0, sqrt(alpha) I_q]]` is compared with centered differences `(r(theta0+eps w) - r(theta0-eps w))/(2 eps)` at `eps=1e-6` for two full-c unit directions (random support on all map coordinates, zero dx), two full-dx unit directions, and one joint random direction; every direction is scored on the whole residual, the data rows, the prior rows, and each per-frequency row slice.
- Directional batch relative Frobenius error `||FD - J W||_F / ||J W||_F = 4.6595e-10`; max directional column error 2.63586e-09.

## Honest caveats

- The full-wave map is nonlinear in c; the Born mode is linear in c but still nonlinear in the pose perturbation dx (poses enter through Green functions/incident fields).  Both are compared against the same first-order prediction, so any mismatch is a measured nonlinear/sampling effect.
- Empirical sample covariance uses 120 completed Born and 120 completed full-wave trials (fewer where the solver did not report success); the most-confounded directions carry huge variance, so Frobenius deviations contain large sampling noise.
- `1/rho` is the classical retention/information-loss inflation factor and equals the variance ratio only when the generalized eigendirections also diagonalize both covariances (for example p=1 or K_IS proportional to I).  For the multi-direction c-space studied here the exact predicted ratio from the linearized covariances is (v^T P_free v)/(v^T P_known v); both are tabulated.
- The Born pose Jacobian helper in this script is locally implemented (no reused module exposes a Born B) and passed one centered finite-difference check; see JSON `finite_difference_self_check`.
- No claim is made that the fitter reached a global minimum; statuses and nfev/cost per trial are recorded in the JSON.

## Artifacts

- Results: `results/family10_online_slam_toy_n120.json`
- Notes: `notes/family10_online_slam_toy_n120.md`
- Figure: `figures/family10_cov_eigenvalues_n120.png` (cov_eigenvalues)
- Figure: `figures/family10_confounded_direction_ratios_n120.png` (direction_ratios)
- Figure: `figures/family10_born_vs_fullwave_deviation_n120.png` (deviation)

## Artifacts and digests

```text
44a976104b0fa329a5699ce8c048e2fa3ed125115465bd6807ad117fb303e45b  src/family10_online_slam_toy.py
ee3b4ac45ad296dba533a8ccefd107a6cc53e4357b7a768d07f7db90f3ee6585  src/helmholtz.py
762af3bdd33af1ab80563c7257b3805e7972f786fa88af1bf5e1b51a768526d7  src/family1_pilot.py
173d2b5d9f66be274c9a8422f4ac2c381bb91a60132b37c734748414cf3eb235  src/family2_algebraic_spine.py
3d4533ca828990ea1021db805b6f68765873d84870abe6c8abd2df960692efd8  results/family10_online_slam_toy_n120.json
39a18e0a5e4d21f22f83fc1de074750d4a56d46c5d3ad05f8003fc6cdd10c474  figures/family10_cov_eigenvalues_n120.png
8fb25b848a674a8b629629a451c791a1067d47d6828d739fce0c6a5d6acce472  figures/family10_confounded_direction_ratios_n120.png
d950510b31456c290d4ed78ca57ecb41ded3ccee7a45cace2a27799f17cb1250  figures/family10_born_vs_fullwave_deviation_n120.png
```
