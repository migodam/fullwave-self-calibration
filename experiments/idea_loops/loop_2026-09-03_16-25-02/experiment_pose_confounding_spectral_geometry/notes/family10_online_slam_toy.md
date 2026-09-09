# Family 10: empirical NLS covariance toy vs linearized spectral predictions

Date: 2026-09-03T13:00:42.742177+00:00 (UTC)

## Scope

finite-dimensional toy only; no continuum, global-nonlinearity, real-SLAM-system, or production claim  
All fits are local (started at the true parameters); the report is a finite-dimensional empirical check and claims nothing about continuum or real SLAM systems.

## What was done

- Scene/config: see `results/family10_online_slam_toy.json` (`config`); coefficients `c0` chosen by least squares (`chi_true = S c0`).
- Linearized predictions computed from whitened/realified multi-frequency stacks: `K_IS`, `K_eff(alpha=1.0)`, P_known and P_free.
- NLS Monte Carlo: 3 trials per mode; known-pose and free-pose (pose-prior) fits via `scipy.optimize.least_squares` with analytic dense Jacobians.
- Born known/free converged trials: 3/3 and 3/3.  Full-wave known/free: 3/3 and 3/3.

## Key comparison numbers

Relative Frobenius deviations ||Cov_emp - P_pred||_F / ||P_pred||_F:

- Born known-pose: 0.899201
- Born free-pose:  42.236
- Full-wave known-pose: 0.826615
- Full-wave free-pose:  5.9198

Predicted retained DOF tr(K_eff K_IS^{-1}) and trace ratios:

- Born retained DOF 14.0265; trace ratio predicted 6.06832, empirical 424.931.
- Full-wave retained DOF 14.3655; trace ratio predicted 3.62905, empirical 21.196.

Three most-confounded direction rows (rho ascending):

| mode | dir | rho | predicted 1/rho | predicted exact ratio | empirical ratio |
|---|---|---|---|---|---|
| born | 0 | 0.000764505 | 1308.04 | 19.1762 | 202.84 |
| born | 1 | 0.00389207 | 256.933 | 7.11919 | 1077.61 |
| born | 2 | 0.0207158 | 48.2724 | 8.69062 | 58.1022 |
| full_wave | 0 | 0.000497834 | 2008.7 | 7.5052 | 3.34177 |
| full_wave | 1 | 0.00792522 | 126.179 | 2.83403 | 44.4166 |
| full_wave | 2 | 0.0273153 | 36.6095 | 2.53753 | 380.796 |

## Honest caveats

- The full-wave map is nonlinear in c; the Born mode is linear in c but still nonlinear in the pose perturbation dx (poses enter through Green functions/incident fields).  Both are compared against the same first-order prediction, so any mismatch is a measured nonlinear/sampling effect.
- Empirical sample covariance is 200-trial (fewer where the solver did not report success) and the most-confounded directions carry huge variance; Frobenius deviations therefore contain large sampling noise.
- `1/rho` is the classical retention/information-loss inflation factor and equals the variance ratio only when the generalized eigendirections also diagonalize both covariances (for example p=1 or K_IS proportional to I).  For the multi-direction c-space studied here the exact predicted ratio from the linearized covariances is (v^T P_free v)/(v^T P_known v); both are tabulated.
- The Born pose Jacobian helper in this script is locally implemented (no reused module exposes a Born B) and passed one centered finite-difference check; see JSON `finite_difference_self_check`.
- No claim is made that the fitter reached a global minimum; statuses and nfev/cost per trial are recorded in the JSON.

## Artifacts

- Results: `results/family10_online_slam_toy.json`
- Notes: `notes/family10_online_slam_toy.md`
- Figure: `figures/family10_cov_eigenvalues.png` (cov_eigenvalues)
- Figure: `figures/family10_confounded_direction_ratios.png` (direction_ratios)
- Figure: `figures/family10_born_vs_fullwave_deviation.png` (deviation)

## Artifacts and digests

```text
16cfe0ee07e7d91a5dfd50db214c0dfaf3996248e14c6dab4cdd9adf02f951a9  src/family10_online_slam_toy.py
ee3b4ac45ad296dba533a8ccefd107a6cc53e4357b7a768d07f7db90f3ee6585  src/helmholtz.py
762af3bdd33af1ab80563c7257b3805e7972f786fa88af1bf5e1b51a768526d7  src/family1_pilot.py
173d2b5d9f66be274c9a8422f4ac2c381bb91a60132b37c734748414cf3eb235  src/family2_algebraic_spine.py
4cd7e119794b20fd059638530ac893344d577901fac917d82c2ca2aaa75b4cc8  results/family10_online_slam_toy.json
9ac10cbbbec48c9faa184675db17cfad747c46d1b3f1deb4647085b2dc5d7d37  figures/family10_cov_eigenvalues.png
77ec7f7941c34e402a1bc29baff7d54e9cbc9dc8394fcfea561242aeadc7b907  figures/family10_confounded_direction_ratios.png
74afca0be3305cf0c1183bb6a89ae65077ed434c864264bf50a3679383405066  figures/family10_born_vs_fullwave_deviation.png
```
