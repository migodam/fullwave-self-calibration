# Family 4c: same-standoff, per-pose-energy-normalized trajectory control

Date: 2026-09-03 (SGT; UTC stamp in
`results/family4c_trajectory_control.json`).  Experiment:
`experiment_pose_confounding_spectral_geometry`.

Supplementary control requested by the Family 4 parent audit rule 7: hold
curved-path standoff fixed (R = 1.6) and equalize each pose's realified
A-block energy before projection, so angular coverage is compared after
removing standoff/radius and per-pose signal-strength scaling.  **No existing
src/results/notes file was modified.**

## Exact command and runtime

```bash
.venv/bin/python src/family4c_trajectory_control.py
```

Wall runtime: 0.147 s.  Platform:
macOS-26.6.2-arm64-arm-64bit, machine arm64,
Python 3.12.13, numpy 2.5.2,
scipy 1.18.1, matplotlib
3.11.1.  Deterministic dense linear algebra; no
RNG used.

## Source SHA-256

- `src/helmholtz.py`: `ee3b4ac45ad296dba533a8ccefd107a6cc53e4357b7a768d07f7db90f3ee6585`
- `src/family1_pilot.py`: `762af3bdd33af1ab80563c7257b3805e7972f786fa88af1bf5e1b51a768526d7`
- `src/family2_algebraic_spine.py`: `173d2b5d9f66be274c9a8422f4ac2c381bb91a60132b37c734748414cf3eb235`
- `src/family4_frequency_trajectory.py`: `ee7cac8b6233bcb0dbbb215df02452746dabf50774760910be5b16a0c0a40242`
- `src/family4c_trajectory_control.py`: `3a5c71eae7f2cf66c7094c321a6e9998760d16e937e4b1caee206978772445bb`

## Scenario and normalization rule

Same scene as Family 4/6: N=16, T=6, n_rx=4, receivers at the four body
offsets, tx at the body origin, two-blob chi0 (amp 0.3/0.5, sigma 0.09/0.07
at (-0.15,-0.12) and (0.18,0.14)), p=24 unit-column smooth Gaussian-RBF basis
(4x6 centres on [-0.3,0.3]^2, sigma_b=0.16), f=1.0 with k_b=2*pi.

for pose t, let block_t be its full realified row block (both Re and Im slices, 2*n_rx rows in A_pix_R/B_R; whiten_realify stacks sqrt(2)*[Re; Im]); n_t = ||A_pix_R[block_t]||_F; scale A_pix_R and B_R block rows by the SAME scalar 1/n_t; then A_s_norm = A_pix_R_norm @ S and B_R_norm = B_R.  Each pose's total A-block energy is normalized to 1 and standoff/signal-strength scaling is removed.

whiten_realify stacks 2*T*n_rx rows as [Re(A_hat) (T*n_rx); Im(A_hat) (T*n_rx)] with complex rows grouped by pose; pose t block rows = pose_block_indices(T,n_rx,t) = concat(arange(t*n_rx,(t+1)*n_rx), arange(T*n_rx+t*n_rx,T*n_rx+(t+1)*n_rx)), matching src/family4_parent_controls.py. A contiguous slice(t*n_rx,(t+1)*n_rx) alone would cover only the real-half rows and was therefore not used.

### Trajectory builders

| trajectory | builder | design length | min/max range |
|---|---|---:|---:|
| straight_same_mid | line p=(x, 1.6), x=linspace(-L/2, L/2, T), L=pi*1.6 (varying standoff; midpoint standoff 1.6) | 5.026548 | 1.677099 / 2.979353 |
| arc90_same_standoff | R=1.6, phi=linspace(-45, 45, T) degrees | 2.513274 | 1.600000 / 1.600000 |
| arc180_same_standoff | R=1.6, phi=linspace(-90, 90, T) degrees | 5.026548 | 1.600000 / 1.600000 |
| arc270_same_standoff | R=1.6, phi=linspace(-135, 135, T) degrees | 7.539822 | 1.600000 / 1.600000 |
| circle360_same_standoff | R=1.6, phi=linspace(0, 360, T, endpoint=False) degrees | 10.053096 | 1.600000 / 1.600000 |

## 1. Normalization verification

For every normalized trajectory, each pose's post-normalization total
A-block Frobenius norm must be 1 within 1e-12.

| trajectory | max |1 - ||A_block||_F| | pass <= 1e-12 |
|---|---:|---|
| straight_same_mid | 2.220e-16 | PASS |
| arc90_same_standoff | 2.220e-16 | PASS |
| arc180_same_standoff | 2.220e-16 | PASS |
| arc270_same_standoff | 1.110e-16 | PASS |
| circle360_same_standoff | 1.110e-16 | PASS |

Max deviation over all trajectories:
2.220e-16; pass =
True.

### Per-pose A/B block Frobenius norms (raw pre, normalized post) and scales

| trajectory | A_fro pre (per pose) | B_fro pre | A_fro post | B_fro post | scale = 1/A_fro |
|---|---|---|---|---|---|
| straight_same_mid | 1.496862e-02, 2.053776e-02, 2.717578e-02, 2.705171e-02, 2.042652e-02, 1.493986e-02 | 3.161497e-02, 6.789246e-02, 5.925286e-02, 6.756181e-02, 7.325302e-02, 5.396108e-02 | 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00 | 2.112083e+00, 3.305739e+00, 2.180355e+00, 2.497506e+00, 3.586172e+00, 3.611887e+00 | 66.806426, 48.690809, 36.797467, 36.966241, 48.955965, 66.935036 |
| arc90_same_standoff | 2.873276e-02, 2.852984e-02, 2.828399e-02, 2.822281e-02, 2.840155e-02, 2.854201e-02 | 9.347072e-02, 5.014030e-02, 3.908120e-02, 8.131162e-02, 1.012286e-01, 1.030301e-01 | 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00 | 3.253106e+00, 1.757469e+00, 1.381743e+00, 2.881060e+00, 3.564191e+00, 3.609768e+00 | 34.803478, 35.051022, 35.355694, 35.432332, 35.209341, 35.036072 |
| arc180_same_standoff | 2.843746e-02, 2.873836e-02, 2.839591e-02, 2.829340e-02, 2.852568e-02, 2.844447e-02 | 1.960986e-02, 9.531817e-02, 2.883473e-02, 9.407642e-02, 9.811892e-02, 2.986247e-02 | 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00 | 6.895783e-01, 3.316757e+00, 1.015453e+00, 3.325031e+00, 3.439669e+00, 1.049851e+00 | 35.164880, 34.796694, 35.216336, 35.343934, 35.056128, 35.156213 |
| arc270_same_standoff | 2.863027e-02, 2.848240e-02, 2.852984e-02, 2.840155e-02, 2.841759e-02, 2.870018e-02 | 9.863442e-02, 2.287733e-02, 5.014030e-02, 1.012286e-01, 5.035127e-02, 8.995543e-02 | 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00 | 3.445110e+00, 8.032092e-01, 1.757469e+00, 3.564191e+00, 1.771834e+00, 3.134315e+00 | 34.928067, 35.109397, 35.051022, 35.209341, 35.189468, 34.842983 |
| circle360_same_standoff | 2.821983e-02, 2.849295e-02, 2.867979e-02, 2.821437e-02, 2.859129e-02, 2.870121e-02 | 6.219971e-02, 9.215466e-02, 9.150715e-02, 5.887011e-02, 8.857697e-02, 8.756079e-02 | 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00 | 2.204114e+00, 3.234297e+00, 3.190649e+00, 2.086529e+00, 3.098041e+00, 3.050770e+00 | 35.436080, 35.096400, 34.867754, 35.442929, 34.975689, 34.841739 |

## 2. Hypothesis check on the four core trajectories (normalized)

Expected retained order (core): circle360 > arc180 > arc90 > straight.
Expected confusable order (ascending): straight < arc90 < arc180 < circle360.

Observed retained desc: arc90_same_standoff, straight_same_mid, arc180_same_standoff, circle360_same_standoff.
Observed confusable asc: arc90_same_standoff, straight_same_mid, arc180_same_standoff, circle360_same_standoff.
retained_consistent = False;
confusable_consistent = False.

**hypothesis_holds_for_four_core = False**
(not forced to pass).

Raw four-core context: retained desc = arc90_same_standoff, straight_same_mid, arc180_same_standoff, circle360_same_standoff;
hypothesis_holds_for_four_core (raw) =
False.

## 3. Ordering flips and the previous counterexample

Raw full retained order: arc90_same_standoff, straight_same_mid, arc180_same_standoff, circle360_same_standoff, arc270_same_standoff.
Normalized full retained order: arc90_same_standoff, straight_same_mid, arc180_same_standoff, circle360_same_standoff, arc270_same_standoff.
retained_order_flipped = False.

Raw full confusable desc: arc270_same_standoff, circle360_same_standoff, arc180_same_standoff, straight_same_mid, arc90_same_standoff.
Normalized full confusable desc:
arc270_same_standoff, circle360_same_standoff, arc180_same_standoff, straight_same_mid, arc90_same_standoff.
confusable_order_flipped = False.

Family 4 raw counterexample retained order:
arc90, straight, circle360, arc180.  Family 4c raw core
order: arc90_same_standoff, straight_same_mid, arc180_same_standoff, circle360_same_standoff.  Family 4c
normalized core order:
arc90_same_standoff, straight_same_mid, arc180_same_standoff, circle360_same_standoff.
Family-4-to-family4c-raw core change =
True; normalized gate =
**False**;
status: persists as a counterexample under the normalized control.

Family 4's raw order is a scenario counterexample under equal design length with different radii; family4c raw removes the radius/standoff difference between curved paths (all R=1.6) while retaining a varying-standoff straight path and changing path length with angular coverage.  Only the normalized four-core comparison is the requested hypothesis gate.

## Raw and normalized metric tables


#### Raw table (same-standoff geometry, un-normalized)

| trajectory | retained | confusable | theta_min deg | log_volume | rho<1e-6 | rho_min | r_A | r_B | design L | open polyline | min range | max range |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| straight_same_mid | 8.353517 | 15.646483 | 9.887569e-05 | -204.0601 | 6 | 2.9781e-12 | 24 | 18 | 5.026548 | 5.026548 | 1.677099 | 2.979353 |
| arc90_same_standoff | 9.407162 | 14.592838 | 4.253070e-04 | -167.7461 | 5 | 5.5101e-11 | 24 | 18 | 2.513274 | 2.502951 | 1.600000 | 1.600000 |
| arc180_same_standoff | 7.999302 | 16.000698 | 3.193125e-04 | -175.8303 | 6 | 3.1059e-11 | 24 | 18 | 5.026548 | 4.944272 | 1.600000 | 1.600000 |
| arc270_same_standoff | 7.536362 | 16.463638 | 4.018091e-04 | -189.3070 | 6 | 4.9181e-11 | 24 | 18 | 7.539822 | 7.263848 | 1.600000 | 1.600000 |
| circle360_same_standoff | 7.703633 | 16.296367 | 3.036807e-04 | -199.5708 | 6 | 2.8092e-11 | 24 | 18 | 10.053096 | 8.000000 | 1.600000 | 1.600000 |


#### Normalized table

| trajectory | retained | confusable | theta_min deg | log_volume | rho<1e-6 | rho_min | r_A | r_B | design L | open polyline | min range | max range |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| straight_same_mid | 8.269820 | 15.730180 | 9.947104e-05 | -203.9637 | 6 | 3.0140e-12 | 24 | 18 | 5.026548 | 5.026548 | 1.677099 | 2.979353 |
| arc90_same_standoff | 9.408661 | 14.591339 | 4.270746e-04 | -167.6930 | 5 | 5.5560e-11 | 24 | 18 | 2.513274 | 2.502951 | 1.600000 | 1.600000 |
| arc180_same_standoff | 7.997790 | 16.002210 | 3.216665e-04 | -175.8197 | 6 | 3.1519e-11 | 24 | 18 | 5.026548 | 4.944272 | 1.600000 | 1.600000 |
| arc270_same_standoff | 7.535343 | 16.464657 | 4.032271e-04 | -189.3062 | 6 | 4.9528e-11 | 24 | 18 | 7.539822 | 7.263848 | 1.600000 | 1.600000 |
| circle360_same_standoff | 7.701223 | 16.298777 | 3.053205e-04 | -199.5498 | 6 | 2.8397e-11 | 24 | 18 | 10.053096 | 8.000000 | 1.600000 | 1.600000 |

## Path geometry record (same in both tables)

| trajectory | continuous design L | open polyline | closed polyline if counted | closure counted for circle | min range | max range |
|---|---:|---:|---:|:---:|---:|---:|
| straight_same_mid | 5.026548 | 5.026548 | 10.053096 | False | 1.677099 | 2.979353 |
| arc90_same_standoff | 2.513274 | 2.502951 | 4.765693 | False | 1.600000 | 1.600000 |
| arc180_same_standoff | 5.026548 | 4.944272 | 8.144272 | False | 1.600000 | 1.600000 |
| arc270_same_standoff | 7.539822 | 7.263848 | 9.526590 | False | 1.600000 | 1.600000 |
| circle360_same_standoff | 10.053096 | 8.000000 | 9.600000 | True | 1.600000 | 1.600000 |

## Smooth A and B norm record (matrix_norm_stats)

| trajectory | ||A_s||_F raw | ||A_s_norm||_F | ||B_R||_F raw | ||B_R_norm||_F | ||A_s_norm||_2 | ||B_R_norm||_2 |
|---|---:|---:|---:|---:|---:|---:|
| straight_same_mid | 2.390104e-02 | 1.131280e+00 | 1.481965e-01 | 7.230951e+00 | 6.205971e-01 | 3.610955e+00 |
| arc90_same_standoff | 3.217209e-02 | 1.130697e+00 | 2.006294e-01 | 7.043792e+00 | 6.403748e-01 | 3.608476e+00 |
| arc180_same_standoff | 3.183231e-02 | 1.117895e+00 | 1.722527e-01 | 6.041265e+00 | 5.677440e-01 | 3.437726e+00 |
| arc270_same_standoff | 3.199692e-02 | 1.121480e+00 | 1.834138e-01 | 6.424127e+00 | 5.408808e-01 | 3.562678e+00 |
| circle360_same_standoff | 3.208995e-02 | 1.126705e+00 | 1.992820e-01 | 6.982671e+00 | 5.978228e-01 | 3.231214e+00 |

## Cannot-establish section

1. All claims are finite-dimensional and model-specific (N=16, p=24 smooth basis, T=6, R=1.6, this two-blob scene); no continuum-limit, universal-trajectory, or estimator claim.
2. Per-pose normalization is a declared gain/whitening control on the linearized blocks; it does not prescribe a physical noise covariance and is not a claim about optimal SNR whitening under correlated noise.
3. Same standoff fixes curved-path range but changes measurement count/budget with angular coverage (continuous design length and sampled polyline lengths differ), and the straight path still varies standoff, so even the normalized five-trajectory set is not a pure angular-coverage isolation.
4. No universal trajectory-ordering theorem is asserted; all ordering statements are observations on this discrete model and normalization/geometry convention.

## Artifacts

- results: `results/family4c_trajectory_control.json`
- figure: `figures/family4c_trajectory_control.png`
- this report: `notes/family4c_trajectory_control.md`
