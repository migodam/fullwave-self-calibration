# Family 13: rank-boundary autopsy (ring-scene c2 residual)

Date: 2026-09-03T16:54:54.920112+00:00.  Experiment: `experiment_pose_confounding_spectral_geometry`.

Family 13 autopsies the Family-9 ring-scene machine-rank boundary: c2 rank identity r_KIS - r_KSL = r_A + r_B - r_AB records residual 1 on the ring scene and 0 on the two-blob reference.  It reuses (read-only) the standard builders and the family-5 algebraic falsifier; no existing source/result/figure/note is modified.

## Exact command and runtime

```bash
.venv/bin/python src/family13_rank_autopsy.py
```

Wall runtime: 0.575 s.  Platform: macOS-26.6.2-arm64-arm-64bit.  Python 3.12.13, numpy 2.5.2, scipy 1.18.1, matplotlib 3.11.1.  CPU-only (Apple Silicon).

Source SHA-256 (this script): `531ef1957dd0e215f17e2ac46a5e662d83d0040b983407761a5691f2d523fdfa`

Reused-module SHA-256:
* `helmholtz.py` `ee3b4ac45ad296dba533a8ccefd107a6cc53e4357b7a768d07f7db90f3ee6585`
* `family1_pilot.py` `762af3bdd33af1ab80563c7257b3805e7972f786fa88af1bf5e1b51a768526d7`
* `family2_algebraic_spine.py` `173d2b5d9f66be274c9a8422f4ac2c381bb91a60132b37c734748414cf3eb235`
* `family4_frequency_trajectory.py` `ee7cac8b6233bcb0dbbb215df02452746dabf50774760910be5b16a0c0a40242`
* `family5_parent_generalized.py` `fbdc8c2f11f6e3503f30f637ec65570062b2c3207d849e8bde1aa62236cd3e4d`
* `family9_replications_ablation.py` `a003d3b3067995d5201a9a892c440880669cc27bf3330980bebc344f8db78803`

Figure SHA-256: `4eddd1fa5a8239ed3304aa74533e016b95bbb488798dd4df12e4348367cc8dbc`.

## Reproduction row (ring scene, f=1.0)

Family-2 machine ranks: r_A=24, r_B=18, r_AB=42, r_KIS=24, r_KSL=23; lhs=1, rhs=0, residual=1.

Stored family-9 ring row: r_KIS=24, r_KSL=23, r_A=24, r_B=18, r_AB=42, residual=1.  Row match: **True**.

The ring chi is the raw `0.5*exp(-((|r|-0.25)/0.05)^2)` grid vector (no L2 rescaling); a global scene scale leaves every c2 rank unchanged, which is why the raw ring row equals the family-9 row.

## A. Raw spectra

Full descending singular-value arrays for A_s, B_R, AB, K_IS, K_SLAM(direct), C, C^T C and R_op, plus ascending/descending eigenvalues of K_SLAM and R_op, are embedded in `results/family13_rank_autopsy.json` (`scenes.<scene>.raw_spectra`).

Key relative minima:

* ring: sigma_24/sigma_1 C=4.74995e-08, K_SLAM=2.25753e-15, R_op=9.96215e-11.
* two_blob: sigma_24/sigma_1 C=3.90533e-07, K_SLAM=1.52533e-13, R_op=5.51003e-11.

## B. Predeclared relative-tolerance grid

tol_rel grid fixed in CONFIG before computing: 1e-16 ... 1e-8 (nine values); rank = count(sigma > tol_rel*sigma_1).

Ring:

| tol_rel | r_A | r_B | r_AB | r_KIS | r_KSL | r_C | r_CtC | r_Rop | c2 res. |
| ------- | --- | --- | ---- | ----- | ----- | --- | ----- | ----- | ------- |
| 1e-16   | 24  | 18  | 42   | 24    | 24    | 24  | 24    | 24    | 0       |
| 1e-15   | 24  | 18  | 42   | 24    | 24    | 24  | 24    | 24    | 0       |
| 1e-14   | 24  | 18  | 42   | 24    | 23    | 24  | 23    | 24    | 1       |
| 1e-13   | 24  | 18  | 42   | 24    | 23    | 24  | 23    | 24    | 1       |
| 1e-12   | 24  | 18  | 42   | 24    | 21    | 24  | 21    | 24    | 3       |
| 1e-11   | 24  | 18  | 42   | 23    | 19    | 24  | 19    | 24    | 4       |
| 1e-10   | 24  | 18  | 42   | 22    | 18    | 24  | 18    | 23    | 4       |
| 1e-09   | 24  | 18  | 42   | 20    | 17    | 24  | 17    | 21    | 3       |
| 1e-08   | 24  | 18  | 41   | 19    | 17    | 24  | 17    | 20    | 1       |

Two-blob reference:

| tol_rel | r_A | r_B | r_AB | r_KIS | r_KSL | r_C | r_CtC | r_Rop | c2 res. |
| ------- | --- | --- | ---- | ----- | ----- | --- | ----- | ----- | ------- |
| 1e-16   | 24  | 18  | 42   | 24    | 24    | 24  | 24    | 24    | 0       |
| 1e-15   | 24  | 18  | 42   | 24    | 24    | 24  | 24    | 24    | 0       |
| 1e-14   | 24  | 18  | 42   | 24    | 24    | 24  | 24    | 24    | 0       |
| 1e-13   | 24  | 18  | 42   | 24    | 24    | 24  | 24    | 24    | 0       |
| 1e-12   | 24  | 18  | 42   | 24    | 23    | 24  | 23    | 24    | 1       |
| 1e-11   | 24  | 18  | 42   | 23    | 22    | 24  | 22    | 24    | 1       |
| 1e-10   | 24  | 18  | 42   | 21    | 20    | 24  | 20    | 23    | 1       |
| 1e-09   | 24  | 18  | 42   | 21    | 19    | 24  | 19    | 22    | 2       |
| 1e-08   | 24  | 18  | 41   | 20    | 18    | 24  | 18    | 20    | 1       |

## C. Stable-rank / effective-rank diagnostics

Ring:

| matrix        | eff rank | stable rank | sigma_24/1        | sigma_23/24 | smallest rel. (top 4 of 6)             |
| ------------- | -------- | ----------- | ----------------- | ----------- | -------------------------------------- |
| A_s           | 6.989    | 2.132       | 1.13e-06          | 4           | 0.000184, 6.11e-05, 2.54e-05, 1.11e-05 |
| B_R           | 6.302    | 5.98        | 1.86e-05 (last/1) | n/a         | 0.000332, 0.000332, 0.000203, 0.000203 |
| AB            | 6.933    | 5.926       | 0.000105          | 1.6         | 2.55e-07, 1.82e-07, 8.3e-08, 6.86e-08  |
| K_IS          | 3.675    | 1.237       | 1.27e-12          | 16          | 3.37e-08, 3.73e-09, 6.45e-10, 1.23e-10 |
| K_SLAM_direct | 1.775    | 1.21        | 2.26e-15          | 58          | 3.19e-11, 8.3e-12, 2.81e-12, 9.57e-13  |
| C             | 2.216    | 1.465       | 4.75e-08          | 7.62        | 5.65e-06, 2.88e-06, 1.68e-06, 9.78e-07 |
| C^T C         | 1.775    | 1.21        | 2.25e-15          | 58.2        | 3.19e-11, 8.3e-12, 2.81e-12, 9.57e-13  |
| R_op          | 10.12    | 9.397       | 9.96e-11          | 4.56        | 2.6e-06, 2.8e-08, 8.08e-09, 7.27e-10   |

Two-blob:

| matrix        | eff rank | stable rank | sigma_24/1       | sigma_23/24 | smallest rel. (top 4 of 6)             |
| ------------- | -------- | ----------- | ---------------- | ----------- | -------------------------------------- |
| A_s           | 7.518    | 3.105       | 1.89e-06         | 2.61        | 0.000172, 0.000128, 3.45e-05, 6.92e-06 |
| B_R           | 6.314    | 3.795       | 0.00655 (last/1) | n/a         | 0.0179, 0.0161, 0.0108, 0.0103         |
| AB            | 7.083    | 3.855       | 0.000552         | 1.43        | 5.21e-07, 1.68e-07, 9.42e-08, 5.6e-08  |
| K_IS          | 5.015    | 1.923       | 3.57e-12         | 6.83        | 2.94e-08, 1.64e-08, 1.19e-09, 4.79e-11 |
| K_SLAM_direct | 2.654    | 1.36        | 1.53e-13         | 17.9        | 1.42e-09, 1.48e-10, 4.66e-11, 1.62e-11 |
| C             | 4.036    | 1.9         | 3.91e-07         | 4.23        | 3.77e-05, 1.22e-05, 6.83e-06, 4.03e-06 |
| C^T C         | 2.654    | 1.36        | 1.53e-13         | 17.9        | 1.42e-09, 1.48e-10, 4.66e-11, 1.62e-11 |
| R_op          | 9.917    | 8.923       | 5.51e-11         | 2.57        | 2.91e-06, 5.53e-08, 5.73e-09, 1.1e-09  |

Smallest eigvalsh eigenvalue of ring K_SLAM: 2.33309e-19 (positive=True).

## D. Backward residuals

Ring:

| quantity                                                           | value             |
| ------------------------------------------------------------------ | ----------------- |
| rel Fro ||K_SLAM - C^T C|| / max(||K_SLAM||F, ||C||F^2)            | 5.00069e-16       |
| rel Fro ||R_op - V^T diag(1/s) K_SLAM V diag(1/s)|| / max(F norms) | 1.51994e-07       |
| R_op definition residual ||R_op - Q_A^T P Q_A||F / ||R_op||F       | 0.0 (same object) |
| ||P^2-P||F / ||P||F (idempotency)                                  | 3.59945e-16       |
| ||Z^T Z - I||F                                                     | 1.94626e-15       |
| c2 lhs (r_KIS - r_KSL)                                             | 1                 |
| c2 rhs (r_A + r_B - r_AB)                                          | 0                 |
| c2 residual                                                        | 1                 |
| family2 rank tol(K_SLAM)                                           | 5.50553e-19       |
| smallest retained sigma(K_SLAM) - tol                              | 1.29794e-17       |
| largest dropped sigma(K_SLAM) below tol                            | 3.17324e-19       |

Two-blob:

| quantity                                                           | value             |
| ------------------------------------------------------------------ | ----------------- |
| rel Fro ||K_SLAM - C^T C|| / max(||K_SLAM||F, ||C||F^2)            | 3.17762e-14       |
| rel Fro ||R_op - V^T diag(1/s) K_SLAM V diag(1/s)|| / max(F norms) | 8.26455e-09       |
| R_op definition residual ||R_op - Q_A^T P Q_A||F / ||R_op||F       | 0.0 (same object) |
| ||P^2-P||F / ||P||F (idempotency)                                  | 3.45535e-16       |
| ||Z^T Z - I||F                                                     | 1.67139e-15       |
| c2 lhs (r_KIS - r_KSL)                                             | 0                 |
| c2 rhs (r_A + r_B - r_AB)                                          | 0                 |
| c2 residual                                                        | 0                 |
| family2 rank tol(K_SLAM)                                           | 1.10159e-20       |
| smallest retained sigma(K_SLAM) - tol                              | 3.04288e-19       |
| largest dropped sigma(K_SLAM) below tol                            | n/a               |

## E. Verdict

### ring

Supported explanation: `b_tolerance_classification_of_a_near_null_singular_value`.

Supported explanation: (b).  On the ring scene C=(I-ZZ^T)A_s has full column rank at every predeclared relative tolerance (sigma_24(C)/sigma_1(C) is far above the grid), while the directly formed K_SLAM=A^T P A has a near-null singular value sigma_24(K_SLAM)/sigma_1(K_SLAM) = [sigma_24(C)/sigma_1(C)]^2 that falls below tolerance on part of the grid.  This is a numerical rank classification of a near-null direction, not a discontinuous exact rank change; K_SLAM and C^T C agree to machine-level relative Frobenius error, so unstable direct subtraction between large separately computed terms is not the mechanism.

Supporting numbers:

* sigma_24_over_sigma_1_C: 4.749949680753936e-08
* sigma_24_over_sigma_1_K_SLAM_direct: 2.2575339767240223e-15
* sigma_24_over_sigma_1_R_op: 9.962146194239543e-11
* C_full_rank_all_grid_tol_rel: True
* K_SLAM_direct_rank_drops_some_grid_tol_rel: True
* residual_nonzero_on_some_grid_tol_rel: True
* rel_fro_K_SLAM_vs_C^T_C: 5.000691004738626e-16
* family2_machine_residual: 1

No theorem is claimed: the verdict above is supported by the computed finite-dimensional numbers only.

### two_blob

Supported explanation: `none_family2_machine_residual_zero`.

No rank boundary under the declared family-2 machine rule: the two-blob c2 residual is 0 and K_SLAM's smallest singular value ratio (1.5e-13) sits well above the family-2 relative threshold 24*eps ~= 5.3e-15.  As in the ring scene, coarse relative tolerances near/above 1e-12 begin shedding the near-null K_SLAM modes (that is tolerance classification, not a discontinuous rank event), but no such shedding occurs at the declared machine rule, and the exact identity is not in question here.

Supporting numbers:

* sigma_24_over_sigma_1_C: 3.90532992971459e-07
* sigma_24_over_sigma_1_K_SLAM_direct: 1.5253257714905237e-13
* sigma_24_over_sigma_1_R_op: 5.510029054500707e-11
* C_full_rank_all_grid_tol_rel: True
* K_SLAM_direct_rank_drops_some_grid_tol_rel: True
* residual_nonzero_on_some_grid_tol_rel: True
* rel_fro_K_SLAM_vs_C^T_C: 3.1776196516009236e-14
* family2_machine_residual: 0

No theorem is claimed: the verdict above is supported by the computed finite-dimensional numbers only.

## F. Engineered algebraic falsifier (two_blob, f=1.0)

B(t) = U diag(s_1,...,s_17,|t| s_18) V^T, reusing `family5_parent_generalized.engineered_rank_event`/`range_projector`/`retention_decomposition` read-only.  This is an **algebraic control** for the discontinuous no-prior projector at exact rank loss; it is not a physical ring event and says nothing about the ring scene.

| t     | rank(B(t)) | retained mass | ||P(t)-P(0)||_2 | rel Fro info jump / ||K_IS||F |
| ----- | ---------- | ------------- | --------------- | ----------------------------- |
| 0     | 17         | 10.40528058   | 0               | 0                             |
| 1e-06 | 18         | 9.40716187    | 1               | 0.139957                      |
| 1e-03 | 18         | 9.40716187    | 1               | 0.139957                      |
| 1e+00 | 18         | 9.40716187    | 1               | 0.139957                      |

Retained mass jumps from 9.40716187 (t>0, full-rank B) to 10.40528058 (t=0) while the no-prior projector jumps by operator norm 1 between t=0 and every t>0 row (family-5 reference with t=1e-8: 1).

## Scope

All numbers are finite-dimensional N=16 dense linear algebra on the whitened/realified smooth p=24 model at f=1.0 with the standard pose arc and receiver set.  The ring residual of 1 is a numerical rank classification boundary and is **not evidence against the exact c2 identity**; C=(I-ZZ^T)A_s remains full rank at every predeclared grid tolerance, and K_SLAM agrees with C^T C at machine-level relative Frobenius error.  No continuum/transversality or theorem claim is made.

## Artifacts

* [results/family13_rank_autopsy.json](results/family13_rank_autopsy.json)
* [figures/family13_rank_autopsy.png](figures/family13_rank_autopsy.png)
* [notes/family13_rank_autopsy.md](notes/family13_rank_autopsy.md)
