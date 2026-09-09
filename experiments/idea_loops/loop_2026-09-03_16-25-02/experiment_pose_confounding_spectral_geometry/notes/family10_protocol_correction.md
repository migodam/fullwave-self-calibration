# Family 10: fixed-pose sandwich (`P_samp`) targets for the nonlinear stacks

Date: 2026-09-03

Follow-up to `results/family10_online_slam_toy_n120.json` and the affine
control in `results/family10_affine_control.json` / `notes/family10_affine_control.md`.

## Protocol statement

The n120 Monte Carlo draws only data noise; the hidden pose stays fixed at
`dx = 0` on every trial. For that replication scheme the correct linearized
covariance of the map component of the penalized free-pose NLS estimator is
the sandwich

    P_samp = [H^-1 H_data H^-1]_{1:p,1:p},
    H_data = [[K_IS, A^T B], [(A^T B)^T, B^T B]],
    H      = H_data + diag(0_p, alpha I_q),        alpha = 1.0,

not the Bayesian marginal `P_free = K_eff^-1`. `K_eff^-1` is the posterior /
prior-consistent map marginal, and is realized empirically only when each
trial draws `dx_i ~ N(0, alpha^-1 I_q)` and generates data at that drawn pose.

This note recomputes the finite-dimensional stacks and both targets for the
nonlinear Born and full-wave models, then compares them with the stored n120
empirics.

## Computation provenance

All stacks come from the import-safe builders in `src/family10_online_slam_toy.py`
(no `main` was called and no heavy nonlinear fit was rerun):

- `build_scene(CONFIG)` and `build_model_blocks(scene, CONFIG, mode)`, mode in
  `{"born", "full_wave"}`;
- same scene/config/frequencies `{1.0, 1.4, 1.8}`, alpha = 1.0, snr = 100.0;
- whitened/realified stacked shapes `A_stack` (144 x 24), `B_stack` (144 x 18);
- theory reuse via `linearized_predictions(A_stack, B_stack, alpha)` and the
  affine-control `_block_quantities` sandwich construction.

Reproducibility cross-check against the stored JSON theory:

| check | born | full_wave |
|---|---:|---:|
| max abs diff, `rho_asc` vs JSON | 0.0 | 0.0 |
| max abs diff, `P_known` eigenvalues vs JSON | 0.0 | 0.0 |
| max abs diff, `P_free` eigenvalues vs JSON | 0.0 | 0.0 |

The rebuilt stacks therefore match the n120 run bit-for-bit at the level of
the stored theory spectra. Runtime of the linear algebra was under one second;
no existing result file was modified.

## Corrected comparison table

| mode | empirical fixed-pose trace ratio (n120) | `P_free` target = tr(P_free)/tr(P_known) | `P_samp` target = tr(P_samp)/tr(P_known) | interpretation |
|---|---:|---:|---:|---|
| born | 82.92 | 6.068 | 3.889 | Corrected target is **36% lower** than the old `P_free` target; nonlinearity still inflates the free-pose covariance by a factor of roughly **21.3** (empirical ratio / P_samp ratio) or **22.5** (tr(Cov_free)/tr(P_samp)) |
| full_wave | 11.019 | 3.629 | 2.755 | Corrected target is **24% lower** than the old target; nonlinearity still inflates the free-pose covariance by a factor of roughly **4.00** (empirical ratio / P_samp ratio) or **5.78** (tr(Cov_free)/tr(P_samp)) |

Trace magnitudes:

| mode | tr(P_known) | tr(P_free) | tr(P_samp) | tr(Cov_known, emp) | tr(Cov_free, emp) | tr(Cov_free)/tr(P_samp) |
|---|---:|---:|---:|---:|---:|---:|
| born | 28930.45 | 175559.33 | 112504.76 | 30585.35 | 2536156.15 | 22.54 |
| full_wave | 6444.68 | 23388.06 | 17755.55 | 9316.03 | 102650.85 | 5.78 |

## Generalized eigen-directions (3 smallest `rho` of (K_eff, K_IS))

Directions `v` are the Euclidean-normalized generalized eigenvectors from the
Cholesky `K_IS = L L^T` construction. `1/rho` is reported only as the naive
spectral inflation; the exact variance ratios are the quadratic forms
`(v^T P v)/(v^T P_known v)`.

### born

| dir | rho | predicted 1/rho | exact `P_free` ratio | exact `P_samp` ratio | empirical n120 ratio |
|---|---:|---:|---:|---:|---:|
| 0 | 7.645e-4 | 1308.0 | 19.18 | 11.10 | 121.83 |
| 1 | 3.892e-3 | 256.9 | 7.119 | 4.417 | 82.34 |
| 2 | 2.072e-2 | 48.27 | 8.691 | 4.675 | 257.09 |

### full_wave

| dir | rho | predicted 1/rho | exact `P_free` ratio | exact `P_samp` ratio | empirical n120 ratio |
|---|---:|---:|---:|---:|---:|
| 0 | 4.978e-4 | 2008.7 | 7.505 | 4.443 | 4.466 |
| 1 | 7.925e-3 | 126.2 | 2.834 | 2.370 | 7.554 |
| 2 | 2.732e-2 | 36.61 | 2.538 | 2.202 | 11.786 |

The JSON does store empirical per-direction ratios (`direction_rows` →
`empirical_inflation_ratio`), reproduced in the last column. For Born, the
largest nonlinear excess is along direction 2 (empirical 257x vs P_samp 4.7x),
while the smallest-rho direction shows a 11x excess (empirical 121.8 vs P_samp
11.1). For full_wave the smallest-rho direction is the only one near target
(4.47 empirical vs 4.44 P_samp), while directions 1-2 exceed P_samp by 3.2x and
5.4x.

## Known-pose match numbers (stored in the n120 JSON)

| mode | known success | free success | known rel_fro | known rel_spectral | tr(Cov_known)/tr(P_known) |
|---|---:|---:|---:|---:|---:|
| born | 120 / 120 | 115 / 120 | 0.0779 | 0.0783 | 1.057 |
| full_wave | 120 / 120 | 119 / 120 | 0.4745 | 0.4535 | 1.446 |

So the fixed-pose map estimator agrees with `P_known` closely for Born and
only moderately for full_wave (whose data model is not affine in `c`; the
stored `y0 - A c0` mismatch is ~0.446 relative norm in that mode).

## Interpretation and uncertainty

1. `K_eff^{-1}` (old target, ratio 6.07 / 3.63) is the Bayesian marginal and
   is valid when true poses are drawn from the prior per trial. Under the
   fixed-true-pose protocol the linearized frequentist target is `P_samp`
   (ratio 3.89 / 2.76), consistent with the Born affine control
   (`trace_ratio_exact_fixed_pose = 3.889`).
2. Correcting the target makes the Born discrepancy worse, not better: the
   stored `free_rel_fro` against `P_free` is 13.49, while `Cov_free` against
   `P_samp` has rel_fro 22.0; the corresponding full_wave numbers are 4.04 vs
   5.40.
3. Nonlinearity therefore still inflates the fixed-pose free covariance:
   Born by ~22.5x in trace (or ~21.3x in the ratio-to-ratio sense) and
   full_wave by ~5.8x (or ~4.0x), after the `P_samp` correction.
4. Born's empirical free covariance is heavily tail-dominated: successful
   free fits include per-trial `c_error_l2` values up to 7.4e3 (several
   trials in the 4e3-7.4e3 range), versus a 4.9e2 max for known-pose fits.
   With only 115-119 free successes, the factor estimates carry wide Monte
   Carlo uncertainty even though the qualitative conclusion (empirical >>
   P_samp) is robust.

No result files were modified. Only this note was created.
