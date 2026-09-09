# Family 4b: block-normalized frequency-diversity control

Date: 2026-09-03 (SGT; UTC stamp in `results/family4b_normalized_control.json`).
Experiment: `experiment_pose_confounding_spectral_geometry`.
Supplementary control for Family 4 (deferred in `results/family4_results.json`
and `notes/family4_report.md`).  **No Family 4 result or report file was
modified.**

## Exact command and runtime

```bash
.venv/bin/python src/family4b_normalized_control.py
```

Wall runtime: 0.17 s.  Platform:
macOS-26.6.2-arm64-arm-64bit, Python 3.12.13,
numpy 2.5.2, scipy 1.18.1.
Deterministic dense linear algebra; no RNG used.

## Scenario and normalization rule

Same scenario as Family 4: N=16, T=6 poses on the 90-degree arc at R=1.6,
n_rx=4, two-blob chi0, p=24 smooth basis, raw blocks f in
{1.0, 1.4, 1.8}; distinct stacks use {1.0, 1.4}; duplicate is block 1
repeated with c=2.0.  Rule:

A_f_n = A_f / ||A_f||_F and B_f_n = B_f / ||A_f||_F (same scalar for A and B within each frequency block; per-block linearized model rescaled consistently; duplicate blocks normalized the same way before c-scaling)

Directions u0..u2 are the three most-confounded generalized eigenvectors of
the **raw** single-frequency K_SLAM_1 w.r.t. K_IS_1
(u = V_A diag(1/s_A) w); they are scale-invariant and are reused unchanged
on the normalized stacks.  Direction agreement with the Family 4 stored u
vectors (abs cosine): 0: 1.000000000000, 1: 1.000000000000, 2: 1.000000000000.
Recomputed raw rho_single differs from the Family 4 reference by
0.00e+00, 0.00e+00, 0.00e+00.

## Claim status

**Does block-normalization change the Family 4 diversity conclusion?  NO.**
The normalized distinct stack still moves all three confounded directions
well above the 1e-4 gate (3/3),
the normalized duplicate movement remains at roundoff level
(max |5.232e-16| < 1e-10), and normalized PSD
monotonicity passes for all four prefix increments.

### C control rows (same u0,u1,u2; movement vs raw rho_single)

Raw Family 4 distinct movement: 1.212691e-01, 2.239175e-01, 4.768305e-01.
Normalized distinct movement: 1.059568e-01, 1.944408e-01, 3.814595e-01.

| dir | rho_single raw | rho_dist_norm | movement_dist_norm | z-resid dist norm | rho_dup_norm | movement_dup_norm | z-resid dup norm | gate >1e-4 | gate <1e-10 |
|---:|---:|---:|---:|---:|---:|---:|---:|:---:|:---:|
| 0 | 5.510027e-11 | 1.059568e-01 | 1.059568e-01 | 3.255100e-01 | 5.510025e-11 | -2.058506e-17 | 7.422945e-06 | PASS | PASS |
| 1 | 1.414337e-10 | 1.944408e-01 | 1.944408e-01 | 4.409545e-01 | 1.414342e-10 | 5.231665e-16 | 1.189259e-05 | PASS | PASS |
| 2 | 1.095607e-09 | 3.814595e-01 | 3.814595e-01 | 6.176241e-01 | 1.095607e-09 | 1.025875e-16 | 3.309995e-05 | PASS | PASS |

### A control: normalized PSD monotonicity rows

| comparison | min eig D | tol | gate |
|---|---:|---:|---|
| K_SLAM_norm F=1->2 | 2.238206e-10 | 1.000e-10 | PASS |
| K_eff_norm F=1->2 | 1.237880e-10 | 1.000e-10 | PASS |
| K_SLAM_norm F=2->3 | 3.094432e-10 | 1.000e-10 | PASS |
| K_eff_norm F=2->3 | 1.346684e-09 | 1.000e-10 | PASS |

### Block-norm record (raw vs normalized)

| block | ||A||_F raw | ||A||_F norm | ||B||_F raw | ||B||_F norm | B/A ratio raw | B/A ratio norm |
|---|---:|---:|---:|---:|---:|---:|
| f=1.0 | 3.21720872e-02 | 1.00000000e+00 | 2.00629404e-01 | 6.23613267e+00 | 6.236133 | 6.236133 |
| f=1.4 | 2.50030130e-02 | 1.00000000e+00 | 2.58631880e-01 | 1.03440285e+01 | 10.344029 | 10.344029 |
| f=1.8 | 1.67621688e-02 | 1.00000000e+00 | 1.65885763e-01 | 9.89643792e+00 | 9.896438 | 9.896438 |

## Cannot-establish section

* All claims are finite-dimensional discrete-model statements (N=16, p=24,
  the specific arc/chi0/blobs); no continuum-limit, universal-trajectory, or
  estimator claim is made.
* Block normalization rescales each raw block by its own ||A_f||_F, which
  equalizes A-block energy across frequencies but does not itself prescribe a
  physical noise covariance; it is one declared energy-matched control, not
  an optimal or correlated-noise whitening claim.
* The shared-z residual and rho are evaluated at only the three raw
  single-frequency confounded directions; other directions or larger stacks
  are not exhaustively checked (same boundary as Family 4 check C).

## Artifacts

- results: `results/family4b_normalized_control.json`
- figure: figures/family4b_normalized_shared_z.png
- this report: notes/family4b_normalized_control.md

Self-cell formula used: self-cell v2 corrected 2026-09-03: complete equal-area disk integral (adds the -1/k_b^2 lower-endpoint term omitted by v1).
