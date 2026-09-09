# Family 5: trajectory sensitivity, robust surrogate, and rank-event stability

Date: 2026-09-03 (SGT; UTC stamp in `results/family5_sensitivity.json`).
Experiment: `experiment_pose_confounding_spectral_geometry`.

Family 5 reuses, without rerunning, the corrected Family 1 forward/Jacobian
builders (`src/helmholtz.py`, `src/family1_pilot.py`), the Family 2 smooth
p=24 RBF basis (`src/family2_algebraic_spine.py`), and the Family 4 K_eff /
whitening conventions.  All numbers are for the discrete N=16 whitened/
realified model; no continuum-limit, recovery, or global-SE(2) claim is made.

## Exact command and runtime

```bash
.venv/bin/python src/family5_sensitivity.py
```

Wall runtime: 18.04 s.  Platform:
macOS-26.6.2-arm64-arm-64bit, Python 3.12.13,
numpy 2.5.2, scipy 1.18.1,
matplotlib 3.11.1.

Scenario: N=16, T=6 arc poses, n_rx=4, q=18 pose parameters, two-blob chi0,
p=24 smooth basis, alpha=1.0.  Pose perturbation norm is Euclidean on R^18
(dual norm Euclidean).  Rank tolerance:
`tol(M) = max(M.shape) * eps_machine * sigma_1(M)`.

## Tolerances used (Family 5 gates)

| gate | value |
|---|---:|
| P1 rel_full at eps=1e-3 | < 1e-4 |
| P1 rel_frozen at eps=1e-3 | > 1e-4 |
| P1 slope of rel_full over clean window | >= 1.0 (expect ~2) |
| P2 rel_full_eig at eps=1e-3 (gapped) | < 1e-3 |
| P2 rel_frozen_eig / rel_full_eig at eps=1e-3 | > 10x |
| P2 chosen gap | any < 1e-6 flags that eigenvalue as non-simple and excludes it from the P2 gate (others asserted separately) |
| P3 err_lam_adv at eps=1e-2 | < 5e-2 and decreasing with eps |
| P3 sampled_min | >= lam_adv - 1e-12*max(1,|lam_r|) |
| P4 violations of empirical Weyl lower bound | == 0 (empirical-only) |
| P5 | recorded; no forced pass |

## Claim-status table

| check | status | executed comparison | key numbers |
|---|---|---|---|
| P1 derivative formula vs FD, full (includes dW) | **PASS** | rel_full at 1e-3, clean slope | rel_full(1e-3)=6.247e-07, slope=1.91 |
| P1 frozen-nuisance approximation is clearly worse | **PASS** | rel_frozen at 1e-3 > 1e-4 | rel_frozen(1e-3)=4.442e-02 |
| P1 dW term dominance | observed | ||A0^T dW A0||_F / ||DK_full||_F | 4.442e-02 |
| P2 simple-eigenvalue first-order prediction (full, gapped) | **PASS** | rel_full_eig(1e-3) and slope for well-gapped chosen | full slopes ~2.00 for idx 23 and 22; rel_full_eig(1e-3) in table below |
| P2 frozen prediction clearly worse (gapped) | **PASS** | ratio frozen/full at 1e-3 for gapped | idx 23 ratio=39.4x, idx 22 ratio=14.8x |
| P2 gap flags (median in dense cluster) | flagged | min adjacent gap < 1e-6 invalidates simple-eigenvalue condition there | {'any_chosen_gap_lt_1e-6': True, 'invalidated_indices': [12]}; gated eigenvalues pass, median idx [12] is not asserted |
| P3 robust first-order surrogate | **PASS** | err at 1e-2, decreasing, sampled_min check | err(1e-2)=4.463e-05, slope=2.00, sampled_ok=True |
| P3 adverse ~ sampled minimum | observed | gap rel to max(1,|lam0|) | eps sweep 1.0e-03, 3.0e-03, 1.0e-02 |
| P4 empirical Weyl lower-bound self-check | **PASS (empirical only)** | violations on same samples | L_emp = 2.108826e-04, L_nominal = 2.241500e-04, violations = 0, L_i mean = 1.247619e-04 |
| P5 projector movement / rank scan | observed | recorded tables | min adjacent eig gap = 6.624e-15 at tau* = 0.0500 (pair (1, 2)); rank event = False; sigma_min(B) < 1e-8 event = False |
| P5 crossing prediction degradation | observed | crossing pair first-order errors | max crossing rel err = 1.316e-06 vs P2 well-gapped 1.128e-07 (ratio 11.67x) |

## P1 raw rows (relative Frobenius error vs DK_fd)

| eps | rel_full | rel_frozen | ||DK_fd||_F |
|---:|---:|---:|---:|
| 1.000e-04 | 1.728e-07 | 4.442e-02 | 1.200e-04 |
| 2.371e-04 | 1.894e-07 | 4.442e-02 | 1.200e-04 |
| 5.623e-04 | 2.990e-07 | 4.442e-02 | 1.200e-04 |
| 1.000e-03 | 6.247e-07 | 4.442e-02 | 1.200e-04 |
| 1.334e-03 | 1.007e-06 | 4.442e-02 | 1.200e-04 |
| 3.162e-03 | 5.096e-06 | 4.442e-02 | 1.200e-04 |
| 7.499e-03 | 2.812e-05 | 4.442e-02 | 1.200e-04 |
| 1.778e-02 | 1.576e-04 | 4.441e-02 | 1.200e-04 |
| 4.217e-02 | 8.852e-04 | 4.438e-02 | 1.200e-04 |
| 1.000e-01 | 4.961e-03 | 4.438e-02 | 1.202e-04 |

The tolerance rows labelled `at eps=1e-3` are evaluated at an explicit
`on_grid=False` sample at eps=1e-3 because `logspace(-4,-1,9)` has no exact
1e-3 node; the 9 on-grid samples are used for slope fits (P2/P3 figures keep
the same convention).

Full clean window: eps in
[5.6e-04,
1.0e-01], slope
1.91 (r^2=0.998).
Frozen clean window slope -0.00.
||DK_full||_F = 1.199524e-04,
||DK_frozen||_F = 1.198848e-04,
||A0^T dW A0||_F = 5.328559e-06
(share 4.442e-02).

## P2 raw rows (chosen eigenvalues)

**chosen idx 23 (largest)**, lam=3.314405e-04, min gap=9.513e-05, full slope=2.00

| idx | eps | rel_full_eig | rel_frozen_eig | overlap+ |
|---:|---:|---:|---:|---:|
| 23 | 1.000e-04 | 8.033e-10 | 3.089e-07 | 1.00000000 |
| 23 | 2.371e-04 | 4.515e-09 | 7.352e-07 | 1.00000000 |
| 23 | 5.623e-04 | 2.539e-08 | 1.758e-06 | 0.99999999 |
| 23 | 1.000e-03 | 8.029e-08 | 3.162e-06 | 0.99999998 |
| 23 | 1.334e-03 | 1.428e-07 | 4.252e-06 | 0.99999996 |
| 23 | 3.162e-03 | 8.036e-07 | 1.055e-05 | 0.99999980 |
| 23 | 7.499e-03 | 4.528e-06 | 2.763e-05 | 0.99999886 |
| 23 | 1.778e-02 | 2.558e-05 | 8.038e-05 | 0.99999350 |
| 23 | 4.217e-02 | 1.455e-04 | 2.754e-04 | 0.99996197 |
| 23 | 1.000e-01 | 8.408e-04 | 1.149e-03 | 0.99976388 |

**chosen idx 12 (median)**, lam=6.164475e-08, min gap=5.936e-08, full slope=1.99

| idx | eps | rel_full_eig | rel_frozen_eig | overlap+ |
|---:|---:|---:|---:|---:|
| 12 | 1.000e-04 | 2.768e-09 | 1.308e-06 | 1.00000000 |
| 12 | 2.371e-04 | 1.558e-08 | 3.112e-06 | 0.99999998 |
| 12 | 5.623e-04 | 8.760e-08 | 7.429e-06 | 0.99999989 |
| 12 | 1.000e-03 | 2.770e-07 | 1.333e-05 | 0.99999966 |
| 12 | 1.334e-03 | 4.924e-07 | 1.790e-05 | 0.99999940 |
| 12 | 3.162e-03 | 2.765e-06 | 4.405e-05 | 0.99999660 |
| 12 | 7.499e-03 | 1.549e-05 | 1.134e-04 | 0.99998078 |
| 12 | 1.778e-02 | 8.643e-05 | 3.186e-04 | 0.99989085 |
| 12 | 4.217e-02 | 4.775e-04 | 1.028e-03 | 0.99937086 |
| 12 | 1.000e-01 | 2.602e-03 | 3.907e-03 | 0.99623189 |

**chosen idx 22 (well_gapped)**, lam=2.363068e-04, min gap=9.275e-05, full slope=2.00

| idx | eps | rel_full_eig | rel_frozen_eig | overlap+ |
|---:|---:|---:|---:|---:|
| 22 | 1.000e-04 | 1.130e-09 | 1.765e-07 | 1.00000000 |
| 22 | 2.371e-04 | 6.346e-09 | 4.149e-07 | 1.00000000 |
| 22 | 5.623e-04 | 3.567e-08 | 9.632e-07 | 0.99999998 |
| 22 | 1.000e-03 | 1.128e-07 | 1.664e-06 | 0.99999993 |
| 22 | 1.334e-03 | 2.005e-07 | 2.168e-06 | 0.99999988 |
| 22 | 3.162e-03 | 1.127e-06 | 4.491e-06 | 0.99999931 |
| 22 | 7.499e-03 | 6.326e-06 | 6.994e-06 | 0.99999609 |
| 22 | 1.778e-02 | 3.545e-05 | 3.866e-06 | 0.99997796 |
| 22 | 4.217e-02 | 1.977e-04 | 1.228e-04 | 0.99987494 |
| 22 | 1.000e-01 | 1.089e-03 | 9.118e-04 | 0.99927873 |

Max gap ratio among chosen eigenvalues (ascending convention):
1.603e+03.

Gate semantics: P2 pass reflects the well-gapped chosen eigenvalues only (tolerance rows at exact eps=1e-3): rel_full_eig < 1e-3 and rel_frozen_eig > 10x rel_full_eig.  Chosen eigenvalues whose min adjacent gap is < 1e-6 are flagged (gap_flags / simple_condition_flagged_indices) and the simple-eigenvalue condition is reported as invalid there, without being asserted

## P3 raw rows

Chosen eigenvalue idx 23 (largest),
lam0 = 3.31440505e-04, ||g_r||_2 = 2.618225e-04,
Ns = 150.

| eps | predicted worst | lam_adv | sampled_min | err_adv | err_sampled | sampled_min check |
|---:|---:|---:|---:|---:|---:|---|
| 1.000e-03 | 3.311787e-04 | 3.311788e-04 | 3.312846e-04 | 4.460e-07 | 3.195e-04 | True |
| 3.000e-03 | 3.306550e-04 | 3.306564e-04 | 3.309733e-04 | 4.015e-06 | 9.602e-04 | True |
| 1.000e-02 | 3.288223e-04 | 3.288371e-04 | 3.298891e-04 | 4.463e-05 | 3.219e-03 | True |

## P4 empirical Lipschitz / Weyl self-check

L_emp = 2.108826e-04, L_nominal = 2.241500e-04, violations = 0, L_i mean = 1.247619e-04.  L_emp and L_nominal are finite-sample / basis-direction diagnostics only: no uniform operator-Lipschitz constant over the full ball is derived, so no bound here is certified.  The Weyl check is self-consistent by construction on the same samples.

## P5 raw rows

### Projector movement (X0 + eps*dX_unit, dX seed 2718)

| eps | ||P(X+eps)-P(X0)||_2 | rank A shifted | sigma_min(B) shifted |
|---:|---:|---:|---:|
| 1.000e-03 | 5.238e-03 | 24 | 6.737e-04 |
| 1.000e-02 | 5.251e-02 | 24 | 6.628e-04 |
| 1.000e-01 | 5.115e-01 | 24 | 5.527e-04 |

### Near-crossing / rank-event scan

min adjacent eig gap = 6.624e-15 at tau* = 0.0500 (pair (1, 2)); rank event = False; sigma_min(B) < 1e-8 event = False.  Full scan table (41 tau values) is in the JSON.

the minimum adjacent gap along the scan occurs in the numerically-zero eigenvalue cluster (tau* pair magnitudes 8.01e-15 and 1.46e-14, well below sigma_min(A)~1e-8); sigma_min(B) stays ~6e-4..8e-4 and rank(B)=18 throughout

### First-order prediction at tau of minimum gap (crossing pair, eps=0.001)

| ascending idx | lam(tau*) | pred(+eps) | lam_obs(+eps) | rel err | overlap |
|---:|---:|---:|---:|---:|---:|
| 1 | 8.012115e-15 | 8.010955e-15 | 8.010954e-15 | 2.113e-07 | 1.000000 |
| 2 | 1.463572e-14 | 1.461050e-14 | 1.461052e-14 | 1.316e-06 | 1.000000 |

P2 well-gapped reference rel_full_eig at eps=1e-3: 1.128e-07.

## Cannot-establish section

* All claims are finite-dimensional and model-specific statements about the
  discrete N=16 whitened/realified smooth-basis Jacobians; no continuum-limit,
  exact-global-SE(2), estimator, or recovery claim is made.
* The P4 Lipschitz numbers are **explicitly empirical only**: L_emp comes from
  150 sampled (dX, eps) pairs and L_nominal from 18 basis directions.  No
  uniform operator-Lipschitz constant over the full ball is derived, and the
  Weyl lower-bound self-check is satisfied by construction on the same
  samples.  Nothing in P4 is certified.
* P1's DK_full is itself a centred finite-difference derivative (delta=1e-4),
  so its numerical floor (~O(delta^2) + roundoff/delta) explains the
  flattening at the smallest eps; slopes are fit only on the recorded clean
  window.
* P2 eigenvector matching assumes the chosen eigenvalue remains identifiable
  by maximum overlap; near-degeneracies or crossings can invalidate that
  label (gap flags are reported).
* P5 reports a near-crossing/rank-event scan along one seed-9191 random path;
  no theorem guarantees a crossing or rank event over that path, and the
  sigma_min(B)/rank events use one stated tolerance.  No uniform continuity
  or rank-event theorem is claimed.
* P5 projector movement uses the machine-rank-tolerance Range(A) basis; the
  observed movement mixes within-range rotation of near-degenerate directions
  and genuine subspace movement.

## Artifacts

- results: `results/family5_sensitivity.json`
- figures: family5_derivative_formula.png, family5_eigenvalue_prediction.png,
  family5_robust_surrogate.png, family5_rank_events.png
- this report: notes/family5_report.md

Self-cell formula used: self-cell v2 corrected 2026-09-03: complete equal-area disk integral (adds the -1/k_b^2 lower-endpoint term omitted by v1).

## Artifacts and digests

```text
ee3b4ac45ad296dba533a8ccefd107a6cc53e4357b7a768d07f7db90f3ee6585  src/helmholtz.py
762af3bdd33af1ab80563c7257b3805e7972f786fa88af1bf5e1b51a768526d7  src/family1_pilot.py
173d2b5d9f66be274c9a8422f4ac2c381bb91a60132b37c734748414cf3eb235  src/family2_algebraic_spine.py
ee7cac8b6233bcb0dbbb215df02452746dabf50774760910be5b16a0c0a40242  src/family4_frequency_trajectory.py
5e7368041cf69a3f23ee6bfcf498cc19f372c37bdf6fe84dfb1c434090758bb7  src/family5_sensitivity.py
782dd079312a5413925d757b70e43c4aa3d0a25680f664aea5ce4fe2bcd8bb6d  results/family5_sensitivity.json
a029e8518e5bc92900c690076bfa436578cd575abf21502ba634a0cd8a95a558  figures/family5_derivative_formula.png
730832c845b8161007afad4c556f4a9b4c0b8211ea424d7d34c83aeb56fac522  figures/family5_eigenvalue_prediction.png
0dd1ccacd9657679a61402f7acd0226bf3d604f5a2decee906754163a57caea0  figures/family5_robust_surrogate.png
86dac85c07ab01eff84d1d631376d41706c7cea7694695dc6b64e4672536e79b  figures/family5_rank_events.png
```
