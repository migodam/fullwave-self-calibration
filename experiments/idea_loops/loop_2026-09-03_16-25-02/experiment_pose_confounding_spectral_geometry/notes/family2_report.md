# Family 2: algebraic spine - claim-status report

Date: 2026-09-03 (UTC run 09:46:36; SGT ~17:46).
Experiment: `experiment_pose_confounding_spectral_geometry`.
Family 1 (forward model + analytic Jacobians) is accepted; this run consumes
only its validated `helmholtz.build_AB` / `whiten_realify` output.

## Scope and method

`results/family2_results.json` was produced by the exact command (from the
experiment root):

```bash
.venv/bin/python src/family2_algebraic_spine.py
```

Wall runtime 4.25 s (recorded in the JSON).  Environment (JSON `platform`):
Apple Silicon CPU (`arm64`), Python 3.12.13, numpy 2.5.2, scipy 1.18.1,
matplotlib 3.11.1; no GPU/MPS/CUDA; no RNG (deterministic dense algebra).

The whitened/realified data space has m = 2*T*n_rx = 48 real rows, pose block
q = 3*T = 18.  Two map parameterizations are analysed:

- pixel basis: A_R = whiten_realify(A_pix, B, None)[0], n = N^2 = 256,
  r_A = 48 = m (Range(A) is the full data space);
- smooth basis: A_smooth_R = A_pix_R @ S, S (256 x 24) unit-2-norm Gaussian
  RBF columns on a 4 x 6 centre lattice of D = [-0.5, 0.5]^2,
  sigma_b = 0.16; n = p = 24, r_A = 24.

B_R is shared by both cases.  Rank tolerance throughout:
`tol(M) = max(M.shape)*eps_machine*sigma_1(M)`.  K_eff is
`K_IS - A^T B (B^T B + alpha I)^{-1} B^T A` with the q x q system solved by
`np.linalg.solve` (or `np.linalg.pinv` for the singular prior).

Retention operators are evaluated as `Q_A^T W Q_A` (W = P_perp for no prior,
W = I - B(B^T B + alpha I)^{-1} B^T for a finite prior), which is algebraically
identical to `diag(1/s_A) V_A^T K V_A diag(1/s_A)` because
`A V_A diag(1/s_A) = Q_A`.  The factored form avoids amplification by
`1/s_A^2` in near-null map directions (`sigma_min(A)` is ~1e-8 in both
cases).  `rank(L_X)` (c4/c5) is counted from
`X = (B^T B + alpha I)^{-1/2} B^T A`, with `rank(L_X) = rank(X)` in exact
arithmetic; the direct difference `K_IS - K_eff(alpha)` is also recorded and
is inflated by a ~`eps ||K_IS||` cancellation floor for large alpha.

Self-cell marker (corrected harness): the runs stamp
`hh.SELF_CELL_FORMULA` and version:
`I_self = (i*pi*a/(2*k_b))*H_1^(1)(k_b*a) - 1/k_b^2, a = h/sqrt(pi)`
(self-cell v2, 2026-09-03).

## Tolerances used

| gate | value |
|---|---:|
| rank tolerance | `max(shape)*eps*sigma_1(M)` |
| c1 pass | subspace distance < `1e-8*max(1, ||A||_2^2)` (backward-scaled) |
| c3 predicted-spectrum pass | max |rho - rho_pred| < `1e-10*max(1,||A||_2^2)` |
| c3 bound check | `-1e-10 <= rho <= 1 + 1e-10` |
| c5 PSD/interlacing | violations above `-1e-10*max(1,||K_IS||_2)` reject |
| c8 invariance | relative Frobenius / spectrum differences < `1e-10` |

`max(1, ||A||_2^2) = 1` for both bases here, so the c1/c3 numerical gates are
1e-8 / 1e-10 absolute.

## Claim-status table

| claim | status | executed check | key numbers (pixel / smooth) |
|---|---|---|---|
| null(K_SLAM) = null(P_perp A) | **PASS (smooth); PASS with numerical caveat (pixel)** | c1 subspace distance, null residuals | smooth dist = 0 (empty nulls, rA=24=rank KSL); pixel dist = 8.47e-7 > 1e-8 gate, but min singular(N1^T N2) = 1 - 4e-13, max ||K_SLAM N2|| = 2e-19, max ||C N1|| = 2e-13, nullities 226 = 226 |
| rank identity r_KIS - r_KS = r_A + r_B - r_AB | **PASS** | c2 exact integer residual | pixel 48-30 = 18 = 48+18-48, residual 0; smooth 24-24 = 0 = 24+18-42, residual 0; d_inter >= max(0, r_A+r_B-m) holds |
| retention spectrum = 1 - cos^2 (padded) | **PASS** | c3 | max |rho-rho_pred| = 2.89e-15 / 1.55e-15; count<1 = 18 = r_B both; rho in [-1e-15, 1+3e-15] pixel, [5.5e-11, 1+2e-15] smooth |
| at most r_B retention values differ from 1 | **PASS** | c3 count | 18 = r_B both |
| rank(L_X) <= r_B, alpha in {1e-6..1e2} | **PASS** | c4 stable factor rank | max factor rank = 18 = r_B both; direct-difference rank is inflated at alpha=100 (pixel 129, smooth 23) by cancellation floor and is documented, not used for the gate |
| K_SLAM <= K_eff(alpha) <= K_IS (PSD) | **PASS** | c5 min eigenvalues | all >= -2.5e-19 (pixel), >= -2.9e-20 (smooth); gate -1e-10 |
| interlacing with p = rank(L_X) | **PASS** | c5 | max positive violation 5.0e-21 (pixel), 6.7e-21 (smooth) <= 1e-10 |
| alpha->0 path approaches K_SLAM | **PASS (as limiting behaviour; not reached at grid's left edge)** | c6 | min d0 on grid 0.209 (pixel) / 0.137 (smooth) at alpha=1e-6; supplementary sweeps show monotone d0 -> 0.0056/0.0036 at alpha=1e-8 |
| alpha->inf regular path approaches K_IS | **PASS** | c6 | dinf(1e4) = 5.17e-7 / 5.50e-7; min over grid at right edge |
| finite prior gives weighted shrinkage in (0,1] | **PASS** | c6 | rho_X in (0.086, 1+3e-15] pixel and (0.0868, 1+1e-15] smooth at alpha=1e-3; max sorted-spectrum change vs no-prior ~0.9995 / 0.9889 (nonzero) |
| singular-support prior never converges to K_IS as alpha->inf | **PASS** | c7 | dinf_sing(1e4) = 4.72e-4 / 2.30e-4, 418-912x the regular-prior dinf at the same alpha. The illustrative `> 0.1` magnitude is not met for the unpenalised pose column q-1 (last pose, theta) in this geometry; the non-convergence is still clear |
| duplicate no prior: K_IS, K_SLAM scale by 1+c^2, rho unchanged | **PASS** | c8(a) | rel K_IS 1.11e-15 / 1.08e-15; rel K_SLAM 2.53e-13 / 2.90e-13; max |rho2-rho| 6.66e-16 / 1.19e-12, all < 1e-10 |
| duplicate fixed finite prior changes rho (data-to-prior weighting) | **PASS** | c8(b) | max |rho_X2-rho_X| = 3.99e-2 / 3.98e-2 at alpha=1; monotonically decreasing over alpha in {1e-3..1e2} |
| duplicate with jointly scaled prior reproduces single system | **PASS** | c8(c) | rel K_eff 1.00e-15 / 1.04e-15; max rho diff 1.67e-15 / 1.22e-15 |
| smooth-basis c1-c3 persist under N=32, N=40 | **PASS** | resolution | c1 distances 0, c2 residuals 0, c3 pred residuals 1.90e-15 / 8.24e-15; rB=18; build+check 1.06 s / 3.04 s |

## Notes on the only gate that is not clean

c1 (pixel basis): the nullspace of K_SLAM is 226-dimensional, but the smallest
positive singular value of K_SLAM is ~1.2e-14 (`sigma_min(P_perp A)^2`, with
`sigma_min(A) ~ 9.7e-9`).  Any orthonormal basis of the zero cluster returned
by `eigh` (or even by SVD of K_SLAM) carries basis mixing at the
`~eps*||K_SLAM||/gap` level, observed here as a subspace distance 8.47e-7
against the full-SVD nullspace of C = P_perp A.  The identity is nevertheless
certified by dimension equality (both nullities 226), by
max ||K_SLAM u|| over N2 = 1.96e-19 (N2 is in the nullspace of K_SLAM), by
min singular(N1^T N2) = 1 - 4e-13, and by rank(K_SLAM) = rank(C) = 30.
`pass_gate` in the JSON is therefore `false` for pixel while
`identity_support` is `true`; the JSON explains the numerical reason.

c6: on the specified grid `logspace(-6,4,25)` the smallest distance to
K_SLAM is at alpha = 1e-6 but is not small (0.209 pixel, 0.137 smooth):
sigma_min(B_R) = 6.75e-4 has sigma^2 = 4.6e-7 < alpha = 1e-6, so the weakest
pose direction is still only weakly removed at the left edge of the grid.
The sweep is monotone and a supplementary check gives d0 = 5.6e-3 (pixel) /
3.6e-3 (smooth) at alpha = 1e-8, consistent with K_eff -> K_SLAM as alpha -> 0.

c7: with the prescribed P[q-1,q-1] = 0 (unpenalised coordinate 18, pose 6
theta), dinf_sing(1e4) = 4.72e-4 / 2.30e-4.  It is clearly nonzero and
418-912x the regular-prior dinf(1e4), so the non-convergence claim is
supported; the prompt's illustrative "e.g. > 0.1" magnitude is not met,
because that residual-gauge direction has modest sensitivity in this
geometry.  This is recorded (`example_magnitude_gate_gt_0_1 = false`).

c8 smooth: max |rho2-rho| = 1.19e-12 passes the 1e-10 gate but is larger than
the pixel value (6.66e-16).  This is consistent with the ~1e-10 bound
`eps*||A_dup||_2/sigma_min(A_smooth)` for singular-subspace accuracy in the
96-row duplicate system; the K-matrix scaling checks are at 1e-15.

## Artifacts and digests

```text
173d2b5d9f66be274c9a8422f4ac2c381bb91a60132b37c734748414cf3eb235  src/family2_algebraic_spine.py
330d0a24a63a8f68b3d197818dbf5f0d01374bd17f4d46b868d5cf6dbf2f5862  results/family2_results.json
ca3d6458f043eab7df9954d5d6a24935950e8323fc080aece4efb0a2cc68385b  figures/family2_retention_spectra.png
34112ae2462ba969b6b8c17c2aa876f55dadf6acff34c0ba85caca2bdc6e1244  figures/family2_prior_sweep.png
a68356562a93de756462dac7e544ca5a2d35b85b26ef07ccbece03a001cffb43  figures/family2_duplicate_prior_effect.png
ee3b4ac45ad296dba533a8ccefd107a6cc53e4357b7a768d07f7db90f3ee6585  src/helmholtz.py
762af3bdd33af1ab80563c7257b3805e7972f786fa88af1bf5e1b51a768526d7  src/family1_pilot.py
```

## Cannot-establish section

- These are finite-dimensional, dense, double-precision linear-algebra checks
  on the single discrete (A_R, B_R) pair (and its duplicate).  They do not
  establish continuum-limit behaviour, convergence under grid refinement of
  the retention spectrum itself, or any claim about noise-whitening other than
  the identity-W (realified) case.
- c1 pixel-basis subspace distance does not meet its absolute 1e-8 gate with
  the prescribed eigh-based null basis because of the ill-conditioned zero
  cluster of K_SLAM (smallest nonzero singular value ~1e-14); the algebraic
  identity is certified by the residual/dimension evidence above, not by the
  eigh-basis distance alone.
- The prior-sweep statement "d0 small at alpha=1e-6" cannot be established on
  this geometry/grid (see c6 note); only the monotone limiting behaviour is
  established.  The singular-prior `> 0.1` illustrative magnitude is not met
  for the prescribed unpenalised coordinate.
- No inverse-problem, estimator, bias, or recovery claim is made anywhere in
  this report; spectra/ranks describe the linearised information geometry
  only.
