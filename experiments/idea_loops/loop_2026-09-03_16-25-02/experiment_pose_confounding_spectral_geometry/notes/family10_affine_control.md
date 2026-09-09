# Family 10 affine control (Born stacks)

Date: 2026-09-03T14:37:29.569910+00:00 (UTC)

## Scope

finite-dimensional affine/linearized control only; no continuum, global, or production claim

## Protocol

- Reused family10 scene and Born stack builders by import (no main() side effects).
- A_stack/B_stack/y0 shapes: [144, 24] / [144, 18] / [144].
- Linear-model consistency ||y0 - A c0||/||y0|| = 8.478e-15.
- Primary affine MC: 500 trials, y = y0 + noise, noise ~ N(0,I), seed=20260903;
  known-pose min_c ||A c - y||^2 and free-pose min_{c,dx} ||[A B][c;dx]-y||^2 + alpha||dx||^2
  solved by lstsq. Covariances are unbiased sample covariances of c_hat - c0.
- alpha = 1.0; generalized directions from family10 (K_eff v = rho K_IS v).

## Primary results (requested fixed-true-pose protocol)

Relative Frobenius deviations ||Cov_emp - P_pred||_F/||P_pred||_F:

- known-pose vs P_known = K_IS^-1: 0.085050
- free-pose vs P_free = K_eff^-1: 0.392003
- free-pose vs exact fixed-pose sampling covariance P_samp: 0.033441

Trace ratio tr(Cov_free)/tr(Cov_known):

- empirical: 3.616413
- predicted P_free/P_known: 6.068324
- exact fixed-pose P_samp/P_known: 3.888801

Three most-confounded generalized directions (predicted exact ratio vs empirical ratio;
the affine exact fixed-pose prediction is shown for comparison):

| dir | rho | predicted (K_eff^-1) | empirical | exact fixed-pose |
|---|---|---|---|---|
| 0 | 7.645052e-04 | 19.176164 | 10.052782 | 11.104848 |
| 1 | 3.892066e-03 | 7.119186 | 3.991712 | 4.416986 |
| 2 | 2.071578e-02 | 8.690621 | 4.293703 | 4.674657 |

## Diagnostic: prior-consistent replication

With a per-trial hidden pose dx_i ~ N(0, alpha^-1 I) and data y_i = y0 + B dx_i + noise_i,
the same affine free estimator is empirically consistent with K_eff^-1:

- free-pose vs P_free rel Frobenius: 0.066088
- trace ratio empirical: 5.446708 vs predicted 6.068324

## Honest interpretation

1. The affine control does **not** match P_free = K_eff^-1 under the prompt-specified protocol
   (y = y0 + noise only, dx always 0 in the data): the exact normal-equation estimator has
   repeated-sampling covariance P_samp = top-left block of H^-1 H_data H^-1, which is smaller.
2. This is an exact linear-algebra statement, not sampling noise and not nonlinearity:
   the prior rows carry zero residual noise in this protocol, so Cov = M J^T diag(I,0) J M
   rather than M.  The 500-trial empirical free covariance matches P_samp to sampling error.
3. P_free = K_eff^-1 is realized by the prior-consistent replication above (hidden dx drawn per
   trial), and the code's W/K_eff/P_free algebra (Schur identity) is verified to roundoff.
4. Therefore family10's reported free-pose mismatch is not explained solely by model
   nonlinearity under the implemented Monte Carlo protocol: even the exact affine model
   would fail the requested P_free comparison by roughly x1.5 (trace) and x1.7-1.9 in the
   most-confounded directions.  Nonlinearity adds a further large inflation on top.
5. To isolate nonlinearity alone against K_eff^-1, the Monte Carlo should draw dx_i per trial
   (prior-consistent generative model).  Alternatively, compare the fixed-dx NLS Monte Carlo
   to P_samp.  Which target is intended is a parent-level decision.

## Artifacts

- Results: `results/family10_affine_control.json`
- Notes: `notes/family10_affine_control.md`
- Figure: `figures/family10_affine_control.png`

## Artifacts and digests

```text
25d2418d133547e93702b482c8363a9dfd490246e0b3ed7893607bee61ba93eb  src/family10_affine_control.py
44a976104b0fa329a5699ce8c048e2fa3ed125115465bd6807ad117fb303e45b  src/family10_online_slam_toy.py
ee3b4ac45ad296dba533a8ccefd107a6cc53e4357b7a768d07f7db90f3ee6585  src/helmholtz.py
762af3bdd33af1ab80563c7257b3805e7972f786fa88af1bf5e1b51a768526d7  src/family1_pilot.py
173d2b5d9f66be274c9a8422f4ac2c381bb91a60132b37c734748414cf3eb235  src/family2_algebraic_spine.py
b36ea33076ca7c5c7b00a760e86113c772a57f95eed163e447a6faf2488b5643  results/family10_affine_control.json
56db948be417a7543680dc63093dc4bc07e8bde4691af6f9adbd67254032ac7f  figures/family10_affine_control.png
```
