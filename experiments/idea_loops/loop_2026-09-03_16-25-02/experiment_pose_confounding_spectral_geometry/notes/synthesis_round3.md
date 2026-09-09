# Round-3 synthesis (read-only)

Date: 2026-09-04 (SGT). Experiment:
`experiment_pose_confounding_spectral_geometry`.

**Read-only note:** this synthesis was produced without rerunning any
experiment and without modifying any source, result, figure, or manuscript
file. Every number below was copied from the cited family reports
(`notes/family10_protocol_correction.md`,
`notes/family10_online_slam_toy_n120.md`,
`notes/family11_snr_diversity.md`, `notes/family12_trajectory_replication.md`,
`notes/family13_rank_autopsy.md`, `notes/family14_robustness_scope.md`,
`notes/family15_born_control.md`, and `notes/synthesis_round2.md`) and from
the corresponding result JSONs, checked only where an exact value or an
aggregate count was needed.

## Round-3 families, roles, and artifacts

| family | role | report path | result JSON path | figure paths | key status |
| --- | --- | --- | --- | --- | --- |
| Family 10 | protocol-corrected online NLS covariance toy (n120; fixed true pose `dx=0`; corrected `P_samp` sandwich target) | `notes/family10_online_slam_toy_n120.md`; `notes/family10_protocol_correction.md` | `results/family10_online_slam_toy_n120.json` (empirics; supporting `results/family10_affine_control.json`) | `figures/family10_cov_eigenvalues_n120.png`; `figures/family10_confounded_direction_ratios_n120.png`; `figures/family10_born_vs_fullwave_deviation_n120.png` | corrected target lowers predicted ratio, but free-pose covariance mismatch persists for Born and full-wave |
| Family 11 | SNR/frequency-set replication across scenes (non-monotonicity, plateau, equal-SNR sets, duplicate invariance, bandwidth) | `notes/family11_snr_diversity.md` | `results/family11_snr_diversity.json` | `figures/family11_frequency_sets.png`; `figures/family11_snr_sweep_scenes.png` | two_blob and ring replicate both f2 values; low_contrast f2=1.8 replicates 2/3 non-monotone |
| Family 12 | trajectory replication and mechanism (finite N=16) | `notes/family12_trajectory_replication.md` | `results/family12_trajectory_replication.json` | `figures/family12_retained_rankings.png`; `figures/family12_principal_angles.png`; `figures/family12_mechanism.png` | trajectory-ordering hypothesis false 0/24; c2 rank identity holds 93/96 with 3 near-null threshold exceptions |
| Family 13 | rank-boundary autopsy of the ring-scene c2 residual (tolerance classification) | `notes/family13_rank_autopsy.md` | `results/family13_rank_autopsy.json` | `figures/family13_rank_autopsy.png` | ring residual 1 is a numerical rank classification of a near-null singular value, not an algebraic failure |
| Family 14 | robustness-scope audit (explicit demotion) | `notes/family14_robustness_scope.md` | `results/family14_robustness_scope.json` | none (report produces no figure) | full nonlinear robustness claims demoted; only the affine tangent certificate is retained |
| Family 15 | Born analytic / literature control | `notes/family15_born_control.md` | `results/family15_born_control.json` | `figures/family15_born_control.png` | Born identity exact; observed slopes ~1; no quantitative external Born benchmark |

## Key numbers (copied from the cited reports)

### Family 10 — protocol-corrected online NLS toy

- Protocol: the n120 Monte Carlo draws only data noise with the hidden pose
  fixed at `dx = 0` on every trial, so the correct linearized frequentist
  target of the map component of the penalized free-pose estimator is the
  sandwich `P_samp = [H^-1 H_data H^-1]_{1:p,1:p}` with `H_data =
  [[K_IS, A^T B], [(A^T B)^T, B^T B]]` and `H = H_data + diag(0_p, alpha I_q)`
  (`alpha = 1.0`), not the Bayesian marginal `P_free = K_eff^-1`.
- Known-pose convergence: Born known/free 120/120 and 115/120; full-wave
  known/free 120/120 and 119/120. Known-pose rel Fro vs `P_known`: Born
  `0.0779`, full-wave `0.4745` (rel spectral `0.0783` / `0.4535`).
- Stored free-pose deviation against the old target: rel Fro `13.49` (Born)
  and `4.04` (full-wave) vs `P_free = K_eff^-1`.
- Trace-ratio comparison: empirical free-pose ratio `82.92` (Born) and
  `11.019` (full-wave); old `P_free` target `6.068` / `3.629`; corrected
  `P_samp` target `3.889` / `2.755` (Born affine control
  `trace_ratio_exact_fixed_pose = 3.889`). The correction lowers the target by
  ~36% (Born) and ~24% (full-wave) but does not remove the mismatch:
  recomputed rel Fro of `Cov_free` vs `P_samp` is `22.0` (Born) and `5.40`
  (full-wave).
- Residual nonlinearity after correction: Born empirical free covariance is
  inflated ~21.3x in the ratio-to-ratio sense or ~22.5x in
  `tr(Cov_free)/tr(P_samp)`; full-wave ~4.0x or ~5.8x.
- Directional FD self-check of the full-wave stacked residual Jacobian: batch
  rel Fro `4.6595e-10`; max directional column error `2.63586e-09` at
  `eps = 1e-6`.
- Caveat copied from the reports: Born free covariance is heavily
  tail-dominated (free-fit `c_error_l2` up to `7.4e3` vs `4.9e2` max for
  known-pose fits), so the factors carry wide Monte Carlo uncertainty even
  though empirical >> `P_samp` is robust.

### Family 11 — SNR / frequency-set replication

- Part A per-(scene, f2) replication, all directions on the 13-point primary
  grid `logspace(-2, 4, 13)`: two_blob 3/3 non-monotone for both f2 values;
  ring 3/3 for both f2 values; low_contrast f2=1.4 is 3/3 but low_contrast
  f2=1.8 is 2/3 non-monotone (`pattern_replicates_all_3 = false`). The
  "replicates both f2 values" flag is True for two_blob and ring and False for
  low_contrast.
- Example peak movements (u0/u1/u2): two_blob f2=1.4 `0.1257` / `0.2335` /
  `0.5756`; two_blob f2=1.8 `0.4343` / `0.1612` / `0.4596`; low_contrast
  f2=1.8 `0.2500` / `0.4881` / `0.3999`. Finite positive high-SNR tail
  plateau medians are recorded (e.g., two_blob f2=1.4 `5.011e-3` / `4.315e-3`
  / `1.261e-2`).
- Part B equal-per-frequency-SNR retained DOF (sum rho) grows with frequency
  count on all scenes: two_blob F1 `11.2929456138` -> F5 `15.9117483745`;
  ring F1 `15.8438774096` -> F5 `20.0568946868`; low_contrast F1
  `16.9803786027` -> F5 `18.9645532086`.
- Duplicate-block (DUP) invariance PASS under the joint-scaled prior
  `alpha = 2` on every scene (all gates `< 1e-10`); max |movement| residuals
  `1.256e-16` (two_blob), `2.682e-16` (ring), `2.187e-16` (low_contrast).
- Fixed-count bandwidth sweep (always three frequencies): BW_mid
  (1-1.4-1.8) has the highest retained DOF on each scene
  (two_blob `15.6364688685`; ring `19.2602533997`; low_contrast
  `18.6183940010`), with BW_narrow and BW_wide lower on every scene.

### Family 12 — trajectory replication (0/24 and 93/96)

- Scope: 4 scenes x controls A/B x 4 trajectories x F1/F2/F3, giving 96
  rank/metric appendix cells and 24 (scene, control, frequency-set) ordering
  rows. Hypothesis ordering compared is
  `['circle360', 'arc180', 'arc90', 'straight']` by retained DOF descending.
- Trajectory-ordering hypothesis holds in **0/24** cells: both the control-A
  and control-B "equals old hypothesis" flags are False in every one of the 24
  rows (no forced pass). `rank(B) = 18` in all 96 cells.
- c2 rank identity (`r_KIS - r_KSL = r_A + r_B - r_AB`) holds in **93/96**
  cells. The 3 exceptions are single-frequency cells whose directly formed
  `K_SLAM` singular values fall below the family-2 machine rank tolerance
  (relative sizes ~1e-13 to 1e-17 of `sigma_1(K_SLAM)`): ring/arc90/F1 control
  A (lhs 3, rhs 0), ring/arc90/F1 control B (lhs 1, rhs 0), and
  offset_edge/arc90/F1 control A (lhs 1, rhs 0). These are recorded as
  numerical rank-threshold observations, not claims of algebraic identity
  failure.
- Cross-check 1 (control B, F1, two_blob vs family4c report) PASS within
  `1e-5`: arc90 `9.408661034`, straight `8.269820239`, arc180 `7.997790261`,
  circle360 `7.701222702` (abs diffs `2.4e-7`..`3.4e-7`).
- Cross-check 2 (control A, F1, two_blob): computed ordering `arc90 >
  straight > circle360 > arc180` matches the `family4_results.json` artifact
  ordering but fails the request's stated ordering (which transposes
  arc180/circle360); recorded without a forced pass.
- Descriptive account only: scene, frequency set, control geometry, and
  normalisation each change the ranking in at least one cell; no universal
  ordering or certified mechanism is claimed.

### Family 13 — tolerance-classification verdict and backward residuals

- Ring verdict (`b_tolerance_classification_of_a_near_null_singular_value`):
  `C = (I - Z Z^T) A_s` has full column rank at every predeclared relative
  tolerance (`sigma_24(C)/sigma_1(C) = 4.74995e-08`), while the directly
  formed `K_SLAM` has a near-null singular value
  `sigma_24/sigma_1 = 2.25753e-15` that drops below tolerance on part of the
  grid. This is a numerical rank classification, not a discontinuous exact
  rank change and not unstable subtraction.
- Backward residuals (ring / two-blob): rel Fro
  `||K_SLAM - C^T C|| / max(||K_SLAM||_F, ||C||_F^2)` = `5.00069e-16` /
  `3.17762e-14`; rel Fro R_op-definition residual = `1.51994e-07` /
  `8.26455e-09`; projector idempotency `||P^2-P||_F / ||P||_F` ~ `3.6e-16`;
  `||Z^T Z - I||_F` ~ `1.9e-15` / `1.7e-15`.
- Two-blob verdict (`none_family2_machine_residual_zero`): c2 residual 0;
  `sigma_24/sigma_1(K_SLAM) = 1.52533e-13` sits above the family-2 relative
  threshold `24*eps ~ 5.3e-15`, so no rank boundary occurs under the declared
  machine rule.
- Tolerance-grid behavior is recorded for both scenes (ring first shows c2
  residual 1 at tol_rel 1e-14..1e-13; coarse tolerances near/above 1e-12 shed
  near-null modes on both scenes); no theorem is claimed.
- Engineered algebraic falsifier (two_blob, f=1.0; `rank(B(t))` 18 -> 17 at
  t=0): retained mass jumps `9.40716187` (t>0) to `10.40528058` (t=0) with
  no-prior projector operator-norm jump `1` and rel Fro information jump
  `0.139957`. This is an algebraic control, not a physical ring event.

### Family 14 — robustness demotion with affine certificate and sampled margins

- Recomputation reproduces every stored affine certificate exactly at IEEE
  double precision (max abs relative diff `0.000e+00`):
  `L_cert_affine = 1.130863e-02` (eps 1e-3), `1.143821e-02` (eps 3e-3),
  `1.189479e-02` (eps 1e-2); pointwise structural constant
  `L_struct_tight(X0) = 1.124399e-02`.
- FD validation supports the implemented Jacobian only (conditional layer):
  J_A/J_B rel Fro `5.990607e-07` / `5.374392e-07` (~6e-7); max per-column rel
  diffs `7.671831e-07` / `7.613361e-07`; consecutive-difference ratios ~4.0
  support O(h^2) for the centered-FD estimates.
- Sampled margins (executed evidence layer): full-nonlinear sampling was 500
  unit directions per eps (seed 9191) and affine sampling 2000 per eps (seed
  7171), with zero sampled violations. At eps=1e-2 the full-nonlinear worst
  sampled quotient is `3.228342e-04` (margin `1.157196e-02`, ratio `36.8449`)
  versus the affine certificate `1.189479e-02`; affine worst quotient
  `3.744008e-04` (affine margin `1.152039e-02`, affine ratio `31.7702`).
- Explicit demotion (verbatim from the report):
  > All full nonlinear execution-error robustness claims are DEMOTED to
  > empirical/structural observations; only the affine tangent certificate is
  > retained as a finite-dimensional conditional certificate.
- No claims: no interval proof, no Hankel bound, no continuum robustness, and
  no "certified nonlinear robustness". The generalized-derivative formulas
  require locally constant rank and are inapplicable at rank events (engineered
  falsifier in Family 13/5); repeated eigenvalues require compressed
  derivatives with only cluster-control evidence.

### Family 15 — Born slopes/errors and literature statement

- Born composition identity is bitwise-exact at every scale (rel Fro `0` vs
  the explicit `G_S diag(E_inc)` stack); the independent `A(chi=0) == A_born`
  check is `7e-17` rel Fro (Family 3b check D).
- No-prior relative deviations grow with scale: rel A `0.000499888` (s=0.001)
  -> `0.504938` (s=1); rel K_IS `0.000701769` -> `0.652688`; rel K_eff
  `0.000701769` -> `0.651679`. Retained DOF full vs Born: `9.15161` vs
  `9.15117` (s=0.001) and `9.40716` vs `9.21936` (s=1).
- Observed log-log slopes (full 6-point grid, reported not forced): A operator
  `1.00261` (r^2 `0.999994`), K_IS `0.994365` (r^2 `0.999841`), K_eff
  `0.994199` (r^2 `0.999837`) - near 1, as expected if the relative Born error
  is dominated by the O(s) first correction.
- Finite-prior (alpha=1.0) K_eff retained DOF stays 24 at low s and is
  `23.9602` (full and Born) at s=1; most-confounded-subspace diagonal squared
  overlap of u1 drops from `0.999965` (s=0.001) to `0.0271854` (s=1), with
  best-pairing sum `2.99997` -> `2.12697`.
- Literature field (verbatim citation): M. L. Diong, A. Roueff, P. Lasaygues,
  A. Litman, "Impact of the Born approximation on the estimation error in 2D
  inverse scattering", Inverse Problems, 2016. Family 15 is a
  formula-level/internal analytic control with the same conceptual structure.
- No-benchmark statement (verbatim): "No quantitative external benchmark was
  possible because the published setup's geometry, normalization, noise model,
  parameterization, and units were not reproduced here; this is not
  independent validation."

## Writeup

`manuscript/main_round3.tex` (120,321 bytes) and its compiled PDF
`manuscript/main_round3.pdf` (5.2 MiB) were produced in this workspace on
2026-09-04 00:46, consistent with the accompanying `manuscript/main_round3.log`
from the same build. The draft keeps every claim aligned to the result JSONs
and imposes the four-layer scope separation: (1) exact finite-dimensional
linear algebra (rank/retention/principal-angle/kernel identities plus the
affine tangent certificate); (2) conditional/discretized full-wave
differentiation (FD-validated Jacobians certify only the implemented discrete
model); (3) executed numerical evidence and preserved counterexamples
(labelled OBS/CTX, never relabelled); and (4) open continuum, global, real
system, and novelty questions. The explicit demotion sentence from Family 14
("All full nonlinear execution-error robustness claims are DEMOTED to
empirical/structural observations; only the affine tangent certificate is
retained as a finite-dimensional conditional certificate") is recorded in the
manuscript, and the family-12 result table reports the trajectory hypothesis
as false in 0/24 cells and the c2 rank identity as 93/96 cells with the three
near-null threshold exceptions.

## Failures and uncertain claims

- Trajectory hypothesis false 0/24.
- c2 rank identity 93/96 with 3 near-null threshold exceptions.
- low_contrast f2=1.8 non-monotone 2/3.
- Online-toy free-pose covariance mismatch.
- Full nonlinear robustness demoted.
- No quantitative external Born benchmark.

## Next round should do

- (i) Only if desired, attempt a genuine validated interval/rigorous bound for
  the full nonlinear operator with a validated arithmetic library (otherwise
  keep the demotion).
- (ii) Pilot 3D scalar or vector Helmholtz extension while keeping G_S/
  current-space vs map-tangent A/K_eff explicit.
- (iii) Larger-N online-toy covariance study or alternative sampler to test
  the free-pose covariance mismatch.
- (iv) Quantitative Born benchmark only if a published setup can be reproduced
  exactly.
- (v) Keep the manuscript aligned to result JSONs.

## Confirmation

Synthesis written to `notes/synthesis_round3.md` (read-only; no experiment
rerun). Failure list: trajectory hypothesis false 0/24; c2 rank identity 93/96
with 3 near-null threshold exceptions; low_contrast f2=1.8 non-monotone 2/3;
online-toy free-pose covariance mismatch; full nonlinear robustness demoted;
no quantitative external Born benchmark.
