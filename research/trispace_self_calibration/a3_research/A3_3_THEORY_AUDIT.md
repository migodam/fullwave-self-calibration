# A3_3 integration and mathematical corrections

2026-09-08. Parent read Q3_3.md and all 1256 lines of A3_3.md. These are
research inputs, not independently verified proofs or publication evidence.
Previous experiments and sampled-error baseline registration remain unchanged.

## What is accepted, and what is not

The useful organizing question is physical attribution: which phase changes
are due to geometry, electronics, material, or forward-model discrepancy?
This strengthens the existing direction; renaming it a spectral factor graph
does not establish novelty. SOM remains a replaceable numerical engine.

| Input proposal | Parent decision and required qualification |
| --- | --- |
| Complex gains can absorb geometry phase | Accepted conditionally; existing exact translation/clock and far-field results specify the geometry and gain sharing. Effective calibration does not identify electronic truth. |
| Incidence complement equals cycle space | For a connected bipartite sum-incidence matrix, a node sign change gives oriented incidence; rank is V-1. With whitening/field factors the relevant projector is the transformed tangent projector, not the unweighted graph projector. Nonzero signals and gain sharing matter. |
| Task-aware generalized spectrum | Useful candidate diagnostic, not yet a performance guarantee. Realify, whiten, scale parameters, profile other nuisances, and use consistent coordinates before comparing spectra. |
| Orthogonal discrepancy can be safely suppressed | Exact statement only for orthogonality to the **profiled** task tangent in the specified metric and a frozen local linear problem. Orthogonality to raw pose columns is insufficient. |
| Aligned discrepancy forces model upgrade | Exact indistinguishability can be proved for unrestricted aligned error. Finite covariance downweighting can still improve risk; near alignment alone does not prove impossibility. A new reference must actually distinguish the ambiguous direction. |
| Low-rank covariance is cheap | The small matrix inversion is cheap; obtaining trustworthy modes may dominate. Our earlier multifidelity timing failures remain evidence against assuming a total saving. |
| Risk equals inverse-information trace plus bias | Correct for ordinary LS with unit noise, or correctly specified GLS. For heuristic weights and unchanged white noise, the variance is a sandwich, not generally inverse weighted information. |
| Spectral policy is certificate-driven | Not yet. Adjacent meshes and sampled covariance do not enclose true error. Until a valid enclosure/probability model is supplied, call this a diagnostic or surrogate-risk policy. |
| MFSOM 2025 detailed equations and numerical claims | Primary v2 Sections II–IV-A independently rechecked this pass; complex frequency/transmitter factors and reported imaging improvement confirmed. Physical-source attribution remains our analysis, not their demonstrated result. See A3_3_PRIOR_LEDGER.md. |
| Plumlee / generalized contrastive PCA / LIS / VarPro parallels | Prior-work leads, not originality evidence. Generic projection, marginalization and spectral ratios are not claimed new. Verify exact scope before citing new claims. |

## 1. Correct weighted local risk

Use real, noise-whitened data, identifiable scaled coordinates, fixed full-column-
rank J, and deterministic SPD precision W. Let y=J theta + d + epsilon,
E epsilon=0, Cov epsilon=I. Define

$$
L_W=(J^T W J)^{-1}J^T W.
$$

For a specified linear target T, the exact linear mean-square error is

$$
E\|T(\hat\theta-\theta)\|^2
=\operatorname{tr}(T L_W L_W^T T^T)+\|T L_W d\|^2.
$$

Proof: L_W J=I, so the error is L_W(d+epsilon); expand the squared norm and
use the zero-mean noise. The cross term vanishes. The inverse-information
shortcut is only justified by the corresponding noise covariance assumption.
For d=E eta with ||eta||<=1, the worst squared bias is ||T L_W E||_2^2.
This is an operational upper bound only if the **actual** discrepancy belongs
to that set. With pilot-dependent W, this fixed-design calculation is not an
automatic conditional risk formula for reused pilot observations.

## 2. Safe suppression: a sufficient exact condition, not an angle slogan

Write J=[B,N], B_v=(I-P_N)B and assume B_v has full column rank. Consider
W=I-alpha u u^T for a unit u and 0<=alpha<1. If u^T B_v=0, then
W B_v=B_v and N^T W B_v=0. Writing B=B_v+N K shows that the profiled
normal equation for the task is unchanged. In particular the entire local
linear task estimator remains B_v^dagger y, for every y. Suppression is
neutral for this task, not an automatic improvement in that task.

Counterexample to using raw B: B=(1,0)^T, N=(1,1)^T and u=(0,1)^T.
Here u is orthogonal to B but not B_v. For error covariance I+t uu^T,
the profiled information is 1/(2+t), below its original 1/2.
With unchanged white noise, the actual estimator variance remains 2 in this
square example, showing separately why weighted inverse information cannot
be substituted for actual variance.

## 3. A precise impossibility statement

If an admissible error family contains B h and both parameter points are
admissible, (x,d)=(h,0) and (0,B h) produce exactly the same data in the
linear model, with the same noise distribution. Every point estimator has
maximum risk across these two worlds at least ||T h||^2/4. This follows by
averaging the two squared distances around their midpoint. It is a local
two-world lower bound, not a nonlinear global guarantee. Known error mean,
restricted amplitude, independent references or distinct observations can
change the conclusion; a generalized eigenvalue alone cannot.

## 4. How the next algorithm must use this

Continue the already registered sampled-error baselines first: they directly
test whether existing approximation-error methods explain the pilot-mode win.
Then evaluate task-specific risk and lost profiled information at frozen
pilots, separately for geometry and material. Use parameter scales / target
losses explicitly. A task spectrum selects candidate actions, while the full
weighted estimator and its sandwich risk score evaluate them; do not equate
an eigenvalue ratio with final recovery error.

Compare retain, downweight and higher-fidelity actions with identical data,
counting mode construction and refinement cost. An acquisition action requires
new rows that separate the identified ambiguity. Register the action policy
and tests before inspecting their outcomes. Include the raw-angle and
inverse-information counterexamples as regression checks.

## 5. New experiment that prevents a premature positive narrative

All 32 registered calibrate-then-image fits completed and passed input/fixed-
geometry auditing. Rank-one-derived geometry followed by low-only imaging
improves material versus direct rank-one in five of eight case/reference
conditions, but versus low-only in only two. These are four reused development
cases, not eight independent trials. In the high-contrast ellipsoid without
reference, material error worsens from 0.222% (direct rank-one) to 3.233%
(two-stage); even oracle geometry gives 3.307%. Thus residual material bias
cannot be eliminated by simply solving geometry first. This motivates joint
error attribution, not a claim that the supplied spectral policy already works.
