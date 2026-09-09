# Mandatory Codex parent audit for development round 3

Date: 2026-09-04.  This file overrides conflicting experiment numbering or
claim wording in `current_idea.json`.  Read it before editing or running code.

## Existing artifacts are immutable inputs

- Families 1--9, Family 4b/4c, Family 5b/5b2, and all `parent_*` supplements
  already exist.  Do not overwrite or silently reinterpret them.
- `family10_*` already denotes the nonlinear/affine covariance toy.  Its
  fixed-true-pose Monte Carlo target is the sandwich covariance `P_samp`, not
  the Bayesian marginal `K_eff^{-1}`.  Read
  `notes/family10_protocol_correction.md`; the original comparison is retained
  as a negative/protocol-correction result.
- `family11_*` already denotes the multi-scene SNR/frequency-set study.  Read
  `notes/family11_snr_diversity.md`; extend it only under a new `family11b_*`
  name and never overwrite the existing source, JSON, figures, or report.
- Read `notes/family5_parent_generalized.md` and
  `results/family5_parent_generalized.json`.  The central trajectory result is
  the generalized derivative
  `dot(rho)=v^T(dot(K_eff)-rho dot(K_IS))v`, not merely the ordinary
  `K_eff`-eigenvalue derivative.  Repeated eigenvalues require compressed
  derivatives; no-prior projector formulas require locally constant rank.

## Round-3 task names and scientific gates

### Family 12 -- trajectory replication and mechanism

Use new files `family12_trajectory_replication.*`.  Compare at least the
two-blob, ring, low-contrast, and offset/edge scenes; at least two frequency
sets; and straight, arc90, arc180, and circle360 paths.  Record both:

1. an equal-continuous-path-length/measurement-count control; and
2. a same-standoff/per-pose-whitened-energy control.

Record retained mass, log-volume, `rho_min`, principal angles, rank(B), shared
pose-compensation residuals, and actual discrete polyline/closed-loop lengths.
Do not seek a universal ordering.  The useful result is a replicated reversal
or an explicit account of which metric, scene, standoff, and frequency changes
the ranking.  Do not tune trajectories after observing the target metric.

### Family 11b -- optional frequency deepening

The required multi-scene 2/3/5-frequency and bandwidth study already exists.
Only add `family11b_*` if it answers a missing mechanistic question with an
equal-total-measurement and equal-total-whitened-energy design.  More points on
the same curve are not sufficient.  Duplicate controls must jointly scale the
prior when claiming finite-prior invariance.

### Family 13 -- rank-boundary autopsy

Use new files `family13_rank_autopsy.*`.  Reproduce the ring-scene machine-rank
boundary from Family 9.  Report raw singular spectra, numerical rank across a
predeclared relative-tolerance grid, stable/effective-rank diagnostics, and
backward residuals for the factorized identities.  Determine whether the
6/7 count was a discontinuous mathematical rank change, a tolerance
classification, or an unstable direct subtraction.  Do not change a tolerance
after seeing the answer and call that a theorem.  Include the engineered rank
loss from `family5_parent_generalized` only as an algebraic falsifier, not as a
physical event.

### Family 14 -- robustness scope, not pseudo-certification

Use `family14_robustness_scope.*` only if new computation is necessary.
Otherwise produce a precise manuscript audit.  The affine tangent bound is a
finite-dimensional structural certificate conditional on the implemented
Jacobian.  The full nonlinear Helmholtz operator has no interval/uniform
certificate here.  Zero sampled violations cannot be called certification.
Integrate the generalized derivative and cluster/rank-event results from the
parent supplement.

### Family 15 -- Born analytic/literature control

Use new files `family15_born_control.*`.  A self-derived Born/weak-scattering
formula checked against the implementation is an analytic control, not an
external benchmark.  A published number may be compared only if geometry,
normalization, noise, parameterization, and units match.  Otherwise provide a
formula-level comparison to a verified source and explicitly state that no
quantitative external benchmark was possible.  Do not invent citations or
claim that internal agreement is independent validation.

## Writeup gate

Produce `manuscript/main_round3.tex` and `main_round3.pdf` without overwriting
Round 2.  The paper must lead with the finite-dimensional theorem/identity and
its correct scope.  Keep these layers separate:

1. exact finite-dimensional linear algebra;
2. conditional/discretized full-wave differentiation;
3. executed numerical evidence and preserved counterexamples;
4. open continuum, nonlinear-global, real-system, and novelty questions.

Do not use “first”, “no prior work”, “certified nonlinear robustness”, or a
universal trajectory/frequency design claim.  The 429-limited novelty search
is a retrieval limitation, not a novelty verdict.  The final stdout summary
must name failures and uncertain claims, not only passes.
