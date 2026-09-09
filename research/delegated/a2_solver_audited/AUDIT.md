# Pre-test correction ledger, 2026-09-06

Original worker output is retained unchanged in `../a2_solver`. No final scene
1001–1020 has been used during these corrections.

1. Correct chain rule: z_x=T x implies gradient_z=T^-1 gradient_x. Both exact
   coherent and exact Rice objectives are now independently checked in optimizer
   coordinates. The old physical-coordinate FD test did not cover this bug.
2. Reduced state/derivative blocks stack by square-sum, not sum of norms. Scaled
   gradient certification uses the same optimizer norm. This improves a valid but
   unnecessarily pessimistic bound, without using true approximation errors.
3. A proposed reduced trial must itself pass the certificate; every accepted
   reduced step is checked for descent in the exact physical objective. No
   unverified trial survives a budget interruption. Basis is frozen in a chart;
   derivatives approximate physical tangent equations, not a moving-proxy map.
4. Extra operator products are charged, conservatively counting reported skinny
   products as full matvec equivalents; raw counts and setup/wall time are
   published. SVD and LU are not free but are charged to wall time, per A2.
5. A single reduced probe, or direct fallback when a probe is unaffordable,
   replaces an exhaustive promotion chain that spent the budget without a step.
   This is a computational safeguard, not evidence of a SOM advantage. No oracle
   singular value or true residual is used to admit a chart.
6. Low-pass width now means Gaussian standard deviation ell_res, not ell_res/4.
   Reference unknown-map pose covariance is the pose block of the joint inverse,
   transformed by the lever-arm metric; it is not known-map pose covariance.
7. Common initialization is explicitly the nonoracle constant alpha=0.5, with
   zero fitting budget. The old 32-RHS low-frequency procedure always fell back
   after spending work, and would leak coherent information into a phaseless
   method if it succeeded. This is a pre-test amendment to the A2 initialization
   procedure, not retrospective selection. Truth is never used in initialization.
8. Stage deadlines are deliberate continuation events. Stopping BEFORE the next
   unaffordable evaluation is not an optimizer-failure criterion. A completed
   budget policy must evaluate the final four-frequency stage and return its last
   accepted state. Missing that stage or violating any task/prediction criterion
   is failure. Status alone can never establish success. Published errors are
   unconditional, including all failures.
9. Final tests require a frozen config and code hash. Bootstrap pairing is unique
   by (seed, initialization); 20 seed clusters, 10,000 resamples. Primary
   Bonferroni one-sided levels are 1-.05/3, with per-aperture normalized margins.
   Secondary comparisons and exploratory larger-budget/high-dimensional studies
   cannot rescue a failed primary result.

Limitations before tests: the implementation is SOM-informed physical reduction,
not Chen's classical independent-current SOM. Its numerical rank is not r_free;
r_free=0 throughout. Its sensing numerical rank is NOT L_det, which needs a
data/noise cutoff. Its local covariance full-rank diagnostic is NOT a calibrated
frequency-admission rule, nor an implementation of the entire A2 controller.
The direct_control ablation shares chart restart cadence, NOT active acquisition.
Source-extension and full nonlinear eight-action acquisition are not implemented.

Frozen tuning will use seeds1–10, one deterministic representative initialization
per seed, four candidates per method; primary200 and secondary800 use that same
selected configuration. This bounded tuning is disclosed, not called exhaustive.
