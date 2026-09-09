# Parent review of the first Family-2 run

Review status: **partially accepted; targeted correction required before
manuscript use**.

The command `.venv/bin/python src/family2_algebraic_spine.py` completed and the
rank identity, principal-angle spectrum, loss-rank bound, PSD ordering,
interlacing, and correctly stacked duplicate controls are useful bounded
evidence.  The following points must remain visible and must be corrected or
qualified before writeup.

1. **Pixel-basis kernel projector gate failed.**  The computed null spaces have
   equal dimension 226, minimum singular value of `N1^T N2` equal to
   `0.999999999999641`, and cross residuals around `1.6e-13`, but the spectral
   projector distance is `8.47e-7`, above the declared `1e-8` gate.  This is a
   real numerical failure of the chosen eigenbasis comparison on the squared
   operator `K_SLAM=C^T C`; do not relabel the gate as passed.  Add a stable
   row-space/null-space comparison and explain the conditioning, while
   retaining the failed original metric.
2. **The small-prior limit was not reached.**  At `alpha=1e-6`,
   `||K_eff-K_SLAM||_F/||K_IS||_F` is about `0.209` (pixel) and `0.137`
   (smooth).  Since the weakest reported pose singular value has square around
   `4.5e-7`, this alpha is not asymptotic.  Extend the sweep downward (with a
   stable SVD-filter expression) or describe the `alpha -> 0` result as an
   analytic conditional limit rather than a numerically demonstrated one.
3. **The singular-prior direction is not automatically a gauge.**  The last
   pose-angle coordinate was merely left unpenalized.  Call it an unpenalized
   nuisance direction, not a residual gauge direction.  Its limiting loss is
   geometry-dependent; compare the large-alpha result to the explicit
   one-dimensional limiting projector/loss instead of using an arbitrary
   magnitude threshold.
4. **Duplicate construction is accepted.**  The implementation correctly uses
   `A_aug=[A;cA]`, `B_aug=[B;cB]`, so the factor is `1+c^2`.  This overrides the
   contradictory wording in the automatically generated child-task prompt.
5. The pixel case has `rank(A)=48=m`, so `range(A)` fills the realified data
   space and the 18-dimensional pose range is wholly confounded.  This is an
   overparameterized finite-data example, not evidence that generic
   full-wave inverse-scattering maps always have this defect.  The 24-column
   smooth-basis case is the nondegenerate comparison.

No Family-2 number proves continuum identifiability, genericity, or a physical
trajectory ordering.
