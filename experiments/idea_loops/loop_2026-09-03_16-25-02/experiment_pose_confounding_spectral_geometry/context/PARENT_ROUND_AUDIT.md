# Parent mathematical audit gates for Families 2--5

Status: binding audit constraints added while development round 0 is active.
They supplement `workshop.md` and `PARENT_CORRECTIONS.md`; they do not replace
the experiment plan.  A numerical pass is inadmissible if it violates one of
these conventions.

## Common representation

All algebra below is performed only after noise whitening and realification:

\[
A_R=\sqrt 2\begin{bmatrix}\Re(WA)\\\Im(WA)\end{bmatrix},\qquad
B_R=\sqrt 2\begin{bmatrix}\Re(WB)\\\Im(WB)\end{bmatrix}.
\]

The map and pose parameters are real.  Drop the subscript only after recording
this conversion in the result metadata.  Complex projectors/ranks are not a
substitute for the realified pose range.

## Family 2: algebra and prior gates

For whitened-real `A,B`, use

\[
K_{IS}=A^T A,\quad K_0=A^T(I-P_B)A,\quad
K_X=A^TA-A^TB(B^TB+J_X)^\dagger B^TA.
\]

1. Check the no-prior identities at a predeclared SVD tolerance:
   `rank(K_IS)-rank(K_0) = dim(range(A) intersect range(B))`, and
   `ker(K_0) = {u: Au in range(B)}` by residuals on both sides.
2. On the observable support of `K_IS`, if singular values of
   `Q_B^T Q_A` are in descending order, report generalized retentions in
   ascending order as `rho_i = 1-sigma_i^2`; pad unmatched directions with
   `rho=1`.  Do not call ratios on ordinary eigenvectors of `K_IS` the
   generalized retention spectrum.
3. Check the PSD differences in the correct directions:
   `K_IS-K_X >= 0` and `K_X-K_0 >= 0`.  Also check
   `rank(K_IS-K_X) <= rank(B)`.
4. With eigenvalues in descending order and
   `p=rank(K_IS-K_X)`, check
   `lambda_i(K_IS) >= lambda_i(K_X) >= lambda_{i+p}(K_IS)` wherever the
   lower index exists.  Record absolute and scale-aware residuals.
5. A pseudoinverse Schur complement must state its range condition.  For
   `J_X >= 0`, explicitly verify numerically that
   `range(B^T A) subset range(B^T B+J_X)` (equivalently its nullspace
   residual), even though this follows analytically from positivity.
6. A scaled duplicate means vertical stacking
   `A_aug=[A;cA]`, `B_aug=[B;cB]`, not replacement by `cA,cB`.  Hence the
   data Fisher factor is `gamma=1+c^2`.  This construction preserves normalized
   retentions only for `J_X=0`, or when the prior is also changed to
   `J_X_aug=gamma J_X`.  With a fixed nonzero prior, show rather than suppress
   the change.  If only `cA,cB` are used, the factor is `c^2`, so a check
   against `1+c^2` is mathematically invalid.
7. A singular-prior residual has no universal lower threshold such as `0.1`.
   Its limiting nonzero loss depends on how the unpenalized pose subspace
   couples to `A`.  Record the limiting operator and compare it with the
   analytically predicted residual using that support; if the selected
   unpenalized coordinate happens to decouple, change to a declared coupled
   direction or report the degeneracy rather than tuning a threshold.

## Family 3: gauge and Born gates

1. Roundoff-level `A xi_chi+B xi_X=0` is expected only from an exactly
   equivariant discrete construction.  Pixel interpolation/cropping errors
   must instead be reported versus grid/basis refinement and cannot be called
   exact gauge verification.
2. Exact gauge implies zero retention; zero retention does not imply that the
   direction is a physical group gauge.
3. In Born at `chi_0=0`, verify `B=0` and first-order `K_X=K_IS`, but also
   evaluate the mixed bilinear remainder
   `(D_X A[delta X]) delta chi`.  Never infer that pose errors are harmless.

## Family 4: frequency and trajectory gates

1. Frequency rows share one physical pose increment.  Stack `B_f` vertically;
   do not introduce an independent nuisance copy per frequency.
2. Adding data must be checked through absolute PSD monotonicity of both
   `K_IS` and the efficient information `K_X`.  Normalized retentions need not
   improve coordinatewise because their denominator changes.
3. Include the exact-duplicate control under (a) no prior, (b) fixed prior,
   and optionally (c) jointly scaled prior.  Claim invariance only in (a)/(c).
4. Equal-budget trajectory comparisons are empirical orderings for the tested
   scenes.  Do not turn them into a universal trajectory theorem.

## Family 5: sensitivity and rank-event gates

On a locally constant-rank no-prior stratum,

\[
\dot P_B=(I-P_B)\dot B B^\dagger+(B^\dagger)^T\dot B^T(I-P_B),
\]

and

\[
\dot K=\dot A^T(I-P_B)A+A^T(I-P_B)\dot A-A^T\dot P_B A.
\]

For finite prior with invertible `C=B^TB+J_X`, use `E=I-BC^{-1}B^T`,

\[
\dot K=\dot A^TEA+A^TE\dot A+A^T\dot E A,
\]

\[
\dot E=-\dot BC^{-1}B^T-BC^{-1}\dot B^T
       +BC^{-1}\dot C C^{-1}B^T,
\quad
\dot C=\dot B^TB+B^T\dot B+\dot J_X.
\]

1. Compare these derivatives with centered finite differences away from rank
   events.  Record the spectral gap and smallest retained singular value.
2. For a simple generalized eigenpair normalized by
   `v^T K_IS v=1`, use
   `dot rho = v^T(dot K_X-rho dot K_IS)v`; for a cluster use the compressed
   perturbation operator, not arbitrary eigenvector matching.
3. Near a rank threshold/gap closure, flag the smooth formula as inapplicable.
   Do not hide a discontinuity behind an SVD tolerance.
4. A sampled perturbation sweep supplies an empirical stress test, not a
   uniform Lipschitz certificate.  A certified robust lower bound requires an
   independently valid uniform constant over the entire uncertainty set.

## Required result discipline

Every family writes machine-readable raw tables, exact commands, seeds,
tolerances, runtime, source hashes, and a statement of what the check cannot
establish.  Preserve failures and counterexamples; do not tune them away.
