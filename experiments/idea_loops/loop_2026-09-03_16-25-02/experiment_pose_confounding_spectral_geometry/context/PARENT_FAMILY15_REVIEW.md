# Mandatory parent review of Family 15

This note controls the interpretation of the Born/full-wave experiment and
overrides any stronger language in generated reports or manuscript drafts.

## 1. Separate known-pose and pose-eliminated comparisons

The comparison of `A_born` with `A_full`, and consequently of
`A_born^T A_born` with `A_full^T A_full`, is a legitimate internal check of
the map tangent and known-pose local information under the declared discrete
model.

If the no-prior or finite-prior calculation uses `A_born` together with the
**full-wave** pose Jacobian `B_full`, it is a mixed-tangent isolation control:

`(A_born | B_full)`.

It isolates the effect of replacing the map tangent while holding the
nuisance tangent fixed.  It is not the Fisher/Schur geometry of a complete
Born forward model and must not be called `Born SLAM information`, `Born
retention`, or a full Born-versus-full-wave comparison without the qualifier
`mixed tangent, fixed B_full`.

A genuine Born pose-elimination comparison requires the analytic or
finite-difference Born pose Jacobian

`B_born = D_X [G_S(X) D_chi e_inc(X)]`,

including all pose-dependent Tx/Rx terms in the implemented model.  If this is
not computed and validated, list it as an open control rather than inferring it
from `B_full`.

This missing control has now been implemented additively as Family 15b.  Read
`notes/family15b_born_pose_control.md` and
`results/family15b_born_pose_control.json`.  The analytic `B_born` agrees with
centered finite differences with an observed step-error slope 2.0000; for all
nonzero contrast scales its range projector is invariant to about `4.6e-14`.
Use Family 15b, not the mixed Family 15 rows, for complete Born-pair claims.

## 2. Born empty-background scaling

The Born map Jacobian is independent of the contrast scale in the linear Born
model, whereas the Born pose Jacobian is proportional to the scene contrast
for the current measurement model and vanishes at `chi=0`.  Thus a weak-
contrast sweep of `(A_born | B_full)` does not by itself reproduce the complete
Born empty-background nuisance geometry.  Keep the already established
bilinear second-order statement separate from any mixed-tangent spectrum.

## 3. Statistical meaning of the Diong et al. neighbor

Verified bibliographic entry:

M. L. Diong, A. Roueff, P. Lasaygues, and A. Litman, "Impact of the Born
approximation on the estimation error in 2D inverse scattering," *Inverse
Problems* 32(6), 065006 (2016), DOI
`10.1088/0266-5611/32/6/065006`.

The published abstract states that it characterizes the mean and variance of
the linear Born maximum-likelihood estimator and compares the resulting
estimation error with the full-model Cramer--Rao bound.  Therefore:

- `K_IS,born^{-1}` is at most the covariance/variance expression for a
  correctly specified, full-rank linear Gaussian Born model under the exact
  parameterization and noise convention.
- When full-wave data are fit by a Born model, approximation-induced bias is
  part of MSE and is not represented by `K_IS,born^{-1}` alone.
- This experiment neither reproduces the paper's geometry, scalar
  parameterization, noise, units, nor numerical values.  The citation is a
  conceptual/statistical neighbor, not an external benchmark or independent
  validation.

## 4. Metric names

- In the no-prior projector case, principal angles may be reported using the
  SVD-based subspace construction and consistent ordering.
- Under a finite pose prior, any `asin(sqrt(rho))` quantity is only an
  angle-like monotone retention transform, not a principal angle.
- If zero/sub-threshold retention values are omitted from a logarithmic sum,
  call the result positive-support log-pseudovolume.  A true log-determinant
  ratio is negative infinity when a retained direction is exactly zero.

## 5. Permissible conclusion

Family 15 may establish finite-dimensional convergence trends for the map
tangent and known-pose information and may document how the mixed control
changes retained metrics.  It cannot validate a physical Born error law, a
complete Born-SLAM Fisher model, an estimator MSE, a continuum asymptotic, or
the numerical conclusions of Diong et al.  Family 15b supplies the previously
missing internal `B_born` control, but no externally matched model or estimator
experiment has been supplied.

Family 15b additionally establishes the following exact finite-dimensional
Born scaling statement.  For `F_B(chi,X)=A_0(X)chi` and `chi=s chi_bar`,
`A_B(s)=A_0` and `B_B(s)=s B_1`.  Consequently the no-prior pose projector and
Born `K_SLAM` are constant for every `s != 0` but jump at `s=0` unless
`Range(A_0)` is orthogonal to `Range(B_1)`.  With `J=alpha I`, `alpha>0`, the
loss is continuous and bounded by
`s^2 ||A_0^T B_1||_2^2 / alpha`.  Preserve the singular-prior caveat and do
not interpret nonzero-s projector invariance as practical information at
vanishing signal.
