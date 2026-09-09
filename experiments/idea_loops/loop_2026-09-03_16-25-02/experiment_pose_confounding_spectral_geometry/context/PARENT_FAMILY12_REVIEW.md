# Mandatory parent review of Family 12

This note is a scientific-meaning override for every later experiment report and
manuscript.  It does not invalidate the numerical arrays in Family 12; it fixes
their interpretation and nomenclature.

## 1. The replicated result is a negative design result

The observed retained-information ordering changes with the scene, frequency
set, and normalization/control.  Therefore Family 12 supports only the scoped
claim that no universal ordering among straight, 90-degree arc, 180-degree arc,
and full-circle paths is justified by these finite-dimensional experiments.
Do not promote any one ordering to a theorem or a scene-independent design
recommendation.

For the original two-blob, raw-control, one-frequency artifact, the recorded
ordering is

`arc90 > straight > circle360 > arc180`,

not `arc90 > straight > arc180 > circle360`.  The latter was a transcription
error in the parent request.  The same-stand-off normalized control instead
records

`arc90 > straight > arc180 > circle360`.

## 2. Finite-prior quantities are not principal angles

For a finite pose prior, the eigenvalues of

`Q_A^T (I + B J_0^{-1} B^T)^{-1} Q_A`

are prior-weighted retention factors.  Applying
`asin(sqrt(rho))` is only an angle-like monotone transform used for plotting; it
is not a principal angle between the map and nuisance subspaces.  In all later
artifacts rename `theta_min_keff` to `retention_angle_transform_keff`, or label
it explicitly as a monotone display transform and never as a principal angle.

Only the no-prior projector case admits the principal-angle identity, with
careful ordering between retention eigenvalues and squared cosines.

## 3. Log volume at destroyed directions

If a retention eigenvalue is zero on the observable map support, the determinant
ratio is zero and its log is negative infinity.  Any implementation that omits
non-positive or sub-threshold eigenvalues computes a positive-support
**log-pseudovolume**, not log-volume.  Rename and qualify this metric.  In
particular, the symmetric ring/circle F1 cells with `rank([A B]) = 40` while
`rank(A) + rank(B) = 42` have a genuine two-dimensional intersection; a finite
sum over remaining positive eigenvalues must not hide the destroyed directions.

## 4. Stable rank identity and the numerical-squaring diagnostic

For the no-prior case, evaluate the loss rank using the factor

`C = (I - P_B) A`

or its whitened/factored equivalent.  The stable identity is

`rank(A) - rank(C) = rank(A) + rank(B) - rank([A B])`

under one declared, scale-aware SVD threshold.  Directly ranking
`K_SLAM = A^T (I-P_B) A = C^T C` squares the condition number and may drop
small but nonzero directions.  Preserve direct-matrix failures as conditioning
diagnostics; do not relax a tolerance until they disappear and do not call them
algebraic counterexamples.

The Family 9 ring/arc case is presently classified as this numerical-squaring
effect.  Keep it distinct from the genuine ring/circle subspace intersection
above.  Family 13 must report raw singular spectra, gaps, tolerance sweeps, and
both direct and factored ranks before making any rank-event claim.

## 5. Stored cosine arrays and ordering

`rho_slam_desc` is stored in descending order, so the elementwise array
`1-rho_slam_desc` is ascending in squared-cosine value.  Do not label it as a
descending principal-angle-cosine spectrum.  Use the independently computed
SVD-based `cos2_desc_principal_angles` array for principal-angle plots and
comparisons.

## 6. Trajectory-control language

Control A equalizes a continuous design-length formula before discretization;
the sampled path lengths and radii are not identical.  Control B fixes standoff
only for the curved paths; the straight path has varying standoff and a matched
midpoint.  Neither control equalizes all of path length, standoff, row count,
total energy, aperture, and conditioning.  State exactly which nuisance factors
are controlled and treat the remaining differences as limitations.

## Required manuscript wording

Family 12 is evidence for protocol sensitivity and against a universal
trajectory ordering.  It is not evidence for a globally optimal path.  Report
principal angles only in the no-prior projector case, prior-weighted retention
otherwise, and separate true subspace intersections from finite-precision rank
loss caused by forming normal matrices.
