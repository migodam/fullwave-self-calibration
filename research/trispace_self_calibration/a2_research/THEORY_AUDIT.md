# Parent scientific audit of A2

Date: 2026-09-05; final update 2026-09-06. Authoring role: Codex scientific orchestrator.
Status: independent algebra audit and bounded numerical checks completed;
algorithmic advantage and submission gates remain unmet. See `FINAL_REVIEW.md`.

## Accepted algebra after independent inspection

The efficient-score contraction proof is valid for the same target/nuisance
family and a parameter-independent observation channel. Its equality condition
is measurability of the coherent efficient score, not of the raw target score.
Use the raw total field as parent when the incident field depends on geometry.

The nested nuisance loss, neutral-admission kernel, bias/variance identity and
shared-map innovation follow from orthogonal projection, least squares and
completion of the square. Their proofs are sound under the stated finite-rank
conditions. These basic tools are prior mathematics; the design translation is
the research candidate, not an unqualified new-information theorem.

For complex additions the admissible real kernel must be invariant under the
complex structure. The intersection ker(F) intersect ker(FJ) has this property;
arbitrary real null vectors do not define a complex SOM subspace.

Physical state elimination has no independent free-current nuisance. Numerical
rank refinement cannot be charged the free-nuisance information loss. Classical
SOM data cutoff and complementary optimization dimension vary oppositely.
The sequential TSOM family P_(S-) V_(D,+) is not generally an exact intersection.

The state-defect bound is local at the candidate material and geometry. It
certifies numerical error, not material-model correctness or the correct phase
branch. Bounds on inverse GN curvature alone do not certify coverage.

## A2-to-algorithm addition: a passive-medium resolvent bound

This is an independently derived, model-specific extension to make A2's state
certificate inexpensive. It is NOT yet a literature novelty claim.

Let the implemented cell-integrated operator be D=k^2 G_D. For equal-area cells
of area h^2 with disk radius a=h/sqrt(pi), its off-diagonal imaginary entries are
k^2 h^2 J0(k|z_i-z_j|)/4. The point-collocation matrix Q with that formula also
on its diagonal is positive semidefinite, because

$$
J_0(k|z_i-z_j|)=\frac1{2\pi}\int_0^{2\pi}
 e^{ik\hat s(\theta)\cdot(z_i-z_j)}\,d\theta
$$

is a Gram kernel. The disk-integrated diagonal gives

$$
\operatorname{Im}D=Q-\delta I,\qquad
\delta=\frac{k^2h^2}{4}-\frac{\pi ka}{2}J_1(ka).
$$

Consequently Im(D) >= -max(delta,0) I. Here Im means the Hermitian imaginary
part (D-D*)/(2i); the matrix is complex symmetric, so it agrees with its
entrywise imaginary part. This conclusion would not hold for an arbitrary
nonreciprocal operator without revisiting the argument.

Suppose the material contrast is chi=(1+i*tau)u, u_i>0 and tau>0, as in the
declared Ohmic law at each frequency. Let U=diag(u), T=U^(1/2), and

$$
A'=(1+i\tau)^{-1}I-TDT,\quad
\alpha=\frac{\tau}{1+\tau^2}-\max(\delta,0)\max_i u_i.
$$

If alpha>0 then Im(A') <= -alpha I. For every v,
alpha||v||^2 <= |v* A'v| <= ||v|| ||A'v||, hence
||A'^(-1)|| <= 1/alpha. Since

$$
T^{-1}MT=I-(1+i\tau)TDT=(1+i\tau)A',
$$

we have the analytic upper bound

$$
\boxed{\|M^{-1}\|_2\le
 \frac{\sqrt{u_{\max}/u_{\min}}}{\sqrt{1+\tau^2}\,\alpha}.}
$$

More usefully, if z=M*j_tilde-b, the data error obeys

$$
\boxed{\|WS(j_{\rm tilde}-j)\|_2\le
 \frac{\|WST\|_2\,\|T^{-1}z\|_2}
 {\sqrt{1+\tau^2}\,\alpha}.}
$$

An inexpensive conservative evaluation replaces ||WST||_2 by its Frobenius
norm; no inverse or state solution is needed for this bound. A corresponding
derivative bound uses the residual identity already proved in A2. For multiple
RHS, square and sum bounds, or use the Frobenius residual norm directly.

Limits: the simple proof needs a homogeneous real background, the stated
quadrature, a common positive loss ratio and positive u on active cells. When
u=0, its current is exactly zero and must be eliminated before applying the
bound. Mixed-sign contrast, nonuniform loss ratio, lossless media or alpha<=0
require a different certificate or direct fallback. A failed bound is not proof
of a singular forward model. Floating-point evaluation is not interval arithmetic;
record a conservative margin and validate it against dense singular-value
calculations as an independent diagnostic, never a hidden estimator oracle.

## Implementation tests executed

1. Compare Im(D) with Q-delta I and verify the Gram lower bound.
2. On physical N8/N16 matrices compare the analytic resolvent upper bound with
   1/sigma_min(M), keeping the latter evaluation-only.
3. For deliberately inaccurate SOM-informed currents test the data-error bound.
4. Sweep tau toward zero, high contrast and coarser grids; alpha<=0 must refuse
   certification and not silently clip to a favorable positive number.
5. Include the overhead of this certificate and all basis construction in timing.

Results: `results/passivity_checks.json` records 12 positive-loss enclosing
checks, six lossless refusals, and an extreme-contrast refusal. The worst state
bound/error ratio is about 210 and derivative ratio about 9,051. `results/
enrichment_probe.json` records 16 nonoracle approximation probes: low-frequency
enrichment can pass at rank 72, but highest-frequency rank 128 remains unadmitted.
The frozen E4 method accepts zero reduced updates. Overconservative skinny-product
charging and stage budgets prevent a useful method comparison; merely charging
the overhead does not establish fair computational discrimination.

`results/final_parent_checks.json` additionally checks efficient physical Fisher
contraction after material elimination and independently rejects the measured
worker's initial propagation/passivity/test-gain implementation. The corrected
measured pilot is a model check only. Original proofs are not expanded into
claims about actual phase-branch coverage, real-array correctness, or global
novelty. Those distinctions are preserved in the manuscript and next Pro prompt.
