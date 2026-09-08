# What the frozen discrepancy weight does and does not mean

Parent mathematical audit, 2026-09-08. These are standard weighted least-squares
identities, not claims of new theory. The registered experiment tests their
usefulness in the declared electromagnetic setting.

## A constrained error fit, not extra information

For a real-whitened residual r and a frozen mode matrix E, let C=I+EE^T.
Completion of the square gives

$$
\min_z\{\|r-Ez\|^2+\|z\|^2\}
=r^T(I+EE^T)^{-1}r.
$$

Indeed z=(I+E^TE)^{-1}E^Tr and the Woodbury identity gives the result.
The proposed weighting therefore permits a penalized error component along E.
It does not create additional phase constraints. It suppresses mismatch along
that component but can also suppress genuine geometry/material signatures.
This is precisely why an isotropic control and a low-only control are necessary.

The rank-two protocol uses E=[delta,R(i d)*sqrt(2)/sigma]/sqrt(2). Realification
gives R(d)^T R(i d)=0 and equal norms. For nonzero d its two covariance increments
are both ||delta||^2/2, while rank one has the single increment ||delta||^2.
Both have the same total trace as the isotropic increment. Equal trace does not
make the directional statistical information equal; that is the intervention.

## Do not report the surrogate Hessian as a physical Fisher matrix

First consider a deterministic design: let true whitened data be y=J theta+b+n,
with zero-mean Cov(n)=I noise, and let P be a deterministic block precision
(identity on low/reference, C^-1 on high). For full-column-rank J,

$$
A_P=(J^TPJ)^{-1}J^TP,\qquad
\widehat\theta-\theta=A_Pb+A_Pn.
$$

In this fixed-design setting the covariance is A_P A_P^T, not generally (J^TPJ)^-1.
For a fixed target scaling T and a separately justified error set b=F eta,
||eta||<=1, the worst fixed-design linear squared risk is

$$
\operatorname{tr}(T A_P A_P^T T^T)+\|T A_P F\|_2^2.
$$

One may use (J^TPJ)^-1 as covariance only under the additional model that the
actual stochastic residual has covariance P^-1. The coarse/fine difference
alone does not establish that model. In our pipeline the pilot, and therefore
P, depends on reused low-band data. Conditioning on that pilot generally changes
the low-block noise mean and covariance; it does not restore the fixed-design
assumptions. Thus the displayed formula provides neither automatic pilot-
conditional nor unconditional uncertainty for the complete data-dependent
estimator. Empirical recovery errors remain the primary outcomes. Independent
pilot data would remove this particular dependence, but would be extra data and
a different acquisition protocol; nonlinear uncertainty would still need care.

The sign of d is immaterial for C but would matter for a deterministic mean
correction; no mean correction is part of this registered experiment. The mode
is frozen, so neither its derivative nor a likelihood log determinant varies
within a fit. Recomputing modes during optimization would be a different method.

## Prior-evidence checkpoint

ScholarQA focused paper search on 2026-09-08 at 08:35:50 UTC used the query
`electromagnetic inverse scattering Bayesian approximation error` with five
requested records and no snippets. It returned HTTP 429, zero candidate IDs,
and exit status 1. No cited paper or novelty conclusion was obtained. Expansion
was stopped rather than treating retrieval failure as absence of prior work.
The exact failure was: Semantic Scholar API error (429): Too Many Requests.
Closest-prior verification remains open; the algebra above is independently
checkable but does not substitute for literature evidence.
