# Why additional coherent data can improve information yet worsen recovery

Parent derivation, 2026-09-08. The matrix principles are established linear
estimation/Schur-complement machinery, not a claimed new theorem. Their purpose
here is to make a falsifiable interpretation of the matched full-wave experiment.

## 1. A precise positive statement, with the right scope

Realify and whiten the low-band and added observations. Let h be the physical
target increment and v a common nuisance increment, with Jacobians (B_L,N_L)
and (B_H,N_H). At a fixed true parameter, with the same correct physical and
electronic model, the profiled information quadratic forms are

$$
h^T I_L h=\min_v\|B_Lh+N_Lv\|^2,
\qquad
h^T I_{L+H}h=\min_v\{\|B_Lh+N_Lv\|^2+\|B_Hh+N_Hv\|^2\}.
$$

Therefore I_(L+H) >= I_L in Loewner order. Additional independent coherent
data cannot decrease this local, correctly specified information. This does
not prove better finite-noise nonconvex optimization, smaller deterministic
bias, global branch coverage, or superiority over a different equal-count
acquisition. Replacing a frequency is not nested data addition.

If H=N_L^T N_L is positive definite after removing gauges, set

$$
Q=B_H-N_HH^{-1}N_L^TB_L,\qquad
W=(I+N_HH^{-1}N_H^T)^{-1}.
$$

Then the exact increment is

$$
I_{L+H}-I_L=Q^TWQ.
$$

Proof: write v=-H^{-1}N_L^TB_Lh+u. The low term becomes h^T I_L h+u^T H u;
minimize u^T H u+||Qh+N_Hu||² and use the Woodbury identity. W is positive
definite. A direction gains strictly positive information exactly when Qh is
nonzero: the added response must not be reproducible by the nuisance change
that already best explains the low-band response.

The variational form handles rank-deficient nuisance without pretending H is
invertible. In that case, strict improvement requires that no low-band
minimizer can also make the added term zero. Newly free per-frequency gains
may absorb the added response and remove strict improvement. Changing the
old-data parameter freedom invalidates the fixed-model comparison.

## 2. Information is not risk under model discrepancy

In a full-column-rank local linear model, let

$$
y=J\theta+E\eta+n,\qquad n\sim N(0,I),\quad\|\eta\|\le1.
$$

E is a specified deterministic error set, not a fitted residual. For a fixed
physically scaled target map T and least squares theta_hat=J^dagger y,

$$
\sup_{\|\eta\|\le1}\mathbb E\|T(\widehat\theta-\theta)\|^2
=\operatorname{tr}\big[T(J^TJ)^{-1}T^T\big]
+\|TJ^\dagger E\|_2^2.
$$

Proof: bias and zero-mean noise separate in expected squared norm. The
covariance is (J^T J)^-1 and the largest squared bias is the top singular value
squared. The same formulas apply to geometry or material components via T;
unscaled sums of meters, radians and relative permittivity are not meaningful.

For a scalar counterexample, the old observation has derivative1 and zero
bias, the added observation derivative a and deterministic error b. Old MSE
is1; joint MSE is

$$
\frac1{1+a^2}+\frac{a^2b^2}{(1+a^2)^2}.
$$

For a !=0 it exceeds1 exactly when b²>1+a². Thus larger derivative/information
can coexist with greater actual recovery error. This is a linear-estimation
counterexample, not by itself an electromagnetic construction or a new result.
The independent full-wave experiments must establish the physical relevance.

## 3. Consequences for this project

- The favorable mechanism is not merely high phase sensitivity. It is added
  **nuisance-distinguishable** response under an adequate full-wave model.
- Frequency-shared electronics and a fixed material family are substantive
  assumptions. Added electronics references supply additional information.
- Underresolved high-frequency fields may produce bias that variance-only
  acquisition scores omit. The relevant error is its projection into estimated
  physical parameters, not only its norm or the final data residual.
- A robust acquisition score based on the displayed risk is an implementable
  algebraic rule **only if E is defensibly available**. An adjacent-grid
  difference is not automatically a verified error set, and using the unknown
  reference truth to choose an online action would be leakage.
- Accordingly, this document does not claim an operational continuum-error
  certificate. The matched low/high/repeat experiment tests the mechanism;
  obtaining an honest, useful error set remains the substantive design problem.

Independent numerical algebra checks are in
`results/frequency_information_checks.json`. They validate the identities, not
the availability of E, uniqueness of full-wave recovery, or scientific novelty.
