# Family 5 parent audit: generalized retention derivatives and rank events

Date: 2026-09-03.  Scope: the discrete, whitened/realified `N=16`, smooth
`p=24`, `q=18` model.  This supplement closes a semantic gap in the original
Family 5: ordinary eigenvalues of `K_eff` were checked there, whereas the
paper's normalized information quantity is the generalized eigenvalue

\[
K_{\mathrm{eff}}(X)v(X)=\rho(X)K_{\mathrm{IS}}(X)v(X),\qquad
v^T K_{\mathrm{IS}}v=1.
\]

## Exact command

```bash
.venv/bin/python src/family5_parent_generalized.py
```

Runtime on the recorded Apple-Silicon environment: 2.68 s.

## Formula and numerical representation

For a simple generalized eigenvalue,

\[
\dot\rho=v^T(\dot K_{\mathrm{eff}}-\rho\dot K_{\mathrm{IS}})v.
\]

The finite-prior derivative uses
`E=I-B(B^T B+alpha I)^{-1}B^T`, including `dot(E)`.  The no-prior derivative
uses the constant-rank projector rule

\[
\dot P_B=(I-P_B)\dot B B^\dagger+(B^\dagger)^T\dot B^T(I-P_B),
\]

and `dot(E)=-dot(P_B)`.  Generalized eigenpairs are computed through the
factored SVD representation `A=Q Sigma V^T`, i.e. the symmetric problem
`Q^T E Q`.  This is essential: `cond(A)=5.30e5`, so direct formation and
testing in coefficient-space `K_IS=A^T A` amplifies cancellation.

## Executed results

| check | finite prior alpha=1 | no prior | status |
|---|---:|---:|---|
| factored generalized backward residual | 4.04e-16 | 9.24e-16 | PASS |
| factored normalization error | 6.18e-12 | 4.07e-12 | PASS |
| worst centered-FD derivative relative error at eps=1e-3, three gapped modes | 4.79e-7 | 1.61e-5 | PASS (<1e-4) |
| first-order prediction-error slopes, three modes | 1.9998--2.0017 | 1.9924--2.0041 | PASS (second order) |

The missing denominator-motion term is not cosmetic.  At `eps=1e-3`, omitting
`-rho dot(K_IS)` produces absolute prediction errors from `1.23e-4` to
`5.91e-4` in the finite-prior modes, versus `1.55e-11` to `2.81e-9` with the
correct formula.  In the no-prior modes, the corresponding wrong-formula
errors are `7.55e-6` to `1.70e-3`, versus `4.79e-7` to `1.47e-6`.

Direct coefficient-space residuals (`1.13e-9` finite prior and `7.97e-6` no
prior) are retained as conditioning diagnostics, not used as the backward
stability gate.  This mirrors the factored-residual policy already used in
Family 2.

## Repeated eigenvalue control

The no-prior spectrum has a structural six-dimensional `rho=1` cluster.
Individual eigenvector derivatives are basis-dependent and therefore are not
reported as invariants.  The compressed derivative matrix

\[
V_c^T(\dot K-\rho\dot K_{\mathrm{IS}})V_c
\]

has stable data-space spectral norm `2.30e-10`; the six perturbed eigenvalues
remain at one to at most `1.44e-15` over `eps=1e-4,1e-3,1e-2`.  The naively
formed coefficient-space compressed norm is `1.26e-6`, another explicit
example of ill-conditioning amplification.

## Engineered rank-event falsifier

With `B=U diag(s) V^T`, the control replaces the weakest singular value by
`|t| s_18`.  For every nonzero tested `t`, `rank(B)=18` and the no-prior
retained mass is `9.40716`.  At `t=0`, `rank(B)` drops to 17, the projector
jump has operator norm `1.0`, retained mass jumps to `10.40528`, and the
relative Frobenius information jump is `0.13996`.  Thus the constant-rank
derivative formula is genuinely inapplicable at a rank event; a sampled smooth
trajectory scan cannot certify otherwise.

## Claim boundary

- The generalized derivative formula and its finite-dimensional implementation
  are validated in three simple, positive-gap modes for each prior case.
- The cluster calculation validates compressed-subspace treatment; it does not
  assign invariant derivatives to individual vectors in a repeated eigenspace.
- The rank event is an engineered algebraic counterexample, not an observed
  event on the physical trajectory family.
- `dA` and `dB` still use centered finite differences, so this is not a
  continuum differentiability proof or a uniform trajectory certificate.

Artifacts: `results/family5_parent_generalized.json`,
`figures/family5_parent_generalized.png`, and
`src/family5_parent_generalized.py`.
