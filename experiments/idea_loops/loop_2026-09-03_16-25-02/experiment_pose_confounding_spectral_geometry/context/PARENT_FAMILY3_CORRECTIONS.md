# Parent corrections for Family 3

Status: binding corrections.  The automatically generated Family-3 worker
prompt contains false expectations; any output following those expectations
without qualification is diagnostic only and must be corrected before use.

## 1. Discrete gauge residuals

For a simultaneous global translation or rotation in free space, the
continuous infinitesimal relation is

\[
A\xi_\chi+B\xi_X=0.
\]

The stated generator signs are consistent with the passive map action
`chi_g(r)=chi(g^{-1}r)`: translation gives
`dchi=-v dot grad(chi), dp=v`; rotation gives
`dchi=-(Jr) dot grad(chi), dp=Jp, dtheta=1`.

There is **no prior ordering** saying the pixel-basis residual must be larger
than the smooth-basis residual.  Pixel samples represent the declared map
tangent exactly at the grid points, while a smooth basis adds representation
error; either can dominate the separate square-domain/quadrature equivariance
error.  Report both raw values and their refinement behavior without tuning to
an expected ordering.  Roundoff-level residuals are not expected from this
fixed square cell discretization unless exact discrete equivariance is proved.

## 2. Radially symmetric scene

For a radial contrast centered at the rotation center,
`dchi_rotation=-(Jr) dot grad(chi)=0`.  This does **not** make simultaneous
global rotation observable.  Instead, the scene has a rotational stabilizer:
rotating every transmitter/receiver pose leaves the data invariant while the
map tangent is zero, so `B dX_rotation=0` in the continuous model.  It is a
pose-only null direction, not a nonzero map gauge direction and not a reason
to expect a large normalized residual.  Avoid dividing by `||dchi||` or
interpreting a nonexistent map coefficient's retention.  Record
`||dchi||`, `||B dX||`, and the symmetry-breaking discretization residual
separately.

## 3. Anchors and masks

The proposed corner/background/boundary masks do not universally increase
`rho_min` by tenfold.  They break a gauge only to the extent that the forbidden
map tangent has non-negligible support there and the remaining parameterization
cannot reproduce it.  A corner in which the two-Gaussian contrast is
essentially zero may do almost nothing.  Preserve such a failure; do not move
the mask or tune a threshold after seeing the result.  Describe these as tested
finite-dimensional constraints, not physical anchors unless the corresponding
world-fixed information is actually in the measurement/prior model.

## 4. Born empty-background bilinear test

At `chi_0=0`, the Born map is

\[
F_B(\chi,X)=A_B(X)\chi,\qquad B_X=D_XF_B(0,X)=0.
\]

The first omitted joint term is

\[
(D_XA_B(X_0)[\delta X])\delta\chi.
\]

Do **not** compare this coefficient with the pose derivative of the full-wave
map at a unit-norm, non-small contrast.  That derivative contains multiple
scattering and is not required to approach the Born coefficient as the pose
finite-difference step shrinks.

An admissible centered coefficient check is

\[
\frac{F_B(\delta\chi,X_0+h\delta X)
      -F_B(\delta\chi,X_0-h\delta X)}{2h}
\to (D_XA_B[\delta X])\delta\chi,
\]

using `hh.born_forward`, with expected centered-FD error `O(h^2)`.  To show
the term's **joint second-order** role, additionally use one common amplitude
`eps`:

\[
F_B(\epsilon\delta\chi,X_0+\epsilon\delta X)
-\epsilon A_B(X_0)\delta\chi
=\epsilon^2(D_XA_B[\delta X])\delta\chi+O(\epsilon^3),
\]

or an equivalent mixed central difference, and report the scaling.  A
full-wave comparison is optional only under a separate weak-contrast
continuation, where multiple-scattering contributions are measured rather
than silently treated as FD error.

Implementation warning for supplementary runs: if `FD_born_h` is computed
from `F_B(chi,X)=A_B(X)chi` with the *same* `A_plus,A_minus` and the same
centered stencil as `BIL_h`, then `FD_born_h == BIL_h` is an algebraic identity
up to arithmetic.  Its observed error/slope is not an independent `O(h^2)`
validation.  Compare the centered coefficient to an independently more
accurate derivative (for example a five-point stencil), and use the common
joint amplitude `epsilon` to demonstrate second-order scaling.  Likewise, at
fixed nonzero contrast, subtracting `h*(D_X A_B)dchi` from a full-wave pose
increment leaves an `O(h ||dchi||^2)` derivative mismatch in general, so its
small-h remainder need not have slope two.

These checks establish a local finite-dimensional Taylor statement only; they
do not show nonlinear global harmlessness or identifiability at an empty
background.
