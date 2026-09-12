# A2 Gaussian parameter-tangent SOM: bounded development protocol

## Scope

This artifact implements Route A from `Theory/MF_TSOMG_DESIGN_ZH.md`: exact
full-wave state elimination followed by a two-operator selection in the
Gaussian **parameter tangent** space.  It does not claim to implement the
explicit-current Route B objective, traditional current-space TSOM, a
certificate, or information creation from a data-null direction.

The material uses two anisotropic additive Gaussian components.  Each has
six real coordinates `(amplitude, cx, cy, log_sx, log_sy, angle)`.  Frequency
dependence is the existing fixed, known lossy law supplied by the A2 physical
module.  There are no calibration nuisance parameters and no neural network.

## Full-wave tangent and metric

For each frequency/illumination, the full-wave discretization solves

\[
 M_l(\theta) j_l = X_l(\theta)e_l, \qquad
 M_l=I-X_lD_l.
\]

The tangent for a scaled real parameter coordinate `q`, with physical update
`dtheta = L dq`, is obtained by differentiating the full state equation:

\[
 M_l\, \partial_q j_l = (\partial_q X_l)(e_l+D_lj_l).
\]

The external and internal operators are respectively

\[
 A=\operatorname{realify}\{W_l S_l K_lL\}_{l}, \qquad
 B=\operatorname{realify}\{Q_lD_lK_lL\}_{l}.
\]

Gaussian material derivatives are analytic (`amplitude`, center, log-width and
angle); only the VIE tangent is solved numerically. `W` is a fixed scalar whitening from the generated complex-noise standard
deviation. `Q` is fixed quadrature normalization `h/sqrt(N*F*Tx)`, so the
second fold measures normalized state/field coupling without changing the
physical parameter metric. `L` consists of frozen positive scales: amplitude
0.25, positions 0.04 m, log-width 0.35, and angle 1 radian.  The same `L`
feeds both operators.  Complex vectors are realified by vertical concatenation
of their real and imaginary parts.

The first fold uses the leading right singular directions of `A`; the remaining
right subspace is `Vd`.  The second fold computes an SVD of `B Vd` and maps its
leading directions back as `V2=Vd Z2`.  Update directions are `[V1,V2]`.
Ranks use a declared relative singular-value cutoff plus an absolute tolerance
in the physical units of each operator, and are recorded in results; for the
compact run a fixed target rank is used only after clipping to the numerical
rank. Each step is clipped to a declared radius in the scaled parameter metric.

`B` is an optimization/selection operator.  The code never appends it to the
data residual, treats it as a new measurement, or reports it as increased
identifiability.

## Explicit-current TSOM variant

The intended Route B objective has one ambient-current coefficient vector per
frequency/illumination and requires a separately constructed current-space
measurement-weak/internal-strong basis and a declared state-residual penalty.
This bounded A2 branch does not yet construct that basis or its penalty
continuation.  Accordingly no `Gaussian + explicit TSOM-current` comparison
is presented.  The result JSON labels this as an implementation limit rather
than treating the tangent-space solver as an equivalent implementation.

## Development experiment

A deterministic two-Gaussian target is generated on the same full-wave A2
geometry.  Data are complex full-wave fields with a single fixed Gaussian-noise
draw.  Half the receivers are used for fitting and the interlaced remainder is
held out.  All methods start from exactly the same perturbed two-Gaussian
initial point and have the same outer-iteration cap.

The methods are:

1. ordinary full-parameter damped Gauss--Newton/LM;
2. first-fold-only parameter SOM (same rank budget as the twofold first fold);
3. twofold parameter TSOMG (`V1` plus selected `V2`).

All residuals are externally measured complex residuals. The run stops on a
noise-discrepancy threshold, a projected-gradient threshold, a rejected trust
step, or the iteration cap. At every accepted
iterate the code records fitting residual, held-out field error against clean
data, material-grid relative error, full-state equation relative residual,
full-wave solve count, wall time, selected ranks, singular values, and internal
operator energy of retained weak directions.  It is one compact development
case, so it is a readiness test rather than a performance or statistical claim.

## Costs and limits

Dense full-wave solves dominate.  With `N` state cells, `P` material parameters,
and `L=F*Tx` right-hand sides, one tangent evaluation stores dense operators
`O(N^2)` and needs a factorization plus `P` tangent right-hand-side solves;
the direct dense cost is approximately `O(N^3 + PN^2 + N^2L)`.  Each SVD works
on the explicit realified tangent matrices and is limited here to `P=12`.

Execution sets BLAS/OpenMP/Accelerate thread limits to two and uses no NN,
training, large sweep, or unrecorded retry.  The use of the same internal
operator for selection does not establish a new observable or prove that the
twofold method beats LM, GSVD, continuation, current-space TSOM/CSI, or a
high-fidelity state model.  Those remain separate gates.
