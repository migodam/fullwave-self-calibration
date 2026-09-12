# Route-B explicit-current development result

## Executed scope

`code/a2/explicit_current.py` ran a deterministic, resumable development set
of 10 two-Gaussian full-wave N32 scenes: five separated/five close, and five
moderate/five strong.  Each scene has one paired 1% complex-noise draw, two
frequencies, six transmitters, a half-aperture of 24 receivers, alternating
train/held receiver split, a common perturbed Gaussian initial point, and a
maximum of 30 outer steps.  The sealed N64/N128 frozen campaign was not opened.

The frozen input/options and source hash are in
`runs/a2/explicit_current/frozen_development_config.json`; per-case histories,
snapshots, conditioning, rank transitions and actual cost ledger are in
`runs/a2/explicit_current/results.json`.

## Route B that was actually tested

For each frequency and transmitter, the deterministic current was

\[
j_{det}=S_{train}^*(S_{train}S_{train}^*)^{-1}y_{train}.
\]

The weak basis was formed as the orthonormal range of

\[
(I-S_{train}^{\dagger}S_{train})V_D,
\]

where `V_D` is the leading matrix-free right singular range of the current
propagator `D`.  Thus added current coordinates do not directly improve the
training data residual.  At fixed Gaussian material, coefficients are obtained
by ridge-regularized least squares on the full VIE state residual.  The weak
rank begins at 8 and doubles to 16 then 32 if the declared relative state
residual stays above 0.08; it is capped at 32.  Gaussian material uses a
profiled finite-difference state Jacobian, a scaled trust radius, and stops on
state discrepancy, gradient, rejected step, or the 30-step cap.

This is a CSI-family state-consistency experiment, not a proof that this basis
is traditional TSOM or an implementation of the full Route-B objective with a
validated state-penalty continuation.

## Observed development result

| method | median material relative L2 | median held scattered error | median wall time | stops |
|---|---:|---:|---:|---|
| explicit-current TSOM/CSI | 70.75% | 7.17% | 3.10 s | 4 gradient, 6 rejected |
| ordinary full-wave Gaussian LM | 0.50% | 0.31% | 1.41 s | 7 discrepancy, 3 rejected |
| adaptive parameter-tangent TSOMG | 0.59% | 0.34% | 2.74 s | 8 discrepancy, 2 rejected |

All ten explicit-current cases expanded to rank 32. Their final state residuals
remained 0.56--0.67, far above the declared 0.012 state-discrepancy threshold.
The rank expansion therefore did not make this explicit-current formulation
physically consistent enough to recover material. It is a negative development
finding, not evidence against current-space TSOM in general.

The adaptive tangent candidate began at first/second rank 4/0, admitted only
when its selected external-gradient capture had a declared deficit above 0.15,
and never appended internal `B` to a data residual. It approached LM quality in
this bounded set but took longer because it constructs all full-wave tangents
and the `B` operator up front. It does not support a speed claim.

## Concrete implementation limits

The explicit scheme fixes `jdet` from noisy partial-aperture data. Its weak
basis is selected from the discrete `D` SVD and an exact measurement-null
projection, but it does not use a jointly optimized data relaxation, a
validated penalty schedule, an adjoint current basis, or independent-grid
state certification. These are specific next implementation questions; the
current negative result must not be hidden by reporting only held data fit.
