# Workshop: Geometry-Lifted TriSpace SOM under Unknown Array Geometry

This workshop is governed by
`research/agentic_input_trispace/RUN_INSTRUCTIONS.md`. Read that file and
`research/trispace_self_calibration/THEORY_SEED.md` completely before choosing
or executing a research action.

The scientific seed is deliberately narrow:

> In a full-wave contrast-source model with an unknown Tx/Rx array pose,
> decompose the receiver-side Green-operator derivative into a canonical
> minimum-norm lift inside a retained SOM current space and an orthogonal
> irreducible data residual. Couple that lift to the linearized state equation
> to expose when a data-equivalent current is physically inconsistent. Test
> whether this graph-based TriSpace construction enables an implementable
> self-calibration algorithm.

The five mandatory experiment families are E1--E5 in `THEORY_SEED.md`:

1. lift identity, rank-truncation behavior, stability, and full-row-rank
   vacuity;
2. pose-dependent SOM projector drift, spectral gaps, and rank events;
3. receiver pseudo-current versus transmitter-induced physical current;
4. data cancellation versus state-equation disambiguation and gauge;
5. an executed reduced nonlinear self-calibration comparison.

The expected useful result is not “all geometry errors are currents.” That
statement is usually trivial. The candidate value is the restricted canonical
lift plus its two witnesses: the irreducible data residual and the state defect.
The run is allowed to narrow or reject this thesis if the experiments do so.

Use the existing corrected Helmholtz implementation and its self-cell/source-
gradient conventions. Run on CPU, retain failures, and keep Phase-I SOM distinct
from an unverified TSOM domain fold.
