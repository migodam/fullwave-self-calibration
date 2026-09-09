# Registered correct-model frequency attribution control

2026-09-08. Registered before generating these new means or fitting them.
This is an intentionally shared-model control, NOT independent Maxwell validation.
Its purpose is to test whether the observed high-frequency degradation persists
when reference/inverse discretization discrepancy is removed. Results will be
descriptive on the same two inspected scenes; no success-rate inference.

## Frozen intervention

Use seeds 8101/8102, saved true parameters, receiver arrays, original absolute
sigma, low_noise, extra_noise and noisy low-band electronics references from
research/delegated/a3_matched_frequency/data_*.npz. Do not regenerate random noise.
Generate means with ExplicitFrequencyModel N32, the SAME model as the estimator,
at low=(3,6,9) and high=18, including the true electronics. Held-out mean fields
also use N32 at the same held-out receivers; structural fields omit electronics.
The ADDA means are the sole intentional substitution relative to the original
matched control. Absolute sigma stays fixed even if the resulting relative SNR
changes slightly. Save newly generated means and data in a new directory only.

Compare low, low+high, low+independent k6 repeat, each with and without the same
saved reference: 12 fits. Exact same initialization (epsilon=2, others zero),
bounds, parameter scales, 35 evaluations and tolerances as matched_frequency.
No warm starts, weighting, gain freedom changes, new material priors, or tuning.
All truth is used only to generate synthetic data and compute frozen metrics,
never supplied to the optimizer or used to select a start.

## Checks and interpretation

Before full collection, check low/reference identity across acquisition choices,
exact use of saved standardized noises and sigma, field shapes/electronic
conventions, and unchanged imported source/data hashes. Use already validated
physical derivative implementation; do not alter it. Save all endpoints before
evaluation, distinguish fit failure from evaluation failure, retain failures,
and use atomic resumable records with exclusive locks and source/data manifests.
Run serially only after the discrepancy-weighting timing collection completes.

Evaluate the same pose/material and held-out low-band sensor/structural metrics.
Compare both within this correct-model control and with the original independent-
code matched experiment. If high addition helps here but hurts there, it supports
model discrepancy as an explanation in these cases; it does not uniquely exclude
nonlinear interactions or prove asymptotic/statistical optimality. If harm persists,
retain that result and inspect optimization/noise rather than assuming bias.
Either way, single-noise realizations cannot refute expected-information monotonicity.
No measured-data or continuum-accuracy claim may use this shared-model control.
