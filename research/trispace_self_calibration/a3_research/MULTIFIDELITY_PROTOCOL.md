# Development protocol: calibration-targeted coarse/fine correction

Frozen before this comparison on 2026-09-08. This is algorithm DEVELOPMENT on
the already inspected ellipsoid seeds 8101–8102, not an unopened final test.

Question: can a locally value/tangent-consistent coarse surrogate reduce fine
Maxwell solves without degrading the final discrete objective and actual
calibration outcome, **beyond simple coarse initialization of a fine solver**?

All methods share N64 ADDA data, N32 target inverse discretization, N16 coarse
model where applicable, the same 13 physical parameters, nominal start,
frequencies (6,9,18), 30 dB complex noise and bounds. No new data/reference is
provided to a preferred method. Compare:

1. fine direct least squares from the nominal start;
2. coarse only (accuracy/cost endpoint, not an equal-fidelity competitor);
3. coarse initialization followed by the same fine direct optimizer;
4. coarse initialization followed by value-only defect correction;
5. coarse initialization followed by value-and-tangent correction.

For anchor z_a, the last surrogate is
mu_c(z) + mu_f(z_a)-mu_c(z_a) + [J_f(z_a)-J_c(z_a)](z-z_a).
It matches the fine value and derivative at the anchor; this is established
multifidelity consistency machinery, not a new general theorem. Inner fits
stay within scaled trust bounds, and trial acceptance uses the actual fine
objective. At most eight outer trials are allowed; failed or stalled trials
remain recorded. Fine derivatives used for acceptance/anchors are fully
charged. No inexpensive rigorous continuum-error certificate is asserted.

Serial timings include model setup, the entire coarse warm start, corrections,
fine checks and rejected trials. Common reference-data generation and final
scientific evaluation are separate. Record all RHS work and wall times, final
fine objective, position/material error, sensor prediction and structural-field
error. A benefit against cold fine optimization alone is insufficient if a
simple coarse warm start performs as well or better. No post-hoc final-population
superiority threshold or significance statement is permitted from two scenes.
