# Frozen discrepancy-mode development ablation

2026-09-08, registered after the matched low/high/repeat result, before running
these new methods. This is development on the same two scenes, not confirmatory
evidence or an originality claim. Model-error covariance/inflation is established
methodology; closest-prior verification and a broader independent test remain needed.

## Question

Can a fixed numerical-discrepancy direction attenuate harmful high-frequency
model error while retaining useful calibration information, beyond ignoring
high frequencies or uniformly downweighting their entire block?

Use the existing matched-frequency data for seeds8101/8102, both reference
conditions, identical low=(3,6,9) plus high18 observations. Start every new
method from its stored low-only pilot (frozen estimate). All low observations
and reference records remain unchanged. Reference/no-reference comparisons are
distinct additional-data conditions, not an algorithm-only information gain.

At the low pilot, compute high-frequency complex predictions with N16 and N32
inverse models, including estimated electronics. Let d be their difference,
delta=sqrt(2)/sigma [Re d; Im d], and n the number of real high-block channels.
No ADDA truth or high-frequency measurement residual enters this construction.
Freeze the discrepancy mode throughout the subsequent optimization. It is a
proxy, not a verified continuum error enclosure or known noise covariance.

Four continuation methods share the same low pilot and low+high observations:

1. `warm_raw`: ordinary white-noise fit (strong warm-start control).
2. `isotropic`: high-block surrogate C=I+(||delta||²/n)I.
3. `rank1`: high-block surrogate C=I+delta delta^T.
4. `rank2`: E=[delta, sqrt(2)/sigma R(i d)]/sqrt(2), C=I+E E^T.

The last three have equal trace(C-I). Apply the symmetric inverse square root
of C to both residual and Jacobian. Only the high-frequency block is modified;
low blocks and the same low-band electronic reference retain original noise
weights. Strength is fixed at1; no parameter sweep or truth-based tuning.
Rank2 is a declared alternative for complex mismatch direction, not a post-hoc
winning substitute. Rank1 versus isotropic is the primary mechanism contrast.

Sixteen fits: two scenes × two reference conditions × four methods. Same
bounds,35 evaluations and tolerances as the original matched-frequency fits.
Retain low-only and independent low-repeat controls from those same data.
The raw warm control should reach the existing raw endpoint if convergence
is not initialization-dependent; discrepancies require inspection, not hiding.

## Checks and accounting

Before full fits: verify W C W^T=I, equal discrepancy trace across methods,
correct real/imag block placement, and the fixed-weight full residual/Jacobian
by finite differences. New methods may cache/reuse their shared pilot and mode
construction physically, but disclose costs: recorded pilot setup+solve,
new discrepancy construction, continuation, and separate frozen evaluation.
Do not label a sum of historical pilot and new continuation times a fresh
end-to-end wall measurement. No speed-superiority claim from this experiment.

Preserve each accepted endpoint before offline evaluation, use safe resumable
records and distinguish evaluation failure from fit failure. Evaluate the same
held-out low-band sensor and structural fields as before, material/pose errors,
and high-block measurement fit separately. Inspect failures and directional
tradeoffs, not just one favorable average. If the mode hides legitimate pose
information or gives no gain over low/repeat/isotropic controls, retain that
negative result and do not promote the method to a core contribution.
