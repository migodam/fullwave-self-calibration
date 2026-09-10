# Development cycle 2 — finite-candidate scalar acquisition, not certificate

2026-09-10, registered after cycle1 and before these new outcomes. Reuse only
the12 independent-DDA development scenes; no final/tuning data generated. No
claim of fresh statistical confirmation. Compare one scalar complex reference
with one added electric-field projection, same one-complex-observation count.
Hardware/reference/computation costs remain distinct.

At each no-reference fit, construct two constrained competing materials with
eps1 ±.25 (whichever remains within the declared domain) and optimize eps2/x
and the allowed gain against training data. A found pair is an upper bound on
minimum separation, NOT a covered lower bound. Candidate grid:24 full-sphere
directions at radii.14,.2m, existing4 illuminations and3 field components.
Select the scalar maximizing the minimum predicted normalized separation from
the available competing candidates. No ground truth enters selection. Reference
score uses corresponding fitted gain differences and reference sigma .01.
Choose reference if its score is at least the selected-EM score, otherwise EM.
This policy is explicitly a surrogate, never a certified ambiguity-breaking rule.

Baselines: fixed scalar index0 and seeded random scalar; one same paired noise
array generated for the entire prospective grid. Electronic gain assumed common
across base and added measurement, an ideal drift-stable acquisition assumption.
Extra-channel noise has same absolute sigma as base channels; no renormalization
of near-field signal. Count all model-selection evaluations and latency.

Optimize from the same3 starts with same limits as cycle1. Report full task
errors, properly whitened residual statistic2*||r||², fixed-world structural
field on radius.15m (do not move that evaluation surface with estimated x).
All scientific acceptance statuses remain unresolved. Save actual noisy base,
reference and selected/random/fixed measurements for replay. Old cycle1 raw
heuristic flags used an unwhitened statistic and old structural metric moved the
surface; preserve them, produce separate corrected audit, never use those flags
or structural values for scientific conclusions.
