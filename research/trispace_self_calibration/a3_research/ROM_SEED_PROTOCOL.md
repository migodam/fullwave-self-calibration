# Matched seed/reuse development comparison

Registered 2026-09-08 before new comparisons. This is development, not final
testing and not a reproduction of the whole classical SOM algorithm family.

Question: do sensing-SOM or Twofold spectral initial spaces yield a practically
useful advantage over physical-RHS or block-Krylov seeds under the same
task enrichment, chart guards, physical parameters and exact trial acceptance?

Use the existing q9/grid16 development scenes 4101–4104 and q49/grid32 scenes
4201–4202. All six methods are rerun serially: direct GN, direct adjoint,
generic-task, sensing-task, Twofold-task, and Krylov-task. Reuse is enabled.
Sensing, Twofold and Krylov seeds have a maximum dimension of 64; all reduced
methods use the same residual/adjoint/material-tangent/pose-RHS enrichment to
dimension 128, at most 14 rounds. The physical-RHS seed has its natural rank.
Sensing rank may be smaller because the stacked measurement matrix has fewer
independent rows. Record actual seed/final ranks rather than force fake rank.

No independent current nuisance variables are introduced. All methods estimate
the same material and pose variables. Frequencies are jointly fitted, without
a method-specific continuation advantage. All have 24 outer iterations and
8,000 full-solve RHS columns as resource ceilings; these are NOT equal wall
budgets. Report wall time including setup separately from offline derivative
audits. A smaller RHS count alone is not a speedup.

The 1e-3 relative state/scattered-output chart guards and exact physical trial
acceptance are unchanged. Every accepted reduced update receives a post-run
high-frequency derivative audit. Activation failure is reported before any
scientific interpretation. Serial timings remain single development runs,
not final-population confidence intervals. Retain all failures. No final
seeds are opened or seed guards disabled.

The new dispatcher must reproduce the original generic/Twofold charts and
verify orthogonality and enriched-chart derivatives before comparisons. The
adapter is process-local; previous source and results remain unchanged.
