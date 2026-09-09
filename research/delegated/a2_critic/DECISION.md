# Decision: venue readiness and next experiments

## Readiness verdict

NOT READY for submission to TAP or TGRS, as either a method/experimental paper
or a theory paper with experimental support.

- The paper's own validation protocol (A2_PRASC_SOM_VALIDATION_PROTOCOL.md,
  line 177) requires, for TAP, a controlled 3D dyadic-Maxwell or
  measured-array validation with known phase-center offsets plus a
  clock/coupling mismatch control; for TGRS, an imaging-relevant task with a
  larger map representation and matched-budget comparative evidence. Neither
  is met.
- The E4 primary endpoint is void (0/240 successes everywhere; zero reduced
  iterations accepted), so the "matched-budget comparative evidence" is
  absent.
- The measured pilot has no position-error ground truth and exhibits eps_r
  non-identifiability under model misspecification — it does not satisfy the
  TAP measured gate.
- The only positive nonlinear result is an exploratory, exact-direct,
  6-scene study; it does not demonstrate any SOM-specific benefit.

The theory core is sound (independently re-verified) and the limits-and-design
content is defensible, but that is a paper fragment, not a submission-ready
manuscript. Conditional route: a theory/diagnostics paper scoped strictly to
the typed rank/bias/acquisition budget and the model-conditional passive bound
COULD become viable, but only if the failed discriminating experiment is
replaced by one that exercises the method and one of the two venue evidence
gates below is genuinely met.

## Three most valuable next experiments (minimum discriminating)

1. Discriminating E4 re-execution with a dimension-aware cost model.
   Re-derive RHS-equivalent units rank-dependently (no skinny-to-dense
   conversion), pre-register stage budgets that guarantee >=1 accepted
   low-frequency step and a completed final-data stage, and re-run the
   compound gate on fresh scenes (no seeds 1-10/1001-1020). The gate is only
   meaningful if >=1 reduced step is accepted; otherwise the protocol itself
   is falsified, which is still a publishable negative if pre-registered.
   Fixes F1/F2.

2. Certificate tightness ablation targeting a useful high-frequency operating
   point. Current bound/true ratios reach 210x (state) and 9,051x
   (derivative) at k=4pi, so rank-128 is unadmitted. Add sharper estimates
   (per-cell zeta, Lanczos/power-iteration ||WST||, block-diagonal inverse
   bounds) with a target of <=10x tightness, and demonstrate an
   admission-vs-rank/cost tradeoff at the top frequency. Without this, guarded
   reduction has no demonstrated operating regime. Fixes F3.

3. Known-offset measured or 3D dyadic-Maxwell calibration validation with a
   mismatch control. A controlled dataset with known phase-center offsets
   plus one clock/timing and one coupling perturbation, matched budgets, and
   joint calibration vs fixed-pose. This is the only route to the TAP gate;
   the Fresnel pilot cannot serve (no offset truth, plane-wave
   misspecification at 4-8 GHz). Fixes F7.

Cheap companion: physical efficient-Fisher check on N=8 tangents (closes F6,
   the abstract's "even after nuisance elimination" claim).

No further runs should be framed as validation until at least experiment 1
yields an accepted reduced iterate under a pre-registered protocol.
