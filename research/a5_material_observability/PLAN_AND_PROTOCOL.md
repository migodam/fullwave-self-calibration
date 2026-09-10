# A5 round 1: integrate A4/A4_2 without changing the target

Registered 2026-09-10 before collecting new results. Status ACTIVE, not TAP-ready.
Target: material permittivity recovery with unknown calibration, not geometry
recovery after deliberately profiling all material information away.

## Source boundary

A4_2 read completely (2050 lines). A4 packaged THEORY_A4, ALGORITHM_A4,
TAP_REVIEW_AND_CLOSURE and PROVENANCE_AND_CORRECTIONS read completely.
A4 original code and frozen results remain under public_release/research/
a4_reliability_v1; they are not modified. A4_2's literature judgments and value
ratings are hypotheses, not a completed independent novelty review.

## Round-1 falsifiable questions

1. Does A4's positive geometry window imply material information? Test the
   opposite exact prediction: all single-electric-mode radial-material changes
   remain hidden by independent complex gains, at every range/frequency.
2. Does multiple scattering always resolve material/gain scaling? No universal
   claim is allowed. Compare the scalar-mode counterexample with two-sphere
   full-Maxwell data under fixed, declared gain sharing.
3. Does a nonzero material tangent mean finite changes can be distinguished?
   Compute both profiled local log-scale sensitivity and finite-pair complex-
   gain residuals; no grid search is a global stability upper bound.

## Frozen small mechanism experiment

Use existing A3 Treams field adapter, plane-wave illumination, two nonoverlapping
spheres centered(-.06,0,0),(.055,.02,0), radii(.035,.025)m. Contrast template
(1+.03i, 2+.05i) multiplied by a real positive scale s. Fixed k18; s in
(.01,.03,.1,.3,1,2); receiver radii(.2,.6,2)m,12 full-sphere receiver positions,
four incident polarization/direction channels, same angular aperture.

Use lmax3 and lmax4 for every baseline, derivative and finite alternative.
Log-scale central derivative step1e-4; alternative scale1.25s. Incident amplitude
is fixed; do not renormalize the forward data. Record absolute sensitivity and
relative sensitivity separately, not a single claimed hardware noise model.
Compare one complex gain shared over the complete block, separate gains for
each illumination, and independent per-entry gains. The last is an exact
negative control. A linear Born scaling control uses s times a fixed reference
field solely to verify homogeneity; it is not labelled a physical Born solver.

No optimization or recovery success claim in this round. No statistics over18
deterministic parameter/range combinations. Multipole differences are reported,
not silently discarded or described as continuum certification. Count complete
runtime; one CPU job, no new services/GPU/cloud infrastructure.

## Next gates, not assumed completed

- Independent A4 regression rerun and source/hash audit.
- Correct theory bridge: scalar-factor counterexample; observable projective
  shape condition; global ambiguity with bounded, nonzero gain families.
- Literature: use ScholarQA for nearest-source verification; unresolved retrieval
  is not a novelty pass. Existing far-field theorem is not a new A5 contribution.
- Then choose extra calibrated modes, shared frequency response or independent
  reference to open material directions; verify actual internal/external material
  recovery in more than a one-parameter known-support example.
- Keep anisotropy/window splitting as a later hypothesis, not an assumed result.

The installed Agentic-AI-Scientist previously ran real loops; DeepSeek later
exhausted balance. Its currently inspected provider router has no GLM branch.
Do not relabel parent scripts as a fresh installed pipeline execution, change the
default provider, or retry the exhausted provider merely to produce activity.
Local scientific work continues; native provider resumption remains a separate
recorded infrastructure task, not scientific completion.
