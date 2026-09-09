# A2 implementation cycle

Started 2026-09-05. Updated 2026-09-06. Status: bounded research cycle and audited
draft delivered; scientific submission gates remain unmet. This status does not
mean all A2 algorithmic conditions or journal acceptance criteria were satisfied.

This is a new research version. The three files under `Theory/Questions/A2*`
are immutable research inputs, not an authority for accepting scientific claims.
The earlier geometry-lift experiments remain historical evidence only.

## Locked scope

Execute A2 E1–E5, a higher-dimensional 2D imaging extension, and a suitability
audit of public measured data. No 3D Maxwell or autonomous-driving validation is
claimed. Codex is the scientific orchestrator; the existing AI Scientist is an
isolated implementation/falsification work package, not the final scientific judge.

## Work packages

- `research/delegated/a2_physics/`: corrected full-wave physical model and checks.
- `research/delegated/a2_theory_validation/`: AI Scientist algebra/likelihood tests.
- `research/delegated/a2_literature/`: retrieval evidence and bounded overlap audit.
- This directory: preregistration, integration, parent audits, final experiment
  analysis and manuscript package.

All data-derived numerical claims must point to raw results, commands, seeds,
configuration hashes, model conventions and an explicit interpretation boundary.
No method superiority is presumed. Final tests must remain inaccessible to tuning.

## Conventions

`L_det` is the data-determined SOM cutoff, `r_num` the numerical state dimension,
and `r_free` the dimension of additional independent current nuisance. The primary
physical estimator has `r_free=0`. A current residual penalty is neither a new
measurement nor a license to label its curvature data Fisher information.

## Current deliverables

- Chinese explanation: `communication/A2_RESEARCH_REPORT_ZH.md` from project root.
- Canonical new manuscript: `PAPER_DRAFT_A2.md` in this directory.
- Rendered review PDF and editable TeX: `output/pdf/trispace_a2_research_draft.*`.
- Chinese next-hard-problems prompt: `GPT_PRO_NEXT_THEORY_AND_DESIGN_ZH.md`.
- Critic disposition and scientific gates: `FINAL_REVIEW.md`.
- Reproduction and exact evidence paths: `REPRODUCIBILITY.md`.
- Audited 2,400-run comparison: `results/E4_RESULTS.md`.

The previous main draft carries a superseded-version link; it is preserved,
not silently combined with A2. Original worker measured fits were rejected;
use `results/measured_corrected_pilot.json` plus the parent review instead.

Bottom line: exploratory exact coherent joint calibration improves images;
SOM-specific reduced performance remains untested because the frozen cost
protocol accepts no reduced iterations. No measured antenna-offset recovery,
universal phase basin, or global novelty claim is accepted.
