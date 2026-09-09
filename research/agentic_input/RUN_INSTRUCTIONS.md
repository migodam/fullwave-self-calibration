# Autonomous run instructions

This is an explicitly requested autonomous research run. Its required endpoint
is a defensible paper draft backed by executable, CPU-only theory-validation
experiments. Preserve the user's existing theory; do not replace it with an
unrelated algorithm or a generic deep-learning project.

## Scientific objective

Develop and test the paper-level thesis that full-wave volumetric
inverse-scattering SLAM admits a pose-confounding effective map-information
geometry based on the whitened, realified map and pose Jacobians. The candidate
selection may sharpen the claim, discover a counterexample, or narrow the scope,
but it must remain within this physical and mathematical program.

## Non-negotiable execution requirements

1. Read and obey the complete workshop material below. It includes the theory
   ledger, exactly five experiment families, and open-question routing.
2. Implement one unified 2D scalar Helmholtz/contrast-source harness with
   Born/full-wave, frequency, trajectory, prior, and perturbation switches.
3. Attempt all five experiment families. A family may end in a documented
   negative result, corrected condition, or computational limitation; it may
   not be silently omitted. Run cheap falsifiers before expensive sweeps.
4. Use Apple Silicon CPU. Use MPS only if PyTorch supports an incidental task
   naturally. Do not install or assume CUDA, GPU Docker, a daemon, scheduler,
   database, web UI, or a new orchestration framework.
5. Keep current-space `G_S`/SOM modes distinct from the map Jacobian `A` and
   map-information eigenmodes. Keep deterministic nuisance pose, random pose
   marginalization, and bounded execution mismatch distinct.
6. Every numerical claim must have exact commands, deterministic seeds,
   tolerances, configs, raw result tables, and figures. Machine-precision rank
   checks and physical dominant-mode thresholds must never be mixed.
7. A passed finite-dimensional experiment supports only the tested model. It
   does not prove continuum transfer, global nonlinear convergence, practical
   online SLAM success, or global novelty.
8. Do not promote unproved continuum, transversality, rank-event, resonance, or
   robust-design claims into the theorem list. Put them in limitations or future
   work unless a valid proof is independently supplied.
9. Literature novelty is retrieval-bounded. Never write `first`, `no prior
   work`, or an exhaustive-novelty claim. Identify the closest overlap and
   separate mature algebra from the proposed physical synthesis.
10. Retain counterexamples and failed hypotheses. If a failure undermines the
    core thesis, prefer abandon or a narrowly justified revision over cosmetic
    metric changes.
11. Produce an English academic manuscript draft in LaTeX/PDF if the pipeline
    reaches writeup, with prominent disclosure that AI tools were used as
    required by the workflow license. The paper must distinguish proved
    finite-dimensional results, conditional derivations, observed numerical
    evidence, and open conjectures.
12. Hard new proofs and scientifically decisive theoretical conflicts remain
    for top-level Codex/GPT Pro review. The coding worker may test them or find
    counterexamples but must not declare them solved from numerical evidence.
13. The repeated-block retention invariance is a no-prior (or jointly scaled
    prior) statement. With fixed finite $J_X$, repeated data change the
    data-to-prior weighting and generally change normalized retention.
14. Do not headline the result as a direct “deformation of the SOM/TSOM
    subspace” unless an explicit $A=G_ST_\chi$ pullback/composition theorem is
    proved and implemented. By default, SOM is a current-space computational
    substrate and the pose-confounding spectrum is a distinct map-tangent/data-
    space layer; show their composition without identifying their modes.
15. Document the 2D Green-function normalization, cell area, singular self-cell
    quadrature/regularization, incident-field convention, and contrast units.
    Centered finite differences validate derivatives of the implemented solver,
    not the physical fidelity of an undocumented discretization; include at
    least one refinement diagnostic and scope the claims accordingly.

## Source priority

The material below is distilled from:

- `Theory/SOM_SLAM_THEORY_CONTEXT.md`
- `Theory/Questions/Q1.md`
- `Theory/Questions/A1.md`
- `research/delegated/context_isolation/summary.md`
- `research/literature/TARGETED_PRIOR_ART_ADDENDUM.md`

If a compressed statement appears ambiguous, inspect those paths directly.
The targeted addendum identifies SAR autofocus and RF-SLAM mapping bounds as
mandatory neighboring work; do not rediscover their ingredients as novel.
