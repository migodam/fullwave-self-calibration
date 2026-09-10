# Material Attribution under Calibration Uncertainty: Independent-Model Diagnostics and Conditional Reference Requirements

## Integrated completion note

The automatic verification has ended without native writeup. Its bounded checks are incorporated below; the candidate-selection claim remains withdrawn. The 53-entry ledger includes input consistency checks, not 53 independent theorems. Corrected worker formulas in Appendix A override the raw report.

### Fixed-loss Born scope

For u=q ell with ell=(0.03,0.05), materials (1.9,2.5) and (2.2,3.0) yield identical fixed-geometry Born data under gains 1 and (30+i)/(40+i). The latter modulus is 0.750182, within the allowed annulus. The q'=45 example instead needs common gain rescaling. Constitutive anchoring is therefore conditional, not a universal removal of ambiguity.

### Verification incorporated into the reference argument

Independent computation reproduces the supplied binary error 0.443423189256. Completed worker calculations yield weaker amplitude-derivative bounds than supplied, without thereby refuting them. No unverified endpoint-infimum or physical-readout claim is promoted into a theorem. Scalar reference-GLS equivalence, gain constraints and complete residual derivatives are implementation requirements, not a novel method. Two raw written derivative errors are corrected in Appendix A (Q5_THEORY_FINAL.md).

The native evaluator's abandon decision concerns treating this bounded verification as a paper contribution; it does not erase the valid checks. The remaining finite-domain, physical and novelty gates are unchanged.


Working manuscript — scientifically incomplete; not submission-ready.

## Abstract

Accurate receiver calibration and small coherent-data residuals need not imply correct quantitative material attribution. We independently reconstruct a two-region electromagnetic inversion experiment with unknown receiver translation and shared complex gain. A weighted noisy reference improves development recovery, but independent forward-model discrepancy remains consequential: joint recovery succeeds in six of twelve independent cases with reference, versus one without. A finite-competitor single-observation selector does not improve these outcomes and is withdrawn. For a separate known-geometry modal class, we derive an implementation-error budget for a conditional material inversion bound. Neither uniform finite separation for the imaging class nor a new electromagnetic impossibility theorem is established. The present results delimit the missing evidence required for a reliable calibration-assisted material claim.

## I. Motivation and scope

Coherent measurements encode material-dependent relative field structure, but unknown electronics and geometry can imitate some changes. Discarding phase avoids one nuisance at the cost of information; retaining it is useful only if correct material attribution survives calibration and model uncertainty. Our scope is two known isotropic supports, not unknown shape, SLAM or a SOM-specific advantage. Nonlinear scattering alone is not a recovery guarantee.

## II. Observation and inverse problem

Let y=g F(epsilon_1,epsilon_2,x)+n, optionally accompanied by z=g+eta. Known losses are 0.03 and 0.05; their constitutive anchoring must not be mistaken for full-wave information. For each bounded material/geometry candidate we minimize the complete weighted residual over an annulus-constrained complex gain. This is standard variable projection, not an original algorithm. Reference uncertainty enters the objective explicitly.

## III. Finite attribution and intervention

Local profiled derivatives diagnose infinitesimal directions; only finite competing worlds address the target material tolerance. We evaluated a finite-competitor action rule using one noisy reference or one shared-gain field scalar. Its score is a heuristic because it lacks a covering lower bound. The independent development results do not support this rule as a superior acquisition method.

## IV. Independent-model development evidence

The protocol and full counts appear in Q5_EXPERIMENT_REPORT.md. Same-principle data yield 11/12 successes without reference and 12/12 with it. Independent DDA data yield 1/12 and 6/12. Selected, random and fixed added EM scalars each yield 1/12; the two-action policy yields 3/12. Thus the additional reference can improve attribution in this development setting, while neither it nor the proposed selector establishes robust recovery. No formal test-set or population-level superiority is claimed.

## V. Conditional reference precision

For separately read known-geometry modal amplitudes h_i with certified derivative lower bounds m_i, a noisy real gain-modulus reference yields material error at most (B_i+H_i B_a)/((a_min-B_a)m_i). The bound makes model bias and reference drift explicit. It is a restricted inverse-stability corollary; its physical readout assumptions and novelty remain to be established. Appendix A is Q5_THEORY_FINAL.md, including the short proof.

## VI. Limitations and decision

The forward discrepancy is not a certified continuum error bound. Full mechanism ablations, electronics-breakdown cases, nearest-work full-text comparison, uniform finite separation and final frozen evaluation are missing. We reject the current selector-superiority claim, not all possible calibration-assisted methods. A mature TAP submission requires closing a physically nontrivial boundary or independent recovery route, rather than promoting these development diagnostics into a completed method paper.

## Appendices and evidence

Theory: Q5_THEORY_FINAL.md. Full development definitions: docs/. Raw results and corrected metrics: results/. Criticism and repair: Q5_REVIEW_AND_REPAIR.md. Bibliographic evidence status: Q5_LITERATURE_LEDGER.md. No unchecked reference is used as support for priority.
