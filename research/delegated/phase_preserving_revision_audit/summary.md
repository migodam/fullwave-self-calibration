# Phase-preserving revision audit — summary

**Verdict: PASS (bounded), with conditions.**

The draft survives the audit as a correctly-hedged theory foundation: no fatal
type/dimension error, no false phaseless-SOM statement, no unjustified novelty or
TGRS/TAP readiness claim. The central method claim (PRASC-SOM beats direct joint
inversion) is honestly labeled open, and the audit does not change that: it is
open, not proven. All defects found are refinements.

## Findings

1. Core math spot-checks pass: state differentiation, realification, lift/leakage
   and full-rank vacuity, Fisher DPI, intensity-Jacobian phase cancellation,
   principal-angle form of ρ, Loewner monotonicity, dimension obstruction,
   stacked shared-map quotient Gram, and ρ(r) non-monotonicity.
2. Eq (4) under-specifies the provenance of `A_χ` across the reduced and
   independent-current formulations; in the independent-current case χ enters data
   only through the state constraint, so the "conservative" claim (line 257)
   should pin down what `A_χ` is to avoid the double-counting the draft itself
   warns about.
3. Theorem 1's assumptions are adequate for the fixed linearization, but the
   statement should explicitly exclude reduced formulations with rank-dependent
   total derivatives `A_χ, B`.
4. Phaseless SOM: zero false statements; the draft explicitly rejects the
   "no current-space constraint" strawman (checked against memory, not a fresh
   primary read).
5. ρ/K_eff vs B_vis/J_x are correctly separated; only the "two sides of the same
   local confounding geometry" phrasing risks implying an unproven duality — one
   sentence pinning the joint-Gram Schur-block relation (or labeling P2.3 open)
   would remove the risk.
6. GPT Pro prompt asks the right decisive question (P4.3, non-avoidable, with a
   downgrade fail-safe), but it is buried as the third subpart of the second
   priority, and much of P1.1 duplicates the already-proved Prop 2. The
   "sufficient to support TGRS/TAP" framing (lines 6, 29) invites confirmation
   bias despite the guardrails.
7. Failure tests against direct joint inversion are substantially present (E4
   baselines, oracles, failure rate, matched cost/accuracy, inherited negative
   control), but there is no pre-registered superiority decision rule, no
   cycle-skipping failure test of the frequency gate (half-circle ≈0.86
   pseudo-minimum is the natural control), and no expected-loss condition
   (non-pose phase corruption).
8. Minor cross-document inconsistency: the paper's gauge-quotiented `p_g` vs the
   review's unquotiented `p = 3` with separate SE(2) handling; harmonize.

## Uncertainty

Audit was bounded to the three files; numeric evidence was read as reported, not
re-derived from raw loop artifacts or figures; phaseless-SOM attribution and the
reference list rely on knowledge/memory, not fresh primary-source or DOI checks;
no nonlinear re-run performed.

## Codex intervention needed

Yes — five small refinements, none fatal to the direction: (1) pin `A_χ`
provenance in eq (4); (2) add the Theorem 1 scope sentence; (3) state or mark
open the §4.4 duality; (4) re-anchor/re-rank the GPT Pro prompt with P4.3 as a
top-level gate; (5) pre-register an E4 superiority decision rule plus
cycle-skipping and expected-loss controls. The method claim itself must stay
labeled open.

## Parent disposition

Applied after the audit on 2026-09-05: the augmented-tangent provenance and
state-feasible nesting caveat are explicit; the joint gauge is no longer split
without an anchor; Section 4.4 now states the block-Schur relation conditionally;
the GPT Pro prompt starts with a SOM-specific pass/fail gate; and E4 now requires a
predeclared matched-budget decision rule, a limited-aperture cycle-skipping control,
and non-pose phase-corruption expected-loss controls. The central performance claim
remains open.
