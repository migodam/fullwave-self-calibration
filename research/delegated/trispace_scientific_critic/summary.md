# Reviewer summary

Date: 2026-09-04. Independent reviewer pass over the corrected
geometry-lifted retained-range SOM study. Retrieval-bounded; no global
novelty claim; no experiment or parent-theory files edited.

## Verdict

**Major revise** (paper draft), not accept, not abandon. Fallback: if the two
gates below fail, reposition as a letters-tier cautionary identifiability
note.

## What survives hostile review

- The full-row-rank vacuity control is honest and correctly framed as a
  boundary, not a result.
- The retained-rank lift/residual decomposition (Prop. 1) and its
  orthogonality (`~1e-18`) are solidly verified, though textbook linear
  algebra.
- The receiver-re-sampling vs transmitter-physical-current split is a real,
  falsifiable modeling distinction (FD checks to `~1e-7`).
- The data-side hiding structure (confounded rank 15/16 -> 2 hidden DOF) is
  exercised with constructed directions.
- The negative state-witness results (E7/E7b/E11) are clean and complete:
  `T_U` is non-discriminating, and neither soft penalties nor exact variable
  projection recover frozen directions.

## What is overclaimed in the current framing

- State-consistency as a mechanism in title/abstract (falsified in every
  tested role).
- "Third space" for `V_P`, which is a pose-dimension subspace inside the
  retained current range.
- Algorithmic benefit from reduced == direct (equivalence only; no
  runtime/DOF evidence).
- `r=M-2` as a structural rule (uniform full circles only; tied to a machine
  degeneracy; no theorem or mode identification).
- Multi-transmitter results as "breaking a pose gauge" (they remove a
  point-source rank deficiency; the global SE(2) gauge is untouched).

## Where top-level Codex judgment is required

1. Final novelty adjudication for the restricted lift + retained-basis
   leakage after full-text checks of receiver-extension FWI - supported
   recombination vs new narrow formulation.
2. Whether to invest in the r=M-2 Fourier theorem, or keep it as an
   observation.
3. Whether Li-Lee-Bresler bilinear identifiability already subsumes
   Proposition 2 algebraically, shrinking the contribution to the physical
   split.
4. Full-length revised paper vs letters-tier negative/identifiability note.

Artifacts: `review.md` (Q1-Q10 answered), `claim_revision.md` (exact edits),
`hard_questions_zh.md` (Chinese questions for GPT Pro).
