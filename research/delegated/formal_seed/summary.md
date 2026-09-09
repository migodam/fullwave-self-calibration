# Formal Seed — Summary

## Findings

1. The source theory is a coherent, self-consistent skeleton: 2D scalar
   Helmholtz contrast-source forward model, local linearization
   $\delta y=A\delta\chi+B\delta X+n$, and the operator family
   $K_{\mathrm{IS}},K_{\mathrm{SLAM}},K_{\mathrm{eff}},L_X$ that converts
   map--pose tangent geometry into absolute and relative information spectra.
2. The algebraic spine (kernel/rank identities, $\rho_i=\sin^2\theta_i$,
   low-rank loss, interlacing, prior ordering, gauge, Born empty-background,
   multi-frequency shared compensation) is preserved verbatim in spirit in
   `seed.md` with accepted/conditional/open/non-claim labels.
3. The Born empty-background caveat is retained as non-negotiable: first-order
   $B=0$ at $\chi_0=0$ must not be read as "pose harmless"; the second-order
   bilinear term governs geometry mismatch there.
4. Whitening + realification, gauge support discipline, and the ban on
   semantic drift between current-space $G_S$/$V_S^\pm$ and the map Jacobian
   $A$/$K_{\mathrm{eff}}$ are encoded as guard rules in the seed.
5. The proposed checks organize cleanly into exactly five experiment families
   (Jacobian consistency; algebraic identities; gauge/Born degeneration;
   frequency and trajectory geometry; sensitivity/robustness/rank events), each
   with hypothesis, controls, metrics, pass/fail rules, and non-claims.
6. Open questions split into four tiers: fatal core-validity, hard
   theory/proof (GPT Pro later), design/construction, and pipeline-solvable.
7. No global novelty is asserted anywhere; the closest-prior-work mapping in the
   sources remains unverified and is scheduled for a fresh systematic check.

## Artifacts (this package)

- `research/delegated/formal_seed/seed.md` — theory-first, falsification-first,
  novelty-first workshop seed.
- `research/delegated/formal_seed/experiments.md` — five experiment families.
- `research/delegated/formal_seed/open_questions.md` — four-tier question
  register.
- `research/delegated/formal_seed/summary.md` — this file.

Absolute root: `/Volumes/migodam's-external-brain/Research/Inv_SLAM/`.

## Uncertainties

- Two-fold SOM $V_D^\pm$ exact definition per a specific Chen source remains
  unresolved.
- Continuity/rank-event theory, prior-weighted spectrum bounds, multi-frequency
  transversality, and Fourier surrogate errors remain unproved open claims.
- Numerical verification has not been run; nothing here upgrades open claims to
  accepted results.
- Novelty/prior-art boundaries are unverified and explicitly not claimed.

## Codex intervention

Not needed for mechanical reading or organization of this package. Needed for:
(1) adjudicating any negative answer to the (a)-tier core-validity questions;
(2) the final novelty judgment after a fresh systematic search; (3) deciding
whether later experiment failures falsify core claims rather than refine
conditions.
