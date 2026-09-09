# Claim revision sheet for the paper draft

Purpose: exact edits needed to bring the draft's claims inside the executed
evidence. All changes preserve the negative results; none requires new code
except the two marked experiment additions (E-mode, E-runtime), which are
small and optional fallback triggers.

## 1. Title

From:

> Geometry-Lifted TriSpace SOM for Self-Calibrating Full-Wave Inverse
> Scattering under Unknown Array Geometry

To (preferred):

> Retained-Rank Pose Lifting in Subspace-Optimization Inversion: When
> Geometry-Error Equivalence Is Vacuous and What Truncation Hides (2D
> Helmholtz Study)

Reasons: "TriSpace" names an object that does not have the standing of the
current/data/state spaces (review Q3); "self-calibrating" without the
state-witness mechanism and without algorithmic gain over direct GN overstates
(Q5, Q6). The replacement names the two results that actually survived:
non-vacuity under restriction and the retained-rank freezing behavior.

## 2. Abstract

Required sentence replacements:

- Delete: "the full-wave state equation supplies a second consistency witness
  that distinguishes merely data-equivalent currents from physically
  admissible changes."
- Replace with: "The full-wave state equation supplies a state-consistency
  witness, which we keep as a pseudo-current diagnostic; in the tested
  regimes it does not discriminate good from bad pose estimates (negative
  result)."
- Delete: any benefit wording for the reduced solver. Replace with: "the
  reduced parameterization coincides with direct joint inversion in
  identifiable pose directions and freezes hidden ones."
- Qualify the rank transition: "for uniform full-circle arrays we observe
  that the retained-rank truncation hides exactly two pose directions at
  r=M-2 across M=8,12,16; the mechanism is conjectured to be azimuthal
  Fourier symmetry and is not proved."
- Replace "break a gauge" (if present) with "remove the point-source
  orientation rank deficiency".

## 3. Contribution list (max four, verbatim)

C1. A model-conditional split of the pose tangent into receiver-side
re-sampling and transmitter/illumination physical-current change under a
world-fixed grid, with the double-counting failure mode made explicit.

C2. A canonical minimum-norm lift of the receiver-side tangent into a
declared retained SOM basis plus the orthogonal retained-basis leakage
residual, with the unrestricted full-row-rank case shown vacuous (negative
control) and per-channel phase/gain identified as the canonical
residual-producing perturbation.

C3. A first-order data-side hiding condition for the retained-rank
parameterization, exercised in confounded and disambiguated regimes, and the
observation that uniform full-circle arrays freeze exactly two pose
directions at r=M-2.

C4. Negative results: the state-consistency witness is non-discriminating in
the tested regimes; exact and soft state-consistent variants fail to recover
frozen directions.

Everything else (SOM/TSOM machinery, joint estimation, blind calibration,
extension FWI, pseudoinverses) is cited as antecedent. The multi-transmitter
result is a sub-result of C1/C3, not a separate contribution.

## 4. Terminology replacements

| Current term | Replacement | Why |
|---|---|---|
| TriSpace | three-channel pose signature (Q, R, D) | V_P is not a third space |
| third space / V_P | pose-equivalent current subspace, explicitly `V_P = Range(U C_U) subseteq Range(U)`, dim <= p | it is subordinate to the retained current space |
| irreducible data residual | retained-basis leakage residual `r_SOM`; state that `r_S = 0` in the full-row-rank regime | "irreducible" is false w.r.t. the full current space |
| lossless | coincides with | equivalence, not gain |
| break the pose gauge | remove the source-symmetry rank deficiency | global SE(d) gauge untouched |
| TSOM | SOM-truncated retained basis (TSOM fold unverified) | no unverified fold presented as canonical |
| verified hiding condition | verified data half of the hiding condition | T_U half never binding |

## 5. Section-level edits

Section 5 ("State-consistent graph TriSpace"): rename to "Three-channel pose
signature"; keep Proposition 2 as algebra; state explicitly that no executed
case has `T_U h = 0`, so the state half is proved but not observed to bind.

Section 4 (lift): add the per-channel phase/gain calibration as the minimal
non-vacuous residual instance (math-sanity C3, residual ratio 1/sqrt(2) for a
single channel), and the structure-preservation negative (C2: rank-deficient
rigid motion still leaves zero residual).

Section 6 (algorithm): either add the runtime/DOF table (E-runtime below) or
delete every benefit sentence. Present the half-circle spurious minimum (pose
error 0.86, all seeds) as a headline limitation of local convergence.

Section 7 (numerics): retitle subsection 7.5 from "state disambiguation" to
"data-side hiding and the state-witness null"; add Table with hidden vs
visible direction residuals (`resid_vs_Sdata` 5.3e-16 vs 0.131) and their
`T_U` values (0.733 vs 0.674) so the non-discrimination is visible in one
place.

Section 2 (related work): add a dedicated Bellomo et al. (TAP 62(5):2450-2462,
2014) paragraph stating they already correct the incident field and
data-equation Green function via measured incident fields and multipolar
expansion; state the difference as joint pose estimation vs pre-calibration.
State that receiver-extension FWI is the closest genus and that the
remaining differentiator is the retained-basis leakage, not the residual
pair.

## 6. Two small experiment additions (fallback gates)

E-mode (for Q4): compute the azimuthal Fourier content of the two hidden
directions and of the two discarded singular modes at r=M-2, and report
alignment; perturb receiver angles slightly and an odd-M case to test
fragility of the s_{M-3}=s_{M-2} degeneracy.

E-runtime (for Q5): record wall-clock, iterations, forward/Jacobian
evaluations, and effective DOF for direct vs reduced at r=4 and r=M-2 across
the settings sweep.

If E-mode does not produce a provable Fourier mechanism (or a theorem), the
r=M-2 sentence becomes "observed in uniform full-circle arrays at three
aperture sizes". If E-runtime shows no advantage, delete all cost/benefit
language.

## 7. Acceptance gates before the revised draft is re-reviewed

1. Title/abstract contain none of: TriSpace, third space, state witness as a
   mechanism, algorithmic benefit, gauge breaking, global novelty.
2. The `T_U` null appears in the abstract and Section 7 as a first-class
   result.
3. Every numerical claim carries the qualifier "deterministic synthetic 2D
   scalar Helmholtz".
4. The three closest extension-FWI papers and Bellomo et al. are full-text
   checked for an equivalent restricted-lift or residual-pair construction;
   the result (found or not) is reported in the originality firewall.
5. Decision recorded: full-length revised paper, or letters-tier cautionary
   note if gates 4/E-mode fail.
