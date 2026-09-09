# Skeptical review of the corrected TriSpace SOM study

Date: 2026-09-04. Independent, hostile-but-fair reviewer pass over the
corrected "geometry-lifted retained-range SOM" study. Scope is the *paper
draft*, not the experiment harness. Every novelty statement here remains
retrieval-bounded; this review makes no global novelty claim and does not edit
experiment or parent-theory files.

Inputs read in full: `research/trispace_self_calibration/THEORY_SEED.md`,
`ORIGINALITY_AUDIT_ZH.md`, `MANUSCRIPT_BLUEPRINT.md`,
`research/delegated/trispace_prior_art/{summary,novelty_overlap,claim_ledger,screened_papers,verified_references}.md`,
`research/delegated/trispace_math_sanity/{summary,counterexamples,derivation}.md`,
`current_idea.json`, `round_00/01/02_outcome.json`,
`consolidated_round3_tables.md`, and the final `experiment_report.json`.
The Bellomo et al. anchor was also verified externally: "An Improved Antenna
Calibration Methodology for Microwave Diffraction Tomography in Limited-Aspect
Configurations", IEEE TAP 62(5):2450-2462, 2014, DOI `10.1109/TAP.2014.2308534`.

## 0. Verdict in one paragraph

**Major revise, not accept, not abandon.** The package has real virtues that
survive hostile scrutiny: an honest full-row-rank vacuity control, a clean
orthogonal lift/residual decomposition, derivative machinery verified against
finite differences, a reproducible retained-rank hiding structure, and
well-kept negative results. But the paper as currently titled and abstracted
overclaims on at least four axes (state-consistency as a mechanism, "third
space", algorithmic benefit, and the r=M-2 "rule"). Its most distinctive
theoretical object, the state witness `T_U`, is a null result in every tested
role. If the two gate items in Q4/Q5 (a theorem for r=M-2 and runtime/DOF
evidence) cannot be produced, the paper should be repositioned as a
cautionary finite-dimensional identifiability note, not a self-calibration
method paper. The internal loop verdicts ("revise" then "lock, lock") judge
*experiment sufficiency*; this review judges *what the evidence licenses the
draft to claim*, and the two disagree.

## 1. Which hypothesis components survived, were narrowed, or were falsified?

| Component | Fate | Evidence |
|---|---|---|
| Unrestricted "geometry error is an equivalent current" | Falsified as a claim; retained as a vacuity control | E1: rank `G_S`=8/8, unrestricted lift residual `5.5e-16` |
| Minimum-norm retained-rank lift `Q_U` + orthogonal residual `R_U` (Prop. 1) | Survived (textbook algebra; see Q2) | r=4 residual ratio 0.784; `Q^H R` orthogonality `8.9e-18` |
| Receiver re-sampling vs transmitter physical-current split (`H_S` vs `M^{-1}H_D`) | Survived as a model-conditional modeling claim | E3: FD errors O(eps^2), `B_l` split identity to `2.6e-18` |
| Pose-dependent spectral drift, rank events, soft filter | Survived | E2: exact degeneracy gap `1.4e-15`; projector jump 1.424 (hard) vs 0.127 (soft) |
| Exact local hiding condition `R_U h=0` and `T_U h=0` (Prop. 2) | Narrowed: only the data half was ever exercised; the state half was never binding | E4: "hidden" direction has `resid_vs_Sdata=5.3e-16` but `resid_vs_Sstate=0.947`, `T_U=0.733` |
| `T_U` as a discriminating consistency witness | Falsified in the tested regimes | E4: `T_U` 0.67-0.80 on hidden and visible directions; E5: 0.81-0.89 across estimators |
| State-consistency as an algorithm ingredient (penalty or variable projection) | Falsified | E7b soft: pose 0.085/map 0.445 at lambda=0, collapse at lambda>=1e-2; E11 soft VP: map 0.838 |
| Reduced parameterization lossless for visible pose DOF | Survived only as equivalence, not benefit | reduced_r4 == direct at solver precision in every setting |
| Freezing of hidden pose DOF | Survived, but it is a demonstrated loss (direct recovers what reduced freezes) | M8 r6: 0.0442 vs direct 0.0158; M12 r10 hidden component 0.072 vs direct 0.00725 |
| r=M-2 hiding transition | Observed, not established (see Q4) | M=8/12/16 uniform full circles; fails for half/three-quarter arcs |
| Multi-transmitter/directional sources remove the point-source theta null | Survived as degeneracy removal, not gauge breaking (see Q7) | L=1 monopole `B_t` rank 2; L=2/L=3 rank 3; 240/240 E6b runs |

The largest silent narrowing: the blueprint's central sentence promises a
state witness that "distinguishes merely data-equivalent currents from
physically admissible changes." The executed study never demonstrated that
distinction in either direction. In E4 the constructed hidden direction is
data-hidden but state-visible (`resid_vs_Sstate=0.947`), i.e. `T_U != 0`, so
it is not a joint hiding direction at all; and `T_U` was never close to zero
for any physical direction, so the state half of Proposition 2 was never
satisfied, let alone shown to be the binding constraint. The audit's own
failure-mode #2 ("T_U is nearly always absorbed by the real contrast tangent
and adds nothing beyond data-only") appears to be exactly what happened.

## 2. Are the minimum-norm lift and orthogonal residual merely textbook linear algebra?

Yes, and the paper should say so. Proposition 1 is the Moore-Penrose
decomposition of `H_S` onto `Ran(K_U)` and its orthogonal complement, plus the
standard minimum-norm property of the pseudoinverse; Corollary 2 is the SVD
formula with the `||Q|| <= ||H_S||/sigma_r(S)` amplification bound. The
vacuity corollary (full row rank implies `Ran(K_U)=C^m`) is elementary. The
`T_U` construction is a projection onto `N_D^perp`, and Proposition 2 is the
joint-Jacobian solvability condition split into two blocks - finite-dimensional
linear algebra throughout. Claiming any of this as new mathematics would be
wrong.

What remains scientifically nontrivial, and should carry the paper:

1. The model-conditional physical split `B h = H_S h + S M^{-1} H_D h` under a
   world-fixed grid with `D_x D = 0`. The claim "lift only `H_S`, keep `H_D`
   in the state equation, or you double-count pose" is a falsifiable modeling
   statement with stated hypotheses (world-fixed domain, homogeneous
   background, transmitter off the rotation center), not linear algebra.
2. The three-residual taxonomy `r_S` (out of `Ran G_S`), `r_SOM` (out of the
   retained-basis image `Ran(G_S P_T)`), `r_Schur` (map-pose confusability).
   In the finite-m full-row-rank regime `r_S = 0` identically and the only
   non-degenerate object is `r_SOM`; the math-sanity counterexamples C1/C2/C3
   show the naive "rank-deficient implies residual" intuition is false and
   per-channel phase/gain is the canonical residual-producing perturbation.
   That correct demarcation has value even if the algebra does not.
3. The honest negative statement: the vacuous unrestricted case, the failure
   of the state witness, and the retained-rank freezing behavior are useful
   to the SOM/TSOM community, which routinely treats `G_S` as known.

Without (1) and (2), the paper collapses to an algebra exercise.

## 3. Does `V_P = Range(Q_U)` deserve the name "third space"?

No. By construction `V_P^{(U)} = Range(U C_U) subseteq Range(U)` with
`dim V_P <= p` (pose dimension 2-3 here). It is a low-dimensional,
scene-dependent tangent object inside the retained current space; it has none
of the standing of the current/data/state spaces or the `V_S^+-/V_D^+-`
folds. The "TriSpace" name invites a false symmetry among three comparable
spaces, and the rejected "eight-cell" `2^3` intersection picture (math-sanity
C8) is the natural misreading it produces.

The graph formulation is strictly more honest and should be adopted: the
object is a single tangent map

```
Gamma_P^{(U)} = { (Q_U h, R_U h, D_U h) : h in R^p }
```

whose three coordinates live in three different ambient spaces (current,
data, state) and are coupled, not three commuting projectors. Recommended
language: "three-channel pose signature (Q, R, D)"; the third channel is a
state-space defect, not a "space" of the current-space reduction. Note also
that `Q_U` itself depends on the scene (`j`, contrast), whitening, current
norm, retained basis, and regularization, so it is not a standing spectral
fold. Retitle the paper accordingly (see claim_revision.md).

## 4. Is `r = M-2` a structural phenomenon or an artifact?

Presently: an artifact-shaped observation tightly coupled to uniform
full-circle Fourier/Hankel symmetry, not yet a publishable structural
phenomenon. The evidence cuts both ways:

- For uniform full circles it is clean and reproducible: M=8 (r6), M=12
  (r10), M=16 (r14), hidden_rank=2, n_vis=1, with `r=M-1,M` collapsing to an
  "empty lift" (`B_red ~ 1e-19`). Across k=8/16 and one alternative scene the
  structure persists.
- But it is not aperture-general: the half-circle and three-quarter arcs
  show r6 with n_vis=3/hid=0, and the half-circle `known_alpha` baseline hits
  a spurious minimum (pose error 0.86) on all three seeds. The rule therefore
  carries a "uniform full circle" qualifier the current abstract does not
  state.
- The transition sits on an exact machine degeneracy (E9: s5-s6 pair
  degenerate to 1e-19 while neighboring gaps are 1e-5 to 1e-4; E11: the hard
  r=6 cut lands exactly on the degenerate s6=s7 pair and produces a
  gauge/numerical null). The rank cut itself is therefore ill-posed at the
  transition: a small physical asymmetry that splits the pair would smear the
  "rule".

The plausible mechanism is symmetry, and the paper should say so: a uniform
circle makes the sensing operator approximately circulant/block-diagonal in
azimuthal Fourier modes, giving paired singular values for `+-k` modes; an
infinitesimal rigid translation/rotation of a uniform array perturbs the rows
primarily through the azimuthal derivative, i.e. through fixed low harmonic
content, so retaining M-2 modes discards exactly the pair that carries pose
and yields a rank-2 hidden subspace. That is a hypothesis, not a theorem.

What is needed to upgrade it:

1. A theorem with explicit hypotheses (uniform full circle, world-fixed grid,
   stated whitening/realification, low-dim real contrast basis) that
   identifies the discarded harmonic pair and proves
   `rank([K_U, contrast-image]) = M-2` at `r=M-2`, with the failure mode under
   aperture truncation stated.
2. A mode-identification experiment: decompose the two hidden directions and
   the two discarded singular modes in azimuthal Fourier components and show
   they align. This is cheap and decisive.
3. Robustness probes: slight angular jitter, odd M, off-center rotation axis,
   finite split of the degenerate pair, to show the rule is not a measure-zero
   coincidence of exact symmetry.

Until (1)-(3) exist, the honest claim is "a dimension-count observation
reproduced at three uniform full-circle apertures", not "a structural
transition at r=M-2".

## 5. Does reduced == direct demonstrate an algorithmic benefit?

No. Equality to solver precision in visible regimes shows only that the
reparameterization spans the same solution set as direct joint GN on the
visible subspace - algebraic losslessness of an in-sample reparameterization
and, implicitly, implementation correctness. The reduced method reuses the
same forward model and derivatives, so it removes only the few discarded
current directions from the optimization variables while the dominant cost
(forward solves and Jacobians) is unchanged. In the hidden regime the reduced
method is strictly worse: it freezes directions direct GN recovers (M8 r6:
0.0442 vs 0.0158; M12 r10 hidden component 0.072 vs direct 0.00725).

The consolidated tables contain pose/map errors and residuals but no
wall-clock times, iteration counts, Jacobian/forward evaluations, or effective
DOF columns, even though the planned metrics list promised runtime and
effective variables. Any algorithmic-benefit sentence is currently
unsupported. To claim a benefit, the paper needs one of: (a) runtime/DOF
tables showing a real reduction with equal accuracy; (b) a regularization case
where the restricted parameterization demonstrably improves basin or noise
behavior; or (c) delete all benefit language and say "coincides with", not
"lossless" or "reduced cost". Recommend (c) as the default.

## 6. How damaging is the negative `T_U`?

Damaging to the paper as currently framed, recoverable with honest
repositioning. The state-consistency witness was the single most distinctive
theoretical object (audit candidate C3, "strongest narrow candidate"), and
every algorithmic use of it failed:

- E4: `T_U` 0.67-0.80 on both hidden and visible directions - no separation.
- E5: 0.81-0.89 across all estimators, good and bad - no scoring power.
- E7: direct `T_U` penalty (`lambda_T=1e-3`) lowers `T_U` to 0.59-0.69 but
  gives pose error 0.52.
- E7b: soft state-consistent variants either land on a meaningless optimum
  (lambda=0) or collapse contrast (map ~1.0 for lambda>=1e-2).
- E11: exact variable projection stalls on the degenerate SVD cut or
  converges to a distant local minimum; the soft version is model-biased
  (map 0.838, discarded-current fraction 0.87), confirming the failure is not
  an artifact of the auxiliary coefficient `c` or a scalar penalty.

So the title/abstract sentence "the state equation supplies a second
consistency witness that distinguishes..." must go; as a mechanism the
witness is falsified in the tested model. But it should not be deleted from
the paper. The right disposition is exactly the one the current idea.json has
already adopted: keep Proposition 2 as exact finite-dimensional algebra, keep
`T_U` as a pseudo-current diagnostic with the null result reported in the
abstract, and add the mechanistic explanation as a tested hypothesis (the
real contrast tangent absorbs `D_U`, so data-side rank alone determines
hiding in this regime; audit failure-mode #2). Two caveats to state: the null
is regime-specific (near-resonance `M^{-1}` amplification could make the
state channel bind elsewhere), and the state half was never shown binding in
any executed case, so "we verified the hiding condition" must read "we
verified its data half".

## 7. Do multi-transmitter/phased-dipole results break a pose gauge?

No - they remove a source-symmetry rank deficiency, which is different. For a
scalar point source, rotation of the single transmitter about its own location
is a null direction of the data: `B_t` has rank 2 with the theta column
structurally collinear with tx/ty for all phi0 (verified at residual 1.7e-18).
A two-element array or a directional element makes the source's orientation a
genuine data DOF: rank 3, theta residual 8.9e-4 (L=2) / 1.2e-3 (L=3), with
SNR-scaled recovery in E6b. This is parameterization physics, not gauge
theory: one cannot "break a gauge" by adding more measurements of the same
symmetric source; one changes the physical model so the formerly null
coordinate becomes observable.

The global SE(2) map-pose gauge (math-sanity C5: `B xi_X = -A xi_chi`,
pose/map interchangeable without an anchor) is untouched by all of this, and
anchors/priors remain necessary. Also temper "escaped": the L=1 phased dipole
reaches rank 3 only weakly - third realified (column-normalized) singular
value 0.0583 against leading 1.41, i.e. roughly 25x worse conditioning for
theta than for the translations at equal noise. Cleanest admissible claim:
"local identifiability of transmitter orientation transitions from null
(isotropic point source) to positive but SNR-limited (multi-element or
directional source); the global pose-map gauge is not affected."

## 8. Proven FD claims vs executed evidence vs interpretation vs open theory

Proven finite-dimensional (subject to whitened/realified inner products,
closed range, stated rank conditions): Prop. 1 (lift-residual identity,
min-norm, rank <= p); Cor. 1 (exact representability iff
`Range(H_S) subseteq Range(K_U)`; full-row-rank vacuity); Cor. 2 (SVD formula
and 1/sigma_r amplification); Prop. 2 (joint hiding iff `R_U h = T_U h = 0`);
the world-fixed pose-tangent split identity; the bridge statement in the
math-sanity derivation. The proofs are short and correct, but the standard
linear algebra content should be labeled as such (Q2).

Executed numerical evidence (deterministic 2D scalar Helmholtz, synthetic):
E1-E11 as summarized in Q1, plus E8 (76/76 replicates), E9 (N=16/24/32,
36/36), E10 (settings sweep), E6b (240/240 multi-Tx runs), and the FD checks
at 1e-7 to 1e-9 relative error. The derivative-chain and decomposition
identities are very solidly verified; the nonlinear results are
initialization/noise-conditioned but replicated.

Plausible interpretation, not proven: the Fourier/circulant mechanism for
r=M-2 (Q4); the two hidden DOF as a translation/rotation pair; contrast-tangent
absorption as the explanation of the `T_U` null; multi-element orientation as
the explanation of rank-3 theta.

Open theory: continuum non-closed-range regularization semantics; metric/
whitening dependence and canonicity of `Q_U`; pose-dependent basis alignment,
rank-event crossing, and global convergence; scene-dependence and possible
amortization of `V_P`; chart dependence under moving grids; the Born/full-wave
separation as chi->0; the unverified TSOM domain fold (so all experiments are
SOM-truncated, not TSOM); whether the state channel binds near resonance or
for richer contrasts; the half-circle local-minimum structure.

## 9. Comparison against the five families; narrowest contribution

- Source/receiver-extension FWI (Huang-Nammour-Symes; Metivier-Brossier;
  da Silva): closest genus. Adding geometry DOFs to absorb misfit and
  estimating corrected geometry is antecedent and must be cited as such. The
  remaining deltas are (i) the contrast-source/current-space substrate, (ii)
  the restricted retained-basis leakage `r_SOM` rather than total misfit, and
  (iii) the explicit physical split. With `T_U` null, the "residual pair"
  differentiator largely collapses to (ii); the honest residual claim is
  data-side only.
- Bellomo et al., IEEE TAP 62(5):2450-2462, 2014 (verified): measured
  incident fields + multipolar expansion correct both the incident field and
  the data-equation Green function, with phase-center treatment, i.e.
  "calibrate G, then invert" as a pre-calibration stage. The current work's
  difference is joint low-dimensional pose estimation from scattered data and
  the identifiability structure in current coordinates - a difference of
  problem formulation, not of the calibration idea. The related-work section
  must draw this boundary explicitly and must not imply Bellomo et al. did
  not already put measured antenna behavior into the Green operator.
- Classical SOM/TSOM: the substrate itself, fixed-geometry. The delta is only
  the pose-dependence of the retained basis and the lift; and because the
  TSOM domain fold is unverified, the paper's "TSOM" is at present a
  truncated SOM basis.
- Blind calibration / bilinear identifiability (Li-Lee-Bresler; radar/SAR
  autofocus): the generic joint-unknowns identifiability-up-to-gauge theory
  likely subsumes the algebra of Proposition 2 as a special joint Jacobian
  rank test. The paper must argue what the Helmholtz specialization adds
  beyond that generic certificate - the physical split of the pose tangent
  and the retained-basis leakage are the only candidates.
- Joint transmitter localization (Karthik-Ghosh DNN): problem-statement
  antecedent in inverse scattering; the current work is model-based and
  analytic, but "joint contrast + transmitter pose estimation" itself is not
  novel.

Narrowest retrieval-bounded contribution the evidence supports: "a first-order
identifiability certificate for antenna-pose absorption in a retained-rank SOM
parameterization of a 2D scalar Helmholtz model: the receiver-side pose
tangent splits into a minimum-norm lift inside the retained current basis
plus an orthogonal retained-basis leakage residual, the unrestricted
full-row-rank case is vacuous, and - in uniform full-circle arrays - the
retained-rank truncation freezes exactly two pose directions; the full-wave
state equation does not add discriminative power in the tested regimes
(negative result)."

## 10. Recommendation

Major revise, with exact edits in `claim_revision.md`. Summary of the
required changes:

1. Retitle and re-abstract; remove "TriSpace", "third space", any
   state-consistency mechanism claim, and any algorithmic-benefit language.
2. Demote `T_U` to a diagnostic with the null result in the abstract; narrow
   "hiding condition verified" to its data half.
3. Replace "break a gauge" with "remove a source-symmetry rank deficiency";
   report the weak theta conditioning of the L=1 dipole.
4. Downgrade r=M-2 to an observation with the Fourier mechanism stated as a
   hypothesis, or produce the theorem + mode-identification experiment of Q4.
5. Add the runtime/DOF table or delete benefit claims; replace "lossless"
   with "coincides with".
6. Promote the half-circle spurious minimum (0.86) to a headline limitation.
7. Add the per-channel phase/gain case (math-sanity C3) as the canonical
   nonzero-residual instance and connect it to Bellomo/phase-center
   discussion.
8. State the boundary against receiver-extension FWI, Bellomo et al., and
   generic bilinear identifiability in the related-work section, and
   full-text check the three closest extension-FWI papers before submission.
9. Mark all results "synthetic 2D scalar Helmholtz, finite-dimensional";
   label every TSOM usage as unverified.
10. Fallback: if gates (4)/(5) fail, target a letters-tier cautionary
    identifiability note instead of a method paper.

Where top-level Codex judgment is required (not decidable by this file-level
review alone):

- Whether the restricted lift + leakage residual is a supported recombination
  versus a genuinely new narrow formulation (final novelty adjudication after
  full-text checks of receiver-extension FWI).
- Whether the Fourier mechanism of Q4 is worth the theorem effort, or the
  paper should stop at "observation".
- Whether the generic bilinear-identifiability literature (Li-Lee-Bresler)
  already subsumes Proposition 2, which would shrink the contribution to the
  physical split alone.
- Whether a null-result state-witness paper is worth full-length form or
  should be a letters/short note.
