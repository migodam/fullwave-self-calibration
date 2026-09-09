# Geometry-Lifted TriSpace SOM for Self-Calibrating Full-Wave Inverse Scattering

Date: 2026-09-04

Purpose: a theory-first, falsification-first seed for the corrected research
direction. This file does not claim global novelty. It distinguishes proved
finite-dimensional algebra, model-conditional physics, numerical hypotheses,
and open problems.

> **Post-experiment note.** Sections 1--8 preserve the theory seed supplied to
> the autonomous research loop. Section 9 is the parent-level adjudication
> after three development rounds and supersedes every unsupported algorithmic
> or state-witness hypothesis above.

## 0. Status vocabulary

- **[proved-fd]**: elementary finite-dimensional result proved below under the
  stated inner products and rank conditions.
- **[model-conditional]**: exact for the stated discretized full-wave model,
  but not automatically for every coordinate convention or hardware model.
- **[hypothesis-to-test]**: a scientific or algorithmic claim requiring the
  specified experiment.
- **[retrieval-bounded candidate contribution]**: not found in the screened
  record in this exact combination; this is not a global novelty claim.
- **[open]**: not established and must not appear as a theorem or accepted
  conclusion.

All complex operators below must be whitened and realified before ranks,
orthogonal projections, pseudoinverses, and pose identifiability are evaluated.
Pose and contrast perturbations are real. This avoids silently replacing a
real nuisance space by a larger complex nuisance space.

## 1. Correct research object

The target is not generic inverse-scattering observability. It is a
self-calibrating extension of SOM/TSOM in which the antenna geometry and hence
the Green operators used to define the SOM coordinates are uncertain:

> Can the data effect of an unknown sensing Green operator be lifted
> canonically into a retained SOM/TSOM current space, while retaining the
> part that cannot be represented by a current and the state-equation defect
> created by that representation, so that array geometry and material contrast
> can be estimated together?

The phrase “geometry error is an equivalent current” is not yet a contribution.
With a full-row-rank unrestricted current-to-data matrix, every data perturbation
has an equivalent current. The useful object must therefore be restricted,
canonical, stable, and coupled back to the full-wave state equation.

## 2. Full-wave independent-current model

At a nominal pose `x`, write

\[
y=S(x)j,\qquad
\Phi(j,\chi,x)=j-D_\chi\{e(x)+D(x)j\}=0,
\]

where `S=G_S` maps induced current to measurements, `D=G_D` propagates current
inside the imaging domain, and `e=E_inc`. Define

\[
M=I-D_\chi D,\qquad E_{\rm tot}=e+Dj.
\]

For a real pose perturbation `h`, an independent current perturbation `delta j`,
and a real contrast perturbation `delta chi`, the linearization is

\[
\delta y=S\,\delta j+H_S h,
\tag{2.1}
\]

\[
M\,\delta j-D_{E_{\rm tot}}\delta\chi-H_Dh=0,
\tag{2.2}
\]

with

\[
H_Sh=(D_xS[h])j,
\qquad
H_Dh=D_\chi\{D_xe[h]+D_xD[h]j\}.
\tag{2.3}
\]

These equations separate two mechanisms that must not be conflated:

1. `H_S h` is a change in how a fixed current is sampled. A current that
   reproduces it is generally a data-equivalent or pseudo-current.
2. `M^{-1}H_Dh`, when it exists, is the actual induced-current change caused
   by a changed illumination or internal propagator at fixed contrast.

**[model-conditional]** On a world-fixed domain grid in a homogeneous fixed
background, moving external antennas changes `S` and `e`, but not the
domain-to-domain propagator `D`; hence `D_xD=0`. A moving/body-fixed domain,
geometry-dependent boundary, or changing background can make `D_xD` nonzero.
The manuscript must state which convention is used instead of writing
`G_D(x)` universally.

If the contrast is held fixed and `M` is invertible, the usual map-eliminated
pose Jacobian is

\[
B h=H_Sh+S M^{-1}H_Dh.
\tag{2.4}
\]

Lifting all of `B h` through `S^\dagger` hides the distinction between the
physical current response and receiver-side re-sampling. The proposed
independent-current construction therefore lifts `H_S` only and keeps `H_D`
in the state equation.

## 3. Canonical SOM-restricted geometry-to-current lift

Let `U` have orthonormal columns spanning a declared retained current space.
For a Phase-I SOM experiment, `U` is a truncated right-singular basis of `S`.
For a verified TSOM implementation it may include the retained
`V_S^-`/domain-supported component, but the exact TSOM definition must be cited
and implemented before that label is used.

Set

\[
K_U=SU,
\qquad
C_U=K_U^\dagger H_S,
\qquad
R_U=(I-K_UK_U^\dagger)H_S,
\qquad
Q_U=UC_U.
\tag{3.1}
\]

`Q_Uh` is the minimum-current-norm element inside `Range(U)` whose data best
approximates the sampling perturbation `H_Sh`. `R_Uh` is the part of that
sampling perturbation that no retained current can reproduce.

### Proposition 1: lift-residual decomposition

**[proved-fd]** For every `H_S`,

\[
H_S=K_UC_U+R_U,
\qquad
K_U^*R_U=0,
\qquad
\operatorname{rank}C_U\le p,
\tag{3.2}
\]

where `p` is the pose dimension. Moreover, `C_Uh` is the unique
minimum-Euclidean-norm coefficient vector minimizing
`||K_Uc-H_Sh||`, and `Q_Uh` is the corresponding minimum-current-norm
representative inside `Range(U)`.

Proof: `K_UK_U^dagger` is the orthogonal projector onto `Range(K_U)`;
the remaining statements are the Moore--Penrose least-squares properties and
rank monotonicity.

### Corollary 1: exact representability and the vacuity test

**[proved-fd]** Exact representation for all pose directions holds iff

\[
\operatorname{Range}(H_S)\subseteq\operatorname{Range}(K_U),
\quad\text{equivalently}\quad R_U=0.
\tag{3.3}
\]

If the unrestricted `S` is full row rank and `U=I`, then `R_U=0` for every
`H_S`. Thus the bare statement “a Green-operator perturbation is equivalent to
a current perturbation” is automatic in the usual underdetermined measurement
geometry and carries no identifiability content.

### Corollary 2: SVD formula and instability

If `S=U_S Sigma V_S^*` and `U=V_r` contains retained nonzero singular modes,
then

\[
Q_r=V_r\Sigma_r^{-1}U_r^*H_S,
\qquad
R_r=(I-U_rU_r^*)H_S,
\tag{3.4}
\]

and

\[
\lVert Q_r\rVert\le
\frac{\lVert H_S\rVert}{\sigma_r(S)}.
\tag{3.5}
\]

The pose-equivalent current can therefore explode when weak singular modes are
admitted. A stable construction needs an explicit truncation/regularization and
rank-event policy; it cannot use an unqualified pseudoinverse.

### Definition: pose-equivalent current subspace

\[
\mathcal V_P^{(U)}=\operatorname{Range}(Q_U)
=\operatorname{Range}(UC_U)\subseteq\operatorname{Range}(U),
\qquad
\dim\mathcal V_P^{(U)}\le p.
\tag{3.6}
\]

This is a current-space object. Intrinsically, exact data-equivalent currents
are equivalence classes modulo `ker S`; the Moore--Penrose choice is a canonical
representative only after the current and data metrics are fixed. In infinite
dimensions, nonclosed range can make the pseudoinverse unbounded; closure and
regularization semantics are **[open]**.

## 4. State-consistent lift and the TriSpace graph

Write a retained current variation as `delta j=U delta c`. From (2.1) and
(3.2), introduce

\[
z=\delta c+C_Uh.
\]

Then the joint linearization becomes

\[
\delta y=K_Uz+R_Uh,
\tag{4.1}
\]

\[
MUz-D_{E_{\rm tot}}\delta\chi-D_Uh=0,
\qquad
D_U=H_D+MUC_U.
\tag{4.2}
\]

`R_Uh` is the irreducible data-space geometry signature. `D_Uh` is the
state-equation price of absorbing the representable sampling change into the
current coordinates. This yields the pose graph

\[
\Gamma_P^{(U)}=
\left\{\big(Q_Uh,R_Uh,D_Uh\big):h\in\mathbb R^p\right\}.
\tag{4.3}
\]

The graph, rather than a forced product of three commuting projectors, is the
proposed TriSpace object:

- sensing/current component `Q_Uh` and its relation to `V_S^+`;
- state/domain component `D_Uh` and its relation to the retained TSOM state
  directions;
- irreducible geometry component `R_Uh`.

### Proposition 2: joint local hiding condition

Let

\[
\mathcal N_D=MU\ker(K_U)+\operatorname{Range}(D_{E_{\rm tot}}),
\qquad
T_U=(I-P_{\mathcal N_D})D_U.
\tag{4.4}
\]

**[proved-fd]** A pose direction `h` can be hidden to first order by some
retained current coefficient and contrast perturbation while satisfying both
the data and state equations iff

\[
R_Uh=0
\quad\text{and}\quad
T_Uh=0.
\tag{4.5}
\]

Proof: because `Range(R_U)` is orthogonal to `Range(K_U)`, zero data variation
requires separately `R_Uh=0` and `z in ker(K_U)`. Equation (4.2) is then
solvable in `z,delta chi` exactly when `D_Uh` belongs to `mathcal N_D`.

Consequently, after declaring the physical gauge subspace `G`, a finite
dimensional local self-calibration condition is

\[
\ker\begin{bmatrix}R_U\\T_U\end{bmatrix}=\mathcal G.
\tag{4.6}
\]

This is an identifiability statement, not a nonlinear convergence theorem.
If `K_U` is full column rank, `ker K_U={0}` and
`mathcal N_D=Range(D_Etot)`.

**[retrieval-bounded candidate contribution]** The exact combination of the
canonical SOM-restricted lift, its orthogonal irreducible residual, and the
state-consistency witness (4.4)--(4.6) is the central candidate contribution.
Its individual ingredients are standard linear algebra; value must be judged
by the full-wave/SOM coupling, the semantic separation it enforces, and whether
the numerical algorithm gains a measurable calibration basin.

## 5. Unknown `G` means a moving coordinate system

At nonlinear iteration `k`, the basis itself depends on the current pose:

\[
U_k=U\{S(x_k),D(x_k),\tau_k\}.
\]

A viable algorithm must therefore do all of the following:

1. recompute or update the sensing operator, illumination, current, and retained
   basis at `x_k`;
2. align successive bases/projectors (for example by orthogonal Procrustes on
   a fixed-rank stratum) rather than comparing arbitrary singular-vector signs
   and rotations;
3. monitor retained singular gaps and threshold crossings;
4. damp or regularize `C_U` when small retained singular values amplify the
   lift;
5. solve map/current/pose updates using both the data residual and the state
   residual, not the data-equivalent current alone;
6. restart or use a soft spectral filter at rank events.

A working local update is a variable-projection or trust-region solve for
`(chi,x)` and the retained ambiguous current coefficients, with the
data-determined SOM component recomputed at each pose. The pipeline may test a
minimal Phase-I version using `U=V_r(S)`; it must not call that implementation
full TSOM unless the exact domain fold has been verified from the cited source.

## 6. Five required falsification experiments

### E1. Canonical lift, exactness, and vacuity

- Verify (3.2) numerically under whitening/realification.
- Sweep retained rank and singular cutoff; report `||R_U||`, `||Q_U||`,
  `sigma_min(K_U)`, and orthogonality residual.
- Demonstrate that unrestricted full-row-rank `S` makes `R=0` for arbitrary
  pose-data perturbations, so unrestricted equivalence is vacuous.
- Create a rank-deficient/limited-aperture or truncated case where `R_U` is
  nonzero.

### E2. Pose-dependent SOM-coordinate drift and rank events

- Perturb receiver pose and compare projectors `P_r(x)` rather than raw singular
  vectors.
- Measure projector distance/leakage, singular gaps, and a valid local
  Davis--Kahan/Wedin-style bound with all conditions stated.
- Deliberately cross a hard threshold and show the rank-event discontinuity;
  compare a soft spectral filter.

### E3. Pseudo-current versus physical induced-current change

- Receiver-only displacement: verify `H_S`, `Q_U`, and `R_U` by centered finite
  differences while `H_D=0`.
- Transmitter-only displacement: verify
  `delta j_phys=M^{-1}H_Dh` against finite differences of the solved current.
- Co-moving array: verify `B h=H_Sh+S delta j_phys` against the full forward
  finite difference.
- Quantify the null/current component lost by replacing the physical current
  with `S^dagger S delta j_phys`.

### E4. Data cancellation versus state-equation disambiguation

- Apply `delta j=-Q_Uh`; measure the residual data signature `R_Uh` and the
  state defect `D_Uh`.
- Project the state defect against admissible real contrast and invisible
  retained-current directions; verify Proposition 2 and the rank of
  `[R_U;T_U]`.
- Include a constructed gauge/confounded case and an anchored or asymmetric
  case. A failed identifiability condition is a retained negative result.

### E5. Implementable reduced self-calibration

- Generate small deterministic full-wave synthetic data at a displaced true
  array geometry.
- Compare at least: wrong-pose fixed SOM; direct joint full-wave
  contrast/pose least squares; and the geometry-lifted reduced update using
  data plus state residuals.
- Report map error, pose error, data residual, state residual, effective degrees
  of freedom, runtime, convergence basin, and failure cases over deterministic
  seeds/noise levels.
- Use small grids and a low-dimensional contrast basis first. No success claim
  is allowed from one favorable initialization.

## 7. Originality firewall

### Mature antecedents; do not sell as new

- contrast sources/equivalent currents and the Lippmann--Schwinger state model;
- classical SOM and twofold SOM as known-geometry current-space reduction;
- pseudoinverse minimum-norm solutions and range/null-space decompositions;
- self-calibration, blind calibration, joint image/calibration optimization;
- source/receiver extension in full-waveform inversion;
- joint inverse-scattering reconstruction and transmitter localization;
- radar/SAR autofocus and motion-error correction;
- variable projection, Schur complements, Fisher information, principal angles,
  and generic subspace perturbation bounds;
- virtual antennas or equivalent contrast-source experiments by themselves.

### Candidate research value, subject to literature and experiment checks

1. A canonical **SOM/TSOM-restricted** geometry-to-current lift, explicitly
   paired with the irreducible residual rather than claiming universal exact
   equivalence.
2. The physical/pseudo-current separation between receiver re-sampling and
   transmitter/internal-state change.
3. The state-consistency witness `T_U` and the exact local hiding condition
   using `[R_U;T_U]`.
4. A graph-based TriSpace organization that keeps current, data, state, map,
   and pose spaces type-correct and avoids fictitious `2^3` intersections.
5. A pose-dependent SOM coordinate algorithm with rank-event monitoring and
   a demonstrated self-calibration benefit or a clearly delimited failure
   regime.

### Not established

- global or exhaustive novelty;
- a unique canonical lift without specified metrics/regularization;
- continuum closed-range or stable-invertibility results;
- a verified universal formula for the TSOM domain fold;
- global nonlinear convergence or large-error calibration;
- real-hardware robustness, 3D Maxwell transfer, broadband antenna effects,
  mutual coupling, clock/phase-center errors, or unknown backgrounds;
- superiority over source/receiver-extension FWI outside the tested model.

## 8. Manuscript claim discipline (pre-experiment hypothesis)

The seed originally proposed claiming that the full-wave state equation
supplied a second discriminating witness. E4/E7/E11 rejected that mechanism,
and E5 rejected an algorithmic-benefit interpretation. The admissible language
after experiment is instead:

> We formulate and test a geometry-lifted SOM diagnostic in which the
> receiver-side Green-operator perturbation is decomposed into a minimum-norm
> retained-current lift and an orthogonal retained-basis leakage. A
> nuisance-space dimension bound identifies pose directions that the retained
> current and map variations can hide; the tested state-consistency variants do
> not recover those directions.

Do not use “first,” “no prior work,” “solves unknown geometry,” “lossless
self-calibration,” or “proves full-wave SLAM convergence.”

## 9. Post-experiment adjudication and corrected rank theorem

### 9.1 What survived

- **[proved-fd]** Proposition 1, Corollaries 1--2, and Proposition 2 remain
  algebraically correct after independent sign, range, and realification
  checks.
- **[model-conditional + numerically verified]** The world-fixed split between
  receiver re-sampling `H_S h` and transmitter/internal physical-current
  response `M^{-1}H_Dh` passed analytic-versus-finite-difference checks. Only
  `H_S` should be pseudo-current lifted in the independent-current model.
- **[numerically verified]** The unrestricted lift is vacuous at full row
  rank, whereas a declared retained basis produces a nonzero retained-basis
  leakage `R_U` and an orthogonal minimum-norm lift.
- **[numerically verified]** Hard singular-vector truncation can jump at rank
  events; a soft spectral filter makes the coordinate change smaller, although
  it does not itself solve self-calibration.

### 9.2 What failed or was narrowed

- **[falsified in tested regimes]** `T_U` did not discriminate good and bad
  pose estimates. Hard state equality and soft state penalties did not recover
  the data-hidden pose directions. Proposition 2 is therefore retained as a
  solvability statement, not advertised as a demonstrated recovery mechanism.
- **[not an algorithmic advantage]** When the visible pose rank is three, the
  tested reduced pose solver is an invertible reparameterization of the direct
  joint objective. Agreement validates implementation and same-basin behavior;
  it is not a speed, accuracy, or basin improvement.
- **[demonstrated loss]** When the visible pose rank is below three, the tested
  reduced solver freezes the complementary pose coordinates. Direct joint
  inversion recovered part of those coordinates and generally performed
  better.
- **[terminology correction]** `V_P^(U)=Range(Q_U)` is a scene- and
  metric-dependent pose-equivalent current tangent inside `Range(U)`, not a
  third independent spectral space. “TriSpace” may be retained only as the
  project/framework name for the typed triple `(Q_Uh,R_Uh,D_Uh)`.

### 9.3 Data-side nuisance-space theorem

Work in the whitened, realified data space `Y_R` of dimension `m`. Let

\[
\mathcal N_U=
\operatorname{Range}_{\mathbb R}\!\big(\mathcal R(K_U)\big)
+\operatorname{Range}(A_\chi)
\]

be the nuisance range generated by complex retained-current coefficients and
real contrast coordinates, and let `B in R^{m x p}` be the total local pose
Jacobian. Define

\[
B_{\rm vis}=(I-P_{\mathcal N_U})B,
\qquad q=\operatorname{rank}(B_{\rm vis}).
\tag{9.1}
\]

Here `\mathcal R` denotes the standard complex-to-real block embedding. The
data-hidden pose subspace is `ker B_vis` and has dimension `p-q`.

**Proposition 3 (dimension obstruction) [proved-fd].**

\[
\dim\ker B_{\rm vis}
\ge
\max\{0,\ p+\dim\mathcal N_U-m\}.
\tag{9.2}
\]

Proof: the range of `B_vis` lies in `\mathcal N_U^\perp`, hence
`q <= m-dim N_U`. Rank--nullity gives
`dim ker B_vis=p-q >= p+dim N_U-m`, together with nonnegativity.

For the executed model, `m=2M`, `p=3`, and the real contrast basis has
dimension `K=3`. Under generic column independence,

\[
\dim\mathcal N_U=\min(2M,2r+3).
\]

Therefore:

- `r <= M-3` leaves enough data codimension for three locally visible pose
  coordinates, but does not guarantee visibility;
- `r=M-2` leaves codimension at most one, so at least two pose directions are
  data-hidden;
- `r>=M-1` generically fills the real data space, so all three pose directions
  are data-hidden.

This is a dimension obstruction, not a special Fourier/Hankel theorem. The
locked autonomous report used a relative-only rank rule after nuisance-range
saturation and mistook roundoff-size singular values for three visible pose
directions. The parent recomputation fixes that interpretation without
changing the locked experiment artifacts.

### 9.4 Corrected scientific object

The defensible paper object is now:

1. a physically typed receiver pseudo-current / transmitter physical-current
   split;
2. a canonical retained-SOM lift and retained-basis leakage, with the
   unrestricted vacuity exposed;
3. a local nuisance-dimension certificate for which pose directions can be
   distinguished from retained current and map changes;
4. an explicit negative result showing that the tested state witness and
   hidden-direction-freezing solver do not supply self-calibration by
   themselves;
5. acquisition guidance: add independent transmitters, directionality,
   frequencies, frames, anchors, or priors until the *stacked* projected pose
   Jacobian has the required rank and conditioning.

This is a limits-and-design theory for self-calibrating SOM, not a claim that
unknown-array full-wave inverse scattering has already been solved.
