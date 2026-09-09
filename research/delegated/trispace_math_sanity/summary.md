# TriSpace SOM-SLAM math sanity review — summary

Scope: isolated mathematical research critic pass over
`fullwave_is_slam_spectral_observability_theory_rep.md`,
the two pasted re-scoping texts, and `Theory/SOM_SLAM_THEORY_CONTEXT.md`.
No novelty claims are made anywhere in this package; every statement carries an
epistemic label: **[theorem candidate]**, **[proposition candidate]**,
**[derivation]**, **[conjecture]**, **[open]**.

## 1. Can $dG_S(X)[\delta X]\,j$ be written as an equivalent current perturbation?

**[proposition candidate]** Exactly when

$$
dG_S(X)[\delta X]\,j \;\in\; \operatorname{Ran} G_S(X),
$$

in which case the minimum-norm equivalent current is
$\delta j_{\mathrm{eq}} = G_S(X)^\dagger\, dG_S(X)[\delta X]\,j$ and every other
representative differs by an element of $\ker G_S(X)$. In whitened realified
data the condition is read against $\operatorname{Ran}(W G_S)$. The residual is

$$
r \;=\; P_{\ker G_S^*}\, dG_S(X)[\delta X]\,j
  \;=\; P_{\operatorname{Ran} G_S^\perp}\, dG_S(X)[\delta X]\,j,
  \qquad \|r\|=\operatorname{dist}\!(v,\operatorname{Ran} G_S).
$$

**[derivation, verified in the harness]** The full pose tangent splits as
$B_l = dG_S[\delta X_l]j + G_S\,\delta j_{X,l}$ with
$\delta j_{X,l}=M^{-1}D_\chi\, dE^{\mathrm{inc}}[\delta X_l]$ (world-fixed grid),
so $B_l \in \operatorname{Ran} G_S \iff$ the sensing term lies in
$\operatorname{Ran} G_S$; the interior term is always in
$\operatorname{Ran} G_S$. Verified numerically at $N=8$, one pose: the split
identity holds to $2.6\times10^{-18}$ and, at rank $4/4$, the equivalent
current reproduces $B_l$ to $9\times10^{-16}$ with zero residual.

**[derivation]** Two degeneracies must be stated up front, or the question is
trivially answered:

1. With $m$ receivers and $S\ge m$ current pixels, generic $G_S$ has full row
   rank, so $\operatorname{Ran} G_S=\mathbb C^m$ and **every** pose-induced
   sensing perturbation is exactly representable at a single snapshot. A
   proper "pose subspace" $V_P$ only becomes non-degenerate across stacked
   poses/shared pose DOF or inside a retained (TSOM) current subspace.
2. If $\operatorname{Ran} G_S$ is not closed (continuum sensing), exact
   representation fails for $v\in\overline{\operatorname{Ran} G_S}\setminus
   \operatorname{Ran} G_S$; only asymptotic compensation with
   $\|\delta j_n\|\to\infty$ exists (the same closed-range pathology as the
   report's $\operatorname{Ran} B$ discussion).

## 2. Interior $G_D$ vs sensing $G_S$ perturbations

They are structurally asymmetric; conflating them is the central risk.

| axis | sensing path $dG_S[\delta X]j$ | interior path $\delta j_X = M^{-1}D_\chi dE^{\mathrm{inc}}[\delta X]$ |
|---|---|---|
| ambient space | data space; needs a pullback | already a current perturbation |
| residual eligibility | may leave $\operatorname{Ran} G_S$ (needs range condition) | automatically in $\operatorname{Ran} G_S$ (zero residual) |
| resolvent coupling | none (linear in fixed $j$) | amplified by $M^{-1}$; near resonance $\|M^{-1}\|$ diverges |
| pose dependence | body-fixed receiver path | Tx/illumination path; vanishes if Tx sits at the rotation center (verified: rotation column has zero interior part) |
| what it changes | which currents are *visible* (can reveal $\ker G_S$ content) | which currents are *physically realizable* |
| identifiability role | autofocus/phase-center bilinear degeneracy | excitation change; distinguishable in principle through multiple-scattering frequency dependence, shrinking to Born degeneracy as $\chi\to0$ |

**[derivation, verified]** Body-fixed rigid motions preserve the row-symmetries
of a rank-deficient $G_S$ (coincident receivers stay coincident), hence
$dG_S[\delta X]j$ stays inside $\operatorname{Ran} G_S$ for such motions
(residual $\approx 10^{-18}$). Per-channel phase/gain calibration breaks the
symmetry and produces a genuine residual (verified: single-channel ratio
$1/\sqrt2$, all-channel $\approx0.89$).

## 3. TriSpace self-calibration algorithm (coherent version)

**[conjecture / proposal]** A defensible version is exact block-coordinate
Gauss–Newton descent on one gauge-anchored objective (not ad hoc projector
products):

1. Anchor the gauge (fix pose 0 / world frame; else the $SE(d)$ gauge makes
   map and pose interchangeable for free).
2. At fixed $X$: SVD of $G_S(X)$ with gap-selected $V_S^+$, TSOM retained basis
   $U_T$, constrained contrast-source solve in $\operatorname{Ran} U_T$.
3. Pose block: minimize the residual with $B(X)$ and prior $J_X$; use the
   residual split $r_S$ (data no current can explain) vs $r_{\mathrm{SOM}}$
   (data the retained subspace cannot explain) to decide which pose
   coordinates are data-identified vs prior-identified.
4. Map block: update $\chi$ from the constrained current via the state
   equation; recompute $U_T$ only outside a trust region around rank events.
5. Certificates each iterate: gauge-anchored, $\sigma_{\min}$ of the reduced
   joint Hessian on the quotient, per-coordinate angles between pose-induced
   and retained-map data subspaces, rank-event flags.

Convergence and identifiability guarantees remain **[open]**: the problem is
non-convex, $V_P$ is scene-dependent ($\delta j_X$ and $j$ depend on $\chi$),
and at finite $m$ with full-row-rank $G_S$ the pose-only residual $r_S$
vanishes, so the algorithm reduces to a SOM-regularized joint inversion rather
than a genuinely new "third space".

## 4. Counterexamples and gauge/identifiability failures found

1. Single-snapshot full-row-rank degeneracy (equivalent current always exists).
2. Structure-preservation lemma: rank-deficient + rigid motion still gives
   zero residual; the naive "rank-deficient ⇒ residual" intuition is false.
3. Per-channel phase/gain calibration: genuine nonzero residual
   ($1/\sqrt2$ ratio), the canonical unrepresentable direction.
4. Compact non-closed range in $\ell^2$: asymptotic-only compensation,
   $\|\delta j_\varepsilon\|\sim(\pi/4\varepsilon)^{1/2}\to\infty$.
5. Global $SE(d)$ gauge: $B\xi_X=-A\xi_\chi$; the equivalent current exists
   with zero residual but equals a physical map change — pose and map are
   data-indistinguishable without an anchor.
6. Pose-only stabilizer (radially symmetric scene): $B\xi_X=0$, a pure pose
   null direction that is not a map confound.
7. Born narrowband/monostatic phase-wrap aliasing; interior vs sensing
   separation collapses as $\chi\to0$ (contrast-dependent identifiability).
8. Retained-subspace leakage: representable in all of $\mathcal H_J$ but not
   in $\operatorname{Ran} U_T$ ($r_{\mathrm{SOM}}\neq0$ while $r_S=0$);
   non-commuting $P_S^\pm, P_P$ falsify any "eight-cell" exact decomposition.
9. Moving-grid / robot-local $\chi$ parameterization makes $V_P$
   coordinate-dependent ($D_\chi(X)$ terms).

## 5. Five minimal falsification experiments (details in experiments.md)

1. E1 — decomposition identity + single-snapshot degeneracy.
2. E2 — rank-deficient leakage: rigid motion vs per-channel calibration.
3. E3 — interior/sensing residual asymmetry and near-resonance amplification.
4. E4 — gauge counterfeiting and pose-only stabilizer.
5. E5 — retained-subspace leakage, non-commutativity, rank event.

All five are implementable from `src/helmholtz.py`
(`build_operators`, `build_AB`, `whiten_realify`, `born_forward`,
`forward_measurements`, `pose_geometry`) plus `family1_pilot.py` /
`family2_algebraic_spine.py`; E4 and parts of E5 are partially covered by the
existing `family3_gauge_born.py` and `family8_som_distinction.py` runs.

## 6. Points that need parent Codex validation

- All proposition/theorem candidates above (esp. the range-condition statement
  and the bridge $B\delta X\in\operatorname{Ran} A \iff G_S^\dagger B\delta X
  \in \operatorname{Ran} T_\chi + \ker G_S$ under $r=0$).
- The "full-row-rank ⇒ single-snapshot $V_P$ is trivial" consequence.
- The structure-preservation lemma and whether per-channel phase calibration
  belongs inside the project's definition of pose $X$.
- The $\ell^2$ asymptotic-compensation growth constant.
- The Born phase-only factorization (derived under the report's far-field
  convention 4.15, not the harness's near-field Hankel setting).
- Algorithm convergence/identifiability status and gauge bookkeeping.
- Whether any of this constitutes novelty: **not judged here**.
