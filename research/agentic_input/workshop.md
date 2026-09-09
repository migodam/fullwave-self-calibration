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
16. For the declared Green function
    $g(r,r')=(i/4)H_0^{(1)}(k_b\lVert r-r'\rVert)$, the equal-area disk
    self-cell integral with $a=h/\sqrt{\pi}$ is
    $I_{\rm self}=\frac{i\pi a}{2k_b}H_1^{(1)}(k_ba)-k_b^{-2}$, so the
    diagonal domain-propagator entry is
    $[G_D]_{nn}=k_b^2I_{\rm self}=\frac{i\pi k_ba}{2}H_1^{(1)}(k_ba)-1$.
    The lower-limit term $-k_b^{-2}$ must not be dropped; without it the
    purported cell integral fails to vanish as $h\to0$. Also use
    $\nabla_s g(z,s)=+(ik_b/4)H_1^{(1)}(k_bR)(z-s)/R=-\nabla_zg(z,s)$ for
    source-position derivatives. A pilot produced before both corrections is
    not admissible evidence and must be rerun.

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

---

# Formal Workshop Seed — Full-Wave Inverse-Scattering SLAM Spectral Observability

> Compact, theory-first, falsification-first, novelty-first seed prepared from
> `Theory/SOM_SLAM_THEORY_CONTEXT.md`, `Theory/Questions/Q1.md`,
> `Theory/Questions/A1.md`, and `research/delegated/context_isolation/summary.md`.
> It preserves the existing theory and the internal (non-peer-reviewed) status
> labels used by those discussions. It solves no new proofs and asserts no
> novelty.

Status legend used throughout this folder:

- **[accepted]** — algebraically derived / adopted inside the prior discussions under the conditions stated in the ledger.
- **[conditional]** — believed correct only under explicitly listed conditions, or only in one regime/model.
- **[open]** — stated conjecture or unresolved question; not to be quoted as a result.
- **[non-claim]** — deliberately excluded or explicitly not claimed.

### Source and exclusion boundary

A read-only check of the user's ChatGPT `电磁逆散射` project on 2026-09-03
confirmed the seven core theory conversations and conversation IDs recorded in
`Theory/SOM_SLAM_THEORY_CONTEXT.md`. The separate recent conversation
`偏微分方程逆散射SLAM描述` restates the same forward-map, tangent-kernel, and
closed-range projection framing; it contributes no additional claim beyond the
local context. Older FAEMIS / learned-imaging experiment conversations in that
ChatGPT project concern a different algorithmic line and are deliberately
excluded from this paper seed to prevent semantic and experimental drift.

## 0. Fixed working assumptions (do not silently relax)

- 2D scalar Helmholtz contrast-source model, discretized, world-fixed map grid,
  uniform known background, static scene, measurements containing scattered
  field.
- Local linearization about a nominal, non-resonant state; all information
  statements are local first-order statements unless a different order is
  stated.
- Whitened measurements; realified operators whenever physical map/pose
  parameters are real; projectors, ranges, ranks, principal angles, and Schur
  complements computed in that whitened, realified space.
- Finite-dimensional setting unless an infinite-dimensional section explicitly
  says otherwise; gauge fixed, or analysis restricted to the identifiable
  support, or Moore--Penrose used without reading gauge variance as ordinary
  covariance.
- Computation target: Apple Silicon CPU. No CUDA/GPU infrastructure, daemon,
  database, web UI, or new orchestration framework.
- Literature mapping in the source documents is internal-discussion context,
  not a verified systematic search. No global novelty claims are permitted.

## 1. Core research statement (preserved)

> In full-wave inverse-scattering SLAM, how does finite-dimensional pose
> uncertainty deform the SOM/TSOM-recoverable subspace; is the deformation a
> low-rank generalized-eigenvalue defect; and can it be tracked stably across
> trajectory-induced rank-changing events?

Geometric equivalent: how close are the map tangent space and the pose tangent
space inside data space?

$$ \operatorname{Range}(A) \cap \operatorname{Range}(B) $$

describes exact local ambiguity; $\sin^2\theta_i$ describes continuous relative
retention; $A^*P_{B^\perp}A$ converts the geometry into absolute effective
information; $DK_{\mathrm{eff}}$ connects it to trajectory sensitivity.

## 2. Preserved theory skeleton

### 2.1 Forward model and Jacobians

Lippmann--Schwinger / contrast-source equations per configuration $t$:

$$
E^{\mathrm{sca}} = G_S J,\qquad
J = D_\chi(E^{\mathrm{inc}} + G_D J),
$$

$$ M_t = I - D_\chi G_{D,t}, \qquad
J_t = M_t^{-1} D_\chi E_t^{\mathrm{inc}},\qquad
F(\chi) = G_S M^{-1} D_\chi E^{\mathrm{inc}}. $$

Joint first-order model about $(\chi_0, X_0)$:

$$ \delta y = A\,\delta\chi + B\,\delta X + n,
\qquad A = D_\chi F,\quad B = D_X F. $$

[conditional] Map Jacobian, when the map is world-fixed and $M_t$ invertible:

$$ A_t = G_{S,t} M_t^{-1}\operatorname{diag}(E_t^{\mathrm{tot}}),
\qquad E_t^{\mathrm{tot}} = E_t^{\mathrm{inc}} + G_D j_t. $$

Basis-parameterized contrast $\chi = S_\chi c$ replaces the last factor by
$D_{E^{\mathrm{tot}}}S_\chi$.

[conditional] Pose Jacobian, full candidate form per pose direction $h$:

$$
B_t h = D_{x_t}G_{S,t}[h]\,j_t
+ G_{S,t}M_t^{-1}
\Big(D_\chi\big(D_{x_t}E^{\mathrm{inc}}[h] + D_{x_t}G_D[h]\,j_t\big)
+ (D_{x_t}D_\chi[h]) E_t^{\mathrm{tot}}\Big)
+ D_{x_t}y^{\mathrm{direct}}[h] + \cdots
$$

where the trailing terms appear only under robot-local map grids, moving
$G_D$, direct-path/calibration channels, and so on. The simplified form

$$ B_t h = D_{x_t}G_{S,t}[h]\,j_t
+ G_{S,t}M_t^{-1} D_\chi D_{x_t}E^{\mathrm{inc}}[h] $$

is only valid when $G_D$ and the map parameterization do not move with pose and
no direct/calibration term exists. Physical meaning: receiver-coordinate
derivative plus re-scattered transmitter-coordinate derivative.

### 2.2 Information operators

Known pose: $K_{\mathrm{IS}} = A^*A$.
Free deterministic pose (no prior):

$$ K_{\mathrm{SLAM}} = A^* P_{B^\perp} A,\qquad P_{B^\perp} = I - BB^\dagger. $$

Pose prior $J_X \succeq 0$ (Schur complement):

$$ K_{\mathrm{eff}} = A^*A - A^*B (B^*B + J_X)^\dagger B^*A. $$

Pose-induced loss:

$$ L_X = K_{\mathrm{IS}} - K_{\mathrm{eff}}
= A^*B(B^*B+J_X)^\dagger B^*A = H_X^* H_X. $$

Map-prior curvature is written separately: $K_{\mathrm{post}} = K_{\mathrm{eff}} +
K_{\mathrm{prior}}$; regularization never creates measurement Fisher
information.

### 2.3 Central identities

Set $\mathcal U = \operatorname{Range}(A)$, $\mathcal V = \operatorname{Range}(B)$.

[accepted]

1. $0 \preceq K_{\mathrm{SLAM}} \preceq K_{\mathrm{IS}}$, hence
   $\lambda_i(K_{\mathrm{SLAM}}) \le \lambda_i(K_{\mathrm{IS}})$.
2. $\ker K_{\mathrm{SLAM}} = \{u : Au \in \mathcal V\}$ (includes $\ker A$).
3. $\operatorname{rank}K_{\mathrm{IS}} - \operatorname{rank}K_{\mathrm{SLAM}}
   = \dim(\mathcal U\cap\mathcal V)$.
4. Equivalence:
   $\mathcal U\cap\mathcal V = \{0\}$ iff
   $\operatorname{rank}[A\ B] = \operatorname{rank}A + \operatorname{rank}B$;
   and $\dim(\mathcal U\cap\mathcal V) \ge \max(0, r_A + r_B - m)$ with $m$ the
   realified data dimension.
5. Directional retention
   $R(u) = \|P_{B^\perp}Au\|^2/\|Au\|^2 = \sin^2\angle(Au, \mathcal V)$.
6. Generalized retention spectrum on the observable support of $K_{\mathrm{IS}}$:
   $K_{\mathrm{SLAM}}v = \rho K_{\mathrm{IS}}v$ has $0\le\rho_i\le1$, and with no
   pose prior $\rho_i = \sin^2\theta_i$, where $\theta_i$ are principal angles
   between $\mathcal U$ and $\mathcal V$. At most $\operatorname{rank}B$
   generalized values differ from 1; ordinary eigenvectors may still rotate
   broadly under a low-rank loss.
7. Continuous confusion dimension
   $d_{\mathrm{conf}} = \sum_i\cos^2\theta_i = \operatorname{tr}(I - R_{\mathrm{op}})
   = \|Q_B^*Q_A\|_F^2$. Distinct quantities: destroyed DoF
   ($\#\{i:\theta_i=0\}$), affected modes ($\#\{i:\theta_i<\pi/2\}$), continuous
   confusion dimension. Volume retention on a task support $\mathcal S$ equals
   $\prod \sin^2\theta_i$.
8. $\operatorname{rank}L_X \le \operatorname{rank}B$ (with $L_0 =
   A^*P_\mathcal V A$ for no prior), Hermitian interlacing
   $\lambda_{i+p}(K_{\mathrm{IS}}) \le \lambda_i(K_{\mathrm{eff}}) \le
   \lambda_i(K_{\mathrm{IS}})$ for $\operatorname{rank}L_X \le p$.
9. Variational forms:
   $u^*K_{\mathrm{SLAM}}u = \min_z\|Au - Bz\|^2$ and
   $u^*K_{\mathrm{eff}}u = \min_z(\|Au-Bz\|^2 + z^*J_Xz)$, giving ordering and
   monotonicity in $J_X$.
10. Exact $SE(2)/SE(3)$ gauge forces $A\xi_\chi + B\xi_X = 0$, so
    $\xi_\chi \in \ker K_{\mathrm{SLAM}}$ and $\rho=0$. The converse is false:
    a local tangent coincidence need not come from a global symmetry.
11. Multi-frequency with one shared physical pose compensation:
    $u\in\ker K_{\mathrm{SLAM}}^{\mathrm{multi}}$ iff $\exists z:\ A_fu = B_fz$
    for every frequency $f$. Consistent stacking never decreases absolute
    $K_{\mathrm{eff}}$ for a fixed pose prior. In the no-prior projector case,
    an exact repeated block $(A_2,B_2)=c(A_1,B_1)$ scales both
    $K_{\mathrm{IS}}$ and $K_{\mathrm{SLAM}}$ by $1+|c|^2$, so every generalized
    retention value $\rho_i$ is unchanged. With a fixed nonzero $J_X$, this
    invariance generally fails because repeated data change the data-to-prior
    weighting; invariance is recovered only under a matching prior scaling.
12. Born model at non-zero nominal scene: $B_\ell = (\partial A/\partial
    X_\ell)\chi_0$; at empty background $\chi_0=0$ the first-order $B=0$ and
    $K_{\mathrm{SLAM}}=K_{\mathrm{IS}}$, but geometry mismatch enters through the
    second-order bilinear term $(D_XA[\delta X])\delta\chi$. Two different
    questions: joint identifiability around non-zero $\chi_0$, and
    sensing-operator uncertainty at empty background.

### 2.4 Conditional results (always state conditions when used)

- Prior limits $J_X\to0$/$J_X\to\infty$ require a fixed-rank regular path
  (e.g. $J_X = \epsilon I$) and prior support covering every data-coupled pose
  direction; unanchored graph gauges never converge to $K_{\mathrm{IS}}$.
- Finite pose prior gives weighted shrinkage, not ordinary principal angles:
  with $D = BJ_X^{-1/2}$,
  $I - B(B^*B+J_X)^{-1}B^* = (I + DD^*)^{-1}$ and
  $R_X = Q_A^*(I+DD^*)^{-1}Q_A$; its eigenvalues lie in $(0,1]$ but are not
  $\sin^2\theta_i$ in general.
- Trajectory derivative formulas
  $\dot K_{\mathrm{eff}} = \dot A^*WA + A^*W\dot A + A^*\dot WA$,
  $\dot W = -\dot BHB^* - BH\dot B^* + BH\dot C HB^*$, require $C =
  B^*B+J_X$ invertible or a consistent generalized-Schur convention. $\dot B$
  needs the pose Hessian $D_X^2F$; omitting it is the explicit
  frozen-nuisance-subspace approximation.
- No-prior projector derivative
  $\dot P_B = P_{B^\perp}\dot B B^\dagger + (P_{B^\perp}\dot B B^\dagger)^*$
  and the bound
  $\|\dot K_{\mathrm{SLAM}}\| \le 2\|A\|\|\dot A\| + 2\|A\|^2\|B^\dagger\|
  \|\dot B\|$ require $B$ of locally constant rank.
- Robust objective surrogate
  $\gamma_r(X,\varepsilon) = \lambda_r(X) - \varepsilon\|\nabla_X\lambda_r\|_*
  + O(\varepsilon^2)$ requires simple $\lambda_r$, positive spectral gap,
  second-order smoothness, and a full norm-ball uncertainty set. A certified
  lower bound $\lambda_r(X+\Delta X)\ge\lambda_r(X)-L_K\|\Delta X\|$ needs a
  valid *uniform* operator-Lipschitz constant for $K$ over the whole uncertainty
  ball. A nominal derivative norm or sampled maximum is only a local/empirical
  diagnostic, not such a certificate.
- Eigenvector sensitivity and $\dot\rho_i$ formulas require simple eigenvalues
  (or compressed-operator treatment) and a phase gauge.
- Full-wave operator bounds
  $\|A_t\| \le \|G_{S,t}\|\,\|M_t^{-1}\|\,\|E_t^{\mathrm{tot}}\|_\infty$ and
  the analogous $B$ bound hold under the stated model; near
  $\sigma_{\min}(M_t)\to0$, absolute information and sensitivity may both grow,
  and linearization radius shrinks.
- Far-field/Born asymptotics ($K_{\mathrm{IS}} = O(R^{-2})$ in 2D monostatic,
  $k_{\perp,\max}$ evanescent reach, phase-modulation vs sampling-location
  terms) are regime-dependent with convention-dependent constants.
- Infinite-dimensional versions require closed ranges / polar or frame
  formulations and stable transversality $\inf_{Au\ne0}
  \mathrm{dist}(Au,\mathcal V)/\|Au\| > 0$; finite-dimensional dimension
  formulas do not transfer automatically.

### 2.5 Open claims (must never be presented as results)

- Generic transversality: for full angular diversity, multiple genuinely
  distinct frequencies, and sufficient MIMO diversity, are map/pose tangents
  separated except for global $SE(2)/SE(3)$ gauge?
- Single-frequency finite-aperture characterization of near-confounded modes and
  missing wedge; whether bandwidth gives a non-trivial lower bound on
  $\theta_{\min}$.
- Whether full-wave resonance information gain outweighs the increase in
  $DK_{\mathrm{eff}}$ and model-mismatch sensitivity.
- Fourier surrogate error for $V_D^+$ (subspace angle / operator norm) and its
  composition with the pose-defect layer.
- Continuity/rank-event theory: smoothness on fixed-rank strata, crossing update
  rules, hysteresis/flag bookkeeping; spectral classification threshold/gap
  stability for $\mathcal E_G/\mathcal E_M/\mathcal E_R$.
- Exact definition of two-fold SOM $V_D^\pm$ in any specific cited Chen version
  and discrete convention; whether hierarchical projector composition (never an
  exact octant decomposition) is valid under which commutativity/gap conditions.
- Prior-weighted spectrum bounds and monotonicity for
  $Q_A^*(I+DD^*)^{-1}Q_A$; joint motion-graph gauge support in generalized Schur
  complements; layered calibration/clock/phase-center defect ranks.
- Nonlinear/global basin questions: cycle skipping, phase ambiguity, and any
  basin-of-attraction bound relative to the local retention spectrum.

## 3. Deliberate non-claims (guard list)

- SOM mode counts describe stable data-resolvable current directions, not
  intrinsic scene dimensionality; a 3D output does not imply every voxel is
  data-supported.
- "Theoretically unique", "large local Fisher information", "algorithm
  converges", and "output looks reasonable" are separate propositions.
- Regularization does not manufacture measurement Fisher information.
- ML modules are outside the core theory until objectives, equivariance,
  calibration, and data-consistency guarantees exist.
- No global novelty claim; closest-prior-work boundaries from the discussions
  are unverified and require a fresh systematic search before submission.

## 4. Non-negotiable caveats

### 4.1 Whitening and realification

All projectors and ranks above are valid in the whitened space. If noise is
colored, first apply $\widehat A = \Sigma_n^{-1/2}A$, $\widehat B =
\Sigma_n^{-1/2}B$, or carry the metric $\Sigma_n^{-1}$ through every projection.
Then realify real physical parameters:

$$ A_{\mathbb R} = \sqrt2\begin{bmatrix}\Re\widehat A\\ \Im\widehat A\end{bmatrix},
\qquad B_{\mathbb R} = \sqrt2\begin{bmatrix}\Re\widehat B\\ \Im\widehat B\end{bmatrix}. $$

Complex contrast requires $(\chi_R,\chi_I)$ as two real blocks
$\begin{bmatrix}\Re\widehat A & -\Im\widehat A\\ \Im\widehat A &
\Re\widehat A\end{bmatrix}$. Using a complex-linear $BB^\dagger$ for real pose
increments is forbidden: it enlarges the nuisance range and overestimates
map--pose confounding.

### 4.2 Gauge discipline

State whether gauge is fixed, analysis is on the identifiable quotient/support,
or Moore--Penrose is used without reading gauge variance as covariance. Exact
gauge implies $\rho=0$; the converse is not asserted. Physical gauge also
depends on anchors, known background, boundaries, and scene symmetries.

### 4.3 Born empty-background caveat

Never write "Born at $\chi_0=0$ proves pose errors are harmless." The correct
statement is that first-order pose Jacobian vanishes there, so pose-map
confusion becomes a second-order bilinear mismatch problem with a different
first-order operator.

### 4.4 Semantic drift ban: current-space $G_S$ vs map Jacobian $A$ (mandatory)

| Object | Where it lives / what it means | What it must never be conflated with |
|---|---|---|
| $G_S$ SVD, $V_S^\pm$ | current space: which induced-current modes reach the array | map-parameter modes, $K_{\mathrm{eff}}$ eigenvectors, map observability spectrum |
| $A = D_\chi F$ (full wave) | $A = G_SM^{-1}\operatorname{diag}(E^{\mathrm{tot}}) \ne G_S$; data-space tangent from map perturbations | $G_S$ itself; $G_S$'s right singular subspace |
| Born $A$ | $A \approx G_S\operatorname{diag}(E^{\mathrm{inc}})$; still a map Jacobian | $G_S$ or current modes |
| $\operatorname{Range}(A)$ | data-space first-order map changes | $\operatorname{Range}(G_S)$, MSR/MUSIC data subspaces, current span |
| $\operatorname{Range}(B)$ | data-space pose-induced changes | the pose coordinate space itself |
| $K_{\mathrm{eff}}$ eigenmodes | map/reduced-coordinate effective information after pose elimination | SOM current modes, $V_S^\pm$, plain per-mode eigenvalue ratios of $K_{\mathrm{IS}}$ |
| $\rho_i$ | Fisher-whitened principal-direction retention | $\lambda_i(K_{\mathrm{SLAM}})/\lambda_i(K_{\mathrm{IS}})$, or a per-$K_{\mathrm{IS}}$-eigenmode ratio |

Current-space and map-space structures connect only through an explicit map
like $T_\chi = \partial J/\partial\chi$ or the pullback
$K_{\chi,\mathrm{pose}} = T_\chi^*K_{\mathrm{pose}}^J T_\chi$, never by
switching names. Rank thresholds used to verify exact algebraic identities must
be machine-precision type, not SOM dominant-subspace thresholds.

### 4.5 Five limits stay separate

Uniqueness/identifiability, local observability, stability, statistical
precision (FIM/CRB), and global multimodality/cycle-skipping are not
interchangeable. Every claim must name which one it is about.

## 5. Falsification-first agenda

Run the five experiment families in `experiments.md`; each must seek evidence
against the named claim before counting it as supported. Priority falsifiers:

1. Rank identity, $\rho=1-\sigma^2(Q_B^*Q_A)$, loss-rank bound, or interlacing
   fails beyond backward-error-scaled tolerance in whitened/realified data.
2. $K_{\mathrm{SLAM}} \not\preceq K_{\mathrm{eff}}(\alpha) \not\preceq
   K_{\mathrm{IS}}$ anywhere in a prior sweep.
3. $SE(2)$ generators do not land in $\ker K_{\mathrm{SLAM}}$ (gauge failure),
   or a zero-retention direction is later found to be non-gauge, contradicting
   the presumed "gauge explains $\rho=0$" reading.
4. Multi-frequency stacking decreases absolute $K_{\mathrm{eff}}$, or a genuine
   diversity block changes nothing while a duplicate block changes $\rho$ --
   either would overturn the stated mechanism.
5. First-order trajectory predictions or the robust surrogate fail by more than
   documented $O(\varepsilon^2)$ terms in their declared validity region, or a
   lower bound backed by a genuinely uniform operator-Lipschitz constant is
   crossed. A bound built only from sampled or nominal derivatives is never
   treated as certified.
6. At nominal (non-resonant) configurations, finite-difference Jacobian checks
   fail, $M$ is effectively singular, or Born/full-wave switching does not
   produce the predicted contrast-dependence.

## 6. Novelty-first discipline

- Mature components (SOM current decomposition, SVD/principal angles,
  Schur/EFIM nuisance elimination, finite-rank interlacing, gauge up to group
  action, Grassmann/projector perturbation, A/D/E-optimal design) are not
  claimed as novel.
- Candidate research space is the *combination*: full-wave volumetric
  contrast-source SLAM + map-pose nuisance tangent geometry + retention
  spectrum + trajectory/robustness consequences. It is a hypothesis, not a
  finding.
- Any "closest work" comparison must be re-verified by a fresh systematic search
  (later Codex-scheduled `scholarqa-research` per project AGENTS.md) before
  submission. No citation in these seeds is an independently verified
  bibliographic claim.
- Do not state "no prior work" or "first" in any artifact or publication draft.

## 7. One-sentence preservation of the user's thesis

SOM/TSOM first compresses current ambiguity through $G_S$ and internal domain
structure; the SLAM Schur/pose defect then characterizes the residual map--pose
ambiguity inside the retained current/map space as a low-rank,
generalized-eigenvalue retention defect whose stability depends on trajectory
geometry, priors, frequency diversity, and rank-event structure.

---

# Experiment Plan — Five Families

> Organized from the checks proposed in `Theory/SOM_SLAM_THEORY_CONTEXT.md`
> ($\S19), `Theory/Questions/A1.md` ($\S19), and
> `research/delegated/context_isolation/summary.md`. Each family is one coherent,
> implementable unit with a falsifiable hypothesis, controlled variables,
> metrics, pass/fail interpretation, and an explicit statement of what it cannot
> establish. All runs assume 2D scalar Helmholtz on Apple Silicon CPU, small
> deterministic matrices, no GPU/database/orchestration infrastructure.

## Shared scenario (reference configuration, reused by every family)

- Domain $D=[-0.5,0.5]^2$; grid $32\times32$ by default with a $40\times40$
  sensitivity run.
- Nominal contrast $\chi_0 = 0.3\,\chi_{\mathrm{disk},1} + 0.5\,\chi_{\mathrm{disk},2}$
  (two smooth-edged blobs), deliberately non-zero. A smooth
  Fourier/spline/basis variant is available whenever infinitesimal translation
  or rotation must remain in the parameter space (gauge checks); piecewise-
  constant pixels are used only where the claim tolerates it.
- Sensors: one Tx plus 4--8 body-fixed Rx per pose; pose $x_t=(p_x,p_y,\theta)$.
- Trajectory set: straight line, $90^\circ$ arc, $180^\circ$ arc, $360^\circ$
  circle; total path length and number of measurement configurations are held
  equal across trajectories within every comparison.
- Physics switches: Born vs full-wave; single vs multiple genuinely distinct
  frequencies; duplicate scaled block $(A_2,B_2)=c(A_1,B_1)$ as a control.
- Nominal points must satisfy invertibility of $M_t$ with a reported
  $\sigma_{\min}(M_t)/\|M_t\|$ margin; resonance cases are labeled separately.
- Noise/whitening: identity covariance default, plus one colored-covariance case
  per family where projections are re-derived with the metric
  $\Sigma_n^{-1}$ (equivalently whiten first). Realification is always applied
  after whitening for real pose/map parameters.
- Rank tolerance for *exact* algebraic checks: machine-precision type,
  $\tau = \max(m,n)\,\epsilon_{\mathrm{mach}}\,\sigma_1$; SOM physics thresholds
  are never used to verify identities.
- Every family returns: status of each claim it touched, exact executed checks,
  pass/fail, artifacts, and a "cannot establish" paragraph. Suggested output
  root: `research/experiments/family-NN/` with configs, logs, and hashes.

---

## Family 1 — Forward model and analytic Jacobians: derivative consistency and linearization validity

**Claims touched:** conditional $A_t$ and $B_t$ formulas; conditional
Born/full-wave regime; resolvent sensitivity caveats.

**Hypothesis (falsifiable).** In the reference scenario (full-wave, world-fixed
map, nominal non-zero $\chi_0$, $M$ well-conditioned), the analytic map
Jacobian $A_t = G_{S,t}M_t^{-1}\operatorname{diag}(E_t^{\mathrm{tot}})$ and the
pose Jacobian $B$ built from the stated coordinate convention agree with
centered finite differences of the forward map to within a documented tolerance
that scales as $O(h^2)$ plus solver error. Along a predeclared weak-scattering
contrast continuation away from resonance, the Born/full-wave *model
discrepancy* should vanish as contrast tends to zero and initially grow in a way
consistent with the Neumann expansion; strict global monotonicity is not
claimed.

**Controlled variables.**
- Fixed domain, grid, nominal contrast, sensor geometry, pose set, and
  measurement count.
- Toggle: Born/full-wave; smooth basis vs pixel map (basis variant only for
  subset of gauge-relevant poses).
- Finite-difference step $h$ swept over a declared range (e.g. $10^{-4}$ to
  $10^{-1}$ relative) to expose truncation vs cancellation.
- State-equation solver tolerance fixed and logged; incident field and Green
  function implementations shared across forward and Jacobian code.

**Metrics.**
- Per-block relative derivative errors:
  $\max_t\|D_hF_t[\delta\chi]-A_t\delta\chi\|/
  \max(\|A_t\delta\chi\|,\epsilon)$, where
  $D_hF_t[\delta\chi]=[F_t(\chi_0+h\delta\chi,X_0)-
  F_t(\chi_0-h\delta\chi,X_0)]/(2h)$, and the analogous centered pose
  derivative with random unit $\delta\chi$, $\delta X$ (several seeds).
- Observed convergence order of FD error vs $h$ (expect near 2 where
  truncation dominates, before roundoff floor).
- Residual of the nonlinear state equation at the nominal point;
  $\sigma_{\min}(M)/\|M\|$ margin; maximum field magnitude.
- Born-vs-full-wave forward and Jacobian discrepancy versus a
  contrast/spectral-radius knob (e.g. scaling the blob amplitudes), with the
  weak-scattering asymptotic region identified separately from resonant or
  cancellation-dominated settings.

**Pass/fail interpretation.**
- Pass: documented relative errors at or below the combined tolerance
  (target $10^{-5}$--$10^{-3}$ depending on $h$, grid, and conditioning) with
  consistent $h$-scaling; Born/full-wave discrepancy tends to zero toward the
  weak-scattering limit and its observed local scaling is reported; $M$ margin
  is reported and stays away from the declared resonance floor.
- Fail: large, $h$-insensitive mismatches; formula errors traceable to a single
  missing term (then the formula is corrected only with the term and the
  condition list updated); Born/full-wave discrepancy fails to approach zero in
  the weak-scattering limit without explanation; $M$ is singular at nominal
  settings with no documented reason.

**What it cannot establish.** Real-world validity beyond the discretized model;
behavior inside resonances/bifurcations; nonlinearity of pose over large
excursions; any claim that matching Jacobians implies matching physics of real
antennas or 3D Maxwell.

---

## Family 2 — Algebraic spine: kernel/rank/retention identities, ordering, prior sweeps

**Claims touched:** accepted identities 1--9 and the conditional prior-limit /
interlacing statements.

**Hypothesis (falsifiable).** In the whitened, realified reference data space
with machine-precision rank tolerance:
$\ker K_{\mathrm{SLAM}} = \{u: Au\in\operatorname{Range}(B)\}$,
$\operatorname{rank}K_{\mathrm{IS}}-\operatorname{rank}K_{\mathrm{SLAM}} =
\dim(\mathcal U\cap\mathcal V)$,
$\operatorname{eig}(R_{\mathrm{op}}) = 1 - \sigma_i^2(Q_B^*Q_A)$,
$0\le\rho_i\le1$, at most $\operatorname{rank}B$ generalized retention values
differ from 1, $\operatorname{rank}L_X\le\operatorname{rank}B$, interlacing
holds, and the prior sweep satisfies
$K_{\mathrm{SLAM}}\preceq K_{\mathrm{eff}}(\alpha)\preceq K_{\mathrm{IS}}$ with
the two documented limits along $J_X=\alpha I$ for
$\alpha\in[10^{-6},10^4]$.

**Controlled variables.**
- One shared nominal scenario set; trajectories and frequency sets as in the
  shared scenario.
- Whitening matrix: identity and one colored case.
- Prior: $J_X = \alpha I$ log-spaced sweep; one $J_X$ with singular support
  (partial/graph-style gauge) to expose limit caveats.
- Rank bookkeeping: thin SVD of $A$ and $B$, $Q_A,Q_B$, observable support of
  $K_{\mathrm{IS}}$; fixed $\tau$ as defined above; every tolerance and basis
  choice logged.

**Metrics.**
- Maximum scaled deviations from each identity (e.g.
  $\|K_{\mathrm{SLAM}} - A^*P_{B^\perp}A\|/\|A^*A\|$,
  $\|R_{\mathrm{op}} - V(I-Q^*ZZ^*Q)V^*\|$, rank differences, count of
  $\rho_i$ differing from 1).
- Generalized-spectrum range and gap; ordinary-vs-generalized eigenvector
  rotation sample statistics (to confirm low-rank loss can rotate many ordinary
  eigenvectors).
- Ordering violations: smallest eigenvalue of each PSD difference
  $K_{\mathrm{eff}}(\alpha)-K_{\mathrm{SLAM}}$ and
  $K_{\mathrm{IS}}-K_{\mathrm{eff}}(\alpha)$ must be $\ge -\delta$ for a
  backward-error-scaled $\delta$.
- Convergence curves of $K_{\mathrm{eff}}(\alpha)$ to both limits; volume-ratio
  $\prod\rho_i$ and $d_{\mathrm{conf}}$ compared to their closed forms.

**Pass/fail interpretation.**
- Pass: all identities hold to documented, conditioning-scaled accuracy; no
  ordering violation; limits and count statements reproduce; in the no-prior
  projector case repeated blocks do not change $\rho$ (checked again in Family
  4).
- Fail: any identity violated beyond tolerance; ordering violation in the
  interior sweep; limit failure explained only by singular-support prior (then
  the caveat, not the unconditional claim, is confirmed).

**What it cannot establish.** Infinite-dimensional transfer; physical
reconstruction quality; nonlinear/global identifiability; continuity across
rank events (Family 5 territory); the claim that eigenvectors of $K_{\mathrm{IS}}$
retain or lose information individually.

---

## Family 3 — Gauge and Born empty-background degeneration

**Claims touched:** accepted gauge item 10 and Born item 12; conditional
basis-representation caveats.

**Hypothesis (falsifiable).** (a) In a symmetry-compatible smooth-basis
representation with uniform background and unanchored scene, the three $SE(2)$
generators (translation $x$, translation $y$, rotation) satisfy
$\|A\delta\chi_g + B\delta X_g\|/(\|A\delta\chi_g\| + \|B\delta X_g\|) \approx 0$
and lie in $\ker K_{\mathrm{SLAM}}$. The expected numerical floor is algebraic
roundoff only when the discrete model is exactly equivariant; otherwise it is a
declared discretization/interpolation error that must decrease under refinement.
(b) In Born with $\chi_0=0$, the first-order pose
Jacobian is numerically zero so $K_{\mathrm{SLAM}}=K_{\mathrm{IS}}$, while the
map--pose response is dominated by the bilinear term
$(D_XA[\delta X])\delta\chi$ over a declared small-perturbation range.
(c) Adding an anchor / known-background / fixed-boundary feature removes the
expected subset of gauge directions.

**Controlled variables.**
- Map basis: smooth (Fourier/spline) vs pixel; only smooth basis is accepted
  for exact infinitesimal gauge checks.
- Background/scene: uniform unanchored; anchored variant; known-background
  variant; one symmetric-scene variant (rotation generator expected to
  degenerate).
- Physics: Born/full-wave; $\chi_0$ non-zero (gauge part) or exactly zero
  (empty-background part).
- Perturbation scales $\|\delta X\|$, $\|\delta\chi\|$ swept over a declared
  range for the bilinear dominance check.

**Metrics.**
- Gauge residual ratio defined above for each generator; distance of generator
  map vectors to $\ker K_{\mathrm{SLAM}}$ (cosine/sine angles).
- Relative sizes: $\|(D_XA[\delta X])\delta\chi\|$ vs $\|A(X_0)\delta\chi\|$ and
  vs the remaining second-order remainder
  $\|F(\delta\chi,X_0+\delta X)-A(X_0)\delta\chi-(D_XA[\delta X])\delta\chi\|$
  as scales shrink; observed exponent of the bilinear term in $\delta X$.
- Rank/kernel dimension of gauge subspace per scene variant; $\rho_{\min}$
  before and after anchoring.

**Pass/fail interpretation.**
- Pass: gauge residuals reach the predeclared floor implied by the independently
  measured forward/Jacobian/discretization errors and decrease under refinement
  unless exact discrete equivariance makes them algebraic-roundoff limited;
  pixel basis shows the documented failure mode (translation leaves the pixel
  space) rather than being used as counterevidence; empty-background behavior
  matches (b); anchors remove expected generators.
- Fail: group generators produce residual well above numerical floor in the
  smooth basis; bilinear term not dominant where predicted; anchoring changes
  nothing.

**What it cannot establish.** That every $\rho=0$ direction is a group gauge
(converse is open); that a finite, bounded, boxed scene possesses exact global
$SE(2)$ invariance in experiments (boundaries/anchors break it); physical
3D/Maxwell gauge behavior; any statement about nonlinear cycle-skipping.

---

## Family 4 — Frequency stacking and trajectory geometry effects on map--pose separation

**Claims touched:** accepted multi-frequency items 11--12; open generic-
transversality conjectures (sampled scenarios only).

**Hypothesis (falsifiable).** (a) Consistent stacking of genuinely distinct
frequency blocks never decreases absolute effective information:
$K_{\mathrm{eff}}^{(1:F+1)}\succeq K_{\mathrm{eff}}^{(1:F)}$. (b) The duplicate
control $(A_2,B_2)=c(A_1,B_1)$ increases absolute information but leaves all
$\rho_i$ unchanged in the no-prior projector case. With a fixed finite prior,
the duplicate instead changes normalized retention because it changes the
data-to-prior weight; scaling the prior by the same $1+|c|^2$ restores the
invariance. (c) A direction fully confounded at single frequency
($\rho\approx0$ because $Au=Bz_f$) is killed in the multi-frequency kernel only
if no single shared compensation $z$ fits all frequencies; adding a genuinely
diverse frequency can move such directions to $\rho>0$. (d) With path length
and measurement count fixed, trajectory geometry orders the retention/defect
statistics as $360^\circ$ circle $\succ$ $180^\circ$ arc $\succ$ $90^\circ$ arc
$\succ$ straight line for the sampled map/sensor configuration, or produces a
documented counterexample.

**Controlled variables.**
- Fixed domain, grid, $\chi_0$, sensors, noise level, and per-configuration
  SNR.
- Frequency sets: $F\in\{1,2,3\}$ genuinely distinct; duplicate scaled set as
  control; identical pose compensation variable $z$ across all blocks by
  construction.
- Trajectories with equal path length and equal number of measurement
  configurations; one shared pose prior setting ($J_X=0$ and one $\alpha>0$).

**Metrics.**
- Absolute spectra $\lambda_i(K_{\mathrm{IS}})$ and $\lambda_i(K_{\mathrm{eff}})$
  per $F$; PSD-difference checks for stacking monotonicity.
- Retention sets $\{\rho_i\}$, $\theta_{\min}$, $d_{\mathrm{conf}}$, defect rank,
  gauge-kernel dimension, before/after each added frequency.
- For each single-frequency-confounded direction: does it survive in
  $\ker K_{\mathrm{SLAM}}^{\mathrm{multi}}$? Quantify movement of $\rho$ with
  diverse vs duplicate addition.
- Duplicate-block control under three settings: $J_X=0$, fixed finite $J_X$,
  and finite $J_X$ scaled by $1+|c|^2$; verify invariance only in the first and
  third settings.
- Trajectory ranking of the four geometries by the above statistics, with
  confidence based on the fixed sampling budgets.

**Pass/fail interpretation.**
- Pass: monotonicity holds; the qualified duplicate-block predictions hold; at least some
  single-frequency confounded directions are separated by genuinely diverse
  frequencies (if none separate, that is a significant negative result refining
  the transversality conjecture); trajectory ordering reproduces the stated
  hypothesis or is replaced by a reproducible counterexample.
- Fail: stacking decreases $K_{\mathrm{eff}}$; a duplicate block alters $\rho$
  in the no-prior case (or fails the matching-scaling prediction); the
  shared-compensation structure is contradicted by construction-internal checks.

**What it cannot establish.** Generic transversality for all scenes (only finite
samples); global optimality of trajectories; real multipath association and
cycle-skipping behavior; engineering observability vs stability of an online
system.

---

## Family 5 — Trajectory sensitivity, robust surrogate, and rank-event stability

**Claims touched:** conditional derivative formulas, robust objective surrogate,
continuity/rank-event open items, frozen-nuisance approximation.

**Hypothesis (falsifiable).** (a) For simple eigenvalues with positive gap, the
first-order prediction $\lambda_i(X+\Delta X)-\lambda_i(X) \approx
v_i^*DK_{\mathrm{eff}}[\Delta X]v_i$ holds with error scaling as
$O(\|\Delta X\|^2)$ over a declared radius. (b) On a full norm ball, the
surrogate $\lambda_r - \varepsilon\|\nabla_X\lambda_r\|_*$ matches the sampled
worst case to first order. If a uniform operator-Lipschitz constant $L_K$ over
the complete ball is actually derived, the certified lower bound
$\lambda_r - \varepsilon L_K$ is never violated (it may be loose); a nominal
derivative norm or sampled maximum is labeled empirical and is not used as a
certificate. (c) While
rank and external spectral gap stay fixed, the defect subspace evolves
smoothly; near gap closure or small $\sigma_{\min}(B)$, fixed-rank tracking
degrades measurably and requires reinitialization/hysteresis. (d) The
frozen-nuisance approximation (dropping $\dot B$) has bounded, quantifiable
error in the sampled regime.

**Controlled variables.**
- Nominal trajectory set from the shared scenario; $\Delta X$ drawn on
  $\|\Delta X\|=\varepsilon$ spheres with fixed seeds, several
  $\varepsilon$ values.
- Spectral regime: simple-eigenvalue/gap cases selected deliberately, plus
  near-crossing/rank-event cases engineered by perturbing trajectory or grid.
- Prior: $J_X=0$ and one $\alpha>0$; norm/dual-norm choice explicit and fixed.
- Derivative implementations: full $DK_{\mathrm{eff}}$ (including $\dot B$ terms
  via pose Hessians) vs frozen-nuisance variant, both compared to finite
  differences.

**Metrics.**
- Relative first-order prediction error and its scaling exponent in
  $\varepsilon$; worst-case sampling minimum vs surrogate; derivation and
  coverage of any uniform $L_K$ used for a certified bound. If only a nominal
  or sampled $L_K$ is available, report its empirical violations without
  calling them violations of a theorem.
- Subspace movement: $\|U(X+\Delta X)U(X+\Delta X)^* - U(X)U(X)^*\|$ per
  perturbation, plus gap and $\sigma_{\min}(B)$ logs along perturbation paths.
- Errors introduced by the frozen-nuisance approximation, expressed in the same
  norm as the sensitivity metrics.

**Pass/fail interpretation.**
- Pass: (a)--(c) hold inside declared validity regions with documented
  tolerances; near rank events the breakdown is visible and correctly flagged
  rather than silently misreported; any genuinely certified Lipschitz bound is
  never crossed, while empirical bounds are labeled as such; frozen-nuisance
  error quantified.
- Fail: first-order predictions fail at $O(\varepsilon)$ instead of
  $O(\varepsilon^2)$ inside the declared gap region; surrogate exceeds sampled
  worst case at first order; a genuinely certified uniform bound is crossed
  beyond numerics; subspace tracking
  "jumps" while rank and gap are still fixed.

**What it cannot establish.** Global robust-trajectory optimality outside the
sampled trajectory family and uncertainty model; combined inner-prior and
outer-bounded-mismatch regimes without an explicit separation design;
non-Gaussian or map-correlated pose-error behavior; online SLAM performance;
real multipath/cycle-skip regimes.

---

## Cross-family discipline

- Every numeric check is tied to a specific claim and status in `seed.md`;
  passing a local check never upgrades an [open] claim to [accepted].
- Report acceptance gates separately: fixture-level experiment passes are not
  production or global-novelty acceptance.
- All runs reproducible from config files; record exact executed commands,
  versions, environment notes (Apple Silicon CPU), and artifact hashes.

---

# Open Questions Register

> Separated by severity and type. (a) questions whose "no" would threaten core
> validity; (b) hard theory/proof questions suitable for later GPT Pro work;
> (c) design/construction questions; (d) ordinary pipeline-solvable tasks.
> Items are carried over from the source discussions with their original open
> status; nothing here is answered in this seed.

## (a) Fatal core-validity questions

These decide whether the research program, as framed, can stand. Codex
intervention is required for any negative adjudication.

1. **Is the local tangent-intersection framing the correct core object for the
   intended claim?** The whole program rests on
   $\delta y = A\delta\chi + B\delta X + n$ near a nominal state. If the
   dominant SLAM failure in real operation is non-local (cycle skipping, wrong
   data association, phase ambiguity), can local retention statements still
   support the intended conclusions, and under which explicitly scoped claims?
2. **Born empty-background / weak-scatterer regime.** At $\chi_0=0$ the
   first-order pose Jacobian vanishes and the leading mismatch is the bilinear
   term $(D_XA[\delta X])\delta\chi$. Does the linearized retention framework
   govern any practical low-contrast operating regime, or is that regime a
   genuinely second-order sensing-operator-uncertainty problem requiring a
   different theory? A "no" here must shrink or relocate the theory's scope.
3. **Does the nuisance model remain faithful when pose dimension grows and
   uncertainties couple?** Per-pose nuisance rank grows as $d_xT$; inner
   $J_X$ (to-be-estimated pose) and outer $\Delta X$ (execution/model
   mismatch) may describe the same physical source. If correlated map--pose
   errors, non-Gaussian pose distributions, or repeated relinearization break
   the Woodbury/Schur equivalence, which accepted formulas no longer correspond
   to any real estimation setting?
4. **Gauge quotient vs bounded-domain reality.** Exact $SE(2)/SE(3)$ invariance
   requires an unbounded/uniform setting. Real experiments have boundaries,
   anchors, and known background that break exact gauge. If exact gauge is
   absent, do $\rho=0$ directions still occur generically as local coincidences,
   and does every claim that leans on gauge still transfer?
5. **Discretization fidelity.** If the infinite-dimensional problem has
   non-closed ranges or accumulating singular values, numerical identities may
   be artifacts of a finite grid. Is stable transversality
   $\inf_{Au\ne0}\mathrm{dist}(Au,\mathcal V)/\|Au\|>0$ verified for the
   intended continuum model, and are finite-dimensional results transferable?
6. **What does local retention actually predict for system success?** If
   $\rho\approx1$ but the basin of attraction is tiny, the metric may be
   locally true and globally misleading. Is a basin/cycle-skip statement needed
   before any practical claim can be made?

## (b) Hard theory/proof questions (later GPT Pro work)

1. Infinite-dimensional operator version: closed-range conditions, polar
   decomposition $A=U_A|A|$, $R_{\mathrm{geom}} = U_A^*P_{B^\perp}U_A$ on
   $\overline{\operatorname{Range}(A^*)}$, frame/polar formulations, and the
   stable-transversality condition.
2. Continuous Helmholtz $SE(2)/SE(3)$ gauge theorem: rigorous group actions on
   contrast and trajectory, Sobolev/domain structure, boundary/anchor/known-
   background symmetry breaking, gauge quotient FIM.
3. Born far-field Fourier theorems: explicit transversality and missing-wedge
   structure for line/circle/partial-aperture trajectories; single-frequency
   near-confounded mode characterization; lower bound on $\theta_{\min}$ from
   bandwidth/aperture.
4. Multi-frequency transversality theorem: under which angle/bandwidth/MIMO
   conditions are map and pose tangents generically separated except for global
   gauge; rigorous version of the shared-compensation kernel condition.
5. Full-wave resolvent perturbation bounds: explicit control of
   $\|DA\|,\|DB\|,\|DK_{\mathrm{eff}}\|$ by $\|M^{-1}\|$, $k$, $R$, $\|\chi\|$,
   domain size, and noise; linearization-radius bounds; whether resonance
   information gain can outweigh sensitivity growth.
6. Prior-weighted spectrum theory: bounds, monotonicity, and interpretation of
   $R_X = Q_A^*(I+BJ_X^{-1}B^*)^{-1}Q_A$; generalized Schur support with
   singular $J_0$ and block-banded motion priors; interaction of motion-graph
   gauges with map gauges.
7. Rank-event and continuity theory: smoothness of defect projectors on fixed-
   rank strata, eigenvalue crossing and cluster-merge update rules, gap
   stability, eigen/subspace derivative bounds with degenerate eigenvalues via
   compressed operators, and validity boundaries for the robust first-order
   surrogate at crossings and under cone-constrained perturbations.
8. Spectral classification and projector-composition theory: threshold/gap
   stability for $\mathcal E_G/\mathcal E_M/\mathcal E_R$; exact conditions
   under which hierarchical $S$/$D$/pose composition is valid and the octant
   decomposition is not; approximation error of intersection-by-projection
   constructions ($P_{S^-}V_D^+$ style).

## (c) Design/construction questions

1. Two-fold SOM exactness: which operator defines $V_D^\pm$ in the intended Chen
   source/version, which projection order and optimization variables are used,
   and whether the pose-defect layer is best placed before or after the
   $S$/$D$ reduction.
2. Fourier/FFT surrogate construction: which basis and frequency selection
   approximate full-wave $V_D^+$ with a bounded Grassmann distance, and whether
   approximate intersection bases are acceptable when projectors do not
   commute.
3. Map parameterization: pixel vs smooth/Fourier/spline basis, the pullback
   $T_\chi = \partial J/\partial\chi$, and which representation keeps gauge
   checks, rank identities, and Born comparisons well-posed.
4. Full model of $B$: world-fixed vs robot-local grids; inclusion of
   $D_xG_D$, $D_xD_\chi$, direct paths, antenna patterns, phase centers, clock,
   and channel gains; how calibration blocks layer on pose blocks and change
   defect rank.
5. Robust trajectory formulation: uncertainty set, norm/dual norm, motion and
   collision constraints/tangent cones, and a CPU-only numerical scheme for
   $\max_X\min_{\|\Delta X\|\le\varepsilon}\lambda_r(K_{\mathrm{eff}})$.
6. Continuous reporting conventions: whitening metric, realification, rank
   tolerance, observable-support definition, principal-angle ordering, and the
   threshold-free metrics used in every artifact.
7. Incremental tracking design: fixed-rank stratum update, rank-event triggers,
   hysteresis windows, and the nested-approximation stopping criterion from the
   small-defect tail formula.

## (d) Ordinary pipeline-solvable tasks

1. Implement the 2D Helmholtz contrast-source forward solver with Born switch,
   $M$ solve/invertibility reports, and analytic $A$/$B$ builders used by the
   experiment families.
2. Build the harness: four equal-budget trajectories, frequency sets (genuine
   vs duplicate control), colored-noise whitening, smooth/pixel basis variants,
   and fixed seeds.
3. Implement the spectral/algebra library: rank-revealing QR/SVD, principal
   angles/CS decomposition, generalized eigenvalue solve restricted to the
   observable support, PSD-ordering checks with scaled tolerances.
4. Run and report prior sweeps, perturbation sweeps, worst-case boundary
   sampling, Lipschitz estimation, and trajectory/frequency comparisons with
   the exact executed checks recorded.
5. Reproducibility packaging: config schemas, versioned outputs, command logs,
   environment notes (Apple Silicon CPU), artifact hashes; no new
   orchestration/daemon/GPU infrastructure.

## Routing note

(a) stays with Codex (core-validity judgment). (b) is offered to later GPT Pro /
deep-theory sessions only after (a) is answered and the algebraic spine is
numerically verified. (c) is resolved by Codex-in-the-loop design review before
building. (d) can be executed by ordinary workers directly.
