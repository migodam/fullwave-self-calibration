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
