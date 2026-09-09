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
