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
