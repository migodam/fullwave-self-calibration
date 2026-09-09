# Five minimal falsification experiments (existing Python harness)

Harness location (read-only reuse):
`experiments/idea_loops/loop_2026-09-03_16-25-02/experiment_pose_confounding_spectral_geometry/`
— `src/helmholtz.py` (`make_grid`, `build_operators`, `build_AB`,
`whiten_realify`, `born_forward`, `forward_measurements`, `pose_geometry`,
`green_grad_first`, `green_grad_source`), `src/family1_pilot.py`
(`build_poses`, `make_chi0`), `src/family2_algebraic_spine.py`
(`build_smooth_basis`), run with `.venv/bin/python` from the experiment root
and `MPLCONFIGDIR=/tmp/mpl` to avoid the font-cache warning.

Each experiment states the claim it targets, the exact computation, and the
**falsification gate** (what observed value would kill the claim). All are
small (N=8–16, T≤3, dense complex128); each should run in seconds to a few
minutes. Two of them (E1, E2) were already executed in this review session
(results cited); E4/E5 partially overlap existing `family3_gauge_born.py` and
`family8_som_distinction.py` runs, noted below.

---

## E1. Decomposition identity + single-snapshot degeneracy

Target claim: (a) $B_l = dG_S[\delta X_l]j + G_S M^{-1}D_\chi
d_Xe^{\mathrm{inc}}[\delta X_l]$ holds exactly; (b) at finite $m\le S$ with
full row rank, every pose-induced sensing perturbation is exactly an
equivalent current (residual 0).

Config: N=8, T=1, 4 Rx (`±0.06` cross), Tx at platform center, Family-1 arc
pose, `make_chi0` contrast.

Compute, per pose coordinate $l\in\{p_x,p_y,\theta\}$:

1. `build_AB` gives $B_l$; rebuild the two terms from `build_operators` +
   explicit Green gradients (as in `build_AB`'s own split).
2. Gate: $\|B_l-\mathrm{sens}-\mathrm{int}\| < 10^{-14}$.
3. Finite-difference cross-check of $B$ vs
   $\frac{F(\chi,X+\delta X)-F(\chi,X-\delta X)}{2\delta X}$
   (gate: relative error $<10^{-6}$).
4. SVD of $G_S$: report $\operatorname{rank}(G_S)/m$ and residual
   $\|P_{\ker G_S^*}\mathrm{sens}\|$; gate for (b): $<10^{-12}$ at full rank,
   while $\|G_S G_S^\dagger\mathrm{sens}-\mathrm{sens}\|=0$ and
   $\|G_S(\delta j_X+G_S^\dagger\mathrm{sens})-B_l\|=0$.
5. Record $\|\delta j_{\mathrm{eq}}\|$ and its location: fraction of norm in
   $V_S^+$ vs $V_S^-$ (foreshadows E5).

Falsifies: (a) any hand-waved "$B$ is not a current object"; (b) the claim
that $V_P$ is a *proper* subspace per snapshot — expected result is residual
$\approx0$ and $V_P$ effectively spanning all visible data directions.

Session evidence: all gates already passed at N=8 (identity $2.6\times10^{-18}$,
residuals $<5\times10^{-18}$, equivalent current reproduces $B_l$ to
$9\times10^{-16}$; rotation column has zero interior part). Codify as a
repeatable `family16_*` script for provenance.

## E2. Rank-deficient leakage: rigid motion vs per-channel calibration

Target claim: the equivalent-current range condition is only exercised by
perturbations that break the null-row symmetries of $G_S$.

Config: N=8, T=1, receivers `[[-0.06,0],[-0.06,0],[0,-0.06],[0,0.06]]`
(pair 0–1 coincident ⇒ rank 3/4; row-nullspace $\ker G_S^*$ is 1-dimensional
and spanned by $e_0-e_1$; $\ker G_S$ is $(S-3)$-dimensional — pull any single
null vector for the test).

Compute:

1. Confirm $\operatorname{rank} G_S=3$ and find $\ker G_S$ via right singular
   vectors; check $\|G_S j_0\|<10^{-12}$ for $j_0\in\ker G_S$.
2. Body-fixed rigid perturbation: $v=dG_S[\delta p_x]j_0$ by central
   differences of $G_S(X\pm\delta p_x)$; gate: residual
   $\|P_{\ker G_S^*}v\|/\|v\|<10^{-10}$ (structure preservation, C2).
3. Per-channel phase calibration: $v=i\delta\varphi\,e_0\,(G_Sj)_0$ (row 0);
   gate: residual ratio $\approx1/\sqrt2>10^{-3}$ (C3). Repeat with an
   all-channel phase vector (expected $\approx0.89$) and a per-channel gain
   vector.
4. Contrast check: interior-path image always has zero residual regardless of
   rank deficiency (gate $<10^{-12}$), proving the interior/sensing residual
   asymmetry at the operator level.

Falsifies: "pose perturbations are always representable as equivalent
currents" (they are not, once calibration DOF are included) and the naive
"rank-deficient ⇒ nonzero residual" claim (rigid motion keeps it zero).

Session evidence: rigid-motion residual ratio $1.3\times10^{-15}$; phase-cal
ratio $7.07\times10^{-1}$ single channel, $8.87\times10^{-1}$ all-channel;
interior residual $\approx7\times10^{-18}$. Codify and extend with a
$\delta\varphi$ sweep to show the ratio is geometry-independent.

## E3. Interior/sensing asymmetry: resolvent amplification and Born collapse

Target claim: the interior path is amplified by $M^{-1}$ (diverges near
resonance) while the sensing path is not; the two must not share one
robustness constant.

Config: N=8, T=1, standard 4-Rx cross. Sweep $k_b$ over a range containing a
resonance of $M=I-D_\chi G_D$ (locate near-resonance by minimizing
$\sigma_{\min}(M)$, e.g. via a coarse $k_b$ scan) and sweep contrast amplitude
$\alpha\chi_0$, $\alpha\in\{0.1,0.5,1\}$.

Compute at each $(k_b,\alpha)$:

1. $\|\mathrm{sens}\|$, $\|\mathrm{int}\|$, ratio
   $\|\mathrm{int}\|/\|\mathrm{sens}\|$ for $\delta p_x$ and $\delta p_y$.
2. $\sigma_{\min}(M)$, $\|M^{-1}\|$.
3. Born switch (`born_forward`): the Born data has no $M^{-1}$; recompute the
   ratio and confirm it collapses to the $M^{-1}$-free baseline as
   $\alpha\to0$.
4. Gate (falsification of "symmetric paths"): the ratio is not constant and
   grows with $\|M^{-1}\|$; if it stays $\approx1$ across the resonance scan
   the claimed asymmetry is numerically unsupported in this model.

Falsifies: treating interior and sensing pose effects as one combined
"geometry" perturbation with a single magnitude budget.

## E4. Gauge counterfeiting and pose-only stabilizer (current level)

Target claim: the $SE(2)$ gauge pair is an exact equivalent-current pair
($B\xi_X=-A\xi_\chi$, zero residual, and the equivalent current is a physical
map change), and radially symmetric scenes give pose-only null directions.

Config: N=16 (smooth basis needed for generators), T=3 poses, standard Rx
cross; map basis from `build_smooth_basis` so the translation/rotation
generators $\xi_\chi=-(v+\Omega x)\cdot\nabla\chi$ remain in span.

Compute:

1. $A,B$ (whitened, realified); SE(2) generators; gate
   $\|A\xi_\chi+B\xi_X\|/(\|A\xi_\chi\|+\|B\xi_X\|)<10^{-8}$ (already the
   `family3_gauge_born.py` check — reuse as control).
2. Current-level identity: $T_\chi=M^{-1}\operatorname{diag}(E^{\mathrm{tot}})$;
   gate $\|G_S(-T_\chi\xi_\chi)-B\xi_X\|/\|B\xi_X\|<10^{-8}$ ⇒ the pose gauge
   data **is** the image of a physical map current change (C5).
3. Radial scene: in this harness $\theta$ rotates the body about the platform
   position $p$, so place $\chi=\varphi(\|x-p\|)$ centered at $p$ with the Tx
   at the platform center; gate $\|B_{\theta}\|\approx0$ (rotation column)
   while translation columns stay nonzero (C6). Repeat with the off-center
   `make_chi0` blobs as the non-symmetric control.
4. Report which retained/current subspace $V_P$ basis vector aligns with
   $-T_\chi\xi_\chi$ (angle between them).

Falsifies: any TriSpace labeling that claims pose and map confounded
directions can be separated without an anchor; and "pose directions always
have data effect".

## E5. Retained-subspace leakage, non-commutativity, rank event

Target claim: pose directions are representable inside the retained TSOM
basis ($r_{\mathrm{SOM}}=0$), and S-fold/pose-fold commute.

Config: N=8–16, T=2–3 poses (stacked, shared pose DOF), standard Rx cross.

Compute:

1. Stacked $G_S$ = block diagonal; SVD; gap-selected $V_S^+$ at rank $r$;
   TSOM retained basis $U_T=[V_S^+, \mathrm{QR}(P_{S^-}V_D^+)]$ with $V_D^+$
   from $G_D$'s SVD (or FFT surrogate per the context file).
2. $r_{\mathrm{SOM},l}=\|P_{(\operatorname{Ran} G_SP_T)^\perp}B_l\|$ for each
   pose coordinate, whitened/realified; gate: generically $>0$ while the
   full-space residual $r_S$ from E1 is $0$ — the falsifying observation for
   "retained basis suffices".
3. Non-commutativity control (report Counterexample 4.73 in the harness):
   dominant current mode of $G_S^*G_S$ vs of $G_S^*WG_S$
   ($W=I-P_B$ or finite-prior contraction); report the angle; gate: nonzero.
4. Leakage metric: perturb pose to $\hat X=X+\delta X$; recompute $V_S^+(\hat
   X)$; report $\|P_S^-(X^\star)P_S^+(\hat X)\|$ as $\|\delta X\|$ sweeps
   upward, and the smallest distance of a singular value of $G_S(\hat X)$ to
   the truncation threshold; find the $\|\delta X\|$ where a rank event
   occurs and record the projector jump (expect $\approx1$ at the crossing,
   report Theorem 4.48).
5. Contrast with `family8_som_distinction.py` outputs (it verified
   $V_S\neq V_A$ and $K_{\mathrm{eff}}$ eigenvector overlap, but not the
   $V_P$/retained-basis residual).

Falsifies: "the pose fold can be analyzed inside $V_S^+$" and any exact
"eight-cell" decomposition built from projector products.

---

## Output contract

Each experiment writes results to
`results/` and a short note to `notes/` under the experiment root, following
the existing family convention; all gates recorded with relative errors and
machine-precision rank tolerances
$\tau_{\mathrm{rank}}=\max(m,n)\epsilon_{\mathrm{mach}}\sigma_1$, never a
physical SOM threshold. Green checks are bounded finite-dimensional
verifications only — none of them establish a continuum theorem, and none
establish novelty.
