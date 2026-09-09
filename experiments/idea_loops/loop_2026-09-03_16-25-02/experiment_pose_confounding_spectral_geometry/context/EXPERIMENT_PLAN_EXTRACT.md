# EXPERIMENT_PLAN_EXTRACT

Read-only extraction (2026-09-03; no code run) from all six listed `context/`
files. Status: **[accepted]** derived under stated conditions;
**[conditional]** only under listed conditions; **[open]** conjecture;
**[non-claim]** excluded. Operators live in the whitened/realified data space;
`U=Range(A)`, `V=Range(B)`; `+` = Moore--Penrose.

# 1. Five experiment families

Shared scenario: `D=[-0.5,0.5]^2`; `32x32` default + `40x40` sensitivity grid;
nonzero two-blob `chi_0`; smooth basis for gauge; 1 Tx + 4-8 body-fixed Rx per
pose `x_t=(p_x,p_y,theta)`; line / 90° / 180° / 360° trajectories, equal path
length/count; Born/full-wave; genuine `F in {1,2,3}` + duplicate control
`(A_2,B_2)=c(A_1,B_1)`; identity + one colored noise; realify after whitening;
report `sigma_min(M_t)/||M_t||`. Exact configs/tolerances: see section 7.

## Family 1 -- Forward model / analytic Jacobians: derivative consistency
Claims: conditional `A_t`,`B_t`; conditional Born/full-wave regime; resolvent
sensitivity. Checks/metrics: (1) centered FD map check
`D_hF_t[dchi]=[F_t(chi_0+h dchi,X_0)-F_t(chi_0-h dchi,X_0)]/(2h)` vs
`A_t dchi`: `max_t||D_hF_t[dchi]-A_t dchi||/max(||A_t dchi||,eps)`, random
unit `dchi` (several seeds); analogous pose check vs `B_t dX`. (2) FD sweep
`h in 10^-4..10^-1` relative; convergence order ~2 until roundoff floor.
(3) state residual; `sigma_min(M)/||M||`; max field. (4) Born-vs-full-wave
discrepancy vs contrast/spectral knob (blob amplitudes), weak-scattering region
separated from resonance/cancellation.
Pass = errors within declared tolerance (target `10^-5..10^-3` by
h/grid/conditioning) with consistent `h` scaling; Born/full-wave discrepancy ->
0 toward weak scattering; `M` away from resonance floor. Fail = large
`h`-insensitive mismatch; missing formula term (correct only with term +
updated conditions); discrepancy fails to vanish; singular `M` unexplained.

## Family 2 -- Algebraic spine: kernel/rank/retention, ordering, prior sweeps
Claims: accepted identities 1-9; conditional prior-limit/interlacing.
Checks/metrics: (1) `ker K_SLAM={u:Au in V}` (both inclusions; includes
`ker A`). (2) `rank K_IS-rank K_SLAM=dim(U∩V)`. (3)
`eig(R_op)=1-sigma_i^2(Q_B^*Q_A)` with `R_op=K_IS^{+/2}K_SLAM K_IS^{+/2}=
V(I-Q^*ZZ^*Q)V^*` (thin SVD `A=Q Sigma V^*`, orthonormal `Z` of `V`);
`0<=rho_i<=1`; at most `rank B` differ from 1. (4) `rank L_X<=rank B`
(`L_0=A^*P_VA`; `L_X=K_IS-K_eff`). (5) Interlacing
`lambda_{i+p}(K_IS)<=lambda_i(K_eff)<=lambda_i(K_IS)` for `rank L_X<=p`;
ordering `K_SLAM<=K_eff(alpha)<=K_IS`. (6) Sweep `J_X=alpha I`, `alpha in
10^-6..10^4` log-spaced; both limits on regular path + one singular-support
`J_X`. (7) PSD order violations: min eigenvalue `>= -delta`
(backward-error-scaled). (8) no-prior duplicate leaves `rho` unchanged (F4).
Additional metrics: max scaled deviations
(`||K_SLAM-A^*P_VA||/||A^*A||`, `||R_op-V(I-Q^*ZZ^*Q)V^*||`), rank
differences, count `rho_i!=1`, spectrum range/gap, ordinary-vs-generalized
eigenvector rotation stats, limit curves, `d_conf`, `prod rho_i`.
Pass = identities hold to conditioning-scaled accuracy; no ordering violation;
limits/counts reproduce. Fail = any violation beyond tolerance; interior
ordering violation; limit failure explainable only by singular-support prior
(caveat confirmed).

## Family 3 -- Gauge and Born empty-background degeneration
Claims: accepted gauge item 10 and Born item 12; conditional basis caveats.
Checks/metrics: (a) smooth-basis `SE(2)` generators (Tx,Ty,rotation): gauge
residual `||A dchi_g+B dX_g||/(||A dchi_g||+||B dX_g||)≈0`; generator-map
distance to `ker K_SLAM`; gauge rank/kernel dim per scene variant.
(b) Born `chi_0=0`: `B=0`, `K_SLAM=K_IS`, bilinear dominance of
`(D_XA[dX])dchi`; compare its norm vs `||A(X_0)dchi||` and vs remainder
`||F(dchi,X_0+dX)-A(X_0)dchi-(D_XA[dX])dchi||` over shrinking scales; observed
exponent in `dX`. (c) anchors/known background/fixed boundaries remove expected
generators; `rho_min` before/after; symmetric-scene variant.
Pass = residuals reach predeclared floor from independent forward/Jacobian/
discretization errors and decrease under refinement (unless exact discrete
equivariance gives roundoff floor); pixel basis shows documented failure;
empty-background/anchor predictions hold. Fail = residuals above floor in
smooth basis; bilinear term not dominant; anchors change nothing. Guard:
converse (`rho=0` => gauge) is [open]; bounded scenes lack exact global `SE(2)`.

## Family 4 -- Frequency stacking and trajectory geometry
Claims: accepted multi-frequency items 11-12; open generic-transversality
(sampled only). Checks/metrics: (a) `K_eff^(1:F+1)>=K_eff^(1:F)` (PSD),
distinct `F=1..3`, fixed prior and `J_X=0`. (b) duplicate control under
`J_X=0`, fixed finite `J_X`, finite `J_X` scaled by `1+|c|^2`; invariance only
in settings 1 and 3; `K^(new)=(1+|c|^2)K^(old)`. (c) shared-`z` kernel:
single-frequency-confounded `u` survives multi-frequency only if one `z` fits
`A_fu=B_f z` for all `f`; quantify `rho` movement, diverse vs duplicate.
(d) trajectory ranking by retention/defect stats, hypothesis `360° > 180° arc >
90° arc > straight line` (counterexample acceptable). Metrics: absolute spectra
per F; PSD checks; `{rho_i}`, `theta_min`, `d_conf`, defect rank, gauge-kernel
dim; multi-kernel survival; ranking.
Pass = monotonicity; qualified duplicate predictions; at least some confounded
directions separated by diverse frequencies (else significant negative result);
trajectory order reproduced or counterexample. Fail = stacking decreases
`K_eff`; duplicate alters `rho` with `J_X=0` or fails matching scaling; shared-
compensation contradicted internally.

## Family 5 -- Trajectory sensitivity, robust surrogate, rank-event stability
Claims: conditional derivative formulas; robust surrogate; open
continuity/rank-event items; frozen-nuisance approximation. Checks/metrics:
(a) simple eigenvalues with positive gap: `lambda_i(X+dX)-lambda_i(X)` vs
`v_i^*DK_eff[dX]v_i` on `||dX||=epsilon` spheres (fixed seeds, several
`epsilon`); error `O(epsilon^2)` in gap regime. (b) full `DK_eff` (with
`dot B` via pose Hessians) vs frozen-nuisance variant, both vs finite
differences. (c) robust surrogate on full norm balls vs sampled worst case;
explicit fixed norm/dual norm. (d) certified bound only with derived uniform
operator-Lipschitz `L_K` over whole ball: `lambda_r(X+dX)>=lambda_r(X)-
L_K||dX||` (Weyl form needs no simple eigenvalue); nominal derivative norm or
sampled max = empirical only; report empirical violations as such. (e)
`||U(X+dX)U(X+dX)^*-U(X)U(X)^*||`, spectral gap, `sigma_min(B)` along paths;
engineer near-crossing/rank events.
Pass = (a)-(c) hold in declared validity regions; rank-event breakdown
visible/flagged; certified bounds never crossed; empirical bounds labeled;
frozen-nuisance error quantified. Fail = `O(epsilon)` errors inside gap regime;
surrogate exceeds sampled worst case; certified bound crossed; tracking jumps
with rank/gap fixed.

# 2. Core finite-dimensional theorem list

(a) `K_IS=A^*A`; `K_SLAM=A^*P_{Range(B)^perp}A` (`P_{Range(B)^perp}=I-BB^+`);
`K_eff=A^*A-A^*B(B^*B+J_X)^+B^*A`, `J_X>=0`;
`L_X=K_IS-K_eff=A^*B(B^*B+J_X)^+B^*A=H_X^*H_X`,
`H_X=(B^*B+J_X)^{+/2}B^*A`.

(b) `0<=K_SLAM<=K_IS`, so `lambda_i(K_SLAM)<=lambda_i(K_IS)`. **[accepted]**

(c) `ker K_SLAM={u:Au in V}` (includes `ker A`). **[accepted]**

(d) `rank K_IS-rank K_SLAM=dim(U∩V)`; **[accepted]** `U∩V={0}` iff
`rank[A B]=rank A+rank B`; `dim(U∩V)>=max(0,r_A+r_B-m)` (`m` realified data
dimension).

(e) `rho_i=sin^2(theta_i)`: on `S=Range(A^*)=ker(A)^perp` solve `K_SLAM v=
rho K_IS v` (`v^*K_ISv=1`) or diagonalize `R_op=K_IS^{+/2}K_SLAM
K_IS^{+/2}=V(Q^*P_{V^perp}Q)V^*`; orthonormal `Z` of `V` gives
`sigma_i(Z^*Q)=cos(theta_i)`, so `eig(R_op)=1-sigma_i^2(Z^*Q)=sin^2(theta_i)`,
`0<=rho_i<=1`. **[accepted]** Conditions: whitened/realified finite-dim, no
pose prior; pad extra angles to `pi/2` if `rank A>rank B`; at most `rank B`
differ from 1; ordinary eigenvectors may rotate broadly. Finite prior gives
weighted shrinkage: `D=BJ_X^{-1/2}`,
`I-B(B^*B+J_X)^{-1}B^*=(I+DD^*)^{-1}`, `R_X=Q_A^*(I+DD^*)^{-1}Q_A`, spectrum
in `(0,1]`, generally not `sin^2`. **[conditional]**

(f) `rank L_X<=rank B` **[accepted]**; `tr L_X=||C^{+/2}B^*A||_F^2`,
`lambda_max(L_X)=||C^{+/2}B^*A||_2^2` (`C=B^*B+J_X`); affected
(`#{theta_i<pi/2}`) vs destroyed (`#{theta_i=0}`) DoF are distinct.

(g) Along fixed-rank regular path `J_X(alpha)=alpha I` with support covering
every data-coupled pose direction: `alpha->0` recovers `K_SLAM`,
`alpha->infinity` recovers `K_IS`. **[conditional]** Singular/unanchored graph
gauges never converge to `K_IS`.

(h) Duplicate `(A_2,B_2)=c(A_1,B_1)` scales `K_IS`,`K_SLAM` by `1+|c|^2` and
leaves every `rho_i` unchanged **[accepted] only no-prior or jointly scaled
prior** (prior scaled by `1+|c|^2`). Caveat (req. 13): fixed finite `J_X`
changes data-to-prior weighting, so normalized retention generally changes.

# 3. Realification / whitening recipe (A1)

For `n~CN(0,Sigma)`: whiten `A_hat=Sigma^{-1/2}A`, `B_hat=Sigma^{-1/2}B`.
Real `chi,X`: `A_R=sqrt(2)[Re A_hat; Im A_hat]`, `B_R=sqrt(2)[Re B_hat;
Im B_hat]`; real FIM `I=[[A_R^T A_R,A_R^T B_R],[B_R^T A_R,B_R^T B_R]]`.
Complex contrast = two real blocks: `A_R=sqrt(2)[[Re A_hat,-Im A_hat],
[Im A_hat,Re A_hat]]`. (Real-param/complex-data FIM note: `2 Re(A^*
Sigma^{-1}A)`.) Compute all projectors/ranges/ranks/angles/Schur complements/
spectra after whitening then realification. Forbidden: complex-linear `BB^+`
for real pose increments.

# 4. Explicit 2D discretization requirements (rule 15 + theory)

Sources give the continuum conventions; the implementer must document
Green-normalization, cell area, singular self-cell quadrature/regularization,
incident convention, contrast units, and include one refinement diagnostic
(`32x32` vs `40x40`). No numeric self-cell rule appears in the sources.

- `D=[-0.5,0.5]^2`, world-fixed uniform `NxN` grid, `N in {32,40}`, step
  `h=1/N`, cell area `h^2` (declare explicitly).
- `g_k(r,r')=(i/4)H_0^{(1)}(k||r-r'||)`; gradient `grad_r g_k=-(ik/4)
  H_1^{(1)}(kR)(r-r')/R`, `R=||r-r'||`.
- `chi=(k^2-k_b^2)/k_b^2`; `G_D` domain propagator, `G_S` domain->receivers;
  `M_t=I-D_chi G_{D,t}`; `E^sca=G_SJ`; `J=D_chi(E^inc+G_DJ)`;
  `F(chi)=G_SM^{-1}D_chi E^inc`.
- Incident point source `[e_t^inc]_n=q_t g_k(z_n,s_t)`; antennas
  `r_{t,a}=p_t+R(theta_t)s_a`, `Dr_{t,a}[dp,dtheta]=dp+dtheta R(theta)J s_a`
  (`J=[[0,-1],[1,0]]`); `[D_xG_S[h]]_{a,n}=grad_1 g_k(r_{t,a},z_n)^T
  Dr_{t,a}[h]`; `[D_xe_t^inc[h]]_n=q_t grad_2 g_k(z_n,s_t)^T D s_t[h]`.
- (22) `A_t=G_{S,t}M_t^{-1}diag(E_t^tot)` (`E_t^tot=e_t^inc+G_{D,t}j_t`;
  basis map uses `diag(E_t^tot)S_chi`). (23) `B_t h=D_{x_t}G_{S,t}[h]j_t+
  G_{S,t}M_t^{-1}D_chi(D_{x_t}e_t^inc[h]+D_{x_t}G_{D,t}[h]j_t)`.
- World-fixed grid: `D_xG_D=0` (also no moving `D_chi`, no direct/calibration
  path) -> `B_th=D_xG_S[h]j_t+G_SM^{-1}D_chi D_xe^inc[h]`.
- Bounds `||A_t||<=||G_{S,t}||||M_t^{-1}||||E_t^tot||_inf`; analogous `B`;
  `sigma_min(M)->0` grows info/sensitivity, shrinks linearization radius.
- Born `A(X)=G_S(X)diag(e^inc(X))`; validity = weak internal multiple
  scattering (e.g. `||D_chiG_D||<1` suitably normed), not merely low contrast.

# 5. Gauge / SE(2) generators and Born-empty-background material

Gauge cancellation check: for each generator pair `(dchi_g,dX_g)`,
`||A dchi_g+B dX_g||/(||A dchi_g||+||B dX_g||)≈0`; floor = roundoff only under
exact discrete equivariance, else declared discretization/interpolation error
decreasing under refinement. Hence `A dchi_g=-B dX_g in V`, `dchi_g in
ker K_SLAM`, `rho=0`.

2D generators: translation `v`: `dchi_v(r)=-v^T grad chi(r)`, `dp_t=v`,
`dtheta_t=0`; rotation `omega`: `dchi_omega(r)=-omega(Jr)^T grad chi(r)`,
`dp_t=omega Jp_t`, `dtheta_t=omega`. At most 3 gauge directions; symmetry can
kill rotation; anchors/known background/boundaries break gauge. Controls:
smooth basis only; uniform unanchored / anchored / known-background /
symmetric variants; nonzero `chi_0` (gauge) vs `chi_0=0` (empty background).
Rules: global gauge => `rho=0` **[accepted]**; converse false **[open]**;
bounded boxed scene lacks exact global `SE(2)` **[non-claim]**.

Born empty background: nonzero `chi_0` -> `B_l=(dA(X)/dX_l)chi_0`,
`V=span{A_{X_1}chi_0,...,A_{X_p}chi_0}`; `chi_0=0` -> `B=0`, `K_SLAM=K_IS`,
expansion `F(dchi,X_0+dX)=A(X_0)dchi+(D_XA(X_0)[dX])dchi+O(||dchi||||dX||^2)`:
leading mismatch is bilinear `(D_XA[dX])dchi` (sensing-operator uncertainty;
do not claim pose errors harmless). Numeric check: bilinear norm vs
`||A(X_0)dchi||` and second-order remainder over shrinking scales; observed
`dX` exponent.

# 6. Robust sensitivity / first-order eigenvalue material

Simple eigenvalue: `delta lambda_i=v_i^*DK_eff[dX]v_i`, error `O(||dX||^2)`
in the declared gap radius; needs simple eigenvalue, positive gap, phase gauge
**[conditional]**. Eigenvector derivative (gauge `v_i^*dot v_i=0`):
`dot v_i=sum_{j!=i} v_j(v_j^*dot K_eff v_i)/(lambda_i-lambda_j)`;
`||dot v_i||<=||dot K_eff||/gap_i`. Normalized generalized retention
(`v_i^*K_ISv_i=1`): `dot rho_i=v_i^*(dot K_eff-rho_i dot K_IS)v_i`.

Full formulas (`C=B^*B+J_X` invertible or consistent convention):
`dot K_eff=dot A^*WA+A^*W dot A+A^*dot W A`, `W=I-BHB^*`, `H=C^{-1}`,
`dot W=-dot BHB^*-BH dot B^*+BH dot C H B^*`, `dot C=dot B^*B+B^*dot B+dot
J_X`; `dot B=D_X^2F[h,.]` needs the pose Hessian; dropping it =
frozen-nuisance-subspace approximation (name it, quantify error). No-prior
projector derivative (locally constant rank):
`dot P_B=P_{B^perp}dot B B^++(P_{B^perp}dot B B^+)^*`;
`||dot K_SLAM||<=2||A||||dot A||+2||A||^2||B^+||||dot B||`.

Worst case: `gamma_r(X,eps)=min_{||dX||<=eps}lambda_r(K_eff(X+dX))`; under
simple `lambda_r`, positive gap, C^2 smoothness, full norm ball:
`gamma_r=lambda_r(X)-eps||grad_X lambda_r||_*+O(eps^2)` (`||.||_*` dual norm;
cone-constrained case uses tangent-cone support function) **[conditional]**.

Certified vs empirical: certified requires a valid uniform operator-Lipschitz
`L_K` over the whole ball: `gamma_r>=lambda_r(K_eff(X))-eps L_K+O(eps^2)`
(Weyl-based, survives crossings, never violated beyond numerics, may be loose).
Nominal derivative norm or sampled maximum = local/empirical diagnostic only;
empirical violations are reported as such, never as theorem violations.

# 7. Explicit tolerances, seeds, configs, resolutions, options

- Domain/grid: `D=[-0.5,0.5]^2`; `32x32` default; `40x40` sensitivity run.
- `chi_0=0.3 chi_disk,1+0.5 chi_disk,2`; smooth basis for gauge.
- 1 Tx + 4-8 body-fixed Rx; line / 90° / 180° / 360°; equal path length and
  measurement count.
- Genuine `F in {1,2,3}`; duplicate `c`-scaled control; one shared `z`.
- Noise identity + colored; whiten then realify.
- `M_t` margin `sigma_min(M_t)/||M_t||`.
- Prior `J_X=alpha I`, `alpha in 10^-6..10^4` + singular-support `J_X`; PSD
  checks `>= -delta` (backward-error-scaled).
- Rank tolerance (exact checks): `tau=max(m,n) eps_mach sigma_1`.
- FD `h in 10^-4..10^-1`; centered `O(h^2)` to roundoff floor; Jacobian error
  target `10^-5..10^-3` by regime.
- Perturbations on `||dX||=epsilon` spheres, fixed seeds, several `epsilon`;
  explicit norm/dual norm.
- Any `L_K`: state derivation/coverage; otherwise empirical-only.
- Every claim: exact commands, versions, environment (Apple Silicon CPU),
  configs, raw tables, figures, artifact hashes; machine rank checks and
  physics thresholds never mixed.
