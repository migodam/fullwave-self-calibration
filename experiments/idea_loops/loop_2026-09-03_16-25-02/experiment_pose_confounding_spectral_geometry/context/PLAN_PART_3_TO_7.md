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
