# CONTEXT_SUMMARY -- experiment_pose_confounding_spectral_geometry

Created: 2026-09-03. Working directory:
`/Volumes/migodam's-external-brain/Research/Inv_SLAM/experiments/idea_loops/loop_2026-09-03_16-25-02/experiment_pose_confounding_spectral_geometry`

Subfolders created: `context/`, `src/`, `results/`, `figures/`, `logs/`, `notes/`.
No solver, experiment, or implementation code has been written yet (setup/read-only phase).

## 0. Important discrepancy and resolution

The task instruction said workshop.md would name **three** source files. The actual
workshop.md names **five** files in its "Source priority" block. All five named
files were located, read, copied into `context/`, and summarized below so no
named source is silently omitted. The three core-theory files are the theory
context plus Q1 and A1; the other two are the delegated context-isolation summary
and the targeted prior-art addendum. All copies were verified byte-identical
(SHA-256 match) to their originals.

## 1. Required files read and copied

Each line: original absolute path -> `context/` copy (filename chosen to stay clear):

1. `/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/agentic_input/workshop.md` -> `context/workshop.md`
2. `/Volumes/migodam's-external-brain/Research/Inv_SLAM/Theory/SOM_SLAM_THEORY_CONTEXT.md` -> `context/SOM_SLAM_THEORY_CONTEXT.md`
3. `/Volumes/migodam's-external-brain/Research/Inv_SLAM/Theory/Questions/Q1.md` -> `context/Q1.md`
4. `/Volumes/migodam's-external-brain/Research/Inv_SLAM/Theory/Questions/A1.md` -> `context/A1.md`
5. `/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/context_isolation/summary.md` -> `context/context_isolation_summary.md`
6. `/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/literature/TARGETED_PRIOR_ART_ADDENDUM.md` -> `context/TARGETED_PRIOR_ART_ADDENDUM.md`

The relative paths in workshop.md ("Theory/...", "research/...") resolve against
the repository root `/Volumes/migodam's-external-brain/Research/Inv_SLAM/`.

## 2. Status-ledger rules from workshop.md (verbatim)

### Status legend (verbatim)

- **[accepted]** -- algebraically derived / adopted inside the prior discussions under the conditions stated in the ledger.
- **[conditional]** -- believed correct only under explicitly listed conditions, or only in one regime/model.
- **[open]** -- stated conjecture or unresolved question; not to be quoted as a result.
- **[non-claim]** -- deliberately excluded or explicitly not claimed.

### Non-negotiable execution requirements (verbatim)

1. Read and obey the complete workshop material below. It includes the theory ledger, exactly five experiment families, and open-question routing.
2. Implement one unified 2D scalar Helmholtz/contrast-source harness with Born/full-wave, frequency, trajectory, prior, and perturbation switches.
3. Attempt all five experiment families. A family may end in a documented negative result, corrected condition, or computational limitation; it may not be silently omitted. Run cheap falsifiers before expensive sweeps.
4. Use Apple Silicon CPU. Use MPS only if PyTorch supports an incidental task naturally. Do not install or assume CUDA, GPU Docker, a daemon, scheduler, database, web UI, or a new orchestration framework.
5. Keep current-space `G_S`/SOM modes distinct from the map Jacobian `A` and map-information eigenmodes. Keep deterministic nuisance pose, random pose marginalization, and bounded execution mismatch distinct.
6. Every numerical claim must have exact commands, deterministic seeds, tolerances, configs, raw result tables, and figures. Machine-precision rank checks and physical dominant-mode thresholds must never be mixed.
7. A passed finite-dimensional experiment supports only the tested model. It does not prove continuum transfer, global nonlinear convergence, practical online SLAM success, or global novelty.
8. Do not promote unproved continuum, transversality, rank-event, resonance, or robust-design claims into the theorem list. Put them in limitations or future work unless a valid proof is independently supplied.
9. Literature novelty is retrieval-bounded. Never write `first`, `no prior work`, or an exhaustive-novelty claim. Identify the closest overlap and separate mature algebra from the proposed physical synthesis.
10. Retain counterexamples and failed hypotheses. If a failure undermines the core thesis, prefer abandon or a narrowly justified revision over cosmetic metric changes.
11. Produce an English academic manuscript draft in LaTeX/PDF if the pipeline reaches writeup, with prominent disclosure that AI tools were used as required by the workflow license. The paper must distinguish proved finite-dimensional results, conditional derivations, observed numerical evidence, and open conjectures.
12. Hard new proofs and scientifically decisive theoretical conflicts remain for top-level Codex/GPT Pro review. The coding worker may test them or find counterexamples but must not declare them solved from numerical evidence.
13. The repeated-block retention invariance is a no-prior (or jointly scaled prior) statement. With fixed finite $J_X$, repeated data change the data-to-prior weighting and generally change normalized retention.
14. Do not headline the result as a direct "deformation of the SOM/TSOM subspace" unless an explicit $A=G_ST_\chi$ pullback/composition theorem is proved and implemented. By default, SOM is a current-space computational substrate and the pose-confounding spectrum is a distinct map-tangent/data-space layer; show their composition without identifying their modes.
15. Document the 2D Green-function normalization, cell area, singular self-cell quadrature/regularization, incident-field convention, and contrast units. Centered finite differences validate derivatives of the implemented solver, not the physical fidelity of an undocumented discretization; include at least one refinement diagnostic and scope the claims accordingly.

### Deliberate non-claims / guard list (verbatim bullets)

- SOM mode counts describe stable data-resolvable current directions, not intrinsic scene dimensionality; a 3D output does not imply every voxel is data-supported.
- "Theoretically unique", "large local Fisher information", "algorithm converges", and "output looks reasonable" are separate propositions.
- Regularization does not manufacture measurement Fisher information.
- ML modules are outside the core theory until objectives, equivariance, calibration, and data-consistency guarantees exist.
- No global novelty claim; closest-prior-work boundaries from the discussions are unverified and require a fresh systematic search before submission.

### Cross-family discipline (verbatim bullets)

- Every numeric check is tied to a specific claim and status in `seed.md`; passing a local check never upgrades an [open] claim to [accepted].
- Report acceptance gates separately: fixture-level experiment passes are not production or global-novelty acceptance.
- All runs reproducible from config files; record exact executed commands, versions, environment notes (Apple Silicon CPU), and artifact hashes.

## 3. Semantic-drift prohibitions from workshop.md

### 4.4 Semantic drift ban: current-space `G_S` vs map Jacobian `A` (mandatory, verbatim table)

| Object | Where it lives / what it means | What it must never be conflated with |
|---|---|---|
| `G_S` SVD, `V_S^±` | current space: which induced-current modes reach the array | map-parameter modes, `K_eff` eigenvectors, map observability spectrum |
| `A = D_chi F` (full wave) | `A = G_S M^-1 diag(E^tot) != G_S`; data-space tangent from map perturbations | `G_S` itself; `G_S`'s right singular subspace |
| Born `A` | `A ≈ G_S diag(E^inc)`; still a map Jacobian | `G_S` or current modes |
| `Range(A)` | data-space first-order map changes | `Range(G_S)`, MSR/MUSIC data subspaces, current span |
| `Range(B)` | data-space pose-induced changes | the pose coordinate space itself |
| `K_eff` eigenmodes | map/reduced-coordinate effective information after pose elimination | SOM current modes, `V_S^±`, plain per-mode eigenvalue ratios of `K_IS` |
| `rho_i` | Fisher-whitened principal-direction retention | `lambda_i(K_SLAM)/lambda_i(K_IS)`, or a per-`K_IS`-eigenmode ratio |

### Caveats (workshop section 4)

- Whitening/realification: all projectors and ranks are valid in the whitened space; after whitening with `Sigma_n^-1/2`, realify real physical parameters using `A_R = sqrt(2)[Re A; Im A]`, `B_R = sqrt(2)[Re B; Im B]`; complex contrast becomes real blocks `[Re A, -Im A; Im A, Re A]`. Using complex-linear `B B^+` for real pose increments is forbidden: it enlarges the nuisance range and overestimates map--pose confounding.
- Gauge discipline: state whether gauge is fixed, analysis is on the identifiable quotient/support, or Moore--Penrose is used without reading gauge variance as covariance. Exact gauge implies `rho=0`; the converse is not asserted.
- Born empty-background caveat: never write "Born at `chi_0=0` proves pose errors are harmless." The correct statement is that first-order pose Jacobian vanishes there, so pose-map confusion becomes a second-order bilinear mismatch problem with a different first-order operator.
- Five limits stay separate: uniqueness/identifiability, local observability, stability, statistical precision (FIM/CRB), and global multimodality/cycle-skipping are not interchangeable. Every claim must name which one it is about.

## 4. Notation definitions (as used in the source files)

Unified notation from the theory context, Q1, A1, and workshop:

- `D` / domain: `D=[-0.5,0.5]^2` (2D scalar Helmholtz, discretized, world-fixed map grid, uniform known background, static scene).
- `chi`, `D_chi = diag(chi)`: medium contrast and its discrete multiplication operator.
- `E^inc`, `E^tot`, `E^sca`: incident, total, and scattered field.
- `J`/`j_t`: contrast source / induced current.
- `G_S` (also `G_S(x_t)`): current-to-receiver data operator; "current space" object.
- `G_D`: internal/domain propagation operator.
- `M_t = I - D_chi G_D,t`: full-wave state operator (must be invertible at the nominal point; report `sigma_min(M)/||M||`).
- `F(chi,X)`: current-eliminated forward map `F(chi) = G_S M^-1 D_chi E^inc`.
- `x_t = (p_x,p_y,theta)`: robot pose at configuration t; `X=(x_1,...,x_T)` the trajectory.
- `A = D_chi F` (map Jacobian): full-wave `A_t = G_S,t M_t^-1 diag(E_t^tot)`; `A` is not `G_S`.
- `B = D_X F` (pose Jacobian): full form has `D_x G_S[h] j + G_S M^-1 D_chi(D_x E^inc[h] + D_x G_D[h] j) + (D_x D_chi[h])E^tot + direct terms`; simplified form used when `G_D`, map parameterization, direct/calibration terms do not move with pose.
- `K_IS = A*A` (known-pose inverse-scattering map information).
- `K_SLAM = A* P_{B^perp} A`, `P_{B^perp}=I-B B^+` (free deterministic pose, no prior).
- `K_eff = A*A - A*B (B*B + J_X)^+ B*A` (pose prior `J_X >= 0`, Schur complement).
- `K_post = K_eff + K_prior`; regularization does not create measurement Fisher information.
- `L_X = K_IS - K_eff = A*B(B*B+J_X)^+ B*A = H_X* H_X` (pose-induced information loss).
- `U = Range(A)`, `V = Range(B)` (data-space tangent subspaces).
- `rho_i`: generalized retention eigenvalues on the observable support of `K_IS`, i.e. solutions of `K_SLAM v = rho K_IS v`; `0 <= rho_i <= 1`; no-prior case `rho_i = sin^2(theta_i)` where `theta_i` are principal angles between `Range(A)` and `Range(B)` in whitened/realified data space.
- `R_geom`/`R_op = K_IS^{+/2} K_SLAM K_IS^{+/2} = Q_A* P_{B^perp} Q_A`; `rho` here means retention (an early "confounding amplitude = cos theta" convention was corrected).
- `R(u) = ||P_{B^perp} A u||^2 / ||A u||^2 = sin^2 angle(Au, V)` (directional retention).
- `d_conf = sum_i cos^2(theta_i) = ||Q_B* Q_A||_F^2 = tr(I - R_op)` (continuous confusion dimension).
- Destroyed DoF `#{theta_i=0}`, affected modes `#{theta_i<pi/2}` are distinct quantities.
- `J_X`: pose/motion prior information matrix (not to be confused with current `J`); `C_X=B*B+J_X`.
- `D = B J_X^{-1/2}` and prior-weighted `R_X = Q_A*(I + D D*)^{-1} Q_A`; eigenvalues in `(0,1]`, generally not `sin^2(theta_i)`.
- `G_S` SOM decomposition: `V_S^+` dominant/data-visible current directions, `V_S^-` complementary current directions (from SVD of `G_S`); `V_D^+`, `V_D^-` are the analogous internal/domain two-fold decomposition whose exact Chen-version definition is unresolved/open.
- `E_G`, `E_M`, `E_R`: gauge-like, mixed, robust spectral classes (threshold-dependent; require gap/stability reporting).
- `SE(2)/SE(3)` generators: translation `v`, rotation `omega`; `delta chi_v = -v^T grad chi`, `delta p_t = v`, `delta theta_t = 0`; `delta chi_omega = -omega (J r)^T grad chi`, `delta p_t = omega J p_t`, `delta theta_t = omega` (2D, with `J=[[0,-1],[1,0]]`).
- Green function (2D free space): `g_k(r,r') = (i/4) H_0^{(1)}(k||r-r'||)`; `grad_r g_k = -(ik/4) H_1^{(1)}(kR) (r-r')/R`.

## 5. Core claims / hypotheses present in the source files

Accepted/derived finite-dimensional statements (always under whitened, realified, locally linearized, gauge-appropriate, closed-range conditions):

1. `0 <= K_SLAM <= K_IS`, so `lambda_i(K_SLAM) <= lambda_i(K_IS)`.
2. `ker K_SLAM = {u : A u in Range(B)}` (includes `ker A`).
3. `rank K_IS - rank K_SLAM = dim(Range(A) cap Range(B))`.
4. `Range(A) cap Range(B) = {0}` iff `rank[A B] = rank A + rank B`; also `dim(U cap V) >= max(0, r_A + r_B - m)` with `m` the realified data dimension.
5. Directional retention identity `R(u) = sin^2 angle(Au, V)`.
6. No-prior generalized retention spectrum `rho_i = sin^2 theta_i`; at most `rank B` values differ from 1; ordinary eigenvectors can rotate broadly under a low-rank loss.
7. `rank L_X <= rank B`; Hermitian interlacing `lambda_{i+p}(K_IS) <= lambda_i(K_eff) <= lambda_i(K_IS)`.
8. Variational forms `u*K_SLAM u = min_z ||Au - Bz||^2`, `u*K_eff u = min_z(||Au-Bz||^2 + z*J_X z)`; monotonicity in `J_X`.
9. Exact `SE(2)/SE(3)` gauge forces `A xi_chi + B xi_X = 0` and `rho=0`; converse is false.
10. Multi-frequency: `u in ker K_SLAM^multi` iff there exists one shared `z` with `A_f u = B_f z` for every frequency f; consistent stacking never decreases absolute `K_eff`; duplicate blocks `(A_2,B_2)=c(A_1,B_1)` scale `K_IS` and `K_SLAM` by `1+|c|^2` and leave every `rho_i` unchanged only in the no-prior (or jointly scaled prior) case.
11. Born at nonzero `chi_0`: `B_l = (dA/dX_l) chi_0`; at `chi_0=0`: `B=0`, `K_SLAM=K_IS`, and leading geometry mismatch is the bilinear term `(D_X A[dX]) dchi`.
12. Volume retention on a task support `S` equals `prod sin^2(theta_i)`; `sum log rho_i` is a log-volume criterion.

Conditional statements carried with conditions (not upgraded to theorems by this phase):

- Prior limits require a fixed-rank regular path (e.g. `J_X = eps I`) and prior support covering every data-coupled pose direction; singular/unanchored graph gauges never converge to `K_IS`.
- `dot K_eff = dot A* W A + A* W dot A + A* dot W A`, `dot W = -dot B H B* - B H dot B* + B H dot C H B*`, `C = B*B+J_X` invertible (or consistent generalized Schur convention); `dot B` needs pose Hessian `D_X^2 F`; dropping it is the frozen-nuisance-subspace approximation.
- No-prior projector derivative `dot P_B = P_{B^perp} dot B B^+ + (P_{B^perp} dot B B^+)*` and `||dot K_SLAM|| <= 2||A||||dot A|| + 2||A||^2||B^+||||dot B||` require B locally constant rank.
- Robust surrogate `gamma_r(X,eps) = lambda_r(X) - eps ||grad_X lambda_r||_* + O(eps^2)` requires simple `lambda_r`, positive spectral gap, second-order smoothness, full norm-ball uncertainty set; a certified bound needs a uniform operator-Lipschitz constant over the whole ball.
- Full-wave bounds `||A_t|| <= ||G_S,t|| ||M_t^-1|| ||E_t^tot||_inf`; near `sigma_min(M_t)->0`, absolute information and sensitivity may both grow and the linearization radius shrinks.
- Far-field/Born asymptotics: `K_IS = O(R^-2)` in 2D monostatic; evanescent `k_perp` reach `k_perp,max ≈ sqrt(k^2 + (log SNR/(2z))^2)`; constants are convention-dependent.
- Infinite-dimensional versions require closed ranges, polar/frame formulations, and stable transversality `inf_{Au != 0} dist(Au,V)/||Au|| > 0`.

Open claims (never presented as results): generic transversality with full angular/frequency/MIMO diversity; single-frequency finite-aperture characterization and `theta_min` lower bounds from bandwidth; full-wave resonance information gain vs sensitivity; Fourier-surrogate error for `V_D^+`; continuity/rank-event theory and spectral class stability; exact two-fold SOM `V_D^±` definition; prior-weighted spectrum bounds; nonlinear/global basin questions (cycle skipping, phase ambiguity).

## 6. Per-source-file summaries

### 6.1 `Theory/SOM_SLAM_THEORY_CONTEXT.md` (3213 lines, Chinese)

Exact path: `/Volumes/migodam's-external-brain/Research/Inv_SLAM/Theory/SOM_SLAM_THEORY_CONTEXT.md`

Role: master consolidated theory ledger and anti-conflation reference, distilled from seven named ChatGPT "电磁逆散射" conversations (S1-S7 with recorded conversation IDs); all internal status labels and conflicts are preserved and are not peer-reviewed.

Technical summary (5-10 lines): It fixes the full-wave 2D scalar Helmholtz contrast-source model and notation, gives explicit Jacobians `A_t = G_S M^-1 diag(E^tot)` and the conditional `B`, then derives the information chain `K_IS`, `K_SLAM = A*P_{B^perp}A`, and Schur `K_eff`. It states and conditions the central identities: kernel/rank/intersection formulas, `rho_i = sin^2 theta_i`, low-rank loss `rank L_X <= rank B`, interlacing, variational/monotonicity forms, exact SE(2)/SE(3) gauge (`A xi_chi + B xi_X = 0`), and multi-frequency shared-compensation conditions. It separates three pose-uncertainty models (deterministic nuisance, random marginalization, bounded execution mismatch), realification rules, and treats trajectory derivatives, robust objectives, Born/far-field degeneracy, evanescent reach, MIMO/MSR/MUSIC/SOM boundaries, and the TSOM/tri-space hierarchy without merging `G_S` current modes with map-information modes. It includes a corrections/conflicts table, an open-questions register, and section 19's minimal 2D experiment recipe, plus a statement that literature mapping is discussion context only.

### 6.2 `Theory/Questions/Q1.md` (968 lines, Chinese)

Exact path: `/Volumes/migodam's-external-brain/Research/Inv_SLAM/Theory/Questions/Q1.md`

Role: original research question/prompt ("位姿不确定性下全波逆散射 SLAM 的谱可观测性极限"), asking for mathematical derivation rather than a deep-learning algorithm.

Technical summary (5-10 lines): It poses the central question of information loss in the observable subspace when pose is unknown, defines the Lippmann-Schwinger/contrast-source model and asks for explicit Jacobians (equations (1)-(3)), and defines `K_IS`, `K_SLAM`, `P_{B^perp}`, retention `R(dchi)`, and the generalized eigenproblem `K_SLAM v = rho K_IS v` together with `R_op`. It requests proofs for ordering, kernel characterization, rank loss, and degree-of-freedom loss in terms of `Range(A) cap Range(B)`, plus continuous spectra without hard thresholds, pose-prior Schur complements, trajectory derivatives and the robust problem `max_X min_{|dX|<=eps} lambda_r(K_eff)`. It lays out regimes (Born, full-wave small-pose, single/multi-frequency, straight/circular/partial-circle/multi-view/active trajectories, far-field asymptotics), asks for identifiability conditions, DoF-loss formulas, retention spectra, stability bounds, and trajectory-design principles, and sets the final deliverable as theorems/approximations with literature comparison, not an algorithm. It explicitly warns not to claim novelty without checking whether `A*P_{B^perp}A` and the retention spectrum already exist.

### 6.3 `Theory/Questions/A1.md` (3324 lines, Chinese)

Exact path: `/Volumes/migodam's-external-brain/Research/Inv_SLAM/Theory/Questions/A1.md`

Role: consolidated, rigorous answer ("总结判断") to the question thread, with proofs, equation-by-equation verdicts, conditions, numerical-experiment design, and literature placement.

Technical summary (5-10 lines): It confirms the core statement `K_eff = A* W_X A` with `W_X = P_{Ran(B)^perp}` (no prior) or `I - B(B*B+J_X)^+ B*` (prior), proves Theorem 1-5: monotonic ordering, `ker K_SLAM = {u:Au in Ran(B)}`, `rank K_IS - rank K_SLAM = dim(U cap V)`, `rho_i = sin^2 theta_i` via thin SVD and CS/principal angles, and monotonic multi-frequency stacking. It gives the variational forms, prior-limit conditions, prior-weighted `R_X = Q*(I+DD*)^-1 Q` and continuous `d_conf`, low-rank `L_X`/interlacing, explicit 2D `SE(2)` generators, the Born `chi_0=0` bilinear degeneration, far-field Fourier sampling and `K_IS=O(R^-2)`, full-wave operator bounds, full `DK_eff` formulas including `dot B = D_X^2 F`, the projector derivative bound, eigenvalue/subspace perturbation formulas, and the robust `gamma_r` surrogate plus conservative Lipschitz bound. It contains a table judging equations (1)-(18), section 19's minimal 2D experiment (32x32 or 40x40 grid, two-blob contrast, one Tx plus 4-8 body-fixed Rx, four equal-budget trajectories, prior sweep `alpha in [10^-6, 10^4]`, machine-precision rank tolerance), and identifies three candidate paper contributions and the next proofs.

### 6.4 `research/delegated/context_isolation/summary.md` (56 lines, English)

Exact path: `/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/context_isolation/summary.md`

Role: delegated structured summary of the theory context (problem, mechanism, assumptions, claims, uncertainties) used to isolate the context for the workshop seed.

Technical summary (5-10 lines): It restates the core problem as full-wave inverse-scattering SLAM with the central risk `Range(A) cap Range(B)`, and the converged question about low-rank generalized-eigenvalue deformation of the SOM/TSOM-recoverable subspace and rank-event tracking. It summarizes forward model, `K_IS`/`K_SLAM`/`K_eff`/`L_X`, the three distinct pose-uncertainty models, generalized eigenproblem `rho_i=sin^2 theta_i`, multi-frequency/trajectory consequences, and the algorithmic TSOM hierarchy. It lists assumptions (local linearization, whitening/realification, world-fixed map, invertible M, finite-dimensional/gauge discipline, prior-limit conditions) and accepted claims (kernel/rank/intersection, interlacing, low-rank loss, gauge one-way implication, Born empty-background bilinear term, multi-frequency monotonicity, robust surrogate conditions). It records uncertainties (two-fold SOM `V_D^±` definition, infinite-dimensional transfer, rank-event theory, prior-weighted spectra, Fourier-surrogate error, novelty unverified) and notes that the minimal 2D numerical verification is proposed but not yet run.

### 6.5 `research/literature/TARGETED_PRIOR_ART_ADDENDUM.md` (120 lines, English)

Exact path: `/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/literature/TARGETED_PRIOR_ART_ADDENDUM.md`

Role: targeted prior-art addendum (checked 2026-09-03) filling gaps left by a partial ScholarQA collection; defines mandatory neighboring work and the bounded novelty statement. Not exhaustive; cannot support a global novelty claim.

Technical summary (5-10 lines): It lists seven closest neighboring lines: (1) RF-SLAM posterior CRB for distributed MIMO SLAM (Deutschmann et al., Asilomar 2025, arXiv:2506.19957); (2) sparse blind deconvolution for distributed radar autofocus imaging with antenna-position ambiguity (Mansour et al., TCI 2018); (3) joint SAR image formation and phase-error estimation (Onhon/Cetin TIP 2012, Scarnati/Gelb JCP 2018); (4) unified bilinear identifiability up to transformation groups (Li/Lee/Bresler, arXiv:1501.06120); (5) mmWave radio-SLAM with geometric multipath maps (2024-2025 representative DOIs); (6) FIM/Jacobian observability for SLAM-based TDOA sensor calibration (Su et al., T-RO 2021); and (7) simultaneous inverse-scattering reconstruction and transmitter localization (Karthik/Ghosh, PIERS 2023), which is a mandatory direct-title neighbor so "joint map + transmitter localization" must not be claimed as new. Its bounded-novelty conclusion is that a potentially differentiating synthesis remains the whitened/realified map-tangent vs physical-pose-tangent geometry for a nonlinear full-wave volumetric contrast model with (i) pose-eliminated finite-rank defect information, (ii) no-prior principal-angle retention, (iii) Born-empty-background second-order operator-uncertainty separation, and (iv) frequency/trajectory/rank-event predictions in one reproducible solver; language must say "we formulate/test" and must not say "first"/"no prior work"/exhaustive.

## 7. Experiment details contained in the sources

Workshop "shared scenario" (also section 19 of the theory context and A1):

- Domain `D=[-0.5,0.5]^2`; grid `32x32` default with `40x40` sensitivity run.
- Nominal contrast `chi_0 = 0.3 chi_disk,1 + 0.5 chi_disk,2` (two smooth-edged blobs, deliberately nonzero). A smooth Fourier/spline/basis variant is used for infinitesimal gauge checks; piecewise-constant pixels only where tolerated.
- Sensors: one Tx plus 4-8 body-fixed Rx per pose; pose `x_t=(p_x,p_y,theta)`.
- Trajectories: straight line, 90-deg arc, 180-deg arc, 360-deg circle; total path length and number of measurement configurations held equal across trajectories in every comparison.
- Switches: Born vs full-wave; single vs multiple genuinely distinct frequencies `F in {1,2,3}`; duplicate scaled block `(A_2,B_2)=c(A_1,B_1)` as control; identity plus one colored covariance (`Sigma_n^-1` metric or whitening first); realification always after whitening.
- Nominal `M_t` must be invertible with reported `sigma_min(M_t)/||M_t||` margin; resonance cases labeled separately.
- Rank tolerance for exact algebraic checks: machine-precision type `tau = max(m,n) eps_mach sigma_1`; SOM physics thresholds never used to verify identities.
- Prior sweep `J_X = alpha I`, `alpha in [10^-6, 10^4]`, plus one singular-support `J_X`; ordering `K_SLAM <= K_eff(alpha) <= K_IS` checked with backward-error-scaled `delta`.

Five experiment families (all mandatory, no silent omission):

- Family 1: Forward model and analytic Jacobians -- derivative consistency and linearization validity (centered FD vs `A dchi`, `B dX`; h sweep; Born/full-wave contrast continuation; M margin).
- Family 2: Algebraic spine -- kernel/rank/retention identities, ordering, prior sweeps (`ker K_SLAM`, rank identity, `eig(R_op)=1-sigma_i^2(Q_B*Q_A)`, `rho_i in [0,1]`, at most rank-B deviations, `rank L_X <= rank B`, interlacing, prior limits).
- Family 3: Gauge and Born empty-background degeneration -- SE(2) generators lie in `ker K_SLAM` (smooth basis), Born `chi_0=0` gives `K_SLAM=K_IS` with bilinear mismatch dominance, anchors remove expected gauge directions.
- Family 4: Frequency stacking and trajectory geometry -- stacking never decreases absolute `K_eff`; duplicate-block invariance only in no-prior/scaled-prior settings; shared pose compensation `z`; trajectory ordering `circle > 180-deg arc > 90-deg arc > straight line` or documented counterexample.
- Family 5: Trajectory sensitivity, robust surrogate, and rank-event stability -- first-order `D lambda_i` predictions at `O(eps^2)`, `gamma_r` surrogate on full norm balls, smooth fixed-rank tracking, breakdown near gap closure/small `sigma_min(B)`, and quantified frozen-nuisance (`dot B`) error.

Every family must return claim status, exact executed checks, pass/fail, artifacts, and a "cannot establish" paragraph. No family result upgrades an [open] claim to [accepted]. Any further solver implementation belongs in a later phase.

## 8. Open questions register highlights (from workshop.md)

- (a) Fatal core-validity: local-tangent framing vs non-local SLAM failure; Born weak-scatterer regime fidelity; nuisance-model faithfulness as pose dimension grows; gauge quotient vs bounded-domain reality; discretization fidelity / stable transversality; what local retention predicts for system success.
- (b) Hard theory/proof work reserved for GPT Pro / deep-theory sessions: infinite-dimensional operator version; continuous Helmholtz SE(2)/SE(3) gauge theorem; Born far-field Fourier transversality; multi-frequency transversality; full-wave resolvent perturbation bounds; prior-weighted spectrum; rank-event/continuity theory; spectral-class/projector-composition theory.
- (c) Design/construction: exact two-fold SOM definition; Fourier/FFT surrogate; map parameterization and pullback `T_chi`; full model of B; robust trajectory formulation; continuous reporting conventions; incremental tracking design.
- (d) Ordinary pipeline tasks: unified 2D Helmholtz contrast-source solver with Born switch and analytic A/B builders; four-trajectory harness with frequency/prior/perturbation switches; spectral/algebra library; sweeps and comparison runs; reproducibility packaging. These are the tasks for later implementation, not this phase.

## 9. Uncertainties / caveats for downstream use

- workshop.md names five source files, not three; this summary covers all five and flags the discrepancy.
- All status labels and "claims" in the sources are internal-discussion labels, not peer-reviewed or externally verified; several are [conditional], [open], or [non-claim].
- Novelty and closest-prior-work boundaries are explicitly unverified; no "first"/"no prior work" statements are licensed.
- No solver, run, numeric check, figure, or result was produced in this phase; no acceptance gate is claimed.
