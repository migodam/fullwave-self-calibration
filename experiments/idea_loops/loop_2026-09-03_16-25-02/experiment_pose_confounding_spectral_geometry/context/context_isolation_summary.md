# SOM/SLAM Theory Context -- Structured Summary

Source: `Theory/SOM_SLAM_THEORY_CONTEXT.md` (consolidated 2026-09-03 local reference over seven prior inverse-scattering/wave-SLAM discussions). Status labels used by the source (derived/accepted, needs conditions, intuition/assumption, open, conflict/correction) are internal to those discussions and are not peer-reviewed or externally verified.

## Problem

- **Core problem**: full-wave inverse-scattering SLAM, i.e. jointly estimating platform trajectory and a volumetric medium contrast field from scattered-field measurements when pose is unknown. Inverse scattering alone assumes known sensor/background geometry; SLAM requires joint pose-map estimation under shared data association, priors, and gauge freedom.
- **Central risk**: local map--pose confusion. A contrast perturbation and a pose perturbation can produce the same first-order change in data space; the document characterizes this through `Range(A) ∩ Range(B)` for map Jacobian `A = D_χ F` and pose Jacobian `B = D_X F`.
- **Converged research question**: how finite-dimensional pose uncertainty deforms the SOM/TSOM-recoverable subspace, whether that deformation is a low-rank generalized-eigenvalue defect, and whether it can be tracked stably across trajectory-induced rank-changing events.
- **Scope boundaries**: distinguishes inverse scattering, imaging, Radio/multipath/Radar SLAM, and full-wave SLAM; also separates uniqueness/identifiability, local observability, stability, statistical precision (FIM/CRB), and global multimodality/cycle-skipping as five non-interchangeable notions.

## Core mechanism

- **Forward model**: scalar 2D Helmholtz with Lippmann--Schwinger equations: data equation `E^sca = G_S J` and state equation `J = D_χ(E^inc + G_D J)`, where `J` is contrast source. Eliminating `J` gives `F(χ) = G_S M^{-1} D_χ E^inc` with `M = I - D_χ G_D`; nonlinearity arises because the unknown medium changes the total field that defines the effective source.
- **First-order joint model**: `δy = A δχ + B δX + n` after noise whitening; `A = G_S M^{-1} diag(E^tot)` and `B` combines receive-geometry and illumination/domain-path terms (full form regime-dependent).
- **Information operators**: pose known gives `K_IS = A*A`; free deterministic pose elimination gives `K_SLAM = A* P_{B⊥} A`; pose prior `J_X` gives Schur-complement `K_eff = A*A - A*B (B*B+J_X)^† B*A`; pose-induced loss `L_X = K_IS - K_eff` has rank at most rank `B`.
- **Pose-uncertainty models are kept separate**: (A) deterministic nuisance pose jointly estimated and eliminated; (B) random pose error marginalized into effective noise (Woodbury-equivalent to `J_X = C_X^{-1}` under Gaussian assumptions); (C) bounded outer execution/model mismatch `‖ΔX‖ ≤ ε` analyzed through trajectory sensitivity of `K_eff`.
- **Map--pose geometry**: generalized eigenproblem `K_SLAM v = ρ K_IS v` on the observable support; no-prior retention satisfies `ρ_i = sin² θ_i` where `θ_i` are principal angles between `Range(A)` and `Range(B)` in data space. Absolute map strength and relative retention are separate axes.
- **Trajectory consequences**: multi-frequency stacking removes confusion only when one shared pose compensation fits all frequencies; adding truly independent blocks monotonically increases absolute `K_eff`; continuous frame/normal operators express sampling-limited vs physics-limited information; evanescent high-frequency reach decays exponentially with distance and improves only logarithmically with SNR.
- **Algorithmic structure**: SOM/TSOM compresses current ambiguity using `G_S` and internal domain structure first; pose defect is a third layer in retained-current/reduced coordinates. FFT-TSOM is only a subspace-level approximation intuition; ML is deliberately excluded from the core theory until objectives, equivariance, calibration, and data-consistency guarantees exist.

## Assumptions

- Local linearization about a nominal solution; whitened measurements and realification of complex data before computing projectors, principal angles, and Schur complements when physical pose/map parameters are real.
- `M` invertible / stable solution branch away from internal resonance or bifurcation in the full-wave `A` expression.
- World-fixed discretization (so `G_D`, `D_χ` do not move with platform) or an explicit alternative parameterization; no undocumented calibration/direct-path terms in simplified `B`.
- Finite-dimensional Hermitian setting with closed ranges and gauge fixed (or analysis restricted to identifiable support / Moore--Penrose with gauge variance not read as ordinary covariance) for rank, kernel, interlacing, and CRB-style statements.
- Deterministic-nuisance vs Gaussian-random vs bounded-mismatch pose models are used in the correct regime and not double-counted; prior limits (`J_X → 0` or `→∞`) need a fixed-rank regular path with the right support.
- Principal-angle statements require the no-prior projector case; a finite pose prior produces weighted shrinkage, not ordinary `sin² θ`.
- Hard spectral classification (gauge/mixed/robust; `V_S^±`, `V_D^±`) requires threshold/gap/stability reporting; an exact `2×2×2` decomposition of three noncommuting projectors is generally invalid.
- Born, far-field, and evanescent-reach formulas are regime-dependent approximations with convention-dependent constants.
- Regularization improves stability/bias--variance but does not create measurement Fisher information.
- Literature mapping reflects the original discussions only; it is not a fresh systematic search, and novelty is not established.

## Claims

- **Derived/accepted within the discussions** (finite-dimensional, whitened, realified, locally linearized, gauge-appropriate):
  - `ker K_SLAM = {u : Au ∈ Range(B)}`; `rank K_IS - rank K_SLAM = dim(Range(A) ∩ Range(B))`; `rank L_X ≤ rank B`; Hermitian interlacing `λ_{i+p}(K_IS) ≤ λ_i(K_eff) ≤ λ_i(K_IS)` for `rank L_X ≤ p`.
  - Generalized retention eigenvalues satisfy `ρ_i = sin² θ_i` with no pose prior; at most `rank B` such values differ from 1 on the observable support, yet ordinary eigenvectors can rotate broadly even under a low-rank loss.
  - Exact global `SE(2)/SE(3)` gauge implies `ρ = 0`; the converse is false (local map--pose tangent coincidence need not come from a group symmetry).
  - In the Born model at empty background `χ_0 = 0`, first-order `B = 0` does not make pose error harmless; confusion enters through the second-order bilinear term `(D_X A[δX]) δχ`.
  - Consistent multi-frequency stacking is monotonically nondecreasing in absolute `K_eff` information; duplicated blocks improve SNR but not relative retention.
  - Robust trajectory objective admits a first-order `λ_r - ε ‖∇λ_r‖_*` surrogate under simplicity, spectral-gap, smoothness, and full-norm-ball assumptions; Lipschitz bounds give more conservative crossing-safe versions.
- **Deliberate non-claims**: SOM mode counts are about stable data-resolvable current directions, not intrinsic scene dimensionality; a 3D volumetric output need not be independently data-supported; "theoretically unique", "large local Fisher information", "algorithm converges", and "output looks reasonable" are separate propositions; `G_S` current spectrum, MSR/MUSIC data subspaces, map Jacobian spectrum, and pose range are different objects.

## Uncertainties

- Two-fold SOM's exact definition of `V_D^±` (which operator, projection order, optimization variables, and Chen-version conventions) is unresolved and must be checked against the specific cited paper.
- Infinite-dimensional versions need extra conditions (closed ranges, polar/frame formulations, stable transversality with positive infimum); finite-dimensional dimension formulas do not transfer automatically.
- Continuity/tracking theory requires fixed rank and positive external spectral gap; rank events, eigenvalue crossings, and Grassmann/Flag bookkeeping lack a complete proof.
- Spectral classification thresholds, noncommuting-projector approximations, and prior-weighted shrinkage behavior need concrete gap/stability analysis and error bounds.
- Full-wave `B`, resolvent sensitivity `‖M^{-1}‖` amplification, Fourier surrogate error for `V_D^+`, prior-weighted spectrum bounds, and multi-frequency transversality conditions are stated as open targets, not proven theorems.
- Regime-dependent claims (far-field amplitude, evanescent reach with `k_{⊥,max}`, resolution/reach limits, physics-vs-sampling split) carry unspecified constants and thresholds; exact constants depend on measurement conventions.
- Robust first-order objectives require smoothness, single eigenvalues with gap, and the correct uncertainty set/dual norm; validity at rank events or under cone-constrained perturbations is open.
- The minimal 2D numerical verification (finite-difference Jacobian checks, rank identities, gauge generators, prior sweeps, trajectory-perturbation checks, Born/full-wave and single/multi-frequency comparisons) is proposed but not yet run.
- Novelty and closest-prior-work boundaries are explicitly unverified; recent literature and "first" claims require fresh systematic re-check before submission.
