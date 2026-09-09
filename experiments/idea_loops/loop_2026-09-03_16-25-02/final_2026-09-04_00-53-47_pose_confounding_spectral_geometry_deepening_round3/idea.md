## Name

pose_confounding_spectral_geometry_deepening_round3

## Title

Pose-Confounding Spectral Geometry for Full-Wave Inverse-Scattering SLAM: Replication, Mechanism, and Robustness Demotion

## Short Hypothesis

After noise whitening and realification, the finite-dimensional pose-eliminated map-information operator is governed by principal angles between map- and pose-induced data tangent spaces. This round tests the currently weak design/robustness claims: trajectory ordering is replicated or mechanistically qualified across scenes and frequency sets; multi-frequency diversity is characterized beyond two frequencies; the machine-rank boundary event is resolved with stable-rank diagnostics; and execution-error robustness is claimed only for the affine tangent layer unless a full nonlinear certification is obtained.

## Related Work

The individual ingredients are mature and must not be claimed as novel: Chen-style subspace optimization, Schur-complement nuisance Fisher information, principal angles and finite-rank Hermitian perturbation, gauge analysis, active sensing, RF-SLAM PCRB/identifiability, blind bilinear inverse problems, and Born/Rytov scattering approximations. Prior completed rounds verified the whitened/realified algebraic spine, the principal-angle retention identity under colored noise and grid refinement, the G_S/SOM versus A/K_eff map-tangent distinction, a preserved trajectory-order counterexample, non-monotone two-frequency SNR diversity, and an affine tangent Lipschitz certificate. The remaining contribution is a finite-dimensional 2D scalar Helmholtz synthesis with replicated negative/design results, stable-rank diagnostics, and explicit robustness demotion where formal certification is unavailable.

## Abstract

We continue a finite-dimensional study of pose-confounded full-wave inverse-scattering SLAM. Prior rounds verified the whitened/realified algebraic spine, principal-angle retention identities under colored noise and grid refinement, and the distinction between map-tangent information (A,K_eff) and current-space G_S/SOM modes. Remaining weaknesses are the single-setting trajectory reversal, narrow frequency-diversity evidence, a non-certified nonlinear robustness bound, and one rank-boundary replication event. This revision therefore (i) replicates and dissects trajectory orderings across four scene classes, multiple frequency sets, and trajectory parameter variations; (ii) tests frequency-set and bandwidth effects with equal-budget multi-frequency designs; (iii) investigates every rank identity with stable-rank/tolerance sweeps and singular-value diagnostics; (iv) either certifies full nonlinear execution-error robustness or demotes it to an explicitly empirical/structural result while keeping the affine tangent certificate; and (v) adds an external literature-derived Born/weak-scattering information comparison. All claims remain scoped to the implemented 2D scalar Helmholtz model and finite discretizations.

## Experiments

### datasets
['Deterministic synthetic 2D scalar Helmholtz contrast-source scenes on D=[-0.5,0.5]^2: two smooth blobs reference, ring scene, offset/edge scene, low-contrast scene, and if CPU permits a rectangular inclusion scene. N=24 pilot sweeps plus N=32 reference and N=40 sensitivity; keep the same solver.', 'Body-fixed 4-8 receiver arrays on straight, arc30, arc60, arc90, arc120, arc180, and circle360 trajectories with matched path length, measurement count, standoff, and per-pose energy where applicable.', 'Frequency sets: single f0, two-frequency f0/f1 with SNR sweep, three- and five-frequency equal-count/equal-energy sets over bandwidths 0.2-1.0, and a duplicate-block control; multiple scenes.', 'Rank-event and robustness scenes: the preserving ring scene and perturbed ring/offset scenes, plus fixed-seed variants used for the machine-rank boundary event.']

### baselines
['Known-pose K_IS vs free-pose K_SLAM vs finite-prior K_eff; analytic A/B vs finite differences; principal-angle closed forms vs generalized eigenvalues.', 'Raw trajectory ordering vs same-standoff/per-pose-energy normalized ordering, now across scenes, frequencies, and trajectory parameter sweeps.', 'Stable-rank/truncated-SVD threshold sweep vs fixed machine-rank tolerance for all rank/kernel identities.', 'Born/Rytov weak-scattering analytic Fisher bound or literature-derived diagonal-Dominant Born bound compared with computed K_eff and K_IS retention eigenvalues.', 'Affine tangent Lipschitz certificate vs sampled full nonlinear perturbations; any full nonlinear bound explicitly labeled certified or empirical/structural.', 'Frequency-set designs with identical total measurement count and energy: genuine diversity vs duplicate blocks; monotonicity and plateau checks.']

### metrics
['Retention eigenvalues, principal angles, exact destroyed DOF, rank(A), rank(B), rank(AB), stable rank curves, singular-value gaps, and condition numbers.', 'Trajectory ordering replicability: fraction of scene/frequency/parameter variants preserving a common ordering; decomposition into B-rank, pose-compensation residual, illumination diversity, and curvature-related pose-Hessian terms.', 'Frequency diversity: retained DOF and log-volume retention for 2-, 3-, and 5-frequency sets; per-frequency SNR/bandwidth plateau; duplicate-block invariance controls.', 'Rank event autopsy: numerical rank vs tolerance, truncation error, singular values, backward error of matrix identities, scene conditioning, and perturbed-scene boundary location.', 'Robustness: affine tangent certificate violations over 2000 directions and epsilon sweep; full nonlinear sampled violations; if formal certification is obtained, its constants and proof conditions; otherwise explicit demotion statement.', 'External-comparison discrepancy: relative difference between computed K_eff/K_IS retention spectra and Born analytic predictions in weak-scattering scene; cite source or label as derived in appendix.', 'Exact commands, configs, seeds, runtime/environment, result tables, figure paths, and SHA-256 hashes.']

### compute_estimate
Apple Silicon CPU only. Reuse Green matrices and factorizations from previous rounds. Run N=24 pilot sweeps for trajectory/frequency/scenes, then selected N=32 reference and N=40 rank/robustness cases. No CUDA, GPU, neural nets, daemons, schedulers, or new orchestration. Document actual runtime and any incomplete family.

### steps
['Before coding, read workshop.md and current synthesis notes; preserve status ledger and semantic drift prohibitions.', 'Family 12: trajectory replication and mechanism, following context/PARENT_ROUND3_AUDIT.md. Use new family12 files; do not overwrite the existing family10 covariance toy. Run straight, arc90, arc180, circle360 on four scene classes and at least two frequency sets under predeclared equal-budget controls; report conditional mechanism, not a universal ordering.', 'Family 11b (optional): deepen the already completed Family 11 frequency study only if an equal-total-measurement and equal-total-whitened-energy design answers a new mechanistic question; never overwrite existing family11 artifacts.', 'Family 13: rank event autopsy. Reproduce the ring-scene machine-rank boundary with predeclared stable-rank/tolerance sweeps and singular-value gaps; determine whether the 6/7 count was a true rank change, tolerance classification, or unstable subtraction.', 'Family 14: robustness scope. Integrate the parent generalized-eigenvalue derivative and cluster/rank-event checks; unless a genuine interval proof is obtained, explicitly demote all full nonlinear execution-error robustness claims to empirical/structural observations.', 'Family 15: Born analytic/literature control. Derive and execute an internal weak-scattering formula check, label it analytic rather than external unless a published setup is quantitatively matched, and never invent or overstate citations.', 'Update synthesis and manuscript: replicated negative/mechanistic trajectory section, multi-frequency section, rank event section, robustness demotion section, and computed external comparison; keep finite-dimensional scope and explicit limitations.']

## Risk Factors And Limitations

- Any trajectory ordering result is likely model/scene/frequency specific; the goal is a replicated counterexample/mechanism, not a universal trajectory design rule.
- Frequency-diversity conclusions remain finite-dimensional and local; no claim of optimal frequency transfer to continuum or real RF systems is made.
- Full nonlinear Helmholtz operator-Lipschitz certification may remain open; all non-certified robustness claims will be demoted to empirical/structural observations.
- The ring-scene rank event may be a numerical tolerance artifact; if so, it must be reported as such, not hidden or resolved into a false stable-rank theorem.
- External comparison is limited to a Born/weak-scattering analytic bound; it is not an independent implementation of the full inverse-scattering model or a real-world dataset.
- Second model/geometry is not added by default; if added it will be a diagnostic, not a claim of generality.
- Finite-dimensional identities do not prove continuum transversality, global nonlinear convergence, or online SLAM performance.
- CPU-only budget may force selected sweeps rather than exhaustive combinations; incomplete family results will be documented explicitly.

