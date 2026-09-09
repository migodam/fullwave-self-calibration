# Revised manuscript blueprint after the autonomous research loop

Date: 2026-09-04

Working title: **Geometry-Lifted SOM under Unknown Array Geometry: Canonical
Equivalent-Current Coordinates, Pose-Hiding Bounds, and Limits of
Self-Calibration**

“TriSpace SOM” remains the project/framework name. In the manuscript it means
the typed three-channel pose signature `(current lift, data leakage, state
defect)`; it is not presented as a third independent spectral space.

## Central paper sentence

Unknown array geometry can be represented non-vacuously in SOM coordinates
only after a retained current space and metric are declared: receiver-side
sampling perturbations then admit a canonical minimum-norm retained-current
lift plus a leakage residual, while a nuisance-space dimension bound states
when pose is necessarily hidden. The tested state-consistency mechanisms do
not recover those hidden directions, so self-calibration requires additional
acquisition diversity, anchors, or priors.

## Supported contribution list

1. A model-conditional split of receiver re-sampling and
   transmitter/illumination-induced physical-current variation on a
   world-fixed grid, including the double-counting failure mode.
2. A canonical retained-SOM lift and orthogonal retained-basis leakage, with
   the unrestricted full-row-rank lift exposed as a vacuous negative control.
3. A finite-dimensional nuisance-space dimension obstruction for local pose
   visibility after retained current and map variations are eliminated.
4. Reproducible positive and negative evidence: spectral-coordinate drift and
   source-diversity rank recovery are observed; state-witness discrimination
   and hidden-direction recovery fail in the tested regimes.

The Moore--Penrose algebra, SOM/TSOM, equivalent currents, joint inversion,
blind calibration, source/receiver extension, and use of data/state residuals
are antecedents rather than contributions.

## Paper structure

### 1. Introduction

- Unknown geometry makes the operator defining SOM coordinates uncertain.
- Bare equivalent-current existence is vacuous at full row rank.
- State the four bounded contributions and headline the negative result.
- Scope: deterministic synthetic 2D scalar Helmholtz, local finite-dimensional
  claims; no global novelty, convergence, hardware, or 3D Maxwell claim.

### 2. Related work and originality firewall

- Classical contrast source, SOM, Twofold SOM, CSI/CC-CSI, virtual experiments.
- Microwave calibration and Bellomo's incident-field/Green-operator/phase-center
  correction.
- Joint inverse scattering plus transmitter localization, blind calibration,
  radar/SAR autofocus.
- Source/receiver-extension FWI.
- End with an antecedent-versus-increment table.

### 3. Typed full-wave model

- Whitened/realified spaces and real physical parameters.
- Independent-current data and state equations.
- `H_S` receiver pseudo-current target versus `M^{-1}H_D` physical current.
- World-fixed versus body-fixed grids.
- State explicitly that numerical `chi` is scattering potential under the
  present Green normalization; fixed coefficients across frequency are not a
  fixed-permittivity multifrequency material experiment.

### 4. Canonical retained-SOM lift

- Define `K_U,C_U,Q_U,R_U`.
- Prove least-squares/minimum-norm/orthogonality and exactness criterion.
- Give the full-row-rank vacuity corollary and `1/sigma_r` stability bound.
- Treat `Range(Q_U)` as a subordinate, scene-dependent tangent.

### 5. Pose hiding after nuisance elimination

- Define the real nuisance range `N_U` and `B_vis=P_Nperp B`.
- Prove `hidden_dim >= max(0,p+dim N_U-m)`.
- Specialize to `m=2M,p=3,K=3` and derive the `M-2/M-1` regimes.
- Give the exact state solvability condition `R_Uh=T_Uh=0`, immediately
  disclosing that the state half never became binding in the experiments.
- Distinguish this data-side branch from the earlier map--pose tangent overlap
  branch based on `A`, `B`, and `Q_A^*Q_B`; merge only in the joint Jacobian.

### 6. Implementable rank-aware update

- Recompute the pose-dependent operator and retained basis.
- Track projectors on fixed-rank strata; use a soft filter near rank events.
- Compute the canonical lift and leakage for diagnosis.
- Solve direct joint GN/LM or its exact Schur/Woodbury pose correction on the
  visible quotient; do not freeze hidden directions and call that recovery.
- If the stacked projected pose Jacobian is deficient or ill-conditioned,
  request anchors/priors or alter source, aperture, frequency, or trajectory.
- Label the executed code SOM-truncated, not TSOM.

### 7. Experiments

- E1: unrestricted vacuity and retained lift/leakage.
- E2: pose-dependent basis drift and hard/soft rank events.
- E3: receiver pseudo-current versus transmitter physical-current derivative.
- E4: data-side hiding and the non-discriminating state witness.
- E5: direct joint versus visible-coordinate restriction, including corrected
  saturation ranks, residuals, seeds, and half-aperture local minima.
- Extensions: multi-transmitter/directional-source rank, robustness, and hard/
  soft state-consistency nulls.
- Never use the old E4 “joint data+state” figure or old E10 rank curve as
  evidence; use the parent-corrected saturation figure.

### 8. Discussion and limitations

- What is physically meaningful versus textbook algebra.
- Why optimization success flags are not recovery certificates.
- Why source diversity removes a source-symmetry deficiency but not the global
  map--pose gauge.
- Open continuum, TSOM, nonlinear basin, Maxwell, coupling, clock, and real-data
  problems.

### 9. Conclusion

- The contribution is a non-vacuous coordinate/diagnostic theory and a hiding
  bound, not a completed universal self-calibration solver.
- Additional independent information is necessary to recover hidden pose.

## Figure policy

Use:

1. `validation/rank_saturation_corrected.png` for the corrected dimension
   obstruction and nonlinear control;
2. `plot_fd_errors.png` for derivative validation;
3. `plot_e6b_multitx_robust.png` for source-diversity/SNR behavior;
4. `plot_e11_vp_state_null.png` for the negative state-consistency result.

Do not use the locked `plot_e10_settings_sweep.png` rank interpretation or
`plot_e4_hiding.png` state label without an explicit invalidation overlay.

## Acceptance status

- Theory seed: independently sanity-checked.
- Five planned experiments: attempted; expanded to E1--E11.
- Autonomous loop: three development rounds, final decision locked.
- Parent rank correction: executed independently and retained separately.
- Scientific status: **major revise / paper draft**, not publication-ready.
- Next decisive gates: exact TSOM primary-source implementation, stacked
  acquisition theorem/design, full-text novelty audit, and higher-fidelity or
  measured-data validation.
