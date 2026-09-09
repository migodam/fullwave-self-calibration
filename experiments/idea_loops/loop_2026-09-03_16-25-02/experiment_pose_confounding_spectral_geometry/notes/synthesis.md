# Pose-confounding spectral geometry: consolidated synthesis (read-only)

Date: 2026-09-03. Working directory:
`/Volumes/migodam's-external-brain/Research/Inv_SLAM/experiments/idea_loops/loop_2026-09-03_16-25-02/experiment_pose_confounding_spectral_geometry`.

This synthesis is a read-only consolidation of the completed experiment
records. No experiment was rerun. Every number below was copied from the JSON
result files and notes named in the section and table rows;
`notes/manuscript_sources.md` records the exact source of each manuscript
table/claim. Executed statuses and preserved failures are reported exactly as
recorded; no recorded failure has been relabelled as a pass.

Binding audit context: `context/PARENT_ROUND_AUDIT.md` and
`context/PARENT_CORRECTIONS.md` (both 2026-09-03), plus per-family parent
reviews `PARENT_FAMILY2_REVIEW.md`, `PARENT_FAMILY3_CORRECTIONS.md`, and
`PARENT_FAMILY4_AUDIT.md`.

## 1. Claim-status tables (as recorded)

Status values use the recorded vocabulary: PASS / FAIL / NOT MET / flagged /
observed / recorded. Where a recorded check is a "PASS with numerical caveat"
or "PASS (as limiting behaviour)", the caveat is preserved.

### 1.1 Family 1 - forward model, corrected self-cell, analytic Jacobians

| claim | recorded status | executed check | key numbers (from JSONs) |
|---|---|---|---|
| corrected equal-area self-cell integral v2 matches direct quadrature | PASS (finite ladder) | 2048-node Gauss-Legendre radial quadrature, `results/self_cell_quadrature_validation.json` | max abs quadrature error 3.13e-16; abs I_self 5.12e-3 (N=8) -> 4.36e-5 (N=128), first/last ratio 117.6; invalidated v1 stays ~2.5e-2 (its 1/k_b^2 limit) |
| `grad_s g = -grad_z g` preserved | PASS | source-coordinate finite-difference unit test | corrected rel err 1.737e-9 vs old-sign 2.000 (discriminating test) |
| map-Jacobian FD consistency, O(h^2) region | PASS | centred FD sweep h=1e-4..1e-1, 3 map seeds | rel errors at h=1e-3: 2.25e-9 / 1.74e-8 / 1.77e-9; fitted slopes 1.976 / 1.968 / 2.001 (r^2 >= 0.9997) |
| pose-Jacobian FD consistency, O(h^2) region | PASS | centred FD sweep, 3 pose seeds | rel errors at h=1e-3: 8.38e-7 / 4.28e-6 / 4.09e-6; fitted slope 2.0000 (seed 10), ~1.999995 (seeds 11, 12) |
| non-resonant state margin and residual | PASS | state-equation residual and M spectrum | sigma_min(M)/||M||_2 = 0.7189, sigma_min(M) = 0.8329, max state residual 9.34e-16, max abs E_tot 8.67e-2 |
| Born/full-wave discrepancy shrinks towards weak scattering | PASS (recorded trend) | contrast sweep s = 0.01..0.2 | rel F discrepancy 1.77e-3 .. 3.54e-2; rel A Frobenius 1.89e-3 .. 3.82e-2; ||D_chi G_D||_2 = 0.326 at s=1 |
| forward-model grid refinement N in {16,24,32,40} | PASS (diagnostic; both planned grids completed) | rel differences vs finest N=40 | 1.45e-3 (N16), 5.03e-4 (N24), 1.60e-4 (N32); not continuum convergence |

### 1.2 Family 2 - algebraic spine (kernel/rank/retention/prior/duplicates)

| claim | recorded status | executed check | key numbers (pixel / smooth unless noted) |
|---|---|---|---|
| `ker K_SLAM = {u : A u in Range(B)}` | PASS (smooth); FAIL on eigenprojector gate, identity certified by residuals (pixel) | c1 subspace distance + null residuals | smooth dist = 0; pixel dist = 8.47e-7 vs 1e-8 gate (`pass_gate=false`, `identity_support=true`): nullity 226 = 226, min singular(N1^T N2) = 1 - 4e-13, max ||K_SLAM N2|| = 1.96e-19, max ||C N1|| = 1.56e-13, rank(K_SLAM) = rank(C) = 30 |
| rank identity `r_KIS - r_KS = r_A + r_B - r_AB` | PASS | c2 exact integer residual | pixel 48-30 = 18 = 48+18-48, residual 0; smooth 24-24 = 0 = 24+18-42, residual 0 |
| no-prior retention = 1 - sigma_i^2(Q_B^T Q_A), padded | PASS | c3 predicted spectrum | max |rho - rho_pred| = 2.89e-15 / 1.55e-15 (gate 1e-10); count < 1 equals r_B = 18 in both |
| at most r_B retention values differ from 1 | PASS | c3 count | 18 = 18 in both bases |
| `rank L_X <= rank B`, alpha in 1e-6..1e2 | PASS | c4 stable factor rank | max factor rank 18 = 18 both; direct-difference rank inflated at alpha=100 (pixel 129, smooth 23), documented as cancellation floor and not used for the gate |
| `K_SLAM <= K_eff(alpha) <= K_IS` and interlacing | PASS | c5 PSD and interlacing residuals | min eigenvalues >= -2.5e-19 / -2.9e-20 (gate -1e-10); max positive interlacing violation 5.0e-21 / 6.7e-21 |
| regular prior alpha -> 0 path approaches K_SLAM | PASS (as limiting behaviour, not at grid left edge) | c6 + parent SVD-filter sweep | d0(alpha=1e-6) = 0.209 / 0.137 (not asymptotic); monotone to d0(1e-14) = 5.70e-9 / 3.67e-9 |
| regular prior alpha -> inf approaches K_IS | PASS | c6 | dinf(1e4) = 5.17e-7 / 5.50e-7 |
| finite prior gives weighted shrinkage in (0,1] | PASS | c6 | rho_X in (0.0862, 1] pixel, (0.0868, 1] smooth at alpha=1e-3; max sorted-spectrum change vs no prior 0.9995 / 0.9889 |
| singular-support prior does not converge to K_IS as alpha -> inf | PASS (non-convergence); illustrative >0.1 magnitude not met | c7 | dinf_sing(1e4) = 4.71e-4 / 2.30e-4 (417-912x regular-prior dinf at same alpha); `example_magnitude_gate_gt_0_1 = false`; direction is an unpenalized nuisance coordinate, not called a gauge |
| stable `C^T C` kernel factorization (parent correction) | PASS | factorized residual | 3.77e-16 (pixel) / 2.70e-16 (smooth), gate 1e-12; does not retro-pass the original eigenprojector gate |
| duplicate no prior: K_IS, K_SLAM scale by 1+c^2, rho unchanged | PASS | c8(a) | rel K_IS 1.11e-15 / 1.08e-15; rel K_SLAM 2.53e-13 / 2.90e-13; max |rho2-rho| 6.66e-16 / 1.19e-12 |
| duplicate with fixed finite prior changes rho | PASS (expected nonzero change) | c8(b) | max |rho_X2 - rho_X| = 3.99e-2 / 3.98e-2 at alpha=1; monotone decrease over alpha in 1e-3..1e2 |
| duplicate with jointly scaled prior reproduces single system | PASS | c8(c) | rel K_eff 1.00e-15 / 1.04e-15; max rho diff 1.67e-15 / 1.22e-15 |
| smooth c1-c3 persist at N=32, N=40 | PASS | resolution sweep | c1 distances 0; c2 residuals 0; c3 residuals 1.90e-15 / 8.24e-15; r_B = 18 |

### 1.3 Family 3 - gauge, Born empty background, anchors (raw record, unchanged)

| claim | recorded status | executed check | key numbers |
|---|---|---|---|
| SE(2) generators cancel in the p=24 smooth basis (residual <= 5e-2, decreasing) | FAIL | gauge residual N=16 and N in {16,32,40} | Tx 1.829e-1 -> 1.828e-1; Ty 8.648e-2 (N16); Rot 3.102e-1 -> 3.093e-1; representation residuals 0.781 / 0.572 / 0.703 |
| smooth generator lies in ker K_SLAM | NOT MET at gauge level | kernel distance ||K_SLAM c_g||/||c_g|| | 1.91e-8 (Tx), 2.31e-9 (Ty), 1.43e-8 (Rot) - nonzero, representation-limited |
| pixel basis (documented larger-residual control) | Expected direction NOT reproduced | pixel gauge residual (N16) | Tx 6.61e-5, Ty 1.58e-5, Rot 4.68e-5 - 1000-20000x smaller than smooth; not an exact-equivariance certificate |
| symmetric scene rotation is not a map-gauge direction | PASS (degenerate / pose-only stabilizer) | dchi_rot, gauge residual | ||dchi_rot|| = 7.8e-17 (numerically zero); absolute ||B dX_rot|| = 1.88e-6 (1.28e-5 of translation scale); normalized residual ~1.0 (ill-scaled), interpreted as pose-only stabilizer |
| smooth retention spectra recorded | PASS (recorded) | R_op spectra | two-blob rho_min = 5.51e-11; symmetric rho_min = 2.20e-11 |
| anchors/known background/fixed boundary remove translation-x gauge (rho_min > 10x base) | FAIL | restricted smooth spectra | ratios vs base 5.51e-11: corner anchor 1.26x, known-background disk 0.011x, outer-ring boundary 0.244x; disk/boundary decrease rho_min |
| B = 0 at empty background (chi0 = 0) | PASS | Frobenius ratio gate < 1e-12 | 0.0 |
| K_SLAM0 = K_IS0 at empty background | PASS | relative Frobenius gate < 1e-12 | 0.0 |
| bilinear term dominates linear Born term | PASS (metric) | ||BIL|| / ||A0_R dchi|| at h=1e-3 | 2.81 |
| full-wave FD matches bilinear term at rel_err < 1e-3 with O(h^2) slope | FAIL | FD of full-wave map at unit-L2 dchi | rel_err(1e-3) = 2.18e-1, h-independent (clean slope 8.6e-8); chi^2-order pose derivatives dominate |
| second-order remainder R(h)-R0 has O(h^2) slope | NOT MET | remainder vs h | slope 1.015 (linear in h); ||R(1e-3)-R0|| = 7.14e-6; ||R(1e-3)|| = R0 = 2.51e-3 |

### 1.4 Family 3b - correctly-conditioned refinement (supplement; raw 3 failures preserved)

| claim | recorded status | executed check | key numbers |
|---|---|---|---|
| A: pixel gauge residual <= 1e-3 and non-increasing N16->N40 | FAIL on non-increasing condition; PASS on 1e-3 magnitude | Tx/Ty/Rot residuals | Tx 6.61e-5 -> 7.95e-5 -> 8.14e-5; Ty 1.58e-5 -> 2.11e-5 -> 2.17e-5; Rot 4.68e-5 -> 6.07e-5 -> 6.24e-5 (all <= 1e-3, none non-increasing) |
| A: pixel kernel distance tiny/decreasing | PASS (recorded) | ||K_SLAM dchi_g||/||dchi_g|| | Tx 3.96e-10 -> 1.22e-10 -> 8.05e-11; Ty 8.73e-11 -> 3.11e-11 -> 2.04e-11; Rot 1.85e-10 -> 6.47e-11 -> 4.27e-11 |
| B: p=24 smooth basis reproduces Family 3 record | PASS (consistency) | side-by-side recomputation | residuals 0.182892 / 0.086484 / 0.310163; max gauge-residual diff vs raw JSON = 0.0 |
| B: p=27 augmented smooth basis represents generators | PASS | representation residual < 1e-10 | 1.42e-15 (Tx), 2.04e-15 (Ty), 1.63e-15 (Rot) |
| B: augmented smooth gauge reaches pixel forward floor | PASS | gauge residual <= 1e-3 | Tx 6.61e-5, Ty 1.58e-5, Rot 4.68e-5 (equal to pixel residuals) |
| C: unanchored translation-x generator retention ~ 0 | PASS | generator-specific rho_g | base rho_g = 9.908e-11 <= 1e-8 |
| C: corner anchor removes translation-x generator (ratio > 10) | PASS | same rho_g under mask | anchor rho_g = 1.182e-9, ratio 11.93 |
| C: known-disk and outer-ring masks (reported, not forced) | PASS on ratio, reported honestly | same generator metric | known-disk rho_g = 2.165e-5 (ratio 2.19e5); outer-ring rho_g = 1.645e-8 (ratio 1.66e2); residual after restriction 0.105 / 7.8e-4 |
| D: A(chi=0) equals Born map | PASS | relative Frobenius | global 7.05e-17; max per-pose 7.09e-17 (gate 1e-12) |
| D: Born FD O(h^2) slope gate | FAIL / not established (degenerate comparison) | rel_err_born vs h | rel_err_born(1e-3) = 2.27e-13 < 1e-4 passes value gate; fitted slope -1.025 because FD_born and BIL are centred differences of the same verified-identical linear map |
| D: full-wave error has small h-independent chi^2 floor | PASS | rel_err_full vs h | 2.166e-4 at every h in 1e-3..1e-1; fitted slope 5.4e-5 (floor dominated) |
| D: second-order remainder drift slope ~2 | PASS | ||R_full(h) - R_full(0)|| | slope 1.986 (r^2 = 0.99999); at h=1e-2 = 1.025e-8; R0 norm = 2.489e-9 |
| D: bilinear-to-linear dominance at small amplitude | PASS (metric) | ||BIL||/||A0_R dchi_small|| | 2.8087 at h=1e-2 (||dchi_small||_2 = 1e-3) |

### 1.5 Family 4 - frequency stacking and trajectory geometry

| claim | recorded status | executed check | key numbers |
|---|---|---|---|
| A: PSD monotonicity of K_SLAM and K_eff under distinct-frequency stacking | PASS | min eig of prefix increments (raw identity-noise metric) | F1->2: K_SLAM 1.517e-13, K_eff 1.288e-14; F2->3: K_SLAM 1.149e-13, K_eff 5.350e-14 (tol 1e-10) |
| B: duplicate, no prior: K scales and rho invariant | PASS | rel Frobenius + rho | rel K_IS = 1.080e-15; rel K_SLAM = 2.893e-13; max |rho_dup - rho_1| = 1.185e-12 |
| B: fixed finite prior breaks duplicate invariance | observed (expected) | max rho / rel Frobenius | max |rho_X dup - single| = 3.984e-2; rel F K_eff = 1.048e-1 |
| B: jointly scaled prior restores invariance | PASS | rel Frobenius + rho | rel F = 1.037e-15; max |rho| = 1.221e-15 |
| C: genuinely distinct frequency moves confounded directions | PASS | rho movement, 3 directions | 3/3 > 1e-4; movements 0.1213 / 0.2239 / 0.4768 |
| C: duplicate frequency leaves rho fixed | PASS | rho movement, 3 directions | max movement 2.90e-16 < 1e-10 |
| D: equal-budget trajectory hypothesis (retained circle360 > arc180 > arc90 > straight) | FAIL (counterexample, not forced) | observed retained_mass order | arc90 > straight > circle360 > arc180 (retained_mass 9.0388 / 8.3535 / 8.1395 / 7.9993); `hypothesis_holds = false`, `forced_pass = false` |
| D: confusable-mass order (inverted hypothesis) | FAIL (counterexample, not forced) | observed confusable_mass order | arc180 > circle360 > straight > arc90 |
| 4b: block-normalized frequency control | PASS (conclusion unchanged) | normalized distinct/duplicate/monotonicity | normalized movements 0.1060 / 0.1944 / 0.3815 (all > 1e-4, 3/3); duplicate max movement 5.23e-16 (< 1e-10); normalized PSD monotonicity passes 4/4 |
| 4-parent: SNR-matched frequency diversity | PASS | per-frequency nominal signal-norm matching | three selected retentions become 0.1170 / 0.2157 / 0.4460 |
| 4-parent: trajectory rankings are metric/standoff-sensitive | observed | equal-length (identity vs SNR) and same-radius curved-path controls | equal-length retained SNR order arc90 > straight > arc180 > circle360; same-radius (R=1.6) retained order arc90 > arc180 > circle360 under both metrics; no universal ordering |

### 1.6 Family 5 - sensitivity, robust surrogate, rank events

| claim | recorded status | executed check | key numbers |
|---|---|---|---|
| P1: full derivative formula (incl. dW/dB pose-Hessian term) vs FD | PASS | rel_full at eps=1e-3 and clean-window slope | rel_full = 6.25e-7 < 1e-4; slope 1.91 (r^2 = 0.998) |
| P1: frozen-nuisance approximation clearly worse | PASS | rel_frozen at eps=1e-3 > 1e-4 | rel_frozen = 4.44e-2 (slope ~ -0.00) |
| P1: dW term dominance | observed | ||A0^T dW A0||_F / ||DK_full||_F | share = 4.44e-2 |
| P2: first-order eigenvalue prediction on well-gapped eigenvalues | PASS | rel_full_eig(1e-3) + slopes | idx23 (largest) 8.03e-8; idx22 (well-gapped) 1.13e-7; full slopes ~2.00 (gated) |
| P2: frozen prediction worse on gapped eigenvalues | PASS | frozen/full ratio at 1e-3 | idx23 39.4x; idx22 14.8x |
| P2: median eigenvalue gap flag | flagged (not asserted) | min adjacent gap < 1e-6 | idx12 min gap 5.94e-8 < 1e-6 -> excluded from P2 gate; `any_chosen_gap_lt_1e-6 = true` |
| P3: robust first-order surrogate | PASS | err at eps=1e-2, decreasing, sampled-min | err(1e-2) = 4.46e-5 < 5e-2; slope 2.00; sampled_min >= adverse - slack at all three eps |
| P4: empirical Weyl lower-bound self-check | PASS (empirical only) | 150 samples, 18 basis directions | violations = 0; L_emp = 2.109e-4; L_nominal = 2.242e-4; L_i mean = 1.248e-4; no uniform operator-Lipschitz certificate |
| P5: projector movement / rank scan | observed (no forced pass) | one random path, tau scan | min adjacent eig gap = 6.62e-15 at tau* = 0.05 (pair (1,2)); rank events false; sigma_min(B) event false (sigma_min(B) ~ 6e-4..8e-4); crossing rel err 1.32e-6 vs gapped 1.13e-7 (ratio 11.7x) |

## 2. Key numbers (exact values from the result JSONs)

| family | quantity | value(s) as recorded in JSON | source |
|---|---|---|---|
| 1 | map-Jacobian FD rel error at h=1e-3 (seeds 0/1/2) | 2.246e-09 / 1.740e-08 / 1.766e-09 | family1_pilot_results.json |
| 1 | pose-Jacobian FD rel error at h=1e-3 (seeds 10/11/12) | 8.382e-07 / 4.279e-06 / 4.087e-06 | family1_pilot_results.json |
| 1 | state margin | sigma_min(M)/||M|| = 0.7189; max state residual 9.34e-16 | family1_pilot_results.json |
| 1 | Born/full-wave discrepancy range (s = 0.01..0.2) | rel F 1.77e-3..3.54e-2; rel A Fro 1.89e-3..3.82e-2 | family1_pilot_results.json |
| 1 | grid refinement rel diff vs N=40 | 1.453e-3 (16), 5.030e-4 (24), 1.601e-4 (32) | family1_grid_refinement.json |
| 1 | self-cell quadrature | max abs quadrature error 3.13e-16 | self_cell_quadrature_validation.json |
| 2 | c1 pixel projector gate | subspace distance 8.4658e-07 vs 1e-8 -> `pass_gate=false`, `identity_support=true` | family2_results.json |
| 2 | c3 retention residuals | 2.8866e-15 (pixel) / 1.5543e-15 (smooth) | family2_results.json |
| 2 | stable kernel factorization residual | 3.773e-16 (pixel) / 2.703e-16 (smooth) | family2_parent_correction.json |
| 2 | small-prior limit (SVD-filter, alpha=1e-14) | d0 = 5.6985e-09 (pixel) / 3.6702e-09 (smooth) | family2_parent_correction.json |
| 2 | singular-prior large-alpha non-convergence | dinf_sing(1e4) = 4.715e-04 / 2.296e-04 (417-912x regular) | family2_results.json |
| 2 | interlacing max positive violation | 5.0e-21 (pixel) / 6.7e-21 (smooth) | family2_results.json |
| 2 | duplicate fixed-prior rho change (alpha=1) | 3.99e-2 (pixel) / 3.98e-2 (smooth) | family2_results.json |
| 3 | smooth gauge residuals (Tx/Ty/Rot, N16) | 0.18289 / 0.08648 / 0.31016 | family3_gauge_born.json |
| 3 | smooth representation residuals | 0.7810 / 0.5719 / 0.7030 | family3_gauge_born.json |
| 3 | pixel gauge residuals (Tx/Ty/Rot, N16) | 6.6131e-05 / 1.5787e-05 / 4.6820e-05 | family3_gauge_born.json |
| 3 | anchor rho_min ratios (corner/disk/boundary) | 1.26x / 0.011x / 0.244x | family3_gauge_born.json |
| 3 | Born empty-background B=0 and K_SLAM0=K_IS0 | both relative Frobenius 0.0 | family3_gauge_born.json |
| 3 | raw full-wave Born rel_err at h=1e-3 | 0.21828 (h-independent) | family3_gauge_born.json |
| 3b | pixel residual plateau N16/32/40 | Tx 6.613e-05/7.955e-05/8.135e-05 (similarly Ty, Rot) | family3b_refinements.json |
| 3b | augmented p27 representation residuals | 1.42e-15 / 2.04e-15 / 1.63e-15 | family3b_refinements.json |
| 3b | corner-anchor generator ratio | 11.93 | family3b_refinements.json |
| 3b | A0 vs Born relative Frobenius | 7.05e-17 global; 7.09e-17 max per pose | family3b_refinements.json |
| 3b | full-wave chi^2 floor | 2.166e-4 (h-independent); remainder slope 1.986 | family3b_refinements.json |
| 4 | PSD increment min eigenvalues (F1->2, F2->3) | K_SLAM 1.517e-13 / 1.149e-13; K_eff 1.288e-14 / 5.350e-14 | family4_results.json |
| 4 | distinct-frequency rho movements | 0.1213 / 0.2239 / 0.4768 (3/3 > 1e-4) | family4_results.json |
| 4 | duplicate-frequency max rho movement | 2.90e-16 | family4_results.json |
| 4 | trajectory retained_mass order | arc90 (9.0388) > straight (8.3535) > circle360 (8.1395) > arc180 (7.9993) | family4_results.json |
| 4 | trajectory confusable_mass order | arc180 (16.0007) > circle360 (15.8605) > straight (15.6465) > arc90 (14.9612) | family4_results.json |
| 4b | normalized distinct movements | 0.10596 / 0.19444 / 0.38146 (3/3 > 1e-4) | family4b_normalized_control.json |
| 4-parent | SNR-matched retentions | 0.1170 / 0.2157 / 0.4460 | family4_parent_controls.json |
| 4-parent | same-radius curved-path retained order | arc90 > arc180 > circle360 (identity and SNR) | family4_parent_controls.json |
| 5 | P1 rel_full vs rel_frozen at eps=1e-3 | 6.2467e-07 vs 4.4422e-02; full slope 1.91 | family5_sensitivity.json |
| 5 | P2 gapped eigenvalue rel errors at 1e-3 | 8.03e-08 (idx23), 1.13e-07 (idx22); frozen ratios 39.4x/14.8x | family5_sensitivity.json |
| 5 | P3 surrogate error at eps=1e-2 | 4.463e-05; slope 2.00; sampled_ok = true | family5_sensitivity.json |
| 5 | P4 empirical Lipschitz | L_emp = 2.1088e-04, L_nominal = 2.2415e-04, violations = 0 (empirical only) | family5_sensitivity.json |
| 5 | P5 minimum adjacent gap | 6.624e-15 at tau*=0.05, pair (1,2); rank events false | family5_sensitivity.json |

## 3. What these checks cannot establish

* All five families are finite-dimensional, dense, double-precision checks on
  the specific discrete N=16 (with N=24/32/40 diagnostics where recorded)
  whitened/realified model built from the implemented 2D scalar Helmholtz
  contrast-source solver. They establish no continuum-limit theorem, no
  infinite-dimensional operator statement, and no closed-range or stable
  transversality property of the continuum inverse problem.
* No global nonlinear convergence or recovery guarantee follows from the
  linearized/Jacobian-level checks; cycle skipping, phase wrapping, posterior
  multimodality, wrong data association, and basin-of-attraction behaviour are
  outside every executed gate.
* No online SLAM performance claim is made. Runtime and conditioning records
  describe offline dense algebra on Apple Silicon CPU; they are not a
  navigation/tracking system test.
* No global novelty claim is made. Related work is cited for retrieval
  context; wording such as "first" or "no prior work" is not used.
* Family 2's pixel-basis c1 projector distance did not meet the predeclared
  1e-8 gate; the algebraic identity is certified by dimension, residual, and
  stable-factorization evidence, not by that eigenbasis distance, and the gate
  failure is preserved.
* Family 3's raw smooth-basis gauge and Born-gate failures are preserved
  records; 3b shows which corrected conditions pass but does not erase or
  relabel those failures. Pixel residuals (~8e-5 at fine N) are a discrete
  forward-model gauge floor, not algebraic zero, and are not monotonically
  decreasing in N.
* Family 4's trajectory ordering is a reproducible scenario counterexample
  under the raw identity-noise metric; 4b and the parent controls demonstrate
  metric/standoff sensitivity, not a universal angular-coverage theorem.
* Family 5's P4 Lipschitz numbers are empirical-only; no uniform
  operator-Lipschitz constant over the full perturbation ball was derived, and
  the Weyl self-check passes by construction on the same samples. Rank events
  were not encountered along the tested path; no rank-event theorem is
  claimed.

## 4. Related-work boundary

The executed ingredients have mature prior art that must not be claimed as
new: Schur-complement (pseudoinverse) elimination of nuisance parameters in
bilinear/self-calibration problems (S. Ling and T. Strohmer,
"Self-Calibration and Bilinear Inverse Problems via Linear Least Squares,"
arXiv:1611.04196, 2016/2017); increasing-stability results for multi-frequency
inverse source scattering in 2D (M. N. Entekhabi and V. Isakov, "On increasing
stability in the two dimensional inverse source scattering problem with many
frequencies," arXiv:1712.08696, 2017); multi-frequency microwave inverse
scattering of realistic volumetric numerical breast phantoms (J. D. Shea,
P. Kosmas, S. C. Hagness, and B. D. Van Veen, "Three-dimensional microwave
imaging of realistic numerical breast phantoms via a multiple-frequency
inverse scattering technique," Medical Physics 37(8), 2010); variational
source conditions and stability estimates for inverse electromagnetic medium
scattering (F. Weidling and T. Hohage, arXiv:1512.06586, 2015); and 4D-radar
SLAM systems where map/pose estimation and loop closure are engineered in
practice (K. Burnett et al., "Are We Ready for Radar to Replace Lidar in
All-Weather Mapping and Localization?," arXiv:2203.10174, 2022; M. Hilger
et al., "Towards Introspective Loop Closure in 4D Radar SLAM,"
arXiv:2404.03940, 2024). The synthesis tested here is narrower: the explicit
finite-dimensional pose-confounding retention theory for a whitened/realified
full-wave contrast-source SLAM model, together with its numerical consequences
(kernel/rank identities, principal-angle retention spectra, prior and duplicate
algebra, gauge and Born-empty-background controls, frequency/trajectory
comparisons, and first-order sensitivity with preserved counterexamples).
