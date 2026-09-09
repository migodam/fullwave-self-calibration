# Family 3: gauge and Born empty-background degeneration - claim-status report

Date: 2026-09-03 (UTC run 09:59:08; SGT ~17:59).
Experiment: `experiment_pose_confounding_spectral_geometry`.
Families 1-2 were not rerun.  Family 3 consumes the validated
`helmholtz.build_AB` / `whiten_realify` / `forward_measurements` outputs and
reuses the Family-2 smooth RBF basis and machine rank-tolerance conventions.

## Scope and method

`results/family3_gauge_born.json` was produced by the exact command (from the
experiment root):

```bash
.venv/bin/python src/family3_gauge_born.py
```

Wall runtime 4.81 s (recorded in the JSON).  Environment (JSON `platform`):
Apple Silicon CPU (`arm64`), Python 3.12.13, numpy 2.5.2, scipy 1.18.1,
matplotlib 3.11.1; no GPU/MPS/CUDA.  The only RNG use is the fixed-seed 31415
pose direction in the Born bilinear check.

The whitened/realified data space is m = 2*T*n_rx = 48 and q = 3*T = 18.
Checks use the analytic two-blob scene, the analytic symmetric scene
`0.5 exp(-|r|^2/(2*0.12^2))`, and the empty (`chi0 = 0`) background.
Self-cell marker: `hh.SELF_CELL_FORMULA` / version are stamped into the JSON
(self-cell v2, corrected 2026-09-03); no self-cell formula is recomputed here.

Rank tolerance throughout: `tol(M) = max(M.shape)*eps_machine*sigma_1(M)`
(Family-2 convention).  Gauge directions are evaluated literally as the
analytic pixel derivative `dchi_g`, with smooth coefficients
`c_g = lstsq(S, dchi_g)` and pixel coefficients `c_g = dchi_g`.

## Tolerances used

| gate | value |
|---|---:|
| smooth two-blob gauge residual (each Tx/Ty/Rot) | <= 5e-2 at N=16 and decreasing to N=40 |
| Born B0 Frobenius ratio | < 1e-12 |
| Born K_SLAM0 == K_IS0 relative Frobenius difference | < 1e-12 |
| Born rel_err at h = 1e-3 | < 1e-3 |
| Born rel_err small-h fitted slope | in [1.8, 2.2] |
| anchor restricted rho_min | > 10 x base rho_min (each of the three restricted masks) |

## Claim-status table

| claim | status | executed check | key numbers |
|---|---|---|---|
| SE(2) generators cancel in smooth basis (gauge residual <= 5e-2, decreasing) | **FAIL** | gauge residual at N=16 and N in {16,32,40} | Tx 1.829e-1 -> 1.828e-1; Ty 8.648e-2; Rot 3.102e-1 -> 3.093e-1 (all well above 5e-2; representation residuals 0.781 / 0.572 / 0.703) |
| smooth generator lies in ker K_SLAM (kernel-distance check) | **NOT MET at a gauge level** | ||K_SLAM c_g||/||c_g|| | 1.91e-8 (Tx), 2.31e-9 (Ty), 1.43e-8 (Rot); nonzero, consistent with the representation-limited c_g |
| pixel basis is a documented larger-residual control | **Expected direction not reproduced** | pixel gauge residual | 6.61e-5 (Tx), 1.58e-5 (Ty), 4.68e-5 (Rot); 1000-20000x SMALLER than smooth residuals.  Pixel does not certify an exact discrete equivariant gauge, but its numerical cancellation is far better because Range(A_pix_R) contains Range(B_R) (Family 2: r_AB = 48, d_inter = 18) |
| symmetric scene rotation is not a gauge direction | **PASS** | dchi_w norm, gauge residual | ||dchi_rot|| = 7.8e-17 (numerically zero), gauge residual = 1.0000 (not small); Tx/Ty residuals 3.66e-2 / 1.33e-1 |
| retention spectra (smooth) recorded | PASS (recorded) | R_op = Q_A^T (I - Z Z^T) Q_A | two-blob rho_min = 5.51e-11; symmetric rho_min = 2.20e-11; 24 eigenvalues each |
| anchoring / known background / fixed boundary removes the translation-x gauge (rho_min > 10 x base) | **FAIL** | restricted smooth A_red spectra | base 5.51e-11; corner anchor 6.97e-11 (x1.26); known-background disk 6.16e-13 (x0.011); outer-ring boundary 1.34e-11 (x0.244).  None exceeds 10x base; disk/boundary rho_min actually decrease |
| B = 0 at empty background | **PASS** | ||B0_R||_F / max(1, ||A0_R||_F) | 0.0 < 1e-12 |
| K_SLAM = K_IS at empty background | **PASS** | ||K_SLAM0 - K_IS0||_F / ||K_IS0||_F | 0.0 < 1e-12 |
| bilinear term dominates the linear Born term | PASS (metric) | ||BIL|| / ||A0_R dchi|| at h=1e-3 | 2.81 |
| full-wave FD matches the bilinear term at rel_err < 1e-3 and slope ~2 | **FAIL** | rel_err(h), clean-window slope | rel_err(1e-3) = 2.18e-1; rel_err is h-independent (clean slope 8.6e-8, first-four slope 1.5e-8), not O(h^2) |
| remaining second-order remainder R(h)-R0 has slope ~2 | **NOT MET as O(h^2)** | ||R(h) - R0|| vs h | slope 1.015 (linear in h); ||R(1e-3)-R0|| = 7.14e-6; ||R(1e-3)|| = 2.51e-3 and R0 = 2.51e-3 |

## Why the anticipated smooth gauge / Born gates did not pass

### Smooth gauge representation

The prescribed smooth basis is the 24-column unit-2-norm Gaussian RBF family
used in Family 2, but the analytic two-blob derivative fields are not in its
span: `||S c_g - dchi_g||/||dchi_g|| = 0.78` (Tx), `0.57` (Ty), `0.70` (Rot).
The map residual therefore cannot be an algebraic roundoff floor.  The Family-2
rank identity also explains why the pixel basis outperforms it here:
Range(A_smooth_R) has dimension 24 and Range(B_R) dimension 18 with
r_AB = 42 (intersection dimension 0), so no nonzero `A_smooth_R c` can exactly
equal `-B_R dX_g`; Range(A_pix_R) fills the 48-dimensional data space
(r_AB = 48, intersection dimension 18), so cancellation is possible up to the
discrete forward-model floor (observed 1.6e-5..6.6e-5).  This is recorded as
the pixel-basis result; it does not certify exact equivariance, and the prompt's
expectation of a larger pixel residual was not reproduced on this scenario.

### Born bilinear check

The bilinear term itself is dominant against the linear term (2.81x at
h=1e-3) and against the R0 remainder (`R0` is 21.9% of the linear norm, and
R(h)/||BIL|| is 7.8% at h=1e-3).  However, `FD_h` is the finite-difference
derivative of the **full-wave** map at the fixed unit-L2 scene `dchi`
(max |dchi| = 0.346), so it contains `chi^2`-order pose derivatives that are
not part of `(D_X A(0)) dchi`.  These produce an h-independent mismatch floor
(rel_err ~0.218) rather than an O(h^2) finite-difference error; accordingly
`R(h)-R0` is dominated by a term linear in h (fitted slope 1.015).  A
Born-linearized proxy would trivially pass the O(h^2) gate, but that is not
the specified `hh.forward_measurements` check.  These raw values establish
that for this O(1) unit-L2 scene amplitude the strict 1e-3 bilinear-dominance
gate is not met; the data do not support claiming pose errors are harmless.

### Anchor variants

Only the upper-right 3x3 corner block (x > 0.30 and y > 0.30) was fixed by the
specified corner anchor, and that region lies in the tails of both blobs.
Its restricted rho_min is essentially unchanged (ratio 1.26).  The disk and
outer-ring masks restrict many pixels with meaningful field sensitivity but
reduce rho_min instead of increasing it, so the expected "anchors remove the
gauge" statement does not hold for these particular masks on this scene.

## Artifacts and digests

```text
0d3252da2c790fa3c841799818f8f49fd4213fce10e2c45ae04633a057c299e7  src/family3_gauge_born.py
464a625b3ade27a75a8ed862d75f8410799f2c1f2b3b0c615e0f6fb29784eb82  results/family3_gauge_born.json
89bb4a8bf21be6e2783e889aceb3de2b0d2bd2ce92c3bec9d2b5cabcb56e4383  figures/family3_gauge_residuals.png
b06303bf3bce45162c780bd8573a7bbd02f565893bec83f26f11d34aa76a1d3d  figures/family3_gauge_refinement.png
b4bbfdf386311426a33baf312b1ec5f04f8e64c458bc851046f1cb6867ac7473  figures/family3_anchor_rho.png
f63a84035a8b93d00bf4ab836056087c54500027e0e1e0b0e6d6a21401475876  figures/family3_born_empty_background.png
ee3b4ac45ad296dba533a8ccefd107a6cc53e4357b7a768d07f7db90f3ee6585  src/helmholtz.py
762af3bdd33af1ab80563c7257b3805e7972f786fa88af1bf5e1b51a768526d7  src/family1_pilot.py
173d2b5d9f66be274c9a8422f4ac2c381bb91a60132b37c734748414cf3eb235  src/family2_algebraic_spine.py
```

## Cannot-establish section

- Finite-dimensional, discrete, dense double-precision checks only; no
  continuum-limit, infinite-dimensional operator, or recovery/estimator claim.
- A bounded, boxed grid has no exact global SE(2) equivariance and this run
  does not construct an exactly equivariant discrete model, so no algebraic
  roundoff-level gauge floor was available.  The pixel gauge residuals
  (1.6e-5..6.6e-5) are the observed forward/discretization floor for the
  analytic pixel derivative, not a proof of exact covariance.
- In the smooth basis, representation residual is part of the floor by
  construction; these runs cannot establish that the SE(2) generators land in
  ker K_SLAM for this 24-column RBF family (kernel distances ~1e-8..1e-9 are
  nonzero).  They also cannot establish the converse (rho_min ~ 1e-11 => a
  physical group gauge), which remains open.
- The symmetric rotation dchi_w is at numerical zero (7.8e-17), so its
  relative representation residual is only meaningful with the recorded
  eps guard; the rotation generator is degenerate, not a small-but-finite
  gauge direction.
- The strict Born bilinear-dominance statement at rel_err < 1e-3 cannot be
  established at this fixed unit-L2 dchi because full-wave chi^2-order terms
  create an h-independent FD mismatch.  The claim should only be evaluated
  at a genuinely small scene amplitude or against the linearized Born map.
- The anchor statement cannot be established for the three specified masks on
  this scene; it is not a statement about other anchor geometries or other
  generators.
