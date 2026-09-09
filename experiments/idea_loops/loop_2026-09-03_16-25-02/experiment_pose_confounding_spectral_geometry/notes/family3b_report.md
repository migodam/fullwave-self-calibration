# Family 3b: correctly-conditioned refinement report

Date: 2026-09-03 (SGT; UTC stamp in `results/family3b_refinements.json`).
Experiment: `experiment_pose_confounding_spectral_geometry`.

This report supplements (and does **not** supersede) the raw Family 3 record.
`results/family3_gauge_born.json` and `notes/family3_report.md` were left
unchanged.  The three Family 3 checks were specified under conditions that do
not match the theory's validity conditions:

1. smooth-basis gauge was tested against generator fields that the p=24 RBF
   basis cannot represent (representation residual 0.57-0.78);
2. the Born bilinear/remainder test used an O(1) `||dchi||_2 = 1` scene, so
   full-wave `chi^2`-order pose derivatives dominated the finite-difference
   comparison;
3. anchor efficacy was measured by the *global* `rho_min` of the full
   retention spectrum instead of generator-specific retention.

`src/family3b_refinements.py` reruns the same claims under corrected
conditions and writes both raw and refined numbers.

## Exact command and runtime

```bash
.venv/bin/python src/family3b_refinements.py
```

Wall runtime: 5.54 s (also in the JSON as `runtime_seconds`).  Platform:
Apple Silicon CPU (`arm64`), Python 3.12.13, numpy 2.5.2, scipy 1.18.1,
matplotlib 3.11.1.  Deterministic dense linear algebra; the only RNG use is
the fixed-seed-31415 pose direction in check D.

Scenario is unchanged: N=16 (A also runs N=32 and N=40), `k_b = 2*pi`, T=6,
n_rx=4, 90-degree arc radius 1.6, two-blob `chi0` from
`family1_pilot.make_chi0`.  Rank tolerance throughout:
`tol(M) = max(M.shape) * eps_machine * sigma_1(M)`.

## Tolerances used (Family 3b gates)

| gate | value |
|---|---:|
| A pixel gauge residual (all generators, all N) | <= 1e-3 and non-increasing from N=16 to N=40 |
| B augmented representation residual | < 1e-10 |
| B augmented gauge residual | <= 1e-3 |
| C base rho_g | <= 1e-8 |
| C corner-anchor rho_g / base rho_g | > 10 |
| D A0 == Born map, max per-pose relative Frobenius | < 1e-12 |
| D Born rel_err slope (clean window) | in [1.8, 2.2] |
| D rel_err_born at h=1e-3 | < 1e-4 |
| D full-wave rel_err at h=1e-2 | < 1e-2 |

## Claim-status table

| claim | status | executed check | key numbers |
|---|---|---|---|
| A: pixel-basis gauge residual small and decreasing/plateauing to roundoff as N increases | **FAIL on non-increasing condition; PASS on 1e-3 magnitude** | `||A_R dchi_g + B_R dX_g||/(||A_R dchi_g||+||B_R dX_g||)` for Tx/Ty/Rot at N=16/32/40 | Tx 6.613e-5 -> 7.955e-5 -> 8.135e-5; Ty 1.579e-5 -> 2.108e-5 -> 2.171e-5; Rot 4.682e-5 -> 6.067e-5 -> 6.239e-5.  All <= 1e-3, but N=32/40 are above the N=16 values, so the strict monotonicity gate is false |
| A: pixel-basis kernel distance tiny/decreasing | **PASS (recorded)** | `||K_SLAM dchi_g||/||dchi_g||` | Tx 3.96e-10 -> 1.22e-10 -> 8.05e-11; Ty 8.73e-11 -> 3.11e-11 -> 2.04e-11; Rot 1.85e-10 -> 6.47e-11 -> 4.27e-11 |
| B: original p=24 smooth basis reproduces Family 3 | **PASS (consistency)** | identical recomputation side by side | gauge residuals 0.182892 / 0.086484 / 0.310163, kernel distances 1.9e-8 / 2.3e-9 / 1.4e-8; max absolute gauge-residual difference vs `family3_gauge_born.json` = 0.0 |
| B: augmented p=27 smooth basis represents the generators | **PASS** | `||S_aug c_g - dchi_g||/||dchi_g||` | 1.42e-15 (Tx), 2.04e-15 (Ty), 1.63e-15 (Rot), all < 1e-10 |
| B: augmented smooth-basis gauge reaches the pixel forward floor | **PASS** | gauge residual with `A_smooth_aug = A_pix_R S_aug`, `B_R` | Tx 6.613e-5, Ty 1.579e-5, Rot 4.682e-5, all <= 1e-3 and equal to the pixel-basis residuals |
| C: unanchored translation-x generator has near-zero generator-specific retention | **PASS** | `rho_g = (c^T K_SLAM_red c)/(c^T K_IS_red c)` | base rho_g = 9.908e-11 <= 1e-8 |
| C: corner anchor removes the translation-x gauge | **PASS** | same rho_g after `free = NOT(x>0.30 & y>0.30)` | anchor rho_g = 1.182e-9, ratio = 11.93 > 10 |
| C: known-background disk and fixed outer ring (reported, not forced) | **PASS on ratio, reported honestly** | same generator-specific metric | known-disk rho_g = 2.165e-5 (ratio 2.19e5); outer-ring rho_g = 1.645e-8 (ratio 1.66e2).  Gauge residuals after restriction: base 6.61e-5, anchor 2.38e-4, known-disk 1.049e-1, outer-ring 7.80e-4 |
| D: A0 (full wave at chi=0) equals Born map | **PASS** | max per-pose-block relative Frobenius difference | 7.087e-17 < 1e-12 (global relative Frobenius 7.054e-17) |
| D: corrected small-amplitude Born rel_err slope ~2 | **FAIL / not established** | `||FD_born_h - BIL_h||/||BIL_h||`, h in [3.16e-3, 1e-1] | fitted slope -1.025 (roundoff-limited; see below).  rel_err_born at h=1e-3 = 2.265e-13 < 1e-4 |
| D: corrected full-wave error has small h-independent chi^2 floor | **PASS** | `||FD_full_h - BIL_h||/||BIL_h||` | rel_err_full = 2.166e-4 at every h from 1e-3 to 1e-1; fitted clean-window slope 5.4e-5 (floor-dominated), value at h=1e-2 = 2.166e-4 < 1e-2 |
| D: corrected second-order remainder drift | **PASS** | `||R_full(h) - R_full(0)||` over h in [3.16e-3, 1e-1] | fitted slope 1.986 (r^2 = 0.99999); at h=1e-2 value = 1.025e-8, R0 norm = 2.489e-9 |
| D: bilinear-to-linear dominance | **PASS (metric)** | `||BIL_h||/||A0_R dchi_small||` at h=1e-2 | 2.8087 (O(1), unchanged from the Family 3 ratio 2.81) |

## Why the raw Family 3 gates remain failed (record unchanged)

`results/family3_gauge_born.json` still records the raw failures:

* smooth-basis gauge residuals 0.183 / 0.086 / 0.310 remain above any
  meaningful tolerance because the p=24 RBF family cannot represent the
  derivative fields (raw representation residuals 0.78 / 0.57 / 0.70);
* the raw global-rho_min anchor comparison remained ratio ~1.26 for the
  corner anchor and decreased for disk/boundary;
* the raw O(1)-amplitude Born full-wave rel_err floor was 0.218, with
  remainder drift slope ~1 rather than ~2.

The corrected tests above do not erase those numbers.  They show which
claims become well-conditioned when (A) the pixel/discrete forward model is
used, (B) the smooth basis is augmented with the generator fields, (C)
retention is evaluated on the specific generator, and (D) the scene
perturbation is genuinely small (`||dchi_small||_2 = 1e-3`,
max |dchi_small| = 3.458e-4).

## Notes on the two gates that did not pass as specified

### A: pixel residual is not monotonically decreasing

The pixel-basis residual is small (all values 1.6e-5-8.1e-5, below 1e-3),
but N=32 and N=40 are each slightly **above** the N=16 value for every
generator, so the observed values do not satisfy the specified non-increasing
condition.  The residual is better interpreted as a small grid-converged
discrete forward-model gauge floor (~8e-5 at the fine resolutions) rather
than an exact-cancellation roundoff plateau.  The kernel-distance version of
the same gauge statement *is* tiny and decreasing with N.

### D: the Born `rel_err_born` slope gate is degenerate

Check D verifies that the full-wave `A` at `chi=0` equals the Born map
(7.1e-17 relative Frobenius).  `BIL_h` is then the centred finite difference
of that map, and `FD_born_h` is the centred finite difference of the same
linear Born map.  Their difference is therefore roundoff-limited (2.3e-13 at
h=1e-3 down to 1.7e-15 at h=1e-1) instead of exhibiting the O(h^2)
truncation slope requested.  The slope gate in [1.8, 2.2] cannot be
established from this pair and is reported as failed; the <1e-4 value gate at
h=1e-3 passes as an exact-equality consistency check.  The meaningful
second-order statement in the corrected regime is the full-wave remainder
drift `||R_full(h)-R_full(0)||`, whose slope is 1.986, and the small
full-wave chi^2 floor of 2.166e-4, which is h-independent as expected.

## Artifacts and digests

```text
6748930de3bb460ab841150e9b3ec61e4c9559ed24909b166985a94e4cc92a47  src/family3b_refinements.py
d4f28a59a578e4968aa20125fec6f652e8fc61716b56020b741851c14032eab4  results/family3b_refinements.json
5df2469564b0435354dee934450f994b5779429dd033c1e28dff4f28049efcc9  figures/family3b_gauge_pixel_refinement.png
18f0c7cd7a14ffcca2a23030b90fb9160a2282bcf2af28746e44489af98cec68  figures/family3b_anchor_generator_retention.png
c10183d208ae92074bc0dda038c27294cf2504db2ee478bd9a36d3bde07d281a  figures/family3b_born_bilinear_refined.png
ee3b4ac45ad296dba533a8ccefd107a6cc53e4357b7a768d07f7db90f3ee6585  src/helmholtz.py
762af3bdd33af1ab80563c7257b3805e7972f786fa88af1bf5e1b51a768526d7  src/family1_pilot.py
173d2b5d9f66be274c9a8422f4ac2c381bb91a60132b37c734748414cf3eb235  src/family2_algebraic_spine.py
0d3252da2c790fa3c841799818f8f49fd4213fce10e2c45ae04633a057c299e7  src/family3_gauge_born.py
```

## Cannot-establish section

* No continuum-limit, exact-global-SE(2), or recovery/estimator claim is made.
  The bounded grid is not exactly equivariant; the pixel gauge floor is
  observed (~8e-5 at fine N), not algebraic zero.
* A decreasing pixel gauge residual to roundoff is **not established**: the
  measured residuals rise slightly and plateau as N goes 16 -> 32 -> 40.
* The p=27 augmented basis contains the three generator fields by
  construction; it is a smooth 27-column tangent-space test, not proof that
  any smooth basis or the full infinite-dimensional problem has the same
  cancellation.
* The requested O(h^2) slope of `rel_err_born` cannot be established from the
  specified comparison because `FD_born_h` and `BIL_h` are centred finite
  differences of the same verified-identical map.  Establishing that slope
  would require a reference derivative that the specified checks do not
  provide.
* The known-background-disk and outer-ring rho_g ratios are reported raw
  (2.19e5 and 1.66e2); no theoretical claim is made that those geometries
  must retain the gauge or predict these exact ratios.
* The full-wave chi^2 floor (2.166e-4) is measured at one amplitude
  `||dchi||_2 = 1e-3`; its exact quadratic scaling in the scene amplitude was
  not swept here.
