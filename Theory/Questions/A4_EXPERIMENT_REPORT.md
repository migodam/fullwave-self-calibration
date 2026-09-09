# A4 executed experiments and numerical audit

This report is generated from preserved raw JSON by `analyze_a4.py`. It does not reproduce A3's frozen 240-fit experiment. Registration was local source/protocol hashing before fresh seeds, not an external preregistration service. The frozen files are never overwritten by the analysis script.

## Inventory

- 21 unit/regression cases pass with `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 PYTHONPATH=src python -m pytest -q`.
- 41 exact-spectrum evaluation points, z from0.01 to100; these are algebraic checks, not independent full-wave scenes.
- 96 DDA receiver records from 48 unique shape/material/frequency/grid state solves. The same state at two receivers must not be charged twice.
- 24 independent branch scenes, six fixed strata with four scenes each, 144 policy rows. One noise realization per scene, no200-replicate inflation.
- 8 new joint-calibration scenes, nine baselines, 72 endpoint rows.
- 54 component/illumination/gain rank controls, and a separate16-scene fresh reference correction audit.
- 602 post-freeze analytical exterior-model identity checks, labelled a supplement, not folded into frozen inference.

Main frozen collection wall time: 10.794197s, in the recorded single-thread environment. This excludes development, literature work, supplementary audit and manuscript generation. It is not a hardware acquisition cost or a large-scale full-wave performance claim.

## 1. Spectrum and exterior fidelity

Maximum vector-relative spectrum discrepancy: 6.163e-15. Near/far log slopes over the declared sampled endpoints: 0.999789, -2.001001, against theoretical+1/-2 limits. The error normalizes the three-singular-value vector; it is not a separate worst relative error for a tiny radial singular value.

![Modal information](../figures/01_modal_information.png)

The exact positive result is conditional on a radial target, controlled modes and scalar frame gain. Its optimum is fixed q and fixed k at kR=1; it is not a universal placement prescription at fixed transmitted power.

The post-freeze supplement verifies exact exterior full/radiative/static identities to6.66e-16 maximum absolute difference. Radiative-only and static-only approximations lose radial information after a free scalar gain. Their errors are analytically computable for this model family, unlike an arbitrary adjacent-grid difference.

## 2. DDA shape/model-class stress

Shapes: sphere axes(.05,.05,.05), near-sphere(.053,.05,.047), ellipsoid(.07,.045,.035), box half extents(.05,.04,.035)m. Permittivities1.5+.02i and6+.1i; k=6,18,36rad/m; grid spacings.025 and.05/3m; receiver radii.12 and.4m along a common oblique direction. This is48 state solves,96 receiver rows. Material differentiation profiles the allowed real-permittivity parameter together with complex gain.

| Shape | Minimum observed / ideal minimum singular value | Median ratio | Maximum ratio | Pattern residual range |
|---|---:|---:|---:|---:|
| box | 0.7685 | 0.9897 | 1.091 | 0.0261--0.3977 |
| ellipsoid | 0.8329 | 0.9692 | 1.105 | 0.0317--0.3063 |
| near_sphere | 0.8954 | 0.9723 | 1.027 | 0.008995--0.1011 |
| sphere | 0.6782 | 0.8861 | 1.013 | 0.0001755--0.0533 |

![DDA shape stress](../figures/02_shape_stress.png)

The ratio is a model-class stress diagnostic, not a validated lower bound for nonspherical targets. Cubic sphere discretizations are not exact radial media; the additional nuisance direction and discretization effects must not be erased. A small discrete linear residual is not a continuum error bound. The code is in-house DDA, not ADDA; treams was not installed/executed in A4. The analytic sphere solution supplies a distinct reference for the spherical case only.

## 3. Nine-baseline quantitative calibration

| Method | Median geometry (mm) | Median material (%) | Mean scaled task loss | Median charged wall time (s) |
|---|---:|---:|---:|---:|
| Low only | 1.978 | 6.29 | 5.355 | 0.01324 |
| Raw all-frequency | 9.947 | 26.82 | 64.006 | 0.01586 |
| Isotropic downweight | 2.447 | 7.397 | 6.8081 | 0.03957 |
| Rank-one downweight | 7.169 | 15.91 | 45.238 | 0.04229 |
| Sampled approximation error | 0.5424 | 0.4274 | 0.071704 | 0.03983 |
| Conditional linear error model | 0.5562 | 0.9181 | 0.18598 | 0.04022 |
| Coarse to fine | 0.6355 | 0.856 | 0.050674 | 0.06107 |
| Full fine | 0.6355 | 0.856 | 0.050674 | 0.04738 |
| Risk selector | 0.6355 | 0.856 | 0.050674 | 0.08468 |

Scaled task loss=(geometry error/.015m)^2+(real-permittivity error/.15)^2. Table geometry/material columns are medians, loss is a mean; do not reconstruct the mean from the medians. Every method uses matched observations, shared initialization, independent pilot, and the same additional noisy electronics reference. Low-only deliberately discards the third-frequency inversion block. The loss.05 is known, not reconstructed.

![Risk comparison](../figures/03_task_risk_baselines.png)

The risk selector chose full fine in 8/8 scenes and produced exactly the same parameter vectors (maximum difference 0.0). Median paired charged time ratio selector/fine is **1.861446**. Its descriptive paired bootstrap interval is 1.744032--2.055453; eight scenes and one serial timing each do not establish a broad population timing law. The structural negative result is more basic: identical endpoints plus additional construction/screening work.

![Charged controller cost](../figures/04_controller_cost.png)

**Controller decision: remove as a main contribution.** Conditional and sampled-error entries refer to these implemented models, not to optimal representatives of all Bayesian approximation-error methods. Sampled-error compensation performs well enough that raw/all-only comparisons would have been misleading. This calibration test uses a Rayleigh coefficient versus exact Mie coefficient with the exact exterior tensor retained; it is not a coarse/fine voxel-tomography benchmark.

## 4. Four distinct branch failures

| Policy | Accepted / 24 | Missing initial bank / 24 | Wrong selected / final covered | Wrong accepted / 24 | Correct rejected / correct estimates |
|---|---:|---:|---:|---:|---:|
| algebraic_multistart | 20/24 | 4/24 | 11/20 | 14/24 | 3/9 |
| branch_acquisition | 20/24 | 4/24 | 0/20 | 4/24 | 4/20 |
| coverage_aware | 16/24 | 4/24 | 0/20 | 4/24 | 8/20 |
| fisher_acquisition | 20/24 | 4/24 | 0/20 | 4/24 | 4/20 |
| naive_multistart | 20/24 | 15/24 | 0/9 | 14/24 | 3/9 |
| random_acquisition | 20/24 | 4/24 | 0/20 | 4/24 | 4/20 |

Definitions: missing bank means no candidate within15mm; wrong selection is conditioned on a final bank containing such a candidate; wrong acceptance is divided by all 24 scenes; correct rejection is divided by the number of correct point estimates, not all scenes. These heterogeneous designed strata are descriptive counts, not IID estimates of one natural scene population.

| Coverage-aware stratum | N | Accepted | Wrong accepted | Correct estimates rejected |
|---|---:|---:|---:|---:|
| antipodal | 4 | 4 | 0 | 0 |
| bank_missing | 4 | 4 | 0 | 0 |
| easy | 4 | 4 | 0 | 0 |
| low_snr | 4 | 0 | 0 | 4 |
| model_discrepancy | 4 | 0 | 0 | 4 |
| tangent_two_world | 4 | 4 | 4 | 0 |

![Coverage strata](../figures/05_coverage_strata.png)

Within the16 exact-modal strata, coverage accepted 12 and wrongly accepted0, rejecting 4. Even under an IID simplification,0/16 has a one-sided95% upper bound 0.1707, not1%; the theoretical noise-event statement is separate and conditional. Low-SNR cases exhausted the6000-cell budget and retained unresolved cells rather than declaring the search complete.

The4 deliberately shifted two-world data sets were **all wrongly accepted** by coverage. Their discrepancy violates beta=0 while matching an alternative admissible physical response across training, held-out and dither measurements. This defeats an unjustified fidelity declaration; it does not contradict the conditional set-membership theorem. Independent validation cannot distinguish identical data laws.

Coverage's overall wrong-accept count did not improve on strong Fisher/random/branch-acquisition controls. It rejected more correct estimates and incurs extra computation. It therefore remains a conditional safety mechanism, not a demonstrated superior branch method. The exact-modal +/- ambiguity is not a test of arbitrary many-cycle phase-wrapping branches.

## 5. Electronics-reference experiment and frozen diagnostic correction

Geometry is known in this separate test. With free gains at3 frequencies, material changes are exactly absorbable before a reference. Three noisy complex reference measurements constrain six real gain components. In the fresh corrected16-scene sample, the no-reference visible norm is at most 4.581e-06 (finite-difference residue), the minimum with-reference visible norm is 18.26, and median material relative error is 0.4573%. The profiled no-reference objective range is at most 1.208e-13.

![Reference diagnostic](../figures/06_reference_rank.png)

**Correction:** original `controls.py` used `J[:-6]` to remove reference rows from a realified array. The reference real/imag rows are noncontiguous; that diagnostic is invalid and withdrawn. Original recovery fits and raw files remain unchanged. `supplements/reference_audit_v2/` contains a corrected row-selection source, separately hashed local protocol, and fresh seed94032. These16 results are not pooled with the original16 to claim32 independent confirmations.

A later default pytest invocation collided with the installed third-party `coverage` package. Its traceback is retained in `postfreeze_default_pytest_import_failure.txt`. The documented test launcher disables auto-loaded plugins and explicitly prepends `src`; all 21 cases then pass. Frozen source hashes are unchanged.

## 6. Cost and evidence limits

Branch acquisition: initial27 complex training and27 validation values; each acquisition method adds9 training and9 fresh validation values. Joint methods use27 (low-only18) fitting values,18 independent pilot values,3 reference and3 independent pilot-reference values. Learned-error methods charge32 fine/coarse pairs per scene. Nominal geometry/support and relative polarization calibration are prior information, not free empirical discoveries.

Raw NPZ archives retain generated noisy arrays. Compute time, acquisition count and reference count are distinct units; no unsupported conversion to money or wall-clock hardware savings is made. No independent external laboratory measurement, hardware displacement recovery, external referee endorsement, full ADDA/treams reproduction, or global general-material certificate was completed.

**Verdict:** B is the only retained main-contribution candidate. A fails its strong-baseline gate; C is supporting. The work has not reached TAP method-paper acceptance.

## Baseline implementation boundary

The A4 rank-one baseline is the leading eigenmode of the sampled high-frequency approximation-error covariance. It is not a reconstruction of the unavailable original A3 pilot coarse/fine rank-one implementation. Its covariance trace need not equal the full sampled covariance trace; isotropic weighting uses the latter trace. The frozen comparison tests these explicitly defined methods, not the entire rank-one or Bayesian approximation-error families.
