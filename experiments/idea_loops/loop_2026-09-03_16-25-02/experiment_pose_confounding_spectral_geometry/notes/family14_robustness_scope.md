# Family 14: robustness-scope audit (demotion audit)

Date: 2026-09-03T16:54:55.609468+00:00 (UTC; also recorded in JSON).  Experiment: `experiment_pose_confounding_spectral_geometry`.  This family audits the scope of previously executed robustness evidence; it computes no new Helmholtz operator and runs no new forward model.

## Exact command and runtime

```bash
.venv/bin/python src/family14_robustness_scope.py
```

Wall runtime: 0.008 s.  Platform: macOS-26.6.2-arm64-arm-64bit.  Python 3.12.13, numpy 2.5.2, scipy 1.18.1, matplotlib 3.11.1.  Apple Silicon CPU only; no GPU/MPS/CUDA.

## Read-only reuse

Loaded artifacts (SHA-256 in the JSON): `results/family5b2_alpha_tight.json`, `results/family5_parent_generalized.json`, plus the two source scripts referenced only for name/structure context (`src/family5b2_alpha_tight.py`, `src/family5_parent_generalized.py`). No existing source, result, figure, or note file was modified.

## Layer A: exact finite-dimensional linear algebra

The affine tangent certificate is a finite-dimensional structural bound conditional on the implemented FD-validated Jacobians; it is NOT a continuum proof and NOT a full nonlinear Helmholtz certificate.

Stored pointwise structural constant: L_struct_tight(X0) = 1.124399e-02.  Stored affine certificates for eps in {1e-3, 3e-3, 1e-2}: 1.130863e-02, 1.143821e-02, 1.189479e-02.  Recomputation from the stored ingredients reproduces every stored value exactly at IEEE double precision (max abs relative diff 0.000e+00); the same formula and operation order as Family 5b2 was used.

| eps | recomputed L_cert_affine | stored L_cert_affine | abs rel diff |
| --- | ------------------------ | -------------------- | ------------ |
| 1e-03 | 1.130863e-02 | 1.130863e-02 | 0.000e+00 |
| 3e-03 | 1.143821e-02 | 1.143821e-02 | 0.000e+00 |
| 1e-02 | 1.189479e-02 | 1.189479e-02 | 0.000e+00 |

## Layer B: conditional/discretized full-wave differentiation

FD validation numbers (J_A/J_B relative Frobenius agreement 5.990607e-07 / 5.374392e-07, i.e. ~6e-7) support the implemented Jacobian but do not certify the continuum derivative.  Max per-column relative differences are 7.671831e-07 (J_A) and 7.613361e-07 (J_B); stored consecutive-difference ratios over h ~ 4.0 provide O(h^2) support for the centered-FD estimates only.

## Layer C: executed numerical evidence

Zero sampled violations over the stored sample counts at the tested eps cannot be called certification; the observed worst quotients and margins are reported below.  Full-nonlinear sampling was 500 unit directions per eps (seed 9191); affine sampling was 2000 unit directions per eps (seed 7171).

| eps | L_cert_affine | full-nonlinear worst ||K(eps u)-K(0)||_F/eps | margin | ratio | affine worst | affine margin |
| --- | ------------- | -------------------------------------------- | ------ | ----- | ------------ | ------------- |
| 1e-03 | 1.130863e-02 | 3.225995e-04 | 1.098603e-02 | 35.0547 | 3.528163e-04 | 1.095582e-02 |
| 3e-03 | 1.143821e-02 | 3.226528e-04 | 1.111556e-02 | 35.4505 | 3.541288e-04 | 1.108408e-02 |
| 1e-02 | 1.189479e-02 | 3.228342e-04 | 1.157196e-02 | 36.8449 | 3.744008e-04 | 1.152039e-02 |

At eps=1e-2: L_cert_affine=1.189479e-02, full-nonlinear worst sampled quotient 3.228342e-04 (max across the three stored full-model eps), margin 1.157196e-02, ratio 36.8449; affine worst sampled quotient 3.744008e-04, affine margin 1.152039e-02, affine ratio 31.7702.

## Layer D: open questions

The full nonlinear Helmholtz operator has NO interval/uniform certificate in this implementation.  The stored generalized-derivative formulas were validated on simple positive-gap modes and require locally constant rank, so they are inapplicable at rank events (see `results/family5_parent_generalized.json`).  Repeated eigenvalues require compressed derivatives; only cluster control evidence is available.

Generalized-derivative audit evidence (stored): finite-prior backward residual (factored) 4.040062e-16, normalization error (factored data) 6.182166e-12, FD relative errors at eps=1e-3 for three simple positive-gap modes 4.792926e-07, 3.919982e-09, 1.466669e-08, prediction-error slopes 2.001706, 2.001466, 1.999819.

No-prior (constant-rank) case: backward residual (factored) 9.238499e-16, normalization error (factored data) 4.074074e-12, FD relative errors at eps=1e-3 2.061912e-06, 1.606920e-05, 6.461263e-06, prediction-error slopes 2.004058, 1.995796, 1.992358.

Repeated-eigenvalue cluster control: the six structural rho=1 modes stay at one to <=1.44e-15 under the stored perturbed rows, with compressed-derivative spectral norm 2.297307e-10; the direct coefficient-space compressed norm is 1.264156e-06.

Engineered rank event (stored): rank(B) 18 -> 17 at t=0, projector jump operator norm 1.000e+00 (1.0), retained mass 9.407162e+00 -> 1.040528e+01 (9.40716 -> 10.40528 at the stored rows), and relative Frobenius information jump 0.13996 (0.13996).  This is an engineered algebraic control, not an observed physical trajectory event, and it falsifies constant-rank derivative applicability at the event.

## Explicit demotion

> All full nonlinear execution-error robustness claims are DEMOTED to empirical/structural observations; only the affine tangent certificate is retained as a finite-dimensional conditional certificate.

## No claims made

* No interval proof is claimed.
* No Hankel bound is claimed.
* No continuum robustness is claimed.
* No 'certified nonlinear robustness' is claimed.

## Artifacts

* [results/family14_robustness_scope.json](results/family14_robustness_scope.json)
* [notes/family14_robustness_scope.md](notes/family14_robustness_scope.md)

Source SHA-256: `074a7131a5853964e84bfe21ef69c6406c668a0816e264b587750d9de020a191`.

Artifact digests:
* `family5b2_json` `5f9e2eca2305106f3f5ab7d22fa0554a193fd8819da71128d18ecfe939a51bc8`
* `parent_generalized_json` `276c595d5dcf2e0ca08e87bbb12499be4c77158f3bab1cf38ec703ccc7ee3b7c`
* `family5b2_script` `a26bf62d177f897b75dfc1578822e464805c06c77db6c35e1bfd7f3b21620f3d`
* `parent_generalized_script` `fbdc8c2f11f6e3503f30f637ec65570062b2c3207d849e8bde1aa62236cd3e4d`
* `source_script` `074a7131a5853964e84bfe21ef69c6406c668a0816e264b587750d9de020a191`
