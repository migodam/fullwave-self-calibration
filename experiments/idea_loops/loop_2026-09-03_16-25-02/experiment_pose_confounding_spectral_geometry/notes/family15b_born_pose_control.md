# Family 15b: complete Born pose-tangent control

This additive parent correction closes the semantic gap in Family 15. Family 15's `(A_born | B_full)` rows remain a useful **mixed-tangent isolation control**; they are not a complete Born-SLAM Fisher model. Family 15b constructs the actual analytic Born pose Jacobian and compares complete tangent pairs.

## Execution

- Command: `.venv/bin/python src/family15b_born_pose_control.py`
- Runtime: 0.991 s on Apple Silicon CPU (macOS-26.6.2-arm64-arm-64bit)
- Dimensions: N=16, 48 real data rows, p=24 smooth map coefficients, q=18 pose coordinates
- Scope: finite-dimensional N=16 scalar Helmholtz toy; local tangents only; no continuum, estimator, hardware, or external-benchmark claim

## Analytic Born pose derivative

For each pose and coordinate, the implemented exact discrete formula is

`B_born,t[:,ell] = (D_X G_S,t[ell]) (chi * E_inc,t) + G_S,t (chi * D_X E_inc,t[ell])`.

At contrast scale s=0.1, centered finite differences over all 18 pose coordinates give:

| h | relative Frobenius error |
|---:|---:|
| 1.0e-03 | 1.756660e-05 |
| 3.0e-04 | 1.581003e-06 |
| 1.0e-04 | 1.756672e-07 |

The observed error-vs-step slope is 2.0000 (R^2=1.000000), supporting the analytic derivative with the expected centered-difference second-order trend on this discrete model.

## Exact finite-dimensional contrast-scaling proposition

For the linear Born map `F_B(chi,X)=A0(X) chi`, set `chi=s chi_bar`. Then `A_B(s)=A0` and `B_B(s)=s B1`. Hence, for every nonzero s, `Range(B_B(s))=Range(B1)`: the no-prior projector and Born `K_SLAM` are independent of contrast amplitude. At s=0, `B_B(0)=0`, so `K_SLAM(0)=K_IS`. The endpoint is discontinuous unless the map tangent is orthogonal to the nonzero Born pose range.

With `J=alpha I`, alpha>0, the finite-prior loss is continuous and obeys

`||K_IS-K_eff(s)||_2 <= s^2 ||A0^T B1||_2^2 / alpha`.

This follows directly from differentiation of the bilinear Born map, invariance of a range under nonzero scalar multiplication, and `(s^2 B1^T B1+alpha I)^-1 <= alpha^-1 I`. A semidefinite prior needs a separate range/kernel analysis.

Executed checks: max `A_B(s)` scale residual 0.000e+00; max `B_B(s)/s` residual 5.530e-16; max nonzero-s Born projector residual 4.521e-14. All finite-prior bound rows hold: True.

The nonzero-s complete-Born retained DOF is 9.15111471 at s=0.001 and 9.15111471 at s=1; the s=0 endpoint is 24.00000000. The tiny nonzero-s variation is numerical. This is a rank-stratum singular limit, not evidence that arbitrarily weak scatterers provide finite practical pose information.

## Complete-pair versus mixed-pair result

| s | rel A | rel B | rel K_IS | rel K_eff complete | projector gap | rdof full | rdof Born complete | rdof mixed |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.001 | 4.999e-04 | 1.765e-04 | 7.018e-04 | 7.018e-04 | 1.658e-04 | 9.151611 | 9.151115 | 9.151173 |
| 0.01 | 5.004e-03 | 1.765e-03 | 7.026e-03 | 7.026e-03 | 1.660e-03 | 9.156332 | 9.151115 | 9.151696 |
| 0.03 | 1.504e-02 | 5.297e-03 | 2.113e-02 | 2.113e-02 | 4.986e-03 | 9.168423 | 9.151115 | 9.152865 |
| 0.1 | 5.047e-02 | 1.767e-02 | 7.101e-02 | 7.100e-02 | 1.671e-02 | 9.224559 | 9.151115 | 9.157027 |
| 0.3 | 1.536e-01 | 5.316e-02 | 2.156e-01 | 2.156e-01 | 5.088e-02 | 9.357869 | 9.151115 | 9.169489 |
| 1 | 5.049e-01 | 1.784e-01 | 6.527e-01 | 6.520e-01 | 1.776e-01 | 9.407162 | 9.151115 | 9.219363 |

Observed log-log slopes are A 1.003, B 1.001, K_IS 0.994, complete-pair K_eff 0.994, and pose-projector gap 1.008. These are finite-grid observations, not continuum asymptotics.

## Scientific boundary

- Only the complete `(A_born | B_born)` rows represent the implemented Born pose-elimination model.
- `(A_born | B_full)` remains explicitly labelled mixed-tangent.
- At s=0, first-order pose information vanishes; for s>0, a no-prior projector ignores the magnitude of B and is therefore a singular statistical idealization as s tends to zero.
- The Diong et al. paper (DOI `10.1088/0266-5611/32/6/065006`) is a conceptual neighbor only. This run does not reproduce its estimator bias, geometry, noise, parameterization, units, or numerical values and is not an external benchmark.

## Artifacts

- `results/family15b_born_pose_control.json`
- `figures/family15b_born_pose_control.png`
- `src/family15b_born_pose_control.py`
