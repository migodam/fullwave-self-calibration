# a2_physics physical-core check summary

Timestamp: 2026-09-06T14:01:27+08:00

## Status

Passed 8/8 bounded physical-core checks.

- **config_validation**: PASS (0.00s)
- **outgoing_convention_far_field**: PASS (0.00s)
- **self_cell_and_k2_factors**: PASS (0.05s)
- **forward_consistency_and_work**: PASS (0.03s)
- **analytic_A_B_vs_fd_both_apertures**: PASS (2.51s)
- **noise_real_variance**: PASS (0.01s)
- **rigid_distance_and_anchor_invariance**: PASS (0.03s)
- **adjoint_gradient_vs_fd**: PASS (0.17s)

## N16/N32 forward+jacobian timing (one call each, full aperture)

| N | rows | wall s | LU factorisations | RHS solves | max state residual |
|---|------|--------|-------------------|------------|--------------------|
| 16 | 288 | 0.169 | 4 | 312 | 9.21e-16 |
| 32 | 288 | 2.421 | 4 | 312 | 1.77e-15 |

## Facts

- State equation residual is at machine precision for every frequency/pose/illumination state (max rel residual 1.8e-15).
- Raw total equals the retained direct-incident plus scattered components by construction and in the checks.
- Analytic material A and pose B match centred finite differences at nonzero shared pose error for both full and limited apertures across all four frequencies (worst normalised errors in checks.json).
- The self-cell diagonal reproduces radial Gauss-Legendre quadrature of the equal-area disk integral to ~1e-15 relative error.
- One LU factorisation per frequency is shared by all 6 poses/illuminations; every LU right-hand side is charged individually.
- The efficient adjoint gradient matches finite differences of the realified objective and charges 12 forward + 12 adjoint RHS for two frequencies, with no additional factorisations.
- Realification sqrt(2)/sigma [Re; Im] makes proper complex noise isotropic with unit real covariance in the check.
- Pose 1 is anchored: its B block is exactly zero; direct-incident derivatives vanish under rigid receiver/transmitter co-motion.

## Artifacts

- `physics.py` (core), `test_physics.py` (executable checks)
- `checks.json`, `environment.json`, `summary.md`

## Uncertainties / assumptions

- 'Width 0.16 m' is interpreted as the Gaussian standard deviation in exp(-|x-c|^2/(2 width^2)), matching the unnormalised Gaussian basis convention used in the surrounding loop code; columns are not normalised.
- The shared unknown x is interpreted as adding the same (dx, dy, dtheta) world transform to the two non-anchored nominal poses.  The documentation states this explicitly; a different interpretation would change B but not the solver core.
- These checks validate the implemented discrete model's internal consistency (finite differences, quadrature, state residuals).  They are not continuum convergence, measurement-truth, or scientific acceptance evidence.
- No seeds 1001+ or optimiser runs were executed.

## Environment

- Python 3.13.13, NumPy 2.5.2, SciPy 1.18.1
- macOS-26.6.2-arm64-arm-64bit-Mach-O
