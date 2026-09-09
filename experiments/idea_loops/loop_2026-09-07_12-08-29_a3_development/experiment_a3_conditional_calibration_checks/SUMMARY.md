# A3 conditional calibration checks - bounded numerical validation

## Scope
Deterministic, single-thread NumPy/SciPy falsification checks for the four A3
algebra statements (fixed-chart derivatives, complex gain quotient structure,
independent-Gaussian conditional branch bound, and dual-residual/resolvent
control). This is an algebra/numerics fixture pass only; it is NOT independent
3D verification, is NOT SOM-performance evidence, and makes no finite-noise
ratio-sufficiency or global-coverage claim.

## What passed
- B1: state-exact counterexample (1A), fixed-chart tangent-span sufficiency
  matching centered finite differences (1B), nonzero-residual coefficient
  derivative including (M_v U)^* r_s (1C), and the moving-U control in which
  U_v must be included (1D).
- B2: cross-ratio gain invariance (2A), factor recovery (2B), complex kernel
  L=row+column gain tangent (2C), the realified whitened quotient projector
  identity under a non-identity SPD covariance with cond 5 (2D), rank-one
  hiding (2E), and the two-component closed form (2F).
- B3: conditional two/three-candidate branch bounds over 10000 independent
  draws each (seeds 3101-3110); empirical errors respect the theorem bounds.
- B4: complex dual-residual identity (4A), SVD-certified resolvent bound
  (4B), and all exact numerical identities in the near-singular example (4C).

## Expected-fail controls (must fail; all reported as expected_fail)
- B2 2G: H with an exactly zero denominator is flagged instead of divided.
- B3 3D: missing candidate (truth 2 e1 not in {0, 4 e1}) gives empirical
  error 1 by definition; theorem bound not applicable.
- B3 3E: refitting a candidate to validation y gives empirical error 1 while
  the plug-in bound evaluated with the refitted distance is invalid.
- B4 4C control: tiny residual (eps=1e-6) with output error 1 shows small
  residual is not a certificate without a verified resolvent gamma.

## Files
- regenerate_checks.py: self-contained regenerator (all cases, fixed seeds,
  deterministic fixture searches where contrast is required).
- checks.json: case-level inputs/metrics/thresholds/status/notes plus counts.
- SUMMARY.md: this summary.
Written to both research/delegated/a3_algebra (canonical) and the inspection
directory of the bounded task.

## Limitations and remaining work
- Finite-difference comparisons use h=1e-6 (5e-4 relative thresholds); they
  validate algebra, not physical-model truth or convergence.
- B2D identity is algebraic; noisy nonlinear cross-ratio statistics remain
  correlated/non-Gaussian and are not claimed sufficient.
- B3 bounds are conditional on branch coverage; no global feasible-set
  coverage or multi-start hypothesis-bank guarantee is claimed.
- Remaining work: independent 3D/forward-model verification, finite-noise
  calibration tests, and A3 manuscript integration by the parent Codex.
