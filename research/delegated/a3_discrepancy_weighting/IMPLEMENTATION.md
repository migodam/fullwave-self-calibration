# Frozen discrepancy-weighting implementation

Status: implemented for parent review on 2026-09-08. Lightweight checks passed.
**No N16/N32 discrepancy modes or 16 continuation fits have been run.** Parent
owns the hypothesis, interpretation and full-run authorization.

## Artifacts

- `discrepancy_weighting.py`: standalone registered experiment importing the
  frozen matched-frequency helpers without changing their files or globals.
- `checks.json`: whitening, equal-trace, block ordering, full residual-Jacobian,
  immutable frozen-input and inherited checkpoint-recovery checks, with hashes.

All 12 prior endpoints are present, converged, and evaluated successfully. The
new runner directly reads `a3_matched_frequency/data_8101.npz`, `data_8102.npz`
and `fits.json`; it does not regenerate observations, reference observations,
held-out means, or pilots. It verifies low-band/reference hashes against all
three prior choices and both reference conditions. Arrays are marked read-only.

## Construction and objective

For each seed/reference condition, use the corresponding stored low-only
endpoint as the common pilot. Evaluate the explicit k18 inverse model at that
pilot on N16 and N32, including its electronics. Define d=fine minus coarse,
delta=`sqrt(2)/sigma * [Re(d); Im(d)]` and its quadrature partner from i*d.
Construction accepts only the pilot, sigma and inverse-grid specification;
neither ADDA truth nor high-frequency observation residuals are arguments.

Four prescribed methods use the same pilot and low+high observations:

| Method | High-block covariance surrogate |
|---|---|
| warm_raw | I |
| isotropic | (1 + delta^T delta / n) I |
| rank1 | I + delta delta^T |
| rank2 | I + E E^T, E=[delta, R(i*d)*sqrt(2)/sigma]/sqrt(2) |

The low-rank symmetric inverse square root uses an SVD of E, retaining the
identity on its orthogonal complement. Isotropic weighting uses its exact
scalar inverse square root. Weights are constructed once from the saved mode
and remain fixed during every continuation. Strength is one. Rank1 versus
isotropic remains the primary registered comparison; rank2 is not selected
after seeing results.

The inherited real residual order is all real training blocks, then all
imaginary training blocks, then the electronics reference. The high-block
indices therefore concatenate `[3*b:4*b]` and `[7*b:8*b]`, where b=144 complex
channels. The same W multiplies those residual rows and all 13 Jacobian
columns. Low and reference rows remain bit-identical to the original objective.

## Bounded verification

Checked all four W matrices against `W C W^T = I`: maximum normalized Frobenius
error **1.61e-15**. The three discrepancy methods have equal trace inflation
within floating-point precision. Symmetry, quadrature-mode orthogonality/equal
norm, and zero-discrepancy identity behavior pass.

Tagged complex data verify high-block real/imag ordering independently of
physics. All 13 columns of the full fixed-weight residual were checked by
centered finite differences for each of four methods with reference absent and
present (eight checks, N8 only). Maximum relative column error **5.29e-9**, below
`1e-7`. Low/reference residuals and Jacobians remained exactly unchanged. These
tests use a fixed synthetic discrepancy vector and the actual frozen matched
observations; they do not construct or tune an N16/N32 discrepancy mode.

Input/source hashes before and after the checks match. No frozen source,
matched estimate, NPZ observation, or matched result artifact is written.

## Fitting, recovery and costs

The prepared experiment has 2 seeds × 2 reference conditions × 4 methods,
N32 continuation, original bounds/scales, `max_nfev=35`, `ftol=xtol=1e-9`,
`gtol=1e-7`. Each continuation starts from the same stored low-only pilot for
its condition. No reference data are added to reference-absent fits.

Each endpoint is atomically saved before parameter or field evaluation.
Frozen low-band evaluation reuses the common k3/k6/k9 held-out receiver means
and the existing sensor/structural metric conventions. High k18 measurement
fit is separately evaluated and reports raw relative field/phase error,
white-noise residual norm/objective, and surrogate-weighted norm/objective.
High measurement fit never enters the common low-band prediction statistic.

The run lock, atomic flush/replace JSON writes, full source/protocol/input hash
manifest, saved mode-vector hashes and pilot hash guard resumability. Completed
fit/evaluation/mode-construction failures are retained and not selectively
retried. Interrupted pre-estimate attempts restart from the pilot with prior
metadata preserved. Frozen estimates resume evaluation without refitting;
completed low/high audit sub-stages are saved and skipped upon evaluation
recovery. An interrupted mode construction restarts with its history retained.
Uncheckpointed interrupted work cost is explicitly unknown, not zero.

Records separate failure stages (mode construction/cache validation,
continuation setup/solve, parameter evaluation, low-band evaluation, and
high-measurement evaluation). Available work counters and attempt times remain
recorded on failures. A failed construction records elapsed time but cannot
recover unreturned internal physical-work counters; no zero-cost claim is made.

The mode JSON preserves N16/N32 per-frequency construction ledgers and setup
cost. Construction is physically shared across the four methods, with a
shared-cost flag; the warm_raw arm does not apply that mode. Pilot setup and
solve costs are explicitly historical. New continuation setup, solve and
attempt wall times are separate from low/high frozen evaluation costs. No sum
is labeled a fresh end-to-end wall measurement and no speed claim is made.

`comparisons.json` will retain signed sensor/structural field and phase,
material and pose differences against low-only, matched low+repeat, matched
raw low+high and each other continuation. For warm_raw against matched raw,
it also reports scaled endpoint distance and raw objective difference for
inspection of initialization dependence. Failed comparisons are explicit.

## Commands

The completed lightweight check command, from the repository root:

```sh
research/trispace_self_calibration/a3_research/.venv3d/bin/python research/delegated/a3_discrepancy_weighting/discrepancy_weighting.py
```

Prepared full command, **not executed; pending parent authorization**:

```sh
research/trispace_self_calibration/a3_research/.venv3d/bin/python research/delegated/a3_discrepancy_weighting/discrepancy_weighting.py --run
```

The script sets single-thread CPU limits before NumPy/SciPy import, uses the
existing environment, and invokes no external models/APIs, installations or
services. The checks establish bounded implementation correctness. Actual
continuation convergence, proxy quality, directional tradeoffs, and scientific
merit remain untested by this implementation step.
