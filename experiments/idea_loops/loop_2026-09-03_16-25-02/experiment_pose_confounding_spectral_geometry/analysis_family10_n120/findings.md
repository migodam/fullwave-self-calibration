# Family10 n120 stored-JSON diagnostics

Source (read-only): `results/family10_online_slam_toy_n120.json`
Scripts: `analysis_family10_n120/diagnostics.py`, `analysis_family10_n120/plots.py`
Full console capture: `analysis_family10_n120/diagnostics_summary.txt`

## Structure

Per trial (both `born` and `full_wave`, 120 completed trials each):

- `known`: `success`, `status`, `nfev`, `cost`, `message`, `c_error_l2`,
  `seconds`
- `free`: the same fields plus `dx_l2`

`c_error_l2 = ||c_hat - c0||_2` (map error norm).  `dx_l2 = ||dx_hat||_2`
(pose MAP displacement; the pose truth is `dx=0`, so this is the free-pose
pose-error norm).  Per-trial error **vectors** (`errs_known`, `errs_free`,
`ok_*`) and `dx_hat` vectors are computed by the experiment but stripped before
JSON serialization (see `src/family10_online_slam_toy.py` around line 1523).
Only norms survive.  Known-pose fits keep `dx=0`, so their pose error is
identically zero and is not stored.

## 1. Per-mode summary (successful trials)

| mode | n | n success | min | median | mean | max of ||e||^2 |
|---|---:|---:|---:|---:|---:|---:|
| born_known | 120 | 120 | 969 | 1.696e4 | 3.048e4 | 2.387e5 |
| born_free | 120 | 115 | 6.247e3 | 2.378e5 | 2.569e6 | 5.549e7 |
| full_wave_known | 120 | 120 | 534.8 | 6.586e3 | 9.296e3 | 4.439e4 |
| full_wave_free | 120 | 119 | 568.4 | 5.798e4 | 1.025e5 | 7.350e5 |

Pose-error norm `||dx_hat||` (free only):

| mode | min | median | mean | max |
|---|---:|---:|---:|---:|
| born_free | 0.765 | 1.623 | 1.695 | 2.851 |
| full_wave_free | 0.874 | 1.671 | 1.737 | 3.028 |

Failures are all `status=0` / `max_nfev=1200`: born free trials {5, 39, 77, 79,
118}, full-wave free trial {35}; none in known-pose fits.

## 2. Born-free outliers (>20x median ||e_free||^2)

Median `||e_free||^2 = 2.378e5`; threshold `4.757e6`.  14 of 115 successful
born-free trials qualify: {1, 2, 31, 38, 45, 51, 57, 67, 69, 74, 80, 93, 94,
112}, with `||e||^2` from `4.76e6` to `5.55e7`.

- All 14 have `success=True`, `status=2` (`ftol`), `nfev` 77-248, cost
  36.2-57.2; they are converged local solutions, not solver failures.
- None is an outlier in the same trial's known-pose fit.  Their known-pose
  `||e||^2` is only 0.21-4.28x the known median (all below 20x), so the huge
  map error is a **free-pose-specific** effect: freeing pose lets the map
  solution escape into low-information directions.
- Free/known `||e||^2` ratios for these trials are 129-5904.
- Log-log correlation between known and free `||e||^2` over successful born
  trials is ~-0.03 (Spearman), so born free error is essentially unrelated to
  the map-noise level seen by the known-pose fit.
- The 14 outlier trials account for 81.8% of the total born-free sum of squared
  map errors.  The empirical `Cov_free` trace is dominated by one eigen-mode:
  top eigenvalue 2.195e6 vs predicted top `P_free` eigenvalue 1.535e5
  (86.6% of empirical trace in the top mode; top three modes 99.3%).

Full-wave free has **no** trials above 20x median (`||e||^2` max 7.35e5 vs
threshold 1.159e6); it is a broad inflation rather than a fat outlier tail.

## 3. Trace ratio excluding outliers

Exact re-computation on the cleaned subset is impossible from this JSON because
the per-trial error vectors and the cleaned-subset mean vector are not stored.
What is exact is the full-sample identity check
`tr(Cov) = (sum||e||^2 - n||mean||^2)/(n-1)`, which reproduces the stored trace.

Born free, full sample: `tr = 2.536e6`; stored trace ratio 82.92 (predicted
6.07).

Born free, cleaned (101/115 trials): `sum||e||^2 = 5.381e7` (18.2% of full);
using the stored full-sample mean as the only available stand-in for the
cleaned mean gives approximate `tr ~ 4.83e5` and approximate trace ratio
`~15.8`; an uncentered upper proxy is `5.38e5` (~17.6 ratio).  So excluding the
14 outliers removes most of the excess, but a ~2.7-3x gap above the predicted
free trace (`1.756e5`) remains.

Full-wave free has no outliers, so no cleaning changes anything:
`tr = 1.027e5`, trace ratio 11.02 (predicted 3.63).

## 4. Pose MAP wandering

Pose MAP displacement norm is 1.6-1.7 on average and never stays near zero
(min 0.77 born / 0.87 full-wave; max 2.85 / 3.03).  Born outlier trials are not
distinguished by larger `||dx||` (outlier dx range 0.91-2.55), so the map blow
up is not simply "larger dx"; it is landing in map directions that are almost
unobservable once pose is free.

## 5. Born-free absolute scale

| quantity | value |
|---|---:|
| predicted `P_free` trace | 1.756e5 |
| predicted `P_free` eigenvalue range | 9.80e-3 .. 1.535e5 |
| empirical `Cov_free` trace | 2.536e6 |
| empirical `Cov_free` diagonal min/median/mean/max | 129 / 2.865e4 / 1.057e5 / 5.317e5 |
| empirical / predicted trace | 14.45x |

The JSON stores only `P_free` eigenvalues (not its full matrix or diagonal), so
the predicted per-coordinate diagonal is unavailable without re-running the
linearized algebra; mean diagonal would be `tr/p = 7315` if it were uniform.

## Bottom line

The free-pose empirical covariance is much larger than the linearized
prediction mainly because of rare-but-converged born trials whose map estimate
lands in near-unobservable directions once pose is free (82% of the born error
sum of squares in 14/115 trials), plus a remaining broad 2.7-3x inflation and a
full-wave 4.4x inflation even without born-style outliers.  Known-pose traces
are close to prediction (born 1.06x, full-wave 1.45x), and per-trial pose MAP
errors stay bounded (~0.8-3.0 norm) but far from the linearization point.
