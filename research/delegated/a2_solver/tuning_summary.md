# A2 E4 tuning summary

- records: 200

| method | runs | success | pose RMSE m | task map RMSE | full map RMSE | wall s | failures |
|---|---:|---:|---:|---:|---:|---:|---:|
| direct | 40 | 0.000 | 0.06508 | 0.1101 | 0.1413 | 0.685 | 40 |
| prasc | 40 | 0.000 | 0.2687 | 0.3421 | 0.3594 | 0.44 | 40 |
| fixed_rank | 40 | 0.000 | 0.4545 | 0.3537 | 0.3725 | 0.385 | 40 |
| phaseless | 40 | 0.000 | 0.05551 | 0.1268 | 0.1566 | 0.717 | 40 |
| direct_control | 40 | 0.000 | 0.08983 | 0.1546 | 0.189 | 0.687 | 40 |

## Bounded tuning observations (worker record, not a scientific verdict)

- Ran tuning seeds 1-10 only, one representative initialization per seed
  (`--inits by_seed`) at budget 800, four equal optimizer candidates per
  method, fixed ranks 16/36/64 cycled for `fixed_rank`. The full 12-init
  grid is reserved for explicit `--mode final`.
- No tuning candidate met the joint success gate; every run's status was
  budget exhaustion under the frozen stage schedule (see `status` and
  `failure` in each record). Direct, phaseless, and direct_control moved
  off the initial pose; fixed_rank moved in 14/40 runs; prasc did not take
  a reduced step before its conservative promotion chain exhausted the cap.
- prasc certificates were available at all examined ranks but the
  conservative passive-medium residual/gradient gates never admitted a
  reduced chart (promotion fallback events are recorded per run).
- All runs' shared material-only initialisation hit the declared
  `init_max_units` cap (100% `init_fallback`) and used the predeclared
  fixed nonzero alpha=0.5; recorded in `init_fallback`/`init_reason`.
- `covariance_surrogate_pass` was true at the final exact Jacobian in all
  200 records while validation success was false; the surrogate-only label
  never served as a coverage certificate. `false_calibrated` records this
  surrogate-vs-validation gap, not an evaluator calibration claim.
- Frozen config: `frozen/frozen_seeds1-10_800d.json` (explicit path
  required by final mode).

## Open items for the parent

1. Under this implementation's stage budget fractions and L-BFGS-B
   combined f/g line search, no method completed the full four-stage
   schedule inside the 200-unit cap, and every budget-800 candidate ended
   in budget exhaustion. Parent should decide whether cap/fraction
   semantics, f/g-per-line-search cost, or budget scale needs revision
   before final generation.
2. Conservative certificates force prasc toward direct fallback at every
   current tuning point; scientific interpretation stays with the parent.
3. Frozen optimizer settings were selected by tuning-loss tie-break because
   all success rates were zero.
