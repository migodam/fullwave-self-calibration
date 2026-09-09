# Family 10 JSON inspection (read-only)

Inspection of `results/family10_online_slam_toy.json` on 2026-09-03.
No modification to the result file; no computation rerun. Read with
`.venv/bin/python` + stdlib `json`.

## Provenance / timing

- `generated_utc`: `2026-09-03T13:00:42.742177+00:00`
- `runner`: `src/family10_online_slam_toy.py`
- `command` recorded in JSON: `.venv/bin/python src/family10_online_slam_toy.py`
- `runtime_seconds`: `41.010709166002925`
- Source SHA-256 fields all match the current source files on disk.

## Config highlights

- `N = 16`, `T = 6`, `n_rx = 4`, `m_c = 24`, `m_real = 144`,
  `p = 24` (coefficient dim), `q_pose = 18` (3 pose params x 6 poses)
- `frequencies = [1.0, 1.4, 1.8]`
- `snr = 100.0`, `alpha = 1.0`
- `monte_carlo.n_trials = 3`, `seed = 20260903`,
  `modes = ["born", "full_wave"]`
- `n_confounded_directions = 3`
- Solver: `scipy least_squares`, `method = "trf"`, `x_scale = "jac"`,
  `max_nfev = 1200`, `xtol = ftol = gtol = 1e-10`
- Fit note: local fits start at the true parameters inside a Fisher-whitened
  coordinate system; `success` means scipy status True.

## Per-mode Monte Carlo (4 effective modes)

Each family (`born`, `full_wave`) stores `n_trials_requested = 3` and two pose
modes, `known` and `free`. All four modes: 3 attempted, 3/3 successful.

| family | pose mode | nfev (min/med/max) | cost (min/med/max) | seconds (min/med/max) |
| --- | --- | --- | --- | --- |
| born | known | 4 / 4 / 4 | 50.94 / 53.44 / 64.74 | 0.066 / 0.066 / 0.067 |
| born | free | 134 / 141 / 199 | 36.55 / 39.96 / 45.63 | 2.27 / 2.36 / 3.37 |
| full_wave | known | 22 / 32 / 32 | 49.19 / 54.16 / 64.21 | 1.29 / 1.84 / 1.84 |
| full_wave | free | 132 / 140 / 198 | 37.69 / 40.95 / 44.73 | 7.59 / 8.05 / 11.37 |

Per-trial records contain `success`, `status`, `nfev`, `cost`, `message`,
`c_error_l2`, `seconds` (plus `dx_l2` for free-pose trials). No
`optimality` field is stored.

## Covariance comparison / trace fields (exact values)

Born:
- `known_rel_fro = 0.8992009874573181`
- `free_rel_fro = 42.23595506090035`
- `known_rel_spectral = 0.8526273023516002`
- `free_rel_spectral = 42.443298138227135`
- `known_trace_emp = 15955.083534668953`
- `known_trace_pred = 28930.450038682666`
- `free_trace_emp = 6779811.499136179`
- `free_trace_pred = 175559.33086239817`
- `trace_ratio_empirical = 424.93111893800256`
- `trace_ratio_predicted = 6.0683235355018414`

Full wave:
- `known_rel_fro = 0.8266145025997169`
- `free_rel_fro = 5.919801677377825`
- `known_rel_spectral = 0.7995421871875931`
- `free_rel_spectral = 6.619427076864533`
- `known_trace_emp = 5198.199429597969`
- `known_trace_pred = 6444.67990263276`
- `free_trace_emp = 110181.28244770705`
- `free_trace_pred = 23388.061528127808`
- `trace_ratio_empirical = 21.196047581465823`
- `trace_ratio_predicted = 3.62904936808008`

## Finite-difference self-checks

- Exactly one self-check is stored:
  `finite_difference_self_check = {mode: "born", eps: 1e-6,
  direction_l2: 1.0, rel_fd_error: 7.698814633667721e-10,
  residual_dim: 162, jacobian_shape: [162, 42]}`.
  Config note says it checks the realified stacked residual Jacobian at the
  true point (Born `B`; duplicated locally because no Born pose Jacobian API is
  reused).
- There is NO stored finite-difference self-check field for the full-wave
  A/B stacks in this JSON. `theory.full_wave` records A/B/y shapes
  (`144 x 24`, `144 x 18`, `144`) but no FD error fields. A recursive scan for
  `fd`/`finite`/`optimality`/`grad` keys inside `theory` and `monte_carlo`
  returned none.

## Caveats

- `config.monte_carlo.n_trials = 3`, but the current source default in
  `src/family10_online_slam_toy.py` is 200 and the stored source hashes match.
  The recorded `command` carries no `--n-trials` flag, so the recorded command
  string appears not to capture the actual invocation flags; the data itself is
  consistent with a 3-trial run.
- "Four modes" in this file are the combinations of two forward models
  (`born`, `full_wave`) x two pose treatments (`known`, `free`), not four
  distinct top-level mode labels.
- Inspected file: `results/family10_online_slam_toy.json`.
