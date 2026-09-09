# a2_highdim exploratory summary

- Generated: 2026-09-06T05:15:49+00:00
- Frozen settings sha256: efa530c3fd1a3ee9d1fba9cae98a4cfd24bf1e171092bef7521a712d0d02dab7

## Scope and epistemic status

EXPLORATORY 2D imaging feasibility only; not the A2 E4 primary endpoint and not publication or production acceptance. Reconstructions are model-conditional within the 49-Gaussian basis and the scalar full-wave physics core. The known-true-pose method is an explicit oracle upper reference and is not a practical baseline.

## Run counts and statuses

| phase | method | n | statuses |
|---|---|---|---|
| tune | coherent_fixedpose | 2 | complete=2, budget_stopped=0, task_error=0, time_stopped=0 |
| tune | coherent_joint | 2 | complete=2, budget_stopped=0, task_error=0, time_stopped=0 |
| tune | intensity_joint | 2 | complete=2, budget_stopped=0, task_error=0, time_stopped=0 |
| tune | oracle | 2 | complete=2, budget_stopped=0, task_error=0, time_stopped=0 |
| test | coherent_fixedpose | 12 | complete=12, budget_stopped=0, task_error=0, time_stopped=0 |
| test | coherent_joint | 12 | complete=12, budget_stopped=0, task_error=0, time_stopped=0 |
| test | intensity_joint | 12 | complete=12, budget_stopped=0, task_error=0, time_stopped=0 |
| test | oracle | 12 | complete=12, budget_stopped=0, task_error=0, time_stopped=0 |
| mismatch | coherent_fixedpose | 3 | complete=3, budget_stopped=0, task_error=0, time_stopped=0 |
| mismatch | coherent_joint | 3 | complete=3, budget_stopped=0, task_error=0, time_stopped=0 |
| smoke | coherent_fixedpose | 1 | complete=1, budget_stopped=0, task_error=0, time_stopped=0 |
| smoke | coherent_joint | 1 | complete=1, budget_stopped=0, task_error=0, time_stopped=0 |
| smoke | intensity_joint | 1 | complete=1, budget_stopped=0, task_error=0, time_stopped=0 |

## Unconditional test metrics (seeds 801-806, radii .125/.5 lambda, n=12 per method)

### Spatial map RMSE

| method | n | mean | median | std | min | max |
|---|---|---|---|---|---|---|
| coherent_fixedpose | 12 | 0.1054 | 0.1010 | 0.0545 | 0.0438 | 0.1892 |
| coherent_joint | 12 | 0.0261 | 0.0264 | 0.0029 | 0.0209 | 0.0304 |
| intensity_joint | 12 | 0.0430 | 0.0425 | 0.0057 | 0.0338 | 0.0502 |
| oracle | 12 | 0.0251 | 0.0257 | 0.0029 | 0.0206 | 0.0285 |

### Coefficient RMSE (49 real coefficients)

| method | n | mean | median | std | min | max |
|---|---|---|---|---|---|---|
| coherent_fixedpose | 12 | 0.1255 | 0.1189 | 0.0390 | 0.0712 | 0.1916 |
| coherent_joint | 12 | 0.0682 | 0.0646 | 0.0123 | 0.0529 | 0.0859 |
| intensity_joint | 12 | 0.0817 | 0.0838 | 0.0155 | 0.0613 | 0.1036 |
| oracle | 12 | 0.0666 | 0.0678 | 0.0118 | 0.0485 | 0.0804 |

### Pose lever-arm error (m)

| method | n | mean | median | std | min | max |
|---|---|---|---|---|---|---|
| coherent_fixedpose | 12 | 0.1562 | 0.1562 | 0.0979 | 0.0625 | 0.2500 |
| coherent_joint | 12 | 0.0124 | 0.0129 | 0.0048 | 0.0058 | 0.0212 |
| intensity_joint | 12 | 0.0262 | 0.0181 | 0.0179 | 0.0097 | 0.0587 |
| oracle | 12 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |

### Masked forward-phase residual (vs clean truth total; |truth|>3 sigma, evaluation-only)

| method | samples | mean abs (rad) | median abs (rad) | rms (rad) | circular std (rad) |
|---|---|---|---|---|---|
| coherent_fixedpose | 3456 | 0.0284 | 0.0106 | 0.0797 | 0.0753 |
| coherent_joint | 3456 | 0.0067 | 0.0032 | 0.0117 | 0.0117 |
| intensity_joint | 3456 | 0.0114 | 0.0056 | 0.0204 | 0.0204 |
| oracle | 3456 | 0.0063 | 0.0031 | 0.0109 | 0.0109 |

## Paired differences, 6-seed cluster bootstrap (descriptive 95% CIs)

| metric | comparison | mean difference | 95% CI |
|---|---|---|---|
| spatial map RMSE | coherent_joint minus coherent_fixedpose | -0.0793 | [-0.0897, -0.0691] |
| coefficient RMSE | coherent_joint minus coherent_fixedpose | -0.0573 | [-0.0700, -0.0454] |
| pose lever-arm error | coherent_joint minus coherent_fixedpose | -0.1439 | [-0.1474, -0.1403] |
| spatial map RMSE | intensity_joint minus coherent_fixedpose | -0.0624 | [-0.0721, -0.0524] |
| coefficient RMSE | intensity_joint minus coherent_fixedpose | -0.0438 | [-0.0517, -0.0338] |
| pose lever-arm error | intensity_joint minus coherent_fixedpose | -0.1301 | [-0.1425, -0.1158] |
| spatial map RMSE | intensity_joint minus coherent_joint | 0.0169 | [0.0121, 0.0212] |
| coefficient RMSE | intensity_joint minus coherent_joint | 0.0135 | [0.0042, 0.0222] |
| pose lever-arm error | intensity_joint minus coherent_joint | 0.0138 | [0.0021, 0.0273] |
| spatial map RMSE | oracle minus coherent_joint | -0.0010 | [-0.0020, -0.0001] |
| coefficient RMSE | oracle minus coherent_joint | -0.0017 | [-0.0050, 0.0017] |
| pose lever-arm error | oracle minus coherent_joint | -0.0124 | [-0.0159, -0.0088] |

Bootstrap convention: resample the six seed clusters (both pose-error radii travel together), 10,000 draws; percentile 2.5-97.5. Exploratory descriptive intervals, not a confirmatory gate.

## Mismatch controls (seed 881)

| kind | method | spatial map RMSE | pose err (m) | fit chi2 p | val chi2 p | false-assurance flag |
|---|---|---|---|---|---|---|
| clock_phase | coherent_fixedpose | 0.1292 | 0.0625 | 0.000e+00 | 0.000e+00 | False |
| clock_phase | coherent_joint | 0.1373 | 0.0800 | 0.000e+00 | 0.000e+00 | False |
| receiver_coupling | coherent_fixedpose | 0.0460 | 0.0625 | 0.000e+00 | 0.000e+00 | False |
| receiver_coupling | coherent_joint | 0.0396 | 0.1559 | 0.000e+00 | 0.000e+00 | False |
| outside_basis | coherent_fixedpose | 0.0405 | 0.0625 | 4.507e-01 | 5.999e-03 | True |
| outside_basis | coherent_joint | 0.0390 | 0.0859 | 9.988e-01 | 7.043e-02 | True |

These controlled mismatches test the coherent geometry-only model against an unmodelled clock phase, receiver coupling and an outside-49-basis inclusion. A fitted residual that looks calibrated under model mismatch is flagged as false-assurance potential; no geometry proof and no calibration guarantee is claimed. The fixed-pose result merely fixes geometry and does not remove mismatch; target artifacts must not be read as validated localisation.

## Figures

- `figures/reconstruction_map.png`
- `figures/error_distributions.png`
- `figures/phase_residuals.png`
- `figures/mismatch_residuals.png`

## Work/wall accounting

Per-run RHS solve columns (forward states, adjoints and any initialisation) are reported in work.units_rhs_columns with kind counts; LU factorisations, operator products and wall seconds are separate raw ledger fields and are never converted into solve-equivalent units. Evaluation-only N32/N20 forward solves are in evaluation_ledger and never feed estimator decisions. This report draws no work-advantage conclusion.

| phase | n runs | median solver wall (s) | total solver wall (s) | median eval wall (s) | total eval wall (s) |
|---|---|---|---|---|---|
| tune | 8 | 4.61 | 37.3 | 0.84 | 6.9 |
| test | 48 | 4.59 | 223.3 | 0.85 | 41.1 |
| mismatch | 6 | 4.89 | 29.4 | 0.92 | 5.6 |
| smoke | 3 | 0.95 | 3.0 | 0.84 | 2.5 |

## Exact commands

```bash
VENV=experiments/idea_loops/loop_2026-09-04_02-58-31/experiment_geometry_lifted_trispace_som/.venv/bin/python
# bounded implementation checks
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 $VENV research/delegated/a2_highdim/test_highdim.py
# basic one-iteration pipeline smoke (seed 71, pre-freeze settings)
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 $VENV research/delegated/a2_highdim/runner.py --mode smoke
# tuning seeds 71-72, then frozen.json written
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 $VENV research/delegated/a2_highdim/runner.py --mode tune
# 48 final exploratory runs, seeds 801-806, frozen settings only
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 $VENV research/delegated/a2_highdim/runner.py --mode test
# 6 controlled mismatch runs, seed 881
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 $VENV research/delegated/a2_highdim/runner.py --mode mismatch
# aggregate summary and figures (reproducible, no reruns)
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 $VENV research/delegated/a2_highdim/runner.py --mode summarize
```

## Limitations

- Exploratory single aperture (full), one noise level (30 dB) and six test scenes with two pose-error radii; no final E4 seeds 1001-1020 were used.
- The inverse N20 grid required a runtime-only extension of the physics Config allowed-N tuple; physics.py itself was not changed and no formula was altered.
- The 49-Gaussian coefficient basis is strongly overlapping (Gram condition ~2.7e4 at N20), so coefficient RMSE is not an interpretable spatial measure; spatial N32 map RMSE is the imaging metric.
- Cluster-bootstrap intervals are descriptive only and cannot establish equivalence or superiority.
- The matched-intensity comparison uses the induced Rice magnitude NLL and exactly the parent |y|; it is not a direct-intensity sensor model.
- Fixed-pose and joint methods share the frozen schedule and identical settings; oracle shares them too but uses x_true, so its advantage is not a calibration claim.
- Mismatch controls are misspecification tests of the coherent geometry-only model; an augmented nuisance estimator was not implemented.