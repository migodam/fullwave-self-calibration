# A3 physical-ROM reconstruction — development summary

Status: numerical foundation implemented, algebra checks pass, rank-error
audit and 56 development nonlinear runs completed.  Development only (seeds
4101–4104); no final seeds and no SOM-falsification claim.

## Artifacts

| File | Purpose |
|---|---|
| `common.py` | fixed grid16/q=9/k=3,6,12 config, scenes, transparent WorkLedger |
| `wave.py` | read-only physical core reuse; fixed-chart LS proxy + residual-term derivative (QR triangular solves) |
| `basis.py` | sensing SVD, sequential Twofold, generic/task/primal-adjoint enrichment, block Krylov charts |
| `audit.py` | exact state/derivative offline audits |
| `solver.py` | direct analytic GN, exact-adjoint baseline, ROM hybrids with shared exact acceptance controller |
| `test_rom.py` | independent FD/adjoint algebra checks |
| `rank_study.py` / `rank_study_raw.json` | fixed-chart rank-error study |
| `run_nonlinear.py` / `records_nonlinear.jsonl` | 56 development nonlinear runs |
| `checks.json`, `environment.json` | check evidence and runtime environment |

## Key numerical evidence

Fixed-chart derivative tests (grid8 smoke, rank 48 chart):

- Fixed-chart coefficient derivative with residual term vs centred FD:
  relative error `4.1e-11`; deliberately omitted-residual version `4.3e-4`.
- Fixed-chart output Jacobian vs FD: `7.7e-10`.
- Reduced loss gradient vs FD: `1.8e-10`.

Rank audit (grid16, three samples, thresholds state/output `1e-3`,
map/pose Jacobian `1e-2`), first-pass target ranks over samples:

| Method | k=3 | k=6 | k=12 |
|---|---:|---:|---:|
| sensing-only | fail | fail | fail |
| Twofold | 96/96/96 | 128/128/128 | 160/192/160 |
| generic task enrichment | 32/48/32 | 64/64/64 | 96/128/96 |
| sensing + task | 96/96/96 | 96/96/96 | 96/128/96 |
| Twofold + task | 96/96/96 | 96/96/96 | 128/128/128 |
| block Krylov | 96/192/96 | 160/256/160 | 256/256/256 |

Development nonlinear runs (4 seeds x methods x continuation schedule):

- ROM methods accepted reduced steps in the all-frequency schedule: on
  average 5–7 per run for generic/twofold-task and 6.8 for block Krylov,
  with every accepted step in the stage that includes k=12.  Generic
  continuation additionally averaged ~6 accepted reduced steps in the
  k=3,6,12 stage across the four seeds.
- ROM runs used fewer full-wave RHS columns than direct analytic GN but more
  wall time (chart construction is fully charged); no speedup claim is made.
- Continuation vs no-continuation on the same final data was not an
  advantage here (mean pose error 0.040 m vs 0.031 m for generic task over
  four seeds).  This is a small development observation, not a claim.
- Exact acceptance controller stops (`no_improving_trial_and_fallback`) are
  preserved in `records_nonlinear.jsonl`; low-loss terminations are not
  relabelled as success gates.

## Correction applied and preserved failures

An initial enrichment issue (near-duplicate directions made appended charts
non-orthonormal) was found by the algebra tests and fixed with two-pass
reorthogonalisation of chart columns.  The fixed-chart residual-term formula
was implemented before the nonlinear runs; its omission is checked in
`test_rom.py` to make the failure mode explicit.  Every nonlinear status,
fallback count, accepted/rejected trial and raw cost ledger is preserved.

## Uncertainties / incomplete comparators

- The chart budgets in `solver.CHART_RANKS` come from the development rank
  audit; they are not independent hold-out tuning.
- The material-weighted domain variant and exact-intersection Twofold were
  not implemented; the Twofold seed is the sequential projection
  `[Vs,(I-Ps)Vd]` as required.
- The reduced state is the physical dependent current on the chart
  (``r_free=0``); there is no separately fitted free-current nuisance.
- Exact state/derivatives are offline audits only; chart construction does
  not call full `M^-1` state solves.
- These checks validate implementation and show the reduced path executes
  with accepted high-frequency steps; they are not physical-model
  verification, a performance comparison, or SOM falsification.

## Exact reproduction commands

```bash
PY="experiments/idea_loops/loop_2026-09-04_02-58-31/experiment_geometry_lifted_trispace_som/.venv/bin/python"
"$PY" research/delegated/a3_rom/test_rom.py
"$PY" research/delegated/a3_rom/rank_study.py
"$PY" research/delegated/a3_rom/run_nonlinear.py
```

Single-thread CPU NumPy/SciPy only, no installs, no GPU, no child workers.
