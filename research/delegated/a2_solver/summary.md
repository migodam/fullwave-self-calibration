# A2 E4 solver implementation summary

Implemented under `research/delegated/a2_solver`:

- `solver.py`: `solve(model,y,sigma,alpha0,x0,method,budget,settings)` with
  methods direct, prasc, fixed_rank, phaseless, direct_control; shared
  four-stage cumulative frequency budget, L-BFGS-B optimizer, box
  constraints, scaled pose coordinates, optional identical quadratic prior.
- `reduced.py`: certified SOM-informed reduced physical states/sensitivities
  (sensing-SVD basis, optional sequential TSOM RRQR block, rank ladder,
  inexact physical Jacobian labelled as such).
- `objective.py`: exact coherent Gaussian objective plus matched
  Rice/noncentral-chi-square intensity NLL with stable i0e/i1e and adjoint
  gradient.
- `cert.py`, `reference.py`, `generate.py`, `calibrate.py`, `common.py`:
  passive certificate integration, reference design/margins, seed scenes,
  material-only initialisation, repeatable tuning timing calibration,
  cost ledger and final-seed gate.
- `run_e4.py`: `--mode tune|final|summarize --budget 200|800`;
  final mode requires `--frozen-config`.

Evidence:

- 6 unit/smoke checks pass (see `checks.json` and `test_solver.py`).
- Tuning seeds 1-10 executed at budgets 200 and 800 (one init/seed;
  no final seed access). Canonical records:
  `records_tune_seeds1-10_800d.jsonl`; `tuning_summary.md` is in this
  directory; no tuning candidate passed the A2 joint success gate.
- Reference design metrics for full and limited apertures are
  full-rank and frozen in `reference_design.json`.
- Timing calibration is frozen in `calibration_tune.json` and
  `timing_summary.md`.

Unresolved issues are listed at the end of `tuning_summary.md`; final
generation remains blocked behind an explicit parent `--mode final`
invocation with the frozen config.
