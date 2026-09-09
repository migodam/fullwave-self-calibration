Implemented the E4 solver package under `research/delegated/a2_solver` without touching the physics core or spawning child workers.

Key evidence:

- API delivered as requested: `solve(model, y, sigma, alpha0, x0, method, budget, settings)` in [solver.py](/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/a2_solver/solver.py) and `run_e4.py --mode tune|final|summarize --budget 200|800` in [run_e4.py](/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/a2_solver/run_e4.py). Final mode requires an explicit `--frozen-config`; tuning cannot access seeds 1001–1020.
- All five methods implemented: direct, prasc (certified reduced states, conservative fallback), fixed_rank (no adaptive certificate), phaseless (exact induced Rice/noncentral-chi-square NLL with stable i0e/i1e adjoint), and direct_control.
- Six unit/smoke checks pass, including reduced-vs-exact consistency at a full basis, phaseless gradient finite differences, certificate behavior, ledger accounting, and the final-seed gate: [checks.json](/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/a2_solver/checks.json).
- Tuning seeds 1–10 ran at budget 800 (one representative init per seed): 200 records in [records_tune_seeds1-10_800d.jsonl](/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/a2_solver/records_tune_seeds1-10_800d.jsonl). No final seed was generated or inspected.
- Reference design is full-rank for both apertures and frozen with margins in [reference_design.json](/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/a2_solver/reference_design.json); timing calibration is frozen in [calibration_tune.json](/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/a2_solver/calibration_tune.json).

Findings and uncertainties to escalate:

1. No tuning candidate passed the joint success gate. Every run’s status was budget exhaustion; direct/phaseless/direct_control moved off the initial pose at budget 800, fixed_rank moved in 14/40 runs, but prasc never took a reduced step before its promotion chain exhausted the cap.
2. The conservative passive-medium certificates were always available but never admitted a reduced chart; prasc therefore recorded only certificate-driven fallback/promotion events. This is reported rather than silently replaced by true-error checks.
3. Under the specified stage fractions and L-BFGS-B combined f/g line search, no method can even complete its first accepted iteration inside 200 units; even 800 ends in budget exhaustion. The parent should decide whether that cap/fraction semantics or the f/g-per-line-search accounting needs revision before final generation.
4. All shared material-only initializations hit the declared `init_max_units` cap and fell back to the predeclared alpha=0.5, with the design change recorded per run.
5. Frozen optimizer settings were selected by loss tie-break because success rates were zero; the frozen path is [frozen_seeds1-10_800d.json](/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/a2_solver/frozen/frozen_seeds1-10_800d.json).
6. Tuning used one init per seed for CPU control; the full 12-init grid remains reserved for explicit final mode.

Full summaries: [summary.md](/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/a2_solver/summary.md), [tuning_summary.md](/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/a2_solver/tuning_summary.md), and [timing_summary.md](/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/a2_solver/timing_summary.md).