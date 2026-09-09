Done. All artifacts and runs are under `research/delegated/a2_highdim/`.

Facts
- Bounded checks: 9/9 pass (N20 adapter, noise reference, corrected scaled-coordinate gradient grad_z = grad_x/scale via FD, Rice NLL gradient FD, exact work charges, budget-preserves-accepted, mismatch controls, E4-seed guard).
- Frozen settings written before test seeds: sha256 `efa530c3…d02dab7`, timestamped in [frozen.json](/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/a2_highdim/frozen.json).
- 48 test runs (seeds 801–806 × radii .125λ/.5λ × 4 methods), 12 per method, all stages visited, no task errors. Spatial map RMSE mean (N32 grid): fixed-pose 0.105, coherent joint 0.026, matched-intensity joint 0.043, oracle (explicitly labelled) 0.025. Coherent-joint pose lever error mean 0.012 m vs 0.026 for intensity joint and 0.156 for fixed pose; coherent joint tracks the oracle closely in this exploratory regime.
- 6 controlled mismatch runs (seed 881): clock phase and receiver coupling leave badly inconsistent residuals (fit chi-square p ≈ 0) for both coherent methods; the outside-49-basis scene gives coherent joint a near-perfect in-sample fit while pose drifts to 0.086 m — flagged as false-assurance potential, not localisation evidence.
- Records/work/wall: [records.jsonl](/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/a2_highdim/records.jsonl) (65 records incl. smoke/tune/test/mismatch), per-run RHS columns ≈ 2,400–2,600 of the 5,000 cap with factorisations/operator products/wall reported separately; test phase ~224 s solver + ~41 s evaluation wall.

Key paths
- [summary.md](/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/a2_highdim/summary.md) and [summary.json](/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/a2_highdim/summary.json)
- [runner.py](/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/a2_highdim/runner.py), [test_highdim.py](/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/a2_highdim/test_highdim.py), [checks.json](/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/a2_highdim/checks.json)
- Figures: [reconstruction_map.png](/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/a2_highdim/figures/reconstruction_map.png), [error_distributions.png](/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/a2_highdim/figures/error_distributions.png), [phase_residuals.png](/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/a2_highdim/figures/phase_residuals.png), [mismatch_residuals.png](/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/a2_highdim/figures/mismatch_residuals.png)

Uncertainties / boundaries
- Exploratory model-conditional feasibility only: scalar full-wave 2D, one aperture/noise level, 49 overlapping Gaussian coefficients, and no E4 seeds 1001–1020.
- N=20 required a runtime-only extension of the physics `Config` allowed-N tuple; `physics.py` was not modified.
- Bootstrap intervals are descriptive (6 seed clusters), not confirmatory; mismatch controls test coherent geometry-only misspecification, with no nuisance estimator or calibration guarantee.