## Name

a3_calibration_v3_plugin_bound_repair_downstream_ablation

## Title

Repairing finite-sample A3 conditional calibration after plug-in bound invalidation, with explicit chart/subspace ablation

## Short Hypothesis

The deterministic A3 identities remain sound, but the prior plug-in B3 branch bound is not calibration-safe under reuse. A corrected cross-validated/empirical bound can provide honest finite-sample calibration, and downstream chart/subspace tradeoffs will be measured as neutral or negative unless regularized alternatives clearly help in tested regimes.

## Related Work

SOM means subspace-based optimization, not self-organizing maps. Generic projection, gain closure, dual residual methods, bootstrap calibration, and linearized/nonlinear least-squares fixtures are prior art. Parent owns literature and scientific judgment; this revision does not claim novelty.

## Abstract

This revision treats the invalid plug-in B3 reused-validation bound as the central finding from the completed round. It keeps the frozen deterministic B1-B4 checks and corrected near-singular identity gate, then diagnoses why the plug-in bound under reuse was invalid in 99.985% of 100,000 repetitions. The revised plan replaces that bound with a corrected leave-one-out/holdout quantile or permutation bound and evaluates it against bootstrap intervals across matched, heavy-tailed, and mismatched-sigma configurations. Downstream work is reframed away from endpoint-accuracy gain: A3 dual-residual correction equals full least squares by identity, so the study instead quantifies chart/subspace error, fixed-chart versus moving-chart behavior, and no-dual-residual ablation across linearized and one small nonlinear 2D pose-graph benchmark. The goal is measured calibration and scope, not a new SLAM or endpoint-gain claim.

## Experiments

### datasets
['Frozen deterministic B1-B4 checks and checks.json from completed round', 'Reused 50-seed Gaussian and heavy-tailed residual fixtures, 2000 draws per config', 'Sharpened B3 reused-validation diagnostic fixtures over 100,000 repetitions', 'Linearized 2D SE(2) pose-graph and 3D pinhole bundle adjustment fixtures from completed round', 'New small nonlinear 2D SE(2) pose-graph trajectories with synthetic odometry and loop closures', 'Condition-number sweep ensemble for near-singular systems']

### baselines
['Plug-in B3 bound vs corrected leave-one-out/holdout quantile bound vs seed-bootstrap percentile', 'Fixed-chart vs moving-chart formula', 'Full A3 dual-residual correction vs no-dual-residual/chart-only correction', 'A3 conditional subspace projection vs unregularized least squares vs ridge-regularized least squares', 'Exact analytic/identity references where available']

### metrics
['Deterministic identity errors and expected-fail control statuses', 'B3 bound validity rate, average invalid bound, miscoverage rate, and interval width for plug-in vs corrected bounds', 'Downstream RMSE mean and standard deviation across Monte Carlo seeds; Wilson failure-rate intervals', 'Chart/subspace error decomposition: correction norm, chart-only error, and full-LS residual relative quantities', 'Nonlinear pose-graph iteration counts, final trajectory/rotation error, and correction identity per linear solve', 'Near-singular residual vs output-error comparisons and conditional correction identities']

### compute_estimate
Single CPU thread, pinned venv, small NumPy/SciPy arrays; frozen checks reused; 50 seeds x 2000 draws, 100k diagnostic, nonlinear pose-graph with 20-50 Monte Carlo trials; roughly 4-8 CPU hours.

### steps
['Re-audit frozen deterministic B1-B4 and prior calibration run as a stable baseline; retain corrected near-singular identity criterion.', 'Diagnose the B3 reused-validation failure: compute empirical distribution of plug-in bound error under reuse and build a corrected leave-one-out/holdout quantile or permutation bound; test validity on independent seeds.', 'Run B3 Gaussian, heavy-tail, and mismatched-sigma configurations with plug-in, corrected, and bootstrap bounds; report validity, width, and average invalid bound.', 'Run explicit fixed-chart vs moving-chart and full vs no-dual-residual ablations on linearized pose-graph/BA fixtures with 50 Monte Carlo seeds and condition sweeps.', 'Implement and run a small nonlinear 2D SE(2) pose-graph benchmark with synthetic noise; compare chart-only, no-dual-residual, full A3 correction, LS, and ridge per linear solve and final trajectory error.', 'Quantify chart/subspace error regimes and downstream tradeoffs against LS/ridge; report negative/neutral outcome where applicable.', 'Audit all non-bool fail flags and expected-fail controls, then produce a calibration/downstream summary with confidence intervals.']

## Risk Factors And Limitations

- The corrected bound may still undercover under heavy tails or model mismatch; claims will be scoped to tested regimes.
- The nonlinear benchmark is small and synthetic, not real SLAM; it only demonstrates behavior under controlled nonlinearities.
- Downstream endpoint improvement remains negative/neutral by identity for linearized fixtures; this revision does not seek or claim practical gain.
- Bootstrap/holdout calibration is finite-sample approximate, not formal proof; the repaired bound may require distributional assumptions.
- Diagnosing the reused-validation invalid bound may be analytically hard; if not fully explained, the corrected empirical bound is presented as measurement rather than theory.
- The prior plug-in B3 bound is not re-introduced as calibration evidence; its invalidation is treated as a core finding.

