# Replicate and close bounded A2 theory-validation checks

Workdir /Volumes/migodam's-external-brain/Research/Inv_SLAM. Isolated DeepSeek Flash worker. No recursive workers or AI Scientist. Do not change physics or solver folders or global configuration. Use apply_patch for edits; CPU single thread. Return <=10 findings, artifact paths, missing work and uncertainties. All detailed outputs under research/delegated/a2_theory_validation/replication/.

Read original TASK.md here, A2_PRASC_SOM_VALIDATION_PROTOCOL.md, and existing package at experiments/idea_loops/loop_2026-09-05_17-55-24_a2_development/experiment_a2_prasc_theory_validation/. Pipeline generated useful E1/E2/E3/E5 but round2 was incomplete. Do not claim round2 replication existed. Reuse code by importing without modifying original outputs. Python experiments/idea_loops/loop_2026-09-04_02-58-31/experiment_geometry_lifted_trispace_som/.venv/bin/python. OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1.

1. Actually rerun the full existing 50 tests and fixed seeds to new output location. Inventory CLI before invoking; scripts may default to historical paths, avoid overwrite by importing modules or copying provided code with apply_patch if needed.
2. Execute optional physical E1 Fisher and E3 limited aperture scripts, capture sample sizes, all pass/fail and non-vacuity (physical nuisance saturation must be explicit). Do not invent a new theorem or assume MC Loewner ordering exact.
3. Execute threshold/negative-control ablations: physical vs free-current envelope; relative-only vs absolute plus bias gates; shared-map vs independent-map; greedy vs pair-lookahead. Existing modules already have exact examples. Vary gate/noise/rank thresholds across 3 values to expose false assurances; record raw data, no selection.
4. Extend E5 physical acquisition with 8 receiver configurations if readily supported by Model and existing physical runner. Use fixed estimated/reference state only for acquisition decision, separate truth evaluation. If unavailable document not run; do not use oracle pose scores or fake nonlinear evidence.
5. summary.md with exact executed commands, environment, run counts, runtimes, metrics, evidence tiers and missing protocol elements. Generate compact 2x2 scientific figure(s) only if informative; do not draft paper.

This is a bounded numerical sanity/negative-control supplement, not E4, not an algorithm superiority claim or journal acceptance. Core theory and final interpretation remain with Codex.
