# Bounded A3 numerical validation worker

Read AGENTS.md, Theory/Questions/A3.md sections Theorems 1–4 and correct fixed-chart derivative (around lines 141–660, 775–820). SOM means subspace-based optimization in electromagnetic inverse scattering, NOT self-organizing maps. Do not start child workers or another pipeline. Codex owns proofs and final claims.

Write code, raw JSON and concise summary in research/delegated/a3_algebra/ only (pipeline internal workdir artifacts also allowed). No changes to A2, source Theory files or main A3 manuscript. CPU single thread, NumPy/SciPy existing Python experiments/idea_loops/loop_2026-09-04_02-58-31/experiment_geometry_lifted_trispace_som/.venv/bin/python. Use apply_patch for source edits. No installations.

Implement and execute deterministic falsification tests:
1. Fixed-chart exact state plus selected tangent span suffices; explicit e1 + theta e2 state-only counterexample. Complex nonzero-residual LS coefficient derivative must include (M_v U)^* r. Compare correct formula and omitted-residual formula with centered finite differences in random well-conditioned systems; moving U derivative control.
2. Nonzero rectangular complex H: row/column gain orbit invariants, realified whitened quotient tangent projector identity for NON-identity covariance, rank-one gain-hiding and two-component formula. Reject/flag near-zero entries rather than divide blindly. Do not claim finite-noise ratio sufficiency.
3. Independent Gaussian branch validation bound. Freeze candidate predictions and acquisition BEFORE drawing validation noise. At least 10000 repetitions for conditional two/three-candidate cases with declared bounded mismatch. Include missing-candidate and reused-validation controls; no global coverage claim. beta is explicitly assumed/known in algebra fixtures, not inferred from fitted residual.
4. Dual residual identity including complex conjugates, certified resolvent bound against exact offline SVD, and near-singular small-residual negative control.

Use reproducible seeds 3101–3110; tolerances scale with conditioning. Include exact failing controls, not only positive tests. Produce checks.json with case-level inputs/metrics/thresholds, a script to regenerate, and SUMMARY.md. Keep stdout <=10 findings, paths, uncertainties. This bounded validation is neither independent 3D verification nor SOM-performance evidence. Do not draft an ACL paper or rename the supplied task on novelty grounds.
