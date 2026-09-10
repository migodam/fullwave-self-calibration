# Reproducibility and execution boundary

## Final native run status and audit

The native process exited 0 after one development round and the 18-turn cap; evaluator decision: abandon. Native writeup was not executed. Worker results/FINAL_REPORT.md and FINAL_SUMMARY.json now exist; use the corrected parent theory/review to interpret them.

src/integration_checks.py checks artifact hashes, a repaired Born counterexample and derivative identities; results/integration_audit.json records the result. The worker validator also reported OVERALL PASS locally for structure, counts and hashes. It references excluded raw logs, so its full success is not claimed for the public clone. PUBLIC_MANIFEST.json instead verifies exported payloads. The correct scientific status remains incomplete.


Run from the Inv_SLAM project root with research/trispace_self_calibration/a3_research/.venv3d/bin/python and one BLAS thread. Required installed packages include NumPy, SciPy and Treams. The development scripts depend on the unchanged A3 maxwell3d.py. This is a workspace reproduction, not yet a self-contained release package.

Order: src/recovery_cycle1.py, src/scalar_cycle2.py, src/audit_development.py, src/parent_checks.py. Output files are intentionally protected against overwrite; use a separate copy/output directory for reruns and preserve the published source/protocol hashes. Cycle-one raw observations are regenerated from its fixed stream; cycle two stores its independent observations explicitly. The corrected audit, not obsolete cycle-one heuristic flags, governs reporting.

Protocols are docs/DEVELOPMENT_PROTOCOL.md and docs/CYCLE2_PROTOCOL.md. Results are results/development_cycle1.json, development_cycle2.json, development_audit.json and parent_checks.json. Three initializations per scene are not independent scenes. No held-out final data exist yet.

The real installed Agentic-AI-Scientist development loop is invoked by src/run_native_checks.py, with task context docs/WORKER_TASK.md and outputs under pipeline/boundary_20260910_171833/. Native paper writing is replaced by a parent-review handoff, not claimed as executed. External calls use only the verified canonical deepseek-flash endpoint for DeepSeek V4.1 Flash; provider_preflight.json preserves the initial unavailable spelling and provider_smoke.json records the successful minimal call. No Pro fallback or default-provider change is authorized.

Worker interval code currently uses an existing cached mpmath package. Portable dependency packaging and complete parent audit remain required. Raw worker output is provisional and is not a scientific conclusion.

An additional independent scalar risk reconstruction runs with src/parent_boundary.py and writes results/parent_boundary.json. It reproduces the 44.3423% supplied binary error example; reference-budget arithmetic is explicitly conditional, not an interval certificate.
