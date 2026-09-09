# Worker summary: Family-1 self-cell correction and revalidation

Status: SUCCESS. Run 2026-09-03T09:26Z (UTC) in
`experiments/idea_loops/loop_2026-09-03_16-25-02/experiment_pose_confounding_spectral_geometry`.
No worker agents were launched; work executed directly per the task constraint.

## Findings

1. Pre-correction snapshot archived verbatim under
   `experiment_pose_confounding_spectral_geometry/archive/INVALIDATED_self_cell_2026-09-03T092339Z/`
   with README + checksums; the live old `logs/environment.md` is superseded,
   not deleted.
2. `src/helmholtz.py` `self_cell_green` corrected to the complete equal-area
   disk integral
   `I_self = (i*pi*a/(2*k_b))*H_1^(1)(k_b*a) - 1/k_b^2` (`a = h/sqrt(pi)`),
   and docs updated; `[G_D]_nn = k_b^2*I_self` now tends to 0 with h.
3. Corrected source-gradient sign preserved via new `green_grad_source`
   helper (`grad_s = -grad_z`); build_AB refactor verified bit-identical to
   the prior validated expression (max abs diff 0.0).
4. 5 new unit tests pass; a mutated old-sign probe fails at relative error
   2.000e+00 vs the finite difference (correct sign 1.737e-09).
5. Radial-quadrature validation (2048-node Gauss-Legendre, N = 8..128): max
   analytic-vs-quadrature abs error 3.13e-16; |I_self| falls monotonically
   5.123e-03 (N=8) to 4.358e-05 (N=128); invalidated rule saturates at
   ~1/k_b^2 ~= 2.53e-02.
6. Family-1 pilot rerun clean: max map FD rel err 1.8e-08 and pose 4.3e-06
   at h_fd=1e-3; slopes ~1.98/2.00; `sigma_min(M)/||M|| = 0.719`; result JSON
   marks corrected formula/version, command, runtime, environment, hashes.
7. Forward-grid refinement (N = 16/24/32/40, same continuous scene and
   receivers) shows decreasing differences (1.45e-03 -> 1.60e-04 vs finest);
   labelled a discrete-model diagnostic, not continuum convergence. N=32 and
   N=40 completed, not deferred.
8. No unresolved mathematical ambiguity remained; all numbers, commands,
   versions, seeds, tolerances, runtimes, and SHA-256 digests are recorded in
   the experiment notes/log and result JSONs.

## Key artifact paths

- Experiment root:
  `experiments/idea_loops/loop_2026-09-03_16-25-02/experiment_pose_confounding_spectral_geometry`
- Corrected source: `src/helmholtz.py`
- Unit tests: `src/test_family1_corrections.py`
- Pilot rerun: `results/family1_pilot_results.json`
- Quadrature validation: `results/self_cell_quadrature_validation.json`
- Grid refinement: `results/family1_grid_refinement.json`
- Report: `notes/family1_self_cell_correction_report.md`
- Environment log: `logs/environment_self_cell_corrected.md`
- Invalidation archive: `archive/INVALIDATED_self_cell_2026-09-03T092339Z/`

## Uncertainties

- Passes are bounded checks of the implemented discrete harness (internal
  derivative consistency plus direct quadrature fidelity of the audited
  self-cell formula); they are not continuum-convergence or physical
  acceptance evidence.

Codex intervention needed: none for this correction gate.
