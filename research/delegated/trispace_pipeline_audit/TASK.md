# Mechanical reproducibility and code audit

Work in `/Volumes/migodam's-external-brain/Research/Inv_SLAM`. This is a
bounded code/results audit. Do not make the final scientific or novelty
judgment. Do not modify the experiment directory.

Target pipeline:

- `experiments/idea_loops/loop_2026-09-04_02-58-31/`
- experiment workdir:
  `experiments/idea_loops/loop_2026-09-04_02-58-31/experiment_geometry_lifted_trispace_som/`
- final report:
  `experiments/idea_loops/loop_2026-09-04_02-58-31/final_2026-09-04_06-31-18_geometry_lifted_retained_range_som/experiment_report.json`

Audit questions:

1. Do the numerical values in the final report and
   `consolidated_round3_tables.md` agree with the raw result JSON files?
2. Are the unrestricted and retained lift formulas implemented as stated,
   including realification, rank threshold, projector, and orthogonality?
3. Are receiver-only, transmitter-only, and co-moving analytic derivatives
   tested against the same forward model without sign or double-count errors?
4. Are `direct`, `reduced_r4/r6/...`, `wrongpose`, `known_pose`, and
   `known_alpha` genuinely distinct baselines, and are oracle variables used
   only in the named oracle baselines?
5. Is “reduced matches direct” a nontrivial outcome or a parameterization
   identity in the tested visible regime? State exactly what it proves.
6. Is the `r=M-2` hidden-rank transition computed from a fixed nominal basis,
   a pose-dependent basis, a threshold artifact, or a symmetry artifact?
7. Do convergence flags hide bad local minima? Check objective/residual and
   parameter errors, especially the half-circle case.
8. Are state-witness, hard variable-projection, and soft variants correctly
   reported as negative, with compatible residual normalizations?
9. Check deterministic seeds, dependency isolation, commands, and whether all
   reported output files can be regenerated. Do not rerun every 240-case sweep;
   rerun cheap core checks and the consolidated-summary generator.
10. Look for metric leakage, use of ground truth inside non-oracle solvers,
    accidental result reuse, silent exception handling, post-selection, or
    formulas evaluated at truth when claimed as implementable diagnostics.

Write only under `research/delegated/trispace_pipeline_audit/`:

- `summary.md`: short verdict, key verified findings, key risks.
- `claim_check.md`: report-claim to raw/code evidence table.
- `reproduction.md`: exact commands run, exit codes, regenerated hashes or
  numerical comparisons.
- `issues.md`: bugs, ambiguities, tautologies, normalization problems, or
  missing evidence, prioritized P0-P3.

Do not edit source/results. Final stdout: success/failure, 3-10 findings,
artifact paths, uncertainties, and whether parent Codex intervention is needed.
