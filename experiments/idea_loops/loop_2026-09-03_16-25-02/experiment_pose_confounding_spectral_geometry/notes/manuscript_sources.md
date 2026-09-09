# Manuscript source mapping (notes/synthesis.md, manuscript/main.tex)

Date: 2026-09-03. This file lists, for every manuscript table/claim group, the
exact result JSON and/or notes files from which the numbers were taken. No
experiment was rerun during manuscript production.

## General environment / audit constraints

- Parent audit gates binding the runs and the write-up:
  `context/PARENT_ROUND_AUDIT.md`, `context/PARENT_CORRECTIONS.md`,
  `context/PARENT_FAMILY2_REVIEW.md`, `context/PARENT_FAMILY3_CORRECTIONS.md`,
  `context/PARENT_FAMILY4_AUDIT.md`.
- Abstract / idea wording and scope boundary: `current_idea.json` in the loop
  root (parent directory) and `context/CONTEXT_SUMMARY.md`,
  `context/EXPERIMENT_PLAN_EXTRACT.md`.

## Section 2 (model/operators) sources

- Equations for g_k, G_D diagonal/self-cell, A_R/B_R, K_IS/K_SLAM/K_eff/L_X,
  rho = 1 - sigma_i^2, prior and duplicate rules:
  `context/PARENT_ROUND_AUDIT.md` (formulas), `context/PARENT_CORRECTIONS.md`
  (self-cell), `context/EXPERIMENT_PLAN_EXTRACT.md` (theorem list 2 and
  realification recipe), `context/CONTEXT_SUMMARY.md` (notation and semantic
  bans).
- Self-cell v2 implementation values (validation ladder, correction report):
  `notes/family1_self_cell_correction_report.md`,
  `results/self_cell_quadrature_validation.json`.

## Family 1 tables and claims

Sources:
- `notes/family1_self_cell_correction_report.md`
- `results/family1_pilot_results.json`
- `results/family1_grid_refinement.json`
- `results/self_cell_quadrature_validation.json`
- `src/test_family1_corrections.py` and `src/validate_self_cell.py` are the
  executed runners referenced by the report; raw numbers live in the JSONs
  above and in the notes tables.

Table/claim mapping:
- Self-cell formula and quadrature: report sections 1-5, validation JSON rows.
- Unit-test status 5/5 pass and sign-discrimination rel error 2.0:
  report section 4.
- Map/pose FD errors, slopes, state margin, Born table: pilot JSON top-level
  fields `map_fd`, `pose_fd`, `convergence_slopes`, `state_check`,
  `born_fullwave`; report sections 6-8.
- Grid refinement rows (N=16/24/32/40, sigma_min(M)/norm, runtime):
  `results/family1_grid_refinement.json` (`per_grid`, `comparisons`,
  `summary`, `metadata`).

## Family 2 tables and claims

Sources:
- `notes/family2_report.md` (claim-status table, cannot-establish)
- `notes/family2_parent_correction.md`
- `results/family2_results.json` (cases.pixel / cases.smooth:
  `check_c1_kernel_identity`, `check_c2_rank_identity`,
  `check_c3_retention_spectrum`, `check_c4_loss_rank`,
  `check_c5_ordering_interlacing`, `check_c6_prior_sweep`,
  `check_c7_singular_prior`, `check_c8_duplicate_invariance`,
  `resolution_smooth_c1c2c3`, `ranks_summary`, metadata)
- `results/family2_parent_correction.json` (stable factorization residuals,
  SVD-filter sweep d0 at alpha=1e-14, explicit singular-prior limit)

Specific numbers:
- c1 pixel gate 8.47e-7 / identity-support evidence:
  family2_results.json `cases.pixel.check_c1_kernel_identity` and
  family2_parent_correction.json `cases.pixel.original_pixel_projector_gate`.
- c3 residuals 2.89e-15 / 1.55e-15:
  family2_results.json `cases.*.check_c3_retention_spectrum.max_abs_diff_rho_vs_pred`.
- Stable factorization 3.77e-16 / 2.70e-16, alpha=1e-14 d0 5.70e-9 / 3.67e-9,
  singular-prior limit 4.71e-4 / 2.29e-4:
  family2_parent_correction.json.
- c6/c7/c8 sweeps and interlacing/c4 ranks:
  family2_results.json check blocks above; report claim table for rendered
  rows.

## Family 3 and 3b tables and claims

Sources:
- `notes/family3_report.md` (raw claim-status table)
- `notes/family3_parent_correction.md`
- `notes/family3b_report.md` (refined claim-status table)
- `results/family3_gauge_born.json`
- `results/family3_parent_correction.json`
- `results/family3b_refinements.json`

Specific numbers:
- Smooth gauge residuals/representation residuals/kernel distances:
  family3_gauge_born.json `gauge_two_blob_smooth.rows` and
  `gauge_refinement_smooth` (N=16/32/40).
- Pixel gauge residuals 6.61e-5 / 1.58e-5 / 4.68e-5:
  family3_gauge_born.json `gauge_two_blob_pixel.rows` (also N rows in
  family3b_refinements.json `check_A_pixel_gauge_vs_N`).
- Symmetric scene rotation: family3_gauge_born.json
  `gauge_symmetric_smooth.rows`; interpretation in
  family3_parent_correction.json `symmetric_rotation_stabilizer`.
- Anchor rho_min ratios: family3_gauge_born.json `anchor_variants.rows`.
- Born empty background B=0 / K_SLAM0=K_IS0 / bilinear dominance / raw FD
  failure: family3_gauge_born.json `born_empty_background`.
- Raw failure preservation flags: family3_parent_correction.json
  `original_failures_preserved` and `pixel_gauge_refinement.rows`;
  `mask_experiment_interpretation` for the mask rows.
- 3b refined checks: family3b_refinements.json `check_A_pixel_gauge_vs_N`,
  `check_B_augmented_smooth_gauge`, `check_C_anchor_generator_retention`,
  `check_D_born_bilinear_corrected`; report narrative in
  notes/family3b_report.md.

## Family 4, 4b, and parent-control tables and claims

Sources:
- `notes/family4_report.md`
- `notes/family4b_normalized_control.md`
- `notes/family4_parent_controls.md`
- `results/family4_results.json`
- `results/family4b_normalized_control.json`
- `results/family4_parent_controls.json`

Specific numbers:
- PSD increments, duplicate controls, shared-z rows:
  family4_results.json `checks.A_psd_monotonicity`,
  `checks.B_duplicate_controls`, `checks.C_shared_pose_compensation`.
- Trajectory counterexample tables (retained/confusable mass, theta_min,
  log_volume, rho_min, row norms): family4_results.json
  `checks.D_equal_budget_trajectories` (trajectories, observed orders,
  hypothesis_holds=false) and report section tables.
- Normalized control rows: family4b_normalized_control.json `checks` and
  `pass_summary`.
- SNR-matched retentions 0.1170/0.2157/0.4460 and trajectory control orders:
  family4_parent_controls.json `frequency_snr_matched`,
  `trajectory_controls.*`; notes/family4_parent_controls.md.

## Family 5 tables and claims

Sources:
- `notes/family5_report.md`
- `results/family5_sensitivity.json`

Specific numbers:
- P1 rows and at-eps-1e-3 values 6.25e-7 / 4.44e-2 and slopes:
  family5_sensitivity.json `checks.part1` (`rows`, `at_eps_1e-3`,
  `clean_window_full`, `clean_window_frozen`, `dW_term_share_of_DK_full`).
- P2 gated eigenvalues/gap flags: `checks.part2` (`gate_rows`, `chosen`,
  `gap_flags`, `simple_condition_flagged_indices`, `pass`).
- P3 rows/slopes/sampled-min: `checks.part3`.
- P4 empirical-only Lipschitz numbers: `checks.part4`.
- P5 projector movement/scan/crossing: `checks.part5`.

## Related-work paragraph

Citation identifiers were verified against local knowledge-bank JSON records
in the loop root (`knowledge_bank/knowledge_bank/`) and web records checked on
2026-09-03:
- Ling & Strohmer arXiv:1611.04196 (local record matches).
- Entekhabi & Isakov arXiv:1712.08696 (local record matches).
- Shea et al., Medical Physics 37(8), 2010, DOI 10.1118/1.3443569 (web check).
- Weidling & Hohage arXiv:1512.06586 (web check).
- Burnett et al. 2022 arXiv:2203.10174 (web check).
- Hilger et al. 2024 arXiv:2404.03940 (web check).

## Appendix claim-status tables

The appendix tables in `manuscript/main.tex` replicate the claim-status tables
of `notes/synthesis.md` sections 1.1-1.6; each row maps to the sources listed
above for its family.
