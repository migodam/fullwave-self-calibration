# Claim check: report/tables vs raw JSON and code

Paths below are relative to
`experiments/idea_loops/loop_2026-09-04_02-58-31/experiment_geometry_lifted_trispace_som/`
unless noted. `report` =
`../final_2026-09-04_06-31-18_geometry_lifted_retained_range_som/experiment_report.json`.

## Bookkeeping-level checks

| Check | Evidence | Status |
|---|---|---|
| `consolidated_round3_tables.md` cells match stored medians/IQRs in the raw JSONs | Rerun of `make_consolidated_summary.py` byte-identical (`cmp` exit 0; hashes in `reproduction.md`) | PASS |
| Report numeric values match the consolidated summary/raw JSON | Every prose number checked below is present in raw JSON at 4-sig precision; exactness claims are the exception | PASS with wording caveats |
| Figures referenced by the report match experiment figures | SHA-256 identical for the four final-figure copies | PASS |

## Report claim-by-claim check

| Report claim / audit question | Raw/code evidence | Status |
|---|---|---|
| Unrestricted lift reproduces perturbation to machine precision | `settings.*.lift_checks.unrestricted.residual_ratio_full` in `results_e10_settings_sweep.json`: ~5e-16 to ~2e-15 | PASS |
| Retained r4/r6 residual ranges 0.474-0.864 / 0.354-0.686; QH R ~1e-18 | `settings.*.lift_checks.retained[*]`: exact ranges and 1.8e-18-1.2e-17; independent recomputation matched to 0.0 | PASS |
| "reduced_r4 matches direct to full precision in every setting" | Half circle: direct med 0.010429502, reduced_r4 med 0.010209086 in `consolidated_round3_summary.json`; M12_k16 full-circle max pose-error difference ~1.6e-7, not machine precision | FAIL / overclaim |
| "reduced_r6 also matches direct everywhere pose is visible at r6" | Half circle has r6 `n_vis=3, hidden_rank=0`; direct med 0.010429502 vs reduced_r6 med 0.027149679; seed 0 pose error 1.286716 | FAIL / overclaim |
| M8 canonical regime: r6 n_vis=1/hidden_rank=2, med 0.0442 vs direct 0.0158 | `results_e10_settings_sweep.json` M8: hid_sv [5.679e-4, 5.37e-19, 1.71e-19], medians 0.044172/0.015824 | PASS |
| M12/M16 transition at r=M-2 with listed hidden medians | `results_e10b_m12_m16_highrank.json` transition rows r10/r14 and run medians (0.075896/0.072033, 0.110962/0.090503) | PASS |
| Trailing r=M-1 and r=M empty-lift: B_red ~1e-19, reduced matches direct exactly | Rank table marks r=M-1 and r=M `vanished_Bred`; **only r=M-1 nonlinear runs exist**; M12 r11 max p diff 1.14e-4, map diff 3.10e-5; M16 r15 max diff ~5e-16 | FAIL / overclaim |
| Oracle baselines and ranges | `known_pose`/`known_alpha` rows are oracle parameterizations (truth only inside named methods); stated map/pose ranges match full-circle and half rows; threequarter `known_pose` map med is 0.006409, so "0.008-0.012" omits one tested setting | PASS with scope note |
| Half-circle local minima recorded | `results_e10_settings_sweep.json`: `known_alpha` pose ~0.857-0.861 on all seeds, residuals 0.266-0.278; reduced_r6 seed0 pose 1.286716, residual 0.039564, all status=1 and optimality<1e-7 | PASS (residuals not shown in md tables) |
| E6b 240/240 runs converged; listed medians | `results_e6b_multitx_robust.json`: 240 runs, 240 success; medians match reported values | PASS; note success criterion is loose (`status>0`), no status stored |
| E6b extension medians, FD 2.169e-7 at eps=1e-4, dipole gauge | `fd_checks_dipole`, `gauge_partA_dipole`: values match; "third realified SV 2.0e-4" is the unnormalized third `realified_sv` 2.0299e-4, while the table prints column-normalized third SV 0.0583 | PASS with unit/label caveat |
| E11 hard null and soft partial-bias numbers | `results_e11_vp_state_null.json` medians: hard 0/2 hidden 0.04408; perturbed hard hidden 0.1959; soft hidden 0.01604, VP residual 0.06229, discard 0.8693, map 0.83799, full-physics residual ~0.83 | PASS |
| T_U remains diagnostic-only | `results_e5_final.json`/`results_e7*.json`: T_U ~0.81-0.89 across good and bad pose solutions with no discriminating separation | PASS |

## Question-by-question audit summary

1. **Do report and tables agree with raw JSONs?** Yes at the rounded table
   level (generator byte-reproducible); report prose exactness claims go
   beyond the raw numbers in the cases above.
2. **Lift formulas implemented as stated?** Yes: retained lift
   `Q=G_s V_r`, `c=pinv(Q)b`, `R=Qc-b`, QH R orthogonality; unrestricted
   `pinv(G_s)`; realified column blocks; 1e-10/1e-8/1e-12 thresholds;
   complement projector; all independently recomputed with zero difference.
3. **Receiver-only/transmitter-only/co-moving derivatives?** Yes: E3 rerun
   byte-identical; E6/E6b FD tables show O(eps^2) agreement on the same
   forward models (`results_e1_e3.json`, `results_e6.json`,
   `results_e6b_multitx_robust.json`).
4. **Distinct baselines / oracle usage?** Yes: parameterizations are distinct;
   truth is used only in observation synthesis, metric evaluation, and the
   `known_pose`/`known_alpha` named oracles.
5. **Is reduced-matches-direct nontrivial?** No: for n_vis=3, V_vis spans R^3
   and reduced minimizes the same full-physics objective under an invertible
   pose reparameterization. It proves implementation consistency and, where
   observed, same-basin convergence; it does not establish a "lossless"
   scientific benefit. The half-circle reduced_r6 counterexample shows basin
   selection is not guaranteed.
6. **Origin of r=M-2 transition?** Fixed nominal basis at
   `(alpha_init,p_init)`; dimension-count/rank bound (orthocomplement
   dimension 1 -> rank(B_red)<=1); robust SVD gap at transition; the
   empty-lift r>=M-1 state is the threshold/numerical artifact and is
   separately flagged.
7. **Do convergence flags hide local minima?** Yes: half-circle records pass
   the gate with residual ratios ~8-12x the noise ratio. Report flags the
   caveat; tables omit final residual columns.
8. **State witness / VP negatives with compatible normalizations?** Negative
   conclusions are supported. VP "data residual" is VP-model residual while
   direct/reduced rows use full physics; JSON notes explain this, but the md
   table lacks the caveat, and hard-row medians include failed runs.
9. **Seeds/deps/commands/regeneration?** Fixed seeds and no global RNG;
   runners write and assert JSON/PNG; `.venv` versions recorded; no
   requirements/command manifest; summary generator and E1/E3 rerun
   byte-identical. Full 240-case sweeps were not rerun per scope.
10. **Leakage/post-selection/silent handling?** None found in audited code:
    no truth in non-oracle residuals/Jacobians, no reuse of prior results
    inside solvers, no post-selection of successful runs in summary medians,
    no silent data exception path in results (failures are explicit records).

