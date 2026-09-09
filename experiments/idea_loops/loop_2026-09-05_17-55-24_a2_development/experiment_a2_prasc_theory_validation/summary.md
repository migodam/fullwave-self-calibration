# A2 PRASC-SOM theory-validation package: summary

Bounded validation supplement. Executed E1/E2/E3/E5 synthetic modules and an optional physical-tangent run (N=8 full aperture). No E4, no nonlinear superiority claim, no publication acceptance claim.

- Test suite: 50/50 passed (E1 9, E2 11, E3 15, E5 14, physical smoke 1); results/test_suite_report.txt.
- Seed/sample ledger: E1 seeds 101-110 x 2000; E2 12 tangents (201-212), 36 checks; E3 seeds 301-310 x 1000; E5 seeds 401-412, 24 budget records; physical 12+12 seeds, declared C1=4/C2=8, no oracle rank.
- Headline results:
  - E1: matched contraction J_ph <= J_coh in all four scalar models; mismatch direct-intensity control reverses ordering (J=1960 vs 2) by design.
  - E2: loss identity residual 0; rho nonmonotone [0.75, 1, 1.55e-15]; 36/36 T5b; max backward-scaled residual 0.2734 (<100).
  - E3: dual spectra exact; K_e<=K_eL<=K0; Theorem-7 analytic risk 9.418 = 3.1997 var + 6.2183 bias^2, MC within 3 SE (max rel 4.1% var, 5.3% bias); confounding crossover at eps=1.
  - E5: shared-map stack>sum (2 vs 0); budget identity 24/24 (max 2.34e-14); random rank enlargements unsafe 0/24 (I<L) and correctly flagged; greedy gap 0.0845; non-submodularity gaps 3.04/7.60; gauge null persists, anchor restores.
  - Physical optional: T5b 24/24, max T5a scaled 0.01006, Theorem-3 residual 3.17e-11; innovation rel 4.48e-16; stack>=sum 12/12. Theorem 7 skipped 12/12 because declared C1 saturates the visible pose direction.
- Protocol ledger (full table in results/A2_quantitative_supplement.md): pass for common conventions and the synthetic E1/E2/E3/E5 core; not run: E1 physical Fisher scene, E3 limited aperture, E5 physical 8-candidate policy comparison, E5 source stabilizer, E4, publication evidence boundary.
- Artifacts: results/*_results.json and *_summary.csv (raw), figures/*.png, results/ci_stats.json, results/test_suite_results.json, results/A2_quantitative_supplement.md.
