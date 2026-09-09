Checked VI.B against `reuse_n16_q9_r128.json`.

| Item | Manuscript | Summary / required correction |
|---|---|---|
| Generic maximum map error: rounding mismatch | \(4.69\times10^{-4}\) | Exact: 0.000468243703994622. At three significant figures: **\(4.68\times10^{-4}\)**. |
| Full forward RHS counts: missing mean qualifier | “2034 for direct GN to 279 for either ROM” | These are `full_rhs_mean` values. Specify **mean counts per scene**. |

All other stated numbers agree: 16 runs (four methods × four scenes), 26 audited updates and passes per ROM, 52 combined passes, mean reuse hits 3.75/3.5, mean losses rounding to 213.416, and online times rounding to 3.15/2.92/1.00 s.

Twofold maximum map error correctly rounds to \(6.18\times10^{-4}\). Maximum pose tangent errors, omitted from the paragraph, are 0.00018041484586684968 (generic) and 0.00022013807766323585 (Twofold); both support the stated \(10^{-2}\) tolerance claim.