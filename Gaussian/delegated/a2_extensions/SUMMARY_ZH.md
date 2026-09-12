# A2 extensions 汇总（自动读取结果）

manifold_v3 已读取 30/30 cases；预算曲线只取 history 中 `rhs <= budget` 的最后点，缺失不插值。

概率 frozen v2 为 30 条；weak/posthoc 为 30 条。固定 RBF LP 可分 118/256。

An optimizer success flag only means the declared numerical optimizer terminated successfully; it does not establish posterior correctness. Weak and initialization results are post-hoc diagnostics and do not replace frozen v2.

Terminal actual RHS and seconds are reported separately because methods often terminate at different RHS; this is not described as equal-RHS comparison.

详细机器可读数据：`runs/a2/extensions/SUMMARY.json`。图：`figures/a2/extensions/probability_summary.png`、`manifold_budget.png`、`manifold_best_worst.png`。
