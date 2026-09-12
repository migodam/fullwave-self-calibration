# A2 ports 工作产物索引

本索引只说明现有 runner、结果和已知范围；不增加理论或评分。

| 主题 | Runner / 主要结果 | 范围与已知缺陷 |
|---|---|---|
| FFT VIE 与 component ports | `code/a2/physics.py`、`ports.py`、`run_ports_probe.py`；`runs/a2/ports/results.json` | 有限网格 2D 标量 VIE。普通 Gaussian 尾部不被当作紧支撑；端口不是自动功率归一 S 参数。 |
| cell quadrature 与 port rank | `run_quadrature_ports.py`；`runs/a2/ports/quadrature_ports_results.json` | 点核默认，cell-integrated 为独立有限网格检查。N16 rank/overlap 结果不是低秩充分性定理。 |
| N24 matrix-free local-T 回归 | `run_ports_n24.py`；`runs/a2/ports/ports_n24_regression.json` | 有完整运行配置、环境和 source hash；约 68% 散射误差是保留的负结果。旧无 manifest 文件在 `legacy_unmanifested/`。 |
| 有限 catalog 证书 | `finite_catalog.py`、`plot_finite_catalog.py`；`runs/a2/finite_catalog/results.json`、`figures/a2/catalog_error_bound.png`、`catalog_decision_counts.png` | 仅冻结的同算子有限字典；不推出连续材料类或独立成像结论。`invalid_v1/` 明确保存错误 covariance/nesting 的旧结果，不能引用。`INTERPRETATION.md` 说明 0.7σ pair 仍单峰且没有 naive false positive。 |
| 显式 current v2 | `explicit_current_v2.py`、`check_explicit_current_v2.py`；`runs/a2/explicit_current_v2/` | 自定义 TSOM 风格诊断；`D Pweak` 组合检查通过，但三例在 rank 128 的 state residual 仍约 0.12–0.13，未修复 closure。不可据此评价传统 TSOM。 |
| matrix-free GN v1 | `benchmark_matrix_free.py`；`runs/a2/matrix_free/benchmark.json` | 固定 Gaussian amplitude 的局部 shared-state 成本控制；v1 dense GN 已快于 range-CG。 |
| matrix-free GN v2 | `benchmark_matrix_free_v2.py`；`runs/a2/matrix_free_v2/benchmark.json`、`time_vs_parameter_rank.png` | 固定中心/宽度，p=784/3136，所有 rank 12/32/64 都报告。强 dual dense GN 仍更快；无 range-preconditioner 时间 crossover，停止数值调参。 |
| Fresnel 实测 | `measured.py`、`fixed_gain_refits.py`；`runs/a2/measured/results.json`、`figures/a2/measured_reconstructions.png` | 单一 2D TM 圆柱、固定几何、gain-profiled nuisance 诊断。主 run 未在开始时保存 source hash，字段明确为 null；十个起点只完整保留训练 objective 和最佳指标，不是配对统计比较。 |

共同边界：绿色检查、有限字典、单一实测对象或同一算子数值控制，均不构成新理论、论文接受、普适优势或生产部署验收。
