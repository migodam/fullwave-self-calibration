# A2 最终数字口径 QA

核对对象：`deliverables/RESEARCH_REPORT_ZH.md` 的非 imaging 数字，以及各已保存 `runs/a2` JSON。此文档只做口径核对，不构成科学验收。

## 逐项核对

| 报告口径 | 结果文件口径 | 结论 |
|---|---|---|
| 每组件 r=24 的三种散射误差约 2.05%、3.85%、6.37% | `ports/quadrature_ports_results.json` 为 2.0540%、3.8484%、6.3732% | 通过 |
| N24、9 维压缩网络误差 68.117% | `ports/ports_n24_regression.json` 为 68.1169% | 通过 |
| 圆柱 N32/N64/N128 约 2.47%、1.45%、0.603% | `ports/quadrature_ports_results.json` 为 2.4729%、1.4516%、0.6026% | 通过 |
| matrix-free 四行时间：p49 2.24/26.00/16.98，p196 9.39/40.02/27.30，p784 2.05/39.67/10.12，p3136 4.67/50.00/14.96 秒 | `matrix_free/benchmark.json` 与 `matrix_free_v2/benchmark.json` 对应值一致；均为报告所列方法及包含的 setup 口径 | 通过 |
| Fresnel gain-profiled 留出百分数 | `measured/results.json`：1G 10.216/15.559/21.788%，2G 10.222/15.304/20.946%，4G 10.586/15.163/20.561% | 与报告四舍五入值一致 |
| Fresnel 十起点范围 | 每模型保留 10 个 `paired_start_objectives`；最佳指标为单一选中分支。1/2/4G 所选分支均为 45 次上限停止；2/4G 的全部成功计数为 0。 | 报告“不作十组完整配对统计、所选结果到上限”一致 |
| finite catalog 旧 9、新含旧并扩展至 27、三真值各十噪声 | `finite_catalog/results.json` 有 3×10 seed、每 stage 36 个候选向量：9 个旧类候选 + 27 个嵌套的新类候选。 | 通过；`36` 是两个类比较时的候选总行数，不与“新类 27”矛盾 |

## 数据索引链接

`data/a2/DATA_INDEX_ZH.md` 的 `../README.md`、`../../runs/a2/measured/results.json` 与 `../../delegated/a2_ports/MEASURED.md` 相对链接均解析到当前 Gaussian 工作树中的存在文件。SHA 与 `data/README.md` 及 measured 结果一致：`476cc9d1cfc98797545ab4adf69302dc5aeb45848a24cf8d7b3d222940cc79eb`。

## 真实矛盾与建议

未发现上述指定数字的真实矛盾。保留报告已有边界：finite catalog 的 36 行候选总数不应被简写成“27 个总候选”；Fresnel 主运行的源码哈希仍为 null，不能事后回填。
