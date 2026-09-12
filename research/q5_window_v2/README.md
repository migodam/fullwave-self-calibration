# Q5 有限窗口与模态读出：可复现审计包

中文汇报与证据表：`MANUSCRIPT_ZH.md`。英文论文修订段落与证明：`MANUSCRIPT.md`。有限割线／径向分裂补充：`WINDOW_PROOF.md`。逐篇文献公式审查：`LITERATURE_AUDIT.md`。

**状态：受限模态结论有条件成立；同时双有损球的 TAP 强结论未完成。** 没有恢复失败选择器，没有宣称新采集策略胜利，没有修改 A5/V1。`PROTOCOL.md`、`EXECUTION_FREEZE.md`、`MODAL_READOUT_PROTOCOL.md` 为执行前注册；提交分支从注册提交 cbebe706 分出，避免与原工作分支后续出现的不同实验报告混合。

## 运行

需要 Python；实际验证环境 Python 3.13.5，其他版本未验证。依赖见 `requirements.txt`。从仓库根目录安装依赖后，运行：

```bash
python -m pip install -r research/q5_window_v2/requirements.txt
bash research/q5_window_v2/reproduce.sh /tmp/q5-new-evidence
```

目标目录必须是新的。脚本复制必要代码和仓库现有、未修改的 `research/trispace_self_calibration/a3_research/maxwell3d.py`，依次生成全部证明记录、独立恢复、机制诊断、有限竞争世界、窗口证书、测试日志与汇总。不要覆盖原始证据。在附带完整结果的归档中，可直接运行：

```bash
OPENBLAS_NUM_THREADS=1 python -m pytest -q research/q5_window_v2/tests
```

25 项测试已通过。主比较每方法每场景 3 个初值、每初值最多 80 次函数评估。原运行主恢复约 11.75 秒、峰值 RSS 288.3 MiB；完整重放还包含更多诊断，时间与机器有关。不需要 GPU、外部 API 或 Treams。

## 文件职责

`interval_certificate.py` 用整数区间和有理输入覆盖 5500 个材料盒；`readout_certificate.py` 输出有限采样泄漏及条件误差预算；`finite_windows.py` 覆盖有限割线和理想径向窗口；`multipole.py` 为标准向量 Mie 多球模型；`run_experiment.py` 为冻结的新 12 场景流；`diagnostics.py` 的 oracle 仅作事后解释；`finite_pair.py` 的候选不是认证反例；`modal_estimator.py` 实现受限读出与单调反演；`tests/test_q5.py` 提供回归检查。

`results/VERIFIED.json`（提交摘要）与重放生成的 `results/REVIEW.json` 是便于审阅的紧凑结果，全部原始复观测、5500 行覆盖账本、主比较初值、失败和日志在完整归档中，亦可由脚本生成。不同机器的计时和 JSON 整体哈希可能变化；物理输入、数组、源文件哈希及确定性数值应分别比较。不要把紧凑结果说成已经包含全部原始数组。

## 边界

模态参考是实幅度高斯测量，独立双球参考是 proper complex IQ 测量，二者不可混用。18 方向定理只在 kR=2；全距离窗口例子采用理想完整模态投影。DDA 网格差、Mie 阶数差、优化距离和局部奇异值都不是连续认证。真实仪器尚未证明达到所列规格。
