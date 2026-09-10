# B1：移动近场接收阵列的定量成像与动态校准

## 提交范围警告

GitHub平台拦截了主实验脚本experiment.py的上传。当前PR仅包含成功写入的文件，不能只靠PR独立复现；完整可运行代码、原始观测和每个初值通过本轮对话ZIP附件交付。没有改编码、换路径或换工具规避该拦截。文档中的复现命令针对完整附件。

## 本轮裁决

**NO — known-scatterer/reference calibration remains the better solution**

此结论限定于本轮选定的机械扫描、已知双区域支撑、共享低维电子学模型及所执行预算。它不是所有移动阵列自校准的不可能性定理。目标无参考的两种候选没有超过强参考方案；增益消元与相同先验的普通联合反演是同一统计目标；独立数据还暴露了显著正演失配。未获得强新 self-calibration 论文的证据，不能标为 TAP 投稿完成。

## 范围和成果

固定两个发射源，两个接收探针组成一个移动接收头；八帧、三个频率；两个静止、已知支撑、未知材料的区域。未知状态为每帧公共复增益、频率线性延迟、有效法向相位中心位置。不分离机械位移和相位中心偏差，不扩展至 UAV、未知形状或稠密 SLAM。

完成三场景比较、近邻检索、物理模型、四类解析歧义/界、两个候选、九主方法、真正两阶段参考校准、因果前缀、自变量/物理模型消融和三维独立求解器检验。36 个开发场景、324 个主方法拟合；72 个两阶段诊断、72 个因果扫描诊断、84 个模型诊断。方法内多个初值是优化检查，不是独立实验样本。

主结果 `07_EXPERIMENTS/SUMMARY.json`、`OUTCOMES.csv`；原始观测与每个初值位于同目录下 S/I/X 等子目录。所有结果不与 Q5 或 A1–A5 合并。完整归档含原始数组；GitHub 中若只提交紧凑摘要，则不声称该摘要是原始数据。

## 阅读顺序

1. `01_ENGINEERING_SCENARIOS.md`：为什么选择机械扫描，为什么不能仅因阵列移动就宣称必须自校准。
2. `02_PRIOR_ART.md`：八项假设维度、全文状态与真实差距。
3. `03_PHYSICAL_MODEL.md`、`04_IDENTIFIABILITY.md`：参数、噪声、sharing、证明和作用域。
4. `06_BASELINES/ALGORITHMS.md`、`07_EXPERIMENTS/REPORT_ZH.md`：实际实现和逐项证据。
5. `08_HARDWARE_PLAN/COST_AND_FEASIBILITY.md`、`09_REVIEWS/ADVERSARIAL_REVIEW.md`、`10_FINAL_VERDICT.md`。

## 复现

Python 3.13.5、NumPy 2.3.5、SciPy 1.17.0；没有学习模型，也不需要 GPU。原始运行使用每进程单 BLAS 线程；多个诊断进程曾并行运行，因此本轮秒数不是干净硬件基准。

```bash
python -m pip install -r requirements.txt
python reproduce.py --output /tmp/B1-smoke --mode smoke
python reproduce.py --output /tmp/B1-full --mode full
```

输出目录必须不存在。脚本把代码复制到新目录后运行，不改动交付的原始数组。smoke 只重建 S/I 第零场景及两种方法，不能说成再次完整复现全部 36 场景。已完成一次干净 smoke：全部原始数组逐元素相同，四个拟合的材料误差差值为零，17 项测试通过。full 运行三层全部主方法与诊断，仍是固定种子重放，不是新的确认性测试。

## 已知限制

模型与原始样本的对照不是实际硬件；DDA 网格差与 VSW 阶数差不是连续 Maxwell 误差界。所有 I/X 相干主方法均触发失配报警，即使材料数值达标也不算经过认证的真实校准。比较计入复标量采集数，没有得到真实停机时间、动态漂移谱或完整硬件成本，所以不能宣称成本 Pareto 最优。

期刊审阅所需的源码归属在 `05_CODE/SOURCE_ATTRIBUTION.md`。A1–A5 只读。
