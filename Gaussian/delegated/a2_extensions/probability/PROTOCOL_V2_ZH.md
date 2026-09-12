# A2 概率成像 v2 协议

旧 `probability_run.py` 结果保留为 pilot，不用于本协议结论。v2 的材料空间是一个**已知**的八个 Gaussian-centre Voronoi patch 有限字典；每个 patch 仅有二元标签。因此它不是任意像素重建，也不是盲 Gaussian 成像。

生成端对 30 个先验预测对象（每个噪声组 10 个）以 N64、`cell_integrated=True`、两频率 1.5/2.75 GHz、4 Tx、24 Rx 作完整 VIE 前向。模型端复用冻结的 256 状态 N32 全波场缓存。偶数接收机为拟合通道，奇数接收机完全留出。先验为全支持的环形 Ising；全空状态没有被删去。噪声为 `fraction * max(truth RMS, all-material N32 RMS)`，fraction 是 0.03、0.1、0.3，故不应称为固定 SNR。

精确基准枚举 256 个状态。Product VI 与 Gaussian-RBF-logit VI 都在同一离散 ELBO 上优化；后者设计矩阵是截距加三个固定空间 Gaussian RBF，而不是多项式。每例保存 KL(q||exact)、ELBO/log-evidence 恒等式、Brier、熵、拟合及留出相对场误差。场误差相对于带噪 N64 数据，不能解释为 clean-manifold 误差。冻结配置、源码哈希和缓存哈希见 `runs/a2/extensions/probability/v2/frozen_config.json`。

`weak/` 是完成 v2 后新增的弱数据诊断：只取第一频率第一发射的两个复接收通道，30 个新对象和 0.3/0.8/1.5 reference-RMS 噪声。它是 post-hoc 诊断，不取代 v2。`optimization_posthoc/` 以 v2 同种子数据重建后比较 RBF 的零、product 边际投影和两个固定随机初始化；这是分离表示能力与优化的诊断，不是修改主结果或真值 oracle。
