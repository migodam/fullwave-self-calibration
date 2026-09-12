# A2 概率成像 v2 结果摘要

主 v2 完成 30/30 个 N64 cell-integrated 生成实例，N64 前向累计 1.97 s，端到端 41.83 s；复用的 N32 256-state 缓存原始构建为 10.48 s。ELBO 与独立 `KL(q||posterior)` 的最大恒等式残差为 `1.82e-12`。所有 product/Gaussian 优化均报告成功；精确枚举没有被计为零成本，而是未单独计时。

在 v2 的带噪奇数接收机留出通道上，噪声 fraction 为 0.03/0.1/0.3 时，exact 的平均相对场误差为 0.0557/0.2571/0.5194；product VI 为 0.0557/0.2571/0.5194，平均 KL 为 `5.53e-5/2.01e-5/4.83e-5`。固定三-RBF Gaussian-logit 的对应留出误差为 0.1851/0.3810/0.5830，KL 为 2925.60/341.06/29.89。这是有限字典、单一固定几何和 N32 近似模型下的结果，不支持泛化或算法优势断言。

弱数据 post-hoc 诊断的 noise fraction 为 0.3/0.8/1.5；对应 exact 平均 patch 熵为 0.369/0.564/0.641，表明该设置不再全部后验饱和。product 的 KL 为 0.840/0.837/0.528；Gaussian-RBF 的 KL 为 3.747/1.520/0.901，且 fraction 0.3 的 Brier 为 0.423。这显示在该受限 RBF 参数化中有模式与因子化误差，不能归因为传统 VI 的一般失败。

v2 的 RBF 初始化 post-hoc 诊断保存了同种子重建的 30 例 clean/noisy fields。选择最高 ELBO 后，v2 noise fraction 0.03/0.1/0.3 的 Gaussian-RBF KL 仍为 1891.71/277.85/16.22；product-marginal 投影初始化分别被选中 6/2/4 次。这是对表示/优化分离的负向诊断，不改变已冻结 v2 结果。

重建图采用统一色标与 cm 坐标：`figures/a2/extensions/probability_v2_reconstructions.png`。详细逐例结果在 `runs/a2/extensions/probability/v2/results.json`、弱数据在 `runs/a2/extensions/probability/weak/results.json`。
