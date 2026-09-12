# 有限类别 full-wave Bayes 成像协议

冻结八个已知 Gaussian-centre Voronoi patch、ROI、两种已知材料、ring Ising prior、频率、半孔径与噪声分组。每个 8-bit label realization 都独立解 full-wave VIE；likelihood 使用缓存的 realization field，绝不使用 `F(E chi)` 代替 `E F(chi)`。

比较精确 256-state posterior、product-Bernoulli mean-field VI 与固定低维 Gaussian-logit VI。ELBO 对枚举 likelihood 精确计算，并用 `log evidence - ELBO = KL(q||posterior)` 检查。30 个对象从同一 prior predictive 分布抽样；reliability 只按该 prior-predictive 语义解释。
