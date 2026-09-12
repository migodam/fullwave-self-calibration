# 有限 patch full-wave Bayes 实验摘要

本实验冻结 8 个已知 Gaussian-centre Voronoi patch（ROI 内）、二元已知材料、ring Ising prior、N32 开发 VIE、1.5/2.75 GHz、4 Tx、半孔径 24 Rx。256 个 label realization 均分别解 full-wave；因此 likelihood 使用缓存的 `F(z)`，没有把 `F(E chi)` 当作 `E F(chi)`。

为避免零散射空对象没有定义的相对 SNR，prior-predictive 抽样条件化为非空 label state；这一条件写入 frozen config。30 个对象分为高/中/低噪声各 10 个。N64 不同网格的三组同物理 state 散射相对差为约 2.55%、2.20%、4.47%，仅是网格敏感性读数。

精确 256-state posterior 是本实验基准。Product Bernoulli VI 和固定四维 Gaussian-logit VI 都对枚举 likelihood 精确优化 ELBO；结果内的 `kl_to_exact` 按 `log evidence - ELBO` 保存，精确 posterior 为零（浮点舍入范围内），未将 ELBO 当作任意下界之外的东西。

高 SNR 时 product VI 与精确 posterior 接近；固定低维 Gaussian-logit 在中/低 SNR 的有限实现中出现较大的 KL、patch Brier 和 posterior-predictive 误差。这是该固定基函数族在本冻结有限任务上的结果，不是概率算法、Gaussian 参数化或 Bayesian inverse scattering 的普适比较。

可靠性图是 30 个 prior-predictive draw 的分箱描述，不能解释成某个固定真值重复加噪声的 calibration。结果不主张泛化、原创性、TAP 优势或真实材料不确定性验收。
