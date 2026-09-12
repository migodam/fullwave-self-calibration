# GPU 表示能力 campaign 协议

此 campaign 比较普通正 Gaussian 表示（K=16、64、144）和 N128 voxel 表示；它不包含 SOM，因此不能作为 SOM 贡献的证据。每个表示的实对比度均乘已知损耗比 `(1+0.15j)`；Gaussian 使用正幅度、受限中心和 Cholesky 参数化的 SPD 协方差，长度单位为米。初始 Gaussian 中心是固定均匀网格，不从真值取得；voxel 从同一 K=64 初始 Gaussian 图开始。

生产配置包含 30 个固定 seed（`20260915+case`），每族 10 个：非网格 16--32 Gaussian 混合、尖边多夹杂、细长曲线并带空洞。数据由 N192 `cell_integrated=True` 的三频（1.25/2/2.75 GHz）、12 Tx、64 半孔径 Rx 完整 VIE 生成；反演为 N128。偶数 Rx 训练，奇数 Rx 留出。每方法 Adam 80 步；每步三频全 12 Tx 前向和隐式伴随。两种表示使用同一归一化场损失和 `1e-4` 的离散相邻像素平方差平滑项，不使用额外 TV。

每 case/method 原子 checkpoint，并拒绝不同配置哈希的恢复；保存真值、重建、轨迹、material L2、带 clean held-field 的相对误差、IoU、实际前向/伴随 RHS 计数、线性迭代、耗时、CUDA 峰值内存。配置与 source hash 在运行前冻结。建议先运行 `--n 32 --cases 0 --steps 2` 验证 GPU 接口，随后才启动默认 120 fits。正则、学习率与步数只是预先声明的发展设置，不能单独证明公平充分或表示优越性。
