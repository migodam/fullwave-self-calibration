# V2 执行细则（生成新场景前冻结）

Treams 与 python-flint 在当前执行环境未安装，外网 pip 解析失败。本轮采用 PROTOCOL 已允许的 SciPy 向量球波／角向投影求解器；不把替代实现称为算法贡献。电／磁球模边界系数、curl 恒等式、平面波重建、小球极限和独立角积分检查均须通过。主模型 lmax=3，Gauss–Legendre 极角14点、方位角28点；高阶检查固定为 lmax=4/5/6。

进一步固定原协议尚未指定的机械细节：附加 EM 候选池为16个 Fibonacci 接收方向×3个笛卡尔分量×4个照射，按 receiver/component/illumination 的 C-order 展平。fixed 方法总取索引0，random 方法用独立 RNG 2026091108 在0…191中均匀抽样；场景和全部噪声仍用 RNG 2026091107。四方法均用相同的基础观测；附加 EM 观测噪声独立、绝对 sigma 与基础观测相同。noisy-reference GLS 使用一个复电子学参考，复 sigma=0.01。不同硬件路径不宣称成本完全相等。

原12个新场景之外的额外计算仅作诊断：前4个场景固定细化 DDA 到 spacing=0.009米、fill_quadrature=4；同一粗网格数据上诊断 exact-gain/oracle-geometry 只用于归因，不进入基线排名。不使用这些诊断选择天线或调整失败选择器。

区间证明使用 Python Fraction 精确有理运算、带严格余项的球贝塞尔级数，材料盒宽1/1000。通过解析消去 sqrt(epsilon) 避免尺寸参数的近似输入；全部输出定向取整。此方法不依赖工作报告的错误导数，也不声称极值必在端点。
