> SUPERSEDED: this is early scaffold documentation. Use protocols/a2/EXTENSIONS_PROTOCOL_ZH.md and manifold_v3 frozen config. Invalid pilot was not accepted.

# 流形切空间开发协议

固定峰值幅度 Gaussian，协方差使用 SPD retraction；协方差导数不含 trace 项。共同目标是训练接收通道的 full-wave 复散射残差，held 通道只作验收。四法为全 SPD-LM、谱 rank continuation、弱方向 secant+可见 LS+全波验收、同候选数随机控制。内部算子只选择探索方向，不作为观测残差。30 例为开发集，不是冻结验收。
