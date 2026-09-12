# A2 实际使用的 Fresnel 数据索引

## 原始文件与来源

- 文件：`Gaussian/data/dielTM_dec8f.exp`
- SHA-256：`476cc9d1cfc98797545ab4adf69302dc5aeb45848a24cf8d7b3d222940cc79eb`
- 来源：Institut Fresnel / IOP 2001 补充实验数据的本地未改动副本；详细出处、获取日期和链接见上级 [README.md](../README.md)。
- 使用边界：该数据按官方说明可用于科学研究；IOP 页面说明权利仍归作者或另有说明者。本仓库不重新授权，也不应重新分发原始文件。
- 获取边界：本工作未绕过 CAPTCHA、注册或访问控制；未获取 Fresnel 3D 数据。3D 下载条款要求额外账户和结果披露许可，故不在 A2 范围内。

## 2001 TM 文件与 A2 取用约定

文件是名义二维、偏心均匀介电圆柱的相干 TM 测量：36 个源视角、每视角 49 个接收样本、1–8 GHz、共 14,112 行。列为 `view receiver frequency_GHz Re(total) Im(total) Re(incident) Im(incident)`；原始论文约定为 `exp(+i omega t)`。

A2 的实测诊断仅取 2、4、6 GHz 和零基视角
`[0,3,6,9,12,15,18,21,24,27,30,33]`。使用既有的共轭转换以对应 A2 的 `exp(-i omega t)` 标量 Green 约定；从背向弧估计每视角 incident/source 因子。固定几何为源半径 0.720 m、接收半径 0.760 m，采用公布的 36/72 角位置；不做几何校准。全局偶数接收索引用于拟合，奇数接收索引保留。

## 增益与解释范围

主拟合还在每个频率的散射训练数据上剖面化一个复增益。它是 incident 因子估计之后的额外 nuisance，不能称作“仅源归一化”，也不能消除材料幅度歧义。固定增益重拟合只是同一预处理下的诊断。

这是一项单一二维圆柱的模型失配检查；它不支持多目标泛化、硬件校准、材料定量恢复、三维验证或神经网络必要性的结论。实际使用和结果口径以 [runs/a2/measured/results.json](../../runs/a2/measured/results.json) 与 [delegated/a2_ports/MEASURED.md](../../delegated/a2_ports/MEASURED.md) 为准。
