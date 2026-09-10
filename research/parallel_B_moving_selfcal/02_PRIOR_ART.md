# 最近邻文献审计：不是按关键词组合宣称原创

访问日期：2026-09-11。F=关键正文/公式已读；A=仅作者/期刊摘要及机构记录，决定性全文未核对；H=所引用书籍或官方技术文档。未读全文的“未显示定理/硬件”不等于断言不存在。下表的8项均按可取得材料报告。

## 八项假设对比

| 文献 | 已知标准/先验 | 静态或离线？ | 几何随采集变化？ | 定量未知材料？ | full-wave？ | joint inversion？ | 可辨识性理论？ | 硬件？/访问 |
|---|---|---|---|---|---|---|---|---|
| R1 Keysight校准注意事项 | 固定夹具及VNA标定 | 校准后应避免电缆移动 | 讨论移动导致误差，不反演几何 | 不研究 | 不适用 | 否 | 否 | 官方说明 H |
| R2 Willi/Guillaume2019 | 光学计量+机器人运动学 | 离线几何标定 | 多工具姿态 | 否 | 不研究散射 | 几何参数拟合 | 摘要未显示 | 实机 A |
| R3 ud Din等2022 | GRL标定、已知样件CAD | 标定后控制扫描 | **移动样件，不是移动天线** | 主任务反射损耗图 | 测量真实EM；非完整材料反演 | 否 | 未给材料唯一性定理 | 实机 F |
| R4 Fikes等2019 | 阵元互耦 | 动态形状自校准 | 阵列弯曲 | 否 | 研究天线耦合，非材料散射反演 | 阵列形状估计 | 摘要未显示 | 8阵元CMOS实机 A |
| R5 Fikes等2021 | 阵元互耦模型 | 形状更新 | 形变阵列 | 否 | 耦合含EM，非目标材料反演 | 形状重建 | 决定性定理未核对 | 多种阵列实机 A |
| R6 Ludvig-Osipov等2025 | 制造阶段S(ε)校准库 | 增益在估计期准静态 | 未知几何不是该模型参数 | 平均复介电常数 | 不简化S的物理响应；数值腔体例 | ε及Rx/Tx增益 | 定义对角等价下可辨识，非本移动域证明 | 数值 F |
| R7 Guo等2022 | **两个均匀校准仿体** | 校准映射后成像 | 非未知移动阵列 | 介电图 | DBIM迭代全场 | 校准+成像两阶段 | 摘要未显示 | 仿体及志愿者 A |
| R8 Repetti等2017 | 图像稀疏、平滑DDE；正文用最亮源先验 | 自校准 | 已知干涉几何；未知方向增益 | 天空强度，不是ε | 傅里叶/RIME，不是多次散射材料模型 | 图像及DDE | 非凸算法收敛，不是所有世界全局唯一 | 数值 F |
| R9 Sob等2021 | 天空模型、增益时间/频率结构 | 求解时间间隔/在线变化 | 增益变化，不联合未知接收位置 | 否 | RIME | 校准/成像流程 | 校准间隔选择，不是Maxwell唯一性 | MeerKAT真数据 F |
| R10 Li等2021 | 已建模波形/路径 | 同步误差分析 | 分布式链路时延等 | 否 | 波形/检测模型 | 检测，不是ε及动态几何联合 | 检测概率解析，非本类辨识定理 | 摘要未显示硬件 A |
| R11 Cheng等2024 | POS与场景高度/点目标辅助 | 运动补偿 | UAV轨迹 | 反射率/SAR，不是ε | 传播相位模型 | 运动补偿成像 | 非本类材料定理 | 实测SAR A |
| R12 Huang等2017 | 源一致性正则 | 迭代FWI | source/receiver-extension，不是任意自由相位许可 | 地震速度，不是ε | 声/地震全波 | 速度及扩展源 | 简单运动学下与旅行时误差联系；全文未核 | 数值 A |
| R13 Takenaka/Moriyama2012 | 已知测量边界的总场 | 无需先验入射场 | 不联合未知边界坐标 | 介质/缺陷 | 场等效原理 | 消去未知入射条件 | 摘要未显示全域辨识定理 | 数值 A |
| R14 Zhang等2023 | 已知测量边界的Cauchy数据 | 源与障碍共同未知 | 未知源位置；接收边界已知 | 障碍形状，不是介电图 | **声波** | 是，分解为两个子问题 | 分解误差及指标渐近 | 数值 A |
| R15 Ji/Liu2019 | 可加入已知磁偶极参考源 | 无相位稀疏远场 | 未知辐射源位置/分布 | 不同时反演一般介质 | 电磁源模型 | 源反演/相位恢复 | 有唯一性类结果，非本移动材料域 | 数值 A |
| R16 Yuyao Chen/Dal Negro2021 | 已知照射/近场数据及PDE | PINN近场参数恢复 | 不联合未知移动校准 | 复ε、µ | 矢量Maxwell PINN | 材料PDE反演 | 收敛实验不等于信息创造 | 摘要未显示实际硬件 A |
| R0 Xudong Chen2018 | 各章节明确观测与源假设 | 多类方法 | 非本动态校准主课题 | 定量ε等 | 含VIE/CSI/SOM等 | 状态与材料一致性 | 各章节不同条件 | 含数值/部分实验 H |

## 关键公式真正重叠在哪里？

R6式(4)为D=RST；R,T分别为对角收发增益，端口匹配、忽略串扰、估计期增益准静态。§2.4.2定义：S(ε1)=R̃S(ε2)T̃只能在ε1=ε2时成立，才称材料可辨识。式(9)–(10)对ε外层搜索、对补偿增益内层拟合。其目标是校准后的S误差，不直接等同于本轮固定噪声原始数据GLS。该文已明确动机包括漂移和运行时无法插入参考。因此“材料和增益联合估计”“in-situ不可插参考”都不能认领为新思想。差别是本轮未知帧几何和delay，但差别本身不是胜利。

R8的y=GFx+b、Y=CXC*+B与双边增益和未知图像联合建模有关，采用块坐标优化及先验。算法收敛保证不等于物理可辨识。R9式(5) Vpq=Gp Cpq Gq^H+噪声+RFI，讨论求解间隔造成的漂移/噪声/模型不全折中。于是时间平滑或短窗校准也不是本轮原创。

R7的两个参考仿体不仅估计电子学，还帮助匹配实际系统和前向模型。本轮只实现公共gain/位置标准模型，未声称逐式复现FD-DBIM；它提醒我们最强传统校准可能比“固定增益一次相除”更强。

R4/R5已有无需外置目标的互耦形状校准。柔性阵列的形状未知不能自动成为新材料自校准的gap。R3证明机械自由度与精确扫描控制可以维持校准条件；必须准确写成机器人移动样件而不是天线。

## 本轮真正尚未闭合的差距

缺的是实际移动头在同等信息/硬件/时间预算下，具有“传统校准无法经济跟踪”的实证条件，以及其剩余材料/校准歧义被可测全波相对结构打破的非局部证明。本文没有以“未找到完全相同标题”代替这些要求。

联合未知Tx位置与定量未知介质的直接电磁最近邻仍未完成全文排除。R14是声学、R15是源反演，不能拿它们替代该检索缺口。此缺口阻止首次性声明，不阻止本轮按实际失败给出NO。

## 可追溯来源

[R0] Xudong Chen, Computational Methods for Electromagnetic Inverse Scattering, Wiley/IEEE, 2018. 用户提供，尤其第2章Green/VIE/偶极、第6章状态—材料一致性和Twofold、第8章无相位。未在交付中复制整书。

[R1] Keysight, Opt.001 Calibration Considerations. https://helpfiles.keysight.com/csg/N1500A/Opt.001_Calibration_Considerations.htm

[R2] D. Willi, S. Guillaume, 2019, Calibration of a Six-Axis Robot for GNSS Antenna Phase Center Estimation. DOI:10.1061/(ASCE)SU.1943-5428.0000291.
https://ascelibrary.org/doi/10.1061/%28ASCE%29SU.1943-5428.0000291

[R3] Salah ud Din et al.,2022, Robotic scanning free-space measurement system for electromagnetic performance evaluation of curved radar absorbing structures. DOI:10.1177/17298806221114554.
https://journals.sagepub.com/doi/10.1177/17298806221114554

[R4] A. C. Fikes et al.,2019, Flexible, Conformal Phased Arrays with Dynamic Array Shape Self-Calibration. DOI:10.1109/MWSYM.2019.8701107.
https://authors.library.caltech.edu/records/27wy6-xp082

[R5] A. C. Fikes, A. Mizrahi, A. Hajimiri,2021, A Framework for Array Shape Reconstruction Through Mutual Coupling. DOI:10.1109/TMTT.2021.3097729.
https://authors.library.caltech.edu/records/cfedq-42j80

[R6] A. Ludvig-Osipov et al.,2025, Auto-Calibration of Near-Field Microwave Measurements for Complex Permittivity Estimation, PIER M136,46–56. DOI:10.2528/PIERM25090302.
https://www.jpier.org/ac_api/download.php?id=25090302

[R7] L. Guo et al.,IEEE TMI41,1087–1103,2022. Calibrated Frequency-Division Distorted Born Iterative Tomography for Real-Life Head Imaging. DOI:10.1109/TMI.2021.3132000.
https://ieeexplore.ieee.org/document/9632616/

[R8] A. Repetti et al.,2017, Non-convex optimization for self-calibration of direction-dependent effects in radio interferometric imaging. DOI:10.1093/mnras/stx1267.
https://arxiv.org/abs/1701.03689

[R9] U. M. Sob et al.,2021, Solution intervals considered harmful: on the optimality of radio interferometric gain solutions. DOI:10.1093/mnras/stab928.
https://academic.oup.com/mnras/article/504/2/1714/6211002

[R10] H. Li et al.,2021, Signal Detection in Distributed MIMO Radar with Non-Orthogonal Waveforms and Sync Errors. DOI:10.1109/TSP.2021.3087897.
https://arxiv.org/abs/2102.09719

[R11] Cheng, Qiu, Meng,2024, Precise Motion Compensation for High-Resolution UAV SAR Imaging Based on Improved PTA. DOI:10.3390/rs16142678.
https://www.mdpi.com/2072-4292/16/14/2678

[R12] G. Huang, R. Nammour, W. Symes,2017, Full-waveform inversion via source-receiver extension. DOI:10.1190/geo2016-0301.1.
https://doi.org/10.1190/geo2016-0301.1

[R13] T. Takenaka,T. Moriyama,2012, Inverse scattering approach based on the field equivalence principle: inversion without a priori information on incident fields. DOI:10.1364/OL.37.003432.

[R14] D. Zhang,Y. Chang,Y. Guo,2023, Jointly determining the point sources and obstacle from Cauchy data.
https://arxiv.org/abs/2306.06665

[R15] X. Ji,X. Liu,2019preprint, Inverse electromagnetic source scattering problems with multi-frequency sparse phased and phaseless far field data. DOI:10.1137/19M1256518.
https://arxiv.org/abs/1906.02187

[R16] Y. Chen,L. Dal Negro,2021, Physics-informed neural networks for imaging and parameter retrieval of photonic nanostructures from near-field data. DOI:10.1063/5.0072969.
https://arxiv.org/abs/2109.12754

## 额外近邻线索，不冒充已阅读全文

10.1002/2017RS006399、10.1002/2017RS006403，以及 Sensors26(11)3517(2026)的天线建模与DBIM工作，需要进一步核对决定性全文。未将其摘要替代完整假设或宣称复现。
