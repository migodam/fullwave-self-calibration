# TriSpace 源文件分节地图（source map）

图例：A = 主报告；B = 研究方向纠偏附件；C = 目标重定位附件；行号均为各文件自身行号。

## 1. Source A（4,748 行）

### 封面与第 0 章

| 内容 | 行范围 |
|---|---|
封面/版本/对象（二维标量 Helmholtz contrast-source） | 1–9 |
§0 证据边界、符号与不可越界语义 | 10–63 |

### §1 一页结论摘要

| 内容 | 行范围 |
|---|---|
§1.1 可作为主结论的命题（13 条） | 66–200 |
§1.2 被反例否定或必须降级的命题 | 201–212 |
§1.3 仍需实质工作的问题 | 213–223 |

### §2–§3 状态与假设

| 内容 | 行范围 |
|---|---|
§2 Feasibility triage 与 P1–P10 状态表 | 224–240 |
§3 假设总表 H1–H19 | 241–267 |

### §4 定理—证明—反例正文

| 包 | 内容 | 行范围 |
|---|---|---|
4.1 | P1 Hilbert-space retention、polar operator、闭包补偿、Friedrichs/稳定横截、紧凑反例、离散收敛 | 270–686 |
4.2 | P2 连续 Helmholtz SE(d) gauge、Lie tangent、stabilizer、gauge 破缺、quotient FIM | 687–892 |
4.3 | P3 Born 远场 Fourier 几何、Ewald/missing wedge、有限采样饱和、frame/coherence 正界、近混淆 | 893–1198 |
4.4 | P4 多频共享位姿补偿、精确维数公式、条件 genericity、单调/反例、submodularity 反例 | 1199–1506 |
4.5 | P5 半正定位姿先验、Schur range condition、短接/增广投影、motion-graph gauge、weighted retention | 1507–1907 |
4.6 | P6 固定秩流形、投影导数、generalized retention 导数、重根分裂、Kato 界、rank-event 分类与 detector | 1908–2265 |
4.7 | P7 full-wave resolvent 界、Hankel 外部界、D^2B=D_X^3F 缺口、内部奇异性、无 uniform constant 的不可能性 | 2266–2695 |
4.8 | P8 robust trajectory 一阶/二阶展开、约束锥、重根 worst problem、rank-event trust region、Pareto 原则 | 2696–2930 |
4.9 | P9 SOM/current subspace 与 map-tangent/pose-defect 的合法桥接、range/等号、奇异值界、不可比语义、不交换反例、two-stage 算法 | 2931–3318 |
4.10 | P10 局部谱→估计误差：三种 protocol、sandwich、Bayesian marginal、LAN/local basin、全局反例、安全表述 | 3319–3589 |

### §5 可直接替换的论文文本（LaTeX 建议）

| 小节 | 内容 | 行范围 |
|---|---|---|
5.1 | Hilbert retention 主定理 | 3594–3661 |
5.2 | 有限网格最小角不能作为 continuum certificate 的反例 | 3662–3683 |
5.3 | rigid-motion gauge 与 quotient information | 3684–3726 |
5.4 | 半正定 pose prior 的正确 Schur 对象 | 3727–3762 |
5.5 | 多频 shared compensation | 3763–3803 |
5.6 | fixed-rank 光滑性与 rank-event 不连续 | 3804–3828 |
5.7 | SOM/current-space 与 map-tangent 的合法桥接 | 3829–3867 |
5.8 | 统计范围 limitation 段落 | 3868–3885 |
5.9 | 建议替换 manuscript 的 limitations | 3886–3926 |
5.10 | 附录：finite-dim generalized Schur condition | 3927–3965 |
5.11 | 附录：local robust eigenvalue remainder | 3966–3993 |

### §6 有数学依据的设计与实验修改

| 小节 | 内容 | 行范围 |
|---|---|---|
6.1 | 先限制到可证 task window | 3998–4024 |
6.2 | Gauge 须显式 quotient/anchor | 4025–4053 |
6.3 | 合法 SOM–pose two-fold reduction（算法 6.1） | 4054–4107 |
6.4 | 多频选择与共享补偿矛盾 | 4108–4156 |
6.5 | 半正定 motion prior：sparse joint solve | 4157–4182 |
6.6 | Rank-event detector：阈值、hysteresis | 4183–4208 |
6.7 | Robust trajectory：margin-limited trust region | 4209–4253 |
6.8 | 轨迹评估为 Pareto 而非预设排序 | 4254–4272 |
6.9 | 三个统计 protocol（F/B/A） | 4273–4295 |
6.10 | 连续/离散验证最低要求 | 4296–4309 |

### §7 文献 claim ledger

| 小节 | 内容 | 行范围 |
|---|---|---|
7.1 | 检索边界与证据强度 A/B/C | 4312–4321 |
7.2 | 逐命题 ledger（C1–C16） | 4322–4342 |
7.3 | 关键文献目录与覆盖对象 | 4343–4385 |
7.4 | Retrieval-bounded novelty 结论 | 4386–4401 |

### §8 尚未解决清单

| 小节 | 内容 | 行范围 |
|---|---|---|
8.1 | P1：真实 Helmholtz 离散化 retention 谱收敛 | 4404–4427 |
8.2 | P2：含 direct path、phase center、标定参数的完整 gauge group | 4428–4443 |
8.3 | P3：line/circle/arc 的 scene-independent exact intersection | 4444–4465 |
8.4 | P4：具体阵列 analytic genericity witness | 4466–4481 |
8.5 | P4/P8：exact shared-Schur weak submodularity | 4482–4501 |
8.6 | P5：无限维 semidefinite motion prior 的 closed-range | 4502–4522 |
8.7 | P6：有噪 Jacobian 下 rank-event detector 概率保证 | 4523–4538 |
8.8 | P7：moving-grid D^2B=D_X^3F continuum bound | 4539–4548 |
8.9 | P7：verified full nonlinear Lipschitz certificate | 4549–4562 |
8.10 | P8：跨 rank strata 的 nonsmooth robust trajectory solver | 4563–4584 |
8.11 | P9：TSOM G_D split 与 pose projection 的统一误差 | 4585–4606 |
8.12 | P10：local basin→observable noise threshold | 4607–4622 |
8.13 | 文献：投稿前系统查重 | 4623–4643 |

### §9–§10 审计与最终状态

| 小节 | 内容 | 行范围 |
|---|---|---|
§9 | 反过度主张审计表（25+ 类禁止表述） | 4644–4678 |
9.1 | 最终语义一致性检查 | 4679–4692 |
9.2 | 最终可主张的论文中心句 | 4693–4700 |
§10 | 最终状态与使用指南 | 4701–4748 |
10.1 | P1–P10 最终标签 | 4703–4717 |
10.2 | 推荐首篇理论稿四条主线 | 4718–4728 |
10.3 | 本报告证据层级 | 4729–4736 |
10.4 | 引用与复用注意 | 4737–4748 |

## 2. Source B（349 行）

| 内容 | 行范围 |
|---|---|
方向纠偏：回到 TriSpace SOM-SLAM；A/B/K_eff 降为辅助层；经典 full-wave contrast-source 与 SOM/Two-fold SOM 回顾；X 未知⇒G_S,G_D,E_inc 未知⇒V_S^±,V_D^± pose-dependent | 1–72 |
“TriSpace 的第三个空间”：Range(D_XF) 在 data space，V_S/V_D 在 current space；V_P 必须为 current-space pose-induced/equivalent-current subspace；不得人为构造三空间交集 | 73–110 |
本项目需要回答的理论问题：pose-dependent SOM decomposition（113–142）；subspace leakage P_S^-(X*)P_S^+(X̂)（143–172）；TriSpace observability/classification（173–184）；rank event 与子空间连续性（185–206）；joint full-wave SOM-SLAM formulation 与先消去 current DOF 再联合恢复 (χ,c,X)（207–244） | 111–244 |
A/B/K_eff 理论定位：local identifiability / Fisher 一致性检查，不得反向解释成 SOM | 245–311 |
最终定位：Pose-Uncertain Full-Wave TriSpace SOM for Inverse-Scattering SLAM；优先解决“Pose space 如何 pull back 到 SOM current space”，不得令 V_P=Range(D_XF) | 312–349 |

## 3. Source C（264 行）

| 内容 | 行范围 |
|---|---|
目标重定位：self-calibrating full-wave inverse scattering；模型 y=F(χ,X)、G_S(X)、D_χ(E_inc(X)+G_D J)；工程动机 5 点；self-calibration 不是 novelty；真正新交叉点 = unknown geometry + nonlinear volumetric full-wave + SOM/TSOM reduction | 1–46 |
为什么用 SOM：joint optimization 高维强非线性；SOM 在优化前对 induced current 做物理降维；unknown X 使 G_S(X) 的谱坐标本身未知 | 47–87 |
TriSpace 真正要解决的问题：每 mode 回答 Sensing/State/Pose 三个物理问题；目标是 J_safe+J_pose-confounded+J_state-recoverable+J_unresolved 分类；不预设正交直和、不制造 2^3 | 88–130 |
理论研究目标：A. pose-dependent SOM subspace theory（135–150）；B. calibration-induced current ambiguity / 合法 pullback（151–162）；C. self-calibrating reduced-order inversion（163–182） | 131–182 |
贡献链：known-geometry SOM→pose-dependent SOM→subspace perturbation/leakage/rank transition→TriSpace 分类→reduced-order self-calibration→与 Born/linearized 等方法比较；实验关注 Born mismatch 的 multiple-scattering 区 | 183–214 |
对此前 A/B/K_eff 的定位：辅助 identifiability 层；不得用 A/K_eff eigenvectors 代替 SOM modes | 215–239 |
最终定位与研究问题英文/中文表述；TriSpace 只有给出 identifiability criterion、subspace theorem、降维或重建 advantage 才构成贡献 | 240–264 |

## 4. Theory/SOM_SLAM_THEORY_CONTEXT.md（按需读取部分）

| 小节 | 内容 | 行范围 |
|---|---|---|
0.1–1.2 | 状态标签、三模型防混淆、从 Chen SOM 到 A/B/K_eff 主线 | 7–135 |
5.4–5.5 | V_D^± 的术语边界与 hard/soft subspace | 463–519 |
9.1–9.6 | G_S/G_D 两折仍在 current space；TSOM 骨架；第三层 pose-fold；map pullback | 992–1238 |
21–22.3 | 冲突/修正总表与最终核心研究问题 | 3028–3141 |
