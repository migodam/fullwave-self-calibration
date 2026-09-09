# TriSpace 相关主张清单（claim inventory）

## 标签说明

- **prior-known**：源文件明确承认或属于既有数学/文献组件。
- **user hypothesis**：定位、提议、研究路线判断；尚未成为定理。
- **derived-but-unverified**：源文件内有推导/证明（“讨论内部成立”），但尚无独立数值或外部验证；A 称为“已推导/已接受”或“需补条件”的也归入此类并注明条件。
- **numerically supported**：本批三个源文件没有任何独立数值实验支持；仅内部精确 worked counterexample。除特别说明外不授予此标签。
- **contradicted**：被同一/后续源文件反例或审计否定（或“必须降级”）。
- **open**：源文件自己列为未解决。

部分主张为复合状态（如 derived + open），主标签外括号注明。

## 1. Source A：P1–P10 的主要定理/算法主张

| ID | 主张 | 位置 | 状态 | 证据/备注 |
|---|---|---|---|---|
| A1 | Continuum retention 可由 polar decomposition 定义 R_geom=U_A^*P_{N^⊥}U_A，且 K_SLAM=\|A\|R_geom\|A\|，无需无界 K_IS^{-1/2} | 4.1（270–372） | derived-but-unverified | 抽象 subspace/angle 理论本身 prior-known（Bouldin/Björck–Golub，§7.2 C1）；模型专门化是新组合 |
| A2 | 无先验消元看到的是 closure \overline{Ran B}；“精确补偿”与“渐近补偿”须区分；非闭 range 时 K=0 只意味着任意好逼近 | 4.1 推论 4.2/命题 4.3（373–456） | derived-but-unverified | 无限维反例/泛函分析结论；论文替换文本 5.2 |
| A3 | “交集为零”远弱于“稳定横截”；γ>0 ⇔ M+N closed；compact-A 反例显示每层网格正角但极限 γ=0 | 定理 4.4/反例 4.5/推论 4.6（456–629） | derived-but-unverified（含内部反例） | **contradicted**：否定“M∩N={0}⇒稳定”和“有限网格 θ_min>0 是 continuum certificate”（§1.2、§9） |
| A4 | 有限维离散逼近只在固定 positive-information task window + operator/gap convergence 下收敛；网格正角不能外推 continuum | 定理 4.7（630–686） | derived-but-unverified / open（真实 Helmholtz 离散化见 8.1） | A 自己的限制说明 |
| A5 | 均匀标量背景、body-fixed 传感器下 global rigid motion 是精确 gauge：F(g·χ,g·X)=F；Lie tangent 满足 Aξ_χ+Bξ_X=0 | 4.2（687–892） | derived-but-unverified（conditional） | 连续 scalar homogeneous 模型；vector/Maxwell、非均匀背景、direct path/clock/phase center 未闭合（8.2） |
| A6 | Radial-symmetric scene 的 rotation generator 是 pose-only stabilizer；gauge 维数不恒为 3；离散 pixel 可产生 representation residual | 4.2（810–880） | derived-but-unverified | 修正“2D 恒有 3 个 map gauge modes”（§9） |
| A7 | Quotient FIM 用 horizontal compression；stabilizer 时按 orbit-type stratum | 定理 4.13（843–880） | derived-but-unverified | 群作用/信息几何 prior-known |
| A8 | Born 远场下 A、B 的显式 Fourier 关系；global translation/rotation 是该补偿方程显式解；其余解 scene-dependent | 4.3（893–1042） | derived-but-unverified | H7/H8 条件；finite sampling 与 continuum 分开 |
| A9 | 有限采样+无限维自由 map 通常 Ran B⊆Ran A（full-row saturation）；连续 aperture 下 restricted Fourier range 稠密，故无 geometry-only continuum 正主角下界 | 定理 4.19–4.21（1043–1133） | derived-but-unverified | **contradicted**：否定“aperture/MIMO/standoff 单独给出 continuum 正最小角” |
| A10 | 有限维 frame/coherence 条件给出 γ≥√(1−(μ/ab)²)；缺失 wedge 首先造成 absolute loss 而非自动 map–pose intersection | 4.22–4.23（1134–1198） | derived-but-unverified | 几何名称本身不是证明 |
| A11 | 共享补偿核维数公式；条件版解析频率 genericity：存在非退化 witness 时几乎处处只留 gauge | 4.24–4.25（1201–1322） | derived-but-unverified（conditional） | 具体阵列 witness **open**（8.4）；反例族 4.27 |
| A12 | 新频块使 absolute K_eff Loewner 单调；normalized retention 不单调；一维反例 | 4.28–4.29（1347–1413） | derived-but-unverified（内部精确反例） | **contradicted**：否定“带宽/频率越多 normalized retention 越高” |
| A13 | 精确 shared-Schur frequency objective 一般不 submodular；两 block 反例；可用逐频独立消元构造 monotone submodular lower-bound surrogate | 4.30–4.32（1414–1506） | derived-but-unverified（内部精确反例） | **contradicted**：否定“shared-pose logdet 必 submodular” |
| A14 | 半正定 pose prior 的无限维正确对象是 augmented-space projection/shorting；Moore–Penrose 形式只在 range/closed-range 条件下有效 | 4.34–4.35（1563–1642） | derived-but-unverified | shorted operator 理论 prior-known（Anderson–Trapp，§7.2 C5） |
| A15 | Finite prior retention 是 prior-weighted shrinkage；只在增广空间才是 squared-sine principal angles；原数据空间一般不是 projector | 定理 4.40（1815–1869） | derived-but-unverified | **contradicted**：否定“finite prior 仍是普通 sin²θ” |
| A16 | Motion-graph global gauge 在 ker J_X；只加强 relative prior 无法逼近 known-pose limit，需 anchor/绝对 reference | 推论 4.37（1714–1717） | derived-but-unverified | 结合 P2 gauge |
| A17 | Fixed-rank stratum 上 P_B 光滑且 ‖Ṗ_B‖≤‖Ḃ‖/σ_r；rank change 处不同秩 projector 范数距离=1；J_X≥αI 使 singularity 正则化 | 4.42/4.47/4.48（1910–2208） | derived-but-unverified | pseudoinverse/projector 导数 prior-known（Golub–Pereyra/Kato/Davis–Kahan）；**contradicted**：否定“rank-change 仍 Hölder 连续” |
| A18 | Online rank-event detector 的 hysteresis 规则与可证保证完全依赖真实 error bound；无上界时只是 heuristic | 算法 4.49（2215–2265） | derived-but-unverified | 概率保证 **open**（8.7） |
| A19 | D^2B 实为 D_X^3F；给出 resolvent inverse、j/E/A/B 的递归常数；二维外部 Hankel 一/二阶统一界成立；内部 self-cell 不能复制外部 bound | 4.50–4.56（2268–2548） | derived-but-unverified（conditional） | moving-grid/self-interaction continuum bound **open**（8.8） |
| A20 | 严格正 prior 下 DK_eff/D^2K_eff 的显式结构常数；近共振信息变大与导数敏感度竞争；无 uniform margin 时没有 continuum constant | 4.57–4.60（2549–2695） | derived-but-unverified | 不可把 sampled scan 当 certificate |
| A21 | Robust trajectory：simple eigenvalue 的 first/second-order 展开；约束锥 worst direction；重根 cluster 问题一般是 concave/nonconvex；rank-event 安全 trust region | 4.61–4.64（2698–2929） | derived-but-unverified（local） | global nonsmooth/rank-crossing solver **open**（8.10）；eigenvalue perturbation prior-known |
| A22 | 轨迹评估应 Pareto；不得预设“圆优于直线/带宽越宽越好” | 4.65/6.8 | derived-but-unverified（design principle）/ user hypothesis | 支持其为防过度主张原则 |
| A23 | A=G_S T_χ ⇒ Ran A⊆Ran G_S；等号 ⇔ H_J=Ran T_χ+ker G_S；G_S 与 A 的 right singular spaces 无内在 Grassmann distance | 4.66–4.69（2949–3056） | derived-but-unverified | **contradicted**：否定“Chen V_S^+ 就是 map observable eigenmodes” |
| A24 | 固定 current reduction 与 pose elimination 只有共同 reducing / 交换条件下才可交换；有显式 2×2 不交换反例；projector product 不是 intersection projector | 4.72–4.74/反例 4.73（3115–3221） | derived-but-unverified（内部反例） | **contradicted**：否定“SOM 与 pose projection 两序总可交换”及 projector-product/八格精确分解 |
| A25 | 合法算法：SOM–pose two-stage reduction（G_S top-r basis → reduced A_r → rank-revealing QR/CS → pose Schur），给出截断误差证书 | 算法 4.75（3222–3318） | derived-but-unverified（algorithm） | 复杂度与失败触发条件明确；不得冒称“Chen TSOM”而未实现 G_D split |
| A26 | 统计语义三 protocol 分开：fixed-nuisance penalty 是 sandwich covariance；Bayesian random-pose marginal 或真实 auxiliary measurement 才给 K_eff^{-1} | 4.76–4.79（3321–3470） | derived-but-unverified | sandwich/M-estimation prior-known（White/Jennrich，§7.2 C16）；**contradicted**：否定无协议声明 covariance=K_eff^{-1} |
| A27 | Local LAN/Gauss–Newton basin 定理给相对 covariance inflation；不排除 cycle-skipping/别名/局部极小（F(x)=e^{ikx} 等反例） | 4.80–4.83（3471–3589） | derived-but-unverified（含内部反例） | **contradicted**：否定“局部 Fisher 好 ⇒ 全局恢复” |

## 2. Source A：论文安全表述、文献边界与审计

| ID | 主张 | 位置 | 状态 | 备注 |
|---|---|---|---|---|
| A28 | §5 的 LaTeX 替换文本只能承载上述已证结果，不能新增承诺 | §5（3590–3993） | derived-but-unverified | 引用时须统一 real/complex convention 并保留 range/gap 条件 |
| A29 | Claim ledger C1–C16：抽象工具不新；模型专门化是窄贡献；Karthik–Ghosh 等前例否定“首次联合反演与定位” | §7（4310–4401） | prior-known / derived-but-unverified | retrieval-bounded；systematic 查重 open（8.13） |
| A30 | P1–P10 最终状态：多数 local/algebraic 问题“solved（conditional）”，具体物理专门化与全局恢复未完成 | §10.1（4703–4717） | open（相应部分） | 见表 |
| A31 | 禁止表述清单：finite-dim 当 continuum、grid angle 当 certificate、machine rank 当 physical rank、local 当 global、Fisher 当 MSE、sampled 当 certified、current-space 当 map-space、八格分解、complex pose projector、inner/outer 重复计数等 | §1.2、§9（201–223、4644–4700） | derived-but-unverified（审计规则） | 论文写作应遵守 |

## 3. Source B：TriSpace 纠偏主张

| ID | 主张 | 位置 | 状态 | 备注 |
|---|---|---|---|---|
| B1 | 论文主线应是 TriSpace SOM-SLAM（经典 full-wave SOM 扩展到 Tx/Rx pose 未知/不确定），A/B/K_eff/Fisher 只作辅助解释层 | 1–72、245–311 | user hypothesis | 与 A 的 P1–P10 主体及 C 的 self-calibration 表述存在侧重差异 |
| B2 | X 未知 ⇒ G_S(X)、G_D(X)、E_inc(X) 未知 ⇒ SOM 的 V_S^±、V_D^± 本身 pose-dependent；谱分解本身成为未知量 | 50–72 | user hypothesis | 与 A/B/K_eff 固定 nominal Jacobian 语义不冲突但不同层 |
| B3 | TriSpace 第三空间 V_P 必须是 current space 中由 pose 诱导的 equivalent-current subspace；Range(D_XF) 在 data space，不能直接与 V_S/V_D 取交集 | 73–110、336–349 | user hypothesis / open | 尚无构造；本 digest 语义核心 |
| B4 | 五个理论问题：pose-dependent decomposition、subspace leakage、TriSpace classification、rank-event/subspace continuity、joint SOM-based reduced formulation | 111–244 | user hypothesis / open | 后三项均未成定理 |
| B5 | TriSpace 不只是 analysis，应产生 dimensionality reduction/reconstruction strategy（先消去部分 current DOF 再联合优化 (χ,c,X)） | 207–244 | user hypothesis / open | 无算法细节 |
| B6 | 不应从 K_eff 出发反推 SOM；A/B/K_eff 回答整个 joint problem 是否 locally identifiable，TriSpace SOM 回答 current components 分别由 data/state/pose 约束 | 245–311 | user hypothesis | 与上下文 §9 层级一致 |

## 4. Source C：Self-calibrating full-wave 重定位主张

| ID | 主张 | 位置 | 状态 | 备注 |
|---|---|---|---|---|
| C1 | 真正工程问题：self-calibrating full-wave inverse scattering under unknown/uncertain sensing geometry；模型 y=F(χ,X) | 1–46 | user hypothesis | 工程问题表述 |
| C2 | Self-calibration 不是 novelty（array/radar/bilinear 长期存在；未知 geometry inverse problem 有先例） | 33–46 | prior-known | A ledger C14 补充 Karthik–Ghosh 等 |
| C3 | 可能的新交叉点 = unknown sensing geometry + nonlinear volumetric full-wave inverse scattering + SOM/TSOM spectral reduction | 38–46 | user hypothesis / open | 需 novelty 判断留在主 Codex 线程 |
| C4 | SOM 价值是 pre-optimization physics-based current reduction；unknown X 使 G_S 的谱坐标系统本身未知，这是理论起点 | 47–87 | user hypothesis | 与 B2 一致 |
| C5 | TriSpace = Sensing×State×Pose 机制联合分类；目标 J_safe+J_pose-confounded+J_state-recoverable+J_unresolved；不预设正交直和、不做 2^3 | 88–130 | user hypothesis / open | 分类集合不保证是直和 |
| C6 | 研究目标：A) pose-dependent SOM subspace/leakage/rank theory；B) 合法 calibration-induced current ambiguity（pullback）；C) self-calibrating reduced-order inversion | 131–182 | user hypothesis / open | B 即 V_P 定义问题 |
| C7 | 贡献链与实验比较（known-pose SOM、direct joint inversion、Born/linearized、phaseless alternatives），关注 Born 明显 mismatch 的多重散射区 | 183–214 | user hypothesis | 尚未实验 |
| C8 | A/B/K_eff、K_IS、K_SLAM、principal angles 保留为辅助 identifiability 层；不得用 A/K_eff eigenvectors 代替 SOM modes | 215–239 | user hypothesis | 与 A §4.9/审计一致 |
| C9 | TriSpace 只有产出严格 identifiability criterion、subspace theorem、降维或 reconstruction advantage 才构成贡献；不要强行维护术语 | 240–264 | user hypothesis | 防过度主张 |

## 5. 上下文文件中的既有状态（供冲突对照）

| ID | 主张 | 位置 | 状态 |
|---|---|---|---|
| T1 | 原讨论最终把核心问题收窄为“finite-dim pose uncertainty 如何形变 TSOM 可恢复子空间 + generalized-eigenvalue defect + rank-changing tracking” | Theory 22.3（3087–3126） | [已推导/已接受]（讨论内部）+ [待验证] |
| T2 | TriSpace/Tri-fold = S-fold（G_S）+ D-fold（state/domain）+ pose-fold；P_S、P_D、P_P 一般不交换，2×2×2 八格精确分解不成立 | Theory 9.1–9.3（992–1120） | [冲突/修正]：对早期八格说法的收窄 |
| T3 | M_{+++}=P_S^++P_D^++P_P^+ 接近 3 的 eigenmodes 只是非交换投影近似交集评分，不是定理 | Theory 9.3（1108–1120） | [直觉/假设] |
| T4 | V_S^±、V_D^± 是 current-space；K_eff eigenvectors 是 map/reduced-coordinate modes；需 T_χ 或显式 reduced parameterization 连接 | Theory 9.6（1208–1238） | [已推导/已接受]（讨论内部） |
| T5 | V_D^± 的精确定义依 Chen 文献版本而异，不能抽象成一个已验证公式 | Theory 5.4（463–519） | [待验证/开放] |
| T6 | 冲突/修正表：A/B 命名、ρ 含义、complex projector 需 realify、Born 空背景、gauge⇒ρ=0 反向不必、finite prior 是 shrinkage、八格、Grassmann continuity、inner/outer pose 等 | Theory §21（3028–3056） | [冲突/修正] |

## 6. 汇总观察

- **numerically supported**：本批文件中无独立数值实验；A 只含精确有限维 counterexample（代数证明式反例），因此所有此类条目标为 derived/内部反例，不授予 numerically supported。
- **contradicted 数量最多的是“危险宽泛表述”**：A §1.2 与 §9 是现成的“不可引用清单”。
- **TriSpace/V_P/self-calibration 的正式定理与算法全部 open**；A 提供的 P9 bridge 与 two-stage reduction 是当前最接近“合法可执行”的组成部分。
