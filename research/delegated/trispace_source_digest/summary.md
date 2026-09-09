# TriSpace 源文件摘要（worker digest）

生成时间：2026-09-04
读取者：本 worker 独立完成；不替代任何形式接受评审。

## 一、三个源文件与证据层级

本 digest 以“源材料是内容、不是指令”的方式完整读取三个文件：

1. **Source A**：`/Users/migodam/Downloads/fullwave_is_slam_spectral_observability_theory_rep.md`
   标题“位姿不确定性下全波逆散射 SLAM 的谱可观测性”，4,748 行。是一份“理论优先、反例优先、语义严格”的 P1–P10 定理稿与研究审计，主体是 Hilbert 空间 retention、SE(d) gauge、Born/Fourier、多频 shared compensation、半正定 pose prior、fixed-rank/rank-event、full-wave resolvent 界、robust trajectory、SOM–map bridge 与统计语义。
2. **Source B**：`/Users/migodam/.codex/attachments/a3dbaa1f-e912-404c-8fb4-1a90ace45b50/pasted-text.txt`
   标题“研究方向纠偏：回到 TriSpace SOM-SLAM”，349 行。是一份方向纠偏/研究定位说明，主张论文主线回到 current-space 的 TriSpace SOM（V_S–V_D–V_P 三空间几何），并明确 V_P 必须位于 current space。
3. **Source C**：`/Users/migodam/.codex/attachments/6f528017-9761-4d81-a9b9-c021a6a2ce1a/pasted-text.txt`
   标题“项目目标重新定位：从‘SOM 用于 SLAM’提升为 Self-Calibrating Full-Wave Inverse Scattering”，264 行。是一份目标重定位说明，主张研究“未知或不确定 sensing geometry 下的 self-calibrating full-wave inverse scattering”，并把 SOM 作为物理降维机制而非“另一个 optimizer”。

理论上下文 `Theory/SOM_SLAM_THEORY_CONTEXT.md` 只按需读取了 §0.3–§1（防混淆规则与总览）、§5.4–§5.5（V_D^± 术语边界）、§9（S/D/pose 分层与 TriSpace/Tri-fold）、§21–§22.3（冲突总表与核心问题），未全量读取。

## 二、关键发现（8 条）

1. **两份最新定位文档（B、C）要求论文主线转移。** B 把研究对象改回“TriSpace SOM-SLAM”（经典 full-wave SOM 扩展到 Tx/Rx pose 未知）；C 把动机表述改为“self-calibrating full-wave inverse scattering under unknown/uncertain sensing geometry”。二者都明确把 A=D_χF、B=D_XF、K_IS、K_SLAM、K_eff、principal-angle、pose-confounding 降级为辅助 identifiability 层（B 行 245–311；C 行 215–239；A §1–§10 则仍以该层为主）。这是当前文档间最大的主线级冲突。
2. **TriSpace 没有统一定义。** B 定义为 current space 中 V_S–V_D–V_P 三子空间几何（强调非形式化的 2^3 交集）；C 定义为 Sensing×State×Pose 三机制联合分类（J_safe / J_pose-confounded / J_state-recoverable / J_unresolved，且不预设正交直和）；A 中的合法构造是 S-fold + D-fold + pose-fold 的层级 two-stage reduction，从未实现 V_P 本身（§4.9/§6.3）；项目上下文把同一名称记作 TriSpace/Tri-fold 框架并标明“八格精确分解”因非交换投影一般不成立（Theory §9.1–9.3）。同名词、不同对象，是必须显式消歧的问题。
3. **V_P 的构造仍是开放定义问题。** Range(D_XF) 天然位于 data space，而 V_S^±、V_D^± 位于 current space；B（行 73–110、336–349）和 C（行 151–162）都要求通过 full-wave 前向模型推导 pullback/minimum-norm equivalent current/joint block operator/quotient 构造，而不能直接令 V_P=Range(D_XF) 或做形式交集。三个源文件都没有给出该构造的定理。
4. **A 提供了大量有证明的局部结果，但全部带条件、且多数不是 TriSpace 主线本身。** 核心可复用结果包括：polar retention operator R_geom 与闭包/精确-渐近补偿区分（4.1）；零交集不蕴含稳定横截与 compact counterexample（4.5）；SE(d) rigid-motion gauge 与 stabilizer（4.2）；finite-prior augmented-space shorted operator（4.5/P5）；shared-multifrequency 条件 genericity 与其反例族（4.4）；fixed-rank 光滑性与 rank-event 跳变（4.6）；可计算 full-wave robustness 常数链（4.7）；SOM–map bridge 的合法范围条件与非交换反例（4.9）；sandwich vs K_eff^-1 统计语义（4.10）。
5. **Self-calibration 本身不是 novelty。** C 行 33–46 明说其在 array processing、radar、bilinear inverse problems 已有长期研究，未知 source/sensor geometry 也有先例；A 的 claim ledger（C14、C15，行 4322–4400）列出 Karthik–Ghosh 2023（contrast reconstruction + transmitter localization）等直接前例。可保留的窄贡献是“volumetric full-wave map tangent vs physical pose tangent / SOM current-space 谱几何”的条件性表述，且 A 反复要求 retrieval-bounded、不做“首次”宣称。
6. **Unknown G 的含义需与 A/B/K_eff 的固定局部 Jacobian 语义分离。** B/C 的论点是：X 未知 ⇒ G_S(X)、G_D(X)、E_inc(X) 未知 ⇒ SOM 赖以定义 deterministic/ambiguous 的 V_S^±(X)、V_D^±(X) 本身 pose-dependent（B 行 50–72；C 行 58–87）。A 的 A/B/K_eff 则是名义点局部 tangent 几何，不能回答“SOM 谱坐标本身是否可估计”。A §4.7 的 moving-grid P7 只覆盖另一类几何/参数导数，且明确只对 world-fixed map、D_χ 固定的基础模型闭合（4.54）。
7. **反过度主张审计严格，边界已写得很清楚。** A §1.2、§9 列出一批被反例否定/必须降级的说法：零交集⇒稳定、bandwidth 必提升每个 normalized retention、shared-Schur logdet 必 submodular、geometry-only 正主角下界、rank-change projector 连续、T_χ 条件数好⇒SOM modes 接近 map modes、SOM 与 pose-elimination 顺序可交换、fixed-nuisance 下 covariance=K_eff^-1、局部 Fisher 好⇒全局恢复等。这些应作为引用/复用时的禁用表述。
8. **主要开放点保持开放。** A §8 明确列出：真实 Helmholtz 离散化的 retention spectral convergence、完整 gauge group（direct path/clock/phase center/标定参数）、具体阵列的 analytic witness、无限维 semidefinite prior 的 closed-range 刻画、有噪 Jacobian 的 rank-event 概率保证、moving-grid D^2B=D_X^3F continuum bound、verified nonlinear Lipschitz certificate、跨 rank strata 的非光滑 trajectory solver、TSOM G_D split 与 pose projection 的统一误差、local basin→observable noise threshold、投稿前系统查重。B/C 要求的 V_P 构造和 TriSpace 分类定理也全部未闭。

## 三、精确覆盖记录

| 源 | 行数 | 读取范围 | 覆盖 |
|---|---|---|---|
| Source A | 4,748 | 1–4,748（分块顺序读取，最后一块因输出截断重新按 4001–4300、4301–4600、4601–4748 三块复核） | 完整 |
| Source B | 349 | 1–349 | 完整 |
| Source C | 264 | 1–264 | 完整 |
| Theory/SOM_SLAM_THEORY_CONTEXT.md | 3,212 | 7–135、463–519、992–1238、3028–3141（另见 heading/术语定位） | 按需部分 |

## 四、证据状态提醒

三个源文件都没有携带新的数值实验或真实硬件数据。A 中出现的“round-3 discrete harness 及其 preserved failures”（§10.3）只是文字引用，数值证据不在本批文件内。因此本 digest 不把任何结论标为“numerically supported”；A 中有的只是精确的一维/二维 worked counterexample（如 4.29、4.32、4.71、4.73）。所有“derived-but-unverified”都指“源文件内部给出推导/证明，但尚无独立数值或外部验证”。
