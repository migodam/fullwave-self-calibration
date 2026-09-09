# [归档：方向修正前] 给网页版 GPT Pro 的高难理论、构造与数学问题总 Prompt

用途：把本轮自动科研后仍未解决、但可能决定论文上限的问题集中交给 GPT Pro。
建议把本文件全文连同以下本地文件上传：

- `THEORY_SEED.md`
- `FINAL_SCIENTIFIC_REVIEW.md`
- `PAPER_DRAFT.md`
- `validation/rank_saturation_results.json`
- `validation/validate_rank_saturation.py`
- 自动实验目录中的 `consolidated_round3_summary.json`
- 自动实验目录中的 `results_e11_vp_state_null.json`
- `Theory/SOM_SLAM_THEORY_CONTEXT.md`

下面横线内是一条完整 prompt，可直接复制。

---

你是一位专攻电磁逆散射、算子理论、可辨识性、数值优化和实验设计的资深理论
研究者。请对下面的研究做 **证明优先、反例优先、语义类型优先** 的深度攻坚。
你的任务不是夸奖想法，也不是把已有概念换名重述，而是尽可能完成严格推导；
如果命题不成立，必须给出最小反例并提出可存活的修正版。

## 一、研究对象与禁止偏移

研究对象是 **未知/漂移天线阵列几何下的 TriSpace SOM 与 full-wave inverse
scattering self-calibration**。核心直觉是：receiver-side Green operator 的变化
可否在一个明确声明的 SOM retained-current space 内，被规范地等效成电流变化；
同时必须保留无法被该 retained space 吸收的数据分量，并区分 transmitter motion
真正引起的 induced-current change。

“TriSpace”这里只允许表示同一 pose tangent 的 typed triple

\[
h\mapsto (Q_Uh,R_Uh,D_Uh),
\]

三项分别属于 current、data、state ambient spaces。禁止把它描述成三个平级、
可交换投影形成的 `2^3` 个交集；`Range(Q_U)` 只是 retained current space 内的
scene-dependent pose-equivalent tangent，不是第三个独立谱空间。

同时必须区分：

1. Chen SOM/Twofold SOM 的 `G_S,G_D` current-space objects；
2. map--pose local data tangents `A=D_chi F`、`B=D_xF`；
3. Fisher/Schur/principal-angle data-space objects；
4. independent-current variable `j` 与 map variable `chi`；
5. operator derivative `D_xG[h]` 与 finite perturbation `Delta G`。

不能把这些空间直接等同。若要比较，必须给出类型正确的 push-forward、pull-back、
quotient 或共同 data-space embedding。

## 二、已固定的模型与已证明有限维事实

在 world-fixed domain grid、固定均匀背景下：

\[
y=S(x)j,
\]

\[
j-D_\chi\{e(x)+Dj\}=0,
\qquad M=I-D_\chi D.
\]

局部线性化为

\[
\delta y=S\delta j+H_Sh,
\qquad H_Sh=(D_xS[h])j,
\]

\[
M\delta j-D_{E_{\rm tot}}\delta\chi-H_Dh=0,
\]

\[
H_Dh=D_\chi D_xe[h]
\]

（本模型 `D_xD=0`）。固定 contrast 且 `M` 可逆时，总 pose Jacobian 是

\[
Bh=H_Sh+SM^{-1}H_Dh.
\]

语义约束：`H_Sh` 是 receiver re-sampling，可被提升为 pseudo-current；
`M^{-1}H_Dh` 是 transmitter/illumination 引起的 physical induced-current change。
禁止把总 `Bh` 整体再提升并同时保留 `H_D`，否则 double count。

给定正交 retained current basis `U`：

\[
K_U=SU,
\qquad C_U=K_U^\dagger H_S,
\qquad Q_U=UC_U,
\qquad R_U=(I-K_UK_U^\dagger)H_S.
\]

有限维中已经证明：

\[
H_S=K_UC_U+R_U,
\qquad K_U^*R_U=0,
\qquad \operatorname{rank}C_U\le p,
\]

且 `Q_Uh` 是声明度量下 retained space 内的最小范数 lift。若 unrestricted `S`
满行秩，则 `R=0` 对所有 `H_S` 自动成立；所以“geometry error 可等效为某个
current”是空洞陈述。

令 `z=delta c+C_Uh`、`D_U=H_D+MUC_U`，有

\[
\delta y=K_Uz+R_Uh,
\]

\[
MUz-D_{E_{\rm tot}}\delta\chi-D_Uh=0.
\]

再令

\[
\mathcal N_D=MU\ker K_U+\operatorname{Range}(D_{E_{\rm tot}}),
\qquad T_U=(I-P_{\mathcal N_D})D_U.
\]

已经证明有限维 iff：一个 pose direction 能被 retained current 与 real contrast
同时隐藏且满足 data/state linearization，当且仅当

\[
R_Uh=0,\qquad T_Uh=0.
\]

注意：这个 iff 只是 solvability statement，实验没有证明 `T_U` 是有效恢复机制。

在白化实化数据空间 `Y_R`，令

\[
\mathcal N_U=
\operatorname{Range}_{\mathbb R}(\mathcal R(K_U))+
\operatorname{Range}(A_\chi),
\qquad
B_{\rm vis}=(I-P_{\mathcal N_U})B.
\]

已经证明

\[
\dim\ker B_{\rm vis}
\ge \max\{0,p+\dim\mathcal N_U-m\}.
\]

在当前实验 `m=2M,p=3,K=3` 且 nuisance columns generic independent 时，
`dim N_U=min(2M,2r+3)`，所以 `r=M-2` 至少隐藏两维，`r>=M-1` 隐藏三维。

## 三、必须接受的负证据，不得重新包装成正结果

1. unrestricted lift residual 约 `5.5e-16`，只证明 full-row-rank vacuity。
2. retained `r=4` 的代表性 leakage ratio 为 `0.784`，正交误差约 `8.9e-18`。
3. hard rank event 的 projector jump 为 `1.424`，soft filter 为 `0.127`；soft filter
   只改善坐标连续性，没有证明 self-calibration。
4. E4 的 data-hidden direction 对 data nuisance residual 为 `5.3e-16`，但 state
   residual 为 `0.947`、`T_U=0.733`。data-visible direction 的 `T_U≈0.674`；好坏
   pose estimator 的 `T_U` 大致都在 `0.81--0.89`。因此现有 `T_U` 不判别。
5. hard state equality 与 exact variable projection 均没有恢复 hidden directions；
   soft VP 的 map error 为 `0.838`、discarded-current fraction 约 `0.87`、完整物理
   residual 约 `0.83`。这是 biased reduced-model fit，不是 recovery。
6. 当 visible rank=3 时，所谓 reduced solver 只是 direct joint objective 的可逆
   正交重参数化；相等不是算法收益。当 visible rank<3 时，它按构造冻结 hidden
   coordinates，且比 direct joint 更差。
7. 旧实验在 `r>=M-1` 仅用 relative singular-value threshold，把 `1e-18` roundoff
   误判为 rank 3。正确结果是 nuisance saturation、visible rank 0、hidden rank 3。
8. half-circle 的 known-map pose-only oracle 三个 seeds 都“收敛”到 pose error 约
   `0.86` 的伪极小；所以 optimizer success flag 不等于 scientific success。
9. multi-transmitter/directional-source 只消除了 isotropic point-source orientation
   rank deficiency；它没有消除 global SE(2) map--pose gauge。
10. 当前代码中的 `chi` 是吸收了频率尺度的 scattering potential；跨 `k` 固定
    coefficient 不是 fixed-permittivity multifrequency experiment。

## 四、请完成的高难任务

### 任务 A：带权实化与 gauge quotient 下的统一局部可辨识性定理

请给出一个完整 theorem，统一包含：

- complex data、complex retained current coefficients、real map/pose variables；
- noise whitening 与任意正定 current metric；
- declared gauge group 的 infinitesimal tangent `G`；
- receiver pseudo-current 与 transmitter physical-current 的类型分裂；
- nuisance elimination 后的 exact local identifiability iff；
- dimension lower bound、equality conditions、rank-deficient counterexamples；
- stacked multi-source/multi-frequency/multi-frame measurements。

不要只写“Jacobian 满秩”。请给出所有 operator domains/codomains、realification
约定、projector metric、closed-range/rank 假设、quotient 结论，并说明它与 generic
bilinear identifiability up to transformation groups 的关系：哪些只是 Li--Lee--
Bresler 已有代数，哪些是 Helmholtz/current-specific physical content。

### 任务 B：无限维非闭值域下的规范 lift

有限维 `K_U^dagger` 的直接使用在 continuum compact sensing operator 下可能无界。
请构造一个可发表的 infinite-dimensional/regularized 版本：

1. 给出 Hilbert spaces、closed-range 或 source condition；
2. 比较 truncated SVD、Tikhonov、spectral filter 的 lift `Q_alpha` 与 leakage；
3. 给出 norm budget 下的最优表示问题；
4. 证明 existence/uniqueness/stability 或给出最弱足够条件；
5. 给出 `alpha -> 0` 时的 graph/closure convergence；
6. 说明什么意义下它仍可称“canonical”，什么情况下只能是 metric-dependent；
7. 提供一个 nonclosed-range counterexample，展示 exact equivalence 与 stable
   equivalence 的差别。

### 任务 C：stacked acquisition 的充分条件与最优设计

对多发射、多频、多帧索引 `ell=1,...,L`，定义各自的 nuisance ranges 与 pose
Jacobians。请推导：

- stacking 后 hidden dimension 的 sharp lower/upper bounds；
- shared pose、shared map、frame-specific current 的 block structure；
- 什么时候每一帧单独 hidden，但 stacking 后 quotient-identifiable；
- 必要条件与可检查的充分条件；
- 最小需要多少独立 source/frequency/frame；
- global gauge 与 source-symmetry stabilizer 如何分别处理；
- 设计目标 `maximize sigma_min,+`、D/E-optimal Fisher criteria 与 worst-case
  principal angle 之间的关系；
- 一个 tractable greedy/gradient algorithm 及其近似或单调性保证。

请提供至少两个反例：增加 measurements 但因为 nuisance 同步扩大而 rank 不增加；
增加多个同构 transmitter 但仍保留 orientation null。

### 任务 D：pose-dependent SOM/TSOM projector 的微分与 rank-event 理论

令 `U(x)` 来自 `S(x)` 的 retained singular subspace，或来自经过原文核验的 Twofold
SOM domain fold。请推导：

- fixed-rank stratum 上 spectral projector 的 Fréchet derivative；
- singular gap 控制的 perturbation bound；
- 与基向量符号/内部 rotation 无关的 Procrustes/parallel-transport update；
- threshold crossing 时 hard projector 不连续的精确条件；
- soft filter 的 bias/stability tradeoff；
- `Q_U(x),R_U(x),B_vis(x)` 对 pose 的完整导数，包括 `D_xU[h]`，而不是把 basis
  冻结；
- rank event 附近可用的 semismooth、manifold-with-strata 或 trust-region 方法。

请判断：是否能给出一个不依赖 arbitrary SVD basis 的局部 convergence theorem？
不能就给出最小反例。

### 任务 E：寻找真正能约束 hidden pose 的 state/physics 机制

当前 `T_U`、hard equality、soft penalty、exact VP 均失败。请不要微调同一个 scalar
penalty。请先解释失败的几何机制，然后构造至少三类本质不同的可检验替代：

1. cross-illumination/shared-contrast state consistency；
2. discarded-current energy 必须低于 noise-consistent bound 的约束；
3. passivity/causality/reciprocity/material-dispersion 等物理约束；
4. near-resonance 或 multiple-scattering sensitivity（如你认为合理）；
5. external calibration target 或 sparse fiducial current。

对每个机制给出：数学 formulation、为什么不是现有 nuisance tangent 的重复、
identifiability 条件、可能失败的反例、最小验证实验。若理论上 state equation 在
当前参数化中必然不能增加独立信息，也请证明这个 negative theorem。

### 任务 F：核验并类型化 Twofold SOM 的 D-fold

请回到 Zhong--Chen Twofold SOM 原始论文及可核验版本，给出准确式号和定义：

- `V_D^+/-` 到底来自 `G_D`、normal operator、材料加权算子还是其他对象；
- 第二折如何作用在 `V_S^-` 或其投影上；
- 变量、投影顺序、阈值、优化目标；
- 它与本文 `D_U` state defect 是否同型；
- 怎样在不混淆 current-space 与 data-space 的前提下加入 P-fold map--pose overlap；
- 最终能否形成一个严格的 S-fold → D-fold → P-fold 层级算法。

如果拿不到全文或式号，请明确停止，不要凭摘要或二手记忆编造公式。

### 任务 G：非线性 basin、有限孔径与 continuation

half-circle 中 known-map pose-only oracle 仍稳定落入伪极小。请：

- 构造该目标在 2D pose chart 中的局部/全局 basin 分析；
- 区分 aperture ambiguity、periodicity、source symmetry、map--pose gauge 与普通
  nonconvex local minima；
- 给出可视化但同时给出数学诊断；
- 研究 multi-start、frequency continuation、homotopy、trust-region radius、
  asymmetry injection 或 anchor 的最小修复；
- 给出一个“优化器收敛但科学失败”的正式 acceptance gate；
- 若可能，证明某种 small-error/local convexity 半径；不能就给反例与经验估计法。

### 任务 H：物理一致的 multifrequency/Maxwell 扩展

请把 2D scalar scattering-potential 公式改写为：

1. fixed material permittivity/conductivity 随频率进入的正确模型；
2. 3D dyadic Maxwell Green tensor、polarization、多分量 antenna response；
3. phase center、mutual coupling、gain/phase/clock nuisance；
4. 哪些 nuisance 可作为低维 pose-like tangent，哪些会使等效-current说法重新
   变得 vacuous；
5. dimension bound 中 real data dimension 与 nuisance count 如何改变；
6. 最小 high-fidelity simulation 和 measured-data calibration protocol。

### 任务 I：原创性与最窄贡献的全文级检验

逐项对照以下族群：Chen SOM、Zhong--Chen Twofold SOM、CSI/CC-CSI、virtual
experiments/virtual antennas、Bellomo 2014 antenna/phase-center calibration、
source/receiver-extension FWI、joint inverse scattering + transmitter
localization、blind gain/phase calibration、radar/SAR autofocus。

重点搜索是否已有以下完全或部分等价对象：

- restricted equivalent-current lift；
- lift plus range-complement leakage；
- receiver pseudo-current / transmitter physical-current split；
- retained-current/map nuisance saturation causing pose hiding；
- state-equation test of geometry-equivalent currents。

每条结论必须给可点击 primary-source citation、对应页码/式号/段落、相同点、不同点
和 novelty impact。不要说“未见报道”而不给检索边界。最后把贡献分成：

- 明确前人工作；
- textbook derivation；
- supported recombination；
- retrieval-bounded novel candidate；
- falsified；
- open。

## 五、必须输出的格式

请输出一份可直接交给论文作者的研究报告，依次包含：

1. **Executive verdict**：哪些问题被解决、哪些被推翻、论文应定位为 method、
   theory/negative-result paper 还是 cautionary letter。
2. **Notation and type table**：每个变量、domain、codomain、real/complex、metric。
3. **Theorem package**：逐条 theorem/lemma/corollary、完整证明、假设必要性、最小
   反例。不能只给 proof sketch。
4. **Relation map**：SOM S-fold、Twofold D-fold、map--pose P-fold、joint Jacobian、
   gauge quotient 的类型正确关系。
5. **Algorithm**：不依赖隐藏真值的 pseudocode、复杂度、rank threshold、失败
   触发器和 acceptance gate。
6. **Experiment package**：每个新实验的 hypothesis、null、设置、metrics、seeds、
   图表、成功/失败判据；优先做能推翻理论的实验。
7. **Claim ledger**：每条 claim 标成 `[proved]`、`[numerically supported]`、
   `[falsified]`、`[open]`、`[prior art]`。
8. **Manuscript-ready replacements**：给出可直接替换到 abstract、contributions、
   limitations 和 conclusion 的英文段落。
9. **Citation audit**：只使用可核验 primary sources；逐条注明证据层级。
10. **Unresolved list**：只保留真正需要后续计算、全文访问或新数据的问题。

## 六、质量与诚信约束

- 不得把 pseudoinverse、projector、rank-nullity、Schur complement 等 textbook
  数学包装成新定理；可以把它们作为物理特化的证明工具。
- 不得把 “full row rank 下可以等效” 写成正贡献。
- 不得把 `T_U` 现有负结果改写成“初步成功”。
- 不得把 reduced=direct 写成加速或精度收益，除非补充并通过 runtime/DOF/basin
  实验。
- 不得把 source-symmetry rank deficiency 称为 global gauge。
- 不得把当前 sensing-SVD 截断代码称为已实现的 Twofold SOM。
- 不得忽略真实参数是 real、current coefficients 是 complex、投影必须 noise-
  whitened 这一事实。
- 不得用优化器 `success/status` 代替 residual、parameter error、basin 与 oracle
  comparison。
- 不得宣称全球新颖性；所有 novelty conclusion 必须附检索范围。
- 如果一个核心命题无法证明，请优先给反例和最窄可存活命题，不要输出车轱辘话。

请现在直接完成尽可能多的证明与设计，不要只返回研究计划或建议我“进一步研究”。

---
