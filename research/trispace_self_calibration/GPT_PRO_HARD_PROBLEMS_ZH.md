# 给网页版 GPT Pro 的核心理论与方法设计总 Prompt

用途：把 TriSpace SOM-SLAM 当前真正决定论文上限的数学、算法与原创性问题，
集中交给 GPT Pro 做一次证明优先的完整攻坚。

本轮不预设这项工作已经足以支撑 TGRS/TAP。请形成当前最强可存活理论与设计，
并对 SOM-specific method claim 给出明确 pass/fail verdict；如果失败，必须给
降级后的最窄论文定位。本轮不要求实际运行新实验。实验应被设计成可执行、可证伪
的方案，但不得把计划写成结果。

建议与本文件一起上传：

- PAPER_DRAFT.md（方向修正后的主稿）
- PAPER_DRAFT_GEOMETRY_LIFTED_2026-09-04.md（旧方向与负结果）
- FINAL_SCIENTIFIC_REVIEW.md
- THEORY_SEED.md
- ORIGINALITY_AUDIT_ZH.md
- validation/rank_saturation_results.json
- experiments/idea_loops/loop_2026-09-04_02-58-31/experiment_geometry_lifted_trispace_som/consolidated_round3_summary.json
- experiments/idea_loops/loop_2026-09-04_02-58-31/experiment_geometry_lifted_trispace_som/results_e11_vp_state_null.json
- Theory/SOM_SLAM_THEORY_CONTEXT.md
- Chen SOM、Zhong--Chen Twofold SOM、Pan 等 phaseless SOM 原文；若无法上传，
  要求 GPT Pro 明确标注无法核验的式号

以下横线内是一条完整 prompt，可以直接复制给 GPT Pro。

---

你是一位专攻电磁逆散射、SOM/Twofold SOM、统计信息几何、算子理论、
可辨识性、数值优化和阵列自校准的资深理论研究者。请给出下面研究当前最强的
核心理论与算法设计，并判断它是否真的足以支撑 IEEE TGRS 或 IEEE TAP 方法
论文；不得预设肯定结论。

工作方式必须是：

1. **证明优先**：先给假设、类型、定理和证明，再给解释。
2. **反例优先**：命题若不成立，给最小反例和最窄可存活修正版。
3. **物理语义优先**：不能把不同 ambient space 的对象靠符号相似强行等同。
4. **原创性克制**：joint pose + image、self-calibration、frequency continuation、
   Schur complement、principal angle 均不是新概念。
5. **方法可实现**：最终算法必须只使用观测、当前迭代量、噪声模型与可声明先验，
   不得使用真值、oracle rank 或事后挑阈值。
6. **不要只返回计划**：请直接完成尽可能多的推导、证明、反例、算法和可替换论文
   段落。

## 一、不可再偏移的论文主线

论文不再以“Green operator 的变化能否等效成 current 变化”为核心终点。
unrestricted current 在 full-row-rank sensing 下可以拟合任意单帧数据扰动，这只是
表达能力冗余，不能增加 non-scattering information，也不能自动改善
\(K_{\mathrm{eff}}\)。

论文的三句核心论点必须保持为：

1. **不要为了规避 geometry sensitivity 而删除 phase；如果 phase corruption
   主要来自低维 pose，就估计并校准这个 pose。**
2. **SOM retained rank 必须同时服务 reconstruction 与 calibration，因为新增
   current freedom 会吸收 pose 信息。**
3. **用低频、低秩获得 calibration basin，再在 pose uncertainty 降低后逐步解锁
   中高频和高秩 current/map 信息。**

候选方法名为：

> phase-preserving rank-adaptive self-calibrating SOM，简称 PRASC-SOM。

它必须保留 complex coherent data。所谓 multifrequency feature fusion 必须指：
在正确材料色散模型、共同 pose/trajectory、共同 map 与 acquisition-specific
current 的物理约束下，做 whitened coherent stacking；不能只是把不同频率特征
随意拼接。

论文价值不能建立在以下说法上：

- “jointly estimate pose and image”；
- “高频 phase 对位置敏感”；
- “用 frequency continuation”；
- “用 Schur complement/Fisher information/principal angle”；
- “phaseless SOM 没有 current-space constraint”。

最后一条是错误的。Pan--Zhong--Chen--Yeo 的 phaseless SOM 仍然使用 spectrum
analysis 以及 deterministic/ambiguous contrast-source partition；只是缺少 phase
后 deterministic recovery 被修改。请按这个事实比较，不得树稻草人。

真正候选的 SOM-specific 科研价值只能是：

\[
\boxed{
\text{retained-rank control}
+\text{dual }\rho/B_{\mathrm{vis}}\text{ observability}
+\text{phase-aware rank/frequency continuation}
}
\]

能否使 self-calibration 比一个充分调优的 direct joint full-wave inversion
具有更大 basin、更低成本、更高精度，或者至少给出可验证的 rank/acquisition
安全策略。如果不能，请明确判定 method claim 不成立。

## 二、类型固定的全波模型

对 acquisition/frequency index \(\ell\)，使用

\[
y_\ell=S_\ell(x)j_\ell+\varepsilon_\ell,
\]

\[
j_\ell-D_{\chi_\ell}\{e_\ell(x)+D_\ell j_\ell\}=0.
\]

变量类型：

- \(y_\ell\in\mathbb C^{M_\ell}\)：complex coherent data；
- \(j_\ell\in\mathbb C^n\)：induced/contrast current；
- \(\chi\in\mathbb R^q\)：real material/map parameter；
- \(x\in\mathbb R^p\)：real pose/trajectory parameter；
- \(\delta c_\ell\in\mathbb C^{r_\ell}\)：retained current coefficient；
- \(h\in\mathbb R^p\)：pose tangent。

在 world-fixed domain grid、固定均匀背景下：

\[
M_\ell=I-D_{\chi_\ell}D_\ell,
\]

\[
\delta y_\ell=S_\ell\delta j_\ell+H_{S,\ell}h,
\qquad
H_{S,\ell}h=(D_xS_\ell[h])j_\ell,
\]

\[
M_\ell\delta j_\ell
-D_{E_{\mathrm{tot},\ell}}\delta\chi_\ell
-H_{D,\ell}h=0,
\]

\[
H_{D,\ell}h=D_{\chi_\ell}D_xe_\ell[h].
\]

必须保持语义：

- \(H_S h\) 是 receiver re-sampling，可 lift 成 pseudo-current；
- \(M^{-1}H_Dh\) 是 transmitter/illumination 引起的 physical induced-current
  change；
- 禁止把 total \(B h\) 整体 lift，同时又保留 \(H_D\)，否则 double count。

局部共同 data-space tangent 写成

\[
\delta y_\ell
=K_{\ell,r}\delta c_\ell
+A_{\chi,\ell}\delta\chi
+B_\ell h,
\qquad
K_{\ell,r}=S_\ell U_{\ell,r}.
\]

这三项首先是一个 conservative augmented envelope，不得自动当成单一物理
parameterization：pure independent-current data block 中
\(A_{\chi,\ell}=0\)，map 通过 state equation 进入；pure reduced-state model
中 \(A_{\chi,\ell},B_\ell\) 是 total derivatives，而 free
\(K_{\ell,r}\delta c_\ell\) 应删掉，除非明确引入额外 current-correction variable
及其 state coupling。请分别推导，禁止 double count。

请明确区分：

1. independent-current SOM formulation；
2. 消去 state 后的 reduced full-wave map/pose formulation；
3. data-only nuisance analysis；
4. state constraint 是否真正减少 nuisance，而不是重复同一数据。

## 三、白化、实化与 gauge

令 \(W_\ell^*W_\ell=\Sigma_{\varepsilon,\ell}^{-1}\)。complex coefficient 的
realification 为

\[
\mathcal R(A)=
\begin{bmatrix}
\Re A&-\Im A\\
\Im A&\Re A
\end{bmatrix},
\]

real coefficient 的 embedding 为

\[
\mathcal E(A)=
\begin{bmatrix}
\Re A\\
\Im A
\end{bmatrix}.
\]

因此

\[
\bar K_{\ell,r}=\mathcal R(W_\ell K_{\ell,r}),\quad
\bar A_\ell=\mathcal E(W_\ell A_{\chi,\ell}),\quad
\bar B_\ell=\mathcal E(W_\ell B_\ell).
\]

一般的 rigid map--pose gauge 是
\(\mathcal G\subset\mathbb R^q\times\mathbb R^p\) 的 joint tangent，不一定能拆成
独立的 \(\mathcal G_x\oplus\mathcal G_\chi\)。请先给 joint quotient；只有在
external anchor 或明确 gauge-fixing constraint 使 map/pose charts 可分后，才可
使用独立的 \(Z_x,Z_\chi\)。所有 rank、projector、Fisher、principal-angle
结论必须在白化、实化、joint quotient 后陈述。

如果使用任意正定 current metric \(M_c\)，请给出 metric pseudoinverse 与
metric projector，不能偷偷退回 Euclidean norm。

## 四、当前已接受的有限维事实

### 4.1 Retained receiver lift

对 receiver tangent：

\[
C_r=(S U_r)^\dagger H_S,\qquad
Q_r=U_rC_r,\qquad
R_r=(I-P_{\operatorname{Range}(S U_r)})H_S.
\]

已接受：

\[
H_S=S U_rC_r+R_r,\qquad
(S U_r)^*R_r=0,\qquad
\operatorname{rank}C_r\le p.
\]

如果 unrestricted \(S\) 满行秩，则 \(R=0\) 对任意 \(H_S\) 自动成立；这不是
positive contribution。

### 4.2 Pose-side visibility

在 quotient real data space 中：

\[
\mathcal C_r=\operatorname{Range}(\bar K_r),
\qquad
\mathcal N_r=\mathcal C_r+\operatorname{Range}(\widehat A),
\]

\[
B_{\mathrm{vis}}(r)=P_{\mathcal N_r^\perp}\widehat B,
\qquad
J_x(r)=B_{\mathrm{vis}}(r)^TB_{\mathrm{vis}}(r).
\]

\(B_{\mathrm{vis}}\) 回答：

> current/map 都可以变化时，还有哪些 pose directions 具有独有 data signature？

### 4.3 Map-side survival

先消去 retained-current nuisance：

\[
\Pi_r=P_{\mathcal C_r^\perp},
\qquad
\widetilde A_r=\Pi_r\widehat A,
\qquad
\widetilde B_r=\Pi_r\widehat B.
\]

\[
K_0(r)=\widetilde A_r^T\widetilde A_r,
\]

\[
K_{\mathrm{eff}}(r)
=\widetilde A_r^T
(I-P_{\operatorname{Range}(\widetilde B_r)})
\widetilde A_r.
\]

定义

\[
K_{\mathrm{eff}}(r)v_i=\rho_i(r)K_0(r)v_i.
\]

\(\rho_i\) 回答：

> known-pose 时可恢复的 map mode，在 pose unknown 后还保留多少相对信息？

无 pose prior 且限制在 \(K_0\) observable support 时，\(\rho_i\) 与
\(\sin^2\theta_i\) 对应，\(\theta_i\) 是
\(\operatorname{Range}(\widetilde A_r)\) 与
\(\operatorname{Range}(\widetilde B_r)\) 的 principal angles。

\(\rho\) 是 generic nuisance-elimination relative spectrum，不是逆散射专属数学；
逆散射特异性来自 \(A,B,S\) 的 full-wave physics 和 SOM current projection。

### 4.4 Nested-rank pose information theorem

固定 linearization、frequency、metric、whitening 和 gauge，若

\[
\operatorname{Range}(U_r)\subseteq\operatorname{Range}(U_{r+1}),
\]

则

\[
\mathcal N_r\subseteq\mathcal N_{r+1},
\]

\[
P_{\mathcal N_{r+1}^\perp}
\preceq
P_{\mathcal N_r^\perp},
\]

\[
J_x(r+1)\preceq J_x(r).
\]

若

\[
\mathcal E_r=\mathcal N_{r+1}\cap\mathcal N_r^\perp,
\]

则

\[
J_x(r)-J_x(r+1)=\widehat B^TP_{\mathcal E_r}\widehat B.
\]

因此 equality 当且仅当
\(P_{\mathcal E_r}\widehat B=0\)。请保留并加强这个 theorem。

注意：\(\rho_i(r)\) 的 numerator 和 denominator 都随 \(r\) 变，一般不保证
monotone。请不要未经证明写成 monotonic。

## 五、必须保留的负证据

这些证据是设计约束，不得包装成新方法成功：

1. unrestricted lift residual 约 \(5.5\times10^{-16}\)，只证明 full-row-rank
   vacuity；
2. retained rank 的 leakage 非零，说明 declared current coordinates 不能吸收
   所有 receiver tangent，但单独不能证明 pose recoverable；
3. hard rank crossing 的 projector jump 很大；soft filter 只改善连续性；
4. \(T_U\)、hard state equality、soft penalty、exact variable projection 在当前
   实验均未恢复 hidden pose；
5. visible rank 完整时，旧 reduced solver 只是 direct objective 的可逆正交
   reparameterization；rank 不完整时则冻结 hidden coordinates；
6. nuisance saturation 时，relative-only threshold 把 roundoff 误判成非零 rank；
   正确结论是 \(r\ge M-1\) 时 tested model 的 visible pose rank 为零；
7. limited aperture 中 known-map pose-only oracle 也会收敛到错误局部极小；
8. multiple/directional transmitters 只消除 source symmetry stabilizer，不消除
   global map--pose gauge；
9. 当前频率 sweep 使用 frequency-scaled scattering potential，不能当作
   fixed-permittivity multifrequency material experiment。

state/passivity/resonance 路线降为次要任务。除非先证明它带来独立信息，否则不要
再用新的 scalar penalty 重复失败路线。

## 六、按优先级完成的核心任务

## G0：先过“为什么必须是 SOM”生死门

在展开大量 textbook 推导之前，先用最小模型回答：

1. PRASC-SOM 相对相同 forward physics、data、initialization、prior 与计算预算的
   direct joint full-wave inversion，究竟多了什么不可被可逆 reparameterization
   消掉的机制？
2. 给一个二者完全等价的模型、一个 adaptive rank 条件性有利的模型、一个
   truncation bias 使其更差的模型。
3. 判断候选增益能否被定理化为 basin、conditioning、complexity 或 acquisition
   guarantee；若都不能，立即把 method claim 判为 fail。
4. 即使 fail，也继续完成有独立价值的 rank-information diagnostic theory，但
   后续全文必须使用降级定位，不能再以“可发方法论文”为目标反推结论。

先输出一段 provisional verdict，再做 P0--P5；完成 P4.3 后必须更新为正式 verdict。
如果输出长度受限，优先顺序为 G0、P1.2、P1.3、P2、P3、P4.3、P5，然后才是
P1.1 的标准证明复核、P6、P7 和次优先级附录。

## P0：先做全局类型审计

请给出完整 notation/type/metric/gauge table，并逐式审计：

- domain/codomain；
- real/complex；
- independent-current 与 reduced-state 语义；
- receiver pseudo-current 与 transmitter physical current；
- shared map、shared pose、frame-specific current；
- whitening 与 metric projector；
- global gauge 与 source symmetry stabilizer。

发现任何类型错误时，先改模型，再继续证明。不要靠维度刚好相同掩盖错误。

## P1：coherent phase 相对 phaseless 的信息定理

其中 P1.2 与 P1.3 是最高优先级；P1.1 已有 draft proof，先做严谨性复核，不要
让标准 Fisher data-processing 重推导挤占核心问题。

### P1.1 Matched observation theorem

设 coherent observation \(Y\sim p_\theta\)，phaseless observation \(Z\) 由
parameter-independent channel \(q(z\mid y)\) 产生，包含
\(Z=|Y|\) 或 \(Z=|Y|^2\)。请：

1. 严格证明 score identity
   \(s_Z=\mathbb E[s_Y\mid Z]\)；
2. 严格证明 joint Fisher data processing
   \(J_Z\preceq J_Y\)；
3. 给 equality 的必要充分条件；
4. 处理 complex Gaussian coherent noise、real intensity noise 与 singular
   Fisher；
5. 给出一个 phase-only pose tangent 被 intensity Jacobian 消去的最小例子；
6. 给出一个 equality 或近 equality 的反例，避免把 inequality 夸大成 strict。

### P1.2 Nuisance-eliminated efficient information

目标参数是 pose 或 map，其他变量是 nuisance。请证明或否定：

\[
J_{\mathrm{eff},Z}^{\mathrm{pose}}
\preceq
J_{\mathrm{eff},Y}^{\mathrm{pose}},
\qquad
J_{\mathrm{eff},Z}^{\mathrm{map}}
\preceq
J_{\mathrm{eff},Y}^{\mathrm{map}},
\]

其中 effective information 用 quotient/Schur/efficient score 定义。

不能仅由 full Fisher 的 blockwise Loewner order 直接跳到 Schur complement。
请用 efficient-score projection、statistical experiment comparison 或明确附加条件
证明；若一般不成立，给最小反例与正确版本。

### P1.3 Phase sensitivity 与 basin

从 2D Hankel Green function 和 3D Helmholtz/dyadic asymptotic 出发，推导：

\[
\delta\phi\approx k\,\delta r,
\]

并明确：

- local pose Jacobian 何时为 \(O(k)\)；
- fixed coherent SNR 下 Fisher 何时为 \(O(k^2)\)；
- amplitude、spreading、noise scaling 如何改变结论；
- phase wrap/cycle skipping 的 ambiguity spacing 为什么是 \(O(k^{-1})\)；
- limited aperture、multi-path、multiple scattering、source symmetry 如何改变
  local basin；
- 一个可以用于 algorithm 的 sufficient local phase-error condition。

最终给出可执行的 frequency gate，例如

\[
k_{\mathrm{next}}
\sqrt{\lambda_{\max}
(L_{\mathrm{path}}\Sigma_xL_{\mathrm{path}}^T)}
\le\gamma_\phi,
\]

并说明它是 theorem、bound、heuristic 还是需要校准的 surrogate。

## P2：统一 \(\rho\)、\(K_{\mathrm{eff}}\)、\(B_{\mathrm{vis}}\) 与 SOM margin

请建立一套双侧 spectral theorem。

### P2.1 Map side

完整处理：

- singular \(K_0\) 的 support restriction；
- generalized eigenvalue multiplicity；
- \(\rho=0,1\) 的几何意义；
- pose prior
  \(\Lambda_x\succeq0\) 下的 regularized Schur form；
- colored noise、arbitrary metric、gauge quotient；
- \(\rho_i=\sin^2\theta_i\) 的精确维数和补一 eigenvalues；
- prior 存在时还能否解释成 modified principal angle。

### P2.2 Pose side

对

\[
B_{\mathrm{vis}}
=P_{(\operatorname{Range}\bar K_r+\operatorname{Range}\widehat A)^\perp}
\widehat B
\]

给出：

- exact local quotient-identifiability iff；
- smallest positive singular value 与 local stability bound；
- dimension lower bound 及 sharpness/equality；
- map prior、pose prior 和 shared-state constraint 如何改变 nuisance projector；
- rank deficiency、near-rank deficiency 与 numerical threshold 的区别。

### P2.3 Duality theorem

请准确说明 \(\rho\) 与 \(B_{\mathrm{vis}}\) 在何种 block Fisher matrix 上构成
map-side/pose-side Schur duality，哪些非零谱有关系，哪些没有。不要说它们
“相等”。如果存在 canonical correlations 或 generalized singular values 的共同
表示，请给定理、证明和维数条件。

最终形成三谱解释：

\[
\sigma_i(S):
\text{ current mode 是否 data-supported},
\]

\[
\rho_i:
\text{ recoverable map mode 是否能在 pose uncertainty 下存活},
\]

\[
\sigma_i(B_{\mathrm{vis}}):
\text{ pose direction 是否能在 current/map nuisance 下存活}.
\]

## P3：retained rank 对 calibration 的单调侵蚀与 Pareto 设计

请把 nested-rank theorem 提升为可投稿版本：

1. 给 finite-dimensional weighted/realified/gauge theorem 和完整证明；
2. 给 equality、strict decrease、rank drop 的必要充分条件；
3. 给 dimension obstruction 的 sharpness 条件；
4. 分析 arbitrary nested basis、SVD basis、soft spectral filter；
5. 区分 hard rank event 与 fixed-rank observability event；
6. 推导 projector derivative 与 singular-gap perturbation bound；
7. 给 rank crossing 处 hard projector 不连续的最小反例；
8. 给 soft filter 的 bias--stability--pose-information tradeoff；
9. 构造 \(\rho(r)\) 非单调的最小反例；
10. 建立 reconstruction utility 与 calibration information 的 Pareto frontier；
11. 把 free retained-current nuisance 与 state-feasible current tangent 分开：
    证明 state elimination 后 effective nuisance 何时仍 nested；若不 nested，
    给反例并给正确的条件性单调 theorem。

请判断以下 rank policy 是否充分：

\[
\mathcal F_k=
\{r:
\eta_S(r)\ge\gamma_S,\;
g_S(r)\ge\gamma_{\mathrm{gap}},\;
\eta_x(r)\ge\gamma_x,\;
\eta_\rho(r)\ge\gamma_\rho\}.
\]

其中

\[
\eta_x(r)
=\frac{\lambda_{\min}(J_x(r))}
{\lambda_{\max}(\widehat B^T\widehat B)+\epsilon},
\qquad
\eta_\rho(r)=\min_{i\in\mathcal T_k}\rho_i(r).
\]

请修正 normalization、singular case、target-mode selection、threshold statistical
calibration、rank hysteresis 与 soft-filter 对应版本。不能使用 ground truth。

## P4：把理论变成 PRASC-SOM，而不是静态 reparameterization

这是第二个最高优先级。

请给出一个完整、可实现的 phase-preserving rank-adaptive self-calibration
algorithm。至少包含：

1. low-frequency/low-rank 初始化；
2. full complex data 与 full-wave state；
3. 当前 pose 下重建 \(S,e,D,U_r\)；
4. \(\sigma(S)\)、gap、\(B_{\mathrm{vis}}\)、\(J_x\)、\(K_0\)、
   \(K_{\mathrm{eff}}\)、\(\rho\) 的计算；
5. calibration-aware rank/filter selection；
6. fixed-rank stratum 上的 joint GN/LM 或 variable projection；
7. rank event 后 projector transport、model rebuild、curvature reset；
8. fixed-rank \(\rho\) crossing 后 map-mode reclassification；
9. phase uncertainty 允许时提高 frequency；
10. pose visibility 不足时触发 acquisition adaptation；
11. full-physics residual、noise consistency、gauge、spectral stability 和 parameter
    change 的 acceptance gate；
12. 明确 failure output，而不是强行返回 calibrated pose。

### P4.1 Fixed-stratum convergence

给出局部 convergence theorem，至少明确：

- objective；
- smoothness；
- fixed-rank projector dependence \(U_r(x)\)；
- Jacobian 是否包含 \(D_xU_r[h]\)；
- gauge quotient；
- trust-region/LM 条件；
- basin 与 singular gap；
- inexact current/map solves 的容差。

如果只能证明 descent 或 subsequence stationarity，请诚实给最强可证版本。

### P4.2 Event-safe dynamics

严格区分：

- \(r\to r\pm1\)：rank event，跨
  \(\mathrm{Gr}(r,n)\) strata；
- rank 不变而 \(\rho_i\) 穿过阈值：observability event；
- singular gap 闭合但 rank 尚未改变：coordinate-instability event。

请设计 semismooth、manifold-with-strata、continuation 或 trust-region event
handling，并证明有限事件、descent 或给出 cycling counterexample 与 hysteresis
修正。

### P4.3 SOM-specific advantage or failure theorem

把 PRASC-SOM 与使用相同 forward model、初值、prior、frequency data、optimizer
预算的 direct joint full-wave inversion 比较。

请给出至少一种可证明的非平凡优势：

- early-stage dimension/conditioning advantage；
- larger local basin sufficient condition；
- reduced forward/adjoint complexity；
- rank policy 避免 pose-hiding failure；
- acquisition decision guarantee。

如果一般不可能证明，请给出：

1. 一个二者完全等价的最小模型；
2. 一个 adaptive rank 严格有利的条件性模型；
3. 一个 adaptive rank 因 truncation bias 更差的模型；
4. 最窄、可实验验证的 method claim。

这是判断论文是否真的有 SOM-specific 科研价值的关键，不能回避。

## P5：stacked acquisition 的充分条件与设计

对 frame-specific current、shared map、shared pose：

\[
\bar K_r=\operatorname{blkdiag}
(\bar K_{1,r_1},\ldots,\bar K_{L,r_L}),
\]

\[
\bar A=
\begin{bmatrix}\widehat A_1\\\vdots\\\widehat A_L\end{bmatrix},
\qquad
\bar B=
\begin{bmatrix}\widehat B_1\\\vdots\\\widehat B_L\end{bmatrix}.
\]

令

\[
\Pi_C=P_{\operatorname{Range}(\bar K_r)^\perp},
\quad
\widetilde A=\Pi_C\bar A,
\quad
\widetilde B=\Pi_C\bar B,
\]

\[
J_{x,\mathrm{stack}}
=\widetilde B^T
(I-P_{\operatorname{Range}(\widetilde A)})
\widetilde B.
\]

请完成：

1. stacked quotient-identifiability 的必要充分条件；
2. hidden dimension 的 sharp bounds；
3. 每帧单独 hidden、stack 后 identifiable 的严格条件；
4. shared-map compensation 不一致为何能产生新信息；
5. 最小 source/frequency/frame 数量的条件性 bound；
6. global gauge 与 source symmetry stabilizer 的分别处理；
7. duplicated/symmetric acquisition 不增加 rank 的反例；
8. nuisance 随 measurement 同步扩张而无收益的反例；
9. E/D/A-optimal、worst principal angle、\(\min\rho\) 的多目标关系；
10. 可实现 greedy/gradient design，并给 monotonicity、submodularity
    approximation 或明确无保证；
11. acquisition 触发条件如何嵌入 PRASC-SOM。

不能把每帧分别消去 shared map 后的 Gram 直接相加，除非证明该操作与共同 map
Schur complement 等价。

## P6：精确接入 SOM/Twofold SOM

### P6.1 SOM-specific current model

回到 Chen SOM 原文，给出：

- deterministic current 与 ambiguous current 的精确定义；
- rank/threshold 的原始选择；
- 当前 \(U_r\)、\(K_r\)、current nuisance 与原 SOM 的对应；
- pose-dependent \(S(x)\) 下 retained projector 的完整 derivative；
- 哪些新增项会被 frozen-basis implementation 漏掉。

### P6.2 Twofold SOM

回到 Zhong--Chen 原文，逐式核验：

- 第二 fold 的 operator、domain/codomain；
- \(V_D^+,V_D^-\) 的来源；
- 第二 fold 如何作用于第一 fold 的 ambiguous current；
- projector 顺序、阈值和优化 objective；
- 它与本项目 state defect \(D_U\) 是否同型。

然后构造类型正确的层级：

\[
\text{S-fold current support}
\longrightarrow
\text{D-fold domain/state refinement}
\longrightarrow
\text{P-side map--pose observability}.
\]

P-side 不是第三个 current-space SVD，不能和 S/D spaces 直接求交。必须通过
共同 data tangent、push-forward 或 Schur complement 连接。

如果无法访问全文或式号，明确写“未核验”，禁止凭摘要补公式。

## P7：原创性审计与 TGRS/TAP 最窄定位

只用可核验 primary sources，逐项对照：

- Chen SOM；
- Zhong--Chen Twofold SOM；
- Pan 等 phaseless SOM；
- CSI/CC-CSI；
- microwave antenna/phase-center calibration；
- joint inverse scattering plus source/transmitter localization；
- source/receiver-extension FWI；
- blind gain/phase calibration；
- movable-array/array-position self-calibration；
- SAR/radar autofocus；
- multifrequency continuation。

重点检索是否已有以下组合：

1. SOM retained rank 作为 pose-calibration information control；
2. nested retained-current space 导致 pose Fisher/Gram Loewner monotonic loss；
3. \(\rho/K_{\mathrm{eff}}\) 与 \(B_{\mathrm{vis}}\) 的双侧 gating；
4. phase-preserving low-frequency/low-rank 到 high-frequency/high-rank dual
   continuation；
5. observability-driven acquisition adaptation。

每条文献必须给：

- DOI/正式链接；
- 可核验页码、式号或原文段落；
- 与本文完全相同、部分相同和不同的部分；
- 对 novelty 的影响；
- 检索范围与访问边界。

最后把所有贡献分成：

- [prior art]
- [textbook]
- [proved recombination]
- [retrieval-bounded novel candidate]
- [falsified]
- [open]

并分别给出 TGRS 和 TAP 两种最窄可投稿定位。不要承诺“能发”，要列出每个期刊
仍缺的理论与证据 gate。

## 七、次优先级附录任务

只有完成 P1--P7 后再做。

### S1：continuum/regularized lift

处理 compact operator、nonclosed range、unbounded pseudoinverse、TSVD/Tikhonov/
spectral filter、source condition、graph convergence。说明 stable equivalence 与
exact equivalence 的区别。它是数学附录，不应抢占主论文 self-calibration 叙事。

### S2：3D Maxwell 与真实阵列 nuisance

扩展到 dyadic Green tensor、polarization、multi-component antennas、phase center、
mutual coupling、gain/phase/clock drift。区分哪些是低维 pose-like nuisance，哪些
会让 current equivalence 再次 vacuous。给 TAP 路线的最小 electromagnetic
validation protocol。

### S3：state/physics constraints

只有在证明 state equation、cross-illumination consistency、passivity、causality、
reciprocity 或 dispersion 增加独立信息后，才把它们加入主算法。先解释旧
\(T_U\) 失败机制，再给 negative theorem 或本质不同的约束。禁止只换一个 penalty
weight。

## 八、必须交付的最终格式

请按下列顺序输出一个完整研究包。

### 1. Executive scientific verdict

用一页以内回答：

- 论文最强可存活 claim；
- 哪些 claim 已证明；
- 哪些只是设计；
- 哪些命题被推翻；
- 当前更像 TGRS method paper、TAP calibration paper，还是仍不足以投稿；
- “为什么必须是 SOM”是否得到非平凡回答。

### 2. Type and assumption ledger

表格列出每个变量、domain/codomain、real/complex、metric、gauge、shared/
frame-specific 属性。每个 theorem 单列最小假设。

### 3. Complete theorem package

逐条给：

- theorem/lemma/corollary；
- 完整证明，不只 proof sketch；
- 假设必要性；
- equality/strictness；
- 最小反例；
- 对算法的具体作用；
- epistemic label。

至少覆盖 P1、P2、P3、P5。

### 4. PRASC-SOM algorithm package

给：

- objective；
- pseudocode；
- rank/filter selection；
- frequency gate；
- rank/observability event handling；
- acquisition trigger；
- complexity；
- initialization；
- stopping and failure certificate；
- fixed-stratum convergence theorem；
- 与 direct joint inversion 的公平比较协议。

### 5. Manuscript-ready English replacement

直接给可替换到论文中的：

- title；
- abstract；
- introduction；
- related work/originality firewall；
- theory sections；
- method section；
- theorem statements；
- limitations；
- conclusion。

不要把实验计划写成实验结果。

### 6. Claim ledger

每条 claim 标成：

- [proved]
- [conditionally proved]
- [prior art]
- [numerically supported by old loop]
- [falsified by old loop]
- [open]

### 7. Theory-validation experiment blueprint

只设计，不运行。必须覆盖：

1. coherent versus phaseless matched information；
2. rank monotonicity/equality；
3. dual \(\rho/B_{\mathrm{vis}}\) predictive power；
4. PRASC-SOM versus direct joint and phaseless baselines；
5. stacked acquisition/symmetry/gauge。

每个实验给 hypothesis、null、data/model、controls、metrics、seeds、plots、成功/
失败判据和最小复现实验。对 experiment 4 必须预注册 primary endpoint、固定
forward/adjoint budget、由 wavelength/resolution/noise 决定的 pose/map
noninferiority margins、seed/grid budget、confidence interval 与 multiplicity
规则；加入 limited-aperture false-basin cycle-skipping control，以及 clock drift、
mutual coupling 或 unmodeled phase noise 这类“phase corruption 并非 pose 主导”
的 expected-loss control。

### 8. Citation audit

只使用可核验 primary sources。每条科研归属给证据位置。不能把搜索摘要当完整
证明来源。

### 9. Unresolved hard-problem list

只保留本次确实无法解决的问题。对每个问题说明：

- 缺少什么；
- 为什么当前无法完成；
- 最小下一步；
- 它是否阻塞主 claim。

## 九、质量红线

- 不得把 pseudoinverse、projector、rank-nullity、Schur complement、Fisher、
  principal angle 本身包装成创新。
- 不得把 unrestricted current equivalence 写成正贡献。
- 不得声称 phaseless SOM 没有 current-space constraint。
- 不得把 joint pose + image 或 frequency continuation 单独写成创新。
- 不得混淆 \(J_X\) 与有限扰动 \(\Delta X\)。
- 不得混淆 current-space SOM 与 map--pose data-space analysis。
- 不得把 \(\rho\) 与 \(B_{\mathrm{vis}}\) 写成同一个谱。
- 不得未经证明声称 \(\rho(r)\) 随 rank 单调。
- 不得把 \(\rho\) threshold crossing 称作 rank event。
- 不得把 source symmetry null 称作 global gauge。
- 不得把当前 truncated sensing SVD 实现称作完整 Twofold SOM。
- 不得把旧 \(T_U\)、state penalty 或 visible-coordinate solver 写成成功方法。
- 不得用 optimizer success/status 代替物理 residual、pose/map error、basin、
  oracle 和 matched baseline。
- 不得使用 ground truth 选择 rank、frequency、threshold 或 acquisition。
- 不得把 2D scalar scattering-potential frequency sweep 说成 fixed-material
  multifrequency Maxwell result。
- 不得宣称全球 novelty；所有 novelty 判断必须附检索边界。
- 如果 PRASC-SOM 无法证明或设计出任何 SOM-specific 优势，请直接说明，并把
  论文降级为 theory/diagnostic result，而不是生成车轱辘话。

请现在直接开始完成上述研究包。优先把 P1--P5 做到可审稿的数学严谨度，再处理
P6--P7；不要只给路线图，也不要建议我另找人完成证明。

---
