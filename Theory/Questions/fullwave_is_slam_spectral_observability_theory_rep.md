# 位姿不确定性下全波逆散射 SLAM 的谱可观测性

## 理论优先、反例优先、语义严格的定理稿与研究审计

**版本日期：2026-09-04**  
**对象：二维标量 Helmholtz contrast-source 模型；必要处给出 Hilbert 空间、有限维参数化与三维推广边界。**

---

## 0. 证据边界、符号与不可越界的语义

本文建立在以下三个层次上，三者不互相冒充：

1. **已锁定的有限维代数**：白化并 realify 后的实矩阵 $A,B$，以及 $K_{\rm IS},K_{\rm SLAM},K_{\rm eff}$ 的 kernel、rank、principal-angle、interlacing、Born rank-stratum 等恒等式。
2. **本文给出的新推导**：Hilbert 空间闭包版本、Friedrichs-angle 稳定横截、$SE(d)$ gauge、半正定先验的短接算子、共享频率 genericity、rank-event 与统计语义等。每条结论都明确假设。
3. **仍未完成的物理专门化**：真实阵列下的非退化 witness、连续 full-wave 三阶微分统一常数、有限孔径的 scene-independent 正角下界、全局非线性恢复与硬件验证。

固定术语如下。

- $j$：contrast source / induced current。
- $J_X$：pose-prior information；绝不与 $j$ 混写。
- $A=D_\chi F$：**map-tangent-to-data** 算子。
- $B=D_XF$：白化、realify 后由**实位姿自由度**生成的 pose-tangent-to-data 算子。
- $G_S$：Chen SOM 的 **current-to-data** 算子。一般 $A=G_ST_\chi$，但 $A\neq G_S$。
- $\delta X$：内层联合估计中的 nuisance pose increment。
- $\Delta X$：外层轨迹执行误差。两者不得重复计数。
- 本文的 Fisher/Schur 结论只描述名义点附近的 local tangent information geometry；不推出全局唯一性、cycle-skipping 消失、Gauss--Newton 收敛或在线 SLAM 成功。

设实 Hilbert 空间

$$
\mathcal H_\chi,\qquad \mathcal H_X,\qquad \mathcal Y,
$$

以及有界线性算子

$$
A:\mathcal H_\chi\to\mathcal Y,\qquad
B:\mathcal H_X\to\mathcal Y.
$$

记

$$
K_{\rm IS}=A^*A,
\qquad
N:=\overline{\operatorname{Ran}B},
\qquad
K_{\rm SLAM}=A^*P_{N^\perp}A.
$$

当 $J_X\succeq0$ 时，最安全的无限维定义不是先写伪逆，而是先定义二次型

$$
q_{\rm eff}(u)
:=\inf_{h\in\mathcal H_X}
\bigl(\|Au-Bh\|^2+\|J_X^{1/2}h\|^2\bigr),
$$

再由 Riesz 表示得到 $K_{\rm eff}$。

---

# 1. 一页结论摘要

## 1.1 可以作为论文主结论的命题

1. **无限维 retention 有合法定义，但必须通过 polar decomposition，而不是无界的 $K_{\rm IS}^{-1/2}$。** 若 $A=U_A|A|$，在
   $\mathcal H_A:=\overline{\operatorname{Ran}A^*}=(\ker A)^\perp$ 上定义

   $$
   R_{\rm geom}:=U_A^*P_{\overline{\operatorname{Ran}B}^{\perp}}U_A.
   $$

   则 $0\preceq R_{\rm geom}\preceq I$，且

   $$
   K_{\rm SLAM}=|A|R_{\rm geom}|A|,
   \qquad
   \frac{\langle K_{\rm SLAM}u,u\rangle}
        {\langle K_{\rm IS}u,u\rangle}
   =\frac{\langle R_{\rm geom}|A|u,|A|u\rangle}{\||A|u\|^2}.
   $$

2. **无限维无先验消元看到的是 $\overline{\operatorname{Ran}B}$。** 若 $\operatorname{Ran}B$ 不闭，$B^\dagger$ 是无界算子，$BB^\dagger$ 不能作为全空间上的有界公式使用。无惩罚最小二乘的零残差可能只是由一列 $\|h_n\|\to\infty$ 的 pose increments 渐近实现。故必须区分

   $$
   Au=Bh
   \quad\text{与}\quad
   Au\in\overline{\operatorname{Ran}B}.
   $$

3. **“交集为零”远弱于“稳定横截”。** 令

   $$
   M:=\overline{\operatorname{Ran}A},\qquad
   \gamma:=\inf_{m\in M,\ \|m\|=1}\|P_{N^\perp}m\|.
   $$

   则 $K_{\rm SLAM}\succeq\gamma^2K_{\rm IS}$，且 $\gamma^2$ 是最优常数。若 $M\cap N=\{0\}$，则

   $$
   \gamma=\sin\theta_F(M,N)>0
   \iff M+N\text{ 闭}.
   $$

   本文给出一个 $M\cap N=\{0\}$、$A$ 紧且非闭值域、但 $\gamma=0$ 的显式 $\ell^2$ 反例；其 $n$ 维截断最小角约为 $1/n$，每个网格都“正”，极限却不稳定。

4. **对连续 Fourier/Born 数据，全空间正角通常根本不存在。** 对开放目标域和连续有限频带/孔径，restricted Fourier operator 常为紧算子且值域稠密；若 $B\neq0$，则 $M=\mathcal Y$，从而全空间 $\gamma=0$。对有限个复采样点和足够丰富的函数地图空间，$A$ 甚至通常满行，从而 $\operatorname{Ran}B\subseteq\operatorname{Ran}A$。因此 aperture、bandwidth、MIMO 数量和 standoff 单独不能给出 scene-、basis-independent 的正主角下界。可证明的正界必须限定到 task subspace、有限维 basis、正则化支撑或 frame-bounded model class。

5. **均匀无限背景中的 global rigid motion 是精确 gauge。** 对同时作用于地图和整条 body-fixed Tx/Rx 轨迹的 $g\in SE(d)$，

   $$
   F(g\cdot\chi,g\cdot X)=F(\chi,X).
   $$

   沿 Lie algebra 求导得到

   $$
   A\xi_\chi+B\xi_X=0.
   $$

   因而在无锚定的 pose prior 下，gauge map tangent 位于 Schur map information 的核中。径向对称场景的旋转 map generator 可等于零，此时产生的是 pose-only stabilizer，而不是非零 map gauge direction。

6. **共享多频补偿的 genericity 可以在有限维、解析依赖和存在一个非退化 witness 的条件下证明。** 若去除固定 gauge 后，存在一个频率组使 stacked matrix $[A_\Omega,-B_\Omega]$ 满列秩，则几乎所有同样大小的频率组也满列秩；异常集合包含在非零实解析 minor 的零集内。没有 witness、数据维数不足、频率重复、尺度复制、窄带退化或场景对称时，命题不成立。

7. **绝对信息单调、归一化 retention 不单调。** 对共享 pose $h$，新增频率块只会增加

   $$
   u^*K_{\rm eff}u
   =\min_h\sum_f\|A_fu-B_fh\|^2+h^*J_Xh,
   $$

   但分母 $u^*K_{\rm IS}u$ 同时改变，所以 $\rho_i$ 无坐标逐项单调性。本文给出一维反例。精确 shared-Schur 频率选择一般也不具 submodularity；两个各自完全可补偿、合起来互相矛盾的频率块构成最小反例。可用“逐频独立消元”的 additive PSD 下界构造 log-det submodular surrogate，但它只是保守下界。

8. **半正定 pose prior 的正确无限维对象是 augmented-space projection / shorted operator。** 令

   $$
   \bar A u=(Au,0),\qquad
   \bar B h=(Bh,J_X^{1/2}h).
   $$

   则

   $$
   q_{\rm eff}(u)=\|P_{\overline{\operatorname{Ran}\bar B}^{\perp}}\bar Au\|^2.
   $$

   这一定义不要求 $B^*B+J_X$ 闭值域。只有在相应 range condition 与闭值域成立时，才能安全写成 Moore--Penrose Schur 公式。有限 prior 的 retention 仍可写成**增广空间**中的 $\sin^2$ principal angles，但一般不能写成原数据空间两个 ordinary subspaces 的 $\sin^2$。

9. **固定秩 stratum 上 projector 平滑；rank event 处无先验 projector 通常范数跳变 1。** 若 $\operatorname{rank}B(t)=r$，则

   $$
   \dot P_B=(I-P_B)\dot BB^\dagger+\bigl((I-P_B)\dot BB^\dagger\bigr)^*,
   \qquad
   \|\dot P_B\|\le
   \frac{\|\dot B\|}{\sigma_r(B)}.
   $$

   不同有限秩的正交投影之差范数等于 1。因此 hard no-prior $P_B$ 在 rank drop 处不可能 Hölder 连续；$K_{\rm SLAM}$ 仅在 $A^*(P_+-P_0)A=0$ 等额外抵消条件下仍可连续。$J_X\succeq\alpha I$ 会把该 rank singularity 正则化为光滑矩阵逆。

10. **真正可计算的 full-wave 鲁棒性常数需要三类 margin：**

    $$
    \|M^{-1}\|\le m_0^{-1},\qquad
    \text{standoff}\ge d_0,
    \qquad
    \sigma_{\min}^+(B)\ge\beta_0
    $$

    或 $J_X\succeq\alpha I$。常数分别以 $m_0^{-1}$ 的高次幂、$d_0^{-1}$ 和 $\beta_0^{-1}$ 爆炸。P7 中的 $D^2B$ 实际是前向映射 $F$ 的三阶 pose 导数；只假设 Green 函数一、二阶导数不足以完成该目标。

11. **robust trajectory 的一阶展开可严格证明，但仅在 fixed-rank、positive-gap 区域。** 若 $\lambda_r$ 简单且 Hessian 有界，

    $$
    \min_{\|\Delta X\|\le\varepsilon}
    \lambda_r(X+\Delta X)
    =\lambda_r(X)-\varepsilon\|\nabla\lambda_r(X)\|_*+R_2,
    \quad |R_2|\le\frac{L_r}{2}\varepsilon^2.
    $$

    重根的一阶分裂由压缩矩阵给出；对凸 uncertainty set 最小化 $\lambda_{\min}$ 一般是 concave minimization，通常非凸。rank event 附近必须缩小 trust region、引入严格正 prior，或直接使用 bounded-nuisance nonsmooth formulation。

12. **Chen SOM 与 map retention 可以严格桥接，但不能同一化。** $A=G_ST_\chi$ 自动给出 $\operatorname{Ran}A\subseteq\operatorname{Ran}G_S$；等号当且仅当

    $$
    \mathcal H_J=\operatorname{Ran}T_\chi+\ker G_S.
    $$

    $G_S$ 的 right singular vectors 位于 current space，而 $A$ 的 right singular vectors位于 map space；在没有 $T_\chi$ 的识别时，两者的 Grassmann distance 甚至没有定义。先做固定 current reduction 再做 pose elimination是合法的；但“重新计算 dominant modes”的两个次序只有在相应 normal operators 有共同 reducing subspace并保持谱排序时才可交换。

13. **$K_{\rm eff}^{-1}$ 的统计语义有限。** 固定真实 pose、只重复采样 data noise，而优化中加入 deterministic pose penalty 时，map covariance 是

    $$
    K_{\rm eff}^{-1}A^*W^2A K_{\rm eff}^{-1},
    \qquad W=I-B(B^*B+J_X)^{-1}B^*,
    $$

    即 sandwich covariance；一般不是 $K_{\rm eff}^{-1}$。后者在正确的 Gaussian random-pose marginal model、或把产生 $J_X$ 的辅助 pose measurements 也纳入重复采样时成立。即使 Fisher 很大，$F(x)=e^{ikx}$ 仍有周期 aliases，故 local spectrum 不能排除 cycle skipping。

## 1.2 被反例否定或必须降级的命题

- “$\operatorname{Ran}A\cap\operatorname{Ran}B=\{0\}$ 就保证稳定”：**假**。
- “增加带宽/频率必然提高每个归一化 retention”：**假**。
- “精确 shared-Schur 频率选择一般是 submodular”：**假**。
- “孔径、MIMO 数量、standoff 单独给出 continuum 正最小角”：**一般假**。
- “rank-changing projector 仍 Hölder 连续”：**假**。
- “$T_\chi$ 条件数好就保证 SOM current modes 接近 map modes”：**假**；还需 range alignment 与谱隙。
- “先做 SOM 再 pose projection 与反序总可交换”：**假**。
- “固定 nuisance 重复采样时 $K_{\rm eff}^{-1}$ 是 covariance”：**一般假**。
- “局部 Fisher 谱好就保证全局恢复”：**假**。

## 1.3 仍需实质工作的问题

- 为具体 full-wave Tx/Rx 阵列构造一个满足 multi-frequency genericity 的解析 witness；
- 在真实连续 Helmholtz 函数类上证明任务限制后的 uniform frame/inf--sup bound；
- 对 moving-grid/self-cell 情形建立 $D^3F$ 的 Sobolev 映射界；
- 将 quotient gauge、motion-graph sparsity 与 full-wave solver 在同一连续模型中闭合；
- 给出非局部 basin/cycle-skipping 条件或明确的全局反例族；
- 完成投稿级系统查重。本次检索只能支持 retrieval-bounded 的窄贡献定位。

---

# 2. Feasibility triage 与 P1--P10 状态

| 包 | 状态 | 本文完成度 | 核心限制 |
|---|---|---|---|
| P1 Hilbert / stable transversality | **solved + counterexample + conditional discretization** | polar retention、闭包、Friedrichs angle、闭和空间、最小反例、task-window 离散收敛均给出 | 全谱底端无 uniform inf--sup 时不能收敛 |
| P2 $SE(2)/SE(3)$ gauge | **solved for scalar homogeneous model; conditional for vector/constraints** | 群作用、不变性、Lie tangent、stabilizer、quotient FIM、破缺分类 | vector Maxwell 需明确极化/传感器表示；非均匀背景仅保留其 symmetry subgroup |
| P3 Born Fourier / missing wedge | **partial theorem + impossibility result + finite-dimensional bound** | 显式 generators、有限采样满行、连续值域稠密、geometry-only 正角不可能、frame/coherence 条件界 | line/circle/arc 的 exact intersection scene-dependent，无统一闭式 |
| P4 多频共享补偿 | **solved in finite dimension under analytic witness; counterexamples otherwise** | 维数公式、genericity、重复/尺度/窄带反例、单调与非单调、submodularity 反例及 surrogate | 阵列非退化 witness 尚未为具体系统证明 |
| P5 半正定 prior | **solved** | range condition、shorted form、kernel/range split、Loewner、等号、增广角、相关 prior | 无限维 Moore--Penrose 公式仍需闭值域；一般应保留 variational definition |
| P6 fixed-rank / rank event | **solved for matrices; conditional operator version** | smoothness、导数界、全 $\dot K$、重根压缩、projector bound、rank jump、detector | infinite-rank operator rank strata 与数据误差界需另设拓扑和 gap |
| P7 resolvent robustness | **conditional; strongest correct revision given** | resolvent 常数链、Hankel 外部界、内部奇性分类、$DK,D^2K$ 常数 | $D^2B=D^3F$ 需三阶 Green/geometry；continuum self-interaction 不可用点态统一界 |
| P8 robust trajectory | **solved locally; nonsmooth/global open** | 一阶/二阶 remainder、constraint cone、重根、trust region、Pareto | rank event 与 global trajectory problem 一般非凸、非光滑 |
| P9 SOM bridge | **solved algebraically + conditional perturbation + algorithm** | range 条件、closure、singular-value 界、mode-distance 条件、不交换反例、合法 two-stage algorithm | Chen TSOM 的 $G_D$ split 与 pose split仍是层级组合，不是三个交换投影 |
| P10 actual error | **solved for local linear/asymptotic semantics; global open** | 三种 covariance、成立条件、LAN/GN local basin、aliasing反例、论文安全表述 | 不提供 global recovery/cycle-skipping certificate |

---

# 3. 假设总表

| 编号 | 假设 | 用于 | 不满足时发生什么 |
|---|---|---|---|
| H1 | $A,B$ 是白化、realify 后的实 Hilbert 空间有界算子 | 全文 | complex pose projector 会扩大不存在的 nuisance range |
| H2 | $N=\overline{\operatorname{Ran}B}$ | P1/P5 | 若直接用非闭 $\operatorname{Ran}B$，正交投影不存在 |
| H3 | $A=U_A|A|$ 的 polar decomposition；工作空间 $\mathcal H_A=(\ker A)^\perp$ | P1 | 避免使用可能无界的 $K_{\rm IS}^{-1/2}$ |
| H4 | $M=\overline{\operatorname{Ran}A}$ 与 $N$ 为闭子空间 | Friedrichs angle | 角度与 closed-sum 定理都针对闭子空间 |
| H5 | 连续 scalar Helmholtz 背景均匀、各向同性；$\chi\in W_c^{1,\infty}$；direct problem 唯一 | P2 | 群不变性可能只剩背景 symmetry subgroup；不可微 map 无 Lie tangent |
| H6 | 传感器位置/方向均 body-fixed，并在 global rigid motion 下共同变换 | P2/P3 | 已知 receiver、anchor 或外部 phase reference 会破 gauge |
| H7 | Born far field：plane-wave incidence、far-field observation、给定 Fourier convention | P3 | 近场轨迹不能直接解释成 Ewald sampling set |
| H8 | finite-dimensional map/pose basis；$A_f,B_f$ 对频率实解析；远离 poles/resonances | P4 genericity | minor 零集论不适用；仅连续依赖不能给 a.e. 结论 |
| H9 | 去除固定 gauge 后存在至少一个 full-rank frequency tuple | P4 genericity | 没有 witness 就不能从“看起来 generic”推出定理 |
| H10 | $J_X$ 为有界 self-adjoint PSD；需要显式伪逆时 $\operatorname{Ran}(B^*B+J_X)$ 闭 | P5 | 只能用 shorted/variational form，不可把无界伪逆当矩阵 |
| H11 | fixed-rank 区间内 $\sigma_{\min}^+(B)\ge\beta_0>0$ | P6/P7 no-prior | projector derivative 爆炸，rank event 处通常跳变 1 |
| H12 | generalized spectrum 所在 support 固定且 $K_{\rm IS}\succeq a_0I$，cluster 外 gap $\delta>0$ | P6/P8 | 单 eigenvector/gradient 无定义或不稳定 |
| H13 | $\|M^{-1}\|\le m_0^{-1}$；参数集合紧；contrast、frequency、domain 有界 | P7 | 接近共振时统一 derivative constant 发散 |
| H14 | sensor-target standoff $d_0>0$ | $G_S$ 的 Hankel 导数 | 外部 Green 点态导数无统一界 |
| H15 | moving-grid/general geometry 若要求 $D^2B$，需 $G_S,G_D,e^{\rm inc}$ 三阶 pose 可微 | P7 | 仅一、二阶 Green 导数不足，因为 $D^2B=D_X^3F$ |
| H16 | robust eigenvalue theorem中 $K(X)$ 为 $C^2$、目标 eigenvalue simple 或 cluster isolated | P8 | rank event/cluster merge时只能用 nonsmooth或压缩算子 |
| H17 | SOM reduction projector固定，且 discarded tangent norm可控 | P9 | adaptive truncation次序不可交换；误差界失效 |
| H18 | local statistical theorem中有限/task维参数、quotient identifiability、truth interior、likelihood正确 | P10 | Fisher covariance可能不是实际 MSE，LAN/GN basin不成立 |
| H19 | 任何离散到连续结论都有 operator-norm/gap convergence 或 uniform frame bound | P1/P3 | 网格上正最小角可能趋零，machine rank不可当 continuum rank |


---

# 4. 定理—证明—反例正文

## 4.1 P1：Hilbert 空间 retention 与稳定横截

### 定理 4.1（polar retention operator）

设 $A\in\mathcal B(\mathcal H_\chi,\mathcal Y)$，其 polar decomposition 为

$$
A=U_A|A|,
\qquad |A|=(A^*A)^{1/2}.
$$

令

$$
\mathcal H_A:=\overline{\operatorname{Ran}|A|}
=\overline{\operatorname{Ran}A^*}
=(\ker A)^\perp,
$$

$$
M:=\overline{\operatorname{Ran}A},
\qquad
N:=\overline{\operatorname{Ran}B}.
$$

则 $U_A:\mathcal H_A\to M$ 是等距满射。定义

$$
\boxed{
R_{\rm geom}:=
U_A^*P_{N^\perp}U_A\big|_{\mathcal H_A}.
}
$$

则：

1. $R_{\rm geom}$ 是 $\mathcal H_A$ 上的 bounded self-adjoint positive contraction：
   $0\preceq R_{\rm geom}\preceq I$；
2. 作为 $\mathcal H_\chi$ 上的 bounded operators，

   $$
   \boxed{
   K_{\rm SLAM}
   =A^*P_{N^\perp}A
   =|A|R_{\rm geom}|A|;
   }
   $$

3. 对任意 $Au\neq0$，令 $x=|A|u\in\operatorname{Ran}|A|\setminus\{0\}$，则

   $$
   \boxed{
   \frac{\langle K_{\rm SLAM}u,u\rangle}
        {\langle K_{\rm IS}u,u\rangle}
   =
   \frac{\langle R_{\rm geom}x,x\rangle}{\|x\|^2}.
   }
   $$

4. 因 $\operatorname{Ran}|A|$ 在 $\mathcal H_A$ 中稠密，所有 generalized Rayleigh quotients 的闭包正是 $R_{\rm geom}$ 的数值域闭包；特别地

   $$
   \inf_{Au\neq0}
   \frac{\langle K_{\rm SLAM}u,u\rangle}
        {\langle K_{\rm IS}u,u\rangle}
   =\inf\sigma(R_{\rm geom}),
   $$

   无需定义无界的 $K_{\rm IS}^{-1/2}$。

#### 证明

Polar decomposition 的标准性质给出：$U_A$ 的 initial space 是
$\overline{\operatorname{Ran}|A|}=\mathcal H_A$，final space 是
$\overline{\operatorname{Ran}A}=M$，且 $U_A$ 在 $\mathcal H_A$ 上为等距满射。

由于 $P_{N^\perp}$ 是 orthogonal projector，

$$
\langle R_{\rm geom}x,x\rangle
=\langle P_{N^\perp}U_Ax,U_Ax\rangle
=\|P_{N^\perp}U_Ax\|^2,
$$

故 $R_{\rm geom}$ self-adjoint、positive，且

$$
0\le \langle R_{\rm geom}x,x\rangle
\le\|U_Ax\|^2=\|x\|^2.
$$

再由 $A=U_A|A|$，

$$
A^*P_{N^\perp}A
=|A|U_A^*P_{N^\perp}U_A|A|
=|A|R_{\rm geom}|A|.
$$

对 $x=|A|u$，分子为 $\langle R_{\rm geom}x,x\rangle$，分母为
$\|Au\|^2=\||A|u\|^2=\|x\|^2$，得到 Rayleigh quotient 公式。
最后，单位球面上连续二次型在稠密子集上的 infimum 等于在闭包上的 infimum；self-adjoint operator 的 Rayleigh infimum等于谱下端。证毕。

### 推论 4.2（核与“精确补偿/渐近补偿”的区别）

在定理 4.1 的条件下，

$$
\ker R_{\rm geom}
=U_A^*(M\cap N),
$$

而

$$
\boxed{
\ker K_{\rm SLAM}
=\{u\in\mathcal H_\chi:Au\in N\}.
}
$$

若 $\operatorname{Ran}B$ 闭，则 $N=\operatorname{Ran}B$，核正好表示存在 $h$ 使 $Au=Bh$。若 $\operatorname{Ran}B$ 不闭，则可能存在

$$
Au\in N\setminus\operatorname{Ran}B,
$$

此时

$$
\inf_h\|Au-Bh\|=0
$$

但没有任何有限 $h$ 精确实现 $Au=Bh$。因此无限维“zero Schur information”表示 **arbitrarily accurate nuisance compensation**，不必表示 attained exact compensation。

#### 证明

$R_{\rm geom}x=0$ 当且仅当其正二次型为零，即
$P_{N^\perp}U_Ax=0$，等价于 $U_Ax\in N$。又 $U_Ax\in M$，故
$x\in U_A^*(M\cap N)$。第二式直接来自

$$
\langle K_{\rm SLAM}u,u\rangle
=\|P_{N^\perp}Au\|^2.
$$

其余结论由距离到非闭子空间的 infimum 定义得到。证毕。

### 命题 4.3（非闭 pose range 时的正则化极限）

对任意 bounded $B$，令

$$
P_\alpha:=B(B^*B+\alpha I)^{-1}B^*,\qquad \alpha>0.
$$

则

$$
\boxed{
P_\alpha\xrightarrow[\alpha\downarrow0]{\rm strong}
P_{\overline{\operatorname{Ran}B}}.
}
$$

该收敛为 operator-norm convergence 当且仅当 $\operatorname{Ran}B$ 闭，等价于 $B$ 的 reduced minimum modulus 为正，亦即 $0$ 与 $B^*B$ 的正谱分离。

#### 证明

由 polar decomposition 或恒等式

$$
B(B^*B+\alpha I)^{-1}B^*
=BB^*(BB^*+\alpha I)^{-1},
$$

并对正算子 $BB^*$ 使用 spectral calculus。标量函数

$$
f_\alpha(t)=\frac{t}{t+\alpha}
$$

在 $t>0$ 时趋于 1，在 $t=0$ 时等于 0，且 $0\le f_\alpha\le1$。由 dominated convergence 得 strong convergence 到 $BB^*$ 的 support projection，即
$P_{\overline{\operatorname{Ran}B}}$。若正谱可任意接近 0，则
$\sup_{t\in\sigma(BB^*)\setminus\{0\}}|1-f_\alpha(t)|=1$，不可能范数收敛；反之若正谱下界为 $c>0$，则误差至多 $\alpha/(c+\alpha)\to0$。闭值域与正 reduced minimum modulus 等价。证毕。

### 定理 4.4（稳定横截、Friedrichs angle 与 form inequality）

设 $M,N$ 如上，并定义

$$
\gamma
:=\inf_{Au\neq0}
\frac{\|P_{N^\perp}Au\|}{\|Au\|}
=\inf_{m\in M,\ \|m\|=1}\operatorname{dist}(m,N).
$$

则：

1. $\gamma^2$ 是使

   $$
   \boxed{
   K_{\rm SLAM}\succeq cK_{\rm IS}
   }
   $$

   成立的最大常数 $c$；
2. 若 $M\cap N\neq\{0\}$，则 $\gamma=0$；
3. 若 $M\cap N=\{0\}$，则

   $$
   \boxed{
   \gamma^2=1-c_0(M,N)^2=
   \sin^2\theta_F(M,N),
   }
   $$

   其中 $c_0$ 是 minimal-angle cosine，而此时与 Friedrichs cosine 相同；
4. 在 $M\cap N=\{0\}$ 时，以下条件等价：

   $$
   \gamma>0,
   \qquad c_0(M,N)<1,
   \qquad M+N\text{ closed}.
   $$

5. 若 $C=M\cap N$ 非零，定义 quotient/transverse constant

   $$
   \gamma_F:=
   \inf_{m\in M\cap C^\perp,\ \|m\|=1}
   \operatorname{dist}(m,N),
   $$

   则 $\gamma_F=\sin\theta_F(M,N)$，且
   $\gamma_F>0\iff M+N$ closed。原始 $\gamma$ 仍为零，因为 exact intersection 已被包含。

#### 证明

$\operatorname{Ran}A$ 在 $M$ 中稠密，且单位向量可由归一化后的 range vectors 逼近，故两种 infimum 相等。对任意 $u$，

$$
\langle K_{\rm SLAM}u,u\rangle
=\|P_{N^\perp}Au\|^2
\ge\gamma^2\|Au\|^2
=\gamma^2\langle K_{\rm IS}u,u\rangle.
$$

若用更大的常数 $c>\gamma^2$，按 infimum 定义可取一列 $u_n$ 使 quotient 趋于 $\gamma^2$，故 inequality 失败；所以 $\gamma^2$ 最优。

若 $M\cap N$ 含单位向量 $m$，则距离为零。若交集为零，对单位 $m\in M$，

$$
\operatorname{dist}(m,N)^2
=1-\|P_Nm\|^2.
$$

因此

$$
\gamma^2
=1-\sup_{\|m\|=1,m\in M}\|P_Nm\|^2
=1-\|P_NP_M\|^2
=1-c_0(M,N)^2.
$$

Friedrichs-angle closed-sum theorem给出
$c(M,N)<1\iff M+N$ closed；在交集为零时 $c=c_0$。一般交集情形在 $C^\perp$ 上应用同一论证，得到第 5 项。证毕。

### 反例 4.5（交集为零但 $\gamma=0$；有限截断虚假稳定）

取

$$
\mathcal Y=\ell^2\oplus\ell^2,
\qquad
M=\ell^2\oplus\{0\},
$$

令 $T:\ell^2\to\ell^2$ 为

$$
Te_n=\frac1n e_n,
$$

并令

$$
N=\operatorname{graph}(T)
=\{(x,Tx):x\in\ell^2\}.
$$

$T$ bounded 且 injective，故 $N$ 是闭子空间并且

$$
M\cap N=\{0\}.
$$

定义 compact injective map operator

$$
Ae_n=\left(\frac1n e_n,0\right),
$$

以及 pose operator

$$
Bx=(x,Tx).
$$

则 $\overline{\operatorname{Ran}A}=M$，$\operatorname{Ran}B=N$ 闭，但

$$
\operatorname{dist}((e_n,0),N)
\le\left\|(e_n,0)-(e_n,n^{-1}e_n)\right\|
=\frac1n,
$$

所以

$$
\boxed{
M\cap N=\{0\}
\quad\text{而}\quad
\gamma=0.
}
$$

在前 $n$ 个坐标截断中，最小 principal angle 满足

$$
\sin\theta_{\min}^{(n)}
=\frac{1/n}{\sqrt{1+1/n^2}}>0,
$$

故每一个有限网格都报告正角，但它趋于零。此例同时证明：

- exact intersection test 不能替代 stable transversality；
- finite-grid $\theta_{\min}>0$ 本身没有 continuum 含义；
- compact $A$ 的小 singular directions 可以逐渐贴近 pose range，而不在任何有限层级精确相交。

### 推论 4.6（dense map-data range 对全空间稳定横截的否定）

若

$$
\overline{\operatorname{Ran}A}=\mathcal Y
$$

且 $B\neq0$，则 $N\neq\{0\}$ 且

$$
M\cap N=N,
\qquad
\gamma=0.
$$

因此对 dense-range compact inverse operator，除非先 quotient 掉 $N$、限制 map task space、加入 map regularity/prior，或改变数据空间 metric，否则不存在全空间的正 stable-transversality constant。

### 定理 4.7（有限维离散逼近的可证明版本）

设 $P_n:\mathcal H_\chi\to\mathcal H_{\chi,n}$、$Q_n:\mathcal Y\to\mathcal Y_n$ 为强收敛到恒等的正交投影，并设
Let $R_n:\mathcal H_X\to\mathcal H_{X,n}$ be the corresponding strongly convergent orthogonal pose projection.

$$
A_n=Q_nAP_n,
\qquad
N_n=\overline{\operatorname{Ran}(Q_nB R_n)}
$$

为离散 map 与 pose spaces。固定 $\tau>0$，令

$$
E_\tau:=\mathbf 1_{[\tau,\infty)}(|A|)\mathcal H_A.
$$

若：

1. $A$ compact，且离散族 collectively compact 并点态一致，从而
   $\|A_nP_n-A\|\to0$；
2. $\tau$ 与 $|A|$ 的谱分离，故 $E_\tau$ 有限维；
3. pose subspaces gap-converge：

   $$
   \|P_{N_n}-P_N\|\to0;
   $$

4. $E_{\tau,n}$ 是 $|A_n|$ 对应 isolated cluster 的 spectral subspace，按 Kato/Davis--Kahan 有
   $\|P_{E_{\tau,n}}-P_{E_\tau}\|\to0$；

则把 $R_{{\rm geom},n}$ transport 到 $E_\tau$ 后有 operator-norm convergence，因而该固定 task window 内的 retention eigenvalues、spectral projectors 与 principal angles 收敛。

反之，若只知道 $A_n\to A$ 点态、或只看扩张 trial spaces 上的最小角，而没有 uniform inf--sup/frame bound

$$
\inf_n\inf_{0\ne u\in E_n}
\frac{\|P_{N_n^\perp}A_nu\|}{\|A_nu\|}>0,
$$

则最底 retention 可趋于零，反例 4.5 已说明这一点。

#### 证明要点

由 collective compactness 加点态收敛得到 compact-operator 的 norm convergence；isolated nonzero singular cluster 因谱隙而 gap-converge。$P_{N_n}$ 的 norm convergence 加上 $U_{A_n}$ 在该 cluster 上的 norm convergence，使

$$
U_{A_n}^*P_{N_n^\perp}U_{A_n}
\to U_A^*P_{N^\perp}U_A
$$

在有限维 $E_\tau$ 上按 operator norm 收敛。Hermitian spectral perturbation 随即给出 eigenvalue 与 projector convergence。由于 compact operator 的 singular values 聚于 0，同一论证不能把 $\tau$ 取到 0。证毕。

> **论文含义。** Continuum claim 应改写成“对固定 positive-information task window，在 operator/gap convergence 与 pose-frame bound 下，离散 retention 收敛”；不能写“网格加密验证了 continuum 最小角”。

---

## 4.2 P2：连续 Helmholtz 的 $SE(d)$ gauge 与破缺

### 模型 4.8

令 $d=2$ 或 $3$，背景为均匀各向同性介质，scalar Green kernel 满足

$$
g_k(Qx+a,Qz+a)=g_k(x,z),
\qquad (Q,a)\in SE(d).
$$

取 $\chi\in W_c^{1,\infty}(\mathbb R^d)$。第 $t$ 个平台 pose 为
$X_t=(R_t,p_t)\in SE(d)$，body-fixed Tx/Rx offsets 分别为 $r_a^{\rm T},r_b^{\rm R}$，物理位置为

$$
x_{t,a}^{\rm T}=R_tr_a^{\rm T}+p_t,
\qquad
x_{t,b}^{\rm R}=R_tr_b^{\rm R}+p_t.
$$

假设 Lippmann--Schwinger state operator 在所考虑参数邻域可逆，且 data 以平台 body frame 中不随 global frame 旋转的 scalar 通道表示。

### 定理 4.9（global rigid-motion invariance）

对 $g=(Q,a)\in SE(d)$ 定义

$$
(g\cdot\chi)(x):=\chi(g^{-1}x)
=\chi(Q^T(x-a)),
$$

$$
g\cdot X_t:=(QR_t,Qp_t+a).
$$

若 incident field 与 receiver sampling 都由上述 body-fixed physical locations 产生，则

$$
\boxed{
F(g\cdot\chi,g\cdot X)=F(\chi,X).
}
$$

#### 证明

以 point-source incidence 为例，名义总场满足

$$
u(x)=u^{\rm inc}(x;X)+
\int g_k(x,z)\chi(z)u(z)\,dz.
$$

定义 $u_g(x):=u(g^{-1}x)$。将 $z=gw$ 代入，并使用 rigid motion 的 Jacobian 为 1、Green kernel 的 Euclidean invariance，以及 transformed Tx 位置，得到

$$
u_g(x)=u_g^{\rm inc}(x;gX)+
\int g_k(x,z)(g\cdot\chi)(z)u_g(z)\,dz.
$$

由 direct problem 唯一性，$u_g$ 就是 transformed configuration 的总场。对 transformed receiver $gx_{t,b}^{\rm R}$，

$$
u_g(gx_{t,b}^{\rm R})=u(x_{t,b}^{\rm R}).
$$

逐通道成立，故 stacked data 不变。plane wave 的方向和 phase origin若也随 body frame 共同旋转/平移，同理。证毕。

### 推论 4.10（Lie algebra gauge tangent）

令 $g(s)=\exp(s\xi)$，其中

$$
\xi=(v,\Omega),
\qquad \Omega^T=-\Omega.
$$

群作用在 map 上的 infinitesimal generator 为

$$
\boxed{
\xi_\chi(x)=-(v+\Omega x)\cdot\nabla\chi(x).
}
$$

pose generator 为

$$
\boxed{
\delta p_t=v+\Omega p_t,
\qquad
\delta R_t=\Omega R_t.
}
$$

若 $F$ Fréchet differentiable，则

$$
\boxed{
A\xi_\chi+B\xi_X=0.
}
$$

因此无 pose prior 时

$$
\xi_\chi\in\ker K_{\rm SLAM}.
$$

若有 prior，则 joint quadratic form沿该 gauge pair 等于

$$
\|A\xi_\chi+B\xi_X\|^2
+\langle J_X\xi_X,\xi_X\rangle
=\langle J_X\xi_X,\xi_X\rangle.
$$

所以该 map gauge 仍在 $K_{\rm eff}$ 核中当且仅当相应 global pose motion 位于 $\ker J_X$；固定 anchor 会把它抬出核。

#### 证明

对定理 4.9 的恒等式关于 $s$ 在 0 处求导并使用 chain rule。map generator由
$g(s)^{-1}x=x-s(v+\Omega x)+o(s)$ 得到。证毕。

### 命题 4.11（径向对称场景的 pose-only stabilizer）

若

$$
\chi(x)=\varphi(\|x-c\|)
$$

且考虑绕 $c$ 的旋转 generator，则

$$
\xi_\chi(x)=-(\Omega(x-c))\cdot\nabla\chi(x)=0.
$$

由 gauge identity 得

$$
B\xi_X=0.
$$

这表示同一场景对该 global rotation 有非平凡 isotropy/stabilizer。它是 pose-only null direction，不应被计数为一个非零 map gauge mode。此时 quotient 不是自由群作用下的普通流形；局部上应按 orbit type 视为 stratified quotient/orbifold。

### 命题 4.12（gauge 的破坏与近似破坏）

以下结论逐项成立。

1. **有限成像窗。** 若窗只是计算截断，而 $\operatorname{supp}\chi$ 与窗边界有正距离，则足够小的 rigid motion 仍在 admissible class 内，local gauge 不被破坏。只有固定物理边界、已知窗外材料或 support 接触边界时，群轨道才离开 admissible set。
2. **已知非均匀背景。** 若背景系数 $b(x)$ 固定在 world frame，则 gauge 只保留满足 $b(gx)=b(x)$ 的 symmetry subgroup；generic background 通常破坏全部 global translations/rotations。
3. **固定 anchor / 已知 receiver。** 若某个 pose、receiver 或 beacon不随 $g$ 变换，则相对几何改变，invariance 被破坏。
4. **外部相位/时钟 reference。** 若 data phase 以固定 world reference 定义，global motion可能留下可测 phase，必须把该 reference 写入 forward model；不能沿用无 reference 的 gauge。
5. **离散 pixel basis。** 连续 generator $-(v+\Omega x)\cdot\nabla\chi$ 一般不在 piecewise-constant pixel space 内，故离散 gauge residual 不为零。这是 representation error，不是 continuum gauge 被物理破坏。
6. **support constraint。** “support 已知”只有在它固定了绝对位置并且不容许 rigid displacement 时才破 gauge；仅知道形状大小或一个宽松包络不一定破坏 local gauge。

### 定理 4.13（quotient parameter manifold 与 quotient FIM）

设联合参数流形

$$
\mathcal Z=\mathcal M_\chi\times\mathcal X,
$$

$G=SE(d)$ 在名义点 $z=(\chi,X)$ 邻域自由且 proper 地作用。令 vertical orbit tangent

$$
\mathcal V_z
=\{(\xi_\chi,\xi_X):\xi\in\mathfrak{se}(d)\},
$$

选定参数空间 metric 后令

$$
\mathcal H_z=\mathcal V_z^\perp.
$$

若 joint information operator 为 $\mathcal I_z\succeq0$，则 quotient
$\mathcal Z/G$ 上的 local information由 horizontal compression 表示：取任意 horizontal isometry
$H_z:\mathbb R^r\to\mathcal H_z$，

$$
\boxed{
\mathcal I_{\rm quot}=H_z^*\mathcal I_zH_z.
}
$$

其谱与 horizontal basis 的选择无关。若群作用有 stabilizer，应按实际 orbit tangent 取 $\mathcal V_z$，并在固定 orbit-type stratum 上使用同一构造。

#### 证明

$DF(z)$ 在 vertical directions 上为零，因此它下降为 quotient tangent 上的线性映射。horizontal lift 给出 quotient tangent 与 $\mathcal H_z$ 的等距同构；Fisher 二次型在该同构下就是上述 compression。换另一个 horizontal orthonormal basis只产生 orthogonal similarity。证毕。

### Remark 4.14（motion-graph gauge 与 map gauge 的关系）

未锚定 odometry/pose-graph Laplacian 的 kernel 是全局 rigid motion。full-wave joint model中的 map gauge不是另一个独立可相加的 gauge；它与该 pose global mode组成同一个 diagonal group orbit

$$
(\xi_\chi,\xi_X).
$$

若只 quotient pose graph 而把 map 固定在 world coordinates，或只固定 map而不给 pose anchor，会造成不一致的 gauge bookkeeping。正确做法是对 joint state quotient，或明确以 anchor/known background定义 world frame。


---

## 4.3 P3：Born 远场 Fourier 几何、missing wedge 与角界

### 约定 4.15（plane-wave/far-field Born model）

采用 Fourier convention

$$
\widehat\chi(q)=\int_D\chi(x)e^{-iq\cdot x}\,dx.
$$

入射 plane wave 的传播方向为 $d\in\mathbb S^{d-1}$，far-field observation方向为 $s\in\mathbb S^{d-1}$。令 Tx/Rx phase origins 为 $p_T,p_R$。在一个固定的标量归一化下，Born datum 写为

$$
\boxed{
y(k,s,d;p_T,p_R)
=C(k,s,d)
 e^{-ikd\cdot p_T}
 e^{iks\cdot p_R}
 \widehat\chi(q),
\qquad q=k(s-d).
}
$$

$C$ 吸收 Green far-field normalization、极化与已知幅度；以下相位符号均相对于此约定。换 Fourier convention 会同时改变两个符号，但 gauge cancellation 不变。

### 定理 4.16（地图与 Tx/Rx 位姿 generators 的显式关系）

对 infinitesimal Tx/Rx translations $a_T,a_R$，

$$
D_{p_T}y[a_T]=-ik(d\cdot a_T)y,
\qquad
D_{p_R}y[a_R]=ik(s\cdot a_R)y.
$$

若两者共同平移 $a_T=a_R=a$，

$$
D_Xy[a]=i(q\cdot a)y.
$$

地图共同平移的 generator 为

$$
\delta\chi_a=-a\cdot\nabla\chi,
$$

且

$$
A\delta\chi_a=-i(q\cdot a)y.
$$

故

$$
\boxed{
A\delta\chi_a+B\delta X_a=0.
}
$$

对 infinitesimal rotation $\Omega^T=-\Omega$，地图 generator

$$
\delta\chi_\Omega(x)=-(\Omega x)\cdot\nabla\chi(x)
$$

满足

$$
\widehat{\delta\chi_\Omega}(q)
=-(\Omega q)\cdot\nabla_q\widehat\chi(q).
$$

若 Tx/Rx positions 与 directions $d,s$ 一同旋转，则 $q$ 的 generator 为
$\delta q=\Omega q$，phase-center dot products保持不变，因而

$$
D_Xy[\Omega]
=C\,e^{\rm phase}
(\Omega q)\cdot\nabla_q\widehat\chi(q),
$$

再次与 map generator cancellation。

#### 证明

translation 导数直接对两个 exponential factors 求导。对地图平移，

$$
\widehat{\chi(\cdot-a)}(q)
=e^{-iq\cdot a}\widehat\chi(q),
$$

在 $a=0$ 求导得到 $-i(q\cdot a)\widehat\chi$。旋转时

$$
\widehat{\chi(Q^T\cdot)}(q)
=\widehat\chi(Q^Tq),
$$

对 $Q(t)=e^{t\Omega}$ 求导即得。共同旋转使 $q(t)=Q(t)q$，两项符号相反。证毕。

### 命题 4.17（一般 pose tangent 的 Fourier 形式）

若某个 pose coordinate $h$ 同时改变已知 prefactor 与 Fourier sampling location，则

$$
\boxed{
D_Xy[h]
=\bigl(D_X\log C_{\rm tot}[h]\bigr)y
+C_{\rm tot}\,
\nabla_q\widehat\chi(q)\cdot D_Xq[h],
}
$$

其中 $C_{\rm tot}$ 包含 amplitude 与 phase factors。故 map--pose exact compensation 在采样点 $q_\ell$ 上满足线性函数方程

$$
\boxed{
\widehat{\delta\chi}(q_\ell)
=\phi_\ell(h)\widehat\chi(q_\ell)
+\nabla_q\widehat\chi(q_\ell)\cdot\delta q_\ell(h),
\quad\forall\ell,
}
$$

其中符号/已知系数已吸收进 $\phi_\ell$。global translations/rotations 是该方程的显式解；其余解高度依赖 nominal scene $\chi$、采样集合与 map class，不能只由轨迹名称决定。

### 命题 4.18（Ewald/Fourier sampling set）

在 fixed $k$ 下：

1. fixed incidence $d_0$、varying observation $s$ 采样 shifted Ewald sphere/circle

   $$
   \Omega_q=k(\mathbb S^{d-1}-d_0);
   $$

2. full bistatic directions $s,d\in\mathbb S^{d-1}$ 的 difference set 满足

   $$
   \{k(s-d)\}=\{q:\|q\|\le2k\};
   $$

3. backscatter/monostatic idealization $s=-d$ 采样半径 $2k$ 的 sphere/circle；
4. limited angular aperture 给出 arc/sector/missing wedge；bandwidth $k\in[k_{\min},k_{\max}]$ 将单一曲线/曲面增厚为 annular sector 或 ball shell。

这些是 plane-wave/far-field 结论。物理空间中的 straight/circle trajectory只有在其 wavefront directions 可近似映射为上述 $d,s$ 时才能转译为 Fourier aperture；近场 spherical-wave 数据不能直接用“圆轨迹 = full Ewald disk”替代。

### 定理 4.19（有限采样下的 map-space saturation）

设 $D\subset\mathbb R^d$ 含非空开集，取互异采样点 $q_1,\ldots,q_m$，定义

$$
A:L^2(D;\mathbb C)\to\mathbb C^m,
\qquad
(Au)_\ell=\int_Du(x)e^{-iq_\ell\cdot x}\,dx.
$$

则 $A$ 满行，即

$$
\boxed{
\operatorname{Ran}A=\mathbb C^m.
}
$$

因此对任意 finite-dimensional pose tangent $B:\mathbb R^p\to\mathbb C^m$，

$$
\operatorname{Ran}B\subseteq\operatorname{Ran}A,
$$

从而所有 pose data directions 都可由某个 unrestricted $L^2$ map tangent 精确复制。

对 real-valued map 与 realified data，若

$$
\{\Re e^{-iq_\ell\cdot x},\Im e^{-iq_\ell\cdot x}\}_{\ell=1}^m
$$

在 $L^2(D;\mathbb R)$ 中线性独立，则相同的 full-row conclusion 对应的 realified operator 成立；若采样包含共轭对称冗余，应先删除依赖行。

#### 证明

$A^*: \mathbb C^m\to L^2(D)$ 为

$$
A^*c=\sum_{\ell=1}^mc_\ell e^{iq_\ell\cdot x}.
$$

若 $A^*c=0$ 在 $D$ 上几乎处处，则有限 exponential sum 是 real-analytic function，因在开集上为零而处处为零。互异频率的 exponentials 线性独立，故 $c=0$。于是 $A^*$ injective，有限维 codomain 下等价于 $A$ surjective。实数版同理。证毕。

> **重要语义。** 该定理不是说 pose 永远无法估计；它说在“有限数据 + 完全自由的无限维地图 tangent”模型中，地图可吸收任何 pose signature。实际可分离性来自多次共享结构、有限/光滑 map parameterization、先验、额外传感器、非线性高阶约束或任务限制，而不是来自有限采样本身的全空间主角。

### 定理 4.20（连续 aperture 下 restricted Fourier range 稠密）

设 $D,\Omega\subset\mathbb R^d$ 均为 bounded open sets，定义

$$
A:L^2(D)\to L^2(\Omega),
\qquad
(Au)(q)=\int_Du(x)e^{-iq\cdot x}\,dx.
$$

则 $A$ 为 Hilbert--Schmidt compact operator；若 $D$ 与 $\Omega$ 含开集，则

$$
\boxed{
\overline{\operatorname{Ran}A}=L^2(\Omega).
}
$$

因而任何非零 closed pose subspace $N\subset L^2(\Omega)$ 都满足全空间
$\gamma=0$。

#### 证明

kernel $e^{-iq\cdot x}$ 在 $\Omega\times D$ 上平方可积，故 $A$ Hilbert--Schmidt。若 $g\in\ker A^*$，则

$$
(A^*g)(x)=\int_\Omega g(q)e^{iq\cdot x}\,dq=0
\quad\text{a.e. on }D.
$$

因为 $g$ 有 bounded support，右侧是 entire/real-analytic function of $x$；它在开集 $D$ 上为零，故处处为零，Fourier uniqueness 给出 $g=0$。于是
$\ker A^*=\{0\}$，由
$\overline{\operatorname{Ran}A}=(\ker A^*)^\perp$ 得结论。证毕。

### 推论 4.21（不存在 geometry-only continuum 正角下界）

即使 aperture、bandwidth、MIMO count 和 standoff 都固定为“良好”数值，只要 map class 足够丰富使定理 4.19 或 4.20 适用，就不存在仅依赖这些几何量的常数 $c>0$ 使

$$
\|P_{N^\perp}Au\|\ge c\|Au\|
\quad\forall u.
$$

更直接地，global $SE(d)$ gauge 已给出零角；即使 quotient 掉 gauge，finite-data full-row map operator仍使其他 pose directions进入 $\operatorname{Ran}A$。因此 P3.4 的“从孔径、带宽、MIMO、standoff 推出统一最小角”在 unrestricted continuum map 上被否定。

### 定理 4.22（有限维 frame/coherence 正角下界）

设有限维、已去 gauge 的 $A\in\mathbb R^{m\times n}$、$B\in\mathbb R^{m\times p}$ 满足

$$
A^TA\succeq a^2I,
\qquad
B^TB\succeq b^2I,
\qquad
\|A^TB\|\le\mu<ab.
$$

则 $\operatorname{Ran}A\cap\operatorname{Ran}B=\{0\}$，且最小 principal angle 满足

$$
\boxed{
\gamma
\ge
\sqrt{1-\left(\frac{\mu}{ab}\right)^2}.
}
$$

#### 证明

令 $Q_A=A(A^TA)^{-1/2}$、$Q_B=B(B^TB)^{-1/2}$。则

$$
\|Q_A^TQ_B\|
\le\|(A^TA)^{-1/2}\|\,\|A^TB\|\,
\|(B^TB)^{-1/2}\|
\le\frac\mu{ab}<1.
$$

最大 canonical correlation 小于 1，故无 exact intersection，且

$$
\gamma^2=1-\|Q_A^TQ_B\|^2
\ge1-(\mu/ab)^2.
$$

证毕。

> 这一定理把 aperture、bandwidth、MIMO 与 standoff 的作用放在正确位置：它们必须先被证明能给出 map frame lower bound $a$、pose frame lower bound $b$ 与 cross-coherence upper bound $\mu$；几何名称本身不是证明。

### 命题 4.23（missing wedge 与 near-confusion 的最强安全表述）

1. missing wedge 首先产生的是 $A$ 的小 singular values，即 **absolute information loss**；它并不自动等价于 map--pose intersection。
2. 在有限维 Fourier/spline subspace $V_n$ 上，若采样 frame operator满足

   $$
   a_n^2\|u\|^2\le\|A u\|^2\le A_n^2\|u\|^2,
   $$

   且 cross coherence 可控，则定理 4.22 给出正角。
3. 在 narrow-band/small-aperture parameter $\eta$ 下，若 $\widehat\chi$、prefactor 和 sampling map 为 $C^2$，并能在中心频率/角度匹配 pose signature 的零阶与一阶 jet，则 Taylor remainder 给出 data residual

   $$
   \|Au_\eta-Bh\|=O(\eta^2)
   $$

   （相对于固定归一化）。这构成 near-confusion 上界；常数依赖 $\sup\|D^2\widehat\chi\|$、phase-center 与参数ization。
4. 不存在对所有 scenes 的 line/circle/partial-arc exact intersection统一闭式。对某些对称 scene，圆周也有 stabilizer；对某些有限 basis，直线采样也可 full rank。轨迹排序必须通过实际 $A,B$ 或经过证明的 frame constants 给出。

---

## 4.4 P4：多频共享位姿补偿

### 定理 4.24（共享补偿核的精确维数）

设有限维

$$
A_f\in\mathbb R^{m_f\times n},
\qquad
B_f\in\mathbb R^{m_f\times p},
$$

并 vertically stack

$$
A_\Omega=\begin{bmatrix}A_{f_1}\\\vdots\\A_{f_F}\end{bmatrix},
\qquad
B_\Omega=\begin{bmatrix}B_{f_1}\\\vdots\\B_{f_F}\end{bmatrix}.
$$

共享补偿 map kernel 为

$$
\mathcal C_\Omega
:=\{u:\exists h,\ A_fu=B_fh\ \forall f\}
=\pi_u\ker[A_\Omega,-B_\Omega].
$$

则

$$
\boxed{
\dim\mathcal C_\Omega
=n+\operatorname{rank}B_\Omega
-\operatorname{rank}[A_\Omega\ B_\Omega].
}
$$

并且

$$
\boxed{
\dim\bigl(\operatorname{Ran}A_\Omega
\cap\operatorname{Ran}B_\Omega\bigr)
=
\operatorname{rank}A_\Omega+
\operatorname{rank}B_\Omega-
\operatorname{rank}[A_\Omega\ B_\Omega].
}
$$

故

$$
\dim\mathcal C_\Omega
=\dim\ker A_\Omega+
\dim(\operatorname{Ran}A_\Omega\cap\operatorname{Ran}B_\Omega).
$$

#### 证明

令 $L=[A_\Omega,-B_\Omega]$。rank-nullity 给出

$$
\dim\ker L=n+p-\operatorname{rank}L.
$$

从 pair kernel投影到 $u$ 的 kernel 是

$$
\ker(\pi_u|_{\ker L})
=\{0\}\times\ker B_\Omega,
$$

维数为 $p-\operatorname{rank}B_\Omega$。因此

$$
\dim\mathcal C_\Omega
=\dim\ker L-\dim\ker B_\Omega
=n+\operatorname{rank}B_\Omega-
\operatorname{rank}L.
$$

第二式是两个有限维 column spaces 的标准维数公式；第三式代数化简。证毕。

### 定理 4.25（解析频率 genericity：条件版）

设频率参数

$$
\omega=(f_1,\ldots,f_F)\in U\subset\mathbb R^F
$$

位于一个连通开集，所有 entries of $A_\omega,B_\omega$ 对 $\omega$ 实解析，并避开 forward resonances。设存在一个固定 gauge pair subspace

$$
\mathcal G\subseteq\ker[A_\omega,-B_\omega]
\quad\forall\omega,
\qquad \dim\mathcal G=g.
$$

取 domain complement matrix $C$ whose columns span $\mathcal G^\perp$，并假设数据维数满足

$$
\sum_fm_f\ge n+p-g.
$$

若存在一个 witness $\omega_0\in U$ 使

$$
[A_{\omega_0},-B_{\omega_0}]C
$$

满列秩 $n+p-g$，则除一个 Lebesgue-measure-zero、empty-interior 的 real-analytic exceptional set 外，

$$
\boxed{
\ker[A_\omega,-B_\omega]=\mathcal G.
}
$$

也就是说，**在这些附加假设下**，几乎所有有限频率选择都只保留 unavoidable gauge。

#### 证明

witness full column rank意味着存在一个 $(n+p-g)\times(n+p-g)$ minor determinant
$\Delta(\omega)$ 在 $\omega_0$ 非零。$\Delta$ 是实解析且不恒为零。非零实解析函数的零集在连通开集内 Lebesgue measure zero 且无内点。对 $\Delta(\omega)\neq0$，reduced matrix满列秩，所以其 kernel 为零；恢复被 quotient 的 $\mathcal G$ 后得到结论。证毕。

### Remark 4.26（为什么该定理尚未闭合具体物理系统）

“阵列 generic”不是 witness。必须对具体 Tx/Rx geometry、map basis、frequency count 与 scene，给出一个可核验的 $\omega_0$ 或证明某个 analytic minor不恒为零。若 map basis维数随频率数量增长、pose blocks随时间独立增长、或 gauge basis未被离散表示，定理的维数和固定-domain假设都要重新检查。

### 反例 4.27（genericity 命题的最小失败族）

1. **重复频率：** $A_{f_2}=A_{f_1}$、$B_{f_2}=B_{f_1}$，只复制同一方程，不减少共享核。
2. **尺度复制：** 对所有 $f$，

   $$
   A_f=c_fA_0,
   \qquad B_f=c_fB_0,
   $$

   则 $A_fu=B_fh$ 对所有 $f$ 等价于一个频率的方程。
3. **窄带极限：** $f_j\to f_0$ 时 stacked matrix趋于 repeated-block matrix；即使每个非零 bandwidth下 exact rank恢复，smallest singular value可趋零，故 generic exact rank不等于 stable transversality。
4. **对称 scene：** scene stabilizer可使 $B_fh=0$ 对所有频率；额外 map/pose gauge也可在每个频率持续。
5. **维数不足：** 若 $\sum m_f<n+p-g$，reduced full column rank不可能。
6. **Born empty background：** 在 $\chi_0=0$，$B_f=0$；一阶问题没有 pose tangent，多频不能用该一阶 $B_f$ 检验二阶 bilinear pose effect。

### 定理 4.28（共享多频的 absolute-information monotonicity）

对频率集合 $S$ 定义

$$
q_S(u)
:=\min_h\left[
\sum_{f\in S}\|A_fu-B_fh\|^2+h^TJ_Xh
\right]
=u^TK_{\rm eff}(S)u.
$$

若 $S\subseteq T$ 且使用同一个 shared pose variable $h$ 和同一个 prior term，则

$$
\boxed{
K_{\rm eff}(T)\succeq K_{\rm eff}(S).
}
$$

#### 证明

对任意固定 $h$，$T$ 的 objective 等于 $S$ 的 objective 加上非负项
$\sum_{f\in T\setminus S}\|A_fu-B_fh\|^2$。对 $h$ 取 infimum 后不等式仍成立。对所有 $u$ 成立即为 Loewner monotonicity。证毕。

### 反例 4.29（归一化 retention 不单调）

取 scalar map 与 scalar pose。第一块

$$
A_1=[1],\qquad B_1=[0]
$$

给出

$$
K_{\rm IS}^{(1)}=1,
\qquad K_{\rm SLAM}^{(1)}=1,
\qquad\rho^{(1)}=1.
$$

加入第二块

$$
A_2=[1],\qquad B_2=[1].
$$

stack 后

$$
A=\begin{bmatrix}1\\1\end{bmatrix},
\qquad
B=\begin{bmatrix}0\\1\end{bmatrix},
$$

故

$$
K_{\rm IS}^{(1,2)}=2,
\qquad
K_{\rm SLAM}^{(1,2)}=1,
\qquad
\rho^{(1,2)}=\frac12.
$$

absolute effective information未下降，但 normalized retention 从 1 降到 $1/2$。

### 定理 4.30（逐频独立消元给 shared compensation 的 additive 下界）

令 $J_X=\sum_{f\in S}J_f$，$J_f\succeq0$，并定义

$$
K_f^{\rm ind}
:=A_f^TA_f-A_f^TB_f(B_f^TB_f+J_f)^\dagger B_f^TA_f
$$

或用其 variational/shorted definition。则

$$
\boxed{
K_{\rm eff}^{\rm shared}(S)
\succeq\sum_{f\in S}K_f^{\rm ind}.
}
$$

#### 证明

shared compensation只允许一个 $h$：

$$
\min_h\sum_f
\bigl(\|A_fu-B_fh\|^2+h^TJ_fh\bigr)
\ge
\sum_f\min_{h_f}
\bigl(\|A_fu-B_fh_f\|^2+h_f^TJ_fh_f\bigr).
$$

右侧即 $u^T(\sum K_f^{\rm ind})u$。证毕。

### 推论 4.31（可证明的 submodular lower-bound surrogate）

取 map ridge $R_0\succ0$，定义

$$
\Phi_{\rm LB}(S)
=\log\det\left(R_0+
\sum_{f\in S}K_f^{\rm ind}\right).
$$

则 $\Phi_{\rm LB}$ monotone submodular。它给出 shared-information log-det 的一个保守 additive surrogate，但一般不是精确 shared-Schur objective。

#### 证明

PSD 增量使 log-det 单调。若 $S\subseteq T$，记
$M_S=R_0+\sum_{f\in S}K_f$，则 $M_T\succeq M_S$。添加 $K_e\succeq0$ 的 marginal gain为

$$
\log\det(M+K_e)-\log\det M
=\int_0^1\operatorname{tr}
\bigl[(M+tK_e)^{-1}K_e\bigr]dt.
$$

$M_T+tK_e\succeq M_S+tK_e$ 导致 inverse 的 Loewner order反向，trace with $K_e\succeq0$ 因而不增，得到 diminishing returns。证毕。

### 反例 4.32（精确 shared-Schur objective 一般不 submodular）

取 scalar map/pose、无 prior，两个频率块

$$
(A_1,B_1)=(1,1),
\qquad
(A_2,B_2)=(1,-1).
$$

单独任一块都可由 pose 完全补偿：

$$
K(\{1\})=K(\{2\})=0.
$$

但合并后

$$
A=\begin{bmatrix}1\\1\end{bmatrix},
\qquad
B=\begin{bmatrix}1\\-1\end{bmatrix},
\qquad A^TB=0,
$$

故

$$
K(\{1,2\})=2.
$$

于是添加第 2 块到空集的 gain 为 0，而添加到 $\{1\}$ 的 gain 为 2，违反 diminishing returns。加入很小的 $J_X=\alpha>0$ 与 map ridge 后该 violation连续保留；当 $\alpha\downarrow0$，任何 weak-submodularity ratio 可趋近 0。因此不能把 exact shared-pose frequency selection一般性地包装成 greedy-submodular theorem。


---

## 4.5 P5：半正定位姿先验、motion-graph gauge 与 weighted retention

### 定理 4.33（finite-dimensional generalized Schur range condition 自动成立）

令

$$
\mathcal I=
\begin{bmatrix}
A^TA&A^TB\\
B^TA&B^TB+J_X
\end{bmatrix},
\qquad J_X\succeq0,
$$

并记 $C=B^TB+J_X$。在有限维中，

$$
\boxed{
\operatorname{Ran}(B^TA)\subseteq\operatorname{Ran}C.
}
$$

因此 Moore--Penrose Schur complement

$$
K_{\rm eff}=A^TA-A^TB C^\dagger B^TA
$$

是合法的 generalized Schur complement，且 $K_{\rm eff}\succeq0$。

#### 证明

若 $h\in\ker C$，则

$$
0=h^TCh=\|Bh\|^2+h^TJ_Xh.
$$

两项非负，故 $Bh=0$ 且 $J_X^{1/2}h=0$。于是

$$
(A^TB)h=A^T(Bh)=0,
$$

即 $\ker C\subseteq\ker(A^TB)$。有限维中取正交补得到

$$
\operatorname{Ran}(B^TA)
=(\ker A^TB)^\perp
\subseteq(\ker C)^\perp
=\operatorname{Ran}C.
$$

PSD block matrix在该 range condition 下的 generalized Schur complement为 PSD。也可直接由下述 variational theorem 得到。证毕。

### 定理 4.34（无限维的短接/增广投影定义）

设 $A,B,J_X^{1/2}$ bounded，定义增广数据空间

$$
\overline{\mathcal Y}:=\mathcal Y\oplus\mathcal H_X
$$

和算子

$$
\bar A u=(Au,0),
\qquad
\bar B h=(Bh,J_X^{1/2}h).
$$

则

$$
\boxed{
q_{\rm eff}(u)
:=\inf_h
\left(\|Au-Bh\|^2+
\|J_X^{1/2}h\|^2\right)
=
\|P_{\overline{\operatorname{Ran}\bar B}^{\perp}}
\bar Au\|^2.
}
$$

因此存在唯一 bounded PSD operator

$$
\boxed{
K_{\rm eff}
=\bar A^*P_{\overline{\operatorname{Ran}\bar B}^{\perp}}
\bar A,
\qquad 0\preceq K_{\rm eff}\preceq K_{\rm IS}.
}
$$

若 $\operatorname{Ran}\bar B$ 闭，则

$$
\bar B^\dagger=(\bar B^*\bar B)^\dagger\bar B^*,
\qquad
\bar B^*\bar B=B^*B+J_X,
$$

且上述定义化为

$$
K_{\rm eff}
=A^*A-A^*B(B^*B+J_X)^\dagger B^*A.
$$

若 $\operatorname{Ran}\bar B$ 不闭，则最后一式可能涉及无界伪逆，必须保留 projection/infimum definition。

#### 证明

最小化式正是 $\bar Au$ 到线性子空间 $\operatorname{Ran}\bar B$ 的距离平方，而到任意线性子空间的距离等于到其闭包的距离。closed subspace 的最近点由 orthogonal projector 给出，故第一式成立。$P^\perp$ 是 positive contraction，给出 boundedness 与 Loewner bounds。闭值域时，orthogonal projector为 $\bar B\bar B^\dagger$，展开

$$
\bar A^*\bar B(\bar B^*\bar B)^\dagger\bar B^*\bar A
=A^*B(B^*B+J_X)^\dagger B^*A.
$$

证毕。

### Remark 4.35（无限维 range condition 的正确强度）

PSD 只自动给出

$$
\operatorname{Ran}(B^*A)
\subseteq\overline{\operatorname{Ran}(B^*B+J_X)^{1/2}},
$$

或等价的 kernel orthogonality。若 $C=B^*B+J_X$ 值域不闭，不能把 closure inclusion升级成 algebraic $\operatorname{Ran}(B^*A)\subseteq\operatorname{Ran}C$。这正是 shorted-operator formulation 比形式伪逆更稳健的原因。

### 定理 4.36（$\ker J_X$ / penalized subspace 分解）

设有限维，或在 Hilbert 空间中假设 $J_X$ closed range。令

$$
\mathcal N=\ker J_X,
\qquad
\mathcal R=\mathcal N^\perp,
\qquad
J_R=J_X|_{\mathcal R}\succ0,
$$

并写

$$
B=[B_0\ B_1],
\qquad B_0=B|_{\mathcal N},\quad B_1=B|_{\mathcal R}.
$$

令

$$
P_0=P_{\overline{\operatorname{Ran}B_0}}.
$$

则

$$
\boxed{
\begin{aligned}
K_{\rm eff}
={}&A^*P_0^\perp A\\
&-A^*P_0^\perp B_1
\bigl(B_1^*P_0^\perp B_1+J_R\bigr)^{-1}
B_1^*P_0^\perp A.
\end{aligned}
}
$$

解释如下。

- $\mathcal N$ 中的 pose directions 是完全自由 nuisance，先被 hard projection 掉；
- $\mathcal R$ 中的 pose directions 在 quotient data space $P_0^\perp\mathcal Y$ 中产生 finite shrinkage；
- 若 $J_R=\alpha J_0$，$\alpha\to\infty$，则

  $$
  K_{\rm eff}\to A^*P_0^\perp A,
  $$

  而不是一般的 $K_{\rm IS}$；只有 $B_0=0$ 时才达到 known-pose limit；
- $\alpha\downarrow0$ 时，在相应闭值域条件下趋于对全部 $\operatorname{Ran}B$ 的 free-pose projection。

#### 证明

对 $h=h_0+h_1$，先固定 $h_1$ 并对 $h_0\in\mathcal N$ 最小化：

$$
\inf_{h_0}\|Au-B_1h_1-B_0h_0\|^2
=\|P_0^\perp(Au-B_1h_1)\|^2.
$$

于是剩余问题为

$$
\inf_{h_1\in\mathcal R}
\bigl[\|P_0^\perp Au-P_0^\perp B_1h_1\|^2
+h_1^*J_Rh_1\bigr].
$$

因 $J_R\succ0$，normal matrix可逆，完成平方即得公式。极限由 inverse monotonicity 得到。证毕。

### 推论 4.37（motion-graph global gauge）

未锚定 pose-graph Laplacian 的 global translation/rotation modes 位于 $\ker J_X$。若相应 $B_0$ 数据方向非零，则即使其余 relative-motion prior趋于无穷强，这些 global modes 仍被 hard projection；known-pose limit不能通过只增强 relative odometry 获得。必须加入 anchor、已知 world feature、known background或其他绝对 reference。

### 定理 4.38（weighted retention 的基本谱界与等号条件）

有限维令

$$
W_J:=I-B(B^TB+J_X)^\dagger B^T,
\qquad
Q_A^TQ_A=I,
\qquad
R_X:=Q_A^TW_JQ_A.
$$

则

$$
\boxed{
0\preceq W_J\preceq I,
\qquad
0\preceq R_X\preceq I.
}
$$

并且对 $y\in\mathcal Y$，

$$
y^TW_Jy
=\min_h\bigl(\|y-Bh\|^2+h^TJ_Xh\bigr).
$$

因此：

1. $W_Jy=y$ 当且仅当 $B^Ty=0$；
2. 在有限维中，$W_Jy=0$ 当且仅当存在 $h\in\ker J_X$ 使 $y=Bh$；无限维中应改为存在 $h_n$ 使 $Bh_n\to y$、$J_X^{1/2}h_n\to0$；
3. $R_Xv=v$ 当且仅当 $B^TQ_Av=0$；
4. $R_Xv=0$ 当且仅当 $Q_Av$ 可由 prior-null pose direction无代价补偿（或其无限维闭包版本）。

#### 证明

variational expression显然位于 $[0,\|y\|^2]$，给出 operator bounds。若 $B^Ty=0$，$h=0$ 是 convex objective 的 minimizer，值为 $\|y\|^2$；反之若最小值等于 $\|y\|^2$，$h=0$ 必为 minimizer，first-order condition给 $B^Ty=0$。最小值为零时两个非负项同时为零，得到有限维条件；闭包版由 minimizing sequence得到。compression 的条件随即成立。证毕。

### 定理 4.39（关于 pose prior 的 Loewner 单调性、严格性与等号）

若

$$
0\preceq J_0\preceq J_1,
$$

则

$$
\boxed{
K_{\rm eff}(J_0)\preceq K_{\rm eff}(J_1).
}
$$

对给定 map direction $u$，令

$$
\mathcal M_0(u)
=\arg\min_h
\bigl(\|Au-Bh\|^2+h^TJ_0h\bigr).
$$

则

$$
u^T[K_{\rm eff}(J_1)-K_{\rm eff}(J_0)]u=0
$$

当且仅当存在 $h_0\in\mathcal M_0(u)$ 同时满足

$$
(J_1-J_0)^{1/2}h_0=0
$$

并且 $h_0$ 也是 $J_1$ 问题的 minimizer。若 $B^TB+J_0\succ0$，旧 minimizer唯一，

$$
h_0=(B^TB+J_0)^{-1}B^TAu,
$$

且等号简化为

$$
\boxed{
(J_1-J_0)^{1/2}h_0=0.
}
$$

因此在一个 map subspace $S$ 上严格 Loewner increase 的充分必要条件是：每个非零 $u\in S$ 的旧最优 compensating pose 都不落在新增 prior 的 kernel 中。

#### 证明

对任意 $h$，$J_1$ objective等于 $J_0$ objective加非负项
$h^T(J_1-J_0)h$，取 infimum得到单调性。若旧 minimizer唯一，取 $h_0$ 代入可见新增项为零时两值相等；若两值相等，new minimizer必须同时达到 old minimum且新增项为零，唯一性迫使其为 $h_0$。非唯一情形用 minimizer sets表述。证毕。

### 定理 4.40（prior-whitened shrinkage 与增广 principal angles）

若 $J_X\succ0$，令

$$
D=BJ_X^{-1/2}.
$$

则

$$
\boxed{
W_J=(I+DD^T)^{-1}.
}
$$

若 $D=U\operatorname{diag}(\eta_\ell)V^T$，则在 $\operatorname{Ran}D$ 的 left singular directions 上，$W_J$ 的 eigenvalues为

$$
\frac1{1+\eta_\ell^2},
$$

在 $\ker D^T$ 上为 1。故 finite-prior retention 是 prior-to-data ratio 决定的 weighted shrinkage。

另一方面，以定理 4.34 的 $\bar A,\bar B$ 为对象，在 $\overline{\mathcal Y}$ 中无惩罚地消元 $\bar B$，可得 generalized retention恰为

$$
\boxed{
\rho_i=\sin^2\bar\theta_i,
}
$$

其中 $\bar\theta_i$ 是 $\operatorname{Ran}\bar A$ 与
$\overline{\operatorname{Ran}\bar B}$ 的 principal angles。

因此：

- finite prior **可以**解释为增广数据空间的 ordinary principal angles；
- 它一般**不能**解释为原数据空间中 $\operatorname{Ran}A$ 与
  $\operatorname{Ran}B$ 的 ordinary $\sin^2\theta_i$；
- 原数据空间 $W_J$ 是 projector 当且仅当其 spectrum只取 $0,1$，这在 genuine finite $J_X\succ0$ 与 $B\neq0$ 时不成立。

#### 证明

Woodbury identity给出

$$
I-B(B^TB+J_X)^{-1}B^T
=I-D(I+D^TD)^{-1}D^T
=(I+DD^T)^{-1}.
$$

SVD 后逐方向得到 shrinkage factors。增广空间部分直接应用定理 4.1 于 $(\bar A,\bar B)$，且
$\bar A^*\bar A=A^*A$。证毕。

### 命题 4.41（block-banded prior、跨时刻相关误差与 map--pose 相关 prior）

1. 若 motion prior $J_X=L^TL$ block-banded/sparse，则使用增广 residual

   $$
   \begin{bmatrix}Au-Bh\\-Lh\end{bmatrix}
   $$

   可保留 sparsity；无需显式形成 dense $J_X$ 或 projector。
2. 若 measurement noise跨时刻相关，必须先用 global whitening factor $W_n$ 变换整个 stacked $[A,B]$；逐帧 Euclidean projector一般错误。
3. 若 joint prior precision为

   $$
   J_{\rm prior}=
   \begin{bmatrix}
   J_{\chi\chi}&J_{\chi X}\\
   J_{X\chi}&J_{XX}
   \end{bmatrix},
   $$

   则 map Schur operator应为

   $$
   \boxed{
   \begin{aligned}
   K_{\rm map}
   ={}&A^TA+J_{\chi\chi}\\
   &-(A^TB+J_{\chi X})
   (B^TB+J_{XX})^\dagger
   (B^TA+J_{X\chi}),
   \end{aligned}
   }
   $$

   并满足相应 range condition/shorted definition。普通 $A^TA-A^TB(\cdots)B^TA$ 在 map--pose correlated prior 下不正确。

---

## 4.6 P6：固定秩流形、广义谱导数、重根与 rank event

### 定理 4.42（constant-rank projector 的光滑性与一阶界）

令 $B(t)\in\mathbb R^{m\times p}$ 为 $C^k$ 曲线，并在开区间内保持
$\operatorname{rank}B(t)=r$。则 Moore--Penrose inverse $B(t)^\dagger$ 与

$$
P_B(t)=B(t)B(t)^\dagger
$$

均为 $C^k$。其一阶导数为

$$
\boxed{
\dot P_B
=(I-P_B)\dot BB^\dagger
+(B^\dagger)^T\dot B^T(I-P_B).
}
$$

令

$$
S=(I-P_B)\dot BB^\dagger.
$$

相对于 $\operatorname{Ran}B\oplus\operatorname{Ran}B^\perp$，
$\dot P_B$ 是 off-diagonal block matrix，因此

$$
\boxed{
\|\dot P_B\|_2=\|S\|_2
\le\frac{\|\dot B\|_2}{\sigma_r(B)}.
}
$$

Frobenius norm满足

$$
\|\dot P_B\|_F\le
\sqrt2\,\frac{\|\dot B\|_F}{\sigma_r(B)}.
$$

#### 证明

constant-rank matrices形成光滑流形，Moore--Penrose inverse在每个 rank stratum 上为 rational/smooth map。由
$P^2=P=P^T$ 求导可知 $P\dot PP=0$ 和
$P^\perp\dot PP^\perp=0$，故 $\dot P$ 只有 off-diagonal blocks。对
$P=BB^\dagger$ 使用标准 pseudoinverse derivative 或直接微分 QR/SVD，得到显示公式。off-diagonal self-adjoint matrix
$\begin{psmallmatrix}0&S^T\\S&0\end{psmallmatrix}$ 的 operator norm等于 $\|S\|$；再用
$\|B^\dagger\|=1/\sigma_r(B)$。证毕。

### 命题 4.43（coarse second-derivative bound）

若沿一条 constant-rank 曲线

$$
\sigma_r(B(t))\ge\beta>0,
\qquad
\|\dot B\|\le b_1,
\qquad
\|\ddot B\|\le b_2,
$$

则一个显式但非 sharp 的 bound 为

$$
\boxed{
\|\ddot P_B\|
\le \frac{2b_2}{\beta}
+\frac{8b_1^2}{\beta^2}.
}
$$

#### 证明

constant-rank pseudoinverse derivative identity给

$$
\|\dot B^\dagger\|
\le3\|B^\dagger\|^2\|\dot B\|
\le3b_1/\beta^2.
$$

对
$S=(I-P)\dot BB^\dagger$ 求导，并使用
$\|\dot P\|\le b_1/\beta$，得到

$$
\|\dot S\|
\le b_2/\beta+4b_1^2/\beta^2.
$$

由 $\ddot P=\dot S+\dot S^T$ 得结论。证毕。

### 定理 4.44（$K_{\rm IS}$、无先验与 finite-prior $K$ 的完整一阶导数）

令所有量沿参数 $t$ 可微。

1. known-pose information：

   $$
   \boxed{
   \dot K_{\rm IS}=\dot A^TA+A^T\dot A.
   }
   $$

2. no-prior、constant-rank $B$：

   $$
   \boxed{
   \begin{aligned}
   \dot K_{\rm SLAM}
   ={}&\dot A^T(I-P_B)A
   +A^T(I-P_B)\dot A\\
   &-A^T\dot P_BA.
   \end{aligned}
   }
   $$

3. 若 $C=B^TB+J_X\succ0$，记

   $$
   H=C^{-1},\qquad W=I-BHB^T,
   $$

   则

   $$
   \boxed{
   \dot K_{\rm eff}
   =\dot A^TWA+A^TW\dot A+A^T\dot WA,
   }
   $$

   其中

   $$
   \boxed{
   \dot C=\dot B^TB+B^T\dot B+\dot J_X,
   }
   $$

   $$
   \boxed{
   \dot W
   =-\dot BHB^T-BH\dot B^T
   +BH\dot C HB^T.
   }
   $$

若 $J_X$ 半正定且 $C$ 奇异，应改用增广 $\bar B$ 的 constant-rank projector derivative；只有在 $C^\dagger$ 的 rank stratum固定时才可逐项微分伪逆。

#### 证明

均由 product rule与
$\dot H=-H\dot C H$ 直接得到。证毕。

### 定理 4.45（simple generalized retention derivative）

设在固定 observable support 上 $G(t)=K_{\rm IS}(t)\succ0$，

$$
K(t)v(t)=\rho(t)G(t)v(t),
\qquad
v(t)^TG(t)v(t)=1,
$$

且 $\rho$ simple。则

$$
\boxed{
\dot\rho
=v^T(\dot K-\rho\dot G)v.
}
$$

#### 证明

微分 generalized eigen-equation：

$$
\dot Kv+K\dot v
=\dot\rho Gv+\rho\dot Gv+\rho G\dot v.
$$

左乘 $v^T$，利用 $v^TK=\rho v^TG$ 与 normalization，包含 $\dot v$ 的两项抵消，得到公式。证毕。

### 定理 4.46（重根簇的一阶分裂）

设 $\rho_0$ 是 multiplicity $s$ 的 isolated generalized eigenvalue，取

$$
V_c^TGV_c=I,
\qquad
KV_c=GV_c\rho_0.
$$

沿给定方向，一阶 eigenvalue splitting 为 Hermitian compression

$$
\boxed{
H_c:=V_c^T(\dot K-\rho_0\dot G)V_c
}
$$

的 $s$ 个 eigenvalues。单个 basis vector 的 Rayleigh derivative 在簇内没有 invariant 意义；只有 $H_c$ 的谱与 cluster projector 有意义。

#### 证明

将 perturbed eigenvectors写成 $V_cz+O(t)$，把 generalized eigen-equation投影到 $V_c$，零阶项抵消，$O(t)$ 项给
$H_cz=\dot\rho z$。证毕。

### 定理 4.47（Kato/Davis--Kahan projector bound）

设 support固定且 $G\succ0$。标准化

$$
T=G^{-1/2}KG^{-1/2}
$$

为 self-adjoint operator。设目标 cluster与其余谱的 separation至少为 $\delta>0$。若 perturbation满足

$$
\|\widetilde T-T\|<\delta/2,
$$

则对应 spectral projectors $P,\widetilde P$ 满足

$$
\boxed{
\|\widetilde P-P\|
\le\frac{2\|\widetilde T-T\|}{\delta}.
}
$$

沿光滑路径，Riesz projector

$$
P(t)=\frac{1}{2\pi i}
\oint_\Gamma(zI-T(t))^{-1}\,dz
$$

在 external gap保持时与 $T$ 同阶光滑。若 $G$ 的 support/rank改变，则该标准化本身不再固定，必须先处理 support event。

### 定理 4.48（rank change 的精确连续性分类）

对 finite-dimensional orthogonal projectors $P,Q$，若

$$
\operatorname{rank}P\neq\operatorname{rank}Q,
$$

则

$$
\boxed{
\|P-Q\|_2=1.
}
$$

因此若 $B(t)$ 连续但在 $t_0$ rank改变，并且 $P_{B(t)}$ 在 $t\neq t_0$ 与 $P_{B(t_0)}$ 秩不同，则

$$
\|P_{B(t)}-P_{B(t_0)}\|=1
$$

对任意足够近但不同的 $t$ 成立。$P_B$ 不仅非 differentiable，而且不 continuous、非 Hölder。

对

$$
B(t)=U\operatorname{diag}
(s_1,\ldots,s_{r-1},|t|s_r,0,\ldots)V^T,
$$

$t\neq0$ 时

$$
P_{B(t)}=
U\operatorname{diag}(1,\ldots,1,1,0,\ldots)U^T,
$$

$t=0$ 时缺少第 $r$ 个方向，故 jump为

$$
P_{B(t)}-P_{B(0)}=u_ru_r^T,
\qquad\|\cdot\|=1.
$$

然而

$$
K_{\rm SLAM}(t)-K_{\rm SLAM}(0)
=-A(0)^Tu_ru_r^TA(0)+o(1)
$$

可能因 $A(0)^Tu_r=0$ 而抵消。因此“$B$ rank event”是 $K$ 非连续的充分警报，不是无条件充分条件。若 $J_X\succeq\alpha I$，
$B(B^TB+J_X)^{-1}B^T$ 对 $B$ 光滑，rank change不再产生 hard jump。

#### 证明

若 $\operatorname{rank}P>\operatorname{rank}Q$，存在单位向量
$x\in\operatorname{Ran}P\cap(\operatorname{Ran}Q)^\perp$，故
$(P-Q)x=x$，得到 norm至少 1；两个 contractions之差 norm至多 1。其余直接展开。证毕。

### 算法 4.49（有条件可证明的 online rank-event detector）

设观测/数值误差有 operator-norm bound

$$
\|\widehat B-B\|\le\eta_B,
$$

一步内真实模型漂移 bound为 $d_B$。选 hysteresis thresholds

$$
0<\tau_{\rm off}<\tau_{\rm on},
\qquad
\tau_{\rm on}-\tau_{\rm off}>2(\eta_B+d_B).
$$

规则：

- inactive singular direction仅在 $\widehat\sigma_i>\tau_{\rm on}$ 时激活；
- active direction仅在 $\widehat\sigma_i<\tau_{\rm off}$ 时移除；
- 中间带保持上一状态。

Weyl inequality给出：

- 若真实 $\sigma_i\ge\tau_{\rm on}+\eta_B$，不会漏检为 inactive；
- 若真实 $\sigma_i\le\tau_{\rm off}-\eta_B$，不会误留为 active；
- 位于 uncertainty band 的事件本质上不可由当前误差模型 certifiably 分类。

对 generalized spectral clusters，若 standardized operator误差
$\|\widehat T-T\|\le\eta_T$，选

$$
g_{\rm split}-g_{\rm merge}>2(\eta_T+d_T),
$$

在 adjacent gap低于 $g_{\rm merge}$ 时合并 cluster，只在高于
$g_{\rm split}$ 时拆分。trust-region step应满足

$$
L_B\|\Delta X\|
+\frac12M_B\|\Delta X\|^2
<\operatorname{dist}
(\sigma(B),\{\tau_{\rm off},\tau_{\rm on}\})-(\eta_B+d_B).
$$

复杂度：dense full SVD约 $O(mp^2)$（$m\ge p$）；incremental rank-$r$ update约
$O(mpr+r^3)$。证明保证完全依赖于 $\eta_B,d_B,L_B,M_B$ 是真实上界；若它们来自 sampling而非 uniform analysis，则 detector只是 heuristic。


---

## 4.7 P7：full-wave resolvent 扰动界与真正的鲁棒性证书

### 先行修正 4.50（P7 原目标的导数阶数缺口）

因为

$$
B=D_XF,
$$

所以

$$
DB=D_X^2F,
\qquad
D^2B=D_X^3F.
$$

因此，若要给 $D^2B$ 的 uniform bound，必须假设 $G_S,G_D,e^{\rm inc}$ 对 pose 至少三阶可微，并控制二维 Hankel Green 函数三阶空间导数，或利用 rigid pullback结构避免直接微分 singular kernel。只控制 Green 一、二阶导数最多闭合到 $DB$，不能闭合 $D^2B$。以下给出满足三阶可微时的最强修正版；若只愿假设二阶，则删除所有下标 3 的项和 $D^2B$ 结论。

### 假设 4.51（compact parameter ball）

在紧参数集合 $\Theta$ 上，令

$$
M(X)=I-D_\chi G_D(X),
\qquad T(X)=M(X)^{-1},
$$

并假设

$$
\|T(X)\|\le t_0:=m_0^{-1},
\qquad
\|D_\chi\|\le c_\chi.
$$

定义 operator derivative suprema

$$
s_r:=\sup_{X\in\Theta}\|D^rG_S(X)\|,
\quad
d_r:=\sup_{X\in\Theta}\|D^rG_D(X)\|,
\quad
e_r:=\sup_{X\in\Theta}\|D^re^{\rm inc}(X)\|,
$$

$r=0,1,2,3$；范数为对应 symmetric multilinear operator norm。令 map injection $S_\chi$ bounded，$\|S_\chi\|=s_\chi$。

### 定理 4.52（resolvent inverse 的显式导数界）

沿任意单位 pose direction，记 $T_r$ 为 $\|D^rT\|$ 的 uniform upper bound。则可取

$$
\boxed{
T_0=t_0,
}
$$

$$
\boxed{
T_1\le t_0^2c_\chi d_1,
}
$$

$$
\boxed{
T_2\le
2t_0^3c_\chi^2d_1^2
+t_0^2c_\chi d_2,
}
$$

$$
\boxed{
T_3\le
6t_0^4c_\chi^3d_1^3
+6t_0^3c_\chi^2d_1d_2
+t_0^2c_\chi d_3.
}
$$

#### 证明

$M'= -D_\chi G_D'$，且

$$
T'=-TM'T=TD_\chi G_D'T.
$$

得到 $T_1$。再次求导：

$$
T''=2TD_\chi G_D'TD_\chi G_D'T
+TD_\chi G_D''T,
$$

得到 $T_2$。三阶求导按 product rule：三个一阶 $G_D'$ 的排列产生 $3!=6$ 项；一个 $G_D''$ 与一个 $G_D'$ 的排列及 inverse factors合计给系数 6；单个 $G_D'''$ 给最后一项。取范数即得。证毕。

### 定理 4.53（state、total field、map Jacobian 与 pose Jacobian 的递归界）

令

$$
j=TD_\chi e^{\rm inc},
\qquad
E=e^{\rm inc}+G_Dj,
\qquad
A=G_STD_ES_\chi,
\qquad
F=G_Sj.
$$

定义

$$
\boxed{
j_r
:=c_\chi\sum_{i=0}^r{r\choose i}T_i e_{r-i},
}
$$

$$
\boxed{
E_r
:=e_r+
\sum_{i=0}^r{r\choose i}d_i j_{r-i}.
}
$$

则 $\|D^rj\|\le j_r$、$\|D^rE\|\le E_r$。对 $r=0,1,2$，

$$
\boxed{
a_r
:=s_\chi
\sum_{i+j+\ell=r}
\frac{r!}{i!j!\ell!}
 s_iT_jE_\ell
}
$$

满足

$$
\|D^rA\|\le a_r.
$$

再定义

$$
\boxed{
f_r:=\sum_{i=0}^r{r\choose i}s_i j_{r-i}.
}
$$

则

$$
\boxed{
\|B\|=\|D_XF\|\le f_1,
\qquad
\|DB\|\le f_2,
\qquad
\|D^2B\|\le f_3.
}
$$

#### 证明

$j,E,A,F$ 都是 bounded multilinear product/composition。对每一阶使用 Leibniz multinomial formula，再代入定理 4.52 的 $T_r$ bounds。证毕。

### Remark 4.54（一般 $B$ 的附加项）

若 moving coordinate grid、direct path、antenna phase center、channel calibration 或 $D_\chi(X)$ 依赖 pose，应把相应 derivative suprema加入 product recursion。上述式子覆盖 world-fixed map 且 $D_\chi$ 固定的基本模型；不能把它声称为所有硬件 parameterization 的完整 $B$。

### 定理 4.55（二维外部 Hankel Green 的一、二阶统一界）

令

$$
g_k(r)=\frac{i}{4}H_0^{(1)}(k\|r\|),
$$

并假设

$$
k\in[k_{\min},k_{\max}],
\qquad
\|r\|\in[d_0,d_1],
\qquad d_0>0.
$$

定义

$$
C_\nu:=
\sup_{z\in[k_{\min}d_0,k_{\max}d_1]}
|H_\nu^{(1)}(z)|,
\qquad \nu=0,1,2.
$$

则

$$
\boxed{
|g_k(r)|\le\frac{C_0}{4},
}
$$

$$
\boxed{
\|\nabla g_k(r)\|
\le\frac{k_{\max}C_1}{4},
}
$$

$$
\boxed{
\|\nabla^2g_k(r)\|
\le
\frac{k_{\max}^2(C_0+C_2)}{8}
+\frac{k_{\max}C_1}{4d_0}.
}
$$

若 $G_S:L^2(D)\to\mathbb C^{m}$ 由 $m$ 个外部 receivers采样，则粗界

$$
\|D^rG_S\|
\le \sqrt{m|D|}\,
\sup_{\rm admissible}|D^rg_k|
$$

成立；连续 receiver curve/surface则将 $m$ 替换为其 measure。

#### 证明

radial derivative identity
$\partial_RH_0^{(1)}(kR)=-kH_1^{(1)}(kR)$ 给 gradient bound。radial Hessian为

$$
\phi''(R)\hat r\hat r^T
+\frac{\phi'(R)}R(I-\hat r\hat r^T),
$$

且

$$
H_1'(z)=\frac12(H_0(z)-H_2(z)).
$$

取 operator norm与 sup 即得。integral operator bound由 Cauchy--Schwarz。证毕。

### Remark 4.56（self-cell、内部奇性与为什么不能复制外部 standoff bound）

二维 Green kernel 在 $R\downarrow0$ 时

$$
H_0^{(1)}(kR)=O(1+|\log R|),
$$

故 $g_k$ 局部属于 $L^1$ 和 $L^2$，$G_D$ 可为 compact/Hilbert--Schmidt 类型。可是

$$
|\nabla g_k|=O(R^{-1})
$$

局部 $L^1$ 但不局部 $L^2$；

$$
|\nabla^2g_k|=O(R^{-2})
$$

不绝对可积，必须按 singular integral、distributional self-term 或 cell-averaged quadrature 处理。因此：

- 外部 sensor derivatives可由 $d_0>0$ 点态统一控制；
- 内部 $G_D$ 的高阶 coordinate derivatives不能直接套同一 sup bound；
- discrete self-cell rule可给有限常数，但常数可能随 mesh size $h$ 爆炸；
- world-fixed $G_D$ 若不随平台 pose变化，则 P7 的 trajectory derivative不需要对其 singular kernel做空间微分，这是最稳妥的基础模型。

若要求 $D^3F$，还需 $H_3$ 与 $d_0^{-2}$ 型外部界，或采用群作用/pullback分析。

### 定理 4.57（严格正 pose prior 下 $DK_{\rm eff}$ 与 $D^2K_{\rm eff}$ 的显式结构常数）

假设

$$
J_X\succeq\alpha I,
\qquad \alpha>0.
$$

令

$$
C=B^TB+J_X,
\qquad H=C^{-1},
\qquad W=I-BHB^T.
$$

给定 uniform bounds

$$
\|A\|\le a_0,
\quad\|DA\|\le a_1,
\quad\|D^2A\|\le a_2,
$$

$$
\|B\|\le b_0,
\quad\|DB\|\le b_1,
\quad\|D^2B\|\le b_2,
$$

以及 $\|DJ_X\|\le j_1$、$\|D^2J_X\|\le j_2$。定义

$$
c_1:=2b_0b_1+j_1,
\qquad
c_2:=2b_1^2+2b_0b_2+j_2,
$$

$$
h_0:=\alpha^{-1},
\qquad
h_1:=h_0^2c_1,
\qquad
h_2:=2h_0^3c_1^2+h_0^2c_2,
$$

$$
w_1:=2b_0b_1h_0+b_0^2h_1,
$$

$$
w_2:=2b_0b_2h_0+2b_1^2h_0
+4b_0b_1h_1+b_0^2h_2.
$$

则

$$
\boxed{
\|DK_{\rm eff}\|
\le 2a_0a_1+a_0^2w_1,
}
$$

$$
\boxed{
\|D^2K_{\rm eff}\|
\le 2a_0a_2+2a_1^2
+4a_0a_1w_1+a_0^2w_2.
}
$$

#### 证明

$\|H\|\le\alpha^{-1}=h_0$。由

$$
DH=-H(DC)H,
$$

$$
D^2H=2H(DC)H(DC)H-H(D^2C)H
$$

得到 $h_1,h_2$。对 $W=I-BHB^T$ 两次 product differentiation得到 $w_1,w_2$。最后对 $K=A^TWA$ 两次求导；因 $0\preceq W\preceq I$，$\|W\|\le1$，逐项取范数即得。证毕。

### 命题 4.58（无先验 fixed-rank 常数）

若 $J_X=0$、$B$ 在参数球内 fixed rank，且

$$
\sigma_{\min}^+(B)\ge\beta_0>0,
$$

则

$$
\|DP_B\|\le b_1/\beta_0,
$$

$$
\|D^2P_B\|
\le2b_2/\beta_0+8b_1^2/\beta_0^2.
$$

代入

$$
K_{\rm SLAM}=A^T(I-P_B)A
$$

即可得到与定理 4.57 同型的 $DK,D^2K$ bounds。若 $\beta_0\downarrow0$，这些常数必然发散；反例 4.48 表明这不是 proof artifact，而是 projector 真正失去连续性。

### Corollary 4.59（参数爆炸率与近共振 tradeoff）

从定理 4.52--4.58 可见：

- $\|A\|$ 至少含一个 $t_0=m_0^{-1}$ factor；
- $\|DA\|$ 含 $t_0^2$ 型项，$\|D^2A\|$ 含 $t_0^3$ 型项；
- $K_{\rm IS}=A^TA$ 的 absolute scale可按 $t_0^2$ 增大；
- $DK$ 常含 $t_0^3$ 或更高次，relative sensitivity通常多出一个 $t_0$；
- 外部 geometry derivatives随 $k$、$k^2$、$d_0^{-1}$ 增大；
- no-prior nuisance geometry额外随 $\beta_0^{-1},\beta_0^{-2}$ 爆炸；finite prior则随 $\alpha^{-1},\alpha^{-2},\alpha^{-3}$ 爆炸。

因此“近共振信息变大”与“模型/轨迹敏感度更快变大”是竞争关系。只比较 $\lambda(K)$ 而不同时报告 derivative constant 与 admissible trust radius，会给出误导性的 trajectory/design conclusion。

### 不可能性结论 4.60（何时没有 continuum uniform constant）

若缺少下列任一项：

- $m_0>0$ 的 uniform resolvent margin；
- $d_0>0$ 的 external standoff；
- fixed-rank $\beta_0>0$ 或 $J_X\succeq\alpha I$；
- 对 map/geometry 的足够 Sobolev regularity；
- internal singular integral/self-cell 的一致处理；

则一般不能在 continuum parameter class上得到有限 uniform $DK$ 或 $D^2K$ 常数。此时可证明的版本必须限制为：

1. 固定离散网格与紧参数球；或
2. 明确的 Sobolev mapping theorem + uniform inverse bound；或
3. strictly positive pose prior + external-only geometry variation。

随机采样方向没有违反只说明 sampled quotient，没有把上述 suprema变成 certified bounds。

---

## 4.8 P8：robust trajectory 的一阶/二阶定理与约束锥

### 定理 4.61（simple eigenvalue 的 robust first-order expansion）

设 $f(X)=\lambda_r(K(X))$，目标 eigenvalue在闭球
$\|H\|\le\varepsilon_0$ 内 simple，并与其余谱有 uniform gap $\delta>0$。设 $K$ 为 $C^2$，且

$$
\|D^2f(X+H)\|\le L_f
\quad\forall\|H\|\le\varepsilon_0.
$$

则对任意参数 norm $\|\cdot\|$ 及其 dual norm $\|\cdot\|_*$，

$$
\boxed{
\min_{\|H\|\le\varepsilon}f(X+H)
=f(X)-\varepsilon\|Df(X)\|_*+R_2,
}
$$

其中

$$
\boxed{
|R_2|\le\frac12L_f\varepsilon^2,
\qquad 0\le\varepsilon\le\varepsilon_0.
}
$$

若 $v$ 是 normalized eigenvector，

$$
Df(X)[H]=v^TDK(X)[H]v.
$$

一个可计算的 Hessian bound 为

$$
\boxed{
L_f
\le M_2+\frac{2M_1^2}{\delta},
}
$$

其中

$$
M_1=\sup\|DK\|,
\qquad
M_2=\sup\|D^2K\|.
$$

对 generalized retention先标准化
$T=K_{\rm IS}^{-1/2}K_{\rm eff}K_{\rm IS}^{-1/2}$ 于固定 positive support，再应用同一定理。

#### 证明

Taylor theorem给

$$
f(X+H)=f(X)+Df(X)[H]+r(H),
\qquad |r(H)|\le\frac12L_f\|H\|^2.
$$

线性泛函在 norm ball 上的最小值为
$-\varepsilon\|Df(X)\|_*$。分别用达到/逼近 dual norm 的方向作上界，并对任意 $H$ 作下界，得到 remainder sandwich。

simple eigenvalue 的二阶 perturbation公式为

$$
D^2f[H,H]
=v^TD^2K[H,H]v
+2\sum_{j\ne r}
\frac{|v_j^TDK[H]v|^2}{\lambda_r-\lambda_j}.
$$

取绝对值、使用 gap $\delta$ 与 Bessel inequality，得到
$M_2+2M_1^2/\delta$。证毕。

### 定理 4.62（motion/collision/equality constraints 的 tangent-cone worst direction）

设 feasible set $\mathcal X$ 在 $X$ 的 Bouligand tangent cone为 $C_X$。对 scaled local uncertainty

$$
H=\varepsilon d+o(\varepsilon),
\qquad d\in C_X,
\qquad\|d\|\le1,
$$

有

$$
\boxed{
\min f(X+H)
=f(X)+\varepsilon
\min_{d\in C_X,\ \|d\|\le1}
Df(X)[d]+o(\varepsilon).
}
$$

该 coefficient等于

$$
\boxed{
-\sigma_{C_X\cap\mathbb B}(-\nabla f),
}
$$

即 $C_X\cap$ unit ball 对 $-\nabla f$ 的 support function取负。若 norm为 Euclidean 且 $C_X$ 为 closed convex cone，则

$$
\boxed{
\min_{d\in C_X,\ \|d\|_2\le1}
\nabla f^Td
=-\|\Pi_{C_X}(-\nabla f)\|_2.
}
$$

对 smooth constraints

$$
c_j(X)=0,
\qquad g_i(X)\le0,
$$

在 regular point，linearized cone为

$$
C_X=
\{d:Dc_j(X)d=0,
\ Dg_i(X)d\le0\ \text{for active }i\}.
$$

#### 证明

由 tangent-cone定义和 directional Taylor expansion得到第一式。第二式是 support function定义。Euclidean convex cone情形用 Moreau decomposition
$z=\Pi_Cz+\Pi_{C^\circ}z$，最大内积由归一化后的 $\Pi_C(-\nabla f)$ 达到。证毕。

### 定理 4.63（重根簇的 worst first-order problem）

设 $\lambda_0$ 是 multiplicity $s$ 的 isolated cluster，$V^TV=I$。若 trajectory坐标为 $X_1,\ldots,X_q$，定义

$$
H_j=V^T\frac{\partial K}{\partial X_j}V.
$$

沿方向 $d\in\mathbb R^q$，cluster的一阶分裂为

$$
\operatorname{eig}
\left(\sum_{j=1}^qd_jH_j\right).
$$

若目标是 cluster 的 lowest branch，则 worst first-order coefficient为

$$
\boxed{
\min_{d\in\mathcal U_1}
\lambda_{\min}
\left(\sum_jd_jH_j\right).
}
$$

$\lambda_{\min}$ 是 affine matrix map 的 concave function，因此在一般 convex uncertainty set上最小化它是 nonconvex concave minimization。特殊情形：

1. 若所有 $H_j$ simultaneously diagonalizable，问题降为有限个 linear function的下包络；
2. 若 $\mathcal U_1$ 是 polytope，global minimum至少有一个 extreme point解，可枚举小规模 vertices；
3. 若要验证给定 $t$ 对所有 polytope vertices满足

   $$
   \sum_jd_jH_j\succeq tI,
   $$

   则是 finite SDP feasibility；
4. 对 ellipsoid/一般 cone，直接把 $d$ 与 eigenvector共同优化会产生 bilinear/nonconvex问题，不能统称为 convex SDP。

### 定理 4.64（rank-event 安全 trust region）

设当前 numerical rank threshold为 $\tau_B$，当前 $B(X)$ 的 singular values与 threshold最小距离为

$$
m_B=\min_i|\sigma_i(B(X))-\tau_B|.
$$

若

$$
\|DB\|\varepsilon+\frac12\|D^2B\|\varepsilon^2<m_B,
$$

则 Weyl inequality保证该 trust region内 threshold rank不变。类似地，若目标 spectral cluster外 gap为 $\delta_K$，且

$$
\|DK\|\varepsilon+\frac12\|D^2K\|\varepsilon^2<\delta_K/2,
$$

则 cluster不与外部谱合并。

在这些条件失败时：

- 不得继续使用 simple-eigenvalue gradient；
- hard no-prior projector value可能不 locally Lipschitz；Clarke gradient也未必适用；
- 可选修正是 $J_X\succeq\alpha I$、缩小 trust region、显式 cluster tracking，或改成 bounded compensating pose game

  $$
  \min_{\|h\|\le H}
  \|Au-Bh\|^2,
  $$

  其 feasible set compact，value对 $(A,B)$ 连续，但它不等同于无界 nuisance Schur projector。

### 设计原则 4.65（无 universal trajectory ordering 的 Pareto formulation）

每条 candidate trajectory 至少同时报告：

1. **absolute map information**：task-restricted $\lambda_r(K_{\rm eff})$、trace 或 regularized log-det；
2. **relative retention**：$R_{\rm geom}$ 或 $K_{\rm IS}$-whitened spectrum；
3. **standoff / safety margin**：最小距离与 collision constraints；
4. **energy / time / path length**；
5. **local robust lower bound**：

   $$
   \lambda_r-\varepsilon\|D\lambda_r\|_*
   -\frac12L_r\varepsilon^2;
   $$

6. **rank/gap margins**：$\sigma_{\min}^+(B)$、target cluster gap；
7. **full-wave state margin**：$\sigma_{\min}(M)$ 或 $\|M^{-1}\|$。

选择应在这些轴上做 Pareto dominance或显式 scalarization，而不是预设“圆优于直线”“带宽越宽越好”。


---

## 4.9 P9：SOM/current subspace 与 map-tangent/pose-defect subspace 的严格桥接

令 current Hilbert space为 $\mathcal H_J$，且

$$
G_S:\mathcal H_J\to\mathcal Y,
\qquad
T_\chi:\mathcal H_\chi\to\mathcal H_J,
\qquad
A=G_ST_\chi.
$$

在 full-wave contrast-source linearization中，典型

$$
T_\chi=M^{-1}D_E S_\chi.
$$

### 定理 4.66（range containment、等号与 closure）

总有

$$
\boxed{
\operatorname{Ran}A
=G_S(\operatorname{Ran}T_\chi)
\subseteq\operatorname{Ran}G_S.
}
$$

精确等号成立当且仅当

$$
\boxed{
\mathcal H_J
=\operatorname{Ran}T_\chi+
\ker G_S.
}
$$

更一般地，closure equality

$$
\overline{\operatorname{Ran}A}
=
\overline{\operatorname{Ran}G_S}
$$

当且仅当 $\operatorname{Ran}T_\chi$ 在 quotient seminorm

$$
\|[j]\|_{G_S}:=\|G_Sj\|,
\qquad [j]\in\mathcal H_J/\ker G_S
$$

中稠密。一个易检验的充分条件是

$$
\overline{\operatorname{Ran}T_\chi+
\ker G_S}^{\,\mathcal H_J}
=\mathcal H_J.
$$

若 $\operatorname{Ran}G_S$ 闭，则 $G_S$ 在
$(\ker G_S)^\perp$ 上 bounded below，quotient norm与原 Hilbert quotient norm等价，上述充分条件同时为必要条件。

#### 证明

containment由 composition直接得到。若 exact equality，任意 $j\in\mathcal H_J$ 的 $G_Sj$ 可写为
$G_ST_\chi u$，故 $j-T_\chi u\in\ker G_S$，即
$j\in\operatorname{Ran}T_\chi+\ker G_S$；反向显然。

closure equality等价于每个 $G_Sj$ 可由 $G_ST_\chi u_n$ 在数据 norm下逼近，即
$\|G_S(j-T_\chi u_n)\|\to0$，正是 quotient seminorm稠密。闭值域时 reduced minimum modulus为正，数据 norm控制到 quotient distance，故与原空间 density等价。证毕。

### 反例 4.67（严格包含）

取 $G_S=I$ 于 $\mathbb R^2$，$T_\chi=P_{e_1}$。则

$$
\operatorname{Ran}A=\operatorname{span}\{e_1\}
\subsetneq\mathbb R^2=
\operatorname{Ran}G_S.
$$

因此 current-to-data 可见的 mode不一定能由任何 admissible map tangent激发。

### 定理 4.68（composition singular-value inequalities）

若 $G_S$ 或 $T_\chi$ compact，则 approximation/singular numbers满足

$$
\boxed{
s_n(G_ST_\chi)
\le\|T_\chi\|s_n(G_S),
\qquad
s_n(G_ST_\chi)
\le\|G_S\|s_n(T_\chi).
}
$$

更一般地，

$$
\boxed{
s_{n+m-1}(G_ST_\chi)
\le s_n(G_S)s_m(T_\chi).
}
$$

这说明 $G_S$ 的 dominant current DoF 只是 $A$ map DoF 的上游瓶颈；$T_\chi$ 的 kernel、ill-conditioning 与 range alignment可进一步减少/旋转 map modes。

### 语义定理 4.69（right singular spaces 不能直接比较）

$G_S$ 的 right singular vectors位于 $\mathcal H_J$，而 $A=G_ST_\chi$ 的 right singular vectors位于 $\mathcal H_\chi$。除非给出一个明确的 isometric identification或将 map modes通过 $T_\chi$ push forward，否则

$$
\operatorname{gap}(V_{G_S},V_A)
$$

没有定义。可合法比较的对象包括：

1. data-left subspaces of $G_S$ and $A$；
2. current vectors $T_\chi V_A$ 与 $V_{G_S}$，但需重新正交化并处理 $\ker T_\chi$；
3. pullback information $T_\chi^*G_S^*G_ST_\chi$ 与 current normal $G_S^*G_S$ 的谱关系。

### 定理 4.70（data-left dominant subspace 的条件性 Davis--Kahan 界）

令

$$
H_0=cG_SG_S^*,
\qquad
H_1=AA^*=G_ST_\chi T_\chi^*G_S^*,
$$

$c>0$。则

$$
\|H_1-H_0\|
\le
\|G_S\|^2
\|T_\chi T_\chi^*-cI\|.
$$

设 $H_0$ 的目标 data-left spectral cluster与其余谱 gap为 $\delta>0$，且右侧小于 $\delta/2$。则对应 projectors满足

$$
\boxed{
\|P_{A,\rm left}-P_{G_S,\rm left}\|
\le
\frac{2\|G_S\|^2
\|T_\chi T_\chi^*-cI\|}{\delta}.
}
$$

该结论要求 $T_\chi T_\chi^*$ 在 relevant current space上近似 scalar identity。仅知道 $T_\chi$ 的 condition number不足。

#### 证明

第一式由 factorization直接取 norm，第二式应用 Davis--Kahan。证毕。

### 反例 4.71（$T_\chi$ condition number 为 1 仍可 mode angle 为 $\pi/2$）

取

$$
G_S=\operatorname{diag}(2,1):\mathbb R^2\to\mathbb R^2,
$$

其 dominant current/right 与 data-left mode均为 $e_1$。令一维 map space通过 isometry

$$
T_\chi u=ue_2.
$$

$T_\chi$ 在其 domain上的 condition number为 1，但

$$
A u=ue_2,
$$

故 map data mode为 $e_2$，与 $G_S$ dominant data mode正交。缺失的条件是 range alignment，而不是 conditioning。

### 定理 4.72（固定 current reduction 与 pose elimination）

令 $P_S$ 是 current space上的固定 orthogonal projector，定义 reduced map tangent

$$
A_S=G_SP_ST_\chi.
$$

对固定 pose weight $W$（无先验时 $W=I-P_B$，finite prior时为相应 contraction），先做 current reduction再 pose elimination得到

$$
K_S=T_\chi^*P_SG_S^*WG_SP_ST_\chi.
$$

若所谓“反序”只是先形成同一个 $W$ 再把同一个固定 $P_S$ 插入，则两者代数上当然相同；真正不交换的是：

- 先按 $G_S^*G_S$ 选择 dominant $P_S$；
- 与先按 pose-eliminated current normal

  $$
  H_W=G_S^*WG_S
  $$

  重新选择 dominant subspace。

固定 $P_S$ 对 $H_W$ 无 cross-coupling的充要条件是

$$
\boxed{
[P_S,H_W]=0.
}
$$

若还要求两个 adaptive dominant projectors相同，则除了 reducing condition外，还必须保证 $P_S$ 内外的 $H_W$ eigenvalue ordering不跨越所用 cutoff，并有正 spectral gap。

在数据空间，若 $Q_S=P_{\operatorname{Ran}(G_SP_S)}$，则 current-data reduction与 hard pose projection按任意顺序相同的充要条件是

$$
\boxed{
[Q_S,P_B]=0.
}
$$

#### 证明

$P_S$ reducing $H_W$ 当且仅当 $P_SH_WP_S^\perp=0$，对 self-adjoint $H_W$ 等价于 commutator为零。两个 orthogonal projectors可交换当且仅当空间分解为四个交叉 intersections，因而投影次序相同。adaptive spectral equality还需 eigenvalue排序不改变。证毕。

### 反例 4.73（SOM split 与 pose-eliminated split 不交换）

取

$$
G_S=\begin{bmatrix}2&0\\0&1\end{bmatrix},
\qquad
\operatorname{Ran}B=\operatorname{span}\{(1,1)^T\}.
$$

原 current normal为

$$
G_S^TG_S=\operatorname{diag}(4,1),
$$

其 top mode是 $e_1$。无先验 pose projector

$$
W=I-\frac12
\begin{bmatrix}1&1\\1&1\end{bmatrix}
=\frac12
\begin{bmatrix}1&-1\\-1&1\end{bmatrix}.
$$

pose-eliminated current normal为

$$
G_S^TWG_S
=\frac12
\begin{bmatrix}4&-2\\-2&1\end{bmatrix},
$$

其唯一非零 eigenvector与 $(2,-1)^T$ 平行，不是 $e_1$。因此先做 $G_S$ SVD truncation与先做 pose elimination再选 dominant current mode得到不同 subspaces。

### 定理 4.74（精确 intersection projector 与近似 CS decomposition）

设 $P,Q$ 为 closed subspaces $M,N$ 的 orthogonal projectors。

1. $PQ$ 是 $M\cap N$ 的 orthogonal projector当且仅当 $P,Q$ commute；一般 $PQ$ 不是 self-adjoint，也不是 projector。
2. 取 orthonormal bases $U_M,U_N$，令

   $$
   C=U_M^TU_N.
   $$

   $\sigma_i(C)=1$ 对应 exact intersection directions；$\sigma_i(C)\approx1$ 对应 small principal angles，而不是 exact intersection。
3. 若 numerical perturbation $\|\widehat C-C\|\le\varepsilon$，目标 singular cluster与其余 singular values gap为 $g>2\varepsilon$，则 Wedin-type bound给

   $$
   \boxed{
   \|\sin\Theta(\widehat U,U)\|
   \le\frac{\varepsilon}{g-\varepsilon}.
   }
   $$

4. alternating projections $(PQP)^k$ 在 closed-sum/positive Friedrichs angle条件下收敛到 $P_{M\cap N}$，其线性速率由 Friedrichs cosine控制；当 angle趋零，收敛任意慢。

这否定了把任意 projector product 或“八象限”直接叫作 exact intersection decomposition 的做法。

### 算法 4.75（数学合法、CPU 可实现的 SOM--pose two-stage reduction）

**目标算子。** 近似 full map information

$$
K_{\rm eff}=A^*WA,
\qquad A=G_ST_\chi,
$$

其中 $W$ 是 hard pose projection或 finite-prior contraction。

**步骤。**

1. 用 Lanczos/randomized SVD求 $G_S$ 的 top-$r$ current basis $V_r$；若实施 Chen TSOM，则用 $G_D$ major basis与 $V_S^-$ 做 CS decomposition，构造

   $$
   U_T=[V_S^+,B_{S^-D^+}],
   $$

   而不是使用未经证明的 projector product。
2. 令 $P_T=U_TU_T^*$，构造 reduced map tangent

   $$
   A_r=G_SP_TT_\chi.
   $$
3. 对 $A_r$ 和 $B$ 做 rank-revealing QR。无 prior时对

   $$
   Z=Q_B^TQ_{A_r}
   $$

   做 SVD；finite prior时改在增广空间

   $$
   \bar A_r=(A_r,0),
   \qquad
   \bar B=(B,J_X^{1/2})
   $$

   做 QR/CS decomposition。
4. 以 reduced operator

   $$
   K_{{\rm eff},r}=A_r^*WA_r
   $$

   求 absolute spectrum与 retention；若只保留 $k$ 个 pose-defect singular modes，则显式记录 tail。

**误差证书。** 令

$$
E=A-A_r=G_S(I-P_T)T_\chi.
$$

因 $\|W\|\le1$，

$$
\boxed{
\|K_{\rm eff}-K_{{\rm eff},r}\|
\le(\|A\|+\|A_r\|)\|E\|.
}
$$

若 $P_T$ 至少包含 $G_S$ 的 top-$r$ right singular space，则

$$
\boxed{
\|E\|
\le\sigma_{r+1}(G_S)\|T_\chi\|.
}
$$

无 prior时，若 defect cross-matrix $Z$ 截断到 rank $k$，则 normalized retention operator误差为

$$
\boxed{
\|Z^TZ-Z_k^TZ_k\|
=\sigma_{k+1}(Z)^2.
}
$$

**复杂度。** 设一次 $G_S/G_S^*$ matvec成本 $C_G$：

- randomized/Lanczos basis约 $O(N_{\rm mv}C_G+(m+n)r^2)$；
- QR/CS cross geometry约 $O(mr^2+mp^2+mrp)$；
- pose solve/eigendecomposition约 $O(p^3+r^3)$，可用低秩结构进一步降低。

**失败触发条件。**

- $\sigma_r(G_S)-\sigma_{r+1}(G_S)$ 太小，dominant current projector不稳定；
- $T_\chi$ 强烈放大 discarded current directions；
- $B$ rank event 或 augmented subspace gap闭合；
- map basis不包含 physical gauge generators；
- 把该算法称为“Chen TSOM”但实际上未实现 $G_D$ split。若只做 $G_S$ + pose 两层，应称 **SOM--pose two-stage reduction**。

---

## 4.10 P10：从局部谱到实际估计误差

### 统计语义 4.76

考虑 linearized whitened model

$$
y=Au_0+Bx_0+n,
\qquad n\sim N(0,I).
$$

三种实验不可混写。

1. **fixed nuisance repeated sampling：** $x_0$ 固定；每次只重采 data noise。
2. **Bayesian random nuisance：** $x\sim N(0,J_X^{-1})$ 也是生成模型的一部分，并被边缘化。
3. **penalized/misspecified M-estimator：** objective含 deterministic penalty或 likelihood与真实生成机制不匹配；covariance为 sandwich/Godambe形式。

### 定理 4.77（fixed true pose + deterministic pose penalty 的 sandwich covariance）

假设 $C=B^TB+J_X\succ0$，令

$$
W=I-BC^{-1}B^T,
\qquad
K_{\rm eff}=A^TWA.
$$

估计器定义为

$$
(\widehat u,\widehat x)
=\arg\min_{u,x}
\bigl(\|y-Au-Bx\|^2+x^TJ_Xx\bigr).
$$

消元 $x$ 后

$$
\widehat u=K_{\rm eff}^{-1}A^TWy
$$

（在 identifiable support上）。若真实 $x_0=0$ 固定，只重采 $n$，则

$$
\boxed{
\operatorname{Cov}_{n}(\widehat u)
=K_{\rm eff}^{-1}A^TW^2A K_{\rm eff}^{-1}.
}
$$

一般不等于 $K_{\rm eff}^{-1}$，因为 finite prior下 $W^2\neq W$。

若真实 fixed $x_0\neq0$，则 map bias为

$$
\boxed{
\operatorname{Bias}(\widehat u)
=K_{\rm eff}^{-1}A^TWBx_0.
}
$$

特殊情形：

- known pose：$W=I$，covariance为 $K_{\rm IS}^{-1}$；
- no-prior hard projection：$W=P_{B^\perp}=W^2$，在 support上 covariance为 $K_{\rm SLAM}^{-1}$；
- genuine finite penalty：通常必须用 sandwich。

#### 证明

对 $x$ 完成平方得到 profiled objective

$$
(y-Au)^TW(y-Au).
$$

normal equation给估计器公式。代入 $y=Au_0+Bx_0+n$，线性变换噪声的 covariance为

$$
K^{-1}A^TW I W A K^{-1}
=K^{-1}A^TW^2AK^{-1}.
$$

取期望得到 bias。证毕。

### 定理 4.78（Bayesian random-pose marginal 下 $K_{\rm eff}^{-1}$ 成立）

Assume $J_X\succ0$.
若真实生成模型为

$$
x\sim N(0,J_X^{-1}),
\qquad n\sim N(0,I),
\qquad x\perp n,
$$

则给定 $u$ 的 marginal data covariance为

$$
\Sigma_y=I+BJ_X^{-1}B^T.
$$

Woodbury identity给

$$
\Sigma_y^{-1}
=I-B(B^TB+J_X)^{-1}B^T=W.
$$

故 marginal Fisher information为

$$
A^T\Sigma_y^{-1}A=K_{\rm eff}.
$$

在 flat map prior、linear Gaussian、full identifiable support下，map posterior covariance / GLS covariance为

$$
\boxed{
K_{\rm eff}^{-1}.
}
$$

同一结论也适用于：$J_X$ 来自真实 auxiliary pose measurements，并且这些辅助 measurement noises在每次 frequentist trial中也一同重采、joint likelihood正确指定。若 pose prior只是 optimizer regularizer而不是生成模型或随机辅助数据，不能使用此结论。

### 定理 4.79（一般 penalized/misspecified M-estimator 的 sandwich）

设 estimating equation为

$$
\sum_{i=1}^N\psi_i(\theta)+\nabla R(\theta)=0,
$$

而 $R$ 是 deterministic penalty。令

$$
H=E\left[-\frac{\partial\psi_i}{\partial\theta}\right]
+\frac1N\nabla^2R,
\qquad
V=\operatorname{Var}(\psi_i).
$$

在标准 M-estimation regularity 下，

$$
\boxed{
\sqrt N(\widehat\theta-\theta_*)
\Rightarrow N(0,H^{-1}VH^{-T}).
}
$$

penalty curvature进入 $H$，但 deterministic penalty没有对应 sampling score variance，不进入 $V$；这正是 Hessian inverse与 covariance分离的机制。

### 定理 4.80（small-noise local asymptotic normality / local NLS validity）

设 quotient 后 finite/task-dimensional参数 $\theta=(u,x)$，观测

$$
y_\varepsilon=F(\theta_0)+\varepsilon n,
\qquad n\sim N(0,\Sigma),
$$

并假设：

1. $F$ 在 $\theta_0$ 邻域 $C^3$；
2. gauge已 quotient/anchored，efficient Jacobian在 horizontal space full rank；
3. truth为 interior point；
4. likelihood正确；
5. 存在随 $\varepsilon\to0$ 落入同一 local basin 的 estimator/initialization；
6. 若使用 prior，其 asymptotic scaling与所声称的生成模型一致。

则 local experiment为 LAN，local MLE/NLS满足

$$
\varepsilon^{-1}(\widehat\theta-\theta_0)
\Rightarrow
N(0,\mathcal I(\theta_0)^{-1}),
$$

map block由相应 Schur/efficient information给出。无 prior时，在
$K_{\rm IS}$-whitened coordinates中，relative covariance inflation eigenvalues为

$$
\boxed{
\rho_i^{-1}
}
$$

对所有 $\rho_i>0$；$\rho_i=0$ 表示一阶不可辨识而非有限方差。

该定理只控制落在正确 local basin的 estimator；不排除远处 aliases、局部极小或错误 data association。

### 定理 4.81（一个可计算的 local injectivity / Gauss--Newton basin lower bound）

设 quotient/horizontal derivative $J_0=DF(\theta_0)$ 满足

$$
\|J_0h\|\ge\mu\|h\|
$$

且 Jacobian在半径 $r$ 球内 Lipschitz：

$$
\|DF(\theta)-DF(\theta_0)\|
\le L\|\theta-\theta_0\|.
$$

则对 $\|h\|\le r$，

$$
\boxed{
\|F(\theta_0+h)-F(\theta_0)\|
\ge
\left(\mu-\frac L2\|h\|\right)\|h\|.
}
$$

因此当 $r<2\mu/L$ 时，名义点相对于自身具有 local metric separation。它是 local basin分析的必要材料，但不保证整个球内任意两点间的 global injectivity，也不保证给定算法从球外进入该 basin。

#### 证明

积分 Taylor formula给

$$
F(\theta_0+h)-F(\theta_0)
=J_0h+
\int_0^1[DF(\theta_0+th)-J_0]h\,dt.
$$

remainder norm至多

$$
\int_0^1Lt\|h\|^2dt
=\frac L2\|h\|^2.
$$

与 $\|J_0h\|\ge\mu\|h\|$ 合并。证毕。

### 反例 4.82（Fisher 很好但 global estimation 失败）

1. **phase alias / cycle skipping：**

   $$
   F(x)=e^{ikx},
   \qquad x\in\mathbb R.
   $$

   local Fisher information为 $k^2/\sigma^2$，可随 $k$ 任意大；但

   $$
   F(x+2\pi n/k)=F(x)
   $$

   对所有整数 $n$ 成立，likelihood具有周期性多峰。高频同时提高 local curvature与 alias密度。
2. **sign ambiguity：**

   $$
   F(x)=x^2.
   $$

   在 $x_0\neq0$ local Fisher为 $4x_0^2/\sigma^2>0$，但 $x_0$ 与 $-x_0$ 全局不可区分。
3. **weakly identified nuisance：** 即使 efficient Fisher正定，若 $\sigma_{\min}B$、spectral gap或 resolvent margin很小，linearization radius可比噪声诱导误差更小，local Gaussian approximation在实际噪声下失效。

### 论文最强安全表述 4.83

可直接使用的边界为：

> 本文的 retention spectrum quantifies the loss of **local tangent Fisher/Gauss--Newton information** after eliminating physical pose increments in a whitened and realified model. Under a correctly specified local Gaussian or small-noise LAN regime, fixed quotient gauge, positive efficient-information gap, and initialization inside the corresponding local basin, the reciprocal retention spectrum describes relative covariance inflation in $K_{\rm IS}$-whitened coordinates. It does not by itself imply global uniqueness, absence of phase aliases or cycle skipping, convergence from arbitrary initialization, or real-system reconstruction performance.


---

# 5. 对当前论文的可直接替换文本

本节不是新的证明层，而是把第 4 节已经证明的结果压缩成可直接放入主文与附录的 LaTeX。为避免把 continuum 与 discretization 混在同一命题中，主文建议只保留 Theorem C.1--C.4；其余细节放入 appendix。以下记号默认全部位于白化、realified 的实 Hilbert 空间。

## 5.1 主文定理：Hilbert-space retention geometry

```latex
\begin{theorem}[Pose-eliminated retention geometry]
\label{thm:continuum-retention}
Let $A\in\mathcal B(\mathcal H_\chi,\mathcal Y)$ and
$B\in\mathcal B(\mathcal H_X,\mathcal Y)$ be the whitened and
realified map and pose tangent operators.  Set
$N:=\overline{\operatorname{Ran}B}$ and let
$A=U_A|A|$ be the polar decomposition.  On
$\mathcal H_A:=\overline{\operatorname{Ran}A^*}=(\ker A)^\perp$, define
\[
R_{\rm geom}:=U_A^*P_{N^\perp}U_A.
\]
Then $R_{\rm geom}$ is a bounded positive contraction and
\[
K_{\rm SLAM}:=A^*P_{N^\perp}A
             =|A|R_{\rm geom}|A|.
\]
Moreover, for every $u$ with $Au\neq0$,
\[
\frac{\langle K_{\rm SLAM}u,u\rangle}
     {\langle A^*Au,u\rangle}
 =\frac{\langle R_{\rm geom}|A|u,|A|u\rangle}
        {\||A|u\|^2}.
\]
Consequently,
\[
\gamma^2:=\inf_{Au\neq0}
\frac{\|P_{N^\perp}Au\|^2}{\|Au\|^2}
=\inf\sigma(R_{\rm geom})
\]
is the largest constant for which
$K_{\rm SLAM}\succeq\gamma^2A^*A$ in quadratic-form order.
\end{theorem}
```

```latex
\begin{remark}[Closure and exact compensation]
\label{rem:closure-pose-range}
The orthogonal projector is taken onto
$\overline{\operatorname{Ran}B}$, not onto a possibly nonclosed range.
Thus $K_{\rm SLAM}u=0$ means that $Au$ can be approximated arbitrarily
well by pose tangents.  It implies the existence of an exact finite-norm
pose increment $h$ satisfying $Au=Bh$ only when
$\operatorname{Ran}B$ is closed and $Au\in\operatorname{Ran}B$.
This distinction disappears in finite dimensions but is essential in the
continuum problem.
\end{remark}
```

```latex
\begin{corollary}[Stable transversality]
\label{cor:stable-transversality}
Let $M:=\overline{\operatorname{Ran}A}$ and
$N:=\overline{\operatorname{Ran}B}$.  If $M\cap N=\{0\}$, then
\[
\gamma=\sin\theta_F(M,N),
\qquad
\gamma>0
\quad\Longleftrightarrow\quad
M+N\ \text{is closed},
\]
where $\theta_F$ is the Friedrichs angle.  Hence trivial exact
intersection is strictly weaker than stable transversality.
\end{corollary}
```

## 5.2 主文反例：有限网格最小角不构成 continuum certificate

```latex
\begin{example}[Zero intersection with vanishing transversality]
\label{ex:vanishing-angle}
Let $\mathcal Y=\ell^2\oplus\ell^2$,
$M=\ell^2\oplus\{0\}$, and
$T e_n=n^{-1}e_n$.  Define
$N=\operatorname{graph}T$,
$Ae_n=(n^{-1}e_n,0)$, and $Bx=(x,Tx)$.
Then $A$ is compact, $\operatorname{Ran}B=N$ is closed, and
$M\cap N=\{0\}$.  Nevertheless,
\[
\operatorname{dist}((e_n,0),N)\le n^{-1},
\qquad \gamma=0.
\]
The $n$-mode truncation has a positive smallest principal angle of order
$n^{-1}$, which converges to zero.  Therefore a positive smallest angle
on every computed grid is not a continuum stability certificate.
\end{example}
```

## 5.3 主文定理：rigid-motion gauge 与 quotient information

```latex
\begin{theorem}[Rigid-motion gauge]
\label{thm:rigid-gauge}
Assume a homogeneous isotropic background and a measurement model that
depends only on relative positions and body-fixed antenna orientations.
For $g=(R,a)\in SE(d)$, let
\[
(g\cdot\chi)(z)=\chi\bigl(R^{-1}(z-a)\bigr),
\qquad
 g\cdot X=(g x_1,\ldots,g x_T).
\]
Whenever the direct scattering problem has a unique outgoing solution,
\[
F(g\cdot\chi,g\cdot X)=F(\chi,X).
\]
If the action and $F$ are differentiable, every Lie-algebra generator
$\xi$ satisfies
\[
A\xi_\chi+B\xi_X=0.
\]
Hence, for an unanchored pose nuisance model,
$\xi_\chi\in\ker K_{\rm SLAM}$.  If the scene is invariant under a
one-parameter subgroup, then the corresponding map generator may vanish;
in that case the generator is a pose-only stabilizer rather than a
nonzero map-gauge direction.
\end{theorem}
```

```latex
\begin{remark}[Gauge breaking]
\label{rem:gauge-breaking}
A finite field of view, a known nonhomogeneous background, fixed receivers,
an absolute phase or position reference, a pose anchor, or a map support
constraint preserves only the subgroup that leaves that additional
structure invariant.  A discrete pixel basis may also break the continuum
gauge approximately because the infinitesimal generator need not belong to
the chosen map space.  All rank and covariance statements are therefore
made either after anchoring or on a chosen horizontal quotient space.
\end{remark}
```

## 5.4 主文定理：半正定 pose prior 的正确 Schur 对象

```latex
\begin{theorem}[Semidefinite pose prior and augmented-space shorting]
\label{thm:semidefinite-prior}
Let $J_X\succeq0$ be bounded and self-adjoint, and define
\[
\bar A u=(Au,0),
\qquad
\bar B h=(Bh,J_X^{1/2}h)
\]
as operators into $\mathcal Y\oplus\mathcal H_X$.  The effective map
information is the closed quadratic form
\[
q_{\rm eff}(u)
 :=\inf_h\bigl(\|Au-Bh\|^2+\|J_X^{1/2}h\|^2\bigr)
 =\|P_{\overline{\operatorname{Ran}\bar B}^{\perp}}\bar A u\|^2.
\]
In finite dimensions, or more generally when the relevant range conditions
and closed-range assumptions hold,
\[
K_{\rm eff}
=A^*A-A^*B(B^*B+J_X)^\dagger B^*A.
\]
If $J_{X,1}\succeq J_{X,0}\succeq0$, then
$K_{\rm eff}(J_{X,1})\succeq K_{\rm eff}(J_{X,0})$.
For $J_X\succ0$, writing $D=BJ_X^{-1/2}$ gives
\[
K_{\rm eff}=A^*(I+DD^*)^{-1}A.
\]
Thus finite-prior retention is a prior-weighted shrinkage in the original
data space; it is an ordinary squared-sine principal-angle spectrum only
after lifting the problem to the augmented data--prior space.
\end{theorem}
```

## 5.5 主文定理：多频 shared compensation

```latex
\begin{proposition}[Shared multifrequency compensation]
\label{prop:shared-frequency}
For a finite frequency set $\Omega$, define the stacked operators
$A_\Omega u=(A_f u)_{f\in\Omega}$ and
$B_\Omega h=(B_f h)_{f\in\Omega}$.  The exactly confounded map space is
\[
\mathcal C_\Omega
=\{u:\exists h,\ A_f u=B_f h\ \text{for every }f\in\Omega\}
=\pi_\chi\ker[A_\Omega,-B_\Omega].
\]
In finite dimensions,
\[
\dim\mathcal C_\Omega
=n_\chi-\operatorname{rank}[A_\Omega,-B_\Omega]
 +\operatorname{rank}B_\Omega.
\]
If, after removing a fixed gauge, the entries of $A_f$ and $B_f$ depend
real-analytically on $f$ and one frequency tuple yields full column rank of
$[A_\Omega,-B_\Omega]$, then almost every tuple of the same cardinality
yields full column rank.  The witness assumption is indispensable.
\end{proposition}
```

```latex
\begin{remark}[No normalized monotonicity claim]
Adding frequency blocks makes the absolute quadratic form
\[
u^*K_{\rm eff}(\Omega)u
=\min_h\left(\sum_{f\in\Omega}\|A_fu-B_fh\|^2+h^*J_Xh\right)
\]
nondecreasing.  It does not make the generalized retention eigenvalues
coordinatewise nondecreasing because $K_{\rm IS}(\Omega)$ changes at the
same time.  Frequency repetition, common scaling, a narrow-band limit, and
scene symmetries provide explicit exceptions to any unconditional
``more bandwidth is better'' statement.
\end{remark}
```

## 5.6 主文定理：fixed-rank smoothness 与 rank-event discontinuity

```latex
\begin{theorem}[Smooth and singular rank strata]
\label{thm:rank-strata}
Let $B(t)$ be $C^1$ and have constant finite rank $r$ on an interval.
Then $P_B(t):=B(t)B(t)^\dagger$ is $C^1$ and
\[
\dot P_B
=(I-P_B)\dot B B^\dagger
 +\bigl((I-P_B)\dot B B^\dagger\bigr)^*.
\]
In particular,
\[
\|\dot P_B\|
\le \frac{\|\dot B\|}{\sigma_r(B)}.
\]
At a genuine rank change, orthogonal projectors on the two sides have
different finite ranks and therefore differ in operator norm by one.
Consequently the no-prior projector is generally not even norm-continuous,
whereas a uniformly positive prior $J_X\succeq\alpha I$ replaces the
pseudoinverse projector by a smooth resolvent.
\end{theorem}
```

## 5.7 主文定理：SOM/current-space 与 map-tangent 的合法桥接

```latex
\begin{theorem}[Current-to-map composition]
\label{thm:som-map-bridge}
Let $G_S:\mathcal H_J\to\mathcal Y$ be the current-to-data operator and
$T_\chi:\mathcal H_\chi\to\mathcal H_J$ the derivative of the induced
current with respect to the material map.  If
$A=G_ST_\chi$, then
\[
\operatorname{Ran}A\subseteq\operatorname{Ran}G_S,
\qquad
\overline{\operatorname{Ran}A}
 \subseteq\overline{\operatorname{Ran}G_S}.
\]
Moreover,
\[
\operatorname{Ran}A=\operatorname{Ran}G_S
\quad\Longleftrightarrow\quad
\mathcal H_J=\operatorname{Ran}T_\chi+\ker G_S.
\]
The right singular vectors of $G_S$ and $A$ belong to different parameter
spaces and therefore have no intrinsic Grassmann distance unless an
explicit identification through $T_\chi$ is supplied.
\end{theorem}
```

```latex
\begin{remark}[Order of reduction and pose elimination]
For a fixed current-space projector $P_r$, replacing $T_\chi$ by
$P_rT_\chi$ and then applying pose elimination is a well-defined reduced
model.  However, recomputing a dominant spectral subspace before and after
pose elimination generally gives different projectors.  Exact commutation
requires a common reducing decomposition for the relevant normal operators;
a product of noncommuting projectors is not the projector onto their
intersection.
\end{remark}
```

## 5.8 主文 statistical limitation 段落

```latex
\paragraph{Statistical scope.}
The operators above describe local tangent information.  If the true pose is
fixed across repeated experiments while a deterministic pose penalty is used
in estimation, the linearized map covariance is a sandwich covariance, not
in general $K_{\rm eff}^{-1}$.  The inverse Schur information is the marginal
map covariance under a correctly specified Gaussian random-pose model, or
when the auxiliary measurements that generate $J_X$ are themselves included
in the repeated-sampling experiment.  Under local asymptotic normality,
positive efficient-information gap, fixed quotient gauge, and initialization
inside the same local basin, reciprocal retention eigenvalues describe
relative covariance inflation in $K_{\rm IS}$-whitened coordinates.  No claim
is made about global uniqueness, phase aliases, cycle skipping, convergence
from arbitrary initialization, or real-system performance.
```

## 5.9 建议替换当前 manuscript 的 limitations

```latex
\section{Limitations}
All continuum statements in this paper are conditional on the explicitly
stated function spaces, closed-range or augmented-space formulations, gauge
quotient, and spectral-gap assumptions.  The exact finite-dimensional rank
identities do not by themselves imply convergence of ranks or principal
angles under grid refinement.  In particular, compact map-to-data operators
may have trivial intersection with the pose range while their stable
transversality constant is zero.  Positive continuum lower bounds therefore
require a task subspace, a regularization window, or a uniform frame/inf--sup
condition.

The rigid-motion result is an invariance theorem only for a homogeneous
background with jointly transformed scene and body-fixed trajectory.  Known
background heterogeneity, anchors, fixed sensors, absolute phase references,
support constraints, and discretization can break this gauge.  The resulting
small but nonzero discrete gauge residual is not evidence of a physical
continuum gauge defect unless a convergent discretization analysis is given.

The multifrequency genericity statement is finite-dimensional and conditional
on an analytically nondegenerate witness.  It does not yield a geometry-only
lower bound on the smallest principal angle, nor does it imply monotonicity of
normalized retention.  Exact shared-pose frequency selection is not generally
submodular; any submodular objective used here is explicitly a surrogate or a
lower bound.

The full-wave sensitivity bounds hold only on compact parameter sets with a
uniform resolvent margin, sensor standoff, and either a positive pose prior or
a fixed-rank pose-tangent gap.  Near a scattering resonance, a collision, or a
rank event, these constants may diverge.  A sampled perturbation test is not a
uniform nonlinear robustness certificate.

Finally, the SOM singular subspaces are current-space objects, whereas the
retention spectrum is defined from map and pose tangents.  Their relation is
only through the explicit composition $A=G_ST_\chi$ and the declared
reduction error.  The literature review is retrieval-bounded and supports no
priority claim.
```

## 5.10 Appendix proof text：finite-dimensional generalized Schur condition

```latex
\begin{lemma}[Range condition for the pseudoinverse Schur complement]
\label{lem:schur-range}
Let
\[
\mathcal I=
\begin{bmatrix}
A^*A & A^*B\\
B^*A & B^*B+J_X
\end{bmatrix},
\qquad J_X\succeq0,
\]
in finite dimensions.  Then
\[
\operatorname{Ran}(B^*A)
\subseteq\operatorname{Ran}(B^*B+J_X).
\]
Consequently the generalized Schur complement is
\[
A^*A-A^*B(B^*B+J_X)^\dagger B^*A\succeq0.
\]
\end{lemma}

\begin{proof}
For positive semidefinite matrices,
\[
\ker(B^*B+J_X)=\ker B\cap\ker J_X.
\]
If $z$ lies in this kernel, then $Bz=0$ and hence
$z^*B^*A=0$.  Therefore
$\ker(B^*B+J_X)\subseteq\ker(A^*B)$, which is equivalent in finite
dimensions to the stated range inclusion.  The Schur complement is the
minimum of
$\|Au-Bh\|^2+h^*J_Xh$ over $h$, and is therefore positive semidefinite.
\end{proof}
```

## 5.11 Appendix proof text：local robust eigenvalue remainder

```latex
\begin{theorem}[Second-order robust eigenvalue expansion]
\label{thm:robust-eigenvalue}
Let $K:X\mapsto K(X)$ be a $C^2$ self-adjoint matrix field.  Suppose
$\lambda_r(K(X_0))$ is simple and separated from the rest of the spectrum by
$\delta>0$ throughout a ball of radius $\varepsilon_0$.  If
$\|DK\|\le M_1$ and $\|D^2K\|\le M_2$ in that ball, then the eigenvalue
Hessian is bounded by
\[
L_r\le M_2+\frac{2M_1^2}{\delta}.
\]
For every norm ball $\{h:\|h\|\le\varepsilon\}$ with
$\varepsilon\le\varepsilon_0$,
\[
\min_{\|h\|\le\varepsilon}\lambda_r(K(X_0+h))
=\lambda_r(K(X_0))
 -\varepsilon\|\nabla\lambda_r(X_0)\|_*
 +R_2,
\qquad
|R_2|\le\frac{L_r}{2}\varepsilon^2.
\]
\end{theorem}
```

---

# 6. 有数学依据的设计与实验修改

本节只列入能够由前述定理支持、并且能明确写出失败触发条件的修改。这里的“证书”均指相应线性化、task space 或 compact parameter ball 内的证书，不指真实硬件端到端保证。

## 6.1 先把 continuum 问题限制到可证 task window

**修改。** 不再报告全像素空间的“最小 retention”作为 continuum 指标。先固定一个有限维 task space $E_r\subset\mathcal H_\chi$，例如：

- 给定 Fourier band；
- 给定 spline/FEM coarse space；
- $K_{\rm IS}$ 的一个与零谱分开的 spectral window；
- 任务导数的 span，例如目标质心、低阶形状或区域平均 contrast。

在 $E_r$ 上报告

$$
\gamma_r
:=\inf_{u\in E_r,\ Au\neq0}
\frac{\|P_{N^\perp}Au\|}{\|Au\|},
$$

$$
a_r:=\inf_{u\in E_r,\|u\|=1}\|Au\|,
\qquad
\lambda_{\min}(K_{\rm eff}|_{E_r}).
$$

**Mathematical basis.** If $a_r>0$ and the discrete operators converge in operator norm on $E_r$, then the retention values $\rho_i$ and the associated spectral projectors converge by ordinary finite-dimensional perturbation theory; the vanishing spectral tail of the compact operator need not be included in the same infimum.

**失败触发。** $a_r$ 接近离散误差或噪声 floor；task window 与其余谱的 gap 闭合；$r$ 随网格增长但没有 uniform frame bound。

## 6.2 Gauge 必须在信息计算之前显式 quotient/anchor

构造 joint gauge matrix

$$
G=\begin{bmatrix}G_\chi\\G_X\end{bmatrix},
\qquad
[A\ B]G\approx0.
$$

有两种合法实现：

1. **anchor formulation：** 固定一个 pose 或加入绝对 reference，重新构造 $B,J_X$；
2. **horizontal formulation：** 选取 $H$ 使 $H^TH=I$、$H^TG=0$，在 coordinates $\delta\theta=Hz$ 中计算 joint FIM，再 short 到 map quotient。

不得把 machine-rank tolerance 自动当作 gauge detector。应分别报告：

$$
r_{\rm gauge}
=\frac{\|[A\ B]G\|}
       {\|A G_\chi\|+\|B G_X\|},
$$

map-basis representation residual，以及 gauge-breaking model terms。

**复杂度。** 若 gauge 维数为 $g\le6$，QR/SVD 成本 $O((n+q)g^2)$，远低于一次 forward/adjoint solve。

**失败触发。** 背景或传感器不满足所声明 symmetry；generator 不在 map basis；anchor 同时被误计入 $J_X$ 与 quotient；径向 scene 的 rotation generator 为零却仍被计作 map gauge。

## 6.3 合法的 SOM--pose two-fold reduction

设 $P_r^J=V_rV_r^*$ 是 Chen $G_S$ current-space dominant projector，且

$$
A=G_ST_\chi,
\qquad
A_r:=G_SP_r^JT_\chi.
$$

采用下列层级算法，而不是三个 projector 的“八象限”乘积。

### 算法 6.1（Current reduction followed by pose defect analysis）

1. 对 $G_S$ 做 truncated/randomized SVD，选择具有外部 singular-value gap 的 $V_r$；
2. 形成 matrix-free $A_r=G_SV_r(V_r^*T_\chi)$；
3. 白化并 realify $A_r,B$；
4. 对 $A_r,B$ 做 rank-revealing QR，得 $Q_A,Q_B$；
5. 对小矩阵 $C=Q_B^TQ_A$ 做 SVD；无 prior retention 为
   $\rho_i=1-\sigma_i(C)^2$；
6. 有 prior 时在 pose-sized system $B^TB+J_X$ 上求解，不显式形成 map-sized dense Schur matrix；
7. 同时报告 current truncation error与 pose-elimination error，不能只报最终重建图。

若

$$
\varepsilon_A:=\|G_S(I-P_r^J)T_\chi\|,
$$

则

$$
\|K_{\rm IS}-K_{{\rm IS},r}\|
\le(2\|A\|+\varepsilon_A)\varepsilon_A.
$$

在无先验且 $\sigma_{\min}^+(B)\ge\beta_0$、$A$ 的 selected data subspace 有 gap $\delta_A$ 时，pose-retention projector 的附加误差由

$$
O\!\left(\frac{\varepsilon_A}{\delta_A}\right)
$$

控制；若同时近似 $B$，再加入 $O(\varepsilon_B/\beta_0)$。

**复杂度。** 设 data dimension $m$、map dimension $n$、pose dimension $q$、current rank $r$：

- randomized SVD：约 $O(s\,C_{G_S}+mr^2)$，$s$ 为少量 power iterations；
- $Q_B$：$O(mq^2)$；
- $Q_A^TQ_B$：$O(mrq)$；
- 小 SVD：$O(\min(r,q)^2\max(r,q))$；
- finite-prior solve：稠密为 $O(q^3)$，motion graph sparse Cholesky 取决于 graph treewidth，通常远小于 $n^3$。

**失败触发。** $G_S$ cutoff 没有 gap；$T_\chi$ 把重要 map directions 映入 $V_r^\perp$；$B$ rank event；current split随 trajectory adaptive recompute而次序被改变；把 $P_SP_D$ 当作 exact intersection projector。

## 6.4 多频选择：优化“共享补偿矛盾”，不使用带宽口号

对候选频率集 $S$，精确 objective 可以取

$$
f_{\rm exact}(S)
=\log\det\bigl(K_{\rm eff}(S)+\lambda I_E\bigr),
$$

其中 $I_E$ 只作用在固定 task space。该 objective 一般不 submodular，因为同一 $h$ 必须同时解释所有 blocks。

可使用两种有明确语义的替代：

1. **exact but non-submodular search：** greedy 只作为 heuristic，每一步重新求 shared Schur；必须与 swap/local-search 后处理比较；
2. **certified lower-bound surrogate：** 独立消元每个 block，

   $$
   \underline K(S)
   :=\sum_{f\in S}A_f^*P_{\overline{\operatorname{Ran}B_f}^{\perp}}A_f
   \preceq K_{\rm SLAM}(S),
   $$

   对

   $$
   \underline f(S)=
   \log\det(\lambda I_E+\underline K(S))
   $$

   使用标准 PSD-logdet submodular greedy。

第二个 objective保守地允许每个频率使用不同 pose compensation，因此不会夸大共享频率价值。它的 greedy $(1-1/e)$ 保证只针对 surrogate，不自动转移给 exact objective。

每个新增频率应报告：

$$
\Delta K_{\rm abs},
\quad
\Delta\rho,
\quad
\operatorname{rank}[A_S,-B_S],
\quad
\sigma_{\min}([A_S,-B_S]|_{H}),
\quad
\text{SNR/energy cost}.
$$

**失败触发。** duplicate blocks；$A_f,B_f$ 共同尺度复制；频率间 whitening 不一致；窄带导致 analytic minors接近零；实际 pose/calibration并非 shared parameter。

## 6.5 半正定 motion prior：保留 sparse joint solve，不先造 dense pseudoinverse

对于 odometry/pose-graph Laplacian $J_X=L_G\otimes Q^{-1}$：

1. 显式保留 global gauge nullspace；
2. 通过 anchor、nullspace QR 或 sparse constrained solve处理；
3. 用 variational action计算

   $$
   K_{\rm eff}u
   =A^*\bigl(Au-Bh_u\bigr),
   $$

   其中

   $$
   (B^*B+J_X)h_u=B^*Au
   $$

   在 quotient 上求最小范数解；
4. 不显式构造 $(B^*B+J_X)^\dagger$。

**复杂度。** 每次 $K_{\rm eff}$ matvec 为一次 $A,A^*$ 与一次 pose graph linear solve。若 graph 为链/固定带宽，factorization近似 $O(Tq_0^3)$；一般稀疏图由 fill-in/treewidth 决定。

**失败触发。** $B^*A$ 在数值上含落在 unhandled gauge nullspace 的分量；map-pose prior correlation 被遗漏；同一 anchor既作为硬约束又作为 prior重复计数。

## 6.6 Rank-event detector：阈值、hysteresis 与可验证条件

令 computed singular values 为 $\widehat\sigma_i$，并有 deterministic perturbation bound

$$
\|\widehat B-B\|\le\eta_B.
$$

使用相对阈值

$$
\tau_{\rm on}=c_{\rm on}\max\{\eta_B,\epsilon_{\rm mach}\|B\|\},
\qquad
\tau_{\rm off}=c_{\rm off}\max\{\eta_B,\epsilon_{\rm mach}\|B\|\},
$$

其中 $c_{\rm on}>c_{\rm off}>1$。一条 singular direction 只有连续 $L$ 帧高于 $\tau_{\rm on}$ 才加入，连续 $L$ 帧低于 $\tau_{\rm off}$ 才删除。由 Weyl inequality：

- 若 true $\sigma_i(B)>\tau_{\rm on}+\eta_B$，则该方向不会漏检；
- 若 true $\sigma_i(B)<\tau_{\rm off}-\eta_B$，则该方向不会误保留；
- 位于中间 uncertainty strip 的方向必须标为 unresolved，不应强制分类。

对 generalized retention cluster，若 estimated perturbation bound为 $\eta_R$，仅当两个 clusters 的 gap $>2\eta_R$ 时分开追踪；否则合并成一个 projector cluster。

**失败触发。** 没有 $\eta_B$ 上界却宣称误报率；阈值用绝对单位且跨频率/whitening比较；rank event附近继续使用 $\dot P_B$；把 repeated zero cluster中的 eigenvectors逐个配对。

## 6.7 Robust trajectory：用 margin-limited trust region

在候选 trajectory $X$ 处同时估计

$$
m_M:=\sigma_{\min}(M),
\qquad
m_B:=\sigma_{\min}^+(B),
\qquad
m_\lambda:=\operatorname{gap}(\lambda_r),
\qquad
m_c:=\operatorname{dist}(X,\partial\mathcal X_{\rm collision}).
$$

选择 trust radius

$$
\boxed{
\Delta_{
m tr}
\le c\min\left
\{
\frac{m_M}{L_M},
\frac{m_B}{L_B},
\frac{m_\lambda}{2L_K},
 m_c
\right\},
\quad 0<c<1.
}
$$

只有在该球内才使用 smooth first-/second-order robust surrogate。若任一 margin无法下界：

- 切换到 $J_X\succeq\alpha I$ 的 regularized objective；或
- 直接优化 variational residual

  $$
  \min_{\Delta X\in\mathcal U}
  \min_h\|A(X+\Delta X)u-B(X+\Delta X)h\|^2+h^*J_Xh;
  $$

- 使用 bundle/trust-region nonsmooth method，并把 event两侧都纳入模型。

不可把随机方向扫描得到的最大 quotient 当作 $L_K$。要称为 certificate，$L_M,L_B,L_K$ 必须来自 interval arithmetic、analytic derivative bounds或覆盖整个 compact ball 的 verified computation。

## 6.8 轨迹评估必须是 Pareto 而非预设排序

每个 trajectory 至少输出六个互不替代的指标：

$$
\begin{aligned}
I_{\rm abs}&=\lambda_r(K_{\rm eff}|_E)
\quad\text{或}\quad \log\det(K_{\rm eff}|_E+\lambda I),\\
I_{\rm ret}&=\lambda_{\min}(R_E)
\quad\text{或}\quad \operatorname{tr}R_E,\\
I_{\rm rob}&=I_{\rm abs}-\varepsilon\|\nabla I_{\rm abs}\|_*,\\
d_{\rm stand}&=\min_{t,z\in D}\|x_t-z\|,\\
C_{\rm energy}&=\sum_t E_t,\\
C_{\rm motion}&=\text{path length/turning/acceleration cost}.
\end{aligned}
$$

报告 nondominated Pareto set，而不是“circle 必然优于 line”。若必须给 scalar score，应先声明权重与单位，并做权重敏感性分析。

## 6.9 统计实验必须拆成三个 protocol

### Protocol F：fixed nuisance frequentist

- 固定真实 $X_0$；
- 只重采样 data noise；
- estimator中的 pose penalty是 deterministic；
- 比较 empirical covariance 与 sandwich target，而非 $K_{\rm eff}^{-1}$。

### Protocol B：Bayesian random nuisance

- 每个 trial 从声明的 $X\sim N(X_0,C_X)$ 生成 pose；
- data conditional model与 prior一致；
- 比较 Bayesian posterior/marginal covariance 与 $K_{\rm eff}^{-1}$，同时报告 posterior bias。

### Protocol A：auxiliary-sensor repeated sampling

- 真实 pose固定或随机均可，但 IMU/odometry measurements在每 trial重采样；
- 把辅助 measurement score纳入 Hessian与 score variance；
- 此时信息恒等式才可能恢复。

三种 protocol 不能在同一 covariance 图中混用。每个 protocol还应分别报告：basin success rate、conditioned-on-success covariance、all-trial error（含 failures）、bias、coverage 与 multimodality diagnostics。

## 6.10 连续/离散验证的最低升级要求

仅做 $N=16,24,32,40$ 的趋势图不够。至少加入：

1. nested 或具有明确 interpolation/restriction 的 approximation spaces；
2. $A_h\to A$、$B_h\to B$ 的 operator-norm或 collectively compact convergence论证；
3. 固定 task window 的 projector gap；
4. consistency error、quadrature error和 map-basis representation error分开；
5. 对 $\theta_{\min,h}$ 同时报告下界证书或承认其可能趋零；
6. forward与inverse使用不同 mesh/quadrature/solver设置，避免 inverse crime；
7. rank theorem用 numerical-algebra tolerance，physical SOM cutoff另行报告。

---

# 7. 文献 claim ledger

## 7.1 检索边界

本节结合所附 targeted prior-art addendum、当前有限维 manuscript 的 reference list，以及针对 P1--P10 的补充检索。检索优先 publisher page、DOI、arXiv 和正式专著；但它不是系统综述，也没有完整读取所有付费全文。下表中的“未检得”只表示本次检索边界内未检得同构命题，不表示文献不存在。

证据强度：

- **A：** 已核对原论文/正式出版元数据，且对象与命题直接相关；
- **B：** 已核对摘要、引言或可靠元数据，但没有完整核对全文全部定理；
- **C：** 仅作为邻近方向线索，不能支撑 novelty conclusion。

## 7.2 逐命题 ledger

| 本工作命题 | 最接近已有结果 | 相同点 | 不同点 | 仍可主张的窄贡献 | 证据 |
|---|---|---|---|---|---|
| C1. 用 polar decomposition 定义 continuum retention，而不用 $K_{\rm IS}^{-1/2}$ | Friedrichs/Dixmier angles；Bouldin 关于 closed range 与 subspace angle；Björck--Golub principal angles | 都研究闭子空间夹角、closed sum、投影与稳定性 | 本文把两子空间指定为 full-wave map tangent closure 与 physical-pose tangent closure，并写成 $|A|R_{\rm geom}|A|$ | 可主张为该 inverse-scattering-SLAM 模型的严格 specialization/synthesis；抽象 angle theorem 不是新结果 | A |
| C2. 交集为零不保证稳定横截，并给 compact-$A$ 最小反例 | Hilbert 子空间 closed-sum theory | 同一 abstract obstruction：angle 可为零而交集平凡 | 本文反例同时包含 compact map operator、closed pose range 与 finite-truncation false stability | 可主张为针对离散 inverse-scattering retention interpretation 的反例与警告 | A |
| C3. $SE(2)/SE(3)$ 同时作用于 map 与 trajectory 产生 gauge tangent | SLAM observability/gauge；bundle-adjustment datum invariance；scattering equation 的 Euclidean covariance | 都有 global rigid-frame redundancy和 Jacobian null directions | 地图为 volumetric contrast field，forward map是 nonlinear Helmholtz/Lippmann--Schwinger，而不是 landmarks/reprojection | 可主张为该连续 full-wave map/trajectory action、stabilizer 与 quotient FIM 的显式 formulation | A/B |
| C4. Schur-eliminated pose loss 与 principal-angle retention | generalized Schur complement、equivalent Fisher information、canonical correlations | nuisance elimination 与 principal angles 均是成熟工具 | 对象是 whitened/realified map tangent与 real physical-pose tangent；强调 complex-pose projector错误 | 可主张为针对这一物理 tangent pair 的统一几何解释，不可把 Schur/principal angles 本身称为新 | A |
| C5. 半正定 motion prior 用 augmented-space shorted operator | Albert generalized Schur；Anderson--Trapp shorted operators | 都处理 PSD block operator、range condition与 variational shorting | 本文明确分离 motion-graph gauge、finite shrinkage与 known-pose limit，并落到 map information | 可主张模型专用分解与语义澄清；shorted operator理论本身不新 | A |
| C6. finite-prior retention 是 weighted shrinkage，只在 lifted space 是 squared-sine angles | ridge/Tikhonov filters、canonical correlations、augmented least squares | 都用 prior whitening 与 singular values解释 shrinkage | 本文将其解释为 map--pose tangent retention并给 Loewner/equality conditions | 可主张专用 theorem package；不能声称发现新的 canonical correlation | A |
| C7. fixed-rank $B\mapsto BB^\dagger$ 可微，rank event跳变 | Golub--Pereyra pseudoinverse/projector derivatives；Kato；Davis--Kahan | constant-rank smoothness、spectral-cluster perturbation均成熟 | 本文把条件显式连接到 pose-defect rank event、online detector和 $K_{\rm eff}$ derivative | 可主张 inverse-scattering-SLAM 的 rank-stratified consequence与 detector design；导数公式本身不新 | A |
| C8. Born empty scene 是二阶 bilinear operator-uncertainty regime；$s\neq0$ no-prior projector scale-invariant而 $s=0$ 奇异 | bilinear inverse-problem identifiability；blind calibration；Born approximation误差分析 | 都涉及 bilinear ambiguity、scaling group或 Born model | 本文固定 physical Born pose tangent $B_B(s)=sB_1$，区分 complete pair 与 mixed isolation control，并给 finite-prior $O(s^2)$ loss | 本次检索未发现同一 full-wave/Born-pose rank-stratum statement；可作为窄命题，仍需投稿前全文查重 | B |
| C9. 多频同一 pose compensation的 exact kernel与 analytic genericity | multifrequency inverse scattering stability；bilinear generic identifiability | 多频常提高稳定性；analytic-minor genericity是成熟论证方式 | 本文约束所有 blocks共享同一 physical $h$，并研究 map/pose tangent intersection | 可主张 conditional theorem；必须保留“存在一个 witness”条件，不能声称宽带普遍 transverse | A/B |
| C10. exact shared-Schur frequency objective一般不 submodular | submodular Gaussian sensor placement与 PSD log-det design | 都研究信息型 set functions | shared nuisance使 blocks产生 complementarity，破坏普通 additive PSD结构 | 可主张最小反例与 conservative submodular lower-bound surrogate；不是对所有 OED 的否定 | A |
| C11. full-wave resolvent margin控制 $DA,DB,D^2K$，近共振同时放大信息与敏感度 | Lippmann--Schwinger resolvent calculus；inverse-scattering stability；operator perturbation | 都使用 resolvent identity与 Green-function regularity | 本文把常数继续传递到 pose-eliminated map information，显式出现 $m_0,d_0,\beta_0$ 或 $\alpha$ | 可主张 conditional bound chain；当前未完成一般 moving-grid $D^2B=D^3F$ continuum bound | A/B |
| C12. robust trajectory 的 first-order dual-norm penalty与 rank-event trust region | Fisher-information trajectory synthesis；active sensing/OED；eigenvalue perturbation | 局部 E-/D-optimal design与 trajectory sensitivity已有大量工作 | 本文 objective同时含 map information、pose elimination、rank strata与 execution perturbation | 可主张组合 formulation、remainder条件和 failure trigger；不能声称信息型 trajectory design本身新 | A |
| C13. Chen SOM current modes 与 map tangent modes只通过 $A=G_ST_\chi$ 联系 | Chen Gs-SOM/TSOM；Zhong--Chen Twofold SOM | 都使用 $G_S,G_D$ 的 current-space SVD和 reduced induced current | 本文研究 map-tangent/pose-tangent geometry，不是 current recovery；给 range/equality/commutation条件 | 可主张桥接定理与不交换反例；SOM decomposition本身不新 | A |
| C14. joint inverse scattering + unknown transmitter/geometry 已有前例 | Karthik--Ghosh 2023；radar autofocus；joint SAR image/phase-error estimation | 都联合估计 image/contrast与未知 geometry/phase nuisance | 现有邻居多为 DNN、SAR phase error、blind deconvolution或有限几何 map，不是本文的 volumetric full-wave tangent retention | 只能主张更窄的 tangent-information synthesis；不得主张“首次联合反演与定位” | A/B |
| C15. RF-SLAM map-vs-pose information已有邻近工作 | Deutschmann et al. distributed MIMO SLAM PCRB；Su et al. TDOA calibration observability；pose-graph CRB | 都使用 FIM/PCRB、joint map/pose或 rank conditions | map为 specular surfaces、landmarks或 sensor calibration，不是 volumetric contrast | 可主张物理 map class与 retention geometry不同；FIM/CRB for SLAM本身不新 | A/B |
| C16. fixed-nuisance penalized estimator需要 sandwich covariance | White misspecified MLE；Jennrich nonlinear least squares | 都区分 Hessian与 score covariance、给 asymptotic sandwich | 本文把该语义用于 pose-penalized full-wave map estimate，并与 Bayesian marginal区分 | 可主张 protocol correction与模型专用 covariance formula；sandwich theory本身不新 | A |

## 7.3 关键文献目录与精确覆盖对象

### Hilbert geometry、Schur 与 perturbation

1. **R. Bouldin**, “The Product of Operators with Closed Range,” *Tohoku Mathematical Journal*, 25(3), 1973, pp. 359--363. DOI: [10.2748/tmj/1178241337](https://doi.org/10.2748/tmj/1178241337).  覆盖 closed range 与 subspace angle；不含 inverse scattering。
2. **A. Björck and G. H. Golub**, “Numerical Methods for Computing Angles Between Linear Subspaces,” *Mathematics of Computation*, 27(123), 1973, pp. 579--594. DOI: [10.1090/S0025-5718-1973-0348991-3](https://doi.org/10.1090/S0025-5718-1973-0348991-3). 覆盖 principal angles/CS computation。
3. **A. Albert**, “Conditions for Positive and Nonnegative Definiteness in Terms of Pseudoinverses,” *SIAM Journal on Applied Mathematics*, 17, 1969, pp. 434--440. DOI: [10.1137/0117041](https://doi.org/10.1137/0117041). 覆盖 generalized Schur range condition。
4. **W. N. Anderson and G. E. Trapp**, “Shorted Operators. II,” *SIAM Journal on Applied Mathematics*, 28(1), 1975, pp. 60--71. DOI: [10.1137/0128007](https://doi.org/10.1137/0128007). 覆盖 Hilbert-space positive shorted operator。
5. **G. H. Golub and V. Pereyra**, “The Differentiation of Pseudo-Inverses and Nonlinear Least Squares Problems Whose Variables Separate,” *SIAM Journal on Numerical Analysis*, 10(2), 1973. DOI: [10.1137/0710036](https://doi.org/10.1137/0710036). 覆盖 constant-rank pseudoinverse/projector derivative与 variable projection。
6. **T. Kato**, *Perturbation Theory for Linear Operators*, Springer. DOI: [10.1007/978-3-642-66282-9](https://doi.org/10.1007/978-3-642-66282-9). 覆盖 analytic operator families、isolated spectral clusters与 Riesz projectors。
7. **C. Davis and W. M. Kahan**, “The Rotation of Eigenvectors by a Perturbation. III,” *SIAM Journal on Numerical Analysis*, 7(1), 1970, pp. 1--46. DOI: [10.1137/0707001](https://doi.org/10.1137/0707001). 覆盖 invariant-subspace perturbation与 sin-theta bounds。

### Inverse scattering、SOM、Born 与 multifrequency

8. **X. Chen**, *Computational Methods for Electromagnetic Inverse Scattering*, Wiley--IEEE Press, 2018, ISBN 9781119311980. 覆盖 Lippmann--Schwinger/contrast-source、Gs-SOM、Twofold SOM、FFT-SOM、MUSIC/MSR、CRB与 regularization；SOM 的子空间对象是 induced-current/current-to-data operator。
9. **X. Chen**, “Subspace-Based Optimization Method for Solving Inverse-Scattering Problems,” *IEEE Transactions on Geoscience and Remote Sensing*, 2010. DOI: [10.1109/TGRS.2009.2025122](https://doi.org/10.1109/TGRS.2009.2025122). 覆盖 Gs-SOM current decomposition；不含 pose nuisance。
10. **Y. Zhong and X. Chen**, “Twofold Subspace-Based Optimization Method for Solving Inverse Scattering Problems,” *Inverse Problems*, 25, 085003, 2009. DOI: [10.1088/0266-5611/25/8/085003](https://doi.org/10.1088/0266-5611/25/8/085003). 覆盖 $G_S/G_D$ 两层 current-space reduction。
11. **M. L. Diong, A. Roueff, P. Lasaygues, and A. Litman**, “Impact of the Born Approximation on the Estimation Error in 2D Inverse Scattering,” *Inverse Problems*, 32, 065006, 2016. DOI: [10.1088/0266-5611/32/6/065006](https://doi.org/10.1088/0266-5611/32/6/065006). 覆盖 Born-based estimator bias/error与 exact-model CRB；不处理 unknown pose tangent。
12. **M. N. Entekhabi and V. Isakov**, "On Increasing Stability in the Two-Dimensional Inverse Source Scattering Problem with Many Frequencies," *Inverse Problems*, 34(5), 055005, 2018. DOI: [10.1088/1361-6420/aab465](https://doi.org/10.1088/1361-6420/aab465); [arXiv:1712.08696](https://arxiv.org/abs/1712.08696). Covers multifrequency increasing stability; it is not a theorem on shared-pose retention genericity.
13. **F. Weidling and T. Hohage**, "Variational Source Conditions and Stability Estimates for Inverse Electromagnetic Medium Scattering Problems," *Inverse Problems and Imaging*, 11(1), 203--220, 2017. DOI: [10.3934/ipi.2017010](https://doi.org/10.3934/ipi.2017010); [arXiv:1512.06586](https://arxiv.org/abs/1512.06586). Covers conditional stability and regularization rates for electromagnetic medium inversion; it does not contain pose Schur geometry.
14. **A. Capozzoli, C. Curcio, and A. Liseno**, "Singular Value Optimization for Multifrequency Multimonostatic Inverse Scattering over Circular Domains under the Born Approximation," *URSI Radio Science Letters*, 4, 1--5, 2022. DOI: [10.46620/22-0065](https://doi.org/10.46620/22-0065). Covers Born acquisition and singular-value design; spectral acquisition design is therefore established prior art.

### Unknown geometry、autofocus、bilinear identifiability 与 SLAM

15. **Y. Li, K. Lee, and Y. Bresler**, "Identifiability in Bilinear Inverse Problems With Applications to Subspace or Sparsity-Constrained Blind Gain and Phase Calibration," *IEEE Transactions on Information Theory*, 63(2), 822--842, 2017. DOI: [10.1109/TIT.2016.2637933](https://doi.org/10.1109/TIT.2016.2637933); precursor [arXiv:1501.06120](https://arxiv.org/abs/1501.06120). Covers transformation-group identifiability; it does not contain Helmholtz pose tangents.
16. **H. Mansour, D. Liu, U. S. Kamilov, and P. T. Boufounos**, “Sparse Blind Deconvolution for Distributed Radar Autofocus Imaging,” *IEEE Transactions on Computational Imaging*, 4(4), 2018. DOI: [10.1109/TCI.2018.2875375](https://doi.org/10.1109/TCI.2018.2875375); arXiv:1805.03269. 覆盖 antenna-position ambiguity 与 multichannel blind deconvolution。
17. **N. Önhon and M. Çetin**, “A Sparsity-Driven Approach for Joint SAR Imaging and Phase Error Correction,” *IEEE Transactions on Image Processing*, 21(4), 2012. DOI: [10.1109/TIP.2011.2179056](https://doi.org/10.1109/TIP.2011.2179056). 覆盖 joint image/phase-error estimation。
18. **T. Scarnati and A. Gelb**, “Joint Image Formation and Two-Dimensional Autofocusing for Synthetic Aperture Radar Data,” *Journal of Computational Physics*, 374, 2018. DOI: [10.1016/j.jcp.2018.07.059](https://doi.org/10.1016/j.jcp.2018.07.059). 覆盖 joint SAR image formation 与 2D phase correction。
19. **Karthik and Ghosh**, “A Scalable Deep Learning Model for Simultaneous Reconstruction and Transmitter Localization in Inverse Scattering,” PIERS 2023. DOI: [10.1109/PIERS59004.2023.10221374](https://doi.org/10.1109/PIERS59004.2023.10221374). 直接覆盖 contrast reconstruction + transmitter localization；因此一般性的“联合 inverse scattering 与 localization”不是可主张 novelty。
20. **C. Deutschmann, S. Li, F. Meyer, and E. Leitinger**, “Posterior Cramér--Rao Bounds on Localization and Mapping Errors in Distributed MIMO SLAM,” Asilomar 2025. DOI: [10.1109/IEEECONF67917.2025.11443670](https://doi.org/10.1109/IEEECONF67917.2025.11443670); [arXiv:2506.19957](https://arxiv.org/abs/2506.19957). 覆盖 mobile transceiver 与 specular-surface map 的 joint PCRB。
21. **S. Su, X. Kong, S. Sukkarieh, and S. Huang**, “Necessary and Sufficient Conditions for Observability of SLAM-Based TDOA Sensor Array Calibration and Source Localization,” *IEEE Transactions on Robotics*, 37(5), 2021. DOI: [10.1109/TRO.2021.3069140](https://doi.org/10.1109/TRO.2021.3069140). 覆盖 Jacobian/FIM rank observability 与 degenerate configurations。
22. **G. P. Huang, A. I. Mourikis, and S. I. Roumeliotis**, “Observability-Based Rules for Designing Consistent EKF SLAM Estimators,” *International Journal of Robotics Research*, 2010. DOI: [10.1177/0278364909353640](https://doi.org/10.1177/0278364909353640). 覆盖 SLAM 的 global translation/rotation unobservable directions与 linearization consistency。
23. **Z. Zhang, G. Gallego, and D. Scaramuzza**, “On the Comparison of Gauge Freedom Handling in Optimization-Based Visual-Inertial State Estimation,” *IEEE Robotics and Automation Letters*, 3(3), 2018. DOI: [10.1109/LRA.2018.2833152](https://doi.org/10.1109/LRA.2018.2833152). 覆盖 anchor/prior/free-gauge handling；不含 wave map。
24. **B. Triggs, P. McLauchlan, R. Hartley, and A. Fitzgibbon**, “Bundle Adjustment -- A Modern Synthesis,” LNCS 1883, 2000. DOI: [10.1007/3-540-44480-7_21](https://doi.org/10.1007/3-540-44480-7_21). 覆盖 structure/pose Schur、sparsity 与 gauge/datum invariance。
25. **Y. Chen, S. Huang, L. Zhao, and G. Dissanayake**, “Cramér--Rao Bounds and Optimal Design Metrics for Pose-Graph SLAM,” *IEEE Transactions on Robotics*, 37(2), 2021. DOI: [10.1109/TRO.2020.3001718](https://doi.org/10.1109/TRO.2020.3001718). 覆盖 pose-graph FIM、Laplacian与 D/T-optimal metrics。

### Active sensing、submodularity 与统计语义

26. **A. D. Wilson, J. A. Schultz, and T. D. Murphey**, “Trajectory Synthesis for Fisher Information Maximization,” *IEEE Transactions on Robotics*, 30(6), 2014. DOI: [10.1109/TRO.2014.2345918](https://doi.org/10.1109/TRO.2014.2345918). 覆盖 nonlinear-system local FIM trajectory optimization；不含 pose-eliminated map FIM。
27. **A. Krause, A. Singh, and C. Guestrin**, “Near-Optimal Sensor Placements in Gaussian Processes: Theory, Efficient Algorithms and Empirical Studies,” *JMLR*, 9, 2008, pp. 235--284. 稳定链接: [JMLR paper](https://jmlr.org/papers/v9/krause08a.html). 覆盖 mutual-information submodularity 与 greedy guarantee；该保证依赖其 Gaussian/additive结构，不能直接迁移到 shared-nuisance Schur objective。
28. **R. I. Jennrich**, "Asymptotic Properties of Non-Linear Least Squares Estimators," *Annals of Mathematical Statistics*, 40(2), 633--643, 1969. DOI: [10.1214/aoms/1177697731](https://doi.org/10.1214/aoms/1177697731). Covers nonlinear least-squares consistency and asymptotics.
29. **H. White**, “Maximum Likelihood Estimation of Misspecified Models,” *Econometrica*, 50(1), 1982, pp. 1--25. DOI: [10.2307/1912526](https://doi.org/10.2307/1912526). 覆盖 misspecified/quasi-MLE sandwich covariance。

## 7.4 Retrieval-bounded novelty conclusion

最安全的 contribution statement 是：

> We formulate and analyze a whitened and realified map-tangent versus physical-pose-tangent geometry for a nonlinear full-wave volumetric contrast model.  The analysis specializes established tools from Schur complementation, principal-angle geometry, group gauge theory, and operator perturbation to this model, and adds model-specific counterexamples and conditional consequences for Born rank strata, shared-frequency compensation, SOM-to-map reduction, and rank-event-aware trajectory design.  The literature review is retrieval-bounded and no priority claim is made.

不得使用：

- “the first inverse-scattering SLAM method”；
- “the first joint reconstruction and transmitter localization”；
- “no prior work studies FIM/principal angles/subspaces”；
- “wide bandwidth guarantees transversality”；
- “SOM modes are the map observability modes”。

---

# 8. 尚未解决清单：最小下一步、所需引理与当前卡点

## 8.1 P1：真实 Helmholtz discretization 的 retention spectral convergence

**未解决对象。** 对实际 volume-integral/FEM/Nyström discretization，证明固定 task window 上

$$
R_{{\rm geom},h}\to R_{\rm geom}
$$

及相应 spectral projectors/gaps收敛。

**最小下一步。** 固定一个 finite-dimensional $E_r\subset W^{s,2}(D)$，证明

$$
\|(A_hI_h-A)|_{E_r}\|\to0,
\qquad
\|B_hJ_h-B\|\to0,
$$

并给 $\beta=\inf\sigma^+(B)>0$ 或直接证明 pose-range projector gap convergence。

**所需工具。** collectively compact approximation；Atkinson theory；Babuška inf--sup；gap convergence of subspaces；Nyström quadrature error for weakly singular Hankel kernels。

**当前卡点。** 已有网格不是一个明示的 nested operator approximation，且全谱底端 compactness排除了 global uniform lower frame bound。

## 8.2 P2：包含 direct path、phase center 与标定参数的完整 gauge group

**未解决对象。** 当数据含 direct Tx--Rx path、unknown clock/phase offset、body-frame phase center、外参和非均匀背景时，确定完整 stabilizer subgroup及其 Lie algebra。

**最小下一步。** 把 forward model写成

$$
F(\chi,X,c)=F_{\rm sca}(\chi,X,c)+F_{\rm dir}(X,c),
$$

逐项计算 group action，求满足 $F(g\cdot\chi,g\cdot X,g\cdot c)=F(\chi,X,c)$ 的最大 subgroup。

**所需工具。** Lie group actions；isotropy/stabilizer theorem；Maxwell/Helmholtz covariance；calibration graph gauge analysis。

**当前卡点。** 不同硬件的 absolute phase reference和 clock convention不同，不能在抽象层面唯一指定 gauge。

## 8.3 P3：line/circle/arc 的 scene-independent exact intersection formula

**结论。** 在无限维 map class中，该目标一般不可实现：有限 Fourier samples可被任意插值，连续 restricted Fourier operator通常 compact/dense-range，因而 exact/near overlap高度 scene-和 function-class-dependent。

**最强可做下一步。** 固定一个明确的 Fourier/spline subspace $E_r$ 与阵列采样 frame，计算显式 matrix

$$
C_{AB}=Q_B^TQ_A
$$

并从 coherence/frame constants推导

$$
\gamma_r^2\ge1-\|C_{AB}\|^2.
$$

对每类 trajectory用 trigonometric-polynomial zero-set 或 Vandermonde determinant给 witness。

**所需工具。** nonharmonic Fourier frames；Ingham inequalities；prolate-spheroidal concentration；Ewald geometry；analytic determinant/nonvanishing minor。

**当前卡点。** 没有指定 map basis、aperture sampling density和 noise metric时，不存在仅依赖“圆/直线/孔径角”的正统一角界。

## 8.4 P4：具体阵列的 analytic genericity witness

**未解决对象。** 对用户实际 MIMO layout与 pose parameterization，给出一个频率 tuple使 gauge quotient上的

$$
[A_\Omega,-B_\Omega]
$$

满列秩。

**最小下一步。** 在 Born far-field finite basis中选择最小 $F$，symbolically/numerically寻找一个 candidate tuple；然后用 interval arithmetic或 exact algebra证明某个 square minor与零分离。

**所需工具。** analytic Fredholm/minor argument；interval determinant or singular-value verification；array-manifold nondegeneracy lemma。

**当前卡点。** 目前只有数值 singular values，没有 rigorous nonzero minor certificate；full-wave resolvent poles也限制全频 analytic continuation。

## 8.5 P4/P8：exact shared-Schur set function 的 weak submodularity

**未解决对象。** 在何种 quantitative coherence条件下，exact frequency objective有 submodularity ratio $\kappa>0$。

**最小下一步。** 对 finite task space，界定每个新增 block的 shared-compensation cross term，并尝试用

$$
\|B_f^TB_g\|,
\quad
\|A_f^TB_g\|,
\quad
\lambda_{\min}(J_X+\sum B_f^TB_f)
$$

控制 marginal-gain ratio。

**所需工具。** weak submodularity/submodularity ratio；matrix inverse update；restricted strong concavity/smoothness。

**当前卡点。** 两频率 sign-flip反例表明没有无条件正 ratio；必须施加强 cross-frequency incoherence或 strictly positive prior。

## 8.6 P5：无限维 semidefinite motion prior 的 closed-range characterization

**未解决对象。** 对 trajectory function space（例如 $H^1([0,T];\mathfrak{se}(d))$）和 odometry differential operator，给出

$$
\operatorname{Ran}\begin{bmatrix}B\\J_X^{1/2}\end{bmatrix}
$$

闭的必要充分条件。

**最小下一步。** 将 $J_X$ 具体化为 elliptic/graph differential form，证明在 quotient gauge上 Poincaré inequality：

$$
\|h\|_{H_X}
\le C\bigl(\|Bh\|_\mathcal Y+\|J_X^{1/2}h\|\bigr).
$$

**所需工具。** closed-range theorem；coercivity modulo finite-dimensional kernel；Poincaré/Korn inequality；compact perturbation of coercive forms。

**当前卡点。** 抽象 $J_X\succeq0$ 太弱；其 kernel可能无限维，且 $B$ 未必在 kernel上提供 inf--sup。

## 8.7 P6：有噪 Jacobian下 rank-event detector 的概率保证

**未解决对象。** 当 $B$ 本身由 noisy calibration/finite differences估计时，给 hysteresis detector的 finite-sample false alarm/miss概率。

**最小下一步。** 建立

$$
\Pr\{\|\widehat B-B\|>\eta_N(\delta)\}\le\delta
$$

的 matrix concentration bound，再代入 Weyl gap rules。

**所需工具。** matrix Bernstein/sub-Gaussian concentration；finite-difference bias bound；sequential change detection。

**当前卡点。** 当前 $B$ 是 deterministic code Jacobian，尚无 Jacobian-estimation noise model。

## 8.8 P7：moving-grid 情形的 $D^2B=D_X^3F$ continuum bound

**未解决对象。** 对 $G_D(X)$、moving map grid、self-cell quadrature与天线 phase center同时依赖 pose时，给完整 $D_X^3F$ bound。

**最小下一步。** 先处理无 self-interaction 的 external-sensor Born model，写出三阶 Hankel radial tensor；再加入 fixed-domain full-wave resolvent；最后单独证明 singular volume potential在 Sobolev spaces中的三阶 shape derivative。

**所需工具。** shape calculus；singular integral mapping；Hankel recurrence；Faà di Bruno/resolvent derivatives；Sobolev multiplication。

**当前卡点。** 目标域内 Green kernel在 diagonal有 logarithmic singularity，不能以简单 pointwise $d_0$ bound处理；moving cells会改变积分域而不仅是 kernel argument。

## 8.9 P7：verified full nonlinear Lipschitz certificate

**未解决对象。** 在一个非平凡 trajectory ball上计算可用而不过度松弛的

$$
\sup_{X\in\mathcal B}\|DK_{\rm eff}(X)\|.
$$

**最小下一步。** 选择很小的 finite-dimensional task/pose model，对 Hankel、matrix inverse和 sparse solve使用 interval arithmetic；以 branch-and-bound subdivision控制 wrapping effect。

**所需工具。** interval linear algebra；verified SVD/eigenvalue bounds；Krawczyk operator；Taylor models。

**当前卡点。** 解析 norm chain在 $m_0^{-1}$、$k$和 domain size上极松；随机方向扫描不提供 supremum certificate。

## 8.10 P8：rank-event 跨越时的 nonsmooth robust trajectory solver

**未解决对象。** 构造对

$$
X\mapsto\lambda_r(A(X)^*P_{B(X)^\perp}A(X))
$$

跨 rank strata仍收敛的算法。

**最小下一步。** 比较三种 formulations：

1. positive-prior smoothing $J_X=\alpha I$ 的 continuation；
2. bounded pose nuisance的 min--max residual；
3. spectral bundle method作用于 cluster sum或 Ky Fan objective。

证明 accumulation point满足 Clarke stationarity或相应 saddle-point condition。

**所需工具。** variational analysis；Clarke subdifferential；spectral bundle/trust-region；epi-convergence as $\alpha\downarrow0$。

**当前卡点。** no-prior projector在 rank change处 norm jump，普通 smooth continuation未必逼近同一 branch。

## 8.11 P9：TSOM $G_D$ split 与 pose projection 的统一误差

**未解决对象。** 对 Chen TSOM 的实际 $V_S^-\cap V_D^+$ approximation，加上 map pullback和 pose elimination后给 end-to-end error。

**最小下一步。** 从原 TSOM明确定义出发，分别界定：

$$
\varepsilon_S=\|(I-P_S)J\|,
\quad
\varepsilon_D=\|(I-P_D)P_S^-J\|,
\quad
\varepsilon_\chi=\|A-A_r\|,
\quad
\varepsilon_P=\|P_B-\widehat P_B\|.
$$

再用 triangle/perturbation bounds合成 $K_{\rm eff}$ error。

**所需工具。** CS decomposition；noncommuting projector estimates；Wedin/Davis--Kahan；randomized SVD residual certificates。

**当前卡点。** 不同 TSOM文献对 $V_D^+$、交集 basis与 Fourier surrogate有具体实现差异，不能先抽象成一个未经核对的 projector公式。

## 8.12 P10：从 local basin 到可观测 noise threshold

**未解决对象。** 给出噪声、波数和 initial pose误差的定量条件，使 estimator以高概率落入正确 phase branch。

**最小下一步。** 从简化 1D/2D Born phase model出发，证明：

- local strong convexity radius $r<2\mu/L$；
- initialization error小于 phase half-period与该 radius；
- noise score不越过 basin boundary的 concentration bound。

随后再加入 full-wave resolvent remainder。

**所需工具。** Newton--Kantorovich theorem；small-noise concentration；periodic likelihood alias analysis；Morse/critical-point theory。

**当前卡点。** 多频、多路径与 full-wave multiple scattering使远处 stationary points的结构高度 scene-dependent；Fisher只提供名义点局部 curvature。

## 8.13 文献：投稿前系统查重

**未解决对象。** full-text 层面确认是否已有相同的“whitened/realified volumetric map tangent vs physical pose tangent principal-angle retention” formulation。

**最小下一步。** 对以下 query clusters做数据库级系统检索并保留 PRISMA-style日志：

- inverse scattering + unknown sensor/source location + Fisher/Schur;
- microwave tomography + antenna position error + CRB;
- autofocus + principal angles/canonical correlations;
- blind calibration + PDE inverse problem + nuisance tangent;
- wave-based SLAM + material map + observability;
- SOM/contrast source + localization/calibration.

对 Karthik--Ghosh、Deutschmann et al.、radar autofocus与 microwave calibration全文逐式核对。

**所需工具。** IEEE Xplore、Scopus/Web of Science、MathSciNet、Google Scholar/Semantic Scholar及出版社全文访问。

**当前卡点。** 现有检索受 429 与部分登录墙限制；“未检得”不能升级成 global novelty。

---

# 9. 反过度主张审计

下面逐类检查本报告和建议 manuscript language。每一行给出禁止表述、审计结果与已经采用的修正版。

| 风险类型 | 禁止/危险表述 | 审计结论 | 已采用的修正版 |
|---|---|---|---|
| 有限维当连续 | “离散 rank identity证明 continuum rank identity” | 已删除 | 有限维 rank公式仅作 matrix theorem；continuum改用 closure、distance、forms与 task-window convergence |
| 有限网格角度当 certificate | “每个网格 $\theta_{\min}>0$，故稳定” | 被反例 4.5 否定 | 必须给 uniform lower bound、closed-sum theorem或固定 task frame bound |
| machine rank当物理 rank | “SVD tolerance检测到 rank change，所以物理 discontinuity” | 已删除 | 区分 exact rank、stable/effective rank与 tolerance classification；需 perturbation margin |
| local 当 global | “$K_{\rm eff}\succ0$ 保证唯一重建和收敛” | 被反例 4.82 否定 | 只称 local tangent identifiability/information；另列 basin、alias与 cycle-skipping条件 |
| Fisher 当 estimator MSE | “covariance = $K_{\rm eff}^{-1}$” | 仅特定生成模型成立 | fixed-nuisance protocol用 sandwich；Bayesian random nuisance或 auxiliary-measurement protocol才用相应 inverse information |
| sampled 当 certified | “500/2000 random directions无 violation，故 nonlinear robust” | 已降级 | 只称 sampled evidence；certificate需要 analytic/interval supremum bound |
| 内部控制当外部验证 | “mixed $(A_{\rm Born},B_{\rm full})$ 验证 Born-SLAM” | 已禁止 | mixed pair仅是 isolation control；complete Born pair必须使用 $B_{\rm Born}$ |
| current-space 当 map-space | “Chen $V_S^+$ 就是 map observable eigenmodes” | 已禁止 | 只通过 $A=G_ST_\chi$、range/pullback与 reduction error联系 |
| $G_D$ 与 pose space混合 | “SOM 的 $V_D$ 是位姿 data space” | 已禁止 | $V_D$ 属 current/domain operator；pose tangent是 data-space $\operatorname{Ran}B$ |
| complex pose projector | “对 complex $B$ 直接用 $BB^\dagger$” | 已禁止 | 先 whitening，再按 real physical parameters realify |
| $j$ 与 $J_X$ 混淆 | “$J$ 同时表示 current 和 pose information” | 已修正 | current用 $j$ 或 $J_{\rm cur}$；prior information固定为 $J_X$ |
| inner/outer uncertainty重复计数 | “Schur消元的 $\delta X$ 再作为同一随机 $\Delta X$ robust perturbation” | 已禁止 | inner nuisance是 estimator-level local variable；outer execution mismatch是 design robustness变量，需声明生成关系 |
| exact intersection与 stable transversality混淆 | “无 exact ambiguity即有稳定 inverse” | 被反例否定 | 分别报告 $M\cap N$ 与 $\gamma/\theta_F$ |
| 闭包语义漂移 | “$K_{\rm SLAM}u=0$ 必存在 $h$ 使 $Au=Bh$” | 无限维不成立 | 只说 $Au\in\overline{\operatorname{Ran}B}$；exact attainment需 closed range |
| ordinary eigenvalue ratio | “$\lambda_i(K_{\rm eff})/\lambda_i(K_{\rm IS})=\rho_i$” | 一般不成立 | 使用 generalized eigenproblem或 polar retention operator |
| finite prior仍称 ordinary angle | “$\rho_i=\sin^2\theta_i$ 对任意 $J_X$” | 一般不成立 | 原数据空间称 weighted shrinkage；只在 augmented space有 angle interpretation |
| low-rank夸大 | “pose只影响 $q$ 个 ordinary eigenvectors” | 不成立 | loss rank/generalized defect至多 $q$；ordinary eigenvectors可全部旋转 |
| rank-event smoothness | “projector在 rank change处仍 Lipschitz/Hölder” | 被定理 4.48 否定 | fixed-rank才 smooth；不同有限 rank projections的 norm distance为 1 |
| 频率单调性 | “更多/更宽频率必提高每个 retention” | 被反例否定 | absolute $K_{\rm eff}$ Loewner单调；normalized spectrum scene/metric dependent |
| submodularity | “shared-pose logdet一定 submodular” | 被两-block反例否定 | 仅 additive independent-elimination lower-bound surrogate有标准 submodularity |
| trajectory universal ordering | “圆周优于圆弧优于直线” | 离散 falsifier已否定且无定理支持 | 使用 Pareto design；排序必须注明 scene、metric、standoff、energy与 bandwidth |
| geometry-only角界 | “aperture/MIMO/standoff给 continuum统一 $\gamma>0$” | 一般不成立 | 必须固定 map class/task basis并给 frame/coherence bound |
| near resonance | “resonance附近信息大，所以更好” | 单边表述错误 | 同时报告 information gain、resolvent derivative、model mismatch与 basin shrinkage |
| gauge维数固定 | “2D 总有恰好3个 map gauge modes” | 不成立 | 最多3个 rigid generators；stabilizer、anchor、background与 basis可降低/破坏 |
| 离散 gauge residual物理化 | “非零 residual证明 continuum gauge被破坏” | 不成立 | 先分解 basis representation、quadrature、window与真实 model-breaking项 |
| novelty过度 | “首次/no prior work” | 已禁止 | 使用 retrieval-bounded claim ledger；只说“we formulate/analyze/test” |
| 联合定位 novelty | “首次联合 inverse scattering 与 transmitter localization” | 被 Karthik--Ghosh直接否定 | 贡献限定为 tangent-information geometry及其条件性 consequences |
| 文献未检得当不存在 | “search未找到，因此没有工作” | 已禁止 | 写“本次检索未发现；full-text coverage不完整” |

## 9.1 最终语义一致性检查

1. 全文 $A$ 始终是 **map tangent to data**；$B$ 始终是 **real physical-pose tangent to data**。
2. $G_S$ 始终是 **current to data**；其 right singular vectors不被称作 map modes。
3. $J_X$ 始终是 pose-prior information；contrast source写作 $j$。
4. $K_{\rm SLAM}$ 的 continuum版本始终使用
   $P_{\overline{\operatorname{Ran}B}^{\perp}}$。
5. 无 prior projector result与 finite-prior weighted result分开。
6. inner nuisance $h$ 与 outer trajectory execution $\Delta X$ 分开。
7. full-wave、complete Born与 mixed-tangent control分开。
8. exact algebra、conditional continuum theorem、numerical observation、counterexample与 open problem分开标注。
9. 所有 covariance陈述都附带 generating protocol。
10. 所有 robust陈述都附带 compact set、rank/gap/resolvent margin或明确标为 sampled evidence。

## 9.2 最终可主张的论文中心句

经过审计后，最强且不过界的中心句是：

> 在白化并按实物理参数 realify 的局部 full-wave inverse-scattering SLAM 模型中，位姿 nuisance 对地图 tangent information 的作用可由 map-data range 与 closed pose-data range 的相对几何描述。无先验时，polar retention operator给出 generalized squared-sine geometry；半正定先验时，正确对象是 augmented-space shorted form。精确交集只描述 exact/asymptotic compensation，而 positive Friedrichs-angle constant才描述 stable transversality。该结构允许在固定 task window、固定 rank和正谱隙下推导 gauge、multi-frequency、rank-event、SOM reduction与 robust-trajectory consequences；它不蕴含 global nonlinear recovery。

---

# 10. 最终状态与使用指南

## 10.1 P1--P10 最终状态

| 包 | 最终标签 | 最适合放置位置 |
|---|---|---|
| P1 | **solved + counterexample + conditional discretization** | 主定理 + appendix |
| P2 | **solved under homogeneous scalar model; conditional extensions** | 主定理 + limitations |
| P3 | **impossibility result + finite-dimensional conditional bounds** | theory/negative result |
| P4 | **finite-dimensional conditional theorem + counterexamples** | main/appendix；具体 witness留 future work |
| P5 | **solved via augmented shorted form** | 主定理 |
| P6 | **solved for matrices/fixed-rank paths** | appendix + algorithm box |
| P7 | **conditional bound chain; original $D^2B$ target narrowed** | appendix + limitations |
| P8 | **local theorem solved; global nonsmooth design open** | methods + limitations |
| P9 | **algebraic bridge solved; quantitative continuum bridge conditional** | main semantic theorem + algorithm |
| P10 | **local statistical semantics solved; global recovery open** | statistical scope + limitations |

## 10.2 推荐的首篇理论稿主线

为保持剃刀原则，首篇稿件建议只保留四条主线：

1. **Hilbert retention + stable transversality：** polar operator、closure、Friedrichs angle、compact counterexample；
2. **physical gauge + semidefinite prior：** $SE(d)$ tangent与 augmented-space shorting；
3. **rank-stratified consequences：** fixed-rank derivative、rank event、Born $s=0$ singularity；
4. **合法 SOM bridge：** $A=G_ST_\chi$、range theorem、two-stage reduction与 noncommutation warning。

P3 的 line/circle显式几何、P4 的具体 array witness、P7 的完整 interval certificate、P8 的大型 online trajectory optimization和 P10 的 global basin theorem宜作为后续工作，而不是在首篇中同时展开。

## 10.3 本报告的证据层级

- **本报告内完整证明：** 抽象 Hilbert/matrix theorem与反例；
- **附加假设下证明：** homogeneous scalar gauge、finite analytic genericity、task-window discretization、resolvent bounds、local robust expansion；
- **来源文件中的 finite-dimensional evidence：** 当前 round-3 discrete harness及其 preserved failures；
- **文献定位：** retrieval-bounded claim ledger；
- **仍开放：** 具体连续 Helmholtz discretization、array witness、full nonlinear uniform certificate与 global estimation。

## 10.4 引用与复用注意

把本报告内容迁入论文时：

- theorem可直接复制，但需统一论文中的 real/complex inner-product convention；
- 任何 “compact task space” 必须在实验中具体定义；
- 若 $J_X$ 是 graph Laplacian，必须说明 anchor或 quotient；
- 若使用 Moore--Penrose inverse，必须保留 range/closed-range condition；
- 若使用 finite-dimensional principal angles，必须注明它们属于 $\operatorname{Ran}A$ 与 $\operatorname{Ran}B$ 的 data-space geometry；
- 引用 Chen SOM 时必须写明 current-space对象；
- novelty段必须保留 retrieval-bounded措辞。

