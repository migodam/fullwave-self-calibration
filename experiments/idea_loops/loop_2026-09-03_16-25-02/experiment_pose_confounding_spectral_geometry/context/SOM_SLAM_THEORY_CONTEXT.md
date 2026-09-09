# 逆散射 SOM / SLAM 理论讨论总 Context

> 用途：长期本地参考；汇总“电磁逆散射”项目中 7 个既有讨论的定义、符号、公式、结论、质疑、冲突与开放问题。本文不是一次新的文献综述，也不把讨论中的猜想升级为已证事实。
>
> 整理日期：2026-09-03

## 0. 阅读说明与证据状态

### 0.1 状态标签

- **[已推导/已接受]**：在原讨论中已经给出代数推导，或作为后续讨论共同采用的工作结论。它表示“讨论内部已成立”，不等同于同行评审或外部文献已经验证。
- **[需补条件]**：核心式子可成立，但必须补充有限维、闭值域、可逆性、谱隙、gauge 固定、噪声白化等条件。
- **[直觉/假设]**：用于组织研究的解释、模型选择或尚未证明的物理直觉。
- **[待验证/开放问题]**：原讨论明确提出，但没有完成严格证明、文献查重或数值验证。
- **[冲突/修正]**：讨论中出现过不同说法，或后续对早期表述作了收窄。

文中的短标签 **[已推导]**、**[已接受]** 均归入“已推导/已接受”，**[待验证]** 归入“待验证/开放问题”；“已推导/需条件”“直觉/需条件”等复合标签表示两个状态同时成立。

### 0.2 七个来源讨论

| 代号 | 用户所列主题 | 项目内实际标题 | Conversation ID |
|---|---|---|---|
| S1 | 全波逆散射 SOM 理论研究 | 全波逆散射SLAM理论研究 | `6a97d0b1-d8cc-83ec-957a-7e472c4854e9` |
| S2 | 总结 SOM 逆散射 SOM 推导 | 总结SOM逆散射SLAM | `6a97b0f1-a934-83ec-9381-2b72ae79556f` |
| S3 | 位置不确定下信息损失 | 推导位姿不确定性下信息损失 | `6a8c0ab7-4194-83ec-8ee7-92fc18cc1d9e` |
| S4 | Regularization 低秩直觉 | Regularization低秩直觉 | `6a9682cb-9588-83ec-981f-b2c335b4c0f6` |
| S5 | 迁移至 SOM 推导 | 迁移至SLAM | `6a969bfe-9eb4-83ec-8427-ac39a935b67d` |
| S6 | 谱可观测性 | 推导谱可观测性 | `6a8c498e-c6cc-83ec-989b-5c2b7cda463b` |
| S7 | 逆散射与 SOM 技术 | 逆散射与SLAM技术 | `6a5726e6-0338-83ec-aa65-a73247c4cc02` |

标题中的 “SOM/SLAM” 混用来自用户列名与项目内实际标题的差异；本文不把两者视为同义词：SOM 是逆散射中的子空间优化方法，SLAM 是位姿与地图联合估计问题。

### 0.3 本文最重要的防混淆规则

原讨论中的“位姿不确定性”至少对应三种不同数学模型，后文始终分开：

1. **未知但要联合估计的确定性 nuisance parameter**：消元后得到投影或 Schur complement，即 $K_{\mathrm{SLAM}}$、$K_{\mathrm{eff}}$。
2. **有概率分布的随机位姿误差**：边缘化后进入有效噪声协方差，或等价地通过 Gaussian prior information 进入 Schur complement。
3. **有界但不设概率的模型/轨迹扰动**：研究 $K(X+\Delta X)$ 的最坏谱与一阶敏感度。

三者回答的问题不同，不能把它们的公式直接互换。

---

## 1. 一页式总览

### 1.1 从 Chen SOM 到逆散射 SLAM 的主线

已知位姿的 contrast-source 逆散射先写成

$$
E^{\mathrm{sca}}=G_SJ,
\qquad
J=D_\chi(E^{\mathrm{inc}}+G_DJ).
$$

Chen SOM 的核心是对已知的 current-to-data operator $G_S$ 做谱分解，把诱导电流 $J$ 分成数据直接稳定确定的部分与需要 state equation、优化或先验补足的部分。

当收发平台位姿 $X=(x_1,\ldots,x_T)$ 也未知时，$G_S$、入射场，乃至坐标表达下的内部算子都可能依赖 $X$。在名义状态附近，全波联合模型的一阶形式是

$$
\boxed{
\delta y=A\,\delta\chi+B\,\delta X+n,
}
$$

其中 $A=D_\chi F$ 是地图 Jacobian，$B=D_XF$ 是位姿 Jacobian。由此得到本文的几何核心：

$$
\boxed{
\operatorname{Range}(A)\cap\operatorname{Range}(B)
}
$$

表示“地图变化”和“位姿变化”在数据空间中能够造成相同变化的部分，即局部地图—位姿混淆。

噪声白化后，已知位姿的地图信息为

$$
K_{\mathrm{IS}}=A^*A.
$$

位姿完全自由、作为确定性 nuisance parameter 消去后：

$$
\boxed{
K_{\mathrm{SLAM}}
=A^*P_{B^\perp}A,
\qquad
P_{B^\perp}=I-BB^\dagger.
}
$$

加入运动/IMU/里程计先验信息 $J_X\succeq0$ 后：

$$
\boxed{
K_{\mathrm{eff}}
=A^*A-A^*B(B^*B+J_X)^\dagger B^*A.
}
$$

所以轨迹具有双重作用：

$$
\boxed{
\text{轨迹改变采样几何并制造信息；轨迹未知又通过地图—位姿耦合消耗信息。}
}
$$

### 1.2 讨论最终收敛到的研究对象

讨论从“dominant/complementary subspace 是否旋转”逐步收敛到更稳健的对象：

- 完整正半定信息算子 $K(X)$，而不是单个基向量；
- 谱投影或连续谱函数 $f(K)$，而不是未经说明的 hard threshold；
- $K_{\mathrm{IS}}$ 与 $K_{\mathrm{SLAM}}$ 的广义特征值/信息保留谱；
- $K_{\mathrm{eff}}(X)$ 的轨迹导数、最弱任务模式及最坏位姿扰动；
- 由实际采样测度产生的 frame/normal operator；
- 在 Born 条件下的 Fourier/phase-space visibility，以及 full-wave 下的局部 Jacobian 与多重散射修正。

建议性的最终极大极小问题是

$$
\boxed{
\max_{X\in\mathcal X}
\min_{\|\Delta X\|\le\varepsilon}
\lambda_r\!\left(K_{\mathrm{eff}}(X+\Delta X)\right),
}
$$

其含义是：在运动约束和最坏位姿误差下，使至少 $r$ 个地图自由度仍尽可能可观测。

---

## 2. 问题边界：逆散射、成像、Radio/Radar SLAM 与全波 SLAM

### 2.1 逆散射并不自动等于 SLAM，也不自动等于完整 3D 可观测

**[已接受]** 原讨论形成的边界是：

- **逆散射**通常假设收发位置、照射和背景已知，反演介电对比度、导电率、折射率、散射势或等效源。
- **SLAM**要求传感器自身状态未知，利用同一批观测联合估计轨迹与地图，并处理 gauge、回环、数据关联和先验。
- **逆散射 SLAM / wave-based SLAM**以 Maxwell/Helmholtz 全波前向算子作为观测模型，地图是物理参数场，位姿与地图共同未知。
- 输出三维体素或 3D 图像不意味着每个体素都由数据独立支持；稳定自由度由显著奇异值、有效信息谱和先验共同决定。

统一观测可写成

$$
y_t=F(\chi,x_t,c;d_t)+n_t,
$$

其中 $c$ 可含外参、时钟、相位中心、通道增益、背景波速/介电常数等标定参数，$d_t$ 可含频率、波形、极化、波束和主动轨迹。

### 2.2 与视觉 Bundle Adjustment 的结构对应

视觉 BA 使用

$$
u_{tj}=\pi(T_t^{-1}P_j),
$$

联合优化位姿 $T_t$ 与路标 $P_j$。逆散射 SLAM 使用

$$
y_t=F(T_t,\chi),
$$

联合优化位姿与连续/离散物性场。两者都有 block Hessian、Schur complement、gauge 与轨迹—地图耦合；区别在于全波算子通常具有更全局的耦合、多次散射、相位敏感性和高维地图块。

### 2.3 与 Radio SLAM、multipath SLAM、Radar SLAM 的边界

- Multipath SLAM 常把镜面反射表示为 virtual anchor、平面或离散散射点；它是 SLAM，但通常不是连续介质的 full-wave inverse scattering。
- Radar SLAM 常在检测/点云层做配准；direct RF/Radio SLAM 更接近直接用 IQ/CSI/波形似然。
- 高阶或 double-bounce 路径既可能增加可观测性，也会带来路径关联与全局多峰问题。
- 全波逆散射 SLAM 的目标更接近“轨迹 + 几何 + 材料/传播参数”的物理地图，而不只是几何点云。

### 2.4 五种不能混叫的“极限”

1. **唯一性/可辨识性**：无噪声下 $F(\theta_1)=F(\theta_2)$ 是否只差已知 gauge。
2. **局部可观测性**：Jacobian 是否存在非 gauge 的零方向。
3. **稳定性**：小奇异值是否把噪声或模型误差无限放大。
4. **统计精度**：FIM/CRB/PCRB/Bayesian posterior 能达到多小误差。
5. **全局阈值与多峰失效**：cycle skipping、相位周跳、错误回环/路径关联等何时使局部 CRB 失效。

**[已接受]** “理论唯一”“局部 Fisher 信息大”“算法能收敛”和“输出视觉上合理”是四个不同命题。

### 2.5 wave-based SLAM 的统一 MAP 形式与工程层次

把 motion、calibration 与全波 data likelihood 同时写入，可用

$$
\boxed{
\begin{aligned}
\min_{\chi,X,c}\quad&
\sum_t
\|y_t-F(\chi,x_t,c;d_t)\|_{\Sigma_t^{-1}}^2\\
&+
\sum_t
\|x_{t+1}-f(x_t,u_t)\|_{Q_t^{-1}}^2\\
&+
\mathcal R_\chi(\chi)
+
\mathcal R_c(c).
\end{aligned}
}
$$

其中第一项是 wave/data residual，第二项是 motion/odometry/IMU prior，后两项分别是 map 与 calibration prior/regularization。局部线性化和 Schur complement 是该 nonlinear MAP 的一轮 Gauss–Newton/Fisher 几何，不等于完整全局求解。

原讨论把一般 SLAM pipeline 拆为：

1. sensing 与时间同步；
2. front end（特征、相关、direct alignment 或 wave residual）；
3. data association / loop closure；
4. back end（pose graph、factor graph、BA 或 joint field inversion）；
5. map representation 与 uncertainty；
6. active sensing / trajectory update。

wave-based 版本可使用 occupancy/grid、point/surface landmarks、TSDF/SDF、mesh、voxel/material contrast field 或 neural implicit map；不同表示改变 $A$ 的参数空间、regularizer 与 gauge 离散误差，不能只换名词而沿用同一 rank 结论。

两个反复出现的工程难点是：

$$
\Delta\phi\approx\frac{2\pi}{\lambda}\Delta R,
$$

即 carrier phase 对微小 range/pose error 极敏感；以及 multipath 既携带额外 geometry/material information，又引入 path association、多峰 posterior 与 model mismatch。本文的 local spectrum 只刻画名义解附近的一阶部分。

---

## 3. 统一符号表

| 符号 | 含义 | 所在空间/备注 |
|---|---|---|
| $D\subset\mathbb R^d$ | 成像/散射区域 | 讨论通常先取 2D scalar Helmholtz |
| $\chi$, $D_\chi=\operatorname{diag}(\chi)$ | 介质 contrast 及其离散乘法算子 | 地图/参数空间 |
| $E_t^{\mathrm{inc}},E_t^{\mathrm{tot}},E_t^{\mathrm{sca}}$ | 入射、总场、散射场 | 可随 Tx、频率、位姿变化 |
| $J_t$ 或 $j_t$ | contrast source / induced current | 目标区域内部 current space |
| $G_D$ | 区域内部传播算子 | current → domain field |
| $G_S(x_t)$ | 目标区域到接收阵列的观测算子 | current → data；依赖 Rx/几何 |
| $M_t=I-D_\chi G_D$ | full-wave state operator | 可逆性需说明；接近共振时敏感 |
| $F(\chi,X)$ | 消去 $J$ 后的完整前向映射 | 地图与轨迹 → data |
| $X=(x_1,\ldots,x_T)$ | 轨迹/全部位姿参数 | 每个位姿可为 2D/3D Lie algebra 坐标 |
| $A=D_\chi F$ | 地图 Jacobian | whitened 后默认吸收 $\Sigma_n^{-1/2}$ |
| $B=D_XF$ | 位姿 Jacobian | 可按时刻、Tx/Rx、外参分块 |
| $K_{\mathrm{IS}}$ | 位姿已知的地图信息 | $A^*A$ |
| $K_{\mathrm{SLAM}}$ | 完全自由位姿消元后的地图信息 | $A^*P_{B^\perp}A$ |
| $J_X$ | 位姿/运动先验信息矩阵 | 非地图 Jacobian；不要与 current $J$ 混淆 |
| $K_{\mathrm{eff}}$ | 有位姿先验时的有效地图信息 | block Schur complement |
| $L_X$ | pose-induced information loss | $K_{\mathrm{IS}}-K_{\mathrm{eff}}$ |
| $G^\dagger$ | Moore–Penrose 伪逆 | 奇异/gauge 情况需限制支撑空间 |
| $P_{\mathcal U}$ | 到子空间 $\mathcal U$ 的正交投影 | $P_{B^\perp}=I-P_{\operatorname{Range}(B)}$ |
| $C_0,C_X,\Sigma_n$ | 地图先验、位姿先验/协方差、噪声协方差 | 信息与协方差不能混写 |
| $\rho_i$ | 信息保留谱/广义特征值 | $0\le\rho_i\le1$（在适当支撑上） |
| $\mu_T$ | 轨迹诱导的经验采样测度 | continuous frame 表述 |

本文默认复数数据内积使用共轭转置 $^{*}$ 或 $^{H}$。原讨论中 $T,H,*$ 有混用；实数模型才可直接用转置。

---

## 4. Full-wave contrast-source 基础

### 4.1 Helmholtz 与 Lippmann–Schwinger

背景波数 $k_b$，介质局部波数 $k(r)$，定义

$$
\chi(r)=\frac{k^2(r)-k_b^2}{k_b^2}.
$$

二维标量模型中

$$
(\nabla^2+k_b^2)E^{\mathrm{tot}}
=-k_b^2\chi E^{\mathrm{tot}},
$$

从背景 Green 函数得到

$$
E^{\mathrm{tot}}(r)
=E^{\mathrm{inc}}(r)
+k_b^2\int_DG_b(r,r')\chi(r')E^{\mathrm{tot}}(r')\,dr'.
$$

定义（常数因子可吸收进算子）

$$
J(r)=\chi(r)E^{\mathrm{tot}}(r),
$$

便得到两条核心方程：

$$
\boxed{E^{\mathrm{sca}}=G_SJ}
\qquad\text{(data equation)},
$$

$$
\boxed{J=D_\chi(E^{\mathrm{inc}}+G_DJ)}
\qquad\text{(state equation)}.
$$

前者对 $J$ 线性，说明外部仪器能看见哪些 current modes；后者使 $J$ 与材料场自洽，承载 full-wave 非线性。

### 4.2 消去 current 后的非线性

离散形式

$$
(I-D_\chi G_D)J=D_\chi E^{\mathrm{inc}},
$$

故

$$
J=M^{-1}D_\chi E^{\mathrm{inc}},
\qquad
M=I-D_\chi G_D,
$$

以及

$$
F(\chi)=G_SM^{-1}D_\chi E^{\mathrm{inc}}.
$$

**[已接受]** 逆散射非线性不是因为给定介质后的 Maxwell/Helmholtz 场方程非线性，而是未知介质改变总场、总场又决定等效源。

### 4.3 Born 近似及适用直觉

若内部再散射算子足够弱，可形式展开

$$
(I-D_\chi G_D)^{-1}
=I+D_\chi G_D+(D_\chi G_D)^2+\cdots.
$$

截断到首项得到

$$
J\approx D_\chi E^{\mathrm{inc}},
\qquad
E^{\mathrm{sca}}\approx G_SD_\chi E^{\mathrm{inc}}.
$$

**[冲突/修正]** Born 条件不能只说“介电常数低”；更相关的是“内部多次散射作用小”，例如某个合适范数/谱半径意义下 $\|D_\chi G_D\|<1$，并共同依赖 contrast、尺度、波长与几何。

### 4.4 两层信息瓶颈

完整信息链为

$$
\chi
\longrightarrow J
\overset{G_S}{\longrightarrow}E^{\mathrm{sca}}
\longrightarrow y.
$$

- **外部观测瓶颈**：$G_S$ 的谱决定哪些内部 current modes 能到达测量空间。
- **物理参数辨识瓶颈**：不同 $\chi$ 是否产生足够可区分的 $J$，受多重散射、照射多样性、内部共振与非线性耦合影响。

因此即使 $J$ 的一部分可被稳定恢复，也不等于 $\chi$ 已被稳定恢复；反之，正则化或学习生成的细节也不自动是测量信息。

---

## 5. Chen 子空间谱系：SOM、MSR、MUSIC 与讨论中的术语边界

### 5.1 原始 SOM 的最小数学骨架

对某次 illumination $p$，data equation 写成

$$
s_p=G_Sw_p,
$$

其中 $w_p$ 是该 illumination 的 contrast source，不是材料 contrast $\chi$ 本身。令

$$
G_S=U\Sigma V^*,
\qquad
\sigma_1\ge\sigma_2\ge\cdots.
$$

选择 dominant/data-resolvable 部分后，

$$
w_p=w_p^{+}+w_p^{-},
$$

$$
w_p^{+}
=\sum_{i\in\mathcal I_+}
\frac{u_i^*s_p}{\sigma_i}v_i,
\qquad
w_p^{-}\in\operatorname{span}\{v_i:i\in\mathcal I_-\}.
$$

原讨论也写作

$$
V_S^{+}=\operatorname{span}\{v_i:i\le L\},
\qquad
V_S^{-}=\operatorname{span}\{v_i:i>L\}.
$$

这里的下标 $S$ 指 sensing/data operator $G_S$；“$+/-$”表示 dominant 与 complementary/ambiguous，而非正负特征值，因为 $G_S^*G_S\succeq0$。

**[已接受]** SOM 的哲学不是把小奇异方向宣告为“没有物理信息”，而是：

$$
\boxed{
\text{数据已经稳定确定的 current component 不再让非线性优化重猜；其余部分由 state equation、优化和先验补足。}
}
$$

“10000 个 pixels → 约 30 个模式”若出现，只表示给定几何、频率、噪声与阈值下约 30 个稳定 data-resolvable current directions；它不表示真实环境本身只有 30 维。

### 5.2 $G_S$ 的 current spectrum 不等于地图 Jacobian 的 observability spectrum

这是后续语义审计中最重要的修正之一：

$$
\boxed{A=\frac{\partial y}{\partial\chi}\neq G_S\quad\text{(一般 full-wave 情形)}.}
$$

在本模型下

$$
A=G_SM^{-1}\operatorname{diag}(E^{\mathrm{tot}}),
$$

Born 下才退化为

$$
A\approx G_S\operatorname{diag}(E^{\mathrm{inc}}).
$$

所以至少要区分两套谱：

1. $G_S$ 的谱：外部阵列能稳定看见哪些 **contrast-source/current modes**；
2. $A=J_\chi$ 或 $K_{\mathrm{eff}}$ 的谱：地图/材料变化在当前状态与位姿耦合下有多可观测。

### 5.3 MSR 与 MUSIC：同样使用 SVD，但子空间所在位置不同

多发多收时，multistatic response (MSR) 矩阵可记为

$$
Y=[y_{r,p}]_{r,p},
$$

其列对应 illuminations/Tx，行对应 receivers/Rx。对小散射体或适当线性近似，$Y$ 可近似由少数 steering vectors 的外积组成，因而近低秩。MUSIC 对 **数据/MSR 矩阵** 的 signal/noise subspace 做投影测试：候选位置的 steering vector 若接近 signal subspace，则其到 noise subspace 的投影较小。

SOM 则主要对 **已知的 current-to-data operator $G_S$** 做谱分解，并用实测数据确定 current 的 dominant component。二者都使用 SVD/投影，但不能把以下对象不加说明地等同：

- MSR 的左/右 singular subspaces（Tx/Rx 数据空间）；
- $G_S$ 的左 singular subspace（data space）；
- $G_S$ 的右 singular subspace $V_S^\pm$（current/domain space）；
- $A=J_\chi$ 的 right singular modes（map space）；
- $B=J_X$ 的 column space（pose-induced data space）。

### 5.4 $V_D^\pm$ 与 Two-fold SOM：保留原讨论的定义边界

原讨论把 $G_D$ 固定为内部/domain propagation operator，并用 $V_D^+\oplus V_D^-$ 表示对其谱或相关 domain operator 作第二次 dominant/complementary 分解；Two-fold SOM 的直觉是同时利用 **外部 data equation 的可见性结构** 与 **内部 state/domain equation 的结构**，进一步减少需要优化的模糊自由度。

但讨论也留下了必须保留的警告：不同文献版本可能对 $V_D^\pm$ 使用 $G_D$、$G_D^*G_D$、经过材料/场加权的组合算子，或两套子空间的交集/组合；因此在没有指定 Chen 原文版本、离散约定与式号前，不能把一个唯一的 Two-fold 公式写成已验证定理。本文只保留安全的术语映射：

$$
\begin{array}{ccl}
V_S^+&:&\text{由外部测量算子稳定支持的 current directions},\\
V_S^-&:&\text{外部数据弱支持的 complementary current directions},\\
V_D^+&:&\text{内部/domain operator 的 dominant directions},\\
V_D^-&:&\text{内部/domain operator 的 complementary directions}.
\end{array}
$$

**[待验证]** Two-fold SOM 的精确投影次序、优化变量和 $V_D^\pm$ 定义必须回到所指具体论文/章节核对；七条讨论中存在概念性映射，但没有形成一个跨版本无歧义的统一公式。

### 5.5 hard split 与 soft subspace

hard spectral projector 可写为

$$
P_\tau(K)=\mathbf 1_{[\tau,\infty)}(K).
$$

若特征值穿过 $\tau$，物理上连续的变化会被人为变成 rank 跳变。讨论因此提出连续替代：

$$
W_\alpha=f_\alpha(K),
\qquad
f_\alpha(\lambda)=\frac{\lambda}{\lambda+\alpha},
$$

或在 prior-whitened Hessian

$$
H=C_0^{1/2}J^*\Sigma_n^{-1}JC_0^{1/2}
$$

上定义

$$
o_i=\frac{\lambda_i}{1+\lambda_i},
\qquad
d_{\mathrm{eff}}=\sum_i\frac{\lambda_i}{1+\lambda_i},
$$

$$
I(m;y)=\frac12\sum_i\log(1+\lambda_i)
=\frac12\log\det(I+H).
$$

解释分别是 continuous observability/data-confidence、有效信息维数和线性 Gaussian information gain。它们把 “physics/data” 与 “prior/learning” 从二元切分改成逐 mode 的连续权重。

---

## 6. 从 full-wave 方程严格线性化到 $A,B$

### 6.1 地图 Jacobian $A$

固定 geometry，令

$$
j=D_\chi(e^{\mathrm{inc}}+G_Dj),
\qquad
E^{\mathrm{tot}}=e^{\mathrm{inc}}+G_Dj.
$$

对 $\chi$ 作一阶扰动：

$$
\delta j
=D_{\delta\chi}E^{\mathrm{tot}}
+D_\chi G_D\delta j.
$$

由于

$$
D_{\delta\chi}E^{\mathrm{tot}}
=\operatorname{diag}(E^{\mathrm{tot}})\delta\chi,
$$

得到

$$
\boxed{
\delta j
=M^{-1}\operatorname{diag}(E^{\mathrm{tot}})\delta\chi,
}
$$

$$
\boxed{
A
=G_SM^{-1}\operatorname{diag}(E^{\mathrm{tot}}).
}
$$

**[已推导/需条件]** 条件是当前模型中 $G_S,G_D,e^{\mathrm{inc}}$ 对 map perturbation 固定，$M$ 在名义点可逆/局部解支稳定。接近内部共振或分支变化时不能无条件使用。

### 6.2 位姿 Jacobian $B$

固定 $\chi$，

$$
\delta y=(\delta G_S)j+G_S\delta j.
$$

若世界坐标中的 $G_D,D_\chi$ 不随平台位姿变化，则

$$
M\delta j=D_\chi\delta e^{\mathrm{inc}},
$$

所以对位姿坐标 $q$

$$
\boxed{
B_q
=\frac{\partial G_S}{\partial q}j
+G_SM^{-1}D_\chi\frac{\partial e^{\mathrm{inc}}}{\partial q}.
}
$$

第一项是 receive-path/measurement geometry 的变化，第二项是 Tx/illumination path 的变化。

更一般的完整候选式为

$$
\begin{aligned}
B_q
=&\;\frac{\partial G_S}{\partial q}j\\
&+G_SM^{-1}\Bigg[
D_\chi\frac{\partial e^{\mathrm{inc}}}{\partial q}
+D_\chi\frac{\partial G_D}{\partial q}j
+\frac{\partial D_\chi}{\partial q}E^{\mathrm{tot}}
\Bigg],
\end{aligned}
$$

并应再加入任何显式 direct-path、校准或接收模型的导数。

**[待验证]** 是否保留 $\partial_qG_D$ 与 $\partial_qD_\chi$，取决于采用世界固定网格还是机器人局部网格、目标是否静止、天线外参如何参数化；原简式不是坐标无关的“完整公式”。

远场相位直觉

$$
G(R)\sim a(R)e^{ikR},
\qquad
\delta G\approx ikG\,(\widehat R\cdot\delta x)
$$

说明短波长可提高几何敏感度，也会放大位姿/相位标定误差。此为 regime-dependent 近似，不是全波普适定理。

### 6.3 堆叠轨迹与 MIMO 分块

对时刻 $t$、发射通道 $p$，

$$
j_{t,p}=D_\chi(e^{\mathrm{inc}}_{t,p}+G_{D,t}j_{t,p}),
\qquad
y_{t,p}=G_{S,t}j_{t,p}+n_{t,p}.
$$

把所有 Rx 放在 $y_{t,p}$ 中，把所有 Tx/频率/极化列堆叠：

$$
Y_t=G_{S,t}J_t+N_t,
\qquad
J_t=[j_{t,1},\ldots,j_{t,P}],
$$

或

$$
\operatorname{vec}(Y_t)
=(I_P\otimes G_{S,t})\operatorname{vec}(J_t)+\operatorname{vec}(N_t).
$$

再按 $t$ 堆叠得到全局

$$
\delta y=A\delta\chi+B\delta X+n.
$$

关键结构：

- 共同静态未知量是 $\chi$，但 $j_{t,p}$ 随 pose/illumination 改变；full-wave 下不能简单把多个 $G_{S,t}$ 直接乘同一个 $\chi$。
- $A=[A_{t,p}]_{t,p}$ 是共同地图的 stacked Jacobian。
- 若每一时刻只有一个刚体位姿 $x_t\in\mathbb R^d$，则 $B$ 对该时刻的列数至多 $d$，即使 MIMO data 很高维；这产生低秩 nuisance 结构。
- 若 Tx/Rx 独立运动、外参/时钟/通道增益也未知，$B$ 应按 pose、Tx、Rx 与 calibration blocks 展开。
- 复数数据而物理 pose/contrast 参数为实数时，可将实部/虚部堆成实系统，或使用一致的 complex-real FIM；不能直接把复参数公式不加说明地套用。

### 6.4 gauge 与坐标选择

若全局平移/旋转同时作用于轨迹和地图而不改变观测，则联合 Jacobian 必有 gauge nullspace。绝对位置、方向、尺度、时钟或背景波速是否可辨识，取决于锚点、IMU/里程计、已知 beacon、外部坐标或材料/几何先验。

**[需补条件]** 所有 rank、inverse 和 CRB 结论必须先说明：

- gauge 已固定；或
- 只在 gauge quotient / identifiable support 上讨论；或
- 明确使用 Moore–Penrose 伪逆，但不把 gauge 方差解释成普通有限协方差。

---

## 7. $K_{\mathrm{IS}},K_{\mathrm{SLAM}},K_{\mathrm{eff}}$ 与信息损失

### 7.1 先白化测量

若 $n\sim\mathcal N_\mathbb C(0,\Sigma_n)$，令

$$
\widetilde A=\Sigma_n^{-1/2}A,
\qquad
\widetilde B=\Sigma_n^{-1/2}B.
$$

后文省略 tilde。若不白化，则所有投影都必须使用 $W=\Sigma_n^{-1}$ 所定义的 metric；不能在 colored noise 下直接用 Euclidean $BB^\dagger$。

### 7.2 已知位姿：ordinary inverse scattering information

$$
\boxed{K_{\mathrm{IS}}=A^*A.}
$$

对地图方向 $v$：

$$
v^*K_{\mathrm{IS}}v=\|Av\|^2.
$$

同一对象可叫 normal operator、Gauss–Newton Hessian 或 Gaussian 局部 Fisher information（差一个 convention-dependent 常数）。$A^*$ 是 adjoint/matched backprojection；它与 time reversal 有结构类比，但不应无条件当成同一个物理操作。

### 7.3 完全自由 pose：投影消元

对固定地图扰动 $z=A\delta\chi$，最能伪装它的 pose perturbation 解

$$
\min_{\delta X}\|z-B\delta X\|^2.
$$

pose 可解释部分为 $P_Bz=BB^\dagger z$，剩余为 $P_{B^\perp}z$。故

$$
\boxed{
K_{\mathrm{SLAM}}
=A^*P_{B^\perp}A.
}
$$

这里的最优 $\delta X^\star$ 是 **compensating/nuisance pose error**，不是要规划的“最优机器人路径”。观测 $y^{\mathrm{obs}}$ 固定；被最小化的是模型 residual。

**[已推导]**

$$
0\preceq K_{\mathrm{SLAM}}\preceq K_{\mathrm{IS}},
$$

因为 $0\preceq P_{B^\perp}\preceq I$。有限维 Hermitian 情况下按同一排序有

$$
\lambda_i(K_{\mathrm{SLAM}})\le\lambda_i(K_{\mathrm{IS}}).
$$

### 7.4 有 pose prior：Schur complement

若 pose deviation 有信息/penalty $J_X\succeq0$，联合信息块为

$$
\mathcal I=
\begin{bmatrix}
A^*A&A^*B\\
B^*A&B^*B+J_X
\end{bmatrix}.
$$

消元得到

$$
\boxed{
K_{\mathrm{eff}}
=A^*A-A^*B(B^*B+J_X)^\dagger B^*A.
}
$$

定义

$$
\boxed{
L_X=K_{\mathrm{IS}}-K_{\mathrm{eff}}
=A^*B(B^*B+J_X)^\dagger B^*A.
}
$$

在 $J_X\succ0$ 或满足 generalized Schur complement 的 range 条件时：

$$
0\preceq K_{\mathrm{SLAM}}
\preceq K_{\mathrm{eff}}
\preceq K_{\mathrm{IS}},
\qquad
L_X\succeq0.
$$

若 $J_X=\alpha J_0,\ J_0\succ0$：

$$
\alpha\to0^+:\ K_{\mathrm{eff}}\to K_{\mathrm{SLAM}},
\qquad
\alpha\to\infty:\ K_{\mathrm{eff}}\to K_{\mathrm{IS}}.
$$

若 prior 只约束部分 pose/gauge directions，极限必须逐子空间解释。

### 7.5 随机 pose 边缘化与 Schur complement 的关系

若 $\delta X\sim\mathcal N(0,C_X)$ 且独立于 $n$，在线性模型中边缘化位姿得到

$$
n_{\mathrm{eff}}=B\delta X+n,
\qquad
\Sigma_{\mathrm{eff}}=\Sigma_n+BC_XB^*.
$$

因此

$$
K_{\mathrm{marg}}=A^*\Sigma_{\mathrm{eff}}^{-1}A.
$$

在 $C_X^{-1}$ 存在时，Woodbury identity 使它与设 $J_X=C_X^{-1}$ 的 Gaussian joint model 等价。若 pose error 是未知固定 bias、与地图相关、非 Gaussian 或被重线性化，这种等价不再自动成立。

### 7.6 精确不可辨识与自由度损失

令

$$
\mathcal A=\operatorname{Range}(A),
\qquad
\mathcal B=\operatorname{Range}(B).
$$

**[已推导]** 有限维下

$$
\boxed{
\ker K_{\mathrm{SLAM}}
=\{v:Av\in\mathcal B\}.
}
$$

其中既含 $Av=0$ 的 ordinary inverse-scattering nullspace，也含 $Av\ne0$ 但完全可被 pose 模拟的新增 null directions。

并且

$$
\boxed{
\operatorname{rank}K_{\mathrm{SLAM}}
=\operatorname{rank}A
-\dim(\mathcal A\cap\mathcal B).
}
$$

无限维情形要处理 range 是否闭、投影到闭包及 Fredholm/compact operator 条件，不能直接照搬有限维维数式。

### 7.7 low-rank loss 的精确含义

$$
\boxed{
\operatorname{rank}L_X\le\operatorname{rank}B.
}
$$

若每个位姿只有 $d$ 个自由度，单时刻 pose-induced loss 至多 rank $d$；多时刻堆叠后上界随独立 pose/calibration 自由度增长。

“低秩”只说明 $L_X$ 在 **地图信息/Hessian 空间** 中由少数 pose-coupled 组合产生；它不保证：

- loss 的幅度小；
- 只改变相同数量的 eigenvalues；
- 重建误差在空间上局部；
- forward operator 本身是低秩；
- 加一个低秩 prior 就能恢复真实数据缺失的信息。

低秩 Hermitian perturbation 仍可旋转很多 eigenvectors，尤其在谱隙小或特征值簇内。

---

## 8. 连续信息保留谱与 principal-angle 几何

### 8.1 单个地图方向的保留率

对 $Av\ne0$：

$$
\boxed{
R(v)
=\frac{v^*K_{\mathrm{SLAM}}v}{v^*K_{\mathrm{IS}}v}
=\frac{\|P_{B^\perp}Av\|^2}{\|Av\|^2}
=\sin^2\angle(Av,\mathcal B).
}
$$

- $R=1$：与所有 pose-induced signatures 正交；
- $R=0$：完全可被 pose 伪装；
- $0<R<1$：部分保留。

angle 在 **data space** 中定义，不是在 map parameter space 中直接比较 $v$ 与 pose vector。

### 8.2 $\rho$ 的历史命名冲突

早期讨论曾定义 confounding amplitude

$$
\rho_{\mathrm{conf}}(v)=\frac{\|P_BAv\|}{\|Av\|}=\cos\theta,
$$

后续又用 $\rho$ 表示 information retention

$$
\rho_{\mathrm{ret}}(v)
=\frac{\|P_{B^\perp}Av\|^2}{\|Av\|^2}
=\sin^2\theta.
$$

二者关系：

$$
\boxed{
\rho_{\mathrm{ret}}=1-\rho_{\mathrm{conf}}^2.
}
$$

本文后续的 $\rho_i\in[0,1]$ 一律表示 **retention**；若引用早期“ρ 越大越混淆”的语句，指的是 $\rho_{\mathrm{conf}}$。

### 8.3 广义特征值，而不是逐项 eigenvalue 比值

直接用

$$
\lambda_i(K_{\mathrm{SLAM}})/\lambda_i(K_{\mathrm{IS}})
$$

并不可靠，因为两套 eigenvectors 未必对应。应在 $K_{\mathrm{IS}}$ 的 observable support 上解

$$
\boxed{
K_{\mathrm{SLAM}}v=\rho K_{\mathrm{IS}}v,
}
$$

或定义

$$
\boxed{
R_{\mathrm{op}}
=K_{\mathrm{IS}}^{\dagger/2}
K_{\mathrm{SLAM}}
K_{\mathrm{IS}}^{\dagger/2}.
}
$$

令

$$
U=A K_{\mathrm{IS}}^{\dagger/2},
$$

则在 observable support 上 $U^*U=I$，columns of $U$ 张成 $\mathcal A$，并有

$$
R_{\mathrm{op}}=U^*P_{B^\perp}U.
$$

故其特征值与 $\mathcal A,\mathcal B$ 的 principal angles 对应：

$$
\boxed{\rho_i=\sin^2\theta_i,\qquad 0\le\rho_i\le1.}
$$

维数不匹配时会有额外的 $1$ 或 $0$；应明确 principal-angle 的排序和支撑空间。

### 8.4 basis-independent 的多 feature 比较

不能任意逐列比较 $A_i$ 与 $B_j$。先作

$$
A=Q_AR_A,
\qquad
B=Q_BR_B,
$$

再对

$$
Q_A^*Q_B
=U_\theta\operatorname{diag}(\cos\theta_i)V_\theta^*
$$

做 SVD/CS decomposition，得到 basis-independent coupling。

### 8.5 absolute observability 与 relative retention 是两条轴

$K_{\mathrm{IS}}$ 衡量地图方向本来能产生多强的数据，$R_{\mathrm{op}}$ 衡量其中有多少在 pose 消元后保留。于是：

| absolute strength | retention 高 | retention 低 |
|---|---|---|
| 强 | 数据强且 pose 可分离，最理想 | 数据强但 pose-confounded，应优先改 geometry/trajectory/pose prior |
| 弱 | 不混 pose，但 SNR/孔径不足 | 弱且混淆，主要由 prior/regularization 决定 |

若 $v_i$ 是 $K_{\mathrm{IS}}$ mode，则方向性的 Rayleigh information 可写成

$$
I_{\mathrm{eff},i}=\sigma_i^2\rho_{\mathrm{ret},i},
$$

但 $v_i$ 通常不会继续是 $K_{\mathrm{SLAM}}$ eigenvector。对前 $r$ 个 SOM modes，真正的 reduced operator 是

$$
K_{\mathrm{SLAM}}^{(r)}
=V_r\Sigma_r
[I-U_r^*P_BU_r]
\Sigma_rV_r^*,
$$

其中 $U_r^*P_BU_r$ 一般非对角，会耦合不同 SOM modes。只有近似共同对角化时，逐 mode 的 $(\sigma_i,\rho_i)$ 分类才精确。

---

## 9. Chen SOM / Two-fold SOM 与 pose-defect 的分层合并

### 9.1 $G_S/G_D$ 的“两折”仍在 current space

设

$$
G_S=U_S\Sigma_SV_S^*,
\qquad
\mathcal J=V_S^+\oplus V_S^-,
$$

并把内部 state/domain 算子所诱导的 current-space 结构记作

$$
\mathcal J=V_D^+\oplus V_D^-.
$$

S2 中整理出的 Two-fold SOM 工作空间为

$$
\boxed{
\mathcal J_{\mathrm{TSOM}}
=
V_S^+
\oplus
(V_S^-\cap V_D^+)
\oplus
(V_S^-\cap V_D^-).
}
$$

其计算上保留的核心通常写成

$$
\boxed{
J
=
V_S^+\alpha^+
+
B_{S^-D^+}\beta,
}
$$

而 $V_S^-\cap V_D^-$ 是 data equation 与 domain/state equation 两边都不利的 double-negative part，因而被丢弃或强正则化。

**[需补条件]** 上式是七条讨论内部形成的术语骨架；$V_D^\pm$ 到底来自 $G_D$、其 normal operator、或另一种材料加权算子，必须按所引用的 Chen 论文版本核对。本文不把跨版本记号统一当成已验证事实。

### 9.2 精确交集与实际投影近似不是一回事

讨论中明确区分：

$$
V_S^-\cap V_D^+
$$

是几何上的精确交集；实际算法常用

$$
\boxed{
B_{S^-D^+}\beta
\approx
V_S^-(V_S^-)^*V_D^+\beta
=P_{S^-}V_D^+\beta
}
$$

来构造近似 basis。若 $P_{S^-}$ 与 $P_{D^+}$ 不交换，投影结果不一定等于精确交集，且顺序会影响结果。

**[冲突/修正]** 曾提出把 $V_S^+\cap V_D^-$ 单独交给机器学习；后续讨论认为这不是原始 TSOM 的稳定语义，最多是新的算法提案，不能反向写成 Chen Two-fold SOM 的定义。

### 9.3 第三层应是 pose defect，而不是把三个投影硬拆成八格

定义 TSOM 保留 basis

$$
U_T=[\,V_S^+,\ B_{S^-D^+}\,],
\qquad
\delta J=U_T\delta c.
$$

在这一 reduced current space 中，联合线性化写成

$$
\boxed{
\delta r=A_T\delta c+B\delta x+n,
}
$$

其中 $A_T$ 把 retained-current coefficient 映到白化 data space，$B$ 是 pose tangent。对应

$$
K_0=A_T^*A_T,
$$

$$
\boxed{
K_P
=K_0-A_T^*B(B^*B+\Lambda_x)^\dagger B^*A_T.
}
$$

因此最稳妥的 hierarchy 是：

1. **S-fold**：由 $G_S$ 区分 data-visible 与 complementary current；
2. **D-fold**：在 complementary current 内利用 state/domain structure；
3. **pose-fold**：在 TSOM 已保留的空间内，研究 finite-dimensional pose nuisance 造成的相对谱缺陷。

这就是讨论所称的 **TriSpace/Tri-fold 研究框架**。前两层主要属于 current recovery，第三层属于 joint map/current–pose distinguishability。

**[冲突/修正]** 对 $P_S^\pm,P_D^\pm,P_P^\pm$ 作对称的 $2\times2\times2$ “八格精确分解”一般不成立，因为

$$
[P_S,P_D]\ne0,
\qquad
[P_D,P_P]\ne0
$$

是常态。只有在相关投影交换或存在共同 reducing subspaces 时，八格才是精确正交分解。

若只是寻找“同时接近三个 plus subspaces”的向量，曾提出看

$$
M_{+++}=P_S^++P_D^++P_P^+
$$

的接近 $3$ 的 eigenmodes。**[直觉/假设]** 这是非交换投影的近似交集评分，不是 TSOM 定理，也不是唯一构造。

### 9.4 reduced pose-defect theorem

在 $K_0$ 的 observable support 上解

$$
\boxed{
K_Pv=\rho K_0v.
}
$$

令

$$
A_T=Q_AR_A,\qquad B=Q_BR_B,
\qquad Z=Q_A^*Q_B,
$$

在无 pose prior 的精确投影情形，

$$
C_P=ZZ^*=Q_A^*P_BQ_A,
\qquad
R=I-C_P.
$$

若 $C_Pu_i=\mu_i u_i$，则

$$
\boxed{\rho_i=1-\mu_i=\sin^2\theta_i.}
$$

这里 $\theta_i$ 是 retained data tangent $\operatorname{Range}(A_T)$ 与 pose tangent $\operatorname{Range}(B)$ 的 principal angles。

由于

$$
\operatorname{rank}C_P\le p:=\operatorname{rank}B,
$$

得到比普通 interlacing 更强的结论：

$$
\boxed{
\#\{i:\rho_i\ne1\}\le p,
}
$$

即在精确算术、无 prior、固定 observable support 下，至少 $r-p$ 个 generalized retention eigenvalues 恰为 $1$。这并不表示普通 $K_P$ 的 eigenvalues 只有 $p$ 个变化；低秩扰动仍可使许多普通 eigenpairs 改变。

定义 pose-defect subspace

$$
\mathcal D_P=\operatorname{Range}(I-R)
$$

后，

$$
\dim\mathcal D_P\le p,
$$

而 $\ker(I-R)$ 是对该一阶 pose nuisance 完全 robust 的 complement。

### 9.5 $G_S/G_D$ 四类与 pose 三类的语义表

在原讨论的工作语义中：

| current-space 类别 | 原始含义 | pose 层可能发生什么 |
|---|---|---|
| $V_S^+\cap V_D^+$ | data 强、state/domain 也支持 | 仍可能与 pose tangent 混淆 |
| $V_S^+\cap V_D^-$ | data 可见但 state/domain 支持较弱或语义依版本 | 不能仅凭标签判为可恢复 |
| $V_S^-\cap V_D^+$ | data 弱，但 state/domain 可用于约束 | TSOM 的主要补偿空间；再检查 pose defect |
| $V_S^-\cap V_D^-$ | 双重不利 | 通常舍弃/强正则化 |

pose 层更适合按连续 retention 分成：

$$
\mathcal E_G:\rho\approx0,\qquad
\mathcal E_M:0<\rho<1,\qquad
\mathcal E_R:\rho\approx1.
$$

- $\mathcal E_G$：gauge-like / strongly confounded；
- $\mathcal E_M$：partially retained；
- $\mathcal E_R$：locally pose-robust。

**[需补条件]** 这些是 spectral clusters/projectors，不应未经阈值敏感性分析就当成三个自然常数维空间。

### 9.6 map-space pullback：不要直接把 current modes 当地图 modes

若 $T_\chi=\partial J/\partial\chi$，current-level information 可写成

$$
K_S^J=A_J^*A_J,
\qquad A_J=G_S
$$

（白化与多帧堆叠省略），pose 消元后为

$$
K_{\mathrm{pose}}^J
=K_S^J-A_J^*P_BA_J.
$$

地图空间需要 pull back：

$$
\boxed{
K_{\chi,S}=T_\chi^*K_S^JT_\chi,
\qquad
K_{\chi,\mathrm{pose}}
=T_\chi^*K_{\mathrm{pose}}^JT_\chi.
}
$$

因此 $V_S^\pm$ 是 current-space basis，$K_{\mathrm{eff}}$ 的 eigenvectors 是 map/reduced-coordinate modes；二者只有通过 $T_\chi$ 或明确 reduced parameterization 才能联系。

---

## 10. 低秩信息损失、谱扰动与 regularization

### 10.1 pose loss 的 factorization

令

$$
C_X=B^*B+J_X.
$$

则

$$
\boxed{
L_X
=K_{\mathrm{IS}}-K_{\mathrm{eff}}
=A^*BC_X^\dagger B^*A
=H_X^*H_X,
}
$$

其中

$$
H_X=C_X^{\dagger/2}B^*A.
$$

于是

$$
\boxed{
\operatorname{rank}L_X
=\operatorname{rank}H_X
\le\operatorname{rank}B,
}
$$

$$
\boxed{
\operatorname{tr}L_X
=\|C_X^{\dagger/2}B^*A\|_F^2,
\qquad
\lambda_{\max}(L_X)
=\|C_X^{\dagger/2}B^*A\|_2^2.
}
$$

普通 $\det L_X$ 在地图维数大于 pose 维数时通常为零；更有意义的是 $\operatorname{pdet}L_X$，或在任务子空间 $\mathcal S$ 上的 determinant ratio：

$$
\frac{\operatorname{pdet}_{\mathcal S}K_{\mathrm{eff}}}
{\operatorname{pdet}_{\mathcal S}K_{\mathrm{IS}}}.
$$

### 10.2 affected directions 与 destroyed degrees of freedom

无 prior 时：

- $\operatorname{rank}L_X$ 等于被 pose 投影触及的独立地图组合数，即 principal angle 小于 $\pi/2$ 的数量；
- 真正新增的 exact nullity / destroyed DoF 为

$$
\boxed{
\operatorname{rank}K_{\mathrm{IS}}
-\operatorname{rank}K_{\mathrm{SLAM}}
=\dim(\operatorname{Range}(A)\cap\operatorname{Range}(B)),
}
$$

即 principal angle 恰为 $0$ 的数量。

所以“受影响”不等于“完全摧毁”。一个 rank-$p$ loss 可以让至多 $p$ 个相对 retention directions 偏离 $1$，但它可能经 eigenvector rotation 影响所有 ordinary modal interpretations。

### 10.3 rank-$p$ interlacing

若 $\operatorname{rank}L_X\le p$，则有限维 Hermitian interlacing 给出

$$
\boxed{
\lambda_{i+p}(K_{\mathrm{IS}})
\le
\lambda_i(K_{\mathrm{eff}})
\le
\lambda_i(K_{\mathrm{IS}}).
}
$$

**[已推导/已接受]** 对固定有限 $p$，finite-rank nuisance 不改变 compact inverse problem 的渐近 eigenvalue-decay class，只能作有限索引位移。

**[需补条件]** 长轨迹若为每个时刻引入独立 pose block，$p\sim d_xT$ 会随数据规模增长，此时不能援引“固定 finite-rank perturbation 不改变渐近类”的简化。

### 10.4 pose prior 的 weighted shrinkage

当 $J_X\succ0$，令

$$
D=BJ_X^{-1/2}.
$$

Woodbury identity 给出

$$
\boxed{
I-B(B^*B+J_X)^{-1}B^*
=(I+DD^*)^{-1}.
}
$$

若

$$
D=Z\operatorname{diag}(\eta_\ell)W^*
$$

是 SVD，则沿 pose-data singular directions 的保留因子是

$$
\boxed{
\frac{1}{1+\eta_\ell^2}.
}
$$

因此 finite prior 下的 normalized operator 是

$$
R_X=Q_A^*(I+DD^*)^{-1}Q_A,
$$

它是 prior-weighted shrinkage，通常不再严格等于某组 ordinary principal angles 的 $\sin^2$。

讨论提出一个连续 confounding dimension：

$$
\boxed{
d_{\mathrm{conf}}(J_X)
=
\sum_\ell
\frac{\eta_\ell^2}{1+\eta_\ell^2}
\|Q_A^*z_\ell\|^2.
}
$$

它同时反映 pose-data coupling 的强度、pose prior 和 map tangent 对应方向。

由变分式

$$
\boxed{
u^*K_{\mathrm{eff}}u
=
\min_z
\bigl(
\|Au-Bz\|^2+z^*J_Xz
\bigr)
}
$$

还可直接看出 prior 单调性：

$$
J_{X,1}\succeq J_{X,0}
\quad\Longrightarrow\quad
K_{\mathrm{eff}}(J_{X,1})
\succeq
K_{\mathrm{eff}}(J_{X,0}).
$$

### 10.5 regularization 不等于“创造信息”

若加入地图先验/regularizer curvature $K_{\mathrm{prior}}$，为避免与 pose Schur complement 的 $K_{\mathrm{eff}}$ 冲突，本文记

$$
\boxed{
K_{\mathrm{post}}=K_{\mathrm{eff}}+K_{\mathrm{prior}}.
}
$$

讨论形成的共识：

- TSVD 产生严格 low-rank estimator；
- nuclear-norm penalty 倾向低秩；
- Tikhonov $\lambda I$ 是 full-rank，但通过 filter factor $\sigma_i^2/(\sigma_i^2+\lambda)$ 形成 effective-rank 截断；
- TV 的 Hessian/局部曲率通常不低秩；
- modal prior $U_r\Gamma U_r^*$ 可以是显式低秩。

**[冲突/修正]** “regularization 本质上是低秩”被否定。更准确的是：某些 regularizers 通过滤波或结构先验减少 effective DoF；它们改善偏差—方差和数值稳定性，但不增加测量本身的 Fisher information。

还需区分两种“低维”：

1. $L_X=A^*BC_X^\dagger B^*A$ 是 Schur elimination 后的低秩 **信息损失矩阵**；
2. 几何误差族

   $$
   \Delta K=\sum_j K_j\,\delta x_j
   $$

   只位于由少数坐标张成的 **低维结构化矩阵族**，每个 $K_j$ 本身未必低秩。

### 10.6 Schur/Gauss–Newton reduced solve

对 joint Gauss–Newton normal equation，消去 pose increment 后：

$$
K_P\delta c=-\widetilde g_c,
$$

$$
\widetilde g_c
=g_c-A_T^*B(B^*B+\Lambda_x)^{-1}g_x,
$$

$$
\delta x
=-(B^*B+\Lambda_x)^{-1}
(g_x+B^*A_T\delta c).
$$

令 $A_T=Q_AR_A$、$w=R_A\delta c$，并设

$$
Z=Q_A^*B(B^*B+\Lambda_x)^{-1/2}.
$$

只需在 pose-sized 小矩阵

$$
H=Z^*Z
$$

上求谱；若 $Hz_j=\mu_jz_j$，则

$$
u_j=\frac{Zz_j}{\sqrt{\mu_j}},
\qquad
\rho_j=1-\mu_j.
$$

无 exact gauge 且加入 relative Levenberg–Marquardt $K_P+\lambda K_0$ 时，白化坐标更新为

$$
\boxed{
w
=-\frac{1}{1+\lambda}b
-
U_P
\operatorname{diag}
\left(
\frac{1}{\rho_j+\lambda}
-
\frac{1}{1+\lambda}
\right)
U_P^*b.
}
$$

其意义是：pose-robust complement 只需统一缩放，而所有额外修正都限制在至多 $p$ 维 defect subspace 中。

若取前 $m$ 个 defect modes 作近似，讨论给出的尾项 operator-norm 指标为

$$
\boxed{
\frac{\mu_{m+1}}
{(1-\mu_{m+1}+\lambda)(1+\lambda)}.
}
$$

可用它而不是任意 mode count 作为 nested approximation 的停止准则。

### 10.7 incremental low-rank basis update

若轨迹更新导致

$$
Z_{k+1}=Z_k+\Delta Z,\qquad Z_k=U_kR_k,
$$

定义

$$
A_k=U_k^*\Delta Z,
\qquad
E_k=(I-U_kU_k^*)\Delta Z=Q_\perp R_\perp.
$$

则有精确的小块表示

$$
Z_{k+1}=U_k(R_k+A_k)+Q_\perp R_\perp.
$$

对右侧小矩阵作 SVD/QR 即可更新 defect basis。fixed-rank、positive-gap 条件下的微分形式为

$$
\dot U=(I-UU^*)\dot Z\,V\Sigma^{-1}.
$$

当最小保留奇异值趋近零时该式变得病态，必须转入 rank-event 处理，而不能继续假定固定维 Grassmann 流。

---

## 11. 谱连续性、Grassmann/Flag 争议与 rank-changing dynamics

### 11.1 应追踪 projector/cluster，不应执着于单个 eigenvector

单个 eigenvector 有相位/符号和簇内旋转不唯一性。若 $\Gamma$ 包围一组与其余谱隔开的 eigenvalues，谱投影可写成

$$
\boxed{
P(X)
=\frac{1}{2\pi i}
\oint_\Gamma
(zI-R(X))^{-1}\,dz.
}
$$

只要该 cluster 与外部谱保持 gap，$P(X)$ 可连续甚至光滑，即使 cluster 内 eigenvalues 重合或交换。其导数为

$$
\dot P=\frac{1}{2\pi i}\oint_\Gamma (zI-R)^{-1}\dot R(zI-R)^{-1}\,dz.
$$

选择 basis $U$ 时存在 $U\sim UQ$ 的 **basis gauge**。平行移动规范

$$
U^*\dot U=0,
\qquad
\dot U=\dot P\,U
$$

用于消除纯 basis rotation。它与 SLAM 的物理 $SE(2)/SE(3)$ gauge 完全不同。

### 11.2 continuity 的真正条件

**[冲突/修正]** “只要 $0<\rho_i<1$，子空间就连续”不成立。真正需要的是：

- operator $R(X)$ 对 $X$ 连续/可微；
- 被追踪 cluster 的 rank 保持；
- cluster 与其余谱之间有正 spectral gap。

$\rho_i=0$ 或 $1$ 本身不会破坏 projector continuity；问题来自 gap closure、rank change 或阈值穿越。特别是 $\rho=1$ 往往有高重数，逐个追踪 robust eigenvectors 没有不变意义；应追踪低秩 defect projector，再取它的正交 complement。

### 11.3 hard classification 的稳定裕量

若用阈值 $\tau_G<\tau_R$ 把谱分成 gauge/mixed/robust 三簇，可定义

$$
\boxed{
m_i
=\min\{|\rho_i-\tau_G|,\ |\rho_i-\tau_R|\}.
}
$$

由 Weyl inequality，若 operator perturbation 满足

$$
\|\Delta R\|<m_i,
$$

则该 eigenvalue 不会跨越分类阈值。若接近阈值，不应把标签当成稳定物理属性。

正常的 Gauss–Newton/Fisher 构造保证 $0\le\rho_i\le1$。若数值上出现 $\rho<0$ 或明显 $\rho>1$，应优先检查：

- whitening/metric 是否一致；
- 是否错误处理 singular support；
- generalized eigenproblem 是否数值病态；
- 使用的是否 indefinite full-Newton Hessian，而非 PSD information operator。

### 11.4 Grassmann 与 Flag 的最终定位

原讨论最终接受的是 **rank-stratified, piecewise-smooth** 图景：

- 在固定 rank 与 gap 的区间，defect/retained subspace 在相应 Grassmann manifold 上平滑演化；
- 在 singular value 穿零、cluster 合并或 hard threshold crossing 时，发生离散 rank event；
- 事件前后需要重新初始化维数、使用 hysteresis 或 tracking window。

Flag manifold 可用于在临界区同时保存 nested subspaces，例如

$$
\mathcal E_G\subset
\mathcal E_G\oplus\mathcal E_M
\subset\mathcal H,
$$

但 **[冲突/修正]** 它不是解决 rank change 的自动定理，也不是主线所必需。最终意见是把 Flag 当作可选 transition buffer / bookkeeping device，而把谱 projector、gap 与 rank events 作为核心数学。

### 11.5 eigenvalue 与 generalized-retention 导数

若

$$
K_{\mathrm{eff}}v_i=\lambda_iv_i,\qquad \|v_i\|=1
$$

且 $\lambda_i$ 单，则

$$
\boxed{
D\lambda_i[h]
=v_i^*DK_{\mathrm{eff}}[h]v_i.
}
$$

取相位规范 $v_i^*\dot v_i=0$，

$$
\dot v_i
=
\sum_{j\ne i}
v_j
\frac{v_j^*\dot K_{\mathrm{eff}}v_i}
{\lambda_i-\lambda_j}.
$$

因此 eigenvector sensitivity 与 inverse gap 成正比。

对 generalized retention mode，若

$$
K_{\mathrm{eff}}v_i=\rho_iK_{\mathrm{IS}}v_i,
\qquad
v_i^*K_{\mathrm{IS}}v_i=1,
$$

则

$$
\boxed{
\dot\rho_i
=
v_i^*
(\dot K_{\mathrm{eff}}-\rho_i\dot K_{\mathrm{IS}})
v_i.
}
$$

重特征值时不能使用单一梯度；一阶 splitting 由 eigenspace $V_E$ 上的压缩矩阵

$$
V_E^*DK[h]V_E
$$

的 eigenvalues 决定。

---

## 12. 位姿不确定性的三种模型与必须补上的实数化

### 12.1 复测量、实位姿：whitening 后仍须 realify

若

$$
n\sim\mathcal{CN}(0,\Sigma),
\qquad
\widehat A=\Sigma^{-1/2}A,\quad
\widehat B=\Sigma^{-1/2}B,
$$

且 $\chi,X$ 都是实参数，则真实 Fisher geometry 应使用

$$
\boxed{
A_{\mathbb R}
=\sqrt2
\begin{bmatrix}
\Re\widehat A\\
\Im\widehat A
\end{bmatrix},
\qquad
B_{\mathbb R}
=\sqrt2
\begin{bmatrix}
\Re\widehat B\\
\Im\widehat B
\end{bmatrix}.
}
$$

若 contrast 是复数，必须把 $\chi_R,\chi_I$ 当成两个实参数块：

$$
\boxed{
A_{\mathbb R}
=\sqrt2
\begin{bmatrix}
\Re\widehat A&-\Im\widehat A\\
\Im\widehat A& \Re\widehat A
\end{bmatrix}.
}
$$

随后所有 projector、range、principal angle 与 Schur complement 都在这个 realified space 中计算。为节省符号，其他章节仍写 $A,B,*$。

**[冲突/修正]** 若直接在 complex-linear space 中使用 $BB^\dagger$，等于允许不物理的 complex pose increment，通常会扩大 nuisance range 并高估 map–pose confounding。

### 12.2 模型 A：确定性 pose nuisance

对

$$
\delta y=A\delta\chi+B\delta X+n
$$

联合估计 $\delta\chi,\delta X$，再消去 $\delta X$。无 prior 得

$$
K_{\mathrm{SLAM}}=A^*P_{B^\perp}A;
$$

有 $J_X$ 得

$$
K_{\mathrm{eff}}
=A^*A-A^*B(B^*B+J_X)^\dagger B^*A.
$$

这回答：“允许 pose 自由补偿后，地图的一阶残差/信息还剩多少？”

### 12.3 模型 B：随机 pose error 的边缘化

若

$$
\delta X\sim\mathcal N(0,C_X),
\qquad
\delta X\perp n,
$$

则条件于 $\delta\chi$ 的有效噪声为

$$
\boxed{
\Sigma_{\mathrm{eff}}
=\Sigma_n+BC_XB^*,
}
$$

地图信息为

$$
\boxed{
K_{\mathrm{marg}}
=A^*\Sigma_{\mathrm{eff}}^{-1}A.
}
$$

当 $C_X^{-1}$ 存在时，Woodbury identity 使它与 Gaussian joint model 中 $J_X=C_X^{-1}$ 的 Schur complement 等价。

**[需补条件]** 若 pose error 是未知固定 bias、与地图相关、非 Gaussian，或每轮 nonlinear relinearization 改变其分布，这一等价不自动成立。

### 12.4 模型 C：外层有界执行/模型失配

若部署时真实轨迹为

$$
X_{\mathrm{true}}=X+\Delta X,
\qquad
\|\Delta X\|\le\varepsilon,
$$

目标是分析

$$
\min_{\|\Delta X\|\le\varepsilon}
\lambda_r(K_{\mathrm{eff}}(X+\Delta X)).
$$

这里 $\Delta X$ 不是内层 Schur complement 中为解释 residual 而优化的 nuisance variable；它是标称轨迹相对真实执行/建模的 perturbation。

**[冲突/修正]** 内层 $J_X$ 与外层 $\Delta X$ 可以同时存在，但必须说明各自代表什么，否则可能把同一来源的不确定性计算两次。

### 12.5 独立 pose blocks 与 motion coupling

若各时刻 pose 完全独立且无 motion prior，

$$
B=\operatorname{blkdiag}(B_1,\ldots,B_T),
$$

$$
P_{B^\perp}
=\operatorname{blkdiag}
(P_{B_1^\perp},\ldots,P_{B_T^\perp}),
$$

从而

$$
\boxed{
K_{\mathrm{SLAM}}
=\sum_{t=1}^T
A_t^*P_{B_t^\perp}A_t.
}
$$

并有

$$
\ker K_{\mathrm{SLAM}}
=
\bigcap_t
\{u:A_tu\in\operatorname{Range}(B_t)\},
$$

$$
\boxed{
\operatorname{rank}K_{\mathrm{SLAM}}
\le
\sum_t
(m_t-\operatorname{rank}B_t),
}
$$

其中 $m_t$ 是 realified 单帧数据维数。

若一帧只有一个复标量、realified 维数为 $2$，且自由 pose Jacobian $B_t$ 对这两个数据分量满行秩，则 $P_{B_t^\perp}=0$：该帧单独不提供一阶有效地图信息。多通道 MIMO、跨时刻 odometry/IMU prior 或共享轨迹参数此时可能是可辨识性的必要条件，而不只是提高精度。

有 motion prior 时，$J_X$ 通常是 block-banded/sparse，不应再把全局 pose projector 写成逐帧独立 block diagonal。

### 12.6 物理 gauge theorem

若均匀背景中同时对地图与全轨迹施加 $g\in SE(2)$ 不改变相对物理几何，

$$
F(g\cdot\chi,g\cdot X)=F(\chi,X).
$$

沿 Lie algebra generator $\xi$ 求导：

$$
\boxed{
A\xi_\chi+B\xi_X=0.
}
$$

所以

$$
A\xi_\chi=-B\xi_X\in\operatorname{Range}(B),
\qquad
\boxed{\xi_\chi\in\ker K_{\mathrm{SLAM}}.}
$$

二维平移 $v\in\mathbb R^2$ 的 map generator 为

$$
\delta\chi_v(r)=-v^T\nabla\chi(r),
\qquad
\delta p_t=v,\quad\delta\theta_t=0.
$$

二维旋转 generator 为

$$
\delta\chi_\omega(r)
=-\omega(Jr)^T\nabla\chi(r),
$$

$$
\delta p_t=\omega Jp_t,
\qquad
\delta\theta_t=\omega.
$$

一般非对称场景最多给出三个 $SE(2)$ gauge directions。已知 anchor、非均匀已知背景或固定世界边界会部分/全部破坏该不变性；场景本身的旋转对称性也可能使某个 map generator 退化。

**[已推导/需条件]** exact gauge 必然给出 $\rho=0$。反过来不成立：$\rho=0$ 也可能是某个局部 map–pose tangent coincidence，并不一定来自全局群不变性。

---

## 13. $DK_{\mathrm{eff}}$、projector derivative 与 robust trajectory

### 13.1 完整一阶导数

令

$$
C=B^*B+J_X,\qquad H=C^{-1},\qquad W=I-BHB^*,
$$

暂假设 $C$ 可逆，则

$$
K_{\mathrm{eff}}=A^*WA.
$$

沿轨迹方向 $h$：

$$
\boxed{
\dot K_{\mathrm{eff}}
=
\dot A^*WA
+
A^*W\dot A
+
A^*\dot WA,
}
$$

$$
\boxed{
\dot W
=-\dot BHB^*
-BH\dot B^*
+BH\dot C HB^*,
}
$$

$$
\boxed{
\dot C
=\dot B^*B+B^*\dot B+\dot J_X.
}
$$

这是把 geometry 对 map Jacobian、pose tangent 与 motion prior 的全部一阶影响分开的总公式。

### 13.2 full-wave $\dot A$ 与 resolvent sensitivity

由

$$
A=G_SM^{-1}D_E
$$

得到

$$
\boxed{
\dot A
=
\dot G_SM^{-1}D_E
+
G_S\dot M^{-1}D_E
+
G_SM^{-1}D_{\dot E}.
}
$$

其中

$$
\dot M^{-1}=-M^{-1}\dot M M^{-1}.
$$

若 $\dot M=-D_\chi\dot G_D$，则

$$
\dot M^{-1}
=M^{-1}D_\chi\dot G_DM^{-1}.
$$

同时

$$
\dot E=\dot e^{\mathrm{inc}}+\dot G_Dj+G_D\dot j.
$$

若 $G_D$ 与 pose 无关，

$$
\dot j=M^{-1}D_\chi\dot e^{\mathrm{inc}},
\qquad
\dot E=(I-G_DD_\chi)^{-1}\dot e^{\mathrm{inc}}.
$$

基本 bound 为

$$
\boxed{
\|A_t\|
\le
\|G_{S,t}\|
\|M_t^{-1}\|
\|E_t^{\mathrm{tot}}\|_\infty.
}
$$

pose Jacobian 则满足

$$
\begin{aligned}
\|B_th\|
\le{}&
\|DG_{S,t}[h]\|\,\|j_t\|\\
&+
\|G_{S,t}\|\|M_t^{-1}\|\|D_\chi\|
\left(
\|De_t^{\mathrm{inc}}[h]\|
+
\|DG_D[h]\|\,\|j_t\|
\right).
\end{aligned}
$$

当 $\sigma_{\min}(M_t)\to0$ 时，map information 可能因场增强而增大，但 $A,B,DK_{\mathrm{eff}}$ 和模型失配敏感度也可能被 resolvent 放大，局部线性化半径反而缩小。

### 13.3 $\dot B$ 需要二阶 forward derivative

因为 $B=D_XF$，

$$
\boxed{
\dot B=D_X^2F[h,\cdot].
}
$$

在 $G_D$ 固定时，若

$$
Bh=DG_S[h]j+G_SDj[h],
$$

则沿另一个方向 $q$

$$
\boxed{
\begin{aligned}
(DB[q])h
={}&D^2G_S[q,h]j\\
&+DG_S[h]Dj[q]\\
&+DG_S[q]Dj[h]\\
&+G_SD^2j[q,h].
\end{aligned}
}
$$

所以只使用 first-order forward Jacobian 而忽略 $\dot B$ 是 **frozen nuisance-subspace approximation**，可以作为近似，但必须明示。

### 13.4 无 prior 时的 projector derivative

若 $B$ 在邻域保持常秩，

$$
P_B=BB^\dagger
$$

的导数为

$$
\boxed{
\dot P_B
=
P_{B^\perp}\dot B B^\dagger
+
(P_{B^\perp}\dot B B^\dagger)^*.
}
$$

从而

$$
\boxed{
\|\dot K_{\mathrm{SLAM}}\|
\le
2\|A\|\|\dot A\|
+
2\|A\|^2\|B^\dagger\|\|\dot B\|.
}
$$

当 $B$ 的最小非零 singular value 很小时，$\|B^\dagger\|$ 变大，pose tangent subspace 对轨迹扰动高度敏感。这给出了 Grassmann 层“子空间突然不稳”的一个具体来源。

### 13.5 robust objective 与适用条件

定义

$$
\gamma_r(X,\varepsilon)
=
\min_{\|\Delta X\|\le\varepsilon}
\lambda_r(K_{\mathrm{eff}}(X+\Delta X)).
$$

若 $\lambda_r$ 是单特征值、与相邻谱有正 gap、$K_{\mathrm{eff}}$ 二阶可微，且不确定集为完整范数球，则

$$
\boxed{
\gamma_r(X,\varepsilon)
=
\lambda_r(X)
-
\varepsilon\|\nabla_X\lambda_r(X)\|_*
+
O(\varepsilon^2).
}
$$

这里 $\|\cdot\|_*$ 是所选轨迹扰动范数的对偶范数。若动力学/避障把扰动限制在 tangent cone $\mathcal T_X$，惩罚项应改成在该 cone 上的 support function，而不是完整对偶范数。

若只知 Lipschitz bound

$$
\|DK(X)[h]\|\le L_K\|h\|,
$$

则不要求 simple eigenvalue 的保守下界是

$$
\boxed{
\gamma_r(X,\varepsilon)
\ge
\lambda_r(K_{\mathrm{eff}}(X))
-
\varepsilon L_K
+
O(\varepsilon^2).
}
$$

这比单 eigenvalue gradient 更保守，但可跨越 eigenvalue crossing。

---

## 14. Born、far-field、多频率与信息“reach”

### 14.1 Born 的空背景退化

Born 模型写成

$$
y=A(X)\chi,
\qquad
A(X)=G_S(X)D_{e^{\mathrm{inc}}(X)}.
$$

在非零名义场景 $\chi_0$ 处，

$$
\boxed{
B_\ell
=\frac{\partial A(X)}{\partial X_\ell}\chi_0,
}
$$

$$
\operatorname{Range}(B)
=
\operatorname{span}
\{A_{X_1}\chi_0,\ldots,A_{X_p}\chi_0\}.
$$

若 $\chi_0=0$，则 $B=0$，于是 joint first-order Jacobian 给出

$$
K_{\mathrm{SLAM}}=K_{\mathrm{IS}}.
$$

**[冲突/修正]** 这不表示空场景附近 pose error 无害。联合展开为

$$
\boxed{
F(\delta\chi,X_0+\delta X)
=
A(X_0)\delta\chi
+
D_XA(X_0)[\delta X]\,\delta\chi
+
O(\|\delta\chi\|\|\delta X\|^2),
}
$$

主导几何失配是 map–pose 双线性项。因此：

- 非零 $\chi_0$ 周围的 $A\delta\chi+B\delta X$ 研究局部联合可辨识；
- $\chi_0=0$ 周围的 $A(X_0+\delta X)-A(X_0)$ 研究 sensing-operator uncertainty。

二者相关但不是同一个一阶问题。

### 14.2 far-field Fourier sampling 结构

plane-wave incidence 与 far-field observation 下，

$$
y_{\alpha,k}
=c_{\alpha,k}\widehat\chi(q_{\alpha,k}),
\qquad
q_{\alpha,k}=k(\widehat r_\alpha-d_\alpha)
$$

（符号取决于 Fourier convention）。对 pose coordinate $X_\ell$：

$$
\boxed{
\frac{\partial y_{\alpha,k}}{\partial X_\ell}
=
\frac{\partial\log c_{\alpha,k}}{\partial X_\ell}y_{\alpha,k}
+
c_{\alpha,k}
\nabla_q\widehat\chi(q_{\alpha,k})^T
\frac{\partial q_{\alpha,k}}{\partial X_\ell}.
}
$$

pose tangent 因此包括：

1. phase/amplitude modulation；
2. Fourier sampling location displacement。

远距离微小平移以前者为主；有限距离或大 aperture motion 会改变 $q$，因此访问 $\nabla_q\widehat\chi$。直线/小角度轨迹产生 missing wedge，圆周与多视角扩大 $q$-coverage。

全局平移的显式验证是

$$
\chi_v(r)=\chi(r-v)
\quad\Longrightarrow\quad
\widehat\chi_v(q)
=e^{-iq^Tv}\widehat\chi(q),
$$

故

$$
\delta y(q)=-i(q^Tv)y(q),
$$

与 sensor reference frame 的反向平移造成同一 data tangent，呼应 $SE(2)$ gauge。

### 14.3 距离、波长与 absolute/relative tradeoff

二维 $kR\gg1$ 时

$$
g_k(R)
\sim
\frac{e^{i(kR+\pi/4)}}{\sqrt{8\pi kR}},
\qquad
\nabla\log g_k(R)
\sim
\left(ik-\frac{1}{2R}\right)\widehat R.
$$

Born Tx–Rx kernel 的幅度近似

$$
|a(z)|\propto(R_{\mathrm{tx}}R_{\mathrm{rx}})^{-1/2}.
$$

monostatic 且两段距离均约为 $R$ 时，

$$
\boxed{
K_{\mathrm{IS}}=O(R^{-2})
\quad\text{in 2D}.
}
$$

同时 $\|\partial_xg\|/\|g\|\sim k$。无 pose prior 时，projector 对 $B$ 的整体缩放不敏感，retention 主要由 subspace orientation 决定；有 prior 时，$B^*B$ 与 $J_X$ 的尺度比也参与。因此远离场景有时可能提高相对 retention，却显著降低 absolute information。

### 14.4 多频率的共同补偿条件

把各频率堆叠：

$$
A=\begin{bmatrix}A_1\\ \vdots\\ A_F\end{bmatrix},
\qquad
B=\begin{bmatrix}B_1\\ \vdots\\ B_F\end{bmatrix}.
$$

若所有频率共享同一个物理 pose perturbation $z$，则

$$
\boxed{
u\in\ker K_{\mathrm{SLAM}}^{\mathrm{multi}}
\iff
\exists z:\ A_fu=B_fz\quad\forall f.
}
$$

每个频率分别存在不同 $z_f$ 并不足够；必须有同一个 compensation 同时拟合所有频率。这是 multi-frequency 打破 map–pose confounding 的关键机制。

由变分式，

$$
u^*K_{\mathrm{eff}}^{(1:F)}u
=
\min_z
\left[
\sum_{f=1}^F\|A_fu-B_fz\|^2+z^*J_Xz
\right],
$$

加入独立的新频率块只增加非负项，故

$$
\boxed{
K_{\mathrm{eff}}^{(1:F+1)}
\succeq
K_{\mathrm{eff}}^{(1:F)}.
}
$$

**[已推导/需条件]** 这是 consistent stacking、同一参数/pose、正确噪声权重与同一 prior 下的 absolute information 单调性。normalized $\rho_i$ 不一定单调，因为 $K_{\mathrm{IS}}$ 也改变。

若新块只是

$$
A_2=cA_1,\qquad B_2=cB_1,
$$

则 absolute information 因独立重复观测增大，但

$$
\rho_i^{\mathrm{new}}=\rho_i^{\mathrm{old}}.
$$

真正改善分离要改变 $\operatorname{Range}(A_f)$ 与 $\operatorname{Range}(B_f)$ 的联合几何，而不只是复制 SNR。

### 14.5 frame/normal operator 与轨迹采样

把一次 measurement configuration 记为

$$
s=(x,\omega,p,\mathrm{Tx},\mathrm{Rx},\ldots),
$$

并写

$$
y(s)=\langle g_s,J\rangle+\eta.
$$

采样测度 $\mu$ 诱导 frame/normal operator

$$
\boxed{
K_\mu
=
\int g_s\otimes g_s\,d\mu(s).
}
$$

连续轨迹 $x(t)$ 的 empirical measure 与 operator 为

$$
\mu_T=\frac1T\int_0^T\delta_{s(t)}\,dt,
\qquad
K_T=\int_0^Tg_{x(t)}\otimes g_{x(t)}\,dt.
$$

有限 aperture 的 Born normal kernel 可写成

$$
\boxed{
K_\Gamma(r,r')
=
k_b^4
\int_\Gamma
G_b^*(r_s,r)G_b(r_s,r')\,dr_s.
}
$$

在 adjoint/backprojection 意义下，point-spread function 就是这一 kernel 的一列/切片。它同时展示：

- physics-limited information：传播、波长、衰减、边界和 aperture；
- sampling-limited information：实际 trajectory/Tx/Rx/frequency 是否充分覆盖。

更多样的独立 measurement blocks 以 PSD 方式累积 information；“测量数多”若只是高度相关重复，不等于 diversity 高。

### 14.6 分辨率、reach 与 evanescent 限制

讨论区分四类 reach：

1. detection reach；
2. localization reach；
3. imaging/reconstruction reach；
4. joint SLAM reach。

它们不能用一个“最远距离”替代。range/cross-range resolution 的经典带宽/孔径近似只是局部成像尺度，不能证明高维材料场可稳定恢复。

对横向空间频率 $k_\perp$，

$$
k_z=\sqrt{k^2-k_\perp^2}.
$$

当 $k_\perp>k$ 时是 evanescent component，随距离 $z$ 衰减为

$$
\exp\!\left[-z\sqrt{k_\perp^2-k^2}\right].
$$

讨论给出的可检测横向频率粗略上限为

$$
\boxed{
k_{\perp,\max}
\approx
\sqrt{
k^2+
\left(\frac{\log\mathrm{SNR}}{2z}\right)^2
}.
}
$$

**[直觉/需条件]** 常数依 measurement convention 和阈值而变，但核心结论是：远场 subwavelength information 随距离指数丢失，而 SNR 对可恢复高频的改善仅呈 logarithmic scaling。

---

## 15. MIMO、MSR、MUSIC 与分块建模

### 15.1 full-wave MIMO 矩阵模型

对时刻 $t$，令：

- $G_T(x_t)$：Tx array 到 DOI 的 incident-field operator；
- $G_R(x_t)$：DOI current 到 Rx array 的 observation operator；
- $S_t$：orthogonal pilots / transmit code matrix；
- $J_t\in\mathbb C^{N\times N_T}$：每个 Tx/pilot 对应的一列 current；
- $Y_t\in\mathbb C^{N_R\times P}$：MIMO data matrix。

full-wave system 为

$$
\boxed{
J_t
=D_\chi\,[G_T(x_t)S_t+G_DJ_t],
\qquad
Y_t=G_R(x_t)J_t+N_t.
}
$$

令

$$
T_\chi=(I-D_\chi G_D)^{-1}D_\chi,
$$

则

$$
\boxed{
Y_t
=G_R(x_t)T_\chi G_T(x_t)S_t+N_t.
}
$$

若 $S_t=I$，可把

$$
H_t=G_R(x_t)T_\chi G_T(x_t)
$$

视作 full-wave multistatic channel / MSR-like matrix。Born 下 $T_\chi\approx D_\chi$。

vectorized data equation 为

$$
\boxed{
\operatorname{vec}(Y_t)
=(I_P\otimes G_R(x_t))
\operatorname{vec}(J_t)
+\operatorname{vec}(N_t).
}
$$

它明确说明 MIMO 的矩阵行/列结构来自 Rx、Tx/pilot；地图 $\chi$ 是跨通道共享的，但 current columns 随 illumination 而变。

### 15.2 MSR / MUSIC / SOM 的连续关系与差别

在点散射或 Born 近似下，MSR 常写成

$$
K_{\mathrm{MSR}}\approx G_SC G_I,
$$

full-wave 延伸为

$$
\boxed{
K_{\mathrm{MSR}}=G_ST_\chi G_I.
}
$$

MUSIC 对 MSR/data matrix 的 **left/right signal/noise subspaces** 测试 steering vector；SOM 则对已知 $G_S$ 的 **right singular current subspace** 分解 $J$，再借助 state equation 求 complementary current 与 contrast。

可把技术谱系概括为：

$$
\boxed{
\text{MSR/MUSIC data-subspace imaging}
\longrightarrow
\text{SOM current-subspace inversion}
\longrightarrow
\text{full-wave MIMO map–pose tangent analysis}.
}
$$

但这个箭头只表示概念/建模关联，不表示三个方法的 subspaces 可直接相等。

### 15.3 “MIMO 低秩”的严格修正

**[冲突/修正]** MIMO/MSR matrix 在代数上并不必然低秩；它完全可以满秩。在成像问题中常见的是 **low numerical/stable rank**，其来源可包括：

- 有限 aperture 与波长决定的有限 electromagnetic DoF；
- 少量有效散射中心或低维参数化；
- compact radiation operator 的 singular-value decay；
- 通道/角度相关性；
- 噪声阈值使小 singular values 不可用。

因此报告 rank 时要给出是 exact rank、numerical rank、effective rank 还是 physics-based DoF。

### 15.4 block diagonal 不是低秩，也不是 Jordan-like

若按时刻/频率排列 measurement operator 得到

$$
\mathcal G
=\operatorname{blkdiag}(G_1,\ldots,G_T),
$$

则

$$
\operatorname{rank}\mathcal G
=\sum_t\operatorname{rank}G_t,
$$

其 singular values 是各 blocks singular values 的并集。block diagonal 本身不意味着低秩。

**[冲突/修正]** 这种结构也不应称为 “Jordan-like”。Jordan form 是 square non-normal operator 在 similarity 变换下的 canonical form；这里适合的工具是 SVD、Hermitian spectral theorem、Schur complement 或 sparse/block linear algebra。

静态地图 inverse problem 也不能因 measurement blocks 独立就写成 map Hessian 的 block diagonal。共同 $\chi$ 给出

$$
\boxed{
A=\begin{bmatrix}A_1\\ \vdots\\ A_T\end{bmatrix},
\qquad
K_{\mathrm{IS}}=A^*A=\sum_tA_t^*A_t.
}
$$

pose Jacobian 在独立 pose 情形可以 block structured；odometry/IMU prior 则使 pose information block-banded。两者经 Schur complement 后都会耦合共同地图变量。

### 15.5 三层“谱”不能混用

讨论最终区分：

1. **radiation observability**：$G_S$ 的 SVD，回答 current 如何辐射到 sensors；
2. **state relevance**：$G_D$ 或相关 state/domain operator 的谱，回答 complementary current 如何受内部物理约束；
3. **pose distinguishability**：$A$ 与 $B$ 在 data space 的几何，回答 map tangent 是否能被 pose tangent 模拟。

GS/GD/pose 三层可以组成算法 pipeline，但其 eigenvalues、thresholds、spaces 和 dimensions 不能直接复用。

---

## 16. FFT-TSOM / NFFT-SOM 与学习模块的定位

### 16.1 讨论中保存下来的算法直觉

原讨论把 FFT-TSOM 理解为：用低频 Fourier basis 近似 $V_D^+$ 的 span，从而避免显式构造/分解大型 domain operator。这个近似是 **subspace-level approximation**，不是要求 Fourier vector 与每个 singular vector 一一对应。

所谓 NFFT-SOM 在这些讨论里指 **New Fast Fourier Transform SOM**，不是通常数值分析中的 non-uniform FFT。其 working form 被概括为

$$
\boxed{
J=J^++F\alpha,
}
$$

其中 $J^+$ 是从 dominant/data-visible part 得到的初始 current，Fourier correction $F\alpha$ 可补偿 $J^+$ 的 bias，并可 coarse-to-fine 地扩张 nested basis。

### 16.2 与 pose defect 的可能组合

**[直觉/假设]** 可先由 TSOM/FFT basis 得到 $U_T$，再只在

$$
\delta J=U_T\delta c
$$

中构造 pose defect $C_P$。这样：

- 大型 current-space reduction 由 SOM/FFT 完成；
- pose 影响由最多 $p$ 维 low-rank update 表达；
- rank/event tracking 只更新小 defect subspace；
- 学习模块若使用，应补足被物理 reduction 留下的 residual，而不是伪装成新增测量信息。

**[待验证]** Fourier approximation 对 full-wave $V_D^+$、不同频率和移动 geometry 的误差，需要 subspace-angle/operator-norm bound；七条讨论没有给出这一证明。

### 16.3 为什么没有把 ML 放入核心理论

讨论中考虑过让 ML 预测：

- $V_S^+\cap V_D^-$ 或其他难解 current components；
- pose-defect correction；
- regularized current/map update。

最终 scope 决定把 ML 从核心理论稿移出，原因不是 ML 无用，而是当前可证明对象已经足够清晰：

$$
\text{full-wave Jacobian}
\rightarrow
\text{Schur loss}
\rightarrow
\text{retention spectrum}
\rightarrow
\text{rank-changing trajectory}.
$$

在没有明确训练目标、equivariance、uncertainty calibration 和 data-consistency guarantee 前，ML 只能列为后续算法扩展。

---

## 17. 轨迹设计的指标层次

### 17.1 已知 pose 时的 classical criteria

对 $K_{\mathrm{IS}}(X)$ 或适当 posterior information：

- E-optimal：最大化 $\lambda_{\min}$ 或任务内第 $r$ 个 eigenvalue；
- D-optimal：最大化 $\log\det$ / pseudo-determinant；
- A-optimal：最小化 $\operatorname{tr}K^\dagger$ 或 posterior covariance trace；
- frame tightness：压缩最大/最小 eigenvalue 差，追求较均匀的 observable directions。

### 17.2 SLAM 情形必须把 nuisance loss 算进去

只优化 $K_{\mathrm{IS}}$ 可能选择“信号很强但 map signature 与 pose signature 几乎平行”的轨迹。讨论建议用

$$
K_{\mathrm{eff}}(X)
=
K_{\mathrm{IS}}(X)-L_X(X)
$$

直接评分，并同时报告：

- absolute spectrum $\lambda_i(K_{\mathrm{eff}})$；
- relative retention $\rho_i$；
- defect rank / $d_{\mathrm{conf}}$；
- execution sensitivity $\|DK_{\mathrm{eff}}\|$ 或 $\|\nabla\lambda_r\|_*$。

一个完整的 robust E-type objective 是

$$
\boxed{
\max_{X\in\mathcal X}
\min_{\|\Delta X\|\le\varepsilon}
\lambda_r(K_{\mathrm{eff}}(X+\Delta X)).
}
$$

其一阶 surrogate 为

$$
\lambda_r(K_{\mathrm{eff}}(X))
-
\varepsilon\|\nabla_X\lambda_r\|_*.
$$

### 17.3 不能把 relative retention 当唯一目标

极弱的 $A$ 可能因几乎不与 $B$ 重叠而得到 $\rho\approx1$，但 absolute information 仍接近零。因此 trajectory objective 至少要同时控制：

$$
\boxed{
\text{absolute map strength}
\quad+\quad
\text{map–pose separation}
\quad+\quad
\text{execution robustness}.
}
$$

这也是为何“远离目标可能更 pose-robust”并不意味着远距离轨迹更适合成像。

---

## 18. 核心定理表、连续指标与条件审计

### 18.1 有限维主结论

在 whitened、realified、局部线性模型

$$
\delta y=A\delta\chi+B\delta X+n
$$

下，令

$$
\mathcal U=\operatorname{Range}(A),
\qquad
\mathcal V=\operatorname{Range}(B).
$$

七条讨论中已经闭合的有限维结论为：

$$
\boxed{
\ker K_{\mathrm{SLAM}}
=\{u:Au\in\mathcal V\},
}
$$

$$
\boxed{
\operatorname{rank}K_{\mathrm{IS}}
-
\operatorname{rank}K_{\mathrm{SLAM}}
=
\dim(\mathcal U\cap\mathcal V),
}
$$

$$
\boxed{
\mathcal U\cap\mathcal V=\{0\}
\iff
\operatorname{rank}[A\ B]
=\operatorname{rank}A+\operatorname{rank}B,
}
$$

$$
\boxed{
\dim(\mathcal U\cap\mathcal V)
\ge
\max(0,r_A+r_B-m),
}
$$

其中 $m$ 是 realified data dimension。

这给出一个纯维数判据：若 $r_A+r_B>m$，exact map–pose confusion 不可避免；物理 gauge 还可能使交集大于 generic subspaces 的维数下界。

### 18.2 三种不同的“混淆维数”

无 pose prior 时：

$$
\boxed{
\begin{aligned}
\text{完全丢失的 DoF}
&=\#\{i:\theta_i=0\},\\
\text{受到任何 pose 影响的 modes}
&=\#\{i:\theta_i<\pi/2\},\\
\text{连续混淆维数}
&=\sum_i\cos^2\theta_i.
\end{aligned}
}
$$

令 $Q_A,Q_B$ 分别为 $\mathcal U,\mathcal V$ 的正交 basis，则

$$
\boxed{
d_{\mathrm{conf}}
=\sum_i\cos^2\theta_i
=\|Q_B^*Q_A\|_F^2
=\operatorname{tr}(I-R_{\mathrm{op}}).
}
$$

这三个量分别回答“完全毁掉多少”“触及多少”“总体混淆强度多大”，不能互相替换。

### 18.3 stable transversality 与 volume retention

若

$$
\rho_{\min}=\sin^2\theta_{\min}>0,
$$

则

$$
\boxed{
\rho_{\min}K_{\mathrm{IS}}
\preceq
K_{\mathrm{SLAM}}
\preceq
\rho_{\max}K_{\mathrm{IS}}
}
$$

在 observable support 上成立。

无 exact confusion 时，

$$
\boxed{
\frac{\det_{\mathcal S}K_{\mathrm{SLAM}}}
{\det_{\mathcal S}K_{\mathrm{IS}}}
=
\prod_{i=1}^r\rho_i
=
\prod_{i=1}^r\sin^2\theta_i.
}
$$

故 $\sum_i\log\rho_i$ 是自然的 log-volume retention criterion；只要存在一个 $\theta_i=0$，volume ratio 就为零。

### 18.4 无限维版本的额外警告

有限维中 $\mathcal U\cap\mathcal V=\{0\}$ 可推出正最小夹角；无限维 compact inverse problem 中不一定，因为 $\operatorname{Range}(A)$ 可能不闭。可能有零交集但存在序列逼近 pose range。稳定条件应写成

$$
\boxed{
\inf_{Au\ne0}
\frac{\operatorname{dist}(Au,\mathcal V)}
{\|Au\|}
>0.
}
$$

此外 $K_{\mathrm{IS}}^{\dagger/2}$ 可能无界。严格 operator 版本可取 polar decomposition

$$
A=U_A|A|
$$

并在 $\overline{\operatorname{Range}(A^*)}$ 上定义

$$
\boxed{
R_{\mathrm{geom}}
=U_A^*P_{B^\perp}U_A.
}
$$

有限维时它与 $K_{\mathrm{IS}}^{\dagger/2}K_{\mathrm{SLAM}}K_{\mathrm{IS}}^{\dagger/2}$ 等价。

### 18.5 公式使用前的条件检查表

| 对象/结论 | 状态 | 必须检查 |
|---|---|---|
| full-wave $A=G_SM^{-1}D_E$ | **[已推导/需条件]** | $M$ 可逆；world-fixed map；map parameterization 明确 |
| 简化 pose $Bh=(DG_S[h])j+G_SM^{-1}D_\chi De^{\mathrm{inc}}[h]$ | **[需补条件]** | $DG_D=0$、$D_XD_\chi=0$、无 direct/calibration path |
| $K_{\mathrm{IS}}=A^*A$ | **[已接受]** | 先 whitening/realification；否则带正确 metric |
| $K_{\mathrm{SLAM}}=A^*P_{B^\perp}A$ | **[已推导/已接受]** | deterministic free nuisance；closed range/finite dimension |
| generalized Schur $K_{\mathrm{eff}}$ | **[已推导/已接受]** | $J_X\succeq0$；一致的 support/range 处理 |
| $J_X\to0\Rightarrow K_{\mathrm{eff}}\to K_{\mathrm{SLAM}}$ | **[需补条件]** | 沿 fixed-rank regular path，例如 $\epsilon I$ |
| $J_X\to\infty\Rightarrow K_{\mathrm{eff}}\to K_{\mathrm{IS}}$ | **[需补条件]** | prior 在所有 data-coupled pose directions 上变强 |
| $\rho_i=\sin^2\theta_i$ | **[已推导/已接受]** | 无 prior projector case；observable support |
| finite-prior $\rho_i$ | **[冲突/修正]** | weighted shrinkage；不再是普通 principal angles |
| rank loss/interlacing | **[已推导/已接受]** | finite dimensional/finite-rank perturbation |
| eigenvector derivative | **[需补条件]** | simple eigenvalue、positive gap、phase gauge |
| projector derivative | **[需补条件]** | $B$ locally constant rank |
| robust first-order objective | **[需补条件]** | smoothness、gap、correct uncertainty set/dual norm |
| hard GS/GD/pose categories | **[直觉/假设]** | threshold、gap、commutativity 与 stability 必须报告 |

---

## 19. 最小二维验证实验（原讨论建议）

目标不是先构建完整 online SLAM，也不是训练 network，而是验证现有恒等式、gauge、频率和轨迹规律。

### 19.1 场景与离散

- $D=[-0.5,0.5]^2$；
- $32\times32$ 或 $40\times40$ grid；
- 非零名义 contrast，例如两个平滑介质体：

$$
\chi_0
=0.3\,\chi_{\mathrm{disk},1}
+0.5\,\chi_{\mathrm{disk},2}.
$$

验证 gauge 时宜用 Fourier/spline/光滑 FEM basis；piecewise-constant pixels 的 infinitesimal translation 不一定精确留在同一离散空间。

### 19.2 sensor 与 trajectory

- 每个 pose：一个 Tx，$4$–$8$ 个固定在机器人载体上的 Rx；
- $x_t=(p_x,p_y,\theta)$；
- 对比 straight line、$90^\circ$ arc、$180^\circ$ arc、$360^\circ$ circle；
- 固定总路径长度和 measurement count，避免把“更多数据”误判为“更好 geometry”。

### 19.3 验证顺序

1. 解 $j_t,E_t^{\mathrm{tot}}$。
2. 用解析式构造 $A,B$。
3. centered finite difference 检查 $A\delta\chi$ 与 $B\delta X$。
4. whitening、realification。
5. 计算 $K_{\mathrm{IS}},K_{\mathrm{SLAM}},K_{\mathrm{eff}}$。
6. 对 $A,B$ 作 rank-revealing QR/SVD，构造 $Q_A,Q_B,Z$。
7. 验证

   $$
   \operatorname{eig}(R_{\mathrm{op}})
   =1-\sigma_i^2(Q_B^*Q_A).
   $$

8. 验证 exact rank identity、loss rank bound 与 interlacing。
9. 构造 $SE(2)$ map/pose generators，检查

   $$
   \frac{\|A\delta\chi+B\delta X\|}
   {\|A\delta\chi\|+\|B\delta X\|}
   \approx0.
   $$

10. 扫描

   $$
    J_X=\alpha I,\qquad \alpha=10^{-6},\ldots,10^4,
   $$

    验证

   $$
    K_{\mathrm{SLAM}}
    \preceq K_{\mathrm{eff}}(\alpha)
    \preceq K_{\mathrm{IS}}.
   $$

11. 对 trajectory 作小扰动，比较实际 $\delta\lambda_i$ 与 $v_i^*DK[\Delta X]v_i$。
12. 在 $\|\Delta X\|=\varepsilon$ 边界采样，比较真实 worst case 与一阶 robust surrogate。
13. 用单频/多频开关验证“共同 $z$”条件，并设置纯重复 block 作为 retention 不变对照。
14. 用 Born/full-wave 开关观察 $M^{-1}$ 对 information 与 sensitivity 的同时放大。

数值 rank tolerance 应取机器精度型标准，例如

$$
\tau_{\mathrm{rank}}
=\max(m,n)\epsilon_{\mathrm{mach}}\sigma_1,
$$

不能用 SOM 的物理 dominant-subspace threshold 来“验证”exact rank theorem。

### 19.4 最小报告内容

每个 geometry 至少报告：

- $A,B$ derivative check relative error；
- realified ranks $r_A,r_B,\operatorname{rank}[A\ B]$；
- $\lambda_i(K_{\mathrm{IS}})$ 与 $\lambda_i(K_{\mathrm{eff}})$；
- $\rho_i$、principal angles、$d_{\mathrm{conf}}$；
- exact/near gauge residual；
- prior sweep；
- trajectory perturbation error vs first-order prediction；
- Born/full-wave 与 single/multi-frequency 的 controlled comparison。

---

## 20. 文献映射与 novelty 边界（只保存原讨论判断）

> 本节是七条讨论中的 literature mapping，不是本次整理重新完成的系统检索。时间敏感的“最新工作”、年份和 novelty 均须在投稿前独立核验。

### 20.1 已知成熟组件

以下单独看都不宜宣称新颖：

- Chen SOM/TSOM 的 dominant/complementary current decomposition；
- MSR/MUSIC 的 signal/noise subspace；
- SVD、CS decomposition、principal angles；
- nuisance elimination 的 Schur complement / equivalent Fisher information；
- gauge 与 identifiability-up-to-group-action；
- finite-rank Hermitian perturbation 与 interlacing；
- Grassmann/projector perturbation；
- active sensing 中的 A/D/E-optimal criteria。

### 20.2 原讨论认为最接近的文献板块

- **Chen SOM / dominant-current / SOM-Net**：研究 $G_S$ 对 current 的压缩以及 complementary current/contrast reconstruction；不直接处理 pose nuisance retention。
- **full-wave inverse-scattering DoF**：以 radiation 或 contrast-modified radiation matrix 的 significant singular values 表达物理 DoF；通常假设 geometry 已知。
- **EFIM / Schur complement**：成熟地消去 channel、clock 或其他 nuisance；说明代数本身不是 novelty。
- **RF-SLAM PCRB / snapshot identifiability / double-bounce SLAM**：处理有限维 landmark/surface/path parameters 与 pose 的联合信息，概念上很近，但地图通常不是 volumetric Helmholtz contrast。
- **coherent direct multipath SLAM**：直接从 coherent RF signal 联合估计而非先提 path features，更接近 wave-level likelihood，但其表示与本文 full-wave tangent spectrum 不同。
- **bilinear inverse problems / blind calibration**：研究 transformation group、gauge 与 identifiability；通常不含 Helmholtz resolvent、trajectory aperture 与 volumetric map。
- **active SLAM / robust Bayesian OED**：会用 information/covariance 选轨迹；本文拟议指标特别把 map strength、pose confounding 和 execution sensitivity 合在一起。

### 20.3 可能的组合性贡献，尚未确立 novelty

原讨论的谨慎表述是：

$$
\boxed{
\text{full-wave physics}
+
\text{nuisance tangent geometry}
+
\text{spectral/trajectory consequences}
}
$$

这一组合可能有研究空间，尤其是：

1. full-wave TSOM-recoverable subspace 内的 pose-induced relative spectral defect；
2. $SE(2)/SE(3)$ gauge 与 $\rho=0$ 的显式联系；
3. Born empty-background degeneration 与二阶 bilinear mismatch；
4. 多频率 shared-compensation transversality；
5. rank-adaptive defect tracking 与 robust trajectory bounds。

**[待验证/开放问题]** 七条讨论没有完成 exhaustive prior-art search，因此不能写“此前无人研究”；最多可说原讨论当时未找到把上述全部组件系统结合的同构工作。

---

## 21. 冲突、修正与尚未接受的说法总表

| 议题 | 早期/冲突说法 | 整理后的状态 |
|---|---|---|
| $A/B$ 命名 | S5 早期用 $A$=pose、$B$=map | **[冲突/修正]** 本文固定 $A=D_\chi F$、$B=D_XF$ |
| $\rho$ 含义 | 早期为 $\cos\theta$ confounding amplitude | **[冲突/修正]** 后续统一为 $\sin^2\theta$ retention |
| $K_{\mathrm{SOM}}$ | 曾表示 raw map information，也曾泛指 pose 后信息 | **[冲突/修正]** 用 $K_{\mathrm{IS}},K_{\mathrm{SLAM}},K_{\mathrm{eff}}$ 分开 |
| $K_{\mathrm{eff}}$ 与 map prior | 曾把 posterior curvature 也叫 $K_{\mathrm{eff}}$ | **[冲突/修正]** map prior 后记 $K_{\mathrm{post}}$ |
| complex projector | 直接用 complex $BB^\dagger$ | **[冲突/修正]** 实 pose 必须先 realify |
| Born $\chi_0=0$ | $B=0$ 被误读为 pose 无害 | **[冲突/修正]** 影响进入 $(D_XA[\delta X])\delta\chi$ 双线性项 |
| $\rho_i$ 对应 modes | 当成普通 $K_{\mathrm{IS}}$ eigenmode 的逐项比值 | **[冲突/修正]** 是 Fisher-whitened principal directions |
| exact gauge | $\rho=0$ 与 global symmetry 等价 | **[冲突/修正]** gauge $\Rightarrow\rho=0$，反向不必成立 |
| $\rho=1$ | 被说成“不是 SLAM mode” | **[冲突/修正]** 只表示当前线性化下一阶 pose-decoupled |
| finite pose prior | 仍套 $\sin^2\theta_i$ | **[冲突/修正]** 一般是 weighted shrinkage |
| Two-fold exact formula | 不同 $V_D^\pm$ 定义被混成一个版本 | **[待验证]** 必须回到具体 Chen 文献/离散约定 |
| 三个投影八格分解 | 默认 $S,D,P$ 可交换 | **[冲突/修正]** 一般不交换，只能 hierarchical 或 approximate |
| Grassmann continuity | $0<\rho<1$ 即保证连续 | **[冲突/修正]** 需要 fixed rank + external spectral gap |
| Flag manifold | 被设想为自动处理 rank change | **[直觉/假设]** 最多作 transition bookkeeping |
| low-rank loss | 被理解为只改变 $p$ 个 ordinary eigenvalues | **[冲突/修正]** generalized defect 至多 $p$ 维；ordinary eigenpairs 可普遍移动 |
| regularization | 被概括为“本质低秩” | **[冲突/修正]** 仅某些形式低秩/effectively low-rank，不创造 data information |
| pose perturbation $\Delta K$ | 因参数少而被称为低秩矩阵 | **[冲突/修正]** 是低维结构化矩阵族，成员可满秩 |
| MIMO matrix | 被视为必然低秩 | **[冲突/修正]** 仅常有 low numerical/stable rank |
| block diagonal | 被称作 low-rank / Jordan-like | **[冲突/修正]** rank 相加；应使用 SVD/block algebra，不是 Jordan 语义 |
| inner/outer pose uncertainty | $J_X$ 与 $\Delta X$ 混写 | **[冲突/修正]** nuisance prior 与 execution mismatch 分开 |
| ML subspace | $V_S^+\cap V_D^-$ 被直接指定给 ML | **[直觉/假设]** 不是原 TSOM 定义，已移出核心 scope |
| 文献“最新/首次” | 讨论中的搜索结果被当成确定 novelty | **[待验证]** 投稿前重做系统查重 |

---

## 22. 待验证问题与研究路线

### 22.1 理论优先级

1. **有限维主定理稿**：把 kernel、rank、principal angles、continuous confounding、determinant、low-rank loss、interlacing 放在统一假设下完整证明。
2. **连续 Helmholtz gauge theorem**：严格定义 $SE(2)/SE(3)$ 对 map/trajectory 的群作用，并说明 boundary、anchor、known background 如何破坏 gauge。
3. **Born far-field trajectory theorem**：对 line/circle/partial aperture 研究 Fourier coverage 与 $A\delta\chi=B\delta X$ 的显式条件。
4. **prior-weighted spectrum**：研究

   $$
   Q_A^*(I+BJ_X^{-1}B^*)^{-1}Q_A
   $$

   的谱界、单调性和与普通 principal angles 的关系。
5. **full-wave resolvent perturbation**：把 $\|DA\|,\|DB\|,\|DK_{\mathrm{eff}}\|$ 显式绑定到 $\|M^{-1}\|,k,R,\|\chi\|$ 与 linearization radius。
6. **multi-frequency transversality**：证明何种角度/带宽条件能除 global gauge 外 generically 分离 map 与 pose。
7. **rank-event theory**：给出 defect projector 在 fixed-rank strata 上的 smoothness 以及 crossing 时的更新规则。
8. **最小数值验证**：先验证理论恒等式，再决定是否扩展 reconstruction/online SLAM/learning。

### 22.2 开放的具体猜想

- 全角、多频、足够 MIMO diversity 时，除 global $SE(2)/SE(3)$ gauge 外，map/pose tangents 是否 generically transverse？
- 单频、有限 aperture 下，额外 near-confounded modes 能否由 missing wedge 明确刻画？
- bandwidth 是否给出 $\theta_{\min}$ 的非平凡下界？
- 接近 full-wave resonance 时，absolute information 增益能否抵消 $DK_{\mathrm{eff}}$ 与 model mismatch 的增长？
- $V_D^+$ 的 Fourier surrogate 在什么波长/contrast/geometry regime 保持小 Grassmann distance？
- motion prior graph 的 global gauge 与 map gauge 如何共同进入 generalized Schur support？
- 若 calibration、clock、phase center、channel gains 同时未知，defect rank/structure 如何分层？
- nonlinear/global ambiguity、cycle skipping 与一阶 retention spectrum 之间可否建立 basin-of-attraction bound？

### 22.3 最终收敛的核心研究问题

原讨论最终把核心问题收窄为：

> **在 full-wave inverse-scattering SLAM 中，有限维位姿不确定性如何形变 TSOM 可恢复子空间；这种形变能否表示成低秩 generalized-eigenvalue defect，并在轨迹引发的 rank-changing events 中被稳定追踪？**

等价的几何表述是：

$$
\boxed{
\text{地图流形与位姿流形在数据空间中的切空间有多接近？}
}
$$

其中：

$$
\operatorname{Range}(A)\cap\operatorname{Range}(B)
$$

描述 exact local ambiguity；

$$
\sin^2\theta_i
$$

描述 continuous relative retention；

$$
A^*P_{B^\perp}A
$$

把几何变成 absolute effective information；

$$
DK_{\mathrm{eff}}
$$

把它连接到 trajectory sensitivity 与 robust reach。

### 22.4 核心论文暂时排除的内容

为保持可证明性，讨论建议首篇理论稿暂不把以下项目设为主贡献：

- 完整 NFFT/FFT-SOM 工程优化；
- 机器学习重建；
- MIMO waveform/array co-design；
- 非交换投影的八格分类；
- 大型 online active-SLAM 系统；
- 未完成查重的“首次”宣称。

它们仍是后续方向，并未被判定为无价值。

---

## 23. 七个来源的主题索引

| 来源 | 本文吸收的主要 context | 主要状态/遗留 |
|---|---|---|
| S1 全波逆散射SLAM理论研究 | TSOM reduced basis；pose low-rank defect；small $p\times p$ solve；LM/Woodbury；incremental basis；rank events；NFFT-SOM 联系 | 数值算法与 Fourier surrogate 仍待验证 |
| S2 总结SOM逆散射SLAM | Chen SOM/Improved/Two-fold/FFT/NFFT 谱系；GS/GD/pose hierarchy；TriSpace；MIMO/MSR/MUSIC；noncommuting projector 与 scope 收敛 | Two-fold 精确版本需对原论文；novelty 未确立 |
| S3 推导位姿不确定性下信息损失 | 严格 $A,B$；realification；kernel/rank/principal-angle theorems；prior；gauge；Born/far-field；multi-frequency；$DK_{\mathrm{eff}}$；robust bound；最小实验 | full local theory 最完整；连续/无限维与实验尚待做 |
| S4 Regularization低秩直觉 | affected vs destroyed；interlacing；prior shrinkage；trace/operator-norm loss；regularization 不创造信息；structured perturbation 区分 | growing pose dimension 与 effective rank 需谨慎 |
| S5 迁移至SLAM | 从 inverse scattering observability 到 joint map–pose tangent；L2 compensation；符号/语义审计；三种 pose models 的起点 | 早期 $A/B,\rho,K$ 命名已在本文修正 |
| S6 推导谱可观测性 | Green-function 字典；absolute vs relative axes；generalized spectrum；trajectory derivative；reach/physics-vs-sampling；研究问题细化 | 部分长回复的接口预览存在截断，后续轮次与其他来源覆盖其主结论 |
| S7 逆散射与SLAM技术 | inverse scattering/SLAM/Radio/Radar 边界；统一 wave-MAP；uniqueness/observability/stability/CRB；frame operator、PSF、evanescent reach、trajectory OED | 文献例子与定量常数需按具体系统再核验 |

### 23.1 读取与完整性说明

- 七个 conversation 均已定位到唯一 ID 并读取到其 conversation page 的末尾；它们不是按标题搜索后拼接的二手摘要。
- S3 的单条超长回答在当前项目已有完整本地副本 Theory/Questions/A1.md，本文据该完整版本提取。
- S2、S7 的长回复通过原 conversation 页面补足了接口预览截断的部分。
- S6 有三条超长历史回复在标准接口预览中出现单条长度截断；其可见正文、后续追问修正及与 S1–S5/S7 重叠的推导均已纳入。无法确认的截断尾句没有被伪装成已读取结论。
- 原讨论中的网页引用和“最新文献”没有在本次整理中重新联网核验；本文保存的是讨论 context，而非重新做 literature review。

---

## 24. 快速引用页

### 24.1 五个核心 operator

$$
\boxed{
\begin{aligned}
A&=D_\chi F,\\
B&=D_XF,\\
K_{\mathrm{IS}}&=A^*A,\\
K_{\mathrm{SLAM}}&=A^*P_{B^\perp}A,\\
K_{\mathrm{eff}}
&=A^*A-A^*B(B^*B+J_X)^\dagger B^*A.
\end{aligned}
}
$$

### 24.2 四个核心谱关系

$$
\boxed{
\begin{aligned}
K_{\mathrm{SLAM}}v_i&=\rho_iK_{\mathrm{IS}}v_i,\\
\rho_i&=\sin^2\theta_i\quad(J_X=0),\\
\operatorname{rank}(K_{\mathrm{IS}}-K_{\mathrm{eff}})
&\le\operatorname{rank}B,\\
\lambda_{i+p}(K_{\mathrm{IS}})
&\le\lambda_i(K_{\mathrm{eff}})
\le\lambda_i(K_{\mathrm{IS}}).
\end{aligned}
}
$$

### 24.3 四条不可忘记的解释

1. $\rho_i$ 是 Fisher-whitened principal directions 的 retention，不是普通 eigenvalue 的逐项比值。
2. $\rho_i\approx1$ 不代表 absolute information 强；必须与 $\lambda_i(K_{\mathrm{IS}})$ 一起看。
3. exact gauge 必然造成 $\rho=0$，但局部 $\rho=0$ 不必来自 global gauge。
4. pose loss 可低秩，而 inverse scattering 的原始 ill-posedness 和 regularization 仍是另一层问题。

### 24.4 一句话主线

$$
\boxed{
\begin{gathered}
\text{SOM/TSOM 先压缩 current ambiguity；}\\
\text{SLAM Schur defect 再刻画保留空间内的 pose ambiguity。}
\end{gathered}
}
$$
