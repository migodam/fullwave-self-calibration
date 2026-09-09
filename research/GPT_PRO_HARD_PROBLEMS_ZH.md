# 给网页版 GPT Pro 的中文总 Prompt：全波逆散射 SLAM 的困难理论、证明与设计问题

> 使用方式：将本文件全文复制给网页版 GPT Pro，并同时上传下列材料中可用的最新版本：
>
> 1. `Theory/SOM_SLAM_THEORY_CONTEXT.md`
> 2. `Theory/Questions/Q1.md`
> 3. `Theory/Questions/A1.md`
> 4. `research/literature/TARGETED_PRIOR_ART_ADDENDUM.md`
> 5. `research/delegated/scholarqa_prior_art/claim_ledger.md`
> 6. 最新论文稿 `experiments/idea_loops/loop_2026-09-03_16-25-02/experiment_pose_confounding_spectral_geometry/manuscript/main_round3.pdf`
> 7. `experiments/idea_loops/loop_2026-09-03_16-25-02/experiment_pose_confounding_spectral_geometry/notes/family5_parent_generalized.md`
> 8. `experiments/idea_loops/loop_2026-09-03_16-25-02/experiment_pose_confounding_spectral_geometry/notes/family15b_born_pose_control.md`
> 9. `research/FINAL_SCIENTIFIC_REVIEW.md`

---

你现在是我的高强度数学科研合作者。请对“位姿不确定性下全波逆散射 SLAM 的谱可观测性”进行一次**理论优先、反例优先、语义严格**的深度研究。你的任务不是重新做已经完成的有限维数值实验，也不是写泛泛的综述，而是解决或严格缩小下面尚未解决的理论、证明、构造和设计问题，最终给出可以直接进入论文定理/附录的数学材料。

如果某个目标在现有假设下不成立，请优先给出最小反例、指出缺失假设并提出可证明的最强修正版；不要为了得到漂亮公式而强行证明错误命题。对无法完成的证明，必须给出明确卡点、可行的下一步和所需引理。不要把数值证据当作证明，不要把有限维结论自动推广到连续模型。

## 一、问题设定与已锁定的有限维事实

从二维频域标量 Helmholtz contrast-source 模型出发。对第 $t$ 个收发配置，

\[
j_t=D_\chi\bigl(e_t^{\mathrm{inc}}+G_{D,t}j_t\bigr),
\qquad y_t=G_{S,t}j_t+n_t,
\]

\[
M_t=I-D_\chi G_{D,t},\qquad
M_tj_t=D_\chi e_t^{\mathrm{inc}},\qquad
E_t^{\mathrm{tot}}=e_t^{\mathrm{inc}}+G_{D,t}j_t.
\]

在 $M_t$ 可逆、名义状态 $(\chi_0,X_0)$ 附近，完整一阶变分为

\[
\delta y=A\,\delta\chi+B\,\delta X+n.
\]

对 world-fixed map、固定 $G_D$ 的基本情形，

\[
A_t=G_{S,t}M_t^{-1}D_{E_t^{\mathrm{tot}}}S_\chi,
\]

而一般位姿方向 $h_t$ 的 $B_t h_t$ 至少包含

\[
DG_{S,t}[h_t]j_t+
G_{S,t}M_t^{-1}D_\chi
\left(De_t^{\mathrm{inc}}[h_t]+DG_{D,t}[h_t]j_t\right),
\]

移动坐标网格、direct path、天线相位中心和标定参数还会增加项。

所有物理参数均按实参数处理：先用噪声协方差白化复数数据，再 realify。后文的 $A,B$ 默认是 whitened、realified 的实矩阵/实算子。禁止直接在复线性位姿空间里用 $BB^\dagger$，因为那会允许不存在的复位姿扰动。

已知位姿、自由位姿和带位姿先验的信息算子分别为

\[
K_{\mathrm{IS}}=A^TA,
\]

\[
K_{\mathrm{SLAM}}=A^T(I-BB^\dagger)A,
\]

\[
K_{\mathrm{eff}}
=A^TA-A^TB(B^TB+J_X)^\dagger B^TA.
\]

以下有限维结论已经有证明并通过数值代数验证，可作为前提，但请检查其连续版本所需条件：

1. $0\preceq K_{\mathrm{SLAM}}\preceq K_{\mathrm{eff}}\preceq K_{\mathrm{IS}}$（中间次序需要对应的 $J_X\succeq0$ 条件）。
2. 
   \[
   \ker K_{\mathrm{SLAM}}
   =\{u:Au\in\operatorname{Ran}(B)\}.
   \]
3. 
   \[
   \operatorname{rank}K_{\mathrm{IS}}-
   \operatorname{rank}K_{\mathrm{SLAM}}
   =\dim(\operatorname{Ran}(A)\cap\operatorname{Ran}(B)).
   \]
4. 在 $K_{\mathrm{IS}}$ 的可观测支撑上，无先验广义保留谱满足
   \[
   K_{\mathrm{SLAM}}v_i=\rho_iK_{\mathrm{IS}}v_i,
   \qquad \rho_i=\sin^2\theta_i,
   \]
   其中 $\theta_i$ 是 $\operatorname{Ran}(A)$ 与
   $\operatorname{Ran}(B)$ 的主夹角。它不是普通
   $K_{\mathrm{IS}}$ eigenmode 的逐项保留率。
5. 信息损失 $K_{\mathrm{IS}}-K_{\mathrm{eff}}$ 的秩不超过
   $\operatorname{rank}(B)$。有限先验的保留谱是 weighted shrinkage，通常不再等于普通 $\sin^2\theta_i$。
6. 多频堆叠会使**绝对**有效信息按 Loewner 次序不减，但归一化保留谱不必单调。完全重复的数据块在无先验时不改变保留谱；有限先验下只有相应共同缩放先验才保持该不变性。
7. Born 模型在 $\chi_0=0$ 时有 $B=0$，所以一阶
   $K_{\mathrm{SLAM}}=K_{\mathrm{IS}}$；位姿影响并未消失，而是以
   $(D_XA[\delta X])\delta\chi$ 的双线性二阶项首先出现。
8. 固定秩区间内，令 $P_B=BB^\dagger$，
   \[
   \dot P_B=(I-P_B)\dot B B^\dagger+(B^\dagger)^T\dot B^T(I-P_B).
   \]
   对简单广义特征值且 $v^TK_{\mathrm{IS}}v=1$，
   \[
   \dot\rho=v^T(\dot K_{\mathrm{eff}}-\rho\dot K_{\mathrm{IS}})v.
   \]
   重根必须用压缩算子；秩事件处不能沿用固定秩导数。
9. 对有限维线性 Born 模型 $F_B(\chi,X)=A_0(X)\chi$ 和
   $\chi=s\bar\chi$，有 $A_B(s)=A_0$、$B_B(s)=sB_1$。因此对每个
   $s\ne0$，无先验投影 $P_{\operatorname{Ran}(B_B(s))}$ 与 Born
   $K_{\mathrm{SLAM}}$ 对振幅 $s$ 不变；但在 $s=0$ 有 $B_B=0$、
   $K_{\mathrm{SLAM}}=K_{\mathrm{IS}}$，除非两数据切向正交，否则这是
   一个 rank-stratum 奇异极限。若 $J_X=\alpha I$ 且 $\alpha>0$，则信息损失
   连续且满足
   \[
   \|K_{\mathrm{IS}}-K_{\mathrm{eff}}(s)\|_2
   \le \frac{s^2}{\alpha}\|A_0^TB_1\|_2^2.
   \]
   该结论不自动覆盖半正定先验，也不表示信号趋零时仍有可用位姿信息。

已经执行的离散实验只说明：上述代数在特定 2D scalar、有限网格、局部 Jacobian 模型中自洽。特别是：广义导数在正谱隙模式上呈二阶有限差分误差；一个人为构造的 $\operatorname{rank}(B):18\to17$ 事件使投影算子跳变范数为 1；多条轨迹的排序会随指标、场景、standoff 和归一化改变。它们都不是连续定理或真实系统验证。

## 二、绝对禁止的语义漂移

请全程遵守：

1. Chen SOM 的 $G_S$ 奇异子空间是 **current-to-data** 空间；本研究的 $A=G_ST_\chi$ 是 **map-tangent-to-data** 算子。二者相关但不相同，除非显式证明 pullback/composition 关系。
2. contrast source $j$ 与 pose-prior information $J_X$ 不得混淆。
3. 内层 nuisance pose $\delta X$ 与外层轨迹执行误差 $\Delta X$ 是不同对象，不得重复计数。
4. local Fisher/Schur 结论不等于全局唯一性、非线性优化收敛、无 cycle skipping、在线 SLAM 成功或真实硬件性能。
5. 不能使用“首次”“无已有工作”等全局新颖性语言。现有检索受到 Semantic Scholar 429 限流，只能做 retrieval-bounded 的条件性定位。已知 Karthik 与 Ghosh（PIERS 2023，DOI 10.1109/PIERS59004.2023.10221374）已经研究过 contrast reconstruction 与 transmitter localization，因此不能声称一般性的“首次联合逆散射与定位”。
6. 不得把采样零违例称为 certified nonlinear robustness。
7. 必须区分完整 Born tangent pair $(A_{\mathrm{Born}},B_{\mathrm{Born}})$
   与 mixed isolation control $(A_{\mathrm{Born}},B_{\mathrm{full}})$；后者不能
   称为完整 Born-SLAM Fisher 模型。

## 三、需要你解决的困难理论包

请先做 feasibility triage，将每题标为：可完整证明、可在附加假设下证明、可给反例、或当前开放。随后按“对论文主结论影响最大”的顺序完成尽可能多的严格结果。

### P1. 无限维 Hilbert 空间版本与稳定横截性

设地图空间 $\mathcal H_\chi$、实位姿空间 $\mathcal H_X$、白化实数据空间 $\mathcal Y$，且

\[
A:\mathcal H_\chi\to\mathcal Y,
\qquad B:\mathcal H_X\to\mathcal Y.
\]

逆散射中的 $A$ 往往紧且 $\operatorname{Ran}(A)$ 不闭，因此
$K_{\mathrm{IS}}^{-1/2}$ 可能无界。请：

1. 用 polar decomposition $A=U_A|A|$ 在
   $\overline{\operatorname{Ran}(A^*)}$ 上严格定义
   \[
   R_{\mathrm{geom}}=U_A^*P_{\operatorname{Ran}(B)^\perp}U_A.
   \]
2. 证明它与有限维广义 Rayleigh quotient 的关系，并精确说明
   $\operatorname{Ran}(B)$ 不闭时应使用闭包还是别的对象。
3. 研究
   \[
   \gamma=inf_{Au\neq0}
   \frac{\|P_{\overline{\operatorname{Ran}(B)}^\perp}Au\|}{\|Au\|}
   \]
   与 Friedrichs angle、闭和空间
   $\overline{\operatorname{Ran}(A)}+\overline{\operatorname{Ran}(B)}$、
   以及 form inequality
   $K_{\mathrm{SLAM}}\succeq\gamma^2K_{\mathrm{IS}}$ 的充要关系。
4. 区分“精确交集为零”和“稳定横截 $\gamma>0$”；构造一个交集为零但
   $\gamma=0$ 的紧算子反例，说明有限网格最小角为什么可能虚假稳定。
5. 给出可写入论文的定理、完整证明、反例和有限维离散逼近的条件（如 collective compactness、gap convergence 或 frame bounds）。

### P2. 连续 Helmholtz 的 $SE(2)/SE(3)$ gauge 定理与破缺

在均匀无限背景中，同时刚体变换场景 contrast 与整条收发轨迹应保持数据不变。请：

1. 在合适的 Sobolev/函数空间中定义群作用
   $g\cdot(\chi,X)$，严格证明前向映射 equivariance/invariance。
2. 对群轨道求导，推出每个 Lie algebra 生成元满足
   \[
   A\xi_\chi+B\xi_X=0,
   \]
   从而得到 gauge tangent 在 Schur-complement FIM 核中。
3. 说明径向对称场景下旋转生成元为何可能成为 pose-only stabilizer，而不是一个非零 map gauge direction。
4. 分别处理有限成像窗、已知非均匀背景、固定 anchor、已知接收器、地图支撑约束和离散像素基如何破坏或近似破坏 gauge。
5. 构造 quotient parameter manifold / horizontal space，并给出 quotient FIM 或投影后的信息算子；说明 motion-graph global gauge 与 map gauge 如何相互作用。

### P3. Born 远场 Fourier 几何、missing wedge 与显式角界

在 plane-wave incidence、far-field observation 和 Born 近似下，把观测写成
$\widehat\chi$ 在 Ewald/Fourier 采样集合上的取值乘相位因子。请：

1. 显式推导地图平移/旋转生成元与 Tx/Rx 位姿生成元在数据空间的关系。
2. 对直线、圆、部分圆弧和有限孔径，刻画
   $\operatorname{Ran}(A)\cap\operatorname{Ran}(B)$ 或其近似交叠。
3. 给出 missing-wedge/低频/有限带宽下的近混淆模式；若不可能得到统一闭式式子，请给可证明的上、下界。
4. 尝试从孔径、带宽、MIMO 数量、最小 standoff 和采样 frame bound 推出最小主夹角或稳定横截常数的下界。
5. 明确哪些结果只在有限维参数化（Fourier/spline 子空间）上成立，哪些能推广到函数空间。

### P4. 多频共享位姿补偿的横截性定理

多频堆叠中，同一个位姿扰动必须同时解释所有频率块：

\[
A_f u=B_f h\quad\text{对所有 }f.
\]

请：

1. 给出共享补偿核
   \[
   \{u:\exists h,\ A_fu=B_fh\ \forall f\}
   \]
   的精确维数/秩表达。
2. 在解析依赖于频率、有限维 map basis 和非退化阵列条件下，证明或反驳“除全局 gauge 外，几乎所有有限频率选择都会使交集横截”的 genericity 命题。
3. 给出频率重复、尺度复制、窄带极限和对称场景的反例。
4. 解释为什么绝对 $K_{\mathrm{eff}}$ 可单调而归一化 $\rho_i$ 非单调；寻找关于频率集合的可证明 surrogate bound，而不是声称“带宽越宽越好”。
5. 若能做到，给出可用于频率选择的 submodularity、weak submodularity 或 frame-potential 条件；若一般不成立，请提供反例。

### P5. 半正定位姿先验、运动图 gauge 与 weighted retention

令 $J_X\succeq0$，允许它是未锚定 odometry/pose-graph Laplacian。请：

1. 严格证明 generalized Schur complement 使用伪逆时的 range condition，并说明本问题中何时自动成立。
2. 在 $\ker J_X$ 与 $\operatorname{Ran}J_X$ 上分解
   $K_{\mathrm{eff}}$，给出自由 gauge、有限收缩和已知位姿极限。
3. 对
   \[
   R_X=Q_A^T\left[I-B(B^TB+J_X)^\dagger B^T\right]Q_A
   \]
   建立谱界、关于 $J_X$ 的 Loewner 单调性、严格单调条件和等号条件。
4. 寻找 prior-whitened canonical-correlation 表达；精确说明何时还能写成某种广义
   $\sin^2\theta$，何时不可以。
5. 处理 block-banded motion prior、跨时刻相关误差、map-pose 相关先验，并指出普通 Schur 公式需要怎样修改。

### P6. 固定秩流形、广义谱导数、重根与秩事件

请把以下离散观察提升为严格矩阵/算子命题：

1. 在固定秩 stratum 上证明 $B\mapsto BB^\dagger$ 的光滑性和投影导数公式，并给出
   $\|\dot P_B\|$ 关于谱隙/$\sigma_{\min}^+(B)$ 的界。
2. 推导 $K_{\mathrm{SLAM}}$、有限先验 $K_{\mathrm{eff}}$ 和广义保留谱的完整一阶导数，包括
   $\dot K_{\mathrm{IS}}$。
3. 对重根簇证明压缩矩阵
   $V_c^T(\dot K-\rho\dot K_{\mathrm{IS}})V_c$ 给出一阶分裂；进一步给 spectral projector 的 Kato/Davis--Kahan 型界。
4. 对 rank change 给出分段光滑、跳变或 Hölder 连续的精确分类。解释人为例子
   $B(t)=U\operatorname{diag}(s_1,\ldots,s_{r-1},|t|s_r)V^T$ 在 $t=0$ 为什么产生投影跳变。
5. 构造适合数值在线跟踪的 rank-event detector：相对阈值、hysteresis、谱簇合并规则和可证明的误报/漏报条件。

### P7. Full-wave resolvent 的显式扰动界与真正的鲁棒性证书

令 $M=I-D_\chi G_D$，并假设目标域与 Tx/Rx 保持最小距离、
$\sigma_{\min}(M)\ge m_0>0$、contrast 和频率有界。请：

1. 利用 resolvent identity 对
   $A=G_SM^{-1}D_E S_\chi$、完整 $B$、
   $DA,DB,D^2A,D^2B$ 给出显式可计算的范数界。
2. 对二维 Hankel Green 函数及其一、二阶空间导数给出在最小距离约束下的统一界；区分 self-cell、域内奇性和外部传感器项。
3. 由这些界推出固定秩且正谱隙区域内的
   $DK_{\mathrm{eff}}$、$D^2K_{\mathrm{eff}}$ 或 operator-Lipschitz 常数。
4. 明确常数如何随 $m_0^{-1}$、波数、contrast、standoff、域大小和
   $\sigma_{\min}^+(B)$ 爆炸；讨论近共振“信息变大”与“敏感度更快变大”的竞争。
5. 若无法在连续模型给出有限统一常数，请证明为什么，并把可证结果限制到离散网格/紧致参数球。不要用随机采样替代统一证明。

### P8. Robust trajectory 的一阶/二阶定理与约束锥

研究

\[
\max_{X\in\mathcal X}
\min_{\Delta X\in\mathcal U(X)}
\lambda_r\bigl(K_{\mathrm{eff}}(X+\Delta X)\bigr)
\]

或归一化保留率版本。请：

1. 对简单特征值和一般范数球证明
   \[
   \min_{\|\Delta X\|\le\varepsilon}
   \lambda_r(X+\Delta X)
   =\lambda_r(X)-\varepsilon\|\nabla_X\lambda_r\|_*+O(\varepsilon^2),
   \]
   并给出可计算的 remainder 常数。
2. 对 motion/collision/equality constraints 的 tangent cone 给出相应支持函数和最坏方向。
3. 对重根簇，将最坏一阶变化写成压缩矩阵线性组合的 min-eigenvalue 问题；说明何时是凸 SDP/非凸问题。
4. 在 rank event 附近给出安全 trust-region 或 nonsmooth formulation，而不是使用失效的光滑梯度。
5. 由于轨迹指标排序已被数值反例否定，请提出不依赖“圆优于直线”之类先验排序的多目标/Pareto 设计：至少同时报告绝对信息、归一化保留、standoff、能量和鲁棒性。

### P9. SOM/current 子空间与 map-tangent/pose-defect 子空间的严格桥接

已知离散全波 Jacobian可写成

\[
A=G_ST_\chi,\qquad T_\chi=M^{-1}D_E S_\chi.
\]

请：

1. 给出 $\operatorname{Ran}(A)\subseteq\operatorname{Ran}(G_S)$ 的精确条件，以及等号、严格包含和闭包版本。
2. 用 singular-value inequalities、Wedin/Davis--Kahan 或 graph-subspace 理论，界定
   $G_S$ dominant current modes 与 $A$ map-tangent modes 的 Grassmann 距离，常数应依赖
   $T_\chi$ 的条件数、核与谱隙。
3. 研究先做 SOM dominant/complementary split、再做 pose projection，与先做 pose elimination、再做 current reduction是否可交换。给出投影可交换的充要条件和不交换反例。
4. 对近似交集构造（例如连续投影、CS decomposition、randomized SVD）给出误差界，禁止用未经证明的“八象限”或 projector 乘积当成交集投影。
5. 最终提出一个数学上合法、CPU 可实现的 two-fold reduction 算法，并明确它近似的是哪一个算子和哪一种误差。

### P10. 从局部谱到实际估计误差：能证明到哪里

现有非线性最小二乘 toy 说明局部 Fisher 预测与有限噪声/非线性估计误差可能相差很大。并且固定真实位姿、只重复采样数据噪声时，带 pose penalty 的 frequentist covariance 应是 sandwich covariance，而不是 Bayesian marginal $K_{\mathrm{eff}}^{-1}$。请：

1. 严格区分三种实验/统计语义：固定 nuisance 的重复采样 covariance、从 prior 生成 nuisance 的 Bayesian marginal、以及 misspecified/penalized M-estimator sandwich covariance。
2. 推导各自的 map covariance，说明 $K_{\mathrm{eff}}^{-1}$ 在哪些生成模型下成立。
3. 给出局部渐近正规性、small-noise LAN 或 Gauss--Newton basin 条件，使 retention spectrum 可以预测实际误差膨胀。
4. 构造 cycle skipping、多峰或弱辨识反例，展示 Fisher 谱很好但全局估计失败的情形。
5. 提出论文中最强但不过界的桥接表述；若没有足够条件，明确把结果限定为局部 tangent information geometry。

## 四、文献与新颖性任务

请进行面向上述 P1--P10 的学术检索，但必须遵守：

1. 优先原始论文、正式期刊/会议和权威专著；给 DOI、arXiv 或稳定链接。
2. 对每条关键文献说明它精确覆盖了哪个对象：inverse scattering with unknown geometry、radar/RF SLAM FIM/PCRB、autofocus、blind calibration、bilinear inverse problem、principal angles、variable projection、operator perturbation、active sensing、SOM/current decomposition 等。
3. 用 claim ledger 输出“本工作命题—最接近已有结果—相同点—不同点—仍可主张的窄贡献—证据强度”。
4. 若找不到并不代表不存在；必须写成“本次检索未发现”，不能写“没有已有工作”。
5. 不要把一般 Schur complement、principal angles、Fisher information、projector derivative 或 SOM SVD 本身当作新颖性。

## 五、最终输出格式

请按以下顺序交付，使用中文解释、英文数学术语可保留：

1. **一页结论摘要**：哪些主命题成立、哪些需缩小、哪些被反例否定。
2. **假设总表**：每个定理所需空间、闭性、秩、谱隙、可逆性、边界、噪声和参数化假设。
3. **定理—证明—反例正文**：每个结果编号；证明不能跳过关键步骤。
4. **困难问题状态表**：P1--P10 分别标注 solved / conditional / counterexample / open。
5. **对当前论文的可直接替换文本**：给出定理、remark、limitations 和 appendix 段落，LaTeX 可直接使用。
6. **设计建议**：只给有数学依据的算法/实验修改，附复杂度和失败触发条件。
7. **文献 claim ledger**：逐项、可核验、有链接。
8. **尚未解决清单**：每项写出最小下一步、所需工具/引理和为什么当前不能完成。

最后请执行一次“反过度主张审计”：逐句检查是否把有限维当连续、把局部当全局、把内部控制当外部验证、把 sampled 当 certified、把 current-space 当 map-tangent-space、或把检索不足当作新颖性。发现任何一项就主动修正。
