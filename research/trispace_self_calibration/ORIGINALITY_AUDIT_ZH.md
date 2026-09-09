# TriSpace SOM 自校准：原创性与“车轱辘话”审计

日期：2026-09-04

> **终稿状态覆盖说明。** 本文件第 1--7 节保留了实验前的原创性假设；第 8 节是
> 三轮自动科研和独立复核后的最终判定。凡二者冲突，以第 8 节为准。尤其是：
> `T_U` 的判别机制与 reduced self-calibration 的算法收益均未成立。

结论先行：这个方向有一个窄而真实的论文候选核心，但绝不是“等效电流”、
“未知阵列自校准”或“联合重建”本身。真正可能形成科研价值的，是把阵列几何
扰动在 **指定的 SOM/TSOM 保留电流空间** 内作规范提升，同时保留无法提升的
数据残差，再用 full-wave 状态方程检查这个“等效电流”是否物理一致。当前判断
只对已检索文献集合成立，不是全球新颖性结论。

## 1. 明确属于前人方法或通用套话的部分

下列内容可以作为背景、工具或动机，但不能写成本文原创贡献：

1. **等效源、等效电流、contrast source 本身。** 用电流/等效源改写非线性
   逆散射是成熟思想。
2. **SOM 与 Twofold SOM 本身。** `G_S` 的奇异子空间、数据确定分量和模糊
   分量，以及进一步使用 domain/state 结构，均属于既有方法。
3. **self-calibration 这个题目本身。** 阵列处理、盲增益/相位标定、雷达
   autofocus、SAR 图像与相位误差联合估计均已有大量先例。
4. **联合估计图像/介质与发射机或接收机位置。** 逆散射里已有“重建 +
   transmitter localization”的直接先例；地震 FWI 也已有 source/receiver
   extension 和 receiver relocalization。
5. **把额外自由度加入优化器。** “把 pose 当变量，交替更新 map 和 pose”只是
   标准联合优化框架，不是贡献。
6. **Moore--Penrose 伪逆、range/null-space 分解、正交投影、主角、Wedin/
   Davis--Kahan 界、variable projection、Gauss--Newton、Schur complement。**
   它们都是标准数学工具。
7. **同时使用 data residual 和 state residual。** CSI/CC-CSI 家族已经这样做；
   新意不能落在“有两个残差”这句话上。
8. **“未知 `G` 会导致 SOM 子空间变化”。** 这是正确且必要的问题陈述，但单独
   看仍是直接推论，不足以构成论文贡献。
9. **“多频、多轨迹有助于标定”。** 若无明确的 shared-pose 条件、反例和量化
   结果，这只是常见直觉。
10. **把实测天线模型写回入射场和 data-equation Green operator。** Bellomo
    等（IEEE TAP 2014, DOI `10.1109/TAP.2014.2308534`）已经用实测 incident
    field 重建天线方向图，并在发射与接收两侧通过多极展开修正入射场和 Green
    operator；其文中还显式讨论 phase-center correction。因此“校准 `G` 后再做
    非线性逆散射”也不是新意。本文必须区别于这种预标定流程：目标是从散射数据
    中联合估计低维阵列 pose，并研究其在 SOM current coordinates 中的可辨识结构。

凡是只把这些内容重新排列、换名为 TriSpace，却没有新的可检验对象、定理或
算法收益，都属于用户所说的“根据前人研究生成出来的车轱辘话”。

## 2. 数学上成立但科研上空洞的核心陷阱

令名义 sensing operator 为 `S=G_S(X)`，阵列微扰产生

\[
v_h=(D_XG_S[h])j.
\]

如果离散问题中 `S` 是满行秩，而 current 像素数通常远大于接收通道数，则

\[
\operatorname{Range}(S)=\mathcal Y,
\]

任何 `v_h` 都能写成 `S delta j_eq`。因此：

> “阵列位置误差可以等效为某个电流变化”在常见单快照离散模型里往往自动
> 成立；它不提供 pose 可辨识性，也不构成新理论。

预实验已在现有 `N=8` Helmholtz harness 中看到该退化：四接收通道时
`rank(G_S)=4/4`，等效电流重建相对误差约 `9e-16`。这个结果应当作为论文的
反过度主张控制，而不是当作成功结果宣传。

## 3. 当前真正有价值的候选贡献

### C1. SOM/TSOM 受限的规范几何电流提升

给定一个明确声明的保留 current basis `U`，定义

\[
C_U=(G_SU)^\dagger H_S,\qquad
Q_U=UC_U,\qquad
R_U=(I-P_{\operatorname{Range}(G_SU)})H_S.
\]

这里 `Q_Uh` 是保留空间中的最小范数等效电流，`R_Uh` 是任何保留电流都无法
解释的数据部分。这一步把空洞的“存在某个电流”改成了一个依赖噪声度量、
截断规则和 SOM 坐标的规范、可计算、可反驳对象。

评价：**检索范围内未发现这个精确对象；是候选贡献，但单独仍主要是标准线性
代数的物理化。** 它必须与 C2--C4 合起来才足够有论文价值。

### C2. receiver pseudo-current 与 transmitter physical-current 的语义分离

在 independent-current 形式中，receiver motion 给出采样项
`H_Sh=(D_XG_S[h])j`；transmitter/illumination motion 给出真实 current 变化

\[
\delta j_{\rm phys}=M^{-1}H_Dh.
\]

只提升前者，后者保留在 state equation，可避免把同一个 pose effect 计算两次。
在 world-fixed homogeneous grid 下 `D_XG_D=0`，不能泛称所有 pose 都使
`G_D` 变化。

评价：**这是重要的语义和建模贡献候选。** 它未必单独新颖，但能把很多表面上
相似、实际上不同的“等效”说法严格拆开。

### C3. 状态一致性见证与精确 hiding condition

重参数化后得到

\[
\delta y=G_SU z+R_Uh,
\]

\[
MUz-D_{E_{\rm tot}}\delta\chi-D_Uh=0,
\qquad D_U=H_D+MUC_U.
\]

再令

\[
\mathcal N_D=MU\ker(G_SU)+\operatorname{Range}(D_{E_{\rm tot}}),
\quad
T_U=(I-P_{\mathcal N_D})D_U.
\]

有限维下，一个 pose 方向可被保留 current 与实 contrast 变化完全隐藏，当且仅当

\[
R_Uh=0,\qquad T_Uh=0.
\]

评价：**这是目前最强的窄理论候选。** 其证明是短的，但它把 SOM 受限提升、
full-wave state consistency 与 pose identifiability 连成一个类型正确的判据。
是否已有高度相似定理仍需对 source/receiver-extension FWI 全文做针对性核查。

### C4. graph-based TriSpace，而不是三个投影硬凑八格

建议把 pose 的三联对象定义为

\[
\Gamma_P^{(U)}=\{(Q_Uh,R_Uh,D_Uh):h\in\mathbb R^p\}.
\]

它同时记录 current-space lift、data-space irreducible residual 与 state-space
defect。TriSpace 因而是一个耦合图结构，而不是假设 `P_S,P_D,P_P` 交换后得到的
`2^3` 个交集。

评价：**这是概念组织与算法设计的候选贡献。** 若 E4/E5 不能证明它产生额外的
可辨识性或算法收益，它会退化成漂亮但无用的新记号。

### C5. 未知 `G` 下可执行的 pose-dependent SOM 坐标更新

每次迭代重建 `G_S(X_k)`、SOM basis、受限提升和状态见证；用 projector/
Procrustes 对齐固定秩子空间，监控 singular gap 与 rank event，并以 soft filter
或 trust region 防止伪逆爆炸。

评价：**是否构成算法贡献完全取决于 E5。** 必须对比 wrong-pose SOM 与直接
joint full-wave inversion，并报告 basin、失败样本、变量数和运行量；只展示一张
成功图不够。

## 4. 当前原创性判断

| 内容 | 当前判断 | 能否作为主贡献 |
|---|---|---|
| 等效电流/contrast source | 明确前人工作 | 否 |
| SOM/TSOM | 明确前人工作 | 否 |
| 自校准、joint map-pose | 明确前人工作 | 否 |
| source/receiver extension | 紧邻前人工作 | 否，必须正面引用 |
| full-row-rank 时可等效 | 标准线性代数；重要负结果 | 只能作边界/反例 |
| SOM-restricted canonical lift `Q_U` + `R_U` | 检索范围内未见精确对应 | 可作候选贡献之一 |
| receiver pseudo / transmitter physical 分离 | 检索范围内未见同样组织 | 可作建模贡献候选 |
| state witness `T_U` 与 hiding iff 判据 | 检索范围内未见精确对应 | 当前最强理论候选 |
| graph TriSpace | 新组合/组织尚待证明有用 | 有 E4/E5 支撑才可 |
| 非线性自校准算法 | 尚未验证 | 由 E5 决定 |

因此目前最诚实的论文定位是：

> 我们不是发明等效源或自校准，而是研究在未知阵列几何下，如何把 Green-operator
> 误差规范地提升到 SOM 保留 current space，并用不可吸收数据残差与 full-wave
> 状态一致性共同判断这种等效何时可用、何时造成不可辨识，以及它是否能形成
> 一个有效的 reduced self-calibration 方法。

## 5. 会使论文主张失败的结果

出现下列任一情况都必须主动收缩或放弃算法主张：

1. 在有意义的 retained rank 下，`R_U` 与 `T_U` 对所有 pose 方向仍同时接近零；
2. `T_U` 几乎总能被实 contrast tangent 吸收，无法比 data-only 提供额外区分；
3. `V_P=Range(Q_U)` 只是在 `V_S^+` 内换了一个名称，没有稳定谱量或预测能力；
4. rank/cutoff 的微小变化使 `Q_U` 完全不稳定，soft regularization 也无法恢复；
5. E5 中 geometry-lifted 方法不优于直接 joint inversion，或只有单一精心初始化
   才成功；
6. source/receiver-extension FWI 全文已经给出等价的受限提升加状态残差判据；
7. 实际 TSOM domain fold 与当前假定不兼容，使 TriSpace 的第二折无法落地。

即使算法主张失败，C1--C3 仍可能形成一篇“何时等效电流表述是空洞、何时由
状态方程恢复可辨识性”的理论/负结果论文；但不能继续声称完整 self-calibration
已经解决。

## 6. 投稿前最低门槛

1. 五组实验全部有可执行记录，尤其是 E5 的多初始化和失败区域；
2. 独立复核 C1 与 C3 的证明、复数 realification、符号和 gauge quotient；
3. 核对 Chen SOM 与 Zhong--Chen Twofold SOM 的原始公式，不能用自造的
   `G_D` split 冒充 TSOM；
4. 对最邻近的三条 source/receiver-extension FWI 与 transmitter-localization
   工作，以及 Bellomo 等的天线方向图/phase-center 标定工作做全文级比对；
5. 明确 nominal basis、pose-dependent basis 和 soft-filter basis 的差别；
6. 在摘要和贡献表中禁用“首次”“此前无人”“解决了未知阵列”等绝对措辞；
7. 若没有真实数据，明确标注为 finite-dimensional synthetic validation，而不是
   practical array calibration validation。

## 7. 证据边界

本审计使用 16 个已核对身份的论文锚点及项目已有文献材料；大多数证据为摘要/
出版页级，Semantic Scholar 匿名检索还受到 429 限流。故“未发现”只表示在当前
检索集合中未发现。最终新颖性判断必须在定理、实验和邻近论文全文核对完成后再
定稿。

## 8. 三轮自动科研后的最终原创性判定

### 8.1 可以保留为论文核心的内容

1. **有条件的物理语义分裂。** 在 world-fixed、均匀固定背景的
   independent-current 模型中，把 receiver re-sampling 的 `H_S h` 与
   transmitter/illumination 引起的真实电流响应 `M^{-1}H_Dh` 分开；只提升前者，
   避免 double counting。该内容有解析导数和有限差分验证，但不能外推到所有
   moving-grid 或 hardware model。
2. **受限而规范的等效电流坐标。** `Q_U=U(G_SU)^\dagger H_S` 与
   `R_U=(I-P_{Range(G_SU)})H_S` 给出声明度量和保留基下的最小范数提升与
   retained-basis leakage。线性代数本身是教材内容，科研价值在于把它放到未知
   `G` 的 SOM 语义中，并同时给出 full-row-rank 空洞性反例。
3. **数据侧 pose-hiding 维数障碍。** 在 `m` 维白化实化数据空间中，若 retained
   current 与 map nuisance 的实 range 为 `N_U`，则
   `hidden_dim >= max(0,p+dim(N_U)-m)`。对本实验的 `m=2M,p=3,K=3`，generic
   情形在 `r=M-2` 至少隐藏两维，在 `r>=M-1` 隐藏三维。这个结论比原先的
   “整圆 Fourier 跃迁”更一般，也纠正了旧的相对阈值 bug。
4. **负结果本身。** `T_U` 在测试区间对好坏位姿没有区分力；hard/soft state
   variants 也没有恢复被丢弃的 pose directions。对已知 `G_S` 的 SOM 文献而言，
   这说明“再加一个 state residual”不会自动得到 self calibration。
5. **设计原则。** 多发射源或有方向性的源可以消除单个各向同性点源的 orientation
   rank deficiency；这不是消除全局 SE(2) gauge。真正目标是让多帧/多频/多源堆叠
   后的 projected pose Jacobian 满秩且具有足够的最小奇异值。

### 8.2 已被否定、必须删除或降级的主张

| 原主张 | 最终状态 | 论文中的处理 |
|---|---|---|
| `T_U` 区分 data-equivalent 与 physical current | **测试中被否定** | 作为负结果和开放问题 |
| graph TriSpace 自身提供额外可辨识性 | **未证明** | 只作 typed bookkeeping，不作贡献定理 |
| `V_P` 是第三个独立空间 | **概念错误** | 改称 retained current 内的 pose-equivalent tangent |
| reduced 与 direct 相等说明算法更好 | **错误解释** | 改称满秩时的重参数化一致性检查 |
| `r=M-2` 是特殊 Fourier/Hankel 结构 | **无需该假设** | 改为一般维数障碍；对称模态解释仍开放 |
| `r>=M-1` 时 pose 又全部可见 | **数值阈值伪影** | 正确为 nuisance 饱和、三维全部隐藏 |
| 多源“打破 pose gauge” | **过度表述** | 仅称消除 source-symmetry rank deficiency |
| Phase-I 实验已经实现 TSOM | **不成立** | 统一称 SOM-truncated；TSOM domain fold 待原文核对 |
| 已解决未知阵列 self calibration | **不成立** | 定位为局部诊断、边界理论与 acquisition guidance |

### 8.3 最窄、最诚实的原创性句子

> 在当前有界检索范围内，我们未发现已有工作把未知阵列几何造成的 receiver-side
> Green-operator tangent，规范地提升到一个明确声明的 SOM 保留电流空间，同时
> 保留其 retained-basis leakage、区分 transmitter-induced physical-current
> response，并给出 retained-current/map nuisance 对 pose 可见维数的显式障碍。
> 本文还给出负证据：所测试的 full-wave state consistency 并不能自动恢复这些
> 被隐藏的位姿方向。

这句话仍是 **retrieval-bounded**，不能改写为“首次”“此前无人”或“证明了全局
新颖性”。Bellomo 等的天线/相位中心预标定、source/receiver-extension FWI、
blind calibration、joint map-pose recovery、SOM/TSOM 与 CSI data/state residual
均属于必须明确承认的前人工作。

### 8.4 论文形态判定

当前材料足以形成一篇 **有限维理论 + 可复现实验 + 负结果的初稿**，但还不够
称为完整 self-calibration 方法论文。若后续无法证明更强的 stacked acquisition
条件或在真实/高保真 Maxwell 数据上获得稳定收益，最合适的投稿形态是
identifiability/cautionary paper，而不是宣称新求解器优于 direct joint inversion。
