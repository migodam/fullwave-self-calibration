# A2 extensions 独立科学与代码评审

## 评审范围与结论

本评审以只读方式核查了：

- `Theory/a2/extensions/CORE_AUDIT_ZH.md`
- `code/a2/extensions/manifold_runner.py`
- `code/a2/extensions/manifold.py`
- `code/a2/extensions/theory_checks.py`
- `code/a2/physics.py` 中被调用的 VIE 正演与材料切向实现
- `code/a2/extensions/probability.py`、`probability_run.py`、`probability_v2.py` 及 v2 后处理/诊断文件
- 已存在的 frozen config、理论检查和运行中结果，仅用于核对实现与证据边界

总评：`CORE_AUDIT_ZH.md` 的主要数学纠偏是可靠而克制的，尤其是静态对象联合似然、峰值 Gaussian 的 SPD 导数、接收散射场二阶矩修正、以及“有限 secant 不能自动升级为证书”的边界。当前不能接受的是由 `manifold_runner.py` 导出的“弱方向有限探索在共同预算下优于基线”类结论：弱方向在满秩时选错了语义对象，且四法实际 RHS 不相等。这两点是该实验主张的致命阻断项，而不是对审计文档全部理论的否定。

本报告截取运行中的中间状态进行代码核查；流形结果文件仍在增长，因此不把其中任何汇总数字当作最终统计结论。概率 v2 主运行已有 30 条记录和匹配源码哈希，但后续初始化、容量和弱数据诊断明确标为 post-hoc，不能并入预冻结的确认性证据。

严重度定义：**致命**表示会使当前核心结论无效；**高**表示必须在对外结论前修正；**中**表示不会推翻局部公式，但会限制解释或复现；**低**表示文档/防御性改进；**通过**表示本次检查未发现实质错误。

## 严重度标注发现

### [致命] M1：满秩时仍把已纳入拟合的方向当作“弱方向探索”

证据：`manifold_runner.py:67` 令 `rank=len(s)` 时可见子空间包含全部 12 个右奇异向量；但 `manifold_runner.py:71` 随后使用

```python
direction = Vh[min(rank, len(s)-1)]
```

当 `rank == len(s)` 时，这等于 `Vh[-1]`，它已经包含在 `V=Vh[:rank].T` 中，并不是 `Range(V)` 的补方向。运行中快照也不是边缘情况：截至本评审快照，`weak_secant` 历史中有 201/391 次记录 `rank=12`；最终数字仍可能变化，但代码错误与比例无关。

后果：这些迭代不能作为“一阶不可见/被截断方向经有限位移显现”的证据。它们只是沿当前完整切空间最弱奇异向量做额外有限步试探。在 `rank=12` 时，`spectral` 也与 `full_lm` 使用同一个全空间，方法标签不再代表四个不同机制。

具体修正：

1. 只有 `rank < p` 时才允许 `direction=Vh[rank]`；`rank == p` 时记录 `no_weak_complement` 并跳过弱方向探索。
2. 若研究对象是“最弱但仍可见”的方向，应明确改名，并从可见拟合空间移除该方向后再构造 profile；不能同时把它放在 `V` 内又称作可见空间外探索。
3. 每次记录 `direction_in_visible_norm=||V^T h||` 与 `direction_complement_norm`，使语义可审计。

### [致命] M2：当前不是相等 RHS 预算比较，不能据此声称同预算加速/优势

RHS 账本的单项计算是正确的：一次双频正演为 `2×6=12` RHS；完整 12 参数双频切向为 `2×12×6=144` RHS。初始化为 12 RHS；普通方法每轮再花 156 RHS；两种 secant 方法每轮最多花 `144+2×12+3×12=204` RHS。代码中的累计计数与这些数一致。

但 `max_outer=20` 与统一的保守 guard 共同造成实际花费不相等：当前每个已完成 case 中 `full_lm`、`spectral` 都是 3132 RHS，而 `weak_secant`、`random_secant` 都是 3480 RHS，后两者多约 11.1%。因此这里最多是“同一上限下不同实际支出”，不是同预算比较。若按迭代数比较，差异更大，因为每轮成本不同。

具体修正：

1. 主比较使用完全相同的实际 RHS 截止点，或报告随累计 RHS 的质量曲线并在共同可达预算处比较。
2. 同时报告墙钟、GMRES 总迭代数和失败数；“一个 RHS”并不保证不同材料点的线性求解工作相同。
3. 不得把提前达到 `max_outer`、拒绝步或未花完预算解释成同精度加速。当前结果可保留为开发轨迹，但不能用于公平性能结论。

### [高] M3：随机对照与弱方向对照的方向分布不匹配

`random_secant` 在全部 12 维坐标中均匀抽取随机方向（`manifold_runner.py:72`），而理论上要比较的是截断可见子空间的补空间方向。随机方向一般含大量可见分量；弱方向在满秩 bug 下甚至完全可见。因此“两者候选数相同”不足以构成机制公平的随机对照。

具体修正：当 `rank<p` 时令随机基线使用 `(I-VV^T)u/||(I-VV^T)u||`，并与弱奇异方向保持相同步长、正负号、retraction、visible refit、裁剪和正演验收；另设全空间随机搜索时应作为第三种不同基线命名。

### [高] M4：有全波单调筛选，但没有可称为信赖域/LM 的充分验收规则

正面结论：每个最终 proposal 都在 `manifold_runner.py:84-85` 重新解完整 VIE，并仅在训练残差平方严格下降时接受；探索 probe 与未选 proposal 的 RHS 也被计入。因此不存在“只用线性模型接受候选”或“只计成功候选成本”的问题。

限制：验收只使用 `val < best[0]`，没有相对下降容差、Armijo 条件或 predicted/actual reduction ratio。运行快照中多个方法的最终目标在约 15 位有效数字内相同，接受次数却不同，说明浮点级变化可以改变 `accepted`。damping 也仅按布尔接受乘 0.7/4，而不是真正的 LM ratio 更新。故应称“damped GN proposal + exact forward monotone selection”，不能据此称标准 trust-region LM。

另一个一致性问题是：secant 的 `delta` 在原始 `h` 上计算，随后 `h+visible(...)` 可能被整体裁剪至半径 0.7；裁剪后的探索分量已不是计算 `delta` 的位置。最终全波验收保证单调性，但局部 profile 的解释已失效。

具体修正：加入尺度相关的最小相对下降与实际/预测下降比；裁剪后重新构造或重新评估 profile；或者保留现实现但收缩名称与理论解释。

### [中] M5：导数实现彼此一致，但扩展自身缺少“材料 chart → VIE 场”端到端有限差分证书

逐式检查结果：

- `a'=a exp(d0)` 对 `d0` 的导数是峰值场 `g`，与 `material():24` 一致。
- `mu'=mu+0.03 dmu` 对中心的导数是 `0.03 g Sigma^{-1}(x-mu)`，符号和尺度与 `material():24` 一致。
- `Sigma'=Sigma^(1/2) exp(H) Sigma^(1/2)` 在零点的物理切向是 `Sigma^(1/2) H Sigma^(1/2)`；`material():25-26` 使用了这一方向，且峰值 Gaussian 导数正确地不含 trace 项。
- `physics.py:195-200` 的材料切向 `dj=M^{-1} diag(dchi)(E+D j)` 与状态方程 `(I-diag(chi)D)j=diag(chi)E` 的微分一致。

已有 `check_manifold_native.py` 检查了 12 个材料列的 chart 有限差分和机械性的 RHS/单调性，但没有把 `material_tangent(...)[scattered]` 与 retraction 后完整 `forward` 的中心差分逐列比较。因此本评审确认“公式/代码一致”，不等于已有独立的端到端场导数证书。

具体修正：在小网格上对 12 列进行复散射场中心差分，分别报告每频率、每参数块误差；测试源码哈希应包含 `manifold_runner.py`、`manifold.py` 和 `physics.py`，而不只哈希检查脚本。

### [高] M6：所谓物理度量实际上是一个未声明权重的乘积坐标度量，弱谱结论依赖该选择

代码用欧氏 `||step||_2` 做裁剪，并对该坐标下的 Jacobian 做 SVD。其隐含度量可写为每个 Gaussian 的

\[
(d\log a)^2+\|d\mu\|^2/L_0^2+
\|\Sigma^{-1/2}(d\Sigma)\Sigma^{-1/2}\|_F^2,
\qquad L_0=0.03\ \mathrm m,
\]

其中 `SYM_BASIS` 的非对角项按 `1/sqrt(2)` 归一化，所以三维系数范数确实对应 Frobenius 范数。该选择对四种方法相同，也把幅度、位置和 SPD 切向无量纲化；若长度单位变换时同时变换 `L0`，位置/SPD 部分可保持单位不变。

问题是协议和 frozen config 没有声明这个度量，`0.03` 只是源码常数；它也不是自动得到的 Gaussian Fisher 度量。奇异值排序、`0.08*s[0]` 截断和“弱方向”都会随 `L0` 及三个参数块的相对权重变化。因此不能笼统写成“物理度量下的弱方向”或“坐标不变”。

具体修正：把 `translation_scale_m` 和各参数块权重冻结到 config，写出上述乘积度量及单位变换规则；对外只称该声明度量下的结果。若要声称度量稳健性，需要预先规定灵敏度分析，而不能事后挑尺度。

### [中] M7：“谱 rank continuation”没有 continuation 机制

`rank=sum(s>0.08*s[0])` 每轮独立硬阈值，没有跨轮秩调度、滞回或误差准则。当前实现是 truncated-SVD damped GN；秩可随谱上下跳动。建议按实际算法改名，或显式冻结 continuation 规则。该问题不影响一次 SVD 的数学正确性，但影响方法归类和新颖性叙述。

### [通过/边界明确] T1：接收散射场二阶修正与余项界在声明的有限维条件下正确

令 `Delta X=diag(delta chi)`，并明确

\[
C_p=\mathbb E[\delta\chi\,\delta\chi^T]
\]

是无共轭的 pseudo-covariance。则

\[
\mathbb E[\Delta X B\Delta X]=C_p\odot B,
\]

且恒等式

\[
j-j_0=U\Delta X(I-B\Delta X)^{-1}E_0
\]

给出审计文档中的

\[
\mathbb E[F]=SME_0+SU(C_p\odot B)E_0+R_F.
\]

若逐 realization 有 `||B Delta X||<=q<1` 和 `||Delta X||<=delta`，则

\[
\|R_F\|\le \|SU\|\,\delta\frac{q^2}{1-q}\|E_0\|
\]

也正确。审计文档准确地区分了这一充分条件与物理可逆性、Gaussian likelihood 或后验误差证书。

### [高] T2：`theory_checks.py` 只验证了一个 2×2 代数例，尚未把矩界变成实际 VIE 的 full-wave 证书

`theory_checks.py:19-27` 的矩阵是手工给定的 2×2 复矩阵；输出也诚实声明“no Maxwell validation claim”。它验证了该例中误差小于保守界，但没有对实际 VIE 状态族验证统一 `q<1`、`delta`、`||SU||`，因此不能从该文件推出所用散射实验满足 Neumann 条件。

此外，四个等权对称状态使三阶中心矩消失；输出误差随 `delta` 呈约四阶缩放。这并不检验一般分布下界的三阶尺度，只是仍被三阶上界覆盖。

具体修正：现阶段只保留“条件有限维定理 + 小矩阵单元测试”的措辞。以后若需要物理应用结论，应对声明的实际有限状态族给出统一常数，或明确报告条件失败；可再加入零均值但三阶矩非零的非对称分布测试。无需为本次正在运行的实验临时补跑。

### [高] T3：静态联合似然的文档语义正确，但所谓数值检查目前是硬编码常数

`CORE_AUDIT_ZH.md` 的

\[
p(y\mid\vartheta)=\sum_zP_\vartheta(z)\prod_\ell p(y_\ell\mid F_\ell(z))
\]

在“同一静态对象 z 被所有频率/照明共同观测，且给定 z 后各噪声块条件独立”时正确；若噪声跨通道相关，乘积应替换为相应联合噪声密度。`probability_run.py` 和 `probability_v2.py` 都把同一个 `z` 用于全部频率、Tx、Rx，再对通道残差求和，符合静态语义，并没有在照明间重抽样物体。

但是 `theory_checks.py:28` 的 `static_likelihood=.09` 与 `resampled_likelihood=.25` 只是直接写入字典，没有从任何状态、先验和条件似然计算，也没有 assertion。它不能充当静态/重抽样语义的可执行检查。

具体修正：构造至少两个静态状态和两个照明，显式计算“先对共同 z 求和”与“每照明分别求和后相乘”，验证两者不同，并检查实现选择前者。

### [中] T4：复随机材料的二阶记号应在公式处消除歧义

审计前文已正确区分 Hermitian covariance `K=E[delta chi delta chi*]` 与 pseudo-covariance `C=E[delta chi delta chi^T]`。但接收场修正式再次写 `C` 时没有就地重申它是后者。对复材料，`E[Delta X B Delta X]` 必须使用无共轭的 `C_p`；若读者误用 `K`，公式即错误。建议把该节全部改记为 `C_p`，并单独说明输出协方差计算仍需要 `K` 与 `C_p` 两者。

### [通过] T5：峰值 Gaussian 导数和保质量分裂修正正确

审计指出峰值参数化不应出现 `-tr(Sigma^{-1}H)/2`，这是正确的；该 trace 项只来自归一化因子。`theory_checks.py` 的有限差分也分别验证了峰值与归一化两种公式。保矩分裂的

\[
a_c=\frac a2\sqrt{\det\Sigma/\det(\Sigma-dd^T)}
\]

正确保持二维连续 Gaussian 的总质量，且必须有 `d^T Sigma^{-1}d<1`。Fourier 检查的四阶缩放对无限域连续 Gaussian 正确；它不是有限方形网格、截断尾部或全波场的验证，现有审计没有越过这一边界。

### [中] P1：概率 v2 的静态 posterior 实现正确，但当前主任务几乎完全后验塌缩，不能验证不确定性质量

主 v2 对 256 个有限状态枚举 N32 模型 likelihood，并用 N64 cell-integrated 数据生成同一静态状态的观测；exact、product VI 和 RBF-logit VI 都用 `E_q[F(z)]` 作预测，没有使用 `F(E_q[chi])`。ELBO/KL 恒等式的最大误差约 `1.82e-12`，实现层面通过。

然而当前 30 个主 v2 case 的 exact posterior patch Brier 在汇总中数值上均为零，entropy/均值也表现为几乎点质量；product VI 因能表示任意点质量而同样几乎精确。这一任务主要检验状态识别和枚举机械性，不能检验相关 posterior、概率校准或“联合分布比边际更重要”的困难情形。图中的 prior-predictive reliability 也只有 30 个对象、240 个相关 patch 观测，不能给固定对象重复噪声覆盖结论。

具体修正：保持 v2 为开发结果并收缩结论。后续确认性协议需要预先定义能产生非平凡多峰/相关 posterior 的信息受限情形，并把对象作为聚类单位报告区间；不要把当前 post-hoc 弱通道诊断回填成原冻结主实验。

### [高] P2：Gaussian-RBF VI 的巨大 KL 同时混合了族容量限制与非凸优化，不能解释为 Gaussian 概率表示的普遍失败

主 v2 的 Gaussian-logit 只有 4 个参数，而 product VI 有 8 个。post-hoc LP 容量诊断显示固定 RBF 设计只能让 256 个二进制状态中的 118 个成为可任意集中的符号模式；因此许多接近点质量的真 posterior 在结构上就不可由该族表示。单次零初始化 L-BFGS 又把族误差与局部优化误差混在一起；`optimizer_success=True` 只表示局部求解器停止成功，不表示找到了全局 ELBO 最大值。post-hoc 多起点确有改善，但文件自己已声明不替代 frozen v2。

具体修正：把结论严格限定为“该固定四维 RBF-logit 族及冻结优化流程在此有限字典上的结果”。若比较表示能力，应给出族内最优/多起点上界和容量覆盖；若比较算法成本，计入全部初始化、likelihood、256-state cache 建造与后处理成本。

### [中] P3：v1、v2 和 post-hoc 材料必须分层引用

- `probability_run.py`/v1 使用 N32 缓存既生成 truth 又做 posterior，是开发性的 inverse-crime 基准；其三个 N64 state 只给网格敏感性读数，并未改变主 likelihood。
- v2 使用独立 N64 cell-integrated 数据和 N32 字典，较好地避免上述问题，但“exact posterior”只对这个可能失配的 N32 有限字典 likelihood 精确，不是连续物理真 posterior。
- `probability_v2_optimization.py`、`check_rbf_capacity.py`、`probability_weak.py` 的标题/配置都明确是 post-hoc 诊断。它们可以解释失败机制，不能当作原冻结假设的独立确认。

建议最终证据表逐项标注 `pilot / frozen development / post-hoc diagnostic / confirmatory`，并固定引用对应源码与结果哈希。

## 三个明确 reviewer 维度

### 1. EM theory（电磁/全波理论）

评价：**条件通过**。

VIE 材料切向、峰值 Gaussian-SPD 复合导数、静态对象全波枚举，以及接收散射场的二阶 resolvent 修正均数学一致。审计也正确拒绝把 `q>=1` 解读为物理不可解，并拒绝从小均值误差跳到 posterior 可靠性。

主要缺口是证据层级：矩界只在任意 2×2 复矩阵上做了单元测试，没有证明当前 VIE/material family 满足统一 Neumann 条件；流形扩展也没有自己的端到端场导数有限差分。故可表述为“离散条件定理与实现公式成立”，不能表述为“实际 Maxwell/VIE 概率矩已获认证”。模型本身仍是二维标量 Helmholtz VIE，不应外推到三维向量 Maxwell。

### 2. Numerical fairness（数值公平性）

评价：**不通过，当前存在致命阻断**。

正向候选确实全部做全波验收、额外探针也计费，这是优点。但四法实际 RHS 不相等；随机方向的采样空间不匹配；满秩时“弱方向”已经在可见空间内；验收对浮点微降过于敏感；“spectral continuation”和“LM”名称均比实际算法更强。运行中的开发图和部分 case 不能修复这些设计问题。当前不得得出同预算速度、有限探索优越性或 chart 失灵结论。

### 3. Novelty（新颖性）

评价：**未建立，不应宣称理论原创**。

审计文档本身在这一点上是诚实的：Neumann/Taylor resolvent 恒等式、Gaussian SPD 流形、变分 Bayes、有限字典枚举、局部投影重拟合都不能仅因组合到全波散射中就称为原创理论。仓库的先验艺术笔记也已记录 Bayesian microwave inversion、Gaussian mixture manifold optimization 与 VarPro 等近邻，但当前检索不是穷尽性文献审查。

现阶段最多可提出待证伪的机制假设：在一个预先声明的物理度量和严格匹配预算下，有限 secant 经旧参数重拟合后是否仍提供超出成熟截断/随机基线的信息；以及联合概率近似误差是否会改变结构决策。当前实现尚未公平验证前者，概率 v2 又几乎是点质量 posterior，尚未验证后者。没有来源与确认性实验前，不得使用“首个”“原创定理”“普适提升”等表述。

## 建议的结论边界

当前可以保留的结论：

1. 峰值 Gaussian 的协方差导数不含归一化 trace 项；runner 中 12 个材料导数与所用 retraction 一致。
2. 在明确有限维、零均值与统一 `q<1` 条件下，接收散射场二阶修正及余项界成立。
3. 静态对象必须在所有照明共享同一状态后再边缘化；概率实现遵守这一语义。
4. 有限候选经过完整正演作了单调筛选，并计入了 probe/candidate RHS。

当前必须否决或暂停的结论：

1. “弱/不可见方向的有限探索已被本 runner 验证”——满秩分支直接违背该语义。
2. “四法在相同 RHS 预算下比较”或任何由此导出的加速结论——实际 RHS 为 3132 对 3480。
3. “当前检查给出实际 VIE 的 full-wave moment certificate”——只有条件定理和 2×2 单元测试。
4. “静态联合似然已由 `theory_checks.py` 数值验证”——对应两个数是硬编码。
5. “Gaussian 概率场一般劣于 product VI”——当前只是四维固定 RBF 族、几乎点质量 posterior 与特定优化流程的结果。
6. 任何理论原创性主张——现有材料只支持有限范围的邻近工作核查和待验证机制。

## 最小修正优先级

在不干扰当前运行的前提下，后续整理应按以下顺序进行：

1. 先修 `rank==p` 的弱方向定义，并匹配随机补空间基线。
2. 以相同实际 RHS/墙钟重新定义汇总方式；现有运行只作为 development trace，不要求现在重启。
3. 显式冻结乘积度量、长度尺度与方法准确名称。
4. 将静态似然硬编码改为真实计算，并补小网格端到端场导数检查。
5. 把 full-wave moment 结果限定为条件定理，除非以后能给实际 VIE 的统一常数。
6. 将概率 v2 主结果与所有 post-hoc 容量、初始化、弱数据诊断严格分栏。

最终科学裁决应由根任务负责人结合运行完成后的冻结哈希和预先声明的统计协议作出；本评审不把运行中的中间文件视为最终验收。
