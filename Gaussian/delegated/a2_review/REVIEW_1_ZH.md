# Gaussian A2 隔离 TAP 审查（Review 1）

日期：2026-09-11  
性质：独立内部审查，不构成期刊接收意见。  
裁决：**大修；当前不得声称参数空间二层 SOM、显式 current 方法或证书驱动方法相对完整 LM/GN 已有优势，也不得声称原创性已闭合。**

## 1. 审查范围与证据边界

本审查只读取任务指定的 A2 理论、三份检查代码、对应结果、先验技术审计和 `runs/a2/explicit_current/results.json`。`delegated/a2_som/EXPLICIT_CURRENT_SUMMARY.md` 不存在。未读取、未运行、也未推断仍在执行的 N64 imaging 或 measured campaigns。

`runs/a2/certificate/results.json` 自身明确标注为 “Independent reconstructed A2 discrete checks; NOT original Pro code rerun”。因此下述 36 项证书结果只能称为本地重构检查。Source Pro 附件仍未取回；现有归档是完整对话正文，不是原附件及其环境的复跑。

## 2. 实际得到支持的结果

1. 在固定有限维二次模型、无约束且 `H=A^T A+lambda I` 正定时，受限解 gap 恒等式成立。30 个随机检查的最大恒等式误差为 `4.88e-15`，坐标协同变换误差最大为 `9.82e-15`。
2. 梯度投影缺失率不是最优步证书。已给出的二维反例 `H=[[1,.9],[.9,1]]`、`g=(1,0)`、`V=span(e1)` 中，梯度缺失率为 0，而完整步为 `(5.263,-4.737)`、受限步为 `(1,0)`，模型 gap 为 `2.1316`。
3. 对二维标量、相同平移方形单元、张量积对称正权 GL 源积分及 `kh<=pi`，角向 form factor 非负，从而参考虚部矩阵是 Gram 半正定矩阵。N8/N12、`kh=0.2,1.0,2.9` 的 12 个浮点谱检查中，扣除自项 shift 后最小余量为 `-1.04e-14`；`kh=3.5` 正确返回“不适用”。这支持代数推导及当前样例的一致性，不是区间算术或连续 Maxwell 证书。
4. 36 个本地重构离散证书样例全部覆盖所测输出误差；嵌套 Galerkin 等价最大相对误差为 `3.48e-14`。这些结果不认证网格误差、模型失配、材料真值或反演结构。
5. 10 个现有 synthetic current 记录中，adaptive tangent 方法相对 ordinary LM 在材料误差和留出散射误差上各胜 7/10，但 wall time 和累计 forward+tangent RHS 各只胜 1/10。显式 current 方法相对 ordinary LM 在材料误差、留出散射误差和 wall time 上均为 0/10 胜，其最终 state relative 位于 `0.564--0.668`。这不支持当前实现优势。

## 3. TAP Reviewer A：电磁条件与证书闭合

### A-1（高严重度）证明对“结构化构造”成立，函数却对任意同字段对象发证

受影响：`Theory/a2/QUADRATURE_STABILITY_EXTENSION_ZH.md` 的实现落地声明；`code/a2/certificate.py::gamma_for_fft_kernel`；`code/a2/check_quadrature_certificate.py`。

角向 Gram 证明本身在所列窄条件下成立：GL 节点成对、权重为正，且 `kh<=pi` 时每个余弦因子非负。问题在于 `gamma_for_fft_kernel` 只读取 `cell_integrated/k/h/quadrature_order` 和一个对角 kernel 值，并不验证所有非对角项确由证明中的同一 GL 平移卷积生成。精确反例是：保留这些元数据和自项不变，在一对非对角位置加入 Hermitian 虚部扰动；该零对角扰动具有一正一负特征值，幅度足够大时会使下界失效。函数仍返回原来的正 `gamma`，但 `Im_H(D)-shift I` 已不再半正定。因而它不是针对“传入 kernel”的自包含证书，只是对受信构造器的条件性公式。

必须修订：将接口绑定到经验证的唯一核构造与版本/参数哈希，或在发证前验证生成结构；若要容纳任意对象，则必须计算有舍入余量的 Hermitian 下界。报告中把“certificate for this kernel”改成“certificate conditional on the trusted FFTGreen construction”。不得用 12 个实例测试替代结构验证。

### A-2（中严重度）电磁号约定、内积和离散未知量仍需成为正式假设

受影响：`Theory/a2/QUADRATURE_STABILITY_EXTENSION_ZH.md` 的命题及推论；`code/a2/certificate.py::certify`。

推论把 `Im(chi)>0` 当作耗散，并在 Euclidean 坐标中使用伴随和范数。这还依赖 Hankel-1 对应的时间号约定、未知量是单元常值/中心配点、所有单元质量相同，或已经严格 mass-whitened。若换时间约定、非均匀单元、接收端再积分或向量 Maxwell，符号和伴随均可能改变。当前文本虽排除了部分泛化，但没有把这些条件全部写进定理假设和运行时元数据。

必须修订：定理前显式声明时间因子、被动材料符号、源/测试离散、质量内积与自项定义；结果文件保存这些元数据。任何条件不匹配时返回“不适用”，不能沿用当前 Gamma。

## 4. TAP Reviewer B：数值公平性、代价与局部论证

### B-1（高严重度）局部 gap 恒等式正确，但不覆盖当前带信赖域/接受拒绝的算法

受影响：`Theory/a2/PARAMETER_SOM_AND_PORTS_ZH.md` 第 3 节；`code/a2/check_parameter_geometry.py`；`runs/a2/explicit_current/results.json` 中 adaptive rank 逻辑。

恒等式只比较同一固定线性化、同一阻尼、无边界/无信赖域时两个精确二次最小点。当前结果包含 `scaled_trust_radius`、多次 trial step、rejected step 和非线性重线性化；因此该恒等式既不给出实际接受步的 objective gap，也不给出全局轨迹排序。更严重的是，实现仍使用已被二维反例否定的 `adaptive_gradient_deficit` 触发增秩：缺失率可为 0 而最优步大量落在子空间外。

必须修订：把该量明确降级为未认证启发式；对实际信赖域/约束子问题记录 KKT 残差、predicted/actual reduction、边界是否激活，并以矩阵自由完整步或可靠误差估计作后验对照。不得把 30 个无约束随机恒等式检查解释为当前 adaptive rank 的正确性验证。

### B-2（致命实验门槛）现有比较不支持净优势，且停止准则与成本预算未匹配

受影响：`Theory/a2/PARAMETER_SOM_AND_PORTS_ZH.md` 第 6 节及任何算法优势表述；`runs/a2/explicit_current/results.json`。

10 个样例显示 adaptive 方法有时精度较好，但代价通常更高。moderate 五例的中位 wall time 为 LM `1.366 s`、adaptive `2.322 s`，累计 RHS 中位数为 `1008` 对 `1680`；strong 五例为 `1.462 s` 对 `3.287 s`，RHS 为 `840` 对 `1848`。两者的迭代次数、trial 次数和停止状态不一，不能把终点精度直接当作同预算优势。显式 current 方法的 state relative 始终大于 0.56，且材料/留出误差对 LM 无一胜出；当前实现尚不是竞争基线。

必须修订：提供同数据、同初值、同先验、同噪声实现下的 operator-equivalent budget 曲线（误差对 RHS、FFT、Jv/J*v、内存和 wall time），重复种子及不确定性，并加入完整矩阵自由 GN/LM、单层截断、等维随机补方向、GSVD/成熟 TSOM。建基、SVD、刷新、失败回退和有限差分 profile 成本必须全计。N64/measured 完成前保持“待验证”。

### B-3（中严重度）离散证书覆盖成立，但尚未证明在算法中有实用紧度

受影响：`runs/a2/certificate/results.json` 及证书作为更新/端口控制工具的表述。

36 项重构检查均覆盖，这是正面结果；但 `bound_rel` 范围为 `0.0129--11.16`，bound/corrected-error 比为 `3.65--19.44`，且最大 corrected relative error 为 `0.874`。覆盖本身不表明证书足够紧到能减少端口、停止 Krylov 或接受材料更新。

必须修订：预先定义证书实际决策（停止、增秩或拒绝）及容差，在相同任务中报告它相对未认证策略节省的工作量和错误接受/拒绝率；继续明确这 36 项是重构检查而非 Pro 附件复跑。

## 5. TAP Reviewer C：信息、原创性与实践含义

### C-1（中严重度）`Fisher=A^T A` 需要更强噪声假设

受影响：`Theory/a2/PARAMETER_SOM_AND_PORTS_ZH.md` 第 2 节。

“协方差固定、参数无关”不足以仅用协方差白化后断言 Fisher 信息等于 `A^T A`。例如同方差 Laplace 位置噪声的 score Fisher 与 covariance-whitened 二次量不同；多通道非同分布时差异甚至不是一个全局常数。

必须修订：假设 proper complex Gaussian/实化 Gaussian 加性噪声，或把 `W` 定义为 score-Fisher whitening 而非 covariance whitening。核心结论可以保留：确定性内部量 `B` 不是新观测，不能增加 observed-data Fisher 信息；但必须限定在局部一阶、给定似然下，不能外推为有限幅度不可辨识性。

### C-2（致命主张门槛）二层内部传播只能是算法选择，不能成为新增测量信息或分辨率证据

受影响：`Theory/A2.md` 主线；`Theory/a2/PARAMETER_SOM_AND_PORTS_ZH.md` 第 1--2 节；所有 SOM 优势表述。

若状态被候选参数精确消去，`B` 是模型导出的确定量。任取一阶零空间方向 `q` 使 `Aq=0` 而 `Bq` 很大，外部似然的一阶导数仍为零。故二层选向最多改变正则、预条件、近似误差或非线性路径；它本身不能证明观测分辨率增加。当前 10 例也没有隔离第二层贡献：最终 rank 常接近 12，总 tangent 成本较大，并缺少等维随机补方向和只增第一层的 matched ablation。

必须修订：删除“增加信息/分辨率”的因果措辞；把假设写成“在给定成本下改善非线性优化或状态近似”，并用 matched ablation 验证。若声称有限幅度可分辨，必须直接报告 `F(theta+tq)-F(theta)` 相对噪声的曲线，而不是引用 `Bq`。

### C-3（致命原创性门槛）最近邻闭合不足，当前材料主动否认核心部件的独创性

受影响：`delegated/a2_audit/PRIOR_ART_ZH.md`；A2 的贡献定位。

本地审计已指出：参数/状态联合自适应与认证降阶、伴随修正误差界、TSOM 内部传播选向、FFT-TSOM 边界改善均已有直接近邻；参数切空间拉回、n-port 代数和局部 gap 目前主要是已知线性代数/隐函数求导的应用。两次 ScholarQA 查询被 429 限流，指定 2023 TAP 论文仅取得题录线索，且 Source Pro 原附件未取回。因此既不能据检索失败声称首次，也不能把重构检查当作原始技术主张的复现。

必须修订：取回并存档原附件；完成逐项 claim chart（假设、算子、选向、证书、复杂度、实测协议）；把可能贡献收缩为可证伪的组合命题，例如“联合材料方向残差与状态证书在大参数重复反演中产生净成本和错误结构控制收益”。在这些证据完成前，原创性只能评为未证。

## 6. 评分与硬门槛

|维度|分数（1--5）|理由|
|---|---:|---|
|正确性|3|GL Gram、二次 gap、残差证书核心代数在窄假设下基本成立；但证书接口未闭合结构假设，Fisher 假设不足，局部恒等式不能认证当前信赖域算法。|
|原创性|2|现有材料把多数核心归入已知方法的具体化；检索和指定近邻对照尚未闭合。|
|实验|2|有小规模重构线性代数/谱检查和 10 个 synthetic current 记录；无已完成的大参数、重复统计或 measured 结论。|
|实践价值|2|证书可能有工具价值，但尚无决策节省；adaptive 成本大多更高，显式 current 当前明显不竞争。|
|表达|4|理论边界和“不增加信息”等限制写得较清楚；仍需统一“证书”与“条件性诊断”的术语。|

**不计算平均分。** 以下三项是独立致命门槛，任一未解除都不能作正面总裁决：

1. 没有同预算、同停止规则、包含成熟基线的净算法优势；
2. 原创性检索和 claim-by-claim 对照未闭合，Source Pro 原附件未取回；
3. 发证接口尚未验证传入离散核满足 Gram 证明的结构假设。

当前最稳妥的科学结论是：A2 已形成若干条件明确的有限维代数工具和反例，但尚未证明新的电磁可辨识性、算法优势或实践收益。
