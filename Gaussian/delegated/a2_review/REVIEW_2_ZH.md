# Gaussian A2 隔离 TAP 第二轮审查（Sol）

日期：2026-09-11  
性质：独立内部复核；根任务保留最终科学判断。本审查不构成期刊接收意见。

## 1. 第二轮裁决

**第一轮指出的窄数学假设、proper-complex Fisher 因子和受信核接口问题已基本修复；实际算法证书、证书决策价值、同成本优势和原创性仍未闭合。新完成的最强数值证据不是支持优势，而是明确支持强 dense/dual GN 基线在已测范围获胜。当前仍应大修，不得声称 A2 已有算法净优势、实用端口压缩优势或原创性闭合。**

本轮只核查指定的已完成记录。没有运行重实验，也没有读取或推断 imaging 的部分结果/检查点。`runs/a2/imaging/frozen_config.json` 与生成器只作为未执行协议审查。Source Pro 原附件仍缺失，这影响来源归档与原环境复跑，不构成本地推导的数学反例；36 项证书继续只能称为本地重构检查，不能称为原附件复跑。

## 2. 第一轮问题逐项复核

|第一轮问题|第二轮状态|复核结论|
|---|---|---|
|EM-GL-01 受信核结构|**已解决（声明范围内）**|`gamma_for_fft_kernel` 先要求精确 `FFTGreen` 类型，再按保存参数重新构造并逐项核对卷积核、正向 FFT 缓冲和伴随 FFT 缓冲。非对角核扰动和仅 FFT 缓冲扰动均返回不适用。该合同防意外构造漂移，不是 Python 安全边界；文本没有作恶意攻击主张。|
|EM-SCOPE-02 号约定/内积/离散|**已解决（声明范围内）**|理论和结果元数据现明确 `exp(-iωt)`、出射 `iH0^(1)/4`、`Im χ≥0`、均匀等面积方形源单元、中心配点、质量白化 Euclidean 坐标及自项定义。改变这些条件需重新推导。|
|INFO-FISHER-01 Fisher 假设|**已解决**|已限定 proper complex Gaussian/实化 Gaussian；非 Gaussian 情形只保留加权 LS 的 GN 曲率。若代码用复噪声 RMS `σ` 除残差后直接堆叠实部/虚部，则实分量协方差为 `I/2`，所以 `I_F=2 A_code^T A_code`；单位实协方差白化需再乘 `√2`。文本中的二倍因子正确。|
|NUM-GAP-01 局部 gap 与实际信赖域|**数学边界已解决；算法证书未检验**|文本已明确恒等式只覆盖同一固定、无约束、精确二次最优步，并把梯度缺失率降为启发式。矩阵自由 Jv/J*v 提供了完整局部 GN 对照，但尚无对实际裁剪、拒绝和重线性化轨迹的 KKT/模型差距证书。|
|NUM-ADV-02 净优势|**未解决，且现有证据反对该主张**|随机谱预条件显著减少部分 CG 迭代，但在 v1/v2 的完整已计成本下始终败给强 dense/dual GN；见第 4 节。|
|NUM-CERT-03 证书实用紧度|**未解决**|36 项重构检查仍全部覆盖，但 `bound_rel=0.01294--11.16135`、bound/corrected-error 比 `3.645--19.443`、最大 corrected relative error `0.87433`；没有预注册的停止、增秩或结构选择决策收益。|
|INFO-SOM-02 内部 B 的信息含义|**不存在待纠正的正面主张**|现理论明确说 B 是确定性内部量，不增加观测 Fisher 信息；它只可能改变选向、正则、预条件或非线性路径。不得虚构“B 增加信息”的主张再批评。|
|NOV-PRIOR-01 最近邻/原创性|**未解决**|标准随机 range、低秩谱 deflation、TSOM 内部传播选向和认证降阶均不能计为本轮原创性。claim chart、指定近邻技术对照及可证伪组合收益仍未闭合。附件缺档应与数学正确性分开记录。|

## 3. 耗散证书的数学与数值来源

修订后的窄命题自洽。对 `g=iH0^(1)/4`，参考虚部为 J0 Gram 核；对称正权张量 GL 的角向 form factor 在 `kh≤π` 时非负。正方形实际自项相对参考自项只产生对角 shift，因此

`-Im_H(X^{-1}-D) ⪰ diag(Imχ/|χ|² + shift)`。

号和材料项均正确；`χ=0` 返回需要先消元，`kh>π` 对积分核返回未知而不是宣称不稳定。证明只涉及声明的二维标量离散核，不覆盖连续模型、非均匀网格、三维 Maxwell 或浮点区间认证。

`runs/a2/quadrature_certificate/results.json` 的检查器哈希与当前 `check_quadrature_certificate.py` 一致，保存的 `physics.py` 构造器哈希也与当前文件一致。N8/N12、点/积分、`kh=0.2/1.0/2.9` 共 12 项的最小 shift 后余量为 `-1.0436e-14`；`kh=3.5` 和两种构造漂移均正确拒绝。它验证当前受信构造与代数的一致性，不是 Source Pro 原附件复跑。

## 4. 矩阵自由基线：正确，但优势主张被否定

小检查支持当前离散 Jv/J*v：伴随相对误差 `1.677e-12`，有限差分相对误差 `1.034e-10`，矩阵自由与稠密 GN 步相对误差 `1.541e-9`。这验证导数/伴随实现，不验证非线性成像成功。

完整成本结果如下；v2 的 dense 基线在 `p>2m` 时用 192 条实化 adjoint rows 组装 J，再用 dual Woodbury 求四个共享状态残差，避免了不必要的 `p³` 求解。

|记录|p|dense/dual GN (s)|plain CG (s)|最佳随机预条件 CG (s)|结论|
|---|---:|---:|---:|---:|---|
|v1|49|2.245|25.998|16.979（rank 12）|预条件优于 plain，但比 dense 慢约 7.56 倍。|
|v1|196|9.388|40.018|27.301（rank 12）|预条件把 capped plain 改为 62--67 次收敛，但仍比 dense 慢约 2.91 倍。|
|v2|784|2.047|39.669|10.115（rank 64）|rank 64 仅 7--8 次 CG，但总时间仍为 dense 的约 4.94 倍。|
|v2|3136|4.672|50.004|14.964（rank 64）|rank 64 为 13 次 CG，但总时间仍为 dense 的约 3.20 倍。|

v2 的 rank-12 方法甚至比 plain CG 更慢：p=784 为 `41.003` 对 `39.669 s`，p=3136 为 `51.880` 对 `50.004 s`。最佳 rank-64 在 p=784 的总 tangent+adjoint RHS 为 `792+824=1616`，dense 为 1536 adjoint RHS；p=3136 则为 `960+992=1952` 对 1536。被测 J 只占 1.20/4.82 MB，因而该范围也没有展示内存优势。

结论必须是：随机 range/deflation 在固定幅值、四个残差的局部共享状态上改善了迭代数和 capped-CG 精度，但没有算法净优势，也没有原创性。结果不能外推到非线性移动/分裂/合并 Gaussian 或端到端反演。

## 5. Fresnel measured 记录

这是一份单个二维 TM 圆柱采集上的模型失配诊断。偶/奇接收器来自同一采集且强相关，不是独立目标或几何外推。主拟合先由 incident back arc 估计每视角源系数，又在散射训练数据上逐频 profile 一个复增益；后者是附加 nuisance，不能称为只用 incident 数据完成的源归一化，也保留材料幅值歧义。

精确负结果应保留：

- gain-profiled 训练目标：1/2/4 Gaussian 分别为 `0.084982/0.080309/0.078513`，known-shape softened-disk oracle 为 `0.058028`。disk 在 2/4/6 GHz 的 held relative error 为 `0.09103/0.13776/0.17242`，均低于所有 Gaussian 结果。2/4 Gaussian 的十个起点没有一个以优化器成功状态结束；所选结果均到 45 次上限。1 Gaussian 虽有 7 个成功起点，最低目标记录本身也到上限。
- gain 固定为一的后验 refit：1/2/4 Gaussian 目标为 `0.174218/0.177243/0.171283`，disk 为 `0.130118`；disk 的 held error `0.19443/0.18576/0.23231` 仍逐频低于三种 Gaussian。disk 是 truth-shape oracle，不能作为公平通用算法基线，但它清楚显示当前 Gaussian 模型与校准并未闭合。
- 现有十起点只保存每模型的十个训练目标和所选状态，没有保存每起点的 held error、时间和参数；不得称为十组完整配对统计。相对拟合 disk 的边界图像差异只是 proxy，不是 Gaussian 逼近下界。

这里有一个真实的文本矛盾和一个真实的来源缺口：结果顶层 `scope` 仍写复增益“handles source normalization only”，与同文件正确的 `gain_profile_interpretation` 相冲突；主长跑的精确执行源码哈希未保存，且当前 `measured.py` 已在运行后修改。名为 `postrun_fixed_gain_source_sha256` 的字段只保存了当前 `measured.py` 哈希，没有保存 `fixed_gain_refits.py` 自身哈希。具体修复是：改正旧 `scope`；将主结果明确标记为执行源码谱系不可完全复核；以后同时冻结 runner、所有导入模块、环境和日志，不能用当前文件哈希回填成历史执行哈希。

## 6. N24 local-T 端口回归

当前 `ports.py` 的大网格分支已把 `T_i` 同时作用于激励 `E` 和每条相互作用边 `DQ_j`，与

`j_i=T_i(E+Σ_{k≠i}Dj_k)`

一致。N24 JSON 的三项完整分量恒等式误差为 `2.41e-9--7.79e-9`，完整 current 状态残差为 `3.98e-10`，支持实现修复。与此同时，9 维压缩网络的 scattered relative error 为 **`0.681169`**，是明确的端口压缩负结果，不能只报告恒等式通过。

该 JSON 没有生成脚本、频率、分量参数、rank/enrichment、容差、随机种子、环境或源哈希；仓库中也找不到重构它的入口。因此数值谱系未闭合。具体修复是保存一个确定性回归 runner 和完整配置/依赖哈希后原样重跑，保留 `0.681169`，不得只补一个无法证明历史来源的哈希。

## 7. 30-case imaging 协议审查（不审结果）

`frozen_config.json` 的 runner 哈希与当前生成器一致，且 N128 cell-integrated 数据/N64 point inverse 的离散失配设计是合理的开发压力测试。但该协议只能回答局部恢复问题：

1. 30 cases 是 3 个固定 morphology templates 各做 10 次平移、噪声和初始化扰动，不是 30 种形状；只有一个总种子序列。
2. 初始化由 `truth_theta.copy()` 产生，只扰动幅度、中心和尺度；模型始终固定为两个 Gaussian。优化函数中的 `truth` 只用于事后误差，不直接参与步计算，但 truth-informed warm start 已使它不是 blind recovery。
3. sharp template 的 `truth_theta` 实际只是两-Gaussian 初始化代理，不是 sharp truth 参数；`fit_parameter_tangent_som` 返回的 `parameter_scaled_error` 对该族没有真值意义。最终材料误差由 runner 用 sharp field 另算才有意义。
4. random-complement 每轮只使用一个固定随机实现，且所有 cases 共用按迭代编号生成的随机系数；它是 matched control，不是随机基线不确定性研究。held receivers 只是同一几何中的奇数相邻通道。
5. `adaptive_tsomg` 仍由梯度捕获率启发式控制，没有 residual/KKT/certificate 的端到端保护。配置只冻结 runner 哈希，未冻结显式 case manifest、`som.py`/`physics.py` 等依赖哈希和环境。

在执行前应把协议正式标为“三模板 truth-informed 两-Gaussian 局部恢复诊断”，冻结 30-case manifest 和依赖环境；sharp 族删除/标无意义的参数真值误差；给 random control 增加预先规定的多种子汇总。即使这 30 cases 完成，也不解除 Frozen900、diverse-real、端到端 certified-adaptive SOM 或 new3D 门槛；这些门槛仍是**未执行/未检验**。

## 8. TAP 硬门槛与逐项评分

|独立硬门槛|状态|依据|
|---|---|---|
|核心无未解决致命数学反例|**通过（仅限声明的有限维标量命题）**|Fisher、二次 gap、local-T 恒等式和 GL Gram 的条件与边界现已明确。|
|受信核发证接口闭合|**通过（理想离散、非区间意义）**|全核/FFT 缓冲重构核对及两项拒绝回归成立。|
|最近邻差异有证据|**未通过**|标准 randomized range/deflation 无原创性；联合贡献 claim chart 与直接近邻对照未完成。|
|冻结测试含完整成本优势|**未通过**|现有完整成本 benchmark 由 dense/dual GN 获胜。30-case development 不能替代 Frozen900。|
|主张对应全波/实测支持|**未检验**|已有全波局部 benchmark 和单圆柱诊断，但没有端到端 SOM 优势的 diverse-real/新三维证据。|
|数据、代码、图可复核|**未通过**|矩阵 benchmark runner 哈希可核对，但 measured 主执行源码缺档、N24 回归无生成入口，多数记录未保存完整依赖/环境哈希。|

|维度|分数（1--5）|理由|
|---|---:|---|
|正确性|3|窄数学和 Jv/J*v 实现已较好闭合；实际信赖域自适应仍无证书，数值来源尚有缺口。|
|原创性|2|未给出超过 TSOM、认证 ROM 和标准 randomized deflation 的可核查差异及收益。|
|实验|2|有强局部负基线和一个 measured 模型失配诊断，但无完成的冻结多场景、端到端或 diverse-real 验证。|
|实践价值|2|dense/dual 基线更快，N24 端口误差 0.681，证书尚未产生决策节省。|
|表达|3|理论边界清楚并保留负结果；measured scope 矛盾和若干来源字段仍妨碍独立复核。|

不计算平均分，也不据此作期刊接收承诺。

## 9. 应升级给 Pro 的三个核心科学瓶颈

1. **可扩展性/基线瓶颈**：推导并验证 dense-primal、dense-dual 与矩阵自由/谱预条件的 break-even 条件，显式包含实数据维数、参数数、照明数、重复残差/反演次数、J/J* RHS 和内存。当前低数据秩区间 dual GN 的胜出不是继续调 randomized rank 能解决的问题。
2. **端到端认证自适应瓶颈**：把局部正规残差、信赖域 KKT/实际-预测下降及状态/输出证书连接成可执行的增秩、拒绝或结构变更规则，预定义 false accept/reject 和工作量指标。否则 36 项覆盖与启发式 adaptive rank 仍是两套未连接工具。
3. **真实问题定义与原创性瓶颈**：在校准 nuisance、blind/model-order 变化、diverse targets 和成熟边界正则/TSOM/认证 ROM 对照下，形成逐 claim 最近邻表。唯一值得继续检验的贡献应收缩为“联合材料方向残差与状态证书是否在某个明确 break-even 区间产生净成本和错误结构控制收益”，而不是把标准 range、端口代数或内部 B 当成新信息。

根任务的最终科学判断应以这些硬门槛为准，而不是以局部检查通过或 CG 迭代数下降代替。
