# A2 运行与证据索引

本目录属于开发和审计，不统一贴成“冻结测试通过”。A1_2 历史结果仍在相邻 a12_* 目录，未覆盖；本索引中所有路径相对 Gaussian/。

|编号/目录|生成入口|实际证据|范围或缺口|
|---|---|---|---|
|a2/source|Theory/a2/sources|Pro 两版完整分析正文|附件原字节缺失；raw CONVERSATION.json 仅本地归档，不入ZIP|
|a2/certificate|code/a2/check_certificate.py|36项修正场覆盖、嵌套/Schur检查|依据正文独立重建；声明离散模型|
|a2/parameter_geometry|code/a2/check_parameter_geometry.py|30项局部二次/尺度检查、交叉耦合反例|固定无约束二次步，不是整个信赖域轨迹|
|a2/quadrature_certificate|code/a2/check_quadrature_certificate.py|12项积分核条件及不适用拒绝|kh充分条件、相同构造，不是区间认证|
|a2/vector_certificate|code/a2/vector_certificate_probe.py|16项三维向量点偶极模型覆盖|不是3D材料反演或连续Maxwell验收|
|a2/ports|run_ports_probe.py、run_quadrature_ports.py、run_ports_n24.py|完整组件、导数/伴随、网格/圆柱、端口压缩|包含压缩失败；N24已新增带配置runner重跑|
|a2/som|code/a2/som.py|早期局部参数SOM与选择诊断|frozen_campaign_call_ledger是预算预估，900任务未执行|
|a2/imaging|code/a2/run_imaging_campaign.py|三模板×十扰动×五方法，逐方法checkpoint和完整轨迹|真值附近初值、固定双Gaussian、开发集；尖边proxy误差无效|
|a2/explicit_current|code/a2/explicit_current.py|早期N32十例显式current分支|状态闭合失败，不充当成熟TSOM|
|a2/explicit_current_v2|code/a2/explicit_current_v2.py|修正D Pweak的三个代表例、rank扫描|rank128仍有11.8%–13.2%状态残差|
|a2/matrix_free|check_matrix_free.py、benchmark_matrix_free.py|导数和GN步一致性；p49/196局部成本|幅度子问题、4共享状态残差|
|a2/matrix_free_v2|benchmark_matrix_free_v2.py|p784/3136，rank12/32/64与强对偶GN|计建基，无端到端非线性优势|
|a2/finite_catalog|code/a2/finite_catalog.py|有限类、30噪声记录、候选误差界与原场数组|invalid_v1排除；未证明防假阳性收益|
|a2/psf|code/a2/local_psf.py|有限切空间投影扰动的局部响应|非全局/任意点PSF|
|a2/measured|measured.py、fixed_gain_refits.py|Fresnel单圆柱实际拟合、图、增益诊断|主runner历史hash缺失；十起点无完整held/time|
|a2/sol_review、sol_review_2|code/codex_stage.py|Sol隔离审查，结构化最终返回与状态|第二轮未审后完成的成像等新项；根任务负责最终裁决|
|a2/package_validation|code/a2/validate_package.py|ZIP CRC、字节、原始来源与入口链接|文件完整性不等于科学复现|
|a2/portable_recheck|独立解压副本快速重跑|五项检查在交付副本重新执行|不重跑长成像、实测或网络训练|

## 来源与冻结状态

input_manifest.json 记录本轮开始时的原始输入哈希。局部runner中的source_sha256、冻结配置、checkpoint各自保留其记录时点。后置 postrun_snapshot.json 只记录封包前当前代码、依赖环境和结果字节，不能证明所有长跑都使用了当前源码。

尤其 Fresnel 主长跑的 executed_primary_source_hash 明确为 null；N24 历史无manifest结果隔离为 legacy_unmanifested，并用新runner重跑；有限目录错误实现隔离为 invalid_v1。这些问题不通过补一个当前hash来掩盖。

成像最初单例预计58分钟，未满足根任务原30分钟调度gate；随后根任务批准同配置继续完整运行。这个时间gate不是科学成功条件，也不允许改写成一开始就通过。随机补充分支只有一条固定控制随机序列，GSVD另改变系数度量，LM计时含额外诊断，均需保留。


## 本轮 A2 扩展（与前表分开）

|目录|入口|已执行内容|证据边界|
|---|---|---|---|
|extensions/source_manifest.json|两段ChatGPT全文归档|三个assistant正文与user提问、字节哈希|对话引用占位符不算已核查文献|
|extensions/theory_checks.json|extensions/theory_checks.py|导数、保矩、联合律、Neumann修正、静态似然、方向错配|有限复矩阵，不是实际Maxwell矩证书|
|extensions/manifold_native_checks.json|extensions/check_manifold_native.py|12参数×2频率场导数与更新机械检查|小网格；与v3同一SPD材料映射|
|extensions/manifold_v2|extensions/manifold_runner.py|30×4开发反演|满秩弱探测语义和随机基线问题，不作优势证据|
|extensions/manifold_v3|extensions/manifold_revised.py|新30×4修订反演、RHS/误差轨迹|同补空间、公用预算上限；未胜完整GN的整体优势|
|extensions/probability/v2|extensions/probability_v2.py|30次先验抽样、N64数据/N32有限字典|25种不同标签；近点质量posterior；带噪held|
|extensions/probability/weak|extensions/probability_weak.py|新30次抽样、仅2复通道|26种不同标签；post-hoc，非原冻结确认|
|extensions/probability/v2/optimization_posthoc|extensions/probability_v2_optimization.py|原v2数据的多初始化|优化诊断，不替换原test结果|
|extensions/rbf_capacity.json|extensions/check_rbf_capacity.py|118/256标签可由固定RBF符号模式表达|有限特征族，不是所有Gaussian表达上限|
|gpu/backend_checks_remote.json、batch_checks_remote.json|gpu/check_torch_vie.py、check_multifrequency.py|4060独立场/伴随/梯度检查及多频布局检查|当前源码哈希匹配；小型测试|
|gpu/scaling_remote.json|gpu/scaling_probe.py|N64/N128/N256全波正演+伴随|单频已知材料、非成像优势|
|gpu/representation_campaign/main128|gpu/representation_campaign.py|N128反演/N192积分数据、3频12Tx64Rx、K16/64/144与voxel|完整性状态以严格SUMMARY与execution_receipt为准；无测量噪声开发比较|

新入口：Theory/a2/extensions/CORE_AUDIT_ZH.md、COMPLEXITY_ZH.md、protocols/a2/EXTENSIONS_PROTOCOL_ZH.md。所有新增分支必须保留frozen development/post-hoc/运行状态区别，不能继承上阶段实测和向量DDA作为新方法验收。

最终GPU新增：90个Gaussian完整几何参数在`gpu/representation_campaign/main128/gaussian_parameters_posthoc.json`；120次N256固定解场检查在`refined_field_posthoc.json`与`REFINEMENT_SUMMARY.json`。后者为明确post-hoc，不替代冻结确认。主运行全部30×4已完成，跨重启边界见执行记录。
