# 从这里读：Gaussian A1＋A1_2＋A2

本包取代尚未发送的旧问题包。**A2 已保存 Pro 两版完整分析正文，并完成新一轮理论、算法、成像和实测诊断；尚未形成通过硬门槛的 TAP 方法论文。** 原始旧Pro附件字节仍缺失，本地代码是独立重建。新增概率与流形对话的完整文字已经纳入 A2；不要只读旧三问题版本。

1. [中文研究报告](RESEARCH_REPORT_ZH.md)：成像任务、实际图片、成本、负结果和判断。
2. [四个 Pro 核心问题](GPT_PRO_QUESTIONS_ZH.md)：决定主贡献能否成立的完整任务。
3. 原始 [A1](../Theory/%20A1.md)、[A1_2](../Theory/A1_2.md)、新 [A2](../Theory/A2.md)；再读 A2 入口下的推导与原 Pro 正文。
4. [路线比较](GAUSSIAN_VS_CALIBRATION_ZH.md)、[研究稿](PAPER_WORKING_DRAFT_ZH.md)、[最终裁决](../reviews/a2/FINAL_VERDICT_ZH.md)。
5. [运行索引](../runs/a2/RESULTS_INDEX_ZH.md) 与 [数据来源](../data/a2/DATA_INDEX_ZH.md)。

## 复核与复现

ZIP 解压后的根目录 START_HERE.md 指向同一入口；图使用包内相对路径。根目录 MANIFEST.json 校验文件字节，不代表预注册或科学通过。需要避免覆盖证据时，先复制解压目录再运行。

A2 的 CPU 代码只需 Python、NumPy、SciPy、Matplotlib，见 [依赖](requirements_a2.txt)。以下为快速数值复核，不需要公开数据、训练或付费 API：

```sh
python Gaussian/code/a2/check_certificate.py
python Gaussian/code/a2/check_parameter_geometry.py
python Gaussian/code/a2/check_quadrature_certificate.py
python Gaussian/code/a2/check_matrix_free.py
python Gaussian/code/a2/vector_certificate_probe.py
python Gaussian/code/a2/extensions/theory_checks.py
python Gaussian/code/a2/extensions/check_manifold_native.py
python Gaussian/code/a2/extensions/check_rbf_capacity.py
```

这些命令实际重新求解小型问题并写结果，不能在唯一证据目录中反复覆盖。首次安装依赖所需网络由使用者自行选择；本轮使用既有本地环境。`code/a12_graph/core.py` 随包提供，是第一条检查所需的轻量物理依赖，不需要 PyTorch。

完整成像是三模板各十扰动、五方法，已经保存逐 case-method 结果和轨迹。重跑前，在复制目录中移走 `Gaussian/runs/a2/imaging`，然后运行：

```sh
python Gaussian/code/a2/run_imaging_campaign.py
python Gaussian/code/a2/benchmark_matrix_free_v2.py
python Gaussian/code/a2/run_ports_n24.py
python Gaussian/code/a2/finite_catalog.py
```

成像 runner 的 resume 会跳过已有 case-method；跳过不是重新复现。成本受硬件和同时运行任务影响，完整耗时、求解次数与记录边界见报告。后置依赖清单不能充当所有旧长跑的历史执行哈希。

Fresnel 原始数据不在 ZIP 中。按数据索引从官方渠道获得对应文件并核对 SHA 后，才可运行 `measured.py`；如果官方入口暂不可用，不能把合成场当成实测复现。包内保存实际拟合参数与图，但原主长跑缺少精确源码快照与完整每起点 held/time 统计，这个缺口未被封包掩盖。

A1_2 的网络训练属于上轮证据，本轮没有重训新网络。旧完整包仍在本地归档；本 A2 包包含相关审计与紧凑历史结果，主要复现入口以上述 A2 代码为准。

本次交付复核：独立解压副本中的上述八组快速检查按最终交付副本执行，结果见 portable_recheck 记录；最终包核对全部Python源码与该副本一致。解压输入包的SHA保存在portable_recheck记录中，最终ZIP因追加该复核记录而具有不同SHA。长成像和实测没有在封包时重复运行。

## 新增概率、流形与4060规模实验

[扩展理论](../Theory/a2/extensions/CORE_AUDIT_ZH.md)和[扩展协议](../protocols/a2/EXTENSIONS_PROTOCOL_ZH.md)对应本轮新增来源。概率v2是冻结有限字典开发；weak、容量LP和多初始化明确是post-hoc。流形v2只保留开发诊断，v3纠正了补空间和预算口径。原A2的实测与3D有限模型不自动成为这两条新方法的验证。

```sh
python Gaussian/code/a2/extensions/manifold_revised.py
python Gaussian/code/a2/extensions/probability_run.py
python Gaussian/code/a2/extensions/probability_v2.py
python Gaussian/code/a2/extensions/probability_weak.py
```

概率第一条会重新构建N32的256状态缓存，v2读取它并生成N64数据。请在副本运行，避免覆盖原证据。

可选GPU后端需要PyTorch，CPU快速检查不要求安装GPU依赖。远端环境与实验范围见[GPU说明](../protocols/a2/gpu/REMOTE_SCOPE_ZH.md)。以下命令均可在包的副本根目录执行，`--device cpu`只适用于有相应资源的小测试：

```sh
python Gaussian/code/a2/gpu/check_torch_vie.py
python Gaussian/code/a2/gpu/check_multifrequency.py
python Gaussian/code/a2/gpu/representation_campaign.py --n 128 --steps 120 --run-name reproduction128
```

大型运行保存逐步checkpoint；续跑必须保持源码与配置哈希一致，且检查点必须可读取。系统重启后本轮出现过损坏文件，已隔离并重算未完成方法。ZIP不包含所有优化器.pt文件；完整Gaussian幅度、中心、协方差另存为 `runs/a2/gpu/representation_campaign/main128/gaussian_parameters_posthoc.json`，供Pro直接检查。复现结果和最初执行结果应分开保存。任何GPU执行状态都不表示论文科学门槛已经通过。
