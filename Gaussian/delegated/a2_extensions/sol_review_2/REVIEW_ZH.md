# A2 extensions / GPU 分支第二轮独立只读评审

## 范围与总评

本轮仅静态核查指定源码、首轮回复、现有 `*checks*.json`、GPU 协议及 `CORE_AUDIT_ZH.md` 第九节；未 SSH、未启动远端任务、未运行重型计算，也不把正在运行的 main128 当作完成证据。

总评：`manifold_revised.py` 已实质修复首轮 M1/M3/M4/M6/M7 中的代码性问题，M5 的导数公式静态核对也一致，但现有端到端检查仍绑定旧 `manifold_runner.py`，没有覆盖 revised 源码。v3 历史足以在共同 RHS 上限处构造预先规定的阶梯轨迹，因此无需强求各法最终恰好消耗同一 RHS；但首轮回复关于“保存每次有效状态图像”的表述不符合源码。GPU 的隐式复伴随和多频布局在公式与现有相对检查中自洽，真实残差也是硬失败条件；当前主要阻断是单频独立检查的源码哈希已过期，以及 run 配置未冻结 device/硬件，可能污染恢复与成本归属。表示比较只能解释为已声明的 bounded development comparison，不能推出一般表示优劣。以下不构成 TAP 接受。

## A. 实际代码/证据缺陷

### [高] G1：单频独立伴随/有限差分检查没有覆盖当前 `torch_vie.py`

`backend_checks.json` 与 `backend_checks_remote.json` 记录的 `torch_vie.py` 哈希均为 `46af...`，当前文件哈希为 `fcc1...`。当前哈希只由 `batch_checks*.json` 覆盖；而该检查是 `TorchMultiVIE` 对同文件内逐频 `TorchVIE` 的相对比较，不能独立排除两者共享的隐式梯度错误。

现有旧哈希检查确实包含：SciPy 场参考、FFT `D` 伴随内积和损失方向导数；数值也通过。但在当前源码重新产生同类检查前，不能把它写成“当前 GPU 后端已获独立梯度验证”。这是轻量重跑门槛，不要求干扰 main128。

### [高] G2：run 配置哈希漏掉 device/硬件，恢复隔离不完整

`representation_campaign.py:153-156` 的 `config_hash` 包含源码、网格、步数和学习率，却不含 `args.device`、GPU 型号或运行环境。相同 run-name 可在配置哈希不变时把 checkpoint 从 CUDA 映射到 CPU/另一 GPU继续运行；随后 `frozen_config.json` 还会被当前调用覆盖。这样累计秒数、求解迭代和数值轨迹可能来自不同设备，却仍被视为同一冻结运行。

main128 若始终只在同一 4060 上恢复，其现有结果不因此自动失效；但正式成本结论前必须由日志/产物证明未发生跨设备恢复。后续应把 device、dtype、求解容差、maxiter 与硬件身份纳入不可覆盖的运行清单，或明确禁止跨环境 resume。

### [中] G3：`actual_optimization_seconds` 混入 checkpoint I/O，且恢复边界的计时口径不严格

每一步都在计时区间内序列化完整模型、Adam 状态和不断增长的 trajectory；不同参数量导致 I/O 成本不同。字段名却称 `actual_optimization_seconds`。最终 checkpoint 自身的写入时间会进入本次最终秒数，但不会进入该 checkpoint 内保存的 `cumulative_seconds`，若恰在写完最终 checkpoint 后中断再恢复，累计口径会发生小幅变化。

该量可以作为“训练加逐步持久化”的端到端成本，但不能无说明地当纯 GPU 优化时间。RHS、批量 BiCGSTAB 迭代与端到端墙钟应分栏解释。

### [中] G4：所谓每 case `cuda_peak_bytes` 实际是进程累计峰值

主循环没有在 case/method 前调用 peak-memory reset；`torch.cuda.max_memory_allocated` 因而是自进程启动以来的峰值，不是该 case 或该方法的独立峰值。除非只报告全程全局峰值，否则该字段标签会误导表示间显存比较。

### [中] G5：`linear_iterations` 对正常收敛的每次 solve 多计 1

`solve_linear` 的 `for it in range(...)` 先更新循环变量，再在循环顶部检查 `active`。若第 k 次迭代末所有 RHS 收敛，下一次进入循环时 `it` 已变成 k+1，随后才 break；`last['iterations']` 和累计 `linear_iterations` 因而记录 k+1，尽管第 k+1 次没有执行 `A`。零 RHS 同样会记 1。该错误不影响解或真实残差门槛，且完整方法的 solve 数相同使相对差值大体保留，但会系统性抬高成本账本；应改用只在实际执行迭代后递增的计数器，并重新解释现有产物。

### [中] M8：首轮回复夸大了 v3 保存的状态内容

`manifold_revised.py:88` 每轮保存 RHS、训练损失、材料相对误差和留出场误差；这足以构造标量质量—RHS 阶梯轨迹。可是 `images.npz` 只在每个 case 完成后保存 truth 与四个方法的最终重建，并未保存“每次有效状态的图像”。因此可以做共同上限标量比较，不能重建任意预算处的材料图。

### [低] M9：若将来修改 `parameter_weights`，当前配置项不会生效

v3 冻结值为 `[1,1,1]`，与源码当前未加权欧氏乘积度量一致，故不影响本次结果。但该字段没有进入 Jacobian、步长或范数计算；将来把它改成非 1 值只会改配置文字，不会改算法。正式开放该参数前应实现映射或删去可配置假象。

### [中] M10：现有端到端场导数检查仍是旧 runner 的证据

`runs/a2/extensions/manifold_native_checks.json` 的 `field_fd_max≈2.0e-8` 数值良好，但其源码账本是 `check_manifold_native.py -> manifold_runner.py`（哈希 `58eb...`）；检查脚本也确实从旧 runner 导入 `material/retract/fit`。它没有哈希或导入 `manifold_revised.py`（当前 `8e7f...`）。因此可以静态确认 revised 的 12 列公式与旧实现同型，不能接受回复中“v3 端到端检查已补足”的证据表述；需让轻量检查直接绑定 revised 源码后才闭合 M5。

## B. 已核对通过的实现点

### [通过] M1/M3/M4：补空间语义、随机基线和有限步组合已修正

- 仅当 `rank < len(s)` 时探测 `Vh[rank]`；满秩明确记录 `no_weak_complement`。
- 随机方向先投影到同一 `V` 的正交补并断言 `||V^T h||<1e-8`。
- `h` 固定为范数 0.5，正交可见修正单独限制到总半径 0.7，不再整体裁剪后冒充原 probe 的 profile。
- 接受条件加入尺度相关的最小下降；方法仍应称“阻尼 GN 提案 + 全波单调验收”，不是标准 trust-region ratio 算法。

### [通过] M2：v3 历史支持共同 RHS 上限分析

每个 history 节点同时记录实际累计 `rhs` 和三项质量量。普通法现有终值为 3600 RHS，secant 法因保守 guard 终值约为 3480--3540 RHS；这不要求人为浪费计算补齐。合法分析是对每 case 先取四法共同可达 cap，再使用 cap 之前最后一个已产生状态的阶梯值/完整轨迹，且预先固定评价量，不能用 truth 或 held-out 误差事后挑选最佳迭代。墙钟没有逐节点历史，因此目前只能构造 RHS 前沿，不能构造严格的时间前沿。

### [通过/证据待补] M5/M6/M7：公式、参数尺度与方法命名边界基本闭合

当前 revised 源码哈希与 v3 frozen config 一致；峰值、中心 0.03 m 尺度和 SPD Frobenius 正交基的材料导数、retraction 及 VIE tangent 组合经静态逐式核对一致。配置冻结了尺度和全 1 权重；谱法应继续称阈值截断 SVD，而非 continuation。`full_lm` 代码标签本身仍略强，论文文字需用实际算法名称。端到端数值证据仍须按 M10 直接绑定 revised 文件。

### [通过] G6：隐式复伴随公式正确

状态方程为 `A(chi)j=diag(chi)E`、`A=I-diag(chi)D`。令 `lambda=A^{-H}S^H grad`，则对共享材料的反传为各发射求和的 `conj(E+Dj)*lambda`，与 `_Implicit.backward` 一致。`A(...,adj=True)` 实现 `I-D^H diag(conj(chi))`；`D` 的裁剪 Toeplitz 伴随使用 `FFT(conj(kernel))`，不是错误地直接共轭未居中 FFT。这里的“通过”是公式/静态实现判断；当前哈希的独立数值证据仍受 G1 限制。

### [通过] G7：多频 batch 布局保持频率独立且次序一致

列布局是 `[frequency, transmitter]`：`E` 按频率拼列，`D` 重排为 `(nf,nt,n,n)`，接收输出为 `(nf,nrx,nt)`，接收伴随再还原为 `N x (nf*nt)`。没有跨频卷积或材料之外的频间耦合。当前源码哈希匹配的 CPU/CUDA batch checks 对逐频场与梯度均通过；该检查验证布局与相对一致性，不替代 G1 的独立梯度证书。

### [通过/硬失败设计] G8：真实残差不是固定传播步的替代品

BiCGSTAB 用递推残差决定 active 集，但返回前重新计算 `b-Ax`，非有限或超过 `1.5*tol` 即抛错；复 64 默认 `tol=2e-5`，与协议一致。成功完成的 production 方法因此意味着每次 solve 都通过运行时门槛。不过结果文件只保留最后一次 solve 的 `last` 于检查脚本，campaign 不保存逐 solve 残差裕量；它支持“门槛已执行”，不支持残差分布或物理模型误差证书。

## C. 已披露限制，不另判为代码 bug

### [限制/高解释风险] R1：表示比较不是一般公平性或表示优越性证明

Gaussian K=16/64/144 与 voxel 共享数据、训练/留出接收机、Adam 步数、学习率及渲染后平滑罚项，但参数映射、有效先验、梯度尺度和初始化不同；voxel 从 K64 初始图开始，而不同 K 的初始材料图不同。协议和 config 已明确披露这些差异。结果最多支持该冻结开发设置下的 quality/cost 点或前沿，不支持“Gaussian/voxel 一般更优”。

GPU trajectory 每步只保存训练 `loss/fit/smooth`。由于每个完整 step 固定一次 36-RHS 正演与一次 36-RHS 伴随，可重建训练 fit 对 RHS 的轨迹；但材料 L2、held-out、IoU 只在终点计算，不能从历史构造这些外部质量量的中途前沿。若四法都完整达到共同 80 steps，可在共同更新数/RHS 点比较终值，同时另报不同 BiCGSTAB 迭代与墙钟；若有方法未完成，不得拿其部分终值与完整终值直接排序。

### [限制] R2：模型与证据范围保持原边界

GPU 是二维标量 VIE、无额外测量噪声且保留 N192/N128 离散失配的开发比较，不是 SOM 测试。求解残差不覆盖离散化、端口或模型误差。已有文献材料只构成有限的 primary-source 邻近检查，不能支撑“首个”、完整 prior-art 排除或原创性裁决。

## D. `CORE_AUDIT_ZH.md` 第九节

### [通过] 有限 posterior 事件概率上下界正确

由 `||F_z-Ftilde_z||<=epsilon_z`，平方范数差满足

`|ell_z-elltilde_z| <= (2||y-Ftilde_z|| epsilon_z + epsilon_z^2)/sigma^2 = d_z`。

所以每个未归一化后验权重位于 `[w_z^-,w_z^+]`。事件概率 `A/(A+B)` 对事件内权重单调增、对事件外权重单调减，文中以事件内下界/事件外上界给下界、反向给上界，公式正确。其成立条件是同一有限状态族、同一 prior、标量 proper-complex Gaussian 噪声约定及逐状态有效的场误差界。

硬门槛保持：`epsilon_z` 必须覆盖求解、离散、端口和物理模型误差；不能用少数状态证书替代全状态/未覆盖质量控制；小噪声会按 `sigma^-2` 放大界；该有限枚举不等于连续场或 NN 后验认证，也不是新贝叶斯定理。

## E. 结果接纳前的硬门槛

1. main128 同步回本地后，逐文件核对 frozen config、`representation_campaign.py`、`torch_vie.py`、`physics.py` 哈希，并确认所有计划 case/method 完成且无 solver exception；远端仍在运行、checkpoint 存在或脚本退出均不等于 TAP 接受。
2. 对当前 `fcc1...` 版 `torch_vie.py` 重新生成单频 SciPy 场、`D` 伴随和独立方向导数检查；多频相对检查不能替代此项。
3. 让小网格 12 列材料/全波场有限差分检查直接导入并哈希 `manifold_revised.py`；旧 runner 的良好误差不能作为 revised 的可执行证书。
4. 证明本次 main run 未跨 device/硬件恢复；否则时间与迭代成本分开作废或分段报告。后续修复 G2 后再建立正式 resume 证据。
5. 流形结论必须按共同 RHS cap 的预定阶梯规则汇总，并保留“条件离散开发比较”措辞；不得从 secant 轨迹推出 chart 失灵、全局优越或理论原创。
6. 表示结论只报告冻结设置下各 family 的终点/可合法重建的训练前沿、实际求解迭代和端到端成本；不得用未保存的中途 held-out/material 指标制造前沿，也不得外推一般 Gaussian 与 voxel 的优劣。

本评审不作 TAP 接受声明。
