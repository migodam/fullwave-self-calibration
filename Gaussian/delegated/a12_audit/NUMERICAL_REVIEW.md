# A12 local-T 数值与公式独立复核

日期：2026-09-11  
范围：只读核查 `Theory/a12/CORE_AUDIT_ZH.md`、`code/a12_theory_checks.py` 与 `runs/a12_theory/results.json`；不改代码，不训练网络，不作新颖性裁决。

## 结论

这组结果可靠地支持一个**有限维构造性结论**：在同一个 16×16 scalar VIE 离散中，把总 contrast 拆成三个可重叠 Gaussian component 后，完整 256 维 local `T_i` 的 component equation 与原始全局 VIE 在机器精度内一致。它还显示：只按原始 `T_i` 的 Frobenius/SVD 能量截断，并不善于保留当前 illumination 与相互作用真正需要的方向；用 plane-wave excitation response 构造输出 range 更有效，而加入两轮 interaction snapshots 可在足够 rank 时突破该 basis 的误差平台。

它不能支持 finite-port full-wave exactness、普适 port-rank law、重叠或强 contrast 定理、端到端加速、Hermite 阶与 multipole 阶对应、GNN 优势或新颖性。尤其“两轮 enrichment 后下降”是这个离散实例、snapshot 集与人为 weighting 下的 basis-design 结果，不是新定理。

## 公式对齐

### 1. 完整 local `T_i` 分解

令 `X_i=diag(xi_i)`，总 current 为 `j=sum_i j_i`，且 `j_i=X_i(E+D j)`。移去本 component 的 self contribution 后：

```text
(I-X_i D) j_i = X_i (E + D sum_{k!=i} j_k)
T_i = (I-X_i D)^(-1) X_i
j_i = T_i (E + D sum_{k!=i} j_k)
```

代码第 31–34 行正是这样构造：`comp[i]=xi_i*E_total`，`T_i=solve(I-X_iD,X_i)`，再检查 `comp[i]-T_i@(E+D@(j-comp[i]))`。四个 case 的最大相对误差为 `0.88e-15` 到 `1.17e-15`。这验证的是**同一离散、同一 Green matrix、同一 reference solve 内的代数恒等式**。

该检查没有验证 continuum、独立 solver、有限 ports、VSWF translation、局部被动 scattering matrix 或截断误差。Gaussian component 的无限尾和 overlap 在 additive-potential algebra 中允许存在，但这不把它变成独立非重叠散射体。

### 2. 三种有限端口近似

代码第 37–52 行把每个完整 `T_i` 写成 `W_i V_i`，并解 reduced component network：

```text
z_i - sum_{k!=i} V_i D W_k z_k = V_i E
j_hat = sum_i W_i z_i
```

三种构造实际为：

| JSON 名称 | 代码实现 | 精确含义 | 主要限制 |
|---|---|---|---|
| `raw_operator_svd` | `W=U_r Sigma_r`, `V=Vh_r` | `T_i` 的 rank-r 最佳 Euclidean/Frobenius operator approximation | 优化所有 256 维 input directions，不针对当前 excitations、interaction 或 receivers；99% `T_i` 能量不保证 network/output 精度。 |
| `plane_wave_response_range` | `Q_i=left_svd(T_i E_b)_r`, `W=Q_i`, `V=Q_i^H T_i` | 左投影 `T_i ~= Q_i Q_i^H T_i`；输出空间由 64 个同频 plane-wave excitations 训练 | 它不是 receiver-observable range，也没有使用 `S`；名称更准确应是 plane-wave-excitation response range。结果依赖角度采样、snapshot weighting 和同频率。 |
| `interaction_enriched_range` | 从 plane-wave `Q_i` 开始，两轮加入 `T_i D Q_j`，每轮重新取 rank-r left SVD | 针对 component-to-component induced fields 的两步 output-range iteration | 两轮是固定 heuristic；每个 snapshot block 被各自 Frobenius-normalize，外部与内部权重人为；不是 TSOM second fold，也不是收敛 Krylov/error-bound 方法。 |

interaction 更新可明确写成：

```text
Q_i^(m+1) = rank-r left range of
  [ normalize(T_i E_b), normalize(T_i D Q_1^m), ..., normalize(T_i D Q_k^m) ]
```

其中省去 `j=i`。实现是同步更新，代数无明显索引错误。三种方法均先显式形成完整 256×256 `T_i`；因此它们目前是 oracle compression diagnostics，不是避免 full local solve 的实际算法。

## JSON 结果复算与解释

原结果中，三个 local `T_i` 的 99% Frobenius-energy rank 为 20–22。下表的 `plain16`、`enrich16` 和 `raw64` 均为 scattered-field relative error；`state16` 为 enrichment rank 16 的 full VIE state-equation residual。每个 component 16 ports 时总 network dimension 为 48。

| 第一对中心距 | amplitude scale | `plain16` | `enrich16` | `state16` | `raw64` |
|---:|---:|---:|---:|---:|---:|
| 0.30 | 0.5 | 0.273% | 0.150% | 0.418% | 0.0688% |
| 0.30 | 2.0 | 1.191% | 0.231% | 0.565% | 0.0824% |
| 0.10 | 0.5 | 1.884% | 0.543% | 1.094% | 0.0680% |
| 0.10 | 2.0 | 8.540% | 0.987% | 1.718% | 0.0921% |

可以支持的具体观察：

- plane-wave range 在 rank 8 到 16 基本平台化：四个 case 的 `field_error(rank16)/field_error(rank8)` 为 0.927、1.034、1.003、1.005。
- rank 16 时，两轮 enrichment 在四个 case 都降低 receiver error；相对 plain range 的 error ratio 为 0.550、0.194、0.288、0.116。改善在 amplitude 2.0 和高 overlap case 最大。
- 改善不是 rank-uniform。rank 8 时，两个 amplitude 0.5 case 反而恶化 153% 和 9.6%；所以不能写成 “interaction enrichment 总是更优”。
- receiver error 可明显小于 current error 与 state residual。最强 overlap/amplitude case 的 receiver error 已到 0.987%，state residual 仍为 1.718%。这支持 CORE_AUDIT 要求 full-state closure，不能只看观测误差。
- raw `T_i` 的 99% energy rank 约 21，但 raw rank 32 的 receiver error仍为 2.16%–3.01%，rank 64 才到约 0.07%–0.09%。这直接说明局部 operator Frobenius-energy threshold 不是 network accuracy certificate。
- raw rank 64 的总 component dimension 为 192，已接近原始 total-current dimension 256；这四个 case 不能证明显著的通用压缩率。

### “overlap”和“strong contrast”的实际范围

我用项目已有 scipy venv 重新构造四个 case，只计算材料幅度、离散条件数和谱半径；没有重跑训练或改产物：

| 距离 | scale | 最大 `abs(chi)` | 最大 component Gram correlation | `cond(I-XD)` | 最大 local `rho(X_iD)` | component `rho(TU)` |
|---:|---:|---:|---:|---:|---:|---:|
| 0.30 | 0.5 | 0.465 | 0.062 | 2.34 | 0.312 | 0.076 |
| 0.30 | 2.0 | 1.859 | 0.062 | 11.17 | 1.248 | 0.229 |
| 0.10 | 0.5 | 0.678 | 0.641 | 3.38 | 0.312 | 0.208 |
| 0.10 | 2.0 | 2.711 | 0.641 | 20.85 | 1.248 | 0.592 |

这里 `cond(I-XD)` 是全局 state matrix 条件数；local block 最大条件数另算为约 2.22 或 10.38，与表中同量级。距离 0.10 的第一对 Gaussian 的离散 correlation 0.641，确实是较强 overlap；距离 0.30 的 correlation 只有 0.062。scale 2.0 比 0.5 更强，但四个全局系统的条件数最高约 21，且 component iteration 的 `rho(TU)` 仍小于 1；没有覆盖 collective resonance、极高 contrast、大电尺寸或 continuum refinement。

两个 scale 2.0 case 的 local `rho(X_iD)` 约 1.248。这是一个有用的独立确认：local inverse 存在，但原文写出的 local Neumann series不收敛。当前直接 `solve` 的 full-space identity仍成立，恰好验证了“resolvent 存在”和“Neumann 路径展开收敛”必须分开。

## 计时边界审计

现有 timing 不允许比较三种 port 方法的端到端速度：

1. `local_setup_svd_seconds` 从代码第 31 行之后开始，包含三个完整 256×256 `T_i` 的 dense solve 与 full SVD，约 0.030–0.034 s；它不含 VIE reference forward、`D` 构造或 geometry。
2. plane-wave 的 `T_i@E_b` 和三个 response SVD 在第 36 行、所有 rank timing 之外完成；其主 basis-construction cost 没有记录。
3. `raw_operator_svd` 与 `plane_wave_response_range` 的 `assembly_solve_seconds` 基本只计 reduced assembly/solve，约 0.2–3.2 ms；sub-ms 单次 wall time 未 warm up、未重复，不能稳定排名。
4. `interaction_enriched_range` 的计时包含每个 rank 从头做两轮 snapshot multiplication/SVD，约 16–28 ms；它与 plain range 的 timing scope 不同，而且没有利用跨 rank cache。
5. 三种方法全都依赖已经显式构造的完整 `T_i`。即使 reduced solve 快，也没有展示怎样以低成本得到实际 local ports。
6. JSON 的全脚本 `wall_seconds=0.915` 在画图和文件写出之前取值，并混合了其余公式 checks；它不是 forward 或 compression benchmark。

可报告的只有“在这次运行中，各代码块耗时为何值”。任何 speedup、复杂度优势或 GNN training throughput claim 都需要统一 setup/solve/adjoint/update 边界、warmup/repeats、hardware/software manifest 和 full-solver matched-error baseline。

## 能支持与不能支持的结论

### 当前可支持

1. 三个 overlapping additive potentials 的 full-space component `T_i` equation 在这个 N=16 scalar discretization 上与 global VIE algebra 一致。
2. task/excitation-aware output range 可在低 rank 下远好于 raw full-operator Frobenius SVD；应按目标误差设计 port basis。
3. 单纯 plane-wave excitation range 会遗漏 multiple-scattering induced directions；两轮 `T_i D Q_j` snapshots 在这四个 case 的 rank 16 上降低了误差平台。
4. overlap 和 amplitude 增加时 plain basis 的平台明显升高；内部 snapshots 在这些实例中有补偿作用。
5. receiver fit、current error和 full-state residual是不同 gate；局部 singular-value energy也不是 network error certificate。

### 当前不能支持

1. “finite ports 严格包含任意阶 full-wave multiple scattering”或任何 continuum truncation theorem。
2. 给定 `k sigma`、contrast、separation 即可得到通用 port rank；这里只固定一个 sigma、一个频率、一个 grid、三个 blobs 与六个 point-source states。
3. enrichment 对任意 overlap/strong contrast 单调更优；低 rank 已有反例，且最强 case仍远离已证 collective resonance。
4. interaction enrichment 是新算法/新定理；它是由完整 `T_i` 驱动、固定两轮、任意 block normalization 的 snapshot-range heuristic，与 block Krylov/reduced-basis 邻域尚未比较。
5. 48-port 网络已实现 universal compression；当前结果只说明 receiver error在四个 case 可低于 1%，state residual最高仍 1.72%。
6. 现有 timing 支持计算加速；basis construction和 full local solve没有公平计入。
7. 该 scalar 2-D grid 结果支持 Maxwell、VSWF、passivity/reciprocity、Hermite-to-multipole mapping 或 NN/GNN physics。

## 最小修复请求

- 把三种方法的对象写清：raw `T_i` operator approximation、plane-wave **excitation-response** range、interaction-enriched response range；不要称后二者 receiver-observable range。
- 统一计时：`D/T_i` 或 local training data构造、basis build、reduced assembly、solve、adjoint与更新全部分别计时，并为每法报告 shared 与 method-specific cost。
- 端口闭合曲线至少扫 grid 16→32、频率、sigma、节点数、随机几何、contrast、loss、overlap/gap与 near-resonance；用独立 finer-grid/full solver作 reference。
- 对 enrichment 做 0/1/2/4 轮和 snapshot-weight ablation；与 block Krylov/POD、randomized range finder、balanced truncation或可用 reduced-basis基线比较。
- port selection 用 full-state/data error或可证 residual estimator，不用 raw `T_i` 99% Frobenius energy单独截断。
- 把 true held-out excitations/frequencies/geometries 从 basis construction 中隔离；当前同频 64 plane waves是任务相关 dictionary，不是泛化证据。
- 报告 total-current dimension 256、component full dimension 768 与 reduced dimension `3r`，避免只选有利分母描述 compression。

## 可追溯性

- 当前 `a12_theory_checks.py` SHA-256：`8b22ba8fca748104b87fa68ae24d24d076d8e80130260f4affd3b5719cc7b84d`。
- `results.json` 内嵌 source hash 与当前代码完全匹配。
- 当前 `results.json` SHA-256：`bd6f5750ccb634caa5a801c3e2928cef159e2f3b8b3a73089c0e21dce9c3570d`。
- JSON 可解析；代码 AST 可解析；所有四个 exact component residual 已逐项读取；三种 rank curves 共 72 行结果，字段完整。
- 独立小检查使用 `research/trispace_self_calibration/a3_research/.venv3d/bin/python`，只重建四个 N=16 case 并计算 correlation、condition number和谱半径；未写入原结果目录。

最终状态：**bounded discrete algebra PASS；finite-port/continuum、universal rank、speed、GNN、novelty 全部 OPEN。**
