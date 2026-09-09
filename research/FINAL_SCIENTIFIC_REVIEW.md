# 最终科学审查：位姿混淆下全波逆散射 SLAM 的谱几何

日期：2026-09-04  
状态：**论文初稿已形成；自动研究回路已锁定；有限维理论验证已复跑；尚不构成连续理论、真实系统或发表级验收。**

## 1. 交付结论

本轮研究形成了一条可辩护、可执行、且经过反例收窄的论文主线：在白化、realify 后的数据空间中，将材料地图切向量记为

\[
A=D_\chi F,
\]

将物理位姿切向量记为

\[
B=D_XF.
\]

消去位姿 nuisance parameter 后，地图信息由 Schur 补

\[
K_{\mathrm{eff}}
=A^TA-A^TB(B^TB+J_X)^\dagger B^TA
\]

控制。无位姿先验时，信息损失可以用两类数据切向子空间
\(\operatorname{Ran}(A)\) 与 \(\operatorname{Ran}(B)\) 的交叠和主夹角解释；有限位姿先验时，它变成 weighted shrinkage，不能继续机械地称为普通主夹角。

论文的价值不在于提出一个普适的“最佳轨迹公式”，而在于建立了以下可验证框架：

1. 用低秩 Schur defect 定位被位姿吸收的地图方向；
2. 用无先验广义保留谱与主夹角描述有限维局部混淆；
3. 用固定秩导数和人为秩事件区分平滑敏感性与非正则跳变；
4. 用轨迹、频率、噪声和 Born/full-wave controls 主动反证过强设计结论；
5. 揭示完整 Born pair 在 \(s\ne0\) 与 \(s=0\) 之间的 rank-stratum 奇异极限。

最终初稿位于：

- `experiments/idea_loops/loop_2026-09-03_16-25-02/experiment_pose_confounding_spectral_geometry/manuscript/main_round3.tex`
- 同目录 `main_round3.pdf`，58 页，Tectonic 编译成功并完成逐页渲染检查。

## 2. 自动科研流水线实际执行情况

### 2.1 使用的自动流程

Agentic-AI-Scientist 不是被 Python 脚本替代，而是作为外层自动研究循环实际运行。恢复命令为：

```bash
PYTHONUNBUFFERED=1 agentic-ai-scientist \
  --resume-loop-dir experiments/idea_loops/loop_2026-09-03_16-25-02 \
  --model deepseek-v4-pro \
  --final-max-turns 60 \
  --min-dev-rounds 3 \
  --max-safety-rounds 8 \
  --num-cite-rounds 20 \
  --writeup-retries 3 \
  --codex-timeout 3600 \
  --worker codex \
  --codex-profile deepseek-flash
```

它完成了 3 个 full development rounds。`loop_state.json` 记录
`full_dev_rounds = 3`、`pending_deepen = false`；最终自动评审在
`development_decision.json` 中给出：

```json
{"decision": "lock", "locked": true}
```

这表示自动评审认为当前多实验族、交叉核验和反证材料已经足以进入论文写作阶段；**它不表示数学证明、真实系统验证或发表验收已经完成。**

自动最终包：

- `experiments/idea_loops/loop_2026-09-03_16-25-02/final_2026-09-04_00-53-47_pose_confounding_spectral_geometry_deepening_round3/`
- 机器可读总结：上述目录中的 `experiment_report.json`
- 锁定证据：上述目录中的 `development_decision.json`

### 2.2 自动写作包装器的边界

自动研究及评审已经完成并锁定；之后通用 post-lock write-up 包装器因两个工具兼容问题未能独立完成最后一步：其客户端工厂不接受 `deepseek-v4-pro`，且通用 LaTeX 模板路径 `ai_scientist/blank_arr_latex` 在当前启动目录下不可解析。该失败发生在研究结果锁定之后，不影响实验产物，但意味着最终论文稿需要由 Codex 父级对已锁定材料进行科学审查、语义修订和编译。本文件与 58 页初稿完成了这一步。

### 2.3 各层职责

- **Agentic-AI-Scientist**：外层假设—实验—反证—评审迭代与 lock 决策。
- **ScholarQA research skill**：邻近先验工作、claim ledger、引用核验和 novelty boundary。
- **项目 Python 环境**：有限维模型、Jacobian、Schur 几何、反例和数值复现。
- **Codex 父级**：核心数学语义、Family 15/15b 区分、证据能否支撑论文 claim、最终初稿和科学审查。

项目内 `Theory/`、各轮 `context/`、research seed，以及个人 ChatGPT“电磁逆散射”文件夹中的相关讨论均只读吸收；现有理论没有被自动流程用一套无关方向替换。

## 3. 最终主张分级

### 3.1 已证明或由精确有限维代数支持

以下结论可在论文中作为有限维命题使用，但必须保留其条件：

1. **Schur 消元与序关系。** 在相应 range condition 和 \(J_X\succeq0\) 条件下，位姿消元产生正半定信息损失；无先验情形为
   \[
   K_{\mathrm{SLAM}}=A^T(I-BB^\dagger)A.
   \]
2. **精确混淆与秩损失。**
   \[
   \ker K_{\mathrm{SLAM}}=\{u:Au\in\operatorname{Ran}(B)\},
   \]
   且有限维秩损失等于两切向 range 的交维数。
3. **无先验保留谱。** 在 \(K_{\mathrm{IS}}\) 的可观测支撑上，广义特征值满足 \(\rho_i=\sin^2\theta_i\)。这不是逐个普通 \(K_{\mathrm{IS}}\) eigenmode 的 retention。
4. **低秩 defect。** \(\operatorname{rank}(K_{\mathrm{IS}}-K_{\mathrm{eff}})\leq\operatorname{rank}(B)\)。
5. **固定秩导数。** 固定秩邻域内可对 \(P_B=BB^\dagger\) 使用投影导数，对简单广义特征值使用 Rayleigh 型导数；重根必须压缩到相应 eigenspace，秩事件必须分开处理。
6. **完整 Born scaling proposition。** 若
   \[
   F_B(\chi,X)=A_0(X)\chi,\qquad \chi=s\bar\chi,
   \]
   则
   \[
   A_B(s)=A_0,\qquad B_B(s)=sB_1.
   \]
   因而每个 \(s\ne0\) 上的无先验 pose projector 和
   \(K_{\mathrm{SLAM}}\) 对 \(|s|\) 不变；在 \(s=0\) 处
   \(B_B=0\)、\(K_{\mathrm{SLAM}}=K_{\mathrm{IS}}\)，通常构成 rank-stratum 奇异端点。若 \(J_X=\alpha I\)、\(\alpha>0\)，则
   \[
   \|K_{\mathrm{IS}}-K_{\mathrm{eff}}(s)\|_2
   \leq \frac{s^2}{\alpha}\|A_0^TB_1\|_2^2.
   \]
   半正定先验不由这个上界自动覆盖。

### 3.2 已由程序数值支持，但不是证明

1. 2D scalar Helmholtz self-cell quadrature 与参考积分吻合到约 \(3.13\times10^{-16}\)。
2. 地图 Jacobian、位姿 Jacobian及完整 Born 位姿 Jacobian 的中心差分误差呈约二阶收敛；Family 15b 的拟合斜率为 1.999997。
3. 一般化 Family 5 的分解残差在约 \(10^{-15}\) 量级；简单谱分支预测误差呈约二阶。
4. 六重 \(\rho=1\) cluster 必须作为整体处理；压缩导数范数约 \(2.30\times10^{-10}\)，不能人为挑单个特征向量作稳定物理解释。
5. 人为 rank event 产生 \(\|\Delta P_B\|_2\approx1\) 和约 0.140 的相对信息跳变，验证固定秩公式不能穿越秩事件。
6. 完整 Born pair 的解析位姿导数通过 18 个位姿坐标的有限差分检查；非零 scale projector 残差约 \(4.52\times10^{-14}\)。

### 3.3 被实验反证或明确降级的命题

1. **不存在当前实验支持的普适轨迹排序。** 原排序假设在 24 个复核单元中为 0/24；raw 与 normalized 排序只在 12 个 scene-frequency 单元中的 2 个一致。
2. **ring 的 6/7 rank 现象不是已证实的真实秩下降。** 它来自对近零奇异值平方后再使用 machine-rank 阈值；分解残差仍约为 \(5\times10^{-16}\)。
3. **完整非线性执行误差没有获得 uniform certificate。** 只保留 affine tangent 层的有限维条件性证书；500 方向 nonlinear samples 只是经验检查。
4. **Family 15 不是完整 Born-SLAM 模型。** 它使用 mixed pair \((A_{\mathrm{Born}},B_{\mathrm{full}})\)，只能隔离地图切向近似误差；完整 pair 由 Family 15b 提供。
5. **外部 Born 论文没有被数值复刻。** 文献中的几何、单位、噪声、参数化和归一化不完整一致，因此这里只能做算子层对照，不能声称数值 benchmark 匹配。
6. **Family 10 的 free-pose covariance mismatch 尚未解决。** 相关结果不能上升为 estimator-MSE 等价结论。

### 3.4 仍然开放、应交给 GPT Pro 或后续研究的问题

- 非闭 range 紧算子下的无限维主夹角、稳定横截性和离散极限；
- 连续 Helmholtz 的 \(SE(2)/SE(3)\) gauge 定理、商流形和 gauge breaking；
- Born 远场 Fourier/Ewald 几何下的 missing-wedge 显式角界；
- 多频共享位姿补偿的 generic transversality，以及频率选择是否具备 weak submodularity；
- 半正定 pose-graph prior、运动图 gauge 和 prior-whitened canonical correlation；
- 重根、固定秩流形与真实秩事件的统一扰动理论；
- 3D vector Maxwell、天线极化和 reciprocal MIMO 约束；
- global identifiability、非线性优化收敛、cycle skipping、真实数据及硬件验证；
- 弱信号下“projector amplitude-invariant”与实际统计可检测性之间的桥接。

这些问题已被整理为可直接交给网页版 GPT Pro 的中文任务包：`research/GPT_PRO_HARD_PROBLEMS_ZH.md`。

## 4. 五条核心理论验证轴及扩展 controls

用户最初整理的五类验证目标已经被覆盖，并在自动回路中扩展为多个 family：

| 验证轴 | 对应程序/证据 | 最终结论 |
|---|---|---|
| 前向离散、自作用项与 Jacobian | `validate_self_cell.py`, Family 1, Family 15b | self-cell 高精度；map/pose/Born pose Jacobian 通过 FD 收敛检查 |
| Schur、主夹角、gauge 与低秩 defect | Families 2, 3, 5/5b/5b2 | 有限维代数成立；重根和 rank event 必须分开处理 |
| 频率、SNR 与轨迹设计 | Families 4/4b/4c, 6, 11, 12 | 绝对信息与归一化 retention 不同；没有普适轨迹排序 |
| 反例、rank autopsy 与鲁棒性边界 | Families 9, 13, 14 | machine-rank 假象已定位；完整 nonlinear robustness claim 降级 |
| Born/full-wave 过渡 | Families 15, 15b | mixed control 与 complete pair 已分离；发现 \(s=0\) rank-stratum 语义陷阱 |

另有 Monte Carlo、online toy 和 SOM 对照用于检查统计及空间语义，但它们不改变上述证据边界。

## 5. 本轮父级实际复跑

所有命令均从
`experiments/idea_loops/loop_2026-09-03_16-25-02/experiment_pose_confounding_spectral_geometry/`
执行，使用项目 `.venv/bin/python` 和 Apple Silicon CPU。

| 命令 | 结果摘要 |
|---|---|
| `.venv/bin/python src/family5_parent_generalized.py` | exit 0；overall gate true；分解残差约 \(10^{-15}\)；rank-event projector jump 约 1 |
| `.venv/bin/python src/family12_trajectory_replication.py` | exit 0；旧排序 0/24；raw/normalized 仅 2/12 一致 |
| `.venv/bin/python src/family13_rank_autopsy.py` | exit 0；ring 现象判为 tolerance classification；人为 rank event 复现 |
| `.venv/bin/python src/family14_robustness_scope.py` | exit 0；affine certificate 保留；full nonlinear claim 明确降级 |
| `.venv/bin/python src/family15_born_control.py` | exit 0；mixed Born/full-wave isolation；composition identity residual 0 |
| `.venv/bin/python src/family15b_born_pose_control.py` | exit 0；complete Born pair；所有 gates true；FD slope 1.999997 |
| `.venv/bin/python -m unittest discover -s src -p 'test_*.py' -v` | 5/5 tests passed |
| `.venv/bin/python src/validate_self_cell.py` | exit 0；最大绝对误差 \(3.13\times10^{-16}\) |
| `.venv/bin/python src/family1_pilot.py` | exit 0；map/pose FD 斜率约 1.976/2.000 |
| `.venv/bin/python src/family2_parent_correction.py` | exit 0；Schur/factorization 残差约 \(10^{-16}\) |
| `.venv/bin/python src/family3_parent_correction.py` | exit 0；Born 分支约二阶；径向旋转 pose-only stabilizer 被保留 |
| `.venv/bin/python src/family4_parent_controls.py` | exit 0；排序依赖 metric、standoff 与 normalization |

JSON 记录、图和解释分别存放于实验目录的 `results/`、`figures/` 和 `notes/`。任何聊天中四舍五入的数字均应让位于这些 JSON 记录。

## 6. 论文语义修订记录

最终父级审查对自动初稿做了以下关键修正：

1. 标题明确为 **Executable Finite-Dimensional Study**，避免暗示连续理论或真实系统已经解决。
2. 摘要前置失败结果：0/24 轨迹排序、2/12 指标一致、rank tolerance autopsy、nonlinear robustness demotion。
3. 将 Family 15 明确标为 mixed-tangent isolation；新增 Family 15b 的完整 Born pose tangent、命题、证明、数值表和图。
4. 无先验 \(\sin^2\theta\) 与有限先验 weighted retention 分开；有限先验的 arccos 可视化不再称作 principal angle。
5. 阈值化 log-volume 改称 **positive-support log-pseudovolume**，避免把伪行列式指标当成完整体积信息。
6. 将“certify derivative”降为“numerically support implemented derivative”。
7. 相关工作加入 radio SLAM PCRB、radar autofocus、joint SAR phase/image、bilinear identifiability、TDOA observability、simultaneous inverse scattering/transmitter localization、SOM 和 multifrequency 邻居。
8. 全文不使用“首次”“没有先验工作”等全局优先权措辞。
9. 添加逐 family claim-status ledger、复现命令、随机种子、SHA-256 来源记录和 AI 使用披露。

## 7. 文献与 novelty 边界

当前可辩护的 novelty 是**组合层和执行层**的，而不是声称每个数学部件都是新的：

- nonlinear full-wave volumetric **map tangent** 与 physical **pose tangent** 在 whitened/realified data space 中的显式并置；
- pose-eliminated low-rank Schur defect 与无先验 principal-angle retention 的统一解释；
- Born \(s\ne0\) projector invariance 与 \(s=0\) rank-stratum 端点的显式命题；
- 一组会主动推翻普适轨迹排序、错误 rank 判定和过强 robustness 叙述的 executable falsifiers。

不得声称“一般性的首次联合逆散射与定位”。Karthik–Ghosh 等邻近工作已经涉及同时 reconstruction 与 transmitter localization。检索还受到 Semantic Scholar 匿名接口 429 限流，因此 novelty 结论必须写成 **retrieval-bounded**，而不是全局穷尽证明。

证据材料：

- `research/literature/TARGETED_PRIOR_ART_ADDENDUM.md`
- `research/delegated/scholarqa_prior_art/claim_ledger.md`
- `research/delegated/scholarqa_prior_art/verified_references.json`

## 8. 验收表

| 层级 | 状态 | 说明 |
|---|---|---|
| Agentic-AI-Scientist 研究循环 | PASS | 3 rounds；自动 evaluator `lock` |
| 有限维程序复现 | PASS | 核心 Python、5/5 unit tests、controls 均成功执行 |
| 理论语义审查 | PASS（有限维范围） | Schur、主夹角、rank event、Born scaling 条件均已区分 |
| 初稿 LaTeX 编译 | PASS | Tectonic 成功；58 页；无 overfull box、fatal 或 unresolved reference |
| PDF 视觉质检 | PASS | 58 页全部渲染；全页 contact sheets 与关键页原尺寸检查无可见裁切/重叠 |
| 连续/无限维定理 | OPEN | 已转入 GPT Pro 困难问题包 |
| 3D vector Maxwell | NOT RUN | 不得从 2D scalar 自动外推 |
| 真实数据/硬件/在线系统 | NOT RUN | 当前结果不是 production validation |
| 外部论文数值 benchmark | NOT ESTABLISHED | 缺少完全一致的公开设置 |
| 发表级最终稿 | NOT YET | 当前为结构完整、证据可追溯的初稿；仍需作者署名、目标期刊格式和理论扩展 |

## 9. 最终建议

当前最有价值的下一步不是再堆一个随意参数扫描，而是二选一：

1. 优先解决 GPT Pro 任务包中的 P1/P2/P5，把有限维谱几何提升到有明确闭包、gauge 和半正定先验条件的定理；或
2. 建立一个可复刻的外部 inverse-scattering/localization benchmark，再评价这一谱诊断能否预测真实重建退化。

在任一步完成前，论文应坚持当前副标题和 claim ledger，不把局部有限维信息几何写成全局可辨识性或真实 SLAM 性能保证。
