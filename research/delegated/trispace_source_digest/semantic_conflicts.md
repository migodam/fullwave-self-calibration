# 语义冲突与消歧记录（semantic conflicts）

本文档重点整理 TriSpace 相关源之间的语义冲突，尤其是：current-space SOM vs map/pose data-space 几何、TriSpace 定义、unknown G、equivalent current、self-calibration。

## 1. current space vs data space：对象的“住所”

### 1.1 固定语义（多源一致的部分）

以下“住所”划分在 A、B、C 与 Theory §9.6 中是一致的：

| 对象 | 所在空间 | 说明 |
|---|---|---|
| J / j（induced current/contrast source） | current space H_J | SOM 的分解对象 |
| G_S（current-to-data） | H_J→Y | Chen SOM sensing operator |
| G_D（internal/domain propagation） | H_J→H_J（current/domain） | state/domain operator |
| V_S^±、V_D^± | current space H_J | SOM/TSOM 的 dominant/complementary 谱子空间 |
| E_inc(X) | domain field（与 current 相关的 domain side） | 受 pose 影响 |
| χ、δχ（map） | map space H_χ | contrast field |
| T_χ=∂J/∂χ | H_χ→H_J | map→current pullback/注入 |
| A=G_S T_χ（map tangent） | H_χ→Y（data） | Ran A⊆Ran G_S⊆Y |
| B=D_XF（pose tangent） | H_X→Y（data） | pose→data |
| K_IS、K_SLAM、K_eff | map/data 二次型（operator on H_χ） | Fisher/Schur retention |
| \operatorname{Ran} A、\operatorname{Ran} B、principal angles | data space Y | P9/上下文都强调这是 data-space 几何 |

### 1.2 已确认的非法比较/混淆（多源共同禁止）

1. 直接把 \operatorname{Range}(D_XF)（data space）与 V_S/V_D（current space）取 intersection/“三空间交集”。B 行 73–110、336–349 明确禁止；C 行 151–162 要求“与 SOM current space 兼容”的构造；A §4.9 与 §9 audit 同步禁止。
2. 把 Chen V_S^+ 称为“map observable eigenmodes”。A 定理 4.66–4.69 证明只有 A=G_S T_χ 的 range containment；V_G_S right singular vectors 与 V_A（map space）没有内禀 Grassmann distance。正确可比较对象：data-left subspaces、T_χ V_A push-forward 再正交化、或 pullback 谱 T_χ^*G_S^*G_ST_χ vs G_S^*G_S。
3. 把 K_eff 的 eigenvectors 当作 SOM modes（A §1.1#12、§9 audit、C 行 231–236）。K_eff 的 eigenvectors 是 map/reduced-coordinate modes；SOM modes 是 current-space objects。
4. 把 G_D/V_D 与“位姿 data space”混写。V_D 属于 current/domain operator；pose tangent 是 data-space Ran B（A §9.1 语义一致性检查第 2 条、Theory §21）。
5. 用 ordinary eigenvalue ratio λ_i(K_eff)/λ_i(K_IS) 代替 generalized retention（A §1.2、§9；Theory §8.3）。应使用 polar retention operator R_geom 或 generalized eigenproblem。

## 2. TriSpace 的定义冲突（同名词不同对象）

| 出处 | TriSpace 所指 | 是否含 V_P | 形式化程度 |
|---|---|---|---|
| Source B | current-space 三空间几何 V_S–V_D–V_P；X 未知下 V_S^±(X)、V_D^±(X) pose-dependent | 是（须构造） | 定位/目标，未形式化 |
| Source C | Sensing×State×Pose 三机制分类；J_safe+J_pose-confounded+J_state-recoverable+J_unresolved；不是形式 intersection | 隐含在 Pose/calibration ambiguity | 定位/分类方案，未形式化 |
| Source A §4.9/§6.3 | S-fold（G_S）→D-fold（state/domain）→pose-fold 的合法 two-stage reduction；不声称三个交换投影 | 否（用 data-space W/P_B 做 pose 层） | 有定理与算法，但对象是 map/reduced-coordinate，不是真正 current-space V_P |
| Theory §9 | TriSpace/Tri-fold = S+D+pose 分层框架；“八格精确分解”不成立 | 概念性 | 讨论内部骨架，V_P 未构造 |

冲突点：

1. **“TriSpace”是否必须是三个 current-space 子空间**：B 是；C 更强调“三种信息机制”；A 的三层实际是 two current folds + one data-space pose fold，只有两折真正在 current space。
2. **2^3 形式交集 vs 最小解释分类**：B/C 都反对为制造 2^3 intersection；C 反对预设四类正交直和；A 定理 4.74 证明 projector product 不是 intersection projector；Theory 9.3 记录八格仅在某些交换条件成立。
3. **V_P 是“pose 诱导的等价 current perturbation 空间”还是“current mode 的 pose-confusion 分类”**：B 指前者（子空间几何），C 指后者（每 mode 回答三个问题）。若要写成一篇论文，必须先消歧或声明只研究其中一个。

## 3. Unknown G 的三层含义（不能混为一谈）

1. **SOM 层的 unknown G（B/C 主线）**：X unknown ⇒ G_S(X)、G_D(X)、E_inc(X) unknown ⇒ V_S^±(X)、V_D^±(X) 和 P_S^±(X)、P_D^±(X) pose-dependent。核心 quantity：\|P_S^+(X+δX)−P_S^+(X)\|、leakage \|P_S^-(X*)P_S^+(X̂)\|、principal angles/spectral gap、rank/threshold events。这是 current-space 谱坐标“本身变成待估量”。
2. **A/B/K_eff 层的 fixed-nominal tangent**：A=D_χF、B=D_XF 在名义点定义；回答“整个 joint problem 是否 locally identifiable”，不回答 SOM 谱坐标。B 行 245–311、C 行 215–239 要求只作辅助层。
3. **A §4.7/P7 的几何导数层**：参数 ball 上 G_S(X)、G_D(X)、e_inc(X) 对 pose 的导数/高阶导数决定 A、B、K_eff 的鲁棒常数；该层只在 world-fixed map、D_χ 固定、外部传感器 standoff 等条件闭合（4.54、4.56、4.60）。moving grid/self-cell 属 **open**（8.8）。

三者回答的问题不同，不应在论文中把 K_eff 的连续性结果当作“SOM 子空间连续性”的证据（A §1.2#9、Theory §11、21）。

## 4. Equivalent current / V_P 的定义状态

问题：数据空间中一个 pose 扰动 BδX 对应 current space 中的哪些扰动，可与 current/data 分解放在同一 ambient space 中分析？

源文件给出的候选构造（都不是已证明定义）：

- **pullback 经 T_χ / map→current 映射**：只对 map 可激发方向有效，Ran A⊆Ran G_S 的严格包含说明并非所有 current modes 都能由 map tangent 激发（A 反例 4.67）。
- **minimum-norm equivalent current**：对 BδX，求极小 \|G_S δJ_eq−BδX\| 的 δJ_eq，或用 G_S 的 pseudoinverse；涉及闭值域/closure 语义（类比 A 推论 4.2 的精确/渐近补偿问题）。
- **joint block operator / quotient construction**：把 (G_S, T_χ, D_XF) 组合成可比较的 block operator；目前无形式化。
- **data-space W 层（A 实际做法）**：A 用 P_B 或 prior contraction W 直接做 pose elimination，不构造 current-space V_P。这合法但**不满足 B/C 的 TriSpace 目标**。

已禁止的捷径：V_P := Range(D_XF)（data space）；P_S P_P 或三 projector 乘积当 exact intersection；2^3 象限直接当精确分解；把 K_eff 谱/defect 直接搬到 current space。

**结论：V_P 的定义问题是 open，且是 B/C 反复指出的“TriSpace 从直觉变成定理”的关键缺口。**

## 5. Self-calibration 的定位冲突

- C 行 33–37：self-calibration 不是 novelty（array/radar/bilinear 长期研究；未知 geometry inverse problem 有先例）；可主张的交叉点是 “unknown sensing geometry + nonlinear volumetric full-wave inverse scattering + SOM/TSOM spectral reduction”。
- C 行 40–46：与 phaseless imaging 互补但不同（保留 coherent phase，估计造成 phase/model mismatch 的 calibration variables）。
- A §7.2 C14/C15：joint inverse scattering + transmitter localization 有 Karthik–Ghosh 2023 等直接前例；RF-SLAM map-vs-pose FIM/CRB 有 Deutschmann 等前例；因此“SLAM 名下的 Fisher 观测性”也不 new。
- 命名冲突：B 坚持“TriSpace SOM-SLAM / Inverse-Scattering SLAM”主线（B 行 312–349），C 反对“将 SOM 方法应用到 SLAM”的表述并建议 self-calibrating inverse scattering（C 行 240–264）。两者实际目标接近，但对外命名/定位不同，且 C 更强调“SOM 降维机制在未知几何下是否存活”。

## 6. A/B/K_eff 层与 TriSpace 层的层级关系（多数源认可的折衷）

折衷表述：

full-wave physics → TriSpace SOM current-space theory → joint SOM pose–map reconstruction，旁边用 A/B/K_eff 作 local identifiability / Fisher-level consistency check（B 行 245–311 的图；C 行 215–239；Theory §9.3 的 hierarchy）。

仍需主 Codex 线程裁决的冲突：A §9.2 的“最终可主张中心句”仍然只落在 data-space retention geometry，而 B/C 要求论文主线前移到 current-space TriSpace/V_P。本 digest 只记录该分歧，不替 Codex 判断哪个是正确主线。

## 7. 最小消歧建议（供后续使用）

1. 写论文/讨论前先固定术语表：TriSpace 当前指三种机制分类还是三子空间几何；V_P 是否必须存在且位于 current space。
2. 永远给对象标注所在空间（map/current/data/pose），禁止无映射跨空间比较。
3. Unknown G 的讨论分开三层：SOM 谱坐标层、nominal tangent/Fisher 层、参数导数鲁棒层。
4. Equivalent current 的任一候选构造都必须以 G_S/G_D/full-wave forward model 出发，并回答闭值域/closure 与精确/渐近补偿语义。
5. 不要用 A/K_eff 的谱替代 V_S/V_D/V_P 的谱；也不要反方向声称 V_S/V_D 直接给出 map observability。
