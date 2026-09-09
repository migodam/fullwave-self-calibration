# A2 后的新一轮 GPT Pro 理论与算法设计总 prompt

你是严格的电磁逆散射理论与数值算法合作者。请继续 TriSpace SOM / self-calibrating full-wave inverse scattering，不要改成一般信息几何论文，不要把 SOM 理解为 self-organizing map。目标是保留 coherent phase，通过低维几何 self-calibration 改善定量成像；必须说明为什么需要 SOM 的结构，而不只是普通 joint inversion。

请先阅读随本 prompt 附上的当前论文、A2 theorem package、THEORY_AUDIT、E4_RESULTS、a2_highdim/summary 及 critic。原始 A2 三个文件是理论输入，不代表全部算法已实现或实验已通过。

## 不能重新包装或否认的事实

1. exact physical state elimination 后没有独立自由电流，$r_{\rm free}=0$。$L_{\rm det}$、$r_{\rm num}$、$r_{\rm free}$ 不同。自由 nuisance 的信息损失定理不适用于提高同一状态方程的数值精度。
2. phaseless SOM 仍有 current 结构。已有 multifrequency SOM 联合复数 transmitter gain/phase calibration（Idriss/Raj, arXiv2503.07316）。不得声称首次 SOM+calibration，也不得说等效本身增加 Fisher 信息。
3. $\rho$ 是 supported map relative retention；绝对 pose singular value 和 bias 要另算。高曲率、小残差不证明在真 branch。Gauge 必须先固定/商掉；先白化/实化；complex-safe kernel 必须对乘 i 不变。
4. 2,400 次冻结 E4 运行在 200/800 工作预算下全部联合失败，约化路径没有接受的迭代。原因包含保守的窄矩阵乘法计价及阶段预算导致未充分执行。不能据此声称方法内在无效或等效，更不能修改旧测试再宣布原 endpoint 通过。
5. 49 参数探索性成像中，direct coherent joint 的空间 RMSE0.0261，固定错误 pose0.1054，matched intensity joint0.0430；这是一般自校准的价值证据，不是 SOM-specific 优势。
6. 新 passive bound 有条件地可计算。18 组合中12正损耗被界包住、6无损拒绝。残差补充让4个低频 case 在 rank72通过，最高频4个到rank128仍拒绝。不能把数值误差证书写成物理模型真确或几何真值证书。
7. 状态界/真误差约 13–210，切线界/真误差约 41–9,051；需要真正降低 conservatism。物理有效 Fisher 的 10 个原场景补查通过，别再把它写成没验证，也不要增加样本计数。
8. 原实测 worker 的波向、吸收符号和 test gain 重拟合有错，已经撤回。修正后的圆柱模型留出归一化平方残差约 0.02162，epsilon 实部约 3.427；这是 object-center/material 拟合，不是阵列偏差恢复。数据没有天线偏差真值；不得沿用原报告“epsilon 不可辨识/高频失效”的错误解释。

## 已有具体物理上界

二维 outgoing $G=iH_0^{(1)}(kr)/4$，等面积网格 h，disk self-cell，$\chi=(1+i\tau)u$，$u_i>0,\tau>0$。令 $T=\operatorname{diag}\sqrt u$，$a=h/\sqrt\pi$，$\delta=k^2h^2/4-\pi kaJ_1(ka)/2$，$\zeta=\tau/(1+\tau^2)-\max(\delta,0)u_{\max}$。当 $\zeta>0$，$\|T^{-1}M^{-1}T\|\le g=[\sqrt{1+\tau^2}\zeta]^{-1}$，输出误差不超过 $\|WST\|_F g\|T^{-1}(M\widetilde j-b)\|$。切线残差可以同样控制。该上界保守，且不是 interval arithmetic。

## 任务 A：真正便宜、且不空的目标导向证书（最高优先级）

能否用输出/伴随残差对 $WSM^{-1}z$ 或目标参数方向给出严格而显著更紧的上界，避免完整 current norm 与 Frobenius 放大？请给可计算条件、所有额外求解/算子成本、失败输出与非空的物理参数区间。不能预先调用精确解、真误差、真实 pose 或未经验证的逆范数估计。如果无法普遍做到，给出不可实现性/成本下界和可行受限版本。

## 任务 B：设计能实际更新 pose 的 SOM-informed 算法

现有 $\min_c\|MUc-b\|$ 的 sensing/domain 基不一定包含 material-weighted currents。请设计 sensing、projected domain、state residual、tangent/adjoint residual 的分工与增广规则。哪些方向只需解状态，哪些对 pose/map 的输出误差必要？如何用 complex-linear RRQR/SVD 处理重复方向？Frozen chart 下近似物理 Jacobian 与移动近似模型的导数必须区分。

给出可实施伪代码、逐次变量、算子形状、停止/拒绝/回退规则。证明的 descent 要与具体步矩阵、约束、非零误差 floor 一致。明确其相对于普通 residual-enriched ROM / block Krylov / direct 的独特收益假设；若不需要 SOM，直接说这是一般 PDE 算法，不要强行贴标签。

## 任务 C：相位 branch 与频率推进

不能用 inverse-GN covariance 假装真实误差 coverage。能否构造可计算且不会轻易接受错误 branch 的 set-based、profile-likelihood、独立验证或多候选频率推进规则？允许输出“仍有多个候选，暂不解锁高频”。需要说明有限孔径、幅值近零、clock 与 geometry 耦合时是否失效。单路径 cosine 的凸区不能直接改名 full-wave basin。

## 任务 D：几何、clock、complex gain、coupling 的边界

对共同 clock delay、每 TX 复数增益、阵列刚体参数，给统一白化 tangent 与 joint gauge。哪些频率/视角能区分，哪些严格混淆？给少量最有诊断力的反例和一组可计算 acquisition 条件。不要随意增加 nuisance 到满秩后宣称几何本质不可辨识。区分独立 channel phase 与低维时钟/几何模型。

## 任务 E：重新设计可区分贡献的公平实验

提出 dimension-aware 的 RHS/operator/basis/QR/LU 计价，支持多个 RHS 的缓存与复用，给与 direct 相同的精度、预处理、continuation 与 tuning。先在新 tuning 集验证每个方法真的执行了有意义的迭代，再冻结新测试。保留旧 primary 的负结果；不得重新拿同一 test seeds 调参数。必须比较 SOM-specific 基、普通 residual enrichment、block Krylov、直接伴随求解，以及匹配时才加入 source extension。明确最小有判别力的 ablation 和预先登记的失败判据。

## 任务 F：实测与投稿叙事

Fresnel calibrated laboratory data 加人工坐标误差，只能证明受控元数据扰动的 robustness。怎样得到足够可信的 phase-center/incident calibration/坐标锚定，避免 calibration nuisance 把所有几何误差吸收？给出不依赖真像的识别方案和必要 metadata。若现成数据不足，明确最低额外测量设计，不要编造未知几何真值。强制检查时间约定、outgoing Hankel、被动介质虚部和 Poynting 通量的一致性；增益必须在训练数据估计后固定到留出预测，不能把测试集重新 profile 当验证。

## 输出格式

先用中文给总体 verdict，再逐任务输出：已知/新推导/条件性/开放；命题与所有假设；证明或反例；可执行算法；计算成本；最小实验；与最近前人的精确区别。最后给 TAP 风格的最小论文主张集及正文/附录布局。不要把已解决的 A2 投影恒等式再次扩成几十页；优先解决非空、非 oracle、可实施、可比较这四个缺口。
