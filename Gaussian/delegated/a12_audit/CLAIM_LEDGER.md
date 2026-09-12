# A1 / A1_2 独立理论来源核查与可证伪清单

日期：2026-09-11  
角色：独立理论审计；不替代根任务的最终理论、方向或新颖性裁决。  
审计对象：`Gaussian/Theory/ A1.md` 与 `Gaussian/Theory/A1_2.md`，重点为 A1_2 新增的 Gaussian/Hermite 局部散射网络与 NN/GNN 路线。

## 审计边界与证据等级

- **代数成立**：在明确的有限维/算子可逆性和空间定义下可直接推出。
- **条件成立**：缺少条件会使公式错误或失去所声称的物理意义。
- **构想**：可以实现，但尚无误差界、数值闭合或材料可识别性证据。
- **来源 A/B/C/U**：A 为本轮检查全文；B 为出版社摘要；C 为元数据；U 为未找到可核文本。
- “未找到等价工作”只登记为未闭合检索，绝不作为新颖性证据。

本轮没有运行数值模型或训练网络。实际执行的检查为：逐行读取 A1_2 全文 4612 行；核对 A1 的既有 claim ledger 和与当前主线有关的原文段落；检查本地 Chen–Zhong–Agarwal 综述、Kanaun、Naik、Baussard 与 Gaussian-mixture split/singularity 证据目录；针对 A1_2 点名的四类 T-matrix/GSM/NN 文献检索出版社或机构页面；检查文件哈希。审计时 SHA-256 为：A1_2 `fb560b5f7951c7e0d6695f5398d1d76d12a811f503b70d5cd033f918495cf0b3`，A1 `027e0811f8f5139d429514ba4b54650da77348053314241ab99fd0d6e9eb0f18`。

## 先行结论

当前网络表达可以作为**受控降阶模型和 learned solver 架构**继续实验，但 A1_2 把四个层次接得过快：

1. 完整函数空间上的 component (t_i)；
2. 有限维 incoming/outgoing modal ports；
3. 电路意义的 power-normalized (S_i)；
4. 外部有限 Tx/Rx 所测得的 (M=C K B)。

这四者不是同一个矩阵。只有第 1 层、完整模态且相关逆存在时，component multiple-scattering 分解可称为 exact full wave；进入有限 ports 后必须加入截断闭合误差或被消去模态的 self-energy。外部 (M) 一般不能识别内部 (K,T_i,U)。因此 `full-wave exact`、`low-port`、`direct de-embedding`、`material recovery` 是四个独立 gate，不能由同一条 resolvent 恒等式一起关闭。

## 行号化公式与假设 ledger

| ID | 原文行号 | 主张或公式 | 独立核查 | 必须修正 / 最小可证伪检查 |
|---|---:|---|---|---|
| C01 | A1_2 122–176 | scalar constitutive consistency 等价于 `rank[u,j]=1` | **错误** | 正确的全局条件是 `rank[u,j]=rank[u]`。当 (u\ne0) 时等价于 rank 1 且 (j\in\mathrm{span}(u))；当 (u=0) 时必须另行要求 (j=0)。单写 `rank≤1` 仍会让 ([0,j\ne0]) 通过。对实数、被动或色散材料，还须把比例系数限制到合法材料集合。 |
| C02 | A1_2 184–217 | (chi^*=u^*j/(u^*u)) 与 angle residual | **条件成立** | 仅在 (u^*u>0) 时定义。零场、极弱场必须 mask/加权，不能用任意正则分母把 (j\ne0) 解释成材料一致。实物参数为实数时应做 real-constrained projection；复材料应包含被动性和频率色散约束。 |
| C03 | A1_2 279–333 | Maxwell 中 (J=-i\omega\epsilon_0\chi E)，`J×E=0` | **必要条件不够** | 复三维向量的叉积为零可表示单照明复线性相关，但仍缺各照明共享同一 (chi)、(E=0\Rightarrow J=0)、各向异性张量材料、不同频率的 (-i\omega\chi(\omega))。应用 stacked rank 条件和合法材料投影。 |
| C04 | A1_2 338–410 | 复相量球面切向电流必有零点（hairy-ball） | **反例推翻** | 对 (r\in S^2)，令 (P_r=I-rr^T)，(J(r)=P_re_x+iP_re_y)。则 (J\cdot r=0)，且 (|J|^2=1+r_z^2\ge1)，处处非零。毛球定理约束实切向场；复 phasor 的实、虚部零点不必重合。Hodge 分解仍可复化使用，但拓扑必零主张必须删除或改成实方向场。 |
| C05 | A1_2 444–554 | Hodge–Gaussian 可表示任意切向电流 | **条件成立** | Hodge 分解要求闭定向曲面、声明 Sobolev 空间与边界条件；genus-(g) harmonic 部分的实维数为 (2g)，复化后为复维 (2g)。有限 Gaussian 势只能近似，不是精确张成。 |
| C06 | A1_2 633–738 | 有限非退化 Gaussian 和不能精确给出非零紧支撑硬物体 | **成立但适用域窄** | 这是解析延拓结论，不是有限分辨率近似失败结论。截断/compactly-supported RBF、soft partition、有限 ROI 会改变前提。 |
| C07 | A1_2 748–887 | Gaussian score softmax 在低温生成 hard cells 并可 gauge-fix | **条件成立** | (\tau\to0) 仅在 away-from-ties 点态成立；界面不可微且导数可随 (1/\tau) 变坏。必须含背景类和有界 ROI。softmax 对所有 logits 加共同**函数**都不变；只约束 bias 和为零未必消除参数族中全部 gauge。 |
| C08 | A1_2 891–1016 | illumination-independent reflectivity 只在 Born/弱散射合理 | **方向正确** | 展开需要相应算子有界并处于 Neumann/局部解析域。零场处不能除。用同一几何在 contrast、频率和视角 sweep 下直接测 (\rho_l=j_l/e_l) 的 illumination dispersion，作为 Born 失效量。 |
| C09 | A1_2 1145–1257 | Gaussian Gram 闭式与 (Q_G=\Phi M_G^{-1/2}) 使谱只依赖 `Ran Φ` | **条件成立** | 闭式是全空间、给定 normalization 下的 scalar (L^2) Gram；有限 ROI、surface、截断 patch 和求积须重算。(M_G) 必须正定；重叠导致秩亏时用 rank-revealing EVD/SVD 与伪逆。换基后 (Q) 只确定到 unitary，奇异值和物理子空间不变，单个坐标向量不唯一；截断阈值会改变结论。 |
| C10 | A1_2 1201–1257, 1862–1923, 2001–2037 | 一次 Gram whitening 足以给出物理 Twofold spectrum | **不充分** | 至少分开四套度量：material tangent Gram、current/state Gram、port flux/power metric、measurement noise covariance。实物参数配复数据时还应对 noise-whitened Jacobian 做 paired-real realification。`C_D v=γ(C_S+εI)v` 是正则 generalized Rayleigh pencil；大比值不表示可检测，必须同时设内部绝对幅度下限和接收端噪声可见下限。 |
| C11 | A1_2 1263–1316 | 高 overlap 的 Gram 条件数可直接作为 merge criterion | **仅诊断表示冗余** | 高条件数说明当前离散表示近线性相关，不证明物理上可合并；merge 会改变非线性材料与局部场。需在 whitening 后同时检查 held-out data change、state residual、材料约束和 refit 后误差。 |
| C12 | A1_2 1320–1454 | Gaussian first fold 与 restricted second fold | **可作为 Galerkin/Rayleigh–Ritz 版本** | 必须把表示误差与 data nullspace 分开。该第二 fold 是先限制在 (Q_N) 后对 (G_DQ_N) 做 SVD，并不等同经典 TSOM 原式中先找 (G_D) 主导方向再投掉 first-fold 部分；需用同一 current metric 做直接 numerical equivalence/ablation。 |
| C13 | A1_2 1476–1604 | (H_0,H_1,H_2,H_4) 是 Gaussian 参数层级 | **局部参数切向解释成立** | rotation 在 isotropic covariance 或重复特征值处不可辨；对称 split 的 (H_4) 是受约束的 rank-one tensor (d^{\otimes4}) 和单侧 birth 参数，不是任意四阶 Hermite coefficient。有限 split 还含 (H_6,H_8,\ldots)。 |
| C14 | A1_2 1608–1784 | (|\widehat H_n|\propto q^ne^{-\sigma^2q^2/2})，峰值 (q\sigma=\sqrt n)，据此判 split 有无意义 | **只是一维/径向理想 Born window** | 方向导数实际含 ((q\cdot d)^n)；(q\perp d) 时完全为零。Ewald shell、vector polarization、sampling、noise、绝对 amplitude 与 full-wave dressing都会改变可见性。0.85% 不是不可恢复证书。用 matched normalization 的 Fisher/held-out response，不得仅凭 (Q\sigma) gate。 |
| C15 | A1_2 1927–2115 | dressed tangent (R_\theta M_u\Phi_\theta) 自动物理 | **局部链式法则成立** | 需要 (I-XD) 可逆且远离不受控共振；它只对当前材料模型、已知 (D,e) 的局部切向物理，不覆盖 geometry/operator/calibration error。它是 parameter-to-state Jacobian，不由命名自动获得新颖性。 |
| C16 | A1_2 2279–2328 | restricted LS 的 exact bias–variance decomposition | **固定子空间时成立** | 要求线性模型、固定非数据选择的 (U)、正确 noise whitening 和声明实/复方差 convention。adaptive split/search 后 (U) 依赖数据，selection bias 和 multiplicity 使等式不再直接适用。 |
| C17 | A1_2 2332–2414 | 少量 Gaussian 捕获 stable subspace 的 principal-angle bound | **上界成立，优势是假设** | 若 (\dim U_G<\dim V_r)，通常会有未覆盖方向/最大角 (\pi/2)。必须报告维数、角度谱、range cutoff 和对 FFT/voxel 的相同成本比较；不能把待验证的 `P_G≪P_F` 当推论。 |
| C18 | A1_2 2657–2782 | (t_i=(I-V_iG_0)^{-1}V_i) 精确吸收 component 内部多散射 | **完整算子层面条件成立** | local resolvent 必须存在。随后写的 (V_i+V_iG_0V_i+\cdots) 还要求例如 (\rho(V_iG_0)<1)；逆存在本身不保证 Neumann 级数收敛。Gaussian 无限尾使 `local` 只是 component label，不是有界独立物体。 |
| C19 | A1_2 2727–2925 | (b=(I-TU)^{-1}Ta=(T^{-1}-U)^{-1}a) 严格 full wave，路径级数包含任意阶 | **第一等式优先；第二等式有额外条件** | ((I-TU)^{-1}T) 只需相关 resolvent；改写为 ((T^{-1}-U)^{-1}) 还要求每个 (T_i) 可逆，而 local/truncated (T_i) 可奇异。路径级数 (T+TUT+\cdots) 要求 (\rho(TU)<1)（或更强的 norm bound）。共振附近可有 resolvent 但该级数发散。 |
| C20 | A1_2 2929–3066 | Gaussian/Hermite 内部阶 (H_n) 可对应 VSWF/multipole ports；H0=monopole、H1=dipole、H2=quadrupole、H4=split | **关键混同** | (H_n) 是 Cartesian material/shape tangent order，multipole (l) 是 Maxwell incoming/outgoing angular-momentum channel。二者没有一一对应：degree (n) 一般含 (l=n,n-2,\ldots)，full-wave dressing和非球对称背景会混合更多 (l)；Maxwell 没有普通电辐射 monopole。必须显式计算 transfer matrix (P_{\rm VSWF}T M_uH_n)，按能量占比选 ports。 |
| C21 | A1_2 2628–3167 | 完整 local (T_i) 可直接变成有限 (r_i)-port，且仍为 exact full wave | **未成立** | 完整 (t_i) 是无限维算子；有限 port 截断后 (P A^{-1}P\ne(PAP)^{-1})。正确的保留空间 inverse 含 discarded-space Schur self-energy (A_{PQ}A_{QQ}^{-1}A_{QP})。必须随 (r_i) 增长展示 full VIE state/data closure，并报告 resolvent amplification。 |
| C22 | A1_2 2782, 2951–2987, 4199–4229 | overlapping Gaussian component 可按传统 single-particle T-matrix + translation network 解释 | **只能二选一地严谨化** | additive-potential/Faddeev 分解允许 overlap，但 local operator 仍在完整场空间，低端口性不自动成立；经典 isolated-object VSWF T-matrix/translation 依赖局部展开域，重叠 circumscribing regions 会破坏标准收敛条件。Gaussian 无限支持使所有组件严格重叠。需使用非重叠 patch/cluster、显式截断误差，或只称 operator component network。 |
| C23 | A1_2 3117–3255 | 若 (R\ll N)，full-wave solve 从 (N\times N) 降到 (R\times R)，(T_i) 可复用且 library 很小 | **计算构想** | 必须计入构造/插值 (T_i)、dense (U)、B/C、梯度、端口增阶和误差控制。平移复用只在均匀平移不变背景、端口随节点移动时成立；边界、substrate、antenna near-field 会改变 local environment。rotation similarity 要保留完整 (l)-shell 或 rotation-closed truncation。`small library` 在 anisotropy/contrast/resonance/frequency 上未证明。 |
| C24 | A1_2 3320–3344 | 2023 已有 HG basis + GSM two-port domain decomposition | **来源 B，身份核实** | Uysal & Akleman, AEU 171 (2023) 154884, DOI `10.1016/j.aeue.2023.154884` 的出版社摘要确实说明 cascaded two-port subdomains 与 HG coefficients。摘要足以阻止 generic `HG+GSM` claim；不足以证明其全文没有 adaptive inverse/material/model-order 机制，负面差异仍需全文。 |
| C25 | A1_2 3352–3467 | full internal access 时 (T^{-1}=K^{-1}+U)，offblock 可作 closure diagnostic | **仅理想完整坐标成立** | 要求已知同一固定 port basis 中的完整 (K,U)，且 (T,K) 可逆。有限端口的正确 effective inverse 可因 discarded modes 非 block-diagonal；噪声和 basis gauge 也会制造 offblock。该量不是外部有限测量可直接计算的 observable。 |
| C26 | A1_2 3471–3635 | 一般 (M=CKB) 不能反演内部网络；特殊 (M=(A^{-1})_{VV}) 给 dark self-energy | **边界判断正确，公式被后文过度推广** | Schur 公式还要求 B/C 选择同一正交 visible block、(A_{DD}) 可逆、该 principal block 可逆。真实 Tx 与 Rx 子空间可以不同，(M) 可矩形。此时没有单一 `VV` 主块，也不能一般声称 (M^{-1}=A_{VV}-\Sigma_D)。 |
| C27 | A1_2 3639–3714 | network Twofold 可用 first fold (CQ)、second fold (UQ_D) 或 (TUQ_D) | **构想，未等价于 classical TSOM** | incoming regular coefficients与 outgoing coefficients是不同 port spaces；(U) 与 (T) 的 domain/codomain和度量必须声明。若目标是从实验识别网络，还需同时考虑 source-excitable 子空间 `Ran(B)`；只对 (C) SVD 只回答接收可见性。second fold 必须验证比 (U) 的普通低秩/SVD、GSVD或正则 GN 更有效。 |
| C28 | A1_2 3718–3820 | parameter-space internal/receiver generalized eigenratio 找到内部强、接收弱方向 | **局部灵敏度工具** | 分母为零时方向严格外部不可辨，不能因 (gamma=\infty) 变成可恢复参数。需同时报告 (\|J_{recv}v\|) 的噪声门槛、(\|J_{int}v\|)、state residual 改善和最终材料误差；还要加入 (\delta C\,b)、(\delta Bx) 及 geometry derivative。 |
| C29 | A1_2 3824–3930 | 单节点材料更新是 (r_i\)-rank Woodbury，measurement 只需小逆 | **固定端口/固定 U,B,C 时代数成立** | 改 shape、位置、port basis/order 会同时改变 (U,B,C)，甚至 (T_i) 的坐标；不能用仅更新 (T_i^{-1}) 的成本代表所有 actions。near-resonance 时小 Schur/Woodbury 块病态，需回退全解和误差监测。 |
| C30 | A1_2 3990–4057 | moment-preserving split 是在原 (A) 上 append 一个 child 的 small-block augmentation | **错误实现语义** | 真 split 是 parent 被两个 covariance-adjusted children **替换**，会改变旧 parent block、两子之间耦合、所有相关 U 行列以及 B/C。原 (A) 保持不变只对应保留 parent 再 birth 一个 child，不是所推导的 split。可仍做低秩/block update，但秩与成本应按 remove-parent + add-two-children 推导。 |
| C31 | A1_2 4047–4057 | (O(d^4)) material/data continuity 自动与 network block augmentation 连续 | **不成立** | (d\to0) 时两个 child 重合，内部 factorization 进入 mixture singularity；translation blocks 可近场病态，端口 gauge 退化。外部 forward data 的四阶连续不意味着内部 (A,T,U) 或其逆条件数连续。必须分别测 data continuity 与 network conditioning。 |
| C32 | A1_2 4061–4099 | 移动节点只改一个 block row/column，perturbation rank (O(2r_i)) | **固定全局端口数时成立** | 同时更新 B/C；若局部 basis 随位置/方向重构或 port order 改变，需把 basis-transfer 纳入。端到端耗时仍需包括生成全部 (U_{ij}) 和求梯度。 |
| C33 | A1_2 4257–4301 | 被动 local scatterer 满足 (S_i^\dagger S_i\preceq I)，无损完整通道 unitary，互易时 (S_i^T=S_i) | **只对特定 S-port convention 条件成立** | A1_2 前文定义的是 transition (T_i)，不是 circuit/power scattering (S_i)。二者不可换名。contractivity/unitarity 只对完整 propagating、power-normalized ports；evanescent/reactive channels使用不同 metric。T-matrix reciprocity通常含 mode metric/parity，不必是普通 transpose symmetry。有限通道的 contractivity也取决于正交投影定义。 |
| C34 | A1_2 4285–4299 | NN enforce passivity + reciprocity + causality 即可保持物理 | **要求正确，方法未给** | 单频 `spectral norm≤1` 与 symmetry 不能保证跨频 causality。需用 passive/causal rational state-space、positive-real impedance再 Cayley变换，或等价可证参数化；跨频共享 poles/material dispersion。还须保留 radiation condition、analytic (U_{ij})、basis-gauge/SE(3) equivariance、complex phase和最终 full-physics residual。 |
| C35 | A1_2 4305–4355 | T-matrix inverse、T+SOM、T+physics-inspired NN 已有 | **来源身份基本核实** | Ishida 2010 DOI `10.1587/transcom.E93.B.2595` 出版社摘要明确“data→T-matrix→SVD radiating/nonradiating currents→dielectric cylinder”；Ye et al. 2013 DOI `10.1109/TAP.2013.2258878` 摘要核实混合 PEC/dielectric T-matrix inverse；Zhang et al. 2023 DOI `10.1049/rsn2.12419` 出版社全文核实 CSI+T-matrix+APU-Net。它们阻断 generic `T-matrix inverse/SVD/NN` claim。 |
| C36 | A1_2 4341–4355 | adaptive Gaussian local modal model、observable de-embedding、network visible/dark、cheap split scoring 可能新 | **未核验主张** | 当前本地证据和点名来源不足以做 priority 判断。至少还需查 multiple-scattering model reduction、network tomography/system identification、adaptive multipole truncation、reduced-basis scattering、graph wave solvers、learned iterative solvers，以及 Gaussian/RBF adaptive inverse media 的全文最近邻。 |
| C37 | A1_2 4413–4479 | 从 (K_V^{est})、dark self-energy 恢复各 (T_i\to(\epsilon_i,\Sigma_i)) | **存在结构不可识别性** | 内部 realization 在 (b'=Lb) 下有 (A'=LAL^{-1},B'=LB,C'=CL^{-1})，外部 (M=CA^{-1}B) 不变。固定解析 U、节点几何、port normalization和 block gauge 可减少自由度，但 partial Tx/Rx 下仍需唯一性/稳定性证明。学习一个 (\Sigma_D) 不会自动补回不可识别信息；最多得到 data-equivalent effective model。 |
| C38 | A1_2 4491–4612 | adaptive low-port Gaussian network 有双重 compression，SOM 选 ports | **核心可证伪假说，未证** | port rank 不只依赖 (k\sigma,\chi,) separation 与 tolerance，还依赖 near-field excitation、gap、loss、distance-to-resonance、background、frequency band、vector polarization、port metric与观测区域。只有 full-state/full-data closure 随 (r_i) 收敛且 end-to-end 成本优于 matched solver，才能建立实际 compression。 |
| C39 | A1 1374–1514；A1 4231–4277；A1_2 2539–2575 | 从 Born/Ewald 主理论切换到 full-wave 主线 | **版本冲突必须显式保留** | A1 前半建议 Born 为理论核心、full wave 只作 mismatch stress test；后半和 A1_2 改为 full-wave 主线。扩实验必须分 matched Born、full-wave generated、full-wave inverted 三类，不得把任一类结果写成另一类证据。 |
| C40 | A1 3188–3426；A1_2 1548–1604 | 保矩 split 经任意 smooth full-wave map 保持四阶 data change | **局部阶次成立，原创性未立** | 四阶来自 material perturbation 的矩匹配，再与可微 forward map 复合；这是 generic smooth composition。真正与观测有关的是 noise-whitened、对 incumbent+nuisance refit 投影后的 signed leading coefficient；若该 projector 后为零，须继续到更高阶。 |
| C41 | A1 2577–2704；A1_2 3526–3635 | dark self-energy 提供 exact elimination | **标准 Schur/Feshbach 恒等式** | exact reduction 必须同时处理 reduced operator、right-hand side 和 output correction；算法贡献只能来自可控近似、复用成本或统计校准，不能来自重新命名恒等式。 |
| C42 | A1 2856–3174, 4101–4118；A1_2 4285–4355 | NN / GRPO / soft partition 是主线增强模块，部分被标绿色 | **支持模块，claim 未清** | generic corrective NN、learned optimizer、soft partition 与 T-matrix+NN 均有近邻。若 GNN 只加速 residual iteration，应按 learned solver 评价；若输出散射块或材料，则另过 physics 与 identifiability gate。A1 的颜色不是证据等级。 |

## NN/GNN 物理保持的最小合同

如果采用 NN-native scattering graph，建议先明确它属于以下哪一类，不能混写：

1. **learned iterative solver**：节点/边仍由 projected full-wave 方程给出，GNN 只预测 damping、preconditioner 或 residual update。每一步都应保留 (r=f-(I-H)b)，并与 block-Jacobi、GMRES、exact reduced solve 作相同步数/容差/时间比较。
2. **local (T_i) surrogate**：输入材料、频率和 shape，输出固定 port convention 下的 (T_i)。必须测 unseen frequency/contrast/aspect ratio/resonance，并用独立 full-wave local solve验真。
3. **dark/self-energy surrogate**：输出被截断模态对保留空间的闭合修正。必须测它是否改善 full-state/data closure；不能用 reduced residual 自证。
4. **inverse GNN**：从外部数据估计材料/节点。必须承认 realization gauge 与有限孔径非唯一性，输出 posterior/uncertainty 或等价类，而非把 latent (T_i) 当材料真值。

四类共同约束：

- node permutation equivariance 不够；还需对局部 port basis 的 unitary/gauge 变换协变，并处理 rigid rotation 下的 Wigner-(D)/vector transformation；
- analytic Green/translation edge 应保留 reciprocity、outgoing phase和几何，不宜由无约束 message MLP 替代；
- paired-real 网络若任意混合实虚通道，会改变复相位结构；需声明允许的 complex nonlinearity和物理等变性；
- local pointwise passive 不等于跨频 causal，也不自动保证截断后的 global interconnection正确；
- 验收量必须同时含 full-state projection error、full VIE receiver error、physics residual、能量/互易/因果缺陷与端到端成本。

## 近邻来源 ledger

| 近邻 | 本轮证据 | 能支持什么 | 仍不能支持什么 |
|---|---|---|---|
| Chen, Zhong & Agarwal, *Subspace Methods for Solving Electromagnetic Inverse Scattering Problems* (2010/2011), DOI `10.4310/MAA.2010.v17.n4.a6` | **A，本地全文** | 点散射体 Foldy–Lax；SOM 的 deterministic/ambiguous current；TSOM 用 internal propagation 强方向进一步限制 ambiguous current。 | 不直接证明 overlapping Gaussian full local (T) 的有限-port closure，也不支持把 second fold 改名为任意 dual Jacobian ratio后主张新意。 |
| Kanaun, PIER B 21 (2010), DOI `10.2528/PIERB10030803` | **A，本地全文** | Gaussian basis 的 3-D dielectric VIE、解析矩阵元、规则网格 Toeplitz/FFT。 | 不提供 adaptive node-local T network、GNN 或 inverse identifiability。 |
| Uysal & Akleman, *Domain decomposition method with generalized scattering matrix for radiowave propagation analysis*, AEU 171 (2023) 154884, DOI `10.1016/j.aeue.2023.154884` | **B，出版社摘要** | HG coefficients、GSM、cascaded two-port subdomains、重复 forward 的矩阵复用。 | 未读全文，不能作“没有 adaptive inverse/material/model-order”这种负面全文断言。 |
| Ishida, *Reconstruction of a Dielectric Cylinder with the Use of the T-Matrix and the Singular Value Decomposition* (2010), DOI `10.1587/transcom.E93.B.2595` | **B，出版社摘要** | measured data 转 T-matrix；SVD 分 radiating/nonradiating currents；恢复 dielectric cylinder。 | 不等价于 Gaussian local-node network，但已经占据 `T-matrix inversion + SVD current split` 的宽 claim。 |
| Ye, Chen, Zhong & Song, IEEE TAP 61 (2013), DOI `10.1109/TAP.2013.2258878` | **B/C，机构摘要与 DOI** | 2-D TM 下同时重建 dielectric 与 PEC，使用 dipole/monopole T-matrix elements。 | 当前摘要不足以逐公式比较 proposed observable-subspace de-embedding。 |
| Zhang et al., IET Radar Sonar Navig. 17 (2023), DOI `10.1049/rsn2.12419` | **A，出版社全文页面** | CSI 与 T-matrix combined parameter model；APU-Net 交替输出 dielectric contrast source 和 T-matrix；含物理 forward module。 | 不是 Gaussian graph，但已否定 generic `T-matrix + physics-inspired NN` 新颖性。 |
| adaptive multipole truncation / reduced T-matrix / GNN wave solver / network tomography | **U，本轮未全文闭合** | 是 port-rank、learned graph 与内部 realization identifiability 的直接邻域。 | 在这些邻域检索完成前，A1_2 4347–4355 的加粗“潜在新东西”只能是研究问题。 |

出版社/机构入口：

- Uysal–Akleman 2023: <https://doi.org/10.1016/j.aeue.2023.154884>
- Ishida 2010: <https://doi.org/10.1587/transcom.E93.B.2595>
- Ye et al. 2013: <https://doi.org/10.1109/TAP.2013.2258878>
- Zhang et al. 2023: <https://doi.org/10.1049/rsn2.12419>

## 最小可证伪矩阵

| Gate | 必须做的反例/实验 | 失败判据 |
|---|---|---|
| G1 完整算子→有限 ports closure | 对固定 truth，以 (r_i) 递增比较 full VIE current、receiver field和能量；扫 overlap/gap、(k\sigma)、contrast、loss和近共振。 | reduced residual 小但 full-state或receiver误差不收敛；误差随 overlap/共振爆炸。 |
| G2 overlap/factorization gauge | 用两种不同重叠 Gaussian decomposition 表示同一 (chi)，比较 full solver与网络输出；另用非重叠 patch/cluster 对照。 | 相同 (chi) 因内部拆分不同产生显著外部预测差，且不能由端口增阶消除。 |
| G3 Hermite↔multipole 映射 | 逐个 (H_0,H_1,H_2,H_4) tangent 计算各 VSWF (l,m,mathrm{TE/TM}) 能量矩阵，随旋转/背景/contrast 复测。 | 能量广泛混阶或 ordering 随场景改变；则禁止以 (H_n) 直接命名 port order。 |
| G4 Gram invariance | 对同一 span 施加随机可逆 rescale/mixing；比较 rank-aware whitening 后奇异值、物理 projector、split decision；扫 near-rank threshold。 | 结论随坐标 scaling 或合理阈值显著变化。 |
| G5 projector detectability | 对 split 的 whitened signed (c_4) 先投影 incumbent+nuisance tangent，再做一侧、多方向搜索校准；同时报告 raw 和 projected。 | raw 响应显著但 projected 响应落入噪声，或 false-birth 失控。 |
| G6 external identifiability | 显式构造内部 similarity/gauge 或不同 (T,U) realization 产生同一 (M)；再逐步固定 geometry/U/port convention。 | 在现实 partial Tx/Rx 下仍有不同材料/网络等价，而方法输出单一“真” (T_i)。 |
| G7 split update correctness | 真正 remove parent + add two children，更新 T/U/B/C；同时比较错误的 append-only。扫 (d\to0)。 | append-only 被当作保矩 split；或外部 (O(d^4)) 同时内部条件数/solver error恶化而未记录。 |
| G8 Twofold 增益 | matched DOF/cost 下比较 classical TSOM、restricted (G_DQ_N)、regularized GN/GSVD、普通 low-rank tangent和随机子空间。 | second fold 只给高 ratio，却不改善 held-out data、state residual、材料误差或稳健性。 |
| G9 NN 物理性 | unseen geometry/frequency/contrast；检查 full-physics residual、passivity metric、reciprocity metric、causal time response/K–K、旋转和 port-basis gauge等变性。 | 图像指标提升但任一独立物理残差系统性失败。 |
| G10 计算优势 | 同精度下比较 setup+forward+adjoint+update 的 wall time、memory、factorizations、端口增长，基线含 full VIE、GMRES/FMM或项目可用等价。 | 仅 (R<N) 或小 inverse 更快，但端到端时间/内存无优势。 |
| G11 新颖性 | 补齐上表 U 类全文；逐项对比 adaptive multipole、model reduction、network inverse/system ID、wave GNN、adaptive RBF inverse。 | 最近邻已经含同一可观测性驱动 port/node增阶机制，或差异只剩 Gaussian 命名。 |

## 独立 hard gates

- **数学定义 gate：未过。** Hairy-ball、rank 条件、(T^{-1}) 可逆性、两处 Neumann 收敛、T/S 区分和真实 split update 必须先修正。
- **finite-port continuum/closure gate：未过。** 尚无 full local operator 到有限 HG/VSWF ports 的截断界或收敛证据。
- **overlap gate：未过。** additive-potential exactness 与 isolated n-port translation 的前提尚未统一。
- **identifiability gate：未过。** (M=CKB) 只给外部 transfer；内部 realization、local (T_i) 和材料的唯一性未建立。
- **NN physics gate：未过。** 当前只有原则，没有满足 passivity/reciprocity/causality/basis-gauge equivariance 的具体参数化和独立验证。
- **evidence/advantage gate：未过。** 无 port-rank law、full VIE closure、matched end-to-end cost或外部分布泛化结果。
- **novelty gate：未过。** HG+GSM、T-matrix inverse/SVD、T+SOM 与 T+physics-inspired NN 均已有明确近邻；GNN/adaptive-port/network-de-embedding 邻域尚未全文闭合。

这些 gate 相互独立。修正公式不关闭数值闭合，数值闭合不关闭可识别性，可识别性也不关闭新颖性。
