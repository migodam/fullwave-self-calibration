# A1_2 近邻文献核查：NN-native scattering matrix / GNN

日期：2026-09-11  
任务性质：有界 prior-art evidence extraction；**不作最终新颖性裁决**。  
主线锚点：`Gaussian/Theory/A1_2.md` 第 11、26–29 节。未修改 `A1.md` 或 `A1_2.md`。

## 检索边界与证据等级

本轮先按 `scholarqa-research` 流程调用已安装 ScholarQA CLI，对 6 个预注册查询做 collection。Semantic Scholar 在 12 个操作中的 11 个返回 HTTP 429；原始输出保存在 [`scholarqa_collect.json`](./scholarqa_collect.json)。依据既有用户授权，没有刷新或重复轰炸该接口，转而核查出版社、arXiv、机构库、PubMed Central 和作者公开稿。因而这里能说明“找到了哪些直接近邻”和“哪些宽泛主张已被占据”，不能证明不存在未检出的论文。

- **Tier A**：读到开放全文或作者公开全文，并在正文中定位方法。
- **Tier B**：只读到出版社/机构的题录与摘要；只用于摘要明确陈述的内容。
- **Tier C**：搜索片段或二手转述；不用于肯定性技术结论。本报告的主要判断不依赖 Tier C。

## 结论先行

`A1_2.md` 对 2023 Hermite–Gaussian + generalized scattering matrix、2010/2013 T-matrix inverse、2022/2023 T-matrix + neural hybrid 的提醒基本准确。更需要补上的近邻有三组：

1. **NN-native 迭代电磁求解器已存在。** Shan 等的 PhiGRL 让 GNN 依据 CFIE 残差迭代修正候选表面电流；Zhu 等 2026 预印本把 MLFMA 的近场块映射为局部 graph message passing，把远场耦合映射为层次模块。因而“用 GNN/message passing 代替末端 U-Net 来做 full-wave scattering solve”本身不能作为新意。
2. **“粒子是节点、物理相互作用是边”的散射代理已存在。** Lamb–Gentine 以相互作用球团簇为图，用 MSTM 数据训练 message-passing GNN 输出光学量；Davidi–Pakizeh 2026 更直接把强耦合 plasmonic nanoparticles 当节点，把极化率和物理相互作用编码进边。它们尚未学习可级联的 local T/S blocks，也没有 adaptive Gaussian move/split/merge inverse，但已明显压缩“scatterer graph”这一层的新颖空间。
3. **学习 scattering matrix 与施加物理约束均有先例。** Jing 等 2023 从几何预测一般散射矩阵；其能量守恒/互易性是训练数据中学到后做的经验误差检查，并非硬保证。Torun 等 2020 已在 NN 中设计 causality/passivity enforcement layers；Lilja 等 2026 则用 QNM 展开构造保证能量守恒与因果性的 S-matrix 网络，并允许加入对称性等约束。故“NN 输出 passive/reciprocal S”也不能单独宣称新。

在这次有界检索中，**没有找到同时具备以下整条机制的论文**：连续可移动、可保持矩量地 split/merge 的 Gaussian material primitives；每个 primitive 的可学习 local modal T/S block；显式 translation-coupled multiple-scattering graph；observable/dark-subspace de-embedding；以及由逆问题证据驱动的 feature-order/model-order growth。这个“组合空缺”只应写成待进一步查新的候选差异，不能写成已证明的新颖性。

## 与 `A1_2.md` 主张逐条核对

| A1_2 主张/方向 | 核查结果 | 原文证据 | 对可写主张的影响 |
|---|---|---|---|
| 2023 年已有 HG basis + GSM + domain decomposition | **支持，Tier B** | Uysal–Akleman 将问题域表示为级联 two-port networks，用 Hermite–Gaussian basis 得到各子域 scattering relations，再以矩阵运算求边界场；出版社正文预览还明确写到子域矩阵可复用、无需逐波追踪多次反射。 | “Gaussian/HG modal basis + scattering network”必须视为已有。该文是 radiowave forward domain decomposition，不是 adaptive inverse material splatting。 |
| 2010 年已有 measured data → T-matrix → SVD radiating/nonradiating currents → dielectric reconstruction | **支持，Tier B** | Ishida 的官方摘要逐步陈述这两个阶段及 SVD 分解。 | “从测量反演 T-matrix / 可见-不可见电流分解”不是空白；Gaussian observable/dark de-embedding 必须说明具体新增结构。 |
| 2013 年 T-matrix + SOM 重建 dielectric/PEC | **支持，Tier A** | Ye 等正文明确把 SOM 扩展到 T-matrix modeling，并以 multipole truncation 在网格子单元上建模；使用 synthetic 与 Fresnel experimental data。 | “T-matrix + SOM”及局部 multipole-order accuracy/cost tradeoff 已有。 |
| 2022/2023 年 T-matrix 与 neural inverse 混合 | **支持，Tier A/B** | Song 等先 BP 得到零阶 T-matrix coefficient image，再由 attention-assisted pix2pix/U-Net 精修；Zhang 等把 CSI 与 T-matrix 合成联合参数模型，以 APU-Net 交替更新 dielectric contrast 和 PEC T-matrix。 | 末端图像网络与 physics-inspired alternating network 都已有。新方案应把可学习量放在 local T/S block 与物理消息传递中，而不是只替换最终图像头。 |
| NN-native scatter solver / learned iterative solver | **已有直接近邻，Tier A** | PhiGRL 模拟 fixed-point iteration：每步显式算 CFIE residual/matrix-vector product，GNN 预测候选解修正直至收敛。SOM-Net 也把 SOM 迭代展开为交替 current/permittivity 更新。 | “learned iteration”或“GNN 求 scattering”不能单列为创新。 |
| GNN 直接表达局部/远程 EM coupling | **已有直接近邻，Tier A/B** | U-PINet 以 surface-mesh graph 的 learned message kernels 表达 MoM near-field block，并以 octree multi-scale module 表达远场；Davidi–Pakizeh 将 plasmonic particles 设为节点、物理相互作用为 edge features。 | Gaussian primitive、modal ports、T/S composition、adaptive topology 与 inverse objective 才可能构成区分。 |
| learned general scattering matrix | **已有，Tier A** | Jing 等以 polygon geometry 输入、complex T-matrix 输出，报告约 7,488× FEM 推理速度；网络为 U-Net-like convolutional architecture。 | “NN 预测一般 scattering matrix”不新。 |
| passive/reciprocal learned local S blocks | **组成部分分别已有；精确组合未检出** | Torun 等以 CEL/PEL 在训练中保证 causal/passive broadband S-parameters；Lilja 等以 QNM expansion 保证 energy conservation/causality并允许手工 symmetry constraints；Jing 等仅经验性观察 reciprocity/energy errors 小。 | 不能声称首次物理约束 S 网络。若主张 local Gaussian block 的 joint passivity+reciprocity，应给出严格参数化、端口归一化和截断通道解释。 |

## 最接近的架构近邻

### 1. PhiGRL：GNN 残差迭代求 CFIE

Shan, Li, Yang, Xu 的期刊版本为 *Solving Combined Field Integral Equations With Physics-Informed Graph Residual Learning for EM Scattering of 3-D PEC Targets*, IEEE TAP 72(1), 733–744, DOI `10.1109/TAP.2023.3331262`。开放的 URSI 2023 全文摘要说明：

- 以 triangular mesh 形成图，处理变化的未知量数；
- 模拟 fixed-point iteration；
- 每一步显式计算 CFIE residual 和 matrix-vector multiplication；
- GNN 预测 candidate solution correction，迭代至 convergence。

这是“物理线性系统 + learnable graph update”的直接先例。它的节点是表面网格未知量，边不是 Gaussian-to-Gaussian translation operator，输出也不是每个节点的 local modal T/S block；任务是 PEC forward solve。

同一作者 2024 年的 ACES-China 短文还直接以 message-passing GNN 从 triangular-mesh graph 预测 3D PEC surface currents（DOI `10.1109/ACES-China62474.2024.10699762`）。这进一步确认 GNN surface-current surrogate 已经是明确先例；本轮只取得题录与摘要，按 Tier B 使用。

### 2. U-PINet：MLFMA 结构映射到 GNN/层次网络

Zhu, Peng, Alexandropoulos, Wang 的 arXiv:2508.03774 v5（2026，标注为投稿中）更接近目标架构：near-field graph encoder 对相邻 surface elements 做多轮 message passing，far-field branch 依照 MLFMA octree 做 coarse-to-fine 耦合，网络输出 induced surface current，并用离散 EFIE residual 训练，不需要 reference current labels。作者还把它作为 MLFMA-GMRES initializer 测试。

它显著削弱“把 EM solver 的局部作用写成 learned message passing、远程作用写成层次传播”的新颖性。仍有的差别是：PEC surface mesh、单次前向 current surrogate、固定对象/重复照射场景；没有 Gaussian local T/S ports、连续 move/split/merge，也没有 inverse de-embedding。

### 3. 粒子图散射代理：Lamb–Gentine 与 Davidi–Pakizeh

Lamb–Gentine 2023 用 MSTM 产生标签，把 black-carbon aggregate 建成 interacting spheres 的图；节点含球坐标与折射率，边含距离，Interaction Network 通过 message passing 后 graph pooling，预测总 extinction/scattering/absorption、asymmetry 和 angle-resolved phase-matrix elements。其 zero-shot 大尺寸测试表明图模型可跨节点数泛化，但论文也明确使用独立 zero-shot test；zero-shot validation 曾参与架构选择。

Davidi–Pakizeh 2026（Optics Express, DOI `10.1364/OE.586635`）从官方摘要可确认：强相互作用 plasmonic nanoparticles 作为 nodes，polarizability 与物理 interaction 进入 edge features，训练数据混合 coupled-dipole approximation 与高精度数值模拟；当前验证集中在 dimer/trimer。该文是当前检索中对“primitive graph + physics-aware edges”最接近的已发表近邻。

二者的输出是 aggregate optical response，而非可解释、可级联并供逆问题 de-embed 的 node-local T/S matrices。这个输出层级差别必须成为方法与实验的一部分，不能只靠命名区分。

### 4. 端口矩阵学习与约束

Jing 等 2023 直接学习低对称二维介质体的 complex scattering matrix。正文明确说物理律没有作为条件放入主网络；其训练损失是 MSE，随后在随机多边形上检查能量守恒与互易误差。因此该文的“inherently satisfies”是数据分布内的经验结论，不是对任意输入的结构性保证。

Torun 等 2020 的 causality enforcement layer 和 passivity enforcement layer 已把 Kramers–Kronig 与 S-parameter singular-value 条件加入端到端训练。Lilja 等 2026 的 QNM-Net 更进一步，把 scattering spectrum 通过 QNM expansion 参数化，给出 energy-conservation 与 causality guarantee，并在可用时加入 losslessness/symmetry/mode-count constraints。联合“local Gaussian + reciprocity + contractivity”仍可能有实现差异，但单项约束没有新颖空间。

## 对 NN-native 路线的证据约束

如果后续架构仍以最终 contrast/permittivity image 为主要 NN 输出，会直接落入 Song 2022、Zhang 2023、SOM-Net 2022 和 Stylianopoulos 2025 的拥挤区域。后者甚至已有 graph-attention backbone，却仍由 U-Net head 做最终重建；这正是用户希望避开的形态。

要把路线与这些近邻实质分开，论文级实现至少需要让下列量成为网络计算图的主体：

1. 节点状态是 incident/outgoing **modal coefficients**，而非像素 latent；
2. 节点模块输出或参数化 Gaussian-specific local `T_i`/`S_i`；
3. 边执行有物理锚点的 translation/coupling，迭代或隐式求解多次散射状态；
4. inverse loss 作用于外部端口响应，并显式处理 observable/dark subspace；
5. move/split/merge 与 Hermite feature-order growth 改变图拓扑和端口阶数；
6. 图像 decoder 仅用于可视化或辅助监督，不能承担核心逆映射。

这六点是从近邻差异推导出的**定位条件**，不是文献已证明的性能优势。

## 科学与表述风险

1. **简单写 `S^T=S` 有归一化风险。** Byrnes–Foreman 对含 evanescent components 的 vector scattering/transfer matrices 推导了更一般的能量、互易、时间反演约束；普通 far-field unitary/symmetric 形式只是特例。Gaussian/Hermite ports 若含非正交、复数或倏逝通道，需要先固定 power normalization、incoming/outgoing pairing 和 reciprocity metric。
2. **reduced S 的 contractivity 不等于 local model 正确。** 截断端口产生的“泄漏”可与材料耗散、未建模远场和基函数误差混淆。passivity projection 只能排除非物理解，不能证明截断精度或可辨识性。
3. **端口物理约束不是训练泛化保证。** Jing 的经验守恒、Torun 的训练期 enforcement、Lilja 的结构参数化是三种不同强度的结论，引用时不能互换。
4. **GNN 的节点数外推已有先例。** Lamb–Gentine 与 Kuhn 等均展示跨图尺寸/连接的泛化动机或实验；“可变 Gaussian 数量”只能作为必要属性，不能独立宣称创新。
5. **learned solver 与 inverse learner 必须分开评价。** PhiGRL/U-PINet 主要解决 forward current；SOM-Net/APU-Net 解决 inverse。新方案若同时做两者，需要分别报告 forward residual/field error、inverse recovery、迭代稳定性和拓扑变化外推。
6. **2026 近邻改变了时间边界。** Davidi–Pakizeh 已正式发表；U-PINet 是 2025 提交、2026 更新的预印本。投稿时必须重新核查其版本与审稿状态。

## 主要来源

1. A. Uysal, F. Akleman, “Domain decomposition method with generalized scattering matrix for radiowave propagation analysis,” *AEU* 171, 154884 (2023), [DOI](https://doi.org/10.1016/j.aeue.2023.154884), [publisher abstract/preview](https://www.sciencedirect.com/science/article/abs/pii/S1434841123003588). Tier B.
2. K. Ishida, “Reconstruction of a Dielectric Cylinder with the Use of the T-Matrix and the Singular Value Decomposition,” *IEICE Trans. Commun.* E93-B(10), 2595–2600 (2010), [official abstract](https://globals.ieice.org/en_transactions/communications/10.1587/transcom.E93.B.2595/_p). Tier B.
3. X. Ye, X. Chen, Y. Zhong, R. Song, “Simultaneous reconstruction of dielectric and perfectly conducting scatterers via T-matrix method,” *IEEE TAP* 61(7), 3774–3781 (2013), [DOI](https://doi.org/10.1109/TAP.2013.2258878), author manuscript in `pdfs/`. Tier A.
4. R. Song et al., “Learning-Based Inversion Method for Solving Electromagnetic Inverse Scattering With Mixed Boundary Conditions,” *IEEE TAP* 70(8), 6218–6228 (2022), [DOI](https://doi.org/10.1109/TAP.2021.3139645), author manuscript in `pdfs/`. Tier A.
5. Q.-q. Zhang et al., “An inverse scattering reconstruction method for perfect electric conductor-dielectric hybrid target based on physics-inspired network,” *IET Radar Sonar Navig.* 17(8), 1286–1298 (2023), [publisher page](https://ietresearch.onlinelibrary.wiley.com/doi/10.1049/rsn2.12419). Tier B for method summary.
6. Y. Liu et al., “SOM-Net: Unrolling the Subspace-based Optimization for Solving Full-wave Inverse Scattering Problems,” *IEEE TGRS* 60 (2022), [arXiv](https://arxiv.org/abs/2209.03567). Tier A.
7. T. Shan et al., “Solving Combined Field Integral Equations With Physics-Informed Graph Residual Learning for EM Scattering of 3-D PEC Targets,” *IEEE TAP* 72(1), 733–744, [DOI](https://doi.org/10.1109/TAP.2023.3331262); [open URSI method summary](https://www.ursi.org/proceedings/procGA23/papers/YSASummaryShanTao.pdf). Tier A for the open summary; journal metadata cross-checked.
8. K. D. Lamb, P. Gentine, “Zero-shot learning of aerosol optical properties with graph neural networks,” *Scientific Reports* 13, 18777 (2023), [open article](https://www.nature.com/articles/s41598-023-45235-8). Tier A.
9. L. Kuhn, T. Repän, C. Rockstuhl, “Exploiting graph neural networks to perform finite-difference time-domain based optical simulations,” *APL Photonics* 8, 036109 (2023), [open institutional copy](https://publikationen.bibliothek.kit.edu/1000157040/150479803). Tier A.
10. S. Bakirtzis et al., “Solving Maxwell's equations with Non-Trainable Graph Neural Network Message Passing,” [arXiv:2405.00814](https://arxiv.org/abs/2405.00814) (2024). Tier A; preprint.
11. R. Zhu et al., “A Physics-Informed Hierarchical Neural Network for Microwave Scattering Analysis of 3D PEC Targets,” [arXiv:2508.03774 v5](https://arxiv.org/abs/2508.03774) (2025/2026). Tier A; submitted preprint.
12. A. A. Davidi, T. Pakizeh, “Physics-aware graph neural networks for optically interacting plasmonic nanoparticle assemblies,” *Optics Express* 34(9), 16790–16801 (2026), [DOI](https://doi.org/10.1364/OE.586635), [PubMed abstract](https://pubmed.ncbi.nlm.nih.gov/42271735/). Tier B in this review because publisher PDF automation was blocked.
13. K. Stylianopoulos et al., “Graph-CNNs for RF Imaging: Learning the Electric Field Integral Equations,” [arXiv:2503.14439](https://arxiv.org/abs/2503.14439) (2025). Tier A.
14. Y. Jing et al., “A deep neural network for general scattering matrix,” *Nanophotonics* 12(13), 2583–2591 (2023), [PMC full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC11501312/). Tier A.
15. H. M. Torun et al., “Causal and Passive Parameterization of S-Parameters Using Neural Networks,” *IEEE T-MTT* 68(10), 4290–4304 (2020), [institutional record/abstract](https://avesis.metu.edu.tr/yayin/147b26d5-9f2d-4e5d-8837-dd3c5e3bfbac/causal-and-passive-parameterization-of-s-parameters-using-neural-networks). Tier B.
16. V. A. Lilja et al., “A general framework for knowledge integration in machine learning for electromagnetic scattering using quasinormal modes,” *Laser & Photonics Reviews* (in press, 2026), [arXiv:2509.06130 v2](https://arxiv.org/abs/2509.06130), [institutional record](https://research.chalmers.se/en/publication/551287). Tier A.
17. N. Byrnes, M. R. Foreman, “Symmetry constraints for vector scattering and transfer matrices containing evanescent components,” *Phys. Rev. Research* 3, 013129 (2021), [open publisher article](https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.3.013129). Tier A.
18. T. Shan, “Solving Electromagnetic Scattering of 3D PEC Targets Based on Graph Neural Networks,” ACES-China (2024), [DOI](https://doi.org/10.1109/ACES-China62474.2024.10699762), [author publication list](https://iemcs-lab.github.io/). Tier B.

机器可读的 claim/source 映射、访问级别、本地文件与 SHA-256 见 [`claim_source_ledger.json`](./claim_source_ledger.json)。
