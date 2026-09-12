# A2 最近邻与主张边界

2026-09-11。本轮 ScholarQA 两个查询均被 Semantic Scholar 429 限流，记录见 scholarqa.json；按授权切换到原始出版源。没有据检索无结果宣称首次。

|候选主张|已核查来源|本轮可作的判断|
|参数空间与状态空间联合自适应并用误差控制|Kartmann 等，2024，[出版全文](https://link.springer.com/article/10.1007/s44207-024-00002-z)|论文已联合参数空间自适应和认证降阶、信赖域。因此这个宽泛组合不是 Gaussian 的独有创新。该文具体假设不能直接迁移到非线性 Gaussian Maxwell。|
|伴随修正、后验误差与自适应信赖域|Keil 等，2021，[原始期刊记录与全文入口](https://www.numdam.org/articles/10.1051/m2an/2021019/)|已有 rigorous NCD-corrected reduced-basis 优化及误差界，不能将 A2 的通用残差恒等式包装为首次。|
|利用内部传播压缩数据弱电流空间|Chen 的 [UCI 官方讲座摘要](https://www.math.uci.edu/taxonomy/term/103?page=1)，并参照 Xu 等作者稿|这正是 TSOM 已有思想，Gaussian 的材料切空间拉回需要说明新增算法收益。|
|TSOM 边界改善、FFT 降成本与 Fresnel 检验|Xu 等，[作者公开稿](https://renchengsong.github.io/papers/TGRS2329032.pdf)，DOI 10.1109/TGRS.2014.2329032|乘性正则 FFT-TSOM 已研究边界保持与实际数据。边界改善不能只比较平滑 Gaussian，需包括成熟边界正则。|
|HG 与广义散射矩阵、T-matrix SOM、原生 GNN|已完成的 [A1_2 近邻表](../a12_literature/PRIOR_ART.md)|保留原始来源等级；本轮 n-port 改写不解除这些近邻约束。|

用户指定的 2023 论文题录线索为 Tiantian Yin、Kai Tan、Xudong Chen，*Far-Field Approximation Learning Method for Millimeter-Wave Short-Range Imaging*，DOI 10.1109/TAP.2023.3242424。本轮检索仅取得题录线索，不据此声称已经完成其技术逐项对照。

值得保留的潜在区别是：用可计算的状态误差和材料方向残差共同决定哪些 Gaussian 更新值得计算、哪些新增结构被数据支持，并在实际大参数问题中取得净成本和错误结构控制收益。这是待验证的贡献目标；当前验证尚未证明它。

补充题录交叉核对：[IEEE 2023年4月官方目录](https://engage.ieee.org/rs/756-GPH-899/images/T-AP_Apr2023.pdf) 列出 Yin、Tan、Chen 的该文。它只加强题录核对，不是论文全文技术证据。
