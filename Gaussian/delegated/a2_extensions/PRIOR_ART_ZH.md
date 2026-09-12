# A2 两条新路线：本轮最近邻核查

检索日期：2026-09-11。ScholarQA 两条 collect 查询都返回 429（同目录 scholarqa.json）；随后核查作者/出版社原始来源。以下是有边界的最近邻核查，不是全领域穷尽检索。对话原文中的引用占位符不算已经核验的文献。

| 来源与等级 | 直接支持什么 | 对本项目的限制 |
|---|---|---|
| [Gharsalli 等，Variational Bayesian inversion for microwave breast imaging](https://cames-old.ippt.pan.pl/index.php/cames/article/view/38)，出版社摘要级 | 非线性微波逆散射、有限均匀材料区、Gauss–Markov–Potts 先验、联合后验、可分变分近似以及 CSI 对照均已有工作 | 不能宣称首次概率材料向量或首次 Bayesian 微波成像。完整 PDF 本次抓取超时；不据摘要断言其所有标签更新细节 |
| [Hosseini 与 Sra，Matrix Manifold Optimization for Gaussian Mixtures，NeurIPS 2015](https://proceedings.neurips.cc/paper/2015/file/dbe272bab69f8e13f14b405e038deb64-Paper.pdf)，出版社全文级 | SPD 流形上的 Gaussian 参数优化、度量和数值对照已有系统研究；直接换流形并不自动更快 | 不是电磁反演论文，但足以排除“Gaussian + manifold”本身的原创性；不能把其统计模型的几何凸性移到全波损失 |
| [Golub 与 Pereyra，Separable nonlinear least squares: the variable projection method and its applications，2003](https://doi.org/10.1088/0266-5611/19/2/201)，本轮只确认 DOI，全文访问被拒 | 作为待精读的 VarPro 近邻 | 本轮不靠未读全文支持具体定理。局部条件 LS 的投影等式在 CORE_AUDIT 中独立给出；全波 Gaussian 幅度一般不是线性可消元参数 |

Gharsalli 网页的卷期是 21(3–4), 2014，页面发布/引用字段显示 2017；论文正式书目年份需在投稿前按出版版本核对，不能悄悄忽略这个元数据冲突。

本轮允许的创新候选是：控制联合概率近似与全波模型误差如何改变结构决策；以及同一物理度量下有限 secant 如何在旧参数重新拟合后仍提供新信息。它们目前是需要证明和比较的机制，不是已成立的首创结论。原 A2 的认证降阶/信赖域近邻仍有效，见原有文献表。
