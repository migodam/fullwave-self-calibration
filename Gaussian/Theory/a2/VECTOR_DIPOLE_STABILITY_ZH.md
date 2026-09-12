# A2 的三维向量扩展：限于声明的点偶极离散模型

这不是把二维 δh 公式移植成 Maxwell 定理，而是对已有三维向量 DDA 模型重新推导。DDA 与辐射反作用修正本身是成熟方法，参见 Draine–Flatau 的[官方作者指南](https://arxiv.org/abs/1202.3424)及其原始文献。以下耗散计算是本地整理的有限系统应用，不宣称原创。

## 模型和结论

采用 e^-iωt、无磁均匀背景。每个位置为三分量电偶极 p_i，标量材料 χ_i、单元体积 V_i>0。K_ij 为自由空间 k² 电场 Green 张量（i≠j），自块为零。采用明确的 Clausius–Mossotti 加辐射修正极化率：

$$
\alpha_i^{-1}=\frac1{3V_i}+\frac1{V_i\chi_i}-i\frac{k^3}{6\pi},
\qquad Lp=E,\quad L=\operatorname{diag}(\alpha_i^{-1}I_3)-K.
$$

此处 α 单位为体积，χ 是相对介电对比度。仅讨论 χ_i≠0；已知空白单元消元。正损耗 Imχ_i>0。写 β=k³/(6π)。出射 Green 张量虚部的角向表示给出

$$
\mathcal G_{ij}=\frac{k^3}{16\pi^2}\int_{\mathbb S^2}
(I-\hat s\hat s^\top)e^{ik\hat s\cdot(x_i-x_j)}d\Omega,
\qquad \mathcal G=\operatorname{Im}_H K+\beta I.
$$

对任意复向量偶极组 z，积分被积项为横向投影后的模平方，所以 $\mathcal G\succeq0$；对角由 $\int(I-\hat s\hat s^\top)d\Omega=8\pi I/3$ 得到 βI。于是辐射自项精确抵消：

$$
-\operatorname{Im}_H L
=\mathcal G+\operatorname{diag}\left(
\frac{\operatorname{Im}\chi_i}{V_i|\chi_i|^2}I_3\right)
\succeq\Gamma\succ0.
$$

由 $|z^*r|\ge-\operatorname{Im}(z^*Lz)\ge\|z\|_\Gamma^2$ 和加权 Cauchy–Schwarz，$z=L^{-1}r$ 满足

$$
\|L^{-1}r\|_\Gamma\le\|r\|_{\Gamma^{-1}}.
$$

因此 A2 的 primal/adjoint 输出残差恒等式及加权上界可用于这个三维**向量离散系统**。三个场分量、横向投影和辐射自项均不可删除。它不保证该点偶极模型接近连续介质的精度。

## 严格边界

1. 极化率、自项与非对角核必须配套；换成其他 LDR、非球形单元或自作用修正，需要重新计算耗散部分。
2. 零损耗时本 Γ 为零，证书失效不表示散射问题无解；辐射 Gram 本身可能提供额外有限系统稳定性，但本轮没有廉价且统一的正下界。
3. 这是标量各向同性材料、三分量向量场；并未处理张量材料、磁性、分层背景或连续核奇性。
4. α 对 χ 非线性。不能把各分量的 α 简单相加，再套用“χ 可加”高斯组件网络；必须从统一状态中的可加材料势重新推导。
5. 该界控制当前材料的离散状态，不认证材料恢复、端口可压缩性或任意形状连续 Maxwell 成像。

## 实际检查

`code/a2/vector_certificate_probe.py` 独立实现三维 dyadic 核；N=3³、4³ 单元（81、192 向量未知数），k=4/8，有损 Gaussian 与两椭球＋已知有损背景，两个局部状态秩，共 16 次输出界覆盖检查。辐射 Gram、耗散不等式与校正场覆盖均在浮点容差内成立。结果在 `runs/a2/vector_certificate/results.json`。

这是新增的向量条件检查。此前 N9³/N11³ 的三维成像/分裂探针保留在 `runs/maxwell`，本轮没有重新执行，也没有由这些小规模数据宣称完成 3D TAP 验收。
