# A2 补充：独立单元积分核的耗散下界

本轮本地推导；尚无原创性主张。用途是消除“Pro 点采样核证书不能自动迁移到新积分核”的实现缺口。二维标量、均匀方形网格、出射 Hankel-1、张量积对称正权 Gauss–Legendre 源单元积分；不适用于任意非均匀网格或三维 Maxwell。

正式号约定为 $e^{-i\omega t}$，被动材料 Imχ≥0；未知量是等面积单元上的电流密度、测试采用中心配点。伴随和范数使用均匀单元质量白化后的 Euclidean 坐标。改变时间号、非均匀质量、测试规则或核定义都需要重新推导。接收白化可以进入 C，但不能改变状态证明所用的内积。

设节点为 $x_i$，单元边长 h，源单元积分节点偏移为 $s_{ab}=(h t_a/2,h t_b/2)$，权重 $w_{ab}=h^2w_aw_b/4$。在全部位移上定义参考虚部核

$$
K_{ij}=\frac{k^2}{4}\sum_{ab}w_{ab}J_0(k|x_i-x_j-s_{ab}|).
$$

实际算法在非对角采用该积分，在对角采用正方形自项 $d_0$；奇异性只在实部的 Y0，自项虚部有限。因此

$$
\operatorname{Im}_H D=K+\delta_{\rm self}I,
\quad \delta_{\rm self}=\operatorname{Im}d_0-K_{ii}.
$$

## 命题：当 kh≤π 时，K 半正定

由 J0 的平面波角向表示和对称积分节点，

$$
K_{ij}=\frac{k^2}{8\pi}\int_0^{2\pi}
e^{ik\hat s(\phi)\cdot(x_i-x_j)}\Psi(\phi)\,d\phi,
$$

其中

$$
\Psi(\phi)=\frac{h^2}{4}
\left[\sum_aw_a\cos\left(\frac{kh}{2}t_a\cos\phi\right)\right]
\left[\sum_bw_b\cos\left(\frac{kh}{2}t_b\sin\phi\right)\right].
$$

若 kh≤π，则每个余弦非负，正权重保证 Ψ≥0。故任意 z 满足

$$
z^*Kz=\frac{k^2}{8\pi}\int_0^{2\pi}
\Psi(\phi)\left|\sum_i z_i e^{-ik\hat s(\phi)\cdot x_i}\right|^2d\phi\ge0.
$$

这是充分条件，不是必要条件。条件失败只返回未知，不能说矩阵必然不稳定。

## 推论：可计算 Γ 不需要全维逆矩阵

在非零材料活动单元上，$L=X^{-1}-D$，于是

$$
-\operatorname{Im}_H L\succeq
\Gamma=\operatorname{diag}\left(\frac{\operatorname{Im}\chi_i}{|\chi_i|^2}+\delta_{\rm self}\right).
$$

Γ 正定时可使用 A2 的加权逆界及 primal/adjoint 输出界。自项修正只需要当前实现的一个对角值及 q² 次 J0 计算；材料项 O(N)。与点核的 Gram＋自项证明一致，但这里明确处理了源单元积分所带来的角向权重。

耗散证明限于声明的精确离散核；浮点特殊函数、积分舍入及 FFT 舍入不是区间算术认证。原始模型误差和有限网格误差也不因此消失。已知 χ=0 单元应消元；本轮函数返回不可用而不偷偷加损耗。

## 实现与验证

- `code/a2/certificate.py::gamma_for_fft_kernel`：先调用 `verify_scalar_kernel`，重新生成受信 FFTGreen 核并核对全部卷积核、正向/伴随 FFT 缓冲区及网格元数据，再计算自项修正。结构不同返回不可用；这个重新生成成本显式计时，不能称为免费。结果保存号约定、内积、自项和构造器源代码哈希。此接口防止意外修改，并非针对 Python monkey-patching 的安全边界。
- `code/a2/check_quadrature_certificate.py`：N8、N12；kh=0.2、1.0、2.9；点/积分核共 12 个小矩阵的 Hermitian 特征值检查。最小余量约 −1.04×10⁻¹⁴，处于浮点误差量级；kh=3.5 的积分核正确返回条件不满足。
- `runs/a2/quadrature_certificate/results.json`：全部结果与源代码哈希。
- Review 1 后追加两项回归：修改非对角核、仅修改 FFT 缓冲区均被拒绝；单纯保留自项/元数据不能再取得正面证书。

这补足了一个实现条件，但没有证明材料反演成功、端口压缩有效，或算法比已有认证降阶方法更快。
