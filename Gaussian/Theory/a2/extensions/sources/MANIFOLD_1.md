可以。你现在接触到的这些东西，其实可以整理成一张很清晰的“逆散射技术树”。最关键的一句话是：

\[
\boxed{\text{最终效果}=
\text{物理模型}
+\text{正向离散}
+\text{未知量表示}
+\text{逆问题算法}
+\text{正则化}}
\]

所以 **SOM、MoM、Born、Gaussian 根本不是四个平级算法**。例如一个完全合理的方法可以叫：

\[
\text{Near-field full-wave}
+\text{MoM}
+\text{Gaussian basis}
+\text{SOM}.
\]

它们可以同时存在。

---

# 1. 先把整个领域分成四层

| 层 | 典型方法 | 它到底决定什么 |
|---|---|---|
| **物理近似** | Born、Rytov、Fresnel、Fraunhofer、Full-wave | 你保留多少真实散射物理 |
| **正向求解器** | MoM/VIE、DDA、FEM、FDTD、T-matrix | 给定材料，怎么算场 |
| **未知量表示** | Pixel、Fourier、Wavelet、Spline、Gaussian、Level-set | 介电常数/电流长什么样 |
| **逆问题求解器** | BP、BIM、DBIM、CSI、SOM、Gauss-Newton | 怎么从测量反推出材料 |

比如经典 SOM 通常实际上是从

\[
J=\chi(E^{inc}+G_DJ),
\qquad
E^s=G_SJ
\]

开始，再用 MoM 把它离散成矩阵。也就是说：

> **MoM 是 SOM 的地基之一，而不是 SOM 的竞争者。**

DDA 本身也可以看成某类 volume-integral/MoM 离散：把连续物体换成许多相互耦合的偶极子。citeturn620853search1turn620853search5

---

# 2. 第一大类：Born / Rytov / 远场近似 —— “把问题先变简单”

这是整个逆散射里面最便宜的一派。

## Born approximation

核心假设：

\[
E(\mathbf r)\approx E^{inc}(\mathbf r).
\]

于是

\[
E^s(\mathbf r_r)
\approx
k^2\int
G(\mathbf r_r,\mathbf r)
\chi(\mathbf r)
E^{inc}(\mathbf r)\,d\mathbf r.
\]

一下从

\[
\chi\rightarrow E\rightarrow E^s
\]

这个 nonlinear 问题，变成：

\[
\boxed{E^s=A\chi}.
\]

### 优点

实现复杂度可以说是 **★☆☆☆☆**。

速度极快，FFT、backprojection、Fourier diffraction theorem 都可以直接上。

如果再加 far-field approximation：

\[
G(\mathbf r_r,\mathbf r)
\sim
\frac{e^{ikr_r}}{r_r}
e^{-ik\hat{\mathbf r}_r\cdot\mathbf r},
\]

数据基本变成物体 Fourier spectrum 的采样。

这就是为什么经典 diffraction tomography / Ewald-circle 方法这么漂亮。

### 缺点

最大的敌人不是“物体一定要小”，而是：

\[
\boxed{\text{物体引起的累积场扰动必须足够弱。}}
\]

高 contrast、尺寸较大、内部 resonant、多次散射明显时就会失败。

经典比较发现 Rytov 的有效范围在一些 diffraction-tomography 条件下能明显宽于 Born，但它又带来 phase-unwrapping 问题。citeturn873517search10

### 我给它的评价

弱散射：

> **几乎是最优方案。**

强散射：

> **不要指望算法技巧把错误物理模型救回来。**

---

# 3. Rytov：比 Born 更“尊重 phase”

它不是假设

\[
E\approx E^{inc},
\]

而通常写成

\[
E=E^{inc}e^{\psi},
\]

然后近似 \(\psi\)。

因此它天然更适合描述：

- phase delay；
- smooth object；
- refractive-index variation；
- transmission imaging。

它通常比 Born 对“累积 phase”更宽容，但 phase wrapping 是实际的大麻烦。已有系统比较也发现其有效范围可以比 first-Born 更大。citeturn873517search10

我的粗略评级：

| | Born | Rytov |
|---|---:|---:|
| 实现 | ★ | ★★ |
| 速度 | ★★★★★ | ★★★★★ |
| 高 contrast | ★ | ★★ |
| phase object | ★★ | ★★★★ |
| multiple scattering | ✗ | ✗ |
| 定量材料恢复 | 弱散射很好 | 弱/中等相位物体很好 |

---

# 4. 远场、Fresnel、近场其实是另一条轴

这里特别容易混淆：

\[
\boxed{\text{far field}\neq\text{Born}}
\]

你完全可以做：

> **full-wave object + far-field receiver**

也可以做：

> **Born object + near-field receiver。**

---

## Far field

进一步把 Green function 的距离关系线性化。

最大的好处：

\[
\text{scattering}\rightarrow \text{Fourier-domain sampling}.
\]

所以：

- 数学特别漂亮；
- inversion 特别快；
- calibration 更简单；
- measurement geometry 很容易分析。

缺点则是近场中的一部分高空间频率和 evanescent 成分不存在了。

因此它通常非常适合：

- radar；
- diffraction tomography；
- weak scattering；
- 大范围 sensing；
- target localization。

---

# 5. Fresnel / paraxial：中间派

保留 quadratic phase：

\[
R \approx
z+
\frac{(x-x')^2+(y-y')^2}{2z}.
\]

你可以把它理解为：

\[
\text{Fraunhofer}
\subset
\text{Fresnel}
\subset
\text{exact Green}.
\]

实现还是非常便宜，但保留了更多 depth / propagation 信息。

光学 tomography、SAR、holography 里面非常常见。

---

# 6. Exact near-field

直接：

\[
G(\mathbf r,\mathbf r')
=
\frac{e^{ik|\mathbf r-\mathbf r'|}}
{4\pi|\mathbf r-\mathbf r'|}
\]

甚至完整 dyadic Green tensor。

这时 geometry information 最丰富。

特别靠近 object 时还有：

\[
1/r,\quad1/r^2,\quad1/r^3
\]

不同场项，以及 evanescent spatial spectrum。

理论上 near field 可以获取 far field 没有的高空间频率信息，但这些成分会随距离迅速衰减，所以现实瓶颈转变成：

\[
\boxed{\text{SNR + calibration + geometry accuracy}}
\]

而不是数学意义上“有没有信息”。

这跟你最近研究的 information window 是完全一致的：**越来越近不一定越来越好，因为你同时会把 calibration/model discrepancy 敏感性放大。**

---

# 7. 第二大类：MoM —— 正统 full-wave 暴力派

Method of Moments 的本质非常简单。

假设

\[
J(\mathbf r)
=
\sum_{n=1}^N a_n f_n(\mathbf r).
\]

代入积分方程：

\[
\mathcal L J = E.
\]

再用 testing functions \(w_m\) 投影：

\[
\langle w_m,\mathcal L f_n\rangle a_n
=
\langle w_m,E\rangle.
\]

最终：

\[
\boxed{\mathbf Za=\mathbf v}.
\]

这就是 MoM。

它是计算电磁里面最经典的 integral-equation 方法之一。citeturn620853search7turn620853search14

---

# 8. MoM 为什么强？

它没有 Born approximation。

所以：

\[
J=\chi E
\]

以及

\[
E=E^{inc}+GJ
\]

里面所有 multiple scattering 都在：

\[
(I-\chi G)^{-1}
\]

里面。

因此它可以处理：

- high contrast；
- resonance；
- internal scattering；
- complex phase；
- arbitrary shapes；
- quantitative permittivity。

这是你现在 full-wave 研究真正的 ground truth 之一。

---

# 9. MoM 为什么贵？

普通 dense VIE MoM：

\[
\text{storage}=O(N^2),
\]

直接求解：

\[
O(N^3).
\]

迭代矩阵乘：

\[
O(KN^2).
\]

这是经典 VIE-MoM 最大的问题。使用 FFT/FMM/MLFMA 后才可以把 matrix-vector multiplication 大幅降下来。citeturn620853search13turn620853search16

规则网格情况下：

\[
GJ
\]

具有 convolution / Toeplitz structure，FFT 后：

\[
O(N\log N).
\]

所以你过去看到 FFT-TSOM，并不是“FFT 让 SOM 有魔法”，而是它利用了这个 Green operator 结构。

---

# 10. DDA：MoM 的一个非常物理直观的亲戚

Discrete Dipole Approximation：

\[
\mathbf p_i
=
\alpha_i
\left(
E_i^{inc}
+
\sum_{j\neq i}
G_{ij}\mathbf p_j
\right).
\]

也就是：

> 整个物体 = \(N\) 个互相散射的小 dipole。

因此特别适合：

- particles；
- nanoparticles；
- arbitrary shape；
- optical scattering；
- heterogeneous dielectric。

DDA 可以严格从 volume integral equation 离散推出来，所以它和 VIE/MoM 在数学上关系非常深。citeturn620853search0turn620853search4

这和你提出的 Gaussian 想法其实已经非常接近：

\[
\boxed{\text{DDA: point dipoles}}
\]

versus

\[
\boxed{\text{你的想法: finite-size Gaussian scatterers}}.
\]

这个比较非常重要。

---

# 11. 第三大类：BIM / DBIM —— 每次只相信局部线性

Born Iterative Method：

先有当前模型 \(\chi_k\)，求 full-wave field，再更新：

\[
\chi_{k+1}
=
\chi_k+\Delta\chi.
\]

DBIM 更进一步，在当前背景附近构造 distorted Green function：

\[
G\rightarrow G[\chi_k].
\]

所以每次都解决：

\[
\delta E
\approx
J(\chi_k)\delta\chi.
\]

它实际上就是：

> **不断更新 Born expansion 的展开点。**

因此可以处理 first-Born 已经失败的问题。DBIM 的经典优势就是能把较强 scattering 拆成一串局部线性问题，但每轮需要 forward solution，因此成本明显增加。citeturn873517search3turn873517search9

我的评级：

\[
\boxed{\text{非常靠谱，但非常传统。}}
\]

如果你今天需要做一个 quantitative microwave tomography baseline：

**DBIM 几乎必放。**

---

# 12. CSI：把“材料”和“电流”一起优化

Contrast Source Inversion 定义：

\[
w=\chi E.
\]

于是同时要求：

数据方程：

\[
E^s=G_Sw,
\]

状态方程：

\[
w=\chi(E^{inc}+G_Dw).
\]

然后最小化两种 mismatch。

它最有意思的一点就是：

> 不必每次 outer iteration 都严格完整求解 forward problem。

因此曾经在 nonlinear inverse scattering 里面非常有影响力。

优点：

- full-wave；
- 非线性；
- 不需要 traditional nested forward solve；
- implementation 比完整 Hessian inversion 更友好。

缺点：

- unknown 从 \(\chi\) 扩充成 \((\chi,w)\)；
- conditioning 仍然困难；
- iteration 可能很多；
- high contrast/local minima 仍然存在。

近年来论文仍然明确把 BIM/DBIM/CSI/SOM 当作 full-wave inverse scattering 的主要传统家族来比较。citeturn873517search2turn873517search11

---

# 13. SOM：CSI 的“spectral compression 版亲戚”

这是你非常熟悉的。

从

\[
E^s=G_SJ
\]

做：

\[
G_S=U\Sigma V^H.
\]

于是把 current 写成：

\[
J=J^{det}+J^{amb}.
\]

大 singular value 对应的部分直接从 data 解：

\[
J^{det}
=
\sum_{i=1}^{r}
\frac{u_i^HE^s}{\sigma_i}v_i.
\]

剩余部分：

\[
J^{amb}
=
\sum_{i>r} c_i v_i
\]

才优化。

所以 SOM 的哲学非常漂亮：

\[
\boxed{
\text{data 已经确定的 current 不要再优化；
只优化 data 看不清楚的 current。}
}
\]

这是 Chen SOM 的核心。citeturn546883search0turn331169search4

---

# 14. SOM 真正比 CSI 多了什么？

不是新的 Maxwell physics。

而是：

\[
\boxed{\text{spectral model reduction}}
\]

也就是优化空间 reduction。

所以论文里传统结果常见：

- convergence 更快；
- noise robustness 好；
- unknown current dimension 更低。

但这一点我要结合你自己的研究特别强调：

> **它并不意味着在现代 matched implementation 下，SOM 必然比 Krylov/VarPro/现代 iterative solver 快。**

你自己的 A3 冻结实验已经很好地说明：

\[
\text{dimension reduction}
\not\Rightarrow
\text{wall-clock superiority}.
\]

这也是为什么我一直赞成你不要再把“SOM 更快”作为现在研究的核心卖点。

---

# 15. TSOM / FFT-TSOM / NFFT-SOM / WT-TSOM

SOM 后来发展出一整个 family。

Twofold SOM 再利用内部 current-to-field operator 的谱结构，进一步限制需要优化的 current subspace。传统文献报告这样可以进一步加快 convergence。citeturn620853search6turn620853search11

FFT 系方法则利用 Green convolution：

\[
GJ
\leftrightarrow
\hat G\hat J.
\]

后来 NFFT-SOM 甚至用完整 Fourier basis 替代传统的小 singular-vector basis，主要目标就是降低 SVD overhead、增加抗噪性并简化实现。citeturn546883search31

Wavelet-TSOM 则把 Fourier basis 换成 wavelet，用 multiresolution/sparsity 当 regularization。已有工作报告它能改善 highly nonlinear inverse scattering 的稳定性。citeturn380530search3

所以你可以把整个 SOM 家族理解成：

\[
\boxed{
\text{SOM}
=
\text{current decomposition}
+
\text{basis choice}
+
\text{regularized optimization}.
}
\]

---

# 16. 第四类：真正决定“物体长什么样”的——基函数

这一层对你的 Gaussian 项目尤其重要。

## Pixel / voxel / pulse basis

\[
\chi(\mathbf r)
=
\sum_i \chi_i
\mathbf 1_{\Omega_i}.
\]

这是最普通的方法。

### 优点

万能。

不做任何 shape assumption。

### 缺点

一个 \(100^3\) grid：

\[
N=10^6.
\]

你还什么 inverse problem 都没算，就先有一百万 unknown。

而且结果天然容易：

- noisy；
- staircase；
- checkerboard；
- ringing。

所以需要 TV/Tikhonov/sparsity。

---

# 17. FEM nodal / polynomial basis

比 pixel smooth。

适合 complicated geometry 和 PDE solver。

但 mesh generation 麻烦。

而且 inverse problem 每一步 mesh / sensitivity / adjoint 都比较重。

---

# 18. Fourier basis

\[
\chi(\mathbf r)
=
\sum_{\mathbf k}
c_{\mathbf k}e^{i\mathbf k\cdot\mathbf r}.
\]

和 Born/far-field 是天生一对。

因为 measurement 本身就在 Fourier domain。

优点：

- FFT；
- cheap；
- global smooth structure 很好。

缺点：

sharp boundary 会 Gibbs：

\[
\text{discontinuity}
\Rightarrow
\text{many Fourier coefficients}.
\]

---

# 19. Wavelet

\[
\chi
=
\sum c_{j,k}\psi_{j,k}.
\]

它最大的优势是：

\[
\boxed{\text{空间定位}+\text{尺度定位}}
\]

同时存在。

因此 piecewise smooth objects 往往很 sparse。

这也是为什么 inverse scattering 已经有大量 wavelet regularization；例如 microwave imaging 中 wavelet projection 被用来增加重建稳定性，而 SOM family 也有 WT-TSOM。citeturn380530search1turn380530search7

所以：

- boundary：比 Fourier 好；
- smooth region：比 pixel 高效；
- multiscale：非常自然。

如果你问我传统 basis 谁最均衡：

\[
\boxed{\text{wavelet 很强。}}
\]

---

# 20. Level-set / shape basis

如果已知物体就是几个 homogeneous region：

\[
\chi(\mathbf r)
=
\chi_1
\mathbf1_{\phi(\mathbf r)>0}
+
\chi_0
\mathbf1_{\phi(\mathbf r)<0}.
\]

那根本没必要估计十万个 pixel。

直接估：

\[
\phi(\mathbf r).
\]

适合：

- sharp boundaries；
- piecewise-constant materials；
- target shape reconstruction。

如果真实 object 很符合假设，它会吊打 pixel。

但复杂 continuous material distribution 就不适合。

---

# 21. 现在来到你的 Gaussian

你当前本质上是在考虑：

\[
\chi(\mathbf r)
=
\sum_{m=1}^{K}
a_m
\exp
\left[
-\frac12
(\mathbf r-\mu_m)^T
\Sigma_m^{-1}
(\mathbf r-\mu_m)
\right].
\]

未知量：

\[
\theta_m=
(a_m,\mu_m,\Sigma_m).
\]

甚至还可以让 Gaussian 带：

- anisotropy；
- complex amplitude；
- dispersion；
- orientation；
- polarization response。

这和普通 voxel 有一个根本差异。

假设：

\[
N_{\rm voxel}=10^6,
\]

而你只要：

\[
K=500
\]

个 Gaussian，每个十几个参数：

\[
N_{\rm param}\sim 5000.
\]

那 optimization dimension 可以下降两个数量级。

---

# 22. Gaussian 真正的优点并不是“看起来像 3DGS”

而是数学上的三个性质。

第一：

\[
\boxed{\text{smooth compact-ish representation}}
\]

一个 Gaussian 可以代表一整个 continuous region。

第二：

\[
\boxed{\text{analytic differentiability}}
\]

比如：

\[
\frac{\partial\chi}{\partial \mu},
\quad
\frac{\partial\chi}{\partial\Sigma},
\quad
\frac{\partial\chi}{\partial a}
\]

全都非常漂亮。

第三：

很多 convolution / Green interaction 可以得到 closed-form、semi-analytic 或低成本 approximation。

这才可能是你的真正 computational advantage。

实际上 inverse scattering 文献里很早就有人使用 Gaussian basis 表示介电分布；例如一项 3-D microwave inversion 工作使用 737 个 3-D Gaussian basis，并与稀疏 regularization 配合，获得比普通 \(l_2\) 重建更少的伪影。citeturn380530search0

所以：

\[
\boxed{\text{Gaussian basis 本身并不是新的。}}
\]

你的创新必须出现在后面。

---

# 23. 你的 Gaussian 真正可能新的地方

你正在试图把它变成：

\[
\boxed{
\text{Gaussian = electromagnetic scattering primitive}
}
\]

而不仅仅是

\[
\boxed{
\text{Gaussian = dielectric interpolation basis}.
}
\]

区别非常大。

普通 Gaussian basis：

\[
\chi=\sum_m a_mg_m
\]

然后还是 voxel/MoM/full solver。

你的设想则更接近：

\[
\text{Gaussian}_i
\xleftrightarrow{S_{ij}}
\text{Gaussian}_j
\]

甚至每个 Gaussian 有 local scattering operator：

\[
T_i
\]

然后 aggregate：

\[
T_{\rm total}
=
T_1\oplus T_2\oplus\cdots
+
\text{multiple-scattering coupling}.
\]

于是：

\[
\boxed{
\text{large voxel Maxwell solve}
\rightarrow
\text{Gaussian interaction network}
}
\]

这就开始接近：

- DDA；
- coupled dipole；
- multiple scattering T-matrix；
- domain decomposition；
- reduced-order network。

这比“把介电常数画成 Gaussian”有价值得多。

---

# 24. 你之前提的 S-matrix 想法放在这里就很自然

如果每个 Gaussian primitive 有：

\[
b_i=S_i a_i
\]

而 primitive 之间 propagation 是：

\[
a_i
=
a_i^{inc}
+
\sum_{j\neq i}
P_{ij}b_j,
\]

那么总体：

\[
b
=
S(a^{inc}+Pb)
\]

所以

\[
\boxed{
b=(I-SP)^{-1}Sa^{inc}.
}
\]

看到没有？

它实际上就是：

\[
J=(I-\chi G)^{-1}\chi E^{inc}
\]

的 **block-scatterer reduced version**。

所以如果这个东西推得成功，它不是 SOM 的小改进，而是：

\[
\boxed{\text{新的 full-wave discretization/model reduction architecture}.}
\]

这比单纯 Gaussian inverse representation 的研究价值高得多。

---

# 25. 但 Gaussian 最大的问题也非常明显：边界

假设真实 object：

\[
\chi(x)=
\begin{cases}
3,&x<0\\
0,&x>0
\end{cases}.
\]

Gaussian 是 analytic smooth function。

有限个 Gaussian 和：

\[
\sum_i a_i e^{-(x-\mu_i)^2/2\sigma_i^2}
\]

天然很难精确表达 discontinuity。

所以你的“长毛问题”本质上就是：

\[
\boxed{
\text{用 smooth primitive 拟合 non-smooth boundary}
}
\]

优化器会不断：

- shrink；
- split；
- pile up；
- elongate Gaussian。

最终边缘出现很多小 Gaussian。

这恰好和 3D Gaussian Splatting 的 densification 类似。

---

# 26. 所以不同 basis 的边界能力可以这么排

| Basis | smooth object | sharp edge | sparsity | 参数效率 |
|---|---:|---:|---:|---:|
| voxel | ★★★ | ★★★★★ | ★ | ★ |
| Fourier | ★★★★★ | ★★ | ★★★ | ★★★ |
| Wavelet | ★★★★ | ★★★★★ | ★★★★★ | ★★★★ |
| B-spline | ★★★★★ | ★★★ | ★★★ | ★★★★ |
| Gaussian | ★★★★★ | ★★～★★★ | ★★★★ | ★★★★★ |
| anisotropic Gaussian | ★★★★★ | ★★★★ | ★★★★★ | ★★★★★ |
| level-set | ★★ | ★★★★★ | ★★★★★ | ★★★★★ |

所以你的下一步绝对不应该只用：

\[
\text{isotropic Gaussian}.
\]

而应该至少：

\[
\boxed{
\Sigma_i
=
R_i
\begin{bmatrix}
\sigma_1^2&&\\
&\sigma_2^2&\\
&&\sigma_3^2
\end{bmatrix}
R_i^T.
}
\]

让 Gaussian 沿 boundary 变成长椭球。

这样边界 approximation efficiency 会明显改善。

---

# 27. Gaussian 和 wavelet 的本质区别

这个非常值得你记。

Wavelet：

> “我的 object 在什么**位置和尺度**上有结构？”

Gaussian：

> “我的 object 可以由哪些**连续空间 primitive** 拼起来？”

因此 wavelet coefficient 是：

\[
c_{j,k},
\]

主要还是 coefficient optimization。

Gaussian 则优化：

\[
a,\mu,\Sigma.
\]

因此 topology / geometry 本身就是变量。

这是 Gaussian 最大的潜力。

---

# 28. Gaussian 和 DDA 的区别更有意思

DDA：

\[
\text{固定位置 point dipoles}
\]

Gaussian：

\[
\text{可移动、可伸缩的 finite-volume scatterers}.
\]

可以粗略理解成：

\[
\boxed{
\text{Gaussian primitive}
\approx
\text{adaptive high-order DDA element}.
}
\]

如果最后能得到好的 interaction law：

\[
K_{ij}
=
\iint
g_i(\mathbf r)
G(\mathbf r,\mathbf r')
g_j(\mathbf r')
\,d\mathbf r\,d\mathbf r',
\]

而这个积分有便宜的 semi-analytic approximation，

那么这才是真正值得研究的地方。

---

# 29. 这里顺便把 FEM / FDTD 放进来

它们主要还是 forward solver。

| Forward solver | 优势 | 最大问题 |
|---|---|---|
| MoM/VIE | open-domain，自然保留 radiation | dense matrix |
| DDA | arbitrary dielectric、物理直观 | 大 N、grid error |
| FEM | complicated material/geometry | mesh + absorbing boundary |
| FDTD | broadband 一次计算 | 3-D memory/time 巨大 |
| T-matrix | repeated/simple particles 极快 | complex arbitrary object 不容易 |
| Gaussian reduced model | 潜在 adaptive + low DoF | interaction law 尚未成熟 |

FDTD 每次 forward 可以直接推进 Maxwell：

\[
E^{n+1},H^{n+1/2}.
\]

优点是一次 broadband。

但如果 inverse optimization 需要几百次 forward：

\[
500\times FDTD
\]

就非常痛苦。

所以 adjoint method 很重要。

---

# 30. Gauss-Newton / adjoint optimization

现代 full-wave inversion 其实大量可以统一成：

\[
\min_\chi
\frac12
\|F(\chi)-y\|^2
+
R(\chi).
\]

然后：

\[
\nabla_\chi L
\]

通过 adjoint state 求。

这个框架最大的优势：

> **你想换什么 forward solver、什么 basis、什么 regularizer 都可以。**

所以从现代 numerical inverse problems 的角度：

SOM 并不是唯一的主角。

一个强 implementation 往往是：

\[
\boxed{
\text{full-wave forward}
+
\text{adjoint}
+
\text{L-BFGS/GN}
+
\text{good parameterization}.
}
\]

这也是为什么你现在 Gaussian 更应该和：

- voxel + L-BFGS；
- wavelet + L-BFGS；
- level-set；
- DBIM；
- SOM；

竞争，而不是只和传统 SOM 比。

---

# 31. 那么“实际效果”到底怎么排？

如果是 **weak scattering quantitative imaging**：

\[
\boxed{
\text{Born/Fourier}
>
\text{其他复杂 full-wave 方法}
}
\]

因为后者白白增加复杂度。

如果是 **moderate scattering**：

\[
\boxed{
\text{DBIM / CSI / SOM / GN}
}
\]

开始体现优势。

如果是 **strong multiple scattering / resonant object**：

\[
\boxed{
\text{真正 full-wave}
}
\]

是原则上的必需。

但此时：

\[
\text{physics accuracy}\uparrow
\]

并不等于：

\[
\text{inverse reliability}\uparrow.
\]

因为 nonlinear landscape、calibration error、model discrepancy 都更严重。

---

# 32. 如果单纯按照“工程实现难度”排

我大致给你这么排：

| 方法 | 实现难度 | 算力 | 强散射 | 定量恢复 |
|---|---:|---:|---:|---:|
| Far-field Born | ★ | ★ | ★ | ★★ |
| Born near-field | ★★ | ★★ | ★ | ★★ |
| Rytov | ★★ | ★ | ★★ | ★★★ |
| BP / diffraction tomography | ★ | ★ | ★ | ★★ |
| BIM | ★★★ | ★★★ | ★★★ | ★★★★ |
| DBIM | ★★★★ | ★★★★ | ★★★★ | ★★★★ |
| CSI | ★★★★ | ★★★ | ★★★★ | ★★★★ |
| SOM | ★★★★ | ★★★ | ★★★★ | ★★★★ |
| TSOM / WT-TSOM | ★★★★★ | ★★★ | ★★★★ | ★★★★ |
| Full-wave adjoint GN | ★★★★ | ★★★★ | ★★★★★ | ★★★★★ |
| voxel MoM brute force | ★★★★ | ★★★★★ | ★★★★★ | — |
| DDA | ★★★ | ★★★★ | ★★★★★ | — |
| Gaussian basis inversion | ★★★ | ★★～★★★ | 取决于 forward | ★★★★ |
| 你的 Gaussian scattering network | ★★★★★研究难度 | **目标★★** | **目标★★★★★** | **目标★★★★★** |

最后一行为什么这么奇怪？

因为如果你的理论最终成功：

\[
\text{开发难度很高}
\]

但是 **部署后的计算复杂度应该很低**。

否则它就没有存在意义。

---

# 33. 如果按“物理保真度”排

大致是：

\[
\text{far-field Born}
<
\text{near-field Born}
<
\text{Rytov / extended Born}
<
\text{BIM}
<
\text{DBIM}
\approx
\text{CSI/SOM}
<
\text{exact full-wave optimization}.
\]

注意 SOM 本身并不会改变 Maxwell physics。

如果 SOM 使用同一个 Lippmann–Schwinger full-wave model，那么：

\[
\boxed{
\text{SOM 与普通 full-wave inversion 的 forward physics 一样。}
}
\]

区别只是 parameterization / optimization。

---

# 34. 如果按速度排

弱散射：

\[
\boxed{
\text{Born/Fourier}
\gg
\text{其他全部}
}
\]

这是没有悬念的。

full-wave 里面真正决定速度的往往不是名字叫：

> SOM / CSI / GN

而是：

\[
\boxed{
\text{一次 }Gx\text{ 有多贵}
}
\]

以及

\[
\boxed{
\text{总共做多少次 }Gx.
}
\]

所以：

- FFT；
- FMM；
- low rank；
- domain decomposition；
- reduced-order basis；
- Gaussian interaction operator；

往往比“优化器名字”重要。

---

# 35. 这恰恰解释了你的 Gaussian 为什么值得继续研究

如果 Gaussian 仅仅：

\[
\chi(x)
=
\sum_i a_i g_i(x)
\]

然后最后每次还是把它 rasterize 到：

\[
256^3
\]

voxel，再做一次 full MoM：

那其实价值有限。

你只是把 regularization 换了。

真正有价值的是：

\[
\boxed{
\text{inverse parameterization}
+
\text{forward discretization}
+
\text{multiple-scattering solver}
}
\]

全部都 Gaussian-native。

即：

\[
\text{Gaussians}
\rightarrow
K_{ij}
\rightarrow
\text{multiple scattering network}
\rightarrow
E^s.
\]

这就可能真正改变复杂度。

---

# 36. 最后我给你一张最重要的“选方法地图”

| 你面对的问题 | 我会首先选 |
|---|---|
| 弱散射、速度第一 | Born + Fourier/FFT |
| smooth phase object | Rytov |
| moderate nonlinear quantitative imaging | DBIM |
| 不想每轮严格 forward solve | CSI |
| current operator 有明显可利用 spectrum | SOM |
| multiscale / piecewise smooth | Wavelet + full-wave inversion |
| 已知目标是明确边界区域 | Level-set |
| 任意复杂 object、baseline ground truth | voxel MoM/FEM/FDTD |
| particle scattering | DDA / T-matrix |
| 很多重复 scatterer | T-matrix / S-matrix |
| smooth/blob-like 3D object | Gaussian basis |
| 希望同时降 forward + inverse DoF | **你的 Gaussian scattering primitive 路线** |

所以站在你现在整个科研路线来看，我对这些方法的定位其实已经比较明确了：

\[
\boxed{
\text{Born/Fourier 是速度上限 baseline}
}
\]

\[
\boxed{
\text{MoM/DDA 是物理可信 baseline}
}
\]

\[
\boxed{
\text{DBIM/GN 是现代 full-wave inversion baseline}
}
\]

\[
\boxed{
\text{SOM 是有趣的 current-space reduction}
}
\]

而你现在真正可能做成新东西的，不是再造一个 **“Gaussian-SOM”**，而是把 Gaussian 提升成：

\[
\boxed{
\textbf{Adaptive finite-size electromagnetic scattering primitive}
}
\]

再得到

\[
\boxed{
\textbf{Gaussian--Gaussian full-wave interaction operator}
}
\]

甚至

\[
\boxed{
\textbf{Gaussian scattering network / block T-matrix}
}
\]

最后让复杂度从“细 voxel 上求 Maxwell”转成“少量 adaptive scatterer 之间求 interaction”。

**如果这一层真的成立，它的理论和计算价值明显高于单纯把 Twofold SOM 的 basis 从 Fourier/wavelet 换成 Gaussian。**这也解释了为什么你前面实验里“Gaussian + Twofold SOM”目前没有出现特别夸张的提升：你改善的主要还是 representation；真正占成本、决定 nonlinear physics 的 full-wave interaction 并没有被彻底改写。
