对，你这个怀疑是对的。我查了一圈以后，结论不是“完全没人做过”，而是：

> **Bayesian / probabilistic inverse scattering 已经存在；甚至“每个像素属于若干材料类别，并得到每类的后验概率”也已经有人做到 nonlinear microwave inverse scattering 里了。**
>
> 但是你现在这个 Gaussian 路线还有一个更有意思、也更可能真正原创的版本：**不要让 Gaussian 直接代表一个模糊的介电常数，而让 Gaussian 去参数化一个 probability-simplex field，然后条件于真实材料 realization 解 full-wave Maxwell。**

而且这样一改，恰好能指出我们之前 Gaussian“概率包”推导里一个很可能存在的理论漏洞。

---

# 1. 先回答：有没有人做过“概率语言的逆散射”？

有，而且有几条已经很接近。

最直接的是 Gharsalli、Ayasso、Duchêne、Mohammad-Djafari 这一系列工作。他们把 nonlinear microwave imaging 写成 Bayesian inverse scattering：

- 每个 pixel 有隐藏材料标签 \(z_i\)；
- 目标由有限个已知/待估材料类别组成；
- contrast 使用 Gaussian mixture；
- 空间标签用 Potts / Markov field，让相邻 pixel 倾向属于同一材料；
- 同时估计 induced current、complex contrast、segmentation labels 和 hyperparameters；
- 最后用 Variational Bayes 近似完整 posterior。

也就是说，算法内部实际上已经会产生类似

\[
q_i=
[q_i(1),q_i(2),\ldots,q_i(K)]
\]

这样的**每像素材料概率向量**。citeturn255194view4turn503690search0

而且更早还有 microwave imaging 的 Gibbs/Markov random field 方法；2012 年也已有 optical diffraction tomography 的 nonlinear variational Bayesian formulation。citeturn503690search7turn503690search1

到现在，Bayesian inverse scattering / uncertainty quantification 已经是相当成熟的一条分支；2026 年仍然有专门讨论 inverse scattering UQ 的文章。citeturn133735search4

另一方面还有另一条完全不同的传统：**random-medium scattering**。这里不是“我不知道材料是什么”，而是 permittivity 本身就是随机场，然后研究 stochastic Maxwell、ensemble-averaged field、stochastic Galerkin 等问题。3D stochastic Maxwell 的 forward/inverse 算法也已经有人做。citeturn133735search13turn155100search67

所以：

**不能 claim：**

> “我们首次用 probability 做 inverse scattering。”

也不能简单 claim：

> “首次让每个点对应一个 material probability vector。”

这个已经非常危险。

---

# 2. 但是这里有一个非常关键的概念分叉

你说

\[
p(\mathbf r)
=
[p_1(\mathbf r),\ldots,p_K(\mathbf r)]
\]

是一个概率向量。

这句话实际上有 **3 种完全不同的物理含义**。

| \(p_k(\mathbf r)\) 是什么 | 真正含义 | forward 怎么做 |
|---|---|---|
| epistemic probability | 我“不知道”这里是哪种材料 | 真材料仍然是某一个 \(k\)，对 latent realization 解 Maxwell |
| physical volume fraction | 一个 voxel 真的由 30% A + 70% B 构成 | 应用 homogenization / effective-medium theory |
| stochastic material | 每次 realization 材料随机变化 | stochastic Maxwell / random-medium scattering |

这三个绝对不能混。

而我认为我们之前 Gaussian probability package 最值得怀疑的地方就在这里。

---

# 3. 我认为之前推导最可能错在哪里

假设

\[
z_i\sim \operatorname{Categorical}(p_{i1},...,p_{iK})
\]

而

\[
\chi_i=\chi_{z_i}.
\]

定义平均 contrast

\[
\bar\chi_i
=
\mathbb E[\chi_i]
=
\sum_kp_{ik}\chi_k.
\]

如果之前做的是：

\[
\boxed{
p
\rightarrow
\bar\chi
\rightarrow
F(\bar\chi)
}
\]

然后说这就是 probability model 的 scattered field——

**在 full-wave 条件下通常是不对的。**

因为：

\[
\boxed{
F(\mathbb E[\chi])
\neq
\mathbb E[F(\chi)].
}
\]

Born approximation 下，

\[
F_B(\chi)=A\chi
\]

是线性的，因此

\[
\mathbb E[F_B(\chi)]
=
A\mathbb E[\chi].
\]

这时没问题。

但 full-wave 是

\[
E
=
E^{inc}+G\chi E
\]

或者离散成

\[
E=(I-GX)^{-1}E^{inc}.
\]

显然 nonlinear。

Random-media 理论本身就长期研究这个问题；ensemble-average wave 并不能简单由 average material 代替，higher-order correlations 会进入传播。citeturn155100search68turn155100search3

---

# 4. 更漂亮的是：我们能把这个误差直接推出来

令

\[
X=\operatorname{diag}(\chi)
=M+\Delta X
\]

其中

\[
M=\operatorname{diag}(\mu),
\qquad
\mu=\mathbb E[\chi],
\]

并且

\[
\mathbb E[\Delta X]=0.
\]

定义 mean-material 下的场：

\[
\bar E=(I-GM)^{-1}E^{inc},
\]

再定义

\[
B=(I-GM)^{-1}G.
\]

那么

\[
E
=
(I-GM-G\Delta X)^{-1}E^{inc}.
\]

作 Neumann expansion：

\[
E
=
\bar E
+
B\Delta X\bar E
+
B\Delta X B\Delta X\bar E
+\cdots
\]

取期望：

\[
\mathbb E[E]
=
\bar E
+
B\,
\mathbb E[\Delta X B\Delta X]
\,\bar E
+
O(\Delta\chi^3).
\]

因为 \(\Delta X\) 是 diagonal，

\[
[\mathbb E(\Delta X B\Delta X)]_{ij}
=
C_{ij}B_{ij},
\]

其中

\[
C_{ij}
=
\mathbb E[
(\chi_i-\mu_i)(\chi_j-\mu_j)
].
\]

所以：

\[
\boxed{
\mathbb E[E]
\approx
\bar E+
B(C\odot B)\bar E
}
\]

这里 \(\odot\) 是 Hadamard product。

这个式子非常有意义。

它直接说：

\[
\boxed{
\text{full-wave probability effect}
\sim
\text{material uncertainty correlation}
\times
\text{multiple-scattering propagator}.
}
\]

Born 下这一项消失。

所以你想要的“概率语言”真正获得物理意义的地方，恰好可能就是 **multiple scattering**。

---

# 5. 对 coherent complex material，还有一个容易漏掉的数学细节

因为

\[
\chi
=
\chi'+j\chi''
\]

一般是 complex。

这里事实上要维护两套二阶统计：

Hermitian covariance：

\[
K_{ij}
=
\mathbb E[
\delta\chi_i
\delta\chi_j^*
],
\]

和 pseudo-covariance：

\[
C_{ij}
=
\mathbb E[
\delta\chi_i
\delta\chi_j
].
\]

它们不是一回事。

full-wave mean correction 上面出现的是后者 \(C\)；

而 prediction uncertainty 一阶近似一般是

\[
\Sigma_y
\approx
\Sigma_n
+
J K J^H.
\]

对于 phase-preserving coherent inverse scattering，这个区分其实挺重要。

---

# 6. 那么我建议你的理论体系这样建立

我暂时叫它

\[
\boxed{
\text{Probability-Simplex Full-Wave Inverse Scattering}
}
\]

不是建议最终论文就叫这个名字。

首先定义：

\[
p(\mathbf r)
\in
\Delta^{K-1}
\]

即

\[
p_k(\mathbf r)\ge0,
\qquad
\sum_kp_k(\mathbf r)=1.
\]

latent material field：

\[
z(\mathbf r)
\sim
\operatorname{Cat}(p(\mathbf r)).
\]

材料：

\[
\chi(\mathbf r,\omega)
=
\chi_{z(\mathbf r)}(\omega).
\]

真正的 physics 永远条件于 \(z\)：

\[
E_z
=
E^{inc}
+
G_\omega[
\chi_zE_z
]
\]

以及

\[
y|z
\sim
\mathcal{CN}
\left(
F_\omega[z],
\Sigma_n
\right).
\]

最终 reconstruction 不是一张普通 contrast image，而是：

\[
\boxed{
p_k(\mathbf r|y)
=
P(z(\mathbf r)=k|y)
}
\]

即 K 张 material probability images。

同时给：

\[
\hat\chi(\mathbf r)
=
\sum_k
p_k(\mathbf r|y)\chi_k
\]

以及 uncertainty map：

\[
H(\mathbf r)
=
-\sum_kp_k\log p_k.
\]

因此一套 reconstruction 会直接输出：

**材料图 + 分类概率图 + uncertainty 图 + boundary uncertainty 图。**

---

# 7. Gaussian 在这里反而找到了一个比“Gaussian 是 scatterer”更合理的位置

这是我现在觉得非常有吸引力的一点。

不要简单写：

\[
\chi(\mathbf r)
=
\sum_ja_j
G_j(\mathbf r).
\]

而可以让 Gaussian 去表示 **probability logits**：

\[
\eta_k(\mathbf r)
=
b_k+
\sum_j
a_{jk}
G(
\mathbf r;\mu_j,\Sigma_j
),
\]

然后

\[
\boxed{
p_k(\mathbf r)
=
\frac{
e^{\eta_k(\mathbf r)}
}{
\sum_\ell e^{\eta_\ell(\mathbf r)}
}.
}
\]

于是：

\[
\mathbf p(\mathbf r)
\in\Delta^{K-1}
\]

天然成立。

---

这样一个 Gaussian 不再宣称：

> “我是这一坨物理介电常数。”

而是：

> **“我是空间材料 posterior / belief field 的一个连续 basis function。”**

这在理论上要干净得多。

而且 Gaussian overlap 也不再那么尴尬——它们是在 logit space 叠加。

---

# 8. 这样你的 Gaussian splitting 突然也有了非常自然的解释

以前最大的疑问之一是：

> Gaussian 为什么应该 split？

现在可以定义 posterior entropy：

\[
H(\mathbf r)
=
-\sum_kp_k(\mathbf r)\log p_k(\mathbf r).
\]

如果一个 Gaussian 覆盖区域中同时存在：

- 高 entropy；
- 高 electromagnetic sensitivity；
- probability field 有明显空间变化；

它就应该 split。

比如定义：

\[
S_j
=
\int
G_j(\mathbf r)
H(\mathbf r)
\,
\|J(\mathbf r)\|^2
\,d\mathbf r.
\]

高 \(S_j\)：

\[
G_j
\rightarrow
G_{j1}+G_{j2}.
\]

而：

- posterior 极其确定；
- 两个 Gaussian 给出相同 material distribution；
- spatially redundant；

则 merge。

这就形成：

\[
\boxed{
\text{physics-aware probabilistic Gaussian splitting}.
}
\]

相比“因为 reconstruction error 高所以 split”，理论故事明显强得多。

---

# 9. 我甚至认为可以先推三个相当漂亮的定理

## 定理 1：Born probability collapse

Born 下：

\[
y=A\chi+n.
\]

因此：

\[
\mathbb E[y]
=
A\mu_\chi.
\]

所以：

\[
\boxed{
p(\mathbf r)
\text{ 的一阶概率信息完全 collapse 成 }
\mu_\chi(\mathbf r).
}
\]

所有产生相同

\[
\mu_\chi
=
\sum_kp_k\chi_k
\]

的概率分布，在 mean coherent Born data 下完全不可区分。

这是一个很漂亮的 negative theorem：

> **Born scattering 本身并不能因为你换成 probability representation 就凭空获得 uncertainty information。**

概率表示只有通过 prior、多频 material signatures、data covariance 或 nonlinear scattering 才真正获得额外可辨识内容。

---

# 10. 定理 2：multi-frequency material-simplex identifiability

令第 \(k\) 种材料的频率 signature 为

\[
s_k
=
[
\chi_k(\omega_1),
\ldots,
\chi_k(\omega_F)
]^T.
\]

则

\[
\bar s
=
\sum_kp_ks_k.
\]

什么时候能够由 \(\bar s\) 唯一确定 \(p\)？

定义实矩阵

\[
\mathcal S
=
\begin{bmatrix}
1&1&\cdots&1\\
\Re s_1&\Re s_2&\cdots&\Re s_K\\
\Im s_1&\Im s_2&\cdots&\Im s_K
\end{bmatrix}.
\]

那么 simplex weight 唯一可辨识的条件就是

\[
\boxed{
\operatorname{rank}_{\mathbb R}\mathcal S=K.
}
\]

换句话说：

> 多频并不是单纯“更多 measurement”。

它可以把不同 material 的 dispersion signature 拉开，从而让**概率向量本身**变得 identifiable。

这个我觉得非常值得继续推。

---

# 11. 定理 3：Full-wave probability lift

一般 nonlinear forward map：

\[
F(\chi).
\]

令

\[
\chi=\mu+\delta\chi.
\]

Taylor 展开：

\[
F(\chi)
=
F(\mu)
+
DF(\mu)\delta\chi
+
\frac12
D^2F(\mu)
[
\delta\chi,\delta\chi
]
+\cdots
\]

于是

\[
\boxed{
\mathbb E[F(\chi)]
=
F(\mu)
+
\frac12
D^2F(\mu):C
+
O(\|\delta\chi\|^3).
}
\]

因此：

- Born：
  \[
  D^2F=0
  \]
  probability higher moments invisible；

- full-wave：
  \[
  D^2F\neq0
  \]
  variance / spatial correlation 会进入 coherent prediction。

可以把这个叫做一种：

\[
\boxed{
\text{nonlinear probability lift}
}
\]

当然名字先别急着包装成 contribution；random-media 文献中肯定有近亲，要进一步做 prior-art audit。citeturn155100search2turn155100search68

---

# 12. 基于这套东西，可以 propose 三档算法

## A. 最简单：Born Probability Tomography

计算

\[
\mu_i
=
\sum_kp_{ik}\chi_k
\]

以及

\[
K_\chi
=
\operatorname{Cov}(\chi).
\]

则：

\[
\mu_y=A\mu,
\]

\[
\Sigma_y
=
\Sigma_n
+
A K_\chi A^H.
\]

然后优化 marginal likelihood：

\[
\boxed{
\mathcal L(p)
=
(y-\mu_y)^H
\Sigma_y^{-1}
(y-\mu_y)
+
\log\det\Sigma_y
+
R(p).
}
\]

其中

\[
p_i\in\Delta^{K-1}.
\]

它已经不只是“把概率平均成介电常数再反演”了，因为 covariance 也进入 likelihood。

可以作为非常干净的 baseline。

---

# 13. B. 我最感兴趣：Moment-Closed Full-Wave Probability Inversion

保持：

\[
\mu_\chi(p),\quad
K_\chi(p),\quad
C_\chi(p).
\]

用 full-wave solver 算

\[
F(\mu_\chi)
\]

以及 Jacobian

\[
J=DF(\mu_\chi).
\]

然后：

\[
\mu_y
\approx
F(\mu_\chi)
+
\frac12D^2F:C_\chi,
\]

\[
\Sigma_y
\approx
\Sigma_n
+
J K_\chi J^H.
\]

求：

\[
\min_{p\in\Delta}
(y-\mu_y)^H
\Sigma_y^{-1}
(y-\mu_y)
+
\log\det\Sigma_y
+
R_{\mathrm spatial}(p).
\]

这可以叫：

> **moment-aware full-wave simplex inversion**

它比 full Monte-Carlo Bayesian solver 便宜很多，又没有犯

\[
F(E[\chi])=E[F(\chi)]
\]

那个错误。

---

# 14. C. 最严格：Full Variational Bayesian Full-Wave Inversion

定义

\[
q(z)
=
\prod_i
\operatorname{Cat}
(z_i;p_i).
\]

优化 ELBO：

\[
\boxed{
\mathcal L_{\rm ELBO}
=
\mathbb E_q[
\log p(y|z)
]
-
D_{\mathrm{KL}}
(q(z)\|p_0(z)).
}
\]

其中每个 sample \(z\) 都真正跑：

\[
F(z)
\]

的 full-wave solver。

这个是 mathematically 最干净的。

但代价最大。

可以：

- Gumbel-softmax；
- stochastic VI；
- adjoint gradient；
- local sampling；
- structured Potts VI。

但这条本身和 2014 那套 Bayesian microwave imaging 太接近，因此单独拿出来原创度不够。citeturn255194view4

---

# 15. 真正属于你 Gaussian 路线的版本，我会这样组合

我认为最有希望的算法不是纯 A/B/C，而是：

\[
\boxed{
\text{Adaptive Gaussian Simplex Full-Wave VI}
}
\]

流程：

### Step 1

用少量 Gaussian 表示

\[
\eta_k(\mathbf r)
\]

而不是直接表示 \(\chi\)。

### Step 2

softmax：

\[
\eta\rightarrow p.
\]

### Step 3

从 \(p\) 得到

\[
\mu_\chi,\ K_\chi,\ C_\chi.
\]

### Step 4

使用 moment-closed full-wave likelihood。

### Step 5

更新：

\[
\{\mu_j,\Sigma_j,a_{jk}\}.
\]

### Step 6

根据

\[
\text{entropy}
\times
\text{EM sensitivity}
\times
\text{model curvature}
\]

决定 Gaussian splitting。

### Step 7

只在高 uncertainty 的局部区域运行 expensive stochastic/full Bayesian sampling。

也就是：

\[
\boxed{
\text{deterministic full-wave}
+
\text{probabilistic representation}
+
\text{local Bayesian refinement}.
}
\]

---

# 16. 这还非常可能改善你一直头疼的 boundary 问题

因为传统 Gaussian fitting 的 boundary 很麻烦：

Gaussian 天生 smooth。

但是 probability representation 中，boundary 不需要 Gaussian 自己具有 sharp edge。

例如两材料：

\[
p_A(\mathbf r)
=
\sigma(
\eta_A-\eta_B
).
\]

真正 boundary 是：

\[
\boxed{
\eta_A(\mathbf r)
=
\eta_B(\mathbf r)
}
\]

也就是一个 level set。

所以：

- Gaussian 可以很 smooth；
- **decision boundary 可以很 sharp。**

最终 MAP：

\[
\hat z(\mathbf r)
=
\arg\max_kp_k(\mathbf r)
\]

自然得到 sharp material boundary。

这其实可能比直接用 Gaussian 的 amplitude 拟合方块、尖角、材料界面合理得多。

---

# 17. 但是有一个非常重要的科研判断

目前我不会把它说成：

> “一种从未有人想到过的 probability inverse scattering theory。”

因为这个说法经不起文献。

最直接的先例已经非常明显：

**Gauss–Markov–Potts + nonlinear full-wave Bayesian microwave inversion + pixel material posterior**，十多年前就有了。citeturn133735search12turn503690search0

真正可能有原创价值的是这三个东西合起来：

\[
\boxed{
\begin{aligned}
&1.\quad
\text{Born probability-collapse / full-wave probability-lift theory}\\
&2.\quad
\text{moment-aware simplex full-wave inversion}\\
&3.\quad
\text{Gaussian-logit adaptive probability representation}
\end{aligned}
}
\]

第三项尤其能重新解释你目前 Gaussian 的作用：

> **Gaussian 不一定应该是“物理 scatterer”。**
>
> 它可能更适合成为**材料 posterior probability field 的自适应连续坐标系**。

这个定位我现在反而觉得比“用 Gaussian 把 full-wave scatterer 本身网络化”更稳。

---

## 我对这个方向目前的评级

**单纯 probability inverse scattering：原创度 2/10。**

**Gaussian + Bayesian segmentation：4/10，已有很多近亲。**

但如果你能真正证明：

\[
\text{Born probability collapse}
\rightarrow
\text{full-wave covariance lift}
\]

并据此设计一个**只在 nonlinear interaction 真正提供额外概率信息时才开启 full-wave refinement**的 Gaussian simplex algorithm，我会给这个研究方向大约 **7.5–8/10 的继续研究价值**。

尤其是它非常适合做一个决定性实验：

同样 mean contrast、不同 probability/correlation structure，构造两个在 Born mean data 下不可区分的对象；然后验证 full-wave 二阶/多重散射项什么时候开始把它们分开。如果这个现象足够稳定、还能形成可计算的 sufficient condition，那就从“Bayesian imaging trick”升级成真正的 **electromagnetic probability-observability theory**。
