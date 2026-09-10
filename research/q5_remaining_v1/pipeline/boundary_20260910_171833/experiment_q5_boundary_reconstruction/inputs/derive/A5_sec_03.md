## 3. 命题 1：固定损耗改变 Born 增益歧义

设几何已知，$B:\mathbb C^2\to\mathbb C^m$ 满复列秩，$u\in\mathbb R^2$ 是两区实对比度，$\ell\in\mathbb R^2$ 精确已知。考虑

$$y=gB(u+i\ell),\qquad g\ne0.$$

若 $D=\det[u,\ell]\ne0$，则 $(u,g)$ 由无噪数据唯一决定。若 $\ell=0$，允许域内的实尺度歧义仍存在。若 $\ell\ne0$ 且 $u=q\ell$，则任何允许的 $q'$ 给出相同数据，只须补偿增益。

### 证明

两世界相等意味着 $gB(u+i\ell)=g'B(u'+i\ell)$。由 $B$ 的单射性，

$$u'+i\ell=c(u+i\ell),\qquad c=g/g'=a+ib.$$

比较虚部得

$$b u+(a-1)\ell=0.$$

当 $u,\ell$ 线性独立，必有 $b=0,a=1$，再比较实部可得 $u'=u$ 和 $g'=g$。当 $\ell=0$，正实 $c$ 可产生允许的尺度歧义。当 $u=q\ell$，选

$$c=\frac{q'+i}{q+i},\qquad g'=g/c,$$

即可保持同一个已知 $\ell$ 并改变实材料。是否构成反例还要检查材料与补偿增益都在允许域内。证毕。

本轮用 $\ell=(0.03,0.05)$、$q=30,q'=45$ 构造了合法两区域反例，两个实介电向量为 $(1.9,2.5)$ 与 $(2.35,3.25)$，增益约为 $0.666831+0.007404i$，相对残差 $1.74\times10^{-16}$。

### 显式有限逆映射

令 $w=B^\dagger y=g(u+i\ell)$，$h=1/g=a+ib$，则

$$M(w)\begin{bmatrix}a\\b\end{bmatrix}=\ell,\qquad M(w)=[\Im w,\Re w].$$

由于 $\det M(w)=-|g|^2D\ne0$，先解这个实 $2\times2$ 系统，再取 $u=\Re(hw)$。这不是仅仅一个切线秩判据。

设压缩数据误差 $\|\delta w\|\le\eta$，$\mu=\sigma_{\min}(M(w))$ 且 $\eta<\mu$。因为 $\|\delta M\|_2\le\eta$，逆矩阵扰动给出

$$|\delta h|\le\frac{|h|\eta}{\mu-\eta},$$
$$\|\delta u\|\le |h|\eta\left[1+\frac{\|w\|+\eta}{\mu-\eta}\right].$$

若原始误差界为 $\|\delta y\|\le\beta$，可用 $\eta=\|B^\dagger\|\beta$。若损耗也有误差 $\delta\ell$，第一式的分子变为 $|h|\eta+\|\delta\ell\|$。因此“已知损耗”是一项可消耗的参考精度，而不是免费信息。

### 精确局部下界

对 $B=I$、$g=1$，投影掉复增益后的实材料 Gram 矩阵为

$$H=I-\frac{uu^T+\ell\ell^T}{S},\qquad S=\|u\|^2+\|\ell\|^2.$$

证明：对任意实 $v$，复直线投影留下的平方范数为 $\|v\|^2-| (u-i\ell)^Tv|^2/S$，展开即得。于是

$$\lambda_{\min}(H)=\frac{1-\sqrt{1-4D^2/S^2}}2.$$

对于一般 $B$，最小化复增益增量前使用 $\|Bz\|\ge\sigma_{\min}(B)\|z\|$，得到

$$\sigma_{\min}(A_{\rm vis})\ge |g|\sigma_{\min}(B)\sqrt{\lambda_{\min}(H)}.$$

需要白化时再乘 $\sqrt2/\sigma$。此处几何已知；增加几何 nuisance 可能继续损失信息，不能把本下界直接搬过去。$D/S$ 很小时，有正秩但稳定性极弱。

