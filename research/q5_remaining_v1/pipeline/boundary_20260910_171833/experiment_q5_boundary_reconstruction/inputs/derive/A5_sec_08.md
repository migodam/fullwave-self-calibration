## 8. 精确 Mie 反应系数与区间导数

令 $m=\sqrt\epsilon$、$x=ka$、$j=j_1$、$y=y_1$，定义

$$D_j(z)=\frac{j_1(z)+zj_1'(z)}z,\qquad
D_y(z)=\frac{y_1(z)+zy_1'(z)}z.$$

电偶极反应系数为

$$q(\epsilon,x)=\frac{m j_1(mx)D_j(x)-j_1(x)D_j(mx)}{m j_1(mx)D_y(x)-y_1(x)D_j(mx)}.$$

本表达式由切向 E/H 连续直接消去内场振幅获得。记分子分母为 $N,D$，则

$$q'=\frac{N'D-ND'}{D^2},\quad m'=\frac1{2m},\quad (mx)'=\frac{x}{2m}.$$

使用

$$j_1(z)=\frac{\sin z}{z^2}-\frac{\cos z}{z},\qquad
D_j(z)=\frac{\sin z}{z}+\frac{\cos z}{z^2}-\frac{\sin z}{z^3},$$
$$D_j'(z)=-j_1'(z)/z-j_1(z)+j_1(z)/z^2,$$

逐项求导即可得到代码中的区间公式。最后

$$h'(\epsilon)=\frac{q'(\epsilon)}{(1+q^2)^{3/2}}.$$

需要控制分母不为零、$q>0$ 及 $q'>0$。所选两个区间的分母下界分别大于 28.14 和 57.13。我们没有声称所有材料、共振或含损耗情形都满足该分支条件。

