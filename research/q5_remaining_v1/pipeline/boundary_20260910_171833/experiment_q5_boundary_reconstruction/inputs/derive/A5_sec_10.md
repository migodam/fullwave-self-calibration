## 10. 命题 7：只有固定 nuisance 秩时才有小扰动保证

令实白化 nuisance 矩阵 $N$ 满列秩，$\sigma_{\min}(N)=\nu>0$；$\widehat N=N+E_N$，$\|E_N\|\le\eta_N<\nu$。则两者同秩，且

$$\|P_{\widehat N}-P_N\|\le\frac{\eta_N}{\nu-\eta_N}.$$

证明：$Q_N\widehat N=Q_NE_N$，而 $\|\widehat N^\dagger\|\le1/(\nu-\eta_N)$，从而 $\|Q_NP_{\widehat N}\|$ 有上述界；两个同秩正交投影的最大主角正弦等于投影差范数。于是对 $\widehat A=A+E_A$，

$$\sigma_{\min}(Q_{\widehat N}\widehat A)
\ge\sigma_{\min}(Q_NA)-\eta_A-\frac{\eta_N}{\nu-\eta_N}\|A\|.$$

这是经典扰动工具，不是新电磁定理。它需要真实导数误差界，而非仅场值网格差。反例 $N_\delta=[e_1,\delta e_2]$、$A=e_2$ 表明：$\delta=0$ 时材料可见，任意非零 $\delta$ 时被 nuisance 完全吸收；场/矩阵扰动虽小，投影仍突变。新增未知模式误差不能不加约束地塞进 nuisance 后仍沿用原谱下界。

