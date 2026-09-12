# 最近邻文献审计：已读公式、适用域与未取得全文

审计日期：2026-09-11。文献不是仅凭标题判定“不一样”。下列“全文已读”指本轮取得可解析正文并核对关键公式；页面截图服务对部分 PDF 返回 InternalError，未据不可见图作结论。其余关键来源明确标为全文不可得，不能据摘要排除原创性冲突。下面讨论与 MANUSCRIPT_ZH.md 的新推导分开。

## 1. 直接重叠的电子学—材料自校准

Ludvig-Osipov, Stenmark, Rylander, McKelvey, “Auto-Calibration of Near-Field Microwave Measurements for Complex Permittivity Estimation,” Progress In Electromagnetics Research M 136, 46–56 (2025)，DOI 10.2528/PIERM25090302。全文：[出版社PDF](https://www.jpier.org/ac_api/download.php?id=25090302)。

已核对：式(4) D=RST，R/T 为收发端对角复增益；式(7) Dtilde=D+RSN^(+)+N^(−)，区分放大前后噪声。假设端口匹配、串扰可忽略、增益在采集中准静态，并已有覆盖均匀材料域的已校准 S(ε) 库。§2.4.1 指出尺度规范并固定 r1=1；§2.4.2 以 S(ε1)=Rtilde S(ε2)Ttilde 仅能在 ε1=ε2 成立定义材料可辨识。式(9)–(10) 是外层材料搜索与内层增益拟合；式(19) 使用 Möbius 近似 M(ε)=(αε+β)/(γε+1)。

因此，联合估计材料和电子学、对增益 profiling、通过相对通道结构排除增益歧义，都不能作为 Q5 未经限定的首创点。本轮不同目标是指定误差集合下的有限材料容差和参考精度边界；该文没有替本轮证明这些边界，也没有替类 C 验证未知接收位移和双球离散误差。该文的事先 S 库不能与一次在线标量参考混成同一成本。

## 2. 单球材料的谱相位和有限先验范围

Romanov, Maltsev, Yurkin, “Retrieving refractive index of single spheres using the phase spectrum of light-scattering pattern,” Optics & Laser Technology 161, 109141 (2023)，DOI 10.1016/j.optlastec.2023.109141。已读作者预印本：[arXiv:2210.00334](https://arxiv.org/pdf/2210.00334)。

决定性式(7)为 F(v)≈exp(i v c0/√m) F_RGD(v/√m)，把谱峰位置与谱相位联系到折射率。§2.3 将 d∈[44,55]、m∈[1.15,1.22] 固定为狭窄工作域，建立位置—相位到尺寸—折射率的插值；相位周期性是限制范围的明确原因。这里的“相位”是散射强度图样的 Fourier 谱相位，不是未知电子学共同相位。RGD 下的归一化图样缺乏该折射率依赖，是文中已经使用的比较。

因此，不得宣称“利用超出弱散射模型的相对谱结构提取材料”“在有限先验范围避开相位歧义”本身为新观点。本轮有限误差参考条件与其观测实验不同，但这种不同尚不足以自动成为 TAP 原创贡献。

## 3. 必须继续取得全文的最近邻

| 文献 | 当前实际访问 | 可确认范围 | 不能据此声称的结论 |
|---|---|---|---|
| Ludlow & Everitt, “Inverse Mie problem,” JOSA A 17(12), 2229–2235 (2000), DOI 10.1364/JOSAA.17.002229 | 出版社全文不可访问；摘要及局部缓存 | 从散射角函数提取 Mie 系数、识别均匀球参数 | 尚未核对其反演公式与稳定性假设，不能排除与受限球体材料识别的重叠 |
| Tsitsas, “A Low-Frequency Electromagnetic Near-Field Inverse Problem for a Spherical Scatterer,” J. Comput. Math. 31(5), 439–448 (2013), DOI 10.4208/jcm.1304-m4388 | 期刊/作者页面摘要可得；完整公式未可靠取得 | 内部电偶极激励、源处近场、低频球介质反演 | 不能以摘要宣称没有有限稳定性或参考条件；网页2018迁移时间不是论文年份 |
| Hislop, Craeye, Gonzalez-Ovejero, TAP 64(4), 1364–1372 (2016), DOI 10.1109/TAP.2016.2526087 | 作者书目及摘要；全文未得 | 以近场已知对象测量校准天线电流模型 | 公共标量 g 不能覆盖此类全部天线模型误差；具体公式比较未完成 |
| Hasar & Simsek, J. Phys. D 42(7), 075403 (2009), DOI 10.1088/0022-3727/42/7/075403 | 大学论文记录及摘要；全文未得 | 与位置/校准不敏感的介电参数测量 | 不把波导/同轴配置与自由空间双球配置视为同一实验，也不凭摘要断言独创 |

全文未取得不等于文献没有相关结论。这四项应保留在 novelty gate，不用已读两篇替代全部最近邻审计。一般集合成员估计、区间法、Chebyshev 中点、经典 Mie/T-matrix 和 GLS 均不是本轮原创。可争取的贡献只能是具体电磁实验域下、经过证明和真实读出验证的有限边界及设计规律。

## 4. 与 Chen 2018 的关系

用户上传的 Xudong Chen《Computational Methods for Electromagnetic Inverse Scattering》用于物理背景：§2.8.3 的 Foldy–Lax 方程区分对象间耦合；§2.9.1 的 CDM/DDA 是独立离散前向模型；§3.2 的模态可读性依赖照射和接收极化。上述背景不为本轮提供新的参考精度结论。本轮没有引用该书作为“全波必定比 Born 更容易自校准”的依据。
