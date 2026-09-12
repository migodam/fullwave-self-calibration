# 最近邻文献：公式、假设与访问状态

检索截止本轮2026-09-11。以下区分“全文可访问且关键节已核读”“仅部分索引全文”“仅摘要”。来源是原论文、作者上传稿或官方出版机构。没有用未读论文的标题推断其定理。各条只概述与本轮最相关的公式，不是穷尽性综述。

## R1. Marx–Mulholland，1983

*Size and Refractive Index Determination of Single Polystyrene Spheres*, J. Res. NBS 88(5), 321–338. DOI: 10.6028/jres.088.016。NIST官方全文：https://nvlpubs.nist.gov/nistpubs/jres/088/jresv88n5p321_A1b.pdf 。

已核读§3.1–3.3：均匀球、已知平面波、双极化角度强度数据，Mie拟合半径及折射率。式(12)为 Q=Σ(Eᵢ−aTᵢ)²，式(13)为 a=ΣEᵢTᵢ/ΣTᵢ²，a处理未知绝对强度。故“球体Mie材料拟合+标量增益消元”明确不是新颖点。与本轮差异是强度比例而非复增益环带、无本轮路径参考误差预算或同时有损双球有限覆盖。

## R2. Francis–Wittmann，2003

*Uncertainty Analysis for Spherical Near-Field Measurements*, AMTA, 43–45。NIST官方全文：https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=31432 。

已核读式(2)、(9)及表1：tₙₘ=Mₙ⁻¹Tₙₘ；Δt=M⁻¹ΔT+Δ(M⁻¹)T，区分测量场与探头误差，采用 m=±1 对称探头。误差表已有模式纯度、角/径向位置、漂移、泄漏、房间散射等。不能把“误差归因”或“模态探头校正”本身当新贡献；该文也没有为本轮装置验证0.1%的容差。

## R3. Epp–Janz，2013预印本

*Spectral approach to the inverse problem for the field of arbitrary changing electric dipole*. arXiv:1308.1662，全文：https://arxiv.org/pdf/1308.1662 。

已核读§3及误差分析：可测磁场相对空间导数 Δₗₘ=(∂ₗHₘ)/Hₘ；式(14)给 I=(kr)²/[1+(kr)²]，式(15)为 r=k⁻¹√[I/(1−I)]；绝对源幅度另需强度。式(38)显示相位梯度距离误差随距离增长。故用近场结构消去源幅度并估距早已存在。本轮G使用三激励张量而非实测梯度，增加有限竞争世界覆盖；不能据此直接声称思想首创。

## R4. Li–Lee–Bresler，2016

*Optimal Sample Complexity for Blind Gain and Phase Calibration*, IEEE TSP 64(21), 5549–5556。作者全文：https://arxiv.org/pdf/1512.07293 。

已核读模型、定理2.1及条件：Y=diag(λ)AX；已知A，X在给定子空间。n>m且(n−1)/(n−m)≤N≤m时，对几乎所有A、λ、X可辨到一个尺度。其generic双线性结论不能替代固定损耗的非线性Maxwell流形、最坏有限噪声裕量或环带约束。公共尺度歧义、共享快照可校准性不是本轮首创。

## R5. Prokopiou–Tsitsas，2018：直接相关，已取得作者全文

*Electromagnetic Excitation of a Spherical Medium by an Arbitrary Dipole and Related Inverse Problems*, Studies in Applied Mathematics 140(4), 438–464. DOI:10.1111/sapm.12206。出版页：https://onlinelibrary.wiley.com/doi/10.1111/sapm.12206 。作者上传的2017稿全文：https://www.researchgate.net/publication/323615456_Electromagnetic_Excitation_of_a_Spherical_Medium_by_an_Arbitrary_Dipole_and_Related_Inverse_Problems 。

已核读§2、§5–6。单球中心/半径、背景和外偶极已知，远场低频系数按已知A归一化。式(5.14) χ₀=3(εᵣ−1)/[2(εᵣ+2)]；式(6.1) εᵣ=(1+2ψ₁)/(1−ψ₁)，ψ₁=2b³m̃₁/(3a₁³)。源参数反演则在材料已知时进行，且须分离不同频率阶。单球材料/源分类、显式材料逆式、不同阶信息早有先例。本轮不能把它们重新命名为新理论；差异只可落在未知电子学、有限误差及实现预算。

## R6. Tsitsas，2013：仍未完成全文公式核验

*A Low-Frequency Electromagnetic Near-Field Inverse Problem for a Spherical Scatterer*, J. Computational Mathematics 31(5),439–448. DOI:10.4208/jcm.1304-m4388。出版页：https://www.global-sci.com/jcm/article/view/12126 。

摘要及作者稿开头可访问：球内偶极、低频近似、在偶极位置使用次级电场反演复介电常数。出版PDF及作者全文抓取失败，关键逆式与误差假设尚未核对，不能断言本轮条件未被涵盖。网页迁移显示2018发布元数据，但卷期/正文为2013，本文采用2013。

## R7. Vasileiou–Koutsoupidou–Kosmas，2026：部分全文索引

*Impact of Increasing Antenna Model Complexity on Microwave Tomography Using DBIM*, Sensors26(11),3517. DOI:10.3390/s26113517。出版页：https://www.mdpi.com/1424-8220/26/11/3517 。

索引全文给出式(26) Sₚq,sc=(Êₜ/Êc)Sₚq,c 的背景校准；CST数据/FDTD反演以及天线模型差异十分接近本轮失配主题。实际页面/PDF抓取失败，未完成整篇核验。不能将“校准有时有效但不能消除模型误差”当新发现，也不能声称该文已有本轮有限覆盖证书。

## R8. Takenaka–Moriyama，2012：摘要

*Inverse scattering approach based on the field equivalence principle: inversion without a priori information on incident fields*, Optics Letters37(16),3432–3434. DOI:10.1364/OL.37.003432。PubMed：https://pubmed.ncbi.nlm.nih.gov/23381281/ 。摘要的入射场消除依赖边界总E/H场与等效原理，不等于任意复增益自校准。全文公式未取得，不作更细原创性排除。

## R9. Bucci–Franceschetti，1987及相关可恢复信息文献：待补全文

*On the spatial bandwidth of scattered fields*, IEEE TAP35(12),1445–1455，DOI:10.1109/TAP.1987.1144024；Bucci–Isernia1997，*Electromagnetic inverse scattering: Retrievable information and measurement strategies*，DOI:10.1029/97RS01826。本轮只核对到元数据/作者稿部分索引，未完成关键定理全文对照。不能把一般“可恢复信息随采集设计变化”说成首次提出。

## 原创性裁决

最接近的已读公式已经覆盖标量消元、近场结构估距、球材料显式反演和多阶信息分类。可能剩余的研究点是：共同/相对电子学约束下的有限任务窗口及其可验证实现界。但本轮没有完成相关文献的排他性核查，也没有完成C类闭环；不宜在摘要中写“首次”“解决了”或宣称优于既有方法。
