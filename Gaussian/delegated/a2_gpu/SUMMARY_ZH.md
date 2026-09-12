# RTX 4060 Gaussian/voxel campaign 候选结果段落

完整性检查通过：30 cases × 4 methods，N128 inverse、N192 cell-integrated data、120 steps、统一 config hash `34aff2c1c886`。

该有界 noiseless development comparison 将普通正 Gaussian component 表示（p=96/384/864）与 voxel（p=16384）比较；它不是 SOM 结果，也不证明正则化或步数设置已达到公平充分。每方法秒数来自包含原子 checkpoint I/O 的优化循环，数据生成与 setup 单列。

报告 material L2、held clean-field error、IoU、RHS/迭代计数及 Gaussian K144 相对 voxel 的 10-case family-paired 差值；训练轨迹只是 fit loss，未保存中途 held 轨迹，故不作 held 前沿结论。全局峰值显存不作为方法间比较。

机器可读汇总：`runs/a2/gpu/representation_campaign/main128/SUMMARY.json`。
