# Final Verdict

**NO — known-scatterer/reference calibration remains the better solution**

限定：本轮机械扫描双接收头、两已知支撑材料区域、共享低维增益/延迟及所测误差域；不是否定所有动态阵列自校准。

Engineering win A 未达到。冻结旧参考在变化状态下会失败，但周期参考、每帧参考及其与普通联合全波反演组合没有被候选超过。实际参考停机和状态创新未测得，不能宣称参考成本使它们不可用。

Physics win B 未达到。Born有准确尺度规范；全波固定几何的相对响应有变化；耦合二点完整响应有一个准确代数解。然而当前观测下未知几何、模型误差和有限读出没有构成全域分离保证。真实已知gain/pose反而暴露材料偏差，不能把11/12目标only开发结果当作已证物理机制。没有从别的Q5模态类搬运成功。

“Exciting”四句审计：
- Conventional calibration fails because **旧表无法预测标定后创新**：在声明类内成立；周期参考经济失效未建立。
- Full-wave measurement contains **相对模态/相互作用结构**：物理上有，但当前可测、不受nuisance模仿的有限保证未完成。
- Therefore measuring **当前三频双探针** makes previous worlds distinguishable：全参数域尚不能推出。
- New sensing rule：本轮只支持**不要丢弃强参考/直接路径、必须同时检验材料与校准以及模型失配**，尚非新颖且已证明胜出的天线设计法则。

可交付成果是完整比较、明确负结果和有限数学约束，不是一篇已可投稿的强方法论文。继续相同参数上的后验调参、把普通联合反演改名或扩展到UAV，不会改变这一裁决。
