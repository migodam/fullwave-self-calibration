# A5 可复现性与证据归属

版本2026-09-10。本包是新的Q5执行检查点，不覆盖附件中已有的A5收尾补充。总科学裁决仍为 **Q5尚未完成**。复现成功只验证代码/记录一致，不改变科学Gate。

## 1. 输入与版本

用户仓库main在此前连接读取时为commit `c9138b63c488ad3422d286c78a114936ed210c56`（2026-09-09归档）。本轮更近的依据是用户上传的2026-09-10 v2交接zip和Q5文本。旧A4报告、原始结果、协议、修正和代码从zip读取；内部A4zip另行解包。全部88个提取输入文件的SHA256前后相同，记录于`results/input_manifest.json`和`results/input_integrity_final.json`。

仓库本次未修改、提交或推送。未自动联系作者、投稿、购买服务、配置常驻设施或调用外部模型provider。没有执行真实Agentic-AI-Scientist pipeline；普通脚本不称为该pipeline。Treams不在本轮执行环境，未运行；旧Treams记录只作为继承比较。

## 2. 环境与依赖

实际：Python3.13.5；NumPy2.3.5；SciPy1.17.0；mpmath1.3.0；pytest9.0.2；Linux x86_64。容器可见5个逻辑CPU，但执行时显式设BLAS/OMP/MKL为1线程。没有GPU计算。完整环境字符串在`results/environment.json`。

`requirements.txt`给出本轮实际版本。复现需要可安装这些依赖的环境；没有保证所有操作系统/旧Python都可直接安装。主要代码不依赖Treams、ADDA、API密钥或网络。`original_maxwell3d.py`保留了一个未调用的Treams函数，不意味着需要安装Treams。

## 3. 安全的完整复现命令

在解压后的包目录中执行：

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python code/reproduce.py --output ../A5_Q5_fresh_run
```

输出目录必须不存在，且不能放在原证据包内部。包装器复制代码和协议，创建新的results/logs，然后逐项执行；若任何命令失败则停止。不要直接在证据包内逐个执行原脚本，否则可能覆盖JSON。原脚本的写出路径由自身位置确定，安全包装器负责隔离。

仅运行核心检查、不运行随机恢复和DDA网格检查：

```sh
python code/reproduce.py --output ../A5_Q5_fresh_quick --quick
```

交付前这两个入口都在全新目录实际跑通。完整新跑得到9个与原包共有的JSON：忽略所有含`seconds`的计时字段后，5,251个数值完全一致，其他非计时内容也一致。该复现没有增加独立场景数，使用同一随机种子。审核记录在`results/clean_reproduction_audit.json`。

包装器不重新生成文献、论文和报告，也不伪称重新做数学证明。它不依赖本机之前的绝对输入路径，因此能独立重跑本包新数值。原A4的21个pytest检查需要用户原交接包中的A4源目录；本包保留本轮运行日志，但不复制整个旧研究归档。

## 4. 文件和入口

| 路径 | 作用 / 执行性质 |
|---|---|
| `code/spherical_maxwell.py` | 本轮新写3Dvector球谐/精确球系数/相互作用正演 |
| `code/original_maxwell3d.py` | 字节保留的附件原DDA/Green kernel/体素化代码 |
| `code/check_solver.py` | 多极/积分/弱对比/能量检查，旧六尺度独立实现复核 |
| `code/mechanism_ablation.py` | 物理Born、孤立、单物体、多极、相互作用、gain结构 |
| `code/check_dda.py` | 四档独立体积分离散与多极比较 |
| `code/recovery.py` | 冻结18场景、60方法场景、180启动；全部失败保留 |
| `code/supplement.py` | 后冻结有限竞争对、完整VarPro导数、GLS等价、多频gain |
| `code/modal_boundary.py` | 精确Mie世界、5500区间盒、风险和一实reference界 |
| `code/model_stress.py` | 距离/两种失配、正确单球固定损耗补充 |
| `code/tests.py` | 19个回归断言，不是19个独立实验场景 |
| `code/reproduce.py` | 新目录安全完整/quick复现入口 |
| `A5_FROZEN_PROTOCOL.md` | 随机恢复前本地注册协议；机制开发在其之前 |
| `results/*.json` | 原始结果、每个初值、场景真值、参数误差和时间 |
| `logs/*.log` | 实际运行输出；`a4_tests.log`为旧21测试本轮执行 |
| `MANIFEST_SHA256.json` | 交付文件内容校验；不包括自身 |

所有相对链接以包根为基准。新文稿不覆盖旧A4稿件。没有附带字体、凭据、文献付费全文或第三方二进制软件。

## 5. 随机性、数据预算和场景划分

主恢复seed202609105：先12个same-model，再6个independent-DDA，保持一个连续随机流。每场景三初值，不产生额外统计样本。相同噪声用于该场景方法配对；有错误reference的场景复用reference随机噪声、增加已声明的幅相偏差。

主观测144复数；一个正确/错误复参考std=.01。validation为16个未用接收位置；structural场为固定36点、无电子学gain的三分量散射场。真值仅用于评估，不供优化。

受限M类噪声诊断seed202609107：两个精确世界，各200次，共400次。一个实gain模长参考std=.001，与主实验的复reference不是同一个信息预算。

失配名义点seed202609108。频率/范围/机制为确定性开发或post-frozen诊断，不称新冻结统计样本。全套重跑复用这些种子，只作可复现性检查。

## 6. 计时定义

代码在导入数值库前设置单线程。`total_seconds`包含一个方法场景的三启动、端点评估等；`data_generation_seconds`包含该场景所需传感器、结构场和验证场生成；脚本`wall_seconds`还包含初始化和写文件。

主运行记录约12.256s；数据生成合计1.837s、全部方法三启动fit合计9.789s。不同运行时间会变，校验明确忽略时间值。没有把代码开发、文献检索、真实硬件、参考/模式校准及测量时间算作零成本，也没有从这些小场景计时推导普遍节时优势。

## 7. 数值证书边界与已知问题

多极阶数和积分差是敏感性检查，不是全材料域连续误差包络。DDA加密非单调；其结果只能称独立离散压力测试，不能冒称认证的连续真值。卡方残差门是heuristic，未给全局覆盖。

类M使用`mpmath.iv`35位区间算术覆盖完整材料区间，而非抽样导数正性。边界端点由十进制有理串构造。解析公式另由差分点检核对。尚未用另一interval backend或形式化证明系统验证，也没有将天线模态误差纳入硬件置信度。

`boundary_check`注释暗示更强的角网格检查，实际是模态界面代数残差，报告已按实际功能降级。开发单球第二物体损耗误用.03，原结果保留，`model_stress.json`另存正确.05的结果。原plug-in reference-first残差门忽略reference噪声；后冻结证明了正确GLS与joint目标等价，撤回基于该弱门的算法优势解释。

回归断言只检查其指定条件，不认证全部论文或真实Maxwell系统。所有参数、噪声及gain范围的变更都应另建开发协议/输出目录，不能反复修改冻结集直到显著。

## 8. 证据标签与归属

原报告：A1–A4统计与旧Treams运行。独立重算：旧六尺度的新多极实现、旧A4测试。新证明：Born固定损耗、明确比例机制、精确Mie有限风险/一实reference、正确GLS/导数和nuisance扰动条件。新数值：本包JSON。待验证：全域主任务误差/分离、现实模态/参考、强设计政策、公平SOTA及硬件。

新代码中复用的原文件保持来源可追溯；本包没有给旧代码或文献重新指派许可证。使用时仍应遵守原仓库及依赖的许可条件。