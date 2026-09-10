# 代码归属与证据边界

`dda.py`逐字复制自本项目基线 `research/trispace_self_calibration/a3_research/maxwell3d.py`，SHA256 bb401b4719e5e74e4dd5614ca566e85294e75e32c9d893540eb75e289958030b。只读复用标准DDA，未更改A1–A5。

`vsw.py`来自用户提供的Q5_V2_VERIFIED_RESEARCH.zip中的经典多球VSW实现，SHA256 8e960c2052fa126427ec203eeefabc955e6ea554369dda7af3cf5fbbf6dfe509；Git blob为546046bdc41ff61603f7104a5468b474b567bdab。二者不作为新方法贡献。冻结协议中该哈希的手抄缺字见09_REVIEWS/DELIVERY_STATUS.md，不修改冻结历史。

`model.py`新增实际点源投影、移动接收坐标、三频本构sharing、准确线性Born和跨频增益；完整附件中的`experiment.py`实现B1冻结流；`safe_run.py`只修正浮点初始边界；`diagnostics.py`新增因果前缀、弱模型及oracle；`reference_only.py`实现事后D2两阶段校准；`checks.py`、`test_core.py`检查实现；`summarize.py`只汇总原结果。

**GitHub因主实验脚本上传被拦截而仅有部分源码，完整执行必须使用对话交付包。** 不声称本PR是附件镜像。

NumPy/SciPy由各自项目授权，未在附件中打包第三方库或书籍。没有宣称本仓库全部内容具有未确认的开源许可证。代码保留在用户项目中用于科研复核。普通profiling、频率continuation、多极T-matrix和DDA都不是本轮首创。

所有正演返回有限数值近似。任何函数都没有实现连续Maxwell误差认证；`reject_diagnostic`仅是声明噪声下的目标残差报警，不是覆盖定理或通用可信度分数。
