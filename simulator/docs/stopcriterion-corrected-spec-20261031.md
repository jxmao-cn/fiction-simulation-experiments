# 修正停止判据：隔离诊断规格

日期：2026-10-31。用户授权按独立审阅建议继续；不改生产，不覆盖旧run，不更改原S2/S6a验收结果。

run-id：`stopcriterion-corrected-20261031-01`；源码位于 `../../experiments/workspaces/fiction-stopcriterion-corrected-20261031-01/`。项目根由fiction.project.json定位。

- 保留raw freeGradientInf字段；新增forceResidualInf，动态为raw/h²，静态为raw；grad与combined按物理力残差[N]判定。
- combined使用gradTolRel配置值，默认step模式保留。
- S6a CLI尺度解释为单体力尺度：1体、4体、独立单体统一使用同一单体总重；本诊断为同质复制，异质体还需逐体归一化，不在本次声称支持范围。
- 单元负例：dt=.005/.01/.02的相同初态力残差为1N；实际更新人为压小而力残差大时必须拒绝combined；改动梯度容差必须改变相应决策；4份独立复制不改变判定。
- Release/CHOLMOD，先CTest；负例失败不运行物理病例。
- M3保持tol=1e-10、maxIter=100、static；S6a保持tol=1e-8、dt=.005、10步、maxIter=100。combined gradTolRel=1e-6、actualStepTol=1e-12、连续2步。材料、边界、惯性、能量和线搜索不改。
- 运行M3、S6a1体、S6a4体（包含独立单体），每次最多600秒、地址空间8GiB，总15分钟；仅诊断允许数值失败后保留现场并运行预登记下一独立病例，资源/执行异常停止。
- 保存原始CSV、完整末态、命令、快照、diff、哈希、CTest、时间/RSS。与旧末态比较最大节点位置差，阈值1e-8m；已有old combined结果仅为数值一致性参考，不是物理真值。
- 测试通过仅说明本次修正和病例有效，不自动合入生产、批准v0.3或宣布S2的M3精度通过。M6与16体不在本批。
