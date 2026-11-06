# S6a 单体主工程重验规格

日期2026-11-02。本批执行单体10步。使用当前生产fiction_core，测试harness复用隔离副本的场景辅助实现，不拷贝旧核心。

固定sphere1K，1760节点/6851tet，63个局部z>=0.45顶部DBC；NH E=1e5 Pa、nu=0.3、rho=1000，重力-z9.81。dt=.005 s、10步、tol=1e-8、maxIter=100、DBC k=1e8；combined，gradTolRel=1e-6、gradScale=单体重约5086.259121 N、actualStepTol=1e-12 m、连续2次。

每步要求收敛、全部状态有限、min J>.9、DBC全部满足且误差<1e-6 m、物理残差/单体重<1e-6，失败立即停止并保存末态，不重试。与原S6a单体及修正隔离单体的完整末态最大节点距离均<1e-8 m。combined仍含原step快路径，但本批单独核查每步力残差。

单批最多600 s、地址空间8GiB。新run-id `s6a-single-main-20261102-01`，保存源码/生产核心哈希快照、编译/测试/命令、逐步CSV、末态、trace及比较结果。实际小步指逐Newton迭代trace中的actual_dx_inf；CSV的state_change_inf是整个时间步节点位移。
