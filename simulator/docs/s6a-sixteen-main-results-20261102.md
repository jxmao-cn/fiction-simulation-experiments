# S6a 16体条件档：主工程短段验证

日期2026-11-02。结论：**本次16体10步条件档通过预登记门槛**。四体前置已通过，规格见 [运行前规格](s6a-sixteen-main-spec-20261102.md)。未修改生产源码、原run或旧验收判定；本批不覆盖接触、异质体和长时稳定。

## 核心来源与配置

通过add_subdirectory编译生产simulator/cpp，链接production/libfiction_core.a，生产include优先。link-command/header-dependencies和核心库/源码/资产/二进制哈希均保存。Release/CHOLMOD，CTest11/11通过（含2项GUI占位）。

16个sphere1K同质实例，4×4平移网格、间隔1.5m；28160节点、109616tet、1008个顶部cap DBC。NH E1e5Pa、nu.3、rho1000，重力-z9.81；dt.005，10步共.05s，tol1e-8、maxIter100、DBC k1e8。combined显式选择，gradTolRel1e-6，单体力尺度5086.259121N，actualStepTol1e-12、连续2次。

## 数值与组织检查

- 10/10步收敛且门槛通过；Newton次数[2, 2, 2, 2, 3, 4, 4, 7, 5, 3]。无maxIter或线搜索失败，状态有限。
- min J最低0.988778770461；DBC误差0，全部约束满足；最大相对物理力残差3.1638742e-10，低于1e-6。已独立按raw IP梯度/dt²复算导出残差一致。
- 初末跨体Hessian triplet、装配后entry及跨体tet均0；初末矩阵非零均3089376，等于16×单体193086。
- 11个保存采样时刻（初态和10步）bbox两两分离下界最小0.5015961m；初态分离在求解前检查。分离检查覆盖这11个采样时刻。
- canonical参考只实际独立求解一次；同质实例仅平移、局部载荷/材料/DBC相同，复用其完整11时刻轨迹作为16体各自局部问题参照。
- 逐体逐时刻176条比较全部完成，最大位置差2.17904198983e-10m，低于1e-8。完整末态另由CSV独立重算，与canonical最大差2.17904198983e-10m，DBC独立复算误差0。176条全时窗对照由程序执行并保存比较表；独立CSV重算覆盖完整末态。

## 成本与证据

进程单次墙钟52.83s，16体十步求解/trace归档合计48.80459s，每步平均4.88046s，中位数4.48409s，最大8.76492s；峰值RSS576420KiB（约562.91MiB），无swap，均低于600s/8GiB预算。线程环境变量设为1，但GNU time显示CPU利用率177%；本次为一次带trace运行，GNU time记录反映实际进程资源使用。

run-id：[s6a-sixteen-main-20261102-01](../output/validation/s6a-sixteen-main-20261102-01/)。包含spec/config、源码快照/hash、真实core链接/header依赖、命令、构建回归、二进制、逐步CSV、全末态、每步独立trace、块对角CSV、bbox间隔、176条对照、summary与time/RSS。无失败重试、未放宽门槛。

本批成功支持主工程可选combined在同质无接触16体短段上的运行、装配及数值等价性；不据此改写原旧run，也不推广到摩擦、接触、异质体或连续体精度。
