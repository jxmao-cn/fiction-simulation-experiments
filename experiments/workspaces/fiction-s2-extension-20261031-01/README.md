# S2 隔离参数实验源码

日期：2026-10-31。基于当前 Archive/fiction 工作区复制，未修改或合并生产源码。cpp/build 排除。

实验修改仅四个文件，详见 experiment-only.patch：

- Config 增加 initialStepScale，默认1.0，构造验证范围 (0,1]。
- 每轮线搜索从 initialStepScale 开始，其余候选步长、能量接受、减半、停止规则不变。
- cantilever_check 添加 --max-iter（默认100）、--tol（默认1e-10）、--initial-alpha（默认1）。实际参数记录在CSV注释；材料、几何、重力、DBC和初始参考形状不变。

initial-alpha=1 保留原始路径；小于1用于受控阻尼对照。这里的参数为每轮线搜索起始上限，不是固定接受步长，亦不是深度学习的梯度下降学习率。未加入初始弯曲形状或加载步进，避免混合变量。

构建目录：/tmp/fiction-s2-extension-20261031-build；Release，CHOLMOD。
构建命令与输出见 configure.log、build.log。物理实验由主会话按固定批次另行执行，本子任务未运行物理批次。

source-hashes.json 为隔离源码哈希；scope-check.json 记录副本相对原代码的修改范围。

高分辨率网格自检：runner 以 long double 独立累加正定向四面体参考体积，保持原1e-12相对自检门槛；不改 builder 的质量或力学数据。原 model.volume 保留原 volume 列作诊断，独立求和写 verified_reference_volume。网格生成器采用 math.fsum 降低几何体积的累加误差。

实际 CLI 同时写 configured_tol / configured_max_iter / configured_initial_alpha 列。
