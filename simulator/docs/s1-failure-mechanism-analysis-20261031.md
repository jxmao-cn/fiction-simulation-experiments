# S1 失败的离线分析：罚约束、旧方向与浮点停滞

日期：2026-10-31。本文结合现有源码、M2 CSV和**标量算术反例**分析罚约束与浮点停滞。

## 1. 判定

| 命题 | 判定 | 依据 |
|---|---|---|
| 罚刚度加倍后仍使用旧方向 | 成立，源码事实 | 更新刚度与能量后，没有在本次线搜索前重新求方向 |
| 旧方向可能在新目标下变成上坡方向 | 成立，数学反例 | 第 4 节给出确实满足增刚度条件的一维例子 |
| 该方向问题就是 M2 此次失败的原因 | 存疑 | 需补充每轮方向、刚度、能量和步长 |
| `<=` 能接受坐标完全不变的 trial | 成立，浮点反例 | 第 5 节核对 binary64 的舍入值 |
| M2 实际反复接受了零位移 | 存疑 | CSV 未记录 alpha、trial 与旧坐标是否相等 |
| 每个载荷调用都从初始罚刚度重新开始 | 成立，源码事实 | DBC 按值复制；runner 没有写回局部增大的刚度 |

阅读顺序建议：先看第 2 节理解程序实际在解什么；第 4 节解释“方向失效”；第 5 节解释“假进展”；第 6–7 节再回到 M2 数据。

## 2. 已知机制：把硬边界暂时变成弹簧，再逐步加硬

### 2.1 目标函数、梯度与 Hessian

S1 没有惯性和重力，真正希望求解的是弹性能在给定端面位移下的平衡。令所有坐标展平为 $x$，约束节点选择矩阵为 $A$，目标端面坐标为 $b$，约束节点质量组成对角矩阵 $W$。当前先解带罚项的问题：

$$
\Phi_k(x)=U(x)+\frac{k}{2}(Ax-b)^TW(Ax-b).
$$

因此

$$
r=Ax-b,\qquad
g_k=\nabla U(x)+kA^TWr,
$$
$$
H_k=\nabla^2 U(x)+kA^TWA.
$$

代码将材料Hessian做PSD投影，求方向使用的近似矩阵记为$B_k$。

| 符号 | 含义 |
|---|---|
| $U,\Phi_k$ | 材料弹性能、含边界罚项的当前目标 |
| $x,A,b,r$ | 坐标、约束选择矩阵、约束目标、边界误差 |
| $W,k$ | 约束节点质量权重、DBC 刚度参数 |
| $g_k,H_k,B_k$ | 目标梯度、精确 Hessian、实际求方向矩阵 |
| $p,\alpha$ | 求出的完整搜索方向、线搜索缩放系数 |
| $\tau,h_{eff}$ | 代码停止容差、有效时间步；静态时后者固定为 1 |

对应源码：[罚能量](../cpp/src/time_integrator.cpp:124)、[罚梯度](../cpp/src/time_integrator.cpp:138)、[罚 Hessian](../cpp/src/time_integrator.cpp:152)。每节点实际罚刚度为$km_i$。静态去惯性见 [目标组装](../cpp/src/time_integrator.cpp:174)，有效步长定义见 [time_integrator.h](../cpp/include/time_integrator.h:122)。静态有效步长1为算法约定。

直观上：端面通过很硬的弹簧被拉向目标，弹簧仍可能伸长；弹簧越硬，位移误差越小，却通常也越难数值求解。

### 2.2 “当前罚问题解得差不多”和“边界准确到位”是两回事

方向求解在未冻结的自由度子空间上近似满足

$$
B_kp=-g_k.
$$

代码使用的方向率与边界满足条件为

$$
\rho_p=\frac{\|p\|_\infty}{h_{eff}},\qquad
\frac{\|x_i-b_i\|_2}{h_{eff}}<\tau.
$$

两者都要满足才退出。如果 $\rho_p\le\tau$，但尚有端面节点未到位，就将 $k$ 加倍。

小反例说明为什么需要增刚度：对

$$
\Phi_k(x)=\frac12x^2+\frac{k}{2}(x-b)^2,
\qquad x_k^*=\frac{k}{1+k}b,
$$

在 $x_k^*$ 处梯度和 Newton 方向都是零，但边界误差仍为 $-b/(1+k)$。例如 $b=1,k=1$ 时，当前软约束问题已经精确收敛在 $x=0.5$，距离硬目标 1 仍差 0.5。增大罚刚度会改变软约束目标及其平衡位置。

一般地，在罚平衡中 $k m_i(x_i-b_i)$ 与材料节点力平衡，误差通常随 $1/(km_i)$ 缩小。网格改变质量分布后，同一 $k$ 的实际边界刚度也会变化，比较网格时需要计入这一变化。

条件位置：[DBC 满足判定](../cpp/src/time_integrator.cpp:315)、[主循环退出条件](../cpp/src/time_integrator.cpp:432)、[增刚度条件](../cpp/src/time_integrator.cpp:444)。这里的residual定义为方向率。因为 $p=-B_k^{-1}g_k$，仅观察小 $p$ 不能脱离矩阵尺度判断 $g_k$ 有多小。

### 2.3 满足后会暂时冻结节点

已满足的整个 DBC 节点会被清零梯度、将矩阵行列消元并置单位对角，令对应方向为零：[梯度清零](../cpp/src/time_integrator.cpp:324)、[行列消元](../cpp/src/time_integrator.cpp:335)。因此它保持当前容差内的近似位置。只要目标未变、节点未再移动，之后仍满足这一判定。

这解释了为什么“50 个边界节点”与“150 个约束自由度”需要区分：CSV 的 `dbcSatisfied` 按节点计数，每节点用三维二范数判定。

## 3. 从代码顺序看：目标换了，方向还没换

当前关键次序如下（省略无关部分）：

```cpp
// 416–418：按当前 k 求能量、方向
E_last = computeIPVal(..., dbc, ...);
sr = computeSearchDirection(..., dbc, ...);
p = sr.p;

// 444–448：改变 k，并更新比较基准
if (res <= cfg_.tol && DBC_not_all_satisfied) {
    dbc.stiffness *= 2.0;
    E_last = computeIPVal(..., dbc, ...);
}

// 451–465：沿旧 p，在新目标下线搜索
xTrial = x + alpha * p;
E_trial = computeIPVal(xTrial, ..., dbc, ...);
if (isfinite(E_trial) && E_trial <= E_last) accepted = true;

// 479–484：先更新坐标，才重新求方向
x += alpha * p;
sr = computeSearchDirection(..., dbc, ...);
p = sr.p;
```

代码更新了能量比较基准，随后沿旧方向搜索新目标。

## 4. 假说 A：旧方向可能变成新目标的上坡方向

### 4.1 推导

固定当前坐标和冻结集合，由 $k$ 改为 $2k$ 后：

$$
g_{2k}=g_k+kA^TWr.
$$

假设旧矩阵在当前可动子空间上正定，则旧方向对旧目标下降：

$$
g_k^Tp=-p^TB_kp<0\quad(p\ne0).
$$

但在新目标下：

$$
\boxed{g_{2k}^Tp=-p^TB_kp+k\,r^TWA p.}
$$

第二项可能为正并超过第一项。于是旧方向不再具有下降保证。仅仅继续减小 alpha，不能保证找到降能步；教材关于回溯存在性的前提本来就是方向对**当前目标**下降，见 [本地教材 Line Search](../../../Archive/phys-sim-book/src/simulation-with-optimization/projected-newton/line_search.md)。

### 4.2 确实触发增刚度分支的数值反例

使用归一化的一维模型，弹性刚度和质量权重均取 1，分析刚度更新与方向之间的关系。令

$$
\Phi_k(x)=\frac12x^2+\frac{k}{2}(x-b)^2,
\quad x=3\times10^{-10},\ b=5\times10^{-10},\ k=1,
\quad h_{eff}=1,\ \tau=10^{-10}.
$$

| 量 | 旧目标 k=1 | 新目标 k=2 |
|---|---:|---:|
| 梯度 | 1e-10 | -1e-10 |
| Hessian | 2 | 3 |
| 正确 Newton 方向 | -5e-11 | +3.333333e-11 |

旧方向大小 5e-11 已不大于容差，但边界误差 2e-10 尚未小于容差，**正好触发代码中的刚度翻倍**。

方向的含义很直接：旧目标希望 x 变小，新目标更强调靠近右侧目标 b，希望 x 变大。两者方向相反。

沿旧方向，新的方向导数和精确二次能量差为

$$
g_2p_{old}=5\times10^{-21}>0,
$$
$$
\boxed{\Phi_2(x+\alpha p_{old})-\Phi_2(x)
=5\times10^{-21}\alpha+3.75\times10^{-21}\alpha^2>0
\quad\forall\alpha>0.}
$$

全步时，能量由 8.5e-20 升到 9.375e-20。**精确算术下没有任何正的可接受降能步**；附带的 Python 标量 binary64 算术脚本对照也在回溯 47 次后碰到 alpha<1e-14 的停止条件，没有接受 trial。

这足以反驳“给旧方向缩步，总能适应新刚度”的普遍断言。本例触发线搜索失败，M2则记录`lineSearchFailed=0`；诊断M2还需要实际方向trace。

## 5. 假说 B：缩步后坐标或能量不再有可辨变化

### 5.1 实际位移的舍入

实际试探增量是 $\alpha p$，而报告记录的是 $\|p\|_\infty/h_{eff}$，没有记录 alpha。binary64 默认就近舍入时，1 向上相邻数的间隔为

$$
\operatorname{ulp}(1)=2^{-52}=2.220446049250313\times10^{-16}.
$$

取一个量级接近 M2 最终方向率的正方向作为**构造值**：

$$
x=1,\quad p=1.011774232477666\times10^{-10},
\quad\alpha=2^{-20}=9.5367431640625\times10^{-7}.
$$

于是

$$
\alpha p=9.649030995155964\times10^{-17}
<\tfrac12\operatorname{ulp}(1).
$$

实际核对结果是 `x + alpha*p == x` 为 True。这是坐标级舍入。

### 5.2 当前接受条件会让这个 trial 通过

构造非驻点目标

$$
f(z)=\frac12(z-1)^2-p(z-1).
$$

在 z=1 处，$f'(1)=-p\ne0$，Newton 方向恰为 p。如果已经试探到上述 alpha，trial 坐标与原坐标完全相同，计算所得两次能量都是 0，因此

```text
trial_equals_x = True
E_trial = E_last = 0
E_trial <= E_last = True
```

随后坐标没有有效变化，下一轮可能算出相同方向，而方向率仍高于 1e-10；如此可消耗迭代上限。源码接受位置为 [第 463 行](../cpp/src/time_integrator.cpp:463)，判断为 [第 465 行](../cpp/src/time_integrator.cpp:465)，实际更新为 [第 479 行](../cpp/src/time_integrator.cpp:479)。代码没有检查更新后的坐标是否相同。

该二次函数的完整Newton步下降。上述计算检查小alpha试探到达后，接受条件如何处理零有效更新；实际回溯路径由逐轮trace检查。

### 5.3 能量相等甚至不要求坐标相等

另取

$$
\widetilde f(z)=1+\frac12(z-1)^2+p(z-1).
$$

在 z=1 沿正 p 取 alpha=1，是严格上坡：精确能量增加 $1.5p^2\approx1.54\times10^{-20}$。但能量在 1 附近的浮点间隔约 2.22e-16，所以本机计算为：

```text
x_trial != x
E_old = 1.0
E_trial = 1.0
```

当前`<=`接受该trial。这个构造例子展示了背景能量尺度与消减误差如何使严格上坡的能量差舍入为零。

浮点背景参考：[Python 官方浮点说明](https://docs.python.org/3/tutorial/floatingpoint.html)、[math.ulp 定义](https://docs.python.org/3/library/math.html#math.ulp)。本次数值全部由本机 elastogen 环境直接计算，见第 9 节。

## 6. M2 数据与待补充诊断

来源：[M2.csv](../output/validation/s1-20261030-01/M2.csv)，只读其 1.02 行：

| 指标 | 数值 | 正确解释 |
|---|---:|---|
| iterations / maxIterHit | 100 / true | 已达到上限，S1 未通过 |
| DBC 满足 | 1/50 | 大部分端面节点仍未到位 |
| DBC 二范数最大误差 | 7.785367e-6 | 高于静态节点判据1e-10约7.8e4倍 |
| 搜索方向率 | 1.011774e-10 | 比 1e-10 高约 1.18%；定义为搜索方向率 |
| 自由梯度无穷范数 | 5.475277e-7 | 自由节点上的物理梯度 |
| 累计回溯 / lineSearchFailed | 3474 / false | 每轮最终找到数值意义下不增能的 trial |
| 端面无穷范数误差 | 7.785222e-6 | 与上面的节点二范数指标略有区别 |
| 能量 / 反力相对误差 | 8.03e-4 / 4.13e-4 | 均超过 1e-5 的验收阈值 |

平均每更新有 34.74 次回溯。如果某轮起始 alpha=1，回溯 j 次对应 $\alpha=2^{-j}$；但是 $2^{-34.74}\approx3.4851\times10^{-11}$ **不是实际观测到的某轮 alpha，也不是所有轮 alpha 的算术平均**。本场景 NH 的初始可行步长返回 1，见 [NH initStepSize](../cpp/include/energies/NeoHookeanTetEnergy.h:134)；即便如此，汇总平均仍不能告诉我们末段每轮的步长。

实际更新需由同一轮的p与alpha计算。

当时停止的直接原因确定是：方向率和 DBC 未同时达标，直到耗尽上限。A、B、罚项尺度及其相互作用由后续逐轮记录区分。

## 7. 另外确认的机制：罚刚度不会跨载荷保留

`stepForward` 接收 const 引用，却在内部执行 [按值复制](../cpp/src/time_integrator.cpp:381)：

```cpp
DBCInfo dbc = dbcIn;
```

倍增只改这个局部副本。runner 在 [第 181 行](../cpp/examples/tension_check.cpp:181) 把外部刚度设为 1e8；每个载荷在 [第 228 行](../cpp/examples/tension_check.cpp:228) 再传入外部 dbc，只取回位置/速度，没有把 `state.report.dbcStiffness` 写回。

因此，真实记录应这样读：

| 调用 | 起始刚度 | 末尾报告刚度 |
|---|---:|---:|
| M2，1.01 | 1e8 | 2.62144e13，即起始值的 2^18 倍 |
| M2，1.02 | **重新从 1e8 开始** | 2e8，即本次增大一次 |

位置采用warm start，罚刚度则跨调用重置；新载荷会重新经历较软约束阶段。后续可比较继承与重置罚刚度的求解路径。

报告更新另有一个时序问题：增刚度后若线搜索立刻失败，在 [第 471 行](../cpp/src/time_integrator.cpp:471) 已退出，而 `report.dbcStiffness` 的更新在 [第 488 行](../cpp/src/time_integrator.cpp:488)，报告可能仍为旧值。本次M2记录lineSearchFailed=false。

## 8. 逐轮诊断设计

下一次受控诊断需要记录每个 Newton 迭代：

- 更新前/后的 k、方向属于哪个 k、方向率、DBC 满足数和最大误差。
- 当前目标下的方向导数 $g_k^Tp$，以及材料/罚能量分量。
- 线搜索起始 alpha、最终 alpha、回溯次数、实际坐标增量和 trial==x。
- 能量差、参考的数值尺度，以及本轮是否触发增刚度。

判别标准：若加倍后 $g_{2k}^Tp_{old}>0$，才能证明本轮方向失效；若 accepted 且 trial==x，才能确认零有效更新；若只看到能量相同而坐标不同，还要检查能量分辨率和原始梯度。

最小修订比较应先只改变“增刚度后重算方向”，其余参数和验收标准固定；再独立比较停滞检测、罚参数继承或明确的硬约束方案。重算方向只在当前矩阵正定等条件下恢复下降性质，不自动保证全局收敛；改成严格 `<` 也不单独解决浮点分辨率问题。

以上为后续诊断设计。

## 9. 可复算材料

- [标量算术脚本](../output/validation/s1-analysis-20261031/check_counterexamples.py)，仅用标准库，不导入仿真器。
- [实际算术输出](../output/validation/s1-analysis-20261031/arithmetic.json)。
- [上一轮完整失败记录](s1-validation-20261030.md)。

```bash
/home/steven/miniconda3/envs/elastogen/bin/python \
  /home/steven/mjx/fiction/simulator/output/validation/s1-analysis-20261031/check_counterexamples.py
```

源码定位以本轮读取的当前版本为准；主要线搜索和 DBC 逻辑与执行副本一致。读取时 `time_integrator.cpp` SHA-256：`84ee9de8ed6118c748d5625b1ddbea7f4b632cb909ac7f8c22978553a9216dc6`；M2.csv SHA-256：`0458e25ade860dba8aa85cc5d0f120a4b8a0189925aa5e951619aec36a95cf95`。
