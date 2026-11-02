# IPC 模型 ↔ fiction C++ 工程连接指南（设计思路与理论支持）

> 面向学习的技术文档，记录 2026-09 这次改动的全部内容：如何让 fiction（mass-spring 增量势能模拟器）
> 直接加载并仿真 IPC（ipc-sim/IPC, SIGGRAPH 2020）`input/tetMeshes/*.msh` 下的四面体模型，
> 实现「重力 + 边弹簧弹力 + 坐标盒区域单一外力」，并在 SSH 无显示环境下导出 OBJ 与 10 秒 mp4 视频。
>
> 阅读顺序建议：第 1 节（总览）→ 第 2 节（网格格式）→ 第 3 节（mass-spring 装配）→ 第 4 节（能量与牛顿法）
> → 第 5 节（场景/边界/外力）→ 第 6 节（可视化与视频）→ 第 7 节（验证）→ 第 8/9 节（文件清单与边界）。
>
> 终端环境无法渲染 LaTeX，本文所有公式均用纯文本/代码块表示。
>
> IPC 全部模型资产的逐文件清单（格式/顶点/边数/可可视化性）见 [`ipc-asset-inventory.md`](ipc-asset-inventory.md)。

---

## 1. 总览：一条数据流水线

核心目标只有一个：**把 IPC 的几何（网格）喂给 fiction 的物理（能量+牛顿求解）**，两端各自保持独立，
中间用一层薄薄的适配模块连接。整条流水线：

```
 IPC/input/tetMeshes/*.msh                 （外部几何资产，MshIO 4.1 ASCII 四面体网格）
        │
        │  ① 解析（纯文件 IO，无第三方依赖）
        ▼
 cpp/include/io/msh_loader.h        ──►  MshMesh { x, tets, surface }
        │
        │  ② 装配（几何 → 物理参数）
        ▼
 cpp/include/scene/mass_spring_builder.h
        │                           ──►  MassSpringModel { x0, masses, edges, l2, k, bbox }
        │
        │  ③ 场景化（命令行 → DBC / 外力 / 初始条件）
        ▼
 cpp/include/scene/ipc_scene.h      ──►  IpcScene { mesh, model, dbc, pinnedNodes, forceNodes }
        │
        │  ④ 求解（复用现有框架，不新造轮子）
        ▼
 TimeIntegrator（隐式欧拉 + projected Newton + CHOLMOD/LDLT）
   ├─ Inertia（内置）
   ├─ ExampleGravityEnergy        ← 重力
   ├─ ExampleSpringEnergy         ← 边弹簧弹力（tet 边）
   ├─ RegionForceEnergy           ← 坐标盒区域恒定外力（对标 IPC NBC）
   └─ DBC 罚能量                  ← pin-box 零速固定
        │
        │  ⑤ 输出
        ├─► cpp/include/io/obj_writer.h          → 表面 OBJ 帧（MeshLab/Blender 可看）
        └─► cpp/include/io/rasterizer.h
            + cpp/include/io/image_writer.h
            + cpp/include/io/video_encoder.h     → PNG/PPM 帧 → ffmpeg → mp4（SSH 可看）
```

**设计原则**（贯穿全部改动）：

1. **不引入 IPC 的依赖树**。IPC 本体依赖 libigl / OSQP / TBB / CHOLMOD / CCD / MshIO……我们只需要
   `.msh` 的几何，就自己写一个几十 KB 的 ASCII 解析器，而不是把 MshIO 拉进来。
   （类比：只需要读一张 Excel 表格，不必把整个 Office 装到项目里。）
2. **物理框架零改动**。`TimeIntegrator`、`IEnergy`、`DBCInfo`、CHOLMOD 求解链路一行没动；
   新东西全是 `IEnergy` 实现 + 头文件工具 + 示例程序。这正是 `energy.h` 设计接口的价值。
3. **几何与物理分离**。`msh_loader` 只管「坐标/单元/表面」，`mass_spring_builder` 只管
   「边、质量、刚度」，`ipc_scene` 只管「命令行与场景语义」。每层可单独测试。
4. **控制台与视频也要无显示环境可用**。SSH 场景没有 OpenGL，所以视频走 CPU 软光栅化，
   仅把 ffmpeg 当编码器。

---

## 2. 读懂 IPC 的模型格式：MshIO 4.1 ASCII

### 2.1 文件长什么样

以 `IPC/input/tetMeshes/cube.msh` 为例（这是最小的真实样本，8 节点 6 四面体）：

```
$MeshFormat
4.1 0 8                 ← 版本 4.1 / fileType=0（ASCII）/ sizeof(size_t)=8
$EndMeshFormat
$Nodes
1 8 1 8                 ← 1 个 entity block，共 8 节点，tag 范围 [1,8]
3 0 0 8                 ← block: entityDim=3, entityTag=0, parametric=0, 8 个节点
1 2 3 4 5 6 7 8         ← 先列出全部节点 tag（每行一个）
0 0 0                   ← 再按相同顺序列出坐标（每行 x y z）
1 0 0
...
$EndNodes
$Elements
1 6 1 6                 ← 1 个 block，共 6 单元
3 0 4 6                 ← entityDim=3, entityTag=0, elementType=4（四面体）, 6 个单元
1 1 5 7 6               ← 单元 tag，4 个节点 tag
...
$EndElements
$Surface
12                      ← IPC 附加段：12 个边界三角面
1 5 6                   ← 每行 3 个节点 tag
...
$EndSurface
```

关键事实（实测 IPC `input/tetMeshes/` 全部 38 个文件）：

| 事实 | 数据 | 设计含义 |
|---|---|---|
| 全部 ASCII | `4.1 0 8` 的第二个数总是 0 | 只需支持 ASCII，二进制直接报错 |
| 单元类型 | 单一 block，`elementType=4`（tet4） | 线性四面体，4 节点/单元 |
| 边界 | 37/38 个文件带 `$Surface` | 优先用它渲染/导出；缺失时自行提取 |
| 规模 | 4 节点（tet.msh）到 22 万节点（Armadillo219K） | 解析器要能应付大文件，演示选小网格 |

### 2.2 解析器的工程要点（`cpp/include/io/msh_loader.h`）

- **tag → 索引映射**：MshIO 的 tag 不保证从 0 连续（一般从 1 开始），解析时用
  `unordered_map<long,int>` 建立映射，单元引用不存在的 tag 立即报错。
- **段结构状态机**：逐行读，遇 `$Xxx` 进入对应解析，遇未知段则跳到配对的 `$EndXxx`；
  行尾 `\r` 统一剥掉（兼容 Windows 换行）。
- **多 block 容错**：按 MshIO 规范处理 `numEntityBlocks > 1` 的情况；
  `parametric != 0` 的节点块（带额外参数坐标）直接报「不支持」而不是默默错位。
- **单元类型分派**：`4 → 四面体`；`2 → 三角面`（有些文件把表面显式写成单元）；
  其余类型跳过（一行一个单元，跳过安全）。
- **表面来源优先级**（这是渲染/导出的关键）：
  1. 显式三角面单元（elementType=2）——最权威；
  2. `$Surface` 段（IPC 自定义）；
  3. **边界提取**：统计所有四面体的 4 个面，只出现一次的即边界。
     - 面用「排序后的三节点索引」哈希计数；计数为 1 时保留其在四面体中的原始走向，
       保证相邻面走向一致（有利于法线方向）。
     - 提取结果排序后输出，保证不同编译器/运行下顺序稳定（OBJ 与测试可复现）。
     - 验证：`cube_no_surface.msh` 提取出 12 个面，每条无向边恰好被 2 个面共享 → 闭合。
- **错误信息带文件名+行号**，方便定位坏文件：`path:line: 期望 $EndElements，实际: ...`。

> 为什么不直接用 MshIO？IPC 的 CMake 会从网络拉取 MshIO 等一堆依赖；fiction 是教学型项目，
> 自研解析器不到 400 行、无依赖、可控、可测，更符合「教学 + 最小实现」的定位。

---

## 3. 从四面体网格到 mass-spring 模型

`cpp/include/scene/mass_spring_builder.h` 做四件事。

### 3.1 唯一边提取（弹簧拓扑）

每个四面体贡献 6 条边：

```
tet (0,1,2,3) → (0,1) (0,2) (0,3) (1,2) (1,3) (2,3)
```

相邻四面体共享边会重复出现，用 `unordered_set<uint64_t>` 去重，key 取
`(min_node << 32) | max_node`。得到 `edges: [E,2]`。

> 小知识：cube.msh 的 6 个四面体只覆盖立方体 28 条几何边中的 **19** 条。
> 弹簧网络取决于四面体剖分方式，这也是 mass-spring 与真实连续体模型差异的来源之一。

### 3.2 集中质量（lumped mass）

真实连续体的质量分布是 ρ·dV；离散化后常用 **lumped mass**（把质量集中到节点）：

```
m_i = ρ · Σ_{tet ∋ i} V_tet / 4
```

- 每个四面体的体积质量 ρ·V 平均分给它的 4 个顶点，节点累加所有相邻四面体。
- 四面体体积（对节点定向不敏感，取绝对值）：

```
V = | (b − a) · ((c − a) × (d − a)) | / 6
```

- 零质量保护：`m_i < 1e-12 · m_total / N` 时抬到下限，避免质量矩阵奇异。
- 默认 ρ = 1000 kg/m³，对齐 IPC 的默认 `density`。

### 3.3 静止长度与刚度

- 静止长度平方：`l2_e = |x0_j − x0_i|²`，零长度保护 `l2 = max(l2, 1e-24)`（防止除零）。
- 刚度两种模式：
  - 显式：`--stiffness k`，所有边统一 k（每边线刚度，单位 N/m）。
  - 自动：`k = E_young · l_avg`，默认 `E_young = 1e5 Pa`（IPC 默认 Young's modulus）。
- 量纲直觉：连续体刚度 k ≈ E·A/l，而边长 l 的网格上 A ≈ l²，于是 k ≈ E·l。
  实测 cube（l_avg≈1.17）auto k≈1.17e5，与 Python 参考实现 1m 正方形取 k=2e4 的量级一致。
- 为什么不能照搬 IPC 的 E、ν？**mass-spring 没有泊松比**，弹簧只沿边传力，
  横向刚度/体积守恒行为与 Neo-Hookean 完全不同。k = E·l_avg 只是「量级对齐」，
  目的是让不同模型的手感差不多，而不是物理等价（见 4.5 节的边界讨论）。

### 3.4 包围盒与坐标盒选点

- `boundingBox(x)` → `[xmin, ymin, zmin, xmax, ymax, zmax]`，供相机取景与选点。
- `nodesInBox(x, box)`：闭区间筛选盒内节点，是 `--pin-box`（固定）与 `--force-box`（外力区域）的基础。
  这直接照搬 IPC 场景脚本里 `DBC/NBC` 的「坐标盒」语义（见 `IPC/input/tutorialExamples/BC/2cubesFall_NBC.txt`）。

---

## 4. 物理内核：IP 框架与三个能量

### 4.1 隐式欧拉 + 增量势能（IP）

一步时间推进从 `x_n, v_n` 开始：

```
x̃ = x_n + h·v_n                      （predictive / 惯性目标位置）
E_IP(x) = ½ Σ_i m_i |x_i − x̃_i|²      （Inertia，内置在 TimeIntegrator）
        + h² · Σ_k E_k(x)             （所有 IEnergy，按 hessianWeight=h² 加权）
        + DBC 罚能量
x_{n+1} = argmin_x E_IP(x)            （projected Newton 迭代）
v_{n+1} = (x_{n+1} − x_n) / h
```

要点：**能量函数统一写成「对位置的能量」，乘 h² 后一起无量纲化**。
这来自 IPC 的增量势能形式，也是为什么每个 `IEnergy` 都有一个 `hessianWeight(h)` 钩子
（默认 h²，重力/弹簧/区域外力都用默认；ZLock 用 1）。

### 4.2 边弹簧能量（弹力）

`cpp/include/energies/SpringEnergy.h`：

```
E_e = ½ k_e l_e² ( |d_e|²/l_e² − 1 )²        d_e = x_i − x_j，l_e = 静止长度
```

- 它惩罚「伸长率的平方」：|d|²/l² → 1 时能量为零。取平方比取一次方更光滑（一次方在 0 处不可导）。
- 梯度（对 i 节点为正、对 j 节点为负）：

```
g_diff = 2 k_e ( |d|²/l² − 1 ) d
```

- 局部 Hessian（差值坐标下的 3×3 块）：

```
H_diff = (2 k_e / l²) · ( 2 d dᵀ + (|d|² − l²) I₃ )
```

  组装成 6×6 边块 `[[H,−H],[−H,H]]`，再 `makePSD` 投影后再进入全局稀疏矩阵。
- **PSD 投影的意义**：精确 Hessian 在弹簧被压缩（|d|<l）时非正定，
  特征值分析给出：

```
沿 d 方向:    λ_∥ = (2k/l²)(3|d|² − l²)     > 0 恒成立
垂直 d 方向:  λ_⊥ = (2k/l²)(|d|² − l²)      < 0 当 |d| < l
```

  投影牛顿法把负特征值截断为 0，保证牛顿方向是下降方向（这也是投影牛顿名称的由来）。
- 这也是 `mesh_check` 第 4 项 FD 校验要在「拉伸位形」（整体放大 1.1 倍）下做的原因：
  拉伸时精确 Hessian 本来就半正定，`makePSD` 是恒等变换，数值差分才能与解析式逐元素对上。

### 4.3 重力（体力）

`cpp/include/energies/Gravity.h`：

```
E_g = −Σ_i m_i · (g · x_i)            g = (0, 0, −9.81)
grad_i = −m_i g                        （常数）
Hess = 0
```

线性势能：能量对位置是一次的，梯度是常数力，Hessian 为零。
牛顿系统不因它变坏，因为 Inertia（对角质量矩阵）已经保证了正定性。

### 4.4 区域恒定外力（对标 IPC NBC）

`cpp/include/energies/RegionForce.h`：

```
E_f = −Σ_{i ∈ region} m_i · (a · x_i)      a = 恒定加速度
grad_i = −m_i a                            （区域外为 0）
Hess = 0，hessianWeight = h²
```

**与 IPC 的对应关系**（这是本次连接最核心的「理论对接」）：

IPC 的 Neumann 边界条件（`IPC/src/TimeStepper/Optimizer.cpp` 的 `energyVal -=` / `gradient -=` 两处）：

```cpp
// 能量值（Optimizer.cpp ~L3242）
energyVal -= dtSq * data.massMatrix.coeff(vI, vI) * data.V.row(vI).dot(NBC.force.transpose());
// 梯度（Optimizer.cpp ~L3453）
gradient.segment<dim>(vI*dim) -= dtSq * data.massMatrix.coeff(vI, vI) * NBC.force;
// 其中 V = x − x̂ 是本步位移，NBC.force 就是脚本里的加速度向量
```

把 `V = x − x̂` 代入、忽略与未知量 x 无关的常数项：

```
IPC:  −h² m a·(x − x̂)  ≡  −h² m a·x + const
本实现: h² · E_f = h² · (−m a·x)
```

梯度完全相同。所以 `--force-box xmin…zmax ax ay az` 就是 IPC 场景文件里
`NBC <box> <ax ay az>` 的等价物（时间窗 t_begin/t_end 暂未实现，见第 9 节）。

> 注意：外力作用在**质量加权**上（grad = −m_i a），等价于「给这些节点一个加速度」，
> 而不是「一个固定的牛顿力」。这与 IPC 一致。

### 4.5 边界讨论：为什么只能定性对照 IPC

- IPC 是 **P1 四面体有限元 + Neo-Hookean（或固定共旋）材料**：每个四面体有连续的变形梯度
  F = ∂x/∂X，能量 ψ(F) 是材料模型，单元刚度块是 12×12。
- fiction 的 mass-spring 是**沿边的离散弹簧**：没有连续极限保证，刚度/变形与网格剖分强相关，
  也无法表达泊松效应。
- 因此本阶段的产物适合「看效果、做实验、演示管线」，**不能**把位移/应力与 IPC 数值直接比较。
- 定量对照的正路在 Phase 2：新增 12×12 局部块的四面体 FEM 能量（`energy.h` 的注释里已经预留了
  「四面体 12×12」的约定），质量矩阵沿用本节的集中质量即可。

---

## 5. 场景、边界与命令行（`cpp/include/scene/ipc_scene.h`）

### 5.1 pin-box → DBC（零速固定）

命令行给出坐标盒 `--pin-box xmin ymin zmin xmax ymax zmax` 后：

1. `nodesInBox` 选出盒内节点并**去重排序**（重复加罚会把刚度放大，必须去重）；
2. 构造 `DBCInfo`：`indices`、`targetVelocities = 0`、`limits = x0`、`isDBC[i]=1`。

TimeIntegrator 的目标位置逻辑（`src/time_integrator.cpp` 主循环）：

```cpp
if ((lim - x_n).dot(vel) > 0) target = x_n + h·vel;   // 速度朝限位方向 → 继续驱动
else                          target = x_n;           // 否则保持当前位置
```

当 `vel = 0` 时恒取 `target = x_n`，即**每个时间步都要求该节点速度为零**——这就是标准的
Dirichlet 速度约束。收敛判据满足后（`|x − target|/h < tol`），该自由度会被「消去」
（Hessian 行列置单位阵、梯度置零），因此实测 pin 节点位移**精确为 0**
（`mesh_check` 检验 6：`maxPinDisp = 0.000000`）。

### 5.2 force-box → RegionForce

`--force-box … ax ay az` 可重复多次；每次选点得到一个节点集合，生成一个 `RegionForceEnergy`。
支持多个盒子 → 多个外力能量叠加，例如「左端向下拉、右端向上拉」的剪切场景。

### 5.3 其余 CLI

| 参数 | 含义 |
|---|---|
| `--mesh` | 必填，.msh 路径 |
| `--density` / `--stiffness` / `--youngs` | 质量与刚度（见 3.2/3.3） |
| `--dt` / `--steps` / `--tol` / `--max-iter` | 时间与牛顿收敛参数 |
| `--gravity gx gy gz` / `--no-gravity` | 重力开关 |
| `--init-vel vx vy vz` | 全体节点初始速度（做纯弹簧振荡演示用） |
| `--pin-box` / `--force-box` | 边界与区域外力（可重复） |
| `--export <dir>` / `--export-every n` | OBJ 帧导出；采样间隔同时被视频复用 |
| `--video` / `--video-fps` / `--video-size` / `--video-static` / `--video-keep-frames` | 视频导出（见 6.3） |
| `--verbose` | 打印 Newton 每步迭代日志 |

参数解析是纯手写、无第三方库（同一风格贯穿全项目），错误信息明确到具体 flag。

---

## 6. 可视化：面片渲染与 SSH 视频导出

### 6.1 实时窗口（`ipc_mesh_visual`）

- `viewer.h/.cpp` 新增 `drawTriangles`：固定管线立即模式，逐面法线 + 定向光
  `shade = 0.35 + 0.65·|n·L|`（取绝对值实现双面可见，模型翻过来也亮）。
- 相机自动按包围盒取景；pin 节点画红点、外力区域画橙点。
- 需要 GLFW/OpenGL 与显示环境；无 DISPLAY 时优雅报错退出。

### 6.2 为什么还要 CPU 软光栅化

SSH 远程没有 X/Wayland，GLFW 窗口开不出来。但要「可视化查看效果」，
最稳的路线是**把渲染也做成纯计算**：自己算投影、光栅化、写图片，再用 ffmpeg 编码。
这样唯一的外部依赖只剩 ffmpeg（连它也只在最后编码时用到）。

`cpp/include/io/rasterizer.h` 的实现要点：

1. **相机**：`SoftCamera {target, radius, yaw, pitch, fov}` 与 viewer 的轨道相机约定一致；
   共用 `camera.h` 的 `perspective` / `lookAt`（从 viewer.cpp 抽出，避免两份数学）。
2. **投影与裁剪**：`clip = P·V·[x,y,z,1]ᵀ`；`w ≤ 0`（相机后）的顶点直接丢弃整个三角形。
   相机按模型框自动取景，不会出现近裁剪面切割问题。
3. **光栅化**：屏幕空间有向边函数

```
area = edge(p0,p1,p2)
w0 = edge(p1,p2,P)，w1 = edge(p2,p0,P)，w2 = edge(p0,p1,P)
P 在三角形内 ⇔ 三个 w 与 area 同号（利用正负号同时处理两种绕向）
重心插值深度：z = (w0·z0 + w1·z1 + w2·z2) / area
```

4. **深度插值为什么正确**：透视除法之后，`z_ndc` 是屏幕坐标的**仿射函数**，
   因此在屏幕空间做线性（重心）插值得到的就是正确的深度值——这就是经典 z-buffer 可以
   只插值 `z/w` 的原因（透视校正只对颜色/UV 等属性必要，纯色平涂不需要）。
5. **超采样抗锯齿**：默认 2×（1920×1440 渲染 → 960×720 输出），盒式降采样。
6. **地面网格**：3D 直线投影后用 DDA 画线、逐像素做深度测试；2×2 覆盖保证降采样后仍可见。
7. **节点高亮**：pin/外力节点画圆点，同样带深度测试。

### 6.3 图像与视频编码

- `io/image_writer.h`：**PNG**（zlib）：8 位真彩，逐行 filter=0，`compress2` 压缩 IDAT，
  CRC32 校验块；无 zlib 时回退 **PPM (P6)**（ffmpeg 同样能读）。
- `io/video_encoder.h`：调用 ffmpeg

```
ffmpeg -y -framerate <fps> -i "<dir>/frame_%06d.png" \
       -c:v libx264 -pix_fmt yuv420p -crf 20 -movflags +faststart <out.mp4>
```

  首选 libx264（兼容性最好），失败自动回退 mpeg4；ffmpeg 不存在时保留帧目录并打印手工命令。
- 帧命名 `frame_000000.png` 与 image2 的 `%06d` 对齐；yuv420p 要求宽高为偶数，程序里自动向下取偶。

### 6.4 时长公式（10 秒视频怎么来）

```
帧数   = steps / export-every + 1
时长   = 帧数 / fps
实时 1x 条件： export-every × dt × fps = 1
```

本次 4 条 10 秒演示视频取：`steps=1000, dt=0.01, export-every=4, fps=25`
→ 251 帧 → 10.04s，且恰为 1 倍速。相机以 18°/s（仿真时间）环绕，10 秒转半圈，
即使模型进入静态平衡也能从旋转中看清三维形变。

---

## 7. 如何证明「连接正确」：验证方法学

`cpp/examples/mesh_check.cpp` 用 6 组独立检验覆盖每一层（全部 PASS）：

| # | 检验 | 验证目标 | 关键思想 |
|---|---|---|---|
| 1 | 解析器 basics | cube.msh: 8 节点/6 四面体/12 面/体积=1.0/bbox/质量=1000 | 拿**解析解已知**的最小网格当标尺 |
| 2 | 无 `$Surface` 回退 + 错误路径 | 边界提取 12 面且每条边被 2 面共享；二进制/2.2 版本/缺段报错 | 闭合性检验等价于「提取的确实是边界」 |
| 3 | 自由落体刚性平移 | 弹簧 + 重力、无 DBC，单步位移精确 = g·h² | **刚性平移是弹簧能量的零空间**：所有边长度不变 → 弹簧力为零；牛顿方程的解可解析写出 |
| 4 | 弹簧 FD | 拉伸位形下 FD 梯度/Hessian vs 解析式 | 位形选取保证精确 Hessian 已 PSD（见 4.2） |
| 5 | 区域外力 FD | 梯度 = −m·a、区域外为 0、Hessian 为空 | 线性势能的解析解直接可比 |
| 6 | pin + force / pin + gravity 仿真 | pin 位移精确 0；外力方向正确；结果有限；重力下垂方向正确 | 端到端行为检验（方向 + 数值健康） |

补充的**覆盖性实验**（一次性跑过，未固化为测试）：

- IPC `input/tetMeshes/` 全部 **38/38** 个 `.msh` 都能加载并建成弹簧网（大步长跑 0 步，只验证解析+装配）。
- 性能实测（本机）：bar-186 约 4–12 ms/步、sphere1K 约 75 ms/步（tol=1e-4）、bunny 约 190 ms/步（tol=1e-2）。
- 视频侧：ffprobe 校验 h264/960×720/251 帧/10.04s；抽帧统计非空（均值/方差非零），防止「黑屏视频」。

FD 校验的通用方法（与 `module_test/gradient_hessian.cpp`、`barrier_check.cpp` 一致）：

```
梯度:   g_i ≈ [E(x+εe_i) − E(x−εe_i)] / (2ε)
海森:   H(:,k) ≈ [g(x+εe_k) − g(x−εe_k)] / (2ε)     （对梯度再做中心差分）
判定:   误差 / max(|解析值|) < 1e-6
```

---

## 8. 改动清单（新增/修改，以及为什么）

### 8.1 新增文件

| 文件 | 层 | 职责 | 设计理由 |
|---|---|---|---|
| `cpp/include/camera.h` | 数学 | `perspective` / `lookAt` | viewer 与软光栅化共用，消灭重复实现 |
| `cpp/include/io/msh_loader.h` | IO | MshIO 4.1 ASCII 解析 + 边界提取 | 无第三方依赖；错误信息带行号 |
| `cpp/include/io/obj_writer.h` | IO | 表面 OBJ 导出 | 最通用的交换格式，MeshLab/Blender 直开 |
| `cpp/include/io/rasterizer.h` | IO | CPU 软光栅化 | SSH 无显示环境可渲染；z-buffer + 2× 超采样 |
| `cpp/include/io/image_writer.h` | IO | PNG（zlib）/PPM | PNG 方便预览，PPM 作为零依赖回退 |
| `cpp/include/io/video_encoder.h` | IO | ffmpeg 编码封装 | 视频是远程查看的最优载体 |
| `cpp/include/scene/mass_spring_builder.h` | 场景 | 边/质量/l2/刚度/选点 | 几何→物理的唯一转换点，可独立测试 |
| `cpp/include/scene/ipc_scene.h` | 场景 | CLI + 场景装配 + 统计 | demo 与 visual 共享，参数与物理不分叉 |
| `cpp/include/energies/RegionForce.h` | 物理 | 区域恒定外力 | IPC NBC 的等价能量，~80 行 |
| `cpp/examples/ipc_mesh_demo.cpp` | 入口 | 控制台仿真 + OBJ + 视频 | 无显示环境主力工具，退出码供 CI |
| `cpp/examples/ipc_mesh_visual.cpp` | 入口 | 交互式窗口 | 有显示环境时的调试/演示 |
| `cpp/examples/mesh_check.cpp` | 测试 | 6 组数值验证 | 分层验证，CI 可跑 |
| `cpp/tests/data/{cube.msh, cube_no_surface.msh, bar-6.msh}` | 测试夹具 | 1KB 级样本 | 测试不依赖仓库外的 IPC 目录 |

### 8.2 修改文件

| 文件 | 改动 | 理由 |
|---|---|---|
| `cpp/src/viewer.cpp` / `cpp/include/viewer.h` | 抽出相机数学到 `camera.h`；新增 `drawTriangles` | 复用 + IPC 网格需要面片渲染 |
| `cpp/CMakeLists.txt` | 新 target：`ipc_mesh_demo` / `ipc_mesh_visual` / `mesh_check`；zlib 探测；新增 ctest | 与既有 `simulator/barrier_check` 的组织方式一致 |
| `.gitignore` | 忽略 `output/` | 帧/视频/OBJ 均为生成物 |
| `AGENTS.md` | 布局、命令、视频配方、注意事项 | 让后续 agent/人类快速找到入口 |
| `cpp/include/energies/README.md` | 补全「现有能量」清单 | 文档与代码同步 |

### 8.3 明确**没有**做的事

- 没改装 `TimeIntegrator` / `IEnergy` 接口 / CHOLMOD 链路——新增全部通过扩展点完成；
- 没改 Python 参考实现（本次是 IO + 场景装配 + 工具链，不属于「算法需双端镜像」的范畴）；
- 没引入 IPC 的任何 C++ 依赖，也没构建 IPC 本体（对照实验留待 Phase 2）。

---

## 9. 已知边界与后续路线

### 9.1 当前边界（务必知晓）

1. **物理定性**：mass-spring 与 IPC 的 Neo-Hookean FEM 不等价；数量级不同、无泊松效应。
2. **无接触/摩擦**：地面网格仅视觉参考，模型会穿地；接触目前只有平面 Barrier 能量（未接入本 demo）。
3. **外力无时间窗**：IPC 的 `NBC … [t_begin] [t_end]` 支持时段开关，本实现是「全程常开」。
4. **模型规模**：>7K 节点（sphere5K、Armadillo 系列）每步几百 ms 起步，演示建议 ≤2K 节点。
5. **仅 tet4**：`elementType` 非 4 的网格不能用于物理（tri 面网格也暂未支持）。
6. **视频离屏渲染是 CPU 逐三角形**：大网格 + 高分辨率的渲染时间会显著增长。

### 9.2 Phase 2 路线（按建议优先级）

| 优先级 | 任务 | 落点 |
|---|---|---|
| 高 | 四面体 FEM 能量（Neo-Hookean / 固定共旋，12×12 局部块） | `cpp/include/energies/`，`energy.h` 已预留约定 |
| 高 | 与 IPC 同场景定量对照（需构建 IPC 跑基准） | 脚本 + 文档，比较位移/能量曲线 |
| 中 | `--force-box` 增加时间窗（t_begin/t_end） | `RegionForceEnergy` 加时间上下文或 demo 层开关 |
| 中 | 解析 IPC 场景 `.txt` 子集（shapes/material/DBC/NBC/ground） | `include/scene/` 新增解析器 |
| 低 | tri 面网格（`input/triMeshes/*.obj`）支持 | loader + 面弹簧拓扑 |
| 低 | 接触/摩擦（把 BarrierEnergy 接入 demo） | 复用 `BarrierEnergy.h` |

---

## 10. 复现命令速查

```bash
# 构建
cmake -S cpp -B cpp/build -DCMAKE_BUILD_TYPE=Release && cmake --build cpp/build -j

# 数值验证（解析器 + 物理 + FD）
./cpp/build/mesh_check
ctest --test-dir cpp/build

# 控制台仿真：左端固定 + 右端外力
./cpp/build/ipc_mesh_demo --mesh /home/steven/mjx/Archive/IPC/input/tetMeshes/bar-186.msh \
    --pin-box -0.1 -0.1 -0.1 0.6 0.6 1.1 \
    --force-box 9.4 -0.1 -0.1 10.1 0.6 1.1 0 0 40 \
    --no-gravity --youngs 1e6 --steps 1000 --export-every 4 --tol 1e-6

# 10 秒视频（1 倍速实时）
./cpp/build/ipc_mesh_demo --mesh /home/steven/mjx/Archive/IPC/input/tetMeshes/bar-186.msh \
    --pin-box -0.1 -0.1 -0.1 0.6 0.6 1.1 --gravity 0 0 -9.81 --youngs 1e6 \
    --steps 1000 --export-every 4 --tol 1e-6 --video output/videos/ipc_models/gravity_bar.mp4

# 全量模型解析兼容性（38/38）
for f in /home/steven/mjx/Archive/IPC/input/tetMeshes/*.msh; do
    ./cpp/build/ipc_mesh_demo --mesh "$f" --no-gravity --steps 0 > /dev/null || echo "FAIL: $f"
done
```

---

## 附：一句话总结

> **IPC 提供「几何」，fiction 提供「物理」，中间的适配层做三件事：
> 解析四面体网格 → 沿边建质量-弹簧系统 → 用与 IPC 同构的能量（重力、弹簧、区域外力）
> 交给既有的 IP 牛顿求解器；再用 CPU 软光栅化 + ffmpeg 把结果变成 SSH 也能看的视频。**
> 全程没有引入 IPC 的依赖，也没有改动 fiction 的求解框架——这正是好接口的价值：
> 新能力是「加文件」，不是「改内核」。
