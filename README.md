# fiction-simulation-experiments

一个以 **C++20 三维超弹性固体仿真器**为主的学习与实验项目。从四面体材料能量、全局装配到隐式时间推进，逐层检查实现，并通过拉伸、弯曲、扭转、圆环与独立多体场景观察离散误差和计算成本。

- [中英文项目介绍与实验图示](https://maojunxiang.divinechime64.chatgpt.site/zh/blog/fiction-solid-simulator/)
- [实验代码与结果导航](EXPERIMENTS.md)
- [仓库检查与已知限制](docs/repository-audit-20261105.md)

## 实现范围

主工程在 [`simulator/cpp/`](simulator/cpp/)，包含：

- Neo-Hookean 四面体有限元：材料能量、梯度、Hessian、单元 Hessian 半正定投影和全局稀疏装配；另保留质量弹簧模型供对照。
- 隐式 Euler、Projected Newton、回溯线搜索、指定位置边界条件，以及迭代次数、收敛状态和物理力残差诊断。
- Eigen 稀疏求解；可选 CHOLMOD 后端。可选 GLFW/OpenGL 查看器，以及无显示环境的 OBJ、图像和视频导出。

**三维网格样例不处理物体间碰撞、自接触或摩擦。** 画面中的地面网格只是视觉参考。主工程虽然有平面 barrier 的单元检查，但它未接入网格样例；隔离的 S6b 接触探索也没有通过物理验收。项目不是通用有限元软件，运行成功或 CTest 通过不代表任意场景的物理精度已验证。

## 快速开始

以下命令均从仓库根目录执行。需要支持 C++20 的编译器、CMake 3.20+、Eigen 3.4+。Linux 上可通过系统包管理器安装；Ubuntu 示例：

```bash
sudo apt-get install build-essential cmake libeigen3-dev
```

最小构建不需要 Python、GPU、图形环境或外部网格：

```bash
cmake -S simulator/cpp -B build -DCMAKE_BUILD_TYPE=Release \
  -DFICTION_USE_VIEWER=OFF -DFICTION_USE_CHOLMOD=OFF
cmake --build build --parallel 3
ctest --test-dir build --output-on-failure
```

### 运行一个三维固体样例

使用仓库自带的单位立方体网格，固定 x=0 端面，在重力下推进 10 步并导出 OBJ：

```bash
./build/ipc_mesh_demo \
  --mesh simulator/cpp/tests/data/cube.msh \
  --material nh --youngs 1e6 \
  --pin-box -0.01 -0.01 -0.01 0.01 1.01 1.01 \
  --dt 0.005 --steps 10 --tol 1e-6 \
  --export output/quickstart --export-every 5
```

正常完成会返回退出码 0，并在 `output/quickstart/` 写出 `0.obj`、`5.obj`、`10.obj`。该小网格只用于确认构建、材料、约束、时间推进与输出链条，不用于定量弯曲精度验证。Newton 未收敛或状态出现非有限值时，控制台样例会报告失败步骤并以非零退出码停止。

**`--material` 默认是 `spring`；运行超弹性有限元必须显式指定 `--material nh`。** 全部参数可通过 `./build/ipc_mesh_demo --help` 查看。

### 可选功能

| 功能 | 依赖与设置 |
|---|---|
| CHOLMOD 稀疏求解 | 安装 SuiteSparse（Ubuntu：`libsuitesparse-dev`），重新配置 `-DFICTION_USE_CHOLMOD=ON`；未找到时回退 Eigen |
| 实时查看器 | GLFW3/OpenGL（Ubuntu：`libglfw3-dev libgl1-mesa-dev`），配置 `-DFICTION_USE_VIEWER=ON`；需可用显示环境 |
| 图像与视频 | zlib 用于 PNG（否则输出 PPM），ffmpeg 用于 MP4；控制台采用 CPU 软光栅化，无需显示环境 |

视频示例可在上述命令中加上 `--video output/quickstart.mp4 --video-static`。若缺少 ffmpeg 或编码失败，保留帧目录并打印手工编码提示；这不等于已经生成 MP4。

## 已开展的实验与证据范围

| 内容 | 已有记录及限制 |
|---|---|
| 材料与时间更新 | 材料有限差分、刚体不变性、装配和线性振子的离散参考检查 |
| 拉伸与静态弯曲 | 指定载荷实验、停止判据分析与空间网格加密；梁理论是近似参考 |
| 动态弯曲与扭转 | 时间步对照及驱动边界实验；原 S3 的 0–2 s、2% 轨迹对照未通过，后续时间收敛趋势不改写原记录 |
| 圆环与复杂网格 | 圆环运动、Armadillo 13K/28K/54K 各 10 步规模记录；不是严格的性能 benchmark 或长时稳定性证明 |
| 无接触多体 | 单体、四体和 16 体短段对照，检查互不耦合实例与独立单体轨迹的一致性 |
| 接触探索（隔离代码） | S6b C0 工具单元检查完成，C1 稳态间隙与反力门槛未通过，后续阶段搁置；不属于已验收主工程能力 |

各阶段的源码版本、报告入口、输入与复现边界见 [EXPERIMENTS.md](EXPERIMENTS.md)。

## 目录与复现边界

| 路径 | 用途 |
|---|---|
| [`simulator/cpp/`](simulator/cpp/) | 当前维护的 C++ 主工程、小网格输入与回归测试 |
| [`simulator/tools/`](simulator/tools/) | 结构化方杆网格生成与部分实验工具 |
| [`simulator/docs/`](simulator/docs/) | 保留原 run-id、参数与结论的历史报告 |
| [`experiments/workspaces/`](experiments/workspaces/) | 各阶段隔离版本；不等同于当前主工程 |
| [`simulator/Python/`](simulator/Python/) | 历史参考实现，仍使用旧导入名 `Python_slovement`，不是可直接运行的快速入口 |
| [`planning/`](planning/) | 历史验证协议与执行设计 |
| [`assets/ipc/`](assets/ipc/) | 随仓库分发的 torus、sphere1K 网格与来源记录 |

此仓库是**源码与部分小型输入的选集，不是完整实验数据归档**。主工程的构建、测试和上述样例可在检出后运行；历史 runner 则可能依赖原机器路径、独立构建、第三方 checkout 或未分发的大网格，不能直接把所有 runner 当作一键复现入口。完整原始轨迹、视频、部分图表和 Armadillo 网格未随仓库分发。历史报告中的本地结果链接因而可能无法在 GitHub 打开。

[`fiction.project.json`](fiction.project.json) 保留原工作区布局，其中 `results`、`figures`、`writing`、`consolidation_manifest` 是未随选集分发的历史路径或输出约定。不要将其视为仓库现有文件清单。

[`SOURCE_MANIFEST.json`](SOURCE_MANIFEST.json) 记录**首次导入时**各文件的来源和 SHA-256。导入包含原工程当时尚未提交的工作文件，不等同于原工程 commit `b3745128119eb5def111159dbd3310a9e31fcc6f` 的原样检出。清单保留不变，后续维护改动以本仓库 Git 历史为准。

## 来源与许可证

学习参考：[Physics-Based Simulation](https://phys-sim-book.github.io/) 及其 [solid-sim-tutorial](https://github.com/phys-sim-book/solid-sim-tutorial)。torus、sphere1K 来自 IPC，来源与哈希见 [`assets/ipc/sources.json`](assets/ipc/sources.json)；该文件也记录了未分发网格的历史信息。rod 测试网格由 ElastoGen 示例转换，说明见 [`simulator/docs/elastogen-rods-experiments.md`](simulator/docs/elastogen-rods-experiments.md)。

本仓库尚未指定项目源码许可证。[`assets/ipc/UPSTREAM_LICENSE`](assets/ipc/UPSTREAM_LICENSE) 仅为上游许可副本，不代表本仓库全部源码或所有第三方资产均采用该许可。
