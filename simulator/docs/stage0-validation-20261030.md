# S0 / T0 阶段零验收记录

日期：2026-10-30。状态：**代码与本轮数值检查完成，停下等待用户理论与工程审查。** 未进入 S1/S2、静态求解、结构化杆生成、大规模、多体、接触或 GPU 实验。

## 1. 授权、范围与改动

用户授权开展首阶段，要求粒度细、适当并发使用低端子 agent，并在阶段零两个编号完成后停止。本轮将两个编号明确解释为 S0、T0；采用三个原生 gpt-6-luna 子 agent，分别负责材料测试、动态测试、独立只读复核，主会话集成与验收。

| 文件 | 本轮改动 |
|---|---|
| [material_validation.cpp](../cpp/examples/material_validation.cpp) | 新增 S0 可执行检查，失败返回非零，输出 FD CSV 与其他指标日志 |
| [dynamic_reference_check.cpp](../cpp/examples/dynamic_reference_check.cpp) | 新增 T0 测试专用线性二次能量、独立离散参考、连续解和逐步 CSV |
| [mesh_check.cpp](../cpp/examples/mesh_check.cpp) | 仅修正旧旋转测试：固定参考，只旋转当前网格，同时检查能量和梯度 |
| [CMakeLists.txt](../cpp/CMakeLists.txt) | 注册两个可执行文件和两个 CTest 条目 |

生产 NH、积分器、质量 builder 和默认求解算法未改动。旧工作区原有修改和未跟踪文件保留，未提交 Git。

## 2. S0 验收

- 静止、当前位形平移、固定参考下刚体旋转，以及变形后叠加旋转的能量与梯度协变检查通过。
- 拉伸、剪切、压缩三态，每态差分步长为 1e-4、1e-5、1e-6、1e-7；要求至少两档相邻步长同时达到归一化梯度误差 <1e-6、原始 Hessian 误差 <1e-4。实际 12/12 档均通过，最大值分别约 **4.61e-8、1.66e-8**。
- PSD 单独验收：三态与独立谱裁剪重构一致，并确认剪切/压缩原始 Hessian 确有负特征值；本次最小值约 -2.47e3 / -4.64e4。投影后有浮点舍入量级负值，按相对容差通过。谱裁剪重构使用 Eigen 特征分解，核对实现路径。
- 固定泊松比时，E 翻倍使能量、梯度范数、原始 Hessian 范数均翻倍。均匀体积加载下，nu=0.20/0.40 的能量导数与解析参考一致，随 nu 增大。
- 两单元加和、平移零模态、总体积 1/3、总质量 400，以及节点质量 [100,100,100,50,50] 均通过。

两单元加和检查装配一致性。S0 的材料状态取正体积、远离截断区域的拉伸、剪切和压缩。

## 3. T0 验收

测试专用势能为线性二次势；未将生产代码的非线性边弹簧误当线性弹簧。刚度 20，求解停止容差 1e-10。离散参考从初态独立累计，不逐步重置为数值解。

| 项目 | 设置与结果 |
|---|---|
| 离散振子 | 质量 1/2/5 × 时间步 0.01/0.02/0.04，共 9 组，各 120 步，三分量位置/速度均比较 |
| 累计位置误差 | 最大约 8.94e-15，阈值 1e-10（位置单位） |
| 累计速度误差 | 最大约 4.02e-14，阈值 1e-10（速度单位） |
| 纯惯性/恒定重力 | 各 120 步，与精确离散更新一致 |
| 连续谐振子 | 质量 2、物理时间 2 s；步长 0.04/0.02/0.01，无量纲状态误差约 0.3240/0.1790/0.09426 |
| 数值耗散 | 同三档末态/初态能量比约 0.4522/0.6709/0.8188，随时间步减小趋近 1；不当作真实材料阻尼 |

共保存 5010 行逐分量轨迹数据。T0 使用无 DBC 小系统和上表输入；当时积分器的失败报告接口待补 StepReport。

## 4. 构建、回归与复现

独立构建目录：`/tmp/fiction-stage0-20261030-build`。Release，GCC 11.4.0、CMake 3.28.0、Eigen 3.4.0，配置确认 CHOLMOD 启用。未验证全新安装环境或本轮 LDLT-only 构建。

```bash
cd /home/steven/mjx/fiction/simulator
cmake -S cpp -B /tmp/fiction-stage0-20261030-build -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/fiction-stage0-20261030-build -j 3
ctest --test-dir /tmp/fiction-stage0-20261030-build --output-on-failure
/tmp/fiction-stage0-20261030-build/material_validation output/validation/stage0-20261030-01/s0-fd.csv
/tmp/fiction-stage0-20261030-build/dynamic_reference_check --csv output/validation/stage0-20261030-01/t0-trajectory.csv
```

CTest **8/8 通过**，其中两个 visual 条目为占位/构建提示，两个新增条目实际执行 S0/T0。旧 simulator 和 mesh_check 另保存完整 stdout，日志无 FAIL 或 Warning；当时 simulator 固定返回 0，验收读取完整日志。

CTest 本轮小检查汇总墙钟为 0.08 s，编译另计。材料和动态 CSV 为最终源码构建后重新执行生成。

## 5. 证据与停点

原始目录：[stage0-20261030-01](../output/validation/stage0-20261030-01)。含 configure/build/CTest 日志、S0/T0 和旧测试日志、CSV、机器可读指标、before/after diff 与状态、相关源码压缩快照、源码和可执行文件 SHA-256。

基础 HEAD 为 `b3745128119eb5def111159dbd3310a9e31fcc6f`；已有核心实现多数未跟踪，复现使用本轮 `source.tar.gz` 和配置记录。before.patch 保存原已跟踪文件的既存差异，未跟踪文件由源码快照保存。

当前停点为 S0/T0 用户审查。后续所需 StepReport、失败语义、静态模式、杆网格生成器及 S1 均未在本轮实施；没有新展示素材、接触、GPU 或大规模结果。
