# energies/ — 能量扩展目录（占位说明）

本目录用于存放你自定义的能量实现，每个能量继承 `include/energy.h:225 IEnergy`。

## 已有接口

- `include/energy.h:1` — 抽象基类 `IEnergy`，含 6 步添加指南与完整伪代码示例 `energy.h:134-190`
- `include/utils.h:1` — `makePSD` 工具，用于非 SPD 局部 Hessian 投影
- `include/time_integrator.h:1` — `TimeIntegrator` 内置 Inertia，其余能量通过 `addEnergy` 注册

## 如何添加一个新能量

> 详见 `include/energy.h:36-196` 顶部大段注释，此处精简版：

1. 新建头文件 `include/energies/my_energy.h`：
```cpp
#pragma once
#include "energy.h"
#include "utils.h"
namespace fiction {
class MyEnergy final : public IEnergy {
 public:
  const char* name() const override { return "MyEnergy"; }
  double val(const Positions& x) const override;
  void addGrad(const Positions& x, Positions& grad) const override;
  void addHess(const Positions& x, std::vector<Triplet>& trips) const override;
  // 可选：double initStepSize(...) const override; // Barrier CCD
  // 可选：double hessianWeight(double h) const override { return 1.0; } // 若不应缩放
};
}
```

2. 新建 `src/energies/my_energy.cpp` 实现上述 3 函数（注意 `addGrad` 用 `+=`，`addHess` 用 `3*node+d` 固定索引，必要时 `makePSD`/`makePSD_dynamic`）。

3. 在 `CMakeLists.txt` 中将 `src/energies/my_energy.cpp` 加入 `fiction_core` 或新 `add_library`，并在 `examples/simulator.cpp` 中：
```cpp
integrator.addEnergy(std::make_shared<MyEnergy>(...));
```

4. 在 `src/time_integrator.cpp` 侧无需改动，`TimeIntegrator` 会自动在 `computeIPVal/Grad/Hess` 循环中调用并按 `hessianWeight(h)` 加权（默认 `h*h`）。

## 现有能量

- `Gravity.h` — 重力（均匀重力场），`Energy.h` 最小骨架示例。
- `SpringEnergy.h` — MassSpring 边弹簧（6×6 局部块 + makePSD），复制后改名适配新单元。
- `BarrierEnergy.h` — Barrier 接触能量（对标 `Python/BarrierEnergy.py`：floor 固定平面 + ceiling 节点锚定平面，含 Hessian 与 initStepSize 数值 CCD），验证见 `examples/barrier_check.cpp`。
- `Gravitation.h` — 万有引力（固定点质量吸引子，Hessian 投影 PSD），演示见 `examples/simulator.cpp`。
- `RegionForce.h` — 坐标盒区域恒定外力（线性势能，对标 IPC Neumann BC），IPC 模型示例见 `examples/ipc_mesh_demo.cpp`。
- `NeoHookeanTetEnergy.h` — Neo-Hookean 四面体 FEM（P1、12×12 局部块、λ/μ 由 E/ν 换算，Elastogen cantilever 复现用）；`--material nh` 启用，FD 校验见 `examples/mesh_check.cpp` 检验 5.5。
- `zlock_energy.h` — 2D 退化 z 锁。
- 其余能量（`FrictionEnergy` 等）请按需自行实现，参考 `Python_slovement/*.py` 对应逻辑与维度约定。

## Inertia 说明

- `InertiaEnergy` 已作为 `TimeIntegrator` 内置项实现（`src/time_integrator.cpp`），不在此目录注册，权重固定为 1，不乘 `h*h`，对标 `time_integrator.py:59`。纯 3D 固定 `[N,3]`。
- 若需替换惯性模型，请直接修改 `TimeIntegrator::inertiaVal/AddGrad/AddHess` 三函数。

## 纯 3D 要点（2D 适配）

- `Positions` 固定 `[N,3]`（RowMajor），`kDim=3`，自由度索引固定 `3*node+d`
- 能量实现直接用 `Vec3` 与 `x.row(i)`，局部块 6×6（边）/9×9（面）/12×12（体）
- 2D 算例：`upgrade2DTo3D` 补零 + `ZLockEnergy`（`zlock_energy.h`）锁 z，`x/y` 运动不受影响
- DBC 使用 `Vec3`（target/limit），不再有 `dbc.dim` 字段
- `utils::makePSD` 提供 3×3/6×6/9×9/12×12 与 `makePSD_dynamic` 通用版（2×2/4×4 已废弃，仅留内联转发）
