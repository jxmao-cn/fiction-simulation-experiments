# AGENTS.md

Guidance for AI coding agents working in this repository.

## Project Overview

**fiction** is a mass-spring solids simulator: implicit Euler time integration where each step minimizes the Incremental Potential (IP) with projected Newton (PSD-projected Hessian, line search, DBC penalty + DOF elimination). 2D scenes are handled by zero-padding into 3D with a z-lock.

The codebase has two parallel implementations:

- `Python/` — the reference implementation (numpy/scipy/pygame), translated into C++ as the project grows.
- `cpp/` — the C++20 rewrite (`fiction_sim`). **The C++ side is the primary development target.**
- `module_test/` — standalone numeric sanity programs (not part of CMake).

Comments throughout the code are in **Chinese**. Keep this convention when editing existing files.

## Repository Layout

```
cpp/
  CMakeLists.txt          # project fiction_sim, C++20, CMake >= 3.20
  include/
    types.h               # [N,3] RowMajor Positions, Vec3, SparseMatrix (kDim=3)
    camera.h              # perspective/lookAt matrices (viewer 与软光栅化共用)
    energy.h              # IEnergy interface + "how to add a new energy" guide
    time_integrator.h     # implicit-Euler + Newton framework
    utils.h               # makePSD Hessian projection, inf-norm helpers
    viewer.h              # GLFW/OpenGL viewer (orbit camera, pause/reset)
    io/                   # msh_loader.h, obj_writer.h, rasterizer.h (CPU 软光栅化),
                          # image_writer.h (PNG/PPM), video_encoder.h (ffmpeg)
    scene/                # mass_spring_builder.h + ipc_scene.h (IPC 场景装配/CLI)
    energies/             # IEnergy implementations + README + example template
                          #   NeoHookeanTetEnergy.h — NH tet FEM（--material nh）
  src/
    time_integrator.cpp   # Newton loop, DBC, line search, CHOLMOD/LDLT solve
    utils.cpp
    viewer.cpp
  examples/
    simulator.cpp         # console validation demo, prints PASS/FAIL
    simulator_visual.cpp  # interactive GLFW window demo
    barrier_check.cpp     # BarrierEnergy FD/CCD/落地数值验证
    ipc_mesh_demo.cpp     # IPC .msh 控制台仿真（重力/弹簧/区域外力 + OBJ 导出）
    ipc_mesh_visual.cpp   # IPC .msh 交互式可视化（面片渲染）
    mesh_check.cpp        # msh 解析 + mass-spring/外力数值验证（CI）
  tests/data/             # 测试夹具：cube.msh, cube_no_surface.msh, bar-6.msh,
                          #   rod16x3x3.msh / rod_z16x3x3.msh（tools/bgeo_to_msh.py 由
                          #   Elastogen rod.bGeo 转换，论文 cantilever 16×3×3）
tools/
  bgeo_to_msh.py          # Elastogen 体素 bgeo → MshIO 4.1 tet（需 elastogen conda 环境）
docs/
  ipc-mesh-integration.md # IPC .msh ↔ fiction 连接的设计思路/理论/验证（学习文档）
  ipc-asset-inventory.md  # IPC/input 全部模型资产扫描存档（格式/规模/可用性）
  elastogen-rods-experiments.md # Elastogen cantilever（bend/twist NH）复现说明
Python/
  simulator.py            # entry point: square mesh + gravity + ground/ceiling DBC
  time_integrator.py      # step_forward, IP assembly, Newton + line search
  *Energy.py              # one file per energy (Inertia, MassSpring, Gravity, ...)
  square_mesh.py, utils.py, readme.md
module_test/
  gradient_hessian.cpp    # finite-difference check of grad/Hessian (compiled independently)
  gravity.cpp
```

## Build (C++)

Dependencies: Eigen3 >= 3.4, optional SuiteSparse/CHOLMOD and GLFW3 (both default ON, graceful fallback if missing).

```bash
cmake -S cpp -B cpp/build -DCMAKE_BUILD_TYPE=Release
cmake --build cpp/build -j
```

Outputs in `cpp/build/`: `simulator`, `simulator_visual`, `ipc_mesh_demo`, `ipc_mesh_visual`, `mesh_check`, `barrier_check`, `cholmod_test`, `libfiction_core.a`, `libfiction_viewer.a`.

Install hints (from CMake warnings): `apt install libeigen3-dev libsuitesparse-dev libglfw3-dev libgl1-mesa-dev` (or brew equivalents).

## Test / Verify

```bash
./cpp/build/simulator        # console validation scenes; expect PASS output
./cpp/build/mesh_check       # msh parser + mass-spring/RegionForce FD checks; expect PASS
ctest --test-dir cpp/build   # simulator_smoke, barrier_check, mesh_check, ipc_mesh_demo_smoke, ...
```

After changing any energy or the integrator, run `simulator` and confirm all scenes pass. The 2D-degradation scene (`upgrade2DTo3D` + z-lock) is a regression test — do not break it.

## Run (IPC .msh models)

`ipc_mesh_demo` / `ipc_mesh_visual` load IPC tet meshes (`include/io/msh_loader.h`, MshIO 4.1 ASCII + `$Surface`), build a mass-spring network on the tet edges with lumped masses `rho * V_tet / 4`, and support gravity, `--pin-box` (DBC) and `--force-box` (constant acceleration region, IPC NBC analog). Spring stiffness defaults to `youngs * avgEdgeLength` (IPC conventions: density 1000, Young's modulus 1e5).

```bash
# 控制台：左端固定 + 右端 +x 拉伸，导出 OBJ 帧
./cpp/build/ipc_mesh_demo --mesh /home/steven/mjx/Archive/IPC/input/tetMeshes/bar-186.msh \
    --pin-box -0.1 -0.1 -0.1 0.6 0.6 1.1 \
    --force-box 9.4 -0.1 -0.1 10.1 0.6 1.1 30 0 0 \
    --no-gravity --steps 100 --export output/bar --export-every 10

# 交互式：重力下垂（需 GLFW；无显示环境会优雅报错退出）
./cpp/build/ipc_mesh_visual --mesh /home/steven/mjx/Archive/IPC/input/tetMeshes/sphere1K.msh \
    --pin-box -2 -2 0.35 2 2 2 --gravity 0 0 -9.81
```

`--help` lists all options (material/poisson/youngs/density/dt/tol/init-vel/pin-box/drive-box/twist-axis/twist-rate/force-box/export/video...). Mass-spring on tet edges is **qualitative only** vs the Neo-Hookean FEM (`--material nh`); quantitative comparison needs the tet FEM energies. Large meshes (sphere5K+, Armadillo*) get slow (seconds/step) — prefer ≤2K nodes for demos.

### Elastogen cantilever experiments (rods.pdf / rods_twist.pdf)

The paper's 16×3×3 cantilever beam is reproduced from `Elastogen/neuralProj/assets/rod*.bgeo` via `tools/bgeo_to_msh.py` → `cpp/tests/data/rod16x3x3.msh`. See `docs/elastogen-rods-experiments.md` for provenance, validation and video recipes (NH bend under gravity + NH twist driven by DBC rotation trajectory).

```bash
# NH 重力弯曲 10s 视频（同 rods.pdf 场景）
./cpp/build/ipc_mesh_demo --mesh cpp/tests/data/rod16x3x3.msh \
    --material nh --youngs 5e7 --poisson 0.32 \
    --pin-box -1 -1 -1 2 3 3 --gravity 0 0 -9.81 \
    --dt 0.002 --steps 5000 --export-every 20 --tol 1e-3 \
    --video output/videos/elastogen_cantilever/bend_rods_nh.mp4
# NH 扭转 10s（rods_twist 场景；DBC 旋转轨迹驱动）
./cpp/build/ipc_mesh_demo --mesh cpp/tests/data/rod16x3x3.msh \
    --material nh --youngs 5e7 --poisson 0.32 \
    --pin-box -1 -1 -1 2 3 3 --drive-box 13 -1 -1 16 3 3 \
    --twist-axis x --twist-rate 36 --no-gravity \
    --dt 0.002 --steps 5000 --export-every 20 --tol 1e-3 \
    --video output/videos/elastogen_cantilever/twist_rods_nh.mp4
```

### Video export (SSH-friendly)

`ipc_mesh_demo --video <path.mp4>` renders frames with a **CPU software rasterizer** (`include/io/rasterizer.h`, no OpenGL/display needed) and encodes with `ffmpeg` (optional; if missing, frames stay in `<path.mp4>.frames/` and the manual command is printed). Frames are sampled every `--export-every` steps, PNG written via zlib (falls back to PPM without zlib). Camera orbits the model by default (`--video-static` to disable).

```bash
# 重力下垂 10s 视频（bar-186 左端固定）
# 时长 = 帧数 / fps，帧数 = steps / export-every + 1；实时 1x 需 export-every * dt * fps = 1
# 即 1000 步 * 0.01s / 每4步一帧 / 25fps -> 251 帧 -> 10s
./cpp/build/ipc_mesh_demo --mesh /home/steven/mjx/Archive/IPC/input/tetMeshes/bar-186.msh \
    --pin-box -0.1 -0.1 -0.1 0.6 0.6 1.1 --gravity 0 0 -9.81 --youngs 1e6 \
    --steps 1000 --export-every 4 --tol 1e-6 --video output/videos/ipc_models/gravity_bar.mp4
```

Dependencies: `ffmpeg` for encoding, zlib for PNG (both optional). Generated demo videos live in `output/videos/<project>/` (gitignored) — **one subdirectory per experiment series**, e.g. `ipc_models/`（IPC .msh 演示）and `elastogen_cantilever/`（论文 cantilever 复现）。New experiment families should get their own subdirectory instead of landing at `output/videos/` top level.

`module_test/` programs are built by hand, e.g.:

```bash
g++ module_test/gradient_hessian.cpp -o module_test/gradient_hessian -O2
```

## Run (Python reference)

```bash
pip install numpy scipy pygame
python Python/simulator.py
```

**Known quirk:** all intra-package imports use the name `Python_slovement` (e.g. `import Python_slovement.square_mesh`) while the folder is named `Python`. If imports fail, either run from a copy/rename of the folder as `Python_slovement`, or fix imports to relative ones. Frames are written to `output/<timestep>.obj`.

## Conventions

- **Adding an energy (C++):** follow the 6-step guide in `cpp/include/energies/README.md`; copy `cpp/include/energies/Gravity.h` or `cpp/include/energies/SpringEnergy.h` as the skeleton. Each energy implements `val` / `addGrad` / `addHess` / `name` / `hessianWeight` / `initStepSize` from `IEnergy` (`cpp/include/energy.h`) and is registered via `addEnergy` in the time integrator.
- **Mesh IO / scenes:** file readers belong in `cpp/include/io/` (header-only, no third-party deps), scene assembly + CLI in `cpp/include/scene/`. `ipc_mesh_demo` and `ipc_mesh_visual` must share `ipc_scene.h` so CLI and physics stay in sync.
- **Hessians** must be PSD-projected via `makePSD` (`cpp/include/utils.h`) before use in the Newton solve; indefinite exact Hessians (e.g. universal gravitation) are projected inside the energy. FD checks must use configurations where the exact Hessian is already PSD (see `mesh_check.cpp` 检验 4).
- **2D-in-3D:** never special-case 2D dimensions; zero-pad positions and use `ZLockEnergy` (`cpp/include/energies/zlock_energy.h`).
- **Linear solve:** CHOLMOD when available, `SimplicialLDLT` fallback — keep both paths working.
- Mirror algorithmic changes in both `cpp/` and `Python/` when the Python file has a counterpart, keeping them numerically equivalent.
- No top-level README exists; `Python/readme.md` describes the original scene.

## Gotchas

- `cpp/CMakeLists.txt` conditionally defines a `time_integrator_demo` target, but `cpp/examples/time_integrator_demo.cpp` does not exist — the target is silently skipped.
- The visual viewer's ground grid is render-only. C++ contact is limited to the plane-barrier `cpp/include/energies/BarrierEnergy.h` (verified numerically by `cpp/examples/barrier_check.cpp`); it is not wired into the viewer demos, and mesh contact/friction exist only in Python. The IPC mesh demos have no contact — models fall through the visual ground plane.
- `cpp/tests/data/*.msh` are small fixtures copied from `IPC/input/tetMeshes/` plus a `cube_no_surface.msh` used to test boundary-face extraction; `rod16x3x3.msh`/`rod_z16x3x3.msh` are **generated** from Elastogen `rod*.bgeo` by `tools/bgeo_to_msh.py` (needs the `elastogen` conda env) — keep regeneration reproducible. `mesh_check` reads the fixtures via the `FICTION_TEST_DATA_DIR` compile definition — keep them small.
- Video frames are written to `<video>.mp4.frames/` and deleted after successful ffmpeg encoding unless `--video-keep-frames`; on failure the directory is kept for manual encoding.
- `cpp/build/`, `__pycache__/`, `output/`, and `module_test/` binaries are gitignored — never commit them.
