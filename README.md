# fiction-simulation-experiments

English | [简体中文](README.zh-CN.md)

A learning and experimentation project built around a **C++20 3D hyperelastic solid simulator**. It follows tetrahedral material energies through global assembly and implicit time integration, using tension, bending, twisting, torus, and independent multibody scenes to study discretization error and computational cost.

- [Project article and experiment figures (Chinese / English)](https://maojunxiang.divinechime64.chatgpt.site/zh/blog/fiction-solid-simulator/)
- [Experiment code and results](EXPERIMENTS.md)
- [Repository maintenance record](docs/repository-audit-20261105.md)

## Implementation

The main project lives in [`simulator/cpp/`](simulator/cpp/):

- Neo-Hookean tetrahedral FEM: material energy, gradient, Hessian, positive semidefinite projection of element Hessians, and global sparse assembly. A mass-spring model is also available for comparison.
- Implicit Euler, Projected Newton, backtracking line search, prescribed-position boundary conditions, and diagnostics for iteration counts, convergence, and physical force residuals.
- Eigen sparse solvers with an optional CHOLMOD backend; an optional GLFW/OpenGL viewer; and headless OBJ, image, and video export.

## Quick start

Run these commands from the repository root. Requirements: a C++20 compiler, CMake 3.20+, and Eigen 3.4+. On Ubuntu:

```bash
sudo apt-get install build-essential cmake libeigen3-dev
```

The minimal build uses the bundled meshes and runs on the CPU:

```bash
cmake -S simulator/cpp -B build -DCMAKE_BUILD_TYPE=Release \
  -DFICTION_USE_VIEWER=OFF -DFICTION_USE_CHOLMOD=OFF
cmake --build build --parallel 3
ctest --test-dir build --output-on-failure
```

### Run a 3D solid example

Use the bundled unit cube, fix its x=0 face, and advance 10 steps under gravity while exporting OBJ files:

```bash
./build/ipc_mesh_demo \
  --mesh simulator/cpp/tests/data/cube.msh \
  --material nh --youngs 1e6 \
  --pin-box -0.01 -0.01 -0.01 0.01 1.01 1.01 \
  --dt 0.005 --steps 10 --tol 1e-6 \
  --export output/quickstart --export-every 5
```

A successful run exits with code 0 and writes `0.obj`, `5.obj`, and `10.obj` to `output/quickstart/`. This example checks the build, material, constraints, time stepping, and export pipeline. If Newton fails to converge or the state contains non-finite values, the console example reports the failed step and exits with a nonzero code.

**The default material is `spring`; pass `--material nh` to run hyperelastic FEM.** Run `./build/ipc_mesh_demo --help` for all options.

### Optional features

| Feature | Dependencies and settings |
|---|---|
| CHOLMOD sparse solver | Install SuiteSparse (`libsuitesparse-dev` on Ubuntu) and configure with `-DFICTION_USE_CHOLMOD=ON`; Eigen is the fallback when CHOLMOD is unavailable |
| Interactive viewer | Install GLFW3/OpenGL (`libglfw3-dev libgl1-mesa-dev` on Ubuntu), configure with `-DFICTION_USE_VIEWER=ON`, and run with a display |
| Images and video | zlib enables PNG output (otherwise PPM); ffmpeg enables MP4 encoding. The console renderer uses CPU rasterization and runs headlessly |

Add `--video output/quickstart.mp4 --video-static` to the example to export video. If ffmpeg is unavailable or encoding fails, the program retains the frame directory and prints a manual encoding command.

## Experiments

| Topic | Recorded work |
|---|---|
| Materials and time integration | Finite-difference material checks, rigid-motion invariance, assembly, and a discrete linear-oscillator reference |
| Tension and static bending | Prescribed-load experiments, stopping-criterion analysis, and spatial mesh refinement, with beam theory as an approximate reference |
| Dynamic bending and twisting | Time-step comparisons and driven-boundary experiments, covering the original S3 0–2 s trajectory comparison with a 2% threshold and subsequent time-step refinement |
| Torus and complex meshes | Torus motion and 10-step computational-cost records for Armadillo 13K/28K/54K |
| Multibody scenes without contact | Short runs with 1, 4, and 16 bodies, comparing uncoupled instances against independent single-body trajectories |

See [EXPERIMENTS.md](EXPERIMENTS.md) for source versions, reports, inputs, and run preparation.

## Layout and run preparation

| Path | Contents |
|---|---|
| [`simulator/cpp/`](simulator/cpp/) | Maintained C++ project, small input meshes, and regression tests |
| [`simulator/tools/`](simulator/tools/) | Structured-bar mesh generation and experiment utilities |
| [`simulator/docs/`](simulator/docs/) | Historical reports with run IDs, parameters, and results |
| [`experiments/workspaces/`](experiments/workspaces/) | Isolated source versions for individual experiment stages |
| [`simulator/Python/`](simulator/Python/) | Historical reference implementation using the original `Python_slovement` import name |
| [`planning/`](planning/) | Validation protocols and execution designs |
| [`assets/ipc/`](assets/ipc/) | Bundled torus and sphere1K meshes with provenance records |

The main build, tests, and cube example run from a checkout. Historical runners require their original machine paths, separate builds, third-party checkouts, or additional meshes. Full raw trajectories, videos, some figures, and Armadillo meshes remain in the original local workspace; historical result links point there.
